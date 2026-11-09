#include "BondVolume.h"
#include "ImagingIO.h"
#include "PlotStyle.h"
#include <Eigen/Cholesky>
#include <Eigen/Eigenvalues>
#include <QCryptographicHash>
#include <QDateTime>
#include <QDir>
#include <QFile>
#include <QFileInfo>
#include <QHash>
#include <QImage>
#include <QPainter>
#include <QPainterPath>
#include <QSaveFile>
#include <QTextStream>
#include <QtEndian>
#include <algorithm>
#include <cstring>
#include <memory>
#include <utility>
#include <numeric>
using namespace ImageUtil;
namespace {
struct RoiData {
    ImagingRoi roi;
    ImageMatrix response, referenceVariance, measured, variance;
};

QJsonObject contrastSourceRevision(const QJsonObject &config) {
    QJsonObject revision;
    const auto input = config["input"].toObject();
    for (const auto &key : {QStringLiteral("sample_nxs"), QStringLiteral("calibration_nxs"),
                            QStringLiteral("roi_table")}) {
        const QFileInfo source(input[key].toString());
        if (!source.isFile())
            fail(QStringLiteral("找不到输入文件：%1").arg(source.filePath()));
        revision[key] = QJsonArray{source.size(), source.lastModified().toMSecsSinceEpoch()};
    }
    return revision;
}

void writeVector(ImagingH5 &file, const QString &path, const QVector<double> &v) {
    file.write(path, {quint64(v.size())}, v.constData(), false);
}
void writeIndex(ImagingH5 &file, const QString &path, const QVector<int> &v,
                const QVector<quint64> &shape = {}) {
    QVector<quint64> data;
    for (int i : v)
        data.append(quint64(i));
    file.writeUInt(path, shape.isEmpty() ? QVector<quint64>{quint64(v.size())} : shape, data.constData(),
                   false);
}
void writeMatrix(ImagingH5 &f, const QString &path, const ImageMatrix &m) {
    f.write(path, {quint64(m.rows()), quint64(m.cols())}, m.data());
}
QJsonObject roiAttributes(const ImagingRoi &r) {
    return {{"label", r.label},
            {"source_id", r.id},
            {"source_label", r.sourceLabel},
            {"detector", r.detector},
            {"placeholder", r.excluded}};
}
void writeGeometry(ImagingH5 &f, const QString &base, const ImagingRoi &r) {
    QVector<quint64> pixels;
    for (auto p : r.pixels) {
        pixels.append(p.y());
        pixels.append(p.x());
    }
    quint64 bounds[] = {quint64(r.rect.y()), quint64(r.rect.y() + r.rect.height()), quint64(r.rect.x()),
                        quint64(r.rect.x() + r.rect.width())};
    f.writeUInt(base + "/bbox_rc", {4}, bounds, false);
    f.writeUInt(base + "/pixels_rc", {quint64(r.pixels.size()), 2}, pixels.constData());
    f.attributes(base, roiAttributes(r));
}
QVector<RoiData> calibration(const QJsonObject &c, const ImagingValidation &v, const QString &directory,
                             const std::atomic_bool &stop, Imaging::Progress progress,
                             const QVector<double> &weights) {
    const auto &g = v.calibration;
    auto acq = c["acquisition"].toObject();
    int nx = int(g.axes[0].size()), ne = int(g.axes[1].size()), width = g.imageSize.width(),
        height = g.imageSize.height(), countPixels = width * height;
    auto rawWeights = Imaging::energyWeights(g.axes[1]);
    quint64 count = 0;
    for (const auto &roi : v.rois)
        count += quint64(ne) * roi.pixels.size() * nx;
    if (count > 134217728 || (!directory.isEmpty() && quint64(nx) * countPixels > 134217728))
        fail("Calibration pixel/grid buffers exceed 1 GiB; use smaller ROIs or grid");
    QVector<std::vector<quint64>> counts;
    for (const auto &roi : v.rois)
        counts.append(std::vector<quint64>(size_t(ne) * roi.pixels.size() * nx, 0));
    QVector<double> rawImage(directory.isEmpty() ? 0 : countPixels, 0),
        band(directory.isEmpty() ? 0 : qsizetype(nx) * countPixels, 0);
    double *rawSum = rawImage.data(), *bandSum = band.data();
    ImagingH5 source(g.file);
    source.frames("/entry/data/" + acq["detector"].toString(),
                  c["runtime"].toObject()["chunk_frames"].toInt(), stop,
                  [&](int f, const quint64 *pixels, int, int) {
                      int x = g.indices[0][f], e = g.indices[1][f];
                      if (!rawImage.isEmpty()) {
                          double factor = rawWeights[e] / g.dose[x * ne + e];
                          auto image = bandSum + qsizetype(x) * countPixels;
                          for (int p = 0; p < countPixels; ++p) {
                              rawSum[p] += double(pixels[p]);
                              image[p] += factor * double(pixels[p]);
                          }
                      }
                      for (int r = 0; r < v.rois.size(); ++r) {
                          const auto &roi = v.rois[r];
                          auto *buffer = counts[r].data();
                          qsizetype size = roi.pixels.size();
                          for (qsizetype p = 0; p < size; ++p) {
                              auto point = roi.pixels[p];
                              auto &target = buffer[(qsizetype(e) * size + p) * nx + x];
                              quint64 value = pixels[point.y() * width + point.x()];
                              if (target > std::numeric_limits<quint64>::max() - value)
                                  fail("Pooled photon count overflow");
                              target += value;
                          }
                      }
                      if (progress && f % 32 == 0)
                          progress(QStringLiteral("薄膜像素提取 %1 / %2").arg(f + 1).arg(g.monitor.size()),
                                   5 + f * 20 / int(g.monitor.size()));
                  });
    std::unique_ptr<ImagingH5> out;
    if (!directory.isEmpty()) {
        out = std::make_unique<ImagingH5>(directory + "/calibration.h5", true);
        out->attributes("/", {{"schema", "XRI.calibration.v1"},
                              {"source_file", g.file},
                              {"roi_file", c["input"].toObject()["roi_table"]},
                              {"monitor", acq["monitor"]},
                              {"exposure_s", acq["calibration_exposure_s"]},
                              {"foil_thickness_mm", acq["foil_thickness_mm"]},
                              {"response_units", "counts / (monitor_current * s) * keV"},
                              {"background_subtracted", false},
                              {"absolute_density_calibrated", false},
                              {"roi_stop_convention", "exclusive"},
                              {"calibration_status", "provisional: inspect ROI tails and background"}});
        writeVector(*out, "x_mm", g.axes[0]);
        writeVector(*out, "energy_kev", g.axes[1]);
        writeVector(*out, "energy_weights_kev", rawWeights);
        writeIndex(*out, "occupancy_xe", g.occupancy, {quint64(nx), quint64(ne)});
        ImageMatrix doses(ne, nx);
        for (int e = 0; e < ne; ++e)
            for (int x = 0; x < nx; ++x)
                doses(e, x) = g.dose[x * ne + e];
        writeMatrix(*out, "dose_ex", doses);
        writeIndex(*out, "frames/x_index", g.indices[0]);
        writeIndex(*out, "frames/energy_index", g.indices[1]);
        writeVector(*out, "frames/x_setpoint_mm", g.coordinates[0]);
        writeVector(*out, "frames/energy_setpoint_kev", g.coordinates[1]);
        writeVector(*out, "frames/i0", g.monitor);
        writeVector(*out, "frames/exposure_s",
                    QVector<double>(g.monitor.size(), acq["calibration_exposure_s"].toDouble()));
        writeVector(*out, "frames/x_readback_mm",
                    source.contains("/entry/instrument/M_Stage_X")
                        ? source.read("/entry/instrument/M_Stage_X")
                        : g.coordinates[0]);
        writeVector(*out, "frames/energy_readback_kev",
                    source.contains("/entry/instrument/M_DCM_B5Energy")
                        ? source.read("/entry/instrument/M_DCM_B5Energy")
                        : g.coordinates[1]);
        out->write("detectors/" + acq["detector"].toString() + "/raw_sum_rc",
                   {quint64(height), quint64(width)}, rawImage.constData());
        out->write("detectors/" + acq["detector"].toString() + "/band_response_xrc",
                   {quint64(nx), quint64(height), quint64(width)}, band.constData());
    }
    QVector<RoiData> selected;
    for (int r = 0; r < v.rois.size(); ++r) {
        cancel(stop);
        const auto &roi = v.rois[r];
        int np = int(roi.pixels.size());
        auto base = "rois/" + roi.key;
        if (out)
            out->writeUInt(base + "/counts_epx", {quint64(ne), quint64(np), quint64(nx)}, counts[r].data());
        if (out)
            writeGeometry(*out, base, roi);
        ImageMatrix raw = ImageMatrix::Zero(np, nx), rawVar = raw, matched = raw, matchedVar = raw;
        const auto *buffer = counts[r].data();
        for (int e = 0; e < ne; ++e)
            for (int p = 0; p < np; ++p)
                for (int x = 0; x < nx; ++x) {
                    double value = double(buffer[(qsizetype(e) * np + p) * nx + x]), d = g.dose[x * ne + e],
                           f = rawWeights[e] / d, m = weights[e] / d;
                    raw(p, x) += f * value;
                    rawVar(p, x) += f * f * value;
                    matched(p, x) += m * value;
                    matchedVar(p, x) += m * m * value;
                }
        if (out) {
            writeMatrix(*out, base + "/R_band_px", raw);
            writeMatrix(*out, base + "/R_band_variance_px", rawVar);
        }
        std::vector<quint64>().swap(counts[r]);
        if (!roi.excluded)
            selected.append({roi, std::move(matched), std::move(matchedVar), {}, {}});
    }
    return selected;
}
void sample(const QJsonObject &c, const ImagingValidation &v, QVector<RoiData> &selected,
            const QString &directory, const std::atomic_bool &stop, Imaging::Progress progress) {
    const auto &g = v.sample;
    auto acq = c["acquisition"].toObject();
    int ny = int(g.axes[0].size()), nz = int(g.axes[1].size()), ne = int(g.axes[2].size()),
        width = g.imageSize.width();
    const auto weights = directory.isEmpty() ? qMakePair(QVector<double>{}, Imaging::energyWeights(g.axes[2]))
                                             : Imaging::matchedWeights(v.calibration.axes[1], g.axes[2]);
    quint64 size = 0;
    for (const auto &r : selected)
        size += r.roi.pixels.size() * quint64(ny) * nz;
    if (size > 134217728)
        fail("Sample ROI buffers exceed memory limit");
    for (auto &r : selected) {
        r.measured = ImageMatrix::Zero(r.roi.pixels.size(), ny * nz);
        r.variance = r.measured;
    }
    ImagingH5 source(g.file);
    source.frames("/entry/data/" + acq["detector"].toString(),
                  c["runtime"].toObject()["chunk_frames"].toInt(), stop,
                  [&](int frame, const quint64 *pixels, int, int) {
                      int yz = g.indices[0][frame] * nz + g.indices[1][frame], e = g.indices[2][frame];
                      double factor = weights.second[e] / g.dose[yz * ne + e];
                      for (auto &r : selected) {
                          for (int p = 0; p < r.roi.pixels.size(); ++p) {
                              auto point = r.roi.pixels[p];
                              double value = double(pixels[point.y() * width + point.x()]);
                              r.measured(p, yz) += factor * value;
                              r.variance(p, yz) += factor * factor * value;
                          }
                      }
                      if (progress && frame % 64 == 0)
                          progress(QStringLiteral("样品积分 %1 / %2").arg(frame + 1).arg(g.monitor.size()),
                                   25 + frame * 20 / int(g.monitor.size()));
                  });
    if (directory.isEmpty())
        return;
    ImagingH5 out(directory + "/tomography_inputs.h5", true);
    bool canonical = true;
    for (int f = 0; f < g.monitor.size(); ++f)
        if ((g.indices[0][f] * nz + g.indices[1][f]) * ne + g.indices[2][f] != f)
            canonical = false;
    out.attributes("/",
                   {{"schema", "XRI.inputs.v1"},
                    {"source_file", g.file},
                    {"calibration_file", directory + "/calibration.h5"},
                    {"exposure_s", acq["sample_exposure_s"]},
                    {"monitor", acq["monitor"]},
                    {"expected_frames", ny * nz * ne},
                    {"actual_frames", g.monitor.size()},
                    {"canonical_order_verified", canonical},
                    {"background_subtracted", false},
                    {"integration",
                     "linear calibration spectral interpolation then sample-grid trapezoidal integration"}});
    writeVector(out, "x_mm", v.calibration.axes[0]);
    writeVector(out, "y_mm", g.axes[0]);
    writeVector(out, "z_mm", g.axes[1]);
    writeVector(out, "energy_kev", g.axes[2]);
    writeVector(out, "energy_weights_kev", weights.second);
    writeVector(out, "reference_effective_energy_weights_kev", weights.first);
    writeIndex(out, "occupancy_yze", g.occupancy, {quint64(ny), quint64(nz), quint64(ne)});
    writeIndex(out, "complete_spectrum_yz", QVector<int>(ny * nz, 1), {quint64(ny), quint64(nz)});
    out.write("dose_yze", {quint64(ny), quint64(nz), quint64(ne)}, g.dose.constData(), false);
    writeIndex(out, "frame_y_index", g.indices[0]);
    writeIndex(out, "frame_z_index", g.indices[1]);
    writeIndex(out, "frame_energy_index", g.indices[2]);
    writeVector(out, "i0", g.monitor);
    QStringList keys;
    for (const auto &r : selected)
        keys.append(r.roi.key);
    out.writeStrings("roi_keys", keys);
    for (const auto &r : selected) {
        auto base = "rois/" + r.roi.key;
        writeMatrix(out, base + "/R_matched_px", r.response);
        writeMatrix(out, base + "/R_matched_variance_px", r.referenceVariance);
        out.write(base + "/D_band_pyz", {quint64(r.measured.rows()), quint64(ny), quint64(nz)},
                  r.measured.data());
        out.write(base + "/D_variance_pyz", {quint64(r.variance.rows()), quint64(ny), quint64(nz)},
                  r.variance.data());
        writeGeometry(out, base, r.roi);
    }
}

struct SpatialSolution {
    QVector<RoiData> accepted;
    QStringList keys;
    QVector<ImageMatrix> response, operators;
    QVector<double> gains, history;
    ImageMatrix h, b, rho;
    ImagingSolverInfo nnlsInfo;
    bool gainConverged;
    int gainIteration;
    double dx;
};

SpatialSolution solveSpatial(const QJsonObject &c, const ImagingValidation &validation,
                             QVector<RoiData> &selected, ImagingResult &result, const std::atomic_bool &stop,
                             const Imaging::Progress &log, qsizetype energyChannels = 1) {
    const auto &progress = log;
    QVector<RoiData> accepted;
    QVector<ImagingPsf> models;
    auto acq = c["acquisition"].toObject(), rec = c["reconstruction"].toObject();
    double pitch = acq["detector_pixel_pitch_mm"].toDouble();
    for (int i = 0; i < selected.size(); ++i) {
        cancel(stop);
        try {
            auto model = Imaging::fitPsf(selected[i].response, selected[i].roi.pixels,
                                         validation.calibration.axes[0], c["psf"].toObject(), stop);
            result.models[selected[i].roi.key] = model.json(pitch);
            models.append(model);
            accepted.append(std::move(selected[i]));
        } catch (const std::exception &e) {
            cancel(stop);
            result.rejected[selected[i].roi.key] = QString::fromUtf8(e.what());
            log("PSF rejected " + selected[i].roi.key + ": " + QString::fromUtf8(e.what()),
                46 + i * 14 / int(selected.size()));
        }
        if (progress)
            progress(QStringLiteral("PSF %1 / %2").arg(i + 1).arg(selected.size()),
                     46 + (i + 1) * 14 / int(selected.size()));
    }
    if (!result.directory.isEmpty())
        saveJson(result.directory + "/psf_models.json",
                 {{"models", result.models}, {"rejected", result.rejected}});
    if (accepted.isEmpty())
        fail("All PSF models failed quality checks");
    QStringList keys;
    for (const auto &r : accepted)
        keys.append(r.roi.key);
    int anchor = keys.indexOf(rec["anchor_roi"].toString());
    if (rec["fit_gains"].toBool() && anchor < 0)
        fail("Gain anchor failed PSF quality checks; choose a verified anchor explicitly: " +
             rec["anchor_roi"].toString());
    QVector<double> lengths;
    for (const auto &m : models)
        lengths.append(m.sampleMmPerPixel);
    auto originalX = validation.calibration.axes[0];
    double dx = c["grid"].toObject()["x_step_mm"].isNull() ? percentile(lengths, 50)
                                                           : c["grid"].toObject()["x_step_mm"].toDouble();
    double cells = std::floor((originalX.last() - originalX.first()) / dx + 1e-9) + 1;
    if (cells < 2 || cells > 4096)
        fail("Sample-space X grid requires 2..4096 cells; adjust x_step_mm");
    int nx = int(cells), ny = int(validation.sample.axes[0].size()),
        nz = int(validation.sample.axes[1].size());
    if (qsizetype(nx) * ny * nz > 16000000 / energyChannels)
        fail(energyChannels > 1 ? QStringLiteral("四维体素超过 1600 万。请增大 X 步长或减少扫描范围。")
                                : QStringLiteral("Reconstructed volume exceeds memory limit"));
    QVector<double> x;
    for (int i = 0; i < nx; ++i)
        x.append((originalX.first() + originalX.last()) / 2 + dx * (i - (nx - 1) / 2.));
    result.axes = {x, validation.sample.axes[0], validation.sample.axes[1]};
    log(QString("Full grid: %1 X x %2 Y x %3 Z; dx=%4 mm; %5 accepted ROIs")
            .arg(nx)
            .arg(ny)
            .arg(nz)
            .arg(dx, 0, 'g', 12)
            .arg(accepted.size()),
        60);
    double scale = dx / acq["foil_thickness_mm"].toDouble();
    QVector<ImageMatrix> response, operators, grams, backs;
    for (int i = 0; i < accepted.size(); ++i) {
        response.append(Imaging::evaluatePsf(models[i], accepted[i].roi.pixels, x));
        operators.append(response.last() * scale);
        grams.append(operators.last().transpose() * operators.last());
        backs.append(operators.last().transpose() * accepted[i].measured);
    }
    QVector<double> gains(accepted.size(), 1), history;
    ImageMatrix h, b, rho;
    auto terms = [&] {
        h = ImageMatrix::Zero(nx, nx);
        b = ImageMatrix::Zero(nx, ny * nz);
        for (int i = 0; i < gains.size(); ++i) {
            h.noalias() += gains[i] * gains[i] * grams[i];
            b.noalias() += gains[i] * backs[i];
        }
    };
    ImagingSolverInfo nnlsInfo;
    bool gainConverged = false;
    int gainIteration = 0;
    log("4/6 Joint NNLS and anchored relative crystal gains", 62);
    if (rec["fit_gains"].toBool()) {
        auto bounds = rec["gain_bounds"].toArray();
        for (int iteration = 1; iteration <= rec["max_gain_iterations"].toInt(); ++iteration) {
            cancel(stop);
            terms();
            rho = Imaging::nnls(h, b, rec["nnls_max_iterations"].toInt(), rec["nnls_tolerance"].toDouble(),
                                stop, &nnlsInfo, rho);
            double difference = 0, norm = 0;
            auto next = gains;
            for (int i = 0; i < gains.size(); ++i) {
                double numerator = (backs[i].array() * rho.array()).sum(),
                       denominator = (rho.array() * (grams[i] * rho).array()).sum();
                if (denominator > 0)
                    next[i] = numerator / denominator;
                next[i] = std::clamp(next[i], bounds[0].toDouble(), bounds[1].toDouble());
                if (i == anchor)
                    next[i] = 1;
                difference += std::pow(next[i] - gains[i], 2);
                norm += gains[i] * gains[i];
            }
            double change = std::sqrt(difference) / std::max(std::sqrt(norm), 1e-15);
            history.append(change);
            gains = std::move(next);
            gainIteration = iteration;
            if (iteration % 20 == 0 || iteration == 1)
                log(QString("Gain iteration %1: relative change %2").arg(iteration).arg(change, 0, 'g', 7),
                    65);
            if (change < rec["gain_tolerance"].toDouble()) {
                gainConverged = true;
                break;
            }
        }
        terms();
        rho = Imaging::nnls(h, b, rec["nnls_max_iterations"].toInt(), rec["nnls_tolerance"].toDouble(), stop,
                            &nnlsInfo, rho);
        gainConverged = gainConverged && nnlsInfo.converged;
    } else {
        terms();
        rho = Imaging::nnls(h, b, rec["nnls_max_iterations"].toInt(), rec["nnls_tolerance"].toDouble(), stop,
                            &nnlsInfo);
        gainConverged = nnlsInfo.converged;
    }
    return {std::move(accepted),
            std::move(keys),
            std::move(response),
            std::move(operators),
            std::move(gains),
            std::move(history),
            std::move(h),
            std::move(b),
            std::move(rho),
            nnlsInfo,
            gainConverged,
            gainIteration,
            dx};
}
void exportVolume(const QString &path, const ImageMatrix &rho, const QVector<QVector<double>> &axes,
                  const QString &directory, const QString &method, const QJsonObject &metadata,
                  const QVector<double> &sensitivity, const QVector<double> &eigenvalues) {
    ImagingH5 out(path, true);
    out.write("reinterp_insect_ck",
              {quint64(axes[0].size()), quint64(axes[1].size()), quint64(axes[2].size())}, rho.data(), true,
              true);
    out.attributes("reinterp_insect_ck",
                   {{"axes", "X,Y,Z"},
                    {"quantity", "relative Kapton-equivalent elastic scattering density"},
                    {"units", "relative"},
                    {"background_subtracted", false},
                    {"energy_integrated", true},
                    {"absolute_density_calibrated", false}});
    const QStringList names{"x_mm", "y_mm", "z_mm"};
    for (int d = 0; d < 3; ++d)
        writeVector(out, names[d], axes[d]);
    out.scales("reinterp_insect_ck", names);
    QJsonObject attributes = metadata;
    attributes["source_inputs"] = directory + "/tomography_inputs.h5";
    attributes["solver"] = method;
    attributes["schema"] = "XRI.result.v1";
    out.attributes("/", attributes);
    writeVector(out, "sensitivity_x", sensitivity);
    out.attributes("sensitivity_x",
                   {{"definition", "sqrt(sum_a,p A_a[p,x]^2); includes fitted gains and dx/foil thickness"}});
    writeVector(out, "gram_eigenvalues", eigenvalues);
}

void exportAvizo(const QString &path, const ImagingResult &result, const std::atomic_bool &stop) {
    const auto &axes = result.axes;
    const qsizetype nx = axes[0].size(), ny = axes[1].size(), nz = axes[2].size();
    QStringList bounds;
    for (const auto &axis : axes) {
        bounds.append(QString::number(axis.first(), 'g', 17));
        bounds.append(QString::number(axis.last(), 'g', 17));
    }
    const auto header = QStringLiteral("# AmiraMesh BINARY-LITTLE-ENDIAN 2.1\n\n"
                                       "define Lattice %1 %2 %3\n\n"
                                       "Parameters {\n"
                                       "    Units { Coordinates \"mm\" },\n"
                                       "    BoundingBox %4,\n"
                                       "    CoordType \"uniform\"\n"
                                       "}\n\n"
                                       "Lattice { float Data } @1\n\n"
                                       "# Data section follows\n"
                                       "@1\n")
                            .arg(nx).arg(ny).arg(nz).arg(bounds.join(' ')).toLatin1();
    cancel(stop);
    QSaveFile file(path);
    if (!file.open(QIODevice::WriteOnly) || file.write(header) != header.size())
        fail(QStringLiteral("无法写入 Avizo 文件 %1。%2").arg(path, file.errorString()));

    // AmiraMesh stores X contiguously. XRI stores Z contiguously within each X plane.
    const qsizetype rowsPerBlock = 16384 / nx;
    QVector<quint32> block(nx * rowsPerBlock);
    for (qsizetype z = 0; z < nz; ++z) {
        for (qsizetype firstY = 0; firstY < ny; firstY += rowsPerBlock) {
            cancel(stop);
            const qsizetype endY = std::min(ny, firstY + rowsPerBlock);
            qsizetype index = 0;
            for (qsizetype y = firstY; y < endY; ++y) {
                for (qsizetype x = 0; x < nx; ++x) {
                    const float value = float(result.primary(x, y * nz + z));
                    quint32 bits;
                    std::memcpy(&bits, &value, sizeof(bits));
                    block[index++] = qToLittleEndian(bits);
                }
            }
            const qint64 bytes = index * qint64(sizeof(quint32));
            if (file.write(reinterpret_cast<const char *>(block.constData()), bytes) != bytes)
                fail(QStringLiteral("无法写入 Avizo 文件 %1。%2").arg(path, file.errorString()));
        }
    }
    cancel(stop);
    if (!file.commit())
        fail(QStringLiteral("无法保存 Avizo 文件 %1。%2").arg(path, file.errorString()));
}

QImage diagnosticImage(const ImagingResult &result, bool bondBand = false) {
    auto color = [](double t) {
        const std::array<QColor, 5> stops{QColor("#440154"), QColor("#3b528b"), QColor("#21918c"),
                                          QColor("#5ec962"), QColor("#fde725")};
        double u = std::clamp(t, 0., 1.) * 4;
        int i = std::min(3, int(u));
        double a = u - i;
        return QColor::fromRgbF(float(stops[i].redF() * (1 - a) + stops[i + 1].redF() * a),
                                float(stops[i].greenF() * (1 - a) + stops[i + 1].greenF() * a),
                                float(stops[i].blueF() * (1 - a) + stops[i + 1].blueF() * a));
    };
    QImage canvas(3000, 1960, QImage::Format_RGB32);
    canvas.fill(Qt::white);
    QPainter p(&canvas);
    p.scale(2, 2);
    p.setRenderHint(QPainter::Antialiasing);
    p.setFont(QFont("Arial", 14));
    int nx = int(result.axes[0].size()), ny = int(result.axes[1].size()), nz = int(result.axes[2].size());
    std::array<ImageMatrix, 3> projections{ImageMatrix::Zero(ny, nz), ImageMatrix::Zero(nx, nz),
                                           ImageMatrix::Zero(nx, ny)};
    for (int x = 0; x < nx; ++x)
        for (int y = 0; y < ny; ++y)
            for (int z = 0; z < nz; ++z) {
                double v = result.primary(x, y * nz + z);
                projections[0](y, z) += v;
                projections[1](x, z) += v;
                projections[2](x, y) += v;
            }
    QStringList titles{"YZ projection (sum X)", "XZ projection (sum Y)", "XY projection (sum Z)"};
    for (int k = 0; k < 3; ++k) {
        auto &m = projections[k];
        QImage image(int(m.cols()), int(m.rows()), QImage::Format_RGB32);
        double maximum = m.maxCoeff();
        for (int row = 0; row < m.rows(); ++row)
            for (int col = 0; col < m.cols(); ++col) {
                double v = maximum > 0 ? m(row, col) / maximum : 0;
                QColor c = color(v);
                image.setPixelColor(col, row, c);
            }
        int rowAxis = k == 0 ? 1 : 0, colAxis = k == 2 ? 1 : 2;
        const auto &xs = result.axes[colAxis], &ys = result.axes[rowAxis];
        double dx = xs.size() > 1 ? xs[1] - xs[0] : 1, dy = ys.size() > 1 ? ys[1] - ys[0] : 1;
        double x0 = xs.first() - .5 * dx, x1 = xs.last() + .5 * dx, y0 = ys.first() - .5 * dy,
               y1 = ys.last() + .5 * dy;
        QRect slot(65 + k * 475, 80, 320, 340);
        double scale = std::min(slot.width() / (x1 - x0), slot.height() / (y1 - y0));
        QSize size(std::max(1, int((x1 - x0) * scale)), std::max(1, int((y1 - y0) * scale)));
        QRect area(QPoint(), size);
        area.moveCenter(slot.center());
        p.drawImage(area, image);
        p.setPen(Qt::black);
        p.drawRect(area);
        for (int y = 0; y < area.height(); ++y) {
            double fraction = 1 - double(y) / std::max(1, area.height() - 1);
            p.setPen(color(fraction));
            p.drawLine(area.right() + 12, area.top() + y, area.right() + 25, area.top() + y);
        }
        p.setPen(Qt::black);
        auto colorTicks = PlotStyle::ticks(0, maximum, 4);
        for (double v : colorTicks.major) {
            double y = area.bottom() - v / maximum * area.height();
            p.drawText(QRectF(area.right() + 30, y - 10, 90, 20), Qt::AlignVCenter,
                       PlotStyle::label(v, colorTicks));
        }
        p.drawText(QRect(area.x() - 20, 35, 420, 30), Qt::AlignCenter, titles[k]);
        PlotStyle::axis(p, area, PlotStyle::ticks(x0, x1, 4), x0, x1, true, false, false, true, 12);
        PlotStyle::axis(p, area, PlotStyle::ticks(y0, y1, 4), y0, y1, false, true, false, true, 12);
        p.drawText(QRect(area.x(), area.bottom() + 33, area.width(), 28), Qt::AlignCenter,
                   QString("%1 (mm)").arg(colAxis == 1 ? "Y" : "Z"));
        p.save();
        p.translate(area.x() - 57, area.center().y());
        p.rotate(-90);
        p.drawText(QRect(-160, -15, 320, 30), Qt::AlignCenter,
                   QString("%1 (mm)").arg(rowAxis == 0 ? "X" : "Y"));
        p.restore();
    }
    auto curve = [&](const QVector<double> &y, const QRectF &area, const QString &title) {
        double max = *std::max_element(y.begin(), y.end());
        p.setPen(Qt::black);
        p.drawRect(area);
        const double ymax = max > 0 ? max * 1.05 : 1;
        const double xlo = result.axes[0].first(), xhi = result.axes[0].last();
        PlotStyle::axis(p, area, PlotStyle::ticks(xlo, xhi, 4), xlo, xhi, true, false, true, true, 12);
        PlotStyle::axis(p, area, PlotStyle::ticks(0, ymax, 4), 0, ymax, false, false, true, true, 12);
        p.setPen(Qt::black);
        p.drawText(QRectF(area.x(), area.y() - 36, area.width(), 28), Qt::AlignCenter, title);
        QPainterPath line;
        for (int i = 0; i < y.size(); ++i) {
            QPointF point(area.left() + double(i) / std::max(1, int(y.size() - 1)) * area.width(),
                          area.bottom() - y[i] / ymax * area.height());
            if (i)
                line.lineTo(point);
            else
                line.moveTo(point);
        }
        p.setPen(QPen(QColor(0, 80, 170), 1.5));
        p.drawPath(line);
        p.setPen(Qt::black);
        p.drawText(QRectF(area.left(), area.bottom() + 33, area.width(), 28), Qt::AlignCenter,
                   QStringLiteral("X (mm)"));
    };
    curve(result.sensitivity, {70, 560, 350, 250}, "Full-X sensitivity (no cropping)");
    QVector<double> integrated;
    for (int x = 0; x < nx; ++x)
        integrated.append(result.primary.row(x).sum());
    curve(integrated, {1020, 560, 350, 250}, "X-integrated density");
    QRectF trajectory(545, 560, 350, 250);
    p.setPen(Qt::black);
    p.drawRect(trajectory);
    p.drawText(QRectF(510, 524, 420, 28), Qt::AlignCenter, "Selected centroid trajectories");
    double lo = INFINITY, hi = -INFINITY, xlo = INFINITY, xhi = -INFINITY;
    int shown = 0;
    for (auto it = result.models.begin(); it != result.models.end() && shown < 4; ++it, ++shown)
        for (auto value : it->toObject()["observations"].toArray()) {
            auto obs = value.toObject();
            lo = std::min(lo, obs["row"].toDouble());
            hi = std::max(hi, obs["row"].toDouble());
            xlo = std::min(xlo, obs["x"].toDouble());
            xhi = std::max(xhi, obs["x"].toDouble());
        }
    if (hi <= lo)
        hi = lo + 1;
    if (xhi <= xlo)
        xhi = xlo + 1;
    PlotStyle::axis(p, trajectory, PlotStyle::ticks(xlo, xhi, 4), xlo, xhi, true, false, true, true, 12);
    PlotStyle::axis(p, trajectory, PlotStyle::ticks(lo, hi, 4), lo, hi, false, false, true, true, 12);
    shown = 0;
    for (auto it = result.models.begin(); it != result.models.end() && shown < 4; ++it, ++shown) {
        auto model = it->toObject();
        auto coefficients = model["row_line"].toArray();
        QColor ink = PlotStyle::curveColor(shown);
        p.setPen(QPen(ink, 1.5));
        QPainterPath line;
        bool first = true;
        for (auto value : model["observations"].toArray()) {
            auto obs = value.toObject();
            double pos = obs["x"].toDouble(), row = obs["row"].toDouble(),
                   prediction = coefficients[0].toDouble() * pos + coefficients[1].toDouble();
            auto point = [&](double r) {
                return QPointF(
                    trajectory.left() + (pos - xlo) / std::max(xhi - xlo, 1e-15) * trajectory.width(),
                    trajectory.bottom() - (r - lo) / std::max(hi - lo, 1e-15) * trajectory.height());
            };
            p.drawEllipse(point(row), 2, 2);
            if (first) {
                line.moveTo(point(prediction));
                first = false;
            } else
                line.lineTo(point(prediction));
        }
        p.drawPath(line);
        p.drawText(QRectF(trajectory.left() + (shown % 2) * 175, trajectory.bottom() + 64 + (shown / 2) * 20,
                          170, 20),
                   it.key());
    }
    p.setPen(Qt::black);
    p.drawText(QRectF(trajectory.left(), trajectory.bottom() + 33, trajectory.width(), 28), Qt::AlignCenter,
               QStringLiteral("Foil X (mm)"));
    p.save();
    p.translate(trajectory.left() - 57, trajectory.center().y());
    p.rotate(-90);
    p.drawText(QRectF(-125, -12, 250, 24), Qt::AlignCenter, "Detector row (pixel)");
    p.restore();
    p.drawText(QRect(40, 925, 1420, 45), Qt::AlignCenter | Qt::TextWordWrap,
               bondBand
                   ? "Integrated-band fit used to determine shared crystal gains. Inspect individual "
                     "energy channels in the four-dimensional HDF5 result."
                   : "Provisional relative Kapton-equivalent density. No sample background subtraction or "
                     "geometric postprocessing. Inspect PSF, sensitivity, rank and residuals together.");
    p.end();
    canvas.setDotsPerMeterX(11811);
    canvas.setDotsPerMeterY(11811);
    return canvas;
}
} // namespace

ImagingValidation BondImaging::validate(const QJsonObject &c, const std::atomic_bool &stop) {
    auto out = Imaging::validate(c, stop, false);
    const auto &foil = out.calibration.axes[1], &sample = out.sample.axes[2];
    const bool matched = sample.first() >= foil.first() - 1e-10 && sample.last() <= foil.last() + 1e-10;
    out.report["psf_energy_band"] = matched ? "sample-matched" : "foil-native";
    return out;
}

BondProjection BondImaging::prepareContrast(const QJsonObject &supplied, double ea, double eb,
                                                const ContrastSpectra &prepareSpectra, const std::atomic_bool &stop,
                                                Imaging::Progress progress) {
    const auto c = Imaging::config(supplied, QDir::currentPath());
    const auto log = [&](const QString &message, int percent) {
        cancel(stop);
        if (progress)
            progress(message, percent);
    };
    log(QStringLiteral("校验扫描与 ROI"), 0);
    const auto revision = contrastSourceRevision(c);
    auto inputs = c;
    auto reconstruction = inputs["reconstruction"].toObject();
    reconstruction["fit_gains"] = false;
    inputs["reconstruction"] = reconstruction;
    auto validation = validate(inputs, stop);
    const auto &g = validation.sample;
    const auto &energy = g.axes[2];
    const auto difference = contrastWeights(energy, ea, eb, {});
    QStringList labels;
    QVector<const ImagingRoi *> rois;
    quint64 elements = 0;
    const qsizetype ny = g.axes[0].size(), nz = g.axes[1].size(), ne = energy.size();
    for (const auto &roi : std::as_const(validation.rois)) {
        if (roi.excluded)
            continue;
        labels.append(roi.label);
        rois.append(&roi);
        elements += quint64(roi.pixels.size()) * ny * nz;
    }
    if (rois.isEmpty())
        fail(QStringLiteral("请至少纳入一个有效晶体。"));
    if (elements > 67108864)
        fail(QStringLiteral("差值投影缓存超过 512 MiB，请缩小 ROI 或减少参与的晶体。"));
    log(QStringLiteral("准备谱处理和双能量差值"), 5);
    const auto spectra = prepareSpectra ? prepareSpectra(energy, labels)
                                        : QVector<BondSpectralSettings>(labels.size());
    if (spectra.size() != labels.size())
        fail(QStringLiteral("谱处理参数与晶体数不匹配。"));
    QVector<QVector<double>> coefficients;
    BondProjection out;
    QVector<ImageMatrix *> projections;
    for (qsizetype i = 0; i < spectra.size(); ++i) {
        const auto &settings = spectra[i];
        auto weights = settings.background.isEmpty() && !settings.compton && !settings.normalization
                           ? difference
                           : contrastWeights(energy, ea, eb, settings);
        if (std::any_of(weights.cbegin(), weights.cend(), [](double v) { return !std::isfinite(v); }))
            fail(QStringLiteral("谱处理权重与实测能量不匹配。"));
        coefficients.append(std::move(weights));
        auto stored = out.differenceByCrystal.try_emplace(rois[i]->key, rois[i]->pixels.size(), ny * nz).first;
        stored->second.setZero();
        projections.append(&stored->second);
    }
    const auto acq = c["acquisition"].toObject();
    ImagingH5 source(g.file);
    source.frames("/entry/data/" + acq["detector"].toString(),
                  c["runtime"].toObject()["chunk_frames"].toInt(), stop,
                  [&](int frame, const quint64 *pixels, int, int) {
                      const qsizetype yz = g.indices[0][frame] * nz + g.indices[1][frame];
                      const int e = g.indices[2][frame];
                      const double dose = g.dose[yz * ne + e];
                      for (qsizetype i = 0; i < rois.size(); ++i) {
                          const double factor = coefficients[i][e] / dose;
                          if (factor == 0.)
                              continue;
                          const auto &points = rois[i]->pixels;
                          auto &projection = *projections[i];
                          for (qsizetype p = 0; p < points.size(); ++p) {
                              const auto pixel = points[p];
                              projection(p, yz) += factor * double(pixels[pixel.y() * g.imageSize.width() + pixel.x()]);
                          }
                      }
                      if (frame % 128 == 0)
                          log(QStringLiteral("计算探测器差值投影 %1 / %2").arg(frame + 1).arg(g.monitor.size()),
                              10 + 89 * frame / int(g.monitor.size()));
                  });
    out.config = c;
    out.ea = ea;
    out.eb = eb;
    if (contrastSourceRevision(c) != revision)
        fail(QStringLiteral("计算期间输入文件已变化，请重新校验后计算差值投影。"));
    out.sourceRevision = revision;
    out.validation = std::move(validation);
    log(QStringLiteral("差值投影完成"), 100);
    return out;
}

BondContrast BondImaging::reconstructContrast(const BondProjection &projection,
                                              const QJsonObject &spatialSettings, const std::atomic_bool &stop,
                                              Imaging::Progress progress) {
    auto settings = projection.config;
    for (const auto &group : {QStringLiteral("grid"), QStringLiteral("reconstruction")})
        settings[group] = spatialSettings.value(group);
    const auto c = Imaging::config(settings, QDir::currentPath());
    const auto rec = c["reconstruction"].toObject();
    if (rec["fit_gains"].toBool() && projection.differenceByCrystal.count(rec["anchor_roi"].toString()) == 0)
        fail(QStringLiteral("增益锚点必须是差值投影中纳入的晶体。"));
    const auto log = [&](const QString &message, int percent) {
        cancel(stop);
        if (progress)
            progress(message, percent);
    };
    log(QStringLiteral("准备薄膜空间标定"), 0);
    if (contrastSourceRevision(c) != projection.sourceRevision)
        fail(QStringLiteral("输入文件自投影计算后已变化，请返回数据与 ROI 重新校验并计算投影。"));
    const auto &validation = projection.validation;
    const auto &energy = validation.sample.axes[2], &foil = validation.calibration.axes[1];
    const bool matched = validation.report.value("psf_energy_band").toString() == "sample-matched";
    const auto band = matched ? Imaging::matchedWeights(foil, energy).first : Imaging::energyWeights(foil);
    auto selected = calibration(c, validation, {}, stop, log, band);
    sample(c, validation, selected, {}, stop, log);
    ImagingResult spatial;
    auto fit = solveSpatial(c, validation, selected, spatial, stop, log);
    selected.clear();
    fit.accepted.clear();
    fit.response.clear();
    fit.rho.resize(0, 0);

    log(QStringLiteral("反演带符号的三维差值"), 80);
    ImageMatrix rhs = std::move(fit.b);
    rhs.setZero();
    for (qsizetype i = 0; i < fit.keys.size(); ++i) {
        cancel(stop);
        const auto &difference = projection.differenceByCrystal.at(fit.keys[i]);
        rhs.noalias() += fit.gains[i] * fit.operators[i].transpose() * difference;
    }
    const auto sensitivity = fit.h.diagonal().eval();
    const double scale = fit.h.cwiseAbs().rowwise().sum().maxCoeff();
    if (!(std::isfinite(scale) && scale > 0.))
        fail(QStringLiteral("无法反演三维差值：空间响应无效。"));
    constexpr double relativeRidge = 1e-6;
    const double ridge = scale * relativeRidge;
    fit.h.diagonal().array() += ridge;
    Eigen::LLT<ImageMatrix> solver(fit.h);
    if (solver.info() != Eigen::Success)
        fail(QStringLiteral("无法分解三维差值的空间响应。"));
    cancel(stop);
    solver.solveInPlace(rhs);
    if (!rhs.allFinite())
        fail(QStringLiteral("三维差值反演产生了无效数值。"));
    const double sensitivityFloor = sensitivity.maxCoeff() * 1e-14;
    for (qsizetype x = 0; x < rhs.rows(); ++x)
        if (sensitivity[x] <= sensitivityFloor)
            rhs.row(x).setConstant(NAN);

    BondContrast out;
    out.axes = std::move(spatial.axes);
    out.values = vector(rhs);
    out.diagnostics = {{"schema", "XRI.bond.contrast.v1"},
                       {"difference", "S(Eb)-S(Ea) before spatial reconstruction"},
                       {"ea_kev", projection.ea},
                       {"eb_kev", projection.eb},
                       {"spectral_processing", projection.spectralProcessing},
                       {"source_revision", projection.sourceRevision},
                       {"solver", "signed ridge least squares"},
                       {"relative_ridge", relativeRidge},
                       {"ridge", ridge},
                       {"crystal_keys", json(fit.keys)},
                       {"relative_gains", json(fit.gains)},
                       {"gain_fit_converged", fit.gainConverged},
                       {"input_validation", validation.report},
                       {"config", c}};
    log(QStringLiteral("三维差值反演完成"), 100);
    return out;
}

BondVolume BondImaging::run(const QJsonObject &supplied, const std::atomic_bool &stop,
                            Imaging::Progress progress) {
    const auto c = Imaging::config(supplied, QDir::currentPath());
    const auto log = [&](const QString &message, int percent) {
        cancel(stop);
        if (progress)
            progress(message, percent);
    };
    log(QStringLiteral("校验样品、薄膜、扫描坐标和 ROI"), 0);
    const auto validation = validate(c, stop);
    const auto &g = validation.sample;
    const auto &energy = g.axes[2], &foil = validation.calibration.axes[1];
    const bool matched = validation.report.value("psf_energy_band").toString() == "sample-matched";
    const auto weights = matched ? Imaging::matchedWeights(foil, energy).first : Imaging::energyWeights(foil);
    auto selected = calibration(c, validation, {}, stop, log, weights);
    sample(c, validation, selected, {}, stop, log);
    ImagingResult spatial;
    auto fit = solveSpatial(c, validation, selected, spatial, stop, log, energy.size());
    const int nx = int(spatial.axes[0].size()), ny = int(g.axes[0].size()), nz = int(g.axes[1].size());
    const int ne = int(energy.size()), nYZ = ny * nz;
    const qsizetype total = qsizetype(nx) * nYZ * ne;
    BondVolume out;
    out.axes = {energy, spatial.axes[0], g.axes[0], g.axes[1]};
    out.intensity.resize(total);
    for (int x = 0; x < nx; ++x)
        out.sensitivity.append(std::sqrt(std::max(0., fit.h(x, x))));
    const double sensitivityFloor =
        *std::max_element(out.sensitivity.cbegin(), out.sensitivity.cend()) * 1e-7;

    // Stack the calibrated detector pixels once. A small batch allows a matrix product
    // even when consecutive frames belong to different energies or spatial cells.
    QVector<qsizetype> offsets;
    for (const auto &item : fit.accepted)
        for (auto pixel : item.roi.pixels)
            offsets.append(qsizetype(pixel.y()) * g.imageSize.width() + pixel.x());
    const qsizetype nPixels = offsets.size();
    const auto output = c["output"].toObject();
    const bool saveBandDiagnostics = output["save_diagnostics"].toBool();
    constexpr int batchSize = 16;
    const qsizetype repeated = g.occupancy.size() - std::count(g.occupancy.cbegin(), g.occupancy.cend(), 1);
    if (nPixels * (qsizetype(nx) + batchSize + repeated) > 32000000)
        fail(QStringLiteral("四维反投影缓冲区过大，请缩小 ROI 或减少重复帧。"));
    ImageMatrix forward(nPixels, nx);
    qsizetype row = 0;
    for (int i = 0; i < fit.accepted.size(); ++i) {
        const auto &op = fit.operators[i];
        forward.middleRows(row, op.rows()) = fit.gains[i] * op;
        row += op.rows();
    }
    std::vector<ImageMatrix> backs(ne, ImageMatrix::Zero(nx, nYZ));
    Eigen::MatrixXd batch(nPixels, batchSize);
    QVector<int> cells;
    QVector<double> norm(ne, 0.);
    QHash<int, Eigen::VectorXd> pooled;
    auto remaining = g.occupancy;
    auto flush = [&] {
        if (cells.isEmpty())
            return;
        Eigen::MatrixXd projected = forward.transpose() * batch.leftCols(cells.size());
        for (int j = 0; j < cells.size(); ++j) {
            int e = cells[j] % ne;
            backs[e].col(cells[j] / ne) = projected.col(j);
            norm[e] += batch.col(j).squaredNorm();
        }
        cells.clear();
    };
    const auto acq = c["acquisition"].toObject();
    {
        ImagingH5 source(g.file);
        source.frames(
            "/entry/data/" + acq["detector"].toString(), c["runtime"].toObject()["chunk_frames"].toInt(),
            stop, [&](int frame, const quint64 *pixels, int, int) {
                int cell = (g.indices[0][frame] * nz + g.indices[1][frame]) * ne + g.indices[2][frame];
                const double factor = 1. / g.dose[cell];
                if (g.occupancy[cell] == 1) {
                    auto column = batch.col(cells.size());
                    for (qsizetype i = 0; i < nPixels; ++i)
                        column[i] = double(pixels[offsets[i]]) * factor;
                } else {
                    auto &counts = pooled[cell];
                    if (!counts.size())
                        counts.setZero(nPixels);
                    for (qsizetype i = 0; i < nPixels; ++i)
                        counts[i] += double(pixels[offsets[i]]);
                    if (--remaining[cell])
                        return;
                    batch.col(cells.size()) = counts * factor;
                    pooled.remove(cell);
                }
                cells.append(cell);
                if (cells.size() == batchSize)
                    flush();
                if (frame % 128 == 0)
                    log(QStringLiteral("按能量累积反投影 %1 / %2").arg(frame + 1).arg(g.monitor.size()),
                        66 + 12 * frame / int(g.monitor.size()));
            });
    }
    flush();
    QJsonArray channels;
    const auto rec = c["reconstruction"].toObject();
    const double bandNorm = fit.b.norm();
    for (int e = 0; e < ne; ++e) {
        log(QStringLiteral("独立求解能量通道 %1 / %2").arg(e + 1).arg(ne), 78 + 21 * e / ne);
        ImagingSolverInfo info;
        ImageMatrix initial =
            bandNorm > 0 ? ImageMatrix(fit.rho * (backs[e].norm() / bandNorm)) : ImageMatrix{};
        auto rho = Imaging::nnls(fit.h, backs[e], rec["nnls_max_iterations"].toInt(),
                                 rec["nnls_tolerance"].toDouble(), stop, &info, initial);
        const double error = norm[e] - 2 * (rho.array() * backs[e].array()).sum() +
                             (rho.array() * (fit.h * rho).array()).sum();
        channels.append(
            QJsonObject{{"incident_kev", energy[e]},
                        {"converged", info.converged},
                        {"iterations", info.iterations},
                        {"projected_kkt_relative", info.kkt},
                        {"relative_l2_residual",
                         norm[e] > 0 ? QJsonValue(std::sqrt(std::max(0., error) / norm[e])) : QJsonValue()}});
        out.converged.append(info.converged);
        std::copy_n(rho.data(), rho.size(), out.intensity.data() + qsizetype(e) * rho.size());
    }
    QStringList crystalLabels;
    for (const auto &item : std::as_const(fit.accepted))
        crystalLabels.append(item.roi.label);
    out.diagnostics = {{"schema", "XRI.bond.v1"},
                       {"config", c},
                       {"input_validation", validation.report},
                       {"psf_band_kev", QJsonArray{matched ? energy.first() : foil.first(),
                                                   matched ? energy.last() : foil.last()}},
                       {"psf_models", spatial.models},
                       {"psf_rejected", spatial.rejected},
                       {"crystal_keys", json(fit.keys)},
                       {"crystal_labels", json(crystalLabels)},
                       {"relative_gains", json(fit.gains)},
                       {"gain_fit_converged", fit.gainConverged},
                       {"gain_fit_iterations", fit.gainIteration},
                       {"x_sensitivity_floor", sensitivityFloor},
                       {"channel_diagnostics", channels}};
    if (saveBandDiagnostics) {
        for (int i = 0; i < fit.accepted.size(); ++i) {
            cancel(stop);
            // Reuse the integrated sample buffers once the spatial fit and channel solves are done.
            auto &op = fit.operators[i];
            auto &predicted = fit.accepted[i].variance;
            auto &residual = fit.accepted[i].measured;
            op *= fit.gains[i];
            predicted.noalias() = op * fit.rho;
            residual = predicted - residual;
            out.bandDiagnostics.append(
                {fit.keys[i], std::move(op), std::move(predicted), std::move(residual)});
        }
    }
    if (output["plots"].toBool()) {
        cancel(stop);
        spatial.primary = fit.rho;
        spatial.sensitivity = out.sensitivity;
        out.diagnosticImage = diagnosticImage(spatial, true);
    }
    log(QStringLiteral("四维求解完成"), 100);
    return out;
}

ImagingResult Imaging::run(const QJsonObject &supplied, const std::atomic_bool &stop, Progress progress) {
    auto c = config(supplied, QDir::currentPath());
    if (progress)
        progress(QStringLiteral("校验扫描计划、坐标、I0 与 ROI"), 0);
    auto validation = validate(c, stop);
    auto directory = c["runtime"].toObject()["output_dir"].toString();
    if (QFileInfo::exists(directory))
        fail("Output directory exists; no overwrite/resume: " + directory);
    if (!QDir().mkpath(directory))
        fail("Cannot create imaging output directory");
    auto log = [&](const QString &message, int percent) {
        cancel(stop);
        QFile file(directory + "/run.log");
        if (!file.open(QIODevice::WriteOnly | QIODevice::Append))
            fail("Cannot append imaging run log");
        auto line = QDateTime::currentDateTime().toString("HH:mm:ss ") + message + "\n";
        file.write(line.toUtf8());
        if (progress)
            progress(message, percent);
    };
    try {
        saveJson(directory + "/config.resolved.json", c);
        saveJson(directory + "/input_validation.json", validation.report);
        QJsonObject sources;
        auto input = c["input"].toObject();
        for (auto it = input.begin(); it != input.end(); ++it) {
            QFileInfo f(it->toString());
            sources[it.key()] = QJsonObject{{"path", f.absoluteFilePath()},
                                            {"size_bytes", double(f.size())},
                                            {"mtime_ms", double(f.lastModified().toMSecsSinceEpoch())}};
        }
        QFile roi(input["roi_table"].toString());
        if (!roi.open(QIODevice::ReadOnly))
            fail("Cannot hash ROI source");
        auto roiSource = sources["roi_table"].toObject();
        roiSource["sha256"] =
            QString(QCryptographicHash::hash(roi.readAll(), QCryptographicHash::Sha256).toHex());
        sources["roi_table"] = roiSource;
        saveJson(
            directory + "/provenance.json",
            {{"versions",
              QJsonObject{
                  {"xri", "0.2.2-native"}, {"qt", qVersion()}, {"eigen", "3.4.0"}, {"hdf5", H5_VERSION}}},
             {"inputs", sources},
             {"started_utc", QDateTime::currentDateTimeUtc().toString(Qt::ISODate)}});
        log("1/6 Extracting calibration pixels and pooled monitor doses", 3);
        auto selected =
            calibration(c, validation, directory, stop, progress,
                        matchedWeights(validation.calibration.axes[1], validation.sample.axes[2]).first);
        log("2/6 Mapping sample YZE and matching energy integrals", 25);
        sample(c, validation, selected, directory, stop, progress);
        log("3/6 Fitting moving elliptical PSFs and affine trajectories", 46);
        ImagingResult result;
        result.directory = directory;
        auto solution = solveSpatial(c, validation, selected, result, stop, log);
        const auto &[accepted, keys, response, operators, gains, history, h, b, rho, nnlsInfo, gainConverged,
                     gainIteration, dx] = solution;
        int nx = int(result.axes[0].size()), ny = int(result.axes[1].size()), nz = int(result.axes[2].size());
        const auto &x = result.axes[0];
        const auto &originalX = validation.calibration.axes[0];
        const auto acq = c["acquisition"].toObject(), rec = c["reconstruction"].toObject();
        const double pitch = acq["detector_pixel_pitch_mm"].toDouble();
        result.nnls = rho;
        result.primary = rho;
        QString method = "joint NNLS";
        std::array<double, 3> spacing{dx, (result.axes[1].last() - result.axes[1].first()) / (ny - 1),
                                      (result.axes[2].last() - result.axes[2].first()) / (nz - 1)};
        auto tvSettings = c["tv"].toObject();
        ImagingSolverInfo tvInfo;
        double beta = 0;
        if (tvSettings["enabled"].toBool()) {
            log("5/6 Fixed-gain nonnegative TV", 75);
            if (tvSettings["lambda"].isNull()) {
                auto absolute = vector(b.cwiseAbs());
                beta = tvSettings["relative_fraction"].toDouble() * percentile(absolute, 95) *
                       std::min({spacing[0], spacing[1], spacing[2]});
            } else
                beta = tvSettings["lambda"].toDouble();
            log(QString("TV lambda = %1").arg(beta, 0, 'g', 12), 75);
            result.tv = tv(h, b, ny, nz, spacing, beta, tvSettings["max_iterations"].toInt(),
                           tvSettings["tolerance"].toDouble(), rho, stop, &tvInfo, log);
            if (tvInfo.converged) {
                result.primary = result.tv;
                method = "joint NNLS gains + fixed-gain TV";
            } else
                log("TV did not converge: primary output uses NNLS; TV iterate saved separately", 95);
        } else
            log("5/6 TV disabled", 95);
        Eigen::SelfAdjointEigenSolver<ImageMatrix> eigen(h);
        if (eigen.info() != Eigen::Success)
            fail("Gram eigensolver failed");
        for (int i = 0; i < nx; ++i) {
            result.sensitivity.append(std::sqrt(std::max(0., h(i, i))));
            result.eigenvalues.append(eigen.eigenvalues()[i]);
        }
        int rank = 0;
        for (double e : result.eigenvalues)
            if (e > std::max(0., result.eigenvalues.last()) * 1e-14)
                ++rank;
        double normData = 0, nnlsError = 0, tvError = 0;
        QJsonObject residuals, gainMap;
        QStringList atBounds;
        auto bounds = rec["gain_bounds"].toArray();
        for (int i = 0; i < accepted.size(); ++i) {
            const auto &data = accepted[i].measured;
            double norm = data.squaredNorm();
            normData += norm;
            ImageMatrix prediction = gains[i] * operators[i] * rho;
            nnlsError += (prediction - data).squaredNorm();
            if (result.tv.size()) {
                ImageMatrix pred = gains[i] * operators[i] * result.tv;
                tvError += (pred - data).squaredNorm();
            }
            ImageMatrix pred = gains[i] * operators[i] * result.primary;
            residuals[keys[i]] = norm > 0 ? QJsonValue((pred - data).norm() / std::sqrt(norm)) : QJsonValue();
            gainMap[keys[i]] = gains[i];
            for (auto limit : bounds)
                if (std::abs(gains[i] - limit.toDouble()) <= 1e-8 + 1e-5 * std::abs(limit.toDouble())) {
                    atBounds.append(keys[i]);
                    break;
                }
        }
        QJsonArray shape{nx, ny, nz}, ranges;
        for (const auto &axis : result.axes)
            ranges.append(QJsonArray{axis.first(), axis.last()});
        result.summary = {
            {"status", "provisional"},
            {"solver", method},
            {"shape_xyz", shape},
            {"spacing_xyz_mm", QJsonArray{spacing[0], spacing[1], spacing[2]}},
            {"ranges_xyz_mm", ranges},
            {"x_grid_policy", "full calibrated span; centered uniform centers; no sensitivity cropping"},
            {"x_step_source", c["grid"].toObject()["x_step_mm"].isNull()
                                  ? "median calibrated sample mm per detector pixel"
                                  : "explicit sample-space step"},
            {"detector_pixel_pitch_mm", pitch},
            {"foil_thickness_mm", acq["foil_thickness_mm"]},
            {"calibration_x_positions", originalX.size()},
            {"energy_range_kev",
             QJsonArray{validation.sample.axes[2].first(), validation.sample.axes[2].last()}},
            {"energy_points", validation.sample.axes[2].size()},
            {"input_roi_count", selected.size()},
            {"accepted_roi_count", accepted.size()},
            {"roi_keys", json(keys)},
            {"rejected_psf_models", result.rejected},
            {"gains_by_roi", gainMap},
            {"gain_anchor", rec["fit_gains"].toBool() ? rec["anchor_roi"] : QJsonValue()},
            {"gain_fit_enabled", rec["fit_gains"]},
            {"gain_fit_converged", gainConverged},
            {"gain_fit_iterations", gainIteration},
            {"gain_relative_change", history.isEmpty() ? 0 : history.last()},
            {"nnls_converged", nnlsInfo.converged},
            {"nnls_final_iterations", nnlsInfo.iterations},
            {"nnls_projected_kkt_relative", nnlsInfo.kkt},
            {"gains_at_bounds", json(atBounds)},
            {"tv_enabled", tvSettings["enabled"]},
            {"tv_lambda", result.tv.size() ? QJsonValue(beta) : QJsonValue()},
            {"tv_converged", result.tv.size() ? QJsonValue(tvInfo.converged) : QJsonValue()},
            {"tv_iterations", tvInfo.iterations},
            {"tv_primal_relative_change", result.tv.size() ? QJsonValue(tvInfo.primal) : QJsonValue()},
            {"tv_dual_relative_change", result.tv.size() ? QJsonValue(tvInfo.dual) : QJsonValue()},
            {"nnls_relative_l2_residual",
             normData > 0 ? QJsonValue(std::sqrt(nnlsError / normData)) : QJsonValue()},
            {"tv_relative_l2_residual",
             result.tv.size() && normData > 0 ? QJsonValue(std::sqrt(tvError / normData)) : QJsonValue()},
            {"numerical_response_rank", rank},
            {"x_cells", nx},
            {"gram_eigenvalue_relative_cutoff", 1e-14},
            {"roi_relative_l2_residual", residuals},
            {"sample_background_subtracted", false},
            {"fitted_reference_offset_in_forward_model", false},
            {"postprocessing_crop", false},
            {"postprocessing_smooth", false},
            {"symmetry_constraint", false},
            {"limitations", QJsonArray{"Grid spacing is not measured spatial resolution; foil positions are "
                                       "observations, not X cell count.",
                                       "Moving Gaussian PSF and affine centroid with constant "
                                       "amplitude/width; finite foil thickness is not deconvolved.",
                                       "Full X includes weak or unobservable depths; inspect sensitivity and "
                                       "Gram eigenvalues; zeros do not prove empty material.",
                                       "TV does not add measurements or improve measured resolution.",
                                       "Relative Kapton-equivalent density; no sample background, absorption "
                                       "or composition correction.",
                                       "Unweighted least squares; counting variances exclude monitor, model, "
                                       "gain and background uncertainty."}}};
        if (!gainConverged)
            log("Gain/NNLS fit not converged: outputs remain provisional", 96);
        log("6/6 Saving volumes, XYZ axes, sensitivity, responses and residuals", 96);
        QJsonObject metadata{{"detector_pixel_pitch_mm", pitch},
                             {"sample_x_spacing_mm", dx},
                             {"foil_thickness_mm", acq["foil_thickness_mm"]},
                             {"crystal_gains", json(gains)},
                             {"crystal_keys", keys.join(',')},
                             {"provisional", true},
                             {"gain_fit_converged", gainConverged},
                             {"nnls_converged", nnlsInfo.converged},
                             {"tv_requested", tvSettings["enabled"]},
                             {"no_geometric_postprocessing", true}};
        if (result.tv.size()) {
            metadata["tv_converged"] = tvInfo.converged;
            metadata["tv_lambda"] = beta;
        }
        exportVolume(directory + "/reconstruction.h5", result.primary, result.axes, directory, method,
                     metadata, result.sensitivity, result.eigenvalues);
        exportVolume(directory + "/reconstruction_nnls.h5", rho, result.axes, directory, "joint NNLS",
                     metadata, result.sensitivity, result.eigenvalues);
        if (result.tv.size())
            exportVolume(directory + "/reconstruction_tv.h5", result.tv, result.axes, directory,
                         "fixed-gain TV", metadata, result.sensitivity, result.eigenvalues);
        {
            ImagingH5 diagnosticsFile(directory + "/solver_diagnostics.h5", true);
            writeVector(diagnosticsFile, "gain_relative_change_history", history);
            writeVector(diagnosticsFile, "x_mm", x);
            writeVector(diagnosticsFile, "sensitivity_x", result.sensitivity);
            writeVector(diagnosticsFile, "gram_eigenvalues", result.eigenvalues);
            diagnosticsFile.writeStrings("roi_keys", keys);
            writeVector(diagnosticsFile, "gains", gains);
            if (c["output"].toObject()["save_diagnostics"].toBool())
                for (int i = 0; i < accepted.size(); ++i) {
                    cancel(stop);
                    auto base = "rois/" + keys[i];
                    ImageMatrix op = gains[i] * operators[i];
                    writeMatrix(diagnosticsFile, base + "/response_px", response[i]);
                    writeMatrix(diagnosticsFile, base + "/operator_px", op);
                    auto volume = [&](const QString &name, const ImageMatrix &m) {
                        diagnosticsFile.write(base + "/" + name,
                                              {quint64(m.rows()), quint64(ny), quint64(nz)}, m.data());
                    };
                    volume("measured_pyz", accepted[i].measured);
                    ImageMatrix pred = op * rho;
                    volume("predicted_nnls_pyz", pred);
                    volume("residual_nnls_pyz", pred - accepted[i].measured);
                    if (result.tv.size()) {
                        pred.noalias() = op * result.tv;
                        volume("predicted_tv_pyz", pred);
                        volume("residual_tv_pyz", pred - accepted[i].measured);
                    }
                    QVector<quint64> pixels;
                    for (auto point : accepted[i].roi.pixels) {
                        pixels.append(point.y());
                        pixels.append(point.x());
                    }
                    diagnosticsFile.writeUInt(
                        base + "/pixels_rc", {quint64(accepted[i].roi.pixels.size()), 2}, pixels.constData());
                }
        }
        if (c["output"].toObject()["plots"].toBool() &&
            !diagnosticImage(result).save(directory + "/diagnostics.png"))
            fail("Cannot save imaging diagnostics PNG");
        if (c["output"].toObject()["save_avizo"].toBool()) {
            log(QStringLiteral("正在保存 Avizo 体数据 reconstruction.am"), 99);
            exportAvizo(directory + "/reconstruction.am", result, stop);
        }
        saveJson(directory + "/summary.json", result.summary);
        log("Completed: " + directory + "/reconstruction.h5", 100);
        return result;
    } catch (const std::exception &e) {
        QFile file(directory + "/run.log");
        if (file.open(QIODevice::WriteOnly | QIODevice::Append))
            file.write((QString("FAILED: ") + QString::fromUtf8(e.what()) + "\n").toUtf8());
        saveJson(directory + "/run_status.json",
                 {{"completed", false}, {"cancelled", stop.load()}, {"error", QString::fromUtf8(e.what())}});
        throw;
    }
}
