#include "BondVolume.h"
#include <QCoreApplication>
#include "ImagingIO.h"
#include <Eigen/QR>
#include <QDir>
#include <QFileInfo>
#include <QJsonDocument>
#include <QTemporaryDir>
#include <algorithm>
#include <cmath>

using namespace ImageUtil;
namespace {
constexpr qsizetype maxSamples = 16000000;
const QString fieldPath = QStringLiteral("/entry/bond/intensity");
const QString diagnosticsPath = QStringLiteral("/entry/bond/diagnostics_json");
const QStringList axisPaths{"/entry/bond/incident_energy_kev", "/entry/bond/x_mm", "/entry/bond/y_mm",
                            "/entry/bond/z_mm"};

qsizetype checkedSize(const QVector<quint64> &shape) {
    qsizetype count = 1;
    if (shape.size() != 4)
        fail(QCoreApplication::translate("BondVolume", "结果数组的维数与格式不一致。"));
    for (auto n : shape) {
        if (n < 2 || n > quint64(maxSamples / count))
            fail(QCoreApplication::translate("BondVolume", "结果尺寸无效或超过 1600 万采样点。"));
        count *= qsizetype(n);
    }
    return count;
}

void checkAxis(const QVector<double> &axis) {
    if (axis.size() < 2)
        fail(QCoreApplication::translate("BondVolume", "各轴至少需要两个坐标。"));
    for (qsizetype i = 0; i < axis.size(); ++i)
        if (!std::isfinite(axis[i]) || (i && axis[i] <= axis[i - 1]))
            fail(QCoreApplication::translate("BondVolume", "坐标轴必须有限且严格递增。"));
}

struct Bracket {
    int first, last;
    double fraction;
};
Bracket bracket(const QVector<double> &axis, double at) {
    constexpr double epsilon = 1e-9;
    if (!std::isfinite(at) || at < axis.first() - epsilon || at > axis.last() + epsilon)
        fail(QCoreApplication::translate("BondVolume", "指定能量超出实测范围。"));
    int hi = int(std::lower_bound(axis.cbegin(), axis.cend(), at) - axis.cbegin());
    if (hi == axis.size())
        hi--;
    if (hi == 0 || std::abs(axis[hi] - at) <= epsilon)
        return {hi, hi, 0};
    if (std::abs(axis[hi - 1] - at) <= epsilon)
        return {hi - 1, hi - 1, 0};
    return {hi - 1, hi, (at - axis[hi - 1]) / (axis[hi] - axis[hi - 1])};
}

double valueAt(const double *data, qsizetype stride, Bracket b) {
    const double first = data[b.first * stride];
    return b.first == b.last ? first : first + b.fraction * (data[b.last * stride] - first);
}

struct SpectralIntegral {
    QVector<double> energy;
    QVector<Bracket> points;

    SpectralIntegral(const QVector<double> &axis, QPair<double, double> range) {
        const auto [lo, hi] = range;
        if (!(lo < hi))
            fail(QCoreApplication::translate("BondVolume", "积分区间必须递增。"));
        energy.append(lo);
        for (double e : axis)
            if (lo < e && e < hi)
                energy.append(e);
        energy.append(hi);
        for (double e : energy)
            points.append(bracket(axis, e));
    }

    bool coveredBy(const QVector<bool> &usable) const {
        return std::all_of(points.cbegin(), points.cend(), [&](Bracket b) {
            return usable[b.first] && usable[b.last];
        });
    }

    double area(const double *data, qsizetype stride = 1) const {
        double sum = 0, prev = valueAt(data, stride, points.first());
        for (qsizetype i = 1; i < points.size(); ++i) {
            const double current = valueAt(data, stride, points[i]);
            sum += .5 * (prev + current) * (energy[i] - energy[i - 1]);
            prev = current;
        }
        return sum;
    }
};

bool supported(const BondVolume &v, qsizetype voxel) {
    const auto yz = v.axes[2].size() * v.axes[3].size();
    return v.sensitivity[voxel / yz] > v.diagnostics.value("x_sensitivity_floor").toDouble();
}
} // namespace

qsizetype BondVolume::voxelCount() const {
    return axes.size() == 4 ? axes[1].size() * axes[2].size() * axes[3].size() : 0;
}

qsizetype BondVolume::voxel(int x, int y, int z) const {
    return (qsizetype(x) * axes[2].size() + y) * axes[3].size() + z;
}

QString BondImaging::diagnosticPngPath(const QString &path) {
    const QFileInfo result(path);
    return result.dir().filePath(result.completeBaseName() + ".diagnostics.png");
}

void BondImaging::save(const QString &path, const BondVolume &volume, const std::atomic_bool &stop) {
    cancel(stop);
    QVector<quint64> shape;
    for (const auto &axis : volume.axes) {
        checkAxis(axis);
        shape.append(quint64(axis.size()));
    }
    if (checkedSize(shape) != volume.intensity.size() || volume.sensitivity.size() != qsizetype(shape[1]) ||
        volume.converged.size() != qsizetype(shape[0]))
        fail(QCoreApplication::translate("BondVolume", "四维结果与坐标或诊断尺寸不一致。"));
    const QFileInfo target(path);
    if (target.exists())
        fail(QCoreApplication::translate("BondVolume", "结果文件已存在，请选择新文件名。"));
    const bool plots = !volume.diagnosticImage.isNull();
    const auto png = diagnosticPngPath(target.absoluteFilePath());
    if (plots && QFileInfo::exists(png))
        fail(QCoreApplication::translate("BondVolume", "诊断 PNG 已存在，请选择新的结果文件名。"));
    QTemporaryDir temp(QDir(target.absolutePath()).filePath(".xri-bond-XXXXXX"));
    if (!temp.isValid())
        fail(QCoreApplication::translate("BondVolume", "无法在输出目录建立临时文件。"));
    const auto staged = temp.filePath("volume.h5");
    const auto stagedPng = temp.filePath("diagnostics.png");
    if (plots && !volume.diagnosticImage.save(stagedPng))
        fail(QCoreApplication::translate("BondVolume", "无法写入诊断 PNG。"));
    {
        ImagingH5 hdf(staged, true);
        auto axis = [&](const QString &name, const QVector<double> &data) {
            hdf.write(name, {quint64(data.size())}, data.constData(), false);
        };
        hdf.write(fieldPath, shape, volume.intensity.constData());
        for (int i = 0; i < axisPaths.size(); ++i)
            axis(axisPaths[i], volume.axes[i]);
        hdf.scales(fieldPath, axisPaths);
        axis("/entry/bond/x_sensitivity", volume.sensitivity);
        QVector<quint64> valid;
        for (bool good : volume.converged)
            valid.append(good ? 1 : 0);
        hdf.writeUInt("/entry/bond/channel_converged", {quint64(valid.size())}, valid.constData(), false);
        hdf.attributes(fieldPath, {{"axis_order", "E,X,Y,Z"},
                                   {"units", "relative intensity"},
                                   {"energy_kind", "incident"},
                                   {"schema", "XRI.bond.v1"}});
        hdf.writeText(diagnosticsPath,
                      QString::fromUtf8(QJsonDocument(volume.diagnostics).toJson(QJsonDocument::Compact)));
        for (const auto &item : volume.bandDiagnostics) {
            cancel(stop);
            const auto base = "/entry/bond/band_diagnostics/" + item.crystal;
            const auto rows = quint64(item.forward.rows());
            hdf.write(base + "/operator_px", {rows, shape[1]}, item.forward.data());
            hdf.write(base + "/predicted_pyz", {rows, shape[2], shape[3]}, item.predicted.data());
            hdf.write(base + "/residual_pyz", {rows, shape[2], shape[3]}, item.residual.data());
        }
    }
    cancel(stop);
    if (plots && !QFile::rename(stagedPng, png))
        fail(QCoreApplication::translate("BondVolume", "无法提交诊断 PNG。"));
    if (!QFile::rename(staged, target.absoluteFilePath())) {
        if (plots)
            QFile::remove(png);
        fail(QCoreApplication::translate("BondVolume", "无法提交四维结果文件。"));
    }
}

BondVolume BondImaging::load(const QString &path) {
    ImagingH5 hdf(path);
    if (!hdf.contains(fieldPath) && hdf.contains("/entry/bond/contrast"))
        fail(QCoreApplication::translate("BondVolume", "此文件是旧版三维差值结果。请在数据浏览器中查看，或选择四维体素 HDF5。"));
    const auto shape = HdfFile::inspect(path, fieldPath).dims;
    checkedSize(shape);
    BondVolume out;
    auto readVector = [&](const QString &name, quint64 expected) {
        if (HdfFile::inspect(path, name).dims != QVector<quint64>{expected})
            fail(QCoreApplication::translate("BondVolume", "四维结果坐标或诊断尺寸不一致。"));
        return hdf.read(name);
    };
    for (int i = 0; i < axisPaths.size(); ++i) {
        auto axis = readVector(axisPaths[i], shape[i]);
        checkAxis(axis);
        out.axes.append(std::move(axis));
    }
    out.sensitivity = readVector("/entry/bond/x_sensitivity", shape[1]);
    for (double value : out.sensitivity)
        if (!std::isfinite(value) || value < 0)
            fail(QCoreApplication::translate("BondVolume", "无效的 X 灵敏度。"));
    for (double value : readVector("/entry/bond/channel_converged", shape[0])) {
        if (value != 0 && value != 1)
            fail(QCoreApplication::translate("BondVolume", "无效的通道收敛标记。"));
        out.converged.append(value == 1);
    }
    out.diagnostics = QJsonDocument::fromJson(hdf.text(diagnosticsPath).toUtf8()).object();
    if (out.diagnostics.value("schema").toString() != "XRI.bond.v1")
        fail(QCoreApplication::translate("BondVolume", "无法识别四维结果版本。"));
    out.intensity = hdf.read(fieldPath);
    return out;
}

QVector<double> BondImaging::process(const BondVolume &v, const QVector<double> &energy,
                                     const BondSpectralSettings &settings, const std::atomic_bool &stop) {
    checkAxis(energy);
    if (energy.size() != v.axes[0].size())
        fail(QCoreApplication::translate("BondVolume", "处理能量轴与四维结果不一致。"));
    auto checkWindow = [&](QPair<double, double> range) {
        if (!(range.first < range.second))
            fail(QCoreApplication::translate("BondVolume", "拟合或归一化能区必须递增。"));
        bracket(energy, range.first);
        bracket(energy, range.second);
    };
    for (auto range : settings.background)
        checkWindow(range);
    const auto *compton = settings.compton ? &*settings.compton : nullptr;
    auto usable = v.converged;
    std::optional<SpectralIntegral> match, normalize;
    auto prepareIntegral = [&](QPair<double, double> range) {
        SpectralIntegral integral(energy, range);
        if (!integral.coveredBy(usable))
            fail(QCoreApplication::translate("BondVolume", "积分区间包含未收敛通道或超出参考轮廓覆盖。"));
        return integral;
    };
    double referenceArea = 0;
    if (compton) {
        if (compton->reference.size() != energy.size() || compton->valence.size() != energy.size() ||
            compton->correction.size() != energy.size())
            fail(QCoreApplication::translate("BondVolume", "康普顿参考与重建能量轴不一致。"));
        for (qsizetype e = 0; e < energy.size(); ++e) {
            if (!(std::isfinite(compton->correction[e]) && compton->correction[e] > 0))
                fail(QCoreApplication::translate("BondVolume", "康普顿截面修正无效。"));
            usable[e] = usable[e] && std::isfinite(compton->valence[e]);
        }
        match = prepareIntegral(compton->matchWindow);
        referenceArea = match->area(compton->reference.constData());
        if (!(std::isfinite(referenceArea) && referenceArea > 0))
            fail(QCoreApplication::translate("BondVolume", "幅度匹配区间需要完整覆盖和正的宽谱积分。"));
    }
    if (settings.normalization)
        normalize = prepareIntegral(*settings.normalization);
    QVector<int> fitPoints;
    for (int e = 0; e < energy.size(); ++e)
        if (usable[e] &&
            std::any_of(settings.background.cbegin(), settings.background.cend(),
                        [&](auto w) { return energy[e] >= w.first - 1e-9 && energy[e] <= w.second + 1e-9; }))
            fitPoints.append(e);
    if (!settings.background.isEmpty() && fitPoints.size() < 2)
        fail(QCoreApplication::translate("BondVolume", "背底能区至少需要两个已收敛且有参考覆盖的能量点。"));
    const auto n = v.voxelCount();
    QVector<double> out(v.intensity.size(), NAN), row(energy.size());
    const auto *input = v.intensity.constData();
    auto *result = out.data();
    for (qsizetype i = 0; i < n; ++i) {
        if (i % 128 == 0)
            cancel(stop);
        if (!supported(v, i))
            continue;
        for (int e = 0; e < energy.size(); ++e)
            row[e] = usable[e] ? input[e * n + i] : NAN;
        if (compton) {
            for (qsizetype e = 0; e < energy.size(); ++e)
                row[e] *= compton->correction[e];
            const double scale = match->area(row.constData()) / referenceArea;
            if (!(std::isfinite(scale) && scale > 0))
                continue;
            for (qsizetype e = 0; e < energy.size(); ++e)
                row[e] -= scale * compton->valence[e];
        }
        if (!fitPoints.isEmpty()) {
            double meanX = 0, meanY = 0;
            int count = 0;
            for (int e : fitPoints)
                if (std::isfinite(row[e])) {
                    meanX += energy[e];
                    meanY += row[e];
                    ++count;
                }
            if (count < 2)
                continue;
            meanX /= count;
            meanY /= count;
            double xx = 0, xy = 0;
            for (int e : fitPoints)
                if (std::isfinite(row[e])) {
                    const double dx = energy[e] - meanX;
                    xx += dx * dx;
                    xy += dx * (row[e] - meanY);
                }
            for (int e = 0; e < energy.size(); ++e)
                row[e] -= meanY + (xy / xx) * (energy[e] - meanX);
        }
        // Return to the reconstruction's intensity scale before any optional area normalization.
        if (compton)
            for (qsizetype e = 0; e < energy.size(); ++e)
                row[e] /= compton->correction[e];
        const double area = normalize ? normalize->area(row.constData()) : 1.;
        if (!(std::isfinite(area) && area > 0))
            continue;
        for (int e = 0; e < energy.size(); ++e)
            result[e * n + i] = row[e] / area;
    }
    cancel(stop);
    if (std::none_of(out.cbegin(), out.cend(), [](double value) { return std::isfinite(value); }))
        fail(QCoreApplication::translate("BondVolume", "所选参数未生成有效体素谱，请检查能区、参考覆盖和积分面积。"));
    return out;
}

QVector<double> BondImaging::energySlice(const BondVolume &v, const QVector<double> &data,
                                         const QVector<double> &energy, double at, bool requireConvergence) {
    checkAxis(energy);
    if (energy.size() != v.axes[0].size() || data.size() != v.intensity.size())
        fail(QCoreApplication::translate("BondVolume", "能量轴与四维结果不一致。"));
    const auto b = bracket(energy, at);
    if (requireConvergence && (!v.converged[b.first] || !v.converged[b.last]))
        fail(QCoreApplication::translate("BondVolume", "指定能量所需的通道未收敛，请选择已收敛通道。"));
    const auto n = v.voxelCount();
    QVector<double> out(n, NAN);
    for (qsizetype i = 0; i < n; ++i)
        if (!requireConvergence || supported(v, i))
            out[i] = valueAt(data.constData() + i, n, b);
    return out;
}

QVector<double> BondImaging::energyIntegral(const BondVolume &v, const QVector<double> &data,
                                           QPair<double, double> boundsKev) {
    const auto n = v.voxelCount();
    if (!n || v.converged.size() != v.axes[0].size() || v.sensitivity.size() != v.axes[1].size() ||
        data.size() != n * v.axes[0].size())
        fail(QCoreApplication::translate("BondVolume", "能量轴与四维结果不一致。"));
    const auto &energy = v.axes[0];
    checkAxis(energy);
    const SpectralIntegral integral(energy, boundsKev);
    if (!integral.coveredBy(v.converged))
        fail(QCoreApplication::translate("BondVolume", "积分区间包含未收敛通道，请选择连续的已收敛能区。"));

    QVector<double> out(n, NAN);
    const auto *input = data.constData();
    for (qsizetype i = 0; i < n; ++i) {
        if (!supported(v, i))
            continue;
        const double area = integral.area(input + i, n);
        if (std::isfinite(area))
            out[i] = area;
    }
    return out;
}

QVector<BondTriangle> BondImaging::surface(const BondVolume &v, const QVector<double> &field, double level,
                                           const std::atomic_bool &stop) {
    if (v.axes.size() != 4)
        fail(QCoreApplication::translate("BondVolume", "等值面输入无效。"));
    return Isosurface::extract(v.axes.sliced(1), field.constData(), field.size(), level, stop);
}

BondReferenceSlice BondImaging::fitReferenceSlice(
    const BondVolume &v, int z, double analyzerKev, const BondReferenceSpectrum &matrix,
    const BondReferenceSpectrum &target, std::optional<QPair<double, double>> window, double threshold,
    const std::atomic_bool &stop) {
    if (v.axes.size() != 4 || z < 0 || z >= v.axes[3].size() || !(analyzerKev > 0) ||
        !std::isfinite(analyzerKev) || !std::isfinite(threshold) || threshold < 0)
        fail(QCoreApplication::translate("BondVolume", "请检查 Z 切片、分析器能量和系数阈值。"));
    for (const auto *ref : {&matrix, &target}) {
        checkAxis(ref->energyEv);
        if (ref->intensity.size() != ref->energyEv.size())
            fail(QCoreApplication::translate("BondVolume", "参考谱的能量与强度长度不一致。"));
    }
    BondReferenceSlice out;
    out.z = z;
    out.threshold = threshold;
    QVector<qsizetype> channels;
    const double lo = std::max(matrix.energyEv.first(), target.energyEv.first());
    const double hi = std::min(matrix.energyEv.last(), target.energyEv.last());
    if (window && (!(window->first < window->second) || window->first < lo || window->second > hi))
        fail(QCoreApplication::translate("BondVolume", "拟合能区须位于两条参考谱共同覆盖的能量范围内。"));
    const auto &incident = v.axes[0];
    for (qsizetype i = 0; i < incident.size(); ++i) {
        const double e = 1000 * (incident[i] - analyzerKev);
        if (!v.converged[i] || e < lo - 1e-9 || e > hi + 1e-9 ||
            (window && (e < window->first - 1e-9 || e > window->second + 1e-9)))
            continue;
        const double a = valueAt(matrix.intensity.constData(), 1, bracket(matrix.energyEv, e));
        const double b = valueAt(target.intensity.constData(), 1, bracket(target.energyEv, e));
        if (!std::isfinite(a) || !std::isfinite(b))
            fail(QCoreApplication::translate("BondVolume", "拟合能区的参考谱含缺失强度。"));
        channels.append(i);
        out.energyEv.append(e);
        out.matrixReference.append(a);
        out.targetReference.append(b);
    }
    if (channels.size() < 5)
        fail(QCoreApplication::translate("BondVolume", "两条参考谱与样品的共同拟合能区内至少需要五个已收敛通道。"));
    for (auto *ref : {&out.matrixReference, &out.targetReference}) {
        double scale = 0.;
        for (double value : std::as_const(*ref))
            scale = std::max(scale, std::abs(value));
        if (!(scale > 0.))
            fail(QCoreApplication::translate("BondVolume", "参考谱在拟合能区内没有有效强度。"));
        for (double &value : *ref)
            value /= scale;
    }

    const auto n = channels.size();
    const double span = out.energyEv.last() - out.energyEv.first();
    const double center = out.energyEv.first() + .5 * span;
    Eigen::MatrixXd design(n, 4);
    for (qsizetype i = 0; i < n; ++i)
        design.row(i) << 1., (out.energyEv[i] - center) / span, out.matrixReference[i], out.targetReference[i];
    const Eigen::MatrixXd matrixOnly = design.leftCols(3);
    const auto matrixFit = matrixOnly.colPivHouseholderQr();
    auto fullFit = design.colPivHouseholderQr();
    fullFit.setThreshold(1e-8);
    if (fullFit.rank() < 4)
        fail(QCoreApplication::translate("BondVolume", "参考谱和线性背底在所选能区内无法独立区分。"));

    const auto nx = v.axes[1].size(), ny = v.axes[2].size(), voxels = v.voxelCount();
    out.total.fill(NAN, nx * ny);
    out.residual.fill(NAN, nx * ny);
    out.coefficient.fill(NAN, nx * ny);
    out.positiveMean.fill(0., n);
    out.negativeMean.fill(0., n);
    Eigen::VectorXd row(n), corrected(n);
    qsizetype fitted = 0;
    for (qsizetype y = 0; y < ny; ++y) {
        for (qsizetype x = 0; x < nx; ++x) {
            const auto pixel = y * nx + x;
            if (pixel % 128 == 0)
                cancel(stop);
            const auto voxel = v.voxel(int(x), int(y), z);
            if (!supported(v, voxel))
                continue;
            for (qsizetype i = 0; i < n; ++i)
                row[i] = v.intensity[channels[i] * voxels + voxel];
            if (!row.allFinite())
                continue;
            const Eigen::Vector3d backgroundFit = matrixFit.solve(row);
            const Eigen::Vector4d coefficients = fullFit.solve(row);
            ++fitted;
            out.total[pixel] = row.sum();
            out.residual[pixel] = (row - matrixOnly * backgroundFit).squaredNorm();
            out.coefficient[pixel] = coefficients[3];
            if (std::abs(coefficients[3]) <= threshold)
                continue;
            corrected = row - design.leftCols(2) * coefficients.head<2>();
            const double peak = corrected.cwiseAbs().maxCoeff();
            if (!(peak > 0.))
                continue;
            auto &mean = coefficients[3] > 0. ? out.positiveMean : out.negativeMean;
            auto &count = coefficients[3] > 0. ? out.positiveCount : out.negativeCount;
            ++count;
            for (qsizetype i = 0; i < n; ++i)
                mean[i] += corrected[i] / peak;
        }
    }
    if (!fitted)
        fail(QCoreApplication::translate("BondVolume", "当前切片没有覆盖完整拟合能区的有效体素。"));
    for (auto [curve, count] : {std::pair{&out.positiveMean, out.positiveCount},
                               std::pair{&out.negativeMean, out.negativeCount}}) {
        if (count == 0)
            curve->clear();
        else
            for (double &value : *curve)
                value /= count;
    }
    cancel(stop);
    return out;
}
