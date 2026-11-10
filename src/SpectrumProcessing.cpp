#include "SpectrumProcessing.h"
#include <QCoreApplication>
#include "ImagingIO.h"
#include <Eigen/Cholesky>
#include <QDir>
#include <QFile>
#include <QFileInfo>
#include <QJsonDocument>
#include <QTextStream>
#include <QUuid>
#include <algorithm>
#include <numeric>

namespace {
using ImageUtil::fail;
constexpr auto group = "/entry/spectrum/";
constexpr double energyTolerance = 1e-7;
constexpr qsizetype maxPoints = 200000, maxCrystals = 256, maxValues = 20000000;
void checkCancelled(const std::atomic_bool &stop) {
    if (stop.load())
        fail(QCoreApplication::translate("SpectrumProcessing", "数据处理已取消。"));
}
ProcessedSpectrum readSpectrum(const QString &file, const SpectrumHdfLayout &layout) {
    auto axis = HdfFile::inspect(file, layout.energy);
    if (!axis.numeric || axis.dims.size() != 1 || axis.dims[0] < 2 || axis.dims[0] > maxPoints)
        fail(QCoreApplication::translate("SpectrumProcessing", "能量转移数据须为含 2 至 200000 点的一维数值数组。"));
    ImagingH5 hdf(file);
    auto labels = hdf.strings(layout.labels, maxCrystals);
    const auto n = axis.dims[0], rows = quint64(labels.size());
    const auto shape = layout.energyFirst ? QVector<quint64>{n, rows} : QVector<quint64>{rows, n};
    auto intensity = HdfFile::inspect(file, layout.intensity);
    if (!intensity.numeric || intensity.dims != shape || n * rows > maxValues)
        fail(QCoreApplication::translate("SpectrumProcessing", "强度矩阵的维度须与所选排列、能量轴和晶体标签数量匹配。"));
    ProcessedSpectrum out;
    out.sources = {file};
    auto unitPath = QString(group) + "intensity_unit";
    out.intensityUnit = hdf.contains(unitPath) ? hdf.text(unitPath) : layout.externalIntensityUnit;
    out.intensityUnit = out.intensityUnit.trimmed();
    if (out.intensityUnit.isEmpty())
        fail(QCoreApplication::translate("SpectrumProcessing", "请指定强度单位。"));
    out.energyEv = hdf.read(layout.energy);
    for (double &value : out.energyEv) {
        if (layout.energyInKev)
            value *= 1000;
        if (!std::isfinite(value))
            fail(QCoreApplication::translate("SpectrumProcessing", "能量轴含非有限值。"));
    }
    bool descending = out.energyEv.last() < out.energyEv.first();
    if (descending)
        std::reverse(out.energyEv.begin(), out.energyEv.end());
    for (qsizetype i = 1; i < out.energyEv.size(); ++i)
        if (out.energyEv[i] <= out.energyEv[i - 1])
            fail(QCoreApplication::translate("SpectrumProcessing", "能量轴须严格单调，且同一文件内的能量点不可重复。"));
    const auto values = hdf.read(layout.intensity);
    for (qsizetype r = 0; r < labels.size(); ++r) {
        const auto &label = labels[r];
        if (label.isEmpty() || label.size() > 256 || out.crystals.contains(label))
            fail(QCoreApplication::translate("SpectrumProcessing", "晶体标签须非空且唯一：%1").arg(label));
        QVector<double> row(n);
        for (qsizetype i = 0; i < row.size(); ++i) {
            auto source = descending ? row.size() - 1 - i : i;
            double value = values[layout.energyFirst ? source * labels.size() + r : r * row.size() + source];
            if (std::isinf(value))
                fail(QCoreApplication::translate("SpectrumProcessing", "晶体 %1 的强度含无穷值。").arg(label));
            row[i] = value;
        }
        out.crystals.insert(label, std::move(row));
    }
    QVector<double> geometry[3];
    const QStringList names{"q_mean_inv_angstrom", "analyzer_energy_kev", "two_theta_deg"};
    for (int i = 0; i < 3; ++i) {
        const auto path = QString(group) + names[i];
        if (!hdf.contains(path))
            continue;
        if (HdfFile::inspect(file, path).dims != QVector<quint64>{rows})
            fail(QCoreApplication::translate("SpectrumProcessing", "晶体几何数据集长度不匹配：%1").arg(path));
        geometry[i] = hdf.read(path);
    }
    for (qsizetype i = 0; i < labels.size(); ++i)
        out.geometry.insert(
            labels[i], {geometry[0].value(i, NAN), geometry[1].value(i, NAN), geometry[2].value(i, NAN)});
    const auto widePath = QString(group) + "compton_wide_intensity";
    const auto sourcePath = QString(group) + "source_index";
    if (layout.intensity == SpectrumHdfLayout{}.intensity && !layout.energyFirst && hdf.contains(widePath)) {
        if (HdfFile::inspect(file, widePath).dims != shape ||
            HdfFile::inspect(file, sourcePath).dims != shape)
            fail(QCoreApplication::translate("SpectrumProcessing", "拼接文件中的宽谱或来源索引维度不匹配。"));
        const auto original = hdf.read(widePath), origin = hdf.read(sourcePath);
        out.sources = hdf.strings(QString(group) + "source_files", 256);
        for (qsizetype r = 0; r < labels.size(); ++r) {
            QVector<double> raw(n);
            QVector<int> codes(n);
            for (qsizetype i = 0; i < qsizetype(n); ++i) {
                const auto index = r * qsizetype(n) + (descending ? qsizetype(n) - 1 - i : i);
                raw[i] = original[index];
                const double code = origin[index];
                if (!std::isfinite(code) || code < 0 || code >= out.sources.size() ||
                    std::floor(code) != code)
                    fail(QCoreApplication::translate("SpectrumProcessing", "拼接文件的来源索引无效。"));
                codes[i] = int(code);
            }
            out.wideCrystals.insert(labels[r], std::move(raw));
            if (std::any_of(codes.cbegin(), codes.cend(), [](int code) { return code != 0; }))
                out.pointSources.insert(labels[r], std::move(codes));
        }
    }
    return out;
}

double inferAngle(const QVector<double> &energy, double e0, double q) {
    if (!(e0 > 0 && q > 0) || !std::isfinite(q) || e0 + energy.first() / 1000. <= 0)
        return NAN;
    constexpr double k = 6.283185307179586 / 12.398419297617678;
    auto mean = [&](double sineSquared) {
        double sum = 0;
        for (double loss : energy) {
            double delta = loss / 1000.;
            sum += std::sqrt(delta * delta + 4 * e0 * (e0 + delta) * sineSquared);
        }
        return k * sum / energy.size();
    };
    if (q <= mean(0) || q >= mean(1))
        return NAN;
    double lo = 0, hi = 1;
    for (int i = 0; i < 42; ++i) {
        double mid = (lo + hi) / 2;
        if (mean(mid) < q)
            lo = mid;
        else
            hi = mid;
    }
    return (360. / 3.14159265358979323846) * std::asin(std::sqrt((lo + hi) / 2));
}

ProcessedSpectrum readReducedTable(const QString &file, const SpectrumHdfLayout &layout) {
    QFile input(file);
    if (!input.open(QIODevice::ReadOnly))
        fail(input.errorString());
    QTextStream stream(&input);
    const auto headers = stream.readLine().split('\t');
    auto axis = headers.indexOf("Energy Transfer (eV)");
    if (axis < 0)
        fail(QCoreApplication::translate("SpectrumProcessing", "请选择含 Energy Transfer (eV) 的逐晶体 *_all_data.txt。"));
    ProcessedSpectrum out;
    out.sources = {file};
    out.intensityUnit = layout.externalIntensityUnit;
    QVector<QPair<qsizetype, QString>> columns;
    for (qsizetype i = 0; i < headers.size(); ++i) {
        auto label = headers[i].trimmed();
        if (i == axis || label.isEmpty() || label.startsWith("Unnamed:"))
            continue;
        if (out.crystals.contains(label))
            fail(QCoreApplication::translate("SpectrumProcessing", "强度表含重复晶体标签：%1").arg(label));
        out.crystals.insert(label, {});
        columns.append({i, label});
    }
    if (columns.isEmpty() || columns.size() > maxCrystals)
        fail(QCoreApplication::translate("SpectrumProcessing", "逐晶体强度列数量无效。"));
    while (!stream.atEnd()) {
        auto line = stream.readLine();
        if (line.trimmed().isEmpty())
            continue;
        const auto fields = line.split('\t');
        bool ok;
        const double e = fields.value(axis).toDouble(&ok);
        if (!ok || !std::isfinite(e) || (!out.energyEv.isEmpty() && e <= out.energyEv.last()) ||
            fields.size() != headers.size() || out.energyEv.size() >= maxPoints ||
            (out.energyEv.size() + 1) * columns.size() > maxValues)
            fail(QCoreApplication::translate("SpectrumProcessing", "逐晶体表的列数、能量顺序或大小无效。"));
        out.energyEv.append(e);
        for (const auto &[i, label] : columns) {
            auto field = fields[i].trimmed();
            double value = field.toDouble(&ok);
            if (!ok && !field.isEmpty() && field.compare("nan", Qt::CaseInsensitive) != 0)
                fail(QCoreApplication::translate("SpectrumProcessing", "晶体 %1 含无法读取的强度。").arg(label));
            out.crystals[label].append(ok && std::isfinite(value) ? value : NAN);
        }
    }
    if (out.energyEv.size() < 2)
        fail(QCoreApplication::translate("SpectrumProcessing", "能量点不足。"));
    if (file.endsWith("_all_data.txt", Qt::CaseInsensitive)) {
        const auto path = file.chopped(13) + "_rois.txt";
        if (QFileInfo::exists(path))
            out.geometry = SpectrumProcessing::loadGeometry(path, out);
    }
    return out;
}

using Parameters = Eigen::Matrix<double, 5, 1>;
double pearsonValue(double x, const Parameters &p, Parameters *jacobian = nullptr) {
    const double width = std::exp(p[3]), m = std::exp(p[4]);
    const double factor = std::expm1(std::log(2.) / m);
    const double u = 2 * (x - p[2]) / width, q = factor * u * u;
    const double logarithm = std::log1p(q), profile = std::exp(-m * logarithm);
    if (jacobian) {
        const double height = p[1] * profile;
        *jacobian << 1., profile, height * 4 * m * factor * u / (width * (1 + q)),
            height * 2 * m * q / (1 + q),
            height * (std::log(2.) * (1 + factor) * u * u / (1 + q) - m * logarithm);
    }
    return p[0] + p[1] * profile;
}
} // namespace

ProcessedSpectrum SpectrumProcessing::fromResult(const SpectrumResult &result, const CrystalRois &rois,
                                                 bool normalized) {
    ProcessedSpectrum out;
    out.energyEv = result.energyTransferEv;
    out.intensityUnit = normalized ? "Intensity / |I0|" : "Counts";
    for (qsizetype i = 0; i < result.included.size(); ++i)
        if (result.included[i]) {
            const auto &label = rois.items()[i].label;
            out.crystals.insert(label, result.crystals[i]);
            out.geometry.insert(label, {result.momentum.value(i).mean, result.analyzerEnergyKev.value(i, NAN),
                                        result.scatteringAngleDeg.value(i, NAN)});
        }
    return out;
}

QMap<QString, SpectrumGeometry> SpectrumProcessing::loadGeometry(const QString &path,
                                                                 const ProcessedSpectrum &spectrum) {
    QFile rois(path);
    if (!rois.open(QIODevice::ReadOnly))
        fail(rois.errorString());
    QTextStream rows(&rois);
    const auto names = rows.readLine().split('\t');
    const auto labelCol = names.indexOf("crystal"), angleCol = names.indexOf("two_theta_deg");
    auto qCol = names.indexOf("q_ave"), centerCol = names.indexOf("center");
    if (qCol < 0)
        qCol = names.indexOf("q_ave_A^-1");
    if (centerCol < 0)
        centerCol = names.indexOf("center_eV");
    if (labelCol < 0 || qCol < 0 || centerCol < 0)
        fail(QCoreApplication::translate("SpectrumProcessing", "晶体几何表须含 crystal、q_ave（Å⁻¹）和 center（eV）列。"));
    const auto enabledCol = names.indexOf("add");
    QSet<QString> seen;
    QMap<QString, SpectrumGeometry> geometry;
    while (!rows.atEnd()) {
        auto fields = rows.readLine().split('\t');
        const auto label = fields.value(labelCol).trimmed();
        if (!spectrum.crystals.contains(label))
            continue;
        if (seen.contains(label)) {
            geometry.remove(label);
            continue;
        }
        seen.insert(label);
        if (enabledCol >= 0 && fields.value(enabledCol).trimmed() == "0")
            continue;
        bool qOk, eOk, aOk;
        const double q = fields.value(qCol).toDouble(&qOk);
        const double e0 = fields.value(centerCol).toDouble(&eOk) / 1000.;
        const double angle = fields.value(angleCol).toDouble(&aOk);
        if (qOk && eOk && std::isfinite(q) && std::isfinite(e0))
            geometry[label] = {q, e0, aOk ? angle : inferAngle(spectrum.energyEv, e0, q)};
    }
    return geometry;
}

ProcessedSpectrum SpectrumProcessing::load(const QString &path, const SpectrumHdfLayout &layout) {
    if (path.endsWith(".txt", Qt::CaseInsensitive) || path.endsWith(".tsv", Qt::CaseInsensitive))
        return readReducedTable(path, layout);
    auto result = readSpectrum(path, layout);
    const auto companion = QFileInfo(path).dir().filePath("spectrum_rois.txt");
    if (QFileInfo(path).fileName() == "spectra.h5" && QFileInfo::exists(companion) &&
        std::none_of(result.geometry.cbegin(), result.geometry.cend(),
                     [](const SpectrumGeometry &g) { return std::isfinite(g.qMean); }))
        result.geometry = loadGeometry(companion, result);
    return result;
}

ProcessedSpectrum SpectrumProcessing::replaceEdges(const ProcessedSpectrum &wide,
                                                   const QVector<ProcessedSpectrum> &edges,
                                                   const std::atomic_bool &stop) {
    if (wide.crystals.isEmpty())
        fail(QCoreApplication::translate("SpectrumProcessing", "请先载入康普顿宽谱。"));
    auto result = wide;
    result.pointSources.clear();
    result.wideCrystals = edges.isEmpty() ? QMap<QString, QVector<double>>{} : wide.crystals;
    result.sources = {wide.sources.isEmpty() ? QCoreApplication::translate("SpectrumProcessing", "当前原始谱") : wide.sources.join(" | ")};
    for (const auto &edge : edges) {
        checkCancelled(stop);
        if (edge.intensityUnit != wide.intensityUnit)
            fail(QCoreApplication::translate("SpectrumProcessing", "宽谱与近边谱的强度单位不一致。"));
        QVector<qsizetype> indices;
        for (double value : edge.energyEv) {
            auto it = std::lower_bound(wide.energyEv.cbegin(), wide.energyEv.cend(), value);
            if (it != wide.energyEv.cbegin() &&
                (it == wide.energyEv.cend() || std::abs(*(it - 1) - value) < std::abs(*it - value)))
                --it;
            if (it == wide.energyEv.cend() || std::abs(*it - value) > 1e-6)
                fail(QCoreApplication::translate("SpectrumProcessing", "近边能量 %1 eV 与宽谱网格不对齐。请按共同能量网格重新导出。")
                         .arg(value, 0, 'g', 12));
            auto index = it - wide.energyEv.cbegin();
            if (!indices.isEmpty() && index <= indices.last())
                fail(QCoreApplication::translate("SpectrumProcessing", "多个近边点对应同一宽谱点。"));
            indices.append(index);
        }
        qsizetype replaced = 0;
        const int source = int(result.sources.size());
        for (auto row = edge.crystals.cbegin(); row != edge.crystals.cend(); ++row) {
            auto target = result.crystals.find(row.key());
            if (target == result.crystals.end())
                continue;
            auto &origin = result.pointSources[row.key()];
            if (origin.isEmpty())
                origin.fill(0, result.energyEv.size());
            for (qsizetype i = 0; i < indices.size(); ++i) {
                if (!std::isfinite(row.value()[i]))
                    continue;
                const auto index = indices[i];
                if (origin[index] != 0)
                    fail(QCoreApplication::translate("SpectrumProcessing", "近边文件在晶体 %1 的能区重叠，请重新选择。").arg(row.key()));
                target.value()[index] = row.value()[i];
                origin[index] = source;
                ++replaced;
            }
        }
        if (!replaced)
            fail(QCoreApplication::translate("SpectrumProcessing", "近边文件没有可替换的同名晶体强度。"));
        result.sources.append(edge.sources.join(" | "));
    }
    return result;
}

ProcessedSpectrum SpectrumProcessing::stitch(const QStringList &files, const SpectrumHdfLayout &layout,
                                             const std::atomic_bool &stop) {
    if (files.isEmpty())
        fail(QCoreApplication::translate("SpectrumProcessing", "请选择谱学 HDF5 文件。"));
    QVector<ProcessedSpectrum> inputs;
    ProcessedSpectrum out;
    QSet<QString> seenFiles;
    QSet<QString> labels;
    qsizetype totalValues = 0;
    for (const auto &file : files) {
        checkCancelled(stop);
        const auto canonical = QFileInfo(file).canonicalFilePath();
        if (seenFiles.contains(canonical))
            continue;
        seenFiles.insert(canonical);
        ProcessedSpectrum input;
        try {
            input = readSpectrum(file, layout);
        } catch (const std::exception &e) {
            fail(QStringLiteral("%1\n%2").arg(QFileInfo(file).fileName(), QString::fromUtf8(e.what())));
        }
        if (!out.intensityUnit.isEmpty() && out.intensityUnit != input.intensityUnit)
            fail(QCoreApplication::translate("SpectrumProcessing", "所选文件的强度单位不一致，请先统一归一化方式。"));
        totalValues += input.crystals.size() * input.energyEv.size();
        if (totalValues > maxValues)
            fail(QCoreApplication::translate("SpectrumProcessing", "待拼接谱学数据超出内存限制，请分批处理。"));
        out.intensityUnit = input.intensityUnit;
        out.energyEv += input.energyEv;
        out.sources += input.sources;
        for (const auto &label : input.crystals.keys())
            labels.insert(label);
        inputs.append(std::move(input));
    }
    std::sort(out.energyEv.begin(), out.energyEv.end());
    qsizetype size = 0;
    for (double energy : std::as_const(out.energyEv))
        if (size == 0 || energy - out.energyEv[size - 1] > energyTolerance)
            out.energyEv[size++] = energy;
    out.energyEv.resize(size);
    if (size < 2 || size > maxPoints || labels.size() > maxCrystals || size * labels.size() > maxValues)
        fail(QCoreApplication::translate("SpectrumProcessing", "拼接后的能量点或晶体数超出支持范围。"));
    for (const auto &label : labels) {
        checkCancelled(stop);
        QVector<double> row(size, NAN);
        QVector<int> counts(size, 0);
        bool present = false;
        for (const auto &input : inputs) {
            auto it = input.crystals.constFind(label);
            if (it == input.crystals.cend())
                continue;
            qsizetype target = 0;
            for (qsizetype i = 0; i < input.energyEv.size(); ++i) {
                double energy = input.energyEv[i], value = it.value()[i];
                while (target + 1 < size && out.energyEv[target + 1] <= energy)
                    ++target;
                if (!std::isfinite(value))
                    continue;
                int count = ++counts[target];
                row[target] = count == 1 ? value : row[target] * (1. - 1. / count) + value / count;
                present = true;
            }
        }
        if (present)
            out.crystals.insert(label, std::move(row));
    }
    if (out.crystals.isEmpty())
        fail(QCoreApplication::translate("SpectrumProcessing", "所选数据集没有有效晶体强度。"));
    return out;
}

CrystalBackground SpectrumProcessing::fit(const QVector<double> &energy, const QVector<double> &intensity,
                                          const BackgroundSettings &settings, const std::atomic_bool &stop) {
    if (energy.size() != intensity.size() || energy.size() < 2 || settings.windows.isEmpty())
        fail(QCoreApplication::translate("SpectrumProcessing", "缺少可拟合的晶体谱或能量区间。"));
    for (auto window : settings.windows)
        if (!std::isfinite(window.first) || !std::isfinite(window.second) || window.first >= window.second)
            fail(QCoreApplication::translate("SpectrumProcessing", "拟合区间的结束值须大于起始值。"));
    const double span = energy.last() - energy.first(), origin = energy.first() + .5 * span;
    if (!(span > 0) || !std::isfinite(span))
        fail(QCoreApplication::translate("SpectrumProcessing", "能量范围无效。"));
    QVector<double> x, y;
    double scale = 0;
    for (qsizetype i = 0; i < energy.size(); ++i) {
        if (!std::isfinite(intensity[i]))
            continue;
        if (std::any_of(settings.windows.begin(), settings.windows.end(), [&](auto window) {
                return energy[i] >= window.first && energy[i] <= window.second;
            })) {
            x.append((energy[i] - origin) / span);
            y.append(intensity[i]);
            scale = std::max(scale, std::abs(intensity[i]));
        }
    }
    if (x.size() < (settings.pearson ? 8 : 2))
        fail(QCoreApplication::translate("SpectrumProcessing", "拟合区间内至少需要 %1 个有效点。").arg(settings.pearson ? 8 : 2));
    scale = scale > 0 ? scale : 1;
    for (double &value : y)
        value /= scale;
    const double meanX = std::accumulate(x.begin(), x.end(), 0.) / x.size();
    const double meanY = std::accumulate(y.begin(), y.end(), 0.) / y.size();
    double xx = 0, xy = 0, yy = 0;
    for (qsizetype i = 0; i < x.size(); ++i) {
        const double dx = x[i] - meanX, dy = y[i] - meanY;
        xx += dx * dx;
        xy += dx * dy;
        yy += dy * dy;
    }
    if (!(xx > 0))
        fail(QCoreApplication::translate("SpectrumProcessing", "拟合区间缺少不同的能量点。"));
    const double slope = xy / xx, intercept = meanY - slope * meanX;
    Parameters best;
    double bestError = INFINITY;
    auto squaredError = [&](const Parameters &p) {
        double error = 0;
        for (qsizetype i = 0; i < x.size(); ++i) {
            const double residual = y[i] - pearsonValue(x[i], p);
            error += residual * residual;
        }
        return error;
    };
    if (settings.pearson) {
        auto extrema = std::minmax_element(y.begin(), y.end());
        const double height = *extrema.second - *extrema.first;
        if (height < 1e-10)
            fail(QCoreApplication::translate("SpectrumProcessing", "所选区间没有可辨识的康普顿峰。"));
        const Parameters lower{-INFINITY, 0., -1.5, std::log(1e-5), std::log(.51)};
        const Parameters upper{INFINITY, 100., 1.5, std::log(5.), std::log(1000.)};
        for (double initialWidth : {.08, .25, .7}) {
            Parameters p;
            p << *extrema.first, height, x[extrema.second - y.begin()], std::log(initialWidth), std::log(2.);
            double error = squaredError(p), damping = 1e-3;
            for (int iteration = 0; iteration < 400; ++iteration) {
                checkCancelled(stop);
                Eigen::Matrix<double, 5, 5> normal = Eigen::Matrix<double, 5, 5>::Zero();
                Parameters gradient = Parameters::Zero();
                for (qsizetype i = 0; i < x.size(); ++i) {
                    Parameters j;
                    const double residual = y[i] - pearsonValue(x[i], p, &j);
                    normal.noalias() += j * j.transpose();
                    gradient.noalias() += j * residual;
                }
                for (int k = 1; k < 5; ++k) {
                    if ((p[k] <= lower[k] && gradient[k] < 0) || (p[k] >= upper[k] && gradient[k] > 0)) {
                        normal.row(k).setZero();
                        normal.col(k).setZero();
                        normal(k, k) = 1.;
                        gradient[k] = 0.;
                    }
                }
                if (gradient.norm() < 1e-10 || error < 1e-16) {
                    if (error < bestError) {
                        best = p;
                        bestError = error;
                    }
                    break;
                }
                normal.diagonal().array() += damping * (normal.diagonal().array() + 1e-12);
                Parameters candidate = p + normal.ldlt().solve(gradient);
                if (!candidate.allFinite())
                    break;
                candidate = candidate.cwiseMax(lower).cwiseMin(upper);
                const double trial = squaredError(candidate);
                if (trial < error) {
                    const double improvement = error - trial;
                    p = candidate;
                    error = trial;
                    damping = std::max(1e-12, damping * .3);
                    if (improvement < 1e-12 * (1 + error)) {
                        if (error < bestError) {
                            best = p;
                            bestError = error;
                        }
                        break;
                    }
                } else {
                    damping *= 10;
                    if (damping > 1e14)
                        break;
                }
            }
        }
        if (!std::isfinite(bestError) || best[1] <= 1e-9)
            fail(QCoreApplication::translate("SpectrumProcessing", "Pearson VII 拟合未收敛或峰高无效。请调整康普顿拟合区间。"));
    } else {
        bestError = 0;
        for (qsizetype i = 0; i < x.size(); ++i) {
            double residual = y[i] - (intercept + slope * x[i]);
            bestError += residual * residual;
        }
    }
    CrystalBackground result;
    result.settings = settings;
    result.rSquared = yy > 0 ? 1 - bestError / yy : NAN;
    if (settings.pearson) {
        result.centerEv = origin + span * best[2];
        result.fwhmEv = span * std::exp(best[3]);
        result.exponent = std::exp(best[4]);
        result.peakHeight = scale * best[1];
    }
    for (qsizetype i = 0; i < energy.size(); ++i) {
        const double pos = (energy[i] - origin) / span;
        const double value = scale * (settings.pearson ? pearsonValue(pos, best) : intercept + slope * pos);
        result.background.append(std::isfinite(intensity[i]) ? value : NAN);
        result.corrected.append(intensity[i] - value);
    }
    return result;
}

void SpectrumProcessing::save(const QString &path, const ProcessedSpectrum &data,
                              const QMap<QString, CrystalBackground> &backgrounds) {
    if (data.crystals.isEmpty() || data.energyEv.size() < 2)
        fail(QCoreApplication::translate("SpectrumProcessing", "没有可保存的谱学结果。"));
    if (QFileInfo::exists(path))
        fail(QCoreApplication::translate("SpectrumProcessing", "输出文件已存在，请使用新文件名。"));
    const QString staging = path + ".partial-" + QUuid::createUuid().toString(QUuid::Id128);
    try {
        {
            ImagingH5 hdf(staging, true);
            const auto labels = data.crystals.keys();
            const auto n = data.energyEv.size();
            QVector<quint64> shape{quint64(labels.size()), quint64(n)};
            hdf.write(QString(group) + "energy_transfer_ev", {quint64(n)}, data.energyEv.constData());
            hdf.writeStrings(QString(group) + "crystal_labels", labels);
            hdf.writeText(QString(group) + "intensity_unit", data.intensityUnit);
            if (!data.sources.isEmpty())
                hdf.writeStrings(QString(group) + "source_files", data.sources);
            auto matrix = [&](const QString &name,
                              const std::function<QVector<double>(const QString &)> &row) {
                QVector<double> values;
                values.reserve(labels.size() * n);
                for (const auto &label : labels) {
                    auto next = row(label);
                    if (next.size() != n)
                        fail(QCoreApplication::translate("SpectrumProcessing", "晶体 %1 的强度长度与能量轴不一致。").arg(label));
                    values += next;
                }
                hdf.write(QString(group) + name, shape, values.constData());
            };
            matrix("intensity", [&](const QString &label) { return data.crystals.value(label); });
            if (!data.wideCrystals.isEmpty())
                matrix("compton_wide_intensity",
                       [&](const QString &label) { return data.wideCrystals.value(label); });
            if (!data.processing.isEmpty())
                hdf.writeText(
                    QString(group) + "processing_json",
                    QString::fromUtf8(QJsonDocument(data.processing).toJson(QJsonDocument::Compact)));

            if (!data.geometry.isEmpty()) {
                QVector<double> q, e0, angles;
                for (const auto &label : labels) {
                    const auto g = data.geometry.value(label);
                    q.append(g.qMean);
                    e0.append(g.analyzerKev);
                    angles.append(g.angleDeg);
                }
                hdf.write(QString(group) + "q_mean_inv_angstrom", {quint64(labels.size())}, q.constData());
                hdf.write(QString(group) + "analyzer_energy_kev", {quint64(labels.size())}, e0.constData());
                hdf.write(QString(group) + "two_theta_deg", {quint64(labels.size())}, angles.constData());
            }
            if (!data.pointSources.isEmpty()) {
                QVector<quint64> sources;
                sources.reserve(labels.size() * n);
                for (const auto &label : labels) {
                    auto it = data.pointSources.constFind(label);
                    for (qsizetype i = 0; i < n; ++i)
                        sources.append(it == data.pointSources.cend() ? 0 : quint64(it.value().value(i)));
                }
                hdf.writeUInt(QString(group) + "source_index", shape, sources.constData());
            }

            if (!backgrounds.isEmpty()) {
                matrix("background", [&](const QString &label) {
                    auto fit = backgrounds.constFind(label);
                    return fit == backgrounds.cend() ? QVector<double>(n, NAN) : fit->background;
                });
                matrix("corrected", [&](const QString &label) {
                    auto fit = backgrounds.constFind(label);
                    return fit == backgrounds.cend() ? QVector<double>(n, NAN) : fit->corrected;
                });
                QJsonObject metadata;
                for (auto fit = backgrounds.cbegin(); fit != backgrounds.cend(); ++fit) {
                    QJsonArray windows;
                    for (auto range : fit->settings.windows)
                        windows.append(QJsonArray{range.first, range.second});
                    metadata[fit.key()] = QJsonObject{
                        {"model", fit->method.isEmpty() ? (fit->settings.pearson ? "Pearson VII" : "linear")
                                                        : fit->method},
                        {"fit_windows_ev", windows},
                        {"r_squared", fit->rSquared},
                        {"center_ev", fit->centerEv},
                        {"fwhm_ev", fit->fwhmEv},
                        {"exponent", fit->exponent}};
                }
                hdf.writeText(QString(group) + "background_fit_json",
                              QString::fromUtf8(QJsonDocument(metadata).toJson(QJsonDocument::Compact)));
            }
            hdf.attributes(QString(group).chopped(1), {{"schema", "xri-spectrum-1"}});
        }
        if (!QFile::rename(staging, path))
            fail(QCoreApplication::translate("SpectrumProcessing", "无法完成 HDF5 文件写入：%1").arg(path));
    } catch (...) {
        QFile::remove(staging);
        throw;
    }
}
