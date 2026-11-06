#include "Spectroscopy.h"
#include <QDateTime>
#include <QDir>
#include <QFileInfo>
#include <QJsonArray>
#include <QJsonDocument>
#include <QJsonObject>
#include <QSaveFile>
#include <QTextStream>
#include <QUuid>
#include <algorithm>
#include <array>
#include <numeric>

namespace {
[[noreturn]] void fail(const QString &message) {
    throw std::runtime_error(message.toUtf8().constData());
}
void checkCancel(const std::atomic_bool &cancel) {
    if (cancel.load())
        fail(QStringLiteral("已取消谱学计算。"));
}
QString detectorPath(int d) {
    return d == 0 ? "/entry/data/D_LAMBDA" : "/entry/data/D_MINIPIX";
}
QVector<double> vectorData(const QString &file, const QString &path, int size = -1) {
    auto data = HdfFile::read({file, path, {}, 0, -1});
    if (!data.error.isEmpty() || data.sampled || data.dims.size() != 1 ||
        (size >= 0 && data.values.size() != size))
        fail(QStringLiteral("向量读取失败或长度不匹配：%1\n%2").arg(path, data.error));
    return data.values;
}
QString energyPath(const QString &file, const QString &requested) {
    if (requested != "auto")
        return requested;
    for (const auto &path : {"/entry/instrument/M_DCM_B5Energy", "/entry/instrument/M_DCM_B5Energy_readback"})
        try {
            if (HdfFile::inspect(file, path).dims.size() == 1)
                return path;
        } catch (const std::exception &) {
        }
    fail(QStringLiteral("找不到实测能量轴，请指定能量数据集。"));
}
DataSlice imageData(const QString &name, const QSize &size) {
    DataSlice image;
    image.path = name;
    image.width = size.width();
    image.height = size.height();
    image.xAxis = 1;
    image.yAxis = 0;
    image.dims = {quint64(size.height()), quint64(size.width())};
    image.starts = {0, 0};
    image.strides = {1, 1};
    image.values.fill(0, size.width() * size.height());
    return image;
}
void imageStats(DataSlice &image) {
    image.minimum = INFINITY;
    image.maximum = -INFINITY;
    long double sum = 0;
    image.finiteCount = 0;
    for (double v : image.values)
        if (std::isfinite(v)) {
            image.minimum = std::min(image.minimum, v);
            image.maximum = std::max(image.maximum, v);
            sum += v;
            ++image.finiteCount;
        }
    image.sum = double(sum);
    image.mean = image.finiteCount ? double(sum / image.finiteCount) : NAN;
    if (!image.finiteCount)
        image.minimum = image.maximum = NAN;
}
void validateSettings(const SpectrumSettings &s) {
    if (!std::isfinite(s.filterFraction) || s.filterFraction < 0 || s.filterFraction > 1 || s.roiSize < 1 ||
        s.roiSize > 2048 || s.fitModel < 0 || s.fitModel > 1 || !(s.centerMax > s.centerMin) ||
        !(s.maxFwhmEv > 0) || s.minRSquared < 0 || s.minRSquared > 1)
        fail(QStringLiteral("掩膜或拟合参数无效。"));
}
SpectrumScan loadScan(const QString &file, const SpectrumSettings &settings, bool elastic,
                      const ElasticCalibration *calibration, bool projection, const std::atomic_bool &cancel,
                      Spectroscopy::Progress progress) {
    checkCancel(cancel);
    if (HdfFile::formatFromFileName(file) == ScanFormat::GridScan)
        fail(QStringLiteral("grid_scan 属于成像数据，不能作为一维谱学扫描：%1").arg(file));
    SpectrumScan scan;
    scan.file = QFileInfo(file).absoluteFilePath();
    scan.energyPath = energyPath(file, elastic ? settings.elasticEnergy : settings.xrsEnergy);
    scan.energyKev = vectorData(file, scan.energyPath);
    int count = int(scan.energyKev.size());
    if (count < 2)
        fail(QStringLiteral("能量扫描至少需要两个点。"));
    int direction = 0;
    for (int i = 0; i < count; ++i) {
        if (!std::isfinite(scan.energyKev[i]) || scan.energyKev[i] <= 0)
            fail(QStringLiteral("能量轴包含非有限值或非正能量。"));
        if (i) {
            double d = scan.energyKev[i] - scan.energyKev[i - 1];
            int sign = d > 0 ? 1 : d < 0 ? -1 : 0;
            if (!sign || (direction && direction != sign))
                fail(QStringLiteral("能量必须严格单调；重复或多维网格数据不能直接作为一维谱处理。"));
            direction = sign;
        }
    }
    const bool normalize = elastic ? settings.normalizeElastic : settings.normalizeXrs;
    scan.i0Path = normalize ? (elastic ? settings.elasticI0 : settings.xrsI0) : QString();
    scan.i0Before = normalize ? vectorData(file, scan.i0Path, count) : QVector<double>(count, 1);
    for (double &v : scan.i0Before)
        v = std::abs(v);
    scan.i0After = scan.i0Before;
    QVector<bool> glitch(count, false);
    if (normalize && !elastic && settings.repairGlitches) {
        for (int i = 1; i < count - 1; ++i) {
            double a = scan.i0After[i - 1], b = scan.i0After[i], c = scan.i0After[i + 1];
            if (std::isfinite(a) && std::isfinite(b) && std::isfinite(c) && a > 0 && c > 0 && b / a < 0.7 &&
                b / c < 0.7) {
                scan.i0After[i] = (a + c) / 2;
                glitch[i] = true;
                scan.glitches.append(i);
            }
        }
    }
    for (double v : scan.i0After)
        if (!std::isfinite(v) || v <= 0)
            fail(QStringLiteral("I0 存在零值或非有限值且无法修复，请检查扫描或关闭归一化。"));
    if (calibration) {
        scan.raw.resize(75);
        scan.filtered.resize(75);
        for (int i = 0; i < 75; ++i)
            if (calibration->adjusted.items()[i].state == RoiState::Valid) {
                scan.raw[i].fill(0, count);
                scan.filtered[i].fill(0, count);
            }
    }
    for (int d = 0; d < 2; ++d) {
        checkCancel(cancel);
        const auto node = HdfFile::inspect(file, detectorPath(d));
        if (node.dims.size() != 3 || node.dims[0] != quint64(count))
            fail(QStringLiteral("探测器帧数与能量轴不一致：%1").arg(detectorPath(d)));
        QSize size(int(node.dims[2]), int(node.dims[1]));
        if (calibration && calibration->sizes[d] != size)
            fail(QStringLiteral("XRS 探测器尺寸与弹性标定不一致。"));
        DataSlice image = imageData(detectorPath(d), size);
        QVector<double> pending, previous;
        auto process = [&](int index, QVector<double> &raw) {
            const double divisor = normalize ? scan.i0After[index] : 1;
            const double *pixels = raw.constData();
            if (projection) {
                double *sum = image.values.data();
                for (qsizetype k = 0; k < raw.size(); ++k) {
                    if (!std::isfinite(pixels[k]))
                        fail(QStringLiteral("探测器包含非有限计数。"));
                    sum[k] += pixels[k] / divisor;
                }
            }
            if (calibration) {
                int start = d == 0 ? 0 : 60, end = d == 0 ? 60 : 75;
                for (int r = start; r < end; ++r) {
                    if (calibration->adjusted.items()[r].state != RoiState::Valid)
                        continue;
                    long double whole = 0, filtered = 0;
                    for (int pixel : calibration->pixels[r]) {
                        double v = pixels[pixel];
                        if (!std::isfinite(v))
                            fail(QStringLiteral("ROI 内出现非有限计数。"));
                        whole += v;
                    }
                    for (int pixel : calibration->filteredPixels[r])
                        filtered += pixels[pixel];
                    scan.raw[r][index] = double(whole / divisor);
                    scan.filtered[r][index] = double(filtered / divisor);
                }
            }
        };
        HdfFile::visitFrames(file, detectorPath(d), cancel,
                             [&](int i, const double *raw, int width, int height) {
                                 if (i == 0)
                                     pending = QVector<double>(raw, raw + width * height);
                                 else {
                                     if (glitch[i - 1]) {
                                         double *values = pending.data();
                                         const double *prev = previous.constData();
                                         for (int k = 0; k < width * height; ++k)
                                             values[k] = (prev[k] + raw[k]) / 2;
                                     }
                                     process(i - 1, pending);
                                     previous = std::move(pending);
                                     pending = QVector<double>(raw, raw + width * height);
                                 }
                                 if (progress && (i % 16 == 0 || i == count - 1))
                                     progress(QStringLiteral("%1 · %2").arg(QFileInfo(file).fileName(),
                                                                            d == 0 ? "LAMBDA" : "MINIPIX"),
                                              (d * 100 + i * 100 / count) / 2);
                             });
        process(count - 1, pending);
        if (projection) {
            imageStats(image);
            scan.images.append(image);
        }
    }
    return scan;
}
bool solve4(double matrix[4][4], double rhs[4], double solution[4]) {
    for (int i = 0; i < 4; ++i) {
        int best = i;
        for (int j = i + 1; j < 4; ++j)
            if (std::abs(matrix[j][i]) > std::abs(matrix[best][i]))
                best = j;
        if (std::abs(matrix[best][i]) < 1e-18)
            return false;
        for (int k = 0; k < 4; ++k)
            std::swap(matrix[i][k], matrix[best][k]);
        std::swap(rhs[i], rhs[best]);
        double scale = matrix[i][i];
        for (int k = i; k < 4; ++k)
            matrix[i][k] /= scale;
        rhs[i] /= scale;
        for (int j = 0; j < 4; ++j)
            if (j != i) {
                double factor = matrix[j][i];
                for (int k = i; k < 4; ++k)
                    matrix[j][k] -= factor * matrix[i][k];
                rhs[j] -= factor * rhs[i];
            }
    }
    for (int i = 0; i < 4; ++i)
        solution[i] = rhs[i];
    return true;
}
double peak(double x, const std::array<double, 4> &p, int model, double *derivative = nullptr) {
    const double width = std::exp(p[2]), z = (x - p[1]) / width;
    const double g = model == 0 ? 1 / (1 + z * z) : std::exp(-0.5 * z * z);
    if (derivative) {
        derivative[0] = g;
        derivative[1] = model == 0 ? p[0] * 2 * z / width * g * g : p[0] * g * z / width;
        derivative[2] = model == 0 ? p[0] * 2 * z * z * g * g : p[0] * g * z * z;
        derivative[3] = 1;
    }
    return p[0] * g + p[3];
}
double interpolate(const QVector<double> &x, const QVector<double> &y, double target) {
    auto it = std::lower_bound(x.begin(), x.end(), target);
    if (it == x.begin())
        return y.first();
    if (it == x.end())
        return y.last();
    int i = int(it - x.begin());
    double t = (target - x[i - 1]) / (x[i] - x[i - 1]);
    return y[i - 1] * (1 - t) + y[i] * t;
}
QJsonArray numbers(const QVector<double> &values) {
    QJsonArray out;
    for (auto v : values)
        out.append(std::isfinite(v) ? QJsonValue(v) : QJsonValue());
    return out;
}
} // namespace

ElasticFit Spectroscopy::fitPeak(const QVector<double> &energy, const QVector<double> &intensity,
                                 const SpectrumSettings &s) {
    ElasticFit result;
    if (energy.size() != intensity.size() || energy.size() < 8) {
        result.reason = QStringLiteral("拟合至少需要 8 个有效能量点。");
        return result;
    }
    for (int i = 0; i < energy.size(); ++i)
        if (!std::isfinite(energy[i]) || !std::isfinite(intensity[i])) {
            result.reason = QStringLiteral("拟合数据含非有限值。");
            return result;
        }
    auto minmax = std::minmax_element(intensity.begin(), intensity.end());
    double scale = *minmax.second;
    if (!(scale > 0) || *minmax.second - *minmax.first <= std::max(1e-12, std::abs(scale) * 1e-10)) {
        result.reason = QStringLiteral("谱线为空、全零或无峰形变化。");
        return result;
    }
    double origin = (energy.first() + energy.last()) * 500,
           span = std::abs(energy.last() - energy.first()) * 1000;
    if (!(span > 0)) {
        result.reason = QStringLiteral("能量范围为零。");
        return result;
    }
    QVector<double> x, y;
    for (int i = 0; i < energy.size(); ++i) {
        x.append(energy[i] * 1000 - origin);
        y.append(intensity[i] / scale);
    }
    int maximum = int(std::max_element(y.begin(), y.end()) - y.begin());
    double baseline = (y.first() + y.last()) / 2;
    std::array<double, 4> p{std::max(0.01, y[maximum] - baseline), x[maximum],
                            std::log(std::min(1.0, span / 5)), baseline};
    auto error = [&](const auto &parameters) {
        double sum = 0;
        for (int i = 0; i < x.size(); ++i) {
            double delta = y[i] - peak(x[i], parameters, s.fitModel);
            sum += delta * delta;
        }
        return sum;
    };
    double current = error(p), lambda = 1e-3;
    bool converged = false;
    for (int iteration = 0; iteration < 300; ++iteration) {
        double normal[4][4]{}, rhs[4]{};
        for (int i = 0; i < x.size(); ++i) {
            double jac[4], residual = y[i] - peak(x[i], p, s.fitModel, jac);
            for (int a = 0; a < 4; ++a) {
                rhs[a] += jac[a] * residual;
                for (int b = 0; b < 4; ++b)
                    normal[a][b] += jac[a] * jac[b];
            }
        }
        double gradient = 0;
        for (double v : rhs)
            gradient += v * v;
        if (gradient < 1e-20) {
            converged = true;
            break;
        }
        for (int a = 0; a < 4; ++a)
            normal[a][a] += lambda * (normal[a][a] + 1e-12);
        double delta[4];
        if (!solve4(normal, rhs, delta)) {
            lambda *= 10;
            continue;
        }
        auto candidate = p;
        for (int a = 0; a < 4; ++a)
            candidate[a] += delta[a];
        candidate[1] = std::clamp(candidate[1], -span, span);
        candidate[2] = std::clamp(candidate[2], std::log(span / 100000), std::log(span * 5));
        double trial = error(candidate);
        if (std::isfinite(trial) && trial < current) {
            double improvement = current - trial;
            p = candidate;
            current = trial;
            lambda = std::max(1e-12, lambda * 0.3);
            if (improvement < 1e-12 * (1 + current)) {
                converged = true;
                break;
            }
        } else {
            lambda *= 10;
            if (lambda > 1e18)
                break;
        }
    }
    result.converged = converged;
    result.centerKev = (origin + p[1]) / 1000;
    result.fwhmEv = std::exp(p[2]) * (s.fitModel == 0 ? 2 : 2.3548200450309493);
    result.amplitude = p[0] * scale;
    result.background = p[3] * scale;
    double mean = std::accumulate(y.begin(), y.end(), 0.0) / y.size(), variance = 0;
    for (double v : y)
        variance += (v - mean) * (v - mean);
    result.rSquared = variance > 0 ? 1 - current / variance : NAN;
    for (double v : x)
        result.fitted.append(peak(v, p, s.fitModel) * scale);
    if (!converged)
        result.reason = QStringLiteral("拟合未收敛。");
    else if (result.amplitude <= 0)
        result.reason = QStringLiteral("拟合峰高非正。");
    else if (result.centerKev < s.centerMin || result.centerKev > s.centerMax)
        result.reason = QStringLiteral("弹性峰中心超出允许范围。");
    else if (result.fwhmEv > s.maxFwhmEv)
        result.reason = QStringLiteral("FWHM 超过阈值。");
    else if (!std::isfinite(result.rSquared) || result.rSquared < s.minRSquared)
        result.reason = QStringLiteral("R² 低于阈值。");
    else
        result.accepted = true;
    return result;
}

ElasticCalibration Spectroscopy::calibrate(const SpectrumSettings &s, const CrystalRois &rois,
                                           const std::atomic_bool &cancel, Progress progress) {
    validateSettings(s);
    if (rois.count(RoiState::Pending))
        fail(QStringLiteral("请先完成 ROI：未设置的晶体应框选或标记失效。"));
    if (!rois.count(RoiState::Valid))
        fail(QStringLiteral("没有有效晶体 ROI。"));
    ElasticCalibration out;
    out.original = rois;
    out.adjusted = rois;
    auto first = loadScan(s.elasticFile, s, true, nullptr, true, cancel, progress);
    out.images = first.images;
    for (const auto &image : out.images)
        out.sizes.append(QSize(image.width, image.height));
    out.pixels.resize(75);
    out.filteredPixels.resize(75);
    out.fits.resize(75);
    QVector<QVector<int>> rawUse(2), filteredUse(2);
    for (int d = 0; d < 2; ++d) {
        rawUse[d].fill(0, out.images[d].values.size());
        filteredUse[d].fill(0, out.images[d].values.size());
    }
    for (int i = 0; i < 75; ++i) {
        checkCancel(cancel);
        const auto &roi = rois.items()[i];
        if (roi.state != RoiState::Valid) {
            out.fits[i].reason = QStringLiteral("人工标记失效。");
            continue;
        }
        int d = i < 60 ? 0 : 1;
        const auto &image = out.images[d];
        QRect bounds(0, 0, image.width, image.height), rect = roi.rect;
        if (rect.isEmpty() || !bounds.contains(rect))
            fail(QStringLiteral("ROI 超出当前探测器：%1").arg(roi.label));
        if (s.autoRoi) {
            double maximum = -INFINITY;
            QPoint location;
            for (int y = rect.top(); y <= rect.bottom(); ++y)
                for (int x = rect.left(); x <= rect.right(); ++x) {
                    double value = image.values[y * image.width + x];
                    if (value > maximum) {
                        maximum = value;
                        location = {x, y};
                    }
                }
            auto shift = s.shifts.value(i), expand = s.expansions.value(i);
            int x1 = int(std::nearbyint(location.x() - s.roiSize / 2.0)) + shift.x() - expand.x(),
                x2 = int(std::nearbyint(location.x() + s.roiSize / 2.0)) + shift.x() + expand.x();
            int y1 = int(std::nearbyint(location.y() - s.roiSize / 2.0)) + shift.y() - expand.y(),
                y2 = int(std::nearbyint(location.y() + s.roiSize / 2.0)) + shift.y() + expand.y();
            rect = QRect(x1, y1, x2 - x1, y2 - y1).intersected(bounds);
            if (!out.adjusted.assign(i, rect, out.sizes[d]))
                fail(QStringLiteral("自动调整后 ROI 为空：%1").arg(roi.label));
        }
        double maximum = -INFINITY;
        for (int y = rect.top(); y <= rect.bottom(); ++y)
            for (int x = rect.left(); x <= rect.right(); ++x)
                maximum = std::max(maximum, image.values[y * image.width + x]);
        for (int y = rect.top(); y <= rect.bottom(); ++y)
            for (int x = rect.left(); x <= rect.right(); ++x) {
                int pixel = y * image.width + x;
                out.pixels[i].append(pixel);
                ++rawUse[d][pixel];
                if (image.values[pixel] > maximum * s.filterFraction) {
                    out.filteredPixels[i].append(pixel);
                    ++filteredUse[d][pixel];
                }
            }
    }
    for (int d = 0; d < 2; ++d) {
        auto mask = imageData(detectorPath(d) + " mask", out.sizes[d]);
        for (int pixel = 0; pixel < mask.values.size(); ++pixel)
            if ((s.useFilter ? filteredUse[d][pixel] : rawUse[d][pixel]) > 1)
                ++out.overlapPixels;
        out.maskImages.append(mask);
    }
    for (int i = 0; i < 75; ++i) {
        int d = i < 60 ? 0 : 1;
        if (s.removeOverlap) {
            auto remove = [&](QVector<int> &pixels, const QVector<int> &use) {
                pixels.erase(
                    std::remove_if(pixels.begin(), pixels.end(), [&](int pixel) { return use[pixel] > 1; }),
                    pixels.end());
            };
            remove(out.pixels[i], rawUse[d]);
            remove(out.filteredPixels[i], filteredUse[d]);
        }
        for (int pixel : s.useFilter ? out.filteredPixels[i] : out.pixels[i])
            out.maskImages[d].values[pixel] += 1;
    }
    for (auto &image : out.maskImages)
        imageStats(image);
    out.scan = loadScan(s.elasticFile, s, true, &out, false, cancel, progress);
    for (int i = 0; i < 75; ++i) {
        checkCancel(cancel);
        if (out.adjusted.items()[i].state != RoiState::Valid)
            continue;
        if ((s.useFilter ? out.filteredPixels[i] : out.pixels[i]).isEmpty()) {
            out.fits[i].reason = QStringLiteral("筛选或剔除重叠后无有效像素。");
            continue;
        }
        out.fits[i] = fitPeak(out.scan.energyKev, s.useFilter ? out.scan.filtered[i] : out.scan.raw[i], s);
    }
    return out;
}
QVector<SpectrumScan> Spectroscopy::extract(const SpectrumSettings &s, const ElasticCalibration &calibration,
                                            const QStringList &files, const std::atomic_bool &cancel,
                                            Progress progress) {
    if (files.isEmpty())
        fail(QStringLiteral("请至少勾选一个 XRS 扫描。"));
    QVector<SpectrumScan> out;
    QSet<QString> seen;
    for (const auto &file : files) {
        auto path = QFileInfo(file).canonicalFilePath();
        if (seen.contains(path))
            fail(QStringLiteral("不能重复叠加同一 XRS 文件。"));
        seen.insert(path);
        out.append(loadScan(file, s, false, &calibration, true, cancel, progress));
    }
    return out;
}
QMap<QString, double> Spectroscopy::readAngles(const QString &file) {
    QMap<QString, double> result;
    for (const auto &module : {"VB", "HL", "VU", "HB", "VD"}) {
        QString key = QString("BL_ID33_XRS_m_%1").arg(module);
        for (const auto &path : {"/entry/instrument/" + key + "_RBV", "/entry/gui/md/" + key + ".RBV"}) {
            auto value = HdfFile::read({file, path, {}, -1, -1});
            if (value.error.isEmpty() && value.values.size() == 1 && std::isfinite(value.values[0])) {
                result[module] = value.values[0];
                break;
            }
        }
    }
    return result;
}
QPair<double, double> Spectroscopy::momentum(const QString &label, double angle, double incident,
                                             double final, const SpectrumSettings &s) {
    const QString module = label.section('-', 0, 0), crystal = label.section('-', 1, 1);
    if (crystal.size() != 2 || crystal[0] < 'A' || crystal[0] > 'E' || crystal[1] < '1' || crystal[1] > '3' ||
        incident <= 0 || final <= 0 || s.crystalDistance <= 0 || s.crystalDiameter <= 0)
        return {NAN, NAN};
    constexpr double pi = 3.14159265358979323846, hc = 12.398419297617678;
    double reverse = (module == "HB" || module == "HL") ? 1 : -1,
           theta = (angle + s.modulePitchDegrees * (crystal[0].unicode() - 'C') * reverse) * pi / 180,
           delta = s.modulePitchDegrees * (crystal[1].unicode() - '2') * pi / 180;
    double scattering = std::acos(std::clamp(std::cos(theta) * std::cos(delta), -1.0, 1.0)),
           ki = 2 * pi * incident / hc, kf = 2 * pi * final / hc;
    double q = std::sqrt(std::max(0.0, ki * ki + kf * kf - 2 * ki * kf * std::cos(scattering)));
    return {q,
            q > 1e-12 ? ki * kf * s.crystalDiameter * std::sin(scattering) / (q * s.crystalDistance) : NAN};
}
SpectrumResult Spectroscopy::combine(const SpectrumSettings &s, const ElasticCalibration &calibration,
                                     const QVector<SpectrumScan> &scans, const std::atomic_bool &cancel) {
    if (scans.isEmpty())
        fail(QStringLiteral("请先提取 XRS 谱。"));
    if (!(s.energyStepEv > 0) || !std::isfinite(s.energyStepEv) || s.qMax < s.qMin)
        fail(QStringLiteral("能量步长或 Q 范围无效。"));
    SpectrumResult out;
    out.angles =
        s.angleSource == 2 ? s.angles : readAngles(s.angleSource == 0 ? s.elasticFile : scans.first().file);
    out.momentum.resize(75);
    out.included.fill(false, 75);
    out.reasons.resize(75);
    out.crystals.resize(75);
    double lower = -INFINITY, upper = INFINITY;
    int selected = 0;
    for (int r = 0; r < 75; ++r) {
        checkCancel(cancel);
        const auto &roi = calibration.adjusted.items()[r];
        const auto &fit = calibration.fits[r];
        auto &reason = out.reasons[r];
        if (roi.state != RoiState::Valid) {
            reason = QStringLiteral("人工失效");
            continue;
        }
        if (!fit.accepted) {
            reason = fit.reason;
            continue;
        }
        QString module = roi.label.section('-', 0, 0);
        if (!out.angles.contains(module))
            fail(QStringLiteral("缺少模组角度 %1，请指定角度来源或手动填写。").arg(module));
        long double sum = 0, dsum = 0;
        int count = 0;
        double min = INFINITY, max = -INFINITY, dmin = INFINITY, dmax = -INFINITY;
        for (const auto &scan : scans)
            for (double energy : scan.energyKev) {
                auto q = momentum(roi.label, out.angles[module], energy, fit.centerKev, s);
                if (!std::isfinite(q.first) || !std::isfinite(q.second))
                    fail(QStringLiteral("Q/dQ 计算失败：%1，检查角度和几何参数。").arg(roi.label));
                sum += q.first;
                dsum += q.second;
                min = std::min(min, q.first);
                max = std::max(max, q.first);
                dmin = std::min(dmin, q.second);
                dmax = std::max(dmax, q.second);
                ++count;
            }
        auto &q = out.momentum[r];
        q = {double(sum / count), double(dsum / count), max - min, dmax - dmin};
        if (s.excludedCrystals.contains(r)) {
            reason = QStringLiteral("手动排除");
            continue;
        }
        if (!s.modules.contains(module)) {
            reason = QStringLiteral("未选模组");
            continue;
        }
        if (q.mean < s.qMin || q.mean > s.qMax) {
            reason = QStringLiteral("Q 超出范围");
            continue;
        }
        out.included[r] = true;
        ++selected;
        for (const auto &scan : scans) {
            double a = (scan.energyKev.first() - fit.centerKev) * 1000,
                   b = (scan.energyKev.last() - fit.centerKev) * 1000;
            lower = std::max(lower, std::min(a, b));
            upper = std::min(upper, std::max(a, b));
        }
    }
    if (!selected)
        fail(QStringLiteral("没有通过 ROI、拟合质量、模组和 Q 筛选的晶体。"));
    double first = std::ceil(lower / s.energyStepEv - 1e-10) * s.energyStepEv,
           last = std::floor(upper / s.energyStepEv + 1e-10) * s.energyStepEv;
    if (!(last >= first))
        fail(QStringLiteral("所选谱没有共同能量转移范围。"));
    double gridSize = std::floor((last - first) / s.energyStepEv + 1e-8) + 1;
    if (!std::isfinite(gridSize) || gridSize > 200000)
        fail(QStringLiteral("插值点超过 200000，请增大步长。"));
    int n = int(gridSize);
    if (n < 2)
        fail(QStringLiteral("共同范围内不足两个插值点。"));
    for (int i = 0; i < n; ++i)
        out.energyTransferEv.append(first + i * s.energyStepEv);
    out.total.fill(0, n);
    for (const auto &m : s.modules)
        out.modules[m].fill(0, n);
    for (int r = 0; r < 75; ++r)
        if (out.included[r]) {
            checkCancel(cancel);
            auto &curve = out.crystals[r];
            curve.fill(0, n);
            for (const auto &scan : scans) {
                QVector<double> x, y = s.useFilter ? scan.filtered[r] : scan.raw[r];
                for (double energy : scan.energyKev)
                    x.append((energy - calibration.fits[r].centerKev) * 1000);
                if (x.first() > x.last()) {
                    std::reverse(x.begin(), x.end());
                    std::reverse(y.begin(), y.end());
                }
                for (int i = 0; i < n; ++i)
                    curve[i] += interpolate(x, y, out.energyTransferEv[i]);
            }
            QString module = calibration.adjusted.items()[r].label.section('-', 0, 0);
            for (int i = 0; i < n; ++i) {
                out.total[i] += curve[i];
                out.modules[module][i] += curve[i];
            }
        }
    return out;
}

bool Spectroscopy::saveResult(const QString &directory, const SpectrumSettings &s,
                              const ElasticCalibration &c, const QVector<SpectrumScan> &scans,
                              const SpectrumResult &result, QString *error) {
    QString staging;
    try {
        if (result.total.isEmpty())
            fail(QStringLiteral("没有可保存的叠加结果。"));
        if (QFileInfo::exists(directory))
            fail(QStringLiteral("输出目录已存在，请选择新名称，避免覆盖。"));
        staging = directory + ".partial-" + QUuid::createUuid().toString(QUuid::Id128);
        if (!QDir().mkpath(staging))
            fail(QStringLiteral("不能创建输出目录。"));
        auto write = [&](const QString &name, const QByteArray &data) {
            QSaveFile file(staging + "/" + name);
            if (!file.open(QIODevice::WriteOnly) || file.write(data) != data.size() || !file.commit())
                fail(QStringLiteral("写入失败：%1").arg(name));
        };
        auto table = [&](const QStringList &headers, const QVector<QVector<double>> &columns) {
            QString text;
            QTextStream stream(&text);
            stream.setRealNumberPrecision(17);
            stream << headers.join('\t') << '\n';
            for (int row = 0; row < result.energyTransferEv.size(); ++row) {
                for (int col = 0; col < columns.size(); ++col) {
                    if (col)
                        stream << '\t';
                    stream << columns[col].value(row, NAN);
                }
                stream << '\n';
            }
            return text.toUtf8();
        };
        QStringList summaryHeaders{"Energy Transfer (eV)", "Intensity"};
        QVector<QVector<double>> summary{result.energyTransferEv, result.total};
        for (auto it = result.modules.begin(); it != result.modules.end(); ++it) {
            summaryHeaders << it.key();
            summary << it.value();
        }
        write("spectrum_data.txt", table(summaryHeaders, summary));
        QStringList allHeaders{"Energy Transfer (eV)"};
        QVector<QVector<double>> all{result.energyTransferEv};
        QString roiText;
        QTextStream rois(&roiText);
        rois.setRealNumberPrecision(17);
        rois << "crystal\tdetector\tx1\tx2\ty1\ty2\tstate\tpixels\tfiltered_pixels\tcenter_eV\tfwhm_"
                "eV\theight\tbackground\tr_squared\tfit_ok\tq_ave_A^-1\tdq_ave_A^-1\tq_range_A^-1\tdq_range_"
                "A^-1\tadd\treason\n";
        for (int r = 0; r < 75; ++r) {
            const auto &roi = c.adjusted.items()[r];
            const auto &fit = c.fits[r];
            const auto &q = result.momentum[r];
            rois << roi.label << '\t' << roi.detector << '\t' << roi.rect.x() << '\t'
                 << roi.rect.x() + roi.rect.width() << '\t' << roi.rect.y() << '\t'
                 << roi.rect.y() + roi.rect.height() << '\t' << int(roi.state) << '\t' << c.pixels[r].size()
                 << '\t' << c.filteredPixels[r].size() << '\t' << fit.centerKev * 1000 << '\t' << fit.fwhmEv
                 << '\t' << fit.amplitude << '\t' << fit.background << '\t' << fit.rSquared << '\t'
                 << fit.accepted << '\t' << q.mean << '\t' << q.dMean << '\t' << q.range << '\t' << q.dRange
                 << '\t' << result.included[r] << '\t' << result.reasons[r] << '\n';
            if (result.included[r]) {
                allHeaders << roi.label;
                all << result.crystals[r];
                write(roi.label + "_data.txt", table({"Energy Transfer (eV)", "Intensity"},
                                                     {result.energyTransferEv, result.crystals[r]}));
            }
        }
        write("spectrum_all_data.txt", table(allHeaders, all));
        write("spectrum_rois.txt", roiText.toUtf8());
        write("original_roi.roi.json", c.original.serialize(s.elasticFile, c.sizes));
        write("adjusted_roi.roi.json", c.adjusted.serialize(s.elasticFile, c.sizes));
        QJsonArray records;
        for (const auto &scan : scans) {
            QJsonArray glitches;
            for (int i : scan.glitches)
                glitches.append(i);
            records.append(QJsonObject{{"file", scan.file},
                                       {"energy_path", scan.energyPath},
                                       {"i0_path", scan.i0Path},
                                       {"glitches", glitches},
                                       {"energy_keV", numbers(scan.energyKev)},
                                       {"i0_before", numbers(scan.i0Before)},
                                       {"i0_after", numbers(scan.i0After)}});
        }
        QJsonArray maskRecords;
        for (int r = 0; r < 75; ++r) {
            QJsonArray raw, filtered;
            for (int p : c.pixels[r])
                raw.append(p);
            for (int p : c.filteredPixels[r])
                filtered.append(p);
            maskRecords.append(QJsonObject{
                {"crystal", c.adjusted.items()[r].label}, {"pixels", raw}, {"filtered_pixels", filtered}});
        }
        write("spectrum_masks.json", QJsonDocument(maskRecords).toJson());
        QJsonObject angles;
        for (auto it = result.angles.begin(); it != result.angles.end(); ++it)
            angles[it.key()] = it.value();
        QJsonArray modules;
        for (const auto &module : s.modules)
            modules.append(module);
        QJsonArray exclusions, corrections;
        for (int r = 0; r < 75; ++r) {
            if (s.excludedCrystals.contains(r))
                exclusions.append(c.adjusted.items()[r].label);
            corrections.append(QJsonObject{{"crystal", c.adjusted.items()[r].label},
                                           {"x_shift", s.shifts.value(r).x()},
                                           {"y_shift", s.shifts.value(r).y()},
                                           {"x_expand", s.expansions.value(r).x()},
                                           {"y_expand", s.expansions.value(r).y()}});
        }
        QJsonObject info{{"version", 1},
                         {"created", QDateTime::currentDateTimeUtc().toString(Qt::ISODate)},
                         {"elastic_file", s.elasticFile},
                         {"elastic_energy_path", c.scan.energyPath},
                         {"elastic_i0_path", c.scan.i0Path},
                         {"normalize_elastic", s.normalizeElastic},
                         {"normalize_xrs", s.normalizeXrs},
                         {"repair_i0_glitches", s.repairGlitches},
                         {"glitch_ratio", 0.7},
                         {"use_filter", s.useFilter},
                         {"filter_fraction", s.filterFraction},
                         {"auto_roi", s.autoRoi},
                         {"roi_size", s.roiSize},
                         {"corrections", corrections},
                         {"remove_overlap", s.removeOverlap},
                         {"overlap_pixels", c.overlapPixels},
                         {"fit_model", s.fitModel == 0 ? "Lorentzian" : "Gaussian"},
                         {"center_min_keV", s.centerMin},
                         {"center_max_keV", s.centerMax},
                         {"max_fwhm_eV", s.maxFwhmEv},
                         {"min_r_squared", s.minRSquared},
                         {"module_angles_degrees", angles},
                         {"angle_source", s.angleSource},
                         {"crystal_diameter_mm", s.crystalDiameter},
                         {"crystal_distance_mm", s.crystalDistance},
                         {"module_pitch_degrees", s.modulePitchDegrees},
                         {"modules", modules},
                         {"excluded_crystals", exclusions},
                         {"q_min_A^-1", s.qMin},
                         {"q_max_A^-1", s.qMax},
                         {"energy_step_eV", s.energyStepEv},
                         {"combination", "linear interpolation on common range; sum over scans and crystals, "
                                         "no averaging or extrapolation"},
                         {"scans", records}};
        auto infoBytes = QJsonDocument(info).toJson();
        write("spectrum_info.json", infoBytes);
        write("spectrum_info.txt", infoBytes);
        if (!QDir().rename(staging, directory))
            fail(QStringLiteral("不能完成输出目录重命名。"));
        return true;
    } catch (const std::exception &e) {
        if (error)
            *error = QString::fromUtf8(e.what()) +
                     (staging.isEmpty() ? QString() : QStringLiteral("\n未完成输出目录：") + staging);
        return false;
    }
}
