#include "CrystalRoi.h"
#include <QCoreApplication>
#include "ImagingIO.h"
#include <QDir>
#include <QFile>
#include <QFileInfo>
#include <QJsonDocument>
#include <QRegularExpression>
#include <QSet>
#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>
using namespace ImageUtil;
namespace {
const QStringList calibrationMotors{"M_Stage_X", "M_DCM_B5Energy"};
const QStringList sampleMotors{"M_Stage_Y", "M_Stage_Z", "M_DCM_B5Energy"};
constexpr double energyCoverageTolerance = 1e-10;

struct BlankScan {
    QVector<double> axis, dose;
    QVector<int> frameBins;
};

QPair<qsizetype, qsizetype> coveredSamples(const QVector<double> &sample, double low, double high) {
    return {std::lower_bound(sample.cbegin(), sample.cend(), low - energyCoverageTolerance) - sample.cbegin(),
            std::upper_bound(sample.cbegin(), sample.cend(), high + energyCoverageTolerance) - sample.cbegin()};
}

QVector<double> interpolateWeights(const QVector<double> &source, const QVector<double> &target,
                                   const QVector<double> &weights) {
    QVector<double> result(source.size(), 0.);
    for (qsizetype k = 0; k < target.size(); ++k) {
        const double energy = std::clamp(target[k], source.first(), source.last());
        const auto it = std::upper_bound(source.cbegin(), source.cend(), energy);
        const auto j = std::clamp(it - source.cbegin() - 1, qsizetype(0), source.size() - 2);
        const double fraction = (energy - source[j]) / (source[j + 1] - source[j]);
        result[j] += weights[k] * (1 - fraction);
        result[j + 1] += weights[k] * fraction;
    }
    return result;
}

BlankScan blankScan(const QString &path, const ImagingGrid &sample, const QJsonObject &acq,
                    double exposure, const std::atomic_bool &stop) {
    cancel(stop);
    if (!QFileInfo(path).isFile() ||
        QFileInfo(path).canonicalFilePath() == QFileInfo(sample.file).canonicalFilePath())
        fail("Select an existing blank scan independent of the sample");
    ImagingH5 source(path);
    const auto energy = source.read("/entry/instrument/M_DCM_B5Energy_setpoint");
    const auto monitor = source.read("/entry/data/" + acq["monitor"].toString());
    const auto detector = HdfFile::inspect(path, "/entry/data/" + acq["detector"].toString());
    if (energy.isEmpty() || energy.size() != monitor.size() || !detector.type.startsWith("uint") ||
        detector.dims != QVector<quint64>{quint64(energy.size()), quint64(sample.imageSize.height()),
                                          quint64(sample.imageSize.width())})
        fail("Blank detector geometry, monitor and energy records must match the acquired frames");
    BlankScan blank;
    blank.axis = energy;
    for (double e : energy)
        if (!std::isfinite(e))
            fail("Nonfinite blank energy");
    std::sort(blank.axis.begin(), blank.axis.end());
    blank.axis.erase(std::unique(blank.axis.begin(), blank.axis.end(), [](double a, double b) {
                         return std::abs(a - b) <= energyCoverageTolerance;
                     }),
                     blank.axis.end());
    if (blank.axis.size() < 2)
        fail(QCoreApplication::translate("ImagingConfig", "空白扫描至少需要两个不同的能量点。"));
    blank.dose.fill(0., blank.axis.size());
    blank.frameBins.reserve(energy.size());
    for (qsizetype f = 0; f < energy.size(); ++f) {
        if (!(monitor[f] > 0. && std::isfinite(monitor[f])))
            fail("Invalid blank monitor");
        const auto bin = std::lower_bound(blank.axis.cbegin(), blank.axis.cend(), energy[f] - energyCoverageTolerance);
        if (bin == blank.axis.cend())
            fail("Cannot match blank energy");
        const int e = int(bin - blank.axis.cbegin());
        blank.frameBins.append(e);
        blank.dose[e] += monitor[f] * exposure;
    }
    return blank;
}

void selectEnergyBand(ImagingValidation &validation, const BlankScan &blank, bool matchFoil) {
    const auto &sample = validation.sample.axes[2];
    const auto &foil = validation.calibration.axes[1];
    qsizetype first = 0, end = sample.size();
    if (matchFoil) {
        double low = foil.first(), high = foil.last();
        if (!blank.axis.isEmpty()) {
            low = std::max(low, blank.axis.first());
            high = std::min(high, blank.axis.last());
        }
        const auto bounds = coveredSamples(sample, low, high);
        first = bounds.first;
        end = bounds.second;
        if (end - first < 2)
            fail(QCoreApplication::translate(
                "ImagingConfig", "共同覆盖的样品能量点不足两个。可关闭能区匹配，用独立薄膜空间响应重建。"));
    }

    auto &band = validation.energyBand;
    band.energies = sample.mid(first, end - first);
    const auto usedWeights = Imaging::energyWeights(band.energies);
    band.sampleWeights.fill(0., sample.size());
    std::copy(usedWeights.cbegin(), usedWeights.cend(), band.sampleWeights.begin() + first);
    if (matchFoil) {
        band.calibrationWeights = interpolateWeights(foil, band.energies, usedWeights);
    } else {
        band.calibrationWeights = Imaging::energyWeights(foil);
        // Keep the operator on the sample integration scale without evaluating the foil at sample energies.
        const double scale = (sample.last() - sample.first()) / (foil.last() - foil.first());
        for (double &weight : band.calibrationWeights)
            weight *= scale;
    }

    const auto range = [](const QVector<double> &axis) { return QJsonArray{axis.first(), axis.last()}; };
    QJsonArray discarded;
    for (qsizetype e = 0; e < first; ++e)
        discarded.append(sample[e]);
    for (qsizetype e = end; e < sample.size(); ++e)
        discarded.append(sample[e]);
    const auto &frameIndices = validation.sample.indices[2];
    const auto usedFrames = std::count_if(frameIndices.cbegin(), frameIndices.cend(),
                                          [&](int e) { return band.sampleWeights[e] > 0.; });
    QJsonObject report{{"calibration_mode", matchFoil ? "sample-matched" : "foil-native-spatial"},
                       {"original_sample_range_kev", range(sample)},
                       {"calibration_range_kev", range(foil)},
                       {"used_range_kev", range(band.energies)},
                       {"used_sample_energies_kev", json(band.energies)},
                       {"discarded_sample_energies_kev", discarded},
                       {"used_sample_frames", qsizetype(usedFrames)},
                       {"discarded_sample_frames", frameIndices.size() - qsizetype(usedFrames)}};
    if (!blank.axis.isEmpty()) {
        const auto [lo, hi] = coveredSamples(band.energies, blank.axis.first(), blank.axis.last());
        const auto covered = band.energies.mid(lo, hi - lo);
        report["blank_corrected_sample_energies_kev"] = json(covered);
        if (!covered.isEmpty()) {
            const auto blankWeights = interpolateWeights(blank.axis, covered, usedWeights.mid(lo, hi - lo));
            band.blankFrameWeights.reserve(blank.frameBins.size());
            for (int bin : blank.frameBins)
                band.blankFrameWeights.append(blankWeights[bin] / blank.dose[bin]);
        }
        report["blank_range_kev"] = range(blank.axis);
    }
    if (!matchFoil) {
        band.notice = QCoreApplication::translate(
            "ImagingConfig", "使用薄膜自身能区确定空间响应，积分全部 %1 个样品能量点。").arg(sample.size());
    } else if (!discarded.isEmpty()) {
        band.notice = QCoreApplication::translate(
                          "ImagingConfig", "按共同覆盖能区积分：%1–%2 keV，保留 %3/%4 个样品能量点，舍弃 %5 帧。")
                          .arg(QString::number(band.energies.first(), 'g', 12),
                               QString::number(band.energies.last(), 'g', 12))
                          .arg(band.energies.size()).arg(sample.size())
                          .arg(frameIndices.size() - qsizetype(usedFrames));
    }
    if (!blank.axis.isEmpty() && report["blank_corrected_sample_energies_kev"].toArray().size() < band.energies.size())
        band.notice += ' ' + QCoreApplication::translate("ImagingConfig", "空白扣除仅覆盖其中 %1 个样品能量点。")
                                 .arg(report["blank_corrected_sample_energies_kev"].toArray().size());
    if (!validation.verification.file.isEmpty()) {
        const auto &axis = validation.verification.axes[1];
        const auto [lo, hi] = coveredSamples(band.energies, axis.first(), axis.last());
        band.verificationEnergies = band.energies.mid(lo, hi - lo);
        report["verification_foil_range_kev"] = range(axis);
        report["verification_sample_energies_kev"] = json(band.verificationEnergies);
        if (!band.verificationEnergies.isEmpty())
            band.verificationWeights = Imaging::matchedWeights(axis, band.verificationEnergies).first;
    }
    validation.report["energy_integration"] = report;
}

QString loopAxis(const QString &motor) {
    return motor == "M_DCM_B5Energy" ? QStringLiteral("DCM") : motor.right(1);
}
QString loopCountKey(const QString &role, const QString &axis) {
    return role + '_' + axis.toLower() + "_points";
}

QJsonObject merge(const QJsonObject &defaults, const QJsonObject &input, const QString &prefix = {}) {
    auto result = defaults;
    for (auto i = input.begin(); i != input.end(); ++i) {
        if (!defaults.contains(i.key()))
            fail("Unknown imaging field: " + prefix + i.key());
        auto value = defaults[i.key()];
        if (value.isObject() && !value.toObject().isEmpty()) {
            if (!i->isObject())
                fail(prefix + i.key() + " must be an object");
            result[i.key()] = merge(value.toObject(), i->toObject(), prefix + i.key() + ".");
        } else
            result[i.key()] = i.value();
    }
    return result;
}
// Safe Python-literal subset for the existing P.grid_scan tuple. No eval/code execution.
class Literal {
    QString s;
    qsizetype p = 0;
    void white() {
        while (p < s.size() && s[p].isSpace())
            ++p;
    }
    bool take(QChar c) {
        white();
        if (p < s.size() && s[p] == c) {
            ++p;
            return true;
        }
        return false;
    }
    QJsonValue parse(int depth) {
        white();
        if (depth > 64 || p >= s.size())
            fail("Invalid scan plan literal");
        QChar c = s[p];
        if (c == '\'' || c == '"') {
            ++p;
            QString out;
            while (p < s.size()) {
                QChar v = s[p++];
                if (v == c)
                    return out;
                if (v == '\\') {
                    if (p >= s.size())
                        break;
                    v = s[p++];
                    if (v == 'n')
                        v = '\n';
                    else if (v == 'r')
                        v = '\r';
                    else if (v == 't')
                        v = '\t';
                }
                out += v;
            }
            fail("Unterminated plan string");
        }
        if (c == '(' || c == '[') {
            ++p;
            QChar end = c == '(' ? ')' : ']';
            QJsonArray out;
            if (take(end))
                return out;
            for (;;) {
                out.append(parse(depth + 1));
                if (take(end))
                    return out;
                if (!take(','))
                    fail("Invalid plan sequence");
                if (take(end))
                    return out;
            }
        }
        if (c == '{') {
            ++p;
            QJsonObject out;
            if (take('}'))
                return out;
            for (;;) {
                auto key = parse(depth + 1);
                if (!key.isString() || !take(':'))
                    fail("Invalid plan dictionary");
                out[key.toString()] = parse(depth + 1);
                if (take('}'))
                    return out;
                if (!take(','))
                    fail("Invalid plan dictionary separator");
                if (take('}'))
                    return out;
            }
        }
        for (const auto &word : {QString("None"), QString("True"), QString("False")})
            if (s.mid(p, word.size()) == word) {
                p += word.size();
                return word == "None" ? QJsonValue() : QJsonValue(word == "True");
            }
        static QRegularExpression numeric("^[+-]?(?:[0-9]+(?:\\.[0-9]*)?|\\.[0-9]+)(?:[eE][+-]?[0-9]+)?");
        auto match = numeric.match(s.mid(p));
        if (!match.hasMatch())
            fail("Non-literal scan plan token");
        p += match.capturedLength();
        bool ok;
        double value = match.captured().toDouble(&ok);
        if (!ok || !std::isfinite(value))
            fail("Nonfinite scan plan coordinate");
        return value;
    }

  public:
    explicit Literal(QString text) : s(std::move(text)) {}
    QJsonValue value() {
        auto out = parse(0);
        white();
        if (p != s.size())
            fail("Unexpected text after scan plan");
        return out;
    }
};
ImagingGrid grid(const QString &file, const QString &role, const QStringList &motors, double exposure,
                 const QJsonObject &acq, const QJsonObject &loops, const std::atomic_bool &stop) {
    cancel(stop);
    ImagingGrid result;
    result.file = file;
    ImagingH5 source(file);
    result.plan = source.text("/entry/gui/plan_cmd");
    auto plan = Literal(result.plan).value().toArray();
    if (plan.size() < 2 || plan[0].toString() != "P.grid_scan" || !plan[1].isArray())
        fail("Expected P.grid_scan metadata");
    auto args = plan[1].toArray();
    if (args.size() < 1 || (args.size() - 1) % 4)
        fail(QCoreApplication::translate("ImagingConfig", "扫描计划中的电机参数不完整：%1").arg(file));
    QStringList actual;
    for (int i = 1; i < args.size(); i += 4) {
        auto name = args[i].toString();
        name.replace('.', '_');
        actual.append(name);
    }
    auto roleLabel = role == "calibration" ? QCoreApplication::translate("ImagingConfig", "薄膜标定") : QCoreApplication::translate("ImagingConfig", "样品成像");
    auto sortedActual = actual;
    auto sortedMotors = motors;
    sortedActual.sort();
    sortedMotors.sort();
    if (sortedActual != sortedMotors) {
        auto hint = actual.contains("M_DCM_B5Energy")
                        ? QCoreApplication::translate("ImagingConfig", "请选择扫描轴与用途匹配的文件。")
                        : QCoreApplication::translate("ImagingConfig", "该扫描计划没有能量轴，不能用于当前需要能量积分的处理流程。");
        fail(QCoreApplication::translate("ImagingConfig", "%1文件的扫描轴不匹配：\n%2\n实际：%3\n所需轴：%4\n%5")
                 .arg(roleLabel, file, actual.join(", "), motors.join(", "), hint));
    }
    for (const auto &motor : actual)
        result.loopOrder.append(loopAxis(motor));
    const auto requestedOrder = loops[role + "_order"].toString();
    if (!requestedOrder.isEmpty() && requestedOrder != result.loopOrder.join(',')) {
        fail(QCoreApplication::translate("ImagingConfig", "%1文件的循环顺序与设置不符。\n%2\n文件（外到内）：%3\n设置（外到内）：%4")
                 .arg(roleLabel, file, result.loopOrder.join(" → "), requestedOrder.split(',').join(" → ")));
    }
    qsizetype cells = 1;
    int frames = -1;
    // Downstream arrays keep their physical X/E or Y/Z/E axes for every scan order.
    for (const auto &motor : motors) {
        auto start = 1 + 4 * actual.indexOf(motor);
        const auto axisName = loopAxis(motor);
        const auto requestedCount = loops[loopCountKey(role, axisName)];
        double lo = args[start + 1].toDouble(NAN), hi = args[start + 2].toDouble(NAN),
               n = requestedCount.isNull() ? args[start + 3].toDouble(NAN) : requestedCount.toDouble();
        if (!std::isfinite(lo) || !std::isfinite(hi) || !(hi > lo) || n < 2 || n > 100000 ||
            std::floor(n) != n)
            fail(QCoreApplication::translate("ImagingConfig", "%1的 %2 轴需要递增的扫描范围及至少两个点。").arg(roleLabel, axisName));
        if (cells > 10000000 / int(n))
            fail("Imaging grid too large");
        cells *= qsizetype(n);
        QVector<double> axis;
        for (int i = 0; i < int(n); ++i)
            axis.append(i == int(n) - 1 ? hi : lo + (hi - lo) * i / (n - 1));
        result.axes.append(axis);
        QVector<quint64> shape;
        auto coordinates = source.read("/entry/instrument/" + motor + "_setpoint", &shape);
        if (shape.size() != 1 || coordinates.isEmpty() || (frames >= 0 && coordinates.size() != frames))
            fail("Frame coordinates must be matching vectors");
        frames = int(coordinates.size());
        result.coordinates.append(coordinates);
        QVector<int> index;
        double tolerance = axisName == "DCM" ? 1e-8 : 1e-6;
        for (double c : coordinates) {
            if (!std::isfinite(c))
                fail("Nonfinite frame coordinate");
            auto it = std::lower_bound(axis.begin(), axis.end(), c);
            int j = it == axis.end() ? int(axis.size() - 1) : int(it - axis.begin());
            if (j > 0 && std::abs(c - axis[j - 1]) <= std::abs(c - axis[j]))
                --j;
            if (std::abs(c - axis[j]) > tolerance)
                fail(QCoreApplication::translate("ImagingConfig", "%1的 %2 轴逐帧坐标与 %3 个扫描点不匹配。请核对扫描范围和该层帧数。")
                         .arg(roleLabel, axisName)
                         .arg(int(n)));
            index.append(j);
        }
        result.indices.append(index);
    }
    result.occupancy.fill(0, cells);
    result.dose.fill(0, cells);
    QVector<quint64> monitorShape;
    result.monitor = source.read("/entry/data/" + acq["monitor"].toString(), &monitorShape);
    if (monitorShape.size() != 1 || result.monitor.size() != frames)
        fail("I0 must match detector frames");
    for (int f = 0; f < frames; ++f) {
        double i0 = result.monitor[f];
        if (!std::isfinite(i0) || i0 <= 0)
            fail("I0 must be positive and finite; imaging does not use abs(I0)");
        qsizetype index = 0;
        for (int k = 0; k < motors.size(); ++k)
            index = index * result.axes[k].size() + result.indices[k][f];
        ++result.occupancy[index];
        result.dose[index] += i0 * exposure;
    }
    if (result.occupancy.contains(0))
        fail(QCoreApplication::translate("ImagingConfig", "%1 扫描不完整：%2 个计划网格单元缺测；不以零填充。")
                 .arg(motors.size() == 2 ? QCoreApplication::translate("ImagingConfig", "薄膜")
                                         : QCoreApplication::translate("ImagingConfig", "样品"))
                 .arg(result.occupancy.count(0)));
    auto node = HdfFile::inspect(file, "/entry/data/" + acq["detector"].toString());
    if (node.dims.size() != 3 || node.dims[0] != quint64(frames) || !node.type.startsWith("uint"))
        fail("Expected unsigned photon counts [frame,row,column] matching motor records");
    if (!node.dims[1] || !node.dims[2] || node.dims[1] > 4096 || node.dims[2] > 4096)
        fail("Unsupported imaging detector dimensions");
    result.imageSize = {int(node.dims[2]), int(node.dims[1])};
    return result;
}
QVector<ImagingRoi> roiTable(const QString &file, const QString &detector, const QSize &size,
                             const QJsonObject &settings) {
    QFile input(file);
    if (!input.open(QIODevice::ReadOnly))
        fail("Cannot open ROI file: " + file);
    QByteArray bytes = input.readAll();
    if (bytes.startsWith("\xEF\xBB\xBF"))
        bytes.remove(0, 3);
    QSet<QString> excluded;
    for (auto value : settings["excluded_ids"].toArray())
        excluded.insert(value.toString());
    auto labels = settings["label_overrides"].toObject();
    QVector<ImagingRoi> result;
    if (bytes.trimmed().startsWith('{')) {
        auto document = QJsonDocument::fromJson(bytes).object();
        QVector<QSize> geometry;
        for (auto value : document["detector_sizes"].toArray()) {
            auto a = value.toArray();
            if (a.size() != 2)
                fail("Invalid shared ROI dimensions");
            geometry.append({a[0].toInt(), a[1].toInt()});
        }
        CrystalRois rois;
        QString error;
        if (!rois.restore(bytes, geometry, &error))
            fail(error);
        int d = detector == "D_LAMBDA" ? 0 : detector == "D_MINIPIX" ? 1 : -1;
        if (d < 0 || geometry.size() != 2 || geometry[d] != size)
            fail("Shared ROI scheme does not match chosen imaging detector");
        for (int i = 0; i < 75; ++i) {
            const auto &r = rois.items()[i];
            if (r.detector != detector)
                continue;
            if (r.state == RoiState::Pending)
                fail("Shared ROI still contains pending crystals for the selected detector");
            auto id = QString::number(i), label = labels.value(id).toString(r.label);
            result.append({id,
                           label,
                           r.label,
                           id + "_" + label,
                           detector,
                           r.rect,
                           r.state == RoiState::Invalid || excluded.contains(id),
                           r.mask,
                           r.state == RoiState::Invalid});
        }
    } else {
        auto lines = QString::fromUtf8(bytes).split('\n');
        if (lines.size() < 2)
            fail("Empty ROI table");
        auto header = lines.takeFirst().trimmed().split('\t'); // Preserve the unlabeled source-ID column.
        if (!header.isEmpty() && header[0] == "roi_label")
            header.prepend("");
        for (const auto &line : lines) {
            if (line.trimmed().isEmpty())
                continue;
            auto text = line;
            if (text.endsWith('\r'))
                text.chop(1);
            auto fields = text.split('\t');
            if (fields.size() != header.size())
                fail("ROI TSV column count mismatch");
            QMap<QString, QString> row;
            for (int c = 0; c < header.size(); ++c)
                row[header[c]] = fields[c].trimmed();
            auto id = row.value("", QString::number(result.size()));
            if (id.isEmpty())
                id = QString::number(result.size());
            auto integer = [&](const QString &key, bool optional = false) {
                if (optional && row.value(key).isEmpty())
                    return 0;
                bool ok;
                int n = row.value(key).toInt(&ok);
                if (!ok || std::abs(double(n)) > 10000000)
                    fail("Invalid ROI integer: " + key);
                return n;
            };
            int dx = integer("x_shift", true), dy = integer("y_shift", true), ex = integer("x_expand", true),
                ey = integer("y_expand", true);
            int x0 = integer("x1") + dx - ex, x1 = integer("x2") + dx + ex, y0 = integer("y1") + dy - ey,
                y1 = integer("y2") + dy + ey;
            if (x0 < 0 || y0 < 0 || x1 <= x0 || y1 <= y0 || x1 > size.width() || y1 > size.height())
                fail("TSV ROI is outside detector: " + id);
            auto original = row["roi_label"], label = labels.value(id).toString(original);
            result.append({id,
                           label,
                           original,
                           id + "_" + label,
                           detector,
                           {x0, y0, x1 - x0, y1 - y0},
                           excluded.contains(id),
                           {}});
        }
    }
    QSet<QString> ids, keys;
    QVector<int> coverage(size.width() * size.height(), 0);
    int active = 0;
    for (auto &r : result) {
        bool integerId = false;
        r.id.toLongLong(&integerId);
        if (!integerId)
            fail("Imaging ROI source IDs must be integers: " + r.id);
        if (ids.contains(r.id) || keys.contains(r.key) || r.label.isEmpty() || r.label.contains('/'))
            fail("Duplicate ROI ID/key or invalid label");
        ids.insert(r.id);
        keys.insert(r.key);
        if (!r.rect.isEmpty() && !QRect(QPoint(), size).contains(r.rect))
            fail("ROI outside detector: " + r.key);
        if (!r.excluded && r.rect.isEmpty())
            fail("Empty valid ROI: " + r.key);
        if (r.pixels.isEmpty())
            for (int y = r.rect.top(); !r.rect.isEmpty() && y <= r.rect.bottom(); ++y)
                for (int x = r.rect.left(); x <= r.rect.right(); ++x)
                    r.pixels.append({x, y});
        if (!r.excluded) {
            for (const auto &p : r.pixels)
                if (++coverage[p.y() * size.width() + p.x()] > 1)
                    fail("Valid imaging ROIs overlap; shared pixels must not be treated independently");
            ++active;
        }
    }
    if (!active)
        fail("No usable imaging ROI remains");
    for (const auto &id : excluded)
        if (!ids.contains(id))
            fail("Excluded ROI ID not in selected detector scheme: " + id);
    for (auto i = labels.begin(); i != labels.end(); ++i)
        if (!ids.contains(i.key()))
            fail("ROI label override ID not found: " + i.key());
    std::sort(result.begin(), result.end(),
              [](const ImagingRoi &a, const ImagingRoi &b) { return a.id.toLongLong() < b.id.toLongLong(); });
    return result;
}
} // namespace
QJsonObject Imaging::defaults() {
    return QJsonObject{
        {"input", QJsonObject{{"calibration_nxs", QJsonValue()},
                              {"sample_nxs", QJsonValue()},
                              {"roi_table", QJsonValue()}}},
        {"acquisition", QJsonObject{{"calibration_exposure_s", QJsonValue()},
                                    {"sample_exposure_s", QJsonValue()},
                                    {"monitor", "D_SiC_Current1"},
                                    {"detector", "D_LAMBDA"},
                                    {"foil_thickness_mm", .05},
                                    {"detector_pixel_pitch_mm", .055}}},
        {"scan_loops", QJsonObject{{"calibration_order", ""},
                                   {"calibration_x_points", QJsonValue()},
                                   {"calibration_dcm_points", QJsonValue()},
                                   {"sample_order", ""},
                                   {"sample_y_points", QJsonValue()},
                                   {"sample_z_points", QJsonValue()},
                                   {"sample_dcm_points", QJsonValue()}}},
        {"roi", QJsonObject{{"excluded_ids", QJsonArray()}, {"label_overrides", QJsonObject()}}},
        {"grid", QJsonObject{{"x_step_mm", QJsonValue()}}},
        {"spatial_calibration", QJsonObject{{"match_sample_band", false}, {"verification_nxs", ""},
                                            {"cache_directory", ""}}},
        {"statistics", QJsonObject{{"variance_weighted", false},
                                   {"sigma_floor", QJsonValue()},
                                   {"blank_nxs", ""},
                                   {"blank_exposure_s", QJsonValue()},
                                   {"background_yz_mm", QJsonArray{}}}},
        {"psf", QJsonObject{{"auto_select_model", true},
                            {"min_peak_fraction", .3},
                            {"min_r_squared", .8},
                            {"edge_margin_pixels", 1.2},
                            {"min_observations", 3},
                            {"max_trajectory_rms_pixels", .75},
                            {"sigma_bounds_pixels", QJsonArray{.3, 12.}},
                            {"max_fit_evaluations", 4000}}},
        {"reconstruction", QJsonObject{{"fit_gains", true},
                                       {"anchor_roi", "36_VU-C1"},
                                       {"gain_bounds", QJsonArray{.05, 20.}},
                                       {"max_gain_iterations", 500},
                                       {"gain_tolerance", 1e-5},
                                       {"nnls_max_iterations", 20000},
                                       {"nnls_tolerance", 1e-7}}},
        {"tv", QJsonObject{{"enabled", true},
                           {"lambda", QJsonValue()},
                           {"relative_fraction", .005},
                           {"max_iterations", 40000},
                           {"tolerance", 1e-6}}},
        {"runtime", QJsonObject{{"chunk_frames", 16}, {"output_dir", QJsonValue()}}},
        {"output", QJsonObject{{"save_diagnostics", true}, {"plots", true}, {"save_avizo", false}}}};
}
QJsonObject Imaging::config(const QJsonObject &input, const QString &directory) {
    auto c = merge(defaults(), input);
    auto spatial = c["spatial_calibration"].toObject();
    if (!spatial["match_sample_band"].isBool() || !spatial["verification_nxs"].isString() ||
        !spatial["cache_directory"].isString())
        fail("Invalid spatial calibration configuration");
    for (const auto &key : {QStringLiteral("verification_nxs"), QStringLiteral("cache_directory")}) {
        const auto path = spatial[key].toString().trimmed();
        spatial[key] = path.isEmpty() ? QString{} : QFileInfo(QDir(directory), path).absoluteFilePath();
    }
    c["spatial_calibration"] = spatial;
    auto statistics = c["statistics"].toObject();
    if (!statistics["variance_weighted"].isBool() || !statistics["blank_nxs"].isString() ||
        !statistics["background_yz_mm"].isArray())
        fail("Invalid counting statistics configuration");
    const auto blank = statistics["blank_nxs"].toString().trimmed();
    statistics["blank_nxs"] =
        blank.isEmpty() ? QString() : QFileInfo(QDir(directory), blank).absoluteFilePath();
    for (const auto &key : {"sigma_floor", "blank_exposure_s"})
        if (!statistics[key].isNull() &&
            !(statistics[key].toDouble() > 0 && std::isfinite(statistics[key].toDouble())))
            fail("The variance floor standard deviation and blank exposure must be positive");
    auto rect = statistics["background_yz_mm"].toArray();
    if (!rect.isEmpty() &&
        (!blank.isEmpty() || rect.size() != 4 ||
         !std::all_of(rect.begin(), rect.end(),
                      [](QJsonValue v) { return v.isDouble() && std::isfinite(v.toDouble()); }) ||
         rect[0].toDouble() > rect[1].toDouble() || rect[2].toDouble() > rect[3].toDouble()))
        fail("Choose one background source: blank scan or [Ymin,Ymax,Zmin,Zmax] in mm");
    c["statistics"] = statistics;
    for (const auto &section : {"input", "runtime"}) {
        auto group = c[section].toObject();
        auto keys = QString(section) == "input" ? QStringList{"calibration_nxs", "sample_nxs", "roi_table"}
                                                : QStringList{"output_dir"};
        for (const auto &key : keys) {
            auto value = group[key];
            if (!value.isString() || value.toString().trimmed().isEmpty())
                fail("Missing imaging path: " + QString(section) + "." + key);
            auto path = value.toString();
            if (path.startsWith("~/"))
                path = QDir::home().filePath(path.mid(2));
            group[key] = QFileInfo(QDir(directory), path).absoluteFilePath();
        }
        c[section] = group;
    }
    auto positive = [&](const QString &section, const QString &key, bool integral = false) {
        const auto v = c.value(section).toObject().value(key);
        double n = v.toDouble(NAN);
        if (!v.isDouble() || !std::isfinite(n) || n <= 0 ||
            (integral && (std::floor(n) != n || n > std::numeric_limits<int>::max() - 1)))
            fail("Invalid positive imaging parameter: " + section + "." + key);
    };
    for (const auto &key :
         {"calibration_exposure_s", "sample_exposure_s", "foil_thickness_mm", "detector_pixel_pitch_mm"})
        positive("acquisition", key);
    for (const auto &key : {"monitor", "detector"}) {
        const auto value = c.value("acquisition").toObject().value(key);
        if (!value.isString() || value.toString().isEmpty() || value.toString().contains('/'))
            fail("Invalid detector/monitor name");
    }
    const auto loops = c["scan_loops"].toObject();
    auto checkLoops = [&](const QString &role, const QStringList &motors, const QString &label) {
        QStringList expected;
        for (const auto &motor : motors) {
            auto axis = loopAxis(motor);
            expected.append(axis);
            auto count = loops[loopCountKey(role, axis)];
            if (count.isNull())
                continue;
            double n = count.toDouble(NAN);
            if (!(n >= 2 && n <= 100000) || std::floor(n) != n)
                fail(QCoreApplication::translate("ImagingConfig", "%1的 %2 层帧数须为 2 至 100000 的整数。").arg(label, axis));
        }
        auto order = loops[role + "_order"];
        if (order.isString() && order.toString().isEmpty())
            return;
        auto axes = order.toString().split(',');
        axes.sort();
        expected.sort();
        if (!order.isString() || axes != expected)
            fail(QCoreApplication::translate("ImagingConfig", "%1循环顺序须包含所需的每个轴，且每个轴只出现一次。").arg(label));
    };
    checkLoops("calibration", calibrationMotors, QCoreApplication::translate("ImagingConfig", "薄膜标定"));
    checkLoops("sample", sampleMotors, QCoreApplication::translate("ImagingConfig", "样品成像"));
    if (!c["grid"].toObject()["x_step_mm"].isNull())
        positive("grid", "x_step_mm");
    for (const auto &key :
         {"min_peak_fraction", "min_r_squared", "edge_margin_pixels", "max_trajectory_rms_pixels"})
        positive("psf", key);
    for (const auto &key : {"min_observations", "max_fit_evaluations"})
        positive("psf", key, true);
    auto psf = c["psf"].toObject();
    if (psf["min_observations"].toInt() < 3 || psf["min_peak_fraction"].toDouble() > 1 ||
        psf["min_r_squared"].toDouble() > 1)
        fail("Invalid PSF thresholds");
    for (const auto &pair : {qMakePair(QString("psf"), QString("sigma_bounds_pixels")),
                             qMakePair(QString("reconstruction"), QString("gain_bounds"))}) {
        auto bounds = c[pair.first].toObject()[pair.second].toArray();
        if (bounds.size() != 2 || !bounds[0].isDouble() || !bounds[1].isDouble() ||
            !(bounds[0].toDouble() > 0) || !(bounds[1].toDouble() > bounds[0].toDouble()))
            fail("Invalid imaging parameter bounds");
        if (pair.first == "reconstruction" && (bounds[0].toDouble() > 1 || bounds[1].toDouble() < 1))
            fail("Gain bounds must contain anchor 1");
    }
    auto rec = c["reconstruction"].toObject();
    if (!rec["anchor_roi"].isString() || rec["anchor_roi"].toString().isEmpty())
        fail("Gain anchor key must be explicit");
    for (const auto &key : {"gain_tolerance", "nnls_tolerance"})
        positive("reconstruction", key);
    for (const auto &key : {"max_gain_iterations", "nnls_max_iterations"})
        positive("reconstruction", key, true);
    positive("runtime", "chunk_frames", true);
    positive("tv", "max_iterations", true);
    positive("tv", "tolerance");
    for (const auto &key : {"lambda", "relative_fraction"}) {
        const auto v = c.value("tv").toObject().value(key);
        if (QString(key) == "lambda" && v.isNull())
            continue;
        if (!v.isDouble() || v.toDouble() < 0 || !std::isfinite(v.toDouble()))
            fail("Invalid TV coefficient");
    }
    for (const auto &pair : {qMakePair(QString("reconstruction"), QString("fit_gains")),
                             qMakePair(QString("psf"), QString("auto_select_model")),
                             qMakePair(QString("tv"), QString("enabled")),
                             qMakePair(QString("output"), QString("save_diagnostics")),
                             qMakePair(QString("output"), QString("plots")),
                             qMakePair(QString("output"), QString("save_avizo"))})
        if (!c[pair.first].toObject()[pair.second].isBool())
            fail("Imaging flags must be booleans");
    auto roi = c["roi"].toObject();
    if (!roi["excluded_ids"].isArray() || !roi["label_overrides"].isObject())
        fail("Invalid ROI exclusions/overrides");
    QJsonArray excluded;
    for (auto id : roi["excluded_ids"].toArray()) {
        if (!id.isString() && !id.isDouble())
            fail("Invalid ROI ID");
        excluded.append(id.isString() ? id.toString() : QString::number(id.toDouble(), 'g', 17));
    }
    roi["excluded_ids"] = excluded;
    for (auto v : roi["label_overrides"].toObject())
        if (!v.isString() || v.toString().isEmpty() || v.toString().contains('/'))
            fail("Invalid ROI label override");
    c["roi"] = roi;
    return c;
}
QJsonObject Imaging::loadConfig(const QString &file, QString *bondOutputFile) {
    QFile f(file);
    if (!f.open(QIODevice::ReadOnly))
        fail("Cannot open imaging config: " + file);
    auto bytes = f.readAll();
    if (bytes.startsWith("\xEF\xBB\xBF"))
        bytes.remove(0, 3);
    QJsonParseError error;
    auto doc = QJsonDocument::fromJson(bytes, &error);
    if (error.error != QJsonParseError::NoError || !doc.isObject())
        fail("Invalid imaging JSON: " + error.errorString());
    auto object = doc.object();
    // A resolved result records the actual forward model used, including pre-selection releases.
    if (QFileInfo(file).fileName() == "config.resolved.json") {
        auto psf = object["psf"].toObject();
        if (!psf.contains("auto_select_model")) {
            psf["auto_select_model"] = false;
            object["psf"] = psf;
        }
        if (!object.contains("spatial_calibration"))
            object["spatial_calibration"] = QJsonObject{{"match_sample_band", true}};
    }
    const auto output = object.take("bond_output_file").toString().trimmed();
    object.remove("bond_contrast_output_file"); // Accept configurations saved by older XRI versions.
    if (bondOutputFile)
        *bondOutputFile = output.isEmpty() ? QString{} : QFileInfo(QFileInfo(file).absoluteDir(), output).absoluteFilePath();
    return config(object, QFileInfo(file).absolutePath());
}
ImagingValidation Imaging::validate(const QJsonObject &c, const std::atomic_bool &stop,
                                    bool selectIntegrationBand) {
    auto validated = config(c, QDir::currentPath());
    auto inputs = validated["input"].toObject(), acq = validated["acquisition"].toObject();
    for (auto i = inputs.begin(); i != inputs.end(); ++i)
        if (!QFileInfo(i->toString()).isFile())
            fail("Input does not exist: " + i->toString());
    ImagingValidation out;
    auto loops = validated["scan_loops"].toObject();
    const auto foilFile = inputs["calibration_nxs"].toString();
    cancel(stop);
    {
        const ImagingH5 source(foilFile);
        if (source.contains("xri_prepared_foil")) {
            auto &foil = out.calibration;
            foil.preparedFoil = QJsonDocument::fromJson(source.text("xri_prepared_foil").toUtf8()).object();
            const auto &preset = foil.preparedFoil;
            const auto image = preset["image_shape"].toArray();
            if (preset["schema"] != "XRI.prepared-foil.v1" || image.size() != 2 ||
                image[0].toInt() < 1 || image[1].toInt() < 1 ||
                preset["source_scan"].toString().isEmpty())
                fail(QCoreApplication::translate("ImagingConfig", "预计算薄膜标定文件无效。"));
            bool compatible = preset["psf"] == validated["psf"] &&
                              !validated["spatial_calibration"].toObject()["match_sample_band"].toBool();
            for (const auto &field : {"detector", "monitor", "calibration_exposure_s",
                                      "foil_thickness_mm", "detector_pixel_pitch_mm"})
                compatible &= preset[field] == acq[field];
            foil.file = foilFile;
            foil.imageSize = {image[1].toInt(), image[0].toInt()};
            foil.axes = {source.read("x_mm"), source.read("energy_kev")};
            for (const auto &axis : foil.axes) {
                if (axis.size() < 2)
                    fail(QCoreApplication::translate("ImagingConfig", "预计算薄膜标定文件无效。"));
                energyWeights(axis);
            }
            const auto order = preset["loop_order"].toArray();
            for (const auto &name : order)
                foil.loopOrder.append(name.toString());
            auto sorted = foil.loopOrder;
            sorted.sort();
            const auto cells = foil.axes[0].size() * foil.axes[1].size();
            if (sorted != QStringList{"DCM", "X"} || preset["actual_frames"].toInteger() < cells)
                fail(QCoreApplication::translate("ImagingConfig", "预计算薄膜标定文件无效。"));
            const auto requestedOrder = loops["calibration_order"].toString();
            compatible &= requestedOrder.isEmpty() || requestedOrder == foil.loopOrder.join(',');
            for (int axis = 0; axis < calibrationMotors.size(); ++axis) {
                const auto count = loops[loopCountKey("calibration", loopAxis(calibrationMotors[axis]))];
                compatible &= count.isNull() || count.toInteger() == foil.axes[axis].size();
            }
            if (!compatible)
                fail(QCoreApplication::translate(
                    "ImagingConfig", "当前参数与预计算标定不一致。请恢复预设参数，或导入薄膜 NeXus 重新标定。"));
        }
    }
    if (out.calibration.preparedFoil.isEmpty())
        out.calibration = grid(foilFile, "calibration", calibrationMotors,
                               acq["calibration_exposure_s"].toDouble(), acq, loops, stop);
    out.sample = grid(inputs["sample_nxs"].toString(), "sample", sampleMotors,
                      acq["sample_exposure_s"].toDouble(), acq, loops, stop);
    if (out.calibration.imageSize != out.sample.imageSize)
        fail("Calibration/sample detector dimensions differ");
    auto report = [](const ImagingGrid &g, const QStringList &motors, double exposure) {
        QJsonObject axes;
        QMap<QString, qsizetype> counts;
        QJsonArray shape, loopCounts;
        for (int i = 0; i < motors.size(); ++i) {
            axes[motors[i]] = json(g.axes[i]);
            shape.append(g.axes[i].size());
            counts[loopAxis(motors[i])] = g.axes[i].size();
        }
        for (const auto &axis : g.loopOrder)
            loopCounts.append(counts.value(axis));
        const auto cells = g.preparedFoil.isEmpty() ? qint64(g.occupancy.size())
                                                    : qint64(g.axes[0].size()) * g.axes[1].size();
        const auto frames = g.preparedFoil.isEmpty() ? qint64(g.monitor.size())
                                                     : g.preparedFoil["actual_frames"].toInteger();
        QJsonObject result{{"source_file", g.file},
                           {"axes", axes},
                           {"grid_shape", shape},
                           {"loop_order_outer_to_inner", json(g.loopOrder)},
                           {"loop_counts_outer_to_inner", loopCounts},
                           {"actual_frames", frames},
                           {"expected_grid_cells", cells},
                           {"repeated_frames", frames - cells},
                           {"image_shape", QJsonArray{g.imageSize.height(), g.imageSize.width()}},
                           {"exposure_s", exposure}};
        if (!g.preparedFoil.isEmpty())
            result["prepared_foil"] = g.preparedFoil;
        return result;
    };
    out.report = {
        {"calibration", report(out.calibration, calibrationMotors, acq["calibration_exposure_s"].toDouble())},
        {"sample", report(out.sample, sampleMotors, acq["sample_exposure_s"].toDouble())}};
    if (selectIntegrationBand) {
        const auto stats = validated["statistics"].toObject();
        const auto rect = stats["background_yz_mm"].toArray();
        if (!rect.isEmpty()) {
            auto count = [&](int d) {
                return std::count_if(out.sample.axes[d].begin(), out.sample.axes[d].end(), [&](double t) {
                    return t >= rect[d * 2].toDouble() && t <= rect[d * 2 + 1].toDouble();
                });
            };
            if (count(0) * count(1) < 2)
                fail("Background YZ rectangle must include at least two scan positions");
        }
        const auto path = stats["blank_nxs"].toString();
        const double exposure = stats["blank_exposure_s"].toDouble(acq["sample_exposure_s"].toDouble());
        const auto spatial = validated["spatial_calibration"].toObject();
        const auto verification = spatial["verification_nxs"].toString();
        if (!verification.isEmpty()) {
            const QFileInfo file(verification);
            if (!file.isFile() || file.canonicalFilePath() == QFileInfo(out.calibration.file).canonicalFilePath())
                fail(QCoreApplication::translate("ImagingConfig", "空间响应验证需要另一份独立的薄膜扫描。"));
            out.verification = grid(verification, "calibration", calibrationMotors, 1., acq,
                                    defaults()["scan_loops"].toObject(), stop);
            if (out.verification.imageSize != out.calibration.imageSize)
                fail("Verification/calibration detector dimensions differ");
            auto verificationReport = report(out.verification, calibrationMotors, 1.);
            verificationReport.remove("exposure_s");
            out.report["spatial_verification_scan"] = verificationReport;
        }
        selectEnergyBand(out, path.isEmpty() ? BlankScan{} : blankScan(path, out.sample, acq, exposure, stop),
                         spatial["match_sample_band"].toBool());
    }
    validateRois(out, validated);
    return out;
}
void Imaging::validateRois(ImagingValidation &data, const QJsonObject &settings) {
    auto rois = roiTable(settings["input"].toObject()["roi_table"].toString(),
                         settings["acquisition"].toObject()["detector"].toString(),
                         data.calibration.imageSize, settings["roi"].toObject());
    QStringList keys, excluded;
    qsizetype pixels = 0;
    std::unique_ptr<ImagingH5> prepared;
    if (!data.calibration.preparedFoil.isEmpty())
        prepared = std::make_unique<ImagingH5>(data.calibration.file);
    for (const auto &r : rois) {
        if (prepared && !r.excluded) {
            const auto path = "rois/" + r.key + "/pixels_rc";
            QVector<quint64> shape;
            const auto positions = prepared->contains(path) ? prepared->read(path, &shape) : QVector<double>{};
            bool same = shape == QVector<quint64>{quint64(r.pixels.size()), 2};
            for (qsizetype i = 0; same && i < r.pixels.size(); ++i)
                same = positions[2 * i] == r.pixels[i].y() && positions[2 * i + 1] == r.pixels[i].x();
            if (!same)
                fail(QCoreApplication::translate(
                         "ImagingConfig", "晶体 %1 的 ROI 与预计算标定不一致。请使用随附的 ROI，或导入薄膜 NeXus 重新标定。")
                         .arg(r.key));
        }
        if (r.excluded)
            excluded.append(r.id);
        else {
            keys.append(r.key);
            pixels += r.pixels.size();
        }
    }
    const auto rec = settings["reconstruction"].toObject();
    if (rec["fit_gains"].toBool() && !keys.contains(rec["anchor_roi"].toString()))
        fail("Gain anchor is not included: " + rec["anchor_roi"].toString());
    data.rois = std::move(rois);
    data.report["roi_keys"] = json(keys);
    data.report["excluded_roi_ids"] = json(excluded);
    data.report["roi_pixels"] = pixels;
}
QVector<double> Imaging::energyWeights(const QVector<double> &axis) {
    if (axis.isEmpty())
        fail("Empty energy axis");
    QVector<double> w(axis.size(), 0);
    for (int i = 0; i < axis.size(); ++i)
        if (!std::isfinite(axis[i]) || (i && axis[i] <= axis[i - 1]))
            fail("Invalid increasing energy axis");
    if (axis.size() == 1) {
        w[0] = 1;
        return w;
    }
    for (int i = 1; i < axis.size(); ++i) {
        double step = (axis[i] - axis[i - 1]) / 2;
        w[i - 1] += step;
        w[i] += step;
    }
    return w;
}
QPair<QVector<double>, QVector<double>> Imaging::matchedWeights(const QVector<double> &source,
                                                                const QVector<double> &target,
                                                                const QString &referenceLabel) {
    if (source.size() < 2)
        fail("Calibration energy grid requires at least two points");
    energyWeights(source);
    auto w = energyWeights(target);
    if (target.first() < source.first() - energyCoverageTolerance ||
        target.last() > source.last() + energyCoverageTolerance) {
        const auto label = referenceLabel.isEmpty()
                               ? QCoreApplication::translate("ImagingConfig", "薄膜标定")
                               : referenceLabel;
        fail(QCoreApplication::translate(
                 "ImagingConfig", "样品能区 %1–%2 keV 超出%3能区 %4–%5 keV。请选择覆盖整个样品能区的%3文件。")
                 .arg(QString::number(target.first(), 'g', 12), QString::number(target.last(), 'g', 12),
                      label, QString::number(source.first(), 'g', 12), QString::number(source.last(), 'g', 12)));
    }
    return {interpolateWeights(source, target, w), w};
}
