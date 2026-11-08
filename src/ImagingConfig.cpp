#include "CrystalRoi.h"
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
using namespace ImageUtil;
namespace {
const QStringList calibrationMotors{"M_Stage_X", "M_DCM_B5Energy"};
const QStringList sampleMotors{"M_Stage_Y", "M_Stage_Z", "M_DCM_B5Energy"};

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
        fail(QStringLiteral("扫描计划中的电机参数不完整：%1").arg(file));
    QStringList actual;
    for (int i = 1; i < args.size(); i += 4) {
        auto name = args[i].toString();
        name.replace('.', '_');
        actual.append(name);
    }
    auto roleLabel = role == "calibration" ? QStringLiteral("薄膜标定") : QStringLiteral("样品成像");
    auto sortedActual = actual;
    auto sortedMotors = motors;
    sortedActual.sort();
    sortedMotors.sort();
    if (sortedActual != sortedMotors) {
        auto hint = actual.contains("M_DCM_B5Energy")
                        ? QStringLiteral("请选择扫描轴与用途匹配的文件。")
                        : QStringLiteral("该扫描计划没有能量轴，不能用于当前需要能量积分的处理流程。");
        fail(QStringLiteral("%1文件的扫描轴不匹配：\n%2\n实际：%3\n所需轴：%4\n%5")
                 .arg(roleLabel, file, actual.join(", "), motors.join(", "), hint));
    }
    for (const auto &motor : actual)
        result.loopOrder.append(loopAxis(motor));
    const auto requestedOrder = loops[role + "_order"].toString();
    if (!requestedOrder.isEmpty() && requestedOrder != result.loopOrder.join(',')) {
        fail(QStringLiteral("%1文件的循环顺序与设置不符。\n%2\n文件（外到内）：%3\n设置（外到内）：%4")
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
            fail(QStringLiteral("%1的 %2 轴需要递增的扫描范围及至少两个点。").arg(roleLabel, axisName));
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
                fail(QStringLiteral("%1的 %2 轴逐帧坐标与 %3 个扫描点不匹配。请核对扫描范围和该层帧数。")
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
        fail(QStringLiteral("%1 扫描不完整：%2 个计划网格单元缺测；不以零填充。")
                 .arg(motors.size() == 2 ? "薄膜" : "样品")
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
                           {}});
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
        for (int y = r.rect.top(); !r.rect.isEmpty() && y <= r.rect.bottom(); ++y)
            for (int x = r.rect.left(); x <= r.rect.right(); ++x) {
                r.pixels.append({x, y});
                if (!r.excluded && ++coverage[y * size.width() + x] > 1)
                    fail("Valid imaging ROIs overlap; shared pixels must not be treated independently");
            }
        if (!r.excluded)
            ++active;
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
        {"psf", QJsonObject{{"min_peak_fraction", .3},
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
        {"output", QJsonObject{{"save_diagnostics", true}, {"plots", true}}}};
}
QJsonObject Imaging::config(const QJsonObject &input, const QString &directory) {
    auto c = merge(defaults(), input);
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
        auto v = c[section].toObject()[key];
        double n = v.toDouble(NAN);
        if (!v.isDouble() || !std::isfinite(n) || n <= 0 ||
            (integral && (std::floor(n) != n || n > std::numeric_limits<int>::max() - 1)))
            fail("Invalid positive imaging parameter: " + section + "." + key);
    };
    for (const auto &key :
         {"calibration_exposure_s", "sample_exposure_s", "foil_thickness_mm", "detector_pixel_pitch_mm"})
        positive("acquisition", key);
    for (const auto &key : {"monitor", "detector"}) {
        auto value = c["acquisition"].toObject()[key];
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
                fail(QStringLiteral("%1的 %2 层帧数须为 2 至 100000 的整数。").arg(label, axis));
        }
        auto order = loops[role + "_order"];
        if (order.isString() && order.toString().isEmpty())
            return;
        auto axes = order.toString().split(',');
        axes.sort();
        expected.sort();
        if (!order.isString() || axes != expected)
            fail(QStringLiteral("%1循环顺序须包含所需的每个轴，且每个轴只出现一次。").arg(label));
    };
    checkLoops("calibration", calibrationMotors, QStringLiteral("薄膜标定"));
    checkLoops("sample", sampleMotors, QStringLiteral("样品成像"));
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
        auto v = c["tv"].toObject()[key];
        if (QString(key) == "lambda" && v.isNull())
            continue;
        if (!v.isDouble() || v.toDouble() < 0 || !std::isfinite(v.toDouble()))
            fail("Invalid TV coefficient");
    }
    for (const auto &pair : {qMakePair(QString("reconstruction"), QString("fit_gains")),
                             qMakePair(QString("tv"), QString("enabled")),
                             qMakePair(QString("output"), QString("save_diagnostics")),
                             qMakePair(QString("output"), QString("plots"))})
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
    const auto output = object.take("bond_output_file");
    if (bondOutputFile) {
        const auto path = output.toString().trimmed();
        *bondOutputFile =
            path.isEmpty() ? QString{} : QFileInfo(QFileInfo(file).absoluteDir(), path).absoluteFilePath();
    }
    return config(object, QFileInfo(file).absolutePath());
}
ImagingValidation Imaging::validate(const QJsonObject &c, const std::atomic_bool &stop,
                                    bool requireMatchedBand) {
    auto validated = config(c, QDir::currentPath());
    auto inputs = validated["input"].toObject(), acq = validated["acquisition"].toObject();
    for (auto i = inputs.begin(); i != inputs.end(); ++i)
        if (!QFileInfo(i->toString()).isFile())
            fail("Input does not exist: " + i->toString());
    ImagingValidation out;
    auto loops = validated["scan_loops"].toObject();
    out.calibration = grid(inputs["calibration_nxs"].toString(), "calibration", calibrationMotors,
                           acq["calibration_exposure_s"].toDouble(), acq, loops, stop);
    out.sample = grid(inputs["sample_nxs"].toString(), "sample", sampleMotors,
                      acq["sample_exposure_s"].toDouble(), acq, loops, stop);
    if (out.calibration.imageSize != out.sample.imageSize)
        fail("Calibration/sample detector dimensions differ");
    if (requireMatchedBand)
        matchedWeights(out.calibration.axes[1], out.sample.axes[2]);
    out.rois = roiTable(inputs["roi_table"].toString(), acq["detector"].toString(), out.calibration.imageSize,
                        validated["roi"].toObject());
    QStringList keys, excluded;
    int pixels = 0;
    for (const auto &r : out.rois) {
        if (r.excluded)
            excluded.append(r.id);
        else {
            keys.append(r.key);
            pixels += int(r.pixels.size());
        }
    }
    auto rec = validated["reconstruction"].toObject();
    if (rec["fit_gains"].toBool() && !keys.contains(rec["anchor_roi"].toString()))
        fail("Gain anchor is not included: " + rec["anchor_roi"].toString());
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
        return QJsonObject{{"source_file", g.file},
                           {"axes", axes},
                           {"grid_shape", shape},
                           {"loop_order_outer_to_inner", json(g.loopOrder)},
                           {"loop_counts_outer_to_inner", loopCounts},
                           {"actual_frames", g.monitor.size()},
                           {"expected_grid_cells", g.occupancy.size()},
                           {"repeated_frames", g.monitor.size() - g.occupancy.size()},
                           {"image_shape", QJsonArray{g.imageSize.height(), g.imageSize.width()}},
                           {"exposure_s", exposure}};
    };
    out.report = {
        {"calibration", report(out.calibration, calibrationMotors, acq["calibration_exposure_s"].toDouble())},
        {"sample", report(out.sample, sampleMotors, acq["sample_exposure_s"].toDouble())},
        {"roi_keys", json(keys)},
        {"excluded_roi_ids", json(excluded)},
        {"roi_pixels", pixels}};
    return out;
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
                                                                const QVector<double> &target) {
    if (source.size() < 2)
        fail("Calibration energy grid requires at least two points");
    energyWeights(source);
    auto w = energyWeights(target);
    if (target.first() < source.first() - 1e-10 || target.last() > source.last() + 1e-10)
        fail("Sample energy range extends outside calibration; extrapolation is not allowed");
    QVector<double> v(source.size(), 0);
    for (int k = 0; k < target.size(); ++k) {
        double e = std::clamp(target[k], source.first(), source.last());
        auto it = std::upper_bound(source.begin(), source.end(), e);
        int j = std::clamp(int(it - source.begin()) - 1, 0, int(source.size()) - 2);
        double t = (e - source[j]) / (source[j + 1] - source[j]);
        v[j] += w[k] * (1 - t);
        v[j + 1] += w[k] * t;
    }
    return {v, w};
}
