#include "CrystalRoi.h"
#include <QCoreApplication>
#include <QFile>
#include <QFileInfo>
#include <QHash>
#include <QJsonArray>
#include <QJsonDocument>
#include <QJsonObject>
#include <QSet>
#include <QTextStream>
#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>

namespace {
struct TextRoi {
    QString id, label, location;
    std::array<int, 4> bounds{};
    bool enabled = true;
    bool empty = false;
};

QStringList textColumns(const QString &line) {
    auto fields = line.contains('\t') ? line.split('\t') : line.simplified().split(' ');
    for (auto &field : fields)
        field = field.trimmed();
    return fields;
}

QVector<TextRoi> readTextRois(const QString &path) {
    const auto fail = [&](int line, const QString &reason) {
        throw std::runtime_error(
            QString("%1:%2 · %3").arg(QFileInfo(path).fileName()).arg(line).arg(reason).toStdString());
    };
    QFile file(path);
    if (!file.open(QIODevice::ReadOnly | QIODevice::Text))
        fail(1, file.errorString());
    QTextStream input(&file);
    const auto header = textColumns(input.readLine());
    for (const auto *column : {"roi_label", "x1", "x2", "y1", "y2"}) {
        if (header.count(column) != 1)
            fail(1, QCoreApplication::translate("CrystalRoi", "缺少或重复的文本 ROI 列 %1。").arg(column));
    }

    QVector<TextRoi> rows;
    int line = 1;
    while (!input.atEnd()) {
        const auto text = input.readLine();
        ++line;
        if (text.trimmed().isEmpty())
            continue;
        auto fields = textColumns(text);
        QString id;
        if (fields.size() == header.size() + 1 && header.first() == "roi_label")
            id = fields.takeFirst();
        else if (header.first().isEmpty() || header.first() == "id" || header.first() == "index")
            id = fields.value(0);
        if (fields.size() != header.size())
            fail(line, QCoreApplication::translate("CrystalRoi", "文本 ROI 列数与表头不一致。"));
        const auto field = [&](const QString &name) { return fields.value(header.indexOf(name)); };
        const auto integer = [&](const QString &name, bool optional = false) {
            const auto value = field(name);
            if (optional && value.isEmpty())
                return 0;
            bool ok;
            const int n = value.toInt(&ok);
            if (!ok || n < -100000 || n > 100000)
                fail(line, QCoreApplication::translate("CrystalRoi", "文本 ROI 列 %1 须为有效整数。").arg(name));
            return n;
        };
        TextRoi roi;
        roi.id = id;
        roi.label = field("roi_label").toUpper();
        roi.location = QString("%1:%2").arg(QFileInfo(path).fileName()).arg(line);
        const int x1 = integer("x1"), x2 = integer("x2"), y1 = integer("y1"), y2 = integer("y2");
        roi.empty = x1 == 0 && x2 == 0 && y1 == 0 && y2 == 0;
        if (!roi.empty && (x1 < 0 || y1 < 0 || x2 <= x1 || y2 <= y1))
            fail(line, QCoreApplication::translate("CrystalRoi", "文本 ROI 原始矩形无效。"));
        const int dx = integer("x_shift", true), dy = integer("y_shift", true);
        const int ex = integer("x_expand", true), ey = integer("y_expand", true);
        roi.bounds = {x1 + dx - ex, x2 + dx + ex, y1 + dy - ey, y2 + dy + ey};
        const auto add = field("add").toLower();
        if (!add.isEmpty() && add != "0" && add != "1" && add != "false" && add != "true")
            fail(line, QCoreApplication::translate("CrystalRoi", "文本 ROI 的 add 值无效。"));
        roi.enabled = add != "0" && add != "false";
        rows.append(std::move(roi));
    }
    if (rows.isEmpty())
        fail(1, QCoreApplication::translate("CrystalRoi", "文本 ROI 没有数据行。"));
    return rows;
}

bool hasLegacyVuLayout(const QVector<TextRoi> &rows) {
    const auto has = [&](const char *id, const char *label, bool empty) {
        return std::any_of(rows.cbegin(), rows.cend(), [&](const TextRoi &r) {
            return r.id == QString::fromLatin1(id) && r.label == QString::fromLatin1(label) && r.empty == empty;
        });
    };
    // This legacy layout is documented in docs/imaging-migration.md.
    return has("32", "VU-B3", false) && has("35", "VU-C3", false) && has("38", "VU-C3", false) &&
           has("41", "VU-E3", false) && (has("31", "VU-A3", true) || has("44", "VU-A3", true));
}
} // namespace

QVector<QPoint> CrystalRoi::pixels() const {
    if (!mask.isEmpty())
        return mask;
    QVector<QPoint> result;
    for (int y = rect.top(); !rect.isEmpty() && y <= rect.bottom(); ++y)
        for (int x = rect.left(); x <= rect.right(); ++x)
            result.append({x, y});
    return result;
}

QVector<QPoint> CrystalRois::decodeMask(const QJsonArray &runs, const QSize &size) {
    if (size.isEmpty())
        return {};
    const qint64 limit = qint64(size.width()) * size.height();
    qint64 previousEnd = 0;
    QVector<QPoint> mask;
    for (const auto &value : runs) {
        const auto pair = value.toArray();
        if (pair.size() != 2)
            return {};
        const qint64 start = pair[0].toInteger(-1), length = pair[1].toInteger(-1);
        if (start < previousEnd || start >= limit || length <= 0 || length > limit - start)
            return {};
        const qint64 end = start + length;
        // Limit a single ROI even if a malformed file encodes a huge rectangle in one run.
        if (mask.size() + length > 1024 * 1024)
            return {};
        for (qint64 i = start; i < end; ++i)
            mask.append({int(i % size.width()), int(i / size.width())});
        previousEnd = end;
    }
    return mask;
}

CrystalRois::CrystalRois() {
    // Notebook order: 60 on LAMBDA (VB, HL, VU, HB), 15 on MINIPIX (VD).
    // Canonical labels are unique; copied notebook output contains duplicated VU labels.
    for (const auto &module : {"VB", "HL", "VU", "HB", "VD"})
        for (QChar letter : QStringLiteral("ABCDE"))
            for (int i = 1; i <= 3; ++i)
                m_items.append({QString("%1-%2%3").arg(module).arg(letter).arg(i),
                                QString(module) == "VD" ? "D_MINIPIX" : "D_LAMBDA",
                                {},
                                RoiState::Pending,
                                {}});
}
bool CrystalRois::assign(int row, const QRect &rect, const QSize &size) {
    if (row < 0 || row >= m_items.size() || rect.width() <= 0 || rect.height() <= 0 ||
        !QRect(QPoint(0, 0), size).contains(rect))
        return false;
    m_items[row].rect = rect;
    m_items[row].state = RoiState::Valid;
    m_items[row].mask.clear();
    return true;
}
bool CrystalRois::assignMask(int row, const QVector<QPoint> &mask, const QSize &size) {
    if (mask.isEmpty())
        return false;
    const QRect detector(QPoint(), size);
    QRect rect;
    qint64 previous = -1;
    for (const auto &p : mask) {
        const qint64 offset = qint64(p.y()) * size.width() + p.x();
        if (!detector.contains(p) || offset <= previous)
            return false;
        rect = rect.united(QRect(p, QSize(1, 1)));
        previous = offset;
    }
    if (!assign(row, rect, size))
        return false;
    m_items[row].mask = mask;
    return true;
}
void CrystalRois::invalidate(int row) {
    if (row >= 0 && row < m_items.size())
        m_items[row].state = RoiState::Invalid;
}
void CrystalRois::clear(int row) {
    if (row >= 0 && row < m_items.size()) {
        m_items[row].state = RoiState::Pending;
        m_items[row].rect = {};
        m_items[row].mask.clear();
    }
}
int CrystalRois::count(RoiState s) const {
    int n = 0;
    for (const auto &r : m_items)
        if (r.state == s)
            ++n;
    return n;
}

QByteArray CrystalRois::convertTextFiles(const QStringList &files, const QVector<QSize> &sizes, QStringList *notes) {
    const auto tr = [](const char *text) { return QCoreApplication::translate("CrystalRoi", text); };
    if (files.isEmpty() || sizes.size() != 2)
        throw std::runtime_error(tr("请选择文本 ROI，并提供两台探测器的尺寸。").toStdString());
    CrystalRois converted;
    QHash<QString, int> indices;
    for (int i = 0; i < converted.items().size(); ++i) {
        indices.insert(converted.items()[i].label, i);
        converted.invalidate(i);
    }

    QStringList messages;
    QSet<int> mapped;
    for (const auto &file : files) {
        const auto rows = readTextRois(file);
        const bool legacyVuLayout = hasLegacyVuLayout(rows);
        for (const auto &roi : rows) {
            auto label = roi.label;
            if (legacyVuLayout && roi.id == "35" && roi.label == "VU-C3") {
                label = "VU-A3";
                messages.append(tr("%1 · 按旧成像映射将编号 35 的 VU-C3 修正为 VU-A3。").arg(roi.location));
            }
            const int index = indices.value(label, -1);
            if (index < 0)
                throw std::runtime_error(tr("%1 · 无法识别晶体 %2。").arg(roi.location, label).toStdString());
            if (roi.empty)
                continue;
            if (mapped.contains(index))
                throw std::runtime_error(tr("%1 · 晶体 %2 有多个非空 ROI，请先核对编号。")
                                             .arg(roi.location, label).toStdString());
            mapped.insert(index);
            const QSize size = sizes[index < 60 ? 0 : 1];
            const auto [x1, sourceX2, y1, sourceY2] = roi.bounds;
            const int x2 = std::min(sourceX2, size.width()), y2 = std::min(sourceY2, size.height());
            if (x1 < 0 || y1 < 0 || x2 <= x1 || y2 <= y1 ||
                !converted.assign(index, QRect(x1, y1, x2 - x1, y2 - y1), size))
                throw std::runtime_error(tr("%1 · 晶体 %2 的 ROI 无可用探测器范围。")
                                             .arg(roi.location, label).toStdString());
            if (x2 != sourceX2 || y2 != sourceY2)
                messages.append(tr("%1 · %2 的矩形上限从 (%3, %4) 截为 (%5, %6)。")
                                    .arg(roi.location, label).arg(sourceX2).arg(sourceY2).arg(x2).arg(y2));
            if (!roi.enabled)
                converted.invalidate(index);
        }
    }
    if (converted.count(RoiState::Valid) == 0)
        throw std::runtime_error(tr("文本 ROI 中没有可用晶体。").toStdString());
    auto obj = QJsonDocument::fromJson(converted.serialize({}, sizes)).object();
    obj["legacy_source_files"] = QJsonArray::fromStringList(files);
    obj["legacy_conversion_notes"] = QJsonArray::fromStringList(messages);
    if (notes)
        *notes = messages;
    return QJsonDocument(obj).toJson();
}
QByteArray CrystalRois::serialize(const QString &source, const QVector<QSize> &sizes) const {
    QJsonArray geometry, items;
    for (const auto &s : sizes)
        geometry.append(QJsonArray{s.width(), s.height()});
    int version = 1;
    for (const auto &r : m_items) {
        QJsonObject row{{"roi_label", r.label},
                                 {"detector", r.detector},
                                 {"state", r.state == RoiState::Valid     ? "valid"
                                           : r.state == RoiState::Invalid ? "invalid"
                                                                          : "pending"},
                                 {"add", r.state == RoiState::Valid},
                                 {"x1", r.rect.x()},
                                 {"x2", r.rect.x() + r.rect.width()},
                                 {"y1", r.rect.y()},
                                 {"y2", r.rect.y() + r.rect.height()}};
        if (!r.mask.isEmpty()) {
            const int width = sizes.value(r.detector == "D_MINIPIX" ? 1 : 0).width();
            QJsonArray runs;
            qint64 start = -1, end = -1;
            for (const auto &p : r.mask) {
                const qint64 offset = qint64(p.y()) * width + p.x();
                if (offset != end) {
                    if (start >= 0)
                        runs.append(QJsonArray{start, end - start});
                    start = offset;
                }
                end = offset + 1;
            }
            runs.append(QJsonArray{start, end - start});
            row["mask_runs"] = runs;
            version = 2;
        }
        items.append(row);
    }
    return QJsonDocument(QJsonObject{{"schema", "xri-crystal-roi"},
                                     {"version", version},
                                     {"source", source},
                                     {"coordinate_convention", "[x1,x2), [y1,y2)"},
                                     {"detector_sizes", geometry},
                                     {"crystals", items}})
        .toJson();
}
bool CrystalRois::restore(const QByteArray &bytes, const QVector<QSize> &sizes, QString *error) {
    auto fail = [error](const QString &s) {
        if (error)
            *error = s;
        return false;
    };
    QJsonParseError parse;
    auto doc = QJsonDocument::fromJson(bytes, &parse);
    auto obj = doc.object();
    const auto version = obj.value("version").toInt();
    if (parse.error != QJsonParseError::NoError || obj.value("schema") != "xri-crystal-roi" ||
        (version != 1 && version != 2))
        return fail(QCoreApplication::translate("CrystalRoi", "不是支持的 XRI ROI 方案。"));
    QJsonArray geometry;
    for (const auto &s : sizes)
        geometry.append(QJsonArray{s.width(), s.height()});
    if (obj.value("detector_sizes").toArray() != geometry)
        return fail(QCoreApplication::translate("CrystalRoi", "探测器尺寸与当前文件不一致。"));
    auto rows = obj.value("crystals").toArray();
    if (rows.size() != 75)
        return fail(QCoreApplication::translate("CrystalRoi", "方案必须包含 75 个晶体。"));
    CrystalRois next;
    for (int i = 0; i < 75; ++i) {
        auto r = rows[i].toObject();
        auto &target = next.m_items[i];
        if (r.value("roi_label") != target.label || r.value("detector") != target.detector)
            return fail(QCoreApplication::translate("CrystalRoi", "晶体编号或探测器对应关系错误。"));
        QString state = r.value("state").toString();
        if (state != "pending" && state != "valid" && state != "invalid")
            return fail(QCoreApplication::translate("CrystalRoi", "ROI 状态无效。"));
        QVector<int> v;
        for (const auto &key : {"x1", "x2", "y1", "y2"}) {
            auto val = r.value(key);
            double d = val.toDouble(-1);
            if (!val.isDouble() || d < 0 || d > 100000 || std::floor(d) != d)
                return fail(QCoreApplication::translate("CrystalRoi", "ROI 坐标必须为非负整数。"));
            v.append(int(d));
        }
        QRect rect(v[0], v[2], v[1] - v[0], v[3] - v[2]);
        if (rect.isNull() && (v[0] || v[1] || v[2] || v[3]))
            return fail(QCoreApplication::translate("CrystalRoi", "空 ROI 必须使用全零坐标。"));
        const auto size = sizes.value(i < 60 ? 0 : 1);
        if (!rect.isNull() && !next.assign(i, rect, size))
            return fail(QCoreApplication::translate("CrystalRoi", "ROI 坐标超出探测器范围或为空。"));
        if (state == "valid" && rect.isNull())
            return fail(QCoreApplication::translate("CrystalRoi", "有效晶体必须具有非空 ROI。"));
        if (r.contains("mask_runs") &&
            (version != 2 || !next.assignMask(i, decodeMask(r["mask_runs"].toArray(), size), size) ||
             target.rect != rect))
            return fail(QCoreApplication::translate("CrystalRoi", "晶体 %1 的掩膜编码或边界无效。").arg(target.label));
        target.rect = rect;
        target.state = state == "valid"     ? RoiState::Valid
                       : state == "invalid" ? RoiState::Invalid
                                            : RoiState::Pending;
    }
    *this = next;
    return true;
}
