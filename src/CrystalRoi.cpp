#include "CrystalRoi.h"
#include <QCoreApplication>
#include <QJsonArray>
#include <QJsonDocument>
#include <QJsonObject>
#include <cmath>

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
