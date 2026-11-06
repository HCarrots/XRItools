#include "CrystalRoi.h"
#include <QJsonArray>
#include <QJsonDocument>
#include <QJsonObject>
#include <cmath>

CrystalRois::CrystalRois() {
    // Notebook order: 60 on LAMBDA (VB, HL, VU, HB), 15 on MINIPIX (VD).
    // Canonical labels are unique; copied notebook output contains duplicated VU labels.
    for (const auto &module : {"VB", "HL", "VU", "HB", "VD"})
        for (QChar letter : QStringLiteral("ABCDE"))
            for (int i = 1; i <= 3; ++i)
                m_items.append({QString("%1-%2%3").arg(module).arg(letter).arg(i),
                                QString(module) == "VD" ? "D_MINIPIX" : "D_LAMBDA",
                                {},
                                RoiState::Pending});
}
bool CrystalRois::assign(int row, const QRect &rect, const QSize &size) {
    if (row < 0 || row >= m_items.size() || rect.width() <= 0 || rect.height() <= 0 ||
        !QRect(QPoint(0, 0), size).contains(rect))
        return false;
    m_items[row].rect = rect;
    m_items[row].state = RoiState::Valid;
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
    for (const auto &r : m_items)
        items.append(QJsonObject{{"roi_label", r.label},
                                 {"detector", r.detector},
                                 {"state", r.state == RoiState::Valid     ? "valid"
                                           : r.state == RoiState::Invalid ? "invalid"
                                                                          : "pending"},
                                 {"add", r.state == RoiState::Valid},
                                 {"x1", r.rect.x()},
                                 {"x2", r.rect.x() + r.rect.width()},
                                 {"y1", r.rect.y()},
                                 {"y2", r.rect.y() + r.rect.height()}});
    return QJsonDocument(QJsonObject{{"schema", "xri-crystal-roi"},
                                     {"version", 1},
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
    if (parse.error != QJsonParseError::NoError || obj.value("schema") != "xri-crystal-roi" ||
        obj.value("version").toInt() != 1)
        return fail(QStringLiteral("不是支持的 XRI ROI 方案。"));
    QJsonArray geometry;
    for (const auto &s : sizes)
        geometry.append(QJsonArray{s.width(), s.height()});
    if (obj.value("detector_sizes").toArray() != geometry)
        return fail(QStringLiteral("探测器尺寸与当前文件不一致。"));
    auto rows = obj.value("crystals").toArray();
    if (rows.size() != 75)
        return fail(QStringLiteral("方案必须包含 75 个晶体。"));
    CrystalRois next;
    for (int i = 0; i < 75; ++i) {
        auto r = rows[i].toObject();
        auto &target = next.m_items[i];
        if (r.value("roi_label") != target.label || r.value("detector") != target.detector)
            return fail(QStringLiteral("晶体编号或探测器对应关系错误。"));
        QString state = r.value("state").toString();
        if (state != "pending" && state != "valid" && state != "invalid")
            return fail(QStringLiteral("ROI 状态无效。"));
        QVector<int> v;
        for (const auto &key : {"x1", "x2", "y1", "y2"}) {
            auto val = r.value(key);
            double d = val.toDouble(-1);
            if (!val.isDouble() || d < 0 || d > 100000 || std::floor(d) != d)
                return fail(QStringLiteral("ROI 坐标必须为非负整数。"));
            v.append(int(d));
        }
        QRect rect(v[0], v[2], v[1] - v[0], v[3] - v[2]);
        if (rect.isNull() && (v[0] || v[1] || v[2] || v[3]))
            return fail(QStringLiteral("空 ROI 必须使用全零坐标。"));
        if (!rect.isNull() && !next.assign(i, rect, sizes.value(i < 60 ? 0 : 1)))
            return fail(QStringLiteral("ROI 坐标超出探测器范围或为空。"));
        if (state == "valid" && rect.isNull())
            return fail(QStringLiteral("有效晶体必须具有非空 ROI。"));
        target.rect = rect;
        target.state = state == "valid"     ? RoiState::Valid
                       : state == "invalid" ? RoiState::Invalid
                                            : RoiState::Pending;
    }
    *this = next;
    return true;
}
