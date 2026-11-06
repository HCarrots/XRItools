#include "ScientificPlot.h"
#include <QMouseEvent>
#include <QPainter>
#include <QPainterPath>
#include <QWheelEvent>
#include <algorithm>
#include <cmath>

namespace {
QString number(double v) {
    return QString::number(v, 'g', 5);
}
QVector<double> ticks(double lo, double hi) {
    QVector<double> result;
    if (!std::isfinite(lo) || !std::isfinite(hi) || hi <= lo)
        return result;
    const double rough = (hi - lo) / 5.0, magnitude = std::pow(10.0, std::floor(std::log10(rough)));
    if (!(magnitude > 0))
        return result;
    const double fraction = rough / magnitude;
    const double step = (fraction <= 1 ? 1 : fraction <= 2 ? 2 : fraction <= 5 ? 5 : 10) * magnitude;
    const double first = std::ceil(lo / step) * step;
    for (int i = 0; i < 12; ++i) {
        double v = first + i * step;
        if (v > hi + step * 1e-8)
            break;
        result.append(std::abs(v) < step * 1e-8 ? 0 : v);
    }
    return result;
}
} // namespace
ScientificPlot::ScientificPlot(QWidget *parent) : QWidget(parent) {
    setMinimumSize(360, 320);
    setMouseTracking(true);
    setToolTip(tr("滚轮缩放 · 双击复位"));
    setAutoFillBackground(true);
}
void ScientificPlot::setData(const DataSlice &data) {
    const bool reset = m_data.path != data.path || m_data.width != data.width ||
                       m_data.height != data.height || m_data.xAxis != data.xAxis ||
                       m_data.yAxis != data.yAxis;
    m_data = data;
    m_hover = {-1, -1};
    if (reset) {
        resetView();
    }
    rebuildImage();
    update();
}
void ScientificPlot::setLogarithmic(bool v) {
    m_log = v;
    rebuildImage();
    update();
}
void ScientificPlot::setPalette(int v) {
    m_palette = v;
    rebuildImage();
    update();
}
void ScientificPlot::setAspectLocked(bool v) {
    m_aspect = v;
    update();
}
void ScientificPlot::resetView() {
    m_view = QRectF(0, 0, std::max(1, m_data.width), std::max(1, m_data.height));
    update();
}
QColor ScientificPlot::color(double t) const {
    t = std::clamp(t, 0.0, 1.0);
    if (m_palette == 2)
        return QColor::fromRgbF(t, t, t);
    const QVector<QColor> viridis{QColor("#440154"), QColor("#3b528b"), QColor("#21918c"), QColor("#5ec962"),
                                  QColor("#fde725")};
    const QVector<QColor> inferno{QColor("#000004"), QColor("#57106e"), QColor("#bc3754"), QColor("#f98e09"),
                                  QColor("#fcffa4")};
    const auto &colors = m_palette == 1 ? inferno : viridis;
    double pos = t * 4;
    int i = std::min(3, int(pos));
    double a = pos - i;
    return QColor::fromRgbF(colors[i].redF() * (1 - a) + colors[i + 1].redF() * a,
                            colors[i].greenF() * (1 - a) + colors[i + 1].greenF() * a,
                            colors[i].blueF() * (1 - a) + colors[i + 1].blueF() * a);
}
void ScientificPlot::rebuildImage() {
    if (m_data.values.isEmpty() || m_data.yAxis < 0) {
        m_image = {};
        return;
    }
    m_image = QImage(m_data.width, m_data.height, QImage::Format_RGB32);
    double lo = m_data.minimum, hi = m_data.maximum;
    // Signed log1p preserves negative values and includes zero-count frames.
    auto transform = [this](double v) { return m_log ? std::copysign(std::log1p(std::abs(v)), v) : v; };
    lo = transform(lo);
    hi = transform(hi);
    QVector<QRgb> lut;
    for (int i = 0; i < 1024; ++i)
        lut.append(color(i / 1023.0).rgb());
    for (int y = 0; y < m_data.height; ++y) {
        auto row = reinterpret_cast<QRgb *>(m_image.scanLine(y));
        for (int x = 0; x < m_data.width; ++x) {
            double v = m_data.values[y * m_data.width + x];
            double t = hi > lo ? (transform(v) - lo) / (hi - lo) : 0;
            row[x] = std::isfinite(v) ? lut[std::clamp(int(t * 1023), 0, 1023)] : qRgb(210, 210, 210);
        }
    }
}
QRectF ScientificPlot::plotRect(const QRectF &b) const {
    QRectF r = b.adjusted(82, 57, -108, -73);
    if (m_data.yAxis >= 0 && m_aspect) {
        double scale = std::min(r.width() / m_view.width(), r.height() / m_view.height());
        QSizeF size(m_view.width() * scale, m_view.height() * scale);
        r = QRectF(r.center() - QPointF(size.width() / 2, size.height() / 2), size);
    }
    return r;
}
void ScientificPlot::paintEvent(QPaintEvent *) {
    QPainter p(this);
    draw(p, rect(), true);
}
void ScientificPlot::draw(QPainter &p, const QRectF &b, bool interaction) {
    p.fillRect(b, Qt::white);
    p.setRenderHint(QPainter::Antialiasing);
    if (m_data.values.isEmpty()) {
        p.setPen(QColor("#647080"));
        p.setFont(QFont("Microsoft YaHei UI", 11));
        p.drawText(b.adjusted(30, 30, -30, -30), Qt::AlignCenter | Qt::TextWordWrap,
                   m_data.error.isEmpty() ? (m_data.text.isEmpty() ? tr("选择数值数据集以预览") : m_data.text)
                                          : m_data.error);
        return;
    }
    const QRectF r = plotRect(b);
    p.setFont(QFont("Arial", 10));
    p.setPen(Qt::black);
    p.drawText(QRectF(b.left() + 20, b.top() + 16, b.width() - 40, 25), Qt::AlignCenter,
               m_data.path.section('/', -1));
    const bool image = m_data.yAxis >= 0;
    double ymin = m_data.minimum, ymax = m_data.maximum;
    auto tx = [this](double v) { return m_log ? std::copysign(std::log1p(std::abs(v)), v) : v; };
    auto inv = [this](double v) { return m_log ? std::copysign(std::expm1(std::abs(v)), v) : v; };
    ymin = tx(ymin);
    ymax = tx(ymax);
    if (!std::isfinite(ymin) || !std::isfinite(ymax)) {
        ymin = 0;
        ymax = 1;
    }
    if (ymin == ymax) {
        double delta = std::max(1.0, std::abs(ymin) * 0.05);
        ymin -= delta;
        ymax += delta;
    } else {
        double pad = (ymax - ymin) * 0.06;
        ymin -= pad;
        ymax += pad;
    }
    p.save();
    p.setClipRect(r);
    if (image) {
        p.setRenderHint(QPainter::SmoothPixmapTransform, false);
        p.drawImage(r, m_image, m_view);
    } else {
        p.setPen(QPen(QColor("#e3e7eb"), 0.6, Qt::DotLine));
        for (int i = 1; i < 5; ++i) {
            double y = r.bottom() - i * r.height() / 5;
            p.drawLine(QPointF(r.left(), y), QPointF(r.right(), y));
        }
        QPainterPath path;
        bool active = false;
        int first = std::max(0, int(m_view.left())),
            last = std::min(m_data.width, int(std::ceil(m_view.right())) + 1);
        for (int x = first; x < last; ++x) {
            double v = m_data.values[x];
            if (!std::isfinite(v)) {
                active = false;
                continue;
            }
            QPointF point(r.left() + (x - m_view.left()) / m_view.width() * r.width(),
                          r.bottom() - (tx(v) - ymin) / (ymax - ymin) * r.height());
            if (active)
                path.lineTo(point);
            else {
                path.moveTo(point);
                active = true;
            }
        }
        p.setPen(QPen(QColor("#195b9a"), 1.25));
        p.drawPath(path);
        if (m_data.width == 1 && std::isfinite(m_data.values[0])) {
            p.setBrush(QColor("#195b9a"));
            p.drawEllipse(QPointF(r.center().x(),
                                  r.bottom() - (tx(m_data.values[0]) - ymin) / (ymax - ymin) * r.height()),
                          3.5, 3.5);
            p.setBrush(Qt::NoBrush);
        }
    }
    if (interaction && m_hover.x() >= 0 && r.contains(m_hover)) {
        p.setPen(QPen(QColor(160, 160, 160, 180), 0.8, Qt::DashLine));
        p.drawLine(QPointF(m_hover.x(), r.top()), QPointF(m_hover.x(), r.bottom()));
        if (image)
            p.drawLine(QPointF(r.left(), m_hover.y()), QPointF(r.right(), m_hover.y()));
    }
    p.restore();
    p.setPen(QPen(Qt::black, 0.9));
    p.setBrush(Qt::NoBrush);
    p.drawRect(r);
    quint64 sx = m_data.xAxis >= 0 ? m_data.strides.value(m_data.xAxis, 1) : 1,
            sy = m_data.yAxis >= 0 ? m_data.strides.value(m_data.yAxis, 1) : 1;
    const double xoffset = double(m_data.starts.value(m_data.xAxis, 0)),
                 yoffset = double(m_data.starts.value(m_data.yAxis, 0));
    for (double label :
         ticks(xoffset + m_view.left() * double(sx), xoffset + (m_view.right() - 1) * double(sx))) {
        double x = r.left() + ((label - xoffset) / double(sx) - m_view.left()) / m_view.width() * r.width();
        p.drawLine(QPointF(x, r.bottom()), QPointF(x, r.bottom() - 5));
        p.drawLine(QPointF(x, r.top()), QPointF(x, r.top() + 5));
        p.drawText(QRectF(x - 40, r.bottom() + 8, 80, 19), Qt::AlignCenter, number(label));
    }
    const auto yticks =
        image ? ticks(yoffset + m_view.top() * double(sy), yoffset + (m_view.bottom() - 1) * double(sy))
              : ticks(ymin, ymax);
    for (double tick : yticks) {
        double y =
            image ? r.top() + ((tick - yoffset) / double(sy) - m_view.top()) / m_view.height() * r.height()
                  : r.bottom() - (tick - ymin) / (ymax - ymin) * r.height();
        p.drawLine(QPointF(r.left(), y), QPointF(r.left() + 5, y));
        p.drawLine(QPointF(r.right(), y), QPointF(r.right() - 5, y));
        double label = image ? tick : inv(tick);
        p.drawText(QRectF(b.left() + 2, y - 10, r.left() - b.left() - 11, 20),
                   Qt::AlignRight | Qt::AlignVCenter, number(label));
    }
    p.drawText(QRectF(r.left(), r.bottom() + 36, r.width(), 22), Qt::AlignCenter,
               m_data.xAxis < 0 ? "Scalar" : QString("Dimension %1 / index").arg(m_data.xAxis));
    p.save();
    p.translate(b.left() + 18, r.center().y());
    p.rotate(-90);
    p.drawText(QRectF(-r.height() / 2, -10, r.height(), 20), Qt::AlignCenter,
               image ? QString("Dimension %1 / index").arg(m_data.yAxis)
                     : (m_log ? "Signal / signed log1p" : "Signal / a.u."));
    p.restore();
    if (image) {
        QRectF bar(r.right() + 25, r.top(), 15, r.height());
        const bool constant = m_data.minimum == m_data.maximum;
        for (int i = 0; i < 256; ++i)
            p.fillRect(
                QRectF(bar.left(), bar.top() + i * bar.height() / 256, bar.width(), bar.height() / 256 + 1),
                m_data.finiteCount ? color(constant ? 0 : 1 - i / 255.0) : QColor(210, 210, 210));
        p.setPen(Qt::black);
        p.drawRect(bar);
        if (constant || !m_data.finiteCount)
            p.drawText(QRectF(bar.right() + 5, bar.center().y() - 10, 65, 20), Qt::AlignVCenter,
                       number(m_data.minimum));
        else
            for (int i = 0; i <= 4; ++i) {
                double t = i / 4.0, y = bar.bottom() - t * bar.height();
                double v = inv(tx(m_data.minimum) + t * (tx(m_data.maximum) - tx(m_data.minimum)));
                p.drawText(QRectF(bar.right() + 5, y - 10, 65, 20), Qt::AlignVCenter, number(v));
            }
        p.drawText(QRectF(bar.left() - 8, bar.bottom() + 10, 80, 20), m_log ? "signed log1p" : "a.u.");
    }
}
QPointF ScientificPlot::toData(const QPointF &pos) const {
    const auto r = plotRect(rect());
    return QPointF(m_view.left() + (pos.x() - r.left()) / r.width() * m_view.width(),
                   m_view.top() + (pos.y() - r.top()) / r.height() * m_view.height());
}
void ScientificPlot::mouseMoveEvent(QMouseEvent *e) {
    m_hover = e->position();
    auto d = toData(m_hover);
    int x = int(std::floor(d.x())), y = m_data.yAxis >= 0 ? int(std::floor(d.y())) : 0;
    if (x >= 0 && x < m_data.width && y >= 0 && y < m_data.height && !m_data.values.isEmpty())
        emit cursorChanged(QString("x = %1   y = %2   value = %3")
                               .arg(quint64(x) * m_data.strides.value(m_data.xAxis, 1))
                               .arg(quint64(y) * m_data.strides.value(m_data.yAxis, 1))
                               .arg(m_data.values[y * m_data.width + x], 0, 'g', 12));
    update();
}
void ScientificPlot::mouseDoubleClickEvent(QMouseEvent *) {
    resetView();
}
void ScientificPlot::wheelEvent(QWheelEvent *e) {
    if (m_data.values.isEmpty() || !plotRect(rect()).contains(e->position()))
        return;
    auto center = toData(e->position());
    double factor = e->angleDelta().y() > 0 ? 0.8 : 1.25;
    double w = std::clamp(m_view.width() * factor, std::min(4.0, double(m_data.width)), double(m_data.width));
    double h = m_data.yAxis >= 0 ? std::clamp(m_view.height() * factor, std::min(4.0, double(m_data.height)),
                                              double(m_data.height))
                                 : 1;
    double left = std::clamp(center.x() - (center.x() - m_view.left()) * w / m_view.width(), 0.0,
                             double(m_data.width) - w);
    double top = std::clamp(center.y() - (center.y() - m_view.top()) * h / m_view.height(), 0.0,
                            double(m_data.height) - h);
    m_view = QRectF(left, top, w, h);
    update();
    e->accept();
}
