#include "ScientificPlot.h"
#include "PlotStyle.h"
#include <QContextMenuEvent>
#include <QKeyEvent>
#include <QMenu>
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
} // namespace
ScientificPlot::ScientificPlot(QWidget *parent) : QWidget(parent) {
    setMinimumSize(360, 320);
    setMouseTracking(true);
    setToolTip(tr("滚轮缩放 · 双击复位"));
    setAutoFillBackground(true);
}
void ScientificPlot::setData(const DataSlice &data) {
    m_profileDrag = ProfileDrag::None;
    unsetCursor();
    m_dragging = false;
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
    m_profileDrag = ProfileDrag::None;
    m_dragging = false;
    unsetCursor();
    m_view = QRectF(0, 0, std::max(1, m_data.width), std::max(1, m_data.height));
    update();
}
void ScientificPlot::setProfilePosition(const QPoint &position) {
    m_profilePosition = position;
    update();
}
void ScientificPlot::setAxisLabels(const QString &x, const QString &y) {
    m_xLabel = x;
    m_yLabel = y;
    update();
}
void ScientificPlot::setImageCoordinates(double xFirst, double xStep, double yFirst, double yStep) {
    m_physicalAxes = std::isfinite(xFirst) && std::isfinite(yFirst) && std::isfinite(xStep) &&
                     std::isfinite(yStep) && xStep > 0 && yStep > 0;
    m_xOrigin = xFirst - .5 * xStep;
    m_yOrigin = yFirst - .5 * yStep;
    m_xStep = xStep;
    m_yStep = yStep;
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
    QRectF r = b.adjusted(82, 57, -130, -73);
    if (m_data.yAxis >= 0 && m_aspect) {
        double w = m_view.width() * (m_physicalAxes ? m_xStep : 1),
               h = m_view.height() * (m_physicalAxes ? m_yStep : 1);
        double scale = std::min(r.width() / w, r.height() / h);
        QSizeF size(w * scale, h * scale);
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
    // Overlays also belong in read-only mask/scan previews; selection permission
    // controls mouse editing independently.
    if (image && (!m_rois.isEmpty() || m_dragging)) {
        auto screenRect = [&](const QRectF &roi) {
            return QRectF(r.left() + (roi.x() - m_view.left()) / m_view.width() * r.width(),
                          r.top() + (roi.y() - m_view.top()) / m_view.height() * r.height(),
                          roi.width() / m_view.width() * r.width(),
                          roi.height() / m_view.height() * r.height());
        };
        for (const auto &roi : m_rois) {
            if (roi.rect.isEmpty() || roi.state == RoiState::Pending)
                continue;
            QColor ink = roi.label == m_selectedRoi       ? QColor("#00ffff")
                         : roi.state == RoiState::Invalid ? QColor("#999999")
                                                          : QColor("#ff655f");
            p.setPen(QPen(ink, roi.label == m_selectedRoi ? 2 : 1,
                          roi.state == RoiState::Invalid ? Qt::DashLine : Qt::SolidLine));
            p.setBrush(Qt::NoBrush);
            auto sr = screenRect(roi.rect);
            p.drawRect(sr);
            p.setFont(QFont("Arial", 9));
            p.drawText(sr.topLeft() + QPointF(2, -3), roi.label);
        }
        if (m_dragging) {
            p.setPen(QPen(QColor("#00ffff"), 1.5, Qt::DashLine));
            p.setBrush(QColor(0, 255, 255, 25));
            p.drawRect(screenRect(QRectF(m_dragStart, m_dragEnd).normalized()));
            p.setBrush(Qt::NoBrush);
        }
    }
    if (image && m_profilePosition.x() >= 0 && m_profilePosition.x() < m_data.width &&
        m_profilePosition.y() >= 0 && m_profilePosition.y() < m_data.height) {
        double x = r.left() + (m_profilePosition.x() + 0.5 - m_view.left()) / m_view.width() * r.width();
        double y = r.top() + (m_profilePosition.y() + 0.5 - m_view.top()) / m_view.height() * r.height();
        p.setPen(QPen(QColor("#00bfff"), 1.3, Qt::DashLine));
        p.drawLine(QPointF(x, r.top()), QPointF(x, r.bottom()));
        p.drawLine(QPointF(r.left(), y), QPointF(r.right(), y));
    }
    p.restore();
    p.setPen(QPen(Qt::black, 0.9));
    p.setBrush(Qt::NoBrush);
    p.drawRect(r);
    double sx = m_physicalAxes && image
                    ? m_xStep
                    : double(m_data.xAxis >= 0 ? m_data.strides.value(m_data.xAxis, 1) : 1),
           sy = m_physicalAxes && image
                    ? m_yStep
                    : double(m_data.yAxis >= 0 ? m_data.strides.value(m_data.yAxis, 1) : 1);
    const double xoffset = m_physicalAxes && image ? m_xOrigin : double(m_data.starts.value(m_data.xAxis, 0)),
                 yoffset = m_physicalAxes && image ? m_yOrigin : double(m_data.starts.value(m_data.yAxis, 0));
    const double xlo = xoffset + m_view.left() * sx, xhi = xoffset + m_view.right() * sx;
    const double ylo = image ? yoffset + m_view.top() * sy : inv(ymin),
                 yhi = image ? yoffset + m_view.bottom() * sy : inv(ymax);
    PlotStyle::axis(p, r, PlotStyle::ticks(xlo, xhi, std::max(2, int(r.width() / 100))), xlo, xhi, true,
                    false, image ? m_imageGrid : m_curveGrid, m_minorTicks);
    PlotStyle::axis(p, r, PlotStyle::ticks(ylo, yhi, std::max(2, int(r.height() / 65)), !image && m_log), ylo,
                    yhi, false, image, image ? m_imageGrid : m_curveGrid, m_minorTicks);
    p.drawText(QRectF(r.left(), r.bottom() + 36, r.width(), 22), Qt::AlignCenter,
               !m_xLabel.isEmpty() ? m_xLabel
               : m_data.xAxis < 0  ? "Scalar"
                                   : QString("Dimension %1 / index").arg(m_data.xAxis));
    p.save();
    p.translate(b.left() + 18, r.center().y());
    p.rotate(-90);
    p.drawText(QRectF(-r.height() / 2, -10, r.height(), 20), Qt::AlignCenter,
               !m_yLabel.isEmpty() ? m_yLabel
               : image             ? QString("Dimension %1 / index").arg(m_data.yAxis)
                                   : (m_log ? "Signal (a.u.; signed log1p scale)" : "Signal (a.u.)"));
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
            p.drawText(QRectF(bar.right() + 5, bar.center().y() - 10, 90, 20), Qt::AlignVCenter,
                       number(m_data.minimum));
        else {
            auto scale = PlotStyle::ticks(m_data.minimum, m_data.maximum, 4, m_log);
            double lo = tx(m_data.minimum), hi = tx(m_data.maximum);
            double previous = INFINITY;
            for (double v : scale.major) {
                double y = bar.bottom() - (tx(v) - lo) / (hi - lo) * bar.height();
                p.drawLine(QPointF(bar.left(), y), QPointF(bar.left() + 4, y));
                if (std::abs(y - previous) >= 24) {
                    p.drawText(QRectF(bar.right() + 5, y - 10, 90, 20), Qt::AlignVCenter,
                               PlotStyle::label(v, scale));
                    previous = y;
                }
            }
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
    if (m_profileDrag != ProfileDrag::None) {
        moveProfile(e->position());
        update();
        return;
    }
    auto hit = profileHit(e->position());
    if (hit == ProfileDrag::None && m_profilePicking &&
        (m_pickAnywhere || e->modifiers().testFlag(Qt::ControlModifier)) && m_data.yAxis >= 0 &&
        !m_data.values.isEmpty() && plotRect(rect()).contains(e->position()))
        hit = ProfileDrag::Both;
    if (hit == ProfileDrag::X)
        setCursor(Qt::SizeHorCursor);
    else if (hit == ProfileDrag::Y)
        setCursor(Qt::SizeVerCursor);
    else if (hit == ProfileDrag::Both)
        setCursor(Qt::SizeAllCursor);
    else
        unsetCursor();
    auto d = toData(m_hover);
    if (m_dragging)
        m_dragEnd = QPointF(std::clamp(d.x(), 0.0, double(m_data.width)),
                            std::clamp(d.y(), 0.0, double(m_data.height)));
    int x = int(std::floor(d.x())), y = m_data.yAxis >= 0 ? int(std::floor(d.y())) : 0;
    if (x >= 0 && x < m_data.width && y >= 0 && y < m_data.height && !m_data.values.isEmpty()) {
        if (m_physicalAxes && m_data.yAxis >= 0)
            emit cursorChanged(QString("%1 = %2 mm   %3 = %4 mm   value = %5")
                                   .arg(m_xLabel.section(' ', 0, 0))
                                   .arg(m_xOrigin + (x + .5) * m_xStep, 0, 'g', 8)
                                   .arg(m_yLabel.section(' ', 0, 0))
                                   .arg(m_yOrigin + (y + .5) * m_yStep, 0, 'g', 8)
                                   .arg(m_data.values[y * m_data.width + x], 0, 'g', 12));
        else
            emit cursorChanged(QString("x = %1   y = %2   value = %3")
                                   .arg(quint64(x) * m_data.strides.value(m_data.xAxis, 1))
                                   .arg(quint64(y) * m_data.strides.value(m_data.yAxis, 1))
                                   .arg(m_data.values[y * m_data.width + x], 0, 'g', 12));
    }
    update();
}
ScientificPlot::ProfileDrag ScientificPlot::profileHit(const QPointF &point) const {
    const auto r = plotRect(rect());
    if (!m_profilePicking || m_data.values.isEmpty() || m_data.yAxis < 0 || !r.contains(point) ||
        m_profilePosition.x() < 0 || m_profilePosition.y() < 0)
        return ProfileDrag::None;
    const double x = r.left() + (m_profilePosition.x() + 0.5 - m_view.left()) / m_view.width() * r.width();
    const double y = r.top() + (m_profilePosition.y() + 0.5 - m_view.top()) / m_view.height() * r.height();
    bool nearX = std::abs(point.x() - x) <= 6, nearY = std::abs(point.y() - y) <= 6;
    return nearX && nearY ? ProfileDrag::Both
           : nearX        ? ProfileDrag::X
           : nearY        ? ProfileDrag::Y
                          : ProfileDrag::None;
}
void ScientificPlot::moveProfile(const QPointF &point) {
    if (!m_profilePicking || m_data.values.isEmpty() || m_data.yAxis < 0)
        return;
    auto data = toData(point);
    QPoint next(std::clamp(int(std::floor(data.x())), 0, m_data.width - 1),
                std::clamp(int(std::floor(data.y())), 0, m_data.height - 1));
    if (m_profileDrag == ProfileDrag::X)
        next.setY(m_profilePosition.y());
    if (m_profileDrag == ProfileDrag::Y)
        next.setX(m_profilePosition.x());
    if (next != m_profilePosition)
        emit profilePointPicked(next);
}
void ScientificPlot::setRoiSelectionEnabled(bool enabled) {
    m_roiEnabled = enabled;
    m_dragging = false;
    setFocusPolicy(enabled ? Qt::StrongFocus : Qt::NoFocus);
    setToolTip(enabled ? tr("拖动框选当前晶体 · 滚轮缩放 · 双击复位 · Esc 取消框选")
                       : tr("滚轮缩放 · 双击复位"));
    update();
}
void ScientificPlot::setRoiOverlays(const QVector<CrystalRoi> &rois, const QString &selected) {
    m_rois = rois;
    m_selectedRoi = selected;
    update();
}
void ScientificPlot::mousePressEvent(QMouseEvent *e) {
    if (e->button() == Qt::LeftButton && m_profilePicking && m_data.yAxis >= 0 && !m_data.values.isEmpty() &&
        plotRect(rect()).contains(e->position())) {
        auto hit =
            e->modifiers().testFlag(Qt::ControlModifier) ? ProfileDrag::Both : profileHit(e->position());
        if (m_pickAnywhere && hit == ProfileDrag::None)
            hit = ProfileDrag::Both;
        if (hit != ProfileDrag::None) {
            setFocus();
            m_dragging = false;
            m_profileDrag = hit;
            setCursor(hit == ProfileDrag::X   ? Qt::SizeHorCursor
                      : hit == ProfileDrag::Y ? Qt::SizeVerCursor
                                              : Qt::SizeAllCursor);
            moveProfile(e->position());
            e->accept();
            return;
        }
    }
    if (m_roiEnabled && !m_data.values.isEmpty() && m_data.yAxis >= 0 && e->button() == Qt::LeftButton &&
        plotRect(rect()).contains(e->position())) {
        setFocus();
        m_dragging = true;
        m_screenStart = e->position();
        m_dragStart = m_dragEnd = toData(e->position());
        update();
    } else
        QWidget::mousePressEvent(e);
}
void ScientificPlot::mouseReleaseEvent(QMouseEvent *e) {
    if (e->button() == Qt::LeftButton && m_profileDrag != ProfileDrag::None) {
        moveProfile(e->position());
        m_profileDrag = ProfileDrag::None;
        unsetCursor();
        e->accept();
        return;
    }
    if (!m_dragging || e->button() != Qt::LeftButton)
        return;
    m_dragging = false;
    auto end = toData(e->position());
    auto area = QRectF(m_dragStart, end).normalized();
    int x1 = std::max(0, int(std::floor(area.left()))), y1 = std::max(0, int(std::floor(area.top())));
    int x2 = std::min(m_data.width, int(std::ceil(area.right()))),
        y2 = std::min(m_data.height, int(std::ceil(area.bottom())));
    if ((e->position() - m_screenStart).manhattanLength() >= 4 && x2 > x1 && y2 > y1)
        emit rectangleSelected(QRect(x1, y1, x2 - x1, y2 - y1));
    update();
}
void ScientificPlot::keyPressEvent(QKeyEvent *e) {
    if (e->key() == Qt::Key_Escape && m_profileDrag != ProfileDrag::None) {
        m_profileDrag = ProfileDrag::None;
        unsetCursor();
        e->accept();
        return;
    }
    if (e->key() == Qt::Key_Escape && m_dragging) {
        m_dragging = false;
        update();
        e->accept();
    } else
        QWidget::keyPressEvent(e);
}
void ScientificPlot::mouseDoubleClickEvent(QMouseEvent *) {
    resetView();
}
void ScientificPlot::wheelEvent(QWheelEvent *e) {
    m_profileDrag = ProfileDrag::None;
    unsetCursor();
    m_dragging = false;
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

void ScientificPlot::contextMenuEvent(QContextMenuEvent *e) {
    QMenu menu(this);
    menu.addAction(tr("复位坐标范围"), this, &ScientificPlot::resetView);
    const bool image = m_data.yAxis >= 0;
    auto grid = menu.addAction(tr("主刻度网格"));
    grid->setCheckable(true);
    grid->setChecked(image ? m_imageGrid : m_curveGrid);
    connect(grid, &QAction::toggled, this, [this, image](bool v) {
        (image ? m_imageGrid : m_curveGrid) = v;
        update();
    });
    auto minor = menu.addAction(tr("次刻度"));
    minor->setCheckable(true);
    minor->setChecked(m_minorTicks);
    connect(minor, &QAction::toggled, this, [this](bool v) {
        m_minorTicks = v;
        update();
    });
    menu.exec(e->globalPos());
}
