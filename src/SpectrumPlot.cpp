#include "SpectrumPlot.h"
#include "PlotStyle.h"
#include <QContextMenuEvent>
#include <QMenu>
#include <QMouseEvent>
#include <QPainter>
#include <QPainterPath>
#include <QToolTip>
#include <QWheelEvent>
#include <algorithm>
#include <cmath>

SpectrumPlot::SpectrumPlot(QWidget *parent) : QWidget(parent) {
    setMinimumSize(300, 220);
    setMouseTracking(true);
}
QRectF SpectrumPlot::area() const {
    int legendRows = std::min(3, (int(m_curves.size()) + 1) / 2);
    int top = 40 + legendRows * 19;
    return QRectF(82, top, std::max(1, width() - 108), std::max(1, height() - top - 65));
}
void SpectrumPlot::setCurves(const QVector<SpectrumCurve> &curves, const QString &title,
                             const QString &xLabel, const QString &yLabel) {
    m_curves = curves;
    m_title = title;
    m_xLabel = xLabel;
    m_yLabel = yLabel;
    double x0 = INFINITY, x1 = -INFINITY, y0 = INFINITY, y1 = -INFINITY;
    for (const auto &c : curves)
        for (qsizetype i = 0; i < std::min(c.x.size(), c.y.size()); ++i)
            if (std::isfinite(c.x[i]) && std::isfinite(c.y[i])) {
                x0 = std::min(x0, c.x[i]);
                x1 = std::max(x1, c.x[i]);
                y0 = std::min(y0, c.y[i]);
                y1 = std::max(y1, c.y[i]);
            }
    m_hasData = std::isfinite(x0);
    if (!m_hasData) {
        x0 = 0;
        x1 = 1;
        y0 = 0;
        y1 = 1;
    }
    if (x1 <= x0) {
        x0 -= .5;
        x1 += .5;
    }
    if (y1 <= y0) {
        y0 -= .5;
        y1 += .5;
    }
    double margin = (y1 - y0) * .06;
    m_extent = {x0, y0 - margin, x1 - x0, y1 - y0 + 2 * margin};
    m_view = m_extent;
    update();
}
void SpectrumPlot::paintEvent(QPaintEvent *) {
    QPainter p(this);
    p.fillRect(rect(), Qt::white);
    p.setRenderHint(QPainter::Antialiasing);
    auto a = area();
    QFont f("Arial", 9);
    p.setFont(f);
    p.setPen(Qt::black);
    p.drawText(QRectF(5, 8, width() - 10, 24), Qt::AlignCenter, m_title);
    if (!m_hasData) {
        p.drawText(a, Qt::AlignCenter,
                   m_curves.isEmpty() ? tr("完成计算后显示曲线")
                                      : tr("所选晶体没有有效曲线，请查看下方状态或选择其他晶体"));
        return;
    }
    auto point = [&](double x, double y) {
        return QPointF(a.left() + (x - m_view.left()) / m_view.width() * a.width(),
                       a.bottom() - (y - m_view.top()) / m_view.height() * a.height());
    };
    PlotStyle::axis(p, a, PlotStyle::ticks(m_view.left(), m_view.right(), std::max(2, int(a.width() / 100))),
                    m_view.left(), m_view.right(), true, false, m_grid, m_minorTicks);
    PlotStyle::axis(p, a, PlotStyle::ticks(m_view.top(), m_view.bottom(), std::max(2, int(a.height() / 60))),
                    m_view.top(), m_view.bottom(), false, false, m_grid, m_minorTicks);
    p.save();
    p.setClipRect(a);
    for (int k = 0; k < m_curves.size(); ++k) {
        const auto &c = m_curves[k];
        QColor color = c.color.isValid() ? c.color : PlotStyle::curveColor(k);
        QPen pen(color, 1.4, c.dashed ? Qt::DashLine : Qt::SolidLine);
        p.setPen(pen);
        QPainterPath path;
        bool begin = true;
        for (qsizetype i = 0; i < std::min(c.x.size(), c.y.size()); ++i) {
            if (!std::isfinite(c.x[i]) || !std::isfinite(c.y[i])) {
                begin = true;
                continue;
            }
            auto pos = point(c.x[i], c.y[i]);
            if (begin) {
                path.moveTo(pos);
                begin = false;
            } else
                path.lineTo(pos);
        }
        p.drawPath(path);
    }
    p.restore();
    p.setPen(Qt::black);
    p.drawRect(a);
    p.drawText(QRectF(a.left(), height() - 28, a.width(), 22), Qt::AlignCenter, m_xLabel);
    p.save();
    p.translate(16, a.center().y());
    p.rotate(-90);
    p.drawText(QRectF(-a.height() / 2, -12, a.height(), 24), Qt::AlignCenter, m_yLabel);
    p.restore();
    int visible = std::min(6, int(m_curves.size()));
    for (int k = 0; k < visible; ++k) {
        const auto &c = m_curves[k];
        QColor color = c.color.isValid() ? c.color : PlotStyle::curveColor(k);
        double x = a.left() + (k % 2) * a.width() / 2;
        double y = 39 + (k / 2) * 19;
        p.setPen(QPen(color, 1.5, c.dashed ? Qt::DashLine : Qt::SolidLine));
        p.drawLine(QPointF(x, y), QPointF(x + 23, y));
        p.setPen(Qt::black);
        if (m_curves.size() > visible && k == visible - 1) {
            p.drawText(QRectF(x, y - 9, a.width() / 2 - 5, 18), Qt::AlignVCenter,
                       tr("另有 %1 条曲线").arg(m_curves.size() - visible + 1));
            continue;
        }
        p.drawText(QRectF(x + 30, y - 9, a.width() / 2 - 35, 18), Qt::AlignVCenter,
                   fontMetrics().elidedText(c.name, Qt::ElideRight, int(a.width() / 2 - 35)));
    }
}
void SpectrumPlot::wheelEvent(QWheelEvent *e) {
    if (!m_hasData || !area().contains(e->position()))
        return;
    auto a = area();
    double fx = (e->position().x() - a.left()) / a.width(),
           fy = (a.bottom() - e->position().y()) / a.height(), factor = e->angleDelta().y() > 0 ? .8 : 1.25;
    double x = m_view.left() + fx * m_view.width(), y = m_view.top() + fy * m_view.height();
    QRectF next(x - fx * m_view.width() * factor, y - fy * m_view.height() * factor, m_view.width() * factor,
                m_view.height() * factor);
    if (next.width() > m_extent.width() * 1e-8 && next.height() > m_extent.height() * 1e-8 &&
        next.width() < m_extent.width() * 1e4 && next.height() < m_extent.height() * 1e4)
        m_view = next;
    update();
    e->accept();
}
void SpectrumPlot::mouseDoubleClickEvent(QMouseEvent *) {
    m_view = m_extent;
    update();
}
void SpectrumPlot::mouseMoveEvent(QMouseEvent *e) {
    auto a = area();
    if (a.contains(e->position()) && m_hasData) {
        double x = m_view.left() + (e->position().x() - a.left()) / a.width() * m_view.width(),
               y = m_view.top() + (a.bottom() - e->position().y()) / a.height() * m_view.height();
        QToolTip::showText(e->globalPosition().toPoint(),
                           QString("x = %1\ny = %2").arg(x, 0, 'g', 9).arg(y, 0, 'g', 9), this);
    } else
        QToolTip::hideText();
}

void SpectrumPlot::contextMenuEvent(QContextMenuEvent *e) {
    QMenu menu(this);
    menu.addAction(tr("复位坐标范围"), this, [this] {
        m_view = m_extent;
        update();
    });
    auto grid = menu.addAction(tr("主刻度网格"));
    grid->setCheckable(true);
    grid->setChecked(m_grid);
    connect(grid, &QAction::toggled, this, [this](bool v) {
        m_grid = v;
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
