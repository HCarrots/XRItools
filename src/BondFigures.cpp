#include "BondFigures.h"
#include "PlotStyle.h"
#include <QMouseEvent>
#include <QPainterPath>
#include <QPolygonF>
#include <QWheelEvent>
#include <algorithm>

namespace {
QColor blend(QColor a, QColor b, double t) {
    return QColor::fromRgbF(a.redF() * (1 - t) + b.redF() * t, a.greenF() * (1 - t) + b.greenF() * t,
                            a.blueF() * (1 - t) + b.blueF() * t);
}
struct SpectralLimits {
    double lo, hi;
    int peakRow;
};
SpectralLimits limits(const QVector<QVector<double>> &rows) {
    double lo = INFINITY, hi = -INFINITY;
    double peak = -1.;
    int peakRow = -1;
    for (int j = 0; j < rows.size(); ++j)
        for (double value : rows[j])
            if (std::isfinite(value)) {
                lo = std::min(lo, value);
                hi = std::max(hi, value);
                if (std::abs(value) > peak) {
                    peak = std::abs(value);
                    peakRow = j;
                }
            }
    return {lo, hi > lo ? hi : lo + std::max(1e-10, std::abs(lo) * 1e-6), peakRow};
}
} // namespace

BondFigures::BondFigures(QWidget *parent) : QWidget(parent) {
    setMinimumSize(390, 300);
}

QColor BondFigures::fieldColor(double fraction, bool signedContrast) {
    const double f = std::clamp(fraction, 0., 1.);
    if (signedContrast)
        return f < .5 ? blend(QColor("#2166ac"), QColor("#f7f7f7"), 2 * f)
                      : blend(QColor("#f7f7f7"), QColor("#b2182b"), 2 * f - 1);
    return f < .5 ? blend(QColor("#440154"), QColor("#21918c"), 2 * f)
                  : blend(QColor("#21918c"), QColor("#fde725"), 2 * f - 1);
}

void BondFigures::setMontage(const QVector<QVector<double>> &axes, const QVector<double> &field,
                             const QVector<int> &zIndices, double minimum, double maximum,
                             bool signedContrast, const QString &scaleLabel) {
    m_kind = Kind::Montage;
    m_axes = axes;
    m_field = field;
    m_z = zIndices;
    m_min = minimum;
    m_max = maximum;
    m_signed = signedContrast;
    m_label = scaleLabel;
    update();
}

void BondFigures::setSpectralSection(SpectralStyle style, const QVector<double> &energy,
                                     const QVector<double> &x, const QVector<QVector<double>> &rows,
                                     const QString &energyLabel, const QString &intensityLabel,
                                     const QString &sectionLabel) {
    m_kind = Kind::Spectral;
    m_spectralStyle = style;
    m_axes = {energy, x};
    m_rows = rows;
    m_label = energyLabel;
    m_caption = intensityLabel;
    m_sectionLabel = sectionLabel;
    setToolTip(tr("左键拖动旋转，滚轮缩放，双击复位"));
    update();
}

void BondFigures::setMesh(const QVector<QVector<double>> &axes, const QVector<SurfaceTriangle> &first,
                          const QVector<SurfaceTriangle> &second, const QString &caption,
                          const std::array<QColor, 2> &colors) {
    m_kind = Kind::Mesh;
    m_axes = axes;
    m_mesh[0] = first;
    m_mesh[1] = second;
    m_meshColors = colors;
    m_caption = caption;
    update();
}

void BondFigures::paintEvent(QPaintEvent *) {
    QPainter p(this);
    p.fillRect(rect(), Qt::white);
    p.setRenderHint(QPainter::Antialiasing);
    p.setFont(QFont("Arial", 9));
    if (m_kind == Kind::Montage)
        drawMontage(p);
    else if (m_kind == Kind::Spectral)
        drawSpectralSection(p);
    else
        drawMesh(p);
}

void BondFigures::drawMontage(QPainter &p) {
    if (m_field.isEmpty() || m_z.isEmpty()) {
        p.drawText(rect(), Qt::AlignCenter, tr("当前能量没有可显示的切片"));
        return;
    }
    const auto &x = m_axes[0], &y = m_axes[1], &z = m_axes[2];
    const double dx = x[1] - x[0], dy = y[1] - y[0];
    const double x0 = x.first() - dx / 2, x1 = x.last() + dx / 2;
    const double y0 = y.first() - dy / 2, y1 = y.last() + dy / 2;
    const int columns = std::min({6, std::max(1, (width() - 104) / 115), int(m_z.size())});
    const int rows = (int(m_z.size()) + columns - 1) / columns;
    const double cellW = std::min(150., (width() - 104.) / columns), cellH = (height() - 32.) / rows;
    const double left = std::max(10., (width() - 104. - columns * cellW) / 2);
    for (int i = 0; i < m_z.size(); ++i) {
        double w = cellW - 42, h = cellH - 58;
        if (w <= 0 || h <= 0)
            return;
        const double ratio = (y1 - y0) / (x1 - x0);
        w = std::min(w, h * ratio);
        h = w / ratio;
        QRectF area(left + (i % columns) * cellW + (cellW - w) / 2,
                    6 + (i / columns) * cellH + (cellH - h) / 2, w, h);
        QImage image(int(y.size()), int(x.size()), QImage::Format_RGB32);
        for (int j = 0; j < y.size(); ++j)
            for (int k = 0; k < x.size(); ++k) {
                const auto v = m_field[(qsizetype(k) * y.size() + j) * z.size() + m_z[i]];
                image.setPixelColor(j, k,
                                    std::isfinite(v) ? fieldColor((v - m_min) / (m_max - m_min), m_signed)
                                                     : QColor("#dddddd"));
            }
        p.drawImage(area, image);
        p.setPen(Qt::black);
        p.drawRect(area);
        PlotStyle::axis(p, area, PlotStyle::ticks(y0, y1, 3), y0, y1, true, false, false, true, 8);
        PlotStyle::axis(p, area, PlotStyle::ticks(x0, x1, 5), x0, x1, false, true, false, true, 8);
        p.drawText(QRectF(area.center().x() - cellW / 2, area.top() - 24, cellW, 24), Qt::AlignCenter,
                   QStringLiteral("Z = %1 mm").arg(z[m_z[i]], 0, 'g', 5));
    }
    QRectF bar(width() - 75., 45, 12, height() - 98.);
    QLinearGradient gradient(bar.bottomLeft(), bar.topLeft());
    for (int i = 0; i <= 10; ++i)
        gradient.setColorAt(i / 10., fieldColor(i / 10., m_signed));
    p.fillRect(bar, gradient);
    p.setPen(Qt::black);
    p.drawRect(bar);
    auto t = PlotStyle::ticks(m_min, m_max, 6);
    for (double value : t.major) {
        double tickY = bar.bottom() - bar.height() * (value - m_min) / (m_max - m_min);
        p.drawLine(QPointF(bar.right(), tickY), QPointF(bar.right() + 4, tickY));
        p.drawText(QRectF(bar.right() + 6, tickY - 8, 56, 16), Qt::AlignVCenter, PlotStyle::label(value, t));
    }
    p.drawText(QRectF(15, height() - 32, width() - 30, 30), Qt::AlignCenter,
               tr("Y（横向）和 X（纵向） / mm    ·    %1").arg(m_label));
}

void BondFigures::drawSpectralSection(QPainter &p) {
    if (m_rows.isEmpty()) {
        p.drawText(rect(), Qt::AlignCenter, tr("载入四维结果后显示三维谱图"));
        return;
    }
    const auto &e = m_axes[0], &x = m_axes[1];
    const auto spectrumLimits = limits(m_rows);
    const double lo = spectrumLimits.lo, hi = spectrumLimits.hi;
    const int peakRow = spectrumLimits.peakRow;
    if (peakRow < 0) {
        p.drawText(rect(), Qt::AlignCenter, tr("当前 Y/Z 位置没有可显示的谱点"));
        return;
    }
    constexpr double ew = 1.65, xw = 1.12, iw = .9;
    const double ca = std::cos(m_yaw), sa = std::sin(m_yaw);
    const double cp = std::cos(m_pitch), sp = std::sin(m_pitch);
    const double projectedW = std::abs(ca) * ew + std::abs(sa) * xw;
    const double projectedH = std::abs(sp) * (std::abs(sa) * ew + std::abs(ca) * xw) + std::abs(cp) * iw;
    const double scale = m_zoom * std::min((width() - 182.) / projectedW, (height() - 140.) / projectedH);
    const QPointF center(width() / 2., height() / 2.);
    struct Projected {
        QPointF screen;
        double depth;
    };
    auto project = [&](double u, double v, double w) -> Projected {
        u = (u - .5) * ew;
        v = (v - .5) * xw;
        w = (w - .5) * iw;
        const double forward = sa * u + ca * v;
        return {center + QPointF(scale * (ca * u - sa * v), scale * (sp * forward - cp * w)),
                cp * forward + sp * w};
    };
    const double farE = sa < 0 ? 1. : 0., farX = ca < 0 ? 1. : 0.;
    const auto et = PlotStyle::ticks(e.first(), e.last(), 5);
    const auto xt = PlotStyle::ticks(x.first(), x.last(), 5);
    const auto it = PlotStyle::ticks(lo, hi, 5);
    auto energyFraction = [&](double value) { return (value - e.first()) / (e.last() - e.first()); };
    auto positionFraction = [&](double value) { return (value - x.first()) / (x.last() - x.first()); };
    auto intensityFraction = [&](double value) { return (value - lo) / (hi - lo); };
    p.setPen(QPen(QColor("#d0d4da"), .8, Qt::DashLine));
    for (double value : et.major) {
        const double u = energyFraction(value);
        p.drawLine(project(u, 0, 0).screen, project(u, 1, 0).screen);
        p.drawLine(project(u, farX, 0).screen, project(u, farX, 1).screen);
    }
    for (double value : xt.major) {
        const double v = positionFraction(value);
        p.drawLine(project(0, v, 0).screen, project(1, v, 0).screen);
        p.drawLine(project(farE, v, 0).screen, project(farE, v, 1).screen);
    }
    for (double value : it.major) {
        const double w = intensityFraction(value);
        p.drawLine(project(0, farX, w).screen, project(1, farX, w).screen);
        p.drawLine(project(farE, 0, w).screen, project(farE, 1, w).screen);
    }

    auto vertex = [&](int j, int i) {
        return project(energyFraction(e[i]), positionFraction(x[j]), intensityFraction(m_rows[j][i]));
    };
    if (m_spectralStyle == SpectralStyle::Surface) {
        struct Face {
            std::array<QPointF, 4> points;
            double depth, color;
        };
        std::vector<Face> faces;
        faces.reserve((e.size() - 1) * (x.size() - 1));
        for (int j = 0; j + 1 < x.size(); ++j)
            for (int i = 0; i + 1 < e.size(); ++i) {
                const double a = m_rows[j][i], b = m_rows[j][i + 1], c = m_rows[j + 1][i + 1],
                             d = m_rows[j + 1][i];
                if (!std::isfinite(a) || !std::isfinite(b) || !std::isfinite(c) || !std::isfinite(d))
                    continue;
                const auto p0 = vertex(j, i), p1 = vertex(j, i + 1), p2 = vertex(j + 1, i + 1),
                           p3 = vertex(j + 1, i);
                faces.push_back({{p0.screen, p1.screen, p2.screen, p3.screen},
                                 p0.depth + p1.depth + p2.depth + p3.depth,
                                 intensityFraction((a + b + c + d) / 4)});
            }
        std::sort(faces.begin(), faces.end(), [](const Face &a, const Face &b) { return a.depth < b.depth; });
        for (const auto &face : faces) {
            const auto color = fieldColor(face.color, false);
            p.setBrush(color);
            p.setPen(QPen(color.darker(120), .42));
            p.drawPolygon(face.points.data(), int(face.points.size()));
        }
        if (faces.empty()) {
            p.setPen(Qt::NoPen);
            for (int j = 0; j < x.size(); ++j)
                for (int i = 0; i < e.size(); ++i)
                    if (std::isfinite(m_rows[j][i])) {
                        p.setBrush(fieldColor(intensityFraction(m_rows[j][i]), false));
                        p.drawEllipse(vertex(j, i).screen, 2, 2);
                    }
        }
        p.setBrush(Qt::NoBrush);
    } else {
        const int count = std::min(24, int(x.size()));
        QVector<int> displayed;
        for (int k = 0; k < count; ++k)
            displayed.append(int(k * (x.size() - 1) / (count - 1)));
        if (!displayed.contains(peakRow)) {
            auto nearest =
                std::min_element(displayed.begin() + 1, displayed.end() - 1, [peakRow](int a, int b) {
                    return std::abs(a - peakRow) < std::abs(b - peakRow);
                });
            *nearest = peakRow;
            std::sort(displayed.begin(), displayed.end());
        }
        if (farX > .5)
            std::reverse(displayed.begin(), displayed.end());
        for (int j : displayed) {
            p.setPen(QPen(fieldColor(positionFraction(x[j]), false), 1.15));
            QPainterPath line;
            bool connected = false;
            for (int i = 0; i < e.size(); ++i) {
                if (!std::isfinite(m_rows[j][i])) {
                    connected = false;
                    continue;
                }
                const auto point = vertex(j, i).screen;
                if (!connected && (i + 1 == e.size() || !std::isfinite(m_rows[j][i + 1])))
                    p.drawPoint(point);
                if (connected)
                    line.lineTo(point);
                else
                    line.moveTo(point);
                connected = true;
            }
            p.drawPath(line);
        }
    }

    auto axis = [&](QPointF start, QPointF end, const PlotStyle::Ticks &ticks, auto fraction) {
        p.setPen(QColor("#333940"));
        p.drawLine(start, end);
        QPointF normal = (start + end) / 2 - center;
        const double length = std::hypot(normal.x(), normal.y());
        if (length > 0)
            normal /= length;
        for (double value : ticks.major) {
            const auto at = start + fraction(value) * (end - start);
            p.drawLine(at, at + 4 * normal);
            if (QLineF(start, end).length() >= 82) {
                const auto text = PlotStyle::label(value, ticks);
                const double offset = 11 + 22 * std::abs(normal.x());
                const auto label = at + offset * normal;
                p.drawText(QRectF(label.x() - 38, label.y() - 8, 76, 16), Qt::AlignCenter, text);
            }
        }
    };
    axis(project(0, 1 - farX, 0).screen, project(1, 1 - farX, 0).screen, et, energyFraction);
    axis(project(1 - farE, 0, 0).screen, project(1 - farE, 1, 0).screen, xt, positionFraction);
    const auto left = project(farE, 1 - farX, 0).screen.x() < project(1 - farE, farX, 0).screen.x();
    const double vertE = left ? farE : 1 - farE, vertX = left ? 1 - farX : farX;
    axis(project(vertE, vertX, 0).screen, project(vertE, vertX, 1).screen, it, intensityFraction);

    p.setPen(Qt::black);
    p.drawText(QRectF(8, 8, width() - 16, 24), Qt::AlignCenter, m_sectionLabel);
    p.drawText(QRectF(8, height() - 45, width() - 16, 21), Qt::AlignCenter,
               tr("%1（能量轴）    X (mm)（位置轴）    %2（高度）").arg(m_label, m_caption));
    p.setPen(QColor("#68717a"));
    p.drawText(QRectF(8, height() - 23, width() - 16, 18), Qt::AlignCenter,
               tr("左键拖动旋转 · 滚轮缩放 · 双击复位"));
}

void BondFigures::drawMesh(QPainter &p) {
    if (m_axes.size() != 3) {
        p.drawText(rect(), Qt::AlignCenter, m_caption);
        return;
    }
    const auto &x = m_axes[0], &y = m_axes[1], &z = m_axes[2];
    const QVector3D lo(x.first(), y.first(), z.first()), hi(x.last(), y.last(), z.last()),
        center = (lo + hi) / 2;
    const double ca = std::cos(m_yaw), sa = std::sin(m_yaw);
    const double cp = std::cos(m_pitch), sp = std::sin(m_pitch);
    const auto span = hi - lo;
    const double projectedW = std::abs(ca) * span.x() + std::abs(sa) * span.y();
    const double projectedH =
        std::abs(sp) * (std::abs(sa) * span.x() + std::abs(ca) * span.y()) + std::abs(cp) * span.z();
    const double scale = m_zoom * std::min((width() - 160.) / projectedW, (height() - 110.) / projectedH);
    auto project = [&](const QVector3D &point) {
        const auto v = point - center;
        double horizontal = v.x() * ca - v.y() * sa, forward = v.x() * sa + v.y() * ca;
        return QVector3D(float(width() * .5 + horizontal * scale),
                         float(height() * .49 - (v.z() * cp - forward * sp) * scale),
                         float((v.z() * sp + forward * cp) * scale));
    };
    auto screen = [](const QVector3D &v) { return QPointF(v.x(), v.y()); };
    std::array<QVector3D, 8> box;
    for (int i = 0; i < 8; ++i)
        box[i] =
            project(QVector3D(i & 1 ? hi.x() : lo.x(), i & 2 ? hi.y() : lo.y(), i & 4 ? hi.z() : lo.z()));
    p.setPen(QPen(QColor("#aeb6be"), 1, Qt::DashLine));
    for (int i = 0; i < 8; ++i)
        for (int bit : {1, 2, 4})
            if (!(i & bit))
                p.drawLine(screen(box[i]), screen(box[i | bit]));
    struct Face {
        std::array<QVector3D, 3> projected;
        double depth;
        int layer;
    };
    std::vector<Face> faces;
    faces.reserve(m_mesh[0].size() + m_mesh[1].size());
    for (int layer = 0; layer < 2; ++layer)
        for (const auto &triangle : m_mesh[layer]) {
            Face f{{project(triangle[0]), project(triangle[1]), project(triangle[2])}, 0, layer};
            f.depth = f.projected[0].z() + f.projected[1].z() + f.projected[2].z();
            faces.push_back(f);
        }
    std::sort(faces.begin(), faces.end(), [](const Face &a, const Face &b) { return a.depth < b.depth; });
    p.setRenderHint(QPainter::Antialiasing, false);
    p.setPen(Qt::NoPen);
    for (const auto &face : faces) {
        const auto &q = face.projected;
        auto normal = QVector3D::normal(q[1] - q[0], q[2] - q[0]);
        const auto base = m_meshColors[face.layer];
        p.setBrush(blend(base.darker(160), base.lighter(155), std::abs(normal.z())));
        QPointF points[]{screen(q[0]), screen(q[1]), screen(q[2])};
        p.drawPolygon(points, 3);
    }
    p.setRenderHint(QPainter::Antialiasing);
    p.setPen(Qt::black);
    const QVector<QVector<double>> axes{x, y, z};
    for (int i = 0; i < 3; ++i) {
        const auto start = screen(box[0]), end = screen(box[std::size_t(1) << i]);
        p.drawLine(start, end);
        const auto &axis = axes[i];
        auto ticks = PlotStyle::ticks(axis.first(), axis.last(), 5);
        for (double v : ticks.major) {
            auto point = start + (end - start) * ((v - axis.first()) / (axis.last() - axis.first()));
            p.drawLine(point, point + QPointF(0, 4));
            if (QLineF(start, end).length() > 95)
                p.drawText(point + QPointF(3, 13), PlotStyle::label(v, ticks));
        }
        p.drawText(end + QPointF(5, -7), QStringLiteral("%1 / mm").arg(QStringLiteral("XYZ").at(i)));
    }
    p.drawText(QRectF(8, height() - 30, width() - 16, 24), Qt::AlignCenter, m_caption);
}

void BondFigures::mousePressEvent(QMouseEvent *event) {
    if (m_kind == Kind::Montage || event->button() != Qt::LeftButton) {
        QWidget::mousePressEvent(event);
        return;
    }
    m_lastMouse = event->position();
    setCursor(Qt::ClosedHandCursor);
    event->accept();
}

void BondFigures::mouseMoveEvent(QMouseEvent *event) {
    if (m_kind != Kind::Montage && (event->buttons() & Qt::LeftButton)) {
        const auto delta = event->position() - m_lastMouse;
        m_yaw += delta.x() * .012;
        m_pitch = std::clamp(m_pitch + delta.y() * .012, -1.5, 1.5);
        m_lastMouse = event->position();
        update();
        event->accept();
    }
}

void BondFigures::mouseReleaseEvent(QMouseEvent *event) {
    if (event->button() == Qt::LeftButton)
        unsetCursor();
    QWidget::mouseReleaseEvent(event);
}

void BondFigures::mouseDoubleClickEvent(QMouseEvent *event) {
    if (m_kind == Kind::Montage || event->button() != Qt::LeftButton) {
        QWidget::mouseDoubleClickEvent(event);
        return;
    }
    m_yaw = -.7;
    m_pitch = .46;
    m_zoom = 1.;
    m_lastMouse = event->position();
    unsetCursor();
    update();
    event->accept();
}

void BondFigures::wheelEvent(QWheelEvent *event) {
    if (m_kind == Kind::Montage || event->angleDelta().y() == 0) {
        QWidget::wheelEvent(event);
        return;
    }
    m_zoom = std::clamp(m_zoom * std::pow(1.15, event->angleDelta().y() / 120.), .55, 2.7);
    update();
    event->accept();
}
