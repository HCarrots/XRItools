#include "BondFigures.h"
#include "PlotStyle.h"
#include <QMouseEvent>
#include <QPainterPath>
#include <QPolygonF>
#include <QWheelEvent>
#include <algorithm>
#include <cmath>

namespace {
QColor blend(QColor a, QColor b, double t) {
    return QColor::fromRgbF(a.redF() * (1 - t) + b.redF() * t, a.greenF() * (1 - t) + b.greenF() * t,
                            a.blueF() * (1 - t) + b.blueF() * t);
}
struct MeshProjection {
    QVector3D center;
    QPointF origin;
    double scale, yaw, pitch;

    auto projector() const {
        const double ca = std::cos(yaw), sa = std::sin(yaw);
        const double cp = std::cos(pitch), sp = std::sin(pitch);
        return [view = *this, ca, sa, cp, sp](const QVector3D &point) {
            const auto v = point - view.center;
            const double horizontal = v.x() * ca - v.y() * sa, forward = v.x() * sa + v.y() * ca;
            return QVector3D(float(view.origin.x() + horizontal * view.scale),
                             float(view.origin.y() - (v.z() * cp - forward * sp) * view.scale),
                             float((v.z() * sp + forward * cp) * view.scale));
        };
    }
};

void drawMeshFaces(QPainter &p, const std::array<QVector<SurfaceTriangle>, 2> &meshes,
                   const std::array<QColor, 2> &colors, const MeshProjection &view,
                   const QVector<double> &depths = {}, const QVector<double> &sensitivity = {},
                   double threshold = 0.) {
    const auto project = view.projector();
    struct Face {
        std::array<QVector3D, 3> projected;
        double depth;
        int layer;
        bool lowResponse = false;
    };
    std::vector<Face> faces;
    faces.reserve(meshes[0].size() + meshes[1].size());
    for (int layer = 0; layer < 2; ++layer)
        for (const auto &triangle : meshes[layer]) {
            Face f{{project(triangle[0]), project(triangle[1]), project(triangle[2])}, 0, layer};
            f.depth = f.projected[0].z() + f.projected[1].z() + f.projected[2].z();
            if (threshold > 0 && depths.size() >= 2 && sensitivity.size() == depths.size()) {
                const double x = (triangle[0].x() + triangle[1].x() + triangle[2].x()) / 3.;
                const auto index = std::clamp(
                    qRound((x - depths.first()) / (depths.last() - depths.first()) * (depths.size() - 1)), 0,
                    int(depths.size() - 1));
                f.lowResponse = sensitivity[index] < threshold;
            }
            faces.push_back(f);
        }
    std::sort(faces.begin(), faces.end(), [](const Face &a, const Face &b) { return a.depth < b.depth; });
    p.setRenderHint(QPainter::Antialiasing, false);
    p.setPen(Qt::NoPen);
    const std::array<QColor, 2> dark{colors[0].darker(160), colors[1].darker(160)};
    const std::array<QColor, 2> light{colors[0].lighter(155), colors[1].lighter(155)};
    for (const auto &face : faces) {
        const auto &q = face.projected;
        const auto normal = QVector3D::normal(q[1] - q[0], q[2] - q[0]);
        auto color = blend(dark[face.layer], light[face.layer], std::abs(normal.z()));
        if (face.lowResponse)
            color.setAlpha(75);
        p.setBrush(color);
        QPointF points[]{q[0].toPointF(), q[1].toPointF(), q[2].toPointF()};
        p.drawPolygon(points, 3);
    }
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
    m_depthProjector = {};
    m_axes = axes;
    m_mesh[0] = first;
    m_mesh[1] = second;
    m_meshColors = colors;
    m_caption = caption;
    update();
}

bool BondFigures::hasMesh() const {
    return m_kind == Kind::Mesh && (!m_mesh[0].isEmpty() || !m_mesh[1].isEmpty());
}

void BondFigures::setDepthProfile(const QVector<double> &profile, int index, double threshold) {
    m_depthProfile = profile;
    m_depthIndex = index;
    m_lowThreshold = threshold;
    update();
}

QRectF BondFigures::depthPlotRect() const {
    return QRectF(85, height() - 87, width() - 142, 47);
}

std::function<QImage(double)> BondFigures::rotationFrames(const QSize &size) const {
    if (!hasMesh() || size.isEmpty())
        return {};
    const double yaw = m_yaw, pitch = m_pitch;
    return [meshes = m_mesh, colors = m_meshColors, size, yaw, pitch,
            scale = 1., center = QVector3D(), palette = QList<QRgb>()](double turn) mutable {
        if (palette.isEmpty()) {
            // Fit the bounding cylinder at the chosen elevation for the whole turn.
            QVector3D lo = (meshes[0].isEmpty() ? meshes[1] : meshes[0]).first()[0], hi = lo;
            for (const auto &mesh : meshes)
                for (const auto &triangle : mesh)
                    for (const auto &vertex : triangle)
                        for (int axis = 0; axis < 3; ++axis) {
                            lo[axis] = std::min(lo[axis], vertex[axis]);
                            hi[axis] = std::max(hi[axis], vertex[axis]);
                        }
            center = (lo + hi) / 2;
            const auto span = hi - lo;
            const double diameter = std::hypot(double(span.x()), double(span.y()));
            const double height = std::abs(std::sin(pitch)) * diameter + std::abs(std::cos(pitch)) * span.z();
            const double fit = std::max(diameter / size.width(), height / size.height());
            if (fit > 0.)
                scale = .9 / fit;

            palette.append(qRgb(255, 255, 255));
            const int shades = !meshes[0].isEmpty() && !meshes[1].isEmpty() ? 127 : 255;
            for (int layer = 0; layer < 2; ++layer) {
                if (meshes[layer].isEmpty())
                    continue;
                const auto dark = colors[layer].darker(160), light = colors[layer].lighter(155);
                for (int i = 0; i < shades; ++i)
                    palette.append(blend(dark, light, double(i) / (shades - 1)).rgb());
            }
            while (palette.size() < 256)
                palette.append(qRgb(255, 255, 255));
        }
        constexpr double fullTurn = 6.283185307179586;
        QImage image(size, QImage::Format_RGB32);
        image.fill(Qt::white);
        QPainter painter(&image);
        drawMeshFaces(painter, meshes, colors,
                      {center, QPointF(size.width() * .5, size.height() * .5), scale, yaw + fullTurn * turn, pitch});
        painter.end();
        return image.convertToFormat(QImage::Format_Indexed8, palette, Qt::ThresholdDither | Qt::AvoidDither);
    };
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
    const bool inspect = m_depthProfile.size() == x.size() && m_depthIndex >= 0 && m_depthIndex < x.size();
    const double modelHeight = height() - (inspect ? 94 : 0);
    const double scale = m_zoom * std::min((width() - 160.) / projectedW, (modelHeight - 110.) / projectedH);
    const MeshProjection view{center, QPointF(width() * .5, modelHeight * .49), scale, m_yaw, m_pitch};
    const auto project = view.projector();
    std::array<QVector3D, 8> box;
    for (int i = 0; i < 8; ++i)
        box[i] =
            project(QVector3D(i & 1 ? hi.x() : lo.x(), i & 2 ? hi.y() : lo.y(), i & 4 ? hi.z() : lo.z()));
    p.setPen(QPen(QColor("#aeb6be"), 1, Qt::DashLine));
    for (int i = 0; i < 8; ++i)
        for (int bit : {1, 2, 4})
            if (!(i & bit))
                p.drawLine(box[i].toPointF(), box[i | bit].toPointF());
    drawMeshFaces(p, m_mesh, m_meshColors, view, x, m_depthProfile, inspect ? m_lowThreshold : 0.);
    p.setRenderHint(QPainter::Antialiasing);
    if (inspect) {
        auto plane = [&](double depth) {
            QPolygonF face;
            for (int corner : {0, 1, 3, 2})
                face << project(QVector3D(depth, corner & 1 ? hi.y() : lo.y(), corner & 2 ? hi.z() : lo.z()))
                            .toPointF();
            return face;
        };
        // Highlight complete low-response depth intervals, including empty parts of the field.
        p.setPen(Qt::NoPen);
        p.setBrush(QColor(214, 140, 39, 34));
        for (int i = 0; i < x.size();) {
            if (m_depthProfile[i] >= m_lowThreshold) {
                ++i;
                continue;
            }
            const int first = i++;
            while (i < x.size() && m_depthProfile[i] < m_lowThreshold)
                ++i;
            const auto a = plane(first ? (x[first - 1] + x[first]) * .5 : x.first());
            const auto b = plane(i < x.size() ? (x[i - 1] + x[i]) * .5 : x.last());
            p.drawPolygon(a);
            p.drawPolygon(b);
            for (int k = 0; k < 4; ++k) {
                const int n = (k + 1) % 4;
                p.drawPolygon(QPolygonF{a[k], a[n], b[n], b[k]});
            }
        }
        p.setPen(QPen(QColor("#176e85"), 1));
        p.setBrush(QColor(40, 157, 177, 27));
        p.drawPolygon(plane(x[m_depthIndex]));
        m_xAxisOnScreen = QLineF(box[0].toPointF(), box[1].toPointF());
        m_depthProjector = project;
    }
    p.setPen(Qt::black);
    const QVector<QVector<double>> axes{x, y, z};
    for (int i = 0; i < 3; ++i) {
        const auto start = box[0].toPointF(), end = box[std::size_t(1) << i].toPointF();
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
    p.drawText(QRectF(8, modelHeight - 30, width() - 16, 24), Qt::AlignCenter, m_caption);
    if (!inspect)
        return;
    const auto area = depthPlotRect();
    auto plotted = [&](int i) {
        return QPointF(area.left() + area.width() * (x[i] - x.first()) / (x.last() - x.first()),
                       area.bottom() - area.height() * m_depthProfile[i]);
    };
    p.fillRect(area, Qt::white);
    for (int i = 0; i < x.size(); ++i) {
        if (m_depthProfile[i] < m_lowThreshold) {
            const double left = i ? (plotted(i - 1).x() + plotted(i).x()) * .5 : area.left();
            const double right = i + 1 < x.size() ? (plotted(i).x() + plotted(i + 1).x()) * .5 : area.right();
            p.fillRect(QRectF(left, area.top(), right - left, area.height()), QColor(214, 140, 39, 42));
        }
    }
    p.setBrush(Qt::NoBrush);
    p.setPen(QColor("#89939a"));
    p.drawRect(area);
    p.drawText(QRectF(0, area.top() - 8, area.left() - 6, 18), Qt::AlignRight, "100%");
    p.drawText(QRectF(0, area.bottom() - 9, area.left() - 6, 18), Qt::AlignRight, "0%");
    QPainterPath curve(plotted(0));
    for (int i = 1; i < x.size(); ++i)
        curve.lineTo(plotted(i));
    p.setPen(QPen(QColor("#176e85"), 1.6));
    p.drawPath(curve);
    const double selected = plotted(m_depthIndex).x();
    p.drawLine(QPointF(selected, area.top()), QPointF(selected, area.bottom()));
    p.setBrush(QColor("#176e85"));
    p.drawEllipse(plotted(m_depthIndex), 3., 3.);
    p.setPen(Qt::black);
    p.drawText(QPointF(area.left() - 6, area.bottom() + 16), QString::number(x.first(), 'g', 5));
    p.drawText(QRectF(area.right() - 72, area.bottom() + 2, 80, 20), Qt::AlignRight,
               QString::number(x.last(), 'g', 5));
    p.drawText(QRectF(area.left(), area.bottom() + 1, area.width(), 23), Qt::AlignCenter,
               tr("相对灵敏度 · X (mm) · 点击曲线选取深度"));
}

void BondFigures::mousePressEvent(QMouseEvent *event) {
    if (m_kind == Kind::Montage || event->button() != Qt::LeftButton) {
        QWidget::mousePressEvent(event);
        return;
    }
    if (m_kind == Kind::Mesh && m_axes.size() == 3 && m_depthIndex >= 0 &&
        m_depthProfile.size() == m_axes[0].size()) {
        const auto pos = event->position();
        const auto graph = depthPlotRect();
        double t = -1.;
        if (graph.adjusted(-3, -3, 3, 3).contains(pos))
            t = std::clamp((pos.x() - graph.left()) / graph.width(), 0., 1.);
        else if (m_xAxisOnScreen.length() > 12) {
            const QPointF d = m_xAxisOnScreen.p2() - m_xAxisOnScreen.p1();
            const double fraction =
                QPointF::dotProduct(pos - m_xAxisOnScreen.p1(), d) / QPointF::dotProduct(d, d);
            if (fraction >= 0 && fraction <= 1 && QLineF(pos, m_xAxisOnScreen.pointAt(fraction)).length() < 6)
                t = fraction;
        }
        if (t >= 0) {
            emit depthSelected(qRound(t * (m_axes[0].size() - 1)));
            event->accept();
            return;
        }
    }
    m_lastMouse = event->position();
    m_pressPosition = m_lastMouse;
    m_dragged = false;
    setCursor(Qt::ClosedHandCursor);
    event->accept();
}

void BondFigures::mouseMoveEvent(QMouseEvent *event) {
    if (m_kind != Kind::Montage && cursor().shape() == Qt::ClosedHandCursor && (event->buttons() & Qt::LeftButton)) {
        if (!m_dragged && QLineF(m_pressPosition, event->position()).length() < 3.)
            return;
        m_dragged = true;
        const auto delta = event->position() - m_lastMouse;
        m_yaw += delta.x() * .012;
        m_pitch = std::clamp(m_pitch + delta.y() * .012, -1.5, 1.5);
        m_lastMouse = event->position();
        update();
        event->accept();
    }
}

void BondFigures::mouseReleaseEvent(QMouseEvent *event) {
    if (event->button() == Qt::LeftButton && cursor().shape() == Qt::ClosedHandCursor && !m_dragged &&
        m_kind == Kind::Mesh && m_depthProjector && m_depthIndex >= 0 && m_axes.size() == 3) {
        // A screen point may hit several faces. Take the frontmost hit and interpolate its physical X.
        auto cross = [](QPointF a, QPointF b) { return a.x() * b.y() - a.y() * b.x(); };
        double front = -INFINITY, selectedX = 0.;
        for (const auto &mesh : m_mesh)
            for (const auto &face : mesh) {
                const auto a = m_depthProjector(face[0]), b = m_depthProjector(face[1]),
                           c = m_depthProjector(face[2]);
                const auto u = (b - a).toPointF(), v = (c - a).toPointF(),
                           d = event->position() - a.toPointF();
                const double det = cross(u, v);
                if (std::abs(det) < 1e-8)
                    continue;
                const double s = cross(d, v) / det, t = cross(u, d) / det;
                const double depth = a.z() + s * (b.z() - a.z()) + t * (c.z() - a.z());
                if (s >= 0 && t >= 0 && s + t <= 1 && depth > front) {
                    front = depth;
                    selectedX =
                        face[0].x() + s * (face[1].x() - face[0].x()) + t * (face[2].x() - face[0].x());
                }
            }
        if (std::isfinite(front)) {
            const auto &x = m_axes[0];
            emit depthSelected(
                std::clamp(qRound((selectedX - x.first()) / (x.last() - x.first()) * (x.size() - 1)), 0,
                           int(x.size() - 1)));
        }
    }
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
