#pragma once
#include <QColor>
#include <QFont>
#include <QPainter>
#include <QStringList>
#include <QVector>
#include <algorithm>
#include <cmath>

// Shared scientific axes. Values stay in measurement units; only positions are
// transformed for signed-log displays. No changes to the underlying samples.
namespace PlotStyle {
struct Ticks {
    QVector<double> major, minor;
    double step = 1;
    bool logarithmic = false;
    bool decades = false;
};
inline double transform(double v, bool logarithmic) {
    return logarithmic ? std::copysign(std::log1p(std::abs(v)), v) : v;
}
inline Ticks ticks(double lo, double hi, int count = 5, bool logarithmic = false) {
    Ticks t;
    t.logarithmic = logarithmic;
    if (!std::isfinite(lo) || !std::isfinite(hi) || !(hi > lo))
        return t;
    const double span = hi - lo;
    if (!std::isfinite(span))
        return t;
    const double rough = span / std::clamp(count, 2, 20);
    const double base = std::pow(10., std::floor(std::log10(rough)));
    if (!(base > 0) || !std::isfinite(base))
        return t;
    const double u = rough / base;
    t.step = (u < 1.5 ? 1 : u < 2.25 ? 2 : u < 3.75 ? 2.5 : u < 7.5 ? 5 : 10) * base;
    if (logarithmic && transform(hi, true) - transform(lo, true) > 2) {
        if (lo <= 0 && hi >= 0)
            t.major.append(0);
        double largest = std::max(std::abs(lo), std::abs(hi));
        int last = std::min(307, int(std::ceil(std::log10(largest))));
        double nearest = lo > 0 ? lo : hi < 0 ? -hi : 1;
        int first = std::max(-307, std::max(last - 12, int(std::floor(std::log10(nearest)))));
        for (int e = first; e <= last; ++e)
            for (int sign : {-1, 1})
                for (int m : {1, 2, 5}) {
                    double v = sign * m * std::pow(10., e);
                    if (v >= lo && v <= hi) {
                        if (m == 1)
                            t.major.append(v);
                        else
                            t.minor.append(v);
                    }
                }
        std::sort(t.major.begin(), t.major.end());
        if (t.major.size() >= 2) {
            t.decades = true;
            return t;
        }
        t.major.clear();
        t.minor.clear();
    }
    double first = std::ceil(lo / t.step - 1e-10) * t.step;
    for (int i = 0; i < 42; ++i) {
        double v = first + i * t.step;
        if (v > hi + t.step * 1e-9)
            break;
        if (v >= lo - t.step * 1e-9)
            t.major.append(std::abs(v) < t.step * 1e-9 ? 0 : v);
    }
    const double small = t.step / 5;
    first = std::ceil(lo / small) * small;
    for (int i = 0; i < 210; ++i) {
        double v = first + i * small;
        if (v > hi + small * 1e-9)
            break;
        if (std::abs(v / t.step - std::round(v / t.step)) > 1e-7)
            t.minor.append(v);
    }
    return t;
}
inline QString label(double value, const Ticks &t) {
    if (!std::isfinite(value))
        return QStringLiteral("—");
    if (value == 0)
        return QStringLiteral("0");
    double largest = 0;
    for (double v : t.major)
        largest = std::max(largest, std::abs(v));
    if (largest >= 1e5 || (largest > 0 && largest < .001)) {
        int precision =
            t.decades ? 0 : std::clamp(int(std::ceil(std::log10(largest) - std::log10(t.step))), 1, 14);
        return QString::number(value, 'e', precision);
    }
    if (t.decades)
        return QString::number(value, 'g', 4);
    int decimals =
        std::clamp(int(std::ceil(-std::log10(t.step))) +
                       (std::abs(t.step / std::pow(10., std::floor(std::log10(t.step))) - 2.5) < .01 ? 1 : 0),
                   0, 12);
    return QString::number(value, 'f', decimals);
}
inline QColor curveColor(int index) {
    static const char *colors[]{"#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9", "#333333"};
    return QColor(colors[index % 7]);
}
inline void axis(QPainter &p, const QRectF &r, const Ticks &t, double lo, double hi, bool horizontal,
                 bool down, bool grid, bool minor = true, int fontPoints = 9, bool fineGrid = false) {
    const double a = transform(lo, t.logarithmic), b = transform(hi, t.logarithmic);
    if (!(b > a))
        return;
    auto position = [&](double v) {
        double f = (transform(v, t.logarithmic) - a) / (b - a);
        return horizontal ? r.left() + f * r.width()
               : down     ? r.top() + f * r.height()
                          : r.bottom() - f * r.height();
    };
    p.save();
    p.setFont(QFont(QStringLiteral("Arial"), fontPoints));
    const QFontMetricsF metrics(p.font());
    for (int major = 0; major <= 1; ++major) {
        if (!major && !minor)
            continue;
        double previousLabel = -INFINITY;
        for (double v : major ? t.major : t.minor) {
            double pos = position(v);
            if (!std::isfinite(pos))
                continue;
            if (grid && (major || fineGrid)) {
                p.setPen(fineGrid ? QPen(major ? QColor(192, 192, 192) : QColor(229, 229, 229), .6)
                                  : QPen(QColor(210, 216, 221), .6, Qt::DashLine));
                if (horizontal)
                    p.drawLine(QPointF(pos, r.top()), QPointF(pos, r.bottom()));
                else
                    p.drawLine(QPointF(r.left(), pos), QPointF(r.right(), pos));
            }
            p.setPen(QPen(QColor("#242424"), .8));
            double length = major ? 5 : 2.5;
            if (horizontal) {
                p.drawLine(QPointF(pos, r.bottom()), QPointF(pos, r.bottom() - length));
                p.drawLine(QPointF(pos, r.top()), QPointF(pos, r.top() + length));
                if (major) {
                    const auto text = label(v, t);
                    const double spacing =
                        fineGrid ? std::max(52., metrics.horizontalAdvance(text) + 10) : 80.;
                    if (std::abs(pos - previousLabel) >= spacing) {
                        p.drawText(QRectF(pos - 48, r.bottom() + 7, 96, 20), Qt::AlignHCenter, text);
                        previousLabel = pos;
                    }
                }
            } else {
                p.drawLine(QPointF(r.left(), pos), QPointF(r.left() + length, pos));
                p.drawLine(QPointF(r.right(), pos), QPointF(r.right() - length, pos));
                if (major && std::abs(pos - previousLabel) >= 24) {
                    p.drawText(QRectF(r.left() - 74, pos - 10, 65, 20), Qt::AlignRight | Qt::AlignVCenter,
                               label(v, t));
                    previousLabel = pos;
                }
            }
        }
    }
    p.restore();
}
} // namespace PlotStyle
