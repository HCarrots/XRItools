#pragma once
#include <QVector>
#include <QVector3D>
#include <array>
#include <atomic>

using SurfaceTriangle = std::array<QVector3D, 3>;

namespace Isosurface {
// Samples are laid out as [X,Y,Z], with Z contiguous.
QVector<SurfaceTriangle> extract(const QVector<QVector<double>> &axes, const double *field,
                                 qsizetype count, double level, const std::atomic_bool &stop);
}
