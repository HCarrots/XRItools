#include "Isosurface.h"
#include "ImagingIO.h"
#include <algorithm>
#include <cmath>

using namespace ImageUtil;

QVector<SurfaceTriangle> Isosurface::extract(const QVector<QVector<double>> &axes, const double *field,
                                             qsizetype count, double level, const std::atomic_bool &stop) {
    cancel(stop);
    if (axes.size() != 3 || !field || !std::isfinite(level))
        fail(QStringLiteral("等值面输入无效。"));
    const qsizetype nx = axes[0].size(), ny = axes[1].size(), nz = axes[2].size();
    if (nx < 2 || ny < 2 || nz < 2 || nx > count / ny / nz || nx * ny * nz != count)
        fail(QStringLiteral("等值面需要与体数据一致的三维坐标轴。"));
    for (const auto &axis : axes)
        for (qsizetype i = 0; i < axis.size(); ++i)
            if (!std::isfinite(axis[i]) || (i > 0 && axis[i] <= axis[i - 1]))
                fail(QStringLiteral("等值面的空间坐标须有限且严格递增。"));
    constexpr int tetrahedra[6][4]{{0, 1, 3, 7}, {0, 3, 2, 7}, {0, 2, 6, 7},
                                   {0, 6, 4, 7}, {0, 4, 5, 7}, {0, 5, 1, 7}};
    QVector<SurfaceTriangle> mesh;
    for (qsizetype x = 0; x + 1 < nx; ++x) {
        cancel(stop);
        for (qsizetype y = 0; y + 1 < ny; ++y)
            for (qsizetype z = 0; z + 1 < nz; ++z) {
                if (z % 1024 == 0)
                    cancel(stop);
                std::array<double, 8> values;
                for (int i = 0; i < 8; ++i) {
                    const qsizetype ix = x + (i & 1), iy = y + ((i >> 1) & 1), iz = z + (i >> 2);
                    values[i] = field[(ix * ny + iy) * nz + iz];
                }
                const auto [min, max] = std::minmax_element(values.cbegin(), values.cend());
                if (*min >= level || *max < level ||
                    !std::all_of(values.cbegin(), values.cend(), [](double f) { return std::isfinite(f); }))
                    continue;
                std::array<QVector3D, 8> points;
                for (int i = 0; i < 8; ++i)
                    points[i] = QVector3D(float(axes[0][x + (i & 1)]), float(axes[1][y + ((i >> 1) & 1)]),
                                          float(axes[2][z + (i >> 2)]));
                auto edge = [&](int a, int b) {
                    return points[a] +
                           float((level - values[a]) / (values[b] - values[a])) * (points[b] - points[a]);
                };
                for (const auto &tetra : tetrahedra) {
                    int above[4], below[4], a = 0, b = 0;
                    for (int i : tetra)
                        if (values[i] >= level)
                            above[a++] = i;
                        else
                            below[b++] = i;
                    if (a == 1)
                        mesh.append(SurfaceTriangle{edge(above[0], below[0]), edge(above[0], below[1]),
                                                 edge(above[0], below[2])});
                    else if (a == 3)
                        mesh.append(SurfaceTriangle{edge(below[0], above[0]), edge(below[0], above[1]),
                                                 edge(below[0], above[2])});
                    else if (a == 2) {
                        const auto p = edge(above[0], below[0]), q = edge(above[0], below[1]),
                                   r = edge(above[1], below[0]), s = edge(above[1], below[1]);
                        mesh.append(SurfaceTriangle{p, q, r});
                        mesh.append(SurfaceTriangle{q, s, r});
                    }
                }
                if (mesh.size() > 500000)
                    fail(QStringLiteral("等值面过于复杂，请调整阈值或使用较粗的 X 网格。"));
            }
    }
    return mesh;
}
