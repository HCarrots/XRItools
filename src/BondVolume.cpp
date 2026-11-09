#include "BondVolume.h"
#include "ImagingIO.h"
#include <QDir>
#include <QFileInfo>
#include <QJsonDocument>
#include <QTemporaryDir>
#include <algorithm>
#include <cmath>

using namespace ImageUtil;
namespace {
constexpr qsizetype maxSamples = 16000000;
const QString fieldPath = QStringLiteral("/entry/bond/intensity");
const QStringList axisPaths{"/entry/bond/incident_energy_kev", "/entry/bond/x_mm", "/entry/bond/y_mm",
                            "/entry/bond/z_mm"};

qsizetype checkedSize(const QVector<quint64> &shape) {
    qsizetype count = 1;
    if (shape.size() != 4)
        fail(QStringLiteral("结果必须为 [E,X,Y,Z] 四维数组。"));
    for (auto n : shape) {
        if (n < 2 || n > quint64(maxSamples / count))
            fail(QStringLiteral("四维结果尺寸无效或超过 1600 万采样点。"));
        count *= qsizetype(n);
    }
    return count;
}

void checkAxis(const QVector<double> &axis) {
    if (axis.size() < 2)
        fail(QStringLiteral("各轴至少需要两个坐标。"));
    for (qsizetype i = 0; i < axis.size(); ++i)
        if (!std::isfinite(axis[i]) || (i && axis[i] <= axis[i - 1]))
            fail(QStringLiteral("四维坐标轴必须有限且严格递增。"));
}

struct Bracket {
    int first, last;
    double fraction;
};
Bracket bracket(const QVector<double> &axis, double at) {
    constexpr double epsilon = 1e-9;
    if (!std::isfinite(at) || at < axis.first() - epsilon || at > axis.last() + epsilon)
        fail(QStringLiteral("指定能量超出实测范围。"));
    int hi = int(std::lower_bound(axis.cbegin(), axis.cend(), at) - axis.cbegin());
    if (hi == axis.size())
        hi--;
    if (hi == 0 || std::abs(axis[hi] - at) <= epsilon)
        return {hi, hi, 0};
    if (std::abs(axis[hi - 1] - at) <= epsilon)
        return {hi - 1, hi - 1, 0};
    return {hi - 1, hi, (at - axis[hi - 1]) / (axis[hi] - axis[hi - 1])};
}

double valueAt(const double *data, qsizetype stride, Bracket b) {
    const double first = data[b.first * stride];
    return b.first == b.last ? first : first + b.fraction * (data[b.last * stride] - first);
}

struct SpectralIntegral {
    QVector<double> energy;
    QVector<Bracket> points;

    SpectralIntegral(const QVector<double> &axis, QPair<double, double> range,
                     const QVector<bool> &usable) {
        const auto [lo, hi] = range;
        if (!(lo < hi))
            fail(QStringLiteral("积分区间必须递增。"));
        energy.append(lo);
        for (double e : axis)
            if (lo < e && e < hi)
                energy.append(e);
        energy.append(hi);
        for (double e : energy) {
            const auto b = bracket(axis, e);
            if (!usable[b.first] || !usable[b.last])
                fail(QStringLiteral("积分区间包含未收敛通道或超出参考轮廓覆盖。"));
            points.append(b);
        }
    }

    double area(const QVector<double> &row) const {
        double sum = 0, prev = valueAt(row.constData(), 1, points.first());
        for (qsizetype i = 1; i < points.size(); ++i) {
            const double current = valueAt(row.constData(), 1, points[i]);
            sum += .5 * (prev + current) * (energy[i] - energy[i - 1]);
            prev = current;
        }
        return sum;
    }
};

bool supported(const BondVolume &v, qsizetype voxel) {
    const auto yz = v.axes[2].size() * v.axes[3].size();
    return v.sensitivity[voxel / yz] > v.diagnostics.value("x_sensitivity_floor").toDouble();
}
} // namespace

qsizetype BondVolume::voxelCount() const {
    return axes.size() == 4 ? axes[1].size() * axes[2].size() * axes[3].size() : 0;
}

qsizetype BondVolume::voxel(int x, int y, int z) const {
    return (qsizetype(x) * axes[2].size() + y) * axes[3].size() + z;
}

QString BondImaging::diagnosticPngPath(const QString &path) {
    const QFileInfo result(path);
    return result.dir().filePath(result.completeBaseName() + ".diagnostics.png");
}

void BondImaging::save(const QString &path, const BondVolume &volume, const std::atomic_bool &stop) {
    cancel(stop);
    QVector<quint64> shape;
    for (const auto &axis : volume.axes) {
        checkAxis(axis);
        shape.append(quint64(axis.size()));
    }
    if (checkedSize(shape) != volume.intensity.size() || volume.sensitivity.size() != qsizetype(shape[1]) ||
        volume.converged.size() != qsizetype(shape[0]))
        fail(QStringLiteral("四维结果与坐标或诊断尺寸不一致。"));
    const QFileInfo target(path);
    if (target.exists())
        fail(QStringLiteral("结果文件已存在，请选择新文件名。"));
    const bool plots = !volume.diagnosticImage.isNull();
    const auto png = diagnosticPngPath(target.absoluteFilePath());
    if (plots && QFileInfo::exists(png))
        fail(QStringLiteral("诊断 PNG 已存在，请选择新的结果文件名。"));
    QTemporaryDir temp(QDir(target.absolutePath()).filePath(".xri-bond-XXXXXX"));
    if (!temp.isValid())
        fail(QStringLiteral("无法在输出目录建立临时文件。"));
    const auto staged = temp.filePath("volume.h5");
    const auto stagedPng = temp.filePath("diagnostics.png");
    if (plots && !volume.diagnosticImage.save(stagedPng))
        fail(QStringLiteral("无法写入诊断 PNG。"));
    {
        ImagingH5 hdf(staged, true);
        auto axis = [&](const QString &name, const QVector<double> &data) {
            hdf.write(name, {quint64(data.size())}, data.constData(), false);
        };
        hdf.write(fieldPath, shape, volume.intensity.constData());
        for (int i = 0; i < axisPaths.size(); ++i)
            axis(axisPaths[i], volume.axes[i]);
        hdf.scales(fieldPath, axisPaths);
        axis("/entry/bond/x_sensitivity", volume.sensitivity);
        QVector<quint64> valid;
        for (bool good : volume.converged)
            valid.append(good ? 1 : 0);
        hdf.writeUInt("/entry/bond/channel_converged", {quint64(valid.size())}, valid.constData(), false);
        hdf.attributes(fieldPath, {{"axis_order", "E,X,Y,Z"},
                                   {"units", "relative intensity"},
                                   {"energy_kind", "incident"},
                                   {"schema", "XRI.bond.v1"}});
        hdf.writeText("/entry/bond/diagnostics_json",
                      QString::fromUtf8(QJsonDocument(volume.diagnostics).toJson(QJsonDocument::Compact)));
        for (const auto &item : volume.bandDiagnostics) {
            cancel(stop);
            const auto base = "/entry/bond/band_diagnostics/" + item.crystal;
            const auto rows = quint64(item.forward.rows());
            hdf.write(base + "/operator_px", {rows, shape[1]}, item.forward.data());
            hdf.write(base + "/predicted_pyz", {rows, shape[2], shape[3]}, item.predicted.data());
            hdf.write(base + "/residual_pyz", {rows, shape[2], shape[3]}, item.residual.data());
        }
    }
    cancel(stop);
    if (plots && !QFile::rename(stagedPng, png))
        fail(QStringLiteral("无法提交诊断 PNG。"));
    if (!QFile::rename(staged, target.absoluteFilePath())) {
        if (plots)
            QFile::remove(png);
        fail(QStringLiteral("无法提交四维结果文件。"));
    }
}

BondVolume BondImaging::load(const QString &path) {
    const auto shape = HdfFile::inspect(path, fieldPath).dims;
    checkedSize(shape);
    ImagingH5 hdf(path);
    BondVolume out;
    auto readVector = [&](const QString &name, quint64 expected) {
        if (HdfFile::inspect(path, name).dims != QVector<quint64>{expected})
            fail(QStringLiteral("四维结果坐标或诊断尺寸不一致。"));
        return hdf.read(name);
    };
    for (int i = 0; i < axisPaths.size(); ++i) {
        auto axis = readVector(axisPaths[i], shape[i]);
        checkAxis(axis);
        out.axes.append(std::move(axis));
    }
    out.sensitivity = readVector("/entry/bond/x_sensitivity", shape[1]);
    for (double value : out.sensitivity)
        if (!std::isfinite(value) || value < 0)
            fail(QStringLiteral("无效的 X 灵敏度。"));
    for (double value : readVector("/entry/bond/channel_converged", shape[0])) {
        if (value != 0 && value != 1)
            fail(QStringLiteral("无效的通道收敛标记。"));
        out.converged.append(value == 1);
    }
    out.diagnostics = QJsonDocument::fromJson(hdf.text("/entry/bond/diagnostics_json").toUtf8()).object();
    if (out.diagnostics.value("schema").toString() != "XRI.bond.v1")
        fail(QStringLiteral("无法识别四维结果版本。"));
    out.intensity = hdf.read(fieldPath);
    return out;
}

QVector<double> BondImaging::process(const BondVolume &v, const QVector<double> &energy,
                                     const BondSpectralSettings &settings, const std::atomic_bool &stop) {
    checkAxis(energy);
    if (energy.size() != v.axes[0].size())
        fail(QStringLiteral("处理能量轴与四维结果不一致。"));
    auto checkWindow = [&](QPair<double, double> range) {
        if (!(range.first < range.second))
            fail(QStringLiteral("拟合或归一化能区必须递增。"));
        bracket(energy, range.first);
        bracket(energy, range.second);
    };
    for (auto range : settings.background)
        checkWindow(range);
    const auto *compton = settings.compton ? &*settings.compton : nullptr;
    auto usable = v.converged;
    std::optional<SpectralIntegral> match, normalize;
    double referenceArea = 0;
    if (compton) {
        if (compton->reference.size() != energy.size() || compton->valence.size() != energy.size() ||
            compton->correction.size() != energy.size())
            fail(QStringLiteral("康普顿参考与重建能量轴不一致。"));
        for (qsizetype e = 0; e < energy.size(); ++e) {
            if (!(std::isfinite(compton->correction[e]) && compton->correction[e] > 0))
                fail(QStringLiteral("康普顿截面修正无效。"));
            usable[e] = usable[e] && std::isfinite(compton->valence[e]);
        }
        match.emplace(energy, compton->matchWindow, usable);
        referenceArea = match->area(compton->reference);
        if (!(std::isfinite(referenceArea) && referenceArea > 0))
            fail(QStringLiteral("幅度匹配区间需要完整覆盖和正的宽谱积分。"));
    }
    if (settings.normalization)
        normalize.emplace(energy, *settings.normalization, usable);
    QVector<int> fitPoints;
    for (int e = 0; e < energy.size(); ++e)
        if (usable[e] &&
            std::any_of(settings.background.cbegin(), settings.background.cend(),
                        [&](auto w) { return energy[e] >= w.first - 1e-9 && energy[e] <= w.second + 1e-9; }))
            fitPoints.append(e);
    if (!settings.background.isEmpty() && fitPoints.size() < 2)
        fail(QStringLiteral("背底能区至少需要两个已收敛且有参考覆盖的能量点。"));
    const auto n = v.voxelCount();
    QVector<double> out(v.intensity.size(), NAN), row(energy.size());
    const auto *input = v.intensity.constData();
    auto *result = out.data();
    for (qsizetype i = 0; i < n; ++i) {
        if (i % 128 == 0)
            cancel(stop);
        if (!supported(v, i))
            continue;
        for (int e = 0; e < energy.size(); ++e)
            row[e] = usable[e] ? input[e * n + i] : NAN;
        if (compton) {
            for (qsizetype e = 0; e < energy.size(); ++e)
                row[e] *= compton->correction[e];
            const double scale = match->area(row) / referenceArea;
            if (!(std::isfinite(scale) && scale > 0))
                continue;
            for (qsizetype e = 0; e < energy.size(); ++e)
                row[e] -= scale * compton->valence[e];
        }
        if (!fitPoints.isEmpty()) {
            double meanX = 0, meanY = 0;
            int count = 0;
            for (int e : fitPoints)
                if (std::isfinite(row[e])) {
                    meanX += energy[e];
                    meanY += row[e];
                    ++count;
                }
            if (count < 2)
                continue;
            meanX /= count;
            meanY /= count;
            double xx = 0, xy = 0;
            for (int e : fitPoints)
                if (std::isfinite(row[e])) {
                    const double dx = energy[e] - meanX;
                    xx += dx * dx;
                    xy += dx * (row[e] - meanY);
                }
            for (int e = 0; e < energy.size(); ++e)
                row[e] -= meanY + (xy / xx) * (energy[e] - meanX);
        }
        // Return to the reconstruction's intensity scale before any optional area normalization.
        if (compton)
            for (qsizetype e = 0; e < energy.size(); ++e)
                row[e] /= compton->correction[e];
        const double area = normalize ? normalize->area(row) : 1.;
        if (!(std::isfinite(area) && area > 0))
            continue;
        for (int e = 0; e < energy.size(); ++e)
            result[e * n + i] = row[e] / area;
    }
    cancel(stop);
    if (std::none_of(out.cbegin(), out.cend(), [](double value) { return std::isfinite(value); }))
        fail(QStringLiteral("所选参数未生成有效体素谱，请检查能区、参考覆盖和积分面积。"));
    return out;
}

QVector<double> BondImaging::energySlice(const BondVolume &v, const QVector<double> &data,
                                         const QVector<double> &energy, double at, bool requireConvergence) {
    checkAxis(energy);
    if (energy.size() != v.axes[0].size() || data.size() != v.intensity.size())
        fail(QStringLiteral("能量轴与四维结果不一致。"));
    const auto b = bracket(energy, at);
    if (requireConvergence && (!v.converged[b.first] || !v.converged[b.last]))
        fail(QStringLiteral("指定能量所需的通道未收敛，请选择已收敛通道。"));
    const auto n = v.voxelCount();
    QVector<double> out(n, NAN);
    for (qsizetype i = 0; i < n; ++i)
        if (!requireConvergence || supported(v, i))
            out[i] = valueAt(data.constData() + i, n, b);
    return out;
}

QVector<BondTriangle> BondImaging::surface(const BondVolume &v, const QVector<double> &field, double level,
                                           const std::atomic_bool &stop) {
    if (field.size() != v.voxelCount() || !std::isfinite(level))
        fail(QStringLiteral("等值面输入无效。"));
    constexpr int tetrahedra[6][4]{{0, 1, 3, 7}, {0, 3, 2, 7}, {0, 2, 6, 7},
                                   {0, 6, 4, 7}, {0, 4, 5, 7}, {0, 5, 1, 7}};
    QVector<BondTriangle> mesh;
    for (int x = 0; x + 1 < v.axes[1].size(); ++x) {
        cancel(stop);
        for (int y = 0; y + 1 < v.axes[2].size(); ++y)
            for (int z = 0; z + 1 < v.axes[3].size(); ++z) {
                std::array<double, 8> values;
                std::array<QVector3D, 8> points;
                for (int i = 0; i < 8; ++i) {
                    const int ix = x + (i & 1), iy = y + ((i >> 1) & 1), iz = z + (i >> 2);
                    values[i] = field[v.voxel(ix, iy, iz)];
                    points[i] = QVector3D(float(v.axes[1][ix]), float(v.axes[2][iy]), float(v.axes[3][iz]));
                }
                const auto [min, max] = std::minmax_element(values.cbegin(), values.cend());
                if (*min >= level || *max < level ||
                    !std::all_of(values.cbegin(), values.cend(), [](double f) { return std::isfinite(f); }))
                    continue;
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
                        mesh.append(BondTriangle{edge(above[0], below[0]), edge(above[0], below[1]),
                                                 edge(above[0], below[2])});
                    else if (a == 3)
                        mesh.append(BondTriangle{edge(below[0], above[0]), edge(below[0], above[1]),
                                                 edge(below[0], above[2])});
                    else if (a == 2) {
                        const auto p = edge(above[0], below[0]), q = edge(above[0], below[1]),
                                   r = edge(above[1], below[0]), s = edge(above[1], below[1]);
                        mesh.append(BondTriangle{p, q, r});
                        mesh.append(BondTriangle{q, s, r});
                    }
                }
                if (mesh.size() > 500000)
                    fail(QStringLiteral("等值面过于复杂，请调整阈值或使用较粗的 X 网格。"));
            }
    }
    return mesh;
}
