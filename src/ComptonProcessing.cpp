#include "ComptonProcessing.h"
#include <QCoreApplication>
#include <QFile>
#include <QFileInfo>
#include <QJsonArray>
#include <QRegularExpression>
#include <QTextStream>
#include <algorithm>
#include <stdexcept>

namespace {
[[noreturn]] void fail(const QString &s) {
    throw std::runtime_error(s.toUtf8().constData());
}
void cancel(const std::atomic_bool &stop) {
    if (stop.load())
        fail(QCoreApplication::translate("ComptonProcessing", "已取消谱学分析。"));
}
bool inWindows(double x, const SpectrumWindows &windows) {
    return std::any_of(windows.cbegin(), windows.cend(),
                       [x](auto w) { return x >= w.first && x <= w.second; });
}
void checkWindows(const SpectrumWindows &windows, const QString &name) {
    if (windows.isEmpty())
        fail(QCoreApplication::translate("ComptonProcessing", "请填写%1。").arg(name));
    for (auto w : windows)
        if (!std::isfinite(w.first) || !std::isfinite(w.second) || w.second <= w.first)
            fail(QCoreApplication::translate("ComptonProcessing", "%1的起止值无效。").arg(name));
}

constexpr double electronKev = 510.99895069, alpha = 0.0072973525643, hartreeEv = 27.211386245988;
constexpr double electronRadius = 2.8179403205e-15;
struct Scattering {
    double e0, halfAngleFactor;
    explicit Scattering(SpectrumGeometry g)
        : e0(g.analyzerKev),
          halfAngleFactor(2 * std::pow(std::sin(g.angleDeg * 3.14159265358979323846 / 360.), 2)) {
        if (!(e0 > 0 && std::isfinite(e0) && g.angleDeg > 0 && g.angleDeg <= 180))
            fail(QCoreApplication::translate("ComptonProcessing", "缺少有效的分析器能量或散射角。"));
    }
    double qKev(double e1) const {
        return std::sqrt(std::pow(e1 - e0, 2) + 2 * e1 * e0 * halfAngleFactor);
    }
    double pz(double e1) const {
        return (.5 * qKev(e1) -
                (e1 - e0) * std::sqrt(.25 + electronKev * electronKev / (2 * e1 * e0 * halfAngleFactor))) /
               (electronKev * alpha);
    }
    double correction(double e1) const {
        const double p = pz(e1) * electronKev * alpha, q = qKev(e1), mass2 = electronKev * electronKev;
        const double electron = std::sqrt(mass2 + p * p);
        const double a = ((e1 - e0) * electron - e1 * e0 * halfAngleFactor) / q;
        const double d = (e1 - e0 + e0 * halfAngleFactor) * a / q;
        const double r = e1 * (electron - d), r2 = r - e1 * e0 * halfAngleFactor;
        const double delta = mass2 * (1 / r - 1 / r2);
        const double chi = r / r2 + r2 / r + 2 * delta + delta * delta;
        return 2 * e1 * q * electron / (mass2 * electronRadius * electronRadius * e0 * chi) *
               (1e-28 * electronKev * alpha);
    }
    double energyAt(double p, double lo, double hi) const {
        for (int i = 0; i < 38; ++i) {
            double mid = (lo + hi) / 2;
            if (pz(mid) > p)
                lo = mid;
            else
                hi = mid;
        }
        return ((lo + hi) / 2 - e0) * 1000;
    }
};

struct Shell {
    double binding = 0, electrons = 0;
    QVector<double> momentum, profile, derivative;

    void prepare() {
        if (momentum.size() < 2 || momentum.size() != profile.size() || momentum.first() != 0)
            fail(QCoreApplication::translate("ComptonProcessing", "HF 壳层动量网格无效。"));
        momentum.append(10000.);
        profile.append(0.);
        derivative.fill(0., profile.size());
        QVector<double> step, slope;
        for (qsizetype i = 1; i < momentum.size(); ++i) {
            auto h = momentum[i] - momentum[i - 1];
            if (!(h > 0))
                fail(QCoreApplication::translate("ComptonProcessing", "HF 壳层动量网格须递增。"));
            step.append(h);
            slope.append((profile[i] - profile[i - 1]) / h);
        }
        // The positive half is mirrored at pz=0, where the derivative is zero.
        for (qsizetype i = 1; i < profile.size() - 1; ++i) {
            if (slope[i - 1] * slope[i] <= 0)
                continue;
            const double a = 2 * step[i] + step[i - 1], b = step[i] + 2 * step[i - 1];
            derivative[i] = (a + b) / (a / slope[i - 1] + b / slope[i]);
        }
    }
    double value(double pz) const {
        const double p = std::abs(pz);
        const auto b = std::upper_bound(momentum.cbegin(), momentum.cend(), p);
        if (b == momentum.cend())
            return 0.;
        const auto i = b - momentum.cbegin() - 1;
        const double h = momentum[i + 1] - momentum[i], t = (p - momentum[i]) / h;
        const double t2 = t * t, t3 = t2 * t;
        return (2 * t3 - 3 * t2 + 1) * profile[i] + (t3 - 2 * t2 + t) * h * derivative[i] +
               (-2 * t3 + 3 * t2) * profile[i + 1] + (t3 - t2) * h * derivative[i + 1];
    }
};

QVector<Shell> readShells(const ComptonSettings &settings) {
    QMap<QString, double> weights;
    for (const auto &item :
         settings.composition.split(QRegularExpression("[,，;；\\s]+"), Qt::SkipEmptyParts)) {
        const auto parts = item.split(':');
        bool ok = true;
        double weight = parts.size() == 1 ? 1. : parts.value(1).toDouble(&ok);
        if (parts.size() > 2 || !ok || !std::isfinite(weight) || weight <= 0 || parts[0].isEmpty())
            fail(QCoreApplication::translate("ComptonProcessing", "组分格式示例 Li:1, Ni:0.9, Co:0.05, Mn:0.05, O:2。"));
        weights[parts[0]] += weight;
    }
    if (weights.isEmpty())
        fail(QCoreApplication::translate("ComptonProcessing", "请填写样品元素与原子计量数。"));
    QFile file(settings.hfTable);
    if (!file.open(QIODevice::ReadOnly))
        fail(QCoreApplication::translate("ComptonProcessing", "请选择可读取的 DABAX ComptonProfiles.dat。"));
    QTextStream stream(&file);
    const QRegularExpression space("\\s+");
    QVector<Shell> shells, element;
    QString symbol;
    auto commit = [&] {
        if (!weights.contains(symbol))
            return;
        if (element.isEmpty())
            fail(QCoreApplication::translate("ComptonProcessing", "HF 表缺少 %1 的壳层数据。").arg(symbol));
        for (auto &shell : element) {
            shell.prepare();
            shells.append(std::move(shell));
        }
        weights.remove(symbol);
    };
    while (!stream.atEnd()) {
        auto line = stream.readLine().trimmed();
        if (line.isEmpty())
            continue;
        auto fields = line.split(space, Qt::SkipEmptyParts);
        const auto tag = fields.first();
        if (tag == "#S") {
            commit();
            element.clear();
            symbol = fields.last();
        } else if (weights.contains(symbol)) {
            if (tag == "#UOCCUP") {
                element.resize(fields.size() - 1);
                for (qsizetype i = 0; i < element.size(); ++i) {
                    bool ok;
                    const double n = fields[i + 1].toDouble(&ok);
                    if (!ok || !std::isfinite(n) || n < 0)
                        fail(QCoreApplication::translate("ComptonProcessing", "HF 壳层占据数无效。"));
                    element[i].electrons = n * weights[symbol];
                }
            } else if (tag == "#UBIND") {
                if (fields.size() != element.size() + 1)
                    fail(QCoreApplication::translate("ComptonProcessing", "HF 壳层结合能列数无效。"));
                for (qsizetype i = 0; i < element.size(); ++i) {
                    bool ok;
                    element[i].binding = fields[i + 1].toDouble(&ok);
                    if (!ok || !std::isfinite(element[i].binding) || element[i].binding < 0)
                        fail(QCoreApplication::translate("ComptonProcessing", "HF 壳层结合能无效。"));
                }
            } else if (!tag.startsWith('#')) {
                if (fields.size() != element.size() + 2)
                    fail(QCoreApplication::translate("ComptonProcessing", "HF 壳层强度列数无效。"));
                bool pOk;
                double p = fields[0].toDouble(&pOk);
                if (!pOk || !std::isfinite(p))
                    fail(QCoreApplication::translate("ComptonProcessing", "HF 动量数据无效。"));
                for (qsizetype i = 0; i < element.size(); ++i) {
                    bool ok;
                    double v = fields[i + 2].toDouble(&ok);
                    if (!ok || !std::isfinite(v) || v < 0)
                        fail(QCoreApplication::translate("ComptonProcessing", "HF 壳层强度无效。"));
                    element[i].momentum.append(p);
                    element[i].profile.append(v);
                }
            }
        }
    }
    commit();
    if (!weights.isEmpty())
        fail(QCoreApplication::translate("ComptonProcessing", "HF 表缺少元素 %1。").arg(weights.keys().join(", ")));
    return shells;
}

double interpolate(const QVector<double> &axis, const QVector<double> &values, double point,
                   double outside = NAN) {
    if (axis.isEmpty() || point < axis.first() || point > axis.last())
        return outside;
    const auto b = std::lower_bound(axis.cbegin(), axis.cend(), point);
    const auto i = b - axis.cbegin();
    if (i == 0 || *b == point)
        return values[i];
    const double t = (point - axis[i - 1]) / (axis[i] - axis[i - 1]);
    return (1 - t) * values[i - 1] + t * values[i];
}

struct Prepared {
    QVector<double> pz, q, normalizedWide, normalizedTarget, core;
    double scale = 1.;
};
Prepared prepare(const QVector<double> &energy, const QVector<double> &wide, const QVector<double> &target,
                 SpectrumGeometry geometry, const QVector<Shell> &shells, const ComptonSettings &settings,
                 const std::atomic_bool &stop) {
    Scattering scattering(geometry);
    const double e0 = geometry.analyzerKev;
    double hi = e0 * 4;
    while (scattering.pz(hi) >= -100 && hi < e0 * 32)
        hi *= 2;
    if (!(scattering.pz(hi) < -100))
        fail(QCoreApplication::translate("ComptonProcessing", "分析器几何超出 HF 动量反算范围。"));
    // Use the reference program's 0.01 a.u. integration grid and shell f-sum convention.
    QVector<double> eGrid, pGrid, qGrid;
    const int first = std::min(9999, int(std::floor(scattering.pz(e0) * 100)));
    for (int index = first; index >= -10000; --index) {
        if (index % 256 == 0)
            cancel(stop);
        const double p = index / 100., loss = scattering.energyAt(p, e0, hi);
        eGrid.append(loss);
        pGrid.append(p);
        qGrid.append(scattering.qKev(e0 + loss / 1000.) / (electronKev * alpha));
    }
    QVector<double> total(eGrid.size(), 0.), core(eGrid.size(), 0.);
    for (const auto &shell : shells) {
        cancel(stop);
        QVector<double> values(eGrid.size(), 0.);
        double integral = 0;
        for (qsizetype i = 0; i < eGrid.size(); ++i) {
            if (eGrid[i] >= shell.binding)
                values[i] = shell.value(pGrid[i]);
            if (i)
                integral += (values[i - 1] * eGrid[i - 1] + values[i] * eGrid[i]) *
                            (eGrid[i] - eGrid[i - 1]) / (2 * hartreeEv * hartreeEv);
        }
        if (!(integral > 0))
            continue;
        for (qsizetype i = 0; i < values.size(); ++i) {
            double v = values[i] * shell.electrons * qGrid[i] * qGrid[i] / (2 * hartreeEv * integral);
            total[i] += v;
            if (shell.binding >= settings.valenceCutoffEv)
                core[i] += v;
        }
    }
    Prepared out;
    QVector<double> hf;
    for (qsizetype i = 0; i < energy.size(); ++i) {
        const double incident = e0 + energy[i] / 1000.;
        if (!(incident > 0))
            fail(QCoreApplication::translate("ComptonProcessing", "入射能量无效。"));
        const double factor = scattering.correction(incident);
        if (!(factor > 0 && std::isfinite(factor)))
            fail(QCoreApplication::translate("ComptonProcessing", "康普顿截面修正无效。"));
        out.pz.append(scattering.pz(incident));
        out.q.append(scattering.qKev(incident) / (electronKev * alpha));
        out.normalizedWide.append(wide[i] * factor);
        out.normalizedTarget.append(target[i] * factor);
        out.core.append(interpolate(eGrid, core, energy[i], energy[i] < 0 ? 0. : NAN));
        hf.append(interpolate(eGrid, total, energy[i], energy[i] < 0 ? 0. : NAN));
    }
    double observed = 0, expected = 0;
    for (qsizetype i = 1; i < energy.size(); ++i) {
        if (!inWindows(energy[i - 1], settings.normalization) ||
            !inWindows(energy[i], settings.normalization) || !std::isfinite(hf[i - 1]) ||
            !std::isfinite(hf[i]) || !std::isfinite(out.normalizedWide[i - 1]) ||
            !std::isfinite(out.normalizedWide[i]))
            continue;
        double step = (energy[i] - energy[i - 1]) / 2;
        observed += step * (out.normalizedWide[i - 1] + out.normalizedWide[i]);
        expected += step * (hf[i - 1] + hf[i]);
    }
    if (!(observed > 0 && expected > 0))
        fail(QCoreApplication::translate("ComptonProcessing", "宽谱归一化区间缺少正的实验或 HF 面积。"));
    out.scale = expected / observed;
    if (!std::isfinite(out.scale))
        fail(QCoreApplication::translate("ComptonProcessing", "HF 归一化系数无效。"));
    for (auto &v : out.normalizedWide)
        v *= out.scale;
    for (auto &v : out.normalizedTarget)
        v *= out.scale;
    return out;
}
struct Reference {
    QVector<double> pz, j, clean;
};
Reference referenceFor(const QVector<double> &e, const Prepared &data, const ComptonSettings &settings,
                       const std::atomic_bool &stop) {
    auto residual = data.normalizedWide;
    for (qsizetype i = 0; i < residual.size(); ++i)
        residual[i] -= data.core[i] * settings.coreScale;
    auto fit = SpectrumProcessing::fit(e, residual, {false, settings.coreTail}, stop);
    Reference out;
    for (qsizetype i = e.size(); i-- > 0;) {
        if (!out.pz.isEmpty() && data.pz[i] <= out.pz.last())
            fail(QCoreApplication::translate("ComptonProcessing", "该晶体的电子动量轴不单调。"));
        out.pz.append(data.pz[i]);
        out.j.append(data.q[i] * fit.corrected[i]);
        out.clean.append(inWindows(e[i], settings.cleanReference) ? out.j.last() : NAN);
    }
    return out;
}
double median(QVector<double> &values, qsizetype required) {
    if (values.size() < required)
        return NAN;
    std::sort(values.begin(), values.end());
    return (values[(values.size() - 1) / 2] + values[values.size() / 2]) / 2;
}
QJsonArray ranges(const SpectrumWindows &windows) {
    QJsonArray result;
    for (auto range : windows)
        result.append(QJsonArray{range.first, range.second});
    return result;
}

void subtractLinear(ComptonResult &out, const ProcessedSpectrum &spectrum, const QStringList &labels,
                    const SpectrumWindows &windows, const QString &method, const std::atomic_bool &stop) {
    for (const auto &label : labels) {
        cancel(stop);
        try {
            ComptonCrystal item;
            item.input = spectrum.crystals.value(label);
            item.fit = SpectrumProcessing::fit(spectrum.energyEv, item.input, {false, windows}, stop);
            item.fit.method = method;
            out.crystals.insert(label, std::move(item));
        } catch (const std::exception &ex) {
            out.failures[label] = QString::fromUtf8(ex.what());
        }
    }
    cancel(stop);
}
QMap<QString, Prepared> fitReference(const ProcessedSpectrum &wide, const ProcessedSpectrum &target,
                                     const QStringList &high, const ComptonSettings &settings,
                                     ComptonResult &out, const std::atomic_bool &stop) {
    const auto &e = wide.energyEv;
    checkWindows(settings.normalization, QCoreApplication::translate("ComptonProcessing", "HF 归一化区间"));
    if (settings.normalization.size() != 1)
        fail(QCoreApplication::translate("ComptonProcessing", "HF 归一化使用一个连续区间。"));
    checkWindows(settings.coreTail, QCoreApplication::translate("ComptonProcessing", "HF 芯谱尾部区间"));
    checkWindows(settings.cleanReference, QCoreApplication::translate("ComptonProcessing", "价电子轮廓拟合区间"));
    if (!std::isfinite(settings.coreScale) || settings.coreScale < 0 ||
        !std::isfinite(settings.valenceCutoffEv) || settings.valenceCutoffEv <= 0)
        fail(QCoreApplication::translate("ComptonProcessing", "HF 芯谱比例或价电子截止能量无效。"));
    const auto shells = readShells(settings);
    QMap<QString, Prepared> prepared;
    QVector<Reference> refs;
    for (const auto &label : high) {
        cancel(stop);
        try {
            auto data = prepare(e, wide.crystals.value(label), target.crystals.value(label),
                                wide.geometry.value(label), shells, settings, stop);
            try {
                refs.append(referenceFor(e, data, settings, stop));
                out.references.append(label);
            } catch (const std::exception &ex) {
                out.referenceFailures[label] = QString::fromUtf8(ex.what());
            }
            prepared.insert(label, std::move(data));
        } catch (const std::exception &ex) {
            out.failures[label] = QString::fromUtf8(ex.what());
        }
    }
    if (refs.isEmpty())
        fail(QCoreApplication::translate("ComptonProcessing", "没有可用于提取价电子轮廓的高 q 宽谱。"));
    QVector<double> clean;
    const auto required = (2 * refs.size() + 2) / 3;
    // As in xrs_ana, take the median across crystals with at least 2/3 reference coverage.
    for (int i = -1600; i < 1000; ++i) {
        if (i % 128 == 0)
            cancel(stop);
        const double p = i * .005;
        QVector<double> measured, selected;
        for (const auto &ref : refs) {
            double raw = interpolate(ref.pz, ref.j, p);
            double fit = interpolate(ref.pz, ref.clean, p);
            if (std::isfinite(raw))
                measured.append(raw);
            if (std::isfinite(fit))
                selected.append(fit);
        }
        out.pz.append(p);
        out.measuredJ.append(median(measured, required));
        clean.append(median(selected, required));
    }
    out.referenceFit =
        SpectrumProcessing::fit(out.pz, clean, {true, {{out.pz.first(), out.pz.last()}}}, stop);
    const auto &peak = out.referenceFit;
    out.parameters["pearson_pz"] = QJsonObject{{"center_au", peak.centerEv},
                                               {"fwhm_au", peak.fwhmEv},
                                               {"exponent", peak.exponent},
                                               {"height", peak.peakHeight},
                                               {"r_squared", peak.rSquared}};
    out.parameters["reference_crystals"] = QJsonArray::fromStringList(out.references);
    for (qsizetype i = 0; i < out.pz.size(); ++i) {
        double u = 2 * (out.pz[i] - peak.centerEv) / peak.fwhmEv;
        out.fittedJ.append(
            std::isfinite(out.measuredJ[i])
                ? peak.peakHeight *
                      std::pow(1 + std::expm1(std::log(2.) / peak.exponent) * u * u, -peak.exponent)
                : NAN);
    }
    return prepared;
}

double valenceAt(const ComptonResult &reference, double pz, double q) {
    return interpolate(reference.pz, reference.fittedJ, pz) / q;
}
} // namespace

ComptonResult analyzeLinear(const ProcessedSpectrum &spectrum, const SpectrumWindows &windows,
                            const std::atomic_bool &stop) {
    if (spectrum.crystals.isEmpty())
        fail(QCoreApplication::translate("ComptonProcessing", "请先选择逐晶体谱。"));
    checkWindows(windows, QCoreApplication::translate("ComptonProcessing", "线性拟合区间"));
    ComptonResult out;
    out.parameters = {{"background_mode", "linear_only"}, {"linear_ev", ranges(windows)}};
    subtractLinear(out, spectrum, spectrum.crystals.keys(), windows, QStringLiteral("Linear"), stop);
    return out;
}

ComptonResult analyzeCompton(const ProcessedSpectrum &wide, const ProcessedSpectrum &stitched,
                             const ComptonSettings &settings, const std::atomic_bool &stop) {
    if (wide.crystals.isEmpty() || wide.energyEv != stitched.energyEv || !std::isfinite(settings.highQMin) ||
        settings.highQMin < 0)
        fail(QCoreApplication::translate("ComptonProcessing", "请载入宽谱并填写有效的高 q 阈值。"));
    ComptonResult out;
    out.parameters = {{"background_mode", "momentum_split"},
                      {"high_q_min_inv_angstrom", settings.highQMin},
                      {"hf_table", QFileInfo(settings.hfTable).fileName()},
                      {"composition", settings.composition},
                      {"hf_core_scale", settings.coreScale},
                      {"valence_cutoff_ev", settings.valenceCutoffEv},
                      {"hf_normalization_ev", ranges(settings.normalization)},
                      {"core_tail_ev", ranges(settings.coreTail)},
                      {"clean_reference_ev", ranges(settings.cleanReference)},
                      {"high_residual_ev", ranges(settings.highResidual)},
                      {"low_linear_ev", ranges(settings.lowLinear)}};
    QStringList high, low;
    for (auto it = stitched.crystals.cbegin(); it != stitched.crystals.cend(); ++it) {
        const double q = wide.geometry.value(it.key()).qMean;
        if (!std::isfinite(q) || q <= 0)
            out.failures[it.key()] = QCoreApplication::translate("ComptonProcessing", "缺少唯一、有效的晶体 q 值。");
        else if (q >= settings.highQMin)
            high.append(it.key());
        else
            low.append(it.key());
    }
    const auto &e = wide.energyEv;
    if (!low.isEmpty())
        checkWindows(settings.lowLinear, QCoreApplication::translate("ComptonProcessing", "低 q 线性区间"));
    subtractLinear(out, stitched, low, settings.lowLinear, QStringLiteral("Low q linear"), stop);
    if (high.isEmpty())
        return out;
    try {
        checkWindows(settings.highResidual, QCoreApplication::translate("ComptonProcessing", "高 q 残余直线区间"));
        const auto prepared = fitReference(wide, stitched, high, settings, out, stop);
        for (auto it = prepared.cbegin(); it != prepared.cend(); ++it) {
            cancel(stop);
            try {
                ComptonCrystal item;
                item.high = true;
                item.scale = it->scale;
                item.input = it->normalizedTarget;
                item.hfCore = it->core;
                auto residual = item.input;
                for (qsizetype i = 0; i < e.size(); ++i) {
                    double valence = valenceAt(out, it->pz[i], it->q[i]);
                    item.valence.append(valence);
                    residual[i] -= valence;
                }
                item.fit = SpectrumProcessing::fit(e, residual, {false, settings.highResidual}, stop);
                for (qsizetype i = 0; i < e.size(); ++i)
                    item.fit.background[i] += item.valence[i];
                item.fit.method = "High q HF-normalized Compton valence + residual linear";
                out.crystals.insert(it.key(), std::move(item));
            } catch (const std::exception &ex) {
                out.failures[it.key()] = QString::fromUtf8(ex.what());
            }
        }
    } catch (const std::exception &ex) {
        for (const auto &label : high)
            if (!out.crystals.contains(label) && !out.failures.contains(label))
                out.failures[label] = QString::fromUtf8(ex.what());
    }
    cancel(stop);
    return out;
}

QMap<QString, ComptonTransfer> transferCompton(const ProcessedSpectrum &wide, const QStringList &labels,
                                                const QVector<double> &incidentKev, const ComptonSettings &settings,
                                                const std::atomic_bool &stop) {
    QStringList high;
    for (auto it = wide.crystals.cbegin(); it != wide.crystals.cend(); ++it) {
        const double q = wide.geometry.value(it.key()).qMean;
        if (std::isfinite(q) && q > 0 && q >= settings.highQMin)
            high.append(it.key());
    }
    ComptonResult reference;
    const auto prepared = fitReference(wide, wide, high, settings, reference, stop);
    QMap<QString, ComptonTransfer> transfers;
    for (const auto &label : labels) {
        const auto geometry = wide.geometry.value(label);
        const Scattering scattering(geometry);
        if (!prepared.contains(label))
            fail(QCoreApplication::translate("ComptonProcessing", "晶体 %1 的宽谱不可用。%2").arg(label, reference.failures.value(label)));

        ComptonTransfer out;
        out.analyzerKev = geometry.analyzerKev;
        const auto &source = prepared.constFind(label)->normalizedWide;
        for (double incident : incidentKev) {
            cancel(stop);
            if (!(std::isfinite(incident) && incident > 0))
                fail(QCoreApplication::translate("ComptonProcessing", "重建通道的入射能量无效。"));
            const double loss = 1000 * (incident - geometry.analyzerKev);
            const double factor = scattering.correction(incident);
            if (!(std::isfinite(factor) && factor > 0))
                fail(QCoreApplication::translate("ComptonProcessing", "重建通道的康普顿截面修正无效。"));
            out.correction.append(factor);
            out.reference.append(interpolate(wide.energyEv, source, loss));
            out.valence.append(valenceAt(reference, scattering.pz(incident),
                                         scattering.qKev(incident) / (electronKev * alpha)));
        }
        transfers.insert(label, std::move(out));
    }
    return transfers;
}

ComptonTransfer transferCompton(const ProcessedSpectrum &wide, const QString &label,
                                const QVector<double> &incidentKev, const ComptonSettings &settings,
                                const std::atomic_bool &stop) {
    return transferCompton(wide, QStringList{label}, incidentKev, settings, stop).value(label);
}
