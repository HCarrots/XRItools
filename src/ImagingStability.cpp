#include "ImagingFit.h"
#include "ImagingIO.h"
#include <QCoreApplication>
#include <algorithm>
#include <cmath>
#include <numeric>
using namespace ImageUtil;
namespace {
QString packed(const QByteArray &v) {
    return QString::fromLatin1(qCompress(v, 6).toBase64());
}
QJsonObject measureAxis(const QVector<double> &axis, const QVector<double> &profile,
                        const QVector<bool> &valid) {
    const int peak = int(std::max_element(profile.begin(), profile.end()) - profile.begin());
    const double max = profile[peak], mass = std::accumulate(profile.begin(), profile.end(), 0.);
    if (!(mass > 0 && max > 0))
        return {};
    QJsonObject result{
        {"center_mm", std::inner_product(axis.begin(), axis.end(), profile.begin(), 0.) / mass}};
    const double half = max * .5;
    int lo = peak, hi = peak;
    while (lo > 0 && valid[lo - 1] && profile[lo - 1] >= half)
        --lo;
    while (hi + 1 < axis.size() && valid[hi + 1] && profile[hi + 1] >= half)
        ++hi;
    auto cross = [&](int i, int j) {
        return axis[i] + (half - profile[i]) * (axis[j] - axis[i]) / (profile[j] - profile[i]);
    };
    if (lo > 0 && valid[lo - 1])
        result["left_mm"] = cross(lo, lo - 1);
    if (hi + 1 < axis.size() && valid[hi + 1])
        result["right_mm"] = cross(hi, hi + 1);
    if (result.contains("left_mm") && result.contains("right_mm"))
        result["thickness_mm"] = result["right_mm"].toDouble() - result["left_mm"].toDouble();
    return result;
}
} // namespace
QJsonObject compareCrystalGroups(const QVector<FitObservation> &all, const QVector<QVector<double>> &axes,
                                 const QJsonObject &c, double lambda, const std::atomic_bool &stop,
                                 Imaging::Progress log) {
    if (all.size() < 2)
        fail("At least two accepted crystals required for independent groups");
    const int nx = int(axes[0].size()), ny = int(axes[1].size()), nz = int(axes[2].size());
    const std::array<double, 3> spacing{axes[0][1] - axes[0][0], axes[1][1] - axes[1][0],
                                        axes[2][1] - axes[2][0]};
    QVector<Eigen::ArrayXd> power;
    Eigen::ArrayXd total = Eigen::ArrayXd::Zero(nx);
    for (const auto &o : all) {
        cancel(stop);
        power.append(o.response->colwise().squaredNorm().transpose().array());
        total += power.last();
    }
    if (!(total.maxCoeff() > 0.))
        fail("No observable depth range");
    const double fullFloor = total.maxCoeff() * 1e-4; // 1% of peak amplitude sensitivity
    QVector<int> order(all.size());
    std::iota(order.begin(), order.end(), 0);
    QVector<double> strength;
    for (const auto &p : power)
        strength.append((p / total.max(fullFloor)).sum());
    std::stable_sort(order.begin(), order.end(), [&](int a, int b) { return strength[a] > strength[b]; });
    std::array<QVector<int>, 2> groups;
    std::array<Eigen::ArrayXd, 2> sums{Eigen::ArrayXd::Zero(nx), Eigen::ArrayXd::Zero(nx)};
    for (int i : order) {
        const auto diff = ((sums[0] - sums[1]) * power[i] / total.max(fullFloor).square()).sum();
        const int g = groups[0].size() >= (all.size() + 1) / 2   ? 1
                      : groups[1].size() >= (all.size() + 1) / 2 ? 0
                                                                 : (diff <= 0 ? 0 : 1);
        groups[g].append(i);
        sums[g] += power[i];
    }
    QVector<bool> common(nx);
    QJsonArray support;
    for (int i = 0; i < nx; ++i) {
        common[i] = total[i] >= fullFloor && std::min(sums[0][i], sums[1][i]) >= .01 * total[i];
        support.append(bool(common[i]));
    }
    bool adjacent = false;
    for (int i = 1; i < nx; ++i)
        adjacent |= common[i] && common[i - 1];
    if (!adjacent)
        fail("Groups have insufficient overlapping depth coverage for boundary comparison");
    const auto tv = c["tv"].toObject(), stats = c["statistics"].toObject();
    const bool fitGains = c["reconstruction"].toObject()["fit_gains"].toBool();
    QJsonArray runs, members;
    QByteArray votes(nx * ny * nz, 0), tentative(votes);
    int converged = 0, measured = 0;
    for (int g = 0; g < 2; ++g) {
        QVector<FitObservation> obs;
        QStringList keys;
        for (int i : groups[g]) {
            obs.append(all[i]);
            keys.append(all[i].key);
        }
        members.append(json(keys));
        auto rec = c["reconstruction"].toObject();
        int anchor = keys.indexOf(rec["anchor_roi"].toString());
        if (anchor < 0)
            anchor = 0;
        QVector<int> anchors{anchor};
        if (fitGains && keys.size() > 1)
            anchors.append(anchor == 0 ? 1 : 0);
        for (int pass = 0; pass < anchors.size(); ++pass) {
            cancel(stop);
            const auto reference = keys[anchors[pass]];
            rec["anchor_roi"] = reference;
            QString label =
                QString("%1 / %2").arg(g == 0 ? "A" : "B", fitGains ? reference : QStringLiteral("g=1"));
            auto notify = [&](const QString &s, int phase) {
                if (log)
                    log(label + " · " + s, std::min(99, g * 50 + pass * 25 + phase / 4));
            };
            notify("NNLS", 0);
            auto fit = fitOperators(obs, rec, stats, stop, notify);
            if (!(lambda > 0.)) {
                lambda = tv["lambda"].toDouble();
                if (!(lambda > 0.))
                    lambda = tv["relative_fraction"].toDouble() * percentile(vector(fit.b.cwiseAbs()), 95) *
                             (*std::min_element(spacing.begin(), spacing.end()));
            }
            QVector<double> factors{0.};
            if (lambda > 0.) {
                if (!pass)
                    factors.append(.5);
                factors.append(1.);
                if (!pass)
                    factors.append(2.);
            }
            QJsonObject gains;
            for (int a = 0; a < keys.size(); ++a)
                gains[keys[a]] = fit.gains[a];
            for (double factor : factors) {
                cancel(stop);
                ImagingSolverInfo info;
                ImageMatrix rho;
                if (factor > 0.) {
                    notify(QString("TV × %1").arg(factor), 70);
                    rho = Imaging::tv(fit.normal(), fit.b, ny, nz, spacing, factor * lambda,
                                      tv["max_iterations"].toInt(), tv["tolerance"].toDouble(), fit.rho, stop,
                                      &info, notify);
                } else
                    rho = fit.rho;
                const bool ok = fit.gainConverged && fit.nnls.converged && (!factor || info.converged);
                QJsonObject row{{"group", g == 0 ? "A" : "B"},
                                {"anchor", fitGains ? QJsonValue(reference) : QJsonValue()},
                                {"gain_converged", fit.gainConverged},
                                {"nnls_converged", fit.nnls.converged},
                                {"gains", gains},
                                {"gain_iterations", fit.history.size()},
                                {"last_gain_change",
                                 fit.history.isEmpty() ? QJsonValue() : QJsonValue(fit.history.last())},
                                {"nnls_iterations", fit.nnls.iterations},
                                {"tv_iterations", info.iterations},
                                {"tv_factor", factor},
                                {"lambda", factor * lambda},
                                {"tv_converged", factor ? QJsonValue(info.converged) : QJsonValue()},
                                {"converged", ok}};
                QVector<QVector<double>> profiles{QVector<double>(nx), QVector<double>(ny),
                                                  QVector<double>(nz)};
                double peak = 0;
                for (int i = 0; i < nx; ++i)
                    if (common[i]) {
                        peak = std::max(peak, rho.row(i).maxCoeff());
                        for (int y = 0; y < ny; ++y)
                            for (int z = 0; z < nz; ++z) {
                                const double v = rho(i, y * nz + z);
                                profiles[0][i] += v;
                                profiles[1][y] += v;
                                profiles[2][z] += v;
                            }
                    }
                QByteArray mask(votes.size(), 0);
                if (peak > 0.) {
                    ++measured;
                    if (ok)
                        ++converged;
                    for (int i = 0; i < nx; ++i)
                        if (common[i])
                            for (int col = 0; col < ny * nz; ++col)
                                if (rho(i, col) >= .5 * peak) {
                                    const int k = i * ny * nz + col;
                                    mask[k] = 1;
                                    ++tentative[k];
                                    if (ok)
                                        ++votes[k];
                                }
                    for (int d = 0; d < 3; ++d)
                        row[QString("XYZ")[d]] = measureAxis(
                            axes[d], profiles[d], d == 0 ? common : QVector<bool>(axes[d].size(), true));
                }
                row["surface_mask"] = packed(mask);
                row["has_signal"] = peak > 0.;
                runs.append(row);
            }
        }
    }
    cancel(stop);
    return {{"schema", "XRI.boundary-stability.v1"},
            {"groups", members},
            {"runs", runs},
            {"common_x", support},
            {"surface_fraction", .5},
            {"reference_lambda", lambda},
            {"shape", QJsonArray{nx, ny, nz}},
            {"converged_runs", converged},
            {"signal_runs", measured},
            {"converged_votes", packed(votes)},
            {"all_votes", packed(tentative)}};
}
