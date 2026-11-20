#include "ImagingFit.h"
#include "ImagingIO.h"
#include <QTemporaryFile>
#include <algorithm>
#include <cmath>
#include <memory>
using namespace ImageUtil;
namespace {
// Rows are depth pairs; columns are physical scan positions. Spill large per-crystal
// caches to temporary storage so gain iterations do not repeatedly integrate detector pixels.
struct WeightedGram {
    QVector<int> pairs;
    ImageMatrix memory;
    std::shared_ptr<QTemporaryFile> disk;
    int sites;
    void visit(const std::atomic_bool &stop, const std::function<void(int, const ImageMatrix &)> &use) const {
        if (disk && !disk->seek(0))
            fail("Cannot rewind weighted normal matrix");
        for (int begin = 0; begin < pairs.size(); begin += 64) {
            cancel(stop);
            int count = std::min(64, int(pairs.size()) - begin);
            ImageMatrix tile(count, sites);
            if (disk) {
                qint64 bytes = tile.size() * qint64(sizeof(double));
                if (disk->read(reinterpret_cast<char *>(tile.data()), bytes) != bytes)
                    fail("Cannot read temporary weighted normal matrix");
            } else
                tile = memory.middleRows(begin, count);
            use(begin, tile);
        }
    }
};
WeightedGram buildGram(const ImageMatrix &op, const ImageMatrix &w, const ImageMatrix &geometry,
                       qint64 &budget, const std::atomic_bool &stop) {
    WeightedGram cache;
    cache.sites = int(w.cols());
    QVector<QPair<int, int>> ij;
    for (int i = 0, k = 0; i < op.cols(); ++i)
        for (int j = i; j < op.cols(); ++j, ++k)
            if (geometry(i, j) > 0.) {
                cache.pairs.append(k);
                ij.append({i, j});
            }
    const qint64 bytes = qint64(ij.size()) * w.cols() * sizeof(double);
    if (bytes <= budget) {
        cache.memory.resize(ij.size(), w.cols());
        budget -= bytes;
    } else {
        cache.disk = std::make_shared<QTemporaryFile>();
        if (!cache.disk->open())
            fail("Cannot create temporary weighted normal matrix");
    }
    for (int begin = 0; begin < ij.size(); begin += 64) {
        cancel(stop);
        const int count = std::min(64, int(ij.size()) - begin);
        ImageMatrix product(count, op.rows());
        for (int s = 0; s < count; ++s) {
            const auto [i, j] = ij[begin + s];
            product.row(s) = (op.col(i).array() * op.col(j).array()).matrix().transpose();
        }
        ImageMatrix tile = product * w;
        if (cache.disk) {
            qint64 size = tile.size() * qint64(sizeof(double));
            if (cache.disk->write(reinterpret_cast<const char *>(tile.data()), size) != size)
                fail("Insufficient temporary disk space for weighted fit");
        } else
            cache.memory.middleRows(begin, count) = tile;
    }
    return cache;
}
} // namespace
OperatorFit fitOperators(const QVector<FitObservation> &obs, const QJsonObject &rec, const QJsonObject &stats,
                         const std::atomic_bool &stop, Imaging::Progress log) {
    if (obs.isEmpty())
        fail("No crystals in reconstruction");
    const int nx = int(obs.first().response->cols()), sites = int(obs.first().data->cols());
    const bool weighted = stats["variance_weighted"].toBool();
    const int pairs = nx * (nx + 1) / 2;
    if (weighted && qint64(pairs) * sites > 67108864)
        fail("Weighted normal matrix exceeds 512 MiB; increase the X step");
    OperatorFit out;
    QVector<ImageMatrix> geometry, backs;
    QVector<WeightedGram> caches;
    QVector<double> means;
    double pooled = 0.;
    qint64 positive = 0;
    if (weighted) {
        for (const auto &o : obs) {
            const auto &v = *o.variance;
            if (v.rows() != o.data->rows() || v.cols() != sites || !v.allFinite() || v.minCoeff() < 0.)
                fail("Counting variance missing, negative or incompatible with data");
            const auto count = (v.array() > 0).count();
            const double sum = v.sum();
            means.append(count ? sum / count : 0.);
            pooled += sum;
            positive += count;
        }
    }
    qint64 budget = 512LL * 1024 * 1024;
    int anchor = -1;
    for (int a = 0; a < obs.size(); ++a) {
        cancel(stop);
        const auto &o = obs[a];
        const auto &op = *o.response, &d = *o.data;
        if (o.key == rec["anchor_roi"].toString())
            anchor = a;
        geometry.append(op.transpose() * op);
        if (weighted) {
            if (log)
                log(QString("Weighted normal matrix %1/%2 · %3").arg(a + 1).arg(obs.size()).arg(o.key), 62);
            const double sigma = stats["sigma_floor"].toDouble();
            const double floor =
                o.savedVarianceFloor > 0. ? o.savedVarianceFloor
                : stats["sigma_floor"].isDouble()
                    ? sigma * sigma
                    : .01 * (means[a] > 0 ? means[a] : pooled / std::max<qint64>(1, positive));
            if (!(floor > 0 && std::isfinite(floor)))
                fail("No positive counting scale; specify a positive sigma floor");
            out.varianceFloors.append(floor);
            ImageMatrix w = (o.variance->array() + floor).inverse();
            backs.append(op.transpose() * (w.array() * d.array()).matrix());
            caches.append(buildGram(op, w, geometry.last(), budget, stop));
        } else
            backs.append(op.transpose() * d);
    }
    const bool gains = rec["fit_gains"].toBool();
    if (gains && anchor < 0)
        fail("Gain anchor is absent from accepted crystals: " + rec["anchor_roi"].toString());
    out.gains.fill(1., obs.size());
    auto solve = [&] {
        out.h = ImageMatrix::Zero(nx, nx);
        out.b = ImageMatrix::Zero(nx, sites);
        if (weighted)
            out.weightedH = ImageMatrix::Zero(pairs, sites);
        for (int a = 0; a < obs.size(); ++a) {
            cancel(stop);
            const double g = out.gains[a];
            out.h.noalias() += g * g * geometry[a];
            out.b.noalias() += g * backs[a];
            if (weighted)
                caches[a].visit(stop, [&](int begin, const ImageMatrix &tile) {
                    for (int k = 0; k < tile.rows(); ++k)
                        out.weightedH.row(caches[a].pairs[begin + k]) += g * g * tile.row(k);
                });
        }
        out.rho = Imaging::nnls(out.normal(), out.b, rec["nnls_max_iterations"].toInt(),
                                rec["nnls_tolerance"].toDouble(), stop, &out.nnls, out.rho);
    };
    const auto limits = rec["gain_bounds"].toArray();
    for (int iter = 0; gains && iter < rec["max_gain_iterations"].toInt(); ++iter) {
        solve();
        auto next = out.gains;
        ImageMatrix dyads;
        if (weighted) {
            dyads.resize(pairs, sites);
            for (int i = 0, k = 0; i < nx; ++i)
                for (int j = i; j < nx; ++j, ++k)
                    dyads.row(k).array() =
                        (i == j ? 1. : 2.) * out.rho.row(i).array() * out.rho.row(j).array();
        }
        double change = 0., norm = 0.;
        for (int a = 0; a < obs.size(); ++a) {
            if (a != anchor) {
                double denominator = 0.;
                if (weighted)
                    caches[a].visit(stop, [&](int begin, const ImageMatrix &tile) {
                        for (int k = 0; k < tile.rows(); ++k)
                            denominator += tile.row(k).dot(dyads.row(caches[a].pairs[begin + k]));
                    });
                else
                    denominator = (out.rho.array() * (geometry[a] * out.rho).array()).sum();
                if (denominator > 0.)
                    next[a] = std::clamp((backs[a].array() * out.rho.array()).sum() / denominator,
                                         limits[0].toDouble(), limits[1].toDouble());
            }
            change += std::pow(next[a] - out.gains[a], 2);
            norm += out.gains[a] * out.gains[a];
        }
        out.gains = next;
        change = std::sqrt(change / std::max(norm, 1e-30));
        out.history.append(change);
        if (log && (iter == 0 || (iter + 1) % 20 == 0))
            log(QString("Gain iteration %1: relative change %2").arg(iter + 1).arg(change), 65);
        if (change < rec["gain_tolerance"].toDouble()) {
            out.gainConverged = true;
            break;
        }
    }
    solve();
    out.gainConverged = (!gains || out.gainConverged) && out.nnls.converged;
    return out;
}
