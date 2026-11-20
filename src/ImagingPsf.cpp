#include "ImagingIO.h"
#include <Eigen/Cholesky>
#include <algorithm>
#include <cmath>
#include <numeric>
using namespace ImageUtil;
namespace {
using Parameters = Eigen::Matrix<double, 7, 1>;
double gaussian(double row, double col, const Parameters &p, Parameters *jac = nullptr) {
    double u = (row - p[1]) / p[3], v = (col - p[2]) / p[4], c = p[5], d = 1 - c * c,
           n = u * u - 2 * c * u * v + v * v, g = std::exp(-n / (2 * d)), peak = p[0] * g;
    if (jac) {
        (*jac) << g, peak * (u - c * v) / (d * p[3]), peak * (v - c * u) / (d * p[4]),
            peak * u * (u - c * v) / (d * p[3]), peak * v * (v - c * u) / (d * p[4]),
            peak * (u * v * d - c * n) / (d * d), 1;
    }
    return peak + p[6];
}
bool fit(const QVector<QPoint> &pixels, const QVector<double> &image, const QJsonObject &settings,
         Parameters *out, double *r2, const std::atomic_bool &stop) {
    const auto bounds = settings["sigma_bounds_pixels"].toArray();
    double scale = *std::max_element(image.begin(), image.end());
    if (!(scale > 0))
        return false;
    QVector<double> y;
    for (auto v : image)
        y.append(v / scale);
    int maximum = int(std::max_element(y.begin(), y.end()) - y.begin());
    double bg = std::max(0., percentile(y, 20));
    double rlo = INFINITY, rhi = -INFINITY, clo = INFINITY, chi = -INFINITY;
    for (auto p : pixels) {
        rlo = std::min(rlo, double(p.y()));
        rhi = std::max(rhi, double(p.y()));
        clo = std::min(clo, double(p.x()));
        chi = std::max(chi, double(p.x()));
    }
    if (rhi <= rlo || chi <= clo)
        return false;
    Parameters p, lo, hi;
    p << std::max(y[maximum] - bg, 1e-10 / scale), pixels[maximum].y(), pixels[maximum].x(),
        std::clamp(2., bounds[0].toDouble(), bounds[1].toDouble()),
        std::clamp(2., bounds[0].toDouble(), bounds[1].toDouble()), 0, bg;
    lo << 0, rlo, clo, bounds[0].toDouble(), bounds[0].toDouble(), -.9, 0;
    hi << INFINITY, rhi, chi, bounds[1].toDouble(), bounds[1].toDouble(), .9, INFINITY;
    auto error = [&](const Parameters &v) {
        double sum = 0;
        for (int i = 0; i < pixels.size(); ++i) {
            double diff = gaussian(pixels[i].y(), pixels[i].x(), v) - y[i];
            sum += diff * diff;
        }
        return sum;
    };
    double current = error(p), lambda = .001;
    bool success = false;
    int evaluations = 1;
    while (evaluations < settings["max_fit_evaluations"].toInt()) {
        if (evaluations % 20 == 1)
            cancel(stop);
        Eigen::Matrix<double, 7, 7> h = Eigen::Matrix<double, 7, 7>::Zero();
        Parameters b = Parameters::Zero();
        for (int i = 0; i < pixels.size(); ++i) {
            Parameters jac;
            double residual = y[i] - gaussian(pixels[i].y(), pixels[i].x(), p, &jac);
            h.noalias() += jac * jac.transpose();
            b.noalias() += jac * residual;
        }
        Parameters projected = b;
        for (int k = 0; k < 7; ++k)
            if ((p[k] <= lo[k] + 1e-12 && b[k] < 0) || (p[k] >= hi[k] - 1e-12 && b[k] > 0))
                projected[k] = 0;
        if (projected.lpNorm<Eigen::Infinity>() < 1e-10) {
            success = true;
            break;
        }
        for (int k = 0; k < 7; ++k)
            h(k, k) += lambda * (h(k, k) + 1e-10);
        Parameters delta = h.ldlt().solve(b);
        if (!delta.allFinite())
            return false;
        Parameters next = (p + delta).cwiseMax(lo).cwiseMin(hi);
        double trial = error(next);
        ++evaluations;
        if (std::isfinite(trial) && trial < current) {
            double improvement = current - trial;
            p = next;
            current = trial;
            lambda = std::max(1e-12, lambda * .3);
            if (improvement < 1e-13 * (1 + current)) {
                success = true;
                break;
            }
        } else {
            lambda *= 10;
            if (lambda > 1e18)
                break;
        }
    }
    if (!success)
        return false;
    double mean = std::accumulate(y.begin(), y.end(), 0.) / y.size(), ss = 0;
    for (double v : y)
        ss += (v - mean) * (v - mean);
    *r2 = ss > 0 ? 1 - current / ss : 0;
    p[0] *= scale;
    p[6] *= scale;
    *out = p;
    return p.allFinite();
}
ImagingPsf baseline(const ImageMatrix &response, const QVector<QPoint> &pixels, const QVector<double> &x,
                    const QJsonObject &settings, const std::atomic_bool &stop,
                    QHash<double, ImagingSpotFit> &cache, const QVector<double> &height, int held) {
    double maxHeight = 0.;
    for (int j = 0; j < x.size(); ++j)
        if (j != held)
            maxHeight = std::max(maxHeight, height[j]);
    const double threshold = maxHeight * settings["min_peak_fraction"].toDouble();
    QVector<Parameters> fits;
    QVector<double> validX, r2s;
    double rlo = INFINITY, rhi = -INFINITY, clo = INFINITY, chi = -INFINITY;
    for (auto pixel : pixels) {
        rlo = std::min(rlo, double(pixel.y()));
        rhi = std::max(rhi, double(pixel.y()));
        clo = std::min(clo, double(pixel.x()));
        chi = std::max(chi, double(pixel.x()));
    }
    for (int j = 0; j < x.size(); ++j) {
        cancel(stop);
        if (j == held || height[j] <= threshold)
            continue;
        ImagingSpotFit spot;
        if (cache.contains(x[j]))
            spot = cache.value(x[j]);
        else {
            QVector<double> image;
            for (int p = 0; p < pixels.size(); ++p)
                image.append(response(p, j));
            Parameters fitted;
            spot.fitted = fit(pixels, image, settings, &fitted, &spot.rSquared, stop);
            if (spot.fitted)
                std::copy_n(fitted.data(), 7, spot.parameters.begin());
            cache.insert(x[j], spot);
        }
        if (!spot.fitted)
            continue;
        const Parameters fitted = Eigen::Map<const Parameters>(spot.parameters.data());
        const double r2 = spot.rSquared;
        double margin = settings["edge_margin_pixels"].toDouble();
        if (r2 >= settings["min_r_squared"].toDouble() && fitted[1] - rlo > margin &&
            rhi - fitted[1] > margin && fitted[2] - clo > margin && chi - fitted[2] > margin) {
            fits.append(fitted);
            validX.append(x[j]);
            r2s.append(r2);
        }
    }
    if (fits.size() < settings["min_observations"].toInt())
        fail(QString("Only %1 usable PSF observations; cannot fit a reliable trajectory").arg(fits.size()));
    // np.polyfit(..., w=sqrt(amplitude)): least-squares weights are amplitude.
    double weight = 0, xMean = 0, rowMean = 0, colMean = 0;
    for (int i = 0; i < fits.size(); ++i) {
        double w = fits[i][0];
        weight += w;
        xMean += w * validX[i];
        rowMean += w * fits[i][1];
        colMean += w * fits[i][2];
    }
    xMean /= weight;
    rowMean /= weight;
    colMean /= weight;
    double xx = 0, xr = 0, xc = 0;
    for (int i = 0; i < fits.size(); ++i) {
        double t = validX[i] - xMean, w = fits[i][0];
        xx += w * t * t;
        xr += w * t * (fits[i][1] - rowMean);
        xc += w * t * (fits[i][2] - colMean);
    }
    if (!(xx > 0))
        fail("No valid X spread for PSF trajectory");
    ImagingPsf m;
    m.rowSlope = xr / xx;
    m.columnSlope = xc / xx;
    m.rowIntercept = rowMean - m.rowSlope * xMean;
    m.columnIntercept = colMean - m.columnSlope * xMean;
    double speed = std::hypot(m.rowSlope, m.columnSlope);
    if (!(speed > 0))
        fail("No depth-dependent detector displacement");
    m.sampleMmPerPixel = 1 / speed;
    double squared = 0;
    QVector<double> amps, rows, cols, corr, offsets;
    for (int i = 0; i < fits.size(); ++i) {
        auto &p = fits[i];
        double dr = m.rowSlope * validX[i] + m.rowIntercept - p[1],
               dc = m.columnSlope * validX[i] + m.columnIntercept - p[2];
        squared += dr * dr + dc * dc;
        amps.append(p[0]);
        rows.append(p[3]);
        cols.append(p[4]);
        corr.append(p[5]);
        offsets.append(p[6]);
        m.observations.append(QJsonObject{{"x", validX[i]},
                                          {"row", p[1]},
                                          {"column", p[2]},
                                          {"r_squared", r2s[i]},
                                          {"calibration_offset", p[6]}});
    }
    m.trajectoryRms = std::sqrt(squared / fits.size());
    if (m.trajectoryRms > settings["max_trajectory_rms_pixels"].toDouble())
        fail(QString("Centroid trajectory residual %1 pixels exceeds %2; possible mixed spots")
                 .arg(m.trajectoryRms)
                 .arg(settings["max_trajectory_rms_pixels"].toDouble()));
    m.amplitude = percentile(amps, 50);
    m.sigmaRow = percentile(rows, 50);
    m.sigmaColumn = percentile(cols, 50);
    m.correlation = percentile(corr, 50);
    m.calibrationOffset = percentile(offsets, 50);
    m.xMin = *std::min_element(validX.begin(), validX.end());
    m.xMax = *std::max_element(validX.begin(), validX.end());
    m.minR2 = *std::min_element(r2s.begin(), r2s.end());
    double projected =
        std::sqrt(std::pow(m.rowSlope * m.sigmaRow, 2) + std::pow(m.columnSlope * m.sigmaColumn, 2) +
                  2 * m.rowSlope * m.columnSlope * m.correlation * m.sigmaRow * m.sigmaColumn) /
        speed;
    m.fwhmMm = 2.354820045 * projected / speed;
    return m;
}

constexpr std::array<int, 3> shapeIndices{0, 3, 4};
const std::array<QString, 3> modelNames{"constant_gaussian", "varying_gaussian", "measured_template"};

void fitTrends(ImagingPsf &m, const QHash<double, ImagingSpotFit> &cache) {
    m.trendCenter = 0.;
    for (const auto &obs : m.observations)
        m.trendCenter += obs.toObject()["x"].toDouble() / m.observations.size();
    for (int i = 0; i < 3; ++i) {
        double mean = 0., xx = 0., xy = 0., low = INFINITY, high = -INFINITY;
        for (const auto &obs : m.observations) {
            const double x = obs.toObject()["x"].toDouble(), t = x - m.trendCenter;
            const double value = std::log(std::max(1e-30, cache.value(x).parameters[shapeIndices[i]]));
            mean += value / m.observations.size();
            xx += t * t;
            xy += t * value;
            low = std::min(low, value);
            high = std::max(high, value);
        }
        m.trends[i] = {mean, xy / xx, low, high};
    }
}

// Missing ROI neighbours contribute zero residual, retaining the analytic Gaussian at masked pixels.
template <class Value>
double interpolate(const QHash<QPoint, int> &index, double r, double c, Value value,
                   double *coverage = nullptr) {
    const int row = int(std::floor(r)), col = int(std::floor(c));
    const double dr = r - row, dc = c - col;
    double sum = 0., weight = 0.;
    for (int i : {0, 1})
        for (int j : {0, 1}) {
            const double w = (i ? dr : 1 - dr) * (j ? dc : 1 - dc);
            if (!w)
                continue;
            const auto it = index.constFind(QPoint(col + j, row + i));
            if (it != index.cend()) {
                sum += w * value(it.value());
                weight += w;
            }
        }
    if (coverage)
        *coverage = weight;
    return sum;
}

// Relative antisymmetric component around the fitted centroid, measured on paired ROI support.
double spotAsymmetry(const ImageMatrix &response, int j, const QVector<QPoint> &pixels,
                     const QHash<QPoint, int> &index, const ImagingSpotFit &spot) {
    const auto &p = spot.parameters;
    auto signal = [&](int k) { return std::max(0., response(k, j) - p[6]); };
    double odd = 0., even = 0.;
    for (int k = 0; k < pixels.size(); ++k) {
        const double r = pixels[k].y(), c = pixels[k].x();
        if (std::pow((r - p[1]) / p[3], 2) + std::pow((c - p[2]) / p[4], 2) > 9.)
            continue;
        double coverage;
        const double mirrored = interpolate(index, 2 * p[1] - r, 2 * p[2] - c, signal, &coverage);
        if (coverage < .999999)
            continue;
        odd += std::pow(signal(k) - mirrored, 2);
        even += std::pow(signal(k) + mirrored, 2);
    }
    return even > 0 ? std::sqrt(odd / even) : 0.;
}

bool prepareTemplates(ImagingPsf &m, const ImageMatrix &response, const QVector<QPoint> &pixels,
                      const QVector<double> &x, const QHash<QPoint, int> &index,
                      QHash<double, ImagingSpotFit> &cache, int minimum) {
    if (m.observations.size() < minimum)
        return false;
    QVector<double> asymmetry;
    QVector<int> columns;
    for (const auto &obs : m.observations) {
        const double pos = obs.toObject()["x"].toDouble();
        int j = int(std::lower_bound(x.cbegin(), x.cend(), pos) - x.cbegin());
        auto &spot = cache[pos];
        if (spot.asymmetry < 0.)
            spot.asymmetry = spotAsymmetry(response, j, pixels, index, spot);
        columns.append(j);
        asymmetry.append(spot.asymmetry);
    }
    const double score = percentile(asymmetry, 50);
    m.selection["median_asymmetry"] = score;
    if (score < .12)
        return false;
    m.pixelIndex = index;
    m.templateResiduals.resize(pixels.size(), columns.size());
    for (int n = 0; n < columns.size(); ++n) {
        const int j = columns[n];
        const double pos = x[j], offset = cache.value(pos).parameters[6];
        m.knots.append(pos);
        Parameters p;
        p << m.amplitude, m.rowSlope * pos + m.rowIntercept, m.columnSlope * pos + m.columnIntercept,
            m.sigmaRow, m.sigmaColumn, m.correlation, 0.;
        for (int k = 0; k < pixels.size(); ++k)
            m.templateResiduals(k, n) =
                std::max(0., response(k, j) - offset) - gaussian(pixels[k].y(), pixels[k].x(), p);
    }
    return true;
}
} // namespace

ImagingPsf Imaging::fitPsf(const ImageMatrix &response, const QVector<QPoint> &pixels,
                           const QVector<double> &x, const QJsonObject &settings,
                           const std::atomic_bool &stop, QHash<double, ImagingSpotFit> *spotCache) {
    if (pixels.isEmpty() || x.isEmpty() || response.rows() != pixels.size() || response.cols() != x.size() ||
        !response.allFinite() || response.minCoeff() < 0 || !std::is_sorted(x.begin(), x.end()))
        fail("Expected ordered foil positions and finite nonnegative PSF response[p,x]");
    QVector<double> height;
    for (int j = 0; j < x.size(); ++j) {
        cancel(stop);
        const auto values = vector(response.col(j));
        height.append(percentile(values, 99.5) - percentile(values, 20));
    }
    QHash<double, ImagingSpotFit> local;
    auto &cache = spotCache ? *spotCache : local;
    auto train = [&](int held) { return baseline(response, pixels, x, settings, stop, cache, height, held); };
    auto model = train(-1);
    if (!settings["auto_select_model"].toBool()) {
        model.selection["status"] = "disabled";
        return model;
    }
    QHash<QPoint, int> index;
    for (int k = 0; k < pixels.size(); ++k)
        index.insert(pixels[k], k);
    fitTrends(model, cache);
    const bool empirical = prepareTemplates(model, response, pixels, x, index, cache, 8);
    std::array<double, 3> squared{};
    double norm = 0.;
    int scored = 0, templateInterpolations = 0, templateFallbacks = 0;
    QJsonArray folds;
    for (int held = 0; held < x.size(); ++held) {
        cancel(stop);
        QJsonObject row{{"x_mm", x[held]}};
        try {
            auto fold = train(held);
            fitTrends(fold, cache);
            const bool hasTemplate =
                empirical && prepareTemplates(fold, response, pixels, x, index, cache, 7);
            if (empirical && !hasTemplate)
                ++templateFallbacks;
            if (hasTemplate && x[held] > fold.knots.first() && x[held] < fold.knots.last())
                ++templateInterpolations;
            const double power = response.col(held).squaredNorm();
            for (int kind = 0; kind <= (empirical ? 2 : 1); ++kind) {
                fold.kind = kind == 2 && !hasTemplate ? ImagingPsf::Constant : ImagingPsf::Kind(kind);
                const ImageMatrix prediction = evaluatePsf(fold, pixels, {x[held]});
                const double error =
                    (prediction.col(0).array() + fold.calibrationOffset - response.col(held).array())
                        .square()
                        .sum();
                squared[kind] += error;
                row[modelNames[kind]] = std::sqrt(error / pixels.size()); // per-image pixel RMSE
            }
            norm += power;
            ++scored;
        } catch (const std::exception &e) {
            cancel(stop);
            row["error"] = QString::fromUtf8(e.what());
        }
        folds.append(row);
    }
    QJsonObject scores;
    if (scored == x.size() && norm > 0.) {
        for (int kind = 0; kind <= (empirical && templateInterpolations >= 3 ? 2 : 1); ++kind) {
            scores[modelNames[kind]] = std::sqrt(squared[kind] / norm);
            if (squared[kind] < squared[model.kind] - 1e-14 * norm)
                model.kind = ImagingPsf::Kind(kind);
        }
    }
    model.selection["status"] = scores.isEmpty() ? "incomplete_folds" : "compared";
    model.selection["relative_l2"] = scores;
    model.selection["scored_positions"] = scored;
    model.selection["template_status"] = model.observations.size() < 8 ? "too_few_positions"
                                         : !empirical                  ? "low_asymmetry"
                                         : scored != x.size()          ? "incomplete_folds"
                                         : templateInterpolations < 3  ? "insufficient_interpolation_tests"
                                                                       : "compared";
    model.selection["template_interior_predictions"] = templateInterpolations;
    model.selection["template_baseline_folds"] = templateFallbacks;
    model.selection["fold_rmse"] = folds;
    if (model.kind != ImagingPsf::Empirical) {
        model.knots.clear();
        model.templateResiduals.resize(0, 0);
        model.pixelIndex.clear();
    }
    return model;
}
ImageMatrix Imaging::evaluatePsf(const ImagingPsf &m, const QVector<QPoint> &pixels,
                                 const QVector<double> &x) {
    ImageMatrix response(pixels.size(), x.size());
    Parameters p;
    p << m.amplitude, 0, 0, m.sigmaRow, m.sigmaColumn, m.correlation, 0;
    for (int j = 0; j < x.size(); ++j) {
        p[1] = m.rowSlope * x[j] + m.rowIntercept;
        p[2] = m.columnSlope * x[j] + m.columnIntercept;
        if (m.kind == ImagingPsf::Varying) {
            const double t = std::clamp(x[j], m.xMin, m.xMax) - m.trendCenter;
            for (int i = 0; i < 3; ++i) {
                const auto &trend = m.trends[i];
                p[shapeIndices[i]] = std::exp(std::clamp(trend[0] + t * trend[1], trend[2], trend[3]));
            }
        }
        for (int k = 0; k < pixels.size(); ++k)
            response(k, j) = gaussian(pixels[k].y(), pixels[k].x(), p);
        if (m.kind != ImagingPsf::Empirical)
            continue;
        const auto &knots = m.knots;
        const int hi = std::clamp(int(std::lower_bound(knots.begin(), knots.end(), x[j]) - knots.begin()), 1,
                                  int(knots.size() - 1));
        const double step = knots[hi] - knots[hi - 1];
        const double t = (x[j] - knots[hi - 1]) / step;
        // Fade empirical corrections over one knot interval beyond the observed support.
        for (int n : {hi - 1, hi}) {
            const double w = std::max(0., 1. - std::abs(t - (n == hi ? 1. : 0.)));
            if (w == 0.)
                continue;
            const double dr = (x[j] - knots[n]) * m.rowSlope, dc = (x[j] - knots[n]) * m.columnSlope;
            for (int k = 0; k < pixels.size(); ++k)
                response(k, j) += w * interpolate(m.pixelIndex, pixels[k].y() - dr, pixels[k].x() - dc,
                                                  [&](int idx) { return m.templateResiduals(idx, n); });
        }
        response.col(j) = response.col(j).cwiseMax(0.);
    }
    return response;
}
QJsonObject ImagingPsf::json(double pitch) const {
    double speed = 1 / sampleMmPerPixel;
    auto diagnostic = selection;
    diagnostic["selected_model"] = modelNames[kind];
    diagnostic["fwhm_reference"] = "baseline_gaussian";
    diagnostic["calibration_offset"] = calibrationOffset;
    if (kind == Varying) {
        QJsonArray values;
        for (const auto &v : trends)
            values.append(QJsonArray{v[0], v[1], v[2], v[3]});
        diagnostic["log_shape_trends"] = values;
        diagnostic["trend_center_x_mm"] = trendCenter;
    } else if (kind == Empirical)
        diagnostic["template_x_mm"] = ImageUtil::json(knots);
    return {{"model_selection", diagnostic},
            {"row_line", QJsonArray{rowSlope, rowIntercept}},
            {"column_line", QJsonArray{columnSlope, columnIntercept}},
            {"sigma_row", sigmaRow},
            {"sigma_col", sigmaColumn},
            {"correlation", correlation},
            {"amplitude", amplitude},
            {"fitted_x_min", xMin},
            {"fitted_x_max", xMax},
            {"detector_pixels_per_sample_mm", speed},
            {"detector_pitch_mm", pitch},
            {"sample_mm_per_detector_pixel", sampleMmPerPixel},
            {"detector_mm_per_sample_mm", pitch * speed},
            {"projected_psf_fwhm_sample_mm", fwhmMm},
            {"centroid_fit_rms_detector_pixels", trajectoryRms},
            {"valid_observations", observations.size()},
            {"fit_r2_min", minR2},
            {"observations", observations}};
}
