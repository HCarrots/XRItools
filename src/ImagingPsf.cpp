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
} // namespace
ImagingPsf Imaging::fitPsf(const ImageMatrix &response, const QVector<QPoint> &pixels,
                           const QVector<double> &x, const QJsonObject &settings,
                           const std::atomic_bool &stop) {
    if (response.rows() != pixels.size() || response.cols() != x.size() || !response.allFinite() ||
        response.minCoeff() < 0 || pixels.isEmpty())
        fail("Expected finite nonnegative PSF response[p,x]");
    QVector<double> height;
    for (int j = 0; j < x.size(); ++j) {
        QVector<double> values;
        for (int p = 0; p < pixels.size(); ++p)
            values.append(response(p, j));
        height.append(percentile(values, 99.5) - percentile(values, 20));
    }
    double threshold =
        *std::max_element(height.begin(), height.end()) * settings["min_peak_fraction"].toDouble();
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
        if (height[j] <= threshold)
            continue;
        QVector<double> image;
        for (int p = 0; p < pixels.size(); ++p)
            image.append(response(p, j));
        Parameters fitted;
        double r2;
        if (!fit(pixels, image, settings, &fitted, &r2, stop))
            continue;
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
    QVector<double> amps, rows, cols, corr;
    for (int i = 0; i < fits.size(); ++i) {
        auto &p = fits[i];
        double dr = m.rowSlope * validX[i] + m.rowIntercept - p[1],
               dc = m.columnSlope * validX[i] + m.columnIntercept - p[2];
        squared += dr * dr + dc * dc;
        amps.append(p[0]);
        rows.append(p[3]);
        cols.append(p[4]);
        corr.append(p[5]);
        m.observations.append(
            QJsonObject{{"x", validX[i]}, {"row", p[1]}, {"column", p[2]}, {"r_squared", r2s[i]}});
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
ImageMatrix Imaging::evaluatePsf(const ImagingPsf &m, const QVector<QPoint> &pixels,
                                 const QVector<double> &x) {
    ImageMatrix response(pixels.size(), x.size());
    Parameters p;
    p << m.amplitude, 0, 0, m.sigmaRow, m.sigmaColumn, m.correlation, 0;
    for (int j = 0; j < x.size(); ++j) {
        p[1] = m.rowSlope * x[j] + m.rowIntercept;
        p[2] = m.columnSlope * x[j] + m.columnIntercept;
        for (int k = 0; k < pixels.size(); ++k)
            response(k, j) = gaussian(pixels[k].y(), pixels[k].x(), p);
    }
    return response;
}
QJsonObject ImagingPsf::json(double pitch) const {
    double speed = 1 / sampleMmPerPixel;
    return {{"row_line", QJsonArray{rowSlope, rowIntercept}},
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
