#include "SpectralMask.h"
#include <Eigen/Dense>
#include <QCoreApplication>
#include <algorithm>
#include <numeric>
#include <random>
#include <stdexcept>

namespace {
using Matrix = Eigen::MatrixXd;
using Vector = Eigen::VectorXd;
[[noreturn]] void fail(const QString &message) {
    throw std::runtime_error(message.toUtf8().constData());
}
void checkCancel(const std::atomic_bool &cancel) {
    if (cancel.load())
        fail(QCoreApplication::translate("SpectralMask", "已取消。"));
}

// Alternating non-negative least squares using hierarchical coordinate updates (HALS).
void updateFactors(Matrix &factor, const Matrix &cross, const Matrix &gram) {
    for (int k = 0; k < factor.cols(); ++k) {
        if (gram(k, k) > 1e-20) {
            const Vector delta = (cross.col(k) - factor * gram.col(k)) / gram(k, k);
            factor.col(k) = (factor.col(k) + delta).cwiseMax(0.0);
        }
    }
}
Vector randomVector(int n, std::mt19937 &rng) {
    std::uniform_real_distribution<double> uniform(0.1, 1.0);
    return Vector::NullaryExpr(n, [&] { return uniform(rng); });
}
} // namespace

SpectralComponents decomposePixels(const PixelSpectra &data, const QVector<int> &rows, UnmixMethod method,
                                   int count, const std::atomic_bool &cancel) {
    const int n = int(rows.size()), m = int(data.energyKev.size());
    if (count < 1 || count > std::min(n, m) || m < 2)
        throw NoSpectralComponent(QCoreApplication::translate("SpectralMask", "组分数须在有效像素数和能量点数允许的范围内。"));
    checkCancel(cancel);
    Matrix x(n, m);
    for (int p = 0; p < n; ++p) {
        if (rows[p] < 0 || rows[p] >= data.intensity.size() || data.intensity[rows[p]].size() != m)
            fail(QCoreApplication::translate("SpectralMask", "像素谱矩阵维度不一致。"));
        x.row(p) = Eigen::Map<const Eigen::RowVectorXd>(data.intensity[rows[p]].constData(), m);
    }
    if (!x.allFinite() || (method == UnmixMethod::NNMF && x.minCoeff() < 0))
        fail(QCoreApplication::translate("SpectralMask", "谱数据须为有限值，NNMF 还要求强度非负。"));
    const double scale = x.cwiseAbs().maxCoeff();
    if (scale <= 0)
        throw NoSpectralComponent(QCoreApplication::translate("SpectralMask", "候选像素在参考扫描中没有信号。"));
    x /= scale;
    Matrix weights(n, count), shapes(m, count);
    SpectralComponents result;
    std::mt19937 rng(17);
    if (method == UnmixMethod::NNMF) {
        for (int k = 0; k < count; ++k) {
            weights.col(k) = randomVector(n, rng);
            shapes.col(k) = randomVector(m, rng).normalized();
        }
        double previous = 1e100;
        result.converged = false;
        for (int iteration = 0; iteration < 600; ++iteration) {
            checkCancel(cancel);
            updateFactors(weights, x * shapes, shapes.transpose() * shapes);
            updateFactors(shapes, x.transpose() * weights, weights.transpose() * weights);
            if (iteration % 5 == 0) {
                const double error = (x - weights * shapes.transpose()).norm() / x.norm();
                if (error < 1e-7 || std::abs(previous - error) < 1e-7 * std::max(error, 1e-4)) {
                    result.converged = true;
                    break;
                }
                previous = error;
            }
        }
    } else {
        // Observations are pixels. Center each energy channel across pixels.
        const Eigen::RowVectorXd mean = x.colwise().mean();
        x.rowwise() -= mean;
        Eigen::JacobiSVD<Matrix> svd(x, Eigen::ComputeThinU | Eigen::ComputeThinV);
        if (count > svd.rank())
            throw NoSpectralComponent(QCoreApplication::translate("SpectralMask", "中心化后可分辨的组分不足，请减少组分数。"));
        weights = svd.matrixU().leftCols(count) * std::sqrt(double(n));
        shapes = svd.matrixV().leftCols(count) * svd.singularValues().head(count).asDiagonal() /
                 std::sqrt(double(n));
        if (method == UnmixMethod::ICA) {
            // Deflationary FastICA on the whitened pixel observations.
            Matrix directions = Matrix::Zero(count, count);
            for (int k = 0; k < count; ++k) {
                Vector w = randomVector(count, rng).normalized();
                bool converged = false;
                for (int iteration = 0; iteration < 400; ++iteration) {
                    checkCancel(cancel);
                    Vector g = (weights * w).array().tanh();
                    Vector next = weights.transpose() * g / n - (1 - g.array().square()).mean() * w;
                    if (k)
                        next -= (directions.leftCols(k) * (directions.leftCols(k).transpose() * next)).eval();
                    if (next.norm() < 1e-12)
                        throw NoSpectralComponent(QCoreApplication::translate("SpectralMask",
                                                         "ICA 分离方向退化，请减少组分数或选择 PCA。"));
                    next.normalize();
                    converged = 1 - std::abs(next.dot(w)) < 1e-7;
                    w = next;
                    if (converged)
                        break;
                }
                result.converged &= converged;
                directions.col(k) = w;
            }
            weights = (weights * directions).eval();
            shapes = (shapes * directions).eval();
        }
    }
    result.relativeResidual = (x - weights * shapes.transpose()).norm() / x.norm();
    QVector<int> order(count);
    std::iota(order.begin(), order.end(), 0);
    if (method == UnmixMethod::NNMF)
        std::stable_sort(order.begin(), order.end(), [&](int a, int b) {
            return weights.col(a).sum() * shapes.col(a).sum() > weights.col(b).sum() * shapes.col(b).sum();
        });
    for (int k : order) {
        Vector h = shapes.col(k);
        if (h.norm() > 0)
            h.normalize();
        Eigen::Index peak;
        h.cwiseAbs().maxCoeff(&peak);
        if (h[peak] < 0)
            h = -h;
        result.spectra.append(QVector<double>(h.data(), h.data() + m));
        const Vector centered = h.array() - h.mean();
        QVector<double> scores(n);
        for (int p = 0; p < n; ++p) {
            const Eigen::Map<const Vector> measured(data.intensity[rows[p]].constData(), m);
            scores[p] = (measured.array() - measured.mean()).matrix().dot(centered) / (m - 1);
        }
        result.covariance.append(scores);
    }
    return result;
}
