#include "ImagingIO.h"
#include <Eigen/Eigenvalues>
#include <algorithm>
#include <cmath>
using namespace ImageUtil;
namespace {
struct Resolved {
    ImageMatrix matrix, basis, target;
    Eigen::VectorXd spectrum;
    double lipschitz;
    bool reduced;
    Resolved(const ImageMatrix &h, const ImageMatrix &b) {
        if (h.rows() != h.cols() || h.rows() != b.rows() || !h.allFinite() || !b.allFinite())
            fail("Invalid finite NNLS sufficient statistics");
        Eigen::SelfAdjointEigenSolver<ImageMatrix> eigen(h);
        if (eigen.info() != Eigen::Success)
            fail("Gram eigensolver failed");
        lipschitz = eigen.eigenvalues().maxCoeff();
        if (!(lipschitz > 0))
            fail("Zero response operator");
        int count = 0;
        for (auto v : eigen.eigenvalues())
            if (v > lipschitz * 1e-14)
                ++count;
        basis.resize(h.rows(), count);
        spectrum.resize(count);
        int k = 0;
        for (int i = 0; i < h.rows(); ++i)
            if (eigen.eigenvalues()[i] > lipschitz * 1e-14) {
                basis.col(k) = eigen.eigenvectors().col(i);
                spectrum[k++] = eigen.eigenvalues()[i] / lipschitz;
            }
        matrix.noalias() = basis * spectrum.asDiagonal() * basis.transpose();
        target.noalias() = basis * (basis.transpose() * b);
        target /= lipschitz;
        reduced = count * 2 < h.rows();
    }
    ImageMatrix product(const ImageMatrix &value) const {
        if (reduced) {
            ImageMatrix small = basis.transpose() * value;
            small.array().colwise() *= spectrum.array();
            return basis * small;
        }
        return matrix * value;
    }
};
} // namespace
ImageMatrix Imaging::nnls(const ImageMatrix &h, const ImageMatrix &b, int maxIterations, double tolerance,
                          const std::atomic_bool &stop, ImagingSolverInfo *info, const ImageMatrix &initial) {
    if (maxIterations < 1 || !std::isfinite(tolerance) || tolerance <= 0)
        fail("Invalid NNLS stopping parameters");
    Resolved terms(h, b);
    ImageMatrix current = initial.size() ? initial : ImageMatrix::Zero(b.rows(), b.cols());
    if (current.rows() != b.rows() || current.cols() != b.cols() || !current.allFinite() ||
        current.minCoeff() < 0)
        fail("Invalid initial NNLS density");
    ImageMatrix extrapolated = current;
    double acceleration = 1;
    *info = {};
    double targetNorm = terms.target.norm();
    for (int iteration = 1; iteration <= maxIterations; ++iteration) {
        if (iteration % 20 == 1)
            cancel(stop);
        ImageMatrix updated = (extrapolated - terms.product(extrapolated) + terms.target).cwiseMax(0.);
        info->iterations = iteration;
        if (iteration % 20 == 0 || iteration == 1 || iteration == maxIterations) {
            ImageMatrix g = terms.product(updated) - terms.target;
            ImageMatrix projected = updated - (updated - g).cwiseMax(0.);
            info->kkt = projected.norm() / std::max({updated.norm(), targetNorm, 1e-15});
            if (info->kkt < tolerance) {
                current = updated;
                info->converged = true;
                break;
            }
        }
        double next = (1 + std::sqrt(1 + 4 * acceleration * acceleration)) / 2;
        ImageMatrix nextExtrapolated = updated + ((acceleration - 1) / next) * (updated - current);
        if (((extrapolated - updated).array() * (updated - current).array()).sum() > 0) {
            next = 1;
            nextExtrapolated = updated;
        }
        current = std::move(updated);
        extrapolated = std::move(nextExtrapolated);
        acceleration = next;
    }
    return current;
}
std::array<ImageMatrix, 3> Imaging::gradient(const ImageMatrix &rho, int ny, int nz,
                                             const std::array<double, 3> &spacing) {
    if (rho.cols() != ny * nz)
        fail("Volume shape does not match gradient axes");
    std::array<ImageMatrix, 3> g;
    for (auto &m : g)
        m = ImageMatrix::Zero(rho.rows(), rho.cols());
    for (int x = 0; x < rho.rows(); ++x)
        for (int y = 0; y < ny; ++y)
            for (int z = 0; z < nz; ++z) {
                int c = y * nz + z;
                if (x + 1 < rho.rows())
                    g[0](x, c) = (rho(x + 1, c) - rho(x, c)) / spacing[0];
                if (y + 1 < ny)
                    g[1](x, c) = (rho(x, c + nz) - rho(x, c)) / spacing[1];
                if (z + 1 < nz)
                    g[2](x, c) = (rho(x, c + 1) - rho(x, c)) / spacing[2];
            }
    return g;
}
ImageMatrix Imaging::adjoint(const std::array<ImageMatrix, 3> &g, int ny, int nz,
                             const std::array<double, 3> &spacing) {
    ImageMatrix out = ImageMatrix::Zero(g[0].rows(), g[0].cols());
    for (int x = 0; x < out.rows(); ++x)
        for (int y = 0; y < ny; ++y)
            for (int z = 0; z < nz; ++z) {
                int c = y * nz + z;
                if (x + 1 < out.rows()) {
                    double v = g[0](x, c) / spacing[0];
                    out(x, c) -= v;
                    out(x + 1, c) += v;
                }
                if (y + 1 < ny) {
                    double v = g[1](x, c) / spacing[1];
                    out(x, c) -= v;
                    out(x, c + nz) += v;
                }
                if (z + 1 < nz) {
                    double v = g[2](x, c) / spacing[2];
                    out(x, c) -= v;
                    out(x, c + 1) += v;
                }
            }
    return out;
}
ImageMatrix Imaging::tv(const ImageMatrix &h, const ImageMatrix &b, int ny, int nz,
                        const std::array<double, 3> &spacing, double lambda, int maxIterations,
                        double tolerance, const ImageMatrix &initial, const std::atomic_bool &stop,
                        ImagingSolverInfo *info, Progress progress) {
    if (lambda < 0 || !std::isfinite(lambda) || maxIterations < 1 || tolerance <= 0 || b.cols() != ny * nz)
        fail("Invalid TV parameters");
    for (auto step : spacing)
        if (!std::isfinite(step) || step <= 0)
            fail("TV requires positive physical spacing");
    Resolved terms(h, b);
    double lam = lambda / terms.lipschitz, norm = 0;
    int dims[] = {int(b.rows()), ny, nz};
    for (int axis = 0; axis < 3; ++axis)
        if (dims[axis] > 1)
            norm += 4 / (spacing[axis] * spacing[axis]);
    double tau = .9, sigma = norm > 0 ? .45 / norm : 1;
    ImageMatrix rho = initial;
    if (rho.rows() != b.rows() || rho.cols() != b.cols() || !rho.allFinite() || rho.minCoeff() < 0)
        fail("Invalid initial TV density");
    std::array<ImageMatrix, 3> dual;
    for (auto &m : dual)
        m = ImageMatrix::Zero(b.rows(), b.cols());
    *info = {};
    for (int iteration = 1; iteration <= maxIterations; ++iteration) {
        if (iteration % 20 == 1)
            cancel(stop);
        ImageMatrix next =
            (rho - tau * (terms.product(rho) - terms.target + adjoint(dual, ny, nz, spacing))).cwiseMax(0.);
        auto g = gradient(2 * next - rho, ny, nz, spacing);
        std::array<ImageMatrix, 3> newDual;
        for (int a = 0; a < 3; ++a)
            newDual[a] = dual[a] + sigma * g[a];
        if (lam == 0) {
            for (auto &m : newDual)
                m.setZero();
        } else
            for (Eigen::Index i = 0; i < rho.size(); ++i) {
                double length = std::sqrt(newDual[0].data()[i] * newDual[0].data()[i] +
                                          newDual[1].data()[i] * newDual[1].data()[i] +
                                          newDual[2].data()[i] * newDual[2].data()[i]),
                       scale = std::max(1., length / lam);
                for (auto &m : newDual)
                    m.data()[i] /= scale;
            }
        info->primal = (next - rho).norm() / std::max(next.norm(), 1e-15);
        double delta = 0, size = 0;
        for (int a = 0; a < 3; ++a) {
            delta += (newDual[a] - dual[a]).squaredNorm();
            size += newDual[a].squaredNorm();
        }
        info->dual = std::sqrt(delta) / std::max(std::sqrt(size), 1e-15);
        rho = std::move(next);
        dual = std::move(newDual);
        info->iterations = iteration;
        if (progress && iteration % 1000 == 0)
            progress(QString("TV %1 · primal=%2 · dual=%3")
                         .arg(iteration)
                         .arg(info->primal, 0, 'g', 5)
                         .arg(info->dual, 0, 'g', 5),
                     75 + iteration * 20 / maxIterations);
        if (iteration > 10 && std::max(info->primal, info->dual) < tolerance) {
            info->converged = true;
            break;
        }
    }
    return rho;
}
