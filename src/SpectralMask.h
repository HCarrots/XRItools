#pragma once
#include <QPoint>
#include <QSize>
#include <QString>
#include <QVector>
#include <atomic>
#include <stdexcept>

class NoSpectralComponent : public std::runtime_error {
  public:
    explicit NoSpectralComponent(const QString &reason) : std::runtime_error(reason.toStdString()) {}
};

enum class UnmixMethod { NNMF, PCA, ICA };
struct PixelSpectra {
    QSize detectorSize;
    QVector<QPoint> pixels;
    QVector<double> energyKev;
    QVector<QVector<double>> intensity; // [pixel][incident energy], normalized once if requested
};
struct SpectralComponents {
    QVector<QVector<double>> spectra, covariance; // one curve and one score per candidate pixel
    double relativeResidual = 0;
    bool converged = true;
};
SpectralComponents decomposePixels(const PixelSpectra &data, const QVector<int> &rows, UnmixMethod method,
                                   int count, const std::atomic_bool &cancel);
