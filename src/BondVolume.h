#pragma once
#include "Imaging.h"
#include <QImage>
#include <QVector3D>
#include <optional>

struct BondBandDiagnostic {
    QString crystal;
    ImageMatrix forward, predicted, residual;
};

struct BondVolume {
    QVector<QVector<double>> axes; // incident keV, X/Y/Z mm
    QVector<double> intensity;     // E, X, Y, Z
    QVector<double> sensitivity;
    QVector<bool> converged;
    QJsonObject diagnostics;
    QVector<BondBandDiagnostic> bandDiagnostics;
    QImage diagnosticImage;
    qsizetype voxelCount() const;
    qsizetype voxel(int x, int y, int z) const;
};

struct BondSpectralSettings {
    QVector<QPair<double, double>> background;
    std::optional<QPair<double, double>> normalization;
};

using BondTriangle = std::array<QVector3D, 3>;

namespace BondImaging {
ImagingValidation validate(const QJsonObject &, const std::atomic_bool &);
BondVolume run(const QJsonObject &, const std::atomic_bool &, Imaging::Progress = {});
void save(const QString &, const BondVolume &, const std::atomic_bool &);
QString diagnosticPngPath(const QString &);
BondVolume load(const QString &);
QVector<double> process(const BondVolume &, const QVector<double> &energy, const BondSpectralSettings &,
                        const std::atomic_bool &);
QVector<double> energySlice(const BondVolume &, const QVector<double> &values, const QVector<double> &energy,
                            double at, bool requireConvergence);
QVector<BondTriangle> surface(const BondVolume &, const QVector<double> &field, double level,
                              const std::atomic_bool &);
} // namespace BondImaging
