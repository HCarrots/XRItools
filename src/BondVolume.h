#pragma once
#include "Imaging.h"
#include "Isosurface.h"
#include <QImage>
#include <memory>
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

struct BondSpectralReference {
    QVector<double> correction, reference, valence;
    QPair<double, double> matchWindow;
};

struct BondSpectralSettings {
    QVector<QPair<double, double>> background;
    std::optional<QPair<double, double>> normalization;
    std::optional<BondSpectralReference> compton;
};

struct BondInputs {
    QJsonObject config;
    ImagingValidation validation;
    QJsonObject sourceRevision;
};

struct BondSpatialModel;
struct BondChannelProjections;

struct BondReferenceSpectrum {
    QVector<double> energyEv, intensity;
};

struct BondReferenceSlice {
    QVector<double> total, residual, coefficient; // Y rows, X columns
    QVector<double> energyEv, matrixReference, targetReference, positiveMean, negativeMean;
    qsizetype positiveCount = 0, negativeCount = 0;
    int z = 0;
    double threshold = 0.;
};

using BondTriangle = SurfaceTriangle;

namespace BondImaging {
ImagingValidation validate(const QJsonObject &, const std::atomic_bool &);
std::shared_ptr<const BondInputs> prepareInputs(const QJsonObject &, const std::atomic_bool &);
std::shared_ptr<const BondSpatialModel> prepareSpatial(std::shared_ptr<const BondInputs>,
                                                       const QJsonObject &spatialSettings,
                                                       const std::atomic_bool &, Imaging::Progress = {});
QJsonObject spatialSummary(const BondSpatialModel &);
std::shared_ptr<const BondChannelProjections> prepareChannels(std::shared_ptr<const BondSpatialModel>,
                                                              const std::atomic_bool &, Imaging::Progress = {});
QVector<double> channelNorms(const BondChannelProjections &);
BondVolume reconstructVolume(const BondChannelProjections &, const QJsonObject &outputOptions,
                             const std::atomic_bool &, Imaging::Progress = {});
BondVolume run(const QJsonObject &, const std::atomic_bool &, Imaging::Progress = {});
BondReferenceSlice fitReferenceSlice(const BondVolume &, int z, double analyzerKev,
                                     const BondReferenceSpectrum &matrix, const BondReferenceSpectrum &target,
                                     std::optional<QPair<double, double>> window, double threshold,
                                     const std::atomic_bool &);
void save(const QString &, const BondVolume &, const std::atomic_bool &);
QString diagnosticPngPath(const QString &);
BondVolume load(const QString &);
QVector<double> process(const BondVolume &, const QVector<double> &energy, const BondSpectralSettings &,
                        const std::atomic_bool &);
QVector<double> energySlice(const BondVolume &, const QVector<double> &values, const QVector<double> &energy,
                            double at, bool requireConvergence);
QVector<double> energyIntegral(const BondVolume &, const QVector<double> &values, QPair<double, double> boundsKev);
QVector<BondTriangle> surface(const BondVolume &, const QVector<double> &field, double level,
                              const std::atomic_bool &);
} // namespace BondImaging
