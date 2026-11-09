#pragma once
#include "Imaging.h"
#include "Isosurface.h"
#include <QImage>
#include <map>
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

struct BondContrast {
    QVector<QVector<double>> axes; // X, Y, Z in mm
    QVector<double> values;
    QJsonObject diagnostics;
};

struct BondProjection {
    QJsonObject config;
    ImagingValidation validation;
    QJsonObject sourceRevision;
    std::map<QString, ImageMatrix> differenceByCrystal; // detector pixel by Y/Z
    double ea = 0., eb = 0.;
    QJsonObject spectralProcessing;
};

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
using ContrastSpectra = std::function<QVector<BondSpectralSettings>(const QVector<double> &, const QStringList &)>;
ImagingValidation validate(const QJsonObject &, const std::atomic_bool &);
BondVolume run(const QJsonObject &, const std::atomic_bool &, Imaging::Progress = {});
BondProjection prepareContrast(const QJsonObject &, double ea, double eb, const ContrastSpectra &,
                               const std::atomic_bool &, Imaging::Progress = {});
BondContrast reconstructContrast(const BondProjection &, const QJsonObject &spatialSettings,
                                 const std::atomic_bool &, Imaging::Progress = {});
BondReferenceSlice fitReferenceSlice(const BondVolume &, int z, double analyzerKev,
                                     const BondReferenceSpectrum &matrix, const BondReferenceSpectrum &target,
                                     std::optional<QPair<double, double>> window, double threshold,
                                     const std::atomic_bool &);
QVector<double> contrastWeights(const QVector<double> &energy, double ea, double eb,
                               const BondSpectralSettings &settings);
void saveContrast(const QString &, const BondContrast &, const std::atomic_bool &);
std::optional<BondContrast> loadContrast(const QString &);
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
