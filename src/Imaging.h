#pragma once
#include "HdfFile.h"
#include <Eigen/Core>
#include <QHash>
#include <QJsonArray>
#include <QJsonObject>
#include <QPoint>
#include <QRect>
#include <array>
#include <functional>

using ImageMatrix = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
struct ImagingRoi {
    QString id, label, sourceLabel, key, detector;
    QRect rect;
    bool excluded = false;
    QVector<QPoint> pixels;
    bool sourceInvalid = false;
};
struct ImagingGrid {
    QString file;
    QStringList loopOrder;
    QVector<QVector<double>> axes, coordinates;
    QVector<QVector<int>> indices;
    QVector<int> occupancy;
    QVector<double> monitor, dose;
    QSize imageSize;
    QString plan;
    QJsonObject preparedFoil;
};
struct ImagingEnergyBand {
    QVector<double> energies, sampleWeights, calibrationWeights, blankFrameWeights;
    QVector<double> verificationEnergies, verificationWeights;
    QString notice;
};
struct ImagingValidation {
    ImagingGrid calibration, sample;
    ImagingGrid verification;
    ImagingEnergyBand energyBand;
    QVector<ImagingRoi> rois;
    QJsonObject report;
};
struct ImagingPsf {
    enum Kind { Constant, Varying, Empirical };
    Kind kind = Constant;
    double calibrationOffset = 0., trendCenter = 0.;
    // Log amplitude, row width and column width: intercept, slope, observed lower/upper bounds.
    std::array<std::array<double, 4>, 3> trends{};
    QVector<double> knots;
    ImageMatrix templateResiduals;
    QHash<QPoint, int> pixelIndex;
    QJsonObject selection;
    double rowSlope = 0, rowIntercept = 0, columnSlope = 0, columnIntercept = 0;
    double sigmaRow = 0, sigmaColumn = 0, correlation = 0, amplitude = 0;
    double sampleMmPerPixel = 0, trajectoryRms = 0, fwhmMm = 0, xMin = 0, xMax = 0, minR2 = 0;
    QJsonArray observations;
    QJsonObject json(double pitch) const;
};
// Independent fit to one measured foil image. Trajectory fitting and quality cuts remain fold-specific.
struct ImagingSpotFit {
    std::array<double, 7> parameters{};
    double rSquared = 0.;
    bool fitted = false;
    double asymmetry = -1.;
};
struct ImagingSolverInfo {
    bool converged = false;
    int iterations = 0;
    double kkt = 0, primal = 0, dual = 0;
};
struct ImagingResult {
    QVector<QVector<double>> axes;
    ImageMatrix nnls, tv, primary;
    QVector<double> sensitivity, eigenvalues;
    QJsonObject summary, models, rejected;
    QString directory;
};
class Imaging {
  public:
    using Progress = std::function<void(const QString &, int)>;
    static QJsonObject defaults();
    static QJsonObject config(const QJsonObject &, const QString &baseDirectory);
    static QJsonObject loadConfig(const QString &file, QString *bondOutputFile = nullptr);
    static ImagingValidation validate(const QJsonObject &, const std::atomic_bool &cancel,
                                      bool selectIntegrationBand = true);
    static void validateRois(ImagingValidation &, const QJsonObject &);
    static ImagingResult run(const QJsonObject &, const std::atomic_bool &, Progress = {});
    static QJsonObject validateDepth(const QString &resultDirectory, const std::atomic_bool &, Progress = {});
    static QJsonObject comparePsfModels(const QString &resultDirectory, const std::atomic_bool &, Progress = {});
    static QJsonObject validateBoundaries(const QString &resultDirectory, const std::atomic_bool &, Progress = {});
    static QVector<double> energyWeights(const QVector<double> &);
    static QPair<QVector<double>, QVector<double>> matchedWeights(const QVector<double> &source,
                                                                  const QVector<double> &target,
                                                                  const QString &referenceLabel = {});
    static ImagingPsf fitPsf(const ImageMatrix &, const QVector<QPoint> &, const QVector<double> &,
                             const QJsonObject &, const std::atomic_bool &,
                             QHash<double, ImagingSpotFit> *spotCache = nullptr);
    static ImageMatrix evaluatePsf(const ImagingPsf &, const QVector<QPoint> &, const QVector<double> &);
    static ImageMatrix nnls(const ImageMatrix &h, const ImageMatrix &b, int maxIterations, double tolerance,
                            const std::atomic_bool &, ImagingSolverInfo *, const ImageMatrix &initial = {});
    static std::array<ImageMatrix, 3> gradient(const ImageMatrix &, int ny, int nz,
                                               const std::array<double, 3> &spacing);
    static ImageMatrix adjoint(const std::array<ImageMatrix, 3> &, int ny, int nz,
                               const std::array<double, 3> &spacing);
    static ImageMatrix tv(const ImageMatrix &, const ImageMatrix &, int ny, int nz,
                          const std::array<double, 3> &spacing, double lambda, int maxIterations,
                          double tolerance, const ImageMatrix &initial, const std::atomic_bool &,
                          ImagingSolverInfo *, Progress = {});
};
