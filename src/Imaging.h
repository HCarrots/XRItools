#pragma once
#include "HdfFile.h"
#include <Eigen/Core>
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
};
struct ImagingGrid {
    QString file;
    QVector<QVector<double>> axes, coordinates;
    QVector<QVector<int>> indices;
    QVector<int> occupancy;
    QVector<double> monitor, dose;
    QSize imageSize;
    QString plan;
};
struct ImagingValidation {
    ImagingGrid calibration, sample;
    QVector<ImagingRoi> rois;
    QJsonObject report;
};
struct ImagingPsf {
    double rowSlope = 0, rowIntercept = 0, columnSlope = 0, columnIntercept = 0;
    double sigmaRow = 0, sigmaColumn = 0, correlation = 0, amplitude = 0;
    double sampleMmPerPixel = 0, trajectoryRms = 0, fwhmMm = 0, xMin = 0, xMax = 0, minR2 = 0;
    QJsonArray observations;
    QJsonObject json(double pitch) const;
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
    static QJsonObject loadConfig(const QString &file);
    static ImagingValidation validate(const QJsonObject &, const std::atomic_bool &cancel);
    static ImagingResult run(const QJsonObject &, const std::atomic_bool &, Progress = {});
    static QVector<double> energyWeights(const QVector<double> &);
    static QPair<QVector<double>, QVector<double>> matchedWeights(const QVector<double> &,
                                                                  const QVector<double> &);
    static ImagingPsf fitPsf(const ImageMatrix &, const QVector<QPoint> &, const QVector<double> &,
                             const QJsonObject &, const std::atomic_bool &);
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
