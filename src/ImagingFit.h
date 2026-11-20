#pragma once
#include "Imaging.h"

struct FitObservation {
    QString key;
    const ImageMatrix *response, *data, *variance;
    double savedVarianceFloor = 0.;
};
struct OperatorFit {
    ImageMatrix h, weightedH, b, rho;
    QVector<double> gains, history, varianceFloors;
    ImagingSolverInfo nnls;
    bool gainConverged = false;
    const ImageMatrix &normal() const {
        return weightedH.size() ? weightedH : h;
    }
};
OperatorFit fitOperators(const QVector<FitObservation> &, const QJsonObject &reconstruction,
                         const QJsonObject &statistics, const std::atomic_bool &, Imaging::Progress = {});
QJsonObject compareCrystalGroups(const QVector<FitObservation> &, const QVector<QVector<double>> &axes,
                                 const QJsonObject &config, double referenceLambda, const std::atomic_bool &,
                                 Imaging::Progress = {});
