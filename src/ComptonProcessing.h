#pragma once
#include "SpectrumProcessing.h"

using SpectrumWindows = QVector<QPair<double, double>>;
struct ComptonSettings {
    double highQMin = 9.5, coreScale = 1., valenceCutoffEv = 20.;
    QString hfTable, composition;
    SpectrumWindows normalization, coreTail, cleanReference, highResidual, lowLinear;
};
struct ComptonCrystal {
    QVector<double> input, hfCore, valence;
    CrystalBackground fit;
    bool high = false;
    double scale = 1.;
};
struct ComptonResult {
    QMap<QString, ComptonCrystal> crystals;
    QMap<QString, QString> failures, referenceFailures;
    QStringList references;
    QVector<double> pz, measuredJ, fittedJ;
    CrystalBackground referenceFit;
    QJsonObject parameters;
};
struct ComptonTransfer {
    QVector<double> correction, reference, valence;
    double analyzerKev = NAN;
};

ComptonResult analyzeCompton(const ProcessedSpectrum &wide, const ProcessedSpectrum &stitched,
                             const ComptonSettings &, const std::atomic_bool &stop);
ComptonResult analyzeLinear(const ProcessedSpectrum &, const SpectrumWindows &, const std::atomic_bool &stop);
ComptonTransfer transferCompton(const ProcessedSpectrum &wide, const QString &label,
                                const QVector<double> &incidentKev, const ComptonSettings &,
                                const std::atomic_bool &stop);
QMap<QString, ComptonTransfer> transferCompton(const ProcessedSpectrum &wide, const QStringList &labels,
                                               const QVector<double> &incidentKev, const ComptonSettings &,
                                               const std::atomic_bool &stop);
