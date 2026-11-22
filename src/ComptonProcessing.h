#pragma once
#include "SpectrumProcessing.h"
#include <array>

using SpectrumWindows = QVector<QPair<double, double>>;
enum class BackgroundComponent { Polynomial, HfCore, ComptonPeak, Asymmetry, Count };
inline constexpr int backgroundComponentCount = int(BackgroundComponent::Count);
struct StagedBackgroundSettings {
    std::array<SpectrumWindows, backgroundComponentCount> windows;
    SpectrumWindows protectedWindows;
    QString hfTable, composition;
    QStringList hfShells;
    int polynomialOrder = 1;
    double asymmetryRise = 1., asymmetryDecay = 1.;
};
struct StagedSpectrumFit {
    QMap<BackgroundComponent, CrystalBackground> terms;
    CrystalBackground total;
};
struct StagedBackgroundResult {
    QMap<QString, StagedSpectrumFit> spectra;
    QMap<QString, QString> failures;
    QJsonObject parameters;
};
QString backgroundComponentKey(BackgroundComponent);
QString effectiveHfTable(const QString &path);
QStringList backgroundHfShellSelection(const QString &text);
struct HfShellInfo {
    double bindingEv;
    QString edge;
};
struct HfShellSelection {
    QStringList identifiers, missing;
};
QMap<QString, HfShellInfo> backgroundHfShells(const QString &table, const QString &composition);
HfShellSelection resolveHfShellSelection(const QMap<QString, HfShellInfo> &, const QStringList &);
CrystalBackground fitBackgroundComponent(BackgroundComponent, const StagedBackgroundSettings &,
                                         const QVector<double> &energy, const QVector<double> &residual,
                                         SpectrumGeometry, const std::atomic_bool &stop);
StagedSpectrumFit fitBackgroundStages(const StagedBackgroundSettings &, const QVector<double> &energy,
                                      const QVector<double> &input, SpectrumGeometry, const std::atomic_bool &stop);
CrystalBackground sumBackground(const QVector<double> &input, const QMap<BackgroundComponent, CrystalBackground> &);
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

ComptonResult analyzeCompton(const ProcessedSpectrum &wide, const ProcessedSpectrum &target,
                             const ComptonSettings &, const std::atomic_bool &stop);
ComptonResult analyzeLinear(const ProcessedSpectrum &, const SpectrumWindows &, const std::atomic_bool &stop);
ComptonTransfer transferCompton(const ProcessedSpectrum &wide, const QString &label,
                                const QVector<double> &incidentKev, const ComptonSettings &,
                                const std::atomic_bool &stop);
QMap<QString, ComptonTransfer> transferCompton(const ProcessedSpectrum &wide, const QStringList &labels,
                                               const QVector<double> &incidentKev, const ComptonSettings &,
                                               const std::atomic_bool &stop);
