#pragma once
#include "Spectroscopy.h"
#include <QJsonObject>

struct SpectrumGeometry {
    double qMean = NAN, analyzerKev = NAN, angleDeg = NAN;
};
struct ProcessedSpectrum {
    QVector<double> energyEv;
    QMap<QString, QVector<double>> crystals;
    QString intensityUnit;
    QStringList sources;
    QMap<QString, SpectrumGeometry> geometry;
    QMap<QString, QVector<int>> pointSources;
    QMap<QString, QVector<double>> wideCrystals;
    QJsonObject processing;
};
struct SpectrumHdfLayout {
    QString energy = "/entry/spectrum/energy_transfer_ev";
    QString labels = "/entry/spectrum/crystal_labels";
    QString intensity = "/entry/spectrum/intensity";
    QString externalIntensityUnit = "a.u.";
    bool energyInKev = false, energyFirst = false;
};
struct BackgroundSettings {
    bool pearson = false;
    QVector<QPair<double, double>> windows;
};
struct CrystalBackground {
    QVector<double> background, corrected;
    BackgroundSettings settings;
    double rSquared = NAN, centerEv = NAN, fwhmEv = NAN, exponent = NAN;
    double peakHeight = NAN;
    QString method;
};

namespace SpectrumProcessing {
ProcessedSpectrum fromResult(const SpectrumResult &, const CrystalRois &, bool normalized);
ProcessedSpectrum load(const QString &, const SpectrumHdfLayout &);
QMap<QString, SpectrumGeometry> loadGeometry(const QString &, const ProcessedSpectrum &);
ProcessedSpectrum replaceEdges(const ProcessedSpectrum &, const QVector<ProcessedSpectrum> &,
                               const std::atomic_bool &);
ProcessedSpectrum stitch(const QStringList &, const SpectrumHdfLayout &, const std::atomic_bool &);
CrystalBackground fit(const QVector<double> &energy, const QVector<double> &intensity,
                      const BackgroundSettings &, const std::atomic_bool &);
void save(const QString &path, const ProcessedSpectrum &,
          const QMap<QString, CrystalBackground> &backgrounds = {});
} // namespace SpectrumProcessing
