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
    // Raw wide scans retain the actual incident samples shifted by each crystal's elastic center.
    QMap<QString, QVector<double>> crystalEnergyEv;
    const QVector<double> &energyForCrystal(const QString &label) const {
        const auto axis = crystalEnergyEv.constFind(label);
        return axis == crystalEnergyEv.cend() ? energyEv : axis.value();
    }
};
struct SpectrumHdfLayout {
    QString energy = "/entry/spectrum/energy_transfer_ev";
    QString labels = "/entry/spectrum/crystal_labels";
    QString intensity = "/entry/spectrum/intensity";
    QString externalIntensityUnit = "a.u.";
    bool energyInKev = false, energyFirst = false;
};
struct WideScanCalibration {
    SpectrumSettings settings;
    ElasticCalibration calibration;
    QStringList labels;
    QMap<QString, double> angles;
};
struct SpectrumSamples {
    QVector<double> energyEv, intensity;
};
struct BackgroundSettings {
    bool pearson = false;
    QVector<QPair<double, double>> windows;
    bool zeroBaseline = false;
};
struct CrystalBackground {
    QVector<double> background, corrected;
    BackgroundSettings settings;
    double rSquared = NAN, centerEv = NAN, fwhmEv = NAN, exponent = NAN;
    double peakHeight = NAN;
    QString method;
    QJsonObject parameters;
};

namespace SpectrumProcessing {
ProcessedSpectrum fromResult(const SpectrumResult &, const CrystalRois &, bool normalized);
ProcessedSpectrum load(const QString &, const SpectrumHdfLayout &);
ProcessedSpectrum loadWideReference(const QString &, const SpectrumHdfLayout &, const WideScanCalibration *,
                                    const std::atomic_bool &);
SpectrumSamples selectSamples(const ProcessedSpectrum &, const QStringList &crystals, double lo, double hi);
ProcessedSpectrum spliceCompton(const ProcessedSpectrum &wide, const ProcessedSpectrum &near,
                                bool averageNear, const std::atomic_bool &);
QMap<QString, SpectrumGeometry> loadGeometry(const QString &, const ProcessedSpectrum &);
ProcessedSpectrum replaceEdges(const ProcessedSpectrum &, const QVector<ProcessedSpectrum> &,
                               const std::atomic_bool &);
ProcessedSpectrum stitch(const QStringList &, const SpectrumHdfLayout &, const std::atomic_bool &);
CrystalBackground fit(const QVector<double> &energy, const QVector<double> &intensity,
                      const BackgroundSettings &, const std::atomic_bool &);
void save(const QString &path, const ProcessedSpectrum &,
          const QMap<QString, CrystalBackground> &backgrounds = {},
          const QMap<QString, QMap<QString, QVector<double>>> &components = {});
void saveText(const QString &path, const ProcessedSpectrum &, const QMap<QString, CrystalBackground> &backgrounds,
              const QMap<QString, QMap<QString, QVector<double>>> &components);
} // namespace SpectrumProcessing
