#pragma once
#include "CrystalRoi.h"
#include "HdfFile.h"
#include <QMap>
#include <QSet>
#include <cmath>
#include <functional>
#include <limits>

struct SpectrumSettings {
    QString elasticFile, elasticEnergy = "auto", xrsEnergy = "auto";
    QString elasticI0 = "/entry/data/D_SiC_Current1", xrsI0 = "/entry/data/D_SiC_Current1";
    bool normalizeElastic = true, normalizeXrs = true, repairGlitches = true;
    bool useFilter = true, autoRoi = false, removeOverlap = true;
    double filterFraction = 0.15;
    int roiSize = 20;
    QVector<QPoint> shifts = QVector<QPoint>(75), expansions = QVector<QPoint>(75);
    int fitModel = 0; // 0: Lorentzian (HWHM), 1: Gaussian (sigma)
    double centerMin = 9.67, centerMax = 9.69, maxFwhmEv = 2, minRSquared = 0.8;
    int angleSource = 0; // elastic file, first XRS file, manual
    QMap<QString, double> angles;
    double crystalDiameter = 100, crystalDistance = 1070, modulePitchDegrees = 7;
    QStringList modules{"VB", "HL", "VU", "HB", "VD"};
    double qMin = 0, qMax = 10, energyStepEv = 0.2;
    QSet<int> excludedCrystals;
};
struct ElasticFit {
    bool converged = false, accepted = false;
    double centerKev = NAN, fwhmEv = NAN, amplitude = NAN, background = NAN, rSquared = NAN;
    QString reason;
    QVector<double> fitted;
};
struct SpectrumScan {
    QString file, energyPath, i0Path;
    QVector<double> energyKev, i0Before, i0After;
    QVector<int> glitches;
    QVector<QVector<double>> raw, filtered;
    QVector<DataSlice> images;
};
struct ElasticCalibration {
    CrystalRois original, adjusted;
    QVector<QSize> sizes;
    QVector<QVector<int>> pixels, filteredPixels;
    QVector<DataSlice> images, maskImages;
    SpectrumScan scan;
    QVector<ElasticFit> fits;
    int overlapPixels = 0;
};
struct MomentumStats {
    double mean = NAN, dMean = NAN, range = NAN, dRange = NAN;
};
struct SpectrumResult {
    QVector<double> energyTransferEv, total;
    QVector<QVector<double>> crystals;
    QMap<QString, QVector<double>> modules;
    QVector<MomentumStats> momentum;
    QVector<bool> included;
    QMap<QString, double> angles;
    QStringList reasons;
};
class Spectroscopy {
  public:
    using Progress = std::function<void(const QString &, int)>;
    static ElasticCalibration calibrate(const SpectrumSettings &, const CrystalRois &,
                                        const std::atomic_bool &, Progress = {});
    static QVector<SpectrumScan> extract(const SpectrumSettings &, const ElasticCalibration &,
                                         const QStringList &, const std::atomic_bool &, Progress = {});
    static SpectrumResult combine(const SpectrumSettings &, const ElasticCalibration &,
                                  const QVector<SpectrumScan> &, const std::atomic_bool &);
    static ElasticFit fitPeak(const QVector<double> &energyKev, const QVector<double> &intensity,
                              const SpectrumSettings &);
    static QPair<double, double> momentum(const QString &label, double angleDeg, double incidentKev,
                                          double finalKev, const SpectrumSettings &);
    static QMap<QString, double> readAngles(const QString &file);
    static bool saveResult(const QString &directory, const SpectrumSettings &, const ElasticCalibration &,
                           const QVector<SpectrumScan> &, const SpectrumResult &, QString *error);
};
