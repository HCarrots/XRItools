#pragma once
#include "Spectroscopy.h"
#include <QDialog>
#include <QJsonArray>
#include <QTemporaryDir>

class QComboBox;
class QSpinBox;
class QDoubleSpinBox;
class QLabel;
class QPushButton;
class QCheckBox;
class QListWidget;
class ScientificPlot;
class SpectrumPlot;

class SpectralMaskDialog : public QDialog {
    Q_OBJECT
  public:
    SpectralMaskDialog(const SpectrumSettings &, const QVector<CrystalRoi> &, const QMap<QString, QSize> &,
                       int initialRow, const QStringList &references, QWidget *parent);
    SpectralSelection selection() const;
    int currentRow() const {
        return row;
    }

    const QMap<int, SpectralSelection> &acceptedMasks() const { return accepted; }

  protected:
    void keyPressEvent(QKeyEvent *event) override;

  private:
    struct MaskOptions {
        int component = 0;
        int thresholdMode = 0;
        double threshold = .2;
        bool flip = false;
    };
    struct Review {
        QVector<QPoint> pixels;
        QVector<double> energyKev;
        SpectralComponents result;
        QVector<int> rows;
        QJsonArray history;
        MaskOptions options;
        QString error;
        bool noComponent = false;
    };
    struct MaskPreview {
        QVector<double> covariance;
        QVector<int> selected;
        double cutoff = 0.;
    };
    void switchCrystal(int index);
    void showReview();
    void refreshDecisions();
    void recordDecision(bool useComponent);
    void applyBatch();
    void resetDecisions();
    void editPreview();
    QMap<int, SpectralSelection> accepted;
    QSet<int> declined;
    bool switching = false;
    QVector<CrystalRoi> crystals;
    QMap<QString, QSize> sizes;
    QMap<int, Review> reviews;
    QTemporaryDir pixelCache;
    QVector<QVector<double>> intensity;
    int row = -1;
    QListWidget *crystal;
    QPushButton *previous, *next;
    QLabel *position, *reviewProgress;
    void calculate(bool refine);
    bool calculateAll();
    bool runDecomposition(const QVector<int> &targets, bool refine = false);
    QString pixelCachePath(int crystalRow) const;
    void cacheIntensity(int crystalRow, const QVector<QVector<double>> &values) const;
    QVector<QVector<double>> loadIntensity(int crystalRow) const;
    void clearPreview();
    void preview();
    MaskOptions maskOptions() const;
    void setThresholdMode(int mode, double value);
    static MaskPreview selectPixels(const Review &, const MaskOptions &);
    QJsonObject roundRecord(const Review &, const MaskOptions &, const MaskPreview &) const;
    SpectralSelection selection(const Review &, const MaskOptions &, const MaskPreview &) const;
    SpectrumSettings settings;
    CrystalRoi roi;
    MaskPreview currentMask;
    QComboBox *reference, *method, *component, *thresholdMode;
    QSpinBox *count;
    QDoubleSpinBox *threshold;
    QCheckBox *flip;
    QLabel *status;
    QPushButton *decomposeAll, *iterate, *apply, *batchApply, *keepOriginal, *updateRoi;
    SpectrumPlot *curves, *scores, *comparison;
    ScientificPlot *map;
};
