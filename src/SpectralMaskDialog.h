#pragma once
#include "Spectroscopy.h"
#include <QCache>
#include <QDialog>
#include <QJsonArray>

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
    struct Review {
        PixelSpectra data;
        SpectralComponents result;
        QVector<int> rows;
        QJsonArray history;
        int round, component, thresholdMode;
        double threshold;
        bool flip;
    };
    void switchCrystal(int index);
    void refreshDecisions();
    void recordDecision(bool useComponent);
    void resetDecisions();
    void editPreview();
    QMap<int, SpectralSelection> accepted;
    QSet<int> declined, computed;
    bool switching = false;
    QVector<CrystalRoi> crystals;
    QMap<QString, QSize> sizes;
    QCache<int, Review> reviews{128 * 1024}; // KiB of recently viewed pixel spectra
    int row = -1;
    bool autoDecompose = false;
    QListWidget *crystal;
    QPushButton *previous, *next;
    QLabel *position, *reviewProgress;
    void calculate(bool nextRound);
    void clearResult();
    void preview();
    QJsonObject roundRecord() const;
    QJsonArray history;
    double cutoffValue = 0;
    SpectrumSettings settings;
    CrystalRoi roi;
    QSize size;
    PixelSpectra data;
    SpectralComponents result;
    QVector<int> rows, selected;
    int round = 0;
    QComboBox *reference, *method, *component, *thresholdMode;
    QSpinBox *count;
    QDoubleSpinBox *threshold;
    QCheckBox *flip;
    QLabel *status;
    QPushButton *iterate, *apply, *keepOriginal, *updateRoi;
    SpectrumPlot *curves, *scores, *comparison;
    ScientificPlot *map;
};
