#pragma once
#include "Spectroscopy.h"
#include <QDialog>
#include <QFutureWatcher>

class QTabWidget;
class QLineEdit;
class QComboBox;
class QCheckBox;
class QDoubleSpinBox;
class QSpinBox;
class QTableWidget;
class QListWidget;
class QLabel;
class QPushButton;
class QProgressBar;
class ScientificPlot;
class SpectrumPlot;
class SpectrumProcessingPage;

class SpectroscopyDialog : public QDialog {
    Q_OBJECT
  public:
    explicit SpectroscopyDialog(const QString &initialFile = {}, QWidget *parent = nullptr);
    ~SpectroscopyDialog() override;
    void selectStep(int step);
    bool isBusy() const {
        return m_busy;
    }
    bool hasCalibration() const {
        return m_haveCalibration;
    }
    bool hasScans() const {
        return m_haveScans;
    }
    bool hasResult() const {
        return m_haveResult;
    }
    QString statusText() const;
    bool loadRois(const QString &path);
    void addScans(const QStringList &files);
    void runStep(int step);
    void setAvailableFiles(const QStringList &files) {
        m_availableFiles = files;
    }
    void useRoiDraft(const QString &file, const QString &path);
    const SpectrumResult &result() const {
        return m_result;
    }

  signals:
    void roiRequested(const QString &file);

  private:
    struct Work {
        ElasticCalibration calibration;
        QVector<SpectrumScan> scans;
        SpectrumResult result;
        QString error;
    };
    SpectrumSettings settings() const;
    void buildUi();
    void invalidate(int stage);
    void refreshState();
    void refreshActions();
    void refreshCalibration();
    void refreshImages();
    void refreshElasticCurve();
    void refreshScanCurves();
    void refreshResult();
    void refreshResultCurve();
    void loadDraft();
    void save();
    void message(const QString &text, bool error = false);
    CrystalRois m_rois;
    ElasticCalibration m_calibration;
    QVector<SpectrumScan> m_scans;
    SpectrumResult m_result;
    SpectrumSettings m_calibrationSettings, m_extractionSettings, m_resultSettings;
    QFutureWatcher<Work> m_watcher;
    std::atomic_bool m_cancel{false};
    bool m_busy = false, m_updating = false, m_haveCalibration = false, m_haveScans = false,
         m_haveResult = false;
    int m_job = 0;
    QStringList m_availableFiles;
    QPair<QString, QString> m_pendingRoi;
    QPushButton *m_nextButton = nullptr, *m_previousButton = nullptr;
    QTabWidget *m_tabs = nullptr;
    QLineEdit *m_elasticFile = nullptr, *m_roiPath = nullptr;
    QComboBox *m_elasticEnergy = nullptr, *m_xrsEnergy = nullptr, *m_elasticI0 = nullptr, *m_xrsI0 = nullptr,
              *m_model = nullptr;
    QCheckBox *m_normalizeElastic = nullptr, *m_normalizeXrs = nullptr, *m_repair = nullptr,
              *m_filter = nullptr, *m_auto = nullptr, *m_overlap = nullptr;
    QDoubleSpinBox *m_fraction = nullptr, *m_centerMin = nullptr, *m_centerMax = nullptr, *m_fwhm = nullptr,
                   *m_r2 = nullptr;
    QSpinBox *m_roiSize = nullptr;
    QListWidget *m_files = nullptr;
    QTableWidget *m_corrections = nullptr, *m_fits = nullptr, *m_inclusion = nullptr;
    QComboBox *m_imageSource = nullptr, *m_elasticView = nullptr, *m_scanChoice = nullptr,
              *m_crystalChoice = nullptr, *m_resultChoice = nullptr;
    ScientificPlot *m_images[2]{};
    SpectrumPlot *m_elasticPlot = nullptr, *m_i0Plot = nullptr, *m_scanPlot = nullptr,
                 *m_resultPlot = nullptr;
    SpectrumProcessingPage *m_processing = nullptr;
    QComboBox *m_angleSource = nullptr;
    QMap<QString, QDoubleSpinBox *> m_angles;
    QVector<QCheckBox *> m_modules;
    QDoubleSpinBox *m_diameter = nullptr, *m_distance = nullptr, *m_pitch = nullptr, *m_qMin = nullptr,
                   *m_qMax = nullptr, *m_step = nullptr;
    QLabel *m_status = nullptr, *m_roiStatus = nullptr;
    QProgressBar *m_progress = nullptr;
    QPushButton *m_calibrateButton = nullptr, *m_extractButton = nullptr, *m_combineButton = nullptr,
                *m_saveButton = nullptr, *m_cancelButton = nullptr;
};
