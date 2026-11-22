#pragma once
#include "Spectroscopy.h"
#include <QDialog>
#include <QFutureWatcher>
#include <QJsonArray>

class QTabWidget;
class QDir;
class QLineEdit;
class QMenu;
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
    bool isBusy() const;
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
    void useGlobalStatus();
    bool importConfiguration(const QString &path);
    bool saveConfiguration(const QString &path);
    bool loadRois(const QString &path);
    bool addEdgeGroup(const QString &element, const QString &shell);
    void addScans(const QStringList &files);
    void runStep(int step);
    void setAvailableFiles(const QStringList &files);
    void useRoiDraft(const QString &file, const QString &path);
    const SpectrumResult &result() const {
        return m_result;
    }

  signals:
    void roiRequested(const QString &file);
    void statusChanged(const QString &text);

  private:
    enum class EdgeStatus { Pending, Combined, Partial, Completed };
    enum class MaskSource { Current, Parameters, Spectral, Original };
    struct EdgeConfiguration {
        QString elasticFile, roiFile, wideFile;
        CrystalRois rois;
        QVector<QSize> sizes;
        QMap<int, SpectralSelection> masks;
        QSet<int> exclusions;
        QSet<int> manualPages;
        QJsonObject parameters;
        QJsonArray corrections;
        QList<QPair<QString, bool>> scans;
        QJsonValue batchTargets;
        bool originalRoiOnly = false;
    };
    struct EdgeGroup {
        QString id, name;
        EdgeConfiguration config;
        EdgeStatus status = EdgeStatus::Pending;
        QJsonObject lastSaved;
        QString targetElement, targetShell;
    };
    struct Work {
        ElasticCalibration masks, calibration;
        QVector<SpectrumScan> scans;
        SpectrumResult result;
        QString error;
        MaskSource maskSource = MaskSource::Current;
    };
    SpectrumSettings settings(MaskSource source = MaskSource::Current) const;
    QMap<QString, QWidget *> configurationFields() const;
    EdgeConfiguration captureEdgeConfiguration() const;
    EdgeConfiguration readEdgeConfiguration(const QJsonObject &, const QDir &) const;
    QJsonObject writeEdgeConfiguration(const EdgeConfiguration &) const;
    void applyEdgeConfiguration(const EdgeConfiguration &);
    bool validEdgeName(const QString &name, int except = -1) const;
    void chooseNewEdgeGroup();
    void switchEdgeGroup(int index);
    void refreshEdgeGroups();
    void chooseEdgeHfShells();
    void populateScanAddMenu(QMenu &menu);
    void chooseBrowserScans(const QStringList &files);
    QSet<QString> listedScanPaths() const;
    void setEdgeStatus(EdgeStatus);
    void refreshExportIdentity();
    QString sampleStem() const;
    QString resultStem() const;
    QJsonObject edgeIdentity() const;
    static QString edgeStatusKey(EdgeStatus);
    void buildUi();
    void invalidate(int stage, bool preserveEdgeStatus = false);
    void calculationParametersEdited(int stage);
    void processingParametersEdited(int page);
    int firstManualPage() const;
    QString manualCalculationHint(int page) const;
    void refreshState();
    void refreshActions();
    void advanceToPage(int target);
    void runStep(int step, MaskSource source);
    void refreshCalibration();
    void editSpectralMask();
    void maskParametersEdited();
    void applySpectralMasks(QMap<int, SpectralSelection> masks);
    void refreshSpectralMaskLabels();
    void refreshImages();
    void refreshElasticCurve();
    void refreshScanCurves();
    void refreshResult();
    void refreshResultCurve();
    void loadDraft();
    void save();
    void message(const QString &text, bool error = false);
    CrystalRois m_rois;
    QVector<QSize> m_roiSizes;
    QSet<int> m_configExclusions;
    QSet<int> m_manualPages;
    QMap<int, SpectralSelection> m_spectralMasks;
    ElasticCalibration m_masks, m_calibration;
    QVector<DataSlice> m_projection;
    QVector<SpectrumScan> m_scans;
    SpectrumResult m_result;
    SpectrumSettings m_calibrationSettings, m_extractionSettings, m_resultSettings;
    QFutureWatcher<Work> m_watcher;
    std::atomic_bool m_cancel{false};
    bool m_busy = false, m_updating = false, m_haveCalibration = false, m_haveScans = false,
         m_haveResult = false;
    bool m_useOriginalRoi = false;
    int m_job = 0;
    int m_nextPage = -1;
    QStringList m_availableFiles;
    QList<EdgeGroup> m_edgeGroups;
    EdgeConfiguration m_edgeDraft;
    int m_activeEdge = -1;
    QPair<QString, QString> m_pendingRoi;
    QPushButton *m_nextButton = nullptr, *m_previousButton = nullptr;
    QPushButton *m_fitElasticButton = nullptr;
    QPushButton *m_calculateMaskButton = nullptr, *m_resetMaskButton = nullptr;
    QDoubleSpinBox *m_glitchDrop = nullptr;
    QTabWidget *m_tabs = nullptr;
    QLineEdit *m_elasticFile = nullptr, *m_roiPath = nullptr, *m_sampleName = nullptr;
    QLineEdit *m_hfTable = nullptr, *m_hfComposition = nullptr, *m_hfShells = nullptr;
    QCheckBox *m_measureCompton = nullptr;
    QComboBox *m_elasticEnergy = nullptr, *m_xrsEnergy = nullptr, *m_elasticI0 = nullptr, *m_xrsI0 = nullptr,
              *m_model = nullptr;
    QCheckBox *m_normalizeElastic = nullptr, *m_normalizeXrs = nullptr, *m_repair = nullptr,
              *m_filter = nullptr, *m_auto = nullptr, *m_overlap = nullptr;
    QDoubleSpinBox *m_fraction = nullptr, *m_centerMin = nullptr, *m_centerMax = nullptr, *m_fwhm = nullptr,
                   *m_r2 = nullptr;
    QSpinBox *m_roiSize = nullptr;
    QListWidget *m_files = nullptr;
    QTableWidget *m_edgeTable = nullptr;
    QPushButton *m_newEdge = nullptr, *m_renameEdge = nullptr, *m_removeEdge = nullptr,
                *m_addScans = nullptr, *m_removeScans = nullptr, *m_saveConfig = nullptr;
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
    QPushButton *m_extractButton = nullptr, *m_combineButton = nullptr,
                *m_saveButton = nullptr, *m_cancelButton = nullptr;
};
