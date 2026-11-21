#pragma once
#include "ComptonProcessing.h"
#include <QFutureWatcher>
#include <QWidget>
#include <functional>
#include <optional>

class QCheckBox;
class QComboBox;
class QDoubleSpinBox;
class QLabel;
class QLineEdit;
class QListWidget;
class QListWidgetItem;
class QPushButton;
class QToolButton;
class SpectrumPlot;

class SpectrumProcessingPage : public QWidget {
    Q_OBJECT
  public:
    SpectrumProcessingPage(QLineEdit &hfTable, QLineEdit &composition, QLineEdit &hfShells, QWidget *parent = nullptr);
    ~SpectrumProcessingPage() override;
    bool isBusy() const {
        return m_busy;
    }
    void setCurrentResult(ProcessedSpectrum result);
    void setAvailableFiles(const QStringList &files) { m_availableFiles = files; }
    void setWideCalibration(std::function<std::optional<WideScanCalibration>()> getCalibration) {
        m_getWideCalibration = std::move(getCalibration);
    }
    QMap<QString, QWidget *> configurationFields() const;
    static QJsonObject backgroundConfigurationDefaults();
    void resetInputs();
    QString wideFile() const { return m_wideFile; }
    QJsonValue batchSelection() const;
    void restoreEdgeInputs(const QString &wideFile, const QJsonValue &batchTargets);
    void restoreWideIfNeeded();
    void setExportIdentity(const QString &stem, const QJsonObject &identity) {
        m_exportStem = stem;
        m_exportIdentity = identity;
    }
    bool prepareAnalysisInput();
    bool ensureAnalysis();
    QWidget *analysisPage() const {
        return m_analysisPage;
    }
    QWidget *resultPage() const {
        return m_resultPage;
    }

  signals:
    void busyChanged();
    void analysisInputReady();
    void analysisFinished(bool usable);
    void analysisParametersEdited();
    void analysisBatchCompleted(bool usable, bool partial);
    void resultSaved(const QString &format, const QString &path);

  private:
    enum Acquisition { Continuous, Separate };
    enum Scope { PerCrystal, ByModule, TotalSpectrum };
    enum class ResultFormat { Hdf5, Text };
    struct Work {
        std::optional<ProcessedSpectrum> wide, stitched;
        std::optional<StagedBackgroundResult> analysis;
        std::optional<CrystalBackground> term;
        QString target;
        QString wideFile, exportFormat, exportPath;
        BackgroundComponent component = BackgroundComponent::Linear;
        quint64 analysisRevision = 0;
        quint64 wideRevision = 0;
        quint64 stitchRevision = 0;
        bool enterAnalysis = false;
        QString message;
        bool error = false;
    };
    SpectrumHdfLayout hdfLayout() const;
    StagedBackgroundSettings settings(int last = backgroundComponentCount - 1) const;
    const ProcessedSpectrum &peakInput() const;
    QStringList previewCrystals(QStringList labels) const;
    void refreshPreviewChoices(const QStringList &labels);
    QString componentName(BackgroundComponent) const;
    QString targetName(const QString &) const;
    bool componentAvailable(BackgroundComponent) const;
    bool canFit(BackgroundComponent) const;
    bool previewComplete() const;
    QStringList batchTargets() const;
    void updateBatchTargets(const std::function<bool(const QListWidgetItem &)> &selected);
    void chooseBatchModules();
    void fitOne(BackgroundComponent);
    void runBatch();
    void loadWide(const QString &path);
    void startSplice(bool enterAnalysis = false);
    void resetSplice();
    void start(const QString &message, std::function<Work()> job, bool cancellable = true);
    void message(const QString &text, bool error = false);
    void resetAnalysis(int first = 0);
    void resetBatch();
    void refreshTargets();
    void refreshInput();
    void refreshAnalysisInput();
    void refreshAnalysis();
    void refreshResult();
    void saveResult(ResultFormat);
    void refreshButtons();
    QStringList resultLabels() const;

    ProcessedSpectrum m_current, m_wide, m_stitched, m_analysisInput;
    QString m_wideFile, m_exportStem;
    QJsonObject m_exportIdentity;
    std::optional<QStringList> m_configuredBatchTargets;
    QStringList m_availableFiles;
    std::function<std::optional<WideScanCalibration>()> m_getWideCalibration;
    QMap<QString, QVector<double>> m_targets;
    QMap<QString, StagedSpectrumFit> m_previews;
    StagedBackgroundResult m_analysis;
    QFutureWatcher<Work> m_watcher;
    std::atomic_bool m_cancel{false};
    bool m_busy = false;
    quint64 m_analysisRevision = 0, m_wideRevision = 0, m_stitchRevision = 0;
    QWidget *m_analysisPage = nullptr, *m_resultPage = nullptr;
    QWidget *m_inputControls = nullptr, *m_fitControls = nullptr, *m_resultControls = nullptr;
    QPushButton *m_loadWide = nullptr, *m_splice = nullptr,
                *m_run = nullptr, *m_cancelButton = nullptr, *m_saveResult = nullptr, *m_saveText = nullptr;
    QPushButton *m_batchModule = nullptr, *m_batchAll = nullptr, *m_batchInvert = nullptr;
    QToolButton *m_previousPreview = nullptr, *m_nextPreview = nullptr;
    QLineEdit *m_energyPath = nullptr, *m_labelsPath = nullptr, *m_intensityPath = nullptr, *m_unit = nullptr;
    QCheckBox *m_kev = nullptr, *m_energyFirst = nullptr, *m_limitPeak = nullptr, *m_normalizeNear = nullptr;
    QDoubleSpinBox *m_peakLo = nullptr, *m_peakHi = nullptr;
    QLineEdit *m_hfTable = nullptr, *m_composition = nullptr, *m_hfShells = nullptr, *m_protected = nullptr;
    QDoubleSpinBox *m_asymmetryRise = nullptr, *m_asymmetryDecay = nullptr;
    std::array<QLineEdit *, backgroundComponentCount> m_windows{};
    std::array<QCheckBox *, backgroundComponentCount> m_enabledTerms{};
    std::array<QWidget *, backgroundComponentCount> m_termControls{};
    std::array<QPushButton *, backgroundComponentCount> m_fitTerms{};
    std::array<QLabel *, backgroundComponentCount> m_termStatus{};
    QDoubleSpinBox *m_resultLo = nullptr, *m_resultHi = nullptr;
    QComboBox *m_inputCrystal = nullptr, *m_fitCrystal = nullptr, *m_resultChoice = nullptr,
              *m_scope = nullptr, *m_resultView = nullptr, *m_acquisition = nullptr;
    QListWidget *m_batchTargets = nullptr;
    QLabel *m_analysisInputInfo = nullptr;
    QLabel *m_inputInfo = nullptr, *m_spliceInfo = nullptr, *m_fitInfo = nullptr, *m_resultInfo = nullptr,
           *m_status[3]{};
    SpectrumPlot *m_inputPlot = nullptr, *m_fitPlot = nullptr,
                 *m_resultPlot = nullptr;
};
