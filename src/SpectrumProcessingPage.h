#pragma once
#include "ComptonProcessing.h"
#include <QFutureWatcher>
#include <QWidget>
#include <optional>

class QCheckBox;
class QComboBox;
class QDoubleSpinBox;
class QLabel;
class QLineEdit;
class QPushButton;
class SpectrumPlot;

class SpectrumProcessingPage : public QWidget {
  public:
    explicit SpectrumProcessingPage(QWidget *parent = nullptr);
    ~SpectrumProcessingPage() override;
    void setCurrentResult(ProcessedSpectrum result);
    QWidget *analysisPage() const {
        return m_analysisPage;
    }
    QWidget *resultPage() const {
        return m_resultPage;
    }

  private:
    enum InputSource { CurrentSpectrum, ExternalSpectrum, StitchingSpectrum };
    struct Work {
        std::optional<ProcessedSpectrum> wide, stitched, external;
        std::optional<ComptonResult> analysis;
        quint64 analysisRevision = 0;
        QString message;
        bool error = false;
    };
    SpectrumHdfLayout hdfLayout() const;
    ComptonSettings settings() const;
    const ProcessedSpectrum &analysisInput() const;
    void start(const QString &message, std::function<Work()> job, bool cancellable = true);
    void message(const QString &text, bool error = false);
    void resetAnalysis();
    void refreshInput();
    void refreshAnalysisInput();
    void refreshAnalysis();
    void refreshResult();
    void refreshButtons();
    QStringList resultLabels() const;

    ProcessedSpectrum m_current, m_wide, m_stitched, m_external;
    ComptonResult m_analysis;
    QFutureWatcher<Work> m_watcher;
    std::atomic_bool m_cancel{false};
    bool m_busy = false;
    quint64 m_analysisRevision = 0;
    QWidget *m_analysisPage = nullptr, *m_resultPage = nullptr;
    QWidget *m_inputControls = nullptr, *m_fitControls = nullptr, *m_resultControls = nullptr;
    QPushButton *m_useCurrent = nullptr, *m_edges = nullptr, *m_geometry = nullptr, *m_saveStitch = nullptr,
                *m_run = nullptr, *m_cancelButton = nullptr, *m_saveResult = nullptr;
    QLineEdit *m_energyPath = nullptr, *m_labelsPath = nullptr, *m_intensityPath = nullptr, *m_unit = nullptr;
    QCheckBox *m_kev = nullptr, *m_energyFirst = nullptr;
    QLineEdit *m_hfTable = nullptr, *m_composition = nullptr, *m_windows[5]{};
    QDoubleSpinBox *m_qMin = nullptr, *m_coreScale = nullptr, *m_cutoff = nullptr;
    QDoubleSpinBox *m_resultLo = nullptr, *m_resultHi = nullptr;
    QPushButton *m_loadDirect = nullptr;
    QComboBox *m_inputCrystal = nullptr, *m_fitCrystal = nullptr, *m_resultChoice = nullptr,
              *m_analysisSource = nullptr, *m_backgroundMode = nullptr;
    QLabel *m_analysisInputInfo = nullptr;
    QLabel *m_inputInfo = nullptr, *m_fitInfo = nullptr, *m_resultInfo = nullptr, *m_status[3]{};
    SpectrumPlot *m_inputPlot = nullptr, *m_referencePlot = nullptr, *m_fitPlot = nullptr,
                 *m_resultPlot = nullptr;
};
