#pragma once
#include "BondVolume.h"
#include <QColor>
#include <QFutureWatcher>
#include <QWidget>

class BondFigures;
class QCheckBox;
class QComboBox;
class QDialog;
class QDoubleSpinBox;
class QLabel;
class QLineEdit;
class QPlainTextEdit;
class QProgressBar;
class QPushButton;
class QSpinBox;
class QTabWidget;
class QTableWidget;
class ScientificPlot;
class SpectrumPlot;

class BondWorkspace : public QWidget {
    Q_OBJECT
  public:
    explicit BondWorkspace(const QString &sample, QWidget *parent = nullptr);
    ~BondWorkspace() override;
    bool importConfiguration(const QString &);
    bool saveConfiguration(const QString &);
    void applyConfiguration(const QJsonObject &);
    void setAvailableFiles(const QStringList &files) {
        m_availableFiles = files;
    }
    void useRoiDraft(const QString &file, const QString &path);
    void start(bool validateOnly);
    void openVolume(const QString &);
    bool isBusy() const;
    bool hasResult() const {
        return m_volume.voxelCount() > 0;
    }
    void selectStep(int);

  signals:
    void roiRequested(const QString &file);

  private:
    struct Work {
        std::optional<ImagingValidation> validation;
        std::optional<BondVolume> volume;
        std::optional<QVector<double>> processed;
        std::optional<std::array<QVector<BondTriangle>, 2>> surfaces;
        std::array<QColor, 2> surfaceColors;
        QString message, error, meshCaption, processingMethod;
        bool normalized = false, subtracted = false;
    };
    void launch(const QString &, std::function<Work()>);
    void status(const QString &, bool error = false);
    QJsonObject configuration() const;
    void showValidation(const ImagingValidation &);
    void installVolume(BondVolume);
    QVector<double> energy() const;
    QString energyLabel() const;
    QString intensityLabel() const;
    bool showingProcessed() const;
    const QVector<double> &values() const;
    void syncViewControls();
    void invalidateInputs();
    void refreshVoxel(bool refreshSection = true, bool refreshSlice = true);
    void refreshMap();
    void refreshMontage();
    void createComptonOptions();
    void applySpectralProcessing();
    void buildSurfaces();
    void changeEnergyUnits();

    QJsonObject m_base;
    QMap<QString, QWidget *> m_inputs;
    QStringList m_availableFiles;
    QByteArray m_roiBytes;
    QPair<QString, QString> m_pendingRoi;
    bool m_refreshing = false;
    std::atomic_bool m_stop{false};
    QFutureWatcher<Work> m_watcher;
    BondVolume m_volume;
    QVector<double> m_processed, m_field;
    double m_referenceKev = 0., m_lo = 0., m_hi = 1.;
    bool m_transfer = false, m_normalized = false, m_subtracted = false;
    QString m_fieldLabel;
    QTabWidget *m_tabs = nullptr, *m_figures = nullptr;
    QLineEdit *m_sample = nullptr, *m_roi = nullptr, *m_output = nullptr, *m_anchor = nullptr,
              *m_background = nullptr, *m_normalization = nullptr;
    QDialog *m_comptonOptions = nullptr;
    QLineEdit *m_wideSpectrum = nullptr, *m_hfTable = nullptr, *m_composition = nullptr;
    std::array<QLineEdit *, 5> m_comptonWindows{};
    QDoubleSpinBox *m_highQ = nullptr, *m_coreScale = nullptr, *m_valenceCutoff = nullptr;
    QCheckBox *m_fitGains = nullptr, *m_loss = nullptr, *m_secondSurface = nullptr;
    QDoubleSpinBox *m_analyzer = nullptr, *m_ea = nullptr, *m_eb = nullptr, *m_level = nullptr,
                   *m_level2 = nullptr, *m_tolerance = nullptr;
    QComboBox *m_mode = nullptr, *m_dataSource = nullptr, *m_spectralStyle = nullptr,
              *m_backgroundMode = nullptr;
    QSpinBox *m_channel = nullptr, *m_x = nullptr, *m_y = nullptr, *m_z = nullptr, *m_firstZ = nullptr,
             *m_countZ = nullptr, *m_strideZ = nullptr, *m_iterations = nullptr;
    QLabel *m_status = nullptr, *m_mapInfo = nullptr, *m_voxelInfo = nullptr, *m_quality = nullptr;
    QLabel *m_selectedEnergy = nullptr;
    QProgressBar *m_progress = nullptr;
    QPushButton *m_cancel = nullptr, *m_run = nullptr;
    QWidget *m_viewBar = nullptr, *m_singleEnergy = nullptr, *m_dualEnergy = nullptr, *m_xControl = nullptr,
            *m_spectralStyleControl = nullptr;
    QTabWidget *m_plots = nullptr;
    QTableWidget *m_crystals = nullptr;
    QPlainTextEdit *m_diagnostics = nullptr;
    ScientificPlot *m_slice = nullptr;
    SpectrumPlot *m_spectrum = nullptr;
    BondFigures *m_montage = nullptr, *m_spectralFigure = nullptr, *m_surface = nullptr;
};
