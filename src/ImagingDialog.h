#pragma once
#include "Imaging.h"
#include "Isosurface.h"
#include <QDialog>
#include <QFutureWatcher>
class QTabWidget;
class QLineEdit;
class QComboBox;
class QCheckBox;
class QTableWidget;
class QPlainTextEdit;
class QSpinBox;
class QLabel;
class QPushButton;
class QProgressBar;
class ScientificPlot;
class SpectrumPlot;
class BondFigures;
class ImagingDialog : public QDialog {
    Q_OBJECT
  public:
    explicit ImagingDialog(const QString &sample = {}, QWidget *parent = nullptr);
    ~ImagingDialog() override;
    void selectStep(int);
    void applyConfiguration(const QJsonObject &);
    QJsonObject configuration() const;
    bool importConfiguration(const QString &);
    bool loadResult(const QString &directory);
    void start(bool validateOnly = false);
    void setAvailableFiles(const QStringList &files) {
        m_availableFiles = files;
    }
    void useRoiDraft(const QString &file, const QString &path);
    bool isBusy() const {
        return m_busy;
    }
    bool hasResult() const {
        return m_haveResult;
    }
    QString statusText() const;
    const ImagingResult &result() const {
        return m_result;
    }
  signals:
    void openFileRequested(const QString &path);
    void roiRequested(const QString &file);

  private:
    struct Work {
        ImagingValidation validation;
        ImagingResult result;
        QVector<SurfaceTriangle> mesh;
        QString error;
    };
    void buildUi();
    void refreshActions();
    void setBusy(bool busy);
    void invalidate();
    void status(const QString &, bool error = false);
    void refreshValidation();
    void refreshResult();
    void refreshImage();
    void refreshModelSource();
    void buildModel();
    const ImageMatrix &volumeData(int method) const;
    void refreshTrajectory();
    void exclusions(bool exclude);
    QJsonObject m_base;
    QStringList m_availableFiles;
    QPair<QString, QString> m_pendingRoi;
    QByteArray m_roiBytes;
    QWidget *m_configTools = nullptr;
    QPushButton *m_nextButton = nullptr, *m_previousButton = nullptr;
    QMap<QString, QWidget *> m_fields;
    QMap<QString, QCheckBox *> m_auto;
    QFutureWatcher<Work> m_watcher;
    ImagingValidation m_validation;
    ImagingResult m_result;
    bool m_busy = false, m_updating = false, m_haveResult = false, m_validateOnly = false;
    bool m_buildingModel = false;
    double m_modelMin = 0., m_modelMax = 0.;
    std::atomic_bool m_cancel{false};
    QTabWidget *m_tabs = nullptr;
    QWidget *m_inputs = nullptr, *m_parameters = nullptr;
    QTableWidget *m_rois = nullptr;
    QPlainTextEdit *m_log = nullptr, *m_summary = nullptr, *m_validationText = nullptr;
    QLabel *m_status = nullptr, *m_info = nullptr, *m_sliceLabel = nullptr, *m_modelRange = nullptr;
    QProgressBar *m_progress = nullptr;
    QPushButton *m_validateButton = nullptr, *m_runButton = nullptr, *m_cancelButton = nullptr,
                *m_importButton = nullptr, *m_exportButton = nullptr, *m_excludeButton = nullptr,
                *m_includeButton = nullptr, *m_openButton = nullptr, *m_buildModelButton = nullptr;
    QComboBox *m_method = nullptr, *m_plane = nullptr, *m_display = nullptr, *m_psfChoice = nullptr;
    QSpinBox *m_slice = nullptr;
    ScientificPlot *m_image = nullptr;
    QComboBox *m_modelSource = nullptr;
    QLineEdit *m_modelLevel = nullptr;
    BondFigures *m_model = nullptr;
    SpectrumPlot *m_sensitivity = nullptr, *m_trajectory = nullptr, *m_spectrum = nullptr;
};
