#pragma once
#include "CrystalRoi.h"
#include "HdfFile.h"
#include <QDialog>
#include <QFutureWatcher>
#include <QTimer>
#include <atomic>

class ScientificPlot;
class QTableWidget;
class QLabel;
class QComboBox;
class QSpinBox;
class QCheckBox;
class QPushButton;
class QGroupBox;
class RoiDialog : public QDialog {
    Q_OBJECT
  public:
    RoiDialog(const QString &file, const QString &workflow, QWidget *parent = nullptr, bool autoDraft = true);
    ~RoiDialog() override;
    bool isLoading() const {
        return m_loading;
    }
    bool imagesReady() const;
    const CrystalRois &rois() const {
        return m_rois;
    }
    bool saveTo(const QString &path);
    bool loadFrom(const QString &path);

  protected:
    void closeEvent(QCloseEvent *) override;
    void reject() override;

  private:
    void loadImages();
    void refresh();
    void selectRow();
    QVector<int> selectedRows() const;
    void assign(int detector, const QRect &rect);
    void nextPending();
    void storeDraft();
    void updateProfiles();
    QString m_file, m_draftPath;
    CrystalRois m_rois;
    QVector<QSize> m_sizes;
    QVector<quint64> m_frames;
    ScientificPlot *m_plots[2]{};
    QTableWidget *m_table = nullptr;
    QLabel *m_progress = nullptr, *m_message = nullptr, *m_current = nullptr;
    QComboBox *m_mode = nullptr;
    QSpinBox *m_frame = nullptr, *m_coords[4]{};
    QCheckBox *m_advance = nullptr;
    QPushButton *m_finish = nullptr;
    QPushButton *m_apply = nullptr, *m_invalid = nullptr, *m_clear = nullptr;
    QCheckBox *m_showProfiles = nullptr;
    QGroupBox *m_profileGroup = nullptr;
    QComboBox *m_profileDetector = nullptr;
    QSpinBox *m_profileX = nullptr, *m_profileY = nullptr;
    ScientificPlot *m_profiles[2]{};
    QPushButton *m_profileCenter = nullptr;
    QFutureWatcher<QVector<DataSlice>> m_watcher;
    std::atomic_bool m_cancel{false};
    std::atomic<int> m_loadPercent{0}, m_activeDetector{0};
    QTimer m_progressTimer;
    bool m_loading = false, m_autoDraft = true;
};
