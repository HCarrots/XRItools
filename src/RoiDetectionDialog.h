#pragma once

#include <QDialog>
#include <QPixmap>
#include <QProcess>
#include <QTemporaryDir>

class QCheckBox;
class QComboBox;
class QDoubleSpinBox;
class QLabel;
class QLineEdit;
class QPushButton;
class QSpinBox;
class QTableWidget;

class RoiDetectionDialog : public QDialog {
    Q_OBJECT
  public:
    explicit RoiDetectionDialog(const QString &source, QWidget *parent = nullptr);
    ~RoiDetectionDialog() override;

  protected:
    bool eventFilter(QObject *watched, QEvent *event) override;

  private:
    void inspectSource();
    void updateFrameRange();
    void run();
    void completed(int code, QProcess::ExitStatus status);
    void exportResult();
    void loadResultForMatching();
    void matchAndSave(const QByteArray &result);
    void setRunning(bool running);
    void resizePreview();

    QTemporaryDir m_work;
    QProcess m_process;
    QWidget *m_parameters = nullptr;
    QLineEdit *m_source = nullptr, *m_python = nullptr, *m_weights = nullptr;
    QComboBox *m_detector = nullptr;
    QCheckBox *m_sum = nullptr;
    QSpinBox *m_frame = nullptr;
    QDoubleSpinBox *m_score = nullptr, *m_maskThreshold = nullptr;
    QLabel *m_message = nullptr, *m_preview = nullptr;
    QTableWidget *m_table = nullptr;
    QPushButton *m_run = nullptr, *m_cancel = nullptr, *m_export = nullptr;
    QPushButton *m_match = nullptr;
    QAction *m_matchCurrent = nullptr;
    QByteArray m_result;
    QPixmap m_overlay;
    bool m_cancelled = false;
};
