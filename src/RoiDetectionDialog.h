#pragma once

#include <QPixmap>
#include <QProcess>
#include <QTemporaryDir>
#include <QWidget>

class QComboBox;
class QDoubleSpinBox;
class QLabel;
class QLineEdit;
class QPushButton;
class QSpinBox;
class QStackedWidget;
class QTableWidget;

class RoiDetectionDialog : public QWidget {
    Q_OBJECT
  public:
    explicit RoiDetectionDialog(const QString &source, QWidget *parent = nullptr);
    ~RoiDetectionDialog() override;
    void setBrowserSource(const QString &source);

  protected:
    bool eventFilter(QObject *watched, QEvent *event) override;

  private:
    void inspectSource();
    bool syncBrowserSource();
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
    QString m_browserSource, m_source;
    QStackedWidget *m_pages = nullptr;
    QWidget *m_detection = nullptr;
    QWidget *m_parameters = nullptr;
    QLineEdit *m_python = nullptr, *m_weights = nullptr;
    QComboBox *m_detector = nullptr, *m_mode = nullptr;
    QSpinBox *m_frame = nullptr;
    QDoubleSpinBox *m_score = nullptr, *m_maskThreshold = nullptr;
    QLabel *m_message = nullptr, *m_preview = nullptr;
    QTableWidget *m_table = nullptr;
    QPushButton *m_run = nullptr, *m_cancel = nullptr, *m_export = nullptr;
    QPushButton *m_matchSaved = nullptr, *m_matchCurrent = nullptr;
    QByteArray m_result;
    QPixmap m_overlay;
    bool m_cancelled = false;
};
