#pragma once
#include "HdfFile.h"
#include <QFutureWatcher>
#include <QHash>
#include <QMainWindow>
#include <QTimer>
#include <QVector>
#include <functional>

class ScientificPlot;
class QTreeWidget;
class QTreeWidgetItem;
class QLineEdit;
class QLabel;
class QComboBox;
class QSpinBox;
class QSlider;
class QTableView;
class QPlainTextEdit;
class QTabWidget;
class QProgressBar;
class QToolButton;
class SliceTableModel;
class QAction;
class SpectroscopyDialog;
class ImagingDialog;
class BondWorkspace;
class RoiDialog;
class QStackedWidget;

class MainWindow : public QMainWindow {
    Q_OBJECT
  public:
    using Relaunch = std::function<bool(const QStringList &)>;
    explicit MainWindow(QWidget *parent = nullptr, Relaunch relaunch = {});
    ~MainWindow() override;
    void openFile(const QString &path);
    void openFiles(const QStringList &paths);
    void activateFile(const QString &path);
    QStringList openedFiles() const;
    QString currentFile() const {
        return m_file;
    }
    void closeCurrentFile();
    void closeAllFiles();
    bool exportFile(const QString &destination);
    bool isExporting() const {
        return m_exporting;
    }
    QString exportError() const {
        return m_exportError;
    }
    void selectPath(const QString &path);
    bool isReading() const {
        return m_reading || m_debounce.isActive() || m_pending;
    }
    const DataSlice &currentData() const {
        return m_slice;
    }
    void setSliceIndex(int dimension, int index);
    SpectroscopyDialog *openSpectroscopy(int step = -1);
    ImagingDialog *openImaging(int step = -1);
    BondWorkspace *openBondImaging();
    RoiDialog *openRoi(const QString &file = {}, const QString &scheme = {});

  protected:
    void dragEnterEvent(QDragEnterEvent *) override;
    void dropEvent(QDropEvent *) override;
    void closeEvent(QCloseEvent *) override;

  private:
    void buildUi();
    bool restartForLanguage(const QString &language);
    void showWorkspace(QWidget *page, const QString &title = {});
    void updateBrowserTitle();
    QTreeWidgetItem *fileRoot(const QString &file) const;
    QTreeWidgetItem *itemAtPath(QTreeWidgetItem *root, const QString &path);
    void rememberFile();
    void syncFiles();
    void clearPreview();
    void populate(QTreeWidgetItem *parent, const QString &path, int depth);
    void selectItem(QTreeWidgetItem *item);
    void configureAxes();
    void buildSliceControls();
    void scheduleRead();
    void startRead();
    void finishRead();
    void updateFilter(const QString &text);
    QString formatLabel() const;
    QString m_file;
    ScanFormat m_format = ScanFormat::Generic;
    HdfNode m_node;
    DataSlice m_slice;
    quint64 m_generation = 0, m_readGeneration = 0;
    bool m_pending = false, m_configuring = false, m_reading = false;
    QFutureWatcher<DataSlice> m_watcher;
    QTimer m_debounce, m_playTimer;
    QTreeWidget *m_tree = nullptr;
    QLineEdit *m_filter = nullptr;
    QLabel *m_path = nullptr, *m_summary = nullptr, *m_stats = nullptr;
    QComboBox *m_xAxis = nullptr, *m_yAxis = nullptr;
    QComboBox *m_files = nullptr;
    QLabel *m_fileCount = nullptr;
    struct FileView {
        QString path;
        QVector<quint64> indices;
        int x = -1, y = -1, tab = 0;
    };
    QHash<QString, FileView> m_fileViews;
    QVector<quint64> m_indices;
    int m_frameDimension = -1;
    ScientificPlot *m_plot = nullptr;
    SliceTableModel *m_tableModel = nullptr;
    QPlainTextEdit *m_text = nullptr;
    QTabWidget *m_tabs = nullptr;
    QSlider *m_slider = nullptr;
    QSpinBox *m_frame = nullptr;
    QToolButton *m_play = nullptr;
    QProgressBar *m_progress = nullptr;
    QAction *m_roiAction = nullptr;
    QAction *m_exportAction = nullptr, *m_closeFileAction = nullptr, *m_closeAllAction = nullptr,
            *m_cancelExportAction = nullptr;
    QFutureWatcher<QString> m_exportWatcher;
    std::atomic_bool m_exportCancel{false};
    bool m_exporting = false;
    QString m_exportError;
    QProgressBar *m_exportProgress = nullptr;
    SpectroscopyDialog *m_spectroscopy = nullptr;
    ImagingDialog *m_imaging = nullptr;
    BondWorkspace *m_bond = nullptr;
    QStackedWidget *m_workspaces = nullptr;
    QWidget *m_browser = nullptr, *m_roiReturn = nullptr, *m_frameControls = nullptr;
    QToolButton *m_backToBrowser = nullptr;
    QHash<QString, RoiDialog *> m_roiEditors;
    bool m_explorerVisible = true;
    Relaunch m_relaunch;
};
