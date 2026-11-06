#pragma once
#include "HdfFile.h"
#include <QFutureWatcher>
#include <QMainWindow>
#include <QTimer>
#include <QVector>

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

class MainWindow : public QMainWindow {
    Q_OBJECT
  public:
    explicit MainWindow(QWidget *parent = nullptr);
    ~MainWindow() override;
    void openFile(const QString &path);
    void selectPath(const QString &path);
    bool isReading() const {
        return m_reading || m_debounce.isActive() || m_pending;
    }
    const DataSlice &currentData() const {
        return m_slice;
    }
    void setSliceIndex(int dimension, int index);

  protected:
    void dragEnterEvent(QDragEnterEvent *) override;
    void dropEvent(QDropEvent *) override;
    void closeEvent(QCloseEvent *) override;

  private:
    void buildUi();
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
};
