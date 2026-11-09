#include "MainWindow.h"
#include "BondWorkspace.h"
#include "ImagingDialog.h"
#include "RoiDialog.h"
#include "ScientificPlot.h"
#include "SpectroscopyDialog.h"
#include "WorkflowUi.h"
#include <QSaveFile>
#include <QtConcurrent/QtConcurrentRun>
#include <QtWidgets>
#include <cmath>
#include <limits>

namespace {
constexpr int PathRole = Qt::UserRole, GroupRole = Qt::UserRole + 1, LoadedRole = Qt::UserRole + 2,
              DepthRole = Qt::UserRole + 3, FileRole = Qt::UserRole + 4;
QString sourceFile(QTreeWidgetItem *item) {
    if (!item)
        return {};
    while (item->parent())
        item = item->parent();
    return item->data(0, FileRole).toString();
}
QString value(double v) {
    return QString::number(v, 'g', 10);
}
} // namespace

class SliceTableModel : public QAbstractTableModel {
  public:
    explicit SliceTableModel(QObject *parent) : QAbstractTableModel(parent) {}
    DataSlice slice;
    void setSlice(const DataSlice &s) {
        beginResetModel();
        slice = s;
        endResetModel();
    }
    int rowCount(const QModelIndex &p = {}) const override {
        return p.isValid() ? 0 : (slice.yAxis >= 0 ? slice.height : slice.width);
    }
    int columnCount(const QModelIndex &p = {}) const override {
        return p.isValid() ? 0 : (slice.yAxis >= 0 ? slice.width : 1);
    }
    QVariant data(const QModelIndex &i, int role) const override {
        if (!i.isValid() || slice.values.isEmpty())
            return {};
        if (role == Qt::TextAlignmentRole)
            return int(Qt::AlignRight | Qt::AlignVCenter);
        if (role == Qt::DisplayRole) {
            int n = slice.yAxis >= 0 ? i.row() * slice.width + i.column() : i.row();
            return value(slice.values.value(n));
        }
        return {};
    }
    QVariant headerData(int s, Qt::Orientation o, int role) const override {
        if (role != Qt::DisplayRole)
            return {};
        if (slice.yAxis < 0 && o == Qt::Horizontal)
            return QStringLiteral("Value");
        int axis = o == Qt::Horizontal ? slice.xAxis : (slice.yAxis >= 0 ? slice.yAxis : slice.xAxis);
        return QString::number(quint64(s) * slice.strides.value(axis, 1));
    }
};

MainWindow::MainWindow(QWidget *parent) : QMainWindow(parent) {
    setObjectName("XRI.MainWindow");
    setWindowTitle(tr("XRI — 非共振非弹性散射成像数据处理工具"));
    resize(1440, 940);
    setMinimumSize(1080, 720);
    setAcceptDrops(true);
    buildUi();
    m_debounce.setSingleShot(true);
    m_debounce.setInterval(90);
    connect(&m_debounce, &QTimer::timeout, this, &MainWindow::startRead);
    connect(&m_watcher, &QFutureWatcher<DataSlice>::finished, this, &MainWindow::finishRead);
    connect(&m_exportWatcher, &QFutureWatcher<QString>::finished, this, [this] {
        m_exporting = false;
        m_exportError = m_exportWatcher.result();
        m_exportProgress->hide();
        syncFiles();
        statusBar()->showMessage(m_exportError.isEmpty() ? tr("HDF5 文件副本已导出") : m_exportError);
    });
    m_playTimer.setInterval(120);
    connect(&m_playTimer, &QTimer::timeout, this, [this] {
        if (!isReading()) {
            if (m_frame->value() >= m_frame->maximum()) {
                m_play->setChecked(false);
                return;
            }
            m_frame->setValue(m_frame->value() + 1);
        }
    });
    QSettings settings;
    restoreGeometry(settings.value("geometry").toByteArray());
    restoreState(settings.value("windowState").toByteArray());
}
MainWindow::~MainWindow() {
    m_exportCancel.store(true);
    m_exportWatcher.waitForFinished();
    m_watcher.waitForFinished();
}

void MainWindow::buildUi() {
    auto fileMenu = menuBar()->addMenu(tr("文件(&F)"));
    fileMenu->setObjectName("FileMenu");
    auto open = fileMenu->addAction(WorkflowUi::icon("folder-open"), tr("导入 HDF5 文件(&O)…"),
                                    QKeySequence::Open);
    open->setObjectName("ImportHdfFiles");
    connect(open, &QAction::triggered, this, [this] {
        auto files = QFileDialog::getOpenFileNames(
            this, tr("导入 HDF5 / NeXus 文件（可多选）"), QSettings().value("lastDirectory").toString(),
            "HDF5 / NeXus (*.nxs *.h5 *.hdf5 *.hdf *.nx5);;All files (*)");
        openFiles(files);
    });
    m_exportAction = fileMenu->addAction(WorkflowUi::icon("download"), tr("导出当前文件副本…"));
    m_exportAction->setObjectName("ExportHdfFile");
    connect(m_exportAction, &QAction::triggered, this, [this] {
        auto path = QFileDialog::getSaveFileName(this, tr("导出完整 HDF5 / NeXus 文件副本"),
                                                 QFileInfo(m_file).absolutePath() + "/" +
                                                     QFileInfo(m_file).completeBaseName() + "_copy." +
                                                     QFileInfo(m_file).suffix(),
                                                 "HDF5 / NeXus (*.nxs *.h5 *.hdf5);;All files (*)");
        if (!path.isEmpty())
            exportFile(path);
    });
    m_cancelExportAction = fileMenu->addAction(WorkflowUi::icon("square"), tr("取消文件导出"));
    m_cancelExportAction->setObjectName("CancelHdfExport");
    connect(m_cancelExportAction, &QAction::triggered, this, [this] { m_exportCancel.store(true); });
    fileMenu->addSeparator();
    m_closeFileAction =
        fileMenu->addAction(tr("关闭当前文件"), QKeySequence::Close, this, &MainWindow::closeCurrentFile);
    m_closeFileAction->setObjectName("CloseBrowserFile");
    m_closeAllAction = fileMenu->addAction(tr("关闭全部文件"), this, &MainWindow::closeAllFiles);
    m_closeAllAction->setObjectName("CloseAllBrowserFiles");
    fileMenu->addSeparator();
    fileMenu->addAction(tr("退出(&X)"), QKeySequence::Quit, this, &QWidget::close);
    auto viewMenu = menuBar()->addMenu(tr("数据浏览器(&B)"));
    viewMenu->setObjectName("DataBrowserMenu");
    auto showBrowser = viewMenu->addAction(tr("显示数据浏览器"));
    showBrowser->setObjectName("ShowDataBrowser");
    connect(showBrowser, &QAction::triggered, this,
            [this] { showWorkspace(m_browser); });
    viewMenu->addSeparator();
    const QStringList modes{tr("图像 / 曲线"), tr("数据表"), tr("文本 / 元数据")};
    for (int i = 0; i < modes.size(); ++i) {
        auto action = viewMenu->addAction(modes[i]);
        action->setObjectName(QString("BrowserMode%1").arg(i));
        connect(action, &QAction::triggered, this, [this, i] {
            showWorkspace(m_browser);
            m_tabs->setCurrentIndex(i);
        });
    }
    viewMenu->addSeparator();
    auto reset = viewMenu->addAction(WorkflowUi::icon("rotate-ccw"), tr("复位图表"));
    reset->setObjectName("ResetBrowserPlot");
    connect(reset, &QAction::triggered, this, [this] { m_plot->resetView(); });
    auto roiMenu = menuBar()->addMenu(tr("ROI(&R)"));
    m_roiAction = roiMenu->addAction(WorkflowUi::icon("scan"), tr("选取分析晶体 ROI…"));
    m_roiAction->setObjectName("SelectCrystalRois");
    m_roiAction->setEnabled(false);
    connect(m_roiAction, &QAction::triggered, this, [this] { openRoi(); });
    auto spectroscopyMenu = menuBar()->addMenu(tr("谱学(&S)"));
    auto spectrum = spectroscopyMenu->addAction(WorkflowUi::icon("chart-line"), tr("打开谱学处理工作区"));
    spectrum->setObjectName("SpectroscopyStep0");
    connect(spectrum, &QAction::triggered, this, [this] { openSpectroscopy(); });
    auto imagingMenu = menuBar()->addMenu(tr("成像(&I)"));
    auto imaging = imagingMenu->addAction(WorkflowUi::icon("image"), tr("打开成像处理工作区"));
    imaging->setObjectName("ImagingStep0");
    connect(imaging, &QAction::triggered, this, [this] { openImaging(); });
    auto chemicalBondMenu = menuBar()->addMenu(tr("化学键衬度分析(&O)"));
    chemicalBondMenu->addAction(WorkflowUi::icon("layers"), tr("打开四维衬度工作区"), this,
                                [this] { openBondImaging(); })
        ->setObjectName("OpenBondWorkspace");
    auto help = menuBar()->addMenu(tr("帮助(&H)"));
    help->addAction(tr("关于 XRI"), this, [this] {
        QMessageBox::about(
            this, tr("关于 XRI"),
            tr("<h3>XRI 0.2.1</h3>"
               "<p>非弹性非共振数据处理工具</p>"
               "<p>用于同步辐射 HDF5 / NeXus 数据浏览、分析晶体 ROI "
               "选取、谱学处理、三维成像与四维衬度分析。</p>"
               "<p>作者：Hu Shiqi<br>"
               "邮箱：<a href=\"mailto:husq@ihep.ac.cn\">husq@ihep.ac.cn</a><br>"
               "Copyright &copy; 2026 Hu Shiqi</p>"
               "<p>本项目源代码采用 <a href=\"https://opensource.org/license/mit\">MIT 开源协议</a>发布。<br>"
               "第三方组件遵循各自许可证，详见随附的 THIRD_PARTY.md。</p>"
               "<p>源代码：<a href=\"https://github.com/HCarrots/XRItools\">HCarrots/XRItools</a></p>"));
    });
    auto explorer = new QDockWidget(tr("文件列表与结构"), this);
    explorer->setObjectName("ExplorerDock");
    auto explorerBody = new QWidget;
    auto explorerLayout = new QVBoxLayout(explorerBody);
    explorerLayout->setContentsMargins(6, 6, 6, 6);
    explorerLayout->setSpacing(6);
    m_filter = new QLineEdit;
    m_filter->setPlaceholderText(tr("筛选已展开节点…"));
    m_filter->setClearButtonEnabled(true);
    explorerLayout->addWidget(m_filter);
    m_tree = new QTreeWidget;
    m_tree->setObjectName("BrowserFileTree");
    m_tree->setHeaderLabels({tr("名称"), tr("维度")});
    m_tree->setAlternatingRowColors(true);
    m_tree->setUniformRowHeights(true);
    m_tree->setAnimated(false);
    m_tree->setColumnWidth(0, 220);
    m_tree->setMinimumWidth(270);
    m_tree->setSelectionMode(QAbstractItemView::SingleSelection);
    explorerLayout->addWidget(m_tree);
    auto hint = new QLabel(tr("双击组以展开 · 选择数据集以预览"));
    hint->setWordWrap(true);
    explorerLayout->addWidget(hint);
    explorer->setWidget(explorerBody);
    addDockWidget(Qt::LeftDockWidgetArea, explorer);
    viewMenu->addAction(explorer->toggleViewAction());
    connect(m_filter, &QLineEdit::textChanged, this, &MainWindow::updateFilter);
    connect(m_tree, &QTreeWidget::itemExpanded, this, [this](QTreeWidgetItem *item) {
        if (!item->data(0, LoadedRole).toBool()) {
            populate(item, item->data(0, PathRole).toString(), item->data(0, DepthRole).toInt() + 1);
            updateFilter(m_filter->text());
        }
    });
    connect(m_tree, &QTreeWidget::currentItemChanged, this, [this](QTreeWidgetItem *item, QTreeWidgetItem *) {
        if (item)
            selectItem(item);
    });

    auto center = new QWidget;
    auto layout = new QVBoxLayout(center);
    layout->setContentsMargins(6, 4, 6, 4);
    layout->setSpacing(4);
    auto fileRow = new QHBoxLayout;
    fileRow->addWidget(new QLabel(tr("当前文件")));
    m_files = new QComboBox;
    m_files->setObjectName("BrowserFiles");
    m_files->setSizeAdjustPolicy(QComboBox::AdjustToMinimumContentsLengthWithIcon);
    m_files->setMinimumContentsLength(35);
    fileRow->addWidget(m_files, 1);
    m_fileCount = new QLabel;
    fileRow->addWidget(m_fileCount);
    layout->addLayout(fileRow);
    connect(m_files, &QComboBox::currentIndexChanged, this, [this](int index) {
        if (index >= 0)
            activateFile(m_files->itemData(index).toString());
    });
    m_path = new QLabel(tr("未选择数据集"));
    QFont pathFont = m_path->font();
    pathFont.setBold(true);
    m_path->setFont(pathFont);
    m_path->setTextInteractionFlags(Qt::TextSelectableByMouse);
    m_path->setSizePolicy(QSizePolicy::Ignored, QSizePolicy::Preferred);
    m_summary = new QLabel(tr("支持标量、字符串、一维信号和多维探测器数据"));
    m_summary->setTextInteractionFlags(Qt::TextSelectableByMouse);
    auto dataset = new QHBoxLayout;
    dataset->addWidget(m_path, 1);
    dataset->addWidget(m_summary);
    layout->addLayout(dataset);
    auto plotToolsBody = new QWidget;
    auto plotTools = new QHBoxLayout(plotToolsBody);
    plotTools->setContentsMargins(0, 0, 0, 0);
    plotTools->addWidget(new QLabel(tr("X 轴")));
    m_xAxis = new QComboBox;
    m_xAxis->setMinimumWidth(100);
    plotTools->addWidget(m_xAxis);
    plotTools->addSpacing(8);
    plotTools->addWidget(new QLabel(tr("Y 轴")));
    m_yAxis = new QComboBox;
    m_yAxis->setMinimumWidth(100);
    plotTools->addWidget(m_yAxis);
    plotTools->addStretch();
    auto palette = new QComboBox;
    palette->addItems({"Viridis", "Inferno", tr("灰度")});
    plotTools->addWidget(palette);
    auto log = new QCheckBox("signed log1p");
    log->setToolTip(tr("符号保留的 log(1+|值|)，支持零及负数"));
    plotTools->addWidget(log);
    auto aspect = new QCheckBox(tr("等比例"));
    aspect->setChecked(true);
    plotTools->addWidget(aspect);
    layout->addWidget(plotToolsBody);
    m_tabs = new QTabWidget;
    m_tabs->setObjectName("BrowserModes");
    m_tabs->setDocumentMode(true);
    m_plot = new ScientificPlot;
    m_plot->setObjectName("ScientificPlot");
    m_tabs->addTab(m_plot, tr("图像 / 曲线"));
    auto table = new QTableView;
    m_tableModel = new SliceTableModel(this);
    table->setModel(m_tableModel);
    table->setAlternatingRowColors(true);
    table->horizontalHeader()->setDefaultSectionSize(105);
    table->verticalHeader()->setDefaultSectionSize(24);
    m_tabs->addTab(table, tr("数据表"));
    m_text = new QPlainTextEdit;
    m_text->setReadOnly(true);
    m_text->setFont(QFont("Consolas", 10));
    m_tabs->addTab(m_text, tr("文本 / 元数据"));
    connect(m_tabs, &QTabWidget::currentChanged, plotToolsBody,
            [plotToolsBody](int tab) { plotToolsBody->setVisible(tab == 0); });
    layout->addWidget(m_tabs, 1);
    m_frameControls = new QWidget;
    auto frames = new QHBoxLayout(m_frameControls);
    frames->setContentsMargins(0, 0, 0, 0);
    m_play = new QToolButton;
    m_play->setText(tr("播放"));
    const auto playIcon = WorkflowUi::icon("play");
    const auto pauseIcon = WorkflowUi::icon("pause");
    m_play->setIcon(playIcon);
    m_play->setIconSize(QSize(16, 16));
    m_play->setToolButtonStyle(Qt::ToolButtonTextBesideIcon);
    m_play->setCheckable(true);
    frames->addWidget(m_play);
    frames->addWidget(new QLabel(tr("帧 / 首个固定维度")));
    m_slider = new QSlider(Qt::Horizontal);
    m_slider->setObjectName("FrameSlider");
    frames->addWidget(m_slider, 1);
    m_frame = new QSpinBox;
    m_frame->setObjectName("FrameIndex");
    m_frame->setMinimumWidth(100);
    frames->addWidget(m_frame);
    layout->addWidget(m_frameControls);
    m_frameControls->hide();
    m_stats = new QLabel;
    m_stats->setTextInteractionFlags(Qt::TextSelectableByMouse);
    m_stats->setWordWrap(true);
    layout->addWidget(m_stats);
    m_browser = center;
    m_backToBrowser = new QToolButton;
    m_backToBrowser->setText(tr("返回数据浏览器"));
    m_backToBrowser->setAutoRaise(true);
    m_backToBrowser->setIcon(WorkflowUi::icon("chevron-left"));
    m_backToBrowser->setIconSize(QSize(16, 16));
    m_backToBrowser->setToolButtonStyle(Qt::ToolButtonTextBesideIcon);
    m_backToBrowser->setObjectName("BackToBrowser");
    menuBar()->setCornerWidget(m_backToBrowser);
    m_backToBrowser->hide();
    m_workspaces = new QStackedWidget;
    m_workspaces->setObjectName("ProcessingWorkspaces");
    m_workspaces->addWidget(center);
    setCentralWidget(m_workspaces);
    connect(m_backToBrowser, &QToolButton::clicked, this,
            [this] { showWorkspace(m_browser); });
    connect(palette, &QComboBox::currentIndexChanged, m_plot, &ScientificPlot::setPalette);
    connect(log, &QCheckBox::toggled, m_plot, &ScientificPlot::setLogarithmic);
    connect(aspect, &QCheckBox::toggled, m_plot, &ScientificPlot::setAspectLocked);
    connect(m_xAxis, &QComboBox::currentIndexChanged, this, [this] {
        if (!m_configuring) {
            if (m_xAxis->currentData() == m_yAxis->currentData()) {
                QSignalBlocker b(m_yAxis);
                m_yAxis->setCurrentIndex(0);
            }
            buildSliceControls();
            scheduleRead();
        }
    });
    connect(m_yAxis, &QComboBox::currentIndexChanged, this, [this] {
        if (!m_configuring) {
            if (m_yAxis->currentData() == m_xAxis->currentData() && m_yAxis->currentData().toInt() >= 0) {
                int other =
                    m_node.dims.size() > 1 ? (m_yAxis->currentData().toInt() + 1) % m_node.dims.size() : -1;
                QSignalBlocker b(m_xAxis);
                m_xAxis->setCurrentIndex(m_xAxis->findData(other));
                if (other < 0) {
                    QSignalBlocker c(m_yAxis);
                    m_yAxis->setCurrentIndex(0);
                    m_xAxis->setCurrentIndex(0);
                }
            }
            buildSliceControls();
            scheduleRead();
        }
    });
    connect(m_slider, &QSlider::valueChanged, m_frame, &QSpinBox::setValue);
    connect(m_frame, &QSpinBox::valueChanged, m_slider, &QSlider::setValue);
    connect(m_frame, &QSpinBox::valueChanged, this, [this](int v) {
        if (m_frameDimension >= 0) {
            m_indices[m_frameDimension] = quint64(v);
            scheduleRead();
        }
    });
    connect(m_play, &QToolButton::toggled, this, [this, playIcon, pauseIcon](bool play) {
        m_play->setText(play ? tr("暂停") : tr("播放"));
        m_play->setIcon(play ? pauseIcon : playIcon);
        if (play)
            m_playTimer.start();
        else
            m_playTimer.stop();
    });
    connect(m_plot, &ScientificPlot::cursorChanged, this,
            [this](const QString &s) { statusBar()->showMessage(s); });
    m_progress = new QProgressBar;
    m_progress->setRange(0, 0);
    m_progress->setMaximumWidth(140);
    m_progress->hide();
    statusBar()->addPermanentWidget(m_progress);
    m_exportProgress = new QProgressBar;
    m_exportProgress->setRange(0, 100);
    m_exportProgress->setMaximumWidth(170);
    m_exportProgress->setFormat(tr("文件导出 %p%"));
    m_exportProgress->hide();
    statusBar()->addPermanentWidget(m_exportProgress);
    statusBar()->addPermanentWidget(new QLabel(tr("只读 · HDF5 / NeXus")));
    statusBar()->showMessage(tr("就绪"));
    resizeDocks({explorer}, {310}, Qt::Horizontal);
    syncFiles();
}

void MainWindow::showWorkspace(QWidget *page, const QString &title) {
    m_play->setChecked(false);
    if (auto editor = qobject_cast<RoiDialog *>(m_workspaces->currentWidget()); editor && editor != page) {
        const auto file = m_roiEditors.key(editor);
        const auto path =
            QStandardPaths::writableLocation(QStandardPaths::AppLocalDataLocation) + "/roi-drafts/" +
            QCryptographicHash::hash(file.toUtf8(), QCryptographicHash::Sha256).toHex() + ".json";
        if (m_spectroscopy)
            m_spectroscopy->useRoiDraft(file, path);
        if (m_imaging)
            m_imaging->useRoiDraft(file, path);
        if (m_bond)
            m_bond->useRoiDraft(file, path);
    }
    auto explorer = findChild<QDockWidget *>("ExplorerDock");
    if (m_workspaces->currentWidget() == m_browser && page != m_browser)
        m_explorerVisible = !explorer->isHidden();
    explorer->setVisible(page == m_browser && m_explorerVisible);
    m_workspaces->setCurrentWidget(page);
    page->show();
    if (page == m_browser)
        updateBrowserTitle();
    else
        setWindowTitle(tr("%1 — XRI").arg(title));
    m_backToBrowser->setVisible(page != m_browser);
    m_roiAction->setEnabled(page != m_browser || !m_file.isEmpty());
    show();
}

void MainWindow::updateBrowserTitle() {
    if (m_workspaces->currentWidget() == m_browser)
        setWindowTitle(m_file.isEmpty() ? tr("XRI — 非共振非弹性散射成像数据处理工具")
                                       : QFileInfo(m_file).fileName() + " [" + formatLabel() + tr("] — XRI"));
}
SpectroscopyDialog *MainWindow::openSpectroscopy(int step) {
    if (!m_spectroscopy) {
        m_spectroscopy = new SpectroscopyDialog(m_file, this);
        m_spectroscopy->setWindowFlags(Qt::Widget);
        m_spectroscopy->setMinimumSize(0, 0);
        m_spectroscopy->findChild<QPushButton *>("CloseWorkflow")->hide();
        m_workspaces->addWidget(m_spectroscopy);
        connect(m_spectroscopy, &SpectroscopyDialog::roiRequested, this, [this](const QString &file) {
            openRoi(file, m_spectroscopy->findChild<QLineEdit *>("SpectrumRoiScheme")->text());
        });
    }
    m_spectroscopy->setAvailableFiles(openedFiles());
    if (step >= 0)
        m_spectroscopy->selectStep(step);
    showWorkspace(m_spectroscopy, tr("谱学 · scan → 弹性标定 → XRS 提取 → 叠加与保存"));
    return m_spectroscopy;
}
ImagingDialog *MainWindow::openImaging(int step) {
    if (!m_imaging) {
        m_imaging = new ImagingDialog(m_file, this);
        m_imaging->setWindowFlags(Qt::Widget);
        m_imaging->setMinimumSize(0, 0);
        m_imaging->findChild<QPushButton *>("CloseWorkflow")->hide();
        m_workspaces->addWidget(m_imaging);
        connect(m_imaging, &ImagingDialog::roiRequested, this, [this](const QString &file) {
            openRoi(file, m_imaging->findChild<QLineEdit *>("input.roi_table")->text());
        });
        connect(m_imaging, &ImagingDialog::openFileRequested, this, [this](const QString &path) {
            openFile(path);
            selectPath("/reinterp_insect_ck");
            showWorkspace(m_browser);
        });
    }
    m_imaging->setAvailableFiles(openedFiles());
    if (step >= 0)
        m_imaging->selectStep(step);
    showWorkspace(m_imaging, tr("成像 · grid_scan → 输入检查 → 三维重建 → 结果与诊断"));
    return m_imaging;
}
BondWorkspace *MainWindow::openBondImaging(int step) {
    if (!m_bond) {
        m_bond = new BondWorkspace(m_file, this);
        m_workspaces->addWidget(m_bond);
        connect(m_bond, &BondWorkspace::roiRequested, this, [this](const QString &file) {
            openRoi(file, m_bond->findChild<QLineEdit *>("input.roi_table")->text());
        });
    }
    m_bond->setAvailableFiles(openedFiles());
    if (step >= 0)
        m_bond->selectStep(step);
    showWorkspace(m_bond, tr("化学键衬度分析 · 四维重建 → 体素谱 → 衬度图"));
    return m_bond;
}
RoiDialog *MainWindow::openRoi(const QString &source, const QString &scheme) {
    if (source.isEmpty())
        if (auto current = qobject_cast<RoiDialog *>(m_workspaces->currentWidget()))
            return current;
    QString roiScheme = scheme;
    QString selected = source.isEmpty() ? m_file : source;
    if (source.isEmpty() && m_workspaces->currentWidget() == m_spectroscopy) {
        selected = m_spectroscopy->findChild<QLineEdit *>("ElasticFile")->text();
        roiScheme = m_spectroscopy->findChild<QLineEdit *>("SpectrumRoiScheme")->text();
    }
    if (source.isEmpty() && m_workspaces->currentWidget() == m_imaging) {
        selected = m_imaging->findChild<QLineEdit *>("input.sample_nxs")->text();
        roiScheme = m_imaging->findChild<QLineEdit *>("input.roi_table")->text();
    }
    if (source.isEmpty() && m_workspaces->currentWidget() == m_bond) {
        selected = m_bond->findChild<QLineEdit *>("input.sample_nxs")->text();
        roiScheme = m_bond->findChild<QLineEdit *>("input.roi_table")->text();
    }
    const QString file = QFileInfo(selected).canonicalFilePath();
    if (file.isEmpty()) {
        statusBar()->showMessage(tr("请先为当前处理步骤选择有效的数据文件。"));
        return nullptr;
    }
    if (!qobject_cast<RoiDialog *>(m_workspaces->currentWidget()))
        m_roiReturn = m_workspaces->currentWidget();
    auto editor = m_roiEditors.value(file);
    if (!editor) {
        editor = new RoiDialog(file, tr("共用 ROI"), this);
        editor->setWindowFlags(Qt::Widget);
        editor->setMinimumSize(0, 0);
        m_workspaces->addWidget(editor);
        m_roiEditors.insert(file, editor);
        connect(editor, &QDialog::finished, this, [this](int) {
            if (m_roiReturn == m_spectroscopy)
                openSpectroscopy(0);
            else if (m_roiReturn == m_imaging)
                openImaging(0);
            else if (m_roiReturn == m_bond)
                openBondImaging(0);
            else
                showWorkspace(m_browser);
        });
    }
    if (roiScheme.endsWith(".json", Qt::CaseInsensitive) && QFileInfo::exists(roiScheme))
        editor->loadFrom(roiScheme);
    showWorkspace(editor, tr("共用 ROI · 75 个分析晶体 — %1").arg(QFileInfo(file).fileName()));
    return editor;
}

void MainWindow::populate(QTreeWidgetItem *parent, const QString &path, int depth) {
    if (depth > 32) {
        parent->setToolTip(0, tr("已达到 32 层浏览深度，避免循环链接。"));
        return;
    }
    try {
        auto nodes = HdfFile::children(sourceFile(parent), path);
        parent->setData(0, LoadedRole, true);
        for (const auto &n : nodes) {
            auto item = new QTreeWidgetItem(parent, {n.name, n.shape});
            item->setData(0, PathRole, n.path);
            item->setData(0, GroupRole, n.group);
            item->setData(0, DepthRole, depth);
            item->setToolTip(0, n.path + "\n" + n.type);
            item->setIcon(0, style()->standardIcon(n.group ? QStyle::SP_DirIcon : QStyle::SP_FileIcon));
            if (n.group)
                item->setChildIndicatorPolicy(QTreeWidgetItem::ShowIndicator);
        }
    } catch (const std::exception &e) {
        statusBar()->showMessage(QString::fromUtf8(e.what()));
    }
}
void MainWindow::openFile(const QString &path) {
    auto file = QFileInfo(path).canonicalFilePath();
    try {
        HdfFile::inspect(path, "/");
    } catch (const std::exception &e) {
        QMessageBox::critical(this, tr("打开失败"), QString::fromUtf8(e.what()));
        return;
    }
    if (auto existing = fileRoot(file)) {
        rememberFile();
        auto name = sourceFile(existing);
        m_fileViews.remove(name);
        m_file.clear();
        m_node = {};
        m_indices.clear();
        activateFile(name);
        return;
    }
    QFileInfo info(file);
    QSettings().setValue("lastDirectory", info.absolutePath());
    auto format = HdfFile::formatFromFileName(file);
    auto label = format == ScanFormat::Scan       ? tr("scan · 一维扫描")
                 : format == ScanFormat::GridScan ? tr("grid_scan · 网格扫描")
                                                  : tr("通用 HDF5 / NeXus");
    auto root = new QTreeWidgetItem(m_tree, {info.fileName(), label});
    root->setToolTip(0, file + "\n" + label);
    root->setToolTip(1, label);
    root->setData(0, PathRole, "/");
    root->setData(0, FileRole, file);
    root->setData(0, GroupRole, true);
    root->setData(0, DepthRole, 0);
    populate(root, "/", 1);
    root->setExpanded(true);
    activateFile(file);
}
void MainWindow::openFiles(const QStringList &paths) {
    if (!paths.isEmpty())
        showWorkspace(m_browser);
    for (const auto &path : paths)
        openFile(path);
}
QStringList MainWindow::openedFiles() const {
    QStringList files;
    for (int i = 0; i < m_tree->topLevelItemCount(); ++i)
        files.append(sourceFile(m_tree->topLevelItem(i)));
    return files;
}
QTreeWidgetItem *MainWindow::fileRoot(const QString &file) const {
    auto name = QFileInfo(file).canonicalFilePath();
    if (name.isEmpty())
        name = QFileInfo(file).absoluteFilePath();
    for (int i = 0; i < m_tree->topLevelItemCount(); ++i) {
        auto root = m_tree->topLevelItem(i);
        if (sourceFile(root).compare(name, Qt::CaseInsensitive) == 0)
            return root;
    }
    return nullptr;
}
void MainWindow::rememberFile() {
    if (!m_file.isEmpty() && !m_node.path.isEmpty() && !m_node.group)
        m_fileViews[m_file] = {m_node.path, m_indices, m_xAxis->currentData().toInt(),
                               m_yAxis->currentData().toInt(), m_tabs->currentIndex()};
}
void MainWindow::syncFiles() {
    if (m_spectroscopy)
        m_spectroscopy->setAvailableFiles(openedFiles());
    if (m_imaging)
        m_imaging->setAvailableFiles(openedFiles());
    if (m_bond)
        m_bond->setAvailableFiles(openedFiles());
    QSignalBlocker block(m_files);
    m_files->clear();
    for (const auto &file : openedFiles()) {
        m_files->addItem(QFileInfo(file).fileName(), file);
        m_files->setItemData(m_files->count() - 1, file, Qt::ToolTipRole);
    }
    m_files->setCurrentIndex(m_files->findData(m_file));
    m_files->setEnabled(m_files->count() > 0);
    m_fileCount->setText(tr("共 %1 个文件").arg(m_files->count()));
    m_roiAction->setEnabled(!m_file.isEmpty());
    m_exportAction->setEnabled(!m_file.isEmpty() && !m_exporting);
    m_closeFileAction->setEnabled(!m_file.isEmpty());
    m_closeAllAction->setEnabled(m_files->count() > 0);
    m_cancelExportAction->setEnabled(m_exporting);
}
void MainWindow::activateFile(const QString &file) {
    auto root = fileRoot(file);
    if (!root)
        return;
    for (int i = 0; i < m_tree->topLevelItemCount(); ++i) {
        auto other = m_tree->topLevelItem(i);
        if (other != root)
            other->setExpanded(false);
    }
    root->setExpanded(true);
    auto name = sourceFile(root);
    auto desired = m_fileViews.value(name).path;
    QTreeWidgetItem *target = nullptr;
    if (!desired.isEmpty())
        target = itemAtPath(root, desired);
    if (!target)
        target = itemAtPath(root, "/entry/data/D_LAMBDA");
    if (!target)
        target = itemAtPath(root, "/reinterp_insect_ck");
    if (!target)
        target = root;
    m_filter->clear();
    if (m_tree->currentItem() == target)
        selectItem(target);
    else
        m_tree->setCurrentItem(target);
    m_tree->scrollToItem(target);
}
QTreeWidgetItem *MainWindow::itemAtPath(QTreeWidgetItem *root, const QString &path) {
    auto current = root;
    if (!current)
        return nullptr;
    QString prefix;
    for (const auto &part : path.split('/', Qt::SkipEmptyParts)) {
        prefix += '/' + part;
        if (!current->data(0, LoadedRole).toBool())
            populate(current, current->data(0, PathRole).toString(), current->data(0, DepthRole).toInt() + 1);
        QTreeWidgetItem *next = nullptr;
        for (int i = 0; i < current->childCount(); ++i)
            if (current->child(i)->data(0, PathRole).toString() == prefix) {
                next = current->child(i);
                break;
            }
        if (!next)
            return nullptr;
        current->setExpanded(true);
        current = next;
    }
    return current;
}
void MainWindow::selectPath(const QString &path) {
    auto item = itemAtPath(fileRoot(m_file), path);
    if (!item)
        return;
    if (m_tree->currentItem() == item)
        selectItem(item);
    else
        m_tree->setCurrentItem(item);
    m_tree->scrollToItem(item);
}
void MainWindow::clearPreview() {
    ++m_generation;
    m_pending = false;
    m_debounce.stop();
    m_play->setChecked(false);
    m_file.clear();
    m_format = ScanFormat::Generic;
    m_node = {};
    m_slice = {};
    m_indices.clear();
    m_plot->setData({});
    m_tableModel->setSlice({});
    m_text->clear();
    m_stats->clear();
    m_path->setText(tr("未选择数据集"));
    m_path->setToolTip({});
    m_summary->setText(tr("通过“文件 → 导入 HDF5 文件”添加一个或多个文件"));
    configureAxes();
    updateBrowserTitle();
    syncFiles();
}
void MainWindow::closeCurrentFile() {
    auto root = fileRoot(m_file);
    if (!root)
        return;
    auto file = sourceFile(root);
    int index = m_tree->indexOfTopLevelItem(root);
    {
        QSignalBlocker block(m_tree);
        delete m_tree->takeTopLevelItem(index);
    }
    m_fileViews.remove(file);
    clearPreview();
    if (m_tree->topLevelItemCount())
        activateFile(sourceFile(m_tree->topLevelItem(std::min(index, m_tree->topLevelItemCount() - 1))));
}
void MainWindow::closeAllFiles() {
    {
        QSignalBlocker block(m_tree);
        m_tree->clear();
    }
    m_fileViews.clear();
    clearPreview();
}
bool MainWindow::exportFile(const QString &destination) {
    if (m_file.isEmpty() || m_exporting)
        return false;
    auto name = QFileInfo(destination).canonicalFilePath();
    if (name.isEmpty())
        name = QFileInfo(destination).absoluteFilePath();
    for (const auto &file : openedFiles())
        if (name.compare(file, Qt::CaseInsensitive) == 0) {
            m_exportError = tr("导出路径不能是任何已导入的原始文件，请选择新文件名。");
            statusBar()->showMessage(m_exportError);
            return false;
        }
    const auto source = m_file;
    m_exportCancel.store(false);
    m_exportError.clear();
    m_exporting = true;
    m_exportProgress->setValue(0);
    m_exportProgress->show();
    syncFiles();
    m_exportWatcher.setFuture(QtConcurrent::run([this, source, destination] {
        QFile input(source);
        if (!input.open(QIODevice::ReadOnly))
            return input.errorString();
        QSaveFile output(destination);
        if (!output.open(QIODevice::WriteOnly))
            return output.errorString();
        qint64 total = input.size(), done = 0;
        int previous = -1;
        while (done < total) {
            if (m_exportCancel.load())
                return tr("文件导出已取消。");
            auto data = input.read(std::min<qint64>(8 * 1024 * 1024, total - done));
            if (data.isEmpty())
                return tr("源文件读取中断，导出未完成。");
            if (output.write(data) != data.size())
                return output.errorString();
            done += data.size();
            int percent = int(done * 100 / std::max<qint64>(total, 1));
            if (percent != previous) {
                previous = percent;
                QMetaObject::invokeMethod(
                    this,
                    [this, percent] {
                        if (m_exporting)
                            m_exportProgress->setValue(percent);
                    },
                    Qt::QueuedConnection);
            }
        }
        if (m_exportCancel.load())
            return tr("文件导出已取消。");
        if (input.size() != total)
            return tr("源文件大小在导出期间发生变化，请等待采集完成后重试。");
        if (!output.commit())
            return output.errorString();
        return QString();
    }));
    return true;
}
void MainWindow::selectItem(QTreeWidgetItem *item) {
    rememberFile();
    const auto file = sourceFile(item);
    if (file != m_file) {
        m_indices.clear();
        m_node = {};
    }
    m_file = file;
    m_format = HdfFile::formatFromFileName(file);
    ++m_generation;
    m_pending = false;
    m_debounce.stop();
    m_play->setChecked(false);
    m_slice = {};
    m_plot->setData({});
    m_tableModel->setSlice({});
    m_stats->clear();
    m_text->clear();
    syncFiles();
    updateBrowserTitle();
    try {
        m_node = HdfFile::inspect(m_file, item->data(0, PathRole).toString());
    } catch (const std::exception &e) {
        statusBar()->showMessage(QString::fromUtf8(e.what()));
        return;
    }
    m_path->setText(m_node.path);
    m_path->setToolTip(m_node.path);
    m_summary->setText(formatLabel() + "    |    " + m_node.type + "    " + m_node.shape);
    QString text = "File: " + m_file + "\n" + m_node.path + "\nScan format: " + formatLabel() +
                   "\nType: " + m_node.type + "\nShape: " + m_node.shape + "\n";
    for (const auto &attr : m_node.attributes)
        text += "\n@" + attr.first + " = " + attr.second;
    m_text->setPlainText(text);
    m_tabs->setCurrentIndex(m_node.numeric ? 0 : 2);
    configureAxes();
    const auto saved = m_fileViews.value(file);
    if (saved.path == m_node.path && !m_node.group) {
        m_configuring = true;
        int x = m_xAxis->findData(saved.x), y = m_yAxis->findData(saved.y);
        if (x >= 0 && y >= 0 && (saved.y < 0 || saved.y != saved.x)) {
            m_xAxis->setCurrentIndex(x);
            m_yAxis->setCurrentIndex(y);
        }
        m_indices = saved.indices;
        buildSliceControls();
        m_configuring = false;
        m_tabs->setCurrentIndex(std::clamp(saved.tab, 0, 2));
    }
    if (!m_node.group)
        scheduleRead();
}
void MainWindow::configureAxes() {
    m_configuring = true;
    m_xAxis->clear();
    m_yAxis->clear();
    m_yAxis->addItem(tr("无 · 一维曲线"), -1);
    for (int i = 0; i < m_node.dims.size(); ++i) {
        auto label = QString("Dim %1 (%2)").arg(i).arg(m_node.dims[i]);
        m_xAxis->addItem(label, i);
        m_yAxis->addItem(label, i);
    }
    int rank = int(m_node.dims.size());
    m_xAxis->setCurrentIndex(rank - 1);
    m_yAxis->setCurrentIndex(rank >= 2 ? rank - 1 : 0);
    m_xAxis->setEnabled(m_node.numeric && rank > 0);
    m_yAxis->setEnabled(m_node.numeric && rank > 1);
    buildSliceControls();
    m_configuring = false;
}
void MainWindow::buildSliceControls() {
    m_play->setChecked(false);
    m_indices.resize(m_node.dims.size());
    m_frameDimension = -1;
    const int x = m_xAxis->currentData().toInt(), y = m_yAxis->currentData().toInt();
    for (int i = 0; i < m_node.dims.size(); ++i) {
        const quint64 maximum = m_node.dims[i] ? m_node.dims[i] - 1 : 0;
        m_indices[i] = std::min(m_indices[i], maximum);
        const bool fixed = m_node.type.contains("string") || (i != x && i != y && m_node.numeric);
        if (fixed && m_frameDimension < 0)
            m_frameDimension = i;
    }
    QSignalBlocker frameBlocker(m_frame), sliderBlocker(m_slider);
    const bool enabled = m_frameDimension >= 0;
    const int maximum =
        enabled ? int(std::min<quint64>(m_node.dims[m_frameDimension] ? m_node.dims[m_frameDimension] - 1 : 0,
                                        std::numeric_limits<int>::max()))
                : 0;
    const int index = enabled ? int(std::min<quint64>(m_indices[m_frameDimension], maximum)) : 0;
    m_frame->setRange(0, maximum);
    m_slider->setRange(0, maximum);
    m_frame->setValue(index);
    m_slider->setValue(index);
    m_frame->setEnabled(enabled);
    m_slider->setEnabled(enabled);
    m_play->setEnabled(enabled);
    m_frameControls->setVisible(enabled);
}
void MainWindow::scheduleRead() {
    if (m_file.isEmpty() || m_node.group)
        return;
    ++m_generation;
    m_pending = true;
    m_debounce.start();
}
void MainWindow::startRead() {
    if (m_reading)
        return;
    if (!m_pending || m_file.isEmpty() || m_node.group)
        return;
    m_pending = false;
    m_readGeneration = m_generation;
    SliceRequest request;
    request.file = m_file;
    request.path = m_node.path;
    request.xAxis = m_node.dims.isEmpty() ? -1 : m_xAxis->currentData().toInt();
    request.yAxis = m_yAxis->currentData().toInt();
    request.indices = m_indices;
    m_progress->show();
    statusBar()->showMessage(tr("正在读取切片…"));
    m_reading = true;
    m_watcher.setFuture(QtConcurrent::run([request] { return HdfFile::read(request); }));
}
void MainWindow::finishRead() {
    m_reading = false;
    m_progress->hide();
    if (m_readGeneration == m_generation) {
        m_slice = m_watcher.result();
        m_plot->setData(m_slice);
        m_tableModel->setSlice(m_slice);
        if (!m_slice.error.isEmpty()) {
            statusBar()->showMessage(m_slice.error);
            m_stats->setText(tr("读取失败：") + m_slice.error);
        } else if (m_slice.values.isEmpty()) {
            m_text->setPlainText(m_node.path + "\nScan format: " + formatLabel() + "\nType: " + m_node.type +
                                 "\nShape: " + m_node.shape + "\n\n" + m_slice.text);
            m_tabs->setCurrentIndex(2);
            statusBar()->showMessage(tr("文本预览已加载"));
        } else {
            m_stats->setText(QString("%1  ·  %2 × %3  ·  min %4   max %5   mean %6   Σ %7  ·  finite %8/%9")
                                 .arg(m_slice.sampled ? tr("步进采样预览") : tr("完整切片"))
                                 .arg(m_slice.width)
                                 .arg(m_slice.height)
                                 .arg(value(m_slice.minimum), value(m_slice.maximum), value(m_slice.mean),
                                      value(m_slice.sum))
                                 .arg(m_slice.finiteCount)
                                 .arg(m_slice.values.size()));
            statusBar()->showMessage(tr("切片已加载 · 滚轮缩放，双击复位"));
        }
    }
    if (m_pending && !m_debounce.isActive())
        startRead();
}
void MainWindow::setSliceIndex(int dimension, int index) {
    if (dimension < 0 || dimension >= m_indices.size())
        return;
    const int x = m_xAxis->currentData().toInt(), y = m_yAxis->currentData().toInt();
    if (!m_node.type.contains("string") && (!m_node.numeric || dimension == x || dimension == y))
        return;
    const quint64 maximum = m_node.dims[dimension] ? m_node.dims[dimension] - 1 : 0;
    const quint64 bounded = std::min<quint64>(quint64(std::max(0, index)), maximum);
    if (dimension == m_frameDimension)
        m_frame->setValue(int(std::min<quint64>(bounded, std::numeric_limits<int>::max())));
    else if (m_indices[dimension] != bounded) {
        m_indices[dimension] = bounded;
        scheduleRead();
    }
}
QString MainWindow::formatLabel() const {
    switch (m_format) {
    case ScanFormat::Scan:
        return tr("scan · 一维扫描");
    case ScanFormat::GridScan:
        return tr("grid_scan · 网格扫描");
    default:
        return tr("通用 HDF5 / NeXus");
    }
}

void MainWindow::updateFilter(const QString &text) {
    std::function<bool(QTreeWidgetItem *)> filter = [&](QTreeWidgetItem *item) {
        bool show = text.isEmpty() || item->text(0).contains(text, Qt::CaseInsensitive);
        for (int i = 0; i < item->childCount(); ++i)
            show = filter(item->child(i)) || show;
        item->setHidden(!show);
        return show;
    };
    for (int i = 0; i < m_tree->topLevelItemCount(); ++i)
        filter(m_tree->topLevelItem(i));
}
void MainWindow::dragEnterEvent(QDragEnterEvent *e) {
    if (e->mimeData()->hasUrls() && e->mimeData()->urls().first().isLocalFile())
        e->acceptProposedAction();
}
void MainWindow::dropEvent(QDropEvent *e) {
    if (e->mimeData()->hasUrls()) {
        QStringList paths;
        for (const auto &url : e->mimeData()->urls())
            if (url.isLocalFile() && QFileInfo(url.toLocalFile()).isFile())
                paths.append(url.toLocalFile());
        openFiles(paths);
        e->acceptProposedAction();
    }
}
void MainWindow::closeEvent(QCloseEvent *e) {
    m_playTimer.stop();
    QSettings settings;
    settings.setValue("geometry", saveGeometry());
    settings.setValue("windowState", saveState());
    QMainWindow::closeEvent(e);
}
