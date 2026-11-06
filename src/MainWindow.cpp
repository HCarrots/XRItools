#include "MainWindow.h"
#include "ScientificPlot.h"
#include <QtConcurrent/QtConcurrentRun>
#include <QtWidgets>
#include <cmath>
#include <limits>

namespace {
constexpr int PathRole = Qt::UserRole, GroupRole = Qt::UserRole + 1, LoadedRole = Qt::UserRole + 2,
              DepthRole = Qt::UserRole + 3;
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
    m_watcher.waitForFinished();
}

void MainWindow::buildUi() {
    auto fileMenu = menuBar()->addMenu(tr("文件(&F)"));
    auto open = fileMenu->addAction(style()->standardIcon(QStyle::SP_DialogOpenButton), tr("打开(&O)…"),
                                    QKeySequence::Open);
    connect(open, &QAction::triggered, this, [this] {
        auto f = QFileDialog::getOpenFileName(this, tr("打开 HDF5 / NeXus 文件"),
                                              QSettings().value("lastDirectory").toString(),
                                              "HDF5 / NeXus (*.nxs *.h5 *.hdf5 *.hdf *.nx5);;All files (*)");
        if (!f.isEmpty())
            openFile(f);
    });
    auto close = fileMenu->addAction(tr("关闭文件"));
    connect(close, &QAction::triggered, this, [this] {
        ++m_generation;
        m_pending = false;
        m_debounce.stop();
        m_play->setChecked(false);
        m_file.clear();
        m_format = ScanFormat::Generic;
        setWindowTitle(tr("XRI — 同步辐射数据浏览器"));
        m_node = {};
        m_slice = {};
        m_tree->clear();
        m_plot->setData({});
        m_tableModel->setSlice({});
        m_text->clear();
        m_path->setText(tr("未选择数据集"));
        m_summary->clear();
        m_stats->clear();
        m_configuring = true;
        m_xAxis->clear();
        m_yAxis->clear();
        buildSliceControls();
        m_configuring = false;
    });
    fileMenu->addSeparator();
    fileMenu->addAction(tr("退出(&X)"), QKeySequence::Quit, this, &QWidget::close);
    auto viewMenu = menuBar()->addMenu(tr("视图(&V)"));
    auto reset = viewMenu->addAction(tr("复位图表"));
    connect(reset, &QAction::triggered, this, [this] { m_plot->resetView(); });
    auto help = menuBar()->addMenu(tr("帮助(&H)"));
    help->addAction(tr("关于 XRI"), this, [this] {
        QMessageBox::about(this, "XRI", tr("XRI 0.1 · 非弹性非共振数据处理工具 Hu Shiqi husq@ihep.ac.cn"));
    });
    auto explorer = new QDockWidget(tr("文件结构"), this);
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
    layout->setContentsMargins(10, 8, 10, 8);
    layout->setSpacing(8);
    m_path = new QLabel(tr("未选择数据集"));
    QFont pathFont = m_path->font();
    pathFont.setPointSize(11);
    pathFont.setBold(true);
    m_path->setFont(pathFont);
    m_path->setTextInteractionFlags(Qt::TextSelectableByMouse);
    layout->addWidget(m_path);
    m_summary = new QLabel(tr("支持标量、字符串、一维信号和多维探测器数据"));
    m_summary->setTextInteractionFlags(Qt::TextSelectableByMouse);
    layout->addWidget(m_summary);
    auto plotTools = new QHBoxLayout;
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
    layout->addLayout(plotTools);
    m_tabs = new QTabWidget;
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
    layout->addWidget(m_tabs, 1);
    auto frames = new QHBoxLayout;
    m_play = new QToolButton;
    m_play->setText(tr("播放"));
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
    layout->addLayout(frames);
    m_stats = new QLabel;
    m_stats->setTextInteractionFlags(Qt::TextSelectableByMouse);
    m_stats->setWordWrap(true);
    layout->addWidget(m_stats);
    setCentralWidget(center);
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
    connect(m_play, &QToolButton::toggled, this, [this](bool play) {
        m_play->setText(play ? tr("暂停") : tr("播放"));
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
    statusBar()->addPermanentWidget(new QLabel(tr("只读 · HDF5 / NeXus")));
    statusBar()->showMessage(tr("就绪"));
    resizeDocks({explorer}, {310}, Qt::Horizontal);
}

void MainWindow::populate(QTreeWidgetItem *parent, const QString &path, int depth) {
    if (depth > 32) {
        parent->setToolTip(0, tr("已达到 32 层浏览深度，避免循环链接。"));
        return;
    }
    try {
        auto nodes = HdfFile::children(m_file, path);
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
    try {
        HdfFile::inspect(path, "/");
    } catch (const std::exception &e) {
        QMessageBox::critical(this, tr("打开失败"), QString::fromUtf8(e.what()));
        return;
    }
    ++m_generation;
    m_pending = false;
    m_debounce.stop();
    m_play->setChecked(false);
    m_file = QFileInfo(path).absoluteFilePath();
    m_format = HdfFile::formatFromFileName(path);
    m_indices.clear();
    m_node = {};
    m_tree->clear();
    m_filter->clear();
    QFileInfo info(path);
    setWindowTitle(info.fileName() + " [" + formatLabel() + tr("] — XRI"));
    QSettings().setValue("lastDirectory", info.absolutePath());
    auto root = new QTreeWidgetItem(m_tree, {info.fileName(), formatLabel()});
    root->setToolTip(0, m_file + "\n" + formatLabel());
    root->setToolTip(1, formatLabel());
    root->setData(0, PathRole, "/");
    root->setData(0, GroupRole, true);
    root->setData(0, DepthRole, 0);
    populate(root, "/", 1);
    root->setExpanded(true);
    m_tree->setCurrentItem(root);
    selectPath("/entry/data/D_LAMBDA");
}
void MainWindow::selectPath(const QString &path) {
    if (m_tree->topLevelItemCount() == 0)
        return;
    auto current = m_tree->topLevelItem(0);
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
            return;
        current->setExpanded(true);
        current = next;
    }
    m_tree->setCurrentItem(current);
    m_tree->scrollToItem(current);
}
void MainWindow::selectItem(QTreeWidgetItem *item) {
    ++m_generation;
    m_pending = false;
    m_debounce.stop();
    m_play->setChecked(false);
    m_slice = {};
    m_plot->setData({});
    m_tableModel->setSlice({});
    m_stats->clear();
    m_text->clear();
    try {
        m_node = HdfFile::inspect(m_file, item->data(0, PathRole).toString());
    } catch (const std::exception &e) {
        statusBar()->showMessage(QString::fromUtf8(e.what()));
        return;
    }
    m_path->setText(m_node.path);
    m_summary->setText(formatLabel() + "    |    " + m_node.type + "    " + m_node.shape);
    QString text = m_node.path + "\nScan format: " + formatLabel() + "\nType: " + m_node.type +
                   "\nShape: " + m_node.shape + "\n";
    for (const auto &attr : m_node.attributes)
        text += "\n@" + attr.first + " = " + attr.second;
    m_text->setPlainText(text);
    m_tabs->setCurrentIndex(m_node.numeric ? 0 : 2);
    configureAxes();
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
    if (e->mimeData()->hasUrls())
        openFile(e->mimeData()->urls().first().toLocalFile());
}
void MainWindow::closeEvent(QCloseEvent *e) {
    m_playTimer.stop();
    QSettings settings;
    settings.setValue("geometry", saveGeometry());
    settings.setValue("windowState", saveState());
    QMainWindow::closeEvent(e);
}
