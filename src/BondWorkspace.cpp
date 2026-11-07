#include "BondWorkspace.h"
#include "BondFigures.h"
#include "ImagingIO.h"
#include "ScientificPlot.h"
#include "SpectrumPlot.h"
#include "WorkflowUi.h"
#include <QtConcurrent/QtConcurrentRun>
#include <QtWidgets>
#include <cmath>

namespace {
const QString scans = QStringLiteral("NeXus / HDF5 (*.nxs *.h5 *.hdf5);;All files (*)");
const QString volumes = QStringLiteral("XRI four-dimensional HDF5 (*.h5 *.hdf5)");

QLabel *hint(const QString &text = {}) {
    auto label = new QLabel(text);
    label->setWordWrap(true);
    label->setTextInteractionFlags(Qt::TextSelectableByMouse);
    return label;
}

QDoubleSpinBox *number(int decimals, double minimum, double maximum, double initial) {
    auto box = new QDoubleSpinBox;
    box->setDecimals(decimals);
    box->setRange(minimum, maximum);
    box->setValue(initial);
    box->setKeyboardTracking(false);
    box->setMinimumWidth(90);
    return box;
}

QSpinBox *index(const QString &name) {
    auto box = new QSpinBox;
    box->setObjectName(name);
    box->setRange(0, 0);
    box->setKeyboardTracking(false);
    return box;
}

QVector<QPair<double, double>> windows(const QString &input) {
    QVector<QPair<double, double>> out;
    for (const auto &part : input.trimmed().split(QRegularExpression("[,，]+"), Qt::SkipEmptyParts)) {
        const auto bounds = part.split(QRegularExpression("[:：]"));
        bool loOk, hiOk;
        double lo = bounds.value(0).toDouble(&loOk), hi = bounds.value(1).toDouble(&hiOk);
        if (bounds.size() != 2 || !loOk || !hiOk || !std::isfinite(lo) || !std::isfinite(hi) || !(hi > lo))
            ImageUtil::fail(QObject::tr("能区格式为 起点:终点，可用逗号分隔多段。"));
        out.append({lo, hi});
    }
    return out;
}

DataSlice image(int width, int height, const QString &title, const std::function<double(int, int)> &pixel) {
    DataSlice result;
    result.path = title;
    result.width = width;
    result.height = height;
    result.xAxis = 1;
    result.yAxis = 0;
    result.dims = {quint64(height), quint64(width)};
    result.strides = {1, 1};
    result.starts = {0, 0};
    result.minimum = INFINITY;
    result.maximum = -INFINITY;
    result.values.reserve(qsizetype(width) * height);
    for (int y = 0; y < height; ++y)
        for (int x = 0; x < width; ++x) {
            const double value = pixel(x, y);
            result.values.append(value);
            if (std::isfinite(value)) {
                result.minimum = std::min(result.minimum, value);
                result.maximum = std::max(result.maximum, value);
                ++result.finiteCount;
            }
        }
    if (!result.finiteCount)
        result.minimum = result.maximum = 0;
    return result;
}

QString compact(double value) {
    return QString::number(value, 'g', 9);
}
} // namespace

BondWorkspace::BondWorkspace(const QString &sample, QWidget *parent)
    : QWidget(parent), m_base(Imaging::defaults()) {
    setObjectName("BondWorkspace");
    auto root = new QVBoxLayout(this);
    m_tabs = new QTabWidget;
    m_tabs->setObjectName("BondSteps");
    m_viewBar = new QWidget;
    auto viewRow = new QHBoxLayout(m_viewBar);
    viewRow->setContentsMargins(0, 0, 0, 0);
    viewRow->addWidget(new QLabel(tr("显示数据")));
    m_dataSource = new QComboBox;
    m_dataSource->setObjectName("BondViewSource");
    m_dataSource->addItem(tr("原始体素"));
    viewRow->addWidget(m_dataSource);
    m_singleEnergy = new QWidget;
    auto energyPicker = new QHBoxLayout(m_singleEnergy);
    energyPicker->setContentsMargins(14, 0, 0, 0);
    energyPicker->addWidget(new QLabel(tr("能量通道")));
    m_channel = index("BondSpectrumChannel");
    energyPicker->addWidget(m_channel);
    m_selectedEnergy = new QLabel;
    energyPicker->addWidget(m_selectedEnergy);
    viewRow->addWidget(m_singleEnergy);
    viewRow->addStretch();
    root->addWidget(m_viewBar);
    root->addWidget(m_tabs, 1);
    auto filePage = new QWidget;
    auto fileLayout = new QVBoxLayout(filePage);
    auto fileActions = new QHBoxLayout;
    auto import = new QPushButton(tr("载入成像配置…"));
    auto reopen = new QPushButton(tr("打开四维结果…"));
    fileActions->addWidget(import);
    fileActions->addWidget(reopen);
    fileActions->addStretch();
    fileLayout->addLayout(fileActions);
    auto files = new QFormLayout;
    auto path = [&](const QString &label, const QString &name, const QString &filter, bool output = false) {
        auto edit = new QLineEdit;
        edit->setObjectName(name);
        auto browse = new QPushButton(tr("浏览…"));
        auto row = new QHBoxLayout;
        row->addWidget(edit, 1);
        row->addWidget(browse);
        files->addRow(label, row);
        connect(browse, &QPushButton::clicked, this, [=] {
            auto selected = output ? QFileDialog::getSaveFileName(this, label, edit->text(), filter, nullptr,
                                                                  QFileDialog::DontConfirmOverwrite)
                                   : QFileDialog::getOpenFileName(this, label, edit->text(), filter);
            if (!selected.isEmpty()) {
                if (output && QFileInfo(selected).suffix().isEmpty())
                    selected += ".h5";
                edit->setText(selected);
            }
        });
        return edit;
    };
    m_sample = path(tr("样品 Y/Z/能量扫描"), "BondSample", scans);
    m_foil = path(tr("薄膜 X/能量扫描"), "BondFoil", scans);
    m_roi = path(tr("分析晶体 ROI 表"), "BondRois", tr("ROI 表 (*.txt *.tsv *.csv *.json);;All files (*)"));
    m_sampleExposure = number(5, 0, 100000, 0);
    m_foilExposure = number(5, 0, 100000, 0);
    m_sampleExposure->setSpecialValueText(tr("待确认"));
    m_foilExposure->setSpecialValueText(tr("待确认"));
    files->addRow(tr("样品曝光 / s"), m_sampleExposure);
    files->addRow(tr("薄膜曝光 / s"), m_foilExposure);
    m_monitor = new QLineEdit("D_SiC_Current1");
    m_detector = new QLineEdit("D_LAMBDA");
    fileLayout->addLayout(files);
    auto acquisition = new QWidget;
    auto acquisitionFields = new QFormLayout(acquisition);
    acquisitionFields->addRow(tr("I₀ 通道"), m_monitor);
    acquisitionFields->addRow(tr("探测器"), m_detector);
    fileLayout->addWidget(WorkflowUi::disclosure(tr("采集通道"), acquisition));
    fileLayout->addWidget(
        hint(tr("成像配置同时带入扫描循环、X 网格、薄膜参数和 ROI 排除项。请核对曝光及所选晶体。")));
    auto next = new QPushButton(tr("校验输入并继续"));
    next->setObjectName("CheckBondInputs");
    fileLayout->addWidget(next, 0, Qt::AlignLeft);
    fileLayout->addStretch();
    m_tabs->addTab(filePage, tr("1 样品与标定"));

    auto runPage = new QWidget;
    auto runLayout = new QVBoxLayout(runPage);
    files = new QFormLayout;
    m_output = path(tr("保存四维结果"), "BondOutput", volumes, true);
    auto solver = new QWidget;
    auto solverForm = new QFormLayout(solver);
    m_fitGains = new QCheckBox(tr("拟合相对晶体增益"));
    m_fitGains->setChecked(true);
    m_anchor = new QLineEdit("36_VU-C1");
    m_anchor->setMaximumWidth(150);
    m_iterations = new QSpinBox;
    m_iterations->setRange(1, 1000000);
    m_iterations->setValue(20000);
    m_tolerance = number(12, 1e-12, 1., 1e-7);
    m_tolerance->setSingleStep(1e-7);
    solverForm->addRow(m_fitGains);
    solverForm->addRow(tr("增益锚点"), m_anchor);
    solverForm->addRow(tr("NNLS 最大迭代"), m_iterations);
    solverForm->addRow(tr("投影 KKT 容差"), m_tolerance);
    m_quality = hint(tr("请先在第 1 步校验输入。"));
    runLayout->addWidget(m_quality);
    auto details = new QTabWidget;
    m_crystals = new QTableWidget(0, 3);
    m_crystals->setHorizontalHeaderLabels({tr("纳入"), tr("晶体"), tr("ROI 像素")});
    m_crystals->horizontalHeader()->setSectionResizeMode(1, QHeaderView::Stretch);
    m_crystals->verticalHeader()->hide();
    m_crystals->setEditTriggers(QAbstractItemView::NoEditTriggers);
    m_diagnostics = new QPlainTextEdit;
    m_diagnostics->setObjectName("BondDiagnostics");
    m_diagnostics->setReadOnly(true);
    details->addTab(m_crystals, tr("晶体 ROI"));
    details->addTab(m_diagnostics, tr("诊断详情"));
    runLayout->addWidget(details, 1);
    runLayout->addWidget(WorkflowUi::disclosure(tr("求解设置"), solver));
    runLayout->addLayout(files);
    m_run = new QPushButton(tr("开始重建并保存"));
    m_run->setObjectName("RunBondVolume");
    m_run->setEnabled(false);
    runLayout->addWidget(m_run, 0, Qt::AlignRight);
    m_tabs->addTab(runPage, tr("2 四维重建"));

    auto voxelPage = new QWidget;
    auto voxelLayout = new QVBoxLayout(voxelPage);
    auto coordinates = new QHBoxLayout;
    coordinates->addWidget(new QLabel(tr("空间索引")));
    m_x = index("BondVoxelX");
    m_y = index("BondVoxelY");
    m_z = index("BondVoxelZ");
    m_xControl = new QWidget;
    auto xRow = new QHBoxLayout(m_xControl);
    xRow->setContentsMargins(0, 0, 0, 0);
    xRow->addWidget(new QLabel("X"));
    xRow->addWidget(m_x);
    coordinates->addWidget(m_xControl);
    coordinates->addWidget(new QLabel("Y"));
    coordinates->addWidget(m_y);
    coordinates->addWidget(new QLabel("Z"));
    coordinates->addWidget(m_z);
    coordinates->addStretch();
    m_spectralStyleControl = new QWidget;
    auto styleRow = new QHBoxLayout(m_spectralStyleControl);
    styleRow->setContentsMargins(0, 0, 0, 0);
    styleRow->addWidget(new QLabel(tr("图形")));
    m_spectralStyle = new QComboBox;
    m_spectralStyle->setObjectName("BondSpectralStyle");
    m_spectralStyle->addItems({tr("三维曲面"), tr("谱线瀑布")});
    styleRow->addWidget(m_spectralStyle);
    coordinates->addWidget(m_spectralStyleControl);
    voxelLayout->addLayout(coordinates);

    auto energyOptions = new QWidget;
    auto energyRow = new QHBoxLayout(energyOptions);
    m_loss = new QCheckBox(tr("显示能量转移"));
    m_analyzer = number(8, 0., 100., 0.);
    m_analyzer->setSpecialValueText(tr("待输入"));
    m_analyzer->setToolTip(tr("输入适用于所有已选晶体的分析器参考能量。"));
    energyRow->addWidget(m_loss);
    energyRow->addWidget(new QLabel(tr("分析器参考 / keV")));
    energyRow->addWidget(m_analyzer);
    energyRow->addStretch();
    voxelLayout->addWidget(WorkflowUi::disclosure(tr("能量坐标设置"), energyOptions));
    auto processing = new QWidget;
    auto fitRow = new QHBoxLayout(processing);
    m_background = new QLineEdit;
    m_background->setObjectName("BondBackground");
    m_background->setPlaceholderText(tr("起点:终点，可填多段"));
    m_normalization = new QLineEdit;
    m_normalization->setPlaceholderText(tr("起点:终点"));
    auto apply = new QPushButton(tr("应用谱处理"));
    apply->setObjectName("ApplyBondSpectralProcessing");
    fitRow->addWidget(new QLabel(tr("线性背底")));
    fitRow->addWidget(m_background, 2);
    fitRow->addWidget(new QLabel(tr("面积归一化")));
    fitRow->addWidget(m_normalization, 1);
    fitRow->addWidget(apply);
    auto processingOptions = WorkflowUi::disclosure(tr("谱处理 · 扣背底与归一化"), processing);
    processingOptions->findChild<QToolButton *>()->setObjectName("BondProcessingOptions");
    voxelLayout->addWidget(processingOptions);
    m_voxelInfo = hint(tr("载入四维结果后可点击 XY 切片查看单体素谱。"));
    voxelLayout->addWidget(m_voxelInfo);
    m_plots = new QTabWidget;
    m_plots->setObjectName("BondSpectralPlots");
    auto pair = new QSplitter;
    m_slice = new ScientificPlot;
    m_slice->setObjectName("BondVoxelSlice");
    m_slice->setAxisLabels("Y (mm)", "X (mm)");
    m_slice->setProfilePickingEnabled(true);
    m_slice->setPickAnywhere(true);
    m_spectrum = new SpectrumPlot;
    m_spectrum->setObjectName("BondVoxelSpectrum");
    m_spectrum->setFineGrid(true);
    pair->addWidget(m_slice);
    pair->addWidget(m_spectrum);
    pair->setSizes({480, 480});
    m_plots->addTab(pair, tr("XY 选点与体素谱"));
    m_spectralFigure = new BondFigures;
    m_spectralFigure->setObjectName("BondEnergyPosition");
    m_plots->addTab(m_spectralFigure, tr("沿 X 三维谱图"));
    voxelLayout->addWidget(m_plots, 1);
    auto toFigures = new QPushButton(tr("查看衬度结果"));
    voxelLayout->addWidget(toFigures, 0, Qt::AlignRight);
    m_tabs->addTab(voxelPage, tr("3 选点与谱处理"));

    auto figurePage = new QWidget;
    auto figureLayout = new QVBoxLayout(figurePage);
    auto mapRow = new QHBoxLayout;
    m_mode = new QComboBox;
    m_mode->setObjectName("BondContrastMode");
    m_mode->addItems({tr("单能量强度"), tr("双能量差值 S(Eb) − S(Ea)")});
    m_dualEnergy = new QWidget;
    auto dualEnergyRow = new QHBoxLayout(m_dualEnergy);
    dualEnergyRow->setContentsMargins(10, 0, 0, 0);
    m_ea = number(10, -1e5, 1e5, 0.);
    m_eb = number(10, -1e5, 1e5, 0.);
    m_ea->setObjectName("BondEa");
    m_eb->setObjectName("BondEb");
    mapRow->addWidget(m_mode);
    dualEnergyRow->addWidget(new QLabel("Ea"));
    dualEnergyRow->addWidget(m_ea);
    dualEnergyRow->addWidget(new QLabel("Eb"));
    dualEnergyRow->addWidget(m_eb);
    mapRow->addWidget(m_dualEnergy);
    mapRow->addStretch();
    figureLayout->addLayout(mapRow);
    m_mapInfo = hint();
    m_mapInfo->setObjectName("BondMapInfo");
    figureLayout->addWidget(m_mapInfo);
    m_figures = new QTabWidget;
    auto montagePage = new QWidget;
    auto montageLayout = new QVBoxLayout(montagePage);
    auto montageRow = new QHBoxLayout;
    m_firstZ = index("BondFirstZ");
    m_countZ = index("BondCountZ");
    m_strideZ = index("BondStrideZ");
    m_countZ->setRange(1, 9);
    m_countZ->setValue(6);
    m_strideZ->setRange(1, 1);
    for (const auto &control : {qMakePair(tr("起始 Z 索引"), m_firstZ), qMakePair(tr("切片数"), m_countZ),
                                qMakePair(tr("间隔"), m_strideZ)}) {
        montageRow->addWidget(new QLabel(control.first));
        montageRow->addWidget(control.second);
    }
    montageRow->addStretch();
    montageLayout->addLayout(montageRow);
    m_montage = new BondFigures;
    m_montage->setObjectName("BondMontage");
    montageLayout->addWidget(m_montage, 1);
    m_figures->addTab(montagePage, tr("XY 切片组"));
    auto surfacePage = new QWidget;
    auto surfaceLayout = new QVBoxLayout(surfacePage);
    auto surfaceRow = new QHBoxLayout;
    m_level = number(12, -1e12, 1e12, 0);
    m_level2 = number(12, -1e12, 1e12, 0);
    m_secondSurface = new QCheckBox(tr("第二等值面"));
    auto mesh = new QPushButton(tr("生成等值面"));
    mesh->setObjectName("BuildBondSurface");
    surfaceRow->addWidget(new QLabel(tr("阈值")));
    surfaceRow->addWidget(m_level);
    surfaceRow->addWidget(m_secondSurface);
    surfaceRow->addWidget(m_level2);
    surfaceRow->addWidget(mesh);
    surfaceRow->addStretch();
    surfaceLayout->addLayout(surfaceRow);
    m_surface = new BondFigures;
    m_surface->setObjectName("BondIsosurface");
    m_surface->setMesh({}, {}, {}, tr("生成等值面后，按住左键拖动旋转。"));
    surfaceLayout->addWidget(m_surface, 1);
    m_figures->addTab(surfacePage, tr("三维等值面"));
    figureLayout->addWidget(m_figures, 1);
    m_tabs->addTab(figurePage, tr("4 衬度结果"));
    m_tabs->setTabEnabled(2, false);
    m_tabs->setTabEnabled(3, false);

    auto footer = new QHBoxLayout;
    m_status = hint();
    m_status->setObjectName("BondStatus");
    m_progress = new QProgressBar;
    m_progress->setMaximumWidth(160);
    m_progress->hide();
    m_cancel = new QPushButton(tr("取消"));
    m_cancel->hide();
    footer->addWidget(m_status, 1);
    footer->addWidget(m_progress);
    footer->addWidget(m_cancel);
    root->addLayout(footer);

    connect(import, &QPushButton::clicked, this, [this] {
        auto path = QFileDialog::getOpenFileName(this, tr("选择成像配置"), {}, "JSON (*.json)");
        if (!path.isEmpty())
            importConfiguration(path);
    });
    connect(reopen, &QPushButton::clicked, this, [this] {
        auto path = QFileDialog::getOpenFileName(this, tr("选择 XRI 四维结果"), {}, volumes);
        if (!path.isEmpty())
            openVolume(path);
    });
    connect(next, &QPushButton::clicked, this, [this] { start(true); });
    connect(m_run, &QPushButton::clicked, this, [this] { start(false); });
    connect(toFigures, &QPushButton::clicked, this, [this] { selectStep(3); });
    connect(m_fitGains, &QCheckBox::toggled, m_anchor, &QLineEdit::setEnabled);
    connect(m_cancel, &QPushButton::clicked, this, [this] {
        m_stop.store(true);
        status(tr("正在取消…"));
    });
    for (auto edit : {m_sample, m_foil, m_roi, m_detector, m_monitor})
        connect(edit, &QLineEdit::textChanged, this, [this] { invalidateInputs(); });
    for (auto spin : {m_sampleExposure, m_foilExposure})
        connect(spin, &QDoubleSpinBox::valueChanged, this, [this] { invalidateInputs(); });
    connect(m_crystals, &QTableWidget::itemChanged, this, [this](QTableWidgetItem *item) {
        if (item->column() != 0)
            return;
        auto config = m_base.value("roi").toObject();
        auto excluded = config.value("excluded_ids").toArray();
        const auto id = item->data(Qt::UserRole).toString();
        for (int i = int(excluded.size()) - 1; i >= 0; --i)
            if (excluded[i].toString() == id)
                excluded.removeAt(i);
        if (item->checkState() == Qt::Unchecked)
            excluded.append(id);
        config["excluded_ids"] = excluded;
        m_base["roi"] = config;
    });
    connect(m_loss, &QCheckBox::toggled, this, [this] { changeEnergyUnits(); });
    connect(m_analyzer, &QDoubleSpinBox::valueChanged, this, [this] { changeEnergyUnits(); });
    connect(m_channel, &QSpinBox::valueChanged, this, [this] {
        if (m_refreshing)
            return;
        syncViewControls();
        refreshVoxel(false);
        if (m_mode->currentIndex() == 0)
            refreshMap();
    });
    connect(m_tabs, &QTabWidget::currentChanged, this, [this] { syncViewControls(); });
    connect(m_plots, &QTabWidget::currentChanged, this, [this] { syncViewControls(); });
    connect(m_spectralStyle, &QComboBox::currentIndexChanged, this, [this] { refreshVoxel(); });
    for (auto pos : {m_x, m_y, m_z})
        connect(pos, &QSpinBox::valueChanged, this, [this, pos] { refreshVoxel(pos != m_x, pos == m_z); });
    connect(m_slice, &ScientificPlot::profilePointPicked, this, [this](QPoint p) {
        const QSignalBlocker blockX(m_x), blockY(m_y);
        m_x->setValue(p.y());
        m_y->setValue(p.x());
        refreshVoxel(true, false);
    });
    connect(apply, &QPushButton::clicked, this, [this] { applySpectralProcessing(); });
    connect(m_dataSource, &QComboBox::currentIndexChanged, this, [this] {
        refreshVoxel();
        refreshMap();
    });
    connect(m_mode, &QComboBox::currentIndexChanged, this, [this] {
        syncViewControls();
        refreshMap();
    });
    for (auto e : {m_ea, m_eb})
        connect(e, &QDoubleSpinBox::valueChanged, this, [this] { refreshMap(); });
    for (auto z : {m_firstZ, m_countZ, m_strideZ})
        connect(z, &QSpinBox::valueChanged, this, [this] { refreshMontage(); });
    connect(mesh, &QPushButton::clicked, this, [this] { buildSurfaces(); });
    m_level2->hide();
    connect(m_secondSurface, &QCheckBox::toggled, m_level2, &QDoubleSpinBox::setVisible);
    connect(&m_watcher, &QFutureWatcher<Work>::finished, this, [this] {
        auto result = m_watcher.result();
        m_tabs->setEnabled(true);
        m_viewBar->setEnabled(true);
        m_cancel->hide();
        m_progress->hide();
        if (!result.error.isEmpty()) {
            status(result.error, true);
            return;
        }
        if (result.validation) {
            showValidation(*result.validation);
            m_run->setEnabled(true);
            selectStep(1);
        }
        if (result.volume)
            installVolume(std::move(*result.volume));
        if (result.processed) {
            m_processed = std::move(*result.processed);
            m_normalized = result.normalized;
            {
                const QSignalBlocker block(m_dataSource);
                if (m_dataSource->count() == 1)
                    m_dataSource->addItem(tr("谱处理后"));
                m_dataSource->setCurrentIndex(1);
            }
            refreshVoxel();
            refreshMap();
        }
        if (result.surfaces)
            m_surface->setMesh(m_volume.axes, (*result.surfaces)[0], (*result.surfaces)[1],
                               result.meshCaption, result.surfaceColors);
        status(result.message);
    });
    if (HdfFile::formatFromFileName(sample) == ScanFormat::GridScan)
        m_sample->setText(sample);
    syncViewControls();
}

BondWorkspace::~BondWorkspace() {
    m_stop.store(true);
    m_watcher.waitForFinished();
}

bool BondWorkspace::isBusy() const {
    return !m_tabs->isEnabled();
}

void BondWorkspace::selectStep(int step) {
    if (m_tabs->isTabEnabled(step))
        m_tabs->setCurrentIndex(step);
}

void BondWorkspace::syncViewControls() {
    const int step = m_tabs->currentIndex();
    const bool voxel = m_plots->currentIndex() == 0, dual = m_mode->currentIndex() == 1;
    m_viewBar->setVisible(hasResult() && step >= 2);
    m_singleEnergy->setVisible((step == 2 && voxel) || (step == 3 && !dual));
    m_dualEnergy->setVisible(dual);
    m_xControl->setVisible(voxel);
    m_voxelInfo->setVisible(voxel);
    m_spectralStyleControl->setVisible(!voxel);
    if (hasResult() && m_channel->value() < m_volume.axes[0].size()) {
        const double incident = m_volume.axes[0][m_channel->value()];
        m_selectedEnergy->setText(
            tr("%1 %2").arg(compact(m_transfer ? 1000 * (incident - m_referenceKev) : incident),
                            m_transfer ? QStringLiteral("eV") : QStringLiteral("keV")));
    }
}

void BondWorkspace::invalidateInputs() {
    m_run->setEnabled(false);
    m_crystals->setRowCount(0);
    m_quality->setText(tr("请先在第 1 步校验输入。"));
    m_diagnostics->clear();
}

void BondWorkspace::status(const QString &text, bool error) {
    m_status->setText(text);
    m_status->setStyleSheet(error ? "color:#a02020" : "");
}

void BondWorkspace::launch(const QString &message, std::function<Work()> fn) {
    m_stop.store(false);
    m_tabs->setEnabled(false);
    m_viewBar->setEnabled(false);
    m_cancel->show();
    m_progress->setRange(0, 0);
    m_progress->show();
    status(message);
    m_watcher.setFuture(QtConcurrent::run([fn = std::move(fn)] {
        try {
            return fn();
        } catch (const std::exception &error) {
            Work result;
            result.error = QString::fromUtf8(error.what());
            return result;
        }
    }));
}

bool BondWorkspace::importConfiguration(const QString &path) {
    if (isBusy())
        return false;
    try {
        applyConfiguration(Imaging::loadConfig(path));
        status(tr("已载入成像配置，请核对文件和曝光。"));
        return true;
    } catch (const std::exception &error) {
        status(QString::fromUtf8(error.what()), true);
        return false;
    }
}

void BondWorkspace::applyConfiguration(const QJsonObject &config) {
    m_base = config;
    const auto input = config.value("input").toObject(), acq = config.value("acquisition").toObject(),
               rec = config.value("reconstruction").toObject();
    m_sample->setText(input.value("sample_nxs").toString());
    m_foil->setText(input.value("calibration_nxs").toString());
    m_roi->setText(input.value("roi_table").toString());
    invalidateInputs();
    m_sampleExposure->setValue(acq.value("sample_exposure_s").toDouble());
    m_foilExposure->setValue(acq.value("calibration_exposure_s").toDouble());
    m_monitor->setText(acq.value("monitor").toString("D_SiC_Current1"));
    m_detector->setText(acq.value("detector").toString("D_LAMBDA"));
    m_fitGains->setChecked(rec.value("fit_gains").toBool(true));
    m_anchor->setText(rec.value("anchor_roi").toString());
    m_iterations->setValue(rec.value("nnls_max_iterations").toInt(20000));
    m_tolerance->setValue(rec.value("nnls_tolerance").toDouble(1e-7));
}

QJsonObject BondWorkspace::configuration() const {
    auto c = m_base;
    c["input"] = QJsonObject{{"sample_nxs", m_sample->text().trimmed()},
                             {"calibration_nxs", m_foil->text().trimmed()},
                             {"roi_table", m_roi->text().trimmed()}};
    auto acq = c.value("acquisition").toObject(), rec = c.value("reconstruction").toObject(),
         runtime = c.value("runtime").toObject();
    acq["sample_exposure_s"] = m_sampleExposure->value();
    acq["calibration_exposure_s"] = m_foilExposure->value();
    acq["monitor"] = m_monitor->text().trimmed();
    acq["detector"] = m_detector->text().trimmed();
    rec["fit_gains"] = m_fitGains->isChecked();
    rec["anchor_roi"] = m_anchor->text().trimmed();
    rec["nnls_max_iterations"] = m_iterations->value();
    rec["nnls_tolerance"] = m_tolerance->value();
    const auto output = m_output->text().trimmed();
    runtime["output_dir"] = output.isEmpty() ? QDir::currentPath() : QFileInfo(output).absolutePath();
    c["acquisition"] = acq;
    c["reconstruction"] = rec;
    c["runtime"] = runtime;
    return Imaging::config(c, QDir::currentPath());
}

void BondWorkspace::start(bool validateOnly) {
    if (isBusy())
        return;
    if (validateOnly)
        invalidateInputs();
    else if (!m_run->isEnabled()) {
        status(tr("请先在第 1 步校验输入。"), true);
        return;
    }
    try {
        auto config = configuration();
        auto path = m_output->text().trimmed();
        if (!validateOnly && path.isEmpty()) {
            path = QFileDialog::getSaveFileName(this, tr("保存四维结果"),
                                                QFileInfo(m_sample->text()).completeBaseName() + "_bond.h5",
                                                volumes, nullptr, QFileDialog::DontConfirmOverwrite);
            if (path.isEmpty())
                return;
        }
        if (!validateOnly) {
            if (QFileInfo(path).suffix().isEmpty())
                path += ".h5";
            if (!QFileInfo(path).absoluteDir().exists())
                ImageUtil::fail(tr("四维结果的保存目录不存在。"));
            if (QFileInfo::exists(path))
                ImageUtil::fail(tr("结果文件已存在，请选择新文件名。"));
            m_output->setText(path);
            config = configuration();
        }
        if (!validateOnly)
            selectStep(1);
        launch(tr("正在检查扫描…"), [this, validateOnly, config, path] {
            Work out;
            if (validateOnly) {
                out.validation = BondImaging::validate(config, m_stop);
                out.message = tr("输入已校验，可在晶体表中调整纳入项。");
                return out;
            }
            out.volume = BondImaging::run(config, m_stop, [this](const QString &text, int percent) {
                QMetaObject::invokeMethod(
                    this,
                    [this, text, percent] {
                        if (isBusy() && !m_stop.load()) {
                            m_progress->setRange(0, 100);
                            m_progress->setValue(percent);
                            status(text);
                        }
                    },
                    Qt::QueuedConnection);
            });
            BondImaging::save(path, *out.volume, m_stop);
            out.message = tr("四维结果已保存到 %1").arg(path);
            return out;
        });
    } catch (const std::exception &error) {
        status(QString::fromUtf8(error.what()), true);
    }
}

void BondWorkspace::openVolume(const QString &path) {
    if (isBusy())
        return;
    launch(tr("正在打开四维结果…"), [this, path] {
        Work out;
        out.volume = BondImaging::load(path);
        ImageUtil::cancel(m_stop);
        out.message = tr("已打开四维结果 %1").arg(path);
        return out;
    });
}

void BondWorkspace::showValidation(const ImagingValidation &validation) {
    const QSignalBlocker block(m_crystals);
    m_crystals->setRowCount(int(validation.rois.size()));
    for (int row = 0; row < validation.rois.size(); ++row) {
        const auto &roi = validation.rois[row];
        auto choice = new QTableWidgetItem(roi.id);
        choice->setData(Qt::UserRole, roi.id);
        choice->setCheckState(roi.excluded ? Qt::Unchecked : Qt::Checked);
        if (roi.rect.isEmpty())
            choice->setFlags(choice->flags() & ~Qt::ItemIsUserCheckable);
        m_crystals->setItem(row, 0, choice);
        m_crystals->setItem(row, 1, new QTableWidgetItem(roi.label));
        m_crystals->setItem(row, 2, new QTableWidgetItem(QString::number(roi.pixels.size())));
    }
    m_diagnostics->setPlainText(QString::fromUtf8(QJsonDocument(validation.report).toJson()));
    m_quality->setText(tr("样品帧 %1 · 独立能量点 %2")
                           .arg(validation.sample.monitor.size())
                           .arg(validation.sample.axes[2].size()));
}

void BondWorkspace::installVolume(BondVolume volume) {
    m_refreshing = true;
    m_volume = std::move(volume);
    m_processed.clear();
    m_normalized = false;
    if (m_dataSource->count() > 1)
        m_dataSource->removeItem(1);
    m_background->clear();
    m_normalization->clear();
    const auto e = energy();
    auto setIndex = [](QSpinBox *box, qsizetype count, int start) {
        box->setMaximum(int(count) - 1);
        box->setValue(start);
    };
    int channel = int(e.size() / 2);
    if (!m_volume.converged[channel]) {
        auto first = std::find(m_volume.converged.cbegin(), m_volume.converged.cend(), true);
        if (first != m_volume.converged.cend())
            channel = int(first - m_volume.converged.cbegin());
    }
    setIndex(m_channel, e.size(), channel);
    const auto n = m_volume.voxelCount(), ny = m_volume.axes[2].size(), nz = m_volume.axes[3].size();
    auto begin = m_volume.intensity.cbegin() + channel * n;
    qsizetype brightest = m_volume.voxel(int(m_volume.axes[1].size() / 2), int(ny / 2), int(nz / 2));
    double maximum = 0;
    for (qsizetype i = 0; i < n; ++i)
        if (std::isfinite(begin[i]) && begin[i] > maximum) {
            maximum = begin[i];
            brightest = i;
        }
    setIndex(m_x, m_volume.axes[1].size(), int(brightest / (ny * nz)));
    setIndex(m_y, ny, int((brightest / nz) % ny));
    setIndex(m_z, nz, int(brightest % nz));
    setIndex(m_firstZ, m_volume.axes[3].size(), 0);
    m_strideZ->setMaximum(int(m_volume.axes[3].size() - 1));
    m_strideZ->setValue(std::max(1, int(m_volume.axes[3].size() / 6)));
    for (auto spin : {m_ea, m_eb}) {
        spin->setRange(e.first(), e.last());
        spin->setSingleStep(e[1] - e[0]);
    }
    int a = std::max(0, channel - 1), b = std::max(channel, 1);
    if (!m_volume.converged[a] || !m_volume.converged[b]) {
        QVector<int> valid;
        for (int i = 0; i < e.size(); ++i)
            if (m_volume.converged[i])
                valid.append(i);
        if (valid.size() >= 2) {
            a = valid.first();
            b = valid.last();
        }
    }
    m_ea->setValue(e[a]);
    m_eb->setValue(e[b]);
    m_mode->setCurrentIndex(0);
    m_diagnostics->setPlainText(QString::fromUtf8(QJsonDocument(m_volume.diagnostics).toJson()));
    const auto converged = std::count(m_volume.converged.cbegin(), m_volume.converged.cend(), true);
    m_quality->setText(
        tr("收敛通道 %1 / %2 · 选用晶体 %3 · 增益拟合 %4")
            .arg(converged)
            .arg(e.size())
            .arg(m_volume.diagnostics.value("crystal_keys").toArray().size())
            .arg(m_volume.diagnostics.value("gain_fit_converged").toBool() ? tr("已收敛") : tr("未收敛")));
    m_tabs->setTabEnabled(2, true);
    m_tabs->setTabEnabled(3, true);
    selectStep(2);
    m_refreshing = false;
    syncViewControls();
    refreshVoxel();
    refreshMap();
}

QVector<double> BondWorkspace::energy() const {
    if (!hasResult())
        return {};
    auto out = m_volume.axes[0];
    if (m_transfer)
        for (auto &e : out)
            e = 1000 * (e - m_referenceKev);
    return out;
}

QString BondWorkspace::energyLabel() const {
    return m_transfer ? tr("能量转移 (eV)") : tr("入射能量 (keV)");
}

QString BondWorkspace::intensityLabel() const {
    if (showingProcessed() && m_normalized)
        return m_transfer ? tr("面积归一强度 (1/eV)") : tr("面积归一强度 (1/keV)");
    return showingProcessed() ? tr("扣背底后的相对强度") : tr("相对强度");
}

bool BondWorkspace::showingProcessed() const {
    return m_dataSource->currentIndex() == 1 && !m_processed.isEmpty();
}

const QVector<double> &BondWorkspace::values() const {
    return showingProcessed() ? m_processed : m_volume.intensity;
}

void BondWorkspace::changeEnergyUnits() {
    const bool transfer = m_loss->isChecked() && m_analyzer->value() > 0;
    if (m_loss->isChecked() && !transfer)
        status(tr("请先填写适用于所选晶体的分析器参考能量。"), true);
    if (transfer == m_transfer && (!transfer || m_referenceKev == m_analyzer->value()))
        return;
    const double ea = m_transfer ? m_referenceKev + m_ea->value() / 1000 : m_ea->value();
    const double eb = m_transfer ? m_referenceKev + m_eb->value() / 1000 : m_eb->value();
    m_referenceKev = m_analyzer->value();
    m_transfer = transfer;
    if (!hasResult())
        return;
    m_refreshing = true;
    const auto e = energy();
    for (auto spin : {m_ea, m_eb}) {
        spin->setRange(e.first(), e.last());
        spin->setSingleStep(e[1] - e[0]);
    }
    m_ea->setValue(transfer ? 1000 * (ea - m_referenceKev) : ea);
    m_eb->setValue(transfer ? 1000 * (eb - m_referenceKev) : eb);
    m_background->clear();
    m_normalization->clear();
    m_processed.clear();
    m_normalized = false;
    if (m_dataSource->count() > 1)
        m_dataSource->removeItem(1);
    m_refreshing = false;
    syncViewControls();
    refreshVoxel();
    refreshMap();
    if (!m_loss->isChecked() || transfer)
        status(tr("能量坐标已更新，拟合区间按当前单位重新输入。"));
}

void BondWorkspace::refreshVoxel(bool refreshSection, bool refreshSlice) {
    if (m_refreshing || !hasResult())
        return;
    const auto e = energy();
    const auto &x = m_volume.axes[1], &y = m_volume.axes[2], &z = m_volume.axes[3];
    const int ix = m_x->value(), iy = m_y->value(), iz = m_z->value(), channel = m_channel->value();
    const auto n = m_volume.voxelCount();
    const auto &samples = values();
    const bool provisional = !m_volume.converged[channel];
    if (refreshSlice) {
        auto slice =
            image(int(y.size()), int(x.size()),
                  tr("XY · Z %1 mm · %2%3")
                      .arg(compact(z[iz]), compact(e[channel]), provisional ? tr(" · 未收敛") : QString{}),
                  [&](int col, int row) { return samples[channel * n + m_volume.voxel(row, col, iz)]; });
        m_slice->setImageCoordinates(y.first(), y[1] - y[0], x.first(), x[1] - x[0]);
        m_slice->setData(slice);
    }
    m_slice->setProfilePosition({iy, ix});
    const auto voxel = m_volume.voxel(ix, iy, iz);
    QVector<double> current, raw;
    for (int i = 0; i < e.size(); ++i) {
        const auto offset = i * n + voxel;
        current.append(samples[offset]);
        raw.append(m_volume.intensity[offset]);
    }
    QVector<SpectrumCurve> curves;
    if (showingProcessed() && !m_normalized)
        curves.append({tr("重建原始谱"), e, raw, QColor("#929292"), true});
    if (!showingProcessed() && std::any_of(m_volume.converged.cbegin(), m_volume.converged.cend(),
                                           [](bool valid) { return !valid; }))
        m_voxelInfo->setText(tr("X %1 · Y %2 · Z %3 mm    原始谱含未收敛通道，详见重建诊断")
                                 .arg(compact(x[ix]), compact(y[iy]), compact(z[iz])));
    else
        m_voxelInfo->setText(tr("X %1 · Y %2 · Z %3 mm    点击切片可选取体素")
                                 .arg(compact(x[ix]), compact(y[iy]), compact(z[iz])));
    curves.append({showingProcessed() ? tr("处理后") : tr("原始谱"), e, current, QColor("#0072b2")});
    m_spectrum->setCurves(curves, tr("体素谱"), energyLabel(), intensityLabel());
    if (!refreshSection)
        return;
    QVector<double> energySamples, positions;
    QVector<QVector<double>> lines;
    const int energyCount = std::min(256, int(e.size())), positionCount = std::min(256, int(x.size()));
    QVector<qsizetype> indices;
    for (int i = 0; i < energyCount; ++i) {
        const auto index = i * (e.size() - 1) / (energyCount - 1);
        energySamples.append(e[index]);
        indices.append(index);
    }
    for (int j = 0; j < positionCount; ++j) {
        const int k = int(j * (x.size() - 1) / (positionCount - 1));
        positions.append(x[k]);
        const auto offset = m_volume.voxel(k, iy, iz);
        QVector<double> row;
        row.reserve(energyCount);
        for (auto i : indices)
            row.append(samples[i * n + offset]);
        lines.append(std::move(row));
    }
    auto title = tr("Y %1 mm · Z %2 mm").arg(compact(y[iy]), compact(z[iz]));
    const bool sampled = energyCount < e.size() || positionCount < x.size();
    const bool waterfall = m_spectralStyle->currentIndex() == 1;
    if (waterfall)
        title += tr(" · %1 条谱线").arg(std::min(24, positionCount));
    if (sampled)
        title += tr(" · 显示抽样");
    m_spectralFigure->setSpectralSection(
        waterfall ? BondFigures::SpectralStyle::Waterfall : BondFigures::SpectralStyle::Surface,
        energySamples, positions, lines, energyLabel(), intensityLabel(), title);
}

void BondWorkspace::refreshMap() {
    if (m_refreshing || !hasResult())
        return;
    const bool dual = m_mode->currentIndex() == 1;
    m_field.clear();
    m_surface->setMesh({}, {}, {}, tr("生成等值面后，按住左键拖动旋转。"));
    try {
        const auto e = energy();
        const auto &samples = values();
        const int channel = m_channel->value();
        if (dual) {
            if (m_ea->value() == m_eb->value())
                ImageUtil::fail(tr("请选择两个不同的能量。"));
            m_field = BondImaging::energySlice(m_volume, samples, e, m_eb->value(), true);
            auto a = BondImaging::energySlice(m_volume, samples, e, m_ea->value(), true);
            for (qsizetype i = 0; i < a.size(); ++i)
                m_field[i] -= a[i];
            m_fieldLabel = tr("双能量衬度 · S(%1) − S(%2) · %3")
                               .arg(compact(m_eb->value()), compact(m_ea->value()), intensityLabel());
        } else {
            m_field = BondImaging::energySlice(m_volume, samples, e, e[channel], showingProcessed());
            m_fieldLabel = tr("单能量 %1 · %2%3")
                               .arg(compact(e[channel]), intensityLabel(),
                                    m_volume.converged[channel] ? QString{} : tr(" · 通道未收敛，仅供诊断"));
        }
        m_lo = INFINITY;
        m_hi = -INFINITY;
        qsizetype valid = 0;
        for (double value : std::as_const(m_field))
            if (std::isfinite(value)) {
                m_lo = std::min(m_lo, value);
                m_hi = std::max(m_hi, value);
                ++valid;
            }
        if (!valid)
            ImageUtil::fail(tr("当前设置下没有有效体素。"));
        const double range = m_hi - m_lo;
        m_level->setValue(m_lo + .65 * range);
        m_level2->setValue(m_lo + .35 * range);
        if (dual || m_lo < 0.) {
            m_hi = std::max({std::abs(m_lo), std::abs(m_hi), 1e-12});
            m_lo = -m_hi;
        } else if (m_lo == m_hi) {
            m_hi = m_lo + std::max(1e-10, std::abs(m_lo) * .001);
        }
        m_mapInfo->setText(tr("%1    ·    %2    ·    有效体素 %3 / %4")
                               .arg(m_fieldLabel, energyLabel())
                               .arg(valid)
                               .arg(m_field.size()) +
                           (valid < m_field.size() ? tr("    灰色为无效") : QString{}));
    } catch (const std::exception &error) {
        m_field.clear();
        m_mapInfo->setText(QString::fromUtf8(error.what()));
    }
    refreshMontage();
}

void BondWorkspace::refreshMontage() {
    if (m_refreshing || !hasResult())
        return;
    QVector<int> slices;
    for (int z = m_firstZ->value(); z < m_volume.axes[3].size() && slices.size() < m_countZ->value();
         z += m_strideZ->value())
        slices.append(z);
    m_montage->setMontage(m_volume.axes, m_field, slices, m_lo, m_hi,
                          m_mode->currentIndex() == 1 || m_lo < 0., intensityLabel());
}

void BondWorkspace::applySpectralProcessing() {
    if (!hasResult() || isBusy())
        return;
    try {
        BondSpectralSettings settings;
        settings.background = windows(m_background->text());
        const auto norm = windows(m_normalization->text());
        if (norm.size() > 1)
            ImageUtil::fail(tr("面积归一化请使用一段连续能区。"));
        if (!norm.isEmpty())
            settings.normalization = norm.first();
        if (settings.background.isEmpty() && !settings.normalization)
            ImageUtil::fail(tr("请填写线性背底或面积归一化能区。"));
        const auto volume = m_volume;
        const auto e = energy();
        launch(tr("正在逐体素处理谱…"), [this, volume, e, settings] {
            Work out;
            out.processed = BondImaging::process(volume, e, settings, m_stop);
            out.normalized = settings.normalization.has_value();
            out.message = tr("逐体素处理完成，未收敛通道和无效归一化体素以灰色显示。");
            return out;
        });
    } catch (const std::exception &error) {
        status(QString::fromUtf8(error.what()), true);
    }
}

void BondWorkspace::buildSurfaces() {
    if (m_field.isEmpty() || isBusy()) {
        status(tr("请先选择有有效体素的强度或衬度场。"), true);
        return;
    }
    const auto volume = m_volume;
    const auto field = m_field;
    const double level = m_level->value(), second = m_level2->value();
    const bool includeSecond = m_secondSurface->isChecked();
    const bool diverging = m_mode->currentIndex() == 1 || m_lo < 0.;
    const auto color = [&](double value) {
        return BondFigures::fieldColor((value - m_lo) / (m_hi - m_lo), diverging);
    };
    const std::array<QColor, 2> colors{color(level), color(second)};
    launch(tr("正在提取三维等值面…"), [this, volume, field, level, second, includeSecond, colors] {
        Work out;
        out.surfaceColors = colors;
        out.surfaces.emplace();
        (*out.surfaces)[0] = BondImaging::surface(volume, field, level, m_stop);
        if (includeSecond)
            (*out.surfaces)[1] = BondImaging::surface(volume, field, second, m_stop);
        out.meshCaption =
            tr("阈值 %1%2    ·    左键拖动旋转")
                .arg(compact(level), includeSecond ? tr("    /    %1").arg(compact(second)) : QString{});
        const auto triangles = (*out.surfaces)[0].size() + (*out.surfaces)[1].size();
        out.message = triangles ? tr("等值面已生成，共 %1 个三角形。").arg(triangles)
                                : tr("当前阈值没有等值面，请调整阈值。");
        return out;
    });
}
