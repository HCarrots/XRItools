#include "BondWorkspace.h"
#include "BondFigures.h"
#include "ComptonProcessing.h"
#include "ImagingIO.h"
#include "ImagingInputUi.h"
#include "ScientificPlot.h"
#include "SpectrumPlot.h"
#include "WorkflowUi.h"
#include "WorkspacePaths.h"
#include <QJsonArray>
#include <QtConcurrent/QtConcurrentRun>
#include <QtWidgets>
#include <cmath>
#include <limits>
#include <utility>

namespace {
using WorkspacePaths::Kind;
using WorkspacePaths::Workflow;

const QString volumes = QStringLiteral("XRI HDF5 (*.h5 *.hdf5)");
enum MapMode { SingleEnergy, EnergyDifference, IntervalIntegral };

QLabel *hint(const QString &text = {}) {
    auto label = new QLabel(text);
    label->setWordWrap(true);
    label->setTextInteractionFlags(Qt::TextSelectableByMouse);
    auto policy = label->sizePolicy();
    policy.setVerticalPolicy(QSizePolicy::Maximum);
    label->setSizePolicy(policy);
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

ProcessedSpectrum wideReference(const QString &path) {
    auto wide = SpectrumProcessing::load(path, {});
    if (!wide.wideCrystals.isEmpty())
        wide.crystals = std::move(wide.wideCrystals);
    return wide;
}

QStringList highMomentum(const ProcessedSpectrum &wide, const QStringList &labels, double minimumQ) {
    QStringList high;
    for (const auto &label : labels) {
        const double q = wide.geometry.value(label).qMean;
        if (!wide.crystals.contains(label) || !std::isfinite(q) || q <= 0)
            ImageUtil::fail(QObject::tr("参考谱缺少晶体 %1 的匹配谱或有效 q 值。").arg(label));
        if (q >= minimumQ)
            high.append(label);
    }
    return high;
}

QPair<double, double> readComptonWindows(ComptonSettings &parameters, const QStringList &ranges) {
    parameters.normalization = windows(ranges[0]);
    parameters.coreTail = windows(ranges[1]);
    parameters.cleanReference = windows(ranges[2]);
    parameters.highResidual = windows(ranges[3]);
    const auto match = windows(ranges[4]);
    if (parameters.highResidual.isEmpty() || match.size() != 1)
        ImageUtil::fail(QObject::tr("请填写高 q 残余直线区间和一段连续的参考幅度匹配区间。"));
    return match.first();
}

BondSpectralSettings comptonSettings(ComptonTransfer ref, const ComptonSettings &parameters,
                                     QPair<double, double> match) {
    const auto incidentWindow = [analyzer = ref.analyzerKev](QPair<double, double> range) {
        return qMakePair(analyzer + range.first / 1000, analyzer + range.second / 1000);
    };
    BondSpectralSettings settings;
    for (auto window : parameters.highResidual)
        settings.background.append(incidentWindow(window));
    settings.compton = BondSpectralReference{std::move(ref.correction), std::move(ref.reference),
                                            std::move(ref.valence), incidentWindow(match)};
    return settings;
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

QStringList contributingCrystals(const BondVolume &volume) {
    const auto keys = volume.diagnostics.value("crystal_keys").toArray();
    const auto labels = volume.diagnostics.value("crystal_labels").toArray();
    QStringList result;
    for (qsizetype i = 0; i < keys.size(); ++i) {
        const auto key = keys[i].toString();
        const auto label = labels.size() == keys.size() ? labels[i].toString() : key.section('_', 1);
        if (label.isEmpty() || result.contains(label))
            ImageUtil::fail(QObject::tr("重建结果的晶体标签缺失或重复，无法匹配逐晶体参考谱。"));
        result.append(label);
    }
    if (result.isEmpty())
        ImageUtil::fail(QObject::tr("重建结果缺少实际参与求解的晶体信息。"));
    return result;
}
} // namespace

BondWorkspace::BondWorkspace(const QString &sample, QWidget *parent)
    : QWidget(parent), m_base(Imaging::defaults()) {
    setObjectName("BondWorkspace");
    auto root = new QVBoxLayout(this);
    WorkflowUi::compact(root);
    m_tabs = new QTabWidget;
    m_tabs->setObjectName("BondSteps");
    m_tabs->setUsesScrollButtons(true);
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
    auto import = WorkflowUi::button("folder-open", tr("载入已有方案…"));
    import->setObjectName("ImportBondConfig");
    auto saveConfig = WorkflowUi::button("save", tr("保存当前配置…"));
    saveConfig->setObjectName("SaveBondConfig");
    auto reopen = WorkflowUi::button("folder-open", tr("查看已有结果…"));
    m_configTools = WorkflowUi::cornerTools({import, saveConfig, reopen});
    m_tabs->setCornerWidget(WorkflowUi::cornerTools({m_configTools, m_viewBar}));
    root->addWidget(m_tabs, 1);
    auto filePage = new QWidget;
    auto fileLayout = new QVBoxLayout(filePage);
    WorkflowUi::compact(fileLayout);
    auto columns = new QHBoxLayout;
    auto column = [&](int stretch) {
        auto body = new QWidget;
        auto layout = new QVBoxLayout(body);
        auto scroll = new QScrollArea;
        scroll->setWidgetResizable(true);
        scroll->setFrameShape(QFrame::NoFrame);
        scroll->setWidget(body);
        columns->addWidget(scroll, stretch);
        return layout;
    };
    auto left = column(3);
    auto right = column(2);
    auto addText = [this](QFormLayout *form, const QString &key, const QString &label, int kind = 0) {
        auto edit = new QLineEdit;
        edit->setObjectName(key);
        edit->setProperty("fieldKind", kind);
        m_inputs[key] = edit;
        form->addRow(label, edit);
        connect(edit, &QLineEdit::textChanged, this, &BondWorkspace::invalidateInputs);
    };
    auto addNumber = [this, &addText](QFormLayout *form, const QString &key, const QString &label,
                                      bool integer = false) {
        if (!integer) {
            addText(form, key, label, 1);
            return;
        }
        auto spin = new QSpinBox;
        spin->setObjectName(key);
        spin->setRange(1, std::numeric_limits<int>::max() - 1);
        spin->setKeyboardTracking(false);
        m_inputs[key] = spin;
        form->addRow(label, spin);
        connect(spin, &QSpinBox::valueChanged, this, &BondWorkspace::invalidateInputs);
    };
    ImagingInputUi::build(
        this, left, right, m_inputs, m_availableFiles, addText, addNumber, [this] { invalidateInputs(); },
        [this](const QString &file) {
            if (!QFileInfo(file).isFile())
                status(tr("请先选择样品文件。"), true);
            else
                emit roiRequested(file);
        });
    m_sample = static_cast<QLineEdit *>(m_inputs["input.sample_nxs"]);
    m_roi = static_cast<QLineEdit *>(m_inputs["input.roi_table"]);
    left->addStretch();

    right->addStretch();
    fileLayout->addLayout(columns, 1);
    auto next = WorkflowUi::button("check", tr("校验数据"));
    next->setObjectName("CheckBondInputs");
    fileLayout->addWidget(next, 0, Qt::AlignRight);
    m_tabs->addTab(filePage, tr("1 数据与 ROI"));

    m_tabs->addTab(createSpatialPage(), tr("2 ROI标定"));
    m_tabs->addTab(createCrystalPage(), tr("3 筛选晶体"));
    m_tabs->addTab(createVolumeBuildPage(), tr("4 四维体素重建"));

    auto voxelPage = new QWidget;
    auto voxelLayout = new QVBoxLayout(voxelPage);
    WorkflowUi::compact(voxelLayout);
    m_volumeQuality = hint({});
    voxelLayout->addWidget(m_volumeQuality);
    m_volumeFileStatus = hint();
    voxelLayout->addWidget(m_volumeFileStatus);
    m_volumeSaveTools = new QWidget;
    auto saveRow = new QHBoxLayout(m_volumeSaveTools);
    saveRow->setContentsMargins(0, 0, 0, 0);
    saveRow->addWidget(new QLabel(tr("新结果文件")));
    m_output = new QLineEdit;
    m_output->setObjectName("BondOutput");
    saveRow->addWidget(m_output, 1);
    m_chooseOutput = new QPushButton(tr("浏览…"));
    saveRow->addWidget(m_chooseOutput);
    m_workspaceOutput = hint();
    m_workspaceOutput->setObjectName("IxsWorkspaceOutput");
    saveRow->addWidget(m_workspaceOutput, 1);
    refreshOutputDestination();
    connect(m_chooseOutput, &QPushButton::clicked, this, [this] {
        auto path = QFileDialog::getSaveFileName(this, tr("选择 HDF5 保存位置"), m_output->text(), volumes, nullptr,
                                                QFileDialog::DontConfirmOverwrite);
        if (!path.isEmpty())
            m_output->setText(QFileInfo(path).suffix().isEmpty() ? path + ".h5" : path);
    });
    auto saveVolume = WorkflowUi::button("save", tr("保存原始四维体素 HDF5"));
    saveVolume->setObjectName("SaveBondVolume");
    saveRow->addWidget(saveVolume);
    connect(saveVolume, &QPushButton::clicked, this, &BondWorkspace::saveResult);
    voxelLayout->addWidget(m_volumeSaveTools);
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
    m_spectralStyleControl = new QWidget;
    auto styleRow = new QHBoxLayout(m_spectralStyleControl);
    styleRow->setContentsMargins(0, 0, 0, 0);
    styleRow->addWidget(new QLabel(tr("图形")));
    m_spectralStyle = new QComboBox;
    m_spectralStyle->setObjectName("BondSpectralStyle");
    m_spectralStyle->addItems({tr("三维曲面"), tr("谱线瀑布")});
    styleRow->addWidget(m_spectralStyle);

    auto energyOptions = new QWidget;
    auto energyRow = new QHBoxLayout(energyOptions);
    energyRow->setContentsMargins(0, 0, 0, 0);
    m_loss = new QCheckBox(tr("显示能量转移"));
    m_loss->setObjectName("BondShowEnergyTransfer");
    m_analyzer = number(8, 0., 100., 0.);
    m_analyzer->setObjectName("BondAnalyzerReference");
    m_analyzer->setSpecialValueText(tr("待输入"));
    m_analyzer->setToolTip(tr("输入适用于所有已选晶体的分析器参考能量。"));
    energyRow->addWidget(m_loss);
    energyRow->addWidget(new QLabel(tr("分析器参考 / keV")));
    energyRow->addWidget(m_analyzer);
    energyRow->addStretch();
    coordinates->addSpacing(6);
    coordinates->addWidget(WorkflowUi::disclosureButton(tr("能量坐标设置"), energyOptions));
    auto processing = new QWidget;
    auto processingLayout = new QVBoxLayout(processing);
    processingLayout->setContentsMargins(0, 0, 0, 0);
    auto fitRow = new QHBoxLayout;
    m_backgroundMode = new QComboBox;
    m_backgroundMode->setObjectName("BondBackgroundMode");
    m_backgroundMode->addItems({tr("逐体素线性扣除"), tr("按 q 分流（高 q 康普顿与低 q 线性）")});
    createComptonOptions();
    auto reference = WorkflowUi::button("sliders-horizontal", tr("宽谱与康普顿参数…"));
    reference->setObjectName("BondComptonParameters");
    reference->hide();
    connect(reference, &QPushButton::clicked, this, [this] { m_comptonOptions->exec(); });
    auto modeRow = new QHBoxLayout;
    modeRow->addWidget(m_backgroundMode);
    modeRow->addWidget(reference);
    modeRow->addStretch();
    processingLayout->addLayout(modeRow);
    m_background = new QLineEdit;
    m_background->setObjectName("BondBackground");
    m_background->setPlaceholderText(tr("起点:终点，可填多段"));
    m_normalization = new QLineEdit;
    m_normalization->setObjectName("BondAreaNormalization");
    m_normalization->setPlaceholderText(tr("起点:终点"));
    m_applyProcessing = WorkflowUi::button("filter", tr("应用谱处理"));
    m_applyProcessing->setObjectName("ApplyBondSpectralProcessing");
    auto linearLabel = new QLabel(tr("线性拟合 / keV"));
    fitRow->addWidget(linearLabel);
    fitRow->addWidget(m_background, 2);
    auto normalizationLabel = new QLabel(tr("面积归一化 / keV"));
    normalizationLabel->setToolTip(tr("可选。在扣背底之后对指定区间积分并归一化。"));
    fitRow->addWidget(normalizationLabel);
    fitRow->addWidget(m_normalization, 1);
    fitRow->addWidget(m_applyProcessing);
    processingLayout->addLayout(fitRow);
    auto processingHint = hint(tr("线性拟合可填多段入射能区。仅做面积归一化时，线性拟合区留空。"));
    processingLayout->addWidget(processingHint);
    connect(m_backgroundMode, &QComboBox::currentIndexChanged, this,
            [this, reference, linearLabel, processingHint](int mode) {
                reference->setVisible(mode != 0);
                linearLabel->setText(mode ? tr("低 q 拟合 / keV") : tr("线性拟合 / keV"));
                processingHint->setText(mode ? tr("纯低 q 体数据使用线性扣除。高 q 需单晶体重建和同一样品的逐晶体宽谱。")
                                             : tr("线性拟合可填多段入射能区。仅做面积归一化时，线性拟合区留空。"));
            });
    auto processingOptions = WorkflowUi::disclosureButton(tr("扣背底与归一化"), processing);
    processingOptions->setObjectName("BondProcessingOptions");
    coordinates->addWidget(processingOptions);
    coordinates->addStretch();
    coordinates->addWidget(m_spectralStyleControl);
    voxelLayout->addLayout(coordinates);
    voxelLayout->addWidget(energyOptions);
    voxelLayout->addWidget(processing);
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
    auto toFigures = WorkflowUi::button("layers", tr("查看三维结果分析"));
    toFigures->setObjectName("AnalyzeBondVolume");
    m_plots->setCornerWidget(toFigures);

    auto figurePage = new QWidget;
    auto figureLayout = new QVBoxLayout(figurePage);
    WorkflowUi::compact(figureLayout);
    m_mapControls = new QWidget;
    auto mapControlsLayout = new QVBoxLayout(m_mapControls);
    mapControlsLayout->setContentsMargins(0, 0, 0, 0);
    auto mapRow = new QHBoxLayout;
    m_mode = new QComboBox;
    m_mode->setObjectName("BondContrastMode");
    m_mode->addItems({tr("单能量强度"), tr("四维体素的双能量差值 S(Eb) − S(Ea)"), tr("区间积分强度")});
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
    m_integralBounds = new QWidget;
    auto boundsRow = new QHBoxLayout(m_integralBounds);
    boundsRow->setContentsMargins(10, 0, 0, 0);
    m_integralLo = number(10, -1e5, 1e5, 0.);
    m_integralHi = number(10, -1e5, 1e5, 0.);
    m_integralLo->setObjectName("BondIntegralLo");
    m_integralHi->setObjectName("BondIntegralHi");
    boundsRow->addWidget(new QLabel(tr("积分下限")));
    boundsRow->addWidget(m_integralLo);
    boundsRow->addWidget(new QLabel(tr("积分上限")));
    boundsRow->addWidget(m_integralHi);
    mapRow->addWidget(m_integralBounds);
    mapRow->addStretch();
    mapControlsLayout->addLayout(mapRow);
    m_mapInfo = hint();
    m_mapInfo->setObjectName("BondMapInfo");
    mapControlsLayout->addWidget(m_mapInfo);
    figureLayout->addWidget(m_mapControls);
    m_figures = new QTabWidget;
    auto montagePage = new QWidget;
    auto montageLayout = new QVBoxLayout(montagePage);
    WorkflowUi::compact(montageLayout);
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
    WorkflowUi::compact(surfaceLayout);
    auto surfaceRow = new QHBoxLayout;
    m_level = number(12, -1e12, 1e12, 0);
    m_level2 = number(12, -1e12, 1e12, 0);
    m_secondSurface = new QCheckBox(tr("第二等值面"));
    m_buildSurface = WorkflowUi::button("layers", tr("生成等值面"));
    m_buildSurface->setObjectName("BuildBondSurface");
    surfaceRow->addWidget(new QLabel(tr("阈值")));
    surfaceRow->addWidget(m_level);
    surfaceRow->addWidget(m_secondSurface);
    surfaceRow->addWidget(m_level2);
    surfaceRow->addWidget(m_buildSurface);
    surfaceRow->addStretch();
    surfaceLayout->addLayout(surfaceRow);
    m_surface = new BondFigures;
    m_surface->setObjectName("BondIsosurface");
    m_surface->setMesh({}, {}, {}, tr("生成等值面后，按住左键拖动旋转。"));
    surfaceLayout->addWidget(m_surface, 1);
    m_figures->addTab(surfacePage, tr("三维等值面"));
    m_referencePage = createReferenceMaps();
    m_figures->addTab(m_referencePage, tr("参考谱拟合四联图"));
    connect(m_figures, &QTabWidget::currentChanged, this, [this] { syncViewControls(); });
    figureLayout->addWidget(m_figures, 1);
    m_tabs->addTab(voxelPage, tr("5 四维结果与体素谱"));
    m_tabs->addTab(figurePage, tr("6 三维结果分析"));
    m_stepHint = hint();
    m_stepHint->setObjectName("BondStepHint");
    root->addWidget(m_stepHint);
    auto footer = new QHBoxLayout;
    m_status = hint();
    m_status->setObjectName("BondStatus");
    m_progress = new QProgressBar;
    m_progress->setMaximumWidth(160);
    m_progress->hide();
    m_cancel = WorkflowUi::button("square", tr("取消"));
    m_cancel->hide();
    footer->addWidget(m_status, 1);
    footer->addWidget(m_progress);
    footer->addWidget(m_cancel);
    root->addLayout(footer);

    connect(import, &QPushButton::clicked, this, [this] {
        auto path = QFileDialog::getOpenFileName(this, tr("载入已有方案"),
                                                 WorkspacePaths::directory(Workflow::Ixs, Kind::Config),
                                                 "JSON (*.json)");
        if (!path.isEmpty())
            importConfiguration(path);
    });
    connect(saveConfig, &QPushButton::clicked, this, [this] {
        const auto path = WorkspacePaths::saveFile(this, Workflow::Ixs, Kind::Config, m_sample->text(),
                                                   tr("保存当前配置"), "ixs.config.json", "JSON (*.json)",
                                                   "bond.config.json");
        if (!path.isEmpty())
            saveConfiguration(path);
    });
    connect(reopen, &QPushButton::clicked, this, [this] {
        auto path = QFileDialog::getOpenFileName(this, tr("选择 XRI 四维体素结果"),
                                                 WorkspacePaths::directory(Workflow::Ixs, Kind::Result), volumes);
        if (!path.isEmpty())
            openVolume(path);
    });
    connect(next, &QPushButton::clicked, this, &BondWorkspace::validateInputs);
    connect(m_run, &QPushButton::clicked, this, &BondWorkspace::startVolume);
    connect(toFigures, &QPushButton::clicked, this, [this] { selectStage(Results); });
    connect(m_fitGains, &QCheckBox::toggled, m_anchor, &QLineEdit::setEnabled);
    connect(m_cancel, &QPushButton::clicked, this, [this] {
        m_stop.store(true);
        status(tr("正在取消…"));
    });
    connect(m_crystals, &QTableWidget::itemChanged, this, [this](QTableWidgetItem *item) {
        if (item->column() != 0)
            return;
        auto field = m_inputs["roi.excluded_ids"];
        auto excluded = ImagingInputUi::read(field, "roi.excluded_ids").toArray();
        const auto id = item->data(Qt::UserRole).toString();
        for (int i = int(excluded.size()) - 1; i >= 0; --i)
            if (excluded[i].toString() == id)
                excluded.removeAt(i);
        if (item->checkState() == Qt::Unchecked)
            excluded.append(id);
        const QSignalBlocker block(field);
        ImagingInputUi::write(field, excluded);
        invalidateSpatial();
    });
    connect(m_loss, &QCheckBox::toggled, this, [this] { changeEnergyUnits(); });
    connect(m_analyzer, &QDoubleSpinBox::valueChanged, this, [this] { changeEnergyUnits(); });
    connect(m_channel, &QSpinBox::valueChanged, this, [this] {
        if (m_refreshing)
            return;
        syncViewControls();
        refreshVoxel(false);
        if (m_mode->currentIndex() == SingleEnergy)
            refreshMap();
    });
    connect(m_tabs, &QTabWidget::currentChanged, this, [this] {
        syncViewControls();
        refreshWorkflow();
    });
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
    connect(m_applyProcessing, &QPushButton::clicked, this, &BondWorkspace::applySpectralProcessing);
    connect(m_dataSource, &QComboBox::currentIndexChanged, this, [this] {
        refreshVoxel();
        refreshMap();
    });
    connect(m_mode, &QComboBox::currentIndexChanged, this, [this] {
        syncViewControls();
        refreshMap();
    });
    for (auto e : {m_ea, m_eb, m_integralLo, m_integralHi})
        connect(e, &QDoubleSpinBox::valueChanged, this, [this] { refreshMap(); });
    for (auto z : {m_firstZ, m_countZ, m_strideZ})
        connect(z, &QSpinBox::valueChanged, this, [this] { refreshMontage(); });
    connect(m_buildSurface, &QPushButton::clicked, this, &BondWorkspace::buildSurfaces);
    m_level2->hide();
    connect(m_secondSurface, &QCheckBox::toggled, m_level2, &QDoubleSpinBox::setVisible);
    connect(&m_watcher, &QFutureWatcher<Work>::finished, this, [this] {
        auto result = m_watcher.result();
        setBusy(false);
        m_cancel->hide();
        m_progress->hide();
        const auto pending = std::exchange(m_pendingRoi, {});
        if (!result.error.isEmpty()) {
            if (!pending.second.isEmpty())
                useRoiDraft(pending.first, pending.second);
            status(result.error, true);
            return;
        }
        if (result.inputs) {
            m_checkedInputs = std::move(result.inputs);
            showValidation(m_checkedInputs->validation);
            refreshWorkflow();
            selectStage(Spatial);
        }
        if (result.spatial) {
            m_spatial = std::move(result.spatial);
            showSpatialSummary();
            refreshWorkflow();
        }
        if (result.channels) {
            m_channels = std::move(result.channels);
            showChannelSummary();
            refreshWorkflow();
        }
        if (result.volume) {
            m_volumePendingSave = result.resultPath.isEmpty();
            m_volumeFileStatus->setText(m_volumePendingSave ? tr("四维体素已求解，尚未保存。")
                                                            : tr("结果文件：%1").arg(result.resultPath));
            installVolume(std::move(*result.volume));
        }
        if (result.savedVolume) {
            m_volumePendingSave = false;
            m_volume.bandDiagnostics.clear();
            m_volume.diagnosticImage = {};
            m_volumeFileStatus->setText(tr("结果文件：%1").arg(result.resultPath));
            syncViewControls();
        }
        if (result.processed) {
            m_processed = std::move(*result.processed);
            m_normalized = result.normalized;
            m_subtracted = result.subtracted;
            {
                const QSignalBlocker block(m_dataSource);
                if (m_dataSource->count() == 1)
                    m_dataSource->addItem(result.processingMethod);
                else
                    m_dataSource->setItemText(1, result.processingMethod);
                m_dataSource->setCurrentIndex(1);
            }
            refreshVoxel();
            refreshMap();
        }
        if (result.referenceSlice)
            showReferenceMaps(*result.referenceSlice);
        if (result.surfaces)
            m_surface->setMesh(spatialAxes(), (*result.surfaces)[0], (*result.surfaces)[1],
                               result.meshCaption, result.surfaceColors);
        status(result.message);
        if (!pending.second.isEmpty())
            useRoiDraft(pending.first, pending.second);
    });
    auto initial = Imaging::defaults();
    auto input = initial.value("input").toObject();
    if (!sample.isEmpty() && HdfFile::formatFromFileName(sample) != ScanFormat::Scan)
        input["sample_nxs"] = QFileInfo(sample).absoluteFilePath();
    initial["input"] = input;
    auto defaults = initial.value("acquisition").toObject();
    defaults["calibration_exposure_s"] = .5;
    defaults["sample_exposure_s"] = .5;
    initial["acquisition"] = defaults;
    applyConfiguration(initial);
    syncViewControls();
}

QWidget *BondWorkspace::createCrystalPage() {
    auto crystalsPage = new QWidget;
    auto crystalsLayout = new QVBoxLayout(crystalsPage);
    WorkflowUi::compact(crystalsLayout);
    m_quality = hint(tr("校验输入后核对参与计算的晶体。"));
    crystalsLayout->addWidget(m_quality);
    auto inputDetails = new QTabWidget;
    m_crystals = new QTableWidget(0, 3);
    m_crystals->setObjectName("BondCrystalSelection");
    m_crystals->setHorizontalHeaderLabels({tr("纳入"), tr("晶体"), tr("ROI 像素")});
    m_crystals->horizontalHeader()->setSectionResizeMode(1, QHeaderView::Stretch);
    m_crystals->verticalHeader()->hide();
    m_crystals->setEditTriggers(QAbstractItemView::NoEditTriggers);
    m_diagnostics = new QPlainTextEdit;
    m_diagnostics->setObjectName("BondDiagnostics");
    m_diagnostics->setReadOnly(true);
    inputDetails->addTab(m_crystals, tr("晶体 ROI"));
    inputDetails->addTab(m_diagnostics, tr("输入校验记录"));
    crystalsLayout->addWidget(inputDetails, 1);
    auto actions = new QHBoxLayout;
    actions->addStretch();
    m_backToSpatial = WorkflowUi::button("chevron-left", tr("返回第 2 步更新标定"));
    m_backToSpatial->setObjectName("UpdateBondSpatial");
    actions->addWidget(m_backToSpatial);
    connect(m_backToSpatial, &QPushButton::clicked, this, [this] { selectStage(Spatial); });
    m_confirmCrystals = WorkflowUi::button("chevron-right", tr("确认晶体并进入四维重建"));
    m_confirmCrystals->setObjectName("ConfirmBondCrystals");
    actions->addWidget(m_confirmCrystals);
    connect(m_confirmCrystals, &QPushButton::clicked, this, &BondWorkspace::confirmCrystals);
    crystalsLayout->addLayout(actions);
    return crystalsPage;
}

QWidget *BondWorkspace::createSpatialPage() {
    auto spatialPage = new QWidget;
    auto spatialLayout = new QVBoxLayout(spatialPage);
    WorkflowUi::compact(spatialLayout);
    spatialLayout->addWidget(hint(tr("用薄膜标定晶体的空间响应，用样品积分谱拟合相对增益。两条路线共用这套模型。")));
    auto solver = new QGroupBox(tr("ROI 标定参数"));
    auto solverForm = new QFormLayout(solver);
    m_xStep = new QLineEdit;
    m_xStep->setObjectName("BondXStep");
    m_xStep->setPlaceholderText(tr("自动（按薄膜标定）"));
    solverForm->addRow(tr("X 网格步长 / mm"), m_xStep);
    m_fitGains = new QCheckBox(tr("拟合相对晶体增益"));
    m_fitGains->setChecked(true);
    m_anchor = new QLineEdit("36_VU-C1");
    m_anchor->setObjectName("BondGainAnchor");
    m_iterations = new QSpinBox;
    m_iterations->setRange(1, 1000000);
    m_iterations->setValue(20000);
    m_tolerance = number(12, 1e-12, 1., 1e-7);
    m_tolerance->setSingleStep(1e-7);
    solverForm->addRow(m_fitGains);
    solverForm->addRow(tr("增益锚点（编号_标签）"), m_anchor);
    solverForm->addRow(tr("NNLS 最大迭代"), m_iterations);
    solverForm->addRow(tr("投影 KKT 容差"), m_tolerance);
    m_calibrate = WorkflowUi::button("play", tr("计算 ROI 标定"));
    m_calibrate->setObjectName("PrepareBondSpatial");
    solverForm->addRow(m_calibrate);
    auto spatialColumns = new QHBoxLayout;
    auto spatialControls = new QVBoxLayout;
    spatialControls->addWidget(solver);
    spatialControls->addStretch();
    spatialColumns->addLayout(spatialControls, 2);
    auto spatialReview = new QVBoxLayout;
    m_spatialStatus = hint(tr("计算后检查体素网格、采用的晶体和增益收敛状态。"));
    m_spatialStatus->setObjectName("BondSpatialSummary");
    spatialReview->addWidget(m_spatialStatus);
    auto spatialDetails = new QTabWidget;
    m_spatialTable = new QTableWidget(0, 3);
    m_spatialTable->setObjectName("BondSpatialCrystals");
    m_spatialTable->setHorizontalHeaderLabels({tr("晶体编号"), tr("相对增益"), tr("空间模型或拒绝原因")});
    m_spatialTable->horizontalHeader()->setStretchLastSection(true);
    m_spatialTable->setEditTriggers(QAbstractItemView::NoEditTriggers);
    spatialDetails->addTab(m_spatialTable, tr("晶体标定"));
    m_spatialDiagnostics = new QPlainTextEdit;
    m_spatialDiagnostics->setReadOnly(true);
    spatialDetails->addTab(m_spatialDiagnostics, tr("标定记录"));
    spatialReview->addWidget(spatialDetails, 1);
    spatialColumns->addLayout(spatialReview, 3);
    spatialLayout->addLayout(spatialColumns, 1);
    m_nextAfterSpatial = WorkflowUi::button("chevron-right", tr("下一步：筛选晶体"));
    m_nextAfterSpatial->setObjectName("ContinueBondSpatial");
    spatialLayout->addWidget(m_nextAfterSpatial, 0, Qt::AlignRight);
    connect(m_calibrate, &QPushButton::clicked, this, &BondWorkspace::prepareSpatial);
    connect(m_nextAfterSpatial, &QPushButton::clicked, this, [this] { selectStage(Crystals); });
    const auto invalidateCalibration = [this] { invalidateSpatial(); };
    for (auto *field : {m_xStep, m_anchor})
        connect(field, &QLineEdit::textChanged, this, invalidateCalibration);
    connect(m_fitGains, &QCheckBox::toggled, this, invalidateCalibration);
    connect(m_iterations, &QSpinBox::valueChanged, this, invalidateCalibration);
    connect(m_tolerance, &QDoubleSpinBox::valueChanged, this, invalidateCalibration);
    return spatialPage;
}

QWidget *BondWorkspace::createVolumeBuildPage() {
    auto channelsPage = new QWidget;
    auto channelsLayout = new QHBoxLayout(channelsPage);
    WorkflowUi::compact(channelsLayout);
    auto channelControls = new QVBoxLayout;
    auto backProject = new QGroupBox(tr("1. 逐能量反投影"));
    auto backLayout = new QVBoxLayout(backProject);
    backLayout->addWidget(hint(tr("汇总各能量、各 YZ 位置的剂量归一化计数，应用共用空间模型得到逐能量反投影。")));
    m_prepareChannels = WorkflowUi::button("play", tr("计算逐能量反投影"));
    m_prepareChannels->setObjectName("PrepareBondChannels");
    backLayout->addWidget(m_prepareChannels);
    channelControls->addWidget(backProject);
    auto volumeSolve = new QGroupBox(tr("2. 四维体素求解"));
    auto solveLayout = new QVBoxLayout(volumeSolve);
    solveLayout->addWidget(hint(tr("对每个入射能量独立求解非负体素。扣背底在求解后的体素谱上设置。")));
    auto flag = [this, solveLayout](const QString &key, const QString &label) {
        auto box = new QCheckBox(label);
        box->setObjectName(key);
        m_inputs[key] = box;
        solveLayout->addWidget(box);
    };
    flag("output.save_diagnostics", tr("保存积分谱的各晶体正演矩阵、预测和残差"));
    flag("output.plots", tr("保存积分谱诊断 PNG"));
    m_run = WorkflowUi::button("play", tr("求解四维体素"));
    m_run->setObjectName("RunBondVolume");
    solveLayout->addWidget(m_run);
    channelControls->addWidget(volumeSolve);
    channelControls->addStretch();
    channelsLayout->addLayout(channelControls, 2);
    auto channelReview = new QVBoxLayout;
    m_channelStatus = hint(tr("计算后检查各能量通道的归一化探测器信号。"));
    m_channelStatus->setObjectName("BondChannelSummary");
    channelReview->addWidget(m_channelStatus);
    m_channelTable = new QTableWidget(0, 2);
    m_channelTable->setObjectName("BondChannelInputs");
    m_channelTable->setHorizontalHeaderLabels({tr("入射能量 / keV"), tr("探测器信号 L2 范数")});
    m_channelTable->horizontalHeader()->setSectionResizeMode(QHeaderView::Stretch);
    m_channelTable->setEditTriggers(QAbstractItemView::NoEditTriggers);
    channelReview->addWidget(m_channelTable, 1);
    channelsLayout->addLayout(channelReview, 3);
    connect(m_prepareChannels, &QPushButton::clicked, this, &BondWorkspace::prepareChannels);
    return channelsPage;
}

BondWorkspace::~BondWorkspace() {
    m_stop.store(true);
    m_watcher.waitForFinished();
}

bool BondWorkspace::isBusy() const {
    return m_busy;
}

void BondWorkspace::selectStage(Stage stage) {
    m_tabs->setCurrentIndex(stage);
}

void BondWorkspace::confirmCrystals() {
    if (isBusy() || !m_spatial)
        return;
    m_crystalsConfirmed = true;
    refreshWorkflow();
    selectStage(VolumeBuild);
}

void BondWorkspace::refreshWorkflow() {
    const bool inputsReady = bool(m_checkedInputs);
    const bool spatialReady = bool(m_spatial);
    const bool reconstructionReady = spatialReady && m_crystalsConfirmed;
    m_calibrate->setEnabled(inputsReady);
    m_nextAfterSpatial->setEnabled(spatialReady);
    m_confirmCrystals->setEnabled(spatialReady);
    m_backToSpatial->setVisible(inputsReady && !spatialReady);
    m_prepareChannels->setEnabled(reconstructionReady);
    m_run->setEnabled(bool(m_channels));
    m_applyProcessing->setEnabled(hasVolume());
    m_buildSurface->setEnabled(hasVolume());
    m_fitReference->setEnabled(hasVolume());

    const int step = m_tabs->currentIndex();
    QString hint;
    if (step >= Volume && !hasVolume())
        hint = tr("请先完成第 4 步四维体素重建，或在第 1 步打开已有结果。");
    else if (step >= Spatial && step <= VolumeBuild && !inputsReady)
        hint = tr("请先在第 1 步校验数据。");
    else if (step >= Crystals && step <= VolumeBuild && !spatialReady)
        hint = tr("请先完成第 2 步 ROI 标定。");
    else if (step == VolumeBuild && !m_crystalsConfirmed)
        hint = tr("请先在第 3 步确认晶体。");
    m_stepHint->setText(hint);
    m_stepHint->setVisible(!hint.isEmpty());
}

void BondWorkspace::syncViewControls() {
    const int step = m_tabs->currentIndex();
    const bool voxel = m_plots->currentIndex() == 0;
    const int mode = m_mode->currentIndex();
    const bool referenceMaps = m_figures->currentWidget() == m_referencePage;
    m_mapControls->setVisible(!referenceMaps);
    const bool showViewTools = hasVolume() && (step == Volume || (step == Results && !referenceMaps));
    m_configTools->setVisible(step == Inputs);
    m_viewBar->setVisible(showViewTools);
    m_tabs->cornerWidget()->setVisible(step == Inputs || showViewTools);
    m_singleEnergy->setVisible((step == Volume && voxel) || (step == Results && mode == SingleEnergy));
    m_dualEnergy->setVisible(mode == EnergyDifference);
    m_integralBounds->setVisible(mode == IntervalIntegral);
    m_xControl->setVisible(voxel);
    m_voxelInfo->setVisible(voxel);
    m_spectralStyleControl->setVisible(!voxel);
    m_volumeSaveTools->setVisible(m_volumePendingSave);
    if (hasVolume() && m_channel->value() < m_volume.axes[0].size()) {
        const double incident = m_volume.axes[0][m_channel->value()];
        m_selectedEnergy->setText(
            tr("%1 %2").arg(compact(m_transfer ? 1000 * (incident - m_referenceKev) : incident),
                            m_transfer ? QStringLiteral("eV") : QStringLiteral("keV")));
    }
}

void BondWorkspace::invalidateInputs() {
    QFile roi(m_roi->text().trimmed());
    m_roiBytes = !roi.fileName().isEmpty() && roi.open(QIODevice::ReadOnly) ? roi.readAll() : QByteArray{};
    const int step = m_tabs->currentIndex();
    if (step >= Spatial && step <= VolumeBuild)
        selectStage(Inputs);
    m_checkedInputs.reset();
    invalidateSpatial();
    m_crystals->setRowCount(0);
    m_quality->setText(tr("校验输入后核对参与计算的晶体。"));
    m_diagnostics->clear();
}

void BondWorkspace::invalidateSpatial() {
    if (m_spatial)
        status(tr("晶体选择或空间参数已变更，请重新计算 ROI 标定。"));
    m_channels.reset();
    m_spatial.reset();
    m_crystalsConfirmed = false;
    showSpatialSummary();
    showChannelSummary();
    refreshWorkflow();
}

void BondWorkspace::useRoiDraft(const QString &file, const QString &path) {
    if (QFileInfo(m_sample->text().trimmed()).canonicalFilePath() != file)
        return;
    if (isBusy()) {
        m_pendingRoi = {file, path};
        return;
    }
    QFile roi(path);
    if (!roi.open(QIODevice::ReadOnly))
        return;
    if (m_roi->text() != path)
        m_roi->setText(path);
    else if (roi.readAll() != m_roiBytes)
        invalidateInputs();
}

void BondWorkspace::status(const QString &text, bool error) {
    m_status->setText(text);
    m_status->setStyleSheet(error ? "color:#a02020" : "");
}

void BondWorkspace::launch(const QString &message, std::function<Work()> fn) {
    m_stop.store(false);
    setBusy(true);
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

void BondWorkspace::setBusy(bool busy) {
    m_busy = busy;
    for (int i = 0; i < m_tabs->count(); ++i)
        m_tabs->widget(i)->setEnabled(!busy);
    m_configTools->setEnabled(!busy);
    m_viewBar->setEnabled(!busy);
}

bool BondWorkspace::importConfiguration(const QString &path) {
    if (isBusy())
        return false;
    try {
        QString output;
        applyConfiguration(Imaging::loadConfig(path, &output));
        if (!output.isEmpty())
            m_output->setText(output);
        status(tr("已载入方案，请核对文件和曝光。"));
        return true;
    } catch (const std::exception &error) {
        status(QString::fromUtf8(error.what()), true);
        return false;
    }
}

bool BondWorkspace::saveConfiguration(const QString &path) {
    try {
        auto config = configuration();
        if (WorkspacePaths::root().isEmpty() && !m_output->text().trimmed().isEmpty())
            config["bond_output_file"] = QFileInfo(m_output->text().trimmed()).absoluteFilePath();
        ImageUtil::saveJson(path, config);
        status(tr("配置已保存到 %1").arg(path));
        return true;
    } catch (const std::exception &error) {
        status(QString::fromUtf8(error.what()), true);
        return false;
    }
}

void BondWorkspace::refreshOutputDestination() {
    const auto folder = WorkspacePaths::directory(Workflow::Ixs, Kind::Result);
    m_output->setVisible(folder.isEmpty());
    m_chooseOutput->setVisible(folder.isEmpty());
    m_workspaceOutput->setVisible(!folder.isEmpty());
    m_workspaceOutput->setText(QDir::toNativeSeparators(folder));
}

void BondWorkspace::applyConfiguration(const QJsonObject &config) {
    if (isBusy())
        return;
    m_base = config;
    for (auto it = m_inputs.cbegin(); it != m_inputs.cend(); ++it) {
        const auto keys = it.key().split('.');
        const QSignalBlocker block(it.value());
        ImagingInputUi::write(it.value(), config.value(keys[0]).toObject().value(keys[1]));
    }
    const auto xStep = config.value("grid").toObject().value("x_step_mm");
    m_xStep->setText(xStep.isDouble() ? QString::number(xStep.toDouble(), 'g', 15) : QString{});
    const auto rec = config.value("reconstruction").toObject();
    m_fitGains->setChecked(rec.value("fit_gains").toBool(true));
    m_anchor->setText(rec.value("anchor_roi").toString());
    m_iterations->setValue(rec.value("nnls_max_iterations").toInt(20000));
    m_tolerance->setValue(rec.value("nnls_tolerance").toDouble(1e-7));
    const auto folder = config.value("runtime").toObject().value("output_dir").toString();
    const auto sample = QFileInfo(m_sample->text()).completeBaseName();
    m_output->setText(folder.isEmpty() ? QString{} : QDir(folder).filePath(sample + "_bond.h5"));
    invalidateInputs();
}

QJsonObject BondWorkspace::configuration(bool spatial) const {
    auto c = m_base;
    for (auto it = m_inputs.cbegin(); it != m_inputs.cend(); ++it) {
        const auto keys = it.key().split('.');
        auto group = c.value(keys[0]).toObject();
        group[keys[1]] = ImagingInputUi::read(it.value(), it.key());
        c[keys[0]] = group;
    }
    auto rec = c.value("reconstruction").toObject(), runtime = c.value("runtime").toObject();
    rec["fit_gains"] = spatial && m_fitGains->isChecked();
    if (spatial) {
        rec["anchor_roi"] = m_anchor->text().trimmed();
        rec["nnls_max_iterations"] = m_iterations->value();
        rec["nnls_tolerance"] = m_tolerance->value();
        auto grid = c.value("grid").toObject();
        const auto text = m_xStep->text().trimmed();
        if (text.isEmpty()) {
            grid["x_step_mm"] = QJsonValue();
        } else {
            bool ok;
            const double step = text.toDouble(&ok);
            if (!ok || !std::isfinite(step) || step <= 0)
                ImageUtil::fail(tr("X 网格步长必须为正数，留空时根据薄膜标定自动选择。"));
            grid["x_step_mm"] = step;
        }
        c["grid"] = grid;
    }
    const auto output = m_output->text().trimmed();
    const auto workspace = WorkspacePaths::directory(Workflow::Ixs, Kind::Result);
    runtime["output_dir"] = !workspace.isEmpty() ? workspace
                          : output.isEmpty()    ? QDir::currentPath()
                                                : QFileInfo(output).absolutePath();
    c["reconstruction"] = rec;
    c["runtime"] = runtime;
    return Imaging::config(c, QDir::currentPath());
}

void BondWorkspace::validateInputs() {
    if (isBusy())
        return;
    invalidateInputs();
    try {
        const auto config = configuration(false);
        launch(tr("正在校验扫描、剂量与 ROI…"), [this, config] {
            Work out;
            out.inputs = BondImaging::prepareInputs(config, m_stop);
            out.message = tr("输入已校验。请计算 ROI 标定。");
            return out;
        });
    } catch (const std::exception &error) {
        status(QString::fromUtf8(error.what()), true);
    }
}

void BondWorkspace::prepareSpatial() {
    if (isBusy() || !m_checkedInputs)
        return;
    try {
        const auto config = configuration();
        invalidateSpatial();
        launch(tr("正在计算 ROI 标定…"), [this, config, inputs = m_checkedInputs] {
            Work out;
            out.spatial = BondImaging::prepareSpatial(
                inputs, config, m_stop, [this](const QString &text, int percent) { reportProgress(text, percent); });
            out.message = tr("ROI 标定完成。检查体素网格和晶体增益后继续。");
            return out;
        });
    } catch (const std::exception &error) {
        status(QString::fromUtf8(error.what()), true);
    }
}

void BondWorkspace::prepareChannels() {
    if (isBusy() || !m_spatial || !m_crystalsConfirmed)
        return;
    m_channels.reset();
    showChannelSummary();
    refreshWorkflow();
    launch(tr("正在计算逐能量反投影…"), [this, spatial = m_spatial] {
        Work out;
        out.channels = BondImaging::prepareChannels(
            spatial, m_stop, [this](const QString &text, int percent) { reportProgress(text, percent); });
        out.message = tr("逐能量反投影已完成。检查通道后可求解四维体素。");
        return out;
    });
}

QString BondWorkspace::outputPath() {
    const auto input = m_volume.diagnostics.value("config").toObject().value("input").toObject();
    const auto automatic = WorkspacePaths::newPath(Workflow::Ixs, Kind::Result,
                                                   input.value("sample_nxs").toString(), "voxels.h5");
    if (!automatic.isEmpty())
        return automatic;
    auto path = m_output->text().trimmed();
    if (path.isEmpty() || QFileInfo::exists(path)) {
        const auto suggested = QFileInfo(m_sample->text()).completeBaseName() + "_bond.h5";
        path = QFileDialog::getSaveFileName(this, tr("保存四维结果"),
                                            QFileInfo(path).dir().filePath(suggested), volumes, nullptr,
                                            QFileDialog::DontConfirmOverwrite);
        if (path.isEmpty())
            return {};
    }
    if (QFileInfo(path).suffix().isEmpty())
        path += ".h5";
    if (QFileInfo::exists(path) || !QDir().mkpath(QFileInfo(path).absolutePath()))
        ImageUtil::fail(tr("无法写入结果，请选择尚未使用的文件名。"));
    m_output->setText(path);
    return path;
}

void BondWorkspace::startVolume() {
    if (isBusy() || !m_channels)
        return;
    try {
        const auto output = configuration(false).value("output").toObject();
        launch(tr("正在逐能量求解四维体素…"), [this, output, channels = m_channels] {
            Work out;
            out.volume = BondImaging::reconstructVolume(
                *channels, output, m_stop, [this](const QString &text, int percent) { reportProgress(text, percent); });
            out.message = tr("四维体素已求解。检查结果和通道收敛情况后保存。");
            return out;
        });
    } catch (const std::exception &error) {
        status(QString::fromUtf8(error.what()), true);
    }
}

void BondWorkspace::saveResult() {
    if (isBusy() || !m_volumePendingSave)
        return;
    try {
        const auto path = outputPath();
        if (path.isEmpty())
            return;
        launch(tr("正在写入 HDF5 结果…"), [this, path] {
            Work out;
            BondImaging::save(path, m_volume, m_stop);
            out.savedVolume = true;
            out.resultPath = path;
            out.message = tr("结果已保存到 %1").arg(path);
            return out;
        });
    } catch (const std::exception &error) {
        status(QString::fromUtf8(error.what()), true);
    }
}

void BondWorkspace::reportProgress(const QString &text, int percent) {
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
}

void BondWorkspace::openVolume(const QString &path) {
    if (isBusy())
        return;
    launch(tr("正在打开四维体素结果…"), [this, path] {
        Work out;
        out.volume = BondImaging::load(path);
        ImageUtil::cancel(m_stop);
        out.resultPath = path;
        out.message = tr("已打开结果 %1").arg(path);
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
        if (roi.sourceInvalid || roi.rect.isEmpty())
            choice->setFlags(choice->flags() & ~Qt::ItemIsUserCheckable);
        m_crystals->setItem(row, 0, choice);
        m_crystals->setItem(row, 1, new QTableWidgetItem(roi.label));
        m_crystals->setItem(row, 2, new QTableWidgetItem(QString::number(roi.pixels.size())));
    }
    m_diagnostics->setPlainText(QString::fromUtf8(QJsonDocument(validation.report).toJson()));
    const auto &incident = validation.sample.axes[2];
    m_quality->setText(tr("样品帧 %1 · 独立能量点 %2")
                           .arg(validation.sample.monitor.size())
                           .arg(incident.size()));
}

void BondWorkspace::showSpatialSummary() {
    if (!m_spatial) {
        m_spatialStatus->setText(tr("计算后检查体素网格、采用的晶体和增益收敛状态。"));
        m_spatialTable->setRowCount(0);
        m_spatialDiagnostics->clear();
        return;
    }
    const auto summary = BondImaging::spatialSummary(*m_spatial);
    const auto grid = summary.value("grid_xyz").toArray();
    const auto keys = summary.value("crystal_keys").toArray();
    const auto gains = summary.value("relative_gains").toArray();
    const auto models = summary.value("psf_models").toObject();
    const auto rejected = summary.value("psf_rejected").toObject();
    const auto gainStatus = !m_fitGains->isChecked() ? tr("未启用")
                          : summary.value("gain_fit_converged").toBool() ? tr("已收敛") : tr("未收敛");
    m_spatialStatus->setText(
        tr("网格 %1 × %2 × %3 · ΔX %4 mm · 采用 %5 块晶体 · 拒绝 %6 块 · 增益 %7 · 积分谱 NNLS %8")
            .arg(grid[0].toInt()).arg(grid[1].toInt()).arg(grid[2].toInt())
            .arg(compact(summary.value("x_step_mm").toDouble()))
            .arg(keys.size())
            .arg(rejected.size())
            .arg(gainStatus)
            .arg(summary.value("integrated_nnls_converged").toBool() ? tr("已收敛") : tr("未收敛")));
    m_spatialTable->setRowCount(int(keys.size() + rejected.size()));
    int row = 0;
    const auto addCrystal = [this, &row](const QString &key, const QString &gain, const QString &detail) {
        m_spatialTable->setItem(row, 0, new QTableWidgetItem(key));
        m_spatialTable->setItem(row, 1, new QTableWidgetItem(gain));
        auto item = new QTableWidgetItem(detail);
        item->setToolTip(detail);
        m_spatialTable->setItem(row++, 2, item);
    };
    const QMap<QString, QString> modelNames{{"constant_gaussian", tr("固定高斯")},
                                            {"varying_gaussian", tr("随 X 变化的高斯")},
                                            {"measured_template", tr("实测模板")}};
    for (qsizetype i = 0; i < keys.size(); ++i) {
        const auto key = keys[i].toString();
        const auto selection = models.value(key).toObject().value("model_selection").toObject();
        const auto model = selection.value("selected_model").toString();
        addCrystal(key, compact(gains[i].toDouble()), modelNames.value(model, model));
    }
    for (auto it = rejected.constBegin(); it != rejected.constEnd(); ++it)
        addCrystal(it.key(), QStringLiteral("—"), it.value().toString());
    m_spatialDiagnostics->setPlainText(QString::fromUtf8(QJsonDocument(summary).toJson()));
}

void BondWorkspace::showChannelSummary() {
    if (!m_channels) {
        m_channelStatus->setText(tr("计算后检查各能量通道的归一化探测器信号。"));
        m_channelTable->setRowCount(0);
        return;
    }
    const auto &energy = m_checkedInputs->validation.sample.axes[2];
    const auto norms = BondImaging::channelNorms(*m_channels);
    m_channelTable->setRowCount(int(energy.size()));
    for (int row = 0; row < energy.size(); ++row) {
        m_channelTable->setItem(row, 0, new QTableWidgetItem(compact(energy[row])));
        m_channelTable->setItem(row, 1, new QTableWidgetItem(compact(norms[row])));
    }
    m_channelStatus->setText(tr("已准备 %1 个能量通道，供下一步非负求解。").arg(energy.size()));
}

void BondWorkspace::installVolume(BondVolume volume) {
    m_refreshing = true;
    m_volume = std::move(volume);
    m_processed.clear();
    m_normalized = false;
    m_subtracted = false;
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
    m_strideZ->setMaximum(std::max(1, int(m_volume.axes[3].size() - 1)));
    m_strideZ->setValue(std::max(1, int(m_volume.axes[3].size() / 6)));
    for (auto spin : {m_ea, m_eb, m_integralLo, m_integralHi}) {
        spin->setRange(e.first(), e.last());
        spin->setSingleStep(e[1] - e[0]);
    }
    for (auto spin : {m_integralLo, m_integralHi})
        spin->setSuffix(m_transfer ? QStringLiteral(" eV") : QStringLiteral(" keV"));
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
    int start = 0, lo = 0, hi = 1;
    double widest = 0.;
    for (int i = 0; i < e.size(); ++i) {
        if (!m_volume.converged[i]) {
            start = i + 1;
            continue;
        }
        if (const double width = e[i] - e[start]; width > widest) {
            widest = width;
            lo = start;
            hi = i;
        }
    }
    m_integralLo->setValue(e[lo]);
    m_integralHi->setValue(e[hi]);
    m_mode->setCurrentIndex(SingleEnergy);
    const auto converged = std::count(m_volume.converged.cbegin(), m_volume.converged.cend(), true);
    const auto recon = m_volume.diagnostics.value("config").toObject().value("reconstruction").toObject();
    const auto gainStatus = !recon.value("fit_gains").toBool(true) ? tr("未启用")
                          : m_volume.diagnostics.value("gain_fit_converged").toBool() ? tr("已收敛") : tr("未收敛");
    m_volumeQuality->setText(
        tr("收敛通道 %1 / %2 · 选用晶体 %3 · 增益拟合 %4")
            .arg(converged)
            .arg(e.size())
            .arg(m_volume.diagnostics.value("crystal_keys").toArray().size())
            .arg(gainStatus));
    m_referenceZ->setRange(0, int(m_volume.axes[3].size()) - 1);
    m_referenceZ->setValue(m_z->value());
    for (auto *map : m_referenceMaps)
        map->setData({});
    m_referenceSpectra->setCurves({}, tr("D · 分区平均谱与参考谱"));
    refreshWorkflow();
    selectStage(Volume);
    m_refreshing = false;
    syncViewControls();
    refreshVoxel();
    refreshMap();
}

QVector<QVector<double>> BondWorkspace::spatialAxes() const {
    return hasVolume() ? m_volume.axes.sliced(1) : QVector<QVector<double>>{};
}

QVector<double> BondWorkspace::energy() const {
    if (!hasVolume())
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

double BondWorkspace::incidentEnergy(double displayedEnergy) const {
    return m_transfer ? m_referenceKev + displayedEnergy / 1000 : displayedEnergy;
}

QString BondWorkspace::intensityLabel() const {
    if (showingProcessed() && m_normalized)
        return tr("面积归一强度 (1/keV)");
    return showingProcessed() ? tr("扣背底后的相对强度") : tr("相对强度");
}

QString BondWorkspace::mapIntensityLabel() const {
    if (m_mode->currentIndex() != IntervalIntegral)
        return intensityLabel();
    return showingProcessed() && m_normalized ? tr("积分面积（无量纲）") : tr("积分强度（相对单位 × keV）");
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
    const std::array spins{m_ea, m_eb, m_integralLo, m_integralHi};
    std::array<double, spins.size()> incident{};
    for (size_t i = 0; i < spins.size(); ++i)
        incident[i] = incidentEnergy(spins[i]->value());
    m_referenceKev = m_analyzer->value();
    m_transfer = transfer;
    if (!hasVolume())
        return;
    m_refreshing = true;
    const auto e = energy();
    for (size_t i = 0; i < spins.size(); ++i) {
        auto spin = spins[i];
        spin->setRange(e.first(), e.last());
        spin->setSingleStep(e[1] - e[0]);
        spin->setValue(transfer ? 1000 * (incident[i] - m_referenceKev) : incident[i]);
    }
    for (auto spin : {m_integralLo, m_integralHi})
        spin->setSuffix(transfer ? QStringLiteral(" eV") : QStringLiteral(" keV"));
    m_refreshing = false;
    syncViewControls();
    refreshVoxel();
    refreshMap();
    if (!m_loss->isChecked() || transfer)
        status(tr("显示能量坐标已更新。谱处理使用设置中标明的单位。"));
}

void BondWorkspace::refreshVoxel(bool refreshSection, bool refreshSlice) {
    if (m_refreshing || !hasVolume())
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
    if (showingProcessed() && !m_normalized && m_subtracted) {
        QVector<double> removed(e.size());
        for (qsizetype i = 0; i < e.size(); ++i)
            removed[i] = raw[i] - current[i];
        curves.append({tr("扣除项"), e, removed, QColor("#b04a00"), true});
    }
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
    if (m_refreshing || !hasVolume())
        return;
    const int mode = m_mode->currentIndex();
    const bool dual = mode == EnergyDifference;
    m_field.clear();
    m_surface->setMesh({}, {}, {}, tr("生成等值面后，按住左键拖动旋转。"));
    try {
        const auto e = energy();
        const auto &samples = values();
        const int channel = m_channel->value();
        if (mode == IntervalIntegral) {
            const double lo = m_integralLo->value(), hi = m_integralHi->value();
            m_field = BondImaging::energyIntegral(m_volume, samples, {incidentEnergy(lo), incidentEnergy(hi)});
            m_fieldLabel = tr("区间积分强度 · %1 至 %2 · %3").arg(compact(lo), compact(hi), mapIntensityLabel());
        } else if (dual) {
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
    {
        const auto nz = int(m_volume.axes[3].size());
        const QSignalBlocker firstBlock(m_firstZ), strideBlock(m_strideZ);
        m_firstZ->setMaximum(nz - 1);
        m_strideZ->setMaximum(std::max(1, nz - 1));
    }
    refreshMontage();
}

void BondWorkspace::refreshMontage() {
    if (m_refreshing || !hasVolume())
        return;
    const auto axes = spatialAxes();
    QVector<int> slices;
    for (int z = m_firstZ->value(); axes.size() == 3 && z < axes[2].size() && slices.size() < m_countZ->value();
         z += m_strideZ->value())
        slices.append(z);
    m_montage->setMontage(axes, m_field, slices, m_lo, m_hi,
                          m_mode->currentIndex() == EnergyDifference || m_lo < 0., mapIntensityLabel());
}

QWidget *BondWorkspace::createReferenceMaps() {
    auto page = new QWidget;
    auto layout = new QVBoxLayout(page);
    WorkflowUi::compact(layout);
    auto files = new QGridLayout;
    for (int i = 0; i < 2; ++i) {
        const auto label = i == 0 ? tr("基体参考谱") : tr("目标参考谱");
        auto edit = new QLineEdit;
        edit->setObjectName(i == 0 ? "BondMatrixReference" : "BondTargetReference");
        edit->setPlaceholderText(tr("单条谱，能量转移以 eV 表示"));
        m_referenceFiles[i] = edit;
        auto browse = WorkflowUi::button("folder-open", tr("浏览…"));
        const int col = 3 * i;
        files->addWidget(new QLabel(label), 0, col);
        files->addWidget(edit, 0, col + 1);
        files->addWidget(browse, 0, col + 2);
        files->setColumnStretch(col + 1, 1);
        connect(browse, &QPushButton::clicked, page, [this, edit, label] {
            const auto path = QFileDialog::getOpenFileName(
                this, label, edit->text(), tr("参考谱 (*.h5 *.hdf5 *.tsv *.txt)"));
            if (!path.isEmpty())
                edit->setText(path);
        });
    }
    layout->addLayout(files);
    auto controls = new QHBoxLayout;
    m_referenceZ = index("BondReferenceZ");
    m_referenceAnalyzer = number(8, 0., 100., 0.);
    m_referenceAnalyzer->setObjectName("BondReferenceAnalyzer");
    m_referenceWindow = new QLineEdit;
    m_referenceWindow->setPlaceholderText(tr("自动取共同能区"));
    m_referenceWindow->setObjectName("BondReferenceWindow");
    m_referenceThreshold = number(8, 0., 1e9, 0.);
    m_referenceThreshold->setObjectName("BondReferenceThreshold");
    controls->addWidget(new QLabel(tr("Z 索引")));
    controls->addWidget(m_referenceZ);
    controls->addWidget(new QLabel(tr("分析器能量 / keV")));
    controls->addWidget(m_referenceAnalyzer);
    controls->addWidget(new QLabel(tr("拟合能区 / eV")));
    controls->addWidget(m_referenceWindow, 1);
    controls->addWidget(new QLabel(tr("系数阈值 ±")));
    controls->addWidget(m_referenceThreshold);
    m_fitReference = WorkflowUi::button("play", tr("拟合当前切片"));
    m_fitReference->setObjectName("FitBondReferenceMaps");
    controls->addWidget(m_fitReference);
    layout->addLayout(controls);
    layout->addWidget(hint(tr("使用四维原始体素谱和线性背底。基体与目标参考谱均须覆盖拟合能区，"
                              "且至少有五个已收敛通道。")));
    auto charts = new QGridLayout;
    const QStringList titles{tr("A · 总强度"), tr("B · 基体拟合残差平方和"), tr("C · 目标谱有符号系数")};
    for (int i = 0; i < 3; ++i) {
        auto plot = new ScientificPlot;
        plot->setMinimumSize(300, 225);
        plot->setAxisLabels(tr("X / mm"), tr("Y / mm"));
        plot->setObjectName(QStringLiteral("BondReferenceMap%1").arg(i));
        auto blank = DataSlice{};
        blank.path = titles[i];
        plot->setData(blank);
        m_referenceMaps[i] = plot;
        charts->addWidget(plot, i / 2, i % 2);
    }
    m_referenceMaps[1]->setPalette(1);
    m_referenceMaps[2]->setPalette(3);
    m_referenceSpectra = new SpectrumPlot;
    m_referenceSpectra->setObjectName("BondReferenceSpectra");
    m_referenceSpectra->setCurves({}, tr("D · 分区平均谱与参考谱"));
    charts->addWidget(m_referenceSpectra, 1, 1);
    for (int i = 0; i < 2; ++i) {
        charts->setRowStretch(i, 1);
        charts->setColumnStretch(i, 1);
    }
    layout->addLayout(charts, 1);
    connect(m_fitReference, &QPushButton::clicked, this, &BondWorkspace::fitReferenceMaps);
    const auto clear = [this, titles] {
        for (int i = 0; i < 3; ++i) {
            DataSlice blank;
            blank.path = titles[i];
            m_referenceMaps[i]->setData(blank);
        }
        m_referenceSpectra->setCurves({}, tr("D · 分区平均谱与参考谱"));
    };
    for (auto *edit : {m_referenceFiles[0], m_referenceFiles[1], m_referenceWindow})
        connect(edit, &QLineEdit::textChanged, page, clear);
    for (auto *spin : {m_referenceAnalyzer, m_referenceThreshold})
        connect(spin, &QDoubleSpinBox::valueChanged, page, clear);
    connect(m_referenceZ, &QSpinBox::valueChanged, page, clear);
    return page;
}

void BondWorkspace::fitReferenceMaps() {
    if (!hasVolume() || isBusy())
        return;
    try {
        const auto matrix = m_referenceFiles[0]->text().trimmed(), target = m_referenceFiles[1]->text().trimmed();
        if (matrix.isEmpty() || target.isEmpty())
            ImageUtil::fail(tr("请选择基体和目标的两条参考谱。"));
        const auto ranges = windows(m_referenceWindow->text());
        if (ranges.size() > 1)
            ImageUtil::fail(tr("参考谱拟合使用一段连续能区。"));
        const auto window = ranges.isEmpty() ? std::optional<QPair<double, double>>{} : std::make_optional(ranges.first());
        const auto analyzer = m_referenceAnalyzer->value(), threshold = m_referenceThreshold->value();
        if (analyzer <= 0.)
            ImageUtil::fail(tr("请填写与四维样品谱匹配的分析器能量。"));
        const int z = m_referenceZ->value();
        const auto volume = m_volume;
        launch(tr("正在拟合当前切片的体素谱…"), [this, volume, z, analyzer, threshold, matrix, target, window] {
            const auto read = [](const QString &path) {
                const auto data = SpectrumProcessing::load(path, {});
                if (data.crystals.size() != 1)
                    ImageUtil::fail(tr("每个参考文件必须恰好包含一条强度谱：%1").arg(path));
                return BondReferenceSpectrum{data.energyEv, data.crystals.first()};
            };
            Work out;
            const auto base = read(matrix);
            const auto component = read(target);
            out.referenceSlice = BondImaging::fitReferenceSlice(volume, z, analyzer, base, component, window,
                                                               threshold, m_stop);
            out.message = tr("切片拟合完成 · %1 个能量通道 · 正系数 %2 · 负系数 %3")
                              .arg(out.referenceSlice->energyEv.size())
                              .arg(out.referenceSlice->positiveCount)
                              .arg(out.referenceSlice->negativeCount);
            return out;
        });
    } catch (const std::exception &error) {
        status(QString::fromUtf8(error.what()), true);
    }
}

void BondWorkspace::showReferenceMaps(const BondReferenceSlice &result) {
    const auto &x = m_volume.axes[1], &y = m_volume.axes[2];
    const auto z = compact(m_volume.axes[3][result.z]);
    const QStringList titles{tr("A · 拟合通道总强度"), tr("B · 基体拟合残差平方和"),
                             tr("C · 目标谱系数（阈值 ±%1）").arg(compact(result.threshold))};
    const std::array fields{&result.total, &result.residual, &result.coefficient};
    for (int i = 0; i < 3; ++i) {
        auto slice = image(int(x.size()), int(y.size()), tr("%1 · Z %2 mm").arg(titles[i], z),
                          [&](int col, int row) {
            const double value = (*fields[i])[row * x.size() + col];
            return i == 2 && std::abs(value) <= result.threshold ? 0. : value;
        });
        if (i == 2 && slice.finiteCount) {
            const double span = std::max({std::abs(slice.minimum), std::abs(slice.maximum), result.threshold, 1e-15});
            slice.minimum = -span;
            slice.maximum = span;
        }
        m_referenceMaps[i]->setImageCoordinates(x.first(), x[1] - x[0], y.first(), y[1] - y[0]);
        m_referenceMaps[i]->setData(slice);
    }
    QVector<SpectrumCurve> curves;
    if (result.positiveCount)
        curves.append({tr("正系数区（%1）").arg(result.positiveCount), result.energyEv, result.positiveMean, QColor("#b22032")});
    if (result.negativeCount)
        curves.append({tr("负系数区（%1）").arg(result.negativeCount), result.energyEv, result.negativeMean, QColor("#285598")});
    curves.append({tr("基体参考"), result.energyEv, result.matrixReference, QColor("#a15b06"), true});
    curves.append({tr("目标参考"), result.energyEv, result.targetReference, QColor("#168264"), true});
    m_referenceSpectra->setCurves(curves, tr("D · 分区平均谱与参考谱"), tr("能量转移 / eV"), tr("归一化强度"));
}

void BondWorkspace::createComptonOptions() {
    m_comptonOptions = new QDialog(this);
    m_comptonOptions->setWindowTitle(tr("按 q 分流的参考与参数"));
    auto form = new QFormLayout(m_comptonOptions);
    form->addRow(hint(tr("请选择同一样品、相同晶体标定且已做 I0 归一化的逐晶体宽谱。晶体按平均 q 分流。"
                         "高 q 的能量转移按宽谱记录的晶体分析器能量计算。")));
    auto fileRow = [&](const QString &label, const QString &name, const QString &filter) {
        auto edit = new QLineEdit;
        edit->setObjectName(name);
        auto browse = new QPushButton(tr("浏览…"));
        auto row = new QHBoxLayout;
        row->addWidget(edit, 1);
        row->addWidget(browse);
        form->addRow(label, row);
        connect(browse, &QPushButton::clicked, m_comptonOptions, [this, edit, label, filter] {
            auto path = QFileDialog::getOpenFileName(m_comptonOptions, label, edit->text(), filter);
            if (!path.isEmpty())
                edit->setText(path);
        });
        return edit;
    };
    m_wideSpectrum = fileRow(tr("逐晶体宽谱"), "BondWideSpectrum", tr("逐晶体谱 (*.h5 *.hdf5 *.txt *.tsv)"));
    m_highQ = number(3, 0, 1000, 9.5);
    form->addRow(tr("高 q 最低值 / Å⁻¹"), m_highQ);
    m_hfTable = fileRow(tr("HF 壳层表"), "BondHFTable", tr("DABAX (*.dat);;所有文件 (*)"));
    const QSettings saved;
    m_hfTable->setText(saved.value("spectroscopy/hfTable").toString());
    m_composition = new QLineEdit(saved.value("spectroscopy/hfComposition").toString());
    m_composition->setObjectName("BondHFComposition");
    m_composition->setPlaceholderText("Li:1, Ni:0.9, Co:0.05, Mn:0.05, O:2");
    form->addRow(tr("元素及原子计量"), m_composition);
    m_coreScale = number(4, 0, 100, 1);
    form->addRow(tr("HF 芯谱比例"), m_coreScale);
    m_valenceCutoff = number(3, .01, 10000, 20);
    form->addRow(tr("价电子结合能上限 / eV"), m_valenceCutoff);
    const QStringList labels{tr("宽谱 HF 归一化"), tr("宽谱芯谱尾部直线"), tr("宽谱价电子轮廓拟合"),
                             tr("高 q 残余直线"), tr("参考幅度匹配")};
    const QStringList names{"HFNormalization", "CoreTail", "CleanReference", "HighResidual", "MatchWindow"};
    for (qsizetype i = 0; i < labels.size(); ++i) {
        m_comptonWindows[i] = new QLineEdit;
        m_comptonWindows[i]->setObjectName("BondCompton" + names[i]);
        m_comptonWindows[i]->setPlaceholderText(i == 0 || i == 4 ? tr("起点:终点")
                                                               : tr("起点:终点，可填多段"));
        form->addRow(labels[i] + tr(" / eV"), m_comptonWindows[i]);
    }
    form->addRow(hint(tr("幅度匹配使用一段与宽谱共有、避开待分析吸收边的能区。"
                         "低 q 只需逐晶体谱中的 q 和主页面的线性拟合区间。")));
    auto buttons = new QDialogButtonBox(QDialogButtonBox::Close);
    connect(buttons, &QDialogButtonBox::rejected, m_comptonOptions, &QDialog::reject);
    form->addRow(buttons);
}

void BondWorkspace::applySpectralProcessing() {
    if (!hasVolume() || isBusy())
        return;
    try {
        BondSpectralSettings settings;
        const auto linear = m_background->text();
        const auto norm = windows(m_normalization->text());
        if (norm.size() > 1)
            ImageUtil::fail(tr("面积归一化请使用一段连续能区。"));
        if (!norm.isEmpty())
            settings.normalization = norm.first();
        const bool split = m_backgroundMode->currentIndex() == 1;
        if (!split)
            settings.background = windows(linear);
        if (!split && settings.background.isEmpty() && !settings.normalization)
            ImageUtil::fail(tr("请填写线性背底或面积归一化能区。"));
        const auto path = m_wideSpectrum->text().trimmed();
        if (split && path.isEmpty())
            ImageUtil::fail(tr("请先选择含晶体 q 的逐晶体宽谱。"));
        ComptonSettings parameters;
        parameters.highQMin = m_highQ->value();
        parameters.coreScale = m_coreScale->value();
        parameters.valenceCutoffEv = m_valenceCutoff->value();
        parameters.hfTable = m_hfTable->text().trimmed();
        parameters.composition = m_composition->text().trimmed();
        QStringList ranges;
        for (auto edit : m_comptonWindows)
            ranges.append(edit->text());
        const auto volume = m_volume;
        auto process = [this, volume, settings, split, path, parameters, ranges, linear]() mutable {
            Work out;
            const auto &e = volume.axes[0];
            QString method = settings.background.isEmpty() ? tr("面积归一化") : tr("线性扣除");
            if (split) {
                const auto wide = wideReference(path);
                const auto labels = contributingCrystals(volume);
                const auto high = highMomentum(wide, labels, parameters.highQMin);
                if (high.isEmpty()) {
                    settings.background = windows(linear);
                    if (settings.background.isEmpty())
                        ImageUtil::fail(tr("请填写低 q 线性拟合区间。"));
                    method = tr("低 q 线性扣除");
                } else {
                    if (labels.size() != 1)
                        ImageUtil::fail(tr("此结果联合了多个晶体并含高 q 晶体。请在谱处理页仅纳入目标高 q 晶体，"
                                           "关闭相对晶体增益拟合并重新重建，再应用康普顿扣除。"));
                    const auto match = readComptonWindows(parameters, ranges);
                    const auto &label = high.first();
                    const auto normalization = settings.normalization;
                    settings = comptonSettings(transferCompton(wide, label, e, parameters, m_stop),
                                               parameters, match);
                    settings.normalization = normalization;
                    method = tr("高 q 康普顿扣除 · %1").arg(label);
                }
            }
            out.processed = BondImaging::process(volume, e, settings, m_stop);
            out.normalized = settings.normalization.has_value();
            out.subtracted = !settings.background.isEmpty();
            out.processingMethod = out.subtracted && out.normalized ? tr("%1 + 面积归一化").arg(method) : method;
            out.message = tr("%1已应用于四维结果。无有效处理值的位置以灰色显示。").arg(out.processingMethod);
            return out;
        };
        launch(split ? tr("正在校验晶体 q、准备参考并处理体素谱…") : tr("正在逐体素处理谱…"),
               std::move(process));
    } catch (const std::exception &error) {
        status(QString::fromUtf8(error.what()), true);
    }
}

void BondWorkspace::buildSurfaces() {
    if (m_field.isEmpty() || isBusy()) {
        status(tr("请先选择有有效体素的强度或衬度场。"), true);
        return;
    }
    const auto axes = spatialAxes();
    const auto field = m_field;
    const double level = m_level->value(), second = m_level2->value();
    const bool includeSecond = m_secondSurface->isChecked();
    const bool diverging = m_mode->currentIndex() == EnergyDifference || m_lo < 0.;
    const auto color = [&](double value) {
        return BondFigures::fieldColor((value - m_lo) / (m_hi - m_lo), diverging);
    };
    const std::array colors{color(level), color(second)};
    launch(tr("正在提取三维等值面…"), [this, axes, field, level, second, includeSecond, colors] {
        Work out;
        out.surfaceColors = colors;
        out.surfaces.emplace();
        (*out.surfaces)[0] = Isosurface::extract(axes, field.constData(), field.size(), level, m_stop);
        if (includeSecond)
            (*out.surfaces)[1] = Isosurface::extract(axes, field.constData(), field.size(), second, m_stop);
        out.meshCaption =
            tr("阈值 %1%2    ·    左键拖动旋转")
                .arg(compact(level), includeSecond ? tr("    /    %1").arg(compact(second)) : QString{});
        const auto triangles = (*out.surfaces)[0].size() + (*out.surfaces)[1].size();
        out.message = triangles ? tr("等值面已生成，共 %1 个三角形。").arg(triangles)
                                : tr("当前阈值没有等值面，请调整阈值。");
        return out;
    });
}
