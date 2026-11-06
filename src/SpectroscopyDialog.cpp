#include "SpectroscopyDialog.h"
#include "RoiDialog.h"
#include "ScientificPlot.h"
#include "SpectrumPlot.h"
#include "WorkflowUi.h"
#include <QtConcurrent/QtConcurrentRun>
#include <QtWidgets>

namespace {
QString number(double v) {
    return std::isfinite(v) ? QString::number(v, 'g', 9) : QStringLiteral("—");
}
QDoubleSpinBox *realBox(double value, double lo, double hi, int decimals = 4) {
    auto b = new QDoubleSpinBox;
    b->setDecimals(decimals);
    b->setRange(lo, hi);
    b->setValue(value);
    b->setKeyboardTracking(false);
    return b;
}
QComboBox *channelBox() {
    auto b = new QComboBox;
    b->setEditable(true);
    b->addItems({"auto", "/entry/instrument/M_DCM_B5Energy", "/entry/instrument/M_DCM_B5Energy_readback"});
    return b;
}
QTableWidgetItem *cell(const QString &text) {
    auto item = new QTableWidgetItem(text);
    item->setFlags(Qt::ItemIsSelectable | Qt::ItemIsEnabled);
    return item;
}
void configureTable(QTableWidget *t, const QStringList &headers) {
    t->setColumnCount(headers.size());
    t->setHorizontalHeaderLabels(headers);
    t->setAlternatingRowColors(true);
    t->setSelectionBehavior(QAbstractItemView::SelectRows);
    t->setSelectionMode(QAbstractItemView::ExtendedSelection);
    t->verticalHeader()->hide();
    t->horizontalHeader()->setStretchLastSection(true);
}
QWidget *scrollPanel(QWidget *body) {
    auto s = new QScrollArea;
    s->setWidgetResizable(true);
    s->setFrameShape(QFrame::NoFrame);
    s->setWidget(body);
    return s;
}
} // namespace

SpectroscopyDialog::SpectroscopyDialog(const QString &initialFile, QWidget *parent) : QDialog(parent) {
    setObjectName("SpectroscopyWorkbench");
    setWindowTitle(tr("谱学处理 — XRI"));
    setWindowFlag(Qt::Window);
    resize(1280, 860);
    setMinimumSize(980, 650);
    buildUi();
    connect(&m_watcher, &QFutureWatcher<Work>::finished, this, [this] {
        auto work = m_watcher.result();
        m_busy = false;
        m_tabs->setEnabled(true);
        m_progress->hide();
        if (!work.error.isEmpty()) {
            message(work.error, true);
            refreshState();
            return;
        }
        if (m_job == 0) {
            m_calibration = std::move(work.calibration);
            m_haveCalibration = true;
            refreshCalibration();
            m_tabs->setCurrentIndex(2);
            int accepted = 0;
            for (auto &fit : m_calibration.fits)
                if (fit.accepted)
                    ++accepted;
            message(tr("弹性标定完成：%1 / %2 个有效 ROI 通过拟合，重叠像素 %3。可在拟合表中查看未通过原因。")
                        .arg(accepted)
                        .arg(m_rois.count(RoiState::Valid))
                        .arg(m_calibration.overlapPixels));
        }
        if (m_job == 1) {
            m_scans = std::move(work.scans);
            m_haveScans = true;
            {
                QSignalBlocker block(m_scanChoice);
                m_scanChoice->clear();
                for (const auto &scan : m_scans)
                    m_scanChoice->addItem(QFileInfo(scan.file).fileName());
            }
            refreshScanCurves();
            refreshImages();
            m_tabs->setCurrentIndex(3);
            int count = 0;
            for (const auto &scan : m_scans)
                count += scan.glitches.size();
            message(tr("提取完成：%1 个扫描，修复 %2 个 I0 异常点及对应探测器帧。")
                        .arg(m_scans.size())
                        .arg(count));
        }
        if (m_job == 2) {
            m_result = std::move(work.result);
            m_haveResult = true;
            refreshResult();
            m_tabs->setCurrentIndex(4);
            int count = 0;
            for (bool v : m_result.included)
                count += v;
            message(tr("叠加完成：%1 个晶体，%2 个能量点，范围 %3 至 %4 eV。多扫描与晶体均求和。")
                        .arg(count)
                        .arg(m_result.total.size())
                        .arg(number(m_result.energyTransferEv.first()),
                             number(m_result.energyTransferEv.last())));
        }
        refreshState();
    });
    if (!initialFile.isEmpty() && HdfFile::formatFromFileName(initialFile) != ScanFormat::GridScan) {
        m_elasticFile->setText(QFileInfo(initialFile).absoluteFilePath());
        loadDraft();
    }
    refreshState();
}
SpectroscopyDialog::~SpectroscopyDialog() {
    m_cancel.store(true);
    m_watcher.waitForFinished();
}
void SpectroscopyDialog::selectStep(int step) {
    m_tabs->setCurrentIndex(std::clamp(step, 0, 4));
}
QString SpectroscopyDialog::statusText() const {
    return m_status->text();
}
void SpectroscopyDialog::message(const QString &text, bool error) {
    m_status->setText(text);
    m_status->setStyleSheet(error ? "color: #a02020;" : "");
}

void SpectroscopyDialog::buildUi() {
    auto root = new QVBoxLayout(this);
    m_tabs = new QTabWidget;
    m_tabs->setObjectName("SpectroscopySteps");
    root->addWidget(m_tabs, 1);
    auto inputs = new QWidget;
    auto inputLayout = new QHBoxLayout(inputs);
    auto forms = new QWidget;
    auto formLayout = new QVBoxLayout(forms);
    auto filesGroup = new QGroupBox(tr("弹性参考与共用 ROI"));
    auto fileForm = new QFormLayout(filesGroup);
    m_elasticFile = new QLineEdit;
    m_elasticFile->setObjectName("ElasticFile");
    auto browse = new QPushButton(tr("选择…"));
    auto fileRow = new QWidget;
    auto fl = new QHBoxLayout(fileRow);
    fl->setContentsMargins(0, 0, 0, 0);
    fl->addWidget(m_elasticFile, 1);
    fl->addWidget(browse);
    fileForm->addRow(tr("弹性扫描"), fileRow);
    connect(browse, &QPushButton::clicked, this, [this] {
        auto path =
            WorkflowUi::selectSource(this, m_availableFiles, m_elasticFile->text(), tr("选择弹性扫描"));
        if (!path.isEmpty()) {
            m_elasticFile->setText(path);
            loadDraft();
        }
    });
    m_roiPath = new QLineEdit;
    m_roiPath->setObjectName("SpectrumRoiScheme");
    m_roiPath->setReadOnly(true);
    m_roiPath->setPlaceholderText(tr("使用 ROI 菜单保存的方案或当前弹性文件的 ROI 草稿"));
    auto load = new QPushButton(tr("载入方案…"));
    auto edit = new QPushButton(tr("选取 ROI…"));
    auto roiButtons = new QWidget;
    auto rl = new QHBoxLayout(roiButtons);
    rl->setContentsMargins(0, 0, 0, 0);
    rl->addWidget(load);
    rl->addWidget(edit);
    fileForm->addRow(tr("ROI 方案"), m_roiPath);
    fileForm->addRow(roiButtons);
    m_roiStatus = new QLabel;
    m_roiStatus->setWordWrap(true);
    fileForm->addRow(m_roiStatus);
    connect(load, &QPushButton::clicked, this, [this] {
        auto path = QFileDialog::getOpenFileName(this, tr("载入共用 ROI 方案"), {}, "XRI ROI (*.json)");
        if (!path.isEmpty())
            loadRois(path);
    });
    connect(edit, &QPushButton::clicked, this, [this] {
        if (!QFileInfo::exists(m_elasticFile->text())) {
            message(tr("请先选择弹性扫描。"), true);
            return;
        }
        if (!isWindow()) {
            emit roiRequested(m_elasticFile->text());
            return;
        }
        RoiDialog dialog(m_elasticFile->text(), tr("谱学数据处理"), this);
        dialog.exec();
        loadDraft();
    });
    formLayout->addWidget(filesGroup);
    auto channels = new QGroupBox(tr("能量与入射光强"));
    auto cf = new QFormLayout(channels);
    m_elasticEnergy = channelBox();
    m_elasticEnergy->setObjectName("ElasticEnergy");
    m_xrsEnergy = channelBox();
    m_elasticI0 = new QComboBox;
    m_elasticI0->setEditable(true);
    m_elasticI0->addItem("/entry/data/D_SiC_Current1");
    m_xrsI0 = new QComboBox;
    m_xrsI0->setEditable(true);
    m_xrsI0->addItem("/entry/data/D_SiC_Current1");
    cf->addRow(tr("弹性能量 / keV"), m_elasticEnergy);
    cf->addRow(tr("XRS 能量 / keV"), m_xrsEnergy);
    cf->addRow(tr("弹性 I0"), m_elasticI0);
    cf->addRow(tr("XRS I0"), m_xrsI0);
    m_normalizeElastic = new QCheckBox(tr("弹性谱除以 |I0|"));
    m_normalizeElastic->setChecked(true);
    m_normalizeXrs = new QCheckBox(tr("XRS 谱除以 |I0|"));
    m_normalizeXrs->setChecked(true);
    m_repair = new QCheckBox(tr("修复孤立 I0 跌落与对应帧（邻点比例 < 0.7）"));
    m_repair->setChecked(true);
    cf->addRow(m_normalizeElastic);
    cf->addRow(m_normalizeXrs);
    cf->addRow(m_repair);
    auto energyNote = new QLabel(tr("auto 使用实测能量通道；请确认能量单位为 keV。能量轴须严格单调。"));
    energyNote->setWordWrap(true);
    cf->addRow(energyNote);
    formLayout->addWidget(WorkflowUi::disclosure(tr("高级：能量通道与 I0 归一化"), channels));
    auto mask = new QGroupBox(tr("ROI 与像素掩膜"));
    auto mf = new QFormLayout(mask);
    m_filter = new QCheckBox(tr("使用弹性叠加图生成的强度掩膜"));
    m_filter->setChecked(true);
    m_fraction = realBox(.15, 0, 1, 3);
    m_fraction->setSingleStep(.01);
    m_auto = new QCheckBox(tr("在原 ROI 内寻找峰值并重定位"));
    m_roiSize = new QSpinBox;
    m_roiSize->setRange(1, 2048);
    m_roiSize->setValue(20);
    m_overlap = new QCheckBox(tr("从所有涉及的 ROI 剔除重叠像素"));
    m_overlap->setChecked(true);
    mf->addRow(m_filter);
    mf->addRow(tr("阈值 / ROI 最大值"), m_fraction);
    mf->addRow(m_auto);
    mf->addRow(tr("自动 ROI 边长 / pixel"), m_roiSize);
    mf->addRow(m_overlap);

    auto fitting = new QGroupBox(tr("弹性峰拟合与质量阈值"));
    auto ff = new QFormLayout(fitting);
    m_model = new QComboBox;
    m_model->addItems({"Lorentzian", "Gaussian"});
    m_centerMin = realBox(9.67, 0, 100, 6);
    m_centerMax = realBox(9.69, 0, 100, 6);
    m_fwhm = realBox(2, .001, 10000, 3);
    m_r2 = realBox(.8, 0, 1, 4);
    ff->addRow(tr("峰形 + 常数背景"), m_model);
    ff->addRow(tr("中心下限 / keV"), m_centerMin);
    ff->addRow(tr("中心上限 / keV"), m_centerMax);
    ff->addRow(tr("最大 FWHM / eV"), m_fwhm);
    ff->addRow(tr("最小 R²"), m_r2);

    formLayout->addStretch();
    inputLayout->addWidget(scrollPanel(forms), 3);
    auto scanGroup = new QGroupBox(tr("XRS 扫描（勾选参与本次提取的文件）"));
    auto sl = new QVBoxLayout(scanGroup);
    m_files = new QListWidget;
    m_files->setObjectName("XrsFiles");
    m_files->setSelectionMode(QAbstractItemView::ExtendedSelection);
    sl->addWidget(m_files, 1);
    auto sb = new QHBoxLayout;
    auto add = new QPushButton(tr("添加扫描…"));
    auto remove = new QPushButton(tr("移除选中"));
    sb->addWidget(add);
    sb->addWidget(remove);
    sl->addLayout(sb);
    auto scanNote = new QLabel(tr("支持多个 scan 文件。取消勾选可剔除整次扫描。grid_scan "
                                  "使用成像流程。所有计算读取原文件，不修改源数据。"));
    scanNote->setWordWrap(true);
    sl->addWidget(scanNote);
    inputLayout->addWidget(scanGroup, 2);
    connect(add, &QPushButton::clicked, this, [this] {
        if (!m_availableFiles.isEmpty()) {
            QMenu menu(this);
            for (const auto &file : m_availableFiles) {
                if (HdfFile::formatFromFileName(file) == ScanFormat::GridScan)
                    continue;
                auto action = menu.addAction(QFileInfo(file).fileName());
                action->setData(file);
            }
            menu.addSeparator();
            auto all = menu.addAction(tr("添加全部已导入的 scan 文件"));
            auto disk = menu.addAction(tr("选择磁盘上的其他文件…"));
            auto selected = menu.exec(QCursor::pos());
            if (!selected)
                return;
            if (selected == all) {
                QStringList scans;
                for (const auto &file : m_availableFiles)
                    if (HdfFile::formatFromFileName(file) == ScanFormat::Scan)
                        scans.append(file);
                addScans(scans);
                return;
            }
            if (selected != disk) {
                addScans({selected->data().toString()});
                return;
            }
        }
        addScans(QFileDialog::getOpenFileNames(this, tr("添加 XRS 扫描"),
                                               QFileInfo(m_elasticFile->text()).absolutePath(),
                                               "NeXus (*.nxs *.h5 *.hdf5)"));
    });
    connect(remove, &QPushButton::clicked, this, [this] {
        for (auto item : m_files->selectedItems())
            delete m_files->takeItem(m_files->row(item));
        invalidate(1);
    });
    connect(m_files, &QListWidget::itemChanged, this, [this] { invalidate(1); });
    m_tabs->addTab(inputs, tr("1 数据与 ROI"));

    auto masks = new QWidget;
    auto ml = new QVBoxLayout(masks);
    ml->addWidget(WorkflowUi::disclosure(tr("像素掩膜参数"), mask));
    auto maskTools = new QHBoxLayout;
    maskTools->addWidget(new QLabel(tr("预览")));
    m_imageSource = new QComboBox;
    m_imageSource->addItems({tr("弹性叠加图"), tr("最终像素掩膜")});
    maskTools->addWidget(m_imageSource);
    auto log = new QCheckBox(tr("signed log1p 显示"));
    log->setChecked(true);
    maskTools->addWidget(log);
    auto saveRoi = new QPushButton(tr("保存调整后 ROI…"));
    maskTools->addStretch();
    maskTools->addWidget(saveRoi);
    ml->addLayout(maskTools);
    auto maskSplit = new QSplitter;
    auto correctionWidget = new QWidget;
    auto cl = new QVBoxLayout(correctionWidget);
    auto correctionNote =
        new QLabel(tr("自动重定位启用时应用偏移与扩展。失效晶体始终跳过；原方案保持不变。"));
    correctionNote->setWordWrap(true);
    cl->addWidget(correctionNote);
    m_corrections = new QTableWidget(75, 5);
    m_corrections->setObjectName("RoiCorrections");
    configureTable(m_corrections, {tr("晶体"), tr("ΔX"), tr("ΔY"), tr("扩 X"), tr("扩 Y")});
    for (int r = 0; r < 75; ++r) {
        m_corrections->setItem(r, 0, cell(m_rois.items()[r].label));
        for (int c = 1; c < 5; ++c) {
            auto item = new QTableWidgetItem("0");
            m_corrections->setItem(r, c, item);
        }
    }
    for (int c = 0; c < 5; ++c)
        m_corrections->setColumnWidth(c, c == 0 ? 80 : 52);
    cl->addWidget(m_corrections);
    maskSplit->addWidget(correctionWidget);
    auto detectorTabs = new QTabWidget;
    for (int d = 0; d < 2; ++d) {
        m_images[d] = new ScientificPlot;
        m_images[d]->setLogarithmic(true);
        detectorTabs->addTab(m_images[d], d == 0 ? "LAMBDA" : "MINIPIX");
    }
    maskSplit->addWidget(detectorTabs);
    maskSplit->setStretchFactor(1, 1);
    maskSplit->setSizes({340, 800});
    connect(m_corrections, &QTableWidget::currentCellChanged, this,
            [this, detectorTabs](int row, int, int, int) {
                if (!m_haveCalibration || row < 0)
                    return;
                detectorTabs->setCurrentIndex(row < 60 ? 0 : 1);
                for (int d = 0; d < 2; ++d) {
                    QVector<CrystalRoi> overlays;
                    for (const auto &roi : m_calibration.adjusted.items())
                        if (roi.detector == (d == 0 ? "D_LAMBDA" : "D_MINIPIX"))
                            overlays.append(roi);
                    m_images[d]->setRoiOverlays(overlays, m_rois.items()[row].label);
                }
            });
    ml->addWidget(maskSplit, 1);
    m_tabs->addTab(masks, tr("2 ROI 与掩膜"));
    connect(m_imageSource, &QComboBox::currentIndexChanged, this, [this] { refreshImages(); });
    connect(log, &QCheckBox::toggled, this, [this](bool on) {
        for (auto image : m_images)
            image->setLogarithmic(on);
    });
    connect(saveRoi, &QPushButton::clicked, this, [this] {
        if (!m_haveCalibration) {
            message(tr("请先完成弹性标定。"), true);
            return;
        }
        auto name =
            QFileDialog::getSaveFileName(this, tr("保存调整后 ROI"), "adjusted.roi.json", "XRI ROI (*.json)");
        if (name.isEmpty())
            return;
        QSaveFile file(name);
        auto bytes = m_calibration.adjusted.serialize(m_calibrationSettings.elasticFile, m_calibration.sizes);
        if (!file.open(QIODevice::WriteOnly) || file.write(bytes) != bytes.size() || !file.commit())
            message(tr("ROI 保存失败。"), true);
        else
            message(tr("已保存调整后 ROI：%1").arg(name));
    });

    auto elastic = new QWidget;
    auto el = new QVBoxLayout(elastic);
    el->addWidget(WorkflowUi::disclosure(tr("峰形拟合与质量阈值"), fitting));
    m_elasticView = new QComboBox;
    m_elasticView->addItems({tr("所选晶体：原 ROI / 掩膜 / 拟合"), tr("全部合格晶体：弹性零点对齐")});
    el->addWidget(m_elasticView);
    auto es = new QSplitter(Qt::Vertical);
    m_elasticPlot = new SpectrumPlot;
    m_elasticPlot->setObjectName("ElasticSpectrum");
    es->addWidget(m_elasticPlot);
    m_fits = new QTableWidget(75, 7);
    m_fits->setObjectName("ElasticFits");
    configureTable(m_fits, {tr("晶体"), tr("像素 / 掩膜"), tr("中心 / eV"), tr("FWHM / eV"), "R²", tr("拟合"),
                            tr("原因")});
    es->addWidget(m_fits);
    es->setSizes({430, 240});
    el->addWidget(es, 1);
    m_tabs->addTab(elastic, tr("3 弹性标定"));
    connect(m_fits, &QTableWidget::itemSelectionChanged, this, &SpectroscopyDialog::refreshElasticCurve);
    connect(m_elasticView, &QComboBox::currentIndexChanged, this, [this] { refreshElasticCurve(); });

    auto xrs = new QWidget;
    auto xl = new QVBoxLayout(xrs);
    auto choices = new QHBoxLayout;
    m_scanChoice = new QComboBox;
    m_crystalChoice = new QComboBox;
    for (const auto &roi : m_rois.items())
        m_crystalChoice->addItem(roi.label);
    choices->addWidget(new QLabel(tr("扫描")));
    choices->addWidget(m_scanChoice, 1);
    choices->addWidget(new QLabel(tr("晶体")));
    choices->addWidget(m_crystalChoice);
    xl->addLayout(choices);
    auto xs = new QSplitter(Qt::Vertical);
    m_i0Plot = new SpectrumPlot;
    m_i0Plot->setObjectName("I0Spectrum");
    m_scanPlot = new SpectrumPlot;
    m_scanPlot->setObjectName("XrsSpectrum");
    xs->addWidget(m_i0Plot);
    xs->addWidget(m_scanPlot);
    xl->addWidget(xs, 1);
    m_tabs->addTab(xrs, tr("4 XRS 与 I0"));
    connect(m_scanChoice, &QComboBox::currentIndexChanged, this, [this] {
        refreshScanCurves();
        refreshImages();
    });
    connect(m_crystalChoice, &QComboBox::currentIndexChanged, this, [this] { refreshScanCurves(); });

    auto resultPage = new QWidget;
    auto resLayout = new QHBoxLayout(resultPage);
    auto options = new QWidget;
    auto ol = new QVBoxLayout(options);
    auto geometry = new QGroupBox(tr("散射几何"));
    auto gf = new QFormLayout(geometry);
    m_angleSource = new QComboBox;
    m_angleSource->addItems({tr("从弹性参考读取"), tr("从首个 XRS 扫描读取"), tr("手动角度（degree）")});
    gf->addRow(tr("模组角度"), m_angleSource);
    for (const auto &module : {"VB", "HL", "VU", "HB", "VD"}) {
        auto b = realBox(90, -360, 360, 4);
        m_angles[module] = b;
        b->setEnabled(false);
        gf->addRow(module, b);
    }
    m_diameter = realBox(100, .001, 10000, 3);
    m_distance = realBox(1070, .001, 100000, 3);
    m_pitch = realBox(7, 0, 90, 4);
    gf->addRow(tr("晶体口径 / mm"), m_diameter);
    gf->addRow(tr("样品距离 / mm"), m_distance);
    gf->addRow(tr("晶体间隔 / degree"), m_pitch);
    ol->addWidget(geometry);
    auto selection = new QGroupBox(tr("筛选与插值"));
    auto sf = new QFormLayout(selection);
    auto moduleRow = new QWidget;
    auto mr = new QHBoxLayout(moduleRow);
    mr->setContentsMargins(0, 0, 0, 0);
    for (const auto &module : {"VB", "HL", "VU", "HB", "VD"}) {
        auto c = new QCheckBox(module);
        c->setChecked(true);
        m_modules.append(c);
        mr->addWidget(c);
    }
    sf->addRow(moduleRow);
    m_qMin = realBox(0, 0, 100, 5);
    m_qMax = realBox(10, 0, 100, 5);
    m_step = realBox(.2, .000001, 10000, 6);
    sf->addRow(tr("Q 下限 / Å⁻¹"), m_qMin);
    sf->addRow(tr("Q 上限 / Å⁻¹"), m_qMax);
    sf->addRow(tr("能量步长 / eV"), m_step);
    auto resultNote = new QLabel(
        tr("能量转移 = 实测能量 − 各晶体弹性峰中心。仅在合格曲线的共同范围内线性插值；扫描与晶体均求和。"));
    resultNote->setWordWrap(true);
    sf->addRow(resultNote);
    ol->addWidget(selection);
    ol->addStretch();
    auto optionScroll = scrollPanel(options);
    optionScroll->setMinimumWidth(300);
    optionScroll->setMaximumWidth(390);
    resLayout->addWidget(optionScroll);
    auto resultBody = new QWidget;
    auto rb = new QVBoxLayout(resultBody);
    m_resultChoice = new QComboBox;
    rb->addWidget(m_resultChoice);
    auto rs = new QSplitter(Qt::Vertical);
    m_resultPlot = new SpectrumPlot;
    m_resultPlot->setObjectName("CombinedSpectrum");
    rs->addWidget(m_resultPlot);
    auto inclusionBody = new QWidget;
    auto il = new QVBoxLayout(inclusionBody);
    auto batch = new QHBoxLayout;
    batch->addWidget(new QLabel(tr("勾选允许叠加的晶体（Ctrl / Shift 可多选）")));
    batch->addStretch();
    auto include = new QPushButton(tr("纳入选中"));
    auto exclude = new QPushButton(tr("排除选中"));
    batch->addWidget(include);
    batch->addWidget(exclude);
    il->addLayout(batch);
    m_inclusion = new QTableWidget(75, 7);
    m_inclusion->setObjectName("SpectrumInclusion");
    configureTable(m_inclusion, {tr("晶体"), tr("Q 均值 / Å⁻¹"), tr("dQ 均值 / Å⁻¹"), tr("Q 变化"),
                                 tr("dQ 变化"), tr("叠加"), tr("原因")});
    il->addWidget(m_inclusion);
    rs->addWidget(inclusionBody);
    rs->setSizes({400, 250});
    rb->addWidget(rs, 1);
    resLayout->addWidget(resultBody, 1);
    m_tabs->addTab(resultPage, tr("5 叠加与结果"));
    auto batchSet = [this](Qt::CheckState state) {
        bool changed = false;
        for (const auto &index : m_inclusion->selectionModel()->selectedRows()) {
            auto item = m_inclusion->item(index.row(), 0);
            if (item && (item->flags() & Qt::ItemIsUserCheckable)) {
                QSignalBlocker b(m_inclusion);
                item->setCheckState(state);
                changed = true;
            }
        }
        if (changed)
            invalidate(2);
    };
    connect(include, &QPushButton::clicked, this, [batchSet] { batchSet(Qt::Checked); });
    connect(exclude, &QPushButton::clicked, this, [batchSet] { batchSet(Qt::Unchecked); });
    connect(m_inclusion, &QTableWidget::itemChanged, this, [this](QTableWidgetItem *item) {
        if (item->column() == 0)
            invalidate(2);
    });
    connect(m_resultChoice, &QComboBox::currentIndexChanged, this, [this] { refreshResultCurve(); });
    m_status = new QLabel(tr("选择弹性扫描，载入共用 ROI 方案，然后开始弹性标定。"));
    m_status->setObjectName("SpectrumStatus");
    m_status->setWordWrap(true);
    m_status->setTextInteractionFlags(Qt::TextSelectableByMouse);
    root->addWidget(m_status);
    m_progress = new QProgressBar;
    m_progress->hide();
    root->addWidget(m_progress);
    auto actions = new QHBoxLayout;
    m_calibrateButton = new QPushButton(tr("计算弹性标定"));
    m_calibrateButton->setObjectName("CalibrateSpectrum");
    m_extractButton = new QPushButton(tr("提取 XRS"));
    m_extractButton->setObjectName("ExtractSpectrum");
    m_combineButton = new QPushButton(tr("计算 Q 并叠加"));
    m_combineButton->setObjectName("CombineSpectrum");
    m_saveButton = new QPushButton(tr("保存谱学结果…"));
    m_saveButton->setObjectName("SaveSpectrum");
    m_cancelButton = new QPushButton(tr("取消计算"));
    m_previousButton = new QPushButton(tr("上一步"));
    m_previousButton->setObjectName("PreviousSpectrumStep");
    m_nextButton = new QPushButton(tr("下一步"));
    m_nextButton->setObjectName("NextSpectrumStep");
    actions->addWidget(m_previousButton);
    actions->addStretch();
    actions->addWidget(m_calibrateButton);
    actions->addWidget(m_extractButton);
    actions->addWidget(m_combineButton);
    actions->addWidget(m_saveButton);
    actions->addWidget(m_nextButton);
    actions->addWidget(m_cancelButton);
    auto close = new QPushButton(tr("关闭"));
    close->setObjectName("CloseWorkflow");
    actions->addWidget(close);
    root->addLayout(actions);
    connect(m_tabs, &QTabWidget::currentChanged, this, [this] { refreshActions(); });
    connect(m_nextButton, &QPushButton::clicked, this, [this] { selectStep(m_tabs->currentIndex() + 1); });
    connect(m_previousButton, &QPushButton::clicked, this,
            [this] { selectStep(m_tabs->currentIndex() - 1); });
    connect(close, &QPushButton::clicked, this, &QWidget::hide);
    connect(m_calibrateButton, &QPushButton::clicked, this, [this] { runStep(0); });
    connect(m_extractButton, &QPushButton::clicked, this, [this] { runStep(1); });
    connect(m_combineButton, &QPushButton::clicked, this, [this] { runStep(2); });
    connect(m_saveButton, &QPushButton::clicked, this, &SpectroscopyDialog::save);
    connect(m_cancelButton, &QPushButton::clicked, this, [this] {
        m_cancel.store(true);
        message(tr("正在取消，等待当前帧读取结束…"));
    });
    // Invalidate only the affected stage. No stale calibration or result may be saved.
    connect(m_elasticFile, &QLineEdit::textChanged, this, [this] {
        m_rois = CrystalRois();
        m_roiPath->clear();
        invalidate(0);
    });
    for (auto c : {m_elasticEnergy, m_elasticI0, m_model})
        connect(c, &QComboBox::currentTextChanged, this, [this] { invalidate(0); });
    for (auto c : {m_xrsEnergy, m_xrsI0})
        connect(c, &QComboBox::currentTextChanged, this, [this] { invalidate(1); });
    for (auto c : {m_normalizeElastic, m_filter, m_auto, m_overlap})
        connect(c, &QCheckBox::toggled, this, [this] { invalidate(0); });
    for (auto c : {m_normalizeXrs, m_repair})
        connect(c, &QCheckBox::toggled, this, [this] { invalidate(1); });
    for (auto b : {m_fraction, m_centerMin, m_centerMax, m_fwhm, m_r2})
        connect(b, &QDoubleSpinBox::valueChanged, this, [this] { invalidate(0); });
    connect(m_roiSize, &QSpinBox::valueChanged, this, [this] { invalidate(0); });
    connect(m_corrections, &QTableWidget::itemChanged, this, [this] { invalidate(0); });
    connect(m_angleSource, &QComboBox::currentIndexChanged, this, [this] {
        for (auto b : m_angles)
            b->setEnabled(m_angleSource->currentIndex() == 2);
        invalidate(2);
    });
    for (auto b : m_angles)
        connect(b, &QDoubleSpinBox::valueChanged, this, [this] { invalidate(2); });
    for (auto b : {m_diameter, m_distance, m_pitch, m_qMin, m_qMax, m_step})
        connect(b, &QDoubleSpinBox::valueChanged, this, [this] { invalidate(2); });
    for (auto c : m_modules)
        connect(c, &QCheckBox::toggled, this, [this] { invalidate(2); });
}

SpectrumSettings SpectroscopyDialog::settings() const {
    SpectrumSettings s;
    s.elasticFile = m_elasticFile->text().trimmed();
    s.elasticEnergy = m_elasticEnergy->currentText().trimmed();
    s.xrsEnergy = m_xrsEnergy->currentText().trimmed();
    s.elasticI0 = m_elasticI0->currentText().trimmed();
    s.xrsI0 = m_xrsI0->currentText().trimmed();
    s.normalizeElastic = m_normalizeElastic->isChecked();
    s.normalizeXrs = m_normalizeXrs->isChecked();
    s.repairGlitches = m_repair->isChecked();
    s.useFilter = m_filter->isChecked();
    s.autoRoi = m_auto->isChecked();
    s.removeOverlap = m_overlap->isChecked();
    s.filterFraction = m_fraction->value();
    s.roiSize = m_roiSize->value();
    s.fitModel = m_model->currentIndex();
    s.centerMin = m_centerMin->value();
    s.centerMax = m_centerMax->value();
    s.maxFwhmEv = m_fwhm->value();
    s.minRSquared = m_r2->value();
    for (int r = 0; r < 75; ++r) {
        int values[4];
        for (int c = 1; c < 5; ++c) {
            bool ok = false;
            values[c - 1] = m_corrections->item(r, c)->text().toInt(&ok);
            if (!ok || std::abs(double(values[c - 1])) > 2048 || (c >= 3 && values[c - 1] < 0))
                throw std::runtime_error(
                    tr("ROI 修正必须为整数，偏移范围 ±2048，扩展范围 0 至 2048。").toUtf8().constData());
        }
        s.shifts[r] = {values[0], values[1]};
        s.expansions[r] = {values[2], values[3]};
        auto item = m_inclusion->item(r, 0);
        if (item && (item->flags() & Qt::ItemIsUserCheckable) && item->checkState() != Qt::Checked)
            s.excludedCrystals.insert(r);
    }
    s.angleSource = m_angleSource->currentIndex();
    for (auto it = m_angles.begin(); it != m_angles.end(); ++it)
        s.angles[it.key()] = it.value()->value();
    s.crystalDiameter = m_diameter->value();
    s.crystalDistance = m_distance->value();
    s.modulePitchDegrees = m_pitch->value();
    s.qMin = m_qMin->value();
    s.qMax = m_qMax->value();
    s.energyStepEv = m_step->value();
    s.modules.clear();
    for (auto c : m_modules)
        if (c->isChecked())
            s.modules.append(c->text());
    return s;
}
void SpectroscopyDialog::invalidate(int stage) {
    if (m_updating || m_busy)
        return;
    if (stage <= 0) {
        m_haveCalibration = false;
        m_calibration = ElasticCalibration();
        m_fits->clearContents();
        m_inclusion->clearContents();
        m_elasticPlot->setCurves({}, {});
        for (auto image : m_images)
            image->setData({});
    }
    if (stage <= 1) {
        m_haveScans = false;
        m_scans.clear();
        {
            QSignalBlocker b(m_scanChoice);
            m_scanChoice->clear();
        }
        m_i0Plot->setCurves({}, {});
        m_scanPlot->setCurves({}, {});
        refreshImages();
    }
    m_haveResult = false;
    m_result = SpectrumResult();
    m_resultPlot->setCurves({}, {});
    m_resultChoice->clear();
    if (stage == 2) {
        QSignalBlocker b(m_inclusion);
        for (int r = 0; r < 75; ++r)
            for (int c = 1; c < 7; ++c)
                if (m_inclusion->item(r, c))
                    m_inclusion->item(r, c)->setText(QStringLiteral("—"));
    }
    refreshState();
    message(stage == 0   ? tr("弹性参数或 ROI 已变更，请重新标定。")
            : stage == 1 ? tr("XRS 扫描或归一化参数已变更，请重新提取。")
                         : tr("叠加条件已变更，请重新计算 Q 并叠加。"));
}
void SpectroscopyDialog::refreshState() {
    if (!m_busy && !m_pendingRoi.second.isEmpty()) {
        const auto pending = m_pendingRoi;
        m_pendingRoi = {};
        useRoiDraft(pending.first, pending.second);
    }
    refreshActions();
    m_calibrateButton->setEnabled(!m_busy);
    m_extractButton->setEnabled(!m_busy && m_haveCalibration);
    m_combineButton->setEnabled(!m_busy && m_haveScans);
    m_saveButton->setEnabled(!m_busy && m_haveResult);
    m_cancelButton->setEnabled(m_busy);
    m_roiStatus->setText(tr("有效 %1 · 失效 %2 · 待设置 %3 / 75")
                             .arg(m_rois.count(RoiState::Valid))
                             .arg(m_rois.count(RoiState::Invalid))
                             .arg(m_rois.count(RoiState::Pending)));
}
bool SpectroscopyDialog::loadRois(const QString &path) {
    if (m_busy)
        return false;
    try {
        QVector<QSize> sizes;
        for (auto detector : {"D_LAMBDA", "D_MINIPIX"}) {
            auto n = HdfFile::inspect(m_elasticFile->text(), QString("/entry/data/") + detector);
            if (n.dims.size() != 3)
                throw std::runtime_error("Detector must be a 3D stack");
            sizes.append({int(n.dims[2]), int(n.dims[1])});
        }
        QFile file(path);
        QString error;
        CrystalRois rois;
        if (!file.open(QIODevice::ReadOnly))
            throw std::runtime_error(file.errorString().toUtf8().constData());
        if (!rois.restore(file.readAll(), sizes, &error))
            throw std::runtime_error(error.toUtf8().constData());
        bool changed = rois.serialize({}, sizes) != m_rois.serialize({}, sizes);
        m_rois = rois;
        m_roiPath->setText(QFileInfo(path).absoluteFilePath());
        if (!changed) {
            refreshState();
            return true;
        }
        invalidate(0);
        message(tr("已载入 ROI：%1；待设置的晶体须先框选或标记失效。").arg(path));
        return true;
    } catch (const std::exception &e) {
        message(QString::fromUtf8(e.what()), true);
        return false;
    }
}
void SpectroscopyDialog::loadDraft() {
    auto canonical = QFileInfo(m_elasticFile->text()).canonicalFilePath();
    if (canonical.isEmpty()) {
        message(tr("弹性文件不存在。"), true);
        return;
    }
    auto path = QStandardPaths::writableLocation(QStandardPaths::AppLocalDataLocation) + "/roi-drafts/" +
                QCryptographicHash::hash(canonical.toUtf8(), QCryptographicHash::Sha256).toHex() + ".json";
    if (QFileInfo::exists(path))
        loadRois(path);
    else
        message(tr("此弹性文件尚无 ROI 草稿。请载入已有方案，或使用“选取 ROI”。"));
}
void SpectroscopyDialog::addScans(const QStringList &files) {
    if (m_busy)
        return;
    for (const auto &file : files) {
        auto path = QFileInfo(file).absoluteFilePath();
        bool duplicate = false;
        for (int i = 0; i < m_files->count(); ++i)
            if (QFileInfo(m_files->item(i)->data(Qt::UserRole).toString()).canonicalFilePath() ==
                QFileInfo(path).canonicalFilePath())
                duplicate = true;
        if (duplicate)
            continue;
        auto item = new QListWidgetItem(QFileInfo(path).fileName());
        item->setToolTip(path);
        item->setData(Qt::UserRole, path);
        item->setFlags(item->flags() | Qt::ItemIsUserCheckable);
        item->setCheckState(Qt::Checked);
        m_files->addItem(item);
    }
    invalidate(1);
}
void SpectroscopyDialog::runStep(int step) {
    if (m_busy || step < 0 || step > 2)
        return;
    if ((step == 1 && !m_haveCalibration) || (step == 2 && !m_haveScans)) {
        message(tr("请按顺序完成前面的处理步骤。"), true);
        return;
    }
    try {
        auto s = settings();
        QStringList files;
        for (int i = 0; i < m_files->count(); ++i)
            if (m_files->item(i)->checkState() == Qt::Checked)
                files.append(m_files->item(i)->data(Qt::UserRole).toString());
        if (step == 1 && files.isEmpty()) {
            message(tr("请先在“文件与参数”中添加并勾选 XRS 扫描。"), true);
            m_tabs->setCurrentIndex(0);
            return;
        }
        invalidate(step);
        m_job = step;
        if (step == 0)
            m_calibrationSettings = s;
        if (step == 1)
            m_extractionSettings = s;
        if (step == 2)
            m_resultSettings = s;
        m_busy = true;
        m_cancel.store(false);
        m_tabs->setEnabled(false);
        m_progress->setRange(0, 100);
        m_progress->setValue(0);
        m_progress->show();
        refreshState();
        message(tr("正在计算，可取消…"));
        auto rois = m_rois;
        auto calibration = m_calibration;
        auto scans = m_scans;
        auto progress = [this](const QString &text, int value) {
            QMetaObject::invokeMethod(
                this,
                [this, text, value] {
                    if (m_busy) {
                        m_progress->setValue(value);
                        message(text);
                    }
                },
                Qt::QueuedConnection);
        };
        m_watcher.setFuture(QtConcurrent::run([this, step, s, rois, calibration, scans, files, progress] {
            Work out;
            try {
                if (step == 0)
                    out.calibration = Spectroscopy::calibrate(s, rois, m_cancel, progress);
                if (step == 1)
                    out.scans = Spectroscopy::extract(s, calibration, files, m_cancel, progress);
                if (step == 2)
                    out.result = Spectroscopy::combine(s, calibration, scans, m_cancel);
            } catch (const std::exception &e) {
                out.error = QString::fromUtf8(e.what());
            }
            return out;
        }));
    } catch (const std::exception &e) {
        message(QString::fromUtf8(e.what()), true);
    }
}
void SpectroscopyDialog::refreshCalibration() {
    m_updating = true;
    {
        QSignalBlocker fitBlock(m_fits), inclusionBlock(m_inclusion);
        for (int r = 0; r < 75; ++r) {
            const auto &roi = m_calibration.adjusted.items()[r];
            const auto &fit = m_calibration.fits[r];
            QStringList values{roi.label,
                               QString("%1 / %2")
                                   .arg(m_calibration.pixels[r].size())
                                   .arg(m_calibration.filteredPixels[r].size()),
                               number(fit.centerKev * 1000),
                               number(fit.fwhmEv),
                               number(fit.rSquared),
                               fit.accepted ? tr("合格") : tr("未通过"),
                               fit.reason};
            for (int c = 0; c < values.size(); ++c)
                m_fits->setItem(r, c, cell(values[c]));
            auto item = cell(roi.label);
            if (fit.accepted) {
                item->setFlags(item->flags() | Qt::ItemIsUserCheckable);
                item->setCheckState(Qt::Checked);
            }
            m_inclusion->setItem(r, 0, item);
            for (int c = 1; c < 7; ++c)
                m_inclusion->setItem(r, c, cell(c == 6 ? fit.reason : QStringLiteral("—")));
        }
        m_fits->resizeColumnsToContents();
        const int minimum[]{75, 105, 135, 120, 115, 75, 260};
        for (int col = 0; col < 7; ++col)
            m_fits->setColumnWidth(col, std::max(minimum[col], m_fits->columnWidth(col)));
        int selected = 0;
        for (int row = 0; row < 75; ++row)
            if (m_calibration.fits[row].accepted) {
                selected = row;
                break;
            }
        m_fits->setCurrentCell(selected, 0);
    }
    m_updating = false;
    refreshImages();
    refreshElasticCurve();
}
void SpectroscopyDialog::refreshImages() {
    if (!m_images[0])
        return;
    int choice = m_imageSource->currentIndex();
    if (m_haveCalibration) {
        if (m_imageSource->count() == 2)
            m_imageSource->addItem(tr("当前 XRS 扫描叠加图"));
        for (int d = 0; d < 2; ++d) {
            DataSlice preview;
            if (choice == 0)
                preview = m_calibration.images[d];
            else if (choice == 1)
                preview = m_calibration.maskImages[d];
            else if (m_haveScans && m_scanChoice->currentIndex() >= 0)
                preview = m_scans[m_scanChoice->currentIndex()].images[d];
            m_images[d]->setData(preview);
            QVector<CrystalRoi> overlays;
            for (const auto &roi : m_calibration.adjusted.items())
                if (roi.detector == (d == 0 ? "D_LAMBDA" : "D_MINIPIX"))
                    overlays.append(roi);
            m_images[d]->setRoiOverlays(overlays, {});
        }
    } else
        for (auto image : m_images)
            image->setData({});
}
void SpectroscopyDialog::refreshElasticCurve() {
    if (!m_haveCalibration)
        return;
    QVector<SpectrumCurve> curves;
    int row = std::max(0, m_fits->currentRow());
    const auto &scan = m_calibration.scan;
    if (m_elasticView->currentIndex() == 0) {
        const auto &fit = m_calibration.fits[row];
        auto label = m_rois.items()[row].label;
        curves.append({tr("原 ROI"), scan.energyKev, scan.raw.value(row), QColor(100, 100, 100), false});
        curves.append({tr("强度掩膜"), scan.energyKev, scan.filtered.value(row), QColor(0, 100, 180), false});
        if (!fit.fitted.isEmpty())
            curves.append({tr("拟合"), scan.energyKev, fit.fitted, QColor(190, 45, 35), true});
        m_elasticPlot->setCurves(curves, label + tr(" · 弹性峰"), "Energy (keV)",
                                 m_calibrationSettings.normalizeElastic ? "Intensity / |I0|" : "Counts");
    } else {
        for (int r = 0; r < 75; ++r)
            if (m_calibration.fits[r].accepted) {
                QVector<double> x;
                for (double energy : scan.energyKev)
                    x.append((energy - m_calibration.fits[r].centerKev) * 1000);
                curves.append({m_rois.items()[r].label,
                               x,
                               m_calibrationSettings.useFilter ? scan.filtered[r] : scan.raw[r],
                               {},
                               false});
            }
        m_elasticPlot->setCurves(curves, tr("合格晶体 · 弹性零点对齐"), "Energy transfer (eV)",
                                 m_calibrationSettings.normalizeElastic ? "Intensity / |I0|" : "Counts");
    }
}
void SpectroscopyDialog::refreshScanCurves() {
    if (!m_haveScans || m_scanChoice->currentIndex() < 0)
        return;
    const auto &scan = m_scans[m_scanChoice->currentIndex()];
    int r = m_crystalChoice->currentIndex();
    m_i0Plot->setCurves({{tr("修复前 |I0|"), scan.energyKev, scan.i0Before, QColor(120, 120, 120), false},
                         {tr("修复后 |I0|"), scan.energyKev, scan.i0After, QColor(185, 45, 35), true}},
                        tr("I0 · 异常点 %1 · %2")
                            .arg(scan.glitches.size())
                            .arg(scan.i0Path.isEmpty() ? tr("未归一化，I0 = 1") : scan.i0Path),
                        "Energy (keV)", "|I0| (a.u.)");
    m_scanPlot->setCurves(
        {{tr("原 ROI"), scan.energyKev, scan.raw.value(r), QColor(100, 100, 100), false},
         {tr("强度掩膜"), scan.energyKev, scan.filtered.value(r), QColor(0, 100, 180), false}},
        m_rois.items()[r].label + " · " + QFileInfo(scan.file).fileName(), "Energy (keV)",
        m_extractionSettings.normalizeXrs ? "Intensity / |I0|" : "Counts");
}
void SpectroscopyDialog::refreshResult() {
    m_updating = true;
    {
        QSignalBlocker block(m_inclusion);
        for (int r = 0; r < 75; ++r) {
            const auto &q = m_result.momentum[r];
            QStringList values{number(q.mean),
                               number(q.dMean),
                               number(q.range),
                               number(q.dRange),
                               m_result.included[r] ? tr("纳入") : tr("排除"),
                               m_result.reasons[r]};
            for (int c = 0; c < 6; ++c)
                m_inclusion->item(r, c + 1)->setText(values[c]);
        }
        m_inclusion->resizeColumnsToContents();
        m_inclusion->setColumnWidth(6, 240);
    }
    {
        QSignalBlocker block(m_resultChoice);
        m_resultChoice->clear();
        m_resultChoice->addItem(tr("总谱与各模组"), -1);
        m_resultChoice->addItem(tr("全部单晶体谱"), -2);
        for (int r = 0; r < 75; ++r)
            if (m_result.included[r])
                m_resultChoice->addItem(m_rois.items()[r].label + tr(" · 各扫描与叠加"), r);
    }
    if (m_angleSource->currentIndex() != 2)
        for (auto it = m_result.angles.begin(); it != m_result.angles.end(); ++it)
            if (m_angles.contains(it.key())) {
                QSignalBlocker block(m_angles[it.key()]);
                m_angles[it.key()]->setValue(it.value());
            }
    m_updating = false;
    refreshResultCurve();
}
void SpectroscopyDialog::refreshResultCurve() {
    if (!m_haveResult)
        return;
    QVector<SpectrumCurve> curves;
    int choice = m_resultChoice->currentData().toInt();
    if (choice == -1) {
        curves.append({tr("Total"), m_result.energyTransferEv, m_result.total, Qt::black, false});
        for (auto it = m_result.modules.begin(); it != m_result.modules.end(); ++it)
            curves.append({it.key(), m_result.energyTransferEv, it.value(), {}, false});
    } else if (choice == -2) {
        for (int r = 0; r < 75; ++r)
            if (m_result.included[r])
                curves.append(
                    {m_rois.items()[r].label, m_result.energyTransferEv, m_result.crystals[r], {}, false});
    } else if (choice >= 0 && choice < 75 && m_result.included[choice]) {
        curves.append(
            {tr("扫描叠加"), m_result.energyTransferEv, m_result.crystals[choice], Qt::black, false});
        for (const auto &scan : m_scans) {
            QVector<double> x;
            for (double energy : scan.energyKev)
                x.append((energy - m_calibration.fits[choice].centerKev) * 1000);
            curves.append({QFileInfo(scan.file).fileName(),
                           x,
                           m_resultSettings.useFilter ? scan.filtered[choice] : scan.raw[choice],
                           {},
                           true});
        }
    }
    m_resultPlot->setCurves(curves, m_resultChoice->currentText(), "Energy transfer (eV)",
                            m_resultSettings.normalizeXrs ? "Summed intensity / |I0|" : "Summed counts");
}
void SpectroscopyDialog::save() {
    if (!m_haveResult || m_busy)
        return;
    auto parent = QFileDialog::getExistingDirectory(this, tr("选择谱学结果的父目录"),
                                                    QFileInfo(m_elasticFile->text()).absolutePath());
    if (parent.isEmpty())
        return;
    auto directory =
        QDir(parent).filePath("spectrum_" + QDateTime::currentDateTime().toString("yyyyMMdd_HHmmss_zzz"));
    QString error;
    if (Spectroscopy::saveResult(directory, m_resultSettings, m_calibration, m_scans, m_result, &error))
        message(tr("已保存总谱、模组谱、单晶体谱、ROI、掩膜和完整参数：%1").arg(directory));
    else
        message(error, true);
}

void SpectroscopyDialog::refreshActions() {
    int step = m_tabs->currentIndex();
    m_previousButton->setVisible(step > 0);
    m_previousButton->setEnabled(!m_busy);
    m_nextButton->setVisible(step < 4);
    m_nextButton->setEnabled(!m_busy && (step < 2 || (step == 2 ? m_haveCalibration : m_haveScans)));
    m_calibrateButton->setVisible(step == 1 || step == 2);
    m_extractButton->setVisible(step == 3);
    m_combineButton->setVisible(step == 4);
    m_saveButton->setVisible(step == 4);
    m_cancelButton->setVisible(m_busy);
    m_nextButton->setDefault(step < 2);
    m_calibrateButton->setDefault((step == 1 || step == 2) && !m_haveCalibration);
    m_extractButton->setDefault(step == 3 && !m_haveScans);
    m_combineButton->setDefault(step == 4 && !m_haveResult);
    m_saveButton->setDefault(step == 4 && m_haveResult);
}
void SpectroscopyDialog::useRoiDraft(const QString &file, const QString &path) {
    if (m_busy) {
        m_pendingRoi = {file, path};
        return;
    }
    if (QFileInfo(m_elasticFile->text()).canonicalFilePath() == file && QFileInfo::exists(path))
        loadRois(path);
}
