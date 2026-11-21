#include "SpectroscopyDialog.h"
#include "ImagingIO.h"
#include "RoiDialog.h"
#include "ScientificPlot.h"
#include "SpectralMaskDialog.h"
#include "SpectrumPlot.h"
#include "SpectrumProcessingPage.h"
#include "WorkflowUi.h"
#include "WorkspacePaths.h"
#include <QJsonArray>
#include <QJsonDocument>
#include <QMetaProperty>
#include <QUuid>
#include <QtConcurrent/QtConcurrentRun>
#include <QtWidgets>

namespace {
using WorkspacePaths::Kind;
using WorkspacePaths::Workflow;

constexpr int prepareMaskStep = -1;
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
QVector<QSize> detectorSizes(const QString &elasticFile) {
    QVector<QSize> sizes;
    for (auto detector : {"D_LAMBDA", "D_MINIPIX"}) {
        auto node = HdfFile::inspect(elasticFile, QString("/entry/data/") + detector);
        if (node.dims.size() != 3)
            throw std::runtime_error("Detector must be a 3D stack");
        sizes.append({int(node.dims[2]), int(node.dims[1])});
    }
    return sizes;
}
QString absolutePath(const QDir &base, const QString &path) {
    return path.isEmpty() ? QString{} : QDir::cleanPath(base.absoluteFilePath(path));
}
QString fileToken(QString text, const QString &fallback) {
    static const QRegularExpression unsupported("[^A-Za-z0-9_-]+");
    static const QRegularExpression trim("_+$|^_+");
    text = text.replace(unsupported, "_").left(40).remove(trim);
    return text.isEmpty() ? fallback : text;
}
QStringList edgeElements() {
    QFile table(effectiveHfTable({}));
    if (!table.open(QIODevice::ReadOnly))
        return {};
    QTextStream stream(&table);
    QStringList elements;
    while (!stream.atEnd()) {
        const auto line = stream.readLine().simplified();
        if (line.startsWith("#S "))
            elements.append(line.section(' ', -1));
    }
    elements.sort(Qt::CaseInsensitive);
    return elements;
}

QString edgeGroupName(const QString &element, const QString &shell) {
    return element + '_' + shell;
}
QMetaProperty configProperty(QWidget *editor) {
    if (auto combo = qobject_cast<QComboBox *>(editor)) {
        auto meta = combo->metaObject();
        return meta->property(meta->indexOfProperty(combo->isEditable() ? "currentText" : "currentIndex"));
    }
    return editor->metaObject()->userProperty();
}
bool validConfigValue(QWidget *editor, const QJsonValue &value) {
    if (value.type() != QJsonValue::fromVariant(configProperty(editor).read(editor)).type())
        return false;
    if (auto spin = qobject_cast<QDoubleSpinBox *>(editor))
        return value.toDouble() >= spin->minimum() && value.toDouble() <= spin->maximum();
    if (auto spin = qobject_cast<QSpinBox *>(editor))
        return value.toInt(-1) >= spin->minimum() && value.toInt(-1) <= spin->maximum();
    if (auto combo = qobject_cast<QComboBox *>(editor); combo && !combo->isEditable())
        return value.toInt(-1) >= 0 && value.toInt(-1) < combo->count();
    return true;
}
void appendScan(QListWidget *list, const QString &path, bool enabled) {
    auto item = new QListWidgetItem(QFileInfo(path).fileName());
    item->setToolTip(path);
    item->setData(Qt::UserRole, path);
    item->setFlags(item->flags() | Qt::ItemIsUserCheckable);
    item->setCheckState(enabled ? Qt::Checked : Qt::Unchecked);
    list->addItem(item);
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
            m_nextPage = -1;
            message(work.error, true);
            refreshState();
            return;
        }
        if (m_job <= 0) {
            m_masks = std::move(work.masks);
            m_projection = m_masks.images;
        }
        if (m_job == prepareMaskStep) {
            m_imageSource->setCurrentIndex(1);
            refreshImages();
            selectStep(1);
            message(tr("掩膜已计算。可修正像素，或点击下一步计算弹性标定。"));
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
            setEdgeStatus(EdgeStatus::Combined);
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
        if (m_nextPage >= 0)
            advanceToPage(m_nextPage);
    });
    if (!initialFile.isEmpty() && HdfFile::formatFromFileName(initialFile) != ScanFormat::GridScan) {
        m_elasticFile->setText(QFileInfo(initialFile).absoluteFilePath());
        loadDraft();
    }
    refreshState();
}
bool SpectroscopyDialog::isBusy() const {
    return m_busy || (m_processing && m_processing->isBusy());
}

SpectroscopyDialog::~SpectroscopyDialog() {
    m_cancel.store(true);
    m_watcher.waitForFinished();
}
void SpectroscopyDialog::selectStep(int step) {
    m_tabs->setCurrentIndex(std::clamp(step, 0, m_tabs->count() - 1));
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
    WorkflowUi::compact(root);
    m_tabs = new QTabWidget;
    m_tabs->setObjectName("SpectroscopySteps");
    m_tabs->setUsesScrollButtons(true);
    auto importConfig = WorkflowUi::button("folder-open", tr("载入已有方案…"));
    importConfig->setObjectName("ImportSpectrumConfig");
    m_saveConfig = WorkflowUi::button("save", tr("保存当前配置…"));
    m_saveConfig->setObjectName("SaveSpectrumConfig");
    m_tabs->setCornerWidget(WorkflowUi::cornerTools({importConfig, m_saveConfig}));
    connect(importConfig, &QPushButton::clicked, this, [this] {
        const auto path = QFileDialog::getOpenFileName(
            this, tr("载入已有方案"), WorkspacePaths::directory(Workflow::Spectroscopy, Kind::Config),
            "JSON (*.json)");
        if (!path.isEmpty())
            importConfiguration(path);
    });
    connect(m_saveConfig, &QPushButton::clicked, this, [this] {
        const auto name = sampleStem() + "_spectroscopy_groups.config.json";
        const auto path = WorkspacePaths::saveFile(this, Workflow::Spectroscopy, Kind::Config, {},
                                                   tr("保存谱学配置"), name, "JSON (*.json)");
        if (!path.isEmpty())
            saveConfiguration(path);
    });
    root->addWidget(m_tabs, 1);
    auto inputs = new QWidget;
    auto inputLayout = new QHBoxLayout(inputs);
    WorkflowUi::compact(inputLayout);
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
    auto load = WorkflowUi::button("folder-open", tr("载入 ROI…"));
    auto edit = WorkflowUi::button("scan", tr("选取 ROI…"));
    auto roiButtons = new QWidget;
    auto rl = new QHBoxLayout(roiButtons);
    rl->setContentsMargins(0, 0, 0, 0);
    rl->addWidget(m_roiPath, 1);
    rl->addWidget(load);
    rl->addWidget(edit);
    fileForm->addRow(tr("ROI 方案"), roiButtons);
    m_roiStatus = new QLabel;
    m_roiStatus->setWordWrap(true);
    fileForm->addRow(m_roiStatus);
    m_hfComposition = new QLineEdit;
    m_hfComposition->setObjectName("BackgroundHfComposition");
    m_hfComposition->setPlaceholderText(tr("例如 Li:1, Ni:0.8, Co:0.1, Mn:0.1, O:2"));
    fileForm->addRow(tr("元素与化学比"), m_hfComposition);
    auto hfRow = new QWidget;
    auto hfLayout = new QHBoxLayout(hfRow);
    hfLayout->setContentsMargins(0, 0, 0, 0);
    m_hfTable = new QLineEdit;
    m_hfTable->setObjectName("BackgroundHfTable");
    m_hfTable->setPlaceholderText(tr("内置默认表 · ComptonProfiles.dat"));
    auto chooseTable = WorkflowUi::button("folder-open", tr("选择…"));
    auto defaultTable = new QPushButton(tr("恢复默认"));
    defaultTable->setObjectName("RestoreDefaultHfTable");
    defaultTable->setAutoDefault(false);
    defaultTable->setEnabled(false);
    hfLayout->addWidget(m_hfTable, 1);
    hfLayout->addWidget(chooseTable);
    hfLayout->addWidget(defaultTable);
    fileForm->addRow(tr("HF 壳层表"), hfRow);
    connect(chooseTable, &QPushButton::clicked, this, [this] {
        const auto path = QFileDialog::getOpenFileName(this, tr("选择 DABAX HF 壳层数据"),
                                                       m_hfTable->text(), "DABAX (*.dat);;All files (*)");
        if (!path.isEmpty())
            m_hfTable->setText(path);
    });
    connect(defaultTable, &QPushButton::clicked, m_hfTable, &QLineEdit::clear);
    connect(m_hfTable, &QLineEdit::textChanged, defaultTable, [defaultTable](const QString &path) {
        defaultTable->setEnabled(!path.trimmed().isEmpty());
    });
    connect(load, &QPushButton::clicked, this, [this] {
        auto path = QFileDialog::getOpenFileName(this, tr("载入共用 ROI 方案"),
                                                 WorkspacePaths::directory(Workflow::Roi, Kind::Config),
                                                 "XRI ROI (*.json)");
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
        if (QFileInfo(dialog.currentSchemePath()).isFile())
            loadRois(dialog.currentSchemePath());
        else
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
    m_repair = new QCheckBox(tr("修复 I0 孤立低点及对应帧"));
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
    m_filter = new QCheckBox(tr("使用像素掩膜（优先采用谱分辨结果）"));
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
    auto scanGroup = new QGroupBox(tr("XRS 扫描与元素分组"));
    auto sl = new QVBoxLayout(scanGroup);
    auto sample = new QFormLayout;
    m_sampleName = new QLineEdit;
    m_sampleName->setObjectName("SpectroscopySampleName");
    m_sampleName->setPlaceholderText(tr("留空时使用弹性扫描名"));
    sample->addRow(tr("样品名称"), m_sampleName);
    sl->addLayout(sample);

    m_edgeTable = new QTableWidget(0, 3);
    m_edgeTable->setObjectName("SpectroscopyEdgeGroups");
    configureTable(m_edgeTable, {tr("元素边"), tr("已勾选 / 总数"), tr("状态")});
    m_edgeTable->setSelectionMode(QAbstractItemView::SingleSelection);
    m_edgeTable->setEditTriggers(QAbstractItemView::NoEditTriggers);
    m_edgeTable->horizontalHeader()->setStretchLastSection(false);
    m_edgeTable->horizontalHeader()->setSectionResizeMode(0, QHeaderView::Stretch);
    m_edgeTable->horizontalHeader()->setSectionResizeMode(1, QHeaderView::ResizeToContents);
    m_edgeTable->horizontalHeader()->setSectionResizeMode(2, QHeaderView::ResizeToContents);
    m_edgeTable->setMinimumHeight(110);
    m_edgeTable->setMaximumHeight(165);
    sl->addWidget(m_edgeTable);
    auto groupActions = new QHBoxLayout;
    m_newEdge = WorkflowUi::button("plus", tr("新建分组…"));
    m_newEdge->setObjectName("AddSpectroscopyEdgeGroup");
    m_renameEdge = new QPushButton(tr("重命名…"));
    m_removeEdge = new QPushButton(tr("删除分组…"));
    for (auto action : {m_newEdge, m_renameEdge, m_removeEdge}) {
        action->setAutoDefault(false);
        groupActions->addWidget(action);
    }
    sl->addLayout(groupActions);
    connect(m_sampleName, &QLineEdit::textChanged, this, &SpectroscopyDialog::refreshExportIdentity);
    connect(m_edgeTable, &QTableWidget::currentCellChanged, this, [this](int row, int, int, int) {
        switchEdgeGroup(row);
    });
    connect(m_newEdge, &QPushButton::clicked, this, &SpectroscopyDialog::chooseNewEdgeGroup);
    connect(m_renameEdge, &QPushButton::clicked, this, [this] {
        if (isBusy() || m_activeEdge < 0)
            return;
        bool accepted;
        auto &group = m_edgeGroups[m_activeEdge];
        const auto name = QInputDialog::getText(this, tr("重命名元素边分组"), tr("元素边名称，例如 O_K 或 Ni_L3"),
                                                QLineEdit::Normal, group.name, &accepted).trimmed();
        if (!accepted)
            return;
        if (!validEdgeName(name, m_activeEdge)) {
            message(tr("请使用非空且未重复的分组名称。"), true);
            return;
        }
        group.name = name;
        refreshExportIdentity();
        refreshEdgeGroups();
    });
    connect(m_removeEdge, &QPushButton::clicked, this, [this] {
        if (isBusy() || m_activeEdge < 0)
            return;
        const auto name = m_edgeGroups[m_activeEdge].name;
        if (QMessageBox::question(this, tr("删除元素边分组"),
                                 tr("删除方案中的分组「%1」及其扫描配置？扫描和已保存的结果文件会保留。").arg(name),
                                 QMessageBox::Yes | QMessageBox::No, QMessageBox::No) != QMessageBox::Yes)
            return;
        m_edgeGroups.removeAt(m_activeEdge);
        m_activeEdge = -1;
        applyEdgeConfiguration(m_edgeDraft);
    });

    auto hfShellForm = new QFormLayout;
    auto hfShellRow = new QWidget;
    auto hfShellLayout = new QHBoxLayout(hfShellRow);
    hfShellLayout->setContentsMargins(0, 0, 0, 0);
    m_hfShells = new QLineEdit;
    m_hfShells->setObjectName("BackgroundHfShells");
    m_hfShells->setPlaceholderText(tr("例如 O:K 或 C:1"));
    m_hfShells->setToolTip(tr("新建分组时自动关联所选壳层。也可使用 HF 表中的编号，例如 C:1。"));
    auto chooseShells = new QPushButton(tr("选壳层…"));
    chooseShells->setObjectName("ChooseSpectrumHfShells");
    chooseShells->setAutoDefault(false);
    hfShellLayout->addWidget(m_hfShells, 1);
    hfShellLayout->addWidget(chooseShells);
    hfShellForm->addRow(tr("当前分组的扣除壳层"), hfShellRow);
    sl->addLayout(hfShellForm);
    connect(chooseShells, &QPushButton::clicked, this, &SpectroscopyDialog::chooseEdgeHfShells);

    sl->addWidget(new QLabel(tr("当前分组的扫描")));
    m_files = new QListWidget;
    m_files->setObjectName("XrsFiles");
    m_files->setSelectionMode(QAbstractItemView::ExtendedSelection);
    sl->addWidget(m_files, 1);
    auto sb = new QHBoxLayout;
    m_addScans = WorkflowUi::button("plus", tr("添加扫描…"));
    m_removeScans = new QPushButton(tr("移除选中"));
    sb->addWidget(m_addScans);
    sb->addWidget(m_removeScans);
    sl->addLayout(sb);
    auto scanNote = new QLabel(tr("选择一个元素边分组，再勾选参与该组计算的扫描。"));
    scanNote->setWordWrap(true);
    sl->addWidget(scanNote);
    inputLayout->addWidget(scanGroup, 2);
    connect(m_addScans, &QPushButton::clicked, this, [this] {
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
    connect(m_removeScans, &QPushButton::clicked, this, [this] {
        if (m_files->selectedItems().isEmpty())
            return;
        for (auto item : m_files->selectedItems())
            delete m_files->takeItem(m_files->row(item));
        invalidate(1);
    });
    connect(m_files, &QListWidget::itemChanged, this, [this] { invalidate(1); });
    m_tabs->addTab(inputs, tr("1 数据与 ROI"));

    auto masks = new QWidget;
    auto ml = new QVBoxLayout(masks);
    WorkflowUi::compact(ml);
    auto maskTools = new QHBoxLayout;
    maskTools->addWidget(new QLabel(tr("预览")));
    m_imageSource = new QComboBox;
    m_imageSource->setObjectName("MaskPreviewSource");
    m_imageSource->addItems({tr("弹性叠加图"), tr("最终像素掩膜")});
    maskTools->addWidget(m_imageSource);
    auto log = new QCheckBox(tr("signed log1p 显示"));
    log->setChecked(true);
    maskTools->addWidget(log);
    maskTools->addWidget(WorkflowUi::disclosureButton(tr("像素掩膜参数"), mask));
    auto saveRoi = WorkflowUi::button("save", tr("保存调整后 ROI…"));
    auto *unmix = WorkflowUi::button("chart-line", tr("像素谱分解…"));
    unmix->setObjectName("OpenSpectralMask");
    maskTools->addWidget(unmix);
    connect(unmix, &QPushButton::clicked, this, &SpectroscopyDialog::editSpectralMask);
    auto *reset = WorkflowUi::button("rotate-ccw", tr("清除此晶体谱掩膜"));
    maskTools->addWidget(reset);
    connect(reset, &QPushButton::clicked, this, [this] {
        const int r = m_corrections->currentRow();
        if (m_spectralMasks.remove(r)) {
            m_corrections->item(r, 0)->setText(m_rois.items()[r].label);
            updateMasks();
        } else
            message(tr("请在表中选择已经应用谱掩膜的晶体。"));
    });
    maskTools->addStretch();
    maskTools->addWidget(saveRoi);
    ml->addLayout(maskTools);
    ml->addWidget(mask);
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
        m_images[d]->setObjectName(d == 0 ? "LambdaMaskPreview" : "MiniPixMaskPreview");
        m_images[d]->setLogarithmic(true);
        detectorTabs->addTab(m_images[d], d == 0 ? "LAMBDA" : "MINIPIX");
    }
    maskSplit->addWidget(detectorTabs);
    maskSplit->setStretchFactor(1, 1);
    maskSplit->setSizes({340, 800});
    connect(m_corrections, &QTableWidget::currentCellChanged, this,
            [this, detectorTabs](int row, int, int, int) {
                if (row < 0)
                    return;
                detectorTabs->setCurrentIndex(row < 60 ? 0 : 1);
                for (int d = 0; d < 2; ++d) {
                    QVector<CrystalRoi> overlays;
                    for (const auto &roi : m_masks.adjusted.items())
                        if (roi.detector == (d == 0 ? "D_LAMBDA" : "D_MINIPIX"))
                            overlays.append(roi);
                    m_images[d]->setRoiOverlays(overlays, m_rois.items()[row].label);
                }
            });
    ml->addWidget(maskSplit, 1);
    m_tabs->addTab(masks, tr("2 掩膜与谱分辨"));
    connect(m_imageSource, &QComboBox::currentIndexChanged, this, [this] { refreshImages(); });
    connect(log, &QCheckBox::toggled, this, [this](bool on) {
        for (auto image : m_images)
            image->setLogarithmic(on);
    });
    connect(saveRoi, &QPushButton::clicked, this, [this] {
        if (m_masks.images.isEmpty()) {
            message(tr("请先计算掩膜。"), true);
            return;
        }
        const auto name = WorkspacePaths::saveFile(this, Workflow::Spectroscopy, Kind::Result, m_elasticFile->text(),
                                                   tr("保存调整后 ROI"), "adjusted.roi.json", "XRI ROI (*.json)");
        if (name.isEmpty())
            return;
        QSaveFile file(name);
        auto bytes = m_masks.adjusted.serialize(m_elasticFile->text(), m_masks.sizes);
        if (!file.open(QIODevice::WriteOnly) || file.write(bytes) != bytes.size() || !file.commit())
            message(tr("ROI 保存失败。"), true);
        else
            message(tr("已保存调整后 ROI：%1").arg(name));
    });

    auto elastic = new QWidget;
    auto el = new QVBoxLayout(elastic);
    WorkflowUi::compact(el);
    auto elasticTools = new QHBoxLayout;
    m_elasticView = new QComboBox;
    m_elasticView->addItems({tr("所选晶体：原 ROI / 掩膜 / 拟合"), tr("全部合格晶体：弹性零点对齐")});
    elasticTools->addWidget(m_elasticView);
    elasticTools->addWidget(WorkflowUi::disclosureButton(tr("峰形拟合与质量阈值"), fitting));
    elasticTools->addStretch();
    el->addLayout(elasticTools);
    el->addWidget(fitting);
    auto es = new QSplitter(Qt::Horizontal);
    es->setObjectName("ElasticResultsSplitter");
    es->setChildrenCollapsible(false);
    m_elasticPlot = new SpectrumPlot(nullptr, SpectrumPlot::Style::Spectroscopy);
    m_elasticPlot->setObjectName("ElasticSpectrum");
    es->addWidget(m_elasticPlot);
    m_fits = new QTableWidget(75, 7);
    m_fits->setObjectName("ElasticFits");
    configureTable(m_fits, {tr("晶体"), tr("像素 / 掩膜"), tr("中心 / eV"), tr("FWHM / eV"), "R²", tr("拟合"),
                            tr("原因")});
    es->addWidget(m_fits);
    es->setStretchFactor(0, 3);
    es->setStretchFactor(1, 2);
    es->setSizes({650, 450});
    el->addWidget(es, 1);
    m_tabs->addTab(elastic, tr("3 弹性标定"));
    connect(m_fits, &QTableWidget::itemSelectionChanged, this, &SpectroscopyDialog::refreshElasticCurve);
    connect(m_elasticView, &QComboBox::currentIndexChanged, this, [this] { refreshElasticCurve(); });

    auto xrs = new QWidget;
    auto xl = new QVBoxLayout(xrs);
    WorkflowUi::compact(xl);
    auto choices = new QHBoxLayout;
    m_scanChoice = new QComboBox;
    m_scanChoice->setObjectName("ExtractedXrsScan");
    m_scanChoice->setSizeAdjustPolicy(QComboBox::AdjustToMinimumContentsLengthWithIcon);
    m_scanChoice->setMinimumContentsLength(28);
    m_scanChoice->setMaximumWidth(400);
    m_crystalChoice = new QComboBox;
    for (const auto &roi : m_rois.items())
        m_crystalChoice->addItem(roi.label);
    choices->addWidget(new QLabel(tr("扫描")));
    choices->addWidget(m_scanChoice, 1);
    choices->addWidget(new QLabel(tr("晶体")));
    choices->addWidget(m_crystalChoice);
    auto *repairHere = new QCheckBox(m_repair->text());
    repairHere->setObjectName("RepairXrsGlitches");
    repairHere->setChecked(m_repair->isChecked());
    choices->addWidget(repairHere);
    choices->addWidget(new QLabel(tr("最小跌幅")));
    m_glitchDrop = realBox(100 * SpectrumSettings{}.glitchDropFraction, 0.1, 99.9, 1);
    m_glitchDrop->setObjectName("GlitchDropPercent");
    m_glitchDrop->setSuffix(" %");
    choices->addWidget(m_glitchDrop);
    connect(repairHere, &QCheckBox::toggled, m_repair, &QCheckBox::setChecked);
    connect(m_repair, &QCheckBox::toggled, repairHere, &QCheckBox::setChecked);
    auto updateRepairControls = [this, repairHere] {
        const bool normalize = m_normalizeXrs->isChecked();
        m_repair->setEnabled(normalize);
        repairHere->setEnabled(normalize);
        m_glitchDrop->setEnabled(normalize && m_repair->isChecked());
    };
    connect(m_normalizeXrs, &QCheckBox::toggled, this, updateRepairControls);
    connect(m_repair, &QCheckBox::toggled, this, updateRepairControls);
    updateRepairControls();
    choices->addStretch(1);
    xl->addLayout(choices);
    auto xs = new QSplitter(Qt::Horizontal);
    xs->setChildrenCollapsible(false);
    m_i0Plot = new SpectrumPlot(nullptr, SpectrumPlot::Style::Spectroscopy);
    m_i0Plot->setObjectName("I0Spectrum");
    m_scanPlot = new SpectrumPlot(nullptr, SpectrumPlot::Style::Spectroscopy);
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
    optionScroll->setObjectName("SpectrumResultParameters");
    optionScroll->setMinimumWidth(300);
    optionScroll->setMaximumWidth(390);
    resLayout->addWidget(optionScroll);
    auto *toggle = new QToolButton;
    toggle->setObjectName("ToggleSpectrumParameters");
    toggle->setCheckable(true);
    toggle->setAutoRaise(true);
    toggle->setIcon(WorkflowUi::icon("chevron-left"));
    toggle->setIconSize({12, 12});
    toggle->setFixedSize(18, 18);
    toggle->setToolTip(tr("收起参数面板"));
    toggle->setAccessibleName(toggle->toolTip());
    resLayout->addWidget(toggle, 0, Qt::AlignTop);
    connect(toggle, &QToolButton::toggled, this, [this, optionScroll, toggle](bool collapsed) {
        optionScroll->setVisible(!collapsed);
        toggle->setIcon(WorkflowUi::icon(collapsed ? "chevron-right" : "chevron-left"));
        toggle->setToolTip(collapsed ? tr("展开参数面板") : tr("收起参数面板"));
        toggle->setAccessibleName(toggle->toolTip());
    });
    auto rs = new QSplitter(Qt::Horizontal);
    rs->setObjectName("SpectrumResultsSplitter");
    rs->setChildrenCollapsible(false);
    resLayout->addWidget(rs, 1);
    auto plotPane = new QWidget;
    auto rb = new QVBoxLayout(plotPane);
    rb->setContentsMargins(0, 0, 0, 0);
    m_resultChoice = new QComboBox;
    m_resultChoice->setMinimumContentsLength(16);
    rb->addWidget(m_resultChoice, 0, Qt::AlignLeft);
    m_resultPlot = new SpectrumPlot(nullptr, SpectrumPlot::Style::Spectroscopy);
    m_resultPlot->setObjectName("CombinedSpectrum");
    rb->addWidget(m_resultPlot, 1);
    rs->addWidget(plotPane);
    auto inclusionBody = new QWidget;
    auto il = new QVBoxLayout(inclusionBody);
    auto batch = new QHBoxLayout;
    il->setContentsMargins(0, 0, 0, 0);
    auto selectionHint = new QLabel(tr("勾选允许叠加的晶体（Ctrl / Shift 可多选）"));
    selectionHint->setWordWrap(true);
    il->addWidget(selectionHint);
    auto include = new QPushButton(tr("纳入选中"));
    auto exclude = new QPushButton(tr("排除选中"));
    batch->addWidget(include);
    batch->addWidget(exclude);
    batch->addStretch();
    il->addLayout(batch);
    m_inclusion = new QTableWidget(75, 7);
    m_inclusion->setObjectName("SpectrumInclusion");
    configureTable(m_inclusion, {tr("晶体"), tr("Q 均值 / Å⁻¹"), tr("dQ 均值 / Å⁻¹"), tr("Q 变化"),
                                 tr("dQ 变化"), tr("叠加"), tr("原因")});
    il->addWidget(m_inclusion);
    rs->addWidget(inclusionBody);
    rs->setStretchFactor(0, 3);
    rs->setStretchFactor(1, 2);
    rs->setSizes({570, 335});
    m_tabs->addTab(resultPage, tr("5 原始谱结果"));
    m_processing = new SpectrumProcessingPage(*m_hfTable, *m_hfComposition, *m_hfShells);
    m_processing->setWideCalibration([this]() -> std::optional<WideScanCalibration> {
        if (!m_haveCalibration)
            return std::nullopt;
        WideScanCalibration context{m_haveResult ? m_resultSettings : settings(), m_calibration, {}, m_result.angles};
        for (qsizetype i = 0; i < m_result.included.size(); ++i)
            if (m_result.included[i])
                context.labels.append(m_calibration.adjusted.items()[i].label);
        return context;
    });
    m_tabs->addTab(m_processing, tr("6 Compton 峰分析"));
    m_tabs->addTab(m_processing->analysisPage(), tr("7 背景分析"));
    m_tabs->addTab(m_processing->resultPage(), tr("8 近边结构结果分析"));
    auto batchSet = [this](Qt::CheckState state) {
        bool changed = false;
        for (const auto &index : m_inclusion->selectionModel()->selectedRows()) {
            auto item = m_inclusion->item(index.row(), 0);
            if (item && (item->flags() & Qt::ItemIsUserCheckable)) {
                QSignalBlocker b(m_inclusion);
                item->setCheckState(state);
                if (state == Qt::Checked)
                    m_configExclusions.remove(item->row());
                else
                    m_configExclusions.insert(item->row());
                changed = true;
            }
        }
        if (changed)
            invalidate(2);
    };
    connect(include, &QPushButton::clicked, this, [batchSet] { batchSet(Qt::Checked); });
    connect(exclude, &QPushButton::clicked, this, [batchSet] { batchSet(Qt::Unchecked); });
    connect(m_inclusion, &QTableWidget::itemChanged, this, [this](QTableWidgetItem *item) {
        if (item->column() == 0 && (item->flags() & Qt::ItemIsUserCheckable)) {
            if (item->checkState() == Qt::Checked)
                m_configExclusions.remove(item->row());
            else
                m_configExclusions.insert(item->row());
            invalidate(2);
        }
    });
    connect(m_resultChoice, &QComboBox::currentIndexChanged, this, [this] { refreshResultCurve(); });
    m_status = new QLabel(tr("选择弹性扫描并载入 ROI，点击下一步计算掩膜。"));
    m_status->setObjectName("SpectrumStatus");
    m_status->setWordWrap(true);
    m_status->setTextInteractionFlags(Qt::TextSelectableByMouse);
    m_progress = new QProgressBar;
    m_progress->setMaximumWidth(140);
    m_progress->hide();
    auto actions = new QHBoxLayout;
    m_extractButton = WorkflowUi::button("chart-line", tr("提取 XRS"));
    m_extractButton->setObjectName("ExtractSpectrum");
    m_combineButton = WorkflowUi::button("layers", tr("计算 Q 并叠加"));
    m_combineButton->setObjectName("CombineSpectrum");
    m_saveButton = WorkflowUi::button("save", tr("保存谱学结果…"));
    m_saveButton->setObjectName("SaveSpectrum");
    m_cancelButton = WorkflowUi::button("square", tr("取消计算"));
    m_previousButton = WorkflowUi::button("chevron-left", tr("上一步"));
    m_previousButton->setObjectName("PreviousSpectrumStep");
    m_nextButton = WorkflowUi::button("chevron-right", tr("下一步"));
    m_nextButton->setObjectName("NextSpectrumStep");
    m_fitElasticButton = WorkflowUi::button("crosshair", tr("计算弹性峰"));
    m_fitElasticButton->setObjectName("FitElasticPeaks");
    actions->addWidget(m_previousButton);
    actions->addWidget(m_status, 1);
    actions->addWidget(m_progress);
    actions->addWidget(m_extractButton);
    actions->addWidget(m_combineButton);
    actions->addWidget(m_saveButton);
    actions->addWidget(m_fitElasticButton);
    actions->addWidget(m_nextButton);
    actions->addWidget(m_cancelButton);
    auto close = new QPushButton(tr("关闭"));
    close->setObjectName("CloseWorkflow");
    actions->addWidget(close);
    root->addLayout(actions);
    connect(m_tabs, &QTabWidget::currentChanged, this, [this](int step) {
        if (m_processing->isBusy())
            m_nextPage = -1;
        if (step == 5 && !m_updating && m_activeEdge >= 0 && m_haveResult)
            m_processing->restoreWideIfNeeded();
        refreshActions();
    });
    connect(m_processing, &SpectrumProcessingPage::busyChanged, this, &SpectroscopyDialog::refreshState);
    connect(m_processing, &SpectrumProcessingPage::analysisInputReady, this, [this] {
        if (m_tabs->currentIndex() == 5)
            selectStep(6);
    });
    connect(m_processing, &SpectrumProcessingPage::analysisFinished, this, [this](bool usable) {
        if (m_nextPage == 7 && usable)
            selectStep(7);
        m_nextPage = -1;
    });
    connect(m_processing, &SpectrumProcessingPage::analysisParametersEdited, this, [this] {
        if (!m_updating && m_activeEdge >= 0 && m_edgeGroups[m_activeEdge].status != EdgeStatus::Pending)
            setEdgeStatus(EdgeStatus::Combined);
    });
    connect(m_processing, &SpectrumProcessingPage::analysisBatchCompleted, this, [this](bool usable, bool partial) {
        setEdgeStatus(!usable ? EdgeStatus::Combined : partial ? EdgeStatus::Partial : EdgeStatus::Completed);
    });
    connect(m_processing, &SpectrumProcessingPage::resultSaved, this, [this](const QString &format, const QString &path) {
        if (m_activeEdge >= 0)
            m_edgeGroups[m_activeEdge].lastSaved[format] = QFileInfo(path).absoluteFilePath();
    });
    connect(m_nextButton, &QPushButton::clicked, this, [this] {
        advanceToPage(m_tabs->currentIndex() + 1);
    });
    connect(m_previousButton, &QPushButton::clicked, this,
            [this] { selectStep(m_tabs->currentIndex() - 1); });
    connect(close, &QPushButton::clicked, this, &QWidget::hide);
    connect(m_extractButton, &QPushButton::clicked, this, [this] { runStep(1); });
    connect(m_fitElasticButton, &QPushButton::clicked, this, [this] { runStep(0); });
    connect(m_combineButton, &QPushButton::clicked, this, [this] { runStep(2); });
    connect(m_saveButton, &QPushButton::clicked, this, &SpectroscopyDialog::save);
    connect(m_cancelButton, &QPushButton::clicked, this, [this] {
        m_cancel.store(true);
        message(tr("正在取消，等待当前帧读取结束…"));
    });
    // Invalidate only the affected stage. No stale calibration or result may be saved.
    connect(m_elasticFile, &QLineEdit::textChanged, this, [this] {
        if (m_updating)
            return;
        m_spectralMasks.clear();
        m_roiSizes.clear();
        m_configExclusions.clear();
        for (int r = 0; r < 75; ++r) {
            QSignalBlocker block(m_corrections);
            m_corrections->item(r, 0)->setText(m_rois.items()[r].label);
        }
        m_rois = CrystalRois();
        m_roiPath->clear();
        m_projection.clear();
        refreshExportIdentity();
        invalidate(prepareMaskStep);
    });
    auto projectionChanged = [this] {
        m_projection.clear();
        invalidate(prepareMaskStep);
    };
    for (auto c : {m_elasticEnergy, m_elasticI0})
        connect(c, &QComboBox::currentTextChanged, this, projectionChanged);
    connect(m_normalizeElastic, &QCheckBox::toggled, this, projectionChanged);
    connect(m_model, &QComboBox::currentTextChanged, this, [this] { invalidate(0); });
    for (auto c : {m_xrsEnergy, m_xrsI0})
        connect(c, &QComboBox::currentTextChanged, this, [this] { invalidate(1); });
    for (auto c : {m_filter, m_auto, m_overlap})
        connect(c, &QCheckBox::toggled, this, &SpectroscopyDialog::updateMasks);
    for (auto c : {m_normalizeXrs, m_repair})
        connect(c, &QCheckBox::toggled, this, [this] { invalidate(1); });
    connect(m_glitchDrop, &QDoubleSpinBox::valueChanged, this, [this] { invalidate(1); });
    for (auto b : {m_centerMin, m_centerMax, m_fwhm, m_r2})
        connect(b, &QDoubleSpinBox::valueChanged, this, [this] { invalidate(0); });
    connect(m_fraction, &QDoubleSpinBox::valueChanged, this, &SpectroscopyDialog::updateMasks);
    connect(m_roiSize, &QSpinBox::valueChanged, this, &SpectroscopyDialog::updateMasks);
    connect(m_corrections, &QTableWidget::itemChanged, this, &SpectroscopyDialog::updateMasks);
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
    s.spectralMasks = m_spectralMasks;
    s.excludedCrystals = m_configExclusions;
    s.elasticFile = m_elasticFile->text().trimmed();
    s.elasticEnergy = m_elasticEnergy->currentText().trimmed();
    s.xrsEnergy = m_xrsEnergy->currentText().trimmed();
    s.elasticI0 = m_elasticI0->currentText().trimmed();
    s.xrsI0 = m_xrsI0->currentText().trimmed();
    s.normalizeElastic = m_normalizeElastic->isChecked();
    s.normalizeXrs = m_normalizeXrs->isChecked();
    s.repairGlitches = m_repair->isChecked();
    s.glitchDropFraction = m_glitchDrop->value() / 100.;
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

QMap<QString, QWidget *> SpectroscopyDialog::configurationFields() const {
    auto fields = m_processing->configurationFields();
    fields.insert({{"elastic_energy", m_elasticEnergy},
                   {"xrs_energy", m_xrsEnergy},
                   {"elastic_i0", m_elasticI0},
                   {"xrs_i0", m_xrsI0},
                   {"normalize_elastic", m_normalizeElastic},
                   {"normalize_xrs", m_normalizeXrs},
                   {"repair_glitches", m_repair},
                   {"glitch_drop_percent", m_glitchDrop},
                   {"use_filter", m_filter},
                   {"auto_roi", m_auto},
                   {"remove_overlap", m_overlap},
                   {"filter_fraction", m_fraction},
                   {"roi_size", m_roiSize},
                   {"fit_model", m_model},
                   {"center_min_kev", m_centerMin},
                   {"center_max_kev", m_centerMax},
                   {"max_fwhm_ev", m_fwhm},
                   {"min_r_squared", m_r2},
                   {"angle_source", m_angleSource},
                   {"crystal_diameter_mm", m_diameter},
                   {"crystal_distance_mm", m_distance},
                   {"module_pitch_degrees", m_pitch},
                   {"q_min", m_qMin},
                   {"q_max", m_qMax},
                   {"energy_step_ev", m_step}});
    for (auto it = m_angles.cbegin(); it != m_angles.cend(); ++it)
        fields.insert("angle_" + it.key(), it.value());
    for (auto module : m_modules)
        fields.insert("include_" + module->text(), module);
    return fields;
}

SpectroscopyDialog::EdgeConfiguration SpectroscopyDialog::captureEdgeConfiguration() const {
    const auto spectrum = settings();
    const auto base = QDir::current();
    EdgeConfiguration out;
    out.elasticFile = absolutePath(base, spectrum.elasticFile);
    out.roiFile = absolutePath(base, m_roiPath->text().trimmed());
    out.rois = m_rois;
    out.sizes = m_roiSizes;
    out.masks = m_spectralMasks;
    out.exclusions = spectrum.excludedCrystals;
    out.wideFile = m_processing->wideFile();
    out.batchTargets = m_processing->batchSelection();
    const auto fields = configurationFields();
    for (auto it = fields.cbegin(); it != fields.cend(); ++it)
        out.parameters[it.key()] = QJsonValue::fromVariant(configProperty(it.value()).read(it.value()));
    out.parameters["hf_table"] = absolutePath(base, out.parameters.value("hf_table").toString().trimmed());
    for (int i = 0; i < m_files->count(); ++i) {
        const auto item = m_files->item(i);
        out.scans.append({item->data(Qt::UserRole).toString(), item->checkState() == Qt::Checked});
    }
    for (int i = 0; i < spectrum.shifts.size(); ++i) {
        const auto shift = spectrum.shifts[i], expansion = spectrum.expansions[i];
        out.corrections.append(QJsonArray{shift.x(), shift.y(), expansion.x(), expansion.y()});
    }
    return out;
}

QJsonObject SpectroscopyDialog::writeEdgeConfiguration(const EdgeConfiguration &edge) const {
    QJsonArray scans, exclusions;
    for (const auto &[file, enabled] : edge.scans)
        scans.append(QJsonObject{{"file", file}, {"enabled", enabled}});
    for (int row = 0; row < m_corrections->rowCount(); ++row)
        if (edge.exclusions.contains(row))
            exclusions.append(row);
    QJsonObject out{{"elastic_file", edge.elasticFile}, {"roi_file", edge.roiFile},
                    {"wide_file", edge.wideFile}, {"background_batch_targets", edge.batchTargets},
                    {"scans", scans}, {"parameters", edge.parameters},
                    {"roi_corrections", edge.corrections}, {"excluded_crystals", exclusions}};
    if (!edge.sizes.isEmpty()) {
        out["roi_snapshot"] = QJsonDocument::fromJson(edge.rois.serialize(edge.elasticFile, edge.sizes)).object();
        QJsonObject masks;
        for (auto it = edge.masks.cbegin(); it != edge.masks.cend(); ++it)
            masks[QString::number(it.key())] = it->provenance;
        out["spectral_selections"] = masks;
    }
    return out;
}

SpectroscopyDialog::EdgeConfiguration SpectroscopyDialog::readEdgeConfiguration(const QJsonObject &config,
                                                                                const QDir &base) const {
    const auto invalid = [this](const QString &key) {
        ImageUtil::fail(tr("谱学配置字段无效：%1").arg(key));
    };
    EdgeConfiguration out;
    for (const auto &key : {"elastic_file", "roi_file"})
        if (!config.value(key).isString())
            invalid(key);
    out.elasticFile = absolutePath(base, config.value("elastic_file").toString());
    out.roiFile = absolutePath(base, config.value("roi_file").toString());
    if (config.contains("wide_file") && !config.value("wide_file").isString())
        invalid("wide_file");
    out.wideFile = absolutePath(base, config.value("wide_file").toString());
    out.batchTargets = config.value("background_batch_targets");
    if (out.batchTargets.isArray()) {
        for (const auto &target : out.batchTargets.toArray())
            if (!target.isString() || target.toString().isEmpty())
                invalid("background_batch_targets");
    } else if (!(out.batchTargets.isNull() || out.batchTargets.isUndefined())) {
        invalid("background_batch_targets");
    }
    if (out.batchTargets.isUndefined())
        out.batchTargets = QJsonValue::Null;
    const auto fields = configurationFields();
    out.parameters = config.value("parameters").toObject();
    if (!out.parameters.contains("normalize_near_edge_by_scan_count"))
        out.parameters["normalize_near_edge_by_scan_count"] = false;
    const auto defaults = SpectrumProcessingPage::backgroundConfigurationDefaults();
    for (auto it = defaults.constBegin(); it != defaults.constEnd(); ++it)
        if (!out.parameters.contains(it.key()))
            out.parameters.insert(it.key(), it.value());
    for (auto it = fields.cbegin(); it != fields.cend(); ++it)
        if (!validConfigValue(it.value(), out.parameters.value(it.key())))
            invalid(it.key());
    out.parameters["hf_table"] = absolutePath(base, out.parameters.value("hf_table").toString());

    if (!config.value("scans").isArray())
        invalid("scans");
    for (const auto &entry : config.value("scans").toArray()) {
        const auto scan = entry.toObject();
        const auto scanPath = scan.value("file").toString();
        if (scanPath.isEmpty() || !scan.value("enabled").isBool())
            invalid("scans");
        out.scans.append({absolutePath(base, scanPath), scan.value("enabled").toBool()});
    }
    out.corrections = config.value("roi_corrections").toArray();
    if (out.corrections.size() != m_corrections->rowCount())
        invalid("roi_corrections");
    for (const auto &entry : out.corrections) {
        const auto row = entry.toArray();
        if (row.size() != 4)
            invalid("roi_corrections");
        for (int col = 0; col < row.size(); ++col) {
            const int v = row[col].toInt(-2049);
            if (v < (col < 2 ? -2048 : 0) || v > 2048)
                invalid("roi_corrections");
        }
    }
    if (!config.value("excluded_crystals").isArray())
        invalid("excluded_crystals");
    for (const auto &entry : config.value("excluded_crystals").toArray()) {
        const int row = entry.toInt(-1);
        if (row < 0 || row >= m_corrections->rowCount())
            invalid("excluded_crystals");
        out.exclusions.insert(row);
    }
    if (config.contains("roi_snapshot")) {
        const auto roi = config.value("roi_snapshot").toObject();
        const auto geometry = roi.value("detector_sizes").toArray();
        if (geometry.size() != 2)
            invalid("roi_snapshot");
        for (const auto &entry : geometry) {
            const auto size = entry.toArray();
            if (size.size() != 2 || size[0].toInt() <= 0 || size[1].toInt() <= 0)
                invalid("roi_snapshot");
            out.sizes.append({size[0].toInt(), size[1].toInt()});
        }
        if (QFileInfo(out.elasticFile).isFile())
            out.sizes = detectorSizes(out.elasticFile);
        QString error;
        if (!out.rois.restore(QJsonDocument(roi).toJson(), out.sizes, &error))
            ImageUtil::fail(error);
        const auto records = config.value("spectral_selections").toObject();
        for (auto it = records.constBegin(); it != records.constEnd(); ++it) {
            bool ok;
            const int row = it.key().toInt(&ok);
            if (!ok || row < 0 || row >= out.rois.items().size() || out.rois.items()[row].mask.isEmpty() ||
                !it.value().isObject())
                invalid("spectral_selections");
            out.masks.insert(row, {out.rois.items()[row].mask, it.value().toObject()});
        }
    }
    return out;
}

void SpectroscopyDialog::applyEdgeConfiguration(const EdgeConfiguration &edge) {
    {
        const QScopedValueRollback updating(m_updating, true);
        const QSignalBlocker blockFiles(m_files), blockCorrections(m_corrections);
        m_processing->resetInputs();
        m_elasticFile->setText(edge.elasticFile);
        m_roiPath->setText(edge.roiFile);
        m_rois = edge.rois;
        m_roiSizes = edge.sizes;
        m_spectralMasks = edge.masks;
        m_configExclusions = edge.exclusions;
        m_files->clear();
        for (const auto &[file, enabled] : edge.scans)
            appendScan(m_files, file, enabled);
        const auto fields = configurationFields();
        for (auto it = fields.cbegin(); it != fields.cend(); ++it)
            configProperty(it.value()).write(it.value(), edge.parameters.value(it.key()).toVariant());
        for (int row = 0; row < edge.corrections.size(); ++row) {
            const auto label = m_rois.items()[row].label;
            m_corrections->item(row, 0)->setText(m_spectralMasks.contains(row) ? label + QStringLiteral(" ✓") : label);
            const auto values = edge.corrections[row].toArray();
            for (int col = 0; col < values.size(); ++col)
                m_corrections->item(row, col + 1)->setText(QString::number(values[col].toInt()));
        }
        m_processing->restoreEdgeInputs(edge.wideFile, edge.batchTargets);
    }
    m_projection.clear();
    m_nextPage = -1;
    selectStep(0);
    invalidate(prepareMaskStep, true);
    refreshExportIdentity();
}

QString SpectroscopyDialog::edgeStatusKey(EdgeStatus status) {
    switch (status) {
    case EdgeStatus::Pending: return "pending";
    case EdgeStatus::Combined: return "raw_ready";
    case EdgeStatus::Partial: return "partial";
    case EdgeStatus::Completed: return "completed";
    }
    return {};
}

bool SpectroscopyDialog::saveConfiguration(const QString &path) {
    if (isBusy() || m_edgeGroups.isEmpty())
        return false;
    try {
        if (m_activeEdge >= 0)
            m_edgeGroups[m_activeEdge].config = captureEdgeConfiguration();
        QJsonArray groups;
        for (const auto &group : std::as_const(m_edgeGroups)) {
            auto snapshot = writeEdgeConfiguration(group.config);
            snapshot["id"] = group.id;
            snapshot["name"] = group.name;
            snapshot["status"] = edgeStatusKey(group.status);
            snapshot["last_saved_results"] = group.lastSaved;
            if (!group.targetElement.isEmpty())
                snapshot["target_edge"] = QJsonObject{{"element", group.targetElement}, {"shell", group.targetShell},
                                                       {"subtraction_initialized", true}};
            groups.append(snapshot);
        }
        const QJsonObject config{{"schema", "xri-spectroscopy-config"}, {"version", 2},
                                 {"sample_name", m_sampleName->text().trimmed()},
                                 {"active_group_id", m_activeEdge >= 0 ? m_edgeGroups[m_activeEdge].id : QString{}},
                                 {"groups", groups}};
        ImageUtil::saveJson(path, config);
        message(tr("配置已保存：%1").arg(path));
        return true;
    } catch (const std::exception &e) {
        message(QString::fromUtf8(e.what()), true);
        return false;
    }
}

bool SpectroscopyDialog::importConfiguration(const QString &path) {
    if (isBusy())
        return false;
    try {
        QFile file(path);
        if (!file.open(QIODevice::ReadOnly))
            ImageUtil::fail(file.errorString());
        const auto config = QJsonDocument::fromJson(file.readAll()).object();
        const auto version = config.value("version").toInt();
        if (config.value("schema") != "xri-spectroscopy-config" || (version != 1 && version != 2))
            ImageUtil::fail(tr("请选择有效的 XRI 谱学配置 JSON。"));
        const QDir base = QFileInfo(path).absoluteDir();
        QList<EdgeGroup> groups;
        int active = -1;
        if (version == 1) {
            groups.append(EdgeGroup{QUuid::createUuid().toString(QUuid::Id128), QStringLiteral("Default"),
                                    readEdgeConfiguration(config, base)});
            active = 0;
        } else {
            const auto records = config.value("groups").toArray();
            if (records.isEmpty() || !config.value("active_group_id").isString() || !config.value("sample_name").isString())
                ImageUtil::fail(tr("谱学方案中的元素分组无效。"));
            const auto selected = config.value("active_group_id").toString();
            QSet<QString> names, ids;
            for (const auto &record : records) {
                const auto snapshot = record.toObject();
                const auto name = snapshot.value("name").toString().trimmed();
                const auto id = snapshot.value("id").toString();
                if (name.isEmpty() || id.isEmpty() || names.contains(name.toCaseFolded()) || ids.contains(id))
                    ImageUtil::fail(tr("谱学方案中的元素分组无效。"));
                names.insert(name.toCaseFolded());
                ids.insert(id);
                EdgeGroup group{id, name, readEdgeConfiguration(snapshot, base)};
                if (snapshot.contains("target_edge")) {
                    const auto target = snapshot.value("target_edge").toObject();
                    group.targetElement = target.value("element").toString().trimmed();
                    group.targetShell = target.value("shell").toString().trimmed();
                    if (group.targetElement.isEmpty() || group.targetShell.isEmpty())
                        ImageUtil::fail(tr("谱学方案中的元素分组无效。"));
                    if (!target.value("subtraction_initialized").toBool() &&
                        group.config.parameters.value("background_hf_shells").toString().trimmed().isEmpty())
                        group.config.parameters["background_hf_shells"] = group.targetElement + ':' + group.targetShell;
                }
                const auto status = snapshot.value("status").toString();
                for (const auto candidate : {EdgeStatus::Combined, EdgeStatus::Partial, EdgeStatus::Completed})
                    if (status == edgeStatusKey(candidate))
                        group.status = candidate;
                const auto outputs = snapshot.value("last_saved_results").toObject();
                for (auto it = outputs.constBegin(); it != outputs.constEnd(); ++it)
                    if (it.value().isString())
                        group.lastSaved[it.key()] = absolutePath(base, it.value().toString());
                if (id == selected)
                    active = groups.size();
                groups.append(std::move(group));
            }
            if (!selected.isEmpty() && active < 0)
                ImageUtil::fail(tr("谱学方案中的元素分组无效。"));
        }
        auto draft = groups.first().config;
        draft.scans.clear();
        draft.wideFile.clear();
        draft.batchTargets = QJsonValue::Null;
        draft.parameters["background_hf_shells"] = "";
        m_edgeDraft = std::move(draft);
        m_edgeGroups = std::move(groups);
        m_activeEdge = active;
        m_sampleName->setText(config.value("sample_name").toString());
        applyEdgeConfiguration(m_activeEdge >= 0 ? m_edgeGroups[m_activeEdge].config : m_edgeDraft);
        message(tr("已载入 %1 个元素分组：%2。").arg(m_edgeGroups.size()).arg(path));
        return true;
    } catch (const std::exception &e) {
        message(QString::fromUtf8(e.what()), true);
        return false;
    }
}

bool SpectroscopyDialog::validEdgeName(const QString &name, int except) const {
    if (name.isEmpty())
        return false;
    for (int i = 0; i < m_edgeGroups.size(); ++i)
        if (i != except && m_edgeGroups[i].name.compare(name, Qt::CaseInsensitive) == 0)
            return false;
    return true;
}

void SpectroscopyDialog::chooseNewEdgeGroup() {
    if (isBusy())
        return;
    const auto elements = edgeElements();
    if (elements.isEmpty()) {
        message(tr("无法从内置表读取元素列表。"), true);
        return;
    }

    QDialog dialog(this);
    dialog.setWindowTitle(tr("新建元素边分组"));
    auto form = new QFormLayout(&dialog);
    auto element = new QComboBox;
    element->setObjectName("NewEdgeGroupElement");
    element->setPlaceholderText(tr("选择元素"));
    element->addItems(elements);
    element->setCurrentIndex(-1);
    auto shell = new QComboBox;
    shell->setObjectName("NewEdgeGroupShell");
    shell->setPlaceholderText(tr("选择目标壳层"));
    // Names for the experimental absorption edge are independent of the HF table's numbered profiles.
    shell->addItems({"K", "L1", "L2", "L3", "L23", "M1", "M2", "M3", "M23", "M4", "M5", "M45",
                     "N1", "N2", "N3", "N23", "N4", "N5", "N45", "N6", "N7", "N67",
                     "O1", "O2", "O3", "O23", "O4", "O5", "O45", "O6", "O7", "O67",
                     "P1", "P2", "P3", "P23", "P4", "P5", "P45"});
    shell->setCurrentIndex(-1);
    auto preview = new QLabel(QStringLiteral("—"));
    preview->setObjectName("NewEdgeGroupName");
    form->addRow(tr("元素"), element);
    form->addRow(tr("目标壳层"), shell);
    form->addRow(tr("分组名称"), preview);
    auto buttons = new QDialogButtonBox(QDialogButtonBox::Ok | QDialogButtonBox::Cancel);
    auto create = buttons->button(QDialogButtonBox::Ok);
    create->setText(tr("新建分组"));
    create->setEnabled(false);
    form->addRow(buttons);
    connect(buttons, &QDialogButtonBox::accepted, &dialog, &QDialog::accept);
    connect(buttons, &QDialogButtonBox::rejected, &dialog, &QDialog::reject);
    auto updateName = [this, element, shell, preview, create] {
        const bool selected = element->currentIndex() >= 0 && shell->currentIndex() >= 0;
        const auto name = selected ? edgeGroupName(element->currentText(), shell->currentText()) : QString{};
        const bool available = selected && validEdgeName(name);
        preview->setText(!selected ? QStringLiteral("—") : available ? name : tr("%1（已存在）").arg(name));
        create->setEnabled(available);
    };
    connect(element, &QComboBox::currentIndexChanged, &dialog, updateName);
    connect(shell, &QComboBox::currentIndexChanged, &dialog, updateName);
    if (dialog.exec() == QDialog::Accepted)
        addEdgeGroup(element->currentText(), shell->currentText());
}

bool SpectroscopyDialog::addEdgeGroup(const QString &element, const QString &shell) {
    if (isBusy())
        return false;
    const auto symbol = element.trimmed(), edgeLabel = shell.trimmed();
    const auto title = edgeGroupName(symbol, edgeLabel);
    if (symbol.isEmpty() || edgeLabel.isEmpty() || !validEdgeName(title)) {
        message(tr("请使用非空且未重复的分组名称。"), true);
        return false;
    }
    try {
        auto edge = captureEdgeConfiguration();
        if (m_activeEdge >= 0)
            m_edgeGroups[m_activeEdge].config = edge;
        else
            m_edgeDraft = edge;
        edge.scans.clear();
        edge.wideFile.clear();
        edge.batchTargets = QJsonValue::Null;
        edge.parameters["background_hf_shells"] = symbol + ':' + edgeLabel;
        m_edgeGroups.append(EdgeGroup{QUuid::createUuid().toString(QUuid::Id128), title, std::move(edge),
                                      EdgeStatus::Pending, {}, symbol, edgeLabel});
        m_activeEdge = m_edgeGroups.size() - 1;
        applyEdgeConfiguration(m_edgeGroups[m_activeEdge].config);
        message(tr("已新建分组「%1」，请添加该组的 XRS 扫描。").arg(title));
        return true;
    } catch (const std::exception &e) {
        message(QString::fromUtf8(e.what()), true);
        return false;
    }
}

void SpectroscopyDialog::switchEdgeGroup(int index) {
    if (index == m_activeEdge || m_updating)
        return;
    if (isBusy() || index < -1 || index >= m_edgeGroups.size()) {
        refreshEdgeGroups();
        return;
    }
    try {
        auto edge = captureEdgeConfiguration();
        if (m_activeEdge >= 0)
            m_edgeGroups[m_activeEdge].config = std::move(edge);
        else
            m_edgeDraft = std::move(edge);
        m_activeEdge = index;
        applyEdgeConfiguration(index >= 0 ? m_edgeGroups[index].config : m_edgeDraft);
        if (index >= 0)
            message(tr("当前分组「%1」。点击下一步按该组方案计算。").arg(m_edgeGroups[index].name));
        else
            message(tr("请先选择元素边分组。"));
    } catch (const std::exception &e) {
        refreshEdgeGroups();
        message(QString::fromUtf8(e.what()), true);
    }
}

void SpectroscopyDialog::refreshEdgeGroups() {
    QSignalBlocker block(m_edgeTable);
    m_edgeTable->setRowCount(m_edgeGroups.size());
    for (int row = 0; row < m_edgeGroups.size(); ++row) {
        const auto &group = m_edgeGroups[row];
        const bool active = row == m_activeEdge;
        qsizetype checked = 0, count = 0;
        if (active) {
            count = m_files->count();
            for (int i = 0; i < m_files->count(); ++i)
                checked += m_files->item(i)->checkState() == Qt::Checked;
        } else {
            count = group.config.scans.size();
            for (const auto &scan : group.config.scans)
                checked += scan.second;
        }
        QString status;
        switch (group.status) {
        case EdgeStatus::Pending: status = tr("未处理"); break;
        case EdgeStatus::Combined: status = tr("已叠加"); break;
        case EdgeStatus::Partial: status = tr("部分完成"); break;
        case EdgeStatus::Completed: status = tr("已完成"); break;
        }
        if (active && isBusy())
            status = tr("处理中");
        const QStringList values{group.name, QString("%1 / %2").arg(checked).arg(count), status};
        for (int column = 0; column < values.size(); ++column) {
            auto item = m_edgeTable->item(row, column);
            if (!item) {
                item = cell({});
                m_edgeTable->setItem(row, column, item);
            }
            item->setText(values[column]);
            if (column == 0)
                item->setToolTip(group.name);
            if (column == 2 && group.status != EdgeStatus::Pending)
                item->setToolTip(tr("记录上次处理状态。切换分组或重新载入方案后，可按该组配置重新计算。"));
            else if (column == 2)
                item->setToolTip({});
        }
    }
    if (m_activeEdge >= 0 && m_edgeTable->currentRow() != m_activeEdge)
        m_edgeTable->setCurrentCell(m_activeEdge, 0);
    else if (m_activeEdge < 0)
        m_edgeTable->setCurrentCell(-1, -1);
}

void SpectroscopyDialog::setEdgeStatus(EdgeStatus status) {
    if (m_activeEdge < 0 || m_edgeGroups[m_activeEdge].status == status)
        return;
    m_edgeGroups[m_activeEdge].status = status;
    refreshEdgeGroups();
}

void SpectroscopyDialog::chooseEdgeHfShells() {
    if (isBusy() || m_activeEdge < 0)
        return;
    try {
        const auto bindings = backgroundHfShells(m_hfTable->text(), m_hfComposition->text().trimmed());
        auto selected = resolveHfShellSelection(bindings, backgroundHfShellSelection(m_hfShells->text()));
        selected.identifiers.sort();
        QDialog dialog(this);
        dialog.setWindowTitle(tr("选择「%1」要扣除的 HF 壳层").arg(m_edgeGroups[m_activeEdge].name));
        auto layout = new QVBoxLayout(&dialog);
        layout->addWidget(new QLabel(tr("勾选当前分组要扣除的壳层。")));
        if (!selected.missing.isEmpty()) {
            auto missing = new QLabel(tr("未在当前组分及 HF 表中找到：%1").arg(selected.missing.join(", ")));
            missing->setWordWrap(true);
            layout->addWidget(missing);
        }
        auto list = new QListWidget;
        list->setObjectName("SpectroscopyHfShellChoices");
        for (auto it = bindings.cbegin(); it != bindings.cend(); ++it) {
            const auto label = it->edge.isEmpty() ? it.key() : it.key() + QStringLiteral(" · ") + it->edge;
            auto item = new QListWidgetItem(tr("%1 · 结合能 %2 eV").arg(label, number(it->bindingEv)), list);
            item->setData(Qt::UserRole, it.key());
            item->setFlags(item->flags() | Qt::ItemIsUserCheckable);
            item->setCheckState(selected.identifiers.contains(it.key()) ? Qt::Checked : Qt::Unchecked);
        }
        layout->addWidget(list);
        auto buttons = new QDialogButtonBox(QDialogButtonBox::Ok | QDialogButtonBox::Cancel);
        connect(buttons, &QDialogButtonBox::accepted, &dialog, &QDialog::accept);
        connect(buttons, &QDialogButtonBox::rejected, &dialog, &QDialog::reject);
        layout->addWidget(buttons);
        if (dialog.exec() == QDialog::Accepted) {
            QStringList shells;
            for (int i = 0; i < list->count(); ++i)
                if (const auto item = list->item(i); item->checkState() == Qt::Checked)
                    shells.append(item->data(Qt::UserRole).toString());
            if (shells != selected.identifiers)
                m_hfShells->setText(shells.join(", "));
        }
    } catch (const std::exception &e) {
        message(QString::fromUtf8(e.what()), true);
    }
}

QString SpectroscopyDialog::sampleStem() const {
    const auto elastic = QFileInfo(m_elasticFile->text().trimmed()).completeBaseName();
    return fileToken(m_sampleName->text().trimmed(), fileToken(elastic, QStringLiteral("sample")));
}

QString SpectroscopyDialog::resultStem() const {
    if (m_activeEdge < 0)
        return {};
    const auto &group = m_edgeGroups[m_activeEdge];
    auto token = fileToken(group.name, "edge_" + fileToken(group.id, "group").left(8));
    for (const auto &other : m_edgeGroups)
        if (other.id != group.id && token.compare(fileToken(other.name, {}), Qt::CaseInsensitive) == 0) {
            token += '_' + fileToken(group.id, "group").left(8);
            break;
        }
    return sampleStem() + '_' + token;
}

QJsonObject SpectroscopyDialog::edgeIdentity() const {
    if (m_activeEdge < 0)
        return {};
    const auto &group = m_edgeGroups[m_activeEdge];
    const auto name = m_sampleName->text().trimmed();
    QJsonObject identity{{"sample_name", name.isEmpty() ? sampleStem() : name},
                         {"edge_group", group.name}, {"edge_group_id", group.id}};
    if (!group.targetElement.isEmpty()) {
        identity["target_element"] = group.targetElement;
        identity["target_shell"] = group.targetShell;
    }
    return identity;
}

void SpectroscopyDialog::refreshExportIdentity() {
    m_processing->setExportIdentity(resultStem(), edgeIdentity());
}

void SpectroscopyDialog::invalidate(int stage, bool preserveEdgeStatus) {
    if (m_updating || m_busy)
        return;
    if (!preserveEdgeStatus)
        setEdgeStatus(EdgeStatus::Pending);
    if (stage == prepareMaskStep)
        m_masks = {};
    if (stage <= 0) {
        m_haveCalibration = false;
        m_calibration = ElasticCalibration();
        m_fits->clearContents();
        m_inclusion->clearContents();
        m_elasticPlot->setCurves({}, {});
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
    m_processing->setCurrentResult({});
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
    message(stage == prepareMaskStep ? tr("点击下一步计算掩膜。")
            : stage == 0 ? tr("弹性参数或 ROI 已变更，请重新标定。")
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
    refreshEdgeGroups();
    const bool idle = !isBusy(), hasEdge = m_activeEdge >= 0;
    m_edgeTable->setEnabled(idle);
    m_hfComposition->setEnabled(idle);
    m_hfTable->parentWidget()->setEnabled(idle);
    m_hfShells->parentWidget()->setEnabled(idle && hasEdge);
    m_sampleName->setEnabled(idle);
    m_newEdge->setEnabled(idle);
    m_renameEdge->setEnabled(idle && hasEdge);
    m_removeEdge->setEnabled(idle && hasEdge);
    m_files->setEnabled(idle && hasEdge);
    m_addScans->setEnabled(idle && hasEdge);
    m_removeScans->setEnabled(idle && hasEdge && m_files->count() > 0);
    m_saveConfig->setEnabled(idle && !m_edgeGroups.isEmpty());
    m_extractButton->setEnabled(idle && hasEdge && m_haveCalibration);
    m_combineButton->setEnabled(idle && hasEdge && m_haveScans);
    m_saveButton->setEnabled(idle && hasEdge && m_haveResult);
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
        const auto sizes = detectorSizes(m_elasticFile->text());
        QFile file(path);
        QString error;
        CrystalRois rois;
        if (!file.open(QIODevice::ReadOnly))
            throw std::runtime_error(file.errorString().toUtf8().constData());
        if (!rois.restore(file.readAll(), sizes, &error))
            throw std::runtime_error(error.toUtf8().constData());
        bool changed = rois.serialize({}, sizes) != m_rois.serialize({}, sizes);
        if (changed) {
            m_spectralMasks.clear();
            m_configExclusions.clear();
            QSignalBlocker block(m_corrections);
            for (int r = 0; r < 75; ++r)
                m_corrections->item(r, 0)->setText(rois.items()[r].label);
        }
        m_rois = rois;
        m_roiSizes = sizes;
        m_roiPath->setText(QFileInfo(path).absoluteFilePath());
        if (!changed) {
            refreshState();
            return true;
        }
        invalidate(prepareMaskStep);
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
void SpectroscopyDialog::setAvailableFiles(const QStringList &files) {
    m_availableFiles = files;
    m_processing->setAvailableFiles(files);
}

void SpectroscopyDialog::addScans(const QStringList &files) {
    if (isBusy() || files.isEmpty())
        return;
    if (m_activeEdge < 0) {
        message(tr("请先选择元素边分组。"), true);
        return;
    }
    const int before = m_files->count();
    for (const auto &file : files) {
        auto path = QFileInfo(file).absoluteFilePath();
        bool duplicate = false;
        const auto canonical = QFileInfo(path).canonicalFilePath();
        for (int i = 0; i < m_files->count(); ++i) {
            const auto existing = m_files->item(i)->data(Qt::UserRole).toString();
            if (existing == path || (!canonical.isEmpty() && QFileInfo(existing).canonicalFilePath() == canonical))
                duplicate = true;
        }
        if (duplicate)
            continue;
        appendScan(m_files, path, true);
    }
    if (m_files->count() != before)
        invalidate(1);
}
void SpectroscopyDialog::runStep(int step) {
    if (isBusy() || step < prepareMaskStep || step > 2)
        return;
    if (m_activeEdge < 0) {
        message(tr("请先选择元素边分组。"), true);
        selectStep(0);
        return;
    }
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
            message(tr("请为当前元素边分组添加并勾选 XRS 扫描。"), true);
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
        auto masks = m_masks;
        auto projection = m_projection;
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
        m_watcher.setFuture(QtConcurrent::run([this, step, s, rois, calibration, masks, projection, scans, files, progress] {
            Work out;
            try {
                if (step == prepareMaskStep || step == 0) {
                    out.masks = masks.images.isEmpty()
                        ? Spectroscopy::prepareMasks(s, rois, m_cancel, progress, projection) : masks;
                    if (step == 0)
                        out.calibration = Spectroscopy::fitElastic(s, out.masks, m_cancel, progress);
                }
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

void SpectroscopyDialog::advanceToPage(int target) {
    if (isBusy())
        return;
    if (m_activeEdge < 0) {
        message(tr("请先选择元素边分组。"), true);
        selectStep(0);
        return;
    }
    m_nextPage = -1;
    if (target <= 5) {
        int job;
        if (m_masks.images.isEmpty())
            job = prepareMaskStep;
        else if (target >= 2 && !m_haveCalibration)
            job = 0;
        else if (target >= 3 && !m_haveScans)
            job = 1;
        else if (target >= 4 && !m_haveResult)
            job = 2;
        else {
            selectStep(target);
            return;
        }
        runStep(job);
        if (m_busy)
            m_nextPage = target;
        return;
    }
    if (target == 6) {
        if (m_processing->prepareAnalysisInput())
            selectStep(6);
        return;
    }
    if (m_processing->ensureAnalysis()) {
        selectStep(target);
    } else if (m_processing->isBusy()) {
        m_nextPage = target;
    }
}
void SpectroscopyDialog::updateMasks() {
    if (m_busy || m_updating)
        return;
    invalidate(prepareMaskStep);
    if (m_projection.isEmpty())
        return;
    try {
        const std::atomic_bool cancel{false};
        m_masks = Spectroscopy::prepareMasks(settings(), m_rois, cancel, {}, m_projection);
        refreshImages();
        message(tr("掩膜已计算。可修正像素，或点击下一步计算弹性标定。"));
    } catch (const std::exception &e) {
        message(QString::fromUtf8(e.what()), true);
    }
}

void SpectroscopyDialog::editSpectralMask() {
    if (m_busy)
        return;
    if (!m_rois.count(RoiState::Valid)) {
        message(tr("请先载入初始 ROI。"), true);
        return;
    }
    try {
        const auto s = settings();
        QMap<QString, QSize> sizes;
        for (const auto &roi : m_rois.items()) {
            if (roi.state != RoiState::Valid || sizes.contains(roi.detector))
                continue;
            const auto node = HdfFile::inspect(s.elasticFile, "/entry/data/" + roi.detector);
            if (node.dims.size() != 3)
                throw std::runtime_error("Expected a detector frame stack.");
            sizes.insert(roi.detector, {int(node.dims[2]), int(node.dims[1])});
        }
        QStringList files;
        for (int i = 0; i < m_files->count(); ++i)
            files.append(m_files->item(i)->data(Qt::UserRole).toString());
        SpectralMaskDialog dialog(s, m_rois.items(), sizes, m_corrections->currentRow(), files, this);
        if (dialog.exec() != QDialog::Accepted)
            return;
        m_corrections->setCurrentCell(dialog.currentRow(), 0);
        const auto &accepted = dialog.acceptedMasks();
        if (accepted.isEmpty()) {
            message(tr("已完成确认，所有晶体均保留原 ROI。"));
            return;
        }
        auto updated = m_rois;
        auto options = s;
        for (auto it = accepted.cbegin(); it != accepted.cend(); ++it) {
            const auto &roi = m_rois.items()[it.key()];
            if (!updated.assignMask(it.key(), it->pixels, sizes.value(roi.detector)))
                throw std::runtime_error(tr("%1 的组分掩膜无效。").arg(roi.label).toStdString());
            options.spectralMasks.insert(it.key(), it.value());
        }
        // Construct the complete preview before publishing any of the reviewed ROIs.
        const std::atomic_bool cancel{false};
        auto masks = Spectroscopy::prepareMasks(options, updated, cancel, {}, m_projection);
        invalidate(prepareMaskStep);
        m_rois = std::move(updated);
        m_spectralMasks = std::move(options.spectralMasks);
        m_masks = std::move(masks);
        m_projection = m_masks.images;
        {
            QSignalBlocker block(m_corrections);
            for (auto it = accepted.cbegin(); it != accepted.cend(); ++it)
                m_corrections->item(it.key(), 0)->setText(m_rois.items()[it.key()].label + QStringLiteral(" ✓"));
        }
        m_imageSource->setCurrentIndex(1);
        refreshImages();
        message(tr("已更新 %1 个晶体的原始 ROI。").arg(accepted.size()));
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
                item->setCheckState(m_configExclusions.contains(r) ? Qt::Unchecked : Qt::Checked);
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
    if (m_haveScans && m_imageSource->count() == 2)
        m_imageSource->addItem(tr("当前 XRS 扫描叠加图"));
    if (!m_haveScans && m_imageSource->count() == 3)
        m_imageSource->removeItem(2);
    int choice = m_imageSource->currentIndex();
    if (!m_masks.images.isEmpty()) {
        for (int d = 0; d < 2; ++d) {
            DataSlice preview;
            if (choice == 0)
                preview = m_masks.images[d];
            else if (choice == 1)
                preview = m_masks.maskImages[d];
            else if (m_haveScans && m_scanChoice->currentIndex() >= 0)
                preview = m_scans[m_scanChoice->currentIndex()].images[d];
            m_images[d]->setData(preview);
            QVector<CrystalRoi> overlays;
            for (const auto &roi : m_masks.adjusted.items())
                if (roi.detector == (d == 0 ? "D_LAMBDA" : "D_MINIPIX"))
                    overlays.append(roi);
            const int row = m_corrections->currentRow();
            m_images[d]->setRoiOverlays(overlays, row >= 0 ? m_rois.items()[row].label : QString());
        }
    } else
        for (auto image : m_images) {
            image->setData({});
            image->setRoiOverlays({}, {});
        }
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
        curves.append({tr("像素掩膜"), scan.energyKev, scan.filtered.value(row), QColor(0, 100, 180), false});
        if (!fit.fitted.isEmpty())
            curves.append({tr("拟合"), scan.energyKev, fit.fitted, QColor(190, 45, 35), true, SpectrumCurve::Draw::Line});
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
    QVector<SpectrumCurve> monitor{{tr("修复前 |I0|"), scan.energyKev, scan.i0Before, QColor(120, 120, 120)},
                                   {tr("修复后 |I0|"), scan.energyKev, scan.i0After, QColor(185, 45, 35), true, SpectrumCurve::Draw::Line}};
    if (!scan.glitches.isEmpty()) {
        SpectrumCurve repaired{tr("已修复点"), {}, {}, QColor(185, 45, 35), false, SpectrumCurve::Draw::Circles};
        for (int i : scan.glitches) {
            repaired.x.append(scan.energyKev[i]);
            repaired.y.append(scan.i0After[i]);
        }
        monitor.append(repaired);
    }
    m_i0Plot->setCurves(monitor, tr("I0 · 已修复 %1 点 · %2")
                            .arg(scan.glitches.size())
                            .arg(scan.i0Path.isEmpty() ? tr("未归一化，I0 = 1") : scan.i0Path),
                        "Energy (keV)", "|I0| (a.u.)");
    m_scanPlot->setCurves(
        {{tr("原 ROI"), scan.energyKev, scan.raw.value(r), QColor(100, 100, 100), false},
         {tr("像素掩膜"), scan.energyKev, scan.filtered.value(r), QColor(0, 100, 180), false}},
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
                m_resultChoice->addItem(m_rois.items()[r].label, r);
    }
    if (m_angleSource->currentIndex() != 2)
        for (auto it = m_result.angles.begin(); it != m_result.angles.end(); ++it)
            if (m_angles.contains(it.key())) {
                QSignalBlocker block(m_angles[it.key()]);
                m_angles[it.key()]->setValue(it.value());
            }
    m_updating = false;
    refreshResultCurve();
    auto spectrum =
        SpectrumProcessing::fromResult(m_result, m_calibration.adjusted, m_resultSettings.normalizeXrs);
    for (const auto &scan : m_scans)
        spectrum.sources.append(scan.file);
    m_processing->setCurrentResult(std::move(spectrum));
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
    if (!m_haveResult || isBusy() || m_activeEdge < 0)
        return;
    QString directory;
    const auto name = resultStem() + "_raw";
    try {
        directory = WorkspacePaths::newPath(Workflow::Spectroscopy, Kind::Result, {}, name);
    } catch (const std::exception &error) {
        message(QString::fromUtf8(error.what()), true);
        return;
    }
    if (directory.isEmpty()) {
        const auto parent = QFileDialog::getExistingDirectory(this, tr("选择谱学结果的父目录"),
                                                               QFileInfo(m_elasticFile->text()).absolutePath());
        if (parent.isEmpty())
            return;
        directory = QDir(parent).filePath(name + '_' + QDateTime::currentDateTime().toString("yyyyMMdd_HHmmss_zzz"));
    }
    QString error;
    if (Spectroscopy::saveResult(directory, m_resultSettings, m_calibration, m_scans, m_result, &error, edgeIdentity())) {
        m_edgeGroups[m_activeEdge].lastSaved["raw"] = QFileInfo(directory).absoluteFilePath();
        message(tr("已保存谱学 HDF5、总谱、模组谱、单晶体谱、ROI、掩膜和完整参数：%1").arg(directory));
    } else
        message(error, true);
}

void SpectroscopyDialog::refreshActions() {
    int step = m_tabs->currentIndex();
    m_tabs->cornerWidget()->setVisible(step == 0);
    m_tabs->cornerWidget()->setEnabled(!isBusy());
    m_previousButton->setVisible(step > 0);
    m_previousButton->setEnabled(!isBusy());
    m_nextButton->setVisible(step < m_tabs->count() - 1);
    m_nextButton->setEnabled(!isBusy() && m_activeEdge >= 0);
    m_extractButton->setVisible(step == 3);
    m_fitElasticButton->setVisible(step == 2);
    m_fitElasticButton->setEnabled(!isBusy() && m_activeEdge >= 0);
    m_combineButton->setVisible(step == 4);
    m_saveButton->setVisible(step == 4);
    m_cancelButton->setVisible(m_busy);
    m_nextButton->setDefault(m_nextButton->isVisibleTo(this));
}
void SpectroscopyDialog::useRoiDraft(const QString &file, const QString &path) {
    if (m_busy) {
        m_pendingRoi = {file, path};
        return;
    }
    if (QFileInfo(m_elasticFile->text()).canonicalFilePath() == file && QFileInfo::exists(path))
        loadRois(path);
}
