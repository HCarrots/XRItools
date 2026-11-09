#include "ImagingDialog.h"
#include "ImagingIO.h"
#include "ImagingInputUi.h"
#include "RoiDialog.h"
#include "ScientificPlot.h"
#include "SpectrumPlot.h"
#include "WorkflowUi.h"
#include <QtConcurrent/QtConcurrentRun>
#include <QtWidgets>
#include <algorithm>
#include <limits>
#include <numeric>
namespace {
QJsonValue fieldValue(const QJsonObject &c, const QString &path) {
    auto keys = path.split('.');
    return c[keys[0]].toObject()[keys[1]];
}
void put(QJsonObject &c, const QString &path, const QJsonValue &value) {
    auto keys = path.split('.');
    auto group = c[keys[0]].toObject();
    group[keys[1]] = value;
    c[keys[0]] = group;
}
QWidget *scrollBody(QWidget *body) {
    auto scroll = new QScrollArea;
    scroll->setWidgetResizable(true);
    scroll->setFrameShape(QFrame::NoFrame);
    scroll->setWidget(body);
    return scroll;
}
QTableWidgetItem *cell(const QString &s) {
    auto item = new QTableWidgetItem(s);
    item->setFlags(Qt::ItemIsEnabled | Qt::ItemIsSelectable);
    return item;
}
QString real(double v) {
    return QString::number(v, 'g', 8);
}
} // namespace
ImagingDialog::ImagingDialog(const QString &sample, QWidget *parent)
    : QDialog(parent), m_base(Imaging::defaults()) {
    setObjectName("ImagingWorkbench");
    setWindowTitle(tr("成像处理 — 薄膜标定与三维联合反演 — XRI"));
    setWindowFlag(Qt::Window);
    resize(1300, 880);
    setMinimumSize(1000, 690);
    buildUi();
    auto initial = Imaging::defaults();
    put(initial, "acquisition.calibration_exposure_s", .5);
    put(initial, "acquisition.sample_exposure_s", .5);
    if (!sample.isEmpty() && HdfFile::formatFromFileName(sample) != ScanFormat::Scan)
        put(initial, "input.sample_nxs", QFileInfo(sample).absoluteFilePath());
    put(initial, "runtime.output_dir",
        QDir(QFileInfo(sample).absolutePath())
            .filePath("imaging_" + QDateTime::currentDateTime().toString("yyyyMMdd_HHmmss_zzz")));
    applyConfiguration(initial);
    connect(&m_watcher, &QFutureWatcher<Work>::finished, this, [this] {
        auto work = m_watcher.result();
        m_busy = false;
        m_inputs->setEnabled(true);
        m_parameters->setEnabled(true);
        for (auto b : {m_validateButton, m_runButton, m_importButton, m_exportButton, m_includeButton,
                       m_excludeButton})
            b->setEnabled(true);
        m_cancelButton->setEnabled(false);
        m_progress->hide();
        if (!work.error.isEmpty()) {
            m_excludeButton->setEnabled(false);
            m_includeButton->setEnabled(false);
            status(work.error, true);
            return;
        }
        m_validation = std::move(work.validation);
        refreshValidation();
        if (m_validateOnly) {
            m_tabs->setCurrentIndex(2);
            status(tr("输入校验通过：薄膜与样品网格完整，能量范围匹配。输出目录须为尚不存在的新目录。"));
        } else {
            m_result = std::move(work.result);
            m_haveResult = true;
            refreshResult();
            m_tabs->setCurrentIndex(3);
            status(tr("重建已保存：%1；主结果采用 %2。请结合残差、灵敏度和收敛状态检查。")
                       .arg(m_result.directory, m_result.summary["solver"].toString()));
        }
    });
}
ImagingDialog::~ImagingDialog() {
    m_cancel.store(true);
    m_watcher.waitForFinished();
}
QString ImagingDialog::statusText() const {
    return m_status->text();
}
void ImagingDialog::status(const QString &message, bool error) {
    m_status->setText(message);
    refreshActions();
    if (!m_busy && !m_pendingRoi.second.isEmpty()) {
        const auto pending = m_pendingRoi;
        m_pendingRoi = {};
        useRoiDraft(pending.first, pending.second);
        return;
    }
    m_status->setStyleSheet(error ? "color:#a02020;" : "");
}
void ImagingDialog::selectStep(int step) {
    m_tabs->setCurrentIndex(std::clamp(step, 0, 4));
}
void ImagingDialog::buildUi() {
    auto root = new QVBoxLayout(this);
    WorkflowUi::compact(root);
    m_configTools = new QWidget;
    auto toolbar = new QHBoxLayout(m_configTools);
    toolbar->setContentsMargins(0, 0, 0, 0);
    m_importButton = new QPushButton(tr("载入已有方案…"));
    m_importButton->setObjectName("ImportImagingConfig");
    m_exportButton = new QPushButton(tr("保存当前配置…"));
    auto load = new QPushButton(tr("查看已有重建…"));
    toolbar->addWidget(m_importButton);
    toolbar->addWidget(m_exportButton);
    toolbar->addWidget(load);
    m_tabs = new QTabWidget;
    m_tabs->setObjectName("ImagingSteps");
    m_tabs->setCornerWidget(m_configTools);
    root->addWidget(m_tabs, 1);
    auto addText = [this](QFormLayout *layout, const QString &key, const QString &label, int kind = 0) {
        auto edit = new QLineEdit;
        edit->setObjectName(key);
        edit->setProperty("fieldKind", kind);
        m_fields[key] = edit;
        layout->addRow(label, edit);
        connect(edit, &QLineEdit::textChanged, this, &ImagingDialog::invalidate);
        return edit;
    };
    auto addFlag = [this](QFormLayout *layout, const QString &key, const QString &label) {
        auto box = new QCheckBox(label);
        box->setObjectName(key);
        m_fields[key] = box;
        layout->addRow(box);
        connect(box, &QCheckBox::toggled, this, &ImagingDialog::invalidate);
        return box;
    };
    auto addNumber = [this](QFormLayout *layout, const QString &key, const QString &label,
                            bool integer = false) {
        if (integer) {
            auto spin = new QSpinBox;
            spin->setRange(1, std::numeric_limits<int>::max() - 1);
            spin->setKeyboardTracking(false);
            spin->setObjectName(key);
            m_fields[key] = spin;
            layout->addRow(label, spin);
            connect(spin, &QSpinBox::valueChanged, this, &ImagingDialog::invalidate);
            return static_cast<QWidget *>(spin);
        }
        auto edit = new QLineEdit;
        edit->setObjectName(key);
        edit->setProperty("fieldKind", 1);
        m_fields[key] = edit;
        layout->addRow(label, edit);
        connect(edit, &QLineEdit::textChanged, this, &ImagingDialog::invalidate);
        return static_cast<QWidget *>(edit);
    };
    auto addNullable = [this, &addNumber](QFormLayout *layout, const QString &key, const QString &label,
                                          const QString &automatic) {
        auto editor = addNumber(layout, key, label);
        auto box = new QCheckBox(automatic);
        m_auto[key] = box;
        layout->addRow(box);
        connect(box, &QCheckBox::toggled, this, [this, editor](bool on) {
            editor->setEnabled(!on);
            invalidate();
        });
    };
    m_inputs = new QWidget;
    auto inputLayout = new QHBoxLayout(m_inputs);
    WorkflowUi::compact(inputLayout);
    auto left = new QWidget;
    auto ll = new QVBoxLayout(left);
    inputLayout->addWidget(scrollBody(left), 3);
    auto right = new QWidget;
    auto rr = new QVBoxLayout(right);
    ImagingInputUi::build(
        this, ll, rr, m_fields, m_availableFiles, addText, addNumber, [this] { invalidate(); },
        [this](const QString &file) {
            if (!QFileInfo::exists(file)) {
                status(tr("请先选择样品文件。"), true);
                return;
            }
            if (!isWindow()) {
                emit roiRequested(file);
                return;
            }
            RoiDialog dialog(file, tr("成像数据处理"), this);
            dialog.exec();
            const auto canonical = QFileInfo(file).canonicalFilePath();
            const auto draft =
                QStandardPaths::writableLocation(QStandardPaths::AppLocalDataLocation) + "/roi-drafts/" +
                QCryptographicHash::hash(canonical.toUtf8(), QCryptographicHash::Sha256).toHex() + ".json";
            if (QFileInfo::exists(draft))
                useRoiDraft(canonical, draft);
            else
                status(tr("此样品尚无共用 ROI 草稿，请先选取或载入已有方案。"), true);
        });
    ll->addStretch();
    auto fileField = [this](QFormLayout *form, const QString &key, const QString &label,
                            bool output = false) {
        auto row = new QWidget;
        auto l = new QHBoxLayout(row);
        l->setContentsMargins(0, 0, 0, 0);
        auto edit = new QLineEdit;
        edit->setObjectName(key);
        m_fields[key] = edit;
        auto choose = new QPushButton(tr("浏览…"));
        l->addWidget(edit, 1);
        l->addWidget(choose);
        form->addRow(label, row);
        connect(edit, &QLineEdit::textChanged, this, &ImagingDialog::invalidate);
        connect(choose, &QPushButton::clicked, this, [this, edit, key, output] {
            QString path;
            if (output) {
                auto parent = QFileDialog::getExistingDirectory(this, tr("选择结果的父目录"),
                                                                QFileInfo(edit->text()).absolutePath());
                if (!parent.isEmpty())
                    path = QDir(parent).filePath(
                        "imaging_" + QDateTime::currentDateTime().toString("yyyyMMdd_HHmmss_zzz"));
            } else if (key != "input.roi_table") {
                path = WorkflowUi::selectSource(this, m_availableFiles, edit->text(), tr("选择 %1").arg(key));
            } else
                path = QFileDialog::getOpenFileName(
                    this, tr("选择 %1").arg(key), QFileInfo(edit->text()).absolutePath(),
                    key == "input.roi_table" ? "ROI (*.json *.txt *.tsv);;All files (*)"
                                             : "NeXus / HDF5 (*.nxs *.h5 *.hdf5)");
            if (!path.isEmpty())
                edit->setText(path);
        });
    };
    auto output = new QGroupBox(tr("结果与诊断"));
    auto of = new QFormLayout(output);
    fileField(of, "runtime.output_dir", tr("新输出目录"), true);
    addFlag(of, "output.save_diagnostics", tr("保存各晶体正演矩阵、预测和残差"));
    addFlag(of, "output.plots", tr("保存诊断 PNG"));
    auto outputNote = new QLabel(tr("完整流程自动保存主结果、NNLS/TV、XYZ "
                                    "坐标、像素缓存、PSF、灵敏度、残差与参数。已存在的输出目录会拒绝写入。取"
                                    "消或失败时保留中间文件与日志。"));
    outputNote->setWordWrap(true);
    of->addRow(outputNote);
    rr->addWidget(output);
    auto boundary = new QLabel(
        tr("主输出：相对 Kapton 等效弹性散射密度 "
           "I(X,Y,Z)。能量已积分；样品背景未扣除。完整坐标保留，不追加几何裁剪、阈值分割、对称清理或平滑。"));
    boundary->setWordWrap(true);
    rr->addWidget(boundary);
    rr->addStretch();
    inputLayout->addWidget(scrollBody(right), 2);
    m_tabs->addTab(m_inputs, tr("1 文件与 ROI"));
    m_parameters = new QWidget;
    auto pl = new QHBoxLayout(m_parameters);
    auto psf = new QGroupBox(tr("连续移动椭圆高斯 PSF"));
    auto pf = new QFormLayout(psf);
    addNumber(pf, "psf.min_peak_fraction", tr("最小相对峰高"));
    addNumber(pf, "psf.min_r_squared", tr("最低拟合 R²"));
    addNumber(pf, "psf.edge_margin_pixels", tr("中心距边界 / pixel"));
    addNumber(pf, "psf.min_observations", tr("最少合格观测"), true);
    addNumber(pf, "psf.max_trajectory_rms_pixels", tr("轨迹最大 RMS / pixel"));
    addText(pf, "psf.sigma_bounds_pixels", tr("σ 边界 [min,max]"), 4);
    addNumber(pf, "psf.max_fit_evaluations", tr("单幅拟合调用上限"), true);
    auto psfNote =
        new QLabel(tr("薄膜位置是标定观测，不是恢复 X "
                      "层数。中间位置由单个移动光斑与线性质心轨迹生成；拟合的常数背景不进入正演核。"));
    psfNote->setWordWrap(true);
    pf->addRow(psfNote);
    auto psfBody = new QWidget;
    auto pbl = new QVBoxLayout(psfBody);
    pbl->addWidget(psf);
    pbl->addStretch();
    pl->addWidget(scrollBody(psfBody));
    auto reconstruction = new QGroupBox(tr("联合非负最小二乘与晶体增益"));
    auto nf = new QFormLayout(reconstruction);
    addNullable(nf, "grid.x_step_mm", tr("样品 X 步长 / mm"), tr("自动使用合格晶体样品长度/像素的中位数"));
    addFlag(nf, "reconstruction.fit_gains", tr("交替拟合相对晶体增益"));
    addText(nf, "reconstruction.anchor_roi", tr("固定增益 = 1 的锚点"));
    addText(nf, "reconstruction.gain_bounds", tr("增益边界 [min,max]"), 4);
    addNumber(nf, "reconstruction.max_gain_iterations", tr("增益迭代上限"), true);
    addNumber(nf, "reconstruction.gain_tolerance", tr("增益相对变化容差"));
    addNumber(nf, "reconstruction.nnls_max_iterations", tr("NNLS 迭代上限"), true);
    addNumber(nf, "reconstruction.nnls_tolerance", tr("投影 KKT 容差"));
    auto nnlsNote =
        new QLabel(tr("A = 增益 × Δx / 薄膜厚度 × 连续响应。全部晶体约束同一个非负密度。数值零模式使用相对 "
                      "1e−14 规则；不按灵敏度裁剪 X。"));
    nnlsNote->setWordWrap(true);
    nf->addRow(nnlsNote);
    auto nb = new QWidget;
    auto nbl = new QVBoxLayout(nb);
    nbl->addWidget(reconstruction);
    nbl->addStretch();
    pl->addWidget(scrollBody(nb));
    auto tv = new QGroupBox(tr("固定增益的非负 TV 反演"));
    auto tf = new QFormLayout(tv);
    addFlag(tf, "tv.enabled", tr("启用 TV 正则化"));
    addNullable(tf, "tv.lambda", tr("显式 λ"), tr("使用相对尺度规则计算 λ"));
    addNumber(tf, "tv.relative_fraction", tr("相对尺度比例"));
    addNumber(tf, "tv.max_iterations", tr("TV 迭代上限"), true);
    addNumber(tf, "tv.tolerance", tr("原始 / 对偶变化容差"));
    auto tvNote =
        new QLabel(tr("自动 λ = 比例 × P95(|AᵀD|) × 最小轴步长。梯度使用实际 mm 间距；以 NNLS "
                      "为初值，增益固定。TV 未收敛时主结果回退 NNLS，保留 TV 最后迭代并明确标记。"));
    tvNote->setWordWrap(true);
    tf->addRow(tvNote);
    auto tb = new QWidget;
    auto tbl = new QVBoxLayout(tb);
    tbl->addWidget(tv);
    tbl->addStretch();
    pl->addWidget(scrollBody(tb));
    m_tabs->addTab(m_parameters, tr("2 重建参数"));
    auto review = new QWidget;
    auto rvl = new QVBoxLayout(review);
    WorkflowUi::compact(rvl);
    auto batch = new QHBoxLayout;
    batch->addWidget(new QLabel(tr("ROI 保持固定像素顺序；支持 Ctrl / Shift 多选")));
    batch->addStretch();
    m_excludeButton = new QPushButton(tr("排除选中编号"));
    m_includeButton = new QPushButton(tr("撤销所选额外排除"));
    batch->addWidget(m_excludeButton);
    batch->addWidget(m_includeButton);
    rvl->addLayout(batch);
    auto split = new QSplitter(Qt::Horizontal);
    split->setChildrenCollapsible(false);
    m_rois = new QTableWidget(0, 8);
    m_rois->setObjectName("ImagingRois");
    m_rois->setHorizontalHeaderLabels({tr("原编号"), tr("标签 / 锚点键"), tr("状态"), tr("像素"), tr("PSF"),
                                       tr("RMS / pixel"), tr("增益"), tr("拒绝原因")});
    m_rois->setSelectionBehavior(QAbstractItemView::SelectRows);
    m_rois->setSelectionMode(QAbstractItemView::ExtendedSelection);
    m_rois->setAlternatingRowColors(true);
    m_rois->verticalHeader()->hide();
    m_rois->horizontalHeader()->setStretchLastSection(true);
    split->addWidget(m_rois);
    m_validationText = new QPlainTextEdit;
    m_validationText->setReadOnly(true);
    split->addWidget(m_validationText);
    split->setSizes({840, 420});
    rvl->addWidget(split, 1);
    m_tabs->addTab(review, tr("3 输入与标定检查"));
    connect(m_excludeButton, &QPushButton::clicked, this, [this] { exclusions(true); });
    connect(m_includeButton, &QPushButton::clicked, this, [this] { exclusions(false); });
    auto volume = new QWidget;
    auto vl = new QVBoxLayout(volume);
    WorkflowUi::compact(vl);
    auto tools = new QHBoxLayout;
    m_method = new QComboBox;
    m_method->addItem(tr("主结果"), 0);
    m_plane = new QComboBox;
    m_plane->setObjectName("ImagingPlane");
    m_plane->addItems({"XY", "XZ", "YZ"});
    m_plane->setCurrentIndex(2);
    m_display = new QComboBox;
    m_display->setObjectName("ImagingDisplay");
    m_display->addItems({tr("切片"), tr("沿法向求和投影")});
    m_slice = new QSpinBox;
    m_slice->setObjectName("ImagingSlice");
    m_slice->setKeyboardTracking(false);
    m_sliceLabel = new QLabel;
    tools->addWidget(m_method);
    tools->addWidget(m_plane);
    tools->addWidget(m_display);
    tools->addWidget(new QLabel(tr("切片索引")));
    tools->addWidget(m_slice);
    tools->addWidget(m_sliceLabel);
    auto log = new QCheckBox(tr("signed log1p 显示"));
    tools->addWidget(log);
    tools->addStretch();
    m_openButton = new QPushButton(tr("在主窗口查看 HDF5"));
    m_openButton->setObjectName("OpenImagingHdf");
    tools->addWidget(m_openButton);
    vl->addLayout(tools);
    m_image = new ScientificPlot;
    m_image->setObjectName("ImagingVolume");
    vl->addWidget(m_image, 1);
    auto cursor = new QLabel(tr("在图像上移动鼠标查看坐标与数值；滚轮缩放，双击复位。"));
    vl->addWidget(cursor);
    connect(m_image, &ScientificPlot::cursorChanged, cursor, &QLabel::setText);
    m_info = new QLabel(tr("完整重建完成后显示三维切片与投影。"));
    m_info->setWordWrap(true);
    vl->addWidget(m_info);
    m_tabs->addTab(volume, tr("4 三维结果"));
    connect(log, &QCheckBox::toggled, m_image, &ScientificPlot::setLogarithmic);
    for (auto c : {m_method, m_plane, m_display})
        connect(c, &QComboBox::currentIndexChanged, this, [this] { refreshImage(); });
    connect(m_slice, &QSpinBox::valueChanged, this, [this] { refreshImage(); });
    connect(m_openButton, &QPushButton::clicked, this, [this] {
        if (m_haveResult)
            emit openFileRequested(m_result.directory + "/reconstruction.h5");
    });
    auto diagnostic = new QTabWidget;
    m_sensitivity = new SpectrumPlot;
    m_sensitivity->setObjectName("ImagingSensitivity");
    diagnostic->addTab(m_sensitivity, tr("X 灵敏度"));
    auto trajectories = new QWidget;
    auto tl = new QVBoxLayout(trajectories);
    WorkflowUi::compact(tl);
    m_psfChoice = new QComboBox;
    m_psfChoice->setMinimumContentsLength(12);
    tl->addWidget(m_psfChoice, 0, Qt::AlignLeft);
    m_trajectory = new SpectrumPlot;
    tl->addWidget(m_trajectory, 1);
    diagnostic->addTab(trajectories, tr("质心轨迹"));
    m_spectrum = new SpectrumPlot;
    diagnostic->addTab(m_spectrum, tr("Gram 特征值"));
    m_log = new QPlainTextEdit;
    m_log->setReadOnly(true);
    m_log->setMaximumBlockCount(5000);
    m_log->setObjectName("ImagingLog");
    diagnostic->addTab(m_log, tr("运行日志"));
    m_summary = new QPlainTextEdit;
    m_summary->setReadOnly(true);
    diagnostic->addTab(m_summary, tr("完整摘要"));
    m_tabs->addTab(diagnostic, tr("5 诊断与日志"));
    connect(m_psfChoice, &QComboBox::currentIndexChanged, this, [this] { refreshTrajectory(); });
    m_status = new QLabel(tr("可载入已有方案，或选择薄膜、样品和 ROI 后校验输入。"));
    m_status->setObjectName("ImagingStatus");
    m_status->setWordWrap(true);
    m_status->setTextInteractionFlags(Qt::TextSelectableByMouse);
    m_progress = new QProgressBar;
    m_progress->setMaximumWidth(140);
    m_progress->hide();
    auto actions = new QHBoxLayout;
    m_validateButton = new QPushButton(tr("仅校验输入"));
    m_validateButton->setObjectName("ValidateImaging");
    m_runButton = new QPushButton(tr("开始完整重建并保存"));
    m_runButton->setObjectName("RunImaging");
    m_cancelButton = new QPushButton(tr("取消计算"));
    m_cancelButton->setObjectName("CancelImaging");
    m_cancelButton->setEnabled(false);
    m_previousButton = new QPushButton(tr("上一步"));
    m_previousButton->setObjectName("PreviousImagingStep");
    m_nextButton = new QPushButton(tr("下一步"));
    m_nextButton->setObjectName("NextImagingStep");
    actions->addWidget(m_previousButton);
    actions->addWidget(m_status, 1);
    actions->addWidget(m_progress);
    actions->addWidget(m_validateButton);
    actions->addWidget(m_runButton);
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
    refreshActions();
    connect(close, &QPushButton::clicked, this, &QWidget::hide);
    connect(m_validateButton, &QPushButton::clicked, this, [this] { start(true); });
    connect(m_runButton, &QPushButton::clicked, this, [this] { start(false); });
    connect(m_cancelButton, &QPushButton::clicked, this, [this] {
        m_cancel.store(true);
        status(tr("正在取消，等待当前读取或迭代结束…"));
    });
    connect(m_importButton, &QPushButton::clicked, this, [this] {
        auto path = QFileDialog::getOpenFileName(this, tr("载入已有方案"), {}, "JSON (*.json)");
        if (!path.isEmpty())
            importConfiguration(path);
    });
    connect(m_exportButton, &QPushButton::clicked, this, [this] {
        try {
            auto config = configuration();
            auto path = QFileDialog::getSaveFileName(this, tr("保存成像配置"), "imaging.config.json",
                                                     "JSON (*.json)");
            if (!path.isEmpty()) {
                ImageUtil::saveJson(path, config);
                status(tr("配置已保存：%1").arg(path));
            }
        } catch (const std::exception &e) {
            status(QString::fromUtf8(e.what()), true);
        }
    });
    connect(load, &QPushButton::clicked, this, [this] {
        if (m_busy)
            return;
        auto directory = QFileDialog::getExistingDirectory(this, tr("选择包含 reconstruction.h5 的结果目录"));
        if (!directory.isEmpty())
            loadResult(directory);
    });
}
QJsonObject ImagingDialog::configuration() const {
    auto c = m_base;
    for (auto i = m_fields.begin(); i != m_fields.end(); ++i) {
        if (m_auto.contains(i.key()) && m_auto[i.key()]->isChecked()) {
            put(c, i.key(), QJsonValue());
            continue;
        }
        put(c, i.key(), ImagingInputUi::read(i.value(), i.key()));
    }
    return Imaging::config(c, QDir::currentPath());
}
void ImagingDialog::applyConfiguration(const QJsonObject &config) {
    if (m_busy)
        return;
    m_updating = true;
    m_base = config;
    for (auto i = m_fields.begin(); i != m_fields.end(); ++i) {
        auto value = fieldValue(config, i.key());
        ImagingInputUi::write(i.value(), value);
        if (m_auto.contains(i.key())) {
            m_auto[i.key()]->setChecked(value.isNull());
            i.value()->setEnabled(!value.isNull());
        }
    }
    m_updating = false;
    invalidate();
}
bool ImagingDialog::importConfiguration(const QString &path) {
    if (m_busy)
        return false;
    try {
        auto c = Imaging::loadConfig(path);
        applyConfiguration(c);
        status(tr("已载入方案 %1。请核对原文件路径并选择新的输出目录。").arg(path));
        return true;
    } catch (const std::exception &e) {
        status(QString::fromUtf8(e.what()), true);
        return false;
    }
}
void ImagingDialog::invalidate() {
    if (m_updating || m_busy)
        return;
    QFile roi(qobject_cast<QLineEdit *>(m_fields["input.roi_table"])->text());
    m_roiBytes = !roi.fileName().isEmpty() && roi.open(QIODevice::ReadOnly) ? roi.readAll() : QByteArray();
    m_haveResult = false;
    m_result = ImagingResult();
    m_image->setData({});
    m_info->hide();
    m_summary->clear();
    m_sensitivity->setCurves({}, {});
    m_spectrum->setCurves({}, {});
    m_trajectory->setCurves({}, {});
    m_psfChoice->clear();
    m_openButton->setEnabled(false);
    m_excludeButton->setEnabled(false);
    m_includeButton->setEnabled(false);
    for (int row = 0; row < m_rois->rowCount(); ++row)
        for (int col = 4; col < 8; ++col)
            if (m_rois->item(row, col))
                m_rois->item(row, col)->setText(QStringLiteral("—"));
    status(tr("参数已变更；按当前设置重新校验或重建。"));
}
void ImagingDialog::start(bool validateOnly) {
    if (m_busy)
        return;
    try {
        auto c = configuration();
        invalidate();
        m_validateOnly = validateOnly;
        m_busy = true;
        m_cancel.store(false);
        m_inputs->setEnabled(false);
        m_parameters->setEnabled(false);
        for (auto b : {m_validateButton, m_runButton, m_importButton, m_exportButton, m_includeButton,
                       m_excludeButton})
            b->setEnabled(false);
        m_cancelButton->setEnabled(true);
        m_progress->setValue(0);
        m_progress->show();
        m_log->clear();
        status(validateOnly ? tr("正在校验输入…") : tr("正在运行完整成像流程…"));
        auto progress = [this](const QString &text, int percent) {
            QMetaObject::invokeMethod(
                this,
                [this, text, percent] {
                    if (m_busy) {
                        m_progress->setValue(percent);
                        status(text);
                        m_log->appendPlainText(text);
                    }
                },
                Qt::QueuedConnection);
        };
        m_watcher.setFuture(QtConcurrent::run([this, c, validateOnly, progress] {
            Work work;
            try {
                work.validation = Imaging::validate(c, m_cancel);
                if (!validateOnly)
                    work.result = Imaging::run(c, m_cancel, progress);
            } catch (const std::exception &e) {
                work.error = QString::fromUtf8(e.what());
            }
            return work;
        }));
    } catch (const std::exception &e) {
        status(QString::fromUtf8(e.what()), true);
    }
}
void ImagingDialog::refreshValidation() {
    m_excludeButton->setEnabled(true);
    m_includeButton->setEnabled(true);
    m_rois->setRowCount(m_validation.rois.size());
    for (int row = 0; row < m_validation.rois.size(); ++row) {
        const auto &roi = m_validation.rois[row];
        QStringList values{roi.id,
                           roi.key,
                           roi.excluded ? tr("排除") : tr("参与标定"),
                           QString::number(roi.pixels.size()),
                           QStringLiteral("—"),
                           QStringLiteral("—"),
                           QStringLiteral("—"),
                           QString()};
        for (int col = 0; col < 8; ++col)
            m_rois->setItem(row, col, cell(values[col]));
    }
    m_rois->resizeColumnsToContents();
    auto describeLoops = [this](const QString &key, const QString &label) {
        auto report = m_validation.report[key].toObject();
        auto order = report["loop_order_outer_to_inner"].toArray();
        auto counts = report["loop_counts_outer_to_inner"].toArray();
        QStringList layers;
        for (qsizetype i = 0; i < order.size(); ++i)
            layers.append(tr("%1（%2 点）").arg(order[i].toString()).arg(counts[i].toInt()));
        return tr("%1（外 → 内） %2 · 共 %3 帧")
            .arg(label, layers.join(QStringLiteral(" → ")))
            .arg(report["actual_frames"].toInt());
    };
    m_validationText->setPlainText(describeLoops("calibration", tr("薄膜标定")) + '\n' +
                                   describeLoops("sample", tr("样品成像")) + "\n\n" +
                                   QString::fromUtf8(QJsonDocument(m_validation.report).toJson()));
}
void ImagingDialog::exclusions(bool exclude) {
    if (m_busy)
        return;
    auto edit = qobject_cast<QLineEdit *>(m_fields["roi.excluded_ids"]);
    QSet<QString> ids;
    for (const auto &id : edit->text().split(QRegularExpression("[,，;；\\s]+"), Qt::SkipEmptyParts))
        ids.insert(id);
    for (const auto &index : m_rois->selectionModel()->selectedRows()) {
        auto id = m_rois->item(index.row(), 0)->text();
        if (exclude)
            ids.insert(id);
        else
            ids.remove(id);
    }
    QStringList values(ids.begin(), ids.end());
    std::sort(values.begin(), values.end(),
              [](const QString &a, const QString &b) { return a.toLongLong() < b.toLongLong(); });
    edit->setText(values.join(','));
    status(tr("额外排除列表已更新。共用 ROI 中的人工失效仍保持排除；请重新校验。"));
}
void ImagingDialog::refreshResult() {
    m_openButton->setEnabled(true);
    m_summary->setPlainText(QString::fromUtf8(QJsonDocument(m_result.summary).toJson()));
    {
        QSignalBlocker b(m_method);
        m_method->clear();
        m_method->addItem(tr("主结果"), 0);
        m_method->addItem("NNLS", 1);
        if (m_result.tv.size())
            m_method->addItem(m_result.summary["tv_converged"].toBool() ? tr("TV（已收敛）")
                                                                        : tr("TV 最后迭代（未收敛）"),
                              2);
    }
    auto axes = m_result.axes;
    {
        QSignalBlocker block(m_slice);
        int normal = m_plane->currentIndex() == 0 ? 2 : m_plane->currentIndex() == 1 ? 1 : 0;
        m_slice->setRange(0, int(axes[normal].size()) - 1);
        m_slice->setValue(int(axes[normal].size()) / 2);
    }
    auto &summary = m_result.summary;
    m_info->setText(tr("%1 × %2 × %3 · %4 · 相对 Kapton 等效密度 · NNLS 残差 %5 · TV %6 · "
                       "样品背景未扣除；完整范围保留，无几何后处理。")
                        .arg(axes[0].size())
                        .arg(axes[1].size())
                        .arg(axes[2].size())
                        .arg(summary["solver"].toString(),
                             real(summary["nnls_relative_l2_residual"].toDouble()),
                             !m_result.tv.size()                ? tr("关闭")
                             : summary["tv_converged"].toBool() ? tr("收敛")
                                                                : tr("未收敛，主结果回退 NNLS")));
    m_info->show();
    m_sensitivity->setCurves({{"Sensitivity", axes[0], m_result.sensitivity, QColor(0, 80, 160), false}},
                             tr("完整 X 范围的灵敏度（含增益和 Δx / 薄膜厚度）"), "X (mm)", "sqrt(sum A^2)");
    QVector<double> index;
    for (int i = 0; i < m_result.eigenvalues.size(); ++i)
        index.append(i);
    m_spectrum->setCurves({{"Gram eigenvalues", index, m_result.eigenvalues, QColor(0, 80, 160), false}},
                          tr("Gram 特征值 · 数值秩 %1 / %2 · 相对阈值 1e−14")
                              .arg(summary["numerical_response_rank"].toInt())
                              .arg(axes[0].size()),
                          "Mode index", "Eigenvalue");
    {
        QSignalBlocker b(m_psfChoice);
        m_psfChoice->clear();
        for (auto i = m_result.models.begin(); i != m_result.models.end(); ++i)
            m_psfChoice->addItem(i.key());
    }
    for (int row = 0; row < m_rois->rowCount(); ++row) {
        auto key = m_rois->item(row, 1)->text();
        if (m_result.models.contains(key)) {
            auto model = m_result.models[key].toObject();
            m_rois->item(row, 4)->setText(tr("合格"));
            m_rois->item(row, 5)->setText(real(model["centroid_fit_rms_detector_pixels"].toDouble()));
            m_rois->item(row, 6)->setText(real(summary["gains_by_roi"].toObject()[key].toDouble()));
        } else if (m_result.rejected.contains(key)) {
            m_rois->item(row, 4)->setText(tr("拒绝"));
            m_rois->item(row, 7)->setText(m_result.rejected[key].toString());
        }
    }
    for (int column : {4, 5, 6})
        m_rois->resizeColumnToContents(column);
    refreshTrajectory();
    refreshImage();
}
void ImagingDialog::refreshImage() {
    if (!m_haveResult)
        return;
    int method = m_method->currentData().toInt(), plane = m_plane->currentIndex();
    const auto &rho = method == 1 ? m_result.nnls : method == 2 ? m_result.tv : m_result.primary;
    if (!rho.size())
        return;
    const auto &axes = m_result.axes;
    int nx = int(axes[0].size()), ny = int(axes[1].size()), nz = int(axes[2].size()),
        normal = plane == 0   ? 2
                 : plane == 1 ? 1
                              : 0;
    {
        QSignalBlocker b(m_slice);
        m_slice->setRange(0, int(axes[normal].size()) - 1);
    }
    int fixed = m_slice->value();
    m_slice->setEnabled(m_display->currentIndex() == 0);
    m_sliceLabel->setText(QString("%1 = %2 mm").arg(QString("XYZ")[normal]).arg(real(axes[normal][fixed])));
    int colAxis = plane == 2 ? 1 : 0, rowAxis = plane == 0 ? 1 : 2;
    DataSlice preview;
    preview.path = QString("%1 · %2").arg(m_method->currentText(), m_plane->currentText());
    preview.width = int(axes[colAxis].size());
    preview.height = int(axes[rowAxis].size());
    preview.xAxis = 1;
    preview.yAxis = 0;
    preview.dims = {quint64(preview.height), quint64(preview.width)};
    preview.strides = {1, 1};
    preview.starts = {0, 0};
    preview.values.fill(0, preview.width * preview.height);
    for (int row = 0; row < preview.height; ++row)
        for (int col = 0; col < preview.width; ++col) {
            int position[3] = {0, 0, 0};
            position[colAxis] = col;
            position[rowAxis] = row;
            double value = 0;
            if (m_display->currentIndex() == 0) {
                position[normal] = fixed;
                value = rho(position[0], position[1] * nz + position[2]);
            } else
                for (int k = 0; k < axes[normal].size(); ++k) {
                    position[normal] = k;
                    value += rho(position[0], position[1] * nz + position[2]);
                }
            preview.values[row * preview.width + col] = value;
        }
    preview.minimum = *std::min_element(preview.values.begin(), preview.values.end());
    preview.maximum = *std::max_element(preview.values.begin(), preview.values.end());
    preview.finiteCount = preview.values.size();
    preview.sum = std::accumulate(preview.values.begin(), preview.values.end(), 0.);
    preview.mean = preview.sum / preview.values.size();
    m_image->setData(preview);
    m_image->setImageCoordinates(
        axes[colAxis].first(), (axes[colAxis].last() - axes[colAxis].first()) / (axes[colAxis].size() - 1),
        axes[rowAxis].first(), (axes[rowAxis].last() - axes[rowAxis].first()) / (axes[rowAxis].size() - 1));
    m_image->setAxisLabels(QString("%1 (mm)").arg(QString("XYZ")[colAxis]),
                           QString("%1 (mm)").arg(QString("XYZ")[rowAxis]));
    Q_UNUSED(nx);
    Q_UNUSED(ny);
}
void ImagingDialog::refreshTrajectory() {
    if (!m_haveResult || m_psfChoice->currentText().isEmpty())
        return;
    auto model = m_result.models[m_psfChoice->currentText()].toObject();
    QVector<double> x, row, column, fitRow, fitCol;
    auto r = model["row_line"].toArray(), c = model["column_line"].toArray();
    for (auto value : model["observations"].toArray()) {
        auto obs = value.toObject();
        double pos = obs["x"].toDouble();
        x.append(pos);
        row.append(obs["row"].toDouble());
        column.append(obs["column"].toDouble());
        fitRow.append(r[0].toDouble() * pos + r[1].toDouble());
        fitCol.append(c[0].toDouble() * pos + c[1].toDouble());
    }
    m_trajectory->setCurves({{"Row centroid", x, row, QColor(0, 80, 160), false},
                             {"Column centroid", x, column, QColor(175, 50, 35), false},
                             {"Row trajectory", x, fitRow, QColor(0, 80, 160), true},
                             {"Column trajectory", x, fitCol, QColor(175, 50, 35), true}},
                            m_psfChoice->currentText() +
                                QString(" · RMS %1 pixel · sample %2 mm/pixel")
                                    .arg(real(model["centroid_fit_rms_detector_pixels"].toDouble()),
                                         real(model["sample_mm_per_detector_pixel"].toDouble())),
                            "Foil X (mm)", "Detector pixel");
}
bool ImagingDialog::loadResult(const QString &directory) {
    if (m_busy)
        return false;
    try {
        ImagingResult result;
        result.directory = QFileInfo(directory).absoluteFilePath();
        QFile state(result.directory + "/run_status.json");
        if (state.open(QIODevice::ReadOnly) &&
            !QJsonDocument::fromJson(state.readAll()).object()["completed"].toBool())
            ImageUtil::fail("This output is marked incomplete; inspect run.log and intermediate files "
                            "instead of treating it as a completed reconstruction");
        ImagingH5 f(result.directory + "/reconstruction.h5");
        for (auto name : {"x_mm", "y_mm", "z_mm"}) {
            auto axis = f.read(name);
            if (axis.isEmpty())
                ImageUtil::fail("Empty result axis");
            result.axes.append(axis);
        }
        auto read = [&](const QString &path) {
            ImagingH5 input(path);
            QVector<quint64> shape;
            auto values = input.read("reinterp_insect_ck", &shape);
            if (shape != QVector<quint64>{quint64(result.axes[0].size()), quint64(result.axes[1].size()),
                                          quint64(result.axes[2].size())})
                ImageUtil::fail("Result XYZ shape mismatch");
            ImageMatrix rho = Eigen::Map<ImageMatrix>(values.data(), result.axes[0].size(),
                                                      result.axes[1].size() * result.axes[2].size());
            if (!rho.allFinite() || rho.minCoeff() < 0)
                ImageUtil::fail("Result density must be finite/nonnegative");
            return rho;
        };
        result.primary = read(result.directory + "/reconstruction.h5");
        result.nnls = read(result.directory + "/reconstruction_nnls.h5");
        if (QFileInfo::exists(result.directory + "/reconstruction_tv.h5"))
            result.tv = read(result.directory + "/reconstruction_tv.h5");
        result.sensitivity = f.read("sensitivity_x");
        result.eigenvalues = f.read("gram_eigenvalues");
        QFile summary(result.directory + "/summary.json");
        if (!summary.open(QIODevice::ReadOnly))
            ImageUtil::fail("Cannot read imaging summary");
        result.summary = QJsonDocument::fromJson(summary.readAll()).object();
        if (result.summary.isEmpty())
            ImageUtil::fail("Invalid imaging summary");
        QFile psf(result.directory + "/psf_models.json");
        if (psf.open(QIODevice::ReadOnly)) {
            auto models = QJsonDocument::fromJson(psf.readAll()).object();
            result.models = models["models"].toObject();
            result.rejected = models["rejected"].toObject();
        }
        const auto savedConfig = result.directory + "/config.resolved.json";
        if (QFileInfo::exists(savedConfig)) {
            auto settings = Imaging::loadConfig(savedConfig);
            auto input = settings["input"].toObject();
            bool same = m_validation.report["sample"].toObject()["source_file"].toString() ==
                            input["sample_nxs"].toString() &&
                        m_validation.report["calibration"].toObject()["source_file"].toString() ==
                            input["calibration_nxs"].toString();
            applyConfiguration(settings);
            if (!same) {
                m_validation = ImagingValidation();
                m_rois->setRowCount(0);
                m_validationText->setPlainText(
                    tr("已载入保存时的配置与结果；可点击“仅校验输入”重新核对原始文件。"));
            }
        }
        m_result = std::move(result);
        m_haveResult = true;
        refreshResult();
        m_tabs->setCurrentIndex(3);
        status(tr("已只读载入已有成像结果：%1").arg(directory));
        return true;
    } catch (const std::exception &e) {
        status(QString::fromUtf8(e.what()), true);
        return false;
    }
}

void ImagingDialog::refreshActions() {
    const int step = m_tabs->currentIndex();
    m_configTools->setVisible(step == 0);
    m_previousButton->setVisible(step > 0);
    m_previousButton->setEnabled(!m_busy);
    m_nextButton->setVisible(step < 4);
    m_nextButton->setEnabled(!m_busy && (step != 2 || m_haveResult));
    m_validateButton->setVisible(step == 0 || step == 2);
    m_runButton->setVisible(step == 2);
    m_cancelButton->setVisible(m_busy);
    m_nextButton->setDefault(step < 2);
    m_runButton->setDefault(step == 2);
}
void ImagingDialog::useRoiDraft(const QString &file, const QString &path) {
    if (m_busy) {
        m_pendingRoi = {file, path};
        return;
    }
    auto sample = qobject_cast<QLineEdit *>(m_fields["input.sample_nxs"]);
    if (QFileInfo(sample->text()).canonicalFilePath() != file || !QFileInfo::exists(path))
        return;
    auto target = qobject_cast<QLineEdit *>(m_fields["input.roi_table"]);
    if (target->text() != path)
        target->setText(path);
    else {
        QFile roi(path);
        if (roi.open(QIODevice::ReadOnly) && roi.readAll() != m_roiBytes)
            invalidate();
    }
}
