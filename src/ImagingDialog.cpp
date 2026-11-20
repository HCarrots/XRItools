#include "ImagingDialog.h"
#include "AnimatedGif.h"
#include "BondFigures.h"
#include "ImagingIO.h"
#include "ImagingInputUi.h"
#include "RoiDialog.h"
#include "ScientificPlot.h"
#include "SpectrumPlot.h"
#include "WorkflowUi.h"
#include "WorkspacePaths.h"
#include <QSaveFile>
#include <QtEndian>
#include <QtConcurrent/QtConcurrentRun>
#include <QtWidgets>
#include <algorithm>
#include <limits>
#include <numeric>
namespace {
using WorkspacePaths::Kind;
using WorkspacePaths::Workflow;

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
    connect(&m_watcher, &QFutureWatcher<Work>::progressValueChanged, this, [this](int value) {
        if (m_exportingModelGif)
            m_progress->setValue(value);
    });
    connect(&m_watcher, &QFutureWatcher<Work>::finished, this, [this] {
        auto work = m_watcher.result();
        setBusy(false);
        m_progress->hide();
        if (m_diagnosticTask != DiagnosticTask::None) {
            const bool depth = m_diagnosticTask == DiagnosticTask::Depth;
            const bool boundaries = m_diagnosticTask == DiagnosticTask::Boundaries;
            m_diagnosticTask = DiagnosticTask::None;
            if (m_cancel.load() || !work.error.isEmpty()) {
                status(m_cancel.load() ? tr("已取消诊断计算。") : work.error, !m_cancel.load());
                return;
            }
            if (depth) {
                m_depthReport = work.diagnosticReport;
                refreshDepth();
            } else if (boundaries) {
                m_boundaryReport = work.diagnosticReport;
            } else {
                m_psfComparison = work.diagnosticReport;
                refreshTrajectory();
            }
            try {
                ImageUtil::saveJson(m_result.directory + (depth        ? "/depth_validation.json"
                                                          : boundaries ? "/boundary_stability.json"
                                                                       : "/psf_model_comparison.json"),
                                    work.diagnosticReport);
                status(depth        ? tr("薄箔留出验证完成，报告已保存。")
                       : boundaries ? tr("边界复算完成，点击「查看边界对照」检查收敛状态和形貌。")
                                    : tr("光斑模型比较完成，已保存逐晶体报告。"));
            } catch (const std::exception &e) {
                status(tr("验证完成，报告可在当前窗口查看。%1").arg(QString::fromUtf8(e.what())), true);
            }
            return;
        }
        if (m_exportingModelGif) {
            m_exportingModelGif = false;
            m_progress->setRange(0, 100);
            m_cancelButton->setText(tr("取消计算"));
            if (!work.error.isEmpty())
                status(work.error, true);
            else if (work.gifPath.isEmpty())
                status(tr("已取消 GIF 导出。"));
            else
                status(tr("旋转 GIF 已保存到 %1").arg(QDir::toNativeSeparators(work.gifPath)));
            return;
        }
        if (m_buildingModel) {
            m_buildingModel = false;
            m_progress->setRange(0, 100);
            if (!work.error.isEmpty() || m_cancel.load()) {
                status(work.error.isEmpty() ? tr("已取消等值面提取。") : work.error, true);
                return;
            }
            if (work.mesh.isEmpty()) {
                m_model->setMesh(m_result.axes, {}, {}, tr("当前阈值没有等值面，请调整阈值。"));
                status(tr("当前阈值没有等值面，请调整阈值。"));
                return;
            }
            m_model->setMesh(m_result.axes, work.mesh, {},
                             tr("%1 · 阈值 %2 · 左键拖动旋转 · 滚轮缩放 · 双击复位")
                                 .arg(m_modelSource->currentText(), m_modelLevel->text().trimmed()),
                             {QColor("#4186ae"), QColor()});
            status(tr("等值面已生成，共 %1 个三角形。").arg(work.mesh.size()));
            return;
        }
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
    m_tabs->setCurrentIndex(std::clamp(step, 0, m_tabs->count() - 1));
}
void ImagingDialog::buildUi() {
    auto root = new QVBoxLayout(this);
    WorkflowUi::compact(root);
    m_importButton = WorkflowUi::button("folder-open", tr("载入已有方案…"));
    m_importButton->setObjectName("ImportImagingConfig");
    m_exportButton = WorkflowUi::button("save", tr("保存当前配置…"));
    m_exportButton->setObjectName("SaveImagingConfig");
    auto load = WorkflowUi::button("folder-open", tr("查看已有重建…"));
    m_configTools = WorkflowUi::cornerTools({m_importButton, m_exportButton, load});
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
            const auto draft = dialog.currentSchemePath();
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
        if (output) {
            m_chooseOutput = choose;
            m_workspaceOutput = new QLabel;
            m_workspaceOutput->setObjectName("ImagingWorkspaceOutput");
            m_workspaceOutput->setWordWrap(true);
            m_workspaceOutput->setTextInteractionFlags(Qt::TextSelectableByMouse);
            l->addWidget(m_workspaceOutput, 1);
        }
        form->addRow(label, row);
        connect(edit, &QLineEdit::textChanged, this, &ImagingDialog::invalidate);
        connect(choose, &QPushButton::clicked, this, [this, edit, output] {
            QString path;
            if (output) {
                auto parent = QFileDialog::getExistingDirectory(this, tr("选择结果的父目录"),
                                                                QFileInfo(edit->text()).absolutePath());
                if (!parent.isEmpty())
                    path = QDir(parent).filePath(
                        "imaging_" + QDateTime::currentDateTime().toString("yyyyMMdd_HHmmss_zzz"));
            } else {
                path = QFileDialog::getOpenFileName(this, tr("选择独立空白扫描"),
                                                   QFileInfo(edit->text()).absolutePath(),
                                                   "NeXus / HDF5 (*.nxs *.h5 *.hdf5)");
            }
            if (!path.isEmpty())
                edit->setText(path);
        });
    };
    auto output = new QGroupBox(tr("结果与诊断"));
    auto of = new QFormLayout(output);
    fileField(of, "runtime.output_dir", tr("新输出目录"), true);
    refreshOutputDestination();
    addFlag(of, "output.save_diagnostics", tr("保存各晶体正演矩阵、预测和残差"));
    addFlag(of, "output.plots", tr("保存诊断 PNG"));
    addFlag(of, "output.save_avizo", tr("同时保存 Avizo 体数据（.am）"));
    auto outputNote = new QLabel(tr("完整流程自动保存主结果、NNLS/TV、XYZ "
                                    "坐标、像素缓存、PSF、灵敏度、残差与参数。已存在的输出目录会拒绝写入。取"
                                    "消或失败时保留中间文件与日志。"));
    outputNote->setWordWrap(true);
    of->addRow(outputNote);
    rr->addWidget(output);
    auto boundary = new QLabel(
        tr("主输出：相对 Kapton 等效弹性散射密度 "
           "I(X,Y,Z)。能量已积分；背景来源可在重建参数中指定。"));
    boundary->setWordWrap(true);
    rr->addWidget(boundary);
    rr->addStretch();
    inputLayout->addWidget(scrollBody(right), 2);
    m_tabs->addTab(m_inputs, tr("1 文件与 ROI"));
    m_parameters = new QWidget;
    auto pl = new QHBoxLayout(m_parameters);
    auto psf = new QGroupBox(tr("晶体光斑模型（PSF）"));
    auto pf = new QFormLayout(psf);
    addFlag(pf, "psf.auto_select_model", tr("用薄箔留出误差自动选择模型"));
    addNumber(pf, "psf.min_peak_fraction", tr("最小相对峰高"));
    addNumber(pf, "psf.min_r_squared", tr("最低拟合 R²"));
    addNumber(pf, "psf.edge_margin_pixels", tr("中心距边界 / pixel"));
    addNumber(pf, "psf.min_observations", tr("最少合格观测"), true);
    addNumber(pf, "psf.max_trajectory_rms_pixels", tr("轨迹最大 RMS / pixel"));
    addText(pf, "psf.sigma_bounds_pixels", tr("σ 边界 [min,max]"), 4);
    addNumber(pf, "psf.max_fit_evaluations", tr("单幅拟合调用上限"), true);
    auto psfNote =
        new QLabel(tr("每次留出一个薄箔位置，比较整幅 ROI 的预测误差。至少 8 "
                      "个合格位置且光斑有明显不对称时，加入实测模板。关闭自动选择可使用恒定高斯。"));
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
    auto statistics = new QGroupBox(tr("计数方差与背景"));
    auto sf = new QFormLayout(statistics);
    addFlag(sf, "statistics.variance_weighted", tr("按传播计数方差加权（NNLS、增益与 TV）"));
    addNullable(sf, "statistics.sigma_floor", tr("σ 下限（积分信号单位）"),
                tr("自动：方差下限取各晶体正方差均值的 1%"));
    fileField(sf, "statistics.blank_nxs", tr("独立空白扫描"));
    addNullable(sf, "statistics.blank_exposure_s", tr("空白曝光 / s"), tr("使用样品曝光时间"));
    addText(sf, "statistics.background_yz_mm", tr("背景区 [Ymin,Ymax,Zmin,Zmax] / mm"), 4);
    auto bgHint = new QLabel(tr("空白扫描按能量汇总整份扫描的计数与剂量。也可指定样品扫描中已知为空白的 YZ "
                                "区域，至少覆盖两个位置。"));
    bgHint->setWordWrap(true);
    sf->addRow(bgHint);
    nbl->addWidget(statistics);
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
    m_tabs->addTab(volume, tr("4 剖面分析"));
    connect(log, &QCheckBox::toggled, m_image, &ScientificPlot::setLogarithmic);
    for (auto c : {m_method, m_plane, m_display})
        connect(c, &QComboBox::currentIndexChanged, this, [this] { refreshImage(); });
    connect(m_slice, &QSpinBox::valueChanged, this, [this] { refreshImage(); });
    connect(m_openButton, &QPushButton::clicked, this, [this] {
        if (m_haveResult)
            emit openFileRequested(m_result.directory + "/reconstruction.h5");
    });

    auto modelPage = new QWidget;
    auto modelLayout = new QVBoxLayout(modelPage);
    WorkflowUi::compact(modelLayout);
    auto modelTools = new QHBoxLayout;
    m_modelSource = new QComboBox;
    m_modelSource->setObjectName("ImagingModelSource");
    m_modelLevel = new QLineEdit;
    m_modelLevel->setObjectName("ImagingModelThreshold");
    m_modelLevel->setMaximumWidth(150);
    m_modelLevel->setToolTip(tr("输入相对密度，可使用科学计数法。"));
    m_modelRange = new QLabel;
    m_buildModelButton = WorkflowUi::button("layers", tr("生成等值面"));
    m_buildModelButton->setObjectName("BuildImagingModel");
    m_saveModelGifButton = WorkflowUi::button("download", tr("保存旋转 GIF…"));
    m_saveModelGifButton->setObjectName("SaveImagingModelGif");
    m_saveModelGifButton->setAutoDefault(false);
    m_saveModelGifButton->setToolTip(tr("从当前视角绕 Z 轴旋转，只导出白底模型。"));
    modelTools->addWidget(new QLabel(tr("体数据")));
    modelTools->addWidget(m_modelSource);
    modelTools->addWidget(new QLabel(tr("密度阈值")));
    modelTools->addWidget(m_modelLevel);
    modelTools->addWidget(m_buildModelButton);
    modelTools->addWidget(m_modelRange);
    modelTools->addStretch();
    modelTools->addWidget(m_saveModelGifButton);
    modelLayout->addLayout(modelTools);
    m_model = new BondFigures;
    m_model->setObjectName("ImagingModel");
    m_model->setMesh({}, {}, {}, tr("完成重建或载入已有结果后，可以生成三维等值面。"));
    modelLayout->addWidget(m_model, 1);
    auto depthTools = new QHBoxLayout;
    m_showLow = new QCheckBox(tr("半透明标记低响应区"));
    m_showLow->setChecked(true);
    m_lowPercent = new QDoubleSpinBox;
    m_lowPercent->setRange(0., 100.);
    m_lowPercent->setDecimals(1);
    m_lowPercent->setValue(1.);
    m_lowPercent->setSuffix("%");
    m_lowPercent->setToolTip(tr("以全 X 峰值的百分比定义低响应区，仅影响显示。"));
    m_validateDepth = new QPushButton(tr("深度分辨率验证"));
    m_validateDepth->setObjectName("ValidateImagingDepth");
    m_validateDepth->setToolTip(tr("逐个留出薄箔位置进行重建，完成后在结果目录保存 depth_validation.json。"));
    m_showDepth = new QPushButton(tr("查看逐点报告"));
    m_showDepth->setObjectName("ImagingDepthReport");
    for (auto button : {m_validateDepth, m_showDepth})
        button->setAutoDefault(false);
    depthTools->addWidget(m_showLow);
    depthTools->addWidget(m_lowPercent);
    depthTools->addStretch();
    depthTools->addWidget(m_validateDepth);
    depthTools->addWidget(m_showDepth);
    modelLayout->addLayout(depthTools);
    m_depthText = new QLabel;
    m_depthText->setObjectName("ImagingDepthReadout");
    m_depthText->setWordWrap(true);
    modelLayout->addWidget(m_depthText);
    connect(m_showLow, &QCheckBox::toggled, this, &ImagingDialog::refreshDepth);
    connect(m_lowPercent, &QDoubleSpinBox::valueChanged, this, &ImagingDialog::refreshDepth);
    connect(m_model, &BondFigures::depthSelected, this, [this](int index) {
        m_depthIndex = index;
        refreshDepth();
    });
    connect(m_validateDepth, &QPushButton::clicked, this, [this] { runDiagnostic(DiagnosticTask::Depth); });
    connect(m_showDepth, &QPushButton::clicked, this, &ImagingDialog::showDepthReport);
    m_tabs->addTab(modelPage, tr("5 三维模型"));
    connect(m_modelSource, &QComboBox::currentIndexChanged, this, &ImagingDialog::refreshModelSource);
    connect(m_buildModelButton, &QPushButton::clicked, this, &ImagingDialog::buildModel);
    connect(m_saveModelGifButton, &QPushButton::clicked, this, &ImagingDialog::saveModelGif);
    connect(m_modelLevel, &QLineEdit::returnPressed, this, &ImagingDialog::buildModel);

    auto diagnostic = new QTabWidget;
    m_sensitivity = new SpectrumPlot;
    m_sensitivity->setObjectName("ImagingSensitivity");
    diagnostic->addTab(m_sensitivity, tr("X 灵敏度"));
    auto trajectories = new QWidget;
    auto tl = new QVBoxLayout(trajectories);
    WorkflowUi::compact(tl);
    m_psfChoice = new QComboBox;
    m_psfChoice->setMinimumContentsLength(12);
    auto modelRow = new QHBoxLayout;
    modelRow->addWidget(m_psfChoice);
    m_comparePsf = new QPushButton(tr("用已有标定比较模型"));
    m_comparePsf->setObjectName("CompareImagingPsf");
    m_comparePsf->setAutoDefault(false);
    modelRow->addWidget(m_comparePsf);
    modelRow->addStretch();
    tl->addLayout(modelRow);
    connect(m_comparePsf, &QPushButton::clicked, this, [this] { runDiagnostic(DiagnosticTask::Models); });
    m_psfSelection = new QLabel;
    m_psfSelection->setObjectName("ImagingPsfSelection");
    m_psfSelection->setWordWrap(true);
    m_psfSelection->setTextInteractionFlags(Qt::TextSelectableByMouse);
    tl->addWidget(m_psfSelection);
    m_trajectory = new SpectrumPlot;
    tl->addWidget(m_trajectory, 1);
    diagnostic->addTab(trajectories, tr("光斑模型与质心轨迹"));
    m_spectrum = new SpectrumPlot;
    diagnostic->addTab(m_spectrum, tr("Gram 特征值"));
    auto boundaries = new QWidget;
    auto bl = new QVBoxLayout(boundaries);
    auto guide = new QLabel(tr("按 X 灵敏度把晶体分成两组，在共同覆盖的深度范围内独立复算。比较锚点和 TV "
                               "强度对质心、半高厚度及局部突起的影响。"));
    guide->setWordWrap(true);
    bl->addWidget(guide);
    m_compareBounds = new QPushButton(tr("运行晶体分组与边界复算"));
    m_compareBounds->setObjectName("ValidateImagingBoundaries");
    m_showBounds = new QPushButton(tr("查看边界对照"));
    m_showBounds->setObjectName("ShowImagingBoundaries");
    for (auto button : {m_compareBounds, m_showBounds}) {
        button->setAutoDefault(false);
        bl->addWidget(button, 0, Qt::AlignLeft);
    }
    bl->addStretch();
    diagnostic->addTab(boundaries, tr("增益与边界稳定性"));
    connect(m_compareBounds, &QPushButton::clicked, this,
            [this] { runDiagnostic(DiagnosticTask::Boundaries); });
    connect(m_showBounds, &QPushButton::clicked, this, &ImagingDialog::showBoundaryReport);
    m_log = new QPlainTextEdit;
    m_log->setReadOnly(true);
    m_log->setMaximumBlockCount(5000);
    m_log->setObjectName("ImagingLog");
    diagnostic->addTab(m_log, tr("运行日志"));
    m_summary = new QPlainTextEdit;
    m_summary->setReadOnly(true);
    diagnostic->addTab(m_summary, tr("完整摘要"));
    m_tabs->addTab(diagnostic, tr("6 诊断与日志"));
    connect(m_psfChoice, &QComboBox::currentIndexChanged, this, [this] { refreshTrajectory(); });
    m_status = new QLabel(tr("可载入已有方案，或选择薄膜、样品和 ROI 后校验输入。"));
    m_status->setObjectName("ImagingStatus");
    m_status->setWordWrap(true);
    m_status->setTextInteractionFlags(Qt::TextSelectableByMouse);
    m_progress = new QProgressBar;
    m_progress->setMaximumWidth(140);
    m_progress->hide();
    auto actions = new QHBoxLayout;
    m_validateButton = WorkflowUi::button("check", tr("仅校验输入"));
    m_validateButton->setObjectName("ValidateImaging");
    m_runButton = WorkflowUi::button("play", tr("开始完整重建并保存"));
    m_runButton->setObjectName("RunImaging");
    m_cancelButton = WorkflowUi::button("square", tr("取消计算"));
    m_cancelButton->setObjectName("CancelImaging");
    m_cancelButton->setEnabled(false);
    m_previousButton = WorkflowUi::button("chevron-left", tr("上一步"));
    m_previousButton->setObjectName("PreviousImagingStep");
    m_nextButton = WorkflowUi::button("chevron-right", tr("下一步"));
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
        status(m_exportingModelGif ? tr("正在取消 GIF 导出…") : tr("正在取消，等待当前读取或迭代结束…"));
    });
    connect(m_importButton, &QPushButton::clicked, this, [this] {
        auto path = QFileDialog::getOpenFileName(this, tr("载入已有方案"),
                                                 WorkspacePaths::directory(Workflow::Imaging, Kind::Config),
                                                 "JSON (*.json)");
        if (!path.isEmpty())
            importConfiguration(path);
    });
    connect(m_exportButton, &QPushButton::clicked, this, [this] {
        try {
            auto config = configuration();
            auto path = WorkspacePaths::saveFile(this, Workflow::Imaging, Kind::Config,
                                                 fieldValue(config, "input.sample_nxs").toString(),
                                                 tr("保存成像配置"), "imaging.config.json", "JSON (*.json)");
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
        auto directory = QFileDialog::getExistingDirectory(
            this, tr("选择包含 reconstruction.h5 的结果目录"),
            WorkspacePaths::directory(Workflow::Imaging, Kind::Result));
        if (!directory.isEmpty())
            loadResult(directory);
    });
}
void ImagingDialog::refreshOutputDestination() {
    const auto folder = WorkspacePaths::directory(Workflow::Imaging, Kind::Result);
    m_fields["runtime.output_dir"]->setVisible(folder.isEmpty());
    m_chooseOutput->setVisible(folder.isEmpty());
    m_workspaceOutput->setVisible(!folder.isEmpty());
    m_workspaceOutput->setText(QDir::toNativeSeparators(folder));
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
    const auto output = WorkspacePaths::newPath(Workflow::Imaging, Kind::Result,
                                                fieldValue(c, "input.sample_nxs").toString(), "imaging");
    if (!output.isEmpty())
        put(c, "runtime.output_dir", output);
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
    refreshModelSource();
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
        setBusy(true);
        m_cancel.store(false);
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
    m_boundaryReport={};
    QFile savedBounds(m_result.directory+"/boundary_stability.json");
    if(savedBounds.open(QIODevice::ReadOnly)){
        const auto report=QJsonDocument::fromJson(savedBounds.readAll()).object();
        if(report["schema"]=="XRI.boundary-stability.v1" &&
           report["shape"].toArray()==m_result.summary["shape_xyz"].toArray()) m_boundaryReport=report;
    }
    m_psfComparison = {};
    QFile comparison(m_result.directory + "/psf_model_comparison.json");
    if (comparison.open(QIODevice::ReadOnly)) {
        const auto report = QJsonDocument::fromJson(comparison.readAll()).object();
        if (report["schema"] == "XRI.psf-comparison.v1")
            m_psfComparison = report;
    }
    m_depthReport = {};
    QFile report(m_result.directory + "/depth_validation.json");
    if (report.open(QIODevice::ReadOnly)) {
        const auto saved = QJsonDocument::fromJson(report.readAll()).object();
        if (saved["schema"] == "XRI.foil-loo.v1" && !saved["folds"].toArray().isEmpty())
            m_depthReport = saved;
    }
    m_relativeSensitivity.clear();
    const auto &s = m_result.sensitivity;
    if (s.size() == m_result.axes[0].size() && s.size() >= 2 &&
        std::all_of(s.cbegin(), s.cend(), [](double v) { return std::isfinite(v) && v >= 0.; })) {
        const double peak = *std::max_element(s.cbegin(), s.cend());
        if (peak > 0.)
            for (double v : s)
                m_relativeSensitivity.append(v / peak);
    }
    m_depthIndex = int(m_relativeSensitivity.size()) / 2;
    m_openButton->setEnabled(true);
    m_summary->setPlainText(QString::fromUtf8(QJsonDocument(m_result.summary).toJson()));
    for (auto *choice : {m_method, m_modelSource}) {
        QSignalBlocker block(choice);
        choice->clear();
        choice->addItem(tr("主结果"), 0);
        choice->addItem("NNLS", 1);
        if (m_result.tv.size())
            choice->addItem(m_result.summary["tv_converged"].toBool() ? tr("TV（已收敛）")
                                                                      : tr("TV 最后迭代（未收敛）"),
                            2);
    }
    refreshModelSource();
    auto axes = m_result.axes;
    {
        QSignalBlocker block(m_slice);
        int normal = m_plane->currentIndex() == 0 ? 2 : m_plane->currentIndex() == 1 ? 1 : 0;
        m_slice->setRange(0, int(axes[normal].size()) - 1);
        m_slice->setValue(int(axes[normal].size()) / 2);
    }
    auto &summary = m_result.summary;
    m_info->setText(tr("%1 × %2 × %3 · %4 · 相对 Kapton 等效密度 · NNLS 残差 %5 · TV %6 · "
                       "%7")
                        .arg(axes[0].size())
                        .arg(axes[1].size())
                        .arg(axes[2].size())
                        .arg(summary["solver"].toString(),
                             real(summary["nnls_relative_l2_residual"].toDouble()),
                             !m_result.tv.size()                ? tr("关闭")
                             : summary["tv_converged"].toBool() ? tr("收敛")
                                                                : tr("未收敛，主结果回退 NNLS"),
                             summary["sample_background_subtracted"].toBool()?tr("已校正独立背景"):tr("未指定样品背景")));
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
    const auto &rho = volumeData(method);
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

const ImageMatrix &ImagingDialog::volumeData(int method) const {
    return method == 1 ? m_result.nnls : method == 2 ? m_result.tv : m_result.primary;
}

void ImagingDialog::refreshModelSource() {
    m_modelMin = m_modelMax = 0.;
    if (!m_haveResult) {
        m_depthIndex = -1;
        m_relativeSensitivity.clear();
        m_depthReport = {};
        refreshDepth();
        m_modelLevel->clear();
        m_modelRange->clear();
        m_model->setMesh({}, {}, {}, tr("完成重建或载入已有结果后，可以生成三维等值面。"));
        return;
    }
    const auto &volume = volumeData(m_modelSource->currentData().toInt());
    m_modelMin = volume.minCoeff();
    m_modelMax = volume.maxCoeff();
    const double threshold = m_modelMin + .25 * (m_modelMax - m_modelMin);
    m_modelLevel->setText(QString::number(threshold, 'g', 17));
    m_modelRange->setText(tr("密度范围 %1 至 %2").arg(real(m_modelMin), real(m_modelMax)));
    m_model->setMesh(m_result.axes, {}, {}, m_modelMin < m_modelMax
                                     ? tr("输入密度阈值后，点击「生成等值面」。")
                                     : tr("当前体数据的密度恒定，无法提取等值面。"));
    refreshDepth();
    refreshActions();
}

QString ImagingDialog::depthMethod(const QJsonObject &row) const {
    const int source = m_modelSource->currentData().toInt();
    return source == 2 || (source == 0 && row["tv_converged"].toBool()) ? "tv" : "nnls";
}

QString ImagingDialog::depthQuality(const QJsonObject &row) const {
    if (row.contains("error"))
        return row["error"].toString();
    const auto method = depthMethod(row);
    const auto peak = row[method].toObject();
    QStringList flags;
    if (m_modelSource->currentData().toInt() == 0 && row.contains("tv") && !row["tv_converged"].toBool())
        flags << tr("TV 未收敛，主结果采用 NNLS");
    if (row["relative_sensitivity"].isDouble() &&
        row["relative_sensitivity"].toDouble() < .01 * m_lowPercent->value())
        flags << tr("已知位置响应低于显示阈值");
    if (!row["gain_converged"].toBool() || !row[method + "_converged"].toBool())
        flags << tr("拟合未收敛");
    if (!peak["peak_x_mm"].isDouble())
        flags << tr("未检出峰");
    else {
        if (!peak["fwhm_mm"].isDouble())
            flags << tr("半高处超出网格");
        if (peak["multiple_lobes"].toBool())
            flags << tr("存在多个半高峰区");
    }
    return flags.join(QStringLiteral(" · "));
}

void ImagingDialog::refreshDepth() {
    m_lowPercent->setEnabled(!m_relativeSensitivity.isEmpty() && m_showLow->isChecked());
    m_model->setDepthProfile(m_relativeSensitivity, m_depthIndex,
                             m_showLow->isChecked() ? .01 * m_lowPercent->value() : 0.);
    if (!m_haveResult || m_depthIndex < 0 || m_depthIndex >= m_relativeSensitivity.size()) {
        m_depthText->setText(m_haveResult ? tr("没有可用的 X 灵敏度数据。") : QString());
        return;
    }
    const auto &x = m_result.axes[0];
    const int i = m_depthIndex;
    QString info =
        tr("X = %1 mm   ·   网格间距 %2 mm   ·   灵敏度 %3（峰值的 %4%）")
            .arg(real(x[i]), real(i ? x[i] - x[i - 1] : x[1] - x[0]), real(m_result.sensitivity[i]),
                 QString::number(100 * m_relativeSensitivity[i], 'g', 4));
    const auto folds = m_depthReport["folds"].toArray();
    if (folds.isEmpty())
        info += "\n" + tr("尚无薄箔定位实测结果。");
    else {
        QJsonObject nearest;
        double distance = std::numeric_limits<double>::infinity();
        for (const auto &item : folds) {
            const auto row = item.toObject();
            const double delta = std::abs(row["foil_x_mm"].toDouble() - x[i]);
            if (delta < distance) {
                nearest = row;
                distance = delta;
            }
        }
        auto metric = [](QJsonValue v) {
            return v.isDouble() ? QString::number(v.toDouble(), 'g', 5) : QStringLiteral("—");
        };
        const auto method = depthMethod(nearest);
        const auto peak = nearest[method].toObject();
        info += "\n" + tr("最近实测薄箔 %1 mm（相距 %2 mm） · %3 定位偏差 %4 mm · FWHM %5 mm")
                           .arg(metric(nearest["foil_x_mm"]), metric(distance), method.toUpper(),
                                metric(peak["bias_mm"]), metric(peak["fwhm_mm"]));
        const auto quality = depthQuality(nearest);
        if (!quality.isEmpty())
            info += " · " + quality;
    }
    m_depthText->setText(info);
}

void ImagingDialog::runDiagnostic(DiagnosticTask task) {
    if (m_busy || !m_haveResult)
        return;
    const QString dir = m_result.directory;
    m_diagnosticTask = task;
    m_cancel.store(false);
    setBusy(true);
    m_progress->setRange(0, 100);
    m_progress->setValue(0);
    m_progress->show();
    status(task==DiagnosticTask::Boundaries ? tr("正在读取分组重建数据…") : tr("正在读取薄箔标定数据…"));
    m_watcher.setFuture(QtConcurrent::run([this, dir, task] {
        Work work;
        try {
            const auto progress = [this](const QString &s, int pct) {
                QMetaObject::invokeMethod(
                    this,
                    [this, s, pct] {
                        if (m_diagnosticTask != DiagnosticTask::None) {
                            m_status->setText(s);
                            m_progress->setValue(pct);
                        }
                    },
                    Qt::QueuedConnection);
            };
            work.diagnosticReport = task == DiagnosticTask::Depth ? Imaging::validateDepth(dir, m_cancel, progress)
                            : task==DiagnosticTask::Boundaries ? Imaging::validateBoundaries(dir,m_cancel,progress)
                                                        : Imaging::comparePsfModels(dir, m_cancel, progress);
        } catch (const std::exception &e) {
            work.error = QString::fromUtf8(e.what());
        }
        return work;
    }));
}

void ImagingDialog::showBoundaryReport() {
    QDialog dialog(this);
    dialog.setWindowTitle(tr("晶体分组、增益与边界对照"));
    auto layout = new QVBoxLayout(&dialog);
    const auto runs = m_boundaryReport["runs"].toArray(), common = m_boundaryReport["common_x"].toArray();
    const auto &axes = m_result.axes;
    const int nx = int(axes[0].size()), ny = int(axes[1].size()), nz = int(axes[2].size());
    if (common.size() != nx || runs.isEmpty())
        return;
    auto note =
        new QLabel(tr("共同覆盖 %1/%2 个 X 位置。表中质心按强度加权，左右边界为轴向投影主峰的半高交点，单位 "
                      "mm。截断的交点留空。图中显示各次复算达到体峰值 50% 的区域。")
                       .arg(std::count(common.begin(), common.end(), QJsonValue(true)))
                       .arg(nx));
    note->setWordWrap(true);
    layout->addWidget(note);
    auto controls = new QHBoxLayout;
    auto membership = new QPushButton(tr("查看分组"));
    controls->addWidget(membership);
    connect(membership, &QPushButton::clicked, &dialog, [&] {
        QStringList lines;
        const auto members = m_boundaryReport["groups"].toArray();
        for (int g = 0; g < members.size(); ++g) {
            QStringList names;
            for (const auto &k : members[g].toArray())
                names << k.toString();
            lines << QString("%1 (%2)\n%3").arg(g ? "B" : "A").arg(names.size()).arg(names.join(", "));
        }
        QMessageBox::information(&dialog, tr("参与复算的晶体"), lines.join("\n\n"));
    });
    auto axis = new QComboBox;
    axis->addItems({"X", "Y", "Z"});
    controls->addWidget(axis);
    auto source = new QComboBox;
    const int trusted = m_boundaryReport["converged_runs"].toInt(),
              all = m_boundaryReport["signal_runs"].toInt();
    int fromA = 0;
    for (const auto &item : runs) {
        const auto row = item.toObject();
        if (row["group"] == "A" && row["converged"].toBool() && row["has_signal"].toBool())
            ++fromA;
    }
    note->setText(note->text() +
                  tr(" 可用于对照的收敛复算：A 组 %1 次，B 组 %2 次。").arg(fromA).arg(trusted - fromA));
    source->addItem(tr("收敛复算出现比例（%1 次）").arg(trusted));
    source->addItem(tr("全部试算出现比例（%1 次）").arg(all));
    if (!trusted)
        source->setCurrentIndex(1);
    for (const auto &v : runs) {
        const auto r = v.toObject();
        source->addItem(QString("%1 · %2 · TV ×%3")
                            .arg(r["group"].toString(), r["anchor"].toString("g=1"))
                            .arg(r["tv_factor"].toDouble()));
    }
    controls->addWidget(source);
    const qsizetype planeSize = qsizetype(ny) * nz;
    auto maskValues = [&](int index) -> QByteArray {
        const auto field = index < 2 ? m_boundaryReport.value(index == 0 ? "converged_votes" : "all_votes")
                                     : runs[index - 2].toObject().value("surface_mask");
        const auto compressed = QByteArray::fromBase64(field.toString().toLatin1());
        if (compressed.size() < 4 ||
            qFromBigEndian<quint32>(compressed.constData()) != quint64(nx) * planeSize)
            return {};
        return qUncompress(compressed);
    };
    auto slice = new QSpinBox;
    slice->setRange(0, nx - 1);
    const auto initialVotes = maskValues(source->currentIndex());
    if (initialVotes.size() == nx * planeSize) {
        qint64 bestSum = -1;
        for (int x = 0; x < nx; ++x) {
            if (!common[x].toBool())
                continue;
            qint64 sum = 0;
            for (qsizetype p = 0; p < planeSize; ++p)
                sum += uchar(initialVotes[x * planeSize + p]);
            if (sum > bestSum) {
                bestSum = sum;
                slice->setValue(x);
            }
        }
    }
    controls->addWidget(new QLabel(tr("X 切片")));
    controls->addWidget(slice);
    auto readout = new QLabel;
    controls->addWidget(readout);
    controls->addStretch();
    layout->addLayout(controls);
    auto table = new QTableWidget(int(runs.size()), 9);
    table->setHorizontalHeaderLabels({tr("组"), tr("锚点"), tr("TV 倍数"), "λ", tr("质心"), tr("左边界"),
                                      tr("右边界"), tr("厚度"), tr("状态")});
    auto deltaLabel = new QLabel;
    auto fill = [&] {
        std::array<QJsonObject, 2> reference;
        for (int i = 0; i < runs.size(); ++i) {
            const auto r = runs[i].toObject(), dims = r[axis->currentText()].toObject();
            QStringList flags;
            if (r["anchor"].isString() && !r["gain_converged"].toBool())
                flags << tr("增益未收敛");
            if (!r["nnls_converged"].toBool())
                flags << tr("NNLS 未收敛");
            if (r["tv_converged"].isBool() && !r["tv_converged"].toBool())
                flags << tr("TV 未收敛");
            if (!r["has_signal"].toBool())
                flags << tr("无信号");
            QStringList values{r["group"].toString(), r["anchor"].toString("g=1"),
                               real(r["tv_factor"].toDouble()), real(r["lambda"].toDouble())};
            for (auto key : {"center_mm", "left_mm", "right_mm", "thickness_mm"})
                values.append(dims[key].isDouble() ? QString::number(dims[key].toDouble(), 'g', 6)
                                                   : QStringLiteral("—"));
            values << (flags.isEmpty() ? tr("已收敛") : flags.join(", "));
            for (int j = 0; j < values.size(); ++j)
                table->setItem(i, j, cell(values[j]));
            int g = r["group"] == "A" ? 0 : 1;
            if (reference[g].isEmpty() && r["tv_factor"].toDouble() == 1.)
                reference[g] = r;
        }
        auto delta = [&](const QString &key) {
            auto a = reference[0][axis->currentText()].toObject().value(key);
            auto b = reference[1][axis->currentText()].toObject().value(key);
            return a.isDouble() && b.isDouble() ? QString::number(b.toDouble() - a.toDouble(), 'g', 5)
                                                : QStringLiteral("—");
        };
        deltaLabel->setText(tr("首选锚点、TV ×1，B 减 A：质心差 %1 mm，厚度差 %2 mm。%3")
                                .arg(delta("center_mm"), delta("thickness_mm"),
                                     reference[0]["converged"].toBool() && reference[1]["converged"].toBool()
                                         ? QString()
                                         : tr("该对照缺少两组均收敛的结果。")));
        table->resizeColumnsToContents();
        table->horizontalHeader()->setStretchLastSection(true);
    };
    layout->addWidget(table);
    layout->addWidget(deltaLabel);
    auto plot = new ScientificPlot;
    layout->addWidget(plot, 1);
    auto cursor = new QLabel;
    layout->addWidget(cursor);
    connect(plot, &ScientificPlot::cursorChanged, cursor, &QLabel::setText);
    auto redraw = [&] {
        const int s = source->currentIndex(), x = slice->value();
        const int n = s == 0 ? trusted : s == 1 ? all : 1;
        readout->setText(QString("X = %1 mm").arg(real(axes[0][x])));
        if (!common[x].toBool() || !n) {
            plot->setData({});
            readout->setText(readout->text() + " · " + tr("无共同覆盖或无可用复算"));
            return;
        }
        const auto values = maskValues(s);
        if (values.size() != qsizetype(nx) * ny * nz) {
            plot->setData({});
            return;
        }
        DataSlice d;
        d.width = ny;
        d.height = nz;
        d.xAxis = 1;
        d.yAxis = 0;
        d.dims = {quint64(nz), quint64(ny)};
        d.starts = {0, 0};
        d.strides = {1, 1};
        d.minimum = 0.;
        d.maximum = 1.;
        d.finiteCount = ny * nz;
        for (int z = 0; z < nz; ++z)
            for (int y = 0; y < ny; ++y)
                d.values.append(double(uchar(values[(x * ny + y) * nz + z])) / n);
        plot->setData(d);
        plot->setImageCoordinates(axes[1][0], axes[1][1] - axes[1][0], axes[2][0], axes[2][1] - axes[2][0]);
        plot->setAxisLabels("Y (mm)", "Z (mm)");
    };
    connect(axis, &QComboBox::currentIndexChanged, &dialog, fill);
    connect(source, &QComboBox::currentIndexChanged, &dialog, redraw);
    connect(slice, &QSpinBox::valueChanged, &dialog, redraw);
    connect(table, &QTableWidget::cellClicked, &dialog,
            [&](int row, int) { source->setCurrentIndex(row + 2); });
    fill();
    redraw();
    dialog.resize(950, 730);
    dialog.exec();
}

void ImagingDialog::showDepthReport() {
    QDialog dialog(this);
    dialog.setWindowTitle(tr("薄箔留出验证 · %1").arg(m_modelSource->currentText()));
    auto layout = new QVBoxLayout(&dialog);
    auto note = new QLabel(tr("每行留出一个薄箔位置，其余位置重新拟合标定。偏差 = 恢复峰位 − "
                              "已知位置。峰宽为单条深度剖面的半高全宽，含薄箔自身厚度（%1 mm）。")
                               .arg(real(m_depthReport["foil_thickness_mm"].toDouble())));
    note->setWordWrap(true);
    layout->addWidget(note);
    auto table = new QTableWidget;
    table->setColumnCount(9);
    table->setHorizontalHeaderLabels({tr("薄箔 X (mm)"), tr("恢复 X (mm)"), tr("偏差 (mm)"), tr("FWHM (mm)"),
                                      tr("网格间距 (mm)"), tr("晶体数"), tr("相对灵敏度"), tr("算法"),
                                      tr("状态")});
    const auto folds = m_depthReport["folds"].toArray();
    table->setRowCount(int(folds.size()));
    for (int i = 0; i < folds.size(); ++i) {
        const auto row = folds[i].toObject();
        const auto method = depthMethod(row);
        const auto peak = row[method].toObject();
        int col = 0;
        for (const auto &v : {row["foil_x_mm"], peak["peak_x_mm"], peak["bias_mm"], peak["fwhm_mm"],
                              row["dx_mm"], row["roi_count"], row["relative_sensitivity"]})
            table->setItem(i, col++,
                           cell(v.isDouble() ? QString::number(v.toDouble(), 'g', 6) : QStringLiteral("—")));
        const QString quality = depthQuality(row);
        table->setItem(i, col++, cell(row.contains("error") ? QStringLiteral("—") : method.toUpper()));
        table->setItem(i, col, cell(quality.isEmpty() ? tr("已收敛") : quality));
    }
    table->resizeColumnsToContents();
    table->horizontalHeader()->setStretchLastSection(true);
    layout->addWidget(table);
    dialog.resize(910, 490);
    dialog.exec();
}

void ImagingDialog::buildModel() {
    if (m_busy || !m_haveResult)
        return;
    bool valid;
    const double level = m_modelLevel->text().trimmed().toDouble(&valid);
    if (!valid || !(m_modelMin < level && level < m_modelMax)) {
        status(tr("阈值须介于 %1 和 %2 之间。").arg(real(m_modelMin), real(m_modelMax)), true);
        return;
    }
    const auto *volume = &volumeData(m_modelSource->currentData().toInt());
    const auto axes = m_result.axes;
    m_buildingModel = true;
    setBusy(true);
    m_cancel.store(false);
    m_progress->setRange(0, 0);
    m_progress->show();
    status(tr("正在提取三维等值面…"));
    m_watcher.setFuture(QtConcurrent::run([this, volume, axes, level] {
        Work work;
        try {
            work.mesh = Isosurface::extract(axes, volume->data(), volume->size(), level, m_cancel);
        } catch (const std::exception &e) {
            work.error = QString::fromUtf8(e.what());
        }
        return work;
    }));
}

void ImagingDialog::saveModelGif() {
    if (m_busy || !m_model->hasMesh())
        return;
    const auto title = tr("保存旋转 GIF");
    bool accepted = false;
    const double seconds = QInputDialog::getDouble(this, title, tr("每转一圈的时间（秒）"), m_rotationPeriod,
                                                   2., 30., 1, &accepted);
    if (!accepted)
        return;

    const auto path = WorkspacePaths::saveFile(
        this, Workflow::Imaging, Kind::Result, m_result.directory, title, "model-rotation.gif",
        tr("GIF 动图 (*.gif)"), QDir(m_result.directory).filePath("model-rotation.gif"));
    if (path.isEmpty())
        return;
    auto render = m_model->rotationFrames(QSize(720, 720));
    m_rotationPeriod = seconds;
    const int period = qRound(seconds * 100);
    const int frames = std::clamp(qRound(seconds * 20), 60, 180);

    m_exportingModelGif = true;
    m_cancel.store(false);
    setBusy(true);
    m_cancelButton->setText(tr("取消导出"));
    m_progress->setRange(0, frames);
    m_progress->setValue(0);
    m_progress->show();
    status(tr("正在导出旋转 GIF…"));
    m_watcher.setFuture(QtConcurrent::run(
        [this, path, render = std::move(render), period, frames](QPromise<Work> &promise) {
            promise.setProgressRange(0, frames);
            QSaveFile file(path);
            auto save = [&] {
                if (m_cancel.load() || !file.open(QIODevice::WriteOnly))
                    return false;
                auto frame = render(0.);
                AnimatedGif gif(file, frame);
                int elapsed = 0;
                for (int i = 0; i < frames; ++i) {
                    if (m_cancel.load())
                        return false;
                    if (i > 0)
                        frame = render(double(elapsed) / period);
                    // Sample the angle at the same centisecond times stored in the GIF.
                    const int next = qRound(double((i + 1) * period) / frames);
                    if (!gif.append(frame, next - elapsed))
                        return false;
                    elapsed = next;
                    promise.setProgressValue(i + 1);
                }
                return gif.finish() && !m_cancel.load() && file.commit();
            };
            Work work;
            try {
                if (save())
                    work.gifPath = path;
                else if (!m_cancel.load())
                    work.error = tr("无法保存旋转 GIF（%1）。").arg(file.errorString());
            } catch (const std::exception &e) {
                work.error = tr("无法保存旋转 GIF（%1）。").arg(QString::fromUtf8(e.what()));
            }
            promise.addResult(std::move(work));
        }));
}

void ImagingDialog::refreshTrajectory() {
    if (!m_haveResult || m_psfChoice->currentText().isEmpty()) {
        m_psfSelection->clear();
        return;
    }
    auto model = m_result.models[m_psfChoice->currentText()].toObject();
    auto cv = model["model_selection"].toObject();
    const QStringList keys{"constant_gaussian", "varying_gaussian", "measured_template"};
    const QStringList names{tr("恒定高斯"), tr("变宽度、变幅度高斯"), tr("实测模板插值")};
    const int selected = std::max(0, int(keys.indexOf(cv["selected_model"].toString())));
    QString detail = tr("采用模型：%1").arg(names[selected]);
    const auto comparison = m_psfComparison["models"].toObject().value(m_psfChoice->currentText()).toObject();
    if (!comparison.isEmpty()) {
        cv = comparison["model_selection"].toObject();
        detail += tr(" · 本次标定比较推荐：%1")
                      .arg(names[std::max(0, int(keys.indexOf(cv["selected_model"].toString())))]);
    }
    const auto scores = cv["relative_l2"].toObject();
    if (!scores.isEmpty()) {
        QStringList shown;
        for (int i = 0; i < keys.size(); ++i)
            if (scores.contains(keys[i]))
                shown << names[i] + " " + QString::number(100 * scores[keys[i]].toDouble(), 'f', 2) + "%";
        detail += tr(" · %1 个留出位置 · 整幅 ROI 相对误差：%2")
                      .arg(cv["scored_positions"].toInt())
                      .arg(shown.join(" / "));
    } else if (cv["status"] == "incomplete_folds")
        detail += tr(" · 部分留出轮次无法标定，保留基线");
    else
        detail += tr(" · 本次重建使用固定基线");
    const auto templateStatus = cv["template_status"].toString();
    if (templateStatus == "too_few_positions")
        detail += tr(" · 模板未参评：合格位置 %1/8").arg(model["valid_observations"].toInt());
    else if (templateStatus == "low_asymmetry")
        detail += tr(" · 模板未参评：不对称指标 %1，门限 0.12")
                      .arg(QString::number(cv["median_asymmetry"].toDouble(), 'f', 3));
    else if (templateStatus == "insufficient_interpolation_tests")
        detail += tr(" · 模板未参评：区间内有效留出预测少于 3 次");
    m_psfSelection->setText(detail);
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
    m_nextButton->setVisible(step + 1 < m_tabs->count());
    m_nextButton->setEnabled(!m_busy && (step != 2 || m_haveResult));
    m_validateButton->setVisible(step == 0 || step == 2);
    m_runButton->setVisible(step == 2);
    m_cancelButton->setVisible(m_busy);
    m_modelSource->setEnabled(!m_busy && m_haveResult);
    m_modelLevel->setEnabled(!m_busy && m_haveResult && m_modelMin < m_modelMax);
    m_buildModelButton->setEnabled(m_modelLevel->isEnabled());
    m_saveModelGifButton->setEnabled(!m_busy && m_model->hasMesh());
    m_validateDepth->setEnabled(!m_busy && m_haveResult &&
                                QFileInfo::exists(m_result.directory + "/config.resolved.json") &&
                                QFileInfo::exists(m_result.directory + "/tomography_inputs.h5"));
    m_comparePsf->setEnabled(m_validateDepth->isEnabled());
    m_compareBounds->setEnabled(m_validateDepth->isEnabled() &&
                                m_result.summary["accepted_roi_count"].toInt() >= 2);
    m_showBounds->setEnabled(!m_busy && m_haveResult && !m_boundaryReport.isEmpty());
    m_showDepth->setEnabled(!m_busy && !m_depthReport.isEmpty());
    m_showLow->setEnabled(!m_relativeSensitivity.isEmpty());
    m_lowPercent->setEnabled(m_showLow->isEnabled() && m_showLow->isChecked());
    m_nextButton->setDefault(step < 2);
    m_runButton->setDefault(step == 2);
}

void ImagingDialog::setBusy(bool busy) {
    m_busy = busy;
    m_inputs->setEnabled(!busy);
    m_parameters->setEnabled(!busy);
    for (auto button : {m_validateButton, m_runButton, m_importButton, m_exportButton})
        button->setEnabled(!busy);
    m_cancelButton->setEnabled(busy);
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
