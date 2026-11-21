#include "SpectrumProcessingPage.h"
#include "SpectrumPlot.h"
#include "WorkflowUi.h"
#include "WorkspacePaths.h"
#include <QJsonArray>
#include <QtConcurrent/QtConcurrentRun>
#include <QtWidgets>
#include <cmath>

namespace {
using WorkspacePaths::Kind;
using WorkspacePaths::Workflow;

const QString inputFilter = QStringLiteral("Spectra (*.h5 *.hdf5 *.nxs *.txt *.tsv);;All files (*)");
constexpr auto linearKey = "low_q_linear_windows";
constexpr auto defaultLinearWindows = "460:510";
constexpr auto totalTarget = "total";
constexpr int batchModuleRole = Qt::UserRole + 1;
const std::array<QColor, backgroundComponentCount> backgroundColors{
    QColor("#d55e00"), QColor("#996025"), QColor("#8e5ab5"), QColor("#b83e78")};
QDoubleSpinBox *numberBox(double min, double max, double value, int precision = 4) {
    auto box = new QDoubleSpinBox;
    box->setDecimals(precision);
    box->setRange(min, max);
    box->setValue(value);
    box->setKeyboardTracking(false);
    return box;
}
SpectrumPlot *spectroscopyPlot(const QString &name) {
    auto plot = new SpectrumPlot(nullptr, SpectrumPlot::Style::Spectroscopy);
    plot->setObjectName(name);
    return plot;
}
QLabel *textLabel(const QString &text = {}) {
    auto label = new QLabel(text);
    label->setWordWrap(true);
    label->setTextInteractionFlags(Qt::TextSelectableByMouse);
    label->setTextFormat(Qt::PlainText);
    return label;
}
QString formatNumber(double value) {
    return std::isfinite(value) ? QString::number(value, 'g', 7) : QStringLiteral("—");
}
struct MeasuredRange {
    double lo = NAN, hi = NAN, maxSpacing = 0., maximum = -INFINITY, peakEnergy = NAN;
    qsizetype count = 0;
};
MeasuredRange measuredRange(const QVector<double> &axis, const QVector<double> &values) {
    MeasuredRange result;
    for (qsizetype i = 0; i < std::min(axis.size(), values.size()); ++i) {
        if (!std::isfinite(values[i]))
            continue;
        if (result.count == 0)
            result.lo = axis[i];
        else
            result.maxSpacing = std::max(result.maxSpacing, axis[i] - result.hi);
        result.hi = axis[i];
        ++result.count;
        if (values[i] > result.maximum) {
            result.maximum = values[i];
            result.peakEnergy = axis[i];
        }
    }
    return result;
}
SpectrumWindows parseWindows(const QString &text, const QString &label) {
    SpectrumWindows result;
    for (const auto &part : text.split(QRegularExpression("[,，;；]+"), Qt::SkipEmptyParts)) {
        auto pair = part.split(QRegularExpression("[:：]"));
        bool loOk, hiOk;
        double lo = pair.value(0).trimmed().toDouble(&loOk), hi = pair.value(1).trimmed().toDouble(&hiOk);
        if (pair.size() != 2 || !loOk || !hiOk || !std::isfinite(lo) || !std::isfinite(hi) || hi <= lo)
            throw std::runtime_error(
                QObject::tr("%1的格式应为 起点:终点，可用逗号分隔多段。").arg(label).toUtf8().constData());
        result.append({lo, hi});
    }
    return result;
}
} // namespace

SpectrumProcessingPage::SpectrumProcessingPage(QLineEdit &hfTable, QLineEdit &composition, QLineEdit &hfShells,
                                               QWidget *parent)
    : QWidget(parent), m_hfTable(&hfTable), m_composition(&composition), m_hfShells(&hfShells) {
    auto input = new QVBoxLayout(this);
    WorkflowUi::compact(input);
    m_inputControls = new QWidget;
    auto controls = new QVBoxLayout(m_inputControls);
    controls->setContentsMargins(0, 0, 0, 0);
    auto acquisitionBox = new QGroupBox(tr("数据源"));
    auto acquisitionLayout = new QVBoxLayout(acquisitionBox);
    WorkflowUi::compact(acquisitionLayout);
    auto commands = new QFormLayout;
    m_acquisition = new QComboBox;
    m_acquisition->setObjectName("ComptonAcquisition");
    m_acquisition->addItems({QStringLiteral("Full Scan"), tr("独立测量")});
    m_loadWide = WorkflowUi::button("folder-open", tr("选择宽谱文件…"));
    m_loadWide->setObjectName("LoadComptonWideSpectrum");
    m_loadWide->hide();
    commands->addRow(tr("测量方式"), m_acquisition);
    commands->addRow(m_loadWide);
    acquisitionLayout->addLayout(commands);
    controls->addWidget(acquisitionBox);
    const SpectrumHdfLayout defaults;
    auto mapping = new QWidget;
    auto form = new QFormLayout(mapping);
    m_energyPath = new QLineEdit(defaults.energy);
    m_labelsPath = new QLineEdit(defaults.labels);
    m_intensityPath = new QLineEdit(defaults.intensity);
    m_unit = new QLineEdit(defaults.externalIntensityUnit);
    m_kev = new QCheckBox(tr("能量以 keV 存储"));
    m_energyFirst = new QCheckBox(tr("矩阵按 [能量, 晶体] 排列"));
    form->addRow(tr("能量转移数据集"), m_energyPath);
    form->addRow(tr("晶体标签数据集"), m_labelsPath);
    form->addRow(tr("强度数据集"), m_intensityPath);
    form->addRow(tr("外部文件强度单位"), m_unit);
    form->addRow(m_kev);
    form->addRow(m_energyFirst);
    acquisitionLayout->addWidget(WorkflowUi::disclosure(tr("外部 HDF5 数据集映射"), mapping));
    auto previewBox = new QGroupBox(tr("谱线与能区"));
    auto previewLayout = new QVBoxLayout(previewBox);
    WorkflowUi::compact(previewLayout);
    auto selection = new QFormLayout;
    m_inputCrystal = new QComboBox;
    m_inputCrystal->setObjectName("ComptonPeakCrystal");
    m_inputCrystal->setMinimumContentsLength(12);
    m_normalizeNear = new QCheckBox(tr("近边谱按扫描次数归一化"));
    m_normalizeNear->setObjectName("NormalizeComptonNearEdge");
    m_normalizeNear->setToolTip(tr("用第 5 步的叠加强度除以参与叠加的扫描数，并将平均强度用于拼接。"));
    m_limitPeak = new QCheckBox(tr("限定预览能区"));
    m_limitPeak->setObjectName("LimitComptonPeakRange");
    m_peakLo = numberBox(-1e8, 1e8, 0., 6);
    m_peakHi = numberBox(-1e8, 1e8, 0., 6);
    m_peakLo->setObjectName("ComptonPeakMinimum");
    m_peakHi->setObjectName("ComptonPeakMaximum");
    for (auto bound : {m_peakLo, m_peakHi})
        bound->setSuffix(" eV");
    m_splice = WorkflowUi::button("layers", tr("Compton 峰拼接"));
    m_splice->setObjectName("SpliceComptonPeak");
    selection->addRow(tr("预览谱"), m_inputCrystal);
    selection->addRow(m_normalizeNear);
    selection->addRow(m_limitPeak);
    selection->addRow(tr("能区下限"), m_peakLo);
    selection->addRow(tr("能区上限"), m_peakHi);
    previewLayout->addLayout(selection);
    controls->addWidget(previewBox);
    controls->addWidget(m_splice);
    controls->addStretch();

    auto inputSplit = new QSplitter(Qt::Horizontal);
    inputSplit->setObjectName("ComptonPeakSplitter");
    inputSplit->setChildrenCollapsible(false);
    auto inputScroll = new QScrollArea;
    inputScroll->setWidget(m_inputControls);
    inputScroll->setWidgetResizable(true);
    inputScroll->setFrameShape(QFrame::NoFrame);
    inputScroll->setMinimumWidth(300);
    inputSplit->addWidget(inputScroll);

    auto inputReview = new QWidget;
    auto inputReviewLayout = new QVBoxLayout(inputReview);
    inputReviewLayout->setContentsMargins(4, 0, 0, 0);
    m_inputInfo = textLabel();
    m_inputInfo->setObjectName("ComptonPeakCoverage");
    inputReviewLayout->addWidget(m_inputInfo);
    m_inputPlot = spectroscopyPlot("ComptonPeakSpectrum");
    inputReviewLayout->addWidget(m_inputPlot, 1);
    m_spliceInfo = textLabel();
    m_spliceInfo->setObjectName("ComptonSpliceInfo");
    inputReviewLayout->addWidget(m_spliceInfo);
    inputSplit->addWidget(inputReview);
    inputSplit->setStretchFactor(1, 1);
    inputSplit->setSizes({330, 720});
    input->addWidget(inputSplit, 1);
    m_status[0] = textLabel();
    m_status[0]->setObjectName("ComptonPeakStatus");
    input->addWidget(m_status[0]);

    m_analysisPage = new QWidget;
    auto analysis = new QVBoxLayout(m_analysisPage);
    WorkflowUi::compact(analysis);
    m_fitControls = new QWidget;
    auto fitLayout = new QVBoxLayout(m_fitControls);
    fitLayout->setContentsMargins(0, 0, 0, 0);
    auto sourceBox = new QGroupBox(tr("背景分析输入"));
    auto sourceFields = new QFormLayout(sourceBox);
    m_analysisInputInfo = textLabel();
    m_analysisInputInfo->setObjectName("BackgroundInputInfo");
    sourceFields->addRow(m_analysisInputInfo);
    m_scope = new QComboBox;
    m_scope->setObjectName("BackgroundScope");
    m_scope->addItems({tr("分晶体"), tr("模组谱"), tr("总谱")});
    sourceFields->addRow(tr("处理对象"), m_scope);
    m_fitCrystal = new QComboBox;
    m_fitCrystal->setObjectName("BackgroundDiagnosticCrystal");
    auto previewChoice = new QHBoxLayout;
    previewChoice->addWidget(m_fitCrystal, 1);
    auto previewArrow = [this, previewChoice](const QString &name, int offset, const QString &icon,
                                               const QString &label) {
        auto button = new QToolButton;
        button->setObjectName(name);
        button->setIcon(WorkflowUi::icon(icon));
        button->setIconSize({16, 16});
        button->setToolTip(label);
        button->setAccessibleName(label);
        connect(button, &QToolButton::clicked, this, [this, offset] {
            const int next = m_fitCrystal->currentIndex() + offset;
            if (next >= 0 && next < m_fitCrystal->count())
                m_fitCrystal->setCurrentIndex(next);
        });
        previewChoice->addWidget(button);
        return button;
    };
    m_previousPreview = previewArrow("PreviousBackgroundPreview", -1, "chevron-left", tr("上一条谱"));
    m_nextPreview = previewArrow("NextBackgroundPreview", 1, "chevron-right", tr("下一条谱"));
    sourceFields->addRow(tr("预览谱"), previewChoice);
    m_protected = new QLineEdit;
    m_protected->setObjectName("BackgroundProtectedWindows");
    m_protected->setPlaceholderText(tr("可选，例如 280:310"));
    m_protected->setToolTip(tr("拟合各项时跳过这里的采样点。多段用逗号分隔。"));
    sourceFields->addRow(tr("保护能区 / eV"), m_protected);
    fitLayout->addWidget(sourceBox);

    m_asymmetryRise = numberBox(.01, 100, 1., 3);
    m_asymmetryDecay = numberBox(.01, 100, 1., 3);

    for (int i = 0; i < backgroundComponentCount; ++i) {
        const auto component = BackgroundComponent(i);
        const auto key = backgroundComponentKey(component);
        auto group = new QGroupBox(componentName(component));
        auto groupLayout = new QVBoxLayout(group);
        WorkflowUi::compact(groupLayout);
        auto command = new QHBoxLayout;
        m_enabledTerms[i] = new QCheckBox(tr("启用"));
        m_enabledTerms[i]->setObjectName("EnableBackground_" + key);
        m_enabledTerms[i]->setChecked(component == BackgroundComponent::Linear);
        m_termStatus[i] = textLabel();
        m_fitTerms[i] = WorkflowUi::button("chart-line", tr("拟合"));
        m_fitTerms[i]->setObjectName("FitBackground_" + key);
        command->addWidget(m_enabledTerms[i]);
        command->addWidget(m_termStatus[i], 1);
        command->addWidget(m_fitTerms[i]);
        groupLayout->addLayout(command);

        m_termControls[i] = new QWidget;
        auto fields = new QFormLayout(m_termControls[i]);
        WorkflowUi::compact(fields);
        m_windows[i] = new QLineEdit;
        m_windows[i]->setObjectName("BackgroundWindows_" + key);
        m_windows[i]->setPlaceholderText(tr("起点:终点，多段用逗号分隔"));
        if (component == BackgroundComponent::Linear) {
            m_windows[i]->setObjectName("LowMomentumFitWindows");
            m_windows[i]->setText(defaultLinearWindows);
        }
        fields->addRow(tr("拟合能区 / eV"), m_windows[i]);
        if (component == BackgroundComponent::HfCore) {
            fields->addRow(textLabel(tr("在第 1 步设置元素与化学比、HF 壳层表和当前分组的扣除壳层。")));
        } else if (component == BackgroundComponent::Asymmetry) {
            fields->addRow(tr("上升宽度 / a.u."), m_asymmetryRise);
            fields->addRow(tr("衰减宽度 / a.u."), m_asymmetryDecay);
        }
        groupLayout->addWidget(m_termControls[i]);
        fitLayout->addWidget(group);
        auto changed = [this, i] {
            resetAnalysis(i);
            emit analysisParametersEdited();
        };
        connect(m_enabledTerms[i], &QCheckBox::toggled, this, changed);
        connect(m_windows[i], &QLineEdit::textChanged, this, changed);
        connect(m_fitTerms[i], &QPushButton::clicked, this, [this, component] { fitOne(component); });
    }
    m_batchTargets = new QListWidget;
    m_batchTargets->setObjectName("BackgroundBatchTargets");
    m_batchTargets->setMaximumHeight(144);
    auto batchControls = new QWidget;
    auto batchLayout = new QVBoxLayout(batchControls);
    batchLayout->setContentsMargins(0, 0, 0, 0);
    auto batchActions = new QHBoxLayout;
    m_batchModule = new QPushButton(tr("按模组选…"));
    m_batchModule->setObjectName("SelectBackgroundBatchModules");
    m_batchAll = new QPushButton(tr("全选"));
    m_batchAll->setObjectName("SelectAllBackgroundTargets");
    m_batchInvert = new QPushButton(tr("反选"));
    m_batchInvert->setObjectName("InvertBackgroundTargets");
    for (auto button : {m_batchModule, m_batchAll, m_batchInvert}) {
        button->setAutoDefault(false);
        batchActions->addWidget(button);
    }
    batchLayout->addLayout(batchActions);
    batchLayout->addWidget(m_batchTargets);
    fitLayout->addWidget(WorkflowUi::disclosure(tr("批量目标"), batchControls));
    m_run = WorkflowUi::button("filter", tr("批量扣除"));
    m_run->setObjectName("RunMomentumBackground");
    fitLayout->addStretch();
    auto scroll = new QScrollArea;
    scroll->setWidget(m_fitControls);
    scroll->setWidgetResizable(true);
    scroll->setMinimumWidth(340);
    scroll->setFrameShape(QFrame::NoFrame);
    auto backgroundSplit = new QSplitter(Qt::Horizontal);
    backgroundSplit->setChildrenCollapsible(false);
    backgroundSplit->addWidget(scroll);
    auto preview = new QWidget;
    auto review = new QVBoxLayout(preview);
    review->setContentsMargins(4, 0, 0, 0);
    m_fitPlot = spectroscopyPlot("MomentumBackgroundDiagnostics");
    review->addWidget(m_fitPlot, 1);
    m_fitInfo = textLabel();
    review->addWidget(m_fitInfo);
    backgroundSplit->addWidget(preview);
    backgroundSplit->setStretchFactor(1, 1);
    backgroundSplit->setSizes({360, 700});
    analysis->addWidget(backgroundSplit, 1);
    auto runStatus = new QHBoxLayout;
    runStatus->addWidget(m_status[1] = textLabel(), 1);
    m_cancelButton = WorkflowUi::button("square", tr("取消分析"));
    m_cancelButton->hide();
    runStatus->addWidget(m_cancelButton);
    runStatus->addWidget(m_run);
    analysis->addLayout(runStatus);

    m_resultPage = new QWidget;
    auto results = new QVBoxLayout(m_resultPage);
    WorkflowUi::compact(results);
    m_resultControls = new QWidget;
    auto resultControls = new QVBoxLayout(m_resultControls);
    resultControls->setContentsMargins(0, 0, 0, 0);
    auto resultOptions = new QGroupBox(tr("谱图设置"));
    auto resultFields = new QFormLayout(resultOptions);
    WorkflowUi::compact(resultFields);
    m_resultChoice = new QComboBox;
    m_resultChoice->setObjectName("BackgroundResultSpectrum");
    m_resultView = new QComboBox;
    m_resultView->setObjectName("BackgroundResultView");
    m_resultView->addItems({tr("扣除后"), tr("原谱与背景分量")});
    m_resultLo = numberBox(-1e8, 1e8, 520, 6);
    m_resultHi = numberBox(-1e8, 1e8, 600, 6);
    m_resultLo->setObjectName("CorrectedEnergyMinimum");
    m_resultHi->setObjectName("CorrectedEnergyMaximum");
    for (auto bound : {m_resultLo, m_resultHi})
        bound->setSuffix(" eV");
    resultFields->addRow(tr("结果谱"), m_resultChoice);
    resultFields->addRow(tr("显示内容"), m_resultView);
    resultFields->addRow(tr("能区下限"), m_resultLo);
    resultFields->addRow(tr("能区上限"), m_resultHi);
    resultControls->addWidget(resultOptions);

    auto exportOptions = new QGroupBox(tr("保存当前能区"));
    auto exportLayout = new QVBoxLayout(exportOptions);
    WorkflowUi::compact(exportLayout);
    m_saveResult = WorkflowUi::button("save", tr("保存为 HDF5…"));
    m_saveResult->setObjectName("SaveNearEdgeHdf5");
    m_saveText = WorkflowUi::button("save", tr("保存为文本…"));
    m_saveText->setObjectName("SaveNearEdgeText");
    exportLayout->addWidget(m_saveResult);
    exportLayout->addWidget(m_saveText);
    resultControls->addWidget(exportOptions);
    resultControls->addStretch();

    auto resultSplit = new QSplitter(Qt::Horizontal);
    resultSplit->setObjectName("NearEdgeResultSplitter");
    resultSplit->setChildrenCollapsible(false);
    m_resultControls->setMinimumWidth(280);
    resultSplit->addWidget(m_resultControls);
    auto resultPreview = new QWidget;
    auto resultPreviewLayout = new QVBoxLayout(resultPreview);
    resultPreviewLayout->setContentsMargins(4, 0, 0, 0);
    m_resultPlot = spectroscopyPlot("CorrectedCrystalSpectrum");
    resultPreviewLayout->addWidget(m_resultPlot, 1);
    resultPreviewLayout->addWidget(m_resultInfo = textLabel());
    resultSplit->addWidget(resultPreview);
    resultSplit->setStretchFactor(1, 1);
    resultSplit->setSizes({300, 760});
    results->addWidget(resultSplit, 1);
    results->addWidget(m_status[2] = textLabel());

    auto spliceSettingsChanged = [this] {
        resetSplice();
        emit analysisParametersEdited();
    };
    connect(m_acquisition, &QComboBox::currentIndexChanged, this, spliceSettingsChanged);
    auto wideMappingChanged = [this] {
        ++m_wideRevision;
        m_wide = {};
        resetSplice();
        emit analysisParametersEdited();
    };
    for (auto edit : {m_energyPath, m_labelsPath, m_intensityPath, m_unit})
        connect(edit, &QLineEdit::textChanged, this, wideMappingChanged);
    for (auto option : {m_kev, m_energyFirst})
        connect(option, &QCheckBox::toggled, this, wideMappingChanged);
    connect(m_loadWide, &QPushButton::clicked, this, [this] {
        auto suggestion = m_wideFile;
        if (suggestion.isEmpty())
            suggestion = QDir(WorkspacePaths::directory(Workflow::Spectroscopy, Kind::Result)).filePath("spectra.h5");
        const auto path = WorkflowUi::selectSource(this, m_availableFiles, suggestion,
                                                   tr("选择 Compton 宽谱文件"), inputFilter);
        if (!path.isEmpty())
            loadWide(QFileInfo(path).absoluteFilePath());
    });
    connect(m_limitPeak, &QCheckBox::toggled, this, &SpectrumProcessingPage::refreshInput);
    connect(m_normalizeNear, &QCheckBox::toggled, this, spliceSettingsChanged);
    for (auto bound : {m_peakLo, m_peakHi})
        connect(bound, &QDoubleSpinBox::valueChanged, this, &SpectrumProcessingPage::refreshInput);
    connect(m_splice, &QPushButton::clicked, this, [this] { startSplice(); });
    connect(m_scope, &QComboBox::currentIndexChanged, this, [this] {
        m_configuredBatchTargets.reset();
        refreshAnalysisInput();
        emit analysisParametersEdited();
    });
    connect(m_run, &QPushButton::clicked, this, &SpectrumProcessingPage::runBatch);
    connect(m_batchTargets, &QListWidget::itemChanged, this, [this] {
        m_configuredBatchTargets = batchTargets();
        resetBatch();
        refreshButtons();
        emit analysisParametersEdited();
    });
    connect(m_batchModule, &QPushButton::clicked, this, &SpectrumProcessingPage::chooseBatchModules);
    connect(m_batchAll, &QPushButton::clicked, this, [this] {
        updateBatchTargets([](const auto &) { return true; });
    });
    connect(m_batchInvert, &QPushButton::clicked, this, [this] {
        updateBatchTargets([](const auto &item) { return item.checkState() != Qt::Checked; });
    });
    connect(m_cancelButton, &QPushButton::clicked, this, [this] {
        m_cancel.store(true);
        message(tr("正在取消…"));
    });
    connect(m_inputCrystal, &QComboBox::currentIndexChanged, this, &SpectrumProcessingPage::refreshInput);
    connect(m_fitCrystal, &QComboBox::currentIndexChanged, this, [this] {
        refreshAnalysis();
        refreshButtons();
    });
    connect(m_protected, &QLineEdit::textChanged, this, [this] {
        resetAnalysis();
        emit analysisParametersEdited();
    });
    for (auto edit : {m_hfTable, m_composition, m_hfShells})
        connect(edit, &QLineEdit::textChanged, this, [this] {
            resetAnalysis(int(BackgroundComponent::HfCore));
            emit analysisParametersEdited();
        });
    for (auto spin : {m_asymmetryRise, m_asymmetryDecay})
        connect(spin, &QDoubleSpinBox::valueChanged, this, [this] {
            resetAnalysis(int(BackgroundComponent::Asymmetry));
            emit analysisParametersEdited();
        });
    connect(m_resultChoice, &QComboBox::currentIndexChanged, this, &SpectrumProcessingPage::refreshResult);
    connect(m_resultView, &QComboBox::currentIndexChanged, this, &SpectrumProcessingPage::refreshResult);
    for (auto spin : {m_resultLo, m_resultHi})
        connect(spin, &QDoubleSpinBox::valueChanged, this, [this] { refreshResult(); });
    connect(m_saveResult, &QPushButton::clicked, this, [this] { saveResult(ResultFormat::Hdf5); });
    connect(m_saveText, &QPushButton::clicked, this, [this] { saveResult(ResultFormat::Text); });
    connect(&m_watcher, &QFutureWatcher<Work>::finished, this, [this] {
        auto result = m_watcher.result();
        m_busy = false;
        if (result.wide && !result.wide->crystalEnergyEv.isEmpty() && result.wideRevision != m_wideRevision) {
            result.wide.reset();
            result.error = true;
            result.message = tr("ROI 或谱学配置已更新，请重新读取原始宽谱。");
        }
        if (result.wide) {
            const bool newSource = m_wideFile != result.wideFile;
            m_wideFile = result.wideFile;
            m_wide = std::move(*result.wide);
            if (!m_wide.wideCrystals.isEmpty()) {
                m_wide.crystals = std::move(m_wide.wideCrystals);
                m_wide.pointSources.clear();
                m_wide.sources = {m_wide.sources.value(0)};
            }
            resetSplice();
            if (newSource)
                emit analysisParametersEdited();
        }
        bool inputReady = false;
        if (result.stitched) {
            if (result.stitchRevision == m_stitchRevision) {
                m_stitched = std::move(*result.stitched);
                refreshInput();
                if (result.enterAnalysis)
                    inputReady = prepareAnalysisInput();
            } else {
                result.message = tr("输入谱已改变，请重新拼接。");
            }
        }
        if (m_cancel.load() || ((result.analysis || result.term) && result.analysisRevision != m_analysisRevision)) {
            result.analysis.reset();
            result.term.reset();
            result.message = m_cancel.load() ? tr("已取消。") : tr("输入谱已更新，请重新拟合。");
            result.error = false;
        }
        if (result.term) {
            auto &preview = m_previews[result.target];
            result.term->corrected.clear();
            preview.terms.insert(result.component, std::move(*result.term));
            preview.total = sumBackground(m_targets.constFind(result.target).value(), preview.terms);
        }
        if (result.analysis) {
            m_analysis = std::move(*result.analysis);
            const QSignalBlocker block(m_resultChoice);
            m_resultChoice->clear();
            if (m_scope->currentIndex() == PerCrystal && m_analysis.spectra.size() > 1) {
                m_resultChoice->addItem(tr("全部晶体等权平均"), m_analysis.spectra.keys());
                QMap<QString, QStringList> groups;
                for (const auto &label : m_analysis.spectra.keys()) {
                    const double q = m_analysisInput.geometry.value(label).qMean;
                    if (!std::isfinite(q) || q <= 0)
                        continue;
                    const double lower = .5 * std::floor(q * 2);
                    groups[tr("%1 · q %2–%3 Å⁻¹").arg(label.section('-', 0, 0))
                               .arg(lower, 0, 'f', 1).arg(lower + .5, 0, 'f', 1)].append(label);
                }
                for (auto it = groups.cbegin(); it != groups.cend(); ++it)
                    if (it->size() > 1)
                        m_resultChoice->addItem(it.key(), it.value());
            }
            for (const auto &label : m_analysis.spectra.keys())
                m_resultChoice->addItem(targetName(label), QStringList{label});
            refreshResult();
        }
        if (result.term || result.analysis)
            refreshAnalysis();
        message(result.message.isEmpty() ? tr("已取消。") : result.message, result.error);
        m_inputControls->setEnabled(true);
        m_fitControls->setEnabled(true);
        m_resultControls->setEnabled(true);
        m_cancelButton->hide();
        refreshButtons();
        emit busyChanged();
        if (inputReady)
            emit analysisInputReady();
        if (result.analysis)
            emit analysisBatchCompleted(!m_analysis.spectra.isEmpty(), !m_analysis.failures.isEmpty());
        if (!result.exportPath.isEmpty())
            emit resultSaved(result.exportFormat, result.exportPath);
        emit analysisFinished(result.analysis.has_value() && !m_analysis.spectra.isEmpty());
    });
    refreshInput();
    refreshAnalysisInput();
}

SpectrumProcessingPage::~SpectrumProcessingPage() {
    m_cancel.store(true);
    m_watcher.waitForFinished();
}

QJsonObject SpectrumProcessingPage::backgroundConfigurationDefaults() {
    QJsonObject defaults{{"background_scope", int(PerCrystal)}, {"background_protected_ev", ""},
                         {"hf_table", ""}, {"hf_composition", ""},
                         {"background_hf_shells", ""}, {"background_asymmetry_rise", 1.},
                         {"background_asymmetry_decay", 1.}, {"compton_acquisition", int(Continuous)},
                         {"compton_limit_preview", false}, {"compton_preview_min_ev", 0.},
                         {"compton_preview_max_ev", 0.}, {"background_result_view", 0}};
    for (int i = 0; i < backgroundComponentCount; ++i) {
        const auto component = BackgroundComponent(i);
        const auto key = backgroundComponentKey(component);
        defaults["background_enable_" + key] = component == BackgroundComponent::Linear;
        if (component != BackgroundComponent::Linear)
            defaults["background_windows_" + key] = "";
    }
    return defaults;
}

QMap<QString, QWidget *> SpectrumProcessingPage::configurationFields() const {
    QMap<QString, QWidget *> fields{{"external_energy_dataset", m_energyPath},
                                  {"external_labels_dataset", m_labelsPath},
                                  {"external_intensity_dataset", m_intensityPath},
                                  {"external_intensity_unit", m_unit},
                                  {"external_energy_kev", m_kev},
                                  {"external_energy_first", m_energyFirst},
                                  {"normalize_near_edge_by_scan_count", m_normalizeNear},
                                  {"compton_acquisition", m_acquisition},
                                  {"compton_limit_preview", m_limitPeak},
                                  {"compton_preview_min_ev", m_peakLo},
                                  {"compton_preview_max_ev", m_peakHi},
                                  {"background_scope", m_scope},
                                  {"background_protected_ev", m_protected},
                                  {"background_hf_shells", m_hfShells},
                                  {"background_asymmetry_rise", m_asymmetryRise},
                                  {"background_asymmetry_decay", m_asymmetryDecay},
                                  {"hf_table", m_hfTable},
                                  {"hf_composition", m_composition},
                                  {"background_result_view", m_resultView},
                                  {"result_energy_min_ev", m_resultLo},
                                  {"result_energy_max_ev", m_resultHi}};
    for (int i = 0; i < backgroundComponentCount; ++i) {
        const auto component = BackgroundComponent(i);
        const auto key = backgroundComponentKey(component);
        fields["background_enable_" + key] = m_enabledTerms[i];
        fields[component == BackgroundComponent::Linear ? QString(linearKey) : "background_windows_" + key] =
            m_windows[i];
    }
    return fields;
}

QJsonValue SpectrumProcessingPage::batchSelection() const {
    return m_configuredBatchTargets.has_value() ? QJsonValue(QJsonArray::fromStringList(*m_configuredBatchTargets))
                                                : QJsonValue(QJsonValue::Null);
}

void SpectrumProcessingPage::restoreEdgeInputs(const QString &wideFile, const QJsonValue &batchTargets) {
    m_wideFile = wideFile;
    if (batchTargets.isArray()) {
        QStringList selected;
        for (const auto &target : batchTargets.toArray())
            selected.append(target.toString());
        m_configuredBatchTargets = selected;
    }
    refreshInput();
}

void SpectrumProcessingPage::resetInputs() {
    ++m_wideRevision;
    m_current = {};
    m_wide = {};
    m_wideFile.clear();
    m_configuredBatchTargets.reset();
    m_acquisition->setCurrentIndex(Continuous);
    m_limitPeak->setChecked(false);
    m_normalizeNear->setChecked(false);
    resetSplice();
    message({});
}

void SpectrumProcessingPage::restoreWideIfNeeded() {
    if (!m_busy && m_acquisition->currentIndex() == Separate && m_wide.crystals.isEmpty() &&
        !m_wideFile.isEmpty() && !m_current.crystals.isEmpty())
        loadWide(m_wideFile);
}

void SpectrumProcessingPage::loadWide(const QString &path) {
    const auto mapping = hdfLayout();
    std::optional<WideScanCalibration> calibration;
    try {
        if (m_getWideCalibration)
            calibration = m_getWideCalibration();
    } catch (const std::exception &error) {
        message(QString::fromUtf8(error.what()), true);
        return;
    }
    start(tr("正在载入宽谱…"),
          [this, path, mapping, calibration, revision = m_wideRevision] {
              Work result;
              result.wide = SpectrumProcessing::loadWideReference(
                  path, mapping, calibration ? &*calibration : nullptr, m_cancel);
              result.wideRevision = revision;
              result.wideFile = path;
              result.message = tr("已载入宽谱 %1。").arg(QFileInfo(path).fileName());
              return result;
          }, false);
}

void SpectrumProcessingPage::setCurrentResult(ProcessedSpectrum result) {
    if (result.crystals.isEmpty()) {
        ++m_wideRevision;
        if (!m_wide.crystalEnergyEv.isEmpty())
            m_wide = {};
    }
    m_current = std::move(result);
    resetSplice();
}

bool SpectrumProcessingPage::prepareAnalysisInput() {
    if (m_busy)
        return false;
    if (!m_analysisInput.crystals.isEmpty())
        return true;
    if (m_current.crystals.isEmpty()) {
        message(tr("请先在第 5 步生成原始谱。"), true);
        return false;
    }
    const bool separate = m_acquisition->currentIndex() == Separate;
    if (separate && m_stitched.crystals.isEmpty()) {
        if (m_wide.crystals.isEmpty()) {
            if (m_wideFile.isEmpty())
                message(tr("请选择 Compton 宽谱文件。"), true);
            else
                restoreWideIfNeeded();
        } else {
            startSplice(true);
        }
        return false;
    }
    m_analysisInput = separate ? m_stitched : m_current;
    refreshAnalysisInput();
    return true;
}

bool SpectrumProcessingPage::ensureAnalysis() {
    if (m_busy)
        return false;
    if (!m_analysis.spectra.isEmpty()) {
        refreshResult();
        return true;
    }
    if (m_run->isEnabled())
        runBatch();
    else if (m_targets.isEmpty())
        message(m_analysisInputInfo->text(), true);
    else if (!previewComplete())
        message(tr("请先逐项拟合当前预览谱中已启用的背景。"), true);
    else
        message(tr("请勾选要批量扣除的谱。"), true);
    return false;
}

const ProcessedSpectrum &SpectrumProcessingPage::peakInput() const {
    return m_acquisition->currentIndex() == Separate ? m_wide : m_current;
}

QStringList SpectrumProcessingPage::previewCrystals(QStringList labels) const {
    const auto choice = m_inputCrystal->currentData().toString();
    if (choice.startsWith("crystal:"))
        return {choice.mid(8)};
    if (choice.startsWith("module:")) {
        const auto prefix = choice.mid(7) + '-';
        labels.removeIf([&](const QString &label) { return !label.startsWith(prefix); });
    }
    return labels;
}

void SpectrumProcessingPage::refreshPreviewChoices(const QStringList &labels) {
    QList<QPair<QString, QString>> entries;
    if (!labels.isEmpty()) {
        entries.append(qMakePair(tr("总谱"), QStringLiteral("total")));
        for (const auto *module : {"VB", "HL", "VU", "HB", "VD"}) {
            const auto name = QString::fromLatin1(module);
            const auto prefix = name + '-';
            if (std::any_of(labels.cbegin(), labels.cend(),
                            [&](const QString &label) { return label.startsWith(prefix); }))
                entries.append(qMakePair(tr("%1 模组谱").arg(name), "module:" + name));
        }
        for (const auto &label : labels)
            entries.append(qMakePair(label, "crystal:" + label));
    }
    QStringList shown, wanted;
    for (int i = 0; i < m_inputCrystal->count(); ++i)
        shown.append(m_inputCrystal->itemData(i).toString());
    for (const auto &entry : entries)
        wanted.append(entry.second);
    if (shown == wanted)
        return;
    const auto choice = m_inputCrystal->currentData();
    const QSignalBlocker block(m_inputCrystal);
    m_inputCrystal->clear();
    for (const auto &entry : entries)
        m_inputCrystal->addItem(entry.first, entry.second);
    if (const auto i = m_inputCrystal->findData(choice); i >= 0)
        m_inputCrystal->setCurrentIndex(i);
}

void SpectrumProcessingPage::resetSplice() {
    ++m_stitchRevision;
    m_stitched = {};
    m_analysisInput = {};
    m_status[0]->clear();
    refreshInput();
    refreshAnalysisInput();
}

void SpectrumProcessingPage::startSplice(bool enterAnalysis) {
    if (m_busy || m_wide.crystals.isEmpty() || m_current.crystals.isEmpty())
        return;
    start(tr("正在拼接 Compton 宽谱与近边谱…"),
          [this, wide = m_wide, near = m_current, average = m_normalizeNear->isChecked(),
           revision = m_stitchRevision, enterAnalysis] {
              Work result;
              result.stitched = SpectrumProcessing::spliceCompton(wide, near, average, m_cancel);
              result.stitchRevision = revision;
              result.enterAnalysis = enterAnalysis;
              result.message = tr("已完成 Compton 峰拼接。");
              return result;
          }, false);
}

SpectrumHdfLayout SpectrumProcessingPage::hdfLayout() const {
    return {m_energyPath->text().trimmed(),
            m_labelsPath->text().trimmed(),
            m_intensityPath->text().trimmed(),
            m_unit->text().trimmed(),
            m_kev->isChecked(),
            m_energyFirst->isChecked()};
}

QString SpectrumProcessingPage::componentName(BackgroundComponent component) const {
    switch (component) {
    case BackgroundComponent::Linear: return tr("线性背景");
    case BackgroundComponent::HfCore: return tr("HF-core 拟合");
    case BackgroundComponent::ComptonPeak: return tr("Compton 峰拟合背景");
    case BackgroundComponent::Asymmetry: return tr("非对称函数修正");
    case BackgroundComponent::Count: return {};
    }
    return {};
}

QString SpectrumProcessingPage::targetName(const QString &key) const {
    if (m_scope->currentIndex() == TotalSpectrum)
        return tr("总谱");
    if (m_scope->currentIndex() == ByModule)
        return tr("%1 模组谱").arg(key);
    return key;
}

bool SpectrumProcessingPage::componentAvailable(BackgroundComponent component) const {
    return m_scope->currentIndex() == PerCrystal ||
           component == BackgroundComponent::Linear || component == BackgroundComponent::ComptonPeak;
}

bool SpectrumProcessingPage::canFit(BackgroundComponent component) const {
    const int index = int(component);
    if (m_busy || m_fitCrystal->currentIndex() < 0 || !componentAvailable(component) ||
        !m_enabledTerms[index]->isChecked())
        return false;
    const auto preview = m_previews.constFind(m_fitCrystal->currentData().toString());
    if (preview != m_previews.cend() && preview->terms.contains(component))
        return false;
    for (int i = 0; i < index; ++i)
        if (componentAvailable(BackgroundComponent(i)) && m_enabledTerms[i]->isChecked() &&
            (preview == m_previews.cend() || !preview->terms.contains(BackgroundComponent(i))))
            return false;
    return true;
}

bool SpectrumProcessingPage::previewComplete() const {
    const auto preview = m_previews.constFind(m_fitCrystal->currentData().toString());
    if (preview == m_previews.cend() || preview->terms.isEmpty())
        return false;
    for (int i = 0; i < backgroundComponentCount; ++i)
        if (componentAvailable(BackgroundComponent(i)) && m_enabledTerms[i]->isChecked() &&
            !preview->terms.contains(BackgroundComponent(i)))
            return false;
    return true;
}

QStringList SpectrumProcessingPage::batchTargets() const {
    QStringList selected;
    for (int i = 0; i < m_batchTargets->count(); ++i)
        if (const auto item = m_batchTargets->item(i); item->checkState() == Qt::Checked)
            selected.append(item->data(Qt::UserRole).toString());
    return selected;
}

void SpectrumProcessingPage::updateBatchTargets(const std::function<bool(const QListWidgetItem &)> &selected) {
    bool changed = false;
    const QSignalBlocker block(m_batchTargets);
    for (int i = 0; i < m_batchTargets->count(); ++i) {
        auto item = m_batchTargets->item(i);
        const auto state = selected(*item) ? Qt::Checked : Qt::Unchecked;
        if (item->checkState() != state) {
            item->setCheckState(state);
            changed = true;
        }
    }
    if (changed) {
        m_configuredBatchTargets = batchTargets();
        resetBatch();
        refreshButtons();
        emit analysisParametersEdited();
    }
}

void SpectrumProcessingPage::chooseBatchModules() {
    QMap<QString, int> modules;
    for (int i = 0; i < m_batchTargets->count(); ++i) {
        const auto module = m_batchTargets->item(i)->data(batchModuleRole).toString();
        if (!module.isEmpty())
            ++modules[module];
    }
    if (modules.isEmpty())
        return;

    QDialog dialog(m_analysisPage);
    dialog.setWindowTitle(tr("按模组选择批量目标"));
    auto layout = new QVBoxLayout(&dialog);
    QList<QCheckBox *> choices;
    for (auto it = modules.cbegin(); it != modules.cend(); ++it) {
        auto check = new QCheckBox(tr("%1 · %2 条谱").arg(it.key()).arg(it.value()));
        check->setProperty("module", it.key());
        choices.append(check);
        layout->addWidget(check);
    }
    auto buttons = new QDialogButtonBox(QDialogButtonBox::Ok | QDialogButtonBox::Cancel);
    layout->addWidget(buttons);
    connect(buttons, &QDialogButtonBox::accepted, &dialog, &QDialog::accept);
    connect(buttons, &QDialogButtonBox::rejected, &dialog, &QDialog::reject);
    if (dialog.exec() != QDialog::Accepted)
        return;
    QSet<QString> selected;
    for (auto choice : choices)
        if (choice->isChecked())
            selected.insert(choice->property("module").toString());
    updateBatchTargets([&](const auto &item) { return selected.contains(item.data(batchModuleRole).toString()); });
}

StagedBackgroundSettings SpectrumProcessingPage::settings(int last) const {
    StagedBackgroundSettings out;
    for (int i = 0; i <= last; ++i) {
        const auto component = BackgroundComponent(i);
        if (!componentAvailable(component) || !m_enabledTerms[i]->isChecked())
            continue;
        out.windows[i] = parseWindows(m_windows[i]->text(), componentName(component));
        if (out.windows[i].isEmpty())
            throw std::runtime_error(tr("请填写「%1」的拟合能区。").arg(componentName(component)).toUtf8().constData());
    }
    out.protectedWindows = parseWindows(m_protected->text(), tr("保护能区"));
    out.hfTable = effectiveHfTable(m_hfTable->text());
    out.composition = m_composition->text().trimmed();
    out.hfShells = backgroundHfShellSelection(m_hfShells->text());
    out.asymmetryRise = m_asymmetryRise->value();
    out.asymmetryDecay = m_asymmetryDecay->value();
    return out;
}

void SpectrumProcessingPage::fitOne(BackgroundComponent component) {
    if (!canFit(component))
        return;
    try {
        const auto parameters = settings(int(component));
        const auto target = m_fitCrystal->currentData().toString();
        const auto preview = m_previews.constFind(target);
        const auto residual = preview == m_previews.cend() ? m_targets.value(target) : preview->total.corrected;
        const auto geometry = m_analysisInput.geometry.value(target);
        const auto name = componentName(component);
        start(tr("正在拟合 %1…").arg(name),
              [this, parameters, target, residual, geometry, name, component, energy = m_analysisInput.energyEv,
               revision = m_analysisRevision] {
                  Work result;
                  result.term = fitBackgroundComponent(component, parameters, energy, residual, geometry, m_cancel);
                  result.target = target;
                  result.component = component;
                  result.analysisRevision = revision;
                  result.message = tr("已拟合 %1。").arg(name);
                  return result;
              });
    } catch (const std::exception &ex) {
        message(QString::fromUtf8(ex.what()), true);
    }
}

void SpectrumProcessingPage::runBatch() {
    if (!m_run->isEnabled())
        return;
    try {
        const auto parameters = settings();
        const auto selected = batchTargets();
        const auto scope = m_scope->currentIndex();
        QJsonObject metadata{{"background_scope", scope},
                             {"fit_hf_table", QFileInfo(parameters.hfTable).fileName()},
                             {"fit_hf_composition", parameters.composition},
                             {"fit_hf_shells", QJsonArray::fromStringList(parameters.hfShells)},
                             {"asymmetry_rise", parameters.asymmetryRise},
                             {"asymmetry_decay", parameters.asymmetryDecay}};
        auto windows = [](const SpectrumWindows &ranges) {
            QJsonArray out;
            for (const auto &[lo, hi] : ranges)
                out.append(QJsonArray{lo, hi});
            return out;
        };
        metadata["protected_energy_windows_ev"] = windows(parameters.protectedWindows);
        QJsonArray order;
        for (int i = 0; i < backgroundComponentCount; ++i)
            if (!parameters.windows[i].isEmpty()) {
                const auto key = backgroundComponentKey(BackgroundComponent(i));
                order.append(key);
                metadata[key + "_fit_windows_ev"] = windows(parameters.windows[i]);
            }
        metadata["background_order"] = order;
        start(tr("正在批量拟合并扣除背景…"),
              [this, parameters, selected, metadata, rows = m_targets, target = m_analysisInput,
               revision = m_analysisRevision] {
                  Work result;
                  StagedBackgroundResult batch;
                  batch.parameters = metadata;
                  for (const auto &label : selected) {
                      if (m_cancel.load())
                          break;
                      try {
                          batch.spectra.insert(label, fitBackgroundStages(parameters, target.energyEv, rows.value(label),
                                                                         target.geometry.value(label), m_cancel));
                      } catch (const std::exception &ex) {
                          batch.failures.insert(label, QString::fromUtf8(ex.what()));
                      }
                  }
                  result.message = tr("批量扣除完成 · 成功 %1 条 · 失败 %2 条")
                                       .arg(batch.spectra.size()).arg(batch.failures.size());
                  result.error = !batch.failures.isEmpty();
                  result.analysis = std::move(batch);
                  result.analysisRevision = revision;
                  return result;
              });
    } catch (const std::exception &ex) {
        message(QString::fromUtf8(ex.what()), true);
    }
}

void SpectrumProcessingPage::message(const QString &text, bool error) {
    for (auto label : m_status) {
        label->setText(text);
        label->setStyleSheet(error ? "color:#a02020" : "");
    }
}

void SpectrumProcessingPage::start(const QString &text, std::function<Work()> job, bool cancellable) {
    m_busy = true;
    emit busyChanged();
    for (auto controls : {m_inputControls, m_fitControls, m_resultControls})
        controls->setEnabled(false);
    refreshButtons();
    m_cancel.store(false);
    m_cancelButton->setVisible(cancellable);
    message(text);
    m_watcher.setFuture(QtConcurrent::run([job = std::move(job)] {
        try {
            return job();
        } catch (const std::exception &ex) {
            Work failed;
            failed.error = true;
            failed.message = QString::fromUtf8(ex.what());
            return failed;
        }
    }));
}

void SpectrumProcessingPage::resetBatch() {
    ++m_analysisRevision;
    m_analysis = {};
    m_resultChoice->clear();
    refreshResult();
    refreshAnalysis();
}

void SpectrumProcessingPage::resetAnalysis(int first) {
    if (first == 0) {
        m_previews.clear();
    } else {
        for (auto it = m_previews.begin(); it != m_previews.end(); ++it) {
            for (int i = first; i < backgroundComponentCount; ++i)
                it->terms.remove(BackgroundComponent(i));
            it->total = sumBackground(m_targets.constFind(it.key()).value(), it->terms);
        }
    }
    m_status[1]->clear();
    m_status[2]->clear();
    resetBatch();
    refreshButtons();
}

void SpectrumProcessingPage::refreshInput() {
    const auto &wide = peakInput();
    const bool separate = m_acquisition->currentIndex() == Separate;
    const bool stitched = separate && !m_stitched.crystals.isEmpty();
    const auto nearScanCount = m_current.sources.size();
    const auto names = stitched ? m_stitched.crystals.keys() : wide.crystals.keys();
    const bool limited = m_limitPeak->isChecked();
    refreshPreviewChoices(names);
    m_inputCrystal->setEnabled(!names.isEmpty());
    m_loadWide->setVisible(separate);
    m_loadWide->setToolTip(m_wideFile);
    m_splice->setVisible(separate);
    m_normalizeNear->setVisible(separate);
    m_normalizeNear->setEnabled(nearScanCount > 1);
    m_limitPeak->setEnabled(!names.isEmpty());
    m_peakLo->setEnabled(limited && !names.isEmpty());
    m_peakHi->setEnabled(limited && !names.isEmpty());
    m_spliceInfo->setVisible(stitched);
    if (stitched)
        m_spliceInfo->setText(tr("拼接结果 · %1 个同名晶体 · %2 个实测能量坐标")
                                  .arg(m_stitched.crystals.size()).arg(m_stitched.energyEv.size()));
    if (names.isEmpty()) {
        m_inputPlot->setCurves({}, {});
        m_inputInfo->setText(separate ? tr("第 5 步原始谱作为近边谱。可选择原始扫描或已处理的逐晶体宽谱。")
                                     : tr("请先在第 5 步生成连续测量的逐晶体谱。"));
        refreshButtons();
        return;
    }
    const auto label = m_inputCrystal->currentText();
    const auto labels = previewCrystals(names);
    const double lo = limited ? m_peakLo->value() : -INFINITY;
    const double hi = limited ? m_peakHi->value() : INFINITY;
    if (!(lo < hi)) {
        m_inputInfo->setText(tr("能区上限须大于下限。"));
        m_inputPlot->setCurves({}, label);
        refreshButtons();
        return;
    }
    auto wideSamples = SpectrumProcessing::selectSamples(wide, labels, lo, hi);
    const auto &axis = wideSamples.energyEv;
    auto &row = wideSamples.intensity;
    if (!limited && !axis.isEmpty()) {
        const QSignalBlocker low(m_peakLo), high(m_peakHi);
        const auto &edge = m_current.energyEv;
        const auto lower = separate && !edge.isEmpty() ? std::min(axis.first(), edge.first()) : axis.first();
        const auto upper = separate && !edge.isEmpty() ? std::max(axis.last(), edge.last()) : axis.last();
        m_peakLo->setValue(std::floor(lower * 1e6) / 1e6);
        m_peakHi->setValue(std::ceil(upper * 1e6) / 1e6);
    }
    const auto points = measuredRange(axis, row);
    QStringList info;
    if (!m_inputCrystal->currentData().toString().startsWith("crystal:"))
        info.append(tr("汇总 %1 个晶体").arg(labels.size()));
    info.append(points.count ? tr("宽谱 %1 至 %2 eV · %3 个点 · 最大间距 %4 eV · 最大采样点 %5 eV")
                                  .arg(formatNumber(points.lo), formatNumber(points.hi))
                                  .arg(points.count)
                                  .arg(formatNumber(points.maxSpacing), formatNumber(points.peakEnergy))
                             : tr("所选能区没有有效宽谱点。"));
    QVector<SpectrumCurve> curves;
    if (stitched && !m_current.energyEv.isEmpty()) {
        for (qsizetype i = 0; i < axis.size(); ++i)
            if (axis[i] >= m_current.energyEv.first() - 1e-7 && axis[i] <= m_current.energyEv.last() + 1e-7)
                row[i] = NAN;
    }
    if (std::any_of(row.cbegin(), row.cend(), [](double value) { return std::isfinite(value); }))
        curves.append({QStringLiteral("full_scan"), axis, row, QColor("#0072b2"), false, SpectrumCurve::Draw::Circles});
    auto coverage = info.join("    ·    ");
    if (separate) {
        info.clear();
        auto nearSamples = SpectrumProcessing::selectSamples(m_current, labels, lo, hi);
        const auto near = measuredRange(nearSamples.energyEv, nearSamples.intensity);
        if (!near.count) {
            info.append(tr("所选晶体和能区没有对应的近边谱。"));
        } else {
            info.append(tr("近边 %1 至 %2 eV · %3 个点")
                            .arg(formatNumber(near.lo), formatNumber(near.hi)).arg(near.count));
            if (nearScanCount > 1) {
                if (m_normalizeNear->isChecked()) {
                    for (auto &intensity : nearSamples.intensity)
                        intensity /= nearScanCount;
                    info.append(tr("近边谱已按 %1 次扫描取平均").arg(nearScanCount));
                } else {
                    info.append(tr("近边谱由 %1 次扫描叠加").arg(nearScanCount));
                }
            }
            if (points.count) {
                const double overlapLo = std::max(points.lo, near.lo), overlapHi = std::min(points.hi, near.hi);
                info.append(overlapHi >= overlapLo
                                ? tr("能区交集 %1 至 %2 eV").arg(formatNumber(overlapLo), formatNumber(overlapHi))
                                : tr("两段能区相离，间隔 %1 eV").arg(formatNumber(overlapLo - overlapHi)));
            }
            if (m_current.intensityUnit == wide.intensityUnit)
                curves.append({tr("近边谱"), nearSamples.energyEv, nearSamples.intensity, QColor("#d55e00"), false,
                               SpectrumCurve::Draw::Circles});
            else
                info.append(tr("强度单位不同，图中仅显示宽谱。"));
        }
        coverage += '\n' + info.join("    ·    ");
    }
    m_inputInfo->setText(coverage);
    m_inputPlot->setCurves(curves, label, "Energy transfer (eV)", wide.intensityUnit);
    refreshButtons();
}

void SpectrumProcessingPage::refreshTargets() {
    const auto &target = m_analysisInput;
    m_targets.clear();
    const int scope = m_scope->currentIndex();
    if (scope == PerCrystal) {
        m_targets = target.crystals;
        return;
    }

    QMap<QString, QStringList> groups;
    for (const auto &label : target.crystals.keys())
        groups[scope == TotalSpectrum ? QString(totalTarget) : label.section('-', 0, 0)].append(label);
    for (auto it = groups.cbegin(); it != groups.cend(); ++it) {
        const auto &labels = it.value();
        auto &result = m_targets[it.key()];
        result.fill(0., target.energyEv.size());
        for (const auto &label : labels) {
            const auto &row = target.crystals.constFind(label).value();
            for (qsizetype i = 0; i < row.size(); ++i)
                result[i] += row[i];
        }

        if (m_acquisition->currentIndex() != Separate || m_wide.crystalEnergyEv.isEmpty() ||
            m_current.energyEv.isEmpty() || labels.size() < 2)
            continue;
        // Wide scans have different native grids. Align the two wide flanks independently on
        // their shared measured coverage, without joining across the near edge or a scan gap.
        const double lo = std::nextafter(m_current.energyEv.first() - 1e-7, -INFINITY);
        const double hi = std::nextafter(m_current.energyEv.last() + 1e-7, INFINITY);
        for (const auto &[first, last] : {QPair<double, double>{-INFINITY, lo}, {hi, INFINITY}}) {
            const auto wide = SpectrumProcessing::selectSamples(m_wide, labels, first, last);
            for (qsizetype i = 0; i < wide.energyEv.size(); ++i) {
                const auto point = std::lower_bound(target.energyEv.cbegin(), target.energyEv.cend(), wide.energyEv[i]);
                if (point != target.energyEv.cend() && *point == wide.energyEv[i])
                    result[point - target.energyEv.cbegin()] = wide.intensity[i];
            }
        }
    }
}

void SpectrumProcessingPage::refreshAnalysisInput() {
    const auto label = m_fitCrystal->currentData().toString();
    refreshTargets();
    {
        const QSignalBlocker blockCombo(m_fitCrystal), blockBatch(m_batchTargets);
        m_fitCrystal->clear();
        m_batchTargets->clear();
        for (const auto &target : m_targets.keys()) {
            const auto name = targetName(target);
            m_fitCrystal->addItem(name, target);
            auto item = new QListWidgetItem(name, m_batchTargets);
            item->setData(Qt::UserRole, target);
            if (m_scope->currentIndex() == PerCrystal)
                item->setData(batchModuleRole, target.section('-', 0, 0));
            else if (m_scope->currentIndex() == ByModule)
                item->setData(batchModuleRole, target);
            item->setFlags(item->flags() | Qt::ItemIsUserCheckable);
            item->setCheckState(!m_configuredBatchTargets || m_configuredBatchTargets->contains(target)
                                    ? Qt::Checked : Qt::Unchecked);
        }
        if (const auto index = m_fitCrystal->findData(label); index >= 0)
            m_fitCrystal->setCurrentIndex(index);
    }
    resetAnalysis();
    const auto &target = m_analysisInput;
    const auto &e = target.energyEv;
    if (!e.isEmpty() && (m_resultLo->value() < e.first() || m_resultHi->value() > e.last() ||
                         m_resultHi->value() <= e.first() || m_resultLo->value() >= e.last())) {
        m_resultLo->setValue(e.first());
        m_resultHi->setValue(e.last());
    }
    if (target.crystals.isEmpty()) {
        m_analysisInputInfo->setText(tr("请在第 6 步点击「下一步」导入谱。"));
    } else {
        const auto source = m_acquisition->currentIndex() == Separate
            ? tr("第 6 步 Compton 拼接谱") : tr("第 6 步 Full Scan");
        m_analysisInputInfo->setText(tr("%1\n%2 个晶体 · %3 个能量点 · %4 至 %5 eV")
                                         .arg(source).arg(target.crystals.size()).arg(e.size())
                                         .arg(formatNumber(e.first()), formatNumber(e.last())));
    }
}

void SpectrumProcessingPage::refreshAnalysis() {
    const auto label = m_fitCrystal->currentData().toString();
    const auto input = m_targets.constFind(label);
    if (input == m_targets.cend()) {
        m_fitPlot->setCurves({}, {});
        m_fitInfo->clear();
        return;
    }
    const auto &energy = m_analysisInput.energyEv;
    const auto draw = m_acquisition->currentIndex() == Separate
        ? SpectrumCurve::Draw::Circles : SpectrumCurve::Draw::Samples;
    const auto backgroundDraw = draw == SpectrumCurve::Draw::Circles ? draw : SpectrumCurve::Draw::Line;
    QVector<SpectrumCurve> curves{{tr("输入谱"), energy, *input, QColor("#0072b2"), false, draw}};
    QStringList info;
    const double q = m_analysisInput.geometry.value(label).qMean;
    if (std::isfinite(q))
        info.append(tr("q %1 Å⁻¹").arg(formatNumber(q)));
    const auto preview = m_previews.constFind(label);
    const auto finished = m_analysis.spectra.constFind(label);
    const auto fit = finished != m_analysis.spectra.cend() ? &*finished
        : preview != m_previews.cend() ? &*preview : nullptr;
    if (fit && !fit->terms.isEmpty()) {
        for (auto it = fit->terms.cbegin(); it != fit->terms.cend(); ++it) {
            const auto name = componentName(it.key());
            curves.append({name, energy, it->background, backgroundColors[int(it.key())], false, backgroundDraw});
            QString details = tr("%1 · R² %2").arg(name, formatNumber(it->rSquared));
            if (std::isfinite(it->centerEv))
                details += tr(" · 峰心 %1 eV · FWHM %2 eV")
                               .arg(formatNumber(it->centerEv), formatNumber(it->fwhmEv));
            if (it->parameters.contains("scale"))
                details += tr(" · 比例 %1").arg(formatNumber(it->parameters.value("scale").toDouble()));
            info.append(details);
        }
        if (fit->terms.size() > 1)
            curves.append({tr("总背景"), energy, fit->total.background, QColor("#555555"), true, backgroundDraw});
        curves.append({tr("扣除后"), energy, fit->total.corrected, QColor("#009e73"), false, draw});
    }
    if (m_analysis.failures.contains(label))
        info.append(tr("批量扣除未完成：%1").arg(m_analysis.failures.value(label)));
    m_fitPlot->setCurves(curves, targetName(label), "Energy transfer (eV)", m_analysisInput.intensityUnit);
    m_fitInfo->setText(info.join('\n'));
}

QStringList SpectrumProcessingPage::resultLabels() const {
    return m_resultChoice->currentData().toStringList();
}

void SpectrumProcessingPage::refreshResult() {
    const auto labels = resultLabels();
    const double lo = m_resultLo->value(), hi = m_resultHi->value();
    const auto &target = m_analysisInput;
    const auto &axis = target.energyEv;
    const bool stitched = m_acquisition->currentIndex() == Separate;
    QVector<SpectrumCurve> curves;
    QVector<double> x, corrected;
    if (hi > lo && !labels.isEmpty()) {
        const auto first = std::lower_bound(axis.cbegin(), axis.cend(), lo) - axis.cbegin();
        const auto last = std::upper_bound(axis.cbegin(), axis.cend(), hi) - axis.cbegin();
        x = axis.mid(first, last - first);
        auto average = [&](const auto &values) {
            QVector<double> y(x.size(), 0.);
            for (const auto &label : labels) {
                const auto &row = values(label);
                for (qsizetype i = 0; i < y.size(); ++i)
                    y[i] += row[first + i];
            }
            for (auto &value : y)
                value /= labels.size();
            return y;
        };
        auto spectrum = [&](const QString &label) -> const StagedSpectrumFit & {
            return m_analysis.spectra.constFind(label).value();
        };
        corrected = average([&](const QString &label) -> const QVector<double> & {
            return spectrum(label).total.corrected;
        });
        const auto draw = stitched ? SpectrumCurve::Draw::Circles : SpectrumCurve::Draw::Samples;
        const auto backgroundDraw = stitched ? draw : SpectrumCurve::Draw::Line;
        if (m_resultView->currentIndex() != 0) {
            curves.append({tr("输入谱"), x, average([&](const QString &label) -> const QVector<double> & {
                               return m_targets.constFind(label).value();
                           }), QColor("#0072b2"), false, draw});
            const auto &terms = spectrum(labels.first()).terms;
            for (auto it = terms.cbegin(); it != terms.cend(); ++it) {
                const auto term = it.key();
                curves.append({componentName(term), x, average([&](const QString &label) -> const QVector<double> & {
                                   return spectrum(label).terms.constFind(term)->background;
                               }), backgroundColors[int(term)], false, backgroundDraw});
            }
            if (terms.size() > 1)
                curves.append({tr("总背景"), x, average([&](const QString &label) -> const QVector<double> & {
                                   return spectrum(label).total.background;
                               }), QColor("#555555"), true, backgroundDraw});
        }
        curves.append({tr("扣背底结果"), x, corrected, QColor("#009e73"), false, draw});
    }
    m_resultPlot->setCurves(curves, m_resultChoice->currentText(), "Energy transfer (eV)", target.intensityUnit);
    qsizetype peak = -1, finite = 0;
    for (qsizetype i = 0; i < corrected.size(); ++i)
        if (std::isfinite(corrected[i])) {
            ++finite;
            if (peak < 0 || corrected[i] > corrected[peak])
                peak = i;
        }
    m_resultInfo->setText(peak < 0 ? tr("当前能区没有有效结果，请运行扣背底或调整能量范围。")
                                   : tr("%1 条谱 · %2 个有效能量点 · 最大采样点 %3 eV，强度 %4")
                                         .arg(labels.size()).arg(finite)
                                         .arg(formatNumber(x[peak]), formatNumber(corrected[peak])));
    for (auto button : {m_saveResult, m_saveText})
        button->setEnabled(!m_busy && finite >= 2);
}

void SpectrumProcessingPage::saveResult(ResultFormat format) {
    const auto labels = resultLabels();
    if (m_busy || labels.isEmpty() || m_resultLo->value() >= m_resultHi->value())
        return;
    const auto &target = m_analysisInput;
    const auto &axis = target.energyEv;
    const auto first = std::lower_bound(axis.cbegin(), axis.cend(), m_resultLo->value()) - axis.cbegin();
    const auto last = std::upper_bound(axis.cbegin(), axis.cend(), m_resultHi->value()) - axis.cbegin();
    if (last - first < 2)
        return;
    const bool plainText = format == ResultFormat::Text;
    const auto exportName = (m_exportStem.isEmpty() ? QStringLiteral("background_corrected") : m_exportStem + "_near_edge") +
                            (plainText ? ".txt" : ".h5");
    const auto path = WorkspacePaths::saveFile(
        m_resultPage, Workflow::Spectroscopy, Kind::Result,
        m_exportStem.isEmpty() ? target.sources.value(0) : QString{},
        tr("保存近边结构结果"), exportName,
        plainText ? tr("制表符分隔文本 (*.txt *.tsv)") : tr("HDF5 (*.h5 *.hdf5)"));
    if (path.isEmpty())
        return;

    ProcessedSpectrum spectrum;
    spectrum.energyEv = axis.mid(first, last - first);
    spectrum.sources = target.sources;
    spectrum.processing = target.processing;
    for (auto it = m_exportIdentity.constBegin(); it != m_exportIdentity.constEnd(); ++it)
        spectrum.processing.insert(it.key(), it.value());
    for (auto it = m_analysis.parameters.constBegin(); it != m_analysis.parameters.constEnd(); ++it)
        spectrum.processing.insert(it.key(), it.value());
    spectrum.intensityUnit = target.intensityUnit;
    QMap<QString, CrystalBackground> fits;
    QMap<QString, QMap<QString, QVector<double>>> components;
    QJsonObject fitDetails;
    for (const auto &label : labels) {
        const auto &item = m_analysis.spectra.constFind(label).value();
        spectrum.crystals[label] = m_targets.constFind(label)->mid(first, last - first);
        if (target.geometry.contains(label))
            spectrum.geometry[label] = target.geometry.value(label);
        if (target.pointSources.contains(label))
            spectrum.pointSources[label] = target.pointSources.value(label).mid(first, last - first);
        auto fit = item.total;
        fit.background = fit.background.mid(first, last - first);
        fit.corrected = fit.corrected.mid(first, last - first);
        fits.insert(label, std::move(fit));
        QJsonObject details;
        for (auto it = item.terms.cbegin(); it != item.terms.cend(); ++it) {
            const auto name = backgroundComponentKey(it.key());
            components[name].insert(label, it->background.mid(first, last - first));
            auto parameters = it->parameters;
            if (std::isfinite(it->rSquared))
                parameters["r_squared"] = it->rSquared;
            if (std::isfinite(it->centerEv)) {
                parameters["center_ev"] = it->centerEv;
                parameters["fwhm_ev"] = it->fwhmEv;
                parameters["exponent"] = it->exponent;
                parameters["height"] = it->peakHeight;
            }
            details[name] = parameters;
        }
        fitDetails[label] = details;
    }
    spectrum.processing["background_fits"] = fitDetails;
    start(tr("正在保存扣除结果…"),
          [path, spectrum, fits, components, plainText] {
              if (plainText)
                  SpectrumProcessing::saveText(path, spectrum, fits, components);
              else
                  SpectrumProcessing::save(path, spectrum, fits, components);
              Work result;
              result.exportPath = path;
              result.exportFormat = plainText ? "near_edge_text" : "near_edge_hdf5";
              result.message = tr("已保存扣除结果 %1。").arg(path);
              return result;
          }, false);
}

void SpectrumProcessingPage::refreshButtons() {
    m_splice->setEnabled(!m_busy && m_stitched.crystals.isEmpty() &&
                        !m_wide.crystals.isEmpty() && !m_current.crystals.isEmpty());
    const auto label = m_fitCrystal->currentData().toString();
    const int previewIndex = m_fitCrystal->currentIndex();
    m_previousPreview->setEnabled(!m_busy && previewIndex > 0);
    m_nextPreview->setEnabled(!m_busy && previewIndex >= 0 && previewIndex + 1 < m_fitCrystal->count());
    const bool hasTargets = m_batchTargets->count() > 0;
    m_batchModule->setEnabled(!m_busy && hasTargets && m_scope->currentIndex() != TotalSpectrum);
    m_batchAll->setEnabled(!m_busy && hasTargets);
    m_batchInvert->setEnabled(!m_busy && hasTargets);
    const auto preview = m_previews.constFind(label);
    for (int i = 0; i < backgroundComponentCount; ++i) {
        const auto component = BackgroundComponent(i);
        const bool available = componentAvailable(component);
        const bool enabled = m_enabledTerms[i]->isChecked();
        const bool fitted = preview != m_previews.cend() && preview->terms.contains(component);
        m_enabledTerms[i]->setEnabled(available);
        m_termControls[i]->setVisible(available && enabled);
        m_fitTerms[i]->setEnabled(canFit(component));
        m_termStatus[i]->setText(!available ? tr("仅分晶体可用") : !enabled ? QString()
                                     : fitted ? tr("已拟合") : tr("待拟合"));
    }
    m_run->setEnabled(!m_busy && previewComplete() && !batchTargets().isEmpty());
}
