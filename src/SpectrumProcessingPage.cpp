#include "SpectrumProcessingPage.h"
#include "SpectrumPlot.h"
#include "WorkflowUi.h"
#include <QtConcurrent/QtConcurrentRun>
#include <QtWidgets>

namespace {
const QString inputFilter = QStringLiteral("Spectra (*.h5 *.hdf5 *.nxs *.txt *.tsv);;All files (*)");
QDoubleSpinBox *numberBox(double min, double max, double value, int precision = 4) {
    auto box = new QDoubleSpinBox;
    box->setDecimals(precision);
    box->setRange(min, max);
    box->setValue(value);
    box->setKeyboardTracking(false);
    return box;
}
SpectrumPlot *finePlot(const QString &name) {
    auto plot = new SpectrumPlot;
    plot->setFineGrid(true);
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
QString savePath(QWidget *parent, const QString &name) {
    auto path = QFileDialog::getSaveFileName(parent, QObject::tr("保存谱学结果"), name, "HDF5 (*.h5 *.hdf5)");
    if (!path.isEmpty() && QFileInfo(path).suffix().isEmpty())
        path += ".h5";
    return path;
}
} // namespace

SpectrumProcessingPage::SpectrumProcessingPage(QWidget *parent) : QWidget(parent) {
    auto input = new QVBoxLayout(this);
    m_inputControls = new QWidget;
    auto controls = new QVBoxLayout(m_inputControls);
    controls->setContentsMargins(0, 0, 0, 0);
    auto commands = new QHBoxLayout;
    m_useCurrent = new QPushButton(tr("用当前原始谱作宽谱"));
    m_useCurrent->setObjectName("LoadCurrentSpectrumForProcessing");
    auto wide = new QPushButton(tr("载入康普顿宽谱…"));
    m_edges = new QPushButton(tr("选择近边谱并替换…"));
    m_geometry = new QPushButton(tr("补充晶体几何…"));
    m_saveStitch = new QPushButton(tr("保存拼接 HDF5…"));
    for (auto button : {m_useCurrent, wide, m_edges, m_geometry})
        commands->addWidget(button);
    commands->addStretch();
    commands->addWidget(m_saveStitch);
    controls->addLayout(commands);
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
    controls->addWidget(WorkflowUi::disclosure(tr("外部 HDF5 数据集映射"), mapping));
    auto selection = new QHBoxLayout;
    m_inputCrystal = new QComboBox;
    m_inputInfo = textLabel(tr("宽谱与近边分开采集时可在此拼接。直接做线性扣除可进入第 7 步。"));
    selection->addWidget(new QLabel(tr("预览晶体")));
    selection->addWidget(m_inputCrystal);
    selection->addWidget(m_inputInfo, 1);
    controls->addLayout(selection);
    input->addWidget(m_inputControls);
    m_inputPlot = finePlot("StitchedSpectrum");
    input->addWidget(m_inputPlot, 1);
    input->addWidget(m_status[0] = textLabel());

    m_analysisPage = new QWidget;
    auto analysis = new QVBoxLayout(m_analysisPage);
    auto columns = new QHBoxLayout;
    m_fitControls = new QWidget;
    auto fitLayout = new QVBoxLayout(m_fitControls);
    fitLayout->setContentsMargins(0, 0, 0, 0);
    auto sourceBox = new QGroupBox(tr("扣背底输入"));
    auto sourceFields = new QFormLayout(sourceBox);
    m_analysisSource = new QComboBox;
    m_analysisSource->setObjectName("BackgroundSpectrumSource");
    m_analysisSource->addItems({tr("第 5 步原始谱"), tr("外部逐晶体谱"), tr("第 6 步宽谱或拼接谱")});
    sourceFields->addRow(tr("输入谱"), m_analysisSource);
    m_loadDirect = new QPushButton(tr("载入外部逐晶体谱…"));
    m_loadDirect->setToolTip(tr("载入已处理的 HDF5 或逐晶体文本谱。自定义 HDF5 路径可在第 6 步设置。"));
    sourceFields->addRow(m_loadDirect);
    m_analysisInputInfo = textLabel();
    sourceFields->addRow(m_analysisInputInfo);
    m_backgroundMode = new QComboBox;
    m_backgroundMode->setObjectName("BackgroundMode");
    m_backgroundMode->addItems({tr("全部晶体用线性背底"), tr("高 q 用康普顿法，低 q 用线性背底")});
    sourceFields->addRow(tr("处理方式"), m_backgroundMode);
    fitLayout->addWidget(sourceBox);

    auto comptonBox = new QGroupBox(tr("康普顿背底参数"));
    auto fields = new QFormLayout(comptonBox);
    m_qMin = numberBox(0, 1000, 9.5);
    m_qMin->setObjectName("HighMomentumMinimum");
    fields->addRow(tr("高 q 最低值 / Å⁻¹"), m_qMin);
    m_qMin->setToolTip(tr("按每个晶体的平均动量转移 q 分流。达到阈值的晶体使用康普顿法。"));
    m_hfTable = new QLineEdit;
    m_hfTable->setPlaceholderText("ComptonProfiles.dat");
    QSettings saved;
    m_hfTable->setText(saved.value("spectroscopy/hfTable").toString());
    auto browse = new QPushButton(tr("选择…"));
    auto pathRow = new QHBoxLayout;
    pathRow->addWidget(m_hfTable, 1);
    pathRow->addWidget(browse);
    fields->addRow(tr("HF 壳层表"), pathRow);
    m_composition = new QLineEdit;
    m_composition->setText(saved.value("spectroscopy/hfComposition").toString());
    m_composition->setPlaceholderText("Li:1, Ni:0.9, Co:0.05, Mn:0.05, O:2");
    fields->addRow(tr("元素及原子计量"), m_composition);
    m_coreScale = numberBox(0, 100, 1);
    m_cutoff = numberBox(.01, 10000, 20);
    fields->addRow(tr("HF 芯谱比例"), m_coreScale);
    fields->addRow(tr("价电子结合能上限 / eV"), m_cutoff);
    const QStringList names{tr("HF 归一化"), tr("芯谱尾部直线"), tr("价电子轮廓拟合"), tr("高 q 残余直线")};
    const QStringList ranges{"80:970", "700:740, 930:970", "120:510, 580:620, 680:740, 920:970", "460:510"};
    for (int i = 0; i < names.size(); ++i) {
        m_windows[i] = new QLineEdit(ranges[i]);
        fields->addRow(names[i] + tr(" / eV"), m_windows[i]);
    }
    auto linearBox = new QGroupBox(tr("线性背底参数"));
    auto linearFields = new QFormLayout(linearBox);
    m_windows[4] = new QLineEdit("460:510");
    m_windows[4]->setObjectName("LowMomentumFitWindows");
    linearFields->addRow(tr("拟合区间 / eV"), m_windows[4]);
    linearFields->addRow(textLabel(
        tr("区间格式为 起点:终点，多段用逗号分隔。初始区间取自参考程序的 NCM 示例，请按实际扫描修改。")));
    m_run = new QPushButton(tr("逐晶体扣背底"));
    m_run->setObjectName("RunMomentumBackground");
    fitLayout->addWidget(linearBox);
    fitLayout->addWidget(comptonBox);
    fitLayout->addWidget(m_run);
    fitLayout->addStretch();
    auto scroll = new QScrollArea;
    scroll->setWidget(m_fitControls);
    scroll->setWidgetResizable(true);
    scroll->setMinimumWidth(380);
    scroll->setMaximumWidth(505);
    columns->addWidget(scroll);
    auto review = new QVBoxLayout;
    auto reviewRow = new QHBoxLayout;
    m_fitCrystal = new QComboBox;
    reviewRow->addWidget(new QLabel(tr("诊断晶体")));
    reviewRow->addWidget(m_fitCrystal);
    reviewRow->addStretch();
    review->addLayout(reviewRow);
    auto diagnostic = new QTabWidget;
    m_fitPlot = finePlot("MomentumBackgroundDiagnostics");
    m_referencePlot = finePlot("ValenceMomentumReference");
    diagnostic->addTab(m_fitPlot, tr("晶体扣除诊断"));
    diagnostic->addTab(m_referencePlot, tr("高 q 价电子轮廓"));
    diagnostic->setTabVisible(1, false);
    comptonBox->hide();
    review->addWidget(diagnostic, 1);
    m_fitInfo = textLabel();
    review->addWidget(m_fitInfo);
    columns->addLayout(review, 1);
    analysis->addLayout(columns, 1);
    auto runStatus = new QHBoxLayout;
    runStatus->addWidget(m_status[1] = textLabel(), 1);
    m_cancelButton = new QPushButton(tr("取消分析"));
    m_cancelButton->hide();
    runStatus->addWidget(m_cancelButton);
    analysis->addLayout(runStatus);

    m_resultPage = new QWidget;
    auto results = new QVBoxLayout(m_resultPage);
    m_resultControls = new QWidget;
    auto resultRow = new QHBoxLayout(m_resultControls);
    resultRow->setContentsMargins(0, 0, 0, 0);
    m_resultChoice = new QComboBox;
    m_resultLo = numberBox(-1e8, 1e8, 520, 6);
    m_resultHi = numberBox(-1e8, 1e8, 600, 6);
    m_resultLo->setObjectName("CorrectedEnergyMinimum");
    m_resultHi->setObjectName("CorrectedEnergyMaximum");
    m_saveResult = new QPushButton(tr("保存当前能区 HDF5…"));
    resultRow->addWidget(new QLabel(tr("结果谱")));
    resultRow->addWidget(m_resultChoice, 1);
    resultRow->addWidget(new QLabel(tr("能量范围 / eV")));
    resultRow->addWidget(m_resultLo);
    resultRow->addWidget(new QLabel(tr("至")));
    resultRow->addWidget(m_resultHi);
    resultRow->addWidget(m_saveResult);
    results->addWidget(m_resultControls);
    m_resultPlot = finePlot("CorrectedCrystalSpectrum");
    results->addWidget(m_resultPlot, 1);
    results->addWidget(m_resultInfo = textLabel());
    results->addWidget(m_status[2] = textLabel());

    connect(m_useCurrent, &QPushButton::clicked, this, [this] {
        m_wide = m_current;
        m_stitched = m_wide;
        refreshInput();
        if (m_analysisSource->currentIndex() == StitchingSpectrum)
            refreshAnalysisInput();
        message(tr("已载入当前原始谱快照。"));
    });
    connect(wide, &QPushButton::clicked, this, [this] {
        auto path = QFileDialog::getOpenFileName(this, tr("选择逐晶体康普顿宽谱"), {}, inputFilter);
        if (path.isEmpty())
            return;
        auto mapping = hdfLayout();
        start(
            tr("正在载入宽谱…"),
            [path, mapping] {
                Work result;
                result.wide = SpectrumProcessing::load(path, mapping);
                result.message = tr("已载入宽谱 %1。").arg(QFileInfo(path).fileName());
                return result;
            },
            false);
    });
    connect(m_edges, &QPushButton::clicked, this, [this] {
        auto files =
            QFileDialog::getOpenFileNames(this, tr("选择要替换宽谱对应能区的近边文件"), {}, inputFilter);
        if (files.isEmpty())
            return;
        auto wideData = m_wide;
        auto mapping = hdfLayout();
        start(tr("正在对齐同名晶体和能量点…"), [this, files, wideData, mapping] {
            QVector<ProcessedSpectrum> edges;
            QSet<QString> seen;
            for (const auto &file : files) {
                if (m_cancel.load())
                    return Work{};
                auto path = QFileInfo(file).canonicalFilePath();
                if (!seen.contains(path))
                    edges.append(SpectrumProcessing::load(file, mapping));
                seen.insert(path);
            }
            Work result;
            result.stitched = SpectrumProcessing::replaceEdges(wideData, edges, m_cancel);
            result.message = tr("已替换 %1 个近边扫描的同名晶体和对应能量点。").arg(edges.size());
            return result;
        });
    });
    connect(m_saveStitch, &QPushButton::clicked, this, [this] {
        auto path = savePath(this, "stitched_spectra.h5");
        if (path.isEmpty())
            return;
        auto spectrum = m_stitched;
        start(
            tr("正在保存拼接谱…"),
            [path, spectrum] {
                SpectrumProcessing::save(path, spectrum);
                Work result;
                result.message = tr("已保存拼接谱 %1。").arg(path);
                return result;
            },
            false);
    });
    connect(m_geometry, &QPushButton::clicked, this, [this] {
        auto path = QFileDialog::getOpenFileName(this, tr("选择宽谱对应的逐晶体 ROI 几何表"), {},
                                                 "ROI geometry (*.txt *.tsv);;All files (*)");
        if (path.isEmpty())
            return;
        auto wideData = m_wide, target = m_stitched;
        start(
            tr("正在读取逐晶体几何…"),
            [path, wideData, target]() mutable {
                wideData.geometry = SpectrumProcessing::loadGeometry(path, wideData);
                if (wideData.geometry.isEmpty())
                    throw std::runtime_error(tr("没有匹配到可用的晶体几何记录。").toUtf8().constData());
                target.geometry = wideData.geometry;
                Work result;
                result.wide = std::move(wideData);
                result.stitched = std::move(target);
                result.message = tr("已更新宽谱和拼接谱的晶体几何。");
                return result;
            },
            false);
    });
    connect(browse, &QPushButton::clicked, this, [this] {
        auto path = QFileDialog::getOpenFileName(m_analysisPage, tr("选择 DABAX HF 壳层数据"),
                                                 m_hfTable->text(), "DABAX (*.dat);;All files (*)");
        if (!path.isEmpty())
            m_hfTable->setText(path);
    });
    connect(m_loadDirect, &QPushButton::clicked, this, [this] {
        const auto path =
            QFileDialog::getOpenFileName(m_analysisPage, tr("选择外部逐晶体谱"), {}, inputFilter);
        if (path.isEmpty())
            return;
        auto mapping = hdfLayout();
        start(
            tr("正在载入逐晶体谱…"),
            [path, mapping] {
                Work result;
                result.external = SpectrumProcessing::load(path, mapping);
                result.message = tr("已载入逐晶体谱 %1。").arg(QFileInfo(path).fileName());
                return result;
            },
            false);
    });
    connect(m_analysisSource, &QComboBox::currentIndexChanged, this, [this] { refreshAnalysisInput(); });
    connect(m_backgroundMode, &QComboBox::currentIndexChanged, this,
            [this, comptonBox, diagnostic](int mode) {
                comptonBox->setVisible(mode != 0);
                diagnostic->setTabVisible(1, mode != 0);
                refreshAnalysisInput();
            });
    connect(m_run, &QPushButton::clicked, this, [this] {
        try {
            const bool hybrid = m_backgroundMode->currentIndex() != 0;
            ComptonSettings parameters;
            if (hybrid) {
                parameters = settings();
                QSettings saved;
                saved.setValue("spectroscopy/hfTable", parameters.hfTable);
                saved.setValue("spectroscopy/hfComposition", parameters.composition);
            } else {
                parameters.lowLinear = parseWindows(m_windows[4]->text(), tr("线性拟合区间"));
            }
            const auto target = analysisInput();
            const auto wideData = hybrid ? m_wide : ProcessedSpectrum{};
            const auto revision = m_analysisRevision;
            start(hybrid ? tr("正在逐晶体计算 HF 归一化、价电子轮廓和背底…") : tr("正在逐晶体拟合线性背底…"),
                  [this, wideData, target, parameters, hybrid, revision] {
                      Work result;
                      result.analysis = hybrid ? analyzeCompton(wideData, target, parameters, m_cancel)
                                               : analyzeLinear(target, parameters.lowLinear, m_cancel);
                      result.analysisRevision = revision;
                      const auto &a = *result.analysis;
                      result.error = !a.failures.isEmpty();
                      result.message = tr("已处理 %1 个晶体，未生成结果 %2 个。切换晶体可查看原因。")
                                           .arg(a.crystals.size())
                                           .arg(a.failures.size());
                      return result;
                  });
        } catch (const std::exception &ex) {
            message(QString::fromUtf8(ex.what()), true);
        }
    });
    connect(m_cancelButton, &QPushButton::clicked, this, [this] {
        m_cancel.store(true);
        message(tr("正在取消…"));
    });
    connect(m_inputCrystal, &QComboBox::currentIndexChanged, this, [this] { refreshInput(); });
    connect(m_fitCrystal, &QComboBox::currentIndexChanged, this, [this] { refreshAnalysis(); });
    connect(m_qMin, &QDoubleSpinBox::valueChanged, this, [this] {
        resetAnalysis();
        refreshInput();
    });
    for (auto edit :
         {m_hfTable, m_composition, m_windows[0], m_windows[1], m_windows[2], m_windows[3], m_windows[4]})
        connect(edit, &QLineEdit::textChanged, this, [this] { resetAnalysis(); });
    for (auto spin : {m_coreScale, m_cutoff})
        connect(spin, &QDoubleSpinBox::valueChanged, this, [this] { resetAnalysis(); });
    connect(m_resultChoice, &QComboBox::currentIndexChanged, this, [this] { refreshResult(); });
    for (auto spin : {m_resultLo, m_resultHi})
        connect(spin, &QDoubleSpinBox::valueChanged, this, [this] { refreshResult(); });
    connect(m_saveResult, &QPushButton::clicked, this, [this] {
        const auto labels = resultLabels();
        if (labels.isEmpty())
            return;
        const auto &target = analysisInput();
        const auto &axis = target.energyEv;
        const auto first = std::lower_bound(axis.cbegin(), axis.cend(), m_resultLo->value()) - axis.cbegin();
        const auto last = std::upper_bound(axis.cbegin(), axis.cend(), m_resultHi->value()) - axis.cbegin();
        if (last - first < 2)
            return;
        ProcessedSpectrum spectrum;
        spectrum.energyEv = axis.mid(first, last - first);
        spectrum.sources = target.sources;
        spectrum.processing = m_analysis.parameters;
        QJsonObject normalizationScales;
        spectrum.intensityUnit =
            m_analysis.crystals.value(labels.first()).high ? "HF-normalized S (1/eV)" : target.intensityUnit;
        QMap<QString, CrystalBackground> fits;
        for (const auto &label : labels) {
            const auto &item = m_analysis.crystals[label];
            spectrum.crystals[label] = item.input.mid(first, last - first);
            if (target.geometry.contains(label))
                spectrum.geometry[label] = target.geometry.value(label);
            normalizationScales[label] = item.scale;
            if (target.pointSources.contains(label))
                spectrum.pointSources[label] = target.pointSources.value(label).mid(first, last - first);
            auto fit = item.fit;
            fit.background = fit.background.mid(first, last - first);
            fit.corrected = fit.corrected.mid(first, last - first);
            fits.insert(label, std::move(fit));
        }
        spectrum.processing["normalization_scales"] = normalizationScales;
        auto path = savePath(m_resultPage, "background_corrected.h5");
        if (path.isEmpty())
            return;
        start(
            tr("正在保存扣除结果…"),
            [path, spectrum, fits] {
                SpectrumProcessing::save(path, spectrum, fits);
                Work result;
                result.message = tr("已保存扣除结果 %1。").arg(path);
                return result;
            },
            false);
    });
    connect(&m_watcher, &QFutureWatcher<Work>::finished, this, [this] {
        auto result = m_watcher.result();
        m_busy = false;
        const bool updatedStitching = result.wide || result.stitched;
        if (result.wide) {
            m_wide = std::move(*result.wide);
            m_stitched = m_wide;
            if (!m_wide.wideCrystals.isEmpty()) {
                m_wide.crystals = std::move(m_wide.wideCrystals);
                m_wide.pointSources.clear();
                m_wide.sources = {m_stitched.sources.value(0)};
            }
        }
        if (result.stitched)
            m_stitched = std::move(*result.stitched);
        if (updatedStitching) {
            refreshInput();
            if (m_analysisSource->currentIndex() == StitchingSpectrum)
                refreshAnalysisInput();
        }
        if (result.external) {
            m_external = std::move(*result.external);
            if (m_analysisSource->currentIndex() == ExternalSpectrum)
                refreshAnalysisInput();
            else
                m_analysisSource->setCurrentIndex(ExternalSpectrum);
        }
        if (result.analysis && result.analysisRevision != m_analysisRevision) {
            result.analysis.reset();
            result.message = tr("输入谱已更新，请重新扣背底。");
            result.error = false;
        }
        if (result.analysis) {
            m_analysis = std::move(*result.analysis);
            refreshAnalysis();
            QSignalBlocker block(m_resultChoice);
            m_resultChoice->clear();
            QStringList high, low;
            QMap<QString, QStringList> groups;
            const auto &geometry = analysisInput().geometry;
            for (auto it = m_analysis.crystals.cbegin(); it != m_analysis.crystals.cend(); ++it) {
                (it->high ? high : low).append(it.key());
                double q = geometry.value(it.key()).qMean;
                if (!std::isfinite(q) || q <= 0)
                    continue;
                double lo = .5 * std::floor(q * 2);
                auto group = tr("%1 · %2 · q %3–%4 Å⁻¹")
                                 .arg(it->high ? tr("康普顿") : tr("线性"), it.key().section('-', 0, 0))
                                 .arg(lo, 0, 'f', 1)
                                 .arg(lo + .5, 0, 'f', 1);
                groups[group].append(it.key());
            }
            if (!high.isEmpty())
                m_resultChoice->addItem(tr("高 q 晶体等权平均"), high);
            if (!low.isEmpty())
                m_resultChoice->addItem(m_backgroundMode->currentIndex() == 0 ? tr("全部晶体等权平均")
                                                                              : tr("低 q 晶体等权平均"),
                                        low);
            for (auto group = groups.cbegin(); group != groups.cend(); ++group)
                m_resultChoice->addItem(group.key(), group.value());
            for (const auto &label : m_analysis.crystals.keys())
                m_resultChoice->addItem(label, QStringList{label});
            refreshResult();
        }
        message(result.message.isEmpty() ? tr("已取消。") : result.message, result.error);
        m_inputControls->setEnabled(true);
        m_fitControls->setEnabled(true);
        m_resultControls->setEnabled(true);
        m_cancelButton->hide();
        refreshButtons();
    });
    refreshInput();
    refreshAnalysisInput();
}

SpectrumProcessingPage::~SpectrumProcessingPage() {
    m_cancel.store(true);
    m_watcher.waitForFinished();
}

void SpectrumProcessingPage::setCurrentResult(ProcessedSpectrum result) {
    m_current = std::move(result);
    if (m_analysisSource->currentIndex() == CurrentSpectrum)
        refreshAnalysisInput();
    else
        refreshButtons();
}

const ProcessedSpectrum &SpectrumProcessingPage::analysisInput() const {
    switch (m_analysisSource->currentIndex()) {
    case ExternalSpectrum:
        return m_external;
    case StitchingSpectrum:
        return m_stitched;
    default:
        return m_current;
    }
}

SpectrumHdfLayout SpectrumProcessingPage::hdfLayout() const {
    return {m_energyPath->text().trimmed(),
            m_labelsPath->text().trimmed(),
            m_intensityPath->text().trimmed(),
            m_unit->text().trimmed(),
            m_kev->isChecked(),
            m_energyFirst->isChecked()};
}

ComptonSettings SpectrumProcessingPage::settings() const {
    ComptonSettings out;
    out.highQMin = m_qMin->value();
    out.coreScale = m_coreScale->value();
    out.valenceCutoffEv = m_cutoff->value();
    out.hfTable = m_hfTable->text().trimmed();
    out.composition = m_composition->text().trimmed();
    SpectrumWindows *windows[]{&out.normalization, &out.coreTail, &out.cleanReference, &out.highResidual,
                               &out.lowLinear};
    for (int i = 0; i < 5; ++i)
        *windows[i] = parseWindows(m_windows[i]->text(), tr("第 %1 项能区").arg(i + 1));
    return out;
}

void SpectrumProcessingPage::message(const QString &text, bool error) {
    for (auto label : m_status) {
        label->setText(text);
        label->setStyleSheet(error ? "color:#a02020" : "");
    }
}

void SpectrumProcessingPage::start(const QString &text, std::function<Work()> job, bool cancellable) {
    m_busy = true;
    for (auto controls : {m_inputControls, m_fitControls, m_resultControls})
        controls->setEnabled(false);
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

void SpectrumProcessingPage::resetAnalysis() {
    ++m_analysisRevision;
    m_analysis = {};
    m_resultChoice->clear();
    refreshAnalysis();
    refreshResult();
    refreshButtons();
}

void SpectrumProcessingPage::refreshInput() {
    if (m_stitched.crystals.isEmpty())
        return;
    auto label = m_inputCrystal->currentText();
    const auto names = m_stitched.crystals.keys();
    QStringList shown;
    for (int i = 0; i < m_inputCrystal->count(); ++i)
        shown.append(m_inputCrystal->itemText(i));
    if (shown != names) {
        QSignalBlocker inputBlock(m_inputCrystal);
        m_inputCrystal->clear();
        m_inputCrystal->addItems(names);
        label = m_inputCrystal->currentText();
    }
    const auto &e = m_stitched.energyEv;
    int high = 0, low = 0;
    for (const auto &name : names) {
        double q = m_wide.geometry.value(name).qMean;
        if (std::isfinite(q) && q > 0)
            (q >= m_qMin->value() ? high : low)++;
    }
    qsizetype count = 0;
    auto origins = m_stitched.pointSources.constFind(label);
    if (origins != m_stitched.pointSources.cend())
        count = std::count_if(origins->cbegin(), origins->cend(), [](int n) { return n != 0; });
    m_inputInfo->setText(tr("%1 个晶体 · %2 个能量点 · 高 q %3 · 低 q %4 · 缺少 q %5 · 当前晶体替换 %6 点")
                             .arg(names.size())
                             .arg(e.size())
                             .arg(high)
                             .arg(low)
                             .arg(names.size() - high - low)
                             .arg(count));
    QVector<SpectrumCurve> curves{{tr("康普顿宽谱"), e, m_wide.crystals.value(label), QColor("#787878")}};
    if (!m_stitched.pointSources.isEmpty())
        curves.append({tr("近边替换后的谱"), e, m_stitched.crystals.value(label), QColor("#0072b2")});
    m_inputPlot->setCurves(curves, label + tr(" · 宽谱与拼接预览"), "Energy transfer (eV)",
                           m_stitched.intensityUnit);
    refreshButtons();
}

void SpectrumProcessingPage::refreshAnalysisInput() {
    const auto &target = analysisInput();
    const auto &e = target.energyEv;
    const auto label = m_fitCrystal->currentText();
    {
        QSignalBlocker block(m_fitCrystal);
        m_fitCrystal->clear();
        m_fitCrystal->addItems(target.crystals.keys());
        if (m_fitCrystal->findText(label) >= 0)
            m_fitCrystal->setCurrentText(label);
    }
    resetAnalysis();
    if (!e.isEmpty() && (m_resultLo->value() < e.first() || m_resultHi->value() > e.last() ||
                         m_resultHi->value() <= e.first() || m_resultLo->value() >= e.last())) {
        m_resultLo->setValue(e.first());
        m_resultHi->setValue(e.last());
    }
    const bool hybrid = m_backgroundMode->currentIndex() != 0;
    if (hybrid && (m_analysisSource->currentIndex() != StitchingSpectrum || m_wide.crystals.isEmpty())) {
        m_analysisInputInfo->setText(tr("康普顿法需要在第 6 步载入宽谱，并在这里选择第 6 步的谱。"));
    } else if (target.crystals.isEmpty()) {
        m_analysisInputInfo->setText(tr("请先在第 5 步生成原始谱、载入外部逐晶体谱，或选择第 6 步的谱。"));
    } else {
        m_analysisInputInfo->setText(tr("%1 个晶体 · %2 个能量点 · %3 至 %4 eV")
                                         .arg(target.crystals.size())
                                         .arg(e.size())
                                         .arg(formatNumber(e.first()), formatNumber(e.last())));
    }
}

void SpectrumProcessingPage::refreshAnalysis() {
    const auto label = m_fitCrystal->currentText();
    const auto &target = analysisInput();
    const auto geo = target.geometry.value(label);
    auto item = m_analysis.crystals.constFind(label);
    if (item == m_analysis.crystals.cend()) {
        m_fitPlot->setCurves({}, tr("晶体扣除诊断"));
        m_fitInfo->setText(
            tr("%1 · q %2 Å⁻¹ · %3")
                .arg(label, formatNumber(geo.qMean), m_analysis.failures.value(label, tr("等待分析"))));
    } else {
        const auto &e = target.energyEv;
        QVector<SpectrumCurve> curves{{tr("输入谱"), e, item->input, QColor("#0072b2")},
                                      {tr("扣除项"), e, item->fit.background, QColor("#d55e00"), true},
                                      {tr("扣除后"), e, item->fit.corrected, QColor("#009e73")}};
        if (item->high) {
            curves.append({tr("价电子轮廓"), e, item->valence, QColor("#895fa2"), true});
            curves.append({tr("HF 芯谱参考"), e, item->hfCore, QColor("#858585"), true});
        }
        m_fitPlot->setCurves(curves, label, "Energy transfer (eV)",
                             item->high ? "HF-normalized S (1/eV)" : target.intensityUnit);
        auto info = tr("%1 · q %2 Å⁻¹ · %3 · 拟合 R² %4")
                        .arg(label, formatNumber(geo.qMean), item->high ? tr("康普顿法") : tr("线性法"),
                             formatNumber(item->fit.rSquared));
        if (item->high)
            info += tr(" · HF 归一化系数 %1").arg(formatNumber(item->scale));
        if (m_analysis.referenceFailures.contains(label))
            info += tr(" · 未参与参考轮廓 %1").arg(m_analysis.referenceFailures.value(label));
        m_fitInfo->setText(info);
    }
    QVector<SpectrumCurve> reference;
    if (!m_analysis.measuredJ.isEmpty())
        reference.append({tr("高 q 中位轮廓"), m_analysis.pz, m_analysis.measuredJ, QColor("#777777")});
    if (!m_analysis.fittedJ.isEmpty())
        reference.append(
            {tr("Pearson VII 价电子轮廓"), m_analysis.pz, m_analysis.fittedJ, QColor("#0072b2")});
    m_referencePlot->setCurves(reference,
                               tr("高 q 参考 · %1 个晶体 · Pearson R² %2")
                                   .arg(m_analysis.references.size())
                                   .arg(formatNumber(m_analysis.referenceFit.rSquared)),
                               "pz (a.u.)", "q S (a.u.)");
}

QStringList SpectrumProcessingPage::resultLabels() const {
    auto names = m_resultChoice->currentData().toStringList();
    if (names.size() < 2 || m_resultHi->value() <= m_resultLo->value())
        return names;
    const auto &axis = analysisInput().energyEv;
    const auto first = std::lower_bound(axis.cbegin(), axis.cend(), m_resultLo->value()) - axis.cbegin();
    const auto last = std::upper_bound(axis.cbegin(), axis.cend(), m_resultHi->value()) - axis.cbegin();
    names.removeIf([&](const QString &label) {
        const auto &row = m_analysis.crystals.constFind(label)->fit.corrected;
        return std::any_of(row.cbegin() + first, row.cbegin() + last,
                           [](double value) { return !std::isfinite(value); });
    });
    return names;
}

void SpectrumProcessingPage::refreshResult() {
    const auto labels = resultLabels();
    const double lo = m_resultLo->value(), hi = m_resultHi->value();
    const auto &target = analysisInput();
    const auto &axis = target.energyEv;
    QVector<double> x, y;
    if (hi > lo && !labels.isEmpty()) {
        for (qsizetype i = 0; i < axis.size(); ++i) {
            if (axis[i] < lo || axis[i] > hi)
                continue;
            double sum = 0;
            int count = 0;
            for (const auto &label : labels) {
                const auto &row = m_analysis.crystals.constFind(label).value().fit.corrected;
                if (std::isfinite(row[i])) {
                    sum += row[i];
                    ++count;
                }
            }
            x.append(axis[i]);
            y.append(count ? sum / count : NAN);
        }
    }
    m_resultPlot->setCurves(x.isEmpty() ? QVector<SpectrumCurve>{}
                                        : QVector<SpectrumCurve>{{tr("扣背底结果"), x, y, QColor("#009e73")}},
                            m_resultChoice->currentText(), "Energy transfer (eV)",
                            !labels.isEmpty() && m_analysis.crystals.value(labels.first()).high
                                ? "HF-normalized S (1/eV)"
                                : target.intensityUnit);
    qsizetype peak = -1, finite = 0;
    for (qsizetype i = 0; i < y.size(); ++i)
        if (std::isfinite(y[i])) {
            ++finite;
            if (peak < 0 || y[i] > y[peak])
                peak = i;
        }
    m_resultInfo->setText(peak < 0 ? tr("当前能区没有有效结果，请运行扣背底或调整能量范围。")
                                   : tr("%1 个晶体 · %2 个有效能量点 · 最大采样点 %3 eV，强度 %4")
                                         .arg(labels.size())
                                         .arg(finite)
                                         .arg(formatNumber(x[peak]), formatNumber(y[peak])));
    m_saveResult->setEnabled(!m_busy && finite >= 2);
}

void SpectrumProcessingPage::refreshButtons() {
    const bool ready = !m_wide.crystals.isEmpty();
    m_useCurrent->setEnabled(!m_busy && !m_current.crystals.isEmpty());
    m_edges->setEnabled(!m_busy && ready);
    m_geometry->setEnabled(!m_busy && ready);
    m_saveStitch->setEnabled(!m_busy && !m_stitched.pointSources.isEmpty());
    const bool hybrid = m_backgroundMode->currentIndex() != 0;
    const bool inputReady = !analysisInput().crystals.isEmpty();
    m_run->setEnabled(!m_busy && inputReady &&
                      (!hybrid || (m_analysisSource->currentIndex() == StitchingSpectrum && ready)));
}
