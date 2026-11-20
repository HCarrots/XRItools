#include "SpectralMaskDialog.h"
#include "ScientificPlot.h"
#include "SpectrumPlot.h"
#include <QtConcurrent/QtConcurrentRun>
#include <QtWidgets>
#include <memory>
#include <numeric>

SpectralMaskDialog::SpectralMaskDialog(const SpectrumSettings &s, const QVector<CrystalRoi> &items,
                                       const QMap<QString, QSize> &dimensions, int initialRow,
                                       const QStringList &references, QWidget *parent)
    : QDialog(parent), crystals(items), sizes(dimensions), settings(s) {
    setObjectName("SpectralMaskDialog");
    resize(1150, 700);
    auto *root = new QVBoxLayout(this);
    auto *body = new QHBoxLayout;
    root->addLayout(body, 1);
    auto *sidebar = new QVBoxLayout;
    body->addLayout(sidebar);
    sidebar->addWidget(new QLabel(tr("晶体")));
    crystal = new QListWidget;
    crystal->setObjectName("SpectralCrystal");
    crystal->setMinimumWidth(180);
    sidebar->addWidget(crystal, 1);
    reviewProgress = new QLabel;
    sidebar->addWidget(reviewProgress);
    auto *layout = new QVBoxLayout;
    body->addLayout(layout, 5);
    auto *navigation = new QHBoxLayout;
    previous = new QPushButton(tr("上一个晶体"));
    next = new QPushButton(tr("下一个晶体"));
    previous->setObjectName("PreviousSpectralCrystal");
    next->setObjectName("NextSpectralCrystal");
    previous->setToolTip(QKeySequence("Alt+Left").toString(QKeySequence::NativeText));
    next->setToolTip(QKeySequence("Alt+Right").toString(QKeySequence::NativeText));
    new QShortcut(QKeySequence("Alt+Left"), this, [this] { previous->click(); });
    new QShortcut(QKeySequence("Alt+Right"), this, [this] { next->click(); });
    int initialIndex = 0;
    for (int r = 0; r < crystals.size(); ++r) {
        if (crystals[r].state != RoiState::Valid)
            continue;
        if (r == initialRow)
            initialIndex = crystal->count();
        auto *item = new QListWidgetItem(crystal);
        item->setData(Qt::UserRole, r);
    }
    position = new QLabel;
    navigation->addWidget(previous);
    navigation->addWidget(next);
    navigation->addWidget(position);
    navigation->addStretch();
    layout->addLayout(navigation);
    auto *sources = new QHBoxLayout;
    sources->addWidget(new QLabel(tr("参考 XRS 扫描")));
    reference = new QComboBox;
    reference->setObjectName("SpectralReference");
    reference->setEditable(true);
    reference->addItems(references);
    sources->addWidget(reference, 1);
    auto *browse = new QPushButton(tr("浏览…"));
    sources->addWidget(browse);
    layout->addLayout(sources);
    auto *controls = new QHBoxLayout;
    method = new QComboBox;
    method->addItems({"NNMF", "PCA", "ICA"});
    count = new QSpinBox;
    count->setRange(1, 16);
    count->setValue(2);
    controls->addWidget(new QLabel(tr("分解方法")));
    controls->addWidget(method);
    controls->addWidget(new QLabel(tr("组分数")));
    controls->addWidget(count);
    auto *calculateButton = new QPushButton(tr("从初始 ROI 分解"));
    calculateButton->setObjectName("DecomposePixels");
    controls->addWidget(calculateButton);
    iterate = new QPushButton(tr("用保留像素再分解"));
    iterate->setObjectName("IterateSpectralMask");
    controls->addWidget(iterate);
    controls->addStretch();
    layout->addLayout(controls);
    auto *selectors = new QHBoxLayout;
    component = new QComboBox;
    component->setObjectName("SampleComponent");
    selectors->addWidget(new QLabel(tr("目标组分")));
    selectors->addWidget(component);
    flip = new QCheckBox(tr("反转组分符号"));
    selectors->addWidget(flip);
    thresholdMode = new QComboBox;
    thresholdMode->setObjectName("CovarianceThresholdMode");
    thresholdMode->addItems({tr("协方差阈值 / 正峰值"), tr("协方差绝对阈值")});
    selectors->addWidget(thresholdMode);
    threshold = new QDoubleSpinBox;
    threshold->setObjectName("CovarianceThreshold");
    threshold->setDecimals(5);
    threshold->setRange(0, 1);
    threshold->setSingleStep(.05);
    threshold->setValue(.2);
    selectors->addWidget(threshold);
    selectors->addStretch();
    layout->addLayout(selectors);
    auto *panels = new QSplitter;
    curves = new SpectrumPlot(nullptr, SpectrumPlot::Style::Spectroscopy);
    curves->setObjectName("SpectralComponents");
    panels->addWidget(curves);
    auto *views = new QTabWidget;
    map = new ScientificPlot;
    map->setAxisLabels("Detector X", "Detector Y");
    scores = new SpectrumPlot;
    comparison = new SpectrumPlot(nullptr, SpectrumPlot::Style::Spectroscopy);
    views->addTab(map, tr("候选与保留像素"));
    views->addTab(scores, tr("逐像素协方差"));
    views->addTab(comparison, tr("积分谱对照"));
    panels->addWidget(views);
    layout->addWidget(panels, 1);
    status = new QLabel;
    status->setObjectName("SpectralMaskStatus");
    status->setWordWrap(true);
    layout->addWidget(status);
    auto *buttons = new QDialogButtonBox(QDialogButtonBox::Close);
    buttons->button(QDialogButtonBox::Close)->setText(tr("关闭"));
    apply = buttons->addButton(tr("采用所选组分"), QDialogButtonBox::ActionRole);
    apply->setObjectName("ApplySpectralMask");
    keepOriginal = buttons->addButton(tr("不采用，保留原 ROI"), QDialogButtonBox::ActionRole);
    keepOriginal->setObjectName("KeepOriginalRoi");
    updateRoi = buttons->addButton(tr("更新原始 ROI"), QDialogButtonBox::ActionRole);
    updateRoi->setObjectName("UpdateOriginalRoi");
    updateRoi->setToolTip(tr("全部晶体确认后，一次性更新已采用的 ROI。"));
    root->addWidget(buttons);
    connect(updateRoi, &QPushButton::clicked, this, [this] {
        if (updateRoi->isEnabled())
            accept();
    });
    connect(apply, &QPushButton::clicked, this, [this] { recordDecision(true); });
    connect(keepOriginal, &QPushButton::clicked, this, [this] { recordDecision(false); });
    connect(buttons, &QDialogButtonBox::rejected, this, &QDialog::reject);
    connect(browse, &QPushButton::clicked, this, [this] {
        const auto file = QFileDialog::getOpenFileName(this, tr("参考 XRS 扫描"), reference->currentText(),
                                                       "NeXus / HDF5 (*.nxs *.h5 *.hdf5);;All (*)");
        if (!file.isEmpty())
            reference->setCurrentText(file);
    });
    connect(reference, &QComboBox::currentTextChanged, this, [this] {
        autoDecompose = false;
        data = {};
        resetDecisions();
    });
    connect(method, &QComboBox::currentIndexChanged, this, &SpectralMaskDialog::resetDecisions);
    connect(count, &QSpinBox::valueChanged, this, &SpectralMaskDialog::resetDecisions);
    connect(calculateButton, &QPushButton::clicked, this, [this] { calculate(false); });
    connect(iterate, &QPushButton::clicked, this, [this] { calculate(true); });
    connect(component, &QComboBox::currentIndexChanged, this, &SpectralMaskDialog::editPreview);
    connect(flip, &QCheckBox::toggled, this, &SpectralMaskDialog::editPreview);
    connect(threshold, &QDoubleSpinBox::valueChanged, this, &SpectralMaskDialog::editPreview);
    connect(thresholdMode, &QComboBox::currentIndexChanged, this, [this](int mode) {
        threshold->setDecimals(mode ? 10 : 5);
        threshold->setRange(mode ? -1e30 : 0, mode ? 1e30 : 1);
        threshold->setValue(mode ? 0 : .2);
        editPreview();
    });
    crystal->setCurrentRow(initialIndex);
    connect(crystal, &QListWidget::currentRowChanged, this, &SpectralMaskDialog::switchCrystal);
    connect(previous, &QPushButton::clicked, this,
            [this] { crystal->setCurrentRow(crystal->currentRow() - 1); });
    connect(next, &QPushButton::clicked, this,
            [this] { crystal->setCurrentRow(crystal->currentRow() + 1); });
    switchCrystal(crystal->currentRow());
    refreshDecisions();
}

void SpectralMaskDialog::keyPressEvent(QKeyEvent *event) {
    if (event->modifiers() == Qt::AltModifier &&
        (event->key() == Qt::Key_Left || event->key() == Qt::Key_Right)) {
        auto *button = event->key() == Qt::Key_Left ? previous : next;
        button->click();
        event->accept();
        return;
    }
    QDialog::keyPressEvent(event);
}

void SpectralMaskDialog::switchCrystal(int index) {
    if (index < 0)
        return;
    QScopedValueRollback restoring(switching, true);
    if (row >= 0 && !result.spectra.isEmpty()) {
        const qsizetype bytes =
            data.pixels.size() *
                (sizeof(QPoint) + sizeof(QVector<double>) + sizeof(double) * data.energyKev.size()) +
            result.spectra.size() * (rows.size() + data.energyKev.size()) * sizeof(double);
        reviews.insert(row,
                       new Review{std::move(data), std::move(result), std::move(rows), std::move(history),
                                  round, component->currentIndex(), thresholdMode->currentIndex(),
                                  threshold->value(), flip->isChecked()},
                       1 + bytes / 1024);
    }
    row = crystal->item(index)->data(Qt::UserRole).toInt();
    roi = crystals[row];
    size = sizes.value(roi.detector);
    data = {};
    clearResult();
    previous->setEnabled(index > 0);
    next->setEnabled(index + 1 < crystal->count());
    position->setText(QStringLiteral("%1 / %2 · %3").arg(index + 1).arg(crystal->count()).arg(roi.detector));
    setWindowTitle(tr("掩膜与谱分辨 · %1").arg(roi.label));
    if (std::unique_ptr<Review> review{reviews.take(row)}) {
        thresholdMode->setCurrentIndex(review->thresholdMode);
        threshold->setValue(review->threshold);
        flip->setChecked(review->flip);
        data = std::move(review->data);
        result = std::move(review->result);
        rows = std::move(review->rows);
        history = std::move(review->history);
        round = review->round;
        QSignalBlocker block(component);
        for (int i = 0; i < result.spectra.size(); ++i)
            component->addItem(QString::number(i + 1));
        component->setCurrentIndex(review->component);
        preview();
    } else if (accepted.contains(row)) {
        auto image = map->data();
        for (const QPoint &p : accepted.value(row).pixels)
            image.values[(p.y() - roi.rect.top()) * image.width + p.x() - roi.rect.left()] = 2;
        image.maximum = 2;
        image.path = tr("候选 = 1，保留 = 2");
        map->setData(image);
        map->setImageCoordinates(roi.rect.left(), 1, roi.rect.top(), 1);
        status->setText(tr("已记录此晶体的组分掩膜。可重新分解以修改选择。"));
    } else if (autoDecompose && !declined.contains(row)) {
        calculate(false);
    }
    refreshDecisions();
}

void SpectralMaskDialog::refreshDecisions() {
    for (int i = 0; i < crystal->count(); ++i) {
        auto *item = crystal->item(i);
        const int r = item->data(Qt::UserRole).toInt();
        QString state = computed.contains(r) ? tr("待确认") : tr("未计算");
        if (declined.contains(r))
            state = tr("保留原 ROI");
        else if (auto it = accepted.constFind(r); it != accepted.cend()) {
            const auto last = it->provenance.value("rounds").toArray().last().toObject();
            state = tr("采用组分 %1").arg(last.value("target_component_1based").toInt());
        }
        item->setText(crystals[r].label + QStringLiteral("  ·  ") + state);
    }
    const int confirmed = accepted.size() + declined.size();
    reviewProgress->setText(tr("已确认 %1 / %2").arg(confirmed).arg(crystal->count()));
    updateRoi->setEnabled(confirmed > 0 && confirmed == crystal->count());
    keepOriginal->setEnabled(computed.contains(row));
}

void SpectralMaskDialog::recordDecision(bool useComponent) {
    if (!computed.contains(row) || (useComponent && selected.isEmpty()))
        return;
    if (useComponent) {
        accepted.insert(row, selection());
        declined.remove(row);
    } else {
        accepted.remove(row);
        declined.insert(row);
    }
    refreshDecisions();
    const int current = crystal->currentRow();
    for (int offset = 1; offset < crystal->count(); ++offset) {
        const int i = (current + offset) % crystal->count();
        const int r = crystal->item(i)->data(Qt::UserRole).toInt();
        if (!accepted.contains(r) && !declined.contains(r)) {
            crystal->setCurrentRow(i);
            break;
        }
    }
}

void SpectralMaskDialog::resetDecisions() {
    reviews.clear();
    computed.clear();
    accepted.clear();
    declined.clear();
    clearResult();
    refreshDecisions();
}

void SpectralMaskDialog::editPreview() {
    if (!switching && !result.spectra.isEmpty()) {
        accepted.remove(row);
        declined.remove(row);
        refreshDecisions();
    }
    preview();
}

void SpectralMaskDialog::clearResult() {
    result = {};
    round = 0;
    history = {};
    rows.clear();
    selected.clear();
    component->clear();
    apply->setEnabled(false);
    iterate->setEnabled(false);
    flip->setChecked(false);
    flip->setEnabled(method->currentIndex() != 0);
    curves->setCurves({}, {});
    scores->setCurves({}, {});
    comparison->setCurves({}, {});
    DataSlice image;
    if (roi.state == RoiState::Valid) {
        image.path = tr("当前原始 ROI");
        image.width = roi.rect.width();
        image.height = roi.rect.height();
        image.dims = {quint64(image.height), quint64(image.width)};
        image.xAxis = 1;
        image.yAxis = 0;
        image.starts = {0, 0};
        image.strides = {1, 1};
        image.values.fill(0, image.width * image.height);
        for (const auto &p : roi.pixels())
            image.values[(p.y() - roi.rect.top()) * image.width + p.x() - roi.rect.left()] = 1;
        image.minimum = 0;
        image.maximum = 1;
        image.finiteCount = image.values.size();
    }
    map->setData(image);
    map->setImageCoordinates(roi.rect.left(), 1, roi.rect.top(), 1);
    status->setText(tr("选择具有样品特征的参考扫描，分解后根据谱形选择目标组分。"));
}

void SpectralMaskDialog::calculate(bool nextRound) {
    const QString file = reference->currentText().trimmed();
    const int k = count->value();
    const auto algorithm = static_cast<UnmixMethod>(method->currentIndex());
    QVector<int> candidates = nextRound ? selected : QVector<int>{};
    SpectralComponents output;
    QString error;
    bool noComponent = false;
    std::atomic_bool cancel{false};
    QProgressDialog progress(tr("正在读取像素谱并分解…"), tr("取消"), 0, 0, this);
    progress.setWindowModality(Qt::WindowModal);
    QFutureWatcher<void> watcher;
    connect(&watcher, &QFutureWatcher<void>::finished, &progress, &QDialog::accept);
    connect(&progress, &QProgressDialog::canceled, this, [&cancel] { cancel.store(true); });
    watcher.setFuture(QtConcurrent::run([&] {
        try {
            if (data.energyKev.isEmpty())
                data = Spectroscopy::pixelSpectra(file, settings, roi, size, cancel);
            if (!nextRound) {
                candidates.resize(data.pixels.size());
                std::iota(candidates.begin(), candidates.end(), 0);
            }
            output = decomposePixels(data, candidates, algorithm, k, cancel);
        } catch (const NoSpectralComponent &e) {
            noComponent = true;
            error = QString::fromUtf8(e.what());
        } catch (const std::exception &e) {
            error = QString::fromUtf8(e.what());
        }
    }));
    if (!watcher.isFinished())
        progress.exec();
    if (!watcher.isFinished())
        cancel.store(true);
    watcher.waitForFinished();
    accepted.remove(row);
    declined.remove(row);
    if (cancel.load() || !error.isEmpty()) {
        clearResult();
        if (noComponent && !cancel.load()) {
            computed.insert(row);
            autoDecompose = true;
        } else {
            computed.remove(row);
        }
        refreshDecisions();
        status->setText(cancel.load() ? tr("已取消。") : error);
        return;
    }
    if (nextRound)
        history.append(roundRecord());
    else
        history = {};
    rows = std::move(candidates);
    autoDecompose = true;
    computed.insert(row);
    result = std::move(output);
    round = nextRound ? round + 1 : 1;
    {
        QSignalBlocker block(component);
        component->clear();
        for (int i = 0; i < k; ++i)
            component->addItem(QString::number(i + 1));
    }
    flip->setChecked(false);
    preview();
    refreshDecisions();
}

void SpectralMaskDialog::preview() {
    const int k = component->currentIndex();
    if (k < 0 || k >= result.spectra.size())
        return;
    const double sign = flip->isChecked() ? -1 : 1;
    QVector<double> covariance = result.covariance[k];
    for (auto &v : covariance)
        v *= sign;
    const double cutoff =
        thresholdMode->currentIndex()
            ? threshold->value()
            : threshold->value() * std::max(0.0, *std::max_element(covariance.begin(), covariance.end()));
    cutoffValue = cutoff;
    QVector<SpectrumCurve> components;
    for (int i = 0; i < result.spectra.size(); ++i) {
        auto shape = result.spectra[i];
        if (i == k && sign < 0)
            for (auto &v : shape)
                v = -v;
        components.append({QString::number(i + 1), data.energyKev, shape, {}, i != k, SpectrumCurve::Draw::Line});
    }
    curves->setCurves(components, tr("组分谱 · 所选组分为实线"), tr("入射能量 (keV)"), tr("单位范数组分"));
    DataSlice image;
    image.width = roi.rect.width();
    image.height = roi.rect.height();
    image.dims = {quint64(image.height), quint64(image.width)};
    image.xAxis = 1;
    image.yAxis = 0;
    image.starts = {0, 0};
    image.strides = {1, 1};
    image.minimum = 0;
    image.maximum = 2;
    image.values.fill(0, image.width * image.height);
    image.finiteCount = image.values.size();
    selected.clear();
    QVector<double> all(data.energyKev.size(), 0), kept(all), indices(rows.size());
    for (int i = 0; i < rows.size(); ++i) {
        const int p = rows[i];
        const bool keep = covariance[i] > cutoff;
        if (keep)
            selected.append(p);
        const QPoint local = data.pixels[p] - roi.rect.topLeft();
        image.values[local.y() * image.width + local.x()] = keep ? 2 : 1;
        indices[i] = p + 1;
        for (int j = 0; j < all.size(); ++j) {
            all[j] += data.intensity[p][j];
            if (keep)
                kept[j] += data.intensity[p][j];
        }
    }
    image.path = tr("候选 = 1，保留 = 2");
    map->setData(image);
    map->setImageCoordinates(roi.rect.left(), 1, roi.rect.top(), 1);
    scores->setCurves({{tr("协方差"), indices, covariance, {}},
                       {tr("阈值"), indices, QVector<double>(indices.size(), cutoff), Qt::red, true}},
                      roi.label, tr("初始 ROI 像素序号"), tr("协方差"));
    comparison->setCurves(
        {{tr("本轮候选"), data.energyKev, all, {}}, {tr("保留像素"), data.energyKev, kept, Qt::darkGreen}},
        roi.label, tr("入射能量 (keV)"), settings.normalizeXrs ? "Counts / |I0|" : "Counts");
    status->setText(tr("第 %1 轮 · 保留 %2 / %3 像素 · 协方差阈值 %4 · 相对重构残差 %5")
                        .arg(round)
                        .arg(selected.size())
                        .arg(rows.size())
                        .arg(cutoff, 0, 'g', 5)
                        .arg(result.relativeResidual, 0, 'g', 4) +
                    (result.converged ? QString() : tr(" · 达到迭代上限")));
    apply->setEnabled(!selected.isEmpty());
    iterate->setEnabled(selected.size() >= count->value() && selected.size() < rows.size());
}

QJsonObject SpectralMaskDialog::roundRecord() const {
    return {{"method", method->currentText()},
            {"components", result.spectra.size()},
            {"target_component_1based", component->currentIndex() + 1},
            {"sign", flip->isChecked() ? -1 : 1},
            {"threshold_mode", thresholdMode->currentIndex() ? "absolute" : "fraction_positive_max"},
            {"threshold_setting", threshold->value()},
            {"covariance_cutoff", cutoffValue},
            {"candidates", rows.size()},
            {"retained", selected.size()},
            {"relative_residual", result.relativeResidual},
            {"converged", result.converged}};
}

SpectralSelection SpectralMaskDialog::selection() const {
    QVector<QPoint> points;
    for (int p : selected)
        points.append(data.pixels[p]);
    auto rounds = history;
    rounds.append(roundRecord());
    return {points,
            {{"reference_scan", QFileInfo(reference->currentText().trimmed()).absoluteFilePath()},
             {"energy_channel", settings.xrsEnergy},
             {"normalize_i0", settings.normalizeXrs},
             {"i0_channel", settings.xrsI0},
             {"repair_glitches", settings.repairGlitches},
             {"glitch_ratio", 1 - settings.glitchDropFraction},
             {"rounds", rounds}}};
}
