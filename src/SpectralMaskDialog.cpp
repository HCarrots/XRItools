#include "SpectralMaskDialog.h"
#include "ScientificPlot.h"
#include "SpectrumPlot.h"
#include <QDataStream>
#include <QSaveFile>
#include <QtConcurrent/QtConcurrentRun>
#include <QtWidgets>
#include <algorithm>
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
    decomposeAll = new QPushButton(tr("分解全部晶体"));
    decomposeAll->setObjectName("DecomposeAllSpectralMasks");
    controls->addWidget(decomposeAll);
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
    batchApply = buttons->addButton(tr("批量采用所选组分"), QDialogButtonBox::ActionRole);
    batchApply->setObjectName("ApplyBatchSpectralMask");
    batchApply->setToolTip(tr("对待确认晶体使用相同的组分编号、阈值和符号。组分编号对应各晶体自己的分解结果。"));
    keepOriginal = buttons->addButton(tr("不采用，保留原 ROI"), QDialogButtonBox::ActionRole);
    keepOriginal->setObjectName("KeepOriginalRoi");
    updateRoi = buttons->addButton(tr("应用谱分解掩膜"), QDialogButtonBox::ActionRole);
    updateRoi->setObjectName("UpdateOriginalRoi");
    updateRoi->setToolTip(tr("全部晶体确认后，一次性应用谱分解掩膜。"));
    root->addWidget(buttons);
    connect(updateRoi, &QPushButton::clicked, this, [this] {
        if (updateRoi->isEnabled())
            accept();
    });
    connect(apply, &QPushButton::clicked, this, [this] { recordDecision(true); });
    connect(batchApply, &QPushButton::clicked, this, &SpectralMaskDialog::applyBatch);
    connect(keepOriginal, &QPushButton::clicked, this, [this] { recordDecision(false); });
    connect(buttons, &QDialogButtonBox::rejected, this, &QDialog::reject);
    connect(browse, &QPushButton::clicked, this, [this] {
        const auto file = QFileDialog::getOpenFileName(this, tr("参考 XRS 扫描"), reference->currentText(),
                                                       "NeXus / HDF5 (*.nxs *.h5 *.hdf5);;All (*)");
        if (!file.isEmpty())
            reference->setCurrentText(file);
    });
    connect(reference, &QComboBox::currentTextChanged, this, &SpectralMaskDialog::resetDecisions);
    connect(method, &QComboBox::currentIndexChanged, this, &SpectralMaskDialog::resetDecisions);
    connect(count, &QSpinBox::valueChanged, this, &SpectralMaskDialog::resetDecisions);
    connect(calculateButton, &QPushButton::clicked, this, [this] { calculate(false); });
    connect(decomposeAll, &QPushButton::clicked, this, &SpectralMaskDialog::calculateAll);
    connect(iterate, &QPushButton::clicked, this, [this] { calculate(true); });
    connect(component, &QComboBox::currentIndexChanged, this, &SpectralMaskDialog::editPreview);
    connect(flip, &QCheckBox::toggled, this, &SpectralMaskDialog::editPreview);
    connect(threshold, &QDoubleSpinBox::valueChanged, this, &SpectralMaskDialog::editPreview);
    connect(thresholdMode, &QComboBox::currentIndexChanged, this, [this](int mode) {
        setThresholdMode(mode, mode ? 0. : .2);
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
    row = crystal->item(index)->data(Qt::UserRole).toInt();
    roi = crystals[row];
    previous->setEnabled(index > 0);
    next->setEnabled(index + 1 < crystal->count());
    position->setText(QStringLiteral("%1 / %2 · %3").arg(index + 1).arg(crystal->count()).arg(roi.detector));
    setWindowTitle(tr("掩膜与谱分辨 · %1").arg(roi.label));
    showReview();
    refreshDecisions();
}

void SpectralMaskDialog::showReview() {
    const QScopedValueRollback restoring(switching, true);
    intensity.clear();
    clearPreview();
    const auto review = reviews.constFind(row);
    const MaskOptions options = review == reviews.cend() ? MaskOptions{} : review->options;
    {
        const QSignalBlocker blockMode(thresholdMode);
        thresholdMode->setCurrentIndex(options.thresholdMode);
        setThresholdMode(options.thresholdMode, options.threshold);
    }
    flip->setChecked(options.flip);
    flip->setEnabled(method->currentIndex() != 0);
    if (review == reviews.cend()) {
        status->setText(tr("请选择参考扫描，然后分解当前晶体或全部晶体。"));
        return;
    }
    if (review->result.spectra.isEmpty()) {
        status->setText(review->error);
        return;
    }
    try {
        intensity = loadIntensity(row);
    } catch (const std::exception &error) {
        status->setText(QString::fromUtf8(error.what()));
        return;
    }
    for (qsizetype i = 0; i < review->result.spectra.size(); ++i)
        component->addItem(QString::number(i + 1));
    component->setCurrentIndex(options.component);
    preview();
}

void SpectralMaskDialog::refreshDecisions() {
    bool pendingDecomposition = false;
    const bool anyResult = std::any_of(reviews.cbegin(), reviews.cend(),
                                       [](const Review &review) { return !review.result.spectra.isEmpty(); });
    for (int i = 0; i < crystal->count(); ++i) {
        auto *item = crystal->item(i);
        const int r = item->data(Qt::UserRole).toInt();
        const auto review = reviews.constFind(r);
        const bool ready = review != reviews.cend() && (!review->result.spectra.isEmpty() || review->noComponent);
        const bool decided = accepted.contains(r) || declined.contains(r);
        pendingDecomposition |= !ready && !decided;
        QString state = tr("未计算");
        if (declined.contains(r))
            state = tr("保留原 ROI");
        else if (auto it = accepted.constFind(r); it != accepted.cend()) {
            const auto last = it->provenance.value("rounds").toArray().last().toObject();
            state = tr("采用组分 %1").arg(last.value("target_component_1based").toInt());
        } else if (ready) {
            state = tr("待确认");
        } else if (review != reviews.cend() && !review->error.isEmpty()) {
            state = tr("计算失败");
        }
        item->setText(crystals[r].label + QStringLiteral("  ·  ") + state);
    }
    const int confirmed = accepted.size() + declined.size();
    reviewProgress->setText(tr("已确认 %1 / %2").arg(confirmed).arg(crystal->count()));
    updateRoi->setEnabled(confirmed > 0 && confirmed == crystal->count());
    decomposeAll->setEnabled(pendingDecomposition && !reference->currentText().trimmed().isEmpty());
    batchApply->setEnabled(component->currentIndex() >= 0 && confirmed < crystal->count());
    const auto current = reviews.constFind(row);
    keepOriginal->setEnabled(current != reviews.cend() &&
                            (current->noComponent || !current->result.spectra.isEmpty() || anyResult));
}

void SpectralMaskDialog::recordDecision(bool useComponent) {
    if (useComponent ? currentMask.selected.isEmpty() : !keepOriginal->isEnabled())
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

void SpectralMaskDialog::applyBatch() {
    const auto options = maskOptions();
    if (options.component < 0 || !calculateAll())
        return;
    int applied = 0, kept = 0, pending = 0;
    for (auto it = reviews.begin(); it != reviews.end(); ++it) {
        const int r = it.key();
        if (accepted.contains(r) || declined.contains(r))
            continue;
        if (it->noComponent) {
            declined.insert(r);
            ++kept;
            continue;
        }
        if (options.component >= it->result.spectra.size()) {
            ++pending;
            continue;
        }
        it->options = options;
        const auto mask = selectPixels(*it, options);
        if (mask.selected.isEmpty()) {
            ++pending;
            continue;
        }
        accepted.insert(r, selection(*it, options, mask));
        ++applied;
    }
    refreshDecisions();
    status->setText(tr("组分 %1 · 已采用 %2 个 · 无可用组分，保留原 ROI %3 个 · 待确认 %4 个")
                        .arg(options.component + 1).arg(applied).arg(kept).arg(pending));
}

void SpectralMaskDialog::resetDecisions() {
    // These files are private to the current dialog and tied to the reference and fit settings.
    for (auto it = reviews.cbegin(); it != reviews.cend(); ++it)
        QFile::remove(pixelCachePath(it.key()));
    reviews.clear();
    accepted.clear();
    declined.clear();
    showReview();
    refreshDecisions();
}

SpectralMaskDialog::MaskOptions SpectralMaskDialog::maskOptions() const {
    return {component->currentIndex(), thresholdMode->currentIndex(), threshold->value(), flip->isChecked()};
}

void SpectralMaskDialog::setThresholdMode(int mode, double value) {
    const QSignalBlocker block(threshold);
    threshold->setDecimals(mode ? 10 : 5);
    threshold->setRange(mode ? -1e30 : 0, mode ? 1e30 : 1);
    threshold->setValue(value);
}

void SpectralMaskDialog::editPreview() {
    if (switching || component->currentIndex() < 0)
        return;
    reviews[row].options = maskOptions();
    accepted.remove(row);
    declined.remove(row);
    preview();
    refreshDecisions();
}

void SpectralMaskDialog::clearPreview() {
    currentMask = {};
    component->clear();
    apply->setEnabled(false);
    iterate->setEnabled(false);
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
}

QString SpectralMaskDialog::pixelCachePath(int crystalRow) const {
    return pixelCache.filePath(QString::number(crystalRow) + QStringLiteral(".bin"));
}

void SpectralMaskDialog::cacheIntensity(int crystalRow, const QVector<QVector<double>> &values) const {
    if (!pixelCache.isValid())
        throw std::runtime_error(tr("无法创建像素谱缓存。").toStdString());
    QSaveFile file(pixelCachePath(crystalRow));
    if (!file.open(QIODevice::WriteOnly))
        throw std::runtime_error(file.errorString().toStdString());
    QDataStream stream(&file);
    stream << values;
    if (stream.status() != QDataStream::Ok || !file.commit())
        throw std::runtime_error(tr("无法写入像素谱缓存。").toStdString());
}

QVector<QVector<double>> SpectralMaskDialog::loadIntensity(int crystalRow) const {
    QFile file(pixelCachePath(crystalRow));
    QVector<QVector<double>> values;
    if (file.open(QIODevice::ReadOnly)) {
        QDataStream stream(&file);
        stream >> values;
        if (stream.status() == QDataStream::Ok)
            return values;
    }
    throw std::runtime_error(tr("无法读取像素谱缓存，请重新分解当前晶体。").toStdString());
}

void SpectralMaskDialog::calculate(bool refine) {
    if (row >= 0)
        runDecomposition({row}, refine);
}

bool SpectralMaskDialog::calculateAll() {
    QVector<int> targets;
    for (int i = 0; i < crystal->count(); ++i) {
        const int r = crystal->item(i)->data(Qt::UserRole).toInt();
        const auto review = reviews.constFind(r);
        if (!declined.contains(r) && (review == reviews.cend() ||
                                     (review->result.spectra.isEmpty() && !review->noComponent)))
            targets.append(r);
    }
    return targets.isEmpty() || runDecomposition(targets);
}

bool SpectralMaskDialog::runDecomposition(const QVector<int> &targets, bool refine) {
    const QString file = reference->currentText().trimmed();
    const int k = count->value(), activeRow = row;
    const auto activeSize = sizes.value(roi.detector);
    const auto algorithm = static_cast<UnmixMethod>(method->currentIndex());
    auto options = maskOptions();
    options.component = 0;
    options.flip = false;
    const auto activeReview = reviews.value(activeRow);
    const auto initialRows = refine ? currentMask.selected : QVector<int>{};
    const auto activeHistory = refine ? selection().provenance.value("rounds").toArray() : QJsonArray{};
    const auto activeIntensity = intensity;
    QMap<int, Review> outcomes;
    std::atomic_bool cancel{false};
    std::atomic_int completed{0};
    QProgressDialog progress(tr("正在读取像素谱并分解…"), tr("取消"), 0, int(targets.size()), this);
    progress.setWindowModality(Qt::WindowModal);
    progress.setAutoClose(false);
    progress.setAutoReset(false);
    QTimer progressTimer;
    connect(&progressTimer, &QTimer::timeout, &progress, [&] { progress.setValue(completed.load()); });
    QFutureWatcher<void> watcher;
    connect(&watcher, &QFutureWatcher<void>::finished, &progress, &QDialog::accept);
    connect(&progress, &QProgressDialog::canceled, &progress, [&] { cancel.store(true); });
    watcher.setFuture(QtConcurrent::run([&] {
        auto store = [&](int r, Review &&review) {
            if (!cancel.load()) {
                outcomes.insert(r, std::move(review));
                ++completed;
            }
        };
        auto decompose = [&](int r, PixelSpectra pixelData, bool fromReference) {
            Review review;
            review.options = r == activeRow ? options : MaskOptions{};
            try {
                review.rows = initialRows;
                if (!refine) {
                    review.rows.resize(pixelData.pixels.size());
                    std::iota(review.rows.begin(), review.rows.end(), 0);
                }
                review.result = decomposePixels(pixelData, review.rows, algorithm, k, cancel);
                if (fromReference)
                    cacheIntensity(r, pixelData.intensity);
                review.pixels = std::move(pixelData.pixels);
                review.energyKev = std::move(pixelData.energyKev);
                review.history = activeHistory;
            } catch (const NoSpectralComponent &error) {
                review.noComponent = true;
                review.error = QString::fromUtf8(error.what());
            } catch (const std::exception &error) {
                review.result = {};
                review.error = QString::fromUtf8(error.what());
            }
            store(r, std::move(review));
        };
        if (targets.size() == 1 && targets.first() == activeRow && !activeIntensity.isEmpty()) {
            PixelSpectra pixelData{activeSize, activeReview.pixels, activeReview.energyKev, activeIntensity};
            decompose(activeRow, std::move(pixelData), false);
            return;
        }

        QMap<QString, QVector<int>> byDetector;
        for (const int r : targets)
            byDetector[crystals[r].detector].append(r);
        for (auto it = byDetector.cbegin(); it != byDetector.cend() && !cancel.load(); ++it) {
            const auto &indices = it.value();
            QVector<CrystalRoi> items;
            for (const int r : indices)
                items.append(crystals[r]);
            try {
                Spectroscopy::visitPixelSpectra(file, settings, items, sizes.value(it.key()), cancel,
                    [&](qsizetype index, PixelSpectra pixelData, const QString &error) {
                        if (error.isEmpty()) {
                            decompose(indices[index], std::move(pixelData), true);
                        } else {
                            Review review;
                            review.error = error;
                            store(indices[index], std::move(review));
                        }
                    });
            } catch (const std::exception &error) {
                for (const int r : indices)
                    if (!outcomes.contains(r)) {
                        Review review;
                        review.error = QString::fromUtf8(error.what());
                        store(r, std::move(review));
                    }
            }
        }
    }));
    if (!watcher.isFinished()) {
        progressTimer.start(50);
        progress.exec();
    }
    if (!watcher.isFinished())
        cancel.store(true);
    watcher.waitForFinished();
    QString activeError;
    for (auto it = outcomes.begin(); it != outcomes.end(); ++it) {
        const int r = it.key();
        if (r == activeRow)
            activeError = it->error;
        const auto prior = reviews.constFind(r);
        if (!it->error.isEmpty() && prior != reviews.cend() && !prior->result.spectra.isEmpty())
            continue;
        accepted.remove(r);
        declined.remove(r);
        reviews.insert(r, std::move(it.value()));
    }
    if (outcomes.contains(activeRow))
        showReview();
    refreshDecisions();
    if (cancel.load())
        status->setText(tr("已取消。"));
    else if (targets.size() == 1 && !activeError.isEmpty())
        status->setText(activeError);
    else if (targets.size() > 1) {
        qsizetype succeeded = 0, noComponent = 0;
        for (const int r : targets) {
            const auto review = reviews.constFind(r);
            if (review != reviews.cend()) {
                succeeded += !review->result.spectra.isEmpty();
                noComponent += review->noComponent;
            }
        }
        status->setText(tr("分解完成 · 成功 %1 个 · 无可用组分 %2 个 · 失败 %3 个")
                            .arg(succeeded).arg(noComponent).arg(targets.size() - succeeded - noComponent));
    }
    return !cancel.load();
}

SpectralMaskDialog::MaskPreview SpectralMaskDialog::selectPixels(const Review &review, const MaskOptions &options) {
    MaskPreview mask;
    mask.covariance = review.result.covariance[options.component];
    if (options.flip)
        for (auto &value : mask.covariance)
            value = -value;
    mask.cutoff = options.thresholdMode ? options.threshold :
        options.threshold * std::max(0.0, *std::max_element(mask.covariance.cbegin(), mask.covariance.cend()));
    for (qsizetype i = 0; i < review.rows.size(); ++i)
        if (mask.covariance[i] > mask.cutoff)
            mask.selected.append(review.rows[i]);
    return mask;
}

void SpectralMaskDialog::preview() {
    const auto options = maskOptions();
    const int k = options.component;
    const auto review = reviews.constFind(row);
    if (review == reviews.cend() || k < 0 || k >= review->result.spectra.size() || intensity.isEmpty())
        return;
    const auto &result = review->result;
    const auto &rows = review->rows;
    const auto &energy = review->energyKev;
    currentMask = selectPixels(*review, options);
    const auto &covariance = currentMask.covariance;
    const double cutoff = currentMask.cutoff;
    QVector<SpectrumCurve> components;
    for (int i = 0; i < result.spectra.size(); ++i) {
        auto shape = result.spectra[i];
        if (i == k && options.flip)
            for (auto &v : shape)
                v = -v;
        components.append({QString::number(i + 1), energy, shape, {}, i != k, SpectrumCurve::Draw::Line});
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
    QVector<double> all(energy.size(), 0), kept(all), indices(rows.size());
    for (int i = 0; i < rows.size(); ++i) {
        const int p = rows[i];
        const bool keep = covariance[i] > cutoff;
        const QPoint local = review->pixels[p] - roi.rect.topLeft();
        image.values[local.y() * image.width + local.x()] = keep ? 2 : 1;
        indices[i] = p + 1;
        for (int j = 0; j < all.size(); ++j) {
            all[j] += intensity[p][j];
            if (keep)
                kept[j] += intensity[p][j];
        }
    }
    image.path = tr("候选 = 1，保留 = 2");
    map->setData(image);
    map->setImageCoordinates(roi.rect.left(), 1, roi.rect.top(), 1);
    scores->setCurves({{tr("协方差"), indices, covariance, {}},
                       {tr("阈值"), indices, QVector<double>(indices.size(), cutoff), Qt::red, true}},
                      roi.label, tr("初始 ROI 像素序号"), tr("协方差"));
    comparison->setCurves(
        {{tr("本轮候选"), energy, all, {}}, {tr("保留像素"), energy, kept, Qt::darkGreen}},
        roi.label, tr("入射能量 (keV)"), settings.normalizeXrs ? "Counts / |I0|" : "Counts");
    status->setText(tr("第 %1 轮 · 保留 %2 / %3 像素 · 协方差阈值 %4 · 相对重构残差 %5")
                         .arg(review->history.size() + 1)
                         .arg(currentMask.selected.size())
                        .arg(rows.size())
                        .arg(cutoff, 0, 'g', 5)
                        .arg(result.relativeResidual, 0, 'g', 4) +
                    (result.converged ? QString() : tr(" · 达到迭代上限")));
    apply->setEnabled(!currentMask.selected.isEmpty());
    iterate->setEnabled(currentMask.selected.size() >= count->value() && currentMask.selected.size() < rows.size());
}

QJsonObject SpectralMaskDialog::roundRecord(const Review &review, const MaskOptions &options,
                                          const MaskPreview &mask) const {
    return {{"method", method->currentText()},
            {"components", review.result.spectra.size()},
            {"target_component_1based", options.component + 1},
            {"sign", options.flip ? -1 : 1},
            {"threshold_mode", options.thresholdMode ? "absolute" : "fraction_positive_max"},
            {"threshold_setting", options.threshold},
            {"covariance_cutoff", mask.cutoff},
            {"candidates", review.rows.size()},
            {"retained", mask.selected.size()},
            {"relative_residual", review.result.relativeResidual},
            {"converged", review.result.converged}};
}

SpectralSelection SpectralMaskDialog::selection() const {
    const auto review = reviews.constFind(row);
    if (review == reviews.cend() || component->currentIndex() < 0)
        return {};
    return selection(*review, maskOptions(), currentMask);
}

SpectralSelection SpectralMaskDialog::selection(const Review &review, const MaskOptions &options,
                                               const MaskPreview &mask) const {
    QVector<QPoint> points;
    points.reserve(mask.selected.size());
    for (int p : mask.selected)
        points.append(review.pixels[p]);
    auto rounds = review.history;
    rounds.append(roundRecord(review, options, mask));
    return {points,
            {{"reference_scan", QFileInfo(reference->currentText().trimmed()).absoluteFilePath()},
             {"energy_channel", settings.xrsEnergy},
             {"normalize_i0", settings.normalizeXrs},
             {"i0_channel", settings.xrsI0},
             {"repair_glitches", settings.repairGlitches},
             {"glitch_ratio", 1 - settings.glitchDropFraction},
             {"rounds", rounds}}};
}
