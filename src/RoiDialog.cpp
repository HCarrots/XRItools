#include "RoiDialog.h"
#include "ScientificPlot.h"
#include <QCryptographicHash>
#include <QtConcurrent/QtConcurrentRun>
#include <QtWidgets>
#include <algorithm>
#include <cmath>
#include <limits>

namespace {
QString detectorPath(int i) {
    return i == 0 ? "/entry/data/D_LAMBDA" : "/entry/data/D_MINIPIX";
}
QString stateText(RoiState s) {
    return s == RoiState::Valid     ? QStringLiteral("有效")
           : s == RoiState::Invalid ? QStringLiteral("失效")
                                    : QStringLiteral("未设置");
}
DataSlice lineProfile(const DataSlice &image, bool horizontal, int x, int y) {
    DataSlice result;
    result.path = QString("%1 · %2=%3")
                      .arg(image.path.section('/', -1), horizontal ? "Y" : "X")
                      .arg(horizontal ? y : x);
    result.xAxis = 0;
    result.yAxis = -1;
    result.width = horizontal ? image.width : image.height;
    result.height = 1;
    result.dims = {quint64(result.width)};
    result.starts = {0};
    result.strides = {1};
    result.minimum = std::numeric_limits<double>::infinity();
    result.maximum = -result.minimum;
    long double sum = 0;
    result.values.reserve(result.width);
    for (int i = 0; i < result.width; ++i) {
        double value = image.values[horizontal ? y * image.width + i : i * image.width + x];
        result.values.append(value);
        if (std::isfinite(value)) {
            result.minimum = std::min(result.minimum, value);
            result.maximum = std::max(result.maximum, value);
            sum += value;
            ++result.finiteCount;
        }
    }
    result.sum = double(sum);
    result.mean = result.finiteCount ? double(sum / result.finiteCount) : NAN;
    if (!result.finiteCount)
        result.minimum = result.maximum = NAN;
    return result;
}
} // namespace
RoiDialog::RoiDialog(const QString &file, const QString &workflow, QWidget *parent, bool autoDraft)
    : QDialog(parent), m_file(file), m_autoDraft(autoDraft) {
    setObjectName("CrystalRoiDialog");
    setWindowTitle(workflow + tr(" · 第一步：75 个分析晶体 ROI"));
    resize(1580, 940);
    setMinimumSize(1200, 760);
    auto main = new QVBoxLayout(this);
    auto heading = new QLabel(
        tr("选择晶体 → 在对应探测器上拖动框选；无法框选时标记失效。坐标采用 [x1, "
           "x2)、[y1, y2)，索引从 0 开始。列表支持 Ctrl 多选、Shift 连选、Ctrl+A 全选后批量标记失效。"));
    heading->setWordWrap(true);
    main->addWidget(heading);
    auto controls = new QHBoxLayout;
    m_mode = new QComboBox;
    m_mode->setObjectName("RoiImageMode");
    m_mode->addItems({tr("全部帧叠加（原始计数）"), tr("单帧（原始计数）")});
    controls->addWidget(m_mode);
    controls->addWidget(new QLabel(tr("帧")));
    m_frame = new QSpinBox;
    m_frame->setObjectName("RoiFrame");
    controls->addWidget(m_frame);
    auto log = new QCheckBox("signed log1p");
    log->setChecked(true);
    controls->addWidget(log);
    auto reset = new QPushButton(tr("复位视图"));
    controls->addWidget(reset);
    m_advance = new QCheckBox(tr("设置后跳到下一个未设置晶体"));
    m_advance->setChecked(true);
    controls->addWidget(m_advance);
    m_showProfiles = new QCheckBox(tr("显示 X/Y 剖面图"));
    m_showProfiles->setObjectName("ShowRoiProfiles");
    controls->addWidget(m_showProfiles);
    controls->addStretch();
    main->addLayout(controls);
    auto body = new QHBoxLayout;
    auto images = new QHBoxLayout;
    for (int i = 0; i < 2; ++i) {
        auto group = new QGroupBox(i == 0 ? tr("LAMBDA · VB / HL / VU / HB · 60 个晶体")
                                          : tr("MINIPIX · VD · 15 个晶体"));
        auto box = new QVBoxLayout(group);
        m_plots[i] = new ScientificPlot;
        m_plots[i]->setObjectName(i == 0 ? "RoiLambdaPlot" : "RoiMinipixPlot");
        m_plots[i]->setLogarithmic(true);
        m_plots[i]->setRoiSelectionEnabled(true);
        box->addWidget(m_plots[i]);
        images->addWidget(group, 1);
        connect(m_plots[i], &ScientificPlot::rectangleSelected, this,
                [this, i](const QRect &rect) { assign(i, rect); });
        connect(m_plots[i], &ScientificPlot::cursorChanged, this,
                [this](const QString &s) { m_message->setText(s); });
        connect(log, &QCheckBox::toggled, m_plots[i], &ScientificPlot::setLogarithmic);
        connect(reset, &QPushButton::clicked, m_plots[i], &ScientificPlot::resetView);
        try {
            auto node = HdfFile::inspect(file, detectorPath(i));
            if (node.dims.size() != 3)
                throw std::runtime_error("Detector is not a frame stack.");
            m_sizes.append(QSize(int(node.dims[2]), int(node.dims[1])));
            m_frames.append(node.dims[0]);
        } catch (const std::exception &) {
            m_sizes.append(QSize(0, 0));
            m_frames.append(0);
        }
    }
    auto imageArea = new QVBoxLayout;
    imageArea->addLayout(images, 1);
    m_profileGroup = new QGroupBox(tr("X / Y 像素剖面 · 基于当前叠加图或单帧的原始计数"));
    m_profileGroup->setObjectName("RoiProfileGroup");
    auto profileLayout = new QVBoxLayout(m_profileGroup);
    auto profileControls = new QHBoxLayout;
    m_profileDetector = new QComboBox;
    m_profileDetector->setObjectName("RoiProfileDetector");
    m_profileDetector->addItems({"LAMBDA", "MINIPIX"});
    profileControls->addWidget(new QLabel(tr("探测器")));
    profileControls->addWidget(m_profileDetector);
    m_profileX = new QSpinBox;
    m_profileX->setObjectName("RoiProfileX");
    m_profileY = new QSpinBox;
    m_profileY->setObjectName("RoiProfileY");
    profileControls->addWidget(new QLabel(tr("X（列）")));
    profileControls->addWidget(m_profileX);
    profileControls->addWidget(new QLabel(tr("Y（行）")));
    profileControls->addWidget(m_profileY);
    m_profileCenter = new QPushButton(tr("定位到当前 ROI 中心"));
    m_profileCenter->setObjectName("CenterRoiProfiles");
    profileControls->addWidget(m_profileCenter);
    auto profileLog = new QCheckBox(tr("剖面对数显示"));
    profileLog->setObjectName("RoiProfileLog");
    profileControls->addWidget(profileLog);
    profileControls->addStretch();
    profileLayout->addLayout(profileControls);
    auto profileHint = new QLabel(tr("拖动蓝色竖线调整 X、横线调整 Y，交点同时调整两者；Ctrl + "
                                     "拖动可在任意位置实时定位。固定 Y 查看行剖面，固定 X 查看列剖面。"));
    profileHint->setWordWrap(true);
    profileLayout->addWidget(profileHint);
    auto curves = new QHBoxLayout;
    for (int i = 0; i < 2; ++i) {
        m_profiles[i] = new ScientificPlot;
        m_profiles[i]->setObjectName(i == 0 ? "RoiHorizontalProfile" : "RoiVerticalProfile");
        m_profiles[i]->setMinimumSize(300, 225);
        m_profiles[i]->setAxisLabels(i == 0 ? "X / pixel" : "Y / pixel", "Intensity / counts");
        curves->addWidget(m_profiles[i], 1);
        connect(profileLog, &QCheckBox::toggled, m_profiles[i], &ScientificPlot::setLogarithmic);
    }
    profileLayout->addLayout(curves);
    imageArea->addWidget(m_profileGroup);
    m_profileGroup->hide();
    body->addLayout(imageArea, 1);
    auto panel = new QVBoxLayout;
    m_progress = new QLabel;
    panel->addWidget(m_progress);
    m_table = new QTableWidget(75, 4);
    m_table->setObjectName("CrystalRoiTable");
    m_table->setHorizontalHeaderLabels({tr("晶体"), tr("探测器"), tr("状态"), tr("ROI [x1,x2) × [y1,y2)")});
    m_table->setSelectionBehavior(QAbstractItemView::SelectRows);
    m_table->setSelectionMode(QAbstractItemView::ExtendedSelection);
    m_table->setToolTip(tr("Ctrl 多选 · Shift 连选 · Ctrl+A 全选；多选后点击“批量标记失效”。"));
    m_table->setEditTriggers(QAbstractItemView::NoEditTriggers);
    m_table->setAlternatingRowColors(true);
    m_table->verticalHeader()->hide();
    m_table->horizontalHeader()->setStretchLastSection(true);
    m_table->setColumnWidth(0, 65);
    m_table->setColumnWidth(1, 75);
    m_table->setColumnWidth(2, 55);
    m_table->setMinimumWidth(365);
    panel->addWidget(m_table, 1);
    m_current = new QLabel;
    m_current->setWordWrap(true);
    panel->addWidget(m_current);
    auto coords = new QGridLayout;
    for (int i = 0; i < 4; ++i) {
        coords->addWidget(new QLabel(QStringList{"x1", "x2", "y1", "y2"}[i]), i / 2, i % 2 * 2);
        m_coords[i] = new QSpinBox;
        m_coords[i]->setObjectName(QString("RoiCoordinate%1").arg(i));
        m_coords[i]->setRange(0, 4096);
        coords->addWidget(m_coords[i], i / 2, i % 2 * 2 + 1);
    }
    panel->addLayout(coords);
    auto apply = new QPushButton(tr("应用坐标 / 设为有效"));
    m_apply = apply;
    apply->setObjectName("ApplyRoiCoordinates");
    apply->setDefault(true);
    panel->addWidget(apply);
    auto rowButtons = new QHBoxLayout;
    auto invalid = new QPushButton(tr("标记失效"));
    m_invalid = invalid;
    invalid->setObjectName("InvalidateCrystalRoi");
    auto clear = new QPushButton(tr("清除重选"));
    m_clear = clear;
    clear->setObjectName("ClearCrystalRoi");
    rowButtons->addWidget(invalid);
    rowButtons->addWidget(clear);
    panel->addLayout(rowButtons);
    auto next = new QPushButton(tr("下一个未设置晶体"));
    panel->addWidget(next);
    body->addLayout(panel);
    main->addLayout(body, 1);
    m_message = new QLabel;
    m_message->setWordWrap(true);
    main->addWidget(m_message);
    auto bottom = new QHBoxLayout;
    auto load = new QPushButton(tr("加载 ROI 方案…"));
    auto save = new QPushButton(tr("保存 ROI 方案…"));
    m_finish = new QPushButton(tr("完成 ROI 选取"));
    m_finish->setObjectName("FinishCrystalRois");
    auto close = new QPushButton(tr("返回（保留草稿）"));
    close->setObjectName("ReturnFromRoi");
    bottom->addWidget(load);
    bottom->addWidget(save);
    bottom->addStretch();
    bottom->addWidget(m_finish);
    bottom->addWidget(close);
    main->addLayout(bottom);
    m_frame->setRange(0, int(std::max<quint64>(1, std::min(m_frames[0], m_frames[1])) - 1));
    m_frame->setEnabled(false);
    connect(m_mode, &QComboBox::currentIndexChanged, this, [this] { loadImages(); });
    connect(m_frame, &QSpinBox::valueChanged, this, [this] { loadImages(); });
    connect(m_table, &QTableWidget::currentCellChanged, this, [this] { selectRow(); });
    connect(m_table, &QTableWidget::itemSelectionChanged, this, &RoiDialog::selectRow);
    connect(m_showProfiles, &QCheckBox::toggled, this, [this](bool show) {
        m_profileGroup->setVisible(show);
        for (auto plot : m_plots)
            plot->setProfilePickingEnabled(show);
        updateProfiles();
    });
    connect(m_profileDetector, &QComboBox::currentIndexChanged, this, &RoiDialog::updateProfiles);
    connect(m_profileX, &QSpinBox::valueChanged, this, &RoiDialog::updateProfiles);
    connect(m_profileY, &QSpinBox::valueChanged, this, &RoiDialog::updateProfiles);
    for (int d = 0; d < 2; ++d)
        connect(m_plots[d], &ScientificPlot::profilePointPicked, this, [this, d](const QPoint &point) {
            QSignalBlocker a(m_profileDetector), b(m_profileX), c(m_profileY);
            m_profileDetector->setCurrentIndex(d);
            m_profileX->setRange(0, std::max(0, m_plots[d]->data().width - 1));
            m_profileY->setRange(0, std::max(0, m_plots[d]->data().height - 1));
            m_profileX->setValue(point.x());
            m_profileY->setValue(point.y());
            updateProfiles();
        });
    connect(m_profileCenter, &QPushButton::clicked, this, [this] {
        const auto rows = selectedRows();
        if (rows.size() != 1)
            return;
        const int row = rows.first();
        const auto &rect = m_rois.items()[row].rect;
        if (rect.isEmpty())
            return;
        int d = row < 60 ? 0 : 1;
        QSignalBlocker a(m_profileDetector), b(m_profileX), c(m_profileY);
        m_profileDetector->setCurrentIndex(d);
        m_profileX->setRange(0, std::max(0, m_plots[d]->data().width - 1));
        m_profileY->setRange(0, std::max(0, m_plots[d]->data().height - 1));
        m_profileX->setValue(rect.x() + (rect.width() - 1) / 2);
        m_profileY->setValue(rect.y() + (rect.height() - 1) / 2);
        updateProfiles();
    });
    connect(apply, &QPushButton::clicked, this, [this] {
        const auto rows = selectedRows();
        if (rows.size() != 1)
            return;
        int row = rows.first();
        assign(row < 60 ? 0 : 1,
               QRect(m_coords[0]->value(), m_coords[2]->value(), m_coords[1]->value() - m_coords[0]->value(),
                     m_coords[3]->value() - m_coords[2]->value()));
    });
    connect(invalid, &QPushButton::clicked, this, [this] {
        const auto rows = selectedRows();
        if (rows.isEmpty())
            return;
        for (int row : rows)
            m_rois.invalidate(row);
        refresh();
        storeDraft();
        if (rows.size() == 1 && m_advance->isChecked())
            nextPending();
    });
    connect(clear, &QPushButton::clicked, this, [this] {
        const auto rows = selectedRows();
        if (rows.size() != 1)
            return;
        m_rois.clear(rows.first());
        refresh();
        storeDraft();
    });
    connect(next, &QPushButton::clicked, this, &RoiDialog::nextPending);
    connect(close, &QPushButton::clicked, this, &RoiDialog::reject);
    connect(save, &QPushButton::clicked, this, [this] {
        auto path = QFileDialog::getSaveFileName(this, tr("保存 ROI 方案"),
                                                 QFileInfo(m_file).completeBaseName() + ".roi.json",
                                                 "ROI (*.roi.json)");
        if (!path.isEmpty())
            saveTo(path);
    });
    connect(load, &QPushButton::clicked, this, [this] {
        auto path = QFileDialog::getOpenFileName(this, tr("加载 ROI 方案"), {}, "ROI (*.roi.json)");
        if (!path.isEmpty())
            loadFrom(path);
    });
    connect(m_finish, &QPushButton::clicked, this, [this] {
        storeDraft();
        if (m_autoDraft && !saveTo(m_draftPath))
            return;
        accept();
    });
    if (m_autoDraft) {
        auto dir = QStandardPaths::writableLocation(QStandardPaths::AppLocalDataLocation) + "/roi-drafts";
        QDir().mkpath(dir);
        m_draftPath =
            dir + "/" +
            QCryptographicHash::hash(QFileInfo(file).canonicalFilePath().toUtf8(), QCryptographicHash::Sha256)
                .toHex() +
            ".json";
        if (QFileInfo::exists(m_draftPath))
            loadFrom(m_draftPath);
    }
    refresh();
    m_table->setCurrentCell(0, 0);
    connect(&m_watcher, &QFutureWatcher<QVector<DataSlice>>::finished, this, [this] {
        m_progressTimer.stop();
        m_loading = false;
        auto images = m_watcher.result();
        QStringList errors;
        for (int i = 0; i < 2; ++i) {
            m_plots[i]->setData(images[i]);
            m_plots[i]->setEnabled(images[i].error.isEmpty() && !images[i].values.isEmpty() &&
                                   !images[i].sampled);
            if (!images[i].error.isEmpty())
                errors << images[i].path + ": " + images[i].error;
            if (images[i].sampled)
                errors << tr("采样预览不能用于精确 ROI，请使用全部帧叠加模式。");
        }
        m_mode->setEnabled(true);
        m_frame->setEnabled(m_mode->currentIndex() == 1);
        m_message->setText(
            errors.isEmpty()
                ? tr("图像已加载。拖动框选当前晶体；滚轮缩放，双击复位。ROI 更改自动保存为本地草稿。 ")
                : errors.join('\n'));
        refresh();
        updateProfiles();
    });
    m_progressTimer.setInterval(250);
    connect(&m_progressTimer, &QTimer::timeout, this, [this] {
        m_message->setText(tr("正在读取 %1：%2% · 按帧读取，关闭窗口可取消。")
                               .arg(m_activeDetector.load() == 0 ? "LAMBDA" : "MINIPIX")
                               .arg(m_loadPercent.load()));
    });
    loadImages();
}
RoiDialog::~RoiDialog() {
    m_cancel = true;
    m_watcher.waitForFinished();
}
bool RoiDialog::imagesReady() const {
    return !m_loading && !m_plots[0]->data().values.isEmpty() && !m_plots[1]->data().values.isEmpty();
}
void RoiDialog::reject() {
    if (isWindow())
        m_cancel = true;
    QDialog::reject();
}
void RoiDialog::closeEvent(QCloseEvent *e) {
    m_cancel = true;
    QDialog::closeEvent(e);
}
void RoiDialog::loadImages() {
    if (m_loading)
        return;
    m_loading = true;
    m_cancel = false;
    m_loadPercent = 0;
    m_activeDetector = 0;
    m_progressTimer.start();
    m_mode->setEnabled(false);
    m_frame->setEnabled(false);
    for (auto plot : m_plots) {
        plot->setEnabled(false);
        plot->setData({});
    }
    updateProfiles();
    m_message->setText(tr("正在按帧读取探测器图像，请稍候；关闭窗口可取消。"));
    const auto file = m_file;
    const int mode = m_mode->currentIndex(), frame = m_frame->value();
    m_watcher.setFuture(QtConcurrent::run([this, file, mode, frame] {
        QVector<DataSlice> result;
        for (int i = 0; i < 2; ++i) {
            m_activeDetector = i;
            m_loadPercent = 0;
            if (mode == 0)
                result.append(HdfFile::sumFrames(file, detectorPath(i), m_cancel, &m_loadPercent));
            else
                result.append(HdfFile::read({file, detectorPath(i), {quint64(frame), 0, 0}, 2, 1}));
        }
        return result;
    }));
}
void RoiDialog::refresh() {
    QSignalBlocker block(m_table);
    for (int i = 0; i < 75; ++i) {
        const auto &r = m_rois.items()[i];
        QString rect = r.rect.isEmpty() ? QStringLiteral("—")
                                        : QString("[%1,%2) × [%3,%4)")
                                              .arg(r.rect.x())
                                              .arg(r.rect.x() + r.rect.width())
                                              .arg(r.rect.y())
                                              .arg(r.rect.y() + r.rect.height());
        QStringList fields{r.label, i < 60 ? "LAMBDA" : "MINIPIX", stateText(r.state), rect};
        for (int c = 0; c < 4; ++c) {
            auto item = m_table->item(i, c);
            if (!item) {
                item = new QTableWidgetItem;
                m_table->setItem(i, c, item);
            }
            item->setText(fields[c]);
            item->setToolTip(fields[c]);
            item->setForeground(r.state == RoiState::Invalid ? QColor("#777777")
                                : r.state == RoiState::Valid ? QColor("#176b35")
                                                             : palette().color(QPalette::Text));
        }
    }
    int pending = m_rois.count(RoiState::Pending);
    m_progress->setText(tr("共 75 个 · 有效 %1 · 失效 %2 · 未设置 %3")
                            .arg(m_rois.count(RoiState::Valid))
                            .arg(m_rois.count(RoiState::Invalid))
                            .arg(pending));
    m_finish->setEnabled(pending == 0 && !m_loading);
    selectRow();
}
QVector<int> RoiDialog::selectedRows() const {
    QVector<int> rows;
    for (const auto &index : m_table->selectionModel()->selectedRows())
        rows.append(index.row());
    std::sort(rows.begin(), rows.end());
    return rows;
}
void RoiDialog::selectRow() {
    const auto rows = selectedRows();
    const bool single = rows.size() == 1;
    m_profileCenter->setEnabled(single && !m_rois.items()[rows.first()].rect.isEmpty() && !m_loading);
    m_apply->setEnabled(single);
    m_clear->setEnabled(single);
    m_invalid->setEnabled(!rows.isEmpty());
    m_invalid->setText(rows.size() > 1 ? tr("批量标记失效（%1）").arg(rows.size()) : tr("标记失效"));
    for (auto coordinate : m_coords)
        coordinate->setEnabled(single);
    if (!single) {
        m_current->setText(
            rows.isEmpty() ? tr("请选择晶体")
                           : tr("已选 %1 个晶体 · 可批量标记失效；框选和坐标编辑需单选。").arg(rows.size()));
        for (int d = 0; d < 2; ++d) {
            QVector<CrystalRoi> overlays;
            for (const auto &roi : m_rois.items())
                if (roi.detector == (d == 0 ? "D_LAMBDA" : "D_MINIPIX"))
                    overlays.append(roi);
            m_plots[d]->setRoiOverlays(overlays, {});
        }
        return;
    }
    int row = rows.first();
    const auto &r = m_rois.items()[row];
    m_current->setText(tr("当前：%1 · %2").arg(r.label, stateText(r.state)));
    int v[] = {r.rect.x(), r.rect.x() + r.rect.width(), r.rect.y(), r.rect.y() + r.rect.height()};
    for (int i = 0; i < 4; ++i) {
        m_coords[i]->setMaximum(i < 2 ? m_sizes[row < 60 ? 0 : 1].width()
                                      : m_sizes[row < 60 ? 0 : 1].height());
        m_coords[i]->setValue(v[i]);
    }
    for (int d = 0; d < 2; ++d) {
        QVector<CrystalRoi> overlays;
        for (const auto &roi : m_rois.items())
            if (roi.detector == (d == 0 ? "D_LAMBDA" : "D_MINIPIX"))
                overlays.append(roi);
        m_plots[d]->setRoiOverlays(overlays, r.label);
    }
}
void RoiDialog::updateProfiles() {
    if (m_loading)
        m_profileCenter->setEnabled(false);
    for (auto plot : m_plots)
        plot->setProfilePosition(QPoint(-1, -1));
    const int detector = m_profileDetector->currentIndex();
    const auto &image = m_plots[detector]->data();
    const bool ready =
        !m_loading && image.error.isEmpty() && !image.values.isEmpty() && image.yAxis >= 0 && !image.sampled;
    m_profileX->setEnabled(ready);
    m_profileY->setEnabled(ready);
    if (!m_showProfiles->isChecked() || !ready) {
        for (auto plot : m_profiles)
            plot->setData({});
        return;
    }
    {
        QSignalBlocker x(m_profileX), y(m_profileY);
        m_profileX->setRange(0, image.width - 1);
        m_profileY->setRange(0, image.height - 1);
    }
    const int x = m_profileX->value(), y = m_profileY->value();
    m_plots[detector]->setProfilePosition(QPoint(x, y));
    m_profiles[0]->setData(lineProfile(image, true, x, y));
    m_profiles[1]->setData(lineProfile(image, false, x, y));
}
void RoiDialog::assign(int detector, const QRect &rect) {
    const auto rows = selectedRows();
    if (rows.size() != 1) {
        m_message->setText(tr("框选 ROI 时请只选择一个晶体。"));
        return;
    }
    int row = rows.first();
    if (row < 0 || m_loading)
        return;
    if ((row < 60 ? 0 : 1) != detector) {
        m_message->setText(tr("当前晶体属于另一台探测器，请在对应图像上框选。"));
        return;
    }
    if (!m_plots[detector]->isEnabled() || !m_rois.assign(row, rect, m_sizes[detector])) {
        m_message->setText(tr("ROI 必须位于探测器内，且 x2 > x1、y2 > y1。"));
        return;
    }
    refresh();
    storeDraft();
    if (m_advance->isChecked())
        nextPending();
}
void RoiDialog::nextPending() {
    int row = m_table->currentRow();
    for (int step = 1; step <= 75; ++step) {
        int next = (std::max(0, row) + step) % 75;
        if (m_rois.items()[next].state == RoiState::Pending) {
            m_table->setCurrentCell(next, 0, QItemSelectionModel::ClearAndSelect | QItemSelectionModel::Rows);
            m_table->scrollToItem(m_table->item(next, 0));
            return;
        }
    }
    m_message->setText(tr("75 个晶体均已处理，可完成选取。"));
}
bool RoiDialog::saveTo(const QString &path) {
    QSaveFile file(path);
    if (!file.open(QIODevice::WriteOnly)) {
        m_message->setText(file.errorString());
        return false;
    }
    auto bytes = m_rois.serialize(m_file, m_sizes);
    if (file.write(bytes) != bytes.size() || !file.commit()) {
        m_message->setText(file.errorString());
        return false;
    }
    m_message->setText(tr("ROI 方案已保存：") + path);
    return true;
}
bool RoiDialog::loadFrom(const QString &path) {
    QFile file(path);
    QString error;
    if (!file.open(QIODevice::ReadOnly)) {
        m_message->setText(file.errorString());
        return false;
    }
    if (file.size() > 1024 * 1024 || !m_rois.restore(file.readAll(), m_sizes, &error)) {
        m_message->setText(error.isEmpty() ? tr("ROI 方案过大。") : error);
        return false;
    }
    refresh();
    storeDraft();
    return true;
}
void RoiDialog::storeDraft() {
    if (m_autoDraft && !m_draftPath.isEmpty())
        saveTo(m_draftPath);
}
