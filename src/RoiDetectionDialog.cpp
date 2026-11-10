#include "RoiDetectionDialog.h"
#include "CrystalRoi.h"
#include "HdfFile.h"
#include "WorkflowUi.h"
#include <QApplication>
#include <QGraphicsRectItem>
#include <QGraphicsScene>
#include <QGraphicsView>
#include <QJsonArray>
#include <QJsonDocument>
#include <QJsonObject>
#include <QMenu>
#include <QMessageBox>
#include <QMouseEvent>
#include <QSaveFile>
#include <QSet>
#include <QSettings>
#include <QStandardPaths>
#include <climits>
#include <cmath>
#include <functional>

namespace {
class MaskView final : public QGraphicsView {
  public:
    explicit MaskView(QGraphicsScene *scene) : QGraphicsView(scene) {
        setDragMode(ScrollHandDrag);
        setTransformationAnchor(AnchorUnderMouse);
        setHorizontalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
        setVerticalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
    }
    std::function<void(const QPoint &)> pick;

    void resetZoom() {
        fit = true;
        fitInView(sceneRect(), Qt::KeepAspectRatio);
    }

  protected:
    void resizeEvent(QResizeEvent *event) override {
        QGraphicsView::resizeEvent(event);
        if (fit)
            resetZoom();
    }
    void wheelEvent(QWheelEvent *event) override {
        if (!event->angleDelta().y())
            return;
        fit = false;
        const double factor = event->angleDelta().y() > 0 ? 1.25 : .8;
        const double zoom = transform().m11() * factor;
        if (zoom >= .05 && zoom <= 64.)
            scale(factor, factor);
        event->accept();
    }
    void mousePressEvent(QMouseEvent *event) override {
        press = event->position();
        QGraphicsView::mousePressEvent(event);
    }
    void mouseReleaseEvent(QMouseEvent *event) override {
        const auto position = mapToScene(event->position().toPoint());
        QGraphicsView::mouseReleaseEvent(event);
        if (event->button() == Qt::LeftButton &&
            (event->position() - press).manhattanLength() < QApplication::startDragDistance() &&
            sceneRect().contains(position) && pick)
            pick({int(std::floor(position.x())), int(std::floor(position.y()))});
    }

  private:
    QPointF press;
    bool fit = true;
};
} // namespace

RoiDetectionDialog::RoiDetectionDialog(const QString &source, QWidget *parent) : QDialog(parent) {
    setWindowTitle(tr("自动识别光斑 ROI"));
    resize(1150, 680);
    auto layout = new QVBoxLayout(this);
    m_parameters = new QWidget;
    auto form = new QFormLayout(m_parameters);
    WorkflowUi::compact(form);
    auto fileRow = [this, form](const QString &label, const QString &value, const QString &filter) {
        auto row = new QHBoxLayout;
        auto edit = new QLineEdit(value);
        auto browse = new QPushButton(tr("浏览…"));
        row->addWidget(edit, 1);
        row->addWidget(browse);
        form->addRow(label, row);
        connect(browse, &QPushButton::clicked, this, [this, edit, label, filter] {
            const auto file = QFileDialog::getOpenFileName(this, label, edit->text(), filter);
            if (!file.isEmpty()) {
                edit->setText(file);
                if (edit == m_source)
                    inspectSource();
            }
        });
        return edit;
    };
    QSettings settings;
    const auto defaultPython = QStandardPaths::findExecutable(
#ifdef Q_OS_WIN
        "python"
#else
        "python3"
#endif
    );
    m_source = fileRow(tr("扫描文件"), source, tr("NeXus / HDF5 (*.nxs *.h5 *.hdf5);;所有文件 (*)"));
    m_python = fileRow(tr("Python 环境"), settings.value("roiDetection/python", defaultPython).toString(),
                       tr("所有文件 (*)"));
    m_python->setObjectName("SpotRoiPython");
    const QDir appDirectory(QCoreApplication::applicationDirPath());
    auto defaultWeights = appDirectory.filePath("models/roi-maskrcnn.pt");
#ifdef Q_OS_LINUX
    if (!QFileInfo::exists(defaultWeights))
        defaultWeights = appDirectory.absoluteFilePath("../share/xri-tools/models/roi-maskrcnn.pt");
#endif
    m_weights = fileRow(tr("模型权重"), settings.value("roiDetection/weights", defaultWeights).toString(),
                        tr("XRI Mask R-CNN (*.pt);;所有文件 (*)"));
    m_weights->setObjectName("SpotRoiWeights");
    auto choices = new QHBoxLayout;
    m_detector = new QComboBox;
    m_sum = new QCheckBox(tr("全部帧求和"));
    m_sum->setChecked(true);
    m_frame = new QSpinBox;
    m_frame->setEnabled(false);
    m_score = new QDoubleSpinBox;
    m_maskThreshold = new QDoubleSpinBox;
    for (auto *box : {m_score, m_maskThreshold}) {
        box->setRange(.01, 1.);
        box->setSingleStep(.05);
    }
    m_score->setValue(.55);
    m_maskThreshold->setValue(.50);
    choices->addWidget(m_detector);
    choices->addWidget(m_sum);
    choices->addWidget(new QLabel(tr("帧索引")));
    choices->addWidget(m_frame);
    choices->addSpacing(12);
    choices->addWidget(new QLabel(tr("置信度")));
    choices->addWidget(m_score);
    choices->addWidget(new QLabel(tr("掩膜阈值")));
    choices->addWidget(m_maskThreshold);
    choices->addStretch();
    form->addRow(tr("探测器"), choices);
    layout->addWidget(m_parameters);
    auto display = new QSplitter(Qt::Horizontal);
    m_preview = new QLabel(tr("识别后显示编号和掩膜叠加图"));
    m_preview->setAlignment(Qt::AlignCenter);
    m_preview->setSizePolicy(QSizePolicy::Ignored, QSizePolicy::Ignored);
    m_preview->setMinimumWidth(320);
    m_preview->installEventFilter(this);
    display->addWidget(m_preview);
    m_table = new QTableWidget(0, 10);
    m_table->setObjectName("SpotRoiResults");
    m_table->setHorizontalHeaderLabels(
        {tr("编号"), tr("置信度"), tr("X 质心"), tr("Y 质心"), tr("掩膜像素"), tr("原始积分"), tr("净积分"),
         tr("长轴 / px"), tr("短轴 / px"), tr("角度 / °")});
    m_table->setEditTriggers(QAbstractItemView::NoEditTriggers);
    m_table->setSelectionBehavior(QAbstractItemView::SelectRows);
    m_table->horizontalHeader()->setSectionResizeMode(QHeaderView::ResizeToContents);
    m_table->verticalHeader()->hide();
    display->addWidget(m_table);
    display->setStretchFactor(0, 4);
    display->setStretchFactor(1, 5);
    layout->addWidget(display, 1);
    auto actions = new QHBoxLayout;
    m_message = new QLabel;
    m_message->setObjectName("SpotRoiStatus");
    m_message->setWordWrap(true);
    m_run = WorkflowUi::button("play", tr("识别光斑"));
    m_run->setObjectName("RunSpotRoiDetection");
    m_cancel = WorkflowUi::button("square", tr("取消"));
    m_export = WorkflowUi::button("download", tr("导出掩膜与参数…"));
    m_export->setObjectName("ExportSpotRoiDetection");
    m_match = WorkflowUi::button("save", tr("手动匹配晶体并保存 ROI…"));
    auto matchMenu = new QMenu(m_match);
    m_matchCurrent = matchMenu->addAction(tr("使用本次识别结果…"), this, [this] { matchAndSave(m_result); });
    auto matchSaved = matchMenu->addAction(tr("载入已保存的识别结果…"), this, [this] { loadResultForMatching(); });
    matchSaved->setObjectName("MatchSavedSpotRois");
    m_match->setMenu(matchMenu);
    actions->addWidget(m_message, 1);
    actions->addWidget(m_run);
    actions->addWidget(m_cancel);
    actions->addWidget(m_export);
    actions->addWidget(m_match);
    layout->addLayout(actions);
    setRunning(false);
    connect(m_source, &QLineEdit::editingFinished, this, [this] { inspectSource(); });
    connect(m_detector, &QComboBox::currentIndexChanged, this, [this] { updateFrameRange(); });
    connect(m_sum, &QCheckBox::toggled, this, [this] { updateFrameRange(); });
    connect(m_run, &QPushButton::clicked, this, [this] { run(); });
    connect(m_cancel, &QPushButton::clicked, this, [this] {
        m_cancelled = true;
        m_process.kill();
    });
    connect(m_export, &QPushButton::clicked, this, [this] { exportResult(); });
    connect(&m_process, &QProcess::finished, this, [this](int code, QProcess::ExitStatus status) {
        completed(code, status);
    });
    connect(&m_process, &QProcess::errorOccurred, this, [this](QProcess::ProcessError error) {
        if (error == QProcess::FailedToStart) {
            setRunning(false);
            m_message->setText(tr("无法启动 Python。请核对环境路径。%1").arg(m_process.errorString()));
        }
    });
    connect(&m_process, &QProcess::readyReadStandardOutput, this, [this] {
        const auto output = QString::fromUtf8(m_process.readAllStandardOutput()).trimmed();
        if (output.endsWith("Reading detector data"))
            m_message->setText(tr("正在读取探测器数据…"));
        else if (output.endsWith("Running Mask R-CNN"))
            m_message->setText(tr("正在识别光斑掩膜…"));
    });
    inspectSource();
}

RoiDetectionDialog::~RoiDetectionDialog() {
    disconnect(&m_process, nullptr, this, nullptr);
    if (m_process.state() != QProcess::NotRunning) {
        m_process.kill();
        m_process.waitForFinished(1000);
    }
}

void RoiDetectionDialog::inspectSource() {
    const auto last = m_detector->currentText();
    const QSignalBlocker blocked(m_detector);
    m_detector->clear();
    const auto path = m_source->text().trimmed();
    if (QFileInfo(path).isFile()) {
        for (const auto &detector : {QStringLiteral("D_LAMBDA"), QStringLiteral("D_MINIPIX")}) {
            const auto shape = HdfFile::inspect(path, "/entry/data/" + detector).dims;
            if (shape.size() == 3 && shape[0] > 0 && shape[0] <= INT_MAX)
                m_detector->addItem(detector, int(shape[0]));
        }
    }
    m_detector->setCurrentIndex(std::max(0, m_detector->findText(last)));
    updateFrameRange();
    m_run->setEnabled(m_detector->count() > 0);
    m_message->setText(m_detector->count() ? tr("选择探测器和帧范围后运行识别。")
                                         : tr("请选择含 LAMBDA 或 MINIPIX 探测器数据的扫描文件。"));
}

bool RoiDetectionDialog::eventFilter(QObject *watched, QEvent *event) {
    if (watched == m_preview && event->type() == QEvent::Resize)
        resizePreview();
    return QDialog::eventFilter(watched, event);
}

void RoiDetectionDialog::resizePreview() {
    if (!m_overlay.isNull())
        m_preview->setPixmap(m_overlay.scaled(m_preview->size(), Qt::KeepAspectRatio, Qt::SmoothTransformation));
}

void RoiDetectionDialog::updateFrameRange() {
    m_frame->setMaximum(std::max(0, m_detector->currentData().toInt() - 1));
    m_frame->setEnabled(!m_sum->isChecked() && m_detector->count() > 0);
}

void RoiDetectionDialog::setRunning(bool running) {
    m_parameters->setEnabled(!running);
    m_run->setEnabled(!running && m_detector->count() > 0);
    m_cancel->setEnabled(running);
    m_export->setEnabled(!running && !m_result.isEmpty());
    m_match->setEnabled(!running);
    m_matchCurrent->setEnabled(!m_result.isEmpty());
}

void RoiDetectionDialog::run() {
    if (m_process.state() != QProcess::NotRunning || !m_detector->count())
        return;
    const auto python = m_python->text().trimmed(), weights = m_weights->text().trimmed();
    if (python.isEmpty() || !QFileInfo(weights).isFile()) {
        m_message->setText(tr("请指定含 PyTorch 和 torchvision 的 Python 环境及已训练的 XRI 模型。"));
        return;
    }
    if (!m_work.isValid()) {
        m_message->setText(tr("无法建立识别工作目录。"));
        return;
    }
    const auto script = m_work.filePath("roi_maskrcnn.py");
    if (!QFile::exists(script) && !QFile::copy(":/scripts/roi_maskrcnn.py", script)) {
        m_message->setText(tr("无法载入 Mask R-CNN 推理脚本。"));
        return;
    }
    QSettings settings;
    settings.setValue("roiDetection/python", python);
    settings.setValue("roiDetection/weights", weights);
    m_result.clear();
    m_overlay = {};
    m_table->setRowCount(0);
    m_preview->setText(tr("正在生成掩膜…"));
    m_cancelled = false;
    m_message->setText(tr("正在启动 Mask R-CNN…"));
    setRunning(true);
    auto environment = QProcessEnvironment::systemEnvironment();
    environment.insert("PYTHONIOENCODING", "utf-8");
    m_process.setProcessEnvironment(environment);
    m_process.start(python, {"-u", script, "infer", "--file", QFileInfo(m_source->text().trimmed()).absoluteFilePath(),
                             "--detector", m_detector->currentText(), "--frame",
                             QString::number(m_sum->isChecked() ? -1 : m_frame->value()), "--weights",
                             QFileInfo(weights).absoluteFilePath(), "--score", QString::number(m_score->value()),
                             "--mask-threshold", QString::number(m_maskThreshold->value()), "--output",
                             m_work.filePath("spots.json"), "--preview", m_work.filePath("spots.png")});
}

void RoiDetectionDialog::completed(int code, QProcess::ExitStatus status) {
    setRunning(false);
    if (m_cancelled) {
        m_preview->setText(tr("识别已取消"));
        m_message->setText(tr("识别已取消。"));
        return;
    }
    if (code != 0 || status != QProcess::NormalExit) {
        const auto detail = QString::fromUtf8(m_process.readAllStandardError()).trimmed().section('\n', -1).left(350);
        m_preview->setText(tr("未生成掩膜"));
        m_message->setText(tr("识别失败。%1").arg(detail));
        return;
    }
    QFile file(m_work.filePath("spots.json"));
    QPixmap preview(m_work.filePath("spots.png"));
    const auto bytes = file.open(QIODevice::ReadOnly) ? file.readAll() : QByteArray{};
    const auto object = QJsonDocument::fromJson(bytes).object();
    if (object.value("schema").toString() != "xri.roi.maskrcnn.result.v1" || preview.isNull()) {
        m_message->setText(tr("无法读取识别结果。"));
        return;
    }
    m_result = bytes;
    m_overlay = std::move(preview);
    resizePreview();
    const auto spots = object.value("spots").toArray();
    m_table->setRowCount(int(spots.size()));
    for (int row = 0; row < spots.size(); ++row) {
        const auto spot = spots.at(row).toObject();
        const auto center = spot.value("centroid_xy").toArray();
        const double values[] = {spot.value("id").toDouble(), spot.value("score").toDouble(),
                                 center[0].toDouble(), center[1].toDouble(), spot.value("mask_area").toDouble(),
                                 spot.value("integrated_counts").toDouble(), spot.value("net_counts").toDouble(),
                                 spot.value("ellipse_major_px").toDouble(), spot.value("ellipse_minor_px").toDouble(),
                                 spot.value("ellipse_angle_deg").toDouble()};
        for (int col = 0; col < 10; ++col) {
            auto item = new QTableWidgetItem(QString::number(values[col], 'g', 8));
            item->setTextAlignment(Qt::AlignRight | Qt::AlignVCenter);
            m_table->setItem(row, col, item);
        }
    }
    const int frame = object.value("frame").toInt(-1);
    m_message->setText(tr("%1 · %2 · 识别到 %3 个光斑")
                           .arg(object.value("detector").toString(),
                                frame < 0 ? tr("全部帧求和") : tr("第 %1 帧").arg(frame))
                           .arg(spots.size()));
    m_export->setEnabled(true);
    m_matchCurrent->setEnabled(!spots.isEmpty());
}

void RoiDetectionDialog::exportResult() {
    if (m_result.isEmpty())
        return;
    const auto path = QFileDialog::getSaveFileName(this, tr("导出掩膜与参数"), "spot-rois.json", tr("JSON (*.json)"));
    if (path.isEmpty())
        return;
    QSaveFile file(path);
    if (file.open(QIODevice::WriteOnly) && file.write(m_result) == m_result.size() && file.commit())
        m_message->setText(tr("已保存掩膜与参数到 %1").arg(path));
    else
        m_message->setText(tr("保存失败。%1").arg(file.errorString()));
}

void RoiDetectionDialog::loadResultForMatching() {
    const auto path = QFileDialog::getOpenFileName(this, tr("载入自动识别结果"), {}, tr("JSON (*.json)"));
    if (path.isEmpty())
        return;
    QFile file(path);
    if (!file.open(QIODevice::ReadOnly)) {
        m_message->setText(tr("无法读取识别结果。%1").arg(file.errorString()));
        return;
    }
    matchAndSave(file.readAll());
}

void RoiDetectionDialog::matchAndSave(const QByteArray &result) {
    const auto report = QJsonDocument::fromJson(result).object();
    const auto detector = report["detector"].toString();
    const auto shape = report["shape_yx"].toArray();
    const auto spots = report["spots"].toArray();
    const int d = detector == "D_LAMBDA" ? 0 : detector == "D_MINIPIX" ? 1 : -1;
    if (report["schema"] != "xri.roi.maskrcnn.result.v1" || d < 0 || shape.size() != 2 || spots.isEmpty()) {
        m_message->setText(tr("请选择含有效光斑的自动识别结果。"));
        return;
    }
    const auto source = report["source"].toString();
    const QSize size(shape[1].toInt(), shape[0].toInt());
    if (size.isEmpty() || size.width() > 4096 || size.height() > 4096) {
        m_message->setText(tr("识别结果的探测器尺寸无效。"));
        return;
    }
    QVector<QSize> sourceSizes(2, QSize(0, 0));
    if (QFileInfo(source).isFile()) {
        for (int i = 0; i < 2; ++i) {
            const auto path = i == 0 ? "/entry/data/D_LAMBDA" : "/entry/data/D_MINIPIX";
            const auto dims = HdfFile::inspect(source, path).dims;
            if (dims.size() == 3 && dims[1] <= 4096 && dims[2] <= 4096)
                sourceSizes[i] = QSize(int(dims[2]), int(dims[1]));
        }
    }
    if (!sourceSizes[d].isEmpty() && sourceSizes[d] != size) {
        m_message->setText(tr("识别结果的尺寸与原始扫描不一致。"));
        return;
    }
    sourceSizes[d] = size;
    auto sizes = sourceSizes;

    QVector<QVector<QPoint>> masks;
    QVector<QPointF> centers;
    QVector<int> ids;
    QMultiHash<QPoint, int> hits;
    for (const auto &value : spots) {
        const auto spot = value.toObject();
        auto mask = CrystalRois::decodeMask(spot["mask_runs"].toArray(), size);
        const auto center = spot["centroid_xy"].toArray();
        const int id = spot["id"].toInt();
        if (mask.isEmpty() || id <= 0 || ids.contains(id) || center.size() != 2 ||
            !center[0].isDouble() || !center[1].isDouble()) {
            m_message->setText(tr("识别结果中的光斑编号、质心或掩膜无效。"));
            return;
        }
        for (const auto &p : mask)
            hits.insert(p, int(masks.size()));
        ids.append(id);
        centers.append({center[0].toDouble(), center[1].toDouble()});
        masks.append(std::move(mask));
    }

    CrystalRois output;
    const int first = d == 0 ? 0 : 60, count = d == 0 ? 60 : 15;
    QVector<int> matches(count, -1);
    QDialog review(this);
    review.setWindowTitle(tr("手动匹配晶体掩膜 · %1").arg(detector));
    review.resize(1050, 660);
    auto layout = new QVBoxLayout(&review);
    auto info = new QLabel(tr("先在左侧选择晶体，再点击右侧对应的光斑掩膜。滚轮缩放，拖动平移。"));
    info->setWordWrap(true);
    layout->addWidget(info);

    auto split = new QSplitter(Qt::Horizontal);
    auto table = new QTableWidget(count, 3);
    table->setObjectName("SpotRoiCrystalMatches");
    table->setHorizontalHeaderLabels({tr("晶体"), tr("状态"), tr("光斑编号")});
    table->setEditTriggers(QAbstractItemView::NoEditTriggers);
    table->setSelectionBehavior(QAbstractItemView::SelectRows);
    table->setSelectionMode(QAbstractItemView::ExtendedSelection);
    table->verticalHeader()->hide();
    table->horizontalHeader()->setSectionResizeMode(QHeaderView::ResizeToContents);
    table->horizontalHeader()->setStretchLastSection(true);
    for (int row = 0; row < count; ++row) {
        table->setItem(row, 0, new QTableWidgetItem(output.items()[first + row].label));
        table->setItem(row, 1, new QTableWidgetItem);
        table->setItem(row, 2, new QTableWidgetItem);
    }
    split->addWidget(table);
    auto scene = new QGraphicsScene(&review);
    scene->setSceneRect(0, 0, size.width(), size.height());
    const bool hasPreview = result == m_result && m_overlay.size() == size;
    QImage base(size, QImage::Format_RGB32);
    if (hasPreview)
        base = m_overlay.toImage().convertToFormat(QImage::Format_RGB32);
    else
        base.fill(QColor(31, 35, 43));
    auto image = scene->addPixmap({});
    QVector<QGraphicsSimpleTextItem *> labels;
    for (const int id : ids) {
        auto label = scene->addSimpleText(QString::number(id));
        label->setFlag(QGraphicsItem::ItemIgnoresTransformations);
        labels.append(label);
    }
    QPen outlinePen(QColor(250, 203, 80), 2);
    outlinePen.setCosmetic(true);
    auto outline = scene->addRect({}, outlinePen);
    auto view = new MaskView(scene);
    view->setObjectName("SpotRoiMaskView");
    split->addWidget(view);
    split->setSizes({280, 730});
    split->setStretchFactor(1, 1);
    layout->addWidget(split, 1);

    auto controls = new QHBoxLayout;
    auto load = WorkflowUi::button("folder-open", tr("载入已有方案…"));
    load->setObjectName("LoadSpotRoiBaseline");
    auto invalid = new QPushButton(tr("标记失效"));
    auto clear = new QPushButton(tr("清除所选 ROI"));
    auto reset = new QPushButton(tr("复位视图"));
    auto advance = new QCheckBox(tr("设置后跳到下一个未设置晶体"));
    advance->setChecked(true);
    controls->addWidget(load);
    controls->addWidget(invalid);
    controls->addWidget(clear);
    controls->addWidget(advance);
    controls->addStretch();
    controls->addWidget(reset);
    layout->addLayout(controls);
    auto actions = new QHBoxLayout;
    auto status = new QLabel(tr("可从空白方案开始，也可载入已有方案继续编辑。"));
    status->setWordWrap(true);
    auto save = WorkflowUi::button("save", tr("保存共用 ROI 方案…"));
    save->setObjectName("SaveMatchedSpotRois");
    auto cancel = new QPushButton(tr("取消"));
    actions->addWidget(status, 1);
    actions->addWidget(save);
    actions->addWidget(cancel);
    layout->addLayout(actions);

    auto selectedRows = [&] { return table->selectionModel()->selectedRows(); };
    auto selectRow = [&](int row) {
        table->setCurrentCell(row, 0, QItemSelectionModel::ClearAndSelect | QItemSelectionModel::Rows);
        table->scrollToItem(table->item(row, 0));
    };
    auto showSelection = [&] {
        const auto selected = selectedRows();
        outline->hide();
        if (selected.size() == 1) {
            const auto &roi = output.items()[first + selected.first().row()];
            if (roi.state == RoiState::Valid) {
                outline->setRect(QRectF(roi.rect));
                outline->show();
            }
        }
    };
    auto refresh = [&] {
        auto picture = base;
        for (int i = 0; i < masks.size(); ++i) {
            const int row = matches.indexOf(i);
            if (!hasPreview || row >= 0) {
                const QColor color = row >= 0 ? QColor(42, 209, 138) : QColor(106, 151, 247);
                for (const auto &pixel : masks[i])
                    picture.setPixelColor(pixel, color);
            }
            labels[i]->setPos(centers[i] + QPointF(.5, .5));
            labels[i]->setBrush(row >= 0 ? QColor(255, 238, 132) : Qt::white);
            labels[i]->setText(row >= 0 ? output.items()[first + row].label : QString::number(ids[i]));
            labels[i]->setVisible(!hasPreview || row >= 0);
        }
        image->setPixmap(QPixmap::fromImage(picture));
        for (int row = 0; row < count; ++row) {
            const auto &roi = output.items()[first + row];
            table->item(row, 1)->setText(roi.state == RoiState::Pending ? tr("未设置") :
                                         roi.state == RoiState::Invalid ? tr("已失效") :
                                         roi.mask.isEmpty() ? tr("矩形 ROI") : tr("掩膜"));
            table->item(row, 2)->setText(matches[row] < 0 ? QString() : QString::number(ids[matches[row]]));
        }
        save->setEnabled(matches.count(-1) < count);
        showSelection();
    };
    auto advanceFrom = [&](int row) {
        if (!advance->isChecked())
            return;
        for (int step = 1; step < count; ++step) {
            const int next = (row + step) % count;
            if (output.items()[first + next].state == RoiState::Pending) {
                selectRow(next);
                break;
            }
        }
    };
    connect(table, &QTableWidget::itemSelectionChanged, &review, showSelection);
    connect(reset, &QPushButton::clicked, &review, [view] { view->resetZoom(); });
    connect(cancel, &QPushButton::clicked, &review, &QDialog::reject);
    connect(invalid, &QPushButton::clicked, &review, [&] {
        const auto rows = selectedRows();
        for (const auto &row : rows) {
            output.invalidate(first + row.row());
            matches[row.row()] = -1;
        }
        refresh();
        if (rows.size() == 1)
            advanceFrom(rows.first().row());
    });
    connect(clear, &QPushButton::clicked, &review, [&] {
        for (const auto &row : selectedRows()) {
            output.clear(first + row.row());
            matches[row.row()] = -1;
        }
        refresh();
    });
    view->pick = [&](const QPoint &pixel) {
        const auto selected = selectedRows();
        if (selected.size() != 1) {
            status->setText(tr("匹配时请选择一个晶体。"));
            return;
        }
        const auto clicked = hits.values(pixel);
        if (clicked.isEmpty()) {
            status->setText(tr("此处没有光斑掩膜。请点击着色区域。"));
            return;
        }
        int spot = clicked.first();
        if (clicked.size() > 1) {
            QMenu choose(&review);
            for (int i : clicked)
                choose.addAction(tr("光斑 %1").arg(ids[i]))->setData(i);
            auto action = choose.exec(QCursor::pos());
            if (!action)
                return;
            spot = action->data().toInt();
        }
        const int row = selected.first().row(), owner = int(matches.indexOf(spot));
        if (owner >= 0 && owner != row) {
            status->setText(tr("光斑 %1 已匹配 %2。请先清除该晶体的匹配。")
                                .arg(ids[spot]).arg(output.items()[first + owner].label));
            return;
        }
        if (!output.assignMask(first + row, masks[spot], size))
            return;
        matches[row] = spot;
        status->setText(tr("%1 已匹配光斑 %2。").arg(output.items()[first + row].label).arg(ids[spot]));
        refresh();
        advanceFrom(row);
    };
    connect(load, &QPushButton::clicked, &review, [&] {
        if (matches.count(-1) < count &&
            QMessageBox::question(&review, tr("载入已有方案"), tr("载入会替换本窗口尚未保存的匹配。继续吗？")) != QMessageBox::Yes)
            return;
        const auto path = QFileDialog::getOpenFileName(&review, tr("载入已有 ROI 方案"), QFileInfo(source).absolutePath(),
                                                       tr("XRI ROI (*.roi.json *.json)"));
        if (path.isEmpty())
            return;
        QFile file(path);
        if (!file.open(QIODevice::ReadOnly)) {
            status->setText(tr("无法读取方案。%1").arg(file.errorString()));
            return;
        }
        const auto bytes = file.readAll();
        const auto geometry = QJsonDocument::fromJson(bytes).object()["detector_sizes"].toArray();
        auto imported = sourceSizes;
        for (int i = 0; i < 2; ++i) {
            if (imported[i].isEmpty() && geometry.size() == 2) {
                const auto pair = geometry[i].toArray();
                if (pair.size() == 2 && pair[0].toInt() >= 0 && pair[1].toInt() >= 0 &&
                    pair[0].toInt() <= 4096 && pair[1].toInt() <= 4096)
                    imported[i] = QSize(pair[0].toInt(), pair[1].toInt());
            }
        }
        QString error;
        if (!output.restore(bytes, imported, &error)) {
            status->setText(error);
            return;
        }
        sizes = imported;
        matches.fill(-1);
        refresh();
        status->setText(tr("已载入原方案。请选择晶体并点击光斑掩膜。"));
    });
    connect(save, &QPushButton::clicked, &review, [&] {
        QHash<QPoint, QString> owners;
        for (int row = 0; row < count; ++row) {
            const auto &roi = output.items()[first + row];
            if (roi.state != RoiState::Valid)
                continue;
            for (const auto &p : roi.pixels()) {
                const auto other = owners.constFind(p);
                if (other != owners.cend()) {
                    selectRow(row);
                    status->setText(tr("晶体 %1 和 %2 共享探测器像素。请调整匹配或标记其中一个失效。")
                                        .arg(roi.label, *other));
                    return;
                }
                owners.insert(p, roi.label);
            }
        }
        QFileDialog picker(&review, tr("保存共用 ROI 方案"), QFileInfo(source).absolutePath(),
                           tr("XRI ROI (*.roi.json)"));
        picker.setAcceptMode(QFileDialog::AcceptSave);
        picker.setDefaultSuffix("roi.json");
        picker.selectFile("matched-spots.roi.json");
        if (picker.exec() != QDialog::Accepted)
            return;
        const auto path = picker.selectedFiles().first();
        if (!path.endsWith(".roi.json", Qt::CaseInsensitive)) {
            status->setText(tr("共用 ROI 方案的文件名应以 .roi.json 结尾。"));
            return;
        }
        const auto scheme = output.serialize(source, sizes);
        QSaveFile destination(path);
        if (!destination.open(QIODevice::WriteOnly) || destination.write(scheme) != scheme.size() ||
            !destination.commit()) {
            status->setText(tr("无法保存方案。%1").arg(destination.errorString()));
            return;
        }
        m_message->setText(tr("已保存共用 ROI 方案：%1（手动匹配 %2 个，待设置 %3 个）")
                               .arg(path).arg(count - matches.count(-1)).arg(output.count(RoiState::Pending)));
        review.accept();
    });
    selectRow(0);
    refresh();
    review.exec();
}
