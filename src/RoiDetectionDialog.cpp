#include "RoiDetectionDialog.h"
#include "CrystalRoi.h"
#include "HdfFile.h"
#include "WorkflowUi.h"
#include "WorkspacePaths.h"
#include <QApplication>
#include <QGraphicsRectItem>
#include <QGraphicsScene>
#include <QGraphicsView>
#include <QJsonArray>
#include <QJsonDocument>
#include <QJsonObject>
#include <QMessageBox>
#include <QMouseEvent>
#include <QSaveFile>
#include <QSet>
#include <QSettings>
#include <QStandardPaths>
#include <climits>
#include <cmath>
#include <functional>
#include <memory>

namespace {
using WorkspacePaths::Kind;
using WorkspacePaths::Workflow;

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

RoiDetectionDialog::RoiDetectionDialog(const QString &source, QWidget *parent)
    : QWidget(parent), m_browserSource(source), m_source(source) {
    setObjectName("AutomaticRoiWorkspace");
    setWindowTitle(tr("自动识别光斑 ROI"));
    resize(1150, 680);
    m_pages = new QStackedWidget;
    auto root = new QVBoxLayout(this);
    root->setContentsMargins(0, 0, 0, 0);
    root->addWidget(m_pages);
    m_detection = new QWidget;
    m_pages->addWidget(m_detection);
    auto layout = new QVBoxLayout(m_detection);
    WorkflowUi::compact(layout);
    m_parameters = new QWidget;
    auto parameters = new QVBoxLayout(m_parameters);
    parameters->setContentsMargins(0, 0, 0, 0);
    parameters->setSpacing(layout->spacing());
    auto sourceForm = new QFormLayout;
    parameters->addLayout(sourceForm);
    auto fileRow = [this](QFormLayout *form, const QString &label, const QString &value,
                          const QString &filter) {
        auto row = new QHBoxLayout;
        auto edit = new QLineEdit(value);
        auto browse = new QPushButton(tr("浏览…"));
        row->addWidget(edit, 1);
        row->addWidget(browse);
        form->addRow(label, row);
        connect(browse, &QPushButton::clicked, this, [this, edit, label, filter] {
            const auto file = QFileDialog::getOpenFileName(this, label, edit->text(), filter);
            if (!file.isEmpty())
                edit->setText(file);
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
    auto modelSettings = new QWidget;
    auto modelForm = new QFormLayout(modelSettings);
    WorkflowUi::compact(modelForm);
    m_python = fileRow(modelForm, tr("Python 环境"),
                       settings.value("roiDetection/python", defaultPython).toString(), tr("所有文件 (*)"));
    m_python->setObjectName("SpotRoiPython");
    const QDir appDirectory(QCoreApplication::applicationDirPath());
    auto defaultWeights = appDirectory.filePath("models/roi-maskrcnn.pt");
#ifdef Q_OS_LINUX
    if (!QFileInfo::exists(defaultWeights))
        defaultWeights = appDirectory.absoluteFilePath("../share/xri-tools/models/roi-maskrcnn.pt");
#endif
    m_weights =
        fileRow(modelForm, tr("模型权重"), settings.value("roiDetection/weights", defaultWeights).toString(),
                tr("XRI Mask R-CNN (*.pt);;所有文件 (*)"));
    m_weights->setObjectName("SpotRoiWeights");
    auto choices = new QHBoxLayout;
    m_detector = new QComboBox;
    m_mode = new QComboBox;
    m_mode->addItems({tr("全部帧叠加（原始计数）"), tr("单帧（原始计数）")});
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
    choices->addWidget(m_mode);
    choices->addWidget(new QLabel(tr("帧")));
    choices->addWidget(m_frame);
    choices->addStretch();
    sourceForm->addRow(tr("探测器"), choices);
    auto thresholds = new QHBoxLayout;
    thresholds->addWidget(m_score);
    thresholds->addWidget(new QLabel(tr("掩膜阈值")));
    thresholds->addWidget(m_maskThreshold);
    thresholds->addStretch();
    modelForm->addRow(tr("置信度"), thresholds);
    const bool needsModel = m_python->text().isEmpty() || !QFileInfo::exists(m_weights->text());
    parameters->addWidget(WorkflowUi::disclosure(tr("识别参数"), modelSettings, needsModel));
    layout->addWidget(m_parameters);
    auto display = new QSplitter(Qt::Horizontal);
    display->setChildrenCollapsible(false);
    auto previewGroup = new QGroupBox(tr("探测器图像"));
    auto previewLayout = new QVBoxLayout(previewGroup);
    WorkflowUi::compact(previewLayout);
    m_preview = new QLabel(tr("识别后显示编号和掩膜叠加图"));
    m_preview->setAlignment(Qt::AlignCenter);
    m_preview->setSizePolicy(QSizePolicy::Ignored, QSizePolicy::Ignored);
    m_preview->setMinimumWidth(320);
    m_preview->installEventFilter(this);
    previewLayout->addWidget(m_preview);
    display->addWidget(previewGroup);
    m_table = new QTableWidget(0, 10);
    m_table->setObjectName("SpotRoiResults");
    m_table->setHorizontalHeaderLabels({tr("编号"), tr("置信度"), tr("X 质心"), tr("Y 质心"), tr("掩膜像素"),
                                        tr("原始积分"), tr("净积分"), tr("长轴 / px"), tr("短轴 / px"),
                                        tr("角度 / °")});
    m_table->setEditTriggers(QAbstractItemView::NoEditTriggers);
    m_table->setSelectionBehavior(QAbstractItemView::SelectRows);
    m_table->setAlternatingRowColors(true);
    m_table->horizontalHeader()->setSectionResizeMode(QHeaderView::ResizeToContents);
    m_table->verticalHeader()->hide();
    auto resultsGroup = new QGroupBox(tr("识别结果"));
    auto resultsLayout = new QVBoxLayout(resultsGroup);
    WorkflowUi::compact(resultsLayout);
    resultsLayout->addWidget(m_table);
    display->addWidget(resultsGroup);
    display->setSizes({690, 430});
    display->setStretchFactor(0, 1);
    display->setStretchFactor(1, 0);
    layout->addWidget(display, 1);
    auto actions = new QHBoxLayout;
    m_message = new QLabel;
    m_message->setObjectName("SpotRoiStatus");
    m_message->setWordWrap(true);
    m_run = WorkflowUi::button("play", tr("识别光斑"));
    m_run->setObjectName("RunSpotRoiDetection");
    m_cancel = WorkflowUi::button("square", tr("取消识别"));
    m_export = WorkflowUi::button("download", tr("导出掩膜与参数…"));
    m_export->setObjectName("ExportSpotRoiDetection");
    m_matchSaved = WorkflowUi::button("folder-open", tr("载入识别结果…"));
    m_matchSaved->setObjectName("MatchSavedSpotRois");
    m_matchCurrent = WorkflowUi::button("check", tr("匹配晶体"));
    m_matchCurrent->setObjectName("MatchCurrentSpotRois");
    actions->addWidget(m_matchSaved);
    actions->addWidget(m_export);
    actions->addWidget(m_message, 1);
    actions->addWidget(m_run);
    actions->addWidget(m_cancel);
    actions->addWidget(m_matchCurrent);
    layout->addLayout(actions);
    setRunning(false);
    connect(m_detector, &QComboBox::currentIndexChanged, this, [this] { updateFrameRange(); });
    connect(m_mode, &QComboBox::currentIndexChanged, this, [this] { updateFrameRange(); });
    connect(m_run, &QPushButton::clicked, this, [this] { run(); });
    connect(m_cancel, &QPushButton::clicked, this, [this] {
        m_cancelled = true;
        m_process.kill();
    });
    connect(m_export, &QPushButton::clicked, this, [this] { exportResult(); });
    connect(m_matchSaved, &QPushButton::clicked, this, [this] { loadResultForMatching(); });
    connect(m_matchCurrent, &QPushButton::clicked, this, [this] { matchAndSave(m_result); });
    connect(&m_process, &QProcess::finished, this,
            [this](int code, QProcess::ExitStatus status) { completed(code, status); });
    connect(&m_process, &QProcess::errorOccurred, this, [this](QProcess::ProcessError error) {
        if (error == QProcess::FailedToStart) {
            setRunning(false);
            m_message->setText(tr("无法启动 Python。请核对环境路径。%1").arg(m_process.errorString()));
            syncBrowserSource();
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

void RoiDetectionDialog::setBrowserSource(const QString &source) {
    m_browserSource = source;
    syncBrowserSource();
}

bool RoiDetectionDialog::syncBrowserSource() {
    if (m_source == m_browserSource || m_process.state() != QProcess::NotRunning ||
        m_pages->currentWidget() != m_detection)
        return false;
    m_source = m_browserSource;
    m_result.clear();
    m_overlay = {};
    m_preview->setText(tr("识别后显示编号和掩膜叠加图"));
    m_table->setRowCount(0);
    m_frame->setValue(0);
    inspectSource();
    setRunning(false);
    return true;
}

void RoiDetectionDialog::inspectSource() {
    const auto last = m_detector->currentText();
    const QSignalBlocker blocked(m_detector);
    m_detector->clear();
    if (QFileInfo(m_source).isFile()) {
        for (const auto &detector : {QStringLiteral("D_LAMBDA"), QStringLiteral("D_MINIPIX")}) {
            try {
                const auto shape = HdfFile::inspect(m_source, "/entry/data/" + detector).dims;
                if (shape.size() == 3 && shape[0] > 0 && shape[0] <= INT_MAX)
                    m_detector->addItem(detector, int(shape[0]));
            } catch (const std::exception &) {
                // A browser file may contain only one supported detector.
            }
        }
    }
    m_detector->setCurrentIndex(std::max(0, m_detector->findText(last)));
    updateFrameRange();
    m_run->setEnabled(m_detector->count() > 0);
    m_message->setText(m_detector->count()
                           ? tr("%1 · 选择探测器和帧范围后运行识别。").arg(QFileInfo(m_source).fileName())
                           : tr("请在数据浏览器中导入并选中含 LAMBDA 或 MINIPIX 数据的文件。"));
}

bool RoiDetectionDialog::eventFilter(QObject *watched, QEvent *event) {
    if (watched == m_preview && event->type() == QEvent::Resize)
        resizePreview();
    return QWidget::eventFilter(watched, event);
}

void RoiDetectionDialog::resizePreview() {
    if (!m_overlay.isNull())
        m_preview->setPixmap(
            m_overlay.scaled(m_preview->size(), Qt::KeepAspectRatio, Qt::SmoothTransformation));
}

void RoiDetectionDialog::updateFrameRange() {
    m_frame->setMaximum(std::max(0, m_detector->currentData().toInt() - 1));
    m_frame->setEnabled(m_mode->currentIndex() == 1 && m_detector->count() > 0);
}

void RoiDetectionDialog::setRunning(bool running) {
    m_parameters->setEnabled(!running);
    m_run->setEnabled(!running && m_detector->count() > 0);
    m_cancel->setEnabled(running);
    m_export->setEnabled(!running && !m_result.isEmpty());
    m_matchSaved->setEnabled(!running);
    m_matchCurrent->setEnabled(!running && m_table->rowCount() > 0);
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
    m_process.start(
        python, {"-u", script, "infer", "--file", m_source, "--detector", m_detector->currentText(),
                 "--frame", QString::number(m_mode->currentIndex() == 0 ? -1 : m_frame->value()), "--weights",
                 QFileInfo(weights).absoluteFilePath(), "--score", QString::number(m_score->value()),
                 "--mask-threshold", QString::number(m_maskThreshold->value()), "--output",
                 m_work.filePath("spots.json"), "--preview", m_work.filePath("spots.png")});
}

void RoiDetectionDialog::completed(int code, QProcess::ExitStatus status) {
    setRunning(false);
    if (syncBrowserSource())
        return;
    if (m_cancelled) {
        m_preview->setText(tr("识别已取消"));
        m_message->setText(tr("识别已取消。"));
        return;
    }
    if (code != 0 || status != QProcess::NormalExit) {
        const auto detail =
            QString::fromUtf8(m_process.readAllStandardError()).trimmed().section('\n', -1).left(350);
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
        const double values[] = {spot.value("id").toDouble(),
                                 spot.value("score").toDouble(),
                                 center[0].toDouble(),
                                 center[1].toDouble(),
                                 spot.value("mask_area").toDouble(),
                                 spot.value("integrated_counts").toDouble(),
                                 spot.value("net_counts").toDouble(),
                                 spot.value("ellipse_major_px").toDouble(),
                                 spot.value("ellipse_minor_px").toDouble(),
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
    const auto path = WorkspacePaths::saveFile(this, Workflow::Roi, Kind::Result, m_source,
                                               tr("导出掩膜与参数"), "spot-rois.json", tr("JSON (*.json)"));
    if (path.isEmpty())
        return;
    QSaveFile file(path);
    if (file.open(QIODevice::WriteOnly) && file.write(m_result) == m_result.size() && file.commit())
        m_message->setText(tr("已保存掩膜与参数到 %1").arg(path));
    else
        m_message->setText(tr("保存失败。%1").arg(file.errorString()));
}

void RoiDetectionDialog::loadResultForMatching() {
    const auto path = QFileDialog::getOpenFileName(this, tr("载入自动识别结果"),
                                                   WorkspacePaths::directory(Workflow::Roi, Kind::Result),
                                                   tr("JSON (*.json)"));
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

    QVector<QVector<QPoint>> masks;
    QVector<QPointF> centers;
    QVector<int> ids;
    QMultiHash<QPoint, int> hits;
    for (const auto &value : spots) {
        const auto spot = value.toObject();
        auto mask = CrystalRois::decodeMask(spot["mask_runs"].toArray(), size);
        const auto center = spot["centroid_xy"].toArray();
        const int id = spot["id"].toInt();
        if (mask.isEmpty() || id <= 0 || ids.contains(id) || center.size() != 2 || !center[0].isDouble() ||
            !center[1].isDouble()) {
            m_message->setText(tr("识别结果中的光斑编号、质心或掩膜无效。"));
            return;
        }
        for (const auto &p : mask)
            hits.insert(p, int(masks.size()));
        ids.append(id);
        centers.append({center[0].toDouble(), center[1].toDouble()});
        masks.append(std::move(mask));
    }

    const int first = d == 0 ? 0 : 60, count = d == 0 ? 60 : 15;
    struct MatchingState {
        CrystalRois output;
        QVector<QSize> sizes;
        QVector<int> matches;
    };
    const auto state =
        std::make_shared<MatchingState>(MatchingState{CrystalRois{}, sourceSizes, QVector<int>(count, -1)});
    auto review = new QWidget;
    auto layout = new QVBoxLayout(review);
    WorkflowUi::compact(layout);
    auto info = new QLabel(tr("先在右侧选择晶体，再点击左侧对应的光斑掩膜。滚轮缩放，拖动平移。"));
    info->setWordWrap(true);

    auto split = new QSplitter(Qt::Horizontal);
    split->setChildrenCollapsible(false);
    auto table = new QTableWidget(count, 3);
    table->setObjectName("SpotRoiCrystalMatches");
    table->setHorizontalHeaderLabels({tr("晶体"), tr("状态"), tr("光斑编号")});
    table->setEditTriggers(QAbstractItemView::NoEditTriggers);
    table->setSelectionBehavior(QAbstractItemView::SelectRows);
    table->setSelectionMode(QAbstractItemView::ExtendedSelection);
    table->setAlternatingRowColors(true);
    table->verticalHeader()->hide();
    table->horizontalHeader()->setSectionResizeMode(QHeaderView::ResizeToContents);
    table->horizontalHeader()->setStretchLastSection(true);
    for (int row = 0; row < count; ++row) {
        table->setItem(row, 0, new QTableWidgetItem(state->output.items()[first + row].label));
        table->setItem(row, 1, new QTableWidgetItem);
        table->setItem(row, 2, new QTableWidgetItem);
    }
    auto scene = new QGraphicsScene(review);
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
    auto previewGroup = new QGroupBox(detector);
    auto previewLayout = new QVBoxLayout(previewGroup);
    WorkflowUi::compact(previewLayout);
    previewLayout->addWidget(view);
    split->addWidget(previewGroup);
    auto panel = new QWidget;
    auto panelLayout = new QVBoxLayout(panel);
    panelLayout->setContentsMargins(0, 0, 0, 0);
    panelLayout->addWidget(table, 1);
    auto candidates = new QComboBox;
    candidates->setObjectName("OverlappingSpotCandidates");
    candidates->setPlaceholderText(tr("选择此处的光斑"));
    candidates->hide();
    panelLayout->addWidget(candidates);
    split->addWidget(panel);
    split->setSizes({690, 430});
    split->setStretchFactor(0, 1);
    split->setStretchFactor(1, 0);

    auto controls = new QHBoxLayout;
    auto load = WorkflowUi::button("folder-open", tr("载入已有方案…"));
    load->setObjectName("LoadSpotRoiBaseline");
    auto invalid = new QPushButton(tr("标记失效"));
    auto clear = new QPushButton(tr("清除所选 ROI"));
    auto reset = WorkflowUi::button("rotate-ccw", tr("复位视图"));
    auto advance = new QCheckBox(tr("设置后跳到下一个未设置晶体"));
    advance->setChecked(true);
    auto rowButtons = new QHBoxLayout;
    rowButtons->addWidget(invalid);
    rowButtons->addWidget(clear);
    panelLayout->addLayout(rowButtons);
    controls->addWidget(reset);
    controls->addWidget(advance);
    controls->addWidget(info, 1);
    layout->addLayout(controls);
    layout->addWidget(split, 1);
    auto actions = new QHBoxLayout;
    auto status = new QLabel(tr("可从空白方案开始，也可载入已有方案继续编辑。"));
    status->setWordWrap(true);
    auto save = WorkflowUi::button("save", tr("保存ROI 方案…"));
    save->setObjectName("SaveMatchedSpotRois");
    auto back = WorkflowUi::button("chevron-left", tr("返回识别结果"));
    back->setObjectName("ReturnToSpotRoiResults");
    actions->addWidget(load);
    actions->addWidget(save);
    actions->addWidget(status, 1);
    actions->addWidget(back);
    layout->addLayout(actions);

    auto selectedRows = [table] { return table->selectionModel()->selectedRows(); };
    auto selectRow = [table](int row) {
        table->setCurrentCell(row, 0, QItemSelectionModel::ClearAndSelect | QItemSelectionModel::Rows);
        table->scrollToItem(table->item(row, 0));
    };
    auto showSelection = [selectedRows, outline, state, first] {
        const auto selected = selectedRows();
        outline->hide();
        if (selected.size() == 1) {
            const auto &roi = state->output.items()[first + selected.first().row()];
            if (roi.state == RoiState::Valid) {
                outline->setRect(QRectF(roi.rect));
                outline->show();
            }
        }
    };
    auto refresh = [this, state, first, count, base, hasPreview, masks, centers, ids, labels, image, table,
                    save, showSelection] {
        const auto &matches = state->matches;
        const auto &output = state->output;
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
            table->item(row, 1)->setText(roi.state == RoiState::Pending   ? tr("未设置")
                                         : roi.state == RoiState::Invalid ? tr("已失效")
                                         : roi.mask.isEmpty()             ? tr("矩形 ROI")
                                                                          : tr("掩膜"));
            table->item(row, 2)->setText(matches[row] < 0 ? QString() : QString::number(ids[matches[row]]));
        }
        save->setEnabled(matches.count(-1) < count);
        showSelection();
    };
    auto advanceFrom = [advance, state, first, count, selectRow](int row) {
        if (!advance->isChecked())
            return;
        for (int step = 1; step < count; ++step) {
            const int next = (row + step) % count;
            if (state->output.items()[first + next].state == RoiState::Pending) {
                selectRow(next);
                break;
            }
        }
    };
    auto finishReview = [this, review] {
        m_pages->setCurrentWidget(m_detection);
        m_pages->removeWidget(review);
        review->deleteLater();
        syncBrowserSource();
    };
    connect(table, &QTableWidget::itemSelectionChanged, review, [showSelection, candidates] {
        candidates->hide();
        showSelection();
    });
    connect(reset, &QPushButton::clicked, review, [view] { view->resetZoom(); });
    connect(back, &QPushButton::clicked, review, finishReview);
    connect(invalid, &QPushButton::clicked, review,
            [selectedRows, state, first, refresh, advanceFrom, candidates] {
                candidates->hide();
                const auto rows = selectedRows();
                for (const auto &row : rows) {
                    state->output.invalidate(first + row.row());
                    state->matches[row.row()] = -1;
                }
                refresh();
                if (rows.size() == 1)
                    advanceFrom(rows.first().row());
            });
    connect(clear, &QPushButton::clicked, review, [selectedRows, state, first, refresh, candidates] {
        candidates->hide();
        for (const auto &row : selectedRows()) {
            state->output.clear(first + row.row());
            state->matches[row.row()] = -1;
        }
        refresh();
    });
    auto matchSpot = [this, selectedRows, state, first, status, ids, masks, size, refresh,
                      advanceFrom](int spot) {
        const auto selected = selectedRows();
        if (selected.size() != 1)
            return;
        auto &matches = state->matches;
        auto &output = state->output;
        const int row = selected.first().row(), owner = int(matches.indexOf(spot));
        if (owner >= 0 && owner != row) {
            status->setText(tr("光斑 %1 已匹配 %2。请先清除该晶体的匹配。")
                                .arg(ids[spot])
                                .arg(output.items()[first + owner].label));
            return;
        }
        if (!output.assignMask(first + row, masks[spot], size))
            return;
        matches[row] = spot;
        status->setText(tr("%1 已匹配光斑 %2。").arg(output.items()[first + row].label).arg(ids[spot]));
        refresh();
        advanceFrom(row);
    };
    connect(candidates, &QComboBox::activated, review, [candidates, matchSpot](int index) {
        candidates->hide();
        matchSpot(candidates->itemData(index).toInt());
    });
    view->pick = [this, selectedRows, status, hits, ids, candidates, matchSpot](const QPoint &pixel) {
        candidates->hide();
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
        if (clicked.size() == 1) {
            matchSpot(clicked.first());
            return;
        }
        candidates->clear();
        for (const int i : clicked)
            candidates->addItem(tr("光斑 %1").arg(ids[i]), i);
        candidates->setCurrentIndex(-1);
        candidates->show();
    };
    connect(load, &QPushButton::clicked, review,
            [this, review, state, count, source, sourceSizes, status, refresh, candidates] {
                auto &matches = state->matches;
                if (matches.count(-1) < count &&
                    QMessageBox::question(review, tr("载入已有方案"),
                                          tr("载入会替换本窗口尚未保存的匹配。继续吗？")) != QMessageBox::Yes)
                    return;
                auto folder = WorkspacePaths::directory(Workflow::Roi, Kind::Config);
                if (folder.isEmpty())
                    folder = QFileInfo(source).absolutePath();
                const auto path = QFileDialog::getOpenFileName(review, tr("载入已有 ROI 方案"), folder,
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
                if (!state->output.restore(bytes, imported, &error)) {
                    status->setText(error);
                    return;
                }
                state->sizes = imported;
                matches.fill(-1);
                candidates->hide();
                refresh();
                status->setText(tr("已载入原方案。请选择晶体并点击光斑掩膜。"));
            });
    connect(save, &QPushButton::clicked, review,
            [this, review, state, first, count, selectRow, status, source, finishReview] {
                const auto &output = state->output;
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
                const auto path = WorkspacePaths::saveFile(
                    review, Workflow::Roi, Kind::Config, source, tr("保存ROI 方案"), "matched-spots.roi.json",
                    tr("XRI ROI (*.roi.json)"),
                    QFileInfo(source).absoluteDir().filePath("matched-spots.roi.json"));
                if (path.isEmpty())
                    return;
                if (!path.endsWith(".roi.json", Qt::CaseInsensitive)) {
                    status->setText(tr("ROI 方案的文件名应以 .roi.json 结尾。"));
                    return;
                }
                const auto scheme = output.serialize(source, state->sizes);
                QSaveFile destination(path);
                if (!destination.open(QIODevice::WriteOnly) || destination.write(scheme) != scheme.size() ||
                    !destination.commit()) {
                    status->setText(tr("无法保存方案。%1").arg(destination.errorString()));
                    return;
                }
                m_message->setText(tr("已保存ROI 方案：%1（手动匹配 %2 个，待设置 %3 个）")
                                       .arg(path)
                                       .arg(count - state->matches.count(-1))
                                       .arg(output.count(RoiState::Pending)));
                finishReview();
            });
    selectRow(0);
    refresh();
    m_pages->addWidget(review);
    m_pages->setCurrentWidget(review);
}
