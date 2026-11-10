#pragma once
#include "ImagingIO.h"
#include "WorkflowUi.h"
#include <QJsonDocument>
#include <QRegularExpression>
#include <algorithm>
#include <cmath>
#include <functional>

namespace ImagingInputUi {
inline QJsonValue read(QWidget *widget, const QString &key) {
    if (auto box = qobject_cast<QCheckBox *>(widget))
        return box->isChecked();
    else if (auto combo = qobject_cast<QComboBox *>(widget))
        return combo->currentData().toString();
    else if (auto spin = qobject_cast<QSpinBox *>(widget)) {
        const bool automatic =
            spin->property("automaticLoopCount").toBool() && spin->value() == spin->minimum();
        return automatic ? QJsonValue() : QJsonValue(spin->value());
    } else if (auto edit = qobject_cast<QLineEdit *>(widget)) {
        auto text = edit->text().trimmed();
        int kind = edit->property("fieldKind").toInt();
        if (kind == 1) {
            bool ok;
            double v = text.toDouble(&ok);
            if (!ok || !std::isfinite(v))
                ImageUtil::fail("Invalid numeric parameter: " + key);
            return v;
        } else if (kind == 2) {
            QJsonArray ids;
            for (const auto &id : text.split(QRegularExpression("[,，;；\\s]+"), Qt::SkipEmptyParts))
                ids.append(id);
            return ids;
        } else if (kind == 3 || kind == 4) {
            QJsonParseError error;
            auto doc = QJsonDocument::fromJson(text.toUtf8(), &error);
            if (error.error != QJsonParseError::NoError || (kind == 3 && !doc.isObject()) ||
                (kind == 4 && !doc.isArray()))
                ImageUtil::fail("Invalid JSON parameter: " + key);
            return kind == 3 ? QJsonValue(doc.object()) : QJsonValue(doc.array());
        } else
            return text;
    }
    return {};
}

inline void write(QWidget *widget, const QJsonValue &value) {
    if (auto box = qobject_cast<QCheckBox *>(widget))
        box->setChecked(value.toBool());
    else if (auto combo = qobject_cast<QComboBox *>(widget))
        combo->setCurrentIndex(std::max(0, combo->findData(value.toString())));
    else if (auto spin = qobject_cast<QSpinBox *>(widget))
        spin->setValue(value.toInt(spin->minimum()));
    else if (auto edit = qobject_cast<QLineEdit *>(widget)) {
        int kind = edit->property("fieldKind").toInt();
        if (kind == 1)
            edit->setText(value.isNull() ? "0" : QString::number(value.toDouble(), 'g', 15));
        else if (kind == 2) {
            QStringList ids;
            for (auto id : value.toArray())
                ids.append(id.toString());
            edit->setText(ids.join(','));
        } else if (kind == 3)
            edit->setText(QString::fromUtf8(QJsonDocument(value.toObject()).toJson(QJsonDocument::Compact)));
        else if (kind == 4)
            edit->setText(QString::fromUtf8(QJsonDocument(value.toArray()).toJson(QJsonDocument::Compact)));
        else
            edit->setText(value.toString());
    }
}

template <class AddText, class AddNumber>
void build(QWidget *context, QVBoxLayout *ll, QVBoxLayout *rr, QMap<QString, QWidget *> &fields,
           const QStringList &availableFiles, const AddText &addText, const AddNumber &addNumber,
           const std::function<void()> &changed, const std::function<void(const QString &)> &roiRequested) {
    auto files = new QGroupBox(QObject::tr("原始文件与共用 ROI"));
    auto ff = new QFormLayout(files);
    auto fileField = [&](const QString &key, const QString &label, bool scan) {
        auto row = new QWidget;
        auto l = new QHBoxLayout(row);
        l->setContentsMargins(0, 0, 0, 0);
        auto edit = new QLineEdit;
        edit->setObjectName(key);
        fields[key] = edit;
        auto choose = new QPushButton(QObject::tr("浏览…"));
        l->addWidget(edit, 1);
        l->addWidget(choose);
        ff->addRow(label, row);
        QObject::connect(edit, &QLineEdit::textChanged, context, changed);
        QObject::connect(
            choose, &QPushButton::clicked, context, [context, &availableFiles, edit, label, scan] {
                const auto title = QObject::tr("选择%1").arg(label);
                const auto path =
                    scan
                        ? WorkflowUi::selectSource(context, availableFiles, edit->text(), title)
                        : QFileDialog::getOpenFileName(context, title, QFileInfo(edit->text()).absolutePath(),
                                                       "ROI (*.json *.txt *.tsv);;All files (*)");
                if (!path.isEmpty())
                    edit->setText(path);
            });
        return l;
    };
    fileField("input.calibration_nxs", QObject::tr("薄膜标定"), true);
    fileField("input.sample_nxs", QObject::tr("样品成像"), true);
    auto roiButtons = fileField("input.roi_table", QObject::tr("ROI 方案 / TSV"), false);
    auto selectRoi = new QPushButton(QObject::tr("选取共用 ROI…"));
    roiButtons->addWidget(selectRoi);
    auto roiNote =
        new QLabel(QObject::tr("按共用方案的逐像素掩膜或矩形取固定像素。旧 TSV 的平移/扩展按表应用。共用 "
                               "ROI 中的失效晶体保持排除；仅使用所选探测器。"));
    roiNote->setWordWrap(true);
    ff->addRow(roiNote);
    ll->addWidget(files);
    selectRoi->setObjectName("SelectImagingInputRoi");
    QObject::connect(selectRoi, &QPushButton::clicked, context,
                     [roiRequested, sample = fields["input.sample_nxs"]] {
                         roiRequested(static_cast<QLineEdit *>(sample)->text().trimmed());
                     });
    auto acq = new QGroupBox(QObject::tr("采集与尺度（沿用原成像归一化）"));
    auto af = new QFormLayout(acq);
    addText(af, "acquisition.detector", QObject::tr("探测器"));
    addText(af, "acquisition.monitor", QObject::tr("入射监测通道"));
    addNumber(af, "acquisition.calibration_exposure_s", QObject::tr("薄膜曝光 / s"));
    addNumber(af, "acquisition.sample_exposure_s", QObject::tr("样品曝光 / s"));
    addNumber(af, "acquisition.foil_thickness_mm", QObject::tr("薄膜厚度 / mm"));
    addNumber(af, "acquisition.detector_pixel_pitch_mm", QObject::tr("探测器像素间距 / mm"));
    auto scanLoops = new QWidget;
    auto loopLayout = new QVBoxLayout(scanLoops);
    loopLayout->setContentsMargins(0, 0, 0, 0);
    auto loopTabs = new QTabWidget;
    loopLayout->addWidget(loopTabs);
    auto addLoops = [context, changed, &fields, loopTabs](const QString &role, const QString &label,
                                                          const QStringList &axes) {
        auto page = new QWidget;
        auto form = new QFormLayout(page);
        auto order = new QComboBox;
        const auto orderKey = "scan_loops." + role + "_order";
        order->setObjectName(orderKey);
        order->addItem(QObject::tr("自动（按文件计划）"), QString());
        auto permutation = axes;
        std::sort(permutation.begin(), permutation.end());
        do {
            order->addItem(permutation.join(QStringLiteral(" → ")), permutation.join(','));
        } while (std::next_permutation(permutation.begin(), permutation.end()));
        fields[orderKey] = order;
        form->addRow(QObject::tr("循环顺序（外 → 内）"), order);
        QObject::connect(order, &QComboBox::currentIndexChanged, context, changed);
        for (const auto &axis : axes) {
            auto points = new QSpinBox;
            const auto key = "scan_loops." + role + '_' + axis.toLower() + "_points";
            points->setObjectName(key);
            points->setRange(1, 100000);
            points->setSpecialValueText(QObject::tr("自动"));
            points->setKeyboardTracking(false);
            points->setProperty("automaticLoopCount", true);
            fields[key] = points;
            form->addRow(QObject::tr("%1 层帧数（扫描点）").arg(axis), points);
            QObject::connect(points, &QSpinBox::valueChanged, context, changed);
        }
        loopTabs->addTab(page, label);
    };
    addLoops("sample", QObject::tr("样品成像"), {"Y", "Z", "DCM"});
    addLoops("calibration", QObject::tr("薄膜标定"), {"X", "DCM"});
    auto loopNote = new QLabel(QObject::tr(
        "手动顺序需与文件的扫描计划一致。各层帧数可单独指定，程序按文件扫描范围和逐帧坐标校验。"));
    loopNote->setWordWrap(true);
    loopLayout->addWidget(loopNote);
    af->addRow(WorkflowUi::disclosure(QObject::tr("扫描循环与各层帧数"), scanLoops));
    addNumber(af, "runtime.chunk_frames", QObject::tr("每批读取帧数（缓存）"), true);
    auto normNote = new QLabel(QObject::tr(
        "每个网格单元使用 Σ计数 / Σ(I0 × 曝光)。I0 必须为正，不取绝对值，也不应用谱学的异常点替换。"));
    normNote->setWordWrap(true);
    af->addRow(normNote);
    ll->addWidget(WorkflowUi::disclosure(QObject::tr("采集与物理尺度"), acq, true));
    auto rois = new QGroupBox(QObject::tr("ROI 编号与标签"));
    auto rf = new QFormLayout(rois);
    addText(rf, "roi.excluded_ids", QObject::tr("排除原编号（逗号分隔）"), 2);
    addText(rf, "roi.label_overrides", QObject::tr("标签修正（JSON）"), 3);
    auto labels =
        new QLabel(QObject::tr("例如排除 31,44，标签修正 {\"35\":\"VU-A3\"}。增益锚点使用 原编号_标签，例如 "
                               "36_VU-C1。默认不猜测占位框或自动更换锚点。"));
    labels->setWordWrap(true);
    rf->addRow(labels);
    rr->addWidget(WorkflowUi::disclosure(QObject::tr("高级：ROI 排除与标签修正"), rois));
}
} // namespace ImagingInputUi
