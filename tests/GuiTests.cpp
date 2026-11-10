#include "BondWorkspace.h"
#include "ImagingDialog.h"
#include "ImagingFixture.h"
#include "MainWindow.h"
#include "PlotStyle.h"
#include "RoiDetectionDialog.h"
#include "RoiDialog.h"
#include "ScientificPlot.h"
#include "SpectroscopyDialog.h"
#include "SpectrumFixture.h"
#include "SpectrumPlot.h"
#include "SpectrumProcessing.h"
#include <QAction>
#include <QApplication>
#include <QCheckBox>
#include <QComboBox>
#include <QCryptographicHash>
#include <QDataStream>
#include <QDir>
#include <QDoubleSpinBox>
#include <QElapsedTimer>
#include <QFileDialog>
#include <QFileInfo>
#include <QFontDatabase>
#include <QGraphicsView>
#include <QJsonArray>
#include <QJsonDocument>
#include <QJsonObject>
#include <QLabel>
#include <QLineEdit>
#include <QMenu>
#include <QMenuBar>
#include <QMouseEvent>
#include <QPointer>
#include <QPushButton>
#include <QSettings>
#include <QSignalSpy>
#include <QSpinBox>
#include <QStackedWidget>
#include <QTabWidget>
#include <QTableWidget>
#include <QTemporaryDir>
#include <QTimer>
#include <QToolButton>
#include <QTranslator>
#include <QTreeWidget>
#include <QWheelEvent>
#include <QtTest>
#include <cmath>

class GuiTests : public QObject {
    Q_OBJECT
  private slots:
    void bilingualUi() {
        const QString key = "ui/language";
        QSettings settings;
        const auto previous = settings.value(key);
        settings.remove(key);
        {
            QTemporaryDir temp;
            const auto file = QFileInfo(SpectrumFixture::create(temp.path(), true)).canonicalFilePath();
            QStringList restartArguments;
            int launches = 0;
            bool allowRestart = false;
            MainWindow chinese(nullptr, [&](const QStringList &arguments) {
                ++launches;
                restartArguments = arguments;
                return allowRestart;
            });
            chinese.show();
            auto *menu = chinese.findChild<QMenu *>("LanguageMenu");
            auto *english = chinese.findChild<QAction *>("Language_en");
            auto *simplified = chinese.findChild<QAction *>("Language_zh_CN");
            QVERIFY(menu && english && simplified);
            QCOMPARE(menu->title(), QStringLiteral("语言(&L)"));
            chinese.openFile(file);
            chinese.selectPath("/entry/data/D_LAMBDA");
            chinese.setSliceIndex(0, 1);
            QVERIFY(chinese.isReading());
            english->trigger();
            QCOMPARE(launches, 0);
            QVERIFY(simplified->isChecked());
            QTRY_VERIFY(!chinese.isReading());
            QTimer::singleShot(0, [] {
                if (auto *message = qobject_cast<QDialog *>(QApplication::activeModalWidget()))
                    message->accept();
            });
            english->trigger();
            QCOMPARE(launches, 1);
            QVERIFY(simplified->isChecked() && chinese.isVisible());
            QVERIFY(!QSettings().contains(key));
            allowRestart = true;
            english->trigger();
            QCOMPARE(launches, 2);
            QVERIFY(!chinese.isVisible());
            QCOMPARE(QSettings().value(key).toString(), QStringLiteral("en"));
            QCOMPARE(restartArguments, (QStringList{"--dataset", "/entry/data/D_LAMBDA", "--frame", "1", "--", file}));
        }
        QTranslator translator;
        QVERIFY(translator.load(":/i18n/xri_en.qm"));
        QVERIFY(QCoreApplication::installTranslator(&translator));
        {
            int launches = 0;
            MainWindow english(nullptr, [&](const QStringList &) {
                ++launches;
                return true;
            });
            QCOMPARE(english.findChild<QMenu *>("FileMenu")->title(), QStringLiteral("&File"));
            QCOMPARE(english.findChild<QMenu *>("LanguageMenu")->title(), QStringLiteral("&Language"));
            QVERIFY(english.findChild<QAction *>("Language_en")->isChecked());
            auto *spectroscopy = english.openSpectroscopy()->findChild<QTabWidget *>("SpectroscopySteps");
            auto *imaging = english.openImaging()->findChild<QTabWidget *>("ImagingSteps");
            auto *bond = english.openBondImaging()->findChild<QTabWidget *>("BondSteps");
            QVERIFY(spectroscopy && imaging && bond);
            QCOMPARE(spectroscopy->tabText(3), QStringLiteral("4 XRS and I0"));
            QCOMPARE(imaging->tabText(2), QStringLiteral("3 Input and Calibration Checks"));
            QCOMPARE(bond->tabText(1), QStringLiteral("2 Spectral Processing and Spatial Reconstruction"));
            RoiDetectionDialog masks({});
            QCOMPARE(masks.windowTitle(), QStringLiteral("Automatic Spot ROI Detection"));
            CrystalRois rois;
            QString error;
            QVERIFY(!rois.restore("{}", {}, &error));
            QCOMPARE(error, QStringLiteral("This is not a supported XRI ROI scheme."));
            auto *chinese = english.findChild<QAction *>("Language_zh_CN");
            chinese->trigger();
            QCOMPARE(launches, 1);
            QCOMPARE(QSettings().value(key).toString(), QStringLiteral("zh_CN"));
        }
        QCoreApplication::removeTranslator(&translator);
        MainWindow chineseAgain;
        QCOMPARE(chineseAgain.findChild<QMenu *>("FileMenu")->title(), QStringLiteral("文件(&F)"));
        if (previous.isValid())
            settings.setValue(key, previous);
        else
            settings.remove(key);
    }
    void automaticSpotRoi() {
        const auto source = qEnvironmentVariable("XRI_ROI_ML_TEST_FILE");
        const auto python = qEnvironmentVariable("XRI_ROI_ML_PYTHON");
        const auto weights = qEnvironmentVariable("XRI_ROI_ML_WEIGHTS");
        if (source.isEmpty() || python.isEmpty() || weights.isEmpty())
            QSKIP("Optional Mask R-CNN inference requires a scan, Python environment and weights");
        MainWindow window;
        window.show();
        window.openFile(source);
        QTRY_VERIFY_WITH_TIMEOUT(!window.isReading(), 10000);
        auto action = window.findChild<QAction *>("DetectSpotRois");
        QVERIFY(action);
        action->trigger();
        auto run = window.findChild<QPushButton *>("RunSpotRoiDetection");
        auto exportButton = window.findChild<QPushButton *>("ExportSpotRoiDetection");
        auto pythonEdit = window.findChild<QLineEdit *>("SpotRoiPython");
        auto weightsEdit = window.findChild<QLineEdit *>("SpotRoiWeights");
        auto results = window.findChild<QTableWidget *>("SpotRoiResults");
        auto status = window.findChild<QLabel *>("SpotRoiStatus");
        QVERIFY(run && exportButton && pythonEdit && weightsEdit && results && status);
        QVERIFY(run->isEnabled() && !exportButton->isEnabled());
        pythonEdit->setText(python);
        weightsEdit->setText(weights);
        run->click();
        QTRY_VERIFY_WITH_TIMEOUT(run->isEnabled(), 60000);
        QVERIFY2(exportButton->isEnabled(), qPrintable(status->text()));
        QVERIFY(results->rowCount() > 0);
        QVERIFY(results->item(0, 4)->text().toInt() >= 8);
        QVERIFY(results->item(0, 5)->text().toDouble() > 0.);
        const auto screenshot = qEnvironmentVariable("XRI_ROI_ML_SCREENSHOT");
        if (!screenshot.isEmpty())
            QVERIFY(run->window()->grab().save(screenshot));
    }
    void matchSavedSpotRois() {
        QTemporaryDir temp;
        const auto source = SpectrumFixture::create(temp.path(), true);
        auto spot = [](int id, int x, int y, const QJsonArray &runs) {
            return QJsonObject{{"id", id}, {"centroid_xy", QJsonArray{x, y}}, {"mask_runs", runs}};
        };
        const auto report = QJsonDocument(QJsonObject{
            {"schema", "xri.roi.maskrcnn.result.v1"}, {"source", source},
            {"detector", "D_LAMBDA"}, {"shape_yx", QJsonArray{8, 8}},
            {"spots", QJsonArray{spot(1, 1, 1, QJsonArray{QJsonArray{9, 1}, QJsonArray{18, 1}}),
                                 spot(2, 5, 5, QJsonArray{QJsonValue(QJsonArray{45, 1})}),
                                 spot(3, 6, 5, QJsonArray{QJsonValue(QJsonArray{46, 1})})}}}).toJson();
        QFile detection(temp.filePath("detected.json"));
        QVERIFY(detection.open(QIODevice::WriteOnly));
        QCOMPARE(detection.write(report), report.size());
        detection.close();
        const auto output = temp.filePath("shared.roi.json");
        const QStringList paths{detection.fileName(), temp.filePath("test.roi.json"), output};
        RoiDetectionDialog dialog(source);
        dialog.show();
        auto action = dialog.findChild<QAction *>("MatchSavedSpotRois");
        QVERIFY(action && action->isEnabled());
        int pickedFiles = 0;
        bool reviewed = false;
        bool startedEmpty = false;
        bool pickedManually = false;
        bool requestedBaseline = false;
        QPointer<QFileDialog> activePicker;
        QTimer driveDialogs;
        connect(&driveDialogs, &QTimer::timeout, &dialog, [&] {
            auto *modal = qobject_cast<QDialog *>(QApplication::activeModalWidget());
            if (auto *picker = qobject_cast<QFileDialog *>(modal)) {
                if (activePicker != picker) {
                    if (pickedFiles >= paths.size()) {
                        picker->reject();
                        return;
                    }
                    activePicker = picker;
                    const QFileInfo file(paths[pickedFiles++]);
                    picker->setDirectory(file.absolutePath());
                    picker->selectFile(file.fileName());
                    if (auto *name = picker->findChild<QLineEdit *>("fileNameEdit"))
                        name->setText(file.fileName());
                    return;
                }
                QMetaObject::invokeMethod(picker, "accept");
            } else if (modal && !reviewed) {
                auto *table = modal->findChild<QTableWidget *>("SpotRoiCrystalMatches");
                auto *save = modal->findChild<QPushButton *>("SaveMatchedSpotRois");
                auto *load = modal->findChild<QPushButton *>("LoadSpotRoiBaseline");
                auto *view = modal->findChild<QGraphicsView *>("SpotRoiMaskView");
                if (!table || !save || !load || !view) {
                    modal->reject();
                    return;
                }
                if (!requestedBaseline) {
                    startedEmpty = table->item(0, 2)->text().isEmpty() && !save->isEnabled();
                    requestedBaseline = true;
                    QTimer::singleShot(0, load, &QPushButton::click);
                    return;
                }
                reviewed = true;
                if (!startedEmpty || pickedFiles != 2 || save->isEnabled() ||
                    !table->item(0, 2)->text().isEmpty() || !table->item(1, 2)->text().isEmpty()) {
                    modal->reject();
                    return;
                }
                table->selectRow(0);
                QTest::mouseClick(view->viewport(), Qt::LeftButton, Qt::NoModifier, view->mapFromScene(6.5, 5.5));
                table->selectRow(1);
                QTest::mouseClick(view->viewport(), Qt::LeftButton, Qt::NoModifier, view->mapFromScene(1.5, 1.5));
                pickedManually = table->item(0, 2)->text() == "3" && table->item(1, 2)->text() == "1";
                if (!pickedManually || !save->isEnabled()) {
                    modal->reject();
                    return;
                }
                QTimer::singleShot(0, save, &QPushButton::click);
            } else if (modal && !modal->findChild<QTableWidget *>("SpotRoiCrystalMatches")) {
                modal->reject();
            }
        });
        QTimer deadline;
        connect(&deadline, &QTimer::timeout, &dialog, [] {
            if (auto *modal = qobject_cast<QDialog *>(QApplication::activeModalWidget()))
                modal->reject();
        });
        const bool native = !QApplication::testAttribute(Qt::AA_DontUseNativeDialogs);
        QApplication::setAttribute(Qt::AA_DontUseNativeDialogs);
        deadline.start(5000);
        driveDialogs.start(20);
        action->trigger();
        driveDialogs.stop();
        deadline.stop();
        QApplication::setAttribute(Qt::AA_DontUseNativeDialogs, !native);
        const auto *status = dialog.findChild<QLabel *>("SpotRoiStatus");
        QVERIFY2(reviewed && pickedManually,
                 qPrintable(QString("pickedFiles=%1, reviewed=%2, status=%3")
                                .arg(pickedFiles).arg(reviewed).arg(status ? status->text() : QString())));
        QCOMPARE(pickedFiles, paths.size());
        QFile saved(output);
        QVERIFY(saved.open(QIODevice::ReadOnly));
        CrystalRois scheme;
        QString error;
        QVERIFY2(scheme.restore(saved.readAll(), {{8, 8}, {6, 6}}, &error), qPrintable(error));
        QCOMPARE(scheme.items()[0].mask, (QVector<QPoint>{{6, 5}}));
        QCOMPARE(scheme.items()[1].mask, (QVector<QPoint>{{1, 1}, {2, 2}}));
        QCOMPARE(scheme.items()[60].rect, QRect(1, 1, 2, 2));
        QVERIFY(scheme.items()[60].mask.isEmpty());
    }
    void bondWorkflow() {
        QTemporaryDir temp;
        MainWindow window;
        window.resize(1200, 800);
        window.show();
        auto open = window.findChild<QAction *>("OpenBondWorkspace");
        QVERIFY(open);
        open->trigger();
        auto work = window.findChild<BondWorkspace *>("BondWorkspace");
        QVERIFY(work && work->isVisible());
        auto tabs = work->findChild<QTabWidget *>("BondSteps");
        QVERIFY(tabs && !tabs->isTabEnabled(BondWorkspace::Spectra));
        const auto dir = qEnvironmentVariable("XRI_BOND_SCREENSHOTS");
        const auto external = qEnvironmentVariable("XRI_BOND_TEST_VOLUME");
        if (external.isEmpty()) {
            work->applyConfiguration(ImagingFixture::create(temp.path(), false, true, true, true));
            auto sampleOrder = work->findChild<QComboBox *>("scan_loops.sample_order");
            auto sampleY = work->findChild<QSpinBox *>("scan_loops.sample_y_points");
            auto diagnosticPng = work->findChild<QCheckBox *>("output.plots");
            auto diagnostics = work->findChild<QCheckBox *>("output.save_diagnostics");
            auto output = work->findChild<QLineEdit *>("BondOutput");
            QVERIFY(sampleOrder && sampleY && diagnosticPng && diagnostics && output);
            sampleOrder->setCurrentIndex(sampleOrder->findData("DCM,Z,Y"));
            sampleY->setValue(2);
            diagnosticPng->setChecked(true);
            diagnostics->setChecked(true);
            output->setText(temp.filePath("result.h5"));
            const auto configPath = temp.filePath("bond-config.json");
            QVERIFY(work->saveConfiguration(configPath));
            sampleOrder->setCurrentIndex(0);
            sampleY->setValue(1);
            diagnosticPng->setChecked(false);
            output->clear();
            QVERIFY(work->importConfiguration(configPath));
            QCOMPARE(sampleOrder->currentData().toString(), QStringLiteral("DCM,Z,Y"));
            QCOMPARE(sampleY->value(), 2);
            QCOMPARE(output->text(), temp.filePath("result.h5"));
            QVERIFY(diagnosticPng->isChecked() && diagnostics->isChecked());
            if (!dir.isEmpty()) {
                QVERIFY(QDir().mkpath(dir));
                QVERIFY(window.grab().save(QDir(dir).filePath("files-and-roi.png")));
            }
            auto validate = work->findChild<QPushButton *>("CheckBondInputs");
            auto run = work->findChild<QPushButton *>("RunBondVolume");
            QVERIFY(validate && run && !run->isVisible());
            validate->click();
            QTRY_VERIFY_WITH_TIMEOUT(!work->isBusy(), 10000);
            QCOMPARE(tabs->currentIndex(), BondWorkspace::Spectra);
            QVERIFY(run->isVisible() && run->isEnabled());
            run->click();
            QTRY_VERIFY_WITH_TIMEOUT(!work->isBusy(), 10000);
            QVERIFY(work->hasResult() && QFile::exists(output->text()));
            QVERIFY(QFile::exists(BondImaging::diagnosticPngPath(output->text())));
            QVERIFY(!HdfFile::inspect(output->text(), "/entry/bond/band_diagnostics/1_test/residual_pyz")
                         .dims.isEmpty());
            work->openVolume(output->text());
        } else {
            work->openVolume(external);
        }
        QTRY_VERIFY_WITH_TIMEOUT(!work->isBusy(), 10000);
        QVERIFY(work->hasResult());
        QCOMPARE(tabs->currentIndex(), BondWorkspace::Volume);
        auto spectrum = work->findChild<SpectrumPlot *>("BondVoxelSpectrum");
        auto slice = work->findChild<ScientificPlot *>("BondVoxelSlice");
        auto figure = work->findChild<QWidget *>("BondEnergyPosition");
        auto style = work->findChild<QComboBox *>("BondSpectralStyle");
        auto source = work->findChild<QComboBox *>("BondViewSource");
        auto channel = work->findChild<QSpinBox *>("BondSpectrumChannel");
        auto plotTabs = work->findChild<QTabWidget *>("BondSpectralPlots");
        QVERIFY(spectrum && slice && figure && style && source && channel && plotTabs);
        QVERIFY(channel->isVisible() && !style->isVisible());
        QVERIFY(slice->data().finiteCount > 0);
        QVERIFY(spectrum->curves().size() == 1);
        if (external.isEmpty()) {
            const auto base = QRectF(slice->rect()).adjusted(82, 57, -130, -73);
            const double scale = std::min(base.width() / 2., base.height() / 2.1);
            const auto pick = base.center() - QPointF(.62 * scale, .8 * scale);
            QTest::mouseClick(slice, Qt::LeftButton, Qt::NoModifier, pick.toPoint());
            QCOMPARE(work->findChild<QSpinBox *>("BondVoxelY")->value(), 0);
            QCOMPARE(work->findChild<QSpinBox *>("BondVoxelX")->value(), 2);
            auto options = work->findChild<QToolButton *>("BondProcessingOptions");
            auto background = work->findChild<QLineEdit *>("BondBackground");
            auto apply = work->findChild<QPushButton *>("ApplyBondSpectralProcessing");
            QVERIFY(options && background && apply && !apply->isVisible());
            options->click();
            if (!dir.isEmpty()) {
                QTest::qWait(60);
                QVERIFY(window.grab().save(QDir(dir).filePath("voxel-processing.png")));
            }
            background->setText("1:2");
            apply->click();
            QTRY_VERIFY_WITH_TIMEOUT(!work->isBusy(), 10000);
            QCOMPARE(source->currentIndex(), 1);
            source->setCurrentIndex(0);
            QCOMPARE(spectrum->curves().size(), 1);
            source->setCurrentIndex(1);
            QCOMPARE(spectrum->curves().size(), 3);
            auto mode = work->findChild<QComboBox *>("BondBackgroundMode");
            auto widePath = work->findChild<QLineEdit *>("BondWideSpectrum");
            auto status = work->findChild<QLabel *>("BondStatus");
            QVERIFY(mode && widePath && status);
            ProcessedSpectrum reference;
            reference.energyEv = {0, 1000};
            reference.intensityUnit = "a.u.";
            reference.crystals = {{"test", {1, 2}}, {"other", {2, 1}}};
            reference.geometry = {{"test", {5, 1, 60}}, {"other", {12, 1, 120}}};
            auto referencePath = temp.filePath("wide.h5");
            SpectrumProcessing::save(referencePath, reference);
            mode->setCurrentIndex(1);
            widePath->setText(referencePath);
            apply->click();
            QTRY_VERIFY_WITH_TIMEOUT(!work->isBusy(), 10000);
            QCOMPARE(source->itemText(1), QStringLiteral("低 q 线性扣除"));
            auto showTransfer = work->findChild<QCheckBox *>("BondShowEnergyTransfer");
            auto analyzer = work->findChild<QDoubleSpinBox *>("BondAnalyzerReference");
            QVERIFY(showTransfer && analyzer);
            analyzer->setValue(1);
            showTransfer->setChecked(true);
            QCOMPARE(source->currentIndex(), 1);
            QCOMPARE(background->text(), QStringLiteral("1:2"));
            showTransfer->setChecked(false);

            const auto original = work->findChild<QLineEdit *>("BondOutput")->text();
            auto mixed = BondImaging::load(original);
            mixed.diagnostics.remove("crystal_labels"); // Exercise the legacy file's ROI keys.
            mixed.diagnostics["crystal_keys"] = QJsonArray{"1_test", "2_other"};
            auto mixedPath = temp.filePath("mixed.h5");
            std::atomic_bool stop{false};
            BondImaging::save(mixedPath, mixed, stop);
            work->openVolume(mixedPath);
            QTRY_VERIFY_WITH_TIMEOUT(!work->isBusy(), 10000);
            apply->click();
            QTRY_VERIFY_WITH_TIMEOUT(!work->isBusy(), 10000);
            QVERIFY(status->text().contains(QStringLiteral("联合了多个晶体")));
            QCOMPARE(source->count(), 1);
            work->openVolume(original);
            QTRY_VERIFY_WITH_TIMEOUT(!work->isBusy(), 10000);
            mode->setCurrentIndex(0);
            source->setCurrentIndex(0);
            options->click();
        }
        if (!dir.isEmpty()) {
            QVERIFY(QDir().mkpath(dir));
            QTest::qWait(60);
            QVERIFY(window.grab().save(QDir(dir).filePath("voxel.png")));
        }
        plotTabs->setCurrentWidget(figure);
        QVERIFY(figure->isVisible() && style->isVisible() && !channel->isVisible());
        for (int representation = 0; representation < style->count(); ++representation) {
            style->setCurrentIndex(representation);
            const QPoint origin = figure->rect().center(), end = origin + QPoint(95, 30);
            const auto initial = figure->grab().toImage();
            QTest::mousePress(figure, Qt::LeftButton, Qt::NoModifier, origin);
            QMouseEvent drag(QEvent::MouseMove, end, figure->mapToGlobal(end), Qt::NoButton, Qt::LeftButton,
                             Qt::NoModifier);
            QApplication::sendEvent(figure, &drag);
            QTest::mouseRelease(figure, Qt::LeftButton, Qt::NoModifier, end);
            const auto rotated = figure->grab().toImage();
            QVERIFY(rotated != initial);
            if (!dir.isEmpty()) {
                const auto filename = representation == 0 ? "energy-position.png" : "waterfall.png";
                QVERIFY(window.grab().save(QDir(dir).filePath(filename)));
            }
            QWheelEvent zoom(origin, figure->mapToGlobal(origin), {}, {0, 120}, Qt::NoButton, Qt::NoModifier,
                             Qt::NoScrollPhase, false);
            QApplication::sendEvent(figure, &zoom);
            QVERIFY(figure->grab().toImage() != rotated);
            QTest::mouseDClick(figure, Qt::LeftButton, Qt::NoModifier, origin);
            QVERIFY(figure->grab().toImage() == initial);
        }
        work->selectStage(BondWorkspace::Results);
        auto mode = work->findChild<QComboBox *>("BondContrastMode");
        auto map = work->findChild<QLabel *>("BondMapInfo");
        auto montage = work->findChild<QWidget *>("BondMontage");
        auto ea = work->findChild<QDoubleSpinBox *>("BondEa");
        QVERIFY(mode && map && montage && ea);
        QVERIFY(channel->isVisible() && !ea->isVisible());
        const int previousChannel = channel->value();
        channel->setValue(0);
        QVERIFY(map->text().contains(QString::number(spectrum->curves().first().x.first(), 'g', 9)));
        channel->setValue(previousChannel);
        mode->setCurrentIndex(1);
        QVERIFY(!channel->isVisible() && ea->isVisible());
        QVERIFY(map->text().contains(QStringLiteral("有效体素")));
        if (!dir.isEmpty())
            QVERIFY(window.grab().save(QDir(dir).filePath("contrast-sections.png")));
        auto surface = work->findChild<QWidget *>("BondIsosurface");
        QVERIFY(surface);
        auto figureTabs = qobject_cast<QTabWidget *>(surface->parentWidget()->parentWidget()->parentWidget());
        QVERIFY(figureTabs);
        figureTabs->setCurrentIndex(1);
        auto build = work->findChild<QPushButton *>("BuildBondSurface");
        QVERIFY(build && build->isVisible());
        build->click();
        QTRY_VERIFY_WITH_TIMEOUT(!work->isBusy(), 10000);
        auto status = work->findChild<QLabel *>("BondStatus");
        QVERIFY(status->text().contains(QStringLiteral("三角形")));
        if (!dir.isEmpty())
            QVERIFY(window.grab().save(QDir(dir).filePath("contrast-isosurface.png")));
        if (external.isEmpty()) {
            work->selectStage(BondWorkspace::Spectra);
            auto prepare = work->findChild<QPushButton *>("PrepareBondProjection");
            auto run = work->findChild<QPushButton *>("RunBondContrast");
            auto contrastEa = work->findChild<QDoubleSpinBox *>("BondContrastEa");
            auto contrastEb = work->findChild<QDoubleSpinBox *>("BondContrastEb");
            auto output = work->findChild<QLineEdit *>("BondOutput");
            auto projection = work->findChild<ScientificPlot *>("BondProjectionPlot");
            auto xStep = work->findChild<QLineEdit *>("BondXStep");
            QVERIFY(prepare && run && contrastEa && contrastEb && output && projection && xStep);
            QVERIFY(run->isVisible() && !run->isEnabled());
            const auto path = temp.filePath("direct-contrast.h5");
            contrastEa->setValue(2);
            contrastEb->setValue(1);
            prepare->click();
            QTRY_VERIFY_WITH_TIMEOUT(!work->isBusy(), 10000);
            QCOMPARE(tabs->currentIndex(), BondWorkspace::Spectra);
            QVERIFY(projection->isVisible() && projection->data().finiteCount == 4);
            QVERIFY(projection->data().maximum < 0.);
            if (!dir.isEmpty()) {
                QVERIFY(window.grab().save(QDir(dir).filePath("detector-spectra-settings.png")));
            }
            work->selectStage(BondWorkspace::Spectra);
            contrastEb->setValue(1.25);
            QVERIFY(!run->isEnabled() && !projection->isVisible());
            QVERIFY(!QFile::exists(path));
            contrastEb->setValue(1);
            prepare->click();
            QTRY_VERIFY_WITH_TIMEOUT(!work->isBusy(), 10000);
            QCOMPARE(tabs->currentIndex(), BondWorkspace::Spectra);
            QVERIFY(run->isVisible() && run->isEnabled() && projection->isVisible());
            output->setText(path);
            xStep->setText("0.125");
            QVERIFY(run->isEnabled() && projection->data().finiteCount == 4);
            run->click();
            QTRY_VERIFY_WITH_TIMEOUT(!work->isBusy(), 10000);
            QVERIFY2(QFile::exists(path), qPrintable(status->text()));
            const auto contrast = BondImaging::loadContrast(path);
            QVERIFY(contrast);
            QCOMPARE(HdfFile::inspect(path, "/entry/bond/contrast").dims, (QVector<quint64>{17, 2, 2}));
            QCOMPARE(contrast->diagnostics.value("ea_kev").toDouble(), 2.);
            QCOMPARE(contrast->diagnostics.value("eb_kev").toDouble(), 1.);
            QVERIFY(contrast->values[32] < 0);
            QCOMPARE(tabs->currentIndex(), BondWorkspace::Results);
            QCOMPARE(mode->currentIndex(), 2);
            QCOMPARE(figureTabs->currentIndex(), 1);
            QVERIFY(map->text().contains(QStringLiteral("S(1) − S(2)")));
            QVERIFY2(status->text().contains(QStringLiteral("三角形")), qPrintable(status->text()));
            if (!dir.isEmpty())
                QVERIFY(window.grab().save(QDir(dir).filePath("direct-contrast.png")));
            work->selectStage(BondWorkspace::Volume);
            plotTabs->setCurrentIndex(1);
            style->setCurrentIndex(1);
            QVERIFY(figure->isVisible() && style->isVisible());
            work->openVolume(path);
            QTRY_VERIFY_WITH_TIMEOUT(!work->isBusy(), 10000);
            QCOMPARE(tabs->currentIndex(), BondWorkspace::Results);
            QCOMPARE(mode->currentIndex(), 2);
            QVERIFY(tabs->isTabEnabled(BondWorkspace::Volume));
        }
    }
    void bondReferenceMaps() {
        QTemporaryDir temp;
        QVERIFY(temp.isValid());
        BondVolume volume;
        constexpr double analyzer = 9.;
        const QVector<double> e{280, 281, 282, 283, 284, 285, 286, 287};
        const QVector<double> matrix{.1, .25, .6, 1, .6, .3, .2, .1};
        const QVector<double> target{.1, .1, .15, .2, .45, 1, .55, .3};
        volume.axes = {{}, {0, 1, 2, 3, 4, 5, 6, 7}, {0, 1, 2, 3, 4, 5}, {0, 1}};
        volume.diagnostics["schema"] = "XRI.bond.v1";
        volume.sensitivity.fill(1., 8);
        volume.converged.fill(true, e.size());
        for (qsizetype i = 0; i < e.size(); ++i) {
            volume.axes[0].append(analyzer + e[i] / 1000.);
            for (int x = 0; x < 8; ++x)
                for (int y = 0; y < 6; ++y)
                    for (int z = 0; z < 2; ++z) {
                        const double coefficient = x == 2 && y > 0 ? -2. : x == 5 && y < 5 ? 2. : 0.;
                        volume.intensity.append(5 + .5 * i + (2 + .1 * y) * matrix[i] + coefficient * target[i]);
                    }
        }
        const auto path = temp.filePath("spectral-volume.h5");
        const std::atomic_bool stop{false};
        BondImaging::save(path, volume, stop);
        const auto matrixPath = temp.filePath("matrix.h5"), targetPath = temp.filePath("target.h5");
        ProcessedSpectrum ref;
        ref.energyEv = e;
        ref.intensityUnit = "a.u.";
        ref.crystals = {{"matrix", matrix}};
        SpectrumProcessing::save(matrixPath, ref);
        ref.crystals = {{"target", target}};
        SpectrumProcessing::save(targetPath, ref);

        MainWindow window;
        window.resize(1200, 840);
        window.show();
        auto work = window.openBondImaging();
        work->openVolume(path);
        QTRY_VERIFY_WITH_TIMEOUT(!work->isBusy(), 10000);
        QVERIFY(work->hasResult());
        work->selectStage(BondWorkspace::Results);
        auto run = work->findChild<QPushButton *>("FitBondReferenceMaps");
        QVERIFY(run);
        auto page = run->parentWidget();
        auto figures = qobject_cast<QTabWidget *>(page->parentWidget()->parentWidget());
        QVERIFY(figures);
        figures->setCurrentWidget(page);
        QVERIFY(run->isVisible());
        work->findChild<QLineEdit *>("BondMatrixReference")->setText(matrixPath);
        work->findChild<QLineEdit *>("BondTargetReference")->setText(targetPath);
        work->findChild<QDoubleSpinBox *>("BondReferenceAnalyzer")->setValue(analyzer);
        work->findChild<QDoubleSpinBox *>("BondReferenceThreshold")->setValue(.5);
        run->click();
        QTRY_VERIFY_WITH_TIMEOUT(!work->isBusy(), 10000);
        const auto residual = work->findChild<ScientificPlot *>("BondReferenceMap1");
        const auto coefficient = work->findChild<ScientificPlot *>("BondReferenceMap2");
        const auto spectra = work->findChild<SpectrumPlot *>("BondReferenceSpectra");
        QVERIFY(residual && coefficient && spectra);
        QCOMPARE(coefficient->data().finiteCount, 48u);
        const auto &map = coefficient->data().values;
        QVERIFY(std::abs(map[10] + 2.) < 1e-10 && std::abs(map[13] - 2.) < 1e-10);
        QCOMPARE(map[0], 0.);
        QVERIFY(residual->data().values[13] > .01 && residual->data().values[0] < 1e-20);
        QCOMPARE(spectra->curves().size(), 4);
        const auto dir = qEnvironmentVariable("XRI_BOND_SCREENSHOTS");
        if (!dir.isEmpty()) {
            QVERIFY(QDir().mkpath(dir));
            QVERIFY(window.grab().save(QDir(dir).filePath("reference-fit-four-panel.png")));
        }
    }
    void academicTicks() {
        for (auto range : {QPair<double, double>{0, 4095},
                           {-0.003, 0.007},
                           {9.680, 9.690},
                           {-3e6, 9e6},
                           {1e-9, 8e-9},
                           {1e8, 1e8 + 10}}) {
            auto ticks = PlotStyle::ticks(range.first, range.second);
            QVERIFY(ticks.major.size() >= 3 && ticks.major.size() <= 12);
            QStringList labels;
            double previous = -INFINITY;
            for (double v : ticks.major) {
                QVERIFY(v > previous && v >= range.first - ticks.step * 1e-8 &&
                        v <= range.second + ticks.step * 1e-8);
                auto text = PlotStyle::label(v, ticks);
                QVERIFY(!labels.contains(text));
                QVERIFY(!text.startsWith("-0.000000"));
                labels.append(text);
                previous = v;
            }
            QVERIFY(!ticks.minor.isEmpty());
        }
        auto log = PlotStyle::ticks(-4095, 4095, 5, true);
        QVERIFY(log.major.contains(0) && log.major.contains(-1000) && log.major.contains(1000));
        QVERIFY(log.major.size() <= 12);
        QVERIFY(PlotStyle::ticks(NAN, 1).major.isEmpty());
        QVERIFY(PlotStyle::ticks(1, 1).major.isEmpty());
    }
    void sharedWorkflowNavigation() {
        QTemporaryDir temp;
        auto file = SpectrumFixture::create(temp.path(), true);
        MainWindow window;
        window.openFile(file);
        QTRY_VERIFY_WITH_TIMEOUT(!window.isReading(), 10000);
        auto stack = window.findChild<QStackedWidget *>("ProcessingWorkspaces");
        QVERIFY(stack);
        auto browser = stack->currentWidget();
        auto spectrum = window.openSpectroscopy(2);
        QVERIFY(!spectrum->isWindow() && stack->currentWidget() == spectrum);
        QVERIFY(window.findChild<QPushButton *>("CloseWorkflow")->isHidden());
        window.findChild<QAction *>("ShowDataBrowser")->trigger();
        QCOMPARE(stack->currentWidget(), browser);
        QCOMPARE(window.openSpectroscopy(), spectrum);
        QCOMPARE(spectrum->findChild<QTabWidget *>("SpectroscopySteps")->currentIndex(), 2);
        auto editor = window.openRoi(file);
        QVERIFY(editor && !editor->isWindow() && stack->currentWidget() == editor);
        QTRY_VERIFY_WITH_TIMEOUT(!editor->isLoading(), 10000);
        QVERIFY(editor->imagesReady());
        QVERIFY(editor->loadFrom(temp.path() + "/test.roi.json"));
        auto saved = editor->rois().items()[0].rect;
        editor->findChild<QPushButton *>("ReturnFromRoi")->click();
        QCOMPARE(stack->currentWidget(), spectrum);
        auto image = window.openImaging();
        image->findChild<QLineEdit *>("input.sample_nxs")->setText(file);
        QCOMPARE(window.openRoi(file), editor);
        QCOMPARE(editor->rois().items()[0].rect, saved);
        editor->findChild<QPushButton *>("ReturnFromRoi")->click();
        QCOMPARE(stack->currentWidget(), image);
        auto roi = image->findChild<QLineEdit *>("input.roi_table");
        QVERIFY(QFileInfo::exists(roi->text()));
        QVERIFY(!image->findChild<QPushButton *>("RunImaging")->isVisible());
        auto bond = window.openBondImaging();
        auto bondSample = bond->findChild<QLineEdit *>("input.sample_nxs");
        auto bondRoi = bond->findChild<QLineEdit *>("input.roi_table");
        auto selectRoi = bond->findChild<QPushButton *>("SelectImagingInputRoi");
        QVERIFY(bondSample && bondRoi && selectRoi);
        bondSample->setText(file);
        selectRoi->click();
        QCOMPARE(stack->currentWidget(), editor);
        editor->findChild<QPushButton *>("ReturnFromRoi")->click();
        QCOMPARE(stack->currentWidget(), bond);
        QCOMPARE(bondRoi->text(), roi->text());
        window.findChild<QAction *>("ShowDataBrowser")->trigger();
        QCOMPARE(stack->currentWidget(), browser);
        QCOMPARE(window.currentFile(), QFileInfo(file).canonicalFilePath());
    }
    void workflowOverview() {
        auto directory = qEnvironmentVariable("XRI_WORKFLOW_SCREENSHOTS");
        if (directory.isEmpty())
            QSKIP("Optional workflow screenshots not requested");
        auto grid = QFileInfo(qEnvironmentVariable("XRI_TEST_FILE")).canonicalFilePath();
        auto scan = QFileInfo(grid).dir().filePath("14516_scan_20260917142117_1.nxs");
        MainWindow window;
        window.openFiles({scan, grid});
        window.resize(1440, 940);
        window.setSliceIndex(0, 1000);
        QTRY_VERIFY_WITH_TIMEOUT(!window.isReading(), 10000);
        QDir().mkpath(directory);
        auto capture = [&](const QString &name) {
            QTest::qWait(60);
            QVERIFY(window.grab().save(directory + "/" + name + ".png"));
        };
        capture("browser");
        auto roi = window.openRoi(scan);
        QTRY_VERIFY_WITH_TIMEOUT(!roi->isLoading(), 15000);
        QVERIFY(roi->imagesReady());
        capture("roi");
        roi->findChild<QCheckBox *>("ShowRoiProfiles")->setChecked(true);
        capture("roi-profiles");
        auto spec = window.openSpectroscopy(0);
        spec->findChild<QLineEdit *>("ElasticFile")->setText(scan);
        if (!qEnvironmentVariable("XRI_REAL_ROI").isEmpty()) {
            QVERIFY(spec->loadRois(qEnvironmentVariable("XRI_REAL_ROI")));
            spec->runStep(0);
            QTRY_VERIFY_WITH_TIMEOUT(!spec->isBusy(), 15000);
            QVERIFY(spec->hasCalibration());
        }
        for (int i = 0; i < 5; ++i) {
            spec->selectStep(i);
            capture(QString("spectrum%1").arg(i));
        }
        auto image = window.openImaging();
        if (!qEnvironmentVariable("XRI_IMAGING_CONFIG").isEmpty())
            QVERIFY(image->importConfiguration(qEnvironmentVariable("XRI_IMAGING_CONFIG")));
        if (!qEnvironmentVariable("XRI_IMAGING_RESULT").isEmpty())
            QVERIFY(image->loadResult(qEnvironmentVariable("XRI_IMAGING_RESULT")));
        const auto imagingSteps = image->findChild<QTabWidget *>("ImagingSteps");
        QVERIFY(imagingSteps);
        for (int i = 0; i < imagingSteps->count(); ++i) {
            image->selectStep(i);
            capture(QString("imaging%1").arg(i));
        }
    }
    void multipleBrowserFiles() {
        auto grid = QFileInfo(qEnvironmentVariable("XRI_TEST_FILE")).canonicalFilePath();
        auto scan = QFileInfo(grid).dir().filePath("14516_scan_20260917142117_1.nxs");
        if (grid.isEmpty() || !QFileInfo::exists(scan))
            QSKIP("Real browser files unavailable");
        MainWindow window;
        window.show();
        QStringList labels;
        for (auto action : window.menuBar()->actions()) {
            auto label = action->text();
            label.remove(QRegularExpression("\\(&[A-Za-z]\\)"));
            labels.append(label);
        }
        QCOMPARE(labels, QStringList({QString::fromUtf8("文件"), QString::fromUtf8("数据浏览器"), "ROI",
                                      QString::fromUtf8("谱学"), QString::fromUtf8("成像"),
                                      QString::fromUtf8("化学键衬度分析"), QString::fromUtf8("语言"),
                                      QString::fromUtf8("帮助")}));
        auto browser = window.findChild<QMenu *>("DataBrowserMenu");
        auto fileMenu = window.findChild<QMenu *>("FileMenu");
        QVERIFY(browser && fileMenu);
        QVERIFY(browser->actions().contains(window.findChild<QAction *>("ResetBrowserPlot")));
        QVERIFY(!fileMenu->actions().contains(window.findChild<QAction *>("ResetBrowserPlot")));
        window.openFile(grid);
        window.setSliceIndex(0, 1000);
        QTRY_VERIFY_WITH_TIMEOUT(!window.isReading(), 10000);
        QCOMPARE(window.currentData().sum, 5274719.);
        window.openFile(scan);
        QTRY_VERIFY_WITH_TIMEOUT(!window.isReading(), 10000);
        QCOMPARE(window.openedFiles().size(), 2);
        QCOMPARE(window.currentData().sum, 37754.);
        auto files = window.findChild<QComboBox *>("BrowserFiles");
        auto tree = window.findChild<QTreeWidget *>("BrowserFileTree");
        QVERIFY(files && tree);
        QCOMPARE(files->count(), 2);
        QCOMPARE(tree->topLevelItemCount(), 2);
        auto child = [](QTreeWidgetItem *parent, const QString &name) {
            for (int i = 0; i < parent->childCount(); ++i)
                if (parent->child(i)->text(0) == name)
                    return parent->child(i);
            return static_cast<QTreeWidgetItem *>(nullptr);
        };
        auto entry = child(tree->topLevelItem(0), "entry");
        QVERIFY(entry);
        tree->topLevelItem(0)->setExpanded(true);
        entry->setExpanded(true);
        auto instrument = child(entry, "instrument");
        QVERIFY(instrument);
        instrument->setExpanded(true);
        auto energy = child(instrument, "M_DCM_B5Energy_setpoint");
        QVERIFY(energy);
        QCOMPARE(energy->text(1), QString("7056"));
        QCOMPARE(window.currentFile(), QFileInfo(scan).canonicalFilePath());
        files->setCurrentIndex(0);
        QTRY_VERIFY_WITH_TIMEOUT(!window.isReading(), 10000);
        QCOMPARE(window.currentFile(), grid);
        QCOMPARE(window.currentData().sum, 5274719.);
        QCOMPARE(window.findChild<QSpinBox *>("FrameIndex")->value(), 1000);
        window.selectPath("/entry/data/D_MINIPIX");
        window.setSliceIndex(0, 100);
        QTRY_VERIFY_WITH_TIMEOUT(!window.isReading(), 10000);
        QCOMPARE(window.currentData().sum, 19113.);
        files->setCurrentIndex(1);
        files->setCurrentIndex(0);
        QTRY_VERIFY_WITH_TIMEOUT(!window.isReading(), 10000);
        QCOMPARE(window.currentData().path, QString("/entry/data/D_MINIPIX"));
        QCOMPARE(window.currentData().sum, 19113.);
        window.openFile(grid);
        QTRY_VERIFY_WITH_TIMEOUT(!window.isReading(), 10000);
        QCOMPARE(window.openedFiles().size(), 2);
        QCOMPARE(window.findChild<QSpinBox *>("FrameIndex")->value(), 0);
        window.closeCurrentFile();
        QTRY_VERIFY_WITH_TIMEOUT(!window.isReading(), 10000);
        QCOMPARE(window.openedFiles().size(), 1);
        QCOMPARE(window.currentFile(), QFileInfo(scan).canonicalFilePath());
        QCOMPARE(window.currentData().sum, 37754.);
        window.setSliceIndex(0, 50);
        window.closeAllFiles();
        QTRY_VERIFY_WITH_TIMEOUT(!window.isReading(), 10000);
        QVERIFY(window.openedFiles().isEmpty() && window.currentFile().isEmpty() &&
                window.currentData().values.isEmpty());
        QVERIFY(!window.findChild<QAction *>("ExportHdfFile")->isEnabled());
    }
    void browserFileExport() {
        QTemporaryDir temp;
        auto config = ImagingFixture::create(temp.path());
        Q_UNUSED(config);
        auto source = temp.path() + "/sample.nxs", other = temp.path() + "/foil.nxs",
             out = temp.path() + "/exported.h5";
        MainWindow window;
        window.openFiles({other, source});
        QTRY_VERIFY_WITH_TIMEOUT(!window.isReading(), 10000);
        QVERIFY(!window.exportFile(source));
        QVERIFY(!window.exportFile(other));
        QVERIFY(!window.exportError().isEmpty());
        QVERIFY(window.exportFile(out));
        window.closeCurrentFile();
        QTRY_VERIFY_WITH_TIMEOUT(!window.isExporting(), 10000);
        QVERIFY2(window.exportError().isEmpty(), qPrintable(window.exportError()));
        auto digest = [](const QString &path) {
            QFile file(path);
            if (!file.open(QIODevice::ReadOnly))
                return QByteArray();
            return QCryptographicHash::hash(file.readAll(), QCryptographicHash::Sha256);
        };
        QVERIFY(!digest(out).isEmpty());
        QCOMPARE(digest(source), digest(out));
        window.openFile(out);
        QTRY_VERIFY_WITH_TIMEOUT(!window.isReading(), 10000);
        QCOMPARE(window.currentData().dims, QVector<quint64>({8, 47, 19}));
        QCOMPARE(window.openedFiles().size(), 2);
        window.findChild<QAction *>("BrowserMode1")->trigger();
        QCOMPARE(window.findChild<QTabWidget *>()->currentIndex(), 1);
    }
    void imagingWorkflow() {
        QTemporaryDir temp;
        auto config = ImagingFixture::create(temp.path(), false, false, false, false, /*spatial=*/true);
        auto output = config["output"].toObject();
        output["plots"] = true;
        config["output"] = output;
        MainWindow window;
        auto action = window.findChild<QAction *>("ImagingStep0");
        QVERIFY(action && action->isEnabled());
        action->trigger();
        auto dialog = window.findChild<ImagingDialog *>();
        QVERIFY(dialog);
        dialog->applyConfiguration(config);
        auto steps = dialog->findChild<QTabWidget *>("ImagingSteps");
        auto avizo = dialog->findChild<QCheckBox *>("output.save_avizo");
        QVERIFY(steps && avizo && !avizo->isChecked());
        QCOMPARE(steps->tabText(3), QStringLiteral("4 剖面分析"));
        QCOMPARE(steps->tabText(4), QStringLiteral("5 三维模型"));
        QCOMPARE(steps->tabText(5), QStringLiteral("6 诊断与日志"));
        auto buildModel = dialog->findChild<QPushButton *>("BuildImagingModel");
        QVERIFY(buildModel && !buildModel->isEnabled());
        avizo->setChecked(true);
        auto order = dialog->findChild<QComboBox *>("scan_loops.sample_order");
        auto points = dialog->findChild<QSpinBox *>("scan_loops.sample_dcm_points");
        QVERIFY(order && points && order->count() == 7);
        order->setCurrentIndex(order->findData(QStringLiteral("Y,Z,DCM")));
        points->setValue(2);
        auto loops = dialog->configuration()["scan_loops"].toObject();
        QCOMPARE(loops["sample_order"].toString(), QStringLiteral("Y,Z,DCM"));
        QCOMPARE(loops["sample_dcm_points"].toInt(), 2);
        QVERIFY(loops["sample_y_points"].isNull());
        auto validate = dialog->findChild<QPushButton *>("ValidateImaging"),
             run = dialog->findChild<QPushButton *>("RunImaging");
        QVERIFY(validate && run);
        QTest::mouseClick(validate, Qt::LeftButton);
        QTRY_VERIFY_WITH_TIMEOUT(!dialog->isBusy(), 10000);
        auto rois = dialog->findChild<QTableWidget *>("ImagingRois");
        QVERIFY(rois && rois->rowCount() == 1);
        QVERIFY(!dialog->hasResult());
        QVERIFY(!QFileInfo::exists(temp.path() + "/output"));
        QTest::mouseClick(run, Qt::LeftButton);
        QTRY_VERIFY_WITH_TIMEOUT(!dialog->isBusy(), 15000);
        QVERIFY2(dialog->hasResult(), qPrintable(dialog->statusText()));
        {
            const auto outputDir = temp.filePath("output");
            ImagingH5 hdf(outputDir + "/reconstruction.h5");
            const auto voxels = hdf.read("reinterp_insect_ck");
            QVERIFY(voxels[40] != voxels[41] && voxels[41] != voxels[42]);
            QFile mesh(outputDir + "/reconstruction.am");
            QVERIFY(mesh.open(QIODevice::ReadOnly));
            QByteArray header;
            while (!mesh.atEnd()) {
                const auto line = mesh.readLine();
                if (line == "@1\n")
                    break;
                header += line;
            }
            QVERIFY(header.startsWith("# AmiraMesh BINARY-LITTLE-ENDIAN 2.1\n"));
            QVERIFY(header.contains("define Lattice 21 2 2\n"));
            QVERIFY(header.contains("BoundingBox -1 1 0 1 -1 0,"));
            QVERIFY(header.contains("Coordinates \"mm\""));
            const auto payload = mesh.readAll();
            QCOMPARE(payload.size(), qsizetype(21 * 2 * 2 * 4));
            QDataStream stream(payload);
            stream.setByteOrder(QDataStream::LittleEndian);
            stream.setFloatingPointPrecision(QDataStream::SinglePrecision);
            for (int z = 0; z < 2; ++z)
                for (int y = 0; y < 2; ++y)
                    for (int x = 0; x < 21; ++x) {
                        float voxel;
                        stream >> voxel;
                        QCOMPARE(double(voxel), voxels[(x * 2 + y) * 2 + z]);
                    }
            QCOMPARE(stream.status(), QDataStream::Ok);
        }
        QImage diagnostic(temp.path() + "/output/diagnostics.png");
        QVERIFY(!diagnostic.isNull() && diagnostic.dotsPerMeterX() >= 11800);
        if (!qEnvironmentVariable("XRI_WORKFLOW_SCREENSHOTS").isEmpty())
            QVERIFY(diagnostic.save(qEnvironmentVariable("XRI_WORKFLOW_SCREENSHOTS") + "/diagnostic.png"));
        auto image = dialog->findChild<ScientificPlot *>("ImagingVolume");
        QVERIFY(image && !image->data().values.isEmpty());
        QVERIFY(dialog->result().summary["tv_converged"].toBool());
        auto slice = dialog->findChild<QSpinBox *>("ImagingSlice");
        QVERIFY(slice);
        slice->setValue(slice->maximum());
        QVERIFY(!image->data().values.isEmpty());
        dialog->selectStep(4);
        auto model = dialog->findChild<QWidget *>("ImagingModel");
        auto modelSource = dialog->findChild<QComboBox *>("ImagingModelSource");
        auto threshold = dialog->findChild<QLineEdit *>("ImagingModelThreshold");
        QVERIFY(model && modelSource && threshold && model->isVisible() && buildModel->isEnabled());
        QCOMPARE(modelSource->count(), 3);
        const auto validThreshold = threshold->text();
        threshold->setText("-1");
        buildModel->click();
        QVERIFY(!dialog->isBusy() && dialog->statusText().contains(QStringLiteral("阈值须介于")));
        threshold->setText(validThreshold);
        buildModel->click();
        QTRY_VERIFY_WITH_TIMEOUT(!dialog->isBusy(), 10000);
        QVERIFY2(dialog->statusText().contains(QStringLiteral("等值面已生成")), qPrintable(dialog->statusText()));
        const QPoint origin = model->rect().center(), end = origin + QPoint(95, 30);
        const auto initial = model->grab().toImage();
        QTest::mousePress(model, Qt::LeftButton, Qt::NoModifier, origin);
        QMouseEvent drag(QEvent::MouseMove, end, model->mapToGlobal(end), Qt::NoButton, Qt::LeftButton,
                         Qt::NoModifier);
        QApplication::sendEvent(model, &drag);
        QTest::mouseRelease(model, Qt::LeftButton, Qt::NoModifier, end);
        QVERIFY(model->grab().toImage() != initial);
        QTest::mouseDClick(model, Qt::LeftButton, Qt::NoModifier, origin);
        QCOMPARE(model->grab().toImage(), initial);
        dialog->hide();
        window.openImaging(3);
        QCOMPARE(window.findChild<ImagingDialog *>(), dialog);
        QVERIFY(dialog->hasResult());
        auto file = dialog->findChild<QLineEdit *>("input.sample_nxs");
        QVERIFY(file);
        file->setText(temp.path() + "/other.nxs");
        QVERIFY(!dialog->hasResult() && image->data().values.isEmpty());
        QVERIFY(dialog->loadResult(temp.path() + "/output"));
        QVERIFY(dialog->hasResult());
        QVERIFY(avizo->isChecked());
        dialog->selectStep(4);
        modelSource->setCurrentIndex(modelSource->findData(1));
        buildModel->click();
        QTRY_VERIFY_WITH_TIMEOUT(!dialog->isBusy(), 10000);
        QVERIFY2(dialog->statusText().contains(QStringLiteral("等值面已生成")), qPrintable(dialog->statusText()));
        dialog->selectStep(3);
        auto open = dialog->findChild<QPushButton *>("OpenImagingHdf");
        QVERIFY(open);
        QTest::mouseClick(open, Qt::LeftButton);
        QTRY_VERIFY_WITH_TIMEOUT(!window.isReading(), 10000);
        QCOMPARE(window.currentData().path, QString("/reinterp_insect_ck"));
        QCOMPARE(window.currentData().dims, QVector<quint64>({21, 2, 2}));
    }
    void imagingRealResult() {
        auto path = qEnvironmentVariable("XRI_IMAGING_RESULT");
        if (path.isEmpty())
            QSKIP("Optional real imaging result not requested");
        ImagingDialog dialog;
        dialog.show();
        auto config = qEnvironmentVariable("XRI_IMAGING_CONFIG");
        if (!config.isEmpty()) {
            QVERIFY(dialog.importConfiguration(config));
            dialog.start(true);
            QTRY_VERIFY_WITH_TIMEOUT(!dialog.isBusy(), 15000);
        }
        QVERIFY2(dialog.loadResult(path), qPrintable(dialog.statusText()));
        auto directory = qEnvironmentVariable("XRI_IMAGING_SCREENSHOTS");
        if (!directory.isEmpty()) {
            QVERIFY(QDir().mkpath(directory));
            auto steps = dialog.findChild<QTabWidget *>("ImagingSteps");
            auto buildModel = dialog.findChild<QPushButton *>("BuildImagingModel");
            QVERIFY(steps && buildModel);
            dialog.selectStep(4);
            buildModel->click();
            QTRY_VERIFY_WITH_TIMEOUT(!dialog.isBusy(), 15000);
            QVERIFY2(dialog.statusText().contains(QStringLiteral("等值面已生成")), qPrintable(dialog.statusText()));
            for (int tab = 0; tab < steps->count(); ++tab) {
                dialog.selectStep(tab);
                QTest::qWait(50);
                QVERIFY(dialog.grab().save(directory + QString("/step%1.png").arg(tab)));
            }
        }
    }
    void spectroscopyRealCalibration() {
        auto roi = qEnvironmentVariable("XRI_REAL_ROI");
        if (roi.isEmpty())
            QSKIP("Optional real ROI calibration not requested");
        auto file = QFileInfo(qEnvironmentVariable("XRI_TEST_FILE"))
                        .dir()
                        .filePath("14516_scan_20260917142117_1.nxs");
        SpectroscopyDialog dialog(file);
        dialog.show();
        QVERIFY(dialog.loadRois(roi));
        dialog.runStep(0);
        QTRY_VERIFY_WITH_TIMEOUT(!dialog.isBusy(), 15000);
        QVERIFY2(dialog.hasCalibration(), qPrintable(dialog.statusText()));
        auto fits = dialog.findChild<QTableWidget *>("ElasticFits");
        QVERIFY(fits);
        fits->setCurrentCell(45, 0);
        auto directory = qEnvironmentVariable("XRI_SPECTROSCOPY_SCREENSHOTS");
        if (!directory.isEmpty()) {
            QVERIFY(QDir().mkpath(directory));
            dialog.selectStep(2);
            QTest::qWait(50);
            QVERIFY(dialog.grab().save(directory + "/real-elastic.png"));
            dialog.selectStep(1);
            QTest::qWait(50);
            QVERIFY(dialog.grab().save(directory + "/real-masks.png"));
        }
    }
    void spectroscopyWorkflow() {
        QTemporaryDir temp;
        auto elastic = SpectrumFixture::create(temp.path(), true);
        auto xrs = SpectrumFixture::create(temp.path(), false);
        auto reverse = SpectrumFixture::create(temp.path(), false, true);
        MainWindow window;
        window.openFile(elastic);
        QTRY_VERIFY_WITH_TIMEOUT(!window.isReading(), 10000);
        auto action = window.findChild<QAction *>("SpectroscopyStep0");
        QVERIFY(action && action->isEnabled());
        action->trigger();
        auto dialog = window.findChild<SpectroscopyDialog *>();
        QVERIFY(dialog);
        QVERIFY(dialog->loadRois(temp.path() + "/test.roi.json"));
        dialog->addScans({xrs, reverse});
        auto calibrate = dialog->findChild<QPushButton *>("CalibrateSpectrum");
        auto extract = dialog->findChild<QPushButton *>("ExtractSpectrum");
        auto combine = dialog->findChild<QPushButton *>("CombineSpectrum");
        auto save = dialog->findChild<QPushButton *>("SaveSpectrum");
        QVERIFY(calibrate && extract && combine && save);
        QVERIFY(!extract->isEnabled() && !combine->isEnabled() && !save->isEnabled());
        dialog->selectStep(1);
        QVERIFY(calibrate->isVisible() && !extract->isVisible());
        QTest::mouseClick(calibrate, Qt::LeftButton);
        QTRY_VERIFY_WITH_TIMEOUT(!dialog->isBusy(), 10000);
        QVERIFY2(dialog->hasCalibration(), qPrintable(dialog->statusText()));
        QVERIFY(extract->isEnabled());
        auto editor = window.openRoi(elastic, temp.path() + "/test.roi.json");
        QTRY_VERIFY_WITH_TIMEOUT(!editor->isLoading(), 10000);
        QVERIFY(editor->imagesReady());
        editor->findChild<QPushButton *>("ReturnFromRoi")->click();
        QVERIFY(dialog->hasCalibration());
        dialog->selectStep(3);
        QVERIFY(extract->isVisible() && !calibrate->isVisible());
        QTest::mouseClick(extract, Qt::LeftButton);
        QTRY_VERIFY_WITH_TIMEOUT(!dialog->isBusy(), 10000);
        QVERIFY2(dialog->hasScans(), qPrintable(dialog->statusText()));
        QTest::mouseClick(combine, Qt::LeftButton);
        QTRY_VERIFY_WITH_TIMEOUT(!dialog->isBusy(), 10000);
        QVERIFY2(dialog->hasResult(), qPrintable(dialog->statusText()));
        QVERIFY(save->isEnabled());
        auto curve = dialog->findChild<SpectrumPlot *>("CombinedSpectrum");
        QVERIFY(curve && !curve->curves().isEmpty());
        QVERIFY(std::abs(dialog->result().total[10] -
                         15 * (10 + .5 * dialog->result().energyTransferEv[10])) < 1e-5);
        auto loadCurrent = dialog->findChild<QPushButton *>("LoadCurrentSpectrumForProcessing");
        auto fitBackground = dialog->findChild<QPushButton *>("RunMomentumBackground");
        auto corrected = dialog->findChild<SpectrumPlot *>("CorrectedCrystalSpectrum");
        auto source = dialog->findChild<QComboBox *>("BackgroundSpectrumSource");
        auto mode = dialog->findChild<QComboBox *>("BackgroundMode");
        auto qMin = dialog->findChild<QDoubleSpinBox *>("HighMomentumMinimum");
        QVERIFY(loadCurrent && fitBackground && corrected && source && mode && qMin);
        dialog->selectStep(6);
        QVERIFY(fitBackground->isEnabled() && !qMin->isVisible());
        dialog->findChild<QLineEdit *>("LowMomentumFitWindows")->setText("0:1.5,8.5:10");
        fitBackground->click();
        QTRY_VERIFY_WITH_TIMEOUT(!corrected->curves().isEmpty(), 10000);
        QCOMPARE(corrected->curves().first().x, dialog->result().energyTransferEv);
        dialog->selectStep(7);
        dialog->findChild<QDoubleSpinBox *>("CorrectedEnergyMinimum")->setValue(2);
        dialog->findChild<QDoubleSpinBox *>("CorrectedEnergyMaximum")->setValue(4);
        QCOMPARE(corrected->curves().first().x.first(), 2.);
        QCOMPARE(corrected->curves().first().x.last(), 4.);
        dialog->selectStep(6);
        mode->setCurrentIndex(1);
        QVERIFY(qMin->isVisible() && !fitBackground->isEnabled() && corrected->curves().isEmpty());
        dialog->selectStep(5);
        QVERIFY(loadCurrent->isVisible());
        loadCurrent->click();
        dialog->selectStep(6);
        source->setCurrentIndex(2);
        qMin->setValue(100);
        QVERIFY(fitBackground->isEnabled());
        fitBackground->click();
        QTRY_VERIFY_WITH_TIMEOUT(!corrected->curves().isEmpty(), 10000);
        if (!qEnvironmentVariable("XRI_SPECTROSCOPY_SCREENSHOTS").isEmpty()) {
            auto directory = qEnvironmentVariable("XRI_SPECTROSCOPY_SCREENSHOTS");
            QVERIFY(QDir().mkpath(directory));
            for (int tab = 0; tab < 8; ++tab) {
                dialog->selectStep(tab);
                QTest::qWait(50);
                QVERIFY(dialog->grab().save(directory + QString("/step%1.png").arg(tab)));
            }
        }
        auto inclusion = dialog->findChild<QTableWidget *>("SpectrumInclusion");
        QVERIFY(inclusion && (inclusion->item(0, 0)->flags() & Qt::ItemIsUserCheckable));
        QVERIFY(!(inclusion->item(2, 0)->flags() & Qt::ItemIsUserCheckable));
        inclusion->item(0, 0)->setCheckState(Qt::Unchecked);
        QVERIFY(!dialog->hasResult() && !save->isEnabled() && dialog->hasScans());
        dialog->runStep(2);
        QTRY_VERIFY_WITH_TIMEOUT(!dialog->isBusy(), 10000);
        QVERIFY(dialog->hasResult() && !dialog->result().included[0]);
        dialog->hide();
        window.openSpectroscopy(4);
        QCOMPARE(window.findChild<SpectroscopyDialog *>(), dialog);
        QVERIFY(dialog->hasResult());
        auto corrections = dialog->findChild<QTableWidget *>("RoiCorrections");
        corrections->item(0, 1)->setText("1");
        QVERIFY(!dialog->hasCalibration() && !dialog->hasScans() && !dialog->hasResult());
        QVERIFY(!save->isEnabled() && !extract->isEnabled());
    }
    void roiLineProfiles() {
        const auto file = QFileInfo(qEnvironmentVariable("XRI_TEST_FILE"))
                              .dir()
                              .filePath("14516_scan_20260917142117_1.nxs");
        if (!QFileInfo::exists(file))
            QSKIP("ROI test file not present");
        RoiDialog dialog(file, QString::fromUtf8("谱学数据处理"), nullptr, false);
        dialog.show();
        QTRY_VERIFY_WITH_TIMEOUT(!dialog.isLoading(), 20000);
        auto show = dialog.findChild<QCheckBox *>("ShowRoiProfiles");
        auto detector = dialog.findChild<QComboBox *>("RoiProfileDetector");
        auto x = dialog.findChild<QSpinBox *>("RoiProfileX"), y = dialog.findChild<QSpinBox *>("RoiProfileY");
        auto horizontal = dialog.findChild<ScientificPlot *>("RoiHorizontalProfile"),
             vertical = dialog.findChild<ScientificPlot *>("RoiVerticalProfile");
        auto lambda = dialog.findChild<ScientificPlot *>("RoiLambdaPlot"),
             mini = dialog.findChild<ScientificPlot *>("RoiMinipixPlot");
        QVERIFY(show && detector && x && y && horizontal && vertical && lambda && mini);
        QTest::mouseClick(show, Qt::LeftButton);
        QVERIFY(horizontal->isVisible());
        auto checkProfiles = [&](ScientificPlot *source) {
            const auto &image = source->data();
            QCOMPARE(horizontal->data().values.size(), image.width);
            QCOMPARE(vertical->data().values.size(), image.height);
            for (int i = 0; i < image.width; ++i)
                QCOMPARE(horizontal->data().values[i], image.values[y->value() * image.width + i]);
            for (int i = 0; i < image.height; ++i)
                QCOMPARE(vertical->data().values[i], image.values[i * image.width + x->value()]);
            QCOMPARE(source->profilePosition(), QPoint(x->value(), y->value()));
        };
        x->setValue(515);
        y->setValue(515);
        checkProfiles(lambda);
        detector->setCurrentIndex(1);
        QCOMPARE(x->maximum(), 255);
        QCOMPARE(y->maximum(), 255);
        QCOMPARE(lambda->profilePosition(), QPoint(-1, -1));
        checkProfiles(mini);
        QSignalSpy picks(lambda, &ScientificPlot::profilePointPicked),
            rois(lambda, &ScientificPlot::rectangleSelected);
        QTest::mouseClick(lambda, Qt::LeftButton, Qt::ControlModifier, lambda->rect().center());
        QCOMPARE(picks.size(), 1);
        QCOMPARE(rois.size(), 0);
        QCOMPARE(dialog.rois().count(RoiState::Valid), 0);
        QCOMPARE(detector->currentIndex(), 0);
        checkProfiles(lambda);
        // Profiles update during the drag, before mouse release, without creating a crystal ROI.
        QPoint dragEnd = lambda->rect().center() + QPoint(60, 40);
        QTest::mousePress(lambda, Qt::LeftButton, Qt::ControlModifier,
                          lambda->rect().center() + QPoint(20, 20));
        int beforeMove = picks.size();
        QTest::mouseMove(lambda, dragEnd);
        QVERIFY(picks.size() > beforeMove);
        checkProfiles(lambda);
        QTest::mouseRelease(lambda, Qt::LeftButton, Qt::NoModifier, dragEnd);
        const int fixedY = y->value(), oldX = x->value();
        QPoint from = dragEnd - QPoint(0, 30), to = from + QPoint(20, -15);
        QTest::mousePress(lambda, Qt::LeftButton, Qt::NoModifier, from);
        QTest::mouseMove(lambda, to);
        QCOMPARE(y->value(), fixedY);
        QVERIFY(x->value() != oldX);
        checkProfiles(lambda);
        QTest::mouseRelease(lambda, Qt::LeftButton, Qt::NoModifier, to);
        const int fixedX = x->value();
        from = dragEnd - QPoint(50, 0);
        to = from + QPoint(15, -20);
        QTest::mousePress(lambda, Qt::LeftButton, Qt::NoModifier, from);
        QTest::mouseMove(lambda, to);
        QCOMPARE(x->value(), fixedX);
        QVERIFY(y->value() != fixedY);
        checkProfiles(lambda);
        QTest::mouseRelease(lambda, Qt::LeftButton, Qt::NoModifier, to);
        QTest::mousePress(lambda, Qt::LeftButton, Qt::ControlModifier, lambda->rect().center());
        QTest::mouseMove(lambda, QPoint(-50, -50));
        QCOMPARE(x->value(), 0);
        QCOMPARE(y->value(), 0);
        checkProfiles(lambda);
        QTest::keyClick(lambda, Qt::Key_Escape);
        QTest::mouseMove(lambda, lambda->rect().center());
        QCOMPARE(x->value(), 0);
        QCOMPARE(y->value(), 0);
        QTest::mouseRelease(lambda, Qt::LeftButton, Qt::NoModifier, lambda->rect().center());
        QCOMPARE(rois.size(), 0);
        QCOMPARE(dialog.rois().count(RoiState::Valid), 0);
        auto mode = dialog.findChild<QComboBox *>("RoiImageMode");
        mode->setCurrentIndex(1);
        QVERIFY(!x->isEnabled());
        QTRY_VERIFY_WITH_TIMEOUT(!dialog.isLoading(), 10000);
        checkProfiles(lambda);
        auto frame = dialog.findChild<QSpinBox *>("RoiFrame");
        frame->setValue(50);
        QTRY_VERIFY_WITH_TIMEOUT(!dialog.isLoading(), 10000);
        checkProfiles(lambda);
        auto log = dialog.findChild<QCheckBox *>("RoiProfileLog");
        auto counts = horizontal->data().values;
        QTest::mouseClick(log, Qt::LeftButton);
        QCOMPARE(horizontal->data().values, counts);
        QPoint center = lambda->rect().center();
        QTest::mousePress(lambda, Qt::LeftButton, Qt::NoModifier, center - QPoint(15, 15));
        QTest::mouseMove(lambda, center + QPoint(15, 15));
        QTest::mouseRelease(lambda, Qt::LeftButton, Qt::NoModifier, center + QPoint(15, 15));
        QCOMPARE(dialog.rois().items()[0].state, RoiState::Valid);
        auto table = dialog.findChild<QTableWidget *>("CrystalRoiTable");
        table->setCurrentCell(0, 0, QItemSelectionModel::ClearAndSelect | QItemSelectionModel::Rows);
        auto centerButton = dialog.findChild<QPushButton *>("CenterRoiProfiles");
        QVERIFY(centerButton->isEnabled());
        QTest::mouseClick(centerButton, Qt::LeftButton);
        QCOMPARE(QPoint(x->value(), y->value()), dialog.rois().items()[0].rect.center());
        checkProfiles(lambda);
        QTest::mouseClick(show, Qt::LeftButton);
        QVERIFY(!horizontal->isVisible());
        QVERIFY(horizontal->data().values.isEmpty());
        QCOMPARE(lambda->profilePosition(), QPoint(-1, -1));
        QCOMPARE(dialog.rois().items()[0].state, RoiState::Valid);
    }
    void sharedCrystalRoiSelection() {
        const auto grid = qEnvironmentVariable("XRI_TEST_FILE");
        const auto file = QFileInfo(grid).dir().filePath("14516_scan_20260917142117_1.nxs");
        if (!QFileInfo::exists(file))
            QSKIP("ROI test file not present");
        RoiDialog dialog(file, QString::fromUtf8("谱学数据处理"), nullptr, false);
        dialog.show();
        QTRY_VERIFY_WITH_TIMEOUT(!dialog.isLoading(), 20000);
        QVERIFY(dialog.imagesReady());
        auto lambda = dialog.findChild<ScientificPlot *>("RoiLambdaPlot");
        auto mini = dialog.findChild<ScientificPlot *>("RoiMinipixPlot");
        QCOMPARE(lambda->data().sum, 54786802.0);
        QCOMPARE(mini->data().sum, 9.0);
        auto table = dialog.findChild<QTableWidget *>("CrystalRoiTable");
        auto invalid = dialog.findChild<QPushButton *>("InvalidateCrystalRoi");
        auto clear = dialog.findChild<QPushButton *>("ClearCrystalRoi");
        auto finish = dialog.findChild<QPushButton *>("FinishCrystalRois");
        QVERIFY(table && invalid && clear && finish);
        QCOMPARE(table->rowCount(), 75);
        QVERIFY(!finish->isEnabled());
        QSignalSpy selected(lambda, &ScientificPlot::rectangleSelected);
        auto drag = [](ScientificPlot *plot) {
            QPoint center = plot->rect().center();
            QTest::mousePress(plot, Qt::LeftButton, Qt::NoModifier, center + QPoint(25, 25));
            QTest::mouseMove(plot, center - QPoint(25, 25));
            QTest::mouseRelease(plot, Qt::LeftButton, Qt::NoModifier, center - QPoint(25, 25));
        };
        drag(lambda);
        QCOMPARE(selected.size(), 1);
        QCOMPARE(dialog.rois().items()[0].state, RoiState::Valid);
        QCOMPARE(dialog.rois().items()[0].rect, selected[0][0].toRect());
        QCOMPARE(table->currentRow(), 1);
        QTest::mouseClick(invalid, Qt::LeftButton);
        QCOMPARE(dialog.rois().items()[1].state, RoiState::Invalid);
        QVERIFY(dialog.rois().items()[1].rect.isEmpty());
        table->setCurrentCell(60, 0);
        drag(lambda);
        QCOMPARE(dialog.rois().items()[60].state, RoiState::Pending);
        drag(mini);
        QCOMPARE(dialog.rois().items()[60].state, RoiState::Valid);
        table->setCurrentCell(60, 0);
        QTest::mouseClick(clear, Qt::LeftButton);
        QCOMPARE(dialog.rois().items()[60].state, RoiState::Pending);
        QTemporaryDir temp;
        auto path = temp.filePath("crystals.roi.json");
        QVERIFY(dialog.saveTo(path));
        table->setCurrentCell(0, 0);
        QTest::mouseClick(clear, Qt::LeftButton);
        QVERIFY(dialog.loadFrom(path));
        QCOMPARE(dialog.rois().items()[0].state, RoiState::Valid);
        QCOMPARE(dialog.rois().items()[1].state, RoiState::Invalid);
        // Ctrl-click disjoint crystals, including an existing valid ROI and both detectors.
        const auto existingRect = dialog.rois().items()[0].rect;
        table->setCurrentCell(0, 0, QItemSelectionModel::ClearAndSelect | QItemSelectionModel::Rows);
        for (int row : {3, 5})
            QTest::mouseClick(table->viewport(), Qt::LeftButton, Qt::ControlModifier,
                              table->visualRect(table->model()->index(row, 0)).center());
        table->selectionModel()->select(table->model()->index(60, 0),
                                        QItemSelectionModel::Select | QItemSelectionModel::Rows);
        QCOMPARE(table->selectionModel()->selectedRows().size(), 4);
        QTest::mouseClick(invalid, Qt::LeftButton);
        QCOMPARE(table->selectionModel()->selectedRows().size(), 4);
        for (int row : {0, 3, 5, 60})
            QCOMPARE(dialog.rois().items()[row].state, RoiState::Invalid);
        QCOMPARE(dialog.rois().items()[0].rect, existingRect);
        QCOMPARE(dialog.rois().items()[2].state, RoiState::Pending);
        table->scrollToItem(table->item(6, 0));
        QTest::mouseClick(table->viewport(), Qt::LeftButton, Qt::NoModifier,
                          table->visualRect(table->model()->index(6, 0)).center());
        table->scrollToItem(table->item(9, 0));
        QTest::mouseClick(table->viewport(), Qt::LeftButton, Qt::ShiftModifier,
                          table->visualRect(table->model()->index(9, 0)).center());
        QCOMPARE(table->selectionModel()->selectedRows().size(), 4);
        QTest::mouseClick(invalid, Qt::LeftButton);
        for (int row = 6; row <= 9; ++row)
            QCOMPARE(dialog.rois().items()[row].state, RoiState::Invalid);
        QVERIFY(dialog.saveTo(path));
        table->setCurrentCell(3, 0, QItemSelectionModel::ClearAndSelect | QItemSelectionModel::Rows);
        QTest::mouseClick(clear, Qt::LeftButton);
        QVERIFY(dialog.loadFrom(path));
        QCOMPARE(dialog.rois().items()[3].state, RoiState::Invalid);
        table->setFocus();
        QTest::keyClick(table, Qt::Key_A, Qt::ControlModifier);
        QCOMPARE(table->selectionModel()->selectedRows().size(), 75);
        QTest::mouseClick(invalid, Qt::LeftButton);
        QCOMPARE(dialog.rois().count(RoiState::Invalid), 75);
        QCOMPARE(dialog.rois().items()[0].rect, existingRect);
        QVERIFY(finish->isEnabled());
        QCOMPARE(dialog.rois().count(RoiState::Valid) + dialog.rois().count(RoiState::Invalid), 75);
        QTest::mouseClick(finish, Qt::LeftButton);
        QCOMPARE(dialog.result(), int(QDialog::Accepted));
    }
    void realFileInteraction() {
        const QString file = qEnvironmentVariable("XRI_TEST_FILE");
        if (!QFileInfo::exists(file))
            QSKIP("Real NeXus test file is not present");
        MainWindow window;
        window.show();
        window.openFile(file);
        QTRY_VERIFY_WITH_TIMEOUT(!window.isReading(), 10000);
        QCOMPARE(window.currentData().width, 516);
        QVERIFY(window.windowTitle().contains("grid_scan"));
        QCOMPARE(window.currentData().maximum, 0.0);
        auto frame = window.findChild<QSpinBox *>("FrameIndex");
        QVERIFY(frame && frame->isEnabled());
        frame->setValue(1000);
        QTRY_VERIFY_WITH_TIMEOUT(!window.isReading(), 10000);
        QCOMPARE(window.currentData().sum, 5274719.0);

        // Keep a selected table visible when the frame changes.
        auto tabs = window.findChild<QTabWidget *>();
        QVERIFY(tabs);
        tabs->setCurrentIndex(1);
        window.setSliceIndex(0, 100);
        QTRY_VERIFY_WITH_TIMEOUT(!window.isReading(), 10000);
        QCOMPARE(tabs->currentIndex(), 1);
        QCOMPARE(window.currentData().sum, 2415629.0);

        // Multiple queued changes must resolve to the last selected dataset.
        window.setSliceIndex(0, 7055);
        window.selectPath("/entry/data/D_MINIPIX");
        window.setSliceIndex(0, 100);
        QTRY_VERIFY_WITH_TIMEOUT(!window.isReading(), 10000);
        QCOMPARE(window.currentData().path, QString("/entry/data/D_MINIPIX"));
        QCOMPARE(window.currentData().sum, 19113.0);

        auto plot = window.findChild<ScientificPlot *>("ScientificPlot");
        QVERIFY(plot);
        plot->setLogarithmic(true);
        plot->setPalette(1);
        window.selectPath("/entry/data/D_LAMBDA_image");
        window.setSliceIndex(0, 7055);
        QTRY_VERIFY_WITH_TIMEOUT(!window.isReading(), 10000);
        QVERIFY(window.currentData().error.isEmpty());
        QVERIFY(window.currentData().text.contains("1 element(s)"));
        QCOMPARE(tabs->currentIndex(), 2);

        const QString scanFile = QFileInfo(file).dir().filePath("14516_scan_20260917142117_1.nxs");
        if (QFileInfo::exists(scanFile)) {
            window.openFile(scanFile);
            QTRY_VERIFY_WITH_TIMEOUT(!window.isReading(), 10000);
            QVERIFY(window.windowTitle().contains(QString::fromUtf8("scan · 一维扫描")));
            QCOMPARE(window.currentData().dims, QVector<quint64>({51, 516, 516}));
            QCOMPARE(frame->maximum(), 50);
            QCOMPARE(frame->value(), 0);
            QCOMPARE(window.currentData().sum, 37754.0);
            frame->setValue(50);
            QTRY_VERIFY_WITH_TIMEOUT(!window.isReading(), 10000);
            QCOMPARE(window.currentData().sum, 65206.0);
            window.selectPath("/entry/instrument/M_DCM_B5Energy_readback");
            QTRY_VERIFY_WITH_TIMEOUT(!window.isReading(), 10000);
            QCOMPARE(window.currentData().width, 51);
            QVERIFY(window.currentData().error.isEmpty());
            window.openFile(file);
            QTRY_VERIFY_WITH_TIMEOUT(!window.isReading(), 10000);
            QCOMPARE(frame->maximum(), 7055);
            QCOMPARE(frame->value(), 0);
            QVERIFY(window.windowTitle().contains(QString::fromUtf8("grid_scan · 网格扫描")));
        }
    }
};

int main(int argc, char **argv) {
    QApplication app(argc, argv);
#ifdef Q_OS_WIN
    if (!qEnvironmentVariable("XRI_SPECTROSCOPY_SCREENSHOTS").isEmpty() ||
        !qEnvironmentVariable("XRI_IMAGING_SCREENSHOTS").isEmpty() ||
        !qEnvironmentVariable("XRI_WORKFLOW_SCREENSHOTS").isEmpty() ||
        !qEnvironmentVariable("XRI_ROI_ML_SCREENSHOT").isEmpty() ||
        !qEnvironmentVariable("XRI_BOND_SCREENSHOTS").isEmpty()) {
        QFontDatabase::addApplicationFont("C:/Windows/Fonts/arial.ttf");
        QFontDatabase::addApplicationFont("C:/Windows/Fonts/msyh.ttc");
        app.setFont(QFont("Microsoft YaHei", 9));
    }
#endif
    QCoreApplication::setOrganizationName("XRI-Tests");
    QCoreApplication::setApplicationName("XRI-Tests");
    GuiTests tests;
    return QTest::qExec(&tests, argc, argv);
}
#include "GuiTests.moc"
