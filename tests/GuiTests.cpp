#include "ImagingDialog.h"
#include "ImagingFixture.h"
#include "MainWindow.h"
#include "PlotStyle.h"
#include "RoiDialog.h"
#include "ScientificPlot.h"
#include "SpectroscopyDialog.h"
#include "SpectrumFixture.h"
#include "SpectrumPlot.h"
#include <QAction>
#include <QApplication>
#include <QCheckBox>
#include <QComboBox>
#include <QCryptographicHash>
#include <QDir>
#include <QElapsedTimer>
#include <QFileInfo>
#include <QFontDatabase>
#include <QLineEdit>
#include <QMenu>
#include <QMenuBar>
#include <QPushButton>
#include <QSignalSpy>
#include <QSpinBox>
#include <QStackedWidget>
#include <QTabWidget>
#include <QTableWidget>
#include <QTemporaryDir>
#include <QTreeWidget>
#include <QtTest>
#include <cmath>

class GuiTests : public QObject {
    Q_OBJECT
  private slots:
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
        for (int i = 0; i < 5; ++i) {
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
        auto config = ImagingFixture::create(temp.path());
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
            for (int tab = 0; tab < 5; ++tab) {
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
        if (!qEnvironmentVariable("XRI_SPECTROSCOPY_SCREENSHOTS").isEmpty()) {
            auto directory = qEnvironmentVariable("XRI_SPECTROSCOPY_SCREENSHOTS");
            QVERIFY(QDir().mkpath(directory));
            for (int tab = 0; tab < 5; ++tab) {
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
        table->setCurrentCell(6, 0, QItemSelectionModel::ClearAndSelect | QItemSelectionModel::Rows);
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
    if (!qEnvironmentVariable("XRI_SPECTROSCOPY_SCREENSHOTS").isEmpty() ||
        !qEnvironmentVariable("XRI_IMAGING_SCREENSHOTS").isEmpty() ||
        !qEnvironmentVariable("XRI_WORKFLOW_SCREENSHOTS").isEmpty()) {
        QFontDatabase::addApplicationFont("C:/Windows/Fonts/arial.ttf");
        QFontDatabase::addApplicationFont("C:/Windows/Fonts/msyh.ttc");
        app.setFont(QFont("Microsoft YaHei", 9));
    }
    QCoreApplication::setOrganizationName("XRI-Tests");
    QCoreApplication::setApplicationName("XRI-Tests");
    GuiTests tests;
    return QTest::qExec(&tests, argc, argv);
}
#include "GuiTests.moc"
