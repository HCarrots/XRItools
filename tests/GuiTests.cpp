#include "BondWorkspace.h"
#include "ImagingDialog.h"
#include "ImagingFixture.h"
#include "MainWindow.h"
#include "PlotStyle.h"
#include "RoiDetectionDialog.h"
#include "RoiDialog.h"
#include "ScientificPlot.h"
#include "SpectroscopyDialog.h"
#include "SpectralMaskDialog.h"
#include "SpectrumFixture.h"
#include "SpectrumPlot.h"
#include "SpectrumProcessing.h"
#include "WorkspacePaths.h"
#include <QAction>
#include <QApplication>
#include <QCheckBox>
#include <QComboBox>
#include <QContextMenuEvent>
#include <QCryptographicHash>
#include <QDataStream>
#include <QDialogButtonBox>
#include <QDir>
#include <QDoubleSpinBox>
#include <QElapsedTimer>
#include <QFileDialog>
#include <QFileInfo>
#include <QFontDatabase>
#include <QGraphicsView>
#include <QGroupBox>
#include <QImageReader>
#include <QInputDialog>
#include <QJsonArray>
#include <QJsonDocument>
#include <QJsonObject>
#include <QLabel>
#include <QLineEdit>
#include <QListWidget>
#include <QMenu>
#include <QMenuBar>
#include <QMouseEvent>
#include <QPlainTextEdit>
#include <QPointer>
#include <QProgressBar>
#include <QPushButton>
#include <QSettings>
#include <QScopeGuard>
#include <QSignalSpy>
#include <QSpinBox>
#include <QSplitter>
#include <QStackedWidget>
#include <QStatusBar>
#include <QTabBar>
#include <QTabWidget>
#include <QTableView>
#include <QTableWidget>
#include <QTemporaryDir>
#include <QTimer>
#include <QToolButton>
#include <QTranslator>
#include <QTreeWidget>
#include <QWheelEvent>
#include <QtTest>
#include <cmath>
#include <type_traits>

class GuiTests : public QObject {
    Q_OBJECT

    template <typename Control>
    static void chooseWorkflowFile(Control *control, const QString &path, double seconds = 0.) {
        const bool native = !QApplication::testAttribute(Qt::AA_DontUseNativeDialogs);
        QApplication::setAttribute(Qt::AA_DontUseNativeDialogs);
        QElapsedTimer deadline;
        deadline.start();
        QPointer<QFileDialog> activePicker;
        QTimer driveDialogs;
        connect(&driveDialogs, &QTimer::timeout, control, [&] {
            auto *modal = qobject_cast<QDialog *>(QApplication::activeModalWidget());
            if (!modal)
                return;
            if (deadline.hasExpired(5000)) {
                modal->reject();
            } else if (auto *period = qobject_cast<QInputDialog *>(modal)) {
                period->setDoubleValue(seconds);
                period->accept();
            } else if (auto *picker = qobject_cast<QFileDialog *>(modal)) {
                if (activePicker != picker) {
                    activePicker = picker;
                    const QFileInfo file(path);
                    picker->setDirectory(file.absolutePath());
                    picker->selectFile(file.fileName());
                    if (auto *name = picker->findChild<QLineEdit *>("fileNameEdit"))
                        name->setText(file.fileName());
                    return;
                }
                QMetaObject::invokeMethod(picker, "accept");
            } else {
                modal->reject();
            }
        });
        driveDialogs.start(20);
        if constexpr (std::is_same_v<Control, QAction>)
            control->trigger();
        else
            control->click();
        QApplication::setAttribute(Qt::AA_DontUseNativeDialogs, !native);
    }

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
            auto *roiMenu = english.findChild<QMenu *>("RoiMenu");
            QVERIFY(roiMenu);
            QCOMPARE(roiMenu->title(), QStringLiteral("&ROI"));
            QCOMPARE(english.findChild<QAction *>("SelectCrystalRois")->text(), QStringLiteral("Manual ROI"));
            QCOMPARE(english.findChild<QAction *>("DetectSpotRois")->text(), QStringLiteral("Automatic ROI"));
            QVERIFY(english.findChild<QAction *>("Language_en")->isChecked());
            auto *spectroscopy = english.openSpectroscopy()->findChild<QTabWidget *>("SpectroscopySteps");
            auto *imaging = english.openImaging()->findChild<QTabWidget *>("ImagingSteps");
            auto *bond = english.openBondImaging()->findChild<QTabWidget *>("BondSteps");
            QVERIFY(spectroscopy && imaging && bond);
            for (auto tabs : {spectroscopy, imaging, bond}) {
                const auto tools = tabs->cornerWidget()->findChildren<QPushButton *>();
                QVERIFY(tools.size() >= 2);
                QCOMPARE(tools[0]->text(), QStringLiteral("Load existing configuration…"));
                QCOMPARE(tools[1]->text(), QStringLiteral("Save current configuration…"));
            }
            QCOMPARE(spectroscopy->tabText(1), QStringLiteral("2 Masks and Spectral Discrimination"));
            QCOMPARE(spectroscopy->tabText(3), QStringLiteral("4 XRS and I0"));
            QCOMPARE(imaging->tabText(2), QStringLiteral("3 Input and Calibration Checks"));
            QCOMPARE(bond->tabText(1), QStringLiteral("2 ROI Calibration"));
            QCOMPARE(bond->tabText(2), QStringLiteral("3 Crystal Selection"));
            QCOMPARE(bond->count(), 6);
            QCOMPARE(bond->tabText(3), QStringLiteral("4 4D Voxel Reconstruction"));
            QCOMPARE(bond->tabText(4), QStringLiteral("5 4D Results and Voxel Spectra"));
            QCOMPARE(bond->tabText(5), QStringLiteral("6 3D Result Analysis"));
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
        auto *match = window.findChild<QPushButton *>("MatchCurrentSpotRois");
        QVERIFY(match && match->isEnabled());
        match->click();
        auto *back = window.findChild<QPushButton *>("ReturnToSpotRoiResults");
        QVERIFY(back && back->isVisible() && !QApplication::activeModalWidget());
        back->click();
        QVERIFY(results->isVisible() && results->rowCount() > 0);
        QTemporaryDir other;
        window.openFile(SpectrumFixture::create(other.path(), true));
        QTRY_VERIFY_WITH_TIMEOUT(!window.isReading(), 10000);
        action->trigger();
        QCOMPARE(results->rowCount(), 0);
        QVERIFY(run->isEnabled() && !exportButton->isEnabled() && !match->isEnabled());
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
                                 spot(2, 5, 5, QJsonArray{QJsonValue(QJsonArray{45, 2})}),
                                 spot(3, 6, 5, QJsonArray{QJsonValue(QJsonArray{46, 1})})}}}).toJson();
        QFile detection(temp.filePath("detected.json"));
        QVERIFY(detection.open(QIODevice::WriteOnly));
        QCOMPARE(detection.write(report), report.size());
        detection.close();
        const auto output = temp.filePath("shared.roi.json");
        RoiDetectionDialog dialog(source);
        dialog.show();
        auto action = dialog.findChild<QPushButton *>("MatchSavedSpotRois");
        QVERIFY(action && action->isEnabled());
        auto pickFile = [&](QPushButton *button, const QString &path) {
            const bool nonNative = QApplication::testAttribute(Qt::AA_DontUseNativeDialogs);
            QApplication::setAttribute(Qt::AA_DontUseNativeDialogs);
            QTimer driveDialog;
            QElapsedTimer deadline;
            deadline.start();
            bool selected = false;
            connect(&driveDialog, &QTimer::timeout, &dialog, [&] {
                auto *modal = qobject_cast<QDialog *>(QApplication::activeModalWidget());
                if (!modal)
                    return;
                auto *picker = qobject_cast<QFileDialog *>(modal);
                if (!picker || deadline.hasExpired(5000)) {
                    modal->reject();
                } else if (selected) {
                    QMetaObject::invokeMethod(picker, "accept");
                } else {
                    const QFileInfo file(path);
                    picker->setDirectory(file.absolutePath());
                    picker->selectFile(file.fileName());
                    if (auto *name = picker->findChild<QLineEdit *>("fileNameEdit"))
                        name->setText(file.fileName());
                    selected = true;
                }
            });
            driveDialog.start(20);
            button->click();
            QApplication::setAttribute(Qt::AA_DontUseNativeDialogs, nonNative);
            return selected;
        };
        QVERIFY(pickFile(action, detection.fileName()));
        auto *table = dialog.findChild<QTableWidget *>("SpotRoiCrystalMatches");
        auto *save = dialog.findChild<QPushButton *>("SaveMatchedSpotRois");
        auto *load = dialog.findChild<QPushButton *>("LoadSpotRoiBaseline");
        auto *view = dialog.findChild<QGraphicsView *>("SpotRoiMaskView");
        auto *candidates = dialog.findChild<QComboBox *>("OverlappingSpotCandidates");
        QVERIFY(table && save && load && view && candidates);
        QVERIFY(table->isVisible() && !QApplication::activeModalWidget());
        QVERIFY(table->item(0, 2)->text().isEmpty() && !save->isEnabled());
        QVERIFY(pickFile(load, temp.filePath("test.roi.json")));
        QVERIFY(!save->isEnabled() && table->item(1, 2)->text().isEmpty());
        table->selectRow(0);
        QTest::mouseClick(view->viewport(), Qt::LeftButton, Qt::NoModifier, view->mapFromScene(6.5, 5.5));
        QVERIFY(candidates->isVisible());
        QCOMPARE(candidates->count(), 2);
        const int candidateIndex = candidates->findText(QStringLiteral("光斑 3"));
        QVERIFY(candidateIndex >= 0);
        candidates->setCurrentIndex(candidateIndex);
        QVERIFY(QMetaObject::invokeMethod(candidates, "activated", Q_ARG(int, candidateIndex)));
        QVERIFY(!candidates->isVisible());
        table->selectRow(1);
        QTest::mouseClick(view->viewport(), Qt::LeftButton, Qt::NoModifier, view->mapFromScene(1.5, 1.5));
        QCOMPARE(table->item(0, 2)->text(), QStringLiteral("3"));
        QCOMPARE(table->item(1, 2)->text(), QStringLiteral("1"));
        QVERIFY(save->isEnabled() && pickFile(save, output));
        QVERIFY(action->isVisible());
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
    void workflowConfigurationTools() {
        QTemporaryDir temp;
        const auto elastic = SpectrumFixture::create(temp.path(), true);
        const auto scan = SpectrumFixture::create(temp.path(), false);
        const auto reverse = SpectrumFixture::create(temp.path(), false, true);
        const auto roiPath = temp.filePath("test.roi.json");
        const auto configPath = temp.filePath("spectroscopy.config.json");
        SpectroscopyDialog spectrum(elastic);
        spectrum.show();
        QVERIFY(spectrum.loadRois(roiPath));
        auto groups = spectrum.findChild<QTableWidget *>("SpectroscopyEdgeGroups");
        auto nextStep = spectrum.findChild<QPushButton *>("NextSpectrumStep");
        auto newGroup = spectrum.findChild<QPushButton *>("AddSpectroscopyEdgeGroup");
        auto tabs = spectrum.findChild<QTabWidget *>("SpectroscopySteps");
        QVERIFY(tabs);
        auto measureCompton = tabs->widget(0)->findChild<QCheckBox *>("MeasureCompton");
        auto advanced = tabs->widget(0)->findChild<QWidget *>("SpectroscopyAdvancedInputs");
        QVERIFY(advanced);
        auto advancedToggle = advanced->findChild<QToolButton *>();
        QVERIFY(advancedToggle && !advancedToggle->isChecked());
        auto enableHf = tabs->widget(6)->findChild<QCheckBox *>("EnableBackground_hf_core");
        auto hfTable = tabs->widget(0)->findChild<QLineEdit *>("BackgroundHfTable");
        auto composition = tabs->widget(0)->findChild<QLineEdit *>("BackgroundHfComposition");
        auto hfShells = tabs->widget(0)->findChild<QLineEdit *>("BackgroundHfShells");
        auto chooseHfShells = tabs->widget(0)->findChild<QPushButton *>("ChooseSpectrumHfShells");
        auto restoreHfTable = tabs->widget(0)->findChild<QPushButton *>("RestoreDefaultHfTable");
        QVERIFY(hfTable && composition && hfShells && chooseHfShells && restoreHfTable && measureCompton && enableHf);
        QVERIFY(!tabs->widget(6)->findChild<QLineEdit *>("BackgroundHfShells"));
        QVERIFY(!tabs->widget(6)->findChild<QLineEdit *>("BackgroundHfTable"));
        QVERIFY(!tabs->widget(6)->findChild<QLineEdit *>("BackgroundHfComposition"));
        QVERIFY(hfTable->text().isEmpty() && !chooseHfShells->isEnabled());
        QVERIFY(groups && nextStep && newGroup && !nextStep->isEnabled());
        auto chooseGroup = [&](const QString &symbol, const QString &edge, bool duplicate = false) {
            bool ok = false;
            QTimer::singleShot(0, &spectrum, [&] {
                auto picker = qobject_cast<QDialog *>(QApplication::activeModalWidget());
                if (!picker)
                    return;
                auto atom = picker->findChild<QComboBox *>("NewEdgeGroupElement");
                auto shell = picker->findChild<QComboBox *>("NewEdgeGroupShell");
                auto name = picker->findChild<QLabel *>("NewEdgeGroupName");
                auto buttons = picker->findChild<QDialogButtonBox *>();
                if (atom && shell && name && buttons) {
                    auto create = buttons->button(QDialogButtonBox::Ok);
                    const bool initiallyDisabled = !create->isEnabled();
                    atom->setCurrentText(symbol);
                    const bool requiresShell = !create->isEnabled();
                    shell->setCurrentText(edge);
                    const auto expected = symbol + '_' + edge;
                    ok = initiallyDisabled && requiresShell &&
                         (duplicate ? !create->isEnabled() && name->text().startsWith(expected)
                                    : create->isEnabled() && name->text() == expected);
                    if (ok && !duplicate) {
                        const auto images = qEnvironmentVariable("XRI_SPECTROSCOPY_SCREENSHOTS");
                        if (!images.isEmpty() && expected == "Ni_L3")
                            ok = picker->grab().save(images + "/new-edge-group.png");
                        create->click();
                        return;
                    }
                }
                picker->reject();
            });
            newGroup->click();
            return ok;
        };
        spectrum.runStep(0);
        QVERIFY(!spectrum.isBusy() && !spectrum.hasCalibration());
        QVERIFY(chooseGroup("O", "K"));
        QCOMPARE(groups->item(0, 0)->text(), QStringLiteral("O_K"));
        QCOMPARE(hfShells->text(), QStringLiteral("O:K"));
        QVERIFY(chooseGroup("O", "K", true));
        QCOMPARE(groups->rowCount(), 1);
        QVERIFY(nextStep->isEnabled());
        QVERIFY(!measureCompton->isChecked() && !tabs->isTabEnabled(5));
        QVERIFY(!measureCompton->isVisible());
        advancedToggle->click();
        QVERIFY(measureCompton->isVisible());
        QVERIFY(!composition->isEnabled() && !hfTable->isEnabled() && !chooseHfShells->isEnabled());
        measureCompton->click();
        QVERIFY(tabs->isTabEnabled(5) && composition->isEnabled() && hfTable->isEnabled());
        QVERIFY(chooseHfShells->isEnabled() && enableHf->isChecked());
        enableHf->setChecked(false);
        composition->setText("C:1, O:2");
        bool usedBundledTable = false;
        QTimer::singleShot(0, &spectrum, [&] {
            auto picker = qobject_cast<QDialog *>(QApplication::activeModalWidget());
            if (!picker)
                return;
            if (auto choices = picker->findChild<QListWidget *>("SpectroscopyHfShellChoices")) {
                bool carbon = false, oxygen = false;
                for (int i = 0; i < choices->count(); ++i) {
                    auto item = choices->item(i);
                    const auto key = item->data(Qt::UserRole).toString();
                    if (key == "C:1") {
                        carbon = item->text().contains("284.7");
                        item->setCheckState(Qt::Checked);
                    }
                    if (key == "O:1") {
                        oxygen = item->text().contains("543.1") && item->checkState() == Qt::Checked;
                        item->setCheckState(Qt::Unchecked);
                    }
                }
                usedBundledTable = carbon && oxygen && picker->windowTitle().contains("O_K");
            }
            picker->done(usedBundledTable ? QDialog::Accepted : QDialog::Rejected);
        });
        chooseHfShells->click();
        QVERIFY(usedBundledTable);
        QCOMPARE(hfShells->text(), QStringLiteral("C:1"));
        auto scans = spectrum.findChild<QListWidget *>("XrsFiles");
        auto addScans = spectrum.findChild<QPushButton *>("AddSpectrumScans");
        QVERIFY(scans && addScans);
        spectrum.setAvailableFiles({scan, elastic, reverse, scan, temp.filePath("4_grid_scan_20261106000000_1.nxs")});
        auto chooseScanAction = [&](const QString &name, bool fromList, const QStringList &selection = {}) {
            QElapsedTimer elapsed;
            elapsed.start();
            bool actionChosen = false, pickerAccepted = selection.isEmpty();
            QTimer drive;
            connect(&drive, &QTimer::timeout, &spectrum, [&] {
                if (auto *picker = qobject_cast<QDialog *>(QApplication::activeModalWidget()); picker != &spectrum &&
                    picker && picker->objectName() == "SpectrumBrowserScanPicker") {
                    auto *list = picker->findChild<QListWidget *>("BrowserSpectrumScans");
                    auto *add = picker->findChild<QPushButton *>("AddSelectedBrowserScans");
                    if (!actionChosen || elapsed.hasExpired(4000) || !list || !add || list->count() != 3 ||
                        add->isEnabled()) {
                        picker->reject();
                        return;
                    }
                    for (int i = 0; i < list->count(); ++i)
                        list->item(i)->setSelected(selection.contains(list->item(i)->data(Qt::UserRole).toString()));
                    pickerAccepted = list->selectedItems().size() == selection.size() && add->isEnabled();
                    if (pickerAccepted)
                        add->click();
                    else
                        picker->reject();
                    return;
                }
                auto *menu = qobject_cast<QMenu *>(QApplication::activePopupWidget());
                if (!menu)
                    return;
                if (elapsed.hasExpired(4000)) {
                    menu->close();
                    return;
                }
                for (auto *action : menu->actions()) {
                    if (action->objectName() == name && action->isEnabled()) {
                        actionChosen = true;
                        QTimer::singleShot(0, menu, [menu, action] {
                            menu->setActiveAction(action);
                            QTest::keyClick(menu, Qt::Key_Return);
                        });
                        return;
                    }
                }
                if (!selection.isEmpty() && !menu->actions().isEmpty() && menu->actions().first()->menu()) {
                    menu->setActiveAction(menu->actions().first());
                    QTest::keyClick(menu, Qt::Key_Right);
                } else {
                    menu->close();
                }
            });
            drive.start(20);
            if (fromList) {
                const auto point = scans->count() ? scans->visualItemRect(scans->item(0)).center() : QPoint(5, 5);
                QContextMenuEvent event(QContextMenuEvent::Mouse, point, scans->viewport()->mapToGlobal(point));
                QApplication::sendEvent(scans->viewport(), &event);
            } else {
                addScans->click();
            }
            return actionChosen && pickerAccepted;
        };
        QVERIFY(chooseScanAction("ChooseBrowserSpectrumScans", false, {scan, reverse}));
        QCOMPARE(scans->count(), 2);
        QCOMPARE(scans->item(0)->checkState(), Qt::Checked);
        QCOMPARE(scans->item(1)->checkState(), Qt::Checked);
        scans->selectAll();
        QVERIFY(chooseScanAction("UncheckSelectedSpectrumScans", true));
        QCOMPARE(groups->item(0, 1)->text(), QStringLiteral("0 / 2"));
        QVERIFY(chooseScanAction("CheckSelectedSpectrumScans", true));
        QCOMPARE(groups->item(0, 1)->text(), QStringLiteral("2 / 2"));
        scans->item(1)->setCheckState(Qt::Unchecked);
        spectrum.findChild<QComboBox *>("ElasticEnergy")->setEditText("/entry/energy_custom");
        spectrum.findChild<QDoubleSpinBox *>("GlitchDropPercent")->setValue(7.5);
        spectrum.findChild<QTableWidget *>("RoiCorrections")->item(0, 1)->setText("2");
        spectrum.findChild<QLineEdit *>("BackgroundWindows_polynomial")->setText("0:1,9:10");
        spectrum.findChild<QComboBox *>("BackgroundScope")->setCurrentIndex(1);
        spectrum.findChild<QLineEdit *>("BackgroundProtectedWindows")->setText("4:6");
        spectrum.findChild<QDoubleSpinBox *>("CorrectedEnergyMinimum")->setValue(2);
        spectrum.findChild<QCheckBox *>("NormalizeComptonNearEdge")->setChecked(true);
        QVERIFY(chooseGroup("Ni", "L3"));
        QCOMPARE(scans->count(), 0);
        QCOMPARE(hfShells->text(), QStringLiteral("Ni:L3"));
        composition->setText("Ni:1, O:2");
        bool linkedNickel = false;
        QTimer::singleShot(0, &spectrum, [&] {
            auto picker = qobject_cast<QDialog *>(QApplication::activeModalWidget());
            if (!picker)
                return;
            if (auto choices = picker->findChild<QListWidget *>("SpectroscopyHfShellChoices")) {
                for (int i = 0; i < choices->count(); ++i) {
                    const auto item = choices->item(i);
                    if (item->checkState() == Qt::Checked)
                        linkedNickel = item->data(Qt::UserRole).toString() == "Ni:3" && item->text().contains("L23");
                }
            }
            picker->accept();
        });
        chooseHfShells->click();
        QVERIFY(linkedNickel);
        QCOMPARE(hfShells->text(), QStringLiteral("Ni:L3"));
        hfShells->setText("O:1");
        const auto customHfTable = temp.filePath("custom-profiles.dat");
        hfTable->setText(customHfTable);
        QVERIFY(restoreHfTable->isEnabled());
        restoreHfTable->click();
        QVERIFY(hfTable->text().isEmpty());
        hfTable->setText(customHfTable);
        QVERIFY(chooseScanAction("ChooseBrowserSpectrumScans", true, {reverse}));
        spectrum.findChild<QDoubleSpinBox *>("GlitchDropPercent")->setValue(3.5);
        spectrum.findChild<QLineEdit *>("BackgroundWindows_polynomial")->setText("4:8");
        spectrum.findChild<QTableWidget *>("RoiCorrections")->item(0, 1)->setText("4");
        QCOMPARE(groups->item(1, 0)->text(), QStringLiteral("Ni_L3"));
        QCOMPARE(groups->item(1, 1)->text(), QStringLiteral("1 / 1"));
        measureCompton->click();
        QVERIFY(!measureCompton->isChecked() && !tabs->isTabEnabled(5));
        groups->setCurrentCell(0, 0);
        QVERIFY(measureCompton->isChecked() && !tabs->isTabEnabled(5) && !enableHf->isChecked());
        QCOMPARE(scans->count(), 2);
        QCOMPARE(groups->item(0, 1)->text(), QStringLiteral("1 / 2"));
        const auto screenshots = qEnvironmentVariable("XRI_SPECTROSCOPY_SCREENSHOTS");
        if (!screenshots.isEmpty())
            QVERIFY(spectrum.grab().save(screenshots + "/spectroscopy-edge-groups.png"));
        auto save = tabs->cornerWidget()->findChild<QPushButton *>("SaveSpectrumConfig");
        QVERIFY(save && save->isVisible() && !save->icon().isNull());
        chooseWorkflowFile(save, configPath);
        QVERIFY2(QFileInfo::exists(configPath), qPrintable(spectrum.statusText()));
        spectrum.selectStep(1);
        QVERIFY(!save->isVisible());

        SpectroscopyDialog restored;
        restored.show();
        auto import = restored.findChild<QPushButton *>("ImportSpectrumConfig");
        QVERIFY(import && import->isVisible() && !import->icon().isNull());
        chooseWorkflowFile(import, configPath);
        QCOMPARE(restored.findChild<QLineEdit *>("ElasticFile")->text(), elastic);
        QCOMPARE(restored.findChild<QLineEdit *>("SpectrumRoiScheme")->text(), roiPath);
        auto restoredScans = restored.findChild<QListWidget *>("XrsFiles");
        QCOMPARE(restoredScans->count(), 2);
        QCOMPARE(restoredScans->item(0)->data(Qt::UserRole).toString(), scan);
        QCOMPARE(restoredScans->item(1)->data(Qt::UserRole).toString(), reverse);
        QCOMPARE(restoredScans->item(0)->checkState(), Qt::Checked);
        QCOMPARE(restoredScans->item(1)->checkState(), Qt::Unchecked);
        QCOMPARE(restored.findChild<QComboBox *>("ElasticEnergy")->currentText(), QStringLiteral("/entry/energy_custom"));
        QCOMPARE(restored.findChild<QDoubleSpinBox *>("GlitchDropPercent")->value(), 7.5);
        QCOMPARE(restored.findChild<QTableWidget *>("RoiCorrections")->item(0, 1)->text(), QStringLiteral("2"));
        QCOMPARE(restored.findChild<QLineEdit *>("BackgroundWindows_polynomial")->text(), QStringLiteral("0:1,9:10"));
        QCOMPARE(restored.findChild<QComboBox *>("BackgroundScope")->currentIndex(), 1);
        QCOMPARE(restored.findChild<QLineEdit *>("BackgroundProtectedWindows")->text(), QStringLiteral("4:6"));
        QCOMPARE(restored.findChild<QDoubleSpinBox *>("CorrectedEnergyMinimum")->value(), 2.);
        auto *normalize = restored.findChild<QCheckBox *>("NormalizeComptonNearEdge");
        auto *restoredMeasureCompton = restored.findChild<QCheckBox *>("MeasureCompton");
        QVERIFY(restoredMeasureCompton && restoredMeasureCompton->isChecked());
        QVERIFY(!restored.findChild<QCheckBox *>("EnableBackground_hf_core")->isChecked());
        QVERIFY(normalize->isChecked());
        auto restoredGroups = restored.findChild<QTableWidget *>("SpectroscopyEdgeGroups");
        QCOMPARE(restoredGroups->rowCount(), 2);
        QCOMPARE(restoredGroups->currentRow(), 0);
        QCOMPARE(restored.findChild<QLineEdit *>("BackgroundHfComposition")->text(), QStringLiteral("C:1, O:2"));
        QCOMPARE(restored.findChild<QLineEdit *>("BackgroundHfShells")->text(), QStringLiteral("C:1"));
        QVERIFY(restored.findChild<QLineEdit *>("BackgroundHfTable")->text().isEmpty());
        restoredGroups->setCurrentCell(1, 0);
        QVERIFY(!restoredMeasureCompton->isChecked());
        QCOMPARE(restored.findChild<QLineEdit *>("BackgroundHfComposition")->text(), QStringLiteral("Ni:1, O:2"));
        QCOMPARE(restored.findChild<QLineEdit *>("BackgroundHfShells")->text(), QStringLiteral("O:1"));
        QCOMPARE(restored.findChild<QLineEdit *>("BackgroundHfTable")->text(), customHfTable);
        QCOMPARE(restoredScans->count(), 1);
        QCOMPARE(restoredScans->item(0)->data(Qt::UserRole).toString(), reverse);
        QCOMPARE(restored.findChild<QDoubleSpinBox *>("GlitchDropPercent")->value(), 3.5);
        QCOMPARE(restored.findChild<QLineEdit *>("BackgroundWindows_polynomial")->text(), QStringLiteral("4:8"));
        QCOMPARE(restored.findChild<QTableWidget *>("RoiCorrections")->item(0, 1)->text(), QStringLiteral("4"));
        restoredGroups->setCurrentCell(0, 0);
        QFile config(configPath);
        QVERIFY(config.open(QIODevice::ReadOnly));
        const auto grouped = QJsonDocument::fromJson(config.readAll()).object();
        config.close();
        QCOMPARE(grouped.value("version").toInt(), 2);
        const auto groupRecords = grouped.value("groups").toArray();
        QCOMPARE(groupRecords.at(0).toObject().value("target_edge").toObject().value("element").toString(),
                 QStringLiteral("O"));
        QCOMPARE(groupRecords.at(1).toObject().value("target_edge").toObject().value("shell").toString(),
                 QStringLiteral("L3"));
        const auto resavedPath = temp.filePath("restored-group-targets.config.json");
        QVERIFY(restored.saveConfiguration(resavedPath));
        QFile resaved(resavedPath);
        QVERIFY(resaved.open(QIODevice::ReadOnly));
        const auto resavedGroups = QJsonDocument::fromJson(resaved.readAll()).object().value("groups").toArray();
        QCOMPARE(resavedGroups.at(1).toObject().value("target_edge"), groupRecords.at(1).toObject().value("target_edge"));
        auto oldGroup = groupRecords.at(1).toObject();
        auto oldTarget = oldGroup.value("target_edge").toObject();
        auto oldParameters = oldGroup.value("parameters").toObject();
        oldTarget.remove("subtraction_initialized");
        oldParameters["background_hf_shells"] = "";
        oldGroup["target_edge"] = oldTarget;
        oldGroup["parameters"] = oldParameters;
        auto oldScheme = grouped;
        oldScheme["groups"] = QJsonArray{groupRecords.at(0), oldGroup};
        oldScheme["active_group_id"] = oldGroup.value("id");
        QVERIFY(config.open(QIODevice::WriteOnly | QIODevice::Truncate));
        QVERIFY(config.write(QJsonDocument(oldScheme).toJson()) > 0);
        config.close();
        QVERIFY(restored.importConfiguration(configPath));
        auto restoredShells = restored.findChild<QLineEdit *>("BackgroundHfShells");
        QCOMPARE(restoredShells->text(), QStringLiteral("Ni:L3"));
        restoredShells->clear();
        const auto clearedPath = temp.filePath("cleared-shells.config.json");
        QVERIFY(restored.saveConfiguration(clearedPath));
        QVERIFY(restored.importConfiguration(clearedPath));
        QVERIFY(restoredShells->text().isEmpty());
        auto legacy = groupRecords.first().toObject();
        legacy["schema"] = grouped.value("schema");
        legacy["version"] = 1;
        auto parameters = legacy.value("parameters").toObject();
        parameters.remove("compton_measured");
        parameters.remove("normalize_near_edge_by_scan_count");
        parameters.remove("background_scope");
        parameters.remove("background_protected_ev");
        parameters.remove("hf_table");
        parameters["background_enable_linear"] = false;
        parameters["low_q_linear_windows"] = parameters.take("background_windows_polynomial");
        for (const auto *key : {"background_enable_polynomial", "background_polynomial_order",
                                "result_smoothing_enabled", "result_smoothing_window_points", "result_smoothing_order"})
            parameters.remove(key);
        parameters["background_source"] = 1;
        legacy["parameters"] = parameters;
        QVERIFY(config.open(QIODevice::WriteOnly | QIODevice::Truncate));
        QVERIFY(config.write(QJsonDocument(legacy).toJson()) > 0);
        config.close();
        QVERIFY(restored.importConfiguration(configPath));
        QCOMPARE(restoredGroups->rowCount(), 1);
        QVERIFY(restoredMeasureCompton->isChecked());
        QVERIFY(restored.findChild<QLineEdit *>("BackgroundHfTable")->text().isEmpty());
        QVERIFY(!normalize->isChecked());
        QVERIFY(!restored.findChild<QCheckBox *>("EnableBackground_polynomial")->isChecked());
        QCOMPARE(restored.findChild<QLineEdit *>("BackgroundWindows_polynomial")->text(), QStringLiteral("0:1,9:10"));
        QCOMPARE(restored.findChild<QSpinBox *>("BackgroundPolynomialOrder")->value(), 1);
        QVERIFY(!restored.findChild<QCheckBox *>("CompareSmoothedSpectrum")->isChecked());
        QCOMPARE(restored.findChild<QComboBox *>("BackgroundScope")->currentIndex(), 0);
        QVERIFY(!restored.importConfiguration(roiPath));
        QCOMPARE(restoredScans->count(), 2);
    }

    void workspaceSaving() {
        QTemporaryDir temp;
        QVERIFY(temp.isValid());
        const QString key = "workspace/root";
        QSettings settings;
        const auto previous = settings.value(key);
        const auto restore = qScopeGuard([&] {
            if (previous.isValid())
                settings.setValue(key, previous);
            else
                settings.remove(key);
            settings.sync();
        });
        settings.remove(key);
        settings.sync();

        MainWindow window;
        window.show();
        const auto config = ImagingFixture::create(temp.filePath("data"), false, false, false, true);
        auto imaging = window.openImaging();
        imaging->applyConfiguration(config);
        auto ixs = window.openBondImaging();
        ixs->applyConfiguration(config);
        auto select = window.findChild<QAction *>("SpecifyWorkspace");
        QVERIFY(select && QDir(temp.path()).mkdir("Work"));
        chooseWorkflowFile(select, temp.filePath("Work"));
        const QDir workspace(temp.filePath("Work"));
        QCOMPARE(WorkspacePaths::root(), workspace.canonicalPath());
        QVERIFY(QDir(temp.path()).mkdir(QStringLiteral("中文")));
        QVERIFY_EXCEPTION_THROWN(WorkspacePaths::setRoot(temp.filePath(QStringLiteral("中文"))), std::exception);
        QCOMPARE(WorkspacePaths::root(), workspace.canonicalPath());
        QCOMPARE(imaging->findChild<QLabel *>("ImagingWorkspaceOutput")->text(),
                 QDir::toNativeSeparators(workspace.filePath("imaging/results")));

        const auto elastic = SpectrumFixture::create(temp.path(), true);
        const auto scan = SpectrumFixture::create(temp.path(), false);
        RoiDialog editor(elastic, "Workspace test", &window, false);
        QTRY_VERIFY_WITH_TIMEOUT(!editor.isLoading(), 10000);
        QVERIFY(editor.loadFrom(temp.filePath("test.roi.json")));
        auto saveRoi = editor.findChild<QPushButton *>("SaveCrystalRois");
        QVERIFY(saveRoi);
        saveRoi->click();
        const auto savedRoi = editor.currentSchemePath();
        QVERIFY(QFileInfo(savedRoi).isFile());
        saveRoi->click();
        QCOMPARE(QDir(workspace.filePath("ROI/config")).entryList({"*.roi.json"}, QDir::Files).size(), 2);
        QVERIFY(QDir(workspace.filePath("ROI/results")).exists());
        const auto name = WorkspacePaths::newPath(WorkspacePaths::Workflow::Roi, WorkspacePaths::Kind::Result,
                                                    QStringLiteral("中文/样品 A 42.nxs"), "spots.json");
        QVERIFY(QFileInfo(name).fileName().startsWith("A_42_spots_"));

        window.openFile(elastic);
        QTRY_VERIFY_WITH_TIMEOUT(!window.isReading(), 10000);
        auto spectrum = window.openSpectroscopy();
        QVERIFY(spectrum->loadRois(savedRoi));
        QVERIFY(spectrum->addEdgeGroup("O", "K"));
        spectrum->findChild<QLineEdit *>("SpectroscopySampleName")->setText("Sample A");
        spectrum->addScans({scan});
        spectrum->findChild<QPushButton *>("SaveSpectrumConfig")->click();
        const auto savedPlans = QDir(workspace.filePath("spectroscopy/config")).entryList({"*.config.json"}, QDir::Files);
        QCOMPARE(savedPlans.size(), 1);
        QVERIFY(savedPlans.first().startsWith("Sample_A_spectroscopy_groups_"));
        for (int step : {0, 1, 2}) {
            spectrum->runStep(step);
            QTRY_VERIFY_WITH_TIMEOUT(!spectrum->isBusy(), 10000);
        }
        QVERIFY2(spectrum->hasResult(), qPrintable(spectrum->statusText()));
        spectrum->findChild<QPushButton *>("SaveSpectrum")->click();
        const QDir spectra(workspace.filePath("spectroscopy/results"));
        const auto spectraRuns = spectra.entryList(QDir::Dirs | QDir::NoDotAndDotDot);
        QCOMPARE(spectraRuns.size(), 1);
        QVERIFY(spectraRuns.first().startsWith("Sample_A_O_K_raw_"));
        QVERIFY(QFileInfo(spectra.filePath(spectraRuns.first() + "/spectrum_data.txt")).isFile());

        window.openImaging(0);
        imaging->findChild<QPushButton *>("SaveImagingConfig")->click();
        QCOMPARE(QDir(workspace.filePath("imaging/config")).entryList({"*.config.json"}, QDir::Files).size(), 1);
        imaging->start();
        QTRY_VERIFY_WITH_TIMEOUT(!imaging->isBusy(), 15000);
        QVERIFY2(imaging->hasResult(), qPrintable(imaging->statusText()));
        QCOMPARE(QFileInfo(imaging->result().directory).absolutePath(), workspace.filePath("imaging/results"));
        QVERIFY(QFileInfo(QDir(imaging->result().directory).filePath("reconstruction.h5")).isFile());
        QVERIFY(!QFileInfo::exists(temp.filePath("data/output")));

        window.openBondImaging();
        ixs->findChild<QPushButton *>("SaveBondConfig")->click();
        QCOMPARE(QDir(workspace.filePath("ixs/config")).entryList({"*.config.json"}, QDir::Files).size(), 1);
        for (const auto *button : {"CheckBondInputs", "PrepareBondSpatial", "ContinueBondSpatial",
                                  "ConfirmBondCrystals", "PrepareBondChannels", "RunBondVolume", "SaveBondVolume"}) {
            auto control = ixs->findChild<QPushButton *>(button);
            QVERIFY(control && control->isVisible() && control->isEnabled());
            control->click();
            QTRY_VERIFY_WITH_TIMEOUT(!ixs->isBusy(), 10000);
        }
        const QDir volumes(workspace.filePath("ixs/results"));
        const auto files = volumes.entryList({"*.h5"}, QDir::Files);
        QCOMPARE(files.size(), 1);
        QCOMPARE(BondImaging::load(volumes.filePath(files.first())).axes[0].size(), 2);
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
        auto stepHint = work->findChild<QLabel *>("BondStepHint");
        auto previous = work->findChild<QPushButton *>("PreviousBondStep");
        auto *spatialSplit = work->findChild<QSplitter *>("BondSpatialSplitter");
        auto *spatialToggle = work->findChild<QToolButton *>("BondSpatialSplitterToggle");
        auto *fileStatus = work->findChild<QLabel *>("BondVolumeFileStatus");
        QVERIFY(tabs && stepHint);
        QVERIFY(previous && spatialSplit && spatialToggle && fileStatus);
        for (int stage = 1; stage < tabs->count(); ++stage) {
            QVERIFY(tabs->isTabEnabled(stage));
            QTest::mouseClick(tabs->tabBar(), Qt::LeftButton, Qt::NoModifier,
                              tabs->tabBar()->tabRect(stage).center());
            QCOMPARE(tabs->currentIndex(), stage);
            QVERIFY(stepHint->isVisible());
        }
        work->selectStage(BondWorkspace::Inputs);
        auto import = tabs->cornerWidget()->findChild<QPushButton *>("ImportBondConfig");
        auto save = tabs->cornerWidget()->findChild<QPushButton *>("SaveBondConfig");
        QVERIFY(import && save && import->isVisible() && save->isVisible());
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
            {
                QFile legacy(configPath);
                QVERIFY(legacy.open(QIODevice::ReadWrite));
                auto object = QJsonDocument::fromJson(legacy.readAll()).object();
                object["bond_contrast_output_file"] = "previous-result.h5";
                const auto bytes = QJsonDocument(object).toJson();
                QVERIFY(legacy.resize(0) && legacy.seek(0));
                QCOMPARE(legacy.write(bytes), bytes.size());
            }
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
            auto confirm = work->findChild<QPushButton *>("ConfirmBondCrystals");
            auto calibrate = work->findChild<QPushButton *>("PrepareBondSpatial");
            auto next = work->findChild<QPushButton *>("ContinueBondSpatial");
            auto channels = work->findChild<QPushButton *>("PrepareBondChannels");
            auto channelTable = work->findChild<QTableWidget *>("BondChannelInputs");
            auto buildSplit = work->findChild<QSplitter *>("BondVolumeBuildSplitter");
            auto saveVolume = work->findChild<QPushButton *>("SaveBondVolume");
            QVERIFY(validate && run && confirm && calibrate && next && channels && channelTable && buildSplit && saveVolume);
            QCOMPARE(buildSplit->orientation(), Qt::Horizontal);
            QVERIFY(buildSplit->widget(0)->isAncestorOf(channels) && buildSplit->widget(0)->isAncestorOf(run));
            QVERIFY(buildSplit->widget(1)->isAncestorOf(channelTable));
            QVERIFY(!run->isVisible());
            QCOMPARE(validate->text(), QStringLiteral("校验数据"));
            work->selectStage(BondWorkspace::Spatial);
            QVERIFY(!calibrate->isEnabled());
            work->selectStage(BondWorkspace::Inputs);
            validate->click();
            QVERIFY(work->isBusy() && tabs->tabBar()->isEnabled() && tabs->isTabEnabled(BondWorkspace::Results));
            work->selectStage(BondWorkspace::Results);
            QCOMPARE(tabs->currentIndex(), BondWorkspace::Results);
            QTRY_VERIFY_WITH_TIMEOUT(!work->isBusy(), 10000);
            QCOMPARE(tabs->currentIndex(), BondWorkspace::Spatial);
            QVERIFY(!import->isVisible() && !save->isVisible());
            QCOMPARE(window.statusBar()->currentMessage(), work->statusText());
            QCOMPARE(spatialSplit->orientation(), Qt::Horizontal);
            QVERIFY(spatialSplit->widget(0)->isAncestorOf(work->findChild<QLineEdit *>("BondXStep")));
            QVERIFY(spatialSplit->widget(1)->isAncestorOf(work->findChild<QTableWidget *>("BondSpatialCrystals")));
            QTRY_VERIFY(spatialSplit->sizes().first() < spatialSplit->width() / 3);
            const auto parameterWidth = spatialSplit->sizes().first();
            auto *handle = spatialSplit->handle(1);
            const auto center = handle->rect().center();
            QTest::mousePress(handle, Qt::LeftButton, Qt::NoModifier, center);
            QTest::mouseMove(handle, center + QPoint(100, 0));
            QTest::mouseRelease(handle, Qt::LeftButton, Qt::NoModifier, center);
            QTRY_VERIFY(spatialSplit->sizes().first() > parameterWidth + 40);
            spatialSplit->setSizes({parameterWidth, spatialSplit->width() - parameterWidth});
            spatialToggle->click();
            QTRY_VERIFY(!spatialSplit->widget(0)->isVisible() && spatialToggle->isVisible());
            spatialToggle->click();
            QTRY_VERIFY(spatialSplit->widget(0)->isVisible());
            QVERIFY(calibrate->isEnabled() && !confirm->isEnabled() && !channels->isEnabled());
            QVERIFY(!next->isEnabled());
            QCOMPARE(calibrate->parentWidget(), next->parentWidget());
            calibrate->click();
            QTRY_VERIFY_WITH_TIMEOUT(!work->isBusy(), 10000);
            QVERIFY2(next->isEnabled(), qPrintable(work->findChild<QLabel *>("BondStatus")->text()));
            work->selectStage(BondWorkspace::VolumeBuild);
            QVERIFY(!channels->isEnabled() && stepHint->text().contains(QStringLiteral("第 3 步")));
            work->selectStage(BondWorkspace::Spatial);
            if (!dir.isEmpty())
                QVERIFY(window.grab().save(QDir(dir).filePath("spatial-calibration.png")));
            next->click();
            QCOMPARE(tabs->currentIndex(), BondWorkspace::Crystals);
            auto crystalTable = work->findChild<QTableWidget *>("BondCrystalSelection");
            QVERIFY(crystalTable && crystalTable->rowCount());
            auto crystal = crystalTable->item(0, 0);
            QVERIFY(crystal && crystal->checkState() == Qt::Checked);
            crystal->setCheckState(Qt::Unchecked);
            QCOMPARE(tabs->currentIndex(), BondWorkspace::Crystals);
            QVERIFY(!confirm->isEnabled() && previous->isVisible());
            crystal->setCheckState(Qt::Checked);
            QVERIFY(!confirm->isEnabled());
            previous->click();
            QCOMPARE(tabs->currentIndex(), BondWorkspace::Spatial);
            QVERIFY(!next->isEnabled());
            calibrate->click();
            QTRY_VERIFY_WITH_TIMEOUT(!work->isBusy(), 10000);
            next->click();
            QCOMPARE(tabs->currentIndex(), BondWorkspace::Crystals);
            confirm->click();
            QCOMPARE(tabs->currentIndex(), BondWorkspace::VolumeBuild);
            work->selectStage(BondWorkspace::Volume);
            QVERIFY(stepHint->text().contains(QStringLiteral("第 4 步")));
            work->selectStage(BondWorkspace::VolumeBuild);
            QVERIFY(run->isVisible() && !run->isEnabled() && channelTable->rowCount() == 0);
            channels->click();
            QTRY_VERIFY_WITH_TIMEOUT(!work->isBusy(), 10000);
            QCOMPARE(channelTable->rowCount(), 2);
            if (!dir.isEmpty())
                QVERIFY(window.grab().save(QDir(dir).filePath("channel-projections.png")));
            QVERIFY(run->isVisible() && run->isEnabled());
            run->click();
            QTRY_VERIFY_WITH_TIMEOUT(!work->isBusy(), 10000);
            QVERIFY(work->hasResult() && !QFile::exists(output->text()));
            QCOMPARE(tabs->currentIndex(), BondWorkspace::Volume);
            QVERIFY(!stepHint->isVisible());
            QVERIFY(saveVolume->isVisible() && saveVolume->isEnabled());
            QVERIFY(fileStatus->isVisible() && fileStatus->text().contains(QStringLiteral("尚未保存")));
            QVERIFY(fileStatus->parentWidget() == work && fileStatus->y() > tabs->geometry().bottom());
            QVERIFY(saveVolume->mapTo(work, QPoint{}).y() > tabs->geometry().bottom());
            if (!dir.isEmpty())
                QVERIFY(window.grab().save(QDir(dir).filePath("voxel-unsaved.png")));
            saveVolume->click();
            QTRY_VERIFY_WITH_TIMEOUT(!work->isBusy(), 10000);
            QVERIFY(work->hasResult() && QFile::exists(output->text()));
            QVERIFY(!saveVolume->isVisible());
            QVERIFY(fileStatus->text().contains(QFileInfo(output->text()).fileName()));
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
        auto voxelSplit = work->findChild<QSplitter *>("BondVoxelSplitter");
        QVERIFY(spectrum && slice && figure && style && source && channel && plotTabs && voxelSplit);
        auto voxelSettings = work->findChild<QGroupBox *>("BondVoxelPlotSettings");
        QVERIFY(voxelSettings && voxelSettings->isAncestorOf(source) && voxelSettings->isAncestorOf(style));
        QVERIFY(!voxelSplit->isAncestorOf(fileStatus) && fileStatus->isVisible());
        QCOMPARE(voxelSplit->orientation(), Qt::Horizontal);
        QVERIFY(voxelSplit->widget(0)->isAncestorOf(source) && voxelSplit->widget(0)->isAncestorOf(style));
        QVERIFY(voxelSplit->widget(1)->isAncestorOf(slice) && voxelSplit->widget(1)->isAncestorOf(spectrum));
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
        auto analyze = work->findChild<QPushButton *>("AnalyzeBondVolume");
        QVERIFY(analyze && analyze->isVisible());
        analyze->click();
        QCOMPARE(tabs->currentIndex(), BondWorkspace::Results);
        auto mode = work->findChild<QComboBox *>("BondContrastMode");
        auto map = work->findChild<QLabel *>("BondMapInfo");
        auto montage = work->findChild<QWidget *>("BondMontage");
        auto ea = work->findChild<QDoubleSpinBox *>("BondEa");
        auto integralLo = work->findChild<QDoubleSpinBox *>("BondIntegralLo");
        auto integralHi = work->findChild<QDoubleSpinBox *>("BondIntegralHi");
        auto resultSplit = work->findChild<QSplitter *>("BondResultsSplitter");
        auto figureTabs = work->findChild<QTabWidget *>("BondResultPlots");
        QVERIFY(mode && map && montage && ea && integralLo && integralHi && resultSplit && figureTabs);
        QCOMPARE(resultSplit->orientation(), Qt::Horizontal);
        QVERIFY(resultSplit->widget(0)->isAncestorOf(mode) && resultSplit->widget(0)->isAncestorOf(channel));
        QVERIFY(resultSplit->widget(1)->isAncestorOf(montage));
        QCOMPARE(mode->count(), 3);
        QVERIFY(channel->isVisible() && !ea->isVisible() && !integralLo->isVisible());
        const int previousChannel = channel->value();
        channel->setValue(0);
        QVERIFY(map->text().contains(QString::number(spectrum->curves().first().x.first(), 'g', 9)));
        channel->setValue(previousChannel);
        mode->setCurrentIndex(1);
        QVERIFY(!channel->isVisible() && ea->isVisible());
        QVERIFY(map->text().contains(QStringLiteral("有效体素")));
        if (!dir.isEmpty())
            QVERIFY(window.grab().save(QDir(dir).filePath("contrast-sections.png")));
        mode->setCurrentIndex(2);
        QCOMPARE(mode->currentText(), QStringLiteral("区间积分强度"));
        QVERIFY(integralLo->isVisible() && integralHi->isVisible() && !ea->isVisible() && !channel->isVisible());
        QVERIFY2(map->text().contains(QStringLiteral("有效体素")), qPrintable(map->text()));
        const double loKev = integralLo->value(), hiKev = integralHi->value();
        integralLo->setValue(hiKev);
        QVERIFY(map->text().contains(QStringLiteral("积分区间必须递增")));
        integralLo->setValue(loKev);
        QVERIFY(map->text().contains(QStringLiteral("有效体素")));
        auto showTransfer = work->findChild<QCheckBox *>("BondShowEnergyTransfer");
        auto analyzer = work->findChild<QDoubleSpinBox *>("BondAnalyzerReference");
        QVERIFY(showTransfer && analyzer);
        if (!showTransfer->isChecked()) {
            analyzer->setValue(loKev);
            showTransfer->setChecked(true);
            QVERIFY(std::abs(integralLo->value() - 1000 * (loKev - analyzer->value())) < 1e-7);
            QVERIFY(std::abs(integralHi->value() - 1000 * (hiKev - analyzer->value())) < 1e-7);
            QVERIFY(map->text().contains(QStringLiteral("有效体素")));
            showTransfer->setChecked(false);
        }
        if (!dir.isEmpty())
            QVERIFY(window.grab().save(QDir(dir).filePath("integral-sections.png")));
        auto surface = work->findChild<QWidget *>("BondIsosurface");
        QVERIFY(surface && resultSplit->widget(1)->isAncestorOf(surface));
        figureTabs->setCurrentWidget(surface);
        auto build = work->findChild<QPushButton *>("BuildBondSurface");
        QVERIFY(build && build->isVisible() && resultSplit->widget(0)->isAncestorOf(build));
        QVERIFY(!work->findChild<QSpinBox *>("BondFirstZ")->isVisible());
        build->click();
        QTRY_VERIFY_WITH_TIMEOUT(!work->isBusy(), 10000);
        auto status = work->findChild<QLabel *>("BondStatus");
        QVERIFY(status->text().contains(QStringLiteral("三角形")));
        if (!dir.isEmpty())
            QVERIFY(window.grab().save(QDir(dir).filePath("integral-isosurface.png")));
        work->selectStage(BondWorkspace::Volume);
        plotTabs->setCurrentIndex(0);
        QVERIFY(voxelSplit->widget(0)->isAncestorOf(channel) && channel->isVisible());
        QCOMPARE(channel->value(), previousChannel);
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
        auto page = work->findChild<QWidget *>("BondReferenceMaps");
        auto figures = work->findChild<QTabWidget *>("BondResultPlots");
        auto split = work->findChild<QSplitter *>("BondResultsSplitter");
        QVERIFY(run && page && figures && split);
        figures->setCurrentWidget(page);
        QVERIFY(run->isVisible());
        QVERIFY(split->widget(0)->isAncestorOf(run) && split->widget(1)->isAncestorOf(page));
        QVERIFY(!work->findChild<QComboBox *>("BondViewSource")->isVisible());
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
        window.show();
        auto stack = window.findChild<QStackedWidget *>("ProcessingWorkspaces");
        QVERIFY(stack);
        auto browser = stack->currentWidget();
        auto *automaticAction = window.findChild<QAction *>("DetectSpotRois");
        auto *manualAction = window.findChild<QAction *>("SelectCrystalRois");
        auto *convertTextRois = window.findChild<QAction *>("ConvertTextRois");
        auto *backToBrowser = window.findChild<QToolButton *>("BackToBrowser");
        auto *menuBar = window.menuBar();
        auto *roiMenu = window.findChild<QMenu *>("RoiMenu");
        QVERIFY(automaticAction && manualAction && convertTextRois && backToBrowser && roiMenu);
        QCOMPARE(roiMenu->title(), QStringLiteral("ROI(&R)"));
        QCOMPARE(manualAction->text(), QStringLiteral("手动选取 ROI"));
        QCOMPARE(automaticAction->text(), QStringLiteral("自动选取 ROI"));
        QVERIFY(menuBar->actions().contains(roiMenu->menuAction()));
        QVERIFY(!menuBar->actions().contains(manualAction) && !menuBar->actions().contains(automaticAction));
        QCOMPARE(roiMenu->actions(), (QList<QAction *>{manualAction, automaticAction, convertTextRois}));
        QVERIFY(!manualAction->menu() && !automaticAction->menu());
        auto chooseRoi = [&](QAction *action) {
            QTest::mouseClick(menuBar, Qt::LeftButton, Qt::NoModifier,
                              menuBar->actionGeometry(roiMenu->menuAction()).center());
            QTRY_VERIFY(roiMenu->isVisible());
            QTest::mouseClick(roiMenu, Qt::LeftButton, Qt::NoModifier, roiMenu->actionGeometry(action).center());
            QTRY_VERIFY(!roiMenu->isVisible());
        };
        QVERIFY(!automaticAction->isEnabled());
        automaticAction->trigger();
        QCOMPARE(stack->currentWidget(), browser);
        window.openFile(file);
        QTRY_VERIFY_WITH_TIMEOUT(!window.isReading(), 10000);
        QVERIFY(automaticAction->isEnabled());
        chooseRoi(automaticAction);
        auto *automatic = qobject_cast<RoiDetectionDialog *>(stack->currentWidget());
        QVERIFY(automatic && !automatic->isWindow() && !QApplication::activeModalWidget());
        QVERIFY(!automatic->findChild<QLineEdit *>("SpotRoiSource"));
        auto *sourceStatus = automatic->findChild<QLabel *>("SpotRoiStatus");
        QVERIFY(sourceStatus->text().contains(QFileInfo(file).fileName()));
        QVERIFY(automatic->findChild<QPushButton *>("RunSpotRoiDetection")->isEnabled());
        const auto second = SpectrumFixture::create(temp.path(), false);
        window.openFile(second);
        QTRY_VERIFY_WITH_TIMEOUT(!window.isReading(), 10000);
        QVERIFY(sourceStatus->text().contains(QFileInfo(second).fileName()));
        auto *files = window.findChild<QComboBox *>("BrowserFiles");
        const int selectedFile = files->currentIndex();
        files->setCurrentIndex(-1);
        QVERIFY(window.currentFile().isEmpty() && !automaticAction->isEnabled());
        files->setCurrentIndex(selectedFile);
        QVERIFY(!window.currentFile().isEmpty() && automaticAction->isEnabled());
        window.closeCurrentFile();
        QVERIFY(automaticAction->isEnabled());
        window.closeCurrentFile();
        QVERIFY(!automaticAction->isEnabled());
        QVERIFY(!automatic->findChild<QPushButton *>("RunSpotRoiDetection")->isEnabled());
        window.openFile(file);
        QTRY_VERIFY_WITH_TIMEOUT(!window.isReading(), 10000);
        QVERIFY(automaticAction->isEnabled());
        chooseRoi(manualAction);
        auto *manual = qobject_cast<RoiDialog *>(stack->currentWidget());
        QVERIFY(manual && !manual->isWindow());
        automaticAction->trigger();
        QCOMPARE(stack->currentWidget(), automatic);
        backToBrowser->click();
        QCOMPARE(stack->currentWidget(), browser);
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
            QVERIFY(spec->addEdgeGroup("C", "K"));
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
        window.openBondImaging();
        capture("ixs-contrast");
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
                                      QString::fromUtf8("IXS衬度成像"), QString::fromUtf8("语言"),
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
    void databaseFolderImport() {
        QTemporaryDir temp;
        QVERIFY(temp.isValid());
        const auto source = SpectrumFixture::create(temp.path(), true);
        const QDir root(temp.filePath("database"));
        QVERIFY(root.mkpath("nested/level2"));
        QVERIFY(root.mkpath("directory.h5"));
        const QStringList paths{"01.nxs", "02.h5", "03.HDF5", "nested/04.hdf", "nested/level2/05.nx5"};
        for (const auto &path : paths)
            QVERIFY(QFile::copy(source, root.filePath(path)));
        QVERIFY(QFile::copy(source, root.filePath("nested/ignored.txt")));

        MainWindow window;
        window.show();
        auto action = window.findChild<QAction *>("ImportDatabaseHdfFiles");
        QVERIFY(action && window.findChild<QMenu *>("FileMenu")->actions().contains(action));
        const bool nonNative = QApplication::testAttribute(Qt::AA_DontUseNativeDialogs);
        QApplication::setAttribute(Qt::AA_DontUseNativeDialogs);
        QSettings settings;
        const auto previous = settings.value("databaseDirectory");
        const auto restore = qScopeGuard([&] {
            QApplication::setAttribute(Qt::AA_DontUseNativeDialogs, nonNative);
            if (previous.isValid())
                settings.setValue("databaseDirectory", previous);
            else
                settings.remove("databaseDirectory");
        });
        QTimer respond;
        QElapsedTimer deadline;
        deadline.start();
        bool pickedFolder = false, testedSelection = false;
        bool progressStarted = false;
        int found = 0, afterShift = 0, afterCtrl = 0;
        QVector<int> progressSteps;
        QString nestedDirectory;
        connect(&respond, &QTimer::timeout, &window, [&] {
            auto *dialog = qobject_cast<QDialog *>(QApplication::activeModalWidget());
            if (!dialog)
                return;
            if (deadline.hasExpired(5000)) {
                dialog->reject();
                return;
            }
            if (auto *picker = qobject_cast<QFileDialog *>(dialog)) {
                if (pickedFolder) {
                    QMetaObject::invokeMethod(picker, "accept");
                } else {
                    picker->setDirectory(temp.path());
                    picker->selectFile("database");
                    if (auto *name = picker->findChild<QLineEdit *>("fileNameEdit"))
                        name->setText("database");
                    pickedFolder = true;
                }
                return;
            }
            if (testedSelection)
                return;
            auto *list = dialog->findChild<QTableView *>("DatabaseHdfFiles");
            auto *import = dialog->findChild<QPushButton *>("ImportSelectedDatabaseFiles");
            auto *progress = dialog->findChild<QProgressBar *>("DatabaseImportProgress");
            if (!list || !import || !progress || list->model()->rowCount() == 0)
                return;
            found = list->model()->rowCount();
            if (found != paths.size()) {
                dialog->reject();
                return;
            }
            const auto click = [list](int row, Qt::KeyboardModifiers modifiers) {
                const auto index = list->model()->index(row, 0);
                list->scrollTo(index);
                QTest::mouseClick(list->viewport(), Qt::LeftButton, modifiers, list->visualRect(index).center());
            };
            nestedDirectory = list->model()->index(4, 1).data().toString();
            click(0, Qt::NoModifier);
            click(2, Qt::ShiftModifier);
            afterShift = list->selectionModel()->selectedRows().size();
            click(1, Qt::ControlModifier);
            click(4, Qt::ControlModifier);
            afterCtrl = list->selectionModel()->selectedRows().size();
            testedSelection = import->isEnabled();
            if (testedSelection) {
                connect(progress, &QProgressBar::valueChanged, &window, [&](int value) {
                    if (value > 0)
                        progressSteps.append(value);
                });
                import->click();
                progressStarted = progress->isVisible() && progress->value() == 0 && progress->maximum() == afterCtrl;
            } else {
                dialog->reject();
            }
        });
        respond.start(20);
        action->trigger();
        respond.stop();
        QVERIFY(pickedFolder && testedSelection);
        QCOMPARE(found, 5);
        QCOMPARE(afterShift, 3);
        QCOMPARE(afterCtrl, 3);
        QVERIFY(progressStarted);
        QCOMPARE(progressSteps, (QVector<int>{1, 2, 3}));
        QCOMPARE(nestedDirectory, QDir::toNativeSeparators("nested/level2"));
        QTRY_VERIFY_WITH_TIMEOUT(!window.isReading(), 10000);
        QCOMPARE(window.openedFiles(),
                 (QStringList{root.filePath(paths[0]), root.filePath(paths[2]), root.filePath(paths[4])}));
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
        auto *localStatus = dialog->findChild<QLabel *>("ImagingStatus");
        QVERIFY(localStatus && localStatus->isHidden());
        QCOMPARE(window.statusBar()->currentMessage(), dialog->statusText());
        auto steps = dialog->findChild<QTabWidget *>("ImagingSteps");
        auto avizo = dialog->findChild<QCheckBox *>("output.save_avizo");
        QVERIFY(steps && avizo && !avizo->isChecked());
        QCOMPARE(steps->tabText(3), QStringLiteral("4 剖面分析"));
        QCOMPARE(steps->tabText(4), QStringLiteral("5 三维模型"));
        QCOMPARE(steps->tabText(5), QStringLiteral("6 日志"));
        auto buildModel = dialog->findChild<QPushButton *>("BuildImagingModel");
        auto saveGif = dialog->findChild<QPushButton *>("SaveImagingModelGif");
        QVERIFY(buildModel && !buildModel->isEnabled() && saveGif && !saveGif->isEnabled());
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
        dialog->selectStep(5);
        auto *logSplit = steps->widget(5)->findChild<QSplitter *>("ImagingLogSplitter");
        auto *logToggle = steps->widget(5)->findChild<QToolButton *>("ImagingLogSplitterToggle");
        auto *logChoices = steps->widget(5)->findChild<QListWidget *>("ImagingLogSections");
        auto *log = steps->widget(5)->findChild<QPlainTextEdit *>("ImagingLog");
        auto *summary = steps->widget(5)->findChild<QPlainTextEdit *>("ImagingSummary");
        QVERIFY(logSplit && logToggle && logChoices && log && summary);
        QCOMPARE(logSplit->orientation(), Qt::Horizontal);
        QVERIFY(logSplit->widget(0)->isAncestorOf(logChoices) && logSplit->widget(1)->isAncestorOf(log));
        QCOMPARE(logChoices->currentRow(), 0);
        QCOMPARE(logChoices->currentItem()->text(), QStringLiteral("运行日志"));
        QVERIFY(log->isVisible() && log->isReadOnly() && log->toPlainText().contains("Completed:"));
        const auto savedLog = log->toPlainText();
        logChoices->setCurrentRow(1);
        QVERIFY(!log->isVisible() && summary->isVisible() && summary->toPlainText().contains("shape_xyz"));
        logToggle->click();
        QTRY_VERIFY(!logChoices->isVisible());
        QVERIFY(summary->isVisible());
        logToggle->click();
        QTRY_VERIFY(logChoices->isVisible());
        logChoices->setCurrentRow(0);
        if (const auto screenshots = qEnvironmentVariable("XRI_WORKFLOW_SCREENSHOTS"); !screenshots.isEmpty())
            QVERIFY(window.grab().save(screenshots + "/imaging-log.png"));
        dialog->selectStep(3);
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
        auto profile = steps->widget(3)->findChild<QSplitter *>("ImagingProfileSplitter");
        auto plane = dialog->findChild<QComboBox *>("ImagingPlane");
        auto display = dialog->findChild<QComboBox *>("ImagingDisplay");
        auto slice = dialog->findChild<QSpinBox *>("ImagingSlice");
        QVERIFY(profile && plane && display && slice);
        QCOMPARE(profile->orientation(), Qt::Horizontal);
        QCOMPARE(profile->count(), 2);
        QVERIFY(profile->widget(0)->isAncestorOf(plane) && profile->widget(0)->isAncestorOf(slice));
        QVERIFY(profile->widget(1)->isAncestorOf(image));
        auto profileToggle = profile->findChild<QToolButton *>("ImagingProfileSplitterToggle");
        QVERIFY(profileToggle);
        profileToggle->click();
        QTRY_VERIFY(!profile->widget(0)->isVisible() && image->isVisible());
        profileToggle->click();
        QTRY_VERIFY(plane->isVisible());
        plane->setCurrentText("XY");
        display->setCurrentIndex(1);
        QVERIFY(!slice->isEnabled() && !image->data().values.isEmpty());
        display->setCurrentIndex(0);
        QVERIFY(slice->isEnabled());
        slice->setValue(slice->maximum());
        QVERIFY(!image->data().values.isEmpty());
        if (const auto screenshots = qEnvironmentVariable("XRI_WORKFLOW_SCREENSHOTS"); !screenshots.isEmpty())
            QVERIFY(dialog->grab().save(screenshots + "/profile-analysis.png"));
        dialog->selectStep(4);
        auto modelSplit = steps->widget(4)->findChild<QSplitter *>("ImagingModelSplitter");
        auto model = dialog->findChild<QWidget *>("ImagingModel");
        auto modelSource = dialog->findChild<QComboBox *>("ImagingModelSource");
        auto threshold = dialog->findChild<QLineEdit *>("ImagingModelThreshold");
        QVERIFY(modelSplit && model && modelSource && threshold && model->isVisible() && buildModel->isEnabled());
        QCOMPARE(modelSplit->orientation(), Qt::Horizontal);
        QCOMPARE(modelSplit->count(), 2);
        QVERIFY(modelSplit->widget(0)->isAncestorOf(modelSource) &&
                modelSplit->widget(0)->isAncestorOf(threshold) &&
                modelSplit->widget(0)->isAncestorOf(saveGif));
        QVERIFY(modelSplit->widget(1)->isAncestorOf(model));
        QCOMPARE(modelSource->count(), 3);
        const auto validThreshold = threshold->text();
        auto depth = dialog->findChild<QLabel *>("ImagingDepthReadout");
        auto validateDepth = dialog->findChild<QPushButton *>("ValidateImagingDepth");
        auto depthReport = dialog->findChild<QPushButton *>("ImagingDepthReport");
        QVERIFY(depth && validateDepth && depthReport && validateDepth->isEnabled());
        const QString initialDepth = depth->text();
        QTest::mouseClick(model, Qt::LeftButton, Qt::NoModifier, QPoint(86, model->height() - 64));
        QVERIFY(depth->text() != initialDepth && depth->text().contains("X = -1 mm"));
        validateDepth->click();
        QTRY_VERIFY_WITH_TIMEOUT(!dialog->isBusy(), 30000);
        QVERIFY(depthReport->isEnabled() && depth->text().contains("FWHM"));
        auto comparePsf = dialog->findChild<QPushButton *>("CompareImagingPsf");
        auto psfSelection = dialog->findChild<QLabel *>("ImagingPsfSelection");
        QVERIFY(comparePsf && comparePsf->isEnabled() && psfSelection);
        dialog->selectStep(5);
        QVERIFY(!comparePsf->isVisible());
        logChoices->setCurrentRow(3);
        QVERIFY(comparePsf->isVisible() && logSplit->widget(0)->isAncestorOf(comparePsf) && psfSelection->isVisible());
        comparePsf->click();
        QTRY_VERIFY_WITH_TIMEOUT(!dialog->isBusy(), 30000);
        QVERIFY(psfSelection->text().contains(QStringLiteral("本次标定比较推荐")));
        QVERIFY(QFileInfo::exists(temp.path() + "/output/psf_model_comparison.json"));
        QVERIFY(QFileInfo::exists(temp.path() + "/output/depth_validation.json"));
        logChoices->setCurrentRow(0);
        dialog->selectStep(4);
        threshold->setText("-1");
        buildModel->click();
        QVERIFY(!dialog->isBusy() && dialog->statusText().contains(QStringLiteral("阈值须介于")));
        threshold->setText(validThreshold);
        buildModel->click();
        QTRY_VERIFY_WITH_TIMEOUT(!dialog->isBusy(), 10000);
        QVERIFY2(dialog->statusText().contains(QStringLiteral("等值面已生成")), qPrintable(dialog->statusText()));
        if (const auto screenshots = qEnvironmentVariable("XRI_WORKFLOW_SCREENSHOTS"); !screenshots.isEmpty())
            QVERIFY(dialog->grab().save(screenshots + "/model-analysis.png"));
        // Start on the background so the gesture cannot select the projected depth axis.
        const QPoint origin(24, 24), end = origin + QPoint(95, 30);
        const auto initial = model->grab().toImage();
        const auto selectedDepth = depth->text();
        QTest::mousePress(model, Qt::LeftButton, Qt::NoModifier, origin);
        QMouseEvent drag(QEvent::MouseMove, end, model->mapToGlobal(end), Qt::NoButton, Qt::LeftButton,
                         Qt::NoModifier);
        QApplication::sendEvent(model, &drag);
        QTest::mouseRelease(model, Qt::LeftButton, Qt::NoModifier, end);
        QVERIFY(model->grab().toImage() != initial);
        QCOMPARE(depth->text(), selectedDepth);
        QTest::mouseDClick(model, Qt::LeftButton, Qt::NoModifier, origin);
        QCOMPARE(model->grab().toImage(), initial);

        QVERIFY(saveGif->isEnabled());
        const auto gifBase = temp.filePath(QStringLiteral("模型旋转"));
        chooseWorkflowFile(saveGif, gifBase, 2.5);
        QTRY_VERIFY_WITH_TIMEOUT(!dialog->isBusy(), 15000);
        QVERIFY2(dialog->statusText().contains(QStringLiteral("旋转 GIF 已保存")), qPrintable(dialog->statusText()));
        const auto viewAfterExport = model->grab().toImage();
        QTest::mouseDClick(model, Qt::LeftButton, Qt::NoModifier, origin);
        QCOMPARE(model->grab().toImage(), viewAfterExport);
        const auto gifPath = gifBase + ".gif";
        QVERIFY2(QFileInfo::exists(gifPath), qPrintable(dialog->statusText()));
        QImageReader gif(gifPath);
        QVERIFY2(gif.canRead(), qPrintable(gif.errorString()));
        const int frames = gif.imageCount();
        QCOMPARE(frames, 60);
        QCOMPARE(gif.loopCount(), -1);
        QImage firstFrame;
        bool moved = false;
        int duration = 0;
        for (int i = 0; i < frames; ++i) {
            const auto frame = gif.read().convertToFormat(QImage::Format_RGB32);
            QVERIFY2(!frame.isNull(), qPrintable(gif.errorString()));
            duration += gif.nextImageDelay();
            if (i == 0) {
                QCOMPARE(frame.size(), QSize(720, 720));
                firstFrame = frame;
                int surfacePixels = 0;
                for (int y = 0; y < frame.height(); ++y) {
                    const auto *line = reinterpret_cast<const QRgb *>(frame.constScanLine(y));
                    for (int x = 0; x < frame.width(); ++x) {
                        const QRgb pixel = line[x];
                        if (pixel == qRgb(255, 255, 255))
                            continue;
                        QVERIFY(qRed(pixel) < qGreen(pixel) && qGreen(pixel) < qBlue(pixel));
                        ++surfacePixels;
                    }
                }
                QVERIFY(surfacePixels > 0);
            } else {
                moved = moved || frame != firstFrame;
            }
        }
        QVERIFY(moved);
        QCOMPARE(duration, 2500);

        const auto canceledGif = temp.filePath("canceled.gif");
        chooseWorkflowFile(saveGif, canceledGif, 6.);
        auto cancel = dialog->findChild<QPushButton *>("CancelImaging");
        QVERIFY(dialog->isBusy() && cancel && cancel->isEnabled());
        cancel->click();
        QTRY_VERIFY_WITH_TIMEOUT(!dialog->isBusy(), 15000);
        QVERIFY2(dialog->statusText().contains(QStringLiteral("已取消 GIF 导出")), qPrintable(dialog->statusText()));
        QVERIFY(!QFileInfo::exists(canceledGif) && saveGif->isEnabled());

        dialog->hide();
        window.openImaging(3);
        QCOMPARE(window.findChild<ImagingDialog *>(), dialog);
        QVERIFY(dialog->hasResult());
        auto file = dialog->findChild<QLineEdit *>("input.sample_nxs");
        QVERIFY(file);
        file->setText(temp.path() + "/other.nxs");
        QVERIFY(!dialog->hasResult() && image->data().values.isEmpty());
        QVERIFY(dialog->loadResult(temp.path() + "/output"));
        QCOMPARE(log->toPlainText(), savedLog);
        QVERIFY(dialog->hasResult());
        QVERIFY(depthReport->isEnabled() && depth->text().contains("FWHM"));
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
    void imagingWeightedBoundaryWorkflow() {
        QTemporaryDir tmp;
        auto c = ImagingFixture::create(tmp.path(), false, false, false, false, true, true);
        auto stats = c["statistics"].toObject();
        stats["variance_weighted"] = true;
        stats["sigma_floor"] = 1.;
        stats["background_yz_mm"] = QJsonArray{0, 0, -1, 0};
        c["statistics"] = stats;
        auto output = c["output"].toObject();
        output["save_diagnostics"] = false;
        c["output"] = output;
        ImagingDialog dialog;
        dialog.show();
        dialog.applyConfiguration(c);
        auto w = dialog.findChild<QCheckBox *>("statistics.variance_weighted");
        QVERIFY(w && w->isChecked());
        QCOMPARE(dialog.configuration()["statistics"].toObject(), stats);
        dialog.start();
        QTRY_VERIFY_WITH_TIMEOUT(!dialog.isBusy(), 30000);
        QVERIFY2(dialog.hasResult(), qPrintable(dialog.statusText()));
        QVERIFY(dialog.result().summary["sample_background_subtracted"].toBool());
        auto run = dialog.findChild<QPushButton *>("ValidateImagingBoundaries");
        auto view = dialog.findChild<QPushButton *>("ShowImagingBoundaries");
        QVERIFY(run && run->isEnabled() && view && !view->isEnabled());
        run->click();
        QTRY_VERIFY_WITH_TIMEOUT(!dialog.isBusy(), 30000);
        QVERIFY2(view->isEnabled(), qPrintable(dialog.statusText()));
        const QString savedDirectory = dialog.result().directory;
        QVERIFY(dialog.loadResult(savedDirectory));
        QVERIFY(view->isEnabled());
        bool inspected = false;
        QTimer::singleShot(0, &dialog, [&] {
            auto report = qobject_cast<QDialog *>(QApplication::activeModalWidget());
            if (!report)
                return;
            auto table = report->findChild<QTableWidget *>();
            auto plot = report->findChild<ScientificPlot *>();
            inspected = table && table->rowCount() == 8 && plot && !plot->data().values.isEmpty();
            report->accept();
        });
        view->click();
        QVERIFY(inspected);
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
        if(QFileInfo::exists(path+"/boundary_stability.json")) {
            auto show=dialog.findChild<QPushButton *>("ShowImagingBoundaries");
            QVERIFY(show && show->isEnabled());
            bool inspected=false;
            QTimer::singleShot(0, &dialog, [&] {
                auto report = qobject_cast<QDialog *>(QApplication::activeModalWidget());
                if (!report)
                    return;
                auto table = report->findChild<QTableWidget *>();
                auto slice = report->findChild<QSpinBox *>();
                auto plot = report->findChild<ScientificPlot *>();
                if (table && table->rowCount() > 0 && slice && plot) {
                    for (int i = 0; i <= slice->maximum() && plot->data().values.isEmpty(); ++i)
                        slice->setValue(i);
                    inspected = !plot->data().values.isEmpty();
                    const auto snapshots = qEnvironmentVariable("XRI_IMAGING_SCREENSHOTS");
                    if (!snapshots.isEmpty() && QDir().mkpath(snapshots))
                        report->grab().save(snapshots + "/boundaries.png");
                }
                report->accept();
            });
            show->click();
            QVERIFY(inspected);
        }
        auto directory = qEnvironmentVariable("XRI_IMAGING_SCREENSHOTS");
        if (!directory.isEmpty()) {
            QVERIFY(QDir().mkpath(directory));
            auto steps = dialog.findChild<QTabWidget *>("ImagingSteps");
            auto buildModel = dialog.findChild<QPushButton *>("BuildImagingModel");
            QVERIFY(steps && buildModel);
            dialog.selectStep(4);
            buildModel->click();
            QTRY_VERIFY_WITH_TIMEOUT(!dialog.isBusy(), 15000);
            QVERIFY2(dialog.statusText().contains(QStringLiteral("等值面已生成")),
                     qPrintable(dialog.statusText()));
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
        QVERIFY(dialog.addEdgeGroup("C", "K"));
        auto *advance = dialog.findChild<QPushButton *>("NextSpectrumStep");
        auto *pages = dialog.findChild<QTabWidget *>("SpectroscopySteps");
        advance->click();
        QTRY_VERIFY_WITH_TIMEOUT(!dialog.isBusy(), 30000);
        QVERIFY2(pages->currentIndex() == 1, qPrintable(dialog.statusText()));
        QVERIFY(!dialog.hasCalibration());
        QVERIFY(dialog.findChild<ScientificPlot *>("LambdaMaskPreview")->data().sum > 0);
        advance->click();
        QTRY_VERIFY_WITH_TIMEOUT(!dialog.isBusy(), 15000);
        QVERIFY2(dialog.hasCalibration(), qPrintable(dialog.statusText()));
        auto fits = dialog.findChild<QTableWidget *>("ElasticFits");
        QVERIFY(fits);
        // Calibration already selects the first accepted crystal.
        const auto elasticPlot = pages->widget(2)->findChild<SpectrumPlot *>();
        QVERIFY(elasticPlot && !elasticPlot->curves().first().x.isEmpty());
        QCOMPARE(elasticPlot->curves().first().draw, SpectrumCurve::Draw::Samples);
        QCOMPARE(elasticPlot->curves().last().draw, SpectrumCurve::Draw::Line);
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
        const auto dir = QFileInfo(file).dir();
        const auto files = dir.entryList({"14517_scan_*.nxs", "14524_scan_*.nxs"}, QDir::Files, QDir::Name);
        if (files.size() == 2) {
            dialog.addScans({dir.filePath(files[0]), dir.filePath(files[1])});
            dialog.selectStep(3);
            dialog.findChild<QPushButton *>("ExtractSpectrum")->click();
            QTRY_VERIFY_WITH_TIMEOUT(!dialog.isBusy(), 30000);
            QVERIFY2(dialog.hasScans(), qPrintable(dialog.statusText()));
            auto *scan = dialog.findChild<QComboBox *>("ExtractedXrsScan");
            auto *plot = dialog.findChild<SpectrumPlot *>("I0Spectrum");
            const int badPoint[] = {120, 36};
            for (int i = 0; i < 2; ++i) {
                scan->setCurrentIndex(i);
                const auto &curves = plot->curves();
                QCOMPARE(curves.size(), 3);
                const int p = badPoint[i];
                QCOMPARE(curves[2].x, QVector<double>{curves[0].x[p]});
                QVERIFY(curves[1].y[p] > curves[0].y[p]);
                qInfo() << files[i] << "point" << p + 1 << curves[0].y[p] << "->" << curves[1].y[p];
            }
            if (!directory.isEmpty())
                QVERIFY(dialog.grab().save(directory + "/real-i0.png"));
        }
    }
    void spectralMaskRealScan() {
        const QString elastic = qEnvironmentVariable("XRI_MASK_ELASTIC");
        const QString reference = qEnvironmentVariable("XRI_MASK_REFERENCE");
        if (elastic.isEmpty() || reference.isEmpty())
            QSKIP("Optional real spectral-mask scans");
        std::atomic_bool cancel{false};
        for (const QString &detector : {QStringLiteral("D_LAMBDA"), QStringLiteral("D_MINIPIX")}) {
            const auto projection = HdfFile::sumFrames(elastic, "/entry/data/" + detector, cancel);
            QVERIFY2(projection.error.isEmpty(), qPrintable(projection.error));
            QVERIFY(!projection.values.isEmpty());
            const int peak = int(std::max_element(projection.values.cbegin(), projection.values.cend()) -
                                 projection.values.cbegin());
            const QSize size(projection.width, projection.height);
            const QRect box = QRect(peak % size.width() - 6, peak / size.width() - 6, 13, 13)
                                  .intersected(QRect(QPoint(), size));
            CrystalRoi roi{detector, detector, box, RoiState::Valid, {}};
            SpectrumSettings settings;
            settings.normalizeXrs = !qEnvironmentVariableIsSet("XRI_MASK_COUNTS_ONLY");
            SpectralMaskDialog dialog(settings, {roi}, {{detector, size}}, 0, {reference}, nullptr);
            dialog.show();
            dialog.findChild<QPushButton *>("DecomposePixels")->click();
            const auto *components = dialog.findChild<QComboBox *>("SampleComponent");
            const auto counts = HdfFile::sumFrames(reference, "/entry/data/" + detector, cancel);
            QVERIFY2(counts.error.isEmpty(), qPrintable(counts.error));
            const auto pixels = roi.pixels();
            const bool hasCounts = std::any_of(pixels.cbegin(), pixels.cend(), [&](const QPoint &p) {
                return counts.values[p.y() * size.width() + p.x()] > 0;
            });
            if (!hasCounts) {
                QCOMPARE(components->count(), 0);
                QVERIFY(!dialog.findChild<QPushButton *>("ApplySpectralMask")->isEnabled());
                QVERIFY(dialog.findChild<QPushButton *>("KeepOriginalRoi")->isEnabled());
                QVERIFY(!dialog.findChild<QPushButton *>("UpdateOriginalRoi")->isEnabled());
                QCOMPARE(dialog.findChild<QLabel *>("SpectralMaskStatus")->text(),
                         QStringLiteral("候选像素在参考扫描中没有信号。"));
                continue;
            }
            QVERIFY2(components->count() == 2,
                     qPrintable(dialog.findChild<QLabel *>("SpectralMaskStatus")->text()));
            const auto selection = dialog.selection();
            QVERIFY(!selection.pixels.isEmpty());
            qInfo() << detector << box << selection.provenance;
            const auto directory = qEnvironmentVariable("XRI_SPECTROSCOPY_SCREENSHOTS");
            if (!directory.isEmpty()) {
                QDir().mkpath(directory);
                QVERIFY(dialog.grab().save(directory + "/spectral-" + detector + ".png"));
            }
        }
    }
    void spectralMaskWorkflow() {
        QTemporaryDir temp;
        const auto elastic = SpectrumFixture::create(temp.path(), true);
        const auto reference = SpectrumFixture::create(temp.path(), false);
        MainWindow window;
        window.openFile(elastic);
        QTRY_VERIFY_WITH_TIMEOUT(!window.isReading(), 10000);
        auto &workbench = *window.openSpectroscopy();
        QVERIFY(workbench.loadRois(temp.filePath("test.roi.json")));
        QVERIFY(workbench.addEdgeGroup("O", "K"));
        workbench.addScans({reference});
        auto *advance = workbench.findChild<QPushButton *>("NextSpectrumStep");
        auto *calculateMask = workbench.findChild<QPushButton *>("CalculateSpectrumMask");
        auto *resetMask = workbench.findChild<QPushButton *>("ResetSpectrumMask");
        auto *maskFraction = workbench.findChild<QDoubleSpinBox *>("MaskFilterFraction");
        auto *localStatus = workbench.findChild<QLabel *>("SpectrumStatus");
        auto *pages = workbench.findChild<QTabWidget *>("SpectroscopySteps");
        auto *source = workbench.findChild<QComboBox *>("MaskPreviewSource");
        QVERIFY(advance && pages && source && calculateMask && resetMask && maskFraction && localStatus);
        QVERIFY(calculateMask->isHidden() && resetMask->isHidden() && localStatus->isHidden());
        advance->click();
        QTRY_VERIFY_WITH_TIMEOUT(!workbench.isBusy(), 10000);
        QCOMPARE(pages->currentIndex(), 1);
        QVERIFY(calculateMask->isVisible() && calculateMask->geometry().right() < advance->geometry().left());
        QVERIFY(resetMask->isVisible() && resetMask->geometry().right() < calculateMask->geometry().left());
        QCOMPARE(window.statusBar()->currentMessage(), workbench.statusText());
        QVERIFY(!workbench.hasCalibration());
        QCOMPARE(source->currentIndex(), 1);
        auto *lambda = pages->widget(1)->findChild<ScientificPlot *>("LambdaMaskPreview");
        auto *mini = pages->widget(1)->findChild<ScientificPlot *>("MiniPixMaskPreview");
        QVERIFY(lambda && mini);
        QCOMPARE(lambda->data().sum, 4.);
        QCOMPARE(mini->data().sum, 2.);
        auto *corrections = workbench.findChild<QTableWidget *>("RoiCorrections");
        QVERIFY(corrections);
        corrections->setCurrentCell(5, 0); // Invalid row falls back directly to first valid crystal.
        auto *open = workbench.findChild<QPushButton *>("OpenSpectralMask");
        QVERIFY(open);
        bool finished = false;
        int revisedCount = 0;
        QTimer::singleShot(0, [&] {
            auto *dialog = workbench.findChild<SpectralMaskDialog *>("SpectralMaskDialog");
            if (!dialog) {
                if (auto *other = qobject_cast<QDialog *>(QApplication::activeModalWidget()))
                    other->reject();
                return;
            }
            auto close = qScopeGuard([dialog] { if (dialog->isVisible()) dialog->close(); });
            auto *previous = dialog->findChild<QPushButton *>("PreviousSpectralCrystal");
            auto *next = dialog->findChild<QPushButton *>("NextSpectralCrystal");
            auto *apply = dialog->findChild<QPushButton *>("ApplySpectralMask");
            auto *decomposeAll = dialog->findChild<QPushButton *>("DecomposeAllSpectralMasks");
            auto *applyBatch = dialog->findChild<QPushButton *>("ApplyBatchSpectralMask");
            auto *threshold = dialog->findChild<QDoubleSpinBox *>("CovarianceThreshold");
            auto *component = dialog->findChild<QComboBox *>("SampleComponent");
            auto *thresholdMode = dialog->findChild<QComboBox *>("CovarianceThresholdMode");
            auto *picker = dialog->findChild<QListWidget *>("SpectralCrystal");
            auto *update = dialog->findChild<QPushButton *>("UpdateOriginalRoi");
            auto *keep = dialog->findChild<QPushButton *>("KeepOriginalRoi");
            QVERIFY(picker && update && keep && decomposeAll && applyBatch && component && thresholdMode);
            QVERIFY(!update->isEnabled() && !keep->isEnabled() && !apply->isEnabled());
            QVERIFY(decomposeAll->isEnabled() && !applyBatch->isEnabled());
            QCOMPARE(picker->count(), 3);
            QVERIFY(!previous->isEnabled());
            QCOMPARE(dialog->currentRow(), 0);
            dialog->findChild<QPushButton *>("DecomposePixels")->click();
            QVERIFY(apply->isEnabled());
            next->click();
            QCOMPARE(dialog->currentRow(), 1);
            QCOMPARE(component->count(), 0);
            QVERIFY(!apply->isEnabled());
            decomposeAll->click();
            QVERIFY(apply->isEnabled() && !decomposeAll->isEnabled());
            const auto movedReference = reference + ".cache-test";
            QVERIFY(QFile::rename(reference, movedReference));
            const auto restoreReference = qScopeGuard([&] { QFile::rename(movedReference, reference); });
            previous->click();
            QCOMPARE(dialog->currentRow(), 0);
            QVERIFY(apply->isEnabled());
            auto *iterate = dialog->findChild<QPushButton *>("IterateSpectralMask");
            if (iterate->isEnabled())
                iterate->click();
            threshold->setValue(.15);
            const auto draft = dialog->selection();
            QVERIFY(!draft.pixels.isEmpty());
            dialog->activateWindow();
            QTRY_VERIFY(dialog->isActiveWindow());
            QTest::keyClick(dialog, Qt::Key_Right, Qt::AltModifier);
            QTRY_COMPARE(dialog->currentRow(), 1);
            QCOMPARE(dialog->findChild<QComboBox *>("SpectralReference")->currentText(), reference);
            QVERIFY(apply->isEnabled());
            threshold->setValue(.3);
            QTest::keyClick(dialog, Qt::Key_Left, Qt::AltModifier);
            QTRY_COMPARE(dialog->currentRow(), 0);
            QCOMPARE(threshold->value(), .15);
            QCOMPARE(dialog->selection().pixels, draft.pixels);
            QCOMPARE(dialog->selection().provenance, draft.provenance);
            apply->click();
            QCOMPARE(dialog->currentRow(), 1);
            QCOMPARE(threshold->value(), .3);
            QVERIFY(!update->isEnabled());
            QVERIFY(!corrections->item(0, 0)->text().endsWith(QStringLiteral("✓")));
            const auto directory = qEnvironmentVariable("XRI_SPECTROSCOPY_SCREENSHOTS");
            if (!directory.isEmpty()) {
                QDir().mkpath(directory);
                QVERIFY(dialog->grab().save(directory + "/crystal-navigation.png"));
            }
            next->click();
            QCOMPARE(dialog->currentRow(), 60); // Cross from LAMBDA to MINIPIX, skipping invalid ROIs.
            QVERIFY(!next->isEnabled());
            QVERIFY(apply->isEnabled());
            QVERIFY(!update->isEnabled());
            component->setCurrentIndex(1);
            thresholdMode->setCurrentIndex(1);
            threshold->setValue(-1e20);
            QVERIFY(applyBatch->isEnabled());
            applyBatch->click();
            QCOMPARE(dialog->currentRow(), 60);
            QVERIFY(update->isEnabled());
            QCOMPARE(dialog->acceptedMasks().size(), 3);
            QCOMPARE(dialog->acceptedMasks()[0].provenance, draft.provenance);
            for (int r : {1, 60}) {
                const auto batchRound = dialog->acceptedMasks()[r].provenance["rounds"].toArray().last().toObject();
                QCOMPARE(batchRound["target_component_1based"].toInt(), 2);
                QCOMPARE(batchRound["threshold_mode"].toString(), QStringLiteral("absolute"));
            }
            picker->setCurrentRow(1);
            QCOMPARE(component->currentIndex(), 1);
            QCOMPARE(threshold->value(), -1e20);
            keep->click();
            QVERIFY(!dialog->acceptedMasks().contains(1));
            picker->setCurrentRow(0);
            QVERIFY(update->isEnabled());
            threshold->setValue(.16);
            QVERIFY(!update->isEnabled());
            component->setCurrentIndex(1);
            thresholdMode->setCurrentIndex(1);
            threshold->setValue(-1e20);
            QVERIFY(apply->isEnabled());
            apply->click();
            QVERIFY(update->isEnabled());
            QCOMPARE(dialog->acceptedMasks().size(), 2);
            QCOMPARE(dialog->acceptedMasks()[0].provenance["rounds"].toArray().last().toObject()
                         .value("target_component_1based").toInt(), 2);
            picker->setCurrentRow(2);
            component->setCurrentIndex(0);
            thresholdMode->setCurrentIndex(0);
            threshold->setValue(.99);
            revisedCount = dialog->selection().pixels.size();
            QCOMPARE(revisedCount, 1);
            apply->click();
            QVERIFY(update->isEnabled());
            update->click();
            finished = true;
        });
        open->click();
        QVERIFY(finished);
        QCOMPARE(mini->data().sum, double(revisedCount));
        QVERIFY(window.statusBar()->currentMessage().contains(QStringLiteral("参数暂停生效")));
        QVERIFY(!corrections->item(1, 0)->text().endsWith(QStringLiteral("✓")));
        for (const int row : {0, 60})
            QVERIFY(corrections->item(row, 0)->text().endsWith(QStringLiteral("✓")));
        advance->click();
        QTRY_VERIFY_WITH_TIMEOUT(!workbench.isBusy(), 10000);
        QVERIFY2(workbench.hasCalibration(), qPrintable(workbench.statusText()));
        QCOMPARE(pages->currentIndex(), 2);
        auto *fits = pages->widget(2)->findChild<QTableWidget *>();
        QCOMPARE(fits->item(60, 1)->text(), QString("4 / %1").arg(revisedCount));
        QCOMPARE(fits->item(1, 1)->text(), QString("4 / 4"));
        workbench.selectStep(1);
        bool reviewedAgain = false;
        QTimer::singleShot(0, [&] {
            auto *dialog = workbench.findChild<SpectralMaskDialog *>();
            auto close = qScopeGuard([dialog] { dialog->reject(); });
            auto *run = dialog->findChild<QPushButton *>("DecomposePixels");
            auto *runAll = dialog->findChild<QPushButton *>("DecomposeAllSpectralMasks");
            auto *keep = dialog->findChild<QPushButton *>("KeepOriginalRoi");
            auto *update = dialog->findChild<QPushButton *>("UpdateOriginalRoi");
            auto *referenceBox = dialog->findChild<QComboBox *>("SpectralReference");
            QVERIFY(!update->isEnabled());
            runAll->click();
            for (int i = 0; i < 3; ++i) {
                QVERIFY(keep->isEnabled());
                keep->click();
            }
            QVERIFY(update->isEnabled());
            QVERIFY(dialog->acceptedMasks().isEmpty());
            referenceBox->setCurrentText(temp.filePath("missing.nxs"));
            run->click();
            QVERIFY(!keep->isEnabled() && !update->isEnabled());
            reviewedAgain = true;
        });
        open->click();
        QVERIFY(reviewedAgain && workbench.hasCalibration());

        maskFraction->setValue(.05);
        QCOMPARE(mini->data().sum, double(revisedCount));
        QVERIFY(workbench.hasCalibration());
        QVERIFY(!advance->isEnabled() && !pages->isTabEnabled(2));
        advance->click();
        QCOMPARE(pages->currentIndex(), 1);
        QVERIFY(window.statusBar()->currentMessage().contains(QStringLiteral("当前使用谱分解掩膜")));
        corrections->item(0, 1)->setText("invalid");
        calculateMask->click();
        QVERIFY(!workbench.isBusy() && workbench.hasCalibration() && !advance->isEnabled());
        QCOMPARE(mini->data().sum, double(revisedCount));
        QVERIFY(window.statusBar()->currentMessage().contains(QStringLiteral("ROI 修正必须为整数")));
        corrections->item(0, 1)->setText("0");

        auto *cancelMask = workbench.findChild<QPushButton *>("CancelSpectrumCalculation");
        QVERIFY(cancelMask);
        calculateMask->click();
        QVERIFY(workbench.isBusy() && cancelMask->isEnabled());
        cancelMask->click();
        QTRY_VERIFY_WITH_TIMEOUT(!workbench.isBusy(), 10000);
        QVERIFY(!advance->isEnabled() && workbench.hasCalibration());
        QCOMPARE(mini->data().sum, double(revisedCount));

        const auto plan = temp.filePath("spectral-mask-plan.json");
        QVERIFY(workbench.saveConfiguration(plan));
        auto verifySwitch = [&](const QString &configuration, bool pending) {
            SpectroscopyDialog restored;
            QVERIFY(restored.importConfiguration(configuration));
            auto *restoredNext = restored.findChild<QPushButton *>("NextSpectrumStep");
            auto *restoredTabs = restored.findChild<QTabWidget *>("SpectroscopySteps");
            restoredNext->click();
            QTRY_VERIFY_WITH_TIMEOUT(!restored.isBusy(), 10000);
            auto *restoredPreview = restored.findChild<ScientificPlot *>("MiniPixMaskPreview");
            QCOMPARE(restoredTabs->currentIndex(), 1);
            QCOMPARE(restoredNext->isEnabled(), !pending);
            if (pending)
                QVERIFY(restoredPreview->data().values.isEmpty() && !restoredTabs->isTabEnabled(2));
            else
                QCOMPARE(restoredPreview->data().sum, double(revisedCount));
            restored.findChild<QPushButton *>("CalculateSpectrumMask")->click();
            QTRY_VERIFY_WITH_TIMEOUT(!restored.isBusy(), 10000);
            QCOMPARE(restoredPreview->data().sum, 3.);
            QVERIFY(restoredNext->isEnabled() && restoredTabs->isTabEnabled(2));
        };
        const auto roiPath = temp.filePath("test.roi.json");
        const auto movedRoi = temp.filePath("moved.roi.json");
        QVERIFY(QFile::rename(roiPath, movedRoi));
        verifySwitch(plan, true);
        QVERIFY(QFile::rename(movedRoi, roiPath));

        QFile savedPlan(plan);
        QVERIFY(savedPlan.open(QIODevice::ReadOnly));
        auto oldPlan = QJsonDocument::fromJson(savedPlan.readAll()).object();
        auto groups = oldPlan["groups"].toArray();
        auto group = groups[0].toObject();
        QVERIFY(group.contains("base_roi_snapshot"));
        group.remove("base_roi_snapshot");
        group.remove("manual_calculation_pending");
        groups[0] = group;
        oldPlan["groups"] = groups;
        QFile legacy(temp.filePath("legacy-spectral-mask-plan.json"));
        QVERIFY(legacy.open(QIODevice::WriteOnly));
        QVERIFY(legacy.write(QJsonDocument(oldPlan).toJson()) > 0);
        legacy.close();
        verifySwitch(legacy.fileName(), false);

        resetMask->click();
        QTRY_VERIFY_WITH_TIMEOUT(!workbench.isBusy(), 10000);
        QVERIFY(!workbench.hasCalibration() && advance->isEnabled());
        QCOMPARE(lambda->data().sum, 8.);
        QCOMPARE(mini->data().sum, 4.);
        QCOMPARE(maskFraction->value(), .05);
        QCOMPARE(source->currentText(), QStringLiteral("原始 ROI 区域"));
        for (const int row : {0, 60})
            QVERIFY(!corrections->item(row, 0)->text().endsWith(QStringLiteral("✓")));
        QVERIFY(window.statusBar()->currentMessage().contains(QStringLiteral("当前使用原始 ROI")));
        if (const auto screenshots = qEnvironmentVariable("XRI_WORKFLOW_SCREENSHOTS"); !screenshots.isEmpty())
            QVERIFY(window.grab().save(screenshots + "/reset-masks.png"));
        const auto resetPlan = temp.filePath("original-roi-plan.json");
        QVERIFY(workbench.saveConfiguration(resetPlan));
        {
            SpectroscopyDialog restored;
            QVERIFY(restored.importConfiguration(resetPlan));
            auto *next = restored.findChild<QPushButton *>("NextSpectrumStep");
            next->click();
            QTRY_VERIFY_WITH_TIMEOUT(!restored.isBusy(), 10000);
            QCOMPARE(restored.findChild<ScientificPlot *>("MiniPixMaskPreview")->data().sum, 4.);
            next->click();
            QTRY_VERIFY_WITH_TIMEOUT(!restored.isBusy(), 10000);
            QVERIFY2(restored.hasCalibration(), qPrintable(restored.statusText()));
        }
        advance->click();
        QTRY_VERIFY_WITH_TIMEOUT(!workbench.isBusy(), 10000);
        QVERIFY2(workbench.hasCalibration(), qPrintable(workbench.statusText()));
        QCOMPARE(fits->item(60, 1)->text(), QStringLiteral("4 / 4"));
        advance->click();
        QTRY_VERIFY_WITH_TIMEOUT(!workbench.isBusy(), 10000);
        QVERIFY2(workbench.hasScans(), qPrintable(workbench.statusText()));
        workbench.selectStep(1);
        calculateMask->click();
        QTRY_VERIFY_WITH_TIMEOUT(!workbench.isBusy(), 10000);
        QVERIFY(!workbench.hasCalibration());
        QCOMPARE(pages->currentIndex(), 1);
        QVERIFY(advance->isEnabled() && pages->isTabEnabled(2));
        QCOMPARE(mini->data().sum, 3.);
        QCOMPARE(lambda->data().sum, 6.);
        QVERIFY(!corrections->item(60, 0)->text().endsWith(QStringLiteral("✓")));
        QVERIFY(window.statusBar()->currentMessage().contains(QStringLiteral("已应用参数掩膜")));
        resetMask->click();
        QTRY_VERIFY_WITH_TIMEOUT(!workbench.isBusy(), 10000);
        QCOMPARE(mini->data().sum, 4.);
        QCOMPARE(lambda->data().sum, 8.);
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
        QVERIFY(dialog->addEdgeGroup("O", "K"));
        dialog->addScans({xrs, reverse});
        auto next = dialog->findChild<QPushButton *>("NextSpectrumStep");
        auto edgeGroups = dialog->findChild<QTableWidget *>("SpectroscopyEdgeGroups");
        auto fitPeak = dialog->findChild<QPushButton *>("FitElasticPeaks");
        auto extract = dialog->findChild<QPushButton *>("ExtractSpectrum");
        auto combine = dialog->findChild<QPushButton *>("CombineSpectrum");
        auto save = dialog->findChild<QPushButton *>("SaveSpectrum");
        QVERIFY(next && fitPeak && extract && combine && save);
        QVERIFY(!dialog->findChild<QPushButton *>("CalibrateSpectrum"));
        QVERIFY(!extract->isEnabled() && !combine->isEnabled() && !save->isEnabled());
        dialog->selectStep(1);
        QVERIFY(next->isVisible() && !fitPeak->isVisible() && !extract->isVisible());
        QTest::mouseClick(next, Qt::LeftButton);
        QTRY_VERIFY_WITH_TIMEOUT(!dialog->isBusy(), 10000);
        QVERIFY2(dialog->hasCalibration(), qPrintable(dialog->statusText()));
        auto *tabs = dialog->findChild<QTabWidget *>("SpectroscopySteps");
        auto *elasticSplit = dialog->findChild<QSplitter *>("ElasticResultsSplitter");
        QVERIFY(elasticSplit && elasticSplit->orientation() == Qt::Horizontal);
        QVERIFY(fitPeak->isVisible() && fitPeak->x() < next->x());
        auto *fitThreshold = dialog->findChild<QDoubleSpinBox *>("ElasticMinimumRSquared");
        QVERIFY(fitThreshold);
        fitThreshold->setValue(.75);
        QVERIFY(!next->isEnabled() && !tabs->isTabEnabled(3));
        dialog->selectStep(3);
        QCOMPARE(tabs->currentIndex(), 2);
        fitPeak->click();
        QTRY_VERIFY_WITH_TIMEOUT(!dialog->isBusy(), 10000);
        QVERIFY(dialog->hasCalibration() && fitPeak->isVisible() && next->isEnabled());
        QVERIFY(extract->isEnabled());
        auto editor = window.openRoi(elastic, temp.path() + "/test.roi.json");
        QTRY_VERIFY_WITH_TIMEOUT(!editor->isLoading(), 10000);
        QVERIFY(editor->imagesReady());
        editor->findChild<QPushButton *>("ReturnFromRoi")->click();
        QVERIFY(dialog->hasCalibration());
        dialog->selectStep(2);
        QTest::mouseClick(next, Qt::LeftButton);
        QTRY_VERIFY_WITH_TIMEOUT(!dialog->isBusy(), 10000);
        QVERIFY2(dialog->hasScans(), qPrintable(dialog->statusText()));
        QCOMPARE(tabs->currentIndex(), 3);
        auto *plot = dialog->findChild<SpectrumPlot *>("I0Spectrum");
        auto *repair = dialog->findChild<QCheckBox *>("RepairXrsGlitches");
        auto *drop = dialog->findChild<QDoubleSpinBox *>("GlitchDropPercent");
        QVERIFY(repair->isVisible() && drop->value() == 5. && !fitPeak->isVisible());
        QCOMPARE(plot->curves().size(), 3);
        repair->click();
        QVERIFY(!dialog->hasScans() && !drop->isEnabled() && !next->isEnabled());
        extract->click();
        QTRY_VERIFY_WITH_TIMEOUT(!dialog->isBusy(), 10000);
        QVERIFY(next->isEnabled());
        QCOMPARE(plot->curves().size(), 2);
        QCOMPARE(plot->curves()[0].y, plot->curves()[1].y);
        repair->click();
        QVERIFY(!dialog->hasScans() && !next->isEnabled());
        extract->click();
        QTRY_VERIFY_WITH_TIMEOUT(!dialog->isBusy(), 10000);
        QVERIFY(next->isEnabled());
        QTest::mouseClick(next, Qt::LeftButton);
        QTRY_VERIFY_WITH_TIMEOUT(!dialog->isBusy(), 10000);
        QVERIFY2(dialog->hasResult(), qPrintable(dialog->statusText()));
        QVERIFY(edgeGroups);
        QCOMPARE(edgeGroups->item(0, 2)->text(), QStringLiteral("已叠加"));
        QCOMPARE(tabs->currentIndex(), 4);
        QCOMPARE(plot->curves().size(), 3);
        QVERIFY(save->isEnabled());
        auto *qMaximum = dialog->findChild<QDoubleSpinBox *>("SpectrumQMaximum");
        QVERIFY(qMaximum);
        qMaximum->setFocus();
        QTest::keyClick(qMaximum, Qt::Key_A, Qt::ControlModifier);
        QTest::keyClicks(qMaximum, "9.9");
        QCOMPARE(qMaximum->value(), 9.9);
        QVERIFY(!next->isEnabled() && !tabs->isTabEnabled(6) && !save->isEnabled());
        combine->click();
        QTRY_VERIFY_WITH_TIMEOUT(!dialog->isBusy(), 10000);
        QVERIFY(dialog->hasResult() && next->isEnabled());
        auto curve = dialog->findChild<SpectrumPlot *>("CombinedSpectrum");
        QVERIFY(curve && !curve->curves().isEmpty());
        auto *rawChoice = dialog->findChild<QComboBox *>("RawResultSpectrum");
        auto *rawPrevious = dialog->findChild<QToolButton *>("PreviousRawResultSpectrum");
        auto *rawNext = dialog->findChild<QToolButton *>("NextRawResultSpectrum");
        QVERIFY(rawChoice && rawPrevious && rawNext);
        QVERIFY(!rawPrevious->isEnabled() && rawNext->isEnabled());
        QCOMPARE(curve->curves().size(), dialog->result().modules.size() + 1);
        rawNext->click();
        QCOMPARE(rawChoice->currentText(), QStringLiteral("各模组谱"));
        QCOMPARE(curve->curves().size(), dialog->result().modules.size());
        for (const auto &module : curve->curves())
            QCOMPARE(module.y, dialog->result().modules.value(module.name));
        rawNext->click();
        QCOMPARE(rawChoice->currentText(), QStringLiteral("全部单晶体谱"));
        QCOMPARE(curve->curves().size(), 3);
        rawPrevious->click();
        QCOMPARE(rawChoice->currentText(), QStringLiteral("各模组谱"));
        rawChoice->setCurrentIndex(rawChoice->findData("VB"));
        QCOMPARE(curve->curves().size(), 1);
        QCOMPARE(curve->curves().first().y, dialog->result().modules.value("VB"));
        rawChoice->setCurrentIndex(rawChoice->count() - 1);
        QVERIFY(rawPrevious->isEnabled() && !rawNext->isEnabled());
        QCOMPARE(curve->curves().first().name, QStringLiteral("扫描叠加"));
        rawChoice->setCurrentIndex(0);
        auto *splitter = dialog->findChild<QSplitter *>("SpectrumResultsSplitter");
        auto *toggle = dialog->findChild<QToolButton *>("ToggleSpectrumParameters");
        auto *parameters = dialog->findChild<QWidget *>("SpectrumResultParameters");
        QVERIFY(splitter && toggle && parameters && parameters->isVisible());
        QCOMPARE(toggle->iconSize(), QSize(12, 12));
        QVERIFY(!toggle->icon().isNull() && toggle->autoRaise());
        QCOMPARE(splitter->orientation(), Qt::Horizontal);
        const int expandedWidth = splitter->width();
        toggle->click();
        QTRY_VERIFY(!parameters->isVisible() && splitter->width() > expandedWidth);
        QVERIFY(dialog->hasResult());
        const auto resultImages = qEnvironmentVariable("XRI_SPECTROSCOPY_SCREENSHOTS");
        if (!resultImages.isEmpty())
            QVERIFY(dialog->grab().save(resultImages + "/result-collapsed.png"));
        toggle->click();
        QTRY_VERIFY(parameters->isVisible());
        if (!resultImages.isEmpty())
            QVERIFY(dialog->grab().save(resultImages + "/result-expanded.png"));
        QVERIFY(std::abs(dialog->result().total[10] -
                         15 * (10 + .5 * dialog->result().energyTransferEv[10])) < 1e-5);
        auto acquisition = dialog->findChild<QComboBox *>("ComptonAcquisition");
        auto *peakSplit = dialog->findChild<QSplitter *>("ComptonPeakSplitter");
        auto loadWide = dialog->findChild<QPushButton *>("LoadComptonWideSpectrum");
        auto peakPlot = dialog->findChild<SpectrumPlot *>("ComptonPeakSpectrum");
        auto peakCrystal = dialog->findChild<QComboBox *>("ComptonPeakCrystal");
        auto peakCoverage = dialog->findChild<QLabel *>("ComptonPeakCoverage");
        auto spliceInfo = dialog->findChild<QLabel *>("ComptonSpliceInfo");
        auto limitPeak = dialog->findChild<QCheckBox *>("LimitComptonPeakRange");
        auto normalizeNear = dialog->findChild<QCheckBox *>("NormalizeComptonNearEdge");
        auto peakLo = dialog->findChild<QDoubleSpinBox *>("ComptonPeakMinimum");
        auto peakHi = dialog->findChild<QDoubleSpinBox *>("ComptonPeakMaximum");
        auto splice = dialog->findChild<QPushButton *>("SpliceComptonPeak");
        auto fullScan = dialog->findChild<QPushButton *>("ApplyComptonFullScan");
        auto fitBackground = dialog->findChild<QPushButton *>("RunMomentumBackground");
        auto backgroundInput = dialog->findChild<QLabel *>("BackgroundInputInfo");
        auto backgroundPlot = dialog->findChild<SpectrumPlot *>("MomentumBackgroundDiagnostics");
        auto backgroundSidebar = dialog->findChild<QWidget *>("BackgroundParameters");
        auto toggleBackground = dialog->findChild<QToolButton *>("ToggleBackgroundParameters");
        auto corrected = dialog->findChild<SpectrumPlot *>("CorrectedCrystalSpectrum");
        auto scope = dialog->findChild<QComboBox *>("BackgroundScope");
        auto diagnostic = dialog->findChild<QComboBox *>("BackgroundDiagnosticCrystal");
        auto previousPreview = dialog->findChild<QToolButton *>("PreviousBackgroundPreview");
        auto nextPreview = dialog->findChild<QToolButton *>("NextBackgroundPreview");
        auto batchList = dialog->findChild<QListWidget *>("BackgroundBatchTargets");
        auto batchModule = dialog->findChild<QPushButton *>("SelectBackgroundBatchModules");
        auto batchAll = dialog->findChild<QPushButton *>("SelectAllBackgroundTargets");
        auto batchInvert = dialog->findChild<QPushButton *>("InvertBackgroundTargets");
        auto fitPolynomial = dialog->findChild<QPushButton *>("FitBackground_polynomial");
        auto fitCompton = dialog->findChild<QPushButton *>("FitBackground_compton_peak");
        auto enableCompton = dialog->findChild<QCheckBox *>("EnableBackground_compton_peak");
        auto enableAsym = dialog->findChild<QCheckBox *>("EnableBackground_asymmetry");
        auto enableHf = dialog->findChild<QCheckBox *>("EnableBackground_hf_core");
        auto enablePolynomial = dialog->findChild<QCheckBox *>("EnableBackground_polynomial");
        auto polynomialOrder = dialog->findChild<QSpinBox *>("BackgroundPolynomialOrder");
        auto measureCompton = dialog->findChild<QCheckBox *>("MeasureCompton");
        auto protectedWindows = dialog->findChild<QLineEdit *>("BackgroundProtectedWindows");
        auto backgroundResults = dialog->findChild<QComboBox *>("BackgroundResultView");
        auto *nearChoice = dialog->findChild<QComboBox *>("BackgroundResultSpectrum");
        auto *nearPrevious = dialog->findChild<QToolButton *>("PreviousNearEdgeSpectrum");
        auto *nearNext = dialog->findChild<QToolButton *>("NextNearEdgeSpectrum");
        auto *smoothSpectrum = dialog->findChild<QCheckBox *>("CompareSmoothedSpectrum");
        auto *smoothingWindow = dialog->findChild<QSpinBox *>("SpectrumSmoothingWindow");
        auto *smoothingOrder = dialog->findChild<QSpinBox *>("SpectrumSmoothingOrder");
        auto comptonWindows = dialog->findChild<QLineEdit *>("BackgroundWindows_compton_peak");
        QVERIFY(acquisition && peakSplit && loadWide && peakPlot && peakCrystal && peakCoverage && spliceInfo &&
                normalizeNear && limitPeak && peakLo && peakHi && splice && fitBackground && corrected &&
                backgroundInput && backgroundPlot && scope && diagnostic && fitPolynomial && fitCompton &&
                enableCompton && enableAsym && enableHf && enablePolynomial && measureCompton && protectedWindows &&
                backgroundResults && comptonWindows);
        QVERIFY(previousPreview && nextPreview && batchList && batchModule && batchAll && batchInvert);
        QVERIFY(nearChoice && nearPrevious && nearNext && fullScan);
        QVERIFY(polynomialOrder && smoothSpectrum && smoothingWindow && smoothingOrder);
        QVERIFY(backgroundSidebar && toggleBackground);
        QCOMPARE(fitBackground->parentWidget(), next->parentWidget());
        QVERIFY(!fitBackground->isVisible());
        QVERIFY(!dialog->findChild<QPushButton *>("FitComptonPeak"));
        QVERIFY(!dialog->findChild<QComboBox *>("BackgroundSpectrumSource"));
        auto finishBackgroundManually = [&] {
            QVERIFY(!next->isEnabled() && !tabs->isTabEnabled(7) && fitBackground->isEnabled());
            fitBackground->click();
            QTRY_VERIFY_WITH_TIMEOUT(!dialog->isBusy(), 10000);
            QCOMPARE(tabs->currentIndex(), 6);
            QVERIFY(next->isEnabled() && tabs->isTabEnabled(7));
            next->click();
        };

        QVERIFY(!measureCompton->isChecked() && !tabs->isTabEnabled(5));
        next->click();
        QCOMPARE(tabs->currentIndex(), 6);
        QVERIFY(backgroundInput->text().contains(QStringLiteral("第 5 步原始谱")));
        QVERIFY(backgroundSidebar->isVisible() && fitBackground->isVisible());
        QTRY_VERIFY(fitBackground->geometry().right() < next->geometry().left());
        QCOMPARE(fitBackground->geometry().center().y(), next->geometry().center().y());
        const int expandedPlotWidth = backgroundPlot->width();
        toggleBackground->click();
        QTRY_VERIFY(!backgroundSidebar->isVisible() && backgroundPlot->width() > expandedPlotWidth);
        QVERIFY(toggleBackground->isVisible() && fitBackground->isVisible());
        toggleBackground->click();
        QTRY_VERIFY(backgroundSidebar->isVisible() && fitBackground->isVisible());
        QVERIFY(!enableHf->isEnabled() && !enableAsym->isEnabled() && !enableCompton->isEnabled());
        enablePolynomial->setChecked(false);
        finishBackgroundManually();
        QCOMPARE(tabs->currentIndex(), 7);
        QVERIFY(!fitBackground->isVisible());
        QCOMPARE(nearChoice->currentText(), QStringLiteral("总谱与各模组"));
        QVERIFY(!nearPrevious->isEnabled() && nearNext->isEnabled() && !backgroundResults->isVisible());
        QCOMPARE(corrected->curves().size(), 3);
        QCOMPARE(corrected->curves().first().name, QStringLiteral("Total"));
        QVERIFY(std::abs(corrected->curves().first().y[10] - dialog->result().total[10]) < 1e-8);
        QVERIFY(!smoothingWindow->isEnabled() && !smoothingOrder->isEnabled());
        const auto total = corrected->curves().first();
        smoothSpectrum->setChecked(true);
        QVERIFY(smoothingWindow->isEnabled() && smoothingOrder->isEnabled());
        QCOMPARE(corrected->curves().size(), 6);
        QCOMPARE(corrected->curves().first().y, total.y);
        QCOMPARE(corrected->curves()[3].x, total.x);
        QCOMPARE(corrected->curves()[3].name, QStringLiteral("Total · 平滑"));
        const auto initiallySmoothed = corrected->curves()[3].y;
        smoothingWindow->setValue(10);
        QCOMPARE(smoothingWindow->value(), 11);
        QVERIFY(corrected->curves()[3].y != initiallySmoothed);
        smoothSpectrum->setChecked(false);
        QCOMPARE(corrected->curves().size(), 3);
        nearNext->click();
        QCOMPARE(nearChoice->currentText(), QStringLiteral("各模组谱"));
        QCOMPARE(corrected->curves().size(), 2);
        for (const auto &module : corrected->curves())
            QVERIFY(std::abs(module.y[10] - dialog->result().modules.value(module.name)[10]) < 1e-8);
        nearNext->click();
        QCOMPARE(nearChoice->currentText(), QStringLiteral("全部单晶体谱"));
        QCOMPARE(corrected->curves().size(), 3);
        QCOMPARE(corrected->curves()[0].name, QStringLiteral("VB-A1"));
        QCOMPARE(corrected->curves()[0].y, dialog->result().crystals[0]);
        nearPrevious->click();
        QCOMPARE(nearChoice->currentText(), QStringLiteral("各模组谱"));
        nearChoice->setCurrentIndex(nearChoice->findText(QStringLiteral("VB 模组谱")));
        QCOMPARE(corrected->curves().size(), 1);
        QCOMPARE(corrected->curves().first().y, dialog->result().modules.value("VB"));
        QVERIFY(backgroundResults->isVisible());
        nearChoice->setCurrentIndex(nearChoice->findText(QStringLiteral("全部晶体等权平均")));
        QCOMPARE(corrected->curves().size(), 1);
        QVERIFY(std::abs(corrected->curves().first().y[10] * 3. - dialog->result().total[10]) < 1e-8);

        dialog->selectStep(6);
        scope->setCurrentIndex(1);
        finishBackgroundManually();
        QCOMPARE(tabs->currentIndex(), 7);
        QCOMPARE(nearChoice->currentText(), QStringLiteral("总谱与各模组"));
        QCOMPARE(corrected->curves().size(), 3);
        QCOMPARE(nearChoice->findText(QStringLiteral("全部单晶体谱")), -1);

        dialog->selectStep(6);
        scope->setCurrentIndex(2);
        QCOMPARE(backgroundPlot->curves().size(), 1);
        const auto raw = backgroundPlot->curves().first();
        protectedWindows->setText("unused");
        finishBackgroundManually();
        QCOMPARE(tabs->currentIndex(), 7);
        QCOMPARE(nearChoice->count(), 1);
        QVERIFY(!nearPrevious->isEnabled() && !nearNext->isEnabled());
        QCOMPARE(corrected->curves().size(), 1);
        QCOMPARE(corrected->curves().first().name, QStringLiteral("原谱"));
        QCOMPARE(corrected->curves().first().x, raw.x);
        QCOMPARE(corrected->curves().first().y, raw.y);
        auto previous = dialog->findChild<QPushButton *>("PreviousSpectrumStep");
        QVERIFY(previous);
        previous->click();
        QCOMPARE(tabs->currentIndex(), 6);
        previous->click();
        QCOMPARE(tabs->currentIndex(), 4);

        measureCompton->setChecked(true);
        QCOMPARE(tabs->isTabEnabled(5), true);
        scope->setCurrentIndex(0);
        QVERIFY(!tabs->isTabEnabled(6));
        QVERIFY(enableHf->isChecked() && enableAsym->isChecked() && enableCompton->isChecked());
        enableHf->setChecked(false);
        enableAsym->setChecked(false);
        enableCompton->setChecked(false);
        enablePolynomial->setChecked(true);
        protectedWindows->clear();

        next->click();
        QCOMPARE(tabs->currentIndex(), 5);
        QCOMPARE(tabs->tabText(5), QStringLiteral("6 Compton 峰分析"));
        QCOMPARE(acquisition->itemText(0), QStringLiteral("Full Scan"));
        QCOMPARE(acquisition->itemText(1), QStringLiteral("独立测量"));
        QCOMPARE(peakSplit->orientation(), Qt::Horizontal);
        QVERIFY(peakSplit->widget(0)->isAncestorOf(acquisition) && peakSplit->widget(1)->isAncestorOf(peakPlot));
        QVERIFY(!normalizeNear->isVisible() && !normalizeNear->isChecked());
        QVERIFY(!loadWide->isVisible() && !splice->isVisible() && peakPlot->curves().size() == 1);
        QCOMPARE(peakPlot->curves().first().name, QStringLiteral("full_scan"));
        QCOMPARE(peakPlot->curves().first().x, dialog->result().energyTransferEv);
        QCOMPARE(peakCrystal->currentText(), QStringLiteral("总谱"));
        QVERIFY(std::abs(peakPlot->curves().first().y[10] - dialog->result().total[10]) < 1e-8);
        const int vb = peakCrystal->findData("module:VB");
        QVERIFY(vb >= 0);
        peakCrystal->setCurrentIndex(vb);
        QVERIFY(std::abs(peakPlot->curves().first().y[10] - dialog->result().modules.value("VB")[10]) < 1e-8);
        const int crystal = peakCrystal->findData("crystal:VB-A1");
        QVERIFY(crystal >= 0);
        peakCrystal->setCurrentIndex(crystal);
        dialog->findChild<QLineEdit *>("BackgroundWindows_polynomial")->setText("0:1.5,8.5:10");
        QVERIFY(!next->isEnabled() && !tabs->isTabEnabled(6) && fullScan->isEnabled());
        fullScan->click();
        QVERIFY(next->isEnabled());
        next->click();
        QCOMPARE(tabs->currentIndex(), 6);
        QCOMPARE(tabs->tabText(6), QStringLiteral("7 背景分析"));
        QVERIFY(!fitBackground->isEnabled() && fitPolynomial->isEnabled() && corrected->curves().isEmpty());
        QVERIFY(enableHf->isEnabled() && enableAsym->isEnabled() && enableCompton->isEnabled());
        QVERIFY(backgroundInput->text().contains("Full Scan"));
        QCOMPARE(backgroundPlot->curves().size(), 1);
        QCOMPARE(backgroundPlot->curves().first().x, dialog->result().energyTransferEv);
        QCOMPARE(diagnostic->currentIndex(), 0);
        QVERIFY(!previousPreview->isEnabled() && nextPreview->isEnabled());
        const auto firstPreview = backgroundPlot->curves().first().y;
        nextPreview->click();
        QCOMPARE(diagnostic->currentIndex(), 1);
        QVERIFY(previousPreview->isEnabled() && backgroundPlot->curves().first().y != firstPreview);
        previousPreview->click();
        QCOMPARE(backgroundPlot->curves().first().y, firstPreview);
        diagnostic->setCurrentIndex(diagnostic->count() - 1);
        QVERIFY(previousPreview->isEnabled() && !nextPreview->isEnabled());
        diagnostic->setCurrentIndex(0);
        QVERIFY(batchModule->isEnabled() && batchAll->isEnabled() && batchInvert->isEnabled());
        batchInvert->click();
        for (int i = 0; i < batchList->count(); ++i)
            QCOMPARE(batchList->item(i)->checkState(), Qt::Unchecked);
        batchAll->click();
        for (int i = 0; i < batchList->count(); ++i)
            QCOMPARE(batchList->item(i)->checkState(), Qt::Checked);
        bool pickedVb = false;
        QTimer::singleShot(0, dialog, [&] {
            if (auto *picker = qobject_cast<QDialog *>(QApplication::activeModalWidget())) {
                for (auto *module : picker->findChildren<QCheckBox *>())
                    if (module->property("module").toString() == "VB") {
                        module->setChecked(true);
                        pickedVb = true;
                    }
                picker->done(pickedVb ? QDialog::Accepted : QDialog::Rejected);
            }
        });
        batchModule->click();
        QVERIFY(pickedVb);
        nextPreview->click();
        for (int i = 0; i < batchList->count(); ++i) {
            auto *target = batchList->item(i);
            QCOMPARE(target->checkState(), target->data(Qt::UserRole).toString().startsWith("VB-")
                                              ? Qt::Checked : Qt::Unchecked);
        }
        if (!resultImages.isEmpty()) {
            for (auto *toggle : tabs->widget(6)->findChildren<QToolButton *>())
                if (toggle->text() == QStringLiteral("批量目标"))
                    toggle->click();
            QVERIFY(batchList->isVisible());
            QVERIFY(dialog->grab().save(resultImages + "/background-targets.png"));
        }
        batchInvert->click();
        for (int i = 0; i < batchList->count(); ++i) {
            auto *target = batchList->item(i);
            QCOMPARE(target->checkState(), target->data(Qt::UserRole).toString().startsWith("VB-")
                                              ? Qt::Unchecked : Qt::Checked);
        }
        dialog->findChild<QLineEdit *>("BackgroundWindows_polynomial")->setText("0:1,9:10");
        scope->setCurrentIndex(1);
        QVERIFY(!enableHf->isEnabled() && !enableAsym->isEnabled() && enableCompton->isEnabled());
        diagnostic->setCurrentIndex(diagnostic->findData("VB"));
        QCOMPARE(backgroundPlot->curves().first().y[10], dialog->result().modules.value("VB")[10]);
        scope->setCurrentIndex(2);
        QCOMPARE(backgroundPlot->curves().first().y[10], dialog->result().total[10]);
        QVERIFY(!previousPreview->isEnabled() && !nextPreview->isEnabled() && !batchModule->isEnabled());
        QVERIFY(!enableHf->isEnabled() && !enableAsym->isEnabled() && enableCompton->isEnabled());
        fitPolynomial->click();
        QTRY_VERIFY_WITH_TIMEOUT(fitBackground->isEnabled(), 10000);
        QCOMPARE(backgroundPlot->curves().size(), 2);
        polynomialOrder->setValue(2);
        QVERIFY(!fitBackground->isEnabled() && fitPolynomial->isEnabled() && !next->isEnabled());
        QCOMPARE(backgroundPlot->curves().size(), 1);
        fitPolynomial->click();
        QTRY_VERIFY_WITH_TIMEOUT(fitBackground->isEnabled(), 10000);
        finishBackgroundManually();
        QCOMPARE(tabs->currentIndex(), 7);
        dialog->findChild<QDoubleSpinBox *>("CorrectedEnergyMinimum")->setValue(2);
        dialog->findChild<QDoubleSpinBox *>("CorrectedEnergyMaximum")->setValue(4);
        QCOMPARE(corrected->curves().first().x.first(), 2.);
        QCOMPARE(corrected->curves().first().x.last(), 4.);
        dialog->selectStep(6);
        polynomialOrder->setValue(1);
        scope->setCurrentIndex(0);
        QVERIFY(enableHf->isEnabled() && enableAsym->isEnabled() && !fitBackground->isEnabled());
        QVERIFY(corrected->curves().isEmpty());
        dialog->selectStep(5);
        ProcessedSpectrum wide;
        wide.intensityUnit = "Intensity / |I0|";
        const auto label = peakCrystal->currentText();
        const auto nearReference = peakPlot->curves().first();
        const auto slope = .5 * (nearReference.y.last() - nearReference.y.first()) /
                           (nearReference.x.last() - nearReference.x.first());
        auto &wideSamples = wide.crystals[label];
        for (int energy = 100; energy <= 300; energy += 20) {
            wide.energyEv.append(energy);
            wideSamples.append(.5 * nearReference.y.first() + slope * (energy - nearReference.x.first()) +
                               50. * std::pow(1. + .41421356237 * std::pow((energy - 200.) / 30., 2.), -2.));
        }
        const auto widePath = temp.filePath("wide.h5");
        SpectrumProcessing::save(widePath, wide);
        acquisition->setCurrentIndex(1);
        QVERIFY(loadWide->isVisible() && splice->isVisible() && !splice->isEnabled() && peakPlot->curves().isEmpty());
        QVERIFY(normalizeNear->isVisible() && normalizeNear->isEnabled());
        dialog->setAvailableFiles({});
        chooseWorkflowFile(loadWide, widePath);
        QTRY_VERIFY_WITH_TIMEOUT(!dialog->isBusy(), 10000);
        QCOMPARE(peakPlot->curves().size(), 2);
        QCOMPARE(peakPlot->curves()[0].x, wide.energyEv);
        QCOMPARE(peakPlot->curves()[1].x, dialog->result().energyTransferEv);
        QCOMPARE(peakPlot->curves()[1].draw, SpectrumCurve::Draw::Circles);
        QVERIFY(peakCoverage->text().contains(QStringLiteral("两段能区相离")));
        QVERIFY(peakCoverage->text().contains(QStringLiteral("近边谱由 2 次扫描叠加")));
        const auto summedNear = peakPlot->curves()[1].y;
        normalizeNear->click();
        QCOMPARE(peakPlot->curves()[0].y, wideSamples);
        QCOMPARE(peakPlot->curves()[1].x, dialog->result().energyTransferEv);
        for (qsizetype i = 0; i < summedNear.size(); ++i)
            QCOMPARE(peakPlot->curves()[1].y[i], summedNear[i] / 2.);
        QVERIFY(peakCoverage->text().contains(QStringLiteral("近边谱已按 2 次扫描取平均")));
        if (!resultImages.isEmpty())
            QVERIFY(dialog->grab().save(resultImages + "/compton-near-edge-average.png"));
        normalizeNear->click();
        QCOMPARE(peakPlot->curves()[1].y, summedNear);
        limitPeak->setChecked(true);
        peakLo->setValue(120);
        peakHi->setValue(280);
        const auto selectedEnergy = wide.energyEv.mid(1, 9);
        QCOMPARE(peakPlot->curves().size(), 1);
        QCOMPARE(peakPlot->curves().first().x, selectedEnergy);
        QVERIFY(peakCoverage->text().contains(QStringLiteral("9 个点")));
        QVERIFY(splice->isEnabled());
        splice->click();
        QTRY_VERIFY_WITH_TIMEOUT(!dialog->isBusy(), 10000);
        QCOMPARE(tabs->currentIndex(), 5);
        QVERIFY(spliceInfo->isVisible() && !splice->isEnabled() && peakPlot->curves().size() == 1);
        limitPeak->setChecked(false);
        QCOMPARE(peakPlot->curves().first().x, wide.energyEv);
        if (!resultImages.isEmpty())
            QVERIFY(dialog->grab().save(resultImages + "/compton-splice.png"));
        normalizeNear->click();
        QVERIFY(!spliceInfo->isVisible() && splice->isEnabled() && !next->isEnabled());
        splice->click();
        QTRY_VERIFY_WITH_TIMEOUT(!dialog->isBusy(), 10000);
        QCOMPARE(tabs->currentIndex(), 5);
        QVERIFY(next->isEnabled());
        next->click();
        QTRY_COMPARE_WITH_TIMEOUT(tabs->currentIndex(), 6, 10000);
        QVERIFY(backgroundInput->text().contains(QStringLiteral("第 6 步 Compton 拼接谱")));
        QCOMPARE(backgroundPlot->curves().size(), 1);
        QCOMPARE(backgroundPlot->curves().first().draw, SpectrumCurve::Draw::Circles);
        auto expectedEnergy = dialog->result().energyTransferEv;
        expectedEnergy += wide.energyEv;
        QCOMPARE(backgroundPlot->curves().first().x, expectedEnergy);
        QVector<double> expectedIntensity;
        for (const auto value : summedNear)
            expectedIntensity.append(value / 2.);
        expectedIntensity += wideSamples;
        QCOMPARE(backgroundPlot->curves().first().y, expectedIntensity);
        if (!resultImages.isEmpty())
            QVERIFY(dialog->grab().save(resultImages + "/background-analysis-input.png"));
        enableCompton->setChecked(true);
        comptonWindows->setText("100:300");
        QVERIFY(!fitCompton->isEnabled());
        fitPolynomial->click();
        QTRY_VERIFY_WITH_TIMEOUT(fitCompton->isEnabled(), 10000);
        QCOMPARE(backgroundPlot->curves().size(), 2);
        fitCompton->click();
        QTRY_VERIFY_WITH_TIMEOUT(fitBackground->isEnabled(), 10000);
        QCOMPARE(backgroundPlot->curves().size(), 4);
        QCOMPARE(backgroundPlot->curves()[1].name, QStringLiteral("多项式背景"));
        QCOMPARE(backgroundPlot->curves()[2].name, QStringLiteral("Compton 峰拟合背景"));
        QCOMPARE(backgroundPlot->curves().last().name, QStringLiteral("总背景"));
        const auto peakIndex = expectedEnergy.indexOf(200.);
        QVERIFY(peakIndex >= 0);
        QVERIFY(std::abs(backgroundPlot->curves().first().y[peakIndex] -
                         backgroundPlot->curves().last().y[peakIndex]) < 1e-6);
        comptonWindows->setText("100:299");
        QVERIFY(!fitBackground->isEnabled() && fitCompton->isEnabled());
        QCOMPARE(backgroundPlot->curves().size(), 2);
        fitCompton->click();
        QTRY_VERIFY_WITH_TIMEOUT(fitBackground->isEnabled(), 10000);
        if (!resultImages.isEmpty())
            QVERIFY(dialog->grab().save(resultImages + "/staged-background.png"));
        fitBackground->click();
        QTRY_VERIFY_WITH_TIMEOUT(!corrected->curves().isEmpty(), 10000);
        QCOMPARE(backgroundPlot->curves().size(), 4);
        QCOMPARE(backgroundPlot->curves().last().name, QStringLiteral("总背景"));
        dialog->selectStep(7);
        backgroundResults->setCurrentIndex(1);
        QCOMPARE(corrected->curves().size(), 5);
        QCOMPARE(corrected->curves()[3].name, QStringLiteral("总背景"));
        QCOMPARE(corrected->curves().last().name, QStringLiteral("扣背底结果"));
        for (qsizetype i = 0; i < corrected->curves().last().y.size(); ++i)
            QVERIFY(std::abs(corrected->curves().last().y[i] -
                             (corrected->curves().first().y[i] - corrected->curves()[3].y[i])) < 1e-9);
        QCOMPARE(tabs->tabText(7), QStringLiteral("8 近边结构结果分析"));
        auto resultColumns = dialog->findChild<QSplitter *>("NearEdgeResultSplitter");
        auto saveNearEdge = dialog->findChild<QPushButton *>("SaveNearEdgeHdf5");
        auto saveNearEdgeText = dialog->findChild<QPushButton *>("SaveNearEdgeText");
        QVERIFY(resultColumns && saveNearEdge && saveNearEdgeText);
        QCOMPARE(resultColumns->orientation(), Qt::Horizontal);
        QVERIFY(resultColumns->widget(0)->isAncestorOf(saveNearEdgeText));
        QVERIFY(resultColumns->widget(1)->isAncestorOf(corrected));
        QVERIFY(saveNearEdge->isEnabled() && saveNearEdgeText->isEnabled());
        const auto unsmoothedResult = corrected->curves().last();
        smoothSpectrum->setChecked(true);
        QCOMPARE(corrected->curves().size(), 6);
        QCOMPARE(corrected->curves().last().name, QStringLiteral("扣背底结果 · 平滑"));
        QCOMPARE(corrected->curves()[4].y, unsmoothedResult.y);
        if (!resultImages.isEmpty())
            QVERIFY(dialog->grab().save(resultImages + "/near-edge-result.png"));
        {
            QTemporaryDir workspace;
            QVERIFY(workspace.isValid());
            QSettings settings;
            const QString key = "workspace/root";
            const auto previous = settings.value(key);
            const auto restore = qScopeGuard([&] {
                if (previous.isValid())
                    settings.setValue(key, previous);
                else
                    settings.remove(key);
                settings.sync();
            });
            WorkspacePaths::setRoot(workspace.path());
            saveNearEdge->click();
            QTRY_VERIFY_WITH_TIMEOUT(!dialog->isBusy(), 10000);
            saveNearEdgeText->click();
            QTRY_VERIFY_WITH_TIMEOUT(!dialog->isBusy(), 10000);
            const QDir exports(workspace.filePath("spectroscopy/results"));
            const auto textFiles = exports.entryList({"*.txt"}, QDir::Files);
            const auto hdfFiles = exports.entryList({"*.h5"}, QDir::Files);
            QCOMPARE(textFiles.size(), 1);
            QCOMPARE(hdfFiles.size(), 1);
            QVERIFY(textFiles.first().contains("_O_K_near_edge_"));
            QVERIFY(hdfFiles.first().contains("_O_K_near_edge_"));
            QFile textResult(exports.filePath(textFiles.first()));
            QVERIFY(textResult.open(QIODevice::ReadOnly));
            const auto header = textResult.readLine();
            QVERIFY(header.contains("VB-A1_background_compton_peak") && header.contains("VB-A1_background_polynomial"));
            const auto firstRow = textResult.readLine().trimmed().split('\t');
            SpectrumHdfLayout savedLayout;
            savedLayout.intensity = "/entry/spectrum/corrected";
            const auto saved = SpectrumProcessing::load(exports.filePath(hdfFiles.first()), savedLayout);
            QCOMPARE(saved.energyEv, unsmoothedResult.x);
            QCOMPARE(firstRow.first().toDouble(), saved.energyEv.first());
            QVERIFY(saved.crystals.contains("VB-A1"));
            QCOMPARE(firstRow.last().toDouble(), saved.crystals.constFind("VB-A1")->first());
            QCOMPARE(saved.crystals.value("VB-A1"), unsmoothedResult.y);
            QCOMPARE(edgeGroups->item(0, 2)->text(), QStringLiteral("已完成"));
            const auto groupedPath = temp.filePath("processed.config.json");
            QVERIFY(dialog->saveConfiguration(groupedPath));
            QFile plan(groupedPath);
            QVERIFY(plan.open(QIODevice::ReadOnly));
            const auto edge = QJsonDocument::fromJson(plan.readAll()).object().value("groups").toArray().first().toObject();
            QCOMPARE(edge.value("status").toString(), QStringLiteral("completed"));
            QCOMPARE(edge.value("wide_file").toString(), widePath);
            const auto outputs = edge.value("last_saved_results").toObject();
            QCOMPARE(outputs.value("near_edge_text").toString(), exports.filePath(textFiles.first()));
            QCOMPARE(outputs.value("near_edge_hdf5").toString(), exports.filePath(hdfFiles.first()));
            SpectroscopyDialog reopened;
            QVERIFY(reopened.importConfiguration(groupedPath));
            QVERIFY(reopened.findChild<QCheckBox *>("CompareSmoothedSpectrum")->isChecked());
            QCOMPARE(reopened.findChild<QSpinBox *>("SpectrumSmoothingWindow")->value(), 11);
            QCOMPARE(reopened.findChild<QTableWidget *>("SpectroscopyEdgeGroups")->item(0, 2)->text(),
                     QStringLiteral("已完成"));
            for (int step : {0, 1, 2}) {
                reopened.runStep(step);
                QTRY_VERIFY_WITH_TIMEOUT(!reopened.isBusy(), 10000);
            }
            QVERIFY(reopened.hasResult());
            reopened.selectStep(5);
            QTRY_VERIFY_WITH_TIMEOUT(!reopened.isBusy(), 10000);
            const auto *restoredWide = reopened.findChild<SpectrumPlot *>("ComptonPeakSpectrum");
            QVERIFY(restoredWide && !restoredWide->curves().isEmpty());
            QCOMPARE(restoredWide->curves().first().x, wide.energyEv);
        }
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
        QCOMPARE(edgeGroups->item(0, 2)->text(), QStringLiteral("未处理"));
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
    void rawComptonWideScan() {
        QTemporaryDir temp;
        const auto external = qEnvironmentVariable("XRI_COMPTON_WIDE_SCAN");
        const auto elastic = external.isEmpty() ? SpectrumFixture::create(temp.path(), true)
                                                 : qEnvironmentVariable("XRI_COMPTON_ELASTIC");
        const auto wide = external.isEmpty() ? SpectrumFixture::create(temp.path(), false, true)
                                               : external;
        const auto rois = external.isEmpty() ? temp.filePath("test.roi.json")
                                               : qEnvironmentVariable("XRI_COMPTON_ROI");
        MainWindow window;
        window.show();
        window.openFiles({wide, elastic});
        auto dialog = window.openSpectroscopy();
        QVERIFY(dialog->loadRois(rois));
        QVERIFY(dialog->addEdgeGroup("C", "K"));
        dialog->findChild<QCheckBox *>("MeasureCompton")->setChecked(true);
        dialog->runStep(0);
        QTRY_VERIFY_WITH_TIMEOUT(!dialog->isBusy(), 30000);
        QVERIFY2(dialog->hasCalibration(), qPrintable(dialog->statusText()));
        const auto nearPath = qEnvironmentVariable("XRI_COMPTON_NEAR_SCAN");
        if (!nearPath.isEmpty()) {
            dialog->addScans({nearPath});
            for (int step : {1, 2}) {
                dialog->runStep(step);
                QTRY_VERIFY_WITH_TIMEOUT(!dialog->isBusy(), 30000);
            }
            QVERIFY2(dialog->hasResult(), qPrintable(dialog->statusText()));
        }
        dialog->selectStep(5);
        auto acquisition = dialog->findChild<QComboBox *>("ComptonAcquisition");
        auto loadWide = dialog->findChild<QPushButton *>("LoadComptonWideSpectrum");
        auto plot = dialog->findChild<SpectrumPlot *>("ComptonPeakSpectrum");
        auto status = dialog->findChild<QLabel *>("ComptonPeakStatus");
        QVERIFY(acquisition && loadWide && plot && status);
        acquisition->setCurrentIndex(1);

        bool pickedFromBrowser = false;
        const auto canonicalWide = QFileInfo(wide).canonicalFilePath();
        QTimer choose;
        connect(&choose, &QTimer::timeout, &window, [&] {
            auto *menu = qobject_cast<QMenu *>(QApplication::activePopupWidget());
            if (!menu)
                return;
            for (auto action : menu->actions()) {
                if (!canonicalWide.isEmpty() && QFileInfo(action->data().toString()).canonicalFilePath() == canonicalWide) {
                    pickedFromBrowser = true;
                    menu->setActiveAction(action);
                    QTest::keyClick(menu, Qt::Key_Return);
                    return;
                }
            }
            menu->close();
        });
        choose.start(30);
        loadWide->click();
        choose.stop();
        QVERIFY(pickedFromBrowser);
        QTRY_VERIFY_WITH_TIMEOUT(!dialog->isBusy(), 30000);
        QVERIFY2(!plot->curves().isEmpty(), qPrintable(status->text()));
        auto crystal = dialog->findChild<QComboBox *>("ComptonPeakCrystal");
        QVERIFY(crystal && crystal->findData("module:VB") >= 0);
        QCOMPARE(crystal->currentData().toString(), QStringLiteral("total"));
        const int firstCrystal = crystal->findData("crystal:VB-A1");
        QVERIFY(firstCrystal >= 0);
        crystal->setCurrentIndex(firstCrystal);
        auto rawAxis = HdfFile::read({wide, "/entry/instrument/M_DCM_B5Energy", {}, 0, -1}).values;
        if (rawAxis.isEmpty())
            rawAxis = HdfFile::read({wide, "/entry/instrument/M_DCM_B5Energy_readback", {}, 0, -1}).values;
        QVERIFY(!rawAxis.isEmpty());
        if (rawAxis.first() > rawAxis.last())
            std::reverse(rawAxis.begin(), rawAxis.end());
        const auto &curve = plot->curves().first();
        QCOMPARE(curve.name, QStringLiteral("full_scan"));
        QCOMPARE(curve.x.size(), rawAxis.size());
        QCOMPARE(curve.y.size(), rawAxis.size());
        for (qsizetype i = 0; i < rawAxis.size(); ++i)
            QVERIFY(std::abs((curve.x[i] - curve.x.first()) - 1000. * (rawAxis[i] - rawAxis.first())) < 1e-8);
        if (!nearPath.isEmpty()) {
            QCOMPARE(plot->curves().size(), 2);
            QCOMPARE(plot->curves().last().x, dialog->result().energyTransferEv);
            const auto native = curve;
            const auto near = plot->curves().last();
            auto *next = dialog->findChild<QPushButton *>("NextSpectrumStep");
            auto *splice = dialog->findChild<QPushButton *>("SpliceComptonPeak");
            auto *tabs = dialog->findChild<QTabWidget *>("SpectroscopySteps");
            auto *diagnostic = dialog->findChild<QComboBox *>("BackgroundDiagnosticCrystal");
            auto *input = dialog->findChild<SpectrumPlot *>("MomentumBackgroundDiagnostics");
            QVERIFY(next && splice && tabs && diagnostic && input);
            QVERIFY(!next->isEnabled());
            splice->click();
            QTRY_VERIFY_WITH_TIMEOUT(!dialog->isBusy(), 30000);
            QVERIFY(next->isEnabled());
            next->click();
            QTRY_COMPARE_WITH_TIMEOUT(tabs->currentIndex(), 6, 30000);
            QVERIFY(diagnostic->findText("VB-A1") >= 0);
            diagnostic->setCurrentText("VB-A1");
            QCOMPARE(input->curves().size(), 1);
            const auto &joined = input->curves().first();
            QVector<double> retained;
            for (qsizetype i = 0; i < joined.x.size(); ++i)
                if (std::isfinite(joined.y[i]))
                    retained.append(joined.x[i]);
            auto expected = near.x;
            for (qsizetype i = 0; i < native.x.size(); ++i)
                if (std::isfinite(native.y[i]) &&
                    (native.x[i] < near.x.first() - 1e-7 || native.x[i] > near.x.last() + 1e-7))
                    expected.append(native.x[i]);
            std::sort(expected.begin(), expected.end());
            QCOMPARE(retained, expected);
        }
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
