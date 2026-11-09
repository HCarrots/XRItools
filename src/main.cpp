#include "BondWorkspace.h"
#include "ImagingDialog.h"
#include "MainWindow.h"
#include "RoiDialog.h"
#include "SpectroscopyDialog.h"
#include <QApplication>
#include <QCheckBox>
#include <QCommandLineParser>
#include <QDir>
#include <QFileInfo>
#include <QIcon>
#include <QPointer>
#include <QStyleFactory>
#include <QTimer>
#include <cstdio>

int main(int argc, char *argv[]) {
    QApplication app(argc, argv);
    QCoreApplication::setOrganizationName("XRI");
    QCoreApplication::setApplicationName("XRI");
    QCoreApplication::setApplicationVersion("0.2.1");
    app.setWindowIcon(QIcon(QStringLiteral(":/assets/xri.png")));
    // Keep Qt Widgets' Windows platform style, native dialogs, menus and frame.
    if (auto style = QStyleFactory::create("windowsvista"))
        app.setStyle(style);
    QCommandLineParser parser;
    parser.setApplicationDescription("XRI - Synchrotron HDF5 / NeXus viewer");
    parser.addHelpOption();
    parser.addVersionOption();
    parser.addPositionalArgument("files", "HDF5 / NeXus files to import", "[files...]");
    parser.addOption({"dataset", "Select dataset path", "path"});
    parser.addOption({"frame", "Fixed dimension 0 index", "index", "0"});
    parser.addOption({"roi-editor", "Open the shared crystal ROI selection window"});
    parser.addOption({"roi-profiles", "Show X/Y line profiles in the ROI selection window"});
    parser.addOption({"spectroscopy", "Open the spectroscopy workbench"});
    parser.addOption({"imaging", "Open the imaging workbench"});
    parser.addOption({"imaging-config", "Load imaging JSON configuration", "path"});
    parser.addOption({"imaging-result", "View an existing imaging result directory", "path"});
    parser.addOption({"bond-result", "Open an existing four-dimensional contrast result", "path"});
    parser.addOption({"smoke-test", "Render current dataset and exit with verification status", "directory"});
    parser.process(app);
    MainWindow window;
    window.show();
    auto args = parser.positionalArguments();
    if (!args.isEmpty())
        window.openFiles(args);
    if (parser.isSet("dataset"))
        window.selectPath(parser.value("dataset"));
    window.setSliceIndex(0, parser.value("frame").toInt());
    QPointer<RoiDialog> roiEditor;
    QPointer<SpectroscopyDialog> spectroscopy;
    QPointer<ImagingDialog> imaging;
    QPointer<BondWorkspace> bond;
    bool startupOk = true;
    if (parser.isSet("imaging") || parser.isSet("imaging-config") || parser.isSet("imaging-result")) {
        imaging = window.openImaging();
        if (parser.isSet("imaging-config"))
            startupOk = imaging->importConfiguration(parser.value("imaging-config")) && startupOk;
        if (parser.isSet("imaging-result"))
            startupOk = imaging->loadResult(parser.value("imaging-result")) && startupOk;
    }
    if (parser.isSet("spectroscopy"))
        spectroscopy = window.openSpectroscopy();
    if (parser.isSet("bond-result")) {
        bond = window.openBondImaging();
        bond->openVolume(parser.value("bond-result"));
    }
    if (parser.isSet("roi-editor") && !args.isEmpty()) {
        roiEditor = window.openRoi();
        if (parser.isSet("roi-profiles"))
            roiEditor->findChild<QCheckBox *>("ShowRoiProfiles")->setChecked(true);
    }
    if (parser.isSet("smoke-test")) {
        auto directory = QFileInfo(parser.value("smoke-test")).absoluteFilePath();
        QDir().mkpath(directory);
        auto timer = new QTimer(&window);
        timer->setInterval(100);
        int *attempts = new int(0);
        QObject::connect(
            timer, &QTimer::timeout, &window,
            [&window, &app, timer, directory, attempts, roiEditor, spectroscopy, imaging, bond, startupOk] {
                if (++*attempts > 1200) {
                    timer->stop();
                    delete attempts;
                    app.exit(2);
                    return;
                }
                if (window.isReading() || (roiEditor && roiEditor->isLoading()) || (bond && bond->isBusy()))
                    return;
                timer->stop();
                delete attempts;
                auto data = window.currentData();
                bool ok = startupOk && (bond ? bond->hasResult()
                                             : ((imaging && imaging->hasResult()) ||
                                                (data.error.isEmpty() && !data.values.isEmpty())));
                if (roiEditor)
                    ok = roiEditor->imagesReady() && window.grab().save(directory + "/window.png") && ok;
                else if (bond) {
                    ok = window.grab().save(directory + "/voxel.png") && ok;
                    bond->selectStep(3);
                    ok = window.grab().save(directory + "/sections.png") && ok;
                } else if (imaging)
                    ok = window.grab().save(directory + "/window.png") && ok;
                else if (spectroscopy)
                    ok = window.grab().save(directory + "/window.png") && ok;
                else
                    ok = window.grab().save(directory + "/window.png") && ok;
                std::printf("%s: %d x %d, finite=%llu, min=%.12g, max=%.12g, mean=%.12g\n",
                            ok ? "PASS" : "FAIL", data.width, data.height,
                            static_cast<unsigned long long>(data.finiteCount), data.minimum, data.maximum,
                            data.mean);
                app.exit(ok ? 0 : 1);
            });
        timer->start();
    }
    return app.exec();
}
