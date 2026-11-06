#include "MainWindow.h"
#include <QApplication>
#include <QCommandLineParser>
#include <QDir>
#include <QFileInfo>
#include <QStyleFactory>
#include <QTimer>
#include <cstdio>

int main(int argc, char *argv[]) {
    QApplication app(argc, argv);
    QCoreApplication::setOrganizationName("XRI");
    QCoreApplication::setApplicationName("XRI");
    QCoreApplication::setApplicationVersion("0.1.0");
    // Keep Qt Widgets' Windows platform style, native dialogs, menus and frame.
    if (auto style = QStyleFactory::create("windowsvista"))
        app.setStyle(style);
    QCommandLineParser parser;
    parser.setApplicationDescription("XRI - Synchrotron HDF5 / NeXus viewer");
    parser.addHelpOption();
    parser.addVersionOption();
    parser.addPositionalArgument("file", "HDF5 / NeXus file to open");
    parser.addOption({"dataset", "Select dataset path", "path"});
    parser.addOption({"frame", "Fixed dimension 0 index", "index", "0"});
    parser.addOption({"smoke-test", "Render current dataset and exit with verification status", "directory"});
    parser.process(app);
    MainWindow window;
    window.show();
    auto args = parser.positionalArguments();
    if (!args.isEmpty())
        window.openFile(args.first());
    if (parser.isSet("dataset"))
        window.selectPath(parser.value("dataset"));
    window.setSliceIndex(0, parser.value("frame").toInt());
    if (parser.isSet("smoke-test")) {
        auto directory = QFileInfo(parser.value("smoke-test")).absoluteFilePath();
        QDir().mkpath(directory);
        auto timer = new QTimer(&window);
        timer->setInterval(100);
        int *attempts = new int(0);
        QObject::connect(timer, &QTimer::timeout, &window, [&window, &app, timer, directory, attempts] {
            if (++*attempts > 300) {
                timer->stop();
                delete attempts;
                app.exit(2);
                return;
            }
            if (window.isReading())
                return;
            timer->stop();
            delete attempts;
            auto data = window.currentData();
            bool ok = data.error.isEmpty() && !data.values.isEmpty();
            ok = window.grab().save(directory + "/window.png") && ok;
            std::printf("%s: %d x %d, finite=%llu, min=%.12g, max=%.12g, mean=%.12g\n", ok ? "PASS" : "FAIL",
                        data.width, data.height, static_cast<unsigned long long>(data.finiteCount),
                        data.minimum, data.maximum, data.mean);
            app.exit(ok ? 0 : 1);
        });
        timer->start();
    }
    return app.exec();
}
