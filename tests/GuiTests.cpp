#include "MainWindow.h"
#include "ScientificPlot.h"
#include <QApplication>
#include <QDir>
#include <QElapsedTimer>
#include <QFileInfo>
#include <QSpinBox>
#include <QTabWidget>
#include <QtTest>
#include <cmath>

class GuiTests : public QObject {
    Q_OBJECT
  private slots:
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
    QCoreApplication::setOrganizationName("XRI-Tests");
    QCoreApplication::setApplicationName("XRI-Tests");
    GuiTests tests;
    return QTest::qExec(&tests, argc, argv);
}
#include "GuiTests.moc"
