#include "Spectroscopy.h"
#include "SpectrumFixture.h"
#include <QCoreApplication>
#include <QJsonDocument>
#include <QJsonObject>
#include <QTemporaryDir>
#include <iomanip>
#include <iostream>

namespace {
int failures = 0;
void expect(bool condition, const char *text) {
    std::cout << (condition ? "PASS " : "FAIL ") << text << '\n';
    if (!condition)
        ++failures;
}
bool close(double a, double b, double tol = 1e-7) {
    return std::isfinite(a) && std::isfinite(b) && std::abs(a - b) < tol;
}
template <class F> bool throws(F f) {
    try {
        f();
        return false;
    } catch (const std::exception &) {
        return true;
    }
}
} // namespace
int main(int argc, char **argv) {
    std::cout.setf(std::ios::unitbuf);
    std::cout << std::setprecision(17);
    QCoreApplication app(argc, argv);
    QTemporaryDir temp;
    std::atomic_bool cancel{false};
    SpectrumSettings s;
    try {
        if (argc == 4 && QString::fromLocal8Bit(argv[1]) == "--real-elastic") {
            std::cout << "Inspect real detectors\n";
            s.elasticFile = QString::fromLocal8Bit(argv[2]);
            QVector<QSize> sizes;
            for (const auto &name : {"D_LAMBDA", "D_MINIPIX"}) {
                auto node = HdfFile::inspect(s.elasticFile, QString("/entry/data/") + name);
                sizes.append({int(node.dims[2]), int(node.dims[1])});
            }
            QFile file(QString::fromLocal8Bit(argv[3]));
            expect(file.open(QIODevice::ReadOnly), "real ROI scheme opens read-only");
            CrystalRois realRois;
            QString error;
            expect(realRois.restore(file.readAll(), sizes, &error),
                   "real ROI scheme matches detector geometry");
            auto calibration =
                Spectroscopy::calibrate(s, realRois, cancel, [](const QString &phase, int progress) {
                    std::cout << progress << " " << phase.toStdString() << '\n';
                });
            int count = 0;
            for (int r = 0; r < 75; ++r)
                if (realRois.items()[r].state == RoiState::Valid) {
                    const auto &fit = calibration.fits[r];
                    count += fit.accepted;
                    std::cout << realRois.items()[r].label.toStdString() << " center=" << fit.centerKev
                              << " FWHM=" << fit.fwhmEv << " R2=" << fit.rSquared
                              << " accepted=" << fit.accepted << " reason=" << fit.reason.toUtf8().constData()
                              << '\n';
                }
            std::cout << "Accepted " << count << " / " << realRois.count(RoiState::Valid) << '\n';
            expect(calibration.images.size() == 2 && calibration.fits.size() == 75,
                   "real file calibration returns images and explicit fit states");
            return failures ? 1 : 0;
        }
        QVector<double> x, l, g;
        for (int i = 0; i < 201; ++i) {
            double e = 9.680 + i * .00005;
            x.append(e);
            l.append(123 / (1 + std::pow((e - 9.6848) * 1000 / .35, 2)) + 4);
            g.append(123 * std::exp(-.5 * std::pow((e - 9.6848) * 1000 / .35, 2)) + 4);
        }
        auto fit = Spectroscopy::fitPeak(x, l, s);
        expect(fit.accepted && close(fit.centerKev, 9.6848, 1e-9) && close(fit.fwhmEv, .7, 1e-6) &&
                   close(fit.rSquared, 1, 1e-10),
               "Lorentzian fit recovers known center and FWHM");
        s.fitModel = 1;
        fit = Spectroscopy::fitPeak(x, g, s);
        expect(fit.accepted && close(fit.centerKev, 9.6848, 1e-9) &&
                   close(fit.fwhmEv, .35 * 2.354820045, 1e-6),
               "Gaussian fit recovers known center and FWHM");
        expect(!Spectroscopy::fitPeak(x, QVector<double>(201, 1), s).accepted,
               "flat/empty elastic curves fail explicitly");
        s.fitModel = 0;
        s.maxFwhmEv = .3;
        expect(!Spectroscopy::fitPeak(x, l, s).accepted, "FWHM quality threshold excludes broad peak");
        s.maxFwhmEv = 2;
        auto elastic = SpectrumFixture::create(temp.path(), true);
        auto xrs = SpectrumFixture::create(temp.path(), false);
        auto reverse = SpectrumFixture::create(temp.path(), false, true);
        s.elasticFile = elastic;
        auto rois = SpectrumFixture::rois();
        auto c = Spectroscopy::calibrate(s, rois, cancel);
        expect(c.fits[0].accepted && c.fits[1].accepted && c.fits[60].accepted && !c.fits[2].accepted,
               "three valid crystals fit, manually invalid crystals skipped");
        expect(c.pixels[0].size() == 4 && c.filteredPixels[0].size() == 2,
               "threshold mask strictly selects pixels above 15 percent");
        expect(close(c.scan.filtered[0][40], 128.75, 1e-6) && close(c.scan.raw[0][40], 139.05, 1e-6),
               "elastic normalization and raw/masked integrals agree with independent values");
        auto scans = Spectroscopy::extract(s, c, {xrs, reverse}, cancel);
        expect(scans.size() == 2 && scans[0].glitches == QVector<int>{30} &&
                   scans[1].glitches == QVector<int>{50},
               "isolated I0 glitch detected in ascending and descending scans");
        expect(close(scans[0].i0Before[30], .2) && close(scans[0].i0After[30], 2) &&
                   close(scans[0].filtered[0][30], 1.25 * (10 + .5 * 3.75)),
               "glitch repairs both frame and I0 before normalization");
        auto result = Spectroscopy::combine(s, c, scans, cancel);
        double error = 0;
        for (int i = 0; i < result.total.size(); ++i)
            error = std::max(error, std::abs(result.total[i] - 15 * (10 + .5 * result.energyTransferEv[i])));
        expect(error < 1e-5 && result.total.size() == 51,
               "two opposite-direction scans interpolate and sum to analytic spectrum");
        expect(result.included[0] && result.included[1] && result.included[60] && !result.included[2],
               "invalid crystals excluded from result");
        auto q = Spectroscopy::momentum("VB-C2", 90, 9.685, 9.685, s);
        constexpr double pi = 3.14159265358979323846;
        double k = 2 * pi * 9.685 / 12.398419297617678;
        expect(close(q.first, k * std::sqrt(2.0)) && close(q.second, k / std::sqrt(2.0) * 100 / 1070),
               "Q and dQ match analytic 90 degree geometry");
        s.excludedCrystals.insert(1);
        auto selected = Spectroscopy::combine(s, c, scans, cancel);
        expect(!selected.included[1] &&
                   close(selected.total[10], 10 * (10 + .5 * selected.energyTransferEv[10]), 1e-5),
               "manual exclusion removes crystal contribution");
        s.modules = {"VB"};
        selected = Spectroscopy::combine(s, c, scans, cancel);
        expect(selected.included[0] && !selected.included[60] && selected.modules.size() == 1,
               "module selection is applied");
        s.qMax = .01;
        expect(throws([&] { Spectroscopy::combine(s, c, scans, cancel); }),
               "empty Q selection fails without invented output");
        s = SpectrumSettings();
        s.elasticFile = elastic;
        auto overlapping = rois;
        overlapping.assign(1, {1, 1, 2, 2}, {8, 8});
        auto ov = Spectroscopy::calibrate(s, overlapping, cancel);
        expect(ov.overlapPixels == 2 && ov.filteredPixels[0].isEmpty() && ov.filteredPixels[1].isEmpty() &&
                   !ov.fits[0].accepted,
               "overlap pixels removed from every affected crystal");
        s.removeOverlap = false;
        ov = Spectroscopy::calibrate(s, overlapping, cancel);
        expect(ov.filteredPixels[0].size() == 2 && ov.filteredPixels[1].size() == 2,
               "overlap removal can be disabled explicitly");
        s.removeOverlap = true;
        s.autoRoi = true;
        s.roiSize = 2;
        s.shifts[0] = {1, 0};
        auto adjusted = Spectroscopy::calibrate(s, rois, cancel);
        expect(adjusted.original.items()[0].rect == QRect(1, 1, 2, 2) &&
                   adjusted.adjusted.items()[0].rect == QRect(1, 0, 2, 2),
               "auto recenter and manual shift preserve original ROI");
        s.autoRoi = false;
        s.shifts[0] = {};
        expect(throws([&] { Spectroscopy::extract(s, c, {xrs, xrs}, cancel); }), "duplicate scans rejected");
        auto grid = temp.path() + "/4_grid_scan_20261106000000_1.nxs";
        expect(QFile::copy(xrs, grid) && throws([&] { Spectroscopy::extract(s, c, {grid}, cancel); }),
               "grid_scan is not accepted as a spectroscopy scan");
        auto badEnergy = scans;
        badEnergy[1].energyKev.fill(20, badEnergy[1].energyKev.size());
        expect(throws([&] { Spectroscopy::combine(s, c, badEnergy, cancel); }),
               "non-overlapping energy-transfer ranges cannot be combined");
        auto noFluxDir = temp.path() + "/without-i0";
        expect(QDir().mkpath(noFluxDir), "create missing-I0 fixture directory");
        auto noFlux = SpectrumFixture::create(noFluxDir, false, false, false);
        expect(throws([&] { Spectroscopy::extract(s, c, {noFlux}, cancel); }),
               "missing required I0 produces an error");
        s.normalizeXrs = false;
        auto unnormalized = Spectroscopy::extract(s, c, {noFlux}, cancel);
        expect(unnormalized[0].glitches.isEmpty() && unnormalized[0].i0Path.isEmpty() &&
                   close(unnormalized[0].filtered[0][0], 25),
               "normalization off reads raw counts without requiring I0");
        s.normalizeXrs = true;
        auto pending = rois;
        pending.clear(2);
        expect(throws([&] { Spectroscopy::calibrate(s, pending, cancel); }), "pending ROI blocks processing");
        cancel.store(true);
        expect(throws([&] { Spectroscopy::calibrate(s, rois, cancel); }),
               "cancellation interrupts processing");
        cancel.store(false);
        QString errorText;
        auto output = temp.path() + "/result";
        expect(Spectroscopy::saveResult(output, s, c, scans, result, &errorText),
               "results and provenance saved as a new directory");
        QFile info(output + "/spectrum_info.json");
        expect(info.open(QIODevice::ReadOnly), "processing provenance can be reopened");
        auto object = QJsonDocument::fromJson(info.readAll()).object();
        expect(object.value("combination").toString().contains("no averaging") &&
                   QFileInfo::exists(output + "/VB-A1_data.txt") &&
                   QFileInfo::exists(output + "/adjusted_roi.roi.json") &&
                   QFileInfo::exists(output + "/spectrum_masks.json"),
               "saved outputs contain physical convention, individual curves, ROIs and masks");
        expect(!Spectroscopy::saveResult(output, s, c, scans, result, &errorText),
               "existing result directory never overwritten");
        if (argc > 1) {
            QDir().mkpath(QString::fromLocal8Bit(argv[1]));
            auto directory = QString::fromLocal8Bit(argv[1]);
            SpectrumFixture::create(directory, true);
            SpectrumFixture::create(directory, false);
            SpectrumFixture::create(directory, false, true);
        }
    } catch (const std::exception &e) {
        std::cerr << "UNEXPECTED " << e.what() << '\n';
        ++failures;
    }
    return failures ? 1 : 0;
}
