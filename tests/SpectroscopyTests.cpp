#include "ComptonProcessing.h"
#include "ImagingIO.h"
#include "Spectroscopy.h"
#include "SpectrumFixture.h"
#include "SpectrumProcessing.h"
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
        if (argc == 3 && QString::fromLocal8Bit(argv[1]) == "--compton-reference") {
            QDir root(QString::fromLocal8Bit(argv[2]));
            auto load = [&](const QString &name) {
                return SpectrumProcessing::load(
                    root.filePath("workspace/data/" + name + "/" + name + "_all_data.txt"), {});
            };
            const auto wide = load("Compton_all_crystal");
            const auto stitched =
                SpectrumProcessing::replaceEdges(wide, {load("O_K_530_550"), load("Ni_L_840_870")}, cancel);
            qsizetype replacements = 0;
            for (const auto &origin : stitched.pointSources)
                replacements +=
                    std::count_if(origin.cbegin(), origin.cend(), [](int source) { return source > 0; });
            expect(replacements == 22788 && stitched.energyEv == wide.energyEv &&
                       stitched.crystals.keys() == wide.crystals.keys(),
                   "reference O K and Ni L scans replace the expected 22788 wide-scan samples");
            ComptonSettings settings;
            settings.hfTable = root.filePath("xrsana/resources/data/ComptonProfiles.dat");
            settings.composition = "Li:1, Ni:0.9, Co:0.05, Mn:0.05, O:2";
            settings.normalization = {{80, 970}};
            settings.coreTail = {{700, 740}, {930, 970}};
            settings.cleanReference = {{120, 510}, {580, 620}, {680, 740}, {920, 970}};
            settings.highResidual = settings.lowLinear = {{460, 510}};
            auto result = analyzeCompton(wide, stitched, settings, cancel);
            int high = 0, low = 0;
            for (auto it = result.crystals.cbegin(); it != result.crystals.cend(); ++it) {
                if (it->high)
                    ++high;
                else
                    ++low;
                expect(it->high == (wide.geometry.value(it.key()).qMean >= settings.highQMin),
                       "successful crystal uses the selected q branch");
            }
            std::cout << "High=" << high << " Low=" << low << " References=" << result.references.size()
                      << " PearsonR2=" << result.referenceFit.rSquared
                      << " PearsonCenter=" << result.referenceFit.centerEv
                      << " PearsonWidth=" << result.referenceFit.fwhmEv << '\n';
            for (auto it = result.failures.cbegin(); it != result.failures.cend(); ++it)
                std::cout << "Excluded " << it.key().toStdString() << " " << it.value().toUtf8().constData()
                          << '\n';
            expect(high > 0 && low > 0 && result.referenceFit.rSquared > .8,
                   "real reduced scans produce high-q valence subtraction and low-q linear results");
            const auto highCrystal = std::find_if(result.crystals.cbegin(), result.crystals.cend(),
                                                  [](const ComptonCrystal &item) { return item.high; });
            if (highCrystal != result.crystals.cend()) {
                const auto geometry = wide.geometry.value(highCrystal.key());
                QVector<double> incident;
                for (double loss : wide.energyEv)
                    incident.append(geometry.analyzerKev + loss / 1000.);
                const auto model = transferCompton(wide, highCrystal.key(), incident, settings, cancel);
                qsizetype compared = 0;
                double maxError = 0;
                for (qsizetype i = 0; i < wide.energyEv.size(); ++i)
                    if (std::isfinite(model.valence[i]) && std::isfinite(highCrystal->valence[i])) {
                        ++compared;
                        maxError = std::max(maxError, std::abs(model.valence[i] - highCrystal->valence[i]));
                    }
                expect(compared > 100 && maxError < 1e-10 && model.analyzerKev == geometry.analyzerKev,
                       "transferred valence matches spectroscopy on the original incident grid");
            }
            return failures ? 1 : 0;
        }
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
        auto saved = SpectrumProcessing::stitch({output + "/spectra.h5"}, {}, cancel);
        expect(saved.energyEv == result.energyTransferEv && saved.crystals["VB-A1"] == result.crystals[0],
               "step five saves an HDF5 spectrum that can be loaded for processing");
        auto first = temp.path() + "/part-a.h5", second = temp.path() + "/part-b.h5";
        SpectrumProcessing::save(first, {{0, 1, 2}, {{"A", {1, 2, 3}}, {"B", {10, 20, 30}}}, "counts", {}});
        {
            ImagingH5 hdf(second, true);
            const double axis[]{.003, .002, .001};
            const double values[]{40, 5, 50, 7, 60, 9};
            hdf.write("/axis", {3}, axis);
            hdf.writeStrings("/labels", {"B", "A"});
            hdf.write("/values", {3, 2}, values);
        }
        SpectrumHdfLayout external{"/axis", "/labels", "/values", "counts", true, true};
        auto converted = SpectrumProcessing::stitch({second}, external, cancel);
        auto convertedPath = temp.path() + "/converted.h5";
        SpectrumProcessing::save(convertedPath, converted);
        auto stitched = SpectrumProcessing::stitch({first, convertedPath, first}, {}, cancel);
        expect(stitched.energyEv == QVector<double>({0, 1, 2, 3}) &&
                   stitched.crystals["A"] == QVector<double>({1, 5.5, 5, 5}) &&
                   stitched.crystals["B"] == QVector<double>({10, 40, 40, 40}) &&
                   stitched.sources.size() == 2,
               "align crystal names and energy points, average overlap and ignore repeated input files");
        auto wide = SpectrumProcessing::load(first, {});
        ProcessedSpectrum edge{{1, 2}, {{"A", {50, 60}}, {"C", {9, 10}}}, "counts", {"edge"}};
        auto replaced = SpectrumProcessing::replaceEdges(wide, {edge}, cancel);
        expect(replaced.energyEv == wide.energyEv && replaced.crystals["A"] == QVector<double>({1, 50, 60}) &&
                   replaced.crystals["B"] == wide.crystals["B"] && !replaced.crystals.contains("C") &&
                   replaced.pointSources["A"] == QVector<int>({0, 1, 1}),
               "near-edge replacement preserves wide axes, unmeasured crystals and source identity");
        expect(throws([&] { SpectrumProcessing::replaceEdges(wide, {edge, edge}, cancel); }),
               "overlapping replacements require a source decision");
        auto stitchedPath = temp.path() + "/stitched.h5";
        SpectrumProcessing::save(stitchedPath, replaced);
        auto reopened = SpectrumProcessing::load(stitchedPath, {});
        expect(reopened.crystals == replaced.crystals && reopened.wideCrystals == wide.crystals &&
                   reopened.pointSources == replaced.pointSources && reopened.sources == replaced.sources,
               "stitched HDF5 retains the original Compton reference and per-point provenance");
        QVector<double> energy, line, compton;
        for (int i = 0; i <= 100; ++i) {
            double e = i * .1, u = 2 * (e - 4.2) / 3.4;
            energy.append(e);
            line.append(8 + .7 * e + (i == 50 ? 13 : 0));
            compton.append(3 + 27 * std::pow(1 + std::expm1(std::log(2.) / 2.5) * u * u, -2.5));
        }
        auto linear = SpectrumProcessing::fit(energy, line, {false, {{0, 2}, {8, 10}}}, cancel);
        expect(close(linear.corrected[50], 13, 1e-9) && close(linear.corrected[10], 0, 1e-9),
               "linear fit subtracts independent sidebands and keeps excluded signal peak");
        auto secondLine = line;
        for (qsizetype i = 0; i < energy.size(); ++i)
            secondLine[i] = 20 - .4 * energy[i] + (i == 50 ? 5 : 0);
        ProcessedSpectrum edgeOnly{energy, {{"A", line}, {"B", secondLine}}, "counts", {}};
        edgeOnly.geometry["B"].qMean = 12;
        const auto edgePath = temp.path() + "/edge_only.h5";
        SpectrumProcessing::save(edgePath, edgeOnly);
        const auto direct = analyzeLinear(SpectrumProcessing::load(edgePath, {}), {{0, 2}, {8, 10}}, cancel);
        expect(direct.failures.isEmpty() && direct.crystals.size() == 2 && !direct.crystals.value("B").high &&
                   close(direct.crystals.value("A").fit.corrected.value(50), 13) &&
                   close(direct.crystals.value("B").fit.corrected.value(50), 5) &&
                   direct.parameters.value("background_mode").toString() == "linear_only",
               "edge-only linear processing fits high-q and unknown-q crystals independently without a wide "
               "scan");
        auto pearson = SpectrumProcessing::fit(energy, compton, {true, {{0, 10}}}, cancel);
        expect(close(pearson.centerEv, 4.2, 1e-5) && close(pearson.fwhmEv, 3.4, 1e-5) &&
                   close(pearson.exponent, 2.5, 1e-4) && close(pearson.corrected[42], 0, 1e-5),
               "Pearson VII fit recovers known Compton center, FWHM and shape");
        auto correctedPath = temp.path() + "/corrected.h5";
        SpectrumProcessing::save(correctedPath, {energy, {{"A", line}, {"B", compton}}, "counts", {}},
                                 {{"A", linear}});
        SpectrumHdfLayout correctedLayout;
        correctedLayout.intensity = "/entry/spectrum/corrected";
        auto corrected = SpectrumProcessing::stitch({correctedPath}, correctedLayout, cancel);
        expect(corrected.crystals.size() == 1 && close(corrected.crystals["A"][50], 13),
               "exported corrected matrix preserves the fitted crystal and omits unprocessed rows");
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
