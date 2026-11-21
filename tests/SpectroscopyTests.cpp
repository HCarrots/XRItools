#include "ComptonProcessing.h"
#include "ImagingIO.h"
#include "Spectroscopy.h"
#include "SpectrumFixture.h"
#include "SpectrumProcessing.h"
#include <QCoreApplication>
#include <QFile>
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
            ProcessedSpectrum independent = stitched;
            independent.energyEv = {wide.energyEv.first() - 5., 461.1, 470.1, 480.1, 490.1, 505.1,
                                    wide.energyEv.last() + 5.};
            for (auto &row : independent.crystals)
                row.fill(1., independent.energyEv.size());
            auto nativeWide = wide;
            for (auto it = nativeWide.crystals.cbegin(); it != nativeWide.crystals.cend(); ++it)
                nativeWide.crystalEnergyEv.insert(it.key(), wide.energyEv);
            nativeWide.energyEv.clear();
            const auto separate = analyzeCompton(nativeWide, independent, settings, cancel);
            const auto independentlyCorrected =
                std::find_if(separate.crystals.cbegin(), separate.crystals.cend(),
                             [](const ComptonCrystal &item) { return item.high; });
            expect(independentlyCorrected != separate.crystals.cend() &&
                       independentlyCorrected->fit.corrected.size() == independent.energyEv.size() &&
                       std::isnan(independentlyCorrected->fit.corrected.first()) &&
                       std::isnan(independentlyCorrected->fit.corrected.last()) &&
                       std::isfinite(independentlyCorrected->fit.corrected[3]),
                   "independent near-edge grids retain target points with no Compton correction outside wide coverage");
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
        {
            const auto reference = Spectroscopy::pixelSpectra(xrs, s, rois.items()[0], {8, 8}, cancel);
            expect(reference.pixels.size() == 4 && close(reference.intensity[0][30], 11.875),
                   "spectral reference applies glitch correction and I0 exactly once");
            PixelSpectra mixture;
            mixture.energyKev = {1, 2, 3, 4};
            mixture.intensity = {{1, 1, 1, 1}, {4, 4, 4, 4},   {0, 0, 3, 3},
                                 {1, 1, 7, 7}, {9, 9, 12, 12}, {2, 2, 11, 11}};
            const QVector<int> rows{0, 1, 2, 3, 4, 5};
            for (auto method : {UnmixMethod::NNMF, UnmixMethod::PCA, UnmixMethod::ICA}) {
                const auto factors = decomposePixels(mixture, rows, method, 2, cancel);
                expect(factors.spectra.size() == 2 && factors.relativeResidual < 1e-3,
                       "two-component model reconstructs rank-two sample/background spectra");
                const int target =
                    std::abs(factors.covariance[0][2]) > std::abs(factors.covariance[1][2]) ? 0 : 1;
                const auto &scores = factors.covariance[target];
                expect(close(scores[0], 0) && std::abs(scores[2]) > 0.1 && close(scores[3], 2 * scores[2]) &&
                           close(scores[2], scores[4]),
                       "sample covariance rejects flat background and is offset invariant");
            }
            expect(throws([&] { decomposePixels(mixture, {0}, UnmixMethod::PCA, 2, cancel); }),
                   "one surviving pixel cannot be decomposed into two components");
            auto options = s;
            options.spectralMasks[0].pixels = {{1, 1}};
            const auto refined = Spectroscopy::calibrate(options, rois, cancel);
            auto target = Spectroscopy::extract(options, refined, {reverse}, cancel);
            expect(refined.pixels[0].size() == 4 && refined.filteredPixels[0] == QVector<int>{9} &&
                       refined.adjusted.items()[0].mask == options.spectralMasks[0].pixels &&
                       close(target[0].filtered[0][40], 12.5) && target[0].raw[0][40] > 12.5,
                   "spectral mask preserves coarse comparison and transfers to independent target scan");
        }
        auto c = Spectroscopy::calibrate(s, rois, cancel);
        expect(c.fits[0].accepted && c.fits[1].accepted && c.fits[60].accepted && !c.fits[2].accepted,
               "three valid crystals fit, manually invalid crystals skipped");
        expect(c.pixels[0].size() == 4 && c.filteredPixels[0].size() == 2,
               "threshold mask strictly selects pixels above 15 percent");
        expect(close(c.scan.filtered[0][40], 128.75, 1e-6) && close(c.scan.raw[0][40], 139.05, 1e-6),
               "elastic normalization and raw/masked integrals agree with independent values");
        {
            auto masked = rois;
            const QVector<QSize> sizes{{8, 8}, {6, 6}};
            const QVector<QPoint> mask{{1, 1}, {1, 2}, {2, 2}};
            CrystalRois reopened;
            QString error;
            expect(masked.assignMask(0, mask, sizes[0]) &&
                       reopened.restore(masked.serialize(elastic, sizes), sizes, &error) &&
                       reopened.items()[0].mask == mask && reopened.items()[1].mask.isEmpty(),
                   "shared scheme round-trips mask and rectangular ROIs together");
            expect(CrystalRois::decodeMask(QJsonArray{QJsonArray{0, 2}, QJsonArray{1, 2}}, sizes[0]).isEmpty() &&
                       CrystalRois::decodeMask(QJsonArray{QJsonValue(QJsonArray{63, 2})}, sizes[0]).isEmpty(),
                   "shared mask decoder rejects repeated and out-of-range pixels");
            auto options = s;
            options.autoRoi = true;
            options.roiSize = 2;
            auto fitted = Spectroscopy::calibrate(options, reopened, cancel);
            auto maskedScans = Spectroscopy::extract(options, fitted, {xrs}, cancel);
            expect(fitted.fits[0].accepted && fitted.adjusted.items()[0].mask == mask &&
                       fitted.pixels[0] == QVector<int>{9, 17, 18} && fitted.filteredPixels[0] == fitted.pixels[0] &&
                       close(fitted.scan.filtered[0][40], 113.3, 1e-6) &&
                       close(maskedScans[0].filtered[0][40], 13.75, 1e-6),
                   "elastic and XRS scans integrate the exact saved mask despite rectangular ROI options");
        }
        auto scans = Spectroscopy::extract(s, c, {xrs, reverse}, cancel);
        expect(scans.size() == 2 && scans[0].glitches == QVector<int>{30} &&
                   scans[1].glitches == QVector<int>{50},
               "isolated I0 glitch detected in ascending and descending scans");
        expect(close(scans[0].i0Before[30], .2) && close(scans[0].i0After[30], 2) &&
                   close(scans[0].filtered[0][30], 1.25 * (10 + .5 * 3.75)),
               "glitch repairs both frame and I0 before normalization");
        {
            QTemporaryDir irregular;
            auto options = s;
            for (bool descending : {false, true}) {
                const auto file = SpectrumFixture::create(irregular.path(), false, descending, true, true);
                const auto repaired = Spectroscopy::extract(s, c, {file}, cancel).first();
                const int i = descending ? 50 : 30;
                const double w = .625;
                const double flux = (1 - w) * 2 + w * 2.04;
                const double counts = (1 - w) * 2 * (10 + .5 * 29 * .125) +
                                      w * 2.04 * (10 + .5 * 31 * .125);
                expect(repaired.glitches == QVector<int>{i} && close(repaired.i0After[i], flux) &&
                           close(repaired.filtered[0][i], 1.25 * counts / flux) &&
                           close(repaired.filtered[60][i], 3.75 * counts / flux),
                       "an isolated 8% dropout repairs I0 and both detectors on an irregular energy grid");
                options.glitchDropFraction = .30;
                const auto strict = Spectroscopy::extract(options, c, {file}, cancel).first();
                expect(strict.glitches.isEmpty() && strict.i0Before == strict.i0After &&
                           strict.filtered[0][i] == 0, "30% threshold leaves shallow drop and lost frame intact");
            }
        }
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
        ProcessedSpectrum sampled;
        sampled.crystalEnergyEv = {{"VB-A1", {0., 1., 2., 3., 4.}}, {"VB-A2", {.5, 1.5, 2.5, 3.5}}};
        sampled.crystals = {{"VB-A1", {1., 2., 3., 4., 5.}}, {"VB-A2", {10., 12., 14., 16.}}};
        const auto sum = SpectrumProcessing::selectSamples(sampled, {"VB-A1", "VB-A2"}, 0., 4.);
        const auto within = SpectrumProcessing::selectSamples(sampled, {"VB-A1", "VB-A2"}, 1., 3.);
        const auto single = SpectrumProcessing::selectSamples(sampled, {"VB-A1"}, 1., 3.);
        expect(sum.energyEv == QVector<double>{1., 2., 3.} && sum.intensity == QVector<double>{13., 16., 19.} &&
                   within.energyEv == QVector<double>{2.} && within.intensity == QVector<double>{16.} &&
                   single.energyEv == QVector<double>{1., 2., 3.} && single.intensity == QVector<double>{2., 3., 4.},
               "preview sums on shared coverage and uses only samples inside the selected energy interval");
        sampled.intensityUnit = "counts";
        sampled.sources = {"compton"};
        ProcessedSpectrum detailed{{1.25, 2.25, 2.75},
                                    {{"VB-A1", {60., 80., 100.}}, {"VB-A2", {120., 140., 160.}}},
                                    "counts", {"near-a", "near-b"}};
        detailed.geometry["VB-A1"].qMean = 4.;
        const auto joined = SpectrumProcessing::spliceCompton(sampled, detailed, true, cancel);
        const auto &joinedA = joined.crystals["VB-A1"];
        expect(joined.energyEv == QVector<double>{0., .5, 1., 1.25, 2.25, 2.75, 3., 3.5, 4.} &&
                   joinedA[0] == 1. && std::isnan(joinedA[1]) && joinedA[3] == 30. && joinedA[4] == 40. &&
                   joinedA[5] == 50. && joinedA[6] == 4. && std::isnan(joinedA[7]) &&
                   joined.pointSources["VB-A1"][3] == 1 && joined.pointSources["VB-A1"][6] == 0 &&
                   joined.geometry.value("VB-A1").qMean == 4.,
               "contained near-edge ranges replace coarse points while retaining distinct per-crystal energies");
        detailed.energyEv = {7., 7.5, 8.};
        const auto disjoint = SpectrumProcessing::spliceCompton(sampled, detailed, false, cancel);
        expect(disjoint.energyEv == QVector<double>{0., .5, 1., 1.5, 2., 2.5, 3., 3.5, 4., 7., 7.5, 8.} &&
                   disjoint.crystals["VB-A1"].last() == 100. && disjoint.crystals["VB-A2"].last() == 160.,
               "disjoint ranges retain only measured energies and unnormalized intensity when requested");
        detailed.energyEv = {3.25, 4., 5.};
        const auto partial = SpectrumProcessing::spliceCompton(sampled, detailed, false, cancel);
        expect(partial.energyEv == QVector<double>{0., .5, 1., 1.5, 2., 2.5, 3., 3.25, 4., 5.} &&
                   partial.crystals["VB-A1"][8] == 80.,
               "partially overlapping scans keep the near-edge points and the unoverlapped wide-scan tail");
        auto pearson = SpectrumProcessing::fit(energy, compton, {true, {{0, 10}}}, cancel);
        expect(close(pearson.centerEv, 4.2, 1e-5) && close(pearson.fwhmEv, 3.4, 1e-5) &&
                   close(pearson.exponent, 2.5, 1e-4) && close(pearson.corrected[42], 0, 1e-5),
               "Pearson VII fit recovers known Compton center, FWHM and shape");
        auto withGap = line;
        withGap[55] = NAN;
        StagedBackgroundSettings plan;
        plan.windows[int(BackgroundComponent::Linear)] = {{0, 10}};
        plan.protectedWindows = {{4.8, 5.2}};
        const auto protectedLine = fitBackgroundStages(plan, energy, withGap, {}, cancel);
        expect(close(protectedLine.total.corrected[50], 13, 1e-9) &&
                   std::isnan(protectedLine.total.background[55]) && std::isnan(protectedLine.total.corrected[55]),
               "positive near-edge signal is excluded from the fit and masked samples stay masked");
        plan.protectedWindows.clear();
        plan.windows[int(BackgroundComponent::ComptonPeak)] = {{0, 10}};
        QVector<double> reduced;
        for (const double value : compton)
            reduced.append(value - 3.);
        const auto peakOnly = fitBackgroundComponent(BackgroundComponent::ComptonPeak, plan, energy, reduced, {}, cancel);
        expect(close(peakOnly.peakHeight, 27, 1e-5) && close(peakOnly.centerEv, 4.2, 1e-5) &&
                   close(peakOnly.corrected[0], 0, 1e-5),
               "the Compton peak component has a fixed zero baseline");

        plan.hfTable = temp.filePath("test-profiles.dat");
        QFile hfTable(plan.hfTable);
        if (!hfTable.open(QIODevice::WriteOnly))
            throw std::runtime_error("could not write the synthetic HF table");
        hfTable.write("#S 6 C\n#UOCCUP 2 2\n#UBIND 8 284\n#L pz total L1 K\n"
                      "0 2 1 1\n1 1 .2 .8\n2 .6 .1 .5\n5 .2 .03 .17\n10 .05 .01 .04\n");
        hfTable.close();
        plan.composition = "C:1";
        expect(backgroundHfShells(plan.hfTable, plan.composition).value("C:2").bindingEv == 284.,
               "the selectable HF shell identifiers preserve the tabulated binding energies");
        const auto builtinShells = backgroundHfShells({}, "O:2, Ni:1, Kr:1, K:1");
        const auto selection = resolveHfShellSelection(builtinShells, {"O:K", "Ni:L2", "Ni:L3", "Kr:L23", "K:M4"});
        expect(selection.identifiers == QStringList{"O:1", "Ni:3", "Kr:3", "Kr:4"} &&
                   selection.missing == QStringList{"K:M4"} && builtinShells.value("Ni:3").edge == "L23",
               "absorption edges resolve to the actual Biggs profiles, including combined and split spin-orbit shells");
        expect(resolveHfShellSelection(builtinShells, {"Ni:3"}).identifiers == QStringList{"Ni:3"},
               "legacy numbered selections remain valid");
        plan.hfShells = {"C:K"};
        plan.windows[int(BackgroundComponent::HfCore)] = {{350, 900}};
        plan.windows[int(BackgroundComponent::Asymmetry)] = {{15, 180}};
        const SpectrumGeometry geometry{5., 9.69, 60.};
        QVector<double> broadEnergy, positive, odd;
        for (int e = 0; e <= 1000; e += 5) {
            broadEnergy.append(e);
            positive.append(e == 500 ? NAN : 50.);
            odd.append(e < 93 ? 10. : -10.);
        }
        const auto hfFit = fitBackgroundComponent(BackgroundComponent::HfCore, plan, broadEnergy, positive,
                                                  geometry, cancel);
        expect(hfFit.parameters.value("scale").toDouble() > 0. && hfFit.background[30] == 0. &&
                   hfFit.background[80] > 0. && std::isnan(hfFit.background[100]) &&
                   hfFit.parameters.value("resolved_shells").toObject().value("C:2") == "K",
               "the custom table's named K shell drives the fit and its actual identifier is recorded");
        const auto oddFit = fitBackgroundComponent(BackgroundComponent::Asymmetry, plan, broadEnergy, odd,
                                                   geometry, cancel);
        expect(oddFit.background[12] * oddFit.background[24] < 0.,
               "the asymmetry correction is signed on opposite sides of the crystal's Compton center");
        plan.hfShells.clear();
        expect(throws([&] { fitBackgroundComponent(BackgroundComponent::HfCore, plan, broadEnergy, positive,
                                                  geometry, cancel); }),
               "HF fitting requires explicit shell selection");
        const auto stagedPath = temp.filePath("staged-background.h5");
        SpectrumProcessing::save(stagedPath, {energy, {{"A", withGap}}, "counts", {}},
                                 {{"A", protectedLine.total}},
                                 {{"linear", {{"A", protectedLine.total.background}}}});
        ImagingH5 stagedFile(stagedPath);
        const auto savedLine = stagedFile.read("/entry/spectrum/background_linear");
        expect(close(savedLine[50], protectedLine.total.background[50]) && std::isnan(savedLine[55]),
               "HDF5 export keeps each background component and its missing samples");
        const auto textPath = temp.filePath("staged-background.txt");
        const ProcessedSpectrum textSpectrum{energy, {{"A", withGap}}, "counts", {}};
        const QMap<QString, CrystalBackground> textFits{{"A", protectedLine.total}};
        const QMap<QString, QMap<QString, QVector<double>>> textComponents{
            {"linear", {{"A", protectedLine.total.background}}}};
        SpectrumProcessing::saveText(textPath, textSpectrum, textFits, textComponents);
        QFile savedText(textPath);
        if (!savedText.open(QIODevice::ReadOnly))
            throw std::runtime_error("could not open the exported text spectrum");
        const auto lines = savedText.readAll().split('\n');
        savedText.close();
        const auto signal = lines.value(51).split('\t'), gap = lines.value(56).split('\t');
        expect(lines.size() == energy.size() + 2 &&
                   lines.first() == "Energy Transfer (eV)\tA_input\tA_background_linear\tA_background_total\tA_corrected" &&
                   signal.size() == 5 && signal.value(0).toDouble() == energy[50] &&
                   signal.value(2).toDouble() == savedLine[50] && close(signal.value(4).toDouble(), 13.) &&
                   gap == QList<QByteArray>{QByteArray::number(energy[55], 'g', 17), "nan", "nan", "nan", "nan"},
               "text export preserves the energy grid, full-precision intensities, background columns and gaps");
        expect(throws([&] { SpectrumProcessing::saveText(textPath, textSpectrum, textFits, textComponents); }),
               "text export never overwrites an existing spectrum");
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
