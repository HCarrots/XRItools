#include "BondVolume.h"
#include "CrystalRoi.h"
#include "ImagingFixture.h"
#include <QApplication>
#include <QJsonDocument>
#include <QTemporaryDir>
#include <iostream>
namespace {
int failures = 0;
void expect(bool condition, const char *text) {
    std::cout << (condition ? "PASS " : "FAIL ") << text << std::endl;
    if (!condition)
        ++failures;
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
    QApplication app(argc, argv);
    std::atomic_bool stop{false};
    try {
        if (argc == 4 && QString::fromLocal8Bit(argv[1]) == "--bond-reference") {
            auto config = Imaging::loadConfig(QString::fromLocal8Bit(argv[2]));
            auto result = BondImaging::run(config, stop, [](const QString &text, int progress) {
                std::cout << progress << " " << text.toUtf8().constData() << std::endl;
            });
            const auto path = QString::fromLocal8Bit(argv[3]);
            BondImaging::save(path, result, stop);
            auto reopened = BondImaging::load(path);
            expect(reopened.axes == result.axes && reopened.converged == result.converged,
                   "real four-dimensional result round-trips");
            std::cout << QJsonDocument(reopened.diagnostics).toJson().constData();
            return failures ? 1 : 0;
        }
        if (argc == 3) {
            auto command = QString::fromLocal8Bit(argv[1]), path = QString::fromLocal8Bit(argv[2]);
            if (command == "--fixture") {
                ImagingFixture::create(path);
                return 0;
            }
            auto config = Imaging::loadConfig(path);
            if (command == "--validate") {
                auto report = Imaging::validate(config, stop).report;
                std::cout << QJsonDocument(report).toJson().constData();
                return 0;
            }
            if (command == "--run") {
                auto result = Imaging::run(config, stop, [](const QString &text, int progress) {
                    std::cout << progress << " " << text.toUtf8().constData() << std::endl;
                });
                std::cout << QJsonDocument(result.summary).toJson().constData();
                return 0;
            }
        }
        auto weights = Imaging::matchedWeights({1, 2, 4}, {1, 1.5, 4});
        expect(std::abs(std::accumulate(weights.first.begin(), weights.first.end(), 0.) - 3) < 1e-12,
               "reference interpolation preserves integrated constant spectrum");
        expect(throws([&] { Imaging::matchedWeights({1, 2}, {.9, 2}); }), "no energy extrapolation");
        {
            QTemporaryDir tmp;
            auto canonical = ImagingFixture::create(tmp.filePath("canonical"), false, false, false, true);
            auto reordered = ImagingFixture::create(tmp.filePath("reordered"), false, true, true, true);
            const auto v = BondImaging::run(canonical, stop);
            const auto w = BondImaging::run(reordered, stop);
            expect(v.axes == w.axes && v.converged == QVector<bool>{true, true} &&
                       w.converged == v.converged && v.voxelCount() == 21 * 2 * 2,
                   "4D reconstruction preserves both independent energy channels with the imaging X grid");
            double agreement = 0;
            for (qsizetype i = 0; i < v.intensity.size(); ++i)
                agreement = std::max(agreement, std::abs(v.intensity[i] - w.intensity[i]));
            const auto a = BondImaging::energySlice(v, v.intensity, v.axes[0], 1, true);
            const auto b = BondImaging::energySlice(v, v.intensity, v.axes[0], 2, true);
            const auto mid = BondImaging::energySlice(v, v.intensity, v.axes[0], 1.5, true);
            const auto peak = std::max_element(b.cbegin(), b.cend()) - b.cbegin();
            expect(agreement < 1e-8 && peak / 4 == 10 && b[peak] > 0 &&
                       std::abs(b[peak] / a[peak] - 3.) < 1e-6 && std::abs(mid[peak] - 2 * a[peak]) < 1e-6,
                   "reordered duplicate frames preserve a known spectral peak and energy interpolation");
            const auto forwardProjection = BondImaging::prepareContrast(canonical, 1, 2, {}, stop);
            const auto reverseProjection = BondImaging::prepareContrast(reordered, 2, 1, {}, stop);
            const auto forward = BondImaging::reconstructContrast(forwardProjection, canonical, stop);
            const auto reverse = BondImaging::reconstructContrast(reverseProjection, reordered, stop);
            bool reverses = forward.values.size() == reverse.values.size();
            for (qsizetype i = 0; reverses && i < forward.values.size(); ++i)
                reverses = (std::isnan(forward.values[i]) && std::isnan(reverse.values[i])) ||
                           std::abs(forward.values[i] + reverse.values[i]) < 1e-8;
            expect(forward.axes == v.axes.sliced(1) && forward.axes == reverse.axes && reverses &&
                       forward.values[peak] > 0 && reverse.values[peak] < 0,
                   "projection differences preserve sign across energy reversal and reordered duplicate frames");
            QFile changedRoi(canonical["input"].toObject()["roi_table"].toString());
            const bool changed = changedRoi.open(QIODevice::Append) && changedRoi.write("\n") == 1;
            changedRoi.close();
            expect(changed && throws([&] { BondImaging::reconstructContrast(forwardProjection, canonical, stop); }),
                   "spatial reconstruction rejects projections after an input file changes");
            expect(throws([&] { BondImaging::energySlice(v, v.intensity, v.axes[0], .9, true); }),
                   "4D energy interpolation rejects extrapolation");
            auto bad = v;
            bad.converged[1] = false;
            expect(throws([&] { BondImaging::energySlice(bad, bad.intensity, bad.axes[0], 1.5, true); }) &&
                       !BondImaging::energySlice(bad, bad.intensity, bad.axes[0], 2, false).isEmpty(),
                   "unconverged channels stay available for raw diagnostics only");
            auto missing = ImagingFixture::create(tmp.filePath("missing"), true);
            expect(throws([&] { BondImaging::validate(missing, stop); }),
                   "4D validation rejects missing cells");
            auto resultPath = tmp.filePath("volume.h5");
            BondImaging::save(resultPath, v, stop);
            auto read = BondImaging::load(resultPath);
            expect(read.axes == v.axes && read.intensity == v.intensity && read.converged == v.converged &&
                       read.sensitivity == v.sensitivity && read.diagnostics == v.diagnostics &&
                       throws([&] { BondImaging::save(resultPath, v, stop); }),
                   "4D HDF5 round-trip retains data and diagnostics without overwriting");
            std::atomic_bool cancelled{true};
            auto cancelPath = tmp.filePath("cancelled.h5");
            expect(throws([&] { BondImaging::run(canonical, cancelled); }) &&
                       throws([&] { BondImaging::save(cancelPath, v, cancelled); }) &&
                       !QFile::exists(cancelPath),
                   "4D cancellation does not publish a partial result");
        }
        {
            BondVolume v;
            v.axes = {{0, 1, 2, 3, 4}, {0, 2}, {0, 1}, {0, 1}};
            v.converged.fill(true, 5);
            v.sensitivity = {1, 1};
            for (int e = 0; e < 5; ++e)
                for (qsizetype i = 0; i < v.voxelCount(); ++i)
                    v.intensity.append(10 + 2 * e + (e == 2 ? (i ? 4 : -4) : 0));
            BondSpectralSettings settings{{{0, 1}, {3, 4}}, QPair<double, double>{1, 3}};
            auto corrected = BondImaging::process(v, v.axes[0], settings, stop);
            expect(std::isnan(corrected[16]) && std::abs(corrected[17] - 1) < 1e-12 &&
                       std::abs(corrected[9]) < 1e-12,
                   "per-voxel subtraction removes a linear baseline and rejects nonpositive areas");
            settings.normalization.reset();
            const auto spectralDifference = [&](double ea, double eb, qsizetype voxel) {
                const auto weights = BondImaging::contrastWeights(v.axes[0], ea, eb, settings);
                double value = 0;
                for (qsizetype e = 0; e < weights.size(); ++e)
                    if (weights[e] != 0.)
                        value += weights[e] * v.intensity[e * v.voxelCount() + voxel];
                return value;
            };
            expect(std::abs(spectralDifference(.5, 2, 0) + 4.) < 1e-12 &&
                       std::abs(spectralDifference(.5, 2, 1) - 4.) < 1e-12,
                   "linear-background difference supports interpolated energies and both signs");

            const QVector<double> factor{1, 2, 4, 2, 1}, valence{3, 4, 6, 4, NAN};
            for (qsizetype e = 0; e < v.axes[0].size(); ++e)
                for (qsizetype i = 0; i < v.voxelCount(); ++i)
                    v.intensity[e * v.voxelCount() + i] =
                        (i + 1) * (valence[e] + 2 + e + (e == 2 ? 6 : 0)) / factor[e];
            settings.background = {{0, 1}};
            settings.compton = BondSpectralReference{factor, {5, 7, 10, 9, NAN}, valence, {0, 1}};
            corrected = BondImaging::process(v, v.axes[0], settings, stop);
            expect(std::abs(corrected[16] - 1.5) < 1e-12 && std::abs(corrected[17] - 3.) < 1e-12 &&
                       std::abs(corrected[24]) < 1e-12 && std::isnan(corrected[32]),
                   "Compton transfer removes the matched valence and residual while retaining spatial amplitude");
            expect(std::abs(spectralDifference(.5, 2.5, 1) - 1.5) < 1e-12 &&
                       throws([&] { spectralDifference(2, 3.5, 1); }),
                   "Compton difference matches preprocessing and rejects energies outside reference coverage");
            settings.normalization = QPair<double, double>{1, 3};
            corrected = BondImaging::process(v, v.axes[0], settings, stop);
            expect(std::abs(corrected[16] - 1.) < 1e-12 && std::abs(corrected[17] - 1.) < 1e-12,
                   "optional area normalization runs after the Compton background has been removed");
            settings.normalization = QPair<double, double>{1, 4};
            expect(throws([&] { BondImaging::process(v, v.axes[0], settings, stop); }),
                   "area normalization rejects gaps in the Compton model coverage");
            const auto mesh = BondImaging::surface(v, {0, 0, 0, 0, 2, 2, 2, 2}, .5, stop);
            expect(!mesh.isEmpty() &&
                       std::all_of(mesh.cbegin(), mesh.cend(),
                                   [](const BondTriangle &triangle) {
                                       return std::all_of(triangle.cbegin(), triangle.cend(), [](auto p) {
                                           return std::abs(p.x() - .5f) < 1e-6;
                                       });
                                   }),
                   "isosurface interpolation recovers the known physical plane X=0.5 mm");
        }
        ImageMatrix h = ImageMatrix::Zero(2, 2), b(2, 2);
        h(0, 0) = 1;
        h(1, 1) = 1e-30;
        b << 1, 2, 1e-10, 2e-10;
        ImagingSolverInfo info;
        auto rho = Imaging::nnls(h, b, 20000, 1e-7, stop, &info);
        expect(info.converged && std::abs(rho(0, 0) - 1) < 1e-12 && rho.row(1).norm() < 1e-12,
               "NNLS projects numerical zero modes out of H and rhs");
        ImageMatrix op(20, 4), truth(4, 3);
        for (int r = 0; r < 20; ++r)
            for (int col = 0; col < 4; ++col)
                op(r, col) = .2 + std::abs(std::sin(r * 1.7 + col * 2.3));
        truth << 1, 2, .2, 3, .2, 1, .5, 1, 2, 1, .6, 3;
        ImageMatrix resolved =
            Imaging::nnls(op.transpose() * op, op.transpose() * op * truth, 20000, 1e-10, stop, &info);
        expect(info.converged && (resolved - truth).norm() < 1e-6,
               "resolved NNLS recovers known nonnegative multicolumn truth");
        ImageMatrix tvTruth(3, 4);
        tvTruth.row(0).setConstant(.2);
        tvTruth.row(1).setConstant(1);
        tvTruth.row(2).setConstant(.4);
        std::array<double, 3> tvSpacing{.5, .5, .5};
        auto objective = [&](const ImageMatrix &value) {
            auto g = Imaging::gradient(value, 2, 2, tvSpacing);
            double variation = 0;
            for (int i = 0; i < value.size(); ++i)
                variation += std::sqrt(g[0].data()[i] * g[0].data()[i] + g[1].data()[i] * g[1].data()[i] +
                                       g[2].data()[i] * g[2].data()[i]);
            return .5 * (value - tvTruth).squaredNorm() + .01 * variation;
        };
        auto tvVolume = Imaging::tv(ImageMatrix::Identity(3, 3), tvTruth, 2, 2, tvSpacing, .01, 10000, 1e-8,
                                    tvTruth, stop, &info);
        expect(info.converged && tvVolume.minCoeff() >= 0 && objective(tvVolume) < objective(tvTruth),
               "nonnegative isotropic TV converges and lowers its original objective");
        ImageMatrix volume(4, 6);
        for (int i = 0; i < volume.size(); ++i)
            volume.data()[i] = std::sin(i);
        std::array<ImageMatrix, 3> dual;
        for (int a = 0; a < 3; ++a) {
            dual[a].resize(4, 6);
            for (int i = 0; i < dual[a].size(); ++i)
                dual[a].data()[i] = std::cos(i + a);
        }
        std::array<double, 3> spacing{.07, .055, .054};
        auto grad = Imaging::gradient(volume, 2, 3, spacing);
        double left = 0;
        for (int a = 0; a < 3; ++a)
            left += (grad[a].array() * dual[a].array()).sum();
        double right = (volume.array() * Imaging::adjoint(dual, 2, 3, spacing).array()).sum();
        expect(std::abs(left - right) < 1e-10,
               "TV gradient and adjoint satisfy boundary inner-product identity");
        QTemporaryDir temp;
        auto config = ImagingFixture::create(temp.path());
        auto validation = Imaging::validate(config, stop);
        expect(validation.calibration.axes[0].size() == 9 && validation.sample.axes[0].size() == 2 &&
                   validation.rois[0].key == "1_test",
               "plan literal, label override and physical frame mapping");
        auto result = Imaging::run(config, stop);
        expect(!QFile::exists(temp.filePath("output/reconstruction.am")),
               "Avizo output is opt-in for existing imaging configurations");
        expect(result.primary.rows() == 21 && result.primary.cols() == 4 && result.axes[2].last() == 0,
               "full XYZ output keeps high-Z boundary");
        Eigen::Index peak;
        result.primary.col(0).maxCoeff(&peak);
        expect(std::abs(result.axes[0][peak]) < 1e-9 && std::abs(result.primary.col(0).sum() - .2) < .01,
               "point location and density recover known physical scale");
        expect(result.summary["tv_converged"].toBool() && result.primary.minCoeff() >= 0,
               "nonnegative TV convergence is reported");
        expect(throws([&] { Imaging::run(config, stop); }), "existing output directory is never overwritten");
        QTemporaryDir repeated;
        auto repeatedConfig = ImagingFixture::create(repeated.path(), false, true);
        auto repeatedValidation = Imaging::validate(repeatedConfig, stop);
        auto repeatedResult = Imaging::run(repeatedConfig, stop);
        expect(repeatedValidation.report["sample"].toObject()["repeated_frames"].toInt() == 1 &&
                   repeatedValidation.report["calibration"].toObject()["repeated_frames"].toInt() == 1 &&
                   (repeatedResult.primary - result.primary).norm() < 1e-9,
               "reversed frame order and repeats pool counts/dose without changing reconstruction");
        QTemporaryDir energyOuter;
        auto alternate = ImagingFixture::create(energyOuter.path(), false, false, true);
        auto detected = Imaging::validate(alternate, stop);
        expect(detected.calibration.loopOrder == QStringList({"DCM", "X"}) &&
                   detected.sample.loopOrder == QStringList({"DCM", "Z", "Y"}) &&
                   detected.calibration.axes == validation.calibration.axes &&
                   detected.sample.axes == validation.sample.axes,
               "automatic scan loops support DCM as the outer axis with canonical physical grids");
        auto loops = alternate["scan_loops"].toObject();
        loops["calibration_order"] = "DCM,X";
        loops["calibration_x_points"] = 9;
        loops["calibration_dcm_points"] = 3;
        loops["sample_order"] = "DCM,Z,Y";
        loops["sample_y_points"] = 2;
        loops["sample_z_points"] = 2;
        loops["sample_dcm_points"] = 2;
        alternate["scan_loops"] = loops;
        auto alternateResult = Imaging::run(alternate, stop);
        expect((alternateResult.primary - result.primary).norm() < 1e-9,
               "manually selected loops and point counts preserve physical reconstruction");
        loops["sample_dcm_points"] = 3;
        alternate["scan_loops"] = loops;
        expect(throws([&] { Imaging::validate(alternate, stop); }),
               "manual point count rejects missing cells");
        loops["sample_dcm_points"] = 2;
        loops["sample_order"] = "Y,Z,DCM";
        alternate["scan_loops"] = loops;
        expect(throws([&] { Imaging::validate(alternate, stop); }),
               "manual order rejects mismatched scan plan");
        CrystalRois shared;
        for (int i = 0; i < 75; ++i)
            shared.invalidate(i);
        shared.assign(0, {0, 0, 19, 47}, {19, 47});
        auto sharedPath = temp.path() + "/shared.roi.json";
        QFile sharedFile(sharedPath);
        expect(sharedFile.open(QIODevice::WriteOnly) &&
                   sharedFile.write(shared.serialize("foil.nxs", {{19, 47}, {9, 9}})) > 0,
               "write shared 75-crystal scheme fixture");
        sharedFile.close();
        auto common = config;
        auto commonInput = common["input"].toObject();
        commonInput["roi_table"] = sharedPath;
        common["input"] = commonInput;
        common["roi"] = QJsonObject{{"excluded_ids", QJsonArray()}, {"label_overrides", QJsonObject()}};
        auto commonRec = common["reconstruction"].toObject();
        commonRec["anchor_roi"] = "0_VB-A1";
        common["reconstruction"] = commonRec;
        auto commonRuntime = common["runtime"].toObject();
        commonRuntime["output_dir"] = temp.path() + "/shared-output";
        common["runtime"] = commonRuntime;
        auto commonValidation = Imaging::validate(common, stop);
        expect(commonValidation.rois.size() == 60 && !commonValidation.rois[0].excluded &&
                   commonValidation.rois[1].excluded,
               "shared ROI adapter preserves canonical ID and invalid crystals for the selected detector");
        auto commonResult = Imaging::run(common, stop);
        expect((commonResult.primary - result.primary).norm() < 1e-9,
               "shared ROI and legacy TSV produce the same physical reconstruction");
        {
            const int omitted = 23 * 19 + 9;
            auto mask = shared.items()[0].pixels();
            mask.removeAt(omitted);
            expect(shared.assignMask(0, mask, {19, 47}), "assign nonrectangular shared ROI");
            QFile maskFile(temp.filePath("masked.roi.json"));
            expect(maskFile.open(QIODevice::WriteOnly) &&
                       maskFile.write(shared.serialize("foil.nxs", {{19, 47}, {9, 9}})) > 0,
                   "write common mask scheme");
            maskFile.close();
            auto masked = common;
            auto input = masked["input"].toObject();
            input["roi_table"] = maskFile.fileName();
            masked["input"] = input;
            auto runtime = masked["runtime"].toObject();
            runtime["output_dir"] = temp.filePath("masked-output");
            masked["runtime"] = runtime;
            const auto maskRun = Imaging::run(masked, stop);
            ImagingH5 original(commonResult.directory + "/solver_diagnostics.h5"),
                maskedDiagnostics(maskRun.directory + "/solver_diagnostics.h5");
            const auto pixelPath = QStringLiteral("rois/0_VB-A1/pixels_rc");
            const auto samplePath = QStringLiteral("rois/0_VB-A1/measured_pyz");
            auto fullMeasured = original.read(samplePath);
            auto fullPixels = original.read(pixelPath);
            const bool nonzero = fullMeasured[omitted * 4] > 0;
            fullMeasured.remove(omitted * 4, 4);
            fullPixels.remove(omitted * 2, 2);
            expect(nonzero && maskedDiagnostics.read(pixelPath) == fullPixels &&
                       maskedDiagnostics.read(samplePath) == fullMeasured,
                   "imaging reconstructs using precisely the saved mask pixels and per-frame counts");

            auto spectral = ImagingFixture::create(temp.filePath("mask-spectral"), false, false, false, true);
            auto full = common;
            auto fullInput = full["input"].toObject();
            fullInput["sample_nxs"] = spectral["input"].toObject()["sample_nxs"];
            full["input"] = fullInput;
            input["sample_nxs"] = fullInput["sample_nxs"];
            masked["input"] = input;
            const auto fullContrast = BondImaging::prepareContrast(full, 1, 2, {}, stop);
            const auto maskContrast = BondImaging::prepareContrast(masked, 1, 2, {}, stop);
            const auto &a = fullContrast.differenceByCrystal.at("0_VB-A1");
            const auto &b = maskContrast.differenceByCrystal.at("0_VB-A1");
            expect(a(omitted, 0) > 0 && b.rows() == mask.size() &&
                       b.topRows(omitted).isApprox(a.topRows(omitted)) &&
                       b.bottomRows(b.rows() - omitted).isApprox(a.bottomRows(b.rows() - omitted)),
                   "chemical contrast computes the energy difference only over saved mask pixels");
        }
        for (bool enabled : {false, true}) {
            auto c = config;
            auto runtime = c["runtime"].toObject();
            runtime["output_dir"] = temp.path() + (enabled ? "/tv-limited" : "/tv-disabled");
            c["runtime"] = runtime;
            auto tv = c["tv"].toObject();
            tv["enabled"] = enabled;
            tv["max_iterations"] = 1;
            c["tv"] = tv;
            auto reconstruction = c["reconstruction"].toObject();
            reconstruction["fit_gains"] = false;
            c["reconstruction"] = reconstruction;
            auto r = Imaging::run(c, stop);
            expect(r.summary["solver"] == "joint NNLS" && (r.primary - r.nnls).norm() == 0,
                   "disabled/unconverged TV uses NNLS primary without concealing status");
        }
        QTemporaryDir broken;
        auto bad = ImagingFixture::create(broken.path(), true);
        expect(throws([&] { Imaging::run(bad, stop); }) && !QFileInfo::exists(broken.path() + "/output"),
               "incomplete sample fails before creating output");
        auto wrong = config;
        auto cancelled = config;
        auto cancelledRuntime = cancelled["runtime"].toObject();
        cancelledRuntime["output_dir"] = temp.path() + "/cancelled";
        cancelled["runtime"] = cancelledRuntime;
        expect(throws([&] {
                   Imaging::run(cancelled, stop, [&](const QString &, int percent) {
                       if (percent == 3)
                           stop.store(true);
                   });
               }) &&
                   QFileInfo::exists(temp.path() + "/cancelled/run_status.json") &&
                   !QFileInfo::exists(temp.path() + "/cancelled/reconstruction.h5"),
               "cancellation during a run preserves failure status without publishing a primary result");
        stop.store(false);
        auto tv = wrong["tv"].toObject();
        tv["smooth"] = true;
        wrong["tv"] = tv;
        expect(throws([&] { Imaging::config(wrong, temp.path()); }),
               "unknown imaging config fields rejected");
        stop.store(true);
        expect(throws([&] { Imaging::validate(config, stop); }),
               "imaging cancellation checked before file processing");
    } catch (const std::exception &e) {
        std::cerr << "ERROR " << e.what() << std::endl;
        ++failures;
    }
    return failures ? 1 : 0;
}
