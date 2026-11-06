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
