#pragma once
#include "ImagingIO.h"
#include <QDir>
#include <QFile>
#include <cmath>
namespace ImagingFixture {
inline QJsonObject create(const QString &directory, bool missing = false, bool reorderedRepeat = false,
                          bool energyOuter = false, bool spectral = false, bool spatial = false,
                          bool twoCrystals = false, int sampleEnergyPoints = 2,
                          std::array<double, 2> foilEnergyRange = {1., 2.}) {
    if (!QDir().mkpath(directory))
        ImageUtil::fail("Cannot create imaging fixture directory");
    int height = 47, width = twoCrystals ? 38 : 19;
    auto spot = [](double x, int row, int col) {
        double u = (row - 23 - 10 * x) / 1.1, v = (col % 19 - 9 - .3 * x) / 1.3;
        return 2000 * std::exp(-(u * u - .2 * u * v + v * v) / (2 * .99));
    };
    auto write = [](ImagingH5 &f, const QString &path, const QVector<double> &data) {
        f.write(path, {quint64(data.size())}, data.constData(), false);
    };
    {
        ImagingH5 f(directory + "/foil.nxs", true);
        f.writeText(
            "entry/gui/plan_cmd",
            QString(energyOuter
                        ? "('P.grid_scan', (['D.LAMBDA'], 'M.DCM_B5Energy', %1, %2, 3, 'M.Stage_X', -1, 1, 9), {})"
                        : "('P.grid_scan', (['D.LAMBDA'], 'M.Stage_X', -1, 1, 9, 'M.DCM_B5Energy', %1, %2, 3), {})")
                .arg(foilEnergyRange[0]).arg(foilEnergyRange[1]));
        QVector<double> x, energy;
        int frames = reorderedRepeat ? 28 : 27;
        std::vector<quint64> image(frames * height * width);
        for (int frame = 0; frame < frames; ++frame) {
            int index = reorderedRepeat ? (frame < 27 ? 26 - frame : 0) : frame;
            if (energyOuter)
                index = (index % 9) * 3 + index / 9;
            int i = index / 3, e = index % 3;
            double pos = -1 + i * .25;
            x.append(pos);
            energy.append(foilEnergyRange[0] + e * (foilEnergyRange[1] - foilEnergyRange[0]) / 2.);
            for (int r = 0; r < height; ++r)
                for (int col = 0; col < width; ++col)
                    image[(frame * height + r) * width + col] =
                        quint64(std::nearbyint(spot(pos, r, col) + 2));
        }
        write(f, "entry/instrument/M_Stage_X_setpoint", x);
        write(f, "entry/instrument/M_DCM_B5Energy_setpoint", energy);
        write(f, "entry/data/D_SiC_Current1", QVector<double>(frames, 2));
        f.writeUInt("entry/data/D_LAMBDA", {quint64(frames), quint64(height), quint64(width)}, image.data());
    }
    {
        ImagingH5 f(directory + "/sample.nxs", true);
        f.writeText("entry/gui/plan_cmd",
                    QString(energyOuter
                                ? "('P.grid_scan', (['D.LAMBDA'], 'M.DCM_B5Energy', 1, 2, %1, "
                                  "'M.Stage_Z', -1, 0, 2, 'M.Stage_Y', 0, 1, 2), {})"
                                : "('P.grid_scan', (['D.LAMBDA'], 'M.Stage_Y', 0, 1, 2, "
                                  "'M.Stage_Z', -1, 0, 2, 'M.DCM_B5Energy', 1, 2, %1), {})")
                        .arg(sampleEnergyPoints));
        QVector<double> y, z, e, i0;
        const int cells = 4 * sampleEnergyPoints;
        const int frames = cells + (reorderedRepeat ? 1 : 0);
        std::vector<quint64> image(frames * height * width);
        for (int frame = 0; frame < frames; ++frame) {
            int index = reorderedRepeat ? (frame < cells ? cells - 1 - frame : 0) : frame;
            if (energyOuter)
                index = ((index % 2) * 2 + (index / 2) % 2) * sampleEnergyPoints + index / 4;
            const int position = index / sampleEnergyPoints, channel = index % sampleEnergyPoints;
            y.append(position / 2);
            z.append(position % 2 - 1);
            e.append(missing && index == cells - 1 ? 1. : 1. + double(channel) / (sampleEnergyPoints - 1));
            i0.append(channel % 2 ? 4 : 2);
            for (int r = 0; r < height; ++r)
                for (int col = 0; col < width; ++col)
                    image[(frame * height + r) * width + col] =
                        quint64(std::nearbyint(spot(0, r, col) * .1 / .05 * .2)) *
                        (channel % 2 ? (spectral ? 6 : 2) : 1) * (spatial ? 1 + position : 1);
        }
        write(f, "entry/instrument/M_Stage_Y_setpoint", y);
        write(f, "entry/instrument/M_Stage_Z_setpoint", z);
        write(f, "entry/instrument/M_DCM_B5Energy_setpoint", e);
        write(f, "entry/data/D_SiC_Current1", i0);
        f.writeUInt("entry/data/D_LAMBDA", {quint64(frames), quint64(height), quint64(width)}, image.data());
    }
    QFile roi(directory + "/roi.tsv");
    if (!roi.open(QIODevice::WriteOnly) ||
        roi.write("\troi_label\tx1\tx2\ty1\ty2\n1\told\t0\t19\t0\t47\n") < 0)
        ImageUtil::fail("Fixture ROI write failed");
    if (twoCrystals && roi.write("2\tsecond\t19\t38\t0\t47\n") < 0)
        ImageUtil::fail("Second fixture ROI write failed");
    roi.close();
    QJsonObject config{
        {"input", QJsonObject{{"calibration_nxs", "foil.nxs"},
                              {"sample_nxs", "sample.nxs"},
                              {"roi_table", "roi.tsv"}}},
        {"acquisition", QJsonObject{{"calibration_exposure_s", .5}, {"sample_exposure_s", .5}}},
        {"roi", QJsonObject{{"label_overrides", QJsonObject{{"1", "test"}}}}},
        {"grid", QJsonObject{{"x_step_mm", .1}}},
        {"reconstruction", QJsonObject{{"anchor_roi", "1_test"}}},
        {"tv", QJsonObject{{"lambda", 0.}, {"max_iterations", 1000}}},
        {"runtime", QJsonObject{{"output_dir", "output"}, {"chunk_frames", 4}}},
        {"output", QJsonObject{{"plots", false}}}};
    ImageUtil::saveJson(directory + "/config.json", config);
    return Imaging::config(config, directory);
}
} // namespace ImagingFixture
