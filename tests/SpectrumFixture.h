#pragma once
#include "CrystalRoi.h"
#include <QDir>
#include <QFile>
#include <cmath>
#include <hdf5.h>
#include <stdexcept>

namespace SpectrumFixture {
inline void dataset(hid_t file, const QString &path, const QVector<hsize_t> &shape,
                    const QVector<double> &values) {
    auto lcpl = H5Pcreate(H5P_LINK_CREATE);
    H5Pset_create_intermediate_group(lcpl, 1);
    auto space = shape.isEmpty() ? H5Screate(H5S_SCALAR)
                                 : H5Screate_simple(int(shape.size()), shape.constData(), nullptr);
    auto ds =
        H5Dcreate2(file, path.toUtf8().constData(), H5T_IEEE_F64LE, space, lcpl, H5P_DEFAULT, H5P_DEFAULT);
    if (ds < 0 || H5Dwrite(ds, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, values.constData()) < 0)
        throw std::runtime_error("Fixture write failed");
    H5Dclose(ds);
    H5Sclose(space);
    H5Pclose(lcpl);
}
inline CrystalRois rois() {
    CrystalRois r;
    for (int i = 0; i < 75; ++i)
        r.invalidate(i);
    r.assign(0, {1, 1, 2, 2}, {8, 8});
    r.assign(1, {5, 5, 2, 2}, {8, 8});
    r.assign(60, {1, 1, 2, 2}, {6, 6});
    return r;
}
inline QString create(const QString &directory, bool elastic, bool descending = false, bool withI0 = true) {
    auto path = QDir(directory).filePath(elastic      ? "1_scan_20261106000000_1.nxs"
                                         : descending ? "3_scan_20261106000000_1.nxs"
                                                      : "2_scan_20261106000000_1.nxs");
    auto file = H5Fcreate(path.toUtf8().constData(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    if (file < 0)
        throw std::runtime_error("Fixture create failed");
    QVector<double> energy, i0, lambda(81 * 64, 0), mini(81 * 36, 0);
    for (int i = 0; i < 81; ++i) {
        int index = descending ? 80 - i : i;
        double e = (elastic ? 9.680 : 9.685) + index * .000125;
        energy.append(e);
        double flux = !elastic && index == 30 ? .2 : 2;
        i0.append(-flux);
        double signal =
            elastic ? 100 / (1 + std::pow((e - 9.685) * 1000 / .35, 2)) + 3 : 10 + .5 * (e - 9.685) * 1000;
        auto fill = [&](QVector<double> &image, int width, int x, int y, double factor) {
            int base = i * width * width + y * width + x;
            image[base] = flux * signal * factor;
            image[base + 1] = flux * signal * factor * .25;
            image[base + width] = flux * signal * factor * .1;
        };
        fill(lambda, 8, 1, 1, 1);
        fill(lambda, 8, 5, 5, 2);
        fill(mini, 6, 1, 1, 3);
    }
    // The provided real scan has only this fallback channel. Missing the first
    // channel must release HDF handles and the process mutex before retrying.
    dataset(file, "/entry/instrument/M_DCM_B5Energy_readback", {81}, energy);
    if (withI0)
        dataset(file, "/entry/data/D_SiC_Current1", {81}, i0);
    dataset(file, "/entry/data/D_LAMBDA", {81, 8, 8}, lambda);
    dataset(file, "/entry/data/D_MINIPIX", {81, 6, 6}, mini);
    for (const auto &module : {"VB", "HL", "VU", "HB", "VD"})
        dataset(file, QString("/entry/instrument/BL_ID33_XRS_m_%1_RBV").arg(module), {}, {90});
    H5Fclose(file);
    if (elastic) {
        QFile json(QDir(directory).filePath("test.roi.json"));
        if (!json.open(QIODevice::WriteOnly) || json.write(rois().serialize(path, {{8, 8}, {6, 6}})) < 0)
            throw std::runtime_error("Fixture ROI write failed");
    }
    return path;
}
} // namespace SpectrumFixture
