#include "CrystalRoi.h"
#include "HdfFile.h"
#include <QCoreApplication>
#include <QFileInfo>
#include <QJsonArray>
#include <QJsonDocument>
#include <QJsonObject>
#include <QSet>
#include <QTemporaryDir>
#include <cmath>
#include <hdf5.h>
#include <iostream>

namespace {
int failures = 0;
void expect(bool condition, const char *description) {
    std::cout << (condition ? "PASS " : "FAIL ") << description << '\n';
    if (!condition)
        ++failures;
}
void fixtures(const QString &path) {
    auto file = H5Fcreate(path.toUtf8().constData(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    hsize_t dims[] = {2, 3, 4};
    auto space = H5Screate_simple(3, dims, nullptr);
    auto ds = H5Dcreate2(file, "cube", H5T_STD_I32LE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    int data[24];
    for (int i = 0; i < 24; ++i)
        data[i] = i;
    H5Dwrite(ds, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
    H5Dclose(ds);
    H5Sclose(space);
    space = H5Screate(H5S_SCALAR);
    ds = H5Dcreate2(file, "scalar", H5T_IEEE_F64LE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    double scalar = 42.5;
    H5Dwrite(ds, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &scalar);
    H5Dclose(ds);
    auto type = H5Tcopy(H5T_C_S1);
    H5Tset_size(type, H5T_VARIABLE);
    ds = H5Dcreate2(file, "text", type, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    const char *text = "synchrotron test";
    H5Dwrite(ds, type, H5S_ALL, H5S_ALL, H5P_DEFAULT, &text);
    H5Dclose(ds);
    H5Tclose(type);
    H5Sclose(space);
    hsize_t empty = 0;
    space = H5Screate_simple(1, &empty, nullptr);
    ds = H5Dcreate2(file, "empty", H5T_IEEE_F64LE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dclose(ds);
    H5Sclose(space);
    hsize_t large = 200001;
    space = H5Screate_simple(1, &large, nullptr);
    ds = H5Dcreate2(file, "large", H5T_IEEE_F64LE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dclose(ds);
    H5Sclose(space);
    hsize_t length = 3;
    space = H5Screate_simple(1, &length, nullptr);
    ds = H5Dcreate2(file, "nonfinite", H5T_IEEE_F64LE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    double nonfinite[] = {NAN, INFINITY, -INFINITY};
    H5Dwrite(ds, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, nonfinite);
    H5Dclose(ds);
    H5Sclose(space);
    H5Lcreate_soft("/absent", file, "broken", H5P_DEFAULT, H5P_DEFAULT);
    H5Lcreate_soft("/", file, "cycle", H5P_DEFAULT, H5P_DEFAULT);
    H5Fclose(file);
}
} // namespace
int main(int argc, char **argv) {
    QCoreApplication app(argc, argv);
    expect(HdfFile::formatFromFileName("14516_scan_20260917142117_1.nxs") == ScanFormat::Scan,
           "recognize scan basename");
    expect(HdfFile::formatFromFileName("15469_grid_scan_20261014231918_1.nxs") == ScanFormat::GridScan,
           "recognize grid_scan without confusing scan suffix");
    expect(HdfFile::formatFromFileName("C:/grid_scan/14516_SCAN_20260917142117_1.NXS") == ScanFormat::Scan,
           "classify case-insensitive basename independently of directory");
    for (const auto &name : {"test.nxs", "14516_other_scan_20260917142117_1.nxs", "14516_scan_20260917_1.nxs",
                             "14516_scan_20260917142117_1.nxs.bak"})
        expect(HdfFile::formatFromFileName(name) == ScanFormat::Generic,
               "unrecognized filename remains generic");
    QTemporaryDir dir;
    auto file = dir.filePath("回归测试.h5");
    fixtures(file);
    CrystalRois rois;
    QSet<QString> labels;
    for (const auto &roi : rois.items())
        labels.insert(roi.label);
    expect(rois.items().size() == 75 && labels.size() == 75 && rois.items()[59].detector == "D_LAMBDA" &&
               rois.items()[60].label == "VD-A1",
           "75 unique crystals, 60 LAMBDA and 15 MINIPIX");
    expect(rois.assign(0, QRect(500, 500, 16, 16), QSize(516, 516)),
           "ROI includes last detector pixel with exclusive upper bounds");
    expect(!rois.assign(1, QRect(500, 500, 17, 16), QSize(516, 516)) &&
               !rois.assign(1, QRect(10, 10, -3, 4), QSize(516, 516)),
           "reject out-of-bounds and reversed ROI");
    rois.invalidate(1);
    rois.invalidate(0);
    expect(rois.items()[0].rect == QRect(500, 500, 16, 16) && rois.count(RoiState::Invalid) == 2,
           "invalidate with or without geometry without losing crystal identity");
    QVector<QSize> sizes{QSize(516, 516), QSize(256, 256)};
    QString roiError;
    CrystalRois restored;
    auto serialized = rois.serialize("scan.nxs", sizes);
    expect(restored.restore(serialized, sizes, &roiError) && restored.count(RoiState::Invalid) == 2 &&
               restored.count(RoiState::Pending) == 73,
           "round-trip invalid and pending states");
    auto malformed = QJsonDocument::fromJson(serialized).object();
    auto crystals = malformed["crystals"].toArray();
    auto first = crystals[0].toObject();
    first["x2"] = 9999;
    crystals[0] = first;
    malformed["crystals"] = crystals;
    expect(!restored.restore(QJsonDocument(malformed).toJson(), sizes, &roiError) &&
               restored.items()[0].rect == QRect(500, 500, 16, 16),
           "reject invalid imported coordinates atomically");
    expect(!restored.restore(serialized, {QSize(256, 256), QSize(256, 256)}, &roiError),
           "reject detector geometry mismatch");
    rois.clear(0);
    expect(rois.items()[0].state == RoiState::Pending && rois.items()[0].rect.isEmpty(),
           "clear ROI for manual reselection");
    std::atomic_bool cancel{false};
    auto projection = HdfFile::sumFrames(file, "/cube", cancel);
    expect(projection.error.isEmpty() && projection.values.first() == 12 && projection.values.last() == 34 &&
               projection.sum == 276,
           "stream frame sum without loading full stack");
    cancel = true;
    expect(!HdfFile::sumFrames(file, "/cube", cancel).error.isEmpty(), "cancel frame sum");
    expect(HdfFile::inspect(file, "/cube").dims == QVector<quint64>({2, 3, 4}),
           "inspect multidimensional shape and Unicode path");
    auto children = HdfFile::children(file, "/");
    bool broken = false, cycle = false;
    for (const auto &child : children) {
        if (child.name == "broken")
            broken = child.type == "Unavailable link";
        if (child.name == "cycle")
            cycle = child.group;
    }
    expect(children.size() == 8 && broken && cycle, "enumerate broken and cyclic links without recursion");
    auto s = HdfFile::read({file, "/cube", {1, 0, 0}, 2, 1});
    expect(s.error.isEmpty() && s.width == 4 && s.height == 3 && s.values.first() == 12 &&
               s.values.last() == 23 && s.mean == 17.5,
           "read exact frame hyperslab and statistics");
    s = HdfFile::read({file, "/cube", {0, 1, 0}, 0, 2});
    expect(s.error.isEmpty() && s.width == 2 && s.height == 4 &&
               s.values == QVector<double>({4, 16, 5, 17, 6, 18, 7, 19}),
           "transpose arbitrary axes into row-major image");
    s = HdfFile::read({file, "/cube", {1, 1, 0}, 2, -1});
    expect(s.values == QVector<double>({16, 17, 18, 19}), "read vector from fixed dimensions");
    s = HdfFile::read({file, "/scalar", {}, -1, -1});
    expect(s.values.size() == 1 && s.values[0] == 42.5, "read scalar");
    s = HdfFile::read({file, "/text", {}, -1, -1});
    expect(s.text == "synchrotron test", "read variable string and reclaim memory");
    s = HdfFile::read({file, "/empty", {}, 0, -1});
    expect(s.values.isEmpty() && s.error.isEmpty(), "empty dataspace");
    s = HdfFile::read({file, "/nonfinite", {}, 0, -1});
    expect(s.finiteCount == 0 && std::isnan(s.mean), "nonfinite statistics");
    s = HdfFile::read({file, "/large", {}, 0, -1});
    expect(s.sampled && s.width <= 100000 && s.strides[0] == 3, "bounded vector preview");
    s = HdfFile::read({file, "/cube", {2, 0, 0}, 2, 1});
    expect(!s.error.isEmpty(), "reject out of range fixed index");
    s = HdfFile::read({file, "/cube", {}, 1, 1});
    expect(!s.error.isEmpty(), "reject overlapping axes");
    s = HdfFile::read({file, "/missing", {}, 0, -1});
    expect(!s.error.isEmpty(), "report missing dataset");
    if (argc > 1 && QFileInfo::exists(QString::fromLocal8Bit(argv[1]))) {
        QString real = QString::fromLocal8Bit(argv[1]);
        auto n = HdfFile::inspect(real, "/entry/data/D_LAMBDA");
        expect(n.dims == QVector<quint64>({7056, 516, 516}), "real NeXus LAMBDA dimensions");
        s = HdfFile::read({real, n.path, {0, 0, 0}, 2, 1});
        expect(s.error.isEmpty() && s.values.size() == 516 * 516 && s.minimum == 0 && s.maximum == 0,
               "real zero-count first frame");
        s = HdfFile::read({real, "/entry/data/D_MINIPIX", {7055, 0, 0}, 2, 1});
        expect(s.error.isEmpty() && s.values.size() == 256 * 256, "real MINIPIX final frame");
        s = HdfFile::read({real, "/entry/instrument/M_DCM_B5Energy", {}, 0, -1});
        expect(s.error.isEmpty() && s.values.size() == 7056 && s.minimum > 9 && s.maximum < 10,
               "real energy vector");
        s = HdfFile::read({real, "/entry/gui/plan_cmd", {}, -1, -1});
        expect(s.text.contains("grid_scan"), "real scan command string");
        s = HdfFile::read({real, "/entry/data/D_LAMBDA_image", {7055}, 0, -1});
        expect(s.error.isEmpty() && s.text.contains("1 element(s)") && s.text.contains("7055"),
               "bounded string array preview at final index");
        s = HdfFile::read({real, "/entry/processed/fit_ok", {}, 0, -1});
        expect(s.error.isEmpty() && s.values.size() == 60, "real boolean enum");
    }
    return failures ? 1 : 0;
}
