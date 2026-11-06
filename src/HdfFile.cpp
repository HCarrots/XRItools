#include "HdfFile.h"
#include <QFileInfo>
#include <QMutex>
#include <QMutexLocker>
#include <QRegularExpression>
#include <algorithm>
#include <cmath>
#include <hdf5.h>
#include <limits>
#include <numeric>
#include <vector>

namespace {
QMutex hdfMutex;
class Handle {
  public:
    hid_t id = -1;
    herr_t (*close)(hid_t) = nullptr;
    Handle(hid_t value, herr_t (*fn)(hid_t)) : id(value), close(fn) {
        if (id < 0)
            throw std::runtime_error(
                "HDF5 operation failed (missing object, unsupported filter, or unreadable file).");
    }
    ~Handle() {
        if (id >= 0)
            close(id);
    }
    Handle(const Handle &) = delete;
    operator hid_t() const {
        return id;
    }
};
void check(herr_t status) {
    if (status < 0)
        throw std::runtime_error(
            "HDF5 read failed. Check dataset type, compression plugins and external files.");
}
Handle openFile(const QString &file) {
    H5Eset_auto2(H5E_DEFAULT, nullptr, nullptr);
    return Handle(H5Fopen(file.toUtf8().constData(), H5F_ACC_RDONLY, H5P_DEFAULT), H5Fclose);
}
QString typeName(hid_t type) {
    const auto cls = H5Tget_class(type);
    if (cls == H5T_INTEGER)
        return QString("%1int%2")
            .arg(H5Tget_sign(type) == H5T_SGN_NONE ? "u" : "")
            .arg(H5Tget_size(type) * 8);
    if (cls == H5T_FLOAT)
        return QString("float%1").arg(H5Tget_size(type) * 8);
    if (cls == H5T_STRING)
        return H5Tis_variable_str(type) > 0 ? "UTF-8 / variable string" : "fixed string";
    if (cls == H5T_ENUM)
        return "enum / bool";
    if (cls == H5T_COMPOUND)
        return "compound";
    if (cls == H5T_REFERENCE)
        return "reference";
    if (cls == H5T_VLEN)
        return "variable-length array";
    if (cls == H5T_ARRAY)
        return "array";
    return "opaque / unsupported";
}
QVector<quint64> dimensions(hid_t space) {
    int rank = H5Sget_simple_extent_ndims(space);
    if (rank < 0)
        throw std::runtime_error("Cannot read dataspace dimensions.");
    std::vector<hsize_t> dims(size_t(rank), 0);
    check(H5Sget_simple_extent_dims(space, dims.data(), nullptr));
    QVector<quint64> result;
    for (auto d : dims)
        result.append(d);
    return result;
}
QString strings(hid_t object, hid_t type, hid_t space, bool attribute, hid_t fileSpace = H5S_ALL) {
    auto n = H5Sget_simple_extent_npoints(space);
    // Bound variable string metadata as well as numeric previews.
    if (n <= 0)
        return "(empty)";
    if (n > 64)
        return QString("String array: %1 elements; select a smaller slice.").arg(n);
    QStringList result;
    if (H5Tis_variable_str(type) > 0) {
        std::vector<char *> values(size_t(n), nullptr);
        check(attribute ? H5Aread(object, type, values.data())
                        : H5Dread(object, type, space, fileSpace == H5S_ALL ? space : fileSpace, H5P_DEFAULT,
                                  values.data()));
        for (auto v : values)
            result.append(v ? QString::fromUtf8(v) : QString());
        check(H5Dvlen_reclaim(type, space, H5P_DEFAULT, values.data()));
    } else {
        const auto size = H5Tget_size(type);
        if (size > 1024 * 1024)
            return "(string exceeds preview limit)";
        std::vector<char> values(size_t(n) * size + 1, 0);
        check(attribute ? H5Aread(object, type, values.data())
                        : H5Dread(object, type, space, fileSpace == H5S_ALL ? space : fileSpace, H5P_DEFAULT,
                                  values.data()));
        for (hssize_t i = 0; i < n; ++i)
            result.append(QString::fromUtf8(values.data() + size_t(i) * size, qsizetype(size))
                              .remove(QChar('\0'))
                              .trimmed());
    }
    return result.join('\n');
}
HdfNode node(hid_t file, const QString &path, bool attrs) {
    Handle obj(H5Oopen(file, path.toUtf8().constData(), H5P_DEFAULT), H5Oclose);
    HdfNode out;
    out.path = path;
    out.name = path == "/" ? "/" : path.section('/', -1);
    out.group = H5Iget_type(obj) == H5I_GROUP;
    if (out.group) {
        out.type = "Group";
    } else if (H5Iget_type(obj) == H5I_DATASET) {
        Handle type(H5Dget_type(obj), H5Tclose);
        Handle space(H5Dget_space(obj), H5Sclose);
        out.dims = dimensions(space);
        out.type = typeName(type);
        auto cls = H5Tget_class(type);
        out.numeric = cls == H5T_INTEGER || cls == H5T_FLOAT || cls == H5T_ENUM;
        QStringList dims;
        for (auto d : out.dims)
            dims << QString::number(d);
        out.shape = H5Sget_simple_extent_type(space) == H5S_NULL
                        ? "null"
                        : (dims.isEmpty() ? "scalar" : dims.join(" × "));
    } else
        out.type = "Named datatype";
    if (attrs) {
        const int count = H5Aget_num_attrs(obj);
        for (int i = 0; i < count; ++i) {
            Handle attr(
                H5Aopen_by_idx(obj, ".", H5_INDEX_NAME, H5_ITER_INC, hsize_t(i), H5P_DEFAULT, H5P_DEFAULT),
                H5Aclose);
            auto size = H5Aget_name(attr, 0, nullptr);
            QByteArray name(int(size) + 1, 0);
            H5Aget_name(attr, name.size(), name.data());
            Handle type(H5Aget_type(attr), H5Tclose);
            Handle space(H5Aget_space(attr), H5Sclose);
            QString value;
            auto cls = H5Tget_class(type);
            auto n = H5Sget_simple_extent_npoints(space);
            if (cls == H5T_STRING)
                value = strings(attr, type, space, true);
            else if ((cls == H5T_INTEGER || cls == H5T_FLOAT) && n > 0 && n <= 64) {
                std::vector<double> values(size_t(n), 0);
                check(H5Aread(attr, H5T_NATIVE_DOUBLE, values.data()));
                QStringList numbers;
                for (double v : values)
                    numbers << QString::number(v, 'g', 12);
                value = numbers.join(", ");
            } else
                value = QString("%1 (%2 elements)").arg(typeName(type)).arg(n);
            out.attributes.append({QString::fromUtf8(name.constData()), value});
        }
    }
    return out;
}
} // namespace

ScanFormat HdfFile::formatFromFileName(const QString &file) {
    // Match the complete basename so grid_scan cannot fall through to scan,
    // and directory names never affect classification.
    static const QRegularExpression pattern(
        QStringLiteral("^[0-9]+_(grid_scan|scan)_[0-9]{14}_[0-9]+\\.nxs$"),
        QRegularExpression::CaseInsensitiveOption);
    const auto match = pattern.match(QFileInfo(file).fileName());
    if (!match.hasMatch())
        return ScanFormat::Generic;
    return match.captured(1).compare("grid_scan", Qt::CaseInsensitive) == 0 ? ScanFormat::GridScan
                                                                            : ScanFormat::Scan;
}

HdfNode HdfFile::inspect(const QString &file, const QString &path) {
    QMutexLocker lock(&hdfMutex);
    auto f = openFile(file);
    return node(f, path, true);
}

QVector<HdfNode> HdfFile::children(const QString &file, const QString &path) {
    QMutexLocker lock(&hdfMutex);
    auto f = openFile(file);
    Handle group(H5Gopen2(f, path.toUtf8().constData(), H5P_DEFAULT), H5Gclose);
    H5G_info_t info{};
    check(H5Gget_info(group, &info));
    QVector<HdfNode> result;
    for (hsize_t i = 0; i < info.nlinks; ++i) {
        const auto len =
            H5Lget_name_by_idx(group, ".", H5_INDEX_NAME, H5_ITER_INC, i, nullptr, 0, H5P_DEFAULT);
        if (len < 0)
            continue;
        QByteArray name(int(len) + 1, 0);
        check(H5Lget_name_by_idx(group, ".", H5_INDEX_NAME, H5_ITER_INC, i, name.data(), name.size(),
                                 H5P_DEFAULT));
        const auto child = (path == "/" ? path : path + '/') + QString::fromUtf8(name.constData());
        try {
            result.append(node(f, child, false));
        } catch (const std::exception &) {
            HdfNode broken;
            broken.name = QString::fromUtf8(name.constData());
            broken.path = child;
            broken.type = "Unavailable link";
            result.append(broken);
        }
    }
    std::stable_sort(result.begin(), result.end(), [](const HdfNode &a, const HdfNode &b) {
        return a.group != b.group ? a.group : a.name < b.name;
    });
    return result;
}

DataSlice HdfFile::read(const SliceRequest &r) {
    QMutexLocker lock(&hdfMutex);
    DataSlice out;
    out.path = r.path;
    out.xAxis = r.xAxis;
    out.yAxis = r.yAxis;
    try {
        auto f = openFile(r.file);
        Handle dataset(H5Dopen2(f, r.path.toUtf8().constData(), H5P_DEFAULT), H5Dclose);
        Handle type(H5Dget_type(dataset), H5Tclose);
        Handle space(H5Dget_space(dataset), H5Sclose);
        out.dims = dimensions(space);
        const int rank = int(out.dims.size());
        auto cls = H5Tget_class(type);
        if (H5Sget_simple_extent_npoints(space) == 0) {
            out.text = "Empty dataset";
            return out;
        }
        if (cls == H5T_STRING) {
            if (rank == 0)
                out.text = strings(dataset, type, space, false);
            else {
                std::vector<hsize_t> offsets(size_t(rank), 0), counts(size_t(rank), 1);
                for (int i = 0; i < rank; ++i) {
                    offsets[i] = r.indices.value(i, 0);
                    if (offsets[i] >= out.dims[i])
                        throw std::runtime_error("String index outside dataset extent.");
                }
                counts.back() = std::min<hsize_t>(64, out.dims.back() - offsets.back());
                check(H5Sselect_hyperslab(space, H5S_SELECT_SET, offsets.data(), nullptr, counts.data(),
                                          nullptr));
                Handle memory(H5Screate_simple(rank, counts.data(), nullptr), H5Sclose);
                out.text = QString("String preview: %1 element(s), last dimension starts at %2\n\n")
                               .arg(counts.back())
                               .arg(offsets.back()) +
                           strings(dataset, type, memory, false, space);
            }
            return out;
        }
        if (cls != H5T_INTEGER && cls != H5T_FLOAT && cls != H5T_ENUM) {
            out.text = QString("Preview unavailable for %1. Metadata and attributes remain accessible.")
                           .arg(typeName(type));
            return out;
        }
        if ((rank > 0 && (r.xAxis < 0 || r.xAxis >= rank)) || r.yAxis >= rank || r.yAxis < -1 ||
            (r.yAxis >= 0 && r.xAxis == r.yAxis))
            throw std::runtime_error("Invalid display axes.");
        std::vector<hsize_t> start(size_t(rank), 0), count(size_t(rank), 1), stride(size_t(rank), 1);
        out.starts.fill(0, rank);
        out.strides.fill(1, rank);
        for (int i = 0; i < rank; ++i) {
            if (i == r.xAxis || i == r.yAxis) {
                const quint64 limit = r.yAxis >= 0 ? 1536 : 100000;
                stride[i] = std::max<quint64>(1, (out.dims[i] + limit - 1) / limit);
                count[i] = (out.dims[i] + stride[i] - 1) / stride[i];
                out.sampled |= stride[i] > 1;
            } else {
                start[i] = i < r.indices.size() ? r.indices[i] : 0;
                if (start[i] >= out.dims[i])
                    throw std::runtime_error("Slice index outside dataset extent.");
            }
            out.starts[i] = start[i];
            out.strides[i] = stride[i];
        }
        out.width = rank == 0 ? 1 : int(count[r.xAxis]);
        out.height = r.yAxis < 0 ? 1 : int(count[r.yAxis]);
        const auto n = size_t(out.width) * size_t(out.height);
        std::vector<double> raw(n, 0);
        if (rank > 0)
            check(H5Sselect_hyperslab(space, H5S_SELECT_SET, start.data(), stride.data(), count.data(),
                                      nullptr));
        hsize_t memorySize = n;
        Handle memory(H5Screate_simple(1, &memorySize, nullptr), H5Sclose);
        if (cls == H5T_ENUM) {
            // HDF5 does not define enum -> floating-point conversions.
            std::vector<long long> integers(n, 0);
            check(H5Dread(dataset, H5T_NATIVE_LLONG, memory, space, H5P_DEFAULT, integers.data()));
            std::copy(integers.begin(), integers.end(), raw.begin());
        } else
            check(H5Dread(dataset, H5T_NATIVE_DOUBLE, memory, space, H5P_DEFAULT, raw.data()));
        out.values.resize(qsizetype(n));
        // File order is C row-major; transpose only when the selected X axis
        // precedes Y. Fixed axes contribute a count of one.
        for (int y = 0; y < out.height; ++y)
            for (int x = 0; x < out.width; ++x) {
                const size_t source = r.yAxis >= 0 && r.xAxis < r.yAxis ? size_t(x) * out.height + y
                                                                        : size_t(y) * out.width + x;
                out.values[y * out.width + x] = raw[source];
            }
        out.minimum = std::numeric_limits<double>::infinity();
        out.maximum = -out.minimum;
        long double sum = 0;
        for (double v : out.values)
            if (std::isfinite(v)) {
                out.minimum = std::min(out.minimum, v);
                out.maximum = std::max(out.maximum, v);
                sum += v;
                ++out.finiteCount;
            }
        out.sum = double(sum);
        out.mean = out.finiteCount ? double(sum / out.finiteCount) : std::numeric_limits<double>::quiet_NaN();
        if (!out.finiteCount)
            out.minimum = out.maximum = std::numeric_limits<double>::quiet_NaN();
    } catch (const std::exception &e) {
        out.error = QString::fromUtf8(e.what());
    }
    return out;
}
