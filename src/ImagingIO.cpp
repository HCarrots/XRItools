#include "ImagingIO.h"
#include <QFile>
#include <QJsonArray>
#include <QJsonDocument>
#include <QMutexLocker>
#include <QSaveFile>
#include <algorithm>
#include <cmath>
#include <hdf5_hl.h>

namespace ImageUtil {
[[noreturn]] void fail(const QString &message) {
    throw std::runtime_error(message.toUtf8().constData());
}
void cancel(const std::atomic_bool &value) {
    if (value.load())
        fail(QStringLiteral("已取消成像计算；已生成的中间文件和日志保留。"));
}
QJsonArray json(const QVector<double> &v) {
    QJsonArray a;
    for (double x : v)
        a.append(std::isfinite(x) ? QJsonValue(x) : QJsonValue());
    return a;
}
QJsonArray json(const QVector<int> &v) {
    QJsonArray a;
    for (int x : v)
        a.append(x);
    return a;
}
QJsonArray json(const QStringList &v) {
    QJsonArray a;
    for (const auto &x : v)
        a.append(x);
    return a;
}
QVector<double> vector(const ImageMatrix &m) {
    return QVector<double>(m.data(), m.data() + m.size());
}
void saveJson(const QString &path, const QJsonObject &object) {
    QSaveFile f(path);
    auto bytes = QJsonDocument(object).toJson();
    if (!f.open(QIODevice::WriteOnly) || f.write(bytes) != bytes.size() || !f.commit())
        fail(QStringLiteral("写入失败：%1").arg(path));
}
double percentile(QVector<double> v, double q) {
    if (v.isEmpty())
        fail("Empty percentile input");
    std::sort(v.begin(), v.end());
    double index = (v.size() - 1) * q / 100;
    int i = int(index);
    return v[i] + (index - i) * (v[std::min(i + 1, int(v.size() - 1))] - v[i]);
}
} // namespace ImageUtil
namespace {
struct Handle {
    hid_t id;
    herr_t (*close)(hid_t);
    Handle(hid_t i, herr_t (*c)(hid_t)) : id(i), close(c) {
        if (i < 0)
            ImageUtil::fail("HDF5 imaging operation failed");
    }
    ~Handle() {
        close(id);
    }
    operator hid_t() const {
        return id;
    }
};
void check(herr_t x) {
    if (x < 0)
        ImageUtil::fail("HDF5 imaging read/write failed");
}
QVector<hsize_t> dimensions(hid_t space) {
    int rank = H5Sget_simple_extent_ndims(space);
    if (rank < 0)
        ImageUtil::fail("Invalid HDF5 dimensions");
    QVector<hsize_t> shape(rank);
    check(H5Sget_simple_extent_dims(space, shape.data(), nullptr));
    return shape;
}
} // namespace
ImagingH5::ImagingH5(const QString &file, bool create) {
    QMutexLocker lock(&HdfFile::ioMutex());
    H5Eset_auto2(H5E_DEFAULT, nullptr, nullptr);
    m_file = create ? H5Fcreate(file.toUtf8().constData(), H5F_ACC_EXCL, H5P_DEFAULT, H5P_DEFAULT)
                    : H5Fopen(file.toUtf8().constData(), H5F_ACC_RDONLY, H5P_DEFAULT);
    if (m_file < 0)
        ImageUtil::fail(QStringLiteral("无法%1 HDF5：%2").arg(create ? "新建" : "只读打开", file));
}
ImagingH5::~ImagingH5() {
    QMutexLocker lock(&HdfFile::ioMutex());
    if (m_file >= 0)
        H5Fclose(m_file);
}
bool ImagingH5::contains(const QString &path) const {
    QMutexLocker lock(&HdfFile::ioMutex());
    return H5Lexists(m_file, path.toUtf8().constData(), H5P_DEFAULT) > 0;
}
QVector<double> ImagingH5::read(const QString &path, QVector<quint64> *shape) const {
    QMutexLocker lock(&HdfFile::ioMutex());
    Handle ds(H5Dopen2(m_file, path.toUtf8().constData(), H5P_DEFAULT), H5Dclose),
        space(H5Dget_space(ds), H5Sclose);
    auto dims = dimensions(space);
    auto count = H5Sget_simple_extent_npoints(space);
    if (count < 0 || count > 500000000)
        ImageUtil::fail("Imaging array exceeds memory limit");
    if (shape) {
        shape->clear();
        for (auto d : dims)
            shape->append(d);
    }
    QVector<double> data(count);
    if (count)
        check(H5Dread(ds, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data.data()));
    return data;
}
QString ImagingH5::text(const QString &path) const {
    QMutexLocker lock(&HdfFile::ioMutex());
    Handle ds(H5Dopen2(m_file, path.toUtf8().constData(), H5P_DEFAULT), H5Dclose),
        type(H5Dget_type(ds), H5Tclose), space(H5Dget_space(ds), H5Sclose);
    if (H5Tget_class(type) != H5T_STRING || H5Sget_simple_extent_npoints(space) != 1)
        ImageUtil::fail("Expected scalar plan string");
    if (H5Tis_variable_str(type) > 0) {
        char *v = nullptr;
        check(H5Dread(ds, type, H5S_ALL, H5S_ALL, H5P_DEFAULT, &v));
        QString text = QString::fromUtf8(v ? v : "");
        H5free_memory(v);
        return text;
    }
    auto size = H5Tget_size(type);
    if (size > 1000000)
        ImageUtil::fail("Plan metadata too long");
    QByteArray bytes(qsizetype(size) + 1, 0);
    check(H5Dread(ds, type, H5S_ALL, H5S_ALL, H5P_DEFAULT, bytes.data()));
    return QString::fromUtf8(bytes.constData());
}
QStringList ImagingH5::strings(const QString &path, int maxCount) const {
    QMutexLocker lock(&HdfFile::ioMutex());
    Handle ds(H5Dopen2(m_file, path.toUtf8().constData(), H5P_DEFAULT), H5Dclose),
        type(H5Dget_type(ds), H5Tclose), space(H5Dget_space(ds), H5Sclose);
    auto shape = dimensions(space);
    if (shape.size() != 1 || !shape[0] || shape[0] > hsize_t(maxCount) || H5Tget_class(type) != H5T_STRING)
        ImageUtil::fail(QStringLiteral("需要一维 HDF5 字符串标签数组：%1").arg(path));
    const int count = int(shape[0]);
    QStringList result;
    if (H5Tis_variable_str(type) > 0) {
        QVector<char *> values(count, nullptr);
        auto status = H5Dread(ds, type, H5S_ALL, H5S_ALL, H5P_DEFAULT, values.data());
        if (status >= 0)
            for (auto value : values)
                result.append(QString::fromUtf8(value ? value : "").trimmed());
        H5Treclaim(type, space, H5P_DEFAULT, values.data());
        check(status);
    } else {
        auto width = H5Tget_size(type);
        if (!width || width > 4096)
            ImageUtil::fail(QStringLiteral("HDF5 晶体标签过长。"));
        QByteArray bytes(qsizetype(width) * count, 0);
        check(H5Dread(ds, type, H5S_ALL, H5S_ALL, H5P_DEFAULT, bytes.data()));
        for (int i = 0; i < count; ++i) {
            auto label = QByteArrayView(bytes.constData() + i * width, qsizetype(width));
            auto zero = label.indexOf('\0');
            result.append(QString::fromUtf8(zero < 0 ? label : label.first(zero)).trimmed());
        }
    }
    return result;
}
void ImagingH5::writeData(const QString &path, const QVector<quint64> &shape, const void *data, hid_t disk,
                          hid_t memory, bool compressed) {
    QMutexLocker lock(&HdfFile::ioMutex());
    QVector<hsize_t> dims;
    quint64 count = 1;
    for (auto d : shape) {
        dims.append(d);
        count *= d;
    }
    Handle space(shape.isEmpty() ? H5Screate(H5S_SCALAR)
                                 : H5Screate_simple(int(dims.size()), dims.constData(), nullptr),
                 H5Sclose),
        links(H5Pcreate(H5P_LINK_CREATE), H5Pclose), props(H5Pcreate(H5P_DATASET_CREATE), H5Pclose);
    check(H5Pset_create_intermediate_group(links, 1));
    if (compressed && count && !dims.isEmpty()) {
        auto chunk = dims;
        quint64 n = count;
        while (n > 131072) {
            auto max = std::max_element(chunk.begin(), chunk.end());
            n /= *max;
            *max = (*max + 1) / 2;
            n *= *max;
        }
        check(H5Pset_chunk(props, int(chunk.size()), chunk.constData()));
        check(H5Pset_shuffle(props));
        check(H5Pset_deflate(props, 4));
    }
    Handle ds(H5Dcreate2(m_file, path.toUtf8().constData(), disk, space, links, props, H5P_DEFAULT),
              H5Dclose);
    if (count)
        check(H5Dwrite(ds, memory, H5S_ALL, H5S_ALL, H5P_DEFAULT, data));
}
void ImagingH5::write(const QString &p, const QVector<quint64> &s, const double *v, bool compressed,
                      bool float32) {
    writeData(p, s, v, float32 ? H5T_IEEE_F32LE : H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, compressed);
}
void ImagingH5::writeUInt(const QString &p, const QVector<quint64> &s, const quint64 *v, bool compressed) {
    writeData(p, s, v, H5T_STD_U64LE, H5T_NATIVE_ULLONG, compressed);
}
void ImagingH5::writeStrings(const QString &p, const QStringList &values) {
    QMutexLocker lock(&HdfFile::ioMutex());
    QVector<QByteArray> utf8;
    std::vector<const char *> pointers;
    for (const auto &v : values)
        utf8.append(v.toUtf8());
    for (const auto &v : utf8)
        pointers.push_back(v.constData());
    hsize_t size = pointers.size();
    Handle type(H5Tcopy(H5T_C_S1), H5Tclose), space(H5Screate_simple(1, &size, nullptr), H5Sclose);
    check(H5Tset_size(type, H5T_VARIABLE));
    check(H5Tset_cset(type, H5T_CSET_UTF8));
    Handle ds(H5Dcreate2(m_file, p.toUtf8().constData(), type, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT),
              H5Dclose);
    if (size)
        check(H5Dwrite(ds, type, H5S_ALL, H5S_ALL, H5P_DEFAULT, pointers.data()));
}
void ImagingH5::writeText(const QString &path, const QString &value) {
    QMutexLocker lock(&HdfFile::ioMutex());
    Handle type(H5Tcopy(H5T_C_S1), H5Tclose), space(H5Screate(H5S_SCALAR), H5Sclose),
        links(H5Pcreate(H5P_LINK_CREATE), H5Pclose);
    check(H5Pset_create_intermediate_group(links, 1));
    check(H5Tset_size(type, H5T_VARIABLE));
    check(H5Tset_cset(type, H5T_CSET_UTF8));
    Handle ds(H5Dcreate2(m_file, path.toUtf8().constData(), type, space, links, H5P_DEFAULT, H5P_DEFAULT),
              H5Dclose);
    auto bytes = value.toUtf8();
    const char *p = bytes.constData();
    check(H5Dwrite(ds, type, H5S_ALL, H5S_ALL, H5P_DEFAULT, &p));
}
void ImagingH5::attributes(const QString &path, const QJsonObject &object) {
    QMutexLocker lock(&HdfFile::ioMutex());
    Handle obj(H5Oopen(m_file, path.toUtf8().constData(), H5P_DEFAULT), H5Oclose),
        space(H5Screate(H5S_SCALAR), H5Sclose);
    for (auto it = object.begin(); it != object.end(); ++it) {
        QByteArray name = it.key().toUtf8();
        if (it->isBool()) {
            Handle type(H5Tenum_create(H5T_STD_U8LE), H5Tclose);
            unsigned char no = 0, yes = 1, value = it->toBool() ? 1 : 0;
            check(H5Tenum_insert(type, "FALSE", &no));
            check(H5Tenum_insert(type, "TRUE", &yes));
            Handle attr(H5Acreate2(obj, name.constData(), type, space, H5P_DEFAULT, H5P_DEFAULT), H5Aclose);
            check(H5Awrite(attr, type, &value));
            continue;
        }
        if (it->isArray()) {
            bool numeric = true;
            for (auto value : it->toArray())
                numeric = numeric && value.isDouble();
            if (numeric) {
                QVector<double> values;
                for (auto value : it->toArray())
                    values.append(value.toDouble());
                hsize_t count = values.size();
                Handle dims(H5Screate_simple(1, &count, nullptr), H5Sclose);
                Handle attr(H5Acreate2(obj, name.constData(), H5T_IEEE_F64LE, dims, H5P_DEFAULT, H5P_DEFAULT),
                            H5Aclose);
                if (count)
                    check(H5Awrite(attr, H5T_NATIVE_DOUBLE, values.constData()));
                continue;
            }
        }
        if (it->isString() || it->isArray() || it->isObject() || it->isNull()) {
            QByteArray value = it->isString()   ? it->toString().toUtf8()
                               : it->isArray()  ? QJsonDocument(it->toArray()).toJson(QJsonDocument::Compact)
                               : it->isObject() ? QJsonDocument(it->toObject()).toJson(QJsonDocument::Compact)
                                                : QByteArray("null");
            Handle type(H5Tcopy(H5T_C_S1), H5Tclose);
            check(H5Tset_size(type, H5T_VARIABLE));
            check(H5Tset_cset(type, H5T_CSET_UTF8));
            Handle attr(H5Acreate2(obj, name.constData(), type, space, H5P_DEFAULT, H5P_DEFAULT), H5Aclose);
            const char *p = value.constData();
            check(H5Awrite(attr, type, &p));
        } else {
            double value = it->isBool() ? (it->toBool() ? 1 : 0) : it->toDouble();
            Handle attr(H5Acreate2(obj, name.constData(), H5T_IEEE_F64LE, space, H5P_DEFAULT, H5P_DEFAULT),
                        H5Aclose);
            check(H5Awrite(attr, H5T_NATIVE_DOUBLE, &value));
        }
    }
}
void ImagingH5::scales(const QString &volume, const QStringList &axes) {
    QMutexLocker lock(&HdfFile::ioMutex());
    Handle ds(H5Dopen2(m_file, volume.toUtf8().constData(), H5P_DEFAULT), H5Dclose);
    for (int i = 0; i < axes.size(); ++i) {
        auto name = axes[i].toUtf8();
        Handle axis(H5Dopen2(m_file, name.constData(), H5P_DEFAULT), H5Dclose);
        check(H5DSset_scale(axis, name.constData()));
        check(H5DSattach_scale(ds, axis, unsigned(i)));
        check(H5DSset_label(ds, unsigned(i), name.constData()));
    }
}
void ImagingH5::frames(const QString &path, int chunk, const std::atomic_bool &cancel,
                       const std::function<void(int, const quint64 *, int, int)> &visitor) const {
    QVector<hsize_t> dims;
    {
        QMutexLocker lock(&HdfFile::ioMutex());
        Handle ds(H5Dopen2(m_file, path.toUtf8().constData(), H5P_DEFAULT), H5Dclose),
            sp(H5Dget_space(ds), H5Sclose);
        dims = dimensions(sp);
    }
    if (dims.size() != 3 || !dims[0] || !dims[1] || !dims[2] || dims[1] > 4096 || dims[2] > 4096 || chunk < 1)
        ImageUtil::fail("Invalid imaging detector stack");
    qsizetype pixels = qsizetype(dims[1] * dims[2]);
    int actualChunk = std::min(chunk, int(std::max<qsizetype>(1, 16777216 / pixels)));
    std::vector<quint64> buffer(size_t(pixels * actualChunk));
    for (hsize_t start = 0; start < dims[0]; start += actualChunk) {
        ImageUtil::cancel(cancel);
        hsize_t n = std::min(hsize_t(actualChunk), dims[0] - start);
        {
            QMutexLocker lock(&HdfFile::ioMutex());
            Handle ds(H5Dopen2(m_file, path.toUtf8().constData(), H5P_DEFAULT), H5Dclose),
                sp(H5Dget_space(ds), H5Sclose);
            hsize_t counts[] = {n, dims[1], dims[2]}, offset[] = {start, 0, 0};
            check(H5Sselect_hyperslab(sp, H5S_SELECT_SET, offset, nullptr, counts, nullptr));
            Handle mem(H5Screate_simple(3, counts, nullptr), H5Sclose);
            check(H5Dread(ds, H5T_NATIVE_ULLONG, mem, sp, H5P_DEFAULT, buffer.data()));
        }
        for (hsize_t local = 0; local < n; ++local) {
            ImageUtil::cancel(cancel);
            visitor(int(start + local), buffer.data() + local * pixels, int(dims[2]), int(dims[1]));
        }
    }
}
