#pragma once
#include "Imaging.h"
#include <hdf5.h>

// All HDF5 calls share the viewer's process mutex, including close/attribute I/O.
// The frame callback runs outside the lock so it can write intermediate files.
class ImagingH5 {
  public:
    explicit ImagingH5(const QString &file, bool create = false);
    ~ImagingH5();
    ImagingH5(const ImagingH5 &) = delete;
    QVector<double> read(const QString &, QVector<quint64> *shape = nullptr) const;
    QString text(const QString &) const;
    QStringList strings(const QString &path, int maxCount) const;
    bool contains(const QString &) const;
    void write(const QString &, const QVector<quint64> &, const double *, bool compressed = true,
               bool float32 = false);
    void writeUInt(const QString &, const QVector<quint64> &, const quint64 *, bool compressed = true);
    void writeStrings(const QString &, const QStringList &);
    void writeText(const QString &, const QString &);
    void attributes(const QString &object, const QJsonObject &);
    void scales(const QString &volume, const QStringList &axes);
    void frames(const QString &, int chunk, const std::atomic_bool &,
                const std::function<void(int, const quint64 *, int, int)> &) const;

  private:
    hid_t m_file = -1;
    void writeData(const QString &, const QVector<quint64> &, const void *, hid_t diskType, hid_t memoryType,
                   bool);
};
namespace ImageUtil {
[[noreturn]] void fail(const QString &);
void cancel(const std::atomic_bool &);
QJsonArray json(const QVector<double> &);
QJsonArray json(const QVector<int> &);
QJsonArray json(const QStringList &);
QVector<double> vector(const ImageMatrix &);
void saveJson(const QString &, const QJsonObject &);
double percentile(QVector<double>, double);
} // namespace ImageUtil
