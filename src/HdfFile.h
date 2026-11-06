#pragma once
#include <QPair>
#include <QString>
#include <QStringList>
#include <QVector>
#include <atomic>
#include <functional>
#include <stdexcept>

struct HdfNode {
    QString name, path, type, shape, link;
    bool group = false, numeric = false;
    QVector<quint64> dims;
    QVector<QPair<QString, QString>> attributes;
};

struct SliceRequest {
    QString file, path;
    QVector<quint64> indices;
    int xAxis = -1, yAxis = -1;
};

struct DataSlice {
    QString path, text, error;
    QVector<double> values;
    QVector<quint64> dims, starts, strides;
    int width = 0, height = 0, xAxis = -1, yAxis = -1;
    quint64 finiteCount = 0;
    double minimum = 0, maximum = 0, mean = 0, sum = 0;
    bool sampled = false;
};

// Each operation owns its HDF5 handles. A process-wide mutex also supports
// HDF5 builds without thread-safety; background reads never share handles.
enum class ScanFormat { Generic, Scan, GridScan };
class QMutex;

class HdfFile {
  public:
    static QMutex &ioMutex();
    static ScanFormat formatFromFileName(const QString &file);
    static HdfNode inspect(const QString &file, const QString &path);
    static QVector<HdfNode> children(const QString &file, const QString &path);
    static DataSlice read(const SliceRequest &request);
    static DataSlice sumFrames(const QString &file, const QString &path, const std::atomic_bool &cancel,
                               std::atomic<int> *progress = nullptr);
    static void visitFrames(const QString &file, const QString &path, const std::atomic_bool &cancel,
                            const std::function<void(int, const double *, int, int)> &visitor);
};
