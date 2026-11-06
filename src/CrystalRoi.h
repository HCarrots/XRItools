#pragma once
#include <QByteArray>
#include <QRect>
#include <QString>
#include <QVector>

enum class RoiState { Pending, Valid, Invalid };
struct CrystalRoi {
    QString label, detector;
    QRect rect; // numpy convention: [x1,x2), [y1,y2)
    RoiState state = RoiState::Pending;
};
class CrystalRois {
  public:
    CrystalRois();
    const QVector<CrystalRoi> &items() const {
        return m_items;
    }
    bool assign(int row, const QRect &rect, const QSize &detectorSize);
    void invalidate(int row);
    void clear(int row);
    int count(RoiState state) const;
    QByteArray serialize(const QString &source, const QVector<QSize> &sizes) const;
    bool restore(const QByteArray &json, const QVector<QSize> &sizes, QString *error);

  private:
    QVector<CrystalRoi> m_items;
};
