#pragma once
#include <QColor>
#include <QRectF>
#include <QVector>
#include <QWidget>

struct SpectrumCurve {
    QString name;
    QVector<double> x, y;
    QColor color;
    bool dashed = false;
};
// Energy coordinates are explicit: irregular and descending scans retain their physical axes.
class SpectrumPlot : public QWidget {
    Q_OBJECT
  public:
    explicit SpectrumPlot(QWidget *parent = nullptr);
    void setFineGrid(bool enabled) {
        m_fineGrid = enabled;
        update();
    }
    void setCurves(const QVector<SpectrumCurve> &, const QString &title,
                   const QString &xLabel = "Energy (keV)", const QString &yLabel = "Intensity (a.u.)");
    const QVector<SpectrumCurve> &curves() const {
        return m_curves;
    }

  protected:
    void paintEvent(QPaintEvent *) override;
    void wheelEvent(QWheelEvent *) override;
    void mouseDoubleClickEvent(QMouseEvent *) override;
    void mouseMoveEvent(QMouseEvent *) override;
    void contextMenuEvent(QContextMenuEvent *) override;

  private:
    QRectF area() const;
    QVector<SpectrumCurve> m_curves;
    QString m_title, m_xLabel, m_yLabel;
    QRectF m_extent, m_view;
    bool m_grid = true, m_minorTicks = true, m_hasData = false;
    bool m_fineGrid = false;
};
