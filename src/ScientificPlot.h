#pragma once
#include "HdfFile.h"
#include <QImage>
#include <QRectF>
#include <QWidget>

class ScientificPlot : public QWidget {
    Q_OBJECT
  public:
    explicit ScientificPlot(QWidget *parent = nullptr);
    void setData(const DataSlice &data);
    void setLogarithmic(bool enabled);
    void setPalette(int palette);
    void setAspectLocked(bool locked);
    void resetView();
    const DataSlice &data() const {
        return m_data;
    }
  signals:
    void cursorChanged(const QString &text);

  protected:
    void paintEvent(QPaintEvent *) override;
    void mouseMoveEvent(QMouseEvent *) override;
    void wheelEvent(QWheelEvent *) override;
    void mouseDoubleClickEvent(QMouseEvent *) override;

  private:
    void draw(QPainter &painter, const QRectF &bounds, bool interaction);
    void rebuildImage();
    QRectF plotRect(const QRectF &bounds) const;
    QPointF toData(const QPointF &point) const;
    QColor color(double fraction) const;
    DataSlice m_data;
    QImage m_image;
    bool m_log = false, m_aspect = true;
    int m_palette = 0;
    QRectF m_view;
    QPointF m_hover{-1, -1};
};
