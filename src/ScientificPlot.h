#pragma once
#include "CrystalRoi.h"
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
    void setImageGrid(bool enabled) {
        m_imageGrid = enabled;
        update();
    }
    void resetView();
    void setRoiSelectionEnabled(bool enabled);
    void setRoiOverlays(const QVector<CrystalRoi> &rois, const QString &selected);
    void setProfilePosition(const QPoint &position);
    void setProfilePickingEnabled(bool enabled) {
        m_profilePicking = enabled;
        if (!enabled) {
            m_profileDrag = ProfileDrag::None;
            unsetCursor();
        }
    }
    void setPickAnywhere(bool enabled) {
        m_pickAnywhere = enabled;
    }
    QPoint profilePosition() const {
        return m_profilePosition;
    }
    void setAxisLabels(const QString &x, const QString &y);
    void setImageCoordinates(double xFirst, double xStep, double yFirst, double yStep);
    const DataSlice &data() const {
        return m_data;
    }
  signals:
    void cursorChanged(const QString &text);
    void rectangleSelected(const QRect &rect);
    void profilePointPicked(const QPoint &position);

  protected:
    void paintEvent(QPaintEvent *) override;
    void mouseMoveEvent(QMouseEvent *) override;
    void mousePressEvent(QMouseEvent *) override;
    void mouseReleaseEvent(QMouseEvent *) override;
    void keyPressEvent(QKeyEvent *) override;
    void wheelEvent(QWheelEvent *) override;
    void mouseDoubleClickEvent(QMouseEvent *) override;
    void contextMenuEvent(QContextMenuEvent *) override;

  private:
    void draw(QPainter &painter, const QRectF &bounds, bool interaction);
    void rebuildImage();
    QRectF plotRect(const QRectF &bounds) const;
    QPointF toData(const QPointF &point) const;
    QColor color(double fraction) const;
    enum class ProfileDrag { None, X, Y, Both };
    ProfileDrag profileHit(const QPointF &point) const;
    void moveProfile(const QPointF &point);
    ProfileDrag m_profileDrag = ProfileDrag::None;
    DataSlice m_data;
    QImage m_image;
    bool m_log = false, m_aspect = true;
    bool m_curveGrid = true, m_imageGrid = false, m_minorTicks = true;
    int m_palette = 0;
    QRectF m_view;
    QPointF m_hover{-1, -1};
    QVector<CrystalRoi> m_rois;
    QString m_selectedRoi;
    bool m_roiEnabled = false, m_dragging = false;
    QPointF m_dragStart, m_dragEnd, m_screenStart;
    QPoint m_profilePosition{-1, -1};
    bool m_profilePicking = false, m_pickAnywhere = false;
    QString m_xLabel, m_yLabel;
    bool m_physicalAxes = false;
    double m_xOrigin = 0, m_yOrigin = 0, m_xStep = 1, m_yStep = 1;
};
