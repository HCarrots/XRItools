#pragma once
#include "Isosurface.h"
#include <QColor>
#include <QImage>
#include <QLineF>
#include <QWidget>
#include <functional>

class QPainter;

class BondFigures : public QWidget {
    Q_OBJECT
  public:
    explicit BondFigures(QWidget *parent = nullptr);
    void setMontage(const QVector<QVector<double>> &axes, const QVector<double> &field,
                    const QVector<int> &zIndices, double minimum, double maximum, bool signedContrast,
                    const QString &scaleLabel);
    enum class SpectralStyle { Surface, Waterfall };
    void setSpectralSection(SpectralStyle style, const QVector<double> &energy, const QVector<double> &x,
                            const QVector<QVector<double>> &rows, const QString &energyLabel,
                            const QString &intensityLabel, const QString &sectionLabel);
    void setMesh(const QVector<QVector<double>> &axes, const QVector<SurfaceTriangle> &first,
                 const QVector<SurfaceTriangle> &second, const QString &caption,
                 const std::array<QColor, 2> &colors = {});
    bool hasMesh() const;
    void setDepthProfile(const QVector<double> &relativeSensitivity, int index, double lowThreshold);
    // Capture the mesh and view so frames can be rendered on a worker thread.
    std::function<QImage(double)> rotationFrames(const QSize &size) const;
    static QColor fieldColor(double fraction, bool signedContrast);

  signals:
    void depthSelected(int index);

  protected:
    void paintEvent(QPaintEvent *) override;
    void mousePressEvent(QMouseEvent *) override;
    void mouseMoveEvent(QMouseEvent *) override;
    void mouseReleaseEvent(QMouseEvent *) override;
    void mouseDoubleClickEvent(QMouseEvent *) override;
    void wheelEvent(QWheelEvent *) override;

  private:
    enum class Kind { Montage, Spectral, Mesh };
    void drawMontage(QPainter &);
    void drawSpectralSection(QPainter &);
    void drawMesh(QPainter &);
    Kind m_kind = Kind::Montage;
    SpectralStyle m_spectralStyle = SpectralStyle::Surface;
    QVector<QVector<double>> m_axes, m_rows;
    QVector<double> m_field;
    QVector<int> m_z;
    std::array<QVector<SurfaceTriangle>, 2> m_mesh;
    std::array<QColor, 2> m_meshColors;
    double m_min = 0, m_max = 1, m_yaw = -.7, m_pitch = .46, m_zoom = 1.;
    bool m_signed = false;
    QString m_label, m_caption, m_sectionLabel;
    QPointF m_lastMouse;
    QVector<double> m_depthProfile;
    int m_depthIndex = -1;
    double m_lowThreshold = 0;
    QLineF m_xAxisOnScreen;
    QPointF m_pressPosition;
    bool m_dragged = false;
    std::function<QVector3D(const QVector3D &)> m_depthProjector;
    QRectF depthPlotRect() const;
};
