#pragma once
#include <QSize>

class QImage;
class QIODevice;

// A looping GIF with a single shared 256-color palette. Frames fill the canvas.
class AnimatedGif {
  public:
    AnimatedGif(QIODevice &file, const QImage &firstFrame);
    bool append(const QImage &frame, int delayCentiseconds);
    bool finish();

  private:
    QIODevice &m_file;
    QSize m_size;
    bool m_ok = false;
};
