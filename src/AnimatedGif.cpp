#include "AnimatedGif.h"
#include <QIODevice>
#include <QImage>
#include <unordered_map>

namespace {
void appendWord(QByteArray &bytes, int value) {
    bytes.append(char(value & 0xff));
    bytes.append(char((value >> 8) & 0xff));
}

// GIF89a Appendix F, https://www.w3.org/Graphics/GIF/spec-gif89a.txt
QByteArray compressedPixels(const QImage &frame) {
    constexpr int clearCode = 256, endCode = 257, tableSize = 4096;
    std::unordered_map<int, int> table;
    table.reserve(tableSize);
    int nextCode = endCode + 1, codeWidth = 9, pendingBits = 0;
    unsigned pending = 0;
    QByteArray bytes;
    auto writeCode = [&](int code) {
        pending |= unsigned(code) << pendingBits;
        pendingBits += codeWidth;
        while (pendingBits >= 8) {
            bytes.append(char(pending & 0xff));
            pending >>= 8;
            pendingBits -= 8;
        }
        // The decoder learns the last inserted entry upon reading this code.
        if (codeWidth < 12 && nextCode == (1 << codeWidth))
            ++codeWidth;
    };
    writeCode(clearCode);
    int prefix = frame.constScanLine(0)[0];
    for (int y = 0; y < frame.height(); ++y) {
        const auto *line = frame.constScanLine(y);
        for (int x = y == 0 ? 1 : 0; x < frame.width(); ++x) {
            const int pixel = line[x];
            const int key = (prefix << 8) | pixel;
            auto found = table.find(key);
            if (found != table.end()) {
                prefix = found->second;
                continue;
            }
            writeCode(prefix);
            if (nextCode < tableSize) {
                table.emplace(key, nextCode++);
            } else {
                writeCode(clearCode);
                table.clear();
                nextCode = endCode + 1;
                codeWidth = 9;
            }
            prefix = pixel;
        }
    }
    writeCode(prefix);
    writeCode(endCode);
    if (pendingBits > 0)
        bytes.append(char(pending));
    return bytes;
}
} // namespace

AnimatedGif::AnimatedGif(QIODevice &file, const QImage &firstFrame) : m_file(file), m_size(firstFrame.size()) {
    if (firstFrame.format() != QImage::Format_Indexed8 || firstFrame.colorCount() != 256 ||
        m_size.isEmpty() || m_size.width() > 65535 || m_size.height() > 65535)
        return;
    QByteArray header("GIF89a");
    appendWord(header, m_size.width());
    appendWord(header, m_size.height());
    header.append("\xf7\x00\x00", 3);
    for (const auto color : firstFrame.colorTable()) {
        header.append(char(qRed(color)));
        header.append(char(qGreen(color)));
        header.append(char(qBlue(color)));
    }
    header.append("\x21\xff\x0bNETSCAPE2.0\x03\x01\x00\x00\x00", 19);
    m_ok = m_file.write(header) == header.size();
}

bool AnimatedGif::append(const QImage &frame, int delayCentiseconds) {
    if (!m_ok || frame.size() != m_size || frame.format() != QImage::Format_Indexed8)
        return false;
    QByteArray block("\x21\xf9\x04\x04", 4);
    appendWord(block, delayCentiseconds);
    block.append("\x00\x00\x2c\x00\x00\x00\x00", 7);
    appendWord(block, m_size.width());
    appendWord(block, m_size.height());
    block.append("\x00\x08", 2);
    const auto compressed = compressedPixels(frame);
    for (qsizetype offset = 0; offset < compressed.size(); offset += 255) {
        const auto count = qMin(qsizetype(255), compressed.size() - offset);
        block.append(char(count));
        block.append(compressed.constData() + offset, count);
    }
    block.append('\0');
    m_ok = m_file.write(block) == block.size();
    return m_ok;
}

bool AnimatedGif::finish() {
    return m_ok && m_file.putChar(';');
}
