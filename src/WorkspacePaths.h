#pragma once

#include <QString>

class QWidget;

namespace WorkspacePaths {
enum class Workflow { Roi, Spectroscopy, Imaging, Ixs };
enum class Kind { Config, Result };

QString root();
void setRoot(const QString &path);
QString directory(Workflow workflow, Kind kind);
// Returns an unused destination, or an empty string if no workspace is set.
QString newPath(Workflow workflow, Kind kind, const QString &source, const QString &name);
QString saveFile(QWidget *parent, Workflow workflow, Kind kind, const QString &source, const QString &title,
                 const QString &name, const QString &filter, const QString &manualSuggestion = {});
} // namespace WorkspacePaths
