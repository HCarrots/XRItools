#include "WorkspacePaths.h"
#include "ImagingIO.h"
#include <QCoreApplication>
#include <QDateTime>
#include <QDir>
#include <QFileDialog>
#include <QFileInfo>
#include <QMessageBox>
#include <QRegularExpression>
#include <QSettings>
#include <QUuid>
#include <exception>

namespace {
using namespace WorkspacePaths;
constexpr auto rootKey = "workspace/root";

QString subdirectory(const QString &root, Workflow workflow, Kind kind) {
    constexpr const char *folders[] = {"ROI", "spectroscopy", "imaging", "ixs"};
    return QDir(root).filePath(QString::fromLatin1(folders[int(workflow)]) +
                               (kind == Kind::Config ? "/config" : "/results"));
}

void ensureDirectory(const QString &path) {
    if (!QDir().mkpath(path))
        ImageUtil::fail(QCoreApplication::translate("WorkspacePaths", "无法创建工作区目录 %1。")
                            .arg(QDir::toNativeSeparators(path)));
}
} // namespace

QString WorkspacePaths::root() {
    return QSettings().value(rootKey).toString();
}

void WorkspacePaths::setRoot(const QString &path) {
    const QFileInfo location(path);
    const auto folder = location.canonicalFilePath();
    if (!location.isDir() || folder.isEmpty())
        ImageUtil::fail(QCoreApplication::translate("WorkspacePaths", "请选择一个可用的工作文件夹。"));
    static const QRegularExpression han(QStringLiteral("\\p{Han}"));
    if (folder.contains(han))
        ImageUtil::fail(QCoreApplication::translate("WorkspacePaths", "工作区路径不能包含中文，请选择其他文件夹。"));
    for (auto workflow : {Workflow::Roi, Workflow::Spectroscopy, Workflow::Imaging, Workflow::Ixs})
        for (auto kind : {Kind::Config, Kind::Result})
            ensureDirectory(subdirectory(folder, workflow, kind));
    QSettings settings;
    settings.setValue(rootKey, folder);
    settings.sync();
    if (settings.status() != QSettings::NoError)
        ImageUtil::fail(QCoreApplication::translate("WorkspacePaths", "无法保存工作区设置。"));
}

QString WorkspacePaths::directory(Workflow workflow, Kind kind) {
    const auto folder = root();
    return folder.isEmpty() ? QString{} : subdirectory(folder, workflow, kind);
}

QString WorkspacePaths::newPath(Workflow workflow, Kind kind, const QString &source, const QString &name) {
    const auto workspace = root();
    if (workspace.isEmpty())
        return {};
    if (!QFileInfo(workspace).isDir())
        ImageUtil::fail(QCoreApplication::translate("WorkspacePaths", "工作文件夹不可用，请重新指定工作区。"));
    const auto folder = subdirectory(workspace, workflow, kind);
    ensureDirectory(folder);
    static const QRegularExpression unsafeName(QStringLiteral("[^A-Za-z0-9_-]+"));
    static const QRegularExpression edgeUnderscores(QStringLiteral("^_+|_+$"));
    auto scan = QFileInfo(source.trimmed()).completeBaseName().replace(unsafeName, "_").left(52).remove(edgeUnderscores);
    if (!scan.isEmpty())
        scan += '_';
    const QFileInfo file(name);
    const auto suffix = file.completeSuffix().isEmpty() ? QString{} : "." + file.completeSuffix();
    const auto time = QDateTime::currentDateTime().toString("yyyyMMdd_HHmmss_zzz");
    QString path;
    do {
        const auto id = QUuid::createUuid().toString(QUuid::Id128).left(12);
        path = QDir(folder).filePath(scan + file.baseName() + '_' + time + '_' + id + suffix);
    } while (QFileInfo::exists(path));
    return path;
}

QString WorkspacePaths::saveFile(QWidget *parent, Workflow workflow, Kind kind, const QString &source,
                                 const QString &title, const QString &name, const QString &filter,
                                 const QString &manualSuggestion) {
    try {
        const auto path = newPath(workflow, kind, source, name);
        if (!path.isEmpty())
            return path;
    } catch (const std::exception &error) {
        QMessageBox::warning(parent, title, QString::fromUtf8(error.what()));
        return {};
    }
    const auto suggestion = manualSuggestion.isEmpty() ? name : manualSuggestion;
    QFileDialog picker(parent, title, suggestion, filter);
    picker.setAcceptMode(QFileDialog::AcceptSave);
    picker.setDefaultSuffix(name.endsWith(".roi.json") ? "roi.json" : QFileInfo(name).suffix());
    return picker.exec() == QDialog::Accepted ? picker.selectedFiles().constFirst() : QString{};
}
