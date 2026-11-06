#pragma once
#include <QtWidgets>

namespace WorkflowUi {
inline QWidget *disclosure(const QString &title, QWidget *body, bool expanded = false) {
    auto group = new QWidget;
    auto layout = new QVBoxLayout(group);
    layout->setContentsMargins(0, 0, 0, 0);
    auto toggle = new QToolButton;
    toggle->setText(title);
    toggle->setCheckable(true);
    toggle->setChecked(expanded);
    toggle->setToolButtonStyle(Qt::ToolButtonTextBesideIcon);
    toggle->setArrowType(expanded ? Qt::DownArrow : Qt::RightArrow);
    toggle->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Fixed);
    layout->addWidget(toggle);
    layout->addWidget(body);
    body->setVisible(expanded);
    QObject::connect(toggle, &QToolButton::toggled, body, [toggle, body](bool on) {
        body->setVisible(on);
        toggle->setArrowType(on ? Qt::DownArrow : Qt::RightArrow);
    });
    return group;
}
inline QString selectSource(QWidget *parent, const QStringList &files, const QString &current,
                            const QString &title,
                            const QString &filter = "NeXus / HDF5 (*.nxs *.h5 *.hdf5)") {
    if (!files.isEmpty()) {
        QMenu menu(parent);
        auto heading = menu.addAction(QObject::tr("已导入数据浏览器的文件"));
        heading->setEnabled(false);
        for (const auto &file : files) {
            auto action = menu.addAction(QFileInfo(file).fileName());
            action->setToolTip(file);
            action->setData(file);
        }
        menu.addSeparator();
        auto browse = menu.addAction(QObject::tr("选择磁盘上的其他文件…"));
        auto chosen = menu.exec(QCursor::pos());
        if (!chosen)
            return {};
        if (chosen != browse)
            return chosen->data().toString();
    }
    return QFileDialog::getOpenFileName(parent, title, QFileInfo(current).absolutePath(), filter);
}
} // namespace WorkflowUi
