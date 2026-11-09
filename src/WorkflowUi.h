#pragma once
#include <QtWidgets>

namespace WorkflowUi {
inline QIcon icon(const QString &name) {
    return QIcon(QStringLiteral(":/assets/lucide/") + name + QStringLiteral(".svg"));
}

inline QPushButton *button(const QString &iconName, const QString &label) {
    auto action = new QPushButton(icon(iconName), label);
    action->setIconSize(QSize(16, 16));
    return action;
}

inline void compact(QLayout *layout) {
    layout->setContentsMargins(6, 4, 6, 4);
    layout->setSpacing(4);
}

inline QToolButton *disclosureButton(const QString &title, QWidget *body, bool expanded = false) {
    auto toggle = new QToolButton;
    toggle->setText(title);
    toggle->setAutoRaise(true);
    toggle->setCheckable(true);
    toggle->setChecked(expanded);
    toggle->setToolButtonStyle(Qt::ToolButtonTextBesideIcon);
    const auto collapsedIcon = icon("chevron-right");
    const auto expandedIcon = icon("chevron-down");
    toggle->setIconSize(QSize(16, 16));
    toggle->setIcon(expanded ? expandedIcon : collapsedIcon);
    body->setVisible(expanded);
    QObject::connect(toggle, &QToolButton::toggled, body, [toggle, body, collapsedIcon, expandedIcon](bool on) {
        body->setVisible(on);
        toggle->setIcon(on ? expandedIcon : collapsedIcon);
    });
    return toggle;
}

inline QWidget *disclosure(const QString &title, QWidget *body, bool expanded = false) {
    auto group = new QWidget;
    group->setSizePolicy(QSizePolicy::Preferred, QSizePolicy::Maximum);
    auto layout = new QVBoxLayout(group);
    layout->setContentsMargins(0, 0, 0, 0);
    layout->setSpacing(4);
    layout->addWidget(body);
    layout->insertWidget(0, disclosureButton(title, body, expanded), 0, Qt::AlignLeft);
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
