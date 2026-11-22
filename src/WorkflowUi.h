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

inline QScrollArea *scrollBody(QWidget *body) {
    auto scroll = new QScrollArea;
    scroll->setWidgetResizable(true);
    scroll->setFrameShape(QFrame::NoFrame);
    scroll->setWidget(body);
    return scroll;
}

inline QWidget *cornerTools(std::initializer_list<QWidget *> controls) {
    auto bar = new QWidget;
    auto layout = new QHBoxLayout(bar);
    layout->setContentsMargins(0, 0, 0, 0);
    for (auto control : controls) {
        if (auto button = qobject_cast<QPushButton *>(control))
            button->setAutoDefault(false);
        layout->addWidget(control);
    }
    return bar;
}

inline QWidget *spectrumChooser(QComboBox *choice, const QString &suffix,
                                const QString &previousText, const QString &nextText) {
    auto bar = new QWidget;
    auto layout = new QHBoxLayout(bar);
    layout->setContentsMargins(0, 0, 0, 0);
    layout->setSpacing(4);
    layout->addWidget(choice, 1);
    const auto arrow = [&](const QString &name, const QString &symbol, const QString &text, int offset) {
        auto button = new QToolButton;
        button->setObjectName(name + suffix);
        button->setIcon(icon(symbol));
        button->setIconSize({16, 16});
        button->setToolTip(text);
        button->setAccessibleName(text);
        QObject::connect(button, &QToolButton::clicked, choice, [choice, offset] {
            const int next = choice->currentIndex() + offset;
            if (next >= 0 && next < choice->count())
                choice->setCurrentIndex(next);
        });
        layout->addWidget(button);
        return button;
    };
    auto previous = arrow("Previous", "chevron-left", previousText, -1);
    auto next = arrow("Next", "chevron-right", nextText, 1);
    const auto refresh = [choice, previous, next] {
        const int index = choice->currentIndex();
        previous->setEnabled(index > 0);
        next->setEnabled(index >= 0 && index + 1 < choice->count());
    };
    QObject::connect(choice, &QComboBox::currentIndexChanged, bar, refresh);
    QObject::connect(choice->model(), &QAbstractItemModel::rowsInserted, bar, refresh);
    QObject::connect(choice->model(), &QAbstractItemModel::rowsRemoved, bar, refresh);
    QObject::connect(choice->model(), &QAbstractItemModel::modelReset, bar, refresh);
    refresh();
    return bar;
}

inline QToolButton *sidePanelToggle(QWidget *panel, const QString &collapseText, const QString &expandText) {
    auto toggle = new QToolButton;
    toggle->setCheckable(true);
    toggle->setAutoRaise(true);
    toggle->setIconSize({12, 12});
    toggle->setFixedSize(18, 18);
    toggle->setIcon(icon("chevron-left"));
    toggle->setToolTip(collapseText);
    toggle->setAccessibleName(collapseText);
    QObject::connect(toggle, &QToolButton::toggled, panel, [panel, toggle, collapseText, expandText](bool collapsed) {
        panel->setVisible(!collapsed);
        toggle->setIcon(icon(collapsed ? "chevron-right" : "chevron-left"));
        const auto &label = collapsed ? expandText : collapseText;
        toggle->setToolTip(label);
        toggle->setAccessibleName(label);
    });
    return toggle;
}

inline QSplitter *analysisSplit(const QString &name, QWidget *controls, QWidget *preview) {
    auto splitter = new QSplitter(Qt::Horizontal);
    splitter->setObjectName(name);
    splitter->setChildrenCollapsible(false);
    splitter->setHandleWidth(6);
    auto sidebar = scrollBody(controls);
    sidebar->setMinimumWidth(240);
    splitter->addWidget(sidebar);

    auto content = new QWidget;
    auto layout = new QHBoxLayout(content);
    layout->setContentsMargins(0, 0, 0, 0);
    layout->setSpacing(4);
    auto toggle = sidePanelToggle(
        sidebar, QCoreApplication::translate("WorkflowUi", "收起参数面板"),
        QCoreApplication::translate("WorkflowUi", "展开参数面板"));
    toggle->setObjectName(name + "Toggle");
    layout->addWidget(toggle, 0, Qt::AlignTop);
    layout->addWidget(preview, 1);
    splitter->addWidget(content);
    splitter->handle(1)->setToolTip(QCoreApplication::translate("WorkflowUi", "拖动调整参数栏宽度"));
    splitter->setStyleSheet(QStringLiteral(
        "QSplitter#%1::handle { background: palette(midlight); }"
        "QSplitter#%1::handle:hover { background: palette(mid); }").arg(name));
    splitter->setStretchFactor(1, 1);
    splitter->setSizes({270, 930});
    return splitter;
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
