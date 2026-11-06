# XRI

C++17 / Qt 6 Widgets / HDF5 编写的 Windows 同步辐射数据查看器。

## 运行

打包后双击 `dist/XRI/XRI.exe`，使用“文件 → 打开”或拖入 `.nxs` / `.h5` 文件。
在项目目录执行以下命令，会自动打开提供的测试文件：

```powershell
.\scripts\run.ps1
```

也可以传入文件与初始数据集：

```powershell
.\dist\XRI\XRI.exe .\Data\15469_grid_scan_20261014231918_1.nxs --dataset /entry/data/D_LAMBDA --frame 1000
```

## 功能

根据完整文件名区分扫描格式，并在窗口标题、文件树根节点和数据摘要中显示：

| 文件名规则 | 识别类型 | 当前示例 |
| --- | --- | --- |
| `<编号>_scan_<14位时间戳>_<序号>.nxs` | scan · 一维扫描 | `14516_scan_20260917142117_1.nxs`（51 帧） |
| `<编号>_grid_scan_<14位时间戳>_<序号>.nxs` | grid_scan · 网格扫描 | `15469_grid_scan_20261014231918_1.nxs`（7056 帧） |

匹配不区分大小写，仅使用文件名，不受目录名影响。不符合规则的文件显示为“通用 HDF5 / NeXus”，仍可浏览。两种格式共用 HDF5 按需读取，维度与帧数以文件内容为准；打开新文件从第 0 帧开始。

- 原生 Windows Qt Widgets 菜单、文件对话框和可停靠面板；中文界面。打开和复位操作集中在顶部菜单。
- 按层展开 HDF5 分组，显示数据类型、维度和属性；筛选当前已加载节点。
- 标量、字符串、布尔、一维信号、二维图像和多维数组切片。
- 自由选择 X/Y 维度；底部滑块与帧输入框控制首个固定维度，支持播放。其他固定维度默认使用索引 0。
- 后台只读读取，快速切换时丢弃过期结果；探测器只读取当前帧。
- 白底学术图表、内向刻度、Viridis / Inferno / 灰度、等比例、保留符号的 log1p 显示。
- 滚轮缩放、双击复位、鼠标数值提示；显示当前切片的 min/max/mean/sum。
- 虚拟数值表，按需查看当前切片的数值。

索引从 **0** 开始。图像原点在左上，Y 向下；坐标使用数据集维度索引。

## 大文件和数据精度

HDF5 通过 hyperslab 按需读取，**不会整体加载约 4.36 GiB 的测试文件**。二维预览每个轴最多 1536 个采样点，一维预览最多 100000 点；超限时按固定步长采样，界面明确标注。统计对应预览中的采样值，并非整个原始数据集。字符串数组每次预览最后一维最多 64 项，可通过切片索引翻页。

数值预览与计算使用 float64，因此绝对值超过 2⁵³ 的 64 位整数可能舍入。目前不展开复合类型、引用、可变长数值数组和嵌套数组，只提供其元数据。外部链接或压缩过滤器缺失时报告读取错误。分组浏览限制为 32 层以防循环链接；本版本不推断扫描网格、物理轴、探测器几何或单位，也不包含拟合和校正流程。

## 编译

需要 MSVC x64 C++ 开发工具、CMake ≥ 3.24、Ninja、Qt ≥ 6.5（Widgets / Concurrent 开发包；回归测试还需要 Qt Test）和 HDF5 C 开发库。本机脚本默认使用已有的 `miniconda3/envs/xrsabre/Library`；匹配的 Qt 开发包是 `qtbase-devel=6.11.0`。运行时不依赖 Python。

```powershell
.\scripts\build.ps1 -Test
.\scripts\build.ps1 -QtRoot C:\Qt\6.8.3\msvc2022_64 -Hdf5Root C:\Libraries\hdf5 -Test
.\scripts\deploy.ps1
```

脚本只对当前进程设置编译和运行路径。也可以通过标准 CMake 指定 `CMAKE_PREFIX_PATH` 与 `HDF5_ROOT`。
Qt 部署使用官方 [windeployqt](https://doc.qt.io/qt-6/windows-deployment.html)，并额外解析 HDF5 与插件的 DLL 依赖。

## 验证

`ctest --test-dir build --output-on-failure` 运行 HDF5 回归测试，涵盖多维 hyperslab、轴交换、标量、字符串、空数据集、非有限值、采样边界、越界、坏链接和循环链接。若测试文件存在，还验证两组探测器、能量、字符串数组和布尔数据。

Qt 界面回归使用离屏平台，验证真实文件的帧切换、连续数据集选择、数据表页保持，无需占用桌面。

界面检查：

```powershell
.\dist\XRI\XRI.exe .\Data\15469_grid_scan_20261014231918_1.nxs --frame 1000 --smoke-test build\smoke
```

该模式加载数据、保存窗口截图，完成后退出，成功返回 0。读取实现参考 HDF5 官方 [dataspaces 与部分 I/O 文档](https://support.hdfgroup.org/documentation/hdf5/latest/_h5_s__u_g.html)。

## 代码结构

| 文件 | 职责 |
| --- | --- |
| `src/HdfFile.*` | HDF5 RAII、元数据、分块切片和统计 |
| `src/MainWindow.*` | 原生桌面 UI、异步读取调度、数据表 |
| `src/ScientificPlot.*` | Qt Painter 科研图表、图像交互 |
| `tests/HdfTests.cpp` | 合成边界数据及真实 NeXus 文件回归 |
| `scripts/` | 构建、运行和 Windows 部署 |
