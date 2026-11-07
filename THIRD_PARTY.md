# Third-party components

XRI's own source code is released under the MIT License; see LICENSE.
Third-party components retain their original licenses. XRI uses the following libraries.

- Qt 6: Qt Core, Gui, Widgets, Concurrent, Svg and PrintSupport; LGPL v3 / GPL / commercial options, depending on the Qt distribution. Qt source and license information: https://www.qt.io/licensing/open-source-lgpl-obligations and https://download.qt.io/archive/qt/ . This project does not modify Qt.
- HDF5: The HDF Group's HDF5 license. https://www.hdfgroup.org/licenses/ . This project does not modify HDF5.
- Eigen 3.4.0: header-only native imaging linear algebra, compiled with `EIGEN_MPL2_ONLY`. Unmodified sources are in `third_party/eigen3/Eigen`, from https://gitlab.com/libeigen/eigen/-/tree/3.4.0 ; Mozilla Public License 2.0 and applicable upstream notices are in `third_party/eigen3/COPYING.*`. It is compiled into XRI, not a runtime DLL.
- Runtime packages copied by deployment may include zlib, zstd, ICU, HarfBuzz, libpng, PCRE2, FreeType and Microsoft Visual C++ runtime. Their licenses are included from the installed dependency packages when available.

The provided measurement file is user data and is not included in the distributable.
