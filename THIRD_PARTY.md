# Third-party components

XRI's own source code is released under the MIT License; see LICENSE.
Third-party components retain their original licenses. XRI uses the following libraries and assets.

- Qt 6: Qt Core, Gui, Widgets, Concurrent, Svg and PrintSupport; LGPL v3 / GPL / commercial options, depending on the Qt distribution. Qt source and license information: https://www.qt.io/licensing/open-source-lgpl-obligations and https://download.qt.io/archive/qt/ . This project does not modify Qt.
- HDF5: The HDF Group's HDF5 license. https://www.hdfgroup.org/licenses/ . This project does not modify HDF5.
- Eigen 3.4.0: header-only native imaging linear algebra, compiled with `EIGEN_MPL2_ONLY`. Unmodified sources are in `third_party/eigen3/Eigen`, from https://gitlab.com/libeigen/eigen/-/tree/3.4.0 ; Mozilla Public License 2.0 and applicable upstream notices are in `third_party/eigen3/COPYING.*`. It is compiled into XRI, not a runtime DLL.
- Lucide 0.468.0: navigation and action SVG icons with a custom stroke color. See `assets/lucide/LICENSE` for the ISC license and `assets/lucide/LICENSE.feather` for the original Feather MIT license. [Lucide source](https://github.com/lucide-icons/lucide/tree/0.468.0/icons)
- XRI uses DABAX `ComptonProfiles.dat` as the built-in HF electron-shell table. The file credits the tabulations to Biggs, Mendelsohn and Mann (1975). The MIT license is in `third_party/dabax/LICENSE` and accompanies the binary package. [DABAX data and license](https://github.com/oasys-kit/DabaxFiles/tree/82de24d6347fea370671e2deace7ae310f940b00)
- The HF shell labels in `third_party/dabax/ComptonShellMap.tsv` use the Biggs electron configurations from xraylib. The BSD 3-Clause notice is in `third_party/dabax/xraylib-LICENSE`. [xraylib Biggs electron configurations](https://github.com/tschoonj/xraylib/blob/17c3e7287a750bdda3e78ba088ad54e56000e69e/data/comptonprofiles/biggs_electron_configuration.dat)
- Runtime packages copied by deployment may include zlib, zstd, ICU, HarfBuzz, libpng, PCRE2, FreeType and Microsoft Visual C++ runtime. Their licenses are included from the installed dependency packages when available.

The provided measurement file is user data and is not included in the distributable.
