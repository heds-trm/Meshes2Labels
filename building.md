# Building
You will need a modern C++ compiler, the project uses the building tool [CMake](https://cmake.org
) (min. version 3.10).

## Prerequisites
* [ITK](https://itk.org/): image manipulation and IO.
* [VTK](https://vtk.org/) (version > 8.2): mesh manipulation and IO.
* [TCLAP](http://tclap.sourceforge.net/): command-line parser.

## CMake
Once `ITK` and `VTK` are built, use CMake to configure the project. You will have to specify the `ITK_BUILD` and `VTK_BUILD` folders, as well the path to the directory containing the `TCLAP` folder with headers (`TCLAP_INCLUDE_DIR` variable).

Also consider changing the `INSTALL_PATH` variable to install the `m2l` executable in your preferred location.