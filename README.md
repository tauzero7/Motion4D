
# libMotion4D - The Motion4D-library

libMotion was first developed by Thomas Müller and Frank Grave at the 
Visualization Research Center, University of Stuttgart, Germany (2009-2015).

Current contact:
Dr. Thomas Müller  
Haus der Astronomie/Max Planck Institute for Astronomy  
69117 Heidelberg, Germany  
Email: tmueller [at] mpia.de

## Brief description:

The Motion4D-library solves the geodesic equation as well as the parallel- and
Fermi-Walker-transport in four-dimensional spacetimes numerically. Initial
conditions are given with respect to natural local tetrads which are adapted to
the symmetries or the coordinates of the spacetime. Beside some already
implemented metrics like the Schwarzschild and Kerr metric, the object oriented
structure of the library permits to implement other metrics or integrators in a
straightforward manner.


The Motion4D-library is published in Computer Physics Communications:

* T. Müller  
  __Motion4D-library extended (New Version Announcement)__  
  Computer Physics Communications 185, 2798-2799 (2014)  
  DOI: [10.1016/j.cpc.2014.04.020](dx.doi.org/10.1016/j.cpc.2014.04.020)

* T. Müller  
  __Motion4D-library extended (New Version Announcement)__  
  Computer Physics Communications 182, 1386-1388 (2011)  
  DOI: [10.1016/j.cpc.2011.02.009](dx.doi.org/10.1016/j.cpc.2011.02.009)  

* T. Müller, F. Grave  
  __An updated version of the Motion4D library (New Version Announcement)__  
  Computer Physics Communications 181, 703 (2010)  
  DOI: [10.1016/j.cpc.2009.10.021](dx.doi.org/10.1016/j.cpc.2009.10.021)  

* T. Müller, F. Grave  
  __Motion4D - A library for lightrays and timelike worldlines in the theory of relativity__  
  Computer Physics Communications 180, 2355-2360 (2009)  
  DOI: [10.1016/j.cpc.2009.07.014](dx.doi.org/10.1016/j.cpc.2009.07.014)  

Details to the implemented spacetimes can be found also in the
["Catalogue of Spacetimes"](http://arxiv.org/abs/0904.4184).


## Prerequisits:

* You need the GNU Scientific Library which is available [here](http://www.gnu.org/software/gsl/)  
  for Linux systems. A cmake-Version which can also be used for Windows/Visual Studio 2015 can 
  be found [here](https://github.com/ampl/gsl).
* If you want to use lua scripting, then you also need [lua](http://www.lua.org/).


## Installation with cmake (linux):  
You have to have cmake 2.6 or higher available.

1. Build two new folders:  
      `mkdir -p build/debug`  
      `mkdir -p build/release`  

2. Enter the folder `build/Debug` or `build/Release`.

3. Call  `ccmake ../..`  and press 'c'.

4. Configure cmake build stuff:  
  If your compiler supports C++11, then set `CMAKE_CXX_FLAGS`  to  `-std=c++11`
    * Set the build type:  `Debug` or `Release`.
    * Set the installation path.
    * Adjust the path to the GNU Scientific Library and its library.  
      For example:  
         `GSL_DIR       /usr/local/gsl/2.5`  
         `GSL_LIB_DIR   /usr/local/gsl/2.5/lib64`
    * If the test programs shall be build, then set  
         `COMPILE_TEST_PRGS   ON`
    * If you have 'lua' available, then set  
         `LUA_AVAILABLE   ON`  
         `LUA_DIR         /usr/local/lua/5.3.5`  
         `LUA_LIB_DIR     /usr/local/lua/5.3.5/lib`  

  Press 'c' and 'g' when you are finished.

5. `make`  
    `make install`

6. Run the test programs:  
    * Change to the `test` directory.  
    * `./testAll.bash [debug]`  `

  Note, if you have only compiled the debug version you have to add `debug` when
  running `testAll`!


## Install with cmake and Visual Studio 2015
You need GSL for Windows. There is a cmake-Version which can be found [here](https://github.com/ampl/gsl).  
If you also want to use LUA, download the sources from [here](https://www.lua.org/). 
You might need 7-zip or something else to extract the .tar.gz file.

Take care that you always choose the same configuration (Debug/Release), (Win32/x64) 
for all libraries.

### Compile this GSL version with cmake/vs 2015 (32 or 64bit):
1. Open "CMakeLists.txt" in cmake-gui and adjust the corresponding paths.
    - Use ".\build\Debug" or ".\build\Release"  as build paths
    - Set CMAKE_BUILD_TYPE  to either Debug or Release
    - Set CMAKE_CONFIGURATION_TYPES  to either Debug or Release
    - Set CMAKE_INSTALL_PREFIX  to e.g. "D:/local/ampl_gsl_vs2015"
    - Run 'Configure'
    - Run 'Generate'
2. Open ".\build\Debug\GSL.sln" or ".\build\Release\GSL.sln" with Visual Studio 2015
3. Run 'Build Solution'
4. Go to the "INSTALL" entry and run 'Build INSTALL' to install into "CMAKE_INSTALL_PREFIX".


### Compile libMotion4D
1. Open "CMakeLists.txt" in cmake-gui and adjust the corresponding paths.
   - Use .\build\Debug  or .\build\Release as build paths.
   - Set CMAKE_CONFIGURATION_TYPES to either Debug or Release
   - Set GSL_DIR to e.g. "D:\local\ampl_gsl_vs2015"
   - Set GSL_LIB_DIR to e.g. "D:\local\ampl_gsl_vs2015\lib"
   - Run 'Configure'
   - Run 'Generate'
2. Open ".\build\Debug\libmotion4d.sln" or ".\build\Release\libmotion4d.sln" with Visual Studio 2015
3. Run 'Build Solution'
4. Go to the "INSTALL" entry and run 'Build INSTALL' to install into "CMAKE_INSTALL_PREFIX".
   

## Add new metric
First, a few preliminary steps which are necessary to find your metric.
1. Open "m4dMetricList.h" and increase 'NUM_METRICS'.
   Append an 'enum' entry at the end of the 'enum_metric' list.
   Append an 'include' entry at the end of the '#include's.
2. Open "m4dMetricList.cpp" and append the name of the new metric at the end
   of the "stl_metric_names" list.
3. Open "m4dMetricDatabase.cpp" and append a 'case' entry for your new metric.

4. You can now start writing your own metric class. 
   For that, have a look e.g. into 'm4dMetricMinkowski' or copy the header and
   source file and rename the class.
   Take care that 'mMetricName' and the name of your metric in "stl_metric_names"
   are the same.
   
   There are a few abstract methods, which you have to overwrite (see m4dMetric.h):
    - calculateMetric
    - calculateChristoffels
    - localToCoord
    - coordToLocal
    - breakCondition

