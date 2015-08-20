<!--
libMotion4D -- The Motion-4D-library.  
Copyright (c) 2009-2015  Thomas Mueller, Frank Grave  
This file is part of the m4d-library.  
-->

# libMotion4D - The Motion4D-library
Copyright (c) 2009-2015 Thomas Mueller, Frank Grave

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
["Catalogue of Spacetimes"](http://go.visus.uni-stuttgart.de/cos).


## Prerequisits:

* You need the GNU Scientific Library which is available [here](http://www.gnu.org/software/gsl/)  
* If you want to use lua scripting, then you also need [lua](http://www.lua.org/).


## Installation with cmake (linux):  
You have to have cmake 2.6 or higher available.

1. Build two new folders:  
      `mkdir -p build/debug`  
      `mkdir -p build/release`  

2. Enter the folder `build/debug` or `build/release`.

3. Call  `ccmake ../..`  and press 'c'.

4. Configure cmake build stuff:  
  If your compiler supports C++11, then set `CMAKE_CXX_FLAGS`  to  `-std=c++11`
    * Set the build type:  `Debug` or `Release`.
    * Set the installation path.
    * Adjust the path to the GNU Scientific Library and its library.  
      For example:  
         `GSL_DIR       /usr/local/gsl/1.15`  
         `GSL_LIB_DIR   /usr/local/gsl/1.15/lib`
    * If the test programs shall be build, then set  
         `COMPILE_TEST_PRGS   ON`
    * If you have 'lua' available, then set  
         `LUA_AVAILABLE   ON`  
         `LUA_DIR         /usr/local/lua/5.3`  
         `LUA_LIB_DIR     /usr/local/lua/5.3/lib`  

  Press 'c' and 'g' when you are finished.

5. `make`  
    `make install`

6. Run the test programs:  
    * Change to the `test` directory.  
    * `./testAll.bash [debug]`  `

  Note, if you have only compiled the debug version you have to add `debug` when
  running `testAll`!


## Installation with QT 4 (linux):
You have to have a qmake somewhere available.

1. Open "common.pri" in a standard text editor (e.g. gedit).  
   Adjust the path to the GNU Scientific Library, e.g.  
      `GSL_DIR     = /usr/local/gsl/1.15`  
      `GSL_LIB_DIR = $$GSL_DIR/lib`

2. Run `qmake`  
   You can also add:  `CONFIG+=debug`  
   and prefix:        `PREFIX=/usr/local/libm4d`  

3. `make`  
   `make install`

4. Compile the test programs:  
    * Change to the `test` directory.  
    * Run the same qmake command as in (2)  
    * `make`  
    * `make install`

5. Run the test programs:  
    * Change to the `test` directory.  
    * `./testAll.bash [debug]`

  Note, if you have only compiled the debug version you have to add `debug` when
  running `testAll`!



## Installation with VisualStudio 2012 (Windows):  
You need a Windows port of the GSL for your Visual Studio version. Store that in
the libMotion4D subfolder, e.g. `C:\Develop\libMotion4D\gsl`.

1. Open the solution file "libMotion4D.sln" in Visual Studio 2012.

2. Go to the property manager page -> libMotion4D -> m4d.  
    * Right click and select "Properties".
    * The user macro "GSL_DIR" points to the gsl library. If you have installed
      GSL in, e.g., `C:\Develop\libMotion4D\gsl-vs1012`, then you have to adjust
      this macro.
    * If you have installed GSL at a completely different location, you have to
      adjust the following properties:  
         - C/C++ -> "Additional Include Directories"  
         - Linker -> "Additional Library Directories  

3. Build the solution in `debug` or `release` configuration.

4. Run the test programs from a command shell:  
  * Change to the 'test' directory.
  * testAll.bat


## Add new metric:  
* Write `m4dMetric...` .h and .cpp file.
* Adjust `m4dMetricList.h`.
* Add metric to `m4dMetricDatabase::initializeMetric()`.
* Run `make` in TOP_DIR.

Take care that in `m4dMetricList.h` the order of the metrics of `stl_metric_names`
and `enum_metric` are the same !!!


## FILES (generated with 'tree'):
    .
    ├── AUTHORS  
    ├── CMakeLists.txt  
    ├── common.pri  
    ├── COPYING  
    ├── doxyfile  
    ├── libMotion4D.sln  
    ├── libMotion4D.vcxproj  
    ├── libMotion4D.vcxproj.filters  
    ├── lua  
    │   ├── examples  
    │   │   └── schwarzschild.lua  
    │   └── m4dlua_main.cpp  
    ├── m4d_lua.pri  
    ├── m4d.pro  
    ├── m4d.props  
    ├── m4d_sources.pri  
    ├── makeDist.sh  
    ├── NEWS  
    ├── README  
    ├── sonst  
    │   ├── m4dicon.png  
    │   ├── m4dicon.xcf
    │   ├── Maple
    │   │   ├── makeBarriolaVilenkin.mws
    │   │   ├── makeBarriolaVilenkin.txt
    │   │   ├── makeGravWave.mws
    │   │   ├── makeHartleThorne.mws
    │   │   ├── makeIsotropicCoords.mws
    │   │   ├── makeKerrBLbardeen.mws
    │   │   ├── makeKerrCarter.mws
    │   │   ├── makeKerrOrig.mws
    │   │   ├── makeKottler.mws
    │   │   ├── makeKruskalSzekeres.mws
    │   │   ├── makeKruskSzekConf.mws
    │   │   ├── makeMorrisThorne.mws
    │   │   ├── makePainleve.mws
    │   │   ├── makeReissner.mws
    │   │   ├── makeTeoWHl.mws
    │   │   └── reissner.mws
    │   ├── Maxima
    │   │   ├── schwarzschild.mac
    │   │   ├── schwarzschild.wxm
    │   │   ├── warp2.mac
    │   │   └── warp.mac
    │   └── SymPy
    │       ├── alcubierreSimple.py
    │       ├── barriolaVilenkin.py
    │       ├── extremeReissnerNordstromDihole.py
    │       ├── kerr.py
    │       ├── m4d.py
    │       ├── rotDihole.py
    │       ├── schwarzschild.py
    │       └── taubnut.py
    ├── src
    │   ├── extra
    │   │   ├── m4dObject.cpp
    │   │   ├── m4dObject.h
    │   │   ├── m4dPlatform.cpp
    │   │   ├── m4dPlatform.h
    │   │   ├── m4dUtilities.cpp
    │   │   └── m4dUtilities.h
    │   ├── lua
    │   │   ├── m4dlua.h
    │   │   ├── m4dlua_metric.cpp
    │   │   ├── m4dlua_metric.h
    │   │   ├── m4dlua_solver.cpp
    │   │   ├── m4dlua_solver.h
    │   │   ├── m4dlua_utils.cpp
    │   │   └── m4dlua_utils.h
    │   ├── m4dGlobalDefs.cpp
    │   ├── m4dGlobalDefs.h
    │   ├── math
    │   │   ├── m4dJacobi.cpp
    │   │   ├── m4dJacobi.h
    │   │   ├── Mat.cpp
    │   │   ├── Mat.h
    │   │   ├── TransCoordinates.cpp
    │   │   ├── TransCoordinates.h
    │   │   ├── TransfMat.cpp
    │   │   ├── TransfMat.h
    │   │   ├── VnD.cpp
    │   │   └── VnD.h
    │   ├── metric
    │   │   ├── m4dMetricAlcubierreAccel.cpp
    │   │   ├── m4dMetricAlcubierreAccel.h
    │   │   ├── m4dMetricAlcubierre.cpp
    │   │   ├── m4dMetricAlcubierre.h
    │   │   ├── m4dMetricAlcubierreSimple.cpp
    │   │   ├── m4dMetricAlcubierreSimple.h
    │   │   ├── m4dMetricBarriolaVilenkin.cpp
    │   │   ├── m4dMetricBarriolaVilenkin.h
    │   │   ├── m4dMetricBertottiKasner.cpp
    │   │   ├── m4dMetricBertottiKasner.h
    │   │   ├── m4dMetricBesselGravWaveCart.cpp
    │   │   ├── m4dMetricBesselGravWaveCart.h
    │   │   ├── m4dMetricBonnor.cpp
    │   │   ├── m4dMetricBonnor.h
    │   │   ├── m4dMetricChazyCurzonRot.cpp
    │   │   ├── m4dMetricChazyCurzonRot.h
    │   │   ├── m4dMetricCosmicStringSchwarzschild.cpp
    │   │   ├── m4dMetricCosmicStringSchwarzschild.h
    │   │   ├── m4dMetric.cpp
    │   │   ├── m4dMetricCurzon.cpp
    │   │   ├── m4dMetricCurzon.h
    │   │   ├── m4dMetricDatabase.cpp
    │   │   ├── m4dMetricDatabase.h
    │   │   ├── m4dMetricDeSitterUnivConf.cpp
    │   │   ├── m4dMetricDeSitterUnivConf.h
    │   │   ├── m4dMetricDeSitterUniv.cpp
    │   │   ├── m4dMetricDeSitterUniv.h
    │   │   ├── m4dMetricEddFinkIn.cpp
    │   │   ├── m4dMetricEddFinkIn.h
    │   │   ├── m4dMetricEinsteinRosenWaveWWB.cpp
    │   │   ├── m4dMetricEinsteinRosenWaveWWB.h
    │   │   ├── m4dMetricErezRosenVar.cpp
    │   │   ├── m4dMetricErezRosenVar.h
    │   │   ├── m4dMetricErnst.cpp
    │   │   ├── m4dMetricErnst.h
    │   │   ├── m4dMetricErnstSchwarzschild.cpp
    │   │   ├── m4dMetricErnstSchwarzschild.h
    │   │   ├── m4dMetricExtremeReissnerNordstromDihole.cpp
    │   │   ├── m4dMetricExtremeReissnerNordstromDihole.h
    │   │   ├── m4dMetricFriedmanNonEmptyNull.cpp
    │   │   ├── m4dMetricFriedmanNonEmptyNull.h
    │   │   ├── m4dMetricGlampedakis.cpp
    │   │   ├── m4dMetricGlampedakis.h
    │   │   ├── m4dMetricGoedelCart.cpp
    │   │   ├── m4dMetricGoedelCart.h
    │   │   ├── m4dMetricGoedel.cpp
    │   │   ├── m4dMetricGoedel.h
    │   │   ├── m4dMetricGoedelScaledCart.cpp
    │   │   ├── m4dMetricGoedelScaledCart.h
    │   │   ├── m4dMetricGoedelScaled.cpp
    │   │   ├── m4dMetricGoedelScaled.h
    │   │   ├── m4dMetric.h
    │   │   ├── m4dMetricHalilsoyWave.cpp
    │   │   ├── m4dMetricHalilsoyWave.h
    │   │   ├── m4dMetricHartleThorneGB.cpp
    │   │   ├── m4dMetricHartleThorneGB.h
    │   │   ├── m4dMetricJaNeWi.cpp
    │   │   ├── m4dMetricJaNeWi.h
    │   │   ├── m4dMetricKasner.cpp
    │   │   ├── m4dMetricKasner.h
    │   │   ├── m4dMetricKastorTraschen.cpp
    │   │   ├── m4dMetricKastorTraschen.h
    │   │   ├── m4dMetricKerrBL.cpp
    │   │   ├── m4dMetricKerrBL.h
    │   │   ├── m4dMetricKottler.cpp
    │   │   ├── m4dMetricKottler.h
    │   │   ├── m4dMetricList.h
    │   │   ├── m4dMetricMinkowskiConformal.cpp
    │   │   ├── m4dMetricMinkowskiConformal.h
    │   │   ├── m4dMetricMinkowski.cpp
    │   │   ├── m4dMetricMinkowski.h
    │   │   ├── m4dMetricMinkRotLattice.cpp
    │   │   ├── m4dMetricMinkRotLattice.h
    │   │   ├── m4dMetricMorrisThorne.cpp
    │   │   ├── m4dMetricMorrisThorne.h
    │   │   ├── m4dMetricPainleveGullstrand.cpp
    │   │   ├── m4dMetricPainleveGullstrand.h
    │   │   ├── m4dMetricPlaneGravWave.cpp
    │   │   ├── m4dMetricPlaneGravWave.h
    │   │   ├── m4dMetricPravda_C_Can.cpp
    │   │   ├── m4dMetricPravda_C_Can.h
    │   │   ├── m4dMetricPravda_C.cpp
    │   │   ├── m4dMetricPravda_C.h
    │   │   ├── m4dMetricPTD_AI.cpp
    │   │   ├── m4dMetricPTD_AI.h
    │   │   ├── m4dMetricPTD_AII.cpp
    │   │   ├── m4dMetricPTD_AII.h
    │   │   ├── m4dMetricPTD_AIII.cpp
    │   │   ├── m4dMetricPTD_AIII.h
    │   │   ├── m4dMetricPTD_BI.cpp
    │   │   ├── m4dMetricPTD_BI.h
    │   │   ├── m4dMetricPTD_BII.cpp
    │   │   ├── m4dMetricPTD_BII.h
    │   │   ├── m4dMetricPTD_BIII.cpp
    │   │   ├── m4dMetricPTD_BIII.h
    │   │   ├── m4dMetricPTD_C.cpp
    │   │   ├── m4dMetricPTD_C.h
    │   │   ├── m4dMetricReissnerNordstrom.cpp
    │   │   ├── m4dMetricReissnerNordstrom.h
    │   │   ├── m4dMetricRotDihole.cpp
    │   │   ├── m4dMetricRotDihole.h
    │   │   ├── m4dMetricSchwarzschildCart.cpp
    │   │   ├── m4dMetricSchwarzschildCart.h
    │   │   ├── m4dMetricSchwarzschild.cpp
    │   │   ├── m4dMetricSchwarzschild.h
    │   │   ├── m4dMetricSchwarzschildIsotropic.cpp
    │   │   ├── m4dMetricSchwarzschildIsotropic.h
    │   │   ├── m4dMetricSchwarzschildTortoise.cpp
    │   │   ├── m4dMetricSchwarzschildTortoise.h
    │   │   ├── m4dMetricStraightSpinningString.cpp
    │   │   ├── m4dMetricStraightSpinningString.h
    │   │   ├── m4dMetricSultanaDyer.cpp
    │   │   ├── m4dMetricSultanaDyer.h
    │   │   ├── m4dMetricTaubNUT.cpp
    │   │   ├── m4dMetricTaubNUT.h
    │   │   ├── m4dMetricTeoWHl.cpp
    │   │   ├── m4dMetricTeoWHl.h
    │   │   ├── m4dMetricTomimatsuSato.cpp
    │   │   └── m4dMetricTomimatsuSato.h
    │   └── motion
    │       ├── m4dFermiWalker.cpp
    │       ├── m4dFermiWalker.h
    │       ├── m4dGeodesicBS.cpp
    │       ├── m4dGeodesicBS.h
    │       ├── m4dGeodesic.cpp
    │       ├── m4dGeodesicDP54.cpp
    │       ├── m4dGeodesicDP54.h
    │       ├── m4dGeodesicDP65.cpp
    │       ├── m4dGeodesicDP65.h
    │       ├── m4dGeodesicGSL.cpp
    │       ├── m4dGeodesicGSL.h
    │       ├── m4dGeodesic.h
    │       ├── m4dGeodesicRK4.cpp
    │       ├── m4dGeodesicRK4.h
    │       ├── m4dMotion.cpp
    │       ├── m4dMotionDatabase.cpp
    │       ├── m4dMotionDatabase.h
    │       ├── m4dMotion.h
    │       └── m4dMotionList.h
    └── test
        ├── calcGeodesic
        │   ├── calcGeodesic.cpp
        │   ├── calcGeodesic.pro
        │   ├── calcGeodesic.vcxproj
        │   ├── calcGeodesic.vcxproj.filters
        │   ├── kerr.ini
        │   ├── plotGeodesic.gnu
        ├── calcParallel
        │   ├── calcParallel.cpp
        │   ├── calcParallel.pro
        │   ├── calcParallel.vcxproj
        │   ├── calcParallel.vcxproj.filters
        │   ├── kerr_timelike.ini
        ├── expected_results
        │   ├── kerr_geod.png
        │   ├── kerr_lens_m4d.dat
        │   ├── kerr_lens.png
        │   ├── m4dTestCircular.txt
        │   ├── m4dTestDB.txt
        │   ├── m4dTestFWwordline.txt
        │   ├── points_circ.dat
        │   ├── points_fw.dat
        │   ├── points_lens.dat
        │   ├── points_par.dat
        │   ├── schw_geods.png
        │   ├── schw_lens_triple.dat
        │   └── tidal.dat
        ├── README.test
        ├── testAll.bash
        ├── testAll.bat
        ├── testCircular
        │   ├── testCircular.cpp
        │   ├── testCircular.pro
        │   ├── testCircular.vcxproj
        │   └── testCircular.vcxproj.filters
        ├── testDatabase
        │   ├── testDatabase.cpp
        │   ├── testDatabase.pro
        │   ├── testDatabase.vcxproj
        │   └── testDatabase.vcxproj.filters
        ├── testFWworldline
        │   ├── testFWworldline.cpp
        │   ├── testFWworldline.pro
        │   ├── testFWworldline.vcxproj
        │   └── testFWworldline.vcxproj.filters
        ├── testGeodesic
        │   ├── plotSchwGeods.gnu
        │   ├── testGeodesic.cpp
        │   ├── testGeodesic.pro
        │   ├── testGeodesic.vcxproj
        │   └── testGeodesic.vcxproj.filters
        ├── testJacobi
        │   ├── kerr_lens.ini
        │   ├── plotJacobi.gnu
        │   ├── testJacobi.cpp
        │   ├── testJacobi.pro
        │   ├── testJacobi.vcxproj
        │   └── testJacobi.vcxproj.filters
        ├── tests.pro
        └── testTidalAccel
            └── testTidalAccel.cpp
