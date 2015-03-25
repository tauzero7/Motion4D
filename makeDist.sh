#!/bin/bash

M4D_DIR=libMotion4D
TODAY=$(date +'%Y%m%d')

cd ..
zip -r libMotion4D_$TODAY.zip \
$M4D_DIR/CMakeLists.txt \
$M4D_DIR/common.pri \
$M4D_DIR/doxyfile \
$M4D_DIR/m4d.pro \
$M4D_DIR/m4d_sources.pri \
$M4D_DIR/AUTHORS \
$M4D_DIR/README \
$M4D_DIR/COPYING \
$M4D_DIR/NEWS \
$M4D_DIR/libMotion4D.sln \
$M4D_DIR/libMotion4D.vcxproj \
$M4D_DIR/libMotion4D.vcxproj.filters \
$M4D_DIR/m4d.props \
$M4D_DIR/src/m4dGlobalDefs.h \
$M4D_DIR/src/extra/*.cpp \
$M4D_DIR/src/extra/*.h \
$M4D_DIR/src/math/*.cpp \
$M4D_DIR/src/math/*.h \
$M4D_DIR/src/metric/m4dMetric.h \
$M4D_DIR/src/metric/m4dMetricAlcubierre.h \
$M4D_DIR/src/metric/m4dMetricAlcubierreSimple.h \
$M4D_DIR/src/metric/m4dMetricBarriolaVilenkin.h \
$M4D_DIR/src/metric/m4dMetricBertottiKasner.h \
$M4D_DIR/src/metric/m4dMetricBesselGravWaveCart.h \
$M4D_DIR/src/metric/m4dMetricBonnor.h \
$M4D_DIR/src/metric/m4dMetricChazyCurzonRot.h \
$M4D_DIR/src/metric/m4dMetricCosmicStringSchwarzschild.h \
$M4D_DIR/src/metric/m4dMetricCurzon.h \
$M4D_DIR/src/metric/m4dMetricDatabase.h \
$M4D_DIR/src/metric/m4dMetricDeSitterUnivConf.h \
$M4D_DIR/src/metric/m4dMetricDeSitterUniv.h \
$M4D_DIR/src/metric/m4dMetricEddFinkIn.h \
$M4D_DIR/src/metric/m4dMetricEinsteinRosenWaveWWB.h \
$M4D_DIR/src/metric/m4dMetricErezRosenVar.h \
$M4D_DIR/src/metric/m4dMetricErnst.h \
$M4D_DIR/src/metric/m4dMetricExtremeReissnerNordstromDihole.h \
$M4D_DIR/src/metric/m4dMetricFriedmanNonEmptyNull.h \
$M4D_DIR/src/metric/m4dMetricGlampedakis.h \
$M4D_DIR/src/metric/m4dMetricGoedel.h \
$M4D_DIR/src/metric/m4dMetricGoedelCart.h \
$M4D_DIR/src/metric/m4dMetricGoedelScaled.h \
$M4D_DIR/src/metric/m4dMetricGoedelScaledCart.h \
$M4D_DIR/src/metric/m4dMetricHalilsoyWave.h \
$M4D_DIR/src/metric/m4dMetricHartleThorneGB.h \
$M4D_DIR/src/metric/m4dMetricJaNeWi.h \
$M4D_DIR/src/metric/m4dMetricKasner.h \
$M4D_DIR/src/metric/m4dMetricKastorTraschen.h \
$M4D_DIR/src/metric/m4dMetricKerrBL.h \
$M4D_DIR/src/metric/m4dMetricKottler.h \
$M4D_DIR/src/metric/m4dMetricList.h \
$M4D_DIR/src/metric/m4dMetricMinkowski.h \
$M4D_DIR/src/metric/m4dMetricMinkowskiConformal.h \
$M4D_DIR/src/metric/m4dMetricMinkRotLattice.h \
$M4D_DIR/src/metric/m4dMetricMorrisThorne.h \
$M4D_DIR/src/metric/m4dMetricPainleveGullstrand.h \
$M4D_DIR/src/metric/m4dMetricPlaneGravWave.h \
$M4D_DIR/src/metric/m4dMetricPravda_C_Can.h \
$M4D_DIR/src/metric/m4dMetricPravda_C.h \
$M4D_DIR/src/metric/m4dMetricPTD_AI.h \
$M4D_DIR/src/metric/m4dMetricPTD_AII.h \
$M4D_DIR/src/metric/m4dMetricPTD_AIII.h \
$M4D_DIR/src/metric/m4dMetricPTD_BI.h \
$M4D_DIR/src/metric/m4dMetricPTD_BII.h \
$M4D_DIR/src/metric/m4dMetricPTD_BIII.h \
$M4D_DIR/src/metric/m4dMetricPTD_C.h \
$M4D_DIR/src/metric/m4dMetricReissnerNordstrom.h \
$M4D_DIR/src/metric/m4dMetricRotDihole.h \
$M4D_DIR/src/metric/m4dMetricSchwarzschild.h \
$M4D_DIR/src/metric/m4dMetricSchwarzschildCart.h \
$M4D_DIR/src/metric/m4dMetricSchwarzschildIsotropic.h \
$M4D_DIR/src/metric/m4dMetricSchwarzschildTortoise.h \
$M4D_DIR/src/metric/m4dMetricStraightSpinningString.h \
$M4D_DIR/src/metric/m4dMetricSultanaDyer.h \
$M4D_DIR/src/metric/m4dMetricTaubNUT.h \
$M4D_DIR/src/metric/m4dMetricTeoWHl.h \
$M4D_DIR/src/metric/m4dMetricTomimatsuSato.h \
$M4D_DIR/src/metric/m4dMetric.cpp \
$M4D_DIR/src/metric/m4dMetricAlcubierre.cpp \
$M4D_DIR/src/metric/m4dMetricAlcubierreSimple.cpp \
$M4D_DIR/src/metric/m4dMetricBarriolaVilenkin.cpp \
$M4D_DIR/src/metric/m4dMetricBertottiKasner.cpp \
$M4D_DIR/src/metric/m4dMetricBesselGravWaveCart.cpp \
$M4D_DIR/src/metric/m4dMetricBonnor.cpp \
$M4D_DIR/src/metric/m4dMetricChazyCurzonRot.cpp \
$M4D_DIR/src/metric/m4dMetricCosmicStringSchwarzschild.cpp \
$M4D_DIR/src/metric/m4dMetricCurzon.cpp \
$M4D_DIR/src/metric/m4dMetricDatabase.cpp \
$M4D_DIR/src/metric/m4dMetricDeSitterUnivConf.cpp \
$M4D_DIR/src/metric/m4dMetricDeSitterUniv.cpp \
$M4D_DIR/src/metric/m4dMetricEddFinkIn.cpp \
$M4D_DIR/src/metric/m4dMetricEinsteinRosenWaveWWB.cpp \
$M4D_DIR/src/metric/m4dMetricErezRosenVar.cpp \
$M4D_DIR/src/metric/m4dMetricErnst.cpp \
$M4D_DIR/src/metric/m4dMetricExtremeReissnerNordstromDihole.cpp \
$M4D_DIR/src/metric/m4dMetricFriedmanNonEmptyNull.cpp \
$M4D_DIR/src/metric/m4dMetricGlampedakis.cpp \
$M4D_DIR/src/metric/m4dMetricGoedel.cpp \
$M4D_DIR/src/metric/m4dMetricGoedelCart.cpp \
$M4D_DIR/src/metric/m4dMetricGoedelScaled.cpp \
$M4D_DIR/src/metric/m4dMetricGoedelScaledCart.cpp \
$M4D_DIR/src/metric/m4dMetricHalilsoyWave.cpp \
$M4D_DIR/src/metric/m4dMetricHartleThorneGB.cpp \
$M4D_DIR/src/metric/m4dMetricJaNeWi.cpp \
$M4D_DIR/src/metric/m4dMetricKasner.cpp \
$M4D_DIR/src/metric/m4dMetricKastorTraschen.cpp \
$M4D_DIR/src/metric/m4dMetricKerrBL.cpp \
$M4D_DIR/src/metric/m4dMetricKottler.cpp \
$M4D_DIR/src/metric/m4dMetricMinkowski.cpp \
$M4D_DIR/src/metric/m4dMetricMinkowskiConformal.cpp \
$M4D_DIR/src/metric/m4dMetricMinkRotLattice.cpp \
$M4D_DIR/src/metric/m4dMetricMorrisThorne.cpp \
$M4D_DIR/src/metric/m4dMetricPainleveGullstrand.cpp \
$M4D_DIR/src/metric/m4dMetricPlaneGravWave.cpp \
$M4D_DIR/src/metric/m4dMetricPravda_C_Can.cpp \
$M4D_DIR/src/metric/m4dMetricPravda_C.cpp \
$M4D_DIR/src/metric/m4dMetricPTD_AI.cpp \
$M4D_DIR/src/metric/m4dMetricPTD_AII.cpp \
$M4D_DIR/src/metric/m4dMetricPTD_AIII.cpp \
$M4D_DIR/src/metric/m4dMetricPTD_BI.cpp \
$M4D_DIR/src/metric/m4dMetricPTD_BII.cpp \
$M4D_DIR/src/metric/m4dMetricPTD_BIII.cpp \
$M4D_DIR/src/metric/m4dMetricPTD_C.cpp \
$M4D_DIR/src/metric/m4dMetricReissnerNordstrom.cpp \
$M4D_DIR/src/metric/m4dMetricRotDihole.cpp \
$M4D_DIR/src/metric/m4dMetricSchwarzschild.cpp \
$M4D_DIR/src/metric/m4dMetricSchwarzschildCart.cpp \
$M4D_DIR/src/metric/m4dMetricSchwarzschildIsotropic.cpp \
$M4D_DIR/src/metric/m4dMetricSchwarzschildTortoise.cpp \
$M4D_DIR/src/metric/m4dMetricStraightSpinningString.cpp \
$M4D_DIR/src/metric/m4dMetricSultanaDyer.cpp \
$M4D_DIR/src/metric/m4dMetricTaubNUT.cpp \
$M4D_DIR/src/metric/m4dMetricTeoWHl.cpp \
$M4D_DIR/src/metric/m4dMetricTomimatsuSato.cpp \
$M4D_DIR/src/motion/*.cpp \
$M4D_DIR/src/motion/*.h \
$M4D_DIR/test/testAll.bash \
$M4D_DIR/test/testAll.bat \
$M4D_DIR/test/README.test \
$M4D_DIR/test/tests.pro \
$M4D_DIR/test/calcGeodesic/calcGeodesic.cpp \
$M4D_DIR/test/calcGeodesic/calcGeodesic.pro \
$M4D_DIR/test/calcGeodesic/calcGeodesic.vcxproj \
$M4D_DIR/test/calcGeodesic/calcGeodesic.vcxproj.filters \
$M4D_DIR/test/calcGeodesic/kerr.ini \
$M4D_DIR/test/calcGeodesic/plotGeodesic.gnu \
$M4D_DIR/test/calcParallel/calcParallel.cpp \
$M4D_DIR/test/calcParallel/calcParallel.pro \
$M4D_DIR/test/calcParallel/calcParallel.vcxproj \
$M4D_DIR/test/calcParallel/calcParallel.vcxproj.filters \
$M4D_DIR/test/calcParallel/kerr_timelike.ini \
$M4D_DIR/test/testCircular/testCircular.cpp \
$M4D_DIR/test/testCircular/testCircular.pro \
$M4D_DIR/test/testCircular/testCircular.vcxproj \
$M4D_DIR/test/testCircular/testCircular.vcxproj.filters \
$M4D_DIR/test/testDatabase/testDatabase.cpp \
$M4D_DIR/test/testDatabase/testDatabase.pro \
$M4D_DIR/test/testDatabase/testDatabase.vcxproj \
$M4D_DIR/test/testDatabase/testDatabase.vcxproj.filters \
$M4D_DIR/test/testFWworldline/testFWworldline.cpp \
$M4D_DIR/test/testFWworldline/testFWworldline.pro \
$M4D_DIR/test/testFWworldline/testFWworldline.vcxproj \
$M4D_DIR/test/testFWworldline/testFWworldline.vcxproj.filters \
$M4D_DIR/test/testGeodesic/testGeodesic.cpp \
$M4D_DIR/test/testGeodesic/testGeodesic.pro \
$M4D_DIR/test/testGeodesic/testGeodesic.vcxproj \
$M4D_DIR/test/testGeodesic/testGeodesic.vcxproj.filters \
$M4D_DIR/test/testGeodesic/plotSchwGeods.gnu \
$M4D_DIR/test/testJacobi/testJacobi.cpp \
$M4D_DIR/test/testJacobi/testJacobi.pro \
$M4D_DIR/test/testJacobi/testJacobi.vcxproj \
$M4D_DIR/test/testJacobi/testJacobi.vcxproj.filters \
$M4D_DIR/test/testJacobi/plotJacobi.gnu \
$M4D_DIR/test/testJacobi/kerr_lens.ini \
$M4D_DIR/test/expected_results/m4dTestDB.txt \
$M4D_DIR/test/expected_results/kerr_geod.png \
$M4D_DIR/test/expected_results/points_par.dat \
$M4D_DIR/test/expected_results/points_circ.dat \
$M4D_DIR/test/expected_results/points_fw.dat \
$M4D_DIR/test/expected_results/schw_geods.png \
$M4D_DIR/test/expected_results/kerr_lens.png \
$M4D_DIR/test/expected_results/m4dTestCircular.txt \
$M4D_DIR/sonst/SymPy/alcubierreSimple.py \
$M4D_DIR/sonst/SymPy/barriolaVilenkin.py \
$M4D_DIR/sonst/SymPy/kerr.py \
$M4D_DIR/sonst/SymPy/m4d.py \
$M4D_DIR/sonst/SymPy/schwarzschild.py \
$M4D_DIR/sonst/Maxima/schwarzschild.mac \
$M4D_DIR/sonst/Maple/makeBarriolaVilenkin.mws \
$M4D_DIR/sonst/Maple/makeBonnor.mws \
$M4D_DIR/sonst/Maple/makeKerrBLbardeen.mws \
$M4D_DIR/sonst/Maple/makeKottler.mws  \
$M4D_DIR/sonst/Maple/makeMorrisThorne.mws \
$M4D_DIR/sonst/Maple/makePainleve.mws \
$M4D_DIR/sonst/Maple/makeReissner.mws

cd -
