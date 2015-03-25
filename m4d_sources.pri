#
#  Note that M4D_SRC_DIR must be defined outside !
#

M4D_GLOBAL_HEADER = $$M4D_SRC_DIR/m4dGlobalDefs.h
M4D_GLOBAL_SOURCE = $$M4D_SRC_DIR/m4dGlobalDefs.cpp

M4D_EXTRA_HEADERS = $$M4D_SRC_DIR/extra/m4dObject.h \
        $$M4D_SRC_DIR/extra/m4dUtilities.h \
        $$M4D_SRC_DIR/extra/m4dPlatform.h 
    
M4D_EXTRA_SOURCES =  $$M4D_SRC_DIR/extra/m4dObject.cpp \
        $$M4D_SRC_DIR/extra/m4dUtilities.cpp \
        $$M4D_SRC_DIR/extra/m4dPlatform.cpp    
                   
                    
M4D_MATH_HEADERS  = $$M4D_SRC_DIR/math/Mat.h \
        $$M4D_SRC_DIR/math/VnD.h \
        $$M4D_SRC_DIR/math/TransCoordinates.h \
        $$M4D_SRC_DIR/math/TransfMat.h

M4D_MATH_SOURCES  = $$M4D_SRC_DIR/math/Mat.cpp \
        $$M4D_SRC_DIR/math/VnD.cpp \
        $$M4D_SRC_DIR/math/TransCoordinates.cpp \
        $$M4D_SRC_DIR/math/TransfMat.cpp



M4D_METRIC_HEADERS = $$M4D_SRC_DIR/metric/m4dMetric.h \
        $$M4D_SRC_DIR/metric/m4dMetricAlcubierre.h \
        $$M4D_SRC_DIR/metric/m4dMetricAlcubierreSimple.h \
        $$M4D_SRC_DIR/metric/m4dMetricBarriolaVilenkin.h \
        $$M4D_SRC_DIR/metric/m4dMetricBertottiKasner.h \
        $$M4D_SRC_DIR/metric/m4dMetricBesselGravWaveCart.h \
        $$M4D_SRC_DIR/metric/m4dMetricBonnor.h \
        $$M4D_SRC_DIR/metric/m4dMetricChazyCurzonRot.h \
        $$M4D_SRC_DIR/metric/m4dMetricCosmicStringSchwarzschild.h \
        $$M4D_SRC_DIR/metric/m4dMetricCurzon.h \
        $$M4D_SRC_DIR/metric/m4dMetricDatabase.h \
        $$M4D_SRC_DIR/metric/m4dMetricDeSitterUniv.h \
        $$M4D_SRC_DIR/metric/m4dMetricDeSitterUnivConf.h \
        $$M4D_SRC_DIR/metric/m4dMetricEddFinkIn.h \
        $$M4D_SRC_DIR/metric/m4dMetricEinsteinRosenWaveWWB.h \
        $$M4D_SRC_DIR/metric/m4dMetricErezRosenVar.h \
        $$M4D_SRC_DIR/metric/m4dMetricErnst.h \
        $$M4D_SRC_DIR/metric/m4dMetricErnstSchwarzschild.h \
        $$M4D_SRC_DIR/metric/m4dMetricExtremeReissnerNordstromDihole.h \
        $$M4D_SRC_DIR/metric/m4dMetricFriedmanNonEmptyNull.h \
        $$M4D_SRC_DIR/metric/m4dMetricGlampedakis.h \
        $$M4D_SRC_DIR/metric/m4dMetricGoedelCart.h \
        $$M4D_SRC_DIR/metric/m4dMetricGoedel.h \
        $$M4D_SRC_DIR/metric/m4dMetricGoedelScaledCart.h \
        $$M4D_SRC_DIR/metric/m4dMetricGoedelScaled.h \
        $$M4D_SRC_DIR/metric/m4dMetricHalilsoyWave.h \
        $$M4D_SRC_DIR/metric/m4dMetricHartleThorneGB.h \
        $$M4D_SRC_DIR/metric/m4dMetricJaNeWi.h \
        $$M4D_SRC_DIR/metric/m4dMetricKasner.h \
        $$M4D_SRC_DIR/metric/m4dMetricKastorTraschen.h \
        $$M4D_SRC_DIR/metric/m4dMetricKerrBL.h \
        $$M4D_SRC_DIR/metric/m4dMetricKottler.h \
        $$M4D_SRC_DIR/metric/m4dMetricList.h \
        $$M4D_SRC_DIR/metric/m4dMetricMinkowski.h \
        $$M4D_SRC_DIR/metric/m4dMetricMinkowskiConformal.h \
        $$M4D_SRC_DIR/metric/m4dMetricMinkRotLattice.h \
        $$M4D_SRC_DIR/metric/m4dMetricMorrisThorne.h \
        $$M4D_SRC_DIR/metric/m4dMetricPainleveGullstrand.h \
        $$M4D_SRC_DIR/metric/m4dMetricPlaneGravWave.h \
        $$M4D_SRC_DIR/metric/m4dMetricPravda_C.h \
        $$M4D_SRC_DIR/metric/m4dMetricPravda_C_Can.h \
        $$M4D_SRC_DIR/metric/m4dMetricPTD_AI.h \
        $$M4D_SRC_DIR/metric/m4dMetricPTD_AII.h \
        $$M4D_SRC_DIR/metric/m4dMetricPTD_AIII.h \
        $$M4D_SRC_DIR/metric/m4dMetricPTD_BI.h \
        $$M4D_SRC_DIR/metric/m4dMetricPTD_BII.h \
        $$M4D_SRC_DIR/metric/m4dMetricPTD_BIII.h \
        $$M4D_SRC_DIR/metric/m4dMetricPTD_C.h \
        $$M4D_SRC_DIR/metric/m4dMetricReissnerNordstrom.h \
        $$M4D_SRC_DIR/metric/m4dMetricRotDihole.h \
        $$M4D_SRC_DIR/metric/m4dMetricSchwarzschild.h \
        $$M4D_SRC_DIR/metric/m4dMetricSchwarzschildCart.h \
        $$M4D_SRC_DIR/metric/m4dMetricSchwarzschildIsotropic.h \
        $$M4D_SRC_DIR/metric/m4dMetricSchwarzschildTortoise.h \
        $$M4D_SRC_DIR/metric/m4dMetricStraightSpinningString.h \
        $$M4D_SRC_DIR/metric/m4dMetricSultanaDyer.h \
        $$M4D_SRC_DIR/metric/m4dMetricTaubNUT.h \
        $$M4D_SRC_DIR/metric/m4dMetricTeoWHl.h \
	$$M4D_SRC_DIR/metric/m4dMetricTomimatsuSato.h

M4D_METRIC_SOURCES = $$M4D_SRC_DIR/metric/m4dMetric.cpp \
        $$M4D_SRC_DIR/metric/m4dMetricAlcubierre.cpp \
        $$M4D_SRC_DIR/metric/m4dMetricAlcubierreSimple.cpp \
        $$M4D_SRC_DIR/metric/m4dMetricBarriolaVilenkin.cpp \
        $$M4D_SRC_DIR/metric/m4dMetricBertottiKasner.cpp \
        $$M4D_SRC_DIR/metric/m4dMetricBesselGravWaveCart.cpp \
        $$M4D_SRC_DIR/metric/m4dMetricBonnor.cpp \
        $$M4D_SRC_DIR/metric/m4dMetricChazyCurzonRot.cpp \
        $$M4D_SRC_DIR/metric/m4dMetricCosmicStringSchwarzschild.cpp \
        $$M4D_SRC_DIR/metric/m4dMetricCurzon.cpp \
        $$M4D_SRC_DIR/metric/m4dMetricDatabase.cpp \
        $$M4D_SRC_DIR/metric/m4dMetricDeSitterUniv.cpp \
        $$M4D_SRC_DIR/metric/m4dMetricDeSitterUnivConf.cpp \
        $$M4D_SRC_DIR/metric/m4dMetricEddFinkIn.cpp \
        $$M4D_SRC_DIR/metric/m4dMetricEinsteinRosenWaveWWB.cpp \
        $$M4D_SRC_DIR/metric/m4dMetricErezRosenVar.cpp \
        $$M4D_SRC_DIR/metric/m4dMetricErnst.cpp \
        $$M4D_SRC_DIR/metric/m4dMetricErnstSchwarzschild.cpp \
        $$M4D_SRC_DIR/metric/m4dMetricExtremeReissnerNordstromDihole.cpp \
        $$M4D_SRC_DIR/metric/m4dMetricFriedmanNonEmptyNull.cpp \
        $$M4D_SRC_DIR/metric/m4dMetricGlampedakis.cpp \
        $$M4D_SRC_DIR/metric/m4dMetricGoedelCart.cpp \
        $$M4D_SRC_DIR/metric/m4dMetricGoedel.cpp \
        $$M4D_SRC_DIR/metric/m4dMetricGoedelScaledCart.cpp \
        $$M4D_SRC_DIR/metric/m4dMetricGoedelScaled.cpp \
        $$M4D_SRC_DIR/metric/m4dMetricHalilsoyWave.cpp \
        $$M4D_SRC_DIR/metric/m4dMetricHartleThorneGB.cpp \
        $$M4D_SRC_DIR/metric/m4dMetricJaNeWi.cpp \
        $$M4D_SRC_DIR/metric/m4dMetricKasner.cpp \
        $$M4D_SRC_DIR/metric/m4dMetricKastorTraschen.cpp \
        $$M4D_SRC_DIR/metric/m4dMetricKerrBL.cpp \
        $$M4D_SRC_DIR/metric/m4dMetricKottler.cpp \
        $$M4D_SRC_DIR/metric/m4dMetricMinkowski.cpp \
        $$M4D_SRC_DIR/metric/m4dMetricMinkowskiConformal.cpp \
        $$M4D_SRC_DIR/metric/m4dMetricMinkRotLattice.cpp \
        $$M4D_SRC_DIR/metric/m4dMetricMorrisThorne.cpp \
        $$M4D_SRC_DIR/metric/m4dMetricPainleveGullstrand.cpp \
        $$M4D_SRC_DIR/metric/m4dMetricPlaneGravWave.cpp \
        $$M4D_SRC_DIR/metric/m4dMetricPravda_C.cpp \
        $$M4D_SRC_DIR/metric/m4dMetricPravda_C_Can.cpp \
        $$M4D_SRC_DIR/metric/m4dMetricPTD_AI.cpp \
        $$M4D_SRC_DIR/metric/m4dMetricPTD_AII.cpp \
        $$M4D_SRC_DIR/metric/m4dMetricPTD_AIII.cpp \
        $$M4D_SRC_DIR/metric/m4dMetricPTD_BI.cpp \
        $$M4D_SRC_DIR/metric/m4dMetricPTD_BII.cpp \
        $$M4D_SRC_DIR/metric/m4dMetricPTD_BIII.cpp \
        $$M4D_SRC_DIR/metric/m4dMetricPTD_C.cpp \
        $$M4D_SRC_DIR/metric/m4dMetricReissnerNordstrom.cpp \
        $$M4D_SRC_DIR/metric/m4dMetricRotDihole.cpp \
        $$M4D_SRC_DIR/metric/m4dMetricSchwarzschild.cpp \
        $$M4D_SRC_DIR/metric/m4dMetricSchwarzschildCart.cpp \
        $$M4D_SRC_DIR/metric/m4dMetricSchwarzschildIsotropic.cpp \
        $$M4D_SRC_DIR/metric/m4dMetricSchwarzschildTortoise.cpp \
        $$M4D_SRC_DIR/metric/m4dMetricStraightSpinningString.cpp \
        $$M4D_SRC_DIR/metric/m4dMetricSultanaDyer.cpp \
        $$M4D_SRC_DIR/metric/m4dMetricTaubNUT.cpp \
        $$M4D_SRC_DIR/metric/m4dMetricTeoWHl.cpp \
	$$M4D_SRC_DIR/metric/m4dMetricTomimatsuSato.cpp


M4D_MOTION_HEADERS = $$M4D_SRC_DIR/motion/m4dFermiWalker.h \
        $$M4D_SRC_DIR/motion/m4dGeodesicGSL.h \
        $$M4D_SRC_DIR/motion/m4dGeodesic.h \
        $$M4D_SRC_DIR/motion/m4dGeodesicBS.h \
        $$M4D_SRC_DIR/motion/m4dGeodesicRK4.h \
        $$M4D_SRC_DIR/motion/m4dGeodesicDP54.h \
        $$M4D_SRC_DIR/motion/m4dGeodesicDP65.h \
        $$M4D_SRC_DIR/motion/m4dMotion.h \
        $$M4D_SRC_DIR/motion/m4dMotionDatabase.h \
        $$M4D_SRC_DIR/motion/m4dMotionList.h

M4D_MOTION_SOURCES = $$M4D_SRC_DIR/motion/m4dFermiWalker.cpp \
        $$M4D_SRC_DIR/motion/m4dGeodesicGSL.cpp \
        $$M4D_SRC_DIR/motion/m4dGeodesic.cpp \
        $$M4D_SRC_DIR/motion/m4dGeodesicBS.cpp \
        $$M4D_SRC_DIR/motion/m4dGeodesicRK4.cpp \
        $$M4D_SRC_DIR/motion/m4dGeodesicDP54.cpp \
        $$M4D_SRC_DIR/motion/m4dGeodesicDP65.cpp \
        $$M4D_SRC_DIR/motion/m4dMotion.cpp \
        $$M4D_SRC_DIR/motion/m4dMotionDatabase.cpp



M4D_HEADERS   = $$M4D_GLOBAL_HEADER \
                $$M4D_EXTRA_HEADERS \
                $$M4D_MATH_HEADERS \
                $$M4D_METRIC_HEADERS \
                $$M4D_MOTION_HEADERS

M4D_SOURCES   = $$M4D_GLOBAL_SOURCE \
                $$M4D_EXTRA_SOURCES \
                $$M4D_MATH_SOURCES \
                $$M4D_METRIC_SOURCES \
                $$M4D_MOTION_SOURCES



M4D_LUA_HEADERS = $$M4D_DIR/src/lua/m4dlua.h  \
                  $$M4D_DIR/src/lua/m4dlua_metric.h \
                  $$M4D_DIR/src/lua/m4dlua_solver.h \
                  $$M4D_DIR/src/lua/m4dlua_utils.h

M4D_LUA_SOURCES = $$M4D_DIR/src/lua/m4dlua_metric.cpp \
                  $$M4D_DIR/src/lua/m4dlua_solver.cpp \
                  $$M4D_DIR/src/lua/m4dlua_utils.cpp

