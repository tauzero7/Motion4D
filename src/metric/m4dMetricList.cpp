/**
 * @file    m4dMetricList.cpp
 * @author  Thomas Mueller
 *
 *  This file is part of libMotion4D.
 */
#include "m4dMetricList.h"

namespace m4d {

#ifdef ALL_METRICS_AVAILABLE
const char* MetricList::stl_metric_names[] = { "unknown", "Minkowski", "MinkowskiConformal", "MinkowskiRotLattice",
    "Schwarzschild", "SchwarzschildGravWave", "SchwarzschildIsotropic", "SchwarzschildTortoise", "SchwarzschildWT",
    "EddFinkIn", "PainleveGullstrand", "AlcubierreWarp", "AlcubierreWarpSimple", "BarriolaVilenkin", "BertottiKasner",
    "BesselGravWaveCart", "Bonnor", "ChazyCurzonRot", "CosmicStringSchwarzschild", "Curzon", "EinsteinRosenWaveWWB",
    "ErezRosenVar", "Ernst", "ExtremeReissnerNordstromDihole", "FriedmanNonEmptyNull", "Glampedakis", "Goedel",
    "GoedelCart", "GoedelScaled", "GoedelScaledCart", "HalilsoyWave", "HartleThorneGB", "JanisNewmanWinicour", "Kasner",
    "KastorTraschen", "KerrBL", "Kottler", "MorrisThorne", "Petrov_Type_D_AI_ES", "Petrov_Type_D_AII_ES",
    "Petrov_Type_D_AIII_ES", "Petrov_Type_D_BI_ES", "Petrov_Type_D_BII_ES", "Petrov_Type_D_BIII_ES",
    "Petrov_Type_D_C_ES", "PlaneGravWave", "Pravda_C-Metric", "Pravda_C-Metric_Canonical_Coords", "ReissnerNordstrom",
    "RotDihole", "DeSitterUniv", "DeSitterUnivConformal", "StraightSpinningString", "SultanaDyerBlackhole", "TaubNUT",
    "TeoSimpleWH", "TeoWHl", "TomimatsuSato", "VaidyaIncRad" };
#else
const char* MetricList::stl_metric_names[] = { "unknown", "Minkowski", "Schwarzschild", "SchwarzschildIsotropic",
    "KerrBL", "MorrisThorne", "CosmicStringSchwarzschild" };
#endif

} // namespace m4d
