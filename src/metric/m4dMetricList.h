/**
 * @file    m4dMetricList.h
 * @author  Thomas Mueller
 *
 * @brief  Base class for each metric class.
 *
 *  This file is part of libMotion4D.
 */
#ifndef M4D_METRIC_LIST_H
#define M4D_METRIC_LIST_H

#include "m4dMetric.h"
#include "m4dMetricCosmicStringSchwarzschild.h"
#include "m4dMetricKerrBL.h"
#include "m4dMetricMinkowski.h"
#include "m4dMetricMorrisThorne.h"
#include "m4dMetricSchwarzschild.h"
#include "m4dMetricSchwarzschildIsotropic.h"

#ifdef ALL_METRICS_AVAILABLE
#include "m4dMetricAlcubierre.h"
#include "m4dMetricAlcubierreSimple.h"
#include "m4dMetricBarriolaVilenkin.h"
#include "m4dMetricBertottiKasner.h"
#include "m4dMetricBesselGravWaveCart.h"
#include "m4dMetricBonnor.h"
#include "m4dMetricChazyCurzonRot.h"
#include "m4dMetricCurzon.h"
#include "m4dMetricDeSitterUniv.h"
#include "m4dMetricDeSitterUnivConf.h"
#include "m4dMetricEddFinkIn.h"
#include "m4dMetricEinsteinRosenWaveWWB.h"
#include "m4dMetricErezRosenVar.h"
#include "m4dMetricErnst.h"
#include "m4dMetricExtremeReissnerNordstromDihole.h"
#include "m4dMetricFriedmanNonEmptyNull.h"
#include "m4dMetricGlampedakis.h"
#include "m4dMetricGoedel.h"
#include "m4dMetricGoedelCart.h"
#include "m4dMetricGoedelScaled.h"
#include "m4dMetricGoedelScaledCart.h"
#include "m4dMetricHalilsoyWave.h"
#include "m4dMetricHartleThorneGB.h"
#include "m4dMetricJaNeWi.h"
#include "m4dMetricKasner.h"
#include "m4dMetricKastorTraschen.h"
#include "m4dMetricKottler.h"
#include "m4dMetricKruskal.h"
#include "m4dMetricMinkRotLattice.h"
#include "m4dMetricMinkowskiConformal.h"
#include "m4dMetricPTD_AI.h"
#include "m4dMetricPTD_AII.h"
#include "m4dMetricPTD_AIII.h"
#include "m4dMetricPTD_BI.h"
#include "m4dMetricPTD_BII.h"
#include "m4dMetricPTD_BIII.h"
#include "m4dMetricPTD_C.h"
#include "m4dMetricPainleveGullstrand.h"
#include "m4dMetricPlaneGravWave.h"
#include "m4dMetricPravda_C.h"
#include "m4dMetricPravda_C_Can.h"
#include "m4dMetricReissnerNordstrom.h"
#include "m4dMetricRotDihole.h"
#include "m4dMetricSchwarzschildGravWave.h"
#include "m4dMetricSchwarzschildTortoise.h"
#include "m4dMetricSchwarzschildWT.h"
#include "m4dMetricStraightSpinningString.h"
#include "m4dMetricSultanaDyer.h"
#include "m4dMetricTaubNUT.h"
#include "m4dMetricTeoSimpleWH.h"
#include "m4dMetricTeoWHl.h"
#include "m4dMetricTomimatsuSato.h"
#include "m4dMetricVaidyaIncRad.h"

#endif // ALL_METRICS_AVAILABLE

namespace m4d {

class API_M4D_EXPORT MetricList
{
public:
#ifdef ALL_METRICS_AVAILABLE
    static const int NUM_METRICS = 60;
#else
    static const int NUM_METRICS = 7;
#endif // ALL_METRICS_AVAILABLE

/* --------------------------------------------------------
 *   List of all metrics currently implemented
 *
 *   When editing this list please take care of the ordering.
 *
 *   The names here must be equal to the names
 *   given in the constructor of the metric: mMetricName
 * -------------------------------------------------------- */
#ifdef ALL_METRICS_AVAILABLE
    static const char* stl_metric_names[];

    enum enum_metric {
        enum_metric_unknown = 0,
        enum_metric_minkowski = 1,
        enum_metric_minkowski_conf,
        enum_metric_minkowski_rotlattice,
        enum_metric_schwarzschild,
        enum_metric_schwarzschild_gravwave,
        enum_metric_schwarzschild_isotropic,
        enum_metric_schwarzschild_tortoise,
        enum_metric_schwarzschild_wt,
        enum_metric_eddfinkin,
        enum_metric_painleve,
        enum_metric_alcubierre,
        enum_metric_alcubierre_simple,
        enum_metric_barriola,
        enum_metric_bertottikasner,
        enum_metric_bessel_grav_wave_cart,
        enum_metric_bonnor,
        enum_metric_chazy_curzon_rot,
        enum_metric_cosmic_string_schwarzschild,
        enum_metric_curzon,
        enum_metric_einstein_rosen_wave_wwb,
        enum_metric_erezrosenvar,
        enum_metric_ernst,
        enum_metric_extreme_reissner_dihole,
        enum_metric_friedman_nonempty_null,
        enum_metric_glampedakis,
        enum_metric_geodel,
        enum_metric_geodel_cart,
        enum_metric_goedelscaled,
        enum_metric_goedelscaled_cart,
        enum_metric_halilsoy_wave,
        enum_metric_hartlethorne_gb,
        enum_metric_janewi,
        enum_metric_kasner,
        enum_metric_kastortraschen,
        enum_metric_kerrbl,
        enum_metric_kottler,
        enum_metric_kruskal,
        enum_metric_morristhorne,
        enum_metric_Petrov_TD_AI,
        enum_metric_Petrov_TD_AII,
        enum_metric_Petrov_TD_AIII,
        enum_metric_Petrov_TD_BI,
        enum_metric_Petrov_TD_BII,
        enum_metric_Petrov_TD_BIII,
        enum_metric_Petrov_TD_C,
        enum_metric_plane_grav_wave,
        enum_metric_Pravda_C_Metric,
        enum_metric_Pravda_C_Can,
        enum_metric_reissner,
        enum_metric_rotdihole,
        enum_metric_desitter_univ,
        enum_metric_desitter_univ_conf,
        enum_metric_spinning_string,
        enum_metric_sultana_dyer,
        enum_metric_taub_nut,
        enum_metric_teoSimpleWH,
        enum_metric_teowhl,
        enum_metric_tomimatsusato,
        enum_metric_vaidyaincrad
    };
#else
    static const char* stl_metric_names[];

    enum enum_metric {
        enum_metric_unknown = 0,
        enum_metric_minkowski = 1,
        enum_metric_schwarzschild,
        enum_metric_schwarzschild_isotropic,
        enum_metric_kerrbl,
        enum_metric_morristhorne,
        enum_metric_cosmic_string_schwarzschild
    };

#endif // ALL_METRICS_AVAILABLE
};

} // end namespace m4d

#endif
