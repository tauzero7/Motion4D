// -------------------------------------------------------------------------------
/*
   m4dMetricDatabase.cpp

  Copyright (c) 2009-2014  Thomas Mueller, Frank Grave


   This file is part of the m4d-library.

   The m4d-library is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   The m4d-library is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with the m4d-library.  If not, see <http://www.gnu.org/licenses/>.

*/
// -------------------------------------------------------------------------------

#include "m4dMetricDatabase.h"

namespace m4d {

/*!
 */
MetricDatabase::MetricDatabase(bool printDatabase) {
    init();
    if (printDatabase) {
        printMetricList();
    }
}

MetricDatabase::~MetricDatabase() {
    if (!mMetricMap.empty()) {
        mMetricMap.clear();
    }
}


// *********************************** public methods ******************************
/*! Get the number of implemented metrics.
 *
 *  \return int: number of implemented metrics.
 */
int MetricDatabase::getNumMetrics() {
    return NUM_METRICS;
}

/*!  Initialize metric 'num' and return the pointer to it.
 *
 *  \param num : metric number.
 *  \return Metric: pointer to metric.
 */
Metric*
MetricDatabase::getMetric(enum_metric num) {
    return initializeMetric(num);
}

/*! Initialize metric 'mName' and return the pointer to it.
 *
 *  \param mName : name of metric.
 *  \return Metric: pointer to metric.
 */
Metric*
MetricDatabase::getMetric(std::string mName) {
    mMetricMapItr = mMetricMap.find(mName);
    if (mMetricMapItr == mMetricMap.end()) {
        fprintf(stderr, "Metric '%s' is not implemented!\n", mName.c_str());
        return NULL;
    }

    return initializeMetric(mMetricMapItr->second);
}

/*! Get the name of metric 'num'.
 *
 *  \param num : number of metric.
 *  \return string : name of metric.
 */
std::string MetricDatabase::getMetricName(enum_metric num) {
    if (int(num) >= 0 && int(num) < NUM_METRICS) {
        return stl_metric_names[num];
    }

    return std::string();
}

/*! Get the number of the 'mName' metric.
 *
 *  \param mName : name of metric.
 *  \return enum_metric : number of metric.
 */
enum_metric MetricDatabase::getMetricNr(std::string mName) {
    mMetricMapItr = mMetricMap.find(mName);
    if (mMetricMapItr == mMetricMap.end()) {
        fprintf(stderr, "Metric '%s' is not implemented!\n", mName.c_str());
        return enum_metric_unknown;
    }

    return mMetricMapItr->second;
}

/*!  Print list of all metrics.
 *
 *  \param fptr : pointer to file.
 */
void MetricDatabase::printMetricList(FILE* fptr) {
    fprintf(fptr, "      Metric name                # params\n");
    fprintf(fptr, "-----------------------------------------------------------\n");

    int      numParams;
    Metric*  metric;
    std::vector<std::string> paramNames;
    double   value;


    for (int i = 0; i < NUM_METRICS; i++) {
        paramNames.clear();
        if ((metric = getMetric(stl_metric_names[i])) != NULL) {
            numParams = metric->getNumParams();
            metric->getParamNames(paramNames);
        } else {
            numParams = 0;
        }
        fprintf(fptr, "%-30s       %d     (", stl_metric_names[i], numParams);

        if (metric != NULL) {
            for (unsigned int j = 0; j < paramNames.size(); j++) {
                metric->getParam(paramNames[j], value);
                fprintf(fptr, "%s=%f ", paramNames[j].c_str(), value);
            }
        }
        fprintf(fptr, ")\n");
    }
}


// ******************************** protected methods ******************************

/*!  Initialize database: store the metrics in a map.
 *
 */
void MetricDatabase::init() {
    for (int i = 0; i < NUM_METRICS; i++) {
        mMetricMap.insert(std::pair<std::string, enum_metric>(stl_metric_names[i], enum_metric(i)));
    }
}


// ---------------------------------------------------
//           Please sort by enum name !!
// ---------------------------------------------------
/*! Initialize metric: generate a new instance.
 *
 *  \param num : number of metric.
 *  \return Metric : pointer to metric.
 */
Metric*
MetricDatabase::initializeMetric(enum_metric  num) {
    Metric*  currMetric;

    switch (num) {
    default:
    case enum_metric_unknown:
        currMetric = NULL;
        break;
    case enum_metric_alcubierre:
        currMetric = new MetricAlcubierre;
        break;
    case enum_metric_alcubierre_simple:
        currMetric = new MetricAlcubierreSimple;
        break;
    case enum_metric_barriola:
        currMetric = new MetricBarriolaVilenkin;
        break;
    case enum_metric_bertottikasner:
        currMetric = new MetricBertottiKasner;
        break;
    case enum_metric_bessel_grav_wave_cart:
        currMetric = new MetricBesselGravWaveCart;
        break;
    case enum_metric_bonnor:
        currMetric = new MetricBonnor;
        break;
    case enum_metric_chazy_curzon_rot:
        currMetric = new MetricChazyCurzonRot;
        break;
    case enum_metric_cosmic_string_schwarzschild:
        currMetric = new MetricCosmicStringSchwarzschild;
        break;
    case enum_metric_curzon:
        currMetric = new MetricCurzon;
        break;
    case enum_metric_desitter_univ:
        currMetric = new MetricDeSitterUniv;
        break;
    case enum_metric_desitter_univ_conf:
        currMetric = new MetricDeSitterUnivConf;
        break;
    case enum_metric_eddfinkin:
        currMetric = new MetricEddFinkIn;
        break;
    case enum_metric_einstein_rosen_wave_wwb:
        currMetric = new MetricEinsteinRosenWaveWWB;
        break;
    case enum_metric_erezrosenvar:
        currMetric = new MetricErezRosenVar;
        break;
    case enum_metric_ernst:
        currMetric = new MetricErnst;
        break;
    case enum_metric_extreme_reissner_dihole:
        currMetric = new MetricExtremeReissnerNordstromDihole;
        break;
    case enum_metric_friedman_nonempty_null:
        currMetric = new MetricFriedmanNonEmptyNull;
        break;
    case enum_metric_glampedakis:
        currMetric = new MetricGlampedakis;
        break;
    case enum_metric_geodel:
        currMetric = new MetricGoedel;
        break;
    case enum_metric_geodel_cart:
        currMetric = new MetricGoedelCart;
        break;
    case enum_metric_goedelscaled:
        currMetric = new MetricGoedelScaled;
        break;
    case enum_metric_goedelscaled_cart:
        currMetric = new MetricGoedelScaledCart;
        break;
    case enum_metric_halilsoy_wave:
        currMetric = new MetricHalilsoyWave;
        break;
    case enum_metric_hartlethorne_gb:
        currMetric = new MetricHartleThorneGB;
        break;
    case enum_metric_janewi:
        currMetric = new MetricJaNeWi;
        break;
    case enum_metric_kasner:
        currMetric = new MetricKasner;
        break;
    case enum_metric_kastortraschen:
        currMetric = new MetricKastorTraschen;
        break;
    case enum_metric_kerrbl:
        currMetric = new MetricKerrBL;
        break;
    case enum_metric_kottler:
        currMetric = new MetricKottler;
        break;
    case enum_metric_minkowski:
        currMetric = new MetricMinkowski;
        break;
    case enum_metric_minkowski_conf:
        currMetric = new MetricMinkowskiConformal;
        break;
    case enum_metric_minkowski_rotlattice:
        currMetric = new MetricMinkRotLattice;
        break;
    case enum_metric_morristhorne:
        currMetric = new MetricMorrisThorne;
        break;
    case enum_metric_painleve:
        currMetric = new MetricPainleveGullstrand;
        break;
    case enum_metric_plane_grav_wave:
        currMetric = new MetricPlaneGravWave;
        break;
    case enum_metric_reissner:
        currMetric = new MetricReissnerNordstrom;
        break;
    case enum_metric_rotdihole:
        currMetric = new MetricRotDihole;
        break;
    case enum_metric_schwarzschild:
        currMetric = new MetricSchwarzschild;
        break;
    case enum_metric_schwarzschild_cart:
        currMetric = new MetricSchwarzschildCart;
        break;
    case enum_metric_schwarzschild_isotropic:
        currMetric = new MetricSchwarzschildIsotropic;
        break;
    case enum_metric_schwarzschild_tortoise:
        currMetric = new MetricSchwarzschildTortoise;
        break;
    case enum_metric_spinning_string:
        currMetric = new MetricStraightSpinningString;
        break;
    case enum_metric_sultana_dyer:
        currMetric = new MetricSultanaDyer;
        break;
    case enum_metric_taub_nut:
        currMetric = new MetricTaubNUT;
        break;
    case enum_metric_teowhl:
        currMetric = new MetricTeoWHl;
        break;
    case enum_metric_tomimatsusato:
        currMetric = new MetricTomimatsuSato;
        break;
    case enum_metric_Petrov_TD_AI:
        currMetric = new MetricPTD_AI;
        break;
    case enum_metric_Petrov_TD_AII:
        currMetric = new MetricPTD_AII;
        break;
    case enum_metric_Petrov_TD_AIII:
        currMetric = new MetricPTD_AIII;
        break;
    case enum_metric_Petrov_TD_BI:
        currMetric = new MetricPTD_BI;
        break;
    case enum_metric_Petrov_TD_BII:
        currMetric = new MetricPTD_BII;
        break;
    case enum_metric_Petrov_TD_BIII:
        currMetric = new MetricPTD_BIII;
        break;
    case enum_metric_Petrov_TD_C:
        currMetric = new MetricPTD_C;
        break;
    case enum_metric_Pravda_C_Metric:
        currMetric = new MetricPravda_C;
        break;
    case enum_metric_Pravda_C_Can:
        currMetric = new MetricPravda_C_Can;
        break;
    }
    return currMetric;
}

} // end namespace m4d
