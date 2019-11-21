/**
 * @file    m4dMetricDatabase.cpp
 * @author  Thomas Mueller
 *
 *  This file is part of libMotion4D.
 */
#include "m4dMetricDatabase.h"

namespace m4d {

MetricDatabase::MetricDatabase()
{
    init();
}

MetricDatabase::~MetricDatabase()
{
    if (!mMetricMap.empty()) {
        mMetricMap.clear();
    }
}

int MetricDatabase::getNumMetrics()
{
    return MetricList::NUM_METRICS;
}

Metric* MetricDatabase::getMetric(MetricList::enum_metric num)
{
    return initializeMetric(num);
}

Metric* MetricDatabase::getMetric(const char* mName)
{
    mMetricMapItr = mMetricMap.find(mName);
    if (mMetricMapItr == mMetricMap.end()) {
        fprintf(stderr, "Metric '%s' is not implemented!\n", mName);
        return nullptr;
    }

    return initializeMetric(mMetricMapItr->second);
}

const char* MetricDatabase::getMetricName(MetricList::enum_metric num)
{
    if (int(num) >= 0 && int(num) < MetricList::NUM_METRICS) {
        return MetricList::stl_metric_names[num];
    }

    return nullptr;
}

MetricList::enum_metric MetricDatabase::getMetricNr(const char* mName)
{
    mMetricMapItr = mMetricMap.find(mName);
    if (mMetricMapItr == mMetricMap.end()) {
        fprintf(stderr, "Metric '%s' is not implemented!\n", mName);
        return MetricList::enum_metric_unknown;
    }

    return mMetricMapItr->second;
}

void MetricDatabase::printMetricList(FILE* fptr)
{
    fprintf(fptr, "      Metric name                # params\n");
    fprintf(fptr, "-----------------------------------------------------------\n");

    int numParams;
    Metric* metric = nullptr;
    std::vector<std::string> paramNames;
    double value;

    for (int i = 0; i < MetricList::NUM_METRICS; i++) {
        paramNames.clear();

        if ((metric = getMetric(MetricList::stl_metric_names[i])) != nullptr) {
            numParams = metric->getNumParams();
            for (int j = 0; j < numParams; j++) {
                const char* pname = metric->getParamName(j);
                if (pname != nullptr) {
                    paramNames.push_back(std::string(pname));
                }
            }
        }
        else {
            numParams = 0;
        }
        fprintf(fptr, "%-30s       %d     (", MetricList::stl_metric_names[i], numParams);

        if (metric != nullptr) {
            for (unsigned int j = 0; j < paramNames.size(); j++) {
                if (metric->getParam(paramNames[j].c_str(), value)) {
                    fprintf(fptr, "%s=%f ", paramNames[j].c_str(), value);
                }
            }
        }
        fprintf(fptr, ")\n");
    }
}

void MetricDatabase::init()
{
    for (int i = 0; i < MetricList::NUM_METRICS; i++) {
        mMetricMap.insert(std::pair<std::string, MetricList::enum_metric>(
            MetricList::stl_metric_names[i], MetricList::enum_metric(i)));
    }
}

// ---------------------------------------------------
//           Please sort by enum name !!
// ---------------------------------------------------
Metric* MetricDatabase::initializeMetric(MetricList::enum_metric num)
{
    Metric* currMetric = nullptr;

    switch (num) {

        case MetricList::enum_metric_unknown:
            currMetric = nullptr;
            break;

        case MetricList::enum_metric_minkowski:
            currMetric = new MetricMinkowski;
            break;

        case MetricList::enum_metric_schwarzschild:
            currMetric = new MetricSchwarzschild;
            break;

        case MetricList::enum_metric_schwarzschild_isotropic:
            currMetric = new MetricSchwarzschildIsotropic;
            break;

        case MetricList::enum_metric_kerrbl:
            currMetric = new MetricKerrBL;
            break;

        case MetricList::enum_metric_morristhorne:
            currMetric = new MetricMorrisThorne;
            break;

        case MetricList::enum_metric_cosmic_string_schwarzschild:
            currMetric = new MetricCosmicStringSchwarzschild;
            break;

#ifdef ALL_METRICS_AVAILABLE
        case MetricList::enum_metric_alcubierre:
            currMetric = new MetricAlcubierre;
            break;
        case MetricList::enum_metric_alcubierre_simple:
            currMetric = new MetricAlcubierreSimple;
            break;
        case MetricList::enum_metric_barriola:
            currMetric = new MetricBarriolaVilenkin;
            break;
        case MetricList::enum_metric_bertottikasner:
            currMetric = new MetricBertottiKasner;
            break;
        case MetricList::enum_metric_bessel_grav_wave_cart:
            currMetric = new MetricBesselGravWaveCart;
            break;
        case MetricList::enum_metric_bonnor:
            currMetric = new MetricBonnor;
            break;
        case MetricList::enum_metric_chazy_curzon_rot:
            currMetric = new MetricChazyCurzonRot;
            break;
        case MetricList::enum_metric_curzon:
            currMetric = new MetricCurzon;
            break;
        case MetricList::enum_metric_desitter_univ:
            currMetric = new MetricDeSitterUniv;
            break;
        case MetricList::enum_metric_desitter_univ_conf:
            currMetric = new MetricDeSitterUnivConf;
            break;
        case MetricList::enum_metric_eddfinkin:
            currMetric = new MetricEddFinkIn;
            break;
        case MetricList::enum_metric_einstein_rosen_wave_wwb:
            currMetric = new MetricEinsteinRosenWaveWWB;
            break;
        case MetricList::enum_metric_erezrosenvar:
            currMetric = new MetricErezRosenVar;
            break;
        case MetricList::enum_metric_ernst:
            currMetric = new MetricErnst;
            break;
        case MetricList::enum_metric_extreme_reissner_dihole:
            currMetric = new MetricExtremeReissnerNordstromDihole;
            break;
        case MetricList::enum_metric_friedman_nonempty_null:
            currMetric = new MetricFriedmanNonEmptyNull;
            break;
        case MetricList::enum_metric_glampedakis:
            currMetric = new MetricGlampedakis;
            break;
        case MetricList::enum_metric_geodel:
            currMetric = new MetricGoedel;
            break;
        case MetricList::enum_metric_geodel_cart:
            currMetric = new MetricGoedelCart;
            break;
        case MetricList::enum_metric_goedelscaled:
            currMetric = new MetricGoedelScaled;
            break;
        case MetricList::enum_metric_goedelscaled_cart:
            currMetric = new MetricGoedelScaledCart;
            break;
        case MetricList::enum_metric_halilsoy_wave:
            currMetric = new MetricHalilsoyWave;
            break;
        case MetricList::enum_metric_hartlethorne_gb:
            currMetric = new MetricHartleThorneGB;
            break;
        case MetricList::enum_metric_janewi:
            currMetric = new MetricJaNeWi;
            break;
        case MetricList::enum_metric_kasner:
            currMetric = new MetricKasner;
            break;
        case MetricList::enum_metric_kastortraschen:
            currMetric = new MetricKastorTraschen;
            break;
        case MetricList::enum_metric_kottler:
            currMetric = new MetricKottler;
            break;
        case MetricList::enum_metric_minkowski_conf:
            currMetric = new MetricMinkowskiConformal;
            break;
        case MetricList::enum_metric_minkowski_rotlattice:
            currMetric = new MetricMinkRotLattice;
            break;
        case MetricList::enum_metric_painleve:
            currMetric = new MetricPainleveGullstrand;
            break;
        case MetricList::enum_metric_plane_grav_wave:
            currMetric = new MetricPlaneGravWave;
            break;
        case MetricList::enum_metric_reissner:
            currMetric = new MetricReissnerNordstrom;
            break;
        case MetricList::enum_metric_rotdihole:
            currMetric = new MetricRotDihole;
            break;
        case MetricList::enum_metric_schwarzschild_gravwave:
            currMetric = new MetricSchwarzschildGravWave;
            break;
        case MetricList::enum_metric_schwarzschild_tortoise:
            currMetric = new MetricSchwarzschildTortoise;
            break;
        case MetricList::enum_metric_schwarzschild_wt:
            currMetric = new MetricSchwarzschildWT;
            break;
        case MetricList::enum_metric_spinning_string:
            currMetric = new MetricStraightSpinningString;
            break;
        case MetricList::enum_metric_sultana_dyer:
            currMetric = new MetricSultanaDyer;
            break;
        case MetricList::enum_metric_taub_nut:
            currMetric = new MetricTaubNUT;
            break;
        case MetricList::enum_metric_teoSimpleWH:
            currMetric = new MetricTeoSimpleWH;
            break;
        case MetricList::enum_metric_teowhl:
            currMetric = new MetricTeoWHl;
            break;
        case MetricList::enum_metric_tomimatsusato:
            currMetric = new MetricTomimatsuSato;
            break;
        case MetricList::enum_metric_Petrov_TD_AI:
            currMetric = new MetricPTD_AI;
            break;
        case MetricList::enum_metric_Petrov_TD_AII:
            currMetric = new MetricPTD_AII;
            break;
        case MetricList::enum_metric_Petrov_TD_AIII:
            currMetric = new MetricPTD_AIII;
            break;
        case MetricList::enum_metric_Petrov_TD_BI:
            currMetric = new MetricPTD_BI;
            break;
        case MetricList::enum_metric_Petrov_TD_BII:
            currMetric = new MetricPTD_BII;
            break;
        case MetricList::enum_metric_Petrov_TD_BIII:
            currMetric = new MetricPTD_BIII;
            break;
        case MetricList::enum_metric_Petrov_TD_C:
            currMetric = new MetricPTD_C;
            break;
        case MetricList::enum_metric_Pravda_C_Metric:
            currMetric = new MetricPravda_C;
            break;
        case MetricList::enum_metric_Pravda_C_Can:
            currMetric = new MetricPravda_C_Can;
            break;
        case MetricList::enum_metric_vaidyaincrad:
            currMetric = new MetricVaidyaIncRad;
            break;

#endif // ALL_METRICS_AVAILABLE
    }

    return currMetric;
}

} // end namespace m4d
