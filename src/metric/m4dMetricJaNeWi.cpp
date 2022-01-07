/**
 * @file    m4dMetricJaNeWi.cpp
 * @author  Thomas Mueller
 *
 * This file is part of the m4d-library.
 */
#include "m4dMetricJaNeWi.h"

double dzdr_janewi(double x, void* params)
{
    struct_janewi_params* par = (struct_janewi_params*)params;
    double rs = par->rs;
    double g = par->gamma;

    double alpha = 1.0 - rs / (x * g);
    double dzdr2 = rs * (4.0 * x * g * g - rs * (1.0 + g) * (1.0 + g)) / (4.0 * x * x * g * g * pow(alpha, g + 1.0));
    return sqrt(fabs(dzdr2));
}

namespace m4d {

MetricJaNeWi::MetricJaNeWi(double mass, double gamma)
{
    mMetricName = "JanisNewmanWinicour";
    setCoordType(enum_coordinate_spherical);

    mPhysicalUnits = enum_physical_constants_geom;
    mSpeedOfLight = 1.0;
    mGravConstant = 1.0;

    addParam("mass", mass);
    mMass = mass;
    rs = 2.0 * mGravConstant * mMass / (mSpeedOfLight * mSpeedOfLight);

    addParam("gamma", gamma);
    mGamma = gamma;

    calcCriticalPoint();
    setStandardValues();

    mDrawTypes.push_back(enum_draw_embedding);
    mDrawTypes.push_back(enum_draw_effpoti);

    /*  parameters for the embedding diagram  */
    if (!mEmbParam.empty()) {
        mEmbParam.clear();
    }
    mHaveEmbedding = true;

    mEmb_rmin = mCritPoint;
    mEmb_rmax = 5.0 * mCritPoint;
    mEmb_r_num = 20.0;
    mEmb_phi_num = 40.0;
    mEmb_rstep = (mEmb_rmax - mEmb_rmin) / mEmb_r_num;
    mEmb_phistep = 2.0 * M_PI / mEmb_phi_num;
    addEmbeddingParam("emb_rmin", mEmb_rmin);
    addEmbeddingParam("emb_rmax", mEmb_rmax);
    addEmbeddingParam("emb_r_num", 20.0);
    addEmbeddingParam("emb_phi_num", 40.0);

    w = gsl_integration_workspace_alloc(1000);
    F.function = &dzdr_janewi;
}

MetricJaNeWi::~MetricJaNeWi()
{
    gsl_integration_workspace_free(w);
}

// *********************************** public methods ******************************

bool MetricJaNeWi::calculateMetric(const double* pos)
{
    double r = pos[1];
    double theta = pos[2];

    double c = mSpeedOfLight;

    double t1 = 1 / mGamma;
    double t2 = rs * t1;
    double t6 = pow(1.0 - t2 / r, mGamma);
    double t7 = c * c;
    double t9 = 1 / t6;
    double t10 = r * r;
    double t11 = t10 * t9;
    double t12 = r * t9;
    double t15 = sin(theta);
    double t16 = t15 * t15;

    g_compts[0][0] = -t6 * t7;
    g_compts[0][1] = 0.0;
    g_compts[0][2] = 0.0;
    g_compts[0][3] = 0.0;
    g_compts[1][0] = 0.0;
    g_compts[1][1] = t9;
    g_compts[1][2] = 0.0;
    g_compts[1][3] = 0.0;
    g_compts[2][0] = 0.0;
    g_compts[2][1] = 0.0;
    g_compts[2][2] = t11 - t12 * t2;
    g_compts[2][3] = 0.0;
    g_compts[3][0] = 0.0;
    g_compts[3][1] = 0.0;
    g_compts[3][2] = 0.0;
    g_compts[3][3] = t11 * t16 - t12 * t16 * rs * t1;

    return true;
}

bool MetricJaNeWi::calculateChristoffels(const double* pos)
{
    double r = pos[1];
    double theta = pos[2];

    double c = mSpeedOfLight;

    double t1 = mGamma * r;
    double t2 = t1 - rs;
    double t3 = 1 / mGamma;
    double t5 = 1 / r;
    double t7 = pow(t2 * t3 * t5, mGamma);
    double t8 = t7 * t7;
    double t11 = c * c;
    double t13 = 1 / t2;
    double t17 = rs * mGamma;
    double t18 = t5 * t13;
    double t20 = t17 * t18 / 2.0;
    double t22 = 2.0 * t1 - t17 - rs;
    double t24 = t18 * t22 / 2.0;
    double t27 = sin(theta);
    double t29 = cos(theta);
    double t30 = 1 / t27 * t29;
    double t31 = t27 * t27;

    christoffel[0][0][0] = 0.0;
    christoffel[0][0][1] = t8 * rs * mGamma * t5 * t11 * t13 / 2.0;
    christoffel[0][0][2] = 0.0;
    christoffel[0][0][3] = 0.0;
    christoffel[0][1][0] = t20;
    christoffel[0][1][1] = 0.0;
    christoffel[0][1][2] = 0.0;
    christoffel[0][1][3] = 0.0;
    christoffel[0][2][0] = 0.0;
    christoffel[0][2][1] = 0.0;
    christoffel[0][2][2] = 0.0;
    christoffel[0][2][3] = 0.0;
    christoffel[0][3][0] = 0.0;
    christoffel[0][3][1] = 0.0;
    christoffel[0][3][2] = 0.0;
    christoffel[0][3][3] = 0.0;
    christoffel[1][0][0] = t20;
    christoffel[1][0][1] = 0.0;
    christoffel[1][0][2] = 0.0;
    christoffel[1][0][3] = 0.0;
    christoffel[1][1][0] = 0.0;
    christoffel[1][1][1] = -t20;
    christoffel[1][1][2] = 0.0;
    christoffel[1][1][3] = 0.0;
    christoffel[1][2][0] = 0.0;
    christoffel[1][2][1] = 0.0;
    christoffel[1][2][2] = t24;
    christoffel[1][2][3] = 0.0;
    christoffel[1][3][0] = 0.0;
    christoffel[1][3][1] = 0.0;
    christoffel[1][3][2] = 0.0;
    christoffel[1][3][3] = t24;
    christoffel[2][0][0] = 0.0;
    christoffel[2][0][1] = 0.0;
    christoffel[2][0][2] = 0.0;
    christoffel[2][0][3] = 0.0;
    christoffel[2][1][0] = 0.0;
    christoffel[2][1][1] = 0.0;
    christoffel[2][1][2] = t24;
    christoffel[2][1][3] = 0.0;
    christoffel[2][2][0] = 0.0;
    christoffel[2][2][1] = -t22 * t3 / 2.0;
    christoffel[2][2][2] = 0.0;
    christoffel[2][2][3] = 0.0;
    christoffel[2][3][0] = 0.0;
    christoffel[2][3][1] = 0.0;
    christoffel[2][3][2] = 0.0;
    christoffel[2][3][3] = t30;
    christoffel[3][0][0] = 0.0;
    christoffel[3][0][1] = 0.0;
    christoffel[3][0][2] = 0.0;
    christoffel[3][0][3] = 0.0;
    christoffel[3][1][0] = 0.0;
    christoffel[3][1][1] = 0.0;
    christoffel[3][1][2] = 0.0;
    christoffel[3][1][3] = t24;
    christoffel[3][2][0] = 0.0;
    christoffel[3][2][1] = 0.0;
    christoffel[3][2][2] = 0.0;
    christoffel[3][2][3] = t30;
    christoffel[3][3][0] = 0.0;
    christoffel[3][3][1] = -t22 * t31 * t3 / 2.0;
    christoffel[3][3][2] = -t27 * t29;
    christoffel[3][3][3] = 0.0;

    return true;
}

bool MetricJaNeWi::calculateChrisD(const double* pos)
{
    double r = pos[1];
    double theta = pos[2];

    double c = mSpeedOfLight;

    double t1 = mGamma * r;
    double t2 = t1 - rs;
    double t3 = 1 / mGamma;
    double t7 = pow(t2 * t3 / r, mGamma);
    double t8 = t7 * t7;
    double t11 = c * c;
    double t12 = rs * mGamma;
    double t14 = 2.0 * t1;
    double t17 = t2 * t2;
    double t18 = 1 / t17;
    double t19 = r * r;
    double t20 = 1 / t19;
    double t29 = t12 * (t14 - rs) * t18 * t20 / 2.0;
    double t30 = mGamma * mGamma;
    double t38 = rs * rs;
    double t43 = (2.0 * t30 * t19 - 2.0 * t1 * rs - 2.0 * t30 * r * rs + t38 * mGamma + t38) * t18 * t20 / 2.0;
    double t44 = cos(theta);
    double t45 = t44 * t44;
    double t46 = sin(theta);
    double t47 = t46 * t46;
    double t50 = (t45 + t47) / t47;

    chrisD[0][0][0][0] = 0.0;
    chrisD[0][0][0][1] = 0.0;
    chrisD[0][0][0][2] = 0.0;
    chrisD[0][0][0][3] = 0.0;
    chrisD[0][0][1][0] = 0.0;
    chrisD[0][0][1][1] = -t8 * rs * mGamma * t11 * (-2.0 * t12 + t14 - rs) * t18 * t20 / 2.0;
    chrisD[0][0][1][2] = 0.0;
    chrisD[0][0][1][3] = 0.0;
    chrisD[0][0][2][0] = 0.0;
    chrisD[0][0][2][1] = 0.0;
    chrisD[0][0][2][2] = 0.0;
    chrisD[0][0][2][3] = 0.0;
    chrisD[0][0][3][0] = 0.0;
    chrisD[0][0][3][1] = 0.0;
    chrisD[0][0][3][2] = 0.0;
    chrisD[0][0][3][3] = 0.0;
    chrisD[0][1][0][0] = 0.0;
    chrisD[0][1][0][1] = -t29;
    chrisD[0][1][0][2] = 0.0;
    chrisD[0][1][0][3] = 0.0;
    chrisD[0][1][1][0] = 0.0;
    chrisD[0][1][1][1] = 0.0;
    chrisD[0][1][1][2] = 0.0;
    chrisD[0][1][1][3] = 0.0;
    chrisD[0][1][2][0] = 0.0;
    chrisD[0][1][2][1] = 0.0;
    chrisD[0][1][2][2] = 0.0;
    chrisD[0][1][2][3] = 0.0;
    chrisD[0][1][3][0] = 0.0;
    chrisD[0][1][3][1] = 0.0;
    chrisD[0][1][3][2] = 0.0;
    chrisD[0][1][3][3] = 0.0;
    chrisD[0][2][0][0] = 0.0;
    chrisD[0][2][0][1] = 0.0;
    chrisD[0][2][0][2] = 0.0;
    chrisD[0][2][0][3] = 0.0;
    chrisD[0][2][1][0] = 0.0;
    chrisD[0][2][1][1] = 0.0;
    chrisD[0][2][1][2] = 0.0;
    chrisD[0][2][1][3] = 0.0;
    chrisD[0][2][2][0] = 0.0;
    chrisD[0][2][2][1] = 0.0;
    chrisD[0][2][2][2] = 0.0;
    chrisD[0][2][2][3] = 0.0;
    chrisD[0][2][3][0] = 0.0;
    chrisD[0][2][3][1] = 0.0;
    chrisD[0][2][3][2] = 0.0;
    chrisD[0][2][3][3] = 0.0;
    chrisD[0][3][0][0] = 0.0;
    chrisD[0][3][0][1] = 0.0;
    chrisD[0][3][0][2] = 0.0;
    chrisD[0][3][0][3] = 0.0;
    chrisD[0][3][1][0] = 0.0;
    chrisD[0][3][1][1] = 0.0;
    chrisD[0][3][1][2] = 0.0;
    chrisD[0][3][1][3] = 0.0;
    chrisD[0][3][2][0] = 0.0;
    chrisD[0][3][2][1] = 0.0;
    chrisD[0][3][2][2] = 0.0;
    chrisD[0][3][2][3] = 0.0;
    chrisD[0][3][3][0] = 0.0;
    chrisD[0][3][3][1] = 0.0;
    chrisD[0][3][3][2] = 0.0;
    chrisD[0][3][3][3] = 0.0;
    chrisD[1][0][0][0] = 0.0;
    chrisD[1][0][0][1] = -t29;
    chrisD[1][0][0][2] = 0.0;
    chrisD[1][0][0][3] = 0.0;
    chrisD[1][0][1][0] = 0.0;
    chrisD[1][0][1][1] = 0.0;
    chrisD[1][0][1][2] = 0.0;
    chrisD[1][0][1][3] = 0.0;
    chrisD[1][0][2][0] = 0.0;
    chrisD[1][0][2][1] = 0.0;
    chrisD[1][0][2][2] = 0.0;
    chrisD[1][0][2][3] = 0.0;
    chrisD[1][0][3][0] = 0.0;
    chrisD[1][0][3][1] = 0.0;
    chrisD[1][0][3][2] = 0.0;
    chrisD[1][0][3][3] = 0.0;
    chrisD[1][1][0][0] = 0.0;
    chrisD[1][1][0][1] = 0.0;
    chrisD[1][1][0][2] = 0.0;
    chrisD[1][1][0][3] = 0.0;
    chrisD[1][1][1][0] = 0.0;
    chrisD[1][1][1][1] = t29;
    chrisD[1][1][1][2] = 0.0;
    chrisD[1][1][1][3] = 0.0;
    chrisD[1][1][2][0] = 0.0;
    chrisD[1][1][2][1] = 0.0;
    chrisD[1][1][2][2] = 0.0;
    chrisD[1][1][2][3] = 0.0;
    chrisD[1][1][3][0] = 0.0;
    chrisD[1][1][3][1] = 0.0;
    chrisD[1][1][3][2] = 0.0;
    chrisD[1][1][3][3] = 0.0;
    chrisD[1][2][0][0] = 0.0;
    chrisD[1][2][0][1] = 0.0;
    chrisD[1][2][0][2] = 0.0;
    chrisD[1][2][0][3] = 0.0;
    chrisD[1][2][1][0] = 0.0;
    chrisD[1][2][1][1] = 0.0;
    chrisD[1][2][1][2] = 0.0;
    chrisD[1][2][1][3] = 0.0;
    chrisD[1][2][2][0] = 0.0;
    chrisD[1][2][2][1] = -t43;
    chrisD[1][2][2][2] = 0.0;
    chrisD[1][2][2][3] = 0.0;
    chrisD[1][2][3][0] = 0.0;
    chrisD[1][2][3][1] = 0.0;
    chrisD[1][2][3][2] = 0.0;
    chrisD[1][2][3][3] = 0.0;
    chrisD[1][3][0][0] = 0.0;
    chrisD[1][3][0][1] = 0.0;
    chrisD[1][3][0][2] = 0.0;
    chrisD[1][3][0][3] = 0.0;
    chrisD[1][3][1][0] = 0.0;
    chrisD[1][3][1][1] = 0.0;
    chrisD[1][3][1][2] = 0.0;
    chrisD[1][3][1][3] = 0.0;
    chrisD[1][3][2][0] = 0.0;
    chrisD[1][3][2][1] = 0.0;
    chrisD[1][3][2][2] = 0.0;
    chrisD[1][3][2][3] = 0.0;
    chrisD[1][3][3][0] = 0.0;
    chrisD[1][3][3][1] = -t43;
    chrisD[1][3][3][2] = 0.0;
    chrisD[1][3][3][3] = 0.0;
    chrisD[2][0][0][0] = 0.0;
    chrisD[2][0][0][1] = 0.0;
    chrisD[2][0][0][2] = 0.0;
    chrisD[2][0][0][3] = 0.0;
    chrisD[2][0][1][0] = 0.0;
    chrisD[2][0][1][1] = 0.0;
    chrisD[2][0][1][2] = 0.0;
    chrisD[2][0][1][3] = 0.0;
    chrisD[2][0][2][0] = 0.0;
    chrisD[2][0][2][1] = 0.0;
    chrisD[2][0][2][2] = 0.0;
    chrisD[2][0][2][3] = 0.0;
    chrisD[2][0][3][0] = 0.0;
    chrisD[2][0][3][1] = 0.0;
    chrisD[2][0][3][2] = 0.0;
    chrisD[2][0][3][3] = 0.0;
    chrisD[2][1][0][0] = 0.0;
    chrisD[2][1][0][1] = 0.0;
    chrisD[2][1][0][2] = 0.0;
    chrisD[2][1][0][3] = 0.0;
    chrisD[2][1][1][0] = 0.0;
    chrisD[2][1][1][1] = 0.0;
    chrisD[2][1][1][2] = 0.0;
    chrisD[2][1][1][3] = 0.0;
    chrisD[2][1][2][0] = 0.0;
    chrisD[2][1][2][1] = -t43;
    chrisD[2][1][2][2] = 0.0;
    chrisD[2][1][2][3] = 0.0;
    chrisD[2][1][3][0] = 0.0;
    chrisD[2][1][3][1] = 0.0;
    chrisD[2][1][3][2] = 0.0;
    chrisD[2][1][3][3] = 0.0;
    chrisD[2][2][0][0] = 0.0;
    chrisD[2][2][0][1] = 0.0;
    chrisD[2][2][0][2] = 0.0;
    chrisD[2][2][0][3] = 0.0;
    chrisD[2][2][1][0] = 0.0;
    chrisD[2][2][1][1] = -1.0;
    chrisD[2][2][1][2] = 0.0;
    chrisD[2][2][1][3] = 0.0;
    chrisD[2][2][2][0] = 0.0;
    chrisD[2][2][2][1] = 0.0;
    chrisD[2][2][2][2] = 0.0;
    chrisD[2][2][2][3] = 0.0;
    chrisD[2][2][3][0] = 0.0;
    chrisD[2][2][3][1] = 0.0;
    chrisD[2][2][3][2] = 0.0;
    chrisD[2][2][3][3] = 0.0;
    chrisD[2][3][0][0] = 0.0;
    chrisD[2][3][0][1] = 0.0;
    chrisD[2][3][0][2] = 0.0;
    chrisD[2][3][0][3] = 0.0;
    chrisD[2][3][1][0] = 0.0;
    chrisD[2][3][1][1] = 0.0;
    chrisD[2][3][1][2] = 0.0;
    chrisD[2][3][1][3] = 0.0;
    chrisD[2][3][2][0] = 0.0;
    chrisD[2][3][2][1] = 0.0;
    chrisD[2][3][2][2] = 0.0;
    chrisD[2][3][2][3] = 0.0;
    chrisD[2][3][3][0] = 0.0;
    chrisD[2][3][3][1] = 0.0;
    chrisD[2][3][3][2] = -t50;
    chrisD[2][3][3][3] = 0.0;
    chrisD[3][0][0][0] = 0.0;
    chrisD[3][0][0][1] = 0.0;
    chrisD[3][0][0][2] = 0.0;
    chrisD[3][0][0][3] = 0.0;
    chrisD[3][0][1][0] = 0.0;
    chrisD[3][0][1][1] = 0.0;
    chrisD[3][0][1][2] = 0.0;
    chrisD[3][0][1][3] = 0.0;
    chrisD[3][0][2][0] = 0.0;
    chrisD[3][0][2][1] = 0.0;
    chrisD[3][0][2][2] = 0.0;
    chrisD[3][0][2][3] = 0.0;
    chrisD[3][0][3][0] = 0.0;
    chrisD[3][0][3][1] = 0.0;
    chrisD[3][0][3][2] = 0.0;
    chrisD[3][0][3][3] = 0.0;
    chrisD[3][1][0][0] = 0.0;
    chrisD[3][1][0][1] = 0.0;
    chrisD[3][1][0][2] = 0.0;
    chrisD[3][1][0][3] = 0.0;
    chrisD[3][1][1][0] = 0.0;
    chrisD[3][1][1][1] = 0.0;
    chrisD[3][1][1][2] = 0.0;
    chrisD[3][1][1][3] = 0.0;
    chrisD[3][1][2][0] = 0.0;
    chrisD[3][1][2][1] = 0.0;
    chrisD[3][1][2][2] = 0.0;
    chrisD[3][1][2][3] = 0.0;
    chrisD[3][1][3][0] = 0.0;
    chrisD[3][1][3][1] = -t43;
    chrisD[3][1][3][2] = 0.0;
    chrisD[3][1][3][3] = 0.0;
    chrisD[3][2][0][0] = 0.0;
    chrisD[3][2][0][1] = 0.0;
    chrisD[3][2][0][2] = 0.0;
    chrisD[3][2][0][3] = 0.0;
    chrisD[3][2][1][0] = 0.0;
    chrisD[3][2][1][1] = 0.0;
    chrisD[3][2][1][2] = 0.0;
    chrisD[3][2][1][3] = 0.0;
    chrisD[3][2][2][0] = 0.0;
    chrisD[3][2][2][1] = 0.0;
    chrisD[3][2][2][2] = 0.0;
    chrisD[3][2][2][3] = 0.0;
    chrisD[3][2][3][0] = 0.0;
    chrisD[3][2][3][1] = 0.0;
    chrisD[3][2][3][2] = -t50;
    chrisD[3][2][3][3] = 0.0;
    chrisD[3][3][0][0] = 0.0;
    chrisD[3][3][0][1] = 0.0;
    chrisD[3][3][0][2] = 0.0;
    chrisD[3][3][0][3] = 0.0;
    chrisD[3][3][1][0] = 0.0;
    chrisD[3][3][1][1] = -t47;
    chrisD[3][3][1][2] = -(t14 - t12 - rs) * t46 * t3 * t44;
    chrisD[3][3][1][3] = 0.0;
    chrisD[3][3][2][0] = 0.0;
    chrisD[3][3][2][1] = 0.0;
    chrisD[3][3][2][2] = -t45 + t47;
    chrisD[3][3][2][3] = 0.0;
    chrisD[3][3][3][0] = 0.0;
    chrisD[3][3][3][1] = 0.0;
    chrisD[3][3][3][2] = 0.0;
    chrisD[3][3][3][3] = 0.0;

    return true;
}

void MetricJaNeWi::localToCoord(const double* pos, const double* ldir, double* dir, enum_nat_tetrad_type)
{
    double r = pos[1];
    double theta = pos[2];
    double alpha = 1.0 - rs / (mGamma * r);
    double ag2 = pow(alpha, 0.5 * mGamma);
    double agm2 = pow(alpha, 0.5 * mGamma - 0.5);

    dir[0] = ldir[0] / ag2 / mSpeedOfLight;
    dir[1] = ldir[1] * ag2;
    dir[2] = ldir[2] * agm2 / r;
    dir[3] = ldir[3] * agm2 / (r * sin(theta));
}

void MetricJaNeWi::coordToLocal(const double* pos, const double* cdir, double* ldir, enum_nat_tetrad_type)
{
    double r = pos[1];
    double theta = pos[2];
    double alpha = 1.0 - rs / (mGamma * r);
    double ag2 = pow(alpha, 0.5 * mGamma);
    double agm2 = pow(alpha, 0.5 * mGamma - 0.5);

    ldir[0] = cdir[0] * ag2 * mSpeedOfLight;
    ldir[1] = cdir[1] / ag2;
    ldir[2] = cdir[2] * r / agm2;
    ldir[3] = cdir[3] * r * sin(theta) / agm2;
}

bool MetricJaNeWi::breakCondition(const double* pos)
{
    bool br = false;

    if ((pos[1] < 0.0) || (mGamma * mGamma * pos[1] * pos[1] <= (1.0 + M4D_METRIC_EPS) * rs * rs)
        || pos[1] < mCritPoint) {
        br = true;
    }
    return br;
}

bool MetricJaNeWi::setParam(const char* pName, double val)
{
    Metric::setParam(pName, val);

    if (strcmp(pName, "mass") == 0) {
        mMass = val;
        rs = 2.0 * mGravConstant * mMass / (mSpeedOfLight * mSpeedOfLight);
    }
    else if (strcmp(pName, "gamma") == 0) {
        mGamma = val;
    }
    calcCriticalPoint();
    return true;
}

bool MetricJaNeWi::transToEmbedding(vec4 p, vec4& ep)
{
    vec4 cp;
    transToPseudoCart(p, cp);

    double r = p[1];
    double x = cp[1];
    double y = cp[2];
    double z;
    if (r >= mCritPoint) {
        calcEmbeddingZ(r, z);
        ep = vec4(p[0], x, y, z);
        return true;
    }
    return false;
}

bool MetricJaNeWi::setEmbeddingParam(const char* name, double val)
{
    Metric::setEmbeddingParam(name, val);

    if (strcmp(name, "emb_rmin") == 0) {
        mEmb_rmin = val;
    }
    else if (strcmp(name, "emb_rmax") == 0) {
        mEmb_rmax = val;
    }
    else if (strcmp(name, "emb_r_num") == 0) {
        mEmb_r_num = val;
    }
    else if (strcmp(name, "emb_phi_num") == 0) {
        mEmb_phi_num = val;
    }
    return testEmbeddingParams();
}

bool MetricJaNeWi::testEmbeddingParams()
{
    bool allOk = true;
    if (mEmb_rmin < mCritPoint) {
        mEmb_rmin = mCritPoint;
        allOk &= false;
    }
    if (mEmb_rmax < mCritPoint) {
        mEmb_rmax = mCritPoint;
        allOk &= false;
    }
    if (mEmb_r_num < 2.0) {
        mEmb_r_num = 2.0;
        allOk &= false;
    }

    if (mEmb_phi_num < 4) {
        mEmb_phi_num = 4;
        allOk &= false;
    }
    return allOk;
}

bool MetricJaNeWi::effPotentialValue(
    const vec4 pos, const vec4 cdir, enum_geodesic_type type, const double x, double& val)
{
    double kappa = 0.0;
    if (type == enum_geodesic_timelike) {
        kappa = -mSign;
    }

    if (pos[1] < mCritPoint + 1e-2 || x < mCritPoint + 1e-2) {
        return false;
    }

    double alphaInit = 1.0 - rs / (pos[1] * mGamma);
    double alpha = 1.0 - rs / (x * mGamma);
    double h = pos[1] * pos[1] * pow(alphaInit, 1.0 - mGamma) * cdir[3];
    val = 0.5 * pow(alpha, mGamma)
        * (h * h / (x * x) * pow(alpha, mGamma - 1.0) - kappa * mSpeedOfLight * mSpeedOfLight);
    return true;
}

bool MetricJaNeWi::totEnergy(const vec4 pos, const vec4 cdir, const double, double& val)
{
    if (pos[1] < mCritPoint + 1e-2) {
        return false;
    }

    // 1/2*k^2/c^2:
    double alphaInit = 1.0 - rs / (pos[1] * mGamma);
    val = 0.5 * pow(alphaInit, 2.0 * mGamma) * mSpeedOfLight * mSpeedOfLight * cdir[0] * cdir[0];
    return true;
}

bool MetricJaNeWi::report(const vec4 pos, const vec4 cdir, char*& text)
{
    std::stringstream ss;
    ss << "Report for JaNeWi metric\n\tcoordinates : (t,r,theta,phi)\n";
    ss << "---------------------------------------------------------------\n";
    ss << "  physical units ......... no\n";
    ss.precision(DEF_FIXED_REPORT_PRECISION);
    ss.setf(std::ios::fixed);
    ss << "  Schwarzschild radius ... r_s  = 2GM/c^2 = " << rs << std::endl;
    ss << "  Photon orbit ........... r_ph = rs*(1+2*gamma)/(2*gamma)\n";
    ss << "                                = " << 0.5 * rs * (1.0 + 2.0 * mGamma) / mGamma << std::endl;
    calcCriticalPoint();
    ss << "  Critical point ......... r_cr = rs*(1+gamma)^2/(2*gamma)^2\n";
    ss << "                                = " << mCritPoint << std::endl;

    double alphaInit = 1.0 - rs / (pos[1] * mGamma);
    double h = pos[1] * pos[1] * pow(alphaInit, 1.0 - mGamma) * cdir[3];
    double k = pow(alphaInit, mGamma) * mSpeedOfLight * mSpeedOfLight * cdir[0];

    ss << "  constant of motion ........ k = " << k << std::endl;
    ss << "  constant of motion ........ h = " << h << std::endl;

    ss << "  Velocity for circular orbits: beta = sqrt[(r/rs)/(2*(r/rs)^2-(1+1/gamma)*(r/rs))]" << std::endl;
    ss << "                                     = " << getCircularVelocity(pos[1]) << std::endl;
    ss << "  Last innermost stable circular orbit: r_iszo=rs/(2gamma)+3/2*rs+rs/(2*gamma)*sqrt(5*gamma^2-1)"
       << std::endl;
    ss << "                                     = "
       << rs / (2.0 * mGamma) + 3.0 / 2.0 * rs + rs / (2.0 * mGamma) * sqrt(5.0 * mGamma * mGamma - 1) << std::endl;

    text = new char[ss.str().length() + 2];
    return CopyString(ss.str().c_str(), text);
}

void MetricJaNeWi::calcCriticalPoint()
{
    double f = 0.5 * (1.0 + mGamma) / mGamma;
    mCritPoint = rs * f * f;
}

bool MetricJaNeWi::calcEmbeddingZ(const double r, double& z)
{
    struct_janewi_params par = { rs, mGamma };
    F.params = &par;

    size_t limit = 1000;
    int key = GSL_INTEG_GAUSS15;
    double error;

    gsl_set_error_handler_off();
    gsl_integration_qag(&F, mCritPoint, r, 0.0, 1e-6, limit, key, w, &z, &error);
    return true;
}

double MetricJaNeWi::getCircularVelocity(const double r, const enum_nat_tetrad_type)
{
    if (r >= rs / 2.0 * ((1.0 / mGamma + 3.0) + sqrt(5.0 * mGamma * mGamma - 1.0) / mGamma)) {
        return sqrt((r / rs) / (2.0 * pow((r / rs), 2.0) - (1.0 + 1.0 / mGamma) * (r / rs)));
    }
    return 0.0;
}

vec4 MetricJaNeWi::getCircularFourVel(const vec4 pos, const enum_nat_tetrad_type)
{
    double beta = getCircularVelocity(pos[1]);
    if (beta > 0.0 && beta < 1.0) {
        double gamma = 1.0 / sqrt(1.0 - beta * beta);
        vec4 e0, e1, e2, e3;
        getNatTetrad(pos, e0, e1, e2, e3);
        return mSpeedOfLight * gamma * (e0 + beta * e3);
    }
    return vec4();
}

// ********************************* protected methods *****************************
void MetricJaNeWi::setStandardValues()
{
    mInitPos[0] = 0.0;
    mInitPos[1] = 3.0 * rs;
    mInitPos[2] = M_PI_2;
    mInitPos[3] = 0.0;
    mInitDir[0] = 1.0;
    mInitDir[1] = 0.0;
    mInitDir[2] = 0.0;

    mCoordNames[0] = std::string("t");
    mCoordNames[1] = std::string("r");
    mCoordNames[2] = std::string("theta");
    mCoordNames[3] = std::string("phi");
}

} // end namespace m4d
