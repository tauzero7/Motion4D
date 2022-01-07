/**
 * @file    m4dMetricKottler.cpp
 * @author  Thomas Mueller
 *
 * This file is part of the m4d-library.
 */
#include "m4dMetricKottler.h"

double dzdr_kottler(double x, void* params)
{
    struct_kottler_params* par = (struct_kottler_params*)params;
    double rs = par->rs;
    double l = par->lambda;

    double w = rs / x + l * x * x / 3.0;
    double dzdr2 = w / (1.0 - w);
    return sqrt(fabs(dzdr2));
}

namespace m4d {
#define eps 1.0e-6

MetricKottler::MetricKottler(double mass, double lambda)
{
    mMetricName = "Kottler";
    setCoordType(enum_coordinate_spherical);

    mPhysicalUnits = enum_physical_constants_geom;
    mSpeedOfLight = 1.0;
    mGravConstant = 1.0;

    addParam("mass");
    setParam("mass", mass);
    mMass = mass;

    rs = 2.0 * mGravConstant * mMass / (mSpeedOfLight * mSpeedOfLight);

    addParam("lambda");
    setParam("lambda", lambda);
    mLambda = lambda;
    calcCriticalPoints();

    mLocTeds.push_back(enum_nat_tetrad_static);

    mDrawTypes.push_back(enum_draw_embedding);
    mDrawTypes.push_back(enum_draw_effpoti);
    mDrawTypes.push_back(enum_draw_custom);

    /*  parameters for the embedding diagram  */
    if (!mEmbParam.empty()) {
        mEmbParam.clear();
    }
    mHaveEmbedding = true;

    mEmb_rmin = rp + eps;
    mEmb_rmax = rm - eps;
    mEmb_r_num = 20.0;
    mEmb_phi_num = 40.0;
    mEmb_rstep = (mEmb_rmax - mEmb_rmin) / mEmb_r_num;
    mEmb_phistep = 2.0 * M_PI / mEmb_phi_num;
    addEmbeddingParam("emb_r_num", 20.0);
    addEmbeddingParam("emb_phi_num", 40.0);
    addEmbeddingParam("emb_rmax", mEmb_rmax);

    w = gsl_integration_workspace_alloc(1000);
    F.function = &dzdr_kottler;

    setStandardValues();
}

MetricKottler::~MetricKottler()
{
    gsl_integration_workspace_free(w);
}

// *********************************** public methods ******************************

bool MetricKottler::calculateMetric(const double* pos)
{
    double r = pos[1];
    double theta = pos[2];

    double c = mSpeedOfLight;

    double t1 = c * c;
    double t3 = 1 / r;
    double t6 = r * r;
    double t15 = sin(theta);
    double t16 = t15 * t15;

    g_compts[0][0] = -t1 + t1 * rs * t3 + t1 * mLambda * t6 / 3.0;
    g_compts[0][1] = 0.0;
    g_compts[0][2] = 0.0;
    g_compts[0][3] = 0.0;
    g_compts[1][0] = 0.0;
    g_compts[1][1] = 1 / (1.0 - rs * t3 - mLambda * t6 / 3.0);
    g_compts[1][2] = 0.0;
    g_compts[1][3] = 0.0;
    g_compts[2][0] = 0.0;
    g_compts[2][1] = 0.0;
    g_compts[2][2] = t6;
    g_compts[2][3] = 0.0;
    g_compts[3][0] = 0.0;
    g_compts[3][1] = 0.0;
    g_compts[3][2] = 0.0;
    g_compts[3][3] = t6 * t16;

    return true;
}

bool MetricKottler::calculateChristoffels(const double* pos)
{
    double r = pos[1];
    double theta = pos[2];

    double c = mSpeedOfLight;

    double t2 = 3.0 * rs;
    double t3 = r * r;
    double t4 = t3 * r;
    double t5 = mLambda * t4;
    double t6 = -3.0 * r + t2 + t5;
    double t9 = c * c;
    double t11 = -t2 + 2.0 * t5;
    double t15 = 1 / r;
    double t19 = t15 / t6 * t11 / 2.0;
    double t22 = sin(theta);
    double t24 = cos(theta);
    double t25 = 1 / t22 * t24;
    double t26 = t22 * t22;

    christoffel[0][0][0] = 0.0;
    christoffel[0][0][1] = t6 / t4 * t9 * t11 / 18.0;
    christoffel[0][0][2] = 0.0;
    christoffel[0][0][3] = 0.0;
    christoffel[0][1][0] = t19;
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
    christoffel[1][0][0] = t19;
    christoffel[1][0][1] = 0.0;
    christoffel[1][0][2] = 0.0;
    christoffel[1][0][3] = 0.0;
    christoffel[1][1][0] = 0.0;
    christoffel[1][1][1] = -t19;
    christoffel[1][1][2] = 0.0;
    christoffel[1][1][3] = 0.0;
    christoffel[1][2][0] = 0.0;
    christoffel[1][2][1] = 0.0;
    christoffel[1][2][2] = t15;
    christoffel[1][2][3] = 0.0;
    christoffel[1][3][0] = 0.0;
    christoffel[1][3][1] = 0.0;
    christoffel[1][3][2] = 0.0;
    christoffel[1][3][3] = t15;
    christoffel[2][0][0] = 0.0;
    christoffel[2][0][1] = 0.0;
    christoffel[2][0][2] = 0.0;
    christoffel[2][0][3] = 0.0;
    christoffel[2][1][0] = 0.0;
    christoffel[2][1][1] = 0.0;
    christoffel[2][1][2] = t15;
    christoffel[2][1][3] = 0.0;
    christoffel[2][2][0] = 0.0;
    christoffel[2][2][1] = -r + rs + t5 / 3.0;
    christoffel[2][2][2] = 0.0;
    christoffel[2][2][3] = 0.0;
    christoffel[2][3][0] = 0.0;
    christoffel[2][3][1] = 0.0;
    christoffel[2][3][2] = 0.0;
    christoffel[2][3][3] = t25;
    christoffel[3][0][0] = 0.0;
    christoffel[3][0][1] = 0.0;
    christoffel[3][0][2] = 0.0;
    christoffel[3][0][3] = 0.0;
    christoffel[3][1][0] = 0.0;
    christoffel[3][1][1] = 0.0;
    christoffel[3][1][2] = 0.0;
    christoffel[3][1][3] = t15;
    christoffel[3][2][0] = 0.0;
    christoffel[3][2][1] = 0.0;
    christoffel[3][2][2] = 0.0;
    christoffel[3][2][3] = t25;
    christoffel[3][3][0] = 0.0;
    christoffel[3][3][1] = t6 * t26 / 3.0;
    christoffel[3][3][2] = -t22 * t24;
    christoffel[3][3][3] = 0.0;

    return true;
}

bool MetricKottler::calculateChrisD(const double* pos)
{
    double r = pos[1];
    double theta = pos[2];

    double c = mSpeedOfLight;

    double t1 = c * c;
    double t2 = r * rs;
    double t4 = r * r;
    double t5 = t4 * t4;
    double t6 = mLambda * t5;
    double t8 = mLambda * mLambda;
    double t11 = 2.0 * t8 * t5 * t4;
    double t12 = rs * rs;
    double t13 = 9.0 * t12;
    double t22 = mLambda * t4 * r;
    double t26 = 1 / t4;
    double t30 = -3.0 * r + 3.0 * rs + t22;
    double t31 = t30 * t30;
    double t34 = (18.0 * t2 + 6.0 * t6 - t13 - 24.0 * t22 * rs + t11) * t26 / t31 / 2.0;
    double t36 = -1.0 + mLambda * t4;
    double t37 = cos(theta);
    double t38 = t37 * t37;
    double t39 = sin(theta);
    double t40 = t39 * t39;
    double t43 = (t38 + t40) / t40;

    chrisD[0][0][0][0] = 0.0;
    chrisD[0][0][0][1] = 0.0;
    chrisD[0][0][0][2] = 0.0;
    chrisD[0][0][0][3] = 0.0;
    chrisD[0][0][1][0] = 0.0;
    chrisD[0][0][1][1] = t1 * (-6.0 * t2 - 2.0 * t6 + t11 + t13) / t5 / 6.0;
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
    chrisD[0][1][0][1] = -t34;
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
    chrisD[1][0][0][1] = -t34;
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
    chrisD[1][1][1][1] = t34;
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
    chrisD[1][2][2][1] = -t26;
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
    chrisD[1][3][3][1] = -t26;
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
    chrisD[2][1][2][1] = -t26;
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
    chrisD[2][2][1][1] = t36;
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
    chrisD[2][3][3][2] = -t43;
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
    chrisD[3][1][3][1] = -t26;
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
    chrisD[3][2][3][2] = -t43;
    chrisD[3][2][3][3] = 0.0;
    chrisD[3][3][0][0] = 0.0;
    chrisD[3][3][0][1] = 0.0;
    chrisD[3][3][0][2] = 0.0;
    chrisD[3][3][0][3] = 0.0;
    chrisD[3][3][1][0] = 0.0;
    chrisD[3][3][1][1] = t36 * t40;
    chrisD[3][3][1][2] = 2.0 / 3.0 * t30 * t39 * t37;
    chrisD[3][3][1][3] = 0.0;
    chrisD[3][3][2][0] = 0.0;
    chrisD[3][3][2][1] = 0.0;
    chrisD[3][3][2][2] = -t38 + t40;
    chrisD[3][3][2][3] = 0.0;
    chrisD[3][3][3][0] = 0.0;
    chrisD[3][3][3][1] = 0.0;
    chrisD[3][3][3][2] = 0.0;
    chrisD[3][3][3][3] = 0.0;

    return true;
}

void MetricKottler::localToCoord(const double* pos, const double* ldir, double* dir, enum_nat_tetrad_type)
{
    double r = pos[1];
    double theta = pos[2];
    double alpha = sqrt(1.0 - rs / r - mLambda / 3.0 * r * r);

    dir[0] = ldir[0] / (mSpeedOfLight * alpha);
    dir[1] = ldir[1] * alpha;
    dir[2] = ldir[2] / r;
    dir[3] = ldir[3] / (r * sin(theta));
}

void MetricKottler::coordToLocal(const double* pos, const double* cdir, double* ldir, enum_nat_tetrad_type)
{
    double r = pos[1];
    double theta = pos[2];
    double alpha = sqrt(1.0 - rs / r - mLambda / 3.0 * r * r);

    ldir[0] = cdir[0] * mSpeedOfLight * alpha;
    ldir[1] = cdir[1] / alpha;
    ldir[2] = cdir[2] * r;
    ldir[3] = cdir[3] * r * sin(theta);
}

bool MetricKottler::breakCondition(const double* pos)
{
    bool br = false;

    double r = pos[1];
    if ((r < 0.0) || (fabs(1.0 - rs / r - mLambda / 3.0 * r * r) <= eps)) {
        br = true;
    }
    return br;
}

double MetricKottler::testConstraint(const double y[], const double kappa)
{
    double r = y[1];
    double theta = y[2];

    double alpha = 1.0 - rs / r - mLambda / 3.0 * r * r;
    double cm = 1.0 / mSpeedOfLight;

    // Scale the directions with the speed of light before doubling them !!
    double dt = y[4];
    double dr = y[5] * cm;
    double dth = y[6] * cm;
    double dph = y[7] * cm;

    double sum = -kappa;
    sum += -alpha * dt * dt + dr * dr / alpha + r * r * (dth * dth + sin(theta) * sin(theta) * dph * dph);
    return sum;
}

bool MetricKottler::setParam(const char* pName, double val)
{
    Metric::setParam(pName, val);

    if (strcmp(pName, "mass") == 0) {
        mMass = val;
        rs = 2.0 * mGravConstant * mMass / (mSpeedOfLight * mSpeedOfLight);
    }
    else if (strcmp(pName, "lambda") == 0) {
        mLambda = val;
    }

    calcCriticalPoints();
    if (mLambda < 0.0) {
        mEmb_rmin = r1 + eps;
    }
    else {
        mEmb_rmin = rp + eps;
    }
    testEmbeddingParams();
    return true;
}

bool MetricKottler::transToEmbedding(vec4 p, vec4& ep)
{
    vec4 cp;
    transToPseudoCart(p, cp);

    double r = p[1];
    double x = cp[1];
    double y = cp[2];
    double z;

    if (mLambda < 0.0) {
        if (r > r1 + eps) {
            calcEmbeddingZ(r, z);
            ep = vec4(p[0], x, y, z);
            return true;
        }
    }
    else {
        if (r > rp + eps && r < rm - eps) {
            calcEmbeddingZ(r, z);
            ep = vec4(p[0], x, y, z);
            return true;
        }
    }
    return false;
}

bool MetricKottler::setEmbeddingParam(const char* name, double val)
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

bool MetricKottler::testEmbeddingParams()
{
    bool allOk = true;

    if (mLambda < 0.0) {
        if (mEmb_rmin < r1 + eps) {
            mEmb_rmin = r1 + eps;
            allOk &= false;
        }
        if (mEmb_rmax < r1 + eps) {
            mEmb_rmax = r1 + eps;
            allOk &= false;
        }
    }
    else {
        if (mEmb_rmin < rp + eps) {
            mEmb_rmin = rp + eps;
            allOk &= false;
        }
        if (mEmb_rmax < rp + eps) {
            mEmb_rmax = rp + eps;
            allOk &= false;
        }
        if (mEmb_rmin > rm - eps) {
            mEmb_rmin = rm - eps;
            allOk &= false;
        }
        if (mEmb_rmax > rm - eps) {
            mEmb_rmax = rm - eps;
            allOk &= false;
        }
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

void MetricKottler::usePhysicalUnits(const enum_physical_constants units)
{
    Metric::usePhysicalUnits(units);
    rs = 2.0 * mGravConstant * mMass / (mSpeedOfLight * mSpeedOfLight);
}

void MetricKottler::setUnits(const double speed_of_light, const double grav_const, const double diel_perm)
{
    Metric::setUnits(speed_of_light, grav_const, diel_perm);
    rs = 2.0 * mGravConstant * mMass / (mSpeedOfLight * mSpeedOfLight);
}

bool MetricKottler::effPotentialValue(
    const vec4 pos, const vec4 cdir, enum_geodesic_type type, const double x, double& val)
{
    double kappa = 0.0;
    if (type == enum_geodesic_timelike) {
        kappa = -mSign;
    }

    if (mLambda < 0.0) {
        if (pos[1] < r1 + eps || x < r1 + eps) {
            return false;
        }
    }
    else {
        if (pos[1] < rp + eps || x < rp + eps || pos[1] > rm - eps || x > rm - eps) {
            return false;
        }
    }

    double h = pos[1] * pos[1] * cdir[3];
    val = 0.5 * (1.0 - rs / x - mLambda * x * x / 3.0) * (h * h / (x * x) - kappa * mSpeedOfLight * mSpeedOfLight);
    return true;
}

bool MetricKottler::totEnergy(const vec4 pos, const vec4 cdir, const double, double& val)
{
    if (mLambda < 0.0) {
        if (pos[1] < r1 + eps) {
            return false;
        }
    }
    else {
        if (pos[1] < rp + eps || pos[1] > rm - eps) {
            return false;
        }
    }

    // 1/2*k^2/c^2:
    val = 0.5 * pow(1.0 - rs / pos[1] - mLambda * pos[1] * pos[1] / 3.0, 2.0) * mSpeedOfLight * mSpeedOfLight * cdir[0]
        * cdir[0];
    return true;
}

bool MetricKottler::report(const vec4 pos, const vec4 cdir, char*& text)
{
    std::stringstream ss;
    ss << "Report for Kottler metric\n\tcoordinate : (t,r,theta,phi)\n";
    ss << "---------------------------------------------------------------\n";
    ss << "  physical units ................................. yes\n";

    ss.precision(DEF_FIXED_REPORT_PRECISION);
    ss.setf(std::ios::fixed);
    ss << "  Schwarzschild radius ........... r_s = 2GM/c^2 = " << rs << std::endl;

    calcCriticalPoints();
    ss << "  Critical points..................r_1 = " << r1 << std::endl;
    ss << "                                   r_p = " << rp << std::endl;
    ss << "                                   r_m = " << rm << std::endl;

    double h = pos[1] * pos[1] * cdir[3];
    double k = (1.0 - rs / pos[1] - mLambda * pos[1] * pos[1] / 3.0) * mSpeedOfLight * mSpeedOfLight * cdir[0];
    ss << "  constant of motion ............... k = " << k << " " << k * k * 0.5 << std::endl;
    ss << "  constant of motion ............... h = " << h << std::endl;

    text = new char[ss.str().length() + 2];
    return CopyString(ss.str().c_str(), text);
}

// ********************************* protected methods *****************************

void MetricKottler::setStandardValues()
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

void MetricKottler::calcCriticalPoints()
{
    if (mLambda < 0.0) {
        double wl = sqrt(-mLambda);
        r1 = 2.0 / wl * sinh(1.0 / 3.0 * asinh(1.5 * rs * wl));
        rp = rm = 0.0;
    }
    else {
        r1 = 0.0;
        double wl = sqrt(mLambda);
        double diskr = 2.25 * rs * rs / (mLambda * mLambda) - pow(mLambda, -3.0);
        if (diskr > 0.0) {
            rp = rm = 0.0;
        }
        else {
            rp = 2.0 / wl * cos(M_PI / 3.0 + acos(1.5 * rs * wl) / 3.0);
            rm = 2.0 / wl * cos(M_PI / 3.0 - acos(1.5 * rs * wl) / 3.0);
        }
    }
}

bool MetricKottler::calcEmbeddingZ(const double r, double& z)
{
    struct_kottler_params par = { rs, mLambda };
    F.params = &par;

    size_t limit = 1000;
    int key = GSL_INTEG_GAUSS15;
    double error;

    gsl_set_error_handler_off();
    if (mLambda < 0.0) {
        gsl_integration_qag(&F, r1, r, 0.0, 1e-6, limit, key, w, &z, &error);
    }
    else {
        gsl_integration_qag(&F, rp, r, 0.0, 1e-6, limit, key, w, &z, &error);
    }
    return true;
}

} // end namespace m4d
