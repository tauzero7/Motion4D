/**
 * @file    m4dMetricCosmicStringSchwarzschild.cpp
 * @author  Thomas Mueller
 *
 *  This file is part of libMotion4D.
 */
#include "m4dMetricCosmicStringSchwarzschild.h"
#include "extra/m4dUtilities.h"

namespace m4d {

#define eps 1.0e-6

MetricCosmicStringSchwarzschild::MetricCosmicStringSchwarzschild(double mass, double beta)
{
    mMetricName = "CosmicStringSchwarzschild";
    setCoordType(enum_coordinate_spherical);

    mPhysicalUnits = enum_physical_constants_geom;
    mSpeedOfLight = 1.0;
    mGravConstant = 1.0;

    addParam("mass", mass);
    addParam("beta", beta);

    mMass = mass;
    mBeta = beta;
    rs = 2.0 * mGravConstant * mMass / (mSpeedOfLight * mSpeedOfLight);

    mDrawTypes.push_back(enum_draw_embedding);
    mDrawTypes.push_back(enum_draw_twoplusone);
    mDrawTypes.push_back(enum_draw_effpoti);

    /*  parameters for the embedding diagram  */
    if (!mEmbParam.empty()) {
        mEmbParam.clear();
    }
    mHaveEmbedding = true;

    mEmb_rmin = rs;
    mEmb_rmax = 5.0 * rs;
    mEmb_r_num = 20.0;
    mEmb_phi_num = 40.0;
    mEmb_rstep = (mEmb_rmax - mEmb_rmin) / mEmb_r_num;
    mEmb_phistep = 2.0 * M_PI / mEmb_phi_num;
    addEmbeddingParam("emb_rmin", mEmb_rmin);
    addEmbeddingParam("emb_rmax", mEmb_rmax);
    addEmbeddingParam("emb_r_num", 20.0);
    addEmbeddingParam("emb_phi_num", 40.0);

    setStandardValues();
}

MetricCosmicStringSchwarzschild::~MetricCosmicStringSchwarzschild() {}

// *********************************** public methods ******************************

bool MetricCosmicStringSchwarzschild::calculateMetric(const double* pos)
{
    double r = pos[1];
    double theta = pos[2];

    double c = mSpeedOfLight;

    double t1 = c * c;
    double t3 = 1 / r;
    double t9 = r * r;
    double t10 = mBeta * mBeta;
    double t12 = sin(theta);
    double t13 = t12 * t12;

    g_compts[0][0] = -t1 + t1 * rs * t3;
    g_compts[0][1] = 0.0;
    g_compts[0][2] = 0.0;
    g_compts[0][3] = 0.0;
    g_compts[1][0] = 0.0;
    g_compts[1][1] = 1 / (1.0 - rs * t3);
    g_compts[1][2] = 0.0;
    g_compts[1][3] = 0.0;
    g_compts[2][0] = 0.0;
    g_compts[2][1] = 0.0;
    g_compts[2][2] = t9;
    g_compts[2][3] = 0.0;
    g_compts[3][0] = 0.0;
    g_compts[3][1] = 0.0;
    g_compts[3][2] = 0.0;
    g_compts[3][3] = t9 * t10 * t13;

    return true;
}

bool MetricCosmicStringSchwarzschild::calculateChristoffels(const double* pos)
{
    double r = pos[1];
    double theta = pos[2];

    double c = mSpeedOfLight;

    double t1 = r - rs;
    double t2 = r * r;
    double t6 = c * c;
    double t10 = 1 / r;
    double t14 = t10 / t1 * rs / 2.0;
    double t15 = sin(theta);
    double t17 = cos(theta);
    double t18 = 1 / t15 * t17;
    double t19 = mBeta * mBeta;
    double t21 = t15 * t15;

    christoffel[0][0][0] = 0.0;
    christoffel[0][0][1] = t1 / t2 / r * t6 * rs / 2.0;
    christoffel[0][0][2] = 0.0;
    christoffel[0][0][3] = 0.0;
    christoffel[0][1][0] = t14;
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
    christoffel[1][0][0] = t14;
    christoffel[1][0][1] = 0.0;
    christoffel[1][0][2] = 0.0;
    christoffel[1][0][3] = 0.0;
    christoffel[1][1][0] = 0.0;
    christoffel[1][1][1] = -t14;
    christoffel[1][1][2] = 0.0;
    christoffel[1][1][3] = 0.0;
    christoffel[1][2][0] = 0.0;
    christoffel[1][2][1] = 0.0;
    christoffel[1][2][2] = t10;
    christoffel[1][2][3] = 0.0;
    christoffel[1][3][0] = 0.0;
    christoffel[1][3][1] = 0.0;
    christoffel[1][3][2] = 0.0;
    christoffel[1][3][3] = t10;
    christoffel[2][0][0] = 0.0;
    christoffel[2][0][1] = 0.0;
    christoffel[2][0][2] = 0.0;
    christoffel[2][0][3] = 0.0;
    christoffel[2][1][0] = 0.0;
    christoffel[2][1][1] = 0.0;
    christoffel[2][1][2] = t10;
    christoffel[2][1][3] = 0.0;
    christoffel[2][2][0] = 0.0;
    christoffel[2][2][1] = -t1;
    christoffel[2][2][2] = 0.0;
    christoffel[2][2][3] = 0.0;
    christoffel[2][3][0] = 0.0;
    christoffel[2][3][1] = 0.0;
    christoffel[2][3][2] = 0.0;
    christoffel[2][3][3] = t18;
    christoffel[3][0][0] = 0.0;
    christoffel[3][0][1] = 0.0;
    christoffel[3][0][2] = 0.0;
    christoffel[3][0][3] = 0.0;
    christoffel[3][1][0] = 0.0;
    christoffel[3][1][1] = 0.0;
    christoffel[3][1][2] = 0.0;
    christoffel[3][1][3] = t10;
    christoffel[3][2][0] = 0.0;
    christoffel[3][2][1] = 0.0;
    christoffel[3][2][2] = 0.0;
    christoffel[3][2][3] = t18;
    christoffel[3][3][0] = 0.0;
    christoffel[3][3][1] = -t1 * t19 * t21;
    christoffel[3][3][2] = -t19 * t15 * t17;
    christoffel[3][3][3] = 0.0;

    return true;
}

bool MetricCosmicStringSchwarzschild::calculateChrisD(const double* pos)
{
    double r = pos[1];
    double theta = pos[2];

    double c = mSpeedOfLight;

    double t1 = c * c;
    double t3 = 2.0 * r;
    double t6 = r * r;
    double t7 = t6 * t6;
    double t14 = 1 / t6;
    double t15 = r - rs;
    double t16 = t15 * t15;
    double t20 = rs * (t3 - rs) * t14 / t16 / 2.0;
    double t21 = cos(theta);
    double t22 = t21 * t21;
    double t23 = sin(theta);
    double t24 = t23 * t23;
    double t27 = (t22 + t24) / t24;
    double t28 = mBeta * mBeta;
    double t29 = t28 * t24;

    chrisD[0][0][0][0] = 0.0;
    chrisD[0][0][0][1] = 0.0;
    chrisD[0][0][0][2] = 0.0;
    chrisD[0][0][0][3] = 0.0;
    chrisD[0][0][1][0] = 0.0;
    chrisD[0][0][1][1] = -t1 * rs * (t3 - 3.0 * rs) / t7 / 2.0;
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
    chrisD[0][1][0][1] = -t20;
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
    chrisD[1][0][0][1] = -t20;
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
    chrisD[1][1][1][1] = t20;
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
    chrisD[1][2][2][1] = -t14;
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
    chrisD[1][3][3][1] = -t14;
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
    chrisD[2][1][2][1] = -t14;
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
    chrisD[2][3][3][2] = -t27;
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
    chrisD[3][1][3][1] = -t14;
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
    chrisD[3][2][3][2] = -t27;
    chrisD[3][2][3][3] = 0.0;
    chrisD[3][3][0][0] = 0.0;
    chrisD[3][3][0][1] = 0.0;
    chrisD[3][3][0][2] = 0.0;
    chrisD[3][3][0][3] = 0.0;
    chrisD[3][3][1][0] = 0.0;
    chrisD[3][3][1][1] = -t29;
    chrisD[3][3][1][2] = -2.0 * t15 * t28 * t23 * t21;
    chrisD[3][3][1][3] = 0.0;
    chrisD[3][3][2][0] = 0.0;
    chrisD[3][3][2][1] = 0.0;
    chrisD[3][3][2][2] = -t28 * t22 + t29;
    chrisD[3][3][2][3] = 0.0;
    chrisD[3][3][3][0] = 0.0;
    chrisD[3][3][3][1] = 0.0;
    chrisD[3][3][3][2] = 0.0;
    chrisD[3][3][3][3] = 0.0;

    return true;
}

void MetricCosmicStringSchwarzschild::localToCoord(
    const double* pos, const double* ldir, double* dir, enum_nat_tetrad_type)
{
    double r = pos[1];
    double theta = pos[2];
    double w = sqrt(1.0 - rs / r);

    dir[0] = ldir[0] / w / mSpeedOfLight;
    dir[1] = ldir[1] * w;
    dir[2] = ldir[2] / r;
    dir[3] = ldir[3] / (r * mBeta * sin(theta));
}

void MetricCosmicStringSchwarzschild::coordToLocal(
    const double* pos, const double* cdir, double* ldir, enum_nat_tetrad_type)
{
    double r = pos[1];
    double theta = pos[2];
    double w = sqrt(1.0 - rs / r);

    ldir[0] = cdir[0] * w * mSpeedOfLight;
    ldir[1] = cdir[1] / w;
    ldir[2] = cdir[2] * r;
    ldir[3] = cdir[3] * r * mBeta * sin(theta);
}

bool MetricCosmicStringSchwarzschild::breakCondition(const double* pos)
{
    bool br = false;

    if ((pos[1] < 0.0) || (pos[1] * pos[1] <= (1.0 + eps) * rs * rs)) {
        br = true;
    }
    return br;
}

bool MetricCosmicStringSchwarzschild::calcDerivs(const double y[], double dydx[])
{
    double r = y[1];
    double theta = y[2];

    double c = mSpeedOfLight;

    double t1 = 1 / r;
    double t2 = r - rs;
    double t4 = t1 / t2;
    double t8 = r * r;
    double t12 = c * c;
    double t14 = y[4] * y[4];
    double t18 = y[5] * y[5];
    double t22 = y[6] * y[6];
    double t24 = mBeta * mBeta;
    double t26 = sin(theta);
    double t27 = t26 * t26;
    double t28 = y[7] * y[7];
    double t36 = cos(theta);

    dydx[0] = y[4];
    dydx[1] = y[5];
    dydx[2] = y[6];
    dydx[3] = y[7];
    dydx[4] = -t4 * rs * y[5] * y[4];
    dydx[5] = -t2 / t8 / r * t12 * rs * t14 / 2.0 + t4 * rs * t18 / 2.0 + t2 * t22 + t2 * t24 * t27 * t28;
    dydx[6] = -2.0 * t1 * y[6] * y[5] + t24 * t26 * t36 * t28;
    dydx[7] = -2.0 * t1 * y[7] * y[5] - 2.0 / t26 * t36 * y[7] * y[6];

    return true;
}

/*! Tests whether the constraint equation is fulfilled.
 *
 *  The constraint equation for lightlike and timelike geodesics reads:
 \verbatim
     sum = g_{\mu\nu} dot(x)^{\mu} dot(x)^{\nu} - kappa c^2 = 0.
 \endverbatim
 *  However, take care of the limited double precision.
 *  \param  y[]   : pointer to position and direction coordinates.
 *  \param  kappa : timelike (-1.0), lightlike (0.0).
 *  \return double : sum.
 */
double MetricCosmicStringSchwarzschild::testConstraint(const double y[], const double kappa)
{
    double r = y[1];
    double theta = y[2];
    double cm = 1.0 / mSpeedOfLight;

    // Scale the directions with the speed of light before doubling them !!
    double dt = y[4];
    double dr = y[5] * cm;
    double dth = y[6] * cm;
    double dph = y[7] * cm;

    double sum = -kappa;
    sum += -(1.0 - rs / r) * dt * dt + dr * dr / (1.0 - rs / r)
        + r * r * (dth * dth + mBeta * mBeta * sin(theta) * sin(theta) * dph * dph);
    return sum;
}

/*! Calculate the scalar product between u and v:  g_{ab}u^a v^b.
 *
 *  \param pos  :  pointer to position.
 *  \param u    :  pointer to vector.
 *  \param v    :  pointer to vector.
 *  \param prod :  reference to scalar product.
 *  \param preCalcMetric : calculate metric coefficients before evaluating the scalar product.
 *  \return true  : success.
 *  \return false : position is not valid.
 */
bool MetricCosmicStringSchwarzschild::calcProduct(
    const double* pos, const double* u, const double* v, double& prod, bool)
{
    prod = 0.0;
    if (breakCondition(pos)) {
        return false;
    }

    double r = pos[1];
    double theta = pos[2];
    prod = -mSpeedOfLight * mSpeedOfLight * (1.0 - rs / r) * u[0] * v[0] + u[1] * v[1] / (1.0 - rs / r)
        + r * r * (u[2] * v[2] + mBeta * mBeta * sin(theta) * sin(theta) * u[3] * v[3]);
    return true;
}

/*! Set parameter 'pName' to 'val'.
 *
 *  Set 'mass' parameter and adjust Schwarzschild radius  rs=2GM/c^2.
 */
bool MetricCosmicStringSchwarzschild::setParam(const char* pName, double val)
{
    Metric::setParam(pName, val);
    if (strcmp(pName, "mass") == 0) {
        mMass = val;
        rs = 2.0 * mGravConstant * mMass / (mSpeedOfLight * mSpeedOfLight);
    }
    else if (strcmp(pName, "beta") == 0) {
        mBeta = val;
    }
    return true;
}

/*! Transform point p to 2+1 coordinates.
 *
 *  \param  p  : point in proper metric coordinates.
 *  \param  cp : reference to transformed point.
 *  \return true : success.
 */
bool MetricCosmicStringSchwarzschild::transToTwoPlusOne(vec4 p, vec4& cp)
{
    vec4 tp;
    TransCoordinates::toCartesianCoord(mCoordType, p, tp);
    cp = vec4(tp[0], tp[1], tp[2], tp[0]);
    return true;
}

/*! Transform point p to embedding coordinates.
 *
 *  \param p  : point to be transformed.
 *  \param ep : reference to 'embedded' point.
 *  \return true : success.
 *  \return false : otherwise.
 */
bool MetricCosmicStringSchwarzschild::transToEmbedding(vec4 p, vec4& ep)
{
    vec4 cp;
    transToPseudoCart(p, cp);

    double r = p[1];
    double x = cp[1];
    double y = cp[2];
    double z;
    if (r >= rs) {
        z = embFunc(r);
        ep = vec4(p[0], x, y, z);
        return true;
    }
    return false;
}

/*! Set embedding parameters.
 *
 *  \param  name : embedding parameter name.
 *  \param  val  : embedding parameter value.
 *  \return true  : success.
 *  \return false : parameter not valid.
 */
bool MetricCosmicStringSchwarzschild::setEmbeddingParam(const char* name, double val)
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

/*! Test embedding parameters
 *  \return  true : all parameters are ok
 *  \return  false : at least one parameter had to be adjusted.
 */
bool MetricCosmicStringSchwarzschild::testEmbeddingParams()
{
    bool allOk = true;
    if (mEmb_rmin < rs) {
        mEmb_rmin = rs;
        allOk &= false;
    }
    if (mEmb_rmax < rs) {
        mEmb_rmax = rs;
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

unsigned int MetricCosmicStringSchwarzschild::getEmbeddingVertices(
    float*& verts, unsigned int*& indices, unsigned int& numElems, unsigned int& counter)
{
    m4d::SafeDelete<float>(verts);
    m4d::SafeDelete<unsigned int>(indices);

    testEmbeddingParams();
    mEmb_rstep = (mEmb_rmax - mEmb_rmin) / mEmb_r_num;
    mEmb_phistep = 2.0 * M_PI / mEmb_phi_num;

    numElems = mEmb_r_num;
    counter = mEmb_phi_num + 1;

    unsigned int numVerts = numElems * counter;
    int numInds = static_cast<int>(numElems * counter * 2);

    verts = new float[numVerts * 3];
    indices = new unsigned int[numInds];

    float* vptr = verts;
    unsigned int* iptr = indices;

    unsigned int vnum;

    double x, y, z, r, phi;
    for (unsigned int k = 0; k < counter; k++) {
        phi = k * mEmb_phistep;

        for (unsigned int j = 0; j < numElems; j++) {
            r = mEmb_rmin + j * mEmb_rstep;
            x = r * cos(phi);
            y = r * sin(phi);

            if (r >= rs) {
                z = embFunc(r);
                *(vptr++) = static_cast<float>(x);
                *(vptr++) = static_cast<float>(y);
                *(vptr++) = static_cast<float>(z);

                vnum = k * numElems + j;
                *(iptr++) = vnum;
                *(iptr++) = vnum + numElems;
            }
        }
    }

    return numVerts;
}

/*! Effective potential.
 *  \param pos : initial position.
 *  \param cdir : initial four-direction.
 *  \param type : geodesic type.
 *  \param x : abscissa value.
 *  \param val : reference to effective potential value.
 *  \return true : effective potential exists at x.
 */
bool MetricCosmicStringSchwarzschild::effPotentialValue(
    const vec4 pos, const vec4 cdir, enum_geodesic_type type, const double x, double& val)
{
    double kappa = 0.0;
    if (type == enum_geodesic_timelike) {
        kappa = -mSign;
    }

    if (pos[1] < rs + 1e-2 || x < rs + 1e-2) {
        return false;
    }

    double h = pos[1] * pos[1] * mBeta * mBeta * cdir[3];
    val = 0.5 * (1.0 - rs / x) * (h * h / (x * x * mBeta * mBeta) - kappa * mSpeedOfLight * mSpeedOfLight);
    return true;
}

/*! Totatl energy.
 *  \param pos : initial position.
 *  \param cdir : initial four-direction.
 *  \param x : abscissa value.
 *  \param val : reference to total energy value.
 *  \return true : effective potential exists at x.
 */
bool MetricCosmicStringSchwarzschild::totEnergy(const vec4 pos, const vec4 cdir, const double, double& val)
{
    if (pos[1] < rs + 1e-2) {
        return false;
    }

    // 1/2*k^2/c^2:
    val = 0.5 * (1.0 - rs / pos[1]) * (1.0 - rs / pos[1]) * mSpeedOfLight * mSpeedOfLight * cdir[0] * cdir[0];
    return true;
}

bool MetricCosmicStringSchwarzschild::report(const vec4 pos, const vec4 cdir, char*& text)
{
    SafeDelete<char>(text);

    std::stringstream ss;
    ss << "Report for SchwarzschildString metric\n\tcoordinates : (t,r,theta,phi)\n";
    ss << "------------------------------       ---------------------------------\n";
    ss << "  physical units ................................. no\n";
    ss.precision(DEF_FIXED_REPORT_PRECISION);
    ss.setf(std::ios::fixed);
    ss << "  Schwarzschild radius ........... r_s = 2GM/c^2 = " << rs << std::endl;
    ss << "  Photon orbit .............................r_po = " << 1.5 * rs << std::endl;
    ss << "  innermost stable circular orbit  r_isco = 3r_s = " << 3.0 * rs << std::endl;
    ss << "                                            beta = " << 0.5 << std::endl;

    double k = (1.0 - rs / pos[1]) * mSpeedOfLight * mSpeedOfLight * cdir[0];
    double h = pos[1] * pos[1] * mBeta * mBeta * cdir[3];
    ss.precision(DEF_FIXED_REPORT_PRECISION);
    ss << "  constant of motion ......................... k = " << k << std::endl;
    ss << "  constant of motion ......................... h = " << h << std::endl;

    ss << "  critical angle for null geodesic ..... ksiCrit = ";
    double ksicrit;
    if (calcKsiCrit(pos, ksicrit)) {
        ss << ksicrit * RAD_TO_DEG << std::endl;
    }
    else {
        ss << "not valid here.";
    }

    text = new char[ss.str().length() + 2];
    return CopyString(ss.str().c_str(), text);
}

// ***************************** specific public methods ***************************
/*!  Calculate the critical angle.
 *
 *   For an observer at distance \f$x_i=r_s/r_i\f$, r_i=pos[1], the black hole has an angular diameter of
 * \f$\xi_{\mbox{crit}}\f$ with \f[ \xi_{\mbox{crit}} = \arcsin\sqrt{\frac{27}{4}x_i^2(1-x_i)} \f] .
 *
 *  \param pos : current position
 *  \param ksicrit : reference to critical angle.
 */
bool MetricCosmicStringSchwarzschild::calcKsiCrit(const vec4 pos, double& ksicrit)
{
    if (pos[1] < rs) {
        return false;
    }

    double xi = rs / pos[1];
    ksicrit = asin(sqrt(6.75 * xi * xi * (1.0 - xi)));
    return true;
}

// ********************************* protected methods *****************************

void MetricCosmicStringSchwarzschild::setStandardValues()
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

double MetricCosmicStringSchwarzschild::embFunc(const double r)
{
    double z = 0.0;
    if (mBeta * mBeta < 1.0 && (r - rs >= 0.0)) {
        double w = sqrt(r / (r - rs) - mBeta * mBeta);
        double q = sqrt(1.0 - mBeta * mBeta);
        z = w * (r - rs) - 0.5 * rs / q * log((w - q) / (w + q));
    }
    return z;
}

} // end namespace m4d
