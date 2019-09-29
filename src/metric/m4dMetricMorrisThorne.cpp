/**
 * @file    m4dMetricMorrisThorne.cpp
 * @author  Thomas Mueller
 *
 *  This file is part of libMotion4D.
 */
#include "m4dMetricMorrisThorne.h"
#include "extra/m4dUtilities.h"

namespace m4d {

#define eps 1.0e-6

MetricMorrisThorne::MetricMorrisThorne(double b0)
{
    mMetricName = "MorrisThorne";
    setCoordType(enum_coordinate_spherical);

    mPhysicalUnits = enum_physical_constants_geom;
    mSpeedOfLight = 1.0;
    mGravConstant = 1.0;

    addParam("b0", b0);
    mb0 = b0;

    mLocTeds.push_back(enum_nat_tetrad_static);

    if (!mEmbParam.empty()) {
        mEmbParam.clear();
    }
    mHaveEmbedding = true;

    mEmb_lmin = -10.0;
    mEmb_lmax = 10.0;
    mEmb_l_num = 40.0;
    mEmb_phi_num = 40.0;
    mEmb_lstep = (mEmb_lmax - mEmb_lmin) / mEmb_l_num;
    mEmb_phistep = 2.0 * M_PI / mEmb_phi_num;
    addEmbeddingParam("emb_lmin", mEmb_lmin);
    addEmbeddingParam("emb_lmax", mEmb_lmax);
    addEmbeddingParam("emb_l_num", mEmb_l_num);
    addEmbeddingParam("emb_phi_num", mEmb_phi_num);

    mDrawTypes.push_back(enum_draw_embedding);
    mDrawTypes.push_back(enum_draw_effpoti);

    setStandardValues();
}

MetricMorrisThorne::~MetricMorrisThorne() {}

bool MetricMorrisThorne::calculateMetric(const double* pos)
{
    double l = pos[1];
    double theta = pos[2];

    double c = mSpeedOfLight;

    double t1 = c * c;
    double t2 = l * l;
    double t3 = mb0 * mb0;
    double t5 = sin(theta);
    double t6 = t5 * t5;

    g_compts[0][0] = -t1;
    g_compts[0][1] = 0.0;
    g_compts[0][2] = 0.0;
    g_compts[0][3] = 0.0;
    g_compts[1][0] = 0.0;
    g_compts[1][1] = 1.0;
    g_compts[1][2] = 0.0;
    g_compts[1][3] = 0.0;
    g_compts[2][0] = 0.0;
    g_compts[2][1] = 0.0;
    g_compts[2][2] = t2 + t3;
    g_compts[2][3] = 0.0;
    g_compts[3][0] = 0.0;
    g_compts[3][1] = 0.0;
    g_compts[3][2] = 0.0;
    g_compts[3][3] = t2 * t6 + t3 * t6;

    return true;
}

bool MetricMorrisThorne::calculateChristoffels(const double* pos)
{
    double l = pos[1];
    double theta = pos[2];

    double t1 = l * l;
    double t2 = mb0 * mb0;
    double t5 = 1 / (t1 + t2) * l;
    double t6 = sin(theta);
    double t8 = cos(theta);
    double t9 = 1 / t6 * t8;
    double t10 = t6 * t6;

    christoffel[0][0][0] = 0.0;
    christoffel[0][0][1] = 0.0;
    christoffel[0][0][2] = 0.0;
    christoffel[0][0][3] = 0.0;
    christoffel[0][1][0] = 0.0;
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
    christoffel[1][0][0] = 0.0;
    christoffel[1][0][1] = 0.0;
    christoffel[1][0][2] = 0.0;
    christoffel[1][0][3] = 0.0;
    christoffel[1][1][0] = 0.0;
    christoffel[1][1][1] = 0.0;
    christoffel[1][1][2] = 0.0;
    christoffel[1][1][3] = 0.0;
    christoffel[1][2][0] = 0.0;
    christoffel[1][2][1] = 0.0;
    christoffel[1][2][2] = t5;
    christoffel[1][2][3] = 0.0;
    christoffel[1][3][0] = 0.0;
    christoffel[1][3][1] = 0.0;
    christoffel[1][3][2] = 0.0;
    christoffel[1][3][3] = t5;
    christoffel[2][0][0] = 0.0;
    christoffel[2][0][1] = 0.0;
    christoffel[2][0][2] = 0.0;
    christoffel[2][0][3] = 0.0;
    christoffel[2][1][0] = 0.0;
    christoffel[2][1][1] = 0.0;
    christoffel[2][1][2] = t5;
    christoffel[2][1][3] = 0.0;
    christoffel[2][2][0] = 0.0;
    christoffel[2][2][1] = -l;
    christoffel[2][2][2] = 0.0;
    christoffel[2][2][3] = 0.0;
    christoffel[2][3][0] = 0.0;
    christoffel[2][3][1] = 0.0;
    christoffel[2][3][2] = 0.0;
    christoffel[2][3][3] = t9;
    christoffel[3][0][0] = 0.0;
    christoffel[3][0][1] = 0.0;
    christoffel[3][0][2] = 0.0;
    christoffel[3][0][3] = 0.0;
    christoffel[3][1][0] = 0.0;
    christoffel[3][1][1] = 0.0;
    christoffel[3][1][2] = 0.0;
    christoffel[3][1][3] = t5;
    christoffel[3][2][0] = 0.0;
    christoffel[3][2][1] = 0.0;
    christoffel[3][2][2] = 0.0;
    christoffel[3][2][3] = t9;
    christoffel[3][3][0] = 0.0;
    christoffel[3][3][1] = -l * t10;
    christoffel[3][3][2] = -t6 * t8;
    christoffel[3][3][3] = 0.0;

    return true;
}

bool MetricMorrisThorne::calculateChrisD(const double* pos)
{
    double l = pos[1];
    double theta = pos[2];

    double t1 = l * l;
    double t2 = mb0 * mb0;
    double t5 = pow(t1 + t2, 2.0);
    double t7 = (t1 - t2) / t5;
    double t8 = cos(theta);
    double t9 = t8 * t8;
    double t10 = sin(theta);
    double t11 = t10 * t10;
    double t14 = (t9 + t11) / t11;

    chrisD[0][0][0][0] = 0.0;
    chrisD[0][0][0][1] = 0.0;
    chrisD[0][0][0][2] = 0.0;
    chrisD[0][0][0][3] = 0.0;
    chrisD[0][0][1][0] = 0.0;
    chrisD[0][0][1][1] = 0.0;
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
    chrisD[0][1][0][1] = 0.0;
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
    chrisD[1][0][0][1] = 0.0;
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
    chrisD[1][1][1][1] = 0.0;
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
    chrisD[1][2][2][1] = -t7;
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
    chrisD[1][3][3][1] = -t7;
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
    chrisD[2][1][2][1] = -t7;
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
    chrisD[2][3][3][2] = -t14;
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
    chrisD[3][1][3][1] = -t7;
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
    chrisD[3][2][3][2] = -t14;
    chrisD[3][2][3][3] = 0.0;
    chrisD[3][3][0][0] = 0.0;
    chrisD[3][3][0][1] = 0.0;
    chrisD[3][3][0][2] = 0.0;
    chrisD[3][3][0][3] = 0.0;
    chrisD[3][3][1][0] = 0.0;
    chrisD[3][3][1][1] = -t11;
    chrisD[3][3][1][2] = -2.0 * l * t10 * t8;
    chrisD[3][3][1][3] = 0.0;
    chrisD[3][3][2][0] = 0.0;
    chrisD[3][3][2][1] = 0.0;
    chrisD[3][3][2][2] = -t9 + t11;
    chrisD[3][3][2][3] = 0.0;
    chrisD[3][3][3][0] = 0.0;
    chrisD[3][3][3][1] = 0.0;
    chrisD[3][3][3][2] = 0.0;
    chrisD[3][3][3][3] = 0.0;

    return true;
}

void MetricMorrisThorne::localToCoord(const double* pos, const double* ldir, double* dir, enum_nat_tetrad_type)
{
    double w = sqrt(mb0 * mb0 + pos[1] * pos[1]);

    dir[0] = ldir[0] / mSpeedOfLight;
    dir[1] = ldir[1];
    dir[2] = ldir[2] / w;
    dir[3] = ldir[3] / (w * sin(pos[2]));
}

void MetricMorrisThorne::coordToLocal(const double* pos, const double* cdir, double* ldir, enum_nat_tetrad_type)
{
    double w = sqrt(mb0 * mb0 + pos[1] * pos[1]);

    ldir[0] = cdir[0] * mSpeedOfLight;
    ldir[1] = cdir[1];
    ldir[2] = cdir[2] * w;
    ldir[3] = cdir[3] * w * sin(pos[2]);
}

bool MetricMorrisThorne::breakCondition(const double*)
{
    bool br = false;
    return br;
}

int MetricMorrisThorne::transToPseudoCart(vec4 p, vec4& cp)
{
    TransCoordinates::toCartesianCoord(mCoordType, p, cp);
    if (p[1] > 0) {
        return 0;
    }
    return 1;
}

bool MetricMorrisThorne::calcDerivs(const double y[], double dydx[])
{
    dydx[0] = y[4];
    dydx[1] = y[5];
    dydx[2] = y[6];
    dydx[3] = y[7];

    double l = y[1];
    double theta = y[2];
    double st = sin(theta);
    double ct = cos(theta);

    dydx[4] = 0.0;
    dydx[5] = l * y[6] * y[6] + l * st * st * y[7] * y[7];
    dydx[6] = -2.0 * l / (mb0 * mb0 + l * l) * y[5] * y[6] + st * ct * y[7] * y[7];
    dydx[7] = -2.0 * l / (mb0 * mb0 + l * l) * y[5] * y[7] - 2.0 * ct / st * y[6] * y[7];

    return true;
}

double MetricMorrisThorne::testConstraint(const double y[], const double kappa)
{
    double l = y[1];
    double cm = 1.0 / mSpeedOfLight;

    // Scale the directions with the speed of light before doubling them !!
    double dt = y[4];
    double dl = y[5] * cm;
    double dth = y[6] * cm;
    double dph = y[7] * cm;
    double st = sin(y[2]);

    double sum = -kappa;
    sum += -dt * dt + dl * dl + (mb0 * mb0 + l * l) * (dth * dth + st * st * dph * dph);

    return sum;
}

bool MetricMorrisThorne::setParam(const char* pName, double val)
{
    if (Metric::setParam(pName, val)) {
        mb0 = val;
    }

    return true;
}

bool MetricMorrisThorne::transToEmbedding(vec4 p, vec4& ep)
{
    if (mb0 <= 0.0) {
        return false;
    }

    vec4 cp;
    vec4 propCoords = p;
    propCoords[1] = sqrt(mb0 * mb0 + p[1] * p[1]);
    transToPseudoCart(propCoords, cp);

    double r = propCoords[1];
    double x = cp[1];
    double y = cp[2];
    double z;

    z = mb0 * log(r / mb0 + sqrt(r * r / (mb0 * mb0) - 1.0));
    if (p[1] < 0.0) {
        z *= -1.0;
    }

    ep = vec4(p[0], x, y, z);
    return true;
}

bool MetricMorrisThorne::transToCustom(vec4, vec4&)
{
    //  TODO
    return false;
}

bool MetricMorrisThorne::setEmbeddingParam(const char* name, double val)
{
    Metric::setEmbeddingParam(name, val);

    if (strcmp(name, "emb_lmin") == 0) {
        mEmb_lmin = val;
    }
    else if (strcmp(name, "emb_lmax") == 0) {
        mEmb_lmax = val;
    }
    else if (strcmp(name, "emb_l_num") == 0) {
        mEmb_l_num = static_cast<unsigned int>(std::max(val, 5.0));
    }
    else if (strcmp(name, "emb_phi_num") == 0) {
        mEmb_phi_num = static_cast<unsigned int>(std::max(val, 4.0));
    }
    return true;
}

unsigned int MetricMorrisThorne::getEmbeddingVertices(
    float*& verts, unsigned int*& indices, unsigned int& numElems, unsigned int& counter)
{
    m4d::SafeDelete<float>(verts);
    m4d::SafeDelete<unsigned int>(indices);

    mEmb_lstep = (mEmb_lmax - mEmb_lmin) / static_cast<double>(mEmb_l_num);
    mEmb_phistep = 2.0 * M_PI / mEmb_phi_num;

    numElems = mEmb_l_num;
    counter = mEmb_phi_num + 1;

    unsigned int numVerts = numElems * counter;
    int numInds = static_cast<int>(numElems * counter * 2);

    verts = new float[numVerts * 3];
    indices = new unsigned int[numInds];

    float* vptr = verts;
    unsigned int* iptr = indices;

    unsigned int vnum;

    double x, y, z, r, phi, l;
    for (unsigned int k = 0; k < counter; k++) {
        phi = k * mEmb_phistep;

        for (unsigned int j = 0; j < numElems; j++) {
            l = mEmb_lmin + j * mEmb_lstep;
            r = sqrt(mb0 * mb0 + l * l);
            x = r * cos(phi);
            y = r * sin(phi);
            z = mb0 * log(r / mb0 + sqrt(r * r / (mb0 * mb0) - 1.0));

            if (l < 0.0) {
                z *= -1.0;
            }

            *(vptr++) = static_cast<float>(x);
            *(vptr++) = static_cast<float>(y);
            *(vptr++) = static_cast<float>(z);

            vnum = k * numElems + j;
            *(iptr++) = vnum;
            *(iptr++) = vnum + numElems;
        }
    }

    return numVerts;
}

bool MetricMorrisThorne::effPotentialValue(
    const vec4 pos, const vec4 cdir, enum_geodesic_type type, const double x, double& val)
{
    double kappa = 0.0;
    if (type == enum_geodesic_timelike) {
        kappa = -mSign;
    }

    double h = (mb0 * mb0 + pos[1] * pos[1]) * cdir[3];
    val = 0.5 * (h * h / (mb0 * mb0 + x * x) - kappa * mSpeedOfLight * mSpeedOfLight);
    return true;
}

bool MetricMorrisThorne::totEnergy(const vec4, const vec4 cdir, const double, double& val)
{
    // 1/2*k^2/c^2:
    val = 0.5 * mSpeedOfLight * mSpeedOfLight * cdir[0] * cdir[0];
    return true;
}

bool MetricMorrisThorne::report(const vec4 pos, const vec4 cdir, char*& text)
{
    std::stringstream ss;
    ss << "Report for Morris-Thorne metric\n\tcoordinate : (t,l,theta,phi)\n";
    ss << "---------------------------------------------------------------\n";
    ss << "  Throat size .............................. b_0 = " << mb0 << std::endl;

    double k, h;
    constsOfMotion(pos, cdir, k, h);
    ss.precision(10);
    ss << "  constant of motion ......................... k = " << k << std::endl;
    ss << "  constant of motion ......................... h = " << h << std::endl;

    double ksicrit;
    calcKsiCrit(pos, ksicrit);
    ss << "  critical angle .................. ksiCrit(deg) = " << ksicrit * RAD_TO_DEG << std::endl;

    text = new char[ss.str().length() + 2];
    return CopyString(ss.str().c_str(), text);
}

bool MetricMorrisThorne::calcKsiCrit(const vec4 pos, double& ksicrit)
{
    double l = pos[1];
    ksicrit = asin(mb0 / sqrt(mb0 * mb0 + l * l));
    return true;
}

bool MetricMorrisThorne::constsOfMotion(const vec4 pos, const vec4 cdir, double& k, double& h)
{
    k = mSpeedOfLight * mSpeedOfLight * cdir[0];
    h = (mb0 * mb0 + pos[1] * pos[1]) * cdir[3];
    return true;
}

void MetricMorrisThorne::setStandardValues()
{
    mInitPos[0] = 0.0;
    mInitPos[1] = 10.0;
    mInitPos[2] = M_PI_2;
    mInitPos[3] = 0.0;
    mInitDir[0] = 1.0;
    mInitDir[1] = 0.0;
    mInitDir[2] = 0.0;

    mCoordNames[0] = std::string("t");
    mCoordNames[1] = std::string("l");
    mCoordNames[2] = std::string("theta");
    mCoordNames[3] = std::string("phi");
}

} // end namespace m4d
