/**
 * @file    m4dMetricSultanaDyer.cpp
 * @author  Thomas Mueller
 *
 * This file is part of the m4d-library.
 */
#include "m4dMetricSultanaDyer.h"

namespace m4d {

MetricSultanaDyer::MetricSultanaDyer(double mass)
{
    mMetricName = "SultanaDyerBlackhole";
    setCoordType(enum_coordinate_spherical);

    mPhysicalUnits = enum_physical_constants_geom;
    mSpeedOfLight = 1.0;
    mGravConstant = 1.0;

    mMass = mass;
    addParam("mass", mMass);

    mSign = -1.0; // !!

    mLocTeds.push_back(enum_nat_tetrad_static);
    mLocTeds.push_back(enum_nat_tetrad_comoving);

    setStandardValues();
}

MetricSultanaDyer::~MetricSultanaDyer() {}

// *********************************** public methods ******************************

bool MetricSultanaDyer::calculateMetric(const double* pos)
{
    double t = pos[0];
    double r = pos[1];
    double theta = pos[2];

    double m = mMass;

    double t1 = t * t;
    double t2 = t1 * t1;
    double t6 = 2.0 * t2 * m / r;
    double t9 = r * r;
    double t10 = t2 * t9;
    double t11 = sin(theta);
    double t12 = t11 * t11;

    g_compts[0][0] = t2 - t6;
    g_compts[0][1] = -t6;
    g_compts[0][2] = 0.0;
    g_compts[0][3] = 0.0;
    g_compts[1][0] = -t6;
    g_compts[1][1] = -t2 - t6;
    g_compts[1][2] = 0.0;
    g_compts[1][3] = 0.0;
    g_compts[2][0] = 0.0;
    g_compts[2][1] = 0.0;
    g_compts[2][2] = -t10;
    g_compts[2][3] = 0.0;
    g_compts[3][0] = 0.0;
    g_compts[3][1] = 0.0;
    g_compts[3][2] = 0.0;
    g_compts[3][3] = -t10 * t12;

    return true;
}

bool MetricSultanaDyer::calculateChristoffels(const double* pos)
{
    double t = pos[0];
    double r = pos[1];
    double theta = pos[2];

    double m = mMass;

    double t1 = r * r;
    double t2 = t1 * r;
    double t3 = m * m;
    double t5 = 4.0 * r * t3;
    double t6 = t3 * t;
    double t8 = 1 / t;
    double t10 = 1 / t2;
    double t13 = 2.0 * m;
    double t19 = (4.0 * r + t) * t8 * t10;
    double t23 = (r + t13) * m * t19;
    double t27 = 2.0 * (t6 - t2 + t5) * t8 * t10;
    double t28 = 2.0 * t8;
    double t29 = t * r;
    double t37 = t * m;
    double t38 = 2.0 * t37;
    double t40 = r * m;
    double t46 = 1 / r;
    double t48 = -t1 - 2.0 * t40 + t37;
    double t52 = 4.0 * t40 + t29 - t38;
    double t54 = sin(theta);
    double t56 = cos(theta);
    double t57 = 1 / t54 * t56;
    double t58 = t54 * t54;

    christoffel[0][0][0] = 2.0 * (t2 + t5 + t6) * t8 * t10;
    christoffel[0][0][1] = m * (r - t13) * t19;
    christoffel[0][0][2] = 0.0;
    christoffel[0][0][3] = 0.0;
    christoffel[0][1][0] = t23;
    christoffel[0][1][1] = -t27;
    christoffel[0][1][2] = 0.0;
    christoffel[0][1][3] = 0.0;
    christoffel[0][2][0] = 0.0;
    christoffel[0][2][1] = 0.0;
    christoffel[0][2][2] = t28;
    christoffel[0][2][3] = 0.0;
    christoffel[0][3][0] = 0.0;
    christoffel[0][3][1] = 0.0;
    christoffel[0][3][2] = 0.0;
    christoffel[0][3][3] = t28;
    christoffel[1][0][0] = t23;
    christoffel[1][0][1] = -t27;
    christoffel[1][0][2] = 0.0;
    christoffel[1][0][3] = 0.0;
    christoffel[1][1][0] = 2.0 * (t29 * m + t2 + 4.0 * t1 * m + t6 + t5) * t8 * t10;
    christoffel[1][1][1] = -m * (t38 + 4.0 * t1 + 8.0 * t40 + t29) * t8 * t10;
    christoffel[1][1][2] = 0.0;
    christoffel[1][1][3] = 0.0;
    christoffel[1][2][0] = 0.0;
    christoffel[1][2][1] = 0.0;
    christoffel[1][2][2] = t46;
    christoffel[1][2][3] = 0.0;
    christoffel[1][3][0] = 0.0;
    christoffel[1][3][1] = 0.0;
    christoffel[1][3][2] = 0.0;
    christoffel[1][3][3] = t46;
    christoffel[2][0][0] = 0.0;
    christoffel[2][0][1] = 0.0;
    christoffel[2][0][2] = t28;
    christoffel[2][0][3] = 0.0;
    christoffel[2][1][0] = 0.0;
    christoffel[2][1][1] = 0.0;
    christoffel[2][1][2] = t46;
    christoffel[2][1][3] = 0.0;
    christoffel[2][2][0] = -2.0 * t48 * t8;
    christoffel[2][2][1] = -t52 * t8;
    christoffel[2][2][2] = 0.0;
    christoffel[2][2][3] = 0.0;
    christoffel[2][3][0] = 0.0;
    christoffel[2][3][1] = 0.0;
    christoffel[2][3][2] = 0.0;
    christoffel[2][3][3] = t57;
    christoffel[3][0][0] = 0.0;
    christoffel[3][0][1] = 0.0;
    christoffel[3][0][2] = 0.0;
    christoffel[3][0][3] = t28;
    christoffel[3][1][0] = 0.0;
    christoffel[3][1][1] = 0.0;
    christoffel[3][1][2] = 0.0;
    christoffel[3][1][3] = t46;
    christoffel[3][2][0] = 0.0;
    christoffel[3][2][1] = 0.0;
    christoffel[3][2][2] = 0.0;
    christoffel[3][2][3] = t57;
    christoffel[3][3][0] = -2.0 * t58 * t48 * t8;
    christoffel[3][3][1] = -t58 * t52 * t8;
    christoffel[3][3][2] = -t54 * t56;
    christoffel[3][3][3] = 0.0;

    return true;
}

bool MetricSultanaDyer::calculateChrisD(const double* pos)
{
    double t = pos[0];
    double r = pos[1];
    double theta = pos[2];

    double m = mMass;

    double t1 = r * r;
    double t2 = 1 / t1;
    double t3 = m * m;
    double t4 = 4.0 * t3;
    double t7 = t * t;
    double t8 = 1 / t7;
    double t15 = 1 / t;
    double t16 = t1 * t1;
    double t18 = t15 / t16;
    double t20 = 2.0 * t3 * (8.0 * r + 3.0 * t) * t18;
    double t21 = 2.0 * m;
    double t24 = t2 * t8;
    double t27 = 2.0 * t1;
    double t28 = t * r;
    double t29 = r * m;
    double t30 = 8.0 * t29;
    double t31 = t * m;
    double t32 = 3.0 * t31;
    double t37 = r + t21;
    double t40 = 4.0 * t37 * m * t24;
    double t44 = 2.0 * m * (t27 + t28 + t30 + t32) * t18;
    double t48 = 2.0 * t2 * (t1 - t4) * t8;
    double t49 = 2.0 * t8;
    double t50 = 4.0 * t29;
    double t64 = r + m;
    double t70 = 4.0 * m + t;
    double t72 = cos(theta);
    double t73 = t72 * t72;
    double t74 = sin(theta);
    double t75 = t74 * t74;
    double t78 = (t73 + t75) / t75;
    double t89 = t15 * t72;

    chrisD[0][0][0][0] = -2.0 * t2 * (t1 + t4) * t8;
    chrisD[0][0][0][1] = -t20;
    chrisD[0][0][0][2] = 0.0;
    chrisD[0][0][0][3] = 0.0;
    chrisD[0][0][1][0] = -4.0 * m * (r - t21) * t24;
    chrisD[0][0][1][1] = -2.0 * m * (t27 + t28 - t30 - t32) * t18;
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
    chrisD[0][1][0][0] = -t40;
    chrisD[0][1][0][1] = -t44;
    chrisD[0][1][0][2] = 0.0;
    chrisD[0][1][0][3] = 0.0;
    chrisD[0][1][1][0] = -t48;
    chrisD[0][1][1][1] = t20;
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
    chrisD[0][2][2][0] = -t49;
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
    chrisD[0][3][3][0] = -t49;
    chrisD[0][3][3][1] = 0.0;
    chrisD[0][3][3][2] = 0.0;
    chrisD[0][3][3][3] = 0.0;
    chrisD[1][0][0][0] = -t40;
    chrisD[1][0][0][1] = -t44;
    chrisD[1][0][0][2] = 0.0;
    chrisD[1][0][0][3] = 0.0;
    chrisD[1][0][1][0] = -t48;
    chrisD[1][0][1][1] = t20;
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
    chrisD[1][1][0][0] = -2.0 * t2 * (t1 + t50 + t4) * t8;
    chrisD[1][1][0][1] = -2.0 * m * (2.0 * t28 + 4.0 * t1 + t30 + t32) * t18;
    chrisD[1][1][0][2] = 0.0;
    chrisD[1][1][0][3] = 0.0;
    chrisD[1][1][1][0] = t40;
    chrisD[1][1][1][1] = t44;
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
    chrisD[1][2][2][1] = -t2;
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
    chrisD[1][3][3][1] = -t2;
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
    chrisD[2][0][2][0] = -t49;
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
    chrisD[2][1][2][1] = -t2;
    chrisD[2][1][2][2] = 0.0;
    chrisD[2][1][2][3] = 0.0;
    chrisD[2][1][3][0] = 0.0;
    chrisD[2][1][3][1] = 0.0;
    chrisD[2][1][3][2] = 0.0;
    chrisD[2][1][3][3] = 0.0;
    chrisD[2][2][0][0] = -2.0 * t37 * r * t8;
    chrisD[2][2][0][1] = 4.0 * t64 * t15;
    chrisD[2][2][0][2] = 0.0;
    chrisD[2][2][0][3] = 0.0;
    chrisD[2][2][1][0] = 4.0 * t29 * t8;
    chrisD[2][2][1][1] = -t70 * t15;
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
    chrisD[2][3][3][2] = -t78;
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
    chrisD[3][0][3][0] = -t49;
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
    chrisD[3][1][3][1] = -t2;
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
    chrisD[3][2][3][2] = -t78;
    chrisD[3][2][3][3] = 0.0;
    chrisD[3][3][0][0] = -2.0 * t75 * t37 * r * t8;
    chrisD[3][3][0][1] = 4.0 * t75 * t64 * t15;
    chrisD[3][3][0][2] = -4.0 * t74 * (-t1 - 2.0 * t29 + t31) * t89;
    chrisD[3][3][0][3] = 0.0;
    chrisD[3][3][1][0] = 4.0 * t75 * r * m * t8;
    chrisD[3][3][1][1] = -t75 * t70 * t15;
    chrisD[3][3][1][2] = -2.0 * t74 * (t50 + t28 - 2.0 * t31) * t89;
    chrisD[3][3][1][3] = 0.0;
    chrisD[3][3][2][0] = 0.0;
    chrisD[3][3][2][1] = 0.0;
    chrisD[3][3][2][2] = -t73 + t75;
    chrisD[3][3][2][3] = 0.0;
    chrisD[3][3][3][0] = 0.0;
    chrisD[3][3][3][1] = 0.0;
    chrisD[3][3][3][2] = 0.0;
    chrisD[3][3][3][3] = 0.0;

    return true;
}

void MetricSultanaDyer::localToCoord(const double* pos, const double* ldir, double* dir, enum_nat_tetrad_type type)
{
    double t = pos[0];
    double r = pos[1];
    double theta = pos[2];

    double dtt = 1.0 / (t * t);

    if (type == enum_nat_tetrad_comoving) {
        double A = sqrt(1.0 + 2.0 * mMass / r) * dtt;
        double B = -2.0 * mMass / (r * sqrt(1.0 + 2.0 * mMass / r)) * dtt;
        double C = 1.0 / sqrt(1.0 + 2.0 * mMass / r) * dtt;

        dir[0] = ldir[0] * A;
        dir[1] = ldir[0] * B + ldir[1] * C;
        dir[2] = ldir[2] / r * dtt;
        dir[3] = ldir[3] / (r * sin(theta)) * dtt;
    }
    else {
        double A = 1.0 / sqrt(1.0 - 2.0 * mMass / r) * dtt;
        double B = 2.0 * mMass / (r * sqrt(1.0 - 2.0 * mMass / r)) * dtt;
        double C = sqrt(1.0 - 2.0 * mMass / r) * dtt;

        dir[0] = ldir[0] * A + ldir[1] * B;
        dir[1] = ldir[1] * C;
        dir[2] = ldir[2] / r * dtt;
        dir[3] = ldir[3] / (r * sin(theta)) * dtt;
    }
}

void MetricSultanaDyer::coordToLocal(const double* pos, const double* cdir, double* ldir, enum_nat_tetrad_type type)
{
    double t = pos[0];
    double r = pos[1];
    double theta = pos[2];

    double dtt = 1.0 / (t * t);

    if (type == enum_nat_tetrad_comoving) {
        double A = sqrt(1.0 + 2.0 * mMass / r) * dtt;
        double B = -2.0 * mMass / (r * sqrt(1.0 + 2.0 * mMass / r)) * dtt;
        double C = 1.0 / sqrt(1.0 + 2.0 * mMass / r) * dtt;

        ldir[0] = cdir[0] / A;
        ldir[1] = (cdir[1] - B / A * cdir[0]) / C;
        ldir[2] = cdir[2] * r * t * t;
        ldir[3] = cdir[3] * r * sin(theta) * t * t;
    }
    else {
        double A = 1.0 / sqrt(1.0 - 2.0 * mMass / r) * dtt;
        double B = 2.0 * mMass / (r * sqrt(1.0 - 2.0 * mMass / r)) * dtt;
        double C = sqrt(1.0 - 2.0 * mMass / r) * dtt;

        ldir[0] = (cdir[0] - cdir[1] * B / C) / A;
        ldir[1] = cdir[1] / C;
        ldir[2] = cdir[2] * r * t * t;
        ldir[3] = cdir[3] * r * sin(theta) * t * t;
    }
}

bool MetricSultanaDyer::breakCondition(const double* pos)
{
    bool br = false;

    double rs = 2.0 * mMass;
    if ((pos[1] < 0.0) || (pos[1] * pos[1] <= (1.0 + M4D_METRIC_EPS) * rs * rs)) {
        br = true;
    }
    return br;
}

bool MetricSultanaDyer::setParam(const char* pName, double val)
{
    Metric::setParam(pName, val);
    if (strcmp(pName, "mass") == 0) {
        mMass = val;
    }
    return true;
}

double MetricSultanaDyer::testConstraint(const double y[], const double kappa)
{
    double t = y[0];
    double r = y[1];
    double theta = y[2];
    double dt = y[4];
    double dr = y[5];
    double dtheta = y[6];
    double dphi = y[7];

    double sum = -kappa * mSign;
    sum += t * t * t * t
        * ((1.0 - 2.0 * mMass / r) * dt * dt - 4.0 * mMass / r * dt * dr - (1.0 + 2.0 * mMass / r) * dr * dr
            - r * r * (dtheta * dtheta + sin(theta) * sin(theta) * dphi * dphi));
    return sum;
}

bool MetricSultanaDyer::report(const vec4, const vec4, char*& text)
{
    std::stringstream ss;
    ss << "Report for SultanaDyer cosmological black hole\n\tcoordinates : (t,r,theta,phi)\n";
    ss << "---------------------------------------------------------------\n";
    ss << "  physical units ................................. no\n";
    ss.precision(DEF_FIXED_REPORT_PRECISION);
    ss.setf(std::ios::fixed);
    ss << "  mass ........................................ = " << mMass << std::endl;

    text = new char[ss.str().length() + 2];
    return CopyString(ss.str().c_str(), text);
}

// ********************************* protected methods *****************************

void MetricSultanaDyer::setStandardValues()
{
    mInitPos[0] = 1.0;
    mInitPos[1] = 6.0;
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
