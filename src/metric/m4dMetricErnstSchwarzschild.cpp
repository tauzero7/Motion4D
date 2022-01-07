/**
 * @file    m4dMetricErnstSchwarzschild.cpp
 * @author  Thomas Mueller
 *
 * This file is part of the m4d-library.
 */

#include "m4dMetricErnstSchwarzschild.h"

namespace m4d {

MetricErnstSchwarzschild::MetricErnstSchwarzschild(double mass, double B)
{
    mMetricName = "ErnstSchwarzschild";
    setCoordType(enum_coordinate_spherical);

    mPhysicalUnits = enum_physical_constants_geom;
    mSpeedOfLight = 1.0;
    mGravConstant = 1.0;
    mDielectricPerm = 1.0;

    addParam("mass", mass);
    mMass = mass;
    addParam("B", B);
    mB = B;

    mSign = -1.0;
    mLocTeds.push_back(enum_nat_tetrad_static);
    setStandardValues();
}

MetricErnstSchwarzschild::~MetricErnstSchwarzschild() {}

// *********************************** public methods ******************************

bool MetricErnstSchwarzschild::calculateMetric(const double* pos)
{
    double r = pos[1];
    double theta = pos[2];

    double m = mMass;
    double B = mB;

    double t3 = 2.0 * m / r;
    double t4 = B * B;
    double t5 = r * r;
    double t6 = t4 * t5;
    double t7 = sin(theta);
    double t8 = t7 * t7;
    double t9 = t6 * t8;
    double t15 = t4 * t4;
    double t16 = t5 * t5;
    double t17 = t15 * t16;
    double t18 = t8 * t8;
    double t27 = 1 / (1.0 - t3);
    double t43 = pow(1.0 + t9, 2.0);

    g_compts[0][0] = 1.0 - t3 + 2.0 * t9 - 4.0 * t4 * r * t8 * m + t17 * t18 - 2.0 * t15 * t5 * r * t18 * m;
    g_compts[0][1] = 0.0;
    g_compts[0][2] = 0.0;
    g_compts[0][3] = 0.0;
    g_compts[1][0] = 0.0;
    g_compts[1][1] = -t27 - 2.0 * t6 * t8 * t27 - t17 * t18 * t27;
    g_compts[1][2] = 0.0;
    g_compts[1][3] = 0.0;
    g_compts[2][0] = 0.0;
    g_compts[2][1] = 0.0;
    g_compts[2][2] = -t5 - 2.0 * t4 * t16 * t8 - t15 * t16 * t5 * t18;
    g_compts[2][3] = 0.0;
    g_compts[3][0] = 0.0;
    g_compts[3][1] = 0.0;
    g_compts[3][2] = 0.0;
    g_compts[3][3] = -t5 * t8 / t43;

    return true;
}

bool MetricErnstSchwarzschild::calculateChristoffels(const double* pos)
{
    double r = pos[1];
    double theta = pos[2];

    double m = mMass;
    double B = mB;

    double t1 = B * B;
    double t2 = r * r;
    double t3 = t2 * r;
    double t5 = sin(theta);
    double t6 = t5 * t5;
    double t8 = 2.0 * t1 * t3 * t6;
    double t9 = t1 * t2;
    double t11 = t9 * t6 * m;
    double t13 = t8 - 3.0 * t11 + m;
    double t15 = r - 2.0 * m;
    double t18 = t9 * t6;
    double t19 = 1.0 + t18;
    double t20 = 1 / t19;
    double t24 = 1 / r;
    double t26 = cos(theta);
    double t28 = t5 * t26 * t20;
    double t31 = 1 / t15;
    double t33 = t24 * t20;
    double t34 = t13 * t31 * t33;
    double t36 = 2.0 * t9 * t28;
    double t48 = 3.0 * t18 + 1.0;
    double t50 = t48 * t24 * t20;
    double t51 = -1.0 + t18;
    double t53 = t51 * t24 * t20;
    double t56 = t51 * t26;
    double t59 = t56 * t20 / t5;
    double t61 = t1 * t1;
    double t62 = t2 * t2;
    double t64 = t6 * t6;
    double t67 = 1 / (1.0 + 2.0 * t18 + t61 * t62 * t64);
    double t70 = t19 * t19;
    double t72 = 1 / t70 / t19;

    christoffel[0][0][0] = 0.0;
    christoffel[0][0][1] = t13 * t15 / t3 * t20;
    christoffel[0][0][2] = 2.0 * t15 * t1 * t24 * t28;
    christoffel[0][0][3] = 0.0;
    christoffel[0][1][0] = t34;
    christoffel[0][1][1] = 0.0;
    christoffel[0][1][2] = 0.0;
    christoffel[0][1][3] = 0.0;
    christoffel[0][2][0] = t36;
    christoffel[0][2][1] = 0.0;
    christoffel[0][2][2] = 0.0;
    christoffel[0][2][3] = 0.0;
    christoffel[0][3][0] = 0.0;
    christoffel[0][3][1] = 0.0;
    christoffel[0][3][2] = 0.0;
    christoffel[0][3][3] = 0.0;
    christoffel[1][0][0] = t34;
    christoffel[1][0][1] = 0.0;
    christoffel[1][0][2] = 0.0;
    christoffel[1][0][3] = 0.0;
    christoffel[1][1][0] = 0.0;
    christoffel[1][1][1] = (t8 - 5.0 * t11 - m) * t31 * t33;
    christoffel[1][1][2] = -2.0 * t1 * r * t5 * t26 * t31 * t20;
    christoffel[1][1][3] = 0.0;
    christoffel[1][2][0] = 0.0;
    christoffel[1][2][1] = t36;
    christoffel[1][2][2] = t50;
    christoffel[1][2][3] = 0.0;
    christoffel[1][3][0] = 0.0;
    christoffel[1][3][1] = 0.0;
    christoffel[1][3][2] = 0.0;
    christoffel[1][3][3] = -t53;
    christoffel[2][0][0] = t36;
    christoffel[2][0][1] = 0.0;
    christoffel[2][0][2] = 0.0;
    christoffel[2][0][3] = 0.0;
    christoffel[2][1][0] = 0.0;
    christoffel[2][1][1] = t36;
    christoffel[2][1][2] = t50;
    christoffel[2][1][3] = 0.0;
    christoffel[2][2][0] = 0.0;
    christoffel[2][2][1] = -t48 * t15 * t20;
    christoffel[2][2][2] = t36;
    christoffel[2][2][3] = 0.0;
    christoffel[2][3][0] = 0.0;
    christoffel[2][3][1] = 0.0;
    christoffel[2][3][2] = 0.0;
    christoffel[2][3][3] = -t59;
    christoffel[3][0][0] = 0.0;
    christoffel[3][0][1] = 0.0;
    christoffel[3][0][2] = 0.0;
    christoffel[3][0][3] = 0.0;
    christoffel[3][1][0] = 0.0;
    christoffel[3][1][1] = 0.0;
    christoffel[3][1][2] = 0.0;
    christoffel[3][1][3] = -t53;
    christoffel[3][2][0] = 0.0;
    christoffel[3][2][1] = 0.0;
    christoffel[3][2][2] = 0.0;
    christoffel[3][2][3] = -t59;
    christoffel[3][3][0] = 0.0;
    christoffel[3][3][1] = t15 * t67 * t6 * t51 * t72;
    christoffel[3][3][2] = t67 * t5 * t56 * t72;
    christoffel[3][3][3] = 0.0;

    return true;
}

bool MetricErnstSchwarzschild::calculateChrisD(const double* pos)
{
    double r = pos[1];
    double theta = pos[2];

    double m = mMass;
    double B = mB;

    double t1 = B * B;
    double t2 = r * r;
    double t3 = t2 * t2;
    double t5 = sin(theta);
    double t6 = t5 * t5;
    double t7 = t1 * t3 * t6;
    double t8 = t1 * t1;
    double t11 = t6 * t6;
    double t12 = t8 * t3 * t2 * t11;
    double t13 = t2 * r;
    double t16 = t13 * m * t1 * t6;
    double t17 = 2.0 * t16;
    double t20 = t11 * m;
    double t21 = t8 * t3 * r * t20;
    double t23 = m * m;
    double t26 = t23 * t1 * t2 * t6;
    double t27 = 2.0 * t26;
    double t28 = t8 * t3;
    double t30 = t28 * t11 * t23;
    double t32 = r * m;
    double t37 = t1 * t2;
    double t38 = t37 * t6;
    double t40 = pow(1.0 + t38, 2.0);
    double t41 = 1 / t40;
    double t44 = 1 / r;
    double t47 = cos(theta);
    double t49 = r - 2.0 * m;
    double t50 = t49 * t49;
    double t55 = t47 * t5;
    double t56 = t55 * t1;
    double t58 = t1 * t13 * t6;
    double t59 = t6 * m;
    double t60 = t37 * t59;
    double t63 = 1 / t2;
    double t69 = t47 * t47;
    double t71 = t37 * t6 * t69;
    double t74 = -t69 + t71 + t6 + t11 * t1 * t2;
    double t84 = 1 / t50;
    double t86 = t63 * t41;
    double t88 = 2.0 * (-t7 + t12 + 6.0 * t16 - 3.0 * t21 - 6.0 * t26 + 3.0 * t30 + t32 - t23) * t84 * t86;
    double t89 = t1 * r;
    double t92 = 4.0 * t89 * t55 * t41;
    double t95 = 2.0 * t37 * t74 * t41;
    double t112 = t28 * t11;
    double t113 = 3.0 * t112;
    double t116 = (t113 + 1.0) * t63 * t41;
    double t120 = (-4.0 * t38 + t112 - 1.0) * t63 * t41;
    double t122 = t89 * t59;
    double t133 = t28 * t11 * t69;
    double t136 = t11 * t6 * t8 * t3;
    double t140 = (-4.0 * t71 + t133 - t6 + t136 - t69) * t41 / t6;
    double t149 = 2.0 * t38;
    double t151 = 1 / (1.0 + t149 + t112);
    double t152 = t40 * t40;
    double t153 = 1 / t152;
    double t154 = t151 * t153;

    chrisD[0][0][0][0] = 0.0;
    chrisD[0][0][0][1] = 0.0;
    chrisD[0][0][0][2] = 0.0;
    chrisD[0][0][0][3] = 0.0;
    chrisD[0][0][1][0] = 0.0;
    chrisD[0][0][1][1] = -2.0 * (-t7 + t12 + t17 - 7.0 * t21 - t27 + 9.0 * t30 + t32 - 3.0 * t23) / t3 * t41;
    chrisD[0][0][1][2] = 4.0 * t1 * t44 * t5 * t47 * t50 * t41;
    chrisD[0][0][1][3] = 0.0;
    chrisD[0][0][2][0] = 0.0;
    chrisD[0][0][2][1] = -4.0 * t56 * (t58 - m - 3.0 * t60) * t63 * t41;
    chrisD[0][0][2][2] = -2.0 * t49 * t1 * t74 * t44 * t41;
    chrisD[0][0][2][3] = 0.0;
    chrisD[0][0][3][0] = 0.0;
    chrisD[0][0][3][1] = 0.0;
    chrisD[0][0][3][2] = 0.0;
    chrisD[0][0][3][3] = 0.0;
    chrisD[0][1][0][0] = 0.0;
    chrisD[0][1][0][1] = -t88;
    chrisD[0][1][0][2] = t92;
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
    chrisD[0][2][0][1] = t92;
    chrisD[0][2][0][2] = -t95;
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
    chrisD[1][0][0][1] = -t88;
    chrisD[1][0][0][2] = t92;
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
    chrisD[1][1][1][1] = -2.0 * (-t7 + t12 + t17 - 5.0 * t21 - t27 + 5.0 * t30 - t32 + t23) * t84 * t86;
    chrisD[1][1][1][2] = t92;
    chrisD[1][1][1][3] = 0.0;
    chrisD[1][1][2][0] = 0.0;
    chrisD[1][1][2][1] = 4.0 * t56 * (t58 + m - t60) * t84 * t41;
    chrisD[1][1][2][2] = 2.0 * t89 * t74 / t49 * t41;
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
    chrisD[1][2][1][1] = t92;
    chrisD[1][2][1][2] = -t95;
    chrisD[1][2][1][3] = 0.0;
    chrisD[1][2][2][0] = 0.0;
    chrisD[1][2][2][1] = -t116;
    chrisD[1][2][2][2] = t92;
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
    chrisD[1][3][3][1] = t120;
    chrisD[1][3][3][2] = -t92;
    chrisD[1][3][3][3] = 0.0;
    chrisD[2][0][0][0] = 0.0;
    chrisD[2][0][0][1] = t92;
    chrisD[2][0][0][2] = -t95;
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
    chrisD[2][1][1][1] = t92;
    chrisD[2][1][1][2] = -t95;
    chrisD[2][1][1][3] = 0.0;
    chrisD[2][1][2][0] = 0.0;
    chrisD[2][1][2][1] = -t116;
    chrisD[2][1][2][2] = t92;
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
    chrisD[2][2][1][1] = -(8.0 * t38 + t113 - 8.0 * t122 + 1.0) * t41;
    chrisD[2][2][1][2] = -4.0 * t37 * t5 * t47 * t49 * t41;
    chrisD[2][2][1][3] = 0.0;
    chrisD[2][2][2][0] = 0.0;
    chrisD[2][2][2][1] = t92;
    chrisD[2][2][2][2] = -t95;
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
    chrisD[2][3][3][1] = -t92;
    chrisD[2][3][3][2] = t140;
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
    chrisD[3][1][3][1] = t120;
    chrisD[3][1][3][2] = -t92;
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
    chrisD[3][2][3][1] = -t92;
    chrisD[3][2][3][2] = t140;
    chrisD[3][2][3][3] = 0.0;
    chrisD[3][3][0][0] = 0.0;
    chrisD[3][3][0][1] = 0.0;
    chrisD[3][3][0][2] = 0.0;
    chrisD[3][3][0][3] = 0.0;
    chrisD[3][3][1][0] = 0.0;
    chrisD[3][3][1][1] = -(7.0 * t112 - 16.0 * t8 * t13 * t20 - 12.0 * t38 + 24.0 * t122 + 1.0) * t6 * t154;
    chrisD[3][3][1][2] = -2.0 * (t113 - 6.0 * t38 + 1.0) * t49 * t5 * t47 * t151 * t153;
    chrisD[3][3][1][3] = 0.0;
    chrisD[3][3][2][0] = 0.0;
    chrisD[3][3][2][1] = -4.0 * (t149 - 3.0) * t6 * t5 * t47 * t89 * t154;
    chrisD[3][3][2][2] = -(t136 + 7.0 * t133 - 12.0 * t71 + t69 - t6) * t151 * t153;
    chrisD[3][3][2][3] = 0.0;
    chrisD[3][3][3][0] = 0.0;
    chrisD[3][3][3][1] = 0.0;
    chrisD[3][3][3][2] = 0.0;
    chrisD[3][3][3][3] = 0.0;

    return true;
}

void MetricErnstSchwarzschild::localToCoord(const double* pos, const double* ldir, double* dir, enum_nat_tetrad_type)
{
    double r = pos[1];
    double theta = pos[2];

    double w = sqrt(1.0 - 2.0 * mMass / r);
    double L = 1.0 + mB * mB * r * r * sin(theta) * sin(theta);

    dir[0] = ldir[0] / L / w;
    dir[1] = ldir[1] * w / L;
    dir[2] = ldir[2] / (L * r);
    dir[3] = ldir[3] * L / (r * sin(theta));
}

void MetricErnstSchwarzschild::coordToLocal(const double* pos, const double* cdir, double* ldir, enum_nat_tetrad_type)
{
    double r = pos[1];
    double theta = pos[2];

    double w = sqrt(1.0 - 2.0 * mMass / r);
    double L = 1.0 + mB * mB * r * r * sin(theta) * sin(theta);

    ldir[0] = cdir[0] * L * w;
    ldir[1] = cdir[1] * L / w;
    ldir[2] = cdir[2] * L * r;
    ldir[3] = cdir[3] * r * sin(theta) / L;
}

bool MetricErnstSchwarzschild::breakCondition(const double* pos)
{
    bool br = false;
    if ((pos[1] < 0.0) || (pos[1] * pos[1] <= (1.0 + M4D_METRIC_EPS) * 4.0 * mMass * mMass)) {
        br = true;
    }

    return br;
}

double MetricErnstSchwarzschild::testConstraint(const double y[], const double kappa)
{
    double r = y[1];
    double theta = y[2];

    double w2 = 1.0 - 2.0 * mMass / r;
    double L = 1.0 + mB * mB * r * r * sin(theta) * sin(theta);

    // Scale the directions with the speed of light before doubling them !!
    double dt = y[4];
    double dr = y[5];
    double dtheta = y[6];
    double dphi = y[7];

    double sum = -kappa * mSign;
    sum += L * L * (w2 * dt * dt - dr * dr / w2 - r * r * dtheta * dtheta)
        - r * r * sin(theta) * sin(theta) / (L * L) * dphi * dphi;
    return sum;
}

bool MetricErnstSchwarzschild::setParam(const char* pName, double val)
{
    Metric::setParam(pName, val);

    if (strcmp(pName, "mass") == 0) {
        mMass = val;
    }
    else if (strcmp(pName, "b") == 0) {
        mB = val;
    }

    return true;
}

bool MetricErnstSchwarzschild::report(const vec4, const vec4, char*& text)
{
    std::stringstream ss;
    ss << "Report for Ernst-Schwarzschild metric\n\tcoordinate : (t,r,theta,phi)\n";
    ss << "---------------------------------------------------------------\n";
    ss.precision(DEF_FIXED_REPORT_PRECISION);
    ss.setf(std::ios::fixed);
    ss << "  mass ........................... " << mMass << std::endl;
    ss << "  B .............................. " << mB << std::endl;
    text = new char[ss.str().length() + 2];
    return CopyString(ss.str().c_str(), text);
}

// ********************************* protected methods *****************************
void MetricErnstSchwarzschild::setStandardValues()
{
    mInitPos[0] = 0.0;
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
