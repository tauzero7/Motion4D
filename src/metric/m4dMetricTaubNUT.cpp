/**
 * @file    m4dMetricTaubNUT.cpp
 * @author  Thomas Mueller
 *
 * This file is part of the m4d-library.
 */
#include "m4dMetricTaubNUT.h"

namespace m4d {

MetricTaubNUT::MetricTaubNUT(double mass, double l)
{
    mMetricName = "TaubNUT";
    setCoordType(enum_coordinate_spherical);

    mPhysicalUnits = enum_physical_constants_geom;
    mSpeedOfLight = 1.0;
    mGravConstant = 1.0;

    addParam("mass");
    setParam("mass", mass);
    mMass = mass;

    addParam("l");
    setParam("l", l);
    mL = l;

    mLocTeds.push_back(enum_nat_tetrad_static);

    mDrawTypes.push_back(enum_draw_effpoti);

    setStandardValues();
}

MetricTaubNUT::~MetricTaubNUT() {}

// *********************************** public methods ******************************

bool MetricTaubNUT::calculateMetric(const double* pos)
{
    // double r     = pos[1];
    double theta = pos[2];
    double l = mL;

    double D, S;
    calcFunctions(pos, D, S);

    g_compts[0][0] = -D / S;
    g_compts[0][1] = 0;
    g_compts[0][2] = 0;
    g_compts[0][3] = -2 * l * D * cos(theta) / S;
    g_compts[1][0] = 0;
    g_compts[1][1] = S / D;
    g_compts[1][2] = 0;
    g_compts[1][3] = 0;
    g_compts[2][0] = 0;
    g_compts[2][1] = 0;
    g_compts[2][2] = S;
    g_compts[2][3] = 0;
    g_compts[3][0] = -2 * l * D * cos(theta) / S;
    g_compts[3][1] = 0;
    g_compts[3][2] = 0;
    g_compts[3][3] = -4 * pow(l, 2) * D * pow(cos(theta), 2) / S + S * pow(sin(theta), 2);

    return true;
}

bool MetricTaubNUT::calculateChristoffels(const double* pos)
{
    double r = pos[1];
    double theta = pos[2];
    double M = mMass;
    double l = mL;

    christoffel[0][0][0] = 0;
    christoffel[0][0][1] = (-2 * r * (-2 * M * r - pow(l, 2) + pow(r, 2)) + (-2 * M + 2 * r) * (pow(l, 2) + pow(r, 2)))
        * (-2 * M * r - pow(l, 2) + pow(r, 2)) / (2 * pow(pow(l, 2) + pow(r, 2), 3));
    christoffel[0][0][2] = 0;
    christoffel[0][0][3] = 0;
    christoffel[0][1][0] = (-2 * r * (-2 * M * r - pow(l, 2) + pow(r, 2)) + (-2 * M + 2 * r) * (pow(l, 2) + pow(r, 2)))
        / (2 * (pow(l, 2) + pow(r, 2)) * (-2 * M * r - pow(l, 2) + pow(r, 2)));
    christoffel[0][1][1] = 0;
    christoffel[0][1][2] = 0;
    christoffel[0][1][3] = 0;
    christoffel[0][2][0]
        = -2 * pow(l, 2) * (-2 * M * r - pow(l, 2) + pow(r, 2)) / (pow(pow(l, 2) + pow(r, 2), 2) * tan(theta));
    christoffel[0][2][1] = 0;
    christoffel[0][2][2] = 0;
    christoffel[0][2][3] = l * (-2 * M * r - pow(l, 2) + pow(r, 2)) / (pow(pow(l, 2) + pow(r, 2), 2) * sin(theta));
    christoffel[0][3][0] = 0;
    christoffel[0][3][1] = l
        * (-2 * r * (-2 * M * r - pow(l, 2) + pow(r, 2)) + (-2 * M + 2 * r) * (pow(l, 2) + pow(r, 2)))
        * (-2 * M * r - pow(l, 2) + pow(r, 2)) * cos(theta) / pow(pow(l, 2) + pow(r, 2), 3);
    christoffel[0][3][2] = -l * (-2 * M * r - pow(l, 2) + pow(r, 2)) * sin(theta) / pow(pow(l, 2) + pow(r, 2), 2);
    christoffel[0][3][3] = 0;
    christoffel[1][0][0] = (-2 * r * (-2 * M * r - pow(l, 2) + pow(r, 2)) + (-2 * M + 2 * r) * (pow(l, 2) + pow(r, 2)))
        / (2 * (pow(l, 2) + pow(r, 2)) * (-2 * M * r - pow(l, 2) + pow(r, 2)));
    christoffel[1][0][1] = 0;
    christoffel[1][0][2] = 0;
    christoffel[1][0][3] = 0;
    christoffel[1][1][0] = 0;
    christoffel[1][1][1] = (2 * r * (-2 * M * r - pow(l, 2) + pow(r, 2)) - (-2 * M + 2 * r) * (pow(l, 2) + pow(r, 2)))
        / (2 * (pow(l, 2) + pow(r, 2)) * (-2 * M * r - pow(l, 2) + pow(r, 2)));
    christoffel[1][1][2] = 0;
    christoffel[1][1][3] = 0;
    christoffel[1][2][0] = 0;
    christoffel[1][2][1] = 0;
    christoffel[1][2][2] = r / (pow(l, 2) + pow(r, 2));
    christoffel[1][2][3] = 0;
    christoffel[1][3][0] = l
        * (-4 * r * (-2 * M * r - pow(l, 2) + pow(r, 2)) + (-2 * M + 2 * r) * (pow(l, 2) + pow(r, 2))) * cos(theta)
        / ((pow(l, 2) + pow(r, 2)) * (-2 * M * r - pow(l, 2) + pow(r, 2)));
    christoffel[1][3][1] = 0;
    christoffel[1][3][2] = 0;
    christoffel[1][3][3] = r / (pow(l, 2) + pow(r, 2));
    christoffel[2][0][0]
        = -2 * pow(l, 2) * (-2 * M * r - pow(l, 2) + pow(r, 2)) / (pow(pow(l, 2) + pow(r, 2), 2) * tan(theta));
    christoffel[2][0][1] = 0;
    christoffel[2][0][2] = 0;
    christoffel[2][0][3] = l * (-2 * M * r - pow(l, 2) + pow(r, 2)) / (pow(pow(l, 2) + pow(r, 2), 2) * sin(theta));
    christoffel[2][1][0] = 0;
    christoffel[2][1][1] = 0;
    christoffel[2][1][2] = r / (pow(l, 2) + pow(r, 2));
    christoffel[2][1][3] = 0;
    christoffel[2][2][0] = 0;
    christoffel[2][2][1] = -r * (-2 * M * r - pow(l, 2) + pow(r, 2)) / (pow(l, 2) + pow(r, 2));
    christoffel[2][2][2] = 0;
    christoffel[2][2][3] = 0;
    christoffel[2][3][0] = l
        * (-4 * pow(l, 2) * (-2 * M * r - pow(l, 2) + pow(r, 2)) * pow(cos(theta), 2)
            + pow(pow(l, 2) + pow(r, 2), 2) * pow(sin(theta), 2) - 2 * pow(pow(l, 2) + pow(r, 2), 2))
        / (pow(pow(l, 2) + pow(r, 2), 2) * sin(theta));
    christoffel[2][3][1] = 0;
    christoffel[2][3][2] = 0;
    christoffel[2][3][3] = (2 * pow(l, 2) * (-2 * M * r - pow(l, 2) + pow(r, 2)) + pow(pow(l, 2) + pow(r, 2), 2))
        / (pow(pow(l, 2) + pow(r, 2), 2) * tan(theta));
    christoffel[3][0][0] = 0;
    christoffel[3][0][1] = l
        * (-2 * r * (-2 * M * r - pow(l, 2) + pow(r, 2)) + (-2 * M + 2 * r) * (pow(l, 2) + pow(r, 2)))
        * (-2 * M * r - pow(l, 2) + pow(r, 2)) * cos(theta) / pow(pow(l, 2) + pow(r, 2), 3);
    christoffel[3][0][2] = -l * (-2 * M * r - pow(l, 2) + pow(r, 2)) * sin(theta) / pow(pow(l, 2) + pow(r, 2), 2);
    christoffel[3][0][3] = 0;
    christoffel[3][1][0] = l
        * (-4 * r * (-2 * M * r - pow(l, 2) + pow(r, 2)) + (-2 * M + 2 * r) * (pow(l, 2) + pow(r, 2))) * cos(theta)
        / ((pow(l, 2) + pow(r, 2)) * (-2 * M * r - pow(l, 2) + pow(r, 2)));
    christoffel[3][1][1] = 0;
    christoffel[3][1][2] = 0;
    christoffel[3][1][3] = r / (pow(l, 2) + pow(r, 2));
    christoffel[3][2][0] = l
        * (-4 * pow(l, 2) * (-2 * M * r - pow(l, 2) + pow(r, 2)) * pow(cos(theta), 2)
            + pow(pow(l, 2) + pow(r, 2), 2) * pow(sin(theta), 2) - 2 * pow(pow(l, 2) + pow(r, 2), 2))
        / (pow(pow(l, 2) + pow(r, 2), 2) * sin(theta));
    christoffel[3][2][1] = 0;
    christoffel[3][2][2] = 0;
    christoffel[3][2][3] = (2 * pow(l, 2) * (-2 * M * r - pow(l, 2) + pow(r, 2)) + pow(pow(l, 2) + pow(r, 2), 2))
        / (pow(pow(l, 2) + pow(r, 2), 2) * tan(theta));
    christoffel[3][3][0] = 0;
    christoffel[3][3][1] = (-2 * M * r - pow(l, 2) + pow(r, 2))
        * (-8 * pow(l, 2) * r * (-2 * M * r - pow(l, 2) + pow(r, 2)) * pow(cos(theta), 2)
            + 4 * pow(l, 2) * (-2 * M + 2 * r) * (pow(l, 2) + pow(r, 2)) * pow(cos(theta), 2)
            + 2 * r * pow(pow(l, 2) + pow(r, 2), 2) * pow(cos(theta), 2) - 2 * r * pow(pow(l, 2) + pow(r, 2), 2))
        / (2 * pow(pow(l, 2) + pow(r, 2), 3));
    christoffel[3][3][2] = (-4 * pow(l, 2) * (-2 * M * r - pow(l, 2) + pow(r, 2)) - pow(pow(l, 2) + pow(r, 2), 2))
        * sin(theta) * cos(theta) / pow(pow(l, 2) + pow(r, 2), 2);
    christoffel[3][3][3] = 0;

    return true;
}

bool MetricTaubNUT::calculateChrisD(const double* pos)
{
    double r = pos[1];
    double theta = pos[2];
    double M = mMass;
    double l = mL;

    double t1 = r * r;
    double t2 = t1 * t1;
    double t3 = l * l;
    double t4 = t2 * t3;
    double t5 = 3.0 * t4;
    double t6 = t3 * t3;
    double t7 = t1 * t6;
    double t9 = t2 * r;
    double t10 = M * t9;
    double t11 = t1 * r;
    double t12 = M * t11;
    double t13 = t12 * t3;
    double t15 = t6 * M;
    double t16 = t15 * r;
    double t18 = M * M;
    double t19 = t18 * t2;
    double t21 = t18 * t6;
    double t22 = t18 * t1;
    double t23 = t22 * t3;
    double t25 = t6 * t3;
    double t27 = t1 + t3;
    double t28 = t27 * t27;
    double t29 = t28 * t28;
    double t30 = 1 / t29;
    double t31 = (t5 - 8.0 * t7 + t10 - 14.0 * t13 + 9.0 * t16 - 3.0 * t19 - t21 + 8.0 * t23 + t25) * t30;
    double t36 = 1 / t28;
    double t38 = M * r;
    double t40 = t1 - 2.0 * t38 - t3;
    double t41 = t40 * t40;
    double t42 = 1 / t41;
    double t44 = 2.0 * (t5 + t10 - 6.0 * t13 + t16 - t19 + 4.0 * t23 + t25 + t21) * t36 * t42;
    double t45 = cos(theta);
    double t47 = M * t1;
    double t49 = r * t3;
    double t51 = M * t3;
    double t52 = t11 - 3.0 * t47 - 3.0 * t49 + t51;
    double t53 = sin(theta);
    double t54 = 1 / t53;
    double t57 = 1 / t28 / t27;
    double t58 = t52 * t54 * t57;
    double t60 = 4.0 * t45 * t3 * t58;
    double t62 = t53 * t53;
    double t63 = t45 * t45;
    double t64 = t62 + t63;
    double t65 = 1 / t62;
    double t69 = 2.0 * t3 * t40 * t64 * t65 * t36;
    double t73 = 2.0 * l * t52 * t54 * t57;
    double t74 = t65 * t36;
    double t77 = t74 * t40 * l * t45;
    double t78 = t45 * l;
    double t80 = 4.0 * t78 * t31;
    double t82 = t53 * l;
    double t87 = 2.0 * t40 * t57 * t82 * (2.0 * t49 + t47 - t51);
    double t90 = 2.0 * t82 * t52 * t57;
    double t92 = t36 * t40 * t78;
    double t94 = (t1 - t3) * t36;
    double t95 = t2 * t1;
    double t109 = 2.0 * t78
        * (t95 - 6.0 * t10 - 9.0 * t4 + 6.0 * t19 + 3.0 * t7 - 12.0 * t23 - 6.0 * t16 - 3.0 * t25 + 20.0 * t13
            - 2.0 * t21)
        * t36 * t42;
    double t116 = 2.0 * t52 * t53 * l / t27 / t40;
    double t117 = t1 * t3;
    double t119 = t51 * r;
    double t121 = 4.0 * t117 - 4.0 * t119 - t6 + t2;
    double t126 = 8.0 * t3 * l * t63 * t58;
    double t130 = t62 * t3;
    double t133 = t62 * t6;
    double t137 = t3 * t63;
    double t142 = t6 * t63;
    double t149 = t78
        * (14.0 * t62 * t1 * t3 - 16.0 * t130 * t38 - 5.0 * t133 + 3.0 * t62 * t2 + 8.0 * t137 * t1 - 8.0 * t137 * t38
            - 2.0 * t142 + 2.0 * t63 * t2)
        * t65 * t36;
    double t151 = t121 * t64 * t74;
    double t161 = t63 * M;
    double t172 = t2 * t6;
    double t177 = t1 * t25;
    double t183 = t2 * t2;
    double t185 = t6 * t6;
    double t196 = -4.0 * t10 * t130 - 8.0 * t12 * t133 - 4.0 * M * t25 * t62 * r + 8.0 * t9 * t3 * t161
        - 112.0 * t11 * t6 * t161 + 72.0 * r * t25 * t161 + 64.0 * t22 * t142 + 8.0 * t172 * t62 + 24.0 * t172 * t63
        - 64.0 * t177 * t63 + 6.0 * t95 * t3 * t62 + t62 * t183 + 8.0 * t185 * t63 - t62 * t185 + 2.0 * t177 * t62
        - 8.0 * t25 * t63 * t18 - 24.0 * t19 * t137;

    chrisD[0][0][0][0] = 0.0;
    chrisD[0][0][0][1] = 0.0;
    chrisD[0][0][0][2] = 0.0;
    chrisD[0][0][0][3] = 0.0;
    chrisD[0][0][1][0] = 0.0;
    chrisD[0][0][1][1] = -2.0 * t31;
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
    chrisD[0][1][0][1] = -t44;
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
    chrisD[0][2][0][1] = t60;
    chrisD[0][2][0][2] = t69;
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
    chrisD[0][2][3][1] = -t73;
    chrisD[0][2][3][2] = -t77;
    chrisD[0][2][3][3] = 0.0;
    chrisD[0][3][0][0] = 0.0;
    chrisD[0][3][0][1] = 0.0;
    chrisD[0][3][0][2] = 0.0;
    chrisD[0][3][0][3] = 0.0;
    chrisD[0][3][1][0] = 0.0;
    chrisD[0][3][1][1] = -t80;
    chrisD[0][3][1][2] = -t87;
    chrisD[0][3][1][3] = 0.0;
    chrisD[0][3][2][0] = 0.0;
    chrisD[0][3][2][1] = t90;
    chrisD[0][3][2][2] = -t92;
    chrisD[0][3][2][3] = 0.0;
    chrisD[0][3][3][0] = 0.0;
    chrisD[0][3][3][1] = 0.0;
    chrisD[0][3][3][2] = 0.0;
    chrisD[0][3][3][3] = 0.0;
    chrisD[1][0][0][0] = 0.0;
    chrisD[1][0][0][1] = -t44;
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
    chrisD[1][2][2][1] = -t94;
    chrisD[1][2][2][2] = 0.0;
    chrisD[1][2][2][3] = 0.0;
    chrisD[1][2][3][0] = 0.0;
    chrisD[1][2][3][1] = 0.0;
    chrisD[1][2][3][2] = 0.0;
    chrisD[1][2][3][3] = 0.0;
    chrisD[1][3][0][0] = 0.0;
    chrisD[1][3][0][1] = t109;
    chrisD[1][3][0][2] = t116;
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
    chrisD[1][3][3][1] = -t94;
    chrisD[1][3][3][2] = 0.0;
    chrisD[1][3][3][3] = 0.0;
    chrisD[2][0][0][0] = 0.0;
    chrisD[2][0][0][1] = t60;
    chrisD[2][0][0][2] = t69;
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
    chrisD[2][0][3][1] = -t73;
    chrisD[2][0][3][2] = -t77;
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
    chrisD[2][1][2][1] = -t94;
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
    chrisD[2][2][1][1] = -t121 * t36;
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
    chrisD[2][3][0][1] = t126;
    chrisD[2][3][0][2] = t149;
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
    chrisD[2][3][3][1] = -t60;
    chrisD[2][3][3][2] = -t151;
    chrisD[2][3][3][3] = 0.0;
    chrisD[3][0][0][0] = 0.0;
    chrisD[3][0][0][1] = 0.0;
    chrisD[3][0][0][2] = 0.0;
    chrisD[3][0][0][3] = 0.0;
    chrisD[3][0][1][0] = 0.0;
    chrisD[3][0][1][1] = -t80;
    chrisD[3][0][1][2] = -t87;
    chrisD[3][0][1][3] = 0.0;
    chrisD[3][0][2][0] = 0.0;
    chrisD[3][0][2][1] = t90;
    chrisD[3][0][2][2] = -t92;
    chrisD[3][0][2][3] = 0.0;
    chrisD[3][0][3][0] = 0.0;
    chrisD[3][0][3][1] = 0.0;
    chrisD[3][0][3][2] = 0.0;
    chrisD[3][0][3][3] = 0.0;
    chrisD[3][1][0][0] = 0.0;
    chrisD[3][1][0][1] = t109;
    chrisD[3][1][0][2] = t116;
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
    chrisD[3][1][3][1] = -t94;
    chrisD[3][1][3][2] = 0.0;
    chrisD[3][1][3][3] = 0.0;
    chrisD[3][2][0][0] = 0.0;
    chrisD[3][2][0][1] = t126;
    chrisD[3][2][0][2] = t149;
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
    chrisD[3][2][3][1] = -t60;
    chrisD[3][2][3][2] = -t151;
    chrisD[3][2][3][3] = 0.0;
    chrisD[3][3][0][0] = 0.0;
    chrisD[3][3][0][1] = 0.0;
    chrisD[3][3][0][2] = 0.0;
    chrisD[3][3][0][3] = 0.0;
    chrisD[3][3][1][0] = 0.0;
    chrisD[3][3][1][1] = -t196 * t30;
    chrisD[3][3][1][2]
        = -2.0 * t40 * t45 * t53 * (9.0 * r * t6 + 4.0 * t47 * t3 - 4.0 * t15 + t9 + 2.0 * t11 * t3) * t57;
    chrisD[3][3][1][3] = 0.0;
    chrisD[3][3][2][0] = 0.0;
    chrisD[3][3][2][1] = 8.0 * t45 * t53 * t3 * t52 * t57;
    chrisD[3][3][2][2] = -(6.0 * t117 - 8.0 * t119 - 3.0 * t6 + t2) * (-t62 + t63) * t36;
    chrisD[3][3][2][3] = 0.0;
    chrisD[3][3][3][0] = 0.0;
    chrisD[3][3][3][1] = 0.0;
    chrisD[3][3][3][2] = 0.0;
    chrisD[3][3][3][3] = 0.0;

    return true;
}

void MetricTaubNUT::localToCoord(const double* pos, const double* ldir, double* dir, enum_nat_tetrad_type)
{
    double r = pos[1];
    double theta = pos[2];
    double sigma = r * r + mL * mL;
    double delta = r * r - 2.0 * mMass * r - mL * mL;

    dir[0] = sqrt(sigma / delta) * ldir[0] - 2.0 * mL * cos(theta) / sin(theta) / sqrt(sigma) * ldir[3];
    dir[1] = sqrt(delta / sigma) * ldir[1];
    dir[2] = ldir[2] / sqrt(sigma);
    dir[3] = ldir[3] / (sqrt(sigma) * sin(theta));
}

void MetricTaubNUT::coordToLocal(const double* pos, const double* cdir, double* ldir, enum_nat_tetrad_type)
{
    double r = pos[1];
    double theta = pos[2];
    double sigma = r * r + mL * mL;
    double delta = r * r - 2.0 * mMass * r - mL * mL;

    ldir[0] = sqrt(delta / sigma) * cdir[0] + 2.0 * mL * cos(theta) * sqrt(delta / sigma) * cdir[3];
    ldir[1] = sqrt(sigma / delta) * cdir[1];
    ldir[2] = sqrt(sigma) * cdir[2];
    ldir[3] = sqrt(sigma) * sin(theta) * cdir[3];
}

bool MetricTaubNUT::breakCondition(const double* pos)
{
    bool br = false;

    double r = pos[1];
    double delta = r * r - 2.0 * mMass * r - mL * mL;
    if (fabs(delta) < 1e-12) {
        br = true;
    }
    return br;
}

bool MetricTaubNUT::calcDerivs(const double y[], double dydx[])
{
    return false;
    dydx[0] = y[4];
    dydx[1] = y[5];
    dydx[2] = y[6];
    dydx[3] = y[7];

    double r = y[1];
    double theta = y[2];

    double st = sin(theta);
    double ct = cos(theta);
    double Sigma = r * r + mL * mL;
    double S2 = Sigma * Sigma;
    double S3 = S2 * Sigma;
    double Delta = r * r - 2.0 * mMass * r - mL * mL;
    double rho = 2.0 * r * mL * mL + mMass * (r * r - mL * mL);

    double G_t_t_r = Delta * rho / S3;
    double G_t_r_t = rho / Delta / Sigma;
    double G_t_th_t = -2.0 * mL * mL * ct / st * Delta / S2;
    double G_th_ph_ph
        = (r * r * r * r - mL * mL * mL * mL + 4.0 * r * r * mL * mL - 4.0 * mMass * r * mL * mL) * ct / st / S2;
    double G_ph_ph_th
        = -(6.0 * r * r * mL * mL - 8.0 * mMass * r * mL * mL - 3.0 * mL * mL * mL * mL + r * r * r * r) * st * ct / S2;
    double G_t_th_ph = mL * Delta / S2 / st;
    double G_t_ph_r = 2.0 * mL * rho * Delta * ct / S3;
    double G_t_ph_th = -mL * Delta * st / S2;
    double G_r_r_r = -rho / Sigma / Delta;
    double G_r_th_th = r / Sigma;
    double G_r_ph_t = -2.0 * mL * (r * r * (r - 3.0 * mMass) - (3.0 * r - mMass) * mL * mL) * ct / Sigma / Delta;
    double G_r_ph_ph = r / Sigma;
    double G_th_th_r = -r * Delta / Sigma;
    double G_th_ph_t = -mL
        * (ct * ct * (6.0 * r * r * mL * mL - 8.0 * mL * mL * mMass * r - 3.0 * mL * mL * mL * mL + r * r * r * r) + S2)
        / S2 / st;
    double G_ph_ph_r = Delta / S3
        * (ct * ct
                * (9.0 * r * mL * mL * mL * mL + 4.0 * mL * mL * mMass * r * r - 4.0 * mL * mL * mL * mL * mMass
                    + r * r * r * (r * r + 2.0 * mL * mL))
            - r * S2);

    dydx[4] = -2.0 * G_t_r_t * y[4] * y[5] - 2.0 * G_r_ph_t * y[5] * y[7] - 2.0 * G_th_ph_t * y[6] * y[7]
        - 2.0 * G_t_th_t * y[4] * y[6];
    dydx[5] = -G_t_t_r * y[4] * y[4] - 2.0 * G_t_ph_r * y[4] * y[7] - G_r_r_r * y[5] * y[5] - G_th_th_r * y[6] * y[6]
        - G_ph_ph_r * y[7] * y[7];
    dydx[6] = -2.0 * G_t_ph_th * y[4] * y[7] - 2.0 * G_r_th_th * y[5] * y[6] - G_ph_ph_th * y[7] * y[7];
    dydx[7] = -2.0 * G_t_th_ph * y[4] * y[6] - 2.0 * G_r_ph_ph * y[5] * y[7] - 2.0 * G_th_ph_ph * y[6] * y[7];
    return true;
}

double MetricTaubNUT::testConstraint(const double y[], const double kappa)
{
    double r = y[1];
    double theta = y[2];

    double dt = y[4];
    double dr = y[5];
    double dth = y[6];
    double dph = y[7];

    double delta = r * r - 2.0 * mMass * r - mL * mL;
    double sigma = r * r + mL * mL;
    double kl = dt + 2.0 * mL * cos(theta) * dph;
    double sum = -kappa;
    sum += -delta / sigma * kl * kl + sigma * (dr * dr / delta + dth * dth + sin(theta) * sin(theta) * dph * dph);
    return sum;
}

bool MetricTaubNUT::setParam(const char* pName, double val)
{
    Metric::setParam(pName, val);

    if (strcmp(pName, "mass") == 0) {
        mMass = val;
    }
    else if (strcmp(pName, "l") == 0) {
        mL = val;
    }
    return true;
}

bool MetricTaubNUT::effPotentialValue(
    const vec4 pos, const vec4 cdir, enum_geodesic_type type, const double x, double& val)
{
    double kappa = 0.0;
    if (type == enum_geodesic_timelike) {
        kappa = -mSign;
    }

    double delta = x * x - 2.0 * mMass * x - mL * mL;
    if (delta <= 0.0 || x < 0.0) {
        return false;
    }

    double h = (pos[1] * pos[1] + mL * mL) * cdir[3];
    double sigma = x * x + mL * mL;
    val = 0.5 * delta / sigma * (h * h / sigma - kappa);
    return true;
}

bool MetricTaubNUT::totEnergy(const vec4 pos, const vec4 cdir, const double, double& val)
{
    double sigma = (pos[1] * pos[1] + mL * mL);
    double delta = (pos[1] * pos[1] - 2.0 * mMass * pos[1] - mL * mL);

    // 1/2*k^2/c^2:
    double k = delta / sigma * cdir[0];
    val = 0.5 * k * k;
    return true;
}

bool MetricTaubNUT::report(const vec4, const vec4, char*& text)
{
    std::stringstream ss;
    ss << "Report for TaubNUT metric\n\tcoordinate : (t,r,theta,phi)\n";
    ss << "---------------------------------------------------------\n";
    ss << "  physical units ................................. no\n";

    calcCriticalRadius();
    calcPhotonOrbit();
    ss.precision(DEF_FIXED_REPORT_PRECISION);
    ss << "  critical avalue ....................... r_crit = " << mCritRadius << std::endl;
    ss << "  photon orbit ............................ r_po = " << mPhotonOrbit << std::endl;

    text = new char[ss.str().length() + 2];
    return CopyString(ss.str().c_str(), text);
}

// ********************************* protected methods *****************************

void MetricTaubNUT::setStandardValues()
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

void MetricTaubNUT::calcCriticalRadius()
{
    mCritRadius = mMass + sqrt(mMass * mMass + mL * mL);
}

void MetricTaubNUT::calcPhotonOrbit()
{
    double w = sqrt(mMass * mMass + mL * mL);
    double psi = acos(mMass / w);

    mPhotonOrbit = 2.0 * w * cos(psi / 3.0) + mMass;
}

void MetricTaubNUT::calcFunctions(const double* pos, double& D, double& S)
{
    double r = pos[1];
    D = r * r - 2.0 * mMass * r - mL * mL;
    S = r * r + mL * mL;
}

} // end namespace m4d
