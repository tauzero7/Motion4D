/**
 * @file    m4dMetricCurzon.cpp
 * @author  Thomas Mueller
 *
 * This file is part of the m4d-library.
 */
#include "m4dMetricCurzon.h"

namespace m4d {

MetricCurzon::MetricCurzon(double mass)
{
    mMetricName = "Curzon";
    setCoordType(enum_coordinate_cylinder);

    mPhysicalUnits = enum_physical_constants_geom;
    mSpeedOfLight = 1.0;
    mGravConstant = 1.0;
    mDielectricPerm = 1.0;

    addParam("mass", mass);
    mMass = mass;

    mSign = 1.0;
    mLocTeds.push_back(enum_nat_tetrad_static);

    setStandardValues();
}

MetricCurzon::~MetricCurzon() {}

// *********************************** public methods ******************************

bool MetricCurzon::calculateMetric(const double* pos)
{
    double r = pos[1];
    double lambda, nu;
    calcLambdaNu(pos, lambda, nu);

    double m = mMass;
    double z = pos[3];

    double t1 = r * r;
    double t2 = z * z;
    double t3 = t1 + t2;
    double t4 = sqrt(t3);
    double t7 = exp(m / t4);
    double t8 = t7 * t7;
    double t10 = m * m;
    double t12 = t3 * t3;
    double t15 = exp(t10 * t1 / t12);
    double t17 = 1 / t15 * t8;

    g_compts[0][0] = -1 / t8;
    g_compts[0][1] = 0.0;
    g_compts[0][2] = 0.0;
    g_compts[0][3] = 0.0;
    g_compts[1][0] = 0.0;
    g_compts[1][1] = t17;
    g_compts[1][2] = 0.0;
    g_compts[1][3] = 0.0;
    g_compts[2][0] = 0.0;
    g_compts[2][1] = 0.0;
    g_compts[2][2] = t1 * t8;
    g_compts[2][3] = 0.0;
    g_compts[3][0] = 0.0;
    g_compts[3][1] = 0.0;
    g_compts[3][2] = 0.0;
    g_compts[3][3] = t17;

    return true;
}

bool MetricCurzon::calculateChristoffels(const double* pos)
{
    double r = pos[1];
    double lambda, nu, dldr, dldz, dndr, dndz;
    calcDLN(pos, lambda, nu, dldr, dldz, dndr, dndz);

    double m = mMass;
    double z = pos[3];

    double t1 = r * r;
    double t2 = z * z;
    double t3 = t1 + t2;
    double t4 = sqrt(t3);
    double t7 = exp(m / t4);
    double t8 = t7 * t7;
    double t9 = t8 * t8;
    double t11 = m * m;
    double t13 = t3 * t3;
    double t16 = exp(t11 * t1 / t13);
    double t17 = 1 / t9 * t16;
    double t18 = t4 * t3;
    double t19 = 1 / t18;
    double t20 = m * t19;
    double t21 = t20 * r;
    double t23 = t20 * z;
    double t26 = m * t18;
    double t27 = t26 * t1;
    double t29 = t1 * t1;
    double t30 = t29 * t1;
    double t32 = 3.0 * t29 * t2;
    double t33 = t2 * t2;
    double t35 = 3.0 * t1 * t33;
    double t36 = t33 * t2;
    double t38 = t13 * t13;
    double t40 = 1 / t4 / t38;
    double t42 = m * r * (-t27 + t26 * t2 + t30 + t32 + t35 + t36) * t40;
    double t47 = m * z * (2.0 * t27 - t30 - t32 - t35 - t36) * t40;
    double t50 = t18 - m * t1;
    double t52 = 1 / r * t50 * t19;

    christoffel[0][0][0] = 0.0;
    christoffel[0][0][1] = t17 * t21;
    christoffel[0][0][2] = 0.0;
    christoffel[0][0][3] = t17 * t23;
    christoffel[0][1][0] = t21;
    christoffel[0][1][1] = 0.0;
    christoffel[0][1][2] = 0.0;
    christoffel[0][1][3] = 0.0;
    christoffel[0][2][0] = 0.0;
    christoffel[0][2][1] = 0.0;
    christoffel[0][2][2] = 0.0;
    christoffel[0][2][3] = 0.0;
    christoffel[0][3][0] = t23;
    christoffel[0][3][1] = 0.0;
    christoffel[0][3][2] = 0.0;
    christoffel[0][3][3] = 0.0;
    christoffel[1][0][0] = t21;
    christoffel[1][0][1] = 0.0;
    christoffel[1][0][2] = 0.0;
    christoffel[1][0][3] = 0.0;
    christoffel[1][1][0] = 0.0;
    christoffel[1][1][1] = -t42;
    christoffel[1][1][2] = 0.0;
    christoffel[1][1][3] = -t47;
    christoffel[1][2][0] = 0.0;
    christoffel[1][2][1] = 0.0;
    christoffel[1][2][2] = t52;
    christoffel[1][2][3] = 0.0;
    christoffel[1][3][0] = 0.0;
    christoffel[1][3][1] = t47;
    christoffel[1][3][2] = 0.0;
    christoffel[1][3][3] = -t42;
    christoffel[2][0][0] = 0.0;
    christoffel[2][0][1] = 0.0;
    christoffel[2][0][2] = 0.0;
    christoffel[2][0][3] = 0.0;
    christoffel[2][1][0] = 0.0;
    christoffel[2][1][1] = 0.0;
    christoffel[2][1][2] = t52;
    christoffel[2][1][3] = 0.0;
    christoffel[2][2][0] = 0.0;
    christoffel[2][2][1] = -t16 * r * t50 * t19;
    christoffel[2][2][2] = 0.0;
    christoffel[2][2][3] = t16 * t1 * t23;
    christoffel[2][3][0] = 0.0;
    christoffel[2][3][1] = 0.0;
    christoffel[2][3][2] = -t23;
    christoffel[2][3][3] = 0.0;
    christoffel[3][0][0] = t23;
    christoffel[3][0][1] = 0.0;
    christoffel[3][0][2] = 0.0;
    christoffel[3][0][3] = 0.0;
    christoffel[3][1][0] = 0.0;
    christoffel[3][1][1] = t47;
    christoffel[3][1][2] = 0.0;
    christoffel[3][1][3] = -t42;
    christoffel[3][2][0] = 0.0;
    christoffel[3][2][1] = 0.0;
    christoffel[3][2][2] = -t23;
    christoffel[3][2][3] = 0.0;
    christoffel[3][3][0] = 0.0;
    christoffel[3][3][1] = t42;
    christoffel[3][3][2] = 0.0;
    christoffel[3][3][3] = t47;

    return true;
}

bool MetricCurzon::calculateChrisD(const double* pos)
{
    double r = pos[1];
    double m = mMass;
    double z = pos[3];

    double t1 = m * m;
    double t2 = r * r;
    double t3 = t1 * t2;
    double t4 = z * z;
    double t5 = t2 + t4;
    double t6 = t5 * t5;
    double t9 = exp(t3 / t6);
    double t10 = t9 * m;
    double t11 = m * t2;
    double t12 = sqrt(t5);
    double t13 = t12 * t6;
    double t14 = t11 * t13;
    double t16 = t2 * t2;
    double t17 = t16 * t2;
    double t18 = t1 * t17;
    double t19 = 2.0 * t18;
    double t20 = t4 * t4;
    double t22 = 2.0 * t3 * t20;
    double t23 = t16 * t16;
    double t24 = 2.0 * t23;
    double t25 = t17 * t4;
    double t26 = 5.0 * t25;
    double t28 = 3.0 * t16 * t20;
    double t29 = t20 * t4;
    double t30 = t29 * t2;
    double t31 = t20 * t20;
    double t35 = exp(m / t12);
    double t36 = t35 * t35;
    double t37 = t36 * t36;
    double t38 = 1 / t37;
    double t40 = t6 * t6;
    double t43 = 1 / t12 / t40 / t5;
    double t46 = t10 * r;
    double t47 = t12 * t40;
    double t50 = t1 * t23;
    double t54 = t16 * t1;
    double t59 = t23 * t2;
    double t69 = t31 * t4;
    double t71 = 4.0 * m * t47 - 4.0 * t50 - 12.0 * t18 * t4 - 12.0 * t54 * t20 - 4.0 * t3 * t29 - 3.0 * t59
        - 15.0 * t23 * t4 - 30.0 * t17 * t20 - 30.0 * t16 * t29 - 15.0 * t2 * t31 - 3.0 * t69;
    double t73 = t6 * t5;
    double t76 = 1 / t12 / t40 / t73;
    double t80 = m * t13;
    double t81 = 4.0 * t80;
    double t85 = 3.0 * t17;
    double t87 = 9.0 * t16 * t4;
    double t88 = t2 * t20;
    double t89 = 9.0 * t88;
    double t90 = 3.0 * t29;
    double t107 = m * t4;
    double t118 = -2.0 * t31 * t20 - 12.0 * t54 * t29 + t23 * t16 - 12.0 * t18 * t20 - 4.0 * t3 * t31 - 4.0 * t50 * t4
        + 4.0 * t107 * t47 + 3.0 * t59 * t4 - 10.0 * t17 * t29 - 15.0 * t16 * t31 - 9.0 * t2 * t69;
    double t125 = 1 / t13;
    double t126 = m * (2.0 * t2 - t4) * t125;
    double t130 = 3.0 * m * t125 * r * z;
    double t134 = m * (-2.0 * t4 + t2) * t125;
    double t137 = t12 * t5;
    double t141 = t2 * t137;
    double t142 = t141 * t107;
    double t146 = m * (3.0 * t14 - t107 * t13 + t24 + t26 + t28 - t30 - t31 - 6.0 * t16 * t137 * m + 6.0 * t142) * t43;
    double t147 = m * r;
    double t148 = m * t137;
    double t149 = t148 * t2;
    double t157 = t147 * z * (6.0 * t149 - 6.0 * t148 * t4 + 2.0 * t80 - t85 - t87 - t89 - t90) * t43;
    double t162 = t147 * z * (-12.0 * t149 + t81 + t85 + t87 + t89 + t90) * t43;
    double t164 = 5.0 * t30;
    double t165 = 2.0 * t31;
    double t169 = m * (2.0 * t14 - t23 - t25 + t28 + t164 + t165 - 12.0 * t142) * t43;
    double t170 = t16 * m;
    double t172 = t11 * t4;
    double t176 = (t13 - 2.0 * t170 + t172) / t2 * t125;
    double t177 = t137 * t4;
    double t181 = t12 * t73;
    double t185 = 3.0 * z * (t13 - t141 - t177 + t170 + t172) / r / t181;
    double t192 = t23 * m;
    double t195 = t1 * m;
    double t205 = t195 * t16;
    double t208 = t17 * m;
    double t215 = 3.0 * t13 * t16 * t4 + 3.0 * t13 * t2 * t20 - 3.0 * t192 * t4 + 2.0 * t195 * t23 + t13 * t17
        + t13 * t29 - 2.0 * t54 * t13 + 2.0 * t3 * t13 * t4 - 2.0 * t205 * t20 - 9.0 * t208 * t20 - 9.0 * t170 * t29
        - 3.0 * t11 * t31;
    double t239 = -4.0 * t54 * t137 - 4.0 * t3 * t177 + 4.0 * t195 * t17 + 4.0 * t205 * t4 + 3.0 * t47 - 3.0 * t181 * t2
        - 3.0 * t181 * t4 + 3.0 * t192 + 9.0 * t208 * t4 + 9.0 * t170 * t20 + 3.0 * t11 * t29;

    chrisD[0][0][0][0] = 0.0;
    chrisD[0][0][0][1] = 0.0;
    chrisD[0][0][0][2] = 0.0;
    chrisD[0][0][0][3] = 0.0;
    chrisD[0][0][1][0] = 0.0;
    chrisD[0][0][1][1] = t10 * (4.0 * t14 - t19 + t22 - t24 - t26 - t28 + t30 + t31) * t38 * t43;
    chrisD[0][0][1][2] = 0.0;
    chrisD[0][0][1][3] = t46 * z * t71 * t38 * t76;
    chrisD[0][0][2][0] = 0.0;
    chrisD[0][0][2][1] = 0.0;
    chrisD[0][0][2][2] = 0.0;
    chrisD[0][0][2][3] = 0.0;
    chrisD[0][0][3][0] = 0.0;
    chrisD[0][0][3][1] = t46 * z * (t81 - 2.0 * t54 + 2.0 * t1 * t20 - t85 - t87 - t89 - t90) * t38 * t43;
    chrisD[0][0][3][2] = 0.0;
    chrisD[0][0][3][3] = t10 * t118 * t38 * t76;
    chrisD[0][1][0][0] = 0.0;
    chrisD[0][1][0][1] = -t126;
    chrisD[0][1][0][2] = 0.0;
    chrisD[0][1][0][3] = -t130;
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
    chrisD[0][3][0][1] = -t130;
    chrisD[0][3][0][2] = 0.0;
    chrisD[0][3][0][3] = t134;
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
    chrisD[1][0][0][1] = -t126;
    chrisD[1][0][0][2] = 0.0;
    chrisD[1][0][0][3] = -t130;
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
    chrisD[1][1][1][1] = t146;
    chrisD[1][1][1][2] = 0.0;
    chrisD[1][1][1][3] = -t157;
    chrisD[1][1][2][0] = 0.0;
    chrisD[1][1][2][1] = 0.0;
    chrisD[1][1][2][2] = 0.0;
    chrisD[1][1][2][3] = 0.0;
    chrisD[1][1][3][0] = 0.0;
    chrisD[1][1][3][1] = -t162;
    chrisD[1][1][3][2] = 0.0;
    chrisD[1][1][3][3] = -t169;
    chrisD[1][2][0][0] = 0.0;
    chrisD[1][2][0][1] = 0.0;
    chrisD[1][2][0][2] = 0.0;
    chrisD[1][2][0][3] = 0.0;
    chrisD[1][2][1][0] = 0.0;
    chrisD[1][2][1][1] = 0.0;
    chrisD[1][2][1][2] = 0.0;
    chrisD[1][2][1][3] = 0.0;
    chrisD[1][2][2][0] = 0.0;
    chrisD[1][2][2][1] = -t176;
    chrisD[1][2][2][2] = 0.0;
    chrisD[1][2][2][3] = t185;
    chrisD[1][2][3][0] = 0.0;
    chrisD[1][2][3][1] = 0.0;
    chrisD[1][2][3][2] = 0.0;
    chrisD[1][2][3][3] = 0.0;
    chrisD[1][3][0][0] = 0.0;
    chrisD[1][3][0][1] = 0.0;
    chrisD[1][3][0][2] = 0.0;
    chrisD[1][3][0][3] = 0.0;
    chrisD[1][3][1][0] = 0.0;
    chrisD[1][3][1][1] = t162;
    chrisD[1][3][1][2] = 0.0;
    chrisD[1][3][1][3] = t169;
    chrisD[1][3][2][0] = 0.0;
    chrisD[1][3][2][1] = 0.0;
    chrisD[1][3][2][2] = 0.0;
    chrisD[1][3][2][3] = 0.0;
    chrisD[1][3][3][0] = 0.0;
    chrisD[1][3][3][1] = t146;
    chrisD[1][3][3][2] = 0.0;
    chrisD[1][3][3][3] = -t157;
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
    chrisD[2][1][2][1] = -t176;
    chrisD[2][1][2][2] = 0.0;
    chrisD[2][1][2][3] = t185;
    chrisD[2][1][3][0] = 0.0;
    chrisD[2][1][3][1] = 0.0;
    chrisD[2][1][3][2] = 0.0;
    chrisD[2][1][3][3] = 0.0;
    chrisD[2][2][0][0] = 0.0;
    chrisD[2][2][0][1] = 0.0;
    chrisD[2][2][0][2] = 0.0;
    chrisD[2][2][0][3] = 0.0;
    chrisD[2][2][1][0] = 0.0;
    chrisD[2][2][1][1] = -t9 * t215 * t43;
    chrisD[2][2][1][2] = 0.0;
    chrisD[2][2][1][3] = -t9 * r * z * t239 * t43;
    chrisD[2][2][2][0] = 0.0;
    chrisD[2][2][2][1] = 0.0;
    chrisD[2][2][2][2] = 0.0;
    chrisD[2][2][2][3] = 0.0;
    chrisD[2][2][3][0] = 0.0;
    chrisD[2][2][3][1] = -t46 * z * (t19 - t22 + t23 + t25 - t28 - t164 - t165) * t43;
    chrisD[2][2][3][2] = 0.0;
    chrisD[2][2][3][3] = t10 * t2 * (-4.0 * t4 * t3 - 3.0 * t88 - 2.0 * t29 + t17) / t47;
    chrisD[2][3][0][0] = 0.0;
    chrisD[2][3][0][1] = 0.0;
    chrisD[2][3][0][2] = 0.0;
    chrisD[2][3][0][3] = 0.0;
    chrisD[2][3][1][0] = 0.0;
    chrisD[2][3][1][1] = 0.0;
    chrisD[2][3][1][2] = 0.0;
    chrisD[2][3][1][3] = 0.0;
    chrisD[2][3][2][0] = 0.0;
    chrisD[2][3][2][1] = t130;
    chrisD[2][3][2][2] = 0.0;
    chrisD[2][3][2][3] = -t134;
    chrisD[2][3][3][0] = 0.0;
    chrisD[2][3][3][1] = 0.0;
    chrisD[2][3][3][2] = 0.0;
    chrisD[2][3][3][3] = 0.0;
    chrisD[3][0][0][0] = 0.0;
    chrisD[3][0][0][1] = -t130;
    chrisD[3][0][0][2] = 0.0;
    chrisD[3][0][0][3] = t134;
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
    chrisD[3][1][1][1] = t162;
    chrisD[3][1][1][2] = 0.0;
    chrisD[3][1][1][3] = t169;
    chrisD[3][1][2][0] = 0.0;
    chrisD[3][1][2][1] = 0.0;
    chrisD[3][1][2][2] = 0.0;
    chrisD[3][1][2][3] = 0.0;
    chrisD[3][1][3][0] = 0.0;
    chrisD[3][1][3][1] = t146;
    chrisD[3][1][3][2] = 0.0;
    chrisD[3][1][3][3] = -t157;
    chrisD[3][2][0][0] = 0.0;
    chrisD[3][2][0][1] = 0.0;
    chrisD[3][2][0][2] = 0.0;
    chrisD[3][2][0][3] = 0.0;
    chrisD[3][2][1][0] = 0.0;
    chrisD[3][2][1][1] = 0.0;
    chrisD[3][2][1][2] = 0.0;
    chrisD[3][2][1][3] = 0.0;
    chrisD[3][2][2][0] = 0.0;
    chrisD[3][2][2][1] = t130;
    chrisD[3][2][2][2] = 0.0;
    chrisD[3][2][2][3] = -t134;
    chrisD[3][2][3][0] = 0.0;
    chrisD[3][2][3][1] = 0.0;
    chrisD[3][2][3][2] = 0.0;
    chrisD[3][2][3][3] = 0.0;
    chrisD[3][3][0][0] = 0.0;
    chrisD[3][3][0][1] = 0.0;
    chrisD[3][3][0][2] = 0.0;
    chrisD[3][3][0][3] = 0.0;
    chrisD[3][3][1][0] = 0.0;
    chrisD[3][3][1][1] = -t146;
    chrisD[3][3][1][2] = 0.0;
    chrisD[3][3][1][3] = t157;
    chrisD[3][3][2][0] = 0.0;
    chrisD[3][3][2][1] = 0.0;
    chrisD[3][3][2][2] = 0.0;
    chrisD[3][3][2][3] = 0.0;
    chrisD[3][3][3][0] = 0.0;
    chrisD[3][3][3][1] = t162;
    chrisD[3][3][3][2] = 0.0;
    chrisD[3][3][3][3] = t169;

    return true;
}

void MetricCurzon::localToCoord(const double* pos, const double* ldir, double* dir, enum_nat_tetrad_type)
{
    double r = pos[1];

    double lambda, nu;
    calcLambdaNu(pos, lambda, nu);

    double el = exp(lambda);
    double eln = exp(lambda - nu);
    dir[0] = 1.0 / el * ldir[0];
    dir[1] = eln * ldir[1];
    dir[2] = el / r * ldir[2];
    dir[3] = eln * ldir[3];
}

void MetricCurzon::coordToLocal(const double*, const double*, double*, enum_nat_tetrad_type)
{
    fprintf(stderr, "uups... not implemented yet!\n");
    // TODO
}

bool MetricCurzon::breakCondition(const double* pos)
{
    bool br = false;

    double r = pos[1];
    double z = pos[3];
    if (r * r + z * z < M4D_METRIC_EPS) {
        br = true;
    }

    return br;
}

double MetricCurzon::testConstraint(const double y[], const double kappa)
{
    double r = y[1];

    // Scale the directions with the speed of light before doubling them !!
    double dt = y[4];
    double dr = y[5];
    double dphi = y[6];
    double dz = y[7];

    double lambda, nu;
    calcLambdaNu(y, lambda, nu);

    double sum = -kappa * mSign;
    sum += -exp(2.0 * lambda) * dt * dt + exp(2.0 * (nu - lambda)) * (dr * dr + dz * dz)
        + r * r * exp(-2.0 * lambda) * dphi * dphi;

    //  double k = exp(2.0*lambda)*dt;
    //  double h = r*r*exp(-2.0*lambda)*dphi;

    return sum;
}

bool MetricCurzon::setParam(const char* pName, double val)
{
    Metric::setParam(pName, val);

    if (strcmp(pName, "mass") == 0) {
        mMass = val;
    }
    return true;
}

bool MetricCurzon::report(const vec4, const vec4, char*& text)
{
    std::stringstream ss;
    ss << "Report for the  Curzon metric\n\tcoordinate : (t,r,phi,z)\n";
    ss << "---------------------------------------------------------------\n";
    ss.precision(DEF_FIXED_REPORT_PRECISION);
    ss.setf(std::ios::fixed);

    text = new char[ss.str().length() + 2];
    return CopyString(ss.str().c_str(), text);
}

// ********************************* protected methods *****************************

void MetricCurzon::setStandardValues()
{
    mInitPos[0] = 0.0;
    mInitPos[1] = 6.0;
    mInitPos[2] = 0.0;
    mInitPos[3] = 0.0;
    mInitDir[0] = 1.0;
    mInitDir[1] = 0.0;
    mInitDir[2] = 0.0;

    mCoordNames[0] = std::string("t");
    mCoordNames[1] = std::string("r");
    mCoordNames[2] = std::string("phi");
    mCoordNames[3] = std::string("z");
}

double MetricCurzon::calcR(const double* pos)
{
    double r = pos[1];
    double z = pos[3];
    return sqrt(r * r + z * z);
}

void MetricCurzon::calcDR(const double* pos, double& R, double& dRdr, double& dRdz)
{
    double r = pos[1];
    double z = pos[3];
    R = calcR(pos);
    dRdr = r / R;
    dRdz = z / R;
}

void MetricCurzon::calcD2R(
    const double* pos, double& R, double& dRdr, double& dRdz, double& dRdrdr, double& dRdrdz, double& dRdzdz)
{
    double r = pos[1];
    double z = pos[3];
    calcDR(pos, R, dRdr, dRdz);

    double edR3 = 1.0 / (R * R * R);
    dRdrdr = z * z * edR3;
    dRdrdz = -r * z * edR3;
    dRdzdz = r * r * edR3;
}

void MetricCurzon::calcLambdaNu(const double* pos, double& lambda, double& nu)
{
    double r = pos[1];
    double R = calcR(pos);

    double R4 = R * R;
    R4 *= R4;

    lambda = -mMass / R;
    nu = -0.5 * mMass * mMass * r * r / R4;
}

void MetricCurzon::calcDLN(const double* pos, double& lambda, double& nu, double&, double&, double& dndr, double& dndz)
{
    double r = pos[1];
    // double z = pos[3];
    double R = calcR(pos);

    double dRdr, dRdz;
    calcDR(pos, R, dRdr, dRdz);
    calcLambdaNu(pos, lambda, nu);

    double edR2 = 1.0 / (R * R);
    // dldr = mMass*dRdr*edR2;
    // dldz = mMass*dRdz*edR2;

    double edR5 = edR2 / (R * R * R);
    dndr = mMass * mMass * r * (-R + 2.0 * r * dRdr) * edR5;
    dndz = 2.0 * mMass * mMass * r * r * dRdz * edR5;
}

void MetricCurzon::calcD2LN(const double* pos, double&, double&, double&, double&, double&, double&, double& dldrdr,
    double& dldrdz, double& dldzdz, double& dndrdr, double& dndrdz, double& dndzdz)
{
    double r = pos[1];

    double R, dRdr, dRdz, dRdrdr, dRdrdz, dRdzdz;
    calcD2R(pos, R, dRdr, dRdz, dRdrdr, dRdrdz, dRdzdz);

    double edR3 = 1.0 / (R * R * R);
    dldrdr = mMass * (-2.0 * dRdr * dRdr + R * dRdrdr) * edR3;
    dldrdz = -mMass * (2.0 * dRdr * dRdz - dRdrdz * R) * edR3;
    dldzdz = -mMass * (2.0 * dRdz * dRdz - dRdzdz * R) * edR3;

    double edR6 = edR3 * edR3;
    dndrdr
        = mMass * mMass * (-R * R + 8.0 * r * dRdr * R - 10.0 * r * r * dRdr * dRdr + 2.0 * r * r * dRdrdr * R) * edR6;
    dndrdz = -2.0 * mMass * mMass * r * (-2 * dRdz * r + 5.0 * r * dRdr * dRdz - r * dRdrdz * R) * edR6;
    dndzdz = -2.0 * mMass * mMass * r * r * (5.0 * dRdz * dRdz - dRdzdz * R) / edR6;
}

} // end namespace m4d
