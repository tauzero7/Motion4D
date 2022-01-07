/**
 * @file    m4dMetricBesselGravWaveCart.cpp
 * @author  Heiko Munz
 *
 * This file is part of the m4d-library.
 */
#include "m4dMetricBesselGravWaveCart.h"

namespace m4d {

MetricBesselGravWaveCart::MetricBesselGravWaveCart(double C)
{
    mMetricName = "BesselGravWaveCart";
    setCoordType(enum_coordinate_cartesian); // cartesian:  [t, x, y, z]

    mPhysicalUnits = enum_physical_constants_geom;
    mSpeedOfLight = 1.0;

    addParam("c", C);
    mC = C;

    // mLocTeds.push_back(enum_nat_tetrad_static);

    mDrawTypes.push_back(enum_draw_twoplusone);

    setStandardValues();
}

MetricBesselGravWaveCart::~MetricBesselGravWaveCart() {}

// *********************************** public methods ******************************

bool MetricBesselGravWaveCart::calculateMetric(const double* pos)
{
    double t = pos[0];
    double x = pos[1];
    double y = pos[2];

    if (x == 0.0 && y == 0.0) {
        double t1 = cos(t);
        double t2 = mC * t1;
        double t3 = exp(t2);
        double t4 = t3 * t3;
        double t6 = 0.2e1 * t2;
        double t7 = exp(-t6);
        double t8 = exp(t6);

        g_compts[0][0] = -0.1e1 / t4;
        g_compts[0][1] = 0.0e0;
        g_compts[0][2] = 0.0e0;
        g_compts[0][3] = 0.0e0;
        g_compts[1][0] = 0.0e0;
        g_compts[1][1] = t7;
        g_compts[1][2] = 0.0e0;
        g_compts[1][3] = 0.0e0;
        g_compts[2][0] = 0.0e0;
        g_compts[2][1] = 0.0e0;
        g_compts[2][2] = t7;
        g_compts[2][3] = 0.0e0;
        g_compts[3][0] = 0.0e0;
        g_compts[3][1] = 0.0e0;
        g_compts[3][2] = 0.0e0;
        g_compts[3][3] = t8;
    }
    else {
        double t1 = mC * mC;
        double t2 = x * x;
        double t3 = y * y;
        double t4 = t2 + t3;
        double t5 = sqrt(t4);
        // double t7 = BesselJ(0, t5);
        double t7 = gsl_sf_bessel_J0(t5);
        double t8 = t7 * t7;
        // double t9 = BesselJ(1, t5);
        double t9 = gsl_sf_bessel_J1(t5);
        double t10 = t9 * t9;
        double t14 = cos(t);
        double t15 = t14 * t14;
        double t21 = exp(t1 * t5 * (t5 * (t8 + t10) - 0.2e1 * t7 * t9 * t15) / 0.2e1);
        double t22 = t21 * t21;
        double t25 = exp(mC * t7 * t14);
        double t26 = t25 * t25;
        double t27 = 0.1e1 / t26;
        double t28 = t22 * t27;
        double t29 = 0.1e1 / t4;
        double t32 = t27 * t29;
        double t40 = t28 * t29 * x * y - t32 * x * y;

        g_compts[0][0] = -t28;
        g_compts[0][1] = 0.0e0;
        g_compts[0][2] = 0.0e0;
        g_compts[0][3] = 0.0e0;
        g_compts[1][0] = 0.0e0;
        g_compts[1][1] = t28 * t29 * t2 + t32 * t3;
        g_compts[1][2] = t40;
        g_compts[1][3] = 0.0e0;
        g_compts[2][0] = 0.0e0;
        g_compts[2][1] = t40;
        g_compts[2][2] = t28 * t29 * t3 + t32 * t2;
        g_compts[2][3] = 0.0e0;
        g_compts[3][0] = 0.0e0;
        g_compts[3][1] = 0.0e0;
        g_compts[3][2] = 0.0e0;
        g_compts[3][3] = t26;
    }

    return true;
}

bool MetricBesselGravWaveCart::calculateChristoffels(const double* pos)
{
    double t = pos[0];
    double x = pos[1];
    double y = pos[2];

    if (x == 0.0 && y == 0.0) {
        double t1 = sin(t);
        double t2 = mC * t1;
        double t3 = cos(t);
        double t5 = exp(mC * t3);
        double t6 = t5 * t5;
        double t7 = t6 * t6;

        christoffel[0][0][0] = t2;
        christoffel[0][0][1] = 0.0e0;
        christoffel[0][0][2] = 0.0e0;
        christoffel[0][0][3] = 0.0e0;
        christoffel[0][1][0] = 0.0e0;
        christoffel[0][1][1] = t2;
        christoffel[0][1][2] = 0.0e0;
        christoffel[0][1][3] = 0.0e0;
        christoffel[0][2][0] = 0.0e0;
        christoffel[0][2][1] = 0.0e0;
        christoffel[0][2][2] = t2;
        christoffel[0][2][3] = 0.0e0;
        christoffel[0][3][0] = 0.0e0;
        christoffel[0][3][1] = 0.0e0;
        christoffel[0][3][2] = 0.0e0;
        christoffel[0][3][3] = -t2;
        christoffel[1][0][0] = 0.0e0;
        christoffel[1][0][1] = t2;
        christoffel[1][0][2] = 0.0e0;
        christoffel[1][0][3] = 0.0e0;
        christoffel[1][1][0] = t2;
        christoffel[1][1][1] = 0.0e0;
        christoffel[1][1][2] = 0.0e0;
        christoffel[1][1][3] = 0.0e0;
        christoffel[1][2][0] = 0.0e0;
        christoffel[1][2][1] = 0.0e0;
        christoffel[1][2][2] = 0.0e0;
        christoffel[1][2][3] = 0.0e0;
        christoffel[1][3][0] = 0.0e0;
        christoffel[1][3][1] = 0.0e0;
        christoffel[1][3][2] = 0.0e0;
        christoffel[1][3][3] = 0.0e0;
        christoffel[2][0][0] = 0.0e0;
        christoffel[2][0][1] = 0.0e0;
        christoffel[2][0][2] = t2;
        christoffel[2][0][3] = 0.0e0;
        christoffel[2][1][0] = 0.0e0;
        christoffel[2][1][1] = 0.0e0;
        christoffel[2][1][2] = 0.0e0;
        christoffel[2][1][3] = 0.0e0;
        christoffel[2][2][0] = t2;
        christoffel[2][2][1] = 0.0e0;
        christoffel[2][2][2] = 0.0e0;
        christoffel[2][2][3] = 0.0e0;
        christoffel[2][3][0] = 0.0e0;
        christoffel[2][3][1] = 0.0e0;
        christoffel[2][3][2] = 0.0e0;
        christoffel[2][3][3] = 0.0e0;
        christoffel[3][0][0] = 0.0e0;
        christoffel[3][0][1] = 0.0e0;
        christoffel[3][0][2] = 0.0e0;
        christoffel[3][0][3] = -t2;
        christoffel[3][1][0] = 0.0e0;
        christoffel[3][1][1] = 0.0e0;
        christoffel[3][1][2] = 0.0e0;
        christoffel[3][1][3] = 0.0e0;
        christoffel[3][2][0] = 0.0e0;
        christoffel[3][2][1] = 0.0e0;
        christoffel[3][2][2] = 0.0e0;
        christoffel[3][2][3] = 0.0e0;
        christoffel[3][3][0] = -t7 * mC * t1;
        christoffel[3][3][1] = 0.0e0;
        christoffel[3][3][2] = 0.0e0;
        christoffel[3][3][3] = 0.0e0;
    }
    else {
        double t1 = x * x;
        double t2 = y * y;
        double t3 = t1 + t2;
        double t4 = sqrt(t3);
        // double t5 = BesselJ(0, t4);
        double t5 = gsl_sf_bessel_J0(t4);
        double t6 = mC * t5;
        double t7 = sin(t);
        double t8 = mC * t4;
        // double t9 = BesselJ(1, t4);
        double t9 = gsl_sf_bessel_J1(t4);
        double t10 = cos(t);
        double t11 = t9 * t10;
        double t18 = t5 * t5;
        double t20 = t9 * t9;
        double t21 = t10 * t10;
        double t27 = 0.1e1 / t4;
        double t28 = (t8 * t18 + t8 * t20 * t21 - t8 * t18 * t21 + t11) * t27;
        double t29 = mC * x * t28;
        double t30 = y * mC;
        double t31 = t28 * t30;
        double t32 = t1 * mC;
        double t37 = t6 * t27;
        double t38 = (t4 + 0.2e1 * t32 * t11) * t7 * t37;
        double t39 = t7 * t5;
        double t41 = mC * mC;
        double t43 = x * t27;
        double t46 = 0.2e1 * t11 * t39 * t41 * y * t43;
        double t47 = mC * t9;
        double t53 = (0.2e1 * t47 * t10 * t2 + t4) * t7 * t37;
        double t54 = t6 * t7;
        double t64 = t41 * (t18 * t1 + t18 * t2 + t20 * t1 + t20 * t2 - 0.2e1 * t4 * t5 * t9 * t21);
        double t65 = exp(-t64);
        double t67 = t65 * mC * t5;
        double t68 = exp(t64);
        double t69 = t1 * t1;
        double t70 = t68 * t69;
        double t71 = t47 * t10;
        double t72 = t70 * t71;
        double t73 = 0.2e1 * t72;
        double t74 = t68 * t1;
        double t76 = t11 * t2;
        double t77 = t74 * mC * t76;
        double t78 = 0.2e1 * t77;
        double t79 = t74 * t4;
        double t80 = t2 * t4;
        double t84 = 0.1e1 / t4 / t3;
        double t87 = t2 * t2;
        double t88 = t87 * t68;
        double t89 = t88 * t71;
        double t90 = 0.2e1 * t89;
        double t92 = t47 * t10 * t87;
        double t93 = t68 * t4;
        double t94 = t93 * t2;
        double t95 = t74 * t41;
        double t96 = t4 * t20;
        double t97 = t21 * t2;
        double t99 = t95 * t96 * t97;
        double t100 = t32 * t76;
        double t101 = t4 * t18;
        double t103 = t95 * t101 * t2;
        double t105 = t95 * t101 * t97;
        double t106 = 0.3e1 * t77;
        double t108 = t41 * t4 * t18;
        double t109 = t70 * t108;
        double t110 = t70 * t41;
        double t111 = t96 * t21;
        double t112 = t110 * t111;
        double t113 = t101 * t21;
        double t114 = t110 * t113;
        double t115 = t90 - t92 + t94 + t99 - t80 - t100 + t103 - t105 + t106 + t72 + t109 + t112 - t114;
        double t117 = t3 * t3;
        double t119 = 0.1e1 / t4 / t117;
        double t120 = x * t119;
        double t122 = t92 - t99 - t94 + t100 + t105 - t103 + t80 + t77 + t72 - t109 + t114 - t112;
        double t124 = y * t119;
        double t136 = t65 * x * t30 * t39 * (0.2e1 * t74 * t71 + 0.2e1 * t68 * mC * t76 + t93 - t4) * t84;
        double t138 = t47 * t10 * t69;
        double t139 = t4 * t1;
        double t140 = t89 + t99 + t100 - t105 + t77 + t103 + t138 - t114 + t139 + t112 + t109 - t79;
        double t142 = t140 * t65 * t124;
        double t143 = t88 * t108;
        double t144 = t88 * t41;
        double t145 = t144 * t113;
        double t146 = t144 * t111;
        double t147 = t92 + t143 - t145 + t146 - t105 + t99 + t80 + t100 + t103 + t77 - t94 + t72;
        double t149 = t147 * t65 * t120;
        double t150 = t43 * t10;
        double t151 = t47 * t150;
        double t156 = t145 - t143 - t146 + t89 - t103 - t99 + t77 + t105 + t100 + t139 + t138 - t79;
        double t159 = -t89 - t143 - t146 + t145 - t99 + t105 + t100 - t103 - t106 + t138 - t73 - t79 + t139;
        double t164 = t47 * t27 * y * t10;
        double t165 = mC * t18;
        double t168 = mC * t20;
        double t179 = exp(
            -mC * (t165 * t1 + t165 * t2 + t168 * t1 + t168 * t2 - 0.2e1 * t8 * t5 * t9 * t21 - 0.4e1 * t5 * t10));
        double t180 = t179 * mC;

        christoffel[0][0][0] = t6 * t7 * (0.2e1 * t8 * t11 + 0.1e1);
        christoffel[0][0][1] = t29;
        christoffel[0][0][2] = t31;
        christoffel[0][0][3] = 0.0e0;
        christoffel[0][1][0] = t29;
        christoffel[0][1][1] = t38;
        christoffel[0][1][2] = t46;
        christoffel[0][1][3] = 0.0e0;
        christoffel[0][2][0] = t31;
        christoffel[0][2][1] = t46;
        christoffel[0][2][2] = t53;
        christoffel[0][2][3] = 0.0e0;
        christoffel[0][3][0] = 0.0e0;
        christoffel[0][3][1] = 0.0e0;
        christoffel[0][3][2] = 0.0e0;
        christoffel[0][3][3] = -t54;
        christoffel[1][0][0] = t29;
        christoffel[1][0][1] = t38;
        christoffel[1][0][2] = t46;
        christoffel[1][0][3] = 0.0e0;
        christoffel[1][1][0] = t67 * t7 * (t73 + t78 + t79 + t80) * t84;
        christoffel[1][1][1] = t115 * t65 * t120;
        christoffel[1][1][2] = -t122 * t65 * t124;
        christoffel[1][1][3] = 0.0e0;
        christoffel[1][2][0] = t136;
        christoffel[1][2][1] = t142;
        christoffel[1][2][2] = t149;
        christoffel[1][2][3] = 0.0e0;
        christoffel[1][3][0] = 0.0e0;
        christoffel[1][3][1] = 0.0e0;
        christoffel[1][3][2] = 0.0e0;
        christoffel[1][3][3] = -t151;
        christoffel[2][0][0] = t31;
        christoffel[2][0][1] = t46;
        christoffel[2][0][2] = t53;
        christoffel[2][0][3] = 0.0e0;
        christoffel[2][1][0] = t136;
        christoffel[2][1][1] = t142;
        christoffel[2][1][2] = t149;
        christoffel[2][1][3] = 0.0e0;
        christoffel[2][2][0] = t67 * t7 * (t78 + t90 + t94 + t139) * t84;
        christoffel[2][2][1] = -t156 * t65 * t120;
        christoffel[2][2][2] = -t159 * t65 * t124;
        christoffel[2][2][3] = 0.0e0;
        christoffel[2][3][0] = 0.0e0;
        christoffel[2][3][1] = 0.0e0;
        christoffel[2][3][2] = 0.0e0;
        christoffel[2][3][3] = -t164;
        christoffel[3][0][0] = 0.0e0;
        christoffel[3][0][1] = 0.0e0;
        christoffel[3][0][2] = 0.0e0;
        christoffel[3][0][3] = -t54;
        christoffel[3][1][0] = 0.0e0;
        christoffel[3][1][1] = 0.0e0;
        christoffel[3][1][2] = 0.0e0;
        christoffel[3][1][3] = -t151;
        christoffel[3][2][0] = 0.0e0;
        christoffel[3][2][1] = 0.0e0;
        christoffel[3][2][2] = 0.0e0;
        christoffel[3][2][3] = -t164;
        christoffel[3][3][0] = -t180 * t39;
        christoffel[3][3][1] = t180 * t9 * t150;
        christoffel[3][3][2] = t179 * y * mC * t11 * t27;
        christoffel[3][3][3] = 0.0e0;
    }

    return true;
}

bool MetricBesselGravWaveCart::calculateChrisD(const double* pos)
{
    double t = pos[0];
    double x = pos[1];
    double y = pos[2];

    if (x == 0.0 && y == 0.0) {
        double t1 = cos(t);
        double t2 = mC * t1;
        double t3 = mC * mC;
        double t4 = t2 / 0.2e1;
        double t5 = t1 * t1;
        double t6 = t3 * t5;
        double t7 = t3 + t4 - t6;
        double t8 = -t4 - t6 + t3;
        double t9 = exp(t2);
        double t10 = t9 * t9;
        double t11 = t10 * t10;
        double t12 = t3 * t11;
        double t17 = t11 * mC * t1;
        double t19 = t17 / 0.2e1;

        chrisD[0][0][0][0] = t2;
        chrisD[0][0][0][1] = 0.0e0;
        chrisD[0][0][0][2] = 0.0e0;
        chrisD[0][0][0][3] = 0.0e0;
        chrisD[0][0][1][0] = 0.0e0;
        chrisD[0][0][1][1] = t7;
        chrisD[0][0][1][2] = 0.0e0;
        chrisD[0][0][1][3] = 0.0e0;
        chrisD[0][0][2][0] = 0.0e0;
        chrisD[0][0][2][1] = 0.0e0;
        chrisD[0][0][2][2] = t7;
        chrisD[0][0][2][3] = 0.0e0;
        chrisD[0][0][3][0] = 0.0e0;
        chrisD[0][0][3][1] = 0.0e0;
        chrisD[0][0][3][2] = 0.0e0;
        chrisD[0][0][3][3] = 0.0e0;
        chrisD[0][1][0][0] = 0.0e0;
        chrisD[0][1][0][1] = t7;
        chrisD[0][1][0][2] = 0.0e0;
        chrisD[0][1][0][3] = 0.0e0;
        chrisD[0][1][1][0] = t2;
        chrisD[0][1][1][1] = 0.0e0;
        chrisD[0][1][1][2] = 0.0e0;
        chrisD[0][1][1][3] = 0.0e0;
        chrisD[0][1][2][0] = 0.0e0;
        chrisD[0][1][2][1] = 0.0e0;
        chrisD[0][1][2][2] = 0.0e0;
        chrisD[0][1][2][3] = 0.0e0;
        chrisD[0][1][3][0] = 0.0e0;
        chrisD[0][1][3][1] = 0.0e0;
        chrisD[0][1][3][2] = 0.0e0;
        chrisD[0][1][3][3] = 0.0e0;
        chrisD[0][2][0][0] = 0.0e0;
        chrisD[0][2][0][1] = 0.0e0;
        chrisD[0][2][0][2] = t7;
        chrisD[0][2][0][3] = 0.0e0;
        chrisD[0][2][1][0] = 0.0e0;
        chrisD[0][2][1][1] = 0.0e0;
        chrisD[0][2][1][2] = 0.0e0;
        chrisD[0][2][1][3] = 0.0e0;
        chrisD[0][2][2][0] = t2;
        chrisD[0][2][2][1] = 0.0e0;
        chrisD[0][2][2][2] = 0.0e0;
        chrisD[0][2][2][3] = 0.0e0;
        chrisD[0][2][3][0] = 0.0e0;
        chrisD[0][2][3][1] = 0.0e0;
        chrisD[0][2][3][2] = 0.0e0;
        chrisD[0][2][3][3] = 0.0e0;
        chrisD[0][3][0][0] = 0.0e0;
        chrisD[0][3][0][1] = 0.0e0;
        chrisD[0][3][0][2] = 0.0e0;
        chrisD[0][3][0][3] = 0.0e0;
        chrisD[0][3][1][0] = 0.0e0;
        chrisD[0][3][1][1] = 0.0e0;
        chrisD[0][3][1][2] = 0.0e0;
        chrisD[0][3][1][3] = 0.0e0;
        chrisD[0][3][2][0] = 0.0e0;
        chrisD[0][3][2][1] = 0.0e0;
        chrisD[0][3][2][2] = 0.0e0;
        chrisD[0][3][2][3] = 0.0e0;
        chrisD[0][3][3][0] = -t2;
        chrisD[0][3][3][1] = 0.0e0;
        chrisD[0][3][3][2] = 0.0e0;
        chrisD[0][3][3][3] = 0.0e0;
        chrisD[1][0][0][0] = 0.0e0;
        chrisD[1][0][0][1] = t7;
        chrisD[1][0][0][2] = 0.0e0;
        chrisD[1][0][0][3] = 0.0e0;
        chrisD[1][0][1][0] = t2;
        chrisD[1][0][1][1] = 0.0e0;
        chrisD[1][0][1][2] = 0.0e0;
        chrisD[1][0][1][3] = 0.0e0;
        chrisD[1][0][2][0] = 0.0e0;
        chrisD[1][0][2][1] = 0.0e0;
        chrisD[1][0][2][2] = 0.0e0;
        chrisD[1][0][2][3] = 0.0e0;
        chrisD[1][0][3][0] = 0.0e0;
        chrisD[1][0][3][1] = 0.0e0;
        chrisD[1][0][3][2] = 0.0e0;
        chrisD[1][0][3][3] = 0.0e0;
        chrisD[1][1][0][0] = t2;
        chrisD[1][1][0][1] = 0.0e0;
        chrisD[1][1][0][2] = 0.0e0;
        chrisD[1][1][0][3] = 0.0e0;
        chrisD[1][1][1][0] = 0.0e0;
        chrisD[1][1][1][1] = t7;
        chrisD[1][1][1][2] = 0.0e0;
        chrisD[1][1][1][3] = 0.0e0;
        chrisD[1][1][2][0] = 0.0e0;
        chrisD[1][1][2][1] = 0.0e0;
        chrisD[1][1][2][2] = t8;
        chrisD[1][1][2][3] = 0.0e0;
        chrisD[1][1][3][0] = 0.0e0;
        chrisD[1][1][3][1] = 0.0e0;
        chrisD[1][1][3][2] = 0.0e0;
        chrisD[1][1][3][3] = 0.0e0;
        chrisD[1][2][0][0] = 0.0e0;
        chrisD[1][2][0][1] = 0.0e0;
        chrisD[1][2][0][2] = 0.0e0;
        chrisD[1][2][0][3] = 0.0e0;
        chrisD[1][2][1][0] = 0.0e0;
        chrisD[1][2][1][1] = 0.0e0;
        chrisD[1][2][1][2] = t4;
        chrisD[1][2][1][3] = 0.0e0;
        chrisD[1][2][2][0] = 0.0e0;
        chrisD[1][2][2][1] = t4;
        chrisD[1][2][2][2] = 0.0e0;
        chrisD[1][2][2][3] = 0.0e0;
        chrisD[1][2][3][0] = 0.0e0;
        chrisD[1][2][3][1] = 0.0e0;
        chrisD[1][2][3][2] = 0.0e0;
        chrisD[1][2][3][3] = 0.0e0;
        chrisD[1][3][0][0] = 0.0e0;
        chrisD[1][3][0][1] = 0.0e0;
        chrisD[1][3][0][2] = 0.0e0;
        chrisD[1][3][0][3] = 0.0e0;
        chrisD[1][3][1][0] = 0.0e0;
        chrisD[1][3][1][1] = 0.0e0;
        chrisD[1][3][1][2] = 0.0e0;
        chrisD[1][3][1][3] = 0.0e0;
        chrisD[1][3][2][0] = 0.0e0;
        chrisD[1][3][2][1] = 0.0e0;
        chrisD[1][3][2][2] = 0.0e0;
        chrisD[1][3][2][3] = 0.0e0;
        chrisD[1][3][3][0] = 0.0e0;
        chrisD[1][3][3][1] = -t4;
        chrisD[1][3][3][2] = 0.0e0;
        chrisD[1][3][3][3] = 0.0e0;
        chrisD[2][0][0][0] = 0.0e0;
        chrisD[2][0][0][1] = 0.0e0;
        chrisD[2][0][0][2] = t7;
        chrisD[2][0][0][3] = 0.0e0;
        chrisD[2][0][1][0] = 0.0e0;
        chrisD[2][0][1][1] = 0.0e0;
        chrisD[2][0][1][2] = 0.0e0;
        chrisD[2][0][1][3] = 0.0e0;
        chrisD[2][0][2][0] = t2;
        chrisD[2][0][2][1] = 0.0e0;
        chrisD[2][0][2][2] = 0.0e0;
        chrisD[2][0][2][3] = 0.0e0;
        chrisD[2][0][3][0] = 0.0e0;
        chrisD[2][0][3][1] = 0.0e0;
        chrisD[2][0][3][2] = 0.0e0;
        chrisD[2][0][3][3] = 0.0e0;
        chrisD[2][1][0][0] = 0.0e0;
        chrisD[2][1][0][1] = 0.0e0;
        chrisD[2][1][0][2] = 0.0e0;
        chrisD[2][1][0][3] = 0.0e0;
        chrisD[2][1][1][0] = 0.0e0;
        chrisD[2][1][1][1] = 0.0e0;
        chrisD[2][1][1][2] = t4;
        chrisD[2][1][1][3] = 0.0e0;
        chrisD[2][1][2][0] = 0.0e0;
        chrisD[2][1][2][1] = t4;
        chrisD[2][1][2][2] = 0.0e0;
        chrisD[2][1][2][3] = 0.0e0;
        chrisD[2][1][3][0] = 0.0e0;
        chrisD[2][1][3][1] = 0.0e0;
        chrisD[2][1][3][2] = 0.0e0;
        chrisD[2][1][3][3] = 0.0e0;
        chrisD[2][2][0][0] = t2;
        chrisD[2][2][0][1] = 0.0e0;
        chrisD[2][2][0][2] = 0.0e0;
        chrisD[2][2][0][3] = 0.0e0;
        chrisD[2][2][1][0] = 0.0e0;
        chrisD[2][2][1][1] = t8;
        chrisD[2][2][1][2] = 0.0e0;
        chrisD[2][2][1][3] = 0.0e0;
        chrisD[2][2][2][0] = 0.0e0;
        chrisD[2][2][2][1] = 0.0e0;
        chrisD[2][2][2][2] = t7;
        chrisD[2][2][2][3] = 0.0e0;
        chrisD[2][2][3][0] = 0.0e0;
        chrisD[2][2][3][1] = 0.0e0;
        chrisD[2][2][3][2] = 0.0e0;
        chrisD[2][2][3][3] = 0.0e0;
        chrisD[2][3][0][0] = 0.0e0;
        chrisD[2][3][0][1] = 0.0e0;
        chrisD[2][3][0][2] = 0.0e0;
        chrisD[2][3][0][3] = 0.0e0;
        chrisD[2][3][1][0] = 0.0e0;
        chrisD[2][3][1][1] = 0.0e0;
        chrisD[2][3][1][2] = 0.0e0;
        chrisD[2][3][1][3] = 0.0e0;
        chrisD[2][3][2][0] = 0.0e0;
        chrisD[2][3][2][1] = 0.0e0;
        chrisD[2][3][2][2] = 0.0e0;
        chrisD[2][3][2][3] = 0.0e0;
        chrisD[2][3][3][0] = 0.0e0;
        chrisD[2][3][3][1] = 0.0e0;
        chrisD[2][3][3][2] = -t4;
        chrisD[2][3][3][3] = 0.0e0;
        chrisD[3][0][0][0] = 0.0e0;
        chrisD[3][0][0][1] = 0.0e0;
        chrisD[3][0][0][2] = 0.0e0;
        chrisD[3][0][0][3] = 0.0e0;
        chrisD[3][0][1][0] = 0.0e0;
        chrisD[3][0][1][1] = 0.0e0;
        chrisD[3][0][1][2] = 0.0e0;
        chrisD[3][0][1][3] = 0.0e0;
        chrisD[3][0][2][0] = 0.0e0;
        chrisD[3][0][2][1] = 0.0e0;
        chrisD[3][0][2][2] = 0.0e0;
        chrisD[3][0][2][3] = 0.0e0;
        chrisD[3][0][3][0] = -t2;
        chrisD[3][0][3][1] = 0.0e0;
        chrisD[3][0][3][2] = 0.0e0;
        chrisD[3][0][3][3] = 0.0e0;
        chrisD[3][1][0][0] = 0.0e0;
        chrisD[3][1][0][1] = 0.0e0;
        chrisD[3][1][0][2] = 0.0e0;
        chrisD[3][1][0][3] = 0.0e0;
        chrisD[3][1][1][0] = 0.0e0;
        chrisD[3][1][1][1] = 0.0e0;
        chrisD[3][1][1][2] = 0.0e0;
        chrisD[3][1][1][3] = 0.0e0;
        chrisD[3][1][2][0] = 0.0e0;
        chrisD[3][1][2][1] = 0.0e0;
        chrisD[3][1][2][2] = 0.0e0;
        chrisD[3][1][2][3] = 0.0e0;
        chrisD[3][1][3][0] = 0.0e0;
        chrisD[3][1][3][1] = -t4;
        chrisD[3][1][3][2] = 0.0e0;
        chrisD[3][1][3][3] = 0.0e0;
        chrisD[3][2][0][0] = 0.0e0;
        chrisD[3][2][0][1] = 0.0e0;
        chrisD[3][2][0][2] = 0.0e0;
        chrisD[3][2][0][3] = 0.0e0;
        chrisD[3][2][1][0] = 0.0e0;
        chrisD[3][2][1][1] = 0.0e0;
        chrisD[3][2][1][2] = 0.0e0;
        chrisD[3][2][1][3] = 0.0e0;
        chrisD[3][2][2][0] = 0.0e0;
        chrisD[3][2][2][1] = 0.0e0;
        chrisD[3][2][2][2] = 0.0e0;
        chrisD[3][2][2][3] = 0.0e0;
        chrisD[3][2][3][0] = 0.0e0;
        chrisD[3][2][3][1] = 0.0e0;
        chrisD[3][2][3][2] = -t4;
        chrisD[3][2][3][3] = 0.0e0;
        chrisD[3][3][0][0] = -0.4e1 * t12 * t5 + 0.4e1 * t12 - t17;
        chrisD[3][3][0][1] = 0.0e0;
        chrisD[3][3][0][2] = 0.0e0;
        chrisD[3][3][0][3] = 0.0e0;
        chrisD[3][3][1][0] = 0.0e0;
        chrisD[3][3][1][1] = t19;
        chrisD[3][3][1][2] = 0.0e0;
        chrisD[3][3][1][3] = 0.0e0;
        chrisD[3][3][2][0] = 0.0e0;
        chrisD[3][3][2][1] = 0.0e0;
        chrisD[3][3][2][2] = t19;
        chrisD[3][3][2][3] = 0.0e0;
        chrisD[3][3][3][0] = 0.0e0;
        chrisD[3][3][3][1] = 0.0e0;
        chrisD[3][3][3][2] = 0.0e0;
        chrisD[3][3][3][3] = 0.0e0;
    }
    else {
        double t1 = x * x;
        double t2 = y * y;
        double t3 = t1 + t2;
        double t4 = sqrt(t3);
        // double t5 = BesselJ(0, t4);
        double t5 = gsl_sf_bessel_J0(t4);
        double t6 = mC * t5;
        double t7 = mC * t4;
        // double t8 = BesselJ(1, t4);
        double t8 = gsl_sf_bessel_J1(t4);
        double t9 = cos(t);
        double t10 = t9 * t9;
        double t11 = t8 * t10;
        double t18 = mC * x;
        double t19 = sin(t);
        double t20 = t8 * t8;
        double t21 = t20 * t9;
        double t22 = t7 * t21;
        double t24 = t5 * t5;
        double t25 = mC * t24;
        double t26 = t4 * t9;
        double t31 = 0.1e1 / t4;
        double t32 = t19 * (-0.2e1 * t22 - t8 + 0.2e1 * t25 * t26) * t31;
        double t33 = t18 * t32;
        double t34 = mC * y;
        double t35 = t34 * t32;
        double t36 = t24 * t1;
        double t37 = t7 * t36;
        double t38 = t24 * t2;
        double t39 = t7 * t38;
        double t40 = t1 * mC;
        double t41 = t4 * t20;
        double t42 = t41 * t10;
        double t43 = t40 * t42;
        double t44 = t20 * t10;
        double t45 = t44 * t2;
        double t46 = t7 * t45;
        double t47 = t24 * t10;
        double t48 = t47 * t1;
        double t49 = t7 * t48;
        double t50 = t47 * t2;
        double t51 = t7 * t50;
        double t52 = t1 * t8;
        double t53 = t52 * t9;
        double t54 = t8 * t9;
        double t55 = t54 * t2;
        double t56 = t1 * t1;
        double t57 = t56 * mC;
        double t58 = t5 * t8;
        double t61 = t2 * t8;
        double t62 = t61 * t5;
        double t64 = 0.2e1 * t40 * t62;
        double t65 = t5 * t10;
        double t66 = t65 * t8;
        double t69 = t40 * t5;
        double t70 = t11 * t2;
        double t72 = 0.4e1 * t69 * t70;
        double t73 = t1 * t9;
        double t74 = t4 * t5;
        double t76 = -t37 - t39 + t43 - t46 + t49 + t51 + t53 - t55 + 0.2e1 * t57 * t58 + t64 - 0.4e1 * t57 * t66 - t72
            - t73 * t74;
        double t79 = 0.1e1 / t4 / t3;
        double t80 = t76 * mC * t79;
        double t92 = t9 * t5;
        double t93 = t92 * t4;
        double t98 = t18 * y
            * (-0.2e1 * t6 * t52 - 0.2e1 * t6 * t61 - 0.2e1 * t7 * t44 + 0.4e1 * t6 * t11 * t1 + 0.4e1 * t6 * t70 + t93
                - 0.2e1 * t54)
            * t79;
        double t99 = t2 * t2;
        double t100 = t8 * t99;
        double t106 = t2 * t9;
        double t107 = t106 * t74;
        double t108 = -t37 - t39 - t43 + t46 + t49 + t51 - t53 + t55 + t64 + 0.2e1 * t6 * t100 - t72
            - 0.4e1 * t6 * t11 * t99 - t107;
        double t110 = t108 * mC * t79;
        double t111 = mC * t8;
        double t112 = t111 * t1;
        double t118 = t6 * (-0.2e1 * t112 + 0.4e1 * t40 * t11 + t26) * t31;
        double t120 = t9 * t24 * t4;
        double t122 = 0.2e1 * t40 * t120;
        double t123 = t40 * t9;
        double t124 = t41 * t123;
        double t125 = 0.2e1 * t124;
        double t127 = 0.4e1 * t6 * t55;
        double t130 = t19 * x;
        double t131 = t130 * t79;
        double t132 = (t122 - t52 - t125 + t127 - t61) * mC * t131;
        double t133 = t1 * t5;
        double t134 = t111 * t9;
        double t136 = 0.4e1 * t133 * t134;
        double t140 = t19 * y * t79;
        double t141 = (-t136 + t122 - t52 - t125 - t61) * mC * t140;
        double t142 = mC * mC;
        double t143 = t58 * t142;
        double t150 = 0.2e1 * t143 * y * x * (-0.1e1 + 0.2e1 * t10) * t31;
        double t151 = t41 * t1;
        double t152 = t36 * t4;
        double t153 = t52 * t5;
        double t158 = t19 * t9 * t79;
        double t160 = 0.2e1 * (-t151 + t152 - t153 + t62) * y * t142 * t158;
        double t161 = t38 * t4;
        double t162 = t2 * t20;
        double t168 = 0.2e1 * (t153 + t161 - t162 * t4 - t62) * x * t142 * t158;
        double t169 = t2 * mC;
        double t176 = t6 * (-0.2e1 * t169 * t8 + 0.4e1 * t169 * t11 + t26) * t31;
        double t177 = mC * t9;
        double t179 = 0.2e1 * t177 * t161;
        double t181 = t41 * t169 * t9;
        double t182 = 0.2e1 * t181;
        double t185 = (-t52 + t179 - t61 - t127 - t182) * mC * t131;
        double t188 = (t136 - t52 - t61 + t179 - t182) * mC * t140;
        double t189 = t6 * t9;
        double t192 = t111 * t31 * x * t19;
        double t195 = t111 * t31 * y * t19;
        double t199 = t36 + t38 + t20 * t1 + t162 - 0.2e1 * t74 * t11;
        double t200 = t142 * t199;
        double t201 = exp(t200);
        double t212 = t1 * t142;
        double t213 = t212 * t5;
        double t214 = t10 * t9;
        double t215 = t8 * t214;
        double t216 = t215 * t2;
        double t218 = 0.4e1 * t213 * t216;
        double t219 = t61 * t201;
        double t223 = 0.4e1 * t213 * t55;
        double t224 = t4 * t201;
        double t226 = t142 * t5;
        double t235 = exp(-t200);
        double t237 = t6 * t79;
        double t239 = t56 * t8;
        double t240 = t239 * t201;
        double t242 = t177 * t201;
        double t244 = 0.2e1 * t41 * t56 * t242;
        double t245 = t57 * t9;
        double t249 = 0.2e1 * t245 * t201 * t24 * t4;
        double t250 = t24 * t5;
        double t251 = t250 * t4;
        double t254 = 0.2e1 * t212 * t251 * t2;
        double t255 = t2 * t201;
        double t258 = 0.2e1 * t151 * t177 * t255;
        double t259 = t52 * t255;
        double t262 = 0.4e1 * t69 * t54 * t255;
        double t263 = t52 * t2;
        double t264 = t24 * t4;
        double t267 = 0.2e1 * t123 * t255 * t264;
        double t268 = t2 * t4;
        double t269 = t44 * t268;
        double t271 = 0.2e1 * t213 * t269;
        double t272 = t142 * t250;
        double t273 = t272 * t1;
        double t274 = t10 * t4;
        double t275 = t274 * t2;
        double t277 = 0.2e1 * t273 * t275;
        double t278 = t74 * t255;
        double t281 = t5 * t2 * t4;
        double t282 = 0.2e1 * t281;
        double t283 = t274 * t99;
        double t285 = 0.2e1 * t272 * t283;
        double t286 = t226 * t20;
        double t288 = 0.2e1 * t286 * t283;
        double t289 = t4 * t99;
        double t291 = 0.2e1 * t272 * t289;
        double t293 = t9 * t201;
        double t295 = t6 * t8 * t293 * t99;
        double t297 = t240 + t244 - t249 + t254 + t258 + t259 - t262 + t263 - t267 + t271 - t277 - 0.2e1 * t278 + t282
            + t100 - t285 + t288 + t291 - 0.4e1 * t295;
        double t300 = t3 * t3;
        double t302 = 0.1e1 / t4 / t300;
        double t303 = t18 * t302;
        double t306 = t56 * t5 * mC;
        double t307 = t54 * t201;
        double t308 = t306 * t307;
        double t310 = t133 * t4;
        double t311 = 0.2e1 * t310;
        double t312 = t133 * t224;
        double t314 = t240 - t249 + t244 + 0.4e1 * t308 + t259 + t258 - t277 + t263 + t254 - t311 + t262 + 0.2e1 * t312
            - t267 + t271 + t288 + t291 + t100 - t285;
        double t317 = t34 * t302;
        double t319 = t41 * mC;
        double t320 = -t142 * t199;
        double t321 = exp(-t320);
        double t324 = 0.2e1 * t319 * t73 * t321;
        double t325 = t264 * t321;
        double t327 = 0.2e1 * t123 * t325;
        double t328 = t52 * t321;
        double t330 = 0.4e1 * t286 * t275;
        double t331 = t61 * t321;
        double t335 = exp(t320);
        double t338 = t19 * mC * t79;
        double t340 = t212 * t99;
        double t342 = 0.2e1 * t340 * t42;
        double t343 = t142 * mC;
        double t344 = t343 * t56;
        double t346 = t54 * t24;
        double t347 = t344 * t99 * t346;
        double t348 = 0.4e1 * t347;
        double t349 = t1 * t99;
        double t350 = t349 * mC;
        double t351 = t350 * t93;
        double t352 = t56 * t1;
        double t353 = t352 * t201;
        double t354 = t142 * t4;
        double t355 = t354 * t24;
        double t356 = t353 * t355;
        double t357 = t56 * t201;
        double t358 = t357 * mC;
        double t359 = t358 * t107;
        double t361 = t57 * t55;
        double t362 = 0.3e1 * t361;
        double t363 = t343 * t1;
        double t364 = t99 * t2;
        double t365 = t363 * t364;
        double t366 = t215 * t24;
        double t367 = t365 * t366;
        double t368 = 0.2e1 * t367;
        double t369 = t343 * t352;
        double t370 = t369 * t8;
        double t371 = t106 * t24;
        double t372 = t370 * t371;
        double t373 = 0.2e1 * t372;
        double t374 = t1 * t201;
        double t375 = t374 * t99;
        double t376 = t375 * t134;
        double t377 = t20 * t8;
        double t378 = t377 * t214;
        double t379 = t378 * t2;
        double t380 = t369 * t379;
        double t381 = 0.2e1 * t380;
        double t382 = t353 * t142;
        double t383 = t65 * t61;
        double t384 = t382 * t383;
        double t386 = -t342 - t348 + t351 - t356 - 0.3e1 * t359 - t362 + t368 - t373 + t376 - t381 - 0.8e1 * t384;
        double t387 = t56 * t56;
        double t388 = t387 * t201;
        double t392 = t142 * t56;
        double t393 = t392 * t2;
        double t394 = t264 * t10;
        double t396 = 0.2e1 * t393 * t394;
        double t399 = t357 * t142;
        double t400 = t5 * t99;
        double t401 = t400 * t8;
        double t402 = t399 * t401;
        double t403 = 0.2e1 * t402;
        double t405 = t375 * t177 * t74;
        double t407 = t374 * t142;
        double t408 = t264 * t99;
        double t409 = t407 * t408;
        double t412 = 0.2e1 * t212 * t408;
        double t416 = 0.2e1 * t392 * t161;
        double t417 = t10 * t2;
        double t418 = t264 * t417;
        double t419 = t399 * t418;
        double t421 = t353 * t134;
        double t422 = t365 * t346;
        double t423 = 0.2e1 * t422;
        double t424 = -0.4e1 * t388 * t142 * t66 + t396 + 0.2e1 * t388 * t143 + t403 - 0.2e1 * t405 - 0.3e1 * t409
            - t412 - t353 * mC * t93 - t416 + 0.4e1 * t419 + t421 - t423;
        double t426 = t10 * t99;
        double t428 = t407 * t264 * t426;
        double t431 = t214 * t2 * t24;
        double t432 = t370 * t431;
        double t433 = 0.2e1 * t432;
        double t434 = t399 * t161;
        double t436 = t382 * t394;
        double t438 = 0.3e1 * t374 * t268;
        double t439 = t399 * t269;
        double t440 = 0.2e1 * t439;
        double t441 = t364 * mC;
        double t442 = t441 * t54;
        double t443 = t344 * t8;
        double t446 = t443 * t214 * t99 * t24;
        double t447 = 0.4e1 * t446;
        double t449 = 0.2e1 * t393 * t42;
        double t450 = t382 * t62;
        double t454 = t344 * t99 * t377 * t214;
        double t455 = 0.4e1 * t454;
        double t456
            = 0.3e1 * t428 + t433 - 0.4e1 * t434 + t436 + t438 - t440 + t442 + t447 - t449 + 0.4e1 * t450 - t455;
        double t460 = 0.3e1 * t1 * t2 * t4;
        double t462 = 0.2e1 * t340 * t394;
        double t464 = t407 * t41 * t426;
        double t465 = 0.3e1 * t464;
        double t466 = t65 * t100;
        double t467 = t399 * t466;
        double t468 = 0.4e1 * t467;
        double t469 = t358 * t55;
        double t471 = t382 * t42;
        double t472 = t245 * t281;
        double t473 = t349 * t134;
        double t474 = 0.2e1 * t473;
        double t475 = t201 * t364;
        double t476 = t475 * t134;
        double t480 = t363 * t364 * t377 * t214;
        double t481 = 0.2e1 * t480;
        double t482
            = -t224 * t99 - t460 + t462 - t465 - t468 + 0.4e1 * t469 + t471 + t472 - t474 - 0.2e1 * t476 + t289 - t481;
        double t488 = 0.1e1 / t4 / t300 / t3;
        double t490 = t142 * t99;
        double t492 = 0.2e1 * t490 * t42;
        double t493 = t358 * t93;
        double t494 = t443 * t431;
        double t495 = 0.2e1 * t494;
        double t499 = t363 * t8 * t9 * t99 * t24;
        double t500 = 0.4e1 * t499;
        double t501 = t443 * t371;
        double t502 = 0.2e1 * t501;
        double t503 = t343 * t99;
        double t505 = t503 * t378 * t1;
        double t506 = 0.4e1 * t505;
        double t507 = t142 * t2;
        double t509 = 0.2e1 * t507 * t152;
        double t510 = t343 * t364;
        double t512 = 0.2e1 * t510 * t366;
        double t513 = t407 * t466;
        double t514 = 0.4e1 * t513;
        double t516 = 0.2e1 * t407 * t418;
        double t517 = t4 * t1;
        double t518 = 0.2e1 * t517;
        double t519 = t268 * t201;
        double t520 = 0.2e1 * t519;
        double t522 = 0.4e1 * t407 * t269;
        double t523 = t99 * mC;
        double t524 = t523 * t54;
        double t525 = 0.2e1 * t524;
        double t527 = 0.2e1 * t353 * t143;
        double t529 = 0.2e1 * t510 * t378;
        double t531 = 0.2e1 * t399 * t394;
        double t532 = t57 * t54;
        double t533 = 0.2e1 * t532;
        double t535 = 0.2e1 * t490 * t264;
        double t536 = t374 * mC;
        double t537 = t536 * t107;
        double t538 = 0.3e1 * t537;
        double t539 = t492 + t493 - t495 + t500 + t502 + t506 + t509 - t512 + t514 + t516 - t518 - t520 - t522 + t525
            - t527 + t529 + t531 - t533 + t535 + t538;
        double t541 = 0.4e1 * t382 * t66;
        double t542 = t523 * t307;
        double t543 = 0.4e1 * t542;
        double t544 = t399 * t383;
        double t545 = 0.8e1 * t544;
        double t546 = 0.2e1 * t268;
        double t547 = t523 * t93;
        double t548 = t201 * t99;
        double t550 = t548 * mC * t93;
        double t552 = t40 * t2;
        double t553 = t552 * t93;
        double t554 = t517 * t201;
        double t555 = 0.2e1 * t554;
        double t556 = t399 * t62;
        double t557 = 0.4e1 * t556;
        double t558 = t507 * t4;
        double t560 = 0.2e1 * t558 * t48;
        double t563 = 0.4e1 * t112 * t106 * t201;
        double t565 = 0.2e1 * t490 * t394;
        double t566 = t44 * t1;
        double t568 = 0.2e1 * t558 * t566;
        double t570 = 0.2e1 * t357 * t355;
        double t571 = t407 * t401;
        double t572 = 0.2e1 * t571;
        double t574 = 0.2e1 * t407 * t161;
        double t575 = t399 * t42;
        double t580 = t503 * t8 * t214 * t24 * t1;
        double t581 = 0.4e1 * t580;
        double t582 = t344 * t379;
        double t583 = 0.2e1 * t582;
        double t585 = 0.2e1 * t510 * t346;
        double t586 = t541 - t543 + t545 + t546 - t547 + 0.2e1 * t550 - t553 + t555 - t557 - t560 - t563 - t565 + t568
            - t570 - t572 - t574 - 0.4e1 * t575 - t581 + t583 + t585;
        double t589 = t235 * x;
        double t590 = t589 * t488;
        double t595 = mC * t79;
        double t596 = t595 * t335;
        double t598 = t548 * t142;
        double t600 = 0.2e1 * t598 * t394;
        double t601 = t598 * t42;
        double t602 = 0.2e1 * t601;
        double t604 = 0.4e1 * t40 * t55;
        double t606 = 0.2e1 * t548 * t355;
        double t609 = -t492 + t493 + t600 + t495 - t500 - t502 - t506 - t602 - t509 + t512 - t514 + t516 - t604 - t606
            + 0.4e1 * t519 - 0.4e1 * t524 + t527 - t529;
        double t610 = 0.2e1 * t542;
        double t612 = t57 * t307;
        double t613 = 0.2e1 * t612;
        double t614 = 0.2e1 * t575;
        double t615 = -t535 + t537 - t541 + t610 - t545 - 0.4e1 * t268 + t547 + t553 - t613 + t557 + t560 + t565 - t568
            + t572 - t574 + t614 + t581 - t583 - t585;
        double t620 = t99 * t321;
        double t621 = t620 * t1;
        double t624 = t364 * t321;
        double t625 = t624 * t1;
        double t627 = 0.2e1 * t625 * t143;
        double t628 = 0.2e1 * t347;
        double t629 = t99 * t99;
        double t630 = t343 * t629;
        double t633 = t226 * t11;
        double t635 = 0.4e1 * t625 * t633;
        double t637 = t142 * t364;
        double t641 = t2 * t321;
        double t643 = t354 * t44;
        double t644 = t641 * t56 * t643;
        double t646 = t352 * mC;
        double t647 = t54 * t321;
        double t648 = t646 * t647;
        double t650 = t2 * t56 * mC;
        double t651 = t650 * t647;
        double t652 = 0.2e1 * t651;
        double t653 = t289 * t321;
        double t655 = t4 * t321;
        double t656 = t655 * t352;
        double t658 = t142 * t20 * t10;
        double t661 = t9 * t321 * t74;
        double t662 = t650 * t661;
        double t663 = t621 * t643;
        double t665 = 0.2e1 * t446;
        double t666
            = 0.2e1 * t644 + t648 - t652 - t412 + t653 - 0.4e1 * t422 - t656 * t658 + t662 - t442 + 0.3e1 * t663 + t665;
        double t672 = 0.2e1 * t454;
        double t673 = t641 * t352;
        double t675 = 0.2e1 * t673 * t143;
        double t677 = t620 * t56;
        double t679 = 0.4e1 * t677 * t143;
        double t680 = t321 * t1;
        double t682 = 0.3e1 * t268 * t680;
        double t683 = t350 * t661;
        double t685 = t352 * t142;
        double t690 = 0.4e1 * t673 * t633;
        double t692 = 0.8e1 * t677 * t633;
        double t694 = t142 * t24 * t10;
        double t696 = t350 * t647;
        double t697 = 0.3e1 * t696;
        double t703 = -t655 * t685 * t24 + t441 * t93 - t690 + t474 - t692 + t656 * t694 - t289 - t697
            - 0.2e1 * t637 * t42 - 0.2e1 * t630 * t346 - 0.4e1 * t480;
        double t708 = t58 * t214;
        double t711 = t58 * t9;
        double t721 = t10 * t321 * t2;
        double t736 = (0.4e1 * t212 * t708 - 0.4e1 * t212 * t711 + 0.2e1 * t40 * t8 * t321 - 0.4e1 * t40 * t11 * t321
                          - 0.4e1 * t111 * t721 - 0.4e1 * t226 * t55 - t26 * t321 + 0.4e1 * t226 * t216 + t26
                          + 0.2e1 * t111 * t641)
            * y * x * t335 * t5 * t595;
        double t738 = 0.2e1 * t392 * t251;
        double t739 = 0.2e1 * t308;
        double t743 = 0.2e1 * t392 * t250 * t10 * t4;
        double t746 = 0.2e1 * t392 * t5 * t42;
        double t747 = 0.2e1 * t295;
        double t748 = t239 - t240 - t244 + t738 - t739 - t743 + t249 + t746 - t277 + t310 + t271 + t254 + t267 + t263
            - t259 - t258 - t312 + t747 - t281 + t278;
        double t752 = t235 * mC * t302;
        double t753 = t748 * t19 * y * t752;
        double t755 = t99 * t24;
        double t761 = t99 * t20 * t4;
        double t764 = t739 - t258 + t271 + t263 - t310 + t267 - t277 + t312 - t259 + t254 + t100 - t100 * t201
            + 0.2e1 * t755 * mC * t293 * t4 - t285 - 0.2e1 * t761 * t242 + t291 + t288 - t747 + t281 - t278;
        double t767 = t764 * t19 * x * t752;
        double t772 = 0.4e1 * t213 * t42;
        double t778 = (0.2e1 * t536 * t41 * t9 + t772 + t136 - 0.2e1 * t536 * t120 + t52 + t219) * y * t235 * t338;
        double t780 = 0.2e1 * t369 * t346;
        double t782 = 0.2e1 * t369 * t366;
        double t783 = 0.4e1 * t494;
        double t784 = 0.2e1 * t499;
        double t785 = 0.4e1 * t501;
        double t786 = 0.2e1 * t505;
        double t790 = 0.2e1 * t354 * t56 * t20 * t10;
        double t791 = t56 * t24;
        double t794 = 0.2e1 * t354 * t791 * t10;
        double t795 = t780 + t600 - t782 - t783 + t784 + t785 + t786 - t602 + t509 - t514 + t516 + t518 - t606 + t520
            + t790 - t794 - t525 + t527 + t533;
        double t796 = t57 * t93;
        double t798 = 0.2e1 * t369 * t378;
        double t799 = 0.2e1 * t580;
        double t800 = 0.4e1 * t582;
        double t802 = 0.2e1 * t354 * t791;
        double t803 = -t537 - t541 + t543 - t545 - t546 - t550 - t553 - t796 - t555 + t557 - t560 + t563 + t568 + t572
            - t574 + t614 + t798 - t799 + t800 + t802;
        double t808 = (t795 + t803) * x * t235 * y * t488;
        double t809 = t475 * t1;
        double t810 = t809 * t633;
        double t812 = 0.2e1 * t361;
        double t817 = 0.4e1 * t384 + t396 - 0.4e1 * t402 + t405 - t409 - t412 - t416 - t423 + t428 + t433 - t436;
        double t820 = t4 * t56;
        double t821 = t820 * t201;
        double t822 = t809 * t143;
        double t827 = 0.3e1 * t473;
        double t831 = t441 * t9 * t201 * t5 * t4;
        double t832 = t646 * t54;
        double t833 = -t465 + 0.8e1 * t467 + 0.3e1 * t469 + t471 + t472 + t820 - t827 - t476 + t831 - t481 + t832;
        double t837 = (-t342 - t348 + t351 + 0.4e1 * t810 + t356 - t812 + t368 - t373 + 0.2e1 * t376 - t381 + t817
                          + t438 - t440 + t447 - t449 - 0.2e1 * t450 - t455 - t821 - t460 + t462 - 0.2e1 * t822 + t833)
            * t235 * t488;
        double t840 = 0.2e1 * t319 * t106 * t321;
        double t842 = t655 * t2;
        double t844 = 0.2e1 * t177 * t24 * t842;
        double t848 = (t328 + t61 + t127 + t330 + t840 - t844) * x * t19 * t596;
        double t849 = t56 * t321;
        double t850 = t849 * t142;
        double t851 = t850 * t418;
        double t853 = t655 * t364;
        double t854 = t853 * t658;
        double t855 = t850 * t161;
        double t858 = t321 * t5 * t4;
        double t860 = 0.3e1 * t644;
        double t861 = t853 * t694;
        double t862 = t854 - t855 + t646 * t9 * t858 - t860 - t648 + t396 - t861 + t652 - t412 - t416 - t653;
        double t864 = 0.2e1 * t663;
        double t867 = t655 * t637 * t24;
        double t868 = t682 - t460 + t462 + t472 + t690 - t474 + t867 + t692 + t289 + t697 - t481;
        double t872 = (-t342 - t627 - t348 + t351 + t635 + t851 - t362 + t368 - t373 - t381 + t862 - t423 + t433 + t662
                          + t442 - t864 + t447 - t449 - t455 - t675 - t679 + t868)
            * t335 * t488;
        double t878 = t620 * t142;
        double t888 = t680 * t142;
        double t891 = t552 * t661 - t492 + t495 - t500 - t502 - 0.4e1 * t680 * mC * t55 - t506 - t509 + t512
            - 0.2e1 * t878 * t42 + 0.2e1 * t850 * t42 + t518 + 0.4e1 * t624 * t142 * t66 + 0.4e1 * t850 * t383 - t525
            - 0.2e1 * t888 * t418 - t529 + t533 - t535;
        double t913 = -0.4e1 * t878 * t153 - 0.2e1 * t624 * t143 - t546 + t547 + 0.2e1 * t842 - 0.2e1 * t680 * t4 + t553
            - 0.4e1 * t849 * t134 + t560 + t565 - t568 + 0.2e1 * t849 * t355 + 0.2e1 * t888 * t161 - 0.2e1 * t850 * t394
            + 0.8e1 * t878 * t65 * t52 + t245 * t858 + t581 - t583 - 0.2e1 * t850 * t62 - t585;
        double t917 = y * t335 * x * (t891 + t913) * t488;
        double t920 = (t310 - t52 + t61) * t9 * t595;
        double t921 = 0.2e1 * t8;
        double t926 = (t74 - t921) * t9 * x * t34 * t79;
        double t947 = 0.2e1 * t123 * t641 * t264;
        double t950 = 0.2e1 * t151 * t177 * t641;
        double t951 = t52 * t641;
        double t954 = 0.4e1 * t69 * t54 * t641;
        double t960 = 0.2e1 * t523 * t9 * t325;
        double t963 = t100 * t321;
        double t966 = 0.2e1 * t761 * t177 * t321;
        double t967 = -t239 + t743 - t746 - t738 + t277 - t263 + t947 - t950 - t951 - t954 - t271 - t254
            - 0.4e1 * t400 * mC * t647 + t960 - 0.2e1 * t641 * t74 - t963 - t966 + t282;
        double t975 = 0.4e1 * t306 * t647 - t239 - t738 - t746 + t743 - t951 - t271 + t954 - t311 + 0.2e1 * t133 * t655
            - t950 - t254 + t947 - t263 + t277 + t960 - t966 - t963;
        double t986 = t343 * t387;
        double t999 = t854 - t855 - t860 - t396 - t861 - 0.2e1 * t986 * t366 + 0.3e1 * t651 + t416 - 0.2e1 * t685 * t394
            - t646 * t93 - 0.4e1 * t432;
        double t1008 = -t683 + 0.2e1 * t986 * t378 + 0.2e1 * t685 * t264 - t472 + t820 + t690 - t827 + t867 + t692
            + 0.2e1 * t696 + t832;
        double t1013 = 0.8e1 * t513;
        double t1016 = 0.2e1 * t475 * t143;
        double t1018 = t780 - t782 - t783 + t784 + t785 + t786 - t602 + t509 + t1013 - t516 + t604 + 0.4e1 * t517
            - t1016 + t790 - t794 - t531 + 0.4e1 * t532 - t537;
        double t1019 = 0.4e1 * t544;
        double t1021 = 0.2e1 * t556;
        double t1022 = 0.4e1 * t571;
        double t1023 = t475 * t142;
        double t1025 = 0.4e1 * t1023 * t66;
        double t1026 = t610 + t1019 - t550 - t553 - t796 - 0.4e1 * t554 - t613 - t1021 - t560 + t568 + t570 - t1022
            + t574 + t614 + t1025 + t798 - t799 + t800 + t802;
        double t1033 = t201 * mC * t4;
        double t1045 = t780 + 0.2e1 * t493 + t600 - t782 - t783 + t784 + t785 + t786 - 0.4e1 * t601 + t509 + t1013
            + t516 + t518 - t1016 - t606 + t520 - t522 + t790 - t794 - t525;
        double t1047 = t533 + t538 + t1019 - t546 + t550 - t553 - t796 - t555 - 0.4e1 * t612 - t1021 - t560 - t563
            + t568 - t1022 - t574 + t1025 + t798 - t799 + t800 + t802;
        double t1055 = t629 * t201;
        double t1058 = -t342 - t348 + t351 - 0.8e1 * t810 - 0.2e1 * t359 - t812 + t368 - t373 + 0.4e1 * t376 - t381
            + 0.2e1 * t1055 * t143;
        double t1065 = t396 + t403 - 0.3e1 * t405 - 0.4e1 * t409 - t412 - t416 + 0.3e1 * t419 - 0.2e1 * t421 - t423
            + 0.4e1 * t428 + t433 - 0.3e1 * t434;
        double t1073 = t438 - 0.3e1 * t439 + t1023 * t394 + t447 - t449 - t455 - t821 - 0.4e1 * t1055 * t142 * t66
            - t460 + t462 + 0.4e1 * t822;
        double t1077
            = -0.2e1 * t464 - t468 + t469 + t1023 * t42 + t472 - t475 * t355 + t820 - t827 + t476 - t831 - t481 + t832;
        double t1084 = (t52 + t281 - t61) * t9 * t595;
        double t1087 = mC * t20;
        double t1095 = exp(-mC * (t25 * t1 + t25 * t2 + t1087 * t1 + t1087 * t2 - 0.2e1 * t7 * t66 - 0.4e1 * t92));
        double t1096 = mC * t1095;
        double t1097 = t4 * t8;
        double t1098 = t1097 * t9;
        double t1128 = t19
            * (-0.2e1 * t273 - 0.2e1 * t272 * t2 - 0.2e1 * t226 * t566 - 0.2e1 * t226 * t45 + 0.2e1 * t272 * t10 * t1
                + 0.2e1 * t272 * t417 - 0.4e1 * t6 * t1098 - t1097)
            / t3;
        double t1130 = t34 * t1095;
        double t1137 = (0.4e1 * t354 * t66 + 0.4e1 * t189 + 0.1e1) * t31;
        double t1140 = t8 * t24;
        double t1143 = t377 * t10;
        double t1146 = t1140 * t10;
        double t1151 = 0.2e1 * t212 * t1140 * t2;
        double t1155 = 0.2e1 * t212 * t8 * t50;
        double t1158 = 0.2e1 * t212 * t1143 * t2;
        double t1161 = t1096 * t79;
        double t1170 = t142 * t8;
        double t1175 = t142 * t377;
        double t1183 = (-0.2e1 * t212 * t1143 - 0.2e1 * t212 * t1140 + 0.2e1 * t212 * t1146 - 0.4e1 * t22 + t74 - t921
                           + 0.2e1 * t1170 * t50 - 0.2e1 * t1170 * t38 - 0.2e1 * t1175 * t417)
            * t9 * x * t1095 * y * t595;

        chrisD[0][0][0][0] = t6 * (0.4e1 * t7 * t11 + t9 - 0.2e1 * t7 * t8);
        chrisD[0][0][0][1] = t33;
        chrisD[0][0][0][2] = t35;
        chrisD[0][0][0][3] = 0.0e0;
        chrisD[0][0][1][0] = t33;
        chrisD[0][0][1][1] = -t80;
        chrisD[0][0][1][2] = t98;
        chrisD[0][0][1][3] = 0.0e0;
        chrisD[0][0][2][0] = t35;
        chrisD[0][0][2][1] = t98;
        chrisD[0][0][2][2] = -t110;
        chrisD[0][0][2][3] = 0.0e0;
        chrisD[0][0][3][0] = 0.0e0;
        chrisD[0][0][3][1] = 0.0e0;
        chrisD[0][0][3][2] = 0.0e0;
        chrisD[0][0][3][3] = 0.0e0;
        chrisD[0][1][0][0] = t33;
        chrisD[0][1][0][1] = -t80;
        chrisD[0][1][0][2] = t98;
        chrisD[0][1][0][3] = 0.0e0;
        chrisD[0][1][1][0] = t118;
        chrisD[0][1][1][1] = t132;
        chrisD[0][1][1][2] = t141;
        chrisD[0][1][1][3] = 0.0e0;
        chrisD[0][1][2][0] = t150;
        chrisD[0][1][2][1] = t160;
        chrisD[0][1][2][2] = t168;
        chrisD[0][1][2][3] = 0.0e0;
        chrisD[0][1][3][0] = 0.0e0;
        chrisD[0][1][3][1] = 0.0e0;
        chrisD[0][1][3][2] = 0.0e0;
        chrisD[0][1][3][3] = 0.0e0;
        chrisD[0][2][0][0] = t35;
        chrisD[0][2][0][1] = t98;
        chrisD[0][2][0][2] = -t110;
        chrisD[0][2][0][3] = 0.0e0;
        chrisD[0][2][1][0] = t150;
        chrisD[0][2][1][1] = t160;
        chrisD[0][2][1][2] = t168;
        chrisD[0][2][1][3] = 0.0e0;
        chrisD[0][2][2][0] = t176;
        chrisD[0][2][2][1] = t185;
        chrisD[0][2][2][2] = t188;
        chrisD[0][2][2][3] = 0.0e0;
        chrisD[0][2][3][0] = 0.0e0;
        chrisD[0][2][3][1] = 0.0e0;
        chrisD[0][2][3][2] = 0.0e0;
        chrisD[0][2][3][3] = 0.0e0;
        chrisD[0][3][0][0] = 0.0e0;
        chrisD[0][3][0][1] = 0.0e0;
        chrisD[0][3][0][2] = 0.0e0;
        chrisD[0][3][0][3] = 0.0e0;
        chrisD[0][3][1][0] = 0.0e0;
        chrisD[0][3][1][1] = 0.0e0;
        chrisD[0][3][1][2] = 0.0e0;
        chrisD[0][3][1][3] = 0.0e0;
        chrisD[0][3][2][0] = 0.0e0;
        chrisD[0][3][2][1] = 0.0e0;
        chrisD[0][3][2][2] = 0.0e0;
        chrisD[0][3][2][3] = 0.0e0;
        chrisD[0][3][3][0] = -t189;
        chrisD[0][3][3][1] = t192;
        chrisD[0][3][3][2] = t195;
        chrisD[0][3][3][3] = 0.0e0;
        chrisD[1][0][0][0] = t33;
        chrisD[1][0][0][1] = -t80;
        chrisD[1][0][0][2] = t98;
        chrisD[1][0][0][3] = 0.0e0;
        chrisD[1][0][1][0] = t118;
        chrisD[1][0][1][1] = t132;
        chrisD[1][0][1][2] = t141;
        chrisD[1][0][1][3] = 0.0e0;
        chrisD[1][0][2][0] = t150;
        chrisD[1][0][2][1] = t160;
        chrisD[1][0][2][2] = t168;
        chrisD[1][0][2][3] = 0.0e0;
        chrisD[1][0][3][0] = 0.0e0;
        chrisD[1][0][3][1] = 0.0e0;
        chrisD[1][0][3][2] = 0.0e0;
        chrisD[1][0][3][3] = 0.0e0;
        chrisD[1][1][0][0] = (0.4e1 * t57 * t11 * t201 - 0.2e1 * t57 * t8 * t201 + 0.4e1 * t112 * t10 * t201 * t2 + t218
                                 - 0.2e1 * t40 * t219 - t223 + t73 * t224 - 0.4e1 * t226 * t54 * t99
                                 + 0.4e1 * t226 * t215 * t99 + t106 * t4)
            * t235 * t237;
        chrisD[1][1][0][1] = -t297 * t19 * t235 * t303;
        chrisD[1][1][0][2] = -t314 * t19 * t235 * t317;
        chrisD[1][1][0][3] = 0.0e0;
        chrisD[1][1][1][0] = (-t324 + t327 - t328 + t330 - 0.2e1 * t331 + t61 + t127) * x * t335 * t338;
        chrisD[1][1][1][1] = -(t386 + t424 + t456 + t482) * t235 * t488;
        chrisD[1][1][1][2] = (t539 + t586) * y * t590;
        chrisD[1][1][1][3] = 0.0e0;
        chrisD[1][1][2][0] = (t328 - t324 + t327 + t330 + t127 + t61) * y * t19 * t596;
        chrisD[1][1][2][1] = -t589 * y * (t609 + t615) * t488;
        chrisD[1][1][2][2] = -(-t342 - t621 * t354 * t47 + t627 - t628 + t351 + 0.2e1 * t630 * t366 - t635 + t362
                                 + 0.4e1 * t367 + 0.2e1 * t637 * t394 + t666 - 0.2e1 * t637 * t264 - 0.2e1 * t630 * t378
                                 - t672 + t675 + t621 * t355 + t679 - t682 + t460 + t462 + t683 + t703)
            * t335 * t488;
        chrisD[1][1][2][3] = 0.0e0;
        chrisD[1][1][3][0] = 0.0e0;
        chrisD[1][1][3][1] = 0.0e0;
        chrisD[1][1][3][2] = 0.0e0;
        chrisD[1][1][3][3] = 0.0e0;
        chrisD[1][2][0][0] = -t736;
        chrisD[1][2][0][1] = t753;
        chrisD[1][2][0][2] = t767;
        chrisD[1][2][0][3] = 0.0e0;
        chrisD[1][2][1][0] = -t778;
        chrisD[1][2][1][1] = -t808;
        chrisD[1][2][1][2] = t837;
        chrisD[1][2][1][3] = 0.0e0;
        chrisD[1][2][2][0] = -t848;
        chrisD[1][2][2][1] = t872;
        chrisD[1][2][2][2] = t917;
        chrisD[1][2][2][3] = 0.0e0;
        chrisD[1][2][3][0] = 0.0e0;
        chrisD[1][2][3][1] = 0.0e0;
        chrisD[1][2][3][2] = 0.0e0;
        chrisD[1][2][3][3] = 0.0e0;
        chrisD[1][3][0][0] = 0.0e0;
        chrisD[1][3][0][1] = 0.0e0;
        chrisD[1][3][0][2] = 0.0e0;
        chrisD[1][3][0][3] = 0.0e0;
        chrisD[1][3][1][0] = 0.0e0;
        chrisD[1][3][1][1] = 0.0e0;
        chrisD[1][3][1][2] = 0.0e0;
        chrisD[1][3][1][3] = 0.0e0;
        chrisD[1][3][2][0] = 0.0e0;
        chrisD[1][3][2][1] = 0.0e0;
        chrisD[1][3][2][2] = 0.0e0;
        chrisD[1][3][2][3] = 0.0e0;
        chrisD[1][3][3][0] = t192;
        chrisD[1][3][3][1] = -t920;
        chrisD[1][3][3][2] = -t926;
        chrisD[1][3][3][3] = 0.0e0;
        chrisD[2][0][0][0] = t35;
        chrisD[2][0][0][1] = t98;
        chrisD[2][0][0][2] = -t110;
        chrisD[2][0][0][3] = 0.0e0;
        chrisD[2][0][1][0] = t150;
        chrisD[2][0][1][1] = t160;
        chrisD[2][0][1][2] = t168;
        chrisD[2][0][1][3] = 0.0e0;
        chrisD[2][0][2][0] = t176;
        chrisD[2][0][2][1] = t185;
        chrisD[2][0][2][2] = t188;
        chrisD[2][0][2][3] = 0.0e0;
        chrisD[2][0][3][0] = 0.0e0;
        chrisD[2][0][3][1] = 0.0e0;
        chrisD[2][0][3][2] = 0.0e0;
        chrisD[2][0][3][3] = 0.0e0;
        chrisD[2][1][0][0] = -t736;
        chrisD[2][1][0][1] = t753;
        chrisD[2][1][0][2] = t767;
        chrisD[2][1][0][3] = 0.0e0;
        chrisD[2][1][1][0] = -t778;
        chrisD[2][1][1][1] = -t808;
        chrisD[2][1][1][2] = t837;
        chrisD[2][1][1][3] = 0.0e0;
        chrisD[2][1][2][0] = -t848;
        chrisD[2][1][2][1] = t872;
        chrisD[2][1][2][2] = t917;
        chrisD[2][1][2][3] = 0.0e0;
        chrisD[2][1][3][0] = 0.0e0;
        chrisD[2][1][3][1] = 0.0e0;
        chrisD[2][1][3][2] = 0.0e0;
        chrisD[2][1][3][3] = 0.0e0;
        chrisD[2][2][0][0]
            = (-0.4e1 * t392 * t711 + 0.4e1 * t392 * t708 + t218 - 0.2e1 * t40 * t331 + 0.4e1 * t112 * t721 - t223
                  + t73 * t4 + t106 * t655 + 0.4e1 * t111 * t426 * t321 - 0.2e1 * t111 * t620)
            * t335 * t237;
        chrisD[2][2][0][1] = t967 * t19 * t335 * t303;
        chrisD[2][2][0][2] = t975 * t19 * t335 * t317;
        chrisD[2][2][0][3] = 0.0e0;
        chrisD[2][2][1][0] = (t52 + t772 + t136 + t844 + t331 - t840) * x * t335 * t338;
        chrisD[2][2][1][1] = (-t820 * t321 - t627 + t628 + 0.2e1 * t685 * t42 + t635 + t851 - t812 + 0.2e1 * t986 * t346
                                 + 0.4e1 * t372 + 0.4e1 * t380 + t999 - t441 * t647 - t662 - t864 - t665 + t449 + t672
                                 - t675 - t679 + t682 - t460 + t1008)
            * t335 * t488;
        chrisD[2][2][1][2] = (t1018 + t1026) * y * t590;
        chrisD[2][2][1][3] = 0.0e0;
        chrisD[2][2][2][0]
            = (t52 - 0.2e1 * t52 * t201 + t772 + t136 + 0.2e1 * t1033 * t371 - 0.2e1 * t1033 * t21 * t2 - t219) * y
            * t235 * t338;
        chrisD[2][2][2][1] = t589 * y * (t1045 + t1047) * t488;
        chrisD[2][2][2][2] = -(t1058 + t1065 + t1073 + t1077) * t235 * t488;
        chrisD[2][2][2][3] = 0.0e0;
        chrisD[2][2][3][0] = 0.0e0;
        chrisD[2][2][3][1] = 0.0e0;
        chrisD[2][2][3][2] = 0.0e0;
        chrisD[2][2][3][3] = 0.0e0;
        chrisD[2][3][0][0] = 0.0e0;
        chrisD[2][3][0][1] = 0.0e0;
        chrisD[2][3][0][2] = 0.0e0;
        chrisD[2][3][0][3] = 0.0e0;
        chrisD[2][3][1][0] = 0.0e0;
        chrisD[2][3][1][1] = 0.0e0;
        chrisD[2][3][1][2] = 0.0e0;
        chrisD[2][3][1][3] = 0.0e0;
        chrisD[2][3][2][0] = 0.0e0;
        chrisD[2][3][2][1] = 0.0e0;
        chrisD[2][3][2][2] = 0.0e0;
        chrisD[2][3][2][3] = 0.0e0;
        chrisD[2][3][3][0] = t195;
        chrisD[2][3][3][1] = -t926;
        chrisD[2][3][3][2] = -t1084;
        chrisD[2][3][3][3] = 0.0e0;
        chrisD[3][0][0][0] = 0.0e0;
        chrisD[3][0][0][1] = 0.0e0;
        chrisD[3][0][0][2] = 0.0e0;
        chrisD[3][0][0][3] = 0.0e0;
        chrisD[3][0][1][0] = 0.0e0;
        chrisD[3][0][1][1] = 0.0e0;
        chrisD[3][0][1][2] = 0.0e0;
        chrisD[3][0][1][3] = 0.0e0;
        chrisD[3][0][2][0] = 0.0e0;
        chrisD[3][0][2][1] = 0.0e0;
        chrisD[3][0][2][2] = 0.0e0;
        chrisD[3][0][2][3] = 0.0e0;
        chrisD[3][0][3][0] = -t189;
        chrisD[3][0][3][1] = t192;
        chrisD[3][0][3][2] = t195;
        chrisD[3][0][3][3] = 0.0e0;
        chrisD[3][1][0][0] = 0.0e0;
        chrisD[3][1][0][1] = 0.0e0;
        chrisD[3][1][0][2] = 0.0e0;
        chrisD[3][1][0][3] = 0.0e0;
        chrisD[3][1][1][0] = 0.0e0;
        chrisD[3][1][1][1] = 0.0e0;
        chrisD[3][1][1][2] = 0.0e0;
        chrisD[3][1][1][3] = 0.0e0;
        chrisD[3][1][2][0] = 0.0e0;
        chrisD[3][1][2][1] = 0.0e0;
        chrisD[3][1][2][2] = 0.0e0;
        chrisD[3][1][2][3] = 0.0e0;
        chrisD[3][1][3][0] = t192;
        chrisD[3][1][3][1] = -t920;
        chrisD[3][1][3][2] = -t926;
        chrisD[3][1][3][3] = 0.0e0;
        chrisD[3][2][0][0] = 0.0e0;
        chrisD[3][2][0][1] = 0.0e0;
        chrisD[3][2][0][2] = 0.0e0;
        chrisD[3][2][0][3] = 0.0e0;
        chrisD[3][2][1][0] = 0.0e0;
        chrisD[3][2][1][1] = 0.0e0;
        chrisD[3][2][1][2] = 0.0e0;
        chrisD[3][2][1][3] = 0.0e0;
        chrisD[3][2][2][0] = 0.0e0;
        chrisD[3][2][2][1] = 0.0e0;
        chrisD[3][2][2][2] = 0.0e0;
        chrisD[3][2][2][3] = 0.0e0;
        chrisD[3][2][3][0] = t195;
        chrisD[3][2][3][1] = -t926;
        chrisD[3][2][3][2] = -t1084;
        chrisD[3][2][3][3] = 0.0e0;
        chrisD[3][3][0][0]
            = -t1096 * t5 * (-0.4e1 * t226 * t1098 + 0.4e1 * t226 * t1097 * t214 - 0.4e1 * t6 + 0.4e1 * t6 * t10 + t9);
        chrisD[3][3][0][1] = -t18 * t1095 * t1128;
        chrisD[3][3][0][2] = -t1130 * t1128;
        chrisD[3][3][0][3] = 0.0e0;
        chrisD[3][3][1][0] = -t1096 * t8 * t130 * t1137;
        chrisD[3][3][1][1] = (-0.2e1 * t392 * t1140 - 0.2e1 * t392 * t1143 + 0.2e1 * t392 * t1146 - t52 - t1151
                                 - 0.4e1 * t124 + t310 + t1155 - t1158 + t61)
            * t9 * t1161;
        chrisD[3][3][1][2] = t1183;
        chrisD[3][3][1][3] = 0.0e0;
        chrisD[3][3][2][0] = -t1130 * t8 * t19 * t1137;
        chrisD[3][3][2][1] = t1183;
        chrisD[3][3][2][2] = (-t1158 + t52 + t1155 - t1151 - 0.2e1 * t1175 * t426 - 0.2e1 * t1170 * t755 - 0.4e1 * t181
                                 + 0.2e1 * t1170 * t47 * t99 + t281 - t61)
            * t9 * t1161;
        chrisD[3][3][2][3] = 0.0e0;
        chrisD[3][3][3][0] = 0.0e0;
        chrisD[3][3][3][1] = 0.0e0;
        chrisD[3][3][3][2] = 0.0e0;
        chrisD[3][3][3][3] = 0.0e0;
    }

    return true;
}

void MetricBesselGravWaveCart::localToCoord(const double* pos, const double* ldir, double* dir, enum_nat_tetrad_type)
{
    calcLTcoeffs(pos);

    dir[0] = ldir[0] / ltT;
    dir[1] = ldir[1] / ltA;
    dir[2] = ldir[2] / ltC - ldir[1] * ltB / (ltA * ltC);
    dir[3] = ldir[3] / ltZ;
}

void MetricBesselGravWaveCart::coordToLocal(const double* pos, const double* cdir, double* ldir, enum_nat_tetrad_type)
{
    calcLTcoeffs(pos);

    ldir[0] = cdir[0] * ltT;
    ldir[1] = cdir[1] * ltA;
    ldir[2] = cdir[2] * ltC + cdir[1] * ltB;
    ldir[3] = cdir[3] * ltZ;
}

bool MetricBesselGravWaveCart::breakCondition(const double*)
{
    bool br = false;

    return br;
}

double MetricBesselGravWaveCart::testConstraint(const double Y[], const double kappa)
{
    double cm = 1.0 / mSpeedOfLight;

    double t = Y[0];
    double x = Y[1];
    double y = Y[2];
    double z = Y[3];
    double pos[4] = { t, x, y, z };

    calculateMetric(pos);

    // Scale the directions with the speed of light before doubling them !!
    double dt = Y[4];
    double dx = Y[5] * cm;
    double dy = Y[6] * cm;
    double dz = Y[7] * cm;

    double sum = -kappa;
    sum += g_compts[0][0] * dt * dt + g_compts[1][1] * dx * dx + g_compts[2][2] * dy * dy
        + 2.0 * g_compts[2][1] * dx * dy + g_compts[3][3] * dz * dz;

    return sum;
}

bool MetricBesselGravWaveCart::setParam(const char* pName, double val)
{
    Metric::setParam(pName, val);

    if (strcmp(pName, "c") == 0) {
        mC = val;
    }

    return true;
}

bool MetricBesselGravWaveCart::transToTwoPlusOne(vec4 p, vec4& cp)
{
    vec4 tp;
    TransCoordinates::toCartesianCoord(mCoordType, p, tp);
    cp = vec4(tp[0], tp[1], tp[2], tp[0]);
    return true;
}

bool MetricBesselGravWaveCart::report(const vec4, const vec4, char*& text)
{
    std::stringstream ss;
    ss << "Report for gravitational Bessel wave metric\n\tcoordinate : (t,x,y,z)\n";
    ss << "---------------------------------------------------------------\n";
    ss << "  physical units ................. no\n";
    ss.precision(DEF_FIXED_REPORT_PRECISION);
    ss.setf(std::ios::fixed);
    ss << "  \"Amplitude\" C ............. = " << mC << std::endl;

    text = new char[ss.str().length() + 2];
    return CopyString(ss.str().c_str(), text);
}

// ********************************* protected methods *****************************

void MetricBesselGravWaveCart::setStandardValues()
{
    mInitPos[0] = 0.0;
    mInitPos[1] = 0.0;
    mInitPos[2] = 0.0;
    mInitPos[3] = 0.0;
    mInitDir[0] = 1.0;
    mInitDir[1] = 0.0;
    mInitDir[2] = 0.0;

    mCoordNames[0] = std::string("t");
    mCoordNames[1] = std::string("x");
    mCoordNames[2] = std::string("y");
    mCoordNames[3] = std::string("z");
}

void MetricBesselGravWaveCart::calcLTcoeffs(const double* pos)
{
    double t = pos[0];
    double x = pos[1];
    double y = pos[2];

    if (x == 0.0 && y == 0.0) {
        double eCos = exp(mC * cos(t));

        ltA = 1.0 / eCos;
        ltB = 0.0;
        ltC = 1.0 / eCos;
        ltD = 1.0 / eCos;
        ltE = 0.0;
        ltF = 1.0 / eCos;
        ltT = 1.0 / eCos;
        ltZ = eCos;
    }
    else {
        double xx = x * x;
        double yy = y * y;
        double rho = sqrt(xx + yy);
        double ct = cos(t);
        double j0 = gsl_sf_bessel_J0(rho);
        double j1 = gsl_sf_bessel_J1(rho);
        double U = mC * j0 * ct;
        double K = 0.5 * mC * mC * rho * (rho * (j0 * j0 + j1 * j1) - 2.0 * j0 * j1 * ct * ct);
        double eU = exp(U);
        double e2K = exp(2.0 * K);
        double eKU = exp(K - U);

        ltA = eKU * rho * sqrt(1.0 / (e2K * yy + xx));
        ltB = x * y * (e2K - 1) / (eU * rho * sqrt(e2K * yy + xx));
        ltC = sqrt(e2K * yy + xx) / (eU * rho);
        ltD = sqrt(e2K * xx + yy) / (eU * rho);
        ltE = x * y * (e2K - 1) / (eU * rho * sqrt(e2K * xx + yy));
        ltF = eKU * rho * sqrt(1.0 / (e2K * xx + yy));
        ltT = eKU;
        ltZ = eU;
    }
}

} // end namespace m4d
