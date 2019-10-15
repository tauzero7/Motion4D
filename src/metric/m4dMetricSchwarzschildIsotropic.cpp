/**
 * @file    m4dMetricSchwarzschildIsotropic.cpp
 * @author  Thomas Mueller
 *
 *  This file is part of libMotion4D.
 */
#include "m4dMetricSchwarzschildIsotropic.h"
#include "extra/m4dUtilities.h"

namespace m4d {

#define eps 1.0e-6

MetricSchwarzschildIsotropic::MetricSchwarzschildIsotropic(double mass)
{
    mMetricName = "SchwarzschildIsotropic";
    mCoordType = enum_coordinate_cartesian;

    mPhysicalUnits = enum_physical_constants_geom;
    mSpeedOfLight = 1.0;
    mGravConstant = 1.0;

    addParam("mass", mass);
    mMass = mass;
    rho_s = 0.5 * mGravConstant * mMass / (mSpeedOfLight * mSpeedOfLight);
    calc_orbits();

    //  Only a static tetrad is defined
    mLocTeds.push_back(enum_nat_tetrad_static);

    // Have embedding available
    mDrawTypes.push_back(enum_draw_embedding);

    // parameters for the embedding diagram
    if (!mEmbParam.empty()) {
        mEmbParam.clear();
    }
    mHaveEmbedding = true;

    mEmb_rmin = rho_s;
    mEmb_rmax = 30.0 * rho_s;
    mEmb_r_num = 20;
    mEmb_phi_num = 40;
    mEmb_rstep = (mEmb_rmax - mEmb_rmin) / static_cast<double>(mEmb_r_num);
    mEmb_phistep = 2.0 * M_PI / mEmb_phi_num;
    addEmbeddingParam("emb_rmin", mEmb_rmin);
    addEmbeddingParam("emb_rmax", mEmb_rmax);
    addEmbeddingParam("emb_r_num", 20.0);
    addEmbeddingParam("emb_phi_num", 40.0);

    setStandardValues();
}

MetricSchwarzschildIsotropic::~MetricSchwarzschildIsotropic() {}

bool MetricSchwarzschildIsotropic::calculateMetric(const double* pos)
{
    double rho = calc_rho(pos);

    double c = mSpeedOfLight;

    double t1 = rho; // rho(x,y,z);
    double t3 = rho_s / t1;
    double t5 = pow(1.0 + t3, 2.0);
    double t7 = c * c;
    double t8 = 1 / t5 * t7;
    double t11 = rho_s * rho_s;
    double t12 = t1 * t1;
    double t14 = t11 / t12;
    double t24 = t11 * t11;
    double t25 = t12 * t12;
    double t28 = 1.0 + 4.0 * t3 + 6.0 * t14 + 4.0 * t11 * rho_s / t12 / t1 + t24 / t25;

    g_compts[0][0] = -t8 + 2.0 * t8 * t3 - t8 * t14;
    g_compts[0][1] = 0.0;
    g_compts[0][2] = 0.0;
    g_compts[0][3] = 0.0;
    g_compts[1][0] = 0.0;
    g_compts[1][1] = t28;
    g_compts[1][2] = 0.0;
    g_compts[1][3] = 0.0;
    g_compts[2][0] = 0.0;
    g_compts[2][1] = 0.0;
    g_compts[2][2] = t28;
    g_compts[2][3] = 0.0;
    g_compts[3][0] = 0.0;
    g_compts[3][1] = 0.0;
    g_compts[3][2] = 0.0;
    g_compts[3][3] = t28;

    return true;
}

bool MetricSchwarzschildIsotropic::calculateChristoffels(const double* pos)
{
    double drdx, drdy, drdz;
    calc_drho(pos, drdx, drdy, drdz);

    double c = mSpeedOfLight;

    double t1 = rho; // rho(x,y,z);
    double t2 = t1 * t1;
    double t3 = t2 * t2;
    double t7 = rho_s * rho_s;
    double t13 = t7 * t7;
    double t17 = c * c;
    double t18 = t3 / (t3 + 4.0 * rho_s * t2 * t1 + 6.0 * t7 * t2 + 4.0 * t7 * rho_s * t1 + t13) * t17;
    double t19 = drdx; // diff(rho(x,y,z),x);
    double t20 = rho_s * t19;
    double t21 = t1 - rho_s;
    double t22 = t1 + rho_s;
    double t23 = t22 * t22;
    double t26 = t21 / t23 / t22;
    double t30 = drdy; // diff(rho(x,y,z),y);
    double t31 = rho_s * t30;
    double t35 = drdz; // diff(rho(x,y,z),z);
    double t36 = rho_s * t35;
    double t40 = 1 / t22;
    double t42 = t40 / t21;
    double t44 = 2.0 * t20 * t42;
    double t46 = 2.0 * t31 * t42;
    double t48 = 2.0 * t36 * t42;
    double t50 = 1 / t1 * t40;
    double t52 = 2.0 * t20 * t50;
    double t54 = 2.0 * t31 * t50;
    double t56 = 2.0 * t36 * t50;

    christoffel[0][0][0] = 0.0;
    christoffel[0][0][1] = 2.0 * t18 * t20 * t26;
    christoffel[0][0][2] = 2.0 * t18 * t31 * t26;
    christoffel[0][0][3] = 2.0 * t18 * t36 * t26;
    christoffel[0][1][0] = t44;
    christoffel[0][1][1] = 0.0;
    christoffel[0][1][2] = 0.0;
    christoffel[0][1][3] = 0.0;
    christoffel[0][2][0] = t46;
    christoffel[0][2][1] = 0.0;
    christoffel[0][2][2] = 0.0;
    christoffel[0][2][3] = 0.0;
    christoffel[0][3][0] = t48;
    christoffel[0][3][1] = 0.0;
    christoffel[0][3][2] = 0.0;
    christoffel[0][3][3] = 0.0;
    christoffel[1][0][0] = t44;
    christoffel[1][0][1] = 0.0;
    christoffel[1][0][2] = 0.0;
    christoffel[1][0][3] = 0.0;
    christoffel[1][1][0] = 0.0;
    christoffel[1][1][1] = -t52;
    christoffel[1][1][2] = t54;
    christoffel[1][1][3] = t56;
    christoffel[1][2][0] = 0.0;
    christoffel[1][2][1] = -t54;
    christoffel[1][2][2] = -t52;
    christoffel[1][2][3] = 0.0;
    christoffel[1][3][0] = 0.0;
    christoffel[1][3][1] = -t56;
    christoffel[1][3][2] = 0.0;
    christoffel[1][3][3] = -t52;
    christoffel[2][0][0] = t46;
    christoffel[2][0][1] = 0.0;
    christoffel[2][0][2] = 0.0;
    christoffel[2][0][3] = 0.0;
    christoffel[2][1][0] = 0.0;
    christoffel[2][1][1] = -t54;
    christoffel[2][1][2] = -t52;
    christoffel[2][1][3] = 0.0;
    christoffel[2][2][0] = 0.0;
    christoffel[2][2][1] = t52;
    christoffel[2][2][2] = -t54;
    christoffel[2][2][3] = t56;
    christoffel[2][3][0] = 0.0;
    christoffel[2][3][1] = 0.0;
    christoffel[2][3][2] = -t56;
    christoffel[2][3][3] = -t54;
    christoffel[3][0][0] = t48;
    christoffel[3][0][1] = 0.0;
    christoffel[3][0][2] = 0.0;
    christoffel[3][0][3] = 0.0;
    christoffel[3][1][0] = 0.0;
    christoffel[3][1][1] = -t56;
    christoffel[3][1][2] = 0.0;
    christoffel[3][1][3] = -t52;
    christoffel[3][2][0] = 0.0;
    christoffel[3][2][1] = 0.0;
    christoffel[3][2][2] = -t56;
    christoffel[3][2][3] = -t54;
    christoffel[3][3][0] = 0.0;
    christoffel[3][3][1] = t52;
    christoffel[3][3][2] = t54;
    christoffel[3][3][3] = -t56;

    return true;
}

bool MetricSchwarzschildIsotropic::calculateChrisD(const double* pos)
{
    double drdx, drdy, drdz;
    calc_drho(pos, drdx, drdy, drdz);

    double drdxdx, drdxdy, drdxdz, drdydy, drdydz, drdzdz;
    calc_d2rho(pos, drdxdx, drdxdy, drdxdz, drdydy, drdydz, drdzdz);

    double c = mSpeedOfLight;

    double t1 = rho_s * rho_s;
    double t2 = drdx; // diff(rho(x,y,z),x);
    double t3 = t2 * t2;
    double t6 = rho; // rho(x,y,z);
    double t7 = t1 * t6;
    double t8 = drdxdx; // diff(diff(rho(x,y,z),x),x);
    double t10 = t6 * t3;
    double t13 = t6 * t6;
    double t14 = t13 * t6;
    double t20 = c * c;
    double t22 = t13 * t13;
    double t30 = t1 * t1;
    double t35 = pow(t6 + rho_s, 2.0);
    double t36 = t35 * t35;
    double t38 = rho_s / (t22 + 4.0 * rho_s * t14 + 6.0 * t1 * t13 + 4.0 * t1 * rho_s * t6 + t30) / t36;
    double t41 = drdxdy; // diff(diff(rho(x,y,z),x),y);
    double t43 = t2 * t1;
    double t44 = drdy; // diff(rho(x,y,z),y);
    double t47 = t6 * t2;
    double t48 = t44 * rho_s;
    double t59
        = 2.0 * (t7 * t41 + 4.0 * t43 * t44 - 8.0 * t47 * t48 + 2.0 * t2 * t44 * t13 - t14 * t41) * t14 * t20 * t38;
    double t60 = drdxdz; // diff(diff(rho(x,y,z),x),z);
    double t62 = drdz; // diff(rho(x,y,z),z);
    double t65 = t62 * rho_s;
    double t76
        = 2.0 * (t7 * t60 + 4.0 * t43 * t62 - 8.0 * t47 * t65 + 2.0 * t2 * t62 * t13 - t14 * t60) * t14 * t20 * t38;
    double t77 = drdydy; // diff(diff(rho(x,y,z),y),y);
    double t79 = t44 * t44;
    double t82 = t6 * t79;
    double t93 = drdydz; // diff(diff(rho(x,y,z),y),z);
    double t98 = t6 * t44;
    double t109 = 2.0 * (t7 * t93 + 4.0 * t1 * t44 * t62 - 8.0 * t98 * t65 + 2.0 * t44 * t62 * t13 - t14 * t93) * t14
        * t20 * t38;
    double t110 = drdzdz; // diff(diff(rho(x,y,z),z),z);
    double t112 = t62 * t62;
    double t115 = t6 * t112;
    double t126 = t13 * t8;
    double t128 = 2.0 * t10;
    double t131 = 1 / t35;
    double t133 = pow(-t6 + rho_s, 2.0);
    double t135 = t131 / t133;
    double t137 = 2.0 * rho_s * (-t126 + t8 * t1 + t128) * t135;
    double t138 = t13 * t41;
    double t141 = 2.0 * t47 * t44;
    double t145 = 2.0 * rho_s * (-t138 + t41 * t1 + t141) * t135;
    double t146 = t13 * t60;
    double t149 = 2.0 * t47 * t62;
    double t153 = 2.0 * rho_s * (-t146 + t60 * t1 + t149) * t135;
    double t154 = t13 * t77;
    double t156 = 2.0 * t82;
    double t160 = 2.0 * rho_s * (-t154 + t77 * t1 + t156) * t135;
    double t161 = t13 * t93;
    double t164 = 2.0 * t98 * t62;
    double t168 = 2.0 * rho_s * (-t161 + t93 * t1 + t164) * t135;
    double t169 = t13 * t110;
    double t171 = 2.0 * t115;
    double t175 = 2.0 * rho_s * (-t169 + t110 * t1 + t171) * t135;
    double t182 = 1 / t13 * t131;
    double t184 = 2.0 * rho_s * (t126 + t6 * t8 * rho_s - t128 - rho_s * t3) * t182;
    double t187 = rho_s * t2;
    double t192 = 2.0 * rho_s * (t138 + t6 * t41 * rho_s - t141 - t187 * t44) * t182;
    double t199 = 2.0 * rho_s * (-t146 - t6 * t60 * rho_s + t149 + t187 * t62) * t182;
    double t206 = 2.0 * rho_s * (t154 + t6 * t77 * rho_s - t156 - rho_s * t79) * t182;
    double t213 = 2.0 * rho_s * (-t161 - t6 * t93 * rho_s + t164 + t48 * t62) * t182;
    double t220 = 2.0 * rho_s * (t169 + t6 * t110 * rho_s - t171 - rho_s * t112) * t182;

    chrisD[0][0][0][0] = 0.0;
    chrisD[0][0][0][1] = 0.0;
    chrisD[0][0][0][2] = 0.0;
    chrisD[0][0][0][3] = 0.0;
    chrisD[0][0][1][0] = 0.0;
    chrisD[0][0][1][1]
        = -2.0 * (4.0 * t1 * t3 + t7 * t8 - 8.0 * t10 * rho_s - t8 * t14 + 2.0 * t3 * t13) * t14 * t20 * t38;
    chrisD[0][0][1][2] = -t59;
    chrisD[0][0][1][3] = -t76;
    chrisD[0][0][2][0] = 0.0;
    chrisD[0][0][2][1] = -t59;
    chrisD[0][0][2][2]
        = -2.0 * (t7 * t77 + 4.0 * t1 * t79 - 8.0 * t82 * rho_s + 2.0 * t79 * t13 - t14 * t77) * t14 * t20 * t38;
    chrisD[0][0][2][3] = -t109;
    chrisD[0][0][3][0] = 0.0;
    chrisD[0][0][3][1] = -t76;
    chrisD[0][0][3][2] = -t109;
    chrisD[0][0][3][3]
        = -2.0 * (t7 * t110 + 4.0 * t1 * t112 - 8.0 * t115 * rho_s + 2.0 * t112 * t13 - t14 * t110) * t14 * t20 * t38;
    chrisD[0][1][0][0] = 0.0;
    chrisD[0][1][0][1] = -t137;
    chrisD[0][1][0][2] = -t145;
    chrisD[0][1][0][3] = -t153;
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
    chrisD[0][2][0][1] = -t145;
    chrisD[0][2][0][2] = -t160;
    chrisD[0][2][0][3] = -t168;
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
    chrisD[0][3][0][1] = -t153;
    chrisD[0][3][0][2] = -t168;
    chrisD[0][3][0][3] = -t175;
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
    chrisD[1][0][0][1] = -t137;
    chrisD[1][0][0][2] = -t145;
    chrisD[1][0][0][3] = -t153;
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
    chrisD[1][1][1][1] = -t184;
    chrisD[1][1][1][2] = -t192;
    chrisD[1][1][1][3] = t199;
    chrisD[1][1][2][0] = 0.0;
    chrisD[1][1][2][1] = t192;
    chrisD[1][1][2][2] = t206;
    chrisD[1][1][2][3] = -t213;
    chrisD[1][1][3][0] = 0.0;
    chrisD[1][1][3][1] = -t199;
    chrisD[1][1][3][2] = -t213;
    chrisD[1][1][3][3] = t220;
    chrisD[1][2][0][0] = 0.0;
    chrisD[1][2][0][1] = 0.0;
    chrisD[1][2][0][2] = 0.0;
    chrisD[1][2][0][3] = 0.0;
    chrisD[1][2][1][0] = 0.0;
    chrisD[1][2][1][1] = -t192;
    chrisD[1][2][1][2] = -t206;
    chrisD[1][2][1][3] = t213;
    chrisD[1][2][2][0] = 0.0;
    chrisD[1][2][2][1] = -t184;
    chrisD[1][2][2][2] = -t192;
    chrisD[1][2][2][3] = t199;
    chrisD[1][2][3][0] = 0.0;
    chrisD[1][2][3][1] = 0.0;
    chrisD[1][2][3][2] = 0.0;
    chrisD[1][2][3][3] = 0.0;
    chrisD[1][3][0][0] = 0.0;
    chrisD[1][3][0][1] = 0.0;
    chrisD[1][3][0][2] = 0.0;
    chrisD[1][3][0][3] = 0.0;
    chrisD[1][3][1][0] = 0.0;
    chrisD[1][3][1][1] = t199;
    chrisD[1][3][1][2] = t213;
    chrisD[1][3][1][3] = -t220;
    chrisD[1][3][2][0] = 0.0;
    chrisD[1][3][2][1] = 0.0;
    chrisD[1][3][2][2] = 0.0;
    chrisD[1][3][2][3] = 0.0;
    chrisD[1][3][3][0] = 0.0;
    chrisD[1][3][3][1] = -t184;
    chrisD[1][3][3][2] = -t192;
    chrisD[1][3][3][3] = t199;
    chrisD[2][0][0][0] = 0.0;
    chrisD[2][0][0][1] = -t145;
    chrisD[2][0][0][2] = -t160;
    chrisD[2][0][0][3] = -t168;
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
    chrisD[2][1][1][1] = -t192;
    chrisD[2][1][1][2] = -t206;
    chrisD[2][1][1][3] = t213;
    chrisD[2][1][2][0] = 0.0;
    chrisD[2][1][2][1] = -t184;
    chrisD[2][1][2][2] = -t192;
    chrisD[2][1][2][3] = t199;
    chrisD[2][1][3][0] = 0.0;
    chrisD[2][1][3][1] = 0.0;
    chrisD[2][1][3][2] = 0.0;
    chrisD[2][1][3][3] = 0.0;
    chrisD[2][2][0][0] = 0.0;
    chrisD[2][2][0][1] = 0.0;
    chrisD[2][2][0][2] = 0.0;
    chrisD[2][2][0][3] = 0.0;
    chrisD[2][2][1][0] = 0.0;
    chrisD[2][2][1][1] = t184;
    chrisD[2][2][1][2] = t192;
    chrisD[2][2][1][3] = -t199;
    chrisD[2][2][2][0] = 0.0;
    chrisD[2][2][2][1] = -t192;
    chrisD[2][2][2][2] = -t206;
    chrisD[2][2][2][3] = t213;
    chrisD[2][2][3][0] = 0.0;
    chrisD[2][2][3][1] = -t199;
    chrisD[2][2][3][2] = -t213;
    chrisD[2][2][3][3] = t220;
    chrisD[2][3][0][0] = 0.0;
    chrisD[2][3][0][1] = 0.0;
    chrisD[2][3][0][2] = 0.0;
    chrisD[2][3][0][3] = 0.0;
    chrisD[2][3][1][0] = 0.0;
    chrisD[2][3][1][1] = 0.0;
    chrisD[2][3][1][2] = 0.0;
    chrisD[2][3][1][3] = 0.0;
    chrisD[2][3][2][0] = 0.0;
    chrisD[2][3][2][1] = t199;
    chrisD[2][3][2][2] = t213;
    chrisD[2][3][2][3] = -t220;
    chrisD[2][3][3][0] = 0.0;
    chrisD[2][3][3][1] = -t192;
    chrisD[2][3][3][2] = -t206;
    chrisD[2][3][3][3] = t213;
    chrisD[3][0][0][0] = 0.0;
    chrisD[3][0][0][1] = -t153;
    chrisD[3][0][0][2] = -t168;
    chrisD[3][0][0][3] = -t175;
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
    chrisD[3][1][1][1] = t199;
    chrisD[3][1][1][2] = t213;
    chrisD[3][1][1][3] = -t220;
    chrisD[3][1][2][0] = 0.0;
    chrisD[3][1][2][1] = 0.0;
    chrisD[3][1][2][2] = 0.0;
    chrisD[3][1][2][3] = 0.0;
    chrisD[3][1][3][0] = 0.0;
    chrisD[3][1][3][1] = -t184;
    chrisD[3][1][3][2] = -t192;
    chrisD[3][1][3][3] = t199;
    chrisD[3][2][0][0] = 0.0;
    chrisD[3][2][0][1] = 0.0;
    chrisD[3][2][0][2] = 0.0;
    chrisD[3][2][0][3] = 0.0;
    chrisD[3][2][1][0] = 0.0;
    chrisD[3][2][1][1] = 0.0;
    chrisD[3][2][1][2] = 0.0;
    chrisD[3][2][1][3] = 0.0;
    chrisD[3][2][2][0] = 0.0;
    chrisD[3][2][2][1] = t199;
    chrisD[3][2][2][2] = t213;
    chrisD[3][2][2][3] = -t220;
    chrisD[3][2][3][0] = 0.0;
    chrisD[3][2][3][1] = -t192;
    chrisD[3][2][3][2] = -t206;
    chrisD[3][2][3][3] = t213;
    chrisD[3][3][0][0] = 0.0;
    chrisD[3][3][0][1] = 0.0;
    chrisD[3][3][0][2] = 0.0;
    chrisD[3][3][0][3] = 0.0;
    chrisD[3][3][1][0] = 0.0;
    chrisD[3][3][1][1] = t184;
    chrisD[3][3][1][2] = t192;
    chrisD[3][3][1][3] = -t199;
    chrisD[3][3][2][0] = 0.0;
    chrisD[3][3][2][1] = t192;
    chrisD[3][3][2][2] = t206;
    chrisD[3][3][2][3] = -t213;
    chrisD[3][3][3][0] = 0.0;
    chrisD[3][3][3][1] = t199;
    chrisD[3][3][3][2] = t213;
    chrisD[3][3][3][3] = -t220;

    return true;
}

void MetricSchwarzschildIsotropic::localToCoord(
    const double* pos, const double* ldir, double* dir, enum_nat_tetrad_type)
{
    double rho = calc_rho(pos);

    double B = (1.0 + rho_s / rho);
    double A = (1.0 - rho_s / rho) / B;
    double edb2 = 1.0 / (B * B);

    dir[0] = ldir[0] / A / mSpeedOfLight;
    dir[1] = ldir[1] * edb2;
    dir[2] = ldir[2] * edb2;
    dir[3] = ldir[3] * edb2;
}

void MetricSchwarzschildIsotropic::coordToLocal(
    const double* pos, const double* cdir, double* ldir, enum_nat_tetrad_type)
{
    double rho = calc_rho(pos);

    double B = (1.0 + rho_s / rho);
    double A = (1.0 - rho_s / rho) / B;
    double b2 = B * B;

    ldir[0] = cdir[0] * A * mSpeedOfLight;
    ldir[1] = cdir[1] * b2;
    ldir[2] = cdir[2] * b2;
    ldir[3] = cdir[3] * b2;
}

bool MetricSchwarzschildIsotropic::breakCondition(const double* pos)
{
    bool br = false;

    double rho = calc_rho(pos);
    if (1.0 - rho_s / rho < eps) {
        br = true;
    }
    return br;
}

double MetricSchwarzschildIsotropic::testConstraint(const double y[], const double kappa)
{
    double rho = calc_rho(y);
    double cm = 1.0 / mSpeedOfLight;

    double dt = y[4];
    double dx = y[5] * cm;
    double dy = y[6] * cm;
    double dz = y[7] * cm;

    double B = (1.0 + rho_s / rho);
    double A = (1.0 - rho_s / rho) / B;
    double b4 = B * B * B * B;

    double sum = -kappa;
    sum += -A * A * dt * dt + b4 * (dx * dx + dy * dy + dz * dz);
    return sum;
}

bool MetricSchwarzschildIsotropic::calcProduct(const double* pos, const double* u, const double* v, double& prod, bool)
{
    prod = 0.0;
    if (breakCondition(pos)) {
        return false;
    }

    double rho = calc_rho(pos);
    double B = (1.0 + rho_s / rho);
    double A = (1.0 - rho_s / rho) / B;
    double b4 = B * B * B * B;

    prod = -mSpeedOfLight * mSpeedOfLight * A * A * u[0] * v[0] + b4 * (u[1] * v[1] + u[2] * v[2] + u[3] * v[3]);
    return true;
}

bool MetricSchwarzschildIsotropic::setParam(const char* pName, double val)
{
    if (Metric::setParam(pName, val)) {
        mMass = val;
        rho_s = 0.5 * mGravConstant * mMass / (mSpeedOfLight * mSpeedOfLight);
        calc_orbits();
    }
    return true;
}

bool MetricSchwarzschildIsotropic::transToEmbedding(vec4 p, vec4& ep)
{
    double x = p[1];
    double y = p[2];
    double z = p[3];
    double rho = sqrt(x * x + y * y + z * z);
    double rs = 4.0 * rho_s;

    if (rho >= rho_s) {
        double r = rho * pow(1.0 + rho_s / rho, 2.0);
        double z = 2.0 * sqrt(rs) * sqrt(r - rs);
        ep = vec4(p[0], x, y, z);
        return true;
    }
    return false;
}

bool MetricSchwarzschildIsotropic::setEmbeddingParam(const char* name, double val)
{
    Metric::setEmbeddingParam(name, val);

    if (strcmp(name, "emb_rmin") == 0) {
        mEmb_rmin = val;
    }
    else if (strcmp(name, "emb_rmax") == 0) {
        mEmb_rmax = val;
    }
    else if (strcmp(name, "emb_r_num") == 0) {
        mEmb_r_num = static_cast<unsigned int>(val);
    }
    else if (strcmp(name, "emb_phi_num") == 0) {
        mEmb_phi_num = static_cast<unsigned int>(val);
    }
    return testEmbeddingParams();
}

/*! Test embedding parameters
 *  \return  true : all parameters are ok
 *  \return  false : at least one parameter had to be adjusted.
 */
bool MetricSchwarzschildIsotropic::testEmbeddingParams()
{
    bool allOk = true;
    if (mEmb_rmin < rho_s) {
        mEmb_rmin = rho_s;
        allOk &= false;
    }
    if (mEmb_rmax < rho_s) {
        mEmb_rmax = rho_s;
        allOk &= false;
    }
    if (mEmb_r_num < 2) {
        mEmb_r_num = 2;
        allOk &= false;
    }

    if (mEmb_phi_num < 4) {
        mEmb_phi_num = 4;
        allOk &= false;
    }
    return allOk;
}

unsigned int MetricSchwarzschildIsotropic::getEmbeddingVertices(
    float*& verts, unsigned int*& indices, unsigned int& numElems, unsigned int& counter)
{
    m4d::SafeDelete<float>(verts);
    m4d::SafeDelete<unsigned int>(indices);

    testEmbeddingParams();
    mEmb_rstep = (mEmb_rmax - mEmb_rmin) / static_cast<double>(mEmb_r_num);
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

    double rs = 4.0 * rho_s;

    for (unsigned int k = 0; k < counter; k++) {
        double phi = k * mEmb_phistep;

        for (unsigned int j = 0; j < numElems; j++) {
            double rho = mEmb_rmin + j * mEmb_rstep;
            double x = rho * cos(phi);
            double y = rho * sin(phi);
            double r = rho * pow(1.0 + rho_s / rho, 2.0);

            if (rho >= rho_s) {
                double z = 2.0 * sqrt(rs) * sqrt(fabs(r - rs));
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

bool MetricSchwarzschildIsotropic::report(const vec4 pos, const vec4, char*& text)
{
    std::stringstream ss;
    ss << "Report for Schwarzschild metric\n\t isotropic coordinates : (t,x,y,z)\n";
    ss << "---------------------------------------------------------------\n";
    ss << "  physical units ................................. no\n";
    ss.precision(DEF_FIXED_REPORT_PRECISION);
    ss.setf(std::ios::fixed);
    ss << "  Schwarzschild radius ...................... r_s = 2GM/c^2 = " << rho_s * 4.0 << std::endl;
    ss << "                                            rho_s = rs/4 = " << rho_s << std::endl;
    ss << "  Photon orbit ........................... rho_po = " << rho_po << std::endl;
    ss << "  innermost stable circular orbit ...... rho_isco = " << rho_lso << std::endl;

    double beta;
    calcBetaOfCircOrbit(pos.data(), beta);
    ss << "  velocity for circular orbit .............. beta = " << beta << std::endl;

    text = new char[ss.str().length() + 2];
    return CopyString(ss.str().c_str(), text);
}

bool MetricSchwarzschildIsotropic::calcBetaOfCircOrbit(const double* pos, double& beta)
{
    calc_orbits();
    double rho = sqrt(pos[1] * pos[1] + pos[2] * pos[2] + pos[3] * pos[3]);
    if (rho >= rho_lso) {
        beta = 1.0 / sqrt(2.0 * (0.25 * rho / rho_s * (1.0 + rho_s / rho) * (1.0 + rho_s / rho) - 1.0));
        return true;
    }
    beta = 0.0;
    return false;
}

void MetricSchwarzschildIsotropic::setStandardValues()
{
    mInitPos[0] = 0.0;
    mInitPos[1] = 6.0;
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

double MetricSchwarzschildIsotropic::calc_rho(const double* pos)
{
    rho = sqrt(pos[1] * pos[1] + pos[2] * pos[2] + pos[3] * pos[3]);
    return rho;
}

void MetricSchwarzschildIsotropic::calc_drho(const double* pos, double& drdx, double& drdy, double& drdz)
{
    rho = calc_rho(pos);
    drdx = pos[1] / rho;
    drdy = pos[2] / rho;
    drdz = pos[3] / rho;
}

void MetricSchwarzschildIsotropic::calc_d2rho(
    const double* pos, double& drdxdx, double& drdxdy, double& drdxdz, double& drdydy, double& drdydz, double& drdzdz)
{
    rho = calc_rho(pos);
    double x = pos[1];
    double y = pos[2];
    double z = pos[3];
    double edrh3 = pow(rho, -3.0);

    drdxdx = (y * y + z * z) * edrh3;
    drdxdy = -x * y * edrh3;
    drdxdz = -x * y * edrh3;
    drdydy = (x * x + z * z) * edrh3;
    drdydz = -y * z * edrh3;
    drdzdz = (x * x + y * y) * edrh3;
}

void MetricSchwarzschildIsotropic::calc_orbits()
{
    rho_po = (2.0 + sqrt(3.0)) * rho_s;
    rho_lso = (5.0 + 2.0 * sqrt(6.0)) * rho_s;
}

} // end namespace m4d
