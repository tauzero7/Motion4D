// -------------------------------------------------------------------------------
/*
   m4dMetricAlcubierreAccel.cpp

  Copyright (c) 2009-2014-2011  Thomas Mueller, Frank Grave


   This file is part of the m4d-library.

   The m4d-library is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   The m4d-library is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with the m4d-library.  If not, see <http://www.gnu.org/licenses/>.

*/
// -------------------------------------------------------------------------------

#include "m4dMetricAlcubierreAccel.h"

namespace m4d {

#define eps 1.0e-6


/*! Standard constructor for the KerrBL metric.
 *
 * \param  sigma : sharpness of warp bubble.
 * \param  R     : size of warp bubble.
 * \param  x0    : position ofwarp  bubble at t=t0.
 * \param  v0    : constant velocity part of warp bubble.
 * \param  alpha : acceleration of warp bubble.
 * \param  t0    : time synchronization.
 */
MetricAlcubierreAccel::MetricAlcubierreAccel(double sigma, double R, double x0, double v0, double alpha, double t0) {
    mMetricName  = "AlcubierreWarpAccel";
    setCoordType(enum_coordinate_cartesian);

    mPhysicalUnits = enum_physical_constants_geom;
    mSpeedOfLight = 1.0;
    mGravConstant = 1.0;

    mSigma = sigma;
    mR     = R;
    mx0    = x0;
    mv0    = v0;
    malpha = alpha;
    mt0    = t0;

    addParam("sigma", sigma);
    addParam("r", R);
    addParam("x0", x0);
    addParam("v0", v0);
    addParam("alpha", alpha);
    addParam("t0", t0);

    mDrawTypes.push_back(enum_draw_twoplusone);

    setStandardValues();

    mLocTeds.push_back(enum_nat_tetrad_comoving);
    mLocTeds.push_back(enum_nat_tetrad_static);
}

MetricAlcubierreAccel::~MetricAlcubierreAccel() {
}


// *********************************** public methods ******************************
/*! Calculate the contravariant metric components at position 'pos'.
 *
 *  \param pos : pointer to position.
 */
bool MetricAlcubierreAccel::calculateMetric(const double* pos) {
    double c = mSpeedOfLight;
    double f = calcF(pos);

    double v = mv0;
    if (pos[0] > mt0) {
        v += malpha * (pos[0] - mt0);
    }

    double t1 = c * c;
    double t2 = v;        // v(t);
    double t3 = t2 * t2;
    double t4 = f;        // f(t,x,y,z);
    double t5 = t4 * t4;
    double t8 = t2 * t4;

    g_compts[0][0] = -t1 + t3 * t5;
    g_compts[0][1] = -t8;
    g_compts[0][2] = 0.0;
    g_compts[0][3] = 0.0;
    g_compts[1][0] = -t8;
    g_compts[1][1] = 1.0;
    g_compts[1][2] = 0.0;
    g_compts[1][3] = 0.0;
    g_compts[2][0] = 0.0;
    g_compts[2][1] = 0.0;
    g_compts[2][2] = 1.0;
    g_compts[2][3] = 0.0;
    g_compts[3][0] = 0.0;
    g_compts[3][1] = 0.0;
    g_compts[3][2] = 0.0;
    g_compts[3][3] = 1.0;

    return true;
}

/*! Calculate the Christoffel symbols of the second kind at position 'pos'.
 *
 *  \param pos : pointer to position.
 */
bool MetricAlcubierreAccel::calculateChristoffels(const double* pos) {
    double c = mSpeedOfLight;
    double f = calcF(pos);

    double v = mv0;
    double dv = 0.0;
    if (pos[0] > mt0) {
        v += malpha * (pos[0] - mt0);
        dv = malpha;
    }
    double ft, fx, fy, fz;
    calcDF(pos, ft, fx, fy, fz);

    double t1 = v;                 // v(t);
    double t2 = t1 * t1;
    double t4 = f;                 // f(t,x,y,z);
    double t5 = t4 * t4;
    double t7 = fx;                // diff(f(t,x,y,z),x);
    double t8 = c * c;
    double t9 = 1 / t8;
    double t11 = t2 * t1 * t5 * t7 * t9;
    double t12 = dv;               // diff(v(t),t);
    double t16 = ft;               // diff(f(t,x,y,z),t);
    double t19 = t4 * t7;
    double t21 = t2 * t2;
    double t27 = t2 * t4;
    double t28 = fy;               // diff(f(t,x,y,z),y);
    double t30 = fz;               // diff(f(t,x,y,z),z);
    double t32 = t9 * t2;
    double t33 = t32 * t19;
    double t34 = t1 * t28;
    double t35 = t34 / 2.0;
    double t36 = t1 * t30;
    double t37 = t36 / 2.0;
    double t40 = t32 * t4 * t28 / 2.0;
    double t43 = (t2 * t5 + t8) * t9;
    double t45 = t34 * t43 / 2.0;
    double t48 = t32 * t4 * t30 / 2.0;
    double t50 = t36 * t43 / 2.0;
    double t51 = t9 * t1;
    double t54 = t51 * t28 / 2.0;
    double t56 = t51 * t30 / 2.0;

    christoffel[0][0][0] = t11;
    christoffel[0][0][1] = -(t8 * t12 * t4 + t8 * t1 * t16 + t8 * t2 * t19 - t21 * t5 * t4 * t7) * t9;
    christoffel[0][0][2] = -t27 * t28;
    christoffel[0][0][3] = -t27 * t30;
    christoffel[0][1][0] = -t33;
    christoffel[0][1][1] = -t11;
    christoffel[0][1][2] = t35;
    christoffel[0][1][3] = t37;
    christoffel[0][2][0] = -t40;
    christoffel[0][2][1] = -t45;
    christoffel[0][2][2] = 0.0;
    christoffel[0][2][3] = 0.0;
    christoffel[0][3][0] = -t48;
    christoffel[0][3][1] = -t50;
    christoffel[0][3][2] = 0.0;
    christoffel[0][3][3] = 0.0;
    christoffel[1][0][0] = -t33;
    christoffel[1][0][1] = -t11;
    christoffel[1][0][2] = t35;
    christoffel[1][0][3] = t37;
    christoffel[1][1][0] = t51 * t7;
    christoffel[1][1][1] = t33;
    christoffel[1][1][2] = 0.0;
    christoffel[1][1][3] = 0.0;
    christoffel[1][2][0] = t54;
    christoffel[1][2][1] = t40;
    christoffel[1][2][2] = 0.0;
    christoffel[1][2][3] = 0.0;
    christoffel[1][3][0] = t56;
    christoffel[1][3][1] = t48;
    christoffel[1][3][2] = 0.0;
    christoffel[1][3][3] = 0.0;
    christoffel[2][0][0] = -t40;
    christoffel[2][0][1] = -t45;
    christoffel[2][0][2] = 0.0;
    christoffel[2][0][3] = 0.0;
    christoffel[2][1][0] = t54;
    christoffel[2][1][1] = t40;
    christoffel[2][1][2] = 0.0;
    christoffel[2][1][3] = 0.0;
    christoffel[2][2][0] = 0.0;
    christoffel[2][2][1] = 0.0;
    christoffel[2][2][2] = 0.0;
    christoffel[2][2][3] = 0.0;
    christoffel[2][3][0] = 0.0;
    christoffel[2][3][1] = 0.0;
    christoffel[2][3][2] = 0.0;
    christoffel[2][3][3] = 0.0;
    christoffel[3][0][0] = -t48;
    christoffel[3][0][1] = -t50;
    christoffel[3][0][2] = 0.0;
    christoffel[3][0][3] = 0.0;
    christoffel[3][1][0] = t56;
    christoffel[3][1][1] = t48;
    christoffel[3][1][2] = 0.0;
    christoffel[3][1][3] = 0.0;
    christoffel[3][2][0] = 0.0;
    christoffel[3][2][1] = 0.0;
    christoffel[3][2][2] = 0.0;
    christoffel[3][2][3] = 0.0;
    christoffel[3][3][0] = 0.0;
    christoffel[3][3][1] = 0.0;
    christoffel[3][3][2] = 0.0;
    christoffel[3][3][3] = 0.0;

    return true;
}

/*! Calculate Jacobi matrix.
 *
 *  \param pos : pointer to position.
 */
bool MetricAlcubierreAccel::calculateChrisD(const double* pos) {
    double c = mSpeedOfLight;
    double v = mv0;
    if (pos[0] > mt0) {
        v += malpha * (pos[0] - mt0);
    }

    double f = calcF(pos);

    double ft, fx, fy, fz;
    calcDF(pos, ft, fx, fy, fz);

    double ftt, ftx, fty, ftz, fxx, fxy, fxz, fyy, fyz, fzz;
    calcD2F(pos, ftt, ftx, fty, ftz, fxx, fxy, fxz, fyy, fyz, fzz);

    double t1 = v;                 // v(t);
    double t2 = t1 * t1;
    double t3 = f;                 // f(t,x,y,z);
    double t4 = t2 * t3;
    double t5 = fx;                // diff(f(t,x,y,z),x);
    double t7 = malpha;            // diff(v(t),t);
    double t8 = t3 * t5 * t7;
    double t11 = ft;               // diff(f(t,x,y,z),t);
    double t12 = t1 * t5 * t11;
    double t14 = t1 * t3;
    double t15 = ftx;              // diff(diff(f(t,x,y,z),t),x);
    double t16 = t14 * t15;
    double t18 = c * c;
    double t19 = 1 / t18;
    double t21 = t4 * (3.0 * t8 + 2.0 * t12 + t16) * t19;
    double t22 = t2 * t1;
    double t23 = t22 * t3;
    double t24 = t5 * t5;
    double t26 = fxx;              // diff(diff(f(t,x,y,z),x),x);
    double t27 = t3 * t26;
    double t30 = t23 * (2.0 * t24 + t27) * t19;
    double t31 = fy;               // diff(f(t,x,y,z),y);
    double t32 = t5 * t31;
    double t34 = fxy;              // diff(diff(f(t,x,y,z),x),y);
    double t35 = t3 * t34;
    double t38 = t23 * (2.0 * t32 + t35) * t19;
    double t39 = fz;               // diff(f(t,x,y,z),z);
    double t40 = t5 * t39;
    double t42 = fxz;              // diff(diff(f(t,x,y,z),x),z);
    double t43 = t3 * t42;
    double t46 = t23 * (2.0 * t40 + t43) * t19;
    double t47 = 0.0;              // diff(diff(v(t),t),t);
    double t50 = t18 * t7;
    double t53 = t18 * t1;
    double t54 = ftt;              // diff(diff(f(t,x,y,z),t),t);
    double t58 = t18 * t2;
    double t59 = t11 * t5;
    double t63 = t3 * t3;
    double t64 = t63 * t3;
    double t66 = t7 * t5;
    double t69 = t2 * t2;
    double t70 = t69 * t63;
    double t73 = t69 * t64;
    double t86 = t50 * t31;
    double t87 = fty;              // diff(diff(f(t,x,y,z),t),y);
    double t88 = t53 * t87;
    double t96 = t50 * t39;
    double t97 = ftz;              // diff(diff(f(t,x,y,z),t),z);
    double t98 = t53 * t97;
    double t106 = t7 * t31;
    double t109 = t2 * t11;
    double t113 = t2 * t5;
    double t117 = t31 * t31;
    double t118 = t2 * t117;
    double t119 = fyy;             // diff(diff(f(t,x,y,z),y),y);
    double t124 = fyz;             // diff(diff(f(t,x,y,z),y),z);
    double t126 = -t2 * t39 * t31 - t4 * t124;
    double t127 = t7 * t39;
    double t136 = t39 * t39;
    double t137 = t2 * t136;
    double t138 = fzz;             // diff(diff(f(t,x,y,z),z),z);
    double t144 = t1 * (2.0 * t8 + t12 + t16) * t19;
    double t147 = t2 * (t24 + t27) * t19;
    double t150 = t2 * (t32 + t35) * t19;
    double t153 = t2 * (t40 + t43) * t19;
    double t155 = t106 + t1 * t87;
    double t157 = t1 * t34 / 2.0;
    double t159 = t1 * t119 / 2.0;
    double t161 = t1 * t124 / 2.0;
    double t163 = t127 + t1 * t97;
    double t165 = t1 * t42 / 2.0;
    double t167 = t1 * t138 / 2.0;
    double t171 = t1 * t11;
    double t177 = t1 * (2.0 * t3 * t31 * t7 + t171 * t31 + t14 * t87) * t19 / 2.0;
    double t178 = t150 / 2.0;
    double t183 = t2 * (t117 + t3 * t119) * t19 / 2.0;
    double t189 = t2 * (t39 * t31 + t3 * t124) * t19 / 2.0;
    double t190 = t2 * t63;
    double t196 = t3 * t11;
    double t201 = (3.0 * t106 * t190 + t86 + t22 * t87 * t63 + t88 + 2.0 * t22 * t31 * t196) * t19 / 2.0;
    double t210 = t1 * (t34 * t2 * t63 + t34 * t18 + 2.0 * t4 * t32) * t19 / 2.0;
    double t219 = t1 * (t119 * t2 * t63 + t119 * t18 + 2.0 * t118 * t3) * t19 / 2.0;
    double t224 = t3 * t39;
    double t230 = t1 * (t124 * t2 * t63 + t124 * t18 + 2.0 * t2 * t31 * t224) * t19 / 2.0;
    double t238 = t1 * (2.0 * t224 * t7 + t171 * t39 + t14 * t97) * t19 / 2.0;
    double t239 = t153 / 2.0;
    double t244 = t2 * (t136 + t3 * t138) * t19 / 2.0;
    double t254 = (3.0 * t127 * t190 + t96 + t22 * t97 * t63 + t98 + 2.0 * t22 * t39 * t196) * t19 / 2.0;
    double t263 = t1 * (t42 * t2 * t63 + t42 * t18 + 2.0 * t4 * t40) * t19 / 2.0;
    double t272 = t1 * (t138 * t2 * t63 + t138 * t18 + 2.0 * t137 * t3) * t19 / 2.0;
    double t276 = t19 * t1;
    double t278 = t276 * t34;
    double t279 = t276 * t42;
    double t281 = t155 * t19 / 2.0;
    double t282 = t278 / 2.0;
    double t284 = t276 * t119 / 2.0;
    double t286 = t276 * t124 / 2.0;
    double t288 = t163 * t19 / 2.0;
    double t289 = t279 / 2.0;
    double t291 = t276 * t138 / 2.0;

    chrisD[0][0][0][0] = t21;
    chrisD[0][0][0][1] = t30;
    chrisD[0][0][0][2] = t38;
    chrisD[0][0][0][3] = t46;
    chrisD[0][0][1][0] = (-t18 * t47 * t3 - 2.0 * t50 * t11 - t53 * t54 - 2.0 * t53 * t8 - t58 * t59 - t58 * t3 * t15 + 4.0 * t22 * t64 * t66 + 3.0 * t70 * t59 + t73 * t15) * t19;
    chrisD[0][0][1][1] = (-t50 * t5 - t53 * t15 - t58 * t24 - t58 * t27 + 3.0 * t70 * t24 + t73 * t26) * t19;
    chrisD[0][0][1][2] = (-t86 - t88 - t58 * t32 - t58 * t35 + 3.0 * t70 * t32 + t73 * t34) * t19;
    chrisD[0][0][1][3] = (-t96 - t98 - t58 * t40 - t58 * t43 + 3.0 * t70 * t40 + t73 * t42) * t19;
    chrisD[0][0][2][0] = -2.0 * t14 * t106 - t109 * t31 - t4 * t87;
    chrisD[0][0][2][1] = -t113 * t31 - t4 * t34;
    chrisD[0][0][2][2] = -t118 - t4 * t119;
    chrisD[0][0][2][3] = t126;
    chrisD[0][0][3][0] = -2.0 * t14 * t127 - t109 * t39 - t4 * t97;
    chrisD[0][0][3][1] = -t113 * t39 - t4 * t42;
    chrisD[0][0][3][2] = t126;
    chrisD[0][0][3][3] = -t137 - t4 * t138;
    chrisD[0][1][0][0] = -t144;
    chrisD[0][1][0][1] = -t147;
    chrisD[0][1][0][2] = -t150;
    chrisD[0][1][0][3] = -t153;
    chrisD[0][1][1][0] = -t21;
    chrisD[0][1][1][1] = -t30;
    chrisD[0][1][1][2] = -t38;
    chrisD[0][1][1][3] = -t46;
    chrisD[0][1][2][0] = t155 / 2.0;
    chrisD[0][1][2][1] = t157;
    chrisD[0][1][2][2] = t159;
    chrisD[0][1][2][3] = t161;
    chrisD[0][1][3][0] = t163 / 2.0;
    chrisD[0][1][3][1] = t165;
    chrisD[0][1][3][2] = t161;
    chrisD[0][1][3][3] = t167;
    chrisD[0][2][0][0] = -t177;
    chrisD[0][2][0][1] = -t178;
    chrisD[0][2][0][2] = -t183;
    chrisD[0][2][0][3] = -t189;
    chrisD[0][2][1][0] = -t201;
    chrisD[0][2][1][1] = -t210;
    chrisD[0][2][1][2] = -t219;
    chrisD[0][2][1][3] = -t230;
    chrisD[0][2][2][0] = 0.0;
    chrisD[0][2][2][1] = 0.0;
    chrisD[0][2][2][2] = 0.0;
    chrisD[0][2][2][3] = 0.0;
    chrisD[0][2][3][0] = 0.0;
    chrisD[0][2][3][1] = 0.0;
    chrisD[0][2][3][2] = 0.0;
    chrisD[0][2][3][3] = 0.0;
    chrisD[0][3][0][0] = -t238;
    chrisD[0][3][0][1] = -t239;
    chrisD[0][3][0][2] = -t189;
    chrisD[0][3][0][3] = -t244;
    chrisD[0][3][1][0] = -t254;
    chrisD[0][3][1][1] = -t263;
    chrisD[0][3][1][2] = -t230;
    chrisD[0][3][1][3] = -t272;
    chrisD[0][3][2][0] = 0.0;
    chrisD[0][3][2][1] = 0.0;
    chrisD[0][3][2][2] = 0.0;
    chrisD[0][3][2][3] = 0.0;
    chrisD[0][3][3][0] = 0.0;
    chrisD[0][3][3][1] = 0.0;
    chrisD[0][3][3][2] = 0.0;
    chrisD[0][3][3][3] = 0.0;
    chrisD[1][0][0][0] = -t144;
    chrisD[1][0][0][1] = -t147;
    chrisD[1][0][0][2] = -t150;
    chrisD[1][0][0][3] = -t153;
    chrisD[1][0][1][0] = -t21;
    chrisD[1][0][1][1] = -t30;
    chrisD[1][0][1][2] = -t38;
    chrisD[1][0][1][3] = -t46;
    chrisD[1][0][2][0] = t155 / 2.0;
    chrisD[1][0][2][1] = t157;
    chrisD[1][0][2][2] = t159;
    chrisD[1][0][2][3] = t161;
    chrisD[1][0][3][0] = t163 / 2.0;
    chrisD[1][0][3][1] = t165;
    chrisD[1][0][3][2] = t161;
    chrisD[1][0][3][3] = t167;
    chrisD[1][1][0][0] = (t66 + t1 * t15) * t19;
    chrisD[1][1][0][1] = t276 * t26;
    chrisD[1][1][0][2] = t278;
    chrisD[1][1][0][3] = t279;
    chrisD[1][1][1][0] = t144;
    chrisD[1][1][1][1] = t147;
    chrisD[1][1][1][2] = t150;
    chrisD[1][1][1][3] = t153;
    chrisD[1][1][2][0] = 0.0;
    chrisD[1][1][2][1] = 0.0;
    chrisD[1][1][2][2] = 0.0;
    chrisD[1][1][2][3] = 0.0;
    chrisD[1][1][3][0] = 0.0;
    chrisD[1][1][3][1] = 0.0;
    chrisD[1][1][3][2] = 0.0;
    chrisD[1][1][3][3] = 0.0;
    chrisD[1][2][0][0] = t281;
    chrisD[1][2][0][1] = t282;
    chrisD[1][2][0][2] = t284;
    chrisD[1][2][0][3] = t286;
    chrisD[1][2][1][0] = t177;
    chrisD[1][2][1][1] = t178;
    chrisD[1][2][1][2] = t183;
    chrisD[1][2][1][3] = t189;
    chrisD[1][2][2][0] = 0.0;
    chrisD[1][2][2][1] = 0.0;
    chrisD[1][2][2][2] = 0.0;
    chrisD[1][2][2][3] = 0.0;
    chrisD[1][2][3][0] = 0.0;
    chrisD[1][2][3][1] = 0.0;
    chrisD[1][2][3][2] = 0.0;
    chrisD[1][2][3][3] = 0.0;
    chrisD[1][3][0][0] = t288;
    chrisD[1][3][0][1] = t289;
    chrisD[1][3][0][2] = t286;
    chrisD[1][3][0][3] = t291;
    chrisD[1][3][1][0] = t238;
    chrisD[1][3][1][1] = t239;
    chrisD[1][3][1][2] = t189;
    chrisD[1][3][1][3] = t244;
    chrisD[1][3][2][0] = 0.0;
    chrisD[1][3][2][1] = 0.0;
    chrisD[1][3][2][2] = 0.0;
    chrisD[1][3][2][3] = 0.0;
    chrisD[1][3][3][0] = 0.0;
    chrisD[1][3][3][1] = 0.0;
    chrisD[1][3][3][2] = 0.0;
    chrisD[1][3][3][3] = 0.0;
    chrisD[2][0][0][0] = -t177;
    chrisD[2][0][0][1] = -t178;
    chrisD[2][0][0][2] = -t183;
    chrisD[2][0][0][3] = -t189;
    chrisD[2][0][1][0] = -t201;
    chrisD[2][0][1][1] = -t210;
    chrisD[2][0][1][2] = -t219;
    chrisD[2][0][1][3] = -t230;
    chrisD[2][0][2][0] = 0.0;
    chrisD[2][0][2][1] = 0.0;
    chrisD[2][0][2][2] = 0.0;
    chrisD[2][0][2][3] = 0.0;
    chrisD[2][0][3][0] = 0.0;
    chrisD[2][0][3][1] = 0.0;
    chrisD[2][0][3][2] = 0.0;
    chrisD[2][0][3][3] = 0.0;
    chrisD[2][1][0][0] = t281;
    chrisD[2][1][0][1] = t282;
    chrisD[2][1][0][2] = t284;
    chrisD[2][1][0][3] = t286;
    chrisD[2][1][1][0] = t177;
    chrisD[2][1][1][1] = t178;
    chrisD[2][1][1][2] = t183;
    chrisD[2][1][1][3] = t189;
    chrisD[2][1][2][0] = 0.0;
    chrisD[2][1][2][1] = 0.0;
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
    chrisD[2][2][1][1] = 0.0;
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
    chrisD[2][3][3][2] = 0.0;
    chrisD[2][3][3][3] = 0.0;
    chrisD[3][0][0][0] = -t238;
    chrisD[3][0][0][1] = -t239;
    chrisD[3][0][0][2] = -t189;
    chrisD[3][0][0][3] = -t244;
    chrisD[3][0][1][0] = -t254;
    chrisD[3][0][1][1] = -t263;
    chrisD[3][0][1][2] = -t230;
    chrisD[3][0][1][3] = -t272;
    chrisD[3][0][2][0] = 0.0;
    chrisD[3][0][2][1] = 0.0;
    chrisD[3][0][2][2] = 0.0;
    chrisD[3][0][2][3] = 0.0;
    chrisD[3][0][3][0] = 0.0;
    chrisD[3][0][3][1] = 0.0;
    chrisD[3][0][3][2] = 0.0;
    chrisD[3][0][3][3] = 0.0;
    chrisD[3][1][0][0] = t288;
    chrisD[3][1][0][1] = t289;
    chrisD[3][1][0][2] = t286;
    chrisD[3][1][0][3] = t291;
    chrisD[3][1][1][0] = t238;
    chrisD[3][1][1][1] = t239;
    chrisD[3][1][1][2] = t189;
    chrisD[3][1][1][3] = t244;
    chrisD[3][1][2][0] = 0.0;
    chrisD[3][1][2][1] = 0.0;
    chrisD[3][1][2][2] = 0.0;
    chrisD[3][1][2][3] = 0.0;
    chrisD[3][1][3][0] = 0.0;
    chrisD[3][1][3][1] = 0.0;
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
    chrisD[3][2][3][2] = 0.0;
    chrisD[3][2][3][3] = 0.0;
    chrisD[3][3][0][0] = 0.0;
    chrisD[3][3][0][1] = 0.0;
    chrisD[3][3][0][2] = 0.0;
    chrisD[3][3][0][3] = 0.0;
    chrisD[3][3][1][0] = 0.0;
    chrisD[3][3][1][1] = 0.0;
    chrisD[3][3][1][2] = 0.0;
    chrisD[3][3][1][3] = 0.0;
    chrisD[3][3][2][0] = 0.0;
    chrisD[3][3][2][1] = 0.0;
    chrisD[3][3][2][2] = 0.0;
    chrisD[3][3][2][3] = 0.0;
    chrisD[3][3][3][0] = 0.0;
    chrisD[3][3][3][1] = 0.0;
    chrisD[3][3][3][2] = 0.0;
    chrisD[3][3][3][3] = 0.0;

    return true;
}

/*! Calculate Riemann tensor R^a_bcd
 * \param pos : pointer to coordinate position where the Riemann tensor have to be evaluated.
 * \return true : successfull
 */
bool MetricAlcubierreAccel::calculateRiemann(const double* pos) {
    double c  = mSpeedOfLight;
    double v  = mv0;
    double dv = 0.0;
    if (pos[0] > mt0) {
        v += malpha * (pos[0] - mt0);
        dv = malpha;
    }
    double ft, fx, fy, fz;
    calcDF(pos, ft, fx, fy, fz);

    double ftt, ftx, fty, ftz, fxx, fxy, fxz, fyy, fyz, fzz;
    calcD2F(pos, ftt, ftx, fty, ftz, fxx, fxy, fxz, fyy, fyz, fzz);

    double t4 = v;               // v(t);
    double t5 = calcF(pos);      // f(t,x,y,z);
    double t6 = t4 * t5;
    double t7 = c * c;
    double t8 = 1 / t7;
    double t9 = dv;              // diff(v(t),t);
    double t10 = fx;             // diff(f(t,x,y,z),x);
    double t13 = ftx;            // diff(diff(f(t,x,y,z),t),x);
    double t16 = t4 * t4;
    double t17 = t10 * t10;
    double t20 = t16 * t5;
    double t21 = fxx;            // diff(diff(f(t,x,y,z),x),x);
    double t24 = fy;             // diff(f(t,x,y,z),y);
    double t25 = t24 * t24;
    double t26 = t16 * t25;
    double t27 = fz;             // diff(f(t,x,y,z),z);
    double t28 = t27 * t27;
    double t29 = t16 * t28;
    double t30 = 4.0 * t9 * t10 + 4.0 * t4 * t13 + 4.0 * t16 * t17 + 4.0 * t20 * t21 - t26 - t29;
    double t31 = t8 * t30;
    double t33 = t6 * t31 / 4.0;
    double t34 = t9 * t24;
    double t35 = fty;            // diff(diff(f(t,x,y,z),t),y);
    double t36 = t4 * t35;
    double t38 = t16 * t24 * t10;
    double t39 = 2.0 * t38;
    double t40 = fxy;            // diff(diff(f(t,x,y,z),x),y);
    double t41 = t20 * t40;
    double t43 = t34 + t36 + t39 + 2.0 * t41;
    double t44 = t8 * t43;
    double t46 = t6 * t44 / 2.0;
    double t47 = t9 * t27;
    double t48 = ftz;            // diff(diff(f(t,x,y,z),t),z);
    double t49 = t4 * t48;
    double t51 = t16 * t27 * t10;
    double t52 = 2.0 * t51;
    double t53 = fxz;            // diff(diff(f(t,x,y,z),x),z);
    double t54 = t20 * t53;
    double t56 = t47 + t49 + t52 + 2.0 * t54;
    double t57 = t8 * t56;
    double t59 = t6 * t57 / 2.0;
    double t62 = t20 * t8 * t40 / 2.0;
    double t65 = t20 * t8 * t53 / 2.0;
    double t66 = t31 / 4.0;
    double t67 = t44 / 2.0;
    double t68 = t57 / 2.0;
    double t69 = t4 * t40;
    double t71 = t69 * t8 / 2.0;
    double t72 = t4 * t53;
    double t74 = t72 * t8 / 2.0;
    double t77 = (t34 + t36 + t39 + t41) * t8 / 2.0;
    double t78 = t16 * t8;
    double t80 = fyy;            // diff(diff(f(t,x,y,z),y),y);
    double t81 = t5 * t80;
    double t82 = 2.0 * t81;
    double t85 = t78 * (3.0 * t25 + t82) / 4.0;
    double t86 = t24 * t27;
    double t88 = fyz;            // diff(diff(f(t,x,y,z),y),z);
    double t89 = t5 * t88;
    double t90 = 2.0 * t89;
    double t93 = t78 * (3.0 * t86 + t90) / 4.0;
    double t96 = t4 * t80 * t8 / 2.0;
    double t99 = t4 * t88 * t8 / 2.0;
    double t102 = (t47 + t49 + t52 + t54) * t8 / 2.0;
    double t104 = fzz;           // diff(diff(f(t,x,y,z),z),z);
    double t105 = t5 * t104;
    double t106 = 2.0 * t105;
    double t109 = t78 * (3.0 * t28 + t106) / 4.0;
    double t112 = t4 * t104 * t8 / 2.0;
    double t113 = t5 * t5;
    double t114 = t16 * t113;
    double t115 = t7 - t114;
    double t118 = t115 * t30 * t8 / 4.0;
    double t121 = t115 * t43 * t8 / 2.0;
    double t124 = t115 * t56 * t8 / 2.0;
    double t125 = t115 * t8;
    double t127 = t125 * t69 / 2.0;
    double t129 = t125 * t72 / 2.0;
    double t130 = t5 * t9;
    double t141 = t4 * (t130 * t24 + t6 * t35 + 2.0 * t20 * t24 * t10 + t114 * t40 + t40 * t7) * t8 / 2.0;
    double t142 = t4 * t8;
    double t143 = t20 * t25;
    double t146 = t80 * t7;
    double t149 = t142 * (2.0 * t143 + t114 * t80 + t146) / 2.0;
    double t150 = t20 * t86;
    double t153 = t88 * t7;
    double t156 = t142 * (2.0 * t150 + t114 * t88 + t153) / 2.0;
    double t159 = t78 * (t82 + t25) / 4.0;
    double t162 = t78 * (t90 + t86) / 4.0;
    double t173 = t4 * (t130 * t27 + t6 * t48 + 2.0 * t20 * t27 * t10 + t114 * t53 + t53 * t7) * t8 / 2.0;
    double t174 = t20 * t28;
    double t177 = t104 * t7;
    double t180 = t142 * (2.0 * t174 + t114 * t104 + t177) / 2.0;
    double t183 = t78 * (t106 + t28) / 4.0;
    double t186 = t34 / 2.0 + t36 / 2.0 + t38 + t41;
    double t195 = t16 * (3.0 * t25 * t7 + 4.0 * t81 * t7 + t26 * t113) * t8 / 4.0;
    double t204 = t16 * (3.0 * t86 * t7 + 4.0 * t89 * t7 + t86 * t114) * t8 / 4.0;
    double t209 = t4 * (2.0 * t146 + t143) * t8 / 4.0;
    double t214 = t4 * (2.0 * t153 + t150) * t8 / 4.0;
    double t215 = t69 / 2.0;
    double t217 = t78 * t25 / 4.0;
    double t219 = t78 * t86 / 4.0;
    double t222 = t47 / 2.0 + t49 / 2.0 + t51 + t54;
    double t231 = t16 * (3.0 * t28 * t7 + 4.0 * t105 * t7 + t29 * t113) * t8 / 4.0;
    double t236 = t4 * (2.0 * t177 + t174) * t8 / 4.0;
    double t237 = t72 / 2.0;
    double t239 = t78 * t28 / 4.0;

    riem[0][0][0][0] = 0.0;
    riem[0][0][0][1] = -t33;
    riem[0][0][0][2] = -t46;
    riem[0][0][0][3] = -t59;
    riem[0][0][1][0] = t33;
    riem[0][0][1][1] = 0.0;
    riem[0][0][1][2] = t62;
    riem[0][0][1][3] = t65;
    riem[0][0][2][0] = t46;
    riem[0][0][2][1] = -t62;
    riem[0][0][2][2] = 0.0;
    riem[0][0][2][3] = 0.0;
    riem[0][0][3][0] = t59;
    riem[0][0][3][1] = -t65;
    riem[0][0][3][2] = 0.0;
    riem[0][0][3][3] = 0.0;
    riem[0][1][0][0] = 0.0;
    riem[0][1][0][1] = t66;
    riem[0][1][0][2] = t67;
    riem[0][1][0][3] = t68;
    riem[0][1][1][0] = -t66;
    riem[0][1][1][1] = 0.0;
    riem[0][1][1][2] = -t71;
    riem[0][1][1][3] = -t74;
    riem[0][1][2][0] = -t67;
    riem[0][1][2][1] = t71;
    riem[0][1][2][2] = 0.0;
    riem[0][1][2][3] = 0.0;
    riem[0][1][3][0] = -t68;
    riem[0][1][3][1] = t74;
    riem[0][1][3][2] = 0.0;
    riem[0][1][3][3] = 0.0;
    riem[0][2][0][0] = 0.0;
    riem[0][2][0][1] = t77;
    riem[0][2][0][2] = t85;
    riem[0][2][0][3] = t93;
    riem[0][2][1][0] = -t77;
    riem[0][2][1][1] = 0.0;
    riem[0][2][1][2] = -t96;
    riem[0][2][1][3] = -t99;
    riem[0][2][2][0] = -t85;
    riem[0][2][2][1] = t96;
    riem[0][2][2][2] = 0.0;
    riem[0][2][2][3] = 0.0;
    riem[0][2][3][0] = -t93;
    riem[0][2][3][1] = t99;
    riem[0][2][3][2] = 0.0;
    riem[0][2][3][3] = 0.0;
    riem[0][3][0][0] = 0.0;
    riem[0][3][0][1] = t102;
    riem[0][3][0][2] = t93;
    riem[0][3][0][3] = t109;
    riem[0][3][1][0] = -t102;
    riem[0][3][1][1] = 0.0;
    riem[0][3][1][2] = -t99;
    riem[0][3][1][3] = -t112;
    riem[0][3][2][0] = -t93;
    riem[0][3][2][1] = t99;
    riem[0][3][2][2] = 0.0;
    riem[0][3][2][3] = 0.0;
    riem[0][3][3][0] = -t109;
    riem[0][3][3][1] = t112;
    riem[0][3][3][2] = 0.0;
    riem[0][3][3][3] = 0.0;
    riem[1][0][0][0] = 0.0;
    riem[1][0][0][1] = t118;
    riem[1][0][0][2] = t121;
    riem[1][0][0][3] = t124;
    riem[1][0][1][0] = -t118;
    riem[1][0][1][1] = 0.0;
    riem[1][0][1][2] = -t127;
    riem[1][0][1][3] = -t129;
    riem[1][0][2][0] = -t121;
    riem[1][0][2][1] = t127;
    riem[1][0][2][2] = 0.0;
    riem[1][0][2][3] = 0.0;
    riem[1][0][3][0] = -t124;
    riem[1][0][3][1] = t129;
    riem[1][0][3][2] = 0.0;
    riem[1][0][3][3] = 0.0;
    riem[1][1][0][0] = 0.0;
    riem[1][1][0][1] = t33;
    riem[1][1][0][2] = t46;
    riem[1][1][0][3] = t59;
    riem[1][1][1][0] = -t33;
    riem[1][1][1][1] = 0.0;
    riem[1][1][1][2] = -t62;
    riem[1][1][1][3] = -t65;
    riem[1][1][2][0] = -t46;
    riem[1][1][2][1] = t62;
    riem[1][1][2][2] = 0.0;
    riem[1][1][2][3] = 0.0;
    riem[1][1][3][0] = -t59;
    riem[1][1][3][1] = t65;
    riem[1][1][3][2] = 0.0;
    riem[1][1][3][3] = 0.0;
    riem[1][2][0][0] = 0.0;
    riem[1][2][0][1] = t141;
    riem[1][2][0][2] = t149;
    riem[1][2][0][3] = t156;
    riem[1][2][1][0] = -t141;
    riem[1][2][1][1] = 0.0;
    riem[1][2][1][2] = -t159;
    riem[1][2][1][3] = -t162;
    riem[1][2][2][0] = -t149;
    riem[1][2][2][1] = t159;
    riem[1][2][2][2] = 0.0;
    riem[1][2][2][3] = 0.0;
    riem[1][2][3][0] = -t156;
    riem[1][2][3][1] = t162;
    riem[1][2][3][2] = 0.0;
    riem[1][2][3][3] = 0.0;
    riem[1][3][0][0] = 0.0;
    riem[1][3][0][1] = t173;
    riem[1][3][0][2] = t156;
    riem[1][3][0][3] = t180;
    riem[1][3][1][0] = -t173;
    riem[1][3][1][1] = 0.0;
    riem[1][3][1][2] = -t162;
    riem[1][3][1][3] = -t183;
    riem[1][3][2][0] = -t156;
    riem[1][3][2][1] = t162;
    riem[1][3][2][2] = 0.0;
    riem[1][3][2][3] = 0.0;
    riem[1][3][3][0] = -t180;
    riem[1][3][3][1] = t183;
    riem[1][3][3][2] = 0.0;
    riem[1][3][3][3] = 0.0;
    riem[2][0][0][0] = 0.0;
    riem[2][0][0][1] = t186;
    riem[2][0][0][2] = t195;
    riem[2][0][0][3] = t204;
    riem[2][0][1][0] = -t186;
    riem[2][0][1][1] = 0.0;
    riem[2][0][1][2] = -t209;
    riem[2][0][1][3] = -t214;
    riem[2][0][2][0] = -t195;
    riem[2][0][2][1] = t209;
    riem[2][0][2][2] = 0.0;
    riem[2][0][2][3] = 0.0;
    riem[2][0][3][0] = -t204;
    riem[2][0][3][1] = t214;
    riem[2][0][3][2] = 0.0;
    riem[2][0][3][3] = 0.0;
    riem[2][1][0][0] = 0.0;
    riem[2][1][0][1] = -t215;
    riem[2][1][0][2] = -t209;
    riem[2][1][0][3] = -t214;
    riem[2][1][1][0] = t215;
    riem[2][1][1][1] = 0.0;
    riem[2][1][1][2] = t217;
    riem[2][1][1][3] = t219;
    riem[2][1][2][0] = t209;
    riem[2][1][2][1] = -t217;
    riem[2][1][2][2] = 0.0;
    riem[2][1][2][3] = 0.0;
    riem[2][1][3][0] = t214;
    riem[2][1][3][1] = -t219;
    riem[2][1][3][2] = 0.0;
    riem[2][1][3][3] = 0.0;
    riem[2][2][0][0] = 0.0;
    riem[2][2][0][1] = 0.0;
    riem[2][2][0][2] = 0.0;
    riem[2][2][0][3] = 0.0;
    riem[2][2][1][0] = 0.0;
    riem[2][2][1][1] = 0.0;
    riem[2][2][1][2] = 0.0;
    riem[2][2][1][3] = 0.0;
    riem[2][2][2][0] = 0.0;
    riem[2][2][2][1] = 0.0;
    riem[2][2][2][2] = 0.0;
    riem[2][2][2][3] = 0.0;
    riem[2][2][3][0] = 0.0;
    riem[2][2][3][1] = 0.0;
    riem[2][2][3][2] = 0.0;
    riem[2][2][3][3] = 0.0;
    riem[2][3][0][0] = 0.0;
    riem[2][3][0][1] = 0.0;
    riem[2][3][0][2] = 0.0;
    riem[2][3][0][3] = 0.0;
    riem[2][3][1][0] = 0.0;
    riem[2][3][1][1] = 0.0;
    riem[2][3][1][2] = 0.0;
    riem[2][3][1][3] = 0.0;
    riem[2][3][2][0] = 0.0;
    riem[2][3][2][1] = 0.0;
    riem[2][3][2][2] = 0.0;
    riem[2][3][2][3] = 0.0;
    riem[2][3][3][0] = 0.0;
    riem[2][3][3][1] = 0.0;
    riem[2][3][3][2] = 0.0;
    riem[2][3][3][3] = 0.0;
    riem[3][0][0][0] = 0.0;
    riem[3][0][0][1] = t222;
    riem[3][0][0][2] = t204;
    riem[3][0][0][3] = t231;
    riem[3][0][1][0] = -t222;
    riem[3][0][1][1] = 0.0;
    riem[3][0][1][2] = -t214;
    riem[3][0][1][3] = -t236;
    riem[3][0][2][0] = -t204;
    riem[3][0][2][1] = t214;
    riem[3][0][2][2] = 0.0;
    riem[3][0][2][3] = 0.0;
    riem[3][0][3][0] = -t231;
    riem[3][0][3][1] = t236;
    riem[3][0][3][2] = 0.0;
    riem[3][0][3][3] = 0.0;
    riem[3][1][0][0] = 0.0;
    riem[3][1][0][1] = -t237;
    riem[3][1][0][2] = -t214;
    riem[3][1][0][3] = -t236;
    riem[3][1][1][0] = t237;
    riem[3][1][1][1] = 0.0;
    riem[3][1][1][2] = t219;
    riem[3][1][1][3] = t239;
    riem[3][1][2][0] = t214;
    riem[3][1][2][1] = -t219;
    riem[3][1][2][2] = 0.0;
    riem[3][1][2][3] = 0.0;
    riem[3][1][3][0] = t236;
    riem[3][1][3][1] = -t239;
    riem[3][1][3][2] = 0.0;
    riem[3][1][3][3] = 0.0;
    riem[3][2][0][0] = 0.0;
    riem[3][2][0][1] = 0.0;
    riem[3][2][0][2] = 0.0;
    riem[3][2][0][3] = 0.0;
    riem[3][2][1][0] = 0.0;
    riem[3][2][1][1] = 0.0;
    riem[3][2][1][2] = 0.0;
    riem[3][2][1][3] = 0.0;
    riem[3][2][2][0] = 0.0;
    riem[3][2][2][1] = 0.0;
    riem[3][2][2][2] = 0.0;
    riem[3][2][2][3] = 0.0;
    riem[3][2][3][0] = 0.0;
    riem[3][2][3][1] = 0.0;
    riem[3][2][3][2] = 0.0;
    riem[3][2][3][3] = 0.0;
    riem[3][3][0][0] = 0.0;
    riem[3][3][0][1] = 0.0;
    riem[3][3][0][2] = 0.0;
    riem[3][3][0][3] = 0.0;
    riem[3][3][1][0] = 0.0;
    riem[3][3][1][1] = 0.0;
    riem[3][3][1][2] = 0.0;
    riem[3][3][1][3] = 0.0;
    riem[3][3][2][0] = 0.0;
    riem[3][3][2][1] = 0.0;
    riem[3][3][2][2] = 0.0;
    riem[3][3][2][3] = 0.0;
    riem[3][3][3][0] = 0.0;
    riem[3][3][3][1] = 0.0;
    riem[3][3][3][2] = 0.0;
    riem[3][3][3][3] = 0.0;

    return true;
}

/*! Transform local 4-direction to coordinate 4-direction.
 *
 *  \param  pos  :  pointer to position array.
 *  \param  ldir :  pointer to local direction array.
 *  \param  dir  :  pointer to calculated coordinate direction array.
 *  \param  type :  type of tetrad.
 */
void MetricAlcubierreAccel::localToCoord(const double* pos, const double* ldir, double* dir,
        enum_nat_tetrad_type  type) {
    double f = calcF(pos);
    double c = mSpeedOfLight;
    double v = mv0;
    if (pos[0] > mt0) {
        v += malpha * (pos[0] - mt0);
    }

    if (type == enum_nat_tetrad_comoving) {
        dir[0] = ldir[0] / c;
        dir[1] = ldir[0] * v * f / c + ldir[1];
        dir[2] = ldir[2];
        dir[3] = ldir[3];
    } else {
        double w = sqrt(c * c - v * v * f * f);

        dir[0] = (ldir[0] + v * f / c * ldir[1]) / w;
        dir[1] = -w / c * ldir[1];
        dir[2] = ldir[2];
        dir[3] = ldir[3];
    }
}

/*! Transform coordinate 4-direction to local 4-direction.
 *
 *  \param  pos  :  pointer to position array.
 *  \param  cdir :  pointer to coordinate direction.
 *  \param  ldir :  pointer to calculated local direction array.
 *  \param  type :  type of tetrad.
 */
void MetricAlcubierreAccel::coordToLocal(const double* pos, const double* cdir, double* ldir,
        enum_nat_tetrad_type  type) {
    double f = calcF(pos);
    double c = mSpeedOfLight;
    double v = mv0;
    if (pos[0] > mt0) {
        v += malpha * (pos[0] - mt0);
    }

    if (type == enum_nat_tetrad_comoving) {
        ldir[0] = c * cdir[0];
        ldir[1] = cdir[1] - v * f * cdir[0];
        ldir[2] = cdir[2];
        ldir[3] = cdir[3];
    } else {
        double w = sqrt(c * c - v * v * f * f);

        ldir[1] = c / w * cdir[1];
        ldir[0] = w * cdir[0] - ldir[1] * v * f / c;
        ldir[2] = cdir[2];
        ldir[3] = cdir[3];
    }
}


/*!
 *  \param pos  :  position.
 *  \return true  : radial position r < 0.0 or ...
 *  \return false : position is valid.
 */
bool MetricAlcubierreAccel::breakCondition(const double*) {
    bool br = false;
    return br;
}


/*! Tests whether the constraint equation is fulfilled.
 *
 *  The constraint equation for lightlike and timelike geodesics reads:
 \verbatim
     sum = g_{\mu\nu} dot(x)^{\mu} dot(x)^{\nu} - kappa c^2 = 0.
 \endverbatim
 *  \param  y[]   : pointer to position and direction coordinates.
 *  \param  kappa : timelike (-1.0), lightlike (0.0).
 *  \return double : sum.
 */
double MetricAlcubierreAccel::testConstraint(const double y[], const double kappa) {
    double c = mSpeedOfLight;
    double f = calcF(y);
    double v = mv0;
    if (y[0] > mt0) {
        v += malpha * (y[0] - mt0);
    }

    double sum = -kappa;
    sum += -c * c * y[4] * y[4] + pow(y[5] - v * f * y[4], 2.0) + y[6] * y[6] + y[7] * y[7];
    return sum;
}


/*! Set parameter 'pName' to 'val'.
 *
 *  Set 'sigma', 'R', and 'vs' parameters.
 */
bool MetricAlcubierreAccel::setParam(std::string pName, double val) {
    Metric::setParam(pName, val);
    if (pName == "sigma") {
        mSigma = val;
    } else if (pName == "r") {
        mR = val;
    } else if (pName == "v0") {
        mv0 = val;
    } else if (pName == "x0") {
        mx0 = val;
    } else if (pName == "alpha") {
        malpha = val;
    } else if (pName == "t0") {
        mt0 = val;
    }
    return true;
}

/*! Transform point p to 2+1 coordinates.
 *
 *  \param  p  : point in proper metric coordinates.
 *  \param  cp : reference to transformed point.
 *  \return true : success.
 */
bool MetricAlcubierreAccel::transToTwoPlusOne(vec4 p, vec4 &cp) {
    cp = vec4(p[0], p[1], p[2], p[0]);
    return true;
}

/*! Generate report.
 */
bool MetricAlcubierreAccel::report(const vec4 , const vec4 , std::string &text) {
    std::stringstream ss;
    ss << "Report for AlcubierreWarp metric\n\tcoordinate : (t,x,y,z)\n";
    ss << "---------------------------------------------------------------\n";
    ss << "  physical units ................................. no\n";
    ss.precision(DEF_FIXED_REPORT_PRECISION);
    ss.setf(std::ios::fixed);


    text = ss.str();
    return true;
}

// *************************** specific  public methods ****************************


// ********************************* protected methods *****************************
/*!
 */
void MetricAlcubierreAccel::setStandardValues() {
    mInitPos[0] = 0.0;
    mInitPos[1] = 0.0;
    mInitPos[2] = -10.0;
    mInitPos[3] = 0.0;
    mInitDir[0] = 1.0;
    mInitDir[1] = 0.0;
    mInitDir[2] = 0.0;

    mCoordNames[0] = std::string("t");
    mCoordNames[1] = std::string("x");
    mCoordNames[2] = std::string("y");
    mCoordNames[3] = std::string("z");
}

/*! Calculate rs function
 *  \param  pos : pointer to position.
 */
double MetricAlcubierreAccel::calcRs(const double* pos) {
    double t = pos[0];
    double x = pos[1];
    double y = pos[2];
    double z = pos[3];

    double xt = mx0 + mv0 * (t - mt0);
    double v  = mv0;
    if (t > mt0) {
        xt += 0.5 * malpha * (t - mt0) * (t - mt0);
        v  += malpha * (t - mt0);
    }
    return sqrt((x - xt) * (x - xt) + y * y + z * z);
}

double MetricAlcubierreAccel::calcF(const double* pos) {
    double rs = calcRs(pos);
#if 1
    double w1 = tanh(mSigma * (rs + mR));
    double w2 = tanh(mSigma * (rs - mR));
    double w3 = tanh(mSigma * mR);
    return 0.5 * (w1 / w3 - w2 / w3);
#else
    double csR = cosh(mSigma * mR);
    double ssR = sinh(mSigma * mR);
    double csr = cosh(mSigma * rs);
    return csR * csR / (csr * csr + ssR * ssR);
#endif
}

void MetricAlcubierreAccel::calcDF(const double* pos, double &ft, double &fx, double &fy, double &fz) {
    double rs = calcRs(pos);
    double xt = mx0 + mv0 * (pos[0] - mt0);
    double v  = mv0;
    if (pos[0] > mt0) {
        xt += 0.5 * malpha * (pos[0] - mt0) * (pos[0] - mt0);
        v  += malpha * (pos[0] - mt0);
    }

    double thp = tanh(mSigma * (rs + mR));
    double thm = tanh(mSigma * (rs - mR));

    double kl = thm * thm - thp * thp;
    double th = 0.5 * mSigma / (rs * tanh(mSigma * mR));
    double df = kl * th;

    double x = pos[1];
    double y = pos[2];
    double z = pos[3];

//std::cerr << pos[0] << " " << pos[1] << " " << pos[2] << " " << pos[3] << std::endl;
//std::cerr << rs << " " << thp << " " << thm  << " " << df << std::endl;

    ft = -v * (x - xt) * df;
    fx = (x - xt) * df;
    fy = y * df;
    fz = z * df;
}

void MetricAlcubierreAccel::calcD2F(const double* pos, double &ftt, double &ftx, double &fty, double &ftz,
                                    double &fxx, double &fxy, double &fxz, double &fyy,
                                    double &fyz, double &fzz) {
    double t = pos[0];
    double x = pos[1];
    double y = pos[2];
    double z = pos[3];

    double sigma = mSigma;
    double R = mR;

    double xt  = mx0 + mv0 * (t - mt0);
    double dxt = mv0;
    double d2xt = 0.0;
    if (t > mt0) {
        xt += 0.5 * malpha * (t - mt0) * (t - mt0);
        dxt += malpha * (t - mt0);
        d2xt = malpha;
    }

    double t1 = xt;       // x(t);
    double t2 = x - t1;
    double t3 = t2 * t2;
    double t4 = y * y;
    double t5 = z * z;
    double t6 = t3 + t4 + t5;
    double t7 = sqrt(t6);
    double t10 = tanh(sigma * (t7 + R));
    double t11 = t10 * t10;
    double t12 = 1.0 - t11;
    double t13 = t10 * t12;
    double t14 = sigma * sigma;
    double t15 = t13 * t14;
    double t16 = 1 / t6;
    double t18 = dxt;      // diff(x(t),t);
    double t19 = t18 * t18;
    double t20 = t16 * t3 * t19;
    double t23 = t12 * sigma;
    double t25 = 1 / t7 / t6;
    double t27 = t25 * t3 * t19;
    double t29 = 1 / t7;
    double t30 = t29 * t19;
    double t33 = d2xt;     // diff(diff(x(t),t),t);
    double t34 = t29 * t2 * t33;
    double t38 = tanh(sigma * (t7 - R));
    double t39 = t38 * t38;
    double t40 = 1.0 - t39;
    double t41 = t38 * t40;
    double t42 = t41 * t14;
    double t45 = t40 * sigma;
    double t51 = tanh(sigma * R);
    double t52 = 1 / t51;
    double t56 = t18 * t2;
    double t57 = 2.0 * t16 * t2 * t56;
    double t59 = t23 * t25;
    double t60 = 2.0 * t56 * t2;
    double t63 = t29 * t18;
    double t66 = t45 * t25;
    double t73 = t16 * y;
    double t74 = t73 * t56;
    double t77 = t56 * y;
    double t85 = t16 * z;
    double t86 = t85 * t56;
    double t89 = t56 * z;
    double t97 = t14 * t16;
    double t98 = 4.0 * t2 * t2;
    double t99 = t97 * t98;
    double t102 = t25 * t98;
    double t105 = t23 * t29;
    double t110 = t45 * t29;
    double t114 = 2.0 * t73 * t2;
    double t116 = 2.0 * t25 * t2;
    double t117 = t116 * y;
    double t126 = 2.0 * t85 * t2;
    double t128 = t116 * z;
    double t137 = t4 * t97;
    double t140 = t25 * t4;
    double t148 = t85 * y;
    double t152 = t25 * y * z;
    double t160 = t97 * t5;
    double t163 = t25 * t5;

    ftt = (-2.0 * t15 * t20 - t23 * t27 + t23 * t30 - t23 * t34 + 2.0 * t42 * t20 + t45 * t27 - t45 * t30 + t45 * t34) * t52 / 2.0;
    ftx = (t15 * t57 + t59 * t60 / 2.0 - t23 * t63 - t42 * t57 - t66 * t60 / 2.0 + t45 * t63) * t52 / 2.0;
    fty = (2.0 * t15 * t74 + t59 * t77 - 2.0 * t42 * t74 - t66 * t77) * t52 / 2.0;
    ftz = (2.0 * t15 * t86 + t59 * t89 - 2.0 * t42 * t86 - t66 * t89) * t52 / 2.0;
    fxx = (-t13 * t99 / 2.0 - t23 * t102 / 4.0 + t105 + t41 * t99 / 2.0 + t45 * t102 / 4.0 - t110) * t52 / 2.0;
    fxy = (-t15 * t114 - t23 * t117 / 2.0 + t42 * t114 + t45 * t117 / 2.0) * t52 / 2.0;
    fxz = (-t15 * t126 - t23 * t128 / 2.0 + t42 * t126 + t45 * t128 / 2.0) * t52 / 2.0;
    fyy = (-2.0 * t13 * t137 - t23 * t140 + t105 + 2.0 * t41 * t137 + t45 * t140 - t110) * t52 / 2.0;
    fyz = (-2.0 * t15 * t148 - t23 * t152 + 2.0 * t42 * t148 + t45 * t152) * t52 / 2.0;
    fzz = (-2.0 * t13 * t160 - t23 * t163 + t105 + 2.0 * t41 * t160 + t45 * t163 - t110) * t52 / 2.0;
}

} // end namespace m4d
