// -------------------------------------------------------------------------------
/*
   m4dMetricAlcubierre.cpp

  Copyright (c) 2009-2014  Thomas Mueller, Frank Grave


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

#include "m4dMetricAlcubierre.h"

namespace m4d {

#define eps 1.0e-6


/*! Standard constructor for the KerrBL metric.
 *
 * \param  sigma : sharpness of warp bubble.
 * \param  R     : size of warp bubble.
 * \param  vs    : velocity of warp bubble.
 */
MetricAlcubierre::MetricAlcubierre(double sigma, double R, double vs) {
    mMetricName  = "AlcubierreWarp";
    mMetricCPPfilename = "m4dMetricAlcubierre.cpp";
    setCoordType(enum_coordinate_cartesian);

    mPhysicalUnits = enum_physical_constants_geom;
    mSpeedOfLight = 1.0;
    mGravConstant = 1.0;

    mSigma = sigma;
    mR     = R;
    mvs    = vs;

    addParam("sigma", sigma);
    addParam("r", R);
    addParam("vs", vs);

    mDrawTypes.push_back(enum_draw_twoplusone);

    setStandardValues();

    mLocTeds.push_back(enum_nat_tetrad_comoving);
    mLocTeds.push_back(enum_nat_tetrad_static);
}

MetricAlcubierre::~MetricAlcubierre() {
}


// *********************************** public methods ******************************
/*! Calculate the contravariant metric components at position 'pos'.
 *
 *  \param pos : pointer to position.
 */
bool MetricAlcubierre::calculateMetric(const double* pos) {
    double c = mSpeedOfLight;
    double v = mvs;

    double t1 = c * c;
    double t2 = v * v;
    double t3 = calcF(pos);  //f(t,x,y,z);
    double t4 = t3 * t3;
    double t7 = v * t3;

    g_compts[0][0] = -t1 + t2 * t4;
    g_compts[0][1] = -t7;
    g_compts[0][2] = 0.0;
    g_compts[0][3] = 0.0;
    g_compts[1][0] = -t7;
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
bool MetricAlcubierre::calculateChristoffels(const double* pos) {
    double c = mSpeedOfLight;
    double v = mvs;

    double ft, fx, fy, fz;
    calcDF(pos, ft, fx, fy, fz);

    double t1 = v * v;
    double t2 = t1 * v;
    double t3 = calcF(pos);    // f(t,x,y,z);
    double t4 = t3 * t3;
    double t6 = fx;            // diff(f(t,x,y,z),x);
    double t7 = c * c;
    double t8 = 1 / t7;
    double t10 = t2 * t4 * t6 * t8;
    double t11 = ft;           // diff(f(t,x,y,z),t);
    double t14 = t3 * t6;
    double t22 = t1 * t3;
    double t23 = fy;           // diff(f(t,x,y,z),y);
    double t25 = fz;           // diff(f(t,x,y,z),z);
    double t27 = t8 * t1;
    double t28 = t27 * t14;
    double t29 = v * t23;
    double t30 = t29 / 2.0;
    double t31 = v * t25;
    double t32 = t31 / 2.0;
    double t35 = t27 * t3 * t23 / 2.0;
    double t38 = (t1 * t4 + t7) * t8;
    double t40 = t29 * t38 / 2.0;
    double t43 = t27 * t3 * t25 / 2.0;
    double t45 = t31 * t38 / 2.0;
    double t46 = t8 * v;
    double t49 = t46 * t23 / 2.0;
    double t51 = t46 * t25 / 2.0;

    christoffel[0][0][0] = t10;
    christoffel[0][0][1] = v * (-t7 * t11 - t7 * v * t14 + t2 * t4 * t3 * t6) * t8;
    christoffel[0][0][2] = -t22 * t23;
    christoffel[0][0][3] = -t22 * t25;
    christoffel[0][1][0] = -t28;
    christoffel[0][1][1] = -t10;
    christoffel[0][1][2] = t30;
    christoffel[0][1][3] = t32;
    christoffel[0][2][0] = -t35;
    christoffel[0][2][1] = -t40;
    christoffel[0][2][2] = 0.0;
    christoffel[0][2][3] = 0.0;
    christoffel[0][3][0] = -t43;
    christoffel[0][3][1] = -t45;
    christoffel[0][3][2] = 0.0;
    christoffel[0][3][3] = 0.0;
    christoffel[1][0][0] = -t28;
    christoffel[1][0][1] = -t10;
    christoffel[1][0][2] = t30;
    christoffel[1][0][3] = t32;
    christoffel[1][1][0] = t46 * t6;
    christoffel[1][1][1] = t28;
    christoffel[1][1][2] = 0.0;
    christoffel[1][1][3] = 0.0;
    christoffel[1][2][0] = t49;
    christoffel[1][2][1] = t35;
    christoffel[1][2][2] = 0.0;
    christoffel[1][2][3] = 0.0;
    christoffel[1][3][0] = t51;
    christoffel[1][3][1] = t43;
    christoffel[1][3][2] = 0.0;
    christoffel[1][3][3] = 0.0;
    christoffel[2][0][0] = -t35;
    christoffel[2][0][1] = -t40;
    christoffel[2][0][2] = 0.0;
    christoffel[2][0][3] = 0.0;
    christoffel[2][1][0] = t49;
    christoffel[2][1][1] = t35;
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
    christoffel[3][0][0] = -t43;
    christoffel[3][0][1] = -t45;
    christoffel[3][0][2] = 0.0;
    christoffel[3][0][3] = 0.0;
    christoffel[3][1][0] = t51;
    christoffel[3][1][1] = t43;
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
bool MetricAlcubierre::calculateChrisD(const double* pos) {
    double c = mSpeedOfLight;
    double v = mvs;

    double ft, fx, fy, fz;
    calcDF(pos, ft, fx, fy, fz);

    double ftt, ftx, fty, ftz, fxx, fxy, fxz, fyy, fyz, fzz;
    calcD2F(pos, ftt, ftx, fty, ftz, fxx, fxy, fxz, fyy, fyz, fzz);


    double t1 = v * v;
    double t2 = t1 * v;
    double t3 = calcF(pos);        // f(t,x,y,z);
    double t4 = t2 * t3;
    double t5 = fx;                // diff(f(t,x,y,z),x);
    double t6 = ft;                // diff(f(t,x,y,z),t);
    double t7 = t5 * t6;
    double t9 = ftx;               // diff(diff(f(t,x,y,z),t),x);
    double t10 = t3 * t9;
    double t12 = c * c;
    double t13 = 1 / t12;
    double t15 = t4 * (2.0 * t7 + t10) * t13;
    double t16 = t5 * t5;
    double t18 = fxx;              // diff(diff(f(t,x,y,z),x),x);
    double t19 = t3 * t18;
    double t22 = t4 * (2.0 * t16 + t19) * t13;
    double t23 = fy;               // diff(f(t,x,y,z),y);
    double t24 = t5 * t23;
    double t26 = fxy;              // diff(diff(f(t,x,y,z),x),y);
    double t27 = t3 * t26;
    double t30 = t4 * (2.0 * t24 + t27) * t13;
    double t31 = fz;               // diff(f(t,x,y,z),z);
    double t32 = t5 * t31;
    double t34 = fxz;              // diff(diff(f(t,x,y,z),x),z);
    double t35 = t3 * t34;
    double t38 = t4 * (2.0 * t32 + t35) * t13;
    double t39 = t3 * t3;
    double t40 = t2 * t39;
    double t44 = t2 * t39 * t3;
    double t46 = ftt;              // diff(diff(f(t,x,y,z),t),t);
    double t48 = t12 * v;
    double t66 = fty;              // diff(diff(f(t,x,y,z),t),y);
    double t67 = t12 * t66;
    double t76 = ftz;              // diff(diff(f(t,x,y,z),t),z);
    double t77 = t12 * t76;
    double t83 = t1 * t6;
    double t85 = t1 * t3;
    double t88 = t1 * t5;
    double t92 = t23 * t23;
    double t93 = t1 * t92;
    double t94 = fyy;              // diff(diff(f(t,x,y,z),y),y);
    double t97 = t1 * t31;
    double t99 = fyz;              // diff(diff(f(t,x,y,z),y),z);
    double t101 = -t97 * t23 - t85 * t99;
    double t108 = t31 * t31;
    double t109 = t1 * t108;
    double t110 = fzz;             // diff(diff(f(t,x,y,z),z),z);
    double t115 = t1 * (t7 + t10) * t13;
    double t118 = t1 * (t16 + t19) * t13;
    double t121 = t1 * (t24 + t27) * t13;
    double t124 = t1 * (t32 + t35) * t13;
    double t126 = v * t66 / 2.0;
    double t128 = v * t26 / 2.0;
    double t130 = v * t94 / 2.0;
    double t132 = v * t99 / 2.0;
    double t134 = v * t76 / 2.0;
    double t136 = v * t34 / 2.0;
    double t138 = v * t110 / 2.0;
    double t144 = t1 * (t6 * t23 + t3 * t66) * t13 / 2.0;
    double t145 = t121 / 2.0;
    double t150 = t1 * (t92 + t3 * t94) * t13 / 2.0;
    double t156 = t1 * (t31 * t23 + t3 * t99) * t13 / 2.0;
    double t159 = t1 * t23;
    double t160 = t3 * t6;
    double t166 = v * (t66 * t1 * t39 + t67 + 2.0 * t159 * t160) * t13 / 2.0;
    double t175 = v * (t26 * t1 * t39 + t26 * t12 + 2.0 * t85 * t24) * t13 / 2.0;
    double t184 = v * (t94 * t1 * t39 + t94 * t12 + 2.0 * t93 * t3) * t13 / 2.0;
    double t194 = v * (t99 * t1 * t39 + t99 * t12 + 2.0 * t159 * t3 * t31) * t13 / 2.0;
    double t200 = t1 * (t6 * t31 + t3 * t76) * t13 / 2.0;
    double t201 = t124 / 2.0;
    double t206 = t1 * (t108 + t3 * t110) * t13 / 2.0;
    double t214 = v * (t76 * t1 * t39 + t77 + 2.0 * t97 * t160) * t13 / 2.0;
    double t223 = v * (t34 * t1 * t39 + t34 * t12 + 2.0 * t85 * t32) * t13 / 2.0;
    double t232 = v * (t110 * t1 * t39 + t110 * t12 + 2.0 * t109 * t3) * t13 / 2.0;
    double t233 = t13 * v;
    double t236 = t233 * t26;
    double t237 = t233 * t34;
    double t239 = t233 * t66 / 2.0;
    double t240 = t236 / 2.0;
    double t242 = t233 * t94 / 2.0;
    double t244 = t233 * t99 / 2.0;
    double t246 = t233 * t76 / 2.0;
    double t247 = t237 / 2.0;
    double t249 = t233 * t110 / 2.0;

    chrisD[0][0][0][0] = t15;
    chrisD[0][0][0][1] = t22;
    chrisD[0][0][0][2] = t30;
    chrisD[0][0][0][3] = t38;
    chrisD[0][0][1][0] = -v * (-3.0 * t40 * t7 - t44 * t9 + t12 * t46 + t48 * t7 + t48 * t10) * t13;
    chrisD[0][0][1][1] = -v * (-3.0 * t40 * t16 - t44 * t18 + t12 * t9 + t48 * t16 + t48 * t19) * t13;
    chrisD[0][0][1][2] = -v * (-3.0 * t40 * t24 - t44 * t26 + t67 + t48 * t24 + t48 * t27) * t13;
    chrisD[0][0][1][3] = -v * (-3.0 * t40 * t32 - t44 * t34 + t77 + t48 * t32 + t48 * t35) * t13;
    chrisD[0][0][2][0] = -t83 * t23 - t85 * t66;
    chrisD[0][0][2][1] = -t88 * t23 - t85 * t26;
    chrisD[0][0][2][2] = -t93 - t85 * t94;
    chrisD[0][0][2][3] = t101;
    chrisD[0][0][3][0] = -t83 * t31 - t85 * t76;
    chrisD[0][0][3][1] = -t88 * t31 - t85 * t34;
    chrisD[0][0][3][2] = t101;
    chrisD[0][0][3][3] = -t109 - t85 * t110;
    chrisD[0][1][0][0] = -t115;
    chrisD[0][1][0][1] = -t118;
    chrisD[0][1][0][2] = -t121;
    chrisD[0][1][0][3] = -t124;
    chrisD[0][1][1][0] = -t15;
    chrisD[0][1][1][1] = -t22;
    chrisD[0][1][1][2] = -t30;
    chrisD[0][1][1][3] = -t38;
    chrisD[0][1][2][0] = t126;
    chrisD[0][1][2][1] = t128;
    chrisD[0][1][2][2] = t130;
    chrisD[0][1][2][3] = t132;
    chrisD[0][1][3][0] = t134;
    chrisD[0][1][3][1] = t136;
    chrisD[0][1][3][2] = t132;
    chrisD[0][1][3][3] = t138;
    chrisD[0][2][0][0] = -t144;
    chrisD[0][2][0][1] = -t145;
    chrisD[0][2][0][2] = -t150;
    chrisD[0][2][0][3] = -t156;
    chrisD[0][2][1][0] = -t166;
    chrisD[0][2][1][1] = -t175;
    chrisD[0][2][1][2] = -t184;
    chrisD[0][2][1][3] = -t194;
    chrisD[0][2][2][0] = 0.0;
    chrisD[0][2][2][1] = 0.0;
    chrisD[0][2][2][2] = 0.0;
    chrisD[0][2][2][3] = 0.0;
    chrisD[0][2][3][0] = 0.0;
    chrisD[0][2][3][1] = 0.0;
    chrisD[0][2][3][2] = 0.0;
    chrisD[0][2][3][3] = 0.0;
    chrisD[0][3][0][0] = -t200;
    chrisD[0][3][0][1] = -t201;
    chrisD[0][3][0][2] = -t156;
    chrisD[0][3][0][3] = -t206;
    chrisD[0][3][1][0] = -t214;
    chrisD[0][3][1][1] = -t223;
    chrisD[0][3][1][2] = -t194;
    chrisD[0][3][1][3] = -t232;
    chrisD[0][3][2][0] = 0.0;
    chrisD[0][3][2][1] = 0.0;
    chrisD[0][3][2][2] = 0.0;
    chrisD[0][3][2][3] = 0.0;
    chrisD[0][3][3][0] = 0.0;
    chrisD[0][3][3][1] = 0.0;
    chrisD[0][3][3][2] = 0.0;
    chrisD[0][3][3][3] = 0.0;
    chrisD[1][0][0][0] = -t115;
    chrisD[1][0][0][1] = -t118;
    chrisD[1][0][0][2] = -t121;
    chrisD[1][0][0][3] = -t124;
    chrisD[1][0][1][0] = -t15;
    chrisD[1][0][1][1] = -t22;
    chrisD[1][0][1][2] = -t30;
    chrisD[1][0][1][3] = -t38;
    chrisD[1][0][2][0] = t126;
    chrisD[1][0][2][1] = t128;
    chrisD[1][0][2][2] = t130;
    chrisD[1][0][2][3] = t132;
    chrisD[1][0][3][0] = t134;
    chrisD[1][0][3][1] = t136;
    chrisD[1][0][3][2] = t132;
    chrisD[1][0][3][3] = t138;
    chrisD[1][1][0][0] = t233 * t9;
    chrisD[1][1][0][1] = t233 * t18;
    chrisD[1][1][0][2] = t236;
    chrisD[1][1][0][3] = t237;
    chrisD[1][1][1][0] = t115;
    chrisD[1][1][1][1] = t118;
    chrisD[1][1][1][2] = t121;
    chrisD[1][1][1][3] = t124;
    chrisD[1][1][2][0] = 0.0;
    chrisD[1][1][2][1] = 0.0;
    chrisD[1][1][2][2] = 0.0;
    chrisD[1][1][2][3] = 0.0;
    chrisD[1][1][3][0] = 0.0;
    chrisD[1][1][3][1] = 0.0;
    chrisD[1][1][3][2] = 0.0;
    chrisD[1][1][3][3] = 0.0;
    chrisD[1][2][0][0] = t239;
    chrisD[1][2][0][1] = t240;
    chrisD[1][2][0][2] = t242;
    chrisD[1][2][0][3] = t244;
    chrisD[1][2][1][0] = t144;
    chrisD[1][2][1][1] = t145;
    chrisD[1][2][1][2] = t150;
    chrisD[1][2][1][3] = t156;
    chrisD[1][2][2][0] = 0.0;
    chrisD[1][2][2][1] = 0.0;
    chrisD[1][2][2][2] = 0.0;
    chrisD[1][2][2][3] = 0.0;
    chrisD[1][2][3][0] = 0.0;
    chrisD[1][2][3][1] = 0.0;
    chrisD[1][2][3][2] = 0.0;
    chrisD[1][2][3][3] = 0.0;
    chrisD[1][3][0][0] = t246;
    chrisD[1][3][0][1] = t247;
    chrisD[1][3][0][2] = t244;
    chrisD[1][3][0][3] = t249;
    chrisD[1][3][1][0] = t200;
    chrisD[1][3][1][1] = t201;
    chrisD[1][3][1][2] = t156;
    chrisD[1][3][1][3] = t206;
    chrisD[1][3][2][0] = 0.0;
    chrisD[1][3][2][1] = 0.0;
    chrisD[1][3][2][2] = 0.0;
    chrisD[1][3][2][3] = 0.0;
    chrisD[1][3][3][0] = 0.0;
    chrisD[1][3][3][1] = 0.0;
    chrisD[1][3][3][2] = 0.0;
    chrisD[1][3][3][3] = 0.0;
    chrisD[2][0][0][0] = -t144;
    chrisD[2][0][0][1] = -t145;
    chrisD[2][0][0][2] = -t150;
    chrisD[2][0][0][3] = -t156;
    chrisD[2][0][1][0] = -t166;
    chrisD[2][0][1][1] = -t175;
    chrisD[2][0][1][2] = -t184;
    chrisD[2][0][1][3] = -t194;
    chrisD[2][0][2][0] = 0.0;
    chrisD[2][0][2][1] = 0.0;
    chrisD[2][0][2][2] = 0.0;
    chrisD[2][0][2][3] = 0.0;
    chrisD[2][0][3][0] = 0.0;
    chrisD[2][0][3][1] = 0.0;
    chrisD[2][0][3][2] = 0.0;
    chrisD[2][0][3][3] = 0.0;
    chrisD[2][1][0][0] = t239;
    chrisD[2][1][0][1] = t240;
    chrisD[2][1][0][2] = t242;
    chrisD[2][1][0][3] = t244;
    chrisD[2][1][1][0] = t144;
    chrisD[2][1][1][1] = t145;
    chrisD[2][1][1][2] = t150;
    chrisD[2][1][1][3] = t156;
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
    chrisD[3][0][0][0] = -t200;
    chrisD[3][0][0][1] = -t201;
    chrisD[3][0][0][2] = -t156;
    chrisD[3][0][0][3] = -t206;
    chrisD[3][0][1][0] = -t214;
    chrisD[3][0][1][1] = -t223;
    chrisD[3][0][1][2] = -t194;
    chrisD[3][0][1][3] = -t232;
    chrisD[3][0][2][0] = 0.0;
    chrisD[3][0][2][1] = 0.0;
    chrisD[3][0][2][2] = 0.0;
    chrisD[3][0][2][3] = 0.0;
    chrisD[3][0][3][0] = 0.0;
    chrisD[3][0][3][1] = 0.0;
    chrisD[3][0][3][2] = 0.0;
    chrisD[3][0][3][3] = 0.0;
    chrisD[3][1][0][0] = t246;
    chrisD[3][1][0][1] = t247;
    chrisD[3][1][0][2] = t244;
    chrisD[3][1][0][3] = t249;
    chrisD[3][1][1][0] = t200;
    chrisD[3][1][1][1] = t201;
    chrisD[3][1][1][2] = t156;
    chrisD[3][1][1][3] = t206;
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
bool MetricAlcubierre::calculateRiemann(const double* pos) {
    double c = mSpeedOfLight;
    double v = mvs;

    double ft, fx, fy, fz;
    calcDF(pos, ft, fx, fy, fz);

    double ftt, ftx, fty, ftz, fxx, fxy, fxz, fyy, fyz, fzz;
    calcD2F(pos, ftt, ftx, fty, ftz, fxx, fxy, fxz, fyy, fyz, fzz);

    double t1 = v * v;
    double t2 = calcF(pos);           // f(t,x,y,z);
    double t3 = t1 * t2;
    double t4 = c * c;
    double t5 = 1 / t4;
    double t6 = ftx;                  // diff(diff(f(t,x,y,z),t),x);
    double t8 = fx;                   // diff(f(t,x,y,z),x);
    double t9 = t8 * t8;
    double t12 = v * t2;
    double t13 = fxx;                 // diff(diff(f(t,x,y,z),x),x);
    double t16 = fy;                  // diff(f(t,x,y,z),y);
    double t17 = t16 * t16;
    double t19 = fz;                  // diff(f(t,x,y,z),z);
    double t20 = t19 * t19;
    double t22 = 4.0 * t6 + 4.0 * v * t9 + 4.0 * t12 * t13 - v * t17 - v * t20;
    double t25 = t3 * t5 * t22 / 4.0;
    double t26 = fty;                 // diff(diff(f(t,x,y,z),t),y);
    double t29 = 2.0 * v * t16 * t8;
    double t30 = fxy;                 // diff(diff(f(t,x,y,z),x),y);
    double t31 = t12 * t30;
    double t33 = t26 + t29 + 2.0 * t31;
    double t36 = t3 * t5 * t33 / 2.0;
    double t37 = ftz;                 // diff(diff(f(t,x,y,z),t),z);
    double t40 = 2.0 * v * t19 * t8;
    double t41 = fxz;                 // diff(diff(f(t,x,y,z),x),z);
    double t42 = t12 * t41;
    double t44 = t37 + t40 + 2.0 * t42;
    double t47 = t3 * t5 * t44 / 2.0;
    double t50 = t3 * t5 * t30 / 2.0;
    double t53 = t3 * t5 * t41 / 2.0;
    double t54 = t5 * v;
    double t56 = t54 * t22 / 4.0;
    double t58 = t54 * t33 / 2.0;
    double t60 = t54 * t44 / 2.0;
    double t62 = t54 * t30 / 2.0;
    double t64 = t54 * t41 / 2.0;
    double t68 = v * (t26 + t29 + t31) * t5 / 2.0;
    double t69 = t1 * t5;
    double t71 = fyy;                 // diff(diff(f(t,x,y,z),y),y);
    double t72 = t2 * t71;
    double t73 = 2.0 * t72;
    double t76 = t69 * (3.0 * t17 + t73) / 4.0;
    double t77 = t16 * t19;
    double t79 = fyz;                 // diff(diff(f(t,x,y,z),y),z);
    double t80 = t2 * t79;
    double t81 = 2.0 * t80;
    double t84 = t69 * (3.0 * t77 + t81) / 4.0;
    double t87 = v * t71 * t5 / 2.0;
    double t90 = v * t79 * t5 / 2.0;
    double t94 = v * (t37 + t40 + t42) * t5 / 2.0;
    double t96 = fzz;                 // diff(diff(f(t,x,y,z),z),z);
    double t97 = t2 * t96;
    double t98 = 2.0 * t97;
    double t101 = t69 * (3.0 * t20 + t98) / 4.0;
    double t104 = v * t96 * t5 / 2.0;
    double t105 = t2 * t2;
    double t106 = t1 * t105;
    double t108 = (t4 - t106) * t5;
    double t111 = t108 * v * t22 / 4.0;
    double t112 = v * t33;
    double t114 = t108 * t112 / 2.0;
    double t115 = v * t44;
    double t117 = t108 * t115 / 2.0;
    double t118 = v * t30;
    double t120 = t108 * t118 / 2.0;
    double t121 = v * t41;
    double t123 = t108 * t121 / 2.0;
    double t133 = v * (t12 * t26 + 2.0 * t3 * t16 * t8 + t106 * t30 + t30 * t4) * t5 / 2.0;
    double t134 = t3 * t17;
    double t137 = t71 * t4;
    double t140 = t54 * (2.0 * t134 + t106 * t71 + t137) / 2.0;
    double t141 = t3 * t77;
    double t144 = t79 * t4;
    double t147 = t54 * (2.0 * t141 + t106 * t79 + t144) / 2.0;
    double t150 = t69 * (t73 + t17) / 4.0;
    double t153 = t69 * (t81 + t77) / 4.0;
    double t163 = v * (t12 * t37 + 2.0 * t3 * t19 * t8 + t106 * t41 + t41 * t4) * t5 / 2.0;
    double t164 = t3 * t20;
    double t167 = t96 * t4;
    double t170 = t54 * (2.0 * t164 + t106 * t96 + t167) / 2.0;
    double t173 = t69 * (t98 + t20) / 4.0;
    double t174 = t112 / 2.0;
    double t184 = t1 * (3.0 * t17 * t4 + 4.0 * t72 * t4 + t17 * t1 * t105) * t5 / 4.0;
    double t193 = t1 * (3.0 * t77 * t4 + 4.0 * t80 * t4 + t77 * t106) * t5 / 4.0;
    double t198 = v * (2.0 * t137 + t134) * t5 / 4.0;
    double t203 = v * (2.0 * t144 + t141) * t5 / 4.0;
    double t204 = t118 / 2.0;
    double t206 = t69 * t17 / 4.0;
    double t208 = t69 * t77 / 4.0;
    double t209 = t115 / 2.0;
    double t219 = t1 * (3.0 * t20 * t4 + 4.0 * t97 * t4 + t20 * t1 * t105) * t5 / 4.0;
    double t224 = v * (2.0 * t167 + t164) * t5 / 4.0;
    double t225 = t121 / 2.0;
    double t227 = t69 * t20 / 4.0;

    riem[0][0][0][0] = 0.0;
    riem[0][0][0][1] = -t25;
    riem[0][0][0][2] = -t36;
    riem[0][0][0][3] = -t47;
    riem[0][0][1][0] = t25;
    riem[0][0][1][1] = 0.0;
    riem[0][0][1][2] = t50;
    riem[0][0][1][3] = t53;
    riem[0][0][2][0] = t36;
    riem[0][0][2][1] = -t50;
    riem[0][0][2][2] = 0.0;
    riem[0][0][2][3] = 0.0;
    riem[0][0][3][0] = t47;
    riem[0][0][3][1] = -t53;
    riem[0][0][3][2] = 0.0;
    riem[0][0][3][3] = 0.0;
    riem[0][1][0][0] = 0.0;
    riem[0][1][0][1] = t56;
    riem[0][1][0][2] = t58;
    riem[0][1][0][3] = t60;
    riem[0][1][1][0] = -t56;
    riem[0][1][1][1] = 0.0;
    riem[0][1][1][2] = -t62;
    riem[0][1][1][3] = -t64;
    riem[0][1][2][0] = -t58;
    riem[0][1][2][1] = t62;
    riem[0][1][2][2] = 0.0;
    riem[0][1][2][3] = 0.0;
    riem[0][1][3][0] = -t60;
    riem[0][1][3][1] = t64;
    riem[0][1][3][2] = 0.0;
    riem[0][1][3][3] = 0.0;
    riem[0][2][0][0] = 0.0;
    riem[0][2][0][1] = t68;
    riem[0][2][0][2] = t76;
    riem[0][2][0][3] = t84;
    riem[0][2][1][0] = -t68;
    riem[0][2][1][1] = 0.0;
    riem[0][2][1][2] = -t87;
    riem[0][2][1][3] = -t90;
    riem[0][2][2][0] = -t76;
    riem[0][2][2][1] = t87;
    riem[0][2][2][2] = 0.0;
    riem[0][2][2][3] = 0.0;
    riem[0][2][3][0] = -t84;
    riem[0][2][3][1] = t90;
    riem[0][2][3][2] = 0.0;
    riem[0][2][3][3] = 0.0;
    riem[0][3][0][0] = 0.0;
    riem[0][3][0][1] = t94;
    riem[0][3][0][2] = t84;
    riem[0][3][0][3] = t101;
    riem[0][3][1][0] = -t94;
    riem[0][3][1][1] = 0.0;
    riem[0][3][1][2] = -t90;
    riem[0][3][1][3] = -t104;
    riem[0][3][2][0] = -t84;
    riem[0][3][2][1] = t90;
    riem[0][3][2][2] = 0.0;
    riem[0][3][2][3] = 0.0;
    riem[0][3][3][0] = -t101;
    riem[0][3][3][1] = t104;
    riem[0][3][3][2] = 0.0;
    riem[0][3][3][3] = 0.0;
    riem[1][0][0][0] = 0.0;
    riem[1][0][0][1] = t111;
    riem[1][0][0][2] = t114;
    riem[1][0][0][3] = t117;
    riem[1][0][1][0] = -t111;
    riem[1][0][1][1] = 0.0;
    riem[1][0][1][2] = -t120;
    riem[1][0][1][3] = -t123;
    riem[1][0][2][0] = -t114;
    riem[1][0][2][1] = t120;
    riem[1][0][2][2] = 0.0;
    riem[1][0][2][3] = 0.0;
    riem[1][0][3][0] = -t117;
    riem[1][0][3][1] = t123;
    riem[1][0][3][2] = 0.0;
    riem[1][0][3][3] = 0.0;
    riem[1][1][0][0] = 0.0;
    riem[1][1][0][1] = t25;
    riem[1][1][0][2] = t36;
    riem[1][1][0][3] = t47;
    riem[1][1][1][0] = -t25;
    riem[1][1][1][1] = 0.0;
    riem[1][1][1][2] = -t50;
    riem[1][1][1][3] = -t53;
    riem[1][1][2][0] = -t36;
    riem[1][1][2][1] = t50;
    riem[1][1][2][2] = 0.0;
    riem[1][1][2][3] = 0.0;
    riem[1][1][3][0] = -t47;
    riem[1][1][3][1] = t53;
    riem[1][1][3][2] = 0.0;
    riem[1][1][3][3] = 0.0;
    riem[1][2][0][0] = 0.0;
    riem[1][2][0][1] = t133;
    riem[1][2][0][2] = t140;
    riem[1][2][0][3] = t147;
    riem[1][2][1][0] = -t133;
    riem[1][2][1][1] = 0.0;
    riem[1][2][1][2] = -t150;
    riem[1][2][1][3] = -t153;
    riem[1][2][2][0] = -t140;
    riem[1][2][2][1] = t150;
    riem[1][2][2][2] = 0.0;
    riem[1][2][2][3] = 0.0;
    riem[1][2][3][0] = -t147;
    riem[1][2][3][1] = t153;
    riem[1][2][3][2] = 0.0;
    riem[1][2][3][3] = 0.0;
    riem[1][3][0][0] = 0.0;
    riem[1][3][0][1] = t163;
    riem[1][3][0][2] = t147;
    riem[1][3][0][3] = t170;
    riem[1][3][1][0] = -t163;
    riem[1][3][1][1] = 0.0;
    riem[1][3][1][2] = -t153;
    riem[1][3][1][3] = -t173;
    riem[1][3][2][0] = -t147;
    riem[1][3][2][1] = t153;
    riem[1][3][2][2] = 0.0;
    riem[1][3][2][3] = 0.0;
    riem[1][3][3][0] = -t170;
    riem[1][3][3][1] = t173;
    riem[1][3][3][2] = 0.0;
    riem[1][3][3][3] = 0.0;
    riem[2][0][0][0] = 0.0;
    riem[2][0][0][1] = t174;
    riem[2][0][0][2] = t184;
    riem[2][0][0][3] = t193;
    riem[2][0][1][0] = -t174;
    riem[2][0][1][1] = 0.0;
    riem[2][0][1][2] = -t198;
    riem[2][0][1][3] = -t203;
    riem[2][0][2][0] = -t184;
    riem[2][0][2][1] = t198;
    riem[2][0][2][2] = 0.0;
    riem[2][0][2][3] = 0.0;
    riem[2][0][3][0] = -t193;
    riem[2][0][3][1] = t203;
    riem[2][0][3][2] = 0.0;
    riem[2][0][3][3] = 0.0;
    riem[2][1][0][0] = 0.0;
    riem[2][1][0][1] = -t204;
    riem[2][1][0][2] = -t198;
    riem[2][1][0][3] = -t203;
    riem[2][1][1][0] = t204;
    riem[2][1][1][1] = 0.0;
    riem[2][1][1][2] = t206;
    riem[2][1][1][3] = t208;
    riem[2][1][2][0] = t198;
    riem[2][1][2][1] = -t206;
    riem[2][1][2][2] = 0.0;
    riem[2][1][2][3] = 0.0;
    riem[2][1][3][0] = t203;
    riem[2][1][3][1] = -t208;
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
    riem[3][0][0][1] = t209;
    riem[3][0][0][2] = t193;
    riem[3][0][0][3] = t219;
    riem[3][0][1][0] = -t209;
    riem[3][0][1][1] = 0.0;
    riem[3][0][1][2] = -t203;
    riem[3][0][1][3] = -t224;
    riem[3][0][2][0] = -t193;
    riem[3][0][2][1] = t203;
    riem[3][0][2][2] = 0.0;
    riem[3][0][2][3] = 0.0;
    riem[3][0][3][0] = -t219;
    riem[3][0][3][1] = t224;
    riem[3][0][3][2] = 0.0;
    riem[3][0][3][3] = 0.0;
    riem[3][1][0][0] = 0.0;
    riem[3][1][0][1] = -t225;
    riem[3][1][0][2] = -t203;
    riem[3][1][0][3] = -t224;
    riem[3][1][1][0] = t225;
    riem[3][1][1][1] = 0.0;
    riem[3][1][1][2] = t208;
    riem[3][1][1][3] = t227;
    riem[3][1][2][0] = t203;
    riem[3][1][2][1] = -t208;
    riem[3][1][2][2] = 0.0;
    riem[3][1][2][3] = 0.0;
    riem[3][1][3][0] = t224;
    riem[3][1][3][1] = -t227;
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
void MetricAlcubierre::localToCoord(const double* pos, const double* ldir, double* dir,
                                    enum_nat_tetrad_type  type) {
    double f = calcF(pos);
    double c = mSpeedOfLight;

    if (type == enum_nat_tetrad_comoving) {
        dir[0] = ldir[0] / c;
        dir[1] = ldir[0] * mvs * f / c + ldir[1];
        dir[2] = ldir[2];
        dir[3] = ldir[3];
    } else {
        double w = sqrt(c * c - mvs * mvs * f * f);

        dir[0] = (ldir[0] - mvs * f / c * ldir[1]) / w;
        dir[1] = w / c * ldir[1];
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
void MetricAlcubierre::coordToLocal(const double* pos, const double* cdir, double* ldir,
                                    enum_nat_tetrad_type  type) {
    double f = calcF(pos);
    double c = mSpeedOfLight;

    if (type == enum_nat_tetrad_comoving) {
        ldir[0] = c * cdir[0];
        ldir[1] = cdir[1] - mvs * f * cdir[0];
        ldir[2] = cdir[2];
        ldir[3] = cdir[3];
    } else {
        double w = sqrt(c * c - mvs * mvs * f * f);

        ldir[1] = c / w * cdir[1];
        ldir[0] = w * cdir[0] - ldir[1] * mvs * f / c;
        ldir[2] = cdir[2];
        ldir[3] = cdir[3];
    }
}


/*!
 *  \param pos  :  position.
 *  \return true  : radial position r < 0.0 or ...
 *  \return false : position is valid.
 */
bool MetricAlcubierre::breakCondition(const double*) {
    bool br = false;

    return br;
}


/*! Calculate right hand side of the geodesic equation in first order form.
 *
 *  \param  y[]   : pointer to position and direction coordinates.
 *  \param  dydx[] : pointer to right side of geodesic equation.
 */
bool MetricAlcubierre::calcDerivs(const double y[], double dydx[]) {
    dydx[0] = y[4];
    dydx[1] = y[5];
    dydx[2] = y[6];
    dydx[3] = y[7];

    double f = calcF(y);
    double f2 = f * f;

    double c = mSpeedOfLight;
    double edc2 = 1.0 / (c * c);
    //fprintf(stderr,"hier...%f %f\n",y[5]-mvs*y[4],f);
    double v = mvs;
    double v2 = v * v;
    double v3 = v2 * v;

    double ft, fx, fy, fz;
    calcDF(y, ft, fx, fy, fz);

    dydx[4] = edc2 * (- v3 * f * f * fx * y[4] * y[4] + 2.0 * v2 * f * fx * y[4] * y[5] + v2 * f * fy * y[4] * y[6] + v2 * f * fz * y[4] * y[7]
                      - v * fx * y[5] * y[5] - v * fy * y[5] * y[6] - v * fz * y[5] * y[7]);
    dydx[5] = v * ft * y[4] * y[4] + v2 * f * fx * y[4] * y[4] + v * fy * y[4] * y[6] + v * fz * y[4] * y[7]
              + edc2 * (-v2 * v2 * f2 * f * fx * y[4] * y[4] + 2.0 * v3 * f2 * fx * y[4] * y[5] + v3 * f2 * fy * y[4] * y[6]
                        + v3 * f2 * fz * y[4] * y[7] - v2 * f * fx * y[5] * y[5] - v2 * f * fy * y[5] * y[6] - v2 * f * fz * y[5] * y[7]);
    dydx[6] = v2 * f * fy * y[4] * y[4] - v * fy * y[4] * y[5];
    dydx[7] = v2 * f * fz * y[4] * y[4] - v * fz * y[4] * y[5];
    return true;
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
double MetricAlcubierre::testConstraint(const double y[], const double kappa) {
    double c = mSpeedOfLight;
    double f = calcF(y);

    double sum = -kappa;
    sum += -c * c * y[4] * y[4] + pow(y[5] - mvs * f * y[4], 2.0) + y[6] * y[6] + y[7] * y[7];

    //double A = 1.0-mvs*mvs*(1.0-f)*(1.0-f);
    //fprintf(stderr,"%e %e %e %e  %e %e %e %e\n",y[4],y[5],y[6],y[7],f,A,mvs*f+1,mvs*f-1);
    return sum;
}

/*! Resize the tangent vector.
 *
 * \return true: have resized.
 * \return false: have not resized.
 */
bool MetricAlcubierre::resize(double* y, double kappa, double factor) {
    double f = calcF(y);
    double z = mvs * f;

    if (kappa < -0.5) {
        double c2 = y[6] * y[6] + y[7] * y[7];
        //fprintf(stderr,"v: %e %e %e %e \n",y[4],y[5],y[6],y[7]);
        if (y[4] > 0.0) {
            y[4] = -(sqrt((kappa - c2) * z * z - kappa + y[5] * y[5] + c2) - y[5] * z) / (z * z - 1.0);
        } else {
            y[4] = (sqrt((kappa - c2) * z * z - kappa + y[5] * y[5] + c2) + y[5] * z) / (z * z - 1.0);
        }
        //fprintf(stderr,"n: %e %e %e %e \n",y[4],y[5],y[6],y[7]);
    } else {
        //  factor = 0.1;
        y[5] *= factor;
        y[6] *= factor;
        y[7] *= factor;
        double c2 = y[6] * y[6] + y[7] * y[7];
        //  fprintf(stderr,"v:%12.8f\n",y[4]);
        if (y[4] > 0.0) {
            y[4] = -(sqrt(c2 + y[5] * y[5] - c2 * z * z) - y[5] * z) / (z * z - 1.0);
        } else {
            y[4] = (sqrt(c2 + y[5] * y[5] - c2 * z * z) + y[5] * z) / (z * z - 1.0);
        }
        // fprintf(stderr,"n:%12.8f %12.8e\n",y[4],sum);
    }
    return true;
}

/*! Set parameter 'pName' to 'val'.
 *
 *  Set 'sigma', 'R', and 'vs' parameters.
 */
bool MetricAlcubierre::setParam(std::string pName, double val) {
    Metric::setParam(pName, val);
    if (pName == "sigma") {
        mSigma = val;
    } else if (pName == "r") {
        mR = val;
    } else if (pName == "vs") {
        mvs = val;
    }
    return true;
}

/*! Transform point p to 2+1 coordinates.
 *
 *  \param  p  : point in proper metric coordinates.
 *  \param  cp : reference to transformed point.
 *  \return true : success.
 */
bool MetricAlcubierre::transToTwoPlusOne(vec4 p, vec4 &cp) {
    cp = vec4(p[0], p[1], p[2], p[0]);
    return true;
}

/*! Generate report.
 */
bool MetricAlcubierre::report(const vec4 , const vec4 , std::string &text) {
    std::stringstream ss;
    ss << "Report for AlcubierreWarp metric\n\tcoordinate : (t,x,y,z)\n";
    ss << "---------------------------------------------------------------\n";
    ss << "  physical units ................................. no\n";
    ss.precision(DEF_FIXED_REPORT_PRECISION);
    ss.setf(std::ios::fixed);


    text = ss.str();
    return true;
}

// ********************************* protected methods *****************************
/*!
 */
void MetricAlcubierre::setStandardValues() {
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
double MetricAlcubierre::calcRs(const double* pos) {
    double t = pos[0];
    double x = pos[1];
    double y = pos[2];
    double z = pos[3];

    return sqrt((x - mvs * t) * (x - mvs * t) + y * y + z * z);
}

double MetricAlcubierre::calcF(const double* pos) {
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

void MetricAlcubierre::calcDF(const double* pos, double &ft, double &fx, double &fy, double &fz) {
    double rs = calcRs(pos);
#if 1
    double thp = tanh(mSigma * (rs + mR));
    double thm = tanh(mSigma * (rs - mR));

    double kl = thm * thm - thp * thp;
    double th = 0.5 * mSigma / (rs * tanh(mSigma * mR));
    double df = kl * th;

    double t = pos[0];
    double x = pos[1];
    double y = pos[2];
    double z = pos[3];

    //std::cerr << pos[0] << " " << pos[1] << " " << pos[2] << " " << pos[3] << std::endl;
    //std::cerr << rs << " " << thp << " " << thm  << " " << df << std::endl;

    ft = -mvs * (x - mvs * t) * df;
    fx = (x - mvs * t) * df;
    fy = y * df;
    fz = z * df;
#else
    double t = pos[0];
    double x = pos[1];
    double y = pos[2];
    double z = pos[3];

    double csR = cosh(mSigma * mR);
    double ssR = sinh(mSigma * mR);
    double csr = cosh(mSigma * rs);
    double f = csR * csR / (csr * csr + ssR * ssR);
    double dfdr = -mSigma * sinh(2.0 * mSigma * rs) * f * f / (csR * csR);

    double drdt = -mvs * (x - mvs * t) / rs;
    double drdx = (x - mvs * t) / rs;
    double drdy = y / rs;
    double drdz = z / rs;

    ft = dfdr * drdt;
    fx = dfdr * drdx;
    fy = dfdr * drdy;
    fz = dfdr * drdz;
#endif
}

void MetricAlcubierre::calcD2F(const double* pos, double &ftt, double &ftx, double &fty, double &ftz,
                               double &fxx, double &fxy, double &fxz, double &fyy,
                               double &fyz, double &fzz) {
    double t = pos[0];
    double x = pos[1];
    double y = pos[2];
    double z = pos[3];

    double vs = mvs;
    double sigma = mSigma;
    double R = mR;

    double t1 = x * x;
    double t2 = x * vs;
    double t5 = vs * vs;
    double t6 = t * t;
    double t8 = y * y;
    double t9 = z * z;
    double t10 = t1 - 2.0 * t2 * t + t5 * t6 + t8 + t9;
    double t11 = sqrt(t10);
    double t14 = tanh(sigma * (t11 + R));
    double t15 = t14 * t14;
    double t16 = 1.0 - t15;
    double t17 = t14 * t16;
    double t18 = sigma * sigma;
    double t19 = 1 / t10;
    double t20 = t18 * t19;
    double t22 = -t2 + t5 * t;
    double t23 = 4.0 * t22 * t22;
    double t24 = t20 * t23;
    double t27 = t16 * sigma;
    double t29 = 1 / t11 / t10;
    double t30 = t29 * t23;
    double t33 = 1 / t11;
    double t34 = t33 * t5;
    double t38 = tanh(sigma * (t11 - R));
    double t39 = t38 * t38;
    double t40 = 1.0 - t39;
    double t41 = t38 * t40;
    double t44 = t40 * sigma;
    double t50 = tanh(sigma * R);
    double t51 = 1 / t50;
    double t54 = t17 * t18;
    double t56 = x - vs * t;
    double t58 = 4.0 * t19 * t56 * t22;
    double t61 = 2.0 * t29 * t22;
    double t62 = 2.0 * t61 * t56;
    double t65 = t33 * vs;
    double t67 = t41 * t18;
    double t76 = t19 * y;
    double t77 = 2.0 * t76 * t22;
    double t79 = t61 * y;
    double t88 = t19 * z;
    double t89 = 2.0 * t88 * t22;
    double t91 = t61 * z;
    double t100 = 4.0 * t56 * t56;
    double t101 = t20 * t100;
    double t104 = t29 * t100;
    double t107 = t27 * t33;
    double t112 = t44 * t33;
    double t116 = 2.0 * t76 * t56;
    double t118 = 2.0 * t29 * t56;
    double t119 = t118 * y;
    double t128 = 2.0 * t88 * t56;
    double t130 = t118 * z;
    double t139 = t20 * t8;
    double t142 = t29 * t8;
    double t150 = t88 * y;
    double t154 = t29 * y * z;
    double t162 = t20 * t9;
    double t165 = t29 * t9;

    ftt = (-t17 * t24 / 2.0 - t27 * t30 / 4.0 + t27 * t34 + t41 * t24 / 2.0 + t44 * t30 / 4.0 - t44 * t34) * t51 / 2.0;
    ftx = (-t54 * t58 / 2.0 - t27 * t62 / 4.0 - t27 * t65 + t67 * t58 / 2.0 + t44 * t62 / 4.0 + t44 * t65) * t51 / 2.0;
    fty = (-t54 * t77 - t27 * t79 / 2.0 + t67 * t77 + t44 * t79 / 2.0) * t51 / 2.0;
    ftz = (-t54 * t89 - t27 * t91 / 2.0 + t67 * t89 + t44 * t91 / 2.0) * t51 / 2.0;
    fxx = (-t17 * t101 / 2.0 - t27 * t104 / 4.0 + t107 + t41 * t101 / 2.0 + t44 * t104 / 4.0 - t112) * t51 / 2.0;
    fxy = (-t54 * t116 - t27 * t119 / 2.0 + t67 * t116 + t44 * t119 / 2.0) * t51 / 2.0;
    fxz = (-t54 * t128 - t27 * t130 / 2.0 + t67 * t128 + t44 * t130 / 2.0) * t51 / 2.0;
    fyy = (-2.0 * t17 * t139 - t27 * t142 + t107 + 2.0 * t41 * t139 + t44 * t142 - t112) * t51 / 2.0;
    fyz = (-2.0 * t54 * t150 - t27 * t154 + 2.0 * t67 * t150 + t44 * t154) * t51 / 2.0;
    fzz = (-2.0 * t17 * t162 - t27 * t165 + t107 + 2.0 * t41 * t162 + t44 * t165 - t112) * t51 / 2.0;
}

} // end namespace m4d
