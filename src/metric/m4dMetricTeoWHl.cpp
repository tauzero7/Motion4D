// -------------------------------------------------------------------------------
/*
   m4dMetricTeoWHl.cpp

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

#include "m4dMetricTeoWHl.h"

namespace m4d {

#define eps 1.0e-6


/*! Standard constructor for the TeoWHl metric.
 *
 * \param  b0 : throat size.
 */
MetricTeoWHl::MetricTeoWHl(double b0) {
    mMetricName  = "TeoWHl";
    setCoordType(enum_coordinate_spherical);

    mPhysicalUnits = enum_physical_constants_geom;
    mSpeedOfLight = 1.0;
    mGravConstant = 1.0;

    setStandardValues();

    addParam("b0", b0);
    mb0 = b0;

    mLocTeds.push_back(enum_nat_tetrad_static);
    mLocTeds.push_back(enum_nat_tetrad_lnrf);
    pm = 1.0;
}

MetricTeoWHl::~MetricTeoWHl() {

}


// *********************************** public methods ******************************
/*! Calculate the contravariant metric components at position 'pos'.
 *
 *  \param pos : pointer to position.
 */
bool MetricTeoWHl::calculateMetric(const double* pos) {
    double l     = pos[1];
    double theta = pos[2];
    double c = 1.0;

    double t1 = calc_N(l);
    double t2 = t1 * t1;
    double t3 = c * c;
    double t5 = calc_r(l);
    double t6 = t5 * t5;
    double t7 = calc_K(l);
    double t8 = t7 * t7;
    double t9 = t6 * t8;
    double t10 = sin(theta);
    double t11 = t10 * t10;
    double t12 = calc_omega(l);
    double t13 = t12 * t12;
    double t20 = t9 * t11 * t12 * c;

    g_compts[0][0] = -t2 * t3 + t9 * t11 * t13 * t3;
    g_compts[0][1] = 0.0;
    g_compts[0][2] = 0.0;
    g_compts[0][3] = -t20;
    g_compts[1][0] = 0.0;
    g_compts[1][1] = 1.0;
    g_compts[1][2] = 0.0;
    g_compts[1][3] = 0.0;
    g_compts[2][0] = 0.0;
    g_compts[2][1] = 0.0;
    g_compts[2][2] = t9;
    g_compts[2][3] = 0.0;
    g_compts[3][0] = -t20;
    g_compts[3][1] = 0.0;
    g_compts[3][2] = 0.0;
    g_compts[3][3] = t9 * t11;

    return true;
}

/*! Calculate the Christoffel symbols of the second kind at position 'pos'.
 *
 *  \param pos : pointer to position.
 */
bool MetricTeoWHl::calculateChristoffels(const double* pos) {
    double l     = pos[1];
    double theta = pos[2];
    double c = 1.0;

    double t1 = calc_N(l);
    double t2 = c * c;
    double t4 = calc_dN(l);  // diff(N(l),l);
    double t6 = calc_r(l);
    double t7 = calc_K(l);
    double t8 = t7 * t7;
    double t9 = t6 * t8;
    double t10 = sin(theta);
    double t11 = t10 * t10;
    double t12 = t9 * t11;
    double t13 = calc_omega(l);
    double t14 = t13 * t13;
    double t15 = t14 * t2;
    double t16 = calc_dr(l);  // diff(r(l),l);
    double t19 = t6 * t6;
    double t20 = t19 * t7;
    double t21 = t20 * t11;
    double t22 = calc_dK(l);  // diff(K(l),l);
    double t25 = t19 * t8;
    double t26 = t25 * t11;
    double t28 = calc_domega(l);  // diff(omega(l),l);
    double t33 = cos(theta);
    double t42 = t1 * t1;
    double t43 = 1 / t42;
    double t45 = (2.0 * t1 * t4 - t25 * t11 * t13 * t28) * t43 / 2.0;
    double t46 = 1 / t6;
    double t48 = 1 / t7;
    double t54 = t19 * t6;
    double t58 = t8 * t7 * t11 * t28;
    double t60 = t42 * t7;
    double t64 = t42 * t6;
    double t74 = c * t46 * t48 * (2.0 * t13 * t6 * t7 * t1 * t4 - t14 * t54 * t58 - 2.0 * t60 * t13 * t16 - 2.0 * t64 * t13 * t22 - t64 * t7 * t28) * t43 / 2.0;
    double t75 = t13 * c;
    double t77 = t33 / t10;
    double t78 = t75 * t77;
    double t87 = t12 * t75 * t16 + t21 * t75 * t22 + t25 * t11 * t28 * c / 2.0;
    double t90 = t10 * t13 * c * t33;
    double t91 = t46 * t48;
    double t95 = t91 * (t7 * t16 + t6 * t22);
    double t100 = t26 * t28 / c * t43 / 2.0;
    double t110 = t91 * (t13 * t54 * t58 + 2.0 * t60 * t16 + 2.0 * t64 * t22) * t43 / 2.0;

    christoffel[0][0][0] = 0.0;
    christoffel[0][0][1] = t1 * t2 * t4 - t12 * t15 * t16 - t21 * t15 * t22 - t26 * t13 * t2 * t28;
    christoffel[0][0][2] = -t10 * t14 * t2 * t33;
    christoffel[0][0][3] = 0.0;
    christoffel[0][1][0] = t45;
    christoffel[0][1][1] = 0.0;
    christoffel[0][1][2] = 0.0;
    christoffel[0][1][3] = t74;
    christoffel[0][2][0] = 0.0;
    christoffel[0][2][1] = 0.0;
    christoffel[0][2][2] = 0.0;
    christoffel[0][2][3] = -t78;
    christoffel[0][3][0] = 0.0;
    christoffel[0][3][1] = t87;
    christoffel[0][3][2] = t90;
    christoffel[0][3][3] = 0.0;
    christoffel[1][0][0] = t45;
    christoffel[1][0][1] = 0.0;
    christoffel[1][0][2] = 0.0;
    christoffel[1][0][3] = t74;
    christoffel[1][1][0] = 0.0;
    christoffel[1][1][1] = 0.0;
    christoffel[1][1][2] = 0.0;
    christoffel[1][1][3] = 0.0;
    christoffel[1][2][0] = 0.0;
    christoffel[1][2][1] = 0.0;
    christoffel[1][2][2] = t95;
    christoffel[1][2][3] = 0.0;
    christoffel[1][3][0] = t100;
    christoffel[1][3][1] = 0.0;
    christoffel[1][3][2] = 0.0;
    christoffel[1][3][3] = t110;
    christoffel[2][0][0] = 0.0;
    christoffel[2][0][1] = 0.0;
    christoffel[2][0][2] = 0.0;
    christoffel[2][0][3] = -t78;
    christoffel[2][1][0] = 0.0;
    christoffel[2][1][1] = 0.0;
    christoffel[2][1][2] = t95;
    christoffel[2][1][3] = 0.0;
    christoffel[2][2][0] = 0.0;
    christoffel[2][2][1] = -t9 * t16 - t20 * t22;
    christoffel[2][2][2] = 0.0;
    christoffel[2][2][3] = 0.0;
    christoffel[2][3][0] = 0.0;
    christoffel[2][3][1] = 0.0;
    christoffel[2][3][2] = 0.0;
    christoffel[2][3][3] = t77;
    christoffel[3][0][0] = 0.0;
    christoffel[3][0][1] = t87;
    christoffel[3][0][2] = t90;
    christoffel[3][0][3] = 0.0;
    christoffel[3][1][0] = t100;
    christoffel[3][1][1] = 0.0;
    christoffel[3][1][2] = 0.0;
    christoffel[3][1][3] = t110;
    christoffel[3][2][0] = 0.0;
    christoffel[3][2][1] = 0.0;
    christoffel[3][2][2] = 0.0;
    christoffel[3][2][3] = t77;
    christoffel[3][3][0] = 0.0;
    christoffel[3][3][1] = -t9 * t11 * t16 - t20 * t11 * t22;
    christoffel[3][3][2] = -t10 * t33;
    christoffel[3][3][3] = 0.0;

    return true;
}

/*! Calculate Jacobi matrix.
 *
 *  \param pos : pointer to position.
 */
bool MetricTeoWHl::calculateChrisD(const double* pos) {
    double l     = pos[1];
    double theta = pos[2];
    double c = 1.0;

    double t1 = calc_dN(l);  //diff(N(l),l);
    double t2 = t1 * t1;
    double t3 = c * c;
    double t5 = calc_N(l);
    double t7 = calc_d2N(l);  // diff(diff(N(l),l),l);
    double t9 = calc_K(l);
    double t10 = t9 * t9;
    double t11 = calc_dr(l);  // diff(r(l),l);
    double t12 = t11 * t11;
    double t13 = t10 * t12;
    double t14 = sin(theta);
    double t15 = t14 * t14;
    double t16 = calc_omega(l);
    double t17 = t16 * t16;
    double t18 = t15 * t17;
    double t19 = t18 * t3;
    double t21 = calc_r(l);
    double t22 = t21 * t9;
    double t23 = t22 * t15;
    double t24 = t17 * t3;
    double t25 = calc_dK(l);  // diff(K(l),l);
    double t26 = t11 * t25;
    double t30 = t21 * t10;
    double t31 = t30 * t15;
    double t32 = t16 * t3;
    double t33 = calc_domega(l);  // diff(omega(l),l);
    double t34 = t11 * t33;
    double t38 = calc_d2r(l);     // diff(diff(r(l),l),l);
    double t41 = t21 * t21;
    double t42 = t25 * t25;
    double t43 = t41 * t42;
    double t45 = t41 * t9;
    double t46 = t45 * t15;
    double t47 = t25 * t33;
    double t51 = calc_d2K(l);     // diff(diff(K(l),l),l);
    double t54 = t41 * t10;
    double t55 = t33 * t33;
    double t60 = calc_d2omega(l); // diff(diff(omega(l),l),l);
    double t63 = t2 * t3 + t5 * t3 * t7 - t13 * t19 - 4.0 * t23 * t24 * t26 - 4.0 * t31 * t32 * t34 - t31 * t24 * t38 - t43 * t19 - 4.0 * t46 * t32 * t47 - t46 * t24 * t51 - t54 * t15 * t55 * t3 - t54 * t15 * t32 * t60;
    double t64 = t30 * t14;
    double t65 = cos(theta);
    double t66 = t11 * t65;
    double t69 = t45 * t14;
    double t70 = t25 * t65;
    double t73 = t54 * t14;
    double t83 = t65 * t65;
    double t89 = t5 * t5;
    double t94 = t15 * t16;
    double t98 = t5 * t41;
    double t107 = t94 * t60;
    double t109 = t1 * t41;
    double t111 = t94 * t33;
    double t115 = t89 * t5;
    double t116 = 1 / t115;
    double t118 = (2.0 * t5 * t2 - 2.0 * t89 * t7 + 2.0 * t5 * t21 * t10 * t94 * t34 + 2.0 * t98 * t9 * t94 * t47 + t98 * t10 * t15 * t55 + t98 * t10 * t107 - 2.0 * t109 * t10 * t111) * t116 / 2.0;
    double t121 = t65 / t89;
    double t123 = t73 * t16 * t33 * t121;
    double t124 = t115 * t16;
    double t131 = t10 * t10;
    double t133 = t11 * t131 * t5;
    double t134 = t41 * t21;
    double t136 = t15 * t33;
    double t145 = t41 * t41;
    double t147 = t25 * t145 * t5;
    double t148 = t10 * t9;
    double t166 = t145 * t131;
    double t167 = t166 * t5;
    double t180 = t1 * t145 * t131;
    double t188 = -2.0 * t13 * t124 - 2.0 * t43 * t124 + t54 * t115 * t60 + 2.0 * t133 * t17 * t134 * t136 + 2.0 * t11 * t10 * t115 * t21 * t33 + 2.0 * t147 * t17 * t148 * t136 + 2.0 * t25 * t41 * t115 * t9 * t33 + 2.0 * t54 * t5 * t16 * t2 - 2.0 * t54 * t89 * t16 * t7 + 2.0 * t167 * t94 * t55 + t167 * t18 * t60 + 2.0 * t30 * t124 * t38 + 2.0 * t45 * t124 * t51 - 2.0 * t180 * t18 * t33 - 2.0 * t109 * t10 * t89 * t33;
    double t190 = 1 / t41;
    double t191 = 1 / t10;
    double t195 = c * t188 * t190 * t191 * t116 / 2.0;
    double t199 = t14 * t33;
    double t201 = c * t41 * t10 * t17 * t199 * t121;
    double t202 = t33 * c;
    double t205 = t202 * t65 / t14;
    double t206 = t16 * c;
    double t209 = (t15 + t83) / t15;
    double t210 = t206 * t209;
    double t211 = t94 * c;
    double t231 = t13 * t211 + 4.0 * t23 * t206 * t26 + 2.0 * t31 * t202 * t11 + t31 * t206 * t38 + t43 * t211 + 2.0 * t46 * t202 * t25 + t46 * t206 * t51 + t54 * t15 * t60 * c / 2.0;
    double t240 = 2.0 * t64 * t206 * t66 + 2.0 * t69 * t206 * t70 + t73 * t202 * t65;
    double t242 = t199 * c * t65;
    double t245 = t83 * t16 * c - t211;
    double t246 = t30 * t38;
    double t247 = t45 * t51;
    double t250 = (t13 + t43 - t246 - t247) * t190 * t191;
    double t265 = 1 / c;
    double t269 = t23 * (-2.0 * t9 * t33 * t11 * t5 - 2.0 * t21 * t33 * t25 * t5 - t22 * t60 * t5 + 2.0 * t22 * t33 * t1) * t265 * t116 / 2.0;
    double t272 = t73 * t33 * t265 * t121;
    double t301 = (-2.0 * t133 * t16 * t134 * t136 + 2.0 * t13 * t115 - 2.0 * t147 * t16 * t148 * t136 + 2.0 * t43 * t115 - t166 * t5 * t55 * t15 - t167 * t107 - 2.0 * t30 * t115 * t38 - 2.0 * t45 * t115 * t51 + 2.0 * t180 * t111) * t190 * t191 * t116 / 2.0;

    chrisD[0][0][0][0] = 0.0;
    chrisD[0][0][0][1] = 0.0;
    chrisD[0][0][0][2] = 0.0;
    chrisD[0][0][0][3] = 0.0;
    chrisD[0][0][1][0] = 0.0;
    chrisD[0][0][1][1] = t63;
    chrisD[0][0][1][2] = -2.0 * t64 * t24 * t66 - 2.0 * t69 * t24 * t70 - 2.0 * t73 * t32 * t33 * t65;
    chrisD[0][0][1][3] = 0.0;
    chrisD[0][0][2][0] = 0.0;
    chrisD[0][0][2][1] = -2.0 * t14 * t16 * t3 * t65 * t33;
    chrisD[0][0][2][2] = -t83 * t17 * t3 + t19;
    chrisD[0][0][2][3] = 0.0;
    chrisD[0][0][3][0] = 0.0;
    chrisD[0][0][3][1] = 0.0;
    chrisD[0][0][3][2] = 0.0;
    chrisD[0][0][3][3] = 0.0;
    chrisD[0][1][0][0] = 0.0;
    chrisD[0][1][0][1] = -t118;
    chrisD[0][1][0][2] = -t123;
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
    chrisD[0][1][3][1] = -t195;
    chrisD[0][1][3][2] = -t201;
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
    chrisD[0][2][3][1] = -t205;
    chrisD[0][2][3][2] = t210;
    chrisD[0][2][3][3] = 0.0;
    chrisD[0][3][0][0] = 0.0;
    chrisD[0][3][0][1] = 0.0;
    chrisD[0][3][0][2] = 0.0;
    chrisD[0][3][0][3] = 0.0;
    chrisD[0][3][1][0] = 0.0;
    chrisD[0][3][1][1] = t231;
    chrisD[0][3][1][2] = t240;
    chrisD[0][3][1][3] = 0.0;
    chrisD[0][3][2][0] = 0.0;
    chrisD[0][3][2][1] = t242;
    chrisD[0][3][2][2] = t245;
    chrisD[0][3][2][3] = 0.0;
    chrisD[0][3][3][0] = 0.0;
    chrisD[0][3][3][1] = 0.0;
    chrisD[0][3][3][2] = 0.0;
    chrisD[0][3][3][3] = 0.0;
    chrisD[1][0][0][0] = 0.0;
    chrisD[1][0][0][1] = -t118;
    chrisD[1][0][0][2] = -t123;
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
    chrisD[1][0][3][1] = -t195;
    chrisD[1][0][3][2] = -t201;
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
    chrisD[1][2][2][1] = -t250;
    chrisD[1][2][2][2] = 0.0;
    chrisD[1][2][2][3] = 0.0;
    chrisD[1][2][3][0] = 0.0;
    chrisD[1][2][3][1] = 0.0;
    chrisD[1][2][3][2] = 0.0;
    chrisD[1][2][3][3] = 0.0;
    chrisD[1][3][0][0] = 0.0;
    chrisD[1][3][0][1] = -t269;
    chrisD[1][3][0][2] = t272;
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
    chrisD[1][3][3][1] = -t301;
    chrisD[1][3][3][2] = t123;
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
    chrisD[2][0][3][1] = -t205;
    chrisD[2][0][3][2] = t210;
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
    chrisD[2][1][2][1] = -t250;
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
    chrisD[2][2][1][1] = -t13 - 4.0 * t22 * t26 - t246 - t43 - t247;
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
    chrisD[2][3][3][2] = -t209;
    chrisD[2][3][3][3] = 0.0;
    chrisD[3][0][0][0] = 0.0;
    chrisD[3][0][0][1] = 0.0;
    chrisD[3][0][0][2] = 0.0;
    chrisD[3][0][0][3] = 0.0;
    chrisD[3][0][1][0] = 0.0;
    chrisD[3][0][1][1] = t231;
    chrisD[3][0][1][2] = t240;
    chrisD[3][0][1][3] = 0.0;
    chrisD[3][0][2][0] = 0.0;
    chrisD[3][0][2][1] = t242;
    chrisD[3][0][2][2] = t245;
    chrisD[3][0][2][3] = 0.0;
    chrisD[3][0][3][0] = 0.0;
    chrisD[3][0][3][1] = 0.0;
    chrisD[3][0][3][2] = 0.0;
    chrisD[3][0][3][3] = 0.0;
    chrisD[3][1][0][0] = 0.0;
    chrisD[3][1][0][1] = -t269;
    chrisD[3][1][0][2] = t272;
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
    chrisD[3][1][3][1] = -t301;
    chrisD[3][1][3][2] = t123;
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
    chrisD[3][2][3][2] = -t209;
    chrisD[3][2][3][3] = 0.0;
    chrisD[3][3][0][0] = 0.0;
    chrisD[3][3][0][1] = 0.0;
    chrisD[3][3][0][2] = 0.0;
    chrisD[3][3][0][3] = 0.0;
    chrisD[3][3][1][0] = 0.0;
    chrisD[3][3][1][1] = -t13 * t15 - 4.0 * t22 * t15 * t11 * t25 - t30 * t15 * t38 - t43 * t15 - t45 * t15 * t51;
    chrisD[3][3][1][2] = -2.0 * t30 * t14 * t11 * t65 - 2.0 * t45 * t14 * t25 * t65;
    chrisD[3][3][1][3] = 0.0;
    chrisD[3][3][2][0] = 0.0;
    chrisD[3][3][2][1] = 0.0;
    chrisD[3][3][2][2] = -t83 + t15;
    chrisD[3][3][2][3] = 0.0;
    chrisD[3][3][3][0] = 0.0;
    chrisD[3][3][3][1] = 0.0;
    chrisD[3][3][3][2] = 0.0;
    chrisD[3][3][3][3] = 0.0;

    return true;
}


/*! Transform local 4-direction to coordinate 4-direction.
 *
 *  \param  pos  :  pointer to position array.
 *  \param  ldir :  pointer to local direction array.
 *  \param  dir  :  pointer to calculated coordinate direction array.
 *  \param  type :  type of tetrad.
 */
void MetricTeoWHl::localToCoord(const double* pos, const double* ldir, double* dir,
                                enum_nat_tetrad_type  type) {
    double l     = pos[1];
    double theta = pos[2];

    calc_N(l);
    calc_r(l);
    calc_K(l);
    calc_omega(l);

    if (type == enum_nat_tetrad_lnrf) {
        dir[0] = ldir[0] / curr_N;
        dir[1] = ldir[1];
        dir[2] = ldir[2] / (curr_r * curr_K);
        dir[3] = ldir[0] * curr_omega / curr_N + ldir[3] / (curr_r * curr_K * sin(theta));
    } else {
        double w = sqrt(curr_N * curr_N - curr_r * curr_r * curr_K * curr_K * sin(theta) * sin(theta) * curr_omega * curr_omega);
        dir[0] = ldir[0] / w - ldir[3] * pm * curr_r * curr_K * sin(theta) * curr_omega / curr_N / w;
        dir[1] = ldir[1];
        dir[2] = ldir[2] / (curr_r * curr_K);
        dir[3] = ldir[3] * pm * w / (curr_N * curr_r * curr_K * sin(theta));
    }
}

/*! Transform coordinate 4-direction to local 4-direction.
 *
 *  \param  pos  :  pointer to position array.
 *  \param  cdir :  pointer to coordinate direction.
 *  \param  ldir :  pointer to calculated local direction array.
 *  \param  type :  type of tetrad.
 */
void MetricTeoWHl::coordToLocal(const double* , const double* , double* ,
                                enum_nat_tetrad_type) {
    // TODO
}


/*! Test break condition.
 *
 *  \param pos  :  position.
 *  \return true  : radial position r < 0.0 or  r^2<=(1.0+eps)*rs^2.
 *  \return false : position is valid.
 */
bool MetricTeoWHl::breakCondition(const double*) {
    bool br = false;
    return br;
}

/** Transform proper coordinates to pseudo Cartesian coordinates
 * \param p
 * \param cp
 * \return chart ID
 */
int MetricTeoWHl::transToPseudoCart(vec4 p, vec4 &cp) {
    TransCoordinates::toCartesianCoord(mCoordType, p, cp);
    if (p[1] > 0) {
        return 0;
    }
    return 1;
}

/*!
 *  Set 'mass' parameter and adjust Schwarzschild radius  rs=2GM/c^2.
 */
bool MetricTeoWHl::setParam(std::string pName, double val) {
    bool ok = Metric::setParam(pName, val);

    return ok;
}

/*! Generate report.
 */
bool MetricTeoWHl::report(const vec4 , const vec4 , std::string &text) {
    std::stringstream ss;
    ss << "Report for Teo wormhole metric\n\tcoordinate : (t,l,theta,phi)\n";
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
void MetricTeoWHl::setStandardValues() {
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

/*!
 *  \param l : proper radial coordinate.
 */
double MetricTeoWHl::calc_N(double) {
    curr_N = 1.0;
    return curr_N;
}

/*!
 *  \param l : proper radial coordinate.
 */
double MetricTeoWHl::calc_dN(double) {
    curr_dNdl = 0.0;
    return curr_dNdl;
}

/*!
 *  \param l : proper radial coordinate.
 */
double MetricTeoWHl::calc_d2N(double) {
    return 0.0;
}

/*!
 *  \param l : proper radial coordinate.
 */
double MetricTeoWHl::calc_K(double) {
    curr_K = 1.0;
    return curr_K;
}

/*!
 *  \param l : proper radial coordinate.
 */
double MetricTeoWHl::calc_dK(double) {
    curr_dKdl = 0.0;
    return curr_dKdl;
}

/*!
 *  \param l : proper radial coordinate.
 */
double MetricTeoWHl::calc_d2K(double) {
    return 0.0;
}

/*!
 *  \param l : proper radial coordinate.
 */
double MetricTeoWHl::calc_r(double l) {
    curr_r = sqrt(l * l + mb0 * mb0);
    return curr_r;
}

/*!
 *  \param l : proper radial coordinate.
 */
double MetricTeoWHl::calc_dr(double l) {
    curr_drdl = l / sqrt(l * l + mb0 * mb0);
    return curr_drdl;
}

/*!
 *  \param l : proper radial coordinate.
 */
double MetricTeoWHl::calc_d2r(double l) {
    return mb0 * mb0 / pow(l * l + mb0 * mb0, 1.5);
}

/*!
 *  \param l : proper radial coordinate.
 */
double MetricTeoWHl::calc_omega(double l) {
    double r = calc_r(l);
    curr_omega = 0.5 * mb0 * mb0 * pow(r, -3.0);
    return curr_omega;
}

/*!
 *  \param l : proper radial coordinate.
 */
double MetricTeoWHl::calc_domega(double l) {
    curr_domegadl = -1.5 * mb0 * mb0 * l * pow(l * l + mb0 * mb0, -2.5);
    return curr_domegadl;
}

/*!
 *  \param l : proper radial coordinate.
 */
double MetricTeoWHl::calc_d2omega(double) {
    return 0.0;
}

} // end namespace m4d
