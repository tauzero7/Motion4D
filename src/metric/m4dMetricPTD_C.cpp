// -------------------------------------------------------------------------------
/*
   m4dMetricPTD_C.cpp

  Copyright (c) 2010-2014  Thomas Mueller, Frank Grave, Felix Beslmeisl


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

#include "m4dMetricPTD_C.h"

namespace m4d {

#define eps 1.0e-6


/*! Standard constructor for the metric.
 *
 * \param  a : parameter "a" for the metric
 * \param  b : parameter "b" for the metric
 */
MetricPTD_C::MetricPTD_C(double a, double b) {
    mMetricName  = "Petrov_Type_D_C_ES";
    setCoordType(enum_coordinate_custom);

    mPhysicalUnits = enum_physical_constants_geom;
    mSpeedOfLight = 1.0;
    mGravConstant = 1.0;


    Par_a = a;
    Par_b = b;

    addParam("a", Par_a);
    addParam("b", Par_b);

    setStandardValues();


}

/*! Standard destructor for the metric.
 *
 */
MetricPTD_C::~MetricPTD_C() {

}


// *********************************** public methods ******************************
/*! Calculate the contravariant metric components at position 'pos'.
 *
 *  \param pos : pointer to position.
 */
bool MetricPTD_C::calculateMetric(const double* pos) {
    double x = pos[2];
    double y = pos[3];
    double a = Par_a;
    double b = Par_b;

    double t1 = y * y;
    double t4 = t1 * y + y * a - b;
    double t6 = pow(x + y, 2.0);
    double t7 = 1 / t6;
    double t9 = x * x;
    double t12 = t9 * x + x * a + b;
    g_compts[0][0] = t4 * t7;
    g_compts[0][1] = 0.0;
    g_compts[0][2] = 0.0;
    g_compts[0][3] = 0.0;
    g_compts[1][0] = 0.0;
    g_compts[1][1] = t7 * t12;
    g_compts[1][2] = 0.0;
    g_compts[1][3] = 0.0;
    g_compts[2][0] = 0.0;
    g_compts[2][1] = 0.0;
    g_compts[2][2] = t7 / t12;
    g_compts[2][3] = 0.0;
    g_compts[3][0] = 0.0;
    g_compts[3][1] = 0.0;
    g_compts[3][2] = 0.0;
    g_compts[3][3] = -t7 / t4;

    return true;
}

/*! Calculate the Christoffel symbols of the second kind at position 'pos'.
 *
 *  \param pos : pointer to position.
 */
bool MetricPTD_C::calculateChristoffels(const double* pos) {
    double x = pos[2];
    double y = pos[3];
    double a = Par_a;
    double b = Par_b;

    double t2 = 1 / (x + y);
    double t3 = x * x;
    double t4 = t3 * x;
    double t5 = x * a;
    double t6 = t4 + t5 + b;
    double t7 = t2 * t6;
    double t8 = y * y;
    double t9 = t8 * y;
    double t10 = y * a;
    double t11 = t9 + t10 - b;
    double t12 = t7 * t11;
    double t13 = t2 * t11;
    double t15 = 3.0 * t8 * x;
    double t16 = 2.0 * b;
    double t17 = t15 + t9 + t5 - t10 + t16;
    double t20 = 1 / t11;
    double t21 = t2 * t20;
    double t23 = t21 * t17 / 2.0;
    double t25 = 3.0 * t3 * y;
    double t26 = -t4 + t5 + t16 - t25 - t10;
    double t29 = 1 / t6;
    double t30 = t2 * t29;
    double t32 = t30 * t26 / 2.0;
    christoffel[0][0][0] = 0.0;
    christoffel[0][0][1] = 0.0;
    christoffel[0][0][2] = t12;
    christoffel[0][0][3] = t13 * t17 / 2.0;
    christoffel[0][1][0] = 0.0;
    christoffel[0][1][1] = 0.0;
    christoffel[0][1][2] = 0.0;
    christoffel[0][1][3] = 0.0;
    christoffel[0][2][0] = -t2;
    christoffel[0][2][1] = 0.0;
    christoffel[0][2][2] = 0.0;
    christoffel[0][2][3] = 0.0;
    christoffel[0][3][0] = t23;
    christoffel[0][3][1] = 0.0;
    christoffel[0][3][2] = 0.0;
    christoffel[0][3][3] = 0.0;
    christoffel[1][0][0] = 0.0;
    christoffel[1][0][1] = 0.0;
    christoffel[1][0][2] = 0.0;
    christoffel[1][0][3] = 0.0;
    christoffel[1][1][0] = 0.0;
    christoffel[1][1][1] = 0.0;
    christoffel[1][1][2] = t7 * t26 / 2.0;
    christoffel[1][1][3] = -t12;
    christoffel[1][2][0] = 0.0;
    christoffel[1][2][1] = -t32;
    christoffel[1][2][2] = 0.0;
    christoffel[1][2][3] = 0.0;
    christoffel[1][3][0] = 0.0;
    christoffel[1][3][1] = -t2;
    christoffel[1][3][2] = 0.0;
    christoffel[1][3][3] = 0.0;
    christoffel[2][0][0] = -t2;
    christoffel[2][0][1] = 0.0;
    christoffel[2][0][2] = 0.0;
    christoffel[2][0][3] = 0.0;
    christoffel[2][1][0] = 0.0;
    christoffel[2][1][1] = -t32;
    christoffel[2][1][2] = 0.0;
    christoffel[2][1][3] = 0.0;
    christoffel[2][2][0] = 0.0;
    christoffel[2][2][1] = 0.0;
    christoffel[2][2][2] = -t30 * (5.0 * t4 + 3.0 * t5 + t16 + t25 + t10) / 2.0;
    christoffel[2][2][3] = -t13 * t29;
    christoffel[2][3][0] = 0.0;
    christoffel[2][3][1] = 0.0;
    christoffel[2][3][2] = -t2;
    christoffel[2][3][3] = -t2;
    christoffel[3][0][0] = t23;
    christoffel[3][0][1] = 0.0;
    christoffel[3][0][2] = 0.0;
    christoffel[3][0][3] = 0.0;
    christoffel[3][1][0] = 0.0;
    christoffel[3][1][1] = -t2;
    christoffel[3][1][2] = 0.0;
    christoffel[3][1][3] = 0.0;
    christoffel[3][2][0] = 0.0;
    christoffel[3][2][1] = 0.0;
    christoffel[3][2][2] = -t2;
    christoffel[3][2][3] = -t2;
    christoffel[3][3][0] = 0.0;
    christoffel[3][3][1] = 0.0;
    christoffel[3][3][2] = -t7 * t20;
    christoffel[3][3][3] = -t21 * (5.0 * t9 + 3.0 * t10 - t16 + t15 + t5) / 2.0;

    return true;
}

/*! Calculate Jacobi matrix.
 *
 *  \param pos : pointer to position.
 */
bool
MetricPTD_C::calculateChrisD(const double* pos) {
    double x = pos[2];
    double y = pos[3];
    double a = Par_a;
    double b = Par_b;

    double t1 = x * x;
    double t2 = t1 * x;
    double t3 = y * y;
    double t4 = t3 * y;
    double t6 = 4.0 * t2 * t4;
    double t8 = t2 * y * a;
    double t9 = 4.0 * t8;
    double t10 = t2 * b;
    double t11 = 4.0 * t10;
    double t12 = x * a;
    double t13 = t12 * t4;
    double t14 = 8.0 * t13;
    double t15 = a * a;
    double t18 = 6.0 * x * t15 * y;
    double t19 = t12 * b;
    double t20 = 8.0 * t19;
    double t21 = b * y;
    double t22 = t21 * a;
    double t23 = 8.0 * t22;
    double t24 = b * b;
    double t25 = 8.0 * t24;
    double t27 = t3 * t1 * a;
    double t28 = 6.0 * t27;
    double t29 = t3 * x;
    double t30 = t29 * b;
    double t31 = 12.0 * t30;
    double t32 = t3 * t3;
    double t33 = t32 * t1;
    double t34 = 9.0 * t33;
    double t36 = t32 * y * x;
    double t37 = 6.0 * t36;
    double t38 = t32 * a;
    double t39 = 2.0 * t38;
    double t40 = t1 * t15;
    double t41 = t3 * t15;
    double t42 = t32 * t3;
    double t43 = -t6 - t9 + t11 - t14 - t18 + t20 - t23 + t25 + t28 + t31 + t34 + t37 - t39 + t40 + t41 + t42;
    double t44 = x + y;
    double t45 = t44 * t44;
    double t46 = 1 / t45;
    double t49 = y * a;
    double t50 = t4 + t49 - b;
    double t51 = 2.0 * t12;
    double t52 = 4.0 * b;
    double t54 = 3.0 * t1 * y;
    double t55 = 2.0 * t49;
    double t56 = 3.0 * t29;
    double t57 = -t2 + t51 + t52 - t54 - t55 + t56 + t4;
    double t60 = t50 * t57 * t46 / 2.0;
    double t61 = 3.0 * t2;
    double t62 = 2.0 * b;
    double t63 = t61 - t62 + t54 + t55 - t56 - t4;
    double t67 = t2 + t12 + b;
    double t68 = t46 * t67;
    double t69 = t56 + t4 + t12 - t49 + t62;
    double t72 = t46 * t50;
    double t75 = 4.0 * t19;
    double t76 = 4.0 * t22;
    double t77 = 24.0 * t30;
    double t78 = t21 * t1;
    double t79 = 12.0 * t78;
    double t80 = b * t4;
    double t82 = 4.0 * t24;
    double t83 = 3.0 * t33;
    double t84 = 6.0 * t38;
    double t86 = t9 - t14 + t18 - t75 + t76 - t28 + t77 + t79 + t6 - t11 + 16.0 * t80 - t82 - t83 - t37 - t84 + t40 +
                 t41 - 3.0 * t42;
    double t87 = t86 * t46;
    double t91 = t67 * t57 * t46 / 2.0;
    double t92 = t63 * t46;
    double t93 = 1 / t67;
    double t95 = t92 * t93 / 2.0;
    double t96 = 1 / t50;
    double t97 = t46 * t96;
    double t99 = t97 * t69 / 2.0;
    double t100 = 5.0 * t2;
    double t101 = 4.0 * t12;
    double t105 = t50 * (t100 + t101 + t52 + t54 + t56 + t4) * t46 / 2.0;
    double t106 = 1 / t44;
    double t109 = 3.0 * t3 + a;
    double t111 = t106 * t50 * t109 / 2.0;
    double t112 = t50 * t50;
    double t113 = 1 / t112;
    double t115 = t87 * t113 / 4.0;
    double t118 = t106 * t67 * t109 / 2.0;
    double t120 = 12.0 * t80;
    double t121 = 12.0 * t22;
    double t122 = 18.0 * t36;
    double t123 = 3.0 * t41;
    double t125 = t6 + t9 - t11 + 16.0 * t13 + t18 - t75 + t120 + t121 - t25 + t28 + t34 + t122 - t39 + t40 - t123
                  + 5.0 * t42;
    double t127 = t125 * t46 / 4.0;
    double t128 = 8.0 * t8;
    double t129 = t1 * t1;
    double t130 = t129 * a;
    double t131 = 2.0 * t130;
    double t133 = t129 * x * y;
    double t134 = 6.0 * t133;
    double t135 = t129 * t3;
    double t136 = 9.0 * t135;
    double t137 = t129 * t1;
    double t138 = 4.0 * t13;
    double t139 = 4.0 * t80;
    double t140 = -t128 - t18 + t20 - t23 + t28 - t79 + t25 + t40 + t41 - t131 + t134 + t136 + t137 - t6 - t138 -
                  t139;
    double t143 = 24.0 * t78;
    double t145 = 6.0 * t130;
    double t146 = 3.0 * t135;
    double t148 = -t128 + t138 + t18 - t75 + t76 - t28 - t31 - t143 + t6 - 16.0 * t10 + t139 - t82 + t40 + t41 -
                  t145 - t134 - t146 - 3.0 * t137;
    double t149 = t148 * t46;
    double t151 = -t2 + t12 + t62 - t54 - t49;
    double t156 = 3.0 * t4;
    double t157 = t156 + t62 + t56 + t51 - t2 - t54;
    double t161 = t67 * t67;
    double t162 = 1 / t161;
    double t164 = t149 * t162 / 4.0;
    double t165 = t46 * t93;
    double t167 = t165 * t151 / 2.0;
    double t169 = 12.0 * t19;
    double t170 = 12.0 * t10;
    double t171 = 3.0 * t40;
    double t172 = 18.0 * t133;
    double t174 = -16.0 * t8 - t18 + t169 - t28 + t170 + t25 + t171 - t41 + t131 - t172 - t136 - 5.0 * t137 - t6
                  - t138 - t139 - t76;
    double t176 = t174 * t46 / 4.0;
    double t178 = 3.0 * t1 + a;
    double t181 = t106 * t178 * t50 / 2.0;
    double t182 = t157 * t46;
    double t184 = t182 * t96 / 2.0;
    double t187 = t67 * t178 * t106 / 2.0;
    double t188 = 4.0 * t49;
    double t189 = 5.0 * t4;
    double t193 = t67 * (t2 - t52 + t54 + t188 + t189 + t56) * t46 / 2.0;
    double t194 = t57 * t46;
    double t197 = 18.0 * t27;
    double t201 = t138 + t18 - t75 + t121 - t197 + t143 + t6 + t170 + t139 - t25 + t40 - t123 + t145 - 30.0 *
                  t133 - 27.0 * t135 - 7.0 * t137;
    double t210 = 32.0 * t8 + t138 + t18 + t23 + t28 + t31 + 36.0 * t78 + t6 + 20.0 * t10 + t139 - t82 + t171 -
                  t41 + 22.0 * t130 + t172 + t146 + 15.0 * t137;
    double t217 = t165 * (t100 + 3.0 * t12 + t62 + t54 + t49) / 2.0;
    double t226 = t106 * t96 * t109 / 2.0;
    double t229 = t106 * t93 * t178 / 2.0;
    double t233 = (t156 + t55 + t62 + t100 + t101 + t54 + t56) * t46 * t96 / 2.0;
    double t237 = (t61 + t51 - t62 + t54 + t188 + t189 + t56) * t46 * t93 / 2.0;
    double t241 = t97 * (t189 + 3.0 * t49 - t62 + t56 + t12) / 2.0;
    double t245 = -t9 - t18 + t169 - t76 + t197 + t77 - t6 + t11 + t120 + t25 + 27.0 * t33 + 30.0 * t36 - t84 +
                  t171 - t41 + 7.0 * t42;
    double t262 = -t9 - 32.0 * t13 - t18 + t20 - t28 + 36.0 * t30 + t79 - t6 + t11 + 20.0 * t80 + t82 - t83 - t122
                  - 22.0 * t38 + t40 - t123 - 15.0 * t42;
    chrisD[0][0][0][0] = -t43 * t46 / 4.0;
    chrisD[0][0][0][1] = 0.0;
    chrisD[0][0][0][2] = 0.0;
    chrisD[0][0][0][3] = 0.0;
    chrisD[0][0][1][0] = 0.0;
    chrisD[0][0][1][1] = -t60;
    chrisD[0][0][1][2] = 0.0;
    chrisD[0][0][1][3] = 0.0;
    chrisD[0][0][2][0] = 0.0;
    chrisD[0][0][2][1] = 0.0;
    chrisD[0][0][2][2] = t50 * t63 * t46 / 2.0;
    chrisD[0][0][2][3] = -t68 * t69 / 2.0;
    chrisD[0][0][3][0] = 0.0;
    chrisD[0][0][3][1] = 0.0;
    chrisD[0][0][3][2] = t72 * t69 / 2.0;
    chrisD[0][0][3][3] = -t87 / 4.0;
    chrisD[0][1][0][0] = 0.0;
    chrisD[0][1][0][1] = t91;
    chrisD[0][1][0][2] = 0.0;
    chrisD[0][1][0][3] = 0.0;
    chrisD[0][1][1][0] = t60;
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
    chrisD[0][2][0][2] = -t95;
    chrisD[0][2][0][3] = t99;
    chrisD[0][2][1][0] = 0.0;
    chrisD[0][2][1][1] = 0.0;
    chrisD[0][2][1][2] = 0.0;
    chrisD[0][2][1][3] = 0.0;
    chrisD[0][2][2][0] = t105;
    chrisD[0][2][2][1] = 0.0;
    chrisD[0][2][2][2] = 0.0;
    chrisD[0][2][2][3] = 0.0;
    chrisD[0][2][3][0] = t111;
    chrisD[0][2][3][1] = 0.0;
    chrisD[0][2][3][2] = 0.0;
    chrisD[0][2][3][3] = 0.0;
    chrisD[0][3][0][0] = 0.0;
    chrisD[0][3][0][1] = 0.0;
    chrisD[0][3][0][2] = t99;
    chrisD[0][3][0][3] = -t115;
    chrisD[0][3][1][0] = 0.0;
    chrisD[0][3][1][1] = 0.0;
    chrisD[0][3][1][2] = 0.0;
    chrisD[0][3][1][3] = 0.0;
    chrisD[0][3][2][0] = t118;
    chrisD[0][3][2][1] = 0.0;
    chrisD[0][3][2][2] = 0.0;
    chrisD[0][3][2][3] = 0.0;
    chrisD[0][3][3][0] = t127;
    chrisD[0][3][3][1] = 0.0;
    chrisD[0][3][3][2] = 0.0;
    chrisD[0][3][3][3] = 0.0;
    chrisD[1][0][0][0] = 0.0;
    chrisD[1][0][0][1] = t91;
    chrisD[1][0][0][2] = 0.0;
    chrisD[1][0][0][3] = 0.0;
    chrisD[1][0][1][0] = t60;
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
    chrisD[1][1][0][0] = -t91;
    chrisD[1][1][0][1] = 0.0;
    chrisD[1][1][0][2] = 0.0;
    chrisD[1][1][0][3] = 0.0;
    chrisD[1][1][1][0] = 0.0;
    chrisD[1][1][1][1] = t140 * t46 / 4.0;
    chrisD[1][1][1][2] = 0.0;
    chrisD[1][1][1][3] = 0.0;
    chrisD[1][1][2][0] = 0.0;
    chrisD[1][1][2][1] = 0.0;
    chrisD[1][1][2][2] = t149 / 4.0;
    chrisD[1][1][2][3] = t68 * t151 / 2.0;
    chrisD[1][1][3][0] = 0.0;
    chrisD[1][1][3][1] = 0.0;
    chrisD[1][1][3][2] = -t72 * t151 / 2.0;
    chrisD[1][1][3][3] = -t67 * t157 * t46 / 2.0;
    chrisD[1][2][0][0] = 0.0;
    chrisD[1][2][0][1] = 0.0;
    chrisD[1][2][0][2] = 0.0;
    chrisD[1][2][0][3] = 0.0;
    chrisD[1][2][1][0] = 0.0;
    chrisD[1][2][1][1] = 0.0;
    chrisD[1][2][1][2] = -t164;
    chrisD[1][2][1][3] = -t167;
    chrisD[1][2][2][0] = 0.0;
    chrisD[1][2][2][1] = t176;
    chrisD[1][2][2][2] = 0.0;
    chrisD[1][2][2][3] = 0.0;
    chrisD[1][2][3][0] = 0.0;
    chrisD[1][2][3][1] = -t181;
    chrisD[1][2][3][2] = 0.0;
    chrisD[1][2][3][3] = 0.0;
    chrisD[1][3][0][0] = 0.0;
    chrisD[1][3][0][1] = 0.0;
    chrisD[1][3][0][2] = 0.0;
    chrisD[1][3][0][3] = 0.0;
    chrisD[1][3][1][0] = 0.0;
    chrisD[1][3][1][1] = 0.0;
    chrisD[1][3][1][2] = -t167;
    chrisD[1][3][1][3] = -t184;
    chrisD[1][3][2][0] = 0.0;
    chrisD[1][3][2][1] = -t187;
    chrisD[1][3][2][2] = 0.0;
    chrisD[1][3][2][3] = 0.0;
    chrisD[1][3][3][0] = 0.0;
    chrisD[1][3][3][1] = -t193;
    chrisD[1][3][3][2] = 0.0;
    chrisD[1][3][3][3] = 0.0;
    chrisD[2][0][0][0] = 0.0;
    chrisD[2][0][0][1] = 0.0;
    chrisD[2][0][0][2] = -t95;
    chrisD[2][0][0][3] = t99;
    chrisD[2][0][1][0] = 0.0;
    chrisD[2][0][1][1] = 0.0;
    chrisD[2][0][1][2] = 0.0;
    chrisD[2][0][1][3] = 0.0;
    chrisD[2][0][2][0] = t105;
    chrisD[2][0][2][1] = 0.0;
    chrisD[2][0][2][2] = 0.0;
    chrisD[2][0][2][3] = 0.0;
    chrisD[2][0][3][0] = t111;
    chrisD[2][0][3][1] = 0.0;
    chrisD[2][0][3][2] = 0.0;
    chrisD[2][0][3][3] = 0.0;
    chrisD[2][1][0][0] = 0.0;
    chrisD[2][1][0][1] = 0.0;
    chrisD[2][1][0][2] = 0.0;
    chrisD[2][1][0][3] = 0.0;
    chrisD[2][1][1][0] = 0.0;
    chrisD[2][1][1][1] = 0.0;
    chrisD[2][1][1][2] = -t164;
    chrisD[2][1][1][3] = -t167;
    chrisD[2][1][2][0] = 0.0;
    chrisD[2][1][2][1] = t176;
    chrisD[2][1][2][2] = 0.0;
    chrisD[2][1][2][3] = 0.0;
    chrisD[2][1][3][0] = 0.0;
    chrisD[2][1][3][1] = -t181;
    chrisD[2][1][3][2] = 0.0;
    chrisD[2][1][3][3] = 0.0;
    chrisD[2][2][0][0] = -t194 * t93 / 2.0;
    chrisD[2][2][0][1] = 0.0;
    chrisD[2][2][0][2] = 0.0;
    chrisD[2][2][0][3] = 0.0;
    chrisD[2][2][1][0] = 0.0;
    chrisD[2][2][1][1] = t201 * t46 * t162 / 4.0;
    chrisD[2][2][1][2] = 0.0;
    chrisD[2][2][1][3] = 0.0;
    chrisD[2][2][2][0] = 0.0;
    chrisD[2][2][2][1] = 0.0;
    chrisD[2][2][2][2] = -t210 * t46 * t162 / 4.0;
    chrisD[2][2][2][3] = -t217;
    chrisD[2][2][3][0] = 0.0;
    chrisD[2][2][3][1] = 0.0;
    chrisD[2][2][3][2] = -t50 * t151 * t46 * t162 / 2.0;
    chrisD[2][2][3][3] = -t182 * t93 / 2.0;
    chrisD[2][3][0][0] = t226;
    chrisD[2][3][0][1] = 0.0;
    chrisD[2][3][0][2] = 0.0;
    chrisD[2][3][0][3] = 0.0;
    chrisD[2][3][1][0] = 0.0;
    chrisD[2][3][1][1] = t229;
    chrisD[2][3][1][2] = 0.0;
    chrisD[2][3][1][3] = 0.0;
    chrisD[2][3][2][0] = 0.0;
    chrisD[2][3][2][1] = 0.0;
    chrisD[2][3][2][2] = -t217;
    chrisD[2][3][2][3] = -t233;
    chrisD[2][3][3][0] = 0.0;
    chrisD[2][3][3][1] = 0.0;
    chrisD[2][3][3][2] = -t237;
    chrisD[2][3][3][3] = -t241;
    chrisD[3][0][0][0] = 0.0;
    chrisD[3][0][0][1] = 0.0;
    chrisD[3][0][0][2] = t99;
    chrisD[3][0][0][3] = -t115;
    chrisD[3][0][1][0] = 0.0;
    chrisD[3][0][1][1] = 0.0;
    chrisD[3][0][1][2] = 0.0;
    chrisD[3][0][1][3] = 0.0;
    chrisD[3][0][2][0] = t118;
    chrisD[3][0][2][1] = 0.0;
    chrisD[3][0][2][2] = 0.0;
    chrisD[3][0][2][3] = 0.0;
    chrisD[3][0][3][0] = t127;
    chrisD[3][0][3][1] = 0.0;
    chrisD[3][0][3][2] = 0.0;
    chrisD[3][0][3][3] = 0.0;
    chrisD[3][1][0][0] = 0.0;
    chrisD[3][1][0][1] = 0.0;
    chrisD[3][1][0][2] = 0.0;
    chrisD[3][1][0][3] = 0.0;
    chrisD[3][1][1][0] = 0.0;
    chrisD[3][1][1][1] = 0.0;
    chrisD[3][1][1][2] = -t167;
    chrisD[3][1][1][3] = -t184;
    chrisD[3][1][2][0] = 0.0;
    chrisD[3][1][2][1] = -t187;
    chrisD[3][1][2][2] = 0.0;
    chrisD[3][1][2][3] = 0.0;
    chrisD[3][1][3][0] = 0.0;
    chrisD[3][1][3][1] = -t193;
    chrisD[3][1][3][2] = 0.0;
    chrisD[3][1][3][3] = 0.0;
    chrisD[3][2][0][0] = t226;
    chrisD[3][2][0][1] = 0.0;
    chrisD[3][2][0][2] = 0.0;
    chrisD[3][2][0][3] = 0.0;
    chrisD[3][2][1][0] = 0.0;
    chrisD[3][2][1][1] = t229;
    chrisD[3][2][1][2] = 0.0;
    chrisD[3][2][1][3] = 0.0;
    chrisD[3][2][2][0] = 0.0;
    chrisD[3][2][2][1] = 0.0;
    chrisD[3][2][2][2] = -t217;
    chrisD[3][2][2][3] = -t233;
    chrisD[3][2][3][0] = 0.0;
    chrisD[3][2][3][1] = 0.0;
    chrisD[3][2][3][2] = -t237;
    chrisD[3][2][3][3] = -t241;
    chrisD[3][3][0][0] = -t245 * t46 * t113 / 4.0;
    chrisD[3][3][0][1] = 0.0;
    chrisD[3][3][0][2] = 0.0;
    chrisD[3][3][0][3] = 0.0;
    chrisD[3][3][1][0] = 0.0;
    chrisD[3][3][1][1] = t194 * t96 / 2.0;
    chrisD[3][3][1][2] = 0.0;
    chrisD[3][3][1][3] = 0.0;
    chrisD[3][3][2][0] = 0.0;
    chrisD[3][3][2][1] = 0.0;
    chrisD[3][3][2][2] = -t92 * t96 / 2.0;
    chrisD[3][3][2][3] = t67 * t69 * t46 * t113 / 2.0;
    chrisD[3][3][3][0] = 0.0;
    chrisD[3][3][3][1] = 0.0;
    chrisD[3][3][3][2] = -t241;
    chrisD[3][3][3][3] = t262 * t46 * t113 / 4.0;


    return true;
}

/*! Transform local 4-direction to coordinate 4-direction.
 *
 *  \param  pos  :  pointer to position array.
 *  \param  ldir :  pointer to local direction array.
 *  \param  dir  :  pointer to calculated coordinate direction array.
 *  \param  type :  type of tetrad.
 */
void MetricPTD_C::localToCoord(const double* pos, const double* ldir, double* dir,
                               enum_nat_tetrad_type) {
    double x = pos[2];
    double y = pos[3];
    double a = Par_a;
    double b = Par_b;

    dir[0] = ldir[0] / (sqrt(fabs(-y * y * y - y * a + b)) / (x + y));
    dir[1] = ldir[1] / (sqrt(fabs(x * x * x + x * a + b)) / (x + y));
    dir[2] = ldir[2] * ((x + y) * sqrt(fabs(x * x * x + x * a + b)));
    dir[3] = ldir[3] * ((x + y) * sqrt(fabs(-y * y * y - y * a + b)));
}

/*! Transform coordinate 4-direction to local 4-direction.
 *
 *  \param  pos  :  pointer to position array.
 *  \param  cdir :  pointer to coordinate direction.
 *  \param  ldir :  pointer to calculated local direction array.
 *  \param  type :  type of tetrad.
 */
void MetricPTD_C::coordToLocal(const double* pos, const double* cdir, double* ldir,
                               enum_nat_tetrad_type) {
    double x = pos[2];
    double y = pos[3];
    double a = Par_a;
    double b = Par_b;

    ldir[0] = cdir[0] * sqrt(fabs(-y * y * y - y * a + b)) / (x + y);
    ldir[1] = cdir[1] * sqrt(fabs(x * x * x + x * a + b)) / (x + y);
    ldir[2] = cdir[2] / ((x + y) * sqrt(fabs(x * x * x + x * a + b)));
    ldir[3] = cdir[3] / ((x + y) * sqrt(fabs(-y * y * y - y * a + b)));
}

/*! Transform point p to custom coordinats.
 *
 * Deleted.
 */
int MetricPTD_C::transToPseudoCart(vec4 p, vec4& cp) {
    cp[0] = p[0];
    cp[1] = p[2];
    cp[2] = p[3];
    cp[3] = p[1];
    return 0;
}

/*! Tests break condition
 *  \param pos  :  position.
 *
 * not implemented.
 */
bool MetricPTD_C::breakCondition(const double*) {
    bool br = false;
    /*
      double x = pos[2];
      double y = pos[3];
      double a = Par_a;
      double b = Par_b;
    */
    /* if (abs(x+y)<eps) { return true; }

     double t = x;
     if (abs(t*t*t+a*t+b)<eps)  { return true; }
     t= -y;
     if (abs(t*t*t+a*t+b)<eps)  { return true; }
    */

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
double MetricPTD_C::testConstraint(const double y[], const double kappa) {
    calculateMetric(y);
    double sum = -mSign * kappa;
    for (int i = 0; i < 4; i++) {
        sum += g_compts[i][i] * y[4 + i] * y[4 + i];
    }
    return sum;
    //return 0;
}

/*! Set parameter 'pName' to 'val'.
 *
 *
 */
bool MetricPTD_C::setParam(std::string pName, double val) {
    Metric::setParam(pName, val);
    if (pName == "b") {
        Par_b = val;
    }
    if (pName == "a") {
        Par_a = val;
    }
    return true;
}




/*! Generate report.
 */
bool MetricPTD_C::report(const vec4 , const vec4 , std::string &text) {
    std::stringstream ss;
    ss << "Report for C metric\n\tcoordinates : (t,phi,x,y)\n";
    ss << "---------------------------------------------------------------\n";
    ss << "Coordinate Ranges:" << std::endl;
    ss << "     x, y: ( 0 < (x+y) )" << std::endl;
    ss << "        x: ( 0 < f( x) )" << std::endl;
    ss << "        y: ( 0 < f(-y) )" << std::endl;
    ss << "      phi: arbitrary" << std::endl;
    ss << "---------------------------------------------------------------\n";
    ss << "Parameter Ranges:" << std::endl;
    ss << "        a: arbitrary (a= -2^(1/3)*m^(-4/3)/6  of Pravda C Parameters)" << std::endl;
    ss << "        b: arbitrary (b= -1/(18*m^2) + A^2    of Pravda C Parameters)" << std::endl;
    ss << "---------------------------------------------------------------\n";
    ss.precision(DEF_FIXED_REPORT_PRECISION);
    ss.setf(std::ios::fixed);


    vec3 roots;

    calculateRoots(roots, -Par_a, Par_b);

    ss << "  Root(s) of Polynom f(x)= x^3 + a*x + b:" << std::endl;
    ss << "    f negative before........... x1 = " << roots[0] << std::endl;
    ss << "    f positive between... x1 and x2 = " << roots[1] << std::endl;
    ss << "    f negative between... x2 and x3 = " << roots[2] << std::endl;
    ss << "    f positive after x3." << std::endl;

    ss << "  Root(s) of Polynom f(-y):" << std::endl;
    ss << "    f positive before........... y1 = " << -roots[2] << std::endl;
    ss << "    f negative between... y1 and y2 = " << -roots[1] << std::endl;
    ss << "    f positive between... y2 and y3 = " << -roots[0] << std::endl;
    ss << "    f negative after y3." << std::endl;
    ss << "---------------------------------------------------------------\n";
    ss << "quad2d    " << roots[1] << " " << -roots[1] << " " << roots[1] << " " << -roots[0] << " " << roots[2] << " " << -roots[0] << " " << roots[2] << " " << -roots[1] << " 3.0   0.5 1.0 0.5" << std::endl;

    text = ss.str();
    return true;
}

// *************************** specific  public methods ****************************
/*! Returns the sign of real number.
 *
 */
double sgn(double value) {
    if (value > 0) {
        return 1.0;
    }
    if (value < 0) {
        return -1.0;
    }
    return 0;
}

/*! Calculates the roots of a polynom of type x^3 - px + q = 0.
 *
 *\param   p     : coefficient of the linear term.
 *\param   q     : y-axes offset
 *\param roots : reference to Roots of the polynom sorted ascending.
 */
void
MetricPTD_C::calculateRoots(vec3 & roots, double p, double q) {
    double d = q * q / 4.0 - p * p * p / 27.0;
    double z1 = 0;
    double z2 = 0;
    double z3 = 0;

    if (d < 0) {
        double u = acos(-q * 0.5 * sqrt(27.0 / p / p / p));
        z1 = sqrt(4.0 / 3.0 * p) * cos((u) / 3.0);
        z2 = -sqrt(4.0 / 3.0 * p) * cos((u + M_PI) / 3.0);
        z3 = -sqrt(4.0 / 3.0 * p) * cos((u - M_PI) / 3.0);
    }
    if (d == 0) {
        if (p == 0) {
            z1 = 0.0;
            z2 = 0.0;
        } else {
            z1 = -3 * q / p;
            z2 = 1.5 * q / p;
        }
        z3 = z2;
    }
    if (d > 0) {
        z1 = sgn(-q * 0.5 + sqrt(d)) * pow(fabs(-q * 0.5 + sqrt(d)), 1 / 3.0) -  sgn(q * 0.5 + sqrt(d)) * pow(fabs(q * 0.5 + sqrt(d)), 1 / 3.0);
        z3 = z2 = z1;
    }

    double tmp;
    //sort
    if (z2 < z1) {
        tmp = z1;
        z1 = z2;
        z2 = tmp;
    }
    if (z3 < z2) {
        tmp = z3;
        z3 = z2;
        z2 = tmp;
    }
    if (z2 < z1) {
        tmp = z1;
        z1 = z2;
        z2 = tmp;
    }
    roots[0] = z1;
    roots[1] = z2;
    roots[2] = z3;
}

// ********************************* protected methods *****************************
/*!
 */
void MetricPTD_C::setStandardValues() {
    mInitPos[0] = 0.0;
    mInitPos[1] = 0.0;
    mInitPos[2] = 1.0;
    mInitPos[3] = -0.5;
    mInitDir[0] = 1.0;
    mInitDir[1] = 0.0;
    mInitDir[2] = 0.0;

    mCoordNames[0] = std::string("t");
    mCoordNames[1] = std::string("phi");
    mCoordNames[2] = std::string("x");
    mCoordNames[3] = std::string("y");
}


} // end namespace m4d

