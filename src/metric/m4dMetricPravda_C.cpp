// -------------------------------------------------------------------------------
/*
   m4dMetricPravda_C.cpp

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

#include "m4dMetricPravda_C.h"

namespace m4d {

#define eps 1.0e-9


/*! Standard constructor for the metric.
 *
 * \param  m : mass of the black holes.
 * \param  A : acceleration constant of the black holes.
 */
MetricPravda_C::MetricPravda_C(double A, double m) {
    mMetricName  = "Pravda_C-Metric";
    //setCoordType(enum_coordinate_custom);
    setCoordType(enum_coordinate_cartesian);

    mPhysicalUnits = enum_physical_constants_geom;
    mSpeedOfLight = 1.0;
    mGravConstant = 1.0;


    Par_A = A;
    Par_m = m;

    addParam("a", Par_A);
    addParam("m", Par_m);

    setStandardValues();


}

/*! Standard destructor for the metric.
 *
 */
MetricPravda_C::~MetricPravda_C() {

}


// *********************************** public methods ******************************
/*! Calculate the contravariant metric components at position 'pos'.
 *
 *  \param pos : pointer to position.
 */
bool MetricPravda_C::calculateMetric(const double* pos) {
    double x = pos[1];
    double y = pos[2];
    double A = Par_A;
    double m = Par_m;

    double t1 = A * A;
    double  t4 = pow(x + y, 2.0);
    double  t6 = 1 / t1 / t4;
    double  t7 = y * y;
    double  t8 = m * A;
    double  t12 = -1.0 + t7 - 2.0 * t8 * t7 * y;
    double  t14 = x * x;
    double  t18 = 1.0 - t14 - 2.0 * t8 * t14 * x;
    g_compts[0][0] = -t6 * t12;
    g_compts[0][1] = 0.0;
    g_compts[0][2] = 0.0;
    g_compts[0][3] = 0.0;
    g_compts[1][0] = 0.0;
    g_compts[1][1] = t6 / t18;
    g_compts[1][2] = 0.0;
    g_compts[1][3] = 0.0;
    g_compts[2][0] = 0.0;
    g_compts[2][1] = 0.0;
    g_compts[2][2] = t6 / t12;
    g_compts[2][3] = 0.0;
    g_compts[3][0] = 0.0;
    g_compts[3][1] = 0.0;
    g_compts[3][2] = 0.0;
    g_compts[3][3] = t6 * t18;





    return true;
}

/*! Calculate the Christoffel symbols of the second kind at position 'pos'.
 *
 *  \param pos : pointer to position.
 */
bool MetricPravda_C::calculateChristoffels(const double* pos) {
    double x = pos[1];
    double y = pos[2];
    double A = Par_A;
    double m = Par_m;

    double  t2 = 1 / (x + y);
    double  t3 = x * x;
    double  t4 = m * A;
    double  t6 = t4 * t3 * x;
    double  t8 = -1.0 + t3 + 2.0 * t6;
    double  t9 = t2 * t8;
    double  t10 = y * y;
    double  t12 = t4 * t10 * y;
    double  t14 = 1.0 - t10 + 2.0 * t12;
    double  t15 = t9 * t14;
    double  t16 = t2 * t14;
    double  t17 = y * x;
    double  t20 = 3.0 * t4 * t10 * x;
    double  t21 = -1.0 + t12 - t17 + t20;
    double  t23 = 1 / t14;
    double  t24 = t2 * t23;
    double  t25 = t24 * t21;
    double  t26 = 1 / t8;
    double  t27 = t2 * t26;
    double  t32 = 3.0 * t4 * y * t3;
    double  t36 = 1.0 + t6 + t17 + t32;
    double  t37 = t27 * t36;
    christoffel[0][0][0] = 0.0;
    christoffel[0][0][1] = -t15;
    christoffel[0][0][2] = t16 * t21;
    christoffel[0][0][3] = 0.0;
    christoffel[0][1][0] = -t2;
    christoffel[0][1][1] = 0.0;
    christoffel[0][1][2] = 0.0;
    christoffel[0][1][3] = 0.0;
    christoffel[0][2][0] = t25;
    christoffel[0][2][1] = 0.0;
    christoffel[0][2][2] = 0.0;
    christoffel[0][2][3] = 0.0;
    christoffel[0][3][0] = 0.0;
    christoffel[0][3][1] = 0.0;
    christoffel[0][3][2] = 0.0;
    christoffel[0][3][3] = 0.0;
    christoffel[1][0][0] = -t2;
    christoffel[1][0][1] = 0.0;
    christoffel[1][0][2] = 0.0;
    christoffel[1][0][3] = 0.0;
    christoffel[1][1][0] = 0.0;
    christoffel[1][1][1] = -t27 * (-1.0 + 2.0 * t3 + 5.0 * t6 + t17 + t32);
    christoffel[1][1][2] = t16 * t26;
    christoffel[1][1][3] = 0.0;
    christoffel[1][2][0] = 0.0;
    christoffel[1][2][1] = -t2;
    christoffel[1][2][2] = -t2;
    christoffel[1][2][3] = 0.0;
    christoffel[1][3][0] = 0.0;
    christoffel[1][3][1] = 0.0;
    christoffel[1][3][2] = 0.0;
    christoffel[1][3][3] = t37;
    christoffel[2][0][0] = t25;
    christoffel[2][0][1] = 0.0;
    christoffel[2][0][2] = 0.0;
    christoffel[2][0][3] = 0.0;
    christoffel[2][1][0] = 0.0;
    christoffel[2][1][1] = -t2;
    christoffel[2][1][2] = -t2;
    christoffel[2][1][3] = 0.0;
    christoffel[2][2][0] = 0.0;
    christoffel[2][2][1] = t9 * t23;
    christoffel[2][2][2] = -t24 * (1.0 - 2.0 * t10 + 5.0 * t12 - t17 + t20);
    christoffel[2][2][3] = 0.0;
    christoffel[2][3][0] = 0.0;
    christoffel[2][3][1] = 0.0;
    christoffel[2][3][2] = 0.0;
    christoffel[2][3][3] = -t2;
    christoffel[3][0][0] = 0.0;
    christoffel[3][0][1] = 0.0;
    christoffel[3][0][2] = 0.0;
    christoffel[3][0][3] = 0.0;
    christoffel[3][1][0] = 0.0;
    christoffel[3][1][1] = 0.0;
    christoffel[3][1][2] = 0.0;
    christoffel[3][1][3] = t37;
    christoffel[3][2][0] = 0.0;
    christoffel[3][2][1] = 0.0;
    christoffel[3][2][2] = 0.0;
    christoffel[3][2][3] = -t2;
    christoffel[3][3][0] = 0.0;
    christoffel[3][3][1] = -t9 * t36;
    christoffel[3][3][2] = t15;
    christoffel[3][3][3] = 0.0;

    return true;
}

/*! Calculate Jacobi matrix.
 *
 *  \param pos : pointer to position.
 */
bool
MetricPravda_C::calculateChrisD(const double* pos) {
    double x = pos[1];
    double y = pos[2];
    double A = Par_A;
    double m = Par_m;

    double t1 = y * y;
    double  t2 = m * A;
    double  t3 = t1 * y;
    double  t4 = t2 * t3;
    double  t5 = 2.0 * t4;
    double  t6 = 1.0 - t1 + t5;
    double  t7 = x * x;
    double  t8 = t7 * x;
    double  t9 = t2 * t8;
    double  t12 = 2.0 * y * x;
    double  t14 = t2 * y * t7;
    double  t15 = 6.0 * t14;
    double  t16 = 1.0 + t7 + 4.0 * t9 + t12 + t15;
    double  t19 = pow(x + y, 2.0);
    double  t20 = 1 / t19;
    double  t21 = t6 * t16 * t20;
    double  t22 = 2.0 * t9;
    double  t23 = -1.0 + t7 + t22;
    double  t26 = t2 * t1 * x;
    double  t27 = 6.0 * t26;
    double  t28 = -1.0 - t1 + 4.0 * t4 - t12 + t27;
    double  t30 = t23 * t28 * t20;
    double  t31 = t6 * t6;
    double  t33 = m * m;
    double  t34 = A * A;
    double  t35 = t33 * t34;
    double  t36 = t1 * t1;
    double  t38 = t35 * t36 * t1;
    double  t39 = 10.0 * t38;
    double  t40 = t36 * y;
    double  t42 = t40 * m * A;
    double  t46 = t36 * m * A * x;
    double  t49 = t35 * t40 * x;
    double  t52 = t2 * t3 * t7;
    double  t55 = t35 * t36 * t7;
    double  t58 = 2.0 * t3 * x;
    double  t59 = t1 * t7;
    double  t60 = 3.0 * t59;
    double  t61 = 1.0 + t1 - t7 - t5 + t15 + t39 + t12 - 4.0 * t42 - 20.0 * t46 + 36.0 * t49 - 20.0 * t52 + 30.0 * t55 + t58 + t60;
    double  t63 = 3.0 * t1;
    double  t65 = 12.0 * t26;
    double  t67 = 8.0 * t46;
    double  t68 = 12.0 * t49;
    double  t69 = 4.0 * t52;
    double  t70 = 6.0 * t55;
    double  t71 = -1.0 + t63 + t7 - 10.0 * t4 - t65 - t15 + 2.0 * t38 + t12 - t67 + t68 - t69 + t70 + t58 + t59;
    double  t73 = 1 / t31;
    double  t74 = t71 * t20 * t73;
    double  t75 = 12.0 * t14;
    double  t76 = t7 * t7;
    double  t77 = t76 * x;
    double  t79 = t77 * m * A;
    double  t82 = t35 * t76 * t7;
    double  t83 = 10.0 * t82;
    double  t86 = y * t76 * t2;
    double  t87 = 8.0 * t86;
    double  t89 = t35 * y * t77;
    double  t90 = 12.0 * t89;
    double  t92 = t2 * t1 * t8;
    double  t93 = 4.0 * t92;
    double  t95 = t35 * t1 * t76;
    double  t96 = 6.0 * t95;
    double  t98 = 2.0 * y * t8;
    double  t99 = 1.0 + t1 - t7 + t22 + t27 + t75 + t12 + 8.0 * t79 + t83 + 2.0 * t76 + t87 + t90 + t93 + t96 + t59 + t98;
    double  t101 = t23 * t23;
    double  t102 = 1 / t101;
    double  t104 = 3.0 * t7;
    double  t115 = -1.0 + t1 + t104 + 10.0 * t9 + t27 + t75 + t12 + 2.0 * t82 + t87 + t90 + t93 + t96 + t59 + t98;
    double  t117 = t115 * t20 * t102;
    double  t128 = 1.0 - t1 + t7 - t5 - t65 - t15 + t39 + t12 - 8.0 * t42 - t67 + t68 - t69 + t70 + t58 + t59 + 2.0 * t36;
    double  t136 = 1.0 - t1 + t7 + t22 - t27 + t12 + 4.0 * t79 + t83 + 20.0 * t86 + 36.0 * t89 + 20.0 * t92 + 30.0 * t95 + t60 + t98;
    chrisD[0][0][0][0] = 0.0;
    chrisD[0][0][0][1] = 0.0;
    chrisD[0][0][0][2] = 0.0;
    chrisD[0][0][0][3] = 0.0;
    chrisD[0][0][1][0] = 0.0;
    chrisD[0][0][1][1] = -t21;
    chrisD[0][0][1][2] = -t30;
    chrisD[0][0][1][3] = 0.0;
    chrisD[0][0][2][0] = 0.0;
    chrisD[0][0][2][1] = t31 * t20;
    chrisD[0][0][2][2] = t61 * t20;
    chrisD[0][0][2][3] = 0.0;
    chrisD[0][0][3][0] = 0.0;
    chrisD[0][0][3][1] = 0.0;
    chrisD[0][0][3][2] = 0.0;
    chrisD[0][0][3][3] = 0.0;
    chrisD[0][1][0][0] = 0.0;
    chrisD[0][1][0][1] = t20;
    chrisD[0][1][0][2] = t20;
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
    chrisD[0][2][0][1] = t20;
    chrisD[0][2][0][2] = -t74;
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
    chrisD[1][0][0][1] = t20;
    chrisD[1][0][0][2] = t20;
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
    chrisD[1][1][1][1] = t99 * t20 * t102;
    chrisD[1][1][1][2] = t20;
    chrisD[1][1][1][3] = 0.0;
    chrisD[1][1][2][0] = 0.0;
    chrisD[1][1][2][1] = -t6 * (-1.0 + t104 + 8.0 * t9 + t12 + t15) * t20 * t102;
    chrisD[1][1][2][2] = t28 * t20 / t23;
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
    chrisD[1][2][1][1] = t20;
    chrisD[1][2][1][2] = t20;
    chrisD[1][2][1][3] = 0.0;
    chrisD[1][2][2][0] = 0.0;
    chrisD[1][2][2][1] = t20;
    chrisD[1][2][2][2] = t20;
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
    chrisD[1][3][3][1] = -t117;
    chrisD[1][3][3][2] = t20;
    chrisD[1][3][3][3] = 0.0;
    chrisD[2][0][0][0] = 0.0;
    chrisD[2][0][0][1] = t20;
    chrisD[2][0][0][2] = -t74;
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
    chrisD[2][1][1][1] = t20;
    chrisD[2][1][1][2] = t20;
    chrisD[2][1][1][3] = 0.0;
    chrisD[2][1][2][0] = 0.0;
    chrisD[2][1][2][1] = t20;
    chrisD[2][1][2][2] = t20;
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
    chrisD[2][2][1][1] = t16 * t20 / t6;
    chrisD[2][2][1][2] = -t23 * (1.0 - t63 + 8.0 * t4 - t12 + t27) * t20 * t73;
    chrisD[2][2][1][3] = 0.0;
    chrisD[2][2][2][0] = 0.0;
    chrisD[2][2][2][1] = t20;
    chrisD[2][2][2][2] = t128 * t20 * t73;
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
    chrisD[2][3][3][1] = t20;
    chrisD[2][3][3][2] = t20;
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
    chrisD[3][1][3][1] = -t117;
    chrisD[3][1][3][2] = t20;
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
    chrisD[3][2][3][1] = t20;
    chrisD[3][2][3][2] = t20;
    chrisD[3][2][3][3] = 0.0;
    chrisD[3][3][0][0] = 0.0;
    chrisD[3][3][0][1] = 0.0;
    chrisD[3][3][0][2] = 0.0;
    chrisD[3][3][0][3] = 0.0;
    chrisD[3][3][1][0] = 0.0;
    chrisD[3][3][1][1] = -t136 * t20;
    chrisD[3][3][1][2] = -t101 * t20;
    chrisD[3][3][1][3] = 0.0;
    chrisD[3][3][2][0] = 0.0;
    chrisD[3][3][2][1] = t21;
    chrisD[3][3][2][2] = t30;
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
void MetricPravda_C::localToCoord(const double* pos, const double* ldir, double* dir,
                                  enum_nat_tetrad_type) {
    double x = pos[1];
    double y = pos[2];
    double A = Par_A;
    double m = Par_m;

    dir[0] = ldir[0] / sqrt(fabs(-2 * m * A * y * y * y + y * y - 1)) * A * (x + y); //q
    dir[1] = ldir[1] * sqrt(fabs(-2 * m * A * x * x * x - x * x + 1)) * A * (x + y); //x
    dir[2] = ldir[2] * (sqrt(fabs(-2 * m * A * y * y * y + y * y - 1))) * A * (x + y); //y
    dir[3] = ldir[3] / (sqrt(fabs(-2 * m * A * x * x * x - x * x + 1))) * A * (x + y); // G^-1 /A^2/(x+y)^2  //p
}

/*! Transform coordinate 4-direction to local 4-direction.
 *
 *  \param  pos  :  pointer to position array.
 *  \param  cdir :  pointer to coordinate direction.
 *  \param  ldir :  pointer to calculated local direction array.
 *  \param  type :  type of tetrad.
 */
void MetricPravda_C::coordToLocal(const double* pos, const double* cdir, double* ldir,
                                  enum_nat_tetrad_type) {
    double x = pos[1];
    double y = pos[2];
    double A = Par_A;
    double m = Par_m;

    ldir[0] = cdir[1] * sqrt(fabs(-2 * m * A * y * y * y + y * y - 1)) / (A * (x + y));
    ldir[1] = cdir[1] / sqrt(fabs(-2 * m * A * x * x * x - x * x + 1)) / (A * (x + y));
    ldir[2] = cdir[2] / (sqrt(fabs(-2 * m * A * y * y * y + y * y - 1)) * A * (x + y));
    ldir[3] = cdir[3] * (sqrt(fabs(-2 * m * A * x * x * x - x * x + 1))) / (A * (x + y));

}

/*! Transform point p to custom coordinates.
 *
 *  \param p  : point to be transformed.
 *  \param cp : reference to customized point.
 *  \return true : always.
 *
 * Function deleted: Could transform metric to cylindrical coordinates.
 */
int MetricPravda_C::transToPseudoCart(vec4 p, vec4& cp) {
    //ToDo
    /* double phi = p[3];
     double x   = p[1];
     double y   = p[2];
     double A = Par_A;
     double m = Par_m;
     double rho = sqrt((-2*m*A*x*x*x-x*x+1)*(-2*m*A*y*y*y+y*y-1))/A/A/(x+y)/(x+y);

     cp[0] = p[0];
     cp[1] = rho*cos(phi);
     cp[2] = rho*sin(phi);
     cp[3] = (1+m*A*x*y*(x-y)+x*y)/A/A/(x+y)/(x+y);*/
    cp[0] = p[0];
    cp[1] = p[1];
    cp[2] = p[2];
    cp[3] = p[3];
    return 0;
}

/*! Tests break condition.
 *  \param pos  :  position.
 *
 *  not implemented.
 */
bool MetricPravda_C::breakCondition(const double*) {
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
double MetricPravda_C::testConstraint(const double y[], const double kappa) {
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
bool MetricPravda_C::setParam(const char* pName, double val) {
    Metric::setParam(pName, val);
    if (pName == "m") {
        Par_m = val;
    }
    if (pName == "a") {
        Par_A = val;
    }
    return true;
}

/*! Generate report.
 */
bool MetricPravda_C::report(const vec4 , const vec4 , std::string &text) {
    std::stringstream ss;
    ss << "Report for Pravda C metric)\n\tcoordinates : (q,x,y,p)\n";
    ss << "---------------------------------------------------------------\n";
    ss << "Coordinate Ranges:" << std::endl;
    ss << "     x, y: ( 0 < (x+y) )" << std::endl;
    ss << "        x: ( 0 < F(y) )" << std::endl;
    ss << "        y: ( 0 < G(x) )" << std::endl;
    ss << "        p: arbitrary" << std::endl;
    ss << "---------------------------------------------------------------\n";
    ss << "Parameter Range:" << std::endl;
    ss << "        A: ??? " << std::endl;
    ss << "        m: ( m > 0 )" << std::endl;
    ss << "     A, m: 27m^2A^2<1????" << std::endl;
    ss << "---------------------------------------------------------------\n";
    ss.precision(DEF_FIXED_REPORT_PRECISION);
    ss.setf(std::ios::fixed);

    double A = Par_A;
    double m = Par_m;

    vec3 roots;

    calculateRoots(roots, 0.5 / m / A, -0.5 / m / A);

    ss << "  Roots of Polynom G=  1 - x^2 - 2mAx^3:" << std::endl;
    ss << "    G positive before........... x1 = " << roots[0] << std::endl;
    ss << "    G negative between... x1 and x2 = " << roots[1] << std::endl;
    ss << "    G positive between... x2 and x3 = " << roots[2] << std::endl;
    ss << "    G negative after x3." << std::endl;

    //calculateRoots(roots, -0.5/m/A, 0.5/m/A);

    ss << "  Roots of Polynom F= -1 + y^2 - 2mAy^3:" << std::endl;
    ss << "    F positive before........... y1 = " << -roots[2] << std::endl;
    ss << "    F negative between... y1 and y2 = " << -roots[1] << std::endl;
    ss << "    F positive between... y2 and y3 = " << -roots[0] << std::endl;
    ss << "    F negative after y3." << std::endl;
    ss << "---------------------------------------------------------------\n";
    ss << "quad2d    " << roots[1] << " " << -roots[1] << " " << roots[1] << " " << -roots[0] << " " << roots[2] << " " << -roots[0] << " " << roots[2] << " " << -roots[1] << " 3.0   0.5 1.0 0.5" << std::endl;


    text = ss.str();
    return true;
}

// *************************** specific  public methods ****************************
/*! Calculates the roots of a polynom of type x^3 + ax^2 + c = 0.
 *
 *\param   a     : coefficient of the quadratic term.
 *\param   c     : y-axes offset
 *\param roots : reference Roots of the polynom sorted ascending.
 */
void
MetricPravda_C::calculateRoots(vec3 & roots, double a, double c) {

    double p = a * a / 3.0;
    double q = 2.0 * a * a * a / 27.0 + c;
    double u = acos(-q * 0.5 * sqrt(27.0 / p / p / p));

    double z1 =  sqrt(4.0 / 3.0 * p) * cos((u) / 3.0) - a / 3.0;
    double z2 = -sqrt(4.0 / 3.0 * p) * cos((u + M_PI) / 3.0) - a / 3.0;
    double z3 = -sqrt(4.0 / 3.0 * p) * cos((u - M_PI) / 3.0) - a / 3.0;

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
void MetricPravda_C::setStandardValues() {
    mInitPos[0] = 0.0;
    mInitPos[1] = -5.0;
    mInitPos[2] = -5.0;
    mInitPos[3] = 0.0;
    mInitDir[0] = 1.0;
    mInitDir[1] = 0.0;
    mInitDir[2] = 0.0;

    mCoordNames[0] = std::string("q");
    mCoordNames[1] = std::string("x");
    mCoordNames[2] = std::string("y");
    mCoordNames[3] = std::string("p");
}


} // end namespace m4d

