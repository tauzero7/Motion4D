// -------------------------------------------------------------------------------
/*
   m4dMetricExtremeReissnerNordstromDihole.cpp

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

#include "m4dMetricExtremeReissnerNordstromDihole.h"

namespace m4d {

#define eps 1.0e-6


/*! Standard constructor for the Kottler metric.
 *
 * \param  mass1 : mass of the first black hole.
 * \param  mass2 : mass of the second black hole.
 */
MetricExtremeReissnerNordstromDihole::MetricExtremeReissnerNordstromDihole(double mass1, double mass2) {
    mMetricName  = "ExtremeReissnerNordstromDihole";
    setCoordType(enum_coordinate_cartesian);

    mPhysicalUnits = enum_physical_constants_geom;
    mSpeedOfLight = 1.0;
    mGravConstant = 1.0;

    addParam("m1", mass1);
    addParam("m2", mass2);
    mM1 = mass1;
    mM2 = mass2;

    mDrawTypes.push_back(enum_draw_twoplusone);

    mLocTeds.push_back(enum_nat_tetrad_static);

    w = gsl_integration_workspace_alloc(1000);
    F.function = &integrandForPeriodicAlongX;

    setStandardValues();
}

MetricExtremeReissnerNordstromDihole::~MetricExtremeReissnerNordstromDihole() {
    gsl_integration_workspace_free(w);
}


// *********************************** public methods ******************************
/*! Calculate the contravariant metric components at position 'pos'.
 *
 *  \param pos : pointer to position.
 */
bool MetricExtremeReissnerNordstromDihole::calculateMetric(const double* pos) {
    double t1 = calc_U(pos);  // U(x,y,z);
    double t2 = t1 * t1;

    g_compts[0][0] = -1 / t2;
    g_compts[0][1] = 0.0;
    g_compts[0][2] = 0.0;
    g_compts[0][3] = 0.0;
    g_compts[1][0] = 0.0;
    g_compts[1][1] = t2;
    g_compts[1][2] = 0.0;
    g_compts[1][3] = 0.0;
    g_compts[2][0] = 0.0;
    g_compts[2][1] = 0.0;
    g_compts[2][2] = t2;
    g_compts[2][3] = 0.0;
    g_compts[3][0] = 0.0;
    g_compts[3][1] = 0.0;
    g_compts[3][2] = 0.0;
    g_compts[3][3] = t2;

    return true;
}

/*! Calculate the Christoffel symbols of the second kind at position 'pos'.
 *
 *  \param pos : pointer to position.
 */
bool MetricExtremeReissnerNordstromDihole::calculateChristoffels(const double* pos) {
    double t1 = calc_U(pos);  // U(x,y,z);
    double Ux, Uy, Uz;
    calc_dU(pos, Ux, Uy, Uz);

    double t2 = t1 * t1;
    double t3 = t2 * t2;
    double t5 = 1 / t3 / t1;
    double t6 = Ux;   // diff(U(x,y,z),x);
    double t8 = Uy;   // diff(U(x,y,z),y);
    double t10 = Uz;  // diff(U(x,y,z),z);
    double t12 = 1 / t1;
    double t13 = t12 * t6;
    double t14 = t12 * t8;
    double t15 = t12 * t10;

    christoffel[0][0][0] = 0.0;
    christoffel[0][0][1] = -t5 * t6;
    christoffel[0][0][2] = -t5 * t8;
    christoffel[0][0][3] = -t5 * t10;
    christoffel[0][1][0] = -t13;
    christoffel[0][1][1] = 0.0;
    christoffel[0][1][2] = 0.0;
    christoffel[0][1][3] = 0.0;
    christoffel[0][2][0] = -t14;
    christoffel[0][2][1] = 0.0;
    christoffel[0][2][2] = 0.0;
    christoffel[0][2][3] = 0.0;
    christoffel[0][3][0] = -t15;
    christoffel[0][3][1] = 0.0;
    christoffel[0][3][2] = 0.0;
    christoffel[0][3][3] = 0.0;
    christoffel[1][0][0] = -t13;
    christoffel[1][0][1] = 0.0;
    christoffel[1][0][2] = 0.0;
    christoffel[1][0][3] = 0.0;
    christoffel[1][1][0] = 0.0;
    christoffel[1][1][1] = t13;
    christoffel[1][1][2] = -t14;
    christoffel[1][1][3] = -t15;
    christoffel[1][2][0] = 0.0;
    christoffel[1][2][1] = t14;
    christoffel[1][2][2] = t13;
    christoffel[1][2][3] = 0.0;
    christoffel[1][3][0] = 0.0;
    christoffel[1][3][1] = t15;
    christoffel[1][3][2] = 0.0;
    christoffel[1][3][3] = t13;
    christoffel[2][0][0] = -t14;
    christoffel[2][0][1] = 0.0;
    christoffel[2][0][2] = 0.0;
    christoffel[2][0][3] = 0.0;
    christoffel[2][1][0] = 0.0;
    christoffel[2][1][1] = t14;
    christoffel[2][1][2] = t13;
    christoffel[2][1][3] = 0.0;
    christoffel[2][2][0] = 0.0;
    christoffel[2][2][1] = -t13;
    christoffel[2][2][2] = t14;
    christoffel[2][2][3] = -t15;
    christoffel[2][3][0] = 0.0;
    christoffel[2][3][1] = 0.0;
    christoffel[2][3][2] = t15;
    christoffel[2][3][3] = t14;
    christoffel[3][0][0] = -t15;
    christoffel[3][0][1] = 0.0;
    christoffel[3][0][2] = 0.0;
    christoffel[3][0][3] = 0.0;
    christoffel[3][1][0] = 0.0;
    christoffel[3][1][1] = t15;
    christoffel[3][1][2] = 0.0;
    christoffel[3][1][3] = t13;
    christoffel[3][2][0] = 0.0;
    christoffel[3][2][1] = 0.0;
    christoffel[3][2][2] = t15;
    christoffel[3][2][3] = t14;
    christoffel[3][3][0] = 0.0;
    christoffel[3][3][1] = -t13;
    christoffel[3][3][2] = -t14;
    christoffel[3][3][3] = t15;

    return true;
}

/*! Calculate Jacobi matrix.
 *
 *  \param pos : pointer to position.
 */
bool MetricExtremeReissnerNordstromDihole::calculateChrisD(const double* pos) {
    double U = calc_U(pos);
    double Ux, Uy, Uz, Uxx, Uyy, Uzz, Uxy, Uxz, Uyz;
    calc_dU(pos, Ux, Uy, Uz);
    calc_ddU(pos, Uxx, Uyy, Uzz, Uxy, Uxz, Uyz);

    double t1 = Ux;      //diff(U(x,y,z),x);
    double t2 = t1 * t1;
    double t4 = Uxx;     // diff(diff(U(x,y,z),x),x);
    double t5 = U;       // U(x,y,z);
    double t6 = t4 * t5;
    double t8 = t5 * t5;
    double t9 = t8 * t8;
    double t11 = 1 / t9 / t8;
    double t13 = Uy;     // diff(U(x,y,z),y);
    double t14 = t1 * t13;
    double t16 = Uxy;    // diff(diff(U(x,y,z),x),y);
    double t17 = t16 * t5;
    double t19 = (-5.0 * t14 + t17) * t11;
    double t20 = Uz;     // diff(U(x,y,z),z);
    double t21 = t1 * t20;
    double t23 = Uxz;    // diff(diff(U(x,y,z),x),z);
    double t24 = t23 * t5;
    double t26 = (-5.0 * t21 + t24) * t11;
    double t27 = t13 * t13;
    double t29 = Uyy;    // diff(diff(U(x,y,z),y),y);
    double t30 = t29 * t5;
    double t33 = t13 * t20;
    double t35 = Uyz;    // diff(diff(U(x,y,z),y),z);
    double t36 = t35 * t5;
    double t38 = (-5.0 * t33 + t36) * t11;
    double t39 = t20 * t20;
    double t41 = Uzz;    // diff(diff(U(x,y,z),z),z);
    double t42 = t41 * t5;
    double t46 = 1 / t8;
    double t47 = (-t2 + t6) * t46;
    double t49 = (-t14 + t17) * t46;
    double t51 = (-t21 + t24) * t46;
    double t53 = (-t27 + t30) * t46;
    double t55 = (-t33 + t36) * t46;
    double t57 = (-t39 + t42) * t46;

    chrisD[0][0][0][0] = 0.0;
    chrisD[0][0][0][1] = 0.0;
    chrisD[0][0][0][2] = 0.0;
    chrisD[0][0][0][3] = 0.0;
    chrisD[0][0][1][0] = 0.0;
    chrisD[0][0][1][1] = -(-5.0 * t2 + t6) * t11;
    chrisD[0][0][1][2] = -t19;
    chrisD[0][0][1][3] = -t26;
    chrisD[0][0][2][0] = 0.0;
    chrisD[0][0][2][1] = -t19;
    chrisD[0][0][2][2] = -(-5.0 * t27 + t30) * t11;
    chrisD[0][0][2][3] = -t38;
    chrisD[0][0][3][0] = 0.0;
    chrisD[0][0][3][1] = -t26;
    chrisD[0][0][3][2] = -t38;
    chrisD[0][0][3][3] = -(-5.0 * t39 + t42) * t11;
    chrisD[0][1][0][0] = 0.0;
    chrisD[0][1][0][1] = -t47;
    chrisD[0][1][0][2] = -t49;
    chrisD[0][1][0][3] = -t51;
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
    chrisD[0][2][0][1] = -t49;
    chrisD[0][2][0][2] = -t53;
    chrisD[0][2][0][3] = -t55;
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
    chrisD[0][3][0][1] = -t51;
    chrisD[0][3][0][2] = -t55;
    chrisD[0][3][0][3] = -t57;
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
    chrisD[1][0][0][1] = -t47;
    chrisD[1][0][0][2] = -t49;
    chrisD[1][0][0][3] = -t51;
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
    chrisD[1][1][1][1] = t47;
    chrisD[1][1][1][2] = t49;
    chrisD[1][1][1][3] = t51;
    chrisD[1][1][2][0] = 0.0;
    chrisD[1][1][2][1] = -t49;
    chrisD[1][1][2][2] = -t53;
    chrisD[1][1][2][3] = -t55;
    chrisD[1][1][3][0] = 0.0;
    chrisD[1][1][3][1] = -t51;
    chrisD[1][1][3][2] = -t55;
    chrisD[1][1][3][3] = -t57;
    chrisD[1][2][0][0] = 0.0;
    chrisD[1][2][0][1] = 0.0;
    chrisD[1][2][0][2] = 0.0;
    chrisD[1][2][0][3] = 0.0;
    chrisD[1][2][1][0] = 0.0;
    chrisD[1][2][1][1] = t49;
    chrisD[1][2][1][2] = t53;
    chrisD[1][2][1][3] = t55;
    chrisD[1][2][2][0] = 0.0;
    chrisD[1][2][2][1] = t47;
    chrisD[1][2][2][2] = t49;
    chrisD[1][2][2][3] = t51;
    chrisD[1][2][3][0] = 0.0;
    chrisD[1][2][3][1] = 0.0;
    chrisD[1][2][3][2] = 0.0;
    chrisD[1][2][3][3] = 0.0;
    chrisD[1][3][0][0] = 0.0;
    chrisD[1][3][0][1] = 0.0;
    chrisD[1][3][0][2] = 0.0;
    chrisD[1][3][0][3] = 0.0;
    chrisD[1][3][1][0] = 0.0;
    chrisD[1][3][1][1] = t51;
    chrisD[1][3][1][2] = t55;
    chrisD[1][3][1][3] = t57;
    chrisD[1][3][2][0] = 0.0;
    chrisD[1][3][2][1] = 0.0;
    chrisD[1][3][2][2] = 0.0;
    chrisD[1][3][2][3] = 0.0;
    chrisD[1][3][3][0] = 0.0;
    chrisD[1][3][3][1] = t47;
    chrisD[1][3][3][2] = t49;
    chrisD[1][3][3][3] = t51;
    chrisD[2][0][0][0] = 0.0;
    chrisD[2][0][0][1] = -t49;
    chrisD[2][0][0][2] = -t53;
    chrisD[2][0][0][3] = -t55;
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
    chrisD[2][1][1][1] = t49;
    chrisD[2][1][1][2] = t53;
    chrisD[2][1][1][3] = t55;
    chrisD[2][1][2][0] = 0.0;
    chrisD[2][1][2][1] = t47;
    chrisD[2][1][2][2] = t49;
    chrisD[2][1][2][3] = t51;
    chrisD[2][1][3][0] = 0.0;
    chrisD[2][1][3][1] = 0.0;
    chrisD[2][1][3][2] = 0.0;
    chrisD[2][1][3][3] = 0.0;
    chrisD[2][2][0][0] = 0.0;
    chrisD[2][2][0][1] = 0.0;
    chrisD[2][2][0][2] = 0.0;
    chrisD[2][2][0][3] = 0.0;
    chrisD[2][2][1][0] = 0.0;
    chrisD[2][2][1][1] = -t47;
    chrisD[2][2][1][2] = -t49;
    chrisD[2][2][1][3] = -t51;
    chrisD[2][2][2][0] = 0.0;
    chrisD[2][2][2][1] = t49;
    chrisD[2][2][2][2] = t53;
    chrisD[2][2][2][3] = t55;
    chrisD[2][2][3][0] = 0.0;
    chrisD[2][2][3][1] = -t51;
    chrisD[2][2][3][2] = -t55;
    chrisD[2][2][3][3] = -t57;
    chrisD[2][3][0][0] = 0.0;
    chrisD[2][3][0][1] = 0.0;
    chrisD[2][3][0][2] = 0.0;
    chrisD[2][3][0][3] = 0.0;
    chrisD[2][3][1][0] = 0.0;
    chrisD[2][3][1][1] = 0.0;
    chrisD[2][3][1][2] = 0.0;
    chrisD[2][3][1][3] = 0.0;
    chrisD[2][3][2][0] = 0.0;
    chrisD[2][3][2][1] = t51;
    chrisD[2][3][2][2] = t55;
    chrisD[2][3][2][3] = t57;
    chrisD[2][3][3][0] = 0.0;
    chrisD[2][3][3][1] = t49;
    chrisD[2][3][3][2] = t53;
    chrisD[2][3][3][3] = t55;
    chrisD[3][0][0][0] = 0.0;
    chrisD[3][0][0][1] = -t51;
    chrisD[3][0][0][2] = -t55;
    chrisD[3][0][0][3] = -t57;
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
    chrisD[3][1][1][1] = t51;
    chrisD[3][1][1][2] = t55;
    chrisD[3][1][1][3] = t57;
    chrisD[3][1][2][0] = 0.0;
    chrisD[3][1][2][1] = 0.0;
    chrisD[3][1][2][2] = 0.0;
    chrisD[3][1][2][3] = 0.0;
    chrisD[3][1][3][0] = 0.0;
    chrisD[3][1][3][1] = t47;
    chrisD[3][1][3][2] = t49;
    chrisD[3][1][3][3] = t51;
    chrisD[3][2][0][0] = 0.0;
    chrisD[3][2][0][1] = 0.0;
    chrisD[3][2][0][2] = 0.0;
    chrisD[3][2][0][3] = 0.0;
    chrisD[3][2][1][0] = 0.0;
    chrisD[3][2][1][1] = 0.0;
    chrisD[3][2][1][2] = 0.0;
    chrisD[3][2][1][3] = 0.0;
    chrisD[3][2][2][0] = 0.0;
    chrisD[3][2][2][1] = t51;
    chrisD[3][2][2][2] = t55;
    chrisD[3][2][2][3] = t57;
    chrisD[3][2][3][0] = 0.0;
    chrisD[3][2][3][1] = t49;
    chrisD[3][2][3][2] = t53;
    chrisD[3][2][3][3] = t55;
    chrisD[3][3][0][0] = 0.0;
    chrisD[3][3][0][1] = 0.0;
    chrisD[3][3][0][2] = 0.0;
    chrisD[3][3][0][3] = 0.0;
    chrisD[3][3][1][0] = 0.0;
    chrisD[3][3][1][1] = -t47;
    chrisD[3][3][1][2] = -t49;
    chrisD[3][3][1][3] = -t51;
    chrisD[3][3][2][0] = 0.0;
    chrisD[3][3][2][1] = -t49;
    chrisD[3][3][2][2] = -t53;
    chrisD[3][3][2][3] = -t55;
    chrisD[3][3][3][0] = 0.0;
    chrisD[3][3][3][1] = t51;
    chrisD[3][3][3][2] = t55;
    chrisD[3][3][3][3] = t57;

    return true;
}


/*! Calculate Riemann tensor.
 *  \todo use symmetries
 *  \param pos : pointer to position.
 */
bool MetricExtremeReissnerNordstromDihole::calculateRiemann(const double* pos) {
    double U = calc_U(pos);
    double Ux, Uy, Uz, Uxx, Uyy, Uzz, Uxy, Uxz, Uyz;
    calc_dU(pos, Ux, Uy, Uz);
    calc_ddU(pos, Uxx, Uyy, Uzz, Uxy, Uxz, Uyz);

    double U2  = U*U;
    double U6  = U2*U2*U2;
    double Ux2 = Ux * Ux;
    double Uy2 = Uy * Uy;
    double Uz2 = Uz * Uz;

    int t = 0;
    int x = 1;
    int y = 2;
    int z = 3;

    riem[t][x][x][t] = -1/U2*(-3*Ux2+Uxx*U+Uy2+Uz2);
    riem[t][x][y][t] = -1/U2*(-4*Ux*Uy+Uxy*U);
    riem[t][x][z][t] = -1/U2*(-4*Ux*Uz+Uxz*U);
    riem[t][y][x][t] = -1/U2*(-4*Ux*Uy+Uxy*U);
    riem[t][y][y][t] = 1/U2*(3*Uy2-Uyy*U-Ux2-Uz2);
    riem[t][y][z][t] = 1/U2*(4*Uy*Uz-Uyz*U);
    riem[t][z][x][t] = -1/U2*(-4*Ux*Uz+Uxz*U);
    riem[t][z][y][t] = 1/U2*(4*Uy*Uz-Uyz*U);
    riem[t][z][z][t] = -1/U2*(-3*Uz2+Uzz*U+Ux2+Uy2);
    riem[x][y][x][y] = 1/U2*(Uy2-Uyy*U+Ux2-Uxx*U-Uz2);
    riem[x][y][x][z] = 1/U2*(-Uyz*U+2*Uy*Uz);
    riem[x][y][y][z] = 1/U2*(-2*Ux*Uz+Uxz*U);
    riem[x][z][x][y] = 1/U2*(-Uyz*U+2*Uy*Uz);
    riem[x][z][x][z] = -(-Uz2+Uzz*U-Ux2+Uxx*U+Uy2)/U2;
    riem[x][z][y][z] = -(Uxy*U-2*Ux*Uy)/U2;
    riem[x][t][x][t] = -1/U6*(-3*Ux2+Uxx*U+Uy2+Uz2);
    riem[x][t][y][t] = -1/U6*(-4*Ux*Uy+Uxy*U);
    riem[x][t][z][t] = -1/U6*(-4*Ux*Uz+Uxz*U);
    riem[y][x][x][y] = (-Uy2+Uyy*U-Ux2+Uxx*U+Uz2)/U2;
    riem[y][x][x][z] = (Uyz*U-2*Uy*Uz)/U2;
    riem[y][x][y][z] = -1/U2*(-2*Ux*Uz+Uxz*U);
    riem[y][z][x][y] = 1/U2*(-2*Ux*Uz+Uxz*U);
    riem[y][z][x][z] = -(Uxy*U-2*Ux*Uy)/U2;
    riem[y][z][y][z] = -(-Uz2+Uzz*U-Uy2+Uyy*U+Ux2)/U2;
    riem[y][t][x][t] = -1/U6*(-4*Ux*Uy+Uxy*U);
    riem[y][t][y][t] = 1/U6*(3*Uy2-Uyy*U-Ux2-Uz2);
    riem[y][t][z][t] = 1/U6*(4*Uy*Uz-Uyz*U);
    riem[z][x][x][y] = (Uyz*U-2*Uy*Uz)/U2;
    riem[z][x][x][z] = (-Uz2+Uzz*U-Ux2+Uxx*U+Uy2)/U2;
    riem[z][x][y][z] = (Uxy*U-2*Ux*Uy)/U2;
    riem[z][y][x][y] = -1/U2*(-2*Ux*Uz+Uxz*U);
    riem[z][y][x][z] = (Uxy*U-2*Ux*Uy)/U2;
    riem[z][y][y][z] = (-Uz2+Uzz*U-Uy2+Uyy*U+Ux2)/U2;
    riem[z][t][x][t] = -1/U6*(-4*Ux*Uz+Uxz*U);
    riem[z][t][y][t] = 1/U6*(4*Uy*Uz-Uyz*U);
    riem[z][t][z][t] = -1/U6*(-3*Uz2+Uzz*U+Ux2+Uy2);

    return true;
}

/*! Transform local 4-direction to coordinate 4-direction.
 *
 *  \param  pos  :  pointer to position array.
 *  \param  ldir :  pointer to local direction array.
 *  \param  dir  :  pointer to calculated coordinate direction array.
 *  \param  type :  type of tetrad.
 */
void MetricExtremeReissnerNordstromDihole::localToCoord(const double* pos, const double* ldir, double* dir,
        enum_nat_tetrad_type) {
    double U = calc_U(pos);

    dir[0] = ldir[0] * U;
    dir[1] = ldir[1] / U;
    dir[2] = ldir[2] / U;
    dir[3] = ldir[3] / U;
}

/*! Transform coordinate 4-direction to local 4-direction.
 *
 *  \param  pos  :  pointer to position array.
 *  \param  cdir :  pointer to coordinate direction.
 *  \param  ldir :  pointer to calculated local direction array.
 *  \param  type :  type of tetrad.
 */
void MetricExtremeReissnerNordstromDihole::coordToLocal(const double* pos, const double* cdir, double* ldir,
        enum_nat_tetrad_type) {
    double U = calc_U(pos);

    ldir[0] = cdir[0] / U;
    ldir[1] = cdir[1] * U;
    ldir[2] = cdir[2] * U;
    ldir[3] = cdir[3] * U;
}


/*! Test break condition.
 *
 *  \param pos    : pointer to position array.
 *  \return true  : radial position r < 0.0 or  r^2<=(1.0+eps)*rs^2.
 *  \return false : position is valid.
 */
bool MetricExtremeReissnerNordstromDihole::breakCondition(const double*) {
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
double MetricExtremeReissnerNordstromDihole::testConstraint(const double y[], const double kappa) {
    double U = calc_U(y);

    double dt = y[4];
    double dx = y[5];
    double dy = y[6];
    double dz = y[7];

    double sum = -kappa;
    sum += -dt * dt / (U * U) + U * U * (dx * dx + dy * dy + dz * dz);
    return sum;
}


/*! Set parameter 'pName' to 'val'.
 *
 *  Set 'mass' or 'lambda' parameter.
 */
bool MetricExtremeReissnerNordstromDihole::setParam(const char* pName, double val) {
    Metric::setParam(pName, val);

    if (pName == "m1") {
        mM1 = val;
    } else if (pName == "m2") {
        mM2 = val;
    }
    return true;
}

/*! Transform point p to 2+1 coordinates.
 *
 *  \param  p  : point in proper metric coordinates.
 *  \param  cp : reference to transformed point.
 *  \return true : success.
 */
bool MetricExtremeReissnerNordstromDihole::transToTwoPlusOne(vec4 p, vec4 &cp) {
    double x = p[1];
    double y = p[2];
    double z = p[3];
    double rho = sqrt(x * x + y * y);

    cp = vec4(p[0], rho, z, p[0]);
    return true;
}

/*! Generate report.
 */
bool MetricExtremeReissnerNordstromDihole::report(const vec4 pos, const vec4 cdir, std::string &text) {
    std::stringstream ss;
    ss << "Report for extreme ReissnerNordstrom Dihole metric\n\tcoordinate : (t,x,y,z)\n";
    ss << "The two singularities are fixed along the z-axis at z=-1 and z=+1. \n";
    ss << "---------------------------------------------------------------------------------\n";
    ss << "  physical units ................................. no \n \n";
    ss.precision(DEF_FIXED_REPORT_PRECISION);
    ss.setf(std::ios::fixed);

    double rho_init_square = pos[1] * pos[1] + pos[2] * pos[2];
    double rho_init = sqrt(rho_init_square);
    double U = 1 + mM1 / sqrt(rho_init_square + (pos[3] - 1) * (pos[3] - 1.)) + mM2 / sqrt(rho_init_square + (pos[3] + 1) * (pos[3] + 1.));
    double k = (mSpeedOfLight * mSpeedOfLight * cdir[0]) / (U * U);
    double Lz = U * U * (pos[1] * cdir[2] - pos[2] * cdir[1]);
    ss << "  function from metric at initial psition .......... U = " << U << std::endl;
    ss << "  constant of motion (energy) ...................... k = " << k << std::endl;
    ss << "  constant of motion (angular momentum) ........... Lz = " << Lz << std::endl;
    ss << "  point of equilibrium ............................. z = ";

    double pointOfEquilibrium;
    calcPointOfEquilibrium(pointOfEquilibrium);
    ss << pointOfEquilibrium << std::endl;

    ss << "\n\n";
    ss << "for additional information about motion along the z-axis make sure that\n";
    ss << "x=y=0 and chi=0 or chi=180 in the standard tetrad direction ! \n";

    if (pos[1] == 0.0 && pos[2] == 0.0 && cdir[1] < 0.00000000000001 && cdir[1] > -0.00000000000001 && cdir[2] < 0.00000000000001 && cdir[2] > -0.00000000000001) {
        if (pos[3] > 1.0) {
            double betaEscapeAlongZ;
            calcEscapeVelocityAlongZ(pos, betaEscapeAlongZ, mM1, mM2);
            ss << "  \n";
            ss << "  escape velocity upward (chi=0) ................ beta = " << betaEscapeAlongZ << std::endl;
        }
        if (pos[3] < -1.0) {
            double betaEscapeAlongZ;
            calcEscapeVelocityAlongZ(pos, betaEscapeAlongZ , mM1, mM2);
            ss << "  \n";
            ss << "  escape velocity downward (chi=180) ............ beta = " << betaEscapeAlongZ << std::endl;
        }
        if (pos[3] == -1.0) {
            ss << "  \n";
            ss << "  initial value is the lower singularity !!!!! \n";
        }
        if (pos[3] == 1.0) {
            ss << "  \n";
            ss << "  initial value is the upper singularity !!!!! \n";
        }
        if (pos[3] > -1.0 && pos[3] < 1.0) {
            double betaEquilibrium = 0.0;

            if (pos[3] == pointOfEquilibrium) {
                betaEquilibrium = 0.0;
                ss << "  \n";
                ss << "  vel. for asymptotic motion to point of equilibrium ... beta = " << betaEquilibrium << std::endl;
            }
            if (pos[3] < pointOfEquilibrium) {
                calcVelocityToEquilibriumPoint(pos, betaEquilibrium, pointOfEquilibrium, mM1, mM2);
                ss << "  \n";
                ss << "  vel. for asymptotic motion to point of equilibrium ... beta = " << betaEquilibrium << std::endl;
                ss << "  appendant orientation ................................. chi = 0.00000000000000 \n";
            }
            if (pos[3] > pointOfEquilibrium) {
                calcVelocityToEquilibriumPoint(pos, betaEquilibrium, pointOfEquilibrium, mM1, mM2);
                ss << "  \n";
                ss << "  vel. for asymptotic motion to point of equilibrium ... beta = " << betaEquilibrium << std::endl;
                ss << "  appendant orientation ................................. chi = 180.000000000000 \n";
            }
        }

        ss << "  lightlike motion is free, but can end in one of the two singularities \n";
    }

    ss << "\n\n";
    ss << "for additional information about motion along the x-axis make sure that \n";
    ss << "m1=m2, y=z=0, chi=90 and ksi=0 or ksi=180 in the standard tetrad direction ! \n";

    if (mM1 == mM2 && pos[2] == 0.0 && pos[3] == 0.0 && cdir[2] < 0.00000000000001 && cdir[2] > -0.00000000000001 && cdir[3] < 0.00000000000001 && cdir[3] > -0.00000000000001) {
        double betaEscapeAlongX;
        calcEscapeVelocityAlongX(pos, betaEscapeAlongX, mM1);
        ss << "  \n";
        ss << "  escape velocity ................................... beta = " << betaEscapeAlongX << std::endl;

        double periodicAlongX;
        calcPeriodicAlongX(pos, periodicAlongX, mM1);
        ss << "  periodic through the origin at rest (proper time) .... T = " << periodicAlongX << std::endl;

        double betaAtOriginAlongX;
        calcVelocityAtOriginAlongX(pos, betaAtOriginAlongX, mM1);
        ss << "  appendant velocity at the origin .................. beta = " << betaAtOriginAlongX << std::endl;
        ss << "  lightlike motion is free \n";
    }

    ss << "\n\n";
    ss << "for circular orbits in the x/y-plane make sure that m1=m2=:m, z=0 \n";
    ss << "and chi=90 in the standard tetrad direction ! \n";

    if (mM1 == mM2 && pos[3] == 0.0 && cdir[3] < 0.00000000000001 && cdir[3] > -0.00000000000001) {
        double mquer = 1.8371173070873836; // sqrt(27./8.)
        double mstern = 1.0294905868123894; // ((13.+sqrt(129.))*sqrt(710.+70.*sqrt(129.)))/(50.*(7.+sqrt(129.)))

        ss << "  \n";
        ss << "  no Photon orbits, only stable timelike orbits .......... m < " << mstern << std::endl;
        ss << "                                            and .......... m = " << mstern << std::endl;
        ss << "  no Photon orbits, one last unstable timelike orbit ..... m < " << mquer << std::endl;
        ss << "  one Photon orbit, one last unstable timelike orbit ..... m = " << mquer << std::endl;
        ss << "  two Photon orbits, one last unstable timelike orbit .... m > " << mquer << std::endl;

        double r_Nm, r_Np;
        calcPhotOrbitsXY(r_Nm, r_Np, mquer);

        ss << "  \n";
        ss << "  inner, stable Photon orbit .......................... r_nm = " << r_Nm << std::endl;
        ss << "  outer, unstable Photon orbit ........................ r_np = " << r_Np << std::endl;

        ss << "  \n";
        ss << "  innermost stable circular orbit (point of equilibrium) ... r_equ = 0.0" << std::endl;
        ss << "  appendant velocity ........................................ beta = 0.0" << std::endl;

        if (mM1 > mstern && mM1 < mquer) {
            double radiusIucoXY, betaCoXYIuco;
            calcRadiusIucoXY(radiusIucoXY, mM1);
            calcBetaCoXY(betaCoXYIuco, radiusIucoXY);

            ss << "  innermost unstable circular orbit ................. r_iuco = " << radiusIucoXY << std::endl;
            ss << "  appendant velocity .................................. beta = " << betaCoXYIuco << std::endl;
        }
        if (mM1 > mstern) {
            double radiusIscoXY, betaCoXYIsco;
            calcRadiusIscoXY(radiusIscoXY, mM1);
            calcBetaCoXY(betaCoXYIsco, radiusIscoXY);

            ss << "  second innermost stable circular orbit ............ r_isco = " << radiusIscoXY << std::endl;
            ss << "  appendant velocity .................................. beta = " << betaCoXYIsco << std::endl;
        }

        ss << "  \n";
        ss << "  critical angle for null geodesic ................. ksiCrit = ";

        double ksicritXY;
        double radius_pos = sqrt(pos[1] * pos[1] + pos[2] * pos[2]);
        if (radius_pos >= r_Nm && mM1 >= mquer) {
            calcKsiCritXY(pos, ksicritXY, r_Np, mM1);
            ss << ksicritXY << std::endl;
        } else {
            ss << "not valid here \n";
        }

        ss << "  velocity for circular orbit at initial position ..... beta = ";

        double betaCoXY;
        calcBetaCoXY(betaCoXY, rho_init);
        if (betaCoXY <= 1.) {
            ss << betaCoXY << std::endl;
        } else {
            ss << "not valid here \n";
        }


        ss << "  \n";
        ss << "  (The terms stable and unstable refer here only to the motion in the x/y-plane.\n";
        ss << "   Globally seen, all of these orbits are unstable.) \n";
    }

    ss << "\n\n";
    ss << "for circular orbits outside the x/y-plane for equal masses make sure that\n";
    ss << "m1=m2:=m, chi=90 in the standard tetrad direction and z unequal 0 ! \n";

    if (mM1 == mM2 && cdir[3] < 0.00000000000001 && cdir[3] > -0.00000000000001 && pos[3] != 0.0) {
        double zPhot, rPhot;
        calcHeightOfOtherPhotonOrbits(zPhot, mM1);
        calcRadiusForOtherCo(rPhot, zPhot);
        ss << "  \n";
        ss << "  Circular orbits exist for ..................... m <  " << 2.598076211;
        ss << "  \n \n";

        if (mM1 < 2.598076211) {
            ss << "  Photon orbit for .............................. z =  " << zPhot << std::endl;
            ss << "               and .............................. z = " << -zPhot << std::endl;
            ss << "       with radius .............................. r =  " << rPhot << std::endl;
            ss << "              ( calculated numerically, works for m >  0.00005 ) \n";

            double rTime, betaTime;
            calcRadiusForOtherCo(rTime, pos[3]);
            calcVelocityForOtherTimelikeCo(pos, betaTime, mM1);
            ss << "  \n";
            ss << "  Timelike orbit at current z with radius ....... r = " << rTime << std::endl;
            ss << "                             and velocity .... beta = " << betaTime << std::endl;
            ss << "           orbit with the same r and beta ....... z = " << -pos[3] << std::endl;
            ss << "                   ( Attention: beta has to be between 0.0 and 1.0 ! ) \n";
        }
    }

    ss << "\n\n";
    ss << "for circular orbits outside the x/y-plane for unequal masses make sure that \n";
    ss << "m1 unequal m2, chi=90 in the standard tetrad direction and z unequal 0 ! \n";

    if (mM1 != mM2 && cdir[3] < 0.00000000000001 && cdir[3] > -0.00000000000001 && pos[3] != 0.0) {
        double zEqu, zSing;
        double minusEins = -1.;
        double plusEins = 1.;
        calcPointOfEquilibrium(zEqu);
        calcASingularityPoint(zSing);

        if (mM1 < mM2) {
            ss << "  \n";
            ss << "  Circular orbits exist only between ............ z = "  << minusEins << std::endl;
            ss << "                                 and ............ z = "  << zSing << std::endl;
            ss << "                         and between ............ z =  " << zEqu << std::endl;
            ss << "                                 and ............ z =  " << plusEins << std::endl;
        }
        if (mM1 > mM2) {
            ss << "  \n";
            ss << "  Circular orbits exist only between ............ z = "  << minusEins << std::endl;
            ss << "                                 and ............ z = "  << zEqu << std::endl;
            ss << "                         and between ............ z =  " << zSing << std::endl;
            ss << "                                 and ............ z =  " << plusEins << std::endl;
        }

        double zPhot1, zPhot2, rPhot1, rPhot2;
        calcHeightOfOtherPhotonOrbitsForUnequalMasses(zPhot1, zPhot2, mM1, mM2);
        calcRadiusForOtherCoForUnequalMasses(rPhot1, zPhot1, mM1, mM2);
        calcRadiusForOtherCoForUnequalMasses(rPhot2, zPhot2, mM1, mM2);
        ss << "  \n";
        ss << "  Photon orbit for .............................. z = "  << zPhot1 << std::endl;
        ss << "       with radius .............................. r =  " << rPhot1 << std::endl;
        ss << "  Photon orbit for .............................. z =  " << zPhot2 << std::endl;
        ss << "       with radius .............................. r =  " << rPhot2 << std::endl;

        if (mM1 < mM2) {
            if ((pos[3] > -1. && pos[3] < zSing) || (pos[3] > zEqu && pos[3] < 1.)) {
                double rTime1, betaTime1;
                calcRadiusForOtherCoForUnequalMasses(rTime1, pos[3], mM1, mM2);
                calcVelocityForOtherTimelikeCoForUnequalMasses(pos, betaTime1, mM1, mM2);

                ss << "  \n";
                ss << "  Timelike orbit at current z with radius ....... r =  " << rTime1 << std::endl;
                ss << "                             and velocity .... beta =  " << betaTime1 << std::endl;
                ss << "                   ( Attention: beta has to be between 0.0 and 1.0 ! ) \n";
            } else {
                ss << "  \n";
                ss << "  Timelike orbit at current z with radius ....... r =  not valid here \n";
                ss << "                             and velocity .... beta =  not valid here \n";
            }
        }
        if (mM1 > mM2) {
            if ((pos[3] > -1. && pos[3] < zEqu) || (pos[3] > zSing && pos[3] < 1.)) {
                double rTime1, betaTime1;
                calcRadiusForOtherCoForUnequalMasses(rTime1, pos[3], mM1, mM2);
                calcVelocityForOtherTimelikeCoForUnequalMasses(pos, betaTime1, mM1, mM2);

                ss << "  \n";
                ss << "  Timelike orbit at current z with radius ....... r =  " << rTime1 << std::endl;
                ss << "                             and velocity .... beta =  " << betaTime1 << std::endl;
                ss << "                   ( Attention: beta has to be between 0.0 and 1.0 ! ) \n";
            } else {
                ss << "  \n";
                ss << "  Timelike orbit at current z with radius ....... r =  not valid here \n";
                ss << "                             and velocity .... beta =  not valid here \n";
            }
        }
    }

    text = ss.str();
    return true;
}


// ****************************** specific public methods **************************
/*!  Calculate the escape velocity for motion only along the z-axis for initial
 *   position outside the two singularities
 *
 *  \param pos : current position
 *  \param betaEscape : reference to the escape velocity
 *  \param m1 : value of the mass at z=+1
 *  \param m2 : value of the mass at z=-1
 */
bool MetricExtremeReissnerNordstromDihole::calcEscapeVelocityAlongZ(const vec4 pos, double &betaEscape, double m1, double m2) {
    double z = pos[3];
    double w1, w2;

    if (z > 1.0) {
        w1 = z + 1.;
        w2 = z - 1.;
    }
    if (z < -1.0) {
        w1 = -z - 1.;
        w2 = 1. - z;
    }

    double w3 = w1 * w1;
    double w4 = w2 * w2;
    double w5 = w1 * w2;
    double w6 = w5 * w5;

    double numerator = m1 * m1 * w3 + m2 * m2 * w4 + 2.*m1 * w2 * w3 + 2.*m2 * w1 * w4 + 2.*m1 * m2 * w5;
    double denominator = w6 + m1 * m1 * w3 + m2 * m2 * w4 + 2.*m1 * w2 * w3 + 2.*m2 * w1 * w4 + 2.*m1 * m2 * w5;

    betaEscape = sqrt(numerator / denominator);

    return true;
}

/*!  Calculate the escape velocity for motion only along the x-axis for equal masses m1=m2
 *
 *  \param pos : current position
 *  \param betaEscape : reference to the escape velocity
 *  \param MASS : value of the two equal masses
 */
bool MetricExtremeReissnerNordstromDihole::calcEscapeVelocityAlongX(const vec4 pos, double &betaEscape, double MASS) {
    double xSquare = pos[1] * pos[1];
    double MSquare = MASS * MASS;
    double w1 = xSquare + 1.;
    double w2 = sqrt(w1);

    double numerator = 4.*MSquare + 4.*MASS * w2;
    double denominator = w1 + numerator;

    betaEscape = sqrt(numerator / denominator);

    return true;
}

/*!  Calculate the point of equilibrium along the z-axis
 *
 *   For this z-value a timelike particle with no initial velocity remains static.
 *
 *  \param pointOfEquilibrium : reference to the point of equilibrium
 */
bool MetricExtremeReissnerNordstromDihole::calcPointOfEquilibrium(double &pointOfEquilibrium) {
    if (mM1 == mM2) {
        pointOfEquilibrium = 0.0;
        return true;
    } else {
        double quot = mM1 / mM2;
        pointOfEquilibrium = (1. + quot - 2.*sqrt(quot)) / (1. - quot);
        return true;
    }
}

/*!  Calculate a singularity along the z-axis needed for some calculations
 *
 *   At this z-value there is a singularity in the rho(z)-function for circular orbits in
 *   planes normal to the z-axis for unequal masses. This function is needed for some
 *   numerical calculations and some case-differentiations.
 *
 *  \param singularityPoint : reference to the singularity
 */
bool MetricExtremeReissnerNordstromDihole::calcASingularityPoint(double &singularityPoint) {
    double q = mM1 / mM2;
    singularityPoint = (q - 1.) / (q + 1.);

    return true;
}

/*!  Calculate the necessary velocity to asymptotically run to the point of equilibrium
 *   between the two singularities along the z-axis
 *
 *   This value only make sense for -1 < z < +1
 *
 *  \param pos : current position
 *  \param betaEquilibrium : reference to the necessary velocity
 *  \param pointOfEquilibrium : location of the point of equilibrium
 *  \param mass1 : value of the mass at z = +1
 *  \param mass2 : value of the mass at z = -1
 */
bool MetricExtremeReissnerNordstromDihole::calcVelocityToEquilibriumPoint(const vec4 pos, double &betaEquilibrium, double pointOfEquilibrium, double mass1, double mass2) {
    double z0 = pos[3];
    double zE = pointOfEquilibrium;
    double Uz0 = 1. + mass1 / sqrt((z0 - 1.) * (z0 - 1.)) + mass2 / sqrt((z0 + 1.) * (z0 + 1.));
    double UzE = 1. + mass1 / sqrt((zE - 1.) * (zE - 1.)) + mass2 / sqrt((zE + 1.) * (zE + 1.));

    betaEquilibrium = sqrt(1. - UzE * UzE / (Uz0 * Uz0));

    return true;
}

/*!  Calculate the periodic time (proper time) for the motion through the center of mass and
 *   back along the x-axis for equal masses m1 and m2 in the x/y-plane for vanishing initial
 *   velocity
 *
 *   For equal masses in the x/y-plane the effective potential has a global minimum at the
 *   center of mass (x=y=z=0). So there is the possibility to make a periodic motion from the
 *   one side of the center of mass to the other through the center of mass and back again.
 *   The initial velocity is set to null. The calculation is done by a numerical integration.
 *
 *  \param pos : current position
 *  \param periodicAlongX : reference to the periodic time
 *  \param MASS : value of the two equal massea m1 and m2
 */
bool MetricExtremeReissnerNordstromDihole::calcPeriodicAlongX(const vec4 pos, double &periodicAlongX, double MASS) {
    double x0 = pos[1];

    if (x0 == 0.0) {
        periodicAlongX = 0.0;
    } else {
        double temp;
        double offset = 0.00000001;

        struct integrandForPeriodicAlongX_params par = {MASS, x0};
        F.params = &par;

        size_t limit = 1000;
        int    key   = GSL_INTEG_GAUSS21;
        double error;

        gsl_set_error_handler_off();
        gsl_integration_qag(&F, -x0 + offset, x0 - offset, 0.0, 1e-6, limit, key, w, &temp, &error);

        periodicAlongX = sqrt(temp * temp);
    }

    return true;
}

/*!  Calculate the velocity at origin for the periodic motion described at the function
 *   calcPeriodicAlongX(...)
 *
 *  \param pos : current position
 *  \param betaAtOriginAlongX : reference to the searched velocity
 *  \param MASS : value of the two equal massea m1 and m2
 */
bool MetricExtremeReissnerNordstromDihole::calcVelocityAtOriginAlongX(const vec4 pos, double &betaAtOriginAlongX, double MASS) {
    double xSquare = pos[1] * pos[1];
    double w1 = xSquare + 1.;
    double w2 = 1. + 2.*MASS;
    double w3 = 2.*MASS + sqrt(w1);

    double numerator = w3 * w3;
    double denominator = w2 * w2 * w1;

    betaAtOriginAlongX = sqrt(1. - numerator / denominator);

    return true;
}

/*!  Calculate the Photon orbits for equal masses m1 und m2 for motion only in the z=0-plane
 *
 *   The number of Photon orbits depends on the mass-value. Null, one or two orbits
 *   are possible.
 *
 *  \param r_Nm : reference to the smaller orbit
 *  \param r_Np : reference to the larger orbit
 *  \param mquer : mass parameter, for which the two orbits have the same value
 */
bool MetricExtremeReissnerNordstromDihole::calcPhotOrbitsXY(double &r_Nm, double &r_Np, double mquer) {
    if (mM1 >= mquer) {
        double m = mM1 * mM1;
        double numerator = m * (27. - 36.*m + 8.*m * m);
        double denominator = 8.*sqrt(m * m * m * (m - 3.) * (m - 3.) * (m - 3.));
        double w1 = acos(numerator / denominator);
        double w2 = sqrt(m * (m - 3.));
        double w3 = (4. / 3.) * (m + 2.*w2 * cos(w1 / 3.));
        double w4 = (4. / 3.) * (m - 2.*w2 * cos((M_PI + w1) / 3.));

        r_Nm = sqrt(w4 - 1);
        r_Np = sqrt(w3 - 1);
    } else {
        r_Nm = 0.0;
        r_Np = 0.0;
    }
    return true;
}

/*!  Calculate the critical angle for null-geodesics in the x/y-plane for equal masses m1 and m2
 *
 *   For an observer at initial distance from the origin, it has an angular diameter,
 *   the critical angle. This angle additionally depends on the mass value and has a behaviour
 *   related to the critical angle in the ReissnerNordstrom-spacetime.
 *   For the existance of this critical angle the existance of Photon-orbits is necessary, so
 *   the mass has to be larger than sqrt(27/8).
 *
 *  \param pos : current position
 *  \param ksicritXY : reference to the critical angle
 *  \param r_Np : radius of the bigger one of the Photon-orbits
 *  \param MASS : value of the equal masses m1 and m2
 */
bool MetricExtremeReissnerNordstromDihole::calcKsiCritXY(const vec4 pos, double &ksicritXY, double r_Np, double MASS) {
    double sqrt_Np = sqrt(1. + r_Np * r_Np);
    double U_Np = 1. + 2.*MASS / sqrt_Np;

    double radius_pos = sqrt(pos[1] * pos[1] + pos[2] * pos[2]);
    double sqrt_pos = sqrt(1. + radius_pos * radius_pos);
    double U_pos = 1. + 2.*MASS / sqrt_pos;

    if (radius_pos >= r_Np) {
        ksicritXY = 180 - (180 / M_PI) * asin((r_Np * U_Np * U_Np) / (radius_pos * U_pos * U_pos));
    } else {
        ksicritXY = (180 / M_PI) * asin((r_Np * U_Np * U_Np) / (radius_pos * U_pos * U_pos));
    }

    return true;
}

/*!  Calculate the second innermost stable circular orbit in the z=0-plane for equal
 *   masses m1 and m2 with the numerical root-finding algorithm by Newton
 *
 *  \param radiusIscoXY : reference to the radius of the second innermost stable
 *                        circular orbit
 *  \param MASS : mass of the two equal masses m1 and m2
 */
bool MetricExtremeReissnerNordstromDihole::calcRadiusIscoXY(double &radiusIscoXY, double MASS) {
    double mstern = 1.0294905868123894;

    if (MASS == mstern) {
        double x = 15.050471684120383;
        radiusIscoXY  = sqrt(x - 1);
        return true;
    } else {
        //int status;
        const gsl_root_fdfsolver_type *T;
        gsl_root_fdfsolver *s;

        //double x0;
        double x;
        double w1 = sqrt(25.*MASS * MASS - 6.);
        double numerator = MASS * (250.*MASS * MASS - 387.);
        double denominator = 2.*w1 * w1 * w1;
        double w2 = acos(numerator / denominator);
        double w3 = cos(w2 / 3.);
        double w4 = (5.*MASS + 2.*w1 * w3);
        x = 100. + w4 * w4 / 9.;               // initial value

        gsl_function_fdf FDF;
        struct inflectionPoints_params params = {MASS};

        FDF.f = &m4d::inflectionPoints;
        FDF.df = &m4d::inflectionPoints_deriv;
        FDF.fdf = &m4d::inflectionPoints_fdf;
        FDF.params = &params;

        T = gsl_root_fdfsolver_newton;
        s = gsl_root_fdfsolver_alloc(T);
        gsl_root_fdfsolver_set(s, &FDF, x);

        for (int iter = 1; iter <= 30; iter++) {
            //status = gsl_root_fdfsolver_iterate (s);
            gsl_root_fdfsolver_iterate(s);
            //x0 = x;
            x = gsl_root_fdfsolver_root(s);
        }

        radiusIscoXY  = sqrt(x - 1);

        gsl_root_fdfsolver_free(s);
        return true;
    }
}

/*!  Calculate the innermost unstable circular orbit in the z=0-plane for equal
 *   masses m1 and m2 with the numerical root-finding algorithm by Newton
 *
 *  \param radiusIucoXY : reference to the radius of the innermost unstable
 *                        circular orbit
 *  \param MASS : mass of the two equal masses m1 and m2
 */
bool MetricExtremeReissnerNordstromDihole::calcRadiusIucoXY(double &radiusIucoXY, double MASS) {
    double mstern = 1.0294905868123894;

    if (MASS == mstern) {
        double x = 15.050471684120383;
        radiusIucoXY  = sqrt(x - 1);
        return true;
    } else {
        //int status;
        const gsl_root_fdfsolver_type *T;
        gsl_root_fdfsolver *s;

        //double x0;
        double x;
        double w1 = 16. / 3.;
        double w2 = sqrt(25.*MASS * MASS - w1);
        double numerator = MASS * (125.*MASS * MASS - 128.);
        double denominator = w2 * w2 * w2;
        double w3 = acos(numerator / denominator);
        double w4 = cos(w3 / 3.);
        double w5 = 1.25 * MASS + 0.5 * w2 * w4;
        x = w5 * w5;                          // initial value

        gsl_function_fdf FDF;
        struct inflectionPoints_params params = {MASS};

        FDF.f = &m4d::inflectionPoints;
        FDF.df = &m4d::inflectionPoints_deriv;
        FDF.fdf = &m4d::inflectionPoints_fdf;
        FDF.params = &params;

        T = gsl_root_fdfsolver_newton;
        s = gsl_root_fdfsolver_alloc(T);
        gsl_root_fdfsolver_set(s, &FDF, x);

        for (int iter = 1; iter <= 30; iter++) {
            //status = gsl_root_fdfsolver_iterate (s);
            gsl_root_fdfsolver_iterate(s);
            //x0 = x;
            x = gsl_root_fdfsolver_root(s);
        }

        radiusIucoXY  = sqrt(x - 1);

        gsl_root_fdfsolver_free(s);
        return true;
    }
}

/*!  Calculate the velocity on a circular orbit in the z=0-plane for equal masses m1 and m2
 *
 *  \param betaCoXY : reference to the velocity on the circular orbit
 *  \param radiusCoXY : radius of the circular orbit
 */
bool MetricExtremeReissnerNordstromDihole::calcBetaCoXY(double &betaCoXY, double radiusCoXY) {
    double w1 = radiusCoXY * radiusCoXY + 1.;
    double numerator = 2.*mM1 * (w1 - 1.);
    double denominator = 2.*mM1 + sqrt(w1) * sqrt(w1) * sqrt(w1);

    betaCoXY = sqrt(numerator / denominator);

    return true;
}

/*!  Calculate the appendant radius of a circular orbit in a plane normal to the z-axis for
 *   a z-value unequal to null and for equal masses m1=m2
 *
 *   One can show, that such orbits, timelike and lightlike, can only exist for -1 < z < +1.
 *   The relevant z-value (height of the orbit) has to be calculated numerically (behold
 *   also other class-methods).
 *
 *  \param radiusForOtherCo : reference to the radius of circular orbits for z unequal to zero
 *  \param z : height of the circular orbit
 */
bool MetricExtremeReissnerNordstromDihole::calcRadiusForOtherCo(double &radiusForOtherCo, double z) {
    double w1 = 1. + z;
    double w2 = 1. - z;
    double w3 = w1 * w1;
    double w4 = w2 * w2;
    double w5 = w4 / w3;

    double denominator = 1. - pow(w5, 1. / 3.);
    double w6 = 4.*z / denominator - w3;

    radiusForOtherCo = sqrt(w6);
    return true;
}

/*!  Calculate the appendant radius of a circular orbit in a plane normal to the z-axis for
 *   a z-value unequal to null and for unequal masses m1 and m2
 *
 *   One can show, that such orbits, timelike and lightlike, can only exist for -1 < z < +1.
 *   The relevant z-value (height of the orbit) has to be calculated numerically (behold
 *   also other class-methods).
 *
 *  \param radiusForOtherCo : reference to the radius of circular orbits for z unequal to zero
 *  \param z : height of the circular orbit
 *  \param M1 : value of the first mass at z=+1
 *  \param M2 : value of the second mass at z=-1
 */
bool MetricExtremeReissnerNordstromDihole::calcRadiusForOtherCoForUnequalMasses(double &radiusForOtherCo, double z, double M1, double M2) {
    double w1 = 1. + z;
    double w2 = 1. - z;
    double w3 = w1 * w1;
    double w4 = w2 * w2;
    double w5 = w4 / w3;

    double denominator = 1. - pow(M1 / M2, 2. / 3.) * pow(w5, 1. / 3.);
    double w6 = 4.*z / denominator - w3;

    radiusForOtherCo = sqrt(w6);
    return true;
}

/*!  Calculate the height of a Photon orbit in a plane normal to the z-axis for a z-value
 *   unequal to null and for equal masses m1=m2 numerically using the root-finding algorithm
 *   by Newton
 *
 *   The only parameter in this Newton-algorithm is the mass-value of the two equal masses.
 *   The necessary function and its derivative are written further down in this data.
 *
 *  \param heightOfOtherPhotonOrbits : reference to the height of the Photon orbit
 *  \param M : value of the two equal masses m1 and m2
 */
bool MetricExtremeReissnerNordstromDihole::calcHeightOfOtherPhotonOrbits(double &heightOfOtherPhotonOrbits, double M) {
    //int status;
    const gsl_root_fdfsolver_type *T;
    gsl_root_fdfsolver *s;

    //double x0;
    double x = 0.99999999999999;        // initial value

    gsl_function_fdf FDF;
    struct heightOfOtherPhotonOrbits_params params = {M};

    FDF.f = &m4d::heightOfOtherPhotonOrbits;
    FDF.df = &m4d::heightOfOtherPhotonOrbits_deriv;
    FDF.fdf = &m4d::heightOfOtherPhotonOrbits_fdf;
    FDF.params = &params;

    T = gsl_root_fdfsolver_newton;
    s = gsl_root_fdfsolver_alloc(T);
    gsl_root_fdfsolver_set(s, &FDF, x);

    for (int iter = 1; iter <= 100; iter++) {
        //status = gsl_root_fdfsolver_iterate (s);
        gsl_root_fdfsolver_iterate(s);
        //x0 = x;
        x = gsl_root_fdfsolver_root(s);
    }

    heightOfOtherPhotonOrbits = x;

    gsl_root_fdfsolver_free(s);

    return true;
}

/*!  Calculate the velocity on a circular timelike orbit in a plane normal to the z-axis for
 *   z unequal to null and equal masses m1=m2 at initial position
 *
 *  \param pos : current position
 *  \param betaTimelike : reference to the necessary velocity
 *  \param M : mass-value of the two equal masses m1 and m2
 */
bool MetricExtremeReissnerNordstromDihole::calcVelocityForOtherTimelikeCo(const vec4 pos, double &betaTimelike, double M) {
    double rho;
    double z = pos[3];
    calcRadiusForOtherCo(rho, z);

    double rhoSquare = rho * rho;
    double zMinusSquare = (z - 1.) * (z - 1.);
    double zPlusSquare = (z + 1.) * (z + 1.);

    double sqrtMinus = sqrt(rhoSquare + zMinusSquare);
    double sqrtPlus = sqrt(rhoSquare + zPlusSquare);

    double U = 1. + M * (1. / sqrtMinus + 1. / sqrtPlus);
    double dU = -M * rho * (1. / (sqrtMinus * sqrtMinus * sqrtMinus) + 1. / (sqrtPlus * sqrtPlus * sqrtPlus));

    double betaSquare = U / (U + rho * dU) - 1.;

    betaTimelike = sqrt(betaSquare);
    return true;
}

/*!  Calculate the height of the two Photon orbits normal to the z-axis for a z-value
 *   unequal to null and for unequal masses m1 and m2 numerically using the root-finding
 *   algorithm by Brent
 *
 *   The parameters in this Brent-algorithm are the mass-values of the two masses.
 *   The necessary function and its derivative are written further down in this data.
 *
 *  \param heightOfOtherPhotonOrbits1 : reference to the height of one of the Photon orbits
 *  \param heightOfOtherPhotonOrbits2 : reference to the height of the other Photon orbit
 *  \param M1 : value of the mass at z=+1
 *  \param M2 : value of the mass at z=-1
 */
bool MetricExtremeReissnerNordstromDihole::calcHeightOfOtherPhotonOrbitsForUnequalMasses(double &heightOfOtherPhotonOrbits1, double &heightOfOtherPhotonOrbits2, double M1, double M2) {
    //int status;
    const gsl_root_fsolver_type *T;
    gsl_root_fsolver *s;
    double r = 0.;

    double zEqu, zSing;
    calcPointOfEquilibrium(zEqu);
    calcASingularityPoint(zSing);

    double zmin1 = -1. + 0.00000000000001;
    double zmax2 = 1. - 0.00000000000001;
    double zmax1, zmin2;

    if (M1 < M2) {
        zmax1 = zSing - 0.00000000000001;
        zmin2 = zEqu + 0.00000000000001;
    } else {
        zmax1 = zEqu - 0.00000000000001;
        zmin2 = zSing + 0.00000000000001;
    }

    gsl_function F;

    struct heightOfOtherPhotonOrbitsForUnequalMasses_params params = {M1, M2};

    F.function = &m4d::heightOfOtherPhotonOrbitsForUnequalMasses;
    F.params = &params;

    T = gsl_root_fsolver_brent;
    s = gsl_root_fsolver_alloc(T);
    gsl_root_fsolver_set(s, &F, zmin1, zmax1);

    for (int iter = 1; iter <= 100; iter++) {
        //status = gsl_root_fsolver_iterate (s);
        gsl_root_fsolver_iterate(s);
        r = gsl_root_fsolver_root(s);
        zmin1 = gsl_root_fsolver_x_lower(s);
        zmax1 = gsl_root_fsolver_x_upper(s);
    }

    heightOfOtherPhotonOrbits1 = r;


    gsl_root_fsolver_set(s, &F, zmin2, zmax2);

    for (int iter = 1; iter <= 100; iter++) {
        //status = gsl_root_fsolver_iterate (s);
        gsl_root_fsolver_iterate(s);
        r = gsl_root_fsolver_root(s);
        zmin1 = gsl_root_fsolver_x_lower(s);
        zmax1 = gsl_root_fsolver_x_upper(s);
    }

    heightOfOtherPhotonOrbits2 = r;

    gsl_root_fsolver_free(s);
    return true;
}

/*!  Calculate the velocity on a circular timelike orbit in a plane normal to the z-axis for
 *   z unequal to null and unequal masses m1 and m2 at initial position
 *
 *  \param pos : current position
 *  \param betaTimelike : reference to the necessary velocity
 *  \param M1 : mass-value of the mass at z=+1
 *  \param M2 : mass-value of the mass at z=-1
 */
bool MetricExtremeReissnerNordstromDihole::calcVelocityForOtherTimelikeCoForUnequalMasses(const vec4 pos, double &betaTimelike, double M1, double M2) {
    double rho;
    double z = pos[3];
    calcRadiusForOtherCoForUnequalMasses(rho, z, M1, M2);

    double rhoSquare = rho * rho;
    double zMinusSquare = (z - 1.) * (z - 1.);
    double zPlusSquare = (z + 1.) * (z + 1.);

    double sqrtMinus = sqrt(rhoSquare + zMinusSquare);
    double sqrtPlus = sqrt(rhoSquare + zPlusSquare);

    double U = 1. + M1 / sqrtMinus + M2 / sqrtPlus;
    double dU = -M1 * rho / (sqrtMinus * sqrtMinus * sqrtMinus) - M2 * rho / (sqrtPlus * sqrtPlus * sqrtPlus);

    double betaSquare = U / (U + rho * dU) - 1.;

    betaTimelike = sqrt(betaSquare);
    return true;
}



// ********************************* protected methods *****************************
/*!
 */
void MetricExtremeReissnerNordstromDihole::setStandardValues() {
    mInitPos[0] = 0.0;
    mInitPos[1] = 10.0;
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

void MetricExtremeReissnerNordstromDihole::calc_r(const double y[]) {
    r1 = sqrt(y[1] * y[1] + y[2] * y[2] + (y[3] - 1.0) * (y[3] - 1.0));
    r2 = sqrt(y[1] * y[1] + y[2] * y[2] + (y[3] + 1.0) * (y[3] + 1.0));
}

double MetricExtremeReissnerNordstromDihole::calc_U(const double y[]) {
    calc_r(y);
    return 1.0 + mM1 / r1 + mM2 / r2;
}

void MetricExtremeReissnerNordstromDihole::calc_dU(const double y[], double &Ux, double &Uy, double &Uz) {
    double xx = y[1];
    double yy = y[2];
    double zz = y[3];

    double r1hm3 = pow(r1, -3.0);
    double r2hm3 = pow(r2, -3.0);

    Ux = -mM1 * xx * r1hm3 - mM2 * xx * r2hm3;
    Uy = -mM1 * yy * r1hm3 - mM2 * yy * r2hm3;
    Uz = -mM1 * (zz - 1.0) * r1hm3 - mM2 * (zz + 1.0) * r2hm3;
}

void MetricExtremeReissnerNordstromDihole::calc_ddU(const double y[], double &Uxx, double &Uyy, double &Uzz,
        double &Uxy, double &Uxz, double &Uyz) {
    double xx = y[1];
    double yy = y[2];
    double zz = y[3];

    double r1hm3 = pow(r1, -3.0);
    double r2hm3 = pow(r2, -3.0);
    double r1hm5 = pow(r1, -5.0);
    double r2hm5 = pow(r2, -5.0);

    Uxx = 3.0 * mM1 * xx * xx * r1hm5 - mM1 * r1hm3 + 3.0 * mM2 * xx * xx * r2hm5 - mM2 * r2hm3;
    Uyy = 3.0 * mM1 * yy * yy * r1hm5 - mM1 * r1hm3 + 3.0 * mM2 * yy * yy * r2hm5 - mM2 * r2hm3;
    Uzz = 3.0 * mM1 * (zz - 1.0) * (zz - 1.0) * r1hm5 - mM1 * r1hm3 + 3.0 * mM2 * (zz + 1.0) * (zz + 1.0) * r2hm5 - mM2 * r2hm3;

    Uxy = 3.0 * mM1 * xx * yy * r1hm5 + 3.0 * mM2 * xx * yy * r2hm5;
    Uxz = 3.0 * mM1 * (zz - 1.0) * xx * r1hm5 + 3.0 * mM2 * (zz + 1.0) * xx * r2hm5;
    Uyz = 3.0 * mM1 * (zz - 1.0) * yy * r1hm5 + 3.0 * mM2 * (zz + 1.0) * yy * r2hm5;
}



// ************************ for some numerical calculations ************************
/*!  The quantities radiusIscoXY and radiusIucoXY in the z=0-plane for equal masses
 *   m1 and m2 follow from the two roots of the function inflectionPoints(x) with
 *       inflectionPoints(x) = x^3 - 6M*x^(5/2) + 3*x^2 + 22M*x^(3/2) + 16M^2
 *   in the interval 1 < x < inf. M is equal m1 and m2 and is the only parameter
 *   from the struct inflectionPoints_params from the header-file.
 *
 *   The function inflectionPoints_deriv(x) is the derivative of inflectionPoints(x),
 *   inflectionPoints_fdf are both in one. The radius rho coresponds with x via
 *   rho = sqrt(x-1), where x came from a substitution.
 *   The name inflectionPoints was chosen, because this function describes the
 *   inflection points of the appendant effective potential.
 */
double
inflectionPoints(double x, void *params) {
    struct inflectionPoints_params *p  = (struct inflectionPoints_params *) params;

    double M = p->M;
    double sqrtx = sqrt(x);

    return x * x * x - 6.0 * M * x * x * sqrtx + 3.0 * x * x + 22.0 * M * x * sqrtx + 16.0 * M * M;
}
double
inflectionPoints_deriv(double x, void *params) {
    struct inflectionPoints_params *p = (struct inflectionPoints_params *) params;

    double M = p->M;
    double sqrtx = sqrt(x);

    return 3.0 * x * x - 15.0 * M * x * sqrtx + 6.0 * x + 33.0 * M * sqrtx;
}
void
inflectionPoints_fdf(double x, void *params, double *y, double *dy) {
    struct inflectionPoints_params *p = (struct inflectionPoints_params *) params;

    double M = p->M;
    double sqrtx = sqrt(x);

    *y = x * x * x - 6.0 * M * x * x * sqrtx + 3.0 * x * x + 22.0 * M * x * sqrtx + 16.0 * M * M;
    *dy = 3.0 * x * x - 15.0 * M * x * sqrtx + 6.0 * x + 33.0 * M * sqrtx;
}

/*!  The quantity periodicAlongX which is calculated in the method calcPeriodicAlongX(...)
 *   has to be calculated by a integration. This integration is executed numerically, the
 *   integrand of the relevant integration is given here. The integration variable is x, M
 *   is the mass of the two equal masses. M is a parameter of the function, declared in a
 *   struct in the header-file.
 */
double integrandForPeriodicAlongX(double x, void *params) {
    struct integrandForPeriodicAlongX_params *par = (struct integrandForPeriodicAlongX_params*)params;

    double M = par->M;
    double x0 = par->x0;

    double Ux = 1. + 2.*M / sqrt(x * x + 1.);
    double Ux0 = 1. + 2.*M / sqrt(x0 * x0 + 1.);
    double Veffx = 0.5 / (Ux * Ux);
    double Veffx0 = 0.5 / (Ux0 * Ux0);
    double denominator = Veffx0 - Veffx;

    return sqrt(2. / denominator);
}

/*!  The height of a Photon orbit for z unequal to zero for equal masses m1 and m2 is given
 *   by the roots of an "ugly" function, defined here in heightOfOtherPhotonOrbits(...).
 *   The only parameter M is the mass-value of the both equal masses and is declared in a
 *   struct in the header-file.
 *
 *   The calculation of the roots is done by the Newton-algorithm, so the derivative of the
 *   function is necessary, too. The derivative is defined in heightOfOtherPhotonOrbits_deriv(...).
 *   heightOfOtherPhotonOrbits_fdf(...) contains both, function and its derivative.
 */
double heightOfOtherPhotonOrbits(double z, void *params) {
    struct heightOfOtherPhotonOrbits_params *p  = (struct heightOfOtherPhotonOrbits_params *) params;

    double M = p->M;

    double w1 = 1. + z;
    double w2 = 1. - z;
    double w3 = w1 * w1;
    double w4 = w2 * w2;
    double w5 = w4 / w3;
    double w6 = pow(w5, 1. / 3.);
    double w7 = 1. / (1. - w6);
    double w8 = z * w7;

    double w9 = (4.*w8 - w3);
    double w9b = M / pow(w8 - z, 3. / 2.);
    double w9c = M / pow(w8, 3. / 2.);
    double w10 = w9b + w9c;
    double w11 = (2. + M / pow(w8 - z, 1. / 2.) + M / pow(w8, 1. / 2.));

    return w9 * w10 - 2.*w11;
}
double heightOfOtherPhotonOrbits_deriv(double z, void *params) {
    struct heightOfOtherPhotonOrbits_params *p  = (struct heightOfOtherPhotonOrbits_params *) params;

    double M = p->M;

    double w1 = 1. + z;
    double w2 = 1. - z;
    double w3 = w1 * w1;
    double w4 = w2 * w2;
    double w5 = w4 / w3;
    double w6 = pow(w5, 1. / 3.);
    double w7 = 1. / (1. - w6);
    double w8 = z * w7;

    double w12 = 1. / w7;
    double w13 = sqrt(w6);
    double w14 = sqrt(w8 - z);

    double summand1 = 6.*(w7 - 4.*z * w7 * w7 / (3.*w13 * w3)) / pow(w8, 3. / 2.);
    double summand2 = (6.*w12 - 8.*z - 6.*z * z * w12) / (z * w1 * w2 * w12 * w14);
    double summand3 = 6.*(4.*w7 - 2.*w1 - 16.*z * w7 * w7 / (3.*w3 * w13)) * (1. / pow(w8, 3. / 2.) + 1. / (w14 * w14 * w14));
    double summand4 = 9.*(4.*z * w7 - w3) * ((3.*w12 - 4.*z - 3.*z * z * w12) / (3.*z * z * (z * z - 1.) * w6 * w14) - (w7 - 4.*z * w7 * w7 / (3.*w13 * w3)) / pow(w8, 5. / 2.));

    return (M / 6.) * (summand1 + summand2 + summand3 + summand4);
}
void heightOfOtherPhotonOrbits_fdf(double z, void *params, double *y, double *dy) {
    struct heightOfOtherPhotonOrbits_params *p  = (struct heightOfOtherPhotonOrbits_params *) params;

    double M = p->M;

    double w1 = 1. + z;
    double w2 = 1. - z;
    double w3 = w1 * w1;
    double w4 = w2 * w2;
    double w5 = w4 / w3;
    double w6 = pow(w5, 1. / 3.);
    double w7 = 1. / (1. - w6);
    double w8 = z * w7;

    double w9 = (4.*w8 - w3);
    double w9b = M / pow(w8 - z, 3. / 2.);
    double w9c = M / pow(w8, 3. / 2.);
    double w10 = w9b + w9c;
    double w11 = (2. + M / pow(w8 - z, 1. / 2.) + M / pow(w8, 1. / 2.));

    double w12 = 1. / w7;
    double w13 = sqrt(w6);
    double w14 = sqrt(w8 - z);

    double summand1 = 6.*(w7 - 4.*z * w7 * w7 / (3.*w13 * w3)) / pow(w8, 3. / 2.);
    double summand2 = (6.*w12 - 8.*z - 6.*z * z * w12) / (z * w1 * w2 * w12 * w14);
    double summand3 = 6.*(4.*w7 - 2.*w1 - 16.*z * w7 * w7 / (3.*w3 * w13)) * (1. / pow(w8, 3. / 2.) + 1. / (w14 * w14 * w14));
    double summand4 = 9.*(4.*z * w7 - w3) * ((3.*w12 - 4.*z - 3.*z * z * w12) / (3.*z * z * (z * z - 1.) * w6 * w14) - (w7 - 4.*z * w7 * w7 / (3.*w13 * w3)) / pow(w8, 5. / 2.));

    *y = w9 * w10 - 2.*w11;
    *dy = (M / 6.) * (summand1 + summand2 + summand3 + summand4);
}

/*!  The height of a Photon orbit for z unequal to zero for unequal masses m1 and m2 is given
 *   by the roots of a "very ugly" function, defined here in
 *   heightOfOtherPhotonOrbitsForUnequalMasses(...). The two parameters M1 and M2 are the
 *   mass-values of the both unequal masses and are declared in a struct in the header-file.
 *
 *   Because a non-defined area and the "ugliness" of the function, the convergence of a
 *   Newton-algorithm can not be guaranteed, so it is used a simple bisektion-algorithm for
 *   root-finding. Also the declaration of the derivative of the function is not necessary.
 */
double heightOfOtherPhotonOrbitsForUnequalMasses(double z, void *params) {
    struct heightOfOtherPhotonOrbitsForUnequalMasses_params *p  = (struct heightOfOtherPhotonOrbitsForUnequalMasses_params *) params;

    double M1 = p->M1;
    double M2 = p->M2;

    double w1 = 1. + z;
    double w2 = 1. - z;
    double w3 = pow(M1 / M2, 2. / 3.);
    double w4 = pow(w2 / w1, 2. / 3.);
    double w5 = z / (1. - w3 * w4);

    return (4.*w5 - w1 * w1) * (M2 / pow(w5, 3. / 2.) + M1 / pow(w5 - z, 3. / 2.)) - 4. - 2.*M2 / sqrt(w5) - 2.*M1 / sqrt(w5 - z);
}

} // end namespace m4d
