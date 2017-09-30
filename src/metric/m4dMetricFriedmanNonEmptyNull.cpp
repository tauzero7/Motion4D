// -------------------------------------------------------------------------------
/*
   m4dMetricFriedmanNonEmptyNull.cpp

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

#include "m4dMetricFriedmanNonEmptyNull.h"

namespace m4d {

#define eps 1.0e-6
#define edr 0.333333333333333

/*! Standard constructor for the Kottler metric.
 *
 * \param  mass : mass of the black hole.
 * \param  k : Friedman parameter.
 */
MetricFriedmanNonEmptyNull::MetricFriedmanNonEmptyNull(double mass, double k) {
    mMetricName  = "FriedmanNonEmptyNull";
    setCoordType(enum_coordinate_spherical);

    mPhysicalUnits = enum_physical_constants_geom;
    mSpeedOfLight = 1.0;
    mGravConstant = 1.0;

    addParam("mass", mass);
    addParam("k", k);

    mC = 4.0 * mGravConstant * mass / (3.0 * M_PI);
    mK = 0;
    if (k >= 1.0) {
        mK = 1;
    } else if (k <= -1.0) {
        mK = -1;
    }

    mDrawTypes.push_back(enum_draw_twoplusone);

    mLocTeds.push_back(enum_nat_tetrad_static);
    setStandardValues();
}

MetricFriedmanNonEmptyNull::~MetricFriedmanNonEmptyNull() {
}


// *********************************** public methods ******************************
/*! Calculate the contravariant metric components at position 'pos'.
 *
 *  \param pos : pointer to position.
 */
bool MetricFriedmanNonEmptyNull::calculateMetric(const double* pos) {
    double chi   = pos[0];
    double r     = pos[1];
    double theta = pos[2];
    double k = (double)mK;

    double t1 = calc_R(chi);
    double t2 = t1 * t1;
    double t3 = r * r;
    double t7 = pow(1.0 + k * t3 / 4.0, 2.0);
    double t9 = t2 / t7;
    double t11 = sin(theta);
    double t12 = t11 * t11;

    g_compts[0][0] = -t2;
    g_compts[0][1] = 0.0;
    g_compts[0][2] = 0.0;
    g_compts[0][3] = 0.0;
    g_compts[1][0] = 0.0;
    g_compts[1][1] = t9;
    g_compts[1][2] = 0.0;
    g_compts[1][3] = 0.0;
    g_compts[2][0] = 0.0;
    g_compts[2][1] = 0.0;
    g_compts[2][2] = t9 * t3;
    g_compts[2][3] = 0.0;
    g_compts[3][0] = 0.0;
    g_compts[3][1] = 0.0;
    g_compts[3][2] = 0.0;
    g_compts[3][3] = t9 * t3 * t12;

    return true;
}

/*! Calculate the Christoffel symbols of the second kind at position 'pos'.
 *
 *  \param pos : pointer to position.
 */
bool MetricFriedmanNonEmptyNull::calculateChristoffels(const double* pos) {
    double chi   = pos[0];
    double r     = pos[1];
    double theta = pos[2];
    double k = (double)mK;

    double t1 = calc_R(chi);
    double t2 = 1 / t1;
    double t3 = calc_dR(chi); // diff(R(chi),chi);
    double t4 = t2 * t3;
    double t5 = r * r;
    double t6 = k * t5;
    double t7 = 4.0 + t6;
    double t8 = t7 * t7;
    double t9 = 1 / t8;
    double t12 = 1 / t7;
    double t18 = t6 - 4.0;
    double t19 = t12 / r * t18;
    double t20 = t2 * t5;
    double t24 = t12 * r;
    double t26 = sin(theta);
    double t28 = cos(theta);
    double t29 = 1 / t26 * t28;
    double t30 = t26 * t26;

    christoffel[0][0][0] = t4;
    christoffel[0][0][1] = 0.0;
    christoffel[0][0][2] = 0.0;
    christoffel[0][0][3] = 0.0;
    christoffel[0][1][0] = 0.0;
    christoffel[0][1][1] = t4;
    christoffel[0][1][2] = 0.0;
    christoffel[0][1][3] = 0.0;
    christoffel[0][2][0] = 0.0;
    christoffel[0][2][1] = 0.0;
    christoffel[0][2][2] = t4;
    christoffel[0][2][3] = 0.0;
    christoffel[0][3][0] = 0.0;
    christoffel[0][3][1] = 0.0;
    christoffel[0][3][2] = 0.0;
    christoffel[0][3][3] = t4;
    christoffel[1][0][0] = 0.0;
    christoffel[1][0][1] = t4;
    christoffel[1][0][2] = 0.0;
    christoffel[1][0][3] = 0.0;
    christoffel[1][1][0] = 16.0 * t4 * t9;
    christoffel[1][1][1] = -2.0 * t12 * k * r;
    christoffel[1][1][2] = 0.0;
    christoffel[1][1][3] = 0.0;
    christoffel[1][2][0] = 0.0;
    christoffel[1][2][1] = 0.0;
    christoffel[1][2][2] = -t19;
    christoffel[1][2][3] = 0.0;
    christoffel[1][3][0] = 0.0;
    christoffel[1][3][1] = 0.0;
    christoffel[1][3][2] = 0.0;
    christoffel[1][3][3] = -t19;
    christoffel[2][0][0] = 0.0;
    christoffel[2][0][1] = 0.0;
    christoffel[2][0][2] = t4;
    christoffel[2][0][3] = 0.0;
    christoffel[2][1][0] = 0.0;
    christoffel[2][1][1] = 0.0;
    christoffel[2][1][2] = -t19;
    christoffel[2][1][3] = 0.0;
    christoffel[2][2][0] = 16.0 * t20 * t9 * t3;
    christoffel[2][2][1] = t24 * t18;
    christoffel[2][2][2] = 0.0;
    christoffel[2][2][3] = 0.0;
    christoffel[2][3][0] = 0.0;
    christoffel[2][3][1] = 0.0;
    christoffel[2][3][2] = 0.0;
    christoffel[2][3][3] = t29;
    christoffel[3][0][0] = 0.0;
    christoffel[3][0][1] = 0.0;
    christoffel[3][0][2] = 0.0;
    christoffel[3][0][3] = t4;
    christoffel[3][1][0] = 0.0;
    christoffel[3][1][1] = 0.0;
    christoffel[3][1][2] = 0.0;
    christoffel[3][1][3] = -t19;
    christoffel[3][2][0] = 0.0;
    christoffel[3][2][1] = 0.0;
    christoffel[3][2][2] = 0.0;
    christoffel[3][2][3] = t29;
    christoffel[3][3][0] = 16.0 * t20 * t30 * t3 * t9;
    christoffel[3][3][1] = t24 * t30 * t18;
    christoffel[3][3][2] = -t26 * t28;
    christoffel[3][3][3] = 0.0;

    return true;
}

/*! Calculate Jacobi matrix.
 *
 *  \param pos : pointer to position.
 */
bool MetricFriedmanNonEmptyNull::calculateChrisD(const double* pos) {
    double chi   = pos[0];
    double r     = pos[1];
    double theta = pos[2];
    double k = (double)mK;

    double t1 = calc_dR(chi);  // diff(R(chi),chi);
    double t2 = t1 * t1;
    double t3 = calc_ddR(chi); // diff(diff(R(chi),chi),chi);
    double t4 = calc_R(chi);
    double t6 = -t2 + t3 * t4;
    double t7 = t4 * t4;
    double t8 = 1 / t7;
    double t9 = t6 * t8;
    double t10 = r * r;
    double t11 = k * t10;
    double t12 = 4.0 + t11;
    double t13 = t12 * t12;
    double t14 = 1 / t13;
    double t15 = t9 * t14;
    double t17 = 1 / t4;
    double t20 = 1 / t13 / t12;
    double t25 = t11 - 4.0;
    double t29 = k * k;
    double t30 = t10 * t10;
    double t31 = t29 * t30;
    double t32 = 16.0 * t11;
    double t36 = (t31 - t32 - 16.0) * t14 / t10;
    double t43 = t25 * t17 * t20;
    double t46 = t31 + t32 - 16.0;
    double t48 = cos(theta);
    double t49 = t48 * t48;
    double t50 = sin(theta);
    double t51 = t50 * t50;
    double t54 = (t49 + t51) / t51;

    chrisD[0][0][0][0] = t9;
    chrisD[0][0][0][1] = 0.0;
    chrisD[0][0][0][2] = 0.0;
    chrisD[0][0][0][3] = 0.0;
    chrisD[0][0][1][0] = 0.0;
    chrisD[0][0][1][1] = 0.0;
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
    chrisD[0][1][0][1] = 0.0;
    chrisD[0][1][0][2] = 0.0;
    chrisD[0][1][0][3] = 0.0;
    chrisD[0][1][1][0] = t9;
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
    chrisD[0][2][2][0] = t9;
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
    chrisD[0][3][3][0] = t9;
    chrisD[0][3][3][1] = 0.0;
    chrisD[0][3][3][2] = 0.0;
    chrisD[0][3][3][3] = 0.0;
    chrisD[1][0][0][0] = 0.0;
    chrisD[1][0][0][1] = 0.0;
    chrisD[1][0][0][2] = 0.0;
    chrisD[1][0][0][3] = 0.0;
    chrisD[1][0][1][0] = t9;
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
    chrisD[1][1][0][0] = 16.0 * t15;
    chrisD[1][1][0][1] = -64.0 * t17 * t1 * t20 * k * r;
    chrisD[1][1][0][2] = 0.0;
    chrisD[1][1][0][3] = 0.0;
    chrisD[1][1][1][0] = 0.0;
    chrisD[1][1][1][1] = 2.0 * k * t25 * t14;
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
    chrisD[1][2][2][1] = t36;
    chrisD[1][2][2][2] = 0.0;
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
    chrisD[1][3][3][1] = t36;
    chrisD[1][3][3][2] = 0.0;
    chrisD[1][3][3][3] = 0.0;
    chrisD[2][0][0][0] = 0.0;
    chrisD[2][0][0][1] = 0.0;
    chrisD[2][0][0][2] = 0.0;
    chrisD[2][0][0][3] = 0.0;
    chrisD[2][0][1][0] = 0.0;
    chrisD[2][0][1][1] = 0.0;
    chrisD[2][0][1][2] = 0.0;
    chrisD[2][0][1][3] = 0.0;
    chrisD[2][0][2][0] = t9;
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
    chrisD[2][1][2][1] = t36;
    chrisD[2][1][2][2] = 0.0;
    chrisD[2][1][2][3] = 0.0;
    chrisD[2][1][3][0] = 0.0;
    chrisD[2][1][3][1] = 0.0;
    chrisD[2][1][3][2] = 0.0;
    chrisD[2][1][3][3] = 0.0;
    chrisD[2][2][0][0] = 16.0 * t10 * t6 * t8 * t14;
    chrisD[2][2][0][1] = -32.0 * r * t1 * t43;
    chrisD[2][2][0][2] = 0.0;
    chrisD[2][2][0][3] = 0.0;
    chrisD[2][2][1][0] = 0.0;
    chrisD[2][2][1][1] = t46 * t14;
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
    chrisD[2][3][3][2] = -t54;
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
    chrisD[3][0][3][0] = t9;
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
    chrisD[3][1][3][1] = t36;
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
    chrisD[3][2][3][2] = -t54;
    chrisD[3][2][3][3] = 0.0;
    chrisD[3][3][0][0] = 16.0 * t10 * t51 * t15;
    chrisD[3][3][0][1] = -32.0 * r * t51 * t1 * t43;
    chrisD[3][3][0][2] = 32.0 * t17 * t10 * t50 * t1 * t14 * t48;
    chrisD[3][3][0][3] = 0.0;
    chrisD[3][3][1][0] = 0.0;
    chrisD[3][3][1][1] = t51 * t46 * t14;
    chrisD[3][3][1][2] = 2.0 / t12 * r * t50 * t25 * t48;
    chrisD[3][3][1][3] = 0.0;
    chrisD[3][3][2][0] = 0.0;
    chrisD[3][3][2][1] = 0.0;
    chrisD[3][3][2][2] = -t49 + t51;
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
void MetricFriedmanNonEmptyNull::localToCoord(const double* pos, const double* ldir, double* dir,
        enum_nat_tetrad_type) {
    double chi   = pos[0];
    double r     = pos[1];
    double theta = pos[2];
    double k = (double)mK;

    double brk = 1.0 + 0.25 * k * r * r;
    double edR = 1.0 / calc_R(chi);

    dir[0] = ldir[0] * edR;
    dir[1] = ldir[1] * edR * brk;
    dir[2] = ldir[2] * edR * brk / r;
    dir[3] = ldir[3] * edR * brk / (r * sin(theta));
}

/*! Transform coordinate 4-direction to local 4-direction.
 *
 *  \param  pos  :  pointer to position array.
 *  \param  cdir :  pointer to coordinate direction.
 *  \param  ldir :  pointer to calculated local direction array.
 *  \param  type :  type of tetrad.
 */
void MetricFriedmanNonEmptyNull::coordToLocal(const double* pos, const double* cdir, double* ldir,
        enum_nat_tetrad_type) {
    double chi   = pos[0];
    double r     = pos[1];
    double theta = pos[2];
    double k = (double)mK;

    double edbrk = 1.0 / (1.0 + 0.25 * k * r * r);
    double R     = calc_R(chi);

    ldir[0] = cdir[0] * R;
    ldir[1] = cdir[1] * R * edbrk;
    ldir[2] = cdir[2] * R * edbrk * r;
    ldir[3] = cdir[3] * R * edbrk * r * sin(theta);
}


/*! Test break condition.
 *
 *  \param pos    : pointer to position array.
 *  \return true  : radial position r < 0.0 or  r^2<=(1.0+eps)*rs^2.
 *  \return false : position is valid.
 */
bool MetricFriedmanNonEmptyNull::breakCondition(const double* pos) {
    bool br = false;

    double r = pos[1];
    double k = (double)mK;
    if (mK == -1)
        if (fabs(1.0 + 0.25 * k * r * r) < 1e-12) {
            br = true;
        }
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
double MetricFriedmanNonEmptyNull::testConstraint(const double y[], const double kappa) {
    double chi   = y[0];
    double r     = y[1];
    double theta = y[2];
    double k = (double)mK;

    double R     = calc_R(chi);
    double edbrk = 1.0 / (1.0 + 0.25 * k * r * r);

    double dchi = y[4];
    double dr   = y[5];
    double dth  = y[6];
    double dph  = y[7];

    double sum = -kappa;
    sum += R * R * (-dchi * dchi + edbrk * edbrk * (dr * dr + r * r * (dth * dth + sin(theta) * sin(theta) * dph * dph)));
    return sum;
}

/*! Set parameter 'pName' to 'val'.
 *
 *  Set 'mass' or 'lambda' parameter.
 */
bool MetricFriedmanNonEmptyNull::setParam(const char* pName, double val) {
    Metric::setParam(pName, val);

    if (strcmp(pName,"k") == 0) {
        mK = 0;
        if (val >= 1.0) {
            mK = 1;
        } else if (val <= -1.0) {
            mK = -1;
        }
    }
    else if (strcmp(pName,"mass") == 0) {
        mC = 4.0 * mGravConstant * val / (3.0 * M_PI);
    }
    return true;
}

/*! Transform point p to 2+1 coordinates.
 *
 *  \param  p  : point in proper metric coordinates.
 *  \param  cp : reference to transformed point.
 *  \return true : success.
 */
bool MetricFriedmanNonEmptyNull::transToTwoPlusOne(vec4 p, vec4 &cp) {
    vec4 tp;
    TransCoordinates::toCartesianCoord(mCoordType, p, tp);
    cp = vec4(tp[0], tp[1], tp[2], tp[0]);
    return true;
}

/*! Generate report.
 */
bool MetricFriedmanNonEmptyNull::report(const vec4 , const vec4 , std::string &text) {
    std::stringstream ss;
    ss << "Report for non-empty Friedman metric\n\tcoordinate : (t,r,theta,phi)\n";
    ss << "---------------------------------------------------------\n";
    ss << "  physical units ................................. no\n";

    text = ss.str();
    return true;
}


// ********************************* protected methods *****************************
/*!
 */
void MetricFriedmanNonEmptyNull::setStandardValues() {
    mInitPos[0] = 1.0;
    mInitPos[1] = 10.0;
    mInitPos[2] = M_PI_2;
    mInitPos[3] = 0.0;
    mInitDir[0] = 1.0;
    mInitDir[1] = 0.0;
    mInitDir[2] = 0.0;

    mCoordNames[0] = std::string("chi");
    mCoordNames[1] = std::string("r");
    mCoordNames[2] = std::string("theta");
    mCoordNames[3] = std::string("phi");
}

/*! Calculate conformal factor R(chi).
 * \param chi: conformal coordinate.
 */
double MetricFriedmanNonEmptyNull::calc_R(const double chi) {
    if (mK == 0) {
        return 2.25 * mC * chi * chi;
    } else if (mK == 1) {
        return 0.5 * mC * (1.0 - cos(chi));
    }

    return 0.5 * mC * (cosh(chi) - 1.0);
}

/*! Calculate first derivative of conformal factor R(chi).
 * \param chi: conformal coordinate.
 */
double MetricFriedmanNonEmptyNull::calc_dR(const double chi) {
    if (mK == 0) {
        return 4.5 * mC * chi;
    } else if (mK == 1) {
        return 0.5 * mC * sin(chi);
    }

    return 0.5 * mC * sinh(chi);
}

/*! Calculate second derivative of conformal factor R(chi).
 * \param chi: conformal coordinate.
 */
double MetricFriedmanNonEmptyNull::calc_ddR(const double chi) {
    if (mK == 0) {
        return 4.5 * mC;
    } else if (mK == 1) {
        return 0.5 * mC * cos(chi);
    }

    return 0.5 * mC * cosh(chi);
}

} // end namespace m4d
