// -------------------------------------------------------------------------------
/*
   m4dMetricEddFinkIn.cpp

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

#include "m4dMetricEddFinkIn.h"

#include <gsl/gsl_errno.h>
#include <gsl/gsl_roots.h>

namespace m4d {

double tau_of_r(double r, void* params)
{
    struct EddFinkParams* par = static_cast<struct EddFinkParams*>(params);
    double rs = par->rs;
    double r0 = par->r0;
    double tau = par->tau;
    return r * r0 / rs * sqrt(rs / r - rs / r0) + r0 * pow(r0 / rs, 0.5) * atan(sqrt(r0 / r - 1.0)) - tau;
}

#define eps 1.0e-6

/*! Standard constructor for the Eddington-Finkelstein metric.
 *
 * \param  mass : mass of the black hole.
 */
MetricEddFinkIn::MetricEddFinkIn(double mass)
{
    mMetricName = "EddFinkIn";
    setCoordType(enum_coordinate_spherical);

    mPhysicalUnits = enum_physical_constants_geom;
    mSpeedOfLight = 1.0;
    mGravConstant = 1.0;

    addParam("mass", mass);
    mMass = mass;
    rs = 2.0 * mGravConstant * mMass / (mSpeedOfLight * mSpeedOfLight);

    setStandardValues();

    mLocTeds.push_back(enum_nat_tetrad_static);
    mLocTeds.push_back(enum_nat_tetrad_freefall);
}

MetricEddFinkIn::~MetricEddFinkIn() {}

// *********************************** public methods ******************************
/*! Calculate the contravariant metric components at position 'pos'.
 *
 *  \param pos : pointer to position.
 */
bool MetricEddFinkIn::calculateMetric(const double* pos)
{
    double r = pos[1];
    double theta = pos[2];

    double c = mSpeedOfLight;

    double t1 = c * c;
    double t6 = r * r;
    double t7 = sin(theta);
    double t8 = t7 * t7;

    g_compts[0][0] = -t1 + t1 * rs / r;
    g_compts[0][1] = c;
    g_compts[0][2] = 0.0;
    g_compts[0][3] = 0.0;
    g_compts[1][0] = c;
    g_compts[1][1] = 0.0;
    g_compts[1][2] = 0.0;
    g_compts[1][3] = 0.0;
    g_compts[2][0] = 0.0;
    g_compts[2][1] = 0.0;
    g_compts[2][2] = t6;
    g_compts[2][3] = 0.0;
    g_compts[3][0] = 0.0;
    g_compts[3][1] = 0.0;
    g_compts[3][2] = 0.0;
    g_compts[3][3] = t6 * t8;

    return true;
}

/*! Calculate the Christoffel symbols of the second kind at position 'pos'.
 *
 *  \param pos : pointer to position.
 */
bool MetricEddFinkIn::calculateChristoffels(const double* pos)
{
    double r = pos[1];
    double theta = pos[2];

    double c = mSpeedOfLight;

    double t2 = r * r;
    double t5 = c * rs / t2 / 2.0;
    double t6 = r - rs;
    double t10 = c * c;
    double t14 = 1 / r;
    double t16 = 1 / c * r;
    double t17 = sin(theta);
    double t19 = cos(theta);
    double t20 = 1 / t17 * t19;
    double t21 = t17 * t17;

    christoffel[0][0][0] = t5;
    christoffel[0][0][1] = t6 / t2 / r * t10 * rs / 2.0;
    christoffel[0][0][2] = 0.0;
    christoffel[0][0][3] = 0.0;
    christoffel[0][1][0] = 0.0;
    christoffel[0][1][1] = -t5;
    christoffel[0][1][2] = 0.0;
    christoffel[0][1][3] = 0.0;
    christoffel[0][2][0] = 0.0;
    christoffel[0][2][1] = 0.0;
    christoffel[0][2][2] = 0.0;
    christoffel[0][2][3] = 0.0;
    christoffel[0][3][0] = 0.0;
    christoffel[0][3][1] = 0.0;
    christoffel[0][3][2] = 0.0;
    christoffel[0][3][3] = 0.0;
    christoffel[1][0][0] = 0.0;
    christoffel[1][0][1] = -t5;
    christoffel[1][0][2] = 0.0;
    christoffel[1][0][3] = 0.0;
    christoffel[1][1][0] = 0.0;
    christoffel[1][1][1] = 0.0;
    christoffel[1][1][2] = 0.0;
    christoffel[1][1][3] = 0.0;
    christoffel[1][2][0] = 0.0;
    christoffel[1][2][1] = 0.0;
    christoffel[1][2][2] = t14;
    christoffel[1][2][3] = 0.0;
    christoffel[1][3][0] = 0.0;
    christoffel[1][3][1] = 0.0;
    christoffel[1][3][2] = 0.0;
    christoffel[1][3][3] = t14;
    christoffel[2][0][0] = 0.0;
    christoffel[2][0][1] = 0.0;
    christoffel[2][0][2] = 0.0;
    christoffel[2][0][3] = 0.0;
    christoffel[2][1][0] = 0.0;
    christoffel[2][1][1] = 0.0;
    christoffel[2][1][2] = t14;
    christoffel[2][1][3] = 0.0;
    christoffel[2][2][0] = -t16;
    christoffel[2][2][1] = -t6;
    christoffel[2][2][2] = 0.0;
    christoffel[2][2][3] = 0.0;
    christoffel[2][3][0] = 0.0;
    christoffel[2][3][1] = 0.0;
    christoffel[2][3][2] = 0.0;
    christoffel[2][3][3] = t20;
    christoffel[3][0][0] = 0.0;
    christoffel[3][0][1] = 0.0;
    christoffel[3][0][2] = 0.0;
    christoffel[3][0][3] = 0.0;
    christoffel[3][1][0] = 0.0;
    christoffel[3][1][1] = 0.0;
    christoffel[3][1][2] = 0.0;
    christoffel[3][1][3] = t14;
    christoffel[3][2][0] = 0.0;
    christoffel[3][2][1] = 0.0;
    christoffel[3][2][2] = 0.0;
    christoffel[3][2][3] = t20;
    christoffel[3][3][0] = -t16 * t21;
    christoffel[3][3][1] = -t6 * t21;
    christoffel[3][3][2] = -t17 * t19;
    christoffel[3][3][3] = 0.0;

    return true;
}

/*! Calculate Jacobi matrix.
 *
 *  \param pos : pointer to position.
 */
bool MetricEddFinkIn::calculateChrisD(const double* pos)
{
    double r = pos[1];
    double theta = pos[2];

    double c = mSpeedOfLight;

    double t2 = r * r;
    double t5 = c * rs / t2 / r;
    double t6 = c * c;
    double t11 = t2 * t2;
    double t16 = 1 / t2;
    double t17 = 1 / c;
    double t18 = cos(theta);
    double t19 = t18 * t18;
    double t20 = sin(theta);
    double t21 = t20 * t20;
    double t24 = (t19 + t21) / t21;

    chrisD[0][0][0][0] = 0.0;
    chrisD[0][0][0][1] = -t5;
    chrisD[0][0][0][2] = 0.0;
    chrisD[0][0][0][3] = 0.0;
    chrisD[0][0][1][0] = 0.0;
    chrisD[0][0][1][1] = -t6 * rs * (2.0 * r - 3.0 * rs) / t11 / 2.0;
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
    chrisD[0][1][1][0] = 0.0;
    chrisD[0][1][1][1] = t5;
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
    chrisD[1][0][0][1] = 0.0;
    chrisD[1][0][0][2] = 0.0;
    chrisD[1][0][0][3] = 0.0;
    chrisD[1][0][1][0] = 0.0;
    chrisD[1][0][1][1] = t5;
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
    chrisD[1][2][2][1] = -t16;
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
    chrisD[1][3][3][1] = -t16;
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
    chrisD[2][1][2][1] = -t16;
    chrisD[2][1][2][2] = 0.0;
    chrisD[2][1][2][3] = 0.0;
    chrisD[2][1][3][0] = 0.0;
    chrisD[2][1][3][1] = 0.0;
    chrisD[2][1][3][2] = 0.0;
    chrisD[2][1][3][3] = 0.0;
    chrisD[2][2][0][0] = 0.0;
    chrisD[2][2][0][1] = -t17;
    chrisD[2][2][0][2] = 0.0;
    chrisD[2][2][0][3] = 0.0;
    chrisD[2][2][1][0] = 0.0;
    chrisD[2][2][1][1] = -1.0;
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
    chrisD[2][3][3][2] = -t24;
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
    chrisD[3][1][3][1] = -t16;
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
    chrisD[3][2][3][2] = -t24;
    chrisD[3][2][3][3] = 0.0;
    chrisD[3][3][0][0] = 0.0;
    chrisD[3][3][0][1] = -t17 * t21;
    chrisD[3][3][0][2] = -2.0 * t17 * r * t20 * t18;
    chrisD[3][3][0][3] = 0.0;
    chrisD[3][3][1][0] = 0.0;
    chrisD[3][3][1][1] = -t21;
    chrisD[3][3][1][2] = -2.0 * (r - rs) * t20 * t18;
    chrisD[3][3][1][3] = 0.0;
    chrisD[3][3][2][0] = 0.0;
    chrisD[3][3][2][1] = 0.0;
    chrisD[3][3][2][2] = -t19 + t21;
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
void MetricEddFinkIn::localToCoord(const double* pos, const double* ldir, double* dir, enum_nat_tetrad_type type)
{
    double r = pos[1];
    double theta = pos[2];
    double c = mSpeedOfLight;

    if (type != enum_nat_tetrad_freefall) {
        double w = sqrt(1.0 - rs / r);
        dir[0] = (ldir[0] + ldir[1]) / w / c;
        dir[1] = ldir[1] * w;
    }
    else {
        double w = sqrt(rs / r);
        dir[0] = (ldir[0] + ldir[1]) / (1.0 + w) / c;
        dir[1] = -w * ldir[0] + ldir[1];
    }
    dir[2] = ldir[2] / r;
    dir[3] = ldir[3] / (r * sin(theta));
}

/*! Transform coordinate 4-direction to local 4-direction.
 *
 *  \param  pos  :  pointer to position array.
 *  \param  cdir :  pointer to coordinate direction.
 *  \param  ldir :  pointer to calculated local direction array.
 *  \param  type :  type of tetrad.
 */
void MetricEddFinkIn::coordToLocal(const double* pos, const double* cdir, double* ldir, enum_nat_tetrad_type type)
{
    double r = pos[1];
    double theta = pos[2];
    double c = mSpeedOfLight;

    if (type != enum_nat_tetrad_freefall) {
        double w = sqrt(1.0 - rs / r);
        ldir[0] = c * cdir[0] * w - cdir[1] / w;
        ldir[1] = cdir[1] / w;
    }
    else {
        double w = sqrt(rs / r);
        ldir[0] = c * cdir[0] - cdir[1] / (1.0 + w);
        ldir[1] = c * w * cdir[0] + cdir[1] / (1.0 + w);
    }
    ldir[2] = cdir[2] * r;
    ldir[3] = cdir[3] * r * sin(theta);
}

/*! Test break condition.
 *
 *  \param pos    : pointer to position array.
 *  \return true  : radial position r < 0.0.
 *  \return false : position is valid.
 */
bool MetricEddFinkIn::breakCondition(const double* pos)
{
    bool br = false;

    if (pos[1] < 0.0) {
        br = true;
    }
    return br;
}

/*! Calculate right hand side of parallel transport and Jocobi equation.
 */
bool MetricEddFinkIn::calcDerivsSachsJacobi(const double y[], double dydx[])
{
    const double* u = &y[DEF_TG_IDX];
    const double* wSA1 = &y[DEF_SA1_IDX];
    const double* wSA2 = &y[DEF_SA2_IDX];
    const double* wJ1 = &y[DEF_JAC1_IDX];
    const double* wJ2 = &y[DEF_JAC2_IDX];
    const double* wJ1d = &y[DEF_DJ1_IDX];
    const double* wJ2d = &y[DEF_DJ2_IDX];

    double gd[4];
    contrChrisVecVec(y, u, u, gd);

    double zSA1[4], zSA2[4], zJ1a[4], zJ1b[4], zJ2a[4], zJ2b[4];
    contrChrisVecVec(y, u, wSA1, zSA1);
    contrChrisVecVec(y, u, wSA2, zSA2);

    contrChrisVecVec(y, u, wJ1d, zJ1a);
    contrChrisVecVec(y, u, wJ2d, zJ2a);

    contrChrDVecVecVec(y, u, u, wJ1, zJ1b);
    contrChrDVecVecVec(y, u, u, wJ2, zJ2b);

    for (int mu = 0; mu < 4; mu++) {
        dydx[mu] = y[DEF_TG_IDX + mu];
        dydx[DEF_TG_IDX + mu] = -gd[mu];

        dydx[DEF_SA1_IDX + mu] = -zSA1[mu];
        dydx[DEF_SA2_IDX + mu] = -zSA2[mu];
        dydx[DEF_JAC1_IDX + mu] = y[DEF_DJ1_IDX + mu];
        dydx[DEF_JAC2_IDX + mu] = y[DEF_DJ2_IDX + mu];

        dydx[DEF_DJ1_IDX + mu] = -2.0 * zJ1a[mu] - zJ1b[mu];
        dydx[DEF_DJ2_IDX + mu] = -2.0 * zJ2a[mu] - zJ2b[mu];
    }
    return true;
}

/*! Tests whether the constraint equation is fulfilled.
 *
 *  The constraint equation for lightlike and timelike geodesics reads:
 \verbatim
     sum = g_{\mu\nu} dot(x)^{\mu} dot(x)^{\nu} - kappa c^2 = 0.
 \endverbatim
 *  However, take care of the limited double precision.
 *  \param  y[]   : pointer to position and direction coordinates.
 *  \param  kappa : timelike (-1.0), lightlike (0.0).
 *  \return double : sum.
 */
double MetricEddFinkIn::testConstraint(const double y[], const double kappa)
{
    double r = y[1];
    double theta = y[2];
    double cm = 1.0 / mSpeedOfLight;

    // Scale the directions with the speed of light before doubling them !!
    double dv = y[4];
    double dr = y[5] * cm;
    double dth = y[6] * cm;
    double dph = y[7] * cm;

    double sum = -kappa;
    sum += -(1.0 - rs / r) * dv * dv + 2.0 * dv * dr + r * r * (dth * dth + sin(theta) * sin(theta) * dph * dph);
    return sum;
}

/*! Set parameter 'pName' to 'val'.
 *
 *  Set 'mass' parameter and adjust Schwarzschild radius  rs=2GM/c^2.
 */
bool MetricEddFinkIn::setParam(const char* pName, double val)
{
    if (Metric::setParam(pName, val)) {
        mMass = val;
        rs = 2.0 * mGravConstant * mMass / (mSpeedOfLight * mSpeedOfLight);
    }
    return true;
}

/*!
 *  \param units : type of physical constants.
 */
void MetricEddFinkIn::usePhysicalUnits(const enum_physical_constants units)
{
    Metric::usePhysicalUnits(units);
    rs = 2.0 * mGravConstant * mMass / (mSpeedOfLight * mSpeedOfLight);
}

/*!
 *  \param speed_of_light : value for speed of light.
 *  \param grav_const : value for gravitational constant.
 *  \param diel_perm : value for dielectric permittivity.
 */
void MetricEddFinkIn::setUnits(const double speed_of_light, const double grav_const, const double diel_perm)
{
    Metric::setUnits(speed_of_light, grav_const, diel_perm);
    rs = 2.0 * mGravConstant * mMass / (mSpeedOfLight * mSpeedOfLight);
}

/*! Generate report.
 */
bool MetricEddFinkIn::report(const vec4, const vec4, char*& text)
{
    std::stringstream ss;
    ss << "Report for Eddington-Finkelstein metric\n\tcoordinate : (v,r,theta,phi)\n";
    ss << "-----------------------------------------------------------------------\n";
    ss << "  physical units ................................. yes\n";
    ss.precision(DEF_FIXED_REPORT_PRECISION);
    ss.setf(std::ios::fixed);
    ss << "  Schwarzschild radius ........... r_s = 2GM/c^2 = " << rs << std::endl;

    text = new char[ss.str().length() + 2];
    return CopyString(ss.str().c_str(), text);
}

/*!  Get the critical angle for a freely falling observer
 *   with initial position r0 and current position r.
 *   See GRG 40,2185 Eq.(14)
 *
 * @param r0  initial observer position.
 * @param r   current observer position.
 * @param ca  critical angle in radians.
 * @return false  There is no critical angle.
 */
bool MetricEddFinkIn::GetCriticalAngle(const double r0, const double r, double& ca)
{
    if (r0 <= rs || r <= 0.0 || r0 < r) {
        return false;
    }

    double x0 = rs / r0;
    double x = rs / r;
    double p = 4.0 / 27.0;

    double a = x * x * sqrt(1.0 - x0) * sqrt(x - x0);
    double b = sqrt(p * (x * x * x - x * x + p));
    double c = x * x * x - x0 * x * x + p;

    if (r >= 1.5 * rs) {
        ca = acos((a + b) / c);
    }
    else {
        ca = acos((a - b) / c);
    }
    return true;
}

bool MetricEddFinkIn::GetCurrPosition(const double r0, const double tau, double& r)
{
    if (tau < 0.0) {
        return false;
    }
    if (tau == 0) {
        r = r0;
        return true;
    }

    struct EddFinkParams params = { rs, r0, tau };

    const gsl_root_fsolver_type* T = gsl_root_fsolver_brent;
    gsl_root_fsolver* s = gsl_root_fsolver_alloc(T);

    gsl_function F;
    F.function = tau_of_r;
    F.params = &params;
    int status = 0;
    int iter = 0, max_iter = 100;

    double r_lo = 0.001;
    double r_hi = r0;
    gsl_root_fsolver_set(s, &F, r_lo, r_hi);
    do {
        iter++;
        status = gsl_root_fsolver_iterate(s);
        r = gsl_root_fsolver_root(s);
        r_lo = gsl_root_fsolver_x_lower(s);
        r_hi = gsl_root_fsolver_x_upper(s);
        status = gsl_root_test_interval(r_lo, r_hi, 0, 1e-8);
    } while (status == GSL_CONTINUE && iter < max_iter);

    gsl_root_fsolver_free(s);
    return true;
}

/*!  Get freefall velocity at starting from r0.
 *   See GRG 40,2185 Eq.(11)
 *
 * @param r0  initial observer position.
 * @param r   current observer position.
 * @param beta  freefall velocity.
 * @return false  There is no freefall velocity.
 */
bool MetricEddFinkIn::GetFreefallVel(const double r0, const double r, double& beta)
{
    if (r0 <= rs || r <= 0.0 || r0 < r) {
        return false;
    }

    double x0 = rs / r0;
    double x = rs / r;
    beta = sqrt((x - x0) / (1.0 - x0));
    return true;
}

bool MetricEddFinkIn::GetTauCrash(const double r0, double& tau)
{
    tau = 0.5 * M_PI * rs * pow(r0 / rs, 1.5);
    return true;
}

// ********************************* protected methods *****************************
/*!
 */
void MetricEddFinkIn::setStandardValues()
{
    mInitPos[0] = 0.0;
    mInitPos[1] = 3.0 * rs;
    mInitPos[2] = M_PI_2;
    mInitPos[3] = 0.0;
    mInitDir[0] = 1.0;
    mInitDir[1] = 0.0;
    mInitDir[2] = 0.0;

    mCoordNames[0] = std::string("v");
    mCoordNames[1] = std::string("r");
    mCoordNames[2] = std::string("theta");
    mCoordNames[3] = std::string("phi");
}

/*! Contract Christoffel symbol with two vectors.
 *  \param y : pointer to full data
 *  \param u : pointer to first vector
 *  \param w : pointer to second vector
 *  \param z : pointer to contraction
 *  \param calc : calculate the Christoffels before using them
 *
 *  The Christoffel symbols do not have to be calculated before this function can be evaluated!
 */
void MetricEddFinkIn::contrChrisVecVec(const double y[], const double u[], const double w[], double* z, bool)
{
    double r = y[1];
    double theta = y[2];

    double c = mSpeedOfLight;

    double t1 = c * rs;
    double t2 = r * r;
    double t3 = 1 / t2;
    double t4 = t3 * u[0];
    double t9 = 1 / c * r;
    double t12 = sin(theta);
    double t13 = t12 * t12;
    double t18 = r - rs;
    double t22 = c * c;
    double t38 = u[3] * w[3];
    double t41 = 1 / r;
    double t44 = t41 * u[1];
    double t46 = cos(theta);
    double t53 = 1 / t12 * t46;

    z[0] = t1 * t4 * w[0] / 2.0 - t9 * u[2] * w[2] - t9 * t13 * u[3] * w[3];
    z[1] = t18 / t2 / r * t22 * rs * u[0] * w[0] / 2.0 - t1 * t3 * u[1] * w[0] / 2.0 - t1 * t4 * w[1] / 2.0
        - t18 * u[2] * w[2] - t18 * t13 * t38;
    z[2] = t41 * u[2] * w[1] + t44 * w[2] - t12 * t46 * t38;
    z[3] = t41 * u[3] * w[1] + t53 * u[3] * w[2] + t44 * w[3] + t53 * u[2] * w[3];
}

/*! Contract partially derived Christoffel symbols with three vectors.
 *  \param y : pointer to full data
 *  \param u : pointer to first vector
 *  \param v : pointer to second vector
 *  \param w : pointer to third vector
 *  \param z : pointer to contraction
 *  \param calc : calculate the Christoffels before using them
 *
 * The partial derivatives of the Christoffel symbols do not have to be calculated before this function can be
 * evaluated!
 */
void MetricEddFinkIn::contrChrDVecVecVec(
    const double y[], const double u[], const double v[], const double w[], double* z, bool)
{
    double r = y[1];
    double theta = y[2];

    double c = mSpeedOfLight;

    double t2 = r * r;
    double t5 = c * rs / t2 / r;
    double t9 = 1 / c;
    double t11 = v[2] * w[1];
    double t13 = sin(theta);
    double t14 = t13 * t13;
    double t16 = u[3] * v[3];
    double t21 = cos(theta);
    double t23 = v[3] * w[2];
    double t28 = c * c;
    double t34 = t2 * t2;
    double t50 = v[3] * w[1];
    double t59 = 1 / t2;
    double t61 = v[1] * w[1];
    double t63 = t59 * u[1];
    double t65 = t21 * t21;
    double t75 = (t65 + t14) / t14;

    z[0] = -t5 * u[0] * v[0] * w[1] - t9 * u[2] * t11 - t9 * t14 * t16 * w[1] - 2.0 * t9 * r * t13 * t21 * u[3] * t23;
    z[1] = -t28 * rs * (2.0 * r - 3.0 * rs) / t34 * u[0] * v[0] * w[1] / 2.0 + t5 * u[1] * v[0] * w[1]
        + t5 * u[0] * v[1] * w[1] - u[2] * v[2] * w[1] - t14 * u[3] * t50 - 2.0 * (r - rs) * t13 * t21 * t16 * w[2];
    z[2] = -t59 * u[2] * t61 - t63 * t11 + (-t65 + t14) * u[3] * t23;
    z[3] = -t59 * u[3] * t61 - t63 * t50 - t75 * u[3] * v[2] * w[2] - t75 * u[2] * v[3] * w[2];
}

} // end namespace m4d
