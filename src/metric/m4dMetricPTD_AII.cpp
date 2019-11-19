// -------------------------------------------------------------------------------
/*
   m4dMetricPTD_AII.cpp

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

#include "m4dMetricPTD_AII.h"

namespace m4d {

#define eps 1.0e-6

/*! Standard constructor for the metric.
 *
 * \param  b : parameter b of the AII metric
 */
MetricPTD_AII::MetricPTD_AII(double b)
{
    mMetricName = "Petrov_Type_D_AII_ES";
    setCoordType(enum_coordinate_cylinder);

    mPhysicalUnits = enum_physical_constants_geom;
    mSpeedOfLight = 1.0;
    mGravConstant = 1.0;

    Par_b = b;

    addParam("b", Par_b);

    setStandardValues();
}

/*! Standard destructor for the metric.
 *
 */

MetricPTD_AII::~MetricPTD_AII() {}

// *********************************** public methods ******************************
/*! Calculate the contravariant metric components at position 'pos'.
 *
 *  \param pos : pointer to position.
 */
bool MetricPTD_AII::calculateMetric(const double* pos)
{
    double r = pos[1];
    double z = pos[3];
    double b = Par_b;

    double t1 = 1;
    double t3 = 1 / z;
    double t6 = z * z;
    double t7 = sinh(r);
    double t8 = t7 * t7;
    g_compts[0][0] = -t1 * b * t3 + t1;
    g_compts[0][1] = 0.0;
    g_compts[0][2] = 0.0;
    g_compts[0][3] = 0.0;
    g_compts[1][0] = 0.0;
    g_compts[1][1] = t6;
    g_compts[1][2] = 0.0;
    g_compts[1][3] = 0.0;
    g_compts[2][0] = 0.0;
    g_compts[2][1] = 0.0;
    g_compts[2][2] = t6 * t8;
    g_compts[2][3] = 0.0;
    g_compts[3][0] = 0.0;
    g_compts[3][1] = 0.0;
    g_compts[3][2] = 0.0;
    g_compts[3][3] = 1 / (b * t3 - 1.0);

    return true;
}

/*! Calculate the Christoffel symbols of the second kind at position 'pos'.
 *
 *  \param pos : pointer to position.
 */
bool MetricPTD_AII::calculateChristoffels(const double* pos)
{
    double r = pos[1];
    double z = pos[3];
    double b = Par_b;

    double t1 = b - z;
    double t2 = z * z;
    double t6 = 1;
    double t10 = 1 / z;
    double t14 = t10 / t1 * b / 2.0;
    double t15 = sinh(r);
    double t17 = cosh(r);
    double t18 = 1 / t15 * t17;
    double t20 = t15 * t15;
    christoffel[0][0][0] = 0.0;
    christoffel[0][0][1] = 0.0;
    christoffel[0][0][2] = 0.0;
    christoffel[0][0][3] = -t1 / t2 / z * t6 * b / 2.0;
    christoffel[0][1][0] = 0.0;
    christoffel[0][1][1] = 0.0;
    christoffel[0][1][2] = 0.0;
    christoffel[0][1][3] = 0.0;
    christoffel[0][2][0] = 0.0;
    christoffel[0][2][1] = 0.0;
    christoffel[0][2][2] = 0.0;
    christoffel[0][2][3] = 0.0;
    christoffel[0][3][0] = -t14;
    christoffel[0][3][1] = 0.0;
    christoffel[0][3][2] = 0.0;
    christoffel[0][3][3] = 0.0;
    christoffel[1][0][0] = 0.0;
    christoffel[1][0][1] = 0.0;
    christoffel[1][0][2] = 0.0;
    christoffel[1][0][3] = 0.0;
    christoffel[1][1][0] = 0.0;
    christoffel[1][1][1] = 0.0;
    christoffel[1][1][2] = 0.0;
    christoffel[1][1][3] = -t1;
    christoffel[1][2][0] = 0.0;
    christoffel[1][2][1] = 0.0;
    christoffel[1][2][2] = t18;
    christoffel[1][2][3] = 0.0;
    christoffel[1][3][0] = 0.0;
    christoffel[1][3][1] = t10;
    christoffel[1][3][2] = 0.0;
    christoffel[1][3][3] = 0.0;
    christoffel[2][0][0] = 0.0;
    christoffel[2][0][1] = 0.0;
    christoffel[2][0][2] = 0.0;
    christoffel[2][0][3] = 0.0;
    christoffel[2][1][0] = 0.0;
    christoffel[2][1][1] = 0.0;
    christoffel[2][1][2] = t18;
    christoffel[2][1][3] = 0.0;
    christoffel[2][2][0] = 0.0;
    christoffel[2][2][1] = -t15 * t17;
    christoffel[2][2][2] = 0.0;
    christoffel[2][2][3] = -t1 * t20;
    christoffel[2][3][0] = 0.0;
    christoffel[2][3][1] = 0.0;
    christoffel[2][3][2] = t10;
    christoffel[2][3][3] = 0.0;
    christoffel[3][0][0] = -t14;
    christoffel[3][0][1] = 0.0;
    christoffel[3][0][2] = 0.0;
    christoffel[3][0][3] = 0.0;
    christoffel[3][1][0] = 0.0;
    christoffel[3][1][1] = t10;
    christoffel[3][1][2] = 0.0;
    christoffel[3][1][3] = 0.0;
    christoffel[3][2][0] = 0.0;
    christoffel[3][2][1] = 0.0;
    christoffel[3][2][2] = t10;
    christoffel[3][2][3] = 0.0;
    christoffel[3][3][0] = 0.0;
    christoffel[3][3][1] = 0.0;
    christoffel[3][3][2] = 0.0;
    christoffel[3][3][3] = t14;

    return true;
}

/*! Calculate Jacobi matrix.
 *
 *  \param pos : pointer to position.
 */
bool MetricPTD_AII::calculateChrisD(const double* pos)
{
    double r = pos[1];
    double z = pos[3];
    double b = Par_b;

    double t1 = 2.0 * z;
    double t5 = z * z;
    double t6 = t5 * t5;
    double t12 = -b + z;
    double t13 = t12 * t12;
    double t15 = 1 / t5;
    double t18 = b * (t1 - b) / t13 * t15 / 2.0;
    double t19 = cosh(r);
    double t20 = t19 * t19;
    double t21 = sinh(r);
    double t22 = t21 * t21;
    double t25 = (-t20 + t22) / t22;
    chrisD[0][0][0][0] = 0.0;
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
    chrisD[0][0][3][3] = -b * (t1 - 3.0 * b) / t6 / 2.0;
    chrisD[0][1][0][0] = 0.0;
    chrisD[0][1][0][1] = 0.0;
    chrisD[0][1][0][2] = 0.0;
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
    chrisD[0][3][0][3] = -t18;
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
    chrisD[1][1][3][3] = 1.0;
    chrisD[1][2][0][0] = 0.0;
    chrisD[1][2][0][1] = 0.0;
    chrisD[1][2][0][2] = 0.0;
    chrisD[1][2][0][3] = 0.0;
    chrisD[1][2][1][0] = 0.0;
    chrisD[1][2][1][1] = 0.0;
    chrisD[1][2][1][2] = 0.0;
    chrisD[1][2][1][3] = 0.0;
    chrisD[1][2][2][0] = 0.0;
    chrisD[1][2][2][1] = t25;
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
    chrisD[1][3][1][3] = -t15;
    chrisD[1][3][2][0] = 0.0;
    chrisD[1][3][2][1] = 0.0;
    chrisD[1][3][2][2] = 0.0;
    chrisD[1][3][2][3] = 0.0;
    chrisD[1][3][3][0] = 0.0;
    chrisD[1][3][3][1] = 0.0;
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
    chrisD[2][1][2][1] = t25;
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
    chrisD[2][2][1][1] = -t20 - t22;
    chrisD[2][2][1][2] = 0.0;
    chrisD[2][2][1][3] = 0.0;
    chrisD[2][2][2][0] = 0.0;
    chrisD[2][2][2][1] = 0.0;
    chrisD[2][2][2][2] = 0.0;
    chrisD[2][2][2][3] = 0.0;
    chrisD[2][2][3][0] = 0.0;
    chrisD[2][2][3][1] = 2.0 * t12 * t21 * t19;
    chrisD[2][2][3][2] = 0.0;
    chrisD[2][2][3][3] = t22;
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
    chrisD[2][3][2][3] = -t15;
    chrisD[2][3][3][0] = 0.0;
    chrisD[2][3][3][1] = 0.0;
    chrisD[2][3][3][2] = 0.0;
    chrisD[2][3][3][3] = 0.0;
    chrisD[3][0][0][0] = 0.0;
    chrisD[3][0][0][1] = 0.0;
    chrisD[3][0][0][2] = 0.0;
    chrisD[3][0][0][3] = -t18;
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
    chrisD[3][1][1][3] = -t15;
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
    chrisD[3][2][2][3] = -t15;
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
    chrisD[3][3][3][3] = t18;
    return true;
}

/*! Transform local 4-direction to coordinate 4-direction.
 *
 *  \param  pos  :  pointer to position array.
 *  \param  ldir :  pointer to local direction array.
 *  \param  dir  :  pointer to calculated coordinate direction array.
 *  \param  type :  type of tetrad.
 */
void MetricPTD_AII::localToCoord(const double* pos, const double* ldir, double* dir, enum_nat_tetrad_type)
{
    double r = pos[1];
    double z = pos[3];
    double w = sqrt(Par_b / pos[3] - 1);

    dir[0] = ldir[0] / w;
    dir[1] = ldir[1] / z;
    dir[2] = ldir[2] / (z * sinh(r));
    dir[3] = ldir[3] * w;
}

/*! Transform coordinate 4-direction to local 4-direction.
 *
 *  \param  pos  :  pointer to position array.
 *  \param  cdir :  pointer to coordinate direction.
 *  \param  ldir :  pointer to calculated local direction array.
 *  \param  type :  type of tetrad.
 */
void MetricPTD_AII::coordToLocal(const double* pos, const double* cdir, double* ldir, enum_nat_tetrad_type)
{
    double r = pos[1];
    double z = pos[3];
    double w = sqrt(Par_b / pos[3] - 1);

    ldir[0] = cdir[0] * w;
    ldir[1] = cdir[1] * z;
    ldir[2] = cdir[2] * (z * sinh(r));
    ldir[3] = cdir[3] / w;
}

/*! Tests break condition.
 *  \param pos  :  position.
 *  \return true  : position z < 0.0 or z reaches b
 *  \return false : position is valid.
 */
bool MetricPTD_AII::breakCondition(const double* pos)
{
    bool br = false;

    if ((pos[3] < 0.0) || (pos[3] * pos[3] >= (1.0 - eps) * Par_b * Par_b)) {
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
double MetricPTD_AII::testConstraint(const double y[], const double kappa)
{
    double r = y[1];
    double z = y[3];
    double cm = 1.0 / mSpeedOfLight;

    // Scale the directions with the speed of light before doubling them !!
    double dt = y[4];
    double dr = y[5] * cm;
    double dph = y[6] * cm;
    double dz = y[7] * cm;

    double sum = -kappa;
    sum += -(Par_b / z - 1.0) * dt * dt + dz * dz / (Par_b / z - 1.0)
        + z * z * (dr * dr + sinh(r) * sinh(r) * dph * dph);
    return sum;
}

/*! Set parameter 'pName' to 'val'.
 *
 */
bool MetricPTD_AII::setParam(const char* pName, double val)
{
    Metric::setParam(pName, val);
    if (strcmp(pName, "b") == 0) {
        Par_b = val;
    }
    return true;
}

/*! Generate report.
 */
bool MetricPTD_AII::report(const vec4 pos, const vec4 cdir, char*& text)
{
    std::stringstream ss;
    ss << "Report for AII metric\n\tcoordinates : (t,r,phi,z)\n";
    ss << "---------------------------------------------------------------\n";
    ss << "Coordinate Ranges:" << std::endl;
    ss << "        r: (0 <= r)" << std::endl;
    ss << "      phi: (phi mod 2*pi)" << std::endl;
    ss << "        z: (0 < z < b)" << std::endl;
    ss << "---------------------------------------------------------------\n";
    ss << "Parameter Range:" << std::endl;
    ss << "        b: (0 < b)" << std::endl;
    ss << "---------------------------------------------------------------\n";
    ss.precision(DEF_FIXED_REPORT_PRECISION);
    ss.setf(std::ios::fixed);

    double sinp2 = sinh(pos[1]);
    double p14 = pos[3] * pos[3] * pos[3] * pos[3];
    double b = Par_b;

    double C0 = cdir[0] * cdir[0] * (pos[3] - b) * (pos[3] - b) / pos[3] / pos[3];
    double C2 = cdir[2] * cdir[2] * p14 * sinp2 * sinp2 * sinp2 * sinp2;
    double K = cdir[1] * cdir[1] * p14 + cdir[2] * cdir[2] * p14 * sinp2 * sinp2;
    //  double m0 = -K/pos[3]/pos[3]  + cdir[3]*cdir[3]*pos[3]/(pos[3] - b) -  cdir[0]*cdir[0]*(pos[3] - b)/pos[3];
    ss << "Constants of motion:" << std::endl;
    ss << "  C_0^2 = " << C0 << std::endl;
    ss << "  C_2^2 = " << C2 << std::endl;
    ss << "      K = " << K << std::endl;
    //  ss << "  m_0^2 = " << m0 << std::endl;
    ss << "---------------------------------------------------------------\n";

    text = new char[ss.str().length() + 2];
    return CopyString(ss.str().c_str(), text);
}

// *************************** specific  public methods ****************************
// None
// ********************************* protected methods *****************************
/*!
 */
void MetricPTD_AII::setStandardValues()
{
    mInitPos[0] = 0.0;
    mInitPos[1] = 1.0;
    mInitPos[2] = 0.0;
    mInitPos[3] = 0.6;
    mInitDir[0] = 1.0;
    mInitDir[1] = 0.0;
    mInitDir[2] = 0.0;

    mCoordNames[0] = std::string("t");
    mCoordNames[1] = std::string("r");
    mCoordNames[2] = std::string("phi");
    mCoordNames[3] = std::string("z");
}

} // end namespace m4d
