// -------------------------------------------------------------------------------
/*
   m4dMetricPT-DBI.cpp

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

#include "m4dMetricPTD_BI.h"

namespace m4d {

#define eps 1.0e-6


/*! Standard constructor for the metric.
 *
 * \param  b : parameter b of the BI metric
 */
MetricPTD_BI::MetricPTD_BI(double b) {
    mMetricName  = "Petrov_Type_D_BI_ES";
    setCoordType(enum_coordinate_spherical);

    mPhysicalUnits = enum_physical_constants_geom;
    mSpeedOfLight = 1.0;
    mGravConstant = 1.0;

    mDrawTypes.push_back(enum_draw_effpoti);

    Par_b = b;

    addParam("b", Par_b);

    setStandardValues();

// mLocTeds.push_back(enum_nat_tetrad_static);


}

/*! Standard destructor for the metric.
 *
 */
MetricPTD_BI::~MetricPTD_BI() {

}


// *********************************** public methods ******************************
/*! Calculate the contravariant metric components at position 'pos'.
 *
 *  \param pos : pointer to position.
 */
bool MetricPTD_BI::calculateMetric(const double* pos) {
    double r     = pos[1];
    double theta = pos[2];
    double b = Par_b;
    double t1 = r * r;
    double  t2 = sin(theta);
    double  t3 = t2 * t2;
    double  t5 = -r + b;
    g_compts[0][0] = -t1 * t3;
    g_compts[0][1] = 0.0;
    g_compts[0][2] = 0.0;
    g_compts[0][3] = 0.0;
    g_compts[1][0] = 0.0;
    g_compts[1][1] = -1 / t5 * r;
    g_compts[1][2] = 0.0;
    g_compts[1][3] = 0.0;
    g_compts[2][0] = 0.0;
    g_compts[2][1] = 0.0;
    g_compts[2][2] = t1;
    g_compts[2][3] = 0.0;
    g_compts[3][0] = 0.0;
    g_compts[3][1] = 0.0;
    g_compts[3][2] = 0.0;
    g_compts[3][3] = -1 / r * t5;

    return true;
}

/*! Calculate the Christoffel symbols of the second kind at position 'pos'.
 *
 *  \param pos : pointer to position.
 */
bool MetricPTD_BI::calculateChristoffels(const double* pos) {
    double r     = pos[1];
    double theta = pos[2];
    double b = Par_b;
    double t1 = -r + b;
    double  t2 = sin(theta);
    double  t3 = t2 * t2;
    double  t5 = cos(theta);
    double  t7 = 1 / r;
    double  t9 = 1 / t2 * t5;
    double  t13 = 1 / t1 * t7 * b / 2.0;
    double  t14 = r * r;
    christoffel[0][0][0] = 0.0;
    christoffel[0][0][1] = -t1 * t3;
    christoffel[0][0][2] = t2 * t5;
    christoffel[0][0][3] = 0.0;
    christoffel[0][1][0] = t7;
    christoffel[0][1][1] = 0.0;
    christoffel[0][1][2] = 0.0;
    christoffel[0][1][3] = 0.0;
    christoffel[0][2][0] = t9;
    christoffel[0][2][1] = 0.0;
    christoffel[0][2][2] = 0.0;
    christoffel[0][2][3] = 0.0;
    christoffel[0][3][0] = 0.0;
    christoffel[0][3][1] = 0.0;
    christoffel[0][3][2] = 0.0;
    christoffel[0][3][3] = 0.0;
    christoffel[1][0][0] = t7;
    christoffel[1][0][1] = 0.0;
    christoffel[1][0][2] = 0.0;
    christoffel[1][0][3] = 0.0;
    christoffel[1][1][0] = 0.0;
    christoffel[1][1][1] = t13;
    christoffel[1][1][2] = 0.0;
    christoffel[1][1][3] = 0.0;
    christoffel[1][2][0] = 0.0;
    christoffel[1][2][1] = 0.0;
    christoffel[1][2][2] = t7;
    christoffel[1][2][3] = 0.0;
    christoffel[1][3][0] = 0.0;
    christoffel[1][3][1] = 0.0;
    christoffel[1][3][2] = 0.0;
    christoffel[1][3][3] = -t13;
    christoffel[2][0][0] = t9;
    christoffel[2][0][1] = 0.0;
    christoffel[2][0][2] = 0.0;
    christoffel[2][0][3] = 0.0;
    christoffel[2][1][0] = 0.0;
    christoffel[2][1][1] = 0.0;
    christoffel[2][1][2] = t7;
    christoffel[2][1][3] = 0.0;
    christoffel[2][2][0] = 0.0;
    christoffel[2][2][1] = t1;
    christoffel[2][2][2] = 0.0;
    christoffel[2][2][3] = 0.0;
    christoffel[2][3][0] = 0.0;
    christoffel[2][3][1] = 0.0;
    christoffel[2][3][2] = 0.0;
    christoffel[2][3][3] = 0.0;
    christoffel[3][0][0] = 0.0;
    christoffel[3][0][1] = 0.0;
    christoffel[3][0][2] = 0.0;
    christoffel[3][0][3] = 0.0;
    christoffel[3][1][0] = 0.0;
    christoffel[3][1][1] = 0.0;
    christoffel[3][1][2] = 0.0;
    christoffel[3][1][3] = -t13;
    christoffel[3][2][0] = 0.0;
    christoffel[3][2][1] = 0.0;
    christoffel[3][2][2] = 0.0;
    christoffel[3][2][3] = 0.0;
    christoffel[3][3][0] = 0.0;
    christoffel[3][3][1] = t1 / t14 / r * b / 2.0;
    christoffel[3][3][2] = 0.0;
    christoffel[3][3][3] = 0.0;

    return true;
}

/*! Calculate Jacobi matrix.
 *
 *  \param pos : pointer to position.
 */
bool
MetricPTD_BI::calculateChrisD(const double* pos) {
    double r     = pos[1];
    double theta = pos[2];
    double b = Par_b;
    double t1 = sin(theta);
    double  t2 = t1 * t1;
    double  t3 = -r + b;
    double  t5 = cos(theta);
    double  t8 = t5 * t5;
    double  t10 = r * r;
    double  t11 = 1 / t10;
    double  t14 = (t8 + t2) / t2;
    double  t15 = 2.0 * r;
    double  t18 = t3 * t3;
    double  t22 = b * (-t15 + b) / t18 * t11 / 2.0;
    double  t26 = t10 * t10;
    chrisD[0][0][0][0] = 0.0;
    chrisD[0][0][0][1] = 0.0;
    chrisD[0][0][0][2] = 0.0;
    chrisD[0][0][0][3] = 0.0;
    chrisD[0][0][1][0] = 0.0;
    chrisD[0][0][1][1] = t2;
    chrisD[0][0][1][2] = -2.0 * t3 * t1 * t5;
    chrisD[0][0][1][3] = 0.0;
    chrisD[0][0][2][0] = 0.0;
    chrisD[0][0][2][1] = 0.0;
    chrisD[0][0][2][2] = t8 - t2;
    chrisD[0][0][2][3] = 0.0;
    chrisD[0][0][3][0] = 0.0;
    chrisD[0][0][3][1] = 0.0;
    chrisD[0][0][3][2] = 0.0;
    chrisD[0][0][3][3] = 0.0;
    chrisD[0][1][0][0] = 0.0;
    chrisD[0][1][0][1] = -t11;
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
    chrisD[0][2][0][2] = -t14;
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
    chrisD[1][0][0][1] = -t11;
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
    chrisD[1][1][1][1] = -t22;
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
    chrisD[1][2][2][1] = -t11;
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
    chrisD[1][3][3][1] = t22;
    chrisD[1][3][3][2] = 0.0;
    chrisD[1][3][3][3] = 0.0;
    chrisD[2][0][0][0] = 0.0;
    chrisD[2][0][0][1] = 0.0;
    chrisD[2][0][0][2] = -t14;
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
    chrisD[2][1][2][1] = -t11;
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
    chrisD[2][3][3][2] = 0.0;
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
    chrisD[3][1][3][1] = t22;
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
    chrisD[3][3][1][1] = -b * (-t15 + 3.0 * b) / t26 / 2.0;
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

/*! Transform local 4-direction to coordinate 4-direction.
 *
 *  \param  pos  :  pointer to position array.
 *  \param  ldir :  pointer to local direction array.
 *  \param  dir  :  pointer to calculated coordinate direction array.
 *  \param  type :  type of tetrad.
 */
void MetricPTD_BI::localToCoord(const double* pos, const double* ldir, double* dir,
                                enum_nat_tetrad_type) {
    double r     = pos[1];
    double theta = pos[2];
    double w = sqrt(1.0 - Par_b / r);

    dir[0] = ldir[0] / (r * sin(theta));
    dir[1] = ldir[1] * w;
    dir[2] = ldir[2] / r;
    dir[3] = ldir[3] / w;
}

/*! Transform coordinate 4-direction to local 4-direction.
 *
 *  \param  pos  :  pointer to position array.
 *  \param  cdir :  pointer to coordinate direction.
 *  \param  ldir :  pointer to calculated local direction array.
 *  \param  type :  type of tetrad.
 */
void MetricPTD_BI::coordToLocal(const double* pos, const double* cdir, double* ldir,
                                enum_nat_tetrad_type) {
    double r     = pos[1];
    double theta = pos[2];
    double w = sqrt(1.0 - Par_b / r);

    ldir[0] = cdir[0] * (r * sin(theta));
    ldir[1] = cdir[1] / w;
    ldir[2] = cdir[2] * r;
    ldir[3] = cdir[3] * w;
}


/*! Tests break condition
 *  \param pos  :  position.
 *  \return true  : radial position r < 0.0 or r reaches b
 *  \return false : position is valid.
 */
bool MetricPTD_BI::breakCondition(const double* pos) {
    bool br = false;

    if ((pos[1] < 0.0) || (pos[1]*pos[1] <= (1.0 + eps)*Par_b * Par_b)) {
        br = true;
    }
    return br;
}

/*! Calculates the constants of Motion.
 *
 *  \param pos : initial position.
 *  \param cdir : initial four-direction.
 */
void MetricPTD_BI::calcConstantsOfMotion(const vec4 pos, const vec4 cdir) {
    double sinp2 = sin(pos[2]);
    double p14 = pos[1] * pos[1] * pos[1] * pos[1];
    double b = Par_b;

    C0 = cdir[3] * cdir[3] * (pos[1] - b) * (pos[1] - b) / pos[1] / pos[1];
    C2 = cdir[0] * cdir[0] * p14 * sinp2 * sinp2 * sinp2 * sinp2;
    K  = (cdir[2] * cdir[2] - cdir[0] * cdir[0] * sinp2 * sinp2) * p14;
    m0 = -K / pos[1] / pos[1]  - cdir[1] * cdir[1] * pos[1] / (pos[1] - b) - cdir[3] * cdir[3] * (pos[1] - b) / pos[1];
}

/*! Effective potential.
 *  \param pos : initial position.
 *  \param cdir : initial four-direction.
 *  \param type : geodesic type.
 *  \param x : abscissa value.
 *  \param val : reference to effective potential value.
 *  \return true : effective potential exists at x.
 */
bool MetricPTD_BI::effPotentialValue(const vec4 pos, const vec4 cdir , enum_geodesic_type , const double x, double &val) {
    /*
      double kappa = 0.0;
    if (type==enum_geodesic_timelike)
      kappa = -mSign;
    */
    if (x <= 0.0) {
        return false;
    }

    double r = x;
    double b = Par_b;
    calcConstantsOfMotion(pos, cdir);

    val = (r - b) / r / r / r * K + (r - b) / r * m0 + C0;
    return true;
}

/*! Total energy.
 *  \param pos : initial position.
 *  \param cdir : initial four-direction.
 *  \param x : abscissa value.
 *  \param val : reference to total energy value.
 *  \return true : effective potential exists at x.
 */
bool MetricPTD_BI::totEnergy(const vec4 , const vec4 , const double , double &val) {
    val = 0.0;
    return true;
}


/*! Calculates the root of the effective Potential left to r_0 for given constants of motion and m0=0 (lightlike)
 *
 *  \param C02 : constant of motion C0^2
 *  \param K   : constant of motion K
 *  \param r0  : start condition r_0
 *  \return -1 if there is no root left to r_0.
 *
 */
double MetricPTD_BI::calculateVeffRoot(double C02, double K, double r0) {
    double r = 0;
    double b = Par_b;

    if (C02 == 0) { // only 1st-Order Polynom
        r = b;
    } else {
        double p = -K / C02;
        double q = -b * K / C02;
        double d = q * q / 4.0 - p * p * p / 27.0;
        double z1 = 0;
        double z2 = 0;
        double z3 = 0;

        if (d < 0) {
            double u = acos(-q * 0.5 * sqrt(27.0 / p / p / p));
            z1 =  sqrt(4.0 / 3.0 * p) * cos((u) / 3.0);
            z2 =  -sqrt(4.0 / 3.0 * p) * cos((u + M_PI) / 3.0);
            z3 =  -sqrt(4.0 / 3.0 * p) * cos((u - M_PI) / 3.0);
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
            z1 = pow(-q * 0.5 + sqrt(d), 1 / 3.0) - pow(q * 0.5 + sqrt(d), 1 / 3.0);
            z3 = z2 = z1;
        }
        r = z1;
        if ((z2 > r0) && (fabs(r0 - z2) < fabs(r0 - r))) {
            r = z2;
        }
        if ((z3 > r0) && (fabs(r0 - z3) < fabs(r0 - r))) {
            r = z3;
        }
    }

    if (r < r0) {
        return -1;
    }
    return r;
}

/*! Set parameter 'pName' to 'val'.
 *
 *
 */
bool MetricPTD_BI::setParam(const char* pName, double val) {
    Metric::setParam(pName, val);
    if (strcmp(pName,"b") == 0) {
        Par_b = val;
    }
    return true;
}


/*! Generate report.
 */
bool MetricPTD_BI::report(const vec4 pos, const vec4 cdir, std::string &text) {
    std::stringstream ss;
    ss << "Report for BI metric\n\tcoordinates : (t,r,theta,phi)\n";
    ss << "---------------------------------------------------------------\n";
    ss << "Coordinate Ranges:" << std::endl;
    ss << "        r: (0 < b < r) or (b < 0 < r)" << std::endl;
    ss << "    theta: (0 < theta < pi)" << std::endl;
    ss << "      phi: arbitrary" << std::endl;
    ss << "---------------------------------------------------------------\n";
    ss << "Parameter Range:" << std::endl;
    ss << "        b: arbitrary" << std::endl;
    ss << "---------------------------------------------------------------\n";
    ss.precision(DEF_FIXED_REPORT_PRECISION);
    ss.setf(std::ios::fixed);

    //double b= Par_b;
    calcConstantsOfMotion(pos, cdir);
    ss << "Constants of motion:" << std::endl;
    ss << "  C_0^2 = " << C0 << std::endl;
    ss << "  C_2^2 = " << C2 << std::endl;
    ss << "      K = " << K << std::endl;
//  ss << "  m_0^2 = " << m0 << std::endl;
    ss << "---------------------------------------------------------------\n";
    double maxrad = calculateVeffRoot(C0, K, pos[1]);
    if (maxrad >= 0) {
        ss << "  maximum distance for light from origin....r_max= " << maxrad << std::endl;
    } else {
        ss << "  maximum distance from origin..............r_max= infinity" << std::endl;
    }

// int r= 30;
// double xi = asin(sqrt((r-b)/r/r/r*pos[1]*pos[1]*pos[1]/(pos[1]-b)))*180/M_PI;
// ss << "  maximum angle xi for light reaching radius r="<< r << ": " << xi << "." ;
    text = ss.str();
    return true;
}

// *************************** specific  public methods ****************************
//None
// ********************************* protected methods *****************************
/*!
 */
void MetricPTD_BI::setStandardValues() {
    mInitPos[0] = 0.0;
    mInitPos[1] = 6.0 * Par_b;
    mInitPos[2] = M_PI_2;
    mInitPos[3] = 0.0;
    mInitDir[0] = 1.0;
    mInitDir[1] = 0.0;
    mInitDir[2] = 0.0;

    mCoordNames[0] = std::string("t");
    mCoordNames[1] = std::string("r");
    mCoordNames[2] = std::string("theta");
    mCoordNames[3] = std::string("phi");
}


} // end namespace m4d

