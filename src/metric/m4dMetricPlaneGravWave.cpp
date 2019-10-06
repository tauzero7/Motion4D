// -------------------------------------------------------------------------------
/*
 m4dMetricPlaneGravWave.cpp

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

#include "m4dMetricPlaneGravWave.h"

namespace m4d {

#define eps 1.0e-6

/*! Standard constructor for the PlaneGravWave metric.
 */
MetricPlaneGravWave::MetricPlaneGravWave(double longExt, double degree)
{
    mMetricName = "PlaneGravWave";
    setCoordType(enum_coordinate_cartesian);

    mPhysicalUnits = enum_physical_constants_geom;
    mSpeedOfLight = 1.0;
    mGravConstant = 1.0;

    uData = dmData = nullptr;

    addParam("long_ext", longExt);
    setParam("long_ext", longExt);
    addParam("degree", degree);
    setParam("degree", degree);

    mLongExt = longExt;
    mDegree = degree;

    dataCalculated = false;

    mDrawTypes.push_back(enum_draw_twoplusone);

    setStandardValues();
}

MetricPlaneGravWave::~MetricPlaneGravWave()
{
    if (uData != nullptr) {
        delete[] uData;
    }
    uData = nullptr;

    if (dmData != nullptr) {
        delete[] dmData;
    }
    dmData = nullptr;
}

// *********************************** public methods ******************************
/*! Calculate the contravariant metric components at position 'pos'.
 *
 *  \param pos : pointer to position.
 */
bool MetricPlaneGravWave::calculateMetric(const double* pos)
{
    double c = mSpeedOfLight;
    double p = getValP(pos);
    double q = getValQ(pos);

    g_compts[0][0] = -c * c;
    g_compts[0][1] = 0;
    g_compts[0][2] = 0;
    g_compts[0][3] = 0;
    g_compts[1][0] = 0;
    g_compts[1][1] = 1;
    g_compts[1][2] = 0;
    g_compts[1][3] = 0;
    g_compts[2][0] = 0;
    g_compts[2][1] = 0;
    g_compts[2][2] = p * p;
    g_compts[2][3] = 0;
    g_compts[3][0] = 0;
    g_compts[3][1] = 0;
    g_compts[3][2] = 0;
    g_compts[3][3] = q * q;

    return true;
}

/*! Calculate the Christoffel symbols of the second kind at position 'pos'.
 *
 *  \param pos : pointer to position.
 */
bool MetricPlaneGravWave::calculateChristoffels(const double* pos)
{
    double c = mSpeedOfLight;
    double p = getValP(pos);
    double q = getValQ(pos);
    double dp = getValDP(pos); // diff(p,u)
    double dq = getValDQ(pos); // diff(q,u)

    double t1 = p;
    double t2 = 1.0 / t1;
    double t3 = c * dp;
    double t4 = (t2 * t3);
    double t5 = q;
    double t6 = 1.0 / t5;
    double t7 = c * dq;
    double t8 = (t6 * t7);
    double t9 = -dp;
    double t10 = (t2 * t9);
    double t11 = -dq;
    double t12 = (t6 * t11);
    double t13 = c * c;
    double t14 = 1 / t13;

    christoffel[0][0][0] = 0;
    christoffel[0][0][1] = 0;
    christoffel[0][0][2] = 0;
    christoffel[0][0][3] = 0;
    christoffel[0][1][0] = 0;
    christoffel[0][1][1] = 0;
    christoffel[0][1][2] = 0;
    christoffel[0][1][3] = 0;
    christoffel[0][2][0] = 0;
    christoffel[0][2][1] = 0;
    christoffel[0][2][2] = t4;
    christoffel[0][2][3] = 0;
    christoffel[0][3][0] = 0;
    christoffel[0][3][1] = 0;
    christoffel[0][3][2] = 0;
    christoffel[0][3][3] = t8;
    christoffel[1][0][0] = 0;
    christoffel[1][0][1] = 0;
    christoffel[1][0][2] = 0;
    christoffel[1][0][3] = 0;
    christoffel[1][1][0] = 0;
    christoffel[1][1][1] = 0;
    christoffel[1][1][2] = 0;
    christoffel[1][1][3] = 0;
    christoffel[1][2][0] = 0;
    christoffel[1][2][1] = 0;
    christoffel[1][2][2] = t10;
    christoffel[1][2][3] = 0;
    christoffel[1][3][0] = 0;
    christoffel[1][3][1] = 0;
    christoffel[1][3][2] = 0;
    christoffel[1][3][3] = t12;
    christoffel[2][0][0] = 0;
    christoffel[2][0][1] = 0;
    christoffel[2][0][2] = t4;
    christoffel[2][0][3] = 0;
    christoffel[2][1][0] = 0;
    christoffel[2][1][1] = 0;
    christoffel[2][1][2] = t10;
    christoffel[2][1][3] = 0;
    christoffel[2][2][0] = t14 * t1 * t3;
    christoffel[2][2][1] = -t1 * t9;
    christoffel[2][2][2] = 0;
    christoffel[2][2][3] = 0;
    christoffel[2][3][0] = 0;
    christoffel[2][3][1] = 0;
    christoffel[2][3][2] = 0;
    christoffel[2][3][3] = 0;
    christoffel[3][0][0] = 0;
    christoffel[3][0][1] = 0;
    christoffel[3][0][2] = 0;
    christoffel[3][0][3] = t8;
    christoffel[3][1][0] = 0;
    christoffel[3][1][1] = 0;
    christoffel[3][1][2] = 0;
    christoffel[3][1][3] = t12;
    christoffel[3][2][0] = 0;
    christoffel[3][2][1] = 0;
    christoffel[3][2][2] = 0;
    christoffel[3][2][3] = 0;
    christoffel[3][3][0] = t14 * t5 * t7;
    christoffel[3][3][1] = -t5 * t11;
    christoffel[3][3][2] = 0;
    christoffel[3][3][3] = 0;

    return true;
}

/*! Calculate Jacobi matrix.
 *
 *  \param pos : pointer to position.
 */
bool MetricPlaneGravWave::calculateChrisD(const double* pos)
{
    //	double c = mSpeedOfLight;

    double p = getValP(pos);
    double q = getValQ(pos);
    double dp = getValDP(pos);
    double dq = getValDQ(pos);
    double ddp = getValDDP(pos);
    double ddq = getValDDQ(pos);

    double t1 = dp * dp;
    // double t2 = t - x;
    double t3 = p;
    double t4 = ddp * t3;
    double t6 = t3 * t3;
    double t8 = ((-t1 + t4) / t6);
    double t9 = dq * dq;
    double t10 = q;
    double t11 = ddq * t10;
    double t13 = t10 * t10;
    double t15 = ((-t9 + t11) / t13);
    double t16 = (t1 + t4);
    double t17 = (t9 + t11);

    chrisD[0][0][0][0] = 0;
    chrisD[0][0][0][1] = 0;
    chrisD[0][0][0][2] = 0;
    chrisD[0][0][0][3] = 0;
    chrisD[0][0][1][0] = 0;
    chrisD[0][0][1][1] = 0;
    chrisD[0][0][1][2] = 0;
    chrisD[0][0][1][3] = 0;
    chrisD[0][0][2][0] = 0;
    chrisD[0][0][2][1] = 0;
    chrisD[0][0][2][2] = 0;
    chrisD[0][0][2][3] = 0;
    chrisD[0][0][3][0] = 0;
    chrisD[0][0][3][1] = 0;
    chrisD[0][0][3][2] = 0;
    chrisD[0][0][3][3] = 0;
    chrisD[0][1][0][0] = 0;
    chrisD[0][1][0][1] = 0;
    chrisD[0][1][0][2] = 0;
    chrisD[0][1][0][3] = 0;
    chrisD[0][1][1][0] = 0;
    chrisD[0][1][1][1] = 0;
    chrisD[0][1][1][2] = 0;
    chrisD[0][1][1][3] = 0;
    chrisD[0][1][2][0] = 0;
    chrisD[0][1][2][1] = 0;
    chrisD[0][1][2][2] = 0;
    chrisD[0][1][2][3] = 0;
    chrisD[0][1][3][0] = 0;
    chrisD[0][1][3][1] = 0;
    chrisD[0][1][3][2] = 0;
    chrisD[0][1][3][3] = 0;
    chrisD[0][2][0][0] = 0;
    chrisD[0][2][0][1] = 0;
    chrisD[0][2][0][2] = 0;
    chrisD[0][2][0][3] = 0;
    chrisD[0][2][1][0] = 0;
    chrisD[0][2][1][1] = 0;
    chrisD[0][2][1][2] = 0;
    chrisD[0][2][1][3] = 0;
    chrisD[0][2][2][0] = t8;
    chrisD[0][2][2][1] = -t8;
    chrisD[0][2][2][2] = 0;
    chrisD[0][2][2][3] = 0;
    chrisD[0][2][3][0] = 0;
    chrisD[0][2][3][1] = 0;
    chrisD[0][2][3][2] = 0;
    chrisD[0][2][3][3] = 0;
    chrisD[0][3][0][0] = 0;
    chrisD[0][3][0][1] = 0;
    chrisD[0][3][0][2] = 0;
    chrisD[0][3][0][3] = 0;
    chrisD[0][3][1][0] = 0;
    chrisD[0][3][1][1] = 0;
    chrisD[0][3][1][2] = 0;
    chrisD[0][3][1][3] = 0;
    chrisD[0][3][2][0] = 0;
    chrisD[0][3][2][1] = 0;
    chrisD[0][3][2][2] = 0;
    chrisD[0][3][2][3] = 0;
    chrisD[0][3][3][0] = t15;
    chrisD[0][3][3][1] = -t15;
    chrisD[0][3][3][2] = 0;
    chrisD[0][3][3][3] = 0;
    chrisD[1][0][0][0] = 0;
    chrisD[1][0][0][1] = 0;
    chrisD[1][0][0][2] = 0;
    chrisD[1][0][0][3] = 0;
    chrisD[1][0][1][0] = 0;
    chrisD[1][0][1][1] = 0;
    chrisD[1][0][1][2] = 0;
    chrisD[1][0][1][3] = 0;
    chrisD[1][0][2][0] = 0;
    chrisD[1][0][2][1] = 0;
    chrisD[1][0][2][2] = 0;
    chrisD[1][0][2][3] = 0;
    chrisD[1][0][3][0] = 0;
    chrisD[1][0][3][1] = 0;
    chrisD[1][0][3][2] = 0;
    chrisD[1][0][3][3] = 0;
    chrisD[1][1][0][0] = 0;
    chrisD[1][1][0][1] = 0;
    chrisD[1][1][0][2] = 0;
    chrisD[1][1][0][3] = 0;
    chrisD[1][1][1][0] = 0;
    chrisD[1][1][1][1] = 0;
    chrisD[1][1][1][2] = 0;
    chrisD[1][1][1][3] = 0;
    chrisD[1][1][2][0] = 0;
    chrisD[1][1][2][1] = 0;
    chrisD[1][1][2][2] = 0;
    chrisD[1][1][2][3] = 0;
    chrisD[1][1][3][0] = 0;
    chrisD[1][1][3][1] = 0;
    chrisD[1][1][3][2] = 0;
    chrisD[1][1][3][3] = 0;
    chrisD[1][2][0][0] = 0;
    chrisD[1][2][0][1] = 0;
    chrisD[1][2][0][2] = 0;
    chrisD[1][2][0][3] = 0;
    chrisD[1][2][1][0] = 0;
    chrisD[1][2][1][1] = 0;
    chrisD[1][2][1][2] = 0;
    chrisD[1][2][1][3] = 0;
    chrisD[1][2][2][0] = -t8;
    chrisD[1][2][2][1] = t8;
    chrisD[1][2][2][2] = 0;
    chrisD[1][2][2][3] = 0;
    chrisD[1][2][3][0] = 0;
    chrisD[1][2][3][1] = 0;
    chrisD[1][2][3][2] = 0;
    chrisD[1][2][3][3] = 0;
    chrisD[1][3][0][0] = 0;
    chrisD[1][3][0][1] = 0;
    chrisD[1][3][0][2] = 0;
    chrisD[1][3][0][3] = 0;
    chrisD[1][3][1][0] = 0;
    chrisD[1][3][1][1] = 0;
    chrisD[1][3][1][2] = 0;
    chrisD[1][3][1][3] = 0;
    chrisD[1][3][2][0] = 0;
    chrisD[1][3][2][1] = 0;
    chrisD[1][3][2][2] = 0;
    chrisD[1][3][2][3] = 0;
    chrisD[1][3][3][0] = -t15;
    chrisD[1][3][3][1] = t15;
    chrisD[1][3][3][2] = 0;
    chrisD[1][3][3][3] = 0;
    chrisD[2][0][0][0] = 0;
    chrisD[2][0][0][1] = 0;
    chrisD[2][0][0][2] = 0;
    chrisD[2][0][0][3] = 0;
    chrisD[2][0][1][0] = 0;
    chrisD[2][0][1][1] = 0;
    chrisD[2][0][1][2] = 0;
    chrisD[2][0][1][3] = 0;
    chrisD[2][0][2][0] = t8;
    chrisD[2][0][2][1] = -t8;
    chrisD[2][0][2][2] = 0;
    chrisD[2][0][2][3] = 0;
    chrisD[2][0][3][0] = 0;
    chrisD[2][0][3][1] = 0;
    chrisD[2][0][3][2] = 0;
    chrisD[2][0][3][3] = 0;
    chrisD[2][1][0][0] = 0;
    chrisD[2][1][0][1] = 0;
    chrisD[2][1][0][2] = 0;
    chrisD[2][1][0][3] = 0;
    chrisD[2][1][1][0] = 0;
    chrisD[2][1][1][1] = 0;
    chrisD[2][1][1][2] = 0;
    chrisD[2][1][1][3] = 0;
    chrisD[2][1][2][0] = -t8;
    chrisD[2][1][2][1] = t8;
    chrisD[2][1][2][2] = 0;
    chrisD[2][1][2][3] = 0;
    chrisD[2][1][3][0] = 0;
    chrisD[2][1][3][1] = 0;
    chrisD[2][1][3][2] = 0;
    chrisD[2][1][3][3] = 0;
    chrisD[2][2][0][0] = t16;
    chrisD[2][2][0][1] = -t16;
    chrisD[2][2][0][2] = 0;
    chrisD[2][2][0][3] = 0;
    chrisD[2][2][1][0] = t16;
    chrisD[2][2][1][1] = -t16;
    chrisD[2][2][1][2] = 0;
    chrisD[2][2][1][3] = 0;
    chrisD[2][2][2][0] = 0;
    chrisD[2][2][2][1] = 0;
    chrisD[2][2][2][2] = 0;
    chrisD[2][2][2][3] = 0;
    chrisD[2][2][3][0] = 0;
    chrisD[2][2][3][1] = 0;
    chrisD[2][2][3][2] = 0;
    chrisD[2][2][3][3] = 0;
    chrisD[2][3][0][0] = 0;
    chrisD[2][3][0][1] = 0;
    chrisD[2][3][0][2] = 0;
    chrisD[2][3][0][3] = 0;
    chrisD[2][3][1][0] = 0;
    chrisD[2][3][1][1] = 0;
    chrisD[2][3][1][2] = 0;
    chrisD[2][3][1][3] = 0;
    chrisD[2][3][2][0] = 0;
    chrisD[2][3][2][1] = 0;
    chrisD[2][3][2][2] = 0;
    chrisD[2][3][2][3] = 0;
    chrisD[2][3][3][0] = 0;
    chrisD[2][3][3][1] = 0;
    chrisD[2][3][3][2] = 0;
    chrisD[2][3][3][3] = 0;
    chrisD[3][0][0][0] = 0;
    chrisD[3][0][0][1] = 0;
    chrisD[3][0][0][2] = 0;
    chrisD[3][0][0][3] = 0;
    chrisD[3][0][1][0] = 0;
    chrisD[3][0][1][1] = 0;
    chrisD[3][0][1][2] = 0;
    chrisD[3][0][1][3] = 0;
    chrisD[3][0][2][0] = 0;
    chrisD[3][0][2][1] = 0;
    chrisD[3][0][2][2] = 0;
    chrisD[3][0][2][3] = 0;
    chrisD[3][0][3][0] = t15;
    chrisD[3][0][3][1] = -t15;
    chrisD[3][0][3][2] = 0;
    chrisD[3][0][3][3] = 0;
    chrisD[3][1][0][0] = 0;
    chrisD[3][1][0][1] = 0;
    chrisD[3][1][0][2] = 0;
    chrisD[3][1][0][3] = 0;
    chrisD[3][1][1][0] = 0;
    chrisD[3][1][1][1] = 0;
    chrisD[3][1][1][2] = 0;
    chrisD[3][1][1][3] = 0;
    chrisD[3][1][2][0] = 0;
    chrisD[3][1][2][1] = 0;
    chrisD[3][1][2][2] = 0;
    chrisD[3][1][2][3] = 0;
    chrisD[3][1][3][0] = -t15;
    chrisD[3][1][3][1] = t15;
    chrisD[3][1][3][2] = 0;
    chrisD[3][1][3][3] = 0;
    chrisD[3][2][0][0] = 0;
    chrisD[3][2][0][1] = 0;
    chrisD[3][2][0][2] = 0;
    chrisD[3][2][0][3] = 0;
    chrisD[3][2][1][0] = 0;
    chrisD[3][2][1][1] = 0;
    chrisD[3][2][1][2] = 0;
    chrisD[3][2][1][3] = 0;
    chrisD[3][2][2][0] = 0;
    chrisD[3][2][2][1] = 0;
    chrisD[3][2][2][2] = 0;
    chrisD[3][2][2][3] = 0;
    chrisD[3][2][3][0] = 0;
    chrisD[3][2][3][1] = 0;
    chrisD[3][2][3][2] = 0;
    chrisD[3][2][3][3] = 0;
    chrisD[3][3][0][0] = t17;
    chrisD[3][3][0][1] = -t17;
    chrisD[3][3][0][2] = 0;
    chrisD[3][3][0][3] = 0;
    chrisD[3][3][1][0] = t17;
    chrisD[3][3][1][1] = -t17;
    chrisD[3][3][1][2] = 0;
    chrisD[3][3][1][3] = 0;
    chrisD[3][3][2][0] = 0;
    chrisD[3][3][2][1] = 0;
    chrisD[3][3][2][2] = 0;
    chrisD[3][3][2][3] = 0;
    chrisD[3][3][3][0] = 0;
    chrisD[3][3][3][1] = 0;
    chrisD[3][3][3][2] = 0;
    chrisD[3][3][3][3] = 0;

    return true;
}

/*! Transform local 4-direction to coordinate 4-direction.
 *
 *  \param  pos  :  pointer to position array.
 *  \param  ldir :  pointer to local direction array.
 *  \param  dir  :  pointer to calculated coordinate direction array.
 *  \param  type :  type of tetrad.
 */
void MetricPlaneGravWave::localToCoord(const double* pos, const double* ldir, double* dir, enum_nat_tetrad_type)
{
    double c = mSpeedOfLight;
    double p = getValP(pos);
    double q = getValQ(pos);

    dir[0] = ldir[0] / c;
    dir[1] = ldir[1];

    dir[2] = ldir[2] / p;
    dir[3] = ldir[3] / q;
}

/*! Transform coordinate 4-direction to local 4-direction.
 *
 *  \param  pos  :  pointer to position array.
 *  \param  cdir :  pointer to coordinate direction.
 *  \param  ldir :  pointer to calculated local direction array.
 *  \param  type :  type of tetrad.
 */
void MetricPlaneGravWave::coordToLocal(const double* pos, const double* cdir, double* ldir, enum_nat_tetrad_type)
{
    double c = mSpeedOfLight;
    double p = getValP(pos);
    double q = getValQ(pos);

    ldir[0] = cdir[0] * c;
    ldir[1] = cdir[1];

    ldir[2] = cdir[2] * p;
    ldir[3] = cdir[3] * q;
}

/*! Test break condition.
 *
 *  \param pos    : pointer to position array.
 *  \return false : position is always valid.
 */
bool MetricPlaneGravWave::breakCondition(const double* pos)
{
    bool br = false;
    double c = mSpeedOfLight;
    double t = pos[0];
    double x = pos[1];

    if ((c * t - x) >= (1.0 - eps)) {
        br = true;
    }
    return br;
}

/*! Set parameter 'pName' to 'val'.
 *
 *	Set the longitudinal extension 'mLongExt' of the gravitational wave or
 *	the degree 'mDegree' of the Fourier polynomial and reset the state of
 *	'dataCalculated'.
 */
bool MetricPlaneGravWave::setParam(const char* pName, double val)
{
    Metric::setParam(pName, val);
    if (strcmp(pName, "long_ext") == 0) {
        if (val < eps) {
            mLongExt = eps;
            dataCalculated = false;
        }
        else {
            mLongExt = val;
            dataCalculated = false;
        }
    }
    else if (strcmp(pName, "degree") == 0) {
        if (val < 1.0) {
            mDegree = 1.0;
            dataCalculated = false;
        }
        else {
            mDegree = val;
            dataCalculated = false;
        }

        if (uData != nullptr) {
            delete[] uData;
        }
        if (dmData != nullptr) {
            delete[] dmData;
        }

        int n = static_cast<int>(mDegree);
        int N = 2 * n;
        uData = new double[N + 1];
        dmData = new double[N + 1];
    }

    return true;
}

bool MetricPlaneGravWave::transToTwoPlusOne(vec4 p, vec4& cp)
{
    vec4 tp;
    TransCoordinates::toCartesianCoord(mCoordType, p, tp);
    cp = vec4(tp[0], tp[1], tp[2], tp[0]);

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
double MetricPlaneGravWave::testConstraint(const double y[], const double kappa)
{
    double c = mSpeedOfLight;
    double cm = 1.0 / c;
    double p = getValP(y);
    double q = getValQ(y);

    // Scale the directions with the speed of light before doubling them !!
    double dt = y[4];
    double dx = y[5] * cm;
    double dy = y[6] * cm;
    double dz = y[7] * cm;

    double sum = -kappa;

    sum += -dt * dt + dx * dx + p * p * dy * dy + q * q * dz * dz;

    return sum;
}

/*! Generate report.
 */
bool MetricPlaneGravWave::report(const vec4, const vec4, std::string& text)
{
    std::stringstream ss;
    ss << "Report for PlaneGravWave metric\n\tcoordinate : (t,x,y,z)\n";
    ss << "---------------------------------------------------------------\n";
    ss << "  physical units ................................. no\n";
    ss.precision(DEF_FIXED_REPORT_PRECISION);
    ss.setf(std::ios::fixed);
    ss << "  Longitudinal Extension ............. a = " << mLongExt << std::endl;
    ss << "  Degree of the Fourier polynomial ... n = " << static_cast<int>(mDegree) << std::endl;

    text = ss.str();
    return true;
}

// ********************************* protected methods *****************************
/*!
 */
void MetricPlaneGravWave::setStandardValues()
{
    mInitPos[0] = 0.0;
    mInitPos[1] = 2.0;
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

// TODO: Methoden getValP usw. noch anschaulicher gestalten.
// ********************************* specific protected methods *****************************
/*! Calculate the value of function p(u) at 'pos'.
 */
double MetricPlaneGravWave::getValP(const double* pos)
{
    double c = mSpeedOfLight;
    double t = pos[0];
    double x = pos[1];
    double u = c * t - x;
    double p = 0.0;

    if (!dataCalculated) {
        calcFourierCoeff();
    }

    if (u < -mLongExt) {
        p = p0;
    }
    if (-mLongExt <= u && u <= 0.0) {
        int n = static_cast<int>(mDegree);
        double L = 0.0;
        double m = intConst;

        L = 1.0 - u + 1.0 / pow(mLongExt, 2) * pow(u, 3) + 1.0 / (2.0 * pow(mLongExt, 3)) * pow(u, 4);

        // Calculate Fourier polynomial for m(u).
        for (int k = 1; k < n; k++) {
            double k_d = static_cast<double>(k);

            m += -fCoeffB[k] * mLongExt / (k_d * M_PI) * cos(k_d * M_PI / mLongExt * u);
        }

        p = L * exp(m);
    }
    if (u > 0.0) {
        p = 1.0 - u;
    }

    return p;
}

/*! Calculate the value of function q(u) at 'pos'.
 */
double MetricPlaneGravWave::getValQ(const double* pos)
{
    double c = mSpeedOfLight;
    double t = pos[0];
    double x = pos[1];
    double u = c * t - x;
    double q = 0.0;

    if (!dataCalculated) {
        calcFourierCoeff();
    }

    if (u < -mLongExt) {
        q = q0;
    }
    if (-mLongExt <= u && u <= 0.0) {
        int n = static_cast<int>(mDegree);
        double L = 0.0;
        double m = intConst;

        L = 1.0 - u + 1.0 / pow(mLongExt, 2) * pow(u, 3) + 1.0 / (2.0 * pow(mLongExt, 3)) * pow(u, 4);

        // Calculate Fourier polynomial for m(u).
        for (int k = 1; k < n; k++) {
            double k_d = static_cast<double>(k);

            m += -fCoeffB[k] * mLongExt / (k_d * M_PI) * cos(k_d * M_PI / mLongExt * u);
        }

        q = L * exp(-m);
    }
    if (u > 0.0) {
        q = 1.0 - u;
    }

    return q;
}

/*! Calculate the value of function diff(p(u),u) at 'pos'.
 */
double MetricPlaneGravWave::getValDP(const double* pos)
{
    double c = mSpeedOfLight;
    double t = pos[0];
    double x = pos[1];
    double u = c * t - x;
    double dp;

    if (!dataCalculated) {
        calcFourierCoeff();
    }

    if (u < -mLongExt) {
        dp = 0.0;
    }
    if (-mLongExt <= u && u <= 0.0) {
        int n = static_cast<int>(mDegree);
        double L = 0.0;
        double dL = 0.0;
        double m = intConst;
        double dm = 0.0;

        L = 1.0 - u + 1.0 / pow(mLongExt, 2) * pow(u, 3) + 1.0 / (2.0 * pow(mLongExt, 3)) * pow(u, 4);
        dL = -1.0 + 3.0 / pow(mLongExt, 2) * pow(u, 2) + 2.0 / pow(mLongExt, 3) * pow(u, 3);

        // Calculate Fourier polynomial for m(u).
        for (int k = 1; k < n; k++) {
            double k_d = static_cast<double>(k);

            m += -fCoeffB[k] * mLongExt / (k_d * M_PI) * cos(k_d * M_PI / mLongExt * u);
        }
        // Calculate Fourier polynomial for dm/du.
        for (int k = 1; k < n; k++) {
            double k_d = static_cast<double>(k);

            dm += fCoeffB[k] * sin(k_d * M_PI / mLongExt * u);
        }

        dp = (dL + L * dm) * exp(m);
    }
    if (u > 0.0) {
        dp = -1.0;
    }

    return dp;
}

/*! Calculate the value of function diff(q(u),u) at 'pos'.
 */
double MetricPlaneGravWave::getValDQ(const double* pos)
{
    double c = mSpeedOfLight;
    double t = pos[0];
    double x = pos[1];
    double u = c * t - x;
    double dq;

    if (!dataCalculated) {
        calcFourierCoeff();
    }

    if (u < -mLongExt) {
        dq = 0.0;
    }
    if (-mLongExt <= u && u <= 0.0) {
        int n = static_cast<int>(mDegree);
        double L = 0.0;
        double dL = 0.0;
        double m = intConst;
        double dm = 0.0;

        L = 1.0 - u + 1.0 / pow(mLongExt, 2) * pow(u, 3) + 1.0 / (2.0 * pow(mLongExt, 3)) * pow(u, 4);
        dL = -1.0 + 3.0 / pow(mLongExt, 2) * pow(u, 2) + 2.0 / pow(mLongExt, 3) * pow(u, 3);

        // Calculate Fourier polynomial for m(u).
        for (int k = 1; k < n; k++) {
            double k_d = static_cast<double>(k);

            m += -fCoeffB[k] * mLongExt / (k_d * M_PI) * cos(k_d * M_PI / mLongExt * u);
        }
        // Calculate Fourier polynomial for dm/du.
        for (int k = 1; k < n; k++) {
            double k_d = static_cast<double>(k);

            dm += fCoeffB[k] * sin(k_d * M_PI / mLongExt * u);
        }

        dq = (dL - L * dm) * exp(-m);
    }
    if (u > 0.0) {
        dq = -1.0;
    }

    return dq;
}

/*! Calculate the value of function diff(diff(p(u),u),u) at 'pos'.
 */
double MetricPlaneGravWave::getValDDP(const double* pos)
{
    double c = mSpeedOfLight;
    double t = pos[0];
    double x = pos[1];
    double u = c * t - x;
    double ddp;

    if (!dataCalculated) {
        calcFourierCoeff();
    }

    if (u < -mLongExt) {
        ddp = 0.0;
    }
    if (-mLongExt <= u && u <= 0.0) {
        int n = static_cast<int>(mDegree);

        double u2 = u * u;
        double u3 = u * u * u;
        double u4 = u * u * u * u;
        double u5 = u * u * u * u * u;
        double mLongExt2 = pow(mLongExt, 2);
        double mLongExt3 = pow(mLongExt, 3);
        double mLongExt4 = pow(mLongExt, 4);

        double L = 0.0;
        double dL = 0.0;
        double ddL = 0.0;
        double m = intConst;
        double dm = 0.0;
        double ddm = 0.0;

        L = 1.0 - u + u3 / mLongExt2 + u4 / (2.0 * mLongExt3);
        dL = -1.0 + 3.0 * u2 / mLongExt2 + 2.0 * u3 / mLongExt3;
        ddL = 6.0 * (u / mLongExt2 + u2 / mLongExt3);
        ddm = -sqrt(3.0 / (u2 + u * mLongExt))
            * (-2.0 * mLongExt3 * u2 - 5.0 * mLongExt * u4 - 2.0 * u5 + 4.0 * mLongExt3 * u - 4.0 * mLongExt2 * u3
                + mLongExt4)
            / pow((2.0 * mLongExt3 * u - 2.0 * mLongExt * u3 - u4 - 2.0 * mLongExt3), 1.5);

        // Calculate Fourier polynomial for m(u).
        for (int k = 1; k < n; k++) {
            double k_d = static_cast<double>(k);

            m += -fCoeffB[k] * mLongExt / (k_d * M_PI) * cos(k_d * M_PI / mLongExt * u);
        }
        // Calculate Fourier polynomial for dm/du.
        for (int k = 1; k < n; k++) {
            double k_d = static_cast<double>(k);

            dm += fCoeffB[k] * sin(k_d * M_PI / mLongExt * u);
        }

        ddp = (ddL + 2.0 * dL * dm + L * ddm + L * dm * dm) * exp(m);
    }
    if (u > 0.0) {
        ddp = 0.0;
    }

    return ddp;
}

/*! Calculate the value of function diff(diff(q(u),u),u) at 'pos'.
 */
double MetricPlaneGravWave::getValDDQ(const double* pos)
{
    double c = mSpeedOfLight;
    double t = pos[0];
    double x = pos[1];
    double u = c * t - x;
    double ddq;

    if (!dataCalculated) {
        calcFourierCoeff();
    }

    if (u < -mLongExt) {
        ddq = 0.0;
    }
    if (-mLongExt <= u && u <= 0.0) {
        int n = static_cast<int>(mDegree);

        double u2 = u * u;
        double u3 = u * u * u;
        double u4 = u * u * u * u;
        double u5 = u * u * u * u * u;
        double mLongExt2 = pow(mLongExt, 2);
        double mLongExt3 = pow(mLongExt, 3);
        double mLongExt4 = pow(mLongExt, 4);

        double L = 0.0;
        double dL = 0.0;
        double ddL = 0.0;
        double m = intConst;
        double dm = 0.0;
        double ddm = 0.0;

        L = 1.0 - u + u3 / mLongExt2 + u4 / (2.0 * mLongExt3);
        dL = -1.0 + 3.0 * u2 / mLongExt2 + 2.0 * u3 / mLongExt3;
        ddL = 6.0 * (u / mLongExt2 + u2 / mLongExt3);
        ddm = -sqrt(3.0 / (u2 + u * mLongExt))
            * (-2.0 * mLongExt3 * u2 - 5.0 * mLongExt * u4 - 2.0 * u5 + 4.0 * mLongExt3 * u - 4.0 * mLongExt2 * u3
                + mLongExt4)
            / pow((2.0 * mLongExt3 * u - 2.0 * mLongExt * u3 - u4 - 2.0 * mLongExt3), 1.5);

        // Calculate Fourier polynomial for m(u).
        for (int k = 1; k < n; k++) {
            double k_d = static_cast<double>(k);

            m += -fCoeffB[k] * mLongExt / (k_d * M_PI) * cos(k_d * M_PI / mLongExt * u);
        }
        // Calculate Fourier polynomial for dm/du.
        for (int k = 1; k < n; k++) {
            double k_d = static_cast<double>(k);

            dm += fCoeffB[k] * sin(k_d * M_PI / mLongExt * u);
        }

        ddq = (ddL - 2.0 * dL * dm - L * ddm + L * dm * dm) * exp(-m);
    }
    if (u > 0.0) {
        ddq = 0.0;
    }

    return ddq;
}

/*! Calculate the Fourier coefficients to approximate the derivation
 * 	\f$ dm/du \f$ by a Fourier polynomial of degree n (here n = mDegree).
 *
 *  mLongExt : Longitudinal extension of the gravitational Wave.
 *  fCoeffB  : Dynamic array (std::vector) for the Fourier coefficients.
 *
 *  <B>Note</B>: Since \f$ dm/du \f$ can be continued as an odd
 *  						 \f$ 2\pi \f$-periodically function, all Fourier
 *  						 coefficients \f$ a_k \f$ vanish and only the coefficients
 *  						 \f$ b_k \f$ have to be calculated.
 */
void MetricPlaneGravWave::calcFourierCoeff()
{
    assert(uData != nullptr);
    assert(dmData != nullptr);
    int n = static_cast<int>(mDegree);
    int N = 2 * n;
    // double uData[N + 1];
    // double dmData[N + 1];

    intConst = 0.0;
    fCoeffB.clear();
    fCoeffB.resize(n + 1, 0.0);

    // Create a discrete set of values u_i on the interval [-1,+1]
    for (int i = 0; i <= n; i++) {
        double i_d = static_cast<double>(i);
        double n_d = static_cast<double>(n);

        uData[n - i] = -i_d / n_d * mLongExt;
        uData[n + i] = i_d / n_d * mLongExt;
    }

    // Create an odd continued data set of dm/du.
    for (int i = 0; i <= n; i++) {
        dmData[i] = -2.0 * sqrt(3.0)
            * sqrt((uData[i] * uData[i] + mLongExt * uData[i])
                / (2.0 * pow(mLongExt, 3) * uData[i] - 2.0 * mLongExt * pow(uData[i], 3) - pow(uData[i], 4)
                    - 2.0 * pow(mLongExt, 3)));

        dmData[n + i] = 2.0 * sqrt(3.0)
            * sqrt((uData[n + i] * uData[n + i] - mLongExt * uData[n + i])
                / (-2.0 * pow(mLongExt, 3) * uData[n + i] + 2.0 * mLongExt * pow(uData[n + i], 3) - pow(uData[n + i], 4)
                    - 2.0 * pow(mLongExt, 3)));
    }

    // Note: The variable u must be multiplied with the factor M_PI/mLongExt
    //			 to ensure that dm/du is continued 2*PI-periodically.
    for (int k = 0; k <= n; k++) {
        for (int i = 1; i <= N; i++) {
            double k_d = static_cast<double>(k);
            double N_d = static_cast<double>(N);

            fCoeffB[k] += 2.0 / N_d * dmData[i] * sin(k_d * M_PI / mLongExt * uData[i]);
        }
    }

    // Calculate constant of integration 'intConst' for m(u).
    for (int k = 1; k < n; k++) {
        double k_d = static_cast<double>(k);

        intConst += fCoeffB[k] * mLongExt / (k_d * M_PI);
    }

    // Calculate constant values p0 and q0. Therefore calculate
    // m0 := m(u = -mLongExt).
    double m0 = intConst;

    for (int k = 1; k < n; k++) {
        double k_d = static_cast<double>(k);

        m0 += -fCoeffB[k] * mLongExt / (k_d * M_PI) * cos(k_d * M_PI);
    }

    p0 = (1.0 + mLongExt / 2.0) * exp(m0);
    q0 = (1.0 + mLongExt / 2.0) * exp(-m0);

    dataCalculated = true;
}

} // end namespace m4d
