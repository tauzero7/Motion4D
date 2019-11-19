// -------------------------------------------------------------------------------
/*
   m4dMetricRotDihole.cpp

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

#include "m4dMetricRotDihole.h"

namespace m4d {
#define eps 1.0e-6

/*! Standard constructor for the Kottler metric.
 *
 * \param  mass1 : mass of black hole 1
 * \param  mass2 : mass of black hole 2
 * \param  omega : angular frequency
 */
MetricRotDihole::MetricRotDihole(double mass1, double mass2, double omega)
{
    mMetricName = "RotDihole";
    setCoordType(enum_coordinate_cartesian);

    mPhysicalUnits = enum_physical_constants_geom;
    mSpeedOfLight = 1.0;
    mGravConstant = 1.0;

    addParam("mass1", mass1);
    mMass1 = mass1;

    addParam("mass2", mass2);
    mMass2 = mass2;

    addParam("omega", omega);
    mOmega = omega;

    setStandardValues();
}

MetricRotDihole::~MetricRotDihole() {}

// *********************************** public methods ******************************
/*! Calculate the contravariant metric components at position 'pos'.
 *
 *  \param pos : pointer to position.
 */
bool MetricRotDihole::calculateMetric(const double* pos)
{
    double U = calc_U(pos);

    g_compts[0][0] = -1 / pow(U, 2);
    g_compts[0][1] = 0;
    g_compts[0][2] = 0;
    g_compts[0][3] = 0;
    g_compts[1][0] = 0;
    g_compts[1][1] = pow(U, 2);
    g_compts[1][2] = 0;
    g_compts[1][3] = 0;
    g_compts[2][0] = 0;
    g_compts[2][1] = 0;
    g_compts[2][2] = pow(U, 2);
    g_compts[2][3] = 0;
    g_compts[3][0] = 0;
    g_compts[3][1] = 0;
    g_compts[3][2] = 0;
    g_compts[3][3] = pow(U, 2);

    return true;
}

/*! Calculate the Christoffel symbols of the second kind at position 'pos'.
 *
 *  \param pos : pointer to position.
 */
bool MetricRotDihole::calculateChristoffels(const double* pos)
{
    double dUdt, dUdx, dUdy, dUdz;
    double U = calc_U(pos);
    calc_dU(pos, dUdt, dUdx, dUdy, dUdz);

    christoffel[0][0][0] = -dUdt / U;
    christoffel[0][0][1] = -dUdx / pow(U, 5);
    christoffel[0][0][2] = -dUdy / pow(U, 5);
    christoffel[0][0][3] = -dUdz / pow(U, 5);
    christoffel[0][1][0] = -dUdx / U;
    christoffel[0][1][1] = dUdt / U;
    christoffel[0][1][2] = 0;
    christoffel[0][1][3] = 0;
    christoffel[0][2][0] = -dUdy / U;
    christoffel[0][2][1] = 0;
    christoffel[0][2][2] = dUdt / U;
    christoffel[0][2][3] = 0;
    christoffel[0][3][0] = -dUdz / U;
    christoffel[0][3][1] = 0;
    christoffel[0][3][2] = 0;
    christoffel[0][3][3] = dUdt / U;
    christoffel[1][0][0] = -dUdx / U;
    christoffel[1][0][1] = dUdt / U;
    christoffel[1][0][2] = 0;
    christoffel[1][0][3] = 0;
    christoffel[1][1][0] = pow(U, 3) * dUdt;
    christoffel[1][1][1] = dUdx / U;
    christoffel[1][1][2] = -dUdy / U;
    christoffel[1][1][3] = -dUdz / U;
    christoffel[1][2][0] = 0;
    christoffel[1][2][1] = dUdy / U;
    christoffel[1][2][2] = dUdx / U;
    christoffel[1][2][3] = 0;
    christoffel[1][3][0] = 0;
    christoffel[1][3][1] = dUdz / U;
    christoffel[1][3][2] = 0;
    christoffel[1][3][3] = dUdx / U;
    christoffel[2][0][0] = -dUdy / U;
    christoffel[2][0][1] = 0;
    christoffel[2][0][2] = dUdt / U;
    christoffel[2][0][3] = 0;
    christoffel[2][1][0] = 0;
    christoffel[2][1][1] = dUdy / U;
    christoffel[2][1][2] = dUdx / U;
    christoffel[2][1][3] = 0;
    christoffel[2][2][0] = pow(U, 3) * dUdt;
    christoffel[2][2][1] = -dUdx / U;
    christoffel[2][2][2] = dUdy / U;
    christoffel[2][2][3] = -dUdz / U;
    christoffel[2][3][0] = 0;
    christoffel[2][3][1] = 0;
    christoffel[2][3][2] = dUdz / U;
    christoffel[2][3][3] = dUdy / U;
    christoffel[3][0][0] = -dUdz / U;
    christoffel[3][0][1] = 0;
    christoffel[3][0][2] = 0;
    christoffel[3][0][3] = dUdt / U;
    christoffel[3][1][0] = 0;
    christoffel[3][1][1] = dUdz / U;
    christoffel[3][1][2] = 0;
    christoffel[3][1][3] = dUdx / U;
    christoffel[3][2][0] = 0;
    christoffel[3][2][1] = 0;
    christoffel[3][2][2] = dUdz / U;
    christoffel[3][2][3] = dUdy / U;
    christoffel[3][3][0] = pow(U, 3) * dUdt;
    christoffel[3][3][1] = -dUdx / U;
    christoffel[3][3][2] = -dUdy / U;
    christoffel[3][3][3] = dUdz / U;

    return true;
}

/*! Calculate Jacobi matrix.
 *
 *  \param pos : pointer to position.
 */
bool MetricRotDihole::calculateChrisD(const double*)
{

    return true;
}

/*! Transform local 4-direction to coordinate 4-direction.
 *
 *  \param  pos  :  pointer to position array.
 *  \param  ldir :  pointer to local direction array.
 *  \param  dir  :  pointer to calculated coordinate direction array.
 *  \param  type :  type of tetrad.
 */
void MetricRotDihole::localToCoord(const double* pos, const double* ldir, double* dir, enum_nat_tetrad_type)
{
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
void MetricRotDihole::coordToLocal(const double* pos, const double* cdir, double* ldir, enum_nat_tetrad_type)
{
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
bool MetricRotDihole::breakCondition(const double*)
{
    bool br = false;
    return br;
}

/*! Set parameter 'pName' to 'val'.
 *
 *  Set 'mass' or 'lambda' parameter.
 */
bool MetricRotDihole::setParam(const char* pName, double val)
{
    Metric::setParam(pName, val);

    if (strcmp(pName, "mass1") == 0) {
        mMass1 = val;
    }
    else if (strcmp(pName, "mass2") == 0) {
        mMass2 = val;
    }
    else if (strcmp(pName, "omega") == 0) {
        mOmega = val;
    }
    return true;
}

/*! Generate report.
 */
bool MetricRotDihole::report(const vec4, const vec4, char*& text)
{
    std::stringstream ss;
    ss << "Report for RotDihole metric\n\tcoordinates : (t,x,y,z)\n";
    ss << "---------------------------------------------------------------\n";
    ss << "  physical units ......... no\n";
    ss.precision(DEF_FIXED_REPORT_PRECISION);
    ss.setf(std::ios::fixed);

    text = new char[ss.str().length() + 2];
    return CopyString(ss.str().c_str(), text);
}

// ********************************* protected methods *****************************
/*!
 */
void MetricRotDihole::setStandardValues()
{
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

void MetricRotDihole::calc_r(const double* pos, double& r1, double& r2)
{
    double t = pos[0];
    double x = pos[1];
    double y = pos[2];
    double z = pos[3];
    double omega = mOmega;

    r1 = sqrt(M4D_SQR(x) + M4D_SQR(y + sin(omega * t)) + M4D_SQR(z - cos(omega * t)));
    r2 = sqrt(M4D_SQR(x) + M4D_SQR(y - sin(omega * t)) + M4D_SQR(z + cos(omega * t)));
}

double MetricRotDihole::calc_U(const double* pos)
{
    double r1, r2;
    calc_r(pos, r1, r2);
    return 1.0 + mMass1 / r1 + mMass2 / r2;
}

void MetricRotDihole::calc_dU(const double* pos, double& dUdt, double& dUdx, double& dUdy, double& dUdz)
{
    double t = pos[0];
    double x = pos[1];
    double y = pos[2];
    double z = pos[3];

    double M1 = mMass1;
    double M2 = mMass2;
    double omega = mOmega;

    double r1, r2;
    calc_r(pos, r1, r2);

    dUdt = -M1 * omega * ((y + sin(omega * t)) * cos(omega * t) + (z - cos(omega * t)) * sin(omega * t)) / pow(r1, 3.0)
        + M2 * omega * ((y - sin(omega * t)) * cos(omega * t) + (z + cos(omega * t)) * sin(omega * t)) / pow(r2, 3.0);
    dUdx = -M1 * x / pow(r1, 3.0) - M2 * x / pow(r2, 3.0);
    dUdy = -M1 * (y + sin(omega * t)) / pow(r1, 3.0) - M2 * (y - sin(omega * t)) / pow(r2, 3.0);
    dUdz = -M1 * (z - cos(omega * t)) / pow(r1, 3.0) - M2 * (z + cos(omega * t)) / pow(r2, 3.0);
}

} // end namespace m4d
