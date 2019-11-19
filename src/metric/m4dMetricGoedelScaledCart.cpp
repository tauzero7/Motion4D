// -------------------------------------------------------------------------------
/*
    m4dMetricGoedelScaledCart.cpp

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

#include "m4dMetricGoedelScaledCart.h"

namespace m4d {

// ---------------------------------------------------
//    constructur/destructor
// ---------------------------------------------------
MetricGoedelScaledCart::MetricGoedelScaledCart(double rG, double zeta)
{
    mMetricName = "Goedel scaled cartesian";
    setCoordType(enum_coordinate_cartesian);

    mPhysicalUnits = enum_physical_constants_geom;
    mSpeedOfLight = 1.0;
    mGravConstant = 1.0;

    addParam("rG");
    setParam("rG", rG);
    addParam("zeta");
    setParam("zeta", zeta);

    mDrawTypes.push_back(enum_draw_twoplusone);

    setStandardValues();
}

MetricGoedelScaledCart::~MetricGoedelScaledCart() {}

// *********************************** public methods ******************************
// ---------------------------------------------------
//    public::calculateMetric
// ---------------------------------------------------
bool MetricGoedelScaledCart::calculateMetric(const double* pos)
{
    double X = pos[1];
    double Y = pos[2];
    double r_G = mRG;
    double c = mSpeedOfLight;

    double t1 = r_G * r_G;
    double t2 = X * X;
    double t3 = Y * Y;
    double t6 = t1 / (1.0 + t2 + t3);
    double t7 = c * c;
    double t14 = sqrt(2.0);
    double t15 = t6 * t14;
    double t19 = t14 * c;
    double t22 = t3 * Y;
    double t25 = t15 * c * Y * t2 + t6 * t19 * Y + t6 * t19 * t22;
    double t26 = t2 * X;
    double t34 = -t6 * t19 * t26 - t15 * c * X * t3 - t6 * t19 * X;
    double t36 = t6 * t3 * t2;
    double t37 = t3 * t3;
    double t44 = t6 * t26 * Y + t6 * t22 * X;
    double t45 = t2 * t2;

    g_compts[0][0] = -t6 * t7 * t2 - t6 * t7 * t3 - t6 * t7;
    g_compts[0][1] = t25;
    g_compts[0][2] = t34;
    g_compts[0][3] = 0.0;
    g_compts[1][0] = t25;
    g_compts[1][1] = -t36 + t6 - t6 * t37;
    g_compts[1][2] = t44;
    g_compts[1][3] = 0.0;
    g_compts[2][0] = t34;
    g_compts[2][1] = t44;
    g_compts[2][2] = -t6 * t45 - t36 + t6;
    g_compts[2][3] = 0.0;
    g_compts[3][0] = 0.0;
    g_compts[3][1] = 0.0;
    g_compts[3][2] = 0.0;
    g_compts[3][3] = t6 * t2 + t6 + t6 * t3;

    return true;
}

// ---------------------------------------------------
//    public::calculateChristoffels
// ---------------------------------------------------
bool MetricGoedelScaledCart::calculateChristoffels(const double* pos)
{
    double X = pos[1];
    double Y = pos[2];
    double c = mSpeedOfLight;

    double t1 = X * X;
    double t2 = Y * Y;
    double t4 = 1.0 / (1.0 + t1 + t2);
    double t5 = X * t4;
    double t6 = 2.0 * t5;
    double t10 = sqrt(2.0);
    double t13 = (2.0 + t1 + t2) * X * Y * t4 * t10 * c;
    double t14 = t2 * t1;
    double t15 = t2 * t2;
    double t19 = t10 * c;
    double t20 = (t14 + 1.0 + t15 + 2.0 * t2) * t4 * t19;
    double t22 = 2.0 * Y * t4;
    double t23 = t1 * t1;
    double t27 = (t23 + t14 + 1.0 + 2.0 * t1) * t4 * t19;
    double t29 = 1.0 / c;
    double t32 = 2.0 * Y * t10 * t5 * t29;
    double t33 = 2.0 * t14;
    double t34 = 3.0 * t2;
    double t35 = 2.0 * t15;
    double t39 = t33 + t35 + 1.0 + t34;
    double t45 = (t1 - t2) * t10 * t4 * t29;
    double t46 = 2.0 * t23;
    double t47 = 3.0 * t1;
    double t48 = t46 + t47 + t33 + 1.0;
    double t50 = Y * t48 * t4;
    double t52 = X * t39 * t4;

    christoffel[0][0][0] = 0.0;
    christoffel[0][0][1] = 0.0;
    christoffel[0][0][2] = 0.0;
    christoffel[0][0][3] = 0.0;
    christoffel[0][1][0] = t6;
    christoffel[0][1][1] = -t13;
    christoffel[0][1][2] = -t20;
    christoffel[0][1][3] = 0.0;
    christoffel[0][2][0] = t22;
    christoffel[0][2][1] = t27;
    christoffel[0][2][2] = t13;
    christoffel[0][2][3] = 0.0;
    christoffel[0][3][0] = 0.0;
    christoffel[0][3][1] = 0.0;
    christoffel[0][3][2] = 0.0;
    christoffel[0][3][3] = 0.0;
    christoffel[1][0][0] = t6;
    christoffel[1][0][1] = -t13;
    christoffel[1][0][2] = -t20;
    christoffel[1][0][3] = 0.0;
    christoffel[1][1][0] = -t32;
    christoffel[1][1][1] = X * (t33 + t34 + t35 - 1.0) * t4;
    christoffel[1][1][2] = Y * t39 * t4;
    christoffel[1][1][3] = 0.0;
    christoffel[1][2][0] = t45;
    christoffel[1][2][1] = -t50;
    christoffel[1][2][2] = -t52;
    christoffel[1][2][3] = 0.0;
    christoffel[1][3][0] = 0.0;
    christoffel[1][3][1] = 0.0;
    christoffel[1][3][2] = 0.0;
    christoffel[1][3][3] = 0.0;
    christoffel[2][0][0] = t22;
    christoffel[2][0][1] = t27;
    christoffel[2][0][2] = t13;
    christoffel[2][0][3] = 0.0;
    christoffel[2][1][0] = t45;
    christoffel[2][1][1] = -t50;
    christoffel[2][1][2] = -t52;
    christoffel[2][1][3] = 0.0;
    christoffel[2][2][0] = t32;
    christoffel[2][2][1] = X * t48 * t4;
    christoffel[2][2][2] = Y * (t46 + t47 + t33 - 1.0) * t4;
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
    christoffel[3][1][3] = 0.0;
    christoffel[3][2][0] = 0.0;
    christoffel[3][2][1] = 0.0;
    christoffel[3][2][2] = 0.0;
    christoffel[3][2][3] = 0.0;
    christoffel[3][3][0] = 0.0;
    christoffel[3][3][1] = 0.0;
    christoffel[3][3][2] = 0.0;
    christoffel[3][3][3] = 0.0;

    return true;
}

// ---------------------------------------------------
//    public::localToCoord
// ---------------------------------------------------
void MetricGoedelScaledCart::localToCoord(const double* pos, const double* ldir, double* dir, enum_nat_tetrad_type)
{
    // TODO

    double X = pos[1];
    double Y = pos[2];
    double R = sqrt(X * X + Y * Y);
    double R2 = R * R;
    double c = mSpeedOfLight;
    double c2 = c * c;

    double A = R2 * (-sqrt(2.0) * c + (1.0 - R2) * mZeta);
    double B = c2 + sqrt(2.0) * R2 * c * mZeta;
    double Gamma = 1.0 / sqrt(c2 + 2.0 * sqrt(2.0) * R2 * c * mZeta - R2 * (1.0 - R2) * mZeta * mZeta);
    double Delta = 1.0 / (R * c * sqrt(1.0 + R2));

    if (R > m4dGoedelScaledCartEps) {
        // local to coord to cylindrical coordinates + transformation to cartesian coordinate (direction)

        dir[0] = Gamma / mRG * ldir[0] + Delta * Gamma * A / mRG * ldir[2];
        dir[1] = 1.0 / mRG * ldir[1] * sqrt(1.0 + R * R) * X / R
            + (Gamma * mZeta / mRG * ldir[0] + Delta * Gamma * B / mRG * ldir[2]) * (-Y);
        dir[2] = 1.0 / mRG * ldir[1] * sqrt(1.0 + R * R) * Y / R
            + (Gamma * mZeta / mRG * ldir[0] + Delta * Gamma * B / mRG * ldir[2]) * X;
        dir[3] = 1.0 / mRG * ldir[3];
    }
    else {
        // metric converges to minkowski metric
        dir[0] = ldir[0];
        dir[1] = ldir[1];
        dir[2] = ldir[2];
        dir[3] = ldir[3];
    }
}
// ---------------------------------------------------
//    public::coordToLocal
// ---------------------------------------------------
void MetricGoedelScaledCart::coordToLocal(const double*, const double*, double*, enum_nat_tetrad_type)
{
    printf("MetricGoedelScaledCart::coordToLocal missing\n");
}

// ---------------------------------------------------
//    public::breakCondition
// ---------------------------------------------------
bool MetricGoedelScaledCart::breakCondition(const double*)
{
    return false;
}

// ---------------------------------------------------
//    public::setParam
// ---------------------------------------------------
bool MetricGoedelScaledCart::setParam(const char* pName, double val)
{
    Metric::setParam(pName, val);
    if (strcmp(pName, "rG") == 0) {
        mRG = val;
    }
    else if (strcmp(pName, "zeta") == 0) {
        mZeta = val;
    }

    // printf("a=%f, zeta=%f\n",mA,mZeta);

    return true;
}

/*! Transform point p to 2+1 coordinates.
 *
 *  \param  p  : point in proper metric coordinates.
 *  \param  cp : reference to transformed point.
 *  \return true : success.
 */
bool MetricGoedelScaledCart::transToTwoPlusOne(vec4 p, vec4& cp)
{
    cp = vec4(p[0], p[1], p[2], p[0]);
    return true;
}

/*! Generate report.
 */
bool MetricGoedelScaledCart::report(const vec4 pos, const vec4 cdir, char*& text)
{
    std::stringstream ss;
    ss << "Report for Goedel metric\n\tcoordinate : (T,X,Y,Z)\n";
    ss << "\tScaled coordinates, geodesical shape independent of rG\n";
    ss << "\talso valid at X=Y=0 (no coordinate singularity\n";
    ss << "---------------------------------------------------------------\n";
    ss << "  physical units ........... yes\n";
    ss.precision(DEF_FIXED_REPORT_PRECISION);
    ss.setf(std::ios::fixed);
    ss.precision(DEF_FIXED_REPORT_PRECISION);
    ss.setf(std::ios::fixed);
    vec4 locStartDir;
    coordToLocal(pos.data(), cdir.data(), locStartDir.data());
    double k0 = -locStartDir[0] / mRG;
    double k2
        = (pos[1] * sqrt(1.0 + pos[1] * pos[1]) * locStartDir[2] - sqrt(2.0) * pos[1] * pos[1] * locStartDir[0]) / mRG;
    double k3 = locStartDir[3] / mRG;
    ss << "  Goedel radius ............ rG  = " << mRG << std::endl;
    ss << "  constant of motion ....... k_0 = " << k0 << std::endl;
    ss << "  constant of motion ....... k_2 = " << k2 << std::endl;
    ss << "  constant of motion ....... k_3 = " << k3;

    text = new char[ss.str().length() + 2];
    return CopyString(ss.str().c_str(), text);
}

// ********************************* protected methods *****************************
// ---------------------------------------------------
//    protected::setViewerVal
// ---------------------------------------------------
void MetricGoedelScaledCart::setStandardValues()
{
    mInitPos[0] = 0.0;
    mInitPos[1] = 0.0;
    mInitPos[2] = 0.0;
    mInitPos[3] = 0.0;
    mInitDir[0] = 1.0;
    mInitDir[1] = 0.0;
    mInitDir[2] = 0.0;

    mCoordNames[0] = std::string("T");
    mCoordNames[1] = std::string("X");
    mCoordNames[2] = std::string("Y");
    mCoordNames[3] = std::string("Z");
}

} // end namespace m4d
