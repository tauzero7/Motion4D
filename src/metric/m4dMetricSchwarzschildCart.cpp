// -------------------------------------------------------------------------------
/*
   m4dMetricSchwarzschildCart.cpp

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

using namespace std;

#include "m4dMetricSchwarzschildCart.h"

namespace m4d {

#define eps 1.0e-6


/*! Standard constructor for the Schwarzschild metric.
 *
 * \param  mass : mass of the black hole.
 */
MetricSchwarzschildCart::MetricSchwarzschildCart(double mass) {
    mMetricName  = "SchwarzschildCart";
    setCoordType(enum_coordinate_cartesian);

    mPhysicalUnits = enum_physical_constants_geom;
    mSpeedOfLight = 1.0;
    mGravConstant = 1.0;

    addParam("mass");
    setParam("mass", mass);

    rs = 2.0 * mass;

    setStandardValues();
}

MetricSchwarzschildCart::~MetricSchwarzschildCart() {
}


// *********************************** public methods ******************************

/*! Calculate the contravariant metric components at position 'pos'.
 *
 *  \param pos : pointer to position.
 */
bool MetricSchwarzschildCart::calculateMetric(const double* pos) {
    double x = pos[1];
    double y = pos[2];
    double z = pos[3];

    double c = 1.0;

    double t1 = c * c;
    double t3 = x * x;
    double t4 = y * y;
    double t5 = z * z;
    double t6 = t3 + t4 + t5;
    double t7 = sqrt(t6);
    double t8 = 1 / t7;
    double t11 = 1 / t6;
    double t12 = t11 * t3;
    double t15 = 1 / (1.0 - rs * t8);
    double t17 = t11 * t4;
    double t18 = t11 * t5;
    double t20 = rs * t11;
    double t22 = 1 / (t7 - rs);
    double t23 = t22 * x;
    double t25 = t20 * t23 * y;
    double t27 = t20 * t23 * z;
    double t32 = t20 * t22 * y * z;

    g_compts[0][0] = -t1 + t1 * rs * t8;
    g_compts[0][1] = 0.0;
    g_compts[0][2] = 0.0;
    g_compts[0][3] = 0.0;
    g_compts[1][0] = 0.0;
    g_compts[1][1] = t12 * t15 + t17 + t18;
    g_compts[1][2] = t25;
    g_compts[1][3] = t27;
    g_compts[2][0] = 0.0;
    g_compts[2][1] = t25;
    g_compts[2][2] = t12 + t17 * t15 + t18;
    g_compts[2][3] = t32;
    g_compts[3][0] = 0.0;
    g_compts[3][1] = t27;
    g_compts[3][2] = t32;
    g_compts[3][3] = t12 + t17 + t18 * t15;

    return true;
}

/*! Calculate the Christoffel symbols of the second kind at position 'pos'.
 *
 *  \param pos : pointer to position.
 */
bool MetricSchwarzschildCart::calculateChristoffels(const double* pos) {
    double x = pos[1];
    double y = pos[2];
    double z = pos[3];

    double c = 1.0;

    double t1 = c * c;
    double t2 = t1 * rs;
    double t3 = x * x;
    double t4 = y * y;
    double t5 = z * z;
    double t6 = t3 + t4 + t5;
    double t7 = sqrt(t6);
    double t8 = -t7 + rs;
    double t10 = t6 * t6;
    double t11 = 1 / t10;
    double t24 = 1 / t8;
    double t25 = 1 / t6 * t24;
    double t26 = rs * x;
    double t28 = t25 * t26 / 2.0;
    double t29 = rs * y;
    double t31 = t25 * t29 / 2.0;
    double t32 = rs * z;
    double t34 = t25 * t32 / 2.0;
    double t35 = t3 * t5;
    double t37 = 6.0 * t35 * t4;
    double t38 = t5 * t5;
    double t39 = t38 * t3;
    double t41 = t4 * t4;
    double t42 = t41 * t5;
    double t44 = t38 * t4;
    double t46 = t38 * t5;
    double t47 = 2.0 * t46;
    double t48 = t41 * t3;
    double t50 = t41 * t4;
    double t51 = 2.0 * t50;
    double t52 = t3 * t3;
    double t53 = t52 * t3;
    double t54 = t7 * t6;
    double t57 = 2.0 * t4 * t54 * rs;
    double t60 = 2.0 * t5 * t54 * rs;
    double t63 = t10 * t10;
    double t64 = 1 / t63;
    double t65 = (-t37 - 3.0 * t39 - 6.0 * t42 - 6.0 * t44 - t47 - 3.0 * t48 - t51 + t53 + t57 + t60) * t24 * t64;
    double t82 = t54 * rs;
    double t84 = -3.0 * t52 - 6.0 * t4 * t3 - 6.0 * t35 - 3.0 * t41 - 6.0 * t4 * t5 - 3.0 * t38 + 2.0 * t82;
    double t86 = t84 * t24 * t64;
    double t88 = y * t3 * rs * t86 / 2.0;
    double t92 = x * t4 * rs * t86 / 2.0;
    double t98 = t26 * y * z * t84 * t24 * t64 / 2.0;
    double t102 = z * t3 * rs * t86 / 2.0;
    double t106 = x * t5 * rs * t86 / 2.0;
    double t107 = 2.0 * t53;
    double t108 = t4 * t52;
    double t110 = t52 * t5;
    double t114 = 2.0 * t82 * t3;
    double t118 = (-t107 - 3.0 * t108 - 6.0 * t110 - t37 - 6.0 * t39 + t114 - 3.0 * t44 - t47 + t50 + t60) * t24 * t64;
    double t128 = t4 * rs * z * t86 / 2.0;
    double t132 = t5 * rs * y * t86 / 2.0;
    double t139 = (-t107 - 6.0 * t108 - 3.0 * t110 - 6.0 * t48 - t37 + t114 - t51 - 3.0 * t42 + t57 + t46) * t24 * t64;

    christoffel[0][0][0] = 0.0;
    christoffel[0][0][1] = -t2 * x * t8 * t11 / 2.0;
    christoffel[0][0][2] = -t2 * y * t8 * t11 / 2.0;
    christoffel[0][0][3] = -t2 * z * t8 * t11 / 2.0;
    christoffel[0][1][0] = -t28;
    christoffel[0][1][1] = 0.0;
    christoffel[0][1][2] = 0.0;
    christoffel[0][1][3] = 0.0;
    christoffel[0][2][0] = -t31;
    christoffel[0][2][1] = 0.0;
    christoffel[0][2][2] = 0.0;
    christoffel[0][2][3] = 0.0;
    christoffel[0][3][0] = -t34;
    christoffel[0][3][1] = 0.0;
    christoffel[0][3][2] = 0.0;
    christoffel[0][3][3] = 0.0;
    christoffel[1][0][0] = -t28;
    christoffel[1][0][1] = 0.0;
    christoffel[1][0][2] = 0.0;
    christoffel[1][0][3] = 0.0;
    christoffel[1][1][0] = 0.0;
    christoffel[1][1][1] = t26 * t65 / 2.0;
    christoffel[1][1][2] = t29 * t65 / 2.0;
    christoffel[1][1][3] = t32 * t65 / 2.0;
    christoffel[1][2][0] = 0.0;
    christoffel[1][2][1] = -t88;
    christoffel[1][2][2] = -t92;
    christoffel[1][2][3] = -t98;
    christoffel[1][3][0] = 0.0;
    christoffel[1][3][1] = -t102;
    christoffel[1][3][2] = -t98;
    christoffel[1][3][3] = -t106;
    christoffel[2][0][0] = -t31;
    christoffel[2][0][1] = 0.0;
    christoffel[2][0][2] = 0.0;
    christoffel[2][0][3] = 0.0;
    christoffel[2][1][0] = 0.0;
    christoffel[2][1][1] = -t88;
    christoffel[2][1][2] = -t92;
    christoffel[2][1][3] = -t98;
    christoffel[2][2][0] = 0.0;
    christoffel[2][2][1] = t26 * t118 / 2.0;
    christoffel[2][2][2] = t29 * t118 / 2.0;
    christoffel[2][2][3] = t32 * t118 / 2.0;
    christoffel[2][3][0] = 0.0;
    christoffel[2][3][1] = -t98;
    christoffel[2][3][2] = -t128;
    christoffel[2][3][3] = -t132;
    christoffel[3][0][0] = -t34;
    christoffel[3][0][1] = 0.0;
    christoffel[3][0][2] = 0.0;
    christoffel[3][0][3] = 0.0;
    christoffel[3][1][0] = 0.0;
    christoffel[3][1][1] = -t102;
    christoffel[3][1][2] = -t98;
    christoffel[3][1][3] = -t106;
    christoffel[3][2][0] = 0.0;
    christoffel[3][2][1] = -t98;
    christoffel[3][2][2] = -t128;
    christoffel[3][2][3] = -t132;
    christoffel[3][3][0] = 0.0;
    christoffel[3][3][1] = t26 * t139 / 2.0;
    christoffel[3][3][2] = t29 * t139 / 2.0;
    christoffel[3][3][3] = t32 * t139 / 2.0;

    //    coutChristoffel();
    //    abort();

    return true;
}

void MetricSchwarzschildCart::coutChristoffel(){
    cout << "\n";
    for(int i = 0; i<4; ++i){
        for(int j = 0; j<4; ++j){
            for(int k = 0; k<4; ++k){
                if(abs(christoffel[i][j][k]) > 1e-14){
                    cout << "christoffel[" << i << "][" << j << "][" << k << "]=" << christoffel[i][j][k] << "\n";
                }
            }
        }
    }
}


/*! Transform local 4-direction to coordinate 4-direction.
 *
 *  \param  pos  :  pointer to position array.
 *  \param  ldir :  pointer to local direction array.
 *  \param  dir  :  pointer to calculated coordinate direction array.
 *  \param  type :  type of tetrad.
 */
void MetricSchwarzschildCart::localToCoord(const double* pos, const double* ldir, double* dir,
        enum_nat_tetrad_type) {
    double x = pos[1];
    double y = pos[2];
    double z = pos[3];

    double c = 1.0;
    double r = sqrt(x * x + y * y + z * z);
    double f = 1.0 - rs / r;

    calcLTcoeffs(pos);

    dir[0] = 1.0 / (c * sqrt(f)) * ldir[0];
    dir[1] = (A * ldir[1] + B * ldir[2] + D * ldir[3]);
    dir[2] = (C * ldir[2] + E * ldir[3]);
    dir[3] = F * ldir[3];
}

/*! Transform coordinate 4-direction to local 4-direction.
 *
 *  \param  pos  :  pointer to position array.
 *  \param  cdir :  pointer to coordinate direction.
 *  \param  ldir :  pointer to calculated local direction array.
 *  \param  type :  type of tetrad.
 */
void MetricSchwarzschildCart::coordToLocal(const double* pos, const double* cdir, double* ldir,
        enum_nat_tetrad_type) {
    double x = pos[1];
    double y = pos[2];
    double z = pos[3];

    double c = 1.0;
    double r = sqrt(x * x + y * y + z * z);
    double f = 1.0 - rs / r;

    calcLTcoeffs(pos);

    ldir[0] = c * sqrt(f) * cdir[0];
    ldir[1] = 1.0 / A * cdir[1] - B / (A * C) * cdir[2] + (B * E - C * D) / (A * C * F) * cdir[3];
    ldir[2] = 1.0 / C * cdir[2] - E / (C * F) * cdir[3];
    ldir[3] = 1.0 / F * cdir[3];
}


/*! Test break condition.
 *
 *  \param pos    : pointer to position array.
 *  \return true  : radial position r < 0.0 or  r^2<=(1.0+eps)*rs^2.
 *  \return false : position is valid.
 */
bool MetricSchwarzschildCart::breakCondition(const double* pos) {
    bool br = false;
    double r2 = pos[1] * pos[1] + pos[2] * pos[2] + pos[3] * pos[3];
//if (r2 < 0.01) {
    //  std::cerr << pos[1] << " " << pos[2] << " " << pos[3] << " " << r2 << " " << rs*rs << std::endl;
//}
    if (r2 <= (1.0 + eps)*rs * rs) {
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
double MetricSchwarzschildCart::testConstraint(const double y[], const double kappa) {
    double xx = y[1];
    double yy = y[2];
    double zz = y[3];
    double r = sqrt(xx * xx + yy * yy + zz * zz);
    double c = 1.0;

    double edrq = 1.0 / (r * r);
    double f    = 1.0 - rs / r;

    double sum = -kappa;

    sum += - c * c * f * y[4] * y[4]
           + (xx * xx / f + yy * yy + zz * zz) * edrq * y[5] * y[5]
           + (xx * xx + yy * yy / f + zz * zz) * edrq * y[6] * y[6]
           + (xx * xx + yy * yy + zz * zz / f) * edrq * y[7] * y[7]
           + 2.0 * rs / (r * r * (r - rs)) * (xx * yy * y[5] * y[6] + xx * zz * y[5] * y[7] + yy * zz * y[6] * y[7]);

    return sum;
}

/*! Set parameter 'pName' to 'val'.
 *
 *  Set 'mass' parameter and adjust Schwarzschild radius  rs=2GM/c^2.
 */
bool MetricSchwarzschildCart::setParam(const char* pName, double val) {
    if (Metric::setParam(pName, val)) {
        rs = 2.0 * mGravConstant * val / (mSpeedOfLight * mSpeedOfLight);
    }
    return true;
}

/*! Generate report.
 */
bool MetricSchwarzschildCart::report(const vec4 , const vec4 , std::string &text) {
    std::stringstream ss;
    ss << "Report for Schwarzschild metric\n\tcoordinate : (t,x,y,z)\n";
    ss << "---------------------------------------------------------------\n";
    ss << "  physical units ................................. no\n";
    ss.precision(DEF_FIXED_REPORT_PRECISION);
    ss.setf(std::ios::fixed);
    ss << "  Schwarzschild radius ........... r_s = 2GM/c^2 = " << rs << std::endl;
    ss << "  Photon orbit ................... r_ph = 3/2*rs = " << 1.5 * rs << std::endl;
    ss << "  innermost stable circular orbit  r_isco = 3r_s = " << 3.0 * rs << std::endl;

    text = ss.str();
    return true;
}

// ********************************* protected methods *****************************
/*!
 */
void MetricSchwarzschildCart::setStandardValues() {
    mInitPos[0] = 0.0;
    mInitPos[1] = 3.0 * rs;
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

/*! Calculate local tetrad coefficients at position 'pos'.
 *
 *  \param pos : pointer to position.
 */
void MetricSchwarzschildCart::calcLTcoeffs(const double* pos) {
    calculateMetric(pos);

    double gxx = g_compts[1][1];
    double gyy = g_compts[2][2];
    double gzz = g_compts[3][3];

    double gxy = g_compts[1][2];
    double gxz = g_compts[1][3];
    double gyz = g_compts[2][3];

    double W = gxx * gyy * gzz - gxz * gxz * gyy + 2.0 * gxz * gxy * gyz - gxy * gxy * gzz - gxx * gyz * gyz;
    double N = gxx * gyy - gxy * gxy;

    double sW = sqrt(W);
    double sN = sqrt(N);

    A = 1.0 / sqrt(gxx);
    C = 1.0 / sqrt(-gxy * gxy / gxx + gyy);
    B = -gxy / gxx * C;
    F = sN / sW;
    E = (gxz * gxy - gxx * gyz) / (sN * sW);
    D = (gxy * gyz - gxz * gyy) / (sN * sW);
}

} // end namespace m4d
