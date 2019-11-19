// -------------------------------------------------------------------------------
/*
    m4dMetricHalilsoyWave.cpp

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

#include "m4dMetricHalilsoyWave.h"

namespace m4d {

#define eps 1.0e-6

/*! Standard constructor for the GravWave metric.
 */
MetricHalilsoyWave::MetricHalilsoyWave(double alpha, double C)
{
    mMetricName = "HalilsoyWave";
    setCoordType(enum_coordinate_cylinder);

    mPhysicalUnits = enum_physical_constants_geom;
    mSpeedOfLight = 1.0;
    mGravConstant = 1.0;

    mAlpha = alpha;
    mC = C;

    addParam("alpha", alpha);
    addParam("c", C);

    setStandardValues();
}

MetricHalilsoyWave::~MetricHalilsoyWave() {}

// *********************************** public methods ******************************
/*! Calculate the contravariant metric components at position 'pos'.
 *
 *  \param pos : pointer to position.
 */
bool MetricHalilsoyWave::calculateMetric(const double* pos)
{
    double rho = pos[1];
    calcMetricFunc(pos);

    double t1 = fV; // V(t,rho);
    double t2 = fK; // K(t,rho);
    double t3 = exp(t2);
    double t4 = t3 * t3;
    double t5 = t1 * t4;
    double t6 = rho * rho;
    double t8 = 1 / t1;
    double t9 = fA; // A(t,rho);
    double t10 = t9 * t9;
    double t13 = t8 * t9;

    g_compts[0][0] = -t5;
    g_compts[0][1] = 0.0;
    g_compts[0][2] = 0.0;
    g_compts[0][3] = 0.0;
    g_compts[1][0] = 0.0;
    g_compts[1][1] = t5;
    g_compts[1][2] = 0.0;
    g_compts[1][3] = 0.0;
    g_compts[2][0] = 0.0;
    g_compts[2][1] = 0.0;
    g_compts[2][2] = t1 * t6 + t8 * t10;
    g_compts[2][3] = t13;
    g_compts[3][0] = 0.0;
    g_compts[3][1] = 0.0;
    g_compts[3][2] = t13;
    g_compts[3][3] = t8;

    return true;
}

/*! Calculate the Christoffel symbols of the second kind at position 'pos'.
 *
 *  \param pos : pointer to position.
 */
bool MetricHalilsoyWave::calculateChristoffels(const double* pos)
{
    double rho = pos[1];
    calcMetricFunc(pos);
    calcDiffMetricFunc(pos);

    double t1 = fV_t; // diff(V(t,rho),t);
    double t2 = fV; // V(t,rho);
    double t3 = fK_t; // diff(K(t,rho),t);
    double t7 = 1 / t2;
    double t9 = (t1 + 2.0 * t2 * t3) * t7 / 2.0;
    double t10 = fV_rho; // diff(V(t,rho),rho);
    double t11 = fK_rho; // diff(K(t,rho),rho);
    double t16 = (t10 + 2.0 * t2 * t11) * t7 / 2.0;
    double t17 = t2 * t2;
    double t18 = 1 / t17;
    double t19 = rho * rho;
    double t20 = t1 * t19;
    double t22 = fA; // A(t,rho);
    double t23 = fA_t; // diff(A(t,rho),t);
    double t24 = t22 * t23;
    double t27 = 1 / t19;
    double t29 = t18 * (t20 * t2 + t24) * t27 / 2.0;
    double t30 = t22 * t1;
    double t31 = t19 * t2;
    double t34 = t22 * t22;
    double t36 = t17 * t19;
    double t41 = t18 * (2.0 * t30 * t31 + t34 * t23 - t36 * t23) * t27 / 2.0;
    double t44 = t23 * t18 * t27 / 2.0;
    double t45 = t10 * t19;
    double t46 = t45 * t2;
    double t49 = fA_rho; // diff(A(t,rho),rho);
    double t50 = t22 * t49;
    double t54 = t18 * (t46 + 2.0 * t17 * rho + t50) * t27 / 2.0;
    double t55 = t22 * t10;
    double t66 = t18 * (2.0 * t55 * t31 + 2.0 * t22 * t17 * rho + t34 * t49 - t36 * t49) * t27 / 2.0;
    double t69 = t49 * t18 * t27 / 2.0;
    double t73 = t18 * (t50 + t46) * t27 / 2.0;
    double t74 = fK; // K(t,rho);
    double t75 = exp(t74);
    double t76 = t75 * t75;
    double t78 = t17 * t2;
    double t80 = 1 / t76 / t78;
    double t100 = t80 * (t30 - t23 * t2) / 2.0;
    double t104 = t80 * (t55 - t49 * t2) / 2.0;

    christoffel[0][0][0] = t9;
    christoffel[0][0][1] = t16;
    christoffel[0][0][2] = 0.0;
    christoffel[0][0][3] = 0.0;
    christoffel[0][1][0] = t16;
    christoffel[0][1][1] = t9;
    christoffel[0][1][2] = 0.0;
    christoffel[0][1][3] = 0.0;
    christoffel[0][2][0] = 0.0;
    christoffel[0][2][1] = 0.0;
    christoffel[0][2][2] = t29;
    christoffel[0][2][3] = -t41;
    christoffel[0][3][0] = 0.0;
    christoffel[0][3][1] = 0.0;
    christoffel[0][3][2] = t44;
    christoffel[0][3][3] = -t29;
    christoffel[1][0][0] = t16;
    christoffel[1][0][1] = t9;
    christoffel[1][0][2] = 0.0;
    christoffel[1][0][3] = 0.0;
    christoffel[1][1][0] = t9;
    christoffel[1][1][1] = t16;
    christoffel[1][1][2] = 0.0;
    christoffel[1][1][3] = 0.0;
    christoffel[1][2][0] = 0.0;
    christoffel[1][2][1] = 0.0;
    christoffel[1][2][2] = t54;
    christoffel[1][2][3] = -t66;
    christoffel[1][3][0] = 0.0;
    christoffel[1][3][1] = 0.0;
    christoffel[1][3][2] = t69;
    christoffel[1][3][3] = -t73;
    christoffel[2][0][0] = 0.0;
    christoffel[2][0][1] = 0.0;
    christoffel[2][0][2] = t29;
    christoffel[2][0][3] = -t41;
    christoffel[2][1][0] = 0.0;
    christoffel[2][1][1] = 0.0;
    christoffel[2][1][2] = t54;
    christoffel[2][1][3] = -t66;
    christoffel[2][2][0] = t80 * (t20 * t17 - t34 * t1 + 2.0 * t24 * t2) / 2.0;
    christoffel[2][2][1] = -t80 * (t45 * t17 + 2.0 * t78 * rho - t34 * t10 + 2.0 * t50 * t2) / 2.0;
    christoffel[2][2][2] = 0.0;
    christoffel[2][2][3] = 0.0;
    christoffel[2][3][0] = -t100;
    christoffel[2][3][1] = t104;
    christoffel[2][3][2] = 0.0;
    christoffel[2][3][3] = 0.0;
    christoffel[3][0][0] = 0.0;
    christoffel[3][0][1] = 0.0;
    christoffel[3][0][2] = t44;
    christoffel[3][0][3] = -t29;
    christoffel[3][1][0] = 0.0;
    christoffel[3][1][1] = 0.0;
    christoffel[3][1][2] = t69;
    christoffel[3][1][3] = -t73;
    christoffel[3][2][0] = -t100;
    christoffel[3][2][1] = t104;
    christoffel[3][2][2] = 0.0;
    christoffel[3][2][3] = 0.0;
    christoffel[3][3][0] = -t80 * t1 / 2.0;
    christoffel[3][3][1] = t80 * t10 / 2.0;
    christoffel[3][3][2] = 0.0;
    christoffel[3][3][3] = 0.0;

    return true;
}

/*! Calculate Jacobi matrix.
 *
 *  \param pos : pointer to position.
 */
bool MetricHalilsoyWave::calculateChrisD(const double*)
{
    return false;
}

/*! Transform local 4-direction to coordinate 4-direction.
 *
 *  \param  pos  :  pointer to position array.
 *  \param  ldir :  pointer to local direction array.
 *  \param  dir  :  pointer to calculated coordinate direction array.
 *  \param  type :  type of tetrad.
 */
void MetricHalilsoyWave::localToCoord(const double* pos, const double* ldir, double* dir, enum_nat_tetrad_type)
{
    double rho = pos[1];
    calcMetricFunc(pos);

    double sV = sqrt(fV);
    dir[0] = ldir[0] * exp(-fK) / sV;
    dir[1] = ldir[1] * exp(-fK) / sV;
    dir[2] = ldir[2] / (rho * sV);
    dir[3] = -ldir[2] * fA / (rho * sV) + ldir[3] * sV;
}

/*! Transform coordinate 4-direction to local 4-direction.
 *
 *  \param  pos  :  pointer to position array.
 *  \param  cdir :  pointer to coordinate direction.
 *  \param  ldir :  pointer to calculated local direction array.
 *  \param  type :  type of tetrad.
 */
void MetricHalilsoyWave::coordToLocal(const double* pos, const double* cdir, double* ldir, enum_nat_tetrad_type)
{
    double rho = pos[1];
    calcMetricFunc(pos);

    double sV = sqrt(fV);
    ldir[0] = cdir[0] * exp(fK) * sV;
    ldir[1] = cdir[1] * exp(fK) * sV;
    ldir[2] = cdir[2] * rho * sV;
    ldir[3] = (cdir[3] + fA * cdir[2]) / sV;
}

/*! Test break condition.
 *
 *  \param pos    : pointer to position array.
 *  \return false : position is always valid.
 */
bool MetricHalilsoyWave::breakCondition(const double*)
{
    bool br = false;
    return br;
}

/*! Set parameter 'pName' to 'val'.
 *
 */
bool MetricHalilsoyWave::setParam(const char* pName, double val)
{
    Metric::setParam(pName, val);
    if (strcmp(pName, "alpha") == 0) {
        mAlpha = val;
    }
    else if (strcmp(pName, "c") == 0) {
        mC = val;
    }
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
double MetricHalilsoyWave::testConstraint(const double y[], const double kappa)
{
    double rho = y[1];
    double dt = y[4];
    double drho = y[5];
    double dphi = y[6];
    double dz = y[7];

    double sum = -kappa;

    calcMetricFunc(y);
    sum += fV * (exp(2.0 * fK) * (drho * drho - dt * dt) + rho * rho * dphi * dphi)
        + (dz + fA * dphi) * (dz + fA * dphi) / fV;
    return sum;
}

/*! Generate report.
 */
bool MetricHalilsoyWave::report(const vec4, const vec4, char*& text)
{
    std::stringstream ss;
    ss << "Report for HalolsoyWave metric\n\tcoordinate : (t,rho,phi,z)\n";
    ss << "---------------------------------------------------------------\n";
    ss << "  physical units ................................. no\n";
    ss.precision(DEF_FIXED_REPORT_PRECISION);
    ss.setf(std::ios::fixed);
    ss << "  alpha ...................................... = " << mAlpha << std::endl;
    ss << "  C .......................................... = " << mC << std::endl;
    text = new char[ss.str().length() + 2];
    return CopyString(ss.str().c_str(), text);
}

// ********************************* protected methods *****************************
/*!
 */
void MetricHalilsoyWave::setStandardValues()
{
    mInitPos[0] = 0.0;
    mInitPos[1] = 5.0;
    mInitPos[2] = 0.0;
    mInitPos[3] = 0.0;
    mInitDir[0] = 1.0;
    mInitDir[1] = 0.0;
    mInitDir[2] = 0.0;

    mCoordNames[0] = std::string("t");
    mCoordNames[1] = std::string("rho");
    mCoordNames[2] = std::string("phi");
    mCoordNames[3] = std::string("z");
}

// ********************************* specific protected methods *****************************
/*! Calculate metric functions at current position.
 *  \param pos : current position.
 */
void MetricHalilsoyWave::calcMetricFunc(const double* pos)
{
    double t = pos[0];
    double rho = pos[1];

    double b0 = gsl_sf_bessel_J0(rho);
    double b1 = gsl_sf_bessel_J1(rho);

    double ct = cos(t);

    fV = cosh(mAlpha) * cosh(mAlpha) * exp(-2.0 * mC * b0 * ct) + sinh(mAlpha) * sinh(mAlpha) * exp(2.0 * mC * b0 * ct);
    fA = -2.0 * mC * sinh(2.0 * mAlpha) * rho * b1 * sin(t);
    fK = 0.5 * mC * mC * (rho * rho * (b0 * b0 + b1 * b1) - 2.0 * rho * b0 * b1 * ct * ct);
}

/*! Calculate derivatives of the metric functions at current position.
 *  \param pos : current position.
 */
void MetricHalilsoyWave::calcDiffMetricFunc(const double* pos)
{
    double t = pos[0];
    double rho = pos[1];

    double b0 = gsl_sf_bessel_J0(rho);
    double b1 = gsl_sf_bessel_J1(rho);

    double ct = cos(t);
    double st = sin(t);
    double cha = cosh(mAlpha);
    double sha = sinh(mAlpha);

    double emb0 = exp(-2.0 * mC * b0 * ct);
    double epb0 = 1.0 / emb0;

    fV_t = 2.0 * cha * cha * mC * b0 * st * emb0 - 2.0 * sha * sha * mC * b0 * st * epb0;
    fV_tt = 2.0 * cha * cha * mC * b0 * ct * emb0 + 4.0 * cha * cha * mC * mC * b0 * b0 * st * st * emb0
        - 2.0 * sha * sha * mC * b0 * ct * epb0 + 4.0 * sha * sha * mC * mC * b0 * b0 * st * st * epb0;

    fV_rho = 2.0 * cha * cha * mC * b1 * ct * emb0 - 2.0 * sha * sha * mC * b1 * ct * epb0;
    fV_rhorho = 2.0 * cha * cha * mC * (b0 - b1 / rho) * ct * emb0
        + 4.0 * cha * cha * mC * mC * b1 * b1 * ct * ct * emb0 - 2.0 * sha * sha * mC * (b0 - b1 / rho) * ct * epb0
        + 4.0 * sha * sha * mC * mC * b1 * b1 * ct * ct * epb0;

    fV_trho = -2.0 * cha * cha * mC * b1 * st * emb0 + 4.0 * cha * cha * mC * mC * b0 * st * b1 * ct * emb0
        + 2.0 * sha * sha * mC * b1 * st * epb0 + 4.0 * sha * sha * mC * mC * b0 * st * b1 * ct * epb0;

    double sh2a = sinh(2.0 * mAlpha);
    fA_t = -2.0 * mC * sh2a * rho * b1 * ct;
    fA_tt = 2.0 * mC * sh2a * rho * b1 * st;

    fA_rho = -2.0 * mC * sh2a * st * rho * b0;
    fA_rhorho = 2.0 * mC * sh2a * st * (-b0 + b1 * rho);

    fA_trho = -2.0 * mC * sh2a * ct * rho * b0;

    fK_t = 2.0 * mC * mC * rho * b0 * b1 * ct * st;
    fK_tt = 2.0 * mC * mC * rho * b0 * b1 * (-1.0 + 2.0 * ct * ct);

    fK_rho = mC * mC * rho * (b0 * b0 + b1 * b1 * ct * ct - b0 * b0 * ct * ct);
    fK_rhorho = -mC * mC
        * (-b0 * b0 + 2.0 * b0 * b1 * rho + b1 * b1 * ct * ct + b0 * b0 * ct * ct - 4.0 * rho * b0 * b1 * ct * ct);

    fK_trho = -2.0 * mC * mC * ct * st * rho * (b1 * b1 - b0 * b0);
}

} // end namespace m4d
