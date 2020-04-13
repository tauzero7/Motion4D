// -------------------------------------------------------------------------------
/*
   m4dMetricAlcubierreSimple.cpp

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

#include "m4dMetricAlcubierreSimple.h"

namespace m4d {

/*! Standard constructor for the KerrBL metric.
 *
 * \param  R     : size of warp bubble.
 * \param  vs    : velocity of warp bubble.
 */
MetricAlcubierreSimple::MetricAlcubierreSimple(double R, double dR, double vs)
{
    mMetricName = "AlcubierreWarpSimple";
    mMetricCPPfilename = "m4dMetricAlcubierreSimple.cpp";
    setCoordType(enum_coordinate_cartesian);

    mPhysicalUnits = enum_physical_constants_geom;
    mSpeedOfLight = 1.0;
    mGravConstant = 1.0;

    mR = R;
    mDR = dR;
    mvs = vs;

    addParam("r", R);
    addParam("dr", dR);
    addParam("vs", vs);

    // mDrawTypes.push_back(enum_draw_twoplusone);

    setStandardValues();

    mLocTeds.push_back(enum_nat_tetrad_comoving);
    mLocTeds.push_back(enum_nat_tetrad_static);
}

MetricAlcubierreSimple::~MetricAlcubierreSimple() {}

// *********************************** public methods ******************************
/*! Calculate the covariant metric components at position 'pos'.
 *
 *  \param pos : pointer to position.
 */
bool MetricAlcubierreSimple::calculateMetric(const double* pos)
{
    double c = mSpeedOfLight;
    double vs = mvs;

    double f = calcF(pos);

    g_compts[0][0] = -pow(c, 2) + pow(vs, 2) * pow(f, 2);
    g_compts[0][1] = -vs * f;
    g_compts[0][2] = 0;
    g_compts[0][3] = 0;
    g_compts[1][0] = -vs * f;
    g_compts[1][1] = 1;
    g_compts[1][2] = 0;
    g_compts[1][3] = 0;
    g_compts[2][0] = 0;
    g_compts[2][1] = 0;
    g_compts[2][2] = 1;
    g_compts[2][3] = 0;
    g_compts[3][0] = 0;
    g_compts[3][1] = 0;
    g_compts[3][2] = 0;
    g_compts[3][3] = 1;

    return true;
}

/*! Calculate the Christoffel symbols of the second kind at position 'pos'.
 *
 *  \param pos : pointer to position.
 */
bool MetricAlcubierreSimple::calculateChristoffels(const double* pos)
{
    double c = mSpeedOfLight;
    double vs = mvs;

    double f = calcF(pos);

    double ft, fx, fy, fz;
    calcDF(pos, ft, fx, fy, fz);

    christoffel[0][0][0] = pow(vs, 3) * pow(f, 2) * fx / pow(c, 2);
    christoffel[0][0][1] = vs * (-pow(c, 2) * vs * f * fx - pow(c, 2) * ft + pow(vs, 3) * pow(f, 3) * fx) / pow(c, 2);
    christoffel[0][0][2] = -pow(vs, 2) * f * fy;
    christoffel[0][0][3] = -pow(vs, 2) * f * fz;
    christoffel[0][1][0] = -pow(vs, 2) * f * fx / pow(c, 2);
    christoffel[0][1][1] = -pow(vs, 3) * pow(f, 2) * fx / pow(c, 2);
    christoffel[0][1][2] = vs * fy / 2;
    christoffel[0][1][3] = vs * fz / 2;
    christoffel[0][2][0] = -pow(vs, 2) * f * fy / (2 * pow(c, 2));
    christoffel[0][2][1] = vs * (-pow(c, 2) - pow(vs, 2) * pow(f, 2)) * fy / (2 * pow(c, 2));
    christoffel[0][2][2] = 0;
    christoffel[0][2][3] = 0;
    christoffel[0][3][0] = -pow(vs, 2) * f * fz / (2 * pow(c, 2));
    christoffel[0][3][1] = vs * (-pow(c, 2) - pow(vs, 2) * pow(f, 2)) * fz / (2 * pow(c, 2));
    christoffel[0][3][2] = 0;
    christoffel[0][3][3] = 0;
    christoffel[1][0][0] = -pow(vs, 2) * f * fx / pow(c, 2);
    christoffel[1][0][1] = -pow(vs, 3) * pow(f, 2) * fx / pow(c, 2);
    christoffel[1][0][2] = vs * fy / 2;
    christoffel[1][0][3] = vs * fz / 2;
    christoffel[1][1][0] = vs * fx / pow(c, 2);
    christoffel[1][1][1] = pow(vs, 2) * f * fx / pow(c, 2);
    christoffel[1][1][2] = 0;
    christoffel[1][1][3] = 0;
    christoffel[1][2][0] = vs * fy / (2 * pow(c, 2));
    christoffel[1][2][1] = pow(vs, 2) * f * fy / (2 * pow(c, 2));
    christoffel[1][2][2] = 0;
    christoffel[1][2][3] = 0;
    christoffel[1][3][0] = vs * fz / (2 * pow(c, 2));
    christoffel[1][3][1] = pow(vs, 2) * f * fz / (2 * pow(c, 2));
    christoffel[1][3][2] = 0;
    christoffel[1][3][3] = 0;
    christoffel[2][0][0] = -pow(vs, 2) * f * fy / (2 * pow(c, 2));
    christoffel[2][0][1] = vs * (-pow(c, 2) - pow(vs, 2) * pow(f, 2)) * fy / (2 * pow(c, 2));
    christoffel[2][0][2] = 0;
    christoffel[2][0][3] = 0;
    christoffel[2][1][0] = vs * fy / (2 * pow(c, 2));
    christoffel[2][1][1] = pow(vs, 2) * f * fy / (2 * pow(c, 2));
    christoffel[2][1][2] = 0;
    christoffel[2][1][3] = 0;
    christoffel[2][2][0] = 0;
    christoffel[2][2][1] = 0;
    christoffel[2][2][2] = 0;
    christoffel[2][2][3] = 0;
    christoffel[2][3][0] = 0;
    christoffel[2][3][1] = 0;
    christoffel[2][3][2] = 0;
    christoffel[2][3][3] = 0;
    christoffel[3][0][0] = -pow(vs, 2) * f * fz / (2 * pow(c, 2));
    christoffel[3][0][1] = vs * (-pow(c, 2) - pow(vs, 2) * pow(f, 2)) * fz / (2 * pow(c, 2));
    christoffel[3][0][2] = 0;
    christoffel[3][0][3] = 0;
    christoffel[3][1][0] = vs * fz / (2 * pow(c, 2));
    christoffel[3][1][1] = pow(vs, 2) * f * fz / (2 * pow(c, 2));
    christoffel[3][1][2] = 0;
    christoffel[3][1][3] = 0;
    christoffel[3][2][0] = 0;
    christoffel[3][2][1] = 0;
    christoffel[3][2][2] = 0;
    christoffel[3][2][3] = 0;
    christoffel[3][3][0] = 0;
    christoffel[3][3][1] = 0;
    christoffel[3][3][2] = 0;
    christoffel[3][3][3] = 0;

    return true;
}

/*! Calculate Jacobi matrix.
 *
 *  \param pos : pointer to position.
 */
bool MetricAlcubierreSimple::calculateChrisD(const double* pos)
{
    double ft, fx, fy, fz;
    calcDF(pos, ft, fx, fy, fz);

    double ftt, ftx, fty, ftz, fxx, fxy, fxz, fyy, fyz, fzz;
    calcD2F(pos, ftt, ftx, fty, ftz, fxx, fxy, fxz, fyy, fyz, fzz);

    return false;
}

/*! Calculate Riemann tensor R^a_bcd
 * \param pos : pointer to coordinate position where the Riemann tensor have to be evaluated.
 * \return true : successfull
 */
bool MetricAlcubierreSimple::calculateRiemann(const double* pos)
{
    double ft, fx, fy, fz;
    calcDF(pos, ft, fx, fy, fz);

    double ftt, ftx, fty, ftz, fxx, fxy, fxz, fyy, fyz, fzz;
    calcD2F(pos, ftt, ftx, fty, ftz, fxx, fxy, fxz, fyy, fyz, fzz);

    return false;
}

/*! Transform local 4-direction to coordinate 4-direction.
 *
 *  \param  pos  :  pointer to position array.
 *  \param  ldir :  pointer to local direction array.
 *  \param  dir  :  pointer to calculated coordinate direction array.
 *  \param  type :  type of tetrad.
 */
void MetricAlcubierreSimple::localToCoord(const double* pos, const double* ldir, double* dir, enum_nat_tetrad_type type)
{
    double f = calcF(pos);
    double c = mSpeedOfLight;

    if (type == enum_nat_tetrad_comoving) {
        dir[0] = ldir[0] / c;
        dir[1] = ldir[0] * mvs * f / c + ldir[1];
        dir[2] = ldir[2];
        dir[3] = ldir[3];
    }
    else {
        double w = sqrt(c * c - mvs * mvs * f * f);

        dir[0] = (ldir[0] - mvs * f / c * ldir[1]) / w;
        dir[1] = w / c * ldir[1];
        dir[2] = ldir[2];
        dir[3] = ldir[3];
    }
}

/*! Transform coordinate 4-direction to local 4-direction.
 *
 *  \param  pos  :  pointer to position array.
 *  \param  cdir :  pointer to coordinate direction.
 *  \param  ldir :  pointer to calculated local direction array.
 *  \param  type :  type of tetrad.
 */
void MetricAlcubierreSimple::coordToLocal(
    const double* pos, const double* cdir, double* ldir, enum_nat_tetrad_type type)
{
    double f = calcF(pos);
    double c = mSpeedOfLight;

    if (type == enum_nat_tetrad_comoving) {
        ldir[0] = c * cdir[0];
        ldir[1] = cdir[1] - mvs * f * cdir[0];
        ldir[2] = cdir[2];
        ldir[3] = cdir[3];
    }
    else {
        double w = sqrt(c * c - mvs * mvs * f * f);

        ldir[1] = c / w * cdir[1];
        ldir[0] = w * cdir[0] - ldir[1] * mvs * f / c;
        ldir[2] = cdir[2];
        ldir[3] = cdir[3];
    }
}

/*!
 *  \param pos  :  position.
 *  \return true  : radial position r < 0.0 or ...
 *  \return false : position is valid.
 */
bool MetricAlcubierreSimple::breakCondition(const double*)
{
    bool br = false;

    return br;
}

// Calculate right hand side of the geodesic equation in first order form.
//
// \param  y[]   : pointer to position and direction coordinates.
//  \param  dydx[] : pointer to right side of geodesic equation.
//
/*
bool MetricAlcubierreSimple::calcDerivs ( const double y[], double dydx[] )
{
    return false;
    dydx[0] = y[4];
    dydx[1] = y[5];
    dydx[2] = y[6];
    dydx[3] = y[7];

    double f = calcF(y);
    double f2 = f*f;

    double c = mSpeedOfLight;
    double edc2 = 1.0/(c*c);
    //fprintf(stderr,"hier...%f %f\n",y[5]-mvs*y[4],f);
    double v = mvs;
    double v2 = v*v;
    double v3 = v2*v;

    double ft,fx,fy,fz;
    calcDF(y,ft,fx,fy,fz);

    return true;
}
*/

/*! Tests whether the constraint equation is fulfilled.
 *
 *  The constraint equation for lightlike and timelike geodesics reads:
 \verbatim
     sum = g_{\mu\nu} dot(x)^{\mu} dot(x)^{\nu} - kappa c^2 = 0.
 \endverbatim
 *  \param  y[]   : pointer to position and direction coordinates.
 *  \param  kappa : timelike (-1.0), lightlike (0.0).
 *  \return sum.
 */
double MetricAlcubierreSimple::testConstraint(const double y[], const double kappa)
{
    double c = mSpeedOfLight;
    double f = calcF(y);

    double sum = -kappa;
    sum += -c * c * y[4] * y[4] + pow(y[5] - mvs * f * y[4], 2.0) + y[6] * y[6] + y[7] * y[7];

    // double k = -c * c * y[4] + mvs * (1 - f) * (y[5] - mvs * f * y[4]);
    // std::cerr << k << std::endl;

    double dxdt = y[5] / y[4];
    double dydt = y[6] / y[4];
    double w = (dxdt - mvs * f) * (dxdt - mvs * f) + dydt * dydt;
    std::cerr << y[4] << " " << y[5] << " " << y[6] << " " << w << std::endl;

    // double A = 1.0-mvs*mvs*(1.0-f)*(1.0-f);
    // fprintf(stderr,"%e %e %e %e  %e %e %e %e\n",y[4],y[5],y[6],y[7],f,A,mvs*f+1,mvs*f-1);
    return sum;
}

bool MetricAlcubierreSimple::resize(double* y, double kappa, double factor)
{
    double a = 1 / fabs(y[4]);
    y[4] *= a;
    y[5] *= a;
    y[6] *= a;
    y[7] *= a;
    return true;
}

/*! Set parameter 'pName' to 'val'.
 *
 *  Set 'sigma', 'R', and 'vs' parameters.
 */
bool MetricAlcubierreSimple::setParam(const char* pName, double val)
{
    Metric::setParam(pName, val);
    if (strcmp(pName, "r") == 0) {
        mR = val;
    }
    else if (strcmp(pName, "dr") == 0) {
        mDR = val;
    }
    else if (strcmp(pName, "vs") == 0) {
        mvs = val;
    }
    return true;
}

/*! Transform point p to 2+1 coordinates.
 *
 *  \param  p  : point in proper metric coordinates.
 *  \param  cp : reference to transformed point.
 *  \return true : success.
 */
bool MetricAlcubierreSimple::transToTwoPlusOne(vec4 p, vec4& cp)
{
    cp = vec4(p[0], p[1], p[2], p[0]);
    return true;
}

/*! Generate report.
 */
bool MetricAlcubierreSimple::report(const vec4, const vec4, char*& text)
{
    std::stringstream ss;
    ss << "Report for AlcubierreSimple metric\n\tcoordinate : (t,x,y,z)\n";
    ss << "---------------------------------------------------------------\n";
    ss << "  physical units ................................. no\n";
    ss.precision(DEF_FIXED_REPORT_PRECISION);
    ss.setf(std::ios::fixed);

    text = new char[ss.str().length() + 2];
    return CopyString(ss.str().c_str(), text);
}

// *************************** specific  public methods ****************************

// ********************************* protected methods *****************************
/*!
 */
void MetricAlcubierreSimple::setStandardValues()
{
    mInitPos[0] = 0.0;
    mInitPos[1] = 0.0;
    mInitPos[2] = -10.0;
    mInitPos[3] = 0.0;
    mInitDir[0] = 1.0;
    mInitDir[1] = 0.0;
    mInitDir[2] = 0.0;

    mCoordNames[0] = std::string("t");
    mCoordNames[1] = std::string("x");
    mCoordNames[2] = std::string("y");
    mCoordNames[3] = std::string("z");
}

/*! Calculate rs function
 *  \param  pos : pointer to position.
 */
double MetricAlcubierreSimple::calcRs(const double* pos)
{
    double t = pos[0];
    double x = pos[1] - mvs * t;
    double y = pos[2];
    double z = pos[3];

    return sqrt(x * x + y * y + z * z);
}

double MetricAlcubierreSimple::calcF(const double* pos)
{
    double rs = calcRs(pos);

#if 0
    if (rs <= mR) {
        return 1.0 - pow(rs / mR, 4.0);
    }
    return 0.0;
#else
    double R1 = mR - 0.5 * mDR;
    double R2 = mR + 0.5 * mDR;
    if (rs <= R1) {
        return 1.0;
    }

    if (rs >= R1 && rs < R2) {
        double w = (rs - R1) / (R2 - R1);
        double y = 1 - w * w;
        return y * y;
    }

    return 0.0;
#endif
}

void MetricAlcubierreSimple::calcDF(const double* pos, double& ft, double& fx, double& fy, double& fz)
{
    double rs = calcRs(pos);

    double t = pos[0];
    double x = pos[1] - mvs * t;
    double y = pos[2];
    double z = pos[3];

#if 0
    double dfdr = 0.0;
    if (rs <= mR) {
        dfdr = -4.0 * pow(rs / mR, 3.0) / mR;
    }
    double df = dfdr / rs;

    // ft = df/dr * dr/dt ...
    ft = -mvs * x * df;
    fx = x * df;
    fy = y * df;
    fz = z * df;
#else
    double R1 = mR - 0.5 * mDR;
    double R2 = mR + 0.5 * mDR;

    double dfdr = 0.0;

    if (rs >= R1 && rs < R2) {
        double w = (rs - R1) / (R2 - R1);
        dfdr = -4.0 * (1 - w * w) * w / (R2 - R1);
    }

    ft = dfdr * x / rs * (-mvs);
    fx = dfdr * x / rs;
    fy = dfdr * y / rs;
    fz = dfdr * z / rs;
#endif
}

/*
void MetricAlcubierreSimple::calcD2F ( const double* pos, double &ftt, double &ftx, double &fty, double &ftz,
                              double &fxx, double &fxy, double &fxz, double &fyy,
                              double &fyz, double &fzz )
*/
void MetricAlcubierreSimple::calcD2F(
    const double*, double&, double&, double&, double&, double&, double&, double&, double&, double&, double&)
{
    /*
    double t = pos[0];
    double x = pos[1];
    double y = pos[2];
    double z = pos[3];

    double vs = mvs;
    double R = mR;
    */
    // TODO
}

} // end namespace m4d
