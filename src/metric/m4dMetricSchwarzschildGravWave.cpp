// -------------------------------------------------------------------------------
/*
   m4dMetricSchwarzschildGravWave.cpp

  Copyright (c) 2017  Thomas Mueller


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

#include "m4dMetricSchwarzschildGravWave.h"
#include <cmath>

namespace m4d {


/*! Standard constructor for the SchwarzschildGravWave metric.
 *
 * \param  mass : mass of the black hole.
 */
MetricSchwarzschildGravWave::MetricSchwarzschildGravWave(double mass) {
    mMetricName  = "SchwarzschildGravWave";
    mMetricCPPfilename = "m4dMetricSchwarzschildGravWave.cpp";
    setCoordType(enum_coordinate_spherical);

    mPhysicalUnits = enum_physical_constants_geom;

    addParam("mass", mass);
    mMass = mass;
    rs = 2.0 * mMass;

    mEpsilon = 0.1;
    addParam("epsilon", mEpsilon);

    ml = 0;
    mSigma = 1.0;
    addParam("sigma", mSigma);

    /*  Only a static tetrad is defined  */
    mLocTeds.push_back(enum_nat_tetrad_static);

    setStandardValues();
}

/*!
 */
MetricSchwarzschildGravWave::~MetricSchwarzschildGravWave() {
}


// *********************************** public methods ******************************

/*! Calculate the contravariant metric components at position 'pos'.
 *
 *  \param pos : pointer to position.
 */
bool MetricSchwarzschildGravWave::calculateMetric(const double* pos) {
    double r     = pos[1];
    double theta = pos[2];
    double f = 1.0 - rs/r;
    calcPerturbations(pos);

    double c = 1.0;
    double t1 = c*c;
    double t2 = f; // f(r);
    double t5 = htt; // htt(t,r);
    double t9 = hrr; // hrr(t,r);
    double t12 = r*r;
    double t13 = hee; // hee(t,r);
    double t16 = sin(theta);
    double t17 = t16*t16;
    double t19 = hpp; // hpp(t,r);

    g_compts[0][0] = -t1*t2+t1*mEpsilon*t5;
    g_compts[0][1] = 0.0;
    g_compts[0][2] = 0.0;
    g_compts[0][3] = 0.0;
    g_compts[1][0] = 0.0;
    g_compts[1][1] = 1/t2+mEpsilon*t9;
    g_compts[1][2] = 0.0;
    g_compts[1][3] = 0.0;
    g_compts[2][0] = 0.0;
    g_compts[2][1] = 0.0;
    g_compts[2][2] = t12+mEpsilon*t13;
    g_compts[2][3] = 0.0;
    g_compts[3][0] = 0.0;
    g_compts[3][1] = 0.0;
    g_compts[3][2] = 0.0;
    g_compts[3][3] = t12*t17+mEpsilon*t19;

    return true;
}


void MetricSchwarzschildGravWave::calcLegendre(int l, double theta,
    double &Pl, double &Plt, double &Pltt, double &Plttt)
{
    if (l == 0) {
        Pl = 1.0;
        Plt = 0.0;
        Pltt = 0.0;
        Plttt = 0.0f;
    }
    else {
        double ct = cos(theta);
        double st = sin(theta);
        double x = ct;

        if (l == 1) {
            Pl = x;
            Plt = -st;
            Pltt = -ct;
            Plttt = st;
        }
        // TODO
    }
}

void MetricSchwarzschildGravWave::calcPerturbations(const double* pos) {
    double t = pos[0];
    double r = pos[1];
    double theta = pos[2];

    double M = 0.5*rs;
    double r2 = r*r;
    double r4 = r2*r2;

    double f = 1.0 - rs/r;
    double p = M - (M*M + mSigma*mSigma*r4) / (r - 2*M);
    double q = sqrt(f) / r2;

    double X = p * q;
    double Y = 3.0 * M * q;
    double Z = (r - 3.0*M) * q;
    double W = r * q;

    double Pl, Plt, Pltt, Plttt;
    calcLegendre(ml, theta, Pl, Plt, Pltt, Plttt);

    double cst = cos(mSigma * t);
    double sth = sin(theta);

    htt = -f * X * Pl * cst;
    hrr = Y / f * Pl * cst;
    hee = r2 * (Z * Pl + W * Pltt) * cst;
    hpp = r2 * sth*sth * (Z * Pl + W * Plt * cos(theta)/sth) * cst;
}


void MetricSchwarzschildGravWave::calcPerturbationsAndDiffs(const double* pos) {
    double t = pos[0];
    double r = pos[1];
    double theta = pos[2];

    double M = 0.5*rs;
    double r2 = r*r;
    double r3 = r2*r;
    double r4 = r2*r2;

    double f = 1.0 - rs/r;
    double df = rs/r2;
    double sf = sqrt(f);
    double p = M - (M*M + mSigma*mSigma*r4) / (r - 2*M);
    double q = sf / r2;

    double X = p * q;
    double Y = 3.0 * M * q;
    double Z = (r - 3.0*M) * q;
    double W = r * q;

    double Pl, Plt, Pltt, Plttt;
    calcLegendre(ml, theta, Pl, Plt, Pltt, Plttt);

    double cst = cos(mSigma * t);
    double sst = sin(mSigma * t);
    double sth = sin(theta);
    double cot = cos(theta) / sth;

    htt = -f * X * Pl * cst;
    hrr = Y / f * Pl * cst;
    hee = r2 * (Z * Pl + W * Pltt) * cst;
    hpp = r2 * sth*sth * (Z * Pl + W * Plt * cot) * cst;


    double dp = -(mSigma*mSigma*r3*(3.0*r - 4.0*rs) - M*M) / ((r-rs)*(r-rs));
    double dq = (0.5*rs/sf - 2.0*sf*r) / r4;
    double DX = dp * q + p * dq;
    double DY = 3.0 * M * dq;
    double DZ = q + (r - 3*M)*dq;
    double DW = q + r * dq;

    htt_t = f * X * Pl * sst * mSigma;
    htt_r = -df * X * Pl * cst - f * DX * Pl * cst;
    // htt_theta =

    hrr_t = -Y * Pl * sst * mSigma / f;
    hrr_r = -Y * Pl * cst * df / (f*f) + DY * Pl * cst / f;
    // hrr_theta =

    hee_t = -r2 * (Z * Pl + W * Pltt) * sst * mSigma;
    hee_r = 2*r*(Z * Pl + W * Pltt) * cst + r2 * (DZ * Pl + DW * Pltt) * cst;
    // hee_theta =

    hpp_t = -r2*sth*sth * (Z * Pl + W * Plt * cot) * sst * mSigma;
    hpp_r = 2*r*sth*sth * (Z * Pl + W * Plt * cot) * cst + r2*sth*sth*(DZ*Pl + DW*Plt*cot)*cst;
    // hpp_theta =
}


/*! Calculate the Christoffel symbols of the second kind at position 'pos'.
 *
 *  \param pos : pointer to position.
 */
bool MetricSchwarzschildGravWave::calculateChristoffels(const double* pos) {
    double r     = pos[1];
    double theta = pos[2];
    double f = 1.0 - rs/r;
    double df = rs/(r*r);
    calcPerturbationsAndDiffs(pos);

    double c = 1.0;
    double t1 = f;
    double t2 = htt; //htt(t,r);
    double t5 = 1/(t1-mEpsilon*t2);
    double t7 = htt_t; //diff(htt(t,r),t);
    double t10 = c*c;
    double t12 = df; // diff(f(r),r);
    double t13 = htt_r; // diff(htt(t,r),r);
    double t15 = t12-mEpsilon*t13;
    double t16 = hrr; //hrr(t,r);
    double t20 = 1/(1.0+mEpsilon*t16*t1);
    double t25 = t15*t5/2.0;
    double t27 = hrr_t; // diff(hrr(t,r),t);
    double t28 = mEpsilon*t27;
    double t30 = t1*t20*t28/2.0;
    double t31 = r*r;
    double t32 = hee; // hee(t,r);
    double t35 = 1/(t31+mEpsilon*t32);
    double t37 = hee_t; // diff(hee(t,r),t);
    double t39 = t35*mEpsilon*t37/2.0;
    double t40 = sin(theta);
    double t41 = t40*t40;
    double t43 = hpp; //hpp(t,r);
    double t46 = 1/(t31*t41+mEpsilon*t43);
    double t48 = hpp_t; // diff(hpp(t,r),t);
    double t50 = t46*mEpsilon*t48/2.0;
    double t52 = t5/t10;
    double t57 = hrr_r; //diff(hrr(t,r),r);
    double t59 = t1*t1;
    double t65 = hee_r; // diff(hee(t,r),r);
    double t67 = 2.0*r+mEpsilon*t65;
    double t69 = t67*t35/2.0;
    double t72 = hpp_r; //diff(hpp(t,r),r);
    double t74 = 2.0*r*t41+mEpsilon*t72;
    double t76 = t74*t46/2.0;
    double t84 = cos(theta);
    double t85 = t40*t84;
    double t86 = t46*t31*t85;

    christoffel[0][0][0] = -t5*mEpsilon*t7/2.0;
    christoffel[0][0][1] = t1*t10*t15*t20/2.0;
    christoffel[0][0][2] = 0.0;
    christoffel[0][0][3] = 0.0;
    christoffel[0][1][0] = t25;
    christoffel[0][1][1] = t30;
    christoffel[0][1][2] = 0.0;
    christoffel[0][1][3] = 0.0;
    christoffel[0][2][0] = 0.0;
    christoffel[0][2][1] = 0.0;
    christoffel[0][2][2] = t39;
    christoffel[0][2][3] = 0.0;
    christoffel[0][3][0] = 0.0;
    christoffel[0][3][1] = 0.0;
    christoffel[0][3][2] = 0.0;
    christoffel[0][3][3] = t50;
    christoffel[1][0][0] = t25;
    christoffel[1][0][1] = t30;
    christoffel[1][0][2] = 0.0;
    christoffel[1][0][3] = 0.0;
    christoffel[1][1][0] = t52*t28/2.0;
    christoffel[1][1][1] = 1/t1*t20*(-t12+mEpsilon*t57*t59)/2.0;
    christoffel[1][1][2] = 0.0;
    christoffel[1][1][3] = 0.0;
    christoffel[1][2][0] = 0.0;
    christoffel[1][2][1] = 0.0;
    christoffel[1][2][2] = t69;
    christoffel[1][2][3] = 0.0;
    christoffel[1][3][0] = 0.0;
    christoffel[1][3][1] = 0.0;
    christoffel[1][3][2] = 0.0;
    christoffel[1][3][3] = t76;
    christoffel[2][0][0] = 0.0;
    christoffel[2][0][1] = 0.0;
    christoffel[2][0][2] = t39;
    christoffel[2][0][3] = 0.0;
    christoffel[2][1][0] = 0.0;
    christoffel[2][1][1] = 0.0;
    christoffel[2][1][2] = t69;
    christoffel[2][1][3] = 0.0;
    christoffel[2][2][0] = t52*mEpsilon*t37/2.0;
    christoffel[2][2][1] = -t1*t67*t20/2.0;
    christoffel[2][2][2] = 0.0;
    christoffel[2][2][3] = 0.0;
    christoffel[2][3][0] = 0.0;
    christoffel[2][3][1] = 0.0;
    christoffel[2][3][2] = 0.0;
    christoffel[2][3][3] = t86;
    christoffel[3][0][0] = 0.0;
    christoffel[3][0][1] = 0.0;
    christoffel[3][0][2] = 0.0;
    christoffel[3][0][3] = t50;
    christoffel[3][1][0] = 0.0;
    christoffel[3][1][1] = 0.0;
    christoffel[3][1][2] = 0.0;
    christoffel[3][1][3] = t76;
    christoffel[3][2][0] = 0.0;
    christoffel[3][2][1] = 0.0;
    christoffel[3][2][2] = 0.0;
    christoffel[3][2][3] = t86;
    christoffel[3][3][0] = t52*mEpsilon*t48/2.0;
    christoffel[3][3][1] = -t1*t74*t20/2.0;
    christoffel[3][3][2] = -t35*t31*t85;
    christoffel[3][3][3] = 0.0;

    return true;
}

/*! Calculate Jacobi matrix.
 *
 *  \param pos : pointer to position.
 */

/*
bool MetricSchwarzschildGravWave::calculateChrisD(const double* pos) {
    double r     = pos[1];
    double theta = pos[2];
    return true;
}
*/


/*! Transform local 4-direction to coordinate 4-direction.
 *
 *  \param  pos  :  pointer to position array.
 *  \param  ldir :  pointer to local direction array.
 *  \param  dir  :  pointer to calculated coordinate direction array.
 *  \param  type :  type of tetrad.
 */
void MetricSchwarzschildGravWave::localToCoord(const double* pos, const double* ldir, double* dir,
                                       enum_nat_tetrad_type)
{
    calculateMetric(pos);

    double a = sqrt(-g_compts[0][0]);
    double b = sqrt(g_compts[1][1]);
    double c = sqrt(g_compts[2][2]);
    double d = sqrt(g_compts[3][3]);

    dir[0] = ldir[0] / a;
    dir[1] = ldir[1] / b;
    dir[2] = ldir[2] / c;
    dir[3] = ldir[3] / d;
}

/*! Transform coordinate 4-direction to local 4-direction.
 *
 *  \param  pos  :  pointer to position array.
 *  \param  cdir :  pointer to coordinate direction.
 *  \param  ldir :  pointer to calculated local direction array.
 *  \param  type :  type of tetrad.
 */
void MetricSchwarzschildGravWave::coordToLocal(const double* pos, const double* cdir, double* ldir,
                                       enum_nat_tetrad_type) {
    calculateMetric(pos);

    double a = sqrt(-g_compts[0][0]);
    double b = sqrt(g_compts[1][1]);
    double c = sqrt(g_compts[2][2]);
    double d = sqrt(g_compts[3][3]);

    ldir[0] = cdir[0] * a;
    ldir[1] = cdir[1] * b;
    ldir[2] = cdir[2] * c;
    ldir[3] = cdir[3] * d;
}


/*! Test break condition.
 *
 *  \param pos    : pointer to position array.
 *  \return true  : radial position r < 0.0 or  r^2<=(1.0+eps)*rs^2.
 *  \return false : position is valid.
 */
bool MetricSchwarzschildGravWave::breakCondition(const double* pos) {
    bool br = false;

    if ((pos[1] < 0.0) || (pos[1]*pos[1] <= (1.0 + 1e-6)*rs * rs)) {
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
 *  However, take care of the limited double precision.
 *  \param  y[]   : pointer to position and direction coordinates.
 *  \param  kappa : timelike (-1.0), lightlike (0.0).
 *  \return double : sum.
 */
double MetricSchwarzschildGravWave::testConstraint(const double y[], const double kappa) {       
    double dt  = y[4];
    double dr  = y[5];
    double dth = y[6];
    double dph = y[7];

    calculateMetric(y);
    double a = g_compts[0][0];
    double b = g_compts[1][1];
    double c = g_compts[2][2];
    double d = g_compts[3][3];

    double sum = -kappa;
    sum += a * dt*dt + b * dr*dr + c * dth*dth + d * dph*dph;
    return sum;
}


/*! Set parameter 'pName' to 'val'.
 *
 *  Set 'mass' parameter and adjust Schwarzschild radius  rs=2GM/c^2.
 */
bool MetricSchwarzschildGravWave::setParam(std::string pName, double val) {
    Metric::setParam(pName, val);
    if (pName == "mass") {
        mMass = val;
        rs = 2*mMass;
    } else if (pName == "epsilon") {
        mEpsilon = val;
    } else if (pName == "sigma") {
        mSigma = val;
    }
    return true;
}

/*! Generate report.
 * \param pos : initial position.
 * \param cdir : initial coordinate direction.
 * \param text : reference to report text.
 */
bool MetricSchwarzschildGravWave::report(const vec4 pos, const vec4 cdir, std::string &text) {
    std::stringstream ss;
    ss << "Report for SchwarzschildGravWave metric\n\tcoordinates : (t,r,theta,phi)\n";
    ss << "---------------------------------------------------------------\n";
    ss << "  physical units ................................. yes\n";
    ss.precision(DEF_FIXED_REPORT_PRECISION);
    ss.setf(std::ios::fixed);
    ss << "  Schwarzschild radius ........... r_s = 2GM/c^2 = " << rs << std::endl;

    text = ss.str();
    return true;
}


// ********************************* protected methods *****************************
/*!
 */
void MetricSchwarzschildGravWave::setStandardValues() {
    mInitPos[0] = 0.0;
    mInitPos[1] = 3.0 * rs;
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
