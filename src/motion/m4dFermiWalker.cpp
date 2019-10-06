// -------------------------------------------------------------------------------
/*
    m4dFermiWalker.cpp

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

#include "m4dFermiWalker.h"

#define DEF_DIFF_H 1e-7
#define DEF_EDZDIFF_H 5e6 //  1/(2*DEF_DIFF_H)
#define DEF_EDDIFF2_H 1e14 //  1/DEF_DIFF_H^2

namespace m4d {

/*! Standard constructor for Fermi-Walker motion.
 *
 *  \param metric : pointer to metric.
 */
FermiWalker::FermiWalker(Metric* metric)
    : Motion(metric)
{
    for (int i = 0; i < 4; i++) {
        mPropAcc[i] = mInitVel[i] = mCurrAcc[i] = 0.0;
    }
    mCalcWithWorldline = false;
    x_tau = nullptr;
    u_tau = nullptr;
    a_tau = nullptr;
    mParams = nullptr;
    mLambda = 0.0;
}

FermiWalker::~FermiWalker() {}

// *********************************** public methods ******************************
/*! Set acceleration with respect to local tetrad.
 *
 *  \param a1 : proper acceleration in e1-direction.
 *  \param a2 : proper acceleration in e2-direction.
 *  \param a3 : proper acceleration in e3-direction.
 */
void FermiWalker::setCurrPropAccel(double a1, double a2, double a3)
{
    mPropAcc[1] = a1;
    mPropAcc[2] = a2;
    mPropAcc[3] = a3;
}

/*! Get acceleration with respect to local tetrad.
 *
 *  \param a1 : reference to proper acceleration in e1-direction.
 *  \param a2 : reference to proper acceleration in e2-direction.
 *  \param a3 : reference to proper acceleration in e3-direction.
 */
void FermiWalker::getCurrPropAccel(double& a1, double& a2, double& a3)
{
    a1 = mPropAcc[1];
    a2 = mPropAcc[2];
    a3 = mPropAcc[3];
}

/*! Set initial velocity with respect to natural local tetrad.
 *
 *  \param fm : time direction (>0 : future, <0 : past).
 *  \param v  : absolute value of velocity with respect to natural local tetrad.
 *  \param theta :  angle
 *  \param phi   :  angle
 *  \param type  : type of natural local tetrad.
 *  \sa enum_nat_tetrad_type.
 */
bool FermiWalker::setInitialVelocity(double fm, double v, double theta, double phi, enum_nat_tetrad_type type)
{
    if (fabs(v) >= 1.0 || fm == 0.0) {
        return false;
    }

    double gamma = 1.0 / sqrt(1.0 - v * v);
    mInitVel[0] = fm / fabs(fm) * gamma;
    mInitVel[1] = gamma * v * sin(theta) * cos(phi);
    mInitVel[2] = gamma * v * sin(theta) * sin(phi);
    mInitVel[3] = gamma * v * cos(theta);

    mVel = v;
    mTheta = theta;
    mPhi = phi;

    double e0[4];
    mMetric->localToCoord(y, mInitVel, e0, type);
    for (int i = 0; i < 4; i++) {
        y[4 + i] = e0[i];
        y[8 + i] = e0[i];
    }
    return true;
}

/*! Get initial velocity with respect to natural local tetrad.
 *
 */
// TODO Parallel gibt es problem !!
void FermiWalker::getInitialVelocity(double& v, double& theta, double& phi)
{
    v = mVel;
    theta = mTheta;
    phi = mPhi;
}

/*! Calculate motion with at most 'maxNumPoints' points.
 *
 *  \param initPos : initial position.
 *  \param fm : time direction (fm>0: future, fm<0: past).
 *  \param v :  velocity.
 *  \param theta_v : theta-direction of velocity-
 *  \param phi_v : phi-direction of velocity.
 *  \param a :  acceleration.
 *  \param theta_a :  theta-direction of acceleration.
 *  \param phi_a :  phi-direction of acceleration.
 *  \param e0 : local tetrad vector.
 *  \param e1 : local tetrad vector.
 *  \param e2 : local tetrad vector.
 *  \param e3 : local tetrad vector.
 *  \param maxNumPoints : maximum number of points.
 *  \param points :  reference to points.
 *  \param base0 : reference to base vector.
 *  \param base1 : reference to base vector.
 *  \param base2 : reference to base vector.
 *  \param base3 : reference to base vector.
 */
enum_break_condition FermiWalker::calculateMotion(const vec4 initPos, double fm, double v, double theta_v, double phi_v,
    double a, double theta_a, double phi_a, const vec4 e0, const vec4 e1, const vec4 e2, const vec4 e3,
    const int maxNumPoints, std::vector<vec4>& points, std::vector<vec4>& base0, std::vector<vec4>& base1,
    std::vector<vec4>& base2, std::vector<vec4>& base3)
{
    if (fm == 0.0) {
        return enum_break_constraint;
    }
    fm = fm / fabs(fm); // normalize to 1.0

    if (!points.empty()) {
        points.clear();
    }
    if (!base0.empty()) {
        base0.clear();
    }
    if (!base1.empty()) {
        base1.clear();
    }
    if (!base2.empty()) {
        base2.clear();
    }
    if (!base3.empty()) {
        base3.clear();
    }

    setInitialPosition(initPos);
    setInitialTetrad(e0, e1, e2, e3);
    setInitialVelocity(fm, v, theta_v, phi_v);
    setCurrPropAccel(a, theta_a, phi_a);

    points.push_back(initPos);
    base0.push_back(e0);
    base1.push_back(e1);
    base2.push_back(e2);
    base3.push_back(e3);

    enum_break_condition breakType = enum_break_none;
    int i = 0;

    while (i < maxNumPoints && (breakType == enum_break_none)) {
        breakType = nextStep();
        points.push_back(getPosition());
        base0.push_back(getE(0));
        base1.push_back(getE(1));
        base2.push_back(getE(2));
        base3.push_back(getE(3));
        i++;
    }
    if ((int)points.size() >= maxNumPoints) {
        breakType = enum_break_num_exceed;
    }

    return breakType;
}

/*! Set 'true' if the Fermi-Walker transport is along a given worldline.
 *
 *  \param withworldline : (true:) calculate Fermi-Walker transport with a predefined wordline, (false:) acceleration
 * with respect to the local tetrad.
 */
void FermiWalker::setCalcWithWorldline(bool withworldline)
{
    mCalcWithWorldline = withworldline;
}

/*! Set the worldline function.
 *
 *  \param func : function pointer to the worldline x=x(tau).
 */
void FermiWalker::set_x_tau(vec4 (*func)(double, void*))
{
    x_tau = func;
}

/*! Set the tangent to the worldline function.
 *
 *  \param func : function pointer to the four-velocity of the worldline x=x(tau).
 */
void FermiWalker::set_u_tau(vec4 (*func)(double, void*))
{
    u_tau = func;
}

/*!  Set the acceleration to the worldline function.
 *
 *  \param func : function pointer to the four-acceleration of the worldline x=x(tau).
 */
void FermiWalker::set_a_tau(vec4 (*func)(double, void*))
{
    a_tau = func;
}

/*!  Set the parameters to the worldline function.
 *
 *  \param params : pointer to parameters.
 */
void FermiWalker::set_params(void* params)
{
    mParams = params;
}

/*!  Evaluate the worldline for proper time tau.
 *
 *  \param tau : proper time.
 *  \param x   : reference to current position.
 *  \return true : worldline x=x(tau) is given.
 *  \return false : worlined is not defined.
 */
bool FermiWalker::get_x_tau(double tau, vec4& x)
{
    if (x_tau == nullptr) {
        return false;
    }

    x = x_tau(tau, mParams);
    return true;
}

/*!   Evaluate tangent to the worldline for proper time tau.
 *
 *  \param tau : proper time.
 *  \param u   : reference to current four-velocity.
 *  \return true : worldline x=x(tau) is given.
 *  \return false : worlined is not defined.
 */
bool FermiWalker::get_u_tau(double tau, vec4& u)
{
    if (x_tau == nullptr) {
        return false;
    }

    if (u_tau == nullptr) {
        vec4 xp, xm;
        get_x_tau(tau + DEF_DIFF_H, xp);
        get_x_tau(tau - DEF_DIFF_H, xm);

        u = (xp - xm) * DEF_EDZDIFF_H;
    }
    else {
        u = u_tau(tau, mParams);
    }

    return true;
}

/*!  Evaluate acceleration to the worldline for proper time tau.
 *
 *  \param tau : proper time.
 *  \param a   : reference to current four-acceleration.
 *  \return true : worldline x=x(tau) is given.
 *  \return false : worlined is not defined.
 */
bool FermiWalker::get_a_tau(double tau, vec4& a)
{
    if (x_tau == nullptr) {
        return false;
    }

    if (a_tau == nullptr) {
        // Calculate christoffel symbols at current position x(tau).
        vec4 pos;
        get_x_tau(tau, pos);
        mMetric->calculateChristoffels(pos);

        if (u_tau == nullptr) {
            vec4 xp, xm, x;
            get_x_tau(tau + DEF_DIFF_H, xp);
            get_x_tau(tau - DEF_DIFF_H, xm);
            get_x_tau(tau, x);
            a = (xp + xm - 2.0 * x) * DEF_EDDIFF2_H;

            for (int mu = 0; mu < 4; mu++)
                for (int i = 0; i < 4; i++)
                    for (int k = 0; k < 4; k++) {
                        a[mu] += mMetric->getChristoffel(i, k, mu) * ((xp[i] - xm[i]) * DEF_EDZDIFF_H)
                            * ((xp[k] - xm[k]) * DEF_EDZDIFF_H);
                    }
        }
        else {
            vec4 up, um, u;
            up = u_tau(tau + DEF_DIFF_H, mParams);
            um = u_tau(tau - DEF_DIFF_H, mParams);
            u = u_tau(tau, mParams);
            a = (up - um) * DEF_EDZDIFF_H;
            for (int mu = 0; mu < 4; mu++)
                for (int i = 0; i < 4; i++)
                    for (int k = 0; k < 4; k++) {
                        a[mu] += mMetric->getChristoffel(i, k, mu) * u[i] * u[k];
                    }
        }
    }
    else {
        a = a_tau(tau, mParams);
    }
    return true;
}

/*! Initialize worldline.
 *
 *  \param tauStart : initial proper time.
 */
bool FermiWalker::initWorldline(double tauStart)
{
    if (x_tau == nullptr) {
        return false;
    }

    mLambda = tauStart;
    return updateWorldline(mLambda);
}

/*!  Update worldline to current time tau.
 *
 * \param tau : current proper time.
 */
bool FermiWalker::updateWorldline(double)
{
    if (x_tau == nullptr) {
        return false;
    }

    vec4 pos = x_tau(mLambda, mParams);
    vec4 vel;
    get_u_tau(mLambda, vel);

    for (int mu = 0; mu < 4; mu++) {
        y[mu] = pos[mu];
        y[mu + 4] = vel[mu];
        y[mu + 8] = vel[mu];
    }
    return true;
}

/*! Calculate the next step of the Fermi-Walker transport.
 *
 *  \return enum_break_condition.
 */
enum_break_condition FermiWalker::nextStep()
{
    if (mMetric->breakCondition(&y[0])) {
        return enum_break_cond;
    }
    int i;
    double yn[DEF_MAX_YS], dydx[DEF_MAX_YS], k1[DEF_MAX_YS], k2[DEF_MAX_YS], k3[DEF_MAX_YS], k4[DEF_MAX_YS];

    vec4 pos, vel;

    //  pos = x_tau(mLambda,mParams);
    //  get_u_tau(mLambda,vel);

    calcDerivs(y, dydx);
    for (i = 0; i < DEF_MAX_YS; i++) {
        k1[i] = mLambdaStep * dydx[i];
        yn[i] = y[i] + 0.5 * k1[i];
    }

    calcDerivs(yn, dydx);
    for (i = 0; i < DEF_MAX_YS; i++) {
        k2[i] = mLambdaStep * dydx[i];
        yn[i] = y[i] + 0.5 * k2[i];
    }

    calcDerivs(yn, dydx);
    for (i = 0; i < DEF_MAX_YS; i++) {
        k3[i] = mLambdaStep * dydx[i];
        yn[i] = y[i] + k3[i];
    }

    calcDerivs(yn, dydx);
    for (i = 0; i < DEF_MAX_YS; i++) {
        k4[i] = mLambdaStep * dydx[i];
        yn[i] = y[i] + 1.0 / 6.0 * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);
        y[i] = yn[i];
    }

    if (fabs(mMetric->testConstraint(y, -1.0)) > mConstraintEpsilon) {
        return enum_break_constraint;
    }

    mLambda += mLambdaStep;
    return enum_break_none;
}

/*! Calculate the next step of the Fermi-Walker transport along the worldline.
 *
 *  \return enum_break_condition.
 */
enum_break_condition FermiWalker::nextStepWL()
{
    if (mMetric->breakCondition(&y[0])) {
        return enum_break_cond;
    }

    if (x_tau == nullptr) {
        return enum_break_other;
    }

    int i;
    double yn[DEF_MAX_YS], dydx[DEF_MAX_YS], k1[DEF_MAX_YS], k2[DEF_MAX_YS], k3[DEF_MAX_YS], k4[DEF_MAX_YS];

    updateWorldline(mLambda);
    get_a_tau(mLambda, mCurrAcc);

    calcDerivsWL(y, dydx);
    for (i = 0; i < DEF_MAX_YS; i++) {
        k1[i] = mLambdaStep * dydx[i];
        yn[i] = y[i] + 0.5 * k1[i];
    }

    calcDerivsWL(yn, dydx);
    for (i = 0; i < DEF_MAX_YS; i++) {
        k2[i] = mLambdaStep * dydx[i];
        yn[i] = y[i] + 0.5 * k2[i];
    }

    calcDerivsWL(yn, dydx);
    for (i = 0; i < DEF_MAX_YS; i++) {
        k3[i] = mLambdaStep * dydx[i];
        yn[i] = y[i] + k3[i];
    }

    calcDerivsWL(yn, dydx);
    for (i = 0; i < DEF_MAX_YS; i++) {
        k4[i] = mLambdaStep * dydx[i];
        yn[i] = y[i] + 1.0 / 6.0 * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);
        y[i] = yn[i];
    }

    mLambda += mLambdaStep;
    return enum_break_none;
}

// ********************************* protected methods *****************************
/*! Calculate right hand side of Fermi-Walker transport.
 *
 * \param y[] : pointer to y.
 * \param dydx[] : pointer to right hand side of geodesic equation.
 */
bool FermiWalker::calcDerivs(const double y[], double dydx[])
{
    int mu, j, k, l, bidx;

    if (!mMetric->calcDerivsFW(mPropAcc, y, dydx)) {
        double c = mMetric->speed_of_light();

        mMetric->calculateChristoffels(y);

        for (mu = 0; mu < 4; mu++) {
            dydx[mu] = y[mu + 4];

            for (j = 0; j < 4; j++) {
                bidx = mu + 4 * (j + 2);
                dydx[bidx] = 0.0;

                for (k = 0; k < 4; k++) {
                    for (l = 0; l < 4; l++) {
                        dydx[bidx] -= mMetric->getChristoffel(k, l, mu) * y[k + 8] * y[8 + 4 * j + l];
                    }
                    dydx[bidx]
                        += (eta(j, k) * mPropAcc[k] * y[mu + 8] - eta(0, j) * mPropAcc[k] * y[8 + 4 * k + mu]) / c;
                }
            }
            dydx[mu + 4] = dydx[mu + 8];
        }
    }
    return true;
}

/*! Calculate right hand side of Fermi-Walker transport for a given worldline.
 *
 * \param y[] : pointer to y.
 * \param dydx[] : pointer to right hand side of geodesic equation.
 */
bool FermiWalker::calcDerivsWL(const double y[], double dydx[])
{
    int mu, j, k, l, bidx;

    //  if (!mMetric->calcDerivsFW(mPropAcc,y,dydx))

    double c = mMetric->speed_of_light();

    mMetric->calculateMetric(y);
    mMetric->calculateChristoffels(y);

    for (mu = 0; mu < 4; mu++) {
        dydx[mu] = 0.0;
        dydx[mu + 4] = 0.0;
        dydx[mu + 8] = 0.0;

        for (j = 1; j < 4; j++) {
            bidx = mu + 4 * (j + 2);
            dydx[bidx] = 0.0;

            for (k = 0; k < 4; k++) {
                for (l = 0; l < 4; l++) {
                    dydx[bidx] -= mMetric->getChristoffel(k, l, mu) * y[k + 4] * y[8 + 4 * j + l];
                    dydx[bidx] -= mMetric->getMetricCoeff(k, l) * (y[k + 4] * mCurrAcc[mu] - mCurrAcc[k] * y[mu + 4])
                        * y[8 + 4 * j + l] / c;
                }
            }
        }
    }
    return true;
}

} // end namespace m4d
