// -------------------------------------------------------------------------------
/*
    m4dMotionChargedParticle.cpp

  Copyright (c) 2015 Thomas Mueller


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

#include "m4dMotionChargedParticle.h"

namespace m4d {

MotionChargedParticle::MotionChargedParticle(Metric* metric)
    : Motion(metric) {

    mLambda = 0.0;
}


MotionChargedParticle::~MotionChargedParticle() {
}

// *********************************** public methods ******************************

/*! Set initial velocity with respect to natural local tetrad.
 *
 *  \param fm : time direction (>0 : future, <0 : past).
 *  \param v  : absolute value of velocity with respect to natural local tetrad.
 *  \param theta :  angle
 *  \param phi   :  angle
 *  \param type  : type of natural local tetrad.
 *  \sa enum_nat_tetrad_type.
 */
bool MotionChargedParticle::setInitialVelocity(double fm, double v, double theta, double phi, double q_over_m, enum_nat_tetrad_type  type) {
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

    mQoverM = q_over_m;
    return true;
}


void MotionChargedParticle::getInitialVelocity(double &v, double &theta, double &phi) {
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
enum_break_condition
MotionChargedParticle::calculateMotion(const vec4 initPos, double fm, double v, double theta_v, double phi_v,
                               const double q_over_m,
                               const vec4 e0, const vec4 e1, const vec4 e2, const vec4 e3,
                               const int maxNumPoints,
                               std::vector<vec4> &points,
                               std::vector<vec4> &base0, std::vector<vec4> &base1, std::vector<vec4> &base2, std::vector<vec4> &base3) {
    if (fm == 0.0) {
        return enum_break_constraint;
    }
    fm = fm / fabs(fm); //normalize to 1.0

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
    setInitialVelocity(fm, v, theta_v, phi_v, q_over_m);

    points.push_back(initPos);
    base0.push_back(e0);
    base1.push_back(e1);
    base2.push_back(e2);
    base3.push_back(e3);


    enum_break_condition  breakType = enum_break_none;
    register int i = 0;

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



/*! Calculate the next step of the Fermi-Walker transport.
 *
 *  \return enum_break_condition.
 */
enum_break_condition
MotionChargedParticle::nextStep() {
    if (mMetric->breakCondition(&y[0])) {
        return enum_break_cond;
    }
    register int i;
    double yn[DEF_MAX_YS], dydx[DEF_MAX_YS], k1[DEF_MAX_YS], k2[DEF_MAX_YS], k3[DEF_MAX_YS], k4[DEF_MAX_YS];

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



// ********************************* protected methods *****************************
/*! Calculate right hand side of
 *
 * \param y[] : pointer to y.
 * \param dydx[] : pointer to right hand side of geodesic equation.
 */
bool
MotionChargedParticle::calcDerivs(const double y[], double dydx[]) {
    register int mu, k, l;

    mMetric->calculateChristoffels(y);
    mMetric->calcFmu_nu(y);

    for (mu = 0; mu < 4; mu++) {
        dydx[mu]     = y[4 + mu];
        dydx[mu + 4] = 0.0;

        for (k = 0; k < 4; k++) {
            for (l = 0; l < 4; l++) {
                dydx[mu + 4] -= mMetric->getChristoffel(k, l, mu) * y[4 + k] * y[4 + l];
            }
            dydx[mu + 4] -= mQoverM * mMetric->getFmu_nu(mu, k) * y[4 + k];
        }
    }

    return true;
}


} // end namespace m4d
