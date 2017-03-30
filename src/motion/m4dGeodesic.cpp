// -------------------------------------------------------------------------------
/*
    m4dGeodesic.cpp

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

#include "m4dGeodesic.h"

namespace m4d {

/*! Standard constructor for geodesic motion.
 *
 *  \param metric : pointer to metric.
 *  \param type   : type of geodesic.
 *  \sa enum_geodesic_type.
 */
Geodesic::Geodesic(Metric* metric, enum_geodesic_type  type)
    : Motion(metric) {
    mType = type;
    setKappa();

    for (int i = 0; i < 4; i++) {
        mBoundBoxMin[i] = -DBL_MAX;
        mBoundBoxMax[i] = DBL_MAX;
    }
    //mBoundBoxMax[0] = 50.0;

    mCalcWithParTransport = false;
    mNumCoords = 8;

    mConstraintEpsilon = DEF_CONSTRAINT_EPSILON;
    resizeEps = DEF_RESIZE_EPSILON;
    resizeFac = DEF_RESIZE_FACTOR;
}


Geodesic::~Geodesic() {
}

// *********************************** public methods ******************************

/*! Set type of geodesic.
 *
 *  \param type   : type of geodesic.
 *  \sa enum_geodesic_type.
 */
void
Geodesic::setGeodesicType(enum_geodesic_type  type) {
    mType = type;
    setKappa();
}

/*! Get type of geodesic.
 *
 * \return type of geodesic.
 */
enum_geodesic_type
Geodesic::type() {
    return mType;
}


bool Geodesic::setParam(std::string paramName, bool val) {
    if (paramName.compare("stepctrl")==0) {
        mStepsizeControlled = val;
        return true;
    }
    return false;
}


bool Geodesic::setParam(std::string paramName, double value) {
    if (paramName.compare("stepctrl") == 0) {
        mStepsizeControlled = (int(value) == 1);
    }
    else if (paramName.compare("eps_a") == 0 || paramName.compare("epsilon_abs") == 0) {
        epsilon_abs = value;
        return true;
    }
    else if (paramName.compare("eps_r") == 0 || paramName.compare("epsilon_rel") == 0) {
        epsilon_rel = value;
        return true;
    }
    return false;
}


bool Geodesic::setParam(std::string paramName, double v0, double v1, double v2, double v3) {
    return Motion::setParam(paramName, v0, v1, v2, v3);
}


void Geodesic::setEpsilons(double eps_a, double eps_r) {
    epsilon_abs = eps_a;
    epsilon_rel = eps_r;
}


void Geodesic::getEpsilons(double &eps_a, double &eps_r) {
    eps_a = epsilon_abs;
    eps_r = epsilon_rel;
}

/*! Set/Unset stepsize control.
 *
 *  \param  control  :  if true then stepsize control is activated.
 */
void
Geodesic::setStepSizeControlled(bool control) {
    mStepsizeControlled = control;
}

/*! Set to 'true' if the local tetrad of the observer should be parallely transported along the null geodesic.
 *
 * \param calcwith : calculate with parallel transport.
 */
void
Geodesic::setCalcWithParTransport(bool calcwith) {
    mCalcWithParTransport = calcwith;
    if (mCalcWithParTransport) {
        mNumCoords = 24;
    } else {
        mNumCoords = 8;
    }
}

/*!
 * \return calculate with parallel transport.
 */
bool
Geodesic::calcWithParTransport() {
    return mCalcWithParTransport;
}

/*!
 * \param eps: epsilon for resize.
 * \param factor: factor for resize.
 */
void
Geodesic::setResize(double eps, double factor) {
    resizeEps = eps;
    resizeFac = factor;
}

/*!
 * \return epsilong for resize.
 */
void
Geodesic::getResize(double &eps, double &fac) {
    eps = resizeEps;
    fac = resizeFac;
}

/*! Initialize geodesic.
 * \param initPos : initial position.
 * \param initDir : initial coordinate direction.
 * \param cstr
 */
enum_break_condition
Geodesic::initializeGeodesic(const vec4 initPos, const vec4 initDir, double &cstr) {
    resetAffineParam();
    resetAffineParamStep();

    setInitialPosition(initPos);
    setInitialDirection(initDir);

    if ((cstr = fabs(testConstraint())) > mConstraintEpsilon) {
        return enum_break_constraint;
    }
    return enum_break_none;
}

/*! Test, whether the constraint condition is fulfilled:  approx=0 (yes)
 *
 * \return constraint value.
 */
double Geodesic::testConstraint() {
    return mMetric->testConstraint(y, mKappa);
}

/*! Print geodesic solver properties.
 * \param fptr : file pointer.
 */
void Geodesic::printF(FILE* fptr) {
    fprintf(fptr, "\nGeodesic:\n");
    fprintf(fptr, "\tstepsize controlled : %s\n", ((mStepsizeControlled) ? "yes" : "no"));
    fprintf(fptr, "\tepsilon_abs         : %12.8e\n", epsilon_abs);
    fprintf(fptr, "\tepsilon_rel         : %12.8e\n", epsilon_rel);
}

// ********************************* protected methods *****************************
void Geodesic::setKappa() {
    switch (mType) {
    default:
    case enum_geodesic_lightlike: {
        mKappa = 0.0;
        break;
    }
    case enum_geodesic_lightlike_sachs: {
        mKappa = 0.0;
        break;
    }
    case enum_geodesic_timelike: {
        mKappa = -1.0;
        break;
    }
    case enum_geodesic_spacelike: {
        mKappa = 1.0;
        break;
    }
    }
}

/*! Calculate the right side of the geodesic equation.
 *
 * \param y[] : pointer to y.
 * \param dydx[] : pointer to right hand side of geodesic equation.
 */
bool
Geodesic::calcDerivs(const double y[], double dydx[]) {
    //register int mu,j,k,l;
    register int mu, k, l;

    if (!mMetric->calcDerivs(y, dydx)) {
        //double ch;
        mMetric->calculateChristoffels(y);

        for (mu = 0; mu < 4; mu++) {
            dydx[mu]   = y[4 + mu];
            dydx[mu + 4] = 0.0;

            for (k = 0; k < 4; k++)
                for (l = 0; l < 4; l++) {
                    dydx[mu + 4] -= mMetric->getChristoffel(k, l, mu) * y[4 + k] * y[4 + l];
                }
        }
    }
    return true;
}

/*! Calculate the right side of the geodesic equation with parallel transport.
 *
 * \param y[] : pointer to y.
 * \param dydx[] : pointer to right hand side of parallel transport equation.
 */
bool
Geodesic::calcDerivsPar(const double y[], double dydx[]) {
    register int mu, j, k, l, bidx;

    if (!mMetric->calcDerivsPar(y, dydx)) {
        mMetric->calculateChristoffels(y);

        for (mu = 0; mu < 4; mu++) {
            dydx[mu]   = y[mu + 4];
            dydx[mu + 4] = 0.0;

            for (k = 0; k < 4; k++)
                for (l = 0; l < 4; l++) {
                    dydx[mu + 4] -= mMetric->getChristoffel(k, l, mu) * y[4 + k] * y[4 + l];
                }

            for (j = 0; j < 4; j++) {
                bidx = mu + 4 * (j + 2);
                dydx[bidx] = 0.0;

                for (k = 0; k < 4; k++)
                    for (l = 0; l < 4; l++) {
                        dydx[bidx] -= mMetric->getChristoffel(k, l, mu) * y[4 + k] * y[8 + 4 * j + l];
                    }
            }
        }
    }
    return true;
}

/*! Calculate the right side of the geodesic equation with parallel transport and Jacobi equation
 *
 * \param y[] : pointer to y.
 * \param dydx[] : pointer to right hand side of parallel transport equation.
 */
bool
Geodesic::calcDerivsSachsJacobi(const double y[], double dydx[]) {
    if (!mMetric->calcDerivsSachsJacobi(y, dydx)) {
        mMetric->calculateChristoffels(y);
        if (!mMetric->calculateChrisD(y)) {
            return false;
        }

        vec4 e[4];
        vec4 pos(y[0], y[1], y[2], y[3]);
        mMetric->getNatTetrad(pos, e[0], e[1], e[2], e[3]);

        for (int mu = 0; mu < 4; mu++) {
            dydx[mu]              = y[DEF_TG_IDX + mu];
            dydx[DEF_TG_IDX + mu]   = 0.0;

            dydx[DEF_JAC1_IDX + mu] = y[DEF_DJ1_IDX + mu];
            dydx[DEF_DJ1_IDX + mu]  = 0.0;
            dydx[DEF_JAC2_IDX + mu] = y[DEF_DJ2_IDX + mu];
            dydx[DEF_DJ2_IDX + mu]  = 0.0;

            dydx[DEF_SA1_IDX + mu]  = 0.0;
            dydx[DEF_SA2_IDX + mu]  = 0.0;

            for (int k = 0; k < 4; k++)
                for (int l = 0; l < 4; l++) {
                    dydx[DEF_TG_IDX + mu]  -= mMetric->getChristoffel(k, l, mu) * y[DEF_TG_IDX + k] * y[DEF_TG_IDX + l];
                    dydx[DEF_DJ1_IDX + mu] -= 2.0 * mMetric->getChristoffel(k, l, mu) * y[DEF_TG_IDX + k] * y[DEF_DJ1_IDX + l];
                    dydx[DEF_DJ2_IDX + mu] -= 2.0 * mMetric->getChristoffel(k, l, mu) * y[DEF_TG_IDX + k] * y[DEF_DJ2_IDX + l];

                    dydx[DEF_SA1_IDX + mu] -= mMetric->getChristoffel(k, l, mu) * y[DEF_TG_IDX + k] * y[DEF_SA1_IDX + l];
                    dydx[DEF_SA2_IDX + mu] -= mMetric->getChristoffel(k, l, mu) * y[DEF_TG_IDX + k] * y[DEF_SA2_IDX + l];

                    for (int n = 0; n < 4; n++) {
                        dydx[DEF_DJ1_IDX + mu] -= mMetric->getChrisD(k, l, mu, n) * y[DEF_TG_IDX + k] * y[DEF_TG_IDX + l] * y[DEF_JAC1_IDX + n];
                        dydx[DEF_DJ2_IDX + mu] -= mMetric->getChrisD(k, l, mu, n) * y[DEF_TG_IDX + k] * y[DEF_TG_IDX + l] * y[DEF_JAC2_IDX + n];
                    }
                }
        }
    }
    return true;
}

/*! Calculate Sachs basis from the initial local direction.
 * \param localNullDir: initial local direction of the light ray.
 * \param locX : local x-direction.
 * \param locY : local y-direction.
 * \param locZ : local z-direction.
 */
void
Geodesic::calcSachsBasis(const vec3 localNullDir, const vec3 locX, const vec3 locY, const vec3 locZ) {
    vec3 b1, b2;

    if ((locZ ^ localNullDir).isZero()) {
        b1 = locY;
        b2 = locX;
    } else {
        b1 = (locZ ^ localNullDir).getNormalized();
        b2 = (b1 ^ localNullDir).getNormalized();
    }
    //k.print(std::cerr,"k:  ");
    //b1.print(std::cerr,"b1: ");
    //b2.print(std::cerr,"b2: ");
    mSachsBasisB1 = vec3(b1[0], b1[1], b1[2]);
    mSachsBasisB2 = vec3(b2[0], b2[1], b2[2]);
}

/*! Set Sachs basis.
 * \param s1 : first sachs vector.
 * \param s2 : second sachs vector.
 */
void
Geodesic::setSachsBasis(const vec4 s1, const vec4 s2) {
    for (int mu = 0; mu < 4; mu++) {
        y[DEF_SA1_IDX + mu] = s1[mu];
        y[DEF_SA2_IDX + mu] = s2[mu];
    }
}

/*! Calculate the Jacobi parameters.
 * \param lambda : affine parameter.
 * \param y[]    : array
 * \param currJacobi : reference to jacobi values.
 */
void
Geodesic::calcJacobiParams(const double lambda, const double y[], vec5 &currJacobi) {
    vec4 pos   = vec4(&y[0]);
    vec4 jdir1 = vec4(&y[DEF_JAC1_IDX]);
    vec4 jdir2 = vec4(&y[DEF_JAC2_IDX]);
    vec4 cb1   = vec4(&y[DEF_SA1_IDX]);
    vec4 cb2   = vec4(&y[DEF_SA2_IDX]);

    // Calculate the Jacobian...
    double sj11, sj12, sj21, sj22;
    mMetric->calcProduct(pos, jdir1, cb1, sj11);
    mMetric->calcProduct(pos, jdir1, cb2, sj12, false);
    mMetric->calcProduct(pos, jdir2, cb1, sj21, false);
    mMetric->calcProduct(pos, jdir2, cb2, sj22, false);
    // fprintf(stdout,"%f %f %f %f   %f %f %f %f\n",y[DEF_JAC1_IDX+0],y[DEF_JAC1_IDX+1],y[DEF_JAC1_IDX+2],y[DEF_JAC1_IDX+3],sj11,sj12,sj21,sj22);

    double j11, j12, j21, j22;
    j11 = -sj11;
    j12 = -sj21;
    j21 = -sj12;
    j22 = -sj22;

    // Determine the shape parameters of a circle that is transformed by the Jacobian...
    double alpha = j11 * j12 + j21 * j22;
    double beta  = j11 * j11 - j12 * j12 + j21 * j21 - j22 * j22;
    double R     = j11 * j11 + j21 * j21;

    double zeta = 0.5 * atan(2.0 * alpha / beta);
    double zeta_plus = zeta + M_PI_2;

    double sz = sin(zeta);
    double cz = cos(zeta);

    double szp = sin(zeta_plus);
    double czp = cos(zeta_plus);

    double dm = sqrt(2.0 * alpha * sz * cz   - beta * sz * sz   + R);
    double dp = sqrt(2.0 * alpha * szp * czp - beta * szp * szp + R);

    double z1p, z2p, z1m, z2m;

    z1p = j11 * czp + j12 * szp;
    z2p = j21 * czp + j22 * szp;
    z1m = j11 * cz + j12 * sz;
    z2m = j21 * cz + j22 * sz;

#if 0
    fprintf(stderr, "%12.8f %12.8f %12.8f %12.8f\n", z1p, z2p, z1m, z2m);
#endif

  //  if (fabs(dm) > fabs(dp)) {
  //      std::swap(dm, dp);
  //      std::swap(z1p, z1m);
  //      std::swap(z2p, z2m);
  //  }

    double angle = atan2(z2p, z1p);
    double elipt = fabs(dp - dm) / fabs(dp + dm);
    double mu = 1.0;
    if (lambda > 0.0) {
        mu = lambda * lambda / (dp * dm);
    }
    currJacobi = vec5(dp, dm, mu, angle, elipt);
}

/*!
 */
void
Geodesic::findMaxJacobi(vec5 &currJacobi, vec5 &maxJacobi) {
    for (unsigned int i = 0; i < 5; i++)
        if (currJacobi[i] > maxJacobi[i]) {
            maxJacobi[i] = currJacobi[i];
        }
}

} // end namespace m4d
