// -------------------------------------------------------------------------------
/*
    m4dGeodesicDP65.cpp

  Copyright (c) 2011  Thomas Mueller


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

#include "m4dGeodesicDP65.h"
#include <algorithm>

#ifdef USE_DP_INT

namespace m4d {

/*! Standard constructor for geodesic motion.
 *
 *  \param  metric : Metric of the spacetime where the geodesic has to be calculated.
 *  \param  type   : Type of geodesic:  kappa=-1 (timelike), kappa=0 (lightlike).
 */
GeodesicDP65::GeodesicDP65(Metric* metric, enum_geodesic_type  type)
    : Geodesic(metric, type) {
    mCalcWithParTransport = false;
    mNumCoords = 8;

    mConstraintEpsilon = DEF_CONSTRAINT_EPSILON;

    mStepsizeControlled = false;
    epsilon_abs = 1.0e-6;
    epsilon_rel = 0.0;

    mMaxLambdaStep = 1.0;

    mOrder         = 5;
    stepFac        = 1.0 / (1.0 - pow(2.0, -(double)mOrder));
    stepSigma      = 0.5;

    a21 = 1.0 / 10.0;

    a31 = -2.0 / 81.0;
    a32 = 20.0 / 81.0;
    //fprintf(stderr,"%12.10e\n",a31+a32);  // = 2/9

    a41 = 615.0 / 1372.0;
    a42 = -270.0 / 343.0;
    a43 = 1053.0 / 1372.0;
    //fprintf(stderr,"%12.10e\n",a41+a42+a43);  // = 3/7

    a51 = 3243.0 / 5500.0;
    a52 = -54.0 / 55.0;
    a53 = 50949.0 / 71500.0;
    a54 = 4998.0 / 17875.0;
    //fprintf(stderr,"%12.10e\n",a51+a52+a53+a54);  // = 3/5

    a61 = -26492.0 / 37125.0;
    a62 = 72.0 / 55.0;
    a63 = 2808.0 / 23375.0;
    a64 = -24206.0 / 37125.0;
    a65 = 338.0 / 459.0;
    //fprintf(stderr,"%12.10e\n",a61+a62+a63+a64+a65);  // = 4/5

    a71 = 5561.0 / 2376.0;
    a72 = -35.0 / 11.0;
    a73 = -24117.0 / 31603.0;
    a74 = 899983.0 / 200772.0;
    a75 = -5225. / 1836.0;
    a76 = 3925.0 / 4056.0;
    //fprintf(stderr,"%12.10e\n",a71+a72+a73+a74+a75+a76);  // = 1

    a81 = 465467.0 / 266112.0;
    a82 = -2945.0 / 1232.0;
    a83 = -5610201.0 / 14158144.0;
    a84 = 10513573.0 / 3212352.0;
    a85 = -424325.0 / 205632.0;
    a86 = 376225.0 / 454272.0;
    a87 = 0.0;
    //fprintf(stderr,"%12.10e\n",a81+a82+a83+a84+a85+a86+a87);  // = 1

    b1  = 821.0 / 10800.0;
    b2  = 0.0;
    b3  = 19683.0 / 71825.0;
    b4  = 175273.0 / 912600.0;
    b5  = 395.0 / 3672.0;
    b6  = 785.0 / 2704.0;
    b7  = 3.0 / 50.0;
    b8  = 0.0;
    //fprintf(stderr,"%12.10e\n",b1+b2+b3+b4+b5+b6+b7+b8);  // = 1

    db1 = 61.0 / 864.0;
    db2 = 0.0;
    db3 = 98415.0 / 321776.0;
    db4 = 16807.0 / 146016.0;
    db5 = 1375.0 / 7344.0;
    db6 = 1375.0 / 5408.0;
    db7 = -37.0 / 1120.0;
    db8 = 1.0 / 10.0;
    //fprintf(stderr,"%12.10e\n",db1+db2+db3+db4+db5+db6+db7+db8);  // = 1
}


GeodesicDP65::~GeodesicDP65() {
}

// *********************************** public methods ******************************

/*! Calculate a geodesic.
 *
 *  \param  initPos      :  initial position of the geodesic in coordinates.
 *  \param  initDir      :  initial direction of the geodesic in coordinates.
 *  \param  maxNumPoints :  maximum number of points.
 *  \param  points       :  reference to the calculated points.
 *  \param  dirs         :  reference to the calculated tangents.
 *  \param  lambda       :  reference to the affine parameters.
 *
 *  \return enum_break_condition : break condition.
 *  \sa enum_break_condition
 */
enum_break_condition
GeodesicDP65::calculateGeodesic(const vec4 initPos, const vec4 initDir, const int maxNumPoints,
                                  std::vector<vec4> &points, std::vector<vec4> &dirs, std::vector<double> &lambda) {
    if (!points.empty()) {
        points.clear();
    }
    if (!dirs.empty()) {
        dirs.clear();
    }
    if (!lambda.empty()) {
        lambda.clear();
    }

    enum_break_condition  breakType = enum_break_none;
    double cstr;
    breakType = initializeGeodesic(initPos, initDir, cstr);
    if (breakType == enum_break_constraint) {
        return enum_break_constraint;
    }

    points.push_back(vec4(y[0], y[1], y[2], y[3]));
    dirs.push_back(vec4(y[4], y[5], y[6], y[7]));
    lambda.push_back(mLambda);

    int status;
    double delta;
    int64_t t1 = get_system_clock();
    while ((int)points.size() < maxNumPoints && breakType == enum_break_none) {
        if (mMetric->breakCondition(&y[0])) {
            return enum_break_cond;
        }

        h = mLambdaStep;
        if (!mStepsizeControlled) {
            nextStep(status);
        } else {
            int maxAttempts = 3;
            for (int n = 0; n < maxAttempts; n++) {
                nextStep(status);
                delta = 0.0;
                for (int i = 0; i < 8; i++) {
                    delta = std::max(delta, fabs(yerr[i]));
                }
                delta *= stepFac / h;

                h = std::min(mMaxLambdaStep, stepSigma * h * pow(epsilon_abs / delta, 1.0 / (double)mOrder));

                if (delta < epsilon_abs) {
                    break;
                }
            }

            for (int i = 0; i < 8; i++) {
                y[i] = yn[i];
            }
            mLambdaStep = h;
        }
        mLambda += mLambdaStep;

        if (!outsideBoundBox()) {
            points.push_back(vec4(y[0], y[1], y[2], y[3]));
            dirs.push_back(vec4(y[4], y[5], y[6], y[7]));
            lambda.push_back(mLambda);

            if ((cstr = fabs(mMetric->testConstraint(y, mKappa))) > mConstraintEpsilon) {
                breakType = enum_break_constraint;
            }

            if (cstr > resizeEps) {
                mMetric->resize(y, mKappa, resizeFac);
            }
        } else {
            breakType = enum_break_outside;
        }
    }
    if ((int)points.size() >= maxNumPoints) {
        breakType = enum_break_num_exceed;
    }

    int64_t t2 = get_system_clock();
    mCalcTime = (t2 - t1) * 1e-6;
    return breakType;
}

enum_break_condition
GeodesicDP65::calculateGeodesic(const vec4 initPos, const vec4 initDir, const int maxNumPoints,
                                vec4 *&points, vec4 *&dirs, int &numPoints) {
    // TODO
    enum_break_condition  breakType = enum_break_none;
    return breakType;
}

/*! Calculate a geodesic and store points, directions, and constraint.
 *
 *  \param  initPos      :  initial position of the geodesic in coordinates.
 *  \param  initDir      :  initial direction of the geodesic in coordinates.
 *  \param  maxNumPoints :  maximum number of points.
 *  \param  points       :  reference to the calculated points.
 *  \param  dirs         :  reference to the calculated tangents.
 *  \param  epsilons     :  reference to the epsilons.
 *  \param  lambda       :  reference to the affine parameters.
 *
 *  \return enum_break_condition : break condition.
 *  \sa enum_break_condition
 */
enum_break_condition
GeodesicDP65::calculateGeodesicData(const vec4 initPos, const vec4 initDir, const int maxNumPoints,
                                      std::vector<vec4> &points, std::vector<vec4> &dirs,
                                      std::vector<double> &epsilons, std::vector<double> &lambda) {
    if (!points.empty()) {
        points.clear();
    }
    if (!dirs.empty()) {
        dirs.clear();
    }
    if (!epsilons.empty()) {
        epsilons.clear();
    }
    if (!lambda.empty()) {
        lambda.clear();
    }

    enum_break_condition  breakType = enum_break_none;
    double cstr;
    breakType = initializeGeodesic(initPos, initDir, cstr);
    if (breakType == enum_break_constraint) {
        return enum_break_constraint;
    }

    vec4 p(y, 4);
    vec4 d(&(y[4]), 4);
    points.push_back(p);
    dirs.push_back(d);
    epsilons.push_back(cstr);

    int64_t t1 = get_system_clock();
    int status;
    double delta;
    while ((int)points.size() < maxNumPoints && breakType == enum_break_none) {
        if (mMetric->breakCondition(&y[0])) {
            return enum_break_cond;
        }

        h = mLambdaStep;
        if (!mStepsizeControlled) {
            nextStep(status);
        } else {
            int maxAttempts = 3;
            for (int n = 0; n < maxAttempts; n++) {
                nextStep(status);
                delta = 0.0;
                for (int i = 0; i < 8; i++) {
                    delta = std::max(delta, fabs(yerr[i]));
                }
                delta *= stepFac / h;

                h = std::min(mMaxLambdaStep, stepSigma * h * pow(epsilon_abs / delta, 1.0 / (double)mOrder));

                if (delta < epsilon_abs) {
                    break;
                }
            }

            for (int i = 0; i < 8; i++) {
                y[i] = yn[i];
            }
            mLambdaStep = h;
        }
        mLambda += mLambdaStep;

        if (!outsideBoundBox()) {
            vec4 p(y, 4);
            vec4 d(&(y[4]), 4);
            points.push_back(p);
            dirs.push_back(d);
            epsilons.push_back(cstr);
            lambda.push_back(mLambda);

            if ((cstr = fabs(mMetric->testConstraint(y, mKappa))) > mConstraintEpsilon) {
                breakType = enum_break_constraint;
            }

            if (cstr > resizeEps) {
                mMetric->resize(y, mKappa, resizeFac);
            }
        } else {
            breakType = enum_break_outside;
        }
    }
    if ((int)points.size() >= maxNumPoints) {
        breakType = enum_break_num_exceed;
    }

    int64_t t2 = get_system_clock();
    mCalcTime = (t2 - t1) * 1e-6;
    return breakType;
}


/*! Calculate a geodesic and the parallel transported local tetrad of the observer.
 *
 *  \param  initPos      :  initial position of the geodesic in coordinates.
 *  \param  initDir      :  initial direction of the geodesic in coordinates.
 *  \param  e0           :  e0 base vectors of local tetrad.
 *  \param  e1           :  e1 base vectors of local tetrad.
 *  \param  e2           :  e2 base vectors of local tetrad.
 *  \param  e3           :  e3 base vectors of local tetrad.
 *  \param  maxNumPoints :  maximum number of points.
 *  \param  points       :  reference to the calculated points.
 *  \param  lambda       :  reference to the affine parameters.
 *  \param  dirs         :  reference to the calculated directions.
 *  \param  base0        :  references to the 0-base vectors of local tetrad.
 *  \param  base1        :  references to the 1-base vectors of local tetrad.
 *  \param  base2        :  references to the 2-base vectors of local tetrad.
 *  \param  base3        :  references to the 3-base vectors of local tetrad.
 *
 *  \return enum_break_condition : break condition.
 *  \sa enum_break_condition
 */
enum_break_condition
GeodesicDP65::calcParTransport(const vec4 initPos, const vec4 initDir,
                                 const vec4 e0, const vec4 e1, const vec4 e2, const vec4 e3,
                                 const int maxNumPoints,
                                 std::vector<vec4> &points,  std::vector<vec4> &dirs,
                                 std::vector<double> &lambda,
                                 std::vector<vec4> &base0, std::vector<vec4> &base1,
                                 std::vector<vec4> &base2, std::vector<vec4> &base3) {
    setCalcWithParTransport(true);

    if (!points.empty()) {
        points.clear();
    }
    if (!dirs.empty()) {
        dirs.clear();
    }
    if (!lambda.empty()) {
        lambda.clear();
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

    resetAffineParam();
    resetAffineParamStep();

    setInitialPosition(initPos);
    setInitialDirection(initDir);
    setInitialTetrad(e0, e1, e2, e3);


//  initPos.print(std::cerr);
//  initDir.print(std::cerr);
//   e0.print(std::cerr);
//   e1.print(std::cerr);
//   e2.print(std::cerr);
//   e3.print(std::cerr);

    enum_break_condition  breakType = enum_break_none;

    double cstr;
    if ((cstr = fabs(testConstraint())) > mConstraintEpsilon) {
        return enum_break_constraint;
    }

    points.push_back(vec4(y[0], y[1], y[2], y[3]));
    dirs.push_back(vec4(y[4], y[5], y[6], y[7]));
    lambda.push_back(mLambda);
    base0.push_back(e0);
    base1.push_back(e1);
    base2.push_back(e2);
    base3.push_back(e3);

    int status;
    int64_t t1 = get_system_clock();
    double delta;
    while ((int)points.size() < maxNumPoints && breakType == enum_break_none) {
        if (mMetric->breakCondition(&y[0])) {
            return enum_break_cond;
        }

        h = mLambdaStep;
        if (!mStepsizeControlled) {
            nextStep(status);
        } else {

            int maxAttempts = 3;
            for (int n = 0; n < maxAttempts; n++) {
                nextStepPar(status);
                delta = 0.0;
                for (int i = 0; i < DEF_MAX_YS_PAR; i++) {
                    delta = std::max(delta, fabs(yerr[i]));
                }
                delta *= stepFac / h;

                h = std::min(mMaxLambdaStep, stepSigma * h * pow(epsilon_abs / delta, 1.0 / (double)mOrder));

                if (delta < epsilon_abs) {
                    break;
                }
            }

            for (int i = 0; i < DEF_MAX_YS_PAR; i++) {
                y[i] = yn[i];
            }
            mLambdaStep = h;
        }
        mLambda += mLambdaStep;

        if (!outsideBoundBox()) {
            points.push_back(vec4(y[0], y[1], y[2], y[3]));
            dirs.push_back(vec4(y[4], y[5], y[6], y[7]));
            base0.push_back(vec4(y[8], y[9], y[10], y[11]));
            base1.push_back(vec4(y[12], y[13], y[14], y[15]));
            base2.push_back(vec4(y[16], y[17], y[18], y[19]));
            base3.push_back(vec4(y[20], y[21], y[22], y[23]));
            lambda.push_back(mLambda);

            if ((cstr = fabs(mMetric->testConstraint(y, mKappa))) > mConstraintEpsilon) {
                breakType = enum_break_constraint;
            }

            if (cstr > resizeEps) {
                mMetric->resize(y, mKappa, resizeFac);
            }
        } else {
            breakType = enum_break_outside;
        }
    }
    if ((int)points.size() >= maxNumPoints) {
        breakType = enum_break_num_exceed;
    }

    int64_t t2 = get_system_clock();
    mCalcTime = (t2 - t1) * 1e-6;
    return breakType;
}

/*! Calculate a geodesic with Sachs-basis transport and Jacobi equation.
 *
 *  \param  initPos      :  initial position of the geodesic in coordinates.
 *  \param  initCoordDir :  initial direction of the geodesic in coordinates.
 *  \param  localNullDir :  local initial direction of the null geodesinc.
 *  \param  locX         :  locX.
 *  \param  locY         :  locY.
 *  \param  locZ         :  locZ.
 *  \param  b0           :  b0 base vectors of local tetrad.
 *  \param  b1           :  b1 base vectors of local tetrad.
 *  \param  b2           :  b2 base vectors of local tetrad.
 *  \param  b3           :  b3 base vectors of local tetrad.
 *  \param  tetrad_type  :  type of tetrad.
 *  \param  maxNumPoints :  maximum number of points.
 *  \param  points       :  reference to the calculated points.
 *  \param  dirs         :  reference to the calculated directions.
 *  \param  lambda       :  reference to the affine parameters.
 *  \param  sachs0       :  reference to the parallel transported sachs basis.
 *  \param  sachs1       :  reference to the parallel transported sachs basis.
 *  \param  jacobi       :  reference to the Jacobi parameters.
 *  \param  maxJacobi    :  reference to maximum Jacobi parameters.
 *
 *  \return enum_break_condition : break condition.
 *  \sa enum_break_condition
 */
enum_break_condition
GeodesicDP65::calcSachsJacobi(const vec4 initPos, const vec4 initCoordDir,
                                const vec3 localNullDir, const vec3 locX, const vec3 locY, const vec3 locZ,
                                const vec4 b0, const vec4 b1, const vec4 b2, const vec4 b3,
                                const enum_nat_tetrad_type  tetrad_type,
                                const int maxNumPoints,
                                std::vector<vec4> &points,  std::vector<vec4> &dirs,
                                std::vector<double> &lambda,
                                std::vector<vec4> &sachs0, std::vector<vec4> &sachs1,
                                std::vector<vec5> &jacobi, vec5 &maxJacobi) {
    setCalcWithParTransport(true);

    if (!points.empty()) {
        points.clear();
    }
    if (!dirs.empty()) {
        dirs.clear();
    }
    if (!lambda.empty()) {
        lambda.clear();
    }
    if (!sachs0.empty()) {
        sachs0.clear();
    }
    if (!sachs1.empty()) {
        sachs1.clear();
    }
    if (!jacobi.empty()) {
        jacobi.clear();
    }

    maxJacobi = vec5(-DBL_MAX, -DBL_MAX, -DBL_MAX, -DBL_MAX, -DBL_MAX);

    resetAffineParam();
    resetAffineParamStep();

    setInitialPosition(initPos);
    setInitialDirection(initCoordDir);

    // calculate Sachs basis -> mSachsBasisB1,mSachsBasisB2
    calcSachsBasis(localNullDir, locX, locY, locZ);

    vec4 e0, e1, e2, e3;
    mMetric->localToCoord(initPos, b0, e0, tetrad_type);
    mMetric->localToCoord(initPos, b1, e1, tetrad_type);
    mMetric->localToCoord(initPos, b2, e2, tetrad_type);
    mMetric->localToCoord(initPos, b3, e3, tetrad_type);

    vec4 bb1, bb2;
    bb1 = mSachsBasisB1[0] * e1 + mSachsBasisB1[1] * e2 + mSachsBasisB1[2] * e3;
    bb2 = mSachsBasisB2[0] * e1 + mSachsBasisB2[1] * e2 + mSachsBasisB2[2] * e3;
    setSachsBasis(bb1, bb2);

    for (int i = 0; i < 4; i++) {
        y[DEF_JAC1_IDX + i] = 0.0;
        y[DEF_DJ1_IDX + i]  = bb1[i];

        y[DEF_JAC2_IDX + i] = 0.0;
        y[DEF_DJ2_IDX + i]  = bb2[i];
    }

    enum_break_condition  breakType = enum_break_none;
    double cstr;
    if ((cstr = fabs(testConstraint())) > mConstraintEpsilon) {
        return enum_break_constraint;
    }

    vec5 currJacobi;

    points.push_back(vec4(&y[0]));
    dirs.push_back(vec4(&y[DEF_TG_IDX]));
    lambda.push_back(mLambda);
    sachs0.push_back(vec4(&y[DEF_SA1_IDX]));
    sachs1.push_back(vec4(&y[DEF_SA2_IDX]));
    jacobi.push_back(vec5(0.0, 0.0, 1.0, 0.0, 0.0));
    maxJacobi = jacobi[0];

    int status;
    int64_t t1 = get_system_clock();
    double delta;
    while ((int)points.size() < maxNumPoints && breakType == enum_break_none) {
        if (mMetric->breakCondition(&y[0])) {
            return enum_break_cond;
        }

        h = mLambdaStep;
        if (!mStepsizeControlled) {
            nextStepSachsJacobi(status);
        } else {
            int maxAttempts = 3;
            for (int n = 0; n < maxAttempts; n++) {
                nextStepSachsJacobi(status);
                delta = 0.0;
                for (int i = 0; i < DEF_MAX_YS_JAC; i++) {
                    delta = std::max(delta, fabs(yerr[i]));
                }

                delta *= stepFac / h;

                h = std::min(mMaxLambdaStep, stepSigma * h * pow(epsilon_abs / delta, 1.0 / (double)mOrder));

                if (h < 1e-6) {
                    h = 1e-6;
                }
                if (delta < epsilon_abs) {
                    break;
                }
            }

            for (int i = 0; i < DEF_MAX_YS_JAC; i++) {
                y[i] = yn[i];
            }
            mLambdaStep = h;
        }
        mLambda += mLambdaStep;

        if (!outsideBoundBox()) {
            points.push_back(vec4(&y[0]));
            dirs.push_back(vec4(&y[DEF_TG_IDX]));
            sachs0.push_back(vec4(&y[DEF_SA1_IDX]));
            sachs1.push_back(vec4(&y[DEF_SA2_IDX]));

            calcJacobiParams(mLambda, y, currJacobi);
            jacobi.push_back(currJacobi);
            lambda.push_back(mLambda);
            findMaxJacobi(currJacobi, maxJacobi);

            if ((cstr = fabs(mMetric->testConstraint(y, mKappa))) > mConstraintEpsilon) {
                breakType = enum_break_constraint;
            }

            if (cstr > resizeEps) {
                mMetric->resize(y, mKappa, resizeFac);
            }
        } else {
            breakType = enum_break_outside;
        }
    }
    if ((int)points.size() >= maxNumPoints) {
        breakType = enum_break_num_exceed;
    }

    int64_t t2 = get_system_clock();
    mCalcTime = (t2 - t1) * 1e-6;
    return breakType;
}

enum_break_condition
GeodesicDP65::calcSachsJacobi(const vec4 initPos, const vec4 initCoordDir,
                              const vec3 localNullDir, const vec3 locX, const vec3 locY, const vec3 locZ,
                              const vec4 b0, const vec4 b1, const vec4 b2, const vec4 b3,
                              const enum_nat_tetrad_type  tetrad_type,
                              const int maxNumPoints,
                              vec4 *&points, vec4 *&dirs,
                              double *&lambda,
                              vec4 *&sachs0, vec4 *&sachs1,
                              vec5 *&jacobi, vec5 &maxJacobi, int &numPoints) {
    // TODO
    enum_break_condition  breakType = enum_break_none;
    return breakType;
}

/*! Calculate next Runge-Kutta step of the geodesic.
 *
 *  \param yo : pointer to old values.
 *  \param yn : pointer to new values.
 *  \param yerr : pointer to error estimate.
 *  \param h : parameter step.
 */
bool
//GeodesicDP65::nextStep ( double* yo, double* yn, double* yerr, double h )
GeodesicDP65::nextStep(int &status) {
    register int i;
    double k1[8], k2[8], k3[8], k4[8], k5[8], k6[8], k7[8], k8[8], yy[8], dydx[8];

    double* yo = y;
    calcDerivs(yo, dydx);
    for (i = 0; i < 8; i++) {
        k1[i] = h * dydx[i];
        yy[i] = yo[i] + a21 * k1[i];
    }

    calcDerivs(yy, dydx);
    for (i = 0; i < 8; i++) {
        k2[i] = h * dydx[i];
        yy[i] = yo[i] + a31 * k1[i]  + a32 * k2[i];
    }

    calcDerivs(yy, dydx);
    for (i = 0; i < 8; i++) {
        k3[i] = h * dydx[i];
        yy[i] = yo[i] + a41 * k1[i] + a42 * k2[i] + a43 * k3[i];
    }

    calcDerivs(yy, dydx);
    for (i = 0; i < 8; i++) {
        k4[i] = h * dydx[i];
        yy[i] = yo[i] + a51 * k1[i] + a52 * k2[i] + a53 * k3[i] + a54 * k4[i];
    }

    calcDerivs(yy, dydx);
    for (i = 0; i < 8; i++) {
        k5[i] = h * dydx[i];
        yy[i] = yo[i] + a61 * k1[i] + a62 * k2[i] + a63 * k3[i] + a64 * k4[i] + a65 * k5[i];
    }

    calcDerivs(yy, dydx);
    for (i = 0; i < 8; i++) {
        k6[i] = h * dydx[i];
        yy[i] = yo[i] + a71 * k1[i] + a72 * k2[i] + a73 * k3[i] + a74 * k4[i] + a75 * k5[i] + a76 * k6[i];
    }

    calcDerivs(yy, dydx);
    for (int i = 0; i < 8; i++) {
        k7[i] = h * dydx[i];
        yy[i] = yo[i] + a81 * k1[i] + a82 * k2[i] + a83 * k3[i] + a84 * k4[i] + a85 * k5[i] + a86 * k6[i] + a87 * k7[i];
    }

    calcDerivs(yy, dydx);
    for (int i = 0; i < 8; i++) {
        k8[i] = h * dydx[i];
        yn[i] = yo[i] + b1 * k1[i] + b2 * k2[i] + b3 * k3[i] + b4 * k4[i] + b5 * k5[i] + b6 * k6[i] + b7 * k7[i] + b8 * k8[i];
        yerr[i] = yn[i] - (yo[i] + db1 * k1[i] + db2 * k2[i] + db3 * k3[i] + db4 * k4[i] + db5 * k5[i] + db6 * k6[i] + db7 * k7[i] + db8 * k8[i]);
    }
    status = 1;
    return true;
}

/*! Calculate next Runge-Kutta step of the geodesic with parallel transport.
 *
 *  \param yo : pointer to old values.
 *  \param yn : pointer to new values.
 *  \param yerr : pointer to error estimate.
 *  \param h : parameter step.
 */
bool
//GeodesicDP65::nextStepPar ( double* yo, double* yn, double* yerr, double h )
GeodesicDP65::nextStepPar(int &status) {
    register int i;
    double k1[DEF_MAX_YS_PAR], k2[DEF_MAX_YS_PAR], k3[DEF_MAX_YS_PAR], k4[DEF_MAX_YS_PAR],
           k5[DEF_MAX_YS_PAR], k6[DEF_MAX_YS_PAR], k7[DEF_MAX_YS_PAR], yy[DEF_MAX_YS_PAR],
           dydx[DEF_MAX_YS_PAR];

    double* yo = y;
    calcDerivsPar(yo, dydx);
    for (i = 0; i < DEF_MAX_YS_PAR; i++) {
        k1[i] = h * dydx[i];
        yy[i] = yo[i] + a21 * k1[i];
    }

    calcDerivsPar(yy, dydx);
    for (i = 0; i < DEF_MAX_YS_PAR; i++) {
        k2[i] = h * dydx[i];
        yy[i] = yo[i] + a31 * k1[i]  + a32 * k2[i];
    }

    calcDerivsPar(yy, dydx);
    for (i = 0; i < DEF_MAX_YS_PAR; i++) {
        k3[i] = h * dydx[i];
        yy[i] = yo[i] + a41 * k1[i] + a42 * k2[i] + a43 * k3[i];
    }

    calcDerivsPar(yy, dydx);
    for (i = 0; i < DEF_MAX_YS_PAR; i++) {
        k4[i] = h * dydx[i];
        yy[i] = yo[i] + a51 * k1[i] + a52 * k2[i] + a53 * k3[i] + a54 * k4[i];
    }

    calcDerivsPar(yy, dydx);
    for (i = 0; i < DEF_MAX_YS_PAR; i++) {
        k5[i] = h * dydx[i];
        yy[i] = yo[i] + a61 * k1[i] + a62 * k2[i] + a63 * k3[i] + a64 * k4[i] + a65 * k5[i];
    }

    calcDerivsPar(yy, dydx);
    for (i = 0; i < DEF_MAX_YS_PAR; i++) {
        k6[i] = h * dydx[i];
        yy[i] = yo[i] + a71 * k1[i] + a72 * k2[i] + a73 * k3[i] + a74 * k4[i] + a75 * k5[i] + a76 * k6[i];
    }

    calcDerivsPar(yy, dydx);
    for (i = 0; i < DEF_MAX_YS_PAR; i++) {
        k7[i] = h * dydx[i];
        yn[i] = yo[i] + b1 * k1[i] + b2 * k2[i] + b3 * k3[i] + b4 * k4[i] + b5 * k5[i] + b6 * k6[i] + b7 * k7[i];
        yerr[i] = yn[i] - (yo[i] + db1 * k1[i] + db2 * k2[i] + db3 * k3[i] + db4 * k4[i] + db5 * k5[i] + db6 * k6[i] + db7 * k7[i]);
    }

    status = 1;
    return true;
}

/*! Calculate next Runge-Kutta step of the geodesic with parallel transport and Jacobi equation.
 *
 *  \param yo : pointer to old values.
 *  \param yn : pointer to new values.
 *  \param yerr : pointer to error estimate.
 *  \param h : parameter step.
 */
bool
//GeodesicDP65::nextStepSachsJacobi ( double* yo, double* yn, double* yerr, double h )
GeodesicDP65::nextStepSachsJacobi(int &status) {
    register int i;
    double k1[DEF_MAX_YS_JAC], k2[DEF_MAX_YS_JAC], k3[DEF_MAX_YS_JAC], k4[DEF_MAX_YS_JAC],
           k5[DEF_MAX_YS_JAC], k6[DEF_MAX_YS_JAC], k7[DEF_MAX_YS_JAC], k8[DEF_MAX_YS_JAC],
           yy[DEF_MAX_YS_JAC], dydx[DEF_MAX_YS_JAC];

    double* yo = y;
    calcDerivsSachsJacobi(yo, dydx);
    for (i = 0; i < DEF_MAX_YS_JAC; i++) {
        k1[i] = h * dydx[i];
        yy[i] = yo[i] + a21 * k1[i];
    }

    calcDerivsSachsJacobi(yy, dydx);
    for (i = 0; i < DEF_MAX_YS_JAC; i++) {
        k2[i] = h * dydx[i];
        yy[i] = yo[i] + a31 * k1[i]  + a32 * k2[i];
    }

    calcDerivsSachsJacobi(yy, dydx);
    for (i = 0; i < DEF_MAX_YS_JAC; i++) {
        k3[i] = h * dydx[i];
        yy[i] = yo[i] + a41 * k1[i] + a42 * k2[i] + a43 * k3[i];
    }

    calcDerivsSachsJacobi(yy, dydx);
    for (i = 0; i < DEF_MAX_YS_JAC; i++) {
        k4[i] = h * dydx[i];
        yy[i] = yo[i] + a51 * k1[i] + a52 * k2[i] + a53 * k3[i] + a54 * k4[i];
    }

    calcDerivsSachsJacobi(yy, dydx);
    for (i = 0; i < DEF_MAX_YS_JAC; i++) {
        k5[i] = h * dydx[i];
        yy[i] = yo[i] + a61 * k1[i] + a62 * k2[i] + a63 * k3[i] + a64 * k4[i] + a65 * k5[i];
    }

    calcDerivsSachsJacobi(yy, dydx);
    for (i = 0; i < DEF_MAX_YS_JAC; i++) {
        k6[i] = h * dydx[i];
        yy[i] = yo[i] + a71 * k1[i] + a72 * k2[i] + a73 * k3[i] + a74 * k4[i] + a75 * k5[i] + a76 * k6[i];
    }

    calcDerivsSachsJacobi(yy, dydx);
    for (int i = 0; i < DEF_MAX_YS_JAC; i++) {
        k7[i] = h * dydx[i];
        yy[i] = yo[i] + a81 * k1[i] + a82 * k2[i] + a83 * k3[i] + a84 * k4[i] + a85 * k5[i] + a86 * k6[i] + a87 * k7[i];
    }

    calcDerivsSachsJacobi(yy, dydx);
    for (int i = 0; i < DEF_MAX_YS_JAC; i++) {
        k8[i] = h * dydx[i];
        yn[i] = yo[i] + b1 * k1[i] + b2 * k2[i] + b3 * k3[i] + b4 * k4[i] + b5 * k5[i] + b6 * k6[i] + b7 * k7[i] + b8 * k8[i];
        yerr[i] = yn[i] - (yo[i] + db1 * k1[i] + db2 * k2[i] + db3 * k3[i] + db4 * k4[i] + db5 * k5[i] + db6 * k6[i] + db7 * k7[i] + db8 * k8[i]);
    }

    status = 1;
    return true;
}

/*! Print geodesic solver properties.
 * \param fptr : file pointer.
 */
void
GeodesicDP65::print(FILE* fptr) {
    fprintf(fptr, "\nGeodesicDP - RK6(5):\n--------------------\n");
    fprintf(fptr, "\tstepsize controlled : %s\n", ((mStepsizeControlled) ? "yes" : "no"));
    fprintf(fptr, "\tstep size           : %12.8e\n", mLambdaStep);
    fprintf(fptr, "\tepsilon             : %12.8e\n", epsilon_abs);
    fprintf(fptr, "\tconstraint epsilon  : %12.8e\n", mConstraintEpsilon);
    fprintf(fptr, "\tbounding box min    : %14.6e %14.6e %14.6e %14.6e\n", mBoundBoxMin[0], mBoundBoxMin[1], mBoundBoxMin[2], mBoundBoxMin[3]);
    fprintf(fptr, "\tbounding box max    : %14.6e %14.6e %14.6e %14.6e\n", mBoundBoxMax[0], mBoundBoxMax[1], mBoundBoxMax[2], mBoundBoxMax[3]);
}

} // end namespace m4d

#endif // USE_DP_INT
