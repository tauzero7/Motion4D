// -------------------------------------------------------------------------------
/*
    m4dGeodesicBS.cpp

  Copyright (c) 2010  Thomas Mueller


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
   aith the m4d-library.  If not, see <http://www.gnu.org/licenses/>.

*/
// -------------------------------------------------------------------------------

#include "m4dGeodesicBS.h"
#include <algorithm>
#include <limits>

namespace m4d {

/*! Standard constructor for geodesic motion.
 *
 *  \param  metric : Metric of the spacetime where the geodesic has to be calculated.
 *  \param  type   : Type of geodesic:  kappa=-1 (timelike), kappa=0 (lightlike).
 */
GeodesicBS::GeodesicBS(Metric* metric, enum_geodesic_type type)
    : Geodesic(metric, type)
{
    mNumCoords = 8;

    mConstraintEpsilon = DEF_CONSTRAINT_EPSILON;
    epsilon_abs = 1e-6;

    dd = new double*[DEF_MAX_YS];
    for (int i = 0; i < DEF_MAX_YS; i++) {
        dd[i] = new double[DEF_BS_MAX_ROW_NUM];
    }

    mhmin = 1e-16;
    mhmax = mMaxLambdaStep = 1e-2;
    for (int k = 0; k < DEF_BS_MAX_ROW_NUMS; k++) {
        nSeq[k] = 2 * (k + 1);
    }
}

GeodesicBS::~GeodesicBS()
{
    for (int i = 0; i < DEF_MAX_YS; i++) {
        delete[] dd[i];
    }
    delete[] dd;
    dd = NULL;
}

// *********************************** public methods ******************************

void GeodesicBS::setNumberOfSteps(int numSteps)
{
    mNumSteps = numSteps;
}

void GeodesicBS::setMaxAffineParamStep(double hmax)
{
    mMaxLambdaStep = mhmax = hmax;
    std::cerr << "BS max step: " << mhmax << std::endl;
}

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
enum_break_condition GeodesicBS::calculateGeodesic(const vec4 initPos, const vec4 initDir, const int maxNumPoints,
    std::vector<vec4>& points, std::vector<vec4>& dirs, std::vector<double>& lambda)
{
    register int i;

    if (!points.empty()) {
        points.clear();
    }
    if (!dirs.empty()) {
        dirs.clear();
    }
    if (!lambda.empty()) {
        lambda.clear();
    }

    resetAffineParam();
    resetAffineParamStep();

    setInitialPosition(initPos);
    setInitialDirection(initDir);

    enum_break_condition breakType = enum_break_none;

    if (fabs(testConstraint()) > mConstraintEpsilon) {
        // cout << "error testConstraint() " << testConstraint() << endl;
        return enum_break_constraint;
    }
    points.push_back(vec4(y[0], y[1], y[2], y[3]));
    dirs.push_back(vec4(y[4], y[5], y[6], y[7]));
    lambda.push_back(mLambda);

    int numOK = 0, numBAD = 0; // count=0;

    double h = mLambdaStep;
    double hdid, hnext;
    double cstr;

    if ((cstr = fabs(testConstraint())) > mConstraintEpsilon) {
        return enum_break_constraint;
    }

    mNumCoords = 8;

    mFirst = 1;
    mExitFlag = 0;
    epsold = -1.0;
    eps = epsilon_abs;
    while ((int)points.size() < maxNumPoints && breakType == enum_break_none) {

        calcDerivs(y, dydx);

        for (i = 0; i < mNumCoords; i++) {
            yscal[i] = fabs(y[i]) + fabs(dydx[i] * h) + DEF_BS_TINY;
        }

        // overshoot ignored

        breakType = nextStep(h, hdid, hnext, cstr);

        if (!outsideBoundBox()) {
            points.push_back(vec4(y[0], y[1], y[2], y[3]));
            dirs.push_back(vec4(y[4], y[5], y[6], y[7]));
            lambda.push_back(mLambda);
        }
        else {
            breakType = enum_break_outside;
        }

        if (hdid == h) {
            ++numOK;
        }
        else {
            ++numBAD;
        }

        if (fabs(hnext) <= mhmin) {
            fprintf(stderr, "Stepsize too small in GeodesicBS: %e\n", hnext);
        }
        h = hnext;
        if (h > mhmax) {
            h = mhmax; // Maximale Schrittweite, um "runde" Geodäten zu bekommen
        }
    }
    if ((int)points.size() >= maxNumPoints) {
        breakType = enum_break_num_exceed;
    }

    return breakType;
}

enum_break_condition GeodesicBS::calculateGeodesic(
    const vec4 initPos, const vec4 initDir, const int maxNumPoints, vec4*& points, vec4*& dirs, int& numPoints)
{
    register int i;

    if (points != NULL) {
        delete[] points;
    }
    if (dirs != NULL) {
        delete[] dirs;
    }

    points = new vec4[maxNumPoints];
    dirs = new vec4[maxNumPoints];

    resetAffineParam();
    resetAffineParamStep();

    setInitialPosition(initPos);
    setInitialDirection(initDir);

    enum_break_condition breakType = enum_break_none;

    if (fabs(testConstraint()) > mConstraintEpsilon) {
        // cout << "error testConstraint() " << testConstraint() << endl;
        return enum_break_constraint;
    }

    points[0] = vec4(y[0], y[1], y[2], y[3]);
    dirs[0] = vec4(y[4], y[5], y[6], y[7]);

    int numOK = 0, numBAD = 0; // count=0;

    double h = mLambdaStep;
    double hdid, hnext;
    double cstr;

    if ((cstr = fabs(testConstraint())) > mConstraintEpsilon) {
        return enum_break_constraint;
    }

    mNumCoords = 8;

    mFirst = 1;
    mExitFlag = 0;
    epsold = -1.0;
    eps = epsilon_abs;

    numPoints = 1;
    while (numPoints < maxNumPoints && breakType == enum_break_none) {
        calcDerivs(y, dydx);

        for (i = 0; i < mNumCoords; i++) {
            yscal[i] = fabs(y[i]) + fabs(dydx[i] * h) + DEF_BS_TINY;
        }

        // overshoot ignored

        breakType = nextStep(h, hdid, hnext, cstr);

        if (!outsideBoundBox()) {
            points[numPoints] = vec4(y[0], y[1], y[2], y[3]);
            dirs[numPoints] = vec4(y[4], y[5], y[6], y[7]);
            numPoints++;
        }
        else {
            breakType = enum_break_outside;
        }

        if (hdid == h) {
            ++numOK;
        }
        else {
            ++numBAD;
        }

        if (fabs(hnext) <= mhmin) {
            fprintf(stderr, "Stepsize too small in GeodesicBS: %e\n", hnext);
        }
        h = hnext;
        if (h > mhmax) {
            h = mhmax; // Maximale Schrittweite, um "runde" Geodäten zu bekommen
        }
    }
    if (numPoints == maxNumPoints) {
        breakType = enum_break_num_exceed;
    }

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
enum_break_condition GeodesicBS::calculateGeodesicData(const vec4 initPos, const vec4 initDir, const int maxNumPoints,
    std::vector<vec4>& points, std::vector<vec4>& dirs, std::vector<double>& epsilons, std::vector<double>& lambda)
{
    register int i;

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

    resetAffineParam();
    resetAffineParamStep();

    setInitialPosition(initPos);
    setInitialDirection(initDir);

    enum_break_condition breakType = enum_break_none;

    if (fabs(testConstraint()) > mConstraintEpsilon) {
        // cout << "error testConstraint() " << testConstraint() << endl;
        return enum_break_constraint;
    }
    points.push_back(vec4(y[0], y[1], y[2], y[3]));
    lambda.push_back(mLambda);

    int numOK = 0, numBAD = 0; // count=0;

    double h = mLambdaStep;
    double hdid, hnext;
    double cstr;

    if ((cstr = fabs(testConstraint())) > mConstraintEpsilon) {
        return enum_break_constraint;
    }

    mFirst = 1;
    mExitFlag = 0;
    epsold = -1.0;
    eps = epsilon_abs;
    while ((int)points.size() < maxNumPoints && breakType == enum_break_none) {

        calcDerivs(y, dydx);

        for (i = 0; i < mNumCoords; i++) {
            yscal[i] = fabs(y[i]) + fabs(dydx[i] * h) + DEF_BS_TINY;
        }

        // overshoot ignored

        breakType = nextStep(h, hdid, hnext, cstr);

        if (!outsideBoundBox()) {
            vec4 p(y, 4);
            vec4 d(&(y[4]), 4);
            points.push_back(p);
            dirs.push_back(d);
            epsilons.push_back(cstr);
            lambda.push_back(mLambda);
        }
        else {
            breakType = enum_break_outside;
        }

        if (hdid == h) {
            ++numOK;
        }
        else {
            ++numBAD;
        }

        if (fabs(hnext) <= mhmin) {
            fprintf(stderr, "Stepsize too small in GeodesicBS.\n");
        }
        h = hnext;
        if (h > mhmax) {
            h = mhmax; // Maximale Schrittweite, um "runde" Geodten zu bekommen
        }
    }
    if ((int)points.size() >= maxNumPoints) {
        breakType = enum_break_num_exceed;
    }

    return breakType;
}

enum_break_condition GeodesicBS::calcParTransport(const vec4, const vec4, const vec4, const vec4, const vec4,
    const vec4, const int, std::vector<vec4>& points, std::vector<vec4>& dirs, std::vector<double>& lambda,
    std::vector<vec4>&, std::vector<vec4>&, std::vector<vec4>&, std::vector<vec4>&)
{
    if (!points.empty()) {
        points.clear();
    }
    if (!dirs.empty()) {
        dirs.clear();
    }
    if (!lambda.empty()) {
        lambda.clear();
    }

    printf("calcParTransport not implemtent in BS\n");
    return enum_break_not_implemented;
}

enum_break_condition GeodesicBS::calcSachsJacobi(const vec4 initPos, const vec4 initCoordDir, const vec3 localNullDir,
    const vec3 locX, const vec3 locY, const vec3 locZ, const vec4 b0, const vec4 b1, const vec4 b2, const vec4 b3,
    const enum_nat_tetrad_type tetrad_type, const int maxNumPoints, std::vector<vec4>& points, std::vector<vec4>& dirs,
    std::vector<double>& lambda, std::vector<vec4>& sachs0, std::vector<vec4>& sachs1, std::vector<vec5>& jacobi,
    vec5& maxJacobi)
{
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

    double DOUBLE_MAX = std::numeric_limits<double>::max();
    maxJacobi = vec5(-DOUBLE_MAX, -DOUBLE_MAX, -DOUBLE_MAX, -DOUBLE_MAX, -DOUBLE_MAX);

    resetAffineParam();
    resetAffineParamStep();

    setInitialPosition(initPos);
    setInitialDirection(initCoordDir);

    calcSachsBasis(localNullDir, locX, locY, locZ);

    vec4 e0, e1, e2, e3;
    mMetric->localToCoord(initPos, b0, e0, tetrad_type);
    mMetric->localToCoord(initPos, b1, e1, tetrad_type);
    mMetric->localToCoord(initPos, b2, e2, tetrad_type);
    mMetric->localToCoord(initPos, b3, e3, tetrad_type);

    vec4 bb1, bb2;
    bb1 = mSachsBasisB1[0] * e1 + mSachsBasisB1[1] * e2 + mSachsBasisB1[2] * e3;
    bb2 = mSachsBasisB2[0] * e1 + mSachsBasisB2[1] * e2 + mSachsBasisB2[2] * e3;

    for (int i = 0; i < 4; i++) {
        y[DEF_JAC1_IDX + i] = 0.0;
        y[DEF_DJ1_IDX + i] = bb1[i];

        y[DEF_JAC2_IDX + i] = 0.0;
        y[DEF_DJ2_IDX + i] = bb2[i];
    }

    setSachsBasis(bb1, bb2);

    enum_break_condition breakType = enum_break_none;

    if (fabs(testConstraint()) > mConstraintEpsilon) {
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

    int64_t t1 = get_system_clock();

    int numOK = 0, numBAD = 0; // count=0;

    double h = mLambdaStep;
    double hdid, hnext;
    double cstr;

    mNumCoords = 32;

    mFirst = 1;
    mExitFlag = 0;
    epsold = -1.0;
    eps = epsilon_abs;
    while ((int)points.size() < maxNumPoints && breakType == enum_break_none) {
        if (!calcDerivsSachsJacobi(y, dydx)) {
            breakType = enum_break_not_implemented;
        }

        for (int i = 0; i < mNumCoords; i++) {
            yscal[i] = fabs(y[i]) + fabs(dydx[i] * h) + DEF_BS_TINY;
        }

        // overshoot ignored
        breakType = nextStepSachsJacobi(h, hdid, hnext, cstr);

        if (!outsideBoundBox()) {
            points.push_back(vec4(y[0], y[1], y[2], y[3]));
            dirs.push_back(vec4(y[4], y[5], y[6], y[7]));

            sachs0.push_back(vec4(&y[DEF_SA1_IDX]));
            sachs1.push_back(vec4(&y[DEF_SA2_IDX]));

            calcJacobiParams(mLambda, y, currJacobi);
            jacobi.push_back(currJacobi);
            lambda.push_back(mLambda);
            findMaxJacobi(currJacobi, maxJacobi);

            if (fabs(mMetric->testConstraint(y, mKappa)) > mConstraintEpsilon) {
                breakType = enum_break_constraint;
            }
        }
        else {
            breakType = enum_break_outside;
        }

        if (hdid == h) {
            ++numOK;
        }
        else {
            ++numBAD;
        }

        if (fabs(hnext) <= mhmin) {
            fprintf(stderr, "Stepsize too small in GeodesicBS.\n");
        }
        h = hnext;
        if (h > mhmax) {
            h = mhmax; // Maximale Schrittweite, um "runde" Geodäten zu bekommen
        }
    }

    if ((int)points.size() >= maxNumPoints) {
        breakType = enum_break_num_exceed;
    }

    int64_t t2 = get_system_clock();
    mCalcTime = (t2 - t1) * 1e-6;
    return breakType;
}

enum_break_condition GeodesicBS::calcSachsJacobi(const vec4 initPos, const vec4 initCoordDir, const vec3 localNullDir,
    const vec3 locX, const vec3 locY, const vec3 locZ, const vec4 b0, const vec4 b1, const vec4 b2, const vec4 b3,
    const enum_nat_tetrad_type tetrad_type, const int maxNumPoints, vec4*& points, vec4*& dirs, double*& lambda,
    vec4*& sachs0, vec4*& sachs1, vec5*& jacobi, vec5& maxJacobi, int& numPoints)
{
    if (points != nullptr) {
        delete[] points;
    }
    if (dirs != nullptr) {
        delete[] dirs;
    }
    if (lambda != nullptr) {
        delete[] lambda;
    }
    if (sachs0 != nullptr) {
        delete[] sachs0;
    }
    if (sachs1 != nullptr) {
        delete[] sachs1;
    }
    if (jacobi != nullptr) {
        delete[] jacobi;
    }

    points = new vec4[maxNumPoints];
    dirs = new vec4[maxNumPoints];
    lambda = new double[maxNumPoints];
    sachs0 = new vec4[maxNumPoints];
    sachs1 = new vec4[maxNumPoints];
    jacobi = new vec5[maxNumPoints];

    double DOUBLE_MAX = std::numeric_limits<double>::max();
    maxJacobi = vec5(-DOUBLE_MAX, -DOUBLE_MAX, -DOUBLE_MAX, -DOUBLE_MAX, -DOUBLE_MAX);

    resetAffineParam();
    resetAffineParamStep();

    setInitialPosition(initPos);
    setInitialDirection(initCoordDir);

    calcSachsBasis(localNullDir, locX, locY, locZ);

    vec4 e0, e1, e2, e3;
    mMetric->localToCoord(initPos, b0, e0, tetrad_type);
    mMetric->localToCoord(initPos, b1, e1, tetrad_type);
    mMetric->localToCoord(initPos, b2, e2, tetrad_type);
    mMetric->localToCoord(initPos, b3, e3, tetrad_type);

    vec4 bb1, bb2;
    bb1 = mSachsBasisB1[0] * e1 + mSachsBasisB1[1] * e2 + mSachsBasisB1[2] * e3;
    bb2 = mSachsBasisB2[0] * e1 + mSachsBasisB2[1] * e2 + mSachsBasisB2[2] * e3;

    for (int i = 0; i < 4; i++) {
        y[DEF_JAC1_IDX + i] = 0.0;
        y[DEF_DJ1_IDX + i] = bb1[i];

        y[DEF_JAC2_IDX + i] = 0.0;
        y[DEF_DJ2_IDX + i] = bb2[i];
    }

    setSachsBasis(bb1, bb2);

    enum_break_condition breakType = enum_break_none;

    if (fabs(testConstraint()) > mConstraintEpsilon) {
        return enum_break_constraint;
    }

    vec5 currJacobi;

    points[0] = vec4(&y[0]);
    dirs[0] = vec4(&y[DEF_TG_IDX]);
    sachs0[0] = vec4(&y[DEF_SA1_IDX]);
    sachs1[0] = vec4(&y[DEF_SA2_IDX]);
    jacobi[0] = vec5(0.0, 0.0, 1.0, 0.0, 0.0);
    lambda[0] = mLambda;
    maxJacobi = jacobi[0];

    int64_t t1 = get_system_clock();

    int numOK = 0, numBAD = 0; // count=0;

    double h = mLambdaStep;
    double hdid, hnext;
    double cstr;

    mNumCoords = 32;

    mFirst = 1;
    mExitFlag = 0;
    epsold = -1.0;
    eps = epsilon_abs;

    numPoints = 1;
    while (numPoints < maxNumPoints && breakType == enum_break_none) {
        if (!calcDerivsSachsJacobi(y, dydx)) {
            breakType = enum_break_not_implemented;
        }

        for (int i = 0; i < mNumCoords; i++) {
            yscal[i] = fabs(y[i]) + fabs(dydx[i] * h) + DEF_BS_TINY;
        }

        // overshoot ignored
        breakType = nextStepSachsJacobi(h, hdid, hnext, cstr);

        if (!outsideBoundBox()) {

            points[numPoints] = vec4(&y[0]);
            dirs[numPoints] = vec4(&y[DEF_TG_IDX]);
            sachs0[numPoints] = vec4(&y[DEF_SA1_IDX]);
            sachs1[numPoints] = vec4(&y[DEF_SA2_IDX]);
            lambda[numPoints] = mLambda;

            calcJacobiParams(mLambda, y, currJacobi);
            jacobi[numPoints] = currJacobi;
            numPoints++;

            findMaxJacobi(currJacobi, maxJacobi);

            if (fabs(mMetric->testConstraint(y, mKappa)) > mConstraintEpsilon) {
                breakType = enum_break_constraint;
            }
        }
        else {
            breakType = enum_break_outside;
        }

        if (hdid == h) {
            ++numOK;
        }
        else {
            ++numBAD;
        }

        if (fabs(hnext) <= mhmin) {
            fprintf(stderr, "Stepsize too small in GeodesicBS.\n");
        }
        h = hnext;
        if (h > mhmax) {
            h = mhmax; // Maximale Schrittweite, um "runde" Geodäten zu bekommen
        }
    }

    if (numPoints == maxNumPoints) {
        breakType = enum_break_num_exceed;
    }

    int64_t t2 = get_system_clock();
    mCalcTime = (t2 - t1) * 1e-6;
    return breakType;
}

/*! Calculate next Runge-Kutta step of the geodesic.
 *
 *  \param  htry : try stepsize.
 *  \param  hdid : did stepsize.
 *  \param  hnext : next stepsize.
 *  \param  constraint : reference to result of constraint equation.
 */
enum_break_condition GeodesicBS::nextStep(double htry, double& hdid, double& hnext, double& constraint)
{
    if (mMetric->breakCondition(&y[0])) {
        return enum_break_cond;
    }
    register int i, k;
    int km = 0;
    double eps1 = 0.0, xest = 0.0, errmax = DEF_BS_TINY, red = 0.0;

    if (eps != epsold) {
        hnext = xnew = -1e29;
        eps1 = DEF_BS_SAFE1 * eps;
        a[0] = nSeq[0] + 1;
        for (k = 0; k < DEF_BS_MAX_ROW_NUM; k++) {
            a[k + 1] = a[k] + nSeq[k + 1];
        }

        for (i = 1; i < DEF_BS_MAX_ROW_NUM; i++)
            for (k = 0; k < i; k++) {
                alf[k][i] = pow(eps1, (a[k + 1] - a[i + 1]) / ((a[i + 1] - a[0] + 1.0) * (2.0 * (k + 1) + 1.0)));
                //  fprintf(stderr,"%d %d %g %f
                //  %f\n",k,i,eps1,(a[k+1]-a[i+1])/((a[i+1]-a[0]+1.0)*(2.0*(k+1)+1.0)),alf[k][i]);
            }
        epsold = eps;
        for (k = 1; k < DEF_BS_MAX_ROW_NUM - 1; k++)
            if (a[k + 1] > a[k] * alf[k - 1][k]) {
                break;
            }
        mKmax = k;
    }

    double h = htry;
    for (i = 0; i < mNumCoords; i++) {
        ysav[i] = y[i];
    }

    if (mLambda != xnew || h != (hnext)) {
        mFirst = 1;
        mKopt = mKmax;
    }

    // std::cerr << mKopt << " " << mKmax << " " << h << std::endl;

    mReduct = 0;
    for (;;) {
        for (k = 0; k < mKmax; k++) {
            xnew = mLambda + h;
            if (xnew == mLambda) {
                fprintf(stderr, "Stepsize underflow\n");
                mExitFlag = 1;
                break;
            }

            modMidPoint(ysav, dydx, h, nSeq[k], yseq);

            xest = M4D_SQR(h / (double)nSeq[k]);
            polyExtrpol(k, xest, yseq, y, yerr);
            // fprintf(stderr,"%f %f %f %f\n",yseq[0],yseq[1],yseq[2],yseq[3]);
            // exit(1);
            if (k != 0) {
                errmax = DEF_BS_TINY;
                for (i = 0; i < mNumCoords; i++) {
                    errmax = std::max(errmax, fabs(yerr[i] / yscal[i]));
                }
                errmax /= eps;
                km = k - 1;
                err[km] = pow(errmax / DEF_BS_SAFE1, 1.0 / (2.0 * km + 1.0));
            }

            if (k != 0 && (k >= mKopt - 1 || mFirst)) {
                if (errmax < 1.0) {
                    mExitFlag = 1;
                    break;
                }

                if (k == mKmax || k == mKopt + 1) {
                    red = DEF_BS_SAFE2 / err[km];
                    break;
                }
                else if (k == mKopt && alf[mKopt - 1][mKopt] < err[km]) {
                    red = 1.0 / err[km];
                    break;
                }
                else if (mKopt == mKmax && alf[km][mKmax - 1] < err[km]) {
                    red = alf[km][mKmax - 1] * DEF_BS_SAFE2 / err[km];
                    break;
                }
                else if (alf[km][mKopt] < err[km]) {
                    red = alf[km][mKopt - 1] / err[km];
                    break;
                }
            }
        }

        if (mExitFlag) {
            break;
        }

        red = std::min(red, DEF_BS_MIN_RED);
        red = std::max(red, DEF_BS_MAX_RED);

        h *= red;
        mReduct = 1;
    }

    mLambda = xnew;
    hdid = h;
    mFirst = 0;

    double wrkmin = 1e35;
    double fact, work, scale = 1.0;
    for (k = 0; k < km; k++) {
        fact = std::max(err[k], DEF_BS_MAX_SCALE);
        work = fact * a[k + 1];
        if (work < wrkmin) {
            scale = fact;
            wrkmin = work;
            mKopt = k + 1;
        }
    }

    hnext = h / scale;

    if (mKopt >= k && mKopt != mKmax && !mReduct) {
        fact = std::max(scale / alf[mKopt - 1][mKopt], DEF_BS_MAX_SCALE);
        if (a[mKopt + 1] * fact <= wrkmin) {
            hnext = h / fact;
            mKopt++;
        }
    }

    if ((constraint = fabs(mMetric->testConstraint(y, mKappa))) > mConstraintEpsilon) {
        return enum_break_constraint;
    }

    // mLambda += mLambdaStep;
    return enum_break_none;
}

/*! Calculate next BS step of the geodesic.
 *
 *  \param  htry : try stepsize.
 *  \param  hdid : did stepsize.
 *  \param  hnext : next stepsize.
 *  \param  constraint : reference to result of constraint equation.
 */
enum_break_condition GeodesicBS::nextStepSachsJacobi(double htry, double& hdid, double& hnext, double& constraint)
{
    if (mMetric->breakCondition(&y[0])) {
        return enum_break_cond;
    }
    register int i, k;
    int km = 0;
    double eps1 = 0.0, xest = 0.0, errmax = DEF_BS_TINY, red = 0.0;

    if (eps != epsold) {
        hnext = xnew = -1e29;
        eps1 = DEF_BS_SAFE1 * eps;
        a[0] = nSeq[0] + 1;
        for (k = 0; k < DEF_BS_MAX_ROW_NUM; k++) {
            a[k + 1] = a[k] + nSeq[k + 1];
        }

        for (i = 1; i < DEF_BS_MAX_ROW_NUM; i++)
            for (k = 0; k < i; k++) {
                alf[k][i] = pow(eps1, (a[k + 1] - a[i + 1]) / ((a[i + 1] - a[0] + 1.0) * (2.0 * (k + 1) + 1.0)));
                //  fprintf(stderr,"%d %d %g %f
                //  %f\n",k,i,eps1,(a[k+1]-a[i+1])/((a[i+1]-a[0]+1.0)*(2.0*(k+1)+1.0)),alf[k][i]);
            }
        epsold = eps;
        for (k = 1; k < DEF_BS_MAX_ROW_NUM - 1; k++)
            if (a[k + 1] > a[k] * alf[k - 1][k]) {
                break;
            }
        mKmax = k;
    }

    double h = htry;
    for (i = 0; i < mNumCoords; i++) {
        ysav[i] = y[i];
    }

    if (mLambda != xnew || h != (hnext)) {
        mFirst = 1;
        mKopt = mKmax;
    }

    // std::cerr << mKopt << " " << mKmax << " " << h << std::endl;

    mReduct = 0;
    for (;;) {
        for (k = 0; k < mKmax; k++) {
            xnew = mLambda + h;
            if (xnew == mLambda) {
                fprintf(stderr, "Stepsize underflow\n");
                mExitFlag = 1;
                break;
            }

            modMidPointJS(ysav, dydx, h, nSeq[k], yseq);

            xest = M4D_SQR(h / (double)nSeq[k]);
            polyExtrpol(k, xest, yseq, y, yerr);
            // fprintf(stderr,"%f %f %f %f\n",yseq[0],yseq[1],yseq[2],yseq[3]);
            // exit(1);
            if (k != 0) {
                errmax = DEF_BS_TINY;
                for (i = 0; i < mNumCoords; i++) {
                    errmax = std::max(errmax, fabs(yerr[i] / yscal[i]));
                }
                errmax /= eps;
                km = k - 1;
                err[km] = pow(errmax / DEF_BS_SAFE1, 1.0 / (2.0 * km + 1.0));
            }

            if (k != 0 && (k >= mKopt - 1 || mFirst)) {
                if (errmax < 1.0) {
                    mExitFlag = 1;
                    break;
                }

                if (k == mKmax || k == mKopt + 1) {
                    red = DEF_BS_SAFE2 / err[km];
                    break;
                }
                else if (k == mKopt && alf[mKopt - 1][mKopt] < err[km]) {
                    red = 1.0 / err[km];
                    break;
                }
                else if (mKopt == mKmax && alf[km][mKmax - 1] < err[km]) {
                    red = alf[km][mKmax - 1] * DEF_BS_SAFE2 / err[km];
                    break;
                }
                else if (alf[km][mKopt] < err[km]) {
                    red = alf[km][mKopt - 1] / err[km];
                    break;
                }
            }
        }

        if (mExitFlag) {
            break;
        }

        red = std::min(red, DEF_BS_MIN_RED);
        red = std::max(red, DEF_BS_MAX_RED);

        // std::cerr << red << std::endl;
        h *= red;
        mReduct = 1;
    }

    mLambda = xnew;
    hdid = h;
    mFirst = 0;

    double wrkmin = 1e35;
    double fact, work, scale = 1.0;
    for (k = 0; k < km; k++) {
        fact = std::max(err[k], DEF_BS_MAX_SCALE);
        work = fact * a[k + 1];
        if (work < wrkmin) {
            scale = fact;
            wrkmin = work;
            mKopt = k + 1;
        }
    }

    hnext = h / scale;

    if (mKopt >= k && mKopt != mKmax && !mReduct) {
        fact = std::max(scale / alf[mKopt - 1][mKopt], DEF_BS_MAX_SCALE);
        if (a[mKopt + 1] * fact <= wrkmin) {
            hnext = h / fact;
            mKopt++;
        }
    }

    if ((constraint = fabs(mMetric->testConstraint(y, mKappa))) > mConstraintEpsilon) {
        return enum_break_constraint;
    }

    // mLambda += mLambdaStep;
    return enum_break_none;
}

bool GeodesicBS::nextStep(int& status)
{
    double h = mLambdaStep;
    double hdid, hnext;
    double cstr;
    nextStep(h, hdid, hnext, cstr);
    status = 1;
    return true;
}

bool GeodesicBS::nextStepPar(int&)
{
    return false;
}

bool GeodesicBS::nextStepSachsJacobi(int& status)
{
    double h = mLambdaStep;
    double hdid, hnext;
    double cstr;
    nextStepSachsJacobi(h, hdid, hnext, cstr);
    status = 1;
    return true;
}

/*! Print geodesic solver properties.
 * \param fptr : file pointer.
 */
void GeodesicBS::print(FILE* fptr)
{
    fprintf(fptr, "\nGeodesicBS:\n------------\n");
    fprintf(fptr, "\tstepsize controlled : %s\n", "yes");
    fprintf(fptr, "\tstep size           : %12.8e\n", mLambdaStep);
    fprintf(fptr, "\tepsilon_abs         : %12.8e\n", epsilon_abs);
    fprintf(fptr, "\tepsilon_rel         : %12.8e\n", epsilon_rel);
    fprintf(fptr, "\tconstraint epsilon  : %12.8e\n", mConstraintEpsilon);
    fprintf(fptr, "\tmax step size       : %12.8e\n", mhmax);
    fprintf(fptr, "\tbounding box min    : %14.6e %14.6e %14.6e %14.6e\n", mBoundBoxMin[0], mBoundBoxMin[1],
        mBoundBoxMin[2], mBoundBoxMin[3]);
    fprintf(fptr, "\tbounding box max    : %14.6e %14.6e %14.6e %14.6e\n", mBoundBoxMax[0], mBoundBoxMax[1],
        mBoundBoxMax[2], mBoundBoxMax[3]);
}

// *********************************** protected methods ******************************

/*! Modified midpoint method.
 *
 * \param yy : yy
 * \param dydx : dydx.
 * \param H : total step size
 * \param numSteps : number of steps
 * \param yout : pointer to yout.
 */
void GeodesicBS::modMidPoint(double* yy, double* dydx, double H, int numSteps, double* yout)
{
    register int i, n;

    long double h = H / (long double)numSteps;

    for (i = 0; i < mNumCoords; i++) { // first step
        ym[i] = yy[i]; // z0
        yn[i] = yy[i] + h * dydx[i]; // z1
    }

    // double lambda = slambda + h;
    calcDerivs(yn, yout);

    double h2 = 2.0 * h;
    double yswap;
    for (n = 1; n < numSteps; n++) {
        for (i = 0; i < mNumCoords; i++) {
            yswap = ym[i] + h2 * yout[i]; // z_n+1
            ym[i] = yn[i]; // z_n-1 := z_n
            yn[i] = yswap; // z_n+1
        }
        // lambda += h;
        calcDerivs(yn, yout);
    }

    for (i = 0; i < mNumCoords; i++) {
        yout[i] = 0.5 * (ym[i] + yn[i] + h * yout[i]); // lambda*y[i]); //geändert
    }
}

/*! Modified midpoint method for jacobi-sachs.
 *
 * \param yy : yy
 * \param dydx : dydx.
 * \param H : total step size
 * \param numSteps : number of steps
 * \param yout : pointer to yout.
 */
void GeodesicBS::modMidPointJS(double* yy, double* dydx, double H, int numSteps, double* yout)
{
    register int i, n;

    long double h = H / (long double)numSteps;

    for (i = 0; i < mNumCoords; i++) { // first step
        ym[i] = yy[i]; // z0
        yn[i] = yy[i] + h * dydx[i]; // z1
    }

    // double lambda = slambda + h;
    calcDerivsSachsJacobi(yn, yout);

    double h2 = 2.0 * h;
    double yswap;
    for (n = 1; n < numSteps; n++) {
        for (i = 0; i < mNumCoords; i++) {
            yswap = ym[i] + h2 * yout[i]; // z_n+1
            ym[i] = yn[i]; // z_n-1 := z_n
            yn[i] = yswap; // z_n+1
        }
        // lambda += h;
        calcDerivsSachsJacobi(yn, yout);
    }

    for (i = 0; i < mNumCoords; i++) {
        yout[i] = 0.5 * (ym[i] + yn[i] + h * yout[i]); // lambda*y[i]); //geändert
    }
}

/*! Polynomial extrapolation
 */
void GeodesicBS::polyExtrpol(int iest, double xest, double* yest, double* yz, double* dy)
{
    register int k, j;
    double q, f2, f1, delta;

    x[iest] = xest;

    for (j = 0; j < mNumCoords; j++) {
        dy[j] = yz[j] = yest[j];
    }

    if (iest) {
        for (j = 0; j < mNumCoords; j++) {
            cc[j] = yest[j];
        }
        for (k = 0; k < iest; k++) {
            delta = 1.0
                / (x[iest - k - 1] - xest); // Teilen durch 0, da x[iest] = xest!!! -> geändert (iest -k) in iest -k-1
            f1 = xest * delta;
            f2 = x[iest - k - 1] * delta;
            for (j = 0; j < mNumCoords; j++) {
                q = dd[j][k]; // Ist dd initialisiert??
                dd[j][k] = dy[j];
                delta = cc[j] - q;
                dy[j] = f1 * delta;
                cc[j] = f2 * delta;
                yz[j] += dy[j];
            }
        }
    }
    for (j = 0; j < mNumCoords; j++) {
        dd[j][iest] = dy[j];
    }
}

} // end namespace m4d
