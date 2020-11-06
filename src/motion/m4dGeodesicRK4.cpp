/**
 * @file    m4dGeodesicRK4.cpp
 * @author  Thomas Mueller
 *
 *  This file is part of libMotion4D.
 */
#include "m4dGeodesicRK4.h"

namespace m4d {

GeodesicRK4::GeodesicRK4(Metric* metric, enum_geodesic_type type)
    : Geodesic(metric, type)
{
    mCalcWithParTransport = false;
    mNumCoords = 8;

    mConstraintEpsilon = DEF_CONSTRAINT_EPSILON;
}

GeodesicRK4::~GeodesicRK4() {}

enum_break_condition GeodesicRK4::calculateGeodesic(const vec4 initPos, const vec4 initDir, const int maxNumPoints,
    std::vector<vec4>& points, std::vector<vec4>& dirs, std::vector<double>& lambda)
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

    enum_break_condition breakType = enum_break_none;
    double cstr;
    breakType = initializeGeodesic(initPos, initDir, cstr);
    if (breakType == enum_break_constraint) {
        return enum_break_constraint;
    }

    points.push_back(vec4(y[0], y[1], y[2], y[3]));
    dirs.push_back(vec4(y[4], y[5], y[6], y[7]));
    lambda.push_back(mLambda);

    int status;
    int64_t t1 = get_system_clock();
    while ((int)points.size() < maxNumPoints && breakType == enum_break_none) {
        if (!nextStep(status)) {
            breakType = enum_break_cond;
        }

        if (!outsideBoundBox()) {
            points.push_back(vec4(y[0], y[1], y[2], y[3]));
            dirs.push_back(vec4(y[4], y[5], y[6], y[7]));
            lambda.push_back(mLambda);

            if ((cstr = fabs(mMetric->testConstraint(y, mKappa))) > mConstraintEpsilon) {
                breakType = enum_break_constraint;
            }
        }
        else {
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

enum_break_condition GeodesicRK4::calculateGeodesic(
    const vec4 initPos, const vec4 initDir, const int maxNumPoints, vec4*& points, vec4*& dirs, int& numPoints)
{

    if (points != nullptr) {
        delete[] points;
    }
    if (dirs != nullptr) {
        delete[] dirs;
    }

    points = new vec4[maxNumPoints];
    dirs = new vec4[maxNumPoints];

    enum_break_condition breakType = enum_break_none;
    double cstr;
    breakType = initializeGeodesic(initPos, initDir, cstr);
    if (breakType == enum_break_constraint) {
        return enum_break_constraint;
    }

    points[0] = vec4(y[0], y[1], y[2], y[3]);
    dirs[0] = vec4(y[4], y[5], y[6], y[7]);

    int status;
    int64_t t1 = get_system_clock();

    numPoints = 1;
    while (numPoints < maxNumPoints && breakType == enum_break_none) {
        if (!nextStep(status)) {
            breakType = enum_break_cond;
        }

        if (!outsideBoundBox()) {
            points[numPoints] = vec4(y[0], y[1], y[2], y[3]);
            dirs[numPoints] = vec4(y[4], y[5], y[6], y[7]);
            numPoints++;

            if ((cstr = fabs(mMetric->testConstraint(y, mKappa))) > mConstraintEpsilon) {
                breakType = enum_break_constraint;
            }
        }
        else {
            breakType = enum_break_outside;
        }
    }

    if (numPoints == maxNumPoints) {
        breakType = enum_break_num_exceed;
    }

    int64_t t2 = get_system_clock();
    mCalcTime = (t2 - t1) * 1e-6;
    return breakType;
}

enum_break_condition GeodesicRK4::calculateGeodesicData(const vec4 initPos, const vec4 initDir, const int maxNumPoints,
    std::vector<vec4>& points, std::vector<vec4>& dirs, std::vector<double>& epsilons, std::vector<double>& lambda)
{
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

    enum_break_condition breakType = enum_break_none;

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

    int status;
    int64_t t1 = get_system_clock();
    while ((int)points.size() < maxNumPoints && breakType == enum_break_none) {
        if (!nextStep(status)) {
            breakType = enum_break_cond;
        }

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
        }
        else {
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

enum_break_condition GeodesicRK4::calcParTransport(const vec4 initPos, const vec4 initDir, const vec4 e0, const vec4 e1,
    const vec4 e2, const vec4 e3, const int maxNumPoints, std::vector<vec4>& points, std::vector<vec4>& dirs,
    std::vector<double>& lambda, std::vector<vec4>& base0, std::vector<vec4>& base1, std::vector<vec4>& base2,
    std::vector<vec4>& base3)
{
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

    enum_break_condition breakType = enum_break_none;

    double cstr;
    breakType = initializeGeodesic(initPos, initDir, cstr);
    if (breakType == enum_break_constraint) {
        return enum_break_constraint;
    }

    setInitialTetrad(e0, e1, e2, e3);

    points.push_back(vec4(y[0], y[1], y[2], y[3]));
    dirs.push_back(vec4(y[4], y[5], y[6], y[7]));
    lambda.push_back(mLambda);
    base0.push_back(e0);
    base1.push_back(e1);
    base2.push_back(e2);
    base3.push_back(e3);

    int status;
    int64_t t1 = get_system_clock();
    while ((int)points.size() < maxNumPoints && breakType == enum_break_none) {
        if (!nextStepPar(status)) {
            breakType = enum_break_cond;
        }

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
        }
        else {
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

enum_break_condition GeodesicRK4::calcSachsJacobi(const vec4 initPos, const vec4 initCoordDir, const vec3 localNullDir,
    const vec3 locX, const vec3 locY, const vec3 locZ, const vec4 b0, const vec4 b1, const vec4 b2, const vec4 b3,
    const enum_nat_tetrad_type tetrad_type, const int maxNumPoints, std::vector<vec4>& points, std::vector<vec4>& dirs,
    std::vector<double>& lambda, std::vector<vec4>& sachs0, std::vector<vec4>& sachs1, std::vector<vec5>& jacobi,
    vec5& maxJacobi)
{
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

    double DOUBLE_MAX = std::numeric_limits<double>::max();
    maxJacobi = vec5(-DOUBLE_MAX, -DOUBLE_MAX, -DOUBLE_MAX, -DOUBLE_MAX, -DOUBLE_MAX);

    enum_break_condition breakType = enum_break_none;

    double cstr;
    breakType = initializeGeodesic(initPos, initCoordDir, cstr);
    if (breakType == enum_break_constraint) {
        return enum_break_constraint;
    }

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

    /*
    double prod;
    mMetric->calcProduct(initPos,initCoordDir,initCoordDir,prod);  fprintf(stderr,"%g\n",prod);
    mMetric->calcProduct(initPos,initCoordDir,bb1,prod);  fprintf(stderr,"%g\n",prod);
    mMetric->calcProduct(initPos,initCoordDir,bb2,prod);  fprintf(stderr,"%g\n",prod);  std::cerr << std::endl;
    */
    setSachsBasis(bb1, bb2);

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
    while ((int)points.size() < maxNumPoints && breakType == enum_break_none) {
        if (!nextStepSachsJacobi(status)) {
            breakType = enum_break_cond;
        }

        if (!outsideBoundBox()) {
            points.push_back(vec4(&y[0]));
            dirs.push_back(vec4(&y[DEF_TG_IDX]));
            sachs0.push_back(vec4(&y[DEF_SA1_IDX]));
            sachs1.push_back(vec4(&y[DEF_SA2_IDX]));

            calcJacobiParams(mLambda, y, currJacobi);
            jacobi.push_back(currJacobi);
            lambda.push_back(mLambda);
            findMaxJacobi(currJacobi, maxJacobi);
        }
        else {
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

enum_break_condition GeodesicRK4::calcSachsJacobi(const vec4 initPos, const vec4 initCoordDir, const vec3 localNullDir,
    const vec3 locX, const vec3 locY, const vec3 locZ, const vec4 b0, const vec4 b1, const vec4 b2, const vec4 b3,
    const enum_nat_tetrad_type tetrad_type, const int maxNumPoints, vec4*& points, vec4*& dirs, double*& lambda,
    vec4*& sachs0, vec4*& sachs1, vec5*& jacobi, vec5& maxJacobi, int& numPoints)
{
    setCalcWithParTransport(true);

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

    enum_break_condition breakType = enum_break_none;

    double cstr;
    breakType = initializeGeodesic(initPos, initCoordDir, cstr);
    if (breakType == enum_break_constraint) {
        return enum_break_constraint;
    }

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

    vec5 currJacobi;

    points[0] = vec4(&y[0]);
    dirs[0] = vec4(&y[DEF_TG_IDX]);
    sachs0[0] = vec4(&y[DEF_SA1_IDX]);
    sachs1[0] = vec4(&y[DEF_SA2_IDX]);
    jacobi[0] = vec5(0.0, 0.0, 1.0, 0.0, 0.0);
    lambda[0] = mLambda;
    maxJacobi = jacobi[0];

    int status;
    int64_t t1 = get_system_clock();

    numPoints = 1;
    while (numPoints < maxNumPoints && breakType == enum_break_none) {
        if (!nextStepSachsJacobi(status)) {
            breakType = enum_break_cond;
        }

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
        }
        else {
            breakType = enum_break_outside;
        }
    }

    if (numPoints == maxNumPoints) {
        breakType = enum_break_num_exceed;
    }

    int64_t t2 = get_system_clock();
    mCalcTime = (t2 - t1) * 1e-6;
    return breakType;
}

bool GeodesicRK4::nextStep(int& status)
{
    if (mMetric->breakCondition(&y[0])) {
        return false;
    }

    register int i;
    double yn[8], dydx[8], k1[8], k2[8], k3[8], k4[8];

    calcDerivs(y, dydx);
    for (i = 0; i < 8; i++) {
        k1[i] = mLambdaStep * dydx[i];
        yn[i] = y[i] + 0.5 * k1[i];
    }

    calcDerivs(yn, dydx);
    for (i = 0; i < 8; i++) {
        k2[i] = mLambdaStep * dydx[i];
        yn[i] = y[i] + 0.5 * k2[i];
    }

    calcDerivs(yn, dydx);
    for (i = 0; i < 8; i++) {
        k3[i] = mLambdaStep * dydx[i];
        yn[i] = y[i] + k3[i];
    }

    calcDerivs(yn, dydx);
    for (i = 0; i < 8; i++) {
        k4[i] = mLambdaStep * dydx[i];
        yn[i] = y[i] + 1.0 / 6.0 * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);
        y[i] = yn[i];
    }

    mLambda += mLambdaStep;
    status = 1;
    return true;
}

bool GeodesicRK4::nextStepPar(int& status)
{
    if (mMetric->breakCondition(&y[0])) {
        return false;
    }

    register int i;
    double yn[DEF_MAX_YS_PAR], dydx[DEF_MAX_YS_PAR], k1[DEF_MAX_YS_PAR], k2[DEF_MAX_YS_PAR], k3[DEF_MAX_YS_PAR],
        k4[DEF_MAX_YS_PAR];

    calcDerivsPar(y, dydx);
    for (i = 0; i < DEF_MAX_YS_PAR; i++) {
        k1[i] = mLambdaStep * dydx[i];
        yn[i] = y[i] + 0.5 * k1[i];
    }

    calcDerivsPar(yn, dydx);
    for (i = 0; i < DEF_MAX_YS_PAR; i++) {
        k2[i] = mLambdaStep * dydx[i];
        yn[i] = y[i] + 0.5 * k2[i];
    }

    calcDerivsPar(yn, dydx);
    for (i = 0; i < DEF_MAX_YS_PAR; i++) {
        k3[i] = mLambdaStep * dydx[i];
        yn[i] = y[i] + k3[i];
    }

    calcDerivsPar(yn, dydx);
    for (i = 0; i < DEF_MAX_YS_PAR; i++) {
        k4[i] = mLambdaStep * dydx[i];
        yn[i] = y[i] + 1.0 / 6.0 * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);
        y[i] = yn[i];
    }

    mLambda += mLambdaStep;
    status = 1;
    return true;
}

bool GeodesicRK4::nextStepSachsJacobi(int& status)
{
    if (mMetric->breakCondition(&y[0])) {
        return false;
    }
    register int i;
    double yn[DEF_MAX_YS_JAC], dydx[DEF_MAX_YS_JAC], k1[DEF_MAX_YS_JAC], k2[DEF_MAX_YS_JAC], k3[DEF_MAX_YS_JAC],
        k4[DEF_MAX_YS_JAC];

    if (!calcDerivsSachsJacobi(y, dydx)) {
        return false; // enum_break_not_implemented;
    }

    for (i = 0; i < DEF_MAX_YS_JAC; i++) {
        k1[i] = mLambdaStep * dydx[i];
        yn[i] = y[i] + 0.5 * k1[i];
    }

    calcDerivsSachsJacobi(yn, dydx);
    for (i = 0; i < DEF_MAX_YS_JAC; i++) {
        k2[i] = mLambdaStep * dydx[i];
        yn[i] = y[i] + 0.5 * k2[i];
    }

    calcDerivsSachsJacobi(yn, dydx);
    for (i = 0; i < DEF_MAX_YS_JAC; i++) {
        k3[i] = mLambdaStep * dydx[i];
        yn[i] = y[i] + k3[i];
    }

    calcDerivsSachsJacobi(yn, dydx);
    for (i = 0; i < DEF_MAX_YS_JAC; i++) {
        k4[i] = mLambdaStep * dydx[i];
        yn[i] = y[i] + 1.0 / 6.0 * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);
        y[i] = yn[i];
    }

    mLambda += mLambdaStep;
    status = 1;
    return true;
}

void GeodesicRK4::printF(FILE* fptr)
{
    fprintf(fptr, "\nGeodesicRK4:\n------------\n");
    fprintf(fptr, "\tstepsize controlled : %s\n", "no");
    fprintf(fptr, "\tstep size           : %12.8e\n", mLambdaStep);
    fprintf(fptr, "\tepsilon_abs         : %12.8e\n", epsilon_abs);
    fprintf(fptr, "\tepsilon_rel         : %12.8e\n", epsilon_rel);
    fprintf(fptr, "\tconstraint epsilon  : %12.8e\n", mConstraintEpsilon);
    fprintf(fptr, "\tbounding box min    : %14.6e %14.6e %14.6e %14.6e\n", mBoundBoxMin[0], mBoundBoxMin[1],
        mBoundBoxMin[2], mBoundBoxMin[3]);
    fprintf(fptr, "\tbounding box max    : %14.6e %14.6e %14.6e %14.6e\n", mBoundBoxMax[0], mBoundBoxMax[1],
        mBoundBoxMax[2], mBoundBoxMax[3]);
}

} // end namespace m4d
