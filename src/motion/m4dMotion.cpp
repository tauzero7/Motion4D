/**
 * @file    m4dMotion.cpp
 * @author  Thomas Mueller
 *
 *  This file is part of libMotion4D.
 */
#include "m4dMotion.h"

namespace m4d {

Motion::Motion(Metric* metric)
{
    assert(metric != nullptr);
    mMetric = metric;

    mLambda = 0.0;
    mLambdaStep = mLambdaStepInit = 0.01;

    double DOUBLE_MAX = std::numeric_limits<double>::max();
    for (int i = 0; i < 4; i++) {
        mBoundBoxMin[i] = -DOUBLE_MAX;
        mBoundBoxMax[i] = DOUBLE_MAX;
    }

    for (int i = 0; i < DEF_MAX_YS; i++) {
        y[i] = 0.0;
    }

    mConstraintEpsilon = DEF_CONSTRAINT_EPSILON;
    mMaxLambdaStep = DEF_MAX_STEPSIZE;
    mMinLambdaStep = DEF_MIN_STEPSIZE;
}

Motion::~Motion() {}

bool Motion::setInitialPosition(double ipos[4])
{
    for (int i = 0; i < 4; i++) {
        y[i] = ipos[i];
    }
    mMetric->calculateMetric(ipos);
    return mMetric->breakCondition(ipos);
}

bool Motion::setInitialPosition(vec4 ipos)
{
    for (int i = 0; i < 4; i++) {
        y[i] = ipos[i];
    }

    mMetric->calculateMetric(ipos);
    return mMetric->breakCondition(ipos);
}

bool Motion::setInitialDirection(double d0, double d1, double d2, double d3)
{
    y[4] = d0;
    y[5] = d1;
    y[6] = d2;
    y[7] = d3;
    if (fabs(testConstraint()) < mConstraintEpsilon) {
        return true;
    }

    return false;
}

bool Motion::setInitialDirection(vec4 dir)
{
    return setInitialDirection(dir.x(0), dir.x(1), dir.x(2), dir.x(3));
}

void Motion::setInitialTetrad(double e0[4], double e1[4], double e2[4], double e3[4])
{
    for (int i = 0; i < 4; i++) {
        y[DEF_EO_IDX + i] = e0[i];
        y[DEF_E1_IDX + i] = e1[i];
        y[DEF_E2_IDX + i] = e2[i];
        y[DEF_E3_IDX + i] = e3[i];
    }
}

void Motion::setInitialTetrad(vec4 e0, vec4 e1, vec4 e2, vec4 e3)
{
    for (int i = 0; i < 4; i++) {
        y[DEF_EO_IDX + i] = e0[i];
        y[DEF_E1_IDX + i] = e1[i];
        y[DEF_E2_IDX + i] = e2[i];
        y[DEF_E3_IDX + i] = e3[i];
    }
}

bool Motion::setLocalInitialDir(double, double, double, double)
{
    fprintf(stderr,
        "Motion::setLocalInitialDir  ( double d1, double d2, double d3, double beta ) ... not implemented yet!\n");
    return false;
}

bool Motion::setLocalInitialDir(vec3, double)
{
    fprintf(stderr, "Motion::setLocalInitialDir  ( vec3 dir, double beta ) ... not implemented yet!\n");
    return false;
}

void Motion::setMetric(Metric* metric)
{
    assert(metric != nullptr);
    mMetric = metric;
}

Metric* Motion::getMetric()
{
    return mMetric;
}

bool Motion::setParam(std::string paramName, double val)
{
    if (paramName.compare("lambda") == 0) {
        mLambda = val;
        return true;
    }
    else if (paramName.compare("minLambdaStep") == 0) {
        mMinLambdaStep = val;
        return true;
    }
    else if (paramName.compare("maxLambdaStep") == 0) {
        mMaxLambdaStep = val;
        return true;
    }
    else if (paramName.compare("constr_eps") == 0) {
        mConstraintEpsilon = val;
        return true;
    }
    return false;
}

bool Motion::setParam(std::string paramName, double v0, double v1, double v2, double v3)
{
    /*
    if (paramName.compare("init_position")==0) {
        return setInitialPosition(m4d::vec4(v0, v1, v2, v3));
    }
    else if (paramName.compare("init_direction")==0) {
        return setInitialDirection(m4d::vec4(v0, v1, v2, v3));
    }
    else if (paramName.compare("local_init_direction")==0) {
        return setLocalInitialDir(v0, v1, v2, v3);
    }
    else
    */
    if (paramName.compare("lower_bb") == 0) {
        mBoundBoxMin[0] = v0;
        mBoundBoxMin[1] = v1;
        mBoundBoxMin[2] = v2;
        mBoundBoxMin[3] = v3;
        return true;
    }
    else if (paramName.compare("upper_bb") == 0) {
        mBoundBoxMax[0] = v0;
        mBoundBoxMax[1] = v1;
        mBoundBoxMax[2] = v2;
        mBoundBoxMax[3] = v3;
        return true;
    }
    return false;
}

void Motion::getPosition(double* p)
{
    for (int i = 0; i < 4; i++) {
        p[i] = y[i];
    }
}

vec4 Motion::getPosition()
{
    return vec4(y[0], y[1], y[2], y[3]);
}

void Motion::getDirection(double* p)
{
    for (int i = 0; i < 4; i++) {
        p[i] = y[i + 4];
    }
}

vec4 Motion::getDirection()
{
    return vec4(y[4], y[5], y[6], y[7]);
}

void Motion::resetAffineParam(double pt)
{
    mLambda = pt;
}

double Motion::getAffineParam()
{
    return mLambda;
}

void Motion::setAffineParamStep(double step)
{
    mLambdaStep = step;
    mLambdaStepInit = step;
}

double Motion::getAffineParamStep()
{
    return mLambdaStep;
}

void Motion::resetAffineParamStep()
{
    mLambdaStep = mLambdaStepInit;
}

void Motion::setMaxAffineParamStep(double step)
{
    mMaxLambdaStep = step;
}

double Motion::getMaxAffineParamStep()
{
    return mMaxLambdaStep;
}

void Motion::setMinAffineParamStep(double step)
{
    mMinLambdaStep = step;
}

double Motion::getMinAffineParamStep()
{
    return mMinLambdaStep;
}

void Motion::setE(unsigned int i, vec4 ee)
{
    for (int j = 0; j < 4; j++) {
        y[i * 4 + 8 + j] = ee[j];
    }
}

vec4 Motion::getE(unsigned int i)
{
    if (i > 3) {
        return vec4(0.0, 0.0, 0.0, 0.0);
    }

    vec4 e;
    for (int j = 0; j < 4; j++) {
        e[j] = y[i * 4 + 8 + j];
    }
    return e;
}

void Motion::getE0(double* p)
{
    for (int i = 0; i < 4; i++) {
        p[i] = y[DEF_EO_IDX + i];
    }
}

void Motion::getE1(double* p)
{
    for (int i = 0; i < 4; i++) {
        p[i] = y[DEF_E1_IDX + i];
    }
}

void Motion::getE2(double* p)
{
    for (int i = 0; i < 4; i++) {
        p[i] = y[DEF_E2_IDX + i];
    }
}

void Motion::getE3(double* p)
{
    for (int i = 0; i < 4; i++) {
        p[i] = y[DEF_E3_IDX + i];
    }
}

void Motion::getTetrad(double* e0, double* e1, double* e2, double* e3)
{
    for (int i = 0; i < 4; i++) {
        e0[i] = y[DEF_EO_IDX + i];
        e1[i] = y[DEF_EO_IDX + i];
        e2[i] = y[DEF_EO_IDX + i];
        e3[i] = y[DEF_EO_IDX + i];
    }
}

void Motion::getTetrad(vec4& e0, vec4& e1, vec4& e2, vec4& e3)
{
    e0.set(&y[8]);
    e1.set(&y[12]);
    e2.set(&y[16]);
    e3.set(&y[20]);
}

void Motion::getTetrad(mat4& m)
{
    vec4 e0, e1, e2, e3;
    getTetrad(e0, e1, e2, e3);
    m.setCol(0, e0);
    m.setCol(1, e1);
    m.setCol(2, e2);
    m.setCol(3, e3);
}

void Motion::getTetradInv(mat4& m)
{
    getTetrad(m);
    m.invert();
}

void Motion::getTetrad(float* m)
{
    assert(m != nullptr);
    vec4 e0, e1, e2, e3;
    getTetrad(e0, e1, e2, e3);

    for (int i = 0; i < 4; i++) {
        m[4 * i + 0] = float(e0[i]);
        m[4 * i + 1] = float(e1[i]);
        m[4 * i + 2] = float(e2[i]);
        m[4 * i + 3] = float(e3[i]);
    }
}

void Motion::getTetradInv(float* m)
{
    assert(m != nullptr);
    mat4 matrix;
    getTetradInv(matrix);

    matrix.getFloatArray(m);
}

vec4 Motion::coordToLocal(vec4 cv)
{
    mat4 imat;
    getTetradInv(imat);

    return (imat * cv);
}

bool Motion::isOrthonormal()
{
    bool itIs = true;
    double prod;

    int i = 0, j = 0;
    while (i < 4 && itIs) {
        mMetric->calcProduct(&y[0], &y[8 + 4 * i], &y[8 + 4 * j], prod);
        if (fabs(prod - eta(i, j)) < 1.0e-6) {
            itIs = true;
        }
        else {
            itIs = false;
        }
        j++;
        if (j >= 4) {
            j = 0;
            i++;
        }
    }

    return itIs;
}

void Motion::setBoundingBox(double p1[4], double p2[4])
{
    for (int i = 0; i < 4; i++) {
        if (p1[i] < p2[i]) {
            mBoundBoxMin[i] = p1[i];
            mBoundBoxMax[i] = p2[i];
        }
        else {
            mBoundBoxMin[i] = p2[i];
            mBoundBoxMax[i] = p1[i];
        }
    }
}

void Motion::setBoundingBox(vec4 p1, vec4 p2)
{
    for (int i = 0; i < 4; i++) {
        if (p1[i] < p2[i]) {
            mBoundBoxMin[i] = p1.x(i);
            mBoundBoxMax[i] = p2.x(i);
        }
        else {
            mBoundBoxMin[i] = p2.x(i);
            mBoundBoxMax[i] = p1.x(i);
        }
    }
}

void Motion::getBoundingBox(vec4& p1, vec4& p2)
{
    for (int i = 0; i < 4; i++) {
        p1.setX(i, mBoundBoxMin[i]);
        p2.setX(i, mBoundBoxMax[i]);
    }
}

bool Motion::outsideBoundBox()
{
    return ((y[0] < mBoundBoxMin[0]) || (y[0] > mBoundBoxMax[0]) || (y[1] < mBoundBoxMin[1]) || (y[1] > mBoundBoxMax[1])
        || (y[2] < mBoundBoxMin[2]) || (y[2] > mBoundBoxMax[2]) || (y[3] < mBoundBoxMin[3])
        || (y[3] > mBoundBoxMax[3]));
}

void Motion::setConstrEps(double eps)
{
    mConstraintEpsilon = fabs(eps);
}

double Motion::getConstrEps()
{
    return mConstraintEpsilon;
}

bool Motion::isRightHanded()
{
    gsl_matrix* m = gsl_matrix_alloc(4, 4);

    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++) {
            gsl_matrix_set(m, i, j, y[8 + i * 4 + j]);
        }

    int signum;
    gsl_permutation* p = gsl_permutation_alloc(4);

    gsl_linalg_LU_decomp(m, p, &signum);

    double det = gsl_linalg_LU_det(m, signum);

    gsl_permutation_free(p);
    gsl_matrix_free(m);

    return ((det > 0.0) ? true : false);
}

bool Motion::gramSchmidtOrth()
{
    vec4 pos = vec4(y[0], y[1], y[2], y[3]);
    vec4 e0 = vec4(y[8], y[9], y[10], y[11]);
    vec4 u1 = vec4(y[12], y[13], y[14], y[15]);
    vec4 u2 = vec4(y[16], y[17], y[18], y[19]);
    vec4 u3 = vec4(y[20], y[21], y[22], y[23]);

    vec4 e1, e2, e3;

    double prod;
    mMetric->calcProduct(pos, e0, u1, prod);
    e1 = u1 + prod * e0;
    e1.normalize();

    mMetric->calcProduct(pos, e0, u2, prod);
    e2 = u2 + prod * e0;
    e2.normalize();

    mMetric->calcProduct(pos, e0, u3, prod);
    e3 = u3 + prod * e0;
    e3.normalize();

    for (int i = 0; i < 4; i++) {
        y[8 + i] = e0[i];
        y[12 + i] = e1[i];
        y[16 + i] = e2[i];
        y[20 + i] = e3[i];
    }

    return true;
}

double Motion::testConstraint()
{
    // fprintf(stderr,"Motion::testConstraint ( ) ... does nothing here! Only implemented in m4dGeodesic...\n");
    return 0.0;
}

void Motion::getCurrentArray(double* cy)
{
    if (cy == nullptr) {
        return;
    }
    for (unsigned int i = 0; i < DEF_MAX_YS; i++) {
        cy[i] = y[i];
    }
}

void Motion::printTetrad(FILE* fptr)
{
    for (int i = 0; i < 4; i++) {
        fprintf(fptr, "e[%d]: ", i);
        for (int j = 0; j < 4; j++) {
            fprintf(fptr, "%14.10f ", y[8 + 4 * i + j]);
        }
        fprintf(fptr, "\n");
    }
}

} // end namespace m4d
