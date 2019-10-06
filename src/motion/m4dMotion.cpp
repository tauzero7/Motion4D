// -------------------------------------------------------------------------------
/*
    m4dMotion.cpp

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

#include "m4dMotion.h"

namespace m4d {

/*!  Standard constructor.
 *
 *  \param  metric : pointer to metric.
 */
Motion::Motion(Metric* metric)
{
    assert(metric != NULL);
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

// *********************************** public methods ******************************
/*! Set initial position in coordinates.
 *
 * \param ipos : pointer to initial position.
 * \return true : position is valid.
 * \return false : position is not valid.
 */
bool Motion::setInitialPosition(double ipos[4])
{
    for (int i = 0; i < 4; i++) {
        y[i] = ipos[i];
    }
    mMetric->calculateMetric(ipos);
    return mMetric->breakCondition(ipos);
}

/*! Set initial position in coordinates.
 *
 * \param ipos : initial position.
 * \return true : position is valid.
 * \return false : position is not valid.
 */
bool Motion::setInitialPosition(vec4 ipos)
{
    for (int i = 0; i < 4; i++) {
        y[i] = ipos[i];
    }

    mMetric->calculateMetric(ipos);
    return mMetric->breakCondition(ipos);
}

/*! Set the initial direction in coordinates.
 *
 *  \param d0: 0-component of initial direction.
 *  \param d1: 1-component of initial direction.
 *  \param d2: 2-component of initial direction.
 *  \param d3: 3-component of initial direction.
 *  \return true : direction fulfills constraint equation.
 *  \return false : direction does not fulfill constraint equation.
 */
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

/*! Set the initial direction in coordinates.
 *
 *  \param dir: initial direction.
 *  \return true : direction fulfills constraint equation.
 *  \return false : direction does not fulfill constraint equation.
 */
bool Motion::setInitialDirection(vec4 dir)
{
    return setInitialDirection(dir.x(0), dir.x(1), dir.x(2), dir.x(3));
}

/*!  Set the initial tetrad of the observer.
 *
 *  \param e0:  pointer to first tetrad vector.
 *  \param e1:  pointer to second tetrad vector.
 *  \param e2:  pointer to third tetrad vector.
 *  \param e3:  pointer to fourth tetrad vector.
 */
void Motion::setInitialTetrad(double e0[4], double e1[4], double e2[4], double e3[4])
{
    for (int i = 0; i < 4; i++) {
        y[DEF_EO_IDX + i] = e0[i];
        y[DEF_E1_IDX + i] = e1[i];
        y[DEF_E2_IDX + i] = e2[i];
        y[DEF_E3_IDX + i] = e3[i];
    }
}

/*! Set the initial tetrad of the observer.
 *
 *  \param e0: first tetrad vector.
 *  \param e1: second tetrad vector.
 *  \param e2: third tetrad vector.
 *  \param e3: fourth tetrad vector.
 */
void Motion::setInitialTetrad(vec4 e0, vec4 e1, vec4 e2, vec4 e3)
{
    for (int i = 0; i < 4; i++) {
        y[DEF_EO_IDX + i] = e0[i];
        y[DEF_E1_IDX + i] = e1[i];
        y[DEF_E2_IDX + i] = e2[i];
        y[DEF_E3_IDX + i] = e3[i];
    }
}

/*! Set the initial direction with respect to the local initial tetrad.
 *
 * \param d1 : 1-component of local initial direction.
 * \param d2 : 2-component of local initial direction.
 * \param d3 : 3-component of local initial direction.
 * \param beta : initial velocity.
 * \return false : always.
 */
bool Motion::setLocalInitialDir(double, double, double, double)
{
    fprintf(stderr,
        "Motion::setLocalInitialDir  ( double d1, double d2, double d3, double beta ) ... not implemented yet!\n");
    return false;
}

/*! Set the initial direction with respect to the local initial tetrad.
 *
 * \param dir : local initial direction.
 * \param beta : initial velocity.
 * \return false : always.
 */
bool Motion::setLocalInitialDir(vec3, double)
{
    fprintf(stderr, "Motion::setLocalInitialDir  ( vec3 dir, double beta ) ... not implemented yet!\n");
    return false;
}

/*! Set the metric.
 *
 * \param metric : pointer to metric.
 */
void Motion::setMetric(Metric* metric)
{
    assert(metric != NULL);
    mMetric = metric;
}

/*! Get the metric.
 *
 * \return pointer to metric.
 */
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

/*! Get the current position.
 *
 * \param p : pointer to position.
 */
void Motion::getPosition(double* p)
{
    for (int i = 0; i < 4; i++) {
        p[i] = y[i];
    }
}

/*! Get the current position.
 *
 *  \return position vector.
 */
vec4 Motion::getPosition()
{
    return vec4(y[0], y[1], y[2], y[3]);
}

/*! Get the current direction.
 *
 * \param p : pointer to direction.
 */
void Motion::getDirection(double* p)
{
    for (int i = 0; i < 4; i++) {
        p[i] = y[i + 4];
    }
}

/*! Get the current direction.
 *
 *  \return direction vector.
 */
vec4 Motion::getDirection()
{
    return vec4(y[4], y[5], y[6], y[7]);
}

/*! Reset the affine parameter (proper time).
 *
 *  \param  pt : reset affine parameter to value pt.
 */
void Motion::resetAffineParam(double pt)
{
    mLambda = pt;
}

/*! Get the affine parameter (proper time).
 *
 * \return affine parameter.
 */
double Motion::getAffineParam()
{
    return mLambda;
}

/*! Set the affine parameter step (proper time step).
 *
 *  \param step : new affine parameter stepsize.
 */
void Motion::setAffineParamStep(double step)
{
    mLambdaStep = step;
    mLambdaStepInit = step;
}

/*! Get the affine parameter step (proper time step).
 *
 *  \return affine parameter stepsize.
 */
double Motion::getAffineParamStep()
{
    return mLambdaStep;
}

/*! Reset affine parameter stepsize.
 */
void Motion::resetAffineParamStep()
{
    mLambdaStep = mLambdaStepInit;
}

/*! Set maxmimum affine parameter stepsize.
 */
void Motion::setMaxAffineParamStep(double step)
{
    mMaxLambdaStep = step;
}

/*! Get maximum affine parameter stepsize.
 */
double Motion::getMaxAffineParamStep()
{
    return mMaxLambdaStep;
}

/*! Set minimum affine parameter stepsize.
 */
void Motion::setMinAffineParamStep(double step)
{
    mMinLambdaStep = step;
}

/*! Get minimum affine parameter stepsize.
 */
double Motion::getMinAffineParamStep()
{
    return mMinLambdaStep;
}

/*! Set i-th tetrad vector.
 *
 *  param i : number of tetrad vector.
 *  param ee : new tetrad vector.
 */
void Motion::setE(unsigned int i, vec4 ee)
{
    for (int j = 0; j < 4; j++) {
        y[i * 4 + 8 + j] = ee[j];
    }
}

/*! Get i-th tetrad vector.
 *
 *  \param i : tetrad index.
 *  \return tetrad vector i.
 */
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

/*! Get 0-th tetrad vector.
 *
 *  \param p : pointer to e0 tetrad vector.
 */
void Motion::getE0(double* p)
{
    for (int i = 0; i < 4; i++) {
        p[i] = y[DEF_EO_IDX + i];
    }
}

/*! Get 1-th tetrad vector.
 *
 *  \param p : pointer to e1 tetrad vector.
 */
void Motion::getE1(double* p)
{
    for (int i = 0; i < 4; i++) {
        p[i] = y[DEF_E1_IDX + i];
    }
}

/*! Get 2-th tetrad vector.
 *
 *  \param p : pointer to e2 tetrad vector.
 */
void Motion::getE2(double* p)
{
    for (int i = 0; i < 4; i++) {
        p[i] = y[DEF_E2_IDX + i];
    }
}

/*! Get 3-th tetrad vector.
 *
 *  \param p : pointer to e3 tetrad vector.
 */
void Motion::getE3(double* p)
{
    for (int i = 0; i < 4; i++) {
        p[i] = y[DEF_E3_IDX + i];
    }
}

/*! Get all tetrad vectors.
 *
 *  \param e0 : pointer to e0 tetrad vector.
 *  \param e1 : pointer to e1 tetrad vector.
 *  \param e2 : pointer to e2 tetrad vector.
 *  \param e3 : pointer to e3 tetrad vector.
 */
void Motion::getTetrad(double* e0, double* e1, double* e2, double* e3)
{
    for (int i = 0; i < 4; i++) {
        e0[i] = y[DEF_EO_IDX + i];
        e1[i] = y[DEF_EO_IDX + i];
        e2[i] = y[DEF_EO_IDX + i];
        e3[i] = y[DEF_EO_IDX + i];
    }
}

/*! Get all tetrad vectors.
 *
 *  \param e0 : reference to e0 tetrad vector.
 *  \param e1 : reference to e1 tetrad vector.
 *  \param e2 : reference to e2 tetrad vector.
 *  \param e3 : reference to e3 tetrad vector.
 */
void Motion::getTetrad(vec4& e0, vec4& e1, vec4& e2, vec4& e3)
{
    e0.set(&y[8]);
    e1.set(&y[12]);
    e2.set(&y[16]);
    e3.set(&y[20]);
}

/*! Get all tetrad vectors as one matrix.
 *
 *  \param m : reference to tetrad matrix.
 */
void Motion::getTetrad(mat4& m)
{
    vec4 e0, e1, e2, e3;
    getTetrad(e0, e1, e2, e3);
    m.setCol(0, e0);
    m.setCol(1, e1);
    m.setCol(2, e2);
    m.setCol(3, e3);
}

/*! Get all tetrad vectors as one matrix (inverse).
 *
 *  \param m : reference to inverse tetrad matrix.
 */
void Motion::getTetradInv(mat4& m)
{
    getTetrad(m);
    m.invert();
}

/*! Get tetrad as matrix.
 *
 *  \param m : pointer to tetrad matrix.
 */
void Motion::getTetrad(float* m)
{
    assert(m != NULL);
    vec4 e0, e1, e2, e3;
    getTetrad(e0, e1, e2, e3);

    for (int i = 0; i < 4; i++) {
        m[4 * i + 0] = float(e0[i]);
        m[4 * i + 1] = float(e1[i]);
        m[4 * i + 2] = float(e2[i]);
        m[4 * i + 3] = float(e3[i]);
    }
}

/*! Get inverse tetrad as matrix.
 *
 *  \param m : pointer to inverse tetrad matrix.
 */
void Motion::getTetradInv(float* m)
{
    assert(m != NULL);
    mat4 matrix;
    getTetradInv(matrix);

    matrix.getFloatArray(m);
}

/*! Transform coordinate vector to local vecor.
 *
 * \param cv : coordinate vector.
 * \return local vector.
 */
vec4 Motion::coordToLocal(vec4 cv)
{
    mat4 imat;
    getTetradInv(imat);

    return (imat * cv);
}

/*! Test whether tetrad is orthonormal.
 *
 *  \return true : local tetrad is orthonormalized.
 *  \return false : local tetrad is not orthonormal.
 */
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

/*! The bounding box represents the domain for integration.
 *
 * \param p1 : pointer to bounding box.
 * \param p2 : pointer to bounding box.
 */
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

/*! The bounding box represents the domain for integration.
 *
 * \param p1 : vector of bounding box.
 * \param p2 : vector of bounding box.
 */
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

/*! The bounding box represents the domain for integration.
 *
 *  \param p1 : reference to bounding box vector.
 *  \param p2 : reference to bounding box vector.
 */
void Motion::getBoundingBox(vec4& p1, vec4& p2)
{
    for (int i = 0; i < 4; i++) {
        p1.setX(i, mBoundBoxMin[i]);
        p2.setX(i, mBoundBoxMax[i]);
    }
}

/*! Checks whether the current position 'y' is outside the bounding box.
 *
 *  \return true : point y is outside bounding box.
 *  \return false : point lies inside bounding box.
 */
bool Motion::outsideBoundBox()
{
    return ((y[0] < mBoundBoxMin[0]) || (y[0] > mBoundBoxMax[0]) || (y[1] < mBoundBoxMin[1]) || (y[1] > mBoundBoxMax[1])
        || (y[2] < mBoundBoxMin[2]) || (y[2] > mBoundBoxMax[2]) || (y[3] < mBoundBoxMin[3])
        || (y[3] > mBoundBoxMax[3]));
}

/*! Set epsilon for the constraint equation.
 *
 *  \param eps : constraint epsilon.
 */
void Motion::setConstrEps(double eps)
{
    mConstraintEpsilon = fabs(eps);
}

/*! Get epsilon for the constraint equation.
 */
double Motion::getConstrEps()
{
    return mConstraintEpsilon;
}

/*! Tests, whether the local tetrad is right handed.
 *
 * \return true : local tetrad is right handed.
 * \return false : local tetrad is not right handed.
 */
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

/*! Do a Gram-Schmidt orthonormalization.
 *
 *  Note that the base vectors will be normalized with respect to the euclidean norm.
 *  Thus, the orthonormalization can be used only for local tetrads that are given
 *  with respect to a natural local tetrad!
 *
 *  \return true : always.
 */
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

/*! Test constraint condition.
 *
 * \return 0.0
 */
double Motion::testConstraint()
{
    // fprintf(stderr,"Motion::testConstraint ( ) ... does nothing here! Only implemented in m4dGeodesic...\n");
    return 0.0;
}

/*! Get current array y.
 *
 * \param cy : pointer to array (must be of size DEF_MAX_YS).
 */
void Motion::getCurrentArray(double* cy)
{
    if (cy == NULL) {
        return;
    }
    for (unsigned int i = 0; i < DEF_MAX_YS; i++) {
        cy[i] = y[i];
    }
}

/*! Print local tetrad.
 *
 * \param fptr : file pointer.
 */
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
