/**
 * @file    m4dMotion.h
 * @author  Thomas Mueller
 *
 *   Base class for motion/geodesics in 4d spacetimes.
 *
             The Motion class cannot be used itself due to the pure virtual
             method calcDerivs().

             The initial position as well as the initial direction have to
             be given in coordinates of the corresponding metric. For the
             initial direction, the localToCoord() method of the corresponding
             Metric class can be used. The initial position is checked whether it
             violates the break condition of the metric.

             For internal calculations, the current position, direction, and
             base vectors are mapped to the double array y:

            \verbatim
               y[ 0...3]  :  current position,
               y[ 4...7]  :  current direction,
               y[ 8..11]  :  e0,
               y[12..15]  :  e1,
               y[16..19]  :  e2.
               y[20..23]  :  e3. \endverbatim

 *  This file is part of libMotion4D.
 */
#ifndef M4D_MOTION_H
#define M4D_MOTION_H

#include <cassert>
#include <cmath>
#include <iostream>

#include <extra/m4dUtilities.h>
#include <m4dGlobalDefs.h>
#include <metric/m4dMetric.h>

#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>

#ifdef _WIN32
#ifndef __GNUC__
#pragma warning(disable : 4244)
#endif
#endif

namespace m4d {

/**
 * @brief The Motion class
 */
class API_M4D_EXPORT Motion
{
public:
    Motion(Metric* metric);
    virtual ~Motion();

    /**
     * @brief Set initial position in coordinates.
     * @return false if position is not valid.
     */
    virtual bool setInitialPosition(double ipos[4]);
    virtual bool setInitialPosition(vec4 ipos);

    /**
     * @brief Set the initial direction in coordinates.
     *  @param d0: 0-component of initial direction.
     *  @param d1: 1-component of initial direction.
     *  @param d2: 2-component of initial direction.
     *  @param d3: 3-component of initial direction.
     *  @return false if direction does not fulfill constraint equation.
     */
    virtual bool setInitialDirection(double d0, double d1, double d2, double d3);

    /**
     * @brief Set the initial direction in coordinates.
     *  @param dir: initial direction.
     *  @return false if direction does not fulfill constraint equation.
     */
    virtual bool setInitialDirection(vec4 dir);

    /**
     * @brief Set the initial tetrad of the observer.
     *  @param e0:  pointer to first tetrad vector.
     *  @param e1:  pointer to second tetrad vector.
     *  @param e2:  pointer to third tetrad vector.
     *  @param e3:  pointer to fourth tetrad vector.
     */
    virtual void setInitialTetrad(double e0[4], double e1[4], double e2[4], double e3[4]);
    virtual void setInitialTetrad(vec4 e0, vec4 e1, vec4 e2, vec4 e3);

    /**
     * @brief Set the initial direction with respect to the local initial tetrad.
     * @param d1 : 1-component of local initial direction.
     * @param d2 : 2-component of local initial direction.
     * @param d3 : 3-component of local initial direction.
     * @param beta : initial velocity.
     */
    virtual bool setLocalInitialDir(double d1, double d2, double d3, double beta = 1.0);

    /**
     * @brief Set the initial direction with respect to the local initial tetrad.
     * @param dir : local initial direction.
     * @param beta : initial velocity.
     */
    virtual bool setLocalInitialDir(vec3 dir, double beta = 1.0);

    virtual void setMetric(Metric* metric);
    virtual Metric* getMetric();

    virtual bool setParam(std::string paramName, double val);
    virtual bool setParam(std::string paramName, double v0, double v1, double v2, double v3);

    /**
     * @brief Get the current position.
     * @param p : pointer to position.
     */
    void getPosition(double* p);

    /// Get the current position.
    vec4 getPosition();

    /**
     * @brief Get the current direction.
     * @param p : pointer to direction.
     */
    void getDirection(double* p);

    /// Get the current direction.
    vec4 getDirection();

    /**
     * @brief Reset the affine parameter (proper time).
     * @param  pt : reset affine parameter to value pt.
     */
    void resetAffineParam(double pt = 0.0);

    /// Get the affine parameter (proper time).
    double getAffineParam();

    /// Set the affine parameter step (proper time step).
    void setAffineParamStep(double step);

    /// Get the affine parameter step (proper time step).
    double getAffineParamStep();

    /// Reset affine parameter stepsize.
    void resetAffineParamStep();

    /// Set maxmimum affine parameter stepsize.
    virtual void setMaxAffineParamStep(double step);

    /// Get maximum affine parameter stepsize.
    double getMaxAffineParamStep();

    /// Set minimum affine parameter stepsize.
    virtual void setMinAffineParamStep(double step);

    /// Get minimum affine parameter stepsize.
    double getMinAffineParamStep();

    void setE(unsigned int i, vec4 ee);
    vec4 getE(unsigned int i);
    void getE0(double* p);
    void getE1(double* p);
    void getE2(double* p);
    void getE3(double* p);
    void getTetrad(double* e0, double* e1, double* e2, double* e3);
    void getTetrad(vec4& e0, vec4& e1, vec4& e2, vec4& e3);
    void getTetrad(mat4& m);
    void getTetradInv(mat4& m);
    void getTetrad(float* m);
    void getTetradInv(float* m);

    vec4 coordToLocal(vec4 cv);

    /// Test whether tetrad is orthonormal.
    bool isOrthonormal();

    /**
     * @brief The bounding box represents the domain for integration.
     * @param p1 : pointer to bounding box.
     * @param p2 : pointer to bounding box.
     */
    void setBoundingBox(double p1[4], double p2[4]);
    void setBoundingBox(vec4 p1, vec4 p2);
    void getBoundingBox(vec4& p1, vec4& p2);

    /// Checks whether the current position 'y' is outside the bounding box.
    bool outsideBoundBox();

    /// Set epsilon for the constraint equation.
    void setConstrEps(double eps);

    /// Get epsilon for the constraint equation.
    double getConstrEps();

    /// Checks whether the local tetrad is right handed.
    bool isRightHanded();

    /**
     * @brief Do a Gram-Schmidt orthonormalization.
     *  Note that the base vectors will be normalized with respect to the euclidean norm.
     *  Thus, the orthonormalization can be used only for local tetrads that are given
     *  with respect to a natural local tetrad!
     */
    bool gramSchmidtOrth();

    /// Test constraint condition.
    virtual double testConstraint();

    void getCurrentArray(double* cy);

    //! Get duration of calculation.
    double getCalcTime() { return mCalcTime; }

    void printTetrad(FILE* fptr = stderr);

protected:
    //! Calculate the right side of the parallel transport.
    virtual bool calcDerivs(const double yn[], double dydx[]) = 0;

protected:
    //! Pointer to the actual metric.
    Metric* mMetric;

    //! Affine parameter.
    double mLambda;
    //! Affine parameter stepsize.
    double mLambdaStep;
    //! Initial affine parameter stepsize.
    double mLambdaStepInit;
    //! Maximum affine parameter stepsize.
    double mMaxLambdaStep;
    //! Minimum affine parameter stepsize.
    double mMinLambdaStep;

    //! This array holds the current position, direction, and all tetrad vectors.
    double y[DEF_MAX_YS];
    //! Epsilon for the constraint equation.
    double mConstraintEpsilon;

    //! Bounding box minimum.
    double mBoundBoxMin[4];
    //! Bounding box maximum.
    double mBoundBoxMax[4];

    //! Time in seconds for calculation of geodesic.
    double mCalcTime;
};

} // end namespace m4d

#endif
