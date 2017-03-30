// --------------------------------------------------------------------------------
/*
    m4dMotion.h

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

/*!  \class  m4d::Motion
     \brief  Base class for motion/geodesics in 4d spacetimes.

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


*/
// --------------------------------------------------------------------------------

#ifndef M4D_MOTION_H
#define M4D_MOTION_H

#include <iostream>
#include <cassert>
#include <cmath>

#include <m4dGlobalDefs.h>
#include <metric/m4dMetric.h>
#include <extra/m4dUtilities.h>

#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>

#ifdef _WIN32
#ifndef __GNUC__
#pragma warning (disable: 4244 )
#endif
#endif

namespace m4d {

// ---------------------------------------------------
//    class definition:   Motion
// ---------------------------------------------------
class API_M4D_EXPORT Motion {
public:
    Motion(Metric* metric);
    virtual ~Motion();

// --------- public methods -----------
public:
    virtual bool  setInitialPosition(double ipos[4]);
    virtual bool  setInitialPosition(vec4 ipos);

    virtual bool  setInitialDirection(double d0, double d1, double d2, double d3);
    virtual bool  setInitialDirection(vec4 dir);

    virtual void  setInitialTetrad(double e0[4], double e1[4], double e2[4], double e3[4]);
    virtual void  setInitialTetrad(vec4 e0, vec4 e1, vec4 e2, vec4 e3);

    virtual bool  setLocalInitialDir(double d1, double d2, double d3, double beta = 1.0);
    virtual bool  setLocalInitialDir(vec3 dir, double beta = 1.0);

    virtual  void        setMetric(Metric* metric);
    virtual  Metric*     getMetric();

    virtual bool  setParam(std::string paramName, double val);
    virtual bool  setParam(std::string paramName, double v0, double v1, double v2, double v3);

    void    getPosition(double* p);
    vec4    getPosition();

    void    getDirection(double* p);
    vec4    getDirection();

    void    resetAffineParam(double pt = 0.0);
    double  getAffineParam();
    void    setAffineParamStep(double step);
    double  getAffineParamStep();
    void    resetAffineParamStep();

    virtual void    setMaxAffineParamStep(double step);
    double  getMaxAffineParamStep();
    virtual void    setMinAffineParamStep(double step);
    double  getMinAffineParamStep();

    void    setE(unsigned int i, vec4 ee);
    vec4    getE(unsigned int i);
    void    getE0(double* p);
    void    getE1(double* p);
    void    getE2(double* p);
    void    getE3(double* p);
    void    getTetrad(double* e0, double* e1, double* e2, double* e3);
    void    getTetrad(vec4 &e0, vec4 &e1, vec4 &e2, vec4 &e3);
    void    getTetrad(mat4 &m);
    void    getTetradInv(mat4 &m);
    void    getTetrad(float* m);
    void    getTetradInv(float* m);

    vec4    coordToLocal(vec4 cv);

    bool    isOrthonormal();

    void    setBoundingBox(double p1[4], double p2[4]);
    void    setBoundingBox(vec4 p1, vec4 p2);
    void    getBoundingBox(vec4 &p1, vec4 &p2);
    bool    outsideBoundBox();

    void    setConstrEps(double eps);
    double  getConstrEps();

    bool    isRightHanded();
    bool    gramSchmidtOrth();

    virtual double  testConstraint();

    void    getCurrentArray(double* cy);

    //! Get duration of calculation.
    double  getCalcTime() {
        return mCalcTime;
    }

    void    printTetrad(FILE* fptr = stderr);

// --------- protected methods -----------
protected:
    //! Calculate the right side of the parallel transport.
    virtual bool calcDerivs(const double yn[], double dydx[]) = 0;

// -------- protected attribute ---------
protected:
    //! Pointer to the actual metric.
    Metric*    mMetric;

    //! Affine parameter.
    double     mLambda;
    //! Affine parameter stepsize.
    double     mLambdaStep;
    //! Initial affine parameter stepsize.
    double     mLambdaStepInit;
    //! Maximum affine parameter stepsize.
    double     mMaxLambdaStep;
    //! Minimum affine parameter stepsize.
    double     mMinLambdaStep;

    //! This array holds the current position, direction, and all tetrad vectors.
    double     y[DEF_MAX_YS];
    //! Epsilon for the constraint equation.
    double     mConstraintEpsilon;

    //! Bounding box minimum.
    double     mBoundBoxMin[4];
    //! Bounding box maximum.
    double     mBoundBoxMax[4];

    //! Time in seconds for calculation of geodesic.
    double     mCalcTime;
};

} // end namespace m4d

#endif

