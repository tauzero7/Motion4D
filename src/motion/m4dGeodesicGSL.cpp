// -------------------------------------------------------------------------------
/*
    m4dGeodesicGSL.cpp

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

#include "m4dGeodesicGSL.h"
#include "m4dMotionList.h"

namespace m4d {

//----------------------------------------------------------------------------
//         func-,jac-adaptor
//----------------------------------------------------------------------------
int func_adaptor_geod(double x, const double y[], double f[], void *params) {
    GeodesicGSL * obj = reinterpret_cast<GeodesicGSL *>(params);
    return obj->func_geod(x, y, f, NULL);
}

int jac_adaptor_geod(double x, const double y[], double *dfdy, double dfdt[], void *params) {
    GeodesicGSL * obj = reinterpret_cast<GeodesicGSL *>(params);
    return obj->jac_geod(x, y, dfdy, dfdt, NULL);
}

int func_adaptor_par(double x, const double y[], double f[], void *params) {
    GeodesicGSL * obj = reinterpret_cast<GeodesicGSL *>(params);
    return obj->func_par(x, y, f, NULL);
}

int func_adaptor_jacobi(double x, const double y[], double f[], void *params) {
    GeodesicGSL * obj = reinterpret_cast<GeodesicGSL *>(params);
    return obj->func_jacobi(x, y, f, NULL);
}


/*! Standard constructor for geodesic motion.
 *
 *  \param metric : pointer to metric.
 *  \param step_type   : integrator step type.
 *  \param solver_type : solver type number.
 *  \param type : type of geodesic.
 *  \sa enum_geodesic_type.
 */
GeodesicGSL::GeodesicGSL(Metric* metric, const gsl_odeiv_step_type*  step_type, int solver_type,
                           enum_geodesic_type  type)
    : Geodesic(metric, type) {
    mStepsizeControlled = false;
    epsilon_abs = 1.0e-6;
    epsilon_rel = 0.0;

    mStep = NULL;
    mControl = NULL;
    mEvolve  = NULL;

    mStepType = step_type;
    mSolverType = solver_type;

    gsl_set_error_handler_off();
}


GeodesicGSL::~GeodesicGSL() {
    freeMemory();
}

// *********************************** public methods ******************************

/*! Initialize geodesic.
 * \param initPos : initial position.
 * \param initDir : initial coordinate direction.
 * \param cstr
 */
enum_break_condition GeodesicGSL::initializeGeodesic(const vec4 initPos, const vec4 initDir, double &cstr) {
    resetAffineParam();
    resetAffineParamStep();

    setInitialPosition(initPos);
    setInitialDirection(initDir);

    if ((cstr = fabs(testConstraint())) > mConstraintEpsilon) {
        return enum_break_constraint;
    }

    if (mCalcWithParTransport) {
        initialize(mStepType, enum_gslint_partrans);
    } else {
        initialize(mStepType);
    }
    allocMemory();

    // initialize dydt_in from system parameters
    GSL_ODEIV_FN_EVAL(&mSys, mLambda, y, dydt_in);
    return enum_break_none;
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
enum_break_condition GeodesicGSL::calculateGeodesic(const vec4 initPos, const vec4 initDir,
    const int maxNumPoints,
    std::vector<vec4> &points , std::vector<vec4> &dirs, std::vector<double> &lambda)
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

    enum_break_condition  breakType = enum_break_none;
    double cstr;
    breakType = initializeGeodesic(initPos, initDir, cstr);
    if (breakType == enum_break_constraint) {
        return enum_break_constraint;
    }

    points.push_back(vec4(y[0], y[1], y[2], y[3]));
    dirs.push_back(vec4(y[4], y[5], y[6], y[7]));
    lambda.push_back(mLambda);

    int64_t t1 = get_system_clock();

    double tc;
    int status;
    while ((int)points.size() < maxNumPoints && breakType == enum_break_none) {
        if (!nextStep(status)) {
            breakType = enum_break_cond;
            if (status == GSL_EBADFUNC) {
                breakType = enum_break_not_implemented;
            }
        } else {
            if (!outsideBoundBox()) {
                vec4 p = vec4(y[0], y[1], y[2], y[3]);
                points.push_back(p);
                dirs.push_back(vec4(y[4], y[5], y[6], y[7]));
                lambda.push_back(mLambda);
                if ((tc = fabs(mMetric->testConstraint(y, mKappa))) > mConstraintEpsilon) {
                    breakType = enum_break_constraint;
                }

                if (tc > resizeEps) {
                    mMetric->resize(y, mKappa, resizeFac);
                    GSL_ODEIV_FN_EVAL(&mSys, mLambda, y, dydt_in);
                }
            } else {
                breakType = enum_break_outside;
            }
        }
    }

    if ((int)points.size() >= maxNumPoints) {
        breakType = enum_break_num_exceed;
    }

    int64_t t2 = get_system_clock();
    mCalcTime = (t2 - t1) * 1e-6;
    freeMemory();
    return breakType;
}


enum_break_condition GeodesicGSL ::calculateGeodesic(const vec4 initPos, const vec4 initDir,
    const int maxNumPoints,
    vec4 *&points, vec4 *&dirs, int &numPoints)
{
    if (points != NULL) {
        delete [] points;
    }
    if (dirs != NULL) {
        delete [] dirs;
    }

    points = new vec4[maxNumPoints];
    dirs   = new vec4[maxNumPoints];

    enum_break_condition  breakType = enum_break_none;
    double cstr;
    breakType = initializeGeodesic(initPos, initDir, cstr);
    if (breakType == enum_break_constraint) {
        return enum_break_constraint;
    }

    points[0] = vec4(y[0], y[1], y[2], y[3]);
    dirs[0]   = vec4(y[4], y[5], y[6], y[7]);

    int64_t t1 = get_system_clock();

    double tc;
    int status;

    numPoints = 1;
    while (numPoints < maxNumPoints && breakType == enum_break_none) {
        if (!nextStep(status)) {
            breakType = enum_break_cond;
            if (status == GSL_EBADFUNC) {
                breakType = enum_break_not_implemented;
            }
        } else {
            if (!outsideBoundBox()) {
                points[numPoints] = vec4(y[0], y[1], y[2], y[3]);
                dirs[numPoints]   = vec4(y[4], y[5], y[6], y[7]);
                numPoints++;
                if ((tc = fabs(mMetric->testConstraint(y, mKappa))) > mConstraintEpsilon) {
                    breakType = enum_break_constraint;
                }

                if (tc > resizeEps) {
                    mMetric->resize(y, mKappa, resizeFac);
                    GSL_ODEIV_FN_EVAL(&mSys, mLambda, y, dydt_in);
                }
            } else {
                breakType = enum_break_outside;
            }
        }
    }

    if (numPoints == maxNumPoints) {
        breakType = enum_break_num_exceed;
    }

    int64_t t2 = get_system_clock();
    mCalcTime = (t2 - t1) * 1e-6;
    freeMemory();
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
GeodesicGSL::calculateGeodesicData(const vec4 initPos, const vec4 initDir, const int maxNumPoints,
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

    setInitialPosition(initPos);
    setInitialDirection(initDir);

    enum_break_condition  breakType = enum_break_none;

    double cstr;
    if ((cstr = fabs(testConstraint())) > mConstraintEpsilon) {
        return enum_break_constraint;
    }

    vec4 p(y, 4);
    vec4 d(&(y[4]), 4);
    points.push_back(p);
    dirs.push_back(d);
    epsilons.push_back(cstr);

    initialize(mStepType);
    allocMemory();

    // initialize dydt_in from system parameters
    resetAffineParam();
    resetAffineParamStep();
    GSL_ODEIV_FN_EVAL(&mSys, mLambda, y, dydt_in);

    int64_t t1 = get_system_clock();

    int status;
    while ((int)points.size() < maxNumPoints && breakType == enum_break_none) {
        if (!nextStep(status)) {
            breakType = enum_break_cond;
            if (status == GSL_EBADFUNC) {
                breakType = enum_break_not_implemented;
            }
        } else {
            if (!outsideBoundBox()) {
                if ((cstr = fabs(mMetric->testConstraint(y, mKappa))) > mConstraintEpsilon) {
                    breakType = enum_break_constraint;
                }

                if (cstr > resizeEps) {
                    mMetric->resize(y, mKappa, resizeFac);
                    GSL_ODEIV_FN_EVAL(&mSys, mLambda, y, dydt_in);
                }

                vec4 p(y, 4);
                vec4 d(&(y[4]), 4);
                points.push_back(p);
                dirs.push_back(d);
                epsilons.push_back(cstr);
                lambda.push_back(mLambda);
            } else {
                breakType = enum_break_outside;
            }
        }
    }
    if ((int)points.size() >= maxNumPoints) {
        breakType = enum_break_num_exceed;
    }

    int64_t t2 = get_system_clock();
    mCalcTime = (t2 - t1) * 1e-6;
    freeMemory();
    return breakType;
}

/*! Calculate a geodesic and the parallel transported local tetrad of the observer.
 *
 */
enum_break_condition
GeodesicGSL::calcParTransport(const vec4 initPos, const vec4 initDir,
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

    enum_break_condition  breakType = enum_break_none;

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

    initialize(mStepType, enum_gslint_partrans);
    allocMemory();

    // initialize dydt_in from system parameters
    resetAffineParam();
    resetAffineParamStep();
    GSL_ODEIV_FN_EVAL(&mSys, mLambda, y, dydt_in);

    int64_t t1 = get_system_clock();

    int status;
    while ((int)points.size() < maxNumPoints && breakType == enum_break_none) {
        if (!nextStepPar(status)) {
            breakType = enum_break_cond;
            if (status == GSL_EBADFUNC) {
                breakType = enum_break_not_implemented;
            }
        } else {
            if (!outsideBoundBox()) {
                if ((cstr = fabs(mMetric->testConstraint(y, mKappa))) > mConstraintEpsilon) {
                    breakType = enum_break_constraint;
                }
                if (cstr > resizeEps) {
                    mMetric->resize(y, mKappa, resizeFac);
                    GSL_ODEIV_FN_EVAL(&mSys, mLambda, y, dydt_in);
                }

                points.push_back(vec4(y[0], y[1], y[2], y[3]));
                dirs.push_back(vec4(y[4], y[5], y[6], y[7]));
                base0.push_back(vec4(y[8], y[9], y[10], y[11]));
                base1.push_back(vec4(y[12], y[13], y[14], y[15]));
                base2.push_back(vec4(y[16], y[17], y[18], y[19]));
                base3.push_back(vec4(y[20], y[21], y[22], y[23]));
                lambda.push_back(mLambda);
            } else {
                breakType = enum_break_outside;
            }
        }
    }
    if ((int)points.size() >= maxNumPoints) {
        breakType = enum_break_num_exceed;
    }

    int64_t t2 = get_system_clock();
    mCalcTime = (t2 - t1) * 1e-6;
    freeMemory();

    return  breakType;
}


/*! Calculate a geodesic and the parallel transported local tetrad of the observer.
 *
 *  \param  initPos      :  initial position of the geodesic in coordinates.
 *  \param  initCoordDir :  initial coordinate direction of the geodesic in coordinates.
 *  \param  localNullDir :  initial local direction of the geodesic in coordinates.
 *  \param  locX
 *  \param  locY
 *  \param  locZ
 *  \param  b0           :  b0 base vectors of local tetrad.
 *  \param  b1           :  b1 base vectors of local tetrad.
 *  \param  b2           :  b2 base vectors of local tetrad.
 *  \param  b3           :  b3 base vectors of local tetrad.
 *  \param  tetrad_type  :  type of local tetrad.
 *  \param  maxNumPoints :  maximum number of points.
 *  \param  points       :  reference to the calculated points.
 *  \param  dirs         :  reference to the calculated directions.
 *  \param  lambda       :  reference to the affine parameters.
 *  \param  sachs0       :  references to Sachs basis 0.
 *  \param  sachs1       :  references to Sachs basis 1.
 *  \param  jacobi       :  references to Jacobi parameters.
 *  \param  maxJacobi    :  references to maximum Jacobi parameter.
 *
 *  \return enum_break_condition : break condition.
 *  \sa enum_break_condition
 */
enum_break_condition
GeodesicGSL::calcSachsJacobi(const vec4 initPos, const vec4 initCoordDir,
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

    initialize(mStepType, enum_gslint_partrans_jacobi);
    allocMemory();

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

#if 0
    double prod;
    mMetric->calcProduct(initPos,initCoordDir,initCoordDir,prod);  fprintf(stderr,"%g\n",prod);
    mMetric->calcProduct(initPos,initCoordDir,bb1,prod);  fprintf(stderr,"p: %g\n",prod);
    mMetric->calcProduct(initPos,initCoordDir,bb2,prod);  fprintf(stderr,"p: %g\n",prod);  std::cerr << std::endl;
#endif

    enum_break_condition  breakType = enum_break_none;
    double cstr;
    if (fabs(testConstraint()) > mConstraintEpsilon) {
        return enum_break_constraint;
    }

    // initialize dydt_in from system parameters
    resetAffineParam();
    resetAffineParamStep();
    GSL_ODEIV_FN_EVAL(&mSys, mLambda, y, dydt_in);


    vec5 currJacobi;

    points.push_back(vec4(&y[0]));
    lambda.push_back(mLambda);
    dirs.push_back(vec4(&y[DEF_TG_IDX]));
    sachs0.push_back(vec4(&y[DEF_SA1_IDX]));
    sachs1.push_back(vec4(&y[DEF_SA2_IDX]));
    jacobi.push_back(vec5(0.0, 0.0, 1.0, 0.0, 0.0));
    maxJacobi = jacobi[0];



    int64_t t1 = get_system_clock();

    int status;
    while ((int)points.size() < maxNumPoints && breakType == enum_break_none) {
        if (!nextStepSachsJacobi(status)) {
            breakType = enum_break_cond;
            if (status == GSL_EBADFUNC) {
                breakType = enum_break_not_implemented;
            }
        } else {
            if (!outsideBoundBox()) {
                vec4 p = vec4(y[0], y[1], y[2], y[3]);
                points.push_back(p);
                dirs.push_back(vec4(&y[DEF_TG_IDX]));
                sachs0.push_back(vec4(&y[DEF_SA1_IDX]));
                sachs1.push_back(vec4(&y[DEF_SA2_IDX]));

#if 0
                // Sachs basis keeps perpendicular to the null geodesic:
                double prod;
                mMetric->calcProduct(p,vec4(&y[DEF_TG_IDX]),vec4(&y[DEF_TG_IDX]),prod);  fprintf(stderr,"%g\n",prod);
                mMetric->calcProduct(p,vec4(&y[DEF_TG_IDX]),vec4(&y[DEF_SA1_IDX]),prod);  fprintf(stderr,"p: %g\n",prod);
                mMetric->calcProduct(p,vec4(&y[DEF_TG_IDX]),vec4(&y[DEF_SA2_IDX]),prod);  fprintf(stderr,"p: %g\n",prod);  std::cerr << std::endl;
#endif
                /*
                fprintf(stdout,"%12.6f %12.6e %12.6e   %12.6e %12.6e %12.6e %12.6e    %12.6e %12.6e %12.6e %12.6e\n",
                y[0],y[1],y[2],
                y[DEF_JAC1_IDX],y[DEF_JAC1_IDX+1],y[DEF_JAC1_IDX+2],y[DEF_JAC1_IDX+3],
                y[DEF_DJ1_IDX],y[DEF_DJ1_IDX+1],y[DEF_DJ1_IDX+2],y[DEF_DJ1_IDX+3]);


                fprintf(stdout,"%12.6e %12.6e %12.6e   %12.6e %12.6e %12.6e %12.6e  \n",
                y[0],y[1],y[2],
                y[DEF_SA1_IDX],y[DEF_SA1_IDX+1],y[DEF_SA1_IDX+2],y[DEF_SA1_IDX+3]);
                */
                calcJacobiParams(mLambda, y, currJacobi);
                jacobi.push_back(currJacobi);
                lambda.push_back(mLambda);
                findMaxJacobi(currJacobi, maxJacobi);

                if ((cstr = fabs(mMetric->testConstraint(y, mKappa))) > mConstraintEpsilon) {
                    breakType = enum_break_constraint;
                }

                if (cstr > resizeEps) {
                    mMetric->resize(y, mKappa, resizeFac);
                    GSL_ODEIV_FN_EVAL(&mSys, mLambda, y, dydt_in);
                }
            } else {
                breakType = enum_break_outside;
            }
        }
    }

    if ((int)points.size() >= maxNumPoints) {
        breakType = enum_break_num_exceed;
    }

    int64_t t2 = get_system_clock();
    mCalcTime = (t2 - t1) * 1e-6;
    freeMemory();
    return breakType;
}


enum_break_condition
GeodesicGSL::calcSachsJacobi(const vec4 initPos, const vec4 initCoordDir,
                             const vec3 localNullDir, const vec3 locX, const vec3 locY, const vec3 locZ,
                             const vec4 b0, const vec4 b1, const vec4 b2, const vec4 b3,
                             const enum_nat_tetrad_type  tetrad_type,
                             const int maxNumPoints,
                             vec4 *&points, vec4 *&dirs,
                             double *&lambda,
                             vec4 *&sachs0, vec4 *&sachs1,
                             vec5 *&jacobi, vec5 &maxJacobi, int &numPoints) {
    setCalcWithParTransport(true);

    if (points != NULL) {
        delete [] points;
    }
    if (dirs != NULL) {
        delete [] dirs;
    }
    if (lambda != NULL) {
        delete [] lambda;
    }
    if (sachs0 != NULL) {
        delete [] sachs0;
    }
    if (sachs1 != NULL) {
        delete [] sachs1;
    }
    if (jacobi != NULL) {
        delete [] jacobi;
    }

    points = new vec4[maxNumPoints];
    dirs   = new vec4[maxNumPoints];
    lambda = new double[maxNumPoints];
    sachs0 = new vec4[maxNumPoints];
    sachs1 = new vec4[maxNumPoints];
    jacobi = new vec5[maxNumPoints];

    maxJacobi = vec5(-DBL_MAX, -DBL_MAX, -DBL_MAX, -DBL_MAX, -DBL_MAX);

    initialize(mStepType, enum_gslint_partrans_jacobi);
    allocMemory();

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
    if (fabs(testConstraint()) > mConstraintEpsilon) {
        return enum_break_constraint;
    }

    // initialize dydt_in from system parameters
    resetAffineParam();
    resetAffineParamStep();
    GSL_ODEIV_FN_EVAL(&mSys, mLambda, y, dydt_in);


    vec5 currJacobi;

    points[0] = vec4(&y[0]);
    dirs[0]   = vec4(&y[DEF_TG_IDX]);
    sachs0[0] = vec4(&y[DEF_SA1_IDX]);
    sachs1[0] = vec4(&y[DEF_SA2_IDX]);
    jacobi[0] = vec5(0.0, 0.0, 1.0, 0.0, 0.0);
    lambda[0] = mLambda;
    maxJacobi = jacobi[0];


    int64_t t1 = get_system_clock();

    numPoints = 1;

    int status;
    while (numPoints < maxNumPoints && breakType == enum_break_none) {
        if (!nextStepSachsJacobi(status)) {
            breakType = enum_break_cond;
            if (status == GSL_EBADFUNC) {
                breakType = enum_break_not_implemented;
            }
        } else {
            if (!outsideBoundBox()) {

                points[numPoints] = vec4(&y[0]);
                dirs[numPoints]   = vec4(&y[DEF_TG_IDX]);
                sachs0[numPoints] = vec4(&y[DEF_SA1_IDX]);
                sachs1[numPoints] = vec4(&y[DEF_SA2_IDX]);
                lambda[numPoints] = mLambda;

                calcJacobiParams(mLambda, y, currJacobi);
                jacobi[numPoints] = currJacobi;
                numPoints++;

                findMaxJacobi(currJacobi, maxJacobi);


                if ((cstr = fabs(mMetric->testConstraint(y, mKappa))) > mConstraintEpsilon) {
                    breakType = enum_break_constraint;
                }

                if (cstr > resizeEps) {
                    mMetric->resize(y, mKappa, resizeFac);
                    GSL_ODEIV_FN_EVAL(&mSys, mLambda, y, dydt_in);
                }
            } else {
                breakType = enum_break_outside;
            }
        }
    }

    if (numPoints == maxNumPoints) {
        breakType = enum_break_num_exceed;
    }

    int64_t t2 = get_system_clock();
    mCalcTime = (t2 - t1) * 1e-6;
    freeMemory();
    return breakType;
}


/*! Right-hand-side of the geodesic equation.
 *  \param x
 *  \param y[]
 *  \param f[]
 *  \param params
 */
int
GeodesicGSL::func_geod(double , const double y[], double f[], void *) {
    if (mMetric->calcDerivs(y, f)) {
        return GSL_SUCCESS;
    }

    f[0] = y[4];
    f[1] = y[5];
    f[2] = y[6];
    f[3] = y[7];

    mMetric->calculateChristoffels(y);

    for (int k = 0; k < 4; k++) {
        f[k + 4] = 0.0;
        for (int i = 0; i < 4; i++)
            for (int j = 0; j < 4; j++) {
                f[k + 4] += - mMetric->getChristoffel(i, j, k) * y[i + 4] * y[j + 4];
            }
    }
    return GSL_SUCCESS;
}

/*! Right-hand-side of the parallel transport.
 *  \param x
 *  \param y[]
 *  \param f[]
 *  \param params
 */
int
GeodesicGSL::func_par(double , const double y[], double f[], void*) {
    if (mMetric->calcDerivsPar(y, f)) {
        return GSL_SUCCESS;
    }

    if (calcDerivsPar(y, f)) {
        return GSL_SUCCESS;
    }
    return 0;
}


/*!
 *  \param x
 *  \param y[]
 *  \param dfdy
 *  \param dfdt
 *  \param params
 */
int
GeodesicGSL::jac_geod(double , const double y[], double *dfdy, double dfdt[], void *) {
    int dimension = 8;

    if (!mMetric->calculateChrisD(y)) {
        int num = dimension * dimension;
        for (int i = 0; i < num; i++) {
            dfdy[i] = 0.0;
        }
        for (int i = 0; i < dimension; i++) {
            dfdt[i] = 0.0;
        }

        fprintf(stderr, "Spacetime do not have a Jacobi matrix!\n");
        return GSL_EBADFUNC;  // gsl_errno.h: "problem with user-supplied function"
    } else {
        gsl_matrix_view   dfdy_mat = gsl_matrix_view_array(dfdy, dimension, dimension);
        gsl_matrix              *m = &dfdy_mat.matrix;

        mMetric->calculateChrisD(y);

        // TODO: incorporarte null/id-matrix from subsets
        double diff1, diff2;
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                // null-matrix
                gsl_matrix_set(m, i, j, 0.0);

                // id-matrix
                if (j + 4 == i) {
                    gsl_matrix_set(m, i, j + 4, 1.0);
                }

                // diff of christoffels
                diff1 = 0.0;
                diff2 = 0.0;
                for (int s = 0; s < 4; s++) {
                    for (int t = 0; t < 4; t++) {
                        diff1 += -mMetric->getChrisD(s, t, i, j) * y[s + 4] * y[t + 4];
                    }
                    diff2 += -(mMetric->getChristoffel(s, j, i)) * y[s + 4];
                }
                diff2 *= 2.0;

                gsl_matrix_set(m, i + 4, j, diff1);
                gsl_matrix_set(m, i + 4, i + 4, diff2);

                dfdy[i * dimension + j] = gsl_matrix_get(m, i, j);
            }
        }
    }
    return GSL_SUCCESS;
}

int
GeodesicGSL::func_jacobi(double , const double y[], double dydx[], void*) {
    if (!mMetric->calcDerivsSachsJacobi(y, dydx)) {
        mMetric->calculateChristoffels(y);
        if (!mMetric->calculateChrisD(y)) {
            return GSL_EBADFUNC;    // gsl_errno.h: "problem with user-supplied function"
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
    return GSL_SUCCESS;
}


/*! Print geodesic solver properties.
 * \param fptr : file pointer.
 */
void GeodesicGSL::printF(FILE* fptr) {
    fprintf(fptr, "\nGeodesicGSL:\n------------\n");
    fprintf(fptr, "\tstep type           : %s\n", stl_solver_names[mSolverType]);
    fprintf(fptr, "\tstepsize controlled : %s\n", ((mStepsizeControlled) ? "yes" : "no"));
    fprintf(fptr, "\tstep size           : %12.8e\n", mLambdaStep);
    fprintf(fptr, "\tepsilon_abs         : %12.8e\n", epsilon_abs);
    fprintf(fptr, "\tepsilon_rel         : %12.8e\n", epsilon_rel);
    fprintf(fptr, "\tconstraint epsilon  : %12.8e\n", mConstraintEpsilon);
    fprintf(fptr, "\tbounding box min    : %14.6e %14.6e %14.6e %14.6e\n", mBoundBoxMin[0], mBoundBoxMin[1], mBoundBoxMin[2], mBoundBoxMin[3]);
    fprintf(fptr, "\tbounding box max    : %14.6e %14.6e %14.6e %14.6e\n", mBoundBoxMax[0], mBoundBoxMax[1], mBoundBoxMax[2], mBoundBoxMax[3]);
    fprintf(fptr, "\tgeodesic type       : %s\n", stl_geodesic_type[mType]);
}

// ********************************* protected methods *****************************

/*! Initialize function and jacobi adaptor for the GSL integration.
 *
 *  \param  step_type :  constant pointer to a GSL step type.
 *  \param  type : type of gsl integration.
 */
void
GeodesicGSL::initialize(const gsl_odeiv_step_type*  step_type,  enum_gslint_type  type) {
    mStepType = step_type;

    switch (type) {
    default:
    case enum_gslint_geodesic:
        mSys.function = func_adaptor_geod;
        mSys.jacobian = jac_adaptor_geod;
        mSys.dimension = 8;
        break;
    case enum_gslint_geodesic_data:
        break;
    case enum_gslint_partrans:
        mSys.function = func_adaptor_par;
        mSys.dimension = DEF_MAX_YS_PAR;
        break;
    case enum_gslint_partrans_jacobi:
        mSys.function  = func_adaptor_jacobi;
        mSys.dimension = DEF_MAX_YS_JAC;
        break;
    }

    mSys.params = this;
}

/*! Allocate memory for the GSL integration.
 */
void
GeodesicGSL::allocMemory() {
    mStep = gsl_odeiv_step_alloc(mStepType, mSys.dimension);

    if (mStepsizeControlled) {
        mControl = gsl_odeiv_control_y_new(epsilon_abs, epsilon_rel);
        mEvolve  = gsl_odeiv_evolve_alloc(mSys.dimension);
    }
}

/*! Free memory after GSL integration.
 */
void
GeodesicGSL::freeMemory() {
    if (mStepsizeControlled) {
        if (mEvolve != NULL) {
            gsl_odeiv_evolve_free(mEvolve);
        }
        mEvolve = NULL;
        if (mControl != NULL) {
            gsl_odeiv_control_free(mControl);
        }
        mControl = NULL;
    }

    if (mStep != NULL) {
        gsl_odeiv_step_free(mStep);
    }
    mStep = NULL;
}

/*! Calculate next step of the geodesic.
 *
 *  \param status : gsl status after this step.
 */
bool GeodesicGSL::nextStep(int &status) {
    if (mMetric->breakCondition(&y[0])) {
        return false;
    }
    register int i;
    double y_err[DEF_MAX_YS];
    double lambda1 = DEF_GSL_LAMBDA_MAX;

    if (mStepsizeControlled) {
        status = gsl_odeiv_evolve_apply(mEvolve, mControl, mStep, &mSys, &mLambda, lambda1, &mLambdaStep, y);
        //std::cerr << "lambda = " << mLambda << " " << mLambdaStep << " " << status <<  std::endl;

        if (mLambdaStep > mMaxLambdaStep) {
            mLambdaStep = mMaxLambdaStep;
        }
    } else {
        status = gsl_odeiv_step_apply(mStep, mLambda, mLambdaStep, y, y_err, dydt_in, dydt_out, &mSys);

        for (i = 0; i < (int)mSys.dimension; i++) {
            dydt_in[i] = dydt_out[i];
        }
        mLambda += mLambdaStep;
    }

#if 0
    // experimental for TeoSimpleWH
    double l = y[1];
    double dt = y[4];
    double dphi = y[7];
    double b0;
    mMetric->getParam("b0",b0);
    double w = sqrt(l*l+b0*b0);
    double c1 = 2.0*dt + b0*b0/w*(dphi - b0*b0*0.5*pow(w,-3.0)*dt);
    double c2 = w*w*(dphi - b0*b0*0.5*pow(w,-3.0)*dt);
    double dl2 = -1.0 + pow(0.5*c1 - b0*b0*c2*0.5*pow(w,-3.0), 2.0) - c2*c2/(w*w);
    double dp = c2/(w*w) + b0*b0*c1*0.25*pow(w,-3.0) - pow(b0,4.0)*c2*0.25*pow(w,-6);
    fprintf(stderr,"c: %f %f %f %f\n",c1,c2,dphi,dp);
#endif
    if (status == GSL_EBADFUNC) {
        return false;
    }
    return true;
}

bool
GeodesicGSL::nextStepPar(int &status) {
    if (mMetric->breakCondition(&y[0])) {
        return false;
    }
    register int i;
    double y_err[DEF_MAX_YS_PAR];
    double lambda1 = DEF_GSL_LAMBDA_MAX;

    if (mStepsizeControlled) {
        status = gsl_odeiv_evolve_apply(mEvolve, mControl, mStep, &mSys, &mLambda, lambda1, &mLambdaStep, y);
        //std::cerr << "lambda = " << mLambda << " " << mLambdaStep << " " << status <<  std::endl;
    } else {
        status = gsl_odeiv_step_apply(mStep, mLambda, mLambdaStep, y, y_err, dydt_in, dydt_out, &mSys);

        for (i = 0; i < (int)mSys.dimension; i++) {
            dydt_in[i] = dydt_out[i];
        }

        mLambda += mLambdaStep;
    }

    if (status == GSL_EBADFUNC) {
        return false;
    }
    return true;
}

/*! Calculate next step of the parallel and Jacobi equation step.
 *
 *  \param status : gsl status after this step.
 */
bool GeodesicGSL::nextStepSachsJacobi(int &status) {
    if (mMetric->breakCondition(&y[0])) {
        return false;
    }
    register int i;
    double  y_err[DEF_MAX_YS_JAC];
    double  lambda1 = DEF_GSL_LAMBDA_MAX;

    if (mStepsizeControlled) {
        status = gsl_odeiv_evolve_apply(mEvolve, mControl, mStep, &mSys, &mLambda, lambda1, &mLambdaStep, y);
        //fprintf(stderr,"error: %s\n", gsl_strerror (status));
        // cerr << "lambda = " << mLambda << " " << mLambdaStep << endl;
    } else {
        status = gsl_odeiv_step_apply(mStep, mLambda, mLambdaStep, y, y_err, dydt_in, dydt_out, &mSys);

        for (i = 0; i < DEF_MAX_YS_JAC; i++) {
            dydt_in[i] = dydt_out[i];
        }

        mLambda += mLambdaStep;
    }

    if (status == GSL_EBADFUNC) {
        return false;
    }
    return true;
}

} // end namespace m4d

