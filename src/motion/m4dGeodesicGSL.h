// --------------------------------------------------------------------------------
/*
    m4dGeodesicGSL.h

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

/*!  \class  m4d::GeodesicGSL
     \brief  Calculate geodesics with the GNU Scientific Library

             GSL:

 */
// --------------------------------------------------------------------------------

#ifndef M4D_GEODESIC_GSL_H
#define M4D_GEODESIC_GSL_H

#include <iostream>
#include <vector>
#include <cassert>
#include <cmath>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>

#include "m4dGeodesic.h"

#define  DEF_GSL_LAMBDA_MAX  1.0e38

namespace m4d {

enum enum_gslint_type {
    enum_gslint_geodesic = 0,
    enum_gslint_geodesic_data,
    enum_gslint_partrans,
    enum_gslint_partrans_jacobi
};

// prototpyes
int func_adaptor_geod(double x, const double y[], double f[], void *params);
int jac_adaptor_geod(double x, const double y[], double *dfdy, double dfdt[], void *params);

int func_adaptor_par(double x, const double y[], double f[], void *params);

int func_adaptor_jacobi(double x, const double y[], double f[], void *params);


// ---------------------------------------------------
//    class definition:   GeodesicGSL
// ---------------------------------------------------
class API_EXPORT GeodesicGSL : public Geodesic {
public:
    GeodesicGSL(Metric* metric, const gsl_odeiv_step_type*  step_type, int solver_type,
                enum_geodesic_type  type = enum_geodesic_lightlike);
    virtual ~GeodesicGSL();

    // --------- public methods -----------
public:
    virtual enum_break_condition  initializeGeodesic(const vec4 initPos, const vec4 initDir, double &cstr);

    virtual enum_break_condition  calculateGeodesic(const vec4 initPos, const vec4 initDir, const int maxNumPoints,
            std::vector<vec4> &points, std::vector<vec4> &dirs, std::vector<double> &lambda);

    virtual enum_break_condition  calculateGeodesic(const vec4 initPos, const vec4 initDir, const int maxNumPoints,
            vec4 *&points, vec4 *&dirs, int &numPoints);

    virtual enum_break_condition  calculateGeodesicData(const vec4 initPos, const vec4 initDir, const int maxNumPoints,
            std::vector<vec4> &points, std::vector<vec4> &dirs, std::vector<double> &epsilons, std::vector<double> &lambda);

    virtual enum_break_condition  calcParTransport(const vec4 initPos, const vec4 initDir,
            const vec4 e0, const vec4 e1, const vec4 e2, const vec4 e3,
            const int maxNumPoints,
            std::vector<vec4> &points, std::vector<vec4> &dirs,
            std::vector<double> &lambda,
            std::vector<vec4> &base0, std::vector<vec4> &base1,
            std::vector<vec4> &base2, std::vector<vec4> &base3);

    virtual enum_break_condition  calcSachsJacobi(const vec4 initPos, const vec4 initCoordDir,
            const vec3 localNullDir, const vec3 locX, const vec3 locY, const vec3 locZ,
            const vec4 e0, const vec4 e1, const vec4 e2, const vec4 e3,
            const enum_nat_tetrad_type  tetrad_type,
            const int maxNumPoints,
            std::vector<vec4> &points, std::vector<vec4> &dirs,
            std::vector<double> &lambda,
            std::vector<vec4> &sachs0, std::vector<vec4> &sachs1,
            std::vector<vec5> &jacobi, vec5 &maxJacobi);

    virtual enum_break_condition  calcSachsJacobi(const vec4 initPos, const vec4 initCoordDir,
            const vec3 localNullDir, const vec3 locX, const vec3 locY, const vec3 locZ,
            const vec4 b0, const vec4 b1, const vec4 b2, const vec4 b3,
            const enum_nat_tetrad_type  tetrad_type,
            const int maxNumPoints,
            vec4 *&points, vec4 *&dirs,
            double *&lambda,
            vec4 *&sachs0, vec4 *&sachs1,
            vec5 *&jacobi, vec5 &maxJacobi, int &numPoints);

    int   func_geod(double x, const double y[], double f[], void *params);
    int   jac_geod(double x, const double y[], double *dfdy, double dfdt[], void *params);

    int   func_par(double x, const double y[], double f[], void *params);

    int   func_jacobi(double x, const double y[], double f[], void *params);

    virtual bool  nextStep(int &status);
    virtual bool  nextStepPar(int &status);
    virtual bool  nextStepSachsJacobi(int &status);

    virtual void  printF(FILE* fptr = stderr);


    // --------- protected methods -----------
protected:
    void initialize(const gsl_odeiv_step_type*  step_type,  enum_gslint_type type = enum_gslint_geodesic);
    void allocMemory();
    void freeMemory();

    // -------- protected attribute ---------
protected:
    const gsl_odeiv_step_type*  mStepType;
    int                         mSolverType;
    gsl_odeiv_step*      mStep;
    gsl_odeiv_control*   mControl;
    gsl_odeiv_evolve*    mEvolve;

    gsl_odeiv_system     mSys;

    double    dydt_in[DEF_MAX_YS], dydt_out[DEF_MAX_YS];
};

} // end namespace m4d

#endif

