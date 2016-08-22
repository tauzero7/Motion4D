// --------------------------------------------------------------------------------
/*
    m4dGeodesicBS.h

  Copyright (c) 2010-2014  Thomas Mueller


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

/*!  \class  m4d::GeodesicBS
     \brief  Calculate geodesics with the Bulirsch-Stoer method.



 */
// --------------------------------------------------------------------------------

#ifndef M4D_GEODESIC_BS_H
#define M4D_GEODESIC_BS_H

#include <iostream>
#include <vector>
#include <cassert>
#include <cmath>

#include <m4dGlobalDefs.h>
#include <metric/m4dMetric.h>
#include <motion/m4dGeodesic.h>

#include "m4dGeodesic.h"

#define  DEF_BS_TINY         1e-30
#define  DEF_BS_MAX_ROW_NUM  8
#define  DEF_BS_MAX_ROW_NUMS (DEF_BS_MAX_ROW_NUM+1)
#define  DEF_BS_SAFE1        0.25
#define  DEF_BS_SAFE2        0.7
#define  DEF_BS_MAX_RED      1e-5
#define  DEF_BS_MIN_RED      0.7
#define  DEF_BS_MAX_SCALE    0.1

namespace m4d {

// ---------------------------------------------------
//    class definition:   GeodesicBS
// ---------------------------------------------------
class API_EXPORT GeodesicBS : public Geodesic {
public:
    GeodesicBS(Metric* metric, enum_geodesic_type  type = enum_geodesic_lightlike);
    virtual ~GeodesicBS();

    // --------- public methods -----------
public:
    virtual void setMaxAffineParamStep(double hmax);
    void    setNumberOfSteps(int numSteps);

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
            const vec4 b0, const vec4 b1, const vec4 b2, const vec4 b3,
            const enum_nat_tetrad_type  tetrad_type,
            const int maxNumPoints,
            std::vector<vec4> &points,  std::vector<vec4> &dirs,
            std::vector<double> &lambda,
            std::vector<vec4> &sachs0, std::vector<vec4> &sachs1,
            std::vector<vec5> &jacobi, vec5 &maxJacobi);

    virtual enum_break_condition  calcSachsJacobi(const vec4 initPos, const vec4 initCoordDir,
            const vec3 localNullDir, const vec3 locX, const vec3 locY, const vec3 locZ,
            const vec4 e0, const vec4 e1, const vec4 e2, const vec4 e3,
            const enum_nat_tetrad_type  tetrad_type,
            const int maxNumPoints,
            vec4 *&points, vec4 *&dirs,
            double *&lambda,
            vec4 *&sachs0, vec4 *&sachs1,
            vec5 *&jacobi, vec5 &maxJacobi, int &numPoints);


    virtual enum_break_condition  nextStep(double htry, double &hdid, double &hnext, double &constraint);
    virtual enum_break_condition  nextStepSachsJacobi(double htry, double &hdid, double &hnext, double &constraint);

    virtual bool    nextStep(int &status);
    virtual bool    nextStepPar(int &status);
    virtual bool    nextStepSachsJacobi(int &status);


    virtual void                  print(FILE* fptr = stderr);


    // --------- protected methods -----------
protected:
    void       modMidPoint(double* yy, double* dydx, double H, int numSteps, double* yout);
    void       modMidPointJS(double* yy, double* dydx, double H, int numSteps, double* yout);
    void       polyExtrpol(int iest, double xest, double* yest, double* yz, double* dy);

    // -------- protected attribute ---------
protected:
    double     ym[DEF_MAX_YS];
    double     yn[DEF_MAX_YS];
    double     dydx[DEF_MAX_YS];
    double     yerr[DEF_MAX_YS];
    double     ysav[DEF_MAX_YS];
    double     yseq[DEF_MAX_YS];
    double     yscal[DEF_MAX_YS];
    double     cc[DEF_MAX_YS];

    double**   dd;
    double     a[DEF_BS_MAX_ROW_NUMS];
    double     alf[DEF_BS_MAX_ROW_NUM][DEF_BS_MAX_ROW_NUM];
    double     err[DEF_BS_MAX_ROW_NUM];
    double     x[DEF_BS_MAX_ROW_NUM];
    double     eps, epsold, xnew;

    int        nSeq[DEF_BS_MAX_ROW_NUMS + 1];

    int        mFirst, mKmax, mKopt;
    int        mReduct, mExitFlag;

    long double     mhmin;
    double     mhmax;
    int        mNumCoords;
    int        mNumSteps;
};

} // end namespace m4d

#endif
