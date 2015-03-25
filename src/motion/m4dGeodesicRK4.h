// --------------------------------------------------------------------------------
/*
    m4dGeodesicRK4.h

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

/*!  \class  m4d::GeodesicRK4
     \brief  Calculate geodesics with a standard Runge-Kutta fourth-order method.



 */
// --------------------------------------------------------------------------------

#ifndef M4D_GEODESIC_RK4_H
#define M4D_GEODESIC_RK4_H

#include <iostream>
#include <vector>
#include <cassert>
#include <cmath>

#include <m4dGlobalDefs.h>
#include <metric/m4dMetric.h>
#include <motion/m4dGeodesic.h>

namespace m4d {

// ---------------------------------------------------
//    class definition:   GeodesicRK4
// ---------------------------------------------------
class MOTION_API GeodesicRK4 : public Geodesic {
public:
    GeodesicRK4(Metric* metric, enum_geodesic_type  type = enum_geodesic_lightlike);
    virtual ~GeodesicRK4();

    // --------- public methods -----------
public:
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

    virtual bool    nextStep(int &status);
    virtual bool    nextStepPar(int &status);
    virtual bool    nextStepSachsJacobi(int &status);

    //   virtual enum_break_condition  nextStep              ( double &constraint );
    //   virtual enum_break_condition  nextStepPar           ( double &constraint );
    //   virtual enum_break_condition  nextStepSachsJacobi   ( double &constraint );

    virtual void                  print(FILE* fptr = stderr);

    // -------- protected attribute ---------
protected:
    bool       mCalcWithParTransport;
    int        mNumCoords;
};

} // end namespace m4d

#endif

