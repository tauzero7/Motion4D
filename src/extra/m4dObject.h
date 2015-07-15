// --------------------------------------------------------------------------------
/*
    m4dObject.h

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

/*!  \class  m4d::Object
     \brief  Master object that stores all relevant data.


     Settings file consist of:
   \verbatim
        METRIC            <string>
        PARAM             <num>  <name>  <double>
        INIT_POS          <double>  <double>  <double>  <double>
        INIT_DIR          <double>  <double>  <double>
        INIT_ANGLE_VEL    <double>  <double>  <double>
        TIME_DIR          <int>
        GEOD_SOLVER_TYPE  <int>
        GEODESIC_TYPE     <int>
        STEPSIZE_CTRL     <int>
        STEPSIZE          <double>
        STEPSIZE_MAX_MIN  <double> <double>
        EPSILONS          <double> <double>
        CONSTR_EPSILON    <double>
        MAX_NUM_POINTS    <int>
        TETRAD_TYPE       <int>
        BASE_0            <double>  <double>  <double>  <double>
        BASE_1            <double>  <double>  <double>  <double>
        BASE_2            <double>  <double>  <double>  <double>
        BASE_3            <double>  <double>  <double>  <double>
        BOOST             <double>  <double>  <double>
        SPEED_OF_LIGHT    <double>
        GRAV_CONSTANT     <double>
        DIELECTRIC_PERM   <double>                                 \endverbatim

    Description:
    <ul>
      <li>METRIC   :  name of the metric defined in each class.
      <li>PARAM    :  continuous number,  name of the parameter,  value.
      <li>INIT_POS :  initial position in coordinates  (x0, x1, x2, x3).
      <li>INIT_DIR :  initial direction in coordinates (n1, n2, n3) with respect to local tetrad.
      <li>INIT_ANGLE_VEL : initial angles and velocity  (ksi, chi, vel) with respect to local tetrad.
      <li>TIME_DIR : time direction (>0: future, <0: past).
      <li>GEOD_SOLVER_TYPE : type of geodesic solver as numerical value.
      <li>GEODESIC_TYPE : type of geodesic (-1: timelike, 0: lightlike, 1: spacelike).
      <li>STEPSIZE_CTRL : stepsize control (1: yes, 0: no).
      <li>STEPSIZE : stepsize.
      <li>STEPSIZE_MAX_MIN : maximum and minimum stepsize.
      <li>EPSILONS : absolute and relative epsilons for geodesic integration.
      <li>CONSTR_EPSILON: epsilon for the constraint condition.
      <li>RESIZE_EPSILON: epsilon when resizing of y[4] shall be done.
      <li>MAX_NUM_POINTS : maximum number of points to be calculated.
      <li>TETRAD_TYPE : type of local tetrad.
      <li>BASE_0 : base vector of local tetrad with respect to natural local tetrad at initial position.
      <li>BASE_1 : base vector of local tetrad with respect to natural local tetrad at initial position.
      <li>BASE_2 : base vector of local tetrad with respect to natural local tetrad at initial position.
      <li>BASE_3 : base vector of local tetrad with respect to natural local tetrad at initial position.
      <li>BOOST  : boost parameters (ksi, chi, vel).
      <li>SPEED_OF_LIGHT : value of the speed of light (geometric units: 1, physical units: 299792458.0)
      <li>GRAV_CONSTANT  : value of the gravitational constant (geometric units: 1, physical units: 6.67...)
      <li>DIELECTRIC_PERM : value of the dielectric permittivity (geometric units: 1, physical units: ...)
    </ul>
    \sa testDatabase.cpp
    \sa enum_nat_tetrad_type
*/
// --------------------------------------------------------------------------------

#ifndef M4D_OBJECT_H
#define M4D_OBJECT_H


#include <iostream>
#include <fstream>
#include <cassert>

#include <m4dGlobalDefs.h>
#include <metric/m4dMetricDatabase.h>
#include <motion/m4dMotionList.h>
#include <motion/m4dMotionDatabase.h>
#include <extra/m4dUtilities.h>

namespace m4d {

// ---------------------------------------------------
//    class definition:   Object
// ---------------------------------------------------
class EXTRA_API Object {
public:
    Object();
    ~Object();

    // --------- public methods -----------
public:
    bool   setMetric(std::string metricName);
    bool   setMetricParam(std::string paramName, double value);

    bool   setSolver(std::string solverName);
    bool   setSolverParam(std::string paramName, bool val);
    bool   setSolverParam(std::string paramName, double value);
    bool   setSolverParam(std::string paramName, double v0, double v1, double v2, double v3);

    bool   setInitialPosition(double x0, double x1, double x2, double x3);
    bool   setInitialDirection(double v0, double v1, double v2, double v3);

    bool   setInitialLocalNullDirection(enum_time_direction tdir,
                                        double l0, double l1, double l2,
                                        enum_nat_tetrad_type nattype = enum_nat_tetrad_default);

    bool   setInitialLocalTimeDirection(enum_time_direction tdir,
                                        double l0, double l1, double l2, double beta,
                                        enum_nat_tetrad_type natType = enum_nat_tetrad_default);

    enum_break_condition  calculateGeodesic(int numPoints);

    void   clearAll();
    void   resetAll();

    bool   getParam(std::string paramName, int &paramValue);
    bool   getParam(std::string paramName, double &paramValue);
    bool   getParam(std::string paramName, m4d::vec3 &paramValue);
    bool   getParam(std::string paramName, m4d::vec4 &paramValue);

    bool   setLorentzTransf(const double chi, const double ksi, const double beta);
    void   resetLorentzTransf();

    bool   loadSettings(std::string filename, bool printset = false);
    bool   saveSettings(std::string filename, std::string dat = std::string());
    void   printSettings(FILE* fptr = stderr);

    bool   makeReport(std::string  &text);

    // --------- public attributes --------
public:
    MetricDatabase*       metricDB;
    Metric*               currMetric;
    IntegratorDatabase*   solverDB;
    Geodesic*             geodSolver;
    enum_integrator       geodSolverType;

    enum_geodesic_type    type;
    bool                  stepsizeControlled;
    double                stepsize;
    double                max_stepsize;
    double                min_stepsize;
    double                epsAbs;
    double                epsRel;
    double                epsConstr;
    double                epsResize;

    vec4                  startPos;
    vec3                  startDir;
    double                ksi;
    double                chi;
    double                vel;
    vec4                  coordDir;
    vec4                  base[4];
    bool                  isBaseInCoords;

    int                   axes_orient;
    double                boost_ksi;
    double                boost_chi;
    double                boost_beta;
    mat4                  lorentz;

    int                   timeDirection;
    enum_nat_tetrad_type  tetradType;
    unsigned int          maxNumPoints;
    std::vector<vec4>     points;
    std::vector<vec4>     dirs;
    std::vector<double>   lambda;
    std::vector<vec4>     sachs1;
    std::vector<vec4>     sachs2;
    std::vector<vec5>     jacobi;
    vec5                  maxJacobi;

    double                speed_of_light;
    double                grav_constant;
    double                dielectric_perm;
};

} // end namespace m4d

#endif
