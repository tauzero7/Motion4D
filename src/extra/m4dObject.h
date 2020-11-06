/**
 * @file    m4dObject.h
 * @author  Thomas Mueller
 *
 * @brief  Master object that stores all relevant data.

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

#ifndef M4D_OBJECT_H
#define M4D_OBJECT_H

#include <cassert>
#include <fstream>
#include <iostream>

#include <extra/m4dUtilities.h>
#include <m4dGlobalDefs.h>
#include <metric/m4dMetricDatabase.h>
#include <motion/m4dMotionDatabase.h>
#include <motion/m4dMotionList.h>

#ifdef _WIN32
#ifndef __GNUC__
#pragma warning(disable : 4244)
#endif
#endif

namespace m4d {

/**
 * @brief The Object class
 */
class API_M4D_EXPORT Object
{
public:
    Object();
    ~Object();

    bool setMetric(const char* metricName);
    bool setMetricParam(const char* paramName, double value);

    bool setSolver(const char* solverName);
    bool setSolverParam(const char* paramName, bool val);
    bool setSolverParam(const char* paramName, double value);
    bool setSolverParam(const char* paramName, double v0, double v1, double v2, double v3);

    bool setParam(const char* paramName, int val);

    bool setInitialPosition(const double* x);
    bool setInitialPosition(double x0, double x1, double x2, double x3);

    bool setInitialDirection(const double* v);
    bool setInitialDirection(double v0, double v1, double v2, double v3);

    bool setInitialLocalNullDirection(
        enum_time_direction tdir, const double* l, enum_nat_tetrad_type nattype = enum_nat_tetrad_default);
    bool setInitialLocalNullDirection(enum_time_direction tdir, double l0, double l1, double l2,
        enum_nat_tetrad_type nattype = enum_nat_tetrad_default);

    bool setInitialLocalTimeDirection(enum_time_direction tdir, double l0, double l1, double l2, double beta,
        enum_nat_tetrad_type natType = enum_nat_tetrad_default);

    enum_break_condition calculateGeodesic();
    enum_break_condition calcSachsJacobi();

    void printStatus();

    unsigned int getNumPoints();
    vec4 getPosition(unsigned int num);
    double getAffineParam(unsigned int num);

    void clearAll();
    void resetAll();

    /**
     * @brief Get parameter value of Object.
     * @param paramName
     * @param paramValue
     * @return true : if parameter was found.\n
     *         false : parameter was not found.
     */
    bool getParam(const char* paramName, int& paramValue);

    /**
     * @brief Get parameter value of Object.
     * @param paramName
     * @param paramValue
     * @return true : if parameter was found.\n
     *         false : parameter was not found.
     */
    bool getParam(const char* paramName, double& paramValue);

    /**
     * @brief Get parameter value of Object.
     * @param paramName
     * @param paramValue
     * @return true : if parameter was found.\n
     *         false : parameter was not found.
     */
    bool getParam(const char* paramName, m4d::vec3& paramValue);

    /**
     * @brief Get parameter value of Object.
     * @param paramName
     * @param paramValue
     * @return true : if parameter was found.\n
     *         false : parameter was not found.
     */
    bool getParam(const char* paramName, m4d::vec4& paramValue);

    /**
     * @brief Set Lorentz transformation.
     *  \param  chi : angle in deg.
     *  \param  ksi : angle in deg.
     *  \param  beta : velocity (v/c).
     */
    bool setLorentzTransf(const double chi, const double ksi, const double beta);
    bool setLorentzTransf(const m4d::vec3 beta);

    /// Reset Lorentz transformation.
    void resetLorentzTransf();

    /**
     * @brief Load settings.
     *  \param filename : name of setting file.
     *  \param printset : print setting.
     *  \return true : success.
     *  \return false : error occured.
     */
    bool loadSettings(const char* filename, bool printset = false);

    /*!  Save settings.
     *  \param filename : name of the settings file.
     *  \param dat      : current date.
     *  \return true : success.
     *  \return false : error occured.
     */
    bool saveSettings(const char* filename, const char* dat = nullptr);

    /*! Print settings to fptr.
     *  \param fptr : pointer to file.
     */
    void printSettings(FILE* fptr = stderr);

    /*! Prepare a report for the current metric.
     *  \param text : reference to string.
     *  \return true : success.
     *  \return false : no metric available.
     */
    bool makeReport(char*& text);
    void printReport(FILE* fptr = stdout);

public:
    MetricDatabase metricDB;
    Metric* currMetric;
    IntegratorDatabase solverDB;
    Geodesic* geodSolver;
    enum_integrator geodSolverType;

    enum_geodesic_type type;
    bool stepsizeControlled;
    double stepsize;
    double max_stepsize;
    double min_stepsize;
    double epsAbs;
    double epsRel;
    double epsConstr;
    double epsResize;

    vec4 startPos;
    vec3 startDir;
    double ksi;
    double chi;
    double vel;
    vec4 coordDir;
    vec4 base[4];
    bool isBaseInCoords;

    int axes_orient;
    double boost_ksi;
    double boost_chi;
    double boost_beta;
    mat4 lorentz;

    int timeDirection;
    enum_nat_tetrad_type tetradType;
    unsigned int maxNumPoints; //!< maximum number of points to be calculated

    std::vector<vec4> points; //!< trajectory points
    std::vector<vec4> dirs; //!< trajectory coordinate directions
    std::vector<double> lambda; //!< trajectory affine parameter
    std::vector<vec4> sachs1;
    std::vector<vec4> sachs2;
    std::vector<vec5> jacobi;
    std::vector<vec4> trans_lt[4]; //!< transported local tetrad vector e0...e3
    vec5 maxJacobi;

    double speed_of_light;
    double grav_constant;
    double dielectric_perm;
};

} // end namespace m4d

#endif
