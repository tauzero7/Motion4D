// -------------------------------------------------------------------------------- 
/*
     m4dGlobalDefs.h     is part of the m4d-library.

   Copyright (c) 2009-2014  Thomas Mueller, Frank Grave

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 3 of the License, or (at
   your option) any later version.
   
   This program is distributed in the hope that it will be useful, but
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   General Public License for more details.
   
   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.

*/

/*!  \file   m4dGlobalDefs.h
     \brief  Global definitions for the m4d library.

           If the constraint equation
         \f[  g_{ab} \dot{x}^a \dot{x}^b - \kappa  > \mbox{DEF\_CONSTRAINT\_EPSILON} \f]
           then the integration stops.

     \mainpage  The Motion-4D-Library.
     
    <H2>Introduction</H2>
    <p>The Motion4D library solves the geodesic equation as well as the parallel- and Fermi-Walker-transport in four-dimensional Lorentzian spacetimes numerically. Initial conditions are given with respect to natural local tetrads which are adapted to the symmetries or the coordinates of the spacetime. Beside some already implemented metrics like the Schwarzschild and Kerr metric, the object oriented structure of the library permits to implement other metrics or integrators in a straight forward manner.</p>
    <p>A more detailed description of the implemented metrics can be found either in the original literature or in the \ref lit_cos "Catalogue of Spacetimes".</p>


    <H2>Add own metric</H2>
    <p>In order to add a metric, one has to do the following:</p>
    <ol>
      <li>Save an old metric under the new metric name:<br> e.g. m4dMetricSchwarzschild.* -> m4dMetricNewName.*
      <li>Write all the necessary methods.
      <li>Add the new metric header into the m4dMetricList.h file and adjust NUM_METRICS, stl_metric_names, and enum_metric.
      <li>Add the new standard metric constructor to the method initializeMetric() in the m4dMetricDatabase.cpp file.
      <li>Add new filenames to the src/metric/Makefile.am  file.
      <li>autoreconf, configure, make
    </ol>
    
    
    <H2>libMotion4D License</H2>
    <p>Copyright &copy; 2009-2014 by Thomas M&uuml;ller and Frank Grave</p><br>
    <p>Permission to use, copy, modify, and distribute this software and its documentation under the terms of the GNU General Public License is hereby granted. No representations are made about the suitability of this software for any purpose. It is provided "as is" without express or implied warranty. See the <a href="http://www.gnu.org/copyleft/gpl.html">GNU General Public License</a> for more details. </p>
    
    <br>
    <H2>Contact</H2>
    Visualisierungsinstitut der Universit&auml;t Stuttgart<br>
    Allmandring 19<br>
    70569 Stuttgart, Germany<br>
    Email: Thomas.Mueller@vis.uni-stuttgart.de<br><br>
    
    Universit&auml;t Stuttgart<br>
    1. Institut für Theoretische Physik<br>
    Pfaffenwaldring 57 // IV<br>
    70550 Stuttgart, Germany<br>
    <br>

    <H2>Books</H2>
    <ul>
      <li>\anchor lit_chandra Chandrasekhar, S., <b>The Mathematical Theory of Black Holes</b> (Oxford University Press)
      <li>\anchor lit_mtw Misner, C.W., and Thorne, K.S., and Wheeler, J.A., <b>Gravitation</b> (W.H.Freeman and Company, New York, 1973)
      <li>\anchor lit_rindler Rindler, W., <b>Relativity - Special, General and Cosmology</b> (Oxford University Press)
      <li>\anchor lit_wald Wald, R., <b>General Relativity</b> (University of Chicago Press)
      <li>\anchor lit_stephani Stephani, H., Kramer, D., MacCallum, M., Hoenselaers, C. and Herlt, E., <b>Exact Solutions of the Einstein Field Equations</b> (Cambridge University Press, 2. edition, 2009)
    </ul>

    <H2>Links</H2>
    <ul>
      <li>M&uuml;ller, T. and Grave, F.<br>
          <i><b>GeodesicViewer - A tool for exploring geodesics in the theory of relativity</b></i><br>
          Comput. Phys. Commun. <b>181</b>, 413-419 (2010).  DOI: <a href="http://dx.doi.org/10.1016/j.cpc.2009.10.010" target="_blank">10.1016/j.cpc.2009.10.010</a>
      <li>M&uuml;ller, T. and Grave, F.<br>
          <i><b>Motion4D - A library for lightrays and timelike worldlines in the theory of relativity</b></i><br>
          Comput. Phys. Commun. <b>180</b>, 2355-2360 (2009).  DOI: <a href="http://dx.doi.org/10.1016/j.cpc.2009.07.014" target="_blank">10.1016/j.cpc.2009.07.014</a>
      <li>M&uuml;ller, T. and Grave, F.<br>
          <i><b>An updated version of the Motion4D library</b></i><br>
          Comput. Phys. Commun. <b>181</b>, 703 (2010). DOI: <a href="http://dx.doi.org/10.1016/j.cpc.2009.10.021" target="_blank">10.1016/j.cpc.2009.10.021</a>
      <li>M&uuml;ller, T.<br>
          <i><b>Motion4D-library extended</b></i><br>
          Comput. Phys. Commun. <b>182</b>, 1386 (2011). DOI: <a href="http://dx.doi.org/10.1016/j.cpc.2011.02.009" target="_blank">10.1016/j.cpc.2011.02.009</a>
      <li>\anchor lit_cos M&uuml;ller, T. and Grave, F.<br>
          <i><b>Catalogue of Spacetimes</b></i><br>
          <a href="http://xxx.lanl.gov/abs/0904.4184v1" target="_blank">arXiv:0904.4184</a> [gr-qc].
      <li><a href="http://www.gnu.org/software/gsl/" target="_blank">GNU Scientific Library</a>
    </ul>
 */

// -------------------------------------------------------------------------------- 
#ifndef M4D_GLOBAL_DEFS_H
#define M4D_GLOBAL_DEFS_H

#include <iostream>
#include <vector>
#include <cmath>
#include <float.h>

#ifdef _WIN32
#include <windows.h>
#endif

#ifdef _WIN32
#ifdef METRIC_EXPORTS
#define METRIC_API __declspec(dllexport)
#else /* METRIC_EXPORTS */
#define METRIC_API __declspec(dllimport)
#endif /* METRIC_EXPORTS */
#else /* _WIN32 */
#define METRIC_API
#endif /* _WIN32 */


#ifdef _WIN32
#ifdef MOTION_EXPORTS
#define MOTION_API __declspec(dllexport)
#else /* METRIC_EXPORTS */
#define MOTION_API __declspec(dllimport)
#endif /* METRIC_EXPORTS */
#else /* _WIN32 */
#define MOTION_API
#endif /* _WIN32 */

#ifdef _WIN32
#ifdef MATH_EXPORTS
#define MATH_API __declspec(dllexport)
#else /* METRIC_EXPORTS */
#define MATH_API __declspec(dllimport)
#endif /* METRIC_EXPORTS */
#else /* _WIN32 */
#define MATH_API
#endif /* _WIN32 */

#ifdef _WIN32
#ifdef EXTRA_EXPORTS
#define EXTRA_API __declspec(dllexport)
#else /* METRIC_EXPORTS */
#define EXTRA_API __declspec(dllimport)
#endif /* METRIC_EXPORTS */
#else /* _WIN32 */
#define EXTRA_API
#endif /* _WIN32 */

#ifdef _WIN32
#define  isnan(x) ((x) != (x))
#define  asinh(x)   (log(x + sqrt(x*x+1)))
#endif

#ifdef __APPLE__
#define  isnan(x) ((x) != (x))
#define  asinh(x)   (log(x + sqrt(x*x+1)))
#endif

#ifdef _WIN32
#pragma warning(disable: 4251 4273)
#endif


// compare /usr/include/math.h
#ifndef M_PI
#define M_PI           3.14159265358979323846  /* pi */
#define M_PI_2         1.57079632679489661923  /* pi/2 */
#define M_PI_4         0.78539816339744830962  /* pi/4 */
#endif

#ifndef RAD_TO_DEG
#define RAD_TO_DEG  57.295779513082320875
#define DEG_TO_RAD  0.017453292519943295770
#endif

#ifndef DBL_MAX
#define DBL_MAX 1.844674407370955616e19
#endif

#define DEF_FIXED_REPORT_PRECISION  14

#define M4D_SIGN(x)    (x>0 ? 1.0 : -1.0)
#define M4D_DELTA(x,y) (x==y ? 1.0 : 0.0)
#define M4D_MAX(x,y)   (x>y ? x : y)
#define M4D_MIN(x,y)   (x>y ? y : x)
#define M4D_SQR(x)     ((x)*(x))

#define DEF_MAX_STEPSIZE        1.0
#define DEF_MIN_STEPSIZE        1.0e-15

#define DEF_CONSTRAINT_EPSILON  1.0e-6
#define DEF_RESIZE_EPSILON      1.0
#define DEF_RESIZE_FACTOR       1.0

#define DEF_MAX_YS_PAR  24
#define DEF_MAX_YS_JAC  32

#define DEF_MAX_YS      32

// Start index position of tangent vector
#define DEF_TG_IDX       4

// Start index positions of local tetrad vectors
#define DEF_EO_IDX   8
#define DEF_E1_IDX  12
#define DEF_E2_IDX  16
#define DEF_E3_IDX  20

// Start index positions of Sachs basis vectors
#define DEF_SA1_IDX   8    
#define DEF_SA2_IDX  12

// Start index positions of Jacobi vectors
#define DEF_JAC1_IDX  16
#define DEF_DJ1_IDX   20
#define DEF_JAC2_IDX  24
#define DEF_DJ2_IDX   28

/* --------------------------------------------------------
 *   mathematical stuff
 * -------------------------------------------------------- */
#include <math/VnD.h>
#include <math/Mat.h>

namespace m4d
{

typedef VnD<double,2>      vec2;
typedef VnD<double,3>      vec3;
typedef VnD<double,4>      vec4;
typedef VnD<double,5>      vec5;
typedef VnD<double,9>      vec9;

typedef Matrix<double,2,2> mat2;
typedef Matrix<double,3,3> mat3;
typedef Matrix<double,4,4> mat4;

typedef VnD<float,3>       vec3f;
typedef VnD<float,4>       vec4f;
typedef VnD<float,9>       vec9f;

typedef Matrix<float,3,3>  mat3f;
typedef Matrix<float,4,4>  mat4f;

typedef VnD<int,2>         ivec2;
typedef VnD<int,3>         ivec3;
typedef VnD<int,4>         ivec4;

typedef enum  _axisID
{
    axis_X,
    axis_Y,
    axis_Z
} enum_axisID;


//  1 greg. year = 365.2425 days 

/* ---------------------------------------------------------------
 *   http://nssdc.gsfc.nasa.gov/planetary/factsheet/sunfact.html
 * --------------------------------------------------------------- */
#define  M4D_MASS_OF_THE_SUN      1.9891e30

/* ---------------------------------------------------------------
 *   physical constants:
 *   http://physics.nist.gov/cuu/Constants/index.html
 * --------------------------------------------------------------- */
#define M4D_SPEED_OF_LIGHT                      299792.4580              //  km s-1
#define M4D_SPEED_OF_LIGHT_IN_PARSEC_PER_YEAR   3.067484662576687e-1     //  1/3.26

#define M4D_GRAV_CONST                            6.67428e-20      //  km^3 kg^-1 s^-2
#define M4D_GRAV_CONST_IN_LS_PER_KG_PER_SEC       2.477093e-36     //  ls^3 kg^-1 s^-2
#define M4D_GRAV_CONST_IN_LS_PER_MSUN_PER_SEC     4.927186e-6      //  ls^3 Msun^-1 s^-2
#define M4D_GRAV_CONST_IN_KM_PER_MSUN_PER_SEC     1.327581e+11     //  km^3 Msun^-1 s^-2
#define M4D_GRAV_CONST_IN_KM_LS_PER_MSUN_PER_SEC  1.477133         //  km ls^2 Msun^-1 s^-2

#define M4D_DIELECTRIC_PERM                       8.854187817e-3   //  C^2 s^2 kg^-1 km^-3

#define M4D_YEAR_TO_SEC       3.1556952e+7            //  365.2425*24*3600 s
#define M4D_SEC_TO_YEAR       3.168873850681143e-8

#define M4D_LIGHTYEAR         9.460536207068016e+15   //  M4D_YEAR_TO_SEC * M4D_SPEED_OF_LIGHT
#define M4D_PARSEC            3.084134803504173e+16   //  M4D_LIGHTYEAR * 3.26

#define M4D_METER_TO_LS       3.33564095198152e-9




/* --------------------------------------------------------
 *   use physical constants:
 *      notset :
 *      geom   :  use geometric units: G = c = 1
 *      real   :  use real physical units defined above
 *      proper :  use proper physical units, e.g. c=30 km/h
 * -------------------------------------------------------- */
enum  enum_physical_constants
{
    enum_physical_constants_notset = 0,
    enum_physical_constants_geom,
    enum_physical_constants_real,
    enum_physical_constants_proper
};


/* --------------------------------------------------------
 *   standard coordinates :
 *
 *    cartesian :           [t, x, y, z]
 *    spherical :           [t, r, theta, phi]
 *    cylinder  :           [t, r, phi, z]
 *    prolate spheroidal:   [t, sigma, tau, phi]
 * -------------------------------------------------------- */
const unsigned int NUM_ENUM_COORDINATE_TYPES = 5;

enum  enum_coordinate_type
{
    enum_coordinate_cartesian = 0,
    enum_coordinate_spherical,
    enum_coordinate_cylinder,
    enum_coordinate_prolatespheroidal,
    enum_coordinate_custom
};

const char stl_coordinate_types[NUM_ENUM_COORDINATE_TYPES][18] =
{  "cartesian",
   "spherical",
   "cylinder",
   "prolatespheroidal",
   "custom"
};

const int NUM_ENUM_SINGLE_COORDINATE_TYPES = 4;	 

enum  enum_single_coordinate_type
{
    enum_scoord_linear = 0,
    enum_scoord_semilinear,
    enum_scoord_periodic,
    enum_scoord_range
};

 const char stl_single_coordinate_types[NUM_ENUM_SINGLE_COORDINATE_TYPES][12] =
{
    "linear", "semilinear", "periodic", "range"
};

const int NUM_ENUM_COORDINATE_CHARACTERS = 4;

enum  enum_coordinate_character
{
    enum_cchar_lightlike = 0,
    enum_cchar_timelike,
    enum_cchar_spacelike,
    enum_cchar_unknown
};

const char stl_coordinate_characters[NUM_ENUM_COORDINATE_CHARACTERS][10] =
{
    "lightlike", "timelike", "spacelike", "unknown"
};

/*! \brief Single coordinate type.
 */
typedef struct _scoord_type
{
    enum_single_coordinate_type  type;
    double                       min;
    double                       max;
    enum_coordinate_character    character;
} struct_scoord_type;

/* --------------------------------------------------------
 *   classification follows  D.Bini and R.T.Jantzen
 *   "Circular orbits in Kerr spacetime...",
 *   CQG 17 (2000) 1637-1647
 * -------------------------------------------------------- */
const unsigned int NUM_ENUM_NAT_TETRAD_TYPES = 7;

enum  enum_nat_tetrad_type
{
    enum_nat_tetrad_default = 0,
    enum_nat_tetrad_static,        // static with respect to coordinates
    enum_nat_tetrad_lnrf,          // locally nonrotating frame
    enum_nat_tetrad_freefall,      // free fall observer from infinity
    enum_nat_tetrad_comoving,      // comoving frame
    enum_nat_tetrad_cylinder,
    enum_nat_tetrad_cartesian
};

const char stl_nat_tetrad_types[][10]  =
{  "default",
   "static",
   "lnrf",
   "freefall",
   "comoving",
   "cylinder",
   "cartesian"
};

/* --------------------------------------------------------
 *   predefined tetrads
 * -------------------------------------------------------- */
const int NUM_PREDEF_TETRADS = 4;

enum  enum_predef_tetrad
{
    enum_predef_lt_standard = 0,
    enum_predef_lt_righthanded,
    enum_predef_lt_radrev,
    enum_predef_lt_invert
};

const char  stl_predef_tetrads[NUM_PREDEF_TETRADS][14] =
{  "standard",
   "right handed",
   "inward",
   "inverted"
};

const double vec_predef_tetrad[NUM_PREDEF_TETRADS][16] =
{
    {1.0,0.0,0.0,0.0, 0.0,1.0,0.0,0.0,  0.0,0.0,1.0,0.0, 0.0,0.0,0.0,1.0},
    {1.0,0.0,0.0,0.0, 0.0,1.0,0.0,0.0,  0.0,0.0,0.0,1.0, 0.0,0.0,-1.0,0.0},
    {1.0,0.0,0.0,0.0, 0.0,-1.0,0.0,0.0, 0.0,0.0,0.0,1.0, 0.0,0.0,1.0,0.0},
    {1.0,0.0,0.0,0.0, 0.0,0.0,0.0,1.0,  0.0,0.0,1.0,0.0, 0.0,1.0,0.0,0.0}
};

/* --------------------------------------------------------
 *   Minkowski metric
 * -------------------------------------------------------- */
#define eta(j,k)  ((j==k) ? ((j==0)?-1.0:1.0) : 0.0)



/* --------------------------------------------------------
 *   break conditions for numerical integration
 * -------------------------------------------------------- */
const int NUM_ENUM_BREAK_CONDITIONS = 8;

enum  enum_break_condition
{
    enum_break_none = 0,
    enum_break_cond,             /*!< geodesic comes into a region where spacetime break down,... */
    enum_break_outside,          /*!< geodesic is outside bounding box */
    enum_break_num_exceed,       /*!< maximum number of geodesic points exceeded */
    enum_break_step_size,        /*!< stepsize too small */
    enum_break_constraint,       /*!< constraint condition is violated */
    enum_break_not_implemented,  /*!< not implemented */
    enum_break_other             /*!< other error occured */
};

const char stl_break_condition[NUM_ENUM_BREAK_CONDITIONS][34] =
{  "no break",
   "spacetime breaks down",
   "outside bounding box",
   "maximum number of points exceeded",
   "stepsize too small",
   "constraint condition violated",
   "not implemented",
   "other error occured"
};


/* --------------------------------------------------------
 *   Type of the geodesic:
 * -------------------------------------------------------- */
const unsigned int NUM_ENUM_GEODESIC_TYPE = 4;

enum  enum_geodesic_type
{
    enum_geodesic_lightlike = 0,
    enum_geodesic_lightlike_sachs,
    enum_geodesic_timelike,
    enum_geodesic_spacelike
};

const char stl_geodesic_type[NUM_ENUM_GEODESIC_TYPE][16] =
{ "lightlike",
  "lightlike_sachs",
  "timelike",
  "spacelike"
};


/* --------------------------------------------------------
 *   Time direction
 * -------------------------------------------------------- */
enum  enum_time_direction
{
    enum_time_backward = -1,
    enum_time_forward = 1
};
const char  str_time_dir[2][9] = {"backward","forward"};


/* --------------------------------------------------------
 *   drawing type
 * -------------------------------------------------------- */
const int NUM_ENUM_DRAW_TYPE = 6;

enum  enum_draw_type
{
    enum_draw_pseudocart = 0,
    enum_draw_coordinates,
    enum_draw_embedding,
    enum_draw_twoplusone,
    enum_draw_effpoti,
    enum_draw_custom
};

const char stl_draw_type[NUM_ENUM_DRAW_TYPE][20] =
{ "pseudo cartesian",
  "coordinates",
  "embedding diagram",
  "2+1 diagram",
  "effective potential",
  "custom"
};




}  // end namespace m4d

#endif
