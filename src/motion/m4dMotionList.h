// --------------------------------------------------------------------------------
/*
    m4dMotionList.h

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
// --------------------------------------------------------------------------------
#ifndef  M4D_MOTION_LIST_H
#define  M4D_MOTION_LIST_H

#include <motion/m4dMotion.h>
#include <motion/m4dGeodesic.h>
#include <motion/m4dGeodesicBS.h>
#include <motion/m4dGeodesicGSL.h>
#include <motion/m4dGeodesicRK4.h>
#include <motion/m4dGeodesicDP54.h>
#include <motion/m4dGeodesicDP65.h>

namespace m4d {

#ifdef USE_DP_INT
const unsigned int NUM_GEOD_SOLVERS = 12;
#else
const unsigned int NUM_GEOD_SOLVERS = 10;
#endif
const unsigned int GSL_INT_OFFSET   = 4;


/* ----------------------------------------------
 *   List of all solvers currently implemented
 *
 *   When editing this list please take care
 *   of the ordering.
 * ---------------------------------------------- */
static const char  stl_solver_names[NUM_GEOD_SOLVERS][60] = {
    "Runge-Kutta (fourth order)",
    "Bulirsch-Stoer",
    "GSL: Embedded 2nd order Runge-Kutta",
    "GSL: 4th order (classical) Runge-Kutta",
    "GSL: Embedded 4th order Runge-Kutta-Fehlberg",
    "GSL: Embedded 4th order Runge-Kutta Cash-Karp",
    "GSL: Embedded 8th order Runge-Kutta Prince-Dormand",
    "GSL: Implicit 2nd order Runge-Kutta at Gaussian points",
    //"GSL: Implicit Burlisch-Stoer of Bader and Deuflhard",
    "GSL: M=1 implicit Gear",
    "GSL: M=2 implicit Gear"
#ifdef USE_DP_INT
    , "Dormand-Prince: RK5(4)",
    "Dormand-Prince: RK6(5)"
#endif
};

static const char  stl_solver_nicknames[NUM_GEOD_SOLVERS][30] = {
    "RK4",
    "BS",
    "GSL_RK2",
    "GSL_RK4",
    "GSL_RK_Fehlberg",
    "GSL_RK_Cash-Karp",
    "GSL_RK_Prince-Dormand",
    "GSL_RK_Gaussian-points",
// "GSL_BS_Bader-Deuflhard",
    "GSL_M1_IG",
    "GSL_M2_IG"
#ifdef USE_DP_INT
    , "DP54",
    "DP65"
#endif
};

enum enum_integrator {
    gsUnknown = -1,
    gsIrk4 = 0,       // standard fourth-order Runge-Kutta
    gsInrbs,          // Bulirsch-Stoer
    gsIgslrk2,        // gsl_odeiv_step_rk2
    gsIgslrk4,        // gsl_odeiv_step_rk4
    gsIgslfehlberg,   // gsl_odeiv_step_rkf45
    gsIgslcash,       // gsl_odeiv_step_rkck
    gsIgslprinc,      // gsl_odeiv_step_rk8pd
    gsIi2,            // gsl_odeiv_step_rk2imp
//  gsIbs,            // gsl_odeiv_step_bsimp
    gsIm1,            // gsl_odeiv_step_gear1
    gsIm2             // gsl_odeiv_step_gear2
#ifdef USE_DP_INT
    , gsIdp54,         // Dormand-Prince: RK5(4)
    gsIdp65           // Dormand-Prince: RK6(5)
#endif
};

} // end namespace m4d

#endif

