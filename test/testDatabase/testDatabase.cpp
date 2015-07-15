// -------------------------------------------------------------------------------- 
/*
    testDatabase.cpp

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

/*! \file  testDatabase.cpp
    \brief Print all implemented metrics and integrators.

     Usage:  \verbatim ./m4dTestDatabase \endverbatim
*/
// -------------------------------------------------------------------------------- 

#include <iostream>
#include <cmath>

#include <m4dGlobalDefs.h>
#include <metric/m4dMetricDatabase.h>
#include <motion/m4dMotionList.h>

/* ----------------------------------------------------------------- 
 *                     m a i n
 * ----------------------------------------------------------------- */
int main( int argc, char* argv[] )
{
    fprintf(stderr,"\n\n=================== Test Database ===================\n\n");

    m4d::MetricDatabase* metricDB = m4d::MetricDatabase::getInstance();
    metricDB->printMetricList(stdout);

    fprintf(stdout,"\nnum   Integrator name\n-----------------------------------------------------------\n");
    for(unsigned int i=0; i<m4d::NUM_GEOD_SOLVERS; i++) {
        fprintf(stdout,"%2d   %s\n",i, m4d::stl_solver_names[i]);
    }

    fprintf(stdout,"\nnum   name of tetrad type\n-----------------------------------------------------------\n");
    for(unsigned int i=0; i<m4d::NUM_ENUM_NAT_TETRAD_TYPES; i++) {
        fprintf(stdout,"%2d   %s\n",i, m4d::stl_nat_tetrad_types[i]);
    }    
    return 1;
}

