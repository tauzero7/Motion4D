// -------------------------------------------------------------------------------- 
/*
    testTidalAccel.cpp
 
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

/*! \file  testTidalAccel.cpp
    \brief Test tidal acceleration for timelike static observer in the
           Schwarzschild metric.

     Usage: \verbatim ./m4dTestTidalAccel  <radius>\endverbatim


     
*/
// -------------------------------------------------------------------------------- 

#include <iostream>
#include <cmath>

#ifdef _WIN32
#include "stdafx.h"
#endif

#include <m4dGlobalDefs.h>
#include <metric/m4dMetricDatabase.h>
#include <motion/m4dGeodesic.h>
#include <motion/m4dGeodesicRK4.h>
#include <math/TransCoordinates.h>
#include <extra/m4dObject.h>


/* ----------------------------------------------------------------- 
 *                     m a i n
 * ----------------------------------------------------------------- */
int main( int argc, char* argv[] )
{
  if (argc!=2)
  {
    fprintf(stderr,"Usage: ./m4dTestTidalAccel  <radius>\n");
    return -1;
  }
  double radius = atof(argv[1]);
  
  fprintf(stderr,"\n\n=================== Test tidal acceleration ===================\n\n");  
  
  
  /* -----------------------------------------
   *   Initialize metric.
   * ----------------------------------------- */
  m4d::MetricDatabase* metricDB = new m4d::MetricDatabase;
  m4d::Metric* metric = metricDB->getMetric("Schwarzschild");
  metric->print(); 

    
  /* -----------------------------------------
   *    Set initial position.
   * ----------------------------------------- */
  m4d::vec4 pos = m4d::vec4(0.0,radius,M_PI_2,0.0);
  

  /* -----------------------------------------
   *  Get natural local tetrad at pos.
   * ----------------------------------------- */
  m4d::vec4 e0,e1,e2,e3;  
  metric->getNatTetrad(pos,e0,e1,e2,e3);
  std::cerr << std::endl;
  std::cerr << "e0: "; e0.printS();
  std::cerr << "e1: "; e1.printS();
  std::cerr << "e2: "; e2.printS();
  std::cerr << "e3: "; e3.printS();
  std::cerr << std::endl;
  
  metric->calcTidalMatrix(pos,e0,e1,e2,e3);
  m4d::mat4 tM;
  metric->getTidalMatrix(tM);
  tM.printS();
  
  delete metric;
  delete metricDB;
  return 1;
}

