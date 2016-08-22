// -------------------------------------------------------------------------------- 
/*
    testCPmotion.cpp
 
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

#include <iostream>
#include <cmath>

#include <m4dGlobalDefs.h>
#include <metric/m4dMetricDatabase.h>
#include <motion/m4dMotionChargedParticle.h>
#include <math/TransCoordinates.h>
#include <extra/m4dObject.h>


/* ----------------------------------------------------------------- 
 *                     m a i n
 * ----------------------------------------------------------------- */
int main( int argc, char* argv[] ) {
    if (argc!=1) {
        fprintf(stderr,"Usage: ./m4dTestCPmotion\n");
        return -1;
    }
  
    m4d::Object  mObject;

    /* ----------------------------------
    *  initialize metric database
    * ---------------------------------- */  
    m4d::MetricDatabase* metricDB = m4d::MetricDatabase::getInstance();

    m4d::Metric* metric = metricDB->getMetric("Ernst");
    metric->setParam("mass", 1.0);
    metric->setParam("b", 0.1);
    metric->printF(); 

    /* ----------------------------------
    *  Define initial position.
    * ---------------------------------- */  
    m4d::vec4 initPos = m4d::vec4(0.0, 4.0, M_PI_2, 0.0);
  
   
    /* ----------------------------------
    *  Set FermiWalker integrator.  
    * ---------------------------------- */     
    m4d::MotionChargedParticle* motion = new m4d::MotionChargedParticle( metric );
    motion->setAffineParamStep(0.01);
  
  
    /* ----------------------------------
    *  Set initial position, direction, and local tetrad.
    * ---------------------------------- */         
    motion->setInitialPosition(initPos);
    motion->setInitialVelocity(1, 0.5, 0.0, M_PI_2, 0.1);
    
  
  /* ----------------------------------
   *  Calculate one orbit.
   * ---------------------------------- */    
   int maxNumPoints = 2000;
    m4d::vec4 pos,ne;  
    int i=0; 
    while (i<maxNumPoints) {
        pos = motion->getPosition();
        fprintf(stderr,"\r%5d  %8.4f",i,motion->getAffineParam());
        fprintf(stdout,"%5d  %8.4f  %12.8f %12.8f %12.8f %12.8f \n",i,motion->getAffineParam(),pos[0],pos[1],pos[2],pos[3]);
        motion->nextStep();
        i++;
  }
  fprintf(stderr,"\n");
    
  delete motion;
  delete metric;
  return 1;
}

