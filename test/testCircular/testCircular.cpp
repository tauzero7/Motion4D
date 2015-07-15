// -------------------------------------------------------------------------------- 
/*
    testCircular.cpp
 
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

/*! \file  testCircular.cpp
    \brief Test parallel transport for timelike stable circular orbits around a
           Schwarzschild black hole. Note that stable circular orbits can only
           exist for r > 1.5rs.

     Usage: \verbatim ./m4dTestCircular  <radius>\endverbatim

     Output:  lambda  pos[0] pos[1] pos[3]  dirX1 dirY1 angle1

     where
     <ul>
       <li>lambda : proper time
       <li>pos[0] : coordinate time
       <li>pos[1] : coordinate position: r
       <li>pos[3] : coordinate position: phi
       <li>dirX1  : pseudo-cartesian x-direction of e1
       <li>dirY1  : pseudo-cartesian y-direction of e1
       <li>angle1 : angle of e1 in pseudo-cartesian coordinates.
     </ul>     
*/
// -------------------------------------------------------------------------------- 

#include <iostream>
#include <cmath>

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
    fprintf(stderr,"Usage: ./m4dTestCircular  <radius>\n");
    return -1;
  }
    
  fprintf(stderr,"\n\n=================== Test parallel transport on circular orbit ===================\n\n");  
  
  // Schwarzschild radius.
  double rs = 2.0;

  // Radius for circular orbit.
  double radius = atof(argv[1]);
  if (radius<3.0*rs)
  {
    fprintf(stderr,"Radius must be greater than 3*rs = %10.6f !\n",3.0*rs);
    return -1;
  }
  
  // Local velocity for circular orbit.
  double vel = sqrt(0.5*rs/(radius-rs));
  fprintf(stderr,"\nlocal velocity:            beta = %10.6f\n",vel);

  // Proper time for one orbit.
  double tau = 2.0*M_PI*sqrt(2.0*pow(radius,3.0)*(1.0-1.5*rs/radius)/rs);
  fprintf(stderr,"proper time for one orbit:  tau = %10.6f\n",tau);

  // Coordinate time for one orbit.
  double ct = 2.0*M_PI*sqrt(2.0*pow(radius,3.0)/rs);
  fprintf(stderr,"coord. time for one orbit:   ct = %10.6f\n\n",ct);

  m4d::Object  mObject;
  /* -----------------------------------------
   *   Initialize metric.
   * ----------------------------------------- */
  m4d::MetricDatabase* metricDB = m4d::MetricDatabase::getInstance();
  m4d::Metric* metric = metricDB->getMetric("Schwarzschild");
  metric->printF(); 

  
  /* -----------------------------------------
   *    Set geodesic integrator.
   * ----------------------------------------- */
  m4d::GeodesicRK4* geod = new m4d::GeodesicRK4( metric, m4d::enum_geodesic_timelike );  
  geod->setBoundingBox(m4d::vec4(-DBL_MAX,-DBL_MAX,-DBL_MAX,-DBL_MAX),m4d::vec4(DBL_MAX,DBL_MAX,DBL_MAX,DBL_MAX));
  geod->setAffineParamStep(0.001);

  
  /* -----------------------------------------
   *    Set initial position and direction.
   * ----------------------------------------- */
  m4d::vec4 pos = m4d::vec4(0.0,radius,M_PI_2,0.0);

  double y0  = 1.0/sqrt(1.0-vel*vel);
  double eta = vel*y0;

  m4d::vec4 coDir;
  m4d::vec4 locDir = m4d::vec4(y0,0.0,0.0,eta);
  metric->localToCoord(pos,locDir,coDir,m4d::enum_nat_tetrad_default);


  /* -----------------------------------------
   *  Boost transformation for initial tetrad.
   * ----------------------------------------- */
  m4d::vec4 e0,e1,e2,e3,u,a;  
  mObject.setLorentzTransf(0.0,90.0,vel);
  mObject.lorentz.printS();

  fprintf(stderr,"gamma: %f    beta/gamma: %f\n",1.0/sqrt(1.0-vel*vel),vel/sqrt(1.0-vel*vel));

  // Determine initial local tetrad.
  metric->localToCoord(pos,mObject.lorentz*m4d::vec4(1,0,0,0),e0);
  metric->localToCoord(pos,mObject.lorentz*m4d::vec4(0,1,0,0),e1);
  metric->localToCoord(pos,mObject.lorentz*m4d::vec4(0,0,1,0),e2);
  metric->localToCoord(pos,mObject.lorentz*m4d::vec4(0,0,0,1),e3);

  /* -----------------------------------------
   *  Set initial position, direction, 
   *  and tetrad.
   * ----------------------------------------- */
  geod->setInitialPosition(pos);
  geod->setInitialDirection(coDir);
  geod->setInitialTetrad(e0,e1,e2,e3);


  /* -----------------------------------------
   *  Calculate the geodesic !
   * ----------------------------------------- */
  m4d::mat4 initTetrad;
  m4d::vec4 n0,n1,n3;
  double dirX1,dirY1,angle1;
  
  double* currY = new double[DEF_MAX_YS];
  double cstr = 0.0, lambda;
  int status;
  while (pos[3]<=2.0*M_PI)
  {
    // Determine natural local tetrad at current position.
    metric->localToCoord(pos,mObject.lorentz*m4d::vec4(1,0,0,0),e0);
    metric->localToCoord(pos,mObject.lorentz*m4d::vec4(0,1,0,0),e1);
    metric->localToCoord(pos,mObject.lorentz*m4d::vec4(0,0,1,0),e2);
    metric->localToCoord(pos,mObject.lorentz*m4d::vec4(0,0,0,1),e3);
      
    // Calculate dual local tetrad for coordToLocal transformation.
    initTetrad.setCol(0,e0);
    initTetrad.setCol(1,e1);
    initTetrad.setCol(2,e2);
    initTetrad.setCol(3,e3);
    initTetrad.invert();
  
    n0  = initTetrad*(geod->getE(0));
    n1  = initTetrad*(geod->getE(1));
    n3  = initTetrad*(geod->getE(3));
	
	
    lambda = geod->getAffineParam();
    fprintf(stderr,"\rtau = %8.3f",lambda);
  //  fprintf(stdout,"%8.3f   %12.8f %12.8f %12.8f %12.8f   %12.8f %12.8f %12.8f %12.8f   %12.8f %12.8f %12.8f %12.8f\n",lambda,pos[0],pos[1],pos[2],pos[3],e1[0],e1[1],e1[2],e1[3],e3[0],e3[1],e3[2],e3[3]);
    
    dirX1 = n1[1]*cos(pos[3]) - n1[3]*sin(pos[3]);
    dirY1 = n1[1]*sin(pos[3]) + n1[3]*cos(pos[3]);    
    angle1 = atan2(dirY1,dirX1);
    
    fprintf(stdout,"%8.3f  %12.8f %12.8f %12.8f   %12.8f %12.8f  %12.8f\n",lambda,pos[0],pos[1],pos[3]*RAD_TO_DEG,dirX1,dirY1,angle1*RAD_TO_DEG);
    
    geod->nextStepPar(status);
    geod->getCurrentArray(currY);
    if ((cstr=fabs(metric->testConstraint(currY,-1))) > geod->getConstrEps())
        break;
    pos = geod->getPosition();
  }
    
  fprintf(stderr,"\nGeodesic precession of e_1 after one orbit: %f deg.\n\n",atan2(n1[3],n1[1])*RAD_TO_DEG);

  delete geod;
  delete metric;
  return 1;
}

