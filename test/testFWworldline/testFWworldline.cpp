// -------------------------------------------------------------------------------- 
/*
    testFWworldline.cpp
 
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

/*! \file  testFWworldline.cpp
    \brief Test Fermi-Walker transport along a given worldline.

          Consider an object on a timelike circular orbit with radius 'R'
          around a Schwarzschild black hole of mass 'mass' with angular
          frequency 'omega'.

          As an example: mass=1, R=6, omega=1/(6*sqrt(3)), which equals
          the last stable orbit (geodesic).

          After one orbit, the tetrad has undergone a geodesic precession
          with angle 'angle'. For the above example, the angle is approximately
          105 degree.

          To be on a timelike geodesic, the four acceleration must vanish identically. 
          Thus,
           \f[a^1 = \frac{3}{2}r_s\omega^2+\frac{c^2r_s}{2R^2}-R\omega^2 \stackrel{!}{=} 0. \f]

          For a given radius, the angular velocity must be
          \f[ \omega = c\sqrt{\frac{r_s}{2R^2\left(R-\frac{3}{2}r_s\right)}}\f]
          and the corresponding three velocity is
          \f[ v = c\sqrt{\frac{r_s}{2\left(R-r_s\right)}}.\f]

     Usage: \verbatim ./m4dTestFWworldline \endverbatim
 
*/
// -------------------------------------------------------------------------------- 

#include <iostream>
#include <cmath>

#include <m4dGlobalDefs.h>
#include <metric/m4dMetricDatabase.h>
#include <motion/m4dFermiWalker.h>
#include <math/TransCoordinates.h>
#include <extra/m4dObject.h>

/*! \brief Parameters for FWworldline test.
 */
typedef struct
{
  double rs;
  double omega;
  double R;
} mystruct;

/*! Worldline for a circular orbit around a Schwarzschild black hole.
 *
 *  \param tau : current proper time.
 *  \param params : pointer to parameters.
 */
m4d::vec4  x_tau( double tau, void* params )
{
  if (params==NULL)
    return m4d::vec4();
  
  mystruct* par = (mystruct*)params;
  double rs    = par->rs;
  double omega = par->omega;
  double R     = par->R;
  
  double f = sqrt((omega*omega*R*R+1.0)/(1.0-rs/R));
  return m4d::vec4( f*tau, R, M_PI_2, omega*tau );
}


/*! Tangent to the worldline for a circular orbit around a Schwarzschild black hole.
 *
 *  \param tau : current proper time.
 *  \param params : pointer to parameters.
 */
m4d::vec4  u_tau( double tau, void* params )
{
  if (params==NULL)
    return m4d::vec4();
  
  mystruct* par = (mystruct*)params;
  double rs    = par->rs;
  double omega = par->omega;
  double R     = par->R;
  
  double f = sqrt((omega*omega*R*R+1.0)/(1.0-rs/R));
  return m4d::vec4( f, 0.0, 0.0, omega );
}

/*! Acceleration needed to stay on the worldline for a circular orbit around a Schwarzschild black hole.
 *
 *  \param tau : current proper time.
 *  \param params : pointer to parameters.
 */
m4d::vec4  a_tau( double tau, void* params )
{
  if (params==NULL)
    return m4d::vec4();
  
  mystruct* par = (mystruct*)params;
  double rs    = par->rs;
  double omega = par->omega;
  double R     = par->R;
  
  //double f = sqrt((omega*omega*R*R+1.0)/(1.0-rs/R));
  return m4d::vec4( 0.0, 1.5*rs*omega*omega + 0.5*rs/(R*R) - R*omega*omega, 0.0, 0.0 );
}


/* ----------------------------------------------------------------- 
 *                     m a i n
 * ----------------------------------------------------------------- */
int main( int argc, char* argv[] )
{
  if (argc!=1)
  {
    fprintf(stderr,"Usage: ./m4dTestFWworldline\n");
    return -1;
  }
  
  fprintf(stderr,"\n\n=================== Test Fermi-Walker worldline ===================\n\n");
  double mass  = 1.0;
  double omega = 1.0/6.0/sqrt(3.0);
  double R     = 6;
  int    maxNumPoints = 100000;


  /* ----------------------------------
   *  local velocity
   * ---------------------------------- */
  double local_v = omega*R/sqrt(omega*omega*R*R+1.0);
  fprintf(stderr,"local velocity: %f\n",local_v);
  
  
  m4d::Object  mObject;

  /* ----------------------------------
   *  initialize metric database
   * ---------------------------------- */  
  m4d::MetricDatabase metricDB;

  m4d::Metric* metric = metricDB.getMetric("Schwarzschild");
  metric->setParam("mass",mass);
  metric->printF(); 

  /* ----------------------------------
   *  Define initial position.
   * ---------------------------------- */  
  m4d::vec4 initPos = m4d::vec4(0.0, R, M_PI_2, 0.0);
  
  /* ----------------------------------
   *  Calculate boost transformation for initial tetrad.
   * ---------------------------------- */ 
  m4d::vec4 e0,e1,e2,e3,u,a;  
  mObject.setLorentzTransf(0.0,90.0,local_v);
  
  /* ----------------------------------
   *  Determine initial local tetrad.
   * ---------------------------------- */    
  metric->localToCoord(initPos,mObject.lorentz*m4d::vec4(1,0,0,0),e0);
  metric->localToCoord(initPos,mObject.lorentz*m4d::vec4(0,1,0,0),e1);
  metric->localToCoord(initPos,mObject.lorentz*m4d::vec4(0,0,1,0),e2);
  metric->localToCoord(initPos,mObject.lorentz*m4d::vec4(0,0,0,1),e3);

  /* ----------------------------------
   *  Calculate dual local tetrad for coordToLocal transformation.
   * ---------------------------------- */     
  m4d::mat4 initTetrad;
  initTetrad.setCol(0,e0);  
  initTetrad.setCol(1,e1);
  initTetrad.setCol(2,e2);
  initTetrad.setCol(3,e3);
  initTetrad.invert();


  /* ----------------------------------
   *  Set FermiWalker integrator.  
   * ---------------------------------- */     
  m4d::FermiWalker* motion = new m4d::FermiWalker( metric );
  motion->setAffineParamStep(0.01);
  motion->setCalcWithWorldline(true);
  motion->set_x_tau( x_tau );
  motion->set_u_tau( u_tau );
    
  mystruct par = {2.0*mass, omega, R};
  motion->set_params( &par );
  
  
  /* ----------------------------------
   *  Set initial position, direction, and local tetrad.
   * ---------------------------------- */     
  motion->setInitialPosition(initPos);
  motion->setInitialDirection(e0);
  motion->setInitialTetrad(e0,e1,e2,e3);
  motion->get_a_tau(0.0,a);
  //a.printS();

  /* ----------------------------------
   *  Test, whether local tetrad is orthonormal.
   * ---------------------------------- */   
#if 0
  double n;
  metric->calcProduct(initPos,e0,e0,n);
  fprintf(stderr,"%f\n",n);
  metric->calcProduct(initPos,e1,e1,n);
  fprintf(stderr,"%f\n",n);
  metric->calcProduct(initPos,e2,e2,n);
  fprintf(stderr,"%f\n",n);
  metric->calcProduct(initPos,e3,e3,n);
  fprintf(stderr,"%f\n",n);
#endif  
  
  
  /* ----------------------------------
   *  Calculate one orbit.
   * ---------------------------------- */    
  m4d::vec4 pos,ne;  
  int i=0; 
  while (i<maxNumPoints && pos[3]<=2.0*M_PI)
  {
    pos = motion->getPosition();
    e1  = motion->getE(1);
    fprintf(stderr,"\r%5d  %8.4f",i,motion->getAffineParam());
    fprintf(stdout,"%5d  %8.4f  %12.8f %12.8f %12.8f %12.8f  ",i,motion->getAffineParam(),pos[0],pos[1],pos[2],pos[3]);
    e1.printS(stdout,"%12.8f ");
    motion->nextStepWL();
    i++;
  }
  fprintf(stderr,"\n");
  
  // After one orbit, the base vector e1 has rotated:
  motion->getTetrad(e0,e1,e2,e3);
  ne = initTetrad*e1;  
  fprintf(stderr,"\nangle: %10.8f\n",atan2(ne[3],ne[1])*RAD_TO_DEG);
  
  delete motion;
  delete metric;
  return 1;
}

