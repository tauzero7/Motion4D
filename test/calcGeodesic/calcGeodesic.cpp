// -------------------------------------------------------------------------------- 
/*
    calcGeodesic.cpp
 
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

/*! \file  calcGeodesic.cpp
    \brief Calculate a geodesic using a parameter file.

     Usage:  \verbatim ./m4dCalcGeodesic  <file.ini> \endverbatim

     The parameter file is handled by the m4d::Object class. 
     
     Output:  i  pos[0] pos[1] pos[2] pos[3]

     where
     <ul>
       <li>i      : number of point
       <li>pos[n] : n-th component of the point in pseudo-cartesian coordinates
     </ul>         
*/
// -------------------------------------------------------------------------------- 

#include <iostream>
#include <cmath>

#include <m4dGlobalDefs.h>
#include <metric/m4dMetricDatabase.h>
#include <motion/m4dGeodesic.h>
#include <motion/m4dGeodesicGSL.h>
#include <math/TransCoordinates.h>
#include <extra/m4dObject.h>
#include <extra/m4dUtilities.h>

/* ----------------------------------------------------------------- 
 *                     m a i n
 * ----------------------------------------------------------------- */
int main( int argc, char* argv[] )
{
  if (argc!=2)
  {
    fprintf(stderr,"Usage: ./m4dCalcGeodesic <file.ini>\n");
    return -1;
  }

  fprintf(stderr,"\n\n=================== Calculate Geodesic ===================\n\n");   
  
  /* -----------------------------------------
   *    Load parameter setting from file.
   * ----------------------------------------- */
  m4d::Object   mObject;
  if (!mObject.loadSettings(argv[1]))
  {
    fprintf(stderr,"Can't load settings file %s  or settings file is invalid!\n",argv[1]);
    return -1;
  }
  mObject.printSettings();
  
  /* -----------------------------------------
   *    Set geodesic integrator.
   * ----------------------------------------- */
  const gsl_odeiv_step_type*  step_type = NULL;
  switch (mObject.geodSolverType)
  {
    default:
    case m4d::gsIrk4:
      break;
    case m4d::gsIgslrk2:
      step_type = gsl_odeiv_step_rk2;
      break;
    case m4d::gsIgslrk4:
      step_type = gsl_odeiv_step_rk4;
      break;
    case m4d::gsIgslfehlberg:
      step_type = gsl_odeiv_step_rkf45;
      break;
    case m4d::gsIgslcash:
      step_type = gsl_odeiv_step_rkck;
      break;
    case m4d::gsIgslprinc:
      step_type = gsl_odeiv_step_rk8pd;
      break;
    case m4d::gsIi2:
      step_type = gsl_odeiv_step_rk2imp;
      break;
      /*
    case m4d::gsIbs:
      step_type = gsl_odeiv_step_bsimp;
      break;
      */
    case m4d::gsIm1:
      step_type = gsl_odeiv_step_gear1;
      break;
    case m4d::gsIm2:
      step_type = gsl_odeiv_step_gear2;
      break;
  }

  m4d::Geodesic*  geod;
  if (mObject.geodSolverType == m4d::gsIrk4)
    geod = new m4d::GeodesicRK4( mObject.currMetric, mObject.type );
  else 
  {
    geod = new m4d::GeodesicGSL( mObject.currMetric, step_type, mObject.type );
  }
  geod->setEpsilons( mObject.epsAbs, mObject.epsRel);
  geod->setInitialPosition(mObject.startPos);
  geod->setStepSizeControlled(mObject.stepsizeControlled);
  geod->setBoundingBox(m4d::vec4(-DBL_MAX,-DBL_MAX,-DBL_MAX,-DBL_MAX),m4d::vec4(DBL_MAX,DBL_MAX,DBL_MAX,DBL_MAX));


  /* -----------------------------------------
   *    Set local initial direction.
   * ----------------------------------------- */
  double eta = 1.0;
  double y0  = 1.0;
    
  if (mObject.type==m4d::enum_geodesic_timelike && fabs(mObject.vel)<1.0)
  {      
    y0  = 1.0/sqrt(1.0-mObject.vel*mObject.vel);
    eta = mObject.vel*y0;
  }
  
  double chi = mObject.chi*DEG_TO_RAD;
  double ksi = mObject.ksi*DEG_TO_RAD;
  m4d::vec4 locDir = M4D_SIGN(mObject.timeDirection)*y0*mObject.base[0]
                      + eta*sin(chi)*cos(ksi)*mObject.base[1]
                      + eta*sin(chi)*sin(ksi)*mObject.base[2]
                      + eta*cos(chi)*mObject.base[3];
  m4d::vec4 coDir;
  mObject.currMetric->localToCoord(mObject.startPos,locDir,coDir,mObject.tetradType);  
  geod->setInitialDirection(coDir);
  locDir.printS();
  coDir.printS();
  mObject.currMetric->coordToLocal(mObject.startPos,coDir,locDir,mObject.tetradType);
  locDir.printS();
  
  /* -----------------------------------------
   *    Calculate the geodesic.
   * ----------------------------------------- */  
  std::cerr << std::endl << "Integration stops: " << m4d::stl_break_condition[geod->calculateGeodesic(mObject.startPos,coDir,mObject.maxNumPoints,mObject.points,mObject.dirs,mObject.lambda)] << std::endl;
  std::cerr << "#points: " << mObject.points.size() << std::endl << std::endl;
  if (mObject.points.size()<=0)
    return -1;
  
  /* -----------------------------------------
   *    Transform geodesic points into
   *    pseudo-cartesian representation
   *    and print them to standard output.
   * ----------------------------------------- */
  unsigned int i = 0;
  bool ok = true;
  m4d::vec4 oldPos, newPos;
  while (i<mObject.points.size() && ok)
  {
    oldPos = m4d::vec4( mObject.points[i][0],mObject.points[i][1],mObject.points[i][2],mObject.points[i][3] );
    //m4d::TransCoordinates::toCartesianCoord(m4d::enum_coordinate_cartesian,oldPos,newPos);
    m4d::TransCoordinates::toCartesianCoord(mObject.currMetric->getCoordType(),oldPos,newPos);
    
    fprintf(stdout,"%5d  %12.8f %12.8f %12.8f %12.8f\n",i,newPos[0],newPos[1],newPos[2],newPos[3]);
    i++;
  }
  fprintf(stdout,"\n");
  
  delete geod;
  return 1;
}

