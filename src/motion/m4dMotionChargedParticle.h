// --------------------------------------------------------------------------------
/*
    m4dMotionChargedParticle.h

  Copyright (c) 2015 Thomas Mueller


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

/*!  \class  m4d::MotionChargedParticle
     \brief  Motion of a charged particle



 */
// --------------------------------------------------------------------------------

#ifndef M4D_MOTION_CHARGED_PARTICLE_H
#define M4D_MOTION_CHARGED_PARTICLE_H

#include <iostream>

#include "m4dMotion.h"

namespace m4d {

// ---------------------------------------------------
//    class definition:   MotionChargedParticle
// ---------------------------------------------------
class API_M4D_EXPORT MotionChargedParticle : public Motion {
public:
    MotionChargedParticle(Metric* metric);
    virtual ~MotionChargedParticle();

// --------- public methods -----------
public:

    bool  setInitialVelocity(double fm, double v, double theta, double phi, double q_over_m,
                             enum_nat_tetrad_type  type = enum_nat_tetrad_default);
    void  getInitialVelocity(double &v, double &theta, double &phi);


    virtual enum_break_condition  calculateMotion(const vec4 initPos, double fm, double v, double theta_v, double phi_v,
            const double q_over_m,
            const vec4 e0, const vec4 e1, const vec4 e2, const vec4 e3,
            const int maxNumPoints,
            std::vector<vec4> &points,
            std::vector<vec4> &base0, std::vector<vec4> &base1, std::vector<vec4> &base2, std::vector<vec4> &base3);


    enum_break_condition  nextStep();


// --------- protected methods -----------
protected:
    bool  calcDerivs(const double y[], double dydx[]);


// -------- protected attribute ---------
protected:

    //! Current four-velocity in coordinates.
    double mInitVel[4];

    double  mVel;
    double  mTheta;
    double  mPhi;

    double mQoverM;
};

} // end namespace m4d

#endif

