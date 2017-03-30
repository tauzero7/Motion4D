// --------------------------------------------------------------------------------
/*
    m4dFermiWalker.h

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

/*!  \class  m4d::FermiWalker
     \brief  The local tetrad of an arbitrarily moving object has to be Fermi-Walker transported.

        The Fermi-Walker transport equation for the local tetrad \f$\mathbf{e}_{(i)}=e_{(i)}^{\mu}\partial_{\mu}\f$ is given by:
         \f[ \frac{de_{(i)}^{\mu}}{d\tau} = -\Gamma_{\alpha\beta}^{\mu}e_{(0)}^{\alpha}e_{(i)}^{\beta} + \frac{1}{c}\left(\eta_{(i)(k)}a^{(k)}e_{(0)}^{\mu}-\eta_{(0)(i)}a^{(k)}e_{(k)}^{\mu}\right) \f]


        There are two possibilities to calculate the Fermi-Walker transport of
        a local tetrad. On the one hand, the proper 3-acceleration with respect
        to the local tetrad can be given for each time step. On the other hand,
        the worldline can be defined as a function \f$x=x(\tau)\f$.

        Because the four-acceleration is always orthogonal to the four-velocity,
        the zero- or time- component of the proper acceleration with respect to
        the local tetrad where the e0- base vector is tangential to the four-velocity
        vanishes.

        The worldline function \f$x=x(\tau)\f$ has to take care of the
        parameter pointer!


        Note that e0 is tangential to four-velocity : y[8]..y[11] = y[4]..y[7].

 */
// --------------------------------------------------------------------------------

#ifndef M4D_MOTION_FERMIWALKER_H
#define M4D_MOTION_FERMIWALKER_H

#include <iostream>

#include "m4dMotion.h"

namespace m4d {

// ---------------------------------------------------
//    class definition:   FermiWalker
// ---------------------------------------------------
class API_M4D_EXPORT FermiWalker : public Motion {
public:
    FermiWalker(Metric* metric);
    virtual ~FermiWalker();

// --------- public methods -----------
public:
    void  setCurrPropAccel(double a1, double a2, double a3);
    void  getCurrPropAccel(double &a1, double &a2, double &a3);

    bool  setInitialVelocity(double fm, double v, double theta, double phi,
                             enum_nat_tetrad_type  type = enum_nat_tetrad_default);
    void  getInitialVelocity(double &v, double &theta, double &phi);

    //! Reset proper time.
    void  resetProperTime() {
        mLambda = 0.0;
    }

    virtual enum_break_condition  calculateMotion(const vec4 initPos, double fm, double v, double theta_v, double phi_v,
            double a, double theta_a, double phi_a,
            const vec4 e0, const vec4 e1, const vec4 e2, const vec4 e3,
            const int maxNumPoints,
            std::vector<vec4> &points,
            std::vector<vec4> &base0, std::vector<vec4> &base1, std::vector<vec4> &base2, std::vector<vec4> &base3);


    void  setCalcWithWorldline(bool withworldline);

    void  set_x_tau(vec4(*func)(double, void*));
    void  set_u_tau(vec4(*func)(double, void*));
    void  set_a_tau(vec4(*func)(double, void*));
    void  set_params(void* params);

    bool  get_x_tau(double tau, vec4 &x);
    bool  get_u_tau(double tau, vec4 &u);
    bool  get_a_tau(double tau, vec4 &a);

    bool  initWorldline(double tauStart);
    bool  updateWorldline(double tau);

    enum_break_condition  nextStep();
    enum_break_condition  nextStepWL();


// --------- protected methods -----------
protected:
    bool  calcDerivs(const double y[], double dydx[]);
    bool  calcDerivsWL(const double y[], double dydx[]);


// -------- protected attribute ---------
protected:
    //! Acceleration with respect to local tetrad.
    double mPropAcc[4];
    //! Current four-velocity in coordinates.
    double mInitVel[4];
    //! Current four-acceleration.
    vec4 mCurrAcc;

    double  mVel;
    double  mTheta;
    double  mPhi;

    bool     mCalcWithWorldline;
    vec4(*x_tau)(double, void*);
    vec4(*u_tau)(double, void*);
    vec4(*a_tau)(double, void*);
    void*   mParams;
};

} // end namespace m4d

#endif

