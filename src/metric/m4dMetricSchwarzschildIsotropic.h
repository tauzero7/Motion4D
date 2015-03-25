// -------------------------------------------------------------------------------
/*
    m4dMetricSchwarzschildIsotropic.h

  Copyright (c) 2010-2014  Thomas Mueller


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

/*!  \class  m4d::MetricSchwarzschildIsotropic
     \brief  Schwarzschild metric in Cartesian isotropic coordinates (t,x,y,z).

             The line element is given by

             \f[ds^2 = -\left(\frac{1-\rho_s/\rho}{1+\rho_s/\rho}\right)^2 c^2dt^2 + \left(1+\frac{\rho_s}{\rho}\right)^4\left(dx^2+dy^2+dz^2\right),\f]
             where \f$\rho^2=x^2+y^2+z^2\f$ and \f$\rho_s=GM/(2c^2)\f$ is the Schwarzschild radius. G is Newton's
             constant, M is the mass of the black hole, and c is the speed of light.

             The natural local tetrad is given by
             \f[ \mathbf{e}_{(0)} = \frac{1+\rho_s/\rho}{1-\rho_s/\rho}\frac{\partial_t}{c},\quad \mathbf{e}_{(1)}=\left(1+\frac{\rho_s}{\rho}\right)^{-2}\partial_x,\quad \mathbf{e}_{(2)}=\left(1+\frac{\rho_s}{\rho}\right)^{-2}\partial_y,\quad \mathbf{e}_{(3)} = \left(1+\frac{\rho_s}{\rho}\right)^{-2}\partial_z.\f]

             Detailed discussions about the Schwarzschild metric can be found
             in the standard literature, e.g. \ref lit_wald "Wald", \ref lit_rindler "Rindler", \ref lit_mtw "MTW".

*/
// -------------------------------------------------------------------------------

#ifndef M4D_METRIC_SCHWARZSCHILD_ISOTROPIC_H
#define M4D_METRIC_SCHWARZSCHILD_ISOTROPIC_H

#include "m4dMetric.h"

#include <gsl/gsl_errno.h>

namespace m4d {

// ---------------------------------------------------
//    class definition:   MetricSchwarzschildIsotropic
// ---------------------------------------------------
class MetricSchwarzschildIsotropic : public Metric {
public:
    MetricSchwarzschildIsotropic(double mass = 1.0);
    virtual ~MetricSchwarzschildIsotropic();

// --------- public methods -----------
public:
    virtual bool   calculateMetric(const double* pos);
    virtual bool   calculateChristoffels(const double* pos);
    virtual bool   calculateChrisD(const double* pos);

    virtual void   localToCoord(const double* pos, const double* ldir, double* dir,
                                enum_nat_tetrad_type  type = enum_nat_tetrad_default);
    virtual void   coordToLocal(const double* pos, const double* cdir, double* ldir,
                                enum_nat_tetrad_type  type = enum_nat_tetrad_default);

    virtual bool   breakCondition(const double* pos);

    virtual double testConstraint(const double y[], const double kappa);

    virtual bool   calcProduct(const double* pos, const double* u, const double* v, double &prod, bool preCalcMetric = true);

    virtual bool   setParam(std::string pName, double val);

    virtual bool   report(const vec4 pos, const vec4 cdir, std::string &text);


    bool           calcBetaOfCircOrbit(const double* pos, double &beta);


// --------- protected methods -----------
protected:
    virtual void  setStandardValues();

    double  calc_rho(const double* pos);
    void    calc_drho(const double* pos, double &drdx, double &drdy, double &drdz);
    void    calc_d2rho(const double* pos, double &drdxdx, double &drdxdy, double &drdxdz,
                       double &drdydy, double &drdydz, double &drdzdz);

    void    calc_orbits();

// -------- protected attribute ---------
protected:
    //! Schwarzschild radius rs=2GM/c^2.
    double rho_s;
    double mMass;

    double rho_po;  // photon orbit;
    double rho_lso; // last stable timelike circular orbit

    double rho;
};

} // end namespace m4d

#endif
