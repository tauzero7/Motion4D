// -------------------------------------------------------------------------------
/*
    m4dMetricSchwarzschildCart.h

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

/*!  \class  m4d::MetricSchwarzschildCart
     \brief  Schwarzschild metric in cartesian coordinates (t,x,y,z).

             The line element is given by

             \f[ ds^2 = -\left(1-\frac{r_s}{r}\right)c^2dt^2 + \left(\frac{x^2}{1-r_s/r}+y^2+z^2\right)\frac{dx^2}{r^2} + \left(x^2+\frac{y^2}{1-r_s/r}+z^2\right)\frac{dy^2}{r^2} + \left(x^2+y^2+\frac{z^2}{1-r_s/r}\right)\frac{dz^2}{r^2} + \frac{2r_s}{r^2(r-r_s)}\left(xy\,dx\,dy+xz\,dx\,dz+yz\,dy\,dz\right), \f]

             with Schwarzschild radius \f$r_s = 2GM/c^2\f$.
*/
// -------------------------------------------------------------------------------

#ifndef M4D_METRIC_SCHWARZSCHILD_CARTNEW_H
#define M4D_METRIC_SCHWARZSCHILD_CARTNEW_H

#include "m4dMetric.h"

namespace m4d {

// ---------------------------------------------------
//    class definition:   MetricSchwarzschildCart
// ---------------------------------------------------
class MetricSchwarzschildCartNew : public Metric {
public:
    //! Standard constructor for the Schwarzschildcart.
    MetricSchwarzschildCartNew(double mass = 1.0);
    virtual ~MetricSchwarzschildCartNew();

// --------- public methods -----------
public:
    virtual bool   calculateMetric(const double* pos);
    virtual bool   calculateChristoffels(const double* pos);

    virtual void   localToCoord(const double* pos, const double* ldir, double* dir,
                                enum_nat_tetrad_type  type = enum_nat_tetrad_default);
    virtual void   coordToLocal(const double* pos, const double* cdir, double* ldir,
                                enum_nat_tetrad_type  type = enum_nat_tetrad_default);


    virtual bool   breakCondition(const double* pos);

    virtual double testConstraint(const double y[], const double kappa);


    virtual bool   setParam(const char* pName, double val);

    virtual bool   report(const vec4 pos, const vec4 cdir, char*&text);


// --------- protected methods -----------
protected:
    virtual void setStandardValues();

    void   calcLTcoeffs(const double* pos);

// -------- protected attribute ---------
protected:
    double rs;

    void coutChristoffel();

    // Tetrad coefficients;
    double A, B, C, D, E, F;
};

} // end namespace m4d

#endif


