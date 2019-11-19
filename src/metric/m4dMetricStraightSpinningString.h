// -------------------------------------------------------------------------------
/*
    m4dMetricStraightSpinningString.h

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

/*!  \class  m4d::MetricStraightSpinningString
     \brief  StraightSpinningString metric in cylindrical coordinates (t,rho,phi,z).

             The line element is given by
             \f[ds^2 = -\left(c\,dt-a\,d\varphi\right)^2+d\rho^2+k^2\rho^2\,d\varphi^2+dz^2\f]
             with constants a and k>0.

             The natural local tetrad is given by
             \f[ \mathbf{e}_{(0)} = \frac{1}{c}\partial_t,\quad \mathbf{e}_{(1)} = \partial_{\rho},\quad \mathbf{e}_{(2)} = \frac{1}{k\rho}\left(\frac{a}{c}\partial_t+\partial_{\varphi}\right),\quad \mathbf{e}_{(3)} = \partial_z. \f]

             From the Euler-Lagrange formalism, we have
             \f[ \dot{\rho}^2 + \frac{1}{k^2\rho^2}\left(h_2-\frac{ah_1}{c}\right)^2 - \kappa c^2 = \frac{h_1^2}{k^2} \f]
             with constants of motion \f$h_1=c(c\dot{t}-a\dot{\varphi})\f$ and \f$h_2=a(c\dot{t}-a\dot{\varphi})+k^2\rho^2\dot{\varphi}\f$.

             The metric is taken from<br>
             Volker Perlick, Living Reviews in Relativity, 7(9), 2004.<br>
             http://www.livingreviews.org/lrr-2004-9


*/
// -------------------------------------------------------------------------------

#ifndef M4D_METRIC_STRAIGHT_SPINNING_STRING_H
#define M4D_METRIC_STRAIGHT_SPINNING_STRING_H

#include "m4dMetric.h"

namespace m4d {

// ---------------------------------------------------
//    class definition:   MetricStraightSpinningString
// ---------------------------------------------------
class MetricStraightSpinningString : public Metric {
public:
    MetricStraightSpinningString(double a = 1.0, double k = 0.5);
    virtual ~MetricStraightSpinningString();

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

    //virtual bool   calcDerivs             ( const double y[], double dydx[] );
    //virtual double testConstraint ( const double y[], const double kappa );

    virtual bool   setParam(const char* pName, double val);

    virtual bool   effPotentialValue(const vec4 pos, const vec4 cdir, enum_geodesic_type type, const double x, double &val);
    virtual bool   totEnergy(const vec4 pos, const vec4 cdir, const double x, double &val);

    virtual bool   report(const vec4 pos, const vec4 cdir, char*&text);

// --------- protected methods -----------
protected:
    virtual void   setStandardValues();
    void   calcClosestApproach(const double* pos, const double* ldir);

// -------- protected attribute ---------
protected:
    double  mA;
    double  mK;
    double  rho_ca;
};

} // end namespace m4d

#endif
