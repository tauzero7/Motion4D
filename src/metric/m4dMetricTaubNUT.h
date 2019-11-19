// -------------------------------------------------------------------------------
/*
    m4dMetricTaubNUT.h

  Copyright (c) 2009-2014  Thomas Mueller


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

/*!  \class  m4d::MetricTaubNUT
     \brief  TaubNUT metric in Boyer-Lindquist like spherical coordinates (t,r,theta,phi).

             The line element is given by

             \f[ds^2 = -\frac{\Delta}{\Sigma}\left(dt+2l\cos\vartheta d\varphi\right)^2 + \Sigma^2\left(\frac{dr^2}{\Delta}+d\vartheta^2 + \sin(\vartheta)^2 d\varphi^2\right),\f]

             where \f$\Delta=r^2-2Mr-l^2\f$ and \f$\Sigma=r^2+l^2\f$.

             The metric is taken from<br>
             Bini et al, <b>Circular holonomy in the Taub-NUT spacetime</b>,  Class. Quantum Grav .<b>19</b>, 5481 (2002).

*/
// -------------------------------------------------------------------------------

#ifndef M4D_METRIC_TAUB_NUT_H
#define M4D_METRIC_TAUB_NUT_H

#include "m4dMetric.h"

namespace m4d {

// ---------------------------------------------------
//    class definition:   MetricTaubNUT
// ---------------------------------------------------
class MetricTaubNUT : public Metric {
public:
    MetricTaubNUT(double mass = 1.0, double l = 0.5);
    virtual ~MetricTaubNUT();

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

    virtual bool   calcDerivs(const double y[], double dydx[]);
    virtual double testConstraint(const double y[], const double kappa);

    virtual bool   setParam(const char* pName, double val);

    virtual bool   effPotentialValue(const vec4 pos, const vec4 cdir, enum_geodesic_type type, const double x, double &val);
    virtual bool   totEnergy(const vec4 pos, const vec4 cdir, const double x, double &val);

    virtual bool   report(const vec4 pos, const vec4 cdir, char*&text);

// --------- protected methods -----------
protected:
    virtual void   setStandardValues();
    void   calcCriticalRadius();
    void   calcPhotonOrbit();
    void   calcFunctions(const double* pos, double &D, double &S);

// -------- protected attribute ---------
protected:
    double  mMass;
    double  mL;

    double  mCritRadius;
    double  mPhotonOrbit;
};

} // end namespace m4d

#endif
