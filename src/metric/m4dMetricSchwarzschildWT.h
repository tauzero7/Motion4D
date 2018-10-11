// -------------------------------------------------------------------------------
/*
    m4dMetricSchwarzschild.h

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

/*!  \class  m4d::MetricSchwarzschildWT
     \brief  Schwarzschild metric in spherical coordinates (t,r,theta,phi).

             The line element is given by

             \f[ds^2 = -\left(1-\frac{r_s}{r}\right) c^2dt^2 + \frac{dr^2}{1-r_s/r} + r^2\left(d\vartheta^2 + \sin^2\vartheta\,d\varphi^2\right),\f]
             where \f$r_s=2GM/c^2\f$ is the Schwarzschild radius. G is Newton's
             constant, M is the mass of the black hole, and c is the speed of light.

             The natural local tetrad is given by
             \f[ \mathbf{e}_{(0)} = \frac{1}{c\,\sqrt{1-r_s/r}}\partial_t,\quad \mathbf{e}_{(1)}=\sqrt{1-\frac{r_s}{r}}\partial_r,\quad \mathbf{e}_{(2)}=\frac{1}{r}\partial_{\vartheta},\quad \mathbf{e}_{(3)} = \frac{1}{r\sin\vartheta}\partial_{\varphi}.\f]


*/
// -------------------------------------------------------------------------------

#ifndef M4D_METRIC_SCHWARZSCHILD_WT_H
#define M4D_METRIC_SCHWARZSCHILD_WT_H

#include "m4dMetric.h"

namespace m4d {

// ---------------------------------------------------
//    class definition:   MetricSchwarzschildWT
// ---------------------------------------------------
class MetricSchwarzschildWT : public Metric {
public:
    MetricSchwarzschildWT(double mass = 1.0);
    virtual ~MetricSchwarzschildWT();

    // --------- public methods -----------
public:
    /**
     * @brief Calculate the contravariant metric components at position 'pos'.
     * @param pos  Position in metric coordinates where coefficients have to be evaluated.
     * @return true if calculation succeeded.
     */
    virtual bool   calculateMetric(const double* pos);
    virtual bool   calculateChristoffels(const double* pos);

    virtual void   localToCoord(const double* pos, const double* ldir, double* dir,
                                enum_nat_tetrad_type  type = enum_nat_tetrad_default);
    virtual void   coordToLocal(const double* pos, const double* cdir, double* ldir,
                                enum_nat_tetrad_type  type = enum_nat_tetrad_default);

    virtual bool   breakCondition(const double* pos);

//    virtual bool   calcDerivs(const double y[], double dydx[]);

    virtual double testConstraint(const double y[], const double kappa);

    virtual bool   setParam(const char* pName, double val);

    virtual bool   report(const vec4 pos, const vec4 cdir, std::string &text);


   
    // --------- protected methods -----------
protected:
    virtual void  setStandardValues();
    
    // -------- protected attribute ---------
protected:
    //! Schwarzschild radius rs=2GM/c^2.
    double rs;
    double mMass;
};

} // end namespace m4d

#endif
