// --------------------------------------------------------------------------------
/*
    m4dMetricPTD_AI.h

  Copyright (c) 2010-2014  Thomas Mueller, Frank Grave, Felix Beslmeisl


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
/*!  \class  m4d::MetricPTD_AI
     \brief  \ref lit_stephani "Exact Solutions:" Petrov type D solutions - Case AI.

             This case equates the Schwarzschild metric in spherical coordinates (t,r,theta,phi).

             The line element is given by

             \f[ds^2 = r^2\left( d\vartheta^2 + \sin^2 \vartheta d\varphi^2 \right) + \frac{r}{r-b}dr^2-\frac{r-b}{r}dt^2. \f]

             The natural local tetrad is given by
             \f[ \mathbf{e}_{(t)} = \sqrt{{\frac {r}{r-b}}}\partial_t, \quad
                 \mathbf{e}_{(r)} = \sqrt{{\frac {r-b}{r}}}\partial_r, \quad
                 \mathbf{e}_{(\vartheta)} = \frac{1}{r}\partial_{\vartheta}, \quad
                 \mathbf{e}_{(\varphi)} = \frac {1}{r\sin{\vartheta}}\partial_{\varphi}.\f]

             From the Euler-Lagrange formalism, we have an effective Potential \f$\frac{1}{2}\dot{r}^2+\frac{1}{2}V_{\mbox{eff}}(r)=\frac{1}{2}C_0^2\f$

              \f[V_{\mbox{eff}}(r)=-C_0^2+K\frac{-b+r}{r^3}-\kappa c^2\frac{-b+r}{r}\f]

             with the following constants of motion:
             \f[  C_0^2 = \dot{t}^2 \frac{r-b}{r}  \f]
             \f[  K     = \dot{\varphi}^2 r^4 + \dot{\vartheta}^2 r^4 \sin^2{\vartheta}\f]
             \f[  -\kappa c^2 = -K\frac1{r^2}-\dot{r}^2 \frac{r}{r - b} +\dot{t}^2\frac{r - b}{r}.\f]


*/

// --------------------------------------------------------------------------------
#ifndef M4DMETRICPTDAI_H
#define M4DMETRICPTDAI_H

#include "m4dMetric.h"

namespace m4d {

// ---------------------------------------------------
//    class definition:   MetricKerrBL
// ---------------------------------------------------
class MetricPTD_AI : public Metric {
public:
    MetricPTD_AI(double b = 1.0);
    virtual ~MetricPTD_AI();

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

    virtual bool   effPotentialValue(const vec4 pos, const vec4 cdir, enum_geodesic_type type, const double x, double &val);
    virtual bool   totEnergy(const vec4 pos, const vec4 cdir, const double x, double &val);

    virtual bool   setParam(const char* pName, double val);

    virtual bool   report(const vec4 pos, const vec4 cdir, char*&text);


// --------- specific public methods ----------
public:


// --------- protected methods -----------
protected:
    virtual void setStandardValues();
    void  calcConstantsOfMotion(const vec4 pos, const vec4 cdir);


// -------- protected attribute ---------
protected:
    double Par_b;
    double K, C0, C2, m0;

};

} // end namespace m4d

#endif // M4DMETRICPTDAI_H
