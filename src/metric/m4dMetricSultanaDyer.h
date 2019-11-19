// -------------------------------------------------------------------------------
/*
    m4dMetricSultanaDyer.h

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

/*!  \class  m4d::MetricSultanaDyer
     \brief  SultanaDyer cosmological black hole in spherical coordinates (t,r,theta,phi).

      The line element is given by
      \f[ds^2 = t^4\left[\left(1-\frac{2m}{r}\right)dt^2-\frac{4m}{r}dt\,dr-\left(1+\frac{2m}{r}\right)dr^2-r^2\left(d\vartheta^2+\sin^2\vartheta\,d\varphi^2\right)\right].\f]

      The comoving natural local tetrad reads
      \f[ \mathbf{e}_{(0)} = \frac{1+2m/r}{t^2\sqrt{1+2m/r}}\partial_t - \frac{2m/r}{t^2\sqrt{1+2m/r}}\partial_r,\quad \mathbf{e}_{(1)} = \frac{1}{t^2\sqrt{1+2m/r}}\partial_r,\quad \mathbf{e}_{(2)} = \frac{1}{t^2r}\partial_{\vartheta},\quad \mathbf{e}_{(3)} = \frac{1}{t^2r\sin\vartheta}\partial_{\varphi}.\f]

      The static natural local tetrad reads
      \f[ \mathbf{e}_{(0)} = \frac{1}{t^2\sqrt{1-2m/r}},\quad \mathbf{e}_{(1)} = \frac{2m/r}{t^2\sqrt{1-2m/r}}\partial_t + \frac{\sqrt{1-2m/r}}{t^2}\partial_r, \quad \mathbf{e}_{(2)} = \frac{1}{t^2r}\partial_{\vartheta},\quad \mathbf{e}_{(3)} = \frac{1}{t^2r\sin\vartheta}\partial_{\varphi}.\f]

      For further information see...<br>
      J. Sultana and C.C. Dyer, "Cosmological black holes: A black hole in the Einstein-de Sitter universe", Gen. Relativ. Gravit. <b>37</b>, 1349--1370 (2005).
 */
// -------------------------------------------------------------------------------

#ifndef M4D_METRIC_SULTANA_DYER_COS_BH_H
#define M4D_METRIC_SULTANA_DYER_COS_BH_H

#include "m4dMetric.h"

namespace m4d {

// ---------------------------------------------------
//    class definition:   MetricSultanaDyer
// ---------------------------------------------------
class MetricSultanaDyer: public Metric {
public:
    MetricSultanaDyer(double mass = 1.0);
    ~MetricSultanaDyer();


// --------- public methods --------------
public:
    virtual bool   calculateMetric(const double* pos);
    virtual bool   calculateChristoffels(const double* pos);
    virtual bool   calculateChrisD(const double* pos);

    virtual void   localToCoord(const double* pos, const double* ldir, double* dir,
                                enum_nat_tetrad_type type = enum_nat_tetrad_default);
    virtual void   coordToLocal(const double* pos, const double* cdir, double* ldir,
                                enum_nat_tetrad_type type = enum_nat_tetrad_default);

    virtual bool   breakCondition(const double* pos);

    virtual double testConstraint(const double y[], const double kappa);

    virtual bool   setParam(const char* pName, double val);

    virtual bool   report(const vec4 pos, const vec4 cdir, char*&text);

// --------- protected methods -----------
protected:
    virtual void setStandardValues();


// -------- protected attribute ----------
protected:
    double mMass;

};

} // end namespace m4d

#endif
