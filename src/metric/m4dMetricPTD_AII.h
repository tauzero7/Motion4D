// --------------------------------------------------------------------------------
/*
    m4dMetricPTD_AII.h

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
/*!  \class  m4d::MetricPTD_AII
     \brief  \ref lit_stephani "Exact Solutions:" Petrov type D solutions - Case AII.

             The line element is given by

             \f[ds^2 = z^2\left( dr^2 + \sinh^2 r d\varphi^2 \right) + \frac{z}{b-z}dz^2-\frac{b-z}{z}dt^2. \f]

             The natural local tetrad is given by
             \f[ \mathbf{e}_{(t)} = \sqrt{{\frac {z}{b-z}}}\partial_t, \quad
                 \mathbf{e}_{(r)} = \frac{1}{z}\partial_r, \quad
                 \mathbf{e}_{(\varphi)} = \frac {1}{z\sinh{r}}\partial_{\varphi}, \quad
                 \mathbf{e}_{(z)} = \sqrt{{\frac {b-z}{z}}}\partial_{z}.\f]


*/
// --------------------------------------------------------------------------------
#ifndef M4DMETRICPTDAII_H
#define M4DMETRICPTDAII_H

#include "m4dMetric.h"

namespace m4d {

// ---------------------------------------------------
//    class definition:   MetricKerrBL
// ---------------------------------------------------
class MetricPTD_AII : public Metric
{
public:
    MetricPTD_AII(double b = 1.0);
    virtual ~MetricPTD_AII();

    // --------- public methods -----------
public:
    virtual bool calculateMetric(const double* pos);
    virtual bool calculateChristoffels(const double* pos);
    virtual bool calculateChrisD(const double* pos);

    virtual void localToCoord(
        const double* pos, const double* ldir, double* dir, enum_nat_tetrad_type type = enum_nat_tetrad_default);
    virtual void coordToLocal(
        const double* pos, const double* cdir, double* ldir, enum_nat_tetrad_type type = enum_nat_tetrad_default);

    virtual bool breakCondition(const double* pos);

    virtual double testConstraint(const double y[], const double kappa);

    virtual bool setParam(const char* pName, double val);

    virtual bool report(const vec4 pos, const vec4 cdir, char*& text);

    // --------- specific public methods ----------
public:
    // --------- protected methods -----------
protected:
    virtual void setStandardValues();

    // -------- protected attribute ---------
protected:
    double Par_b;
};

} // end namespace m4d

#endif // M4DMETRICPTDAII_H
