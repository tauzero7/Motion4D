// -------------------------------------------------------------------------------
/*
    m4dMetricMinkowski.h

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

/*!  \class  m4d::MetricMinkowski
     \brief  Minkowski metric in cartesian coordinates (t,x,y,z).

             The line element is given by
             \f[ ds^2 = -c^2dt^2 + dx^2 + dy^2 + dz^2. \f]

             The natural local tetrad is given by
             \f[ \mathbf{e}_{(0)} = \frac{1}{c}\partial_t,\quad \mathbf{e}_{(1)}=\partial_x,\quad
   \mathbf{e}_{(2)}=\partial_y,\quad \mathbf{e}_{(3)} = \partial_z.\f]

             Detailed discussions about the Minkowski metric can be found
             in the standard literature, e.g. \ref lit_wald "Wald", \ref lit_rindler "Rindler", \ref lit_mtw "MTW".
*/
// -------------------------------------------------------------------------------

#ifndef M4D_METRIC_MINKOWSKI_H
#define M4D_METRIC_MINKOWSKI_H

#include "m4dMetric.h"

namespace m4d {

// ---------------------------------------------------
//    class definition:   MetricMinkowski
// ---------------------------------------------------
class MetricMinkowski : public Metric
{
public:
    MetricMinkowski();
    virtual ~MetricMinkowski();

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

    virtual bool calcDerivs(const double y[], double dydx[]);
    virtual bool calcDerivsPar(const double y[], double dydx[]);

    virtual double testConstraint(const double y[], const double kappa);

    virtual bool transToTwoPlusOne(vec4 p, vec4& cp);

    virtual bool report(const vec4 pos, const vec4 cdir, char*& text);

    // --------- protected methods -----------
protected:
    virtual void setStandardValues();

    // -------- protected attribute ---------
protected:
};

} // end namespace m4d

#endif
