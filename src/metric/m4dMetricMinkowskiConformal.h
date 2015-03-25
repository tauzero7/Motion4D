// -------------------------------------------------------------------------------
/*
    m4dMetricMinkowskiConformal.h

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

/*!  \class  m4d::MetricMinkowskiConformal
     \brief  Minkowski metric in cartesian coordinates \f$(\psi,\xi,\vartheta,\varphi)\f$.

             The line element is given by
             \f[ ds^2 = -d\psi^2 + d\xi^2 + \sin^2\xi\left(d\vartheta^2+\sin^2\vartheta d\varphi^2\right). \f]

             The natural local tetrad is given by
             \f[ \mathbf{e}_{(\psi)} = \partial_{\psi},\quad \mathbf{e}_{(\xi)}=\partial_{\xi},\quad \mathbf{e}_{(\vartheta)}=\frac{1}{\sin\xi}\partial_{\vartheta},\quad \mathbf{e}_{(\varphi)} = \frac{1}{\sin\xi\sin\vartheta}\partial_{\varphi}.\f]

             Detailed discussions about the conformal Minkowski metric can be found
             in the standard literature.
*/
// -------------------------------------------------------------------------------

#ifndef M4D_METRIC_MINKOWSKI_CONFORMAL_H
#define M4D_METRIC_MINKOWSKI_CONFORMAL_H

#include "m4dMetric.h"

namespace m4d {

// ---------------------------------------------------
//    class definition:   MetricMinkowskiConformal
// ---------------------------------------------------
class MetricMinkowskiConformal : public Metric {
public:
    MetricMinkowskiConformal();
    virtual ~MetricMinkowskiConformal();


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
    //virtual bool   calcDerivsPar          ( const double y[], double dydx[] );

    virtual double testConstraint(const double y[], const double kappa);

    virtual bool   transToTwoPlusOne(vec4 p, vec4 &cp);

    virtual bool   report(const vec4 pos, const vec4 cdir, std::string &text);

// --------- protected methods -----------
protected:
    virtual void setStandardValues();

// -------- protected attribute ---------
protected:
};

} // end namespace m4d

#endif
