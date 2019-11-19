// -------------------------------------------------------------------------------
/*
    m4dMetricChazyCurzonRot.h

  Copyright (c) 2011-2014  Thomas Mueller


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

/*!  \class  m4d::MetricChazyCurzonRot
     \brief  ChazyCurzonRot metric in Weyl coordinates (t,rho,phi,z).

             The line element is given by

             \f[ds^2 = e^{-2U}\left[e^{2k}\left(d\rho^2+dz^2\right)+\rho^2d\varphi^2\right]-e^{2U}\left(dt+A\,d\varphi\right)^2,\f]
             where \f$e^{-2U}=\cosh\frac{2m}{r}-p\sinh\frac{2m}{r}\f$, \f$2k=-m^2\rho^2/r^4\f$, \f$A=2qmz/r\f$, and \f$r^2=\rho^2+z^2\f$.
             The parameters \f$p,q\f$ are related via \f$p^2+q^2=1\f$.

             The natural local tetrad reads:


             See also
             Stephani, H., Kramer, D., MacCallum, M., Hoenselaers, C. and Herlt, E.,<br>
             Exact Solutions of the Einstein Field Equations<br>
             (Cambridge University Press, 2. edition, 2009)

*/
// -------------------------------------------------------------------------------

#ifndef M4D_METRIC_CHAZY_CURZON_ROT_H
#define M4D_METRIC_CHAZY_CURZON_ROT_H

#include "m4dMetric.h"

namespace m4d {

// ---------------------------------------------------
//    class definition:   MetricChazyCurzonRot
// ---------------------------------------------------
class MetricChazyCurzonRot : public Metric {
public:
    MetricChazyCurzonRot(double mass = 1.0, double p = 0.5);
    virtual ~MetricChazyCurzonRot();

// --------- public methods -----------
public:
    virtual bool   calculateMetric(const double* pos);
    virtual bool   calculateChristoffels(const double* pos);
    virtual bool   calculateChrisD(const double* pos);

    virtual void   localToCoord(const double* pos, const double* ldir, double* dir,
                                enum_nat_tetrad_type  type = enum_nat_tetrad_cylinder);
    virtual void   coordToLocal(const double* pos, const double* cdir, double* ldir,
                                enum_nat_tetrad_type  type = enum_nat_tetrad_cylinder);

    virtual bool   breakCondition(const double* pos);

    virtual double testConstraint(const double y[], const double kappa);

    virtual bool   setParam(const char* pName, double val);

    virtual bool   report(const vec4 pos, const vec4 cdir, char*&text);


// --------- protected methods -----------
protected:
    virtual void setStandardValues();

    void calcUkA(const double* pos, double &em2U, double &k, double &A);
    void calcDUka(const double* pos, double &dUdrho, double &dUdz,
                  double &dkdrho, double &dkdz,
                  double &dAdrho, double &dAdz);

// -------- protected attribute ---------
protected:
    double mMass;
    double m_p;
    double m_q;

};

} // end namespace m4d

#endif
