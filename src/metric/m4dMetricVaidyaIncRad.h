// -------------------------------------------------------------------------------
/*
    m4dMetricVaidyaIncRad.h

  Copyright (c) 2016  Thomas Mueller


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

/*!  \class  m4d::MetricVaidyaIncRad
     \brief  Incoming radiation Vaidya metric in spherical coordinates (v,r,theta,phi).

             The line element is given by

             \f[ds^2 = 2dvdr - \left(1-\frac{2m(v)}{r}\right)dv^2 + r^2(d\theta^2+\sin^2\theta d\phi^2).\f]

             Detailed discussions about the Vaidya metric can be found in Griffiths,Podolsky.

*/
// -------------------------------------------------------------------------------

#ifndef M4D_METRIC_VAIDYA_INC_RAD_H
#define M4D_METRIC_VAIDYA_INC_RAD_H

#include "m4dMetric.h"

namespace m4d {

// ---------------------------------------------------
//    class definition:   MetricVaidyaIncRad
// ---------------------------------------------------
class MetricVaidyaIncRad : public Metric {
public:
    MetricVaidyaIncRad(double p = 1.0);
    virtual ~MetricVaidyaIncRad();

// --------- public methods -----------
public:
    virtual bool   calculateMetric(const double* pos);
    virtual bool   calculateChristoffels(const double* pos);
    
    virtual void   localToCoord(const double* pos, const double* ldir, double* dir,
                                enum_nat_tetrad_type  type = enum_nat_tetrad_default);
    virtual void   coordToLocal(const double* pos, const double* cdir, double* ldir,
                                enum_nat_tetrad_type  type = enum_nat_tetrad_default);

    virtual bool   breakCondition(const double* pos);

    virtual bool   report(const vec4 pos, const vec4 cdir, std::string &text);

    virtual bool   setParam(std::string pName, double val);

// --------- protected methods -----------
protected:
    virtual void setStandardValues();

    void calcMassFunc(const double v, double &m);    
    void calcMassFunc(const double v, double &m, double &dmdv);

// -------- protected attribute ---------
protected:
    double mP;  // constant factor
};

} // end namespace m4d

#endif // M4D_METRIC_VAIDYA_INC_RAD_H

