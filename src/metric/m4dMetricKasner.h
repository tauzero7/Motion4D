// -------------------------------------------------------------------------------
/*
    m4dMetricKasner.h

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

/*!  \class  m4d::MetricKasner
     \brief  Kasner metric in cartesian coordinates (t,x,y,z).

             The line element is given by

             \f[ds^2 = -dt^2 + t^{2p_1}dx^2 + t^{2p_2}dy^2 + t^{2p_3}dz^2.\f]

             The parameters \f$ p_1,p_2,p_3\f$ can be represented by the Khalatnikov-Lifshitz parameter \f$u\f$:
             \f[ p_1 = -\frac{u}{1+u+u^2},\quad p_2 = \frac{1+u}{1+u+u^2}, \quad p_3 = \frac{u(1+u)}{1+u+u^2}.\f]

             Detailed discussions about the Kasner metric can be found in MTW.

*/
// -------------------------------------------------------------------------------

#ifndef M4D_METRIC_KASNER_H
#define M4D_METRIC_KASNER_H

#include "m4dMetric.h"

namespace m4d {

// ---------------------------------------------------
//    class definition:   MetricKasner
// ---------------------------------------------------
class MetricKasner : public Metric {
public:
    MetricKasner(double u = 0.0);
    virtual ~MetricKasner();

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

    virtual bool   setParam(const char* pName, double val);

    virtual bool   report(const vec4 pos, const vec4 cdir, char*&text);

// --------- specific public methods ----------
public:
    void  calc_parameters();

// --------- protected methods -----------
protected:
    virtual void setStandardValues();

// -------- protected attribute ---------
protected:
    double  p1, p2, p3, mU;
};

} // end namespace m4d

#endif

