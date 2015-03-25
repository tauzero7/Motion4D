// -------------------------------------------------------------------------------
/*
    m4dMetricDeSitterUniv.h

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

/*!  \class  m4d::MetricDeSitterUniv
     \brief  DeSitter metric in cartesian coordinates (t,x,y,z).

             The line element is given by

             \f[ds^2 = -c^2dt^2+e^{2Ht}\left(dx^2+dy^2+dz^2\right),\f]

             where H is the Hubble parameter.

*/
// -------------------------------------------------------------------------------

#ifndef M4D_METRIC_DESITTER_UNIV_H
#define M4D_METRIC_DESITTER_UNIV_H

#include "m4dMetric.h"

namespace m4d {

// ---------------------------------------------------
//    class definition:   MetricDeSitterUniv
// ---------------------------------------------------
class MetricDeSitterUniv : public Metric {
public:
    MetricDeSitterUniv(double h = 0.1);
    virtual ~MetricDeSitterUniv();

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

    virtual double testConstraint(const double y[], const double kappa);

    virtual bool   transToTwoPlusOne(vec4 p, vec4 &cp);

    virtual bool   setParam(std::string pName, double val);

    virtual bool   report(const vec4 pos, const vec4 cdir, std::string &text);

// --------- protected methods -----------
protected:
    virtual void setStandardValues();

// -------- protected attribute ---------
protected:
    double  mHubble;
};

} // end namespace m4d

#endif
