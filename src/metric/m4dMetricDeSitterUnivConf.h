// -------------------------------------------------------------------------------
/*
    m4dMetricDeSitterUnivConf.h

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

/*!  \class  m4d::MetricDeSitterUnivConf
     \brief  DeSitter metric in conformal, cartesian coordinates (T,x,y,z).

             The line element is given by

             \f[ds^2 = \frac{c^2}{H^2T^2}\left[-dT^2+dx^2+dy^2+dz^2\right],\f]

             where H is the Hubble parameter and \f$c^2/(H^2T^2)\f$ is called the conformal factor.

             If the metric parameter p < 0.0, then the conformal factor is set to 1.
*/
// -------------------------------------------------------------------------------

#ifndef M4D_METRIC_DESITTER_UNIV_CONF_H
#define M4D_METRIC_DESITTER_UNIV_CONF_H

#include "m4dMetric.h"

namespace m4d {

// ---------------------------------------------------
//    class definition:   MetricDeSitterUnivConf
// ---------------------------------------------------
class MetricDeSitterUnivConf : public Metric {
public:
    MetricDeSitterUnivConf(double h = 0.1, double p = 1.0);
    virtual ~MetricDeSitterUnivConf();

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

    virtual bool   setParam(const char* pName, double val);

    virtual bool   report(const vec4 pos, const vec4 cdir, std::string &text);

// --------- protected methods -----------
protected:
    virtual void setStandardValues();

// -------- protected attribute ---------
protected:
    double  mHubble;
    double  mP;
};

} // end namespace m4d

#endif
