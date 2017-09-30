// --------------------------------------------------------------------------------
/*
    m4dMetricGoedelScaled.h

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

/*!  \class  m4d::MetricGoedelScaled
 *   \brief  goedel metric in cylindrical coordinates
 *
 *
 *           source: (metric cylindrical) Endre Kajari, Uni Ulm, scaled coordinates
 *
 *           parameter:  goedel radius rG
 */
// --------------------------------------------------------------------------------

#ifndef M4D_METRIC_GOEDEL_SCALED_H
#define M4D_METRIC_GOEDEL_SCALED_H

#include "m4dMetric.h"

// #define m4dGoedelEps 1.0e-10

namespace m4d {

// ---------------------------------------------------
//    class definition:   MetricGoedelScaled
// ---------------------------------------------------
class MetricGoedelScaled : public Metric {
public:
    //! Standard constructor for the GoedelScaled metric.
    MetricGoedelScaled(double rG = 1.0, double zeta = 0.0);
    virtual ~MetricGoedelScaled();

// --------- public methods -----------
public:
    virtual bool calculateMetric(const double* pos);
    virtual bool calculateChristoffels(const double* pos);

    virtual void localToCoord(const double* pos, const double* ldir, double* dir,
                              enum_nat_tetrad_type  type = enum_nat_tetrad_default);
    virtual void coordToLocal(const double* pos, const double* cdir, double* ldir,
                              enum_nat_tetrad_type  type = enum_nat_tetrad_default);


    virtual bool breakCondition(const double* pos);

    virtual bool  transToTwoPlusOne(vec4 p, vec4 &cp);


    virtual bool setParam(const char* pName, double val);

    virtual bool   report(const vec4 pos, const vec4 cdir, std::string &text);

// --------- protected methods -----------
protected:
    virtual void setStandardValues();

// -------- protected attribute ---------
protected:
    double mRG;
    double mZeta;

    //helpers
//     double A,B,Gamma,Xi;
};

} // end namespace m4d

#endif
