// -------------------------------------------------------------------------------
/*
    m4dMetricErezRosenVar.h

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

/*!  \class  m4d::MetricErezRosenVar

*/
// -------------------------------------------------------------------------------

#ifndef M4D_METRIC_EREZ_ROSEN_VAR_H
#define M4D_METRIC_EREZ_ROSEN_VAR_H

#include "m4dMetric.h"

namespace m4d {

// ---------------------------------------------------
//    class definition:   MetricErezRosenVar
// ---------------------------------------------------
class MetricErezRosenVar : public Metric {
public:
    MetricErezRosenVar(double mass = 1.0, double q = 0.1);
    virtual ~MetricErezRosenVar();

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

    virtual bool   setParam(const char* pName, double val);

    virtual double getCircularVelocity(const double r, const enum_nat_tetrad_type  tedType = enum_nat_tetrad_default);
    virtual vec4   getCircularFourVel(const vec4 pos, const enum_nat_tetrad_type  tedType = enum_nat_tetrad_default);

    virtual bool   report(const vec4 pos, const vec4 cdir, std::string &text);

    // --------- protected methods -----------
protected:
    virtual void setStandardValues();

    bool calcPotentials(const double* pos, double &psi, double &g, double &delta);
    bool calcDiffPots(const double* pos, double &psi, double &g, double &delta,
                      double &dpsidr, double &dpsidtheta,
                      double &dgdr, double &dgdtheta,
                      double &dDdr, double &dDdtheta);

    // -------- protected attribute ---------
protected:
    double mMass;
    double mQ;

};

} // end namespace m4d

#endif
