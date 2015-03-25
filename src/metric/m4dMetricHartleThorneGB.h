// -------------------------------------------------------------------------------
/*
    m4dMetricHartleThorneGB.h

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

/*!  \class  m4d::MetricHartleThorneGB
     \brief  Hartle-Thorne metric in spherical coordinates (t,r,theta,phi).

             The line element is given by

             \f[ds^2 = -dt^2 + .\f]


             Hartle-Thorne metric following

             Kostas Glampedakis and Stanislav Babak,
             "Mapping spacetimes with LISA: inspiral of a test body in a ‘quasi-Kerr’ field",
             Class. Quantum Grav. 23 (2006) 4167–4188

             see also ApJ 753,175 (2012)
*/
// -------------------------------------------------------------------------------

#ifndef M4D_METRIC_HARTLE_THORNE_GB_H
#define M4D_METRIC_HARTLE_THORNE_GB_H

#include "m4dMetric.h"

namespace m4d {

// ---------------------------------------------------
//    class definition:   MetricHartleThorneGB
// ---------------------------------------------------
class MetricHartleThorneGB : public Metric {
public:
    MetricHartleThorneGB(double mass = 1.0, double angmom = 0.0, double eta = 0.0);
    virtual ~MetricHartleThorneGB();

public:
    virtual bool   calculateMetric(const double* pos);
    virtual bool   calculateChristoffels(const double* pos);
    virtual bool   calculateChrisD(const double* pos);

    virtual void   localToCoord(const double* pos, const double* ldir, double* dir,
                                enum_nat_tetrad_type  type = enum_nat_tetrad_default);
    virtual void   coordToLocal(const double* pos, const double* cdir, double* ldir,
                                enum_nat_tetrad_type  type = enum_nat_tetrad_default);

    virtual bool   breakCondition(const double* pos);

    virtual bool   setParam(std::string pName, double val);

    virtual bool   report(const vec4 pos, const vec4 cdir, std::string &text);

public:
    void  calcFunc(const double* pos);
    void  calcFuncDiff(const double* pos);

protected:
    virtual void setStandardValues();

protected:
    double  mMass, mAngmom, mEta;
    double sigma, delta;
    double dsigmadr, dsigmadth, deltadr;
    double htt, hrr, hthth, hphph;
    double dhttdr, dhttdth, dhrrdr, dhrrdth;
    double dhaadr, dhaadth, dhppdr, dhppdth;
};

} // end namespace m4d

#endif

