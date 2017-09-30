// -------------------------------------------------------------------------------
/*
    m4dMetricKastorTraschen.h

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

/*!  \class  m4d::MetricKastorTraschen
     \brief  Kastor-Traschen metric in cartesian coordinates (t,x,y,z).

             The line element is given by


*/
// -------------------------------------------------------------------------------

#ifndef M4D_METRIC_KASTOR_TRASCHEN_H
#define M4D_METRIC_KASTOR_TRASCHEN_H

#include "m4dMetric.h"

namespace m4d {

/**
 * @brief The MetricKastorTraschen class
 */
class MetricKastorTraschen : public Metric {
public:
    MetricKastorTraschen(double H = 0.0);
    virtual ~MetricKastorTraschen();


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

    virtual bool   report(const vec4 pos, const vec4 cdir, std::string &text);


public:
    void  calcPotentials(const double* pos, double &Omega, double &a);
    void  calcPotDiffs(const double* pos, double &Omega, double &a,
                       double &dOdx, double &dOdy, double &dOdz,
                       double &dOdt, double &dadt);

protected:
    virtual void setStandardValues();


protected:
    double m1, z1; // mass and position of first black hole
    double m2, z2; // mass and position of second black hole
    double mH;
};

} // end namespace m4d

#endif

