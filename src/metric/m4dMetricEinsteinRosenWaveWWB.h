// -------------------------------------------------------------------------------
/*
    m4dMetricEinsteinRosenWaveWWB.h

  Copyright (c) 2013-2014  Thomas Mueller


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

/*!  \class  m4d::MetricEinsteinRosenWaveWWB
     \brief  Einstein-Rosen wave with a Weber-Wheeler-Bonnor pulse  in cylindrical coordinates (t,rho,phi,z).

     Detailed discussions can be found in<br><br>

     J.B. Griffiths and S. Micciche<br>
     "The Weber-Wheeler-Bonnor pulse and phase shifts in gravitational soliton interactions"<br>
     Physics Letters A 223, 37-42 (1997)

*/
// -------------------------------------------------------------------------------

#ifndef M4D_METRIC_EINSTEIN_ROSEN_WAVE_WWB_H
#define M4D_METRIC_EINSTEIN_ROSEN_WAVE_WWB_H

#include "m4dMetric.h"

namespace m4d {

class MetricEinsteinRosenWaveWWB : public Metric {
public:
    MetricEinsteinRosenWaveWWB(double c = 1.0, double a = 1.0);
    virtual ~MetricEinsteinRosenWaveWWB();


public:
    virtual bool   calculateMetric(const double* pos);
    virtual bool   calculateChristoffels(const double* pos);
    virtual bool   calculateChrisD(const double* pos);

    virtual void   localToCoord(const double* pos, const double* ldir, double* dir,
                                enum_nat_tetrad_type  type = enum_nat_tetrad_default);
    virtual void   coordToLocal(const double* pos, const double* cdir, double* ldir,
                                enum_nat_tetrad_type  type = enum_nat_tetrad_default);

    virtual bool   breakCondition(const double* pos);
    virtual double testConstraint(const double* y, const double kappa);

    virtual bool   setParam(std::string pName, double val);

    virtual bool   report(const vec4 pos, const vec4 cdir, std::string &text);



protected:
    virtual void setStandardValues();

    void  calcPotentials(const double *pos, double &gam, double &psi);
    void  calcDiffPoti(const double *pos, double &gam, double &psi,
                       double &gamt, double &gamr, double &psit, double &psir);

protected:
    double m_c;
    double m_a;
};

} // end namespace m4d

#endif // M4D_METRIC_EINSTEIN_ROSEN_WAVE_WWB_H
