// -------------------------------------------------------------------------------
/*
    m4dMetricHalilsoyWave.h

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

/*!  \class  m4d::MetricHalilsoyWave
     \brief  Gravitational wave solution given by Halilsoy.

       The Halilsoy standing gravitational wavein cylindrical coordinates (t,rho,phi,z) is given by
       \f[ ds^2 = e^{-2U}\left[e^{2k}\left(d\rho^2-dt^2\right) + \rho^2d\varphi^2\right] + e^{2U}\left(dz+A\,d\varphi\right)^2,\f]
       where \f$e^{-2U}=\cosh^2\alpha\,e^{-2CJ_0(\rho)\cos t}+\sinh^2\alpha\,e^{2CJ_0(\rho)\cos t}\f$, \f$A=-2C\sinh(2\alpha)\rho J_1(\rho)\sin t\f$,
       \f$k=\frac{1}{2}C^2\left[\rho^2\left(J_o(\rho)^2+J_1(\rho)^2\right)-2\rho J_0(\rho)J_1(\rho)\cos^2t\right]\f$, and the Bessel functions J_i.

     For further information see<br>
     M. Halilsoy, <b>Cross-Polarized Cylindrical Gravitational Waves of Einstein and Rosen</b>,<br>
     Il Nuovo Cimento <b>102 B</b>, 563 (1988).<br>
     or<br>
     Exact solutions (p.353).
 */
// -------------------------------------------------------------------------------

#ifndef M4D_METRIC_HALILSOY_WAVE_H
#define M4D_METRIC_HALILSOY_WAVE_H

#include "m4dMetric.h"
#include "gsl/gsl_sf_bessel.h"

namespace m4d {

// ---------------------------------------------------
//    class definition:   MetricHalilsoyWave
// ---------------------------------------------------
class MetricHalilsoyWave: public Metric {
public:
    MetricHalilsoyWave(double alpha = 0.0, double C = 1.0);
    ~MetricHalilsoyWave();


// --------- public methods --------------
public:
    virtual bool   calculateMetric(const double* pos);
    virtual bool   calculateChristoffels(const double* pos);
    virtual bool   calculateChrisD(const double* pos);

    virtual void   localToCoord(const double* pos, const double* ldir, double* dir,
                                enum_nat_tetrad_type type = enum_nat_tetrad_default);
    virtual void   coordToLocal(const double* pos, const double* cdir, double* ldir,
                                enum_nat_tetrad_type type = enum_nat_tetrad_default);

    virtual bool   breakCondition(const double* pos);

    virtual double testConstraint(const double y[], const double kappa);

    virtual bool   setParam(const char* pName, double val);

    virtual bool   report(const vec4 pos, const vec4 cdir, std::string &text);

// --------- protected methods -----------
protected:
    virtual void setStandardValues();

// ----- specific protected methods ------
    void  calcMetricFunc(const double* pos);
    void  calcDiffMetricFunc(const double* pos);

// -------- protected attribute ----------
protected:
    double mAlpha;
    double mC;

    // metric functions
    double fV, fA, fK;

    // derivatives of the metric functions
    double fV_t, fV_tt, fV_rho, fV_rhorho, fV_trho;
    double fA_t, fA_tt, fA_rho, fA_rhorho, fA_trho;
    double fK_t, fK_tt, fK_rho, fK_rhorho, fK_trho;
};

} // end namespace m4d

#endif
