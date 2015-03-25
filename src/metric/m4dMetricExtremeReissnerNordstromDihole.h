// -------------------------------------------------------------------------------
/*
    m4dMetricExtremeReissnerNordstromDihole.h

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

/*!  \class  m4d::MetricExtremeReissnerNordstromDihole
     \brief  Extreme Reissner-Nordstrom Dihole in Cartesian coordinates (t,x,y,z).

             The line element is given by

             \f[ds^2 = -\frac{dt}{U^2}+U^2\left(dx^2+dy^2+dz^2\right),\f]

             where \f$U=1+\frac{M_1}{r_1}+\frac{M_2}{r_2}\f$, \f$r_1=\sqrt{x^2+y^2+(z-1)^2}\f$, and \f$r_2=\sqrt{x^2+y^2+(z+1)^2}\f$.

             The metric is taken from<br>
             S. Chandrasekhar,<br>
             <b>The two-centre problem in general relativity: the scattering of radiation by two extreme Reissner-Nordstrom black-holes,</b><br>
             Proc. R. Soc. Lond. A <b>421</b>, 227-258 (1989).

*/
// -------------------------------------------------------------------------------

#ifndef M4D_METRIC_EXTREME_REISSNER_DIHOLE_H
#define M4D_METRIC_EXTREME_REISSNER_DIHOLE_H

#include "m4dMetric.h"
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_integration.h>

namespace m4d {

// ---------------------------------------------------
//    class definition:   MetricExtremeReissnerNordstromDihole
// ---------------------------------------------------
class MetricExtremeReissnerNordstromDihole : public Metric {
public:
    MetricExtremeReissnerNordstromDihole(double mass1 = 1.0, double mass2 = 1.0);
    virtual ~MetricExtremeReissnerNordstromDihole();

// --------- public methods -----------
public:
    virtual bool   calculateMetric(const double* pos);
    virtual bool   calculateChristoffels(const double* pos);
    virtual bool   calculateChrisD(const double* pos);
    virtual bool   calculateRiemann(const double* pos);

    virtual void   localToCoord(const double* pos, const double* ldir, double* dir,
                                enum_nat_tetrad_type  type = enum_nat_tetrad_default);
    virtual void   coordToLocal(const double* pos, const double* cdir, double* ldir,
                                enum_nat_tetrad_type  type = enum_nat_tetrad_default);

    virtual bool   breakCondition(const double* pos);

    virtual double testConstraint(const double y[], const double kappa);

    virtual bool   setParam(std::string pName, double val);

    virtual bool   transToTwoPlusOne(vec4 p, vec4 &cp);

    virtual bool   report(const vec4 pos, const vec4 cdir, std::string &text);

// --------- specific public methods ----------
public:
    bool  calcEscapeVelocityAlongZ(const vec4 pos, double &betaEscape, double m1, double m2);
    bool  calcEscapeVelocityAlongX(const vec4 pos, double &betaEscape, double MASS);
    bool  calcPointOfEquilibrium(double &pointOfEquilibrium);
    bool  calcASingularityPoint(double &singularityPoint);
    bool  calcVelocityToEquilibriumPoint(const vec4 pos, double &betaEquilibrium, double pointOfEquilibrium, double mass1, double mass2);
    bool  calcPeriodicAlongX(const vec4 pos, double &periodicAlongX, double MASS);
    bool  calcVelocityAtOriginAlongX(const vec4 pos, double &betaAtOriginAlongX, double MASS);
    bool  calcPhotOrbitsXY(double &r_Nm, double &r_Np, double mquer);
    bool  calcKsiCritXY(const vec4 pos, double &ksicritXY, double r_Np, double MASS);
    bool  calcRadiusIscoXY(double &radiusIscoXY, double MASS);
    bool  calcRadiusIucoXY(double &radiusIucoXY, double MASS);
    bool  calcBetaCoXY(double &betaCoXY, double radiusCoXY);
    bool  calcRadiusForOtherCo(double &radiusForOtherCo, double z);
    bool  calcRadiusForOtherCoForUnequalMasses(double &radiusForOtherCo, double z, double M1, double M2);
    bool  calcHeightOfOtherPhotonOrbits(double &heightOfOtherPhotonOrbits, double M);
    bool  calcVelocityForOtherTimelikeCo(const vec4 pos, double &betaTimelike, double M);
    bool  calcHeightOfOtherPhotonOrbitsForUnequalMasses(double &heightOfOtherPhotonOrbits1, double &heightOfOtherPhotonOrbits2, double M1, double M2);
    bool  calcVelocityForOtherTimelikeCoForUnequalMasses(const vec4 pos, double &betaTimelike, double M1, double M2);

// --------- protected methods -----------
protected:
    virtual void   setStandardValues();
    void   calc_r(const double y[]);
    double calc_U(const double y[]);
    void   calc_dU(const double y[], double &Ux, double &Uy, double &Uz);
    void   calc_ddU(const double y[], double &Uxx, double &Uyy, double &Uzz,
                    double &Uxy, double &Uxz, double &Uyz);

// -------- protected attribute ---------
protected:
    double  r1, r2;
    double  mM1;
    double  mM2;

    gsl_integration_workspace* w;
    gsl_function F;
};

// ------- for some numerical calculations --------
/*!           ( no class members ! )           */

struct inflectionPoints_params {
    double M;
};
double inflectionPoints(double x, void *params);
double inflectionPoints_deriv(double x, void *params);
void   inflectionPoints_fdf(double x, void *params, double *y, double *dy);

struct integrandForPeriodicAlongX_params {
    double M;
    double x0;
};
double integrandForPeriodicAlongX(double x, void *params);

struct heightOfOtherPhotonOrbits_params {
    double M;
};
double heightOfOtherPhotonOrbits(double z, void *params);
double heightOfOtherPhotonOrbits_deriv(double z, void *params);
void heightOfOtherPhotonOrbits_fdf(double z, void *params, double *y, double *dy);

struct heightOfOtherPhotonOrbitsForUnequalMasses_params {
    double M1, M2;
};
double heightOfOtherPhotonOrbitsForUnequalMasses(double z, void *params);

} // end namespace m4d

#endif
