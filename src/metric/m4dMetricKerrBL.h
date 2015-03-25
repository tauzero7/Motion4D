// --------------------------------------------------------------------------------
/*
    m4dMetricKerrBL.h

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

/*!  \class  m4d::MetricKerrBL
     \brief  Kerr metric in spherical Boyer-Lindquist coordinates (t,r,theta,phi).

             The line element is given by

            \f[ ds^2 = -\left(1-\frac{2mr}{\Sigma}\right) dt^2 - \frac{4mar \sin(\vartheta)^2}{\Sigma} dt\,d\varphi
                    + \frac{\Sigma}{\Delta} dr^2 + \Sigma d\vartheta^2
                    + \frac{r^2+a^2+(2ma^2r \sin(\vartheta)^2)}{\Sigma} \sin(\vartheta)^2 d\varphi^2,\f]

             with \f$\Delta = r^2-2mr+a^2\f$  and  \f$\Sigma = r^2+a^2 \cos(\vartheta)^2\f$. 'm' is
             the mass of the black hole and 'a' the angular momentum per unit rest mass.

             The natural local tetrad of the Kerr spacetime is given by
             \f[ \mathbf{e}_{(0)} = \Gamma\left(\partial_t+\zeta\partial_{\varphi}\right),\quad \mathbf{e}_{(1)}=\sqrt{\frac{\Delta}{\Sigma}}\partial_r,\quad \mathbf{e}_{(2)}=\frac{1}{\sqrt{\Sigma}}\partial_{\vartheta},\quad \mathbf{e}_{(3)} = \Gamma\left(\pm\frac{g_{t\varphi}+\zeta g_{\varphi\varphi}}{\sqrt{\Delta}\sin\vartheta}\partial_t\mp\frac{g_{tt}+\zeta g_{t\varphi}}{\sqrt{\Delta}\sin\vartheta}\partial_{\varphi}\right),\f]
             with the free parameter \f$\zeta\f$.
             In this metric, we use the following two entities:
             \f[ \omega = \frac{2mar}{\Sigma\left(r^2+a^2\right)+2ma^2r\sin^2\vartheta} \f]
             and
             \f[ \Gamma^{-2} = \left(1-\frac{2mr}{\Sigma}\right)+\frac{2mar\sin^2\vartheta}{\Sigma}\zeta-\left(r^2+a^2+\frac{2ma^2r\sin^2\vartheta}{\Sigma}\right)\zeta^2\sin^2\vartheta. \f]


             The event horizon is given by the radius \f$r_{+}=m+\sqrt{m^2-a^2}\f$, where \f$m\geq a\f$.

             The static limit (outer boundary of the ergosphere) is given by \f$r_0=m+\sqrt{m^2-a^2\cos^2\vartheta}\f$.

             This metric is discussed e.g. in <br><br>

             Bardeen et al,<br><b>"Rotating black holes: Locally nonrotating frames, energy extraction, and scalar synchrotron
             radiation,"</b><br>The Astrophysical Journal <b>178</b>, 347--359 (1972).<br>

             The default local tetrad is the locally nonrotating frame
             (enum_nat_tetrad_default = enum_nat_tetrad_lnrf).
*/
// --------------------------------------------------------------------------------
#ifndef M4D_METRIC_KERRBL_H
#define M4D_METRIC_KERRBL_H

#include "m4dMetric.h"

namespace m4d {

/**
 * @brief The MetricKerrBL class
 */
class MetricKerrBL : public Metric {
public:
    MetricKerrBL(double m = 1.0, double a = 0.5);
    virtual ~MetricKerrBL();

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

    virtual bool   calcDerivs(const double y[], double dydx[]);
    virtual double testConstraint(const double y[], const double kappa);

    virtual bool   setParam(std::string pName, double val);

    virtual bool   transToTwoPlusOne(vec4 p, vec4 &cp);

    virtual bool   effPotentialValue(const vec4 pos, const vec4 cdir, enum_geodesic_type type, const double x, double &val);
    virtual bool   totEnergy(const vec4 pos, const vec4 cdir, const double x, double &val);

    virtual double getCircularVelocity(const double r, const enum_nat_tetrad_type  tedType = enum_nat_tetrad_default);
    virtual vec4   getCircularFourVel(const vec4 pos, const enum_nat_tetrad_type  tedType = enum_nat_tetrad_default);

    virtual bool   report(const vec4 pos, const vec4 cdir, std::string &text);


    // --------- specific public methods ----------
public:
    bool  calcEventHorizon(double &rp);
    bool  calcErgosphere(double theta, double &r0);
    bool  calcPhotonOrbit(double &rphDirect, double &rphRetrograd);
    bool  calcMarginStabOrbit(double &rmsDirect, double &rmsRetrograd);
    bool  calcMSOvelocities(double &velDirect, double &velRetrograd);

    void  getCircularVelocities(const double r, double &b1, double &b2);

    // --------- protected methods -----------
protected:
    virtual void setStandardValues();
    virtual bool calcSigmaAndDelta(double r, double theta);
    virtual bool calcA(double r, double theta);
    virtual bool calcGamma(double r, double theta, double zeta);

    // -------- protected attribute ---------
protected:
    double mass;
    double angmom;

    double sigma;
    double delta;
    double A;
    double ga;

    double pm;
};

} // end namespace m4d

#endif
