// --------------------------------------------------------------------------------
/*
    m4dMetricPravda_C_Can.h

  Copyright (c) 2010-2014  Thomas Mueller, Frank Grave, Felix Beslmeisl


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

/*!  \class  m4d::MetricPravda_C_Can
     \brief  This is the C-Metric as given in class MetricPTD_C in "canonical coordinates adapted to the boost-rotation symmetry".

             The line element is given by

             \f[ds^2 = \frac{1}{z^2-t^2}\left(e^{\rho}r^2(z\,dt-t\,dz)^2 -e^{\lambda}(z\,dz-t\,dt)^2 \right) - e^{\lambda}\,dr^2-r^2e^{-\rho}\,d\varphi^2\f]
             with \f[e^{\rho}=\frac{R_3+R+Z_3-r^2}{4\alpha^2\left(R_1+R+Z_1-r^2\right)},\f]
                  \f[e^{\lambda}=\frac{2\alpha^2\left(R(R+R_1+Z_1)-Z_1r^2\right) \left(R_1R_3+(R+Z_1)(R+Z_3)-(Z_1+Z_3)r^2) \right)}{R_iR_3\left( R(R +R_3+Z_3) - Z_3r^2\right)},\f]
                  \f[R=\frac12\left(z^2-t^2+r^2\right),\f]
                  \f[R_i=\sqrt{(R + Z_i)^2 - 2Z_i r^2},\f]
                  \f[Z_i=z_i-z_2,\f]
                  \f[\alpha^2=\frac14 \frac{m^2}{A^6(z_2-z_1)^2(z3-z1)^2},\f]
                  \f[q=\frac1{4\alpha^2},\f]
             and \f$z_3<z_1<z2\f$ the roots of \f$2A^4z^3-A^2z^2+m^2.\f$
             A and m are real parameters.

             The natural local tetrad for \f$z^2-t^2>0\f$ is given by
             \f[  \mathbf{e}_{(t)} = \frac{1}{\sqrt{z^2-t^2}}\left( qze^{-\rho/2}\,\partial_t, + te^{-\lambda/2}\,\partial_z,\right)\quad
                  \mathbf{e}_{(r)} = e^{-\lambda/2}\,\partial_r, \quad
                  \mathbf{e}_{(\varphi)} = re^{\rho/2}\,\partial_{\varphi}, \quad
                  \mathbf{e}_{(z)} = \frac{1}{\sqrt{z^2-t^2}}\left( qte^{-\rho/2}\,\partial_t, + ze^{-\lambda/2}\,\partial_z,\right)\f]

             and for \f$z^2-t^2<0\f$ by
             \f[  \mathbf{e}_{(t)} = \frac{1}{\sqrt{z^2-t^2}}\left( qte^{-\rho/2}\,\partial_t, + ze^{-\lambda/2}\,\partial_z,\right) \quad
                  \mathbf{e}_{(r)} = e^{-\lambda/2}\,\partial_r, \quad
                  \mathbf{e}_{(\varphi)} = re^{\rho/2}\,\partial_{\varphi}, \quad
                  \mathbf{e}_{(z)} = \frac{1}{\sqrt{z^2-t^2}}\left( qze^{-\rho/2}\,\partial_t, + te^{-\lambda/2}\,\partial_z,\right).\f]

             This metric ist dicussed in<br>
             Pravda, V. and Pravdov&aacute; ,A., <b>Co-accelerated particles in the C-metric</b>  Class. Quantum Gravit. <b>18</b>, 1205 (2001).
*/

// --------------------------------------------------------------------------------
#ifndef M4DMETRICPRAVDA_C_CAN_H
#define M4DMETRICPRAVDA_C_CAN_H

#include "m4dMetric.h"

namespace m4d {

// ---------------------------------------------------
//    class definition:   MetricPravda_C_Can
// ---------------------------------------------------
class MetricPravda_C_Can : public Metric {
public:
    MetricPravda_C_Can(double A = 0.01,  double m = 1.0);
    virtual ~MetricPravda_C_Can();

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

    virtual void   calculateRoots(vec3 & roots, double p, double q);

    virtual bool   report(const vec4 pos, const vec4 cdir, std::string &text);


// --------- specific public methods ----------
public:

    double rho(double tau, double zeta, double eta);
    double rho_tau(double tau, double zeta, double eta);
    double rho_zeta(double tau, double zeta, double eta);
    double rho_eta(double tau, double zeta, double eta);

    double lambda(double tau, double zeta, double eta);
    double lambda_tau(double tau, double zeta, double eta);  //diff(lambda(tau,zeta,eta),tau)
    double lambda_zeta(double tau, double zeta, double eta); //diff(lambda(tau,zeta,eta),zeta)
    double lambda_eta(double tau, double zeta, double eta);  //diff(lambda(tau,zeta,eta),eta)

// --------- protected methods -----------
protected:
    virtual void setStandardValues();


// -------- protected attribute ---------
protected:
    double Par_A;
    double Par_m;
    vec3 z_i;   //z_i Roots of Polynom
    double Z1; //z_1-z_2
    double Z3; //z_3-2_2
    double alpha2;
    double q;

};

} // end namespace m4d

#endif //
