// -------------------------------------------------------------------------------
/*
   m4dMetricPravda_C.cpp

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
// -------------------------------------------------------------------------------

#include "m4dMetricPravda_C_Can.h"

namespace m4d {

#define eps 1.0e-9


/*! Standard constructor for the metric.
*
* \param  m : mass of the black holes.
* \param  A : acceleration constant of the black holes.
*/
MetricPravda_C_Can::MetricPravda_C_Can(double A, double m) {
    mMetricName  = "Pravda_C-Metric_Canonical_Coords";
    setCoordType(enum_coordinate_cylinder);
    //setCoordType(enum_coordinate_cartesian);

    mPhysicalUnits = enum_physical_constants_geom;
    mSpeedOfLight = 1.0;
    mGravConstant = 1.0;

    mSign = -1.0;

    Par_A = A;
    Par_m = m;

    addParam("a", Par_A);
    addParam("m", Par_m);

    setStandardValues();
    setParam("a", Par_A);

    // mLocTeds.push_back(enum_nat_tetrad_static);

}

/*! Standard destructor for the metric.
 *
 */

MetricPravda_C_Can::~MetricPravda_C_Can() {

}


// *********************************** public methods ******************************
/*! Calculate the contravariant metric components at position 'pos'.
*
*  \param pos : pointer to position.
*/
bool
MetricPravda_C_Can::calculateMetric(const double* pos) {
    double tau = pos[0];
    double eta = pos[1];
    double zeta = pos[3];

    double  t1 = zeta * zeta;
    double  t2 = tau * tau;
    double  t4 = 1 / (t1 - t2);
    double  t5 = 2.0 * t1;
    double  t6 = 2.0 * t2;
    double  t7 = eta * eta;
    double  t8 = 2.0 * t7;
    double  t10 = Z3 * t7;
    double  t13 = sqrt(t5 - t6 + t8 + 4.0 * Z3 - 8.0 * t10);
    double  t14 = t13 / 2.0;
    double  t15 = t1 / 2.0;
    double  t16 = t2 / 2.0;
    double  t17 = t7 / 2.0;
    double  t18 = t14 + t15 - t16 - t17 + Z3;
    double  t20 = alpha2;
    double  t22 = t4 * t18 / t20;
    double  t24 = Z1 * t7;
    double  t27 = sqrt(t5 - t6 + t8 + 4.0 * Z1 - 8.0 * t24);
    double  t28 = t27 / 2.0;
    double  t29 = t28 + t15 - t16 - t17 + Z1;
    double  t31 = q * q;
    double  t33 = 1 / t29 / t31;
    double  t38 = t1 - t2 + t7;
    double  t41 = t38 * (t15 - t16 + t17 + t28 + Z1) / 2.0 - t24;
    double  t49 = t27 * t13 / 4.0 + (t15 - t16 + t17 + Z1) * (t15 - t16 + t17 + Z3) - (Z1 + Z3) * t7;
    double  t51 = t4 * t20 * t41 * t49;
    double  t54 = 1 / t27 / t13;
    double  t58 = 1 / (t38 * (t15 - t16 + t17 + t14 + Z3) / 2.0 - t10);
    double  t73 = -t22 * t33 * zeta * tau / 4.0 + 8.0 * t51 * t54 * t58 * zeta * tau;
    g_compts[0][0] = t22 * t33 * t1 / 4.0 - 8.0 * t51 * t54 * t58 * t2;
    g_compts[0][1] = 0.0;
    g_compts[0][2] = 0.0;
    g_compts[0][3] = t73;
    g_compts[1][0] = 0.0;
    g_compts[1][1] = -8.0 * t20 * t41 * t49 * t54 * t58;
    g_compts[1][2] = 0.0;
    g_compts[1][3] = 0.0;
    g_compts[2][0] = 0.0;
    g_compts[2][1] = 0.0;
    g_compts[2][2] = -4.0 * t7 / t18 * t20 * t29;
    g_compts[2][3] = 0.0;
    g_compts[3][0] = t73;
    g_compts[3][1] = 0.0;
    g_compts[3][2] = 0.0;
    g_compts[3][3] = t22 * t33 * t2 / 4.0 - 8.0 * t51 * t54 * t58 * t1;

    return true;
}

/*! Calculate the Christoffel symbols of the second kind at position 'pos'.
*
*  \param pos : pointer to position.
*/
bool
MetricPravda_C_Can::calculateChristoffels(const double* pos) {
    double tau = pos[0];
    double eta = pos[1];
    double zeta = pos[3];

    double  t1 = rho(tau, zeta, eta);
    double  t2 = lambda_tau(tau, zeta, eta);
    double  t3 = t1 * t2;
    double  t4 = q * q;
    double  t5 = tau * tau;
    double  t6 = t5 * t5;
    double  t7 = t4 * t6;
    double  t9 = t5 * tau;
    double  t10 = t9 * zeta;
    double  t11 = lambda(tau, zeta, eta);
    double  t12 = t4 * t4;
    double  t13 = t11 * t12;
    double  t14 = lambda_zeta(tau, zeta, eta);
    double  t15 = t13 * t14;
    double  t16 = t10 * t15;
    double  t17 = t1 * t14;
    double  t18 = t17 * t4;
    double  t19 = t10 * t18;
    double  t20 = t5 * t11;
    double  t21 = zeta * zeta;
    double  t24 = t20 * t21 * t12 * t2;
    double  t25 = t5 * t1;
    double  t26 = rho_tau(tau, zeta, eta);
    double  t28 = t25 * t26 * t21;
    double  t29 = t5 * t21;
    double  t30 = t4 * t11;
    double  t32 = t29 * t30 * t26;
    double  t34 = t2 * t4;
    double  t36 = t25 * t34 * t21;
    double  t38 = tau * t1;
    double  t41 = 2.0 * t38 * t30 * t21;
    double  t42 = t21 * zeta;
    double  t43 = tau * t42;
    double  t44 = rho_zeta(tau, zeta, eta);
    double  t46 = t43 * t30 * t44;
    double  t47 = t1 * t44;
    double  t48 = t43 * t47;
    double  t49 = t1 * t1;
    double  t52 = 2.0 * tau * t49 * t21;
    double  t53 = t21 * t21;
    double  t54 = t11 * t53;
    double  t55 = t4 * t26;
    double  t57 = -t3 * t7 - t16 + t19 - t24 - t28 + 2.0 * t32 + 2.0 * t36 + t41 + t46 - t48 - t52 - t54 * t55;
    double  t58 = 1 / t4;
    double  t60 = -t21 + t5;
    double  t61 = t60 * t60;
    double  t63 = 1 / t1;
    double  t65 = 1 / t11;
    double  t66 = 1 / t61 * t63 * t65;
    double  t69 = rho_eta(tau, zeta, eta);
    double  t71 = lambda_eta(tau, zeta, eta);
    double  t76 = 1 / t60;
    double  t77 = t76 * t58;
    double  t82 = t1 * zeta;
    double  t86 = zeta * t26;
    double  t88 = t30 * t86 * t9;
    double  t90 = zeta * t2;
    double  t93 = t4 * t5;
    double  t95 = t17 * t93 * t21;
    double  t98 = t30 * t44 * t21 * t5;
    double  t99 = t1 * t26;
    double  t103 = t3 * tau * t4 * t42;
    double  t115 = -t13 * t14 * t6 - t82 * t2 * t9 * t4 + 2.0 * t88 - t13 * t90 * t9 + t95 + t98 - t99 * t43 + 2.0 * t103 - t30 * t26 * tau * t42 + 2.0 * t1 * t11 * t42 * t4 - 2.0 * t49 * t42 - t47 * t53;
    double  t119 = t1 * t71;
    double  t121 = t21 * t11;
    double  t125 = t63 * t76;
    double  t127 = (-t119 * t5 + t121 * t69) * t65 * t125 / 2.0;
    double  t129 = t65 * t2 / 2.0;
    double  t134 = t65 * t63;
    double  t137 = (-t119 + t11 * t69) * zeta * tau * t134 * t76 / 2.0;
    double  t139 = t63 * t26 / 2.0;
    double  t142 = t82 * t26 * t9;
    double  t145 = t13 * t14 * t5 * t21;
    double  t147 = t11 * zeta * t4;
    double  t149 = 2.0 * t25 * t147;
    double  t150 = t47 * t29;
    double  t153 = 2.0 * t5 * t49 * zeta;
    double  t156 = t13 * t2 * tau * t42;
    double  t162 = (t17 * t7 - t142 + t88 - t145 + t149 - t150 - t153 + t103 - t156 + t30 * t44 * t53) * t58 * t66 / 2.0;
    double  t170 = t65 * zeta * tau * (t69 - t71 * t4) * t76 * t58 / 2.0;
    double  t173 = t53 * t1;
    double  t178 = (t30 * t26 * t6 - t16 + t19 - t28 - t24 + t41 + t46 - t48 - t52 + t173 * t34) * t58 * t66 / 2.0;
    double  t183 = t14 * zeta;
    double  t186 = tau * t11 * t4;
    double  t207 = 2.0 * t1 - eta * t69;
    double  t209 = t63 / eta * t207 / 2.0;
    double  t211 = t65 * t14 / 2.0;
    double  t218 = (t20 * t69 - t21 * t1 * t71) * t65 * t125 / 2.0;
    double  t219 = eta * eta;
    double  t222 = t44 * zeta;
    double  t230 = t65 / t49 / t1 * t76;
    double  t248 = t63 * t44 / 2.0;
    double  t256 = t9 * t44;
    double  t264 = -t99 * t6 - 2.0 * t9 * t49 + 2.0 * t19 + 2.0 * t9 * t1 * t30 - t256 * t147 - t256 * t82 + t36 + t32 - t43 * t15 - t43 * t18 + 2.0 * t46 - t54 * t12 * t2;
    double  t281 = -t30 * t44 * t6 - t142 + t88 - t145 + 2.0 * t95 - t153 + t149 + 2.0 * t98 - t150 + t103 - t156 - t173 * t14 * t4;
    christoffel[0][0][0] = -t57 * t58 * t66 / 2.0;
    christoffel[0][0][1] = t65 * (-t69 * t21 + t71 * t5 * t4) * t77 / 2.0;
    christoffel[0][0][2] = 0.0;
    christoffel[0][0][3] = -t115 * t58 * t66 / 2.0;
    christoffel[0][1][0] = -t127;
    christoffel[0][1][1] = t129;
    christoffel[0][1][2] = 0.0;
    christoffel[0][1][3] = -t137;
    christoffel[0][2][0] = 0.0;
    christoffel[0][2][1] = 0.0;
    christoffel[0][2][2] = -t139;
    christoffel[0][2][3] = 0.0;
    christoffel[0][3][0] = t162;
    christoffel[0][3][1] = t170;
    christoffel[0][3][2] = 0.0;
    christoffel[0][3][3] = t178;
    christoffel[1][0][0] = -t127;
    christoffel[1][0][1] = t129;
    christoffel[1][0][2] = 0.0;
    christoffel[1][0][3] = -t137;
    christoffel[1][1][0] = -(-t3 * t5 + t2 * t11 * t21 * t4 - t183 * t38 + t183 * t186) * t76 * t134 / 2.0;
    christoffel[1][1][1] = t65 * t71 / 2.0;
    christoffel[1][1][2] = 0.0;
    christoffel[1][1][3] = -(-t90 * t38 + t90 * t186 - t17 * t21 + t14 * t11 * t93) * t76 * t134 / 2.0;
    christoffel[1][2][0] = 0.0;
    christoffel[1][2][1] = 0.0;
    christoffel[1][2][2] = t209;
    christoffel[1][2][3] = 0.0;
    christoffel[1][3][0] = t137;
    christoffel[1][3][1] = t211;
    christoffel[1][3][2] = 0.0;
    christoffel[1][3][3] = t218;
    christoffel[2][0][0] = 0.0;
    christoffel[2][0][1] = 0.0;
    christoffel[2][0][2] = -t139;
    christoffel[2][0][3] = 0.0;
    christoffel[2][1][0] = 0.0;
    christoffel[2][1][1] = 0.0;
    christoffel[2][1][2] = t209;
    christoffel[2][1][3] = 0.0;
    christoffel[2][2][0] = t219 * (-t99 * t5 + t121 * t55 - t222 * t38 + t222 * t186) * t230 / 2.0;
    christoffel[2][2][1] = -t65 * eta * t207 / t49 / 2.0;
    christoffel[2][2][2] = 0.0;
    christoffel[2][2][3] = t219 * (-t86 * t38 + t86 * t186 - t47 * t21 + t44 * t11 * t93) * t230 / 2.0;
    christoffel[2][3][0] = 0.0;
    christoffel[2][3][1] = 0.0;
    christoffel[2][3][2] = -t248;
    christoffel[2][3][3] = 0.0;
    christoffel[3][0][0] = t162;
    christoffel[3][0][1] = t170;
    christoffel[3][0][2] = 0.0;
    christoffel[3][0][3] = t178;
    christoffel[3][1][0] = t137;
    christoffel[3][1][1] = t211;
    christoffel[3][1][2] = 0.0;
    christoffel[3][1][3] = t218;
    christoffel[3][2][0] = 0.0;
    christoffel[3][2][1] = 0.0;
    christoffel[3][2][2] = -t248;
    christoffel[3][2][3] = 0.0;
    christoffel[3][3][0] = -t264 * t58 * t66 / 2.0;
    christoffel[3][3][1] = -t65 * (t69 * t5 - t71 * t21 * t4) * t77 / 2.0;
    christoffel[3][3][2] = 0.0;
    christoffel[3][3][3] = -t281 * t58 * t66 / 2.0;

    return true;
}

/*! Calculate Jacobi matrix.
*   This function is not implemented in this metric.
*   Derivative terms become too extensive.
*
*  \param pos : pointer to position.
*/
bool
MetricPravda_C_Can::calculateChrisD(const double*) {
    // fprintf(stderr,"MetricPravda_C_Can::calculateChrisD ( const double* pos ) ... not implemented yet\n");
    return false;
}

/*! Transform local 4-direction to coordinate 4-direction.
*
*  \param  pos  :  pointer to position array.
*  \param  ldir :  pointer to local direction array.
*  \param  dir  :  pointer to calculated coordinate direction array.
*  \param  type :  type of tetrad.
*/
void
MetricPravda_C_Can::localToCoord(const double* pos, const double* ldir, double* dir,
                                 enum_nat_tetrad_type) {
    double tau = pos[0];
    double zeta = pos[3];
    double eta = pos[1];
    double lambda1 = sqrt(lambda(tau, zeta, eta));
    double rho1 = sqrt(rho(tau, zeta, eta));
    double k = 1 / sqrt(fabs(zeta * zeta - tau * tau));

    if (zeta * zeta - tau * tau > 0) {
        dir[0] = (ldir[0] * (q * zeta / (rho1)) + ldir[3] * (tau / (lambda1))) * k;
        dir[3] = (ldir[0] * (q * tau / (rho1)) + ldir[3] * (zeta / (lambda1))) * k;
    } else {
        dir[0] = (ldir[3] * (q * zeta / (rho1)) + ldir[0] * (tau / (lambda1))) * k;
        dir[3] = (ldir[3] * (q * tau / (rho1))   + ldir[0] * (zeta / (lambda1))) * k;
    }
    dir[1] = ldir[1] / (lambda1);
    dir[2] = ldir[2] / eta * (rho1);
}

/*! Transform coordinate 4-direction to local 4-direction.
*
*  \param  pos  :  pointer to position array.
*  \param  cdir :  pointer to coordinate direction.
*  \param  ldir :  pointer to calculated local direction array.
*  \param  type :  type of tetrad.
*/
void
MetricPravda_C_Can::coordToLocal(const double* pos, const double* cdir, double* ldir,
                                 enum_nat_tetrad_type) {
    double tau = pos[0];
    double zeta = pos[3];
    double eta = pos[1];
    double lambda1 = sqrt(lambda(tau, zeta, eta));
    double rho1 = sqrt(rho(tau, zeta, eta));
    double k = 1 / sqrt(fabs(zeta * zeta - tau * tau));

    if (zeta * zeta - tau * tau > 0) {
        ldir[0] = (cdir[0] * rho1 * zeta / q - cdir[3] * (rho1 * tau / q)) * k;
        ldir[3] = (cdir[0] * lambda1 * tau - cdir[3] * (lambda1 * zeta)) * k;
    } else {
        ldir[0] = (cdir[0] * lambda1 * tau - cdir[3] * (lambda1 * zeta)) * k;
        ldir[3] = (cdir[0] * rho1 * zeta / q - cdir[3] * (rho1 * tau / q)) * k;
    }
    ldir[1] = cdir[1] * lambda1 / k;
    ldir[2] = cdir[2] * eta / rho1 / k;


}


/*! Tests break condition
*  \param pos  :  position.
*
* Not implemented yet.
*/
bool
MetricPravda_C_Can::breakCondition(const double*) {
    bool br = false;
    return br;
}


/*! Tests whether the constraint equation is fulfilled.
*
*  The constraint equation for lightlike and timelike geodesics reads:
\verbatim
 sum = g_{\mu\nu} dot(x)^{\mu} dot(x)^{\nu} - kappa c^2 = 0.
\endverbatim
*  \param  y[]   : pointer to position and direction coordinates.
*  \param  kappa : timelike (-1.0), lightlike (0.0).
*  \return double : sum.
*/
double
MetricPravda_C_Can::testConstraint(const double y[], const double kappa) {
    calculateMetric(y);
    double sum = -mSign * kappa;
    for (int i = 0; i < 4; i++) {
        sum += g_compts[i][i] * y[4 + i] * y[4 + i];
    }
    sum += 2 * g_compts[0][3] * y[4] * y[4 + 3];
    return sum;
    //return 0;
}

/*! Set parameter 'pName' to 'val' and calculates the resulting constants \f$Z_i\f$, \f$\alpha^2\f$ and q.
*
*
*/
bool
MetricPravda_C_Can::setParam(const char* pName, double val) {
    Metric::setParam(pName, val);
    if (pName == "m") {
        Par_m = val;
    }
    if (pName == "a") {
        Par_A = val;
    }
    double A = Par_A;
    double m = Par_m;
    calculateRoots(z_i, -0.5 / A / A, 0.5 * m * m / A / A / A / A);
    Z1 = z_i[0] - z_i[1];
    Z3 = z_i[2] - z_i[1];
    alpha2 = 0.25 * m * m / A / A / A / A / A / A / (z_i[1] - z_i[0]) / (z_i[1] - z_i[0]) / (z_i[2] - z_i[0]) / (z_i[2] - z_i[0]);
    q = 0.25 / alpha2;
    return true;
}

/*! Calculates the roots of a polynom of type x^3 + ax^2 + c = 0.
 *
 *\param   a     : coefficient of the quadratic term.
 *\param   c     : y-axes offset
 *\param roots : reference Roots of the polynom sorted roots[2]<roots[0]<roots[1].
*/
void
MetricPravda_C_Can::calculateRoots(vec3 & roots, double a, double c) {

    double p = a * a / 3.0;
    double q = 2.0 * a * a * a / 27.0 + c;
    double u = acos(-q * 0.5 * sqrt(27.0 / p / p / p));

    double z1 =  sqrt(4.0 / 3.0 * p) * cos((u) / 3.0) - a / 3.0;
    double z2 = -sqrt(4.0 / 3.0 * p) * cos((u + M_PI) / 3.0) - a / 3.0;
    double z3 = -sqrt(4.0 / 3.0 * p) * cos((u - M_PI) / 3.0) - a / 3.0;

    double tmp;
    //sort
    if (z2 < z1) {
        tmp = z1;
        z1 = z2;
        z2 = tmp;
    }
    if (z3 < z2) {
        tmp = z3;
        z3 = z2;
        z2 = tmp;
    }
    if (z2 < z1) {
        tmp = z1;
        z1 = z2;
        z2 = tmp;
    }
    roots[0] = z2;
    roots[1] = z3;
    roots[2] = z1;
}

/*! Generate report.
*/
bool
MetricPravda_C_Can::report(const vec4 pos, const vec4 , std::string &text) {
    std::stringstream ss;
    ss << "Report for Pravda C metric\n\tcanonical coordinates : (t,r,phi,z)\n";
    ss << "---------------------------------------------------------------\n";
    ss << "Coordinate Ranges:" << std::endl;
    ss << "        t: arbitrary" << std::endl;
    ss << "        r: ( 0 < r )" << std::endl;
    ss << "      phi: ( phi mod 2*pi )" << std::endl;
    ss << "        z: arbitrary" << std::endl;
    ss << "---------------------------------------------------------------\n";
    ss << "Parameter Ranges:" << std::endl;
    ss << "        A: ??? " << std::endl;
    ss << "        m: ( m > 0 )" << std::endl;
    ss << "     A, m: 27m^2A^2<1" << std::endl;
    ss << "---------------------------------------------------------------\n";
    ss.precision(DEF_FIXED_REPORT_PRECISION);
    ss.setf(std::ios::fixed);

    double tau = pos[0];
    double eta = pos[1];
    double zeta = pos[3];

    ss << " e^rho = " << rho(tau, zeta, eta) << std::endl;
    ss << " e^lambda = " << lambda(tau, zeta, eta) << std::endl;
    ss << " alpha^2 = " << alpha2 << std::endl;
    ss << " q = " << q << std::endl;
    ss << "---------------------------------------------------------------\n";

    ss << "  Roots of Polynom z^3 - 1/(2*A^2)z^2 + m^2/(2*A^4):" << std::endl;
    ss << "    z1 = " << z_i[0] << std::endl;
    ss << "    z2 = " << z_i[1] << std::endl;
    ss << "    z3 = " << z_i[2] << std::endl;
    ss << "  Z1 = z1-z2 = " << Z1 << std::endl;
    ss << "  Z3 = z3-z2 = " << Z3 << std::endl;

    ss << "---------------------------------------------------------------\n";


    text = ss.str();
    return true;
}

// *************************** specific  public methods ****************************
/*! Calculates \f$e^{\rho(t, z, r)} \f$
* \param tau : t-coordinate
* \param zeta: z-coordinate
* \param eta : r-coordinate
*/
double MetricPravda_C_Can::rho(double tau, double zeta, double eta) {
    double   t1 = zeta * zeta;
    double  t2 = 2.0 * t1;
    double  t3 = tau * tau;
    double  t4 = 2.0 * t3;
    double  t5 = eta * eta;
    double  t6 = 2.0 * t5;
    double  t11 = sqrt(t2 - t4 + t6 + 4.0 * Z3 - 8.0 * Z3 * t5);
    double  t13 = t1 / 2.0;
    double  t14 = t3 / 2.0;
    double  t15 = t5 / 2.0;
    double  t17 = alpha2;
    double  t24 = sqrt(t2 - t4 + t6 + 4.0 * Z1 - 8.0 * Z1 * t5);
    return (t11 / 2.0 + t13 - t14 - t15 + Z3) / t17 / (t24 / 2.0 + t13 - t14 - t15 + Z1) / 4.0;


}

/*! Calculates \f$\frac{\partial}{\partial t}e^{\rho(t, z, r)} \f$
* \param tau : t-coordinate
* \param zeta: z-coordinate
* \param eta : r-coordinate
*/
double MetricPravda_C_Can::rho_tau(double tau, double zeta, double eta) {
    double t1 = zeta * zeta;
    double  t2 = 2.0 * t1;
    double  t3 = tau * tau;
    double  t4 = 2.0 * t3;
    double  t5 = eta * eta;
    double  t6 = 2.0 * t5;
    double  t11 = sqrt(t2 - t4 + t6 + 4.0 * Z3 - 8.0 * Z3 * t5);
    double  t15 = alpha2;
    double  t16 = 1 / t15;
    double  t22 = sqrt(t2 - t4 + t6 + 4.0 * Z1 - 8.0 * Z1 * t5);
    double  t24 = t1 / 2.0;
    double  t25 = t3 / 2.0;
    double  t26 = t5 / 2.0;
    double  t27 = t22 / 2.0 + t24 - t25 - t26 + Z1;
    double  t33 = t27 * t27;
    return (-1 / t11 * tau - tau) * t16 / t27 / 4.0 - (t11 / 2.0 + t24 - t25 - t26 + Z3) * t16 / t33 * (-1 /
            t22 * tau - tau) / 4.0;

}
/*! Calculates \f$\frac{\partial}{\partial z}e^{\rho(t, z, r)} \f$
* \param tau : t-coordinate
* \param zeta: z-coordinate
* \param eta : r-coordinate
*/
double MetricPravda_C_Can::rho_zeta(double tau, double zeta, double eta) {
    double  t1 = zeta * zeta;
    double  t2 = 2.0 * t1;
    double  t3 = tau * tau;
    double  t4 = 2.0 * t3;
    double  t5 = eta * eta;
    double  t6 = 2.0 * t5;
    double  t11 = sqrt(t2 - t4 + t6 + 4.0 * Z3 - 8.0 * Z3 * t5);
    double  t15 = alpha2;
    double  t16 = 1 / t15;
    double  t22 = sqrt(t2 - t4 + t6 + 4.0 * Z1 - 8.0 * Z1 * t5);
    double  t24 = t1 / 2.0;
    double  t25 = t3 / 2.0;
    double  t26 = t5 / 2.0;
    double  t27 = t22 / 2.0 + t24 - t25 - t26 + Z1;
    double  t33 = t27 * t27;
    return (1 / t11 * zeta + zeta) * t16 / t27 / 4.0 - (t11 / 2.0 + t24 - t25 - t26 + Z3) * t16 / t33 * (1 /
            t22 * zeta + zeta) / 4.0;
}

/*! Calculates \f$\frac{\partial}{\partial r}e^{\rho(t, z, r)} \f$
* \param tau : t-coordinate
* \param zeta: z-coordinate
* \param eta : r-coordinate
*/
double MetricPravda_C_Can::rho_eta(double tau, double zeta, double eta) {
    double  t1 = zeta * zeta;
    double  t2 = 2.0 * t1;
    double  t3 = tau * tau;
    double  t4 = 2.0 * t3;
    double  t5 = eta * eta;
    double  t6 = 2.0 * t5;
    double  t11 = sqrt(t2 - t4 + t6 + 4.0 * Z3 - 8.0 * Z3 * t5);
    double  t13 = 4.0 * eta;
    double  t20 = alpha2;
    double  t21 = 1 / t20;
    double  t27 = sqrt(t2 - t4 + t6 + 4.0 * Z1 - 8.0 * Z1 * t5);
    double  t29 = t1 / 2.0;
    double  t30 = t3 / 2.0;
    double  t31 = t5 / 2.0;
    double  t32 = t27 / 2.0 + t29 - t30 - t31 + Z1;
    double  t38 = t32 * t32;
    return (1 / t11 * (t13 - 16.0 * Z3 * eta) / 4.0 - eta) * t21 / t32 / 4.0 - (t11 / 2.0 + t29 - t30 - t31 +
            Z3) * t21 / t38 * (1 / t27 * (t13 - 16.0 * Z1 * eta) / 4.0 - eta) / 4.0;
}
/*! Calculates \f$e^{\lambda(t, z, r)} \f$
* \param tau : t-coordinate
* \param zeta: z-coordinate
* \param eta : r-coordinate
*/
double MetricPravda_C_Can::lambda(double tau, double zeta, double eta) {
    double t1 = alpha2;
    double  t2 = zeta * zeta;
    double  t3 = tau * tau;
    double  t4 = eta * eta;
    double  t5 = t2 - t3 + t4;
    double  t6 = t2 / 2.0;
    double  t7 = t3 / 2.0;
    double  t8 = t4 / 2.0;
    double  t9 = 2.0 * t2;
    double  t10 = 2.0 * t3;
    double  t11 = 2.0 * t4;
    double  t13 = Z1 * t4;
    double  t16 = sqrt(t9 - t10 + t11 + 4.0 * Z1 - 8.0 * t13);
    double  t23 = Z3 * t4;
    double  t26 = sqrt(t9 - t10 + t11 + 4.0 * Z3 - 8.0 * t23);
    return 8.0 * t1 * (t5 * (t6 - t7 + t8 + t16 / 2.0 + Z1) / 2.0 - t13) * (t16 * t26 / 4.0 + (t6 - t7 + t8 + Z3
                                                                                              ) * (t6 - t7 + t8 + Z1) - (Z1 + Z3) * t4) / t16 / t26 / (t5 * (t6 - t7 + t8 + t26 / 2.0 + Z3) / 2.0 - t23);

}

/*! Calculates \f$\frac{\partial}{\partial t}e^{\lambda(t, z, r)} \f$
* \param tau : t-coordinate
* \param zeta: z-coordinate
* \param eta : r-coordinate
*/
double MetricPravda_C_Can::lambda_tau(double tau, double zeta, double eta) { //diff(lambda(tau,zeta,eta),tau)
    double  t1 = alpha2;
    double  t2 = zeta * zeta;
    double  t3 = t2 / 2.0;
    double  t4 = tau * tau;
    double  t5 = t4 / 2.0;
    double  t6 = eta * eta;
    double  t7 = t6 / 2.0;
    double  t8 = 2.0 * t2;
    double  t9 = 2.0 * t4;
    double  t10 = 2.0 * t6;
    double  t12 = Z1 * t6;
    double  t14 = t8 - t9 + t10 + 4.0 * Z1 - 8.0 * t12;
    double  t15 = sqrt(t14);
    double  t17 = t3 - t5 + t7 + t15 / 2.0 + Z1;
    double  t19 = t2 - t4 + t6;
    double  t20 = 1 / t15;
    double  t27 = Z3 * t6;
    double  t29 = t8 - t9 + t10 + 4.0 * Z3 - 8.0 * t27;
    double  t30 = sqrt(t29);
    double  t33 = t3 - t5 + t7 + Z1;
    double  t34 = t3 - t5 + t7 + Z3;
    double  t38 = t15 * t30 / 4.0 + t33 * t34 - (Z1 + Z3) * t6;
    double  t40 = 1 / t30;
    double  t41 = t20 * t40;
    double  t43 = t3 - t5 + t7 + t30 / 2.0 + Z3;
    double  t45 = t19 * t43 / 2.0 - t27;
    double  t46 = 1 / t45;
    double  t47 = t41 * t46;
    double  t52 = t1 * (t19 * t17 / 2.0 - t12);
    double  t65 = t52 * t38;
    double  t69 = t46 * tau;
    double  t79 = t45 * t45;
    return 8.0 * t1 * (-tau * t17 + t19 * (-t20 * tau - tau) / 2.0) * t38 * t47 + 8.0 * t52 * (-t20 * t30 *
            tau / 2.0 - t15 * t40 * tau / 2.0 - tau * t34 - t33 * tau) * t47 + 16.0 * t65 / t15 / t14 * t40 * t69 + 16.0 * t65 *
           t20 / t30 / t29 * t69 - 8.0 * t65 * t41 / t79 * (-tau * t43 + t19 * (-t40 * tau - tau) / 2.0);
}

/*! Calculates \f$\frac{\partial}{\partial z}e^{\lambda(t, z, r)} \f$
*
* \param tau : t-coordinate
* \param zeta: z-coordinate
* \param eta : r-coordinate
*/
double MetricPravda_C_Can::lambda_zeta(double tau, double zeta, double eta) {
    double  t1 = alpha2;
    double  t2 = zeta * zeta;
    double  t3 = t2 / 2.0;
    double  t4 = tau * tau;
    double  t5 = t4 / 2.0;
    double  t6 = eta * eta;
    double  t7 = t6 / 2.0;
    double  t8 = 2.0 * t2;
    double  t9 = 2.0 * t4;
    double  t10 = 2.0 * t6;
    double  t12 = Z1 * t6;
    double  t14 = t8 - t9 + t10 + 4.0 * Z1 - 8.0 * t12;
    double  t15 = sqrt(t14);
    double  t17 = t3 - t5 + t7 + t15 / 2.0 + Z1;
    double  t19 = t2 - t4 + t6;
    double  t20 = 1 / t15;
    double  t27 = Z3 * t6;
    double  t29 = t8 - t9 + t10 + 4.0 * Z3 - 8.0 * t27;
    double  t30 = sqrt(t29);
    double  t33 = t3 - t5 + t7 + Z1;
    double  t34 = t3 - t5 + t7 + Z3;
    double  t38 = t15 * t30 / 4.0 + t33 * t34 - (Z1 + Z3) * t6;
    double  t40 = 1 / t30;
    double  t41 = t20 * t40;
    double  t43 = t3 - t5 + t7 + t30 / 2.0 + Z3;
    double  t45 = t19 * t43 / 2.0 - t27;
    double  t46 = 1 / t45;
    double  t47 = t41 * t46;
    double  t52 = t1 * (t19 * t17 / 2.0 - t12);
    double  t65 = t52 * t38;
    double  t69 = t46 * zeta;
    double  t79 = t45 * t45;
    return  8.0 * t1 * (zeta * t17 + t19 * (t20 * zeta + zeta) / 2.0) * t38 * t47 + 8.0 * t52 * (t20 * t30 *
            zeta / 2.0 + t15 * t40 * zeta / 2.0 + zeta * t34 + t33 * zeta) * t47 - 16.0 * t65 / t15 / t14 * t40 * t69 - 16.0 *
            t65 * t20 / t30 / t29 * t69 - 8.0 * t65 * t41 / t79 * (zeta * t43 + t19 * (t40 * zeta + zeta) / 2.0);
}

/*! Calculates \f$\frac{\partial}{\partial r}e^{\lambda(t, z, r)} \f$
* \param tau : t-coordinate
* \param zeta: z-coordinate
* \param eta : r-coordinate
*/
double MetricPravda_C_Can::lambda_eta(double tau, double zeta, double eta) {
    double t1 = alpha2;
    double  t2 = zeta * zeta;
    double  t3 = t2 / 2.0;
    double  t4 = tau * tau;
    double  t5 = t4 / 2.0;
    double  t6 = eta * eta;
    double  t7 = t6 / 2.0;
    double  t8 = 2.0 * t2;
    double  t9 = 2.0 * t4;
    double  t10 = 2.0 * t6;
    double  t12 = Z1 * t6;
    double  t14 = t8 - t9 + t10 + 4.0 * Z1 - 8.0 * t12;
    double  t15 = sqrt(t14);
    double  t17 = t3 - t5 + t7 + t15 / 2.0 + Z1;
    double  t19 = t2 - t4 + t6;
    double  t20 = 1 / t15;
    double  t21 = 4.0 * eta;
    double  t22 = Z1 * eta;
    double  t24 = t21 - 16.0 * t22;
    double  t33 = Z3 * t6;
    double  t35 = t8 - t9 + t10 + 4.0 * Z3 - 8.0 * t33;
    double  t36 = sqrt(t35);
    double  t39 = t3 - t5 + t7 + Z1;
    double  t40 = t3 - t5 + t7 + Z3;
    double  t42 = Z1 + Z3;
    double  t44 = t15 * t36 / 4.0 + t39 * t40 - t42 * t6;
    double  t46 = 1 / t36;
    double  t47 = t20 * t46;
    double  t49 = t3 - t5 + t7 + t36 / 2.0 + Z3;
    double  t51 = t19 * t49 / 2.0 - t33;
    double  t52 = 1 / t51;
    double  t53 = t47 * t52;
    double  t58 = t1 * (t19 * t17 / 2.0 - t12);
    double  t63 = Z3 * eta;
    double  t65 = t21 - 16.0 * t63;
    double  t76 = t58 * t44;
    double  t91 = t51 * t51;
    return 8.0 * t1 * (eta * t17 + t19 * (eta + t20 * t24 / 4.0) / 2.0 - 2.0 * t22) * t44 * t53 + 8.0 * t58
           * (t20 * t36 * t24 / 8.0 + t15 * t46 * t65 / 8.0 + eta * t40 + t39 * eta - 2.0 * t42 * eta) * t53 - 4.0 * t76 / t15 /
           t14 * t46 * t52 * t24 - 4.0 * t76 * t20 / t36 / t35 * t52 * t65 - 8.0 * t76 * t47 / t91 * (eta * t49 + t19 * (eta +
                   t46 * t65 / 4.0) / 2.0 - 2.0 * t63);
}

// ********************************* protected methods *****************************
/*!
*/
void
MetricPravda_C_Can::setStandardValues() {
    mInitPos[0] = 1.0;
    mInitPos[1] = 1.0;
    mInitPos[2] = 0.0;
    mInitPos[3] = 0.0;
    mInitDir[0] = 1.0;
    mInitDir[1] = 0.0;
    mInitDir[2] = 0.0;

    mCoordNames[0] = std::string("t");
    mCoordNames[1] = std::string("r");
    mCoordNames[2] = std::string("phi");
    mCoordNames[3] = std::string("z");
}


} // end namespace m4d

