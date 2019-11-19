// -------------------------------------------------------------------------------
/*
   m4dMetricKastorTraschen.cpp

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
// -------------------------------------------------------------------------------

#include "m4dMetricKastorTraschen.h"

namespace m4d {

MetricKastorTraschen::MetricKastorTraschen(double H)
{
    mMetricName = "KastorTraschen";
    setCoordType(enum_coordinate_cartesian);

    mPhysicalUnits = enum_physical_constants_geom;
    mSpeedOfLight = 1.0;
    mGravConstant = 1.0;

    /*  Only a static tetrad is defined  */
    mLocTeds.push_back(enum_nat_tetrad_static);

    m1 = 1.0;
    z1 = 1.0;
    m2 = 1.0;
    z2 = -1.0;
    mH = H;

    addParam("m1", m1);
    addParam("z1", z1);
    addParam("m2", m2);
    addParam("z2", z2);
    addParam("h", mH);

    setStandardValues();
}

MetricKastorTraschen::~MetricKastorTraschen() {}

bool MetricKastorTraschen::calculateMetric(const double* pos)
{
    double Omega, a;
    calcPotentials(pos, Omega, a);

    double t1 = Omega; // Omega(t,x,y,z);
    double t2 = t1 * t1;
    double t4 = a; // a(t);
    double t5 = t4 * t4;
    double t6 = t2 * t5;

    g_compts[0][0] = -1 / t2;
    g_compts[0][1] = 0.0;
    g_compts[0][2] = 0.0;
    g_compts[0][3] = 0.0;
    g_compts[1][0] = 0.0;
    g_compts[1][1] = t6;
    g_compts[1][2] = 0.0;
    g_compts[1][3] = 0.0;
    g_compts[2][0] = 0.0;
    g_compts[2][1] = 0.0;
    g_compts[2][2] = t6;
    g_compts[2][3] = 0.0;
    g_compts[3][0] = 0.0;
    g_compts[3][1] = 0.0;
    g_compts[3][2] = 0.0;
    g_compts[3][3] = t6;

    return true;
}

bool MetricKastorTraschen::calculateChristoffels(const double* pos)
{
    double Omega, a;
    double dOdx, dOdy, dOdz, dOdt, dadt;
    calcPotDiffs(pos, Omega, a, dOdx, dOdy, dOdz, dOdt, dadt);

    double t1 = Omega; // Omega(t,x,y,z);
    double t2 = 1 / t1;
    double t3 = dOdt; // diff(Omega(t,x,y,z),t);
    double t5 = t1 * t1;
    double t6 = t5 * t5;
    double t9 = a; // a(t);
    double t10 = t9 * t9;
    double t12 = 1 / t6 / t1 / t10;
    double t13 = dOdx; // diff(Omega(t,x,y,z),x);
    double t15 = dOdy; // diff(Omega(t,x,y,z),y);
    double t17 = dOdz; // diff(Omega(t,x,y,z),z);
    double t19 = t2 * t13;
    double t23 = dadt; // diff(a(t),t);
    double t25 = t9 * t3 + t1 * t23;
    double t26 = t2 / t9 * t25;
    double t27 = t2 * t15;
    double t28 = t2 * t17;
    double t31 = t5 * t1 * t9 * t25;

    christoffel[0][0][0] = -t2 * t3;
    christoffel[0][0][1] = -t12 * t13;
    christoffel[0][0][2] = -t12 * t15;
    christoffel[0][0][3] = -t12 * t17;
    christoffel[0][1][0] = -t19;
    christoffel[0][1][1] = t26;
    christoffel[0][1][2] = 0.0;
    christoffel[0][1][3] = 0.0;
    christoffel[0][2][0] = -t27;
    christoffel[0][2][1] = 0.0;
    christoffel[0][2][2] = t26;
    christoffel[0][2][3] = 0.0;
    christoffel[0][3][0] = -t28;
    christoffel[0][3][1] = 0.0;
    christoffel[0][3][2] = 0.0;
    christoffel[0][3][3] = t26;
    christoffel[1][0][0] = -t19;
    christoffel[1][0][1] = t26;
    christoffel[1][0][2] = 0.0;
    christoffel[1][0][3] = 0.0;
    christoffel[1][1][0] = t31;
    christoffel[1][1][1] = t19;
    christoffel[1][1][2] = -t27;
    christoffel[1][1][3] = -t28;
    christoffel[1][2][0] = 0.0;
    christoffel[1][2][1] = t27;
    christoffel[1][2][2] = t19;
    christoffel[1][2][3] = 0.0;
    christoffel[1][3][0] = 0.0;
    christoffel[1][3][1] = t28;
    christoffel[1][3][2] = 0.0;
    christoffel[1][3][3] = t19;
    christoffel[2][0][0] = -t27;
    christoffel[2][0][1] = 0.0;
    christoffel[2][0][2] = t26;
    christoffel[2][0][3] = 0.0;
    christoffel[2][1][0] = 0.0;
    christoffel[2][1][1] = t27;
    christoffel[2][1][2] = t19;
    christoffel[2][1][3] = 0.0;
    christoffel[2][2][0] = t31;
    christoffel[2][2][1] = -t19;
    christoffel[2][2][2] = t27;
    christoffel[2][2][3] = -t28;
    christoffel[2][3][0] = 0.0;
    christoffel[2][3][1] = 0.0;
    christoffel[2][3][2] = t28;
    christoffel[2][3][3] = t27;
    christoffel[3][0][0] = -t28;
    christoffel[3][0][1] = 0.0;
    christoffel[3][0][2] = 0.0;
    christoffel[3][0][3] = t26;
    christoffel[3][1][0] = 0.0;
    christoffel[3][1][1] = t28;
    christoffel[3][1][2] = 0.0;
    christoffel[3][1][3] = t19;
    christoffel[3][2][0] = 0.0;
    christoffel[3][2][1] = 0.0;
    christoffel[3][2][2] = t28;
    christoffel[3][2][3] = t27;
    christoffel[3][3][0] = t31;
    christoffel[3][3][1] = -t19;
    christoffel[3][3][2] = -t27;
    christoffel[3][3][3] = t28;

    return true;
}

bool MetricKastorTraschen::calculateChrisD(const double*)
{
    return false;
    /*
      double t1 = diff(Omega(t,x,y,z),t);
      double t2 = t1*t1;
      double t3 = diff(diff(Omega(t,x,y,z),t),t);
      double t4 = Omega(t,x,y,z);
      double t7 = t4*t4;
      double t8 = 1/t7;
      double t10 = diff(Omega(t,x,y,z),x);
      double t11 = t1*t10;
      double t12 = diff(diff(Omega(t,x,y,z),t),x);
      double t13 = t12*t4;
      double t15 = (-t11+t13)*t8;
      double t16 = diff(Omega(t,x,y,z),y);
      double t17 = t1*t16;
      double t18 = diff(diff(Omega(t,x,y,z),t),y);
      double t19 = t18*t4;
      double t21 = (-t17+t19)*t8;
      double t22 = diff(Omega(t,x,y,z),z);
      double t23 = t1*t22;
      double t24 = diff(diff(Omega(t,x,y,z),t),z);
      double t25 = t24*t4;
      double t27 = (-t23+t25)*t8;
      double t28 = a(t);
      double t31 = diff(a(t),t);
      double t32 = t10*t31;
      double t37 = t7*t7;
      double t39 = 1/t37/t7;
      double t41 = t28*t28;
      double t43 = 1/t41/t28;
      double t45 = t10*t10;
      double t47 = diff(diff(Omega(t,x,y,z),x),x);
      double t48 = t47*t4;
      double t51 = 1/t41;
      double t53 = t10*t16;
      double t55 = diff(diff(Omega(t,x,y,z),x),y);
      double t56 = t55*t4;
      double t59 = (-5.0*t53+t56)*t39*t51;
      double t60 = t10*t22;
      double t62 = diff(diff(Omega(t,x,y,z),x),z);
      double t63 = t62*t4;
      double t66 = (-5.0*t60+t63)*t39*t51;
      double t69 = t16*t31;
      double t76 = t16*t16;
      double t78 = diff(diff(Omega(t,x,y,z),y),y);
      double t79 = t78*t4;
      double t83 = t16*t22;
      double t85 = diff(diff(Omega(t,x,y,z),y),z);
      double t86 = t85*t4;
      double t89 = (-5.0*t83+t86)*t39*t51;
      double t92 = t22*t31;
      double t99 = t22*t22;
      double t101 = diff(diff(Omega(t,x,y,z),z),z);
      double t102 = t101*t4;
      double t107 = (-t45+t48)*t8;
      double t109 = (-t53+t56)*t8;
      double t111 = (-t60+t63)*t8;
      double t113 = t31*t31;
      double t118 = diff(diff(a(t),t),t);
      double t122 = (-t41*t2-t7*t113+t4*t41*t3+t7*t28*t118)*t8*t51;
      double t124 = (-t76+t79)*t8;
      double t126 = (-t83+t86)*t8;
      double t128 = (-t99+t102)*t8;
      double t129 = t7*t41;
      double t132 = t7*t4;
      double t133 = t132*t28;
      double t138 = t132*t41;
      double t142 = 3.0*t129*t2+6.0*t133*t1*t31+t37*t113+t138*t3+t37*t28*t118;
      double t148 = 3.0*t129*t11+4.0*t133*t32+t138*t12;
      double t154 = 3.0*t129*t17+4.0*t133*t69+t138*t18;
      double t160 = 3.0*t129*t23+4.0*t133*t92+t138*t24;
      chrisD[0][0][0][0] = -(-t2+t3*t4)*t8;
      chrisD[0][0][0][1] = -t15;
      chrisD[0][0][0][2] = -t21;
      chrisD[0][0][0][3] = -t27;
      chrisD[0][0][1][0] = (5.0*t11*t28+2.0*t32*t4-t13*t28)*t39*t43;
      chrisD[0][0][1][1] = -(-5.0*t45+t48)*t39*t51;
      chrisD[0][0][1][2] = -t59;
      chrisD[0][0][1][3] = -t66;
      chrisD[0][0][2][0] = (5.0*t17*t28+2.0*t69*t4-t19*t28)*t39*t43;
      chrisD[0][0][2][1] = -t59;
      chrisD[0][0][2][2] = -(-5.0*t76+t79)*t39*t51;
      chrisD[0][0][2][3] = -t89;
      chrisD[0][0][3][0] = (5.0*t23*t28+2.0*t92*t4-t25*t28)*t39*t43;
      chrisD[0][0][3][1] = -t66;
      chrisD[0][0][3][2] = -t89;
      chrisD[0][0][3][3] = -(-5.0*t99+t102)*t39*t51;
      chrisD[0][1][0][0] = -t15;
      chrisD[0][1][0][1] = -t107;
      chrisD[0][1][0][2] = -t109;
      chrisD[0][1][0][3] = -t111;
      chrisD[0][1][1][0] = t122;
      chrisD[0][1][1][1] = t15;
      chrisD[0][1][1][2] = t21;
      chrisD[0][1][1][3] = t27;
      chrisD[0][1][2][0] = 0.0;
      chrisD[0][1][2][1] = 0.0;
      chrisD[0][1][2][2] = 0.0;
      chrisD[0][1][2][3] = 0.0;
      chrisD[0][1][3][0] = 0.0;
      chrisD[0][1][3][1] = 0.0;
      chrisD[0][1][3][2] = 0.0;
      chrisD[0][1][3][3] = 0.0;
      chrisD[0][2][0][0] = -t21;
      chrisD[0][2][0][1] = -t109;
      chrisD[0][2][0][2] = -t124;
      chrisD[0][2][0][3] = -t126;
      chrisD[0][2][1][0] = 0.0;
      chrisD[0][2][1][1] = 0.0;
      chrisD[0][2][1][2] = 0.0;
      chrisD[0][2][1][3] = 0.0;
      chrisD[0][2][2][0] = t122;
      chrisD[0][2][2][1] = t15;
      chrisD[0][2][2][2] = t21;
      chrisD[0][2][2][3] = t27;
      chrisD[0][2][3][0] = 0.0;
      chrisD[0][2][3][1] = 0.0;
      chrisD[0][2][3][2] = 0.0;
      chrisD[0][2][3][3] = 0.0;
      chrisD[0][3][0][0] = -t27;
      chrisD[0][3][0][1] = -t111;
      chrisD[0][3][0][2] = -t126;
      chrisD[0][3][0][3] = -t128;
      chrisD[0][3][1][0] = 0.0;
      chrisD[0][3][1][1] = 0.0;
      chrisD[0][3][1][2] = 0.0;
      chrisD[0][3][1][3] = 0.0;
      chrisD[0][3][2][0] = 0.0;
      chrisD[0][3][2][1] = 0.0;
      chrisD[0][3][2][2] = 0.0;
      chrisD[0][3][2][3] = 0.0;
      chrisD[0][3][3][0] = t122;
      chrisD[0][3][3][1] = t15;
      chrisD[0][3][3][2] = t21;
      chrisD[0][3][3][3] = t27;
      chrisD[1][0][0][0] = -t15;
      chrisD[1][0][0][1] = -t107;
      chrisD[1][0][0][2] = -t109;
      chrisD[1][0][0][3] = -t111;
      chrisD[1][0][1][0] = t122;
      chrisD[1][0][1][1] = t15;
      chrisD[1][0][1][2] = t21;
      chrisD[1][0][1][3] = t27;
      chrisD[1][0][2][0] = 0.0;
      chrisD[1][0][2][1] = 0.0;
      chrisD[1][0][2][2] = 0.0;
      chrisD[1][0][2][3] = 0.0;
      chrisD[1][0][3][0] = 0.0;
      chrisD[1][0][3][1] = 0.0;
      chrisD[1][0][3][2] = 0.0;
      chrisD[1][0][3][3] = 0.0;
      chrisD[1][1][0][0] = t142;
      chrisD[1][1][0][1] = t148;
      chrisD[1][1][0][2] = t154;
      chrisD[1][1][0][3] = t160;
      chrisD[1][1][1][0] = t15;
      chrisD[1][1][1][1] = t107;
      chrisD[1][1][1][2] = t109;
      chrisD[1][1][1][3] = t111;
      chrisD[1][1][2][0] = -t21;
      chrisD[1][1][2][1] = -t109;
      chrisD[1][1][2][2] = -t124;
      chrisD[1][1][2][3] = -t126;
      chrisD[1][1][3][0] = -t27;
      chrisD[1][1][3][1] = -t111;
      chrisD[1][1][3][2] = -t126;
      chrisD[1][1][3][3] = -t128;
      chrisD[1][2][0][0] = 0.0;
      chrisD[1][2][0][1] = 0.0;
      chrisD[1][2][0][2] = 0.0;
      chrisD[1][2][0][3] = 0.0;
      chrisD[1][2][1][0] = t21;
      chrisD[1][2][1][1] = t109;
      chrisD[1][2][1][2] = t124;
      chrisD[1][2][1][3] = t126;
      chrisD[1][2][2][0] = t15;
      chrisD[1][2][2][1] = t107;
      chrisD[1][2][2][2] = t109;
      chrisD[1][2][2][3] = t111;
      chrisD[1][2][3][0] = 0.0;
      chrisD[1][2][3][1] = 0.0;
      chrisD[1][2][3][2] = 0.0;
      chrisD[1][2][3][3] = 0.0;
      chrisD[1][3][0][0] = 0.0;
      chrisD[1][3][0][1] = 0.0;
      chrisD[1][3][0][2] = 0.0;
      chrisD[1][3][0][3] = 0.0;
      chrisD[1][3][1][0] = t27;
      chrisD[1][3][1][1] = t111;
      chrisD[1][3][1][2] = t126;
      chrisD[1][3][1][3] = t128;
      chrisD[1][3][2][0] = 0.0;
      chrisD[1][3][2][1] = 0.0;
      chrisD[1][3][2][2] = 0.0;
      chrisD[1][3][2][3] = 0.0;
      chrisD[1][3][3][0] = t15;
      chrisD[1][3][3][1] = t107;
      chrisD[1][3][3][2] = t109;
      chrisD[1][3][3][3] = t111;
      chrisD[2][0][0][0] = -t21;
      chrisD[2][0][0][1] = -t109;
      chrisD[2][0][0][2] = -t124;
      chrisD[2][0][0][3] = -t126;
      chrisD[2][0][1][0] = 0.0;
      chrisD[2][0][1][1] = 0.0;
      chrisD[2][0][1][2] = 0.0;
      chrisD[2][0][1][3] = 0.0;
      chrisD[2][0][2][0] = t122;
      chrisD[2][0][2][1] = t15;
      chrisD[2][0][2][2] = t21;
      chrisD[2][0][2][3] = t27;
      chrisD[2][0][3][0] = 0.0;
      chrisD[2][0][3][1] = 0.0;
      chrisD[2][0][3][2] = 0.0;
      chrisD[2][0][3][3] = 0.0;
      chrisD[2][1][0][0] = 0.0;
      chrisD[2][1][0][1] = 0.0;
      chrisD[2][1][0][2] = 0.0;
      chrisD[2][1][0][3] = 0.0;
      chrisD[2][1][1][0] = t21;
      chrisD[2][1][1][1] = t109;
      chrisD[2][1][1][2] = t124;
      chrisD[2][1][1][3] = t126;
      chrisD[2][1][2][0] = t15;
      chrisD[2][1][2][1] = t107;
      chrisD[2][1][2][2] = t109;
      chrisD[2][1][2][3] = t111;
      chrisD[2][1][3][0] = 0.0;
      chrisD[2][1][3][1] = 0.0;
      chrisD[2][1][3][2] = 0.0;
      chrisD[2][1][3][3] = 0.0;
      chrisD[2][2][0][0] = t142;
      chrisD[2][2][0][1] = t148;
      chrisD[2][2][0][2] = t154;
      chrisD[2][2][0][3] = t160;
      chrisD[2][2][1][0] = -t15;
      chrisD[2][2][1][1] = -t107;
      chrisD[2][2][1][2] = -t109;
      chrisD[2][2][1][3] = -t111;
      chrisD[2][2][2][0] = t21;
      chrisD[2][2][2][1] = t109;
      chrisD[2][2][2][2] = t124;
      chrisD[2][2][2][3] = t126;
      chrisD[2][2][3][0] = -t27;
      chrisD[2][2][3][1] = -t111;
      chrisD[2][2][3][2] = -t126;
      chrisD[2][2][3][3] = -t128;
      chrisD[2][3][0][0] = 0.0;
      chrisD[2][3][0][1] = 0.0;
      chrisD[2][3][0][2] = 0.0;
      chrisD[2][3][0][3] = 0.0;
      chrisD[2][3][1][0] = 0.0;
      chrisD[2][3][1][1] = 0.0;
      chrisD[2][3][1][2] = 0.0;
      chrisD[2][3][1][3] = 0.0;
      chrisD[2][3][2][0] = t27;
      chrisD[2][3][2][1] = t111;
      chrisD[2][3][2][2] = t126;
      chrisD[2][3][2][3] = t128;
      chrisD[2][3][3][0] = t21;
      chrisD[2][3][3][1] = t109;
      chrisD[2][3][3][2] = t124;
      chrisD[2][3][3][3] = t126;
      chrisD[3][0][0][0] = -t27;
      chrisD[3][0][0][1] = -t111;
      chrisD[3][0][0][2] = -t126;
      chrisD[3][0][0][3] = -t128;
      chrisD[3][0][1][0] = 0.0;
      chrisD[3][0][1][1] = 0.0;
      chrisD[3][0][1][2] = 0.0;
      chrisD[3][0][1][3] = 0.0;
      chrisD[3][0][2][0] = 0.0;
      chrisD[3][0][2][1] = 0.0;
      chrisD[3][0][2][2] = 0.0;
      chrisD[3][0][2][3] = 0.0;
      chrisD[3][0][3][0] = t122;
      chrisD[3][0][3][1] = t15;
      chrisD[3][0][3][2] = t21;
      chrisD[3][0][3][3] = t27;
      chrisD[3][1][0][0] = 0.0;
      chrisD[3][1][0][1] = 0.0;
      chrisD[3][1][0][2] = 0.0;
      chrisD[3][1][0][3] = 0.0;
      chrisD[3][1][1][0] = t27;
      chrisD[3][1][1][1] = t111;
      chrisD[3][1][1][2] = t126;
      chrisD[3][1][1][3] = t128;
      chrisD[3][1][2][0] = 0.0;
      chrisD[3][1][2][1] = 0.0;
      chrisD[3][1][2][2] = 0.0;
      chrisD[3][1][2][3] = 0.0;
      chrisD[3][1][3][0] = t15;
      chrisD[3][1][3][1] = t107;
      chrisD[3][1][3][2] = t109;
      chrisD[3][1][3][3] = t111;
      chrisD[3][2][0][0] = 0.0;
      chrisD[3][2][0][1] = 0.0;
      chrisD[3][2][0][2] = 0.0;
      chrisD[3][2][0][3] = 0.0;
      chrisD[3][2][1][0] = 0.0;
      chrisD[3][2][1][1] = 0.0;
      chrisD[3][2][1][2] = 0.0;
      chrisD[3][2][1][3] = 0.0;
      chrisD[3][2][2][0] = t27;
      chrisD[3][2][2][1] = t111;
      chrisD[3][2][2][2] = t126;
      chrisD[3][2][2][3] = t128;
      chrisD[3][2][3][0] = t21;
      chrisD[3][2][3][1] = t109;
      chrisD[3][2][3][2] = t124;
      chrisD[3][2][3][3] = t126;
      chrisD[3][3][0][0] = t142;
      chrisD[3][3][0][1] = t148;
      chrisD[3][3][0][2] = t154;
      chrisD[3][3][0][3] = t160;
      chrisD[3][3][1][0] = -t15;
      chrisD[3][3][1][1] = -t107;
      chrisD[3][3][1][2] = -t109;
      chrisD[3][3][1][3] = -t111;
      chrisD[3][3][2][0] = -t21;
      chrisD[3][3][2][1] = -t109;
      chrisD[3][3][2][2] = -t124;
      chrisD[3][3][2][3] = -t126;
      chrisD[3][3][3][0] = t27;
      chrisD[3][3][3][1] = t111;
      chrisD[3][3][3][2] = t126;
      chrisD[3][3][3][3] = t128;
    */
    return true;
}

/*! Transform local 4-direction to coordinate 4-direction.
 *
 *  \param  pos  :  pointer to position array.
 *  \param  ldir :  pointer to local direction array.
 *  \param  dir  :  pointer to calculated coordinate direction array.
 *  \param  type :  type of tetrad.
 */
void MetricKastorTraschen::localToCoord(const double* pos, const double* ldir, double* dir, enum_nat_tetrad_type)
{
    double Omega, a;
    calcPotentials(pos, Omega, a);

    double edOa = 1.0 / (Omega * a);
    dir[0] = Omega * ldir[0];
    dir[1] = edOa * ldir[1];
    dir[2] = edOa * ldir[2];
    dir[3] = edOa * ldir[3];
}

/*! Transform coordinate 4-direction to local 4-direction.
 *
 *  \param  pos  :  pointer to position array.
 *  \param  cdir :  pointer to coordinate direction.
 *  \param  ldir :  pointer to calculated local direction array.
 *  \param  type :  type of tetrad.
 */
void MetricKastorTraschen::coordToLocal(const double* pos, const double* cdir, double* ldir, enum_nat_tetrad_type)
{
    double Omega, a;
    calcPotentials(pos, Omega, a);

    double Oa = Omega * a;

    ldir[0] = cdir[0] / Omega;
    ldir[1] = cdir[1] * Oa;
    ldir[2] = cdir[2] * Oa;
    ldir[3] = cdir[3] * Oa;
}

/*! Test break condition.
 *
 *  \param pos    : pointer to position array.
 *  \return false : position is valid.
 */
bool MetricKastorTraschen::breakCondition(const double*)
{
    bool br = false;
    return br;
}

double MetricKastorTraschen::testConstraint(const double y[], const double kappa)
{
    double Omega, a;
    calcPotentials(y, Omega, a);

    double dt = y[4];
    double dx = y[5];
    double dy = y[6];
    double dz = y[7];

    double sum = -kappa;
    sum += -dt * dt / (Omega * Omega) + Omega * Omega * a * a * (dx * dx + dy * dy + dz * dz);
    // std::cerr << Omega << " " << a << " " << sum << std::endl;
    return sum;
}

/*! Set parameter 'pName' to 'val'.
 *
 *  Set 'mass' or 'lambda' parameter.
 */
bool MetricKastorTraschen::setParam(const char* pName, double val)
{
    Metric::setParam(pName, val);
    if (strcmp(pName, "m1") == 0) {
        m1 = val;
    }
    else if (strcmp(pName, "z1") == 0) {
        z1 = val;
    }
    else if (strcmp(pName, "m2") == 0) {
        m2 = val;
    }
    else if (strcmp(pName, "z2") == 0) {
        z2 = val;
    }
    else if (strcmp(pName, "h") == 0) {
        mH = val;
    }
    return true;
}

/*! Generate report.
 */
bool MetricKastorTraschen::report(const vec4, const vec4, char*& text)
{
    std::stringstream ss;
    ss << "Report for Kastor-Traschen metric\n\tcoordinate : (t,x,y,z)\n";
    ss << "---------------------------------------------------------------\n";
    ss << "  physical units ..................... no\n";
    ss.precision(DEF_FIXED_REPORT_PRECISION);
    ss.setf(std::ios::fixed);

    text = new char[ss.str().length() + 2];
    return CopyString(ss.str().c_str(), text);
}

void MetricKastorTraschen::calcPotentials(const double* pos, double& Omega, double& a)
{
    double t = pos[0];
    double x = pos[1];
    double y = pos[2];
    double z = pos[3];

    double r1 = sqrt(x * x + y * y + (z - z1) * (z - z1));
    double r2 = sqrt(x * x + y * y + (z - z2) * (z - z2));

    a = exp(mH * t);
    Omega = 1.0 + m1 / (r1 * a) + m2 / (r2 * a);
}

void MetricKastorTraschen::calcPotDiffs(
    const double* pos, double& Omega, double& a, double& dOdx, double& dOdy, double& dOdz, double& dOdt, double& dadt)
{
    double t = pos[0];
    double x = pos[1];
    double y = pos[2];
    double z = pos[3];

    double r1 = sqrt(x * x + y * y + (z - z1) * (z - z1));
    double r2 = sqrt(x * x + y * y + (z - z2) * (z - z2));

    a = exp(mH * t);
    Omega = 1.0 + m1 / (r1 * a) + m2 / (r2 * a);

    double dr1dx = x / r1;
    double dr1dy = y / r1;
    double dr1dz = (z - z1) / r1;

    double dr2dx = x / r2;
    double dr2dy = y / r2;
    double dr2dz = (z - z2) / r2;

    dOdx = -m1 / (r1 * r1 * a) * dr1dx - m2 / (r2 * r2 * a) * dr2dx;
    dOdy = -m1 / (r1 * r1 * a) * dr1dy - m2 / (r2 * r2 * a) * dr2dy;
    dOdz = -m1 / (r1 * r1 * a) * dr1dz - m2 / (r2 * r2 * a) * dr2dz;

    dadt = mH * exp(mH * t);
    dOdt = -m1 * dadt / (r1 * a * a) - m2 * dadt / (r2 * a * a);
}

// ********************************* protected methods *****************************
/*!
 */
void MetricKastorTraschen::setStandardValues()
{
    mInitPos[0] = 0.0;
    mInitPos[1] = 30.0;
    mInitPos[2] = 0.0;
    mInitPos[3] = 0.0;
    mInitDir[0] = 1.0;
    mInitDir[1] = 0.0;
    mInitDir[2] = 0.0;

    mCoordNames[0] = std::string("t");
    mCoordNames[1] = std::string("x");
    mCoordNames[2] = std::string("y");
    mCoordNames[3] = std::string("z");
}

} // end namespace m4d
