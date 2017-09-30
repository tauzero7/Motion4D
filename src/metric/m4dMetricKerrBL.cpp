// -------------------------------------------------------------------------------
/*
   m4dMetricKerrBL.cpp

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

#include <cmath>
#include "m4dMetricKerrBL.h"

namespace m4d {

#define eps 1.0e-6


/*! Standard constructor for the KerrBL metric.
 *
 * \param  m : mass of the black hole.
 * \param  a : angular momentum of the black hole.
 */
MetricKerrBL::MetricKerrBL(double m, double a) {
    mMetricName  = "KerrBL";
    setCoordType(enum_coordinate_spherical);

    mPhysicalUnits = enum_physical_constants_geom;
    mSpeedOfLight = 1.0;
    mGravConstant = 1.0;

    mass = m;
    angmom = a;

    addParam("mass", mass);
    addParam("angmom", angmom);

    setStandardValues();

    mLocTeds.push_back(enum_nat_tetrad_static);
    mLocTeds.push_back(enum_nat_tetrad_lnrf);
    pm = 1.0;

    mDrawTypes.push_back(enum_draw_twoplusone);
    mDrawTypes.push_back(enum_draw_effpoti);
}

MetricKerrBL::~MetricKerrBL() {
}


// *********************************** public methods ******************************
/*! Calculate the contravariant metric components at position 'pos'.
 *
 *  \param pos : pointer to position.
 */
bool MetricKerrBL::calculateMetric(const double* pos) {
    double r     = pos[1];
    double theta = pos[2];
    double m = mass;
    double a = angmom;

    double t1 = m * r;
    double t2 = r * r;
    double t3 = a * a;
    double t4 = cos(theta);
    double t5 = t4 * t4;
    double t7 = t2 + t3 * t5;
    double t8 = 1 / t7;
    double t13 = sin(theta);
    double t14 = t13 * t13;
    double t18 = 2.0 * m * a * r * t14 * t8;
    double t21 = 1 / (t2 - 2.0 * t1 + t3);
    double t28 = t14 * t14;

    g_compts[0][0] = -1.0 + 2.0 * t1 * t8;
    g_compts[0][1] = 0.0;
    g_compts[0][2] = 0.0;
    g_compts[0][3] = -t18;
    g_compts[1][0] = 0.0;
    g_compts[1][1] = t21 * t2 + t21 * t3 * t5;
    g_compts[1][2] = 0.0;
    g_compts[1][3] = 0.0;
    g_compts[2][0] = 0.0;
    g_compts[2][1] = 0.0;
    g_compts[2][2] = t7;
    g_compts[2][3] = 0.0;
    g_compts[3][0] = -t18;
    g_compts[3][1] = 0.0;
    g_compts[3][2] = 0.0;
    g_compts[3][3] = t14 * t2 + t14 * t3 + 2.0 * t28 * m * t3 * r * t8;

    return true;
}

/*! Calculate the Christoffel symbols of the second kind at position 'pos'.
 *
 *  \param pos : pointer to position.
 */
bool MetricKerrBL::calculateChristoffels(const double* pos) {
    double r     = pos[1];
    double theta = pos[2];
    double m = mass;
    double a = angmom;

    double t1 = r * r;
    double t2 = m * r;
    double t4 = a * a;
    double t5 = t1 - 2.0 * t2 + t4;
    double t6 = cos(theta);
    double t7 = t6 * t6;
    double t8 = t4 * t7;
    double t9 = t1 + t8;
    double t10 = t9 * t9;
    double t12 = 1 / t10 / t9;
    double t14 = -t1 + t8;
    double t20 = sin(theta);
    double t21 = t4 * t6 * t20;
    double t24 = t4 + t1;
    double t26 = 1 / t9;
    double t28 = t4 * t4;
    double t29 = t28 * t7;
    double t30 = t1 * t4;
    double t31 = t30 * t7;
    double t32 = t4 * m;
    double t35 = 2.0 * t32 * r * t7;
    double t36 = t1 * t1;
    double t37 = t1 * r;
    double t38 = m * t37;
    double t41 = 1 / (t29 + t31 + t30 - t35 + t36 - 2.0 * t38);
    double t42 = t14 * t26 * t41;
    double t43 = t24 * m * t42;
    double t44 = m * a;
    double t45 = t44 * t42;
    double t46 = 1 / t10;
    double t49 = 2.0 * t2 * t46 * t21;
    double t50 = t44 * r;
    double t51 = 1 / t20;
    double t55 = 2.0 * t50 * t6 * t51 * t46;
    double t58 = -1.0 + t7;
    double t61 = t5 * m * a * t58 * t14 * t12;
    double t62 = t20 * t6;
    double t66 = 2.0 * t50 * t62 * t24 * t12;
    double t68 = 1 / t5 * t26;
    double t69 = m * t1;
    double t77 = t26 * t4 * t62;
    double t78 = t26 * r;
    double t85 = (t29 - t31 - t30 - 3.0 * t36) * m * a * t58 * t26 * t41;
    double t87 = t7 * t7;
    double t88 = r * t28 * t87;
    double t89 = m * t28;
    double t90 = t89 * t7;
    double t91 = t89 * t87;
    double t92 = t69 * t8;
    double t95 = 2.0 * t37 * t4 * t7;
    double t96 = t32 * t1;
    double t97 = t36 * r;
    double t102 = (t88 + t90 - t91 - t92 + t95 - t96 + t97 - 2.0 * m * t36) * t26 * t41;
    double t112 = 2.0 * t20 * t58 * m * t4 * a * r * t6 * t46;
    double t113 = t87 * t28;
    double t120 = t6 * (t113 + 2.0 * t32 * r - t35 + 2.0 * t31 + t36) * t51 * t46;
    double t126 = t36 * t4;
    double t129 = t1 * t28;

    christoffel[0][0][0] = 0.0;
    christoffel[0][0][1] = -t5 * t12 * m * t14;
    christoffel[0][0][2] = -2.0 * t12 * m * r * t21;
    christoffel[0][0][3] = 0.0;
    christoffel[0][1][0] = -t43;
    christoffel[0][1][1] = 0.0;
    christoffel[0][1][2] = 0.0;
    christoffel[0][1][3] = -t45;
    christoffel[0][2][0] = -t49;
    christoffel[0][2][1] = 0.0;
    christoffel[0][2][2] = 0.0;
    christoffel[0][2][3] = -t55;
    christoffel[0][3][0] = 0.0;
    christoffel[0][3][1] = -t61;
    christoffel[0][3][2] = t66;
    christoffel[0][3][3] = 0.0;
    christoffel[1][0][0] = -t43;
    christoffel[1][0][1] = 0.0;
    christoffel[1][0][2] = 0.0;
    christoffel[1][0][3] = -t45;
    christoffel[1][1][0] = 0.0;
    christoffel[1][1][1] = -t68 * (t69 - r * t4 + t8 * r - t8 * m);
    christoffel[1][1][2] = t68 * t21;
    christoffel[1][1][3] = 0.0;
    christoffel[1][2][0] = 0.0;
    christoffel[1][2][1] = -t77;
    christoffel[1][2][2] = t78;
    christoffel[1][2][3] = 0.0;
    christoffel[1][3][0] = -t85;
    christoffel[1][3][1] = 0.0;
    christoffel[1][3][2] = 0.0;
    christoffel[1][3][3] = t102;
    christoffel[2][0][0] = -t49;
    christoffel[2][0][1] = 0.0;
    christoffel[2][0][2] = 0.0;
    christoffel[2][0][3] = -t55;
    christoffel[2][1][0] = 0.0;
    christoffel[2][1][1] = -t77;
    christoffel[2][1][2] = t78;
    christoffel[2][1][3] = 0.0;
    christoffel[2][2][0] = 0.0;
    christoffel[2][2][1] = -t5 * t26 * r;
    christoffel[2][2][2] = -t77;
    christoffel[2][2][3] = 0.0;
    christoffel[2][3][0] = -t112;
    christoffel[2][3][1] = 0.0;
    christoffel[2][3][2] = 0.0;
    christoffel[2][3][3] = t120;
    christoffel[3][0][0] = 0.0;
    christoffel[3][0][1] = -t61;
    christoffel[3][0][2] = t66;
    christoffel[3][0][3] = 0.0;
    christoffel[3][1][0] = -t85;
    christoffel[3][1][1] = 0.0;
    christoffel[3][1][2] = 0.0;
    christoffel[3][1][3] = t102;
    christoffel[3][2][0] = -t112;
    christoffel[3][2][1] = 0.0;
    christoffel[3][2][2] = 0.0;
    christoffel[3][2][3] = t120;
    christoffel[3][3][0] = 0.0;
    christoffel[3][3][1] = t5 * t58 * (t97 + t95 + t88 - t96 + t92 + t90 - t91) * t12;
    christoffel[3][3][2] = -t62 * (t36 * t1 + 2.0 * t126 * t7 + t129 * t87 + t126 + 2.0 * t129 * t7 + t28 * t4 * t87 + 4.0 * t38 * t4 - 4.0 * t38 * t8 - 2.0 * t2 * t113 + 2.0 * t89 * r) * t12;
    christoffel[3][3][3] = 0.0;

    return true;
}

/*! Calculate Jacobi matrix.
 *
 *  \param pos : pointer to position.
 */
bool MetricKerrBL::calculateChrisD(const double* pos) {
    double r     = pos[1];
    double theta = pos[2];
    double m = mass;
    double a = angmom;

    double t1 = r * r;
    double t2 = t1 * t1;
    double t3 = t2 * r;
    double t4 = a * a;
    double t5 = t4 * t4;
    double t6 = r * t5;
    double t7 = cos(theta);
    double t8 = t7 * t7;
    double t9 = t8 * t8;
    double t10 = t6 * t9;
    double t11 = m * t2;
    double t12 = 3.0 * t11;
    double t13 = m * t5;
    double t14 = t13 * t9;
    double t15 = t1 * r;
    double t16 = t15 * t4;
    double t17 = t16 * t8;
    double t19 = m * t1;
    double t20 = t4 * t8;
    double t21 = t19 * t20;
    double t24 = t6 * t8;
    double t26 = t3 + t10 - t12 - t14 - 4.0 * t17 + 8.0 * t21 + 2.0 * t16 - 4.0 * t24;
    double t28 = t1 + t20;
    double t29 = t28 * t28;
    double t30 = t29 * t29;
    double t31 = 1 / t30;
    double t34 = m * r;
    double t36 = t1 - 2.0 * t34 + t4;
    double t37 = t36 * m;
    double t39 = sin(theta);
    double t40 = t7 * t39;
    double t47 = m * t4;
    double t56 = t9 * t4;
    double t57 = 4.0 * t56;
    double t58 = t8 * t1;
    double t59 = 2.0 * t58;
    double t60 = -5.0 * t20 + t57 + t1 - t59;
    double t65 = t5 * t4;
    double t66 = t65 * t9;
    double t67 = t66 * m;
    double t68 = t65 * t8;
    double t70 = 3.0 * t68 * r;
    double t71 = t5 * t1;
    double t72 = t9 * m;
    double t73 = t71 * t72;
    double t74 = m * t8;
    double t75 = t71 * t74;
    double t76 = 6.0 * t75;
    double t77 = t5 * t15;
    double t78 = t77 * t8;
    double t80 = t11 * t20;
    double t81 = 6.0 * t80;
    double t82 = t3 * t4;
    double t83 = t82 * t8;
    double t85 = t4 * t2;
    double t86 = t85 * m;
    double t87 = 3.0 * t86;
    double t88 = 2.0 * t82;
    double t89 = t2 * t1;
    double t90 = t89 * m;
    double t91 = t2 * t15;
    double t92 = t67 - t70 - t73 + t76 - 6.0 * t78 + t77 + t81 - 3.0 * t83 - t87 + t88 - t90 + t91;
    double t94 = 1 / t28;
    double t95 = t5 * t8;
    double t96 = t1 * t4;
    double t97 = t96 * t8;
    double t99 = t47 * r * t8;
    double t101 = m * t15;
    double t103 = t95 + t97 + t96 - 2.0 * t99 + t2 - 2.0 * t101;
    double t104 = t103 * t103;
    double t105 = 1 / t104;
    double t108 = 2.0 * t92 * m * t94 * t105;
    double t110 = t4 + t1;
    double t112 = t20 - 3.0 * t1;
    double t114 = 1 / t29;
    double t116 = t114 / t103;
    double t119 = 2.0 * t40 * t47 * t110 * t112 * t116;
    double t122 = 6.0 * t21;
    double t129 = 2.0 * (-t14 + t10 + 3.0 * t24 + 3.0 * t17 - t122 - t16 + t12 - 2.0 * t3) * m * a * t94 * t105;
    double t130 = t4 * a;
    double t132 = m * t112;
    double t135 = 2.0 * t40 * t130 * t132 * t116;
    double t137 = t39 * t4 * t7;
    double t139 = 1 / t29 / t28;
    double t142 = 2.0 * t137 * t132 * t139;
    double t143 = 3.0 * t20;
    double t149 = 2.0 * t47 * r * (-t143 + 2.0 * t56 + t1 - t59) * t139;
    double t150 = a * m;
    double t156 = 2.0 * t150 * t7 * t112 / t39 * t139;
    double t157 = t150 * r;
    double t160 = -1.0 + t8;
    double t161 = 1 / t160;
    double t164 = 2.0 * t157 * (t1 - t143 + t57) * t139 * t161;
    double t168 = 2.0 * t150 * t160 * t26 * t31;
    double t170 = t9 * t5;
    double t178 = 2.0 * t37 * a * t40 * (t2 + t170 + 4.0 * t96 - 2.0 * t95 - 4.0 * t97) * t31;
    double t180 = 3.0 * t2;
    double t181 = 3.0 * t97;
    double t187 = 2.0 * t150 * t39 * t7 * (-t180 + t181 - 5.0 * t96 + t95) * t31;
    double t191 = 2.0 * t157 * t110 * t60 * t31;
    double t192 = m * m;
    double t193 = t192 * t1;
    double t196 = t34 * t95;
    double t200 = t71 * t9;
    double t201 = t85 * t8;
    double t202 = 3.0 * t201;
    double t204 = t34 * t170;
    double t206 = t101 * t20;
    double t208 = t101 * t4;
    double t214 = -8.0 * t193 * t20 + 4.0 * t196 - t68 + 2.0 * t192 * t2 + t71 + t66 - t200 - t202 + 3.0 * t85 + 2.0 * t204 + 8.0 * t206 - 4.0 * t208 - 2.0 * m * t3 - 2.0 * t170 * t192;
    double t215 = t36 * t36;
    double t216 = 1 / t215;
    double t223 = 2.0 * t4 * t7 * r * t39 * t114;
    double t235 = t4 * (t20 - t1 + t59);
    double t239 = t235 * t114;
    double t241 = (-t1 + t20) * t114;
    double t243 = r * t65 * t9;
    double t248 = -t67 + 2.0 * t243 + t70 - t73 + 8.0 * t78 - t76 - t77 - t81 + t83 + t87 - t88 + 3.0 * t90 - 3.0 * t91;
    double t254 = 2.0 * t248 * m * a * t160 * t94 * t105;
    double t265 = 2.0 * t40 * a * m * (t68 + 2.0 * t71 * t8 - 3.0 * t71 + t201 - 6.0 * t85 - 3.0 * t89) * t116;
    double t266 = t5 * t5;
    double t267 = t9 * t8;
    double t268 = t266 * t267;
    double t269 = t65 * t267;
    double t270 = t269 * t1;
    double t272 = 2.0 * t269 * t192;
    double t273 = t269 * t34;
    double t274 = 2.0 * t273;
    double t276 = t1 * t65 * t9;
    double t279 = 2.0 * t66 * t192;
    double t280 = t68 * t34;
    double t282 = t2 * t5;
    double t283 = t282 * t9;
    double t285 = t170 * t101;
    double t287 = t170 * t193;
    double t289 = t95 * t193;
    double t291 = t268 - t270 - t272 + t274 + 3.0 * t276 + t279 - 6.0 * t280 - 3.0 * t283 + 6.0 * t285 - 8.0 * t287 + 12.0 * t289;
    double t292 = t77 * t74;
    double t293 = 12.0 * t292;
    double t294 = t282 * t8;
    double t296 = t13 * t15;
    double t297 = 2.0 * t296;
    double t298 = t89 * t4;
    double t299 = t298 * t8;
    double t303 = 6.0 * t85 * t8 * t192;
    double t306 = 6.0 * t4 * t192 * t2;
    double t307 = t82 * m;
    double t308 = 4.0 * t307;
    double t309 = t2 * t2;
    double t314 = -t293 + 3.0 * t294 + t297 - 3.0 * t299 + t303 + t298 - t306 + t308 - t309 + 4.0 * t91 * m - 4.0 * t192 * t89;
    double t317 = (t291 + t314) * t94 * t105;
    double t323 = 2.0 * t137 * m * (t95 + t97 - 3.0 * t96 - t180) * t116;
    double t337 = 2.0 * t39 * t160 * m * t130 * t7 * t112 * t139;
    double t347 = 2.0 * m * t130 * r * (-5.0 * t58 + 3.0 * t56 + 4.0 * t1 * t9 + t1 - t143) * t139;
    double t348 = t15 * t9;
    double t351 = r * t267;
    double t361 = (4.0 * t348 * t47 - 4.0 * t13 * t351 + 10.0 * t204 - 6.0 * t206 + 2.0 * t208 + t202 + t89 + 3.0 * t200 - 6.0 * t196 + t269) * t139 * t161;
    double t362 = t82 * t74;
    double t364 = t66 * t34;
    double t371 = 3.0 * t270;
    double t376 = 2.0 * t307 + t371 + t272 + t268 - t279 - 16.0 * t287 + 16.0 * t289 + t303 - 6.0 * t273 - 8.0 * t280 - t298;
    double t393 = t82 + t91 + t243 + t77 * t9 + 2.0 * t80 - 4.0 * t73 - 2.0 * t67 + 2.0 * t83 + 8.0 * t75 + 2.0 * t78 - 2.0 * t86 + 2.0 * t68 * m - 4.0 * t13 * t1;
    double t412 = t8 * t3 + 2.0 * t348 * t4 + t351 * t5 - t3 - 2.0 * t17 - t10 - 6.0 * t11 + t122 + 6.0 * t74 * t2 - t72 * t96 - t267 * t5 * m - 5.0 * t47 * t1 + t13 * t8;
    double t438 = 24.0 * t292 - 18.0 * t285 + 4.0 * t9 * t89 * t4 + t299 + 5.0 * t283 + t371 + 2.0 * t267 * t2 * t5 + t268 - t297 - t308 - t298;

    chrisD[0][0][0][0] = 0.0;
    chrisD[0][0][0][1] = 0.0;
    chrisD[0][0][0][2] = 0.0;
    chrisD[0][0][0][3] = 0.0;
    chrisD[0][0][1][0] = 0.0;
    chrisD[0][0][1][1] = -2.0 * m * t26 * t31;
    chrisD[0][0][1][2] = -4.0 * t37 * t4 * t40 * (-2.0 * t1 + t20) * t31;
    chrisD[0][0][1][3] = 0.0;
    chrisD[0][0][2][0] = 0.0;
    chrisD[0][0][2][1] = -2.0 * t47 * t7 * t39 * (-5.0 * t1 + t20) * t31;
    chrisD[0][0][2][2] = 2.0 * t47 * r * t60 * t31;
    chrisD[0][0][2][3] = 0.0;
    chrisD[0][0][3][0] = 0.0;
    chrisD[0][0][3][1] = 0.0;
    chrisD[0][0][3][2] = 0.0;
    chrisD[0][0][3][3] = 0.0;
    chrisD[0][1][0][0] = 0.0;
    chrisD[0][1][0][1] = -t108;
    chrisD[0][1][0][2] = -t119;
    chrisD[0][1][0][3] = 0.0;
    chrisD[0][1][1][0] = 0.0;
    chrisD[0][1][1][1] = 0.0;
    chrisD[0][1][1][2] = 0.0;
    chrisD[0][1][1][3] = 0.0;
    chrisD[0][1][2][0] = 0.0;
    chrisD[0][1][2][1] = 0.0;
    chrisD[0][1][2][2] = 0.0;
    chrisD[0][1][2][3] = 0.0;
    chrisD[0][1][3][0] = 0.0;
    chrisD[0][1][3][1] = t129;
    chrisD[0][1][3][2] = -t135;
    chrisD[0][1][3][3] = 0.0;
    chrisD[0][2][0][0] = 0.0;
    chrisD[0][2][0][1] = -t142;
    chrisD[0][2][0][2] = t149;
    chrisD[0][2][0][3] = 0.0;
    chrisD[0][2][1][0] = 0.0;
    chrisD[0][2][1][1] = 0.0;
    chrisD[0][2][1][2] = 0.0;
    chrisD[0][2][1][3] = 0.0;
    chrisD[0][2][2][0] = 0.0;
    chrisD[0][2][2][1] = 0.0;
    chrisD[0][2][2][2] = 0.0;
    chrisD[0][2][2][3] = 0.0;
    chrisD[0][2][3][0] = 0.0;
    chrisD[0][2][3][1] = -t156;
    chrisD[0][2][3][2] = -t164;
    chrisD[0][2][3][3] = 0.0;
    chrisD[0][3][0][0] = 0.0;
    chrisD[0][3][0][1] = 0.0;
    chrisD[0][3][0][2] = 0.0;
    chrisD[0][3][0][3] = 0.0;
    chrisD[0][3][1][0] = 0.0;
    chrisD[0][3][1][1] = -t168;
    chrisD[0][3][1][2] = -t178;
    chrisD[0][3][1][3] = 0.0;
    chrisD[0][3][2][0] = 0.0;
    chrisD[0][3][2][1] = t187;
    chrisD[0][3][2][2] = -t191;
    chrisD[0][3][2][3] = 0.0;
    chrisD[0][3][3][0] = 0.0;
    chrisD[0][3][3][1] = 0.0;
    chrisD[0][3][3][2] = 0.0;
    chrisD[0][3][3][3] = 0.0;
    chrisD[1][0][0][0] = 0.0;
    chrisD[1][0][0][1] = -t108;
    chrisD[1][0][0][2] = -t119;
    chrisD[1][0][0][3] = 0.0;
    chrisD[1][0][1][0] = 0.0;
    chrisD[1][0][1][1] = 0.0;
    chrisD[1][0][1][2] = 0.0;
    chrisD[1][0][1][3] = 0.0;
    chrisD[1][0][2][0] = 0.0;
    chrisD[1][0][2][1] = 0.0;
    chrisD[1][0][2][2] = 0.0;
    chrisD[1][0][2][3] = 0.0;
    chrisD[1][0][3][0] = 0.0;
    chrisD[1][0][3][1] = t129;
    chrisD[1][0][3][2] = -t135;
    chrisD[1][0][3][3] = 0.0;
    chrisD[1][1][0][0] = 0.0;
    chrisD[1][1][0][1] = 0.0;
    chrisD[1][1][0][2] = 0.0;
    chrisD[1][1][0][3] = 0.0;
    chrisD[1][1][1][0] = 0.0;
    chrisD[1][1][1][1] = -t214 * t216 * t114;
    chrisD[1][1][1][2] = t223;
    chrisD[1][1][1][3] = 0.0;
    chrisD[1][1][2][0] = 0.0;
    chrisD[1][1][2][1] = -2.0 * t137 * (t20 * r + r * t4 - t20 * m + 2.0 * t15 - 3.0 * t19) * t216 * t114;
    chrisD[1][1][2][2] = t235 / t36 * t114;
    chrisD[1][1][2][3] = 0.0;
    chrisD[1][1][3][0] = 0.0;
    chrisD[1][1][3][1] = 0.0;
    chrisD[1][1][3][2] = 0.0;
    chrisD[1][1][3][3] = 0.0;
    chrisD[1][2][0][0] = 0.0;
    chrisD[1][2][0][1] = 0.0;
    chrisD[1][2][0][2] = 0.0;
    chrisD[1][2][0][3] = 0.0;
    chrisD[1][2][1][0] = 0.0;
    chrisD[1][2][1][1] = t223;
    chrisD[1][2][1][2] = -t239;
    chrisD[1][2][1][3] = 0.0;
    chrisD[1][2][2][0] = 0.0;
    chrisD[1][2][2][1] = t241;
    chrisD[1][2][2][2] = t223;
    chrisD[1][2][2][3] = 0.0;
    chrisD[1][2][3][0] = 0.0;
    chrisD[1][2][3][1] = 0.0;
    chrisD[1][2][3][2] = 0.0;
    chrisD[1][2][3][3] = 0.0;
    chrisD[1][3][0][0] = 0.0;
    chrisD[1][3][0][1] = t254;
    chrisD[1][3][0][2] = t265;
    chrisD[1][3][0][3] = 0.0;
    chrisD[1][3][1][0] = 0.0;
    chrisD[1][3][1][1] = 0.0;
    chrisD[1][3][1][2] = 0.0;
    chrisD[1][3][1][3] = 0.0;
    chrisD[1][3][2][0] = 0.0;
    chrisD[1][3][2][1] = 0.0;
    chrisD[1][3][2][2] = 0.0;
    chrisD[1][3][2][3] = 0.0;
    chrisD[1][3][3][0] = 0.0;
    chrisD[1][3][3][1] = t317;
    chrisD[1][3][3][2] = t323;
    chrisD[1][3][3][3] = 0.0;
    chrisD[2][0][0][0] = 0.0;
    chrisD[2][0][0][1] = -t142;
    chrisD[2][0][0][2] = t149;
    chrisD[2][0][0][3] = 0.0;
    chrisD[2][0][1][0] = 0.0;
    chrisD[2][0][1][1] = 0.0;
    chrisD[2][0][1][2] = 0.0;
    chrisD[2][0][1][3] = 0.0;
    chrisD[2][0][2][0] = 0.0;
    chrisD[2][0][2][1] = 0.0;
    chrisD[2][0][2][2] = 0.0;
    chrisD[2][0][2][3] = 0.0;
    chrisD[2][0][3][0] = 0.0;
    chrisD[2][0][3][1] = -t156;
    chrisD[2][0][3][2] = -t164;
    chrisD[2][0][3][3] = 0.0;
    chrisD[2][1][0][0] = 0.0;
    chrisD[2][1][0][1] = 0.0;
    chrisD[2][1][0][2] = 0.0;
    chrisD[2][1][0][3] = 0.0;
    chrisD[2][1][1][0] = 0.0;
    chrisD[2][1][1][1] = t223;
    chrisD[2][1][1][2] = -t239;
    chrisD[2][1][1][3] = 0.0;
    chrisD[2][1][2][0] = 0.0;
    chrisD[2][1][2][1] = t241;
    chrisD[2][1][2][2] = t223;
    chrisD[2][1][2][3] = 0.0;
    chrisD[2][1][3][0] = 0.0;
    chrisD[2][1][3][1] = 0.0;
    chrisD[2][1][3][2] = 0.0;
    chrisD[2][1][3][3] = 0.0;
    chrisD[2][2][0][0] = 0.0;
    chrisD[2][2][0][1] = 0.0;
    chrisD[2][2][0][2] = 0.0;
    chrisD[2][2][0][3] = 0.0;
    chrisD[2][2][1][0] = 0.0;
    chrisD[2][2][1][1] = -(t2 + t181 - 4.0 * t99 - t96 + t95) * t114;
    chrisD[2][2][1][2] = -2.0 * t36 * t114 * r * t137;
    chrisD[2][2][1][3] = 0.0;
    chrisD[2][2][2][0] = 0.0;
    chrisD[2][2][2][1] = t223;
    chrisD[2][2][2][2] = -t239;
    chrisD[2][2][2][3] = 0.0;
    chrisD[2][2][3][0] = 0.0;
    chrisD[2][2][3][1] = 0.0;
    chrisD[2][2][3][2] = 0.0;
    chrisD[2][2][3][3] = 0.0;
    chrisD[2][3][0][0] = 0.0;
    chrisD[2][3][0][1] = -t337;
    chrisD[2][3][0][2] = -t347;
    chrisD[2][3][0][3] = 0.0;
    chrisD[2][3][1][0] = 0.0;
    chrisD[2][3][1][1] = 0.0;
    chrisD[2][3][1][2] = 0.0;
    chrisD[2][3][1][3] = 0.0;
    chrisD[2][3][2][0] = 0.0;
    chrisD[2][3][2][1] = 0.0;
    chrisD[2][3][2][2] = 0.0;
    chrisD[2][3][2][3] = 0.0;
    chrisD[2][3][3][0] = 0.0;
    chrisD[2][3][3][1] = t142;
    chrisD[2][3][3][2] = t361;
    chrisD[2][3][3][3] = 0.0;
    chrisD[3][0][0][0] = 0.0;
    chrisD[3][0][0][1] = 0.0;
    chrisD[3][0][0][2] = 0.0;
    chrisD[3][0][0][3] = 0.0;
    chrisD[3][0][1][0] = 0.0;
    chrisD[3][0][1][1] = -t168;
    chrisD[3][0][1][2] = -t178;
    chrisD[3][0][1][3] = 0.0;
    chrisD[3][0][2][0] = 0.0;
    chrisD[3][0][2][1] = t187;
    chrisD[3][0][2][2] = -t191;
    chrisD[3][0][2][3] = 0.0;
    chrisD[3][0][3][0] = 0.0;
    chrisD[3][0][3][1] = 0.0;
    chrisD[3][0][3][2] = 0.0;
    chrisD[3][0][3][3] = 0.0;
    chrisD[3][1][0][0] = 0.0;
    chrisD[3][1][0][1] = t254;
    chrisD[3][1][0][2] = t265;
    chrisD[3][1][0][3] = 0.0;
    chrisD[3][1][1][0] = 0.0;
    chrisD[3][1][1][1] = 0.0;
    chrisD[3][1][1][2] = 0.0;
    chrisD[3][1][1][3] = 0.0;
    chrisD[3][1][2][0] = 0.0;
    chrisD[3][1][2][1] = 0.0;
    chrisD[3][1][2][2] = 0.0;
    chrisD[3][1][2][3] = 0.0;
    chrisD[3][1][3][0] = 0.0;
    chrisD[3][1][3][1] = t317;
    chrisD[3][1][3][2] = t323;
    chrisD[3][1][3][3] = 0.0;
    chrisD[3][2][0][0] = 0.0;
    chrisD[3][2][0][1] = -t337;
    chrisD[3][2][0][2] = -t347;
    chrisD[3][2][0][3] = 0.0;
    chrisD[3][2][1][0] = 0.0;
    chrisD[3][2][1][1] = 0.0;
    chrisD[3][2][1][2] = 0.0;
    chrisD[3][2][1][3] = 0.0;
    chrisD[3][2][2][0] = 0.0;
    chrisD[3][2][2][1] = 0.0;
    chrisD[3][2][2][2] = 0.0;
    chrisD[3][2][2][3] = 0.0;
    chrisD[3][2][3][0] = 0.0;
    chrisD[3][2][3][1] = t142;
    chrisD[3][2][3][2] = t361;
    chrisD[3][2][3][3] = 0.0;
    chrisD[3][3][0][0] = 0.0;
    chrisD[3][3][0][1] = 0.0;
    chrisD[3][3][0][2] = 0.0;
    chrisD[3][3][0][3] = 0.0;
    chrisD[3][3][1][0] = 0.0;
    chrisD[3][3][1][1] = t160 * (t309 - t293 - 6.0 * t362 + 10.0 * t364 + t276 + 5.0 * t299 + 7.0 * t283 - t294 + 4.0 * t296 - t306 + t376) * t31;
    chrisD[3][3][1][2] = -2.0 * t36 * t7 * t39 * t393 * t31;
    chrisD[3][3][1][3] = 0.0;
    chrisD[3][3][2][0] = 0.0;
    chrisD[3][3][2][1] = -2.0 * t40 * t4 * t412 * t31;
    chrisD[3][3][2][2] = -(-t309 + 2.0 * t8 * t309 - t274 + 10.0 * t280 - t294 + t276 - 16.0 * t9 * t3 * t47 - 4.0 * t267 * t15 * t13 + 20.0 * t362 - 8.0 * t364 + t438) * t31;
    chrisD[3][3][2][3] = 0.0;
    chrisD[3][3][3][0] = 0.0;
    chrisD[3][3][3][1] = 0.0;
    chrisD[3][3][3][2] = 0.0;
    chrisD[3][3][3][3] = 0.0;

    return true;
}


/*! Transform local 4-direction to coordinate 4-direction.
 *
 *  \param  pos  :  pointer to position array.
 *  \param  ldir :  pointer to local direction array.
 *  \param  dir  :  pointer to calculated coordinate direction array.
 *  \param  type :  type of tetrad.
 */
void MetricKerrBL::localToCoord(const double* pos, const double* ldir, double* dir,
                                enum_nat_tetrad_type  type) {
    double r     = pos[1];
    double theta = pos[2];

    calcSigmaAndDelta(r, theta);
    calcA(r, theta);
    double rs = 2.0 * mass;

    if (type == enum_nat_tetrad_static) {
        double w = sqrt(1.0 - rs * r / sigma);
        /*
        dir[0] = ldir[0]/w + ldir[3]*pm*rs*angmom*r*sin(theta)/(w*sqrt(delta)*sigma);
        dir[1] = ldir[1]*sqrt(delta/sigma);
        dir[2] = ldir[2]/sqrt(sigma);
        dir[3] = - ldir[3]*pm*w/sqrt(delta)/sin(theta);
        */
        double e00 = 1.0 / w;
        double e11 = sqrt(delta / sigma);
        double e22 = 1.0 / sqrt(sigma);
        double e30 = rs * angmom * r * sin(theta) / (w * sqrt(delta) * sigma);
        double e33 = -w / (sqrt(delta) * sin(theta));

        dir[0] = ldir[0] * e00 + ldir[3] * e30;
        dir[1] = ldir[1] * e11;
        dir[2] = ldir[2] * e22;
        dir[3] = ldir[3] * e33;
    } else {
        double omega = rs * angmom * r / A;
        ga = sqrt(A / (sigma * delta));
        dir[0] = ga * ldir[0];
        dir[1] = ldir[1] * sqrt(delta / sigma);
        dir[2] = ldir[2] / sqrt(sigma);
        dir[3] = ldir[0] * omega * ga + ldir[3] * sqrt(sigma / A) / sin(theta);
    }
}

/*! Transform coordinate 4-direction to local 4-direction.
 *
 *  \param  pos  :  pointer to position array.
 *  \param  cdir :  pointer to coordinate direction.
 *  \param  ldir :  pointer to calculated local direction array.
 *  \param  type :  type of tetrad.
 */
void MetricKerrBL::coordToLocal(const double* pos, const double* cdir, double* ldir,
                                enum_nat_tetrad_type  type) {
    double r     = pos[1];
    double theta = pos[2];

    calcSigmaAndDelta(r, theta);
    calcA(r, theta);
    double rs = 2.0 * mass;

    if (type == enum_nat_tetrad_static) {
        double w = sqrt(1.0 - rs * r / sigma);
        double st = sin(theta);
        /*
        ldir[0] = w*cdir[0] + rs*angmom*r*sin(theta)*sin(theta)/(w*sigma)*cdir[3];
        ldir[1] = sqrt(sigma/delta)*cdir[1];
        ldir[2] = sqrt(sigma)*cdir[2];
        ldir[3] = - pm/w*sqrt(delta)*sin(theta)*cdir[3];
        */
        double t00 = w;
        double t03 = rs * angmom * r * st * st / (sigma * w);
        double t11 = sqrt(sigma / delta);
        double t22 = sqrt(sigma);
        double t33 = -sqrt(delta) * st / w;

        ldir[0] = t00 * cdir[0] + t03 * cdir[3];
        ldir[1] = t11 * cdir[1];
        ldir[2] = t22 * cdir[2];
        ldir[3] = t33 * cdir[3];
    } else {
        double omega = rs * angmom * r / A;
        ldir[0] = sqrt(delta * sigma / A) * cdir[0];
        ldir[1] = sqrt(sigma / delta) * cdir[1];
        ldir[2] = sqrt(sigma) * cdir[2];
        ldir[3] = sqrt(A / sigma) * sin(theta) * (cdir[3] - omega * cdir[0]);
    }
}


/*!
 *  \param pos  :  position.
 *  \return true  : radial position r < 0.0 or ...
 *  \return false : position is valid.
 */
bool MetricKerrBL::breakCondition(const double* pos) {
    bool br = false;

    double r = pos[1];

    double eh = mass + sqrt(mass * mass - angmom * angmom);
    if ((r < 0.0) || (r <= (1.0 + eps)*eh)) {
        br = true;
    }
    return br;
}


/*! Calculate right hand side of the geodesic equation in first order form.
 *
 *  \param  y[]   : pointer to position and direction coordinates.
 *  \param  dydx[] : pointer to right side of geodesic equation.
 */
bool MetricKerrBL::calcDerivs(const double y[], double dydx[]) {
    dydx[0] = y[4];
    dydx[1] = y[5];
    dydx[2] = y[6];
    dydx[3] = y[7];

    double r = y[1];
    double theta = y[2];
    double rs = 2.0 * mass;

    double st = sin(theta);
    double ct = cos(theta);
    double st2 = st * st;
    double ct2 = ct * ct;

    double r2 = r * r;
    double a2 = angmom * angmom;
    double Sigma = r2 + a2 * ct2;
    double Delta = r2 - rs * r + a2;
    double A = (r2 + a2) * Sigma + rs * r * a2 * st2;

    double w = r2 - a2 * ct2;
    double S2 = Sigma * Sigma;
    double S3 = Sigma * S2;

    double G_t_t_r    = 0.5 * rs * Delta * w / S3;
    double G_t_r_t    = 0.5 * rs * (r2 + a2) * w / S2 / Delta;
    double G_t_th_t   = -rs * r * a2 * st * ct / S2;
    double G_t_ph_r   = -0.5 * Delta * rs * angmom * st2 * w / S3;
    double G_r_r_r    = 0.5 * (2.0 * r * a2 * st2 - rs * w) / Sigma / Delta;
    double G_r_th_r   = -a2 * st * ct / Sigma;
    double G_r_ph_t   = 0.5 * rs * angmom * st2 * (a2 * ct2 * (a2 - r2) - r2 * (a2 + 3.0 * r2)) / S2 / Delta;
    double G_r_ph_ph  = 0.5 * (2.0 * r * S2 + rs * (a2 * a2 * st2 * ct2 - r2 * (Sigma + r2 + a2))) / S2 / Delta;
    double G_th_ph_ph = ct / st / S2 * (S2 + rs * a2 * r * st2);
    double G_ph_ph_r  = 0.5 * Delta * st2 / S3 * (-2.0 * r * S2 + rs * a2 * st2 * w);
    double G_ph_ph_th = -st * ct / S3 * (A * Sigma + (r2 + a2) * rs * a2 * r * st2);
    double G_t_t_th   = -rs * r * a2 * st * ct / S3;
    double G_t_r_ph   = 0.5 * rs * angmom * w / S2 / Delta;
    double G_t_th_ph  = -rs * r * angmom * ct / st / S2;
    double G_t_ph_th  = rs * r * angmom * (r2 + a2) * st * ct / S3;
    double G_r_r_th   = a2 * st * ct / Sigma / Delta;
    double G_r_th_th  = r / Sigma;
    double G_th_th_r  = -r * Delta / Sigma;
    double G_th_th_th = -a2 * st * ct / Sigma;
    double G_th_ph_t  = rs * r * a2 * angmom * st2 * st * ct / S2;

    dydx[4] = -2.0 * G_t_r_t*y[4] * y[5] - 2.0 * G_t_th_t*y[4] * y[6] - 2.0 * G_r_ph_t*y[5] * y[7] - 2.0 * G_th_ph_t*y[6] * y[7];
    dydx[5] = -G_t_t_r * y[4] * y[4] - 2.0 * G_t_ph_r * y[4] * y[7] - G_r_r_r * y[5] * y[5] - 2.0 * G_r_th_r * y[5] * y[6] - G_ph_ph_r * y[7] * y[7] - G_th_th_r * y[6] * y[6];
    dydx[6] = -G_ph_ph_th * y[7] * y[7] - G_t_t_th * y[4] * y[4] - 2.0 * G_t_ph_th * y[4] * y[7] - G_r_r_th * y[5] * y[5] - 2.0 * G_r_th_th * y[5] * y[6] - G_th_th_th * y[6] * y[6];
    dydx[7] = -2.0 * G_r_ph_ph * y[5] * y[7] - 2.0 * G_th_ph_ph * y[6] * y[7] - 2.0 * G_t_r_ph * y[4] * y[5] - 2.0 * G_t_th_ph * y[4] * y[6];

    return true;
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
double MetricKerrBL::testConstraint(const double y[], const double kappa) {
    double r = y[1];
    double theta = y[2];

    calcSigmaAndDelta(r, theta);

    double st2 = sin(theta) * sin(theta);
    double sum = -kappa;
    sum += -(1.0 - 2.0 * mass * r / sigma) * y[4] * y[4] - 4.0 * mass * angmom * r * st2 / sigma * y[4] * y[7] + sigma / delta * y[5] * y[5] + sigma * y[6] * y[6] + (r * r + angmom * angmom + 2.0 * mass * angmom * angmom * r * st2 / sigma) * st2 * y[7] * y[7];
    return sum;
}

/*! Set parameter 'pName' to 'val'.
 *
 *  Set 'mass' and 'angmom' parameter.
 */
bool MetricKerrBL::setParam(const char* pName, double val) {
    Metric::setParam(pName, val);
    if (strcmp(pName,"mass") == 0) {
        mass = val;
    }
    else if (strcmp(pName,"angmom") == 0) {
        angmom = val;
    }
    return true;
}

/*!  Transform point p to 2+1 coordinates.
 *
 *  \param  p  : point in proper metric coordinates.
 *  \param  cp : reference to transformed point.
 *  \return true : success.
 */
bool MetricKerrBL::transToTwoPlusOne(vec4 p, vec4 &cp) {
    vec4 tp;
    TransCoordinates::toCartesianCoord(mCoordType, p, tp);
    cp = vec4(tp[0], tp[1], tp[2], tp[0]);
    return true;
}

/*! Effective potential.
 *  \param pos : initial position.
 *  \param cdir : initial four-direction.
 *  \param type : geodesic type.
 *  \param x : abscissa value.
 *  \param val : reference to effective potential value.
 *  \return true : effective potential exists at x.
 */
bool MetricKerrBL::effPotentialValue(const vec4 pos, const vec4 cdir , enum_geodesic_type type, const double x, double &val) {
    double kappa = 0.0;
    if (type == enum_geodesic_timelike) {
        kappa = -mSign;
    }

    double eh = mass + sqrt(mass * mass - angmom * angmom);
    if ((pos[1] <= (1.0 + eps)*eh) || (x <= (1.0 + eps)*eh)) {
        return false;
    }

    double rs = 2.0 * mass;
    double k = (1.0 - rs / pos[1]) * cdir[0] + rs * angmom / pos[1] * cdir[3];
    double h = (pos[1] * pos[1] + angmom * angmom + rs * angmom * angmom / pos[1]) * cdir[3] - rs * angmom / pos[1] * cdir[0];
    //fprintf(stderr,"%f %f  %f %f\n",cdir[0],cdir[3],k,h);
    val = 0.5 * pow(x, -3.0) * (h * h * (x - rs) + 2.0 * angmom * h * k * rs - k * k * (x * x * x + angmom * angmom * (x + rs))) - 0.5 * kappa * (x * x - rs * x + angmom * angmom) / (x * x);
    return true;
}

/*! Total energy.
 *  \param pos : initial position.
 *  \param cdir : initial four-direction.
 *  \param x : abscissa value.
 *  \param val : reference to total energy value.
 *  \return true : effective potential exists at x.
 */
bool MetricKerrBL::totEnergy(const vec4 , const vec4 , const double , double &val) {
    val = 0.0;
    return true;
}

double MetricKerrBL::getCircularVelocity(const double r, const enum_nat_tetrad_type) {
    double rs = 2.0 * mass;
    double a = angmom;
    double a2 = a * a;
    double delta = r * r - rs * r + a2;
    double A = (r * r + a2) * r * r + rs * a2 * r;

    double beta = (rs * a * pow(r, 4.0) * (3.0 * r * r + a2) - A * sqrt(2.0 * pow(r, 7.0) * rs)) / (r * r * sqrt(delta) * (-2.0 * pow(r, 5.0) + rs * a2 * r * r));
    return beta;
}

void MetricKerrBL::getCircularVelocities(const double r, double &b1, double &b2) {
    double rs = 2.0 * mass;
    double a = angmom;
    double a2 = a * a;
    double delta = r * r - rs * r + a2;
    double A = (r * r + a2) * r * r + rs * a2 * r;

    b1 = (rs * a * pow(r, 4.0) * (3.0 * r * r + a2) - A * sqrt(2.0 * pow(r, 7.0) * rs)) / (r * r * sqrt(delta) * (-2.0 * pow(r, 5.0) + rs * a2 * r * r));
    b2 = (rs * a * pow(r, 4.0) * (3.0 * r * r + a2) + A * sqrt(2.0 * pow(r, 7.0) * rs)) / (r * r * sqrt(delta) * (-2.0 * pow(r, 5.0) + rs * a2 * r * r));
}

vec4
MetricKerrBL::getCircularFourVel(const vec4 pos, const enum_nat_tetrad_type) {
    double beta = getCircularVelocity(pos[1]);
    if (beta > 0.0 && beta < 1.0) {
        double gamma = 1.0 / sqrt(1.0 - beta * beta);
        vec4 e0, e1, e2, e3;
        getNatTetrad(pos, e0, e1, e2, e3);
        return mSpeedOfLight * gamma * (e0 + beta * e3);
    }
    return vec4();
}

/*! Generate report.
 */
bool MetricKerrBL::report(const vec4 pos, const vec4 cdir, std::string &text) {
    std::stringstream ss;
    ss << "Report for Kerr metric in Boyer-Lindquist form\n\tcoordinate : (t,r,theta,phi)\n";
    ss << "---------------------------------------------------------------\n";
    ss << "  physical units ................................ no\n";
    ss.precision(DEF_FIXED_REPORT_PRECISION);
    ss.setf(std::ios::fixed);

    double rp, rphD, rphR, rmsD, rmsR;
    double vD, vR;
    calcMSOvelocities(vD, vR);

    if (calcEventHorizon(rp)) {
        ss << "  Event horizon .......................... r_p = " << rp << std::endl;
    }
    if (calcPhotonOrbit(rphD, rphR)) {
        ss << "  Photon orbit (prograd) ............. r_ph(d) = " << rphD << std::endl;
        ss << "  Photon orbit (retrograd) ........... r_ph(r) = " << rphR << std::endl;
    }
    if (calcMarginStabOrbit(rmsD, rmsR)) {
        ss << "  Margin. stabil orbit (prograd) ..... r_ms(d) = " << rmsD << std::endl;
        ss << "                                          beta = " << vD << std::endl;
        ss << "  Margin. stabil orbit (retrograd) ... r_ms(r) = " << rmsR << std::endl;
    }

    double rs = 2.0 * mass;
    double k = (1.0 - rs / pos[1]) * cdir[0] + rs * angmom / pos[1] * cdir[3];
    double h = (pos[1] * pos[1] + angmom * angmom + rs * angmom * angmom / pos[1]) * cdir[3] - rs * angmom / pos[1] * cdir[0];
    ss << "  constant of motion ....................... k = " << k << std::endl;
    ss << "  constant of motion ....................... h = " << h << std::endl;

    double cvel1, cvel2;
    getCircularVelocities(pos[1], cvel1, cvel2);
    if (fabs(cvel1) < 1.0) {
        ss << "  timelike circular velocity .......pro beta = " << cvel1 << std::endl;
    }
    if (fabs(cvel2) < 1.0) {
        ss << "  timelike circular velocity .....retro beta = " << cvel2 << std::endl;
    }
    text = ss.str();
    return true;
}

// *************************** specific  public methods ****************************
/*! Calculate the radius of the event horizon.
 *
 *  \param  rp  : reference to radius.
 */
bool MetricKerrBL::calcEventHorizon(double &rp) {
    if (mass * mass - angmom * angmom < 0.0) {
        return false;
    }

    rp = mass + sqrt(mass * mass - angmom * angmom);
    return true;
}

/*! Calculate the radius of the outer boundary of the ergosphere.
 *
 *  \param theta : angle.
 *  \param  r0  : reference to radius.
 */
bool MetricKerrBL::calcErgosphere(double theta, double &r0) {
    double w = mass * mass - angmom * angmom * cos(theta) * cos(theta);
    if (w < 0.0) {
        return false;
    }

    r0 = mass + sqrt(w);
    return true;
}

/*! Calculate the radii of the photon orbits for direct and retrograd direction.
 *
 *  \param  rphDirect  : reference to direct radius.
 *  \param  rphRetrograd  : reference to retrograd radius.
 */
bool MetricKerrBL::calcPhotonOrbit(double &rphDirect, double &rphRetrograd) {
    if (angmom > mass) {
        return false;
    }

    rphDirect    = 2.0 * mass * (1.0 + cos(2.0 / 3.0 * acos(-angmom / mass)));
    rphRetrograd = 2.0 * mass * (1.0 + cos(2.0 / 3.0 * acos(angmom / mass)));
    return true;
}

/*! Calculate the radii of the marginally stable orbits for direct and retrograd direction.
 *
 *  \param  rmsDirect  : reference to direct radius.
 *  \param  rmsRetrograd  : reference to retrograd radius.
 */
bool MetricKerrBL::calcMarginStabOrbit(double &rmsDirect, double &rmsRetrograd) {
    if (angmom > mass) {
        return false;
    }

    double adm = angmom / mass;
    double Z1 = 1.0 + pow(1.0 - adm * adm, 1.0 / 3.0) * (pow(1.0 + adm, 1.0 / 3.0) + pow(1.0 - adm, 1.0 / 3.0));
    double Z2 = sqrt(3.0 * adm * adm + Z1 * Z1);

    rmsDirect    = mass * (3.0 + Z2 - sqrt((3.0 - Z1) * (3.0 + Z1 + 2.0 * Z2)));
    rmsRetrograd = mass * (3.0 + Z2 + sqrt((3.0 - Z1) * (3.0 + Z1 + 2.0 * Z2)));
    return true;
}

/*! Calculate velocities for the marginally stable orbits.
 *  \param velDirect : reference to direct velocity.
 *  \param velRetrograd : reference to retrograd velocity.
 */
bool MetricKerrBL::calcMSOvelocities(double &velDirect, double &velRetrograd) {
    double rD, rR;
    calcMarginStabOrbit(rD, rR);
    double ut = (rD * sqrt(rD) + angmom * sqrt(mass)) / sqrt(rD) / sqrt(rD * rD - 3.0 * mass * rD + 2.0 * angmom * sqrt(mass) * sqrt(rD));
    double dt = sqrt(((rD * rD + angmom * angmom) * rD * rD + 2.0 * mass * angmom * angmom * rD) / (rD * rD * (rD * rD - 2.0 * mass * rD + angmom * angmom)));
    double g = ut / dt;
    velDirect = sqrt(1.0 - 1.0 / (g * g));

    if (std::isnan(velDirect)) {
        return false;
    }
    velRetrograd = 0.0;  // TODO
    return true;
}

// ********************************* protected methods *****************************
/*!
 */
void MetricKerrBL::setStandardValues() {
    mInitPos[0] = 0.0;
    mInitPos[1] = 6.0 * mass;
    mInitPos[2] = M_PI_2;
    mInitPos[3] = 0.0;
    mInitDir[0] = 1.0;
    mInitDir[1] = 0.0;
    mInitDir[2] = 0.0;

    mCoordNames[0] = std::string("t");
    mCoordNames[1] = std::string("r");
    mCoordNames[2] = std::string("theta");
    mCoordNames[3] = std::string("phi");
}

/*! Calculate parameters sigma and delta
 *  \param  r : radial position.
 *  \param  theta : theta position.
 */
bool MetricKerrBL::calcSigmaAndDelta(double r, double theta) {
    double ct = cos(theta);
    sigma = r * r + angmom * angmom * ct * ct;
    delta = r * r - 2.0 * mass * r + angmom * angmom;
    return true;
}

/*! Calculate parameter A
 *  \param  r : radial position.
 *  \param  theta : theta position.
 */
bool MetricKerrBL::calcA(double r, double theta) {
    double ct = cos(theta);
    double st = sin(theta);
    sigma = r * r + angmom * angmom * ct * ct;
    A = (r * r + angmom * angmom) * sigma + 2.0 * mass * angmom * angmom * r * st * st;
    return true;
}

/*!
 *  Note that sigma and delta have to be calculated in advance.
 *  \param  r : radial position.
 *  \param  theta : theta position.
 *  \param  zeta  : ...
 */
bool MetricKerrBL::calcGamma(double r, double theta, double zeta) {
    double st2 = sin(theta) * sin(theta);
    double ginv2 = 1.0 - 2.0 * mass * r / sigma + 4.0 * mass * angmom * r * st2 / sigma * zeta - (r * r + angmom * angmom + 2.0 * mass * angmom * angmom * r * st2 / sigma) * st2 * zeta * zeta;
    ga = 1.0 / sqrt(ginv2);
    return true;
}

} // end namespace m4d
