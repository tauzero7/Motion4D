// -------------------------------------------------------------------------------
/*
   m4dMetricHartleThorneGB.cpp

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

#include "m4dMetricHartleThorneGB.h"

namespace m4d {


/*! Standard constructor for the Kottler metric.
 *
 * \param  u : Khalatnikov-Lifshitz parameter.
 */
MetricHartleThorneGB::MetricHartleThorneGB(double mass, double angmom, double eta) {
    mMetricName  = "HartleThorneGB";
    setCoordType(enum_coordinate_spherical);

    mPhysicalUnits = enum_physical_constants_geom;
    mSpeedOfLight = 1.0;
    mGravConstant = 1.0;

    addParam("mass", mass);
    mMass = mass;
    addParam("angmom", angmom);
    mAngmom = angmom;
    addParam("eta", eta);
    mEta = eta;

    setStandardValues();
}

MetricHartleThorneGB::~MetricHartleThorneGB() {
}


// *********************************** public methods ******************************
/*! Calculate the contravariant metric components at position 'pos'.
 *
 *  \param pos : pointer to position.
 */
bool MetricHartleThorneGB::calculateMetric(const double* pos) {
    double M = mMass;
    double a = mAngmom;
    double eta = mEta;

    double r = pos[1];
    double theta = pos[2];
    calcFunc(pos);

    double t1 = M * r;
    double t2 = sigma;      // Sigma(r,theta);
    double t3 = 1 / t2;
    double t6 = a * a;
    double t7 = eta * t6;
    double t8 = htt;        // htt(r,theta);
    double t9 = t2 * t2;
    double t11 = t8 * t2;
    double t14 = M * M;
    double t16 = r * r;
    double t19 = t14 * t6;
    double t20 = sin(theta);
    double t21 = t20 * t20;
    double t22 = t21 * t21;
    double t24 = hphph;       // hpp(r,theta);
    double t29 = 1 / t9;
    double t34 = r * t21;
    double t43 = t21 * t24;
    double t48 = t22 * t24;
    double t57 = -M * a * t34 * t3 + eta * t6 * a * M * t34 * (t11 - 2.0 * t8 * M * r - t43 * t16 * t2 - t43 * t6 * t2 - 2.0 * t48 * M * t6 * r) * t29;
    double t58 = delta;     // Delta(r);
    double t61 = t58 * t58;
    double t64 = hrr;       // hrr(r,theta);
    double t68 = hthth;       // haa(r,theta);
    double t82 = t16 * t16;
    double t92 = t21 * M;
    double t96 = t6 * t6;
    double t97 = t24 * t96;

    g_compts[0][0] = -1.0 + 2.0 * t1 * t3 + t7 * (t8 * t9 - 4.0 * t11 * t1 + 4.0 * t8 * t14 * t16 + 4.0 * t19 * t16 * t22 * t24) * t29;
    g_compts[0][1] = 0.0;
    g_compts[0][2] = 0.0;
    g_compts[0][3] = 2.0 * t57;
    g_compts[1][0] = 0.0;
    g_compts[1][1] = t2 / t58 + t7 * t9 / t61 * t64;
    g_compts[1][2] = 0.0;
    g_compts[1][3] = 0.0;
    g_compts[2][0] = 0.0;
    g_compts[2][1] = 0.0;
    g_compts[2][2] = t2 + t7 * t9 * t68;
    g_compts[2][3] = 0.0;
    g_compts[3][0] = 2.0 * t57;
    g_compts[3][1] = 0.0;
    g_compts[3][2] = 0.0;
    g_compts[3][3] = t21 * t16 + t21 * t6 + 2.0 * t22 * M * t6 * r * t3 + t7 * t22 * (4.0 * t19 * t16 * t8 + t24 * t82 * t9 + 2.0 * t24 * t16 * t9 * t6 + 4.0 * t24 * t16 * r * t2 * t92 * t6 + t97 * t9 + 4.0 * t97 * t2 * t92 * r + 4.0 * t48 * t14 * t96 * t16) * t29;
    return true;
}

/*! Calculate the Christoffel symbols of the second kind at position 'pos'.
 *
 *  \param pos : pointer to position.
 */
bool MetricHartleThorneGB::calculateChristoffels(const double* pos) {
    double M = mMass;
    double a = mAngmom;
    double eta = mEta;

    double r = pos[1];
    double theta = pos[2];
    calcFuncDiff(pos);

    double t1 = sigma;      // Sigma(r,theta);
    double t2 = t1 * t1;
    double t3 = t2 * t2;
    double t4 = 1 / t3;
    double t5 = delta;      // Delta(r);
    double t6 = t5 * t5;
    double t7 = t4 * t6;
    double t8 = a * a;
    double t9 = eta * t8;
    double t10 = hrr;       // hrr(r,theta);
    double t11 = t1 * t10;
    double t14 = 1 / (t5 + t9 * t11);
    double t15 = M * t2;
    double t17 = M * r;
    double t18 = dsigmadr;      // diff(Sigma(r,theta),r);
    double t19 = t18 * t1;
    double t22 = t2 * t1;
    double t23 = dhttdr;        // diff(htt(r,theta),r);
    double t24 = t22 * t23;
    double t26 = t9 * t2;
    double t27 = t23 * M;
    double t28 = t27 * r;
    double t31 = t9 * t1;
    double t32 = htt;           // htt(r,theta);
    double t33 = t32 * t18;
    double t34 = t33 * t17;
    double t37 = t2 * t32;
    double t41 = M * M;
    double t42 = t23 * t41;
    double t43 = r * r;
    double t44 = t42 * t43;
    double t47 = t32 * t41;
    double t48 = t47 * r;
    double t51 = t8 * t8;
    double t52 = eta * t51;
    double t53 = t52 * t1;
    double t54 = t41 * r;
    double t55 = sin(theta);
    double t56 = t55 * t55;
    double t57 = t56 * t56;
    double t58 = hphph;           // hpp(r,theta);
    double t59 = t57 * t58;
    double t60 = t54 * t59;
    double t63 = t41 * t43;
    double t64 = dhppdr;        // diff(hpp(r,theta),r);
    double t65 = t57 * t64;
    double t70 = t47 * t43;
    double t74 = t63 * t59;
    double t77 = 2.0 * t15 - 2.0 * t17 * t19 + t9 * t24 - 4.0 * t26 * t28 + 4.0 * t31 * t34 - 4.0 * t9 * t37 * M + 4.0 * t31 * t44 + 8.0 * t31 * t48 + 8.0 * t53 * t60 + 4.0 * t53 * t63 * t65 - 8.0 * t9 * t18 * t70 - 8.0 * t52 * t18 * t74;
    double t81 = hthth;           // haa(r,theta);
    double t82 = t1 * t81;
    double t85 = 1 / (1.0 + t9 * t82);
    double t86 = t4 * t85;
    double t87 = dsigmadth;     // diff(Sigma(r,theta),theta);
    double t88 = t87 * t1;
    double t89 = t17 * t88;
    double t91 = dhttdth;       // diff(htt(r,theta),theta);
    double t92 = t22 * t91;
    double t94 = t91 * M;
    double t95 = t94 * r;
    double t98 = t32 * t87;
    double t99 = t98 * t17;
    double t102 = t91 * t41;
    double t103 = t102 * t43;
    double t106 = t1 * t41;
    double t107 = t52 * t106;
    double t108 = t56 * t55;
    double t109 = t43 * t108;
    double t110 = cos(theta);
    double t111 = t58 * t110;
    double t112 = t109 * t111;
    double t115 = dhppdth;      // diff(hpp(r,theta),theta);
    double t116 = t57 * t115;
    double t129 = 1 / t2;
    double t130 = t8 * t22;
    double t133 = t43 * t22;
    double t134 = t133 * M;
    double t136 = t56 * t58;
    double t137 = t52 * t136;
    double t138 = t43 * r;
    double t139 = t138 * t2;
    double t140 = M * t18;
    double t141 = t139 * t140;
    double t144 = eta * eta;
    double t145 = t51 * t8;
    double t146 = t144 * t145;
    double t147 = t146 * t136;
    double t148 = t33 * M;
    double t152 = t138 * t22;
    double t153 = t152 * t27;
    double t158 = 16.0 * t147 * t139 * t47;
    double t159 = t8 * t2;
    double t160 = t17 * t18;
    double t163 = t146 * t56;
    double t164 = t58 * t43;
    double t165 = t3 * t23;
    double t170 = t8 * t32;
    double t172 = t133 * eta * t170 * M;
    double t174 = t51 * t22;
    double t175 = t174 * eta;
    double t179 = t51 * t51;
    double t180 = t144 * t179;
    double t181 = t180 * t56;
    double t182 = t58 * t22;
    double t183 = t32 * M;
    double t187 = t51 * t2;
    double t188 = t187 * eta;
    double t191 = eta * t32;
    double t195 = eta * t145;
    double t196 = t195 * t136;
    double t197 = r * t18;
    double t201 = t58 * t3;
    double t204 = t57 * t56;
    double t205 = t58 * t58;
    double t206 = t204 * t205;
    double t207 = t180 * t206;
    double t208 = t43 * t43;
    double t209 = t208 * t1;
    double t210 = t18 * t41;
    double t213 = 8.0 * t207 * t209 * t210;
    double t214 = t180 * t136;
    double t218 = t146 * t59;
    double t219 = t208 * t2;
    double t220 = t219 * t42;
    double t221 = t218 * t220;
    double t223 = t47 * t18;
    double t224 = t209 * t223;
    double t226 = 8.0 * t218 * t224;
    double t228 = t51 * t3;
    double t229 = eta * t23;
    double t231 = -4.0 * t181 * t182 * t183 + 8.0 * t188 * t48 - 4.0 * t174 * t191 * M - 2.0 * t196 * t15 * t197 + t181 * t201 * t23 - t213 + 4.0 * t214 * t37 * t160 - 12.0 * t221 + t226 - 2.0 * t141 + t228 * t229;
    double t233 = t208 * r;
    double t234 = t233 * t145;
    double t235 = t41 * M;
    double t236 = t235 * t57;
    double t237 = t234 * t236;
    double t238 = t144 * t32;
    double t239 = t1 * t64;
    double t240 = t238 * t239;
    double t242 = 8.0 * t237 * t240;
    double t245 = t235 * t32;
    double t248 = 8.0 * t218 * t209 * t245;
    double t249 = t144 * t51;
    double t250 = t249 * t136;
    double t251 = t208 * t43;
    double t252 = t251 * t1;
    double t256 = t152 * eta;
    double t257 = t8 * t23;
    double t258 = t257 * M;
    double t261 = t219 * eta;
    double t265 = t139 * eta;
    double t268 = 8.0 * t265 * t170 * t41;
    double t273 = t146 * t204;
    double t274 = t205 * t233;
    double t278 = 4.0 * t273 * t274 * t2 * t41;
    double t280 = t233 * t2;
    double t282 = t250 * t280 * t47;
    double t284 = t146 * t206;
    double t287 = 4.0 * t284 * t252 * t210;
    double t288 = t204 * t58;
    double t289 = t180 * t288;
    double t290 = t235 * t138;
    double t291 = t290 * t33;
    double t293 = 8.0 * t289 * t291;
    double t294 = t9 * t136;
    double t298 = t249 * t56;
    double t299 = t58 * t208;
    double t302 = t209 * eta;
    double t307 = t33 * t63;
    double t310 = t233 * t22;
    double t314 = t145 * t1;
    double t319 = t56 * t41;
    double t320 = t319 * t51;
    double t321 = t43 * t2;
    double t322 = t321 * t229;
    double t325 = t56 * M;
    double t326 = t325 * t51;
    double t327 = r * t22;
    double t331 = 8.0 * t282 - t287 - t293 - 2.0 * t294 * t280 * t140 + t298 * t299 * t165 - 8.0 * t302 * t8 * t18 * t47 - 8.0 * t53 * t307 - 4.0 * t250 * t310 * t27 - 8.0 * t314 * eta * t18 * t74 - 4.0 * t320 * t322 + 2.0 * t326 * t327 * t229;
    double t334 = t43 * t145;
    double t335 = t334 * t319;
    double t336 = t32 * t32;
    double t337 = t144 * t336;
    double t340 = 4.0 * t335 * t337 * t19;
    double t341 = t52 * t56;
    double t342 = t22 * M;
    double t343 = t164 * t342;
    double t344 = t341 * t343;
    double t346 = t138 * t179;
    double t347 = t235 * t204;
    double t348 = t346 * t347;
    double t350 = 8.0 * t348 * t240;
    double t351 = t208 * t145;
    double t352 = t41 * t57;
    double t353 = t351 * t352;
    double t354 = t2 * t64;
    double t355 = t238 * t354;
    double t356 = t353 * t355;
    double t358 = t195 * t56;
    double t359 = t182 * M;
    double t362 = t43 * t179;
    double t363 = t362 * t352;
    double t364 = t363 * t355;
    double t366 = t346 * t236;
    double t368 = 8.0 * t366 * t240;
    double t369 = t179 * t8;
    double t370 = t144 * t369;
    double t371 = t57 * t57;
    double t372 = t370 * t371;
    double t373 = t205 * t235;
    double t374 = t138 * t18;
    double t377 = 8.0 * t372 * t373 * t374;
    double t380 = t56 * t144;
    double t384 = 4.0 * r * t145 * t41 * t380 * t336 * t2;
    double t385 = t43 * t3;
    double t386 = t9 * t23;
    double t389 = t208 * t22;
    double t391 = t250 * t389 * t183;
    double t393 = t251 * t2;
    double t397 = t205 * t1;
    double t398 = t235 * t43;
    double t401 = 8.0 * t372 * t397 * t398;
    double t402 = t180 * t59;
    double t403 = t1 * t235;
    double t405 = t403 * t138 * t23;
    double t407 = 8.0 * t402 * t405;
    double t408 = t19 * t70;
    double t411 = t370 * t204;
    double t412 = t205 * t2;
    double t415 = 4.0 * t411 * t412 * t54;
    double t416 = t37 * t54;
    double t421 = t170 * t140;
    double t425 = 8.0 * t289 * t405;
    double t427 = 8.0 * t402 * t291;
    double t428 = -4.0 * t391 + 4.0 * t250 * t393 * t42 + t401 + t407 - 8.0 * t214 * t408 + t415 + 8.0 * t214 * t416 - 16.0 * t147 * t224 + 4.0 * t265 * t421 - t425 - t427;
    double t431 = 8.0 * t402 * t408;
    double t433 = t2 * t23 * t63;
    double t437 = 8.0 * t402 * t416;
    double t438 = t209 * t52;
    double t442 = t1 * t32;
    double t443 = t398 * t442;
    double t445 = 8.0 * t289 * t443;
    double t446 = t402 * t433;
    double t451 = 8.0 * t402 * t443;
    double t452 = t370 * t206;
    double t455 = 4.0 * t452 * t19 * t63;
    double t456 = t233 * t1;
    double t460 = 8.0 * t218 * t456 * t235 * t23;
    double t461 = t233 * t235;
    double t464 = 8.0 * t218 * t461 * t33;
    double t465 = t431 + 4.0 * t214 * t433 - t437 - 8.0 * t438 * t210 * t59 + t445 - 12.0 * t446 + 4.0 * t289 * t433 + t451 - t455 + t460 - t464;
    double t466 = t180 * t371;
    double t467 = t205 * t208;
    double t470 = 8.0 * t466 * t467 * t403;
    double t474 = 8.0 * t466 * t274 * t235 * t18;
    double t475 = t24 * t17;
    double t481 = t147 * t133 * t183;
    double t488 = t145 * t2;
    double t489 = t488 * eta;
    double t492 = t9 * t56;
    double t494 = t492 * t299 * t342;
    double t496 = t336 * t1;
    double t499 = 8.0 * t163 * t398 * t496;
    double t500 = t336 * t18;
    double t503 = 8.0 * t163 * t290 * t500;
    double t504 = -t470 - t474 - 4.0 * t214 * t475 + 4.0 * t402 * t475 - 8.0 * t481 + 4.0 * t218 * t153 + 4.0 * t250 * t280 * t148 + 8.0 * t489 * t60 + 2.0 * t494 + t499 - t503;
    double t509 = M * t138;
    double t510 = t509 * t1;
    double t514 = t52 * t57;
    double t515 = t58 * t138;
    double t516 = t1 * M;
    double t520 = t146 * t32;
    double t521 = t41 * t208;
    double t525 = t180 * t32;
    double t529 = M * t233;
    double t538 = t249 * t32;
    double t546 = t51 * t57;
    double t557 = t52 * t32;
    double t558 = t1 * t56;
    double t562 = t146 * t442;
    double t573 = 2.0 * t510 - t159 - 4.0 * t341 * t70 - 4.0 * t514 * t515 * t516 - t321 + 8.0 * t520 * t521 * t136 + 4.0 * t525 * t63 * t136 + 2.0 * t529 * t1 * t294 + 4.0 * t9 * t47 * t208 + 4.0 * t52 * t70 + 4.0 * t538 * t41 * t251 * t136 + 4.0 * t510 * t137 + 4.0 * t521 * eta * t546 * t58 + 4.0 * t63 * eta * t145 * t57 * t58 + t9 * t37 * t43 + 2.0 * t557 * t558 * t17 - 8.0 * t562 * t509 * t136 - 4.0 * t249 * t442 * t529 * t136 - 8.0 * t520 * t521 * t59;
    double t574 = t180 * t442;
    double t583 = t9 * t32;
    double t586 = t516 * r;
    double t589 = t136 * t2;
    double t591 = t195 * t57;
    double t592 = t58 * t1;
    double t596 = t195 * t204;
    double t598 = t58 * t41 * t43;
    double t604 = t59 * t17;
    double t609 = t164 * t2;
    double t616 = t2 * t56;
    double t617 = t616 * t164;
    double t625 = t52 * t37;
    double t629 = -4.0 * t574 * t17 * t136 - 8.0 * t525 * t74 + 4.0 * t525 * t288 * t63 - 4.0 * t583 * t510 - 4.0 * t557 * t586 - t195 * t589 - 4.0 * t591 * t592 * t17 - 4.0 * t596 * t598 + 4.0 * t562 * t59 * t509 + 4.0 * t574 * t604 - t492 * t299 * t2 - 2.0 * t341 * t609 - 2.0 * t325 * t8 * r * t1 + 2.0 * t520 * t617 + t538 * t616 * t299 + t525 * t589 + 2.0 * t586 * t196 + t625 + 2.0 * t17 * t1 * t8;
    double t631 = 1 / (t573 + t629);
    double t633 = t129 * (4.0 * t261 * t257 * t41 + 2.0 * t130 * M - 8.0 * t250 * t252 * t223 - t242 - 4.0 * t172 - 4.0 * t175 * t28 + 8.0 * t147 * t139 * t148 + 2.0 * t163 * t164 * t165 - 4.0 * t256 * t258 - t248 - t278 + t331 + t158 - 2.0 * t159 * t160 + 2.0 * t358 * t359 - 8.0 * t147 * t153 + t350 + 4.0 * t356 + t340 + 4.0 * t344 + 8.0 * t147 * t220 - t384 + t385 * t386 - t368 + t428 - 4.0 * t137 * t141 + t465 + t504 + t268 + 4.0 * t188 * t44 + 4.0 * t188 * t34 + 2.0 * t134 + 4.0 * t364 - t377 + t231) * t631 / 2.0;
    double t634 = a * M;
    double t635 = t515 * t24;
    double t637 = M * t43;
    double t639 = t56 * t64;
    double t640 = t52 * t639;
    double t643 = t456 * t42;
    double t646 = t249 * t47;
    double t650 = M * t208;
    double t652 = t9 * t639;
    double t658 = t467 * t106;
    double t665 = t274 * t210;
    double t668 = t180 * t204;
    double t670 = t41 * t138;
    double t671 = t205 * t18 * t670;
    double t674 = t515 * t15;
    double t677 = t146 * t47;
    double t678 = t138 * t1;
    double t682 = t136 * t22;
    double t684 = t1 * t23;
    double t685 = t670 * t684;
    double t688 = t298 * t635 - 2.0 * t637 * t2 * t640 + 4.0 * t250 * t643 - 4.0 * t646 * t456 * t639 - 2.0 * t650 * t2 * t652 - 4.0 * t250 * t209 * t47 - 4.0 * t273 * t658 - 4.0 * t250 * t233 * t18 * t47 - 4.0 * t273 * t665 - 4.0 * t668 * t671 - 4.0 * t492 * t674 - 4.0 * t677 * t678 * t639 + t52 * t682 - 4.0 * t218 * t685;
    double t689 = t321 * t27;
    double t692 = t43 * t18;
    double t697 = t146 * t37;
    double t698 = r * t56;
    double t699 = t58 * t18;
    double t700 = t698 * t699;
    double t702 = t43 * t57;
    double t703 = t64 * M;
    double t704 = t702 * t703;
    double t707 = t249 * t336;
    double t708 = t1 * t43;
    double t709 = t708 * t140;
    double t715 = t327 * t639;
    double t720 = t146 * t183;
    double t721 = t321 * t639;
    double t724 = t249 * t183;
    double t725 = t138 * t56;
    double t726 = t58 * t2;
    double t732 = t219 * t27;
    double t737 = t33 * t670;
    double t740 = -4.0 * t147 * t689 + 4.0 * t562 * t692 * t59 * M + t697 * t700 - 2.0 * t697 * t704 + 4.0 * t707 * t709 - 4.0 * t677 * t374 * t59 - t520 * t715 + 4.0 * t677 * t708 * t59 + 4.0 * t720 * t721 + 8.0 * t724 * t725 * t726 + 2.0 * t218 * t689 - 4.0 * t250 * t732 - 4.0 * t697 * t604 - 4.0 * t147 * t737;
    double t742 = t183 * r;
    double t745 = t32 * t22;
    double t748 = t219 * t639;
    double t751 = t63 * t1;
    double t760 = t9 * t138;
    double t761 = t22 * t56;
    double t762 = t761 * t64;
    double t769 = t761 * t164;
    double t772 = t52 * r;
    double t774 = t152 * t639;
    double t776 = t249 * t37;
    double t779 = t52 * t43;
    double t780 = t2 * t57;
    double t784 = 4.0 * t26 * t742 - 2.0 * t9 * t745 + 4.0 * t724 * t748 + 4.0 * t707 * t751 - 4.0 * t707 * t670 * t18 + 4.0 * t677 * t678 * t65 + t760 * t762 + 4.0 * t188 * t604 + 3.0 * t492 * t164 * t22 - 3.0 * t538 * t769 + t772 * t762 - t538 * t774 + t776 * t725 * t699 + 2.0 * t779 * t780 * t703;
    double t785 = t197 * t2;
    double t788 = t397 * t63;
    double t791 = t15 * r;
    double t794 = t616 * t699;
    double t796 = t9 * r;
    double t802 = t182 * r * t23;
    double t806 = t442 * t63;
    double t809 = t1 * eta;
    double t810 = t51 * t43;
    double t813 = t58 * M;
    double t814 = t18 * t57 * t813;
    double t817 = t692 * t183;
    double t820 = t336 * t22;
    double t822 = -t707 * t785 - t520 * t682 + 4.0 * t668 * t788 - 4.0 * t707 * t791 - t772 * t794 + 2.0 * t796 * t37 * t18 - t760 * t794 + t163 * t802 + 4.0 * t147 * t685 + 4.0 * t147 * t806 - 4.0 * t809 * t810 * t814 - 4.0 * t31 * t817 + t22 + t249 * t820 - t785;
    double t827 = t634 * t129 * (t688 + t740 + t784 + t822) * t631;
    double t828 = t41 * t108;
    double t830 = t2 * t110;
    double t831 = t830 * t58;
    double t832 = t238 * t831;
    double t835 = t57 * t108;
    double t836 = t835 * t205;
    double t837 = t370 * t836;
    double t838 = t138 * t110;
    double t846 = t1 * t91;
    double t847 = t290 * t846;
    double t850 = t2 * t91;
    double t851 = t63 * t850;
    double t854 = t371 * t55;
    double t855 = t370 * t854;
    double t856 = t41 * t41;
    double t858 = t208 * t110;
    double t864 = t98 * t63;
    double t871 = t9 * t91;
    double t873 = M * t87;
    double t874 = t139 * t873;
    double t876 = t290 * t98;
    double t879 = t98 * t41;
    double t888 = t98 * M;
    double t889 = t139 * t888;
    double t892 = t152 * t94;
    double t896 = t43 * t8 * t41;
    double t904 = 24.0 * t351 * t828 * t832 - 16.0 * t837 * t403 * t838 - 8.0 * t372 * t373 * t138 * t87 - 8.0 * t289 * t847 + 4.0 * t289 * t851 - 32.0 * t855 * t205 * t856 * t858 + 8.0 * t402 * t847 - 8.0 * t53 * t864 - 8.0 * t314 * eta * t87 * t74 + t385 * t871 - 2.0 * t874 - 8.0 * t402 * t876 - 8.0 * t250 * t252 * t879 + t181 * t201 * t91 + 4.0 * t250 * t393 * t102 + 8.0 * t147 * t889 - 8.0 * t147 * t892 - 8.0 * t896 * t55 * t110 * t2 + 4.0 * t250 * t280 * t888;
    double t905 = t219 * t102;
    double t912 = t17 * t87;
    double t921 = t87 * t41;
    double t929 = t209 * t879;
    double t937 = t1 * t115;
    double t938 = t238 * t937;
    double t941 = t2 * t115;
    double t942 = t238 * t941;
    double t945 = t41 * t55;
    double t954 = t351 * t856;
    double t956 = t336 * t110;
    double t963 = -12.0 * t218 * t905 + 8.0 * t218 * t456 * t235 * t91 - 2.0 * t159 * t912 + 4.0 * t188 * t99 + 4.0 * t188 * t103 - 8.0 * t289 * t876 - 8.0 * t207 * t209 * t921 + 4.0 * t218 * t892 + 8.0 * t147 * t905 - 16.0 * t147 * t929 - 4.0 * t284 * t252 * t921 - 4.0 * t137 * t874 + 8.0 * t348 * t938 + 4.0 * t353 * t942 - 8.0 * t334 * t945 * t337 * t830 - 8.0 * t366 * t938 + 4.0 * t363 * t942 - 32.0 * t954 * t55 * t144 * t956 - 8.0 * t438 * t921 * t59;
    double t972 = eta * t91;
    double t976 = t57 * t55;
    double t979 = t208 * t32;
    double t984 = r * t87;
    double t988 = t92 * t17;
    double t1009 = t3 * t91;
    double t1013 = t170 * t873;
    double t1020 = t8 * t91;
    double t1024 = t43 * t32;
    double t1026 = t106 * t1024 * t87;
    double t1029 = -8.0 * t237 * t938 - 4.0 * t250 * t310 * t94 + 8.0 * t218 * t929 + 2.0 * t326 * t327 * t972 - 64.0 * t180 * t976 * t856 * t979 * t111 + t228 * t972 - 2.0 * t196 * t15 * t984 - 4.0 * t214 * t988 - 4.0 * t320 * t321 * t972 - 4.0 * t452 * t88 * t63 + 4.0 * t402 * t988 - 12.0 * t402 * t851 - 8.0 * t219 * t52 * t828 * t111 - 8.0 * t302 * t8 * t87 * t47 + 2.0 * t163 * t164 * t1009 + 4.0 * t265 * t1013 - 8.0 * t466 * t274 * t235 * t87 + 4.0 * t261 * t1020 * t41 + 8.0 * t402 * t1026;
    double t1030 = t976 * t235;
    double t1031 = t145 * t138;
    double t1033 = t809 * t111;
    double t1036 = t336 * t87;
    double t1037 = t290 * t1036;
    double t1045 = t235 * t108;
    double t1047 = t110 * t1;
    double t1048 = t1047 * t58;
    double t1049 = t238 * t1048;
    double t1052 = t810 * t945;
    double t1053 = eta * t110;
    double t1054 = t1053 * t37;
    double t1057 = t138 * t51;
    double t1058 = t235 * t55;
    double t1060 = t1053 * t442;
    double t1063 = t180 * t836;
    double t1071 = t346 * t1030;
    double t1077 = t1020 * M;
    double t1094 = t37 * t912;
    double t1103 = -32.0 * t1030 * t1031 * t1033 - 8.0 * t163 * t1037 - 2.0 * t294 * t280 * t873 + t298 * t299 * t1009 - 48.0 * t346 * t1045 * t1049 + 16.0 * t1052 * t1054 - 32.0 * t1057 * t1058 * t1060 - 16.0 * t1063 * t456 * t235 * t110 + 4.0 * t335 * t337 * t88 + 64.0 * t1071 * t1049 + 24.0 * t362 * t828 * t832 - 4.0 * t256 * t1077 - 8.0 * t214 * t1026 - 8.0 * t488 * eta * t41 * t112 + 4.0 * t214 * t851 - 48.0 * t234 * t1045 * t1049 - 8.0 * t218 * t461 * t98 + 4.0 * t214 * t1094 + 32.0 * t1031 * t1058 * t337 * t1047 - 4.0 * t175 * t95;
    double t1108 = t129 * (t904 + t963 + t1029 + t1103) * t631 / 2.0;
    double t1110 = t55 * t87;
    double t1111 = t1110 * t2;
    double t1112 = r * t110;
    double t1116 = t63 * t1047;
    double t1122 = t110 * t22;
    double t1126 = t110 * t57;
    double t1127 = t813 * r;
    double t1131 = t976 * t115;
    double t1132 = t1131 * t17;
    double t1139 = t976 * t87;
    double t1143 = t2 * t108;
    double t1144 = t58 * t87;
    double t1147 = t22 * t108;
    double t1148 = t115 * t43;
    double t1151 = t108 * t58;
    double t1152 = t43 * t87;
    double t1158 = t52 * t110;
    double t1161 = t249 * t1151;
    double t1162 = t139 * t94;
    double t1165 = -t1111 + 16.0 * t720 * t1112 * t589 + 24.0 * t707 * t1116 + 16.0 * t724 * t838 * t589 - 4.0 * t9 * t1122 * t32 - 12.0 * t697 * t1126 * t1127 - 2.0 * t697 * t1132 + 4.0 * t249 * t496 * t1110 * t17 + 4.0 * t562 * t1139 * t1127 + t520 * t1143 * t1144 - t538 * t1147 * t1148 + t776 * t1151 * t1152 - 4.0 * t677 * t109 * t937 + 4.0 * t1158 * t682 - 4.0 * t1161 * t1162;
    double t1166 = t1147 * t115;
    double t1168 = t9 * t110;
    double t1172 = t138 * t57;
    double t1176 = t146 * t1151;
    double t1177 = t846 * t63;
    double t1180 = t108 * t115;
    double t1181 = t52 * t1180;
    double t1184 = t110 * t32;
    double t1189 = t9 * t1180;
    double t1192 = t2 * eta;
    double t1195 = t51 * t110 * t136;
    double t1198 = t209 * t102;
    double t1201 = t43 * t110;
    double t1202 = t558 * t58;
    double t1213 = t8 * t110;
    double t1214 = t1213 * t136;
    double t1217 = t850 * t17;
    double t1220 = t1 * t57;
    double t1225 = t146 * t108;
    double t1226 = t182 * t91;
    double t1228 = -t520 * t1166 + 4.0 * t1168 * t769 - 32.0 * t146 * t245 * t1172 * t111 + 4.0 * t1176 * t1177 - 2.0 * t791 * t1181 - 16.0 * t751 * t9 * t1184 - 2.0 * t509 * t2 * t1189 - 8.0 * t17 * t1192 * t1195 + 4.0 * t1161 * t1198 - 16.0 * t677 * t1201 * t1202 - 16.0 * t646 * t858 * t1202 - 16.0 * t707 * t290 * t110 - 8.0 * t509 * t1192 * t1214 - 4.0 * t1176 * t1217 + 40.0 * t677 * t1201 * t1220 * t58 + t1225 * t1226;
    double t1237 = t976 * t58;
    double t1238 = t146 * t1237;
    double t1242 = t205 * t87;
    double t1243 = t1242 * t63;
    double t1246 = t37 * t17;
    double t1251 = t373 * t838;
    double t1258 = t51 * t976;
    double t1259 = t809 * t1258;
    double t1260 = t1144 * t17;
    double t1263 = t138 * t108;
    double t1272 = t43 * t976;
    double t1276 = r * t108;
    double t1283 = t17 * t830;
    double t1285 = -16.0 * t107 * t702 * t111 - 4.0 * t1161 * t208 * t87 * t47 + 2.0 * t1238 * t1217 - 4.0 * t180 * t835 * t1243 + 16.0 * t1168 * t1246 - 4.0 * t1176 * t864 - 16.0 * t466 * t1251 - 4.0 * t809 * t8 * t55 * t99 - 4.0 * t1259 * t1260 + 4.0 * t724 * t1263 * t941 - 4.0 * t1238 * t1177 - 4.0 * t707 * t63 * t1110 - 4.0 * t677 * t1272 * t1144 + 4.0 * t720 * t1276 * t941 + 4.0 * t677 * t1272 * t937 - 4.0 * t1283;
    double t1287 = t467 * t921;
    double t1290 = t208 * t108;
    double t1298 = t110 * t56;
    double t1305 = t9 * t55;
    double t1309 = t9 * t108;
    double t1310 = t22 * t115;
    double t1311 = t1310 * t43;
    double t1315 = t52 * t108;
    double t1318 = t52 * t976;
    double t1319 = t941 * t17;
    double t1331 = t164 * t92;
    double t1333 = -4.0 * t146 * t835 * t1287 - 4.0 * t646 * t1290 * t937 - 4.0 * t520 * t1122 * t136 - 4.0 * t249 * t745 * t1298 * t164 + 2.0 * t249 * t820 * t110 + 2.0 * t1305 * t37 * t87 + t1309 * t1311 - t1309 * t726 * t1152 - t1315 * t726 * t87 + 2.0 * t1318 * t1319 + 2.0 * t1122 + 12.0 * t52 * t830 * t604 - 12.0 * t707 * t1283 - t707 * t1111 + t52 * t1166 + t249 * t108 * t1331;
    double t1338 = t631 / t55;
    double t1340 = t17 * a * t129 * (t1165 + t1228 + t1285 + t1333) * t1338;
    double t1342 = t197 * t1;
    double t1351 = r * t2;
    double t1356 = t9 * t43;
    double t1360 = t616 * t64;
    double t1362 = t558 * t699;
    double t1373 = -t2 + t1342 + t9 * t37 - 4.0 * t583 * t586 - 3.0 * t492 * t609 - t52 * t589 - 4.0 * t53 * t604 + t9 * t1351 * t23 - t796 * t442 * t18 - 2.0 * t1356 * t684 * M - t760 * t1360 + t760 * t1362 - t772 * t1360 + t772 * t1362 - 2.0 * t779 * t1220 * t703 + 4.0 * t1356 * t148 + 4.0 * t779 * t814;
    double t1376 = t7 * t14 * t634 * t56 * t1373;
    double t1381 = t830 * t32;
    double t1408 = t592 * t87;
    double t1417 = -2.0 * t830 + t1110 * t1 + 2.0 * t9 * t1381 - 4.0 * t1168 * t442 * t17 - 4.0 * t1168 * t617 - 4.0 * t1158 * t589 - 12.0 * t52 * t1047 * t604 + t9 * t55 * t2 * t91 - t1305 * t442 * t87 - 2.0 * t1305 * t846 * t17 - t1309 * t941 * t43 + t1309 * t592 * t1152 - t52 * t1143 * t115 + t1315 * t1408 - 2.0 * t1318 * t937 * t17 + 4.0 * t1305 * t99 + 4.0 * t1318 * t1260;
    double t1420 = t86 * M * a * r * t55 * t1417;
    double t1421 = 1 / t1;
    double t1425 = deltadr;         // diff(Delta(r),r);
    double t1436 = dhrrdr;          // diff(hrr(r,theta),r);
    double t1449 = dhrrdth;         // diff(hrr(r,theta),theta);
    double t1452 = t87 * t5 + 2.0 * t9 * t11 * t87 + t9 * t2 * t1449;
    double t1457 = t1421 * t14;
    double t1459 = t1457 * t1452 / 2.0;
    double t1463 = dhaadr;        // diff(haa(r,theta),r);
    double t1466 = t18 + 2.0 * t9 * t82 * t18 + t9 * t2 * t1463;
    double t1469 = t1466 * t85 * t1421 / 2.0;
    double t1470 = t634 * t56;
    double t1472 = t146 * t352;
    double t1473 = t233 * t32;
    double t1492 = t370 * t57;
    double t1513 = t1492 * t412 * t197 + 4.0 * t214 * t708 * t148 - t214 * t1351 * t33 + t181 * t802 - 4.0 * t411 * t412 * t17 + 4.0 * t466 * t665 + t181 * t182 * t32 - 4.0 * t214 * t1246 + 4.0 * t218 * t732 + t174 * t191 - t130;
    double t1515 = t180 * t352;
    double t1522 = t138 * t32;
    double t1530 = t8 * t56;
    double t1559 = r * t32;
    double t1567 = t43 * t23;
    double t1575 = -t370 * t57 * t205 * t22 - t133 * t583 + 4.0 * t265 * t546 * t813 + t152 * t386 - t265 * t170 * t18 - t188 * t1559 * t18 + t174 * eta * r * t23 - 4.0 * t489 * t604 - 2.0 * t188 * t1567 * M + 4.0 * t438 * t814 + 4.0 * t302 * t421;
    double t1580 = t145 * t22;
    double t1585 = t180 * t183;
    double t1590 = t204 * t64;
    double t1603 = t761 * t64 * t233;
    double t1609 = t180 * t37;
    double t1612 = t251 * t56;
    double t1613 = t1612 * t354;
    double t1626 = t51 * t56;
    double t1641 = -4.0 * t697 * t65 * t650 - 4.0 * t1609 * t704 + 2.0 * t724 * t1613 - r * t144 * t179 * t745 * t639 - 2.0 * t261 * t258 + 2.0 * t280 * eta * t1530 * t699 + 4.0 * t265 * t1626 * t699 + 4.0 * t466 * t658 + 8.0 * t147 * t209 * t148 + 8.0 * t207 * t209 * t140 - 2.0 * t147 * t139 * t33;
    double t1651 = t180 * t57;
    double t1652 = t205 * t43;
    double t1661 = t2 * t18;
    double t1668 = t249 * t57;
    double t1685 = t146 * t57;
    double t1695 = t205 * t251;
    double t1706 = t298 * t58 * t233 * t24 - t250 * t280 * t33 - 2.0 * t250 * t393 * t27 + 3.0 * t1685 * t274 * t1661 - 3.0 * t298 * t299 * t745 + 4.0 * t250 * t280 * t183 + t1668 * t1695 * t22 + t1685 * t467 * t22 + 4.0 * t163 * t670 * t500 + t159 * t197 + 2.0 * t489 * t700;
    double t1712 = t1470 * t129 * (-8.0 * t218 * t139 * t183 + 4.0 * t250 * t252 * t148 - t1651 * t1652 * t22 + t133 - 2.0 * t163 * t164 * t745 + t1668 * t205 * t208 * t138 * t1661 + 4.0 * t1472 * t1473 * t699 + 4.0 * t1515 * t1522 * t239 + 4.0 * t1515 * t1522 * t699 - 4.0 * t574 * t1590 * t670 + 3.0 * t1651 * t205 * t138 * t1661 - 4.0 * t1515 * t1024 * t592 + 4.0 * t1472 * t979 * t592 + 4.0 * t1472 * t1473 * t239 + 4.0 * t284 * t252 * t140 + t1706 + 4.0 * t273 * t274 * t15 - 4.0 * t163 * t63 * t496 - 4.0 * t147 * t732 + 4.0 * t402 * t689 - 4.0 * t188 * t742 + 4.0 * t720 * t748 + 4.0 * t372 * t671 - 4.0 * t372 * t788 + 2.0 * t326 * t322 + 2.0 * t389 * eta * t1530 * t58 + 4.0 * t314 * eta * t43 * t814 - 2.0 * t1580 * eta * t56 * t58 + 4.0 * t452 * t709 - 4.0 * t402 * t685 + 4.0 * t289 * t685 - 2.0 * t214 * t689 + 4.0 * t289 * t737 - 4.0 * t289 * t806 - t538 * t1603 + t1641 + 2.0 * t163 * t635 + 4.0 * t53 * t817 - 4.0 * t218 * t643 + t1575 - 2.0 * t520 * t774 + t1513 + 2.0 * t1585 * t721 + t139 * t18) * t631;
    double t1718 = t146 * t745;
    double t1723 = t3 * t56;
    double t1724 = t1723 * t64;
    double t1727 = t509 * t22;
    double t1731 = -4.0 * t596 * t354 * t63 - 2.0 * r * t3 + 4.0 * t1718 * t65 * t509 + 4.0 * t134 + t158 + t525 * t1724 - 8.0 * t172 + 4.0 * t1727 * t640 + t213 + 4.0 * t221 + t226;
    double t1736 = -12.0 * t514 * t343 + t242 + t248 + t268 + t278 + 16.0 * t282 + t287 + t293 - t340 + 8.0 * t344 - t350;
    double t1748 = t180 * t745;
    double t1752 = -12.0 * t356 - 12.0 * t364 + t368 + t377 + t384 - 16.0 * t391 + 4.0 * t525 * t22 * t57 * t813 + t538 * t1723 * t64 * t208 - t401 + 4.0 * t1748 * t65 * t17 - t407;
    double t1778 = t208 * t51;
    double t1779 = t1778 * t41;
    double t1780 = t57 * eta;
    double t1781 = t1780 * t354;
    double t1784 = t334 * t41;
    double t1791 = -8.0 * t720 * t774 - 32.0 * t677 * t1172 * t726 - 4.0 * t1609 * t59 * t160 + 4.0 * t538 * t1723 * t515 + 2.0 * t520 * t1723 * t64 * t43 + 2.0 * t557 * t761 * M + 8.0 * t1609 * t288 * t54 - 8.0 * t574 * t204 * t18 * t598 + 4.0 * t1779 * t1781 + 4.0 * t1784 * t1781 - 4.0 * t1585 * t715 + 8.0 * t677 * t748;
    double t1797 = t325 * t8;
    double t1800 = t180 * t47;
    double t1805 = t699 * t63;
    double t1808 = t17 * t22;
    double t1814 = t529 * t22;
    double t1821 = -t415 + 8.0 * t809 * t1626 * t307 + 2.0 * t1797 * t785 + 4.0 * t1800 * t721 + t425 + t427 + 8.0 * t809 * t145 * t204 * t1805 + 2.0 * t1808 * t195 * t639 + 4.0 * t646 * t1613 + 2.0 * t1814 * t652 + 16.0 * t1057 * t41 * t1780 * t726;
    double t1840 = t3 * t64;
    double t1843 = t22 * t64;
    double t1853 = 4.0 * t1609 * t1590 * t63 - 4.0 * t724 * t1603 + 2.0 * t9 * t32 * t3 * r + t431 + 12.0 * t1718 * t59 * t637 - 4.0 * t697 * t59 * t374 * M - 8.0 * t341 * t416 - t492 * t1840 * t208 - 4.0 * t591 * t1843 * t17 - 8.0 * t596 * t726 * t54 - 4.0 * t514 * t1843 * t509;
    double t1856 = t515 * t140;
    double t1879 = t699 * t17;
    double t1882 = 4.0 * t52 * t780 * t1856 - t195 * t1724 - 4.0 * t492 * t201 * t138 - 2.0 * t341 * t1840 * t43 - 4.0 * t591 * t359 - 4.0 * t341 * t201 * r + 4.0 * t520 * t1723 * t58 * r - 2.0 * t625 * t325 * t197 + 4.0 * t195 * t780 * t1879 - t437 - t445;
    double t1888 = 4.0 * t446 - t451 + t455 - t460 + t464 + t470 + t474 - 2.0 * t325 * t130 - 16.0 * t481 + 8.0 * t494 - t499 + t503;
    double t1894 = t129 * (t1731 + t1736 + t1752 + t1791 + t1821 + t1853 + t1882 + t1888) * t631 / 2.0;
    double t1901 = dhaadth;         // diff(haa(r,theta),theta);
    double t1917 = t43 * t336;
    double t1942 = t146 * t828;
    double t1946 = M * t110;
    double t1947 = t1946 * t32;
    double t1948 = t139 * t1947;
    double t1955 = t180 * t1237;
    double t1968 = -t321 * t87 - t159 * t87 - t181 * t1226 + t525 * t761 * t115 + 4.0 * t147 * t1162 + 8.0 * t146 * t945 * t1917 * t1047 + 4.0 * t697 * t116 * t509 - 2.0 * t1585 * t698 * t941 - 4.0 * t720 * t725 * t941 - 4.0 * t163 * t63 * t1036 - 4.0 * t289 * t1177 - 4.0 * t289 * t864 - 4.0 * t372 * t1243 - 4.0 * t678 * eta * t1013 - 32.0 * t1942 * t979 * t1048 + 16.0 * t1176 * t1948 - 32.0 * t180 * t828 * t1024 * t1048 + 24.0 * t1955 * t63 * t1047 * t32 + 4.0 * t574 * t204 * t115 * t63 - 2.0 * t724 * t233 * t56 * t941;
    double t1971 = t139 * t52;
    double t1977 = t321 * eta;
    double t1980 = t976 * t41;
    double t2003 = t55 * t58;
    double t2005 = t830 * t742;
    double t2008 = t976 * t205;
    double t2015 = t55 * M;
    double t2017 = t1351 * t110;
    double t2022 = t314 * t1780;
    double t2034 = -4.0 * t466 * t1287 + 2.0 * t180 * t55 * t182 * t1184 - 4.0 * t180 * t2003 * t2005 - 4.0 * t370 * t2008 * t1283 + 2.0 * t520 * t761 * t1148 - 4.0 * t2015 * t8 * t2017 - 4.0 * t53 * t99 - 4.0 * t2022 * t1260 - 2.0 * t489 * t136 * t87 - 8.0 * t488 * t1053 * t1151 * t17 + t187 * t191 * t87;
    double t2046 = t51 * r;
    double t2069 = t146 * t2003;
    double t2089 = t2 * t87;
    double t2102 = -4.0 * t1977 * t1626 * t1144 + 2.0 * t188 * t95 - 16.0 * t855 * t1251 - 4.0 * t1472 * t979 * t937 + t250 * t219 * t98 + 2.0 * t250 * t280 * t94 - t1668 * t1695 * t2089 - 3.0 * t1685 * t467 * t2089 - t298 * t299 * t92 - 2.0 * t163 * t1331 + 2.0 * t147 * t321 * t98;
    double t2117 = t115 * t208;
    double t2133 = t249 * t2003;
    double t2162 = 2.0 * t2133 * t389 * t1184 - 4.0 * t214 * t88 * t742 - 4.0 * t452 * t89 + 16.0 * t180 * t1151 * t2005 + t181 * t726 * t98 + 2.0 * t214 * t1217 - 2.0 * t326 * t1351 * t972 - 4.0 * t2133 * t280 * t1947 - 4.0 * t146 * t2008 * t280 * t1946 + 2.0 * t265 * t1077 - t174 * t972;
    double t2168 = t1470 * r * t129 * (-8.0 * t1971 * t110 * t108 * t813 + t1977 * t170 * t87 - 16.0 * t1980 * t334 * t1033 - 4.0 * t250 * t456 * t888 - 3.0 * t1651 * t1652 * t2089 - 8.0 * t180 * t2008 * t139 * t1946 - 32.0 * t180 * t1030 * t1522 * t111 - 16.0 * t146 * t55 * t290 * t956 - 16.0 * t837 * t1116 - 16.0 * t1063 * t209 * t41 * t110 + t538 * t761 * t2117 - 8.0 * t207 * t678 * t873 - 8.0 * t147 * t678 * t888 + 4.0 * t2069 * t133 * t1184 - 4.0 * t284 * t456 * t873 - t1492 * t412 * t87 + 4.0 * t1609 * t116 * t17 + 4.0 * t2015 * t2046 * t1054 - 4.0 * t1515 * t1024 * t1144 - 4.0 * t1515 * t1024 * t937 - 4.0 * t1472 * t979 * t1144 - 2.0 * t261 * t1530 * t1144 + 4.0 * t218 * t1198 + t2034 + t1968 - 8.0 * t2069 * t1948 - 16.0 * t1052 * t1060 - 4.0 * t678 * t52 * t57 * t87 * t813 + t2102 - t133 * t871 + t2162 - 4.0 * t402 * t1217 - 4.0 * t218 * t1162 + 4.0 * t402 * t1177) * t631;
    double t2169 = t3 * t108;
    double t2180 = t108 * M;
    double t2181 = t2180 * t8;
    double t2190 = t835 * t41;
    double t2199 = t385 * eta;
    double t2203 = t976 * eta;
    double t2204 = t2203 * t941;
    double t2212 = 2.0 * t520 * t2169 * t1148 + 4.0 * t181 * t201 * t1184 - 2.0 * t385 * t110 + t538 * t2169 * t2117 + 2.0 * t2181 * t1351 * t87 - 4.0 * t1580 * eta * t1132 + 8.0 * t1225 * t1037 - 4.0 * t2190 * t145 * t321 * eta * t115 - 8.0 * t175 * t1184 * t17 + 2.0 * t2199 * t1213 * t32 + 4.0 * t1779 * t2204 + 2.0 * t1808 * t195 * t1180 + 4.0 * t1727 * t1181;
    double t2226 = t208 * t3;
    double t2227 = t2226 * eta;
    double t2238 = t138 * t369;
    double t2241 = t854 * t144 * t1242;
    double t2258 = 32.0 * t954 * t380 * t956 + 32.0 * t208 * t369 * t856 * t371 * t56 * t144 * t205 * t110 + 4.0 * t509 * t1122 - 4.0 * t2227 * t1214 - 8.0 * t2199 * t1195 - 4.0 * t256 * t1258 * t115 * M + 2.0 * t1814 * t1189 + 8.0 * t2238 * t235 * t2241 - 2.0 * t8 * t3 * t110 - 8.0 * t256 * t1213 * t183 + 8.0 * t583 * t521 * t830 + 8.0 * t557 * t63 * t830 + 4.0 * t1784 * t2204;
    double t2260 = t233 * t179;
    double t2274 = t208 * t179;
    double t2284 = t1053 * t726;
    double t2299 = t342 * t1112 * t32;
    double t2304 = 8.0 * t2260 * t235 * t2241 - 2.0 * t2199 * t51 * t108 * t115 + 8.0 * t896 * t1298 * t2 - t2227 * t8 * t108 * t115 + 64.0 * t2274 * t856 * t204 * t238 * t111 - 4.0 * t1942 * t1917 * t88 + 32.0 * t1778 * t352 * t2284 + 4.0 * t1238 * t905 - 8.0 * t720 * t1263 * t1310 + 4.0 * t1800 * t109 * t941 - 4.0 * t1585 * t1276 * t1310 + 24.0 * t402 * t2299 + 4.0 * t1748 * t1132;
    double t2314 = t180 * t835 * t58;
    double t2322 = t152 * t1947;
    double t2346 = 4.0 * t1609 * t835 * t115 * t63 - 4.0 * t724 * t233 * t108 * t1310 - 8.0 * t2314 * t1026 - 72.0 * t1515 * t1024 * t831 + 8.0 * t2314 * t847 + 24.0 * t218 * t2322 - 16.0 * t214 * t2299 - 4.0 * t1238 * t889 + 8.0 * t146 * t319 * t1917 * t830 + 32.0 * t677 * t208 * t56 * t831 + 8.0 * t2190 * t334 * t809 * t1144 + 4.0 * t1718 * t1131 * t509 + t525 * t2169 * t115;
    double t2366 = t234 * t1030;
    double t2367 = t238 * t1144;
    double t2378 = t22 * eta;
    double t2391 = -72.0 * t1472 * t979 * t831 - 4.0 * t1955 * t1094 + 32.0 * t289 * t63 * t1381 + 16.0 * t646 * t1612 * t831 + 32.0 * t1031 * t347 * t1033 + 8.0 * t677 * t1290 * t941 + 8.0 * t2366 * t2367 - 24.0 * t152 * t52 * t1126 * t813 - 32.0 * t319 * t810 * t1054 + 8.0 * t325 * t2046 * t2378 * t1184 + 4.0 * t488 * t2203 * t1260 - 24.0 * t1580 * t1053 * t604 + 4.0 * t1971 * t1139 * t813;
    double t2392 = t56 * t235;
    double t2398 = t145 * t3;
    double t2409 = t1213 * t22;
    double t2432 = t235 * t371;
    double t2434 = t144 * t205;
    double t2435 = t2434 * t1047;
    double t2442 = -32.0 * t146 * t2392 * t138 * t336 * t1047 - t2398 * eta * t108 * t115 - 8.0 * t1797 * t327 * t110 - 4.0 * t2398 * eta * t1298 * t58 + 4.0 * t17 * t2409 + 2.0 * t228 * t1053 * t32 + 8.0 * t2366 * t938 + 8.0 * t17 * t2378 * t145 * t56 * t111 + 16.0 * t509 * t2378 * t1195 + 8.0 * t529 * t2378 * t1214 + 16.0 * t1800 * t43 * t56 * t831 + 16.0 * t2238 * t2432 * t2435 + 4.0 * t646 * t251 * t108 * t941;
    double t2448 = t144 * t58 * t846;
    double t2463 = t2434 * t88;
    double t2485 = 32.0 * t1057 * t2392 * t1060 - 8.0 * t2366 * t2448 + 16.0 * t2260 * t2432 * t2435 + 8.0
                   * t1071 * t2367 + 8.0 * t1071 * t938 + 32.0 * t334 * t352 * t2284 + 4.0 * t43 * t369 * t2190 * t2463 - 8.0 * t1071 * t2448 + 8.0 * t2274 * t2190 * t2463 + 4.0 * t251 * t145 * t2190 * t2463 + 4.0 * t1955 * t851 + 8.0 * t828 * t810 * t809 * t98 - 2.0 * t2180 * t2046 * t1192 * t98;
    double t2506 = t180 * t835 * t235;
    double t2513 = t180 * t1980;
    double t2520 = t146 * t1980;
    double t2531 = -32.0 * t204 * t41 * t334 * t2284 + 48.0 * t366 * t1049 + 48.0 * t237 * t1049 - 32.0 * t147 * t2322 - 16.0 * t250 * t310 * t1947 + 8.0 * t147 * t385 * t1184 + 4.0 * t250 * t2226 * t1184 + 8.0 * t2506 * t1522 * t1144 - 8.0 * t2506 * t1522 * t937 + 8.0 * t2513 * t1024 * t1408 - 12.0 * t2513 * t1024 * t941 + 8.0 * t2520 * t979 * t1408 - 12.0 * t2520 * t979 * t941 - 64.0 * t180 * t347 * t1522 * t1048;
    double t2537 = t129 * (t2212 + t2258 + t2304 + t2346 + t2391 + t2442 + t2485 + t2531) * t1338 / 2.0;
    double t2589 = 2.0 * t327 + 2.0 * t325 * t159 - 2.0 * t1797 * t1342 + 8.0 * t341 * t106 * t1559 + 4.0 * t341 * t106 * t1567 + t492 * t1843 * t208 + 4.0 * t492 * t182 * t138 + 2.0 * t341 * t1843 * t43 + 4.0 * t341 * t182 * r + 4.0 * t514 * t354 * t509 + 12.0 * t514 * t726 * t637 - 4.0 * t52 * t1220 * t1856 + t195 * t762 + 4.0 * t591 * t354 * t17 - 4.0 * t2022 * t1879 + 4.0 * t591 * t726 * M + 4.0 * t596 * t239 * t63 + 8.0 * t596 * t592 * t54 - 8.0 * t341 * t307 - 8.0 * t596 * t1805;
    double t2625 = 2.0 * t1201 * t22 + 2.0 * t2409 + 8.0 * t1797 * t2017 - 2.0 * t2181 * t984 * t1 + 16.0 * t52 * t1298 * t806 + 4.0 * t492 * t1122 * t299 + 8.0 * t341 * t1122 * t164 + 24.0 * t52 * t1126 * t674 + 4.0 * t358 * t1122 * t58 + 24.0 * t195 * t1126 * t726 * t17 + 32.0 * t195 * t204 * t110 * t592 * t63;
    double t2646 = t195 * t835;
    double t2655 = 4.0 * t1315 * t1177 + t1309 * t1310 * t208 + 2.0 * t1315 * t1311 + 4.0 * t1318 * t941 * t509 - 4.0 * t1259 * t515 * t873 + t195 * t1166 + 4.0 * t195 * t976 * t1319 - 4.0 * t195 * t976 * t1 * t1260 +
                   4.0 * t2646 * t937 * t63 - 8.0 * t1315 * t864 - 8.0 * t2646 * t1144 * t63;
    christoffel[0][0][0] = 0.0;
    christoffel[0][0][1] = -t7 * t14 * t77 / 2.0;
    christoffel[0][0][2] = -t86 * (-2.0 * t89 + t9 * t92 - 4.0 * t26 * t95 + 4.0 * t31 * t99 + 4.0 * t31 * t103 + 16.0 * t107 * t112 + 4.0 * t53 * t63 * t116 - 8.0 * t9 * t87 * t70 - 8.0 * t52 * t87 * t74) / 2.0;
    christoffel[0][0][3] = 0.0;
    christoffel[0][1][0] = t633;
    christoffel[0][1][1] = 0.0;
    christoffel[0][1][2] = 0.0;
    christoffel[0][1][3] = t827;
    christoffel[0][2][0] = t1108;
    christoffel[0][2][1] = 0.0;
    christoffel[0][2][2] = 0.0;
    christoffel[0][2][3] = t1340;
    christoffel[0][3][0] = 0.0;
    christoffel[0][3][1] = -t1376;
    christoffel[0][3][2] = -t1420;
    christoffel[0][3][3] = 0.0;
    christoffel[1][0][0] = t633;
    christoffel[1][0][1] = 0.0;
    christoffel[1][0][2] = 0.0;
    christoffel[1][0][3] = t827;
    christoffel[1][1][0] = 0.0;
    christoffel[1][1][1] = t1421 / t5 * t14 * (t18 * t6 - t1 * t1425 * t5 + 2.0 * t31 * t10 * t18 * t5 - 2.0 * t9 * t2 * t10 * t1425 + t9 * t2 * t1436 * t5) / 2.0;
    christoffel[1][1][2] = -t1421 * t85 * t1452 / t6 / 2.0;
    christoffel[1][1][3] = 0.0;
    christoffel[1][2][0] = 0.0;
    christoffel[1][2][1] = t1459;
    christoffel[1][2][2] = t1469;
    christoffel[1][2][3] = 0.0;
    christoffel[1][3][0] = t1712;
    christoffel[1][3][1] = 0.0;
    christoffel[1][3][2] = 0.0;
    christoffel[1][3][3] = t1894;
    christoffel[2][0][0] = t1108;
    christoffel[2][0][1] = 0.0;
    christoffel[2][0][2] = 0.0;
    christoffel[2][0][3] = t1340;
    christoffel[2][1][0] = 0.0;
    christoffel[2][1][1] = t1459;
    christoffel[2][1][2] = t1469;
    christoffel[2][1][3] = 0.0;
    christoffel[2][2][0] = 0.0;
    christoffel[2][2][1] = -t6 * t1466 * t1457 / 2.0;
    christoffel[2][2][2] = (t87 + 2.0 * t9 * t82 * t87 + t9 * t2 * t1901) * t1421 * t85 / 2.0;
    christoffel[2][2][3] = 0.0;
    christoffel[2][3][0] = -t2168;
    christoffel[2][3][1] = 0.0;
    christoffel[2][3][2] = 0.0;
    christoffel[2][3][3] = t2537;
    christoffel[3][0][0] = 0.0;
    christoffel[3][0][1] = -t1376;
    christoffel[3][0][2] = -t1420;
    christoffel[3][0][3] = 0.0;
    christoffel[3][1][0] = t1712;
    christoffel[3][1][1] = 0.0;
    christoffel[3][1][2] = 0.0;
    christoffel[3][1][3] = t1894;
    christoffel[3][2][0] = -t2168;
    christoffel[3][2][1] = 0.0;
    christoffel[3][2][2] = 0.0;
    christoffel[3][2][3] = t2537;
    christoffel[3][3][0] = 0.0;
    christoffel[3][3][1] = -t7 * t14 * t56 * t2589 / 2.0;
    christoffel[3][3][2] = -t86 * t55 * (t2625 + t2655) / 2.0;
    christoffel[3][3][3] = 0.0;
    return true;
}

/*! Calculate Jacobi matrix.
 *
 *  \param pos : pointer to position.
 */
bool MetricHartleThorneGB::calculateChrisD(const double*) {
    return false;
}

/*! Transform local 4-direction to coordinate 4-direction.
 *
 *  \param  pos  :  pointer to position array.
 *  \param  ldir :  pointer to local direction array.
 *  \param  dir  :  pointer to calculated coordinate direction array.
 *  \param  type :  type of tetrad.
 */
void MetricHartleThorneGB::localToCoord(const double* pos, const double* ldir, double* dir,
                                        enum_nat_tetrad_type) {
    calculateMetric(pos);

    double gtt = g_compts[0][0];
    double grr = g_compts[1][1];
    double gaa = g_compts[2][2];
    double gtp = g_compts[0][3];
    double gpp = g_compts[3][3];

    //double omega = -gtp/gpp;
    double zeta = 0.0;

    double gam = 1.0 / sqrt(-(gtt + 2.0 * zeta * gtp + zeta * zeta * gpp));
    double dlt = 1.0 / sqrt(gtp * gtp - gtt * gpp);
    double w1 = gtp + zeta * gpp;
    double w2 = gtt + zeta * gtp;

    dir[0] = gam * (ldir[0] + dlt * w1 * ldir[3]);
    dir[1] = ldir[1] / sqrt(grr);
    dir[2] = ldir[2] / sqrt(gaa);
    dir[3] = gam * (ldir[0] * zeta - dlt * w2 * ldir[3]);

    //double prod;
    //calcProduct(pos,dir,dir,prod);
    //fprintf(stderr,"hier: %e\n",prod);
}

/*! Transform coordinate 4-direction to local 4-direction.
 *
 *  \param  pos  :  pointer to position array.
 *  \param  cdir :  pointer to coordinate direction.
 *  \param  ldir :  pointer to calculated local direction array.
 *  \param  type :  type of tetrad.
 */
void MetricHartleThorneGB::coordToLocal(const double* , const double* , double* ,
                                        enum_nat_tetrad_type) {
    // TODO
}


/*! Test break condition.
 *
 *  \param pos    : pointer to position array.
 *  \return true  : radial position r < 0.0 or  r^2<=(1.0+eps)*rs^2.
 *  \return false : position is valid.
 */
bool MetricHartleThorneGB::breakCondition(const double* pos) {
    bool br = false;
    if (pos[1] <= 2.0 * mMass || isnan(pos[0]) || isnan(pos[1]) || isnan(pos[2]) || isnan(pos[3])) {
        br = true;
    } else {
        calculateMetric(pos);

        double gtt = g_compts[0][0];
        double grr = g_compts[1][1];
        double gaa = g_compts[2][2];
        double gtp = g_compts[0][3];
        double gpp = g_compts[3][3];

        double dd = gtp * gtp - gtt * gpp;

        if (grr <= 0.0 || gaa <= 0.0 || dd <= 0.0) {
            br = true;
        }
    }
    return br;
}


/*! Set parameter 'pName' to 'val'.
 *
 *  Set 'mass' or 'lambda' parameter.
 */
bool MetricHartleThorneGB::setParam(std::string pName, double val) {
    Metric::setParam(pName, val);
    if (pName == "mass") {
        mMass = val;
    } else if (pName == "angmom") {
        mAngmom = val;
    } else if (pName == "eta") {
        mEta = val;
    }
    return true;
}

/*! Generate report.
 */
bool MetricHartleThorneGB::report(const vec4 , const vec4 , std::string &text) {
    std::stringstream ss;
    ss << "Report for HartleThorne metric\n\tcoordinate : (t,x,y,z)\n";
    ss << "---------------------------------------------------------------\n";
    ss << "  physical units ..................... no\n";
    ss.precision(DEF_FIXED_REPORT_PRECISION);
    ss.setf(std::ios::fixed);
    ss << "  param .............................. mass   = " << mMass << std::endl;
    ss << "  param .............................. angmom = " << mAngmom << std::endl;
    ss << "  param .............................. eta    = " << mEta << std::endl;

    text = ss.str();
    return true;
}

// ********************************* protected methods *****************************
void MetricHartleThorneGB::calcFunc(const double* pos) {
    double M = mMass;
    double a = mAngmom;

    double r = pos[1];
    double theta = pos[2];

    double t1 = cos(theta);
    double t2 = t1 * t1;
    double t4 = 1.0 - 3.0 * t2;
    double t6 = 1 / M;
    double t8 = 1 / r;
    double t10 = r - 2.0 * M;
    double t11 = 1 / t10;
    double t13 = M * M;
    double t14 = 2.0 * t13;
    double t15 = M * r;
    double t17 = r * r;
    double t18 = 3.0 * t17;
    double t24 = 1 / t13;
    double t26 = log(r * t11);
    double t30 = -5.0 / 8.0 * (r - M) * t6 * t8 * t11 * (t14 + 6.0 * t15 - t18) - 15.0 / 16.0 * r * t10 * t24 * t26;
    double t34 = 1.0 - 2.0 * M * t8;
    double t39 = 1 / t17;
    double t50 = 5.0 / 8.0 * t6 * t8 * (t14 - 3.0 * t15 - t18) + 15.0 / 16.0 * t24 * (t17 - t14) * t26;
    double t52 = sin(theta);
    double t53 = t52 * t52;

    htt = t4 * t30 / t34;
    hrr = t34 * t4 * t30;
    hthth = -t39 * t4 * t50;
    hphph = -1 / t53 * t39 * t4 * t50;

    sigma = t17 + a * a * t2;
    delta = t17 - 2.0 * t15 + a * a;
}

void MetricHartleThorneGB::calcFuncDiff(const double* pos) {
    double M = mMass;
    double a = mAngmom;

    double r = pos[1];
    double theta = pos[2];

    double t1 = cos(theta);
    double t2 = t1 * t1;
    double t4 = 1.0 - 3.0 * t2;
    double t6 = 1 / M;
    double t8 = 1 / r;
    double t10 = r - 2.0 * M;
    double t11 = 1 / t10;
    double t13 = M * M;
    double t14 = 2.0 * t13;
    double t15 = M * r;
    double t17 = r * r;
    double t18 = 3.0 * t17;
    double t24 = 1 / t13;
    double t26 = log(r * t11);
    double t30 = -5.0 / 8.0 * (r - M) * t6 * t8 * t11 * (t14 + 6.0 * t15 - t18) - 15.0 / 16.0 * r * t10 * t24 * t26;
    double t34 = 1.0 - 2.0 * M * t8;
    double t39 = 1 / t17;
    double t50 = 5.0 / 8.0 * t6 * t8 * (t14 - 3.0 * t15 - t18) + 15.0 / 16.0 * t24 * (t17 - t14) * t26;
    double t52 = sin(theta);
    double t53 = t52 * t52;
    double t58 = t17 * t17;
    double t60 = 3.0 * t26 * t58;
    double t61 = t13 * t13;
    double t63 = t13 * t17;
    double t65 = t13 * M;
    double t66 = t65 * r;
    double t68 = t17 * r;
    double t70 = 6.0 * M * t68;
    double t71 = t26 * t68;
    double t72 = t71 * M;
    double t74 = t26 * r;
    double t77 = t26 * t17;
    double t78 = t77 * t13;
    double t82 = t10 * t10;
    double t88 = t1 * t52;
    double t91 = 4.0 * t61;
    double t94 = -8.0 * t66 + 18.0 * t63 - t70 - t91 + t60 - 12.0 * t72 + 12.0 * t78;
    double t105 = 1 / t68;
    double t118 = -t4 * (-2.0 * t15 + t14 + t77 - 2.0 * t74 * M);
    double t119 = 1 / t58;
    double t131 = 4.0 * t65 - 6.0 * t13 * r - 6.0 * M * t17 + 3.0 * t71 - 6.0 * t74 * t13;
    double t137 = 1 / (-1.0 + t2);

    htt = t4 * t30 / t34;
    hrr = t34 * t4 * t30;
    hthth = -t39 * t4 * t50;
    hphph = -1 / t53 * t39 * t4 * t50;

    dhttdr  = -5.0 / 8.0 * t4 * (t60 + 12.0 * t61 + 30.0 * t63 - 44.0 * t66 - t70 - 18.0 * t72 - 24.0 * t74 * t65 + 36.0 * t78) * t24 / t82 / t10;
    dhttdth = -15.0 / 8.0 * t88 * t94 * t24 / t82;
    dhrrdr  = -5.0 / 8.0 * (t60 - 6.0 * t72 - t70 + 6.0 * t63 + 4.0 * t66 + t91) * t4 * t24 * t105;
    dhrrdth = -15.0 / 8.0 * t94 * t52 * t1 * t24 * t39;
    dhaadr  = 15.0 / 4.0 * t118 * t119 * t11;
    dhaadth = -15.0 / 8.0 * t88 * t131 * t105 * t24;
    dhppdr  = -15.0 / 4.0 * t118 * t137 * t119 * t11;
    dhppdth = 5.0 / 4.0 * t1 * t131 * t137 / t52 * t105 * t24;

    sigma = t17 + a * a * t2;
    delta = t17 - 2.0 * t15 + a * a;

    dsigmadr = 2.0 * r;
    dsigmadth = -2.0 * a * a * cos(theta) * sin(theta);
    deltadr = 2.0 * r - 2.0 * M;
}

/*!
 */
void MetricHartleThorneGB::setStandardValues() {
    mInitPos[0] = 0.0;
    mInitPos[1] = 10.0;
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

} // end namespace m4d
