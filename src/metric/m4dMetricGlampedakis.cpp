// -------------------------------------------------------------------------------
/*
   m4dMetricGlampedakis.cpp

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
#include "m4dMetricGlampedakis.h"

namespace m4d {


/*! Standard constructor for the Glampedakis metric.
 *
 */
MetricGlampedakis::MetricGlampedakis(double mass, double angmom, double epsilon) {
    mMetricName  = "Glampedakis";
    setCoordType(enum_coordinate_spherical);

    mPhysicalUnits = enum_physical_constants_geom;
    mSpeedOfLight = 1.0;
    mGravConstant = 1.0;

    addParam("mass", mass);
    mMass = mass;
    M = mMass;
    addParam("angmom", angmom);
    mAngmom = angmom;
    a   = mAngmom;
    a2  = a * a;
    addParam("epsilon", epsilon);
    mEpsilon = epsilon;

    rs = 2.0 * mMass;

    setStandardValues();
    initToZero();

    metricskip = 0;
    metriccalled = 0;

    christoffelskip = 0;
    christoffelcalled = 0;

    chrisDskip = 0;
    chrisDcalled = 0;

    kerrskip = 0;
    kerrcalled = 0;

    kerrdiffskip = 0;
    kerrdiffcalled = 0;

    kerrdiff2skip = 0;
    kerrdiff2called = 0;

    glampskip = 0;
    glampcalled = 0;

    glampdiffskip = 0;
    glampdiffcalled = 0;

    glampdiff2skip = 0;
    glampdiff2called = 0;

}

MetricGlampedakis::~MetricGlampedakis() {
    std::cout << " kerr: " << kerrskip << "/" << kerrcalled << ", diff: " <<  kerrdiffskip << "/" << kerrdiffcalled
              << ", diff2: " <<  kerrdiff2skip << "/" << kerrdiff2called << "\n";
    std::cout << " glamp: " << glampskip << "/" << glampcalled
              << ", diff: " <<  glampdiffskip << "/" << glampdiffcalled
              << ", diff2: " <<  glampdiff2skip << "/" << glampdiff2called << "\n";
    std::cout << "metric: " << metricskip << "/" << metriccalled << "\n";
    std::cout << "christoffels: " << christoffelskip << "/" << christoffelcalled << "\n";
    std::cout << "chrisD: " << chrisDskip << "/" << chrisDcalled << "\n";
}

#include <iostream>
//#include <iomanip>      // std::setprecision


// *********************************** public methods ******************************
/*! Calculate the contravariant metric components at position 'pos'.
 *
 *  \param pos : pointer to position.
 */
bool MetricGlampedakis::calculateMetric(const double* pos) {

    metriccalled++;
    if ((oldmetricpos[0] == pos[1]) && (oldmetricpos[1] == pos[2])) {
        metricskip++;
        return true;
    }
    oldmetricpos[0] = pos[1];
    oldmetricpos[0] = pos[2];

    calcgComps(pos);

    g_compts[0][0] = gtt;
    g_compts[0][3] = gtph;
    g_compts[1][1] = grr;
    g_compts[2][2] = gthth;
    g_compts[3][0] = gtph;
    g_compts[3][3] = gphph;
    /*
        std::cout.precision(12);
        std::cout << "Metrik:\n";
        std::cout << "tt:\t "   << g_compts[0][0] << "\n";
        std::cout << "tph:\t "  << g_compts[0][3] << "\n";
        std::cout << "rr:\t "   << g_compts[1][1] << "\n";
        std::cout << "thth:\t " << g_compts[2][2] << "\n";
        std::cout << "phph:\t " << g_compts[3][3] << "\n";
    */
    return true;
}

/*! Calculate the Christoffel symbols of the second kind at position 'pos'.
 *
 *  \param pos : pointer to position.
 */
bool MetricGlampedakis::calculateChristoffels(const double* pos) {
    christoffelcalled++;
    if ((oldchristoffelpos[0] == pos[1]) && (oldchristoffelpos[1] == pos[2])) {
        christoffelskip++;
        return true;
    }
    oldchristoffelpos[0] = pos[1];
    oldchristoffelpos[1] = pos[2];


    calcgComps(pos);
    calcgCompsDiff(pos);
    double gtph2 = gtph * gtph;
    double fac = (gphph * gtt - gtph2);
    //111  1
    christoffel[1][1][1] = 0.5 * grrdr / grr;
    //112  2
    christoffel[1][1][2] = -0.5 * grrdth / gthth;
    //121  3
    christoffel[1][2][1] = 0.5 * grrdth / grr;
    christoffel[2][1][1] = christoffel[1][2][1];
    //122  4
    christoffel[1][2][2] = 0.5 * gththdr / gthth;
    christoffel[2][1][2] = christoffel[1][2][2];
    //133  5
    christoffel[1][3][3] = 0.5 * (gtt * gphphdr - gtph * gtphdr) / fac;
    christoffel[3][1][3] = christoffel[1][3][3];
    //130  6
    christoffel[1][3][0] = 0.5 * (-gtph * gphphdr + gphph * gtphdr) / fac;
    christoffel[3][1][0] = christoffel[1][3][0];
    //103  7
    christoffel[1][0][3] = 0.5 * (gtt * gtphdr - gtph * gttdr) / fac;
    christoffel[0][1][3] = christoffel[1][0][3];
    //100  8
    christoffel[1][0][0] = 0.5 * (-gtph * gtphdr + gphph * gttdr) / fac;
    christoffel[0][1][0] = christoffel[1][0][0];
    //221  9
    christoffel[2][2][1] = -0.5 * gththdr / grr;
    //222  10
    christoffel[2][2][2] = 0.5 * gththdth / gthth;
    //233  11
    christoffel[2][3][3] = 0.5 * (gtt * gphphdth - gtph * gtphdth) / fac;
    christoffel[3][2][3] = christoffel[2][3][3];
    //230  12
    christoffel[2][3][0] = 0.5 * (-gtph * gphphdth + gphph * gtphdth) / fac;
    christoffel[3][2][0] = christoffel[2][3][0];
    //203  13
    christoffel[2][0][3] = 0.5 * (gtt * gtphdth - gtph * gttdth) / fac;
    christoffel[0][2][3] = christoffel[2][0][3];
    //200  14
    christoffel[2][0][0] = 0.5 * (-gtph * gtphdth + gphph * gttdth) / fac;
    christoffel[0][2][0] = christoffel[2][0][0];
    //331  15
    christoffel[3][3][1] = -0.5 * gphphdr / grr;
    //332  16
    christoffel[3][3][2] = -0.5 * gphphdth / gthth;
    //301  17
    christoffel[3][0][1] = -0.5 * gtphdr / grr;
    christoffel[0][3][1] = christoffel[3][0][1];
    //302  18
    christoffel[3][0][2] = -0.5 * gtphdth / gthth;
    christoffel[0][3][2] = christoffel[3][0][2];
    //001  19
    christoffel[0][0][1] = -0.5 * gttdr / grr;
    //002  20
    christoffel[0][0][2] = -0.5 * gttdth / gthth;
    /*
        std::cout << "Christoffel:\n";
        std::cout << "rrr:\t " << christoffel[1][1][1] << "\n";
        std::cout << "rrth:\t " << christoffel[1][1][2] << "\n";
        std::cout << "rthr:\t " << christoffel[1][2][1] << "\n";
        std::cout << "rthth:\t " << christoffel[1][2][2] << "\n";
        std::cout << "rphph:\t " << christoffel[1][3][3] << "\n";
        std::cout << "rpht:\t " << christoffel[1][3][0] << "\n";
        std::cout << "rtph:\t " << christoffel[1][0][3] << "\n";
        std::cout << "rtt:\t " << christoffel[1][0][0] << "\n";
        std::cout << "ththr:\t " << christoffel[2][2][1] << "\n";
        std::cout << "ththth:\t " << christoffel[2][2][2] << "\n";
        std::cout << "thphph:\t " << christoffel[2][3][3] << "\n";
        std::cout << "thpht:\t " << christoffel[2][3][0] << "\n";
        std::cout << "thtph:\t " << christoffel[2][0][3] << "\n";
        std::cout << "thtt:\t " << christoffel[2][0][0] << "\n";
        std::cout << "phphr:\t " << christoffel[3][3][1] << "\n";
        std::cout << "phphth:\t " << christoffel[3][3][2] << "\n";
        std::cout << "phtr:\t " << christoffel[3][0][1] << "\n";
        std::cout << "phtth:\t " << christoffel[3][0][2] << "\n";
        std::cout << "ttr:\t " << christoffel[0][0][1] << "\n";
        std::cout << "ttth:\t " << christoffel[0][0][2] << "\n";
    */
    return true;
}

/*! Calculate Jacobi matrix.
 *
 *  \param pos : pointer to position.
 */
bool MetricGlampedakis::calculateChrisD(const double* pos) {
    chrisDcalled++;
    if ((oldchrisDpos[0] == pos[1]) && (oldchrisDpos[1] == pos[2])) {
        chrisDskip++;
        return true;
    }
    oldchrisDpos[0] = pos[1];
    oldchrisDpos[1] = pos[2];

    calcgCompsDiff2(pos);

    double gtt2   = gtt * gtt;
    double gtph2    = gtph * gtph;
    double grr2     = grr * grr;
    double gthth2   = gthth * gthth;
    double gphph2   = gphph * gphph;

    double gttdr2   = gttdr * gttdr;
    double gtphdr2  = gtphdr * gtphdr;
    double grrdr2   = grrdr * grrdr;
    double grrdth2  = grrdth * grrdth;
    double gththdr2 = gththdr * gththdr;
    double gphphdr2 = gphphdr * gphphdr;

    double gttdth2  = gttdth * gttdth;
    double gtphdth2 = gtphdth * gtphdth;
    double gththdth2 = gththdth * gththdth;
    double gphphdth2 = gphphdth * gphphdth;

    double gtph3    = gtph2 * gtph;

    double fac2 = (-gphph * gtt + gtph2) * (-gphph * gtt + gtph2);


    chrisD[1][1][1][1] = 0.5 * (-grrdr2 + grrd2r * grr) / grr2;

    chrisD[1][1][1][2] = 0.5 * (-grrdr * grrdth + grrdrdth * grr) / grr2;

    chrisD[1][1][2][1] = 0.5 * (grrdth * gththdr - grrdrdth * gthth) / gthth2;

    chrisD[1][1][2][2] = 0.5 * (grrdth * gththdth - grrd2th * gthth) / gthth2;

    chrisD[1][2][1][1] = 0.5 * (-grrdr * grrdth + grrdrdth * grr) / grr2;
    chrisD[2][1][1][1] = chrisD[1][2][1][1];

    chrisD[1][2][1][2] = 0.5 * (-grrdth2 + grrd2th * grr) / grr2;
    chrisD[2][1][1][2] = chrisD[1][2][1][2];

    chrisD[1][2][2][1] = 0.5 * (-gththdr2 + gththd2r * gthth) / gthth2;
    chrisD[2][1][2][1] = chrisD[1][2][2][1];

    chrisD[1][2][2][2] = 0.5 * (-gththdr * gththdth + gththdrdth * gthth) / gthth2;
    chrisD[2][1][2][2] = chrisD[1][2][2][2];

    chrisD[1][3][3][1] = -0.5 * (gttdr * gphphdr * gtph2 - gtt2 * gphphd2r * gphph + gtt * gphphd2r * gtph2 + gtphdr2 * gphph * gtt + gtphdr2 * gtph2
                                 + gtph * gtphd2r * gphph * gtt - gtph3 * gtphd2r + gtt2 * gphphdr2 - 3 * gtt * gphphdr * gtph * gtphdr - gtph * gtphdr * gphph * gttdr) / fac2;
    chrisD[3][1][3][1] = chrisD[1][3][3][1];

    chrisD[1][3][3][2] = -0.5 * (gttdth * gphphdr * gtph2 - gtt2 *  gphphdrdth * gphph + gtt * gphphdrdth * gtph2 + gtphdth * gtphdr * gphph * gtt + gtphdth * gtphdr * gtph2
                                 + gtph * gtphdrdth * gphph * gtt - gtph3 * gtphdrdth + gtt2 * gphphdr * gphphdth - 2 * gtt * gphphdr * gtph * gtphdth - gtph * gtphdr * gtt * gphphdth
                                 - gtph * gtphdr * gphph * gttdth) / fac2;
    chrisD[3][1][3][2] = chrisD[1][3][3][2];

    chrisD[1][3][0][1] = 0.5 * (-gtph * gphphd2r * gphph * gtt + gtph3 * gphphd2r + gphph2 * gtphd2r * gtt - gphph * gtphd2r * gtph2 + gtph * gphphdr2 * gtt + gtph * gphphdr * gphph * gttdr
                                - 2 * gtph2 * gphphdr * gtphdr - gphph * gtphdr * gtt * gphphdr - gphph2 * gtphdr * gttdr + 2 * gphph * gtphdr2 * gtph) / fac2;
    chrisD[3][1][0][1] = chrisD[1][3][0][1];


    chrisD[1][3][0][2] = -0.5 * (gtphdth * gphphdr * gphph * gtt + gtphdth * gphphdr * gtph2 + gtph * gphphdrdth * gphph * gtt - gtph3 * gphphdrdth + gphphdth * gtphdr * gtph2
                                 - gphph2 * gtphdrdth * gtt + gphph *  gtphdrdth * gtph2 - gtph * gphphdr * gtt * gphphdth - gtph * gphphdr * gphph * gttdth + gphph2 * gtphdr * gttdth
                                 - 2 * gphph * gtphdr * gtph * gtphdth) / fac2;
    chrisD[3][1][0][2] = chrisD[1][3][0][2];

    chrisD[1][0][3][1] = 0.5 * (gtt2 * gtphd2r * gphph - gtt * gtphd2r * gtph2 - gtph * gttd2r * gphph * gtt + gtph3 * gttd2r - gtt2 * gtphdr * gphphdr - gtt * gtphdr * gphph * gttdr
                                + 2 * gtt * gtphdr2 * gtph + gtph * gttdr * gtt * gphphdr + gtph * gttdr2 * gphph - 2 * gtph2 * gttdr * gtphdr) / fac2;
    chrisD[0][1][3][1] = chrisD[1][0][3][1];

    chrisD[1][0][3][2] = -0.5 * (gttdth * gtphdr * gtph2 - gtt2 * gtphdrdth * gphph + gtt * gtphdrdth * gtph2 + gtphdth * gttdr * gphph * gtt + gtphdth * gttdr * gtph2
                                 + gtph * gttdrdth * gphph * gtt - gtph3 * gttdrdth + gtt2 * gtphdr * gphphdth - 2 * gtt * gtphdr * gtph * gtphdth - gtph * gttdr * gtt * gphphdth
                                 - gtph * gttdr * gphph * gttdth) / fac2;
    chrisD[0][1][3][2] = chrisD[1][0][3][2];

    chrisD[1][0][0][1] = -0.5 * (gtphdr2 * gphph * gtt + gtphdr2 * gtph2 + gtph * gtphd2r * gphph * gtt - gtph3 * gtphd2r + gttdr * gphphdr * gtph2 - gphph2 * gttd2r * gtt + gphph * gttd2r * gtph2
                                 - gtt * gphphdr * gtph * gtphdr - 3 * gtph * gtphdr * gphph * gttdr + gphph2 * gttdr2) / fac2;
    chrisD[0][1][0][1] = chrisD[1][0][0][1];

    chrisD[1][0][0][2] = -0.5 * (gtphdth * gtphdr * gphph * gtt + gtphdth * gtphdr * gtph2 + gtph * gtphdrdth * gphph * gtt - gtph3 * gtphdrdth + gphphdth * gttdr * gtph2
                                 - gphph2 * gttdrdth * gtt + gphph * gttdrdth * gtph2 - gtph * gtphdr * gtt * gphphdth - gtph * gtphdr * gphph * gttdth + gphph2 * gttdr * gttdth
                                 - 2 * gphph * gttdr * gtph * gtphdth) / fac2;
    chrisD[0][1][0][2] = chrisD[1][0][0][2];

    chrisD[2][2][1][1] = -0.5 * (-gththdr * grrdr + gththd2r * grr) / grr2;

    chrisD[2][2][1][2] = -0.5 * (-grrdth * gththdr + gththdrdth * grr) / grr2;

    chrisD[2][2][2][1] = 0.5 * (-gththdr * gththdth + gththdrdth * gthth) / gthth2;

    chrisD[2][2][2][2] = 0.5 * (-gththdth2 + gththd2th * gthth) / gthth2;

    chrisD[2][3][3][1] = -0.5 * (gphphdth * gttdr * gtph2 - gtt2 * gphphdrdth * gphph + gtt * gphphdrdth * gtph2 + gtphdth * gtphdr * gphph * gtt + gtphdth * gtphdr * gtph2
                                 + gtph * gtphdrdth * gphph * gtt - gtph3 * gtphdrdth + gtt2 * gphphdr * gphphdth - 2 * gtph * gtphdr * gtt * gphphdth - gtt * gphphdr * gtph * gtphdth
                                 - gphph * gttdr * gtph * gtphdth) / fac2;
    chrisD[3][2][3][1] = chrisD[2][3][3][1];

    chrisD[2][3][3][2] = 0.5 * (-gttdth * gphphdth * gtph2 + gtt2 * gphphd2th * gphph - gtt * gphphd2th * gtph2 - gtphdth2 * gphph * gtt - gtphdth2 * gtph2
                                - gtph * gtphd2th * gphph * gtt + gtph3 * gtphd2th - gtt2 * gphphdth2 + 3 * gtt * gphphdth * gtph * gtphdth
                                + gtph * gtphdth * gphph * gttdth) / fac2;
    chrisD[3][2][3][2] = chrisD[2][3][3][2];

    chrisD[2][3][0][1] = 0.5 * (-gphphdth * gtphdr * gphph * gtt - gphphdth * gtphdr * gtph2 - gtph * gphphdrdth * gphph * gtt + gtph3 * gphphdrdth - gtphdth * gphphdr * gtph2
                                + gphph2 * gtphdrdth * gtt - gphph * gtphdrdth * gtph2 + gtph * gphphdr * gtt * gphphdth + gtph * gphphdth * gphph * gttdr - gphph2 * gtphdth * gttdr
                                + 2 * gphph * gtphdr * gtph * gtphdth) / fac2;
    chrisD[3][2][0][1] = chrisD[2][3][0][1];


    chrisD[2][3][0][2] = 0.5 * (-gtph * gphphd2th * gphph * gtt + gtph3 * gphphd2th + gphph2 * gtphd2th * gtt - gphph * gtphd2th * gtph2 + gtph * gphphdth2 * gtt
                                + gtph * gphphdth * gphph * gttdth - 2 * gtph2 * gphphdth * gtphdth - gphph * gtphdth * gtt * gphphdth - gphph2 * gtphdth * gttdth
                                + 2 * gphph * gtphdth2 * gtph) / fac2;
    chrisD[3][2][0][2] = chrisD[2][3][0][2];

    chrisD[2][0][3][1] = -0.5 * (gtphdth * gttdr * gtph2 - gtt2 * gtphdrdth * gphph + gtt * gtphdrdth * gtph2 + gttdth * gtphdr * gphph * gtt + gttdth * gtphdr * gtph2
                                 + gtph * gttdrdth * gphph * gtt - gtph3 * gttdrdth + gtt2 * gtphdth * gphphdr - 2 * gtt * gtphdr * gtph * gtphdth - gtph * gttdth * gtt * gphphdr
                                 - gtph * gttdr * gphph * gttdth) / fac2;
    chrisD[0][2][3][1] = chrisD[2][0][3][1];


    chrisD[2][0][3][2] = 0.5 * (gtt2 * gtphd2th * gphph - gtt * gtphd2th * gtph2 - gtph * gttd2th * gphph * gtt + gtph3 * gttd2th - gtt2 * gtphdth * gphphdth
                                - gtt * gtphdth * gphph * gttdth + 2 * gtt * gtphdth2 * gtph + gtph * gttdth * gtt * gphphdth + gtph * gttdth2 * gphph
                                - 2 * gtph2 * gttdth * gtphdth) / fac2;
    chrisD[0][2][3][2] = chrisD[2][0][3][2];

    chrisD[2][0][0][1] = 0.5 * (-gtphdth * gtphdr * gphph * gtt - gtphdth * gtphdr * gtph2 - gtph * gtphdrdth * gphph * gtt + gtph3 * gtphdrdth - gttdth * gphphdr * gtph2
                                + gphph2 * gttdrdth * gtt - gphph * gttdrdth * gtph2 + gtt * gphphdr * gtph * gtphdth + gphph * gttdr * gtph * gtphdth - gphph2 * gttdr * gttdth
                                + 2 * gtph * gtphdr * gphph * gttdth) / fac2;
    chrisD[0][2][0][1] = chrisD[2][0][0][1];

    chrisD[2][0][0][2] = 0.5 * (-gtphdth2 * gphph * gtt - gtphdth2 * gtph2 - gtph * gtphd2th * gphph * gtt + gtph3 * gtphd2th - gttdth * gphphdth * gtph2
                                + gphph2 * gttd2th * gtt - gphph * gttd2th * gtph2 + gtt * gphphdth * gtph * gtphdth + 3 * gtph * gtphdth * gphph * gttdth
                                - gphph2 * gttdth2) / fac2;
    chrisD[0][2][0][2] = chrisD[2][0][0][2];


    chrisD[3][3][1][1] = 0.5 * (gphphdr * grrdr - gphphd2r * grr) / grr2;

    chrisD[3][3][1][2] = 0.5 * (gphphdr * grrdth - gphphdrdth * grr) / grr2;

    chrisD[3][3][2][1] = 0.5 * (gphphdth * gththdr - gphphdrdth * gthth) / gthth2;

    chrisD[3][3][2][2] = 0.5 * (gphphdth * gththdth - gphphd2th * gthth) / gthth2;

    chrisD[3][0][1][1] = 0.5 * (gtphdr * grrdr - gtphd2r * grr) / grr2;
    chrisD[0][3][1][1] = chrisD[3][0][1][1];


    chrisD[3][0][1][2] = 0.5 * (gtphdr * grrdth - gtphdrdth * grr) / grr2;
    chrisD[0][3][1][2] = chrisD[3][0][1][2];

    chrisD[3][0][2][1] = 0.5 * (gtphdth * gththdr - gtphdrdth * gthth) / gthth2;
    chrisD[0][3][2][1] = chrisD[3][0][2][1];

    chrisD[3][0][2][2] = 0.5 * (gtphdth * gththdth - gtphd2th * gthth) / gthth2;
    chrisD[0][3][2][2] = chrisD[3][0][2][2];

    chrisD[0][0][1][1] = 0.5 * (gttdr * grrdr - gttd2r * grr) / grr2;

    chrisD[0][0][1][2] = 0.5 * (gttdr * grrdth - gttdrdth * grr) / grr2;

    chrisD[0][0][2][1] = 0.5 * (gttdth * gththdr - gttdrdth * gthth) / gthth2;

    chrisD[0][0][2][2] = 0.5 * (gttdth * gththdth - gttd2th * gthth) / gthth2;
    /*

        std::cout << "ChrisD:\n";
        std::cout << "chrisD[1][1][1][1]\t " << chrisD[1][1][1][1] << "\n";
        std::cout << "chrisD[1][1][1][2]\t " << chrisD[1][1][1][2] << "\n";
        std::cout << "chrisD[1][1][2][1]\t " << chrisD[1][1][2][1] << "\n";
        std::cout << "chrisD[1][1][2][2]\t " << chrisD[1][1][2][2] << "\n";
        std::cout << "chrisD[1][2][1][1]\t " << chrisD[1][2][1][1] << "\n";
        std::cout << "chrisD[1][2][1][2]\t " << chrisD[1][2][1][2] << "\n";
        std::cout << "chrisD[1][2][2][1]\t " << chrisD[1][2][2][1] << "\n";
        std::cout << "chrisD[1][2][2][2]\t " << chrisD[1][2][2][2] << "\n";
        std::cout << "chrisD[1][3][3][1]\t " << chrisD[1][3][3][1] << "\n";
        std::cout << "chrisD[1][3][3][2]\t " << chrisD[1][3][3][2] << "\n";
        std::cout << "chrisD[1][3][0][1]\t " << chrisD[1][3][0][1] << "\n";
        std::cout << "chrisD[1][3][0][2]\t " << chrisD[1][3][0][2] << "\n";
        std::cout << "chrisD[1][0][3][1]\t " << chrisD[1][0][3][1] << "\n";
        std::cout << "chrisD[1][0][3][2]\t " << chrisD[1][0][3][2] << "\n";
        std::cout << "chrisD[1][0][0][1]\t " << chrisD[1][0][0][1] << "\n";
        std::cout << "chrisD[1][0][0][2]\t " << chrisD[1][0][0][2] << "\n";
        std::cout << "chrisD[2][2][1][1]\t " << chrisD[2][2][1][1] << "\n";
        std::cout << "chrisD[2][2][1][2]\t " << chrisD[2][2][1][2] << "\n";
        std::cout << "chrisD[2][2][2][1]\t " << chrisD[2][2][2][1] << "\n";
        std::cout << "chrisD[2][2][2][2]\t " << chrisD[2][2][2][2] << "\n";
        std::cout << "chrisD[2][3][3][1]\t " << chrisD[2][3][3][1] << "\n";
        std::cout << "chrisD[2][3][3][2]\t " << chrisD[2][3][3][2] << "\n";
        std::cout << "chrisD[2][3][0][1]\t " << chrisD[2][3][0][1] << "\n";
        std::cout << "chrisD[2][3][0][2]\t " << chrisD[2][3][0][2] << "\n";
        std::cout << "chrisD[2][0][3][1]\t " << chrisD[2][0][3][1] << "\n";
        std::cout << "chrisD[2][0][3][2]\t " << chrisD[2][0][3][2] << "\n";
        std::cout << "chrisD[2][0][0][1]\t " << chrisD[2][0][0][1] << "\n";
        std::cout << "chrisD[2][0][0][2]\t " << chrisD[2][0][0][2] << "\n";
        std::cout << "chrisD[3][3][1][1]\t " << chrisD[3][3][1][1] << "\n";
        std::cout << "chrisD[3][3][1][2]\t " << chrisD[3][3][1][2] << "\n";
        std::cout << "chrisD[3][3][2][1]\t " << chrisD[3][3][2][1] << "\n";
        std::cout << "chrisD[3][3][2][2]\t " << chrisD[3][3][2][2] << "\n";
        std::cout << "chrisD[3][0][1][1]\t " << chrisD[3][0][1][1] << "\n";
        std::cout << "chrisD[3][0][1][2]\t " << chrisD[3][0][1][2] << "\n";
        std::cout << "chrisD[3][0][2][1]\t " << chrisD[3][0][2][1] << "\n";
        std::cout << "chrisD[3][0][2][2]\t " << chrisD[3][0][2][2] << "\n";
        std::cout << "chrisD[0][0][1][1]\t " << chrisD[0][0][1][1] << "\n";
        std::cout << "chrisD[0][0][1][2]\t " << chrisD[0][0][1][2] << "\n";
        std::cout << "chrisD[0][0][2][1]\t " << chrisD[0][0][2][1] << "\n";
        std::cout << "chrisD[0][0][2][2]\t " << chrisD[0][0][2][2] << "\n";
    */
    return true;
}


bool MetricGlampedakis::calculateRiemann(const double* pos) {
    calcgComps(pos);
    calcgCompsDiff(pos);
    calcgCompsDiff2(pos);

    double gtt2      = gtt * gtt;
    double grr2      = grr * grr;
    double gthth2    = gthth * gthth;
    double gtph2     = gtph * gtph;
    double gtph3     = gtph * gtph * gtph;
    double gphph2    = gphph * gphph;

    double gttdr2    = gttdr * gttdr;
    double gttdth2   = gttdth * gttdth;
    double gtphdth2  = gtphdth * gtphdth;
    double grrdth2   = grrdth * grrdth;
    double gphphdr2  = gphphdr * gphphdr;
    double gththdr2  = gththdr * gththdr;
    double gphphdth2 = gphphdth * gphphdth;
    double gtphdr2   = gtphdr * gtphdr;


    double denom1    = (gphph*gtt-gtph2);
    double denom12   = denom1 * denom1;

    int ti     = 0;
    int ri     = 1;
    int thetai = 2;
    int phii   = 3;

    riem[ri][thetai][ri][thetai] = 0.25*(-2*grrd2th*grr*gthth-2*gththd2r*grr*gthth+gthth*grrdth2+gthth*gththdr*grrdr+grr*gththdr2+grr*grrdth*gththdth)/gthth/grr2;

    riem[ri][thetai][phii][ti] = 0.25/grr*(gtt*gtphdr*gphphdth-gtt*gphphdr*gtphdth+gtph*gphphdr*gttdth-gtph*gttdr*gphphdth+gphph*gttdr*gtphdth-gphph*gtphdr*gttdth)/denom1;

    riem[ri][phii][ri][phii] = 0.25*(-2*gphphd2r*grr*gthth*gphph*gtt+2*gphphd2r*grr*gthth*gtph2+gphphdr*grrdr*gthth*gphph*gtt-gphphdr*grrdr*gthth*gtph2-grrdth*gphphdth*grr*gphph*gtt+grrdth*gphphdth*grr*gtph2+gtt*gphphdr2*grr*gthth-2*gtph*gphphdr*gtphdr*grr*gthth+gphph*gtphdr2*grr*gthth)/gthth/denom1/grr2;

    riem[ri][phii][ri][ti] = -0.25*(2*gtphd2r*grr*gthth*gphph*gtt-2*gtphd2r*grr*gthth*gtph2-gtphdr*grrdr*gthth*gphph*gtt+gtphdr*grrdr*gthth*gtph2+grrdth*gtphdth*grr*gphph*gtt-grrdth*gtphdth*grr*gtph2-gtt*gphphdr*gtphdr*grr*gthth+gtph*gtphdr2*grr*gthth+gtph*gttdr*gphphdr*grr*gthth-gphph*gttdr*gtphdr*grr*gthth)/gthth/denom1/grr2;

    riem[ri][phii][thetai][phii] = 0.25*(-2*gphphdrdth*grr*gthth*gphph*gtt+2*gphphdrdth*grr*gthth*gtph2+gphphdr*grrdth*gthth*gphph*gtt-gphphdr*grrdth*gthth*gtph2+gphphdth*gththdr*grr*gphph*gtt-gphphdth*gththdr*grr*gtph2+gtt*gphphdr*gphphdth*grr*gthth-gtph*gphphdr*gtphdth*grr*gthth-gtph*gtphdr*gphphdth*grr*gthth+gphph*gtphdr*gtphdth*grr*gthth)/gthth/denom1/grr2;

    riem[ri][phii][thetai][ti] = 0.25*(-2*gtphdrdth*grr*gthth*gphph*gtt+2*gtphdrdth*grr*gthth*gtph2+gtphdr*grrdth*gthth*gphph*gtt-gtphdr*grrdth*gthth*gtph2+gtphdth*gththdr*grr*gphph*gtt-gtphdth*gththdr*grr*gtph2+gtt*gtphdr*gphphdth*grr*gthth-gtph*gtphdr*gtphdth*grr*gthth-gtph*gttdr*gphphdth*grr*gthth+gphph*gttdr*gtphdth*grr*gthth)/gthth/denom1/grr2;

    riem[ri][ti][ri][phii] = -0.25*(2*gtphd2r*grr*gthth*gphph*gtt-2*gtphd2r*grr*gthth*gtph2-gtphdr*grrdr*gthth*gphph*gtt+gtphdr*grrdr*gthth*gtph2+grrdth*gtphdth*grr*gphph*gtt-grrdth*gtphdth*grr*gtph2-gtt*gphphdr*gtphdr*grr*gthth+gtph*gtphdr2*grr*gthth+gtph*gttdr*gphphdr*grr*gthth-gphph*gttdr*gtphdr*grr*gthth)/gthth/denom1/grr2;

    riem[ri][ti][ri][ti] = 0.25*(-2*gttd2r*grr*gthth*gphph*gtt+2*gttd2r*grr*gthth*gtph2+gttdr*grrdr*gthth*gphph*gtt-gttdr*grrdr*gthth*gtph2-grrdth*gttdth*grr*gphph*gtt+grrdth*gttdth*grr*gtph2+gtt*gtphdr2*grr*gthth-2*gtph*gttdr*gtphdr*grr*gthth+gphph*gttdr2*grr*gthth)/gthth/denom1/grr2;

    riem[ri][ti][thetai][phii] = 0.25*(-2*gtphdrdth*grr*gthth*gphph*gtt+2*gtphdrdth*grr*gthth*gtph2+gtphdr*grrdth*gthth*gphph*gtt-gtphdr*grrdth*gthth*gtph2+gtphdth*gththdr*grr*gphph*gtt-gtphdth*gththdr*grr*gtph2+gtt*gphphdr*gtphdth*grr*gthth-gtph*gphphdr*gttdth*grr*gthth-gtph*gtphdr*gtphdth*grr*gthth+gphph*gtphdr*gttdth*grr*gthth)/gthth/denom1/grr2;

    riem[ri][ti][thetai][ti] = -0.25*(2*gttdrdth*grr*gthth*gphph*gtt-2*gttdrdth*grr*gthth*gtph2-gttdr*grrdth*gthth*gphph*gtt+gttdr*grrdth*gthth*gtph2-gttdth*gththdr*grr*gphph*gtt+gttdth*gththdr*grr*gtph2-gtt*gtphdr*gtphdth*grr*gthth+gtph*gtphdr*gttdth*grr*gthth+gtph*gttdr*gtphdth*grr*gthth-gphph*gttdr*gttdth*grr*gthth)/gthth/denom1/grr2;

    riem[thetai][ri][ri][thetai] = 0.25/grr/gthth2*(2*grrd2th*grr*gthth+2*gththd2r*grr*gthth-gthth*grrdth2-gthth*gththdr*grrdr-grr*gththdr2-grr*grrdth*gththdth);

    riem[thetai][ri][phii][ti] = 0.25*(-gtt*gtphdr*gphphdth+gtt*gphphdr*gtphdth-gtph*gphphdr*gttdth+gtph*gttdr*gphphdth-gphph*gttdr*gtphdth+gphph*gtphdr*gttdth)/denom1/gthth;

    riem[thetai][phii][ri][phii] = 0.25*(-2*gphphdrdth*grr*gthth*gphph*gtt+2*gphphdrdth*grr*gthth*gtph2+gphphdr*grrdth*gthth*gphph*gtt-gphphdr*grrdth*gthth*gtph2+gphphdth*gththdr*grr*gphph*gtt-gphphdth*gththdr*grr*gtph2+gtt*gphphdr*gphphdth*grr*gthth-gtph*gphphdr*gtphdth*grr*gthth-gtph*gtphdr*gphphdth*grr*gthth+gphph*gtphdr*gtphdth*grr*gthth)/grr/denom1/gthth2;

    riem[thetai][phii][ri][ti] = 0.25*(-2*gtphdrdth*grr*gthth*gphph*gtt+2*gtphdrdth*grr*gthth*gtph2+gtphdr*grrdth*gthth*gphph*gtt-gtphdr*grrdth*gthth*gtph2+gtphdth*gththdr*grr*gphph*gtt-gtphdth*gththdr*grr*gtph2+gtt*gphphdr*gtphdth*grr*gthth-gtph*gphphdr*gttdth*grr*gthth-gtph*gtphdr*gtphdth*grr*gthth+gphph*gtphdr*gttdth*grr*gthth)/grr/denom1/gthth2;

    riem[thetai][phii][thetai][phii] = 0.25*(-2*gphphd2th*grr*gthth*gphph*gtt+2*gphphd2th*grr*gthth*gtph2-gphphdr*gththdr*gthth*gphph*gtt+gphphdr*gththdr*gthth*gtph2+gphphdth*gththdth*grr*gphph*gtt-gphphdth*gththdth*grr*gtph2+gtt*gphphdth2*grr*gthth-2*gtph*gphphdth*gtphdth*grr*gthth+gphph*gtphdth2*grr*gthth)/grr/denom1/gthth2;

    riem[thetai][phii][thetai][ti] =  -0.25*(2*gtphd2th*grr*gthth*gphph*gtt-2*gtphd2th*grr*gthth*gtph2+gtphdr*gththdr*gthth*gphph*gtt-gtphdr*gththdr*gthth*gtph2-gtphdth*gththdth*grr*gphph*gtt+gtphdth*gththdth*grr*gtph2-gtt*gphphdth*gtphdth*grr*gthth+gtph*gtphdth2*grr*gthth+gtph*gttdth*gphphdth*grr*gthth-gphph*gttdth*gtphdth*grr*gthth)/grr/denom1/gthth2;

    riem[thetai][ti][ri][phii]= 0.25*(-2*gtphdrdth*grr*gthth*gphph*gtt+2*gtphdrdth*grr*gthth*gtph2+gtphdr*grrdth*gthth*gphph*gtt-gtphdr*grrdth*gthth*gtph2+gtphdth*gththdr*grr*gphph*gtt-gtphdth*gththdr*grr*gtph2+gtt*gtphdr*gphphdth*grr*gthth-gtph*gtphdr*gtphdth*grr*gthth-gtph*gttdr*gphphdth*grr*gthth+gphph*gttdr*gtphdth*grr*gthth)/grr/denom1/gthth2;

    riem[thetai][ti][ri][ti] = -0.25*(2*gttdrdth*grr*gthth*gphph*gtt-2*gttdrdth*grr*gthth*gtph2-gttdr*grrdth*gthth*gphph*gtt+gttdr*grrdth*gthth*gtph2-gttdth*gththdr*grr*gphph*gtt+gttdth*gththdr*grr*gtph2-gtt*gtphdr*gtphdth*grr*gthth+gtph*gtphdr*gttdth*grr*gthth+gtph*gttdr*gtphdth*grr*gthth-gphph*gttdr*gttdth*grr*gthth)/grr/denom1/gthth2;

    riem[thetai][ti][thetai][phii] = -0.25*(2*gtphd2th*grr*gthth*gphph*gtt-2*gtphd2th*grr*gthth*gtph2+gtphdr*gththdr*gthth*gphph*gtt-gtphdr*gththdr*gthth*gtph2-gtphdth*gththdth*grr*gphph*gtt+gtphdth*gththdth*grr*gtph2-gtt*gphphdth*gtphdth*grr*gthth+gtph*gtphdth2*grr*gthth+gtph*gttdth*gphphdth*grr*gthth-gphph*gttdth*gtphdth*grr*gthth)/grr/denom1/gthth2;

    riem[thetai][ti][thetai][ti] = 0.25*(-2*gttd2th*grr*gthth*gphph*gtt+2*gttd2th*grr*gthth*gtph2-gttdr*gththdr*gthth*gphph*gtt+gttdr*gththdr*gthth*gtph2+gttdth*gththdth*grr*gphph*gtt-gttdth*gththdth*grr*gtph2+gtt*gtphdth2*grr*gthth-2*gtph*gttdth*gtphdth*grr*gthth+gphph*gttdth2*grr*gthth)/grr/denom1/gthth2;

    riem[phii][ri][ri][phii] = -0.25/grr*(-2*gphphd2r*grr*gthth*gphph*gtt2+2*gtt*gphphd2r*grr*gthth*gtph2+gphphdr*grrdr*gthth*gphph*gtt2-gtt*gphphdr*grrdr*gthth*gtph2-grrdth*gphphdth*grr*gphph*gtt2+gtt*grrdth*gphphdth*grr*gtph2+gtt2*gphphdr2*grr*gthth-3*gtt*gtph*gphphdr*gtphdr*grr*gthth+gtt*gphph*gtphdr2*grr*gthth+2*gtph*gtphd2r*grr*gthth*gphph*gtt-2*gtphd2r*grr*gthth*gtph3-gtph*gtphdr*grrdr*gthth*gphph*gtt+gtphdr*grrdr*gthth*gtph3+gtph*grrdth*gtphdth*grr*gphph*gtt-grrdth*gtphdth*grr*gtph3+gtph2*gtphdr2*grr*gthth+gtph2*gttdr*gphphdr*grr*gthth-gtph*gphph*gttdr*gtphdr*grr*gthth)/denom12/gthth;

    riem[phii][ri][ri][ti] = 0.25*(2*gtphd2r*grr*gthth*gphph*gtt2-2*gtt*gtphd2r*grr*gthth*gtph2-gtphdr*grrdr*gthth*gphph*gtt2+gtt*gtphdr*grrdr*gthth*gtph2+grrdth*gtphdth*grr*gphph*gtt2-gtt*grrdth*gtphdth*grr*gtph2-gtt2*gphphdr*gtphdr*grr*gthth+2*gtt*gtph*gtphdr2*grr*gthth+gtt*gtph*gttdr*gphphdr*grr*gthth-gtt*gphph*gttdr*gtphdr*grr*gthth-2*gtph*gttd2r*grr*gthth*gphph*gtt+2*gttd2r*grr*gthth*gtph3+gtph*gttdr*grrdr*gthth*gphph*gtt-gttdr*grrdr*gthth*gtph3-gtph*grrdth*gttdth*grr*gphph*gtt+grrdth*gttdth*grr*gtph3-2*gtph2*gttdr*gtphdr*grr*gthth+gtph*gphph*gttdr2*grr*gthth)/grr/denom12/gthth;

    riem[phii][ri][thetai][phii] = 0.25*(2*gphphdrdth*grr*gthth*gphph*gtt2-2*gtt*gphphdrdth*grr*gthth*gtph2-gphphdr*grrdth*gthth*gphph*gtt2+gtt*gphphdr*grrdth*gthth*gtph2-gphphdth*gththdr*grr*gphph*gtt2+gtt*gphphdth*gththdr*grr*gtph2-gtt2*gphphdr*gphphdth*grr*gthth+2*gtt*gtph*gphphdr*gtphdth*grr*gthth+gtt*gtph*gtphdr*gphphdth*grr*gthth-gtt*gphph*gtphdr*gtphdth*grr*gthth-2*gtph*gtphdrdth*grr*gthth*gphph*gtt+2*gtphdrdth*grr*gthth*gtph3+gtph*gtphdr*grrdth*gthth*gphph*gtt-gtphdr*grrdth*gthth*gtph3+gtph*gtphdth*gththdr*grr*gphph*gtt-gtphdth*gththdr*grr*gtph3-gtph2*gphphdr*gttdth*grr*gthth-gtph2*gtphdr*gtphdth*grr*gthth+gtph*gphph*gtphdr*gttdth*grr*gthth)/grr/denom12/gthth;

    riem[phii][ri][thetai][ti] = -0.25/grr*(-2*gtphdrdth*grr*gthth*gphph*gtt2+2*gtt*gtphdrdth*grr*gthth*gtph2+gtphdr*grrdth*gthth*gphph*gtt2-gtt*gtphdr*grrdth*gthth*gtph2+gtphdth*gththdr*grr*gphph*gtt2-gtt*gtphdth*gththdr*grr*gtph2+gtt2*gtphdr*gphphdth*grr*gthth-2*gtt*gtph*gtphdr*gtphdth*grr*gthth-gtt*gtph*gttdr*gphphdth*grr*gthth+gtt*gphph*gttdr*gtphdth*grr*gthth+2*gtph*gttdrdth*grr*gthth*gphph*gtt-2*gttdrdth*grr*gthth*gtph3-gtph*gttdr*grrdth*gthth*gphph*gtt+gttdr*grrdth*gthth*gtph3-gtph*gttdth*gththdr*grr*gphph*gtt+gttdth*gththdr*grr*gtph3+gtph2*gtphdr*gttdth*grr*gthth+gtph2*gttdr*gtphdth*grr*gthth-gtph*gphph*gttdr*gttdth*grr*gthth)/denom12/gthth;

    riem[phii][thetai][ri][phii] = 0.25*(2*gphphdrdth*grr*gthth*gphph*gtt2-2*gtt*gphphdrdth*grr*gthth*gtph2-gphphdr*grrdth*gthth*gphph*gtt2+gtt*gphphdr*grrdth*gthth*gtph2-gphphdth*gththdr*grr*gphph*gtt2+gtt*gphphdth*gththdr*grr*gtph2-gtt2*gphphdr*gphphdth*grr*gthth+gtt*gtph*gphphdr*gtphdth*grr*gthth+2*gtt*gtph*gtphdr*gphphdth*grr*gthth-gtt*gphph*gtphdr*gtphdth*grr*gthth-2*gtph*gtphdrdth*grr*gthth*gphph*gtt+2*gtphdrdth*grr*gthth*gtph3+gtph*gtphdr*grrdth*gthth*gphph*gtt-gtphdr*grrdth*gthth*gtph3+gtph*gtphdth*gththdr*grr*gphph*gtt-gtphdth*gththdr*grr*gtph3-gtph2*gtphdr*gtphdth*grr*gthth-gtph2*gttdr*gphphdth*grr*gthth+gtph*gphph*gttdr*gtphdth*grr*gthth)/grr/denom12/gthth;

    riem[phii][thetai][ri][ti] = -0.25/gthth*(-2*gtphdrdth*grr*gthth*gphph*gtt2+2*gtt*gtphdrdth*grr*gthth*gtph2+gtphdr*grrdth*gthth*gphph*gtt2-gtt*gtphdr*grrdth*gthth*gtph2+gtphdth*gththdr*grr*gphph*gtt2-gtt*gtphdth*gththdr*grr*gtph2+gtt2*gphphdr*gtphdth*grr*gthth-gtt*gtph*gphphdr*gttdth*grr*gthth-2*gtt*gtph*gtphdr*gtphdth*grr*gthth+gtt*gphph*gtphdr*gttdth*grr*gthth+2*gtph*gttdrdth*grr*gthth*gphph*gtt-2*gttdrdth*grr*gthth*gtph3-gtph*gttdr*grrdth*gthth*gphph*gtt+gttdr*grrdth*gthth*gtph3-gtph*gttdth*gththdr*grr*gphph*gtt+gttdth*gththdr*grr*gtph3+gtph2*gtphdr*gttdth*grr*gthth+gtph2*gttdr*gtphdth*grr*gthth-gtph*gphph*gttdr*gttdth*grr*gthth)/denom12/grr;

    riem[phii][thetai][thetai][phii] = -0.25/gthth*(-2*gphphd2th*grr*gthth*gphph*gtt2+2*gtt*gphphd2th*grr*gthth*gtph2-gphphdr*gththdr*gthth*gphph*gtt2+gtt*gphphdr*gththdr*gthth*gtph2+gphphdth*gththdth*grr*gphph*gtt2-gtt*gphphdth*gththdth*grr*gtph2+gtt2*gphphdth2*grr*gthth-3*gtt*gtph*gphphdth*gtphdth*grr*gthth+gtt*gphph*gtphdth2*grr*gthth+2*gtph*gtphd2th*grr*gthth*gphph*gtt-2*gtphd2th*grr*gthth*gtph3+gtph*gtphdr*gththdr*gthth*gphph*gtt-gtphdr*gththdr*gthth*gtph3-gtph*gtphdth*gththdth*grr*gphph*gtt+gtphdth*gththdth*grr*gtph3+gtph2*gtphdth2*grr*gthth+gtph2*gttdth*gphphdth*grr*gthth-gtph*gphph*gttdth*gtphdth*grr*gthth)/denom12/grr;

    riem[phii][thetai][thetai][ti] = 0.25*(2*gtphd2th*grr*gthth*gphph*gtt2-2*gtt*gtphd2th*grr*gthth*gtph2+gtphdr*gththdr*gthth*gphph*gtt2-gtt*gtphdr*gththdr*gthth*gtph2-gtphdth*gththdth*grr*gphph*gtt2+gtt*gtphdth*gththdth*grr*gtph2-gtt2*gphphdth*gtphdth*grr*gthth+2*gtt*gtph*gtphdth2*grr*gthth+gtt*gtph*gttdth*gphphdth*grr*gthth-gtt*gphph*gttdth*gtphdth*grr*gthth-2*gtph*gttd2th*grr*gthth*gphph*gtt+2*gttd2th*grr*gthth*gtph3-gtph*gttdr*gththdr*gthth*gphph*gtt+gttdr*gththdr*gthth*gtph3+gtph*gttdth*gththdth*grr*gphph*gtt-gttdth*gththdth*grr*gtph3-2*gtph2*gttdth*gtphdth*grr*gthth+gtph*gphph*gttdth2*grr*gthth)/grr/denom12/gthth;

    riem[phii][phii][ri][thetai] = -0.25*gtph*(-gtt*gtphdr*gphphdth+gtt*gphphdr*gtphdth-gtph*gphphdr*gttdth+gtph*gttdr*gphphdth-gphph*gttdr*gtphdth+gphph*gtphdr*gttdth)/denom12;

    riem[phii][phii][phii][ti] = 0.25*gtph*(-gthth*gttdr*gphphdr+gthth*gtphdr2-grr*gttdth*gphphdth+grr*gtphdth2)/grr/gthth/denom1;

    riem[phii][ti][ri][thetai] = -0.25*gtt*(-gtt*gtphdr*gphphdth+gtt*gphphdr*gtphdth-gtph*gphphdr*gttdth+gtph*gttdr*gphphdth-gphph*gttdr*gtphdth+gphph*gtphdr*gttdth)/denom12;

    riem[phii][ti][phii][ti] = 0.25*gtt*(-gthth*gttdr*gphphdr+gthth*gtphdr2-grr*gttdth*gphphdth+grr*gtphdth2)/grr/gthth/denom1;

    riem[ti][ri][ri][phii] = 0.25/grr*(-2*gtph*gphphd2r*grr*gthth*gphph*gtt+2*gphphd2r*grr*gthth*gtph3+gtph*gphphdr*grrdr*gthth*gphph*gtt-gphphdr*grrdr*gthth*gtph3-gtph*grrdth*gphphdth*grr*gphph*gtt+grrdth*gphphdth*grr*gtph3+gtph*gtt*gphphdr2*grr*gthth-2*gtph2*gphphdr*gtphdr*grr*gthth+2*gtph*gphph*gtphdr2*grr*gthth+2*gtphd2r*grr*gthth*gphph2*gtt-2*gphph*gtphd2r*grr*gthth*gtph2-gtphdr*grrdr*gthth*gphph2*gtt+gphph*gtphdr*grrdr*gthth*gtph2+grrdth*gtphdth*grr*gphph2*gtt-gphph*grrdth*gtphdth*grr*gtph2-gphph*gtt*gphphdr*gtphdr*grr*gthth+gphph*gtph*gttdr*gphphdr*grr*gthth-gphph2*gttdr*gtphdr*grr*gthth)/denom12/gthth;

    riem[ti][ri][ri][ti] = 0.25/grr*(-2*gtph*gtphd2r*grr*gthth*gphph*gtt+2*gtphd2r*grr*gthth*gtph3+gtph*gtphdr*grrdr*gthth*gphph*gtt-gtphdr*grrdr*gthth*gtph3-gtph*grrdth*gtphdth*grr*gphph*gtt+grrdth*gtphdth*grr*gtph3+gtt*gtph*gphphdr*gtphdr*grr*gthth-gtph2*gtphdr2*grr*gthth-gtph2*gttdr*gphphdr*grr*gthth+3*gtph*gphph*gttdr*gtphdr*grr*gthth+2*gttd2r*grr*gthth*gphph2*gtt-2*gphph*gttd2r*grr*gthth*gtph2-gttdr*grrdr*gthth*gphph2*gtt+gphph*gttdr*grrdr*gthth*gtph2+grrdth*gttdth*grr*gphph2*gtt-gphph*grrdth*gttdth*grr*gtph2-gtt*gphph*gtphdr2*grr*gthth-gphph2*gttdr2*grr*gthth)/denom12/gthth;

    riem[ti][ri][thetai][phii] = 0.25/grr*(-2*gtph*gphphdrdth*grr*gthth*gphph*gtt+2*gphphdrdth*grr*gthth*gtph3+gtph*gphphdr*grrdth*gthth*gphph*gtt-gphphdr*grrdth*gthth*gtph3+gtph*gphphdth*gththdr*grr*gphph*gtt-gphphdth*gththdr*grr*gtph3+gtph*gtt*gphphdr*gphphdth*grr*gthth-gtph2*gphphdr*gtphdth*grr*gthth-gtph2*gtphdr*gphphdth*grr*gthth+2*gtph*gphph*gtphdr*gtphdth*grr*gthth+2*gtphdrdth*grr*gthth*gphph2*gtt-2*gphph*gtphdrdth*grr*gthth*gtph2-gtphdr*grrdth*gthth*gphph2*gtt+gphph*gtphdr*grrdth*gthth*gtph2-gtphdth*gththdr*grr*gphph2*gtt+gphph*gtphdth*gththdr*grr*gtph2-gphph*gtt*gphphdr*gtphdth*grr*gthth+gphph*gtph*gphphdr*gttdth*grr*gthth-gphph2*gtphdr*gttdth*grr*gthth)/denom12/gthth;

    riem[ti][ri][thetai][ti] = 0.25/grr*(-2*gtph*gtphdrdth*grr*gthth*gphph*gtt+2*gtphdrdth*grr*gthth*gtph3+gtph*gtphdr*grrdth*gthth*gphph*gtt-gtphdr*grrdth*gthth*gtph3+gtph*gtphdth*gththdr*grr*gphph*gtt-gtphdth*gththdr*grr*gtph3+gtt*gtph*gtphdr*gphphdth*grr*gthth-gtph2*gtphdr*gtphdth*grr*gthth-gtph2*gttdr*gphphdth*grr*gthth+2*gtph*gphph*gttdr*gtphdth*grr*gthth+2*gttdrdth*grr*gthth*gphph2*gtt-2*gphph*gttdrdth*grr*gthth*gtph2-gttdr*grrdth*gthth*gphph2*gtt+gphph*gttdr*grrdth*gthth*gtph2-gttdth*gththdr*grr*gphph2*gtt+gphph*gttdth*gththdr*grr*gtph2-gtt*gphph*gtphdr*gtphdth*grr*gthth+gtph*gphph*gtphdr*gttdth*grr*gthth-gphph2*gttdr*gttdth*grr*gthth)/denom12/gthth;

    riem[ti][thetai][ri][phii] = 0.25/gthth*(-2*gtph*gphphdrdth*grr*gthth*gphph*gtt+2*gphphdrdth*grr*gthth*gtph3+gtph*gphphdr*grrdth*gthth*gphph*gtt-gphphdr*grrdth*gthth*gtph3+gtph*gphphdth*gththdr*grr*gphph*gtt-gphphdth*gththdr*grr*gtph3+gtph*gtt*gphphdr*gphphdth*grr*gthth-gtph2*gphphdr*gtphdth*grr*gthth-gtph2*gtphdr*gphphdth*grr*gthth+2*gtph*gphph*gtphdr*gtphdth*grr*gthth+2*gtphdrdth*grr*gthth*gphph2*gtt-2*gphph*gtphdrdth*grr*gthth*gtph2-gtphdr*grrdth*gthth*gphph2*gtt+gphph*gtphdr*grrdth*gthth*gtph2-gtphdth*gththdr*grr*gphph2*gtt+gphph*gtphdth*gththdr*grr*gtph2-gphph*gtt*gtphdr*gphphdth*grr*gthth+gphph*gtph*gttdr*gphphdth*grr*gthth-gphph2*gttdr*gtphdth*grr*gthth)/denom12/grr;

    riem[ti][thetai][ri][ti] = 0.25/gthth*(-2*gtph*gtphdrdth*grr*gthth*gphph*gtt+2*gtphdrdth*grr*gthth*gtph3+gtph*gtphdr*grrdth*gthth*gphph*gtt-gtphdr*grrdth*gthth*gtph3+gtph*gtphdth*gththdr*grr*gphph*gtt-gtphdth*gththdr*grr*gtph3+gtt*gtph*gphphdr*gtphdth*grr*gthth-gtph2*gphphdr*gttdth*grr*gthth-gtph2*gtphdr*gtphdth*grr*gthth+2*gtph*gphph*gtphdr*gttdth*grr*gthth+2*gttdrdth*grr*gthth*gphph2*gtt-2*gphph*gttdrdth*grr*gthth*gtph2-gttdr*grrdth*gthth*gphph2*gtt+gphph*gttdr*grrdth*gthth*gtph2-gttdth*gththdr*grr*gphph2*gtt+gphph*gttdth*gththdr*grr*gtph2-gtt*gphph*gtphdr*gtphdth*grr*gthth+gtph*gphph*gttdr*gtphdth*grr*gthth-gphph2*gttdr*gttdth*grr*gthth)/denom12/grr;

    riem[ti][thetai][thetai][phii] = 0.25/gthth*(-2*gtph*gphphd2th*grr*gthth*gphph*gtt+2*gphphd2th*grr*gthth*gtph3-gtph*gphphdr*gththdr*gthth*gphph*gtt+gphphdr*gththdr*gthth*gtph3+gtph*gphphdth*gththdth*grr*gphph*gtt-gphphdth*gththdth*grr*gtph3+gtph*gtt*gphphdth2*grr*gthth-2*gtph2*gphphdth*gtphdth*grr*gthth+2*gtph*gphph*gtphdth2*grr*gthth+2*gtphd2th*grr*gthth*gphph2*gtt-2*gphph*gtphd2th*grr*gthth*gtph2+gtphdr*gththdr*gthth*gphph2*gtt-gphph*gtphdr*gththdr*gthth*gtph2-gtphdth*gththdth*grr*gphph2*gtt+gphph*gtphdth*gththdth*grr*gtph2-gphph*gtt*gphphdth*gtphdth*grr*gthth+gphph*gtph*gttdth*gphphdth*grr*gthth-gphph2*gttdth*gtphdth*grr*gthth)/denom12/grr;

    riem[ti][thetai][thetai][ti] = 0.25/gthth*(-2*gtph*gtphd2th*grr*gthth*gphph*gtt+2*gtphd2th*grr*gthth*gtph3-gtph*gtphdr*gththdr*gthth*gphph*gtt+gtphdr*gththdr*gthth*gtph3+gtph*gtphdth*gththdth*grr*gphph*gtt-gtphdth*gththdth*grr*gtph3+gtt*gtph*gphphdth*gtphdth*grr*gthth-gtph2*gtphdth2*grr*gthth-gtph2*gttdth*gphphdth*grr*gthth+3*gtph*gphph*gttdth*gtphdth*grr*gthth+2*gttd2th*grr*gthth*gphph2*gtt-2*gphph*gttd2th*grr*gthth*gtph2+gttdr*gththdr*gthth*gphph2*gtt-gphph*gttdr*gththdr*gthth*gtph2-gttdth*gththdth*grr*gphph2*gtt+gphph*gttdth*gththdth*grr*gtph2-gtt*gphph*gtphdth2*grr*gthth-gphph2*gttdth2*grr*gthth)/denom12/grr;

    riem[ti][phii][ri][thetai] = 0.25*gphph*(-gtt*gtphdr*gphphdth+gtt*gphphdr*gtphdth-gtph*gphphdr*gttdth+gtph*gttdr*gphphdth-gphph*gttdr*gtphdth+gphph*gtphdr*gttdth)/denom12;

    riem[ti][phii][phii][ti] = -0.25*gphph*(-gthth*gttdr*gphphdr+gthth*gtphdr2-grr*gttdth*gphphdth+grr*gtphdth2)/grr/gthth/denom1;

    riem[ti][ti][ri][thetai] = 0.25*gtph*(-gtt*gtphdr*gphphdth+gtt*gphphdr*gtphdth-gtph*gphphdr*gttdth+gtph*gttdr*gphphdth-gphph*gttdr*gtphdth+gphph*gtphdr*gttdth)/denom12;

    riem[ti][ti][phii][ti] = -0.25*gtph*(-gthth*gttdr*gphphdr+gthth*gtphdr2-grr*gttdth*gphphdth+grr*gtphdth2)/grr/gthth/denom1;

    return true;

}


/*! Transform local 4-direction to coordinate 4-direction.
 *
 *  \param  pos  :  pointer to position array.
 *  \param  ldir :  pointer to local direction array.
 *  \param  dir  :  pointer to calculated coordinate direction array.
 *  \param  type :  type of tetrad.
 */
void MetricGlampedakis::localToCoord(const double* pos, const double* ldir, double* dir,
                                     enum_nat_tetrad_type) {
    calcgComps(pos);
    double omega = -gtph / gphph;
    double zeta = 0.0;

    double gam = 1.0 / sqrt(-(gtt + 2.0 * zeta * gtph + zeta * zeta * gphph));
    double dlt = 1.0 / sqrt(gtph * gtph - gtt * gphph);
    double w1 = gtph + zeta * gphph;
    double w2 = gtt + zeta * gtph;

    dir[0] = gam * (ldir[0] + dlt * w1 * ldir[3]);
    dir[1] = ldir[1] / sqrt(grr);
    dir[2] = ldir[2] / sqrt(gthth);
    dir[3] = gam * (ldir[0] * zeta - dlt * w2 * ldir[3]);

}

/*! Transform coordinate 4-direction to local 4-direction.
 *
 *  \param  pos  :  pointer to position array.
 *  \param  cdir :  pointer to coordinate direction.
 *  \param  ldir :  pointer to calculated local direction array.
 *  \param  type :  type of tetrad.
 */
void MetricGlampedakis::coordToLocal(const double* , const double* , double* ,
                                     enum_nat_tetrad_type) {
    // TODO
}


/*! Test break condition.
 *
 *  \param pos    : pointer to position array.
 *  \return true  : radial position r < 0.0 or  r^2<=(1.0+eps)*rs^2.
 *  \return false : position is valid.
 */
bool MetricGlampedakis::breakCondition(const double* pos) {
    bool br = false;
    if (pos[1] <= 2.0 * mMass || std::isnan(pos[0]) || std::isnan(pos[1]) || std::isnan(pos[2]) || std::isnan(pos[3])) {
        br = true;
    } else {
        calcgComps(pos);

        double dd = gtph * gtph - gtt * gphph;


        if (grr <= 0.0 || gthth <= 0.0 || dd <= 0.0) {
            br = true;
        }
    }
    return br;
}




/*! Set parameter 'pName' to 'val'.
 *
 *  Set 'mass' or 'lambda' parameter.
 */
bool MetricGlampedakis::setParam(std::string pName, double val) {
    Metric::setParam(pName, val);
    if (pName == "mass") {
        mMass = val;
        M = mMass;
        rs = 2.0 * mMass;
    } else if (pName == "angmom") {
        mAngmom = val;
        a   = mAngmom;
        a2  = a * a;
    } else if (pName == "epsilon") {
        mEpsilon = val;
    }
    return true;
}

/*! Generate report.
 */
bool MetricGlampedakis::report(const vec4 , const vec4 , std::string &text) {
    std::stringstream ss;
    ss << "Report for HartleThorne metric\n\tcoordinate : (t,x,y,z)\n";
    ss << "---------------------------------------------------------------\n";
    ss << "  physical units ..................... no\n";
    ss.precision(DEF_FIXED_REPORT_PRECISION);
    ss.setf(std::ios::fixed);
    ss << "  param .............................. mass   = " << mMass << std::endl;
    ss << "  param .............................. angmom = " << mAngmom << std::endl;
    ss << "  param .............................. eta    = " << mEpsilon << std::endl;

    text = ss.str();
    return true;
}

// ********************************* protected methods *****************************

// ********************************* Kerr routines *********************************
void MetricGlampedakis::calcKerr(const double* pos) {
    kerrcalled++;
    if ((oldkerrpos[0] == pos[1]) && (oldkerrpos[1] == pos[2])) {
        kerrskip++;
        return;
    }
    oldkerrpos[0] = pos[1];
    oldkerrpos[1] = pos[2];


    r = pos[1];
    r2 = r * r;

    theta = pos[2];
    sth = sin(theta);
    cth = cos(theta);

    sth2 = sth * sth;
    cth2 = cth * cth;

    sigma   = r2 + a2 * cth2;

    delta   = r2 - rs * r + a2;

    ktt     = -1.0 + rs * r / sigma;
    ktph    = -rs * a * r * sth2 / sigma;
    krr     = sigma / delta;
    kthth   = sigma;
    kphph   = (r2 + a2 + rs * a2 * r * sth2 / sigma) * sth2;

}

void MetricGlampedakis::calcKerrDiff(const double* pos) {
    kerrdiffcalled++;
    if ((oldkerrdiffpos[0] == pos[1]) && (oldkerrdiffpos[1] == pos[2])) {
        kerrdiffskip++;
        return;
    }
    oldkerrdiffpos[0] = pos[1];
    oldkerrdiffpos[1] = pos[2];

    calcKerr(pos);

    sigmadr     = 2.0 * r;
    sigmadth    = -2.0 * a2 * cth * sth;

    deltadr     = 2.0 * r - rs;

    sigma2      = sigma * sigma;
    delta2      = delta * delta;

    sigma2Inv   = 1.0 / sigma2;

    kttdr       = rs * (sigma - r * sigmadr) * sigma2Inv;
    kttdth      = -rs * r * sigmadth * sigma2Inv;

    ktphdr      = -rs * a * sth2 * (sigma - r * sigmadr) * sigma2Inv;
    ktphdth     = -rs * a * r * sth * (2.0 * cth * sigma - sth * sigmadth) * sigma2Inv;

    krrdr       = (sigmadr * delta - sigma * deltadr) / delta2;
    krrdth      = sigmadth / delta;

    kththdr     = sigmadr;
    kththdth    = sigmadth;

    kphphdr     = (2.0 * r * sigma2 + rs * a2 * sth2 * (sigma - r * sigmadr)) * sth2 * sigma2Inv;
    kphphdth    = sth * (rs * a2 * r * sth2 * (4.0 * cth * sigma - sth * sigmadth) + 2.0 * cth * sigma2 * (r2 + a2)) * sigma2Inv;

}

void MetricGlampedakis::calcKerrDiff2(const double* pos) {
    kerrdiff2called++;
    if ((oldkerrdiff2pos[0] == pos[1]) && (oldkerrdiff2pos[1] == pos[2])) {
        kerrdiff2skip++;
        return;
    }
    oldkerrdiff2pos[0] = pos[1];
    oldkerrdiff2pos[1] = pos[2];

    calcKerrDiff(pos);


    sigmad2r     = 2.0;
    sigmad2th    = 2.0 * a2 * (sth2 - cth2);

    deltad2r     = 2.0;
    sigma3  = sigma2 * sigma;
    delta3  = delta2 * delta;

    sth3 = sth2 * sth;
    cth3 = cth2 * cth;
    cth4 = cth3 * cth;

    double sigmadr2     = sigmadr * sigmadr;
    double sigmadth2    = sigmadth * sigmadth;
    double deltadr2     = deltadr * deltadr;

    kttd2r    = -rs * (2.0 * sigmadr * sigma - 2.0 * r * sigmadr2 + r * sigmad2r * sigma) / sigma3;
    kttd2th   = -rs * r * (-2.0 * sigmadth2 + sigmad2th * sigma) / sigma3;
    kttdrdth  = -rs * (sigmadth * sigma - 2.0 * r * sigmadr  * sigmadth) / sigma3;

    ktphd2r    = rs * a * sth2 * (2.0 * sigmadr * sigma - 2.0 * r * sigmadr2 + r * sigmad2r * sigma) / sigma3;
    ktphd2th   = rs * a * r * (sigma2 * (2.0 - 4.0 * cth2) + sigma * sigmadth * 4.0 * sth * cth
                               + sth2 * (sigma * sigmad2th - 2.0 * sigmadth2)) / sigma3;

    ktphdrdth  = rs * a * sth * (-2.0 * cth * sigma2 + sth * sigmadth * sigma + 2.0 * r * sigmadr * cth * sigma
                                 - 2.0 * r * sth * sigmadr * sigmadth) / sigma3;

    krrd2r    = (sigmad2r * delta2 - 2.0 * sigmadr * deltadr * delta + sigma * (2.0 * deltadr2 - deltad2r * delta)) / delta3;
    krrd2th   = sigmad2th / delta;
    krrdrdth  = -sigmadth * deltadr / delta2;

    kththd2r    = sigmad2r;
    kththd2th   = sigmad2th;
    kththdrdth  = 0.0;

    kphphd2r   = -sth2 * (-2.0 * sigma3 + rs * a2 * (sth2 * (2.0 * sigmadr * sigma + r * (sigmad2r * sigma - 2.0 * sigmadr2))))
                 / sigma3;
    kphphd2th   = (r2 + a2) * (4.0 * cth2 - 2.0) + rs * a2 * r * ((1.0 - 2.0 * cth2 + cth4) * (2.0 * sigmadth2 - sigmad2th * sigma)
                  + sigma2 * (20.0 * cth2 - 16.0 * cth4 - 4.0)
                  - sigmadth * sigma * 8.0 * sth3 * cth) / sigma3;

    kphphdrdth  = 4.0 * r * sth * cth + rs * a2 * sth * (4.0 * cth * sth2 * (sigma2 - r * sigma * sigmadr)
                  + sigmadth * sth3 * (2.0 * r * sigmadr - sigma)) / sigma3;
    /*
        std::cout.precision(12);

        std::cout << "ktt: " << kttd2r << " " << kttd2th << " " << kttdrdth << "\n";
        std::cout << "ktph: " << ktphd2r << " " << ktphd2th << " " << ktphdrdth << "\n";
        std::cout << "krr: " << krrd2r << " " << krrd2th << " " << krrdrdth << "\n";
        std::cout << "kthth: " << kththd2r << " " << kththd2th << " " << kththdrdth << "\n";
        std::cout << "kphph: " << kphphd2r << " " << kphphd2th << " " << kphphdrdth << "\n";
    */
}

// ********************************* Glampedakis routines **************************
void MetricGlampedakis::calcGlampedakis(const double* pos) {
    glampcalled++;
    if ((oldglamppos[0] == pos[1]) && (oldglamppos[1] == pos[2])) {
        glampskip++;
        return;
    }
    oldglamppos[0] = pos[1];
    oldglamppos[1] = pos[2];
    M2   = M * M;
    r    = pos[1];
    r2   = r * r;

    Mr = r * M;

    theta = pos[2];
    sth  = sin(theta);
    cth  = cos(theta);
    sth2 = sth * sth;
    cth2 = cth * cth;

    lgr = log(r / (r - 2.0 * M));

    F1 = -5.0 * (r - M) / (8.0 * Mr * (r - 2.0 * M)) * (2.0 * M2 + 6.0 * Mr - 3.0 * r2)
         - 15.0 * r * (r - 2 * M) / (16.0 * M2) * lgr;
    F2 = 5.0 / (8.0 * Mr) * (2.0 * M2 - 3.0 * Mr - 3.0 * r2) + 15.0 / (16.0 * M2) * (r2 - 2.0 * M2) * lgr;

    cf = (1.0 - 3.0 * cth2);

    ff1 = cf * F1;
    ff2 = cf * F2;

    fSchw = (1.0 - rs / r);

    htt     = ff1 / fSchw;
    hrr     = ff1 * fSchw;
    hthth   = -ff2 / r2;
    hphph   = -ff2 / r2 / sth2;

    ktt2     = ktt * ktt;
    ktph2    = ktph * ktph;
    krr2     = krr * krr;
    kthth2   = kthth * kthth;
    kphph2   = kphph * kphph;

    gltt    = hphph * ktph2 + htt * ktt2;
    gltph   = ktph * (ktt * htt + kphph * hphph);
    glrr    = hrr * krr2;
    glthth  = hthth * kthth2;
    glphph  = htt * ktph2 + hphph * kphph2;



}

void MetricGlampedakis::calcGlampedakisDiff(const double* pos) {
    glampdiffcalled++;
    if ((oldglampdiffpos[0] == pos[1]) && (oldglampdiffpos[1] == pos[2])) {
        glampdiffskip++;
        return;
    }
    oldglampdiffpos[0] = pos[1];
    oldglampdiffpos[1] = pos[2];
    calcKerr(pos);
    calcGlampedakis(pos);

    M3   = M2 * M;
    M4   = M3 * M;
    M5   = M4 * M;
    r3   = r2 * r;
    r4   = r3 * r;
    r5   = r4 * r;
    sth3 = sth2 * sth;

    F1dr    = 5.0 / 8.0 * (26.0 * M3 * r2 - 4.0 * M4 * r - 24.0 * M2 * r3 + 6.0 * M * r4 + 4.0 * M5 + lgr *
                           (-3.0 * r5 + 15.0 * r4 * M - 24.0 * r3 * M2 + 12.0 * r2 * M3)) / (M2 * r2 * (-r + 2.0 * M)
                                   * (-r + 2.0 * M));

    F2dr    = -(5.0 / 8.0) * (4.0 * M3 * r + 4.0 * M4 + 6.0 * M2 * r2 - 6.0 * M * r3 + lgr * (3.0 * r4 - 6.0 * r3 * M))
              / (-r + 2.0 * M) / M2 / r2;

    cfdth = 6.0 * cth * sth;
    fSchwdr  = rs / r2;

    httdr   = cf / fSchw * (F1dr - F1 * fSchwdr / fSchw);
    httdth  = cfdth * F1 / fSchw;

    hrrdr   = cf * (F1dr * fSchw + F1 * fSchwdr);
    hrrdth  = cfdth * F1 * fSchw;

    hththdr = -cf * (F2dr * r - 2.0 * F2) / r3;
    hththdth = -cfdth * F2 / r2;

    hphphdr = -cf * (F2dr * r - 2.0 * F2) / r3 / sth2;
    hphphdth = F2 * (-cfdth * sth + 2.0 * cf * cth) / r2 / sth3;

    double ttphph   = ktt * htt + kphph * hphph;
    double ttphphdr = kttdr * htt + ktt * httdr + kphphdr * hphph + kphph * hphphdr;
    double ttphphdth = kttdth * htt + ktt * httdth + kphphdth * hphph + kphph * hphphdth;


    glttdr      = hphphdr * ktph2 + 2.0 * hphph * ktph * ktphdr + httdr * ktt2 + 2.0 * htt * kttdr * ktt;
    glttdth     = hphphdth * ktph2 + 2.0 * hphph * ktph * ktphdth + httdth * ktt2 + 2.0 * htt * kttdth * ktt;;

    gltphdr     = ktphdr * ttphph + ktph * ttphphdr;
    gltphdth    = ktphdth * ttphph + ktph * ttphphdth;

    glrrdr      = hrrdr * krr2 + hrr * 2.0 * krr * krrdr;
    glrrdth     = hrrdth * krr2 + hrr * 2.0 * krr * krrdth;

    glththdr    = hththdr * kthth2 + hthth * 2.0 * kthth * kththdr;
    glththdth   = hththdth * kthth2 + hthth * 2.0 * kthth * kththdth;

    glphphdr    = httdr * ktph2 + htt * 2.0 * ktph * ktphdr + hphphdr * kphph2 + hphph * 2.0 * kphph * kphphdr;
    glphphdth   = httdth * ktph2 + htt * 2.0 * ktph * ktphdth + hphphdth * kphph2 + hphph * 2.0 * kphph * kphphdth;

}

void MetricGlampedakis::calcGlampedakisDiff2(const double* pos) {
    glampdiff2called++;
    if ((oldglampdiff2pos[0] == pos[1]) && (oldglampdiff2pos[1] == pos[2])) {
        glampdiff2skip++;
        return;
    }
    oldglampdiff2pos[0] = pos[1];
    oldglampdiff2pos[1] = pos[2];
    calcKerrDiff2(pos);
    calcGlampedakisDiff(pos);

    M6   = M5 * M;
    r6   = r5 * r;

    sth3 = sth2 * sth;
    sth4 = sth3 * sth;

    double lgr = log(r / (r - 2.0 * M));

    F1d2r   = -5.0 / 8.0 * (44.0 * M3 * r3 - 12.0 * M4 * r2 - 24.0 * M5 * r - 30.0 * M2 * r4 + 6.0 * M * r5 + 16.0 * M6
                            + lgr * (-3.0 * r6 + 18.0 * r5 * M - 36.0 * r4 * M2 + 24.0 * r3 * M3))
              / (M2 * r3 * pow((-r + 2.0 * M), 3.0));
    F2d2r   = 5.0 / 8.0 * (-8.0 * M3 * r2 - 4.0 * M4 * r + 16.0 * M5 + 18.0 * M2 * r3 - 6.0 * M * r4
                           + lgr * (3.0 * r5 - 12.0 * r4 * M + 12.0 * r3 * M2)) / (M2 * r3 * pow((-r + 2.0 * M), 2.0));

    cfd2th      = 6.0 * (cth2 - sth2);
    fSchwd2r    = -2.0 * rs / r3;

    httd2r      = -cf * (-F1d2r * fSchw * fSchw + 2.0 * F1dr * fSchwdr * fSchw - 2.0 * F1 * fSchwdr * fSchwdr + F1 * fSchwd2r * fSchw) / pow(fSchw, 3.0);
    httd2th     = cfd2th * F1 / fSchw;
    httdrdth    =  cfdth * (F1dr * fSchw - F1 * fSchwdr) / fSchw / fSchw;

    hrrd2r      = cf * (F1d2r * fSchw + 2.0 * F1dr * fSchwdr + F1 * fSchwd2r);
    hrrd2th     = cfd2th * F1 * fSchw;
    hrrdrdth    = cfdth * (F1dr * fSchw + F1 * fSchwdr);

    hththd2r      = cf * (-F2d2r * r2 + 4.0 * F2dr * r - 6.0 * F2) / r4;
    hththd2th     = -cfd2th * F2 / r2;
    hththdrdth    = -cfdth * (F2dr * r - 2.0 * F2) / r3;

    hphphd2r      = hththd2r / sth2;
    hphphd2th     = -F2 * (cfd2th * sth2 - 4.0 * cfdth * cth * sth + 6.0 * cf * cth2 + 2.0 * cf * sth2) / r2 / sth4;
    hphphdrdth    = (F2dr * r - 2.0 * F2) * (-cfdth * sth + 2.0 * cf * cth) / r3 / sth3;

    //std::cout.precision(12);
    //std::cout << "F1: " << F1 << " F2: " << F2 << " F1dr: " << F1dr << " F2dr: " << F2dr << " F1d2r: " << F1d2r << " F2d2r: " << F2d2r << "\n";
    //std::cout << "htt: " << htt << " hrr: " << hrr << " hthth: " << hthth << " hphph: " << hphph <<  "\n";
    //std::cout << "htt: " << httd2r << " " << httd2th << " " << httdrdth <<  "\n";
    //std::cout << "hphph: " << hphphd2r << " " << hphphd2th << " " << hphphdrdth <<  "\n";


    double kttdr2   = kttdr * kttdr;
    double ktphdr2  = ktphdr * ktphdr;

    double kttdth2  = kttdth * kttdth;
    double ktphdth2 = ktphdth * ktphdth;

    double krrdr2   = krrdr * krrdr;
    double krrdth2  = krrdth * krrdth;

    double kththdr2 = kththdr * kththdr;
    double kththdth2 = kththdth * kththdth;

    double kphphdr2 = kphphdr * kphphdr;
    double kphphdth2 = kphphdth * kphphdth;

    /*
    glttd2r     = hphphd2r * ktph2 + 4.0 * hphphdr * ktph * ktphdr + 2.0 * hphph * (ktphdr2 + ktph * ktphd2r)
            + httd2r * ktt2 + 4.0 * httdr * ktt * kttdr + 2.0 * htt * (kttdr2 + ktt * kttd2r);
    glttd2th    = hphphd2th * ktph2 + 4.0 * hphphdth * ktph * ktphdth + 2.0 * hphph * (ktphdth2 + ktph * ktphd2th)
            + httd2th * ktt2 + 4.0 * httdth * ktt * kttdth + 2.0 * htt * (kttdth2 + ktt * kttd2th);
    glttdrdth   = hphphdrdth * ktph2 + 2.0 * hphphdr * ktph * ktphdth + 2.0 * hphphdth * ktph * ktphdr
            + 2.0 * hphph * (ktphdth * ktphdr + ktph * ktphdrdth) + httdrdth * ktt2 + 2.0 * httdr * ktt * kttdth + 2.0 * httdth * ktt * kttdr
            + 2.0 * htt * (kttdth * kttdr + ktt * kttdrdth);

    */

    glttd2r     = hphphd2r * ktph2 + 4.0 * hphphdr * ktph * ktphdr + 2.0 * hphph * (ktphdr2 + ktph * ktphd2r)
                  + httd2r * ktt2 + 4.0 * httdr * ktt * kttdr + 2.0 * htt * (kttdr2 + ktt * kttd2r);
    glttd2th    = hphphd2th * ktph2 + 4.0 * hphphdth * ktph * ktphdth + 2.0 * hphph * (ktphdth2 + ktph * ktphd2th)
                  + httd2th * ktt2 + 4.0 * httdth * ktt * kttdth + 2.0 * htt * (kttdth2 + ktt * kttd2th);
    glttdrdth   = hphphdrdth * ktph2 + 2.0 * hphphdr * ktph * ktphdth + 2.0 * hphphdth * ktph * ktphdr
                  + 2.0 * hphph * (ktphdth * ktphdr + ktph * ktphdrdth) + httdrdth * ktt2 + 2.0 * httdr * ktt * kttdth + 2.0 * httdth * ktt * kttdr
                  + 2.0 * htt * (kttdth * kttdr + ktt * kttdrdth);

    gltphd2r     = ktphd2r * ktt * htt + ktphd2r * kphph * hphph + 2 * ktphdr * kttdr * htt + 2 * ktphdr * ktt * httdr + 2 * ktphdr * kphphdr * hphph + 2 * ktphdr * kphph * hphphdr + ktph * kttd2r * htt + 2 * ktph * kttdr * httdr + ktph * ktt * httd2r + ktph * kphphd2r * hphph + 2 * ktph * kphphdr * hphphdr + ktph * kphph * hphphd2r;
    gltphd2th    = ktphd2th * ktt * htt + ktphd2th * kphph * hphph + 2 * ktphdth * kttdth * htt + 2 * ktphdth * ktt * httdth + 2 * ktphdth * kphphdth * hphph + 2 * ktphdth * kphph * hphphdth + ktph * kttd2th * htt + 2 * ktph * kttdth * httdth + ktph * ktt * httd2th + ktph * kphphd2th * hphph + 2 * ktph * kphphdth * hphphdth + ktph * kphph * hphphd2th;
    gltphdrdth   = ktphdrdth * ktt * htt + ktphdrdth * kphph * hphph + ktphdr * kttdth * htt + ktphdr * ktt * httdth + ktphdr * kphphdth * hphph + ktphdr * kphph * hphphdth + ktphdth * kttdr * htt + ktphdth * ktt * httdr + ktphdth * kphphdr * hphph + ktphdth * kphph * hphphdr + ktph * kttdrdth * htt + ktph * kttdr * httdth + ktph * kttdth * httdr + ktph * ktt * httdrdth + ktph * kphphdrdth * hphph + ktph * kphphdr * hphphdth + ktph * kphphdth * hphphdr + ktph * kphph * hphphdrdth;

    glrrd2r      = hrrd2r * krr2 + 4.0 * hrrdr * krr * krrdr + 2 * hrr * krrdr2 + 2.0 * hrr * krr * krrd2r;
    glrrd2th     = hrrd2th * krr2 + 4 * hrrdth * krr * krrdth + 2 * hrr * krrdth2 + 2 * hrr * krr * krrd2th;
    glrrdrdth    = hrrdrdth * krr2 + 2 * hrrdr * krr * krrdth + 2 * hrrdth * krr * krrdr + 2 * hrr * krrdth * krrdr + 2 * hrr * krr * krrdrdth;

    glththd2r    = hththd2r * kthth2 + 4 * hththdr * kthth * kththdr + 2 * hthth * kththdr2 + 2 * hthth * kthth * kththd2r;
    glththd2th   = hththd2th * kthth2 + 4 * hththdth * kthth * kththdth + 2 * hthth * kththdth2 + 2 * hthth * kthth * kththd2th;
    glththdrdth  = hththdrdth * kthth2 + 2 * hththdr * kthth * kththdth + 2 * hththdth * kthth * kththdr + 2 * hthth * kththdth * kththdr; //+2*hthth*kthth*kththdrdth;kththdrdth null

    glphphd2r    = httd2r * ktph2 + 4 * httdr * ktph * ktphdr + 2 * htt * ktphdr2 + 2 * htt * ktph * ktphd2r + hphphd2r * kphph2 + 4 * hphphdr * kphph * kphphdr + 2 * hphph * kphphdr2 + 2 * hphph * kphph * kphphd2r;
    glphphd2th   = httd2th * ktph2 + 4 * httdth * ktph * ktphdth + 2 * htt * ktphdth2 + 2 * htt * ktph * ktphd2th + hphphd2th * kphph2 + 4 * hphphdth * kphph * kphphdth + 2 * hphph * kphphdth2 + 2 * hphph * kphph * kphphd2th;
    glphphdrdth  = httdrdth * ktph2 + 2 * httdr * ktph * ktphdth + 2 * httdth * ktph * ktphdr + 2 * htt * ktphdth * ktphdr + 2 * htt * ktph * ktphdrdth + hphphdrdth * kphph2 + 2 * hphphdr * kphph * kphphdth + 2 * hphphdth * kphph * kphphdr + 2 * hphph * kphphdth * kphphdr + 2 * hphph * kphph * kphphdrdth;


}

// ********************************* Metric routines **************************
void MetricGlampedakis::calcgComps(const double* pos) {
    if ((compoldpos[0] == pos[1]) && (compoldpos[1] == pos[2])) {
        return;
    }
    compoldpos[0] = pos[1];
    compoldpos[1] = pos[2];

    calcKerr(pos);
    calcGlampedakis(pos);

    gtt     = ktt + mEpsilon * gltt;
    gtph    = ktph + mEpsilon * gltph;
    grr     = krr + mEpsilon * glrr;
    gthth   = kthth + mEpsilon * glthth;
    gphph   = kphph + mEpsilon * glphph;

}

void MetricGlampedakis::calcgCompsDiff(const double* pos) {

    calcKerrDiff(pos);
    calcGlampedakisDiff(pos);

    gttdr   = kttdr + mEpsilon * glttdr;
    gttdth  = kttdth + mEpsilon * glttdth;

    gtphdr  = ktphdr + mEpsilon * gltphdr;
    gtphdth = ktphdth + mEpsilon * gltphdth;

    grrdr   = krrdr + mEpsilon * glrrdr;
    grrdth  = krrdth + mEpsilon * glrrdth;

    gththdr   = kththdr + mEpsilon * glththdr;
    gththdth  = kththdth + mEpsilon * glththdth;

    gphphdr   = kphphdr + mEpsilon * glphphdr;
    gphphdth  = kphphdth + mEpsilon * glphphdth;

}

void MetricGlampedakis::calcgCompsDiff2(const double* pos) {

    calcKerrDiff2(pos);
    calcGlampedakisDiff2(pos);

    gttd2r   = kttd2r + mEpsilon * glttd2r;
    gttd2th  = kttd2th + mEpsilon * glttd2th;
    gttdrdth  = kttdrdth + mEpsilon * glttdrdth;

    gtphd2r  = ktphd2r + mEpsilon * gltphd2r;
    gtphd2th = ktphd2th + mEpsilon * gltphd2th;
    gtphdrdth = ktphdrdth + mEpsilon * gltphdrdth;

    grrd2r   = krrd2r + mEpsilon * glrrd2r;
    grrd2th  = krrd2th + mEpsilon * glrrd2th;
    grrdrdth  = krrdrdth + mEpsilon * glrrdrdth;

    gththd2r   = kththd2r + mEpsilon * glththd2r;
    gththd2th  = kththd2th + mEpsilon * glththd2th;
    gththdrdth  = /*kththdrdth + //+2*hthth*kthth*kththdrdth;kththdrdth null  */ mEpsilon * glththdrdth;

    gphphd2r   = kphphd2r + mEpsilon * glphphd2r;
    gphphd2th  = kphphd2th + mEpsilon * glphphd2th;
    gphphdrdth  = kphphdrdth + mEpsilon * glphphdrdth;
}



/*!
 */
void MetricGlampedakis::setStandardValues() {
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


/*!
 */
void MetricGlampedakis::initToZero() {
    for (size_t i = 0; i < 4; ++i) {
        for (size_t j = 0; j < 4; ++j) {
            g_compts[i][j] = 0.0;
            for (size_t k = 0; k < 4; ++k) {
                christoffel[i][j][k] = 0.0;
                for (size_t l = 0; l < 4; ++l) {
                    chrisD[i][j][k][l] = 0.0;
                    riem[i][j][k][l]   = 0.0;
                }
            }
        }
    }

}

} // end namespace m4d
