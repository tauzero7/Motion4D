// -------------------------------------------------------------------------------
/*
   m4dMetricTomimatsuSato.cpp

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

#include "m4dMetricTomimatsuSato.h"

namespace m4d {


/*! Standard constructor for the TomimatsuSato metric.
 *
 */
MetricTomimatsuSato::MetricTomimatsuSato(double k0, double p0, double q0) {
    mMetricName  = "TomimatsuSato";
    setCoordType(enum_coordinate_prolatespheroidal);

    mPhysicalUnits = enum_physical_constants_geom;
    mSpeedOfLight = 1.0;
    mGravConstant = 1.0;

    addParam("k", mk);
    addParam("p", mp);
    addParam("q", mq);


    setStandardValues();
    initToZero();

    k = k0;
    p = p0;
    q = q0;

    mprolSphfac = k;
    mk = k0;
    mp = p0;
    mq = q0;

    k2 = k * k;
    p2 = p * p;
    q2 = q * q;

    metricskip = 0;
    metriccalled = 0;

    christoffelskip = 0;
    christoffelcalled = 0;

    chrisDskip = 0;
    chrisDcalled = 0;


}

MetricTomimatsuSato::~MetricTomimatsuSato() {
    std::cout << " comp: " << compskip << "/" << compcalled << ", diff: " <<  compdiffskip << "/" << compdiffcalled
              << ", diff2: " <<  compdiff2skip << "/" << compdiff2called << "\n";
}

#include <iostream>

// *********************************** public methods ******************************
/*! Calculate the contravariant metric components at position 'pos'.
 *
 *  \param pos : pointer to position.
 */
bool MetricTomimatsuSato::calculateMetric(const double* pos) {

    metriccalled++;
    if ((oldmetricpos[0] == pos[1]) && (oldmetricpos[1] == pos[2])) {
        metricskip++;
        return true;
    }
    oldmetricpos[0] = pos[1];
    oldmetricpos[0] = pos[2];

    calcTomimatsuSato(pos);

    g_compts[0][0] = g44;
    g_compts[0][3] = g34;
    g_compts[1][1] = g11;
    g_compts[2][2] = g22;
    g_compts[3][0] = g34;
    g_compts[3][3] = g33;

    return true;
}

/*! Calculate the Christoffel symbols of the second kind at position 'pos'.
 *
 *  \param pos : pointer to position.
 */
bool MetricTomimatsuSato::calculateChristoffels(const double* pos) {
    christoffelcalled++;
    if ((oldchristoffelpos[0] == pos[1]) && (oldchristoffelpos[1] == pos[2])) {
        christoffelskip++;
        return true;
    }
    oldchristoffelpos[0] = pos[1];
    oldchristoffelpos[1] = pos[2];


    calcTomimatsuSato(pos);
    calcTomimatsuSatoDiff(pos);

    double fac = g33 * g44 - g34 * g34;
    //111  1
    christoffel[1][1][1] = 0.5 / g11 * g11dx;
    //112  2
    christoffel[1][1][2] = -0.5 / g22 * g11dy;
    //121  3
    christoffel[1][2][1] = 0.5 / g11 * g11dy;
    christoffel[2][1][1] = christoffel[1][2][1];
    //122  4
    christoffel[1][2][2] = 0.5 / g22 * g22dx;
    christoffel[2][1][2] = christoffel[1][2][2];
    //133  5
    christoffel[1][3][3] = 0.5 * (g44 * g33dx - g34 * g34dx) / fac;
    christoffel[3][1][3] = christoffel[1][3][3];
    //130  6
    christoffel[1][3][0] = 0.5 * (-g34 * g33dx + g33 * g34dx) / fac;
    christoffel[3][1][0] = christoffel[1][3][0];
    //103  7
    christoffel[1][0][3] = 0.5 * (g44 * g34dx - g34 * g44dx) / fac;
    christoffel[0][1][3] = christoffel[1][0][3];
    //100  8
    christoffel[1][0][0] = 0.5 * (-g34 * g34dx + g33 * g44dx) / fac;
    christoffel[0][1][0] = christoffel[1][0][0];
    //221  9
    christoffel[2][2][1] = -0.5 / g11 * g22dx;
    //222  10
    christoffel[2][2][2] = 0.5 / g22 * g22dy;
    //233  11
    christoffel[2][3][3] = 0.5 * (g44 * g33dy - g34 * g34dy) / fac;
    christoffel[3][2][3] = christoffel[2][3][3];
    //230  12
    christoffel[2][3][0] = 0.5 * (-g34 * g33dy + g33 * g34dy) / fac;
    christoffel[3][2][0] = christoffel[2][3][0];
    //203  13
    christoffel[2][0][3] = 0.5 * (g44 * g34dy - g34 * g44dy) / fac;
    christoffel[0][2][3] = christoffel[2][0][3];
    //200  14
    christoffel[2][0][0] = 0.5 * (-g34 * g34dy + g33 * g44dy) / fac;
    christoffel[0][2][0] = christoffel[2][0][0];
    //331  15
    christoffel[3][3][1] = -0.5 / g11 * g33dx;
    //332  16
    christoffel[3][3][2] = -0.5 / g22 * g33dy;
    //301  17
    christoffel[3][0][1] = -0.5 / g11 * g34dx;
    christoffel[0][3][1] = christoffel[3][0][1];
    //302  18
    christoffel[3][0][2] = -0.5 / g22 * g34dy;
    christoffel[0][3][2] = christoffel[3][0][2];
    //001  19
    christoffel[0][0][1] = -0.5 / g11 * g44dx;
    //002  20
    christoffel[0][0][2] = -0.5 / g22 * g44dy;

    /*
        std::cout << "Christoffel:\n";
        std::cout.precision(10);
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
bool MetricTomimatsuSato::calculateChrisD(const double* pos) {
    chrisDcalled++;
    if ((oldchrisDpos[0] == pos[1]) && (oldchrisDpos[1] == pos[2])) {
        chrisDskip++;
        return true;
    }
    oldchrisDpos[0] = pos[1];
    oldchrisDpos[1] = pos[2];

    calcTomimatsuSatoDiff2(pos);

    double g112 = g11 * g11;
    double g222 = g22 * g22;
    double g332 = g33 * g33;
    double g442 = g44 * g44;
    double g342 = g34 * g34;
    double g343 = g342 * g34;

    double fac2 = (g33 * g44 - g34 * g34) * (g33 * g44 - g34 * g34);

    double g11dx2 = g11dx * g11dx;
    double g22dx2 = g22dx * g22dx;
    double g33dx2 = g33dx * g33dx;
    double g34dx2 = g34dx * g34dx;
    double g44dx2 = g44dx * g44dx;

    double g11dy2 = g11dy * g11dy;
    double g22dy2 = g22dy * g22dy;
    double g33dy2 = g33dy * g33dy;
    double g34dy2 = g34dy * g34dy;
    double g44dy2 = g44dy * g44dy;

    chrisD[1][1][1][1] =  0.5 * (-g11dx2 + g11d2x * g11) / g112;


    chrisD[1][1][1][2] = -0.5 * (g11dx * g11dy - g11dxdy * g11) / g112;

    chrisD[1][1][2][1] = 0.5 * (g11dy * g22dx - g11dxdy * g22) / g222;

    chrisD[1][1][2][2] = 0.5 * (g11dy * g22dy - g11d2y * g22) / g222;

    chrisD[1][2][1][1] = -0.5 * (g11dx * g11dy - g11dxdy * g11) / g112;
    chrisD[2][1][1][1] = chrisD[1][2][1][1];

    chrisD[1][2][1][2] =  -0.5 * (g11dy2 - g11d2y * g11) / g112;
    chrisD[2][1][1][2] = chrisD[1][2][1][2];

    chrisD[1][2][2][1] = 0.5 * (-g22dx2 + g22d2x * g22) / g222;
    chrisD[2][1][2][1] = chrisD[1][2][2][1];

    chrisD[1][2][2][2] = 0.5 * (-g22dx * g22dy + g22dxdy * g22) / g222;
    chrisD[2][1][2][2] = chrisD[1][2][2][2];

    chrisD[1][3][3][1] = -0.5 * (g44dx * g33dx * g342 - g442 * g33d2x * g33 + g44 * g33d2x * g342 + g34dx2 * g33 * g44 + g34dx2 * g342
                                 + g34 * g34d2x * g33 * g44 - g343 * g34d2x + g442 * g33dx2 - 3.0 * g44 * g33dx * g34 * g34dx
                                 - g34 * g34dx * g33 * g44dx) / fac2;
    chrisD[3][1][3][1] = chrisD[1][3][3][1];

    chrisD[1][3][3][2] = -0.5 * (g44dy * g33dx * g342 - g442 * g33dxdy * g33 + g44 * g33dxdy * g342 + g34dy * g34dx * g33 * g44 + g34dy * g34dx * g342
                                 + g34 * g34dxdy * g33 * g44 - g343 * g34dxdy + g442 * g33dx * g33dy - 2.0 * g44 * g33dx * g34 * g34dy - g34 * g34dx
                                 * g44 * g33dy - g34 * g34dx * g33 * g44dy) / fac2;
    chrisD[3][1][3][2] = chrisD[1][3][3][2];

    chrisD[1][3][0][1] = 0.5 * (-g34 * g33d2x * g33 * g44 + g343 * g33d2x + g332 * g34d2x * g44 - g33 * g34d2x * g342 + g34 * g33dx2 * g44 + g34 * g33dx * g33 * g44dx
                                - 2.0 * g342 * g33dx * g34dx - g33 * g34dx * g44 * g33dx - g332 * g34dx * g44dx + 2.0 * g33 * g34dx2 * g34) / fac2;
    chrisD[3][1][0][1] = chrisD[1][3][0][1];


    chrisD[1][3][0][2] =  0.5 * (-g34dy * g33dx * g33 * g44 - g34dy * g33dx * g342 - g34 * g33dxdy * g33 * g44 + g343 * g33dxdy - g33dy * g34dx * g342 + g332 * g34dxdy * g44 - g33 * g34dxdy * g342
                                 + g34 * g33dx * g44 * g33dy + g34 * g33dx * g33 * g44dy - g332 * g34dx * g44dy + 2 * g33 * g34dx * g34 * g34dy) / fac2;

    chrisD[3][1][0][2] = chrisD[1][3][0][2];

    chrisD[1][0][3][1] = -0.5 * (-g442 * g34d2x * g33 + g44 * g34d2x * g342 + g34 * g44d2x * g33 * g44 - g343 * g44d2x + g442 * g34dx * g33dx + g44 * g34dx * g33 * g44dx - 2 * g44 * g34dx2 * g34
                                 - g34 * g44dx * g44 * g33dx - g34 * g44dx2 * g33 + 2 * g342 * g44dx * g34dx) / fac2;
    chrisD[0][1][3][1] = chrisD[1][0][3][1];


    chrisD[1][0][3][2] = -0.5 * (g44dy * g34dx * g342 - g442 * g34dxdy * g33 + g44 * g34dxdy * g342 + g34dy * g44dx * g33 * g44 + g34dy * g44dx * g342 + g34 * g44dxdy * g33 * g44 - g343 * g44dxdy
                                 + g442 * g34dx * g33dy - 2 * g44 * g34dx * g34 * g34dy - g34 * g44dx * g44 * g33dy - g34 * g44dx * g33 * g44dy) / fac2;
    chrisD[0][1][3][2] = chrisD[1][0][3][2];

    chrisD[1][0][0][1] = 0.5 * (-g34dx2 * g33 * g44 - g34dx2 * g342 - g34 * g34d2x * g33 * g44 + g343 * g34d2x - g44dx * g33dx * g342 + g332 * g44d2x * g44 - g33 * g44d2x * g342 + g44 * g33dx * g34 * g34dx
                                + 3 * g34 * g34dx * g33 * g44dx - g332 * g44dx2) / fac2;
    chrisD[0][1][0][1] = chrisD[1][0][0][1];

    chrisD[1][0][0][2] = 0.5 * (-g34dy * g34dx * g33 * g44 - g34dy * g34dx * g342 - g34 * g34dxdy * g33 * g44 + g343 * g34dxdy - g33dy * g44dx * g342 + g332 * g44dxdy * g44 - g33 * g44dxdy * g342
                                + g34 * g34dx * g44 * g33dy + g34 * g34dx * g33 * g44dy - g332 * g44dx * g44dy + 2 * g33 * g44dx * g34 * g34dy) / fac2;
    chrisD[0][1][0][2] = chrisD[1][0][0][2];

    chrisD[2][2][1][1] = -0.5 * (-g22dx * g11dx + g22d2x * g11) / g112;

    chrisD[2][2][1][2] = -0.5 * (-g11dy * g22dx + g22dxdy * g11) / g112;

    chrisD[2][2][2][1] = 0.5 * (-g22dx * g22dy + g22dxdy * g22) / g222;

    chrisD[2][2][2][2] = 0.5 * (-g22dy2 + g22d2y * g22) / g222;

    chrisD[2][3][3][1] = 0.5 * (-g33dy * g44dx * g342 + g442 * g33dxdy * g33 - g44 * g33dxdy * g342 - g34dy * g34dx * g33 * g44 - g34dy * g34dx * g342 - g34 * g34dxdy * g33 * g44 + g343 * g34dxdy
                                - g442 * g33dx * g33dy + 2 * g34 * g34dx * g44 * g33dy + g44 * g33dx * g34 * g34dy + g33 * g44dx * g34 * g34dy) / fac2;
    chrisD[3][2][3][1] = chrisD[2][3][3][1];

    chrisD[2][3][3][2] = 0.5 * (-g44dy * g33dy * g342 + g442 * g33d2y * g33 - g44 * g33d2y * g342 - g34dy2 * g33 * g44 - g34dy2 * g342 - g34 * g34d2y * g33 * g44 + g343 * g34d2y
                                - g442 * g33dy2 + 3 * g44 * g33dy * g34 * g34dy + g34 * g34dy * g33 * g44dy) / fac2;
    chrisD[3][2][3][2] = chrisD[2][3][3][2];

    chrisD[2][3][0][1] = -0.5 * (g33dy * g34dx * g33 * g44 + g33dy * g34dx * g342 + g34 * g33dxdy * g33 * g44 - g343 * g33dxdy + g34dy * g33dx * g342 - g332 * g34dxdy * g44 + g33 * g34dxdy * g342
                                 - g34 * g33dx * g44 * g33dy - g34 * g33dy * g33 * g44dx + g332 * g34dy * g44dx - 2 * g33 * g34dx * g34 * g34dy) / fac2;
    chrisD[3][2][0][1] = chrisD[2][3][0][1];


    chrisD[2][3][0][2] = -0.5 * (g34 * g33d2y * g33 * g44 - g343 * g33d2y - g332 * g34d2y * g44 + g33 * g34d2y * g342 - g34 * g33dy2 * g44 - g34 * g33dy * g33 * g44dy + 2 * g342 * g33dy * g34dy
                                 + g33 * g34dy * g44 * g33dy + g332 * g34dy * g44dy - 2 * g33 * g34dy2 * g34) / fac2;
    chrisD[3][2][0][2] = chrisD[2][3][0][2];

    chrisD[2][0][3][1] = -0.5 * (g34dy * g44dx * g342 - g442 * g34dxdy * g33 + g44 * g34dxdy * g342 + g44dy * g34dx * g33 * g44 + g44dy * g34dx * g342 + g34 * g44dxdy * g33 * g44 - g343 * g44dxdy
                                 + g442 * g34dy * g33dx - 2 * g44 * g34dx * g34 * g34dy - g34 * g44dy * g44 * g33dx - g34 * g44dx * g33 * g44dy) / fac2;
    chrisD[0][2][3][1] = chrisD[2][0][3][1];


    chrisD[2][0][3][2] = -0.5 * (-g442 * g34d2y * g33 + g44 * g34d2y * g342 + g34 * g44d2y * g33 * g44 - g343 * g44d2y + g442 * g34dy * g33dy + g44 * g34dy * g33 * g44dy - 2 * g44 * g34dy2 * g34
                                 - g34 * g44dy * g44 * g33dy - g34 * g44dy2 * g33 + 2 * g342 * g44dy * g34dy) / fac2;
    chrisD[0][2][3][2] = chrisD[2][0][3][2];

    chrisD[2][0][0][1] = 0.5 * (-g34dy * g34dx * g33 * g44 - g34dy * g34dx * g342 - g34 * g34dxdy * g33 * g44 + g343 * g34dxdy - g44dy * g33dx * g342 + g332 * g44dxdy * g44 - g33 * g44dxdy * g342
                                + g44 * g33dx * g34 * g34dy + g33 * g44dx * g34 * g34dy - g332 * g44dx * g44dy + 2 * g34 * g34dx * g33 * g44dy) / fac2;
    chrisD[0][2][0][1] = chrisD[2][0][0][1];

    chrisD[2][0][0][2] = 0.5 * (-g34dy2 * g33 * g44 - g34dy2 * g342 - g34 * g34d2y * g33 * g44 + g343 * g34d2y - g44dy * g33dy * g342 + g332 * g44d2y * g44 - g33 * g44d2y * g342 + g44 * g33dy * g34 * g34dy
                                + 3 * g34 * g34dy * g33 * g44dy - g332 * g44dy2) / fac2;
    chrisD[0][2][0][2] = chrisD[2][0][0][2];


    chrisD[3][3][1][1] = 0.5 * (g33dx * g11dx - g33d2x * g11) / g112;

    chrisD[3][3][1][2] = -0.5 * (-g33dx * g11dy + g33dxdy * g11) / g112;

    chrisD[3][3][2][1] = -0.5 * (-g33dy * g22dx + g33dxdy * g22) / g222;

    chrisD[3][3][2][2] = -0.5 * (-g33dy * g22dy + g33d2y * g22) / g222;

    chrisD[3][0][1][1] = 0.5 * (g34dx * g11dx - g34d2x * g11) / g112;
    chrisD[0][3][1][1] = chrisD[3][0][1][1];


    chrisD[3][0][1][2] = -0.5 * (-g34dx * g11dy + g34dxdy * g11) / g112;
    chrisD[0][3][1][2] = chrisD[3][0][1][2];

    chrisD[3][0][2][1] = 0.5 * (g34dy * g22dx - g34dxdy * g22) / g222;
    chrisD[0][3][2][1] = chrisD[3][0][2][1];

    chrisD[3][0][2][2] = -0.5 * (-g34dy * g22dy + g34d2y * g22) / g222;
    chrisD[0][3][2][2] = chrisD[3][0][2][2];

    chrisD[0][0][1][1] = -0.5 * (-g44dx * g11dx + g44d2x * g11) / g112;

    chrisD[0][0][1][2] = -0.5 * (-g44dx * g11dy + g44dxdy * g11) / g112;

    chrisD[0][0][2][1] = -0.5 * (-g44dy * g22dx + g44dxdy * g22) / g222;

    chrisD[0][0][2][2] = -0.5 * (-g44dy * g22dy + g44d2y * g22) / g222;

    /*
        std::cout.precision(10);
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

/*! Transform local 4-direction to coordinate 4-direction.
 *
 *  \param  pos  :  pointer to position array.
 *  \param  ldir :  pointer to local direction array.
 *  \param  dir  :  pointer to calculated coordinate direction array.
 *  \param  type :  type of tetrad.
 */
void MetricTomimatsuSato::localToCoord(const double* pos, const double* ldir, double* dir,
                                       enum_nat_tetrad_type) {
    calcTomimatsuSato(pos);

    //double omega = -g34/g33;
    double zeta = 0.0;

    double gam = 1.0 / sqrt(-(g44 + 2.0 * zeta * g34 + zeta * zeta * g33));
    double dlt = 1.0 / sqrt(g34 * g34 - g44 * g33);
    double w1 = g34 + zeta * g33;
    double w2 = g44 + zeta * g34;

    dir[0] = gam * (ldir[0] + dlt * w1 * ldir[3]);
    dir[1] = ldir[1] / sqrt(g11);
    dir[2] = ldir[2] / sqrt(g22);
    dir[3] = gam * (ldir[0] * zeta - dlt * w2 * ldir[3]);

}

/*! Transform coordinate 4-direction to local 4-direction.
 *
 *  \param  pos  :  pointer to position array.
 *  \param  cdir :  pointer to coordinate direction.
 *  \param  ldir :  pointer to calculated local direction array.
 *  \param  type :  type of tetrad.
 */
void MetricTomimatsuSato::coordToLocal(const double* , const double* , double* ,
                                       enum_nat_tetrad_type) {
    // TODO
}


/*! Test break condition.
 *
 *  \param pos    : pointer to position array.
 *  \return true  : radial position r < 0.0 or  r^2<=(1.0+eps)*rs^2.
 *  \return false : position is valid.
 */
bool MetricTomimatsuSato::breakCondition(const double* pos) {
    bool br = false;
    // x<1 oder y<1 als abbruchbedingung??
    if (isnan(pos[0]) || isnan(pos[1]) || isnan(pos[2]) || isnan(pos[3])) {
        br = true;
    } else {
        calcTomimatsuSato(pos);
        /*
                double dd = gtph*gtph - gtt*gphph;


                if (grr<=0.0 || gthth<=0.0 || dd<=0.0 ){
                    br = true;
                }
                */
    }
    return br;
}




/*! Set parameter 'pName' to 'val'.
 *
 *
 */
bool MetricTomimatsuSato::setParam(std::string pName, double val) {
    Metric::setParam(pName, val);
    if (pName == "k") {
        k = val;
        mk = k;
    } else if (pName == "q") {
        q = val;
        p = sqrt(1.0 - q * q);
        mq = q;
        mp = p;
    } else if (pName == "p") {
        p = val;
        q = sqrt(1.0 - p * p);
        mq = q;
        mp = p;
    }
    q2 = q * q;
    p2 = p * p;
    k2 = k * k;
    return true;
}

/*! Generate report.
 */
bool MetricTomimatsuSato::report(const vec4 , const vec4 , std::string &text) {
    std::stringstream ss;
    ss << "Report for Tomimatsu-Sato metric\n\tcoordinate : (t,x,y,phi)\n";
    ss << "---------------------------------------------------------------\n";
    ss << "  physical units ..................... no\n";
    ss.precision(DEF_FIXED_REPORT_PRECISION);
    ss.setf(std::ios::fixed);
    ss << "  param .............................. k   = " << k << std::endl;
    ss << "  param .............................. p = " << p << std::endl;
    ss << "  param .............................. q    = " << q << std::endl;

    text = ss.str();
    return true;
}

// ********************************* protected methods *****************************


// ********************************* TomimatsuSato routines **************************
void MetricTomimatsuSato::calcTomimatsuSato(const double* pos) {
    compcalled++;
    if ((oldcomppos[0] == pos[1]) && (oldcomppos[1] == pos[2])) {
        compskip++;
        return;
    }
    oldcomppos[0] = pos[1];
    oldcomppos[1] = pos[2];

    x = pos[1];
    y = pos[2];

    x2 = x * x;
    omy2 = 1.0 - y * y;
    x2mo = x * x - 1.0;
    x2my2 = x * x - y * y;

    mu      = p * p * x2mo * x2mo + q * q * omy2 * omy2;
    mu2     = mu * mu;

    sigma   = 2.0 * p * q * x2my2;
    sigma2  = sigma * sigma;

    nu      = 4.0 * x * (p * x2 + 2.0 * x + p);
    tau     = -4.0 * q * omy2 * (p * x + 1.0) / p;

    A       = mu2 - x2mo * omy2 * sigma2;
    f       = A / (mu2 + mu * nu - omy2 * (x2mo * sigma2 - sigma * tau));
    g       = A / (pow(p, 4) * pow(x2my2, 4));

    omega   = -k * omy2 * (x2mo * sigma * nu + mu * tau) / A ;
    omega2  = omega * omega;

    g44     = -f;
    g34     = f * omega;
    g11     = k2 / f * g * x2my2 / x2mo;
    g22     = k2 / f * g * x2my2 / omy2;
    g33     = k2 / f * x2mo * omy2 - f * omega2;
    /*
        std::cout.precision(10);
        std::cout << " k: " << k << " p: " << p << " q: " << q << " k2: " << k2 <<  " p2: " << p2 <<  " q2: " << q2 <<  "\n";
        std::cout << " x: " << x << " y = " << y << " x2mo = " << x2mo << " x2my2 = " << x2my2 <<  "\n";
        std::cout << " A: " << A << " mu = " << mu << " nu = " << nu << " sigma = " << sigma << " tau = " << tau << " omega =  " << omega << "\n";
        std::cout << "\nf: " << f << " g = " << g << "\n";
        std::cout << " g44: " << g44 << " g34 = " << g34 << " g11 = " << g11 << " g22 = " << g22 << " g33 = " << g33 << "\n\n\n";

    */

}

void MetricTomimatsuSato::calcTomimatsuSatoDiff(const double* pos) {
    compdiffcalled++;
    if ((oldcompdiffpos[0] == pos[1]) && (oldcompdiffpos[1] == pos[2])) {
        compdiffskip++;
        return;
    }
    oldcompdiffpos[0] = pos[1];
    oldcompdiffpos[1] = pos[2];

    calcTomimatsuSato(pos);


    mudx    = 4.0 * p * p * x2mo * x;
    mudy    = -4.0 * q * q * omy2 * y;
    sigmadx = 4.0 * p * q * x;
    sigmady = -4.0 * p * q * y;
    nudx    = 4.0 * p * x2 + 8.0 * x + 4.0 * p + 4.0 * x * 2.0 * (p * x + 1.0);

    taudx   = -4.0 * q * omy2;
    taudy   = 8.0 * q * y * (p * x + 1.0) / p;

    Adx     = 2.0 * (mu * mudx - x * omy2 * sigma2 - x2mo * omy2 * sigma * sigmadx);
    Ady     = 2.0 * (mu * mudy + x2mo * y * sigma2 - x2mo * omy2 * sigma * sigmady);

    B =  mu2 + mu * nu - omy2 * (x2mo * sigma2 - sigma * tau) ;

    fdx     = Adx / B - A * (2.0 * mu * mudx + mudx * nu + mu * nudx - omy2 *
                             (2.0 * x * sigma2 + 2.0 * x2mo * sigma * sigmadx - sigmadx * tau - sigma * taudx)) / B / B;
    fdy     = Ady / B - A * (2.0 * mu * mudy + mudy * nu + 2.0 * y *
                             (x2mo * sigma2 - sigma * tau) - omy2 *
                             (2.0 * x2mo * sigma * sigmady - sigmady * tau - sigma * taudy)) / B / B;

    double pow4 = pow(p * x2my2, 4);

    gdx     = Adx / pow4 - 8.0 * A * x / pow4 / x2my2;
    gdy     = Ady / pow4 + 8.0 * A * y / pow4 / x2my2;

    omegadx = -k * omy2 * (2.0 * x * sigma * nu + x2mo * sigmadx * nu + x2mo * sigma * nudx + mudx * tau + mu * taudx) / A
              + k * omy2 * (x2mo * sigma * nu + mu * tau) * Adx / A / A;
    omegady = 2.0 * k * y * (x2mo * sigma * nu + mu * tau) / A - k * omy2 * (x2mo * sigmady * nu + mudy * tau + mu * taudy) / A
              + k * omy2 * (x2mo * sigma * nu + mu * tau) * Ady / A / A;

    x4 = x2 * x2;
    y2 = y * y;
    y4 = y2 * y2;

    g44dx = -fdx;
    g44dy = -fdy;

    g34dx = fdx * omega + f * omegadx;
    g34dy = fdy * omega + f * omegady;

    g11dx = k2 * ((g * fdx - gdx * f) * (-x4 + x2 + x2 * y2 - y2) - 2.0 * g * f * x * omy2) / pow(f * x2mo, 2.0) ;
    g11dy = -k2 * ((g * fdy - gdy * f) * (x2 - y2) + 2.0 * g * f * y) / f / f / x2mo;

    g22dx = k2 * (-(g * fdx - gdx * f) * x2my2 + 2.0 * g * f * x)  / f / f / omy2 ;
    g22dy = k2 * ((g * fdy - gdy * f) * (-x2 + x2 * y2 + y2 - y4) + 2.0 * g * f * y * x2mo) / pow(f * omy2, 2.0);

    g33dx = -k2 * (fdx * (x2 - x2 * y2 - 1.0 + y2) - 2.0 * f * x * omy2)  / f / f - fdx * omega2 - 2.0 * f * omega * omegadx;
    g33dy = -k2 * (fdy * (x2 - x2 * y2 - 1.0 + y2) + 2.0 * f * y * x2mo)  / f / f - fdy * omega2 - 2.0 * f * omega * omegady;

    /*
        std::cout.precision(10);
        std::cout << " A: dx = " << Adx << " dy = " << Ady << "\n";
        std::cout << " f: dx = " << fdx << " dy = " << fdy << "\n";
        std::cout << " g: dx = " << gdx << " dy = " << gdy << "\n";
        std::cout << " omega: dx = " << omegadx << " dy = " << omegady << "\n";
        std::cout << " g44: dx = " << g44dx << " dy = " << g44dy << "\n";
        std::cout << " g34: dx = " << g34dx << " dy = " << g34dy << "\n";
        std::cout << " g11: dx = " << g11dx << " dy = " << g11dy << "\n";
        std::cout << " g22: dx = " << g22dx << " dy = " << g22dy << "\n";
        std::cout << " g33: dx = " << g33dx << " dy = " << g33dy << "\n";

    */
}

void MetricTomimatsuSato::calcTomimatsuSatoDiff2(const double* pos) {
    compdiff2called++;
    if ((oldcompdiff2pos[0] == pos[1]) && (oldcompdiff2pos[1] == pos[2])) {
        compdiff2skip++;
        return;
    }
    oldcompdiff2pos[0] = pos[1];
    oldcompdiff2pos[1] = pos[2];

    calcTomimatsuSato(pos);
    calcTomimatsuSatoDiff(pos);

    mud2x       = 8.0 * p2 * x2 + 4.0 * p2 * x2mo;
    mud2y       = 8.0 * q2 * y2 - 4.0 * q2 * omy2;

    sigmad2x    = 4.0 * p * q;
    sigmad2y    = -sigmad2x;

    nud2x       = 24.0 * p * x + 16.0;

    taud2y      = 8.0 * q * (p * x + 1.0) / p;
    taudxdy     = 8.0 * q * y;

    double mudx2 = mudx * mudx;
    double mudy2 = mudy * mudy;
    double sigmadx2 = sigmadx * sigmadx;
    double sigmady2 = sigmady * sigmady;
    double omegadx2 = omegadx * omegadx;
    double omegady2 = omegady * omegady;

    double fdx2 = fdx * fdx;
    double fdy2 = fdy * fdy;

    Ad2x        = 2.0 * (mudx2) + 2.0 * mu * mud2x - 2.0 * omy2 * sigma2 - 8.0 * x * omy2 * sigma * sigmadx
                  - 2.0 * x2mo * omy2 * (sigmadx2 + sigma * sigmad2x);
    Ad2y        = 2.0 * (mudy2) + 2.0 * mu * mud2y + 2.0 * x2mo * sigma2 + 8.0 * x2mo * y * sigma * sigmady
                  - 2.0 * x2mo * omy2 * (sigmady2 + sigma * sigmad2y);
    Adxdy       = 2.0 * mudx * mudy /*+ 2.0 * mu * mudxdy */ + 4.0 * x * y * sigma2 - 4.0 * x * omy2 * sigma * sigmady
                  + 4.0 * x2mo * y * sigma * sigmadx - 2.0 * x2mo * omy2 * sigmadx * sigmady /*- 2.0 * x2mo * omy2 * sigma * sigmadxdy */;

    //B =  mu2 + mu * nu - omy2 * (x2mo * sigma2 - sigma * tau) ;

    double tmp1 = 2.0 * mu * mudx + mudx * nu + mu * nudx - omy2
                  * (2.0 * (x * sigma2 + x2mo * sigma * sigmadx) - sigmadx * tau - sigma * taudx);
    double tmp2 = 2.0 * mu * mudy + mudy * nu + 2.0 * y * (x2mo * sigma2 - sigma * tau)
                  - omy2 * (2.0 * x2mo * sigma * sigmady - sigmady * tau - sigma * taudy);

    fd2x        = Ad2x / B - 2.0 * Adx * tmp1 / B / B + 2.0 * A * tmp1 * tmp1 / B / B / B
                  - A * (2.0 * mudx2 + 2.0 * mu * mud2x + mud2x * nu + 2.0 * mudx * nudx + mu * nud2x - omy2 *
                         (2.0 * sigma2 + 8.0 * x * sigma * sigmadx + 2.0 * x2mo * (sigmadx2 + sigma * sigmad2x) - sigmad2x * tau
                          - 2.0 * sigmadx * taudx)) / B / B;

    fd2y        = Ad2y / B - 2.0 * Ady * tmp2 / B / B + 2.0 * A * tmp2 * tmp2 / B / B / B
                  - A * (2.0 * mudy2 + 2.0 * mu * mud2y + mud2y * nu + 2.0 * x2mo * sigma2 - 2.0 * sigma * tau + 4.0 * y * (2.0 * x2mo * sigma * sigmady
                          - sigmady * tau - sigma * taudy) - omy2 * (2.0 * x2mo * (sigmady2 + sigma * sigmad2y) - sigmad2y * tau - 2.0 * sigmady * taudy - sigma * taud2y)) / B / B;

    fdxdy       = Adxdy / B  - (Adx * tmp2 + Ady * tmp1) / B / B + 2.0 * A * tmp1 * tmp2 / B / B / B
                  - A * (2.0 * mudx * mudy /*+ 2.0 * mu * mudxdy + mudxdy * nu*/ + mudy * nudx + 2.0 * y * (2.0 * x * sigma2 + 2.0 * x2mo * sigma * sigmadx
                          - sigmadx * tau - sigma * taudx) - omy2 * (4.0 * x * sigma * sigmady + 2.0 * x2mo * /*(*/ sigmadx * sigmady /*+ sigma * sigmadxdy )*/
                                  /*- sigmadxdy * tau */ - sigmadx * taudy - sigmady * taudx - sigma * taudxdy)) / B / B;

    double pow4 = pow(p * x2my2, 4.0);
    double pow5 = pow4 * x2my2;
    double pow6 = pow5 * x2my2;

    gd2x        = Ad2x / pow4 - (16.0 * Adx * x + 8.0 * A) / pow5 + 80.0 * A * x2 / pow6;
    gd2y        = Ad2y / pow4 + (16.0 * Ady * y + 8.0 * A) / pow5 + 80.0 * A * y2 / pow6;
    gdxdy       = Adxdy / pow4 +  8.0 * (Adx * y - Ady * x) / pow5 - 80.0 * A * x * y / pow6;

    double tmp3 = x2mo * sigma * nu + mu * tau;

    omegad2x = -k * omy2 * (2.0 * sigma * nu + 4.0 * x * (sigmadx * nu + sigma * nudx) + x2mo * (sigmad2x * nu + 2.0 * sigmadx * nudx + sigma * nud2x)
                            + mud2x * tau + 2.0 * mudx * taudx /*+ mu * taud2x*/) / A
               + 2.0 * k * omy2 * (2.0 * x * sigma * nu + x2mo * (sigmadx * nu + sigma * nudx) + mudx * tau + mu * taudx) * Adx / A / A
               - 2.0 * k * omy2 * tmp3 * Adx * Adx / A / A / A
               + k * omy2 * tmp3 * Ad2x / A / A;


    omegad2y = 2.0 * k * tmp3 / A + 4.0 * k * y * (x2mo * sigmady * nu + mudy * tau + mu * taudy) / A
               - 4.0 * k * y * tmp3 * Ady / A / A
               - k * omy2 * (x2mo * sigmad2y * nu + mud2y * tau + 2.0 * mudy * taudy + mu * taud2y) / A
               + 2.0 * k * omy2 * (x2mo * sigmady * nu + mudy * tau + mu * taudy) * Ady / A / A
               - 2.0 * k * omy2 * tmp3 * Ady * Ady / A / A / A
               + k * omy2 * tmp3 * Ad2y / A / A;

    omegadxdy = 2.0 * k * y * (2.0 * x * sigma * nu + x2mo * (sigmadx * nu + sigma * nudx) + mudx * tau + mu * taudx) / A
                - k * omy2 * (2.0 * x * sigmady * nu /*+ x2mo * sigmadxdy * nu */ + x2mo * sigmady * nudx /*+ mudxdy * tau */ + mudx * taudy + mudy * taudx + mu * taudxdy) / A
                + k * omy2 * (2.0 * x * sigma * nu + x2mo * (sigmadx * nu +  sigma * nudx) + mudx * tau + mu * taudx) * Ady / A / A
                - 2.0 * k * y * tmp3 * Adx / A / A
                + k * omy2 * (x2mo * sigmady * nu + mudy * tau + mu * taudy) * Adx / A / A
                - 2.0 * k * omy2 * tmp3 * Adx * Ady / A / A / A
                + k * omy2 * tmp3 * Adxdy / A / A;

    g44d2x     = -fd2x;
    g44d2y     = -fd2y;
    g44dxdy    = -fdxdy;

    g34d2x     = fd2x * omega + 2.0 * fdx * omegadx  + f * omegad2x;
    g34d2y     = fd2y * omega + 2.0 * fdy * omegady + f * omegad2y;
    g34dxdy    = fdxdy * omega + fdx * omegady + fdy * omegadx + f * omegadxdy;


    g11d2x     = k2 * (2.0 * g * x2my2 * fdx2 / (f * f * f * x2mo)
                       - (x2my2 * (2.0 * gdx * fdx + g * fd2x) + 4.0 * g * x * fdx) / (f * f * x2mo)
                       + (gd2x * x2my2 + 4.0 * gdx * x + 2.0 * g) / (f * x2mo)
                       - (x2my2 * (4.0 * gdx * x + 2.0 * g) + 8.0 * g * x2) / (f * x2mo * x2mo)
                       + 8.0 * g * x2my2 * x2 / (f * x2mo * x2mo * x2mo)
                       + 4.0 * g * x2my2 * fdx * x / (f * f * x2mo * x2mo)
                      );


    g11d2y     = k2 * (2.0 * g * x2my2 * fdy2 / (f * f * f * x2mo)
                       + (-x2my2 * (2.0 * gdy * fdy + g * fd2y) + 4.0 * g * y * fdy) / (f * f * x2mo)
                       + (gd2y * x2my2 - 4.0 * gdy * y - 2.0 * g) / (f * x2mo)
                      );
    g11dxdy    = k2 * (2.0 * g * x2my2 * fdx * fdy / (f * f * f * x2mo)
                       + (-x2my2 * (gdy * fdx + g * fdxdy + gdx * fdy) + 2.0 * g * (y * fdx - x * fdy)) / (f * f * x2mo)
                       + (gdxdy * x2my2 + 2.0 * (-gdx * y + gdy * x)) / (f * x2mo)
                       + (-2.0 * gdy * x * x2my2 + 4.0 * g * x * y) / (f * x2mo * x2mo)
                       + 2.0 * g * x2my2 * fdy * x / (f * f * x2mo * x2mo)
                      );


    g22d2x     = k2 * (2.0 * g * x2my2 * fdx2 / (f * f * f * omy2)
                       - (x2my2 * (2.0 * gdx * fdx + g * fd2x) + 4.0 * g * x * fdx) / (f * f * omy2)
                       + (gd2x * x2my2 + 4.0 * gdx * x + 2.0 * g) / (f * omy2)
                      );
    g22d2y     = k2 * (2.0 * g * x2my2 * fdy2 / (f * f * f * omy2)
                       + (-x2my2 * (2.0 * gdy * fdy + g * fd2y) + 4.0 * g * y * fdy) / (f * f * omy2)
                       + (gd2y * x2my2 - 4.0 * gdy * y - 2.0 * g) / (f * omy2)
                       + (x2my2 * (4.0 * gdy * y + 2.0 * g) - 8.0 * g * y2) / (f * omy2 * omy2)
                       + 8.0 * g * x2my2 * y2 / (f * omy2 * omy2 * omy2)
                       - 4.0 * g * x2my2 * fdy * y / (f * f * omy2 * omy2)
                      );
    g22dxdy    = k2 * (2.0 * g * x2my2 * fdx * fdy / (f * f * f * omy2)
                       + (-x2my2 * (gdy * fdx + g * fdxdy + gdx * fdy) + 2.0 * g * (y * fdx - x * fdy)) / (f * f * omy2)
                       + (gdxdy * x2my2 + 2.0 * (-gdx * y + gdy * x)) / (f * omy2)
                       + (2.0 * gdx * y * x2my2 + 4.0 * g * x * y) / (f * omy2 * omy2)
                       - 2.0 * g * x2my2 * fdx * y / (f * f * omy2 * omy2)
                      );

    g33d2x     = k2 * (2.0 * x2mo * omy2 * fdx2 / (f * f * f)
                       - omy2 * (4.0 * x * fdx + x2mo * fd2x) / (f * f)
                       + 2.0 * omy2 / f)
                 - fd2x * omega2 - 4.0 * fdx * omega * omegadx - 2.0 * f * (omegadx2 + omega * omegad2x);

    g33d2y     = k2 * (2.0 * x2mo * omy2 * fdy2 / (f * f * f)
                       + x2mo * (4.0 * y * fdy - omy2 * fd2y) / (f * f)
                       - 2.0 * x2mo / f)
                 - fd2y * omega2 - 4.0 * fdy * omega * omegady - 2.0 * f * (omegady2 + omega * omegad2y);
    g33dxdy    = k2 * (2.0 * x2mo * omy2 * fdx * fdy / (f * f * f)
                       + (x2mo * (2.0 * y * fdx - omy2 * fdxdy) - 2.0 * omy2 * x * fdy) / (f * f)
                       - 4.0 * x * y / f)
                 - fdxdy * omega2 - 2.0 * (omega * (fdx * omegady + fdy * omegadx) + f * (omegadx * omegady + omega * omegadxdy));

    /*
        std::cout.precision(10);
        std::cout << " mu: d2x = " << mud2x << " d2y = " << mud2y << "\n";
        std::cout << " sigma: d2x = " << sigmad2x << " d2y = " << sigmad2y << "\n";
        std::cout << " nud2x = " << nud2x << "\n";
        std::cout << " tau: d2y = " << taud2y << " dxdy = " << taudxdy << "\n";

        std::cout << " A: d2x = " << Ad2x << " d2y = " << Ad2y << " dxdy = " << Adxdy << "\n";
        std::cout << " f: d2x = " << fd2x << " d2y = " << fd2y << " dxdy = " << fdxdy << "\n";
        std::cout << " g: d2x = " << gd2x << " d2y = " << gd2y << " dxdy = " << gdxdy << "\n";
        std::cout << " omega: d2x = " << omegad2x << " d2y = " << omegad2y << " dxdy = " << omegadxdy << "\n";

        std::cout << " g11: d2x = " << g11d2x << " d2y = " << g11d2y << " dxdy = " << g11dxdy << "\n";
        std::cout << " g22: d2x = " << g22d2x << " d2y = " << g22d2y << " dxdy = " << g22dxdy << "\n";
        std::cout << " g33: d2x = " << g33d2x << " d2y = " << g33d2y << " dxdy = " << g33dxdy << "\n";
    */
    //std::cout << " g34: d2x = " << g34d2x << " d2y = " << g34d2y << " dxdy = " << g34dxdy << "\n";

}


/*!
 */
void MetricTomimatsuSato::setStandardValues() {
    mInitPos[0] = 0.0;
    mInitPos[1] = 5.0;
    mInitPos[2] = 0.5;
    mInitPos[3] = 0.0;
    mInitDir[0] = 1.0;
    mInitDir[1] = 0.0;
    mInitDir[2] = 0.0;

    mCoordNames[0] = std::string("t");
    mCoordNames[1] = std::string("x");
    mCoordNames[2] = std::string("y");
    mCoordNames[3] = std::string("phi");

}


/*!
 */
void MetricTomimatsuSato::initToZero() {
    for (size_t i = 0; i < 4; ++i) {
        for (size_t j = 0; j < 4; ++j) {
            g_compts[i][j] = 0.0;
            for (size_t k = 0; k < 4; ++k) {
                christoffel[i][j][k] = 0.0;
                for (size_t l = 0; l < 4; ++l) {
                    chrisD[i][j][k][l] = 0.0;
                }
            }
        }
    }

}

} // end namespace m4d
