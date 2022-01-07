/**
 * @file    testJacobi.cpp
 * @author  Thomas Mueller
 *
 * @brief Test Jacobi equation.

     Usage: \verbatim ./m4dTestJacobi <file.ini> \endverbatim

 * This file is part of the m4d-library
 */
#include <cmath>
#include <iostream>

#include <extra/m4dObject.h>
#include <m4dGlobalDefs.h>
#include <math/TransCoordinates.h>
#include <metric/m4dMetricDatabase.h>
#include <motion/m4dGeodesic.h>
#include <motion/m4dGeodesicGSL.h>
#include <motion/m4dGeodesicRK4.h>

#define DEF_USE_KERR

/* -----------------------------------------------------------------
 *                     m a i n
 * ----------------------------------------------------------------- */
int main(int argc, char* argv[])
{
    if (argc != 2) {
        fprintf(stderr, "Usage: ./m4dTestJacobi <file.ini>\n");
        return -1;
    }

    fprintf(stderr, "\n\n=================== Test Jacobi ===================\n\n");

    /* -----------------------------------------
     *    Load parameter setting from file.
     * ----------------------------------------- */
    m4d::Object mObject;
    if (!mObject.loadSettings(argv[1])) {
        fprintf(stderr, "Can't load settings file %s  or settings file is invalid!\n", argv[1]);
        return -1;
    }
    mObject.printSettings();

    mObject.currMetric->printF();
    mObject.geodSolver->printF();

    mObject.geodSolver->setInitialPosition(mObject.startPos);

    /* -----------------------------------------
     *    Set local initial direction.
     * ----------------------------------------- */
    if (mObject.type == m4d::enum_geodesic_timelike) {
        fprintf(stderr, "Error: geodesic type is timelike!\n");
        return -1;
    }

    m4d::vec4 locDir = M4D_SIGN(mObject.timeDirection) * mObject.base[0] + mObject.startDir[0] * mObject.base[1]
        + mObject.startDir[1] * mObject.base[2] + mObject.startDir[2] * mObject.base[3];
    m4d::vec4 coDir;
    mObject.currMetric->localToCoord(mObject.startPos, locDir, coDir, mObject.tetradType);
    mObject.geodSolver->setInitialDirection(coDir);
    m4d::vec3 nullDir = m4d::vec3(mObject.startDir[0], mObject.startDir[1], mObject.startDir[2]);

    m4d::vec3 locX = m4d::vec3(1.0, 0.0, 0.0);
    m4d::vec3 locY = m4d::vec3(0.0, 1.0, 0.0);
    m4d::vec3 locZ = m4d::vec3(0.0, 0.0, 1.0);
    m4d::vec3 lX, lY, lZ;

    lX = locX;
    lY = locY;
    lZ = locZ;

    // m4d::enum_break_condition  breakCond = m4d::enum_break_none;

    //  int64_t t1,t2;
    //  t1 = m4d::get_system_clock();
    mObject.geodSolver->calcSachsJacobi(mObject.startPos, coDir, nullDir, lX, lY, lZ, mObject.base[0], mObject.base[1],
        mObject.base[2], mObject.base[3], mObject.tetradType, mObject.maxNumPoints, mObject.points, mObject.dirs,
        mObject.lambda, mObject.sachs1, mObject.sachs2, mObject.jacobi, mObject.maxJacobi);
    //  t2 = m4d::get_system_clock();

    // fprintf(stderr,"#points: %d  calculated in %f seconds.\n",(int)mObject.points.size(),(t2-t1)*1e-6);

    double t, r, phi, theta, dp, dm, mu, angle, elipt;

    unsigned int j = 0;
    while (j < mObject.points.size()) {
        t = mObject.points[j][0];
        r = mObject.points[j][1];
        theta = mObject.points[j][2];
        phi = mObject.points[j][3];

        dp = mObject.jacobi[j][0];
        dm = mObject.jacobi[j][1];
        mu = mObject.jacobi[j][2];
        angle = mObject.jacobi[j][3] * RAD_TO_DEG;
        elipt = mObject.jacobi[j][4];

        // fprintf(stdout,"%5d  %12.8f %12.8f %12.8f   %12.8e %12.8e %12.8e
        // %12.8e\n",j,t,r,theta,mObject.sachs1[j][0],mObject.sachs1[j][1],mObject.sachs1[j][2],mObject.sachs1[j][3]);
        fprintf(stdout, "%5d %12.8e   %12.8e %12.8e %12.8e %12.12e   %12.8e %12.8e %12.8e %12.8e %12.8e\n", j,
            mObject.lambda[j], t, r, theta, phi, dp, dm, mu, angle, elipt);
        j++;
    }

    return 1;
}
