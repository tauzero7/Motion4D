/**
 * @file    testSimpleGeodesic.h
 * @author  Thomas Mueller
 *
 * @brief Test geodesic calculation within the Schwarzschild spacetime
           in cartesian coordinates.

     Usage: \verbatim ./m4dTestSimpleGeodesic \endverbatim

 * This file is part of the m4d-library.
 */
#include <cmath>
#include <iostream>

#include <extra/m4dObject.h>
#include <m4dGlobalDefs.h>
#include <math/TransCoordinates.h>
#include <metric/m4dMetricDatabase.h>
#include <motion/m4dGeodesic.h>
#include <motion/m4dGeodesicGSL.h>

/* -----------------------------------------------------------------
 *                     m a i n
 * ----------------------------------------------------------------- */
int main(int argc, char* argv[])
{
    if (argc != 1) {
        fprintf(stderr, "Usage: ./m4dTestSimpleGeodesic\n");
        return -1;
    }

    fprintf(stderr, "\n\n=================== Test simple Geodesic ===================\n\n");

    m4d::Object mObject;

    // Set metric
    // mObject.setMetric("SchwarzschildCart");
    mObject.setMetric("TeoSimpleWH");
    mObject.currMetric->printF();

    // Set geodesic integrator.
    mObject.setSolver("GSL_RK_Cash-Karp");
    mObject.setSolverParam("eps_a", 1e-8);
    mObject.setSolverParam("stepctrl", true);
    double boxSize = 100.0;
    mObject.setSolverParam("lower_bb", -DBL_MAX, -boxSize, -boxSize, -boxSize);
    mObject.setSolverParam("upper_bb", DBL_MAX, boxSize, boxSize, boxSize);
    mObject.geodSolver->printF();

    mObject.setInitialPosition(0.0, 10.0, M_PI_2, 0.0);
    mObject.setInitialLocalNullDirection(m4d::enum_time_forward, -1.0, 0.0, 1.0);
    mObject.maxNumPoints = 1000;
    m4d::enum_break_condition bcd;
    bcd = mObject.calculateGeodesic();
    std::cerr << m4d::stl_break_condition[(int)bcd] << std::endl;
    exit(1);

    /* -----------------------------------------
     *    Loop over several initial positions.
     * ----------------------------------------- */
    const int NUM_INIT_POS = 15;
    double y[NUM_INIT_POS] = { 0.0, 1.0, 2.0, 3.0, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 9.0, 10.0 };

    FILE* fptr;
    char fname[32];

    for (int i = 0; i < NUM_INIT_POS; i++) {
#ifdef _WIN32
        sprintf_s(fname, "points_%d.dat", i);
#else
        sprintf(fname, "points_%d.dat", i);
#endif

#ifdef _WIN32
        fopen_s(&fptr, fname, "w");
#else
        fptr = fopen(fname, "w");
#endif
        if (fptr == NULL) {
            fprintf(stderr, "Cannot open file for output!\n");
            return -1;
        }

        mObject.setInitialPosition(0.0, 10.0, y[i], 0.0);
        mObject.setInitialLocalNullDirection(m4d::enum_time_forward, -1.0, 0.0, 0.0);

        m4d::enum_break_condition bc;
        bc = mObject.calculateGeodesic();
        std::cerr << m4d::stl_break_condition[(int)bc];
        std::cerr << "#points: " << mObject.points.size() << std::endl;

        unsigned int j = 0;
        while (j < mObject.points.size()) {
            fprintf(fptr, "%5d  %12.8f %12.8f %12.8f %12.8f\n", j, mObject.points[j][0], mObject.points[j][1],
                mObject.points[j][2], mObject.points[j][3]);
            j++;
        }
        fprintf(fptr, "\n");
        fclose(fptr);
    }

    mObject.resetAll();
    return 1;
}
