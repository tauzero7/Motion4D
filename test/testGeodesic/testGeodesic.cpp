/**
 * @file    testGeodesic.h
 * @author  Thomas Mueller
 *
 * @brief  Test geodesic calculation within the Schwarzschild spacetime in cartesian coordinates.

     Usage: \verbatim ./m4dTestGeodesic \endverbatim

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
        fprintf(stderr, "Usage: ./m4dTestGeodesic\n");
        return -1;
    }

    fprintf(stderr, "\n\n=================== Test Geodesic ===================\n\n");

    m4d::Object mObject;

    /* -----------------------------------------
     *   Initialize metric.
     * ----------------------------------------- */
    m4d::MetricDatabase metricDB;
    m4d::Metric* metric = metricDB.getMetric("SchwarzschildIsotropic");
    metric->printF();

    /* -----------------------------------------
     *    Set geodesic integrator.
     * ----------------------------------------- */
    m4d::GeodesicGSL* geod = new m4d::GeodesicGSL(metric, gsl_odeiv_step_rkf45, 4, m4d::enum_geodesic_lightlike);
    // m4d::GeodesicRK4* geod = new m4d::GeodesicRK4( metric, m4d::enum_geodesic_lightlike );
    geod->setEpsilons(1.0e-8, 0.0);
    geod->setStepSizeControlled(true);
    double boxSize = 20.0;
    geod->setBoundingBox(
        m4d::vec4(-DBL_MAX, -boxSize, -boxSize, -boxSize), m4d::vec4(DBL_MAX, boxSize, boxSize, boxSize));

    m4d::vec4 pos, locDir, coDir;
    locDir = m4d::vec4(1.0, -1.0, 0.0, 0.0);

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
        pos = m4d::vec4(0.0, 10.0, y[i], 0.0);
        metric->localToCoord(pos, locDir, coDir, m4d::enum_nat_tetrad_default);

        geod->setInitialPosition(pos);
        geod->setInitialDirection(coDir);

        // Calculate the geodesic !
        std::vector<m4d::vec4> points;
        std::vector<m4d::vec4> dirs;
        std::vector<double> lambda;
        std::cerr << m4d::stl_break_condition[geod->calculateGeodesic(pos, coDir, 1000, points, dirs, lambda)]
                  << std::endl;
        std::cerr << "#points: " << points.size() << std::endl;

        unsigned int j = 0;
        while (j < points.size()) {
            fprintf(
                fptr, "%5d  %12.8f %12.8f %12.8f %12.8f\n", j, points[j][0], points[j][1], points[j][2], points[j][3]);
            j++;
        }
        fprintf(fptr, "\n");
        fclose(fptr);
    }

    delete geod;
    delete metric;
    return 1;
}
