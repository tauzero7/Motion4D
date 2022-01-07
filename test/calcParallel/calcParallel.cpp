/**
 * @file    calcParallel.cpp
 * @author  Thomas Mueller
 *
 * @brief Calculate parallel transport of local tetrad  with parameters
           from a setting file.

     Usage: \verbatim ./CalcParallel  <file.ini> \endverbatim

     The parameter file is handled by the m4d::Object class.

     Output:  i  pos[0] pos[1] pos[2] pos[3]  base0[0] base0[1] base0[2] base0[3]  base1[0] base1[1]....

     where
     <ul>
       <li>i        : number of point
       <li>pos[n]   : n-th component of the point
       <li>base0[n] : n-th component of the first tetrad vector in coordinates
       <li>base1[n] : n-th component of the second tetrad vector in coordinates
       <li>etc.
     </ul>

 * This file is part of the m4d-library.
 */

#include <cmath>
#include <iostream>

#include <extra/m4dObject.h>
#include <extra/m4dUtilities.h>
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
    if (argc != 2) {
        fprintf(stderr, "Usage: ./m4dCalcParallel <file.ini>\n");
        return -1;
    }

    fprintf(stderr, "\n\n=================== Calculate parallel transport ===================\n\n");

    /* -----------------------------------------
     *    Load parameter setting from file.
     * ----------------------------------------- */
    m4d::Object mObject;
    if (!mObject.loadSettings(argv[1])) {
        fprintf(stderr, "Can't load settings file %s  or settings file is invalid!\n", argv[1]);
        return -1;
    }
    mObject.printSettings();

    /* -----------------------------------------
     *    Set geodesic integrator.
     * ----------------------------------------- */
    const gsl_odeiv_step_type* step_type = NULL;
    switch (mObject.geodSolverType) {
        default:
        case m4d::gsIrk4:
            break;
        case m4d::gsIgslrk2:
            step_type = gsl_odeiv_step_rk2;
            break;
        case m4d::gsIgslrk4:
            step_type = gsl_odeiv_step_rk4;
            break;
        case m4d::gsIgslfehlberg:
            step_type = gsl_odeiv_step_rkf45;
            break;
        case m4d::gsIgslcash:
            step_type = gsl_odeiv_step_rkck;
            break;
        case m4d::gsIgslprinc:
            step_type = gsl_odeiv_step_rk8pd;
            break;
        case m4d::gsIi2:
            step_type = gsl_odeiv_step_rk2imp;
            break;
        /*
        case m4d::gsIbs:
          step_type = gsl_odeiv_step_bsimp;
          break;
        */
        case m4d::gsIm1:
            step_type = gsl_odeiv_step_gear1;
            break;
        case m4d::gsIm2:
            step_type = gsl_odeiv_step_gear2;
            break;
    }

    m4d::Geodesic* geod;
    if (mObject.geodSolverType == m4d::gsIrk4) {
        geod = new m4d::GeodesicRK4(mObject.currMetric, mObject.type);
    }
#ifdef USE_DP_INT
    else if (mObject.geodSolverType == m4d::gsIdp54) {
        geod = new m4d::GeodesicDP54(mObject.currMetric, mObject.type);
    }
    else if (mObject.geodSolverType == m4d::gsIdp65) {
        geod = new m4d::GeodesicDP65(mObject.currMetric, mObject.type);
    }
#endif
    else {
        geod = new m4d::GeodesicGSL(mObject.currMetric, step_type, mObject.type);
    }
    geod->setEpsilons(mObject.epsAbs, mObject.epsRel);
    geod->setInitialPosition(mObject.startPos);
    geod->setStepSizeControlled(mObject.stepsizeControlled);
    geod->setBoundingBox(
        m4d::vec4(-DBL_MAX, -DBL_MAX, -DBL_MAX, -DBL_MAX), m4d::vec4(DBL_MAX, DBL_MAX, DBL_MAX, DBL_MAX));

    /* -----------------------------------------
     *    Set local initial direction.
     * ----------------------------------------- */
    double eta = 1.0;
    double y0 = 1.0;

    if (mObject.type == m4d::enum_geodesic_timelike && fabs(mObject.vel) < 1.0) {
        y0 = 1.0 / sqrt(1.0 - mObject.vel * mObject.vel);
        eta = mObject.vel * y0;
    }

    double chi = mObject.chi * DEG_TO_RAD;
    double ksi = mObject.ksi * DEG_TO_RAD;
    m4d::vec4 locDir = M4D_SIGN(mObject.timeDirection) * y0 * mObject.base[0]
        + eta * sin(chi) * cos(ksi) * mObject.base[1] + eta * sin(chi) * sin(ksi) * mObject.base[2]
        + eta * cos(chi) * mObject.base[3];
    m4d::vec4 coDir;
    mObject.currMetric->localToCoord(mObject.startPos, locDir, coDir, mObject.tetradType);
    geod->setInitialDirection(coDir);

    /* ---------------------------------
     *   Transform local tetrad to
     *   coordinate representation.
     * --------------------------------- */
    m4d::vec4 e0, e1, e2, e3;
    mObject.currMetric->localToCoord(mObject.startPos, mObject.base[0], e0, mObject.tetradType);
    mObject.currMetric->localToCoord(mObject.startPos, mObject.base[1], e1, mObject.tetradType);
    mObject.currMetric->localToCoord(mObject.startPos, mObject.base[2], e2, mObject.tetradType);
    mObject.currMetric->localToCoord(mObject.startPos, mObject.base[3], e3, mObject.tetradType);

    /* ---------------------------------
     *   Calculate parallel transport.
     * --------------------------------- */
    std::vector<m4d::vec4> base[4];
    std::cerr << m4d::stl_break_condition[geod->calcParTransport(mObject.startPos, coDir, e0, e1, e2, e3,
        mObject.maxNumPoints, mObject.points, mObject.dirs, mObject.lambda, base[0], base[1], base[2], base[3])]
              << std::endl;
    std::cerr << "#points: " << mObject.points.size() << std::endl;
    if (mObject.points.size() <= 0)
        return -1;

    unsigned int i = 0;
    m4d::vec4 newBase;
    // double prod;
    while (i < mObject.points.size()) {
        fprintf(stdout, "%5d  %12.8f %12.8f %12.8f %12.8f  ", i, mObject.points[i][0], mObject.points[i][1],
            mObject.points[i][2], mObject.points[i][3]);
        for (int j = 0; j < 4; j++) {
#if 1
            fprintf(
                stdout, "%12.8f %12.8f %12.8f %12.8f  ", base[j][i][0], base[j][i][1], base[j][i][2], base[j][i][3]);
#else
            mObject.currMetric->coordToLocal(mObject.points[i], base[j][i], newBase, mObject.tetradType);
            fprintf(stdout, "%12.8f %12.8f %12.8f %12.8f  ", newBase[0], newBase[1], newBase[2], newBase[3]);
#endif
        }
        fprintf(stdout, "\n");
        i++;
    }
    std::cerr << std::endl;

    delete geod;
    return 1;
}
