/**
 * @file    testCPmotion.cpp
 * @author  Thomas Mueller
 *
 * This file is part of the m4d-library.
 */

#include <cmath>
#include <iostream>

#include <extra/m4dObject.h>
#include <m4dGlobalDefs.h>
#include <math/TransCoordinates.h>
#include <metric/m4dMetricDatabase.h>
#include <motion/m4dMotionChargedParticle.h>

/* -----------------------------------------------------------------
 *                     m a i n
 * ----------------------------------------------------------------- */
int main(int argc, char* argv[])
{
    if (argc != 1) {
        fprintf(stderr, "Usage: ./m4dTestCPmotion\n");
        return -1;
    }

    m4d::Object mObject;

    /* ----------------------------------
     *  initialize metric database
     * ---------------------------------- */
    m4d::MetricDatabase metricDB;

    m4d::Metric* metric = metricDB.getMetric("Ernst");
    metric->setParam("mass", 1.0);
    metric->setParam("b", 0.1);
    metric->printF();

    /* ----------------------------------
     *  Define initial position.
     * ---------------------------------- */
    m4d::vec4 initPos = m4d::vec4(0.0, 4.0, M_PI_2, 0.0);

    /* ----------------------------------
     *  Set FermiWalker integrator.
     * ---------------------------------- */
    m4d::MotionChargedParticle* motion = new m4d::MotionChargedParticle(metric);
    motion->setAffineParamStep(0.01);

    /* ----------------------------------
     *  Set initial position, direction, and local tetrad.
     * ---------------------------------- */
    motion->setInitialPosition(initPos);
    motion->setInitialVelocity(1, 0.5, 0.0, M_PI_2, 0.1);

    /* ----------------------------------
     *  Calculate one orbit.
     * ---------------------------------- */
    int maxNumPoints = 2000;
    m4d::vec4 pos, ne;
    int i = 0;
    while (i < maxNumPoints) {
        pos = motion->getPosition();
        fprintf(stderr, "\r%5d  %8.4f", i, motion->getAffineParam());
        fprintf(stdout, "%5d  %8.4f  %12.8f %12.8f %12.8f %12.8f \n", i, motion->getAffineParam(), pos[0], pos[1],
            pos[2], pos[3]);
        motion->nextStep();
        i++;
    }
    fprintf(stderr, "\n");

    delete motion;
    delete metric;
    return 1;
}
