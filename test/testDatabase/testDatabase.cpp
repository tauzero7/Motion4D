/**
 * @file    testDatabase.cpp
 * @author  Thomas Mueller
 *
 * @brief Print all implemented metrics and integrators.

     Usage:  \verbatim ./m4dTestDatabase \endverbatim

 * This file is part of the m4d-library.
 */
#include <cmath>
#include <iostream>

#include <m4dGlobalDefs.h>
#include <metric/m4dMetricDatabase.h>
#include <motion/m4dMotionList.h>

/* -----------------------------------------------------------------
 *                     m a i n
 * ----------------------------------------------------------------- */
int main(int argc, char* argv[])
{
    fprintf(stderr, "\n\n=================== Test Database ===================\n\n");

    m4d::MetricDatabase metricDB;
    metricDB.printMetricList(stdout);

    fprintf(stdout, "\nnum   Integrator name\n-----------------------------------------------------------\n");
    for (unsigned int i = 0; i < m4d::NUM_GEOD_SOLVERS; i++) {
        fprintf(stdout, "%2d   %s\n", i, m4d::stl_solver_names[i]);
    }

    fprintf(stdout, "\nnum   name of tetrad type\n-----------------------------------------------------------\n");
    for (unsigned int i = 0; i < m4d::NUM_ENUM_NAT_TETRAD_TYPES; i++) {
        fprintf(stdout, "%2d   %s\n", i, m4d::stl_nat_tetrad_types[i]);
    }

    // for (int i = 1; i <= m4d::MetricList::NUM_METRICS; i++) {
    //     m4d::Metric* metric = metricDB.getMetric(static_cast<m4d::MetricList::enum_metric>(i));
    //     std::cerr << i << " " << metric->getMetricName() << std::endl;
    // }

    return 1;
}
