/**
 * @file    m4dMetricDatabase.h
 * @author  Thomas Mueller
 *
 *   @brief  database for all metrics
 *
 *    The getMetric() method generates a new instance of the metric.
 *    The user has to delete this metric himself!
 *
 *  This file is part of libMotion4D.
 */
#ifndef M4D_METRIC_DATABASE_H
#define M4D_METRIC_DATABASE_H

#include "m4dMetricList.h"
#include <iostream>
#include <map>

#ifdef _WIN32
#ifndef __GNUC__
#pragma warning(disable : 4244)
#endif
#endif

namespace m4d {

class API_M4D_EXPORT MetricDatabase
{
public:
    MetricDatabase();
    ~MetricDatabase();

    bool allMetricsAvailable();

    /**
     * @brief Get number of implementred metrics.
     * @return
     */
    int getNumMetrics();

    /**
     * @brief Initialize metric and get pointer to it.
     * @param num
     * @return
     */
    Metric* getMetric(MetricList::enum_metric num);

    /**
     * @brief Initialize metric and get pointer to it.
     * @param mName
     * @return
     */
    Metric* getMetric(const char* mName);

    /**
     * @brief Get name of metric by id.
     * @param num
     * @return
     */
    const char* getMetricName(MetricList::enum_metric num);

    /**
     * @brief Get metric id by name.
     * @param mName
     * @return
     */
    MetricList::enum_metric getMetricNr(const char* mName);

    /**
     * @brief Print list of all available metrics
     * @param fptr
     */
    void printMetricList(FILE* fptr = stderr);

protected:
    void init();
    Metric* initializeMetric(MetricList::enum_metric num);

private:
    std::map<std::string, MetricList::enum_metric> mMetricMap;
    std::map<std::string, MetricList::enum_metric>::iterator mMetricMapItr;
};

} // end namespace m4d

#endif
