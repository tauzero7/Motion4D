// -------------------------------------------------------------------------------
/*
   m4dMetricDatabase.h

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

/*!  \class  m4d::MetricDatabase
 *   \brief  database for all metrics
 *
 *           The getMetric() method generates a new instance of the metric.
 *           The user has to delete this metric himself!
 *
 */
// -------------------------------------------------------------------------------
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

class API_M4D_EXPORT MetricDatabase {
public:
    MetricDatabase();
    ~MetricDatabase();

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
