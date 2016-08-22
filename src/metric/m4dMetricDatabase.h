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

#include <iostream>
#include <map>
#include "m4dMetricList.h"

#ifdef _WIN32
#ifndef __GNUC__
#pragma warning (disable: 4244 )
#endif
#endif

namespace m4d {

class API_EXPORT MetricDatabase {
public:
    static MetricDatabase* M4D_CALL getInstance() {
        static CGuard g;
        if (m_instance == nullptr) {
            m_instance = new MetricDatabase();
        }
        return m_instance;
    }

    int          getNumMetrics();
    Metric*      getMetric(enum_metric  num);
    Metric*      getMetric(const char* mName);
    std::string  getMetricName(enum_metric  num);
	enum_metric  getMetricNr(const char* mName);

    void         printMetricList(FILE* fptr = stderr);

protected:
    void      init();
    Metric*   initializeMetric(enum_metric  num);


    // make the Database class a singleton
private:
    static MetricDatabase* m_instance;
    MetricDatabase();
    MetricDatabase(const MetricDatabase&) {}
    ~MetricDatabase();

    class CGuard {
    public:
        ~CGuard() {
            if (MetricDatabase::m_instance != nullptr) {
                delete MetricDatabase::m_instance;
                MetricDatabase::m_instance = nullptr;
            }
        }
    };

    friend class CGuard;

private:
    std::map<std::string, enum_metric>            mMetricMap;
    std::map<std::string, enum_metric>::iterator  mMetricMapItr;
};

} // end namespace m4d

#endif


