// -------------------------------------------------------------------------------
/*
   m4dMotionDatabase.h

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

/*!  \class  m4d::IntegratorDatabase
 *   \brief  database for all integrators
 *
 *
 */
// -------------------------------------------------------------------------------
#ifndef M4D_MOTION_DATABASE_H
#define M4D_MOTION_DATABASE_H

#include <iostream>
#include <map>
#include "m4dMotionList.h"

#ifdef _WIN32
#ifndef __GNUC__
#pragma warning (disable: 4244 )
#endif
#else
#define MOTION_API
#endif

namespace m4d {

class EXTRA_API IntegratorDatabase {
public:
    static IntegratorDatabase* getInstance() {
        static CGuard g;
        if (m_instance == nullptr) {
            m_instance = new IntegratorDatabase();
        }
        return m_instance;
    }

    int              getNumIntegrators();
    Geodesic*        getIntegrator(Metric* cMetric, enum_integrator  num);
    Geodesic*        getIntegrator(Metric* cMetric, std::string mName);
    std::string      getIntegratorName(enum_integrator  num);
    enum_integrator  getIntegratorNr(std::string mName);

    void         printIntegratorList(FILE* fptr = stderr);

protected:
    void        init();
    Geodesic*   initializeIntegrator(Metric* cMetric, enum_integrator  num);

private:
    static IntegratorDatabase* m_instance;
    IntegratorDatabase();
    IntegratorDatabase(const IntegratorDatabase &) {}
    ~IntegratorDatabase();

    class CGuard {
    public:
        ~CGuard() {
            if (IntegratorDatabase::m_instance != nullptr) {
                delete IntegratorDatabase::m_instance;
                IntegratorDatabase::m_instance = nullptr;
            }
        }
    };

    friend class CGuard;


private:
    std::map<std::string, enum_integrator>            mIntegratorMap;
    std::map<std::string, enum_integrator>            mIntegratorNickMap;
    std::map<std::string, enum_integrator>::iterator  mIntegratorMapItr;
};

} // end namespace m4d

#endif


