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
#endif

namespace m4d {

class API_M4D_EXPORT IntegratorDatabase {
public:
    IntegratorDatabase();
    ~IntegratorDatabase();
    
    int              getNumIntegrators();
    Geodesic*        getIntegrator(Metric* cMetric, enum_integrator  num);
    Geodesic*        getIntegrator(Metric* cMetric, const char* mName);
    const char*      getIntegratorName(enum_integrator  num);
    enum_integrator  getIntegratorNr(const char* mName);

    void         printIntegratorList(FILE* fptr = stderr);

protected:
    void        init();
    Geodesic*   initializeIntegrator(Metric* cMetric, enum_integrator  num);

private:
    std::map<std::string, enum_integrator>            mIntegratorMap;
    std::map<std::string, enum_integrator>            mIntegratorNickMap;
    std::map<std::string, enum_integrator>::iterator  mIntegratorMapItr;
};

} // end namespace m4d

#endif


