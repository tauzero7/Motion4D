// -------------------------------------------------------------------------------
/*
   m4dMetricDatabase.cpp

  Copyright (c) 2009-2014-2011  Thomas Mueller, Frank Grave


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
// -------------------------------------------------------------------------------

#include "m4dMotionDatabase.h"

namespace m4d {

//#ifndef _WIN32
IntegratorDatabase* IntegratorDatabase::m_instance = nullptr;
//#endif // _WIN32

/*!
 */
IntegratorDatabase::IntegratorDatabase() {
    init();
}

IntegratorDatabase::~IntegratorDatabase() {
    if (!mIntegratorMap.empty()) {
        mIntegratorMap.clear();
    }
    if (!mIntegratorNickMap.empty()) {
        mIntegratorNickMap.clear();
    }
}


// *********************************** public methods ******************************
/*! Get the number of implemented metrics.
 *
 *  \return int: number of implemented metrics.
 */
int
IntegratorDatabase::getNumIntegrators() {
    return NUM_GEOD_SOLVERS;
}

/*!  Initialize metric 'num' and return the pointer to it.
 *
 *  \param cMetric : pointer to metric.
 *  \param num     : metric number.
 */
Geodesic*
IntegratorDatabase::getIntegrator(Metric* cMetric, enum_integrator num) {
    return initializeIntegrator(cMetric, num);
}

/*! Initialize metric 'mName' and return the pointer to it.
 *
 *  \param cMetric : pointer to metric.
 *  \param mName   : name of metric.
 */
Geodesic*
IntegratorDatabase::getIntegrator(Metric* cMetric, const char* mName) {
    mIntegratorMapItr = mIntegratorNickMap.find(mName);
    if (mIntegratorMapItr == mIntegratorNickMap.end()) {
        mIntegratorMapItr = mIntegratorMap.find(mName);
        if (mIntegratorMapItr == mIntegratorMap.end()) {
            fprintf(stderr, "Integrator '%s' is not implemented!\n", mName);
            return NULL;
        }
    }

    return initializeIntegrator(cMetric, mIntegratorMapItr->second);
}

/*! Get the name of metric 'num'.
 *
 *  \param num : number of metric.
 *  \return string : name of metric.
 */
const char* IntegratorDatabase::getIntegratorName(enum_integrator num) {
    if (num >= 0 && num < (int)NUM_GEOD_SOLVERS) {
        return stl_solver_names[num];
    }
    return nullptr;
}

/*! Get the number of the 'mName' metric.
 *
 *  \param mName : name of metric.
 *  \return enum_metric : number of metric.
 */
enum_integrator IntegratorDatabase::getIntegratorNr(const char* mName) {
    mIntegratorMapItr = mIntegratorNickMap.find(mName);
    if (mIntegratorMapItr == mIntegratorNickMap.end()) {
        mIntegratorMapItr = mIntegratorMap.find(mName);
        if (mIntegratorMapItr == mIntegratorMap.end()) {
            fprintf(stderr, "Integrator '%s' is not implemented!\n", mName);
            return gsUnknown;
        }
    }

    return mIntegratorMapItr->second;
}

/*!  Print list of all metrics.
 *
 *  \param fptr : pointer to file.
 */
void IntegratorDatabase::printIntegratorList(FILE* fptr) {
    fprintf(fptr, "  Solver nicknanme            Solver full name\n");
    fprintf(fptr, "-----------------------------------------------------------\n");
    for (unsigned int i = 0; i < NUM_GEOD_SOLVERS; i++) {
        fprintf(fptr, "%-23s : %s\n", stl_solver_nicknames[i], stl_solver_names[i]);
    }
}


// ******************************** protected methods ******************************

/*!  Initialize database: store the metrics in a map.
 *
 */
void IntegratorDatabase::init() {
    for (unsigned int i = 0; i < NUM_GEOD_SOLVERS; i++) {
        mIntegratorMap.insert(std::pair<std::string, enum_integrator>(stl_solver_names[i], enum_integrator(i)));
        mIntegratorNickMap.insert(std::pair<std::string, enum_integrator>(stl_solver_nicknames[i], enum_integrator(i)));
    }
}


// ---------------------------------------------------
//           Please sort by enum name !!
// ---------------------------------------------------
/*! Initialize metric: generate a new instance.
 *
 *  \param cMetric : pointer to metric.
 *  \param num     : enumeration of metric.
 */
Geodesic*
IntegratorDatabase::initializeIntegrator(Metric* cMetric,  enum_integrator  num) {

    Geodesic*  currGeo = NULL;

    const gsl_odeiv_step_type*  step_type = gsl_odeiv_step_rk4;
    switch (num) {
    case m4d::gsIrk4:
#ifdef USE_DP_INT
    case m4d::gsIdp54:
    case m4d::gsIdp65:
#endif
    case gsInrbs:
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
    default:
        break;
    }

    if (num == gsIrk4) {
        currGeo = new GeodesicRK4(cMetric);
    } else if (num == gsInrbs) {
        currGeo = new GeodesicBS(cMetric);
    }
#ifdef USE_DP_INT
    else if (num == gsIdp54) {
        currGeo = new GeodesicDP54(cMetric);
    } else if (num == gsIdp65) {
        currGeo = new GeodesicDP65(cMetric);
    }
#endif
    else {
        currGeo = new GeodesicGSL(cMetric, step_type, num);
    }

    // mData->geodSolver->setBoundingBox( m4d::vec4(-DBL_MAX,-DBL_MAX,-DBL_MAX,-DBL_MAX), m4d::vec4(DBL_MAX,DBL_MAX,DBL_MAX,DBL_MAX) );
    // mData->geodSolver->setGeodesicType( mData->type );

    return currGeo;
}

} // end namespace m4d
