/**
 * @file    m4dMotionDatabase.cpp
 * @author  Thomas Mueller
 *
 * This file is part of the m4d-library.
 *
 */

#include "m4dMotionDatabase.h"

namespace m4d {

IntegratorDatabase::IntegratorDatabase()
{
    init();
}

IntegratorDatabase::~IntegratorDatabase()
{
    if (!mIntegratorMap.empty()) {
        mIntegratorMap.clear();
    }
    if (!mIntegratorNickMap.empty()) {
        mIntegratorNickMap.clear();
    }
}

// *********************************** public methods ******************************

int IntegratorDatabase::getNumIntegrators()
{
    return NUM_GEOD_SOLVERS;
}

Geodesic* IntegratorDatabase::getIntegrator(Metric* cMetric, enum_integrator num)
{
    return initializeIntegrator(cMetric, num);
}

Geodesic* IntegratorDatabase::getIntegrator(Metric* cMetric, const char* mName)
{
    mIntegratorMapItr = mIntegratorNickMap.find(mName);
    if (mIntegratorMapItr == mIntegratorNickMap.end()) {
        mIntegratorMapItr = mIntegratorMap.find(mName);
        if (mIntegratorMapItr == mIntegratorMap.end()) {
            fprintf(stderr, "Integrator '%s' is not implemented!\n", mName);
            return nullptr;
        }
    }

    return initializeIntegrator(cMetric, mIntegratorMapItr->second);
}

const char* IntegratorDatabase::getIntegratorName(enum_integrator num)
{
    if (num >= 0 && num < (int)NUM_GEOD_SOLVERS) {
        return stl_solver_names[num];
    }
    return nullptr;
}

enum_integrator IntegratorDatabase::getIntegratorNr(const char* mName)
{
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

void IntegratorDatabase::printIntegratorList(FILE* fptr)
{
    fprintf(fptr, "  Solver nicknanme            Solver full name\n");
    fprintf(fptr, "-----------------------------------------------------------\n");
    for (unsigned int i = 0; i < NUM_GEOD_SOLVERS; i++) {
        fprintf(fptr, "%-23s : %s\n", stl_solver_nicknames[i], stl_solver_names[i]);
    }
}

// ******************************** protected methods ******************************

void IntegratorDatabase::init()
{
    for (unsigned int i = 0; i < NUM_GEOD_SOLVERS; i++) {
        mIntegratorMap.insert(std::pair<std::string, enum_integrator>(stl_solver_names[i], enum_integrator(i)));
        mIntegratorNickMap.insert(std::pair<std::string, enum_integrator>(stl_solver_nicknames[i], enum_integrator(i)));
    }
}

// ---------------------------------------------------
//           Please sort by enum name !!
// ---------------------------------------------------
Geodesic* IntegratorDatabase::initializeIntegrator(Metric* cMetric, enum_integrator num)
{

    Geodesic* currGeo = nullptr;

    const gsl_odeiv_step_type* step_type = gsl_odeiv_step_rk4;
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
    }
    else if (num == gsInrbs) {
        currGeo = new GeodesicBS(cMetric);
    }
#ifdef USE_DP_INT
    else if (num == gsIdp54) {
        currGeo = new GeodesicDP54(cMetric);
    }
    else if (num == gsIdp65) {
        currGeo = new GeodesicDP65(cMetric);
    }
#endif
    else {
        currGeo = new GeodesicGSL(cMetric, step_type, num);
    }

    // mData->geodSolver->setBoundingBox( m4d::vec4(-DBL_MAX,-DBL_MAX,-DBL_MAX,-DBL_MAX),
    // m4d::vec4(DBL_MAX,DBL_MAX,DBL_MAX,DBL_MAX) ); mData->geodSolver->setGeodesicType( mData->type );

    return currGeo;
}

} // end namespace m4d
