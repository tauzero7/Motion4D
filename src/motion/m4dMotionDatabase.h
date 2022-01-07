/**
 * @file    m4dMotionDatabase.h
 * @author  Thomas Mueller
 *
 * @brief   database for all integrators
 *
 * This file is part of the m4d-library.
 *
 */
#ifndef M4D_MOTION_DATABASE_H
#define M4D_MOTION_DATABASE_H

#include "m4dMotionList.h"
#include <iostream>
#include <map>

#ifdef _WIN32
#ifndef __GNUC__
#pragma warning(disable : 4244)
#endif
#endif

namespace m4d {

class API_M4D_EXPORT IntegratorDatabase
{
public:
    IntegratorDatabase();
    ~IntegratorDatabase();

    int getNumIntegrators();
    Geodesic* getIntegrator(Metric* cMetric, enum_integrator num);
    Geodesic* getIntegrator(Metric* cMetric, const char* mName);
    const char* getIntegratorName(enum_integrator num);
    enum_integrator getIntegratorNr(const char* mName);

    void printIntegratorList(FILE* fptr = stderr);

protected:
    void init();
    Geodesic* initializeIntegrator(Metric* cMetric, enum_integrator num);

private:
    std::map<std::string, enum_integrator> mIntegratorMap;
    std::map<std::string, enum_integrator> mIntegratorNickMap;
    std::map<std::string, enum_integrator>::iterator mIntegratorMapItr;
};

} // end namespace m4d

#endif
