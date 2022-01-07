/**
 * @file    m4dMetricRotDihole.h
 * @author  Thomas Mueller
 *
 * @brief  Rotating extreme Reissner-Nordstrom dihole metric in spherical coordinates (t,r,theta,phi).
 *
 *  This file is part of libMotion4D.
 */
#ifndef M4D_METRIC_ROT_DIHOLE_H
#define M4D_METRIC_ROT_DIHOLE_H

#include "m4dMetric.h"
#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>

namespace m4d {

class MetricRotDihole : public Metric
{
public:
    MetricRotDihole(double mass1 = 1.0, double mass2 = 1.0, double omega = 0.0);
    virtual ~MetricRotDihole();

    // --------- public methods -----------
public:
    virtual bool calculateMetric(const double* pos);
    virtual bool calculateChristoffels(const double* pos);
    virtual bool calculateChrisD(const double* pos);

    virtual void localToCoord(
        const double* pos, const double* ldir, double* dir, enum_nat_tetrad_type type = enum_nat_tetrad_default);
    virtual void coordToLocal(
        const double* pos, const double* cdir, double* ldir, enum_nat_tetrad_type type = enum_nat_tetrad_default);

    virtual bool breakCondition(const double* pos);

    virtual bool setParam(const char* pName, double val);

    virtual bool report(const vec4 pos, const vec4 cdir, char*& text);

    // --------- protected methods -----------
protected:
    virtual void setStandardValues();
    void calc_r(const double* pos, double& r1, double& r2);
    double calc_U(const double* pos);
    void calc_dU(const double* pos, double& dUdt, double& dUdx, double& dUdy, double& dUdz);

    // -------- protected attribute ---------
protected:
    double mMass1;
    double mMass2;
    double mOmega;
};

} // end namespace m4d

#endif
