/**
 * @file    m4dMetricErezRosenVar.h
 * @author  Thomas Mueller
 *
 * This file is part of the m4d-library.
 */
#ifndef M4D_METRIC_EREZ_ROSEN_VAR_H
#define M4D_METRIC_EREZ_ROSEN_VAR_H

#include "m4dMetric.h"

namespace m4d {

class MetricErezRosenVar : public Metric
{
public:
    MetricErezRosenVar(double mass = 1.0, double q = 0.1);
    virtual ~MetricErezRosenVar();

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

    virtual double testConstraint(const double y[], const double kappa);

    virtual bool setParam(const char* pName, double val);

    virtual double getCircularVelocity(const double r, const enum_nat_tetrad_type tedType = enum_nat_tetrad_default);
    virtual vec4 getCircularFourVel(const vec4 pos, const enum_nat_tetrad_type tedType = enum_nat_tetrad_default);

    virtual bool report(const vec4 pos, const vec4 cdir, char*& text);

    // --------- protected methods -----------
protected:
    virtual void setStandardValues();

    bool calcPotentials(const double* pos, double& psi, double& g, double& delta);
    bool calcDiffPots(const double* pos, double& psi, double& g, double& delta, double& dpsidr, double& dpsidtheta,
        double& dgdr, double& dgdtheta, double& dDdr, double& dDdtheta);

    // -------- protected attribute ---------
protected:
    double mMass;
    double mQ;
};

} // end namespace m4d

#endif
