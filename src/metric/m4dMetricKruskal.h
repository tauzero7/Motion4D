/**
 * @file    m4dMetricKruskal.h
 * @author  Thomas Mueller
 *
 * @brief  Kruskal metric in spherical coordinates (t,r,theta,phi).

             The line element is given by
 *
 *  This file is part of libMotion4D.
 */
#ifndef M4D_METRIC_KRUSKAL_H
#define M4D_METRIC_KRUSKAL_H

#include "m4dMetric.h"

namespace m4d {

class MetricKruskal : public Metric
{
public:
    MetricKruskal(double mass = 1.0);
    virtual ~MetricKruskal();

public:
    virtual bool calculateMetric(const double* pos);
    virtual bool calculateChristoffels(const double* pos);

    virtual void localToCoord(
        const double* pos, const double* ldir, double* dir, enum_nat_tetrad_type type = enum_nat_tetrad_default);
    virtual void coordToLocal(
        const double* pos, const double* cdir, double* ldir, enum_nat_tetrad_type type = enum_nat_tetrad_default);

    virtual bool breakCondition(const double* pos);
    virtual bool calcDerivs(const double y[], double dydx[]);

    virtual double testConstraint(const double y[], const double kappa);

    virtual bool setParam(const char* pName, double val);

    virtual bool report(const vec4 pos, const vec4 cdir, char*& text);

    virtual int transToPseudoCart(vec4 p, vec4& cp);

protected:
    virtual void setStandardValues();

    double get_r(double T, double X);

protected:
    double mMass;
    double rs;
};

} // end namespace m4d

#endif
