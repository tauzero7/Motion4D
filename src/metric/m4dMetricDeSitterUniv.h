/**
 * @file    m4dMetricDeSitterUniv.h
 * @author  Thomas Mueller
 *
 * @brief  DeSitter metric in cartesian coordinates (t,x,y,z).

     The line element is given by

     \f[ds^2 = -c^2dt^2+e^{2Ht}\left(dx^2+dy^2+dz^2\right),\f]

     where H is the Hubble parameter.

 * This file is part of the m4d-library.
 */

#ifndef M4D_METRIC_DESITTER_UNIV_H
#define M4D_METRIC_DESITTER_UNIV_H

#include "m4dMetric.h"

namespace m4d {

// ---------------------------------------------------
//    class definition:   MetricDeSitterUniv
// ---------------------------------------------------
class MetricDeSitterUniv : public Metric
{
public:
    MetricDeSitterUniv(double h = 0.1);
    virtual ~MetricDeSitterUniv();

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

    virtual bool transToTwoPlusOne(vec4 p, vec4& cp);

    virtual bool setParam(const char* pName, double val);

    virtual bool report(const vec4 pos, const vec4 cdir, char*& text);

    // --------- protected methods -----------
protected:
    virtual void setStandardValues();

    // -------- protected attribute ---------
protected:
    double mHubble;
};

} // end namespace m4d

#endif
