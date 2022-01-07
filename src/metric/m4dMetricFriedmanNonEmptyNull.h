/**
 * @file    m4dMetricFriedmanNonEmptyNull.h
 * @author  Thomas Mueller
 *
 * @brief  The three non-empty Friedman models with vanishing cosmological constant.

     The line element is given in conformal spherical coordinates (chi,r,theta,phi),

     \f[ds^2 =
 R(\chi)^2\left(-d\chi^2+\frac{dr^2+r^2(d\vartheta^2+\sin^2\vartheta\,d\varphi^2)}{\left(1+\frac{k}{4}r^2\right)^2}\right),\f]

     where k=0,+1,-1.

     The metric is taken from Rindler.

 * This file is part of the m4d-library.
 */
#ifndef M4D_METRIC_FRIEDMAN_NON_EMPTY_NULL_H
#define M4D_METRIC_FRIEDMAN_NON_EMPTY_NULL_H

#include "m4dMetric.h"

namespace m4d {

// ---------------------------------------------------
//    class definition:   MetricFriedmanNonEmptyNull
// ---------------------------------------------------
class MetricFriedmanNonEmptyNull : public Metric
{
public:
    MetricFriedmanNonEmptyNull(double mass = 1.0, double k = 0.0);
    virtual ~MetricFriedmanNonEmptyNull();

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

    virtual bool transToTwoPlusOne(vec4 p, vec4& cp);

    virtual bool report(const vec4 pos, const vec4 cdir, char*& text);

    // --------- protected methods -----------
protected:
    virtual void setStandardValues();
    double calc_R(const double chi);
    double calc_dR(const double chi);
    double calc_ddR(const double chi);

    // -------- protected attribute ---------
protected:
    double mC;
    int mK;
};

} // end namespace m4d

#endif
