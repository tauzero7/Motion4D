/**
 * @file    m4dMetricPTD_C.h
 * @author  Felix Beslmeisl
 *
 * @brief  \ref lit_stephani "Exact Solutions:" Petrov type D solutions - Case C.

     The line element is given by

     \f[ds^2 = \frac{1}{(x+y)^2}\left( \frac1{f(x)} dx^2 + f(x) d\varphi^2 - \frac1{f(-y)}dy^2 +f(-y)dt^2 \right) \f]
     with \f$f(u) := \pm (u^3 +au +b)\f$.

     The natural local tetrad is given by
     \f[  \mathbf{e}_{(t)} = (x+y)\frac{1}{\sqrt{-y^3-ay+b}}\,\partial_t, \quad
          \mathbf{e}_{(x)} = (x+y)\sqrt{x^3+ax+b}\,\partial_x, \quad
          \mathbf{e}_{(y)} = (x+y)\sqrt{-y^3-ay+b}\,\partial_y, \quad
          \mathbf{e}_{(\varphi)} = (x+y)\frac{1}{\sqrt{x^3+ax+b}}\,\partial_{\varphi}.\f]

* This file is part of the m4d-library.
*/
#ifndef M4DMETRICPTDC_H
#define M4DMETRICPTDC_H

#include "m4dMetric.h"

namespace m4d {

class MetricPTD_C : public Metric
{
public:
    MetricPTD_C(double a = -1.0, double b = 1.0);
    virtual ~MetricPTD_C();

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
    virtual int transToPseudoCart(vec4 p, vec4& cp);

    virtual double testConstraint(const double y[], const double kappa);

    virtual bool setParam(const char* pName, double val);

    virtual bool report(const vec4 pos, const vec4 cdir, char*& text);

    // --------- specific public methods ----------
public:
    virtual void calculateRoots(vec3& roots, double p, double q);

    // --------- protected methods -----------
protected:
    virtual void setStandardValues();

    // -------- protected attribute ---------
protected:
    double Par_a;
    double Par_b;
};

} // end namespace m4d

#endif // M4DMETRICPTDC_H
