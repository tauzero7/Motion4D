/**
 * @file    m4dMetricPravda_C.h
 * @author  Felix Beslmeisl
 *
 * @brief  This is the C-Metric as given in class MetricPTD_C in coordinates eliminating the linear term in the polynom
 \f$u^3 +au +b\f$.

     The line element is given by

     \f[ds^2 = \frac{1}{A(x+y)^2}\left( \frac1{f(x)} dx^2 + f(x) dp^2 - \frac1{f(-y)}dy^2 +f(-y)dq^2 \right) \f]
     with \f$f(u) := \pm (-2mAu^3 -u^2 + 1)\f$.

     This metric ist dicussed in<br>
     Pravda, V. and Pravdov&aacute; ,A., <b>Co-accelerated particles in the C-metric</b>  Class. Quantum Gravit.
 <b>18</b>, 1205 (2001).

 * This file is part of the m4d-library.
 */
#ifndef M4DMETRICPRAVDA_C_H
#define M4DMETRICPRAVDA_C_H

#include "m4dMetric.h"

namespace m4d {

class MetricPravda_C : public Metric
{
public:
    MetricPravda_C(double A = 0.1, double m = 1.0);
    virtual ~MetricPravda_C();

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
    double Par_A;
    double Par_m;
};

} // end namespace m4d

#endif // M4DMETRICPTDC_H
