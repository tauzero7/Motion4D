/**
 * @file    m4dMetricPTD_BI.h
 * @author  Felix Beslmeisl
 *
 * @brief  \ref lit_stephani "Exact Solutions:" Petrov type D solutions - Case BI.

     The line element is given by

     \f[ds^2 = r^2\left( d\vartheta^2 - \sin^2 \vartheta dt^2 \right) + \frac{r}{r-b}dr^2+\frac{r-b}{r}d\varphi^2. \f]

     The natural local tetrad is given by
     \f[  \mathbf{e}_{(t)} = \frac {1}{r\sin{\vartheta}}\partial_t, \quad
          \mathbf{e}_{(r)} = \sqrt{{\frac {r-b}{r}}}\partial_r, \quad
          \mathbf{e}_{(\vartheta)} = \frac{1}{r}\partial_{\vartheta}, \quad
          \mathbf{e}_{(\varphi)} = \sqrt{{\frac {r}{r-b}}}\partial_{\varphi}.\f]

     From the Euler-Lagrange formalism, we have an effective Potential
 \f$\frac{1}{2}\dot{r}^2+\frac{1}{2}V_{\mbox{eff}}(r)=\frac{1}{2}C_0^2\f$

      \f[V_{\mbox{eff}}(r)=+C_0^2+K\frac{-b+r}{r^3}-\kappa c^2\frac{-b+r}{r}\f]

     with the following constants of motion:

     \f[  C_0^2 = \dot{\varphi}^2 \frac{r-b}{r}  \f]
     \f[  K     = \dot{t}^2 r^4 - \dot{\vartheta}^2 r^4 \sin^2{\vartheta}\f]
     \f[  -\kappa c^2 = -K\frac1{r^2}-\dot{r}^2 \frac{r}{r - b} -\dot{\varphi}^2\frac{r - b}{r},\f]

     but this spacetime is not spherical symmetric, so additionally particles out of the
 \f$\vartheta=\frac{\pi}{2}\f$-plane fall into one of the poles.


 * This file is part of the m4d-library.
 */
#ifndef M4DMETRICPTDBI_H
#define M4DMETRICPTDBI_H

#include "m4dMetric.h"

namespace m4d {

class MetricPTD_BI : public Metric
{
public:
    MetricPTD_BI(double b = 1.0);
    virtual ~MetricPTD_BI();

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

    virtual bool effPotentialValue(
        const vec4 pos, const vec4 cdir, enum_geodesic_type type, const double x, double& val);
    virtual bool totEnergy(const vec4 pos, const vec4 cdir, const double x, double& val);
    virtual double calculateVeffRoot(double C02, double K, double r0);

    virtual bool setParam(const char* pName, double val);

    virtual bool report(const vec4 pos, const vec4 cdir, char*& text);

    // --------- protected methods -----------
protected:
    virtual void setStandardValues();
    void calcConstantsOfMotion(const vec4 pos, const vec4 cdir);

    // -------- protected attribute ---------
protected:
    double Par_b;
    double K, C0, C2, m0;
};

} // end namespace m4d

#endif // M4DMETRICPTDBI_H
