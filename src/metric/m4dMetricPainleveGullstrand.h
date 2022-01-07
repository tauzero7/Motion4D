/**
 * @file    m4dMetricPainleveGullstrand.h
 * @author  Thomas Mueller
 *
 * @brief  Schwarzschild spacetime Painlev{\'e}-Gullstrand coordinates (T,r,theta,phi).

     The line element is given by
         \f[ds^2 = -c^2 dT^2 + \left(dr+\sqrt{\frac{r_s}{r}}c\,dT\right)^2 + r^2 \left(d\vartheta^2 + \sin(\vartheta)^2
 d\varphi^2\right)\f] where \f$r_s=2GM/c^2\f$ is the Schwarzschild radius. G is Newton's constant, M is the mass of the
 black hole, and c is the speed of light.

         There are two natural local tetrads:
              - Static tetrad:  (enum_nat_tetrad_static)
              \f[ \mathbf{e}_{(0)}^s = \frac{1}{c\,\sqrt{1-r_s/r}}\partial_T,\quad \mathbf{e}_{(1)}^s =
 \frac{1}{c\,\sqrt{r/r_s-1}}\partial_T+\sqrt{1-\frac{r_s}{r}}\partial_r,\quad \mathbf{e}_{(2)}^s =
 \frac{1}{r}\partial_{\vartheta},\quad \mathbf{e}_{(3)}^s = \frac{1}{r\sin\vartheta}\partial_{\varphi}. \f]
              - Freely falling tetrad: (enum_nat_tetrad_freefall)
              \f[ \mathbf{e}_{(0)}^f = \frac{1}{c}\partial_T-\sqrt{\frac{r_s}{r}}\partial_r,\quad \mathbf{e}_{(1)}^f =
 \partial_r,\quad \mathbf{e}_{(2)}^f = \frac{1}{r}\partial_{\vartheta},\quad \mathbf{e}_{(3)}^f =
 \frac{1}{r\sin\vartheta}\partial_{\varphi}.\f]

         Detailed discussions about the Painlev{\'e}-Gullstrand metric can be found
         in ...

 * This file is part of the m4d-library.
 */
#ifndef M4D_METRIC_PAINLEVE_GULLSTRAND_H
#define M4D_METRIC_PAINLEVE_GULLSTRAND_H

#include "m4dMetric.h"

namespace m4d {

class MetricPainleveGullstrand : public Metric
{
public:
    MetricPainleveGullstrand(double mass = 1.0);
    virtual ~MetricPainleveGullstrand();

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

    virtual bool report(const vec4 pos, const vec4 cdir, char*& text);

    // --------- protected methods -----------
protected:
    virtual void setStandardValues();

    // -------- protected attribute ---------
protected:
    double rs;
    double mMass;
};

} // end namespace m4d

#endif
