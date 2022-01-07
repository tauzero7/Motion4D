/**
 * @file    m4dMetricBertottiKasner.h
 * @author  Thomas Mueller
 *
 * @brief  Bertotti-Kasner metric in spherical coordinates (t,r,theta,phi).

     The line element is given by

     \f[ds^2 = \Lambda^2\left[\left(1-\frac{2m}{r}\right) dt^2 - \frac{dr^2}{1-2m/r} - r^2d\vartheta^2\right] -
 \frac{r^2\sin^2\vartheta}{\Lambda^2}d\varphi^2,\f] where \f$\Lambda = 1+B^2r^2\sin^2\vartheta\f$.

     The natural local tetrad reads:
     \f[ \mathbf{e}_{(t)} = \frac{1}{\Lambda\sqrt{1-2m/r}}\partial_t,\quad
 \mathbf{e}_{(r)}=\frac{\sqrt{1-2m/r}}{\Lambda}\partial_r,\quad \mathbf{e}_{(\vartheta)}=\frac{1}{\Lambda
 r}\partial_{\vartheta},\quad \mathbf{e}_{(\varphi)}=\frac{\Lambda}{r\sin\vartheta}\partial_{\varphi}.\f]

     Detailed discussions about the BertottiKasner metric can be found in

 * This file is part of the m4d-library.
 */
#ifndef M4D_METRIC_BERTOTTI_KASNER_H
#define M4D_METRIC_BERTOTTI_KASNER_H

#include "m4dMetric.h"

namespace m4d {

class MetricBertottiKasner : public Metric
{
public:
    MetricBertottiKasner(double lambda = 1.0);
    virtual ~MetricBertottiKasner();

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
    double mLambda;
};

} // end namespace m4d

#endif
