/**
 * @file    m4dMetricMinkRotLattice.h
 * @author  Thomas Mueller
 *
 * @brief  Minkowski spacetime in rotating cylindrical coordinates (t,r,phi,z).

     The line element is given by

    \f[ ds^2 = -\left(1-\frac{r^2\omega^2}{c^2}\right) c^2 dt^2 + 2r^2 \omega dt\, d\varphi
           + r^2 d\varphi^2 + dr^2 + dz^2 \f].


     The local tetrad of the comoving observer is
  \f[ \mathbf{e}_{(t)} = \frac{1}{c}\partial_t - \frac{\omega}{c}\partial_{\varphi},\quad
 \mathbf{e}_{(r)}=\partial_r,\quad \mathbf{e}_{(z)}=\partial_z,\quad \mathbf{e}_{(\varphi)} =
 \frac{1}{r}\partial_{\varphi}, \f]

     whereas the static observer has the local tetrad
  \f[ \mathbf{e}_{(t)} = \frac{1}{c\,\sqrt{1-\omega^2r^2/c^2}}\partial_t,\quad \mathbf{e}_{(r)}=\partial_r,\quad
 \mathbf{e}_{(z)}=\partial_z,\quad \mathbf{e}_{(\varphi)} = \frac{\omega r}{c^2\sqrt{1-\omega^2r^2/c^2}}\partial_t +
 \frac{\sqrt{1-\omega^2r^2/c^2}}{r}\partial_{\varphi}.\f]

   natLocTetrad:  default = static

 * This file is part of the m4d-library.
 */
#ifndef M4D_METRIC_MINK_ROTLATTICE_H
#define M4D_METRIC_MINK_ROTLATTICE_H

#include "m4dMetric.h"

namespace m4d {

class MetricMinkRotLattice : public Metric
{
public:
    MetricMinkRotLattice(double omega = 0.0);
    virtual ~MetricMinkRotLattice();

    // --------- public methods -----------
public:
    virtual bool calculateMetric(const double* pos);
    virtual bool calculateChristoffels(const double* pos);

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
    double mOmega;
};

} // end namespace m4d

#endif
