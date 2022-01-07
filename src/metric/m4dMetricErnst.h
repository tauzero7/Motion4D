/**
 * @file    m4dMetricErnst.h
 * @author  Thomas Mueller
 *
 * @brief  Ernst metric in spherical Schwarzschild-like coordinates (t,r,theta,phi).

     The line element is given by

     \f[ds^2 = \Lambda^2\left[-\left(1-\frac{2m}{r}\right) dt^2 + \frac{dr^2}{1-2m/r} + r^2d\vartheta^2\right] +
\frac{r^2\sin^2\vartheta}{\Lambda^2}d\varphi^2,\f] where \f$\Lambda = 1+B^2r^2\sin^2\vartheta\f$.

     The natural local tetrad reads:
     \f[ \mathbf{e}_{(t)} = \frac{1}{\Lambda\sqrt{1-2m/r}}\partial_t,\quad
\mathbf{e}_{(r)}=\frac{\sqrt{1-2m/r}}{\Lambda}\partial_r,\quad \mathbf{e}_{(\vartheta)}=\frac{1}{\Lambda
r}\partial_{\vartheta},\quad \mathbf{e}_{(\varphi)}=\frac{\Lambda}{r\sin\vartheta}\partial_{\varphi}.\f]

     Detailed discussions about the Ernst metric can be found in
     <ul>
       <li> Frederick J. Ernst, "Black holes in a magnetic universe," J. Math. Phys. <b>17</b>, 54--56 (1976).</li>
       <li> R.A. Konoplya, "Magnetised black hole as a gravitational lens," Phys. Lett. B <b>644</b>, 219--223
(2007).</li>
     </ul>

* This file is part of the m4d-library.
*/

#ifndef M4D_METRIC_ERNST_H
#define M4D_METRIC_ERNST_H

#include "m4dMetric.h"

namespace m4d {

// ---------------------------------------------------
//    class definition:   MetricErnst
// ---------------------------------------------------
class MetricErnst : public Metric
{
public:
    MetricErnst(double mass = 1.0, double B = 0.1);
    virtual ~MetricErnst();

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

    virtual bool effPotentialValue(
        const vec4 pos, const vec4 cdir, enum_geodesic_type type, const double x, double& val);
    virtual bool totEnergy(const vec4 pos, const vec4 cdir, const double x, double& val);

    virtual void calcFmu_nu(const double* pos);

    virtual bool report(const vec4 pos, const vec4 cdir, char*& text);

    // --------- protected methods -----------
protected:
    virtual void setStandardValues();

    // -------- protected attribute ---------
protected:
    double mMass;
    double mB;
};

} // end namespace m4d

#endif
