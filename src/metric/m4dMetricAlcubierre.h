/**
 * @file    m4dMetricAlcubierre.h
 * @author  Thomas Mueller
 *
 * @brief  Alcubierre warp metric in Cartesian coordinates (t,x,y,z).

     The line element is given by

    \f[ ds^2 = -c^2 dt^2 + \left[dx-v_s(t)f(r_s(t))\right]^2+dy^2+dz^2\f]

     with \f$v_s=\frac{dx_s(t)}{dt}\f$, \f$r_s(t)=\sqrt{(x-x_s(t)^2+y^2+z^2}\f$, and

    \f[ f(r_s)=\frac{\tanh(\sigma(r_s+R))-\tanh(\sigma(r_s-R))}{2\tanh(\sigma R)}. \f]

     The natural comoving local tetrad  is given by
     \f[ \mathbf{e}_{(0)} = \frac{1}{c}\partial_t+\frac{v_s f(r_s(t))}{c}\partial_x,\quad
 \mathbf{e}_{(1)}=\partial_x,\quad \mathbf{e}_{(2)}=\partial_y,\quad \mathbf{e}_{(3)} = \partial_z.\f] The natural
 static local tetrad is given by \f[ \mathbf{e}_{(0)} = \frac{1}{\sqrt{c^2-v_s^2f(r_s(t))^2}}\partial_t,\quad
 \mathbf{e}_{(1)}=\frac{v_sf(r_s(t))}{c\sqrt{c^2-v_s^2f(r_s(t))^2}}\partial_t+\frac{\sqrt{c^2-v_s^2f(r_s(t))^2}}{c}\partial_x,\quad
 \mathbf{e}_{(2)}=\partial_y,\quad \mathbf{e}_{(3)} = \partial_z.\f]

     Miguel Alcubierre,<br><b>"The warp drive: hyper-fast travel within general relativity,"</b><br>
     Classical Quantum Gravity <b>11</b>, L73--L77 (1994).<br>

     The default local tetrad is the comoving one
     (enum_nat_tetrad_default = enum_nat_tetrad_lnrf).

 * This file is part of the m4d-library.
 */
#ifndef M4D_METRIC_ALCUBIERRE_H
#define M4D_METRIC_ALCUBIERRE_H

#include "m4dMetric.h"

namespace m4d {

// ---------------------------------------------------
//    class definition:   MetricAlcubierre
// ---------------------------------------------------
class MetricAlcubierre : public Metric
{
public:
    MetricAlcubierre(double sigma = 8.0, double R = 1.0, double vs = 1.0);
    virtual ~MetricAlcubierre();

    // --------- public methods -----------
public:
    virtual bool calculateMetric(const double* pos);
    virtual bool calculateChristoffels(const double* pos);
    virtual bool calculateChrisD(const double* pos);

    virtual bool calculateRiemann(const double* pos);

    virtual void localToCoord(
        const double* pos, const double* ldir, double* dir, enum_nat_tetrad_type type = enum_nat_tetrad_default);
    virtual void coordToLocal(
        const double* pos, const double* cdir, double* ldir, enum_nat_tetrad_type type = enum_nat_tetrad_default);

    virtual bool breakCondition(const double* pos);

    virtual bool calcDerivs(const double y[], double dydx[]);

    virtual double testConstraint(const double y[], const double kappa);
    virtual bool resize(double* y, double kappa, double factor = DEF_RESIZE_FACTOR);

    virtual bool setParam(const char* pName, double val);

    virtual bool transToTwoPlusOne(vec4 p, vec4& cp);

    virtual bool report(const vec4 pos, const vec4 cdir, char*& text);

    // --------- specific public methods ----------
public:
    virtual double calcRs(const double* pos);
    virtual double calcF(const double* pos);

    // --------- protected methods -----------
protected:
    virtual void setStandardValues();
    virtual void calcDF(const double* pos, double& ft, double& fx, double& fy, double& fz);
    virtual void calcD2F(const double* pos, double& ftt, double& ftx, double& fty, double& ftz, double& fxx,
        double& fxy, double& fxz, double& fyy, double& fyz, double& fzz);

    // -------- protected attribute ---------
protected:
    double mSigma;
    double mR;
    double mvs;
};

} // end namespace m4d

#endif
