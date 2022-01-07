/**
 * @file    m4dMetricKasner.h
 * @author  Thomas Mueller
 *
 * @brief  Kasner metric in cartesian coordinates (t,x,y,z).

     The line element is given by

     \f[ds^2 = -dt^2 + t^{2p_1}dx^2 + t^{2p_2}dy^2 + t^{2p_3}dz^2.\f]

     The parameters \f$ p_1,p_2,p_3\f$ can be represented by the Khalatnikov-Lifshitz parameter \f$u\f$:
     \f[ p_1 = -\frac{u}{1+u+u^2},\quad p_2 = \frac{1+u}{1+u+u^2}, \quad p_3 = \frac{u(1+u)}{1+u+u^2}.\f]

     Detailed discussions about the Kasner metric can be found in MTW.

 * This file is part of the m4d-library.
 */
#ifndef M4D_METRIC_KASNER_H
#define M4D_METRIC_KASNER_H

#include "m4dMetric.h"

namespace m4d {

// ---------------------------------------------------
//    class definition:   MetricKasner
// ---------------------------------------------------
class MetricKasner : public Metric
{
public:
    MetricKasner(double u = 0.0);
    virtual ~MetricKasner();

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

    virtual bool setParam(const char* pName, double val);

    virtual bool report(const vec4 pos, const vec4 cdir, char*& text);

    // --------- specific public methods ----------
public:
    void calc_parameters();

    // --------- protected methods -----------
protected:
    virtual void setStandardValues();

    // -------- protected attribute ---------
protected:
    double p1, p2, p3, mU;
};

} // end namespace m4d

#endif
