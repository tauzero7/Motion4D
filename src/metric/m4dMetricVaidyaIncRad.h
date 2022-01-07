/**
 * @file    m4dMetricVaidyaIncRad.h
 * @author  Thomas Mueller
 *
 * @brief  Incoming radiation Vaidya metric in spherical coordinates (v,r,theta,phi).

     The line element is given by

     \f[ds^2 = 2dvdr - \left(1-\frac{2m(v)}{r}\right)dv^2 + r^2(d\theta^2+\sin^2\theta d\phi^2).\f]

     Detailed discussions about the Vaidya metric can be found in Griffiths,Podolsky.

     The mass function m(v) is taken from Piesnack and Kassner, https://arxiv.org/pdf/2103.08340v2.pdf
     with a1 = a2 = 0.

 * This file is part of the m4d-library.
 */
#ifndef M4D_METRIC_VAIDYA_INC_RAD_H
#define M4D_METRIC_VAIDYA_INC_RAD_H

#include "m4dMetric.h"

namespace m4d {

/**
 * @brief The MetricVaidyaIncRad class
 */
class MetricVaidyaIncRad : public Metric
{
public:
    MetricVaidyaIncRad(double k = 1.0);
    virtual ~MetricVaidyaIncRad();

    // --------- public methods -----------
public:
    virtual bool calculateMetric(const double* pos);
    virtual bool calculateChristoffels(const double* pos);

    virtual void localToCoord(
        const double* pos, const double* ldir, double* dir, enum_nat_tetrad_type type = enum_nat_tetrad_default);
    virtual void coordToLocal(
        const double* pos, const double* cdir, double* ldir, enum_nat_tetrad_type type = enum_nat_tetrad_default);

    virtual bool breakCondition(const double* pos);

    virtual bool report(const vec4 pos, const vec4 cdir, char*& text);

    virtual bool setParam(const char* pName, double val);

    // --------- protected methods -----------
protected:
    virtual void setStandardValues();

    void calcMassFunc(const double v, double& m);
    void calcMassFunc(const double v, double& m, double& dmdv);

    // -------- protected attribute ---------
protected:
    double m_k; // constant factor
    double m_vl;
    double m_sigma;
};

} // end namespace m4d

#endif // M4D_METRIC_VAIDYA_INC_RAD_H
