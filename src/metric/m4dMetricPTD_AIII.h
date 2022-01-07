/**
 * @file    m4dMetricPTD_AIII.h
 * @author  Felix Beslmeisl
 *
 * @brief  \ref lit_stephani "Exact Solutions:" Petrov type D solutions - Case AIII.

     The line element is given by

     \f[ds^2 =z^2\left( dr^2 + r^2 d\varphi^2 \right) + z dz^2-\frac{1}{z}dt^2.\f]

     The natural local tetrad is given by
     \f[ \mathbf{e}_{(t)} = \sqrt{z}\partial_t, \quad
         \mathbf{e}_{(r)} = \frac{1}{z}\partial_r, \quad
         \mathbf{e}_{(\varphi)} = \frac{1}{zr}\partial_{\varphi}, \quad
         \mathbf{e}_{(z)} = \frac{1}{\sqrt{z}}\partial_{z}.\f]

 * This file is part of the m4d-library.
 */
#ifndef M4DMETRICPTDAIII_H
#define M4DMETRICPTDAIII_H

#include "m4dMetric.h"

namespace m4d {

class MetricPTD_AIII : public Metric
{
public:
    MetricPTD_AIII();
    virtual ~MetricPTD_AIII();

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
};

} // end namespace m4d

#endif // M4DMETRICPTDAIII_H
