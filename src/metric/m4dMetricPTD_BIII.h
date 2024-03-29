/**
 * @file    m4dMetricPTD_BIII.h
 * @author  Felix Beslmeisl
 *
 * @brief  \ref lit_stephani "Exact Solutions:" Petrov type D solutions - Case BIII.

     The line element is given by

     \f[ds^2 = z^2\left( dr^2 - r^2 dt^2 \right) + z dz^2+\frac{1}{z}d\varphi^2.\f]

     The natural local tetrad is given by
     \f[  \mathbf{e}_{(t)} = \frac{1}{zr}\partial_t, \quad
          \mathbf{e}_{(r)} = \frac{1}{z}\partial_r, \quad
          \mathbf{e}_{(\varphi)} = \frac{1}{\sqrt{z}}\partial_{\varphi}, \quad
          \mathbf{e}_{(z)} = \sqrt{z}\partial_{z}.\f]

* This file is part of the m4d-library.
*/
#ifndef M4DMETRICPTDBIII_H
#define M4DMETRICPTDBIII_H

#include "m4dMetric.h"

namespace m4d {

// ---------------------------------------------------
//    class definition:   MetricKerrBL
// ---------------------------------------------------
class MetricPTD_BIII : public Metric
{
public:
    MetricPTD_BIII();
    virtual ~MetricPTD_BIII();

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

    // --------- protected methods -----------
protected:
    virtual void setStandardValues();
};

} // end namespace m4d

#endif // M4DMETRICPTDBIII_H
