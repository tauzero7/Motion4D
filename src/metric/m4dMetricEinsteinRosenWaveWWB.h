/**
 * @file    m4dMetricEinsteinRosenWaveWWB.h
 * @author  Thomas Mueller
 *
 * @brief  Einstein-Rosen wave with a Weber-Wheeler-Bonnor pulse  in cylindrical coordinates (t,rho,phi,z).

     Detailed discussions can be found in<br><br>

     J.B. Griffiths and S. Micciche<br>
     "The Weber-Wheeler-Bonnor pulse and phase shifts in gravitational soliton interactions"<br>
     Physics Letters A 223, 37-42 (1997)

 * This file is part of the m4d-library.
 */
#ifndef M4D_METRIC_EINSTEIN_ROSEN_WAVE_WWB_H
#define M4D_METRIC_EINSTEIN_ROSEN_WAVE_WWB_H

#include "m4dMetric.h"

namespace m4d {

class MetricEinsteinRosenWaveWWB : public Metric
{
public:
    MetricEinsteinRosenWaveWWB(double c = 1.0, double a = 1.0);
    virtual ~MetricEinsteinRosenWaveWWB();

public:
    virtual bool calculateMetric(const double* pos);
    virtual bool calculateChristoffels(const double* pos);
    virtual bool calculateChrisD(const double* pos);

    virtual void localToCoord(
        const double* pos, const double* ldir, double* dir, enum_nat_tetrad_type type = enum_nat_tetrad_default);
    virtual void coordToLocal(
        const double* pos, const double* cdir, double* ldir, enum_nat_tetrad_type type = enum_nat_tetrad_default);

    virtual bool breakCondition(const double* pos);
    virtual double testConstraint(const double* y, const double kappa);

    virtual bool setParam(const char* pName, double val);

    virtual bool report(const vec4 pos, const vec4 cdir, char*& text);

protected:
    virtual void setStandardValues();

    void calcPotentials(const double* pos, double& gam, double& psi);
    void calcDiffPoti(
        const double* pos, double& gam, double& psi, double& gamt, double& gamr, double& psit, double& psir);

protected:
    double m_c;
    double m_a;
};

} // end namespace m4d

#endif // M4D_METRIC_EINSTEIN_ROSEN_WAVE_WWB_H
