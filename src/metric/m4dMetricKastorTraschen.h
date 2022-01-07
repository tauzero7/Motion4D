/**
 * @file    m4dMetricKastorTraschen.h
 * @author  Thomas Mueller
 *
 * @brief  Kastor-Traschen metric in cartesian coordinates (t,x,y,z).
 *
 * This file is part of the m4d-library.
 */
#ifndef M4D_METRIC_KASTOR_TRASCHEN_H
#define M4D_METRIC_KASTOR_TRASCHEN_H

#include "m4dMetric.h"

namespace m4d {

/**
 * @brief The MetricKastorTraschen class
 */
class MetricKastorTraschen : public Metric
{
public:
    MetricKastorTraschen(double H = 0.0);
    virtual ~MetricKastorTraschen();

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

public:
    void calcPotentials(const double* pos, double& Omega, double& a);
    void calcPotDiffs(const double* pos, double& Omega, double& a, double& dOdx, double& dOdy, double& dOdz,
        double& dOdt, double& dadt);

protected:
    virtual void setStandardValues();

protected:
    double m1, z1; // mass and position of first black hole
    double m2, z2; // mass and position of second black hole
    double mH;
};

} // end namespace m4d

#endif
