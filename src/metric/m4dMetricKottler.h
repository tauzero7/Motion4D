/**
 * @file    m4dMetricKottler.h
 * @author  Thomas Mueller
 *
 * @brief  Kottler metric in spherical coordinates (t,r,theta,phi).

     The line element is given by

     \f[ds^2 = -\left(1-\frac{r_s}{r}-\frac{\Lambda r^2}{3}\right) dt^2 + \frac{dr^2}{1-r_s/r-\Lambda r^2/3} +
 r^2\left(d\vartheta^2 + \sin(\vartheta)^2 d\varphi^2\right),\f]

     where \f$r_s=2GM/c^2\f$ is the Schwarzschild radius. G is Newton's
     constant, M is the mass of the black hole, Lambda is the cosmological
     constant, and c is the speed of light.

     Detailed discussions about the Kottler metric can be found
     in e.g. Phys.Rev.D 76, 043006 (2007) and in the original work by
     Kottler...

 * This file is part of the m4d-library.
 */
#ifndef M4D_METRIC_KOTTLER_H
#define M4D_METRIC_KOTTLER_H

#include "m4dMetric.h"
#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>

typedef struct {
    double rs;
    double lambda;
} struct_kottler_params;

namespace m4d {

class MetricKottler : public Metric
{
public:
    MetricKottler(double mass = 1.0, double lambda = 0.05);
    virtual ~MetricKottler();

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

    virtual bool transToEmbedding(vec4 p, vec4& ep);

    virtual bool setEmbeddingParam(const char* name, double val);
    virtual bool testEmbeddingParams();

    virtual void usePhysicalUnits(const enum_physical_constants units);
    virtual void setUnits(const double speed_of_light, const double grav_const, const double diel_perm);

    virtual bool effPotentialValue(
        const vec4 pos, const vec4 cdir, enum_geodesic_type type, const double x, double& val);
    virtual bool totEnergy(const vec4 pos, const vec4 cdir, const double x, double& val);

    virtual bool report(const vec4 pos, const vec4 cdir, char*& text);

    // --------- protected methods -----------
protected:
    virtual void setStandardValues();

    void calcCriticalPoints();
    bool calcEmbeddingZ(const double r, double& z);

    // -------- protected attribute ---------
protected:
    double rs;
    double mMass;
    double mLambda;

    // critical points
    double r1, rp, rm;

    double mEmb_rmin;
    double mEmb_rmax;
    double mEmb_rstep;
    double mEmb_phistep;
    double mEmb_r_num;
    double mEmb_phi_num;

    gsl_integration_workspace* w;
    gsl_function F;
};

} // end namespace m4d

#endif
