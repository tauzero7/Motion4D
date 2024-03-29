/**
 * @file    m4dMetricJaNeWi.h
 * @author  Thomas Mueller
 *
 * @brief  JaNeWi metric in spherical coordinates (t,r,theta,phi).

     The line element is given by

     \f[ds^2 = -\alpha^{\gamma} c^2dt^2 + \alpha^{-\gamma}dr^2 + r^2\alpha^{-\gamma+1}\left(d\vartheta^2+\sin^2\vartheta
 d\varphi^2\right).\f]


     Detailed discussions about the Janis-Newman-Winicour metric can be found in<br><br>
     A. I. Janis, E. T. Newman, and J. Winicour.<br>
     <b>Reality of the Schwarzschild singularity</b><br>
     Phys. Rev. Lett. <b>20</b>, 878–880 (1968).

 * This file is part of the m4d-library.
 */
#ifndef M4D_METRIC_JANEWI_H
#define M4D_METRIC_JANEWI_H

#include "m4dMetric.h"
#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>

typedef struct {
    double rs;
    double gamma;
} struct_janewi_params;

namespace m4d {

// ---------------------------------------------------
//    class definition:   MetricJaNeWi
// ---------------------------------------------------
class MetricJaNeWi : public Metric
{
public:
    MetricJaNeWi(double mass = 1.0, double gamma = 1.0);
    virtual ~MetricJaNeWi();

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

    virtual bool transToEmbedding(vec4 p, vec4& ep);

    virtual bool setEmbeddingParam(const char* name, double val);
    virtual bool testEmbeddingParams();
    //    virtual int    getEmbeddingVertices(std::vector<vec3> &verts,
    //                                        std::vector<int> &indices, unsigned int &numElems, unsigned int &counter);

    virtual double getCircularVelocity(const double r, const enum_nat_tetrad_type tedType = enum_nat_tetrad_default);
    virtual vec4 getCircularFourVel(const vec4 pos, const enum_nat_tetrad_type tedType = enum_nat_tetrad_default);

    virtual bool effPotentialValue(
        const vec4 pos, const vec4 cdir, enum_geodesic_type type, const double x, double& val);
    virtual bool totEnergy(const vec4 pos, const vec4 cdir, const double x, double& val);

    virtual bool report(const vec4 pos, const vec4 cdir, char*& text);

    // --------- specific public methods ----------
public:
    void calcCriticalPoint();
    bool calcEmbeddingZ(const double r, double& z);

    // --------- protected methods -----------
protected:
    virtual void setStandardValues();

    // -------- protected attribute ---------
protected:
    double rs;
    double mMass;
    double mGamma;
    double mCritPoint;

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
