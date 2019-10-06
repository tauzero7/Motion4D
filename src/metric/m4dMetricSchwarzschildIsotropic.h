/**
 * @file    m4dMetricSchwarzschildIsotropic.h
 * @author  Thomas Mueller
 *
 * @brief  Schwarzschild metric in Cartesian isotropic coordinates (t,x,y,z).

    The line element is given by

    \f[ds^2 = -\left(\frac{1-\rho_s/\rho}{1+\rho_s/\rho}\right)^2 c^2dt^2 +
 \left(1+\frac{\rho_s}{\rho}\right)^4\left(dx^2+dy^2+dz^2\right),\f] where \f$\rho^2=x^2+y^2+z^2\f$ and
 \f$\rho_s=GM/(2c^2)\f$ is the Schwarzschild radius. G is Newton's constant, M is the mass of the black hole, and c is
 the speed of light.

    The natural local tetrad is given by
    \f[ \mathbf{e}_{(0)} = \frac{1+\rho_s/\rho}{1-\rho_s/\rho}\frac{\partial_t}{c},\quad
 \mathbf{e}_{(1)}=\left(1+\frac{\rho_s}{\rho}\right)^{-2}\partial_x,\quad
 \mathbf{e}_{(2)}=\left(1+\frac{\rho_s}{\rho}\right)^{-2}\partial_y,\quad \mathbf{e}_{(3)} =
 \left(1+\frac{\rho_s}{\rho}\right)^{-2}\partial_z.\f]

    Detailed discussions about the Schwarzschild metric can be found
    in the standard literature, e.g. \ref lit_wald "Wald", \ref lit_rindler "Rindler", \ref lit_mtw "MTW".
 *
 *  This file is part of libMotion4D.
 */
#ifndef M4D_METRIC_SCHWARZSCHILD_ISOTROPIC_H
#define M4D_METRIC_SCHWARZSCHILD_ISOTROPIC_H

#include "m4dMetric.h"

#include <gsl/gsl_errno.h>

namespace m4d {

/**
 * @brief The MetricSchwarzschildIsotropic class
 */
class MetricSchwarzschildIsotropic : public Metric
{
public:
    MetricSchwarzschildIsotropic(double mass = 1.0);
    virtual ~MetricSchwarzschildIsotropic();

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

    virtual bool calcProduct(
        const double* pos, const double* u, const double* v, double& prod, bool preCalcMetric = true);

    virtual bool setParam(const char* pName, double val);

    virtual bool transToEmbedding(vec4 p, vec4& ep);
    virtual bool setEmbeddingParam(const char* name, double val);
    virtual bool testEmbeddingParams();

    virtual unsigned int getEmbeddingVertices(
        float*& verts, unsigned int*& indices, unsigned int& numElems, unsigned int& counter);

    virtual bool report(const vec4 pos, const vec4 cdir, char*& text);

    bool calcBetaOfCircOrbit(const double* pos, double& beta);

protected:
    virtual void setStandardValues();

    double calc_rho(const double* pos);
    void calc_drho(const double* pos, double& drdx, double& drdy, double& drdz);
    void calc_d2rho(const double* pos, double& drdxdx, double& drdxdy, double& drdxdz, double& drdydy, double& drdydz,
        double& drdzdz);

    void calc_orbits();

protected:
    //! Schwarzschild radius rs=2GM/c^2.
    double rho_s;
    double mMass;

    double rho_po; // photon orbit;
    double rho_lso; // last stable timelike circular orbit

    double rho;

    double mEmb_rmin;
    double mEmb_rmax;
    double mEmb_rstep;
    double mEmb_phistep;
    unsigned int mEmb_r_num;
    unsigned int mEmb_phi_num;
};

} // end namespace m4d

#endif
