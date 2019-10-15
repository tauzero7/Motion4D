/**
 * @file    m4dMetricCosmicStringSchwarzschild.h
 * @author  Thomas Mueller
 *
 * @brief  Cosmic string within the Schwarzschild metric in spherical coordinates (t,r,theta,phi).

     The line element is given by

     \f[ds^2 = -\left(1-\frac{r_s}{r}\right) c^2dt^2 + \frac{dr^2}{1-r_s/r} + r^2\left(\vartheta^2 +
 \beta^2\sin(\vartheta)^2 \varphi^2\right),\f] where \f$r_s=2GM/c^2\f$ is the Schwarzschild radius. G is Newton's
     constant, M is the mass of the black hole, c is the speed of light,
     and \f$\beta\f$ is the string parameter.

     The natural local tetrad is given by
     \f[ \mathbf{e}_{(0)} = \frac{1}{c\,\sqrt{1-r_s/r}}\partial_t,\quad
 \mathbf{e}_{(1)}=\sqrt{1-\frac{r_s}{r}}\partial_r,\quad \mathbf{e}_{(2)}=\frac{1}{r}\partial_{\vartheta},\quad
 \mathbf{e}_{(3)} = \frac{1}{r\beta\sin\vartheta}\partial_{\varphi}.\f]

     Detailed discussions about the Schwarzschild metric with a cosmic string can be found in<br>
     M. Aryal, L. H. Ford, and A. Vilenkin,<br>
     <b>Cosmic strings and black holes</b><br>
     Phys. Rev. D, 34(8):2263–2266, Oct 1986.
 *
 *  This file is part of libMotion4D.
 */
#ifndef M4D_METRIC_COSMIC_STRING_SCHWARZSCHILD_H
#define M4D_METRIC_COSMIC_STRING_SCHWARZSCHILD_H

#include "m4dMetric.h"

namespace m4d {

/**
 * @brief The MetricCosmicStringSchwarzschild class
 */
class MetricCosmicStringSchwarzschild : public Metric
{
public:
    MetricCosmicStringSchwarzschild(double mass = 1.0, double beta = 0.5);
    virtual ~MetricCosmicStringSchwarzschild();

public:
    virtual bool calculateMetric(const double* pos);
    virtual bool calculateChristoffels(const double* pos);
    virtual bool calculateChrisD(const double* pos);

    virtual void localToCoord(
        const double* pos, const double* ldir, double* dir, enum_nat_tetrad_type type = enum_nat_tetrad_default);
    virtual void coordToLocal(
        const double* pos, const double* cdir, double* ldir, enum_nat_tetrad_type type = enum_nat_tetrad_default);

    virtual bool breakCondition(const double* pos);

    virtual bool calcDerivs(const double y[], double dydx[]);

    virtual double testConstraint(const double y[], const double kappa);

    virtual bool calcProduct(
        const double* pos, const double* u, const double* v, double& prod, bool preCalcMetric = true);

    virtual bool setParam(const char* pName, double val);

    virtual bool transToTwoPlusOne(vec4 p, vec4& cp);

    virtual bool transToEmbedding(vec4 p, vec4& ep);
    virtual bool setEmbeddingParam(const char* name, double val);
    virtual bool testEmbeddingParams();

    virtual unsigned int getEmbeddingVertices(
        float*& verts, unsigned int*& indices, unsigned int& numElems, unsigned int& counter);

    virtual bool effPotentialValue(
        const vec4 pos, const vec4 cdir, enum_geodesic_type type, const double x, double& val);
    virtual bool totEnergy(const vec4 pos, const vec4 cdir, const double x, double& val);

    virtual bool report(const vec4 pos, const vec4 cdir, char*& text);

public:
    bool calcKsiCrit(const vec4 pos, double& ksicrit);

protected:
    virtual void setStandardValues();
    double embFunc(const double r);

protected:
    //! Schwarzschild radius rs=2GM/c^2.
    double rs;
    double mMass;
    double mBeta;

    double mEmb_rmin;
    double mEmb_rmax;
    double mEmb_rstep;
    double mEmb_phistep;
    unsigned int mEmb_r_num;
    unsigned int mEmb_phi_num;
};

} // end namespace m4d

#endif
