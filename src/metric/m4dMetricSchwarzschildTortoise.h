// -------------------------------------------------------------------------------
/*
    m4dMetricSchwarzschildTortoise.h

  Copyright (c) 2010-2014  Thomas Mueller


   This file is part of the m4d-library.

   The m4d-library is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   The m4d-library is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with the m4d-library.  If not, see <http://www.gnu.org/licenses/>.

*/

/*!  \class  m4d::MetricSchwarzschildTortoise
     \brief  Schwarzschild metric in tortoise coordinates (t,rho,theta,phi).

             The line element is given by

             \f[ds^2 = -\left(1-\frac{r_s}{r(\rho)}\right) c^2dt^2 + \left(1-\frac{r_s}{r(\rho)}\right) d\rho^2 + r(\rho)^2\left(\vartheta^2 + \sin(\vartheta)^2 \varphi^2\right),\f]
             where \f$r_s=2GM/c^2\f$ is the Schwarzschild radius. G is Newton's
             constant, M is the mass of the black hole, and c is the speed of light.

             The tortoise and Schwarzschild radial coordinates are related by
	     \f[ \rho=r+r_s\ln\left(\frac{r}{r_s}-1\right),\quad\mbox{or}\quad r=r_s\left\{1+\mathcal{W}\left[\exp\left(\frac{\rho}{r_s}-1\right)\right]\right\}. \f]

             The natural local tetrad is given by
             \f[ \mathbf{e}_{(0)} = \frac{1}{c\,\sqrt{1-r_s/r(\rho)}}\partial_t,\quad \mathbf{e}_{(1)}=\frac{1}{c\,\sqrt{1-r_s/r(\rho)}}\partial_{\rho},\quad \mathbf{e}_{(2)}=\frac{1}{r(\rho)}\partial_{\vartheta},\quad \mathbf{e}_{(3)} = \frac{1}{r(\rho)\sin\vartheta}\partial_{\varphi}.\f]

             The embedding diagram (Flamm's paraboloid) is defined by the function
             \f[ z = 2\sqrt{r_s}\sqrt{r(\rho)-r_s}. \f]

             Detailed discussions about the Schwarzschild metric can be found
             in the standard literature, e.g. \ref lit_wald "Wald", \ref lit_rindler "Rindler", \ref lit_mtw "MTW".

*/
// -------------------------------------------------------------------------------

#ifndef M4D_METRIC_SCHWARZSCHILD_TORTOISE_H
#define M4D_METRIC_SCHWARZSCHILD_TORTOISE_H

#include "m4dMetric.h"

#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_lambert.h>

namespace m4d {

// ---------------------------------------------------
//    class definition:   MetricSchwarzschildTortoise
// ---------------------------------------------------
class MetricSchwarzschildTortoise : public Metric {
public:
    MetricSchwarzschildTortoise(double mass = 1.0);
    virtual ~MetricSchwarzschildTortoise();

// --------- public methods -----------
public:
    virtual bool   calculateMetric(const double* pos);
    virtual bool   calculateChristoffels(const double* pos);
    virtual bool   calculateChrisD(const double* pos);

    virtual void   localToCoord(const double* pos, const double* ldir, double* dir,
                                enum_nat_tetrad_type  type = enum_nat_tetrad_default);
    virtual void   coordToLocal(const double* pos, const double* cdir, double* ldir,
                                enum_nat_tetrad_type  type = enum_nat_tetrad_default);

    virtual bool   breakCondition(const double* pos);

    //virtual bool   calcDerivs             ( const double y[], double dydx[] );

    virtual double testConstraint(const double y[], const double kappa);

    virtual bool   calcProduct(const double* pos, const double* u, const double* v, double &prod, bool preCalcMetric = true);

    virtual bool   setParam(const char* pName, double val);

    virtual int    transToPseudoCart(vec4 p, vec4 &cp);
    virtual bool   transToEmbedding(vec4 p, vec4 &ep);

    virtual bool   setEmbeddingParam(const char* name, double val);
    virtual bool   testEmbeddingParams();
//    virtual int    getEmbeddingVertices(std::vector<vec3> &verts,
//                                        std::vector<int> &indices, unsigned int &numElems, unsigned int &counter);
    virtual bool   report(const vec4 pos, const vec4 cdir, std::string &text);


// --------- protected methods -----------
protected:
    virtual void  setStandardValues();

    double  calc_r(const double rho);

// -------- protected attribute ---------
protected:
    //! Schwarzschild radius rs=2GM/c^2.
    double rs;
    double mMass;

    double mEmb_rmin;
    double mEmb_rmax;
    double mEmb_rstep;
    double mEmb_phistep;
    double mEmb_r_num;
    double mEmb_phi_num;
};

} // end namespace m4d

#endif
