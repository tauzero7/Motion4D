// -------------------------------------------------------------------------------
/*
    m4dMetricKottler.h

  Copyright (c) 2009-2014  Thomas Mueller, Frank Grave


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

/*!  \class  m4d::MetricKottler
     \brief  Kottler metric in spherical coordinates (t,r,theta,phi).

             The line element is given by

             \f[ds^2 = -\left(1-\frac{r_s}{r}-\frac{\Lambda r^2}{3}\right) dt^2 + \frac{dr^2}{1-r_s/r-\Lambda r^2/3} + r^2\left(d\vartheta^2 + \sin(\vartheta)^2 d\varphi^2\right),\f]

             where \f$r_s=2GM/c^2\f$ is the Schwarzschild radius. G is Newton's
             constant, M is the mass of the black hole, Lambda is the cosmological
             constant, and c is the speed of light.

             Detailed discussions about the Kottler metric can be found
             in e.g. Phys.Rev.D 76, 043006 (2007) and in the original work by
             Kottler...

*/
// -------------------------------------------------------------------------------

#ifndef M4D_METRIC_KOTTLER_H
#define M4D_METRIC_KOTTLER_H

#include "m4dMetric.h"
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>

/*! \brief Parameters for the Kottler metric.

  Is used to numerical calculate the embedding function.
 */
typedef struct {
    double rs;
    double lambda;
} struct_kottler_params;

namespace m4d {

// ---------------------------------------------------
//    class definition:   MetricKottler
// ---------------------------------------------------
class MetricKottler : public Metric {
public:
    MetricKottler(double mass = 1.0, double lambda = 0.05);
    virtual ~MetricKottler();

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

    virtual double testConstraint(const double y[], const double kappa);

    virtual bool   setParam(const char* pName, double val);

    virtual bool   transToEmbedding(vec4 p, vec4 &ep);

    virtual bool   setEmbeddingParam(const char *name, double val);
    virtual bool   testEmbeddingParams();
//    virtual int    getEmbeddingVertices(std::vector<vec3> &verts,
//                                        std::vector<int> &indices, unsigned int &numElems, unsigned int &counter);

    virtual void   usePhysicalUnits(const enum_physical_constants  units);
    virtual void   setUnits(const double speed_of_light, const double grav_const, const double diel_perm);

    virtual bool   effPotentialValue(const vec4 pos, const vec4 cdir, enum_geodesic_type type, const double x, double &val);
    virtual bool   totEnergy(const vec4 pos, const vec4 cdir, const double x, double &val);

    virtual bool   report(const vec4 pos, const vec4 cdir, std::string &text);

// --------- protected methods -----------
protected:
    virtual void   setStandardValues();

    void   calcCriticalPoints();
    bool   calcEmbeddingZ(const double r, double &z);

// -------- protected attribute ---------
protected:
    double  rs;
    double  mMass;
    double  mLambda;

    // critical points
    double  r1, rp, rm;

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

