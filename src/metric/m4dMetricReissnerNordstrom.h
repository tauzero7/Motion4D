// -------------------------------------------------------------------------------
/*
    m4dMetricReissnerNordstrom.h

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

/*!  \class  m4d::MetricReissnerNordstrom
     \brief  Reissner-Nordstrom metric in spherical coordinates (t,r,theta,phi).

             The line element is given by

             \f[ds^2 = -\left(1-\frac{r_s}{r}+\frac{\rho Q^2}{r^2}\right) dt^2 + \frac{dr^2}{1-r_s/r+\rho Q^2/r^2} + r^2\left(d\vartheta^2 + \sin(\vartheta)^2 d\varphi^2\right),\f]

             where rs=2GM/c^2 is the Schwarzschild radius. G is Newton's
             constant, M is the mass and Q the charge of the black hole,
             c is the speed of light, and rho=G/(epsilon_0*c^4).

             Detailed discussions about the Reissner-Nordstrom metric can be found
             in the standard literature, e.g. MTW.

*/
// -------------------------------------------------------------------------------

#ifndef M4D_METRIC_REISSNER_NORDSTROM_H
#define M4D_METRIC_REISSNER_NORDSTROM_H

#include "m4dMetric.h"
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>

/*! \brief Parameters for the ReissnerNordstrom metric.

  Is used to numerical calculate the embedding function.
 */
typedef struct {
    double mass;
    double rho;
    double q;
} struct_reissner_params;

namespace m4d {

// ---------------------------------------------------
//    class definition:   MetricReissnerNordstrom
// ---------------------------------------------------
class MetricReissnerNordstrom : public Metric {
public:
    MetricReissnerNordstrom(double mass = 1.0, double q = 0.1);
    virtual ~MetricReissnerNordstrom();

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

    virtual bool   calcDerivs(const double yn[], double dydx[]);

    virtual double testConstraint(const double y[], const double kappa);

    virtual bool   setParam(const char* pName, double val);

    virtual bool   transToEmbedding(vec4 p, vec4 &ep);

    virtual bool   setEmbeddingParam(const char* name, double val);
    virtual bool   testEmbeddingParams();
//    virtual int    getEmbeddingVertices(std::vector<vec3> &verts,
//                                        std::vector<int> &indices, unsigned int &numElems, unsigned int &counter);

    virtual bool   effPotentialValue(const vec4 pos, const vec4 cdir, enum_geodesic_type type, const double x, double &val);
    virtual bool   totEnergy(const vec4 pos, const vec4 cdir, const double x, double &val);

    virtual void   usePhysicalUnits(const enum_physical_constants  units);
    virtual void   setUnits(const double speed_of_light, const double grav_const, const double diel_perm);

    virtual bool   report(const vec4 pos, const vec4 cdir, char*&text);

// --------- specific public methods ----------
public:
    bool  calcKsiCrit(const vec4 pos, double &ksicrit, double B);
    bool  calcRadiusIsco(double &radiusIsco, double B);
    bool  calcRadiusIsco2(double &radiusIsco2, double B);
    bool  calcRadiusIuco(double &radiusIsco2, double B);
    bool  calcBetaCo(double &betaCo, double radiusCo, double B);

// --------- protected methods -----------
protected:
    virtual void   setStandardValues();
    void   calcDiskr();
    void   calcCritical();
    void   calcExtremals();

    bool   calcEmbeddingZ(const double r, double &z);

// -------- protected attribute ---------
protected:
    //! Schwarzschild radius rs=2GM/c^2.
    double rs;
    double mMass;

    //! k = G/(epsilon0*c^4);
    double mK;
    //! Charge of the black hole.
    double mQ;
    //! 1.0-4*mK*Q^2/rs^2.
    double mDiskr;

    // critical points
    double rp, rm;

    // extremal points for null geodesics
    double rNp, rNm;

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
