// -------------------------------------------------------------------------------
/*
    m4dMetricSchwarzschild.h

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

/*!  \class  m4d::MetricSchwarzschild
     \brief  Schwarzschild metric in spherical coordinates (t,r,theta,phi).

             The line element is given by

             \f[ds^2 = -\left(1-\frac{r_s}{r}\right) c^2dt^2 + \frac{dr^2}{1-r_s/r} + r^2\left(d\vartheta^2 + \sin^2\vartheta\,d\varphi^2\right),\f]
             where \f$r_s=2GM/c^2\f$ is the Schwarzschild radius. G is Newton's
             constant, M is the mass of the black hole, and c is the speed of light.

             The natural local tetrad is given by
             \f[ \mathbf{e}_{(0)} = \frac{1}{c\,\sqrt{1-r_s/r}}\partial_t,\quad \mathbf{e}_{(1)}=\sqrt{1-\frac{r_s}{r}}\partial_r,\quad \mathbf{e}_{(2)}=\frac{1}{r}\partial_{\vartheta},\quad \mathbf{e}_{(3)} = \frac{1}{r\sin\vartheta}\partial_{\varphi}.\f]

             The embedding diagram (Flamm's paraboloid) is defined by the function
             \f[ z = 2\sqrt{r_s}\sqrt{r-r_s}. \f]

             From the Euler-Lagrange formalism, we have
             \f[ \frac{1}{2}\dot{r}^2+V_{\mbox{eff}} = \frac{1}{2}\frac{k^2}{c^2},\qquad\mbox{with}\qquad V_{\mbox{eff}}=\frac{1}{2}\left(\frac{h^2}{r^2}-\kappa c^2\right)\left(1-\frac{r_s}{r}\right)\f]
             and the constants of motion \f$k=(1-r_s/r)c^2\dot{t}\f$ and \f$h=r^2\dot{\varphi}\f$.

             Detailed discussions about the Schwarzschild metric can be found
             in the standard literature, e.g. \ref lit_wald "Wald", \ref lit_rindler "Rindler", \ref lit_mtw "MTW".

*/
// -------------------------------------------------------------------------------

#ifndef M4D_METRIC_SCHWARZSCHILD_H
#define M4D_METRIC_SCHWARZSCHILD_H

#include "m4dMetric.h"

namespace m4d {

// ---------------------------------------------------
//    class definition:   MetricSchwarzschild
// ---------------------------------------------------
class MetricSchwarzschild : public Metric {
public:
    MetricSchwarzschild(double mass = 1.0);
    virtual ~MetricSchwarzschild();

    // --------- public methods -----------
public:
    virtual bool   calculateMetric(const double* pos);
    virtual bool   calculateChristoffels(const double* pos);
    virtual bool   calculateChrisD(const double* pos);

    virtual bool   calculateRiemann(const double* pos);
    virtual bool   calculateWeyl(const double* pos);
    virtual bool   calculateRicRotCoeffs(const double* pos);
    virtual bool   calculateContrRRC(const double* pos);

    virtual void   localToCoord(const double* pos, const double* ldir, double* dir,
                                enum_nat_tetrad_type  type = enum_nat_tetrad_default);
    virtual void   coordToLocal(const double* pos, const double* cdir, double* ldir,
                                enum_nat_tetrad_type  type = enum_nat_tetrad_default);

    virtual bool   breakCondition(const double* pos);

    virtual bool   calcDerivs(const double y[], double dydx[]);
    virtual bool   calcDerivsPar(const double y[], double dydx[]);
    virtual bool   calcDerivsSachsJacobi(const double y[], double dydx[]);
    virtual bool   calcDerivsFW(const double a[], const double y[], double dydx[]);

    virtual double testConstraint(const double y[], const double kappa);

    virtual bool   calcProduct(const double* pos, const double* u, const double* v, double &prod, bool preCalcMetric = true);

    virtual bool   setParam(std::string pName, double val);

    virtual bool   transToEmbedding(vec4 p, vec4 &ep);
    virtual bool   transToTwoPlusOne(vec4 p, vec4 &cp);
    virtual bool   transToCustom(vec4 p, vec4 &cp);

    virtual bool   setEmbeddingParam(std::string name, double val);
    virtual bool   testEmbeddingParams();
    virtual int    getEmbeddingVertices(std::vector<vec3> &verts,
                                        std::vector<int> &indices, unsigned int &numElems, unsigned int &counter);

    virtual bool   effPotentialValue(const vec4 pos, const vec4 cdir, enum_geodesic_type type, const double x, double &val);
    virtual bool   totEnergy(const vec4 pos, const vec4 cdir, const double x, double &val);

    virtual double getCircularVelocity(const double r, const enum_nat_tetrad_type  tedType = enum_nat_tetrad_default);
    virtual vec4   getCircularFourVel(const vec4 pos, const enum_nat_tetrad_type  tedType = enum_nat_tetrad_default);

    virtual void   usePhysicalUnits(const enum_physical_constants  units);
    virtual void   setUnits(const double speed_of_light, const double grav_const, const double diel_perm);

    virtual bool   report(const vec4 pos, const vec4 cdir, std::string &text);


    // --------- specific public methods ----------
public:
    bool  calcKsiCrit(const vec4 pos, double &ksicrit);

    double  calcKsiCrit(const double r);
    double  getConstOfMotion_k(const double r, const double beta, const int kappa = 0);
    double  getConstOfMotion_h(const double r, const double beta, const double ksi, const int kappa = 0);

    double param_aqua(const double k, const double h, const int kappa = 0);
    double param_b(const double k, const double h, const int kappa = 0);
    void   param_sq(const double k, const double h, const int kappa,
                    double &p, double &q, double &D1);
    double param_psi(const double p, const double q);
    double param_eps(const double r, const double ksi);

    bool  minDist_r(const double r, const double ksi, double &rMin);

    // --------- protected methods -----------
protected:
    virtual void  setStandardValues();
    virtual void  contrChrisVecVec(const double y[], const double v[], const double w[], double* z, bool calc = true);
    virtual void  contrChrDVecVecVec(const double y[], const double u[], const double v[], const double w[], double* z, bool calc = true);

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
