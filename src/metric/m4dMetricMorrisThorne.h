// -------------------------------------------------------------------------------
/*
    m4dMetricMorrisThorne.h

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

/*!  \class  m4d::MetricMorrisThorne
     \brief  The Morris-Thorne wormhole metric in spherical coordinates
             (t,l,theta,phi) with proper radial coordinate l.

             The line element is given by

             \f[ds^2 = -c^2 dt^2 + dl^2 + (l^2+b_0^2)\left(d\vartheta^2+\sin(\vartheta)^2 d\varphi^2\right),\f]

             where \f$b_0\f$ is the throat radius and c is the speed of light.

             The natural local tetrad is given by
             \f[ \mathbf{e}_{(0)} = \frac{1}{c}\partial_t,\quad \mathbf{e}_{(1)}=\partial_l,\quad \mathbf{e}_{(2)}=\frac{1}{\sqrt{l^2+b_0^2}}\partial_{\vartheta},\quad \mathbf{e}_{(3)} = \frac{1}{\sqrt{l^2+b_0^2}\sin\vartheta}\partial_{\varphi}.\f]

             The embedding diagram is defined by the function
             \f[z = \pm b_0\log\left[\frac{r}{b_0}+\sqrt{\left(\frac{r}{b_0}\right)-1}\right],\f]

             where \f$r=\sqrt{b_0^2+l^2}\f$.

             Detailed discussions about this metric can be found in <br><br>
             M.S. Morris and K.S. Thorne,<br><b>"Wormholes in spacetime and their use
             for interstellar travel: A tool for teaching general relativity,"</b><br>
             American Journal of Physics <b>56</b>, 395 (1988).
*/
// -------------------------------------------------------------------------------

#ifndef M4D_METRIC_MORRISTHORNE_H
#define M4D_METRIC_MORRISTHORNE_H

#include "m4dMetric.h"

namespace m4d {

// ---------------------------------------------------
//    class definition:   MetricMorrisThorne
// ---------------------------------------------------
class MetricMorrisThorne : public Metric {
public:
    MetricMorrisThorne(double b0 = 1.0);
    virtual ~MetricMorrisThorne();

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
    virtual int    transToPseudoCart(vec4 p, vec4 &cp);

    virtual bool   calcDerivs(const double y[], double dydx[]);
    virtual double testConstraint(const double y[], const double kappa);

    virtual bool   setParam(const char* pName, double val);

    virtual bool   transToEmbedding(vec4 p, vec4 &ep);
    virtual bool   transToCustom(vec4 p, vec4 &cp);

    virtual bool   setEmbeddingParam(const char* name, double val);
    virtual int    getEmbeddingVertices(std::vector<vec3> &verts,
                                        std::vector<int> &indices, unsigned int &numElems, unsigned int &counter);

    virtual bool   effPotentialValue(const vec4 pos, const vec4 cdir, enum_geodesic_type type, const double x, double &val);
    virtual bool   totEnergy(const vec4 pos, const vec4 cdir, const double x, double &val);

    virtual bool   report(const vec4 pos, const vec4 cdir, std::string &text);


// --------- specific public methods ----------
public:
    bool  calcKsiCrit(const vec4 pos, double &ksicrit);
    bool  constsOfMotion(const vec4 pos, const vec4 cdir, double &k, double &h);


// --------- protected methods -----------
protected:
    virtual void setStandardValues();

// -------- protected attribute ---------
protected:
    double mb0;

    double mEmb_lmin;
    double mEmb_lmax;
    double mEmb_lstep;
    double mEmb_phistep;
    double mEmb_l_num;
    double mEmb_phi_num;
};

} // end namespace m4d
#endif
