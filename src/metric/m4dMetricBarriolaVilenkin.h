// -------------------------------------------------------------------------------
/*
    m4dMetricBarriolaVilenkin.h

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

/*!  \class  m4d::MetricBarriolaVilenkin
     \brief  Barriola-Vilenkin metric in spherical coordinates (t,r,theta,phi).

             The line element is given by
             \f[ ds^2 = -c^2dt^2 + dr^2 + k^2r^2\left(d\vartheta^2+\sin^2\vartheta\,d\varphi^2\right). \f]

             The natural local tetrad is given by
             \f[ \mathbf{e}_{(0)} = \frac{1}{c}\partial_t,\quad \mathbf{e}_{(1)}=\partial_r,\quad \mathbf{e}_{(2)}=\frac{1}{kr}\partial_{\vartheta},\quad \mathbf{e}_{(3)} = \frac{1}{kr\sin\vartheta}\partial_{\varphi}.\f]

             Detailed discussions about the BarriolaVilenkin metric can be found in <br><br>
             Manuel Barriola and Alexander Vilenkin<br><b>"Gravitational Field of a Global Monopole"</b><br>Physical Review Letters <b>63</b>, 341--343 (1989).
*/
// -------------------------------------------------------------------------------

#ifndef M4D_METRIC_BARRIOLA_VILENKIN_H
#define M4D_METRIC_BARRIOLA_VILENKIN_H

#include "m4dMetric.h"

namespace m4d {

// ---------------------------------------------------
//    class definition:   BarriolaVilenkin
// ---------------------------------------------------
class MetricBarriolaVilenkin : public Metric {
public:
    MetricBarriolaVilenkin(double k = 0.8);
    virtual ~MetricBarriolaVilenkin();


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

    virtual bool   calcDerivs(const double y[], double dydx[]);
    virtual double testConstraint(const double y[], const double kappa);

    virtual bool   setParam(const char* pName, double val);

    virtual bool   transToEmbedding(vec4 p, vec4 &ep);

    virtual bool   setEmbeddingParam(const char *name, double val);
    virtual bool   testEmbeddingParams();
//    virtual int    getEmbeddingVertices(std::vector<vec3> &verts,
//                                        std::vector<int> &indices, unsigned int &numElems, unsigned int &counter);

    virtual bool   effPotentialValue(const vec4 pos, const vec4 cdir, enum_geodesic_type type, const double x, double &val);
    virtual bool   totEnergy(const vec4 pos, const vec4 cdir, const double x, double &val);

    virtual bool   report(const vec4 pos, const vec4 cdir, char*&text);

// --------- protected methods -----------
protected:
    virtual void    setStandardValues();
    double  closestApproach(const vec4 pos, const vec4 cdir);

// -------- protected attribute ---------
protected:
    double mk;

    double mEmb_rmin;
    double mEmb_rmax;
    double mEmb_rstep;
    double mEmb_phistep;
    double mEmb_r_num;
    double mEmb_phi_num;
};

} // end namespace m4d

#endif

