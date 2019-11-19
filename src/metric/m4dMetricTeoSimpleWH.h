// -------------------------------------------------------------------------------
/*
    m4dMetricTeoSimpleWH.h

  Copyright (c) 2015  Thomas Mueller


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

/*!  \class  m4d::MetricTeoSimpleWH
     \brief  Teo metric of an axisymmetric rotating wormhole in spherical coordinates (t,l,theta,phi).

*/
// -------------------------------------------------------------------------------

#ifndef M4D_METRIC_TEO_SIMPLE_WH_H
#define M4D_METRIC_TEO_SIMPLE_WH_H

#include "m4dMetric.h"

namespace m4d {

class MetricTeoSimpleWH : public Metric {
public:
    MetricTeoSimpleWH(double b0 = 1.0);
    virtual ~MetricTeoSimpleWH();

    // --------- public methods -----------
public:
    virtual bool   calculateMetric(const double* pos);
    virtual bool   calculateChristoffels(const double* pos);


    virtual void   localToCoord(const double* pos, const double* ldir, double* dir,
                                enum_nat_tetrad_type  type = enum_nat_tetrad_default);
    virtual void   coordToLocal(const double* pos, const double* cdir, double* ldir,
                                enum_nat_tetrad_type  type = enum_nat_tetrad_default);

    virtual bool   breakCondition(const double* pos);
    virtual int    transToPseudoCart(vec4 p, vec4 &cp);

    virtual bool   setParam(const char* pName, double val);

    virtual bool   report(const vec4 pos, const vec4 cdir, char*&text);

    virtual bool   transToEmbedding(vec4 p, vec4 &ep);
    virtual bool   setEmbeddingParam(const char *name, double val);
//    virtual int    getEmbeddingVertices(std::vector<vec3> &verts,
//                                        std::vector<int> &indices, unsigned int &numElems, unsigned int &counter);


    // --------- protected methods -----------
protected:
    virtual void setStandardValues();

    double getShapeVal(double r);

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

#endif // M4D_METRIC_TEO_SIMPLE_WH_H
