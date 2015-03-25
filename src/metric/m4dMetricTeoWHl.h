// -------------------------------------------------------------------------------
/*
    m4dMetricTeoWHl.h

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

/*!  \class  m4d::MetricTeoWHl
     \brief  Teo metric of an axisymmetric rotating wormhole in spherical coordinates (t,l,theta,phi).

             The line element is given by

             \f[ds^2 = -N(l,\vartheta) c^2dt^2 + dl^2 + r(l)^2K(l,\vartheta)^2 \left[d\vartheta^2 + \sin^2(\vartheta)\left(d\varphi-\omega(l,\vartheta)c dt\right)^2\right].\f]

             The several potentials can only be changed in the source code. The standard potentials are taken from [Fechtig]
             \f[ N(l,\vartheta)=1, \qquad K(l,\vartheta)=1, \qquad r(l)=\sqrt{l^2+b_0^2}, \qquad \omega(l,\vartheta)=\frac{b_0^2}{2\left(l^2+b_0^2\right)^{3/2}}\f]

             The default local tetrad is the locally nonrotating frame (enum_nat_tetrad_default = enum_nat_tetrad_lnrf).

             A detailed discussion about this metric can be found in the original article<br>

             Edward Teo,<br><b>"Rotating traversable wormholes"</b>,<br> Phys. Rev. D <b>58</b>, 024014 (1998).<br>

             and in the diploma thesis<br>

             Oliver Fechtig,<br><b>"Physikalische Aspekte und Visualisierung von stationaeren Wurmloechern",</b><br>
             <a href="http://itp1.uni-stuttgart.de/publikationen/?A=1&T=1#F">Click here</a> (in german).

*/
// -------------------------------------------------------------------------------

#ifndef M4D_METRIC_TEOWHL_H
#define M4D_METRIC_TEOWHL_H

#include "m4dMetric.h"

namespace m4d {

// ---------------------------------------------------
//    class definition:   MetricTeoWHl
// ---------------------------------------------------
class MetricTeoWHl : public Metric {
public:
    MetricTeoWHl(double b0 = 1.0);
    virtual ~MetricTeoWHl();

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

    virtual bool   setParam(std::string pName, double val);

    virtual bool   report(const vec4 pos, const vec4 cdir, std::string &text);

    // --------- protected methods -----------
protected:
    virtual void setStandardValues();

    double  calc_N(double l);
    double  calc_dN(double l);
    double  calc_d2N(double l);
    double  calc_K(double l);
    double  calc_dK(double l);
    double  calc_d2K(double l);
    double  calc_r(double l);
    double  calc_dr(double l);
    double  calc_d2r(double l);
    double  calc_omega(double l);
    double  calc_domega(double l);
    double  calc_d2omega(double l);

    // -------- protected attribute ---------
protected:
    double mb0;
    double pm;

    double curr_N, curr_dNdl;
    double curr_K, curr_dKdl;
    double curr_r, curr_drdl;
    double curr_omega, curr_domegadl;
};

} // end namespace m4d

#endif
