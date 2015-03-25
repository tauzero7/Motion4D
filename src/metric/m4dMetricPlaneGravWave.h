// -------------------------------------------------------------------------------
/*
   m4dMetricPlaneGravWave.h

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

/*!  \class  m4d::MetricPlaneGravWave
     \brief  Metric for a plane gravitational wave with finite logitudinal and
             infinite transversal extension (sandwich wave).

             With \f$ u=ct-x \f$ the line element is given by
             \f[ ds^2 = -c^2dt^2 + dx^2 + p^2\left(u\right)dy^2 + q^2\left(u\right)dz^2. \f]
             The functions \f$ p\left(u\right) \f$ and \f$ q\left(u\right) \f$
             are represented through
             \f[ p\left(u\right)=\left\{\begin{array}{lr}
                                          p_0 = \mathrm{const.} \quad &-a>u  \\
                                          L\left(u\right)\mathrm{e}^{m\left(u\right)} \quad  &-a\leq u \leq 0  \\
                                          1-u \quad &0<a
                                        \end{array}\right. \f]
             \f[ q\left(u\right)=\left\{\begin{array}{lr}
                                          q_0 = \mathrm{const.} \quad &-a>u  \\
                                          L\left(u\right)\mathrm{e}^{-m\left(u\right)} \quad  &-a\leq u \leq 0 \\
                                          1-u \quad &0<a
                                        \end{array}\right., \f]
             where the parameter \f$ a \f$ characterize the longitudinal
             extension of the wave. The Functions \f$ L\left(u\right) \f$ and
             \f$ m\left(u\right) \f$ are given by
             \f[ L\left(u\right) = 1 - u + \frac{1}{a^3}u^3 + \frac{1}{2a^3}u^4 \f]
             \f[ m\left(u\right) = \pm2\sqrt{3} \int \sqrt{\frac{u^{2} + au}{2a^{3}u - 2au^{3} - u^{4} - 2a^{3}}} \, \mathrm{d} u.\f]

             The natural local tetrad is given by
             \f[ \mathbf{e}_{(0)} = \frac{1}{c}\partial_t,\quad \mathbf{e}_{(1)}=\partial_x,\quad \mathbf{e}_{(2)}= \frac{1}{p\left(u\right)}\partial_y,\quad \mathbf{e}_{(3)} = \frac{1}{p\left(u\right)}\partial_z.\f]

             Detailed discussions about the metric for a plane gravitational wave can be found
             in the standard literature, \ref lit_rindler "Rindler".
 */
// -------------------------------------------------------------------------------

#ifndef M4D_METRIC_PLANE_GRAV_WAVE_H
#define M4D_METRIC_PLANE_GRAV_WAVE_H

#include "m4dMetric.h"

namespace m4d {

// ---------------------------------------------------
//    class definition:   MetricPlaneGraveWave
// ---------------------------------------------------
class MetricPlaneGravWave: public Metric {
public:
    MetricPlaneGravWave(double longExt = 1.0, double degree = 30.0);
    virtual ~MetricPlaneGravWave();

// ----------- public methods ------------
public:
    virtual bool   calculateMetric(const double* pos);
    virtual bool   calculateChristoffels(const double* pos);
    virtual bool   calculateChrisD(const double* pos);

    virtual void   localToCoord(const double* pos, const double* ldir, double* dir,
                                enum_nat_tetrad_type type = enum_nat_tetrad_default);
    virtual void   coordToLocal(const double* pos, const double* cdir, double* ldir,
                                enum_nat_tetrad_type type = enum_nat_tetrad_default);

    virtual bool   breakCondition(const double* pos);

    virtual double testConstraint(const double y[], const double kappa);

    virtual bool   setParam(std::string pName, double val);
    virtual bool   transToTwoPlusOne(vec4 p, vec4 & cp);

    virtual bool   report(const vec4 pos, const vec4 cdir, std::string &text);

// --------- protected methods -----------
protected:
    virtual void setStandardValues();

// ----- specific protected methods ------
    double getValP(const double* pos);
    double getValQ(const double* pos);
    double getValDP(const double* pos);
    double getValDQ(const double* pos);
    double getValDDP(const double* pos);
    double getValDDQ(const double* pos);
    void calcFourierCoeff();

// --------- protected attribute ----------
protected:
    bool dataCalculated;
    double mDegree;
    double mLongExt;
    double intConst;
    double p0;
    double q0;
    std::vector<double> fCoeffB;

    double *uData, *dmData;
};

} // end namespace m4d

#endif
