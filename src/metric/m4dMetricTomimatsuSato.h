// -------------------------------------------------------------------------------
/*
    m4dMetricTomimatsuSato.h

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

/*!  \class  m4d::MetricGTomimatsuSato
     \brief
             V.S. Manko: Progress of Theoretical Physics 127, 1057 (2012)
*/
// -------------------------------------------------------------------------------

#ifndef M4D_METRIC_TOMIMATSUSATO_H
#define M4D_METRIC_TOMIMATSUSATO_H

#include "m4dMetric.h"

namespace m4d {

// ---------------------------------------------------
//    class definition:   MetricTomimatsuSato
// ---------------------------------------------------
class MetricTomimatsuSato : public Metric {
public:
    // standard values such that M=1
    MetricTomimatsuSato(double k0 = 0.4330127019, double p0 = 0.5, double q0 = 0.8660254038);
    virtual ~MetricTomimatsuSato();

public:
    virtual bool   calculateMetric(const double* pos);
    virtual bool   calculateChristoffels(const double* pos);
    virtual bool   calculateChrisD(const double* pos);

    virtual void   localToCoord(const double* pos, const double* ldir, double* dir,
                                enum_nat_tetrad_type  type = enum_nat_tetrad_default);
    virtual void   coordToLocal(const double* pos, const double* cdir, double* ldir,
                                enum_nat_tetrad_type  type = enum_nat_tetrad_default);

    virtual bool   breakCondition(const double* pos);

    virtual bool   setParam(std::string pName, double val);

    virtual bool   report(const vec4 pos, const vec4 cdir, std::string &text);

public:
    void  calcTomimatsuSato(const double* pos);
    void  calcTomimatsuSatoDiff(const double* pos);
    void  calcTomimatsuSatoDiff2(const double* pos);

protected:
    virtual void setStandardValues();
    void initToZero();

protected:
    double  p, q, k;
    double mp, mq, mk;
    double p2, q2, k2; // squared


    // the metric components in grtensor notation
    double g44, g34, g11, g22, g33;


    // their derivatives
    double g44dx, g44dy;
    double g34dx, g34dy;
    double g11dx, g11dy;
    double g22dx, g22dy;
    double g33dx, g33dy;


    // their second derivatives
    double g44d2x, g44d2y, g44dxdy;
    double g34d2x, g34d2y, g34dxdy;
    double g11d2x, g11d2y, g11dxdy;
    double g22d2x, g22d2y, g22dxdy;
    double g33d2x, g33d2y, g33dxdy;



    // functions defined by V.S. Manko with some additional definitions
    double f, g, omega, A, B; // g=exp(2*gamma) in the paper, A, B not as in the paper, A is the nominator of f,g, denominator of omega
    double mu, sigma, nu, tau;

    double mu2, sigma2, omega2;


    // their derivatives
    double fdx, fdy;
    double gdx, gdy;
    double omegadx, omegady;
    double Adx, Ady;
    double mudx, mudy;
    double sigmadx, sigmady;
    double nudx; // vanishes: nudy
    double taudx, taudy;


    // their second derivatives
    double fd2x, fd2y, fdxdy;
    double gd2x, gd2y, gdxdy;
    double omegad2x, omegad2y, omegadxdy;
    double Ad2x, Ad2y, Adxdy;
    double mud2x, mud2y; // vanish: mudxdy
    double sigmad2x, sigmad2y; //, vanish: sigmadxdy
    double nud2x; // vanish: nud2y, nudxdy
    double taud2y, taudxdy; // vanish: taud2x

    double x, y;
    double x2, y2; //x*x, y*y
    double x4, y4;
    double omy2; // 1-y^2
    double x2mo; // x^2 - 1
    double x2my2; //x^2 - y^2




    // positions where the corresponding functions were last evaluated
    // only r and theta coordinates are relevant for this metric
    double oldcomppos[2];
    double oldcompdiffpos[2];
    double oldcompdiff2pos[2];

    double oldmetricpos[2];
    double oldchristoffelpos[2];
    double oldchrisDpos[2];

    int metricskip;
    int metriccalled;

    int christoffelskip;
    int christoffelcalled;

    int chrisDskip;
    int chrisDcalled;

    int compskip;
    int compcalled;

    int compdiffskip;
    int compdiffcalled;

    int compdiff2skip;
    int compdiff2called;

};

} // end namespace m4d

#endif

