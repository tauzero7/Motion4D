// -------------------------------------------------------------------------------
/*
    m4dMetricHartleThorneGB.h

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

/*!  \class  m4d::MetricGlampedakis
     \brief  Hartle-Thorne metric in spherical coordinates (t,r,theta,phi).

             The line element is given by

             \f[ds^2 = -dt^2 + .\f]


             Hartle-Thorne metric following

             Kostas Glampedakis and Stanislav Babak,
             "Mapping spacetimes with LISA: inspiral of a test body in a ‘quasi-Kerr’ field",
             Class. Quantum Grav. 23 (2006) 4167–4188

             see also ApJ 753,175 (2012)
*/
// -------------------------------------------------------------------------------

#ifndef M4D_METRIC_GLAMPEDAKIS_H
#define M4D_METRIC_GLAMPEDAKIS_H

#include "m4dMetric.h"

namespace m4d {

// ---------------------------------------------------
//    class definition:   MetricGlampedakis
// ---------------------------------------------------
class MetricGlampedakis : public Metric {
public:
    MetricGlampedakis(double mass = 1.0, double angmom = 0.0, double epsilon = 0.0);
    virtual ~MetricGlampedakis();

public:
    virtual bool   calculateMetric(const double* pos);
    virtual bool   calculateChristoffels(const double* pos);
    virtual bool   calculateChrisD(const double* pos);
    virtual bool   calculateRiemann(const double* pos);

    virtual void   localToCoord(const double* pos, const double* ldir, double* dir,
                                enum_nat_tetrad_type  type = enum_nat_tetrad_default);
    virtual void   coordToLocal(const double* pos, const double* cdir, double* ldir,
                                enum_nat_tetrad_type  type = enum_nat_tetrad_default);

    virtual bool   breakCondition(const double* pos);

    virtual bool   setParam(const char* pName, double val);

    virtual bool   report(const vec4 pos, const vec4 cdir, std::string &text);

public:
    void  calcGlampedakis(const double* pos);
    void  calcGlampedakisDiff(const double* pos);
    void  calcGlampedakisDiff2(const double* pos);

    void calcKerr(const double* pos);
    void calcKerrDiff(const double* pos);
    void calcKerrDiff2(const double* pos);

    void calcgComps(const double* pos);
    void calcgCompsDiff(const double* pos);
    void calcgCompsDiff2(const double* pos);

protected:
    virtual void setStandardValues();
    void initToZero();

protected:
    double  mMass, mAngmom, mEpsilon;
    double rs, a, a2;

    double r, r2, r3, r4, r5, r6;
    double M, M2, M3, M4, M5, M6;
    double Mr;
    // for Glampedakis
    double cf, fSchw, lgr, ff1, ff2;
    double cfdth, fSchwdr;
    double cfd2th, fSchwd2r;

    double theta;
    double sth, cth; // sin(theta)
    //double cth = cos(theta);
    double sth2, cth2; //squared
    double sth3, sth4;
    double cth3, cth4;

    // Kerr parameter and derivatives
    double sigma, delta;
    double sigma2, delta2; // sigma*sigma, delta*delta
    double sigma3, delta3;
    double sigma2Inv; // 1/sigma2

    double sigmadr, sigmadth, deltadr;
    double sigmad2r, deltad2r, sigmad2th;
    //Kerr metric components
    double ktt, ktph, krr, kthth, kphph;

    double ktt2, ktph2, krr2, kthth2, kphph2;

    // their derivatives
    double kttdr, kttdth;
    double ktphdr, ktphdth;
    double krrdr, krrdth;
    double kththdr, kththdth;
    double kphphdr, kphphdth;

    //their second derivatives
    double kttd2r, kttd2th, kttdrdth;
    double ktphd2r, ktphd2th, ktphdrdth;
    double krrd2r, krrd2th, krrdrdth;
    double kththd2r, kththd2th, kththdrdth;
    double kphphd2r, kphphd2th, kphphdrdth;


    // additional parameters by Glampedakis
    double htt, hrr, hthth, hphph;
    // their derivatives
    double httdr, httdth, hrrdr, hrrdth;
    double hththdr, hththdth, hphphdr, hphphdth;
    //their second derivatives
    double httd2r, httd2th, httdrdth;
    double hrrd2r, hrrd2th, hrrdrdth;
    double hththd2r, hththd2th, hththdrdth;
    double hphphd2r, hphphd2th, hphphdrdth;

    // the metric corrections
    double gltt, gltph, glrr, glthth, glphph;

    // their derivatives
    double glttdr, glttdth;
    double gltphdr, gltphdth;
    double glrrdr, glrrdth;
    double glththdr, glththdth;
    double glphphdr, glphphdth;


    // their second derivatives
    double glttd2r, glttd2th, glttdrdth;
    double gltphd2r, gltphd2th, gltphdrdth;
    double glrrd2r, glrrd2th, glrrdrdth;
    double glththd2r, glththd2th, glththdrdth;
    double glphphd2r, glphphd2th, glphphdrdth;

    // functions defined by glampedakis
    double F1, F2;
    // their derivatives
    double F1dr, F2dr;
    //their second derivatives
    double F1d2r, F2d2r;

    //the full metric components
    double gtt, gtph, grr, gthth, gphph;
    // their derivatives
    double gttdr, gttdth;
    double gtphdr, gtphdth;
    double grrdr, grrdth;
    double gththdr, gththdth;
    double gphphdr, gphphdth;

    // their second derivatives
    double gttd2r, gttd2th, gttdrdth;
    double gtphd2r, gtphd2th, gtphdrdth;
    double grrd2r, grrd2th, grrdrdth;
    double gththd2r, gththd2th, gththdrdth;
    double gphphd2r, gphphd2th, gphphdrdth;

    // positions where the corresponding functions were last evaluated
    // only r and theta coordinates are relevant for this metric
    double compoldpos[2];
    double compdiffoldpos[2];
    double oldchristoffelpos[2];
    double oldchrisDpos[2];
    double oldglamppos[2];
    double oldglampdiffpos[2];
    double oldglampdiff2pos[2];
    double oldkerrpos[2];
    double oldkerrdiffpos[2];
    double oldkerrdiff2pos[2];
    double oldmetricpos[2];

    int metricskip;
    int metriccalled;

    int christoffelskip;
    int christoffelcalled;

    int chrisDskip;
    int chrisDcalled;

    int kerrskip;
    int kerrcalled;

    int kerrdiffskip;
    int kerrdiffcalled;

    int kerrdiff2skip;
    int kerrdiff2called;

    int glampskip;
    int glampcalled;

    int glampdiffskip;
    int glampdiffcalled;

    int glampdiff2skip;
    int glampdiff2called;

};

} // end namespace m4d

#endif

