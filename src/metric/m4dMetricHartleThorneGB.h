/**
 * @file    m4dMetricHartleThorneGB.h
 * @author  Thomas Mueller
 *
 * @brief  Hartle-Thorne metric in spherical coordinates (t,r,theta,phi).

     The line element is given by

     \f[ds^2 = -dt^2 + .\f]


     Hartle-Thorne metric following

     Kostas Glampedakis and Stanislav Babak,
     "Mapping spacetimes with LISA: inspiral of a test body in a ‘quasi-Kerr’ field",
     Class. Quantum Grav. 23 (2006) 4167–4188

     see also ApJ 753,175 (2012)

 * This file is part of the m4d-library.
 */
#ifndef M4D_METRIC_HARTLE_THORNE_GB_H
#define M4D_METRIC_HARTLE_THORNE_GB_H

#include "m4dMetric.h"

namespace m4d {

// ---------------------------------------------------
//    class definition:   MetricHartleThorneGB
// ---------------------------------------------------
class MetricHartleThorneGB : public Metric
{
public:
    MetricHartleThorneGB(double mass = 1.0, double angmom = 0.0, double eta = 0.0);
    virtual ~MetricHartleThorneGB();

public:
    virtual bool calculateMetric(const double* pos);
    virtual bool calculateChristoffels(const double* pos);
    virtual bool calculateChrisD(const double* pos);

    virtual void localToCoord(
        const double* pos, const double* ldir, double* dir, enum_nat_tetrad_type type = enum_nat_tetrad_default);
    virtual void coordToLocal(
        const double* pos, const double* cdir, double* ldir, enum_nat_tetrad_type type = enum_nat_tetrad_default);

    virtual bool breakCondition(const double* pos);

    virtual bool setParam(const char* pName, double val);

    virtual bool report(const vec4 pos, const vec4 cdir, char*& text);

public:
    void calcFunc(const double* pos);
    void calcFuncDiff(const double* pos);

protected:
    virtual void setStandardValues();

protected:
    double mMass, mAngmom, mEta;
    double sigma, delta;
    double dsigmadr, dsigmadth, deltadr;
    double htt, hrr, hthth, hphph;
    double dhttdr, dhttdth, dhrrdr, dhrrdth;
    double dhaadr, dhaadth, dhppdr, dhppdth;
};

} // end namespace m4d

#endif
