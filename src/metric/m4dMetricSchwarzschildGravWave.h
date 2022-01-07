/**
 * @file    m4dMetricSchwarzschildGravWave.h
 * @author  Thomas Mueller
 *
 * @brief  Schwarzschild gravitational wave metric in spherical coordinates (t,r,theta,phi)
             with perturbations

 *
 * This file is part of the m4d-library.
 */
#ifndef M4D_METRIC_SCHWARZSCHILD_GRAVWAVE_H
#define M4D_METRIC_SCHWARZSCHILD_GRAVWAVE_H

#include "m4dMetric.h"

namespace m4d {

// ---------------------------------------------------
//    class definition:   MetricSchwarzschildGravWave
// ---------------------------------------------------
class MetricSchwarzschildGravWave : public Metric
{
public:
    MetricSchwarzschildGravWave(double mass = 1.0);
    virtual ~MetricSchwarzschildGravWave();

    // --------- public methods -----------
public:
    /**
     * @brief Calculate the contravariant metric components at position 'pos'.
     * @param pos  Position in metric coordinates where coefficients have to be evaluated.
     * @return true if calculation succeeded.
     */
    virtual bool calculateMetric(const double* pos);
    virtual bool calculateChristoffels(const double* pos);
    // virtual bool   calculateChrisD(const double* pos);

    virtual void localToCoord(
        const double* pos, const double* ldir, double* dir, enum_nat_tetrad_type type = enum_nat_tetrad_default);
    virtual void coordToLocal(
        const double* pos, const double* cdir, double* ldir, enum_nat_tetrad_type type = enum_nat_tetrad_default);

    virtual bool breakCondition(const double* pos);

    virtual double testConstraint(const double y[], const double kappa);

    virtual bool setParam(const char* pName, double val);

    virtual bool report(const vec4 pos, const vec4 cdir, char*& text);

    // --------- protected methods -----------
protected:
    virtual void setStandardValues();

    void calcPerturbations(const double* pos);
    void calcPerturbationsAndDiffs(const double* pos);

    void calcLegendre(int l, double theta, double& Pl, double& Plt, double& Pltt, double& Plttt);

    // -------- protected attribute ---------
protected:
    //! Schwarzschild radius rs=2GM/c^2.
    double rs;
    double mMass;
    double mEpsilon;
    double mSigma;
    int ml;

    double htt, hrr, hee, hpp;
    double htt_t, htt_r, htt_theta;
    double hrr_t, hrr_r, hrr_theta;
    double hee_t, hee_r, hee_theta;
    double hpp_t, hpp_r, hpp_theta;
};

} // end namespace m4d

#endif
