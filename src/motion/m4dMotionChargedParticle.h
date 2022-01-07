/**
 * @file    m4dMotionChargedParticle.h
 * @author  Thomas Mueller
 *
 * @brief  Motion of a charged particle
 *
 *
 * This file is part of the m4d-library.
 */
#ifndef M4D_MOTION_CHARGED_PARTICLE_H
#define M4D_MOTION_CHARGED_PARTICLE_H

#include <iostream>

#include "m4dMotion.h"

namespace m4d {

class API_M4D_EXPORT MotionChargedParticle : public Motion
{
public:
    explicit MotionChargedParticle(Metric* metric);
    virtual ~MotionChargedParticle();

public:
    bool setInitialVelocity(double fm, double v, double theta, double phi, double q_over_m,
        enum_nat_tetrad_type type = enum_nat_tetrad_default);
    void getInitialVelocity(double& v, double& theta, double& phi);

    virtual enum_break_condition calculateMotion(const vec4 initPos, double fm, double v, double theta_v, double phi_v,
        const double q_over_m, const vec4 e0, const vec4 e1, const vec4 e2, const vec4 e3, const int maxNumPoints,
        std::vector<vec4>& points, std::vector<vec4>& base0, std::vector<vec4>& base1, std::vector<vec4>& base2,
        std::vector<vec4>& base3);

    enum_break_condition nextStep();

protected:
    bool calcDerivs(const double y[], double dydx[]);

protected:
    //! Current four-velocity in coordinates.
    double mInitVel[4];

    double mVel;
    double mTheta;
    double mPhi;

    double mQoverM;
};

} // end namespace m4d

#endif
