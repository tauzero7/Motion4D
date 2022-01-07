/**
 * @file    m4dMetricBonnor.h
 * @author  Thomas Mueller
 *
 * @brief  Bonnor metric in spherical coordinates (t,r,phi,z).

     The line element is given by

     \f[ds^2 = -\left(\frac{P}{Y}\right)^2 dt^2 + \frac{Y^2P^2}{Q^3Z}\left[dr^2+Z\,d\vartheta^2\right] +
Z\left(\frac{Y}{P}\right)^2\sin^2\vartheta\,d\varphi^2,\f] where P=P(r,theta),...


     The natural local tetrad reads:


     Detailed discussions about the Bonnor metric can be found in
     <ul>
       <li> W.B. Bonnor, "An exact Solution of the Einstein-Maxwell Equations Referring to a Magnetic Dipole,"
Z. Phys. <b>190</b>, 444--445 (1966).</li>
     </ul>

 * This file is part of the m4d-library.
 *
 */
#ifndef M4D_METRIC_BONNOR_H
#define M4D_METRIC_BONNOR_H

#include "m4dMetric.h"

namespace m4d {

// ---------------------------------------------------
//    class definition:   MetricBonnor
// ---------------------------------------------------
class MetricBonnor : public Metric
{
public:
    MetricBonnor(double mass = 1.0, double b = 0.1);
    virtual ~MetricBonnor();

    // --------- public methods -----------
public:
    virtual bool calculateMetric(const double* pos);
    virtual bool calculateChristoffels(const double* pos);
    virtual bool calculateChrisD(const double* pos);

    virtual void localToCoord(
        const double* pos, const double* ldir, double* dir, enum_nat_tetrad_type type = enum_nat_tetrad_cylinder);
    virtual void coordToLocal(
        const double* pos, const double* cdir, double* ldir, enum_nat_tetrad_type type = enum_nat_tetrad_cylinder);

    virtual bool breakCondition(const double* pos);

    // virtual double testConstraint         ( const double y[], const double kappa );

    virtual bool setParam(const char* pName, double val);

    virtual bool report(const vec4 pos, const vec4 cdir, char*& text);

    // --------- protected methods -----------
protected:
    virtual void setStandardValues();
    void calcPotentials(const double* pos);
    void calcPotiDiffs(const double* pos);

    // -------- protected attribute ---------
protected:
    double mMass;
    double mB;
    double P, Q, Y, Z;
    double dPdr, dPdtheta, dQdr, dQdtheta, dYdr, dYdtheta, dZdr;
};

} // end namespace m4d

#endif
