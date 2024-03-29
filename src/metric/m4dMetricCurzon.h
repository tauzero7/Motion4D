/**
 * @file    m4dMetricCurzon.h
 * @author  Thomas Mueller
 *
 * @brief  Curzon metric in cylindrical coordinates (t,r,phi,z).

     The line element is given by

     \f[ds^2 = ,\f]
     where

     The natural local tetrad reads:


     Detailed discussions about the Curzon metric can be found in
     <ul>
       <li> Susan M. Scott and P. Szekeres, "The Curzon Singularity. I: Spatial Sections," Gen. Relativ. Gravit.
<b>18</b>, 557--570 (1986).</li>
     </ul>

 *
 * This file is part of the m4d-library.
 */
#ifndef M4D_METRIC_CURZON_H
#define M4D_METRIC_CURZON_H

#include "m4dMetric.h"

namespace m4d {

// ---------------------------------------------------
//    class definition:   MetricCurzon
// ---------------------------------------------------
class MetricCurzon : public Metric
{
public:
    MetricCurzon(double mass = 1.0);
    virtual ~MetricCurzon();

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

    virtual double testConstraint(const double y[], const double kappa);

    virtual bool setParam(const char* pName, double val);

    virtual bool report(const vec4 pos, const vec4 cdir, char*& text);

    // --------- protected methods -----------
protected:
    virtual void setStandardValues();

    double calcR(const double* pos);
    void calcDR(const double* pos, double& R, double& dRdr, double& dRdz);
    void calcD2R(
        const double* pos, double& R, double& dRdr, double& dRdz, double& dRdrdr, double& dRdrdz, double& dRdzdz);

    void calcLambdaNu(const double* pos, double& lambda, double& nu);
    void calcDLN(const double* pos, double& lambda, double& nu, double& dldr, double& dldz, double& dndr, double& dndz);
    void calcD2LN(const double* pos, double& lambda, double& nu, double& dldr, double& dldz, double& dndr, double& dndz,
        double& dldrdr, double& dldrdz, double& dldzdz, double& dndrdr, double& dndrdz, double& dndzdz);

    // -------- protected attribute ---------
protected:
    double mMass;
};

} // end namespace m4d

#endif
