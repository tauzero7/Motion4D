/**
 * @file    m4dMetricEinsteinRosenWaveWWB.cpp
 * @author  Thomas Mueller
 *
 * This file is part of the m4d-library.
 */
#include "metric/m4dMetricEinsteinRosenWaveWWB.h"

namespace m4d {

MetricEinsteinRosenWaveWWB::MetricEinsteinRosenWaveWWB(double c, double a)
{
    mMetricName = "EinsteinRosenWaveWWB";
    setCoordType(enum_coordinate_cylinder);

    mPhysicalUnits = enum_physical_constants_geom;
    mSpeedOfLight = 1.0;
    mGravConstant = 1.0;

    addParam("c", c);
    m_c = c;
    addParam("a", a);
    m_a = a;

    setStandardValues();
}

MetricEinsteinRosenWaveWWB::~MetricEinsteinRosenWaveWWB() {}

bool MetricEinsteinRosenWaveWWB::calculateMetric(const double* pos)
{
    double rho = pos[1];

    double gam, psi;
    calcPotentials(pos, gam, psi);

    double t1 = gam; // g(t,rho);
    double t2 = exp(t1);
    double t3 = t2 * t2;
    double t4 = psi; // psi(t,rho);
    double t5 = exp(t4);
    double t6 = t5 * t5;
    double t7 = 1 / t6;
    double t8 = t3 * t7;
    double t9 = rho * rho;

    g_compts[0][0] = -t8;
    g_compts[0][1] = 0.0;
    g_compts[0][2] = 0.0;
    g_compts[0][3] = 0.0;
    g_compts[1][0] = 0.0;
    g_compts[1][1] = t8;
    g_compts[1][2] = 0.0;
    g_compts[1][3] = 0.0;
    g_compts[2][0] = 0.0;
    g_compts[2][1] = 0.0;
    g_compts[2][2] = t9 * t7;
    g_compts[2][3] = 0.0;
    g_compts[3][0] = 0.0;
    g_compts[3][1] = 0.0;
    g_compts[3][2] = 0.0;
    g_compts[3][3] = t6;

    return true;
}

bool MetricEinsteinRosenWaveWWB::calculateChristoffels(const double* pos)
{
    double rho = pos[1];

    double gam, psi, gamt, gamr, psit, psir;
    calcDiffPoti(pos, gam, psi, gamt, gamr, psit, psir);

    double t1 = gamt; // diff(g(t,rho),t);
    double t2 = psit; // diff(psi(t,rho),t);
    double t3 = t1 - t2;
    double t4 = gamr; // diff(g(t,rho),rho);
    double t5 = psir; // diff(psi(t,rho),rho);
    double t6 = t4 - t5;
    double t9 = -1.0 + rho * t5;
    double t10 = t9 / rho;
    double t11 = gam; // g(t,rho);
    double t12 = exp(t11);
    double t13 = t12 * t12;
    double t14 = 1 / t13;
    double t15 = rho * rho;
    double t20 = psi; // psi(t,rho);
    double t21 = exp(t20);
    double t22 = t21 * t21;
    double t23 = t22 * t22;
    double t24 = t14 * t23;

    christoffel[0][0][0] = t3;
    christoffel[0][0][1] = t6;
    christoffel[0][0][2] = 0.0;
    christoffel[0][0][3] = 0.0;
    christoffel[0][1][0] = t6;
    christoffel[0][1][1] = t3;
    christoffel[0][1][2] = 0.0;
    christoffel[0][1][3] = 0.0;
    christoffel[0][2][0] = 0.0;
    christoffel[0][2][1] = 0.0;
    christoffel[0][2][2] = -t2;
    christoffel[0][2][3] = 0.0;
    christoffel[0][3][0] = 0.0;
    christoffel[0][3][1] = 0.0;
    christoffel[0][3][2] = 0.0;
    christoffel[0][3][3] = t2;
    christoffel[1][0][0] = t6;
    christoffel[1][0][1] = t3;
    christoffel[1][0][2] = 0.0;
    christoffel[1][0][3] = 0.0;
    christoffel[1][1][0] = t3;
    christoffel[1][1][1] = t6;
    christoffel[1][1][2] = 0.0;
    christoffel[1][1][3] = 0.0;
    christoffel[1][2][0] = 0.0;
    christoffel[1][2][1] = 0.0;
    christoffel[1][2][2] = -t10;
    christoffel[1][2][3] = 0.0;
    christoffel[1][3][0] = 0.0;
    christoffel[1][3][1] = 0.0;
    christoffel[1][3][2] = 0.0;
    christoffel[1][3][3] = t5;
    christoffel[2][0][0] = 0.0;
    christoffel[2][0][1] = 0.0;
    christoffel[2][0][2] = -t2;
    christoffel[2][0][3] = 0.0;
    christoffel[2][1][0] = 0.0;
    christoffel[2][1][1] = 0.0;
    christoffel[2][1][2] = -t10;
    christoffel[2][1][3] = 0.0;
    christoffel[2][2][0] = -t14 * t15 * t2;
    christoffel[2][2][1] = t14 * rho * t9;
    christoffel[2][2][2] = 0.0;
    christoffel[2][2][3] = 0.0;
    christoffel[2][3][0] = 0.0;
    christoffel[2][3][1] = 0.0;
    christoffel[2][3][2] = 0.0;
    christoffel[2][3][3] = 0.0;
    christoffel[3][0][0] = 0.0;
    christoffel[3][0][1] = 0.0;
    christoffel[3][0][2] = 0.0;
    christoffel[3][0][3] = t2;
    christoffel[3][1][0] = 0.0;
    christoffel[3][1][1] = 0.0;
    christoffel[3][1][2] = 0.0;
    christoffel[3][1][3] = t5;
    christoffel[3][2][0] = 0.0;
    christoffel[3][2][1] = 0.0;
    christoffel[3][2][2] = 0.0;
    christoffel[3][2][3] = 0.0;
    christoffel[3][3][0] = t24 * t2;
    christoffel[3][3][1] = -t24 * t5;
    christoffel[3][3][2] = 0.0;
    christoffel[3][3][3] = 0.0;

    return true;
}

bool MetricEinsteinRosenWaveWWB::calculateChrisD(const double*)
{
    return false;
}

void MetricEinsteinRosenWaveWWB::localToCoord(const double* pos, const double* ldir, double* dir, enum_nat_tetrad_type)
{

    double gam, psi;
    calcPotentials(pos, gam, psi);

    double rho = pos[1];
    double epg = exp(psi - gam);

    dir[0] = epg * ldir[0];
    dir[1] = epg * ldir[1];
    dir[2] = exp(psi) / rho * ldir[2];
    dir[3] = exp(-psi) * ldir[3];
}

void MetricEinsteinRosenWaveWWB::coordToLocal(const double* pos, const double* cdir, double* ldir, enum_nat_tetrad_type)
{
    double gam, psi;
    calcPotentials(pos, gam, psi);

    double rho = pos[1];
    double epg = exp(gam - psi);

    ldir[0] = epg * cdir[0];
    ldir[1] = epg * cdir[1];
    ldir[2] = rho * exp(-psi) * cdir[2];
    ldir[3] = exp(psi) * cdir[3];
}

bool MetricEinsteinRosenWaveWWB::breakCondition(const double*)
{
    return false;
}

double MetricEinsteinRosenWaveWWB::testConstraint(const double* y, const double kappa)
{
    double rho = y[1];
    double dt = y[4];
    double drho = y[5];
    double dphi = y[6];
    double dz = y[7];

    double gam, psi;
    calcPotentials(y, gam, psi);

    double sum = -kappa;
    sum += exp(2.0 * (gam - psi)) * (-dt * dt + drho * drho) + rho * rho * exp(-2.0 * psi) * dphi * dphi
        + exp(2.0 * psi) * dz * dz;
    return sum;
}

bool MetricEinsteinRosenWaveWWB::setParam(const char* pName, double val)
{
    Metric::setParam(pName, val);

    if (strcmp(pName, "c") == 0) {
        m_c = val;
    }
    else if (strcmp(pName, "a") == 0) {
        m_a = val;
    }
    return true;
}

bool MetricEinsteinRosenWaveWWB::report(const vec4, const vec4, char*& text)
{
    std::stringstream ss;
    ss << "Report for MetricEinsteinRosenWaveWWB metric\n\tcoordinates : (t,r,theta,phi)\n";
    ss << "---------------------------------------------------------------\n";
    ss << "  physical units ......... no\n";
    ss.precision(DEF_FIXED_REPORT_PRECISION);
    ss.setf(std::ios::fixed);
    text = new char[ss.str().length() + 2];
    return CopyString(ss.str().c_str(), text);
}

void MetricEinsteinRosenWaveWWB::setStandardValues()
{
    mInitPos[0] = 0.0;
    mInitPos[1] = 10.0;
    mInitPos[2] = 0.0;
    mInitPos[3] = 0.0;
    mInitDir[0] = 1.0;
    mInitDir[1] = 0.0;
    mInitDir[2] = 0.0;

    mCoordNames[0] = std::string("t");
    mCoordNames[1] = std::string("rho");
    mCoordNames[2] = std::string("phi");
    mCoordNames[3] = std::string("z");
}

void MetricEinsteinRosenWaveWWB::calcPotentials(const double* pos, double& gam, double& psi)
{
    double t = pos[0];
    double rho = pos[1];
    double a = m_a;
    double c = m_c;

    double t1 = sqrt(2.0);
    double t3 = a * a;
    double t4 = t3 * t3;
    double t5 = rho * rho;
    double t6 = t3 * t5;
    double t8 = t * t;
    double t9 = t3 * t8;
    double t11 = t5 * t5;
    double t14 = t8 * t8;
    double t16 = sqrt(t4 + 2.0 * t6 + 2.0 * t9 + t11 - 2.0 * t5 * t8 + t14);
    double t19 = pow(t3 + t5 - t8, 2.0);
    double t20 = 4.0 * t9;
    double t21 = t19 + t20;
    double t24 = sqrt((t16 + t3 + t5 - t8) / t21);
    double t26 = c * c;
    double t30 = t21 * t21;

    psi = t1 * c * t24;
    gam = t26 / t3 * (1.0 - 2.0 * t6 * (t19 - t20) / t30 + (t5 - t3 - t8) / t16) / 2.0;
}

void MetricEinsteinRosenWaveWWB::calcDiffPoti(
    const double* pos, double& gam, double& psi, double& gamt, double& gamr, double& psit, double& psir)
{
    double t = pos[0];
    double rho = pos[1];

    double a = m_a;
    double c = m_c;
    calcPotentials(pos, gam, psi);

    double t1 = sqrt(2.0);
    double t2 = t1 * c;
    double t3 = a * a;
    double t4 = t3 * t3;
    double t5 = rho * rho;
    double t6 = t3 * t5;
    double t8 = t * t;
    double t9 = t3 * t8;
    double t11 = t5 * t5;
    double t14 = t8 * t8;
    double t15 = t4 + 2.0 * t6 + 2.0 * t9 + t11 - 2.0 * t5 * t8 + t14;
    double t16 = sqrt(t15);
    double t17 = t16 + t3 + t5 - t8;
    double t18 = t3 + t5 - t8;
    double t19 = t18 * t18;
    double t20 = 4.0 * t9;
    double t21 = t19 + t20;
    double t22 = 1 / t21;
    double t24 = sqrt(t17 * t22);
    double t25 = 1 / t24;
    double t26 = 1 / t16;
    double t27 = t3 * t;
    double t30 = t27 - t5 * t + t8 * t;
    double t36 = t21 * t21;
    double t37 = 1 / t36;
    double t38 = t17 * t37;
    double t40 = 4.0 * t18 * t;
    double t41 = 8.0 * t27;
    double t42 = -t40 + t41;
    double t48 = t3 * rho;
    double t49 = t5 * rho;
    double t51 = t48 + t49 - rho * t8;
    double t64 = c * c;
    double t66 = t64 / t3;
    double t71 = t19 - t20;
    double t74 = t71 / t36 / t21;
    double t83 = (t5 - t3 - t8) / t16 / t15;
    double t92 = t3 * t49;

    psit = t2 * t25 * ((2.0 * t26 * t30 - 2.0 * t) * t22 - t38 * t42) / 2.0;
    psir = t2 * t25 * ((2.0 * t26 * t51 + 2.0 * rho) * t22 - 4.0 * t38 * t18 * rho) / 2.0;
    gamt = t66 * (-2.0 * t6 * (-t40 - t41) * t37 + 4.0 * t6 * t74 * t42 - 2.0 * t * t26 - 2.0 * t83 * t30) / 2.0;
    gamr = t66
        * (-4.0 * t48 * t71 * t37 - 8.0 * t92 * t18 * t37 + 16.0 * t92 * t74 * t18 + 2.0 * rho * t26 - 2.0 * t83 * t51)
        / 2.0;
}

} // end namespace m4d
