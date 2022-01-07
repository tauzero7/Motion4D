/**
 * @file    m4dMetricKastorTraschen.cpp
 * @author  Thomas Mueller
 *
 * This file is part of the m4d-library.
 */
#include "m4dMetricKastorTraschen.h"

namespace m4d {

MetricKastorTraschen::MetricKastorTraschen(double H)
{
    mMetricName = "KastorTraschen";
    setCoordType(enum_coordinate_cartesian);

    mPhysicalUnits = enum_physical_constants_geom;
    mSpeedOfLight = 1.0;
    mGravConstant = 1.0;

    /*  Only a static tetrad is defined  */
    mLocTeds.push_back(enum_nat_tetrad_static);

    m1 = 1.0;
    z1 = 1.0;
    m2 = 1.0;
    z2 = -1.0;
    mH = H;

    addParam("m1", m1);
    addParam("z1", z1);
    addParam("m2", m2);
    addParam("z2", z2);
    addParam("h", mH);

    setStandardValues();
}

MetricKastorTraschen::~MetricKastorTraschen() {}

bool MetricKastorTraschen::calculateMetric(const double* pos)
{
    double Omega, a;
    calcPotentials(pos, Omega, a);

    double t1 = Omega; // Omega(t,x,y,z);
    double t2 = t1 * t1;
    double t4 = a; // a(t);
    double t5 = t4 * t4;
    double t6 = t2 * t5;

    g_compts[0][0] = -1 / t2;
    g_compts[0][1] = 0.0;
    g_compts[0][2] = 0.0;
    g_compts[0][3] = 0.0;
    g_compts[1][0] = 0.0;
    g_compts[1][1] = t6;
    g_compts[1][2] = 0.0;
    g_compts[1][3] = 0.0;
    g_compts[2][0] = 0.0;
    g_compts[2][1] = 0.0;
    g_compts[2][2] = t6;
    g_compts[2][3] = 0.0;
    g_compts[3][0] = 0.0;
    g_compts[3][1] = 0.0;
    g_compts[3][2] = 0.0;
    g_compts[3][3] = t6;

    return true;
}

bool MetricKastorTraschen::calculateChristoffels(const double* pos)
{
    double Omega, a;
    double dOdx, dOdy, dOdz, dOdt, dadt;
    calcPotDiffs(pos, Omega, a, dOdx, dOdy, dOdz, dOdt, dadt);

    double t1 = Omega; // Omega(t,x,y,z);
    double t2 = 1 / t1;
    double t3 = dOdt; // diff(Omega(t,x,y,z),t);
    double t5 = t1 * t1;
    double t6 = t5 * t5;
    double t9 = a; // a(t);
    double t10 = t9 * t9;
    double t12 = 1 / t6 / t1 / t10;
    double t13 = dOdx; // diff(Omega(t,x,y,z),x);
    double t15 = dOdy; // diff(Omega(t,x,y,z),y);
    double t17 = dOdz; // diff(Omega(t,x,y,z),z);
    double t19 = t2 * t13;
    double t23 = dadt; // diff(a(t),t);
    double t25 = t9 * t3 + t1 * t23;
    double t26 = t2 / t9 * t25;
    double t27 = t2 * t15;
    double t28 = t2 * t17;
    double t31 = t5 * t1 * t9 * t25;

    christoffel[0][0][0] = -t2 * t3;
    christoffel[0][0][1] = -t12 * t13;
    christoffel[0][0][2] = -t12 * t15;
    christoffel[0][0][3] = -t12 * t17;
    christoffel[0][1][0] = -t19;
    christoffel[0][1][1] = t26;
    christoffel[0][1][2] = 0.0;
    christoffel[0][1][3] = 0.0;
    christoffel[0][2][0] = -t27;
    christoffel[0][2][1] = 0.0;
    christoffel[0][2][2] = t26;
    christoffel[0][2][3] = 0.0;
    christoffel[0][3][0] = -t28;
    christoffel[0][3][1] = 0.0;
    christoffel[0][3][2] = 0.0;
    christoffel[0][3][3] = t26;
    christoffel[1][0][0] = -t19;
    christoffel[1][0][1] = t26;
    christoffel[1][0][2] = 0.0;
    christoffel[1][0][3] = 0.0;
    christoffel[1][1][0] = t31;
    christoffel[1][1][1] = t19;
    christoffel[1][1][2] = -t27;
    christoffel[1][1][3] = -t28;
    christoffel[1][2][0] = 0.0;
    christoffel[1][2][1] = t27;
    christoffel[1][2][2] = t19;
    christoffel[1][2][3] = 0.0;
    christoffel[1][3][0] = 0.0;
    christoffel[1][3][1] = t28;
    christoffel[1][3][2] = 0.0;
    christoffel[1][3][3] = t19;
    christoffel[2][0][0] = -t27;
    christoffel[2][0][1] = 0.0;
    christoffel[2][0][2] = t26;
    christoffel[2][0][3] = 0.0;
    christoffel[2][1][0] = 0.0;
    christoffel[2][1][1] = t27;
    christoffel[2][1][2] = t19;
    christoffel[2][1][3] = 0.0;
    christoffel[2][2][0] = t31;
    christoffel[2][2][1] = -t19;
    christoffel[2][2][2] = t27;
    christoffel[2][2][3] = -t28;
    christoffel[2][3][0] = 0.0;
    christoffel[2][3][1] = 0.0;
    christoffel[2][3][2] = t28;
    christoffel[2][3][3] = t27;
    christoffel[3][0][0] = -t28;
    christoffel[3][0][1] = 0.0;
    christoffel[3][0][2] = 0.0;
    christoffel[3][0][3] = t26;
    christoffel[3][1][0] = 0.0;
    christoffel[3][1][1] = t28;
    christoffel[3][1][2] = 0.0;
    christoffel[3][1][3] = t19;
    christoffel[3][2][0] = 0.0;
    christoffel[3][2][1] = 0.0;
    christoffel[3][2][2] = t28;
    christoffel[3][2][3] = t27;
    christoffel[3][3][0] = t31;
    christoffel[3][3][1] = -t19;
    christoffel[3][3][2] = -t27;
    christoffel[3][3][3] = t28;

    return true;
}

bool MetricKastorTraschen::calculateChrisD(const double*)
{
    return false;
}

void MetricKastorTraschen::localToCoord(const double* pos, const double* ldir, double* dir, enum_nat_tetrad_type)
{
    double Omega, a;
    calcPotentials(pos, Omega, a);

    double edOa = 1.0 / (Omega * a);
    dir[0] = Omega * ldir[0];
    dir[1] = edOa * ldir[1];
    dir[2] = edOa * ldir[2];
    dir[3] = edOa * ldir[3];
}

void MetricKastorTraschen::coordToLocal(const double* pos, const double* cdir, double* ldir, enum_nat_tetrad_type)
{
    double Omega, a;
    calcPotentials(pos, Omega, a);

    double Oa = Omega * a;

    ldir[0] = cdir[0] / Omega;
    ldir[1] = cdir[1] * Oa;
    ldir[2] = cdir[2] * Oa;
    ldir[3] = cdir[3] * Oa;
}

bool MetricKastorTraschen::breakCondition(const double*)
{
    bool br = false;
    return br;
}

double MetricKastorTraschen::testConstraint(const double y[], const double kappa)
{
    double Omega, a;
    calcPotentials(y, Omega, a);

    double dt = y[4];
    double dx = y[5];
    double dy = y[6];
    double dz = y[7];

    double sum = -kappa;
    sum += -dt * dt / (Omega * Omega) + Omega * Omega * a * a * (dx * dx + dy * dy + dz * dz);
    // std::cerr << Omega << " " << a << " " << sum << std::endl;
    return sum;
}

bool MetricKastorTraschen::setParam(const char* pName, double val)
{
    Metric::setParam(pName, val);
    if (strcmp(pName, "m1") == 0) {
        m1 = val;
    }
    else if (strcmp(pName, "z1") == 0) {
        z1 = val;
    }
    else if (strcmp(pName, "m2") == 0) {
        m2 = val;
    }
    else if (strcmp(pName, "z2") == 0) {
        z2 = val;
    }
    else if (strcmp(pName, "h") == 0) {
        mH = val;
    }
    return true;
}

bool MetricKastorTraschen::report(const vec4, const vec4, char*& text)
{
    std::stringstream ss;
    ss << "Report for Kastor-Traschen metric\n\tcoordinate : (t,x,y,z)\n";
    ss << "---------------------------------------------------------------\n";
    ss << "  physical units ..................... no\n";
    ss.precision(DEF_FIXED_REPORT_PRECISION);
    ss.setf(std::ios::fixed);

    text = new char[ss.str().length() + 2];
    return CopyString(ss.str().c_str(), text);
}

void MetricKastorTraschen::calcPotentials(const double* pos, double& Omega, double& a)
{
    double t = pos[0];
    double x = pos[1];
    double y = pos[2];
    double z = pos[3];

    double r1 = sqrt(x * x + y * y + (z - z1) * (z - z1));
    double r2 = sqrt(x * x + y * y + (z - z2) * (z - z2));

    a = exp(mH * t);
    Omega = 1.0 + m1 / (r1 * a) + m2 / (r2 * a);
}

void MetricKastorTraschen::calcPotDiffs(
    const double* pos, double& Omega, double& a, double& dOdx, double& dOdy, double& dOdz, double& dOdt, double& dadt)
{
    double t = pos[0];
    double x = pos[1];
    double y = pos[2];
    double z = pos[3];

    double r1 = sqrt(x * x + y * y + (z - z1) * (z - z1));
    double r2 = sqrt(x * x + y * y + (z - z2) * (z - z2));

    a = exp(mH * t);
    Omega = 1.0 + m1 / (r1 * a) + m2 / (r2 * a);

    double dr1dx = x / r1;
    double dr1dy = y / r1;
    double dr1dz = (z - z1) / r1;

    double dr2dx = x / r2;
    double dr2dy = y / r2;
    double dr2dz = (z - z2) / r2;

    dOdx = -m1 / (r1 * r1 * a) * dr1dx - m2 / (r2 * r2 * a) * dr2dx;
    dOdy = -m1 / (r1 * r1 * a) * dr1dy - m2 / (r2 * r2 * a) * dr2dy;
    dOdz = -m1 / (r1 * r1 * a) * dr1dz - m2 / (r2 * r2 * a) * dr2dz;

    dadt = mH * exp(mH * t);
    dOdt = -m1 * dadt / (r1 * a * a) - m2 * dadt / (r2 * a * a);
}

// ********************************* protected methods *****************************
void MetricKastorTraschen::setStandardValues()
{
    mInitPos[0] = 0.0;
    mInitPos[1] = 30.0;
    mInitPos[2] = 0.0;
    mInitPos[3] = 0.0;
    mInitDir[0] = 1.0;
    mInitDir[1] = 0.0;
    mInitDir[2] = 0.0;

    mCoordNames[0] = std::string("t");
    mCoordNames[1] = std::string("x");
    mCoordNames[2] = std::string("y");
    mCoordNames[3] = std::string("z");
}

} // end namespace m4d
