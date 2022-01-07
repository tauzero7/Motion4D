/**
 * @file    m4dMetricSchwarzschildWT.cpp
 * @author  Thomas Mueller
 *
 * This file is part of the m4d-library.
 */
#include "m4dMetricSchwarzschildWT.h"
#include <cmath>

#define sign(x) ((x >= 0) ? 1.0 : -1.0)

namespace m4d {

MetricSchwarzschildWT::MetricSchwarzschildWT(double mass)
{
    mMetricName = "SchwarzschildWT";
    mMetricCPPfilename = "m4dMetricSchwarzschildWT.cpp";
    setCoordType(enum_coordinate_spherical);

    mPhysicalUnits = enum_physical_constants_geom;
    mSpeedOfLight = 1.0;
    mGravConstant = 1.0;

    addParam("mass", mass);
    mMass = mass;
    rs = 2.0 * mGravConstant * mMass / (mSpeedOfLight * mSpeedOfLight);

    /*  Only a static tetrad is defined  */
    mLocTeds.push_back(enum_nat_tetrad_static);

    setStandardValues();
}

MetricSchwarzschildWT::~MetricSchwarzschildWT() {}

// *********************************** public methods ******************************

bool MetricSchwarzschildWT::calculateMetric(const double* pos)
{
    double r = pos[1];
    double theta = pos[2];

    double c = mSpeedOfLight;
    double c2 = c * c;

    double t1 = c2;
    double t3 = 1.0 / r;
    double t9 = r * r;
    double t10 = sin(theta);
    double t11 = t10 * t10;

    g_compts[0][0] = -t1;
    g_compts[0][1] = 0.0;
    g_compts[0][2] = 0.0;
    g_compts[0][3] = 0.0;
    g_compts[1][0] = 0.0;
    g_compts[1][1] = 1.0 / (1.0 - rs * t3);
    g_compts[1][2] = 0.0;
    g_compts[1][3] = 0.0;
    g_compts[2][0] = 0.0;
    g_compts[2][1] = 0.0;
    g_compts[2][2] = t9;
    g_compts[2][3] = 0.0;
    g_compts[3][0] = 0.0;
    g_compts[3][1] = 0.0;
    g_compts[3][2] = 0.0;
    g_compts[3][3] = t9 * t11;

    return true;
}

bool MetricSchwarzschildWT::calculateChristoffels(const double* pos)
{
    double r = pos[1];
    double theta = pos[2];

    // double c = mSpeedOfLight;

    double t1 = r - rs;
    //  double t2 = r * r;
    // double t6 = c * c;
    double t10 = 1.0 / r;
    //  double t14 = t10/t1*rs/2.0;
    double t14 = t10 / t1 * rs * 0.5;
    double t15 = sin(theta);
    double t17 = cos(theta);
    double t18 = 1.0 / t15 * t17;
    double t19 = t15 * t15;

    christoffel[0][0][0] = 0.0;
    christoffel[0][0][1] = 0.0;
    christoffel[0][0][2] = 0.0;
    christoffel[0][0][3] = 0.0;
    christoffel[0][1][0] = 0.0;
    christoffel[0][1][1] = 0.0;
    christoffel[0][1][2] = 0.0;
    christoffel[0][1][3] = 0.0;
    christoffel[0][2][0] = 0.0;
    christoffel[0][2][1] = 0.0;
    christoffel[0][2][2] = 0.0;
    christoffel[0][2][3] = 0.0;
    christoffel[0][3][0] = 0.0;
    christoffel[0][3][1] = 0.0;
    christoffel[0][3][2] = 0.0;
    christoffel[0][3][3] = 0.0;
    christoffel[1][0][0] = 0.0;
    christoffel[1][0][1] = 0.0;
    christoffel[1][0][2] = 0.0;
    christoffel[1][0][3] = 0.0;
    christoffel[1][1][0] = 0.0;
    christoffel[1][1][1] = -t14;
    christoffel[1][1][2] = 0.0;
    christoffel[1][1][3] = 0.0;
    christoffel[1][2][0] = 0.0;
    christoffel[1][2][1] = 0.0;
    christoffel[1][2][2] = t10;
    christoffel[1][2][3] = 0.0;
    christoffel[1][3][0] = 0.0;
    christoffel[1][3][1] = 0.0;
    christoffel[1][3][2] = 0.0;
    christoffel[1][3][3] = t10;
    christoffel[2][0][0] = 0.0;
    christoffel[2][0][1] = 0.0;
    christoffel[2][0][2] = 0.0;
    christoffel[2][0][3] = 0.0;
    christoffel[2][1][0] = 0.0;
    christoffel[2][1][1] = 0.0;
    christoffel[2][1][2] = t10;
    christoffel[2][1][3] = 0.0;
    christoffel[2][2][0] = 0.0;
    christoffel[2][2][1] = -t1;
    christoffel[2][2][2] = 0.0;
    christoffel[2][2][3] = 0.0;
    christoffel[2][3][0] = 0.0;
    christoffel[2][3][1] = 0.0;
    christoffel[2][3][2] = 0.0;
    christoffel[2][3][3] = t18;
    christoffel[3][0][0] = 0.0;
    christoffel[3][0][1] = 0.0;
    christoffel[3][0][2] = 0.0;
    christoffel[3][0][3] = 0.0;
    christoffel[3][1][0] = 0.0;
    christoffel[3][1][1] = 0.0;
    christoffel[3][1][2] = 0.0;
    christoffel[3][1][3] = t10;
    christoffel[3][2][0] = 0.0;
    christoffel[3][2][1] = 0.0;
    christoffel[3][2][2] = 0.0;
    christoffel[3][2][3] = t18;
    christoffel[3][3][0] = 0.0;
    christoffel[3][3][1] = -t1 * t19;
    christoffel[3][3][2] = -t15 * t17;
    christoffel[3][3][3] = 0.0;
    return true;
}

void MetricSchwarzschildWT::localToCoord(const double* pos, const double* ldir, double* dir, enum_nat_tetrad_type)
{
    double r = pos[1];
    double theta = pos[2];
    double w = sqrt(1.0 - rs / r);

    dir[0] = ldir[0] / mSpeedOfLight;
    dir[1] = ldir[1] * w;
    dir[2] = ldir[2] / r;
    dir[3] = ldir[3] / (r * sin(theta));
}

void MetricSchwarzschildWT::coordToLocal(const double* pos, const double* cdir, double* ldir, enum_nat_tetrad_type)
{
    double r = pos[1];
    double theta = pos[2];
    double w = sqrt(1.0 - rs / r);

    ldir[0] = cdir[0] * mSpeedOfLight;
    ldir[1] = cdir[1] / w;
    ldir[2] = cdir[2] * r;
    ldir[3] = cdir[3] * r * sin(theta);
}

bool MetricSchwarzschildWT::breakCondition(const double* pos)
{
    bool br = false;

    if ((pos[1] < 0.0) || (pos[1] * pos[1] <= (1.0 + M4D_METRIC_EPS) * rs * rs)) {
        br = true;
    }
    return br;
}

double MetricSchwarzschildWT::testConstraint(const double y[], const double kappa)
{
    double r = y[1];
    double theta = y[2];
    double cm = 1.0 / mSpeedOfLight;

    // Scale the directions with the speed of light before doubling them !!
    double dt = y[4];
    double dr = y[5] * cm;
    double dth = y[6] * cm;
    double dph = y[7] * cm;

    double sum = -kappa;
    sum += -dt * dt + dr * dr / (1.0 - rs / r) + r * r * (dth * dth + sin(theta) * sin(theta) * dph * dph);
    return sum;
}

bool MetricSchwarzschildWT::setParam(const char* pName, double val)
{
    if (Metric::setParam(pName, val)) {
        mMass = val;
        rs = 2.0 * mGravConstant * mMass / (mSpeedOfLight * mSpeedOfLight);
        return true;
    }

    return false;
}

bool MetricSchwarzschildWT::report(const vec4, const vec4, char*& text)
{
    std::stringstream ss;
    ss << "Report for Schwarzschild metric\n\tcoordinates : (t,r,theta,phi)\n";
    ss << "---------------------------------------------------------------\n";
    ss << "  physical units ................................. yes\n";
    ss.precision(DEF_FIXED_REPORT_PRECISION);
    ss.setf(std::ios::fixed);

    text = new char[ss.str().length() + 2];
    return CopyString(ss.str().c_str(), text);
}

// ********************************* protected methods *****************************

void MetricSchwarzschildWT::setStandardValues()
{
    mInitPos[0] = 0.0;
    mInitPos[1] = 3.0 * rs;
    mInitPos[2] = M_PI_2;
    mInitPos[3] = 0.0;
    mInitDir[0] = 1.0;
    mInitDir[1] = 0.0;
    mInitDir[2] = 0.0;

    mCoordNames[0] = std::string("t");
    mCoordNames[1] = std::string("r");
    mCoordNames[2] = std::string("theta");
    mCoordNames[3] = std::string("phi");
}

} // end namespace m4d
