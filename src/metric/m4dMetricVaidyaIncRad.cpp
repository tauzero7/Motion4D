/**
 * @file    m4dMetricVaidyaIncRad.cpp
 * @author  Thomas Mueller
 *
 * This file is part of the m4d-library.
 */
#include "m4dMetricVaidyaIncRad.h"

namespace m4d {

MetricVaidyaIncRad::MetricVaidyaIncRad(double k)
{
    mMetricName = "VaidyaIncRad";
    setCoordType(enum_coordinate_spherical);

    mPhysicalUnits = enum_physical_constants_geom;
    mSpeedOfLight = 1.0;
    mGravConstant = 1.0;

    addParam("k", k);
    addParam("vl", 0.0);
    addParam("sigma", 0.5);
    m_k = k;
    m_vl = 0.0;
    m_sigma = 0.5;

    mLocTeds.push_back(enum_nat_tetrad_static);

    setStandardValues();
}

MetricVaidyaIncRad::~MetricVaidyaIncRad() {}

// *********************************** public methods ******************************
bool MetricVaidyaIncRad::calculateMetric(const double* pos)
{
    double r = pos[1];
    double theta = pos[2];

    double m;
    calcMassFunc(pos[0], m);

    double t1 = m;
    double t6 = r * r;
    double t7 = sin(theta);
    double t8 = t7 * t7;

    g_compts[0][0] = -1.0 + 2.0 * t1 / r;
    g_compts[0][1] = 1.0;
    g_compts[0][2] = 0.0;
    g_compts[0][3] = 0.0;
    g_compts[1][0] = 1.0;
    g_compts[1][1] = 0.0;
    g_compts[1][2] = 0.0;
    g_compts[1][3] = 0.0;
    g_compts[2][0] = 0.0;
    g_compts[2][1] = 0.0;
    g_compts[2][2] = t6;
    g_compts[2][3] = 0.0;
    g_compts[3][0] = 0.0;
    g_compts[3][1] = 0.0;
    g_compts[3][2] = 0.0;
    g_compts[3][3] = t6 * t8;

    return true;
}

bool MetricVaidyaIncRad::calculateChristoffels(const double* pos)
{
    double r = pos[1];
    double theta = pos[2];

    double m, dmdv;
    calcMassFunc(pos[0], m, dmdv);

    double t1 = m;
    double t2 = r * r;
    double t4 = t1 / t2;
    double t5 = dmdv; // diff(m(v),v);
    double t8 = t1 * t1;
    double t14 = 1 / r;
    double t16 = -r + 2.0 * t1;
    double t17 = sin(theta);
    double t19 = cos(theta);
    double t20 = 1 / t17 * t19;
    double t21 = t17 * t17;

    christoffel[0][0][0] = t4;
    christoffel[0][0][1] = -(-t5 * t2 - t1 * r + 2.0 * t8) / t2 / r;
    christoffel[0][0][2] = 0.0;
    christoffel[0][0][3] = 0.0;
    christoffel[0][1][0] = 0.0;
    christoffel[0][1][1] = -t4;
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
    christoffel[1][0][1] = -t4;
    christoffel[1][0][2] = 0.0;
    christoffel[1][0][3] = 0.0;
    christoffel[1][1][0] = 0.0;
    christoffel[1][1][1] = 0.0;
    christoffel[1][1][2] = 0.0;
    christoffel[1][1][3] = 0.0;
    christoffel[1][2][0] = 0.0;
    christoffel[1][2][1] = 0.0;
    christoffel[1][2][2] = t14;
    christoffel[1][2][3] = 0.0;
    christoffel[1][3][0] = 0.0;
    christoffel[1][3][1] = 0.0;
    christoffel[1][3][2] = 0.0;
    christoffel[1][3][3] = t14;
    christoffel[2][0][0] = 0.0;
    christoffel[2][0][1] = 0.0;
    christoffel[2][0][2] = 0.0;
    christoffel[2][0][3] = 0.0;
    christoffel[2][1][0] = 0.0;
    christoffel[2][1][1] = 0.0;
    christoffel[2][1][2] = t14;
    christoffel[2][1][3] = 0.0;
    christoffel[2][2][0] = -r;
    christoffel[2][2][1] = t16;
    christoffel[2][2][2] = 0.0;
    christoffel[2][2][3] = 0.0;
    christoffel[2][3][0] = 0.0;
    christoffel[2][3][1] = 0.0;
    christoffel[2][3][2] = 0.0;
    christoffel[2][3][3] = t20;
    christoffel[3][0][0] = 0.0;
    christoffel[3][0][1] = 0.0;
    christoffel[3][0][2] = 0.0;
    christoffel[3][0][3] = 0.0;
    christoffel[3][1][0] = 0.0;
    christoffel[3][1][1] = 0.0;
    christoffel[3][1][2] = 0.0;
    christoffel[3][1][3] = t14;
    christoffel[3][2][0] = 0.0;
    christoffel[3][2][1] = 0.0;
    christoffel[3][2][2] = 0.0;
    christoffel[3][2][3] = t20;
    christoffel[3][3][0] = -r * t21;
    christoffel[3][3][1] = t16 * t21;
    christoffel[3][3][2] = -t17 * t19;
    christoffel[3][3][3] = 0.0;

    return true;
}

void MetricVaidyaIncRad::localToCoord(const double* pos, const double* ldir, double* dir, enum_nat_tetrad_type)
{
    double r = pos[1];
    double theta = pos[2];

    double m;
    calcMassFunc(pos[0], m);

    // static tetrad
    double w = sqrt(1.0 - 2.0 * m / r);
    dir[0] = (ldir[0] + ldir[1]) / w;
    dir[1] = ldir[1] * w;

    dir[2] = ldir[2] / r;
    dir[3] = ldir[3] / (r * sin(theta));
}

void MetricVaidyaIncRad::coordToLocal(const double* pos, const double* cdir, double* ldir, enum_nat_tetrad_type)
{
    // fprintf(stderr, "uups...  MetricVaidyaIncRad::coordToLocal() not implemented yet!\n");
    double r = pos[1];
    double theta = pos[2];

    double m;
    calcMassFunc(pos[0], m);

    // static tetrad
    double w = sqrt(1.0 - 2.0 * m / r);
    ldir[0] = w * cdir[0] - cdir[1] / w;
    ldir[1] = cdir[1] / w;

    ldir[2] = cdir[2] * r;
    ldir[3] = cdir[3] * r * sin(theta);
}

bool MetricVaidyaIncRad::breakCondition(const double* pos)
{
    double m;
    calcMassFunc(pos[0], m);
    double r = pos[1];
    if (r * r <= (1.0 + M4D_METRIC_EPS) * 4.0 * m * m) {
        return true;
    }
    return false;
}

bool MetricVaidyaIncRad::report(const vec4, const vec4, char*& text)
{
    std::stringstream ss;
    ss << "Report for Vaidya metric\n\tcoordinate : (v,r,theta,phi)\n";
    ss << "---------------------------------------------------------------\n";
    ss << "  physical units ..................... no\n";
    ss.precision(DEF_FIXED_REPORT_PRECISION);
    ss.setf(std::ios::fixed);

    text = new char[ss.str().length() + 2];
    return CopyString(ss.str().c_str(), text);
}

bool MetricVaidyaIncRad::setParam(const char* pName, double val)
{
    Metric::setParam(pName, val);
    if (strcmp(pName, "k") == 0) {
        m_k = val;
        return true;
    }
    else if (strcmp(pName, "vl") == 0) {
        m_vl = val;
        return true;
    }
    else if (strcmp(pName, "sigma") == 0) {
        m_sigma = val;
        return true;
    }
    return false;
}

// ********************************* protected methods *****************************

void MetricVaidyaIncRad::setStandardValues()
{
    mInitPos[0] = 0.0;
    mInitPos[1] = 10.0;
    mInitPos[2] = M_PI_2;
    mInitPos[3] = 0.0;
    mInitDir[0] = 1.0;
    mInitDir[1] = 0.0;
    mInitDir[2] = 0.0;

    mCoordNames[0] = std::string("v");
    mCoordNames[1] = std::string("r");
    mCoordNames[2] = std::string("theta");
    mCoordNames[3] = std::string("phi");
}

void MetricVaidyaIncRad::calcMassFunc(const double v, double& m)
{
    if (v >= m_vl) {
        m = 0.0;
        return;
    }

    double th = tanh(m_sigma * (m_vl - v));
    m = m_k * pow(m_vl - v, 1.0 / 3.0) * th * th;
}

void MetricVaidyaIncRad::calcMassFunc(const double v, double& m, double& dmdv)
{

    if (v >= m_vl) {
        m = dmdv = 0.0;
        return;
    }

    double th = tanh(m_sigma * (m_vl - v));
    double ch = cosh(m_sigma * (m_vl - v));

    m = m_k * pow(m_vl - v, 1.0 / 3.0) * th * th;

    dmdv = -m_k / 3.0 * pow(m_vl - v, -2.0 / 3.0) * th * th;
    dmdv -= m_k * pow(m_vl - v, 1.0 / 3.0) * 2.0 * m_sigma * th / (ch * ch);
}

} // end namespace m4d
