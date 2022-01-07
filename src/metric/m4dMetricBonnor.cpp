/**
 * @file    m4dMetricBonnor.cpp
 * @author  Thomas Mueller
 *
 * This file is part of the m4d-library.
 */
#include "m4dMetricBonnor.h"

namespace m4d {

MetricBonnor::MetricBonnor(double mass, double b)
{
    mMetricName = "Bonnor";
    setCoordType(enum_coordinate_spherical);

    mPhysicalUnits = enum_physical_constants_geom;
    mSpeedOfLight = 1.0;
    mGravConstant = 1.0;
    mDielectricPerm = 1.0;

    addParam("mass", mass);
    mMass = mass;
    addParam("b", b);
    mB = b;

    mSign = 1.0;
    mLocTeds.push_back(enum_nat_tetrad_static);

    setStandardValues();
}

MetricBonnor::~MetricBonnor() {}

// *********************************** public methods ******************************

bool MetricBonnor::calculateMetric(const double* pos)
{
    double theta = pos[2];
    calcPotentials(pos);

    double t1 = P; // P(r,theta);
    double t2 = t1 * t1;
    double t3 = Y; // Y(r,theta);
    double t4 = t3 * t3;
    double t7 = t4 * t2;
    double t8 = Q; // Q(r,theta);
    double t9 = t8 * t8;
    double t11 = 1 / t9 / t8;
    double t12 = Z; // Z(r);
    double t19 = sin(theta);
    double t20 = t19 * t19;

    g_compts[0][0] = -t2 / t4;
    g_compts[0][1] = 0.0;
    g_compts[0][2] = 0.0;
    g_compts[0][3] = 0.0;
    g_compts[1][0] = 0.0;
    g_compts[1][1] = t7 * t11 / t12;
    g_compts[1][2] = 0.0;
    g_compts[1][3] = 0.0;
    g_compts[2][0] = 0.0;
    g_compts[2][1] = 0.0;
    g_compts[2][2] = t7 * t11;
    g_compts[2][3] = 0.0;
    g_compts[3][0] = 0.0;
    g_compts[3][1] = 0.0;
    g_compts[3][2] = 0.0;
    g_compts[3][3] = t12 * t4 / t2 * t20;
    return true;
}

bool MetricBonnor::calculateChristoffels(const double* pos)
{
    double theta = pos[2];
    calcPotiDiffs(pos);

    double t1 = P; // P(r,theta);
    double t2 = 1 / t1;
    double t3 = Y; // Y(r,theta);
    double t4 = t3 * t3;
    double t5 = t4 * t4;
    double t8 = t2 / t5 / t3;
    double t9 = Q; // Q(r,theta);
    double t10 = t9 * t9;
    double t11 = t10 * t9;
    double t12 = Z; // Z(r);
    double t14 = dPdr; // diff(P(r,theta),r);
    double t15 = t14 * t3;
    double t16 = dYdr; // diff(Y(r,theta),r);
    double t17 = t1 * t16;
    double t18 = -t15 + t17;
    double t21 = dPdtheta; // diff(P(r,theta),theta);
    double t22 = t21 * t3;
    double t23 = dYdtheta; // diff(Y(r,theta),theta);
    double t24 = t1 * t23;
    double t25 = -t22 + t24;
    double t28 = 1 / t3;
    double t29 = t28 * t2;
    double t30 = t29 * t18;
    double t31 = t29 * t25;
    double t32 = 1 / t9;
    double t33 = 1 / t12;
    double t35 = t9 * t12;
    double t40 = t3 * t1;
    double t41 = dQdr; // diff(Q(r,theta),r);
    double t45 = dZdr; // diff(Z(r),r);
    double t56 = dQdtheta; // diff(Q(r,theta),theta);
    double t60 = t32 * (-2.0 * t24 * t9 - 2.0 * t22 * t9 + 3.0 * t40 * t56);
    double t65 = t29 * t60 / 2.0;
    double t72 = -2.0 * t17 * t9 - 2.0 * t15 * t9 + 3.0 * t40 * t41;
    double t75 = t29 * t32 * t72 / 2.0;
    double t84 = t45 * t3 * t1 + 2.0 * t12 * t16 * t1 - 2.0 * t12 * t3 * t14;
    double t87 = t29 * t33 * t84 / 2.0;
    double t92 = sin(theta);
    double t98 = cos(theta);
    double t101 = t92 * t23 * t1 - t3 * t92 * t21 + t3 * t98 * t1;
    double t103 = t29 / t92 * t101;
    double t104 = t1 * t1;
    double t105 = t104 * t104;
    double t109 = 1 / t105 / t1 * t28 * t11;
    double t110 = t92 * t92;

    christoffel[0][0][0] = 0.0;
    christoffel[0][0][1] = -t8 * t11 * t12 * t18;
    christoffel[0][0][2] = -t8 * t11 * t25;
    christoffel[0][0][3] = 0.0;
    christoffel[0][1][0] = -t30;
    christoffel[0][1][1] = 0.0;
    christoffel[0][1][2] = 0.0;
    christoffel[0][1][3] = 0.0;
    christoffel[0][2][0] = -t31;
    christoffel[0][2][1] = 0.0;
    christoffel[0][2][2] = 0.0;
    christoffel[0][2][3] = 0.0;
    christoffel[0][3][0] = 0.0;
    christoffel[0][3][1] = 0.0;
    christoffel[0][3][2] = 0.0;
    christoffel[0][3][3] = 0.0;
    christoffel[1][0][0] = -t30;
    christoffel[1][0][1] = 0.0;
    christoffel[1][0][2] = 0.0;
    christoffel[1][0][3] = 0.0;
    christoffel[1][1][0] = 0.0;
    christoffel[1][1][1]
        = -t29 * t32 * t33 * (-2.0 * t17 * t35 - 2.0 * t15 * t35 + 3.0 * t40 * t41 * t12 + t40 * t45 * t9) / 2.0;
    christoffel[1][1][2] = t29 * t60 * t33 / 2.0;
    christoffel[1][1][3] = 0.0;
    christoffel[1][2][0] = 0.0;
    christoffel[1][2][1] = -t65;
    christoffel[1][2][2] = -t75;
    christoffel[1][2][3] = 0.0;
    christoffel[1][3][0] = 0.0;
    christoffel[1][3][1] = 0.0;
    christoffel[1][3][2] = 0.0;
    christoffel[1][3][3] = t87;
    christoffel[2][0][0] = -t31;
    christoffel[2][0][1] = 0.0;
    christoffel[2][0][2] = 0.0;
    christoffel[2][0][3] = 0.0;
    christoffel[2][1][0] = 0.0;
    christoffel[2][1][1] = -t65;
    christoffel[2][1][2] = -t75;
    christoffel[2][1][3] = 0.0;
    christoffel[2][2][0] = 0.0;
    christoffel[2][2][1] = t29 * t32 * t12 * t72 / 2.0;
    christoffel[2][2][2] = -t65;
    christoffel[2][2][3] = 0.0;
    christoffel[2][3][0] = 0.0;
    christoffel[2][3][1] = 0.0;
    christoffel[2][3][2] = 0.0;
    christoffel[2][3][3] = t103;
    christoffel[3][0][0] = 0.0;
    christoffel[3][0][1] = 0.0;
    christoffel[3][0][2] = 0.0;
    christoffel[3][0][3] = 0.0;
    christoffel[3][1][0] = 0.0;
    christoffel[3][1][1] = 0.0;
    christoffel[3][1][2] = 0.0;
    christoffel[3][1][3] = t87;
    christoffel[3][2][0] = 0.0;
    christoffel[3][2][1] = 0.0;
    christoffel[3][2][2] = 0.0;
    christoffel[3][2][3] = t103;
    christoffel[3][3][0] = 0.0;
    christoffel[3][3][1] = -t109 * t12 * t110 * t84 / 2.0;
    christoffel[3][3][2] = -t109 * t12 * t92 * t101;
    christoffel[3][3][3] = 0.0;

    return true;
}

bool MetricBonnor::calculateChrisD(const double*)
{
    return false;
}

void MetricBonnor::localToCoord(const double* pos, const double* ldir, double* dir, enum_nat_tetrad_type)
{
    calcPotentials(pos);
    double theta = pos[2];

    dir[0] = ldir[0] * Y / P;
    dir[1] = ldir[1] * pow(Q, 1.5) * sqrt(Z) / (Y * P);
    dir[2] = ldir[2] * pow(Q, 1.5) / (Y * P);
    dir[3] = ldir[3] * P / (Y * sqrt(Z) * sin(theta));
}

void MetricBonnor::coordToLocal(const double*, const double*, double*, enum_nat_tetrad_type) {}

bool MetricBonnor::breakCondition(const double* pos)
{
    bool br = false;
    double r = pos[1];
    if (r <= mMass + sqrt(mMass * mMass + mB * mB)) {
        br = true;
    }
    return br;
}

bool MetricBonnor::setParam(const char* pName, double val)
{
    Metric::setParam(pName, val);

    if (strcmp(pName, "mass") == 0) {
        mMass = val;
    }
    else if (strcmp(pName, "b") == 0) {
        mB = val;
    }
    return true;
}

bool MetricBonnor::report(const vec4, const vec4, char*& text)
{
    std::stringstream ss;
    ss << "Report for the  Bonnor metric\n\tcoordinate : (t,r,theta,phi)\n";
    ss << "---------------------------------------------------------------\n";
    ss.precision(DEF_FIXED_REPORT_PRECISION);
    ss.setf(std::ios::fixed);

    text = new char[ss.str().length() + 2];
    return CopyString(ss.str().c_str(), text);
}

// ********************************* protected methods *****************************
void MetricBonnor::setStandardValues()
{
    mInitPos[0] = 0.0;
    mInitPos[1] = 10.0;
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

void MetricBonnor::calcPotentials(const double* pos)
{
    double r = pos[1];
    double theta = pos[2];
    double ct = cos(theta);

    P = r * r - 2 * mMass * r - mB * mB * ct * ct;
    Q = (r - mMass) * (r - mMass) - (mB * mB + mMass * mMass) * ct * ct;
    Y = r * r - mB * mB * ct * ct;
    Z = r * r - 2 * mMass * r - mB * mB;
}

void MetricBonnor::calcPotiDiffs(const double* pos)
{
    double r = pos[1];
    double theta = pos[2];
    double ct = cos(theta);
    double st = sin(theta);

    calcPotentials(pos);

    dPdr = 2.0 * (r - mMass);
    dPdtheta = 2.0 * mB * mB * ct * st;

    dQdr = 2.0 * (r - mMass);
    dQdtheta = 2.0 * (mB * mB + mMass * mMass) * ct * st;

    dYdr = 2.0 * r;
    dYdtheta = 2.0 * mB * mB * ct * st;

    dZdr = 2.0 * (r - mMass);
}

} // end namespace m4d
