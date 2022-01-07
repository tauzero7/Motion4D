/**
 * @file    m4dMetricMinkRotLattice.cpp
 * @author  Thomas Mueller
 *
 * This file is part of the m4d-library.
 */
#include "m4dMetricMinkRotLattice.h"

namespace m4d {

MetricMinkRotLattice::MetricMinkRotLattice(double omega)
{
    mMetricName = "RotLattice";
    setCoordType(enum_coordinate_cylinder);

    mOmega = omega;

    mPhysicalUnits = enum_physical_constants_geom;
    mSpeedOfLight = 1.0;
    mGravConstant = 1.0;

    mLocTeds.push_back(enum_nat_tetrad_static);
    mLocTeds.push_back(enum_nat_tetrad_lnrf);

    addParam("omega", mOmega);

    setStandardValues();
}

MetricMinkRotLattice::~MetricMinkRotLattice() {}

// *********************************** public methods ******************************

bool MetricMinkRotLattice::calculateMetric(const double* pos)
{
    double c = mSpeedOfLight;
    double r = pos[1];

    double t1 = c * c;
    double t2 = r * r;
    double t3 = mOmega * mOmega;
    double t6 = mOmega * t2;

    g_compts[0][0] = -t1 + t2 * t3;
    g_compts[0][1] = 0.0;
    g_compts[0][2] = t6;
    g_compts[0][3] = 0.0;
    g_compts[1][0] = 0.0;
    g_compts[1][1] = 1.0;
    g_compts[1][2] = 0.0;
    g_compts[1][3] = 0.0;
    g_compts[2][0] = t6;
    g_compts[2][1] = 0.0;
    g_compts[2][2] = t2;
    g_compts[2][3] = 0.0;
    g_compts[3][0] = 0.0;
    g_compts[3][1] = 0.0;
    g_compts[3][2] = 0.0;
    g_compts[3][3] = 1.0;

    return true;
}

bool MetricMinkRotLattice::calculateChristoffels(const double* pos)
{
    double r = pos[1];

    double t1 = mOmega * mOmega;
    double t3 = 1 / r;
    double t4 = mOmega * t3;
    double t5 = mOmega * r;

    christoffel[0][0][0] = 0.0;
    christoffel[0][0][1] = -r * t1;
    christoffel[0][0][2] = 0.0;
    christoffel[0][0][3] = 0.0;
    christoffel[0][1][0] = 0.0;
    christoffel[0][1][1] = 0.0;
    christoffel[0][1][2] = t4;
    christoffel[0][1][3] = 0.0;
    christoffel[0][2][0] = 0.0;
    christoffel[0][2][1] = -t5;
    christoffel[0][2][2] = 0.0;
    christoffel[0][2][3] = 0.0;
    christoffel[0][3][0] = 0.0;
    christoffel[0][3][1] = 0.0;
    christoffel[0][3][2] = 0.0;
    christoffel[0][3][3] = 0.0;
    christoffel[1][0][0] = 0.0;
    christoffel[1][0][1] = 0.0;
    christoffel[1][0][2] = t4;
    christoffel[1][0][3] = 0.0;
    christoffel[1][1][0] = 0.0;
    christoffel[1][1][1] = 0.0;
    christoffel[1][1][2] = 0.0;
    christoffel[1][1][3] = 0.0;
    christoffel[1][2][0] = 0.0;
    christoffel[1][2][1] = 0.0;
    christoffel[1][2][2] = t3;
    christoffel[1][2][3] = 0.0;
    christoffel[1][3][0] = 0.0;
    christoffel[1][3][1] = 0.0;
    christoffel[1][3][2] = 0.0;
    christoffel[1][3][3] = 0.0;
    christoffel[2][0][0] = 0.0;
    christoffel[2][0][1] = -t5;
    christoffel[2][0][2] = 0.0;
    christoffel[2][0][3] = 0.0;
    christoffel[2][1][0] = 0.0;
    christoffel[2][1][1] = 0.0;
    christoffel[2][1][2] = t3;
    christoffel[2][1][3] = 0.0;
    christoffel[2][2][0] = 0.0;
    christoffel[2][2][1] = -r;
    christoffel[2][2][2] = 0.0;
    christoffel[2][2][3] = 0.0;
    christoffel[2][3][0] = 0.0;
    christoffel[2][3][1] = 0.0;
    christoffel[2][3][2] = 0.0;
    christoffel[2][3][3] = 0.0;
    christoffel[3][0][0] = 0.0;
    christoffel[3][0][1] = 0.0;
    christoffel[3][0][2] = 0.0;
    christoffel[3][0][3] = 0.0;
    christoffel[3][1][0] = 0.0;
    christoffel[3][1][1] = 0.0;
    christoffel[3][1][2] = 0.0;
    christoffel[3][1][3] = 0.0;
    christoffel[3][2][0] = 0.0;
    christoffel[3][2][1] = 0.0;
    christoffel[3][2][2] = 0.0;
    christoffel[3][2][3] = 0.0;
    christoffel[3][3][0] = 0.0;
    christoffel[3][3][1] = 0.0;
    christoffel[3][3][2] = 0.0;
    christoffel[3][3][3] = 0.0;

    return true;
}

void MetricMinkRotLattice::localToCoord(const double* pos, const double* ldir, double* dir, enum_nat_tetrad_type type)
{
    double r = pos[1];
    double wdc = mOmega / mSpeedOfLight;

    if (type == enum_nat_tetrad_default || type == enum_nat_tetrad_static) {
        double sq = sqrt(1.0 - r * r * wdc * wdc);

        dir[0] = ldir[0] / sq / mSpeedOfLight + ldir[2] * wdc * r / sq / mSpeedOfLight;
        dir[1] = ldir[1];
        dir[2] = ldir[2] * sq / r;
        dir[3] = ldir[3];
    }
    else if (type == enum_nat_tetrad_lnrf) {
        dir[0] = ldir[0] / mSpeedOfLight;
        dir[1] = ldir[1];
        dir[2] = ldir[2] / r - wdc * ldir[0];
        dir[3] = ldir[3];
    }
}

void MetricMinkRotLattice::coordToLocal(const double* pos, const double* cdir, double* ldir, enum_nat_tetrad_type type)
{
    double r = pos[1];
    double wdc = mOmega / mSpeedOfLight;

    if (type == enum_nat_tetrad_default || type == enum_nat_tetrad_static) {
        double sq = sqrt(1.0 - r * r * wdc * wdc);

        ldir[0] = cdir[0] * mSpeedOfLight * sq - cdir[2] * wdc * r * r / sq;
        ldir[1] = cdir[1];
        ldir[2] = r / sq * cdir[2];
        ldir[3] = cdir[3];
    }
    else if (type == enum_nat_tetrad_lnrf) {
        ldir[0] = cdir[0] * mSpeedOfLight;
        ldir[1] = cdir[1];
        ldir[2] = r * (cdir[2] + mOmega * cdir[0]);
        ldir[3] = cdir[3];
    }
}

bool MetricMinkRotLattice::breakCondition(const double* pos)
{
    if (pos[1] >= fabs(mSpeedOfLight / mOmega)) {
        return true;
    }

    return false;
}

bool MetricMinkRotLattice::setParam(const char* pName, double val)
{
    if (Metric::setParam(pName, val)) {
        mOmega = val;
    }
    return true;
}

double MetricMinkRotLattice::testConstraint(const double y[], const double kappa)
{
    double r = y[1];
    double cm = 1.0 / mSpeedOfLight;

    // Scale the directions with the speed of light before doubling them !!
    double dt = y[4];
    double dr = y[5] * cm;
    double dph = y[6] * cm;
    double dz = y[7] * cm;

    double wdc = mOmega * cm;

    double sum = -kappa;
    sum += -(1.0 - wdc * wdc * r * r) * dt * dt + 2.0 * r * r * cm * mOmega * dt * dph + r * r * dph * dph + dr * dr
        + dz * dz;

    return sum;
}

bool MetricMinkRotLattice::report(const vec4, const vec4, char*& text)
{
    std::stringstream ss;
    ss << "Report for Minkowski as rotating lattice metric\n\tcoordinate : (t,r,phi,z)\n";
    ss << "--------------------------------------------------------\n";
    ss << "  physical units ................................. no\n";
    ss << "  points outside  r=c/omega are not allowed ...... " << mSpeedOfLight / mOmega << std::endl;
    text = new char[ss.str().length() + 2];
    return CopyString(ss.str().c_str(), text);
}

// ********************************* protected methods *****************************

void MetricMinkRotLattice::setStandardValues()
{

    mInitPos[0] = 0.0;
    mInitPos[1] = 1.0;
    mInitPos[2] = 0.0;
    mInitPos[3] = 0.0;
    mInitDir[0] = 1.0;
    mInitDir[1] = 0.0;
    mInitDir[2] = 0.0;

    mCoordNames[0] = std::string("t");
    mCoordNames[1] = std::string("r");
    mCoordNames[2] = std::string("phi");
    mCoordNames[3] = std::string("z");
}

} // end namespace m4d
