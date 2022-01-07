/**
 * @file    m4dMetricAlcubierreSimple.cpp
 * @author  Thomas Mueller
 *
 *  This file is part of the m4d-library.
 */
#include "m4dMetricAlcubierreSimple.h"

namespace m4d {

MetricAlcubierreSimple::MetricAlcubierreSimple(double R, double dR, double vs)
{
    mMetricName = "AlcubierreWarpSimple";
    mMetricCPPfilename = "m4dMetricAlcubierreSimple.cpp";
    setCoordType(enum_coordinate_cartesian);

    mPhysicalUnits = enum_physical_constants_geom;
    mSpeedOfLight = 1.0;
    mGravConstant = 1.0;

    mR = R;
    mDR = dR;
    mvs = vs;

    addParam("r", R);
    addParam("dr", dR);
    addParam("vs", vs);

    // mDrawTypes.push_back(enum_draw_twoplusone);

    setStandardValues();

    mLocTeds.push_back(enum_nat_tetrad_comoving);
    mLocTeds.push_back(enum_nat_tetrad_static);
}

MetricAlcubierreSimple::~MetricAlcubierreSimple() {}

// *********************************** public methods ******************************

bool MetricAlcubierreSimple::calculateMetric(const double* pos)
{
    double c = mSpeedOfLight;
    double vs = mvs;

    double f = calcF(pos);

    g_compts[0][0] = -pow(c, 2) + pow(vs, 2) * pow(f, 2);
    g_compts[0][1] = -vs * f;
    g_compts[0][2] = 0;
    g_compts[0][3] = 0;
    g_compts[1][0] = -vs * f;
    g_compts[1][1] = 1;
    g_compts[1][2] = 0;
    g_compts[1][3] = 0;
    g_compts[2][0] = 0;
    g_compts[2][1] = 0;
    g_compts[2][2] = 1;
    g_compts[2][3] = 0;
    g_compts[3][0] = 0;
    g_compts[3][1] = 0;
    g_compts[3][2] = 0;
    g_compts[3][3] = 1;

    return true;
}

bool MetricAlcubierreSimple::calculateChristoffels(const double* pos)
{
    double c = mSpeedOfLight;
    double vs = mvs;

    double f = calcF(pos);

    double ft, fx, fy, fz;
    calcDF(pos, ft, fx, fy, fz);

    christoffel[0][0][0] = pow(vs, 3) * pow(f, 2) * fx / pow(c, 2);
    christoffel[0][0][1] = vs * (-pow(c, 2) * vs * f * fx - pow(c, 2) * ft + pow(vs, 3) * pow(f, 3) * fx) / pow(c, 2);
    christoffel[0][0][2] = -pow(vs, 2) * f * fy;
    christoffel[0][0][3] = -pow(vs, 2) * f * fz;
    christoffel[0][1][0] = -pow(vs, 2) * f * fx / pow(c, 2);
    christoffel[0][1][1] = -pow(vs, 3) * pow(f, 2) * fx / pow(c, 2);
    christoffel[0][1][2] = vs * fy / 2;
    christoffel[0][1][3] = vs * fz / 2;
    christoffel[0][2][0] = -pow(vs, 2) * f * fy / (2 * pow(c, 2));
    christoffel[0][2][1] = vs * (-pow(c, 2) - pow(vs, 2) * pow(f, 2)) * fy / (2 * pow(c, 2));
    christoffel[0][2][2] = 0;
    christoffel[0][2][3] = 0;
    christoffel[0][3][0] = -pow(vs, 2) * f * fz / (2 * pow(c, 2));
    christoffel[0][3][1] = vs * (-pow(c, 2) - pow(vs, 2) * pow(f, 2)) * fz / (2 * pow(c, 2));
    christoffel[0][3][2] = 0;
    christoffel[0][3][3] = 0;
    christoffel[1][0][0] = -pow(vs, 2) * f * fx / pow(c, 2);
    christoffel[1][0][1] = -pow(vs, 3) * pow(f, 2) * fx / pow(c, 2);
    christoffel[1][0][2] = vs * fy / 2;
    christoffel[1][0][3] = vs * fz / 2;
    christoffel[1][1][0] = vs * fx / pow(c, 2);
    christoffel[1][1][1] = pow(vs, 2) * f * fx / pow(c, 2);
    christoffel[1][1][2] = 0;
    christoffel[1][1][3] = 0;
    christoffel[1][2][0] = vs * fy / (2 * pow(c, 2));
    christoffel[1][2][1] = pow(vs, 2) * f * fy / (2 * pow(c, 2));
    christoffel[1][2][2] = 0;
    christoffel[1][2][3] = 0;
    christoffel[1][3][0] = vs * fz / (2 * pow(c, 2));
    christoffel[1][3][1] = pow(vs, 2) * f * fz / (2 * pow(c, 2));
    christoffel[1][3][2] = 0;
    christoffel[1][3][3] = 0;
    christoffel[2][0][0] = -pow(vs, 2) * f * fy / (2 * pow(c, 2));
    christoffel[2][0][1] = vs * (-pow(c, 2) - pow(vs, 2) * pow(f, 2)) * fy / (2 * pow(c, 2));
    christoffel[2][0][2] = 0;
    christoffel[2][0][3] = 0;
    christoffel[2][1][0] = vs * fy / (2 * pow(c, 2));
    christoffel[2][1][1] = pow(vs, 2) * f * fy / (2 * pow(c, 2));
    christoffel[2][1][2] = 0;
    christoffel[2][1][3] = 0;
    christoffel[2][2][0] = 0;
    christoffel[2][2][1] = 0;
    christoffel[2][2][2] = 0;
    christoffel[2][2][3] = 0;
    christoffel[2][3][0] = 0;
    christoffel[2][3][1] = 0;
    christoffel[2][3][2] = 0;
    christoffel[2][3][3] = 0;
    christoffel[3][0][0] = -pow(vs, 2) * f * fz / (2 * pow(c, 2));
    christoffel[3][0][1] = vs * (-pow(c, 2) - pow(vs, 2) * pow(f, 2)) * fz / (2 * pow(c, 2));
    christoffel[3][0][2] = 0;
    christoffel[3][0][3] = 0;
    christoffel[3][1][0] = vs * fz / (2 * pow(c, 2));
    christoffel[3][1][1] = pow(vs, 2) * f * fz / (2 * pow(c, 2));
    christoffel[3][1][2] = 0;
    christoffel[3][1][3] = 0;
    christoffel[3][2][0] = 0;
    christoffel[3][2][1] = 0;
    christoffel[3][2][2] = 0;
    christoffel[3][2][3] = 0;
    christoffel[3][3][0] = 0;
    christoffel[3][3][1] = 0;
    christoffel[3][3][2] = 0;
    christoffel[3][3][3] = 0;

    return true;
}

bool MetricAlcubierreSimple::calculateChrisD(const double* pos)
{
    double ft, fx, fy, fz;
    calcDF(pos, ft, fx, fy, fz);

    double ftt, ftx, fty, ftz, fxx, fxy, fxz, fyy, fyz, fzz;
    calcD2F(pos, ftt, ftx, fty, ftz, fxx, fxy, fxz, fyy, fyz, fzz);

    return false;
}

bool MetricAlcubierreSimple::calculateRiemann(const double* pos)
{
    double ft, fx, fy, fz;
    calcDF(pos, ft, fx, fy, fz);

    double ftt, ftx, fty, ftz, fxx, fxy, fxz, fyy, fyz, fzz;
    calcD2F(pos, ftt, ftx, fty, ftz, fxx, fxy, fxz, fyy, fyz, fzz);

    return false;
}

void MetricAlcubierreSimple::localToCoord(const double* pos, const double* ldir, double* dir, enum_nat_tetrad_type type)
{
    double f = calcF(pos);
    double c = mSpeedOfLight;

    if (type == enum_nat_tetrad_comoving) {
        dir[0] = ldir[0] / c;
        dir[1] = ldir[0] * mvs * f / c + ldir[1];
        dir[2] = ldir[2];
        dir[3] = ldir[3];
    }
    else {
        double w = sqrt(c * c - mvs * mvs * f * f);

        dir[0] = (ldir[0] - mvs * f / c * ldir[1]) / w;
        dir[1] = w / c * ldir[1];
        dir[2] = ldir[2];
        dir[3] = ldir[3];
    }
}

void MetricAlcubierreSimple::coordToLocal(
    const double* pos, const double* cdir, double* ldir, enum_nat_tetrad_type type)
{
    double f = calcF(pos);
    double c = mSpeedOfLight;

    if (type == enum_nat_tetrad_comoving) {
        ldir[0] = c * cdir[0];
        ldir[1] = cdir[1] - mvs * f * cdir[0];
        ldir[2] = cdir[2];
        ldir[3] = cdir[3];
    }
    else {
        double w = sqrt(c * c - mvs * mvs * f * f);

        ldir[1] = c / w * cdir[1];
        ldir[0] = w * cdir[0] - ldir[1] * mvs * f / c;
        ldir[2] = cdir[2];
        ldir[3] = cdir[3];
    }
}

bool MetricAlcubierreSimple::breakCondition(const double*)
{
    bool br = false;

    return br;
}

double MetricAlcubierreSimple::testConstraint(const double y[], const double kappa)
{
    double c = mSpeedOfLight;
    double f = calcF(y);

    double sum = -kappa;
    sum += -c * c * y[4] * y[4] + pow(y[5] - mvs * f * y[4], 2.0) + y[6] * y[6] + y[7] * y[7];

    // double k = -c * c * y[4] + mvs * (1 - f) * (y[5] - mvs * f * y[4]);
    // std::cerr << k << std::endl;

    double dxdt = y[5] / y[4];
    double dydt = y[6] / y[4];
    double w = (dxdt - mvs * f) * (dxdt - mvs * f) + dydt * dydt;
    std::cerr << y[4] << " " << y[5] << " " << y[6] << " " << w << std::endl;

    // double A = 1.0-mvs*mvs*(1.0-f)*(1.0-f);
    // fprintf(stderr,"%e %e %e %e  %e %e %e %e\n",y[4],y[5],y[6],y[7],f,A,mvs*f+1,mvs*f-1);
    return sum;
}

bool MetricAlcubierreSimple::resize(double* y, double, double)
{
    double a = 1 / fabs(y[4]);
    y[4] *= a;
    y[5] *= a;
    y[6] *= a;
    y[7] *= a;
    return true;
}

bool MetricAlcubierreSimple::setParam(const char* pName, double val)
{
    Metric::setParam(pName, val);
    if (strcmp(pName, "r") == 0) {
        mR = val;
    }
    else if (strcmp(pName, "dr") == 0) {
        mDR = val;
    }
    else if (strcmp(pName, "vs") == 0) {
        mvs = val;
    }
    return true;
}

bool MetricAlcubierreSimple::transToTwoPlusOne(vec4 p, vec4& cp)
{
    cp = vec4(p[0], p[1], p[2], p[0]);
    return true;
}

/*! Generate report.
 */
bool MetricAlcubierreSimple::report(const vec4, const vec4, char*& text)
{
    std::stringstream ss;
    ss << "Report for AlcubierreSimple metric\n\tcoordinate : (t,x,y,z)\n";
    ss << "---------------------------------------------------------------\n";
    ss << "  physical units ................................. no\n";
    ss.precision(DEF_FIXED_REPORT_PRECISION);
    ss.setf(std::ios::fixed);

    text = new char[ss.str().length() + 2];
    return CopyString(ss.str().c_str(), text);
}

// ********************************* protected methods *****************************
void MetricAlcubierreSimple::setStandardValues()
{
    mInitPos[0] = 0.0;
    mInitPos[1] = 0.0;
    mInitPos[2] = -10.0;
    mInitPos[3] = 0.0;
    mInitDir[0] = 1.0;
    mInitDir[1] = 0.0;
    mInitDir[2] = 0.0;

    mCoordNames[0] = std::string("t");
    mCoordNames[1] = std::string("x");
    mCoordNames[2] = std::string("y");
    mCoordNames[3] = std::string("z");
}

double MetricAlcubierreSimple::calcRs(const double* pos)
{
    double t = pos[0];
    double x = pos[1] - mvs * t;
    double y = pos[2];
    double z = pos[3];

    return sqrt(x * x + y * y + z * z);
}

double MetricAlcubierreSimple::calcF(const double* pos)
{
    double rs = calcRs(pos);

#if 0
    if (rs <= mR) {
        return 1.0 - pow(rs / mR, 4.0);
    }
    return 0.0;
#else
    double R1 = mR - 0.5 * mDR;
    double R2 = mR + 0.5 * mDR;
    if (rs <= R1) {
        return 1.0;
    }

    if (rs >= R1 && rs < R2) {
        double w = (rs - R1) / (R2 - R1);
        double y = 1 - w * w;
        return y * y;
    }

    return 0.0;
#endif
}

void MetricAlcubierreSimple::calcDF(const double* pos, double& ft, double& fx, double& fy, double& fz)
{
    double rs = calcRs(pos);

    double t = pos[0];
    double x = pos[1] - mvs * t;
    double y = pos[2];
    double z = pos[3];

#if 0
    double dfdr = 0.0;
    if (rs <= mR) {
        dfdr = -4.0 * pow(rs / mR, 3.0) / mR;
    }
    double df = dfdr / rs;

    // ft = df/dr * dr/dt ...
    ft = -mvs * x * df;
    fx = x * df;
    fy = y * df;
    fz = z * df;
#else
    double R1 = mR - 0.5 * mDR;
    double R2 = mR + 0.5 * mDR;

    double dfdr = 0.0;

    if (rs >= R1 && rs < R2) {
        double w = (rs - R1) / (R2 - R1);
        dfdr = -4.0 * (1 - w * w) * w / (R2 - R1);
    }

    ft = dfdr * x / rs * (-mvs);
    fx = dfdr * x / rs;
    fy = dfdr * y / rs;
    fz = dfdr * z / rs;
#endif
}

void MetricAlcubierreSimple::calcD2F(
    const double*, double&, double&, double&, double&, double&, double&, double&, double&, double&, double&)
{
    fprintf(stderr, "uups... not implemented yet!\n");
    // TODO
}

} // end namespace m4d
