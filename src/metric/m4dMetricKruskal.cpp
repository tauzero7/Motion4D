/**
 * @file    m4dMetricKruskal.cpp
 * @author  Thomas Mueller
 *
 *  This file is part of libMotion4D.
 */
#include "m4dMetricKruskal.h"

#include <gsl/gsl_sf.h>
#include <gsl/gsl_sf_lambert.h>

namespace m4d {

#define eps 1.0e-6

MetricKruskal::MetricKruskal(double mass)
{
    mMetricName = "Kruskal";
    setCoordType(enum_coordinate_spherical);

    mPhysicalUnits = enum_physical_constants_geom;
    mSpeedOfLight = 1.0;
    mGravConstant = 1.0;
    mDielectricPerm = 1.0;

    addParam("mass", mass);
    mMass = mass;
    rs = 2.0 * mMass;

    mSign = 1.0;
    mLocTeds.push_back(enum_nat_tetrad_static);

    setStandardValues();
}

MetricKruskal::~MetricKruskal() {}

bool MetricKruskal::calculateMetric(const double* pos)
{
    double r = get_r(pos[0], pos[1]);
    double fac = 4.0 * rs * rs * rs / r * exp(-r / rs);
    double theta = pos[2];

    g_compts[0][0] = -fac;
    g_compts[0][1] = 0.0;
    g_compts[0][2] = 0.0;
    g_compts[0][3] = 0.0;
    g_compts[1][0] = 0.0;
    g_compts[1][1] = fac;
    g_compts[1][2] = 0.0;
    g_compts[1][3] = 0.0;
    g_compts[2][0] = 0.0;
    g_compts[2][1] = 0.0;
    g_compts[2][2] = r * r;
    g_compts[2][3] = 0.0;
    g_compts[3][0] = 0.0;
    g_compts[3][1] = 0.0;
    g_compts[3][2] = 0.0;
    g_compts[3][3] = r * r * sin(theta) * sin(theta);
    return true;
}

bool MetricKruskal::calculateChristoffels(const double* pos)
{
    double T = pos[0];
    double X = pos[1];
    double r = get_r(T, X);
    double theta = pos[2];

    double GTT_T = T * rs * (r + rs) / (r * r) * exp(-r / rs);
    double GTX_X = GTT_T;
    double GXX_T = GTT_T;

    double GTT_X = -X * rs * (r + rs) / (r * r) * exp(-r / rs);
    double GTX_T = GTT_X;
    double GXX_X = GTT_X;

    double GTt_t = -2.0 * rs * rs / (r * r) * T * exp(-r / rs);
    double GTp_p = GTt_t;

    double GXt_t = 2.0 * rs * rs / (r * r) * X * exp(-r / rs);
    double GXp_p = GXt_t;

    double Gtt_T = -r / (2 * rs) * T;
    double Gtt_X = -r / (2 * rs) * X;
    double Gpp_T = -r / (2 * rs) * T * sin(theta) * sin(theta);
    double Gpp_X = -r / (2 * rs) * X * sin(theta) * sin(theta);
    double Gtp_p = cos(theta) / sin(theta);
    double Gpp_t = -sin(theta) * cos(theta);

    christoffel[0][0][0] = GTT_T;
    christoffel[0][0][1] = GTT_X;
    christoffel[0][0][2] = 0.0;
    christoffel[0][0][3] = 0.0;
    christoffel[0][1][0] = GTX_T;
    christoffel[0][1][1] = GTX_X;
    christoffel[0][1][2] = 0.0;
    christoffel[0][1][3] = 0.0;
    christoffel[0][2][0] = 0.0;
    christoffel[0][2][1] = 0.0;
    christoffel[0][2][2] = GTt_t;
    christoffel[0][2][3] = 0.0;
    christoffel[0][3][0] = 0.0;
    christoffel[0][3][1] = 0.0;
    christoffel[0][3][2] = 0.0;
    christoffel[0][3][3] = GTp_p;
    christoffel[1][0][0] = GTX_T;
    christoffel[1][0][1] = GTX_X;
    christoffel[1][0][2] = 0.0;
    christoffel[1][0][3] = 0.0;
    christoffel[1][1][0] = GXX_T;
    christoffel[1][1][1] = GXX_X;
    christoffel[1][1][2] = 0.0;
    christoffel[1][1][3] = 0.0;
    christoffel[1][2][0] = 0.0;
    christoffel[1][2][1] = 0.0;
    christoffel[1][2][2] = GXt_t;
    christoffel[1][2][3] = 0.0;
    christoffel[1][3][0] = 0.0;
    christoffel[1][3][1] = 0.0;
    christoffel[1][3][2] = 0.0;
    christoffel[1][3][3] = GXp_p;
    christoffel[2][0][0] = 0.0;
    christoffel[2][0][1] = 0.0;
    christoffel[2][0][2] = GTt_t;
    christoffel[2][0][3] = 0.0;
    christoffel[2][1][0] = 0.0;
    christoffel[2][1][1] = 0.0;
    christoffel[2][1][2] = GXt_t;
    christoffel[2][1][3] = 0.0;
    christoffel[2][2][0] = Gtt_T;
    christoffel[2][2][1] = Gtt_X;
    christoffel[2][2][2] = 0.0;
    christoffel[2][2][3] = 0.0;
    christoffel[2][3][0] = 0.0;
    christoffel[2][3][1] = 0.0;
    christoffel[2][3][2] = 0.0;
    christoffel[2][3][3] = Gtp_p;
    christoffel[3][0][0] = 0.0;
    christoffel[3][0][1] = 0.0;
    christoffel[3][0][2] = 0.0;
    christoffel[3][0][3] = GTp_p;
    christoffel[3][1][0] = 0.0;
    christoffel[3][1][1] = 0.0;
    christoffel[3][1][2] = 0.0;
    christoffel[3][1][3] = GXp_p;
    christoffel[3][2][0] = 0.0;
    christoffel[3][2][1] = 0.0;
    christoffel[3][2][2] = 0.0;
    christoffel[3][2][3] = Gtp_p;
    christoffel[3][3][0] = Gpp_T;
    christoffel[3][3][1] = Gpp_X;
    christoffel[3][3][2] = Gpp_t;
    christoffel[3][3][3] = 0.0;

    return true;
}

void MetricKruskal::localToCoord(const double* pos, const double* ldir, double* dir, enum_nat_tetrad_type)
{
    double r = get_r(pos[0], pos[1]);
    double theta = pos[2];

    double fac = sqrt(r) / (2.0 * rs * sqrt(rs)) * exp(r / (2.0 * rs));

    dir[0] = ldir[0] * fac;
    dir[1] = ldir[1] * fac;
    dir[2] = ldir[2] / r;
    dir[3] = ldir[3] / (r * sin(theta));
}

void MetricKruskal::coordToLocal(const double* pos, const double* cdir, double* ldir, enum_nat_tetrad_type)
{
    double r = get_r(pos[0], pos[1]);
    double theta = pos[2];

    double fac = sqrt(r) / (2.0 * rs * sqrt(rs)) * exp(r / (2.0 * rs));

    ldir[0] = cdir[0] / fac;
    ldir[1] = cdir[1] / fac;
    ldir[2] = cdir[2] * r;
    ldir[3] = cdir[3] * r * sin(theta);
}

bool MetricKruskal::breakCondition(const double* pos)
{
    bool br = false;

    double T = pos[0];
    double X = pos[1];

    if (X * X - T * T < -1.0) {
        br = true;
    }

    return br;
}

bool MetricKruskal::calcDerivs(const double y[], double dydx[])
{
    double T = y[0];
    double X = y[1];
    double r = get_r(T, X);
    double theta = y[2];

    double GTT_T = T * rs * (r + rs) / (r * r) * exp(-r / rs);
    double GTX_X = GTT_T;
    double GXX_T = GTT_T;

    double GTT_X = -X * rs * (r + rs) / (r * r) * exp(-r / rs);
    double GTX_T = GTT_X;
    double GXX_X = GTT_X;

    double GTt_t = -2.0 * rs * rs / (r * r) * T * exp(-r / rs);
    double GTp_p = GTt_t;

    double GXt_t = 2.0 * rs * rs / (r * r) * X * exp(-r / rs);
    double GXp_p = GXt_t;

    double Gtt_T = -r / (2 * rs) * T;
    double Gtt_X = -r / (2 * rs) * X;
    double Gpp_T = -r / (2 * rs) * T * sin(theta) * sin(theta);
    double Gpp_X = -r / (2 * rs) * X * sin(theta) * sin(theta);
    double Gtp_p = cos(theta) / sin(theta);
    double Gpp_t = -sin(theta) * cos(theta);

    dydx[0] = y[4];
    dydx[1] = y[5];
    dydx[2] = y[6];
    dydx[3] = y[7];

    dydx[4] = -GTT_T * y[4] * y[4] - GXX_T * y[5] * y[5] - 2.0 * GTX_T * y[4] * y[5] - Gtt_T * y[6] * y[6]
        - Gpp_T * y[7] * y[7];
    dydx[5] = -2.0 * GTX_X * y[4] * y[5] - GTT_X * y[4] * y[4] - GXX_X * y[5] * y[5] - Gtt_X * y[6] * y[6]
        - Gpp_X * y[7] * y[7];
    dydx[6] = -2.0 * GTt_t * y[4] * y[6] - 2.0 * GXt_t * y[5] * y[6] - Gpp_t * y[7] * y[7];
    dydx[7] = -2.0 * GTp_p * y[4] * y[7] - 2.0 * GXp_p * y[5] * y[7] - 2.0 * Gtp_p * y[6] * y[7];

    return true;
}

double MetricKruskal::testConstraint(const double y[], const double kappa)
{
    double dT = y[4];
    double dX = y[5];
    double dtheta = y[6];
    double dphi = y[7];

    double r = get_r(y[0], y[1]);
    double fac = 4.0 * rs * rs * rs / r * exp(-r / rs);
    double theta = y[2];

    double sum = -kappa * mSign;
    sum += fac * (-dT * dT + dX * dX) + r * r * (dtheta * dtheta + sin(theta) * sin(theta) * dphi * dphi);
    return sum;
}

bool MetricKruskal::setParam(const char* pName, double val)
{
    Metric::setParam(pName, val);

    if (strcmp(pName, "mass") == 0) {
        mMass = val;
        rs = 2.0 * mMass;
    }

    return true;
}

bool MetricKruskal::report(const vec4, const vec4, char*& text)
{
    std::stringstream ss;
    ss << "Report for the  Kruskal metric\n\tcoordinate : (t,r,theta,phi)\n";
    ss << "---------------------------------------------------------------\n";
    ss.precision(DEF_FIXED_REPORT_PRECISION);
    ss.setf(std::ios::fixed);

    text = new char[ss.str().length() + 2];
    return CopyString(ss.str().c_str(), text);
}

int MetricKruskal::transToPseudoCart(vec4 p, vec4& cp)
{
    double T = p[0];
    double X = p[1];
    double r = get_r(T, X);
    // std::cerr << T << " " << X << " " << r << std::endl;
    double t = 0.0;
    if (r > rs) {
        t = 2.0 * rs * atanh(T / X);
    }
    else {
        t = 2.0 * rs * atanh(X / T);
    }

    vec4 np = vec4(t, r, p[2], p[3]);
    TransCoordinates::toCartesianCoord(mCoordType, np, cp);
    return 0;
}

void MetricKruskal::setStandardValues()
{
    mInitPos[0] = 0.0;
    mInitPos[1] = 6.0;
    mInitPos[2] = M_PI_2;
    mInitPos[3] = 0.0;
    mInitDir[0] = 1.0;
    mInitDir[1] = 0.0;
    mInitDir[2] = 0.0;

    mCoordNames[0] = std::string("T");
    mCoordNames[1] = std::string("X");
    mCoordNames[2] = std::string("theta");
    mCoordNames[3] = std::string("phi");
}

double MetricKruskal::get_r(double T, double X)
{
    // double a = (X * X - T * T) / exp(1.0);
    double a = (X + T) * (X - T) / exp(1.0);

    gsl_sf_result result;
    int status = gsl_sf_lambert_W0_e(a, &result);
    if (status != GSL_SUCCESS) {
        fprintf(stderr, "Error calculating LambertW func for a = %f\n", a);
    }

    double W = result.val;
    return rs * (W + 1.0);
}

} // end namespace m4d
