/**
 * @file    m4dMetricChazyCurzonRot.cpp
 * @author  Thomas Mueller
 *
 * This file is part of the m4d-library.
 */
#include "m4dMetricChazyCurzonRot.h"

namespace m4d {

MetricChazyCurzonRot::MetricChazyCurzonRot(double mass, double p)
{
    mMetricName = "ChazyCurzonRot";
    setCoordType(enum_coordinate_cylinder);

    mPhysicalUnits = enum_physical_constants_geom;
    mSpeedOfLight = 1.0;
    mGravConstant = 1.0;
    mDielectricPerm = 1.0;

    addParam("mass", mass);
    mMass = mass;
    addParam("p", p);
    m_p = p;
    m_q = sqrt(1.0 - m_p * m_p);

    mSign = 1.0;
    mLocTeds.push_back(enum_nat_tetrad_comoving);
    mLocTeds.push_back(enum_nat_tetrad_static);

    setStandardValues();
}

MetricChazyCurzonRot::~MetricChazyCurzonRot() {}

// *********************************** public methods ******************************

bool MetricChazyCurzonRot::calculateMetric(const double* pos)
{
    double rho = pos[1];
    double em2U, k, A;
    calcUkA(pos, em2U, k, A);

    // t1 = U(rho,z);
    // t2 = exp(t1);
    // t3 = t2*t2;
    double t3 = 1.0 / em2U;
    double t4 = A; // A(rho,z);
    double t5 = t3 * t4;
    double t6 = 1 / t3;
    double t7 = k; // k(rho,z);
    double t8 = exp(t7);
    double t9 = t8 * t8;
    double t10 = t6 * t9;
    double t11 = rho * rho;
    double t13 = t4 * t4;

    g_compts[0][0] = -t3;
    g_compts[0][1] = 0.0;
    g_compts[0][2] = -t5;
    g_compts[0][3] = 0.0;
    g_compts[1][0] = 0.0;
    g_compts[1][1] = t10;
    g_compts[1][2] = 0.0;
    g_compts[1][3] = 0.0;
    g_compts[2][0] = -t5;
    g_compts[2][1] = 0.0;
    g_compts[2][2] = t6 * t11 - t3 * t13;
    g_compts[2][3] = 0.0;
    g_compts[3][0] = 0.0;
    g_compts[3][1] = 0.0;
    g_compts[3][2] = 0.0;
    g_compts[3][3] = t10;

    return true;
}

bool MetricChazyCurzonRot::calculateChristoffels(const double* pos)
{
    double rho = pos[1];

    double em2U, k, A;
    calcUkA(pos, em2U, k, A);
    double U = -0.5 * log(em2U);

    double dUdrho, dUdz, dkdrho, dkdz, dAdrho, dAdz;
    calcDUka(pos, dUdrho, dUdz, dkdrho, dkdz, dAdrho, dAdz);

    double t1 = k; // k(rho,z);
    double t2 = exp(t1);
    double t3 = t2 * t2;
    double t4 = 1 / t3;
    double t5 = U; // U(rho,z);
    double t6 = exp(t5);
    double t7 = t6 * t6;
    double t8 = t7 * t7;
    double t9 = t4 * t8;
    double t10 = dUdrho; // diff(U(rho,z),rho);
    double t12 = dUdz; // diff(U(rho,z),z);
    double t14 = rho * rho;
    double t15 = t14 * t10;
    double t16 = 2.0 * t15;
    double t17 = A; // A(rho,z);
    double t18 = t8 * t17;
    double t19 = dAdrho; // diff(A(rho,z),rho);
    double t20 = t18 * t19;
    double t22 = 1 / t14;
    double t24 = (t16 + t20) * t22 / 2.0;
    double t27 = t8 * t19 * t22 / 2.0;
    double t33 = t8 * (2.0 * t17 * t10 + t19) * t4 / 2.0;
    double t36 = dAdz; // diff(A(rho,z),z);
    double t40 = t8 * (2.0 * t17 * t12 + t36) * t4 / 2.0;
    double t41 = t14 * t12;
    double t43 = t18 * t36;
    double t46 = (2.0 * t41 + t43) * t22 / 2.0;
    double t49 = t8 * t36 * t22 / 2.0;
    double t50 = dkdrho; // diff(k(rho,z),rho);
    double t51 = -t10 + t50;
    double t52 = dkdz; // diff(k(rho,z),z);
    double t53 = t12 - t52;
    double t54 = t14 * t17;
    double t58 = t17 * t17;
    double t59 = t8 * t58;
    double t65 = (4.0 * t54 * t10 + t14 * t19 + t59 * t19 - 2.0 * t17 * rho) * t22 / 2.0;
    double t69 = (t20 + t16 - 2.0 * rho) * t22 / 2.0;
    double t82 = (4.0 * t54 * t12 + t14 * t36 + t59 * t36) * t22 / 2.0;

    christoffel[0][0][0] = 0.0;
    christoffel[0][0][1] = t9 * t10;
    christoffel[0][0][2] = 0.0;
    christoffel[0][0][3] = t9 * t12;
    christoffel[0][1][0] = t24;
    christoffel[0][1][1] = 0.0;
    christoffel[0][1][2] = -t27;
    christoffel[0][1][3] = 0.0;
    christoffel[0][2][0] = 0.0;
    christoffel[0][2][1] = t33;
    christoffel[0][2][2] = 0.0;
    christoffel[0][2][3] = t40;
    christoffel[0][3][0] = t46;
    christoffel[0][3][1] = 0.0;
    christoffel[0][3][2] = -t49;
    christoffel[0][3][3] = 0.0;
    christoffel[1][0][0] = t24;
    christoffel[1][0][1] = 0.0;
    christoffel[1][0][2] = -t27;
    christoffel[1][0][3] = 0.0;
    christoffel[1][1][0] = 0.0;
    christoffel[1][1][1] = t51;
    christoffel[1][1][2] = 0.0;
    christoffel[1][1][3] = t53;
    christoffel[1][2][0] = t65;
    christoffel[1][2][1] = 0.0;
    christoffel[1][2][2] = -t69;
    christoffel[1][2][3] = 0.0;
    christoffel[1][3][0] = 0.0;
    christoffel[1][3][1] = -t53;
    christoffel[1][3][2] = 0.0;
    christoffel[1][3][3] = t51;
    christoffel[2][0][0] = 0.0;
    christoffel[2][0][1] = t33;
    christoffel[2][0][2] = 0.0;
    christoffel[2][0][3] = t40;
    christoffel[2][1][0] = t65;
    christoffel[2][1][1] = 0.0;
    christoffel[2][1][2] = -t69;
    christoffel[2][1][3] = 0.0;
    christoffel[2][2][0] = 0.0;
    christoffel[2][2][1] = t4 * (t15 - rho + t59 * t10 + t20);
    christoffel[2][2][2] = 0.0;
    christoffel[2][2][3] = t4 * (t41 + t59 * t12 + t43);
    christoffel[2][3][0] = t82;
    christoffel[2][3][1] = 0.0;
    christoffel[2][3][2] = -t46;
    christoffel[2][3][3] = 0.0;
    christoffel[3][0][0] = t46;
    christoffel[3][0][1] = 0.0;
    christoffel[3][0][2] = -t49;
    christoffel[3][0][3] = 0.0;
    christoffel[3][1][0] = 0.0;
    christoffel[3][1][1] = -t53;
    christoffel[3][1][2] = 0.0;
    christoffel[3][1][3] = t51;
    christoffel[3][2][0] = t82;
    christoffel[3][2][1] = 0.0;
    christoffel[3][2][2] = -t46;
    christoffel[3][2][3] = 0.0;
    christoffel[3][3][0] = 0.0;
    christoffel[3][3][1] = -t51;
    christoffel[3][3][2] = 0.0;
    christoffel[3][3][3] = -t53;

    return true;
}

bool MetricChazyCurzonRot::calculateChrisD(const double*)
{
    return true;
}

void MetricChazyCurzonRot::localToCoord(
    const double* pos, const double* ldir, double* dir, enum_nat_tetrad_type tetradType)
{
    double rho = pos[1];
    double em2U, k, A, eU, emk;
    calcUkA(pos, em2U, k, A);
    eU = sqrt(1.0 / em2U);
    emk = exp(-k);

    // fprintf(stderr,"%8.4f %8.4f %8.4f %8.4f\n",ldir[0],ldir[1],ldir[2],ldir[3]);

    if (tetradType == enum_nat_tetrad_static) {
        dir[0] = ldir[0] * sqrt(em2U) - A * eU / rho * ldir[2];
        dir[1] = ldir[1] * eU * emk;
        dir[2] = ldir[2] * eU / rho;
        dir[3] = ldir[3] * eU * emk;
    }
    else { // comoving
        double e2U = 1.0 / em2U;
        double g_phiphi = -A * A * e2U + rho * rho * em2U;
        double g_tphi = -A * e2U;
        double g_tt = -e2U;
        double Alpha = sqrt(g_phiphi / (g_tphi * g_tphi - g_tt * g_phiphi));
        double omega = -g_tphi / g_phiphi;

        dir[0] = Alpha * ldir[0];
        dir[1] = ldir[1] * eU * emk;
        dir[2] = Alpha * omega * ldir[0] + ldir[2] / sqrt(g_phiphi);
        dir[3] = ldir[3] * eU * emk;
    }
    fprintf(stderr, "%8.4g %8.4g %8.4g %8.4g\n\n", dir[0], dir[1], dir[2], dir[3]);
}

void MetricChazyCurzonRot::coordToLocal(const double*, const double*, double*, enum_nat_tetrad_type)
{
    fprintf(stderr, "uups... not implemented yet!\n");
    // TODO
}

bool MetricChazyCurzonRot::breakCondition(const double*)
{
    return false;
}

double MetricChazyCurzonRot::testConstraint(const double y[], const double kappa)
{
    double rho = y[1];

    double dt = y[4];
    double drho = y[5];
    double dphi = y[6];
    double dz = y[7];

    double em2U, k, A;
    calcUkA(y, em2U, k, A);

    double sum = -kappa * mSign;
    sum += em2U * (exp(2.0 * k) * (drho * drho + dz * dz) + rho * rho * dphi * dphi)
        - 1.0 / em2U * (dt + A * dphi) * (dt + A * dphi);

    return sum;
}

bool MetricChazyCurzonRot::setParam(const char* pName, double val)
{
    Metric::setParam(pName, val);

    if (strcmp(pName, "mass") == 0) {
        mMass = val;
    }
    else if (strcmp(pName, "p")) {
        m_p = val;
        m_q = sqrt(1.0 - m_p * m_p);
    }
    return true;
}

bool MetricChazyCurzonRot::report(const vec4, const vec4, char*& text)
{
    std::stringstream ss;
    ss << "Report for the  rotating ChazyCurzon metric\n\tWeyl coordinate : (t,rho,phi,z)\n";
    ss << "---------------------------------------------------------------\n";
    ss.precision(DEF_FIXED_REPORT_PRECISION);
    ss.setf(std::ios::fixed);

    text = new char[ss.str().length() + 2];
    return CopyString(ss.str().c_str(), text);
}

// ********************************* protected methods *****************************

void MetricChazyCurzonRot::setStandardValues()
{
    mInitPos[0] = 0.0;
    mInitPos[1] = 6.0;
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

void MetricChazyCurzonRot::calcUkA(const double* pos, double& em2U, double& k, double& A)
{
    double rho = pos[1];
    double z = pos[3];
    double r2 = rho * rho + z * z;
    double r = sqrt(r2);

    em2U = cosh(2.0 * mMass / r) - m_p * sinh(2.0 * mMass / r);
    k = -0.5 * mMass * mMass * rho * rho / (r2 * r2);
    A = 2.0 * m_q * mMass * z / r;
}

void MetricChazyCurzonRot::calcDUka(
    const double* pos, double& dUdrho, double& dUdz, double& dkdrho, double& dkdz, double& dAdrho, double& dAdz)
{
    double rho = pos[1];
    double z = pos[3];
    double r2 = rho * rho + z * z;
    double r = sqrt(r2);

    double ch = cosh(2.0 * mMass / r);
    double sh = sinh(2.0 * mMass / r);
    double w1 = r * r * r;

    dUdrho = mMass * rho * (sh - m_p * ch) / (w1 * (ch - m_p * sh));
    dUdz = mMass * z * (sh - m_p * ch) / (w1 * (ch - m_p * sh));

    dkdrho = mMass * mMass * rho * (rho * rho - z * z) / (r2 * r2 * r2);
    dkdz = 2.0 * mMass * mMass * rho * rho * z / (r2 * r2 * r2);

    dAdrho = -2.0 * mMass * m_q * rho * z / (r * r * r);
    dAdz = 2.0 * mMass * m_q * rho * rho / (r * r * r);
}

} // end namespace m4d
