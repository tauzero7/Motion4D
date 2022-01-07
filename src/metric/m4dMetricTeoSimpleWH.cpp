/**
 * @file    m4dMetricTeoSimpleWH.cpp
 * @author  Thomas Mueller
 *
 * This file is part of the m4d-library.
 */
#include "m4dMetricTeoSimpleWH.h"

namespace m4d {

MetricTeoSimpleWH::MetricTeoSimpleWH(double b0)
{
    mMetricName = "TeoSimpleWH";
    setCoordType(enum_coordinate_spherical);

    mPhysicalUnits = enum_physical_constants_geom;
    mSpeedOfLight = 1.0;
    mGravConstant = 1.0;

    setStandardValues();

    addParam("b0", b0);
    mb0 = b0;

    mLocTeds.push_back(enum_nat_tetrad_static);
    mLocTeds.push_back(enum_nat_tetrad_lnrf);

    mDrawTypes.push_back(enum_draw_embedding);

    if (!mEmbParam.empty()) {
        mEmbParam.clear();
    }
    mHaveEmbedding = true;

    mEmb_lmin = -10.0;
    mEmb_lmax = 10.0;
    mEmb_l_num = 40.0;
    mEmb_phi_num = 40.0;
    mEmb_lstep = (mEmb_lmax - mEmb_lmin) / mEmb_l_num;
    mEmb_phistep = 2.0 * M_PI / mEmb_phi_num;
    addEmbeddingParam("emb_lmin", mEmb_lmin);
    addEmbeddingParam("emb_lmax", mEmb_lmax);
    addEmbeddingParam("emb_l_num", mEmb_l_num);
    addEmbeddingParam("emb_phi_num", mEmb_phi_num);
}

MetricTeoSimpleWH::~MetricTeoSimpleWH() {}

// *********************************** public methods ******************************

bool MetricTeoSimpleWH::calculateMetric(const double* pos)
{
    double l = pos[1];
    double theta = pos[2];
    double c = 1.0;
    double b0 = mb0;
    double ch2 = 1.0;
    double b0h2 = b0 * b0;
    double b0h4 = b0h2 * b0h2;
    double l2 = l * l;
    double st = sin(theta);
    double st2 = st * st;
    double t4 = pow(b0h2 + l2, -1.5);
    double t5 = pow(b0h2 + l2, -2.0);

    g_compts[0][0] = 0.25 * b0h4 * ch2 * t5 * st2 - ch2;
    g_compts[0][1] = 0;
    g_compts[0][2] = 0;
    g_compts[0][3] = b0h2 * c * (-0.5 * b0h2 - 0.5 * l2) * t4 * st2;
    g_compts[1][0] = 0;
    g_compts[1][1] = 1.0;
    g_compts[1][2] = 0;
    g_compts[1][3] = 0;
    g_compts[2][0] = 0;
    g_compts[2][1] = 0;
    g_compts[2][2] = b0h2 + l2;
    g_compts[2][3] = 0;
    g_compts[3][0] = b0h2 * c * (-0.5 * b0h2 - 0.5 * l2) * t4 * st2;
    g_compts[3][1] = 0;
    g_compts[3][2] = 0;
    g_compts[3][3] = (b0h2 + l2) * st2;

    return true;
}

bool MetricTeoSimpleWH::calculateChristoffels(const double* pos)
{
    double l = pos[1];
    double theta = pos[2];
    double c = 1.0;
    double b0 = mb0;
    double b0h2 = b0 * b0;
    double b0h4 = b0h2 * b0h2;
    double ch3 = 1.0;
    double ch2 = 1.0;
    double ch4 = 1.0;
    double ch5 = 1.0;
    double l2 = l * l;
    double st = sin(theta);
    double st2 = st * st;
    double st4 = st2 * st2;
    double st3 = st * st2;
    double ct = cos(theta);
    double r2 = pow(b0h2 + l2, 2.0);
    double invr3 = pow(b0h2 + l2, -3.0);
    double t1 = pow(b0h2 + l2, -0.5);
    double t2 = pow(b0h2 + l2, -3.5);
    double t3 = pow(b0h2 + l2, -2.5);
    double t4 = pow(b0h2 + l2, -1.5);
    double t5 = pow(b0h2 + l2, -2.0);
    double t6 = pow(-0.5 * b0h2 - 0.5 * l2, 2);
    double b0h6 = pow(b0, 6.0);

    christoffel[0][0][0] = 0;
    christoffel[0][0][1] = 0.5 * b0h4 * ch2 * l * invr3 * st2;
    christoffel[0][0][2] = -0.25 * b0h4 * ch2 * invr3 * st * ct;
    christoffel[0][0][3] = 0;
    christoffel[0][1][0] = -1.0L / 2.0L * b0h2 * ch3 * (-0.5 * b0h2 - 0.5 * l2) * t1 * (0.25 * b0h4 * st2 - r2)
            * (-3.0 * b0h2 * c * l * (-0.5 * b0h2 - 0.5 * l2) * t3 * st2 - 1.0 * b0h2 * c * l * t4 * st2) * st2
            / ((0.25 * b0h4 * ch2 * t5 * st2 - ch2)
                * (-0.25 * b0h4 * ch2 * r2 * st4 + ch2 * r2 * (0.25 * b0h4 * st2 - r2) * st2))
        - 1.0 * b0h4 * ch2 * l * invr3
            * ((1.0L / 2.0L) * b0h4 * ch4 * t6 * t5 * (0.25 * b0h4 * st2 - r2) * st4
                    / (pow(0.25 * b0h4 * ch2 * t5 * st2 - ch2, 2)
                        * (-0.25 * b0h4 * ch2 * r2 * st4 + ch2 * r2 * (0.25 * b0h4 * st2 - r2) * st2))
                + (1.0L / 2.0L) / (0.25 * b0h4 * ch2 * t5 * st2 - ch2))
            * st2;
    christoffel[0][1][1] = 0;
    christoffel[0][1][2] = 0;
    christoffel[0][1][3] = 0.5 * b0h6 * ch5 * l * (-0.5 * b0h2 - 0.5 * l2) * t2 * (0.25 * b0h4 * st2 - r2) * st4
            / ((0.25 * b0h4 * ch2 * t5 * st2 - ch2)
                * (-0.25 * b0h4 * ch2 * r2 * st4 + ch2 * r2 * (0.25 * b0h4 * st2 - r2) * st2))
        + (1.0L / 2.0L) * ch2 * pow(b0h2 + l2, 1.0) * (0.25 * b0h4 * st2 - r2)
            * (-3.0 * b0h2 * c * l * (-0.5 * b0h2 - 0.5 * l2) * t3 * st2 - 1.0 * b0h2 * c * l * t4 * st2)
            / (-0.25 * b0h4 * ch2 * r2 * st4 + ch2 * r2 * (0.25 * b0h4 * st2 - r2) * st2);
    christoffel[0][2][0] = 0.5 * b0h4 * ch2 * t5
            * ((1.0L / 2.0L) * b0h4 * ch4 * t6 * t5 * (0.25 * b0h4 * st2 - r2) * st4
                    / (pow(0.25 * b0h4 * ch2 * t5 * st2 - ch2, 2)
                        * (-0.25 * b0h4 * ch2 * r2 * st4 + ch2 * r2 * (0.25 * b0h4 * st2 - r2) * st2))
                + (1.0L / 2.0L) / (0.25 * b0h4 * ch2 * t5 * st2 - ch2))
            * st * ct
        - 1.0 * b0h4 * ch4 * t6 * t5 * (0.25 * b0h4 * st2 - r2) * st3 * ct
            / ((0.25 * b0h4 * ch2 * t5 * st2 - ch2)
                * (-0.25 * b0h4 * ch2 * r2 * st4 + ch2 * r2 * (0.25 * b0h4 * st2 - r2) * st2));
    christoffel[0][2][1] = 0;
    christoffel[0][2][2] = 0;
    christoffel[0][2][3] = 1.0 * b0h2 * ch3 * (-0.5 * b0h2 - 0.5 * l2) * t1 * (0.25 * b0h4 * st2 - r2) * st * ct
            / (-0.25 * b0h4 * ch2 * r2 * st4 + ch2 * r2 * (0.25 * b0h4 * st2 - r2) * st2)
        - 0.25 * b0h6 * ch5 * (-0.5 * b0h2 - 0.5 * l2) * t3 * (0.25 * b0h4 * st2 - r2) * st3 * ct
            / ((0.25 * b0h4 * ch2 * t5 * st2 - ch2)
                * (-0.25 * b0h4 * ch2 * r2 * st4 + ch2 * r2 * (0.25 * b0h4 * st2 - r2) * st2));
    christoffel[0][3][0] = 0;
    christoffel[0][3][1] = 1.5 * b0h2 * c * l * (-0.5 * b0h2 - 0.5 * l2) * t3 * st2 + 0.5 * b0h2 * c * l * t4 * st2;
    christoffel[0][3][2] = -1.0 * b0h2 * c * (-0.5 * b0h2 - 0.5 * l2) * t3 * st * ct;
    christoffel[0][3][3] = 0;
    christoffel[1][0][0] = -1.0L / 2.0L * b0h2 * ch3 * (-0.5 * b0h2 - 0.5 * l2) * t1 * (0.25 * b0h4 * st2 - r2)
            * (-3.0 * b0h2 * c * l * (-0.5 * b0h2 - 0.5 * l2) * t3 * st2 - 1.0 * b0h2 * c * l * t4 * st2) * st2
            / ((0.25 * b0h4 * ch2 * t5 * st2 - ch2)
                * (-0.25 * b0h4 * ch2 * r2 * st4 + ch2 * r2 * (0.25 * b0h4 * st2 - r2) * st2))
        - 1.0 * b0h4 * ch2 * l * invr3
            * ((1.0L / 2.0L) * b0h4 * ch4 * t6 * t5 * (0.25 * b0h4 * st2 - r2) * st4
                    / (pow(0.25 * b0h4 * ch2 * t5 * st2 - ch2, 2)
                        * (-0.25 * b0h4 * ch2 * r2 * st4 + ch2 * r2 * (0.25 * b0h4 * st2 - r2) * st2))
                + (1.0L / 2.0L) / (0.25 * b0h4 * ch2 * t5 * st2 - ch2))
            * st2;
    christoffel[1][0][1] = 0;
    christoffel[1][0][2] = 0;
    christoffel[1][0][3] = 0.5 * b0h6 * ch5 * l * (-0.5 * b0h2 - 0.5 * l2) * t2 * (0.25 * b0h4 * st2 - r2) * st4
            / ((0.25 * b0h4 * ch2 * t5 * st2 - ch2)
                * (-0.25 * b0h4 * ch2 * r2 * st4 + ch2 * r2 * (0.25 * b0h4 * st2 - r2) * st2))
        + (1.0L / 2.0L) * ch2 * pow(b0h2 + l2, 1.0) * (0.25 * b0h4 * st2 - r2)
            * (-3.0 * b0h2 * c * l * (-0.5 * b0h2 - 0.5 * l2) * t3 * st2 - 1.0 * b0h2 * c * l * t4 * st2)
            / (-0.25 * b0h4 * ch2 * r2 * st4 + ch2 * r2 * (0.25 * b0h4 * st2 - r2) * st2);
    christoffel[1][1][0] = 0;
    christoffel[1][1][1] = 0;
    christoffel[1][1][2] = 0;
    christoffel[1][1][3] = 0;
    christoffel[1][2][0] = 0;
    christoffel[1][2][1] = 0;
    christoffel[1][2][2] = 1.0 * l / (b0h2 + l2);
    christoffel[1][2][3] = 0;
    christoffel[1][3][0] = -1.0 * b0h2 * ch3 * l * (-0.5 * b0h2 - 0.5 * l2) * t1 * (0.25 * b0h4 * st2 - r2) * st4
            / ((0.25 * b0h4 * ch2 * t5 * st2 - ch2)
                * (-0.25 * b0h4 * ch2 * r2 * st4 + ch2 * r2 * (0.25 * b0h4 * st2 - r2) * st2))
        + (-3.0 * b0h2 * c * l * (-0.5 * b0h2 - 0.5 * l2) * t3 * st2 - 1.0 * b0h2 * c * l * t4 * st2)
            * ((1.0L / 2.0L) * b0h4 * ch4 * t6 * t5 * (0.25 * b0h4 * st2 - r2) * st4
                    / (pow(0.25 * b0h4 * ch2 * t5 * st2 - ch2, 2)
                        * (-0.25 * b0h4 * ch2 * r2 * st4 + ch2 * r2 * (0.25 * b0h4 * st2 - r2) * st2))
                + (1.0L / 2.0L) / (0.25 * b0h4 * ch2 * t5 * st2 - ch2));
    christoffel[1][3][1] = 0;
    christoffel[1][3][2] = 0;
    christoffel[1][3][3] = -1.0L / 2.0L * b0h2 * ch3 * (-0.5 * b0h2 - 0.5 * l2) * t1 * (0.25 * b0h4 * st2 - r2)
            * (-3.0 * b0h2 * c * l * (-0.5 * b0h2 - 0.5 * l2) * t3 * st2 - 1.0 * b0h2 * c * l * t4 * st2) * st2
            / ((0.25 * b0h4 * ch2 * t5 * st2 - ch2)
                * (-0.25 * b0h4 * ch2 * r2 * st4 + ch2 * r2 * (0.25 * b0h4 * st2 - r2) * st2))
        + 1.0 * ch2 * l * pow(b0h2 + l2, 1.0) * (0.25 * b0h4 * st2 - r2) * st2
            / (-0.25 * b0h4 * ch2 * r2 * st4 + ch2 * r2 * (0.25 * b0h4 * st2 - r2) * st2);
    christoffel[2][0][0] = 0.5 * b0h4 * ch2 * t5
            * ((1.0L / 2.0L) * b0h4 * ch4 * t6 * t5 * (0.25 * b0h4 * st2 - r2) * st4
                    / (pow(0.25 * b0h4 * ch2 * t5 * st2 - ch2, 2)
                        * (-0.25 * b0h4 * ch2 * r2 * st4 + ch2 * r2 * (0.25 * b0h4 * st2 - r2) * st2))
                + (1.0L / 2.0L) / (0.25 * b0h4 * ch2 * t5 * st2 - ch2))
            * st * ct
        - 1.0 * b0h4 * ch4 * t6 * t5 * (0.25 * b0h4 * st2 - r2) * st3 * ct
            / ((0.25 * b0h4 * ch2 * t5 * st2 - ch2)
                * (-0.25 * b0h4 * ch2 * r2 * st4 + ch2 * r2 * (0.25 * b0h4 * st2 - r2) * st2));
    christoffel[2][0][1] = 0;
    christoffel[2][0][2] = 0;
    christoffel[2][0][3] = 1.0 * b0h2 * ch3 * (-0.5 * b0h2 - 0.5 * l2) * t1 * (0.25 * b0h4 * st2 - r2) * st * ct
            / (-0.25 * b0h4 * ch2 * r2 * st4 + ch2 * r2 * (0.25 * b0h4 * st2 - r2) * st2)
        - 0.25 * b0h6 * ch5 * (-0.5 * b0h2 - 0.5 * l2) * t3 * (0.25 * b0h4 * st2 - r2) * st3 * ct
            / ((0.25 * b0h4 * ch2 * t5 * st2 - ch2)
                * (-0.25 * b0h4 * ch2 * r2 * st4 + ch2 * r2 * (0.25 * b0h4 * st2 - r2) * st2));
    christoffel[2][1][0] = 0;
    christoffel[2][1][1] = 0;
    christoffel[2][1][2] = 1.0 * l / (b0h2 + l2);
    christoffel[2][1][3] = 0;
    christoffel[2][2][0] = 0;
    christoffel[2][2][1] = -1.0 * l;
    christoffel[2][2][2] = 0;
    christoffel[2][2][3] = 0;
    christoffel[2][3][0] = 2.0 * b0h2 * c * (-0.5 * b0h2 - 0.5 * l2) * t4
            * ((1.0L / 2.0L) * b0h4 * ch4 * t6 * t5 * (0.25 * b0h4 * st2 - r2) * st4
                    / (pow(0.25 * b0h4 * ch2 * t5 * st2 - ch2, 2)
                        * (-0.25 * b0h4 * ch2 * r2 * st4 + ch2 * r2 * (0.25 * b0h4 * st2 - r2) * st2))
                + (1.0L / 2.0L) / (0.25 * b0h4 * ch2 * t5 * st2 - ch2))
            * st * ct
        - 1.0 * b0h2 * ch3 * (-0.5 * b0h2 - 0.5 * l2) * sqrt(b0h2 + l2) * (0.25 * b0h4 * st2 - r2) * st3 * ct
            / ((0.25 * b0h4 * ch2 * t5 * st2 - ch2)
                * (-0.25 * b0h4 * ch2 * r2 * st4 + ch2 * r2 * (0.25 * b0h4 * st2 - r2) * st2));
    christoffel[2][3][1] = 0;
    christoffel[2][3][2] = 0;
    christoffel[2][3][3] = -1.0 * b0h4 * ch4 * t6 * t5 * (0.25 * b0h4 * st2 - r2) * st3 * ct
            / ((0.25 * b0h4 * ch2 * t5 * st2 - ch2)
                * (-0.25 * b0h4 * ch2 * r2 * st4 + ch2 * r2 * (0.25 * b0h4 * st2 - r2) * st2))
        + 1.0 * ch2 * r2 * (0.25 * b0h4 * st2 - r2) * st * ct
            / (-0.25 * b0h4 * ch2 * r2 * st4 + ch2 * r2 * (0.25 * b0h4 * st2 - r2) * st2);
    christoffel[3][0][0] = 0;
    christoffel[3][0][1] = 1.5 * b0h2 * c * l * (-0.5 * b0h2 - 0.5 * l2) * t3 * st2 + 0.5 * b0h2 * c * l * t4 * st2;
    christoffel[3][0][2] = -1.0 * b0h2 * c * (-0.5 * b0h2 - 0.5 * l2) * t3 * st * ct;
    christoffel[3][0][3] = 0;
    christoffel[3][1][0] = -1.0 * b0h2 * ch3 * l * (-0.5 * b0h2 - 0.5 * l2) * t1 * (0.25 * b0h4 * st2 - r2) * st4
            / ((0.25 * b0h4 * ch2 * t5 * st2 - ch2)
                * (-0.25 * b0h4 * ch2 * r2 * st4 + ch2 * r2 * (0.25 * b0h4 * st2 - r2) * st2))
        + (-3.0 * b0h2 * c * l * (-0.5 * b0h2 - 0.5 * l2) * t3 * st2 - 1.0 * b0h2 * c * l * t4 * st2)
            * ((1.0L / 2.0L) * b0h4 * ch4 * t6 * t5 * (0.25 * b0h4 * st2 - r2) * st4
                    / (pow(0.25 * b0h4 * ch2 * t5 * st2 - ch2, 2)
                        * (-0.25 * b0h4 * ch2 * r2 * st4 + ch2 * r2 * (0.25 * b0h4 * st2 - r2) * st2))
                + (1.0L / 2.0L) / (0.25 * b0h4 * ch2 * t5 * st2 - ch2));
    christoffel[3][1][1] = 0;
    christoffel[3][1][2] = 0;
    christoffel[3][1][3] = -1.0L / 2.0L * b0h2 * ch3 * (-0.5 * b0h2 - 0.5 * l2) * t1 * (0.25 * b0h4 * st2 - r2)
            * (-3.0 * b0h2 * c * l * (-0.5 * b0h2 - 0.5 * l2) * t3 * st2 - 1.0 * b0h2 * c * l * t4 * st2) * st2
            / ((0.25 * b0h4 * ch2 * t5 * st2 - ch2)
                * (-0.25 * b0h4 * ch2 * r2 * st4 + ch2 * r2 * (0.25 * b0h4 * st2 - r2) * st2))
        + 1.0 * ch2 * l * pow(b0h2 + l2, 1.0) * (0.25 * b0h4 * st2 - r2) * st2
            / (-0.25 * b0h4 * ch2 * r2 * st4 + ch2 * r2 * (0.25 * b0h4 * st2 - r2) * st2);
    christoffel[3][2][0] = 2.0 * b0h2 * c * (-0.5 * b0h2 - 0.5 * l2) * t4
            * ((1.0L / 2.0L) * b0h4 * ch4 * t6 * t5 * (0.25 * b0h4 * st2 - r2) * st4
                    / (pow(0.25 * b0h4 * ch2 * t5 * st2 - ch2, 2)
                        * (-0.25 * b0h4 * ch2 * r2 * st4 + ch2 * r2 * (0.25 * b0h4 * st2 - r2) * st2))
                + (1.0L / 2.0L) / (0.25 * b0h4 * ch2 * t5 * st2 - ch2))
            * st * ct
        - 1.0 * b0h2 * ch3 * (-0.5 * b0h2 - 0.5 * l2) * sqrt(b0h2 + l2) * (0.25 * b0h4 * st2 - r2) * st3 * ct
            / ((0.25 * b0h4 * ch2 * t5 * st2 - ch2)
                * (-0.25 * b0h4 * ch2 * r2 * st4 + ch2 * r2 * (0.25 * b0h4 * st2 - r2) * st2));
    christoffel[3][2][1] = 0;
    christoffel[3][2][2] = 0;
    christoffel[3][2][3] = -1.0 * b0h4 * ch4 * t6 * t5 * (0.25 * b0h4 * st2 - r2) * st3 * ct
            / ((0.25 * b0h4 * ch2 * t5 * st2 - ch2)
                * (-0.25 * b0h4 * ch2 * r2 * st4 + ch2 * r2 * (0.25 * b0h4 * st2 - r2) * st2))
        + 1.0 * ch2 * r2 * (0.25 * b0h4 * st2 - r2) * st * ct
            / (-0.25 * b0h4 * ch2 * r2 * st4 + ch2 * r2 * (0.25 * b0h4 * st2 - r2) * st2);
    christoffel[3][3][0] = 0;
    christoffel[3][3][1] = -1.0 * l * st2;
    christoffel[3][3][2] = -1.0 * st * ct;
    christoffel[3][3][3] = 0;

    return true;
}

void MetricTeoSimpleWH::localToCoord(const double* pos, const double* ldir, double* dir, enum_nat_tetrad_type type)
{
    double l = pos[1];
    double theta = pos[2];
    double c = 1.0;
    double r = sqrt(l * l + mb0 * mb0);
    double st = sin(theta);

    if (type == enum_nat_tetrad_static) {
    }
    else {
        dir[0] = ldir[0] / c;
        dir[1] = ldir[1];
        dir[2] = ldir[2] / r;
        dir[3] = ldir[0] * 0.5 * mb0 * mb0 * pow(r, -3.0) + ldir[3] / (st * r);
    }
}

void MetricTeoSimpleWH::coordToLocal(const double*, const double*, double*, enum_nat_tetrad_type)
{
    fprintf(stderr, "uups... not implemented yet!\n");
    // TODO
}

bool MetricTeoSimpleWH::breakCondition(const double*)
{
    return false;
}

int MetricTeoSimpleWH::transToPseudoCart(vec4 p, vec4& cp)
{
    TransCoordinates::toCartesianCoord(mCoordType, p, cp);
    if (p[1] > 0) {
        return 0;
    }

    return 1;
}

bool MetricTeoSimpleWH::setParam(const char* pName, double val)
{
    bool ok = Metric::setParam(pName, val);
    if (ok) {
        mb0 = val;
    }
    return ok;
}

bool MetricTeoSimpleWH::report(const vec4, const vec4, char*& text)
{
    std::stringstream ss;
    ss << "Report for Teo wormhole metric\n\tcoordinate : (t,l,theta,phi)\n";
    ss << "---------------------------------------------------------------\n";
    ss << "  physical units ................................. no\n";
    ss.precision(DEF_FIXED_REPORT_PRECISION);
    ss.setf(std::ios::fixed);

    text = new char[ss.str().length() + 2];
    return CopyString(ss.str().c_str(), text);
}

bool MetricTeoSimpleWH::transToEmbedding(vec4 p, vec4& ep)
{
    if (mb0 <= 0.0) {
        return false;
    }

    vec4 cp;
    vec4 propCoords = p;
    propCoords[1] = sqrt(mb0 * mb0 + p[1] * p[1]);
    transToPseudoCart(propCoords, cp);

    double r = propCoords[1];
    double x = cp[1];
    double y = cp[2];
    double z;

    z = getShapeVal(r);
    if (p[1] < 0.0) {
        z *= -1.0;
    }

    ep = vec4(p[0], x, y, z);
    return true;
}

double MetricTeoSimpleWH::getShapeVal(double r)
{
    return mb0 * log((sqrt(r * r - mb0 * mb0) + r) / mb0);
}

bool MetricTeoSimpleWH::setEmbeddingParam(const char* name, double val)
{
    Metric::setEmbeddingParam(name, val);

    if (strcmp(name, "emb_lmin") == 0) {
        mEmb_lmin = val;
    }
    else if (strcmp(name, "emb_lmax") == 0) {
        mEmb_lmax = val;
    }
    else if (strcmp(name, "emb_l_num") == 0) {
        mEmb_l_num = val;
        if (mEmb_l_num < 5.0) {
            mEmb_l_num = 5.0;
        }
    }
    else if (strcmp(name, "emb_phi_num") == 0) {
        mEmb_phi_num = val;
        if (mEmb_phi_num < 4) {
            mEmb_phi_num = 4;
        }
    }
    return true;
}

// ********************************* protected methods *****************************

void MetricTeoSimpleWH::setStandardValues()
{
    mInitPos[0] = 0.0;
    mInitPos[1] = 10.0;
    mInitPos[2] = M_PI_2;
    mInitPos[3] = 0.0;
    mInitDir[0] = 1.0;
    mInitDir[1] = 0.0;
    mInitDir[2] = 0.0;

    mCoordNames[0] = std::string("t");
    mCoordNames[1] = std::string("l");
    mCoordNames[2] = std::string("theta");
    mCoordNames[3] = std::string("phi");
}

} // end namespace m4d
