// -------------------------------------------------------------------------------
/*
   m4dMetricSchwarzschildTortoise.cpp

  Copyright (c) 2010-2014  Thomas Mueller


   This file is part of the m4d-library.

   The m4d-library is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   The m4d-library is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with the m4d-library.  If not, see <http://www.gnu.org/licenses/>.

*/
// -------------------------------------------------------------------------------

#include "m4dMetricSchwarzschildTortoise.h"

namespace m4d {

#define eps 1.0e-6


/*! Standard constructor for the Schwarzschild metric.
 *
 * \param  mass : mass of the black hole.
 */
MetricSchwarzschildTortoise::MetricSchwarzschildTortoise(double mass) {
    mMetricName  = "SchwarzschildTortoise";
    mCoordType   = enum_coordinate_custom;

    mPhysicalUnits = enum_physical_constants_geom;
    mSpeedOfLight = 1.0;
    mGravConstant = 1.0;

    addParam("mass", mass);
    mMass = mass;
    rs = 2.0 * mGravConstant * mMass / (mSpeedOfLight * mSpeedOfLight);

    /*  Only a static tetrad is defined  */
    mLocTeds.push_back(enum_nat_tetrad_static);

    mDrawTypes.push_back(enum_draw_embedding);

    /*  parameters for the embedding diagram  */
    if (!mEmbParam.empty()) {
        mEmbParam.clear();
    }
    mHaveEmbedding = true;

    mEmb_rmin    = rs;
    mEmb_rmax    = 5.0 * rs;
    mEmb_r_num   = 20.0;
    mEmb_phi_num = 40.0;
    mEmb_rstep = (mEmb_rmax - mEmb_rmin) / mEmb_r_num;
    mEmb_phistep = 2.0 * M_PI / mEmb_phi_num;
    addEmbeddingParam("emb_rmin", mEmb_rmin);
    addEmbeddingParam("emb_rmax", mEmb_rmax);
    addEmbeddingParam("emb_r_num", 20.0);
    addEmbeddingParam("emb_phi_num", 40.0);

    setStandardValues();
    gsl_set_error_handler_off();
}

/*!
 */
MetricSchwarzschildTortoise::~MetricSchwarzschildTortoise() {
}


// *********************************** public methods ******************************

/*! Calculate the contravariant metric components at position 'pos'.
 *
 *  \param pos : pointer to position.
 */
bool MetricSchwarzschildTortoise::calculateMetric(const double* pos) {
    double rho   = pos[1];
    double theta = pos[2];

    double c = mSpeedOfLight;

    double t4 = exp((rho - rs) / rs);
    double t5 = gsl_sf_lambert_W0(t4); // LambertW(t4);
    double t6 = 1.0 + t5;
    double t8 = t5 / t6;
    double t9 = c * c;
    double t11 = rs * rs;
    double t12 = t6 * t6;
    double t13 = t11 * t12;
    double t14 = sin(theta);
    double t15 = t14 * t14;

    g_compts[0][0] = -t8 * t9;
    g_compts[0][1] = 0.0;
    g_compts[0][2] = 0.0;
    g_compts[0][3] = 0.0;
    g_compts[1][0] = 0.0;
    g_compts[1][1] = t8;
    g_compts[1][2] = 0.0;
    g_compts[1][3] = 0.0;
    g_compts[2][0] = 0.0;
    g_compts[2][1] = 0.0;
    g_compts[2][2] = t13;
    g_compts[2][3] = 0.0;
    g_compts[3][0] = 0.0;
    g_compts[3][1] = 0.0;
    g_compts[3][2] = 0.0;
    g_compts[3][3] = t13 * t15;

    return true;
}

/*! Calculate the Christoffel symbols of the second kind at position 'pos'.
 *
 *  \param pos : pointer to position.
 */
bool MetricSchwarzschildTortoise::calculateChristoffels(const double* pos) {
    double rho   = pos[1];
    double theta = pos[2];

    double c = mSpeedOfLight;

    double t2 = 1 / rs;
    double t4 = exp((rho - rs) * t2);
    double t5 = gsl_sf_lambert_W0(t4); // LambertW(t4);
    double t6 = 1.0 + t5;
    double t7 = t6 * t6;
    double t8 = 1 / t7;
    double t9 = c * c;
    double t14 = t8 * t2 / 2.0;
    double t16 = t5 * t8 * t2;
    double t17 = rs * t6;
    double t18 = sin(theta);
    double t20 = cos(theta);
    double t21 = 1 / t18 * t20;
    double t22 = t18 * t18;

    christoffel[0][0][0] = 0.0;
    christoffel[0][0][1] = t8 * t9 * t2 / 2.0;
    christoffel[0][0][2] = 0.0;
    christoffel[0][0][3] = 0.0;
    christoffel[0][1][0] = t14;
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
    christoffel[1][0][0] = t14;
    christoffel[1][0][1] = 0.0;
    christoffel[1][0][2] = 0.0;
    christoffel[1][0][3] = 0.0;
    christoffel[1][1][0] = 0.0;
    christoffel[1][1][1] = t14;
    christoffel[1][1][2] = 0.0;
    christoffel[1][1][3] = 0.0;
    christoffel[1][2][0] = 0.0;
    christoffel[1][2][1] = 0.0;
    christoffel[1][2][2] = t16;
    christoffel[1][2][3] = 0.0;
    christoffel[1][3][0] = 0.0;
    christoffel[1][3][1] = 0.0;
    christoffel[1][3][2] = 0.0;
    christoffel[1][3][3] = t16;
    christoffel[2][0][0] = 0.0;
    christoffel[2][0][1] = 0.0;
    christoffel[2][0][2] = 0.0;
    christoffel[2][0][3] = 0.0;
    christoffel[2][1][0] = 0.0;
    christoffel[2][1][1] = 0.0;
    christoffel[2][1][2] = t16;
    christoffel[2][1][3] = 0.0;
    christoffel[2][2][0] = 0.0;
    christoffel[2][2][1] = -t17;
    christoffel[2][2][2] = 0.0;
    christoffel[2][2][3] = 0.0;
    christoffel[2][3][0] = 0.0;
    christoffel[2][3][1] = 0.0;
    christoffel[2][3][2] = 0.0;
    christoffel[2][3][3] = t21;
    christoffel[3][0][0] = 0.0;
    christoffel[3][0][1] = 0.0;
    christoffel[3][0][2] = 0.0;
    christoffel[3][0][3] = 0.0;
    christoffel[3][1][0] = 0.0;
    christoffel[3][1][1] = 0.0;
    christoffel[3][1][2] = 0.0;
    christoffel[3][1][3] = t16;
    christoffel[3][2][0] = 0.0;
    christoffel[3][2][1] = 0.0;
    christoffel[3][2][2] = 0.0;
    christoffel[3][2][3] = t21;
    christoffel[3][3][0] = 0.0;
    christoffel[3][3][1] = -t17 * t22;
    christoffel[3][3][2] = -t18 * t20;
    christoffel[3][3][3] = 0.0;

    return true;
}

/*! Calculate Jacobi matrix.
 *
 *  \param pos : pointer to position.
 */
bool MetricSchwarzschildTortoise::calculateChrisD(const double* pos) {
    double rho   = pos[1];
    double theta = pos[2];

    double c = mSpeedOfLight;

    double t4 = exp((rho - rs) / rs);
    double t5 = gsl_sf_lambert_W0(t4); // LambertW(t4);
    double t6 = 1.0 + t5;
    double t7 = t6 * t6;
    double t8 = t7 * t7;
    double t9 = 1 / t8;
    double t10 = c * c;
    double t12 = rs * rs;
    double t13 = 1 / t12;
    double t16 = t9 * t13;
    double t17 = t16 * t5;
    double t20 = t5 * (-1.0 + t5) * t16;
    double t21 = 1 / t6;
    double t23 = cos(theta);
    double t24 = t23 * t23;
    double t25 = -1.0 + t24;
    double t26 = 1 / t25;
    double t30 = sin(theta);

    chrisD[0][0][0][0] = 0.0;
    chrisD[0][0][0][1] = 0.0;
    chrisD[0][0][0][2] = 0.0;
    chrisD[0][0][0][3] = 0.0;
    chrisD[0][0][1][0] = 0.0;
    chrisD[0][0][1][1] = -t9 * t10 * t13 * t5;
    chrisD[0][0][1][2] = 0.0;
    chrisD[0][0][1][3] = 0.0;
    chrisD[0][0][2][0] = 0.0;
    chrisD[0][0][2][1] = 0.0;
    chrisD[0][0][2][2] = 0.0;
    chrisD[0][0][2][3] = 0.0;
    chrisD[0][0][3][0] = 0.0;
    chrisD[0][0][3][1] = 0.0;
    chrisD[0][0][3][2] = 0.0;
    chrisD[0][0][3][3] = 0.0;
    chrisD[0][1][0][0] = 0.0;
    chrisD[0][1][0][1] = -t17;
    chrisD[0][1][0][2] = 0.0;
    chrisD[0][1][0][3] = 0.0;
    chrisD[0][1][1][0] = 0.0;
    chrisD[0][1][1][1] = 0.0;
    chrisD[0][1][1][2] = 0.0;
    chrisD[0][1][1][3] = 0.0;
    chrisD[0][1][2][0] = 0.0;
    chrisD[0][1][2][1] = 0.0;
    chrisD[0][1][2][2] = 0.0;
    chrisD[0][1][2][3] = 0.0;
    chrisD[0][1][3][0] = 0.0;
    chrisD[0][1][3][1] = 0.0;
    chrisD[0][1][3][2] = 0.0;
    chrisD[0][1][3][3] = 0.0;
    chrisD[0][2][0][0] = 0.0;
    chrisD[0][2][0][1] = 0.0;
    chrisD[0][2][0][2] = 0.0;
    chrisD[0][2][0][3] = 0.0;
    chrisD[0][2][1][0] = 0.0;
    chrisD[0][2][1][1] = 0.0;
    chrisD[0][2][1][2] = 0.0;
    chrisD[0][2][1][3] = 0.0;
    chrisD[0][2][2][0] = 0.0;
    chrisD[0][2][2][1] = 0.0;
    chrisD[0][2][2][2] = 0.0;
    chrisD[0][2][2][3] = 0.0;
    chrisD[0][2][3][0] = 0.0;
    chrisD[0][2][3][1] = 0.0;
    chrisD[0][2][3][2] = 0.0;
    chrisD[0][2][3][3] = 0.0;
    chrisD[0][3][0][0] = 0.0;
    chrisD[0][3][0][1] = 0.0;
    chrisD[0][3][0][2] = 0.0;
    chrisD[0][3][0][3] = 0.0;
    chrisD[0][3][1][0] = 0.0;
    chrisD[0][3][1][1] = 0.0;
    chrisD[0][3][1][2] = 0.0;
    chrisD[0][3][1][3] = 0.0;
    chrisD[0][3][2][0] = 0.0;
    chrisD[0][3][2][1] = 0.0;
    chrisD[0][3][2][2] = 0.0;
    chrisD[0][3][2][3] = 0.0;
    chrisD[0][3][3][0] = 0.0;
    chrisD[0][3][3][1] = 0.0;
    chrisD[0][3][3][2] = 0.0;
    chrisD[0][3][3][3] = 0.0;
    chrisD[1][0][0][0] = 0.0;
    chrisD[1][0][0][1] = -t17;
    chrisD[1][0][0][2] = 0.0;
    chrisD[1][0][0][3] = 0.0;
    chrisD[1][0][1][0] = 0.0;
    chrisD[1][0][1][1] = 0.0;
    chrisD[1][0][1][2] = 0.0;
    chrisD[1][0][1][3] = 0.0;
    chrisD[1][0][2][0] = 0.0;
    chrisD[1][0][2][1] = 0.0;
    chrisD[1][0][2][2] = 0.0;
    chrisD[1][0][2][3] = 0.0;
    chrisD[1][0][3][0] = 0.0;
    chrisD[1][0][3][1] = 0.0;
    chrisD[1][0][3][2] = 0.0;
    chrisD[1][0][3][3] = 0.0;
    chrisD[1][1][0][0] = 0.0;
    chrisD[1][1][0][1] = 0.0;
    chrisD[1][1][0][2] = 0.0;
    chrisD[1][1][0][3] = 0.0;
    chrisD[1][1][1][0] = 0.0;
    chrisD[1][1][1][1] = -t17;
    chrisD[1][1][1][2] = 0.0;
    chrisD[1][1][1][3] = 0.0;
    chrisD[1][1][2][0] = 0.0;
    chrisD[1][1][2][1] = 0.0;
    chrisD[1][1][2][2] = 0.0;
    chrisD[1][1][2][3] = 0.0;
    chrisD[1][1][3][0] = 0.0;
    chrisD[1][1][3][1] = 0.0;
    chrisD[1][1][3][2] = 0.0;
    chrisD[1][1][3][3] = 0.0;
    chrisD[1][2][0][0] = 0.0;
    chrisD[1][2][0][1] = 0.0;
    chrisD[1][2][0][2] = 0.0;
    chrisD[1][2][0][3] = 0.0;
    chrisD[1][2][1][0] = 0.0;
    chrisD[1][2][1][1] = 0.0;
    chrisD[1][2][1][2] = 0.0;
    chrisD[1][2][1][3] = 0.0;
    chrisD[1][2][2][0] = 0.0;
    chrisD[1][2][2][1] = -t20;
    chrisD[1][2][2][2] = 0.0;
    chrisD[1][2][2][3] = 0.0;
    chrisD[1][2][3][0] = 0.0;
    chrisD[1][2][3][1] = 0.0;
    chrisD[1][2][3][2] = 0.0;
    chrisD[1][2][3][3] = 0.0;
    chrisD[1][3][0][0] = 0.0;
    chrisD[1][3][0][1] = 0.0;
    chrisD[1][3][0][2] = 0.0;
    chrisD[1][3][0][3] = 0.0;
    chrisD[1][3][1][0] = 0.0;
    chrisD[1][3][1][1] = 0.0;
    chrisD[1][3][1][2] = 0.0;
    chrisD[1][3][1][3] = 0.0;
    chrisD[1][3][2][0] = 0.0;
    chrisD[1][3][2][1] = 0.0;
    chrisD[1][3][2][2] = 0.0;
    chrisD[1][3][2][3] = 0.0;
    chrisD[1][3][3][0] = 0.0;
    chrisD[1][3][3][1] = -t20;
    chrisD[1][3][3][2] = 0.0;
    chrisD[1][3][3][3] = 0.0;
    chrisD[2][0][0][0] = 0.0;
    chrisD[2][0][0][1] = 0.0;
    chrisD[2][0][0][2] = 0.0;
    chrisD[2][0][0][3] = 0.0;
    chrisD[2][0][1][0] = 0.0;
    chrisD[2][0][1][1] = 0.0;
    chrisD[2][0][1][2] = 0.0;
    chrisD[2][0][1][3] = 0.0;
    chrisD[2][0][2][0] = 0.0;
    chrisD[2][0][2][1] = 0.0;
    chrisD[2][0][2][2] = 0.0;
    chrisD[2][0][2][3] = 0.0;
    chrisD[2][0][3][0] = 0.0;
    chrisD[2][0][3][1] = 0.0;
    chrisD[2][0][3][2] = 0.0;
    chrisD[2][0][3][3] = 0.0;
    chrisD[2][1][0][0] = 0.0;
    chrisD[2][1][0][1] = 0.0;
    chrisD[2][1][0][2] = 0.0;
    chrisD[2][1][0][3] = 0.0;
    chrisD[2][1][1][0] = 0.0;
    chrisD[2][1][1][1] = 0.0;
    chrisD[2][1][1][2] = 0.0;
    chrisD[2][1][1][3] = 0.0;
    chrisD[2][1][2][0] = 0.0;
    chrisD[2][1][2][1] = -t20;
    chrisD[2][1][2][2] = 0.0;
    chrisD[2][1][2][3] = 0.0;
    chrisD[2][1][3][0] = 0.0;
    chrisD[2][1][3][1] = 0.0;
    chrisD[2][1][3][2] = 0.0;
    chrisD[2][1][3][3] = 0.0;
    chrisD[2][2][0][0] = 0.0;
    chrisD[2][2][0][1] = 0.0;
    chrisD[2][2][0][2] = 0.0;
    chrisD[2][2][0][3] = 0.0;
    chrisD[2][2][1][0] = 0.0;
    chrisD[2][2][1][1] = -t5 * t21;
    chrisD[2][2][1][2] = 0.0;
    chrisD[2][2][1][3] = 0.0;
    chrisD[2][2][2][0] = 0.0;
    chrisD[2][2][2][1] = 0.0;
    chrisD[2][2][2][2] = 0.0;
    chrisD[2][2][2][3] = 0.0;
    chrisD[2][2][3][0] = 0.0;
    chrisD[2][2][3][1] = 0.0;
    chrisD[2][2][3][2] = 0.0;
    chrisD[2][2][3][3] = 0.0;
    chrisD[2][3][0][0] = 0.0;
    chrisD[2][3][0][1] = 0.0;
    chrisD[2][3][0][2] = 0.0;
    chrisD[2][3][0][3] = 0.0;
    chrisD[2][3][1][0] = 0.0;
    chrisD[2][3][1][1] = 0.0;
    chrisD[2][3][1][2] = 0.0;
    chrisD[2][3][1][3] = 0.0;
    chrisD[2][3][2][0] = 0.0;
    chrisD[2][3][2][1] = 0.0;
    chrisD[2][3][2][2] = 0.0;
    chrisD[2][3][2][3] = 0.0;
    chrisD[2][3][3][0] = 0.0;
    chrisD[2][3][3][1] = 0.0;
    chrisD[2][3][3][2] = t26;
    chrisD[2][3][3][3] = 0.0;
    chrisD[3][0][0][0] = 0.0;
    chrisD[3][0][0][1] = 0.0;
    chrisD[3][0][0][2] = 0.0;
    chrisD[3][0][0][3] = 0.0;
    chrisD[3][0][1][0] = 0.0;
    chrisD[3][0][1][1] = 0.0;
    chrisD[3][0][1][2] = 0.0;
    chrisD[3][0][1][3] = 0.0;
    chrisD[3][0][2][0] = 0.0;
    chrisD[3][0][2][1] = 0.0;
    chrisD[3][0][2][2] = 0.0;
    chrisD[3][0][2][3] = 0.0;
    chrisD[3][0][3][0] = 0.0;
    chrisD[3][0][3][1] = 0.0;
    chrisD[3][0][3][2] = 0.0;
    chrisD[3][0][3][3] = 0.0;
    chrisD[3][1][0][0] = 0.0;
    chrisD[3][1][0][1] = 0.0;
    chrisD[3][1][0][2] = 0.0;
    chrisD[3][1][0][3] = 0.0;
    chrisD[3][1][1][0] = 0.0;
    chrisD[3][1][1][1] = 0.0;
    chrisD[3][1][1][2] = 0.0;
    chrisD[3][1][1][3] = 0.0;
    chrisD[3][1][2][0] = 0.0;
    chrisD[3][1][2][1] = 0.0;
    chrisD[3][1][2][2] = 0.0;
    chrisD[3][1][2][3] = 0.0;
    chrisD[3][1][3][0] = 0.0;
    chrisD[3][1][3][1] = -t20;
    chrisD[3][1][3][2] = 0.0;
    chrisD[3][1][3][3] = 0.0;
    chrisD[3][2][0][0] = 0.0;
    chrisD[3][2][0][1] = 0.0;
    chrisD[3][2][0][2] = 0.0;
    chrisD[3][2][0][3] = 0.0;
    chrisD[3][2][1][0] = 0.0;
    chrisD[3][2][1][1] = 0.0;
    chrisD[3][2][1][2] = 0.0;
    chrisD[3][2][1][3] = 0.0;
    chrisD[3][2][2][0] = 0.0;
    chrisD[3][2][2][1] = 0.0;
    chrisD[3][2][2][2] = 0.0;
    chrisD[3][2][2][3] = 0.0;
    chrisD[3][2][3][0] = 0.0;
    chrisD[3][2][3][1] = 0.0;
    chrisD[3][2][3][2] = t26;
    chrisD[3][2][3][3] = 0.0;
    chrisD[3][3][0][0] = 0.0;
    chrisD[3][3][0][1] = 0.0;
    chrisD[3][3][0][2] = 0.0;
    chrisD[3][3][0][3] = 0.0;
    chrisD[3][3][1][0] = 0.0;
    chrisD[3][3][1][1] = t5 * t25 * t21;
    chrisD[3][3][1][2] = -2.0 * t6 * rs * t30 * t23;
    chrisD[3][3][1][3] = 0.0;
    chrisD[3][3][2][0] = 0.0;
    chrisD[3][3][2][1] = 0.0;
    chrisD[3][3][2][2] = -2.0 * t24 + 1.0;
    chrisD[3][3][2][3] = 0.0;
    chrisD[3][3][3][0] = 0.0;
    chrisD[3][3][3][1] = 0.0;
    chrisD[3][3][3][2] = 0.0;
    chrisD[3][3][3][3] = 0.0;

    return true;
}

/*! Transform local 4-direction to coordinate 4-direction.
 *
 *  \param  pos  :  pointer to position array.
 *  \param  ldir :  pointer to local direction array.
 *  \param  dir  :  pointer to calculated coordinate direction array.
 *  \param  type :  type of local tetrad.
 */
void MetricSchwarzschildTortoise::localToCoord(const double* pos, const double* ldir, double* dir,
        enum_nat_tetrad_type) {
    double rho   = pos[1];
    double theta = pos[2];

    double r = calc_r(rho);
    double w = sqrt(1.0 - rs / r);

    dir[0] = ldir[0] / w / mSpeedOfLight;
    dir[1] = ldir[1] / w;
    dir[2] = ldir[2] / r;
    dir[3] = ldir[3] / (r * sin(theta));
}

/*! Transform coordinate 4-direction to local 4-direction.
 *
 *  \param  pos  :  pointer to position array.
 *  \param  cdir :  pointer to coordinate direction.
 *  \param  ldir :  pointer to calculated local direction array.
 *  \param  type :  type of local tetrad.
 */
void MetricSchwarzschildTortoise::coordToLocal(const double* pos, const double* cdir, double* ldir,
        enum_nat_tetrad_type) {
    double rho   = pos[1];
    double theta = pos[2];

    double r = calc_r(rho);
    double w = sqrt(1.0 - rs / r);

    ldir[0] = cdir[0] * w * mSpeedOfLight;
    ldir[1] = cdir[1] * w;
    ldir[2] = cdir[2] * r;
    ldir[3] = cdir[3] * r * sin(theta);
}


/*! Test break condition.
 *
 *  \param pos    : pointer to position array.
 *  \return true  : radial position r < 0.0 or  r^2<=(1.0+eps)*rs^2.
 *  \return false : position is valid.
 */
bool MetricSchwarzschildTortoise::breakCondition(const double*) {
    bool br = false;
    return br;
}


/* Calculate right hand side of the geodesic equation in first order form.
 *
 *  \param  y[]   : pointer to position and direction coordinates.
 *  \param  dydx[] : pointer to right side of geodesic equation.
 *
bool MetricSchwarzschildTortoise::calcDerivs ( const double y[], double dydx[] )
{
  dydx[0] = y[4];
  dydx[1] = y[5];
  dydx[2] = y[6];
  dydx[3] = y[7];

  double r     = y[1];
  double theta = y[2];

  dydx[4] = -rs/(r*(r-rs))*y[4]*y[5];
  dydx[5] = -0.5*mSpeedOfLight*mSpeedOfLight*rs*(r-rs)/pow(r,3.0)*y[4]*y[4] + 0.5*rs/(r*(r-rs))*y[5]*y[5] + (r-rs)*(y[6]*y[6]+ sin(theta)*sin(theta)*y[7]*y[7]);
  dydx[6] = -2.0/r*y[5]*y[6] + sin(theta)*cos(theta)*y[7]*y[7];
  dydx[7] = -2.0/r*y[5]*y[7] - 2.0*cos(theta)/sin(theta)*y[6]*y[7];

  return true;
}
*/

/*! Tests whether the constraint equation is fulfilled.
 *
 *  The constraint equation for lightlike and timelike geodesics reads:
 \verbatim
     sum = g_{\mu\nu} dot(x)^{\mu} dot(x)^{\nu} - kappa c^2 = 0.
 \endverbatim
 *  However, take care of the limited double precision.
 *  \param  y[]   : pointer to position and direction coordinates.
 *  \param  kappa : timelike (-1.0), lightlike (0.0).
 *  \return double : sum.
 */
double MetricSchwarzschildTortoise::testConstraint(const double y[], const double kappa) {
    double rho   = y[1];
    double theta = y[2];
    double cm = 1.0 / mSpeedOfLight;

    double r = calc_r(rho);
    // Scale the directions with the speed of light before doubling them !!
    double dt   = y[4];
    double drho = y[5] * cm;
    double dth  = y[6] * cm;
    double dph  = y[7] * cm;

    double sum = -kappa;
    sum += -(1.0 - rs / r) * dt * dt + (1.0 - rs / r) * drho * drho + r * r * (dth * dth + sin(theta) * sin(theta) * dph * dph);
    return sum;
}

/*! Calculate the scalar product between u and v:  g_{ab}u^a v^b.
 *
 *  \param pos  :  pointer to position.
 *  \param u    :  pointer to vector.
 *  \param v    :  pointer to vector.
 *  \param prod :  reference to scalar product.
 *  \param preCalcMetric : calculate metric coefficients before evaluating the scalar product.
 *  \return true  : success.
 *  \return false : position is not valid.
 */
bool MetricSchwarzschildTortoise::calcProduct(const double* pos, const double* u, const double* v, double &prod, bool) {
    prod = 0.0;
    if (breakCondition(pos)) {
        return false;
    }

    double rho   = pos[1];
    double theta = pos[2];

    double r = calc_r(rho);
    prod = -mSpeedOfLight * mSpeedOfLight * (1.0 - rs / r) * u[0] * v[0] + (1.0 - rs / r) * u[1] * v[1] + r * r * (u[2] * v[2] + sin(theta) * sin(theta) * u[3] * v[3]);
    return true;
}

/*! Set parameter 'pName' to 'val'.
 *
 *  Set 'mass' parameter and adjust Schwarzschild radius  rs=2GM/c^2.
 */
bool MetricSchwarzschildTortoise::setParam(const char* pName, double val) {
    if (Metric::setParam(pName, val)) {
        mMass = val;
        rs = 2.0 * mGravConstant * mMass / (mSpeedOfLight * mSpeedOfLight);
    }
    return true;
}

/*!  Transform point p to pseudo cartesian coordinates.
 *
 *  \param p  : point to be transformed.
 *  \param cp : reference to customized point.
 *  \return ID of chart
 */
int MetricSchwarzschildTortoise::transToPseudoCart(vec4 p, vec4 &cp) {
    double rho   = p[1];
    double theta = p[2];
    double phi   = p[3];
    double r = calc_r(rho);

    cp[0] = p[0];
    cp[1] = r * sin(theta) * cos(phi);
    cp[2] = r * sin(theta) * sin(phi);
    cp[3] = r * cos(theta);
    return 0;
}


/*! Transform point p to embedding coordinates.
 *
 *  \param p  : point to be transformed.
 *  \param ep : reference to 'embedded' point.
 *  \return true : success.
 *  \return false : otherwise.
 */
bool MetricSchwarzschildTortoise::transToEmbedding(vec4 p, vec4 &ep) {
    vec4 cp;
    transToPseudoCart(p, cp);

    double rho = p[1];
    double r = calc_r(rho);

    double x = cp[1];
    double y = cp[2];
    double z;
    if (r >= rs) {
        z = 2.0 * sqrt(rs) * sqrt(r - rs);
        ep = vec4(p[0], x, y, z);
        return true;
    }
    return false;
}


/*! Set embedding parameters.
 *
 *  \param  name : embedding parameter name.
 *  \param  val  : embedding parameter value.
 *  \return true  : success.
 *  \return false : parameter not valid.
 */
bool MetricSchwarzschildTortoise::setEmbeddingParam(const char* name, double val) {
    Metric::setEmbeddingParam(name, val);

    if (strcmp(name,"emb_rmin") == 0) {
        mEmb_rmin = val;
    }
    else if (strcmp(name,"emb_rmax") == 0) {
        mEmb_rmax = val;
    }
    else if (strcmp(name,"emb_r_num") == 0) {
        mEmb_r_num = val;
    }
    else if (strcmp(name,"emb_phi_num") == 0) {
        mEmb_phi_num = val;
    }
    return testEmbeddingParams();
}

/*! Test embedding parameters
 *  \return  true : all parameters are ok
 *  \return  false : at least one parameter had to be adjusted.
 */
bool MetricSchwarzschildTortoise::testEmbeddingParams() {
    bool allOk = true;
    if (mEmb_rmin < rs) {
        mEmb_rmin = rs;
        allOk &= false;
    }
    if (mEmb_rmax < rs) {
        mEmb_rmax = rs;
        allOk &= false;
    }
    if (mEmb_r_num < 2.0) {
        mEmb_r_num = 2.0;
        allOk &= false;
    }

    if (mEmb_phi_num < 4) {
        mEmb_phi_num = 4;
        allOk &= false;
    }
    return allOk;
}

/*! Generate vertices for the embedding diagram.
 *
 *  \param verts : reference to vector of vertices.
 *  \param indices : reference to vector of indices.
 *  \param numElems : number of elements in a strip.
 *  \param counter  : number of strips.
 *  \return int : number of vertices.
 */
//int MetricSchwarzschildTortoise::getEmbeddingVertices(std::vector<vec3> &verts,
//        std::vector<int> &indices, unsigned int &numElems, unsigned int &counter) {
//    if (!verts.empty()) {
//        verts.clear();
//    }

//    if (!indices.empty()) {
//        indices.clear();
//    }

//    testEmbeddingParams();
//    mEmb_rstep = (mEmb_rmax - mEmb_rmin) / mEmb_r_num;
//    mEmb_phistep = 2.0 * M_PI / mEmb_phi_num;

//    numElems = int(mEmb_r_num);
//    counter  = int(mEmb_phi_num) + 1;

//    int vnum;

//    double x, y, z, r, phi;
//    for (unsigned int k = 0; k < counter; k++) {
//        phi = k * mEmb_phistep;
//        for (unsigned int j = 0; j < numElems; j++) {
//            r = mEmb_rmin + j * mEmb_rstep;
//            x = r * cos(phi);
//            y = r * sin(phi);
//            if (r >= rs) {
//                z = 2.0 * sqrt(rs) * sqrt(fabs(r - rs));
//                verts.push_back(vec3(x, y, z));

//                vnum = k * numElems + j;

//                indices.push_back(vnum);
//                indices.push_back(vnum + numElems);
//            }
//        }
//    }

//    int numVerts = (int)verts.size();
//    int numInds  = (int)indices.size();

//    if (2 * numVerts == numInds) {
//        return numVerts;
//    }

//    return 0;
//}



/*! Generate report.
 */
bool MetricSchwarzschildTortoise::report(const vec4 , const vec4 , std::string &text) {
    std::stringstream ss;
    ss << "Report for Schwarzschild metric\n\t tortoise coordinates : (t,rho,theta,phi)\n";
    ss << "---------------------------------------------------------------\n";
    ss << "  physical units ................................. yes\n";
    ss.precision(DEF_FIXED_REPORT_PRECISION);
    ss.setf(std::ios::fixed);
    ss << "  Schwarzschild radius ...................... r_s = 2GM/c^2 = " << rs << std::endl;
    ss << "  Photon orbit .................... rho_ph = (3/2-ln(2))*rs = " << (1.5 - log(2.0))*rs << std::endl;
    ss << "  innermost stable circular orbit  rho_isco = (3+ln(2))*r_s = " << (3.0 + log(2.0))*rs << std::endl;
    ss << "\n  Note that pseudo-cartesian coordinates use the Schwarzschild radial coordinate\n  r = rs*{1+W[exp(rho/rs-1)]}." << std::endl;
    text = ss.str();
    return true;
}



// ********************************* protected methods *****************************
/*!
 */
void MetricSchwarzschildTortoise::setStandardValues() {
    mInitPos[0] = 0.0;
    mInitPos[1] = 1.0;
    mInitPos[2] = M_PI_2;
    mInitPos[3] = 0.0;
    mInitDir[0] = 1.0;
    mInitDir[1] = 0.0;
    mInitDir[2] = 0.0;

    mCoordNames[0] = std::string("t");
    mCoordNames[1] = std::string("rho");
    mCoordNames[2] = std::string("theta");
    mCoordNames[3] = std::string("phi");
}

/*! Calculate r from rho.
 * \param rho
 */
double MetricSchwarzschildTortoise::calc_r(const double rho) {
    return rs * (1.0 + gsl_sf_lambert_W0(exp(rho / rs - 1.0)));
}

} // end namespace m4d
