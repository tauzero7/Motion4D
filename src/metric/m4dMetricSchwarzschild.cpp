// -------------------------------------------------------------------------------
/*
   m4dMetricSchwarzschild.cpp

  Copyright (c) 2009-2014  Thomas Mueller, Frank Grave


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

#include "m4dMetricSchwarzschild.h"
#include <cmath>

#define sign(x) ( (x>=0) ? 1.0 : -1.0 )

namespace m4d {

#define eps 1.0e-6


/*! Standard constructor for the Schwarzschild metric.
 *
 * \param  mass : mass of the black hole.
 */
MetricSchwarzschild::MetricSchwarzschild(double mass) {
    mMetricName  = "Schwarzschild";
    mMetricCPPfilename = "m4dMetricSchwarzschild.cpp";
    setCoordType(enum_coordinate_spherical);

    mPhysicalUnits = enum_physical_constants_geom;
    mSpeedOfLight = 1.0;
    mGravConstant = 1.0;

    addParam("mass", mass);
    mMass = mass;
    rs = 2.0 * mGravConstant * mMass / (mSpeedOfLight * mSpeedOfLight);

    /*  Only a static tetrad is defined  */
    mLocTeds.push_back(enum_nat_tetrad_static);


    mDrawTypes.push_back(enum_draw_embedding);
    mDrawTypes.push_back(enum_draw_twoplusone);
    mDrawTypes.push_back(enum_draw_effpoti);
    mDrawTypes.push_back(enum_draw_custom);


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
}

/*!
 */
MetricSchwarzschild::~MetricSchwarzschild() {
}


// *********************************** public methods ******************************

/*! Calculate the contravariant metric components at position 'pos'.
 *
 *  \param pos : pointer to position.
 */
bool MetricSchwarzschild::calculateMetric(const double* pos) {
    double r     = pos[1];
    double theta = pos[2];

    double c = mSpeedOfLight;
    double c2 = c * c;

    double t1 = c2;
    double t3 = 1.0 / r;
    double t9 = r * r;
    double t10 = sin(theta);
    double t11 = t10 * t10;

    g_compts[0][0] = -t1 + t1 * rs * t3;
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

/*! Calculate the Christoffel symbols of the second kind at position 'pos'.
 *
 *  \param pos : pointer to position.
 */
bool MetricSchwarzschild::calculateChristoffels(const double* pos) {
    double r     = pos[1];
    double theta = pos[2];

    double c = mSpeedOfLight;

    double t1 = r - rs;
    double t2 = r * r;
    double t6 = c * c;
    double t10 = 1.0 / r;
    //  double t14 = t10/t1*rs/2.0;
    double t14 = t10 / t1 * rs * 0.5;
    double t15 = sin(theta);
    double t17 = cos(theta);
    double t18 = 1.0 / t15 * t17;
    double t19 = t15 * t15;

    christoffel[0][0][0] = 0.0;
    christoffel[0][0][1] = t1 / t2 / r * t6 * rs / 2.0;
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

/*! Calculate Jacobi matrix.
 *
 *  \param pos : pointer to position.
 */
bool MetricSchwarzschild::calculateChrisD(const double* pos) {
    double r     = pos[1];
    double theta = pos[2];

    double c = mSpeedOfLight;

    double t1 = c * c;
    double t3 = 2.0 * r;
    double t6 = r * r;
    double t7 = t6 * t6;
    double t14 = 1 / t6;
    double t15 = r - rs;
    double t16 = t15 * t15;
    double t20 = rs * (t3 - rs) * t14 / t16 / 2.0;
    double t21 = cos(theta);
    double t22 = t21 * t21;
    double t23 = sin(theta);
    double t24 = t23 * t23;
    double t27 = (t22 + t24) / t24;

    chrisD[0][0][0][0] = 0.0;
    chrisD[0][0][0][1] = 0.0;
    chrisD[0][0][0][2] = 0.0;
    chrisD[0][0][0][3] = 0.0;
    chrisD[0][0][1][0] = 0.0;
    chrisD[0][0][1][1] = -t1 * rs * (t3 - 3.0 * rs) / t7 / 2.0;
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
    chrisD[0][1][0][1] = -t20;
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
    chrisD[1][0][0][1] = -t20;
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
    chrisD[1][1][1][1] = t20;
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
    chrisD[1][2][2][1] = -t14;
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
    chrisD[1][3][3][1] = -t14;
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
    chrisD[2][1][2][1] = -t14;
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
    chrisD[2][2][1][1] = -1.0;
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
    chrisD[2][3][3][2] = -t27;
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
    chrisD[3][1][3][1] = -t14;
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
    chrisD[3][2][3][2] = -t27;
    chrisD[3][2][3][3] = 0.0;
    chrisD[3][3][0][0] = 0.0;
    chrisD[3][3][0][1] = 0.0;
    chrisD[3][3][0][2] = 0.0;
    chrisD[3][3][0][3] = 0.0;
    chrisD[3][3][1][0] = 0.0;
    chrisD[3][3][1][1] = -t24;
    chrisD[3][3][1][2] = -2.0 * t15 * t23 * t21;
    chrisD[3][3][1][3] = 0.0;
    chrisD[3][3][2][0] = 0.0;
    chrisD[3][3][2][1] = 0.0;
    chrisD[3][3][2][2] = -t22 + t24;
    chrisD[3][3][2][3] = 0.0;
    chrisD[3][3][3][0] = 0.0;
    chrisD[3][3][3][1] = 0.0;
    chrisD[3][3][3][2] = 0.0;
    chrisD[3][3][3][3] = 0.0;

    return true;
}


/*! Calculate the Riemann tensor R^a_bcd at position 'pos'.
 *
 *  \param  pos  :  pointer to coordinate position where the Weyl tensor have to be evaluated.
 *  \return true :  successfull
 */
bool MetricSchwarzschild::calculateRiemann(const double* pos) {
    double r     = pos[1];
    double theta = pos[2];

    double c = mSpeedOfLight;

    double t1 = r * r;
    double t4 = r - rs;
    double t6 = 1 / t1 * rs / t4;
    double t7 = 1 / r;
    double t8 = rs * t7;
    double t9 = t8 / 2.0;
    double t10 = sin(theta);
    double t11 = t10 * t10;
    double t13 = t7 * t11 * rs;
    double t14 = t13 / 2.0;
    double t15 = c * c;
    double t17 = t1 * t1;
    double t20 = t4 * t15 * rs / t17;
    double t21 = t20 / 2.0;
    double t22 = t6 / 2.0;

    riem[0][0][0][0] = 0.0;
    riem[0][0][0][1] = 0.0;
    riem[0][0][0][2] = 0.0;
    riem[0][0][0][3] = 0.0;
    riem[0][0][1][0] = 0.0;
    riem[0][0][1][1] = 0.0;
    riem[0][0][1][2] = 0.0;
    riem[0][0][1][3] = 0.0;
    riem[0][0][2][0] = 0.0;
    riem[0][0][2][1] = 0.0;
    riem[0][0][2][2] = 0.0;
    riem[0][0][2][3] = 0.0;
    riem[0][0][3][0] = 0.0;
    riem[0][0][3][1] = 0.0;
    riem[0][0][3][2] = 0.0;
    riem[0][0][3][3] = 0.0;
    riem[0][1][0][0] = 0.0;
    riem[0][1][0][1] = t6;
    riem[0][1][0][2] = 0.0;
    riem[0][1][0][3] = 0.0;
    riem[0][1][1][0] = -t6;
    riem[0][1][1][1] = 0.0;
    riem[0][1][1][2] = 0.0;
    riem[0][1][1][3] = 0.0;
    riem[0][1][2][0] = 0.0;
    riem[0][1][2][1] = 0.0;
    riem[0][1][2][2] = 0.0;
    riem[0][1][2][3] = 0.0;
    riem[0][1][3][0] = 0.0;
    riem[0][1][3][1] = 0.0;
    riem[0][1][3][2] = 0.0;
    riem[0][1][3][3] = 0.0;
    riem[0][2][0][0] = 0.0;
    riem[0][2][0][1] = 0.0;
    riem[0][2][0][2] = -t9;
    riem[0][2][0][3] = 0.0;
    riem[0][2][1][0] = 0.0;
    riem[0][2][1][1] = 0.0;
    riem[0][2][1][2] = 0.0;
    riem[0][2][1][3] = 0.0;
    riem[0][2][2][0] = t9;
    riem[0][2][2][1] = 0.0;
    riem[0][2][2][2] = 0.0;
    riem[0][2][2][3] = 0.0;
    riem[0][2][3][0] = 0.0;
    riem[0][2][3][1] = 0.0;
    riem[0][2][3][2] = 0.0;
    riem[0][2][3][3] = 0.0;
    riem[0][3][0][0] = 0.0;
    riem[0][3][0][1] = 0.0;
    riem[0][3][0][2] = 0.0;
    riem[0][3][0][3] = -t14;
    riem[0][3][1][0] = 0.0;
    riem[0][3][1][1] = 0.0;
    riem[0][3][1][2] = 0.0;
    riem[0][3][1][3] = 0.0;
    riem[0][3][2][0] = 0.0;
    riem[0][3][2][1] = 0.0;
    riem[0][3][2][2] = 0.0;
    riem[0][3][2][3] = 0.0;
    riem[0][3][3][0] = t14;
    riem[0][3][3][1] = 0.0;
    riem[0][3][3][2] = 0.0;
    riem[0][3][3][3] = 0.0;
    riem[1][0][0][0] = 0.0;
    riem[1][0][0][1] = t20;
    riem[1][0][0][2] = 0.0;
    riem[1][0][0][3] = 0.0;
    riem[1][0][1][0] = -t20;
    riem[1][0][1][1] = 0.0;
    riem[1][0][1][2] = 0.0;
    riem[1][0][1][3] = 0.0;
    riem[1][0][2][0] = 0.0;
    riem[1][0][2][1] = 0.0;
    riem[1][0][2][2] = 0.0;
    riem[1][0][2][3] = 0.0;
    riem[1][0][3][0] = 0.0;
    riem[1][0][3][1] = 0.0;
    riem[1][0][3][2] = 0.0;
    riem[1][0][3][3] = 0.0;
    riem[1][1][0][0] = 0.0;
    riem[1][1][0][1] = 0.0;
    riem[1][1][0][2] = 0.0;
    riem[1][1][0][3] = 0.0;
    riem[1][1][1][0] = 0.0;
    riem[1][1][1][1] = 0.0;
    riem[1][1][1][2] = 0.0;
    riem[1][1][1][3] = 0.0;
    riem[1][1][2][0] = 0.0;
    riem[1][1][2][1] = 0.0;
    riem[1][1][2][2] = 0.0;
    riem[1][1][2][3] = 0.0;
    riem[1][1][3][0] = 0.0;
    riem[1][1][3][1] = 0.0;
    riem[1][1][3][2] = 0.0;
    riem[1][1][3][3] = 0.0;
    riem[1][2][0][0] = 0.0;
    riem[1][2][0][1] = 0.0;
    riem[1][2][0][2] = 0.0;
    riem[1][2][0][3] = 0.0;
    riem[1][2][1][0] = 0.0;
    riem[1][2][1][1] = 0.0;
    riem[1][2][1][2] = -t9;
    riem[1][2][1][3] = 0.0;
    riem[1][2][2][0] = 0.0;
    riem[1][2][2][1] = t9;
    riem[1][2][2][2] = 0.0;
    riem[1][2][2][3] = 0.0;
    riem[1][2][3][0] = 0.0;
    riem[1][2][3][1] = 0.0;
    riem[1][2][3][2] = 0.0;
    riem[1][2][3][3] = 0.0;
    riem[1][3][0][0] = 0.0;
    riem[1][3][0][1] = 0.0;
    riem[1][3][0][2] = 0.0;
    riem[1][3][0][3] = 0.0;
    riem[1][3][1][0] = 0.0;
    riem[1][3][1][1] = 0.0;
    riem[1][3][1][2] = 0.0;
    riem[1][3][1][3] = -t14;
    riem[1][3][2][0] = 0.0;
    riem[1][3][2][1] = 0.0;
    riem[1][3][2][2] = 0.0;
    riem[1][3][2][3] = 0.0;
    riem[1][3][3][0] = 0.0;
    riem[1][3][3][1] = t14;
    riem[1][3][3][2] = 0.0;
    riem[1][3][3][3] = 0.0;
    riem[2][0][0][0] = 0.0;
    riem[2][0][0][1] = 0.0;
    riem[2][0][0][2] = -t21;
    riem[2][0][0][3] = 0.0;
    riem[2][0][1][0] = 0.0;
    riem[2][0][1][1] = 0.0;
    riem[2][0][1][2] = 0.0;
    riem[2][0][1][3] = 0.0;
    riem[2][0][2][0] = t21;
    riem[2][0][2][1] = 0.0;
    riem[2][0][2][2] = 0.0;
    riem[2][0][2][3] = 0.0;
    riem[2][0][3][0] = 0.0;
    riem[2][0][3][1] = 0.0;
    riem[2][0][3][2] = 0.0;
    riem[2][0][3][3] = 0.0;
    riem[2][1][0][0] = 0.0;
    riem[2][1][0][1] = 0.0;
    riem[2][1][0][2] = 0.0;
    riem[2][1][0][3] = 0.0;
    riem[2][1][1][0] = 0.0;
    riem[2][1][1][1] = 0.0;
    riem[2][1][1][2] = t22;
    riem[2][1][1][3] = 0.0;
    riem[2][1][2][0] = 0.0;
    riem[2][1][2][1] = -t22;
    riem[2][1][2][2] = 0.0;
    riem[2][1][2][3] = 0.0;
    riem[2][1][3][0] = 0.0;
    riem[2][1][3][1] = 0.0;
    riem[2][1][3][2] = 0.0;
    riem[2][1][3][3] = 0.0;
    riem[2][2][0][0] = 0.0;
    riem[2][2][0][1] = 0.0;
    riem[2][2][0][2] = 0.0;
    riem[2][2][0][3] = 0.0;
    riem[2][2][1][0] = 0.0;
    riem[2][2][1][1] = 0.0;
    riem[2][2][1][2] = 0.0;
    riem[2][2][1][3] = 0.0;
    riem[2][2][2][0] = 0.0;
    riem[2][2][2][1] = 0.0;
    riem[2][2][2][2] = 0.0;
    riem[2][2][2][3] = 0.0;
    riem[2][2][3][0] = 0.0;
    riem[2][2][3][1] = 0.0;
    riem[2][2][3][2] = 0.0;
    riem[2][2][3][3] = 0.0;
    riem[2][3][0][0] = 0.0;
    riem[2][3][0][1] = 0.0;
    riem[2][3][0][2] = 0.0;
    riem[2][3][0][3] = 0.0;
    riem[2][3][1][0] = 0.0;
    riem[2][3][1][1] = 0.0;
    riem[2][3][1][2] = 0.0;
    riem[2][3][1][3] = 0.0;
    riem[2][3][2][0] = 0.0;
    riem[2][3][2][1] = 0.0;
    riem[2][3][2][2] = 0.0;
    riem[2][3][2][3] = t13;
    riem[2][3][3][0] = 0.0;
    riem[2][3][3][1] = 0.0;
    riem[2][3][3][2] = -t13;
    riem[2][3][3][3] = 0.0;
    riem[3][0][0][0] = 0.0;
    riem[3][0][0][1] = 0.0;
    riem[3][0][0][2] = 0.0;
    riem[3][0][0][3] = -t21;
    riem[3][0][1][0] = 0.0;
    riem[3][0][1][1] = 0.0;
    riem[3][0][1][2] = 0.0;
    riem[3][0][1][3] = 0.0;
    riem[3][0][2][0] = 0.0;
    riem[3][0][2][1] = 0.0;
    riem[3][0][2][2] = 0.0;
    riem[3][0][2][3] = 0.0;
    riem[3][0][3][0] = t21;
    riem[3][0][3][1] = 0.0;
    riem[3][0][3][2] = 0.0;
    riem[3][0][3][3] = 0.0;
    riem[3][1][0][0] = 0.0;
    riem[3][1][0][1] = 0.0;
    riem[3][1][0][2] = 0.0;
    riem[3][1][0][3] = 0.0;
    riem[3][1][1][0] = 0.0;
    riem[3][1][1][1] = 0.0;
    riem[3][1][1][2] = 0.0;
    riem[3][1][1][3] = t22;
    riem[3][1][2][0] = 0.0;
    riem[3][1][2][1] = 0.0;
    riem[3][1][2][2] = 0.0;
    riem[3][1][2][3] = 0.0;
    riem[3][1][3][0] = 0.0;
    riem[3][1][3][1] = -t22;
    riem[3][1][3][2] = 0.0;
    riem[3][1][3][3] = 0.0;
    riem[3][2][0][0] = 0.0;
    riem[3][2][0][1] = 0.0;
    riem[3][2][0][2] = 0.0;
    riem[3][2][0][3] = 0.0;
    riem[3][2][1][0] = 0.0;
    riem[3][2][1][1] = 0.0;
    riem[3][2][1][2] = 0.0;
    riem[3][2][1][3] = 0.0;
    riem[3][2][2][0] = 0.0;
    riem[3][2][2][1] = 0.0;
    riem[3][2][2][2] = 0.0;
    riem[3][2][2][3] = -t8;
    riem[3][2][3][0] = 0.0;
    riem[3][2][3][1] = 0.0;
    riem[3][2][3][2] = t8;
    riem[3][2][3][3] = 0.0;
    riem[3][3][0][0] = 0.0;
    riem[3][3][0][1] = 0.0;
    riem[3][3][0][2] = 0.0;
    riem[3][3][0][3] = 0.0;
    riem[3][3][1][0] = 0.0;
    riem[3][3][1][1] = 0.0;
    riem[3][3][1][2] = 0.0;
    riem[3][3][1][3] = 0.0;
    riem[3][3][2][0] = 0.0;
    riem[3][3][2][1] = 0.0;
    riem[3][3][2][2] = 0.0;
    riem[3][3][2][3] = 0.0;
    riem[3][3][3][0] = 0.0;
    riem[3][3][3][1] = 0.0;
    riem[3][3][3][2] = 0.0;
    riem[3][3][3][3] = 0.0;

    return true;
}

/*! Calculate the Weyl tensor at position 'pos'.
 *
 *  \param  pos  :  pointer to coordinate position where the Weyl tensor have to be evaluated.
 *  \return true :  successfull
 */
bool MetricSchwarzschild::calculateWeyl(const double* pos) {
    double r     = pos[1];
    double theta = pos[2];

    double c = mSpeedOfLight;

    double t1 = c * c;
    double t2 = t1 * rs;
    double t3 = r * r;
    double t6 = t2 / t3 / r;
    double t7 = r - rs;
    double t9 = t7 / t3;
    double t11 = t9 * t2 / 2.0;
    double t12 = sin(theta);
    double t13 = t12 * t12;
    double t16 = t9 * t2 * t13 / 2.0;
    double t18 = 1 / t7 * rs;
    double t19 = t18 / 2.0;
    double t21 = t18 * t13 / 2.0;
    double t23 = r * t13 * rs;

    weyl[0][0][0][0] = 0.0;
    weyl[0][0][0][1] = 0.0;
    weyl[0][0][0][2] = 0.0;
    weyl[0][0][0][3] = 0.0;
    weyl[0][0][1][0] = 0.0;
    weyl[0][0][1][1] = 0.0;
    weyl[0][0][1][2] = 0.0;
    weyl[0][0][1][3] = 0.0;
    weyl[0][0][2][0] = 0.0;
    weyl[0][0][2][1] = 0.0;
    weyl[0][0][2][2] = 0.0;
    weyl[0][0][2][3] = 0.0;
    weyl[0][0][3][0] = 0.0;
    weyl[0][0][3][1] = 0.0;
    weyl[0][0][3][2] = 0.0;
    weyl[0][0][3][3] = 0.0;
    weyl[0][1][0][0] = 0.0;
    weyl[0][1][0][1] = -t6;
    weyl[0][1][0][2] = 0.0;
    weyl[0][1][0][3] = 0.0;
    weyl[0][1][1][0] = t6;
    weyl[0][1][1][1] = 0.0;
    weyl[0][1][1][2] = 0.0;
    weyl[0][1][1][3] = 0.0;
    weyl[0][1][2][0] = 0.0;
    weyl[0][1][2][1] = 0.0;
    weyl[0][1][2][2] = 0.0;
    weyl[0][1][2][3] = 0.0;
    weyl[0][1][3][0] = 0.0;
    weyl[0][1][3][1] = 0.0;
    weyl[0][1][3][2] = 0.0;
    weyl[0][1][3][3] = 0.0;
    weyl[0][2][0][0] = 0.0;
    weyl[0][2][0][1] = 0.0;
    weyl[0][2][0][2] = t11;
    weyl[0][2][0][3] = 0.0;
    weyl[0][2][1][0] = 0.0;
    weyl[0][2][1][1] = 0.0;
    weyl[0][2][1][2] = 0.0;
    weyl[0][2][1][3] = 0.0;
    weyl[0][2][2][0] = -t11;
    weyl[0][2][2][1] = 0.0;
    weyl[0][2][2][2] = 0.0;
    weyl[0][2][2][3] = 0.0;
    weyl[0][2][3][0] = 0.0;
    weyl[0][2][3][1] = 0.0;
    weyl[0][2][3][2] = 0.0;
    weyl[0][2][3][3] = 0.0;
    weyl[0][3][0][0] = 0.0;
    weyl[0][3][0][1] = 0.0;
    weyl[0][3][0][2] = 0.0;
    weyl[0][3][0][3] = t16;
    weyl[0][3][1][0] = 0.0;
    weyl[0][3][1][1] = 0.0;
    weyl[0][3][1][2] = 0.0;
    weyl[0][3][1][3] = 0.0;
    weyl[0][3][2][0] = 0.0;
    weyl[0][3][2][1] = 0.0;
    weyl[0][3][2][2] = 0.0;
    weyl[0][3][2][3] = 0.0;
    weyl[0][3][3][0] = -t16;
    weyl[0][3][3][1] = 0.0;
    weyl[0][3][3][2] = 0.0;
    weyl[0][3][3][3] = 0.0;
    weyl[1][0][0][0] = 0.0;
    weyl[1][0][0][1] = t6;
    weyl[1][0][0][2] = 0.0;
    weyl[1][0][0][3] = 0.0;
    weyl[1][0][1][0] = -t6;
    weyl[1][0][1][1] = 0.0;
    weyl[1][0][1][2] = 0.0;
    weyl[1][0][1][3] = 0.0;
    weyl[1][0][2][0] = 0.0;
    weyl[1][0][2][1] = 0.0;
    weyl[1][0][2][2] = 0.0;
    weyl[1][0][2][3] = 0.0;
    weyl[1][0][3][0] = 0.0;
    weyl[1][0][3][1] = 0.0;
    weyl[1][0][3][2] = 0.0;
    weyl[1][0][3][3] = 0.0;
    weyl[1][1][0][0] = 0.0;
    weyl[1][1][0][1] = 0.0;
    weyl[1][1][0][2] = 0.0;
    weyl[1][1][0][3] = 0.0;
    weyl[1][1][1][0] = 0.0;
    weyl[1][1][1][1] = 0.0;
    weyl[1][1][1][2] = 0.0;
    weyl[1][1][1][3] = 0.0;
    weyl[1][1][2][0] = 0.0;
    weyl[1][1][2][1] = 0.0;
    weyl[1][1][2][2] = 0.0;
    weyl[1][1][2][3] = 0.0;
    weyl[1][1][3][0] = 0.0;
    weyl[1][1][3][1] = 0.0;
    weyl[1][1][3][2] = 0.0;
    weyl[1][1][3][3] = 0.0;
    weyl[1][2][0][0] = 0.0;
    weyl[1][2][0][1] = 0.0;
    weyl[1][2][0][2] = 0.0;
    weyl[1][2][0][3] = 0.0;
    weyl[1][2][1][0] = 0.0;
    weyl[1][2][1][1] = 0.0;
    weyl[1][2][1][2] = -t19;
    weyl[1][2][1][3] = 0.0;
    weyl[1][2][2][0] = 0.0;
    weyl[1][2][2][1] = t19;
    weyl[1][2][2][2] = 0.0;
    weyl[1][2][2][3] = 0.0;
    weyl[1][2][3][0] = 0.0;
    weyl[1][2][3][1] = 0.0;
    weyl[1][2][3][2] = 0.0;
    weyl[1][2][3][3] = 0.0;
    weyl[1][3][0][0] = 0.0;
    weyl[1][3][0][1] = 0.0;
    weyl[1][3][0][2] = 0.0;
    weyl[1][3][0][3] = 0.0;
    weyl[1][3][1][0] = 0.0;
    weyl[1][3][1][1] = 0.0;
    weyl[1][3][1][2] = 0.0;
    weyl[1][3][1][3] = -t21;
    weyl[1][3][2][0] = 0.0;
    weyl[1][3][2][1] = 0.0;
    weyl[1][3][2][2] = 0.0;
    weyl[1][3][2][3] = 0.0;
    weyl[1][3][3][0] = 0.0;
    weyl[1][3][3][1] = t21;
    weyl[1][3][3][2] = 0.0;
    weyl[1][3][3][3] = 0.0;
    weyl[2][0][0][0] = 0.0;
    weyl[2][0][0][1] = 0.0;
    weyl[2][0][0][2] = -t11;
    weyl[2][0][0][3] = 0.0;
    weyl[2][0][1][0] = 0.0;
    weyl[2][0][1][1] = 0.0;
    weyl[2][0][1][2] = 0.0;
    weyl[2][0][1][3] = 0.0;
    weyl[2][0][2][0] = t11;
    weyl[2][0][2][1] = 0.0;
    weyl[2][0][2][2] = 0.0;
    weyl[2][0][2][3] = 0.0;
    weyl[2][0][3][0] = 0.0;
    weyl[2][0][3][1] = 0.0;
    weyl[2][0][3][2] = 0.0;
    weyl[2][0][3][3] = 0.0;
    weyl[2][1][0][0] = 0.0;
    weyl[2][1][0][1] = 0.0;
    weyl[2][1][0][2] = 0.0;
    weyl[2][1][0][3] = 0.0;
    weyl[2][1][1][0] = 0.0;
    weyl[2][1][1][1] = 0.0;
    weyl[2][1][1][2] = t19;
    weyl[2][1][1][3] = 0.0;
    weyl[2][1][2][0] = 0.0;
    weyl[2][1][2][1] = -t19;
    weyl[2][1][2][2] = 0.0;
    weyl[2][1][2][3] = 0.0;
    weyl[2][1][3][0] = 0.0;
    weyl[2][1][3][1] = 0.0;
    weyl[2][1][3][2] = 0.0;
    weyl[2][1][3][3] = 0.0;
    weyl[2][2][0][0] = 0.0;
    weyl[2][2][0][1] = 0.0;
    weyl[2][2][0][2] = 0.0;
    weyl[2][2][0][3] = 0.0;
    weyl[2][2][1][0] = 0.0;
    weyl[2][2][1][1] = 0.0;
    weyl[2][2][1][2] = 0.0;
    weyl[2][2][1][3] = 0.0;
    weyl[2][2][2][0] = 0.0;
    weyl[2][2][2][1] = 0.0;
    weyl[2][2][2][2] = 0.0;
    weyl[2][2][2][3] = 0.0;
    weyl[2][2][3][0] = 0.0;
    weyl[2][2][3][1] = 0.0;
    weyl[2][2][3][2] = 0.0;
    weyl[2][2][3][3] = 0.0;
    weyl[2][3][0][0] = 0.0;
    weyl[2][3][0][1] = 0.0;
    weyl[2][3][0][2] = 0.0;
    weyl[2][3][0][3] = 0.0;
    weyl[2][3][1][0] = 0.0;
    weyl[2][3][1][1] = 0.0;
    weyl[2][3][1][2] = 0.0;
    weyl[2][3][1][3] = 0.0;
    weyl[2][3][2][0] = 0.0;
    weyl[2][3][2][1] = 0.0;
    weyl[2][3][2][2] = 0.0;
    weyl[2][3][2][3] = t23;
    weyl[2][3][3][0] = 0.0;
    weyl[2][3][3][1] = 0.0;
    weyl[2][3][3][2] = -t23;
    weyl[2][3][3][3] = 0.0;
    weyl[3][0][0][0] = 0.0;
    weyl[3][0][0][1] = 0.0;
    weyl[3][0][0][2] = 0.0;
    weyl[3][0][0][3] = -t16;
    weyl[3][0][1][0] = 0.0;
    weyl[3][0][1][1] = 0.0;
    weyl[3][0][1][2] = 0.0;
    weyl[3][0][1][3] = 0.0;
    weyl[3][0][2][0] = 0.0;
    weyl[3][0][2][1] = 0.0;
    weyl[3][0][2][2] = 0.0;
    weyl[3][0][2][3] = 0.0;
    weyl[3][0][3][0] = t16;
    weyl[3][0][3][1] = 0.0;
    weyl[3][0][3][2] = 0.0;
    weyl[3][0][3][3] = 0.0;
    weyl[3][1][0][0] = 0.0;
    weyl[3][1][0][1] = 0.0;
    weyl[3][1][0][2] = 0.0;
    weyl[3][1][0][3] = 0.0;
    weyl[3][1][1][0] = 0.0;
    weyl[3][1][1][1] = 0.0;
    weyl[3][1][1][2] = 0.0;
    weyl[3][1][1][3] = t21;
    weyl[3][1][2][0] = 0.0;
    weyl[3][1][2][1] = 0.0;
    weyl[3][1][2][2] = 0.0;
    weyl[3][1][2][3] = 0.0;
    weyl[3][1][3][0] = 0.0;
    weyl[3][1][3][1] = -t21;
    weyl[3][1][3][2] = 0.0;
    weyl[3][1][3][3] = 0.0;
    weyl[3][2][0][0] = 0.0;
    weyl[3][2][0][1] = 0.0;
    weyl[3][2][0][2] = 0.0;
    weyl[3][2][0][3] = 0.0;
    weyl[3][2][1][0] = 0.0;
    weyl[3][2][1][1] = 0.0;
    weyl[3][2][1][2] = 0.0;
    weyl[3][2][1][3] = 0.0;
    weyl[3][2][2][0] = 0.0;
    weyl[3][2][2][1] = 0.0;
    weyl[3][2][2][2] = 0.0;
    weyl[3][2][2][3] = -t23;
    weyl[3][2][3][0] = 0.0;
    weyl[3][2][3][1] = 0.0;
    weyl[3][2][3][2] = t23;
    weyl[3][2][3][3] = 0.0;
    weyl[3][3][0][0] = 0.0;
    weyl[3][3][0][1] = 0.0;
    weyl[3][3][0][2] = 0.0;
    weyl[3][3][0][3] = 0.0;
    weyl[3][3][1][0] = 0.0;
    weyl[3][3][1][1] = 0.0;
    weyl[3][3][1][2] = 0.0;
    weyl[3][3][1][3] = 0.0;
    weyl[3][3][2][0] = 0.0;
    weyl[3][3][2][1] = 0.0;
    weyl[3][3][2][2] = 0.0;
    weyl[3][3][2][3] = 0.0;
    weyl[3][3][3][0] = 0.0;
    weyl[3][3][3][1] = 0.0;
    weyl[3][3][3][2] = 0.0;
    weyl[3][3][3][3] = 0.0;

    return true;
}

/*! Calculate the Ricci rotation coefficients at position 'pos'.
 *
 *  \param  pos  :  pointer to coordinate position where the Ricci rotation coefficients have to be evaluated.
 *  \return true :  successfull
 */
bool MetricSchwarzschild::calculateRicRotCoeffs(const double* pos) {
    double r     = pos[1];
    double theta = pos[2];

    double t1 = r - rs;
    double t3 = 1 / r;
    double t6 = sqrt(t1 * t3);
    double t9 = 1 / t1 * t3 * t6 * rs / 2.0;
    double t10 = r * r;
    double t14 = 1 / t10 * t1 / t6;
    double t15 = sin(theta);
    double t18 = cos(theta);
    double t19 = t3 / t15 * t18;

    rrc[0][0][0] = 0.0;
    rrc[0][0][1] = 0.0;
    rrc[0][0][2] = 0.0;
    rrc[0][0][3] = 0.0;
    rrc[0][1][0] = -t9;
    rrc[0][1][1] = 0.0;
    rrc[0][1][2] = 0.0;
    rrc[0][1][3] = 0.0;
    rrc[0][2][0] = 0.0;
    rrc[0][2][1] = 0.0;
    rrc[0][2][2] = 0.0;
    rrc[0][2][3] = 0.0;
    rrc[0][3][0] = 0.0;
    rrc[0][3][1] = 0.0;
    rrc[0][3][2] = 0.0;
    rrc[0][3][3] = 0.0;
    rrc[1][0][0] = t9;
    rrc[1][0][1] = 0.0;
    rrc[1][0][2] = 0.0;
    rrc[1][0][3] = 0.0;
    rrc[1][1][0] = 0.0;
    rrc[1][1][1] = 0.0;
    rrc[1][1][2] = 0.0;
    rrc[1][1][3] = 0.0;
    rrc[1][2][0] = 0.0;
    rrc[1][2][1] = 0.0;
    rrc[1][2][2] = -t14;
    rrc[1][2][3] = 0.0;
    rrc[1][3][0] = 0.0;
    rrc[1][3][1] = 0.0;
    rrc[1][3][2] = 0.0;
    rrc[1][3][3] = -t14;
    rrc[2][0][0] = 0.0;
    rrc[2][0][1] = 0.0;
    rrc[2][0][2] = 0.0;
    rrc[2][0][3] = 0.0;
    rrc[2][1][0] = 0.0;
    rrc[2][1][1] = 0.0;
    rrc[2][1][2] = t14;
    rrc[2][1][3] = 0.0;
    rrc[2][2][0] = 0.0;
    rrc[2][2][1] = 0.0;
    rrc[2][2][2] = 0.0;
    rrc[2][2][3] = 0.0;
    rrc[2][3][0] = 0.0;
    rrc[2][3][1] = 0.0;
    rrc[2][3][2] = 0.0;
    rrc[2][3][3] = -t19;
    rrc[3][0][0] = 0.0;
    rrc[3][0][1] = 0.0;
    rrc[3][0][2] = 0.0;
    rrc[3][0][3] = 0.0;
    rrc[3][1][0] = 0.0;
    rrc[3][1][1] = 0.0;
    rrc[3][1][2] = 0.0;
    rrc[3][1][3] = t14;
    rrc[3][2][0] = 0.0;
    rrc[3][2][1] = 0.0;
    rrc[3][2][2] = 0.0;
    rrc[3][2][3] = t19;
    rrc[3][3][0] = 0.0;
    rrc[3][3][1] = 0.0;
    rrc[3][3][2] = 0.0;
    rrc[3][3][3] = 0.0;

    return true;
}

/*! Calculate the contractions of the Ricci rotation coefficients at position 'pos'.
 *
 *  \param  pos  :  pointer to coordinate position where the contractions of the Ricci rotation coefficients have to be evaluated.
 *  \return true :  successfull
 */
bool MetricSchwarzschild::calculateContrRRC(const double* pos) {
    double r     = pos[1];
    double theta = pos[2];

    crrc[0] = crrc[3] = 0.0;
    crrc[1] = 0.5 * (4.0 * r - 3.0 * rs) / (r * r * sqrt(1.0 - rs / r));
    crrc[2] = cos(theta) / (r * sin(theta));
    return true;
}


/*! Transform local 4-direction to coordinate 4-direction.
 *
 *  \param  pos  :  pointer to position array.
 *  \param  ldir :  pointer to local direction array.
 *  \param  dir  :  pointer to calculated coordinate direction array.
 *  \param  type :  type of tetrad.
 */
void MetricSchwarzschild::localToCoord(const double* pos, const double* ldir, double* dir,
                                       enum_nat_tetrad_type) {
    double r     = pos[1];
    double theta = pos[2];
    double w = sqrt(1.0 - rs / r);

    dir[0] = ldir[0] / w / mSpeedOfLight;
    dir[1] = ldir[1] * w;
    dir[2] = ldir[2] / r;
    dir[3] = ldir[3] / (r * sin(theta));
}

/*! Transform coordinate 4-direction to local 4-direction.
 *
 *  \param  pos  :  pointer to position array.
 *  \param  cdir :  pointer to coordinate direction.
 *  \param  ldir :  pointer to calculated local direction array.
 *  \param  type :  type of tetrad.
 */
void MetricSchwarzschild::coordToLocal(const double* pos, const double* cdir, double* ldir,
                                       enum_nat_tetrad_type) {
    double r     = pos[1];
    double theta = pos[2];
    double w = sqrt(1.0 - rs / r);

    ldir[0] = cdir[0] * w * mSpeedOfLight;
    ldir[1] = cdir[1] / w;
    ldir[2] = cdir[2] * r;
    ldir[3] = cdir[3] * r * sin(theta);
}


/*! Test break condition.
 *
 *  \param pos    : pointer to position array.
 *  \return true  : radial position r < 0.0 or  r^2<=(1.0+eps)*rs^2.
 *  \return false : position is valid.
 */
bool MetricSchwarzschild::breakCondition(const double* pos) {
    bool br = false;

    if ((pos[1] < 0.0) || (pos[1]*pos[1] <= (1.0 + eps)*rs * rs)) {
        br = true;
    }
    return br;
}

/*! Calculate right hand side of the geodesic equation in first order form.
 *
 *  \param  y[]   : pointer to position and direction coordinates.
 *  \param  dydx[] : pointer to right side of geodesic equation.
 */
bool MetricSchwarzschild::calcDerivs(const double y[], double dydx[]) {
    dydx[0] = y[4];
    dydx[1] = y[5];
    dydx[2] = y[6];
    dydx[3] = y[7];

    double r     = y[1];
    double theta = y[2];

    dydx[4] = -rs / (r * (r - rs)) * y[4] * y[5];
    dydx[5] = -0.5 * mSpeedOfLight * mSpeedOfLight * rs * (r - rs) / pow(r, 3.0) * y[4] * y[4] + 0.5 * rs / (r * (r - rs)) * y[5] * y[5] + (r - rs) * (y[6] * y[6] + sin(theta) * sin(theta) * y[7] * y[7]);
    dydx[6] = -2.0 / r * y[5] * y[6] + sin(theta) * cos(theta) * y[7] * y[7];
    dydx[7] = -2.0 / r * y[5] * y[7] - 2.0 * cos(theta) / sin(theta) * y[6] * y[7];

    return true;
}

/*! Calculate right hand side of the geodesic equation in first order form with parallel transport.
 *
 *  \param  y[]   : pointer to position and direction coordinates.
 *  \param  dydx[] : pointer to right side of parallel transport equation.
 */
bool MetricSchwarzschild::calcDerivsPar(const double y[], double dydx[]) {
    calcDerivs(y, dydx);

    double r     = y[1];
    double theta = y[2];

    int bidx;
    for (int n = 0; n < 4; n++) {
        bidx = 4 * (n + 2);
        dydx[bidx + 0] = -0.5 * rs / (r * (r - rs)) * (y[4] * y[bidx + 1] + y[5] * y[bidx + 0]);
        dydx[bidx + 1] = -0.5 * rs * (r - rs) / (r * r * r) * mSpeedOfLight * mSpeedOfLight * y[4] * y[bidx + 0] + 0.5 * rs / (r * (r - rs)) * y[5] * y[bidx + 1] + (r - rs) * (y[6] * y[bidx + 2] + pow(sin(theta), 2.0) * y[7] * y[bidx + 3]);
        dydx[bidx + 2] = -y[5] / r * y[bidx + 2] - y[6] / r * y[bidx + 1] + sin(theta) * cos(theta) * y[7] * y[bidx + 3];
        dydx[bidx + 3] = -y[5] * y[bidx + 3] / r - y[7] * y[bidx + 1] / r - cos(theta) / sin(theta) * (y[6] * y[bidx + 3] + y[7] * y[bidx + 2]);
    }

    return true;
}

/*! Calculate right hand side of parallel transport and Jocobi equation.
 */
bool MetricSchwarzschild::calcDerivsSachsJacobi(const double y[], double dydx[]) {
    const double* u    = &y[DEF_TG_IDX];
    const double* wSA1 = &y[DEF_SA1_IDX];
    const double* wSA2 = &y[DEF_SA2_IDX];
    const double* wJ1  = &y[DEF_JAC1_IDX];
    const double* wJ2  = &y[DEF_JAC2_IDX];
    const double* wJ1d = &y[DEF_DJ1_IDX];
    const double* wJ2d = &y[DEF_DJ2_IDX];

    double gd[4];
    contrChrisVecVec(y, u, u, gd);

    double zSA1[4], zSA2[4], zJ1a[4], zJ1b[4], zJ2a[4], zJ2b[4];
    contrChrisVecVec(y, u, wSA1, zSA1, false);
    contrChrisVecVec(y, u, wSA2, zSA2, false);

    contrChrisVecVec(y, u, wJ1d, zJ1a, false);
    contrChrisVecVec(y, u, wJ2d, zJ2a, false);

    contrChrDVecVecVec(y, u, u, wJ1, zJ1b);
    contrChrDVecVecVec(y, u, u, wJ2, zJ2b, false);

    for (int mu = 0; mu < 4; mu++) {
        dydx[mu]              = y[DEF_TG_IDX + mu];
        dydx[DEF_TG_IDX + mu]   = -gd[mu];

        dydx[DEF_SA1_IDX + mu]  = -zSA1[mu];
        dydx[DEF_SA2_IDX + mu]  = -zSA2[mu];
        dydx[DEF_JAC1_IDX + mu] = y[DEF_DJ1_IDX + mu];
        dydx[DEF_JAC2_IDX + mu] = y[DEF_DJ2_IDX + mu];

        dydx[DEF_DJ1_IDX + mu]  = -2.0 * zJ1a[mu] - zJ1b[mu];
        dydx[DEF_DJ2_IDX + mu]  = -2.0 * zJ2a[mu] - zJ2b[mu];
    }
    return true;
}


/*! Calculate right hand side of the Fermi-Walker transport equation.
 *
 *  \param  a[]  : pointer to proper acceleration.
 *  \param  y[]  : pointer to position, direction and tetrad.
 *  \param dydx[] : pointer to right side of Fermi-Walker transport equation.
 *  \return true : always.
 */
bool MetricSchwarzschild::calcDerivsFW(const double a[], const double y[], double dydx[]) {
    //return false;
    dydx[0] = y[4];
    dydx[1] = y[5];
    dydx[2] = y[6];
    dydx[3] = y[7];

    double r     = y[1];
    double theta = y[2];

    // double edc = 1.0/mSpeedOfLight;
    double edr = 1.0 / r;

    int bidx;  // tetrad index
    for (int n = 0; n < 4; n++) {
        bidx = 4 * (n + 2);
        dydx[bidx + 0] = -0.5 * rs / (r * (r - rs)) * (y[4] * y[bidx + 1] + y[5] * y[bidx + 0])  +  a[n] * y[4];
        dydx[bidx + 1] = -0.5 * rs * (r - rs) / (r * r * r) * mSpeedOfLight * mSpeedOfLight * y[4] * y[bidx + 0] + 0.5 * rs / (r * (r - rs)) * y[5] * y[bidx + 1] + (r - rs) * (y[6] * y[bidx + 2] + pow(sin(theta), 2.0) * y[7] * y[bidx + 3])  +  a[n] * y[5];
        dydx[bidx + 2] = -y[5] * y[bidx + 2] * edr - y[6] * y[bidx + 1] * edr + sin(theta) * cos(theta) * y[7] * y[bidx + 3]  +  a[n] * y[6];
        dydx[bidx + 3] = -y[5] * y[bidx + 3] * edr - y[7] * y[bidx + 1] * edr - cos(theta) / sin(theta) * (y[6] * y[bidx + 3] + y[7] * y[bidx + 2])  +  a[n] * y[7];
    }

    bidx = 8;
    for (int mu = 0; mu < 4; mu++) {
        dydx[bidx + mu] += a[1] * y[DEF_E1_IDX + mu] + a[2] * y[DEF_E2_IDX + mu] + a[3] * y[DEF_E3_IDX + mu];
        dydx[mu + 4] = dydx[bidx + mu];
    }
    return true;
}

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
double MetricSchwarzschild::testConstraint(const double y[], const double kappa) {
    double r     = y[1];
    double theta = y[2];
    double cm = 1.0 / mSpeedOfLight;

    // Scale the directions with the speed of light before doubling them !!
    double dt = y[4];
    double dr = y[5] * cm;
    double dth = y[6] * cm;
    double dph = y[7] * cm;

    double sum = -kappa;
    sum += -(1.0 - rs / r) * dt * dt + dr * dr / (1.0 - rs / r) + r * r * (dth * dth + sin(theta) * sin(theta) * dph * dph);
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
bool MetricSchwarzschild::calcProduct(const double* pos, const double* u, const double* v, double &prod, bool) {
    prod = 0.0;
    if (breakCondition(pos)) {
        return false;
    }

    double r     = pos[1];
    double theta = pos[2];
    prod = -mSpeedOfLight * mSpeedOfLight * (1.0 - rs / r) * u[0] * v[0] + u[1] * v[1] / (1.0 - rs / r) + r * r * (u[2] * v[2] + sin(theta) * sin(theta) * u[3] * v[3]);
    return true;
}

/*! Set parameter 'pName' to 'val'.
 *
 *  Set 'mass' parameter and adjust Schwarzschild radius  rs=2GM/c^2.
 */
bool MetricSchwarzschild::setParam(std::string pName, double val) {
    if (Metric::setParam(pName, val)) {
        mMass = val;
        rs = 2.0 * mGravConstant * mMass / (mSpeedOfLight * mSpeedOfLight);
    }
    return true;
}

/*! Transform point p to embedding coordinates.
 *
 *  \param p  : point to be transformed.
 *  \param ep : reference to 'embedded' point.
 *  \return true : success.
 *  \return false : otherwise.
 */
bool MetricSchwarzschild::transToEmbedding(vec4 p, vec4 &ep) {
    vec4 cp;
    transToPseudoCart(p, cp);

    double r = p[1];
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

/*! Transform point p to 2+1 coordinates.
 *
 *  \param  p  : point in proper metric coordinates.
 *  \param  cp : reference to transformed point.
 *  \return true : success.
 */
bool MetricSchwarzschild::transToTwoPlusOne(vec4 p, vec4 &cp) {
    vec4 tp;
    TransCoordinates::toCartesianCoord(mCoordType, p, tp);
    cp = vec4(tp[0], tp[1], tp[2], tp[0]);
    return true;
}

/*!  Transform point p to custom coordinats.
 *
 *  \param p  : point to be transformed.
 *  \param cp : reference to customized point.
 *  \return true : always.
 */
bool MetricSchwarzschild::transToCustom(vec4 p, vec4 &cp) {
    double r = p[1];
    double phi = p[3];

    cp[0] = p[0];
    cp[1] = (10.0 - 9.0 * rs / r) * cos(phi);
    cp[2] = (10.0 - 9.0 * rs / r) * sin(phi);
    cp[3] = 0.0;
    return true;
}

/*! Set embedding parameters.
 *
 *  \param  name : embedding parameter name.
 *  \param  val  : embedding parameter value.
 *  \return true  : success.
 *  \return false : parameter not valid.
 */
bool MetricSchwarzschild::setEmbeddingParam(std::string name, double val) {
    Metric::setEmbeddingParam(name, val);

    if (name == "emb_rmin") {
        mEmb_rmin = val;
    } else if (name == "emb_rmax") {
        mEmb_rmax = val;
    } else if (name == "emb_r_num") {
        mEmb_r_num = val;
    } else if (name == "emb_phi_num") {
        mEmb_phi_num = val;
    }
    return testEmbeddingParams();
}

/*! Test embedding parameters
 *  \return  true : all parameters are ok
 *  \return  false : at least one parameter had to be adjusted.
 */
bool MetricSchwarzschild::testEmbeddingParams() {
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
int MetricSchwarzschild::getEmbeddingVertices(std::vector<vec3> &verts,
        std::vector<int> &indices, unsigned int &numElems, unsigned int &counter) {
    if (!verts.empty()) {
        verts.clear();
    }

    if (!indices.empty()) {
        indices.clear();
    }

    testEmbeddingParams();
    mEmb_rstep = (mEmb_rmax - mEmb_rmin) / mEmb_r_num;
    mEmb_phistep = 2.0 * M_PI / mEmb_phi_num;

    numElems = int(mEmb_r_num);
    counter  = int(mEmb_phi_num) + 1;

    int vnum;

    double x, y, z, r, phi;
    for (unsigned int k = 0; k < counter; k++) {
        phi = k * mEmb_phistep;
        for (unsigned int j = 0; j < numElems; j++) {
            r = mEmb_rmin + j * mEmb_rstep;
            x = r * cos(phi);
            y = r * sin(phi);
            if (r >= rs) {
                z = 2.0 * sqrt(rs) * sqrt(fabs(r - rs));
                verts.push_back(vec3(x, y, z));

                vnum = k * numElems + j;

                indices.push_back(vnum);
                indices.push_back(vnum + numElems);
            }
        }
    }

    int numVerts = (int)verts.size();
    int numInds  = (int)indices.size();

    if (2 * numVerts == numInds) {
        return numVerts;
    }

    return 0;
}

/*! Effective potential.
 *  \param pos : initial position.
 *  \param cdir : initial four-direction.
 *  \param type : geodesic type.
 *  \param x : abscissa value.
 *  \param val : reference to effective potential value.
 *  \return true : effective potential exists at x.
 */
bool MetricSchwarzschild::effPotentialValue(const vec4 pos, const vec4 cdir, enum_geodesic_type type, const double x, double &val) {
    double kappa = 0.0;
    if (type == enum_geodesic_timelike) {
        kappa = -mSign;
    }

    if (pos[1] < rs + 1e-2 || x < rs + 1e-2) {
        return false;
    }

    double h = pos[1] * pos[1] * cdir[3];
    val = 0.5 * (1.0 - rs / x) * (h * h / (x * x) - kappa * mSpeedOfLight * mSpeedOfLight);
    return true;
}

/*! Total energy.
 *  \param pos : initial position.
 *  \param cdir : initial four-direction.
 *  \param x : abscissa value.
 *  \param val : reference to total energy value.
 *  \return true : effective potential exists at x.
 */
bool MetricSchwarzschild::totEnergy(const vec4 pos, const vec4 cdir, const double , double &val) {
    if (pos[1] < rs + 1e-2) {
        return false;
    }

    // 1/2*k^2/c^2:
    val = 0.5 * (1.0 - rs / pos[1]) * (1.0 - rs / pos[1]) * mSpeedOfLight * mSpeedOfLight * cdir[0] * cdir[0];
    return true;
}

/*! Determine the velocity for a closed circular orbit if it exists.
 *   A circular timelike geodesic with respect to r-coordinate does exist
 *   only for r>=3rs (last timelike circular orbit).
 * \param r  Radial coordinate.
 * \param tedType type of tetrad.
 */
double MetricSchwarzschild::getCircularVelocity(const double r, const enum_nat_tetrad_type) {
    if (r >= 3.0 * rs) {
        return 1.0 / sqrt(2.0 * (r / rs - 1.0));
    }
    return 0.0;
}

/**
 * @brief MetricSchwarzschild::getCircularFourVel
 * @param pos
 * @param tedType type of tetrad
 * @return
 */
vec4
MetricSchwarzschild::getCircularFourVel(const vec4 pos, const enum_nat_tetrad_type) {
    double beta = getCircularVelocity(pos[1]);
    if (beta > 0.0 && beta < 1.0) {
        double gamma = 1.0 / sqrt(1.0 - beta * beta);
        vec4 e0, e1, e2, e3;
        getNatTetrad(pos, e0, e1, e2, e3);
        return mSpeedOfLight * gamma * (e0 + beta * e3);
    }
    return vec4();
}

/*!
 *  \param units : type of physical constants.
 */
void MetricSchwarzschild::usePhysicalUnits(const enum_physical_constants  units) {
    Metric::usePhysicalUnits(units);
    rs = 2.0 * mGravConstant * mMass / (mSpeedOfLight * mSpeedOfLight);
}

/*!
 *  \param speed_of_light : value for speed of light.
 *  \param grav_const : value for gravitational constant.
 *  \param diel_perm : value for dielectric permittivity.
 */
void MetricSchwarzschild::setUnits(const double speed_of_light, const double grav_const, const double diel_perm) {
    Metric::setUnits(speed_of_light, grav_const, diel_perm);
    rs = 2.0 * mGravConstant * mMass / (mSpeedOfLight * mSpeedOfLight);
}

/*! Generate report.
 * \param pos : initial position.
 * \param cdir : initial coordinate direction.
 * \param text : reference to report text.
 */
bool MetricSchwarzschild::report(const vec4 pos, const vec4 cdir, std::string &text) {
    std::stringstream ss;
    ss << "Report for Schwarzschild metric\n\tcoordinates : (t,r,theta,phi)\n";
    ss << "---------------------------------------------------------------\n";
    ss << "  physical units ................................. yes\n";
    ss.precision(DEF_FIXED_REPORT_PRECISION);
    ss.setf(std::ios::fixed);
    ss << "  Schwarzschild radius ........... r_s = 2GM/c^2 = " << rs << std::endl;
    ss << "  Photon orbit ................... r_ph = 3/2*rs = " << 1.5 * rs << std::endl;
    ss << "  innermost stable circular orbit  r_isco = 3r_s = " << 3.0 * rs << std::endl;
    ss << "                                            beta = " << getCircularVelocity(pos[1]) << std::endl;

    double k = (1.0 - rs / pos[1]) * mSpeedOfLight * mSpeedOfLight * cdir[0];
    double h = pos[1] * pos[1] * cdir[3];
    ss.precision(DEF_FIXED_REPORT_PRECISION);
    ss << "  constant of motion ..................... k = " << k << std::endl;
    ss << "  constant of motion ..................... h = " << h << std::endl;
    ss << "  critical angle for null geodesic . ksiCrit = ";

    double ksicrit;
    if (calcKsiCrit(pos, ksicrit)) {
        ss << ksicrit*RAD_TO_DEG << std::endl;
    } else {
        ss << "not valid here.";
    }

    ss << "\nThe effective potential is only valid for geodesics in the theta=pi/2 hypersurface." << std::endl;
    text = ss.str();
    return true;
}

// ***************************** specific public methods ***************************
/*!  Calculate the critical angle.
 *
 *   For an observer at distance \f$x_i=r_s/r_i\f$, r_i=pos[1], the black hole has an angular diameter of \f$\xi_{\mbox{crit}}\f$ with
 *    \f[ \xi_{\mbox{crit}} = \arcsin\sqrt{\frac{27}{4}x_i^2(1-x_i)} \f] .
 *
 *  \param pos : current position
 *  \param ksicrit : reference to critical angle.
 */
bool MetricSchwarzschild::calcKsiCrit(const vec4 pos, double &ksicrit) {
    if (pos[1] < rs) {
        return false;
    }

    double xi = rs / pos[1];

    if (pos[1] >= 1.5 * rs) {
        ksicrit = asin(sqrt(6.75 * xi * xi * (1.0 - xi)));
    } else {
        ksicrit = M_PI - asin(sqrt(6.75 * xi * xi * (1.0 - xi)));
    }

    return true;
}

/*! Calculate constant of motion k
 *
 * \param r
 * \param beta
 * \param kappa
 */
double MetricSchwarzschild::getConstOfMotion_k(const double r, const double beta, const int kappa) {
    double gamma = 1.0 / sqrt(1.0 - beta * beta);
    double eta = 1.0;

    if (kappa == -1) {
        eta = beta * gamma;
    } else if (kappa == 0) {
        eta = 1.0;
    }
    return sqrt(eta * eta - kappa) * sqrt(1.0 - 2.0 * mMass / r);
}

double MetricSchwarzschild::calcKsiCrit(const double r) {
    double rs = 2.0 * mMass;
    double rsdri = rs / r;
    double sk = sqrt(6.75 * rsdri * rsdri * (1.0 - rsdri));
    if (r >= 3.0 * mMass) {
        return M_PI - asin(sk);
    }
    return asin(sk);
}

/*! Calculate constant of motion h
 *
 * \param r
 * \param beta
 * \param ksi   in rad
 * \param kappa
 */
double MetricSchwarzschild::getConstOfMotion_h(const double r, const double beta, const double ksi, const int kappa) {
    double gamma = 1.0 / sqrt(1.0 - beta * beta);
    double eta = 1.0;

    if (kappa == -1) {
        eta = beta * gamma;
    } else if (kappa == 0) {
        eta = 1.0;
    }
    return r * eta * sin(ksi);
}


double MetricSchwarzschild::param_aqua(const double k, const double h, const int kappa) {
    double rs = 2.0 * mMass;
    return rs * rs * k * k / (h * h) + param_b(k, h, kappa);
}

double MetricSchwarzschild::param_b(const double , const double h, const int kappa) {
    if (kappa == 0) {
        return 0.0;
    }
    double rs = 2.0 * mMass;
    return rs * rs * double(kappa) / (h * h);
}

void MetricSchwarzschild::param_sq(const double k, const double h, const int kappa,
                                   double &p, double &q, double &D1) {
    p = -param_b(k, h, kappa) / 3.0 - 1.0 / 9.0;
    q = 0.5 * param_aqua(k, h, kappa) - 1.0 / 27.0 - param_b(k, h, kappa) / 6.0;
    D1 = q * q + p * p * p;
}

double MetricSchwarzschild::param_psi(const double p, const double q) {
    double psi = 0.0;
    double D1 = q * q + p * p * p;

    double rho3 = pow(sign(q) * sqrt(fabs(p)), 3.0);
    if ((p < 0.0) && (D1 <= 0.0)) {
        psi = acos(q / rho3);
    } else if ((p < 0.0) && (D1 > 0.0)) {
        psi = acosh(q / rho3);
        //    cerr << "psi=" << psi << endl;
    } else if (p > 0.0) {
        psi = asinh(q / rho3);
    } else if (fabs(p) < 1.0e-10) {
        psi = 0.0;
    }
    return psi;
}

double MetricSchwarzschild::param_eps(const double r, const double ksi) {
    return r * sin(ksi) / sqrt(1.0 - 2.0 * mMass / r);
}


/*! calculate minimal distance to black hole
 *
 * @param r    initial radial distance of geodesic
 * @param ksi  initial direction of geodesic in rad
 * @param rMin  minimal distance of geodesic
 * @return true : minimal distance exists\n
 *         false : minimal distance does not exist
 */
bool MetricSchwarzschild::minDist_r(const double r, const double ksi, double &rMin) {
    double b = param_eps(r, ksi);
    double sq3 = sqrt(3.0);
    double val = -3.0 * sq3 * mMass / b;
    if (fabs(val) > 1.0) {
        //std::cerr << "rMin does not exist!" << std::endl;
        return false;
    }

    double psidr = acos(val) / 3.0;
    rMin = 2.0 * b / sq3 * cos(psidr);
    return true;
}

// ********************************* protected methods *****************************
/*!
 */
void MetricSchwarzschild::setStandardValues() {
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

/*! Contract Christoffel symbol with two vectors.
 *  \param y : pointer to full data
 *  \param u : pointer to first vector
 *  \param w : pointer to second vector
 *  \param z : pointer to contraction
 *  \param calc : calculate the Christoffels before using them
 *
 *  The Christoffel symbols do not have to be calculated before this function can be evaluated!
 */
void MetricSchwarzschild::contrChrisVecVec(const double y[], const double u[], const double w[], double* z, bool calc) {
    double r     = y[1];
    double theta = y[2];

    if (calc) {
        christoffel[0][0][1] = 0.5 * (r - rs) * mSpeedOfLight * mSpeedOfLight * rs * pow(r, -3.0);
        christoffel[0][1][0] = 0.5 * rs / (r * (r - rs));
        christoffel[1][1][1] = -christoffel[0][1][0];
        christoffel[1][2][2] = christoffel[1][3][3] = 1.0 / r;
        christoffel[2][2][1] = -r + rs;
        christoffel[2][3][3] = cos(theta) / sin(theta);
        christoffel[3][3][1] = -(r - rs) * sin(theta) * sin(theta);
        christoffel[3][3][2] = -sin(theta) * cos(theta);
    }

    z[0] = christoffel[0][1][0] * (u[1] * w[0] + u[0] * w[1]);
    z[1] = christoffel[0][0][1] * u[0] * w[0] + christoffel[1][1][1] * u[1] * w[1] + christoffel[2][2][1] * u[2] * w[2] + christoffel[3][3][1] * u[3] * w[3];
    z[2] = christoffel[1][2][2] * (u[2] * w[1] + u[1] * w[2]) + christoffel[3][3][2] * u[3] * w[3];
    z[3] = christoffel[1][3][3] * (u[3] * w[1] + u[1] * w[3]) + christoffel[2][3][3] * (u[3] * w[2] + u[2] * w[3]);
}

/*! Contract partially derived Christoffel symbols with three vectors.
 *  \param y : pointer to full data
 *  \param u : pointer to first vector
 *  \param v : pointer to second vector
 *  \param w : pointer to third vector
 *  \param z : pointer to contraction
 *  \param calc : calculate the Christoffels before using them
 *
 * The partial derivatives of the Christoffel symbols do not have to be calculated before this function can be evaluated!
 */
void MetricSchwarzschild::contrChrDVecVecVec(const double y[], const double u[], const double v[], const double w[], double* z, bool calc) {
    double r     = y[1];
    double theta = y[2];

    if (calc) {
        chrisD[0][0][1][1] = -0.5 * mSpeedOfLight * mSpeedOfLight * rs * (2.0 * r - 3.0 * rs) * pow(r, -4.0);
        chrisD[0][1][0][1] = -0.5 * rs * (2.0 * r - rs) * pow(r * (r - rs), -2.0);
        chrisD[1][1][1][1] = -chrisD[0][1][0][1];
        chrisD[1][2][2][1] = chrisD[1][3][3][1] = -1.0 / (r * r);
        chrisD[2][2][1][1] = -1.0;
        chrisD[3][3][1][1] = -sin(theta) * sin(theta);
        chrisD[2][3][3][2] = 1.0 / chrisD[3][3][1][1];
        chrisD[3][3][1][2] = -(r - rs) * sin(2.0 * theta);
        chrisD[3][3][2][2] = -cos(2.0 * theta);
    }

    z[0] = chrisD[0][1][0][1] * (u[1] * v[0] * w[1] + u[0] * v[1] * w[1]);
    z[1] = chrisD[0][0][1][1] * u[0] * v[0] * w[1] + chrisD[1][1][1][1] * u[1] * v[1] * w[1] + chrisD[2][2][1][1] * u[2] * v[2] * w[1] + chrisD[3][3][1][1] * u[3] * v[3] * w[1] + chrisD[3][3][1][2] * u[3] * v[3] * w[2];
    z[2] = chrisD[1][2][2][1] * (u[2] * v[1] * w[1] + u[1] * v[2] * w[1]) + chrisD[3][3][2][2] * u[3] * v[3] * w[2];
    z[3] = chrisD[1][3][3][1] * (u[3] * v[1] * w[1] + u[1] * v[3] * w[1]) + chrisD[2][3][3][2] * (u[3] * v[2] * w[2] + u[2] * v[3] * w[2]);
}

} // end namespace m4d
