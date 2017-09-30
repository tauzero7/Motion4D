// -------------------------------------------------------------------------------
/*
   m4dMetricReissnerNordstrom.cpp

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

#include "m4dMetricReissnerNordstrom.h"

double dzdr_reissner(double x, void* params) {
    struct_reissner_params* par = (struct_reissner_params*)params;
    double rs  = 2.0 * par->mass;
    double rho = par->rho;
    double q   = par->q;

    double ARN = 1.0 - rs / x + rho * q * q / (x * x);
    double dzdr2 = 1.0 / ARN - 1.0;
    return sqrt(fabs(dzdr2));
}

namespace m4d {

#define eps 1.0e-6


/*! Standard constructor for the Schwarzschild metric.
 *
 * \param  mass : mass of the black hole.
 * \param  q : charge
 */
MetricReissnerNordstrom::MetricReissnerNordstrom(double mass, double q) {
    mMetricName  = "ReissnerNordstrom";
    setCoordType(enum_coordinate_spherical);

    mPhysicalUnits = enum_physical_constants_geom;
    mSpeedOfLight   = 1.0;
    mGravConstant   = 1.0;
    mDielectricPerm = 1.0;
    mK = mGravConstant / (mDielectricPerm * pow(mSpeedOfLight, 4.0));

    addParam("mass", mass);
    mMass = mass;
    rs = 2.0 * mGravConstant * mMass / (mSpeedOfLight * mSpeedOfLight);

    addParam("charge");
    setParam("charge", q);
    mQ = q;
    calcDiskr();
    calcCritical();
    calcExtremals();

    mLocTeds.push_back(enum_nat_tetrad_static);

    mDrawTypes.push_back(enum_draw_embedding);
    mDrawTypes.push_back(enum_draw_effpoti);

    /*  parameters for the embedding diagram  */
    if (!mEmbParam.empty()) {
        mEmbParam.clear();
    }
    mHaveEmbedding = true;

    mEmb_rmin    = rp;
    mEmb_rmax    = 5.0 * rp;
    mEmb_r_num   = 20.0;
    mEmb_phi_num = 40.0;
    mEmb_rstep = (mEmb_rmax - mEmb_rmin) / mEmb_r_num;
    mEmb_phistep = 2.0 * M_PI / mEmb_phi_num;
    addEmbeddingParam("emb_rmin", mEmb_rmin);
    addEmbeddingParam("emb_rmax", mEmb_rmax);
    addEmbeddingParam("emb_r_num", 20.0);
    addEmbeddingParam("emb_phi_num", 40.0);

    w = gsl_integration_workspace_alloc(1000);
    F.function = &dzdr_reissner;

    setStandardValues();
}

/*!
 */
MetricReissnerNordstrom::~MetricReissnerNordstrom() {
    gsl_integration_workspace_free(w);
}


// *********************************** public methods ******************************

/*! Calculate the contravariant metric components at position 'pos'.
 *
 *  \param pos : pointer to position.
 */
bool MetricReissnerNordstrom::calculateMetric(const double* pos) {
    double r     = pos[1];
    double theta = pos[2];

    double c = mSpeedOfLight;

    double t1 = c * c;
    double t3 = 1 / r;
    double t6 = mQ * mQ;
    double t7 = r * r;
    double t8 = 1 / t7;
    double t17 = sin(theta);
    double t18 = t17 * t17;

    g_compts[0][0] = -t1 + t1 * rs * t3 - t1 * mK * t6 * t8;
    g_compts[0][1] = 0.0;
    g_compts[0][2] = 0.0;
    g_compts[0][3] = 0.0;
    g_compts[1][0] = 0.0;
    g_compts[1][1] = 1 / (1.0 - rs * t3 + mK * t6 * t8);
    g_compts[1][2] = 0.0;
    g_compts[1][3] = 0.0;
    g_compts[2][0] = 0.0;
    g_compts[2][1] = 0.0;
    g_compts[2][2] = t7;
    g_compts[2][3] = 0.0;
    g_compts[3][0] = 0.0;
    g_compts[3][1] = 0.0;
    g_compts[3][2] = 0.0;
    g_compts[3][3] = t7 * t18;

    return true;
}

/*! Calculate the Christoffel symbols of the second kind at position 'pos'.
 *
 *  \param pos : pointer to position.
 */
bool MetricReissnerNordstrom::calculateChristoffels(const double* pos) {
    double r     = pos[1];
    double theta = pos[2];

    double c = mSpeedOfLight;

    double t1 = r * r;
    double t2 = rs * r;
    double t3 = mQ * mQ;
    double t4 = mK * t3;
    double t5 = t1 - t2 + t4;
    double t6 = t1 * t1;
    double t10 = c * c;
    double t12 = t2 - 2.0 * t4;
    double t16 = 1 / r;
    double t20 = t16 / t5 * t12 / 2.0;
    double t21 = t5 * t16;
    double t22 = sin(theta);
    double t24 = cos(theta);
    double t25 = 1 / t22 * t24;
    double t26 = t22 * t22;

    christoffel[0][0][0] = 0.0;
    christoffel[0][0][1] = t5 / t6 / r * t10 * t12 / 2.0;
    christoffel[0][0][2] = 0.0;
    christoffel[0][0][3] = 0.0;
    christoffel[0][1][0] = t20;
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
    christoffel[1][0][0] = t20;
    christoffel[1][0][1] = 0.0;
    christoffel[1][0][2] = 0.0;
    christoffel[1][0][3] = 0.0;
    christoffel[1][1][0] = 0.0;
    christoffel[1][1][1] = -t20;
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
    christoffel[2][2][1] = -t21;
    christoffel[2][2][2] = 0.0;
    christoffel[2][2][3] = 0.0;
    christoffel[2][3][0] = 0.0;
    christoffel[2][3][1] = 0.0;
    christoffel[2][3][2] = 0.0;
    christoffel[2][3][3] = t25;
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
    christoffel[3][2][3] = t25;
    christoffel[3][3][0] = 0.0;
    christoffel[3][3][1] = -t21 * t26;
    christoffel[3][3][2] = -t22 * t24;
    christoffel[3][3][3] = 0.0;

    return true;
}

/*! Calculate Jacobi matrix.
 *
 *  \param pos : pointer to position.
 */
bool MetricReissnerNordstrom::calculateChrisD(const double* pos) {
    double r     = pos[1];
    double theta = pos[2];

    double c = mSpeedOfLight;

    double t1 = c * c;
    double t2 = r * r;
    double t5 = 2.0 * rs * t2 * r;
    double t7 = mQ * mQ;
    double t9 = 6.0 * t2 * mK * t7;
    double t10 = rs * rs;
    double t11 = t10 * t2;
    double t13 = rs * r;
    double t14 = mK * t7;
    double t15 = t13 * t14;
    double t17 = mK * mK;
    double t18 = t7 * t7;
    double t19 = t17 * t18;
    double t23 = t2 * t2;
    double t31 = 1 / t2;
    double t33 = t2 - t13 + t14;
    double t34 = t33 * t33;
    double t37 = (t5 - t9 - t11 + 4.0 * t15 - 2.0 * t19) * t31 / t34 / 2.0;
    double t38 = t2 - t14;
    double t40 = cos(theta);
    double t41 = t40 * t40;
    double t42 = sin(theta);
    double t43 = t42 * t42;
    double t46 = (t41 + t43) / t43;

    chrisD[0][0][0][0] = 0.0;
    chrisD[0][0][0][1] = 0.0;
    chrisD[0][0][0][2] = 0.0;
    chrisD[0][0][0][3] = 0.0;
    chrisD[0][0][1][0] = 0.0;
    chrisD[0][0][1][1] = -t1 * (t5 - t9 - 3.0 * t11 + 12.0 * t15 - 10.0 * t19) / t23 / t2 / 2.0;
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
    chrisD[0][1][0][1] = -t37;
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
    chrisD[1][0][0][1] = -t37;
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
    chrisD[1][1][1][1] = t37;
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
    chrisD[1][2][2][1] = -t31;
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
    chrisD[1][3][3][1] = -t31;
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
    chrisD[2][1][2][1] = -t31;
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
    chrisD[2][2][1][1] = -t38 * t31;
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
    chrisD[2][3][3][2] = -t46;
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
    chrisD[3][1][3][1] = -t31;
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
    chrisD[3][2][3][2] = -t46;
    chrisD[3][2][3][3] = 0.0;
    chrisD[3][3][0][0] = 0.0;
    chrisD[3][3][0][1] = 0.0;
    chrisD[3][3][0][2] = 0.0;
    chrisD[3][3][0][3] = 0.0;
    chrisD[3][3][1][0] = 0.0;
    chrisD[3][3][1][1] = -t43 * t38 * t31;
    chrisD[3][3][1][2] = -2.0 * t33 / r * t42 * t40;
    chrisD[3][3][1][3] = 0.0;
    chrisD[3][3][2][0] = 0.0;
    chrisD[3][3][2][1] = 0.0;
    chrisD[3][3][2][2] = -t41 + t43;
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
 *  \param  type :  type of tetrad.
 */
void MetricReissnerNordstrom::localToCoord(const double* pos, const double* ldir, double* dir,
        enum_nat_tetrad_type) {
    double r     = pos[1];
    double theta = pos[2];

    double As = sqrt(1.0 - rs / pos[1] + mK * mQ * mQ / (pos[1] * pos[1]));

    dir[0] = ldir[0] / As / mSpeedOfLight;
    dir[1] = ldir[1] * As;
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
void MetricReissnerNordstrom::coordToLocal(const double* pos, const double* cdir, double* ldir,
        enum_nat_tetrad_type) {
    double r     = pos[1];
    double theta = pos[2];

    double As = sqrt(1.0 - rs / pos[1] + mK * mQ * mQ / (pos[1] * pos[1]));

    ldir[0] = cdir[0] * As * mSpeedOfLight;
    ldir[1] = cdir[1] / As;
    ldir[2] = cdir[2] * r;
    ldir[3] = cdir[3] * r * sin(theta);
}


/*! Test break condition.
 *
 *  \param pos    : pointer to position array.
 *  \return true  : radial position r < 0.0 or  r^2<=(1.0+eps)*rs^2.
 *  \return false : position is valid.
 */
bool MetricReissnerNordstrom::breakCondition(const double* pos) {
    if (pos[1] <= 0.0) {
        return true;
    }

    if (mDiskr > 0.0)
        if (pos[1] >= rm && pos[1] <= rp) {
            return true;
        }

    return false;
}

/*! Calculate right hand side of the geodesic equation in first order form.
 *
 *  \param  y[]   : pointer to position and direction coordinates.
 *  \param  dydx[] : pointer to right side of geodesic equation.
 */
bool MetricReissnerNordstrom::calcDerivs(const double y[], double dydx[]) {
    dydx[0] = y[4];
    dydx[1] = y[5];
    dydx[2] = y[6];
    dydx[3] = y[7];

    double r = y[1];
    double theta = y[2];

    double A = 1.0 - rs / r + mK * mQ * mQ / (r * r);
    double w = rs * r - 2.0 * mK * mQ * mQ;
    double rh3 = r * r * r;

    dydx[4] = -w / A / rh3 * y[4] * y[5];
    dydx[5] = -0.5 * mSpeedOfLight * mSpeedOfLight * A * w / rh3 * y[4] * y[4] + 0.5 * w / A / rh3 * y[5] * y[5] + r * A * (y[6] * y[6] + sin(theta) * sin(theta) * y[7] * y[7]);
    dydx[6] = -2.0 / r * y[5] * y[6] + sin(theta) * cos(theta) * y[7] * y[7];
    dydx[7] = -2.0 / r * y[5] * y[7] - 2.0 * cos(theta) / sin(theta) * y[6] * y[7];

    return true;
}


/*! Tests whether the constraint equation is fulfilled.
 *
 *  The constraint equation for lightlike and timelike geodesics reads:
 \verbatim
     sum = g_{\mu\nu} dot(x)^{\mu} dot(x)^{\nu} - kappa c^2 = 0.
 \endverbatim
 *  \param  y[]   : pointer to position and direction coordinates.
 *  \param  kappa : timelike (-1.0), lightlike (0.0).
 *  \return double : sum.
 */
double MetricReissnerNordstrom::testConstraint(const double y[], const double kappa) {
    double r     = y[1];
    double theta = y[2];
    double st    = sin(theta);
    double cm = 1.0 / mSpeedOfLight;

    // Scale the directions with the speed of light before doubling them !!
    double dt = y[4];
    double dr = y[5] * cm;
    double dth = y[6] * cm;
    double dph = y[7] * cm;

    double sum = -kappa;
    double A = 1.0 - rs / r + mK * mQ * mQ / (r * r);
    sum += -A * dt * dt + dr * dr / A + r * r * (dth * dth + st * st * dph * dph);
    return sum;
}

/*! Set parameter 'pName' to 'val'.
 *
 *  Set 'mass' parameter and adjust Schwarzschild radius  rs=2GM/c^2.
 *  'charge' represents the charge of the black hole.
 */
bool MetricReissnerNordstrom::setParam(const char* pName, double val) {
    Metric::setParam(pName, val);

    if (strcmp(pName, "mass") == 0) {
        mMass = val;
        rs = 2.0 * mGravConstant * mMass / (mSpeedOfLight * mSpeedOfLight);
    }
    else if (strcmp(pName,"charge") == 0) {
        mQ = val;
    }

    calcDiskr();
    calcCritical();
    calcExtremals();
    return true;
}


/*! Transform point p to embedding coordinates.
 *
 *  \param p  : point to be transformed.
 *  \param ep : reference to 'embedded' point.
 *  \return true : success.
 *  \return false : otherwise.
 */
bool MetricReissnerNordstrom::transToEmbedding(vec4 p, vec4 &ep) {
    vec4 cp;
    transToPseudoCart(p, cp);

    double r = p[1];
    double x = cp[1];
    double y = cp[2];
    double z;
    if (r >= rp) {
        calcEmbeddingZ(r, z);
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
bool MetricReissnerNordstrom::setEmbeddingParam(const char *name, double val) {
    Metric::setEmbeddingParam(name, val);

    if (strcmp(name, "emb_rmin") == 0) {
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
bool MetricReissnerNordstrom::testEmbeddingParams() {
    bool allOk = true;
    if (mEmb_rmin < rp) {
        mEmb_rmin = rp;
        allOk &= false;
    }
    if (mEmb_rmax < rp) {
        mEmb_rmax = rp;
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
//int MetricReissnerNordstrom::getEmbeddingVertices(std::vector<vec3> &verts,
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
//            if (r >= rp) {
//                calcEmbeddingZ(r, z);
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

/*! Effective potential.
 *  \param pos : initial position.
 *  \param cdir : initial four-direction.
 *  \param type : geodesic type.
 *  \param x : abscissa value.
 *  \param val : reference to effective potential value.
 *  \return true : effective potential exists at x.
 */
bool MetricReissnerNordstrom::effPotentialValue(const vec4 pos, const vec4 cdir , enum_geodesic_type type, const double x, double &val) {
    double kappa = 0.0;
    if (type == enum_geodesic_timelike) {
        kappa = -mSign;
    }

    if ((pos[1] > rm - eps && pos[1] < rp + eps) || (x > rm - eps && x < rp + eps)) {
        return false;
    }

    double h = pos[1] * pos[1] * cdir[3];
    val = 0.5 * (1.0 - rs / x + mK * mQ * mQ / (x * x)) * (h * h / (x * x) - kappa * mSpeedOfLight * mSpeedOfLight);
    return true;
}

/*! Total energy.
 *  \param pos : initial position.
 *  \param cdir : initial four-direction.
 *  \param x : abscissa value.
 *  \param val : reference to total energy value.
 *  \return true : effective potential exists at x.
 */
bool MetricReissnerNordstrom::totEnergy(const vec4 pos, const vec4 cdir, const double , double &val) {
    if (pos[1] > rm - eps && pos[1] < rp + eps) {
        return false;
    }

    // 1/2*k^2/c^2:
    double Arn = 1.0 - rs / pos[1] + mK * mQ * mQ / (pos[1] * pos[1]);
    val = 0.5 * Arn * Arn * mSpeedOfLight * mSpeedOfLight * cdir[0] * cdir[0];
    return true;
}

/*!
 *  \param units : type of physical constants.
 */
void MetricReissnerNordstrom::usePhysicalUnits(const enum_physical_constants  units) {
    Metric::usePhysicalUnits(units);
    rs = 2.0 * mGravConstant * mMass / (mSpeedOfLight * mSpeedOfLight);
    mK = mGravConstant / (mDielectricPerm * pow(mSpeedOfLight, 4.0));
    calcDiskr();
}

/*!
 *  \param speed_of_light : value for speed of light.
 *  \param grav_const : value for gravitational constant.
 *  \param diel_perm : value for dielectric permittivity.
 */
void MetricReissnerNordstrom::setUnits(const double speed_of_light, const double grav_const, const double diel_perm) {
    Metric::setUnits(speed_of_light, grav_const, diel_perm);
    rs = 2.0 * mGravConstant * mMass / (mSpeedOfLight * mSpeedOfLight);
    mK = mGravConstant / (mDielectricPerm * pow(mSpeedOfLight, 4.0));
    calcDiskr();
}

/*! Generate report.
 */
bool MetricReissnerNordstrom::report(const vec4 pos, const vec4 cdir, std::string &text) {
    std::stringstream ss;
    ss << "Report for Reissner-Nordstrom metric\n\tcoordinate : (t,r,theta,phi)\n";
    ss << "--------------------------------------------------------------------------------\n";
    ss << "  physical units ................................. yes\n\n";
    ss.precision(DEF_FIXED_REPORT_PRECISION);
    ss.setf(std::ios::fixed);

    double temp = mGravConstant / (mSpeedOfLight * mSpeedOfLight * sqrt(mK));
    ss << "  Schwarzschild black hole .................... charge/mass = 0" << std::endl;
    ss << "  black hole region ........................... charge/mass < " << temp << std::endl;
    ss << "  extreme black hole .......................... charge/mass = " << temp << std::endl;
    ss << "  naked singularity (2 P-orbits) .............. charge/mass < " << 1.5 * temp / sqrt(2.) << std::endl;
    ss << "  naked singularity (1 P-orbit) ............... charge/mass = " << 1.5 * temp / sqrt(2.) << std::endl;
    ss << "  naked sing. (1 last unstable circ. orbit) ... charge/mass < " << temp*sqrt(5.) / 2. << std::endl;
    ss << "  naked sing. (only stable circ. orbits) ...... charge/mass = " << temp*sqrt(5.) / 2. << std::endl;
    ss << "                                     and ...... charge/mass > " << temp*sqrt(5.) / 2. << std::endl;
    ss << "\n";

    ss << "  Schwarzschild radius ..................... r_s = " << rs << std::endl;
    ss << "  outer RN-Horizon ......................... r_p = " << rp << std::endl;
    ss << "  inner RN-Horizon ......................... r_m = " << rm << std::endl;
    ss << "  unstable outer Photon orbit ............. r_np = " << rNp << std::endl;
    ss << "  stable inner Photon orbit ............... r_nm = " << rNm << std::endl;
    ss << "\n";

    ss << "  innermost stable circular orbit ........... r_isco = ";

    double B = mK * mQ * mQ / (rs * rs);
    double radiusIsco;
    calcRadiusIsco(radiusIsco, B);
    ss << radiusIsco << std::endl;

    ss << "  appendant velocity .......................... beta = ";

    double betaIsco;
    calcBetaCo(betaIsco, radiusIsco, B);
    ss << betaIsco << std::endl;

    if (B > 0.25 && B < 0.3125) {
        ss << "  second innermost stable circular orbit ... r_isco2 = ";

        double radiusIsco2;
        calcRadiusIsco2(radiusIsco2, B);
        ss << radiusIsco2 << std::endl;

        ss << "  appendant velocity .......................... beta = ";

        double betaIsco2;
        calcBetaCo(betaIsco2, radiusIsco2, B);
        ss << betaIsco2 << std::endl;
    }

    if (B > 0.28125 && B < 0.3125) {
        ss << "  innermost unstable circular orbit ......... r_iuco = ";

        double radiusIuco;
        calcRadiusIuco(radiusIuco, B);
        ss << radiusIuco << std::endl;

        ss << "  appendant velocity .......................... beta = ";

        double betaIuco;
        calcBetaCo(betaIuco, radiusIuco, B);
        ss << betaIuco << std::endl;
    }

    ss << "\n";
    double h   = pos[1] * pos[1] * cdir[3];
    double Arn = 1.0 - rs / pos[1] + mK * mQ * mQ / (pos[1] * pos[1]);
    double k   = Arn * mSpeedOfLight * mSpeedOfLight * cdir[0];
    ss << "  RN parameter ............................. A_RN = " << Arn << std::endl;
    ss << "  constant of motion .......................... k = " << k << std::endl;
    ss << "  constant of motion .......................... h = " << h << std::endl;
    ss << "  vel. for circ. orbit at initial pos. ..... beta = ";

    double betaCo;
    calcBetaCo(betaCo, pos[1], B);
    if (betaCo <= 1.) {
        ss << betaCo << std::endl;
    } else {
        ss << "not valid here" << std::endl;
    }

    ss << "  critical angle for null geodesic ...... ksiCrit = ";

    double ksicrit;
    if (calcKsiCrit(pos, ksicrit, B)) {
        ss << ksicrit*RAD_TO_DEG << std::endl;
    } else {
        ss << "not valid here.";
    }

    ss << "\n\nThe effective potential is only valid for geodesics in the theta=pi/2 hypersurface." << std::endl;
    text = ss.str();
    return true;
}

// ***************************** specific public methods ***************************
/*!  Calculate the critical angle in the black hole region (0<=B<=1/4) and in the case
 *   of a naked singularity for 1/4<B<=9/32. For B>9/32 there is no more Photon-orbit
 *   and ksicrit doesn't exist.
 *
 *   For an observer at distance r_i=pos[1] (initial radius), the black hole has an
 *   angular diameter, the critical angle. This angle additionally depends on the
 *   charge-mass-ratio. [ B=(rho*Q^2)/(r_s^2) ]
 *   For 1/4<B<=9/32 we have a naked singularity, but nevertheless there is an area
 *   around the singularity, you can't see. This angular diameter is described by
 *   ksicrit, too.
 *
 *  \param pos : current position
 *  \param ksicrit : reference to critical angle.
 *  \param B : B=(rho*Q^2)/(r_s^2)
 */
bool MetricReissnerNordstrom::calcKsiCrit(const vec4 pos, double &ksicrit, double B) {
    if (pos[1] < rp || B > 0.28125 || pos[1] < rNm) {
        return false;
    }

    double xi = rs / pos[1];
    double w1 = 3. + sqrt(9. - 32.*B);
    double numerator = w1 * w1 * w1 * w1 * xi * xi * (1. - xi + B * xi * xi);
    double denominator = 32.*(w1 - 8.*B);

    if (pos[1] >= rNp) {
        ksicrit = asin(sqrt(numerator / denominator));
    } else {
        ksicrit = M_PI - asin(sqrt(numerator / denominator));
    }

    return true;
}

/*! Calculate the innermost stable circular orbit

 *  \param radiusIsco : reference to the innermost stable circular orbit
 *  \param B : B=(rho*Q^2)/(r_s^2)
 */
bool MetricReissnerNordstrom::calcRadiusIsco(double &radiusIsco, double B) {
    if (B <= 0.25) {
        double w1 = B * sqrt(64.*B * B - 36.*B + 5.);
        double w2 = 8.*B * B - 9.*B + 2;

        double w3 = pow((w2 - w1), 1. / 3.);
        double w4 = pow((w2 + w1), 1. / 3.);
        double w5 = pow(0.5, 1. / 3.);

        radiusIsco = rs * (1 + w5 * (w3 + w4));
    } else {
        radiusIsco = 2 * B * rs;
    }

    return true;
}

/*! Calculate the second innermost stable circular orbit

 *  \param radiusIsco2 : reference to the second innermost stable circular orbit
 *  \param B : B=(rho*Q^2)/(r_s^2)
 */
bool MetricReissnerNordstrom::calcRadiusIsco2(double &radiusIsco2, double B) {
    if (B > 0.25 && B < 0.3125) {
        double w1 = sqrt(1. - 3.*B);

        double numerator = 2. - 9.*B + 8.*B * B;
        double denominator = 2.*w1 * w1 * w1;

        double w3 = acos(numerator / denominator);
        double w4 = cos(w3 / 3.);

        radiusIsco2 = rs * (1 + 2.*w1 * w4);
    }

    return true;
}

/*! Calculate the innermost unstable circular orbit

 *  \param radiusIuco : reference to the innermost unstable circular orbit
 *  \param B : B=(rho*Q^2)/(r_s^2)
 */
bool MetricReissnerNordstrom::calcRadiusIuco(double &radiusIuco, double B) {
    if (B > 0.28125 && B < 0.3125) {
        double w1 = sqrt(1. - 3.*B);

        double numerator = 2. - 9.*B + 8.*B * B;
        double denominator = 2.*w1 * w1 * w1;

        double w3 = acos(numerator / denominator);
        double w4 = cos((M_PI + w3) / 3.);

        radiusIuco = rs * (1 - 2.*w1 * w4);
    }

    return true;
}

/*! Calculate the velocity on a circular orbit

 *  \param betaCo : reference to the velocity on the innermost stable circular orbit
 *  \param radiusCo : radius of the circular orbit
 *  \param B : B=(rho*Q^2)/(r_s^2)
 */
bool MetricReissnerNordstrom::calcBetaCo(double &betaCo, double radiusCo, double B) {
    double xi = rs / radiusCo;

    double numerator = xi - 2.*B * xi * xi;
    double denominator = 2.*(1 - xi + B * xi * xi);

    betaCo = sqrt(numerator / denominator);

    return true;
}

// ********************************* protected methods *****************************
/*!
 */
void MetricReissnerNordstrom::setStandardValues() {
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

/*!
 */
void MetricReissnerNordstrom::calcDiskr() {
    mDiskr = 1.0 - 4.0 * mK * mQ * mQ / (rs * rs);
}

/*! Calculate critical points.
 */
void MetricReissnerNordstrom::calcCritical() {
    if (mDiskr < 0.0) {
        rp = rm = 0.0;
    } else {
        double w = sqrt(mDiskr);
        rp = 0.5 * rs * (1.0 + w);
        rm = 0.5 * rs * (1.0 - w);
    }
}

/*! Calculate extremal points for null geodesics.
 */
void MetricReissnerNordstrom::calcExtremals() {
    double w2 = 1.0 - 32.0 / 9.0 * mK * mQ * mQ / (rs * rs);
    if (w2 < 0.0) {
        rNp = rNm = 0.0;
    } else {
        rNp = 0.75 * rs * (1.0 + sqrt(w2));
        rNm = 0.75 * rs * (1.0 - sqrt(w2));
    }
}

/*! Calculate embedding function
 *  \param r : radial coordinate.
 *  \param z : reference to z-value.
 */
bool MetricReissnerNordstrom::calcEmbeddingZ(const double r, double &z) {
    struct_reissner_params par = {mMass, mK, mQ};
    F.params = &par;

    size_t limit = 1000;
    int    key   = GSL_INTEG_GAUSS15;
    double error;

    gsl_set_error_handler_off();
    gsl_integration_qag(&F, rp, r, 0.0, 1e-6, limit, key, w, &z, &error);
    return true;
}

} // end namespace m4d
