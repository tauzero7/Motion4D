// -------------------------------------------------------------------------------
/*
   m4dMetricSchwarzschildCart.cpp

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


using namespace std;

#include "m4dMetricSchwarzschildCartNew.h"

namespace m4d {

#define eps 1.0e-6


/*! Standard constructor for the Schwarzschild metric.
 *
 * \param  mass : mass of the black hole.
 */
MetricSchwarzschildCartNew::MetricSchwarzschildCartNew(double mass) {
    mMetricName  = "SchwarzschildCartNew";
    setCoordType(enum_coordinate_cartesian);

    mPhysicalUnits = enum_physical_constants_geom;
    mSpeedOfLight = 1.0;
    mGravConstant = 1.0;

    addParam("mass");
    setParam("mass", mass);

    rs = 2.0 * mass;

    for(int i = 0; i<4; ++i){
        for(int j = 0; j<4; ++j){
            for(int k = 0; k<4; ++k){
                christoffel[i][j][k]  = 0.0;
            }
        }
    }

    setStandardValues();
}

MetricSchwarzschildCartNew::~MetricSchwarzschildCartNew() {
}


// *********************************** public methods ******************************

/*! Calculate the contravariant metric components at position 'pos'.
 *
 *  \param pos : pointer to position.
 */
bool MetricSchwarzschildCartNew::calculateMetric(const double* pos) {
    double x = pos[1];
    double y = pos[2];
    double z = pos[3];

    double c = 1.0;

    double t1 = c * c;
    double t3 = x * x;
    double t4 = y * y;
    double t5 = z * z;
    double t6 = t3 + t4 + t5;
    double t7 = sqrt(t6);
    double t8 = 1 / t7;
    double t11 = 1 / t6;
    double t12 = t11 * t3;
    double t15 = 1 / (1.0 - rs * t8);
    double t17 = t11 * t4;
    double t18 = t11 * t5;
    double t20 = rs * t11;
    double t22 = 1 / (t7 - rs);
    double t23 = t22 * x;
    double t25 = t20 * t23 * y;
    double t27 = t20 * t23 * z;
    double t32 = t20 * t22 * y * z;

    g_compts[0][0] = -t1 + t1 * rs * t8;
    g_compts[0][1] = 0.0;
    g_compts[0][2] = 0.0;
    g_compts[0][3] = 0.0;
    g_compts[1][0] = 0.0;
    g_compts[1][1] = t12 * t15 + t17 + t18;
    g_compts[1][2] = t25;
    g_compts[1][3] = t27;
    g_compts[2][0] = 0.0;
    g_compts[2][1] = t25;
    g_compts[2][2] = t12 + t17 * t15 + t18;
    g_compts[2][3] = t32;
    g_compts[3][0] = 0.0;
    g_compts[3][1] = t27;
    g_compts[3][2] = t32;
    g_compts[3][3] = t12 + t17 + t18 * t15;

    return true;
}

/*! Calculate the Christoffel symbols of the second kind at position 'pos'.
 *
 *  \param pos : pointer to position.
 */
bool MetricSchwarzschildCartNew::calculateChristoffels(const double* pos) {
    double x = pos[1];
    double y = pos[2];
    double z = pos[3];
    //real x = yc[1];
    //real y = yc[2];
    //real z = yc[3];

    double x2 = x*x;
    double y2 = y*y;
    double z2 = z*z;

    double xy = x*y;
    double xz = x*z;
    double yz = y*z;

    // Warum das?
    //dydx[0] = yc[4];
    //dydx[1] = yc[5];
    //dydx[2] = yc[6];
    //dydx[3] = yc[7];

    double REAL_HALF = 0.5;
    double REAL_ONE = 1.0;
    double REAL_TWO = 2.0;


    //double rs = REAL_TWO*d_metricParams[0];

    double r = sqrt(x2+y2+z2);
    double rMinusrs = r - rs;
    double rMinusrs2 = rMinusrs*rMinusrs;
    double r2 = r*r;
    double r3 = r*r*r;
    double r4 = r2*r2;

    double OneMinusRsOverR = (REAL_ONE - rs/r);
    double OneMinusRsOverR2 = OneMinusRsOverR*OneMinusRsOverR;

    double rxDeviation2 = x2/OneMinusRsOverR + y2 + z2;
    double ryDeviation2 = x2 + y2/OneMinusRsOverR + z2;
    double rzDeviation2 = x2 + y2 + z2/OneMinusRsOverR;

    double DGttFac = -rs/r2;
    double rDx = x/r;
    double rDy = y/r;
    double rDz = z/r;

    double gtt = -OneMinusRsOverR;
    double gxx = rxDeviation2/r2;
    double gxy = rs/r2/rMinusrs*x*y;
    double gxz = rs/r2/rMinusrs*x*z;
    double gyz = rs/r2/rMinusrs*y*z;
    double gyy = ryDeviation2/r2;
    double gzz = rzDeviation2/r2;

    double gxySquare =  gxy*gxy;
    double gxzSquare =  gxz*gxz;
    double gyzSquare =  gyz*gyz;

    double giiFac1  = REAL_TWO / OneMinusRsOverR / r2;
    double giiFac1b = REAL_TWO / OneMinusRsOverR2 / r4;
    double giiFac2  = rs / OneMinusRsOverR2 / r4;
    double giiFac3  = REAL_TWO / r3;

    double gttDx = rDx*DGttFac;
    double gttDy = rDy*DGttFac;
    double gttDz = rDz*DGttFac;

    double gxxDx = x*giiFac1  - x2*rDx*giiFac2 - rxDeviation2*rDx*giiFac3;
    double gxxDy = y*giiFac1b - x2*rDy*giiFac2 - rxDeviation2*rDy*giiFac3;
    double gxxDz = z*giiFac1b - x2*rDz*giiFac2 - rxDeviation2*rDz*giiFac3;

    double gyyDx = x*giiFac1b - y2*rDx*giiFac2 - ryDeviation2*rDx*giiFac3;
    double gyyDy = y*giiFac1  - y2*rDy*giiFac2 - ryDeviation2*rDy*giiFac3;
    double gyyDz = z*giiFac1b - y2*rDz*giiFac2 - ryDeviation2*rDz*giiFac3;

    double gzzDx = x*giiFac1b - z2*rDx*giiFac2 - rzDeviation2*rDx*giiFac3;
    double gzzDy = y*giiFac1b - z2*rDy*giiFac2 - rzDeviation2*rDy*giiFac3;
    double gzzDz = z*giiFac1  - z2*rDz*giiFac2 - rzDeviation2*rDz*giiFac3;

    double gijFac1 = REAL_TWO*rs/r3/rMinusrs;
    double gijFac2 = rs / r2 / rMinusrs2;
    double gijFac3 = rs / r2 / rMinusrs;

    double gxyDx = -xy*rDx*(gijFac1 + gijFac2) + y*gijFac3;
    double gxyDy = -xy*rDy*(gijFac1 + gijFac2) + x*gijFac3;
    double gxyDz = -xy*rDz*(gijFac1 + gijFac2);

    double gxzDx = -xz*rDx*(gijFac1 + gijFac2) + z*gijFac3;
    double gxzDy = -xz*rDy*(gijFac1 + gijFac2);
    double gxzDz = -xz*rDz*(gijFac1 + gijFac2) + x*gijFac3;

    double gyzDx = -yz*rDx*(gijFac1 + gijFac2);
    double gyzDy = -yz*rDy*(gijFac1 + gijFac2) + z*gijFac3;
    double gyzDz = -yz*rDz*(gijFac1 + gijFac2) + y*gijFac3;

    double GDenominator = gxx*gyy*gzz+REAL_TWO*gxy*gxz*gyz-gxx*gyzSquare-gxySquare*gzz-gxzSquare*gyy;

    double G_x_t_t = REAL_HALF/gtt*gttDx;
    christoffel[1][0][0] = G_x_t_t;
    christoffel[0][1][0] = G_x_t_t;

    double G_y_t_t = REAL_HALF/gtt*gttDy;
    christoffel[2][0][0] = G_y_t_t;
    christoffel[0][2][0] = G_y_t_t;

    double G_z_t_t = REAL_HALF/gtt*gttDz;
    christoffel[3][0][0] = G_z_t_t;
    christoffel[0][3][0] = G_z_t_t;

    double G_t_t_x = -REAL_HALF*(
                gttDx*gyy*gzz-gttDx*gyzSquare+gttDy*gxz*gyz-gttDy*gxy*gzz+gttDz*gxy*gyz-gttDz*gxz*gyy)/GDenominator;
    christoffel[0][0][1] = G_t_t_x;

    double G_x_x_x = REAL_HALF*(
                gxxDx*gyy*gzz-gxxDx*gyzSquare+REAL_TWO*gxz*gyz*gxyDx-gxz*gyz*gxxDy-REAL_TWO*gxy*gzz*gxyDx+gxy*gzz*gxxDy
                +REAL_TWO*gxy*gyz*gxzDx-gxy*gyz*gxxDz-REAL_TWO*gxz*gyy*gxzDx+gxz*gyy*gxxDz)/GDenominator;
    christoffel[1][1][1] = G_x_x_x;

    double G_x_y_x = -REAL_HALF*(
                -gxxDy*gyy*gzz+gxxDy*gyzSquare-gyyDx*gxz*gyz+gyyDx*gxy*gzz-gxy*gyz*gxzDy-gxy*gyz*gyzDx+gxy*gyz*gxyDz
                +gxz*gyy*gxzDy+gxz*gyy*gyzDx-gxz*gyy*gxyDz)/GDenominator;
    christoffel[1][2][1] = G_x_y_x;
    christoffel[2][1][1] = G_x_y_x;

    double G_x_z_x = REAL_HALF*(
                gxxDz*gyy*gzz-gxxDz*gyzSquare+gxz*gyz*gxyDz+gxz*gyz*gyzDx-gxz*gyz*gxzDy-gxy*gzz*gxyDz-gxy*gzz*gyzDx
                +gxy*gzz*gxzDy+gzzDx*gxy*gyz-gzzDx*gxz*gyy)/GDenominator;
    christoffel[1][3][1] = G_x_z_x;
    christoffel[3][1][1] = G_x_z_x;

    double G_y_y_x = -REAL_HALF*(
                -REAL_TWO*gyy*gzz*gxyDy+gyy*gzz*gyyDx+REAL_TWO*gyzSquare*gxyDy-gyzSquare*gyyDx-gyyDy*gxz*gyz
                +gyyDy*gxy*gzz-REAL_TWO*gxy*gyz*gyzDy+gxy*gyz*gyyDz+REAL_TWO*gxz*gyy*gyzDy-gxz*gyy*gyyDz)/GDenominator;
    christoffel[2][2][1] = G_y_y_x;

    double G_y_z_x = REAL_HALF*(
                gyy*gzz*gxyDz+gyy*gzz*gxzDy-gyy*gzz*gyzDx-gyzSquare*gxyDz-gyzSquare*gxzDy+gyzSquare*gyzDx
                +gyyDz*gxz*gyz-gyyDz*gxy*gzz+gzzDy*gxy*gyz-gzzDy*gxz*gyy)/GDenominator;
    christoffel[2][3][1] = G_y_z_x;
    christoffel[3][2][1] = G_y_z_x;

    double G_z_z_x = -REAL_HALF*(
                -REAL_TWO*gyy*gzz*gxzDz+gyy*gzz*gzzDx+REAL_TWO*gyzSquare*gxzDz-gyzSquare*gzzDx-REAL_TWO*gxz*gyz*gyzDz
                +gxz*gyz*gzzDy+REAL_TWO*gxy*gzz*gyzDz-gxy*gzz*gzzDy-gzzDz*gxy*gyz+gzzDz*gxz*gyy)/GDenominator;
    christoffel[3][3][1] = G_z_z_x;

    double G_t_t_y = REAL_HALF*(
                -gttDx*gxz*gyz+gttDx*gxy*gzz-gttDy*gxx*gzz+gttDy*gxzSquare-gttDz*gxy*gxz+gttDz*gxx*gyz)/GDenominator;
    christoffel[0][0][2] = G_t_t_y;

    double G_x_x_y = -REAL_HALF*(
                -gxxDx*gxz*gyz+gxxDx*gxy*gzz-REAL_TWO*gxx*gzz*gxyDx+gxx*gzz*gxxDy+REAL_TWO*gxzSquare*gxyDx-gxzSquare*gxxDy
                -REAL_TWO*gxy*gxz*gxzDx+gxy*gxz*gxxDz+REAL_TWO*gxx*gyz*gxzDx-gxx*gyz*gxxDz)/GDenominator;
    christoffel[1][1][2] = G_x_x_y;

    double G_x_y_y = REAL_HALF*(
                gxz*gyz*gxxDy-gxy*gzz*gxxDy+gyyDx*gxx*gzz-gyyDx*gxzSquare+gxy*gxz*gxzDy+gxy*gxz*gyzDx-gxy*gxz*gxyDz
                -gxx*gyz*gxzDy-gxx*gyz*gyzDx+gxx*gyz*gxyDz)/GDenominator;
    christoffel[1][2][2] = G_x_y_y;
    christoffel[2][1][2] = G_x_y_y;

    double G_x_z_y = -REAL_HALF*(
                -gxxDz*gxz*gyz+gxxDz*gxy*gzz-gxx*gzz*gxyDz-gxx*gzz*gyzDx+gxx*gzz*gxzDy+gxzSquare*gxyDz+gxzSquare*gyzDx
                -gxzSquare*gxzDy-gzzDx*gxy*gxz+gzzDx*gxx*gyz)/GDenominator;
    christoffel[1][3][2] = G_x_z_y;
    christoffel[3][1][2] = G_x_z_y;

    double G_y_y_y = REAL_HALF*(
                REAL_TWO*gxz*gyz*gxyDy-gyyDx*gxz*gyz-REAL_TWO*gxy*gzz*gxyDy+gyyDx*gxy*gzz+gyyDy*gxx*gzz
                -gyyDy*gxzSquare+REAL_TWO*gxy*gxz*gyzDy-gxy*gxz*gyyDz-REAL_TWO*gxx*gyz*gyzDy+gxx*gyz*gyyDz)/GDenominator;
    christoffel[2][2][2] = G_y_y_y;

    double G_y_z_y = -REAL_HALF*(
                -gxz*gyz*gxyDz-gxz*gyz*gxzDy+gxz*gyz*gyzDx+gxy*gzz*gxyDz+gxy*gzz*gxzDy-gxy*gzz*gyzDx-gyyDz*gxx*gzz
                +gyyDz*gxzSquare-gzzDy*gxy*gxz+gzzDy*gxx*gyz)/GDenominator;
    christoffel[2][3][2] = G_y_z_y;
    christoffel[3][2][2] = G_y_z_y;

    double G_z_z_y = REAL_HALF*(
                REAL_TWO*gxz*gyz*gxzDz-gxz*gyz*gzzDx-REAL_TWO*gxy*gzz*gxzDz+gxy*gzz*gzzDx+REAL_TWO*gxx*gzz*gyzDz
                -gxx*gzz*gzzDy-REAL_TWO*gxzSquare*gyzDz+gxzSquare*gzzDy+gzzDz*gxy*gxz-gzzDz*gxx*gyz)/GDenominator;
    christoffel[3][3][2] = G_z_z_y;

    double G_t_t_z = REAL_HALF*(
                -gttDx*gxy*gyz+gttDx*gxz*gyy-gttDy*gxy*gxz+gttDy*gxx*gyz-gttDz*gxx*gyy+gttDz*gxySquare)/GDenominator;
    christoffel[0][0][3] = G_t_t_z;

    double G_x_x_z = -REAL_HALF*(
                -gxxDx*gxy*gyz+gxxDx*gxz*gyy-REAL_TWO*gxy*gxz*gxyDx+gxy*gxz*gxxDy+REAL_TWO*gxx*gyz*gxyDx-gxx*gyz*gxxDy
                -REAL_TWO*gxx*gyy*gxzDx+gxx*gyy*gxxDz+REAL_TWO*gxySquare*gxzDx-gxySquare*gxxDz)/GDenominator;
    christoffel[1][1][3] = G_x_x_z;

    double G_x_y_z = -REAL_HALF*(
                -gxxDy*gxy*gyz+gxxDy*gxz*gyy-gyyDx*gxy*gxz+gyyDx*gxx*gyz-gxx*gyy*gxzDy-gxx*gyy*gyzDx+gxx*gyy*gxyDz
                +gxySquare*gxzDy+gxySquare*gyzDx-gxySquare*gxyDz)/GDenominator;
    christoffel[1][2][3] = G_x_y_z;
    christoffel[2][1][3] = G_x_y_z;

    double G_x_z_z = REAL_HALF*(
                gxy*gyz*gxxDz-gxz*gyy*gxxDz+gxy*gxz*gxyDz+gxy*gxz*gyzDx-gxy*gxz*gxzDy-gxx*gyz*gxyDz-gxx*gyz*gyzDx
                +gxx*gyz*gxzDy+gzzDx*gxx*gyy-gzzDx*gxySquare)/GDenominator;
    christoffel[1][3][3] = G_x_z_z;
    christoffel[3][1][3] = G_x_z_z;

    double G_y_y_z = REAL_HALF*(
                REAL_TWO*gxy*gyz*gxyDy-gxy*gyz*gyyDx-REAL_TWO*gxz*gyy*gxyDy+gxz*gyy*gyyDx+gyyDy*gxy*gxz-gyyDy*gxx*gyz
                +REAL_TWO*gxx*gyy*gyzDy-gxx*gyy*gyyDz-REAL_TWO*gxySquare*gyzDy+gxySquare*gyyDz)/GDenominator;
    christoffel[2][2][3] = G_y_y_z;

    double G_y_z_z = -REAL_HALF*(
                -gxy*gyz*gxyDz-gxy*gyz*gxzDy+gxy*gyz*gyzDx+gxz*gyy*gxyDz+gxz*gyy*gxzDy-gxz*gyy*gyzDx-gxy*gxz*gyyDz
                +gxx*gyz*gyyDz-gzzDy*gxx*gyy+gzzDy*gxySquare)/GDenominator;
    christoffel[2][3][3] = G_y_z_z;
    christoffel[3][2][3] = G_y_z_z;


    double G_z_z_z = REAL_HALF*(
                REAL_TWO*gxy*gyz*gxzDz-gzzDx*gxy*gyz-REAL_TWO*gxz*gyy*gxzDz+gzzDx*gxz*gyy+REAL_TWO*gxy*gxz*gyzDz
                -gzzDy*gxy*gxz-REAL_TWO*gxx*gyz*gyzDz+gzzDy*gxx*gyz+gzzDz*gxx*gyy-gzzDz*gxySquare)/GDenominator;

    christoffel[3][3][3] = G_z_z_z;



   // coutChristoffel();
    //abort();
    return true;
}

void MetricSchwarzschildCartNew::coutChristoffel(){
    cout << "\n";
    for(int i = 0; i<4; ++i){
        for(int j = 0; j<4; ++j){
            for(int k = 0; k<4; ++k){
                if(abs(christoffel[i][j][k]) > 1e-14){
                    cout << "christoffel[" << i << "][" << j << "][" << k << "]=" << christoffel[i][j][k] << "\n";
                }
            }
        }
    }
}


/*! Transform local 4-direction to coordinate 4-direction.
 *
 *  \param  pos  :  pointer to position array.
 *  \param  ldir :  pointer to local direction array.
 *  \param  dir  :  pointer to calculated coordinate direction array.
 *  \param  type :  type of tetrad.
 */
void MetricSchwarzschildCartNew::localToCoord(const double* pos, const double* ldir, double* dir,
        enum_nat_tetrad_type) {
    double x = pos[1];
    double y = pos[2];
    double z = pos[3];

    double c = 1.0;
    double r = sqrt(x * x + y * y + z * z);
    double f = 1.0 - rs / r;

    calcLTcoeffs(pos);

    dir[0] = 1.0 / (c * sqrt(f)) * ldir[0];
    dir[1] = (A * ldir[1] + B * ldir[2] + D * ldir[3]);
    dir[2] = (C * ldir[2] + E * ldir[3]);
    dir[3] = F * ldir[3];
}

/*! Transform coordinate 4-direction to local 4-direction.
 *
 *  \param  pos  :  pointer to position array.
 *  \param  cdir :  pointer to coordinate direction.
 *  \param  ldir :  pointer to calculated local direction array.
 *  \param  type :  type of tetrad.
 */
void MetricSchwarzschildCartNew::coordToLocal(const double* pos, const double* cdir, double* ldir,
        enum_nat_tetrad_type) {
    double x = pos[1];
    double y = pos[2];
    double z = pos[3];

    double c = 1.0;
    double r = sqrt(x * x + y * y + z * z);
    double f = 1.0 - rs / r;

    calcLTcoeffs(pos);

    ldir[0] = c * sqrt(f) * cdir[0];
    ldir[1] = 1.0 / A * cdir[1] - B / (A * C) * cdir[2] + (B * E - C * D) / (A * C * F) * cdir[3];
    ldir[2] = 1.0 / C * cdir[2] - E / (C * F) * cdir[3];
    ldir[3] = 1.0 / F * cdir[3];
}


/*! Test break condition.
 *
 *  \param pos    : pointer to position array.
 *  \return true  : radial position r < 0.0 or  r^2<=(1.0+eps)*rs^2.
 *  \return false : position is valid.
 */
bool MetricSchwarzschildCartNew::breakCondition(const double* pos) {
    bool br = false;
    double r2 = pos[1] * pos[1] + pos[2] * pos[2] + pos[3] * pos[3];
//if (r2 < 0.01) {
    //  std::cerr << pos[1] << " " << pos[2] << " " << pos[3] << " " << r2 << " " << rs*rs << std::endl;
//}
    if (r2 <= (1.0 + eps)*rs * rs) {
        br = true;
    }
    return br;
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
double MetricSchwarzschildCartNew::testConstraint(const double y[], const double kappa) {
    double xx = y[1];
    double yy = y[2];
    double zz = y[3];
    double r = sqrt(xx * xx + yy * yy + zz * zz);
    double c = 1.0;

    double edrq = 1.0 / (r * r);
    double f    = 1.0 - rs / r;

    double sum = -kappa;

    sum += - c * c * f * y[4] * y[4]
           + (xx * xx / f + yy * yy + zz * zz) * edrq * y[5] * y[5]
           + (xx * xx + yy * yy / f + zz * zz) * edrq * y[6] * y[6]
           + (xx * xx + yy * yy + zz * zz / f) * edrq * y[7] * y[7]
           + 2.0 * rs / (r * r * (r - rs)) * (xx * yy * y[5] * y[6] + xx * zz * y[5] * y[7] + yy * zz * y[6] * y[7]);

    return sum;
}

/*! Set parameter 'pName' to 'val'.
 *
 *  Set 'mass' parameter and adjust Schwarzschild radius  rs=2GM/c^2.
 */
bool MetricSchwarzschildCartNew::setParam(const char* pName, double val) {
    if (Metric::setParam(pName, val)) {
        rs = 2.0 * mGravConstant * val / (mSpeedOfLight * mSpeedOfLight);
        return true;
    }
    return false;
}

/*! Generate report.
 */
bool MetricSchwarzschildCartNew::report(const vec4 , const vec4 , std::string &text) {
    std::stringstream ss;
    ss << "Report for Schwarzschild metric\n\tcoordinate : (t,x,y,z)\n";
    ss << "---------------------------------------------------------------\n";
    ss << "  physical units ................................. no\n";
    ss.precision(DEF_FIXED_REPORT_PRECISION);
    ss.setf(std::ios::fixed);
    ss << "  Schwarzschild radius ........... r_s = 2GM/c^2 = " << rs << std::endl;
    ss << "  Photon orbit ................... r_ph = 3/2*rs = " << 1.5 * rs << std::endl;
    ss << "  innermost stable circular orbit  r_isco = 3r_s = " << 3.0 * rs << std::endl;

    text = ss.str();
    return true;
}

// ********************************* protected methods *****************************
/*!
 */
void MetricSchwarzschildCartNew::setStandardValues() {
    mInitPos[0] = 0.0;
    mInitPos[1] = 3.0 * rs;
    mInitPos[2] = 0.0;
    mInitPos[3] = 0.0;
    mInitDir[0] = 1.0;
    mInitDir[1] = 0.0;
    mInitDir[2] = 0.0;

    mCoordNames[0] = std::string("t");
    mCoordNames[1] = std::string("x");
    mCoordNames[2] = std::string("y");
    mCoordNames[3] = std::string("z");
}

/*! Calculate local tetrad coefficients at position 'pos'.
 *
 *  \param pos : pointer to position.
 */
void MetricSchwarzschildCartNew::calcLTcoeffs(const double* pos) {
    calculateMetric(pos);

    double gxx = g_compts[1][1];
    double gyy = g_compts[2][2];
    double gzz = g_compts[3][3];

    double gxy = g_compts[1][2];
    double gxz = g_compts[1][3];
    double gyz = g_compts[2][3];

    double W = gxx * gyy * gzz - gxz * gxz * gyy + 2.0 * gxz * gxy * gyz - gxy * gxy * gzz - gxx * gyz * gyz;
    double N = gxx * gyy - gxy * gxy;

    double sW = sqrt(W);
    double sN = sqrt(N);

    A = 1.0 / sqrt(gxx);
    C = 1.0 / sqrt(-gxy * gxy / gxx + gyy);
    B = -gxy / gxx * C;
    F = sN / sW;
    E = (gxz * gxy - gxx * gyz) / (sN * sW);
    D = (gxy * gyz - gxz * gyy) / (sN * sW);
}

} // end namespace m4d
