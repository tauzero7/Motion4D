// -------------------------------------------------------------------------------
/*
    Metric.cpp

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

#include "m4dMetric.h"

namespace m4d {

#define eps 1.0e-6

/*! Standard constructor with no arguments.
 */
Metric::Metric() {
    mMetricName = "noName";
    mMetricCPPfilename = "unknown.cpp";

    mPhysicalUnits = enum_physical_constants_geom;
    mSpeedOfLight   = 1.0;
    mGravConstant   = 1.0;
    mDielectricPerm = 1.0;

    // Note that mSign=1.0 for sign(g)=+2 and mSign=-1.0 for sign(g)=-2 !
    mSign = 1.0;

    //init all to 0.0
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            g_compts[i][j] = 0.0;
            for (int k = 0; k < 4; k++) {
                christoffel[i][j][k] = 0.0;
                for (int l = 0; l < 4; l++) {
                    chrisD[i][j][k][l] = 0.0;
                }
            }
        }
    }
    //init minkowski
    g_compts[0][0] = -1.0;
    g_compts[1][1] = g_compts[2][2] = g_compts[3][3] = 1.0;

    mNumParam = 0;
    if (!mLocTeds.empty()) {
        mLocTeds.clear();
    }
    mLocTeds.push_back(enum_nat_tetrad_default);

    if (!mDrawTypes.empty()) {
        mDrawTypes.clear();
    }
    mDrawTypes.push_back(enum_draw_pseudocart);
    mDrawTypes.push_back(enum_draw_coordinates);

    tidalMatrix.setNull();
    coordTidalMatrix.setNull();

    setStandardValues();

    if (!mEmbParam.empty()) {
        mEmbParam.clear();
    }
    mHaveEmbedding    = false;

    mHaveEffPotential = false;

    mEnableResize = false;
}

/*! Standard destructor.
 */
Metric::~Metric() {
    if (!mEmbParam.empty()) {
        mEmbParam.clear();
    }
    if (!mParam.empty()) {
        mParam.clear();
    }
    if (!mLocTeds.empty()) {
        mLocTeds.clear();
    }
    if (!mDrawTypes.empty()) {
        mDrawTypes.clear();
    }
}


// *********************************** public methods ******************************
/*! Get the name of the metric.
 *
 */
std::string Metric::getMetricName() {
    return mMetricName;
}

/*! Get the filename of the metric.
 *
 */
std::string Metric::getMetricCPPfilename() {
    return mMetricCPPfilename;
}

/*! Get the type of the coordinates.
 *
 *  \return coordinate type of the metric.
 *  \sa enum_coordinate_type
 */
enum_coordinate_type
Metric::getCoordType() {
    return mCoordType;
}

/*! Get the name of the coordinate 'num'.
 *
 *  \param  num  :  number of coordinate (0,1,2,3).
 */
std::string Metric::getCoordName(int num) {
    if (num >= 0 && num <= 3) {
        return mCoordNames[num];
    }

    return std::string();
}

/*!  Get all coordinate names.
 *
 *  \param names[] : reference to coordinate names.
 */
void Metric::getCoordNames(std::string names[4]) {
    for (int i = 0; i < 4; i++) {
        names[i] = mCoordNames[i];
    }
}

/*! Return the signum of the signature of the metric.
 *  \return +1.0 : sign(g)=+2
 *  \return -1.0 : sign(g)=-2
 */
double Metric::sign() {
    return mSign;
}

/*! Calculate the contravariant metric components g_ab at position 'pos'.
 *
 *  \param  pos  :  coordinate position where the metric coefficients have to be evaluated.
 *  \return true : successfull
 */
bool Metric::calculateMetric(const vec4 pos) {
    double p[4] = { pos[0], pos[1], pos[2], pos[3] };
    return calculateMetric(p);
}

/*! Calculate the Christoffel symbols of the second kind at position 'pos'.
 *
 *  \param  pos  :  coordinate position where the Christoffels have to be evaluated.
 *  \return true :  successfull
 */
bool Metric::calculateChristoffels(const vec4 pos) {
    double p[4] = { pos[0], pos[1], pos[2], pos[3] };
    return calculateChristoffels(p);
}

/*! Calculate partial derivatives of the Christoffel symbols.
 * \param pos : pointer to position.
 */
bool Metric::calculateChrisD(const double*) {
    fprintf(stderr, "Metric::calculateChrisD ( const double* pos ) ... not implemented yet\n");
    return false;
}

/*! Calculate partial derivatives of the Christoffel symbols.
 * \param pos : coordinate position.
 */
bool Metric::calculateChrisD(const vec4 pos) {
    double p[4] = { pos[0], pos[1], pos[2], pos[3] };
    return calculateChrisD(p);
}

/*! Calculate Riemann tensor R^a_bcd
 * \param pos : pointer to coordinate position where the Riemann tensor have to be evaluated.
 * \return true : successfull
 */
bool Metric::calculateRiemann(const double*) {
    return false;
}

/*! Calculate Riemann tensor R^a_bcd
 * \param pos : coordinate position where the Riemann tensor have to be evaluated.
 * \return true : successfull
 */
bool Metric::calculateRiemann(const vec4 pos) {
    double p[4] = { pos[0], pos[1], pos[2], pos[3] };
    return calculateRiemann(p);
}

/*! Calculate the Weyl tensor at position 'pos'.
 *
 *  \param  pos  :  pointer to coordinate position where the Weyl tensor have to be evaluated.
 *  \return true :  successfull
 */
bool Metric::calculateWeyl(const double*) {
    return false;
}

/*! Calculate the Weyl tensor at position 'pos'.
 *
 *  \param  pos  :  coordinate position where the Weyl tensor have to be evaluated.
 *  \return true :  successfull
 */
bool Metric::calculateWeyl(const vec4 pos) {
    double p[4] = { pos[0], pos[1], pos[2], pos[3] };
    return calculateWeyl(p);
}

/*! Calculate the Ricci at position 'pos'.
 *
 *  \param  pos  :  pointer to coordinate position where the Ricci tensor have to be evaluated.
 *  \return true :  successfull
 */
bool Metric::calculateRicci(const double*) {
    return false;
}

/*! Calculate the Ricci at position 'pos'.
 *
 *  \param  pos  :  coordinate position where the Ricci tensor have to be evaluated.
 *  \return true :  successfull
 */
bool Metric::calculateRicci(const vec4 pos) {
    double p[4] = { pos[0], pos[1], pos[2], pos[3] };
    return calculateRicci(p);
}

/*! Calculate the Ricci rotation coefficients at position 'pos'.
 *
 *  \param  pos  :  pointer to coordinate position where the Ricci rotation coefficients have to be evaluated.
 *  \return true :  successfull
 */
bool Metric::calculateRicRotCoeffs(const double*) {
    return false;
}

/*! Calculate the Ricci rotation coefficients at position 'pos'.
 *
 *  \param  pos  :  coordinate position where the Ricci rotation coefficients have to be evaluated.
 *  \return true :  successfull
 */
bool Metric::calculateRicRotCoeffs(const vec4 pos) {
    double p[4] = { pos[0], pos[1], pos[2], pos[3] };
    return calculateRicRotCoeffs(p);
}

/*! Calculate the contractions of the Ricci rotation coefficients at position 'pos'.
 *
 *  \param  pos  :  pointer to coordinate position where the contractions of the Ricci rotation coefficients have to be evaluated.
 *  \return true :  successfull
 */
bool Metric::calculateContrRRC(const double*) {
    return false;
}

/*! Calculate the contractions of the Ricci rotation coefficients at position 'pos'.
 *
 *  \param  pos  :  coordinate position where the contractions of the Ricci rotation coefficients have to be evaluated.
 *  \return true :  successfull
 */
bool Metric::calculateContrRRC(const vec4 pos) {
    double p[4] = { pos[0], pos[1], pos[2], pos[3] };
    return calculateContrRRC(p);
}

/*!  Return the contravariant metric coefficient g_{mu nu}.
 *
 *   The coefficients have to be calculated with \sa calculateMetric().
 *   \param  mu : component index.
 *   \param  nu : component index.
 */
double Metric::getMetricCoeff(const int mu, const int nu) {
    return g_compts[mu][nu];
}

/*! Return the Christoffel symbol Chr_{mu nu}^{kappa} of the second kind.
 *
 *  The Christoffel symbols have to be calculated with \sa calculateChristoffels().
 *  \param  mu : lower index of Christoffel symbol.
 *  \param  nu : lower index of Christoffel symbol.
 *  \param  kappa : upper index of Christoffel symbol.
 */
double Metric::getChristoffel(const int mu, const int nu, const int kappa) {
    return christoffel[mu][nu][kappa];
}

/*!   Return the first derivative of the Christoffel symbol.
 *
 *  \param  mu : lower index of Christoffel symbol.
 *  \param  nu : lower index of Christoffel symbol.
 *  \param  kappa : upper index of Christoffel symbol.
 *  \param  tau : partial derivative index of Christoffel symbol.
 */
double Metric::getChrisD(const int mu, const int nu, const int kappa, const int tau) {
    return chrisD[mu][nu][kappa][tau];
}


/*!   Return a Riemann tensor component.
 *
 *  \param  mu :  first index.
 *  \param  nu :  second index.
 *  \param  rho : third index.
 *  \param  sigma :  fourth index.
 */
double Metric::getRiemCoeff(const int mu, const int nu, const int rho, const int sigma) {
    return riem[mu][nu][rho][sigma];
}

/*!   Return a Weyl tensor component.
 *
 *  \param  mu :  first index.
 *  \param  nu :  second index.
 *  \param  rho : third index.
 *  \param  sigma :  fourth index.
 */
double Metric::getWeylCoeff(const int mu, const int nu, const int rho, const int sigma) {
    return weyl[mu][nu][rho][sigma];
}

/*!   Return a Ricci tensor component.
 *
 *  \param  mu :  first index.
 *  \param  nu :  second index.
 */
double Metric::getRicciCoeff(const int mu, const int nu) {
    return ric[mu][nu];
}

/*!   Return a Ricci rotation coefficient.
 *
 *  \param  i : first index.
 *  \param  j : second index.
 *  \param  k : third index.
 */
double Metric::getRicRotCoeff(const int i, const int j, const int k) {
    return rrc[i][j][k];
}

/*!   Return a contracted Ricci rotation coefficient.
 *
 *  \param  j :  index.
 */
double Metric::getContrRRCCoeff(const int j) {
    return crrc[j];
}

/*! Get natural local tetrad.
 *
 *  \param  pos  :  position.
 *  \param  e0   :  reference to base vector 0.
 *  \param  e1   :  reference to base vector 1.
 *  \param  e2   :  reference to base vector 2.
 *  \param  e3   :  reference to base vector 3.
 *  \param  type :  type of local tetrad.
 *  \sa enum_nat_tetrad_type
 */
void Metric::getNatTetrad(const vec4 pos, vec4 &e0, vec4 &e1, vec4 &e2, vec4 &e3,
                          enum_nat_tetrad_type  type) {
    vec4 ldir0 = vec4(1.0, 0.0, 0.0, 0.0);
    vec4 ldir1 = vec4(0.0, 1.0, 0.0, 0.0);
    vec4 ldir2 = vec4(0.0, 0.0, 1.0, 0.0);
    vec4 ldir3 = vec4(0.0, 0.0, 0.0, 1.0);

    localToCoord(pos, ldir0, e0, type);
    localToCoord(pos, ldir1, e1, type);
    localToCoord(pos, ldir2, e2, type);
    localToCoord(pos, ldir3, e3, type);
}

void Metric::getNatDualTetrad(const vec4 pos, vec4 &t0, vec4 &t1, vec4 &t2, vec4 &t3,
                              enum_nat_tetrad_type  type) {
    vec4 cdir0 = vec4(1.0, 0.0, 0.0, 0.0);
    vec4 cdir1 = vec4(0.0, 1.0, 0.0, 0.0);
    vec4 cdir2 = vec4(0.0, 0.0, 1.0, 0.0);
    vec4 cdir3 = vec4(0.0, 0.0, 0.0, 1.0);

    coordToLocal(pos, cdir0, t0, type);
    coordToLocal(pos, cdir1, t1, type);
    coordToLocal(pos, cdir2, t2, type);
    coordToLocal(pos, cdir3, t3, type);
}


/*! Transform local 4-direction to coordinate 4-direction.
 *
 *  \param  pos  :  position.
 *  \param  ldir :  local direction.
 *  \param  cdir :  reference to coordinate direction.
 *  \param  type :  type of local tetrad.
 *  \sa enum_nat_tetrad_type
 */
void Metric::localToCoord(const vec4 pos, const vec4 ldir, vec4 &cdir,
                          enum_nat_tetrad_type type) {
    localToCoord(pos.data(), ldir.data(), cdir.data(), type);
}


/*! Transform coordinate 4-direction to local 4-direction.
 *
 *  \param  pos  :  position.
 *  \param  cdir :  coordinate direction.
 *  \param  ldir :  reference to local direction.
 *  \param  type :  type of local tetrad.
 *  \sa enum_nat_tetrad_type
 */
void Metric::coordToLocal(const vec4 pos, const vec4 cdir, vec4 &ldir,
                          enum_nat_tetrad_type  type) {
    coordToLocal(pos.data(), cdir.data(), ldir.data(), type);
}

/*! Test break condition.
 *
 *  \param pos :  position.
 */
bool Metric::breakCondition(const vec4 pos) {
    return breakCondition(pos.data());
}

/*! Calculate right hand side of the geodesic equation in first order form.
 *
 *  \param  y[]    : pointer to position and direction coordinates.
 *  \param  dydx[] : pointer to right side of geodesic equation.
 */
bool Metric::calcDerivs(const double* y, double* dydx) {
    //register int mu,j,k,l;
    register int mu, k, l;

    //double ch;
    calculateChristoffels(y);

    for (mu = 0; mu < 4; mu++) {
        dydx[mu] = y[4 + mu];
        dydx[mu + 4] = 0.0;

        for (k = 0; k < 4; k++)
            for (l = 0; l < 4; l++) {
                // ch += christoffel[k][l][mu]*yn[4+k]*yn[4+4*j+l];
                dydx[mu + 4] -= getChristoffel(k, l, mu) * y[4 + k] * y[4 + l];
            }
    }
    return true;
}

/*! Calculate right hand side of the geodesic equation in first order form with parallel transport.
 *
 *  \param  y[]    : pointer to position and direction coordinates.
 *  \param  dydx[] : pointer to right side of parallel transport equation.
 */
bool Metric::calcDerivsPar(const double* y, double* dydx) {
    register int mu, j, k, l;

    double ch;
    calculateChristoffels(y);

    for (mu = 0; mu < 4; mu++) {
        dydx[mu]   = y[4 + mu];
        dydx[mu + 4] = 0.0;

        for (j = 0; j < 4; j++) {
            dydx[8 + 4 * j + mu] = 0.0;

            ch = 0.0;
            for (k = 0; k < 4; k++) {
                dydx[mu + 4] -= getChristoffel(j, k, mu) * y[4 + j] * y[4 + k];
                for (l = 0; l < 4; l++) {
                    ch += getChristoffel(k, l, mu) * y[4 + k] * y[8 + 4 * j + l];
                }
            }
            dydx[8 + 4 * j + mu] -= ch;
        }
    }
    return true;
}

/*! Calculate right hand side of parallel transport and Jocobi equation.
 */
bool Metric::calcDerivsSachsJacobi(const double* , double*) {
    return false;
}

/*! Calculate right hand side of the Fermi-Walker transport equation.
 *
 *  This method has to be overwritten by the metric child classes. Otherwise, the standard form of
 *  the Fermi-Walker transport equation will be used.
 *  \sa FermiWalker::calcDerivs().
 *  \param  a[]    : pointer to proper acceleration.
 *  \param  y[]    : pointer to position and direction coordinates.
 *  \param  dydx[] : pointer to right hand side of Fermi-Walker transport equation.
 *  \return false : always.
 */
bool Metric::calcDerivsFW(const double* , const double* , double*) {
    return false;
}



/*! Tests whether the constraint equation is fulfilled.
 *
 *  The constraint equation for lightlike and timelike geodesics reads:
 \verbatim
       sum = g_{\mu\nu} dot(x)^{\mu} dot(x)^{\nu} - kappa c^2 = 0.
 \endverbatim
 *  This cause problems because of the limited precision of double.
 *  \param  y[]   : pointer to position and direction coordinates.
 *  \param  kappa : timelike (-1.0), lightlike (0.0).
 *  \return double : sum.
 */
double Metric::testConstraint(const double y[], const double kappa) {
    calculateMetric(y);
    double sum = -mSign * kappa * mSpeedOfLight * mSpeedOfLight;
    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++) {
            sum += g_compts[i][j] * y[4 + i] * y[4 + j];
        }

    return sum;
}

/*! Resize the tangent vector.
 *
 * \return true: have resized.
 * \return false: have not resized.
 */
bool Metric::resize(double*, double , double) {
    return false;
}

/*!
 * \return  enum_physical_constants : type of physical constants.
 */
enum_physical_constants Metric::physicalUnits() {
    return mPhysicalUnits;
}

/*!
 *  \param units : type of physical constants.
 */
void Metric::usePhysicalUnits(const enum_physical_constants  units) {
    mPhysicalUnits = units;

    if (mPhysicalUnits == enum_physical_constants_real) {
        mSpeedOfLight   = M4D_SPEED_OF_LIGHT;
        mGravConstant   = M4D_GRAV_CONST;
        mDielectricPerm = M4D_DIELECTRIC_PERM;
    } else if (mPhysicalUnits == enum_physical_constants_geom) {
        mSpeedOfLight   = 1.0;
        mGravConstant   = 1.0;
        mDielectricPerm = 1.0;
    }
}

/*!
 *  \param speed_of_light : value for speed of light.
 *  \param grav_const : value for gravitational constant.
 *  \param diel_perm : value for dielectric permittivity.
 */
void Metric::setUnits(const double speed_of_light, const double grav_const, const double diel_perm) {
    mPhysicalUnits = enum_physical_constants_proper;
    mSpeedOfLight   = speed_of_light;
    mGravConstant   = grav_const;
    mDielectricPerm = diel_perm;
}

/*! Return speed of light.
 */
double Metric::speed_of_light() {
    return mSpeedOfLight;
}

/*! Return gravitational constant.
 */
double Metric::grav_constant() {
    return mGravConstant;
}

/*! Return dielectric permittivity.
 */
double Metric::dielectric_perm() {
    return mDielectricPerm;
}


/*!  Generate local null direction out of a three direction.
 *
 *  The local null direction \f$\mathbf{k}=(k_0,k_1,k_2,k_3)^T\f$ follows from the three direction \f$\vec{d} = (d_1,d_2,d_3)^T\f$ via
 *   \f[ \mathbf{k} = \pm \mathbf{e}_{(0)} + \frac{1}{|\vec{d}|}\left(d_1\mathbf{e}_{(1)} + d_2\mathbf{e}_{(2)} + d_3\mathbf{e}_{(3)}\right).\f]
 *  Thus, \f$k_0 = \pm 1\f$, \f$k_1 = d_1/|\vec{d}|,\ldots\f$.
 *
 *  \param  pm : time direction.
 *  \param  l3dir : pointer to local direction which will be normalized.
 *  \param  l4dir : pointer to coordinate direction.
 *  \return true : success.
 *  \return false : l3dir cannot be normalized.
 */
bool Metric::localNullDir(const double pm, const double* l3dir, double* l4dir) {
    double norm = (l3dir[0] * l3dir[0] + l3dir[1] * l3dir[1] + l3dir[3] * l3dir[1]);
    if (norm < eps) {
        fprintf(stderr, "Metric::localNullDir():  can't normalize l3dir!\n");
        return false;
    }

    norm = 1.0 / norm;
    l4dir[0] = M4D_SIGN(pm);
    l4dir[1] = l3dir[0] * norm;
    l4dir[2] = l3dir[1] * norm;
    l4dir[3] = l3dir[2] * norm;
    return true;
}

/*!  Generate local null direction out of a three direction.
 *
 *  \param  pm : time direction.
 *  \param  l3dir : pointer to local direction which will be normalized.
 *  \param  l4dir : pointer to coordinate direction.
 *  \return true : success.
 *  \return false : l3dir cannot be normalized.
 */
bool Metric::localNullDir(const double pm, const vec3 l3dir, vec4 &l4dir) {
    double a3dir[3] = {l3dir[0], l3dir[1], l3dir[2]};
    double a4dir[4];
    bool ok = localNullDir(pm, a3dir, a4dir);
    if (ok) {
        l4dir = vec4(a4dir[0], a4dir[1], a4dir[2], a4dir[3]);
    }
    return ok;
}

/*!  Generate local four-velocity out of local three velocity.
 *
 *  The local four-velocity \f$\mathbf{u}=(u_0,u_1,u_2,u_3)^T\f$ follows from the three velocity \f$\vec{\beta}=(\beta_1,\beta_2,\beta_3)^T\f$ via
 *   \f[ \mathbf{u} = \pm c\gamma\left(\mathbf{e}_{(0)} + \beta_1\mathbf{e}_{(1)} + \beta_2\mathbf{e}_{(2)} + \beta_3\mathbf{e}_{(3)}\right), \f]
 *  where \f$\gamma=\left(\beta_1^2+\beta_2^2+\beta_3^2\right)^{-1/2}\f$. Thus, \f$u_0=\pm c\gamma\f$, \f$u_1=\pm c\gamma\beta_1,\ldots\f$.
 *
 *  \param  pm : time direction.
 *  \param  l3vel : pointer to local three velocity.
 *  \param  l4vel : pointer to local four velocity.
 *  \return true : success.
 *  \return false : l3vel exceeds the speed of light.
 */
bool Metric::localFourVel(const double pm, const double* l3vel, double* l4vel) {
    double vel = sqrt(l3vel[0] * l3vel[0] + l3vel[1] * l3vel[1] + l3vel[2] * l3vel[2]);
    if (vel >= mSpeedOfLight) {
        fprintf(stderr, "Metric::localFourVel():  too fast!\n");
        return false;
    }

    double beta  = vel / mSpeedOfLight;
    double gamma = 1.0 / sqrt(1.0 - beta * beta);

    l4vel[0] = M4D_SIGN(pm) * mSpeedOfLight * gamma;
    l4vel[1] = l3vel[0] * gamma;
    l4vel[2] = l3vel[1] * gamma;
    l4vel[3] = l3vel[2] * gamma;
    return true;
}

/*!  Generate local four-velocity out of local three velocity.
 *
 *  \param  pm : time direction.
 *  \param  l3dir : pointer to local three velocity.
 *  \param  l4dir : pointer to local four velocity.
 *  \return true : success.
 *  \return false : l3vel exceeds the speed of light.
 */
bool Metric::localFourVel(const double pm, const vec3 l3dir, vec4 &l4dir) {
    double a3dir[3] = {l3dir[0], l3dir[1], l3dir[2]};
    double a4dir[4];
    bool ok = localFourVel(pm, a3dir, a4dir);
    if (ok) {
        l4dir = vec4(a4dir[0], a4dir[1], a4dir[2], a4dir[3]);
    }
    return ok;
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
bool Metric::calcProduct(const double* pos, const double* u, const double* v, double &prod, bool preCalcMetric) {
    prod = 0.0;
    if (breakCondition(pos)) {
        return false;
    }

    if (preCalcMetric) {
        calculateMetric(pos);
    }

    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++) {
            prod += g_compts[i][j] * u[i] * v[j];
        }

    return true;
}

/*! Calculate the scalar product between u and v:  g_{ab}u^a v^b.
 *
 *  \param pos  :  position.
 *  \param u    :  vector.
 *  \param v    :  vector.
 *  \param prod :  reference to scalar product.
 *  \param preCalcMetric : calculate metric coefficients before evaluating the scalar product.
 *  \return true  : success.
 *  \return false : position is not valid.
 */
bool Metric::calcProduct(const vec4 pos, const vec4 u, const vec4 v, double &prod, bool preCalcMetric) {
    prod = 0.0;
    if (breakCondition(pos)) {
        return false;
    }

    if (preCalcMetric) {
        calculateMetric(pos);
    }

    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++) {
            prod += g_compts[i][j] * u[i] * v[j];
        }

    return true;
}

/*! Set Type of coordinates
 *
 * \param coord_type : Type of coordinates
 */
void Metric::setCoordType(enum_coordinate_type coord_type) {
    mCoordType  = coord_type;
    for (unsigned int i = 0; i < 4; i++) {
        mScoordType[i].type = enum_scoord_linear;
        mScoordType[i].min  = -DBL_MAX;
        mScoordType[i].max  =  DBL_MAX;
        mScoordType[i].character = enum_cchar_spacelike;
    }
    mScoordType[0].character = enum_cchar_timelike;
    switch (mCoordType) {
    // case enum_coordinate_cartesian: {}
    case enum_coordinate_spherical: {
        mScoordType[2].type = enum_scoord_periodic;
        mScoordType[2].min = 0;
        mScoordType[2].max = 1.0 * M_PI;
        mScoordType[3].type = enum_scoord_periodic;
        mScoordType[3].min = 0;
        mScoordType[3].max = 2.0 * M_PI;
        break;
    }
    case enum_coordinate_cylinder: {
        mScoordType[2].type = enum_scoord_periodic;
        mScoordType[2].min = 0;
        mScoordType[2].max = 2.0 * M_PI;
        break;
    }
    default:
    {}
    }
}

/*! Calculate difference between two coordinates.
 *
 * \param coordNum : coordinate number (0<=coordNum<4)
 * \param p1 : coordinate of point 1
 * \param p2 : coordinate of point 2
 * \return difference
 */
double Metric::coordDiff(const unsigned int coordNum, const double p1, const double p2) {
    assert(coordNum < 4);
    double diff;

    if (mScoordType[coordNum].type == enum_scoord_periodic) {
        diff = fmod(p2 - p1, mScoordType[coordNum].max);

        if (diff > 0.5 * mScoordType[coordNum].max) {
            diff = diff - mScoordType[coordNum].max;
        } else if (diff < -0.5 * mScoordType[coordNum].max) {
            diff = diff + mScoordType[coordNum].max;
        }
    } else {
        diff = p2 - p1;
    }
    return diff;
}

/*! Calculate separate space- and time-distance, if possible. The metric is evaluated as point p1.
 * \param p1 : pointer to point 1
 * \param p2 : pointer to point 2
 * \param spaceDist : reference to spacelike distance.
 * \param timeDist  : reference to timelike distance.
 * \return true : success.
 */
bool Metric::calcSepDist(const vec4 p1, const vec4 p2, double &spaceDist, double &timeDist) {
    spaceDist = 0.0;
    timeDist  = 0.0;

    if (!calculateMetric(p1)) {
        return false;
    }

    for (unsigned int i = 0; i < 4; i++) {
        for (unsigned int j = 0; j < 4; j++) {
            if ((mScoordType[i].character == enum_cchar_spacelike) ||
                    (mScoordType[j].character == enum_cchar_spacelike)) {
                spaceDist += g_compts[i][j] * coordDiff(i, p1[i], p2[i]) * coordDiff(j, p1[j], p2[j]);
                //fprintf(stderr,"%f %f ... ",coordDiff(i,p1[i],p2[i]),coordDiff(j,p1[j],p2[j]));
            }

            if ((mScoordType[i].character == enum_cchar_timelike) ||
                    (mScoordType[j].character == enum_cchar_timelike)) {
                timeDist  -= g_compts[i][j] * coordDiff(i, p1[i], p2[i]) * coordDiff(j, p1[j], p2[j]);
            }
        }
    }

    spaceDist = sqrt(fabs(spaceDist));
    timeDist  = sqrt(fabs(timeDist));
    return true;
}

bool Metric::calcTidalMatrix(const double* pos, const double* e0, const double* e1, const double* e2, const double* e3) {
    return calcTidalMatrix(vec4(pos), vec4(e0), vec4(e1), vec4(e2), vec4(e3));
}

bool Metric::calcTidalMatrix(const vec4 pos, const vec4 e0, const vec4 e1, const vec4 e2, const vec4 e3) {
    tidalMatrix.setNull();
    if (!calculateRiemann(pos)) {
        return false;
    }

    mat4 ee;
    ee.setRow(0, e0);
    ee.setRow(1, e1);
    ee.setRow(2, e2);
    ee.setRow(3, e3);
    mat4 tt = mat4(ee);
    tt.invert();

    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            double sum = 0.0;
            for (int mu = 0; mu < 4; mu++) {
                for (int nu = 0; nu < 4; nu++) {
                    for (int rho = 0; rho < 4; rho++) {
                        for (int sigma = 0; sigma < 4; sigma++) {
                            sum += tt.getElem(mu, i) * riem[mu][nu][rho][sigma] * e0[nu] * e0[rho] * ee.getElem(j, sigma);
                        }
                    }
                }
            }
            tidalMatrix.setElem(i, j, sum * mSpeedOfLight * mSpeedOfLight);
        }
    }
    return true;
}

double Metric::getTidalMatrixCoeff(const int i, const int j) {
    return tidalMatrix.getElem(i, j);
}

void Metric::getTidalMatrix(mat4 &m) {
    m = tidalMatrix;
}


bool Metric::calcCoordTidalMatrix(const double* pos, const double* u) {
    return calcCoordTidalMatrix(vec4(pos),vec4(u));
}

bool Metric::calcCoordTidalMatrix(const vec4 pos, const vec4 u) {
    coordTidalMatrix.setNull();
    if (!calculateRiemann(pos)) {
        return false;
    }

    for(int mu = 0; mu < 4; mu++) {
        for (int sigma = 0; sigma < 4; sigma++) {
            double sum = 0.0;
            for (int nu = 0; nu < 4; nu++) {
                for (int rho = 0; rho < 4; rho++) {
                    sum += riem[mu][nu][rho][sigma] * u[nu] * u[rho];
                }
            }
            coordTidalMatrix.setElem(mu,sigma,sum);
        }
    }
    return true;
}

double Metric::getCoordTidalMatrixCoeff(const int i, const int j) {
    return coordTidalMatrix.getElem(i, j);
}

void Metric::getCoordTidalMatrix(mat4 &m) {
    m = coordTidalMatrix;
}

void Metric::calcFmu_nu(const double* ) {
    // has to be implemented in the corresponding metric
}

double Metric::getFmu_nu(const int mu, const int nu) {
    return fmu_nu[mu][nu];
}


/*! Add parameter 'pName'.
 *
 *  \param pName  : new parameter name.
 *  \param val    : new parameter value.
 *  \return true  : parameter was added.
 *  \return false : parameter already exists.
 */
bool Metric::addParam(std::string pName, double val) {
    lowCase(pName);

    mParamItr = mParam.find(pName);
    if (mParamItr == mParam.end()) {
        mParam.insert(std::pair<std::string, double>(pName, val));
        mNumParam++;
        return true;
    } else {
        fprintf(stderr, "Parameter %s already exists!\n", pName.c_str());
        return false;
    }
}


/*! Set parameter 'pName' to 'val'.
 *
 *  \param pName : parameter name.
 *  \param val   : parameter value.
 *  \return true  : parameter was set.
 *  \return false : parameter was not set.
 */
bool Metric::setParam(std::string pName, double val) {
    lowCase(pName);
    mParamItr = mParam.find(pName);
    if (mParamItr == mParam.end()) {
        fprintf(stderr, "Parameter %s do no exist!\n", pName.c_str());
        return false;
    } else {
        mParamItr->second = val;
        return true;
    }
}

/*! Get parameter 'pName'.
 *
 *  \param pName : parameter name.
 *  \param val   : parameter value.
 *  \return true  : parameter exists.
 *  \return false : parameter do not exist.
 */
bool Metric::getParam(std::string pName, double &val) {
    lowCase(pName);
    mParamItr = mParam.find(pName);
    if (mParamItr == mParam.end()) {
        fprintf(stderr, "Parameter %s do no exist!\n", pName.c_str());
        return false;
    } else {
        val = mParamItr->second;
        return true;
    }
}

/*! Get parameter number 'pNr'.
 *
 *  \param pNr   : parameter number.
 *  \param pName : name of parameter.
 *  \param val   : value of parameter.
 *  \return true  : parameter exists.
 *  \return false : parameter do not exist.
 */
bool Metric::getParam(int pNr, std::string& pName, double& val) {
    if (pNr >= 0 && pNr < mNumParam) {
        std::map<std::string, double>::iterator p_itr = mParam.begin();
        for (int i = 0; i < pNr; i++) {
            p_itr++;
        }
        pName = (*p_itr).first;
        val   = (*p_itr).second;

        return true;
    }
    return false;
}

bool Metric::setParam(int pNr, double val) {
    if (pNr >= 0 && pNr < mNumParam) {
        std::map<std::string, double>::iterator p_itr = mParam.begin();
        for (int i = 0; i < pNr; i++) {
            p_itr++;
        }
        std::string pName = p_itr->first;
        return setParam(pName, val);
    }
    return false;
}

/*! Get number of parameters.
 *
 * \return number of parameters.
 */
int Metric::getNumParams() {
    return mNumParam;
}

/*! Get list of parameter names.
 *
 *  \param  names : reference to vector of strings.
 */
void Metric::getParamNames(std::vector<std::string> &names) {
    if (!names.empty()) {
        names.clear();
    }

    mParamItr = mParam.begin();
    while (mParamItr != mParam.end()) {
        names.push_back(mParamItr->first);
        mParamItr++;
    }
}

int Metric::getParamNum(std::string name) {
    int count = 0;
    mParamItr = mParam.begin();
    while (mParamItr != mParam.end()) {
        if (mParamItr->first.compare(name) == 0) {
            return count;
        }
        mParamItr++;
        count++;
    }
    return -1;
}

/*! Get types of local tetrads defined for the metric.
 *
 *  All metrics have a default local tetrads which should be adapted to
 *  the symmetries and the coordinates.
 *  \param locted : reference to vector of types of local tetrads.
 *  \return int : number of local tetrad types defined for the metric.
 *  \sa enum_nat_tetrad_type
 */
int Metric::getLocTedTypes(std::vector<enum_nat_tetrad_type> &locted) {
    locted = mLocTeds;
    return (int)mLocTeds.size();
}

/*! Get type of local tetrad 'num'.
 *
 * \param  num :  number of local tetrad type.
 */
enum_nat_tetrad_type
Metric::getCurrLTtype(int num) {
    if (num < 0 || num >= (int)mLocTeds.size()) {
        return enum_nat_tetrad_default;
    }

    return mLocTeds[num];
}

/*! Get types of drawing methods defined for the metric.
 *
 *  \param drawTypes : reference to vector of draw types.
 *  \sa enum_draw_type.
 */
int Metric::getDrawTypes(std::vector<enum_draw_type> &drawTypes) {
    drawTypes = mDrawTypes;
    return (int)mDrawTypes.size();
}

/*! Get type of drawing method 'num'.
 *
 *  \param num : number of draw type.
 */
enum_draw_type Metric::getCurrDrawType(int num) {
    if (num < 0 || num >= (int)mDrawTypes.size()) {
        return enum_draw_pseudocart;
    }

    return mDrawTypes[num];
}

/*! Get type of drawing method 'num'.
 *
 *  \param name : name of draw type.
 */
enum_draw_type Metric::getCurrDrawType(std::string name) {
    int num = 0;
    while (num < (int)mDrawTypes.size()) {
        if (name == stl_draw_type[mDrawTypes[num]]) {
            return enum_draw_type(mDrawTypes[num]);
        }
        num++;
    }
    return enum_draw_pseudocart;
}

/*! Transform point p to pseudo-cartesian coordinates.
 *
 *  The transformation between proper metric coordinates and pseudo-cartesian
 *  coordinates can also be overloaded by each metric.
 *  \param  p  : point in proper metric coordinates.
 *  \param  cp : reference to transformed point.
 *  \return ID of chart
 */
int Metric::transToPseudoCart(vec4 p, vec4 &cp) {
    if (mCoordType == enum_coordinate_prolatespheroidal) {
        TransCoordinates::transCoordProlSphCart(p, cp, mprolSphfac);
    } else {
        TransCoordinates::toCartesianCoord(mCoordType, p, cp);
    }
    return 0;
}

/*! Transform point p to embedding coordinates.
 *
 *  The transformation between proper metric coordinates and embedding
 *  coordinates should be overloaded by the metric which do have an
 *  embedding.
 *  \param  p  : point in proper metric coordinates.
 *  \param  ep : reference to transformed point.
 *  \return true  : success.
 *  \return false : embedding not available.
 */
bool Metric::transToEmbedding(vec4 , vec4&) {
    return false;
}

/*! Transform point p to 2+1 coordinates.
 *
 *  \param  p  : point in proper metric coordinates.
 *  \param  cp : reference to transformed point.
 *  \return true : success.
 */
bool Metric::transToTwoPlusOne(vec4 p, vec4 &cp) {
    vec4 tp;
    TransCoordinates::toCartesianCoord(mCoordType, p, tp);
    cp = vec4(tp[0], tp[1], tp[2], tp[0]);
    return true;
}

/*! Transform point p to custom coordinates.
 *
 *  \param  p  : point in proper metric coordinates.
 *  \param  cp : reference to transformed point.
 *  \return true : success.
 */
bool Metric::transToCustom(vec4 p, vec4 &cp) {
    cp = p;
    return true;
}

/*! Transform from pseudo-Cartesian coordinates into proper metric coordinates.
 *
 * \param cp : point in pseudo-Cartesian coordinates
 * \param p  : reference to transformed point.
 * \return true : success.
 */
void Metric::transFromPseudoCart(vec4 cp, vec4 &p) {
    TransCoordinates::coordTransf(enum_coordinate_cartesian, cp, mCoordType, p);
}

/*! Add embedding parameters.
 *
 *  \param  name : new embedding parameter name.
 *  \param  val  : new embedding parameter value.
 *  \return true  : parameter was added.
 *  \return false : parameter already exists.
 */
bool Metric::addEmbeddingParam(std::string name, double val) {
    lowCase(name);

    mEmbParamItr = mEmbParam.find(name);
    if (mEmbParamItr == mEmbParam.end()) {
        mEmbParam.insert(std::pair<std::string, double>(name, val));
        return true;
    } else {
#if DEF_SHOW_EMB_WARN==1
        fprintf(stderr, "Embedding parameter %s already exists!\n", name.c_str());
#endif
        return false;
    }
}

/*! Set embedding parameters.
 *
 *  \param  name : embedding parameter name.
 *  \param  val  : embedding parameter value.
 *  \return true  : success.
 *  \return false : parameter not valid.
 */
bool Metric::setEmbeddingParam(std::string name, double val) {
    lowCase(name);
    mEmbParamItr = mEmbParam.find(name);
    if (mEmbParamItr == mEmbParam.end()) {
#if DEF_SHOW_EMB_WARN==1
        fprintf(stderr, "Embedding parameter %s do no exist!\n", name.c_str());
#endif
        return false;
    } else {
        mEmbParamItr->second = val;
        return true;
    }
}

/*! Get embedding parameters.
 *
 *  \param  name : embedding parameter name.
 *  \param  val  : reference to embedding parameter value.
 *  \return true  : success.
 *  \return false : parameter not valid.
 */
bool Metric::getEmbeddingParam(std::string name, double &val) {
    lowCase(name);
    mEmbParamItr = mEmbParam.find(name);
    if (mEmbParamItr == mEmbParam.end()) {
#if DEF_SHOW_EMB_WARN==1
        fprintf(stderr, "Embedding parameter %s do no exist!\n", name.c_str());
#endif
        return false;
    } else {
        val = mEmbParamItr->second;
        return true;
    }
}

/*! Get embedding parameter names.
 *
 *  \param names : reference to vector of strings.
 */
void Metric::getEmbeddingNames(std::vector<std::string> &names) {
    if (!names.empty()) {
        names.clear();
    }

    mEmbParamItr = mEmbParam.begin();
    while (mEmbParamItr != mEmbParam.end()) {
        names.push_back(mEmbParamItr->first);
        mEmbParamItr++;
    }
}

/*! Get all embedding parameters.
 * \param names : reference to parameter name list.
 * \param params : reference to parameter value list.
 * \return true : both lists have the same size.
 */
bool Metric::getAllEmbeddingParams(std::vector<std::string> &names, std::vector<double> &params) {
    if (!params.empty()) {
        params.clear();
    }
    getEmbeddingNames(names);

    double val;
    for (unsigned int i = 0; i < names.size(); i++) {
        if (getEmbeddingParam(names[i], val)) {
            params.push_back(val);
        }
    }

    if (params.size() == names.size()) {
        return true;
    }

    names.clear();
    params.clear();
    return false;
}

/*! Get embedding parameter map.
 *  \param params : reference to embedding parameter map.
 *  \return true: successfull.
 */
bool Metric::getEmbeddingMap(std::map<std::string, double> &params) {
    if (mEmbParam.size() == 0) {
        return false;
    }

    if (!params.empty()) {
        params.clear();
    }

    mEmbParamItr = mEmbParam.begin();
    while (mEmbParamItr != mEmbParam.end()) {
        params.insert(std::pair<std::string, double>(mEmbParamItr->first, mEmbParamItr->second));
        mEmbParamItr++;
    }
    return true;
}

/*! Generate vertices for the embedding diagram.
 *
 *  \param verts : reference to vector of vertices.
 *  \param indices : reference to vector of indices.
 *  \param numElems : number of elements in a strip.
 *  \param counter  : number of strips.
 *  \return int : 0
 */
int Metric::getEmbeddingVertices(std::vector<vec3>&,
                                 std::vector<int>&, unsigned int&, unsigned int&) {
    return 0;
}

/*! Have embedding diagram.
 */
bool Metric::haveEmbedding() {
    return mHaveEmbedding;
}

/*! Get effective potential value.
 */
bool Metric::effPotentialValue(const vec4 , const vec4 , enum_geodesic_type , const double , double &) {
    return false;
}

bool Metric::totEnergy(const vec4 , const vec4 , const double , double &) {
    return false;
}

/*! Have effective potential.
 */
bool Metric::haveEffPotential() {
    return mHaveEffPotential;
}

/*! Determine the velocity for a closed circular orbit if it exists.
 * \param r  Radial coordinate.
 * \param tedType  Velocity with respect to natural local tetrad type.
 * \return  beta = v/c
 */
double Metric::getCircularVelocity(const double , const enum_nat_tetrad_type) {
    fprintf(stderr, "Metric::getCircularVelocity() not available for this metric.\n");
    return 0.0;
}

vec4
Metric::getCircularFourVel(const vec4 , const enum_nat_tetrad_type) {
    fprintf(stderr, "Metric::getCircularFourVel() not available for this metric.\n");
    return vec4();
}

/*! Get standard initial position.
 *
 *  \return vec4 : initial position.
 */
vec4
Metric::getStandardInitialPosition() {
    return vec4(mInitPos[0], mInitPos[1], mInitPos[2], mInitPos[3]);
}

/*! Get standard initial direction with respect to the local tetrad.
 *
 *  \return vec3 : initial direction with respect to the local tetrad.
 */
vec3
Metric::getStandardInitialDirection() {
    return vec3(mInitDir[0], mInitDir[1], mInitDir[2]);
}

bool Metric::isResizeEnabled() {
    return mEnableResize;
}

/*!  Generate report.
 *
 *  \param  pos  : initial position.
 *  \param  cdir : initial coordinate direction.
 *  \param  text : reference to string.
 *  \return true  : success.
 *  \return false : no report available.
 */
bool Metric::report(const vec4 , const vec4 , std::string &text) {
    text = std::string();
    return false;
}


/*!
 *  \param fptr : pointer to FILE.
 */
void Metric::printF(FILE* fptr) {
    fprintf(fptr, "\nMetricName : %s\n", mMetricName.c_str());
    for (int i = 0; i < 13 + (int)mMetricName.length(); i++) {
        fprintf(fptr, "-");
    }
    fprintf(fptr, "\n   sign(g) : %s\n", (mSign > 0.0 ? "+2" : "-2"));
    fprintf(fptr, "\n   #params : %d\n", (int)mParam.size());
    mParamItr = mParam.begin();
    while (mParamItr != mParam.end()) {
        fprintf(fptr, " %12s : %7.4f\n", mParamItr->first.c_str(), mParamItr->second);
        mParamItr++;
    }
}

// ********************************* protected methods *****************************
/*!
 * \param  s : reference to string.
 */
void Metric::lowCase(std::string &s) {
    int l = (int)s.length();
    for (int i = 0; i < l; i++) {
        s[i] = tolower(s[i]);
    }
}

/*! Set coordinate names, initial coordinate position and initial local direction.
 *
 */
void Metric::setStandardValues() {
    mCoordNames[0] = std::string("x0");
    mCoordNames[1] = std::string("x1");
    mCoordNames[2] = std::string("x2");
    mCoordNames[3] = std::string("x3");

    for (unsigned int i = 0; i < 4; i++) {
        mScoordType[i].type = enum_scoord_linear;
        mScoordType[i].min  = -DBL_MAX;
        mScoordType[i].max  =  DBL_MAX;
        mScoordType[i].character = enum_cchar_spacelike;
    }
    mScoordType[0].character = enum_cchar_timelike;
}


/*! Contract Christoffel symbol with two vectors.
 *  \param y : full data
 *  \param v : first vector
 *  \param w : second vector
 *  \param z : contraction
 *  \param calc : calculate the Christoffels before using them (only is child classes).
 *
 * Note that the Christoffel symbols have to be calculated before this function can be evaluated!
 *
 *   \f[z^{\mu} = \Gamma_{\rho\sigma}^{\mu} v^{\rho}w^{\sigma} \f]
 */
void Metric::contrChrisVecVec(const double* , const double* v, const double* w, double* z, bool) {
    for (int mu = 0; mu < 4; mu++) {
        z[mu] = 0.0;
        for (int i = 0; i < 4; i++)
            for (int j = 0; j < 4; j++) {
                z[mu] += christoffel[i][j][mu] * v[i] * w[j];
            }
    }
}

/*! Contract partially derived Christoffel symbols with three vectors.
 *  \param y : full data
 *  \param u : first vector
 *  \param v : second vector
 *  \param w : third vector
 *  \param z : contraction
 *  \param calc : calculate the Christoffels before using them (only is child classes).
 *
 * Note that the partial derivatives of the Christoffel symbols have to be calculated before this function can be evaluated!
 *
 *   \f[z^{\mu} = \Gamma_{\rho\sigma;\lambda}^{\mu} u^{\rho}v^{\sigma}w^{\lambda} \f]
 */
void Metric::contrChrDVecVecVec(const double* , const double* u, const double* v, const double* w, double* z, bool) {
    for (int mu = 0; mu < 4; mu++) {
        z[mu] = 0.0;
        for (int i = 0; i < 4; i++)
            for (int j = 0; j < 4; j++)
                for (int l = 0; l < 4; l++) {
                    z[mu] += chrisD[i][j][mu][l] * u[i] * v[j] * w[l];
                }
    }
}

} // end namespace m4d

