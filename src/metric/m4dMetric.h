// --------------------------------------------------------------------------------
/*
    m4dMetric.h

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

/*!  \class  m4d::Metric
     \brief  Base class for each metric class.

             To determine the metric coefficients and the Christoffel
             symbols one can use e.g. Maple/GrTensorII or maxima:
             <ul>
               <li>\f$g_{\mu\nu}\f$  : grcalc(g(dn,dn));
               <li>\f$\Gamma_{\mu\nu}^{\kappa}\f$ : grcalc(Chr(dn,dn,up)) = grcalc(Chr2);
               <li>\f$\partial_{\lambda}\Gamma_{\mu\nu}^{\kappa}\f$ : grcalc(Chr(dn,dn,up,pdn));
             </ul>

             If physical constants (speed of light, Newton's gravitational
             constant,...) are NOT included in the metric coefficients and the
             Christoffel symbols, the attribute mPhysicalUnits has to be set to
             enum_physical_constants_notset. Otherwise, there are several units
             that can be chosen:

             <ul>
               <li>enum_physical_constants_geom:  G=c=1<br>
                   Distances are measured in units of \f$M\f$, where \f$1M \approx 4.926938\cdot 10^{-6}~ls\f$.
               <li>enum_physical_constants_real:<br>
                   Here, the speed of light as well as the gravitational constant are given in SI units: \f$c=299792458~m\f$ and \f$G=6.67428\cdot 10^{-11}~\frac{m^3}{kg\cdot s^2}\f$.
               <li>enum_physical_constants_proper:<br>
                   The numerical values of all constants can be chosen arbitrarily.
             </ul>

             calcDerivs(), calcDerivsPar(), and testConstraint() can be
             overloaded in the child classes to perform a faster calculation
             (mHaveRightSide must be set true!), otherwise they are calculated
             as in the motion class. For the Fermi-Walker transport, the
             calcDerivsFW() method can be overloaded in the child classes to
             perform a faster calculation.

             Each metric must have a natural local tetrad (default).

             Each metric must have a pseudo-cartesian representation.
             Additionally, the coordinates can be drawn with respect
             to each other.

             The signature of the metric must be adjusted by hand. If not otherwise specified,
             \f[ \textrm{sign}(\mathbf{g})=+2. \f]
             Then mSign=1.0, otherwise mSign=-1.0.
*/
// --------------------------------------------------------------------------------

#ifndef M4D_METRIC_H
#define M4D_METRIC_H

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstring>
#include <iostream>
#include <map>
#include <sstream>
#include <vector>

#include <m4dGlobalDefs.h>
#include <math/TransCoordinates.h>

#ifdef _WIN32
#ifndef __GNUC__
#pragma warning(disable : 4244)
#endif
#endif

#define DEF_SHOW_EMB_WARN 0

namespace m4d {

// ---------------------------------------------------
//    class definition:   Metric
// ---------------------------------------------------
class API_M4D_EXPORT Metric {
public:
    Metric();
    virtual ~Metric();

    // --------- public methods -----------
public:
    const char* getMetricName();
    const char* getMetricCPPfilename();
    enum_coordinate_type getCoordType();
    const char* getCoordName(int num);
    double sign();

    virtual bool calculateMetric(const double* pos) = 0;
    virtual bool calculateMetric(const vec4 pos);
    virtual bool calculateChristoffels(const double* pos) = 0;
    virtual bool calculateChristoffels(const vec4 pos);
    virtual bool calculateChrisD(const double* pos);
    virtual bool calculateChrisD(const vec4 pos);

    virtual bool calculateRiemann(const double* pos);
    virtual bool calculateRiemann(const vec4 pos);
    virtual bool calculateWeyl(const double* pos);
    virtual bool calculateWeyl(const vec4 pos);
    virtual bool calculateRicci(const double* pos);
    virtual bool calculateRicci(const vec4 pos);
    virtual bool calculateRicRotCoeffs(const double* pos);
    virtual bool calculateRicRotCoeffs(const vec4 pos);
    virtual bool calculateContrRRC(const double* pos);
    virtual bool calculateContrRRC(const vec4 pos);

    double getMetricCoeff(const int mu, const int nu);
    double getChristoffel(const int mu, const int nu, const int kappa);
    double getChrisD(const int mu, const int nu, const int kappa, const int tau);

    double getRiemCoeff(const int mu, const int nu, const int rho, const int sigma);
    double getWeylCoeff(const int mu, const int nu, const int rho, const int sigma);
    double getRicciCoeff(const int mu, const int nu);
    double getRicRotCoeff(const int i, const int j, const int k);
    double getContrRRCCoeff(const int j);

    virtual void getNatTetrad(const vec4 pos, vec4& e0, vec4& e1, vec4& e2, vec4& e3,
        enum_nat_tetrad_type type = enum_nat_tetrad_default);

    virtual void getNatDualTetrad(const vec4 pos, vec4& t0, vec4& t1, vec4& t2, vec4& t3,
        enum_nat_tetrad_type type = enum_nat_tetrad_default);

    virtual void localToCoord(const double* pos, const double* ldir, double* dir,
        enum_nat_tetrad_type type = enum_nat_tetrad_default)
        = 0;
    virtual void localToCoord(const vec4 pos, const vec4 ldir, vec4& cdir,
        enum_nat_tetrad_type type = enum_nat_tetrad_default);
    virtual void coordToLocal(const double* pos, const double* cdir, double* ldir,
        enum_nat_tetrad_type type = enum_nat_tetrad_default)
        = 0;
    virtual void coordToLocal(const vec4 pos, const vec4 cdir, vec4& ldir,
        enum_nat_tetrad_type type = enum_nat_tetrad_default);

    virtual bool breakCondition(const double* pos) = 0;
    virtual bool breakCondition(const vec4 pos);

    virtual bool calcDerivs(const double* y, double* dydx);
    virtual bool calcDerivsPar(const double* y, double* dydx);
    virtual bool calcDerivsSachsJacobi(const double* y, double* dydx);
    virtual bool calcDerivsFW(const double* a, const double* y, double* dydx);

    virtual double testConstraint(const double y[], const double kappa);
    virtual bool resize(double* y, double kappa, double factor = DEF_RESIZE_FACTOR);

    bool localNullDir(const double pm, const double* l3dir, double* l4dir);
    bool localNullDir(const double pm, const vec3 l3dir, vec4& l4dir);
    bool localFourVel(const double pm, const double* l3vel, double* l4vel);
    bool localFourVel(const double pm, const vec3 l3dir, vec4& l4dir);

    virtual bool calcProduct(const double* pos, const double* u, const double* v, double& prod, bool preCalcMetric = true);
    virtual bool calcProduct(const vec4 pos, const vec4 u, const vec4 v, double& prod, bool preCalcMetric = true);

    virtual double coordDiff(const unsigned int coordNum, const double p1, const double p2);
    virtual bool calcSepDist(const vec4 p1, const vec4 p2, double& spaceDist, double& timeDist);

    bool calcTidalMatrix(const double* pos, const double* e0, const double* e1, const double* e2, const double* e3);
    bool calcTidalMatrix(const vec4 pos, const vec4 e0, const vec4 e1, const vec4 e2, const vec4 e3);
    double getTidalMatrixCoeff(const int i, const int j);
    void getTidalMatrix(mat4& m);

    bool calcCoordTidalMatrix(const double* pos, const double* u);
    bool calcCoordTidalMatrix(const vec4 pos, const vec4 u);
    double getCoordTidalMatrixCoeff(const int i, const int j);
    void getCoordTidalMatrix(mat4& m);

    virtual void calcFmu_nu(const double* pos);
    double getFmu_nu(const int mu, const int nu);

    virtual enum_physical_constants physicalUnits();
    virtual void usePhysicalUnits(const enum_physical_constants units);
    virtual void setUnits(const double speed_of_light, const double grav_const, const double diel_perm);

    double speed_of_light();
    double grav_constant();
    double dielectric_perm();

    virtual bool addParam(const char* pName, double val = 0.0);
    virtual bool setParam(const char* pName, double val);
    virtual bool getParam(const char* pName, double& val);
    //virtual bool    getParam(int pNr, std::string& pName, double& val);
    virtual bool getParam(int pNr, char*& pName, double& val);
    virtual bool setParam(int pNr, double val);

    int getNumParams();
    void getParamNames(std::vector<std::string>& names);
    int getParamNum(const char* name);
    const char* getParamName(int pNr);

    int getLocTedTypes(std::vector<enum_nat_tetrad_type>& locted);
    enum_nat_tetrad_type getCurrLTtype(int num);

    int getDrawTypes(std::vector<enum_draw_type>& drawTypes);
    enum_draw_type getCurrDrawType(int num);
    enum_draw_type getCurrDrawType(std::string name);

    virtual int transToPseudoCart(vec4 p, vec4& cp);
    virtual bool transToEmbedding(vec4 p, vec4& ep);
    virtual bool transToTwoPlusOne(vec4 p, vec4& cp);
    virtual bool transToCustom(vec4 p, vec4& cp);

    virtual void transFromPseudoCart(vec4 cp, vec4& p);

    virtual bool addEmbeddingParam(const char* name, double val = 0.0);
    virtual bool setEmbeddingParam(const char* name, double val);
    virtual bool getEmbeddingParam(const char* name, double& val);

    virtual void getEmbeddingNames(std::vector<std::string>& names);
    virtual bool getAllEmbeddingParams(std::vector<std::string>& names, std::vector<double>& params);
    virtual bool getEmbeddingMap(std::map<std::string, double>& params);
    virtual int getEmbeddingVertices(std::vector<vec3>& verts,
        std::vector<int>& indices, unsigned int& numElems, unsigned int& counter);

    virtual bool haveEmbedding();

    virtual bool effPotentialValue(const vec4 pos, const vec4 cdir, enum_geodesic_type type, const double x, double& val);
    virtual bool totEnergy(const vec4 pos, const vec4 cdir, const double x, double& val);
    virtual bool haveEffPotential();

    virtual double getCircularVelocity(const double r, const enum_nat_tetrad_type tedType = enum_nat_tetrad_default);
    virtual vec4 getCircularFourVel(const vec4 pos, const enum_nat_tetrad_type tedType = enum_nat_tetrad_default);

    virtual vec4 getStandardInitialPosition();
    virtual vec3 getStandardInitialDirection();

    bool isResizeEnabled();

    bool isChrisDAvailable();

    // virtual bool    report(const vec4 pos, const vec4 cdir, std::string &text);
    virtual void printF(FILE* fptr = stderr);

protected:
    void lowCase(std::string& s);
    void setCoordType(enum_coordinate_type coord_type);

    /**
     * @brief Set coordinate names, initial coordinate position and initial local direction.
     */
    virtual void setStandardValues();
    virtual void contrChrisVecVec(const double* y, const double* v, const double* w, double* z, bool calc = true);
    virtual void contrChrDVecVecVec(const double* y, const double* u, const double* v, const double* w, double* z, bool calc = true);

    // -------- protected attribute ---------
protected:
    std::string mMetricName; //!< Calling name of the metric.
    std::string mMetricCPPfilename; //!< File name of the metric.
    enum_coordinate_type mCoordType; //!< Type of coordinates.
    std::string mCoordNames[4]; //!< Name of coordinates.
    struct_scoord_type mScoordType[4]; //!< Single coordinate type.
    double mSign; //!< Signature of the metric:

    double g_compts[4][4]; //!< Metric coefficients.
    double christoffel[4][4][4]; //!< Christoffel symbols of the second kind.
    double chrisD[4][4][4][4]; //!< First derivative of the Christoffel symbols.

    double riem[4][4][4][4]; //!< Riemann tensor R^a_bcd
    double weyl[4][4][4][4]; //!< Weyl tensor.
    double ric[4][4]; //!< Ricci tensor.
    double rrc[4][4][4]; //!< Ricci rotation coefficients.
    double crrc[4]; //!< Contractions of the Ricci rotation coefficients.

    mat4 tidalMatrix;
    mat4 coordTidalMatrix;

    double fmu_nu[4][4]; //!< electromagnetic field tensor F^mu_nu

    bool inPhysicalUnits;
    enum_physical_constants mPhysicalUnits;

    std::map<std::string, double> mParam; //!< Map of parameters.
    std::map<std::string, double>::iterator mParamItr;
    int mNumParam;

    std::vector<enum_nat_tetrad_type> mLocTeds; //!< Local tetrads defined for the metric.
    std::vector<enum_draw_type> mDrawTypes; //!< Drawing types defined for the metric.

    double mSpeedOfLight; //!< Speed of light (SI: m/s).
    double mGravConstant; //!< Newton's gravitational constant (SI: m^3/(kg*s^2)).
    double mDielectricPerm; //!< Dielectric vacuum permittivity (SI: As/(Vm)).

    double mInitPos[4]; //!< Initial position in coordinates.
    double mInitDir[3]; //!< Initial direction of light ray with respect to local tetrad.

    std::map<std::string, double> mEmbParam; //!< Map of embedding parameters.
    std::map<std::string, double>::iterator mEmbParamItr;
    bool mHaveEmbedding;
    bool mHaveEffPotential;

    //! Enable resize capability. Set to yes if the tangent of the geodesic shall be resized below a certain threshold.
    bool mEnableResize;

    //! Factor for prolate spheroidal Coordinates
    double mprolSphfac;

    bool mHaveChrisD;
};

} // end namespace m4d

#endif
