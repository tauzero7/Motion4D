/**
 * @file    m4dGeodesic.h
 * @author  Thomas Mueller
 *
     \brief  Base class for all geodesics.

       The Geodesic class cannot be used itself due to the pure virtual
       methods calculateGeodesic() and calculateGeodesicData().

       The geodesic equation reads
       \f[ \frac{d^2x^{\mu}}{d\lambda^2} +
 \Gamma_{\alpha\beta}^{\mu}\frac{dx^{\alpha}}{d\lambda}\frac{dx^{\beta}}{d\lambda} = 0, \f] where \f$\lambda\f$ is an
 affine parameter and \f$\Gamma_{\alpha\beta}^{\mu}\f$ are the Christoffell symbols of the second kind defined by \f[
 \Gamma_{\alpha\beta}^{\mu} =
 \frac{1}{2}g^{\mu\rho}\left(g_{\rho\alpha,\beta}+g_{\rho\beta,\alpha}-g_{\alpha\beta,\rho}\right). \f] In the case of a
 timelike geodesic, this affine parameter equals the proper time. A timelike \f$(\kappa=-1)\f$ or a lightlike
 \f$(\kappa=0)\f$ geodesic has to fulfill the constraint equation \f[
 g_{\mu\nu}\frac{dx^{\mu}}{d\lambda}\frac{dx^{\nu}}{d\lambda}=\kappa c^2.\f]

       The parallel transport equation for a vector X along the geodesic with tangent u reads
       \f[ \frac{dX^{\mu}}{d\lambda} + \Gamma_{\alpha\beta}^{\mu}u^{\alpha}X^{\beta} = 0. \f]

 *
 *  This file is part of libMotion4D.
 */
#ifndef M4D_GEODESIC_H
#define M4D_GEODESIC_H

#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>

#include <m4dGlobalDefs.h>
#include <metric/m4dMetric.h>
#include <motion/m4dMotion.h>

namespace m4d {

/**
 * @brief The Geodesic class
 */
class API_M4D_EXPORT Geodesic : public Motion
{
public:
    Geodesic(Metric* metric, enum_geodesic_type type = enum_geodesic_lightlike);
    virtual ~Geodesic();

    /*!
     * @brief  Set type of geodesic.
     * @param type  Type of geodesic: lightlike, timelike,
     */
    void setGeodesicType(enum_geodesic_type type);

    /// Get type of geodesic.
    enum_geodesic_type type();

    /*!
     * @brief Set parameter via key-value.
     * @param paramName   Name of parameter.
     * @param val         Parameter value to be set.
     */
    virtual bool setParam(std::string paramName, bool val);

    /*!
     * @brief Set parameter via key-value.
     * @param paramName   Name of parameter.
     * @param value       Parameter value to be set.
     */
    virtual bool setParam(std::string paramName, double value);

    /*!
     * @brief Set four parameters via key-values.
     * @param paramName   Name of parameter.
     * @param v0          Parameter v0.
     * @param v1          Parameter v1.
     * @param v2          Parameter v2.
     * @param v3          Parameter v3.
     */
    virtual bool setParam(std::string paramName, double v0, double v1, double v2, double v3);

    /*!
     * @brief Set absolute and relative epsilon.
     * @param  eps_a   epsilon absolute.
     * @param  eps_r   epsilon relative.
     */
    void setEpsilons(double eps_a, double eps_r);

    /*!
     * @brief Get absolute and relative epsilons.
     * @param eps_a : reference to epsilon absolute.
     * @param eps_r : reference to epsilon relative.
     */
    void getEpsilons(double& eps_a, double& eps_r);

    /// Set/Unset stepsize control.
    void setStepSizeControlled(bool control = true);

    /// Set if the local tetrad of the observer should be parallely transported along the null geodesic.
    void setCalcWithParTransport(bool calcwith = false);
    bool calcWithParTransport();

    void setResize(double eps, double factor);
    void getResize(double& eps, double& factor);

    /**
     * @brief Initialize geodesic.
     * @param initPos : initial position.
     * @param initDir : initial coordinate direction.
     * @param cstr
     */
    virtual enum_break_condition initializeGeodesic(const vec4 initPos, const vec4 initDir, double& cstr);

    /**
     * @brief Calculate a geodesic.
     *  @param  initPos      :  initial position of the geodesic in coordinates.
     *  @param  initDir      :  initial direction of the geodesic in coordinates.
     *  @param  maxNumPoints :  maximum number of points.
     *  @param  points       :  reference to the calculated points.
     *  @param  dirs         :  reference to the calculated tangents.
     *  @param  lambda       :  reference to the affine parameters.
     *
     *  \return enum_break_condition : break condition.
     *  \sa enum_break_condition
     */
    virtual enum_break_condition calculateGeodesic(const vec4 initPos, const vec4 initDir, const int maxNumPoints,
        std::vector<vec4>& points, std::vector<vec4>& dirs, std::vector<double>& lambda)
        = 0;

    virtual enum_break_condition calculateGeodesic(
        const vec4 initPos, const vec4 initDir, const int maxNumPoints, vec4*& points, vec4*& dirs, int& numPoints)
        = 0;

    /**
     * @brief Calculate a geodesic and store points, directions, and constraint.
     *  @param  initPos      :  initial position of the geodesic in coordinates.
     *  @param  initDir      :  initial direction of the geodesic in coordinates.
     *  @param  maxNumPoints :  maximum number of points.
     *  @param  points       :  reference to the calculated points.
     *  @param  dirs         :  reference to the calculated tangents.
     *  @param  epsilons     :  reference to the epsilons.
     *  @param  lambda       :  reference to the affine parameters.
     *
     *  \return enum_break_condition : break condition.
     *  \sa enum_break_condition
     */
    virtual enum_break_condition calculateGeodesicData(const vec4 initPos, const vec4 initDir, const int maxNumPoints,
        std::vector<vec4>& points, std::vector<vec4>& dirs, std::vector<double>& epsilons, std::vector<double>& lambda)
        = 0;

    virtual enum_break_condition calcParTransport(const vec4 initPos, const vec4 initDir, const vec4 e0, const vec4 e1,
        const vec4 e2, const vec4 e3, const int maxNumPoints, std::vector<vec4>& points, std::vector<vec4>& dirs,
        std::vector<double>& lambda, std::vector<vec4>& base0, std::vector<vec4>& base1, std::vector<vec4>& base2,
        std::vector<vec4>& base3)
        = 0;

    /**
     * @brief Calculate a geodesic and the parallel transported local tetrad of the observer.
     *  @param  initPos      :  initial position of the geodesic in coordinates.
     *  @param  initCoordDir :  initial coordinate direction of the geodesic in coordinates.
     *  @param  localNullDir :  initial local direction of the geodesic in coordinates.
     *  @param  locX
     *  @param  locY
     *  @param  locZ
     *  @param  b0           :  b0 base vectors of local tetrad.
     *  @param  b1           :  b1 base vectors of local tetrad.
     *  @param  b2           :  b2 base vectors of local tetrad.
     *  @param  b3           :  b3 base vectors of local tetrad.
     *  @param  tetrad_type  :  type of local tetrad.
     *  @param  maxNumPoints :  maximum number of points.
     *  @param  points       :  reference to the calculated points.
     *  @param  dirs         :  reference to the calculated directions.
     *  @param  lambda       :  reference to the affine parameters.
     *  @param  sachs0       :  references to Sachs basis 0.
     *  @param  sachs1       :  references to Sachs basis 1.
     *  @param  jacobi       :  references to Jacobi parameters.
     *  @param  maxJacobi    :  references to maximum Jacobi parameter.
     *
     *  \return enum_break_condition : break condition.
     *  \sa enum_break_condition
     */
    virtual enum_break_condition calcSachsJacobi(const vec4 initPos, const vec4 initCoordDir, const vec3 localNullDir,
        const vec3 locX, const vec3 locY, const vec3 locZ, const vec4 e0, const vec4 e1, const vec4 e2, const vec4 e3,
        const enum_nat_tetrad_type tetrad_type, const int maxNumPoints, std::vector<vec4>& points,
        std::vector<vec4>& dirs, std::vector<double>& lambda, std::vector<vec4>& sachs0, std::vector<vec4>& sachs1,
        std::vector<vec5>& jacobi, vec5& maxJacobi)
        = 0;

    virtual enum_break_condition calcSachsJacobi(const vec4 initPos, const vec4 initCoordDir, const vec3 localNullDir,
        const vec3 locX, const vec3 locY, const vec3 locZ, const vec4 e0, const vec4 e1, const vec4 e2, const vec4 e3,
        const enum_nat_tetrad_type tetrad_type, const int maxNumPoints, vec4*& points, vec4*& dirs, double*& lambda,
        vec4*& sachs0, vec4*& sachs1, vec5*& jacobi, vec5& maxJacobi, int& numPoints)
        = 0;

    virtual bool nextStep(int& status) = 0;
    virtual bool nextStepPar(int& status) = 0;
    virtual bool nextStepSachsJacobi(int& status) = 0;

    /**
     * @brief Test, whether the constraint condition is fulfilled.
     * @return constraint value.
     */
    virtual double testConstraint();

    virtual void printF(FILE* fptr = stderr);

protected:
    void setKappa();

    /**
     * @brief Calculate the right side of the geodesic equation.
     * @param y[] : pointer to y.
     * @param dydx[] : pointer to right hand side of geodesic equation.
     */
    virtual bool calcDerivs(const double y[], double dydx[]);

    /**
     * @brief Calculate the right side of the geodesic equation with parallel transport.
     * @param y[] : pointer to y.
     * @param dydx[] : pointer to right hand side of parallel transport equation.
     */
    virtual bool calcDerivsPar(const double y[], double dydx[]);

    /**
     * @brief Calculate the right side of the geodesic equation with parallel transport and Jacobi equation
     * @param y[] : pointer to y.
     * @param dydx[] : pointer to right hand side of parallel transport equation.
     */
    virtual bool calcDerivsSachsJacobi(const double y[], double dydx[]);

    /**
     * @brief Calculate Sachs basis from the initial local direction.
     * @param localNullDir: initial local direction of the light ray.
     * @param locX : local x-direction.
     * @param locY : local y-direction.
     * @param locZ : local z-direction.
     */
    void calcSachsBasis(const vec3 localNullDir, const vec3 locX, const vec3 loxY, const vec3 locZ);

    /**
     * @brief Set Sachs basis.
     * @param s1 : first sachs vector.
     * @param s2 : second sachs vector.
     */
    void setSachsBasis(const vec4 s1, const vec4 s2);

    /**
     * @brief Calculate the Jacobi parameters.
     * @param lambda : affine parameter.
     * @param y[]    : array
     * @param currJacobi : reference to jacobi values.
     */
    void calcJacobiParams(const double lambda, const double y[], vec5& currJacobi);

    void findMaxJacobi(vec5& currJacobi, vec5& maxJacobi);

protected:
    //! Type of geodesic.
    enum_geodesic_type mType;
    double mKappa;

    bool mCalcWithParTransport;
    int mNumCoords;

    bool mStepsizeControlled;
    //! Absolute epsilon.
    double epsilon_abs;
    //! Relative epsilon.
    double epsilon_rel;
    //! Epsilon for resize of y[4].
    double resizeEps;
    double resizeFac;

    vec3 mSachsBasisB1;
    vec3 mSachsBasisB2;
};

} // end namespace m4d

#endif
