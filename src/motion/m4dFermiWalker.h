/**
 * @file    m4dFermiWalker.h
 * @author  Thomas Mueller

     \brief  The local tetrad of an arbitrarily moving object has to be Fermi-Walker transported.

        The Fermi-Walker transport equation for the local tetrad \f$\mathbf{e}_{(i)}=e_{(i)}^{\mu}\partial_{\mu}\f$ is
 given by: \f[ \frac{de_{(i)}^{\mu}}{d\tau} = -\Gamma_{\alpha\beta}^{\mu}e_{(0)}^{\alpha}e_{(i)}^{\beta} +
 \frac{1}{c}\left(\eta_{(i)(k)}a^{(k)}e_{(0)}^{\mu}-\eta_{(0)(i)}a^{(k)}e_{(k)}^{\mu}\right) \f]


        There are two possibilities to calculate the Fermi-Walker transport of
        a local tetrad. On the one hand, the proper 3-acceleration with respect
        to the local tetrad can be given for each time step. On the other hand,
        the worldline can be defined as a function \f$x=x(\tau)\f$.

        Because the four-acceleration is always orthogonal to the four-velocity,
        the zero- or time- component of the proper acceleration with respect to
        the local tetrad where the e0- base vector is tangential to the four-velocity
        vanishes.

        The worldline function \f$x=x(\tau)\f$ has to take care of the
        parameter pointer!


        Note that e0 is tangential to four-velocity : y[8]..y[11] = y[4]..y[7].
 *
 *  This file is part of libMotion4D.
 */
#ifndef M4D_MOTION_FERMIWALKER_H
#define M4D_MOTION_FERMIWALKER_H

#include <iostream>

#include "m4dMotion.h"

namespace m4d {

/**
 * @brief The FermiWalker class
 */
class API_M4D_EXPORT FermiWalker : public Motion
{
public:
    explicit FermiWalker(Metric* metric);
    virtual ~FermiWalker();

    /*! Set acceleration with respect to local tetrad.
     *
     *  \param a1 : proper acceleration in e1-direction.
     *  \param a2 : proper acceleration in e2-direction.
     *  \param a3 : proper acceleration in e3-direction.
     */
    void setCurrPropAccel(double a1, double a2, double a3);

    /*! Get acceleration with respect to local tetrad.
     *
     *  \param a1 : reference to proper acceleration in e1-direction.
     *  \param a2 : reference to proper acceleration in e2-direction.
     *  \param a3 : reference to proper acceleration in e3-direction.
     */
    void getCurrPropAccel(double& a1, double& a2, double& a3);

    /*! Set initial velocity with respect to natural local tetrad.
     *
     *  \param fm : time direction (>0 : future, <0 : past).
     *  \param v  : absolute value of velocity with respect to natural local tetrad.
     *  \param theta :  angle
     *  \param phi   :  angle
     *  \param type  : type of natural local tetrad.
     *  \sa enum_nat_tetrad_type.
     */
    bool setInitialVelocity(
        double fm, double v, double theta, double phi, enum_nat_tetrad_type type = enum_nat_tetrad_default);

    /// Get initial velocity with respect to natural local tetrad.
    void getInitialVelocity(double& v, double& theta, double& phi);

    //! Reset proper time.
    void resetProperTime() { mLambda = 0.0; }

    /*! Calculate motion with at most 'maxNumPoints' points.
     *
     *  \param initPos : initial position.
     *  \param fm : time direction (fm>0: future, fm<0: past).
     *  \param v :  velocity.
     *  \param theta_v : theta-direction of velocity-
     *  \param phi_v : phi-direction of velocity.
     *  \param a :  acceleration.
     *  \param theta_a :  theta-direction of acceleration.
     *  \param phi_a :  phi-direction of acceleration.
     *  \param e0 : local tetrad vector.
     *  \param e1 : local tetrad vector.
     *  \param e2 : local tetrad vector.
     *  \param e3 : local tetrad vector.
     *  \param maxNumPoints : maximum number of points.
     *  \param points :  reference to points.
     *  \param base0 : reference to base vector.
     *  \param base1 : reference to base vector.
     *  \param base2 : reference to base vector.
     *  \param base3 : reference to base vector.
     */
    virtual enum_break_condition calculateMotion(const vec4 initPos, double fm, double v, double theta_v, double phi_v,
        double a, double theta_a, double phi_a, const vec4 e0, const vec4 e1, const vec4 e2, const vec4 e3,
        const int maxNumPoints, std::vector<vec4>& points, std::vector<vec4>& base0, std::vector<vec4>& base1,
        std::vector<vec4>& base2, std::vector<vec4>& base3);

    /*! Set 'true' if the Fermi-Walker transport is along a given worldline.
     *
     *  \param withworldline : (true:) calculate Fermi-Walker transport with a predefined wordline, (false:)
     * acceleration with respect to the local tetrad.
     */
    void setCalcWithWorldline(bool withworldline);

    /*! Set the worldline function.
     *  \param func : function pointer to the worldline x=x(tau).
     */
    void set_x_tau(vec4 (*func)(double, void*));

    /*! Set the tangent to the worldline function.
     *  \param func : function pointer to the four-velocity of the worldline x=x(tau).
     */
    void set_u_tau(vec4 (*func)(double, void*));

    /*!  Set the acceleration to the worldline function.
     *  \param func : function pointer to the four-acceleration of the worldline x=x(tau).
     */
    void set_a_tau(vec4 (*func)(double, void*));

    /*!  Set the parameters to the worldline function.
     *  \param params : pointer to parameters.
     */
    void set_params(void* params);

    /*!  Evaluate the worldline for proper time tau.
     *  \param tau : proper time.
     *  \param x   : reference to current position.
     *  \return true : worldline x=x(tau) is given.
     *  \return false : worlined is not defined.
     */
    bool get_x_tau(double tau, vec4& x);

    /*!   Evaluate tangent to the worldline for proper time tau.
     *  \param tau : proper time.
     *  \param u   : reference to current four-velocity.
     *  \return true : worldline x=x(tau) is given.
     *  \return false : worlined is not defined.
     */
    bool get_u_tau(double tau, vec4& u);

    /*!  Evaluate acceleration to the worldline for proper time tau.
     *  \param tau : proper time.
     *  \param a   : reference to current four-acceleration.
     *  \return true : worldline x=x(tau) is given.
     *  \return false : worlined is not defined.
     */
    bool get_a_tau(double tau, vec4& a);

    /*! Initialize worldline.
     *  \param tauStart : initial proper time.
     */
    bool initWorldline(double tauStart);

    /*!  Update worldline to current time tau.
     * \param tau : current proper time.
     */
    bool updateWorldline(double tau);

    /*! Calculate the next step of the Fermi-Walker transport.
     *  \return enum_break_condition.
     */
    enum_break_condition nextStep();

    /*! Calculate the next step of the Fermi-Walker transport along the worldline.
     *  \return enum_break_condition.
     */
    enum_break_condition nextStepWL();

protected:
    /*! Calculate right hand side of Fermi-Walker transport.
     * \param y[] : pointer to y.
     * \param dydx[] : pointer to right hand side of geodesic equation.
     */
    bool calcDerivs(const double y[], double dydx[]);

    /*! Calculate right hand side of Fermi-Walker transport for a given worldline.
     * \param y[] : pointer to y.
     * \param dydx[] : pointer to right hand side of geodesic equation.
     */
    bool calcDerivsWL(const double y[], double dydx[]);

protected:
    //! Acceleration with respect to local tetrad.
    double mPropAcc[4];
    //! Current four-velocity in coordinates.
    double mInitVel[4];
    //! Current four-acceleration.
    vec4 mCurrAcc;

    double mVel;
    double mTheta;
    double mPhi;

    bool mCalcWithWorldline;
    vec4 (*x_tau)(double, void*);
    vec4 (*u_tau)(double, void*);
    vec4 (*a_tau)(double, void*);
    void* mParams;
};

} // end namespace m4d

#endif
