/**
 * @file    m4dJacobi.cpp
 * @author  Thomas Mueller
 *
 * @brief  Jacobi Elliptic functions.
 *
 * This file is part of the m4d-library.
 */
#include "m4dJacobi.h"
#include <cmath>

namespace m4d {

void sncndn(const double u, const double m, const double eps, double& sn, double& cn, double& dn)
{
    double a[15];
    double b[15];
    double c[15];
    double phi[15];

    double mu, v, hn;
    if (m > 1) {
        mu = 1.0 / m;
        v = u * sqrt(m);
    }
    else {
        mu = m;
        v = u;
    }
    a[0] = 1.0;
    b[0] = sqrt(1.0 - mu);
    c[0] = sqrt(mu);
    int N = 0;

    do {
        N++;
        a[N] = 0.5 * (a[N - 1] + b[N - 1]);
        b[N] = sqrt(a[N - 1] * b[N - 1]);
        c[N] = 0.5 * (a[N - 1] - b[N - 1]);
    } while (fabs(c[N]) > eps);

    phi[N] = pow(2.0, N) * a[N] * v;
    for (int i = N; i > 0; i--) {
        phi[i - 1] = 0.5 * (asin(c[i] / a[i] * sin(phi[i])) + phi[i]);
    }

    sn = sin(phi[0]);
    cn = cos(phi[0]);
    dn = cn / cos(phi[1] - phi[0]);
    if (m > 1) {
        sn *= sqrt(mu);
        hn = dn;
        dn = cn;
        cn = hn;
    }
}

} // end namespace m4d
