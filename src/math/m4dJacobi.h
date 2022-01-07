/**
 * @file    m4dJacobi.h
 * @author  Thomas Mueller
 *
 * @brief  Jacobi Elliptic functions.
 *
 * This file is part of the m4d-library.
 */
#ifndef M4D_JACOBI_H
#define M4D_JACOBI_H

namespace m4d {

void sncndn(const double u, const double m, const double eps, double& sn, double& cn, double& dn);

} // end namespace m4d

#endif
