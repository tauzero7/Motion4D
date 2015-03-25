// --------------------------------------------------------------------------------
/*
    m4dJacobi.h

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

/*!  \file   m4dJacobi.h
     \brief  Jacobi Elliptic functions.

*/
// --------------------------------------------------------------------------------

#ifndef M4D_JACOBI_H
#define M4D_JACOBI_H


namespace m4d {

void  sncndn(const double u, const double m, const double eps,
             double &sn, double &cn, double &dn);

} // end namespace m4d

#endif


























