// --------------------------------------------------------------------------------
/*
    TransCoordinates.h

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

/*!  \class  m4d::TransCoordinates
     \brief  Coordinate transformations cartesian <-> spherical <-> ...

             each possible tansformation has been implemented in 3 ways:
               - transform point
               - transform point and single direction
               - transform point and four directions (e.g. tetrad)
             it is possible to transform single points/dirs or whole
             arrays of points/diections:
               - directly (transCoordAB)
               - to a specific coordiante system (coordTransf)
               - to cartesian coordiantes (toCartesianCoord)

             Code is duplicated often to minimize function calls
*/
// --------------------------------------------------------------------------------

#ifndef M4D_TRANS_COORDINATES_H
#define M4D_TRANS_COORDINATES_H

#include <m4dGlobalDefs.h>

namespace m4d {

class API_EXPORT TransCoordinates {
public:
    TransCoordinates();
    ~TransCoordinates();

    //! transform position to cartesian coordinates
    static void toCartesianCoord(enum_coordinate_type fromCoord,
                                 const vec4& oldPos, vec4& newPos);

    static void toCartesianCoord(enum_coordinate_type fromCoord,
                                 const vec4& oldPos, const vec4& oldDir,
                                 vec4& newPos, vec4& newDir);

    static void toCartesianCoord(enum_coordinate_type fromCoord,
                                 const vec4& oldPos,
                                 const vec4& oldDir0, const vec4& oldDir1, const vec4& oldDir2, const vec4& oldDir3,
                                 vec4& newPos,
                                 vec4& newDir0, vec4& newDir1, vec4& newDir2, vec4& newDir3);

    //! transform position to cartesian coordinates (ray)
    static void toCartesianCoord(enum_coordinate_type fromCoord, const std::vector<vec4>& oldPos,
                                 std::vector<vec4>& newPos);
    static void toCartesianCoord(enum_coordinate_type fromCoord, const std::vector<vec4>& oldPos, const std::vector<vec4>& oldDir,
                                 std::vector<vec4>& newPos, std::vector<vec4>& newDir);
    static void toCartesianCoord(enum_coordinate_type fromCoord, const std::vector<vec4>& oldPos,
                                 const std::vector<vec4>& oldDir0, const std::vector<vec4>& oldDir1,
                                 const std::vector<vec4>& oldDir2, const std::vector<vec4>& oldDir3,
                                 std::vector<vec4>& newPos,
                                 std::vector<vec4>& newDir0, std::vector<vec4>& newDir1,
                                 std::vector<vec4>& newDir2, std::vector<vec4>& newDir3);

    //! coordinate transformation from 'fromCoord' to 'toCoord' coordinates
    static void coordTransf(enum_coordinate_type fromCoord, const vec4& oldPos,
                            enum_coordinate_type toCoord,   vec4& newPos);
    static void coordTransf(enum_coordinate_type fromCoord, const vec4& oldPos, const vec4& oldDir,
                            enum_coordinate_type toCoord,   vec4& newPos, vec4& newDir);
    static void coordTransf(enum_coordinate_type fromCoord, const vec4& oldPos,
                            const vec4& oldDir0, const vec4& oldDir1, const vec4& oldDir2, const vec4& oldDir3,
                            enum_coordinate_type toCoord, vec4& newPos,
                            vec4& newDir0, vec4& newDir1, vec4& newDir2, vec4& newDir3);
    //! coordinate transformation from 'fromCoord' to 'toCoord' coordinates (ray)
    static void coordTransf(enum_coordinate_type fromCoord, const std::vector<vec4>& oldPos,
                            enum_coordinate_type toCoord,   std::vector<vec4>& newPos);
    static void coordTransf(enum_coordinate_type fromCoord, const std::vector<vec4>& oldPos, const std::vector<vec4>& oldDir,
                            enum_coordinate_type toCoord, std::vector<vec4>& newPos, std::vector<vec4>& newDir);
    static void coordTransf(enum_coordinate_type fromCoord, const std::vector<vec4>& oldPos,
                            const std::vector<vec4>& oldDir0, const std::vector<vec4>& oldDir1,
                            const std::vector<vec4>& oldDir2, const std::vector<vec4>& oldDir3,
                            enum_coordinate_type toCoord, std::vector<vec4>& newPos,
                            std::vector<vec4>& newDir0, std::vector<vec4>& newDir1,
                            std::vector<vec4>& newDir2, std::vector<vec4>& newDir3);


    //! coordinate transformation cart->sph
    static void transCoordCartSph(const vec4& oldPos, vec4& newPos);
    static void transCoordCartSph(const vec4& oldPos, const vec4& oldDir,
                                  vec4& newPos, vec4& newDir);
    static void transCoordCartSph(const vec4& oldPos, const vec4& oldDir0, const vec4& oldDir1, const vec4& oldDir2, const vec4& oldDir3,
                                  vec4& newPos, vec4& newDir0, vec4& newDir1, vec4& newDir2, vec4& newDir3);
    //! coordinate transformation sph->cart
    static void transCoordSphCart(const vec4& oldPos, vec4& newPos);
    static void transCoordSphCart(const vec4& oldPos, const vec4& oldDir,
                                  vec4& newPos, vec4& newDir);
    static void transCoordSphCart(const vec4& oldPos, const vec4& oldDir0, const vec4& oldDir1, const vec4& oldDir2, const vec4& oldDir3,
                                  vec4& newPos, vec4& newDir0, vec4& newDir1, vec4& newDir2, vec4& newDir3);
    //! coordinate transformation cart->cyl
    static void transCoordCartCyl(const vec4& oldPos, vec4& newPos);
    static void transCoordCartCyl(const vec4& oldPos, const vec4& oldDir,
                                  vec4& newPos, vec4& newDir);
    static void transCoordCartCyl(const vec4& oldPos, const vec4& oldDir0, const vec4& oldDir1, const vec4& oldDir2, const vec4& oldDir3,
                                  vec4& newPos, vec4& newDir0, vec4& newDir1, vec4& newDir2, vec4& newDir3);
    //! coordinate transformation cyl->cart
    static void transCoordCylCart(const vec4& oldPos, vec4& newPos);
    static void transCoordCylCart(const vec4& oldPos, const vec4& oldDir,
                                  vec4& newPos, vec4& newDir);
    static void transCoordCylCart(const vec4& oldPos, const vec4& oldDir0, const vec4& oldDir1, const vec4& oldDir2, const vec4& oldDir3,
                                  vec4& newPos, vec4& newDir0, vec4& newDir1, vec4& newDir2, vec4& newDir3);


    //! coordinate transformation prolate spheroidal->cart
    static void transCoordProlSphCart(const vec4& oldPos, vec4& newPos, double a);
    static void transCoordProlSphCart(const vec4& oldPos, const vec4& oldDir,
                                      vec4& newPos, vec4& newDir);
    static void transCoordProlSphCart(const vec4& oldPos, const vec4& oldDir0, const vec4& oldDir1, const vec4& oldDir2, const vec4& oldDir3,
                                      vec4& newPos, vec4& newDir0, vec4& newDir1, vec4& newDir2, vec4& newDir3);


};

} // end namespace m4d

#endif
