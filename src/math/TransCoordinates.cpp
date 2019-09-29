// -------------------------------------------------------------------------------
/*
    TransCoordinates.cpp

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

#include "TransCoordinates.h"
#include <cassert>
#include <cmath>
#include <iostream>

namespace m4d {

// ---------------------------------------------------
//    constructur/destructor
// ---------------------------------------------------
TransCoordinates::TransCoordinates() {}

TransCoordinates::~TransCoordinates() {}

/*!
 * transform POINT to cartesian coordinates
 * \param fromCoord : old coordinate system.
 * \param oldPos    : old position in old coordinate system.
 * \param newPos    : new position in cartesian coordinates.
 */
void TransCoordinates::toCartesianCoord(enum_coordinate_type fromCoord,
    const vec4& oldPos, vec4& newPos)
{
    switch (fromCoord) {
    case (enum_coordinate_spherical):
        TransCoordinates::transCoordSphCart(oldPos, newPos);
        break;
    case (enum_coordinate_cylinder):
        TransCoordinates::transCoordCylCart(oldPos, newPos);
        break;
    case (enum_coordinate_cartesian):
        newPos = oldPos;
        break;
    case (enum_coordinate_custom):
        newPos = oldPos;
        break;
    default:
        break;
    }
}

/*!
 * transform POINT to cartesian coordinates
 * \param fromCoord : old coordinate system
 * \param oldPos : old position in old coordinate system
 * \param oldDir : old direction in old coordinates system
 * \param newPos : new position in cartesian coordinates
 * \param newDir : new direction in cartesian coordinates
 */
void TransCoordinates::toCartesianCoord(enum_coordinate_type fromCoord,
    const vec4& oldPos, const vec4& oldDir,
    vec4& newPos, vec4& newDir)
{
    switch (fromCoord) {
    case (enum_coordinate_spherical):
        TransCoordinates::transCoordSphCart(oldPos, oldDir, newPos, newDir);
        break;
    case (enum_coordinate_cylinder):
        TransCoordinates::transCoordCylCart(oldPos, oldDir, newPos, newDir);
        break;
    case (enum_coordinate_cartesian):
        newPos = oldPos;
        newDir = oldDir;
        break;
    case (enum_coordinate_custom):
        newPos = oldPos;
        newDir = oldDir;
        break;
    default:
        break;
    }
}

/*!
 * transform POINT to cartesian coordinates
 * \param fromCoord : old coordinate system.
 * \param oldPos    : old position in old coordinate system.
 * \param oldDir0   : old direction in old coordinate system.
 * \param oldDir1   : old direction in old coordinate system.
 * \param oldDir2   : old direction in old coordinate system.
 * \param oldDir3   : old direction in old coordinate system.
 * \param newPos    : new position in cartesian coordinates.
 * \param newDir0   : new direction in cartesian coordinates.
 * \param newDir1   : new direction in cartesian coordinates.
 * \param newDir2   : new direction in cartesian coordinates.
 * \param newDir3   : new direction in cartesian coordinates.
 */
void TransCoordinates::toCartesianCoord(
    enum_coordinate_type fromCoord, const vec4& oldPos, const vec4& oldDir0,
    const vec4& oldDir1, const vec4& oldDir2, const vec4& oldDir3, vec4& newPos,
    vec4& newDir0, vec4& newDir1, vec4& newDir2, vec4& newDir3)
{
    switch (fromCoord) {
    case (enum_coordinate_spherical):
        TransCoordinates::transCoordSphCart(oldPos, oldDir0, oldDir1, oldDir2,
            oldDir3, newPos, newDir0, newDir1,
            newDir2, newDir3);
        break;
    case (enum_coordinate_cylinder):
        TransCoordinates::transCoordCylCart(oldPos, oldDir0, oldDir1, oldDir2,
            oldDir3, newPos, newDir0, newDir1,
            newDir2, newDir3);
        break;
    case (enum_coordinate_cartesian):
        newPos = oldPos;
        newDir0 = oldDir0;
        newDir1 = oldDir1;
        newDir2 = oldDir2;
        newDir3 = oldDir3;
        break;
    case (enum_coordinate_custom):
        newPos = oldPos;
        newDir0 = oldDir0;
        newDir1 = oldDir1;
        newDir2 = oldDir2;
        newDir3 = oldDir3;
        break;
    default:
        break;
    }
}

/*!
 * transform RAY to cartesian coordinates
 * \param fromCoord : old coordinate system
 * \param oldPos : old positions in old coordinate system
 * \param newPos : new positions in cartesian coordinates
 */
void TransCoordinates::toCartesianCoord(enum_coordinate_type fromCoord,
    const std::vector<vec4>& oldPos,
    std::vector<vec4>& newPos)
{
    unsigned int i = 0;

    // resize if necessary (initialization of newPos not required)
    if (newPos.size() != oldPos.size()) {
        newPos.clear();
        newPos.resize(oldPos.size());
    }

    // implementierung ohne aufruf von' toCartesianCoord (fromCoord,vec4& oldPos,
    // ...) damit nicht dauernd der switch abgearbeotet wird. sollte so schneller
    // sein.

    switch (fromCoord) {
    case (enum_coordinate_spherical):
        for (i = 0; i < oldPos.size(); i++) {
            TransCoordinates::transCoordSphCart(oldPos[i], newPos[i]);
        }
        break;
    case (enum_coordinate_cylinder):
        for (i = 0; i < oldPos.size(); i++) {
            TransCoordinates::transCoordCylCart(oldPos[i], newPos[i]);
        }
        break;
    case (enum_coordinate_cartesian):
        for (i = 0; i < oldPos.size(); i++) {
            newPos[i] = oldPos[i];
        }
        break;
    case (enum_coordinate_custom):
        for (i = 0; i < oldPos.size(); i++) {
            newPos[i] = oldPos[i];
        }
        break;
    default:
        break;
    }
}

/*!
 * transform RAY to cartesian coordinates
 * \param fromCoord : old coordinate system.
 * \param oldPos : old positions in old coordinate system.
 * \param oldDir : old directions in old coordinate system.
 * \param newPos : new positions in cartesian coordinates.
 * \param newDir : new directions in cartesian coordiantes.
 */
void TransCoordinates::toCartesianCoord(enum_coordinate_type fromCoord,
    const std::vector<vec4>& oldPos,
    const std::vector<vec4>& oldDir,
    std::vector<vec4>& newPos,
    std::vector<vec4>& newDir)
{
    unsigned int i = 0;

    // resize if necessary (initialization of newPos not required)
    if (newPos.size() != oldPos.size()) {
        newPos.clear();
        newPos.resize(oldPos.size());
    }
    if (newDir.size() != oldDir.size()) {
        newDir.clear();
        newDir.resize(oldDir.size());
    }

    // implementierung ohne aufruf von' toCartesianCoord (fromCoord,vec4& oldPos,
    // ...) damit nicht dauernd der switch abgearbeotet wird. sollte so schneller
    // sein.

    switch (fromCoord) {
    case (enum_coordinate_spherical):
        for (i = 0; i < oldPos.size(); i++) {
            TransCoordinates::transCoordSphCart(oldPos[i], oldDir[i], newPos[i],
                newDir[i]);
        }
        break;
    case (enum_coordinate_cylinder):
        for (i = 0; i < oldPos.size(); i++) {
            TransCoordinates::transCoordCylCart(oldPos[i], oldDir[i], newPos[i],
                newDir[i]);
        }
        break;
    case (enum_coordinate_cartesian):
        for (i = 0; i < oldPos.size(); i++) {
            newPos[i] = oldPos[i];
            newDir[i] = oldDir[i];
        }
        break;
    case (enum_coordinate_custom):
        for (i = 0; i < oldPos.size(); i++) {
            newPos[i] = oldPos[i];
            newDir[i] = oldDir[i];
        }
        break;
    default:
        break;
    }
}

/**
 * transform RAY to cartesian coordinates
 * \param fromCoord : old coordinate system.
 * \param oldPos    : old positions in old coordinate system.
 * \param oldDir0   : old directions in old coordinate system.
 * \param oldDir1   : old directions in old coordinate system.
 * \param oldDir2   : old directions in old coordinate system.
 * \param oldDir3   : old directions in old coordinate system.
 * \param newPos    : new positions in cartesian coordinates.
 * \param newDir0   : new directions in cartesian coordiantes.
 * \param newDir1   : new directions in cartesian coordiantes.
 * \param newDir2   : new directions in cartesian coordiantes.
 * \param newDir3   : new directions in cartesian coordiantes.
 */
void TransCoordinates::toCartesianCoord(
    enum_coordinate_type fromCoord, const std::vector<vec4>& oldPos,
    const std::vector<vec4>& oldDir0, const std::vector<vec4>& oldDir1,
    const std::vector<vec4>& oldDir2, const std::vector<vec4>& oldDir3,
    std::vector<vec4>& newPos, std::vector<vec4>& newDir0,
    std::vector<vec4>& newDir1, std::vector<vec4>& newDir2,
    std::vector<vec4>& newDir3)
{
    unsigned int i = 0;

    // resize if necessary (initialization of newPos not required)
    if (newPos.size() != oldPos.size()) {
        newPos.clear();
        newPos.resize(oldPos.size());
    }
    if (newDir0.size() != oldDir0.size()) {
        newDir0.clear();
        newDir0.resize(oldDir0.size());
    }
    if (newDir1.size() != oldDir1.size()) {
        newDir1.clear();
        newDir1.resize(oldDir1.size());
    }
    if (newDir2.size() != oldDir2.size()) {
        newDir2.clear();
        newDir2.resize(oldDir2.size());
    }
    if (newDir3.size() != oldDir3.size()) {
        newDir3.clear();
        newDir3.resize(oldDir3.size());
    }

    // implementierung ohne aufruf von' toCartesianCoord (fromCoord,vec4& oldPos,
    // ...) damit nicht dauernd der switch abgearbeotet wird. sollte so schneller
    // sein.

    switch (fromCoord) {
    case (enum_coordinate_spherical):
        for (i = 0; i < oldPos.size(); i++)
            TransCoordinates::transCoordSphCart(
                oldPos[i], oldDir0[i], oldDir1[i], oldDir2[i], oldDir3[i], newPos[i],
                newDir0[i], newDir1[i], newDir2[i], newDir3[i]);
        break;
    case (enum_coordinate_cylinder):
        for (i = 0; i < oldPos.size(); i++)
            TransCoordinates::transCoordCylCart(
                oldPos[i], oldDir0[i], oldDir1[i], oldDir2[i], oldDir3[i], newPos[i],
                newDir0[i], newDir1[i], newDir2[i], newDir3[i]);
        break;
    case (enum_coordinate_cartesian):
        for (i = 0; i < oldPos.size(); i++) {
            newPos[i] = oldPos[i];
            newDir0[i] = oldDir0[i];
            newDir1[i] = oldDir1[i];
            newDir2[i] = oldDir2[i];
            newDir3[i] = oldDir3[i];
        }
        break;
    case (enum_coordinate_custom):
        for (i = 0; i < oldPos.size(); i++) {
            newPos[i] = oldPos[i];
            newDir0[i] = oldDir0[i];
            newDir1[i] = oldDir1[i];
            newDir2[i] = oldDir2[i];
            newDir3[i] = oldDir3[i];
        }
        break;
    default:
        break;
    }
}

/*!
 * transform POINT between coordinate systems.
 * \param fromCoord old coordinate system.
 * \param oldPos point in old coordinates.
 * \param toCoord new coordinate system.
 * \param newPos point in new coordinates.
 */
void TransCoordinates::coordTransf(enum_coordinate_type fromCoord,
    const vec4& oldPos,
    enum_coordinate_type toCoord, vec4& newPos)
{
    switch (fromCoord) {
    default:
        break;
    case (enum_coordinate_cartesian): {
        switch (toCoord) {
        case (enum_coordinate_spherical): {
            TransCoordinates::transCoordCartSph(oldPos, newPos);
            break;
        }
        case (enum_coordinate_cylinder): {
            TransCoordinates::transCoordCartCyl(oldPos, newPos);
            break;
        }
        case (enum_coordinate_cartesian): {
            newPos = oldPos;
            break;
        }
        default:
            break;
        }
        break;
    }
    case (enum_coordinate_spherical): {
        switch (toCoord) {
        case (enum_coordinate_cartesian): {
            TransCoordinates::transCoordSphCart(oldPos, newPos);
            break;
        }
        case (enum_coordinate_cylinder): {
            TransCoordinates::transCoordSphCart(oldPos, newPos);
            TransCoordinates::transCoordCartCyl(newPos, newPos);
            break;
        }
        case (enum_coordinate_spherical): {
            newPos = oldPos;
            break;
        }
        default:
            break;
        }
        break;
    }
    case (enum_coordinate_cylinder): {
        switch (toCoord) {
        case (enum_coordinate_cartesian): {
            TransCoordinates::transCoordCylCart(oldPos, newPos);
            break;
        };
        case (enum_coordinate_spherical): {
            TransCoordinates::transCoordCylCart(oldPos, newPos);
            TransCoordinates::transCoordCartSph(newPos, newPos);
            break;
        };
        case (enum_coordinate_cylinder): {
            newPos = oldPos;
            break;
        }
        default:
            break;
        }
        break;
    }
    }
}

/*!
 * transform POINT between coordinate systems.
 * \param fromCoord : old coordinate system.
 * \param oldPos : point in old coordinates.
 * \param oldDir : direction in old coordinates.
 * \param toCoord : new coordinate system.
 * \param newPos : point in new coordinates.
 * \param newDir : direction in new coordinates.
 */
void TransCoordinates::coordTransf(enum_coordinate_type fromCoord,
    const vec4& oldPos, const vec4& oldDir,
    enum_coordinate_type toCoord, vec4& newPos,
    vec4& newDir)
{
    switch (fromCoord) {
    case (enum_coordinate_cartesian): {
        switch (toCoord) {
        case (enum_coordinate_spherical): {
            TransCoordinates::transCoordCartSph(oldPos, oldDir, newPos, newDir);
            break;
        };
        case (enum_coordinate_cylinder): {
            TransCoordinates::transCoordCartCyl(oldPos, oldDir, newPos, newDir);
            break;
        }
        case (enum_coordinate_cartesian): {
            newPos = oldPos;
            newDir = oldDir;
            break;
        }
        default:
            break;
        }
        break;
    }
    case (enum_coordinate_spherical): {
        switch (toCoord) {
        case (enum_coordinate_cartesian): {
            TransCoordinates::transCoordSphCart(oldPos, oldDir, newPos, newDir);
            break;
        }
        case (enum_coordinate_cylinder): {
            TransCoordinates::transCoordSphCart(oldPos, oldDir, newPos, newDir);
            TransCoordinates::transCoordCartCyl(oldPos, oldDir, newPos, newDir);
            break;
        }
        case (enum_coordinate_spherical): {
            newPos = oldPos;
            newDir = oldDir;
            break;
        }
        default:
            break;
        }
        break;
    }
    case (enum_coordinate_cylinder): {
        switch (toCoord) {
        case (enum_coordinate_cartesian): {
            TransCoordinates::transCoordCylCart(oldPos, oldDir, newPos, newDir);
            break;
        }
        case (enum_coordinate_spherical): {
            TransCoordinates::transCoordCylCart(oldPos, oldDir, newPos, newDir);
            TransCoordinates::transCoordCartSph(oldPos, oldDir, newPos, newDir);
            break;
        }
        case (enum_coordinate_cylinder): {
            newPos = oldPos;
            newDir = oldDir;
            break;
        }
        default:
            break;
        }
        break;
    }
    default:
        break;
    }
}

/*!
 * transform POINT between coordinate systems
 * \param fromCoord : old coordinate system
 * \param oldPos    : reference to point in old coordinates
 * \param oldDir0   : reference to direction in old coordinates
 * \param oldDir1   : reference to direction in old coordinates
 * \param oldDir2   : reference to direction in old coordinates
 * \param oldDir3   : reference to direction in old coordinates
 * \param toCoord   : new coordinate system
 * \param newPos    : reference to point in new coordinates
 * \param newDir0   : reference to direction in new coordinates
 * \param newDir1   : reference to direction in new coordinates
 * \param newDir2   : reference to direction in new coordinates
 * \param newDir3   : reference to direction in new coordinates
 */
void TransCoordinates::coordTransf(enum_coordinate_type fromCoord,
    const vec4& oldPos, const vec4& oldDir0,
    const vec4& oldDir1, const vec4& oldDir2,
    const vec4& oldDir3,
    enum_coordinate_type toCoord, vec4& newPos,
    vec4& newDir0, vec4& newDir1, vec4& newDir2,
    vec4& newDir3)
{
    switch (fromCoord) {
    case (enum_coordinate_cartesian): {
        switch (toCoord) {
        case (enum_coordinate_spherical): {
            TransCoordinates::transCoordCartSph(oldPos, oldDir0, oldDir1, oldDir2,
                oldDir3, newPos, newDir0, newDir1,
                newDir2, newDir3);
            break;
        }
        case (enum_coordinate_cylinder): {
            TransCoordinates::transCoordCartCyl(oldPos, oldDir0, oldDir1, oldDir2,
                oldDir3, newPos, newDir0, newDir1,
                newDir2, newDir3);
            break;
        }
        case (enum_coordinate_cartesian): {
            newPos = oldPos;
            newDir0 = oldDir0;
            newDir1 = oldDir1;
            newDir2 = oldDir2;
            newDir3 = oldDir3;
            break;
        }
        default:
            break;
        }
        break;
    }
    case (enum_coordinate_spherical): {
        switch (toCoord) {
        case (enum_coordinate_cartesian): {
            TransCoordinates::transCoordSphCart(oldPos, oldDir0, oldDir1, oldDir2,
                oldDir3, newPos, newDir0, newDir1,
                newDir2, newDir3);
            break;
        }
        case (enum_coordinate_cylinder): {
            TransCoordinates::transCoordSphCart(oldPos, oldDir0, oldDir1, oldDir2,
                oldDir3, newPos, newDir0, newDir1,
                newDir2, newDir3);
            TransCoordinates::transCoordCartCyl(oldPos, oldDir0, oldDir1, oldDir2,
                oldDir3, newPos, newDir0, newDir1,
                newDir2, newDir3);
            break;
        }
        case (enum_coordinate_spherical): {
            newPos = oldPos;
            newDir0 = oldDir0;
            newDir1 = oldDir1;
            newDir2 = oldDir2;
            newDir3 = oldDir3;
            break;
        }
        default:
            break;
        }
        break;
    }
    case (enum_coordinate_cylinder): {
        switch (toCoord) {
        case (enum_coordinate_cartesian): {
            TransCoordinates::transCoordCylCart(oldPos, oldDir0, oldDir1, oldDir2,
                oldDir3, newPos, newDir0, newDir1,
                newDir2, newDir3);
            break;
        }
        case (enum_coordinate_spherical): {
            TransCoordinates::transCoordCylCart(oldPos, oldDir0, oldDir1, oldDir2,
                oldDir3, newPos, newDir0, newDir1,
                newDir2, newDir3);
            TransCoordinates::transCoordCartSph(oldPos, oldDir0, oldDir1, oldDir2,
                oldDir3, newPos, newDir0, newDir1,
                newDir2, newDir3);
            break;
        }
        case (enum_coordinate_cylinder): {
            newPos = oldPos;
            newDir0 = oldDir0;
            newDir1 = oldDir1;
            newDir2 = oldDir2;
            newDir3 = oldDir3;
            break;
        }
        default:
            break;
        }
        break;
    }
    default:
        break;
    }
}

/*!
 * transform RAY between coordinate systems.
 * \param fromCoord : old coordinate system.
 * \param oldPos : points in old coordinates.
 * \param toCoord : new coordinate system.
 * \param newPos : points in new coordinates.
 */
void TransCoordinates::coordTransf(enum_coordinate_type fromCoord,
    const std::vector<vec4>& oldPos,
    enum_coordinate_type toCoord,
    std::vector<vec4>& newPos)
{
    unsigned int i = 0;

    // resize if necessary (initialization of newPos not required)
    if (newPos.size() != oldPos.size()) {
        newPos.clear();
        newPos.resize(oldPos.size());
    }

    // implementierung ohne aufruf von' toCartesianCoord (fromCoord,vec4& oldPos,
    // ...) damit nicht dauernd der switch abgearbeotet wird. sollte so schneller
    // sein.

    switch (fromCoord) {
    case (enum_coordinate_cartesian): {
        switch (toCoord) {
        case (enum_coordinate_spherical): {
            for (i = 0; i < oldPos.size(); i++) {
                TransCoordinates::transCoordCartSph(oldPos[i], newPos[i]);
            }
            break;
        }
        case (enum_coordinate_cylinder): {
            for (i = 0; i < oldPos.size(); i++) {
                TransCoordinates::transCoordCartCyl(oldPos[i], newPos[i]);
            }
            break;
        }
        case (enum_coordinate_cartesian): {
            for (i = 0; i < oldPos.size(); i++) {
                newPos[i] = oldPos[i];
            }
            break;
        }
        default:
            break;
        }
        break;
    }
    case (enum_coordinate_spherical): {
        switch (toCoord) {
        case (enum_coordinate_cartesian): {
            for (i = 0; i < oldPos.size(); i++) {
                TransCoordinates::transCoordSphCart(oldPos[i], newPos[i]);
            }
            break;
        }
        case (enum_coordinate_cylinder): {
            for (i = 0; i < oldPos.size(); i++) {
                TransCoordinates::transCoordSphCart(oldPos[i], newPos[i]);
                TransCoordinates::transCoordCartCyl(newPos[i], newPos[i]);
            }
            break;
        }
        case (enum_coordinate_spherical): {
            for (i = 0; i < oldPos.size(); i++) {
                newPos[i] = oldPos[i];
            }
            break;
        }
        default:
            break;
        }
        break;
    }
    case (enum_coordinate_cylinder): {
        switch (toCoord) {
        case (enum_coordinate_cartesian): {
            for (i = 0; i < oldPos.size(); i++) {
                TransCoordinates::transCoordCylCart(oldPos[i], newPos[i]);
            }
            break;
        }
        case (enum_coordinate_spherical): {
            for (i = 0; i < oldPos.size(); i++) {
                TransCoordinates::transCoordCylCart(oldPos[i], newPos[i]);
                TransCoordinates::transCoordCartSph(newPos[i], newPos[i]);
            }
            break;
        }
        case (enum_coordinate_cylinder): {
            for (i = 0; i < oldPos.size(); i++) {
                newPos[i] = oldPos[i];
            }
            break;
        }
        default:
            break;
        }
        break;
    }
    default:
        break;
    }
}

/*!
 * transform RAY between coordinate systems
 * \param fromCoord : old coordinate system
 * \param oldPos : points in old coordinates
 * \param oldDir : directions in old coordinates
 * \param toCoord : new coordinate system
 * \param newPos : points in new coordinates
 * \param newDir : directions in new coordinates
 */
void TransCoordinates::coordTransf(enum_coordinate_type fromCoord,
    const std::vector<vec4>& oldPos,
    const std::vector<vec4>& oldDir,
    enum_coordinate_type toCoord,
    std::vector<vec4>& newPos,
    std::vector<vec4>& newDir)
{
    unsigned int i = 0;

    // resize if necessary (initialization of newPos not required)
    if (newPos.size() != oldPos.size()) {
        newPos.clear();
        newPos.resize(oldPos.size());
    }
    if (newDir.size() != oldDir.size()) {
        newDir.clear();
        newDir.resize(oldDir.size());
    }

    // implementierung ohne aufruf von' toCartesianCoord (fromCoord,vec4& oldPos,
    // ...) damit nicht dauernd der switch abgearbeotet wird. sollte so schneller
    // sein.

    switch (fromCoord) {
    case (enum_coordinate_cartesian): {
        switch (toCoord) {
        case (enum_coordinate_spherical): {
            for (i = 0; i < oldPos.size(); i++) {
                TransCoordinates::transCoordCartSph(oldPos[i], oldDir[i], newPos[i],
                    newDir[i]);
            }
            break;
        }
        case (enum_coordinate_cylinder): {
            for (i = 0; i < oldPos.size(); i++) {
                TransCoordinates::transCoordCartCyl(oldPos[i], oldDir[i], newPos[i],
                    newDir[i]);
            }
            break;
        }
        case (enum_coordinate_cartesian): {
            for (i = 0; i < oldPos.size(); i++) {
                newPos[i] = oldPos[i];
                newDir[i] = oldDir[i];
            }
            break;
        }
        default:
            break;
        }
        break;
    }
    case (enum_coordinate_spherical): {
        switch (toCoord) {
        case (enum_coordinate_cartesian): {
            for (i = 0; i < oldPos.size(); i++) {
                TransCoordinates::transCoordSphCart(oldPos[i], oldDir[i], newPos[i],
                    newDir[i]);
            }
            break;
        }
        case (enum_coordinate_cylinder): {
            for (i = 0; i < oldPos.size(); i++) {
                TransCoordinates::transCoordSphCart(oldPos[i], oldDir[i], newPos[i],
                    newDir[i]);
                TransCoordinates::transCoordCartCyl(oldPos[i], oldDir[i], newPos[i],
                    newDir[i]);
            }
            break;
        }
        case (enum_coordinate_spherical): {
            for (i = 0; i < oldPos.size(); i++) {
                newPos[i] = oldPos[i];
                newDir[i] = oldDir[i];
            }
            break;
        }
        default:
            break;
        }
        break;
    }
    case (enum_coordinate_cylinder): {
        switch (toCoord) {
        case (enum_coordinate_cartesian): {
            for (i = 0; i < oldPos.size(); i++) {
                TransCoordinates::transCoordCylCart(oldPos[i], oldDir[i], newPos[i],
                    newDir[i]);
            }
            break;
        }
        case (enum_coordinate_spherical): {
            for (i = 0; i < oldPos.size(); i++) {
                TransCoordinates::transCoordCylCart(oldPos[i], oldDir[i], newPos[i],
                    newDir[i]);
                TransCoordinates::transCoordCartSph(oldPos[i], oldDir[i], newPos[i],
                    newDir[i]);
            }
            break;
        }
        case (enum_coordinate_cylinder): {
            for (i = 0; i < oldPos.size(); i++) {
                newPos[i] = oldPos[i];
                newDir[i] = oldDir[i];
            }
            break;
        }
        default:
            break;
        }
        break;
    }
    default:
        break;
    }
}

/**
 * transform RAY between coordinate systems
 * \param fromCoord  : old coordinate system
 * \param oldPos     : reference to points in old coordinates
 * \param oldDir0    : reference to directions in old coordinates
 * \param oldDir1    : reference to directions in old coordinates
 * \param oldDir2    : reference to directions in old coordinates
 * \param oldDir3    : reference to directions in old coordinates
 * \param toCoord    : new coordinate system
 * \param newPos     : reference to points in new coordinates
 * \param newDir0    : reference to directions in new coordinates
 * \param newDir1    : reference to directions in new coordinates
 * \param newDir2    : reference to directions in new coordinates
 * \param newDir3    : reference to directions in new coordinates
 */
void TransCoordinates::coordTransf(
    enum_coordinate_type fromCoord, const std::vector<vec4>& oldPos,
    const std::vector<vec4>& oldDir0, const std::vector<vec4>& oldDir1,
    const std::vector<vec4>& oldDir2, const std::vector<vec4>& oldDir3,
    enum_coordinate_type toCoord, std::vector<vec4>& newPos,
    std::vector<vec4>& newDir0, std::vector<vec4>& newDir1,
    std::vector<vec4>& newDir2, std::vector<vec4>& newDir3)
{
    unsigned int i = 0;

    // resize if necessary (initialization of newPos not required)
    if (newPos.size() != oldPos.size()) {
        newPos.clear();
        newPos.resize(oldPos.size());
    }
    if (newDir0.size() != oldDir0.size()) {
        newDir0.clear();
        newDir0.resize(oldDir0.size());
    }
    if (newDir1.size() != oldDir1.size()) {
        newDir1.clear();
        newDir1.resize(oldDir1.size());
    }
    if (newDir2.size() != oldDir2.size()) {
        newDir2.clear();
        newDir2.resize(oldDir2.size());
    }
    if (newDir3.size() != oldDir3.size()) {
        newDir3.clear();
        newDir3.resize(oldDir3.size());
    }

    // implementierung ohne aufruf von' toCartesianCoord (fromCoord,vec4& oldPos,
    // ...) damit nicht dauernd der switch abgearbeotet wird. sollte so schneller
    // sein.

    switch (fromCoord) {
    case (enum_coordinate_cartesian): {
        switch (toCoord) {
        case (enum_coordinate_spherical): {
            for (i = 0; i < oldPos.size(); i++)
                TransCoordinates::transCoordCartSph(
                    oldPos[i], oldDir0[i], oldDir1[i], oldDir2[i], oldDir3[i],
                    newPos[i], newDir0[i], newDir1[i], newDir2[i], newDir3[i]);
            break;
        }
        case (enum_coordinate_cylinder): {
            for (i = 0; i < oldPos.size(); i++)
                TransCoordinates::transCoordCartCyl(
                    oldPos[i], oldDir0[i], oldDir1[i], oldDir2[i], oldDir3[i],
                    newPos[i], newDir0[i], newDir1[i], newDir2[i], newDir3[i]);
            break;
        }
        case (enum_coordinate_cartesian): {
            for (i = 0; i < oldPos.size(); i++) {
                newPos[i] = oldPos[i];
                newDir0[i] = oldDir0[i];
                newDir1[i] = oldDir1[i];
                newDir2[i] = oldDir2[i];
                newDir3[i] = oldDir3[i];
            }
            break;
        }
        default:
            break;
        }
        break;
    }
    case (enum_coordinate_spherical): {
        switch (toCoord) {
        case (enum_coordinate_cartesian): {
            for (i = 0; i < oldPos.size(); i++)
                TransCoordinates::transCoordSphCart(
                    oldPos[i], oldDir0[i], oldDir1[i], oldDir2[i], oldDir3[i],
                    newPos[i], newDir0[i], newDir1[i], newDir2[i], newDir3[i]);
            break;
        }
        case (enum_coordinate_cylinder): {
            for (i = 0; i < oldPos.size(); i++) {
                TransCoordinates::transCoordSphCart(
                    oldPos[i], oldDir0[i], oldDir1[i], oldDir2[i], oldDir3[i],
                    newPos[i], newDir0[i], newDir1[i], newDir2[i], newDir3[i]);
                TransCoordinates::transCoordCartCyl(
                    oldPos[i], oldDir0[i], oldDir1[i], oldDir2[i], oldDir3[i],
                    newPos[i], newDir0[i], newDir1[i], newDir2[i], newDir3[i]);
            }
            break;
        }
        case (enum_coordinate_spherical): {
            for (i = 0; i < oldPos.size(); i++) {
                newPos[i] = oldPos[i];
                newDir0[i] = oldDir0[i];
                newDir1[i] = oldDir1[i];
                newDir2[i] = oldDir2[i];
                newDir3[i] = oldDir3[i];
            }
            break;
        }
        default:
            break;
        }
        break;
    }
    case (enum_coordinate_cylinder): {
        switch (toCoord) {
        case (enum_coordinate_cartesian): {
            for (i = 0; i < oldPos.size(); i++)
                TransCoordinates::transCoordCylCart(
                    oldPos[i], oldDir0[i], oldDir1[i], oldDir2[i], oldDir3[i],
                    newPos[i], newDir0[i], newDir1[i], newDir2[i], newDir3[i]);
            break;
        }
        case (enum_coordinate_spherical): {
            for (i = 0; i < oldPos.size(); i++) {
                TransCoordinates::transCoordCylCart(
                    oldPos[i], oldDir0[i], oldDir1[i], oldDir2[i], oldDir3[i],
                    newPos[i], newDir0[i], newDir1[i], newDir2[i], newDir3[i]);
                TransCoordinates::transCoordCartSph(
                    oldPos[i], oldDir0[i], oldDir1[i], oldDir2[i], oldDir3[i],
                    newPos[i], newDir0[i], newDir1[i], newDir2[i], newDir3[i]);
            }
            break;
        }
        case (enum_coordinate_cylinder): {
            for (i = 0; i < oldPos.size(); i++) {
                newPos[i] = oldPos[i];
                newDir0[i] = oldDir0[i];
                newDir1[i] = oldDir1[i];
                newDir2[i] = oldDir2[i];
                newDir3[i] = oldDir3[i];
            }
            break;
        }
        default:
            break;
        }
        break;
    }
    default:
        break;
    }
}

/*!
 * single coordinate transformation from cartesian to spherical coordinates.
 * \param oldPos : point in cartesian coordinates.
 * \param newPos : point in spherical coordinates.
 */
void TransCoordinates::transCoordCartSph(const vec4& oldPos, vec4& newPos)
{
    //  cerr << "TransCoordinates::transCartSph...\n";
    double t = oldPos[0];
    double x = oldPos[1];
    double y = oldPos[2];
    double z = oldPos[3];
    double r = sqrt(x * x + y * y + z * z);

    double theta = acos(z / r);
    double phi = atan2(y, x);
    newPos = vec4(t, r, theta, phi);
}

/*!
 * single coordinate transformation from cartesian to spherical coordinates.
 * \param oldPos : point in cartesian coordinates.
 * \param oldDir : direction in cartesian coordinates.
 * \param newPos : point in spherical coordinates.
 * \param newDir : direction in spherical coordinates.
 */
void TransCoordinates::transCoordCartSph(const vec4& oldPos, const vec4& oldDir,
    vec4& newPos, vec4& newDir)
{
    transCoordCartSph(oldPos, newPos);

    double r = newPos[1];
    double theta = newPos[2];
    double phi = newPos[3];

    double dirX = oldDir[1];
    double dirY = oldDir[2];
    double dirZ = oldDir[3];

    newDir[0] = oldDir[0];
    newDir[1] = dirX * sin(theta) * cos(phi) + dirY * sin(theta) * sin(phi) + dirZ * cos(theta); // dirR
    newDir[2] = dirX * cos(theta) * cos(phi) / r + dirY * cos(theta) * sin(phi) / r - dirZ * sin(theta) / r; // dirTheta
    newDir[3] = -dirX * sin(phi) / r / sin(theta) + dirY * cos(phi) / r / sin(theta); // dirPhi
}

/*!
 * single coordinate transformation from cartesian to spherical coordinates.
 * \param oldPos : point in cartesian coordinates.
 * \param oldDir0 : direction in cartesian coordinates.
 * \param oldDir1 : direction in cartesian coordinates.
 * \param oldDir2 : direction in cartesian coordinates.
 * \param oldDir3 : direction in cartesian coordinates.
 * \param newPos : point in spherical coordinates.
 * \param newDir0 : direction in spherical coordinates.
 * \param newDir1 : direction in spherical coordinates.
 * \param newDir2 : direction in spherical coordinates.
 * \param newDir3 : direction in spherical coordinates.
 */
void TransCoordinates::transCoordCartSph(
    const vec4& oldPos, const vec4& oldDir0, const vec4& oldDir1,
    const vec4& oldDir2, const vec4& oldDir3, vec4& newPos, vec4& newDir0,
    vec4& newDir1, vec4& newDir2, vec4& newDir3)
{
    transCoordCartSph(oldPos, newPos);

    double r = newPos[1];
    double theta = newPos[2];
    double phi = newPos[3];

    double dirX0 = oldDir0[1];
    double dirY0 = oldDir0[2];
    double dirZ0 = oldDir0[3];
    double dirX1 = oldDir1[1];
    double dirY1 = oldDir1[2];
    double dirZ1 = oldDir1[3];
    double dirX2 = oldDir2[1];
    double dirY2 = oldDir2[2];
    double dirZ2 = oldDir2[3];
    double dirX3 = oldDir3[1];
    double dirY3 = oldDir3[2];
    double dirZ3 = oldDir3[3];

    newDir0[0] = oldDir0[0];
    newDir0[1] = dirX0 * sin(theta) * cos(phi) + dirY0 * sin(theta) * sin(phi) + dirZ0 * cos(theta); // dirR
    newDir0[2] = dirX0 * cos(theta) * cos(phi) / r + dirY0 * cos(theta) * sin(phi) / r - dirZ0 * sin(theta) / r; // dirTheta
    newDir0[3] = -dirX0 * sin(phi) / r / sin(theta) + dirY0 * cos(phi) / r / sin(theta); // dirPhi
    newDir1[0] = oldDir1[0];
    newDir1[1] = dirX1 * sin(theta) * cos(phi) + dirY1 * sin(theta) * sin(phi) + dirZ1 * cos(theta); // dirR
    newDir1[2] = dirX1 * cos(theta) * cos(phi) / r + dirY1 * cos(theta) * sin(phi) / r - dirZ1 * sin(theta) / r; // dirTheta
    newDir1[3] = -dirX1 * sin(phi) / r / sin(theta) + dirY1 * cos(phi) / r / sin(theta); // dirPhi
    newDir2[0] = oldDir2[0];
    newDir2[1] = dirX2 * sin(theta) * cos(phi) + dirY2 * sin(theta) * sin(phi) + dirZ2 * cos(theta); // dirR
    newDir2[2] = dirX2 * cos(theta) * cos(phi) / r + dirY2 * cos(theta) * sin(phi) / r - dirZ2 * sin(theta) / r; // dirTheta
    newDir2[3] = -dirX2 * sin(phi) / r / sin(theta) + dirY2 * cos(phi) / r / sin(theta); // dirPhi
    newDir3[0] = oldDir3[0];
    newDir3[1] = dirX3 * sin(theta) * cos(phi) + dirY3 * sin(theta) * sin(phi) + dirZ3 * cos(theta); // dirR
    newDir3[2] = dirX3 * cos(theta) * cos(phi) / r + dirY3 * cos(theta) * sin(phi) / r - dirZ3 * sin(theta) / r; // dirTheta
    newDir3[3] = -dirX3 * sin(phi) / r / sin(theta) + dirY3 * cos(phi) / r / sin(theta); // dirPhi
}

/*!
 * single coordinate transformation from spherical to cartesian coordinates.
 * \param oldPos : point in spherical coordinates.
 * \param newPos : point in cartesian coordinates.
 */
void TransCoordinates::transCoordSphCart(const vec4& oldPos, vec4& newPos)
{
    double t = oldPos[0];
    double r = oldPos[1];
    double theta = oldPos[2];
    double phi = oldPos[3];

    newPos = vec4(t, r * sin(theta) * cos(phi), r * sin(theta) * sin(phi),
        r * cos(theta));
}

/*!
 * single coordinate transformation from spherical to cartesian coordinates.
 * \param oldPos : point in spherical coordinates.
 * \param oldDir : direction in spherical coordinates.
 * \param newPos : point in cartesian coordinates.
 * \param newDir : direction in cartesian coordinates.
 */
void TransCoordinates::transCoordSphCart(const vec4& oldPos, const vec4& oldDir,
    vec4& newPos, vec4& newDir)
{
    double r = oldPos[1];
    double theta = oldPos[2];
    double phi = oldPos[3];

    transCoordSphCart(oldPos, newPos);

    double dirR = oldDir[1];
    double dirTheta = oldDir[2];
    double dirPhi = oldDir[3];

    newDir[0] = oldDir[0];
    newDir[1] = dirR * sin(theta) * cos(phi) + dirTheta * r * cos(theta) * cos(phi) - dirPhi * r * sin(theta) * sin(phi); // dirX
    newDir[2] = dirR * sin(theta) * sin(phi) + dirTheta * r * cos(theta) * sin(phi) + dirPhi * r * sin(theta) * cos(phi); // dirY
    newDir[3] = dirR * cos(theta) - dirTheta * r * sin(theta); // dirZ
}

/*!
 * single coordinate transformation from spherical to cartesian coordinates.
 * \param oldPos : point in spherical coordinates.
 * \param oldDir0 : direction in spherical coordinates.
 * \param oldDir1 : direction in spherical coordinates.
 * \param oldDir2 : direction in spherical coordinates.
 * \param oldDir3 : direction in spherical coordinates.
 * \param newPos : point in cartesian coordinates.
 * \param newDir0 : direction in cartesian coordinates.
 * \param newDir1 : direction in cartesian coordinates.
 * \param newDir2 : direction in cartesian coordinates.
 * \param newDir3 : direction in cartesian coordinates.
 */
void TransCoordinates::transCoordSphCart(
    const vec4& oldPos, const vec4& oldDir0, const vec4& oldDir1,
    const vec4& oldDir2, const vec4& oldDir3, vec4& newPos, vec4& newDir0,
    vec4& newDir1, vec4& newDir2, vec4& newDir3)
{
    double r = oldPos[1];
    double theta = oldPos[2];
    double phi = oldPos[3];

    transCoordSphCart(oldPos, newPos);

    double dirR0 = oldDir0[1];
    double dirTheta0 = oldDir0[2];
    double dirPhi0 = oldDir0[3];
    double dirR1 = oldDir1[1];
    double dirTheta1 = oldDir1[2];
    double dirPhi1 = oldDir1[3];
    double dirR2 = oldDir2[1];
    double dirTheta2 = oldDir2[2];
    double dirPhi2 = oldDir2[3];
    double dirR3 = oldDir3[1];
    double dirTheta3 = oldDir3[2];
    double dirPhi3 = oldDir3[3];

    newDir0[0] = oldDir0[0];
    newDir0[1] = dirR0 * sin(theta) * cos(phi) + dirTheta0 * r * cos(theta) * cos(phi) - dirPhi0 * r * sin(theta) * sin(phi); // dirX
    newDir0[2] = dirR0 * sin(theta) * sin(phi) + dirTheta0 * r * cos(theta) * sin(phi) + dirPhi0 * r * sin(theta) * cos(phi); // dirY
    newDir0[3] = dirR0 * cos(theta) - dirTheta0 * r * sin(theta); // dirZ
    newDir1[0] = oldDir1[0];
    newDir1[1] = dirR1 * sin(theta) * cos(phi) + dirTheta1 * r * cos(theta) * cos(phi) - dirPhi1 * r * sin(theta) * sin(phi); // dirX
    newDir1[2] = dirR1 * sin(theta) * sin(phi) + dirTheta1 * r * cos(theta) * sin(phi) + dirPhi1 * r * sin(theta) * cos(phi); // dirY
    newDir1[3] = dirR1 * cos(theta) - dirTheta1 * r * sin(theta); // dirZ
    newDir2[0] = oldDir2[0];
    newDir2[1] = dirR2 * sin(theta) * cos(phi) + dirTheta2 * r * cos(theta) * cos(phi) - dirPhi2 * r * sin(theta) * sin(phi); // dirX
    newDir2[2] = dirR2 * sin(theta) * sin(phi) + dirTheta2 * r * cos(theta) * sin(phi) + dirPhi2 * r * sin(theta) * cos(phi); // dirY
    newDir2[3] = dirR2 * cos(theta) - dirTheta2 * r * sin(theta); // dirZ
    newDir3[0] = oldDir3[0];
    newDir3[1] = dirR3 * sin(theta) * cos(phi) + dirTheta3 * r * cos(theta) * cos(phi) - dirPhi3 * r * sin(theta) * sin(phi); // dirX
    newDir3[2] = dirR3 * sin(theta) * sin(phi) + dirTheta3 * r * cos(theta) * sin(phi) + dirPhi3 * r * sin(theta) * cos(phi); // dirY
    newDir3[3] = dirR3 * cos(theta) - dirTheta3 * r * sin(theta); // dirZ
}

/*!
 * single coordinate transformation from cartesian to cylindrical coordinates.
 * \param oldPos : point in cartesian coordinates.
 * \param newPos : point in cylindrical coordinates.
 */
void TransCoordinates::transCoordCartCyl(const vec4& oldPos, vec4& newPos)
{
    double t = oldPos[0];
    double x = oldPos[1];
    double y = oldPos[2];
    double z = oldPos[3];
    double r = sqrt(x * x + y * y);

    // Winkel phi
    double phi = 0.0, tan_phi = 0.0;
    if (x > 0.0 || x < 0.0) {
        tan_phi = y / x;
        phi = 0.0;
    }

    // ????

    if (y > 0.0 && x > 0.0) {
        phi = atan(tan_phi);
    } else if (y > 0.0 && x < 0.0) {
        phi = M_PI + atan(tan_phi);
    } else if (y < 0.0 && x < 0.0) {
        phi = M_PI + atan(tan_phi);
    } else if (y < 0.0 && x > 0.0) {
        phi = 2.0 * M_PI + atan(tan_phi);
    } else if (x == 0.0 && y > 0.0) {
        phi = 0.5 * M_PI;
    } else if (x == 0.0 && y < 0.0) {
        phi = 1.5 * M_PI;
    } else if (x > 0.0 && y == 0.0) {
        phi = 0.0;
    } else if (x < 0.0 && y == 0.0) {
        phi = M_PI;
    }

    newPos = vec4(t, r, phi, z);
}

/*!
 * single coordinate transformation from cartesian to cylindrical coordinates.
 * \param oldPos : point in cartesian coordinates.
 * \param oldDir : direction in cartesian coordinates.
 * \param newPos : point in cylindrical coordinates.
 * \param newDir : direction in cylindrical coordinates.
 */
void TransCoordinates::transCoordCartCyl(const vec4& oldPos, const vec4& oldDir,
    vec4& newPos, vec4& newDir)
{
    transCoordCartCyl(oldPos, newPos);
    double r = newPos[1];
    double phi = newPos[2];

    double dirX = oldDir[1];
    double dirY = oldDir[2];
    double dirZ = oldDir[3];

    newDir[0] = oldDir.x(0);
    newDir[1] = dirX * cos(phi) + dirY * sin(phi); // dirR
    newDir[2] = -dirX * sin(phi) / r + dirY * cos(phi) / r; // dirPhi
    newDir[3] = dirZ;
}

/*!
 * single coordinate transformation from cartesian to cylindrical coordinates.
 * \param oldPos : point in cartesian coordinates.
 * \param oldDir0 : direction in cartesian coordinates.
 * \param oldDir1 : direction in cartesian coordinates.
 * \param oldDir2 : direction in cartesian coordinates.
 * \param oldDir3 : direction in cartesian coordinates.
 * \param newPos : point in cylindrical coordinates.
 * \param newDir0 : direction in cylindrical coordinates.
 * \param newDir1 : direction in cylindrical coordinates.
 * \param newDir2 : direction in cylindrical coordinates.
 * \param newDir3 : direction in cylindrical coordinates.
 */
void TransCoordinates::transCoordCartCyl(
    const vec4& oldPos, const vec4& oldDir0, const vec4& oldDir1,
    const vec4& oldDir2, const vec4& oldDir3, vec4& newPos, vec4& newDir0,
    vec4& newDir1, vec4& newDir2, vec4& newDir3)
{
    transCoordCartCyl(oldPos, newPos);
    double r = newPos[1];
    double phi = newPos[2];

    double dirX0 = oldDir0[1];
    double dirY0 = oldDir0[2];
    double dirZ0 = oldDir0[3];
    double dirX1 = oldDir1[1];
    double dirY1 = oldDir1[2];
    double dirZ1 = oldDir1[3];
    double dirX2 = oldDir2[1];
    double dirY2 = oldDir2[2];
    double dirZ2 = oldDir2[3];
    double dirX3 = oldDir3[1];
    double dirY3 = oldDir3[2];
    double dirZ3 = oldDir3[3];

    newDir0[0] = oldDir0.x(0);
    newDir0[1] = dirX0 * cos(phi) + dirY0 * sin(phi); // dirR
    newDir0[2] = -dirX0 * sin(phi) / r + dirY0 * cos(phi) / r; // dirPhi
    newDir0[3] = dirZ0;
    newDir1[0] = oldDir1.x(0);
    newDir1[1] = dirX1 * cos(phi) + dirY1 * sin(phi); // dirR
    newDir1[2] = -dirX1 * sin(phi) / r + dirY1 * cos(phi) / r; // dirPhi
    newDir1[3] = dirZ1;
    newDir2[0] = oldDir2.x(0);
    newDir2[1] = dirX2 * cos(phi) + dirY2 * sin(phi); // dirR
    newDir2[2] = -dirX2 * sin(phi) / r + dirY2 * cos(phi) / r; // dirPhi
    newDir2[3] = dirZ2;
    newDir3[0] = oldDir3.x(0);
    newDir3[1] = dirX3 * cos(phi) + dirY3 * sin(phi); // dirR
    newDir3[2] = -dirX3 * sin(phi) / r + dirY3 * cos(phi) / r; // dirPhi
    newDir3[3] = dirZ3;
}

/*!
 * single coordinate transformation from cylindrical to cartesian coordinates.
 * \param oldPos : point in cylindrical coordinates.
 * \param newPos : point in cartesian coordinates.
 */
void TransCoordinates::transCoordCylCart(const vec4& oldPos, vec4& newPos)
{
    double t = oldPos[0];
    double r = oldPos[1];
    double phi = oldPos[2];
    double z = oldPos[3];

    newPos = vec4(t, r * cos(phi), r * sin(phi), z);
}

/*!
 * single coordinate transformation from cylindrical to cartesian coordinates.
 * \param oldPos : point in cylindrical coordinates.
 * \param oldDir : direction in cylindrical coordinates.
 * \param newPos : point in cartesian coordinates.
 * \param newDir : direction in cartesian coordinates.
 */
void TransCoordinates::transCoordCylCart(const vec4& oldPos, const vec4& oldDir,
    vec4& newPos, vec4& newDir)
{
    double r = oldPos[1];
    double phi = oldPos[2];
    transCoordCylCart(oldPos, newPos);

    double dirR = oldDir[1];
    double dirPhi = oldDir[2];

    newDir[0] = oldDir.x(0);
    newDir[1] = dirR * cos(phi) - dirPhi * r * sin(phi); // dirX
    newDir[2] = dirR * sin(phi) + dirPhi * r * cos(phi); // dirY
    newDir[3] = oldDir.x(3);
}

/**
 * single coordinate transformation from cylindrical to cartesian coordinates.
 * \param oldPos : point in cylindrical coordinates.
 * \param oldDir0 : direction in cylindrical coordinates.
 * \param oldDir1 : direction in cylindrical coordinates.
 * \param oldDir2 : direction in cylindrical coordinates.
 * \param oldDir3 : direction in cylindrical coordinates.
 * \param newPos : point in cartesian coordinates.
 * \param newDir0 : direction in cartesian coordinates.
 * \param newDir1 : direction in cartesian coordinates.
 * \param newDir2 : direction in cartesian coordinates.
 * \param newDir3 : direction in cartesian coordinates.
 */
void TransCoordinates::transCoordCylCart(
    const vec4& oldPos, const vec4& oldDir0, const vec4& oldDir1,
    const vec4& oldDir2, const vec4& oldDir3, vec4& newPos, vec4& newDir0,
    vec4& newDir1, vec4& newDir2, vec4& newDir3)
{
    double r = oldPos[1];
    double phi = oldPos[2];
    transCoordCylCart(oldPos, newPos);

    double dirR0 = oldDir0[1];
    double dirPhi0 = oldDir0[2];
    double dirR1 = oldDir1[1];
    double dirPhi1 = oldDir1[2];
    double dirR2 = oldDir2[1];
    double dirPhi2 = oldDir2[2];
    double dirR3 = oldDir3[1];
    double dirPhi3 = oldDir3[2];

    newDir0[0] = oldDir0.x(0);
    newDir0[1] = dirR0 * cos(phi) - dirPhi0 * r * sin(phi); // dirX
    newDir0[2] = dirR0 * sin(phi) + dirPhi0 * r * cos(phi); // dirY
    newDir0[3] = oldDir0.x(3);
    newDir1[0] = oldDir1.x(0);
    newDir1[1] = dirR1 * cos(phi) - dirPhi1 * r * sin(phi); // dirX
    newDir1[2] = dirR1 * sin(phi) + dirPhi1 * r * cos(phi); // dirY
    newDir1[3] = oldDir1.x(3);
    newDir2[0] = oldDir2.x(0);
    newDir2[1] = dirR2 * cos(phi) - dirPhi2 * r * sin(phi); // dirX
    newDir2[2] = dirR2 * sin(phi) + dirPhi2 * r * cos(phi); // dirY
    newDir2[3] = oldDir2.x(3);
    newDir3[0] = oldDir3.x(0);
    newDir3[1] = dirR3 * cos(phi) - dirPhi3 * r * sin(phi); // dirX
    newDir3[2] = dirR3 * sin(phi) + dirPhi3 * r * cos(phi); // dirY
    newDir3[3] = oldDir3.x(3);
}

/*!
 * single coordinate transformation from cylindrical to cartesian coordinates.
 * \param oldPos : point in cylindrical coordinates.
 * \param newPos : point in cartesian coordinates.
 * \param a : Free parameter in prolate spheroidal coordinates (similar to a
 * radius)
 */
void TransCoordinates::transCoordProlSphCart(const vec4& oldPos, vec4& newPos,
    double a)
{
    double t = oldPos[0];
    double sigma = oldPos[1];
    double tau = oldPos[2];
    double phi = oldPos[3];

    assert(sigma >= 1.0);
    assert(std::abs(tau) <= 1.0);

    double rho = a * sqrt((sigma * sigma - 1.0) * (1.0 - tau * tau));

    newPos = vec4(t, rho * cos(phi), rho * sin(phi), a * sigma * tau);
}

/*!
 * single coordinate transformation from cylindrical to cartesian coordinates.
 * \param oldPos : point in cylindrical coordinates.
 * \param oldDir : direction in cylindrical coordinates.
 * \param newPos : point in cartesian coordinates.
 * \param newDir : direction in cartesian coordinates.
 */
void TransCoordinates::transCoordProlSphCart(const vec4& oldPos,
    const vec4& oldDir, vec4& newPos,
    vec4& newDir)
{
    double r = oldPos[1];
    double phi = oldPos[2];
    transCoordCylCart(oldPos, newPos);

    double dirR = oldDir[1];
    double dirPhi = oldDir[2];

    newDir[0] = oldDir.x(0);
    newDir[1] = dirR * cos(phi) - dirPhi * r * sin(phi); // dirX
    newDir[2] = dirR * sin(phi) + dirPhi * r * cos(phi); // dirY
    newDir[3] = oldDir.x(3);
}

/**
 * single coordinate transformation from cylindrical to cartesian coordinates.
 * \param oldPos : point in cylindrical coordinates.
 * \param oldDir0 : direction in cylindrical coordinates.
 * \param oldDir1 : direction in cylindrical coordinates.
 * \param oldDir2 : direction in cylindrical coordinates.
 * \param oldDir3 : direction in cylindrical coordinates.
 * \param newPos : point in cartesian coordinates.
 * \param newDir0 : direction in cartesian coordinates.
 * \param newDir1 : direction in cartesian coordinates.
 * \param newDir2 : direction in cartesian coordinates.
 * \param newDir3 : direction in cartesian coordinates.
 */
void TransCoordinates::transCoordProlSphCart(
    const vec4& oldPos, const vec4& oldDir0, const vec4& oldDir1,
    const vec4& oldDir2, const vec4& oldDir3, vec4& newPos, vec4& newDir0,
    vec4& newDir1, vec4& newDir2, vec4& newDir3)
{
    double r = oldPos[1];
    double phi = oldPos[2];
    transCoordCylCart(oldPos, newPos);

    double dirR0 = oldDir0[1];
    double dirPhi0 = oldDir0[2];
    double dirR1 = oldDir1[1];
    double dirPhi1 = oldDir1[2];
    double dirR2 = oldDir2[1];
    double dirPhi2 = oldDir2[2];
    double dirR3 = oldDir3[1];
    double dirPhi3 = oldDir3[2];

    newDir0[0] = oldDir0.x(0);
    newDir0[1] = dirR0 * cos(phi) - dirPhi0 * r * sin(phi); // dirX
    newDir0[2] = dirR0 * sin(phi) + dirPhi0 * r * cos(phi); // dirY
    newDir0[3] = oldDir0.x(3);
    newDir1[0] = oldDir1.x(0);
    newDir1[1] = dirR1 * cos(phi) - dirPhi1 * r * sin(phi); // dirX
    newDir1[2] = dirR1 * sin(phi) + dirPhi1 * r * cos(phi); // dirY
    newDir1[3] = oldDir1.x(3);
    newDir2[0] = oldDir2.x(0);
    newDir2[1] = dirR2 * cos(phi) - dirPhi2 * r * sin(phi); // dirX
    newDir2[2] = dirR2 * sin(phi) + dirPhi2 * r * cos(phi); // dirY
    newDir2[3] = oldDir2.x(3);
    newDir3[0] = oldDir3.x(0);
    newDir3[1] = dirR3 * cos(phi) - dirPhi3 * r * sin(phi); // dirX
    newDir3[2] = dirR3 * sin(phi) + dirPhi3 * r * cos(phi); // dirY
    newDir3[3] = oldDir3.x(3);
}

} // end namespace m4d
