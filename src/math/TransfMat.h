// --------------------------------------------------------------------------------
/*
    TransfMat.h

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

/*!
     \class  m4d::TranslateMat3D
     \brief  Translation matrix in 3 dimensions.

                    TranslateMat3D ( float, float, float )
                                   ( vec3 )

     \class  m4d::TranslateMat2D
     \brief  Translation matrix in 2 dimensions.


     \class  m4d::RotateMat3D
     \brief  Rotation matrix in 3 dimensions as 3x4 matrix

                    RotateMatd3D   ( vec3, float )

     \class  m4d::RotateMat3Df
     \brief  Rotation matrix in 3 dimensions as 3x3 matrix.

     \class  m4d::RotateMat2D
     \brief  Rotation matrix in 2 dimensions.

     \class  m4d::ScaleMat3D
     \brief  Scaling matrix in 3 dimensions.

                    ScaleMat3D ( float float float )
                               ( vec3 )

     \class  m4d::ScaleMat2D
     \brief  Scaling matrix in 2 dimensions.

           Transformation matrices have to be Matrix<float,3,4> or
           Matrix<float,2,3> !
*/
// --------------------------------------------------------------------------------

#ifndef M4D_TRANSF_MAT_H
#define M4D_TRANSF_MAT_H

#include <string>
#include <typeinfo>
#include <cassert>

#include <m4dGlobalDefs.h>

namespace m4d {

//----------------------------------------------------------------------------
//         TranslateMat3D
//----------------------------------------------------------------------------
//
//                    (  1  0  0  tx )
//  TranslateMat3D =  (  0  1  0  ty )
//                    (  0  0  0  tz )
//
class MATH_API TranslateMat3D : public Matrix<double,3,4> {
public:
    TranslateMat3D(double tx, double ty, double tz);
    TranslateMat3D(const vec3 &translat);
};

//----------------------------------------------------------------------------
//         RotateMat3D
//
//  angles in radiant!
//----------------------------------------------------------------------------
class MATH_API RotateMat3D : public Matrix<double,3,4> {
public:
    RotateMat3D(const vec3 &rotAxis, double rotAngle);
    RotateMat3D(enum_axisID  mainAxis, double rotAngle);
};

class MATH_API RotateMat3Dd : public Matrix<double,3,3> {
public:
    RotateMat3Dd();
    RotateMat3Dd(const vec3  &rotAxis,  double rotAngle);
    RotateMat3Dd(enum_axisID  mainAxis, double rotAngle);
};


class MATH_API RotateMat3Df : public Matrix<float,3,3> {
public:
    RotateMat3Df();
    RotateMat3Df(const vec3f &rotAxis,  float rotAngle);
    RotateMat3Df(enum_axisID  mainAxis, float rotAngle);
};


//----------------------------------------------------------------------------
//         ScaleMat3D
//----------------------------------------------------------------------------
//
//                ( sx  0  0  0  )
//  ScaleMat3D =  (  0 sy  0  0  )
//                (  0  0 sz  0  )
//
class MATH_API ScaleMat3D : public Matrix<double,3,4> {
public:
    ScaleMat3D(double sx, double sy, double sz);
    ScaleMat3D(const vec3 &scale);
};


//----------------------------------------------------------------------------
//         TranslateMat2D
//----------------------------------------------------------------------------
class MATH_API TranslateMat2D : public Matrix<double,2,3> {
public:
    TranslateMat2D(double tx, double ty);
};

//----------------------------------------------------------------------------
//         RotateMat2D
//----------------------------------------------------------------------------
class MATH_API RotateMat2D : public Matrix<double,2,3> {
public:
    RotateMat2D(double rotAngle);
    RotateMat2D(double rotCenterX, double rotCenterY, double rotAngle);
};

//----------------------------------------------------------------------------
//         ScaleMat2D
//----------------------------------------------------------------------------
class MATH_API ScaleMat2D : public Matrix<double,2,3> {
public:
    ScaleMat2D(double sx, double sy);
};

//----------------------------------------------------------------------------
//         Lorentz matrix
//----------------------------------------------------------------------------
class MATH_API LorentzTransf : public Matrix<double,4,4> {
public:
    LorentzTransf(double beta, VnD<double,3> n);
};

} // end namespace m4d

#endif

