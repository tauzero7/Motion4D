/**
 * @file    TransfMat.h
 * @author  Thomas Mueller
 *
 * @brief  Translation matrix in 3 dimensions.

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

#ifndef M4D_TRANSF_MAT_H
#define M4D_TRANSF_MAT_H

#include <cassert>
#include <string>
#include <typeinfo>

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
class API_M4D_EXPORT TranslateMat3D : public Matrix<double, 3, 4>
{
public:
    TranslateMat3D(double tx, double ty, double tz);
    explicit TranslateMat3D(const vec3& translat);
};

//----------------------------------------------------------------------------
//         RotateMat3D
//
//  angles in radiant!
//----------------------------------------------------------------------------
class API_M4D_EXPORT RotateMat3D : public Matrix<double, 3, 4>
{
public:
    RotateMat3D(const vec3& rotAxis, double rotAngle);
    RotateMat3D(enum_axisID mainAxis, double rotAngle);
};

class API_M4D_EXPORT RotateMat3Dd : public Matrix<double, 3, 3>
{
public:
    RotateMat3Dd();
    RotateMat3Dd(const vec3& rotAxis, double rotAngle);
    RotateMat3Dd(enum_axisID mainAxis, double rotAngle);
};

class API_M4D_EXPORT RotateMat3Df : public Matrix<float, 3, 3>
{
public:
    RotateMat3Df();
    RotateMat3Df(const vec3f& rotAxis, float rotAngle);
    RotateMat3Df(enum_axisID mainAxis, float rotAngle);
};

//----------------------------------------------------------------------------
//         ScaleMat3D
//----------------------------------------------------------------------------
//
//                ( sx  0  0  0  )
//  ScaleMat3D =  (  0 sy  0  0  )
//                (  0  0 sz  0  )
//
class API_M4D_EXPORT ScaleMat3D : public Matrix<double, 3, 4>
{
public:
    ScaleMat3D(double sx, double sy, double sz);
    explicit ScaleMat3D(const vec3& scale);
};

//----------------------------------------------------------------------------
//         TranslateMat2D
//----------------------------------------------------------------------------
class API_M4D_EXPORT TranslateMat2D : public Matrix<double, 2, 3>
{
public:
    TranslateMat2D(double tx, double ty);
};

//----------------------------------------------------------------------------
//         RotateMat2D
//----------------------------------------------------------------------------
class API_M4D_EXPORT RotateMat2D : public Matrix<double, 2, 3>
{
public:
    explicit RotateMat2D(double rotAngle);
    RotateMat2D(double rotCenterX, double rotCenterY, double rotAngle);
};

//----------------------------------------------------------------------------
//         ScaleMat2D
//----------------------------------------------------------------------------
class API_M4D_EXPORT ScaleMat2D : public Matrix<double, 2, 3>
{
public:
    ScaleMat2D(double sx, double sy);
};

//----------------------------------------------------------------------------
//         Lorentz matrix
//----------------------------------------------------------------------------
class API_M4D_EXPORT LorentzTransf : public Matrix<double, 4, 4>
{
public:
    LorentzTransf(double beta, VnD<double, 3> n);
};

} // end namespace m4d

#endif
