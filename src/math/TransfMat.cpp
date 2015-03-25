// -------------------------------------------------------------------------------
/*
    TransfMat.cpp

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

#include "TransfMat.h"

#ifdef _WIN32
#ifndef __GNUC__
#pragma warning (disable: 4244 )
#endif
#endif


namespace m4d {

//----------------------------------------------------------------------------
//     TranslateMat3D:    constructor, destructor
//----------------------------------------------------------------------------
TranslateMat3D::TranslateMat3D(double tx, double ty, double tz) {
    setIdent();

    mat[0][3] = tx;
    mat[1][3] = ty;
    mat[2][3] = tz;

    matIsUnit = (tx == 0.0) && (ty == 0.0) && (tz == 0.0);
}

TranslateMat3D::TranslateMat3D(const vec3& translat) {
    setIdent();

    mat[0][3] = translat[0];
    mat[1][3] = translat[1];
    mat[2][3] = translat[2];

    matIsUnit = (mat[0][3] == 0.0) && (mat[1][3] == 0.0) && (mat[2][3] == 0.0);
}

//----------------------------------------------------------------------------
//     RotateMat3D:    constructor, destructor
//----------------------------------------------------------------------------
//   algorithm by:  Michael E. Pique, Rotation Tools, pp. 465 - 469
//             in:  A. Glassner (ed.), Graphics Gems, Academic Press 1990
//-------------------------------------------------------------------------
RotateMat3D::RotateMat3D(const vec3 &rotAxis, double rotAngle) {
    vec3 axis = rotAxis.getNormalized();

    double sinAngle = sin(rotAngle);
    double cosAngle = cos(rotAngle);
    double one_cosAngle = 1.0 - cosAngle;

    setIdent();

    mat[0][0] = one_cosAngle * axis.x(0) * axis.x(0) + cosAngle;
    mat[0][1] = one_cosAngle * axis.x(0) * axis.x(1) - sinAngle * axis.x(2);
    mat[0][2] = one_cosAngle * axis.x(0) * axis.x(2) + sinAngle * axis.x(1);

    mat[1][0] = one_cosAngle * axis.x(1) * axis.x(0) + sinAngle * axis.x(2);
    mat[1][1] = one_cosAngle * axis.x(1) * axis.x(1) + cosAngle;
    mat[1][2] = one_cosAngle * axis.x(1) * axis.x(2) - sinAngle * axis.x(0);

    mat[2][0] = one_cosAngle * axis.x(2) * axis.x(0) - sinAngle * axis.x(1);
    mat[2][1] = one_cosAngle * axis.x(2) * axis.x(1) + sinAngle * axis.x(0);
    mat[2][2] = one_cosAngle * axis.x(2) * axis.x(2) + cosAngle;

    matIsUnit = (mat[0][0]==1.0) && (mat[0][1]==0.0) && (mat[0][2]==0.0) &&
                (mat[1][0]==0.0) && (mat[1][1]==1.0) && (mat[1][2]==0.0) &&
                (mat[2][0]==0.0) && (mat[2][1]==0.0) && (mat[2][2]==1.0);
}

RotateMat3D::RotateMat3D(enum_axisID mainAxis, double rotAngle) {
    setIdent();

    switch (mainAxis) {
        case axis_X:
            mat[1][1] =  mat[2][2] = cos(rotAngle) ;
            mat[1][2] = - (mat[2][1] = sin(rotAngle)) ;
            matIsUnit = (mat[1][1] == 1.0) && (mat[1][2] == 0.0);
            break;
        case axis_Y:
            mat[0][0] = mat[2][2] = cos(rotAngle) ;
            mat[2][0] = - (mat[0][2] = sin(rotAngle)) ;
            matIsUnit = (mat[0][0] == 1.0) && (mat[2][0] == 0.0);
            break;
        case axis_Z:
            mat[0][0] = mat[1][1] = cos(rotAngle) ;
            mat[0][1] = - (mat[1][0] = sin(rotAngle)) ;
            matIsUnit = (mat[0][0] == 1.0) && (mat[0][1] == 0.0);
            break;
    }
}

//----------------------------------------------------------------------------
//     RotateMat3Dd:    constructor, destructor
//----------------------------------------------------------------------------
RotateMat3Dd::RotateMat3Dd() {
    setIdent();
}

RotateMat3Dd::RotateMat3Dd(const vec3 &rotAxis, double rotAngle) {
    vec3 axis = rotAxis.getNormalized();

    double sinAngle = sin(rotAngle);
    double cosAngle = cos(rotAngle);
    double one_cosAngle = 1.0 - cosAngle;

    setIdent();

    mat[0][0] = one_cosAngle * axis.x(0) * axis.x(0) + cosAngle;
    mat[0][1] = one_cosAngle * axis.x(0) * axis.x(1) - sinAngle * axis.x(2);
    mat[0][2] = one_cosAngle * axis.x(0) * axis.x(2) + sinAngle * axis.x(1);

    mat[1][0] = one_cosAngle * axis.x(1) * axis.x(0) + sinAngle * axis.x(2);
    mat[1][1] = one_cosAngle * axis.x(1) * axis.x(1) + cosAngle;
    mat[1][2] = one_cosAngle * axis.x(1) * axis.x(2) - sinAngle * axis.x(0);

    mat[2][0] = one_cosAngle * axis.x(2) * axis.x(0) - sinAngle * axis.x(1);
    mat[2][1] = one_cosAngle * axis.x(2) * axis.x(1) + sinAngle * axis.x(0);
    mat[2][2] = one_cosAngle * axis.x(2) * axis.x(2) + cosAngle;

    matIsUnit = (mat[0][0]==1.0) && (mat[0][1]==0.0) && (mat[0][2]==0.0) &&
                (mat[1][0]==0.0) && (mat[1][1]==1.0) && (mat[1][2]==0.0) &&
                (mat[2][0]==0.0) && (mat[2][1]==0.0) && (mat[2][2]==1.0);
}

RotateMat3Dd::RotateMat3Dd(enum_axisID mainAxis, double rotAngle) {
    setIdent();
    switch (mainAxis) {
        case axis_X: {
            mat[1][1] =  mat[2][2] = cos(rotAngle);
            mat[1][2] = - (mat[2][1] = sin(rotAngle));
            matIsUnit = (mat[1][1] == 1.0) && (mat[1][2] == 0.0);
            break;
        }
        case axis_Y: {
            mat[0][0] = mat[2][2] = cos(rotAngle);
            mat[2][0] = - (mat[0][2] = sin(rotAngle));
            matIsUnit = (mat[0][0] == 1.0) && (mat[2][0] == 0.0);
            break;
        }
        case axis_Z: {
            mat[0][0] = mat[1][1] = cos(rotAngle);
            mat[0][1] = - (mat[1][0] = sin(rotAngle));
            matIsUnit = (mat[0][0] == 1.0) && (mat[0][1] == 0.0);
            break;
        }
    }
}


//----------------------------------------------------------------------------
//     RotateMat3Df:    constructor, destructor
//----------------------------------------------------------------------------
RotateMat3Df::RotateMat3Df() {
    setIdent();
}

RotateMat3Df::RotateMat3Df(const vec3f &rotAxis, float rotAngle) {
    vec3f axis = rotAxis.getNormalized();

    float sinAngle = sinf(rotAngle);
    float cosAngle = cosf(rotAngle);
    float one_cosAngle = 1.0 - cosAngle;

    setIdent();

    mat[0][0] = one_cosAngle * axis.x(0) * axis.x(0) + cosAngle;
    mat[0][1] = one_cosAngle * axis.x(0) * axis.x(1) - sinAngle * axis.x(2);
    mat[0][2] = one_cosAngle * axis.x(0) * axis.x(2) + sinAngle * axis.x(1);

    mat[1][0] = one_cosAngle * axis.x(1) * axis.x(0) + sinAngle * axis.x(2);
    mat[1][1] = one_cosAngle * axis.x(1) * axis.x(1) + cosAngle;
    mat[1][2] = one_cosAngle * axis.x(1) * axis.x(2) - sinAngle * axis.x(0);

    mat[2][0] = one_cosAngle * axis.x(2) * axis.x(0) - sinAngle * axis.x(1);
    mat[2][1] = one_cosAngle * axis.x(2) * axis.x(1) + sinAngle * axis.x(0);
    mat[2][2] = one_cosAngle * axis.x(2) * axis.x(2) + cosAngle;

    matIsUnit = (mat[0][0]==1.0) && (mat[0][1]==0.0) && (mat[0][2]==0.0) &&
                (mat[1][0]==0.0) && (mat[1][1]==1.0) && (mat[1][2]==0.0) &&
                (mat[2][0]==0.0) && (mat[2][1]==0.0) && (mat[2][2]==1.0);
}

RotateMat3Df::RotateMat3Df(enum_axisID mainAxis, float rotAngle) {
    setIdent();

    switch (mainAxis) {
        case axis_X:
            mat[1][1] =  mat[2][2] = cos(rotAngle) ;
            mat[1][2] = - (mat[2][1] = sin(rotAngle)) ;
            matIsUnit = (mat[1][1] == 1.0) && (mat[1][2] == 0.0);
            break;
        case axis_Y:
            mat[0][0] = mat[2][2] = cos(rotAngle) ;
            mat[2][0] = - (mat[0][2] = sin(rotAngle)) ;
            matIsUnit = (mat[0][0] == 1.0) && (mat[2][0] == 0.0);
            break;
        case axis_Z:
            mat[0][0] = mat[1][1] = cos(rotAngle) ;
            mat[0][1] = - (mat[1][0] = sin(rotAngle)) ;
            matIsUnit = (mat[0][0] == 1.0) && (mat[0][1] == 0.0);
            break;
    }
}

//----------------------------------------------------------------------------
//       ScaleMat3D:   constructor, destructor
//----------------------------------------------------------------------------
ScaleMat3D::ScaleMat3D(double sx, double sy, double sz) {
    setIdent();

    mat[0][0] = sx;
    mat[1][1] = sy;
    mat[2][2] = sz;

    matIsUnit = (sx == 1.0) && (sy == 1.0) && (sz == 1.0);
}

ScaleMat3D::ScaleMat3D(const vec3 &scale) {
    setIdent();

    mat[0][0] = scale[0];
    mat[1][1] = scale[1];
    mat[2][2] = scale[2];

    matIsUnit = (mat[0][0] == 1.0) && (mat[1][1] == 1.0) && (mat[2][2] == 1.0);
}



//----------------------------------------------------------------------------
//     TranslateMat2D:    constructor, destructor
//----------------------------------------------------------------------------
TranslateMat2D::TranslateMat2D(double tx, double ty) {
    setIdent();

    mat[0][2] = tx;
    mat[1][2] = ty;

    matIsUnit = (tx == 0.0) && (ty == 0.0);
}

//-------------------------------------------------------------------------
//    RotateMat2D::constructor, destructor
//-------------------------------------------------------------------------

RotateMat2D::RotateMat2D(double rotAngle) {
    setIdent();

    mat[0][0] = mat[1][1] = cos(rotAngle) ;
    mat[0][1] = - (mat[1][0] = sin(rotAngle)) ;         // ?????
}

RotateMat2D::RotateMat2D(double rotCenterX, double rotCenterY, double rotAngle) {
    Matrix<double,2,3>  trt = TranslateMat2D(rotCenterX, rotCenterY) *
                              RotateMat2D(rotAngle) *
                              TranslateMat2D(-rotCenterX, -rotCenterY);

    int  i ;
    int  j ;
    for (i = 0; i < 2; i++)
        for (j = 0; j < 3; j++) {
            mat[i][j] = trt.getElem(i,j);
        }

}


//-------------------------------------------------------------------------
//   ScaleMat2D::constructor, destructor
//-------------------------------------------------------------------------

ScaleMat2D::ScaleMat2D(double sx, double sy) {
    setIdent();

    mat[0][0] = sx;
    mat[1][1] = sy;
}



//----------------------------------------------------------------------------
//
//----------------------------------------------------------------------------
LorentzTransf::LorentzTransf(double beta, VnD<double,3> n) {
    assert(beta>-1.0 && beta<1.0);

    double gamma = 1.0/sqrt(1.0-beta*beta);
    VnD<double,3> ndir = n;
    ndir = ndir.getNormalized();

    setIdent();
    mat[0][0] = gamma;
    for (int row=1; row<4; row++) {
        for (int col=1; col<4; col++) {
            mat[row][col] = (gamma-1.0)*ndir[row-1]*ndir[col-1] + M4D_DELTA(row,col);
        }
        mat[0][row] = beta*gamma*ndir[row-1];
        mat[row][0] = beta*gamma*ndir[row-1];
    }
}

} // end namespace m4d

