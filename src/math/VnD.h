// --------------------------------------------------------------------------------
/*
    VnD.h

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

/*!  \class  m4d::VnD
     \brief  Template class for n-dimensional vectors.


*/
// --------------------------------------------------------------------------------

#ifndef M4D_VnD_H
#define M4D_VnD_H

#ifndef _WIN32
#ifndef DBL_MAX
#define DBL_MAX 1.844674407370955616e19
#endif
#endif

#include <cstdio>
#include <iostream>
#include <string>
#include <typeinfo>
#include <cassert>
#include <cmath>

namespace m4d {

static double epsilon=1.0e-8;

//---------------------------------------------------------------------------
//    class-template  vType
//---------------------------------------------------------------------------
template <class vType, int size> class VnD {
    vType v[size];
    std::string classType;


public:
    //! Standard constructor initialized to zero.
    VnD() {
        for (int i=0; i<size; i++) {
            v[i] = (vType)0;
        }
    }

    VnD(const VnD<vType,size> &vec) {
        for (int i=0; i<size; i++) {
            v[i] = vec[i];
        }
    }

    //! Constructor for vec2.
    VnD(vType x1, vType x2) {
        assert(size==2);
        v[0]=x1;
        v[1]=x2;
    }

    //! Constructor for vec3.
    VnD(vType x1, vType x2, vType x3) {
        assert(size==3);
        v[0]=x1;
        v[1]=x2;
        v[2]=x3;
    }

    //! Constructor for vec4.
    VnD(vType x1, vType x2, vType x3, vType x4) {
        assert(size==4);
        v[0]=x1;
        v[1]=x2;
        v[2]=x3;
        v[3]=x4;
    }

    //! Constructor for vec5.
    VnD(vType x1, vType x2, vType x3, vType x4, vType x5) {
        assert(size==5);
        v[0]=x1;
        v[1]=x2;
        v[2]=x3;
        v[3]=x4;
        v[4]=x5;
    }

    VnD(const vType* x, unsigned int s = size) {
        unsigned int i=0;
        while (i<s && i<size) {
            v[i] = x[i];
            i++;
        }
    }
    //! Standard destructor.
    ~VnD() {}

    //! Return pointer to const array.
    const vType*  data() const {
        return v;
    }
    //! Return pointer to array.
    vType*        data() {
        return v;
    }

    //! Set vector by array.
    void set(vType y[size]);
    //! Set x-th component.
    void setX(int nr, vType val) {
        v[nr] = val;
    }


    // define operators
    const vType     &operator[](int i) const;
    vType           &operator[](int i);
    void            operator= (const VnD<vType,size> &vec);
    VnD<vType,size> operator+ (const VnD<vType,size> &vec) const;
    void            operator+=(const VnD<vType,size> &vec);
    VnD<vType,size> operator- (const VnD<vType,size> &vec) const;
    VnD<vType,size> operator- () const;

    //! Multiplication with scalar.
    VnD<vType,size> operator* (const vType a) const;
    //! Multiplication with scalar.
    void            operator*=(const vType a);
    //! Multiplication with scale vector.
    void            operator*=(const VnD<vType,size> &vec);

    //! Division by scalar.
    VnD<vType,size> operator/ (const vType a) const;
    //! Division by scalar.
    void            operator/=(const vType a);

    //! Dot product.
    vType           operator| (const VnD<vType,size> &vec) const;
    //! Cross product.
    VnD<vType,size> operator^ (const VnD<vType,size> &vec) const;

    //! Logical operators.
    int             operator==(const VnD<vType,size> &vec) const;
    //! Logical operators.
    int             operator!=(const VnD<vType,size> &vec) const;

    // define return values
    vType            x(int nr) const {
        if (nr>=0 && nr < size) {
            return v[nr];
        } else {
            return (vType)0;
        }
    }
    /*
    vType&           x(int nr)       {
        if (nr>=0 && nr < size) {
            return v[nr];
        } else {
            return (vType)0;
        }
    }
    */

    //! Get norm of vector.
    vType            getNorm() const;
    //! Normalize vector and return norm.
    vType            normalize();
    //! Normalize vector and return vector.
    VnD<vType,size>  getNormalized() const;

    void             vabs();
    VnD<vType,size>  getVabs() const;

    int              getSize() { return size; }

    int              mostDominantCoord() const;
    int              leastDominantCoord() const;

    bool             isZero(void) const;

    void             type(void) const;
    //! Print to stdout.
    void             print(int nks = 4) const;
    void             printS(FILE* fptr = stderr, const std::string format = std::string("%12.8f ")) const;
    void             print(std::ostream & os, std::string text = std::string(), int nks = 4) const;
    void             printF(std::ostream &os=std::cerr) const;

    VnD<vType,3>     getAsV3D() const;
    VnD<vType,4>     getAsV4D() const;
    VnD<vType,4>     get3As4() const;


    // --------- friend functions ---------
    friend VnD<vType,size> operator*(double a, const VnD<vType,size> &vec)  {
        VnD<vType,size> q;
        for (int i=0; i<size; i++) {
            q[i] = a*vec.x(i);
        }
        return q;
    }

    friend std::ostream& operator<<(std::ostream& s, const VnD<vType,size> &vec) {
        s << "( ";
        for (int i=0; i<size; i++) {
            s << vec.x(i) << " ";
        }
        s << ")";
        return s;
    }

//  friend istream& operator>>(istream& s, V3D& v);

};


//---------------------------------------------------------------------------
//      set
//---------------------------------------------------------------------------
template <class vType, int size> void VnD<vType, size>::set(vType y[size]) {
    for (int i=0; i<size; i++) {
        v[i]=y[i];
    }
}

//---------------------------------------------------------------------------
//      operator[]
//---------------------------------------------------------------------------
template <class vType, int size> const vType &VnD<vType,size>::operator[](int i) const {
    return v[i];
}

template <class vType, int size> vType &VnD<vType,size>::operator[](int i) {
    return v[i];
}

//---------------------------------------------------------------------------
//      operator=
//---------------------------------------------------------------------------
template <class vType, int size> void VnD<vType, size>::operator=(const VnD<vType,size> &vec) {
    for (int i=0; i<size; i++) {
        v[i] = vec.x(i);
    }
}

//---------------------------------------------------------------------------
//      operator+
//---------------------------------------------------------------------------
template <class vType, int size> VnD<vType,size> VnD<vType, size>::operator+(const VnD<vType,size> &vec) const {
    VnD<vType,size> q;
    for (int i=0; i<size; i++) {
        q[i] = v[i] + vec.x(i);
    }
    return q;
}

//---------------------------------------------------------------------------
//      operator+=
//---------------------------------------------------------------------------
template <class vType, int size> void VnD<vType, size>::operator+=(const VnD<vType,size> &vec) {
    for (int i=0; i<size; i++) {
        v[i] += vec[i];
    }
}

//---------------------------------------------------------------------------
//      operator-
//---------------------------------------------------------------------------
template <class vType, int size> VnD<vType,size> VnD<vType, size>::operator-(const VnD<vType,size> &vec) const {
    VnD<vType,size> q;
    for (int i=0; i<size; i++) {
        q[i] = v[i] - vec.x(i);
    }
    return q;
}

template <class vType, int size>  VnD<vType, size> VnD<vType, size>::operator-() const {
    VnD<vType,size> q;
    for (int i=0; i<size; i++) {
        q[i] = -v[i];
    }
    return q;
}

//---------------------------------------------------------------------------
//      operator* vType
//---------------------------------------------------------------------------
template <class vType, int size> VnD<vType,size> VnD<vType,size>::operator*(const vType a) const {
    VnD<vType,size> q;
    for (int i=0; i<size; i++) {
        q[i] = v[i]*a;
    }
    return q;
}

template <class vType, int size> vType VnD<vType,size>::operator|(const VnD<vType,size> &vec) const {
    vType q=0;
    for (int i=0; i<size; i++) {
        q += v[i]*vec.v[i];
    }
    return q;
}

template <class vType, int size> void VnD<vType,size>::operator*=(const vType a) {
    for (int i=0; i<size; i++) {
        v[i] *= a;
    }
}

template <class vType, int size> void VnD<vType,size>::operator*=(const VnD<vType,size> &vec) {
    for (int i=0; i<size; i++) {
        v[i] *= vec[i];
    }
}

template <class vType, int size> VnD<vType,size> VnD<vType,size>::operator/(const vType a) const {
    VnD<vType,size> q;
    for (int i=0; i<size; i++) {
        q[i] = v[i]/a;
    }
    return q;
}

template <class vType, int size> void VnD<vType,size>::operator/=(const vType a) {
    for (int i=0; i<size; i++) {
        v[i] /= a;
    }
}

//---------------------------------------------------------------------------
//      operator^ vType
//---------------------------------------------------------------------------

template <class vType, int size> VnD<vType,size> VnD<vType,size>::operator^(const VnD<vType,size> &vec) const {
    if (size == 3) {
        VnD<vType,size> q;
        q[0] = v[1]*vec.v[2]-v[2]*vec.v[1];
        q[1] = v[2]*vec.v[0]-v[0]*vec.v[2];
        q[2] = v[0]*vec.v[1]-v[1]*vec.v[0];
        return q;
    } else {
        VnD<vType,size> q;
        return q;
    }
}

//---------------------------------------------------------------------------
//      operator==
//---------------------------------------------------------------------------
template <class vType, int size> int VnD<vType,size>::operator==(const VnD<vType,size> &vec) const {
    bool isOkay = true;
    for (int i=0; i<size; i++) {
        isOkay &= (fabs(v[i]-vec.x(i))<=epsilon);
    }
    return isOkay;
}
//---------------------------------------------------------------------------
//      operator!=
//---------------------------------------------------------------------------
template <class vType, int size> int VnD<vType,size>::operator!=(const VnD<vType,size> &vec) const {
    return !(*this == vec);
}

//---------------------------------------------------------------------------
//     getNorm()
//---------------------------------------------------------------------------
template <class vType, int size> vType VnD<vType,size>::getNorm() const {
    vType sum=(vType)0;
    for (int i=0; i<size; i++) {
        sum+= v[i]*v[i];
    }
    return (vType)sqrt((double)sum);
}

//---------------------------------------------------------------------------
//     normalize()
//---------------------------------------------------------------------------
template <class vType, int size> vType VnD<vType,size>::normalize() {
    vType n = getNorm();
    if (n <= epsilon) {
        //fprintf(stderr,"VnD::normalize() ... cannot normalize; length <= epsilon.\n");
    } else {
        vType a = (vType)1/n;
        for (int i=0; i<size; i++) {
            v[i] = a*v[i];
        }
    }
    return n;
}

//---------------------------------------------------------------------------
//     getNormalized()
//---------------------------------------------------------------------------
template <class vType, int size> VnD<vType,size> VnD<vType,size>::getNormalized() const {
    VnD<vType,size> q;

    vType n = getNorm();
    if (n<=epsilon) {
        fprintf(stderr,"VnD::getNormalized() ... cannot be normalized!\n");
    } else {
        for (int i=0; i<size; i++) {
            q[i] = v[i]/n;
        }
    }

    return q;
}

//---------------------------------------------------------------------------
//     vabs(), getVabs()
//---------------------------------------------------------------------------
template <class vType, int size> void VnD<vType,size>::vabs() {
    for (int i=0; i<size; i++) {
        v[i] = fabs((double)v[i]);
    }
}

template <class vType, int size> VnD<vType,size> VnD<vType,size>::getVabs() const {
    VnD<vType,size> q;
    for (int i=0; i<size; i++) {
        q[i] = fabs((double)v[i]);
    }
    return q;
}

//---------------------------------------------------------------------------
//     DominantCoord
//---------------------------------------------------------------------------
template <class vType, int size> int VnD<vType,size>::mostDominantCoord() const {
    VnD<vType,size> q;
    int dominant = -1;
    double max = -DBL_MAX;
    for (int i=0; i<size; i++) {
        q[i] = fabs((double)v[i]);
        if (q[i]>max) {
            max=q[i];
            dominant = i;
        }
    }
    return dominant;
}

template <class vType, int size> int VnD<vType,size>::leastDominantCoord() const {
    VnD<vType,size> q;
    int dominant = -1;
    double min = DBL_MAX;
    for (int i=0; i<size; i++) {
        q[i] = fabs((double)v[i]);
        if (q[i]<min) {
            min=q[i];
            dominant = i;
        }
    }
    return dominant;
}

//---------------------------------------------------------------------------
//      isZero()
//---------------------------------------------------------------------------
template <class vType, int size> bool VnD<vType,size>::isZero() const {
    double val = 0.0;
    for (int i=0; i<size; i++) {
        val += (double)(v[i]*v[i]);
    }
    if (val<1e-8) {
        return true;
    }

    return false;
}

//---------------------------------------------------------------------------
//      type()
//---------------------------------------------------------------------------
template <class vType, int size> void VnD<vType,size>::type() const {
    std::cout << "Type: VnD { " << typeid(vType).name() << " of size " << size << " }"<< std::endl;
}


//---------------------------------------------------------------------------
//      print()
//---------------------------------------------------------------------------
template <class vType, int size> void VnD<vType,size>::print(int nks) const {
    // nks = Zahl der Nachkommastellen
    std::cout.precision(nks);
    std::cout << "( ";
    for (int i=0; i<size; i++) {
        std::cout << v[i] << " ";
    }
    std::cout << ")\n";

    // std::cout << "Class if of type " <<  classType << std::endl;
}

template <class vType, int size> void VnD<vType, size>::printS(FILE* fptr, const std::string format) const {
    for (int i=0; i<size; i++) {
        fprintf(fptr,format.c_str(),v[i]);
    }
    fprintf(fptr,"\n");
}


template <class vType, int size> void VnD<vType,size>::print(std::ostream &os, std::string text, int nks) const {
    // nks = Zahl der Nachkommastellen
    os.precision(nks);
    os << text << "( ";
    for (int i=0; i<size; i++) {
        os << v[i] << " ";
    }
    os << ")\n";

    // std::cout << "Class if of type " <<  classType << std::endl;
}


template <class vType, int size> void VnD<vType, size>::printF(std::ostream &os) const {
    os.setf(std::ios::showpoint);
    for (int i=0; i<size; i++) {
        os.width(10);
        os << v[i] << " ";
    }
    os << "\n";
}

//---------------------------------------------------------------------------
//      conversions
//---------------------------------------------------------------------------
template <class vType, int size> VnD<vType,3> VnD<vType, size>::getAsV3D() const {
    VnD<vType,3> nv;
    int k=0;
    int count = (size < 4) ? (size-1) : (3);

    while (k < count) { // copy 'spacelike' components
        nv[k] = v[k+1];
        k++;
    }

    while (k<3) {
        nv[k++] = 0.0;
    }

    return nv;

}

template <class vType, int size> VnD<vType,4> VnD<vType, size>::getAsV4D() const {
    VnD<vType,4> nv;
    int k=0;
    int count = (size < 5) ? (size) : (4);

    while (k < count) {
        nv[k] = v[k];
        k++;
    }

    while (k<4) {
        nv[k++] = 0.0;
    }
    return nv;
}

template <class vType, int size> VnD<vType,4> VnD<vType, size>::get3As4() const {
    assert(size==3);

    VnD<vType,4> nv;

    nv[0] = (vType)0;
    for (int k=0; k<3; k++) {
        nv[k+1] = v[k];
    }

    return nv;
}

} // end namespace m4d

#endif

