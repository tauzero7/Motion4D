// --------------------------------------------------------------------------------
/*
    Mat.h

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

/*!  \class  m4d::Matrix
     \brief  Template class for n x m  matrices.

             n: number of rows, m: number of columns.
*/
// --------------------------------------------------------------------------------

#ifndef M4D_MAT_H
#define M4D_MAT_H

#include <cstdio>
#include <iostream>
#include <cstring>
#include <typeinfo>
#include <cassert>

#include "VnD.h"

#ifdef _WIN32
#ifndef __GNUC__
#pragma warning (disable: 4244 )
#endif
#endif

#ifdef _WIN32
#ifdef MATH_EXPORTS
#define MATH_API __declspec(dllexport)
#else /* METRIC_EXPORTS */
#define MATH_API __declspec(dllimport)
#endif /* METRIC_EXPORTS */
#else /* _WIN32 */
#define MATH_API
#endif /* _WIN32 */

namespace m4d {

//---------------------------------------------------------------------------
//    class-template  mType
//---------------------------------------------------------------------------
template <class mType, int n, int m> class  Matrix {
protected:
    mType**      mat;
    std::string  classType;
    int          nr,nc;  // number of rows, number of columns
    bool         matIsUnit;

public:
    Matrix();
    Matrix(double val);
    Matrix(const Matrix &M);
    Matrix(mType* field);
    ~Matrix();

    void         setAll(mType val);
    void         setElem(int row, int col, mType val);
    mType        getElem(int row, int col) const;

    void         setRow(int row, const VnD<mType,m> &vec);        // set row-vector
    VnD<mType,m> getRow(int row) const;

    void         setCol(int col, const VnD<mType,n> &vec);        // set col-vector
    VnD<mType,n> getCol(int col) const;

    int          numRows() const {
        return nr;
    }
    int          numCols() const {
        return nc;
    }

    void         setNull();   // make 0-matrix;
    void         setIdent();  // make identity-matrix

    bool         isIdentMat() const;

    void         transpose();
    void         invert();

    void         getDoubleArray(double val[]);
    void         getFloatArray(float val[]);

    void         copyComponents(const Matrix<mType,n,m> &M);

    VnD<mType,m>       operator[](int z);
    VnD<mType,n>       operator*(const VnD<mType,m> &v);

    void               operator=(const Matrix<mType,n,m> &M);
    Matrix<mType,n,m>  operator+(const Matrix<mType,n,m> &M) const;
    Matrix<mType,n,m>  operator-(const Matrix<mType,n,m> &M) const;
    Matrix<mType,n,m>  operator|(const Matrix<mType,n,m> &M) const;
    Matrix<mType,n,m>  operator*(const mType a) const;

    Matrix<mType,n,m>  operator^(const int l);

    void               operator*=(const Matrix<mType,n,m> &B);

    void   print(std::ostream& ostr = std::cerr) const;
    void   printS(FILE* fptr = stderr, const std::string format = "%12.8f ") const;

    // --------------------- friend functions ---------------------

    // double * Matrix
    friend Matrix<mType,n,m> operator* (const double a, Matrix<mType,n,m> &M) {
        Matrix<mType,n,m> prod;
        for (int i=0; i<n; i++) {
            for (int j=0; j<m; j++) {
                prod.setElem(i,j,a*M.getElem(i,j));
                prod.matIsUnit = prod.matIsUnit && (prod.getElem(i,j) == double(i==j));
            }
        }
        return prod;
    }

    //  Matrix * VnD (specific matrix vektor multiplication)
    /*
    friend VnD<mType,n> operator* (const Matrix<mType,n,m> &M, const VnD<mType,n> &vec)  {
      VnD<mType,n> q;
      for(int i=0; i<n; i++) {
      q[i] = (mType)0;
      for(int j=0; j<n; j++) {
          q[i] += M.getElem(i,j)*vec.x(j);
      }
      }
      return q;
    }
    */
    friend VnD<mType,n> operator* (const Matrix<mType,n,m> &M, const VnD<mType,n> &vec)  {
        VnD<mType,n> q;
        for (int i=0; i<n; i++) {
            q[i] = (mType)0;
            for (int j=0; j<n; j++) {
                q[i] += M.getElem(i,j)*vec.x(j);
            }
            q[i]+=M.getElem(i,m-1);
        }
        return q;
    }

    // transposeMult  (specific multiplication)
    friend VnD<mType,n> transposeMult(const Matrix<mType,n,m> &M, const VnD<mType,n> &vec) {
        VnD<mType,n> q;
        for (int i=0; i<n; i++) {
            q[i] = (mType)0;
            for (int j=0; j<n; j++) {
                q[i] += M.getElem(j,i)*vec.x(j);
            }
        }
        return q;
    }

    // Matrix * Matrix  ((specific) matrix matrix multiplication)
    friend Matrix<mType,n,m>  operator*(const Matrix<mType,n,m> &m1, const Matrix<mType,n,m> &m2) {
        assert(n+1==m || n==m);
        Matrix<mType,n,m> prod;
        mType sum;
        for (int i=0; i<n; i++) {
            for (int j=0; j<n; j++) {
                sum = 0;
                for (int k=0; k<n; k++) {
                    sum += m1.getElem(i,k)*m2.getElem(k,j);
                }
                prod.setElem(i,j,sum);
            }
            // hier geht die spezielle M-M-Multi weiter
            if ((n+1)==m) {
                sum = 0;
                for (int k=0; k<n; k++) {
                    sum += m1.getElem(i,k)*m2.getElem(k,m-1);
                }
                prod.setElem(i,m-1, sum + m1.getElem(i,m-1));
            }
        }
        return prod;
    }
};



/**
 *  Standard constructor results in null matrix
 */
template <class mType, int n, int m> Matrix<mType,n,m>::Matrix() {
    nr = n;
    nc = m;

    // initialize matrix
    mat = new mType*[n];
    for (int i=0; i<n; i++) {
        mat[i] = new mType[m];
    }
    setNull();  // write null-matrix
}

/**
 *  Standard constructor with fixed entry
 *
 * @param val value
 */
template <class mType, int n, int m> Matrix<mType,n,m>::Matrix(double val) {
    nr = n;
    nc = m;

    // initialize matrix
    mat = new mType*[n];
    for (int i=0; i<n; i++) {
        mat[i] = new mType[m];
    }
    setAll(val);
}

/**
 *  Constructor given a 1d array
 *
 * @param field  1D array
 */
template <class mType, int n, int m> Matrix<mType,n,m>::Matrix(mType* field) {
    nr = n;
    nc = m;

    // initialize matrix
    mat = new mType*[n];
    for (int i=0; i<n; i++) {
        mat[i] = new mType[m];
    }

    for (int i=0; i<n; i++) {
        for (int j=0; j<m; j++) {
            mat[i][j] = field[nr*i+j];
        }
    }
}

/**
 *  Copy-constructor
 */
template <class mType, int n, int m> Matrix<mType,n,m>::Matrix(const Matrix<mType,n,m> &M) {
    mat=new mType*[n];
    for (int i=0; i<n; i++) {
        mat[i]=new mType[m];
    }

    int size=sizeof(mType)*m;
    for (int i=0; i<n; i++) {
        memcpy(mat[i],M.mat[i],size);
    }
}

/**
 *  Destructor
 */
template <class mType, int n, int m> Matrix<mType,n,m>::~Matrix() {
    for (int i=0; i<n; i++) {
        delete [] mat[i];
    }
    delete [] mat;

    mat = NULL;
}


/**
 *  Set all matrix elements to val.
 *
 * @param val  value
 */
template <class mType, int n, int m> void Matrix<mType,n,m>::setAll(mType val) {
    for (int i=0; i<n; i++) {
        for (int j=0; j<m; j++) {
            mat[i][j] = val;
        }
    }
}

/**
 *  Set matrix element.
 *
 * @param row  row index
 * @param col  column index
 * @param val  value
 */
template <class mType, int n, int m> void Matrix<mType,n,m>::setElem(int row, int col, mType val) {
    mat[row][col] = val;
}

/**
 *  Get matrix element.
 *
 * @param row  row index
 * @param col  column index
 * @return value
 */
template <class mType, int n, int m> mType Matrix<mType,n,m>::getElem(int row, int col) const {
    return mat[row][col];
}

//---------------------------------------------------------------------------
//      set/get Row
//---------------------------------------------------------------------------
template <class mType, int n, int m> void  Matrix<mType,n,m>::setRow(int row, const VnD<mType,m> &vec) {
    for (int j=0; j<m; j++) {
        mat[row][j] = vec.x(j);
    }
}

template <class mType, int n, int m> VnD<mType,m> Matrix<mType,n,m>::getRow(int row) const {
    VnD<mType,m> vec;
    for (int j=0; j<m; j++) {
        vec[j] = mat[row][j];
    }
    return vec;
}

//---------------------------------------------------------------------------
//      set/get Col
//---------------------------------------------------------------------------
template <class mType, int n, int m> void  Matrix<mType,n,m>::setCol(int col, const VnD<mType,n> &vec) {
    for (int j=0; j<n; j++) {
        mat[j][col] = vec.x(j);
    }
}

template <class mType, int n, int m> VnD<mType,n> Matrix<mType,n,m>::getCol(int col) const {
    VnD<mType,n> vec;
    for (int j=0; j<n; j++) {
        vec[j] = mat[j][col];
    }
    return vec;
}


//---------------------------------------------------------------------------
//      clear()
//---------------------------------------------------------------------------
template <class mType, int n, int m> void Matrix<mType,n,m>::setNull() {
    for (int i=0; i<n; i++) {
        for (int j=0; j<m; j++) {
            mat[i][j] = (mType)0;
        }
    }
}

//---------------------------------------------------------------------------
//      setIdent()
//---------------------------------------------------------------------------
template <class mType, int n, int m> void Matrix<mType,n,m>::setIdent() {
    for (int i=0; i<n; i++) {
        for (int j=0; j<m; j++) {
            if (i==j) {
                mat[i][j] = (mType)1;
            } else {
                mat[i][j] = (mType)0;
            }
        }
    }
}

//---------------------------------------------------------------------------
//      isIdentMat()
//---------------------------------------------------------------------------
template <class mType, int n, int m> bool Matrix<mType,n,m>::isIdentMat() const {
    for (int i=0; i<n; i++) {
        for (int j=0; j<n; j++) {
            if (mat[i][j] != (mType)(i==j)) {
                return false;
            }
        }
    }
    if (n+1==m) {
        for (int j=0; j<n; j++) {
            if (mat[j][m-1] != (mType)(0)) {
                return false;
            }
        }
    }
    return true;
}

//---------------------------------------------------------------------------
//      transpose()
//---------------------------------------------------------------------------
template <class mType, int n, int m> void Matrix<mType,n,m>::transpose() {
    mType tmp;
    for (int i=0; i<n; i++) {
        for (int j=i; j<m; j++) {
            tmp = mat[i][j];
            mat[i][j] = mat[j][i];
            mat[j][i] = tmp;
        }
    }
}

//---------------------------------------------------------------------------
//      invert()
//---------------------------------------------------------------------------
template <class mType, int n, int m> void Matrix<mType,n,m>::invert() {
    if ((n==m) || (n+1==m)) {
        int    i, j, k, r, s, pr, ps, rowj, colj;
        int* row = new int[m];
        int* col = new int[m];
        double Pivot ;

        // --- initialize ---
        double** A;
        A = new double*[m];
        for (int i=0; i<m; i++) {
            A[i] = new double[m];
        };

        for (i = 0 ; i < n ; i++) {
            for (j = 0 ; j < m ; j++) {
                A[i][j] = mat[i][j];
            }
        }

        if (n+1==m) {
            for (j=0; j<m; j++) {
                A[m-1][j] = 0.0;
            }
            A[m-1][m-1] = 1.0;
        }

        for (i = 0 ; i < m ; i++) {
            row[i] = col[i] = i ;
        }

        //
        for (j = 0 ; j < m ; j++) {
            pr = j ;
            ps = j ;

            for (r = j ; r < m ; r++) {
                for (s = j ; s < m ; s++) {
                    if (fabs(A[row[r]][col[s]]) > fabs(A[row[pr]][col[ps]])) {
                        pr = r ;
                        ps = s ;
                    }
                }
            }

            //
            if (pr > j) {
                int temp = row[j];
                row[j] = row[pr];
                row[pr] = temp;
            }

            if (ps > j) {
                int temp = col[j];
                col[j] = col[ps];
                col[ps] = temp;
            }

            rowj  = row[j] ;
            colj  = col[j] ;
            Pivot = A[rowj][colj] ;

            for (i = 0 ; i < m ; i++) {
                if (i != rowj) {
                    for (k = 0 ; k < m ; k++) {
                        if (k != colj) {
                            A[i][k] -= (A[i][colj] * A[rowj][k] / Pivot) ;
                        }
                    }
                }
            }


            for (k = 0 ; k < m ; k++) {
                if (k != colj) {
                    A[rowj][k] /= (-Pivot) ;
                }
            }


            for (i = 0 ; i < m ; i++) {
                if (i != rowj) {
                    A[i][colj] /= Pivot ;
                }
            }

            A[rowj][colj] = 1.0 / Pivot ;

        } /*for j */

        for (i = 0 ; i < m ; i++) {
            for (k = 0 ; k < m ; k++) {
                if (col[i]<n) {
                    mat[col[i]][row[k]] = A[row[i]][col[k]] ;
                }
            }
        }

        for (i=0; i<m; i++) {
            delete [] A[i];
        }
        delete [] A;
        delete [] row;
        delete [] col;
    }
}


/**
 *  Get matrix as double array
 *
 * @param val  pointer to double array
 */
template <class mType, int n, int m> void Matrix<mType,n,m>::getDoubleArray(double val[]) {
    for (int i=0; i<n; i++) {
        for (int j=0; j<m; j++) {
            val[i*m+j] = mat[i][j];
        }
    }
}

/**
 *  Get matrix as float array
 *
 * @param val  pointer to float array
 */
template <class mType, int n, int m> void Matrix<mType,n,m>::getFloatArray(float val[]) {
    for (int i=0; i<n; i++) {
        for (int j=0; j<m; j++) {
            val[i*m+j] = float(mat[i][j]);
        }
    }
}

/**
 *  Copy matrix components
 */
template <class mType, int n, int m> void Matrix<mType,n,m>::copyComponents(const Matrix<mType,n,m> &M) {
    for (int i=0; i<n; i++) {
        for (int j=0; j<m; j++) {
            mat[i][j] = M.getElem(i,j);
        }
    }
}


//---------------------------------------------------------------------------
//      operator[](int i)  //access row vector
//---------------------------------------------------------------------------
template <class mType, int n, int m> VnD<mType,m> Matrix<mType,n,m>::operator[](int z) {
    assert(z<nr);
    VnD<mType,m> vec;
    for (int i=0; i<nc; i++) {
        vec.setX(i,mat[z][i]);
    }
    return vec;
}

//---------------------------------------------------------------------------
//      operator=
//---------------------------------------------------------------------------
template <class mType, int n, int m> void Matrix<mType,n,m>::operator=(const Matrix<mType,n,m> &M) {
    /*
    if (this != &M) {
        for (int i=0; i<nr; i++) {
            delete [] mat[i];
        }
        delete [] mat;
    }

    mat = new mType*[nr];
    for (int i=0; i<nr; i++) {
        mat[i] = new mType[nc];
    }
    */

    int size=sizeof(mType)*nc;
    for (int i=0; i<nr; i++) {
        memcpy(mat[i],M.mat[i],size);
    }
}


/**
 *  Operator+
 */
template <class mType, int n, int m> Matrix<mType,n,m> Matrix<mType,n,m>::operator+(const Matrix<mType,n,m> &M) const {
    Matrix<mType,n,m> Q;
    for (int i=0; i<n; i++) {
        for (int j=0; j<m; j++) {
            Q.mat[i][j] = mat[i][j] + M.mat[i][j];
        }
    }
    return Q;
}

/**
 *  Operator-
 */
template <class mType, int n, int m> Matrix<mType,n,m> Matrix<mType,n,m>::operator-(const Matrix<mType,n,m> &M) const {
    Matrix<mType,n,m> Q;
    for (int i=0; i<n; i++) {
        for (int j=0; j<m; j++) {
            Q.mat[i][j] = mat[i][j] - M.mat[i][j];
        }
    }
    return Q;
}

/**
 *  Operator|
 */
template <class mType, int n, int m> Matrix<mType,n,m> Matrix<mType,n,m>::operator|(const Matrix<mType,n,m> &M) const {
    assert(n==m);
    Matrix<mType,n,m> Q;

    for (int i=0; i<n; i++) {
        for (int j=0; j<m; j++) {
            Q.mat[i][j] = (mType)0;
            for (int k=0; k<m; k++) {
                Q.mat[i][j] += mat[i][k]*M.mat[k][j];
            }
        }
    }

    return Q;
}

/**
 *  Scalar multiplication
 *
 * @param a  scalar value
 */
template <class mType, int n, int m> Matrix<mType,n,m> Matrix<mType,n,m>::operator*(const mType a) const {
    Matrix<mType,n,m> Q;

    for (int i=0; i<n; i++) {
        for (int j=0; j<m; j++) {
            for (int k=0; k<m; k++) {
                Q.mat[i][j] = a*mat[i][k];
            }
        }
    }
    return Q;
}


/**
 *  Matrix-Vector multiplication
 *
 * @param v  vector
 */
template <class mType, int n, int m> VnD<mType,n> Matrix<mType,n,m>::operator*(const VnD<mType,m> &v) {
    VnD<mType,n> q;
    for (int i=0; i<n; i++) {
        q[i] = (mType)0;
        for (int j=0; j<m; j++) {
            q[i] += mat[i][j]*v.x(j);
        }
    }
    return q;
}

//---------------------------------------------------------------------------
//      operator = matrix
//---------------------------------------------------------------------------
/*
template <class mType, int n, int m> const Matrix<mType,n,m>& Matrix<mType,n,m>::operator=( const Matrix<mType,n,m> &B )
{
   if(this != &B)
   {
      int  i ;
      int  j ;
      for(  i = 0; i < 4; i++ )
         for(  j = 0; j < 4; j++ )
            mat[i][j] = B.mat[i][j];

      matIsUnit = B.matIsUnit;
   }

   return *this;
}
*/

//---------------------------------------------------------------------------
//      operator*= matrix
//---------------------------------------------------------------------------
template <class mType, int n, int m> void Matrix<mType,n,m>::operator*=(const Matrix<mType,n,m> &B) {
    Matrix<mType,n,m> prod(0.0);
    prod.matIsUnit = true;

    register int  i ;
    register int  j ;
    register int  k ;
    //int limit;

    assert((n==m) || (n+1 == m));

    if (n==m) {
        for (i=0; i<n; i++) {
            for (j=0; j<n; j++) {
                for (k=0; k<n; k++) {
                    prod.mat[i][j] += mat[i][k] * B.mat[k][j];
                }

                prod.matIsUnit = prod.matIsUnit && (prod.mat[i][j] == mType(i==j));
            }
        }
    } else if (n+1 == m) {
        for (i=0; i<n; i++) {
            for (j=0; j<m; j++) {
                for (k=0; k<n; k++) {
                    prod.mat[i][j] += mat[i][k] * B.mat[k][j];
                }

                prod.matIsUnit = prod.matIsUnit && (prod.mat[i][j] == mType(i==j));
            }
            prod.mat[i][m-1] += mat[i][m-1];
        }
    }

    *this = prod;
}


//---------------------------------------------------------------------------
//      operator^ =matrix
//---------------------------------------------------------------------------
template <class mType, int n, int m> Matrix<mType,n,m> Matrix<mType,n,m>::operator^(const int l) {
    assert(n==m);
    Matrix<mType,n,m> pow;
    pow.setIdent();

    if (l==0) {
        return pow;
    }

    if (l>0) {
        for (int i=0; i<l; i++) {
            pow = pow*mat;
        }
        pow.matIsUnit = isIdentMat();
        return pow;
    }

    return mat;
}


/**
 *  Print matrix
 *
 * @param ostr   output stream
 */
template <class mType, int n, int m> void Matrix<mType,n,m>::print(std::ostream& ostr) const {
    for (int i=0; i<n; i++) {
        ostr << "( ";
        for (int j=0; j<m; j++) {
            ostr << mat[i][j] << "\t";
        }
        ostr << ")" << std::endl;
    }
}

/**
 *  Print matrix
 *
 * @param fptr   file pointer
 * @param format  output format
 */
template <class mType, int n, int m> void  Matrix<mType,n,m>::printS(FILE* fptr, const std::string format) const {
    for (int i=0; i<n; i++) {
        fprintf(fptr,"(");
        for (int j=0; j<m; j++) {
            fprintf(fptr,format.c_str(),mat[i][j]);
        }
        fprintf(fptr,")\n");
    }
}

} // end namespace m4d

#endif

