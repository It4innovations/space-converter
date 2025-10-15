#ifndef _OPENGADGET3_MATRIX_HPP_
#define _OPENGADGET3_MATRIX_HPP_


#include <iostream>
#include "../System/assert.h"
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <math.h>
#include "Eos.hpp"
#include "FluidVector.hpp"
#include "Vector.hpp"
#include "../CodeBase/precision.h"


//-------------------------------------------------------------------------------------------------
namespace OpenGadget3
{

  //===============================================================================================
  /// \brief   General matrix class for storing data in column-major order (similar to C/C++ arrays)
  /// \author  S. Heigl & D. A. Hubber
  /// \date    13/02/2018
  //===============================================================================================
  template <int M, int N>
  class Matrix
  {
  protected:

    MyFloat data[M][N];                          ///< Matrix data


  public:

    // Constructors and destructor
    //---------------------------------------------------------------------------------------------
    Matrix()
    {
      for (int i=0; i<M; i++) {
        for (int j=0; j<N; j++) {
          data[i][j] = MyFloat(0.0);
        }
      }
    }

    Matrix(const Matrix<M,N> &other)
    {
      for (int i=0; i<M; i++) {
        for (int j=0; j<N; j++) {
          data[i][j] = other(i,j);
        }
      }
    }

    Matrix(const MyFloat value)
    {
      assert(std::isfinite(value));
      for (int i=0; i<M; i++) {
        for (int j=0; j<N; j++) {
          data[i][j] = value;
        }
      }
    }

    virtual ~Matrix() {};


    // Static functions for generating matrices
    //---------------------------------------------------------------------------------------------
    static inline Matrix<M,N> GenerateRandomMatrix()
    {
      Matrix<M,N> randMatrix;
      for (int i=0; i<M; i++) {
        for (int j=0; j<N; j++) {
          randMatrix(i,j) = MyFloat(rand()%RAND_MAX)/MyFloat(RAND_MAX);
        }
      }
      return randMatrix;
    }


    // Operators
    //---------------------------------------------------------------------------------------------
    void operator= (const Matrix<M,N> &other)
    {
      for (int i=0; i<M; i++) {
        for (int j=0; j<N; j++) {
          data[i][j] = other(i,j);
        }
      }
    }

    void operator= (const MyFloat value)
    {
      assert(std::isfinite(value));
      for (int i=0; i<M; i++) {
        for (int j=0; j<N; j++) {
          data[i][j] = value;
        }
      }
    }

    template <int P>
    Matrix<M,P> operator* (const Matrix<N,P> &other) const
    {
      Matrix<M,P> product;
      for (int i=0; i<M; i++) {
        for (int j=0; j<P; j++) {
          for (int k=0; k<N; k++) {
            product(i,j) += data[i][k]*other(k,j);
          }
        }
      }
      return product;
    }

    MyFloat& operator() (const int rindex, const int cindex)
    {
      assert(rindex >= 0 && rindex < M);
      assert(cindex >= 0 && cindex < N);
      return data[rindex][cindex];
    }

    const MyFloat& operator() (const int rindex, const int cindex) const
    {
      assert(rindex >= 0 && rindex < M);
      assert(cindex >= 0 && cindex < N);
      return data[rindex][cindex];
    }

    Matrix<M,N> operator+ (const Matrix<M,N> &other)
    {
      Matrix<M,N> result;
      for (int i=0; i<M; i++) {
        for (int j=0; j<N; j++) {
          result(i,j) = this->data[i][j] + other(i,j);
        }
      }
      return result;
    }

    Matrix<M,N> operator- (const Matrix<M,N> &other)
    {
      Matrix<M,N> result;
      for (int i=0; i<M; i++) {
        for (int j=0; j<N; j++) {
          result(i,j) = this->data[i][j] - other(i,j);
        }
      }
      return result;
    }

    Matrix<M,N> operator* (const MyFloat scalar)
    {
      Matrix<M,N> result;
      for (int i=0; i<M; i++) {
        for (int j=0; j<N; j++) {
          result(i,j) = this->data[i][j]*scalar;
        }
      }
      return result;
    }

    friend std::ostream& operator<< (std::ostream &out, const Matrix<M,N> &matrix)
    {
      for (int i=0; i<M; i++) {
        for (int j=0; j<N; j++) {
          out << matrix(i,j) << " ";
        }
      }
      return out;
    }

    // Other functions
    //---------------------------------------------------------------------------------------------
    inline MyFloat GetSum() const
    {
      MyFloat sum = MyFloat(0.0);
      for (int i=0; i<M; i++) {
        for (int j=0; j<N; j++) {
          sum += data[i][j];
        }
      }
      return sum;
    }

    inline Matrix<N,M> GetTranspose() const
    {
      Matrix<N,M> transpose;
      for (int i=0; i<M; i++) {
        for (int j=0; j<N; j++) {
          transpose(j,i) = data[i][j];
        }
      }
      return transpose;
    }

    //virtual MyFloat GetDeterminant() const = 0;
    //virtual void InvertMatrix() = 0;

    virtual MyFloat GetDeterminant() const
    {
      std::cout << "Can't get determinant for general Matrix object" << std::endl;
      exit(0);
      return MyFloat(0.0);
    }

    virtual void InvertMatrix()
    {
      std::cout << "Can't invert general Matrix object" << std::endl;
      // PANIC?
      exit(0);
    }

  };


  //===============================================================================================
  /// \brief   General matrix class for storing data in column-major order (similar to C/C++ arrays)
  /// \author  S. Heigl & D. A. Hubber
  /// \date    12/09/2018
  //===============================================================================================
  template <int N>
  class SquareMatrix : public Matrix<N,N>
  {
  private:

    using Matrix<N,N>::data;


  public:

    // Constructors and destructor
    //---------------------------------------------------------------------------------------------
    SquareMatrix() : Matrix<N,N>() {};
    SquareMatrix(const SquareMatrix<N> &other) : Matrix<N,N>(other) {};
    SquareMatrix(const Matrix<N,N> &other) : Matrix<N,N>(other) {};
    SquareMatrix(const MyFloat value) : Matrix<N,N>(value) {};
    virtual ~SquareMatrix() {};


    // Operators
    //---------------------------------------------------------------------------------------------
    friend SquareMatrix<N> operator* (const SquareMatrix<N> &mat1, const SquareMatrix<N> &mat2)
    {
      SquareMatrix<N> product;
      for (int i=0; i<N; i++) {
        for (int j=0; j<N; j++) {
          for (int k=0; k<N; k++) {
            product(i,j) += mat1(i,k)*mat2(k,j);
          }
        }
      }
      return product;
    }

    friend Vector<N> operator* (const SquareMatrix<N> &mat, const Vector<N> &vec)
    {
      Vector<N> result;
      for (int i=0; i<N; i++) {
        for (int j=0; j<N; j++) {
          result[i] += mat(i,j)*vec[j];
        }
      }
      return result;
    }

    // Other inlined functions for generating inverse matrices
    //---------------------------------------------------------------------------------------------
    virtual MyFloat GetDeterminant() const override
    {
      if (N == 1) {
        return data[0][0];
      }
      else if (N == 2) {
        return data[0][0]*data[1][1] - data[0][1]*data[1][0];
      }
      else if (N == 3) {
        return data[0][0]*(data[1][1]*data[2][2] - data[2][1]*data[1][2]) -
               data[0][1]*(data[1][0]*data[2][2] - data[1][2]*data[2][0]) +
               data[0][2]*(data[1][0]*data[2][1] - data[1][1]*data[2][0]);
      }
    }

    virtual void InvertMatrix() override
    {
      MyFloat invdet;
      const MyFloat determinant = SquareMatrix<N>::GetDeterminant();

      if (fabs(determinant) > MyFloat(1.0e-30)) {
        invdet = MyFloat(1.0)/determinant;
      }
      else {
        invdet = MyFloat(0.0);
        std::cout << "Error!  Zero determinant for SquareMatrix" << std::endl;
        //exit(0);
        assert(std::isnormal(determinant));
      }

      SquareMatrix<N> invMatrix;
      if (N == 1) {
        invMatrix(0,0) = invdet;
      }
      else if (N == 2) {
        invMatrix(0,0) = invdet*data[1][1];
        invMatrix(0,1) = -MyFloat(1.0)*invdet*data[0][1];
        invMatrix(1,0) = -MyFloat(1.0)*invdet*data[1][0];
        invMatrix(1,1) = invdet*data[0][0];
      }
      else if (N == 3) {
        invMatrix(0,0) = (data[1][1]*data[2][2] - data[2][1]*data[1][2])*invdet;
        invMatrix(0,1) = (data[0][2]*data[2][1] - data[0][1]*data[2][2])*invdet;
        invMatrix(0,2) = (data[0][1]*data[1][2] - data[0][2]*data[1][1])*invdet;
        invMatrix(1,0) = (data[1][2]*data[2][0] - data[1][0]*data[2][2])*invdet;
        invMatrix(1,1) = (data[0][0]*data[2][2] - data[0][2]*data[2][0])*invdet;
        invMatrix(1,2) = (data[1][0]*data[0][2] - data[0][0]*data[1][2])*invdet;
        invMatrix(2,0) = (data[1][0]*data[2][1] - data[2][0]*data[1][1])*invdet;
        invMatrix(2,1) = (data[2][0]*data[0][1] - data[0][0]*data[2][1])*invdet;
        invMatrix(2,2) = (data[0][0]*data[1][1] - data[1][0]*data[0][1])*invdet;
      }
      for (int i=0; i<N; i++) {
        for (int j=0; j<N; j++) {
          data[i][j] = invMatrix(i,j);
        }
      }
      return;
    }

    inline SquareMatrix<N> GetInverseMatrix()
    {
      SquareMatrix<N> inverseMatrix(*this);
      inverseMatrix.InvertMatrix();
      return inverseMatrix;
    }

    inline MyFloat GetTrace() const
    {
      MyFloat trace = MyFloat(0.0);
      for (int i=0; i<N; i++) trace += data[i][i];
      return trace;
    }

    inline SquareMatrix<N> GetTranspose() const
    {
      SquareMatrix<N> transpose;
      for (int i=0; i<N; i++) {
        for (int j=0; j<N; j++) {
          transpose(j,i) = data[i][j];
        }
      }
      return transpose;
    }

    // Static functions
    //---------------------------------------------------------------------------------------------
    inline static SquareMatrix<N> GetInverseMatrixFromMatrix(const SquareMatrix<N> matrix)
    {
      SquareMatrix<N> inverseMatrix(matrix);
      inverseMatrix.InvertMatrix();
      return inverseMatrix;
    }

    static inline SquareMatrix<N> GenerateIdentitySquareMatrix()
    {
      SquareMatrix<N> identityMatrix;
      for (int i=0; i<N; i++) identityMatrix(i,i) = MyFloat(1.0);
      return identityMatrix;
    }

    static inline SquareMatrix<N> GenerateRandomMatrix()
    {
      SquareMatrix<N> randMatrix;
      for (int i=0; i<N; i++) {
        for (int j=0; j<N; j++) {
          randMatrix(i,j) = MyFloat(rand()%RAND_MAX)/MyFloat(RAND_MAX);
        }
      }
      return randMatrix;
    }

  };


  //===============================================================================================
  /// \brief   Specialised child matrix class for rotation operations on vectors only.
  /// \author  S. Heigl & D. A. Hubber
  /// \date    12/09/2018
  //===============================================================================================
  template <int ndim>
  class RotationMatrix : public SquareMatrix<ndim>
  {
  private:

    using Matrix<ndim,ndim>::data;


  public:

    // Constructors and destructor
    //---------------------------------------------------------------------------------------------
    RotationMatrix() : SquareMatrix<ndim>() {};
    RotationMatrix(const RotationMatrix<ndim> &other) : SquareMatrix<ndim>(other) {};
    RotationMatrix(const SquareMatrix<ndim> &other) : SquareMatrix<ndim>(other) {};
    RotationMatrix(const MyFloat value) : SquareMatrix<ndim>(value) {};
    virtual ~RotationMatrix() {};


    // Operators
    //---------------------------------------------------------------------------------------------
    FluidVector<ndim> operator* (const FluidVector<ndim> &fvOrig) const
    {
      FluidVector<ndim> fv(fvOrig);
      for (int k=0; k<ndim; k++) fv[k] = MyFloat(0.0);
      for (int i=0; i<ndim; i++) {
        for (int j=0; j<ndim; j++) {
          fv[i] += data[i][j]*fvOrig[j];
        }
      }
      return fv;
    }


    // Other functions
    //---------------------------------------------------------------------------------------------
    inline RotationMatrix<ndim> GetTranspose() const
    {
      RotationMatrix<ndim> transpose;
      for (int i=0; i<ndim; i++) {
        for (int j=0; j<ndim; j++) {
          transpose(j,i) = data[i][j];
        }
      }
      return transpose;
    }

    inline RotationMatrix<ndim> GetInverseMatrix() const
    {
      return RotationMatrix<ndim>::GetTranspose();
    }


    // Static functions
    //---------------------------------------------------------------------------------------------
    inline static RotationMatrix<1> GetRotationMatrixFromVector(const Vector<1> &v)
    {
      assert(std::isnormal(v.GetMagnitudeSquared()));
      RotationMatrix<1> mat;
      mat(0,0) = v[0];
      return mat;
    }

    inline static RotationMatrix<2> GetRotationMatrixFromVector(const Vector<2> &v)
    {
      assert(std::isnormal(v.GetMagnitudeSquared()));
      RotationMatrix<2> mat;
      mat(0,0) = v[0];
      mat(0,1) = -v[1];
      mat(1,0) = v[1];
      mat(1,1) = v[0];
      return mat;
    }

    inline static RotationMatrix<3> GetRotationMatrixFromVector(const Vector<3> &v)
    {
      assert(std::isnormal(v.GetMagnitudeSquared()));
      RotationMatrix<3> mat;
      const MyFloat radphi = sqrt(v[0]*v[0] + v[1]*v[1]);
      MyFloat cosphi;
      MyFloat sinphi;
      if (radphi == MyFloat(0.0)) {
        cosphi = MyFloat(1.0);
        sinphi = MyFloat(0.0);
      }
      else {
        cosphi = v[0]/radphi;
        sinphi = v[1]/radphi;
      }
      mat(0,0) = radphi*cosphi;
      mat(0,1) = -sinphi;
      mat(0,2) = -v[2]*cosphi;
      mat(1,0) = radphi*sinphi;
      mat(1,1) = cosphi;
      mat(1,2) = -v[2]*sinphi;
      mat(2,0) = v[2];
      mat(2,1) = MyFloat(0.0);
      mat(2,2) = radphi;
      return mat;
    }

    inline static RotationMatrix<1> GetInverseRotationMatrixFromVector(const Vector<1> &v)
    {
      assert(std::isnormal(v.GetMagnitudeSquared()));
      RotationMatrix<1> mat;
      mat(0,0) = v[0];
      return mat;
    }

    inline static RotationMatrix<2> GetInverseRotationMatrixFromVector(const Vector<2> &v)
    {
      assert(std::isnormal(v.GetMagnitudeSquared()));
      RotationMatrix<2> mat;
      mat(0,0) = v[0];
      mat(0,1) = v[1];
      mat(1,0) = -v[1];
      mat(1,1) = v[0];
      return mat;
    }

    inline static RotationMatrix<3> GetInverseRotationMatrixFromVector(const Vector<3> &v)
    {
      assert(std::isnormal(v.GetMagnitudeSquared()));
      RotationMatrix<3> mat;
      const MyFloat radphi = sqrt(v[0]*v[0] + v[1]*v[1]);
      MyFloat cosphi;
      MyFloat sinphi;
      if (radphi == MyFloat(0.0)) {
        cosphi = MyFloat(1.0);
        sinphi = MyFloat(0.0);
      }
      else {
        cosphi = v[0]/radphi;
        sinphi = v[1]/radphi;
      }
      mat(0,0) = radphi*cosphi;
      mat(0,1) = radphi*sinphi;
      mat(0,2) = v[2];
      mat(1,0) = -sinphi;
      mat(1,1) = cosphi;
      mat(1,2) = MyFloat(0.0);
      mat(2,0) = -v[2]*cosphi;
      mat(2,1) = -v[2]*sinphi;
      mat(2,2) = radphi;
      return mat;
    }

  };

}
//-------------------------------------------------------------------------------------------------
#endif
