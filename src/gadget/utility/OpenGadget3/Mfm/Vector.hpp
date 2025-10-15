#ifndef _OPENGADGET3_VECTOR_HPP_
#define _OPENGADGET3_VECTOR_HPP_


#include "../System/assert.h"
#include <cmath>
#include <cstdlib>
#include <ostream>
#include "../CodeBase/precision.h"


//-------------------------------------------------------------------------------------------------
namespace OpenGadget3
{

  //===============================================================================================
  /// \brief   Basic 1d/2d/3d vector class for simplifying storage and calculations with vectors.
  /// \author  D. A. Hubber
  /// \date    27/10/2017
  //===============================================================================================
  template <int ndim>
  class Vector
  {
  private:

    MyFloat data[ndim];                        ///< Basic array containing vector components


  public:

    // Constructors and destructor
    //---------------------------------------------------------------------------------------------
    Vector()
    {
      for (int k=0; k<ndim; k++) data[k] = MyFloat(0.0);
    }

    Vector(const MyFloat value)
    {
      for (int k=0; k<ndim; k++) data[k] = value;
    }

    Vector(const MyFloat v0, const MyFloat v1)
    {
      assert(ndim == 2);
      data[0] = v0;
      data[1] = v1;
    }

    Vector(const MyFloat v0, const MyFloat v1, const MyFloat v2)
    {
      assert(ndim == 3);
      data[0] = v0;
      data[1] = v1;
      data[2] = v2;
    }

    Vector(const Vector<ndim> &other)
    {
      for (int k=0; k<ndim; k++) data[k] = other[k];
    }

    ~Vector() {};


    // Operators
    //---------------------------------------------------------------------------------------------
    inline void operator= (const Vector<ndim> &other)
    {
      for (int k=0; k<ndim; k++) data[k] = other[k];
    }

    inline void operator= (const MyFloat value)
    {
      assert(std::isfinite(value));
      for (int k=0; k<ndim; k++) data[k] = value;
    }

    inline MyFloat& operator[] (const int index)
    {
      assert(index >= 0 && index < ndim);
      return data[index];
    }

    inline const MyFloat& operator[] (const int index) const
    {
      assert(index >= 0 && index < ndim);
      return data[index];
    }

    inline Vector<ndim> operator+ (const Vector<ndim> &otherVector) const
    {
      Vector<ndim> result;
      for (int k=0; k<ndim; k++) result[k] = this->data[k] + otherVector[k];
      return result;
    }

    inline Vector<ndim> operator- (const Vector<ndim> &otherVector) const
    {
      Vector<ndim> result;
      for (int k=0; k<ndim; k++) result[k] = this->data[k] - otherVector[k];
      return result;
    }

    inline void operator+= (const Vector<ndim> &otherVector)
    {
      for (int k=0; k<ndim; k++) data[k] += otherVector[k];
    }

    inline void operator-= (const Vector<ndim> &otherVector)
    {
      for (int k=0; k<ndim; k++) data[k] -= otherVector[k];
    }

    inline Vector<ndim> operator* (const MyFloat multFac)
    {
      assert(std::isfinite(multFac));
      Vector<ndim> result;
      for (int k=0; k<ndim; k++) result[k] = this->data[k]*multFac;
      return result;
    }

    friend inline Vector<ndim> operator* (const MyFloat multFac, const Vector<ndim> &vec)
    {
      assert(std::isfinite(multFac));
      Vector<ndim> result;
      for (int k=0; k<ndim; k++) result[k] = vec[k]*multFac;
      return result;
    }

    friend inline Vector<ndim> operator* (const Vector<ndim> &vec, const MyFloat multFac)
    {
      assert(std::isfinite(multFac));
      Vector<ndim> result;
      for (int k=0; k<ndim; k++) result[k] = vec[k]*multFac;
      return result;
    }

    void operator*= (const MyFloat multFac)
    {
      assert(std::isfinite(multFac));
      for (int k=0; k<ndim; k++) data[k] *= multFac;
    }

    void operator/= (const MyFloat divFac)
    {
      assert(std::fabs(divFac) > MyFloat(0.0));
      const MyFloat invFac = MyFloat(1.0) / divFac;
      for (int k=0; k<ndim; k++) data[k] *= invFac;
    }

    inline friend std::ostream& operator<< (std::ostream &out, const Vector<ndim> &vector)
    {
      for (int k=0; k<ndim; k++) out << vector[k] << " ";
      //out << vector[ndim-1];
      return out;
    }


    // Other functions
    //---------------------------------------------------------------------------------------------
    inline MyFloat GetMagnitudeSquared() const
    {
      MyFloat magSqd = MyFloat(0.0);
      for (int k=0; k<ndim; k++) magSqd += data[k]*data[k];
      return magSqd;
    }

    inline MyFloat GetMagnitude() const
    {
      return sqrt(GetMagnitudeSquared());
    }

    inline Vector<ndim> GetUnitVector() const
    {
      Vector<ndim> unitVector(*this);
      unitVector.Normalise();
      return unitVector;
    }

    inline void Normalise()
    {
      const MyFloat magSqd = GetMagnitudeSquared();
      if (magSqd > MyFloat(0.0)) {
        const MyFloat invMag = MyFloat(1.0) / sqrt(magSqd);
        for (int k=0; k<ndim; k++) data[k] *= invMag;
      }
    }


    // Static functions between two vectors
    //---------------------------------------------------------------------------------------------
    inline static MyFloat DotProduct(const Vector<1> &v1, const Vector<1> &v2)
    {
      return v1[0]*v2[0];
    }

    inline static MyFloat DotProduct(const Vector<2> &v1, const Vector<2> &v2)
    {
      return v1[0]*v2[0] + v1[1]*v2[1];
    }

    inline static MyFloat DotProduct(const Vector<3> &v1, const Vector<3> &v2)
    {
      return v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2];
    }

    inline static MyFloat DistanceSquared(const Vector<1> &v1, const Vector<1> &v2)
    {
      return (v2[0] - v1[0])*(v2[0] - v1[0]);
    }

    inline static MyFloat DistanceSquared(const Vector<2> &v1, const Vector<2> &v2)
    {
      return (v2[0] - v1[0])*(v2[0] - v1[0]) + (v2[1] - v1[1])*(v2[1] - v1[1]);
    }

    inline static MyFloat DistanceSquared(const Vector<3> &v1, const Vector<3> &v2)
    {
      return (v2[0] - v1[0])*(v2[0] - v1[0]) +
       (v2[1] - v1[1])*(v2[1] - v1[1]) + (v2[2] - v1[2])*(v2[2] - v1[2]);
    }

    inline static MyFloat Distance(const Vector<ndim> &v1, const Vector<ndim> &v2)
    {
      return sqrt(DistanceSquared(v1, v2));
    }


    // Static utilty function for generating random vectors (for testing)
    //---------------------------------------------------------------------------------------------
    static inline Vector<ndim> GenerateRandomVector()
    {
      Vector<ndim> v;
      for (int k=0; k<ndim; k++) v[k] = MyFloat(rand()%RAND_MAX)/MyFloat(RAND_MAX);
      return v;
    }

    static inline Vector<ndim> GenerateRandomUnitVector()
    {
      Vector<ndim> v;
      for (int k=0; k<ndim; k++) v[k] = MyFloat(2.0)*MyFloat(rand()%RAND_MAX)/MyFloat(RAND_MAX) - MyFloat(1.0);
      v.Normalise();
      return v;
    }

  };


}
//-------------------------------------------------------------------------------------------------
#endif
