#ifndef _OPENGADGET3_FLUID_VECTOR_HPP_
#define _OPENGADGET3_FLUID_VECTOR_HPP_


#include "../System/assert.h"
#include <cmath>
#include <cstdlib>
#include <math.h>
#include <ostream>
//#include "../CodeBase/allvars.h"
#include "../CodeBase/precision.h"
#include "Vector.hpp"


//-------------------------------------------------------------------------------------------------
namespace OpenGadget3
{

  //===============================================================================================
  /// \brief   General fluid vector class.
  /// \author  S. Heigl & D. A. Hubber
  /// \date    13/02/2018
  //===============================================================================================
  template <int ndim>
  class FluidVector
  {
  protected:

    MyFloat data[ndim+2];                        ///< Fluid vector in various forms (e.g. W, Q, U)


  public:

    static const int ivx = 0;
    static const int ivy = 1;
    static const int ivz = 2;
    static const int irho = ndim;
    static const int ipress = ndim + 1;
    static const int nvar = ndim + 2;


    // Constructors and destructor
    //---------------------------------------------------------------------------------------------
    FluidVector()
    {
      for (int k=0; k<ndim+2; k++) data[k] = MyFloat(0.0);
    }
    FluidVector(const FluidVector<ndim> &other)
    {
      for (int k=0; k<ndim+2; k++) data[k] = other[k];
    }
    FluidVector(const Vector<ndim> vec, const MyFloat s1, const MyFloat s2)
    {
      for (int k=0; k<ndim; k++) data[k] = vec[k];
      data[irho]   = s1;
      data[ipress] = s2;
    }
    ~FluidVector() {};


    // Operators
    //---------------------------------------------------------------------------------------------
    inline void operator= (const FluidVector<ndim> &other)
    {
      for (int k=0; k<ndim+2; k++) data[k] = other[k];
    }

    inline void operator= (const MyFloat value)
    {
      assert(std::isfinite(value));
      for (int k=0; k<ndim+2; k++) data[k] = value;
    }

    inline MyFloat& operator[] (const int index)
    {
      assert(index >= 0 && index < ndim+2);
      return data[index];
    }

    inline const MyFloat& operator[] (const int index) const
    {
      assert(index >= 0 && index < ndim+2);
      return data[index];
    }

    inline FluidVector<ndim> operator+ (const FluidVector<ndim> &other)
    {
      FluidVector<ndim> result;
      for (int k=0; k<ndim+2; k++) result[k] = this->data[k] + other[k];
      return result;
    }

    inline FluidVector<ndim> operator- (const FluidVector<ndim> &otherFluidVector)
    {
      FluidVector<ndim> result;
      for (int k=0; k<ndim+2; k++) result[k] = this->data[k] - otherFluidVector[k];
      return result;
    }

    inline FluidVector<ndim> operator* (const MyFloat multFac)
    {
      assert(std::isfinite(multFac));
      FluidVector<ndim> result;
      for (int k=0; k<ndim+2; k++) result[k] = this->data[k]*multFac;
      return result;
    }

    inline void operator+= (const FluidVector<ndim> &other)
    {
      for (int k=0; k<ndim+2; k++) data[k] += other[k];
    }

    inline void operator-= (const FluidVector<ndim> &other)
    {
      for (int k=0; k<ndim+2; k++) data[k] -= other[k];
    }

    inline friend std::ostream& operator<< (std::ostream &out, const FluidVector<ndim> &FluidVector)
    {
      out << "FluidVector: (" << FluidVector[0];
      for (int k=1; k<ndim+2; k++) out << ", " << FluidVector[k];
      out << ")";
      return out;
    }


    // Static functions
    //---------------------------------------------------------------------------------------------
    static inline FluidVector<ndim> GenerateRandomFluidVector()
    {
      FluidVector<ndim> fv;
      for (int k=0; k<ndim; k++) fv[k] = MyFloat(2.0)*MyFloat(rand()%RAND_MAX)/MyFloat(RAND_MAX) - MyFloat(1.0);
      fv[irho] = MyFloat(rand()%RAND_MAX)/MyFloat(RAND_MAX);
      fv[ipress] = MyFloat(rand()%RAND_MAX)/MyFloat(RAND_MAX);
      return fv;
    }

  };


  template <int ndim>
  class QFluidVector;

  template <int ndim>
  class UFluidVector;


  //===============================================================================================
  /// \brief   Fluid vector of primitive quantities, W = (velocity, density, pressure)
  /// \author  S. Heigl & D. A. Hubber
  /// \date    14/02/2018
  //===============================================================================================
  template <int ndim>
  class WFluidVector final : public FluidVector<ndim>
  {
  protected:

    using FluidVector<ndim>::data;


  public:

    static const int ivx = 0;
    static const int ivy = 1;
    static const int ivz = 2;
    static const int irho = ndim;
    static const int ipress = ndim + 1;
    static const int nvar = ndim + 2;


    // Constructors
    //---------------------------------------------------------------------------------------------
    WFluidVector() : FluidVector<ndim>() {};
    WFluidVector(const WFluidVector<ndim> &other)
    {
      for (int k=0; k<ndim+2; k++) data[k] = other[k];
    }
    WFluidVector(const FluidVector<ndim> &other)
    {
      for (int k=0; k<ndim+2; k++) data[k] = other[k];
    }
    WFluidVector(const Vector<ndim> v, const MyFloat rho, const MyFloat p)
    {
      for (int k=0; k<ndim; k++) data[k] = v[k];
      data[irho]   = rho;
      data[ipress] = p;
    }
    WFluidVector(const MyFloat *v, const MyLongDouble rho, const MyFloat p);
    WFluidVector(const QFluidVector<ndim> &, const MyFloat, const MyFloat, const MyFloat);
    WFluidVector(const UFluidVector<ndim> &, const MyFloat);
    

    // Operators
    //---------------------------------------------------------------------------------------------
    inline void operator= (const MyFloat value)
    {
      assert(std::isfinite(value));
      for (int k=0; k<ndim+2; k++) data[k] = value;
    }


    // Getter functions
    //---------------------------------------------------------------------------------------------
    inline MyFloat GetDensity() const {return data[irho];}
    inline MyFloat GetPressure() const {return data[ipress];}
    inline MyFloat GetSpecificInternalEnergy(const MyFloat gamma) const /* in physical units */
    {
      return data[ipress]/(gamma - MyFloat(1.0))/data[irho];
    }
    inline MyFloat GetSpecificKineticEnergy() const
    {
      MyFloat keTot = MyFloat(0.0);
      for (int k=0; k<ndim; k++) keTot += data[k]*data[k];
      return MyFloat(0.5)*keTot;
    }
    inline MyFloat GetSpecificTotalEnergy(const MyFloat gamma) const
    {
      return GetSpecificKineticEnergy() + GetSpecificInternalEnergy(gamma);
    }
    inline Vector<ndim> GetVelocity() const
    {
      Vector<ndim> v;
      for (int k=0; k<ndim; k++) v[k] = data[k];
      return v;
    }


    // Static functions
    //---------------------------------------------------------------------------------------------
    static WFluidVector<ndim> ComputePrimitiveTimeDerivative(const MyFloat, const WFluidVector<ndim> &,
                                                             const Vector<ndim> *, const Vector<ndim>);

    static inline WFluidVector<ndim> GenerateRandomWFluidVector()
    {
      WFluidVector<ndim> W;
      for (int k=0; k<ndim; k++) W[k] = MyFloat(2.0)*MyFloat(rand()%RAND_MAX)/MyFloat(RAND_MAX) - MyFloat(1.0);
      W[irho] = MyFloat(rand()%RAND_MAX)/MyFloat(RAND_MAX);
      W[ipress] = MyFloat(rand()%RAND_MAX)/MyFloat(RAND_MAX);
      return W;
    }

  };


  //===============================================================================================
  /// \brief   Fluid vector of conserved quantities, Q = (momentum, mass, total energy)
  /// \author  S. Heigl & D. A. Hubber
  /// \date    14/02/2018
  //===============================================================================================
  template <int ndim>
  class UFluidVector final : public FluidVector<ndim>
  {
  protected:

    using FluidVector<ndim>::data;


  public:

    static const int ivx = 0;
    static const int ivy = 1;
    static const int ivz = 2;
    static const int irho = ndim;
    static const int ipress = ndim + 1;
    static const int nvar = ndim + 2;


    // Constructors
    //---------------------------------------------------------------------------------------------
    UFluidVector() : FluidVector<ndim>() {};
    UFluidVector(const UFluidVector<ndim> &other)
    {
      for (int k=0; k<ndim+2; k++) data[k] = other[k];
    }
    UFluidVector(const FluidVector<ndim> &other)
    {
      for (int k=0; k<ndim+2; k++) data[k] = other[k];
    }
    UFluidVector(const QFluidVector<ndim> &Q, const MyFloat vol)
    {
      assert(std::isnormal(vol));
      const MyFloat invVol = MyFloat(1.0) / vol;
      for (int k=0; k<ndim+2; k++) data[k] = Q[k]*invVol;
    }
    UFluidVector(const WFluidVector<ndim> &, const MyFloat);


    // Operators
    //---------------------------------------------------------------------------------------------
    inline void operator= (const MyFloat value)
    {
      assert(std::isfinite(value));
      for (int k=0; k<ndim+2; k++) data[k] = value;
    }


    // Other getter and setter functions
    //---------------------------------------------------------------------------------------------
    MyFloat GetSpecificInternalEnergy() const;

    inline MyFloat GetEnergyDensity() const {return data[ipress];}
    inline MyFloat GetDensity() const {return data[irho];}
    inline Vector<ndim> GetMomentumDensity() const
    {
      Vector<ndim> p;
      for (int k=0; k<ndim; k++) p[k] = data[k];
      return p;
    }
    inline Vector<ndim> GetVelocity() const
    {
      const MyFloat invRho = MyFloat(1.0) / data[irho];
      Vector<ndim> v;
      for (int k=0; k<ndim; k++) v[k] = data[k]*invRho;
      return v;
    }


    // Static functions
    //---------------------------------------------------------------------------------------------
    static inline UFluidVector<ndim> GenerateRandomUFluidVector()
    {
      UFluidVector<ndim> U;
      for (int k=0; k<ndim; k++) U[k] = MyFloat(2.0)*MyFloat(rand()%RAND_MAX)/MyFloat(RAND_MAX) - MyFloat(1.0);
      U[irho] = MyFloat(rand()%RAND_MAX)/MyFloat(RAND_MAX);
      U[ipress] = MyFloat(rand()%RAND_MAX)/MyFloat(RAND_MAX);
      return U;
    }

  };


  //===============================================================================================
  /// \brief   Fluid vector of conserved quantities, Q = (momentum, mass, total energy)
  /// \author  S. Heigl & D. A. Hubber
  /// \date    14/02/2018
  //===============================================================================================
  template <int ndim>
  class QFluidVector final : public FluidVector<ndim>
  {
  protected:

    using FluidVector<ndim>::data;


  public:

    static const int ivx = 0;
    static const int ivy = 1;
    static const int ivz = 2;
    static const int irho = ndim;
    static const int ipress = ndim + 1;
    static const int nvar = ndim + 2;


    // Constructors
    //---------------------------------------------------------------------------------------------
    QFluidVector() : FluidVector<ndim>() {};
    QFluidVector(const QFluidVector<ndim> &other)
    {
      for (int k=0; k<ndim+2; k++) data[k] = other[k];
    }
    QFluidVector(const FluidVector<ndim> &other)
    {
      for (int k=0; k<ndim+2; k++) data[k] = other[k];
    }
    QFluidVector(const Vector<ndim> p, const MyFloat m, const MyFloat E)
    {
      for (int k=0; k<ndim; k++) data[k] = p[k];
      data[irho]   = m;
      data[ipress] = E;
    }
    QFluidVector(const UFluidVector<ndim> &U, const MyFloat vol)
    {
      assert(std::isnormal(vol));
      for (int k=0; k<ndim+2; k++) data[k] = U[k]*vol;
    }
    QFluidVector(const WFluidVector<ndim> &, const MyFloat, const MyFloat, const MyFloat);


    // Operators
    //---------------------------------------------------------------------------------------------
    inline void operator= (const MyFloat value)
    {
      assert(std::isfinite(value));
      for (int k=0; k<ndim+2; k++) data[k] = value;
    }


    // Other getter and setter functions
    //---------------------------------------------------------------------------------------------
    MyFloat GetSpecificInternalEnergy() const;

    inline MyFloat GetEnergy() const {return data[ipress];}
    inline MyFloat GetMass() const {return data[irho];}
    inline Vector<ndim> GetMomentum() const
    {
      Vector<ndim> p;
      for (int k=0; k<ndim; k++) p[k] = data[k];
      return p;
    }
    inline Vector<ndim> GetVelocity() const
    {
      assert(data[irho] > MyFloat(0.0));
      const MyFloat invm = MyFloat(1.0) / data[irho];
      Vector<ndim> v;
      for (int k=0; k<ndim; k++) v[k] = data[k]*invm;
      return v;
    }


    // Static functions
    //---------------------------------------------------------------------------------------------
    static inline QFluidVector<ndim> GenerateRandomQFluidVector()
    {
      QFluidVector<ndim> Q;
      for (int k=0; k<ndim; k++) Q[k] = MyFloat(2.0)*MyFloat(rand()%RAND_MAX)/MyFloat(RAND_MAX) - MyFloat(1.0);
      Q[irho] = MyFloat(rand()%RAND_MAX)/MyFloat(RAND_MAX);
      Q[ipress] = MyFloat(rand()%RAND_MAX)/MyFloat(RAND_MAX);
      return Q;
    }

  };


  template <int ndim>
  using FluxVector = FluidVector<ndim>;

}
//-------------------------------------------------------------------------------------------------
#endif
