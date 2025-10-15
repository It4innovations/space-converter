#ifndef _OPENGADGET3_RIEMANN_SOLVER_HPP_
#define _OPENGADGET3_RIEMANN_SOLVER_HPP_


#include "../System/assert.h"
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <math.h>
#include "Eos.hpp"
#include "FluidVector.hpp"
#include "Vector.hpp"
#include "Matrix.hpp"
#include "../CodeBase/precision.h"


//-------------------------------------------------------------------------------------------------
namespace OpenGadget3
{

  //===============================================================================================
  /// \brief   Base Riemann solver class.
  /// \author  D. A. Hubber
  /// \date    13/02/2018
  //===============================================================================================
  template <int ndim>
  class RiemannSolver
  {
  protected:

    const MyFloat gamma;
    const MyFloat invgamma;
    const MyFloat g1;
    const MyFloat g2;
    const MyFloat g3;
    const MyFloat g4;
    const MyFloat g5;
    const MyFloat g6;
    const MyFloat g7;
    const MyFloat g8;
    const MyFloat g9;
    Eos<ndim>* const eos;


  public:

    static const int ivx = 0;
    static const int ivy = 1;
    static const int ivz = 2;
    static const int irho = ndim;
    static const int ipress = ndim + 1;
    static const int nvar = ndim + 2;

    //---------------------------------------------------------------------------------------------
    RiemannSolver(Eos<ndim>* const);
    virtual ~RiemannSolver() {};

    //---------------------------------------------------------------------------------------------
    virtual bool ComputeFluxes(const bool, const Vector<ndim> &, WFluidVector<ndim> &,
                               WFluidVector<ndim> &, Vector<ndim> &, FluxVector<ndim> &) const = 0;

    // Some useful inline functions
    //---------------------------------------------------------------------------------------------
    inline FluxVector<ndim> ComputeFlux1d(const WFluidVector<ndim> &W) const
    {
      FluxVector<ndim> flux;
      for (int k=0; k<ndim; k++) flux[k] = W[irho]*W[ivx]*W[k];
      flux[ivx]    = W[irho]*W[ivx]*W[ivx] + W[ipress];
      flux[irho]   = W[irho]*W[ivx];
      flux[ipress] = W[ivx]*(W[irho]*W.GetSpecificTotalEnergy(gamma) + W[ipress]);
      return flux;
    }

    inline FluxVector<ndim> ComputeFlux1dAtFace
     (const WFluidVector<ndim> &W,
      const MyFloat vFace) const
    {
      FluxVector<ndim> flux;
      for (int k=0; k<ndim; k++) flux[k] = W[irho]*(W[ivx] - vFace)*W[k];
      flux[ivx]    = W[irho]*(W[ivx] - vFace)*W[ivx] + W[ipress];
      flux[irho]   = W[irho]*(W[ivx] - vFace);
      flux[ipress] = (W[ivx] - vFace)*W[irho]*W.GetSpecificTotalEnergy(gamma) + W[ivx]*W[ipress];
      return flux;
    }

  };



  //===============================================================================================
  /// \brief   Exact Riemann solver, based on algorithm described in Toro (1999).
  /// \author  D. A. Hubber
  /// \date    13/02/2018
  //===============================================================================================
  template <int ndim>
  class ExactRiemannSolver final : public RiemannSolver<ndim>
  {
  protected:

    using RiemannSolver<ndim>::gamma;
    using RiemannSolver<ndim>::invgamma;
    using RiemannSolver<ndim>::g1;
    using RiemannSolver<ndim>::g2;
    using RiemannSolver<ndim>::g3;
    using RiemannSolver<ndim>::g4;
    using RiemannSolver<ndim>::g5;
    using RiemannSolver<ndim>::g6;
    using RiemannSolver<ndim>::g7;
    using RiemannSolver<ndim>::g8;
    using RiemannSolver<ndim>::g9;
    using RiemannSolver<ndim>::eos;

    const int maxNumIterations;
    const MyFloat tolerance;


  public:

    static const int ivx = 0;
    static const int ivy = 1;
    static const int ivz = 2;
    static const int irho = ndim;
    static const int ipress = ndim + 1;
    static const int nvar = ndim + 2;

    //---------------------------------------------------------------------------------------------
    ExactRiemannSolver(Eos<ndim>* const _eos, const int _maxNumIterations, const MyFloat _tolerance) :
     RiemannSolver<ndim>(_eos), maxNumIterations(_maxNumIterations), tolerance(_tolerance) {};
    ~ExactRiemannSolver() {};

    //---------------------------------------------------------------------------------------------
    virtual bool ComputeFluxes(const bool, const Vector<ndim> &, WFluidVector<ndim> &,
                               WFluidVector<ndim> &, Vector<ndim> &,
                               FluxVector<ndim> &) const override;

    void ComputeStarRegion(const MyFloat, const MyFloat, const MyFloat, const MyFloat,
                           const MyFloat, const MyFloat, const MyFloat, const MyFloat,
                           MyFloat &, MyFloat &) const;

    void SampleExactSolution(const MyFloat, const MyFloat, const MyFloat, const MyFloat,
                             const MyFloat, const MyFloat, const MyFloat, const MyFloat,
                             const MyFloat, const MyFloat, const MyFloat,
                             MyFloat &, MyFloat &, MyFloat &) const;
  };



  //===============================================================================================
  /// \brief   HLL approximate Riemann solver
  /// \author  D. A. Hubber
  /// \date    13/02/2018
  //===============================================================================================
  template <int ndim>
  class HllRiemannSolver final : public RiemannSolver<ndim>
  {
  protected:

    using RiemannSolver<ndim>::gamma;
    using RiemannSolver<ndim>::invgamma;
    using RiemannSolver<ndim>::g1;
    using RiemannSolver<ndim>::g2;
    using RiemannSolver<ndim>::g3;
    using RiemannSolver<ndim>::g4;
    using RiemannSolver<ndim>::g5;
    using RiemannSolver<ndim>::g6;
    using RiemannSolver<ndim>::g7;
    using RiemannSolver<ndim>::g8;
    using RiemannSolver<ndim>::g9;
    using RiemannSolver<ndim>::eos;


  public:

    static const int ivx = 0;
    static const int ivy = 1;
    static const int ivz = 2;
    static const int irho = ndim;
    static const int ipress = ndim + 1;
    static const int nvar = ndim + 2;

    //---------------------------------------------------------------------------------------------
    HllRiemannSolver(Eos<ndim>* const _eos) : RiemannSolver<ndim>(_eos) {};
    ~HllRiemannSolver() {};

    //---------------------------------------------------------------------------------------------
    virtual bool ComputeFluxes(const bool, const Vector<ndim> &, WFluidVector<ndim> &,
                               WFluidVector<ndim> &, Vector<ndim> &,
                               FluxVector<ndim> &) const override;

    void ComputeWaveSpeedEstimates(const WFluidVector<ndim> &, const WFluidVector<ndim> &,
                                   MyFloat &, MyFloat &) const;
    MyFloat ComputeStarVelocity(const WFluidVector<ndim> &, const WFluidVector<ndim> &,
                                const MyFloat, const MyFloat) const;
    UFluidVector<ndim> ComputeHllStateVector(const WFluidVector<ndim> &, const WFluidVector<ndim> &,
                                             const MyFloat, const MyFloat) const;
    FluxVector<ndim> ComputeHllFlux1d(const WFluidVector<ndim> &, const WFluidVector<ndim> &,
                                      const MyFloat, const MyFloat) const;

  };

  //===============================================================================================
  /// \brief   HLLC approximate Riemann solver
  /// \author  F. Groth, based on HLL solver of D. A. Hubber
  /// \date    02/12/2022
  //===============================================================================================
  template <int ndim>
  class HllcRiemannSolver final : public RiemannSolver<ndim>
  {
  protected:

    using RiemannSolver<ndim>::gamma;
    using RiemannSolver<ndim>::invgamma;
    using RiemannSolver<ndim>::g1;
    using RiemannSolver<ndim>::g2;
    using RiemannSolver<ndim>::g3;
    using RiemannSolver<ndim>::g4;
    using RiemannSolver<ndim>::g5;
    using RiemannSolver<ndim>::g6;
    using RiemannSolver<ndim>::g7;
    using RiemannSolver<ndim>::g8;
    using RiemannSolver<ndim>::g9;
    using RiemannSolver<ndim>::eos;


  public:

    static const int ivx = 0;
    static const int ivy = 1;
    static const int ivz = 2;
    static const int irho = ndim;
    static const int ipress = ndim + 1;
    static const int nvar = ndim + 2;

    //---------------------------------------------------------------------------------------------
    HllcRiemannSolver(Eos<ndim>* const _eos) : RiemannSolver<ndim>(_eos) {};
    ~HllcRiemannSolver() {};

    //---------------------------------------------------------------------------------------------
    virtual bool ComputeFluxes(const bool, const Vector<ndim> &, WFluidVector<ndim> &,
                               WFluidVector<ndim> &, Vector<ndim> &,
                               FluxVector<ndim> &) const override;

    void ComputeWaveSpeedEstimates(const WFluidVector<ndim> &, const WFluidVector<ndim> &,
                                   MyFloat &, MyFloat &) const;
    MyFloat ComputeStarVelocity(const WFluidVector<ndim> &, const WFluidVector<ndim> &,
                                const MyFloat, const MyFloat) const;
    UFluidVector<ndim> ComputeHllcStateVector(const WFluidVector<ndim> &, const WFluidVector<ndim> &,
                                             const MyFloat, const MyFloat) const;
    FluxVector<ndim> ComputeHllcFlux1d(const WFluidVector<ndim> &, const WFluidVector<ndim> &,
				       const MyFloat, const MyFloat, const MyFloat) const;

  };


  //===============================================================================================
  /// \brief   Roe approximate Riemann solver
  /// \author  F. Groth
  /// \date    02/12/2022
  //===============================================================================================
  template <int ndim>
  class RoeRiemannSolver final : public RiemannSolver<ndim>
  {
  protected:

    using RiemannSolver<ndim>::gamma;
    using RiemannSolver<ndim>::invgamma;
    using RiemannSolver<ndim>::g1;
    using RiemannSolver<ndim>::g2;
    using RiemannSolver<ndim>::g3;
    using RiemannSolver<ndim>::g4;
    using RiemannSolver<ndim>::g5;
    using RiemannSolver<ndim>::g6;
    using RiemannSolver<ndim>::g7;
    using RiemannSolver<ndim>::g8;
    using RiemannSolver<ndim>::g9;
    using RiemannSolver<ndim>::eos;


  public:

    static const int ivx = 0;
    static const int ivy = 1;
    static const int ivz = 2;
    static const int irho = ndim;
    static const int ipress = ndim + 1;
    static const int nvar = ndim + 2;

    //---------------------------------------------------------------------------------------------
    RoeRiemannSolver(Eos<ndim>* const _eos) : RiemannSolver<ndim>(_eos) {};
    ~RoeRiemannSolver() {};

    //---------------------------------------------------------------------------------------------
    virtual bool ComputeFluxes(const bool, const Vector<ndim> &, WFluidVector<ndim> &,
                               WFluidVector<ndim> &, Vector<ndim> &,
                               FluxVector<ndim> &) const override;

    void ComputeAveragedValues(const WFluidVector<ndim> &, const WFluidVector<ndim> &,
			       Vector<nvar> &) const;
    void ComputeEigenvalues(const Vector<nvar> &,
                                Vector<nvar> &) const;
    void ComputeRightEigenvectors(const Vector<nvar> &,
					  SquareMatrix<nvar> &) const;
    void ComputeWaveStrength(const WFluidVector<ndim> &, const WFluidVector<ndim> &, const Vector<nvar>,
			     Vector<nvar> &) const;
    void ComputeStarVelocities(const WFluidVector<ndim> &, const WFluidVector<ndim> &,
			       const Vector<nvar> &, const Vector<nvar> &, const SquareMatrix<nvar> &,
			       MyFloat &, MyFloat &, MyFloat &, MyFloat &) const;
    FluxVector<ndim> ComputeRoeFlux1d(const WFluidVector<ndim> &, const WFluidVector<ndim> &,
				       const Vector<nvar> &, const Vector<nvar> &, const SquareMatrix<nvar> &) const;
  };


}
//-------------------------------------------------------------------------------------------------
#endif
