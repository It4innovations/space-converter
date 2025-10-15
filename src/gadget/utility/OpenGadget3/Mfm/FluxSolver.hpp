#ifndef _OPENGADGET3_FLUX_SOLVER_HPP_
#define _OPENGADGET3_FLUX_SOLVER_HPP_


#include <cstdio>
#include <cstdlib>
#include <string>
#include "FluidVector.hpp"
#include "RiemannSolver.hpp"
#include "../CodeBase/precision.h"


//-------------------------------------------------------------------------------------------------
namespace OpenGadget3
{

  //===============================================================================================
  /// \brief   Class to compute the Godunov flux using the chosen Riemann solver algorithm.
  /// \author  S. Heigl & D. A. Hubber
  /// \date    13/02/2018
  //===============================================================================================
  template <int ndim>
  class FluxSolver
  {
  private:

    Eos<ndim>* const eos;
    ExactRiemannSolver<ndim> *exactRiemann;
    HllRiemannSolver<ndim> *hllRiemann;
    HllcRiemannSolver<ndim> *hllcRiemann;
    RoeRiemannSolver<ndim> *roeRiemann;


  public:

    static const int ivx = 0;
    static const int ivy = 1;
    static const int ivz = 2;
    static const int irho = ndim;
    static const int ipress = ndim + 1;

    FluxSolver(Eos<ndim>* const);
    ~FluxSolver();

    FluxVector<ndim> ComputeFluxes(const bool, const Vector<ndim> &, Vector<ndim> &,
                                   WFluidVector<ndim> &, WFluidVector<ndim> &);

  };

}
//-------------------------------------------------------------------------------------------------
#endif
