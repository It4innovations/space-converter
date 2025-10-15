#include <map>
#include <string>
#include "FluxSolver.hpp"
#include "Matrix.hpp"
#include "Vector.hpp"
#include "../CodeBase/allvars.h"


//-------------------------------------------------------------------------------------------------
namespace OpenGadget3
{

  //===============================================================================================
  /// Constructor for FluxSolver class.
  //===============================================================================================
  template <int ndim>
  FluxSolver<ndim>::FluxSolver(Eos<ndim>* const _eos) : eos(_eos)
  {
    exactRiemann = new ExactRiemannSolver<ndim>(_eos, 8, (MyFloat) 1.0e-4);
    roeRiemann   = new RoeRiemannSolver<ndim>(_eos);
    hllcRiemann  = new HllcRiemannSolver<ndim>(_eos);
    hllRiemann   = new HllRiemannSolver<ndim>(_eos);
  }


  //===============================================================================================
  /// Destructor for FluxSolver class.
  //===============================================================================================
  template <int ndim>
  FluxSolver<ndim>::~FluxSolver()
  {
    delete hllRiemann;
    delete hllcRiemann;
    delete roeRiemann;
    delete exactRiemann;
  }


  //===============================================================================================
  /// Compute the numerical Godunov fluxes for a given Riemann problem with left and right states.
  //===============================================================================================
  template <int ndim>
  FluxVector<ndim> FluxSolver<ndim>::ComputeFluxes
   (const bool zeroMassFlux,                     ///< [in] Compute fluxes in zero-mass flux frame
    const Vector<ndim> &rUnit,                   ///< [in] Unit vector on face
    Vector<ndim> &vFace,                         ///< [inout] Face velocity vector
    WFluidVector<ndim> &Wleft,                   ///< [inout] LHS primitive state
    WFluidVector<ndim> &Wright)                  ///< [inout] RHS primitive state

  {
    FluxVector<ndim> flux;

    // Compute rotation and inverse rotation matrices
    const RotationMatrix<ndim> rotMatrix = RotationMatrix<ndim>::GetRotationMatrixFromVector(rUnit);
    const RotationMatrix<ndim> invRotMatrix = rotMatrix.GetInverseMatrix();

    // Boost velocities to frame of moving frame
    for (int k=0; k<ndim; k++) Wleft[k] -= vFace[k];
    for (int k=0; k<ndim; k++) Wright[k] -= vFace[k];

    // Rotate velocity component so Riemann problem is now a 1d problem set along the x-axis
    Wleft  = invRotMatrix*Wleft;
    Wright = invRotMatrix*Wright;
    
    assert(std::isnormal(Wleft[0]));
    assert(std::isnormal(Wright[0]));

    // First attempt to use the HLL solver; if an error is generated, then use the exact solver
#ifdef RIEMANN_SOLVER_HLL
    if (!hllRiemann->ComputeFluxes(zeroMassFlux, rUnit, Wleft, Wright, vFace, flux)) {
      assert(std::isnormal(Wleft[0]));
      assert(std::isnormal(Wright[0]));
      exactRiemann->ComputeFluxes(zeroMassFlux, rUnit, Wleft, Wright, vFace, flux);
    }
#elif defined(RIEMANN_SOLVER_HLLC)
    if (!hllcRiemann->ComputeFluxes(zeroMassFlux, rUnit, Wleft, Wright, vFace, flux)) {
      assert(std::isnormal(Wleft[0]));
      assert(std::isnormal(Wright[0]));
      exactRiemann->ComputeFluxes(zeroMassFlux, rUnit, Wleft, Wright, vFace, flux);
    }
#elif defined(RIEMANN_SOLVER_ROE)
    if (!roeRiemann->ComputeFluxes(zeroMassFlux, rUnit, Wleft, Wright, vFace, flux)) {
      assert(std::isnormal(Wleft[0]));
      assert(std::isnormal(Wright[0]));
      exactRiemann->ComputeFluxes(zeroMassFlux, rUnit, Wleft, Wright, vFace, flux);
    }
#else
    assert(std::isnormal(Wright[0]));
    exactRiemann->ComputeFluxes(zeroMassFlux, rUnit, Wleft, Wright, vFace, flux);
#endif

    // Rotate velocity component of flux back to orientation of original frame
    flux = rotMatrix*flux;

    //flux[irho] = 0; //explicetely set to zero for MFM

    // Add corrections to flux for boosting back to original lab frame
    Vector<ndim> fv;
    for (int k=0; k<ndim; k++) fv[k] = flux[k];
    flux[ipress] += (MyFloat)0.5*vFace.GetMagnitudeSquared()*flux[irho]; // this term is zero for MFM
    flux[ipress] += Vector<ndim>::DotProduct(vFace, fv);
    //for (int k=0; k<ndim; k++) flux[k] = flux[ipress]*rUnit[k];
    for (int k=0; k<ndim; k++) flux[k] += vFace[k]*flux[irho]; // zero for MFM

    return flux;
  }


#if defined(ONEDIM)
  template class FluxSolver<1>;
#elif defined(TWODIMS)
  template class FluxSolver<2>;
#else
  template class FluxSolver<3>;
#endif

}
//-------------------------------------------------------------------------------------------------
