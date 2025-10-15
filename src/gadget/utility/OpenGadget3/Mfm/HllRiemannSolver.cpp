#include <iostream>
#include "RiemannSolver.hpp"

//-------------------------------------------------------------------------------------------------
namespace OpenGadget3
{

  //===============================================================================================
  /// Compute the numerical fluxes based on the HLL approximate Riemann solver.
  //===============================================================================================
  template <int ndim>
  bool HllRiemannSolver<ndim>::ComputeFluxes
   (const bool zeroMassFlux,                     ///< [in] Compute fluxes in zero-mass flux frame
    const Vector<ndim> &rUnit,                   ///< [in] Unit vector on face
    WFluidVector<ndim> &Wleft,                   ///< [in] LHS primitive state
    WFluidVector<ndim> &Wright,                  ///< [in] RHS primitive state
    Vector<ndim> &vFace,                         ///< [in] Face velocity vector
    FluxVector<ndim> &flux) const                ///< [out] Flux vector
  {
    MyFloat Sleft;                               // Left-travelling wave speed of star region
    MyFloat Sright;                              // Right-travelling wave speed of star region

    // Compute left and right wave speeds, and the speed of the star region
    ComputeWaveSpeedEstimates(Wleft, Wright, Sleft, Sright);

    // If wave speed estimates are invalid (e.g. Sleft and Sright are converging), then signify failure of
    // approximate HLL solver (by returning false) to allow the use of the ExactRiemannSolver instead
    if (Sleft >= Sright) {
#ifndef RIEMANN_SOLVER_HLL_IGNORE_OUTPUT
      std::cout << "SLEFT : " << Sleft << "    Sright : " << Sright << std::endl;
#endif
      return false;
    }

    // If selecting zero-mass flux fame (e.g. for MFM), then move to frame of contact discontinuity
    if (zeroMassFlux) {
      const MyFloat Sstar = ComputeStarVelocity(Wleft, Wright, Sleft, Sright);
      Wleft[ivx] -= Sstar;
      Wright[ivx] -= Sstar;
      vFace += Sstar*rUnit;
    }

    // Compute flux depending on the different wave speeds in the simplified HLL shock structure
    if (Sleft >= MyFloat(0.0)) {
      flux = this->ComputeFlux1d(Wleft);
    }
    else if (Sright <= MyFloat(0.0)) {
      flux = this->ComputeFlux1d(Wright);
    }
    else {
      flux = ComputeHllFlux1d(Wleft, Wright, Sleft, Sright);
    }

    // Guarantee mass flux is zero (should be zero but in case of floating point rounding error)
    if (zeroMassFlux) {
      flux[irho] = MyFloat(0.0);
    }

    return true;
  }


  //===============================================================================================
  /// Exact Riemann solver, based on approach outlined by Toro (1999).
  /// Variable names often follow Toro naming conventions so read chapter for further details.
  //===============================================================================================
  template <int ndim>
  void HllRiemannSolver<ndim>::ComputeWaveSpeedEstimates
   (const WFluidVector<ndim> &Wleft,             ///< [in] LHS fluid state vector
    const WFluidVector<ndim> &Wright,            ///< [in] RHS fluid state vector
    MyFloat &Sleft,                              ///< [out] Left travelling signal speed
    MyFloat &Sright) const                       ///< [out] Right travelling signal speed
  {
    Sleft  = Wleft[ivx] - eos->SoundSpeed(Wleft.GetDensity(), Wleft.GetSpecificInternalEnergy(gamma), 1); //physical
    Sright = Wright[ivx] + eos->SoundSpeed(Wright.GetDensity(), Wright.GetSpecificInternalEnergy(gamma), 1); //physical    
  }


  //===============================================================================================
  /// Exact Riemann solver, based on approach outlined by Toro (1999).
  /// Variable names often follow Toro naming conventions so read chapter for further details.
  //===============================================================================================
  template <int ndim>
  MyFloat HllRiemannSolver<ndim>::ComputeStarVelocity
   (const WFluidVector<ndim> &Wleft,             ///< [in] LHS fluid state vector
    const WFluidVector<ndim> &Wright,            ///< [in] RHS fluid state vector
    const MyFloat Sleft,                         ///< [out] Left travelling signal speed
    const MyFloat Sright) const                  ///< [out] Right travelling signal speed
  {
    const MyFloat jLeft = Wleft[irho]*(Sleft - Wleft[ivx]);
    const MyFloat jRight = Wright[irho]*(Sright - Wright[ivx]);
    MyFloat Sstar = (jRight*Wright[ivx] - jLeft*Wleft[ivx] + Wleft[ipress] - Wright[ipress]);
    Sstar /= (jRight - jLeft);
    return Sstar;
  }


  //===============================================================================================
  /// Exact Riemann solver, based on approach outlined by Toro (1999).
  /// Variable names often follow Toro naming conventions so read chapter for further details.
  //===============================================================================================
  template <int ndim>
  UFluidVector<ndim> HllRiemannSolver<ndim>::ComputeHllStateVector
   (const WFluidVector<ndim> &Wleft,             ///< [in] LHS fluid state vector
    const WFluidVector<ndim> &Wright,            ///< [in] RHS fluid state vector
    const MyFloat Sleft,                         ///< [out] Left travelling signal speed
    const MyFloat Sright) const                  ///< [out] Right travelling signal speed
  {
    const FluxVector<ndim> fluxLeft = this->ComputeFlux1d(Wleft);
    const FluxVector<ndim> fluxRight = this->ComputeFlux1d(Wright);
    const UFluidVector<ndim> Uleft(Wleft, gamma);
    const UFluidVector<ndim> Uright(Wright, gamma);
    const MyFloat invSdiff = MyFloat(1.0) / (Sright - Sleft);
    UFluidVector<ndim> Ustar;
    for (int k=0; k<nvar; k++) {
      Ustar[k] = Uright[k]*Sright - Uleft[k]*Sleft + fluxLeft[k] - fluxRight[k];
      Ustar[k] *= invSdiff;
    }
    return Ustar;
  }


  //===============================================================================================
  /// Exact Riemann solver, based on approach outlined by Toro (1999).
  /// Variable names often follow Toro naming conventions so read chapter for further details.
  //===============================================================================================
  template <int ndim>
  FluxVector<ndim> HllRiemannSolver<ndim>::ComputeHllFlux1d
   (const WFluidVector<ndim> &Wleft,             ///< [in] LHS fluid state vector
    const WFluidVector<ndim> &Wright,            ///< [in] RHS fluid state vector
    const MyFloat Sleft,                         ///< [out] Left travelling signal speed
    const MyFloat Sright) const                  ///< [out] Right travelling signal speed
  {
    const FluxVector<ndim> fluxLeft = this->ComputeFlux1d(Wleft);
    const FluxVector<ndim> fluxRight = this->ComputeFlux1d(Wright);
    const UFluidVector<ndim> Uleft(Wleft, gamma);
    const UFluidVector<ndim> Uright(Wright, gamma);
    const MyFloat invSdiff = MyFloat(1.0) / (Sright - Sleft);
    FluxVector<ndim> flux;
    for (int k=0; k<nvar; k++) {
      flux[k] = Sright*fluxLeft[k] - Sleft*fluxRight[k] + Sleft*Sright*(Uright[k] - Uleft[k]);
      flux[k] *= invSdiff;
    }
    return flux;
  }


#if defined(ONEDIM)
  template class HllRiemannSolver<1>;
#elif defined(TWODIMS)
  template class HllRiemannSolver<2>;
#else
  template class HllRiemannSolver<3>;
#endif


}
//-------------------------------------------------------------------------------------------------
