#include <iostream>
#include "RiemannSolver.hpp"


//-------------------------------------------------------------------------------------------------
namespace OpenGadget3
{

  //===============================================================================================
  /// Compute the numerical fluxes based on the HLLC approximate Riemann solver.
  //===============================================================================================
  template <int ndim>
  bool HllcRiemannSolver<ndim>::ComputeFluxes
   (const bool zeroMassFlux,                     ///< [in] Compute fluxes in zero-mass flux frame
    const Vector<ndim> &rUnit,                   ///< [in] Unit vector on face
    WFluidVector<ndim> &Wleft,                   ///< [in] LHS primitive state
    WFluidVector<ndim> &Wright,                  ///< [in] RHS primitive state
    Vector<ndim> &vFace,                         ///< [in] Face velocity vector
    FluxVector<ndim> &flux) const                ///< [out] Flux vector
  {
    MyFloat Sleft;                               // Left-travelling wave speed of star region
    MyFloat Sright;                              // Right-travelling wave speed of star region
    //MyFloat pstar;
    
    // Compute left and right wave speeds, and the speed of the star region
    ComputeWaveSpeedEstimates(Wleft, Wright, Sleft, Sright);

    // If wave speed estimates are invalid (e.g. Sleft and Sright are converging), then signify failure of
    // approximate HLLC solver (by returning false) to allow the use of the ExactRiemannSolver instead
    if (Sleft >= Sright) {
#ifndef RIEMANN_SOLVER_HLL_IGNORE_OUTPUT
      std::cout << "SLEFT : " << Sleft << "    Sright : " << Sright << std::endl;
#endif
      return false;
    }

    // If selecting zero-mass flux fame (e.g. for MFM), then move to frame of contact discontinuity
    const MyFloat Sstar = ComputeStarVelocity(Wleft, Wright, Sleft, Sright/*, pstar*/);
    if (!std::isfinite(Sstar)) {
      std::cout << "SLEFT : " << Sleft << " Wleft: " << Wleft[ivx] << "    Sright : " << Sright << " Wright: " << Wright[ivx] << std::endl;
    }
    /*if (pstar < 0) {
      return false;
      }*/
    if((Wleft[ivx]==Sleft)||(Wright[ivx]==Sright))
      {
	flux[ivx] = (MyFloat)0.5 * (Wleft[ipress] + Wright[ipress]); // not correct yet...
	if(ndim > 1) flux[ivy] = (MyFloat)0.0;
	if(ndim > 2) flux[ivz] = (MyFloat)0.0;
	flux[irho] = (MyFloat)0.0;
	flux[ipress] = (MyFloat)0.0; // check again
	return false;
      }
    if (!std::isfinite(Sstar)) {
      return false;
    }
    
    if (zeroMassFlux) {
      Wleft[ivx] -= Sstar;
      Wright[ivx] -= Sstar;
      vFace += Sstar*rUnit;
    }

    // Compute flux depending on the different wave speeds in the simplified HLLC shock structure
    if (Sleft >= (MyFloat)0.0) {
      flux = this->ComputeFlux1d(Wleft);
    }
    else if ((Sleft <= (MyFloat)0.0)&&((MyFloat)0.0 <= Sstar)) {
      flux = ComputeHllcFlux1d(Wleft, Wright, Sleft, Sright, Sstar);
      //flux = ComputeHllcFlux1d(Wright, Wleft, Sright, Sleft, Sstar);
    }
    else if ((Sstar <= (MyFloat)0.0)&&((MyFloat)0.0 <= Sright)) {
      //flux = ComputeHllcFlux1d(Wleft, Wright, Sleft, Sright, Sstar);
      flux = ComputeHllcFlux1d(Wright, Wleft, Sright, Sleft, Sstar);
    }
    else if (Sright <= (MyFloat)0.0) {
      flux = this->ComputeFlux1d(Wright);
    }

    // Guarantee mass flux is zero (should be zero but in case of floating point rounding error)
    if (zeroMassFlux) {
      flux[irho] = (MyFloat)0.0;
    }

    return true;
  }


  //===============================================================================================
  /// Hllc Riemann solver, based on approach outlined by Toro (1999).
  /// Variable names often follow Toro naming conventions so read chapter for further details.
  //===============================================================================================
  template <int ndim>
  void HllcRiemannSolver<ndim>::ComputeWaveSpeedEstimates
   (const WFluidVector<ndim> &Wleft,             ///< [in] LHS fluid state vector
    const WFluidVector<ndim> &Wright,            ///< [in] RHS fluid state vector
    MyFloat &Sleft,                              ///< [out] Left travelling signal speed
    MyFloat &Sright) const                       ///< [out] Right travelling signal speed
  {
    Sleft  = Wleft[ivx] - eos->SoundSpeed(Wleft.GetDensity(), Wleft.GetSpecificInternalEnergy(gamma), 1); //physical
    Sright = Wright[ivx] + eos->SoundSpeed(Wright.GetDensity(), Wright.GetSpecificInternalEnergy(gamma), 1); //physical    
  }


  //===============================================================================================
  /// Hllc Riemann solver, based on approach outlined by Toro (1999).
  /// Variable names often follow Toro naming conventions so read chapter for further details.
  //===============================================================================================
  template <int ndim>
  MyFloat HllcRiemannSolver<ndim>::ComputeStarVelocity
   (const WFluidVector<ndim> &Wleft,             ///< [in] LHS fluid state vector
    const WFluidVector<ndim> &Wright,            ///< [in] RHS fluid state vector
    const MyFloat Sleft,                         ///< [in] Left travelling signal speed
    const MyFloat Sright/*,                        ///< [in] Right travelling signal speed
			  const MyFloat &pstar*/) const                  ///< [out] pressure in star region
  {
    const MyFloat jLeft = Wleft[irho]*(Sleft - Wleft[ivx]);
    const MyFloat jRight = Wright[irho]*(Sright - Wright[ivx]);
    MyFloat Sstar = (jRight*Wright[ivx] - jLeft*Wleft[ivx] + Wleft[ipress] - Wright[ipress]);
    Sstar /= (jRight - jLeft);

    //pstar = (MyFloat)0.5 * (Wleft[ipress] + Wright[ipress]) - (MyFloat)0.125 * (Wleft[ivx]+Wright[ivx]) * (Wleft[irho]+Wright[irho]) * (eos->SoundSpeed(Wleft.GetDensity(), Wleft.GetSpecificInternalEnergy(gamma), 1) + eos->SoundSpeed(Wright.GetDensity(), Wright.GetSpecificInternalEnergy(gamma), 1));
    
    return Sstar;
  }


  //===============================================================================================
  /// Hllc Riemann solver, based on approach outlined by Toro (1999).
  /// Variable names often follow Toro naming conventions so read chapter for further details.
  //===============================================================================================
  template <int ndim>
  UFluidVector<ndim> HllcRiemannSolver<ndim>::ComputeHllcStateVector
   (const WFluidVector<ndim> &Wleft,             ///< [in] LHS fluid state vector
    const WFluidVector<ndim> &Wright,            ///< [in] RHS fluid state vector
    const MyFloat Sleft,                         ///< [in] Left travelling signal speed
    const MyFloat Sright) const                  ///< [in] Right travelling signal speed
  {
    const FluxVector<ndim> fluxLeft = this->ComputeFlux1d(Wleft);
    const FluxVector<ndim> fluxRight = this->ComputeFlux1d(Wright);
    const UFluidVector<ndim> Uleft(Wleft, gamma);
    const UFluidVector<ndim> Uright(Wright, gamma);
    const MyFloat invSdiff = (MyFloat)1.0 / (Sright - Sleft);
    UFluidVector<ndim> Ustar;
    for (int k=0; k<nvar; k++) {
      Ustar[k] = Uright[k]*Sright - Uleft[k]*Sleft + fluxLeft[k] - fluxRight[k];
      Ustar[k] *= invSdiff;
    }
    return Ustar;
  }


  //===============================================================================================
  /// Hllc Riemann solver, based on approach outlined by Toro (1999).
  /// Variable names often follow Toro naming conventions so read chapter for further details.
  //===============================================================================================
  template <int ndim>
  FluxVector<ndim> HllcRiemannSolver<ndim>::ComputeHllcFlux1d
   (const WFluidVector<ndim> &Wleft,             ///< [in] LHS fluid state vector
    const WFluidVector<ndim> &Wright,            ///< [in] RHS fluid state vector
    const MyFloat Sleft,                         ///< [in] Left travelling signal speed
    const MyFloat Sright,                        ///< [in] Right travelling signal speed
    const MyFloat Sstar) const                   ///< [in] signal speed of contact disconstinuity
  {
    const FluxVector<ndim> fluxLeft = this->ComputeFlux1d(Wleft);
    //const FluxVector<ndim> fluxRight = this->ComputeFlux1d(Wright);
    const UFluidVector<ndim> Uleft(Wleft, gamma);
    //const UFluidVector<ndim> Uright(Wright, gamma);
    const MyFloat invSdiff = (MyFloat)1.0 / (Sleft - Sstar);
    FluxVector<ndim> flux;
#if HLLC_RIEMANN_VERSION == 3
    MyFloat PLR = (MyFloat)0.5 * (Wleft[ipress] + Wright[ipress] + Wleft[irho] * (Sleft - Wleft[ivx])*(Sstar - Wleft[ivx]) + Wright[irho]*(Sright-Wright[ivx])*(Sstar - Wright[ivx]));
#endif
    for (int k=0; k<nvar; k++) {
#if HLLC_RIEMANN_VERSION == 1 // Toro Eqn. (10.38) and (10.39)
      MyFloat fac;
      if (k == irho) {
	fac = 1;
      } else if (k == ivx) {
	fac = Sstar;
      } else if ((k == ivy) && (ndim > 1)) {
	fac = Wleft[ivy];
      } else if ((k == ivz) && (ndim > 2)) {
	fac = Wleft[ivz];
      } else if (k == ipress) {
	fac = Uleft[ipress] / Wleft[irho] + (Sstar - Wleft[ivx]) * (Sstar + Wleft[ipress]/(Wleft[irho]*(Sleft - Wleft[ivx])));
      } 
      MyFloat Ustarleftk = Wleft[irho] * (Sleft - Wleft[ivx]) * invSdiff * fac;
      flux[k] = fluxLeft[k] + Sleft * (Ustarleftk - Uleft[k]);
#elif HLLC_RIEMANN_VERSION == 2 // Toro Eqn. (10.41)
      MyFloat Dstar;
      if (k == irho) {
	Dstar = 0;
      } else if (k == ivx) {
        Dstar = 1;
      } else if ((k == ivy) && (ndim > 1)) {
        Dstar = 0;
      } else if ((k == ivz) && (ndim > 2)) {
        Dstar = 0;
      } else if (k == ipress) {
        Dstar = Sstar;
      }
      flux[k] = Sstar*(Sleft*Uleft[k] - fluxLeft[k]) + Sleft*(Wleft[ipress] + Wleft[irho]*(Sleft - Wleft[ivx])*(Sstar - Wleft[ivx])) * Dstar ;
      flux[k] *= invSdiff;
#else // HLLC_RIEMANN_VERSION == 3 // Toro Eqn. (10.44)
      MyFloat Dstar;
      if (k == irho) {
	Dstar = 0;
      } else if (k == ivx) {
        Dstar = 1;
      } else if ((k == ivy) && (ndim > 1)) {
        Dstar = 0;
      } else if ((k == ivz) && (ndim > 2)) {
        Dstar = 0;
      } else if (k == ipress) {
        Dstar = Sstar;
      }
      flux[k] = Sstar*(Sleft*Uleft[k] - fluxLeft[k]) + Sleft*PLR * Dstar;
      flux[k] *= invSdiff;
#endif
    }
    return flux;
  }


#if defined(ONEDIM)
  template class HllcRiemannSolver<1>;
#elif defined(TWODIMS)
  template class HllcRiemannSolver<2>;
#else
  template class HllcRiemannSolver<3>;
#endif


}
//-------------------------------------------------------------------------------------------------
