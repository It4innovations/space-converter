#include "RiemannSolver.hpp"
#include "Matrix.hpp"
//#include "../CodeBase/allvars.h"

//-------------------------------------------------------------------------------------------------
namespace OpenGadget3
{

  //===============================================================================================
  /// Compute the numerical fluxes based on the Roe approximate Riemann solver.
  //===============================================================================================
  template <int ndim>
  bool RoeRiemannSolver<ndim>::ComputeFluxes
   (const bool zeroMassFlux,                     ///< [in] Compute fluxes in zero-mass flux frame
    const Vector<ndim> &rUnit,                   ///< [in] Unit vector on face
    WFluidVector<ndim> &Wleft,                   ///< [in] LHS primitive state
    WFluidVector<ndim> &Wright,                  ///< [in] RHS primitive state
    Vector<ndim> &vFace,                         ///< [in] Face velocity vector
    FluxVector<ndim> &flux) const                ///< [out] Flux vector
  {
    Vector<nvar> valtilde;                       // averaged values: soundspeed, velocity, entalpy
    Vector<nvar> lambda;                         // eigenvalues
    SquareMatrix<nvar> Ktilde;                   // eigenvectors
    Vector<nvar> alpha;                          // wave strength
    
    /*if (zeroMassFlux) {
      Wleft[ivx] -= Sstar;
      Wright[ivx] -= Sstar;
      vFace += Sstar*rUnit;
      }*/

    // Compute averaged values
    ComputeAveragedValues(Wleft, Wright, valtilde);
    
    ComputeEigenvalues(valtilde, lambda);
    ComputeRightEigenvectors(valtilde, Ktilde);
    ComputeWaveStrength(Wleft, Wright, valtilde, alpha);
    
    // Compute flux depending on the different wave speeds in the simplified shock structure
    if (lambda[irho] >= (MyFloat)0.0) {
      flux = this->ComputeFlux1d(Wleft);
    }
    //else if ((valtilde[irho] <= (MyFloat)0.0)&&((MyFloat)0.0 <= valtilde[ivx])) {
    //  flux = ComputeRoeFlux1d(Wleft, Wright, alpha, lambda, Ktilde);    
    //} else if ((valtilde[ivx] <= (MyFloat)0.0)&&((MyFloat)0.0 <= valtilde[ipress])) {
    //  flux = ComputeRoeFlux1d(Wright, Wleft, alpha, lambda, Ktilde);    
    else if ((lambda[irho] < (MyFloat)0.0)&&((MyFloat)0.0 < lambda[ipress])) {
      flux = ComputeRoeFlux1d(Wleft, Wright, alpha, lambda, Ktilde);
    } else if (lambda[ipress] <= (MyFloat)0.0) {
      flux = this->ComputeFlux1d(Wright);
    }
    
    // Guarantee mass flux is zero (should be zero but in case of floating point rounding error)
    if (zeroMassFlux) {
      flux[irho] = (MyFloat)0.0;
    }
    
    return true;
  }


  //===============================================================================================
  /// Roe Riemann solver, based on approach outlined by Toro (1999).
  /// Variable names often follow Toro naming conventions so read chapter for further details.
  //===============================================================================================
  template <int ndim>
  void RoeRiemannSolver<ndim>::ComputeAveragedValues
   (const WFluidVector<ndim> &Wleft,             ///< [in] LHS fluid state vector
    const WFluidVector<ndim> &Wright,            ///< [in] RHS fluid state vector
    Vector<nvar> &valtilde) const                ///< [out] average values
  {
    for (int k = 0; k < nvar; k++)
      valtilde[k] = (MyFloat)0.0;
    MyFloat sqrtrhol = sqrt(Wleft[irho]);
    MyFloat sqrtrhor = sqrt(Wright[irho]);
    MyFloat invsqrtrhosum = 1.0 / (sqrtrhol + sqrtrhor);
    valtilde[ivx] = (sqrtrhol*Wleft[ivx] + sqrtrhor*Wright[ivx]) * invsqrtrhosum;
    if(ndim > 1) valtilde[ivy] = (sqrtrhol*Wleft[ivy] + sqrtrhor*Wright[ivy]) * invsqrtrhosum;
    if(ndim > 2) valtilde[ivz] = (sqrtrhol*Wleft[ivz] + sqrtrhor*Wright[ivz]) * invsqrtrhosum;
    MyFloat Hleft = Wleft.GetSpecificInternalEnergy(gamma) + Wleft[ipress]/Wleft[irho];
    MyFloat Hright = Wright.GetSpecificInternalEnergy(gamma) + Wright[ipress]/Wright[irho];
    for(int i = 0; i < ndim; i++)
      {
	Hleft += (MyFloat)0.5 * Wleft[i]*Wleft[i];
	Hright += (MyFloat)0.5 * Wright[i]*Wright[i];
      }
    valtilde[ipress] = (sqrtrhol*Hleft + sqrtrhor*Hright) * invsqrtrhosum; // Htilde
    valtilde[irho] = (gamma-(MyFloat)1.0) * valtilde[ipress];
    for(int i = 0; i < ndim; i++)
      valtilde[irho] -= (MyFloat)0.5 * (gamma-(MyFloat)1.0) * valtilde[i]*valtilde[i]; // atilde
    valtilde[irho] = sqrt(valtilde[irho]);
  }


  //===============================================================================================
  /// Roe Riemann solver, based on approach outlined by Toro (1999).
  /// Variable names often follow Toro naming conventions so read chapter for further details.
  //===============================================================================================
  template <int ndim>
  void RoeRiemannSolver<ndim>::ComputeEigenvalues
   (const Vector<nvar> &val,                     ///< [in] values
    Vector<nvar> &lambda) const                  ///< [out] eigenvalues
  {
    for (int k = 0; k < nvar; k++)
      lambda[k] = (MyFloat)0.0;
    lambda[ivx] = val[ivx];
    if(ndim > 1) lambda[ivy] = val[ivx];
    if(ndim > 2) lambda[ivz] = val[ivx];
    lambda[irho] = val[ivx] - val[irho];
    lambda[ipress] = val[ivx] + val[irho];
  }


  //===============================================================================================
  /// Roe Riemann solver, based on approach outlined by Toro (1999).
  /// Variable names often follow Toro naming conventions so read chapter for further details.
  //===============================================================================================
  template <int ndim>
  void RoeRiemannSolver<ndim>::ComputeRightEigenvectors
   (const Vector<nvar> &val,                     ///< [in] values
    SquareMatrix<nvar> &K) const                 ///< [out] eigenvectors
  {
    for(int i = 0; i < nvar; i++)
      for(int j = 0; j < nvar; j++)
	K(i,j) = (MyFloat)0.0;
    
    K(ivx,ivx) = val[ivx];
    if(ndim > 1) K(ivy,ivx) = val[ivy];
    if(ndim > 2) K(ivz,ivx) = val[ivz];
    K(irho,ivx) = (MyFloat)1.0;
    for(int i = 0; i < ndim; i++)
      K(ipress,ivx) += (MyFloat)0.5 * val[i]*val[i];

    if(ndim > 1)
      {
	K(ivy,ivy) = (MyFloat)1.0;
	K(ipress,ivy) = val[ivy];
      }
    if(ndim > 2)
      {
	K(ivz,ivz) = (MyFloat)1.0;
	K(ipress,ivz) = val[ivz];
      }
    
    K(ivx,irho) = val[ivx] - val[irho];
    if(ndim > 1) K(ivy,irho) = val[ivy];
    if(ndim > 2) K(ivz,irho) = val[ivz];
    K(irho,irho) = (MyFloat)1.0;
    K(ipress,irho) = val[ipress] - val[ivx]*val[irho];

    K(ivx,ipress) = val[ivx] + val[irho];
    if(ndim > 1) K(ivy,ipress) = val[ivy];
    if(ndim > 2) K(ivz,ipress) = val[ivz];
    K(irho,ipress) = (MyFloat)1.0;
    K(ipress,ipress) = val[ipress] + val[ivx]*val[irho];
  }


  //===============================================================================================
  /// Roe Riemann solver, based on approach outlined by Toro (1999).
  /// Variable names often follow Toro naming conventions so read chapter for further details.
  //===============================================================================================
  template <int ndim>
  void RoeRiemannSolver<ndim>::ComputeWaveStrength
   (const WFluidVector<ndim> &Wleft,             ///< [in] LHS fluid state vector
    const WFluidVector<ndim> &Wright,            ///< [in] RHS fluid state vector
    const Vector<nvar> val,                 ///< [in] average velocity
    Vector<nvar> &alpha) const                   ///< [out] wave strength values
  {
    const UFluidVector<ndim> Uleft(Wleft, gamma);
    const UFluidVector<ndim> Uright(Wright, gamma);
    Vector<nvar> deltau;
    for(int i = 0; i< nvar; i++)
      deltau[i] = (MyFloat)0.0;
    for(int i =0; i< nvar; i++)
      deltau[i] = Uright[i] - Uleft[i];
    MyFloat bardeltau5 = deltau[ipress]; //Uright[ipress] - Uleft[ipress]; // other terms are zero
    for(int k = 0; k < nvar; k++)
      alpha[k] = (MyFloat)0.0;
    alpha[ivx] = (gamma-(MyFloat)1.0) / (val[irho]*val[irho]) * (deltau[irho]*(val[ipress]-val[ivx]*val[ivx]) + val[ivx]*deltau[ivx] - bardeltau5);
    if(ndim > 1) alpha[ivy] = deltau[ivy] - deltau[irho]*val[ivy];
    if(ndim > 2) alpha[ivz] = deltau[ivz] - deltau[irho]*val[ivz];
    alpha[irho] = (MyFloat)0.5/val[irho] * (deltau[irho]*(val[ivx] + val[irho]) - deltau[ivx] - val[irho]*alpha[ivx]);
    alpha[ipress] = deltau[irho] - (alpha[irho] + alpha[ivx]);
  }

  
  //===============================================================================================
  /// Roe Riemann solver, based on approach outlined by Toro (1999).
  /// Variable names often follow Toro naming conventions so read chapter for further details.
  //===============================================================================================
  template <int ndim>
  FluxVector<ndim> RoeRiemannSolver<ndim>::ComputeRoeFlux1d
   (const WFluidVector<ndim> &Wleft,             ///< [in] LHS fluid state vector
    const WFluidVector<ndim> &Wright,            ///< [in] RHS fluid state vector
    const Vector<nvar> &alpha,                   ///< [in] wave strength
    const Vector<nvar> &lambda,                  ///< [in] eigenvalues
    const SquareMatrix<nvar> &Ktilde) const      ///< [in] eigenvectors
  {
    MyFloat uStarLeft, uStarRight, aStarLeft, aStarRight;
    ComputeStarVelocities(Wleft, Wright, alpha, lambda, Ktilde, uStarLeft, uStarRight, aStarLeft, aStarRight);
    
#ifndef ROE_RIGHT_FLUXES
    // left
    const FluxVector<ndim> fluxLeft = this->ComputeFlux1d(Wleft);
    MyFloat lambda1Left = Wleft[ivx] - eos->SoundSpeed(Wleft.GetDensity(), Wleft.GetSpecificInternalEnergy(gamma), 1);
    MyFloat lambda1Right = uStarLeft - aStarLeft;
    bool LeftIsTranssonic = ((lambda1Left < 0.0) && (0 < lambda1Right));
    MyFloat barlambda1;
    if(LeftIsTranssonic)
      barlambda1 = lambda1Left * (lambda1Right - lambda[irho]) / (lambda1Right - lambda1Left);
#endif
#ifndef ROE_LEFT_FLUXES
    // right
    const FluxVector<ndim> fluxRight = this->ComputeFlux1d(Wright);
    MyFloat lambda5Left = uStarRight + aStarRight;
    MyFloat lambda5Right = Wright[ivx] + eos->SoundSpeed(Wright.GetDensity(), Wright.GetSpecificInternalEnergy(gamma), 1);
    bool RightIsTranssonic = ((lambda5Left < 0.0) && (0 < lambda5Right));
    MyFloat barlambda5;
    if(RightIsTranssonic)
      barlambda5 = lambda5Right * (lambda[ipress] - lambda5Left) / (lambda5Right - lambda5Left);
#endif
    
    FluxVector<ndim> flux;
    for(int k = 0; k < nvar; k++)
      flux[k] = (MyFloat)0.0;

#ifndef ROE_RIGHT_FLUXES
    flux += fluxLeft;
#endif
#ifndef ROE_LEFT_FLUXES
    flux += fluxRight;
#endif
#ifndef ROE_RIGHT_FLUXES
    // left-sided flux
    if(LeftIsTranssonic) {
      for(int kk = 0; kk < nvar; kk++)
	flux[kk] += alpha[irho]*barlambda1*Ktilde(kk,irho);
    } else {
      for(int k = 0; k < nvar; k++)
	if(lambda[k] <= (MyFloat)0.0)
	  for(int kk = 0; kk < nvar; kk++)
	    flux[kk] += alpha[k]*lambda[k]*Ktilde(kk,k);
    }
#endif
#ifndef ROE_LEFT_FLUXES
    // right-sided flux
    if(RightIsTranssonic) {
      for(int kk = 0; kk < nvar; kk++)
	flux[kk] += alpha[ipress]*barlambda5*Ktilde(kk,ipress);
    } else {
      for(int k = 0; k < nvar; k++)
	if(lambda[k] >= (MyFloat)0.0)
	  for(int kk = 0; kk < nvar; kk++)
	    flux[kk] -= alpha[k]*lambda[k]*Ktilde(kk,k);
    }
#endif

#if(!defined(ROE_LEFT_FLUXES) && !defined(ROE_LEFT_FLUXES))
    for(int k = 0; k < nvar; k++)
      flux[k] = (MyFloat)0.5 * flux[k];
#endif
    return flux;
  }

  //===============================================================================================
  /// Roe Riemann solver, based on approach outlined by Toro (1999).
  /// Variable names often follow Toro naming conventions so read chapter for further details.
  //===============================================================================================
  template <int ndim>
  void RoeRiemannSolver<ndim>::ComputeStarVelocities
   (const WFluidVector<ndim> &Wleft,             ///< [in] LHS fluid state vector
    const WFluidVector<ndim> &Wright,            ///< [in] RHS fluid state vector
    const Vector<nvar> &alpha,                   ///< [in] wave strength
    const Vector<nvar> &lambda,                  ///< [in] eigenvalues
    const SquareMatrix<nvar> &Ktilde,            ///< [in] eigenvectors
    MyFloat &uStarLeft,                          ///< [out] primitive fluid variables in star region, derived from left-sided values
    MyFloat &uStarRight,                         ///< [out] primitive fluid variables in star region, derived from right-sided values
    MyFloat &aStarLeft,                          ///< [out] left sound speed in star region
    MyFloat &aStarRight) const                   ///< [out] right sound speed in star region
  {
#ifdef ROE_AVERAGE
    uStarLeft = (Wleft[irho]*Wleft[ivx] + alpha[irho]*lambda[irho]) / (Wleft[irho] + alpha[irho]);
    MyFloat rhoStar = Wleft[irho] + alpha[irho];
    MyFloat pStar = (gamma - (MyFloat)1.0) * (Wleft[irho]*Wleft.GetSpecificTotalEnergy(gamma) + alpha[irho]*Ktilde(ipress,irho) - (MyFloat)0.5*rhoStar*uStarLeft*uStarLeft );
    aStarLeft = sqrt(gamma*pStar/rhoStar); // replace with sound speed using eos???

    uStarRight = (Wright[irho]*Wright[ivx] - alpha[ipress]*lambda[ipress]) / (Wright[irho] - alpha[ipress]);
    rhoStar = Wright[irho] - alpha[ipress];
    pStar = (gamma - (MyFloat)1.0) * (Wright[irho]*Wright.GetSpecificTotalEnergy(gamma) + alpha[ipress]*Ktilde(ipress,ipress) - (MyFloat)0.5*rhoStar*uStarRight*uStarRight );
    aStarRight = sqrt(gamma*pStar/rhoStar); // replace with sound speed using eos???
#elif defined(PVRS_APPROXIMATION)
    MyFloat barrho = (MyFloat)0.5 * (Wleft[irho] + Wright[irho]);
    MyFloat bara = (MyFloat)0.5 * (eos->SoundSpeed(Wleft.GetDensity(), Wleft.GetSpecificInternalEnergy(gamma), 1) + eos->SoundSpeed(Wright.GetDensity(), Wright.GetSpecificInternalEnergy(gamma), 1));
    MyFloat pStar = (MyFloat)0.5 * ((Wleft[ipress]+Wright[iopress]) + (Wleft[ivx] - Wright[ivx])*barrho*bara);
    pstar = max(0, pStar);
    uStarLeft = (MyFloat)0.5 * ((Wleft[ivx] + Wright[ivx]) + (Wleft[ipress] - Wright[ipress])/(rhobar*abar));
    uStarRight = uStarLeft;
    MyFloat rhoStarLeft = Wleft[irho] + (Wleft[ivx] - uStar)*barrho/bara;
    MyFloat rhoStarRight = Wright[irho] + (uStar - Wright[ivx])*barrho/bara;
    aStarLeft = sqrt(gamma*pStar/rhoStarLeft); // replace with sound speed using eos???
    aStarRight = sqrt(gamma*pStar/rhoStarRight); // replace with sound speed using eos???
#elif defined(TRRS_APPROXIMATION)
    MyFloat z = (gamma-(MyFloat)1.0) / ((MyFloat)2.0 * gamma);
    MyFloat pStar = pow(eos->SoundSpeed(Wleft.GetDensity(), Wleft.GetSpecificInternalEnergy(gamma), 1) + eos->SoundSpeed(Wright.GetDensity(), Wright.GetSpecificInternalEnergy(gamma), 1) - (MyFloat)0.5*(gamma-(MyFloat)1.0)*(Wright[ivx] - Wleft[ivx]) / pow(aLeft/Wleft[ipress],z) + pow(aRight/Wright[ipress],z), 1/z);
    aStarLeft = eos->SoundSpeed(Wleft.GetDensity(), Wleft.GetSpecificInternalEnergy(gamma), 1) * pow(pStar/Wleft[ipress], z);
    uStarLeft = Wleft[ipress] + (MyFloat)2.0/(gamma - (MyFloat)1.0)*(eos->SoundSpeed(Wleft.GetDensity(), Wleft.GetSpecificInternalEnergy(gamma), 1) - aStarLeft);
    aStarRight = eos->SoundSpeed(Wright.GetDensity(), Wright.GetSpecificInternalEnergy(gamma), 1) * pow(pStar/Wright[ipress], z);
    uStarRight = Wright[ipress] + (MyFloat)2.0/(gamma - (MyFloat)1.0)*(aStarRight - eos->SoundSpeed(Wright.GetDensity(), Wright.GetSpecificInternalEnergy(gamma), 1));    
#endif
  }



#if defined(ONEDIM)
  template class RoeRiemannSolver<1>;
#elif defined(TWODIMS)
  template class RoeRiemannSolver<2>;
#else
  template class RoeRiemannSolver<3>;
#endif


}
//-------------------------------------------------------------------------------------------------
