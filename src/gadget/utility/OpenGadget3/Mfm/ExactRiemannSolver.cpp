#include "RiemannSolver.hpp"

//-------------------------------------------------------------------------------------------------
namespace OpenGadget3
{

  //===============================================================================================
  /// Exact Riemann solver, based on approach outlined by Toro (1999).
  //===============================================================================================
  template <int ndim>
  bool ExactRiemannSolver<ndim>::ComputeFluxes
   (const bool zeroMassFlux,                     ///< [in] Compute fluxes in zero-mass flux frame
    const Vector<ndim> &rUnit,                   ///< [in] Unit vector on face
    WFluidVector<ndim> &Wleft,                   ///< [in] LHS primitive state
    WFluidVector<ndim> &Wright,                  ///< [in] RHS primitive state
    Vector<ndim> &vFace,                         ///< [inout] Face velocity vector
    FluxVector<ndim> &flux) const                ///< [out] Flux vector
  {
    const MyFloat cl = eos->SoundSpeed(Wleft.GetDensity(), Wleft.GetSpecificInternalEnergy(gamma), 1); //physical
    const MyFloat cr = eos->SoundSpeed(Wright.GetDensity(), Wright.GetSpecificInternalEnergy(gamma), 1); //physical
    MyFloat uStar;                               // Velocity in star region
    MyFloat pStar;                               // Pressure in star region
    MyFloat u;                                   // Riemann solver state velocity
    MyFloat d;                                   // Riemann solver state density
    MyFloat p;                                   // Riemann solver state pressure
    WFluidVector<ndim> Wface;                    // Riemann solver state vector

    // Compute p and u values at interface (in star region)
    ComputeStarRegion(Wleft[irho], Wright[irho], Wleft[ivx], Wright[ivx],
                      Wleft[ipress], Wright[ipress], cl, cr, uStar, pStar);

    // If not a vacuum state (i.e. pStar > 0), then sample solution
    //---------------------------------------------------------------------------------------------
    if (pStar > MyFloat(0.0)) {
      SampleExactSolution(MyFloat(0.0), pStar, uStar, Wleft[ipress], Wright[ipress], Wleft[irho],
                          Wright[irho], cl, cr, Wleft[ivx], Wright[ivx], p, d, u);

      // Copy Riemann solver state
      Wface[ivx]    = u;
      Wface[irho]   = d;
      Wface[ipress] = p;

      // If selecting face that moves with the star region velocity (i.e. so mass flux is zero),
      // then re-set face velocity to match star region (in x-direction)
      if (zeroMassFlux) {
        Wface[ivx] = MyFloat(0.0);
        for (int k=0; k<ndim; k++) vFace[k] += u*rUnit[k];
      }
      // Otherwise, calculate transverse velocity fluxes from upwind direction
      else {
        if (u > MyFloat(0.0)) {
          for (int kv=1; kv<ndim; kv++) Wface[kv] = Wleft[kv];
        }
        else {
          for (int kv=1; kv<ndim; kv++) Wface[kv] = Wright[kv];
        }
      }

      flux = this->ComputeFlux1d(Wface);
      Wface[ivx] += vFace[ivx];

      return true;
    }
    // Otherwise assume vacuum state conditions
    //---------------------------------------------------------------------------------------------
    else {
      for (int k=0; k<nvar; k++) flux[k] = MyFloat(0.0);
      //exit(0);
      return false;
    }
    //---------------------------------------------------------------------------------------------

  }


  //===============================================================================================
  /// Exact Riemann solver, based on approach outlined by Toro (1999).
  /// Variable names often follow Toro naming conventions so read chapter for further details.
  //===============================================================================================
  template <int ndim>
  void ExactRiemannSolver<ndim>::ComputeStarRegion
   (const MyFloat dl,                            ///< [in] LHS density
    const MyFloat dr,                            ///< [in] RHS density
    const MyFloat ul,                            ///< [in] LHS velocity
    const MyFloat ur,                            ///< [in] RHS velocity
    const MyFloat pl,                            ///< [in] LHS pressure
    const MyFloat pr,                            ///< [in] RHS pressure
    const MyFloat cl,                            ///< [in] LHS sound speed
    const MyFloat cr,                            ///< [in] RHS sound speed
    MyFloat &uStar,                              ///< [out] Velocity of intermediate state
    MyFloat &pStar) const                        ///< [out] Intermediate pressure state
  {
    const MyFloat Al = g5/dl;                                // ..
    const MyFloat Ar = g5/dr;                                // ..
    const MyFloat Bl = pl*g6;                                // ..
    const MyFloat Br = pr*g6;                                // ..
    const MyFloat quser = (MyFloat) 2.0;                     // First guess of solution
    const MyFloat cup = (MyFloat) 0.25*(dl + dr)*(cl + cr);  // ..
    const MyFloat pmin = std::min(pl, pr);                   // ..
    const MyFloat pmax = std::max(pl, pr);                   // ..
    int iteration = 0;                                       // Iteration counter
    MyFloat fl;                                              // Left-state variable
    MyFloat flprime;                                         // Left-state variable gradient
    MyFloat fr;                                              // Right-state variable
    MyFloat frprime;                                         // Right-state variable gradient
    MyFloat pold;                                            // Old intermediate pressure
    MyFloat ppv;                                             // Primitive Variable pressure

    // Check all input variables are valid numbers
    assert(std::isnormal(dl));
    assert(std::isnormal(dr));
    assert(std::isnormal(pl));
    assert(std::isnormal(pr));
    assert(std::isnormal(cl));
    assert(std::isnormal(cr));

    // Check for vacuum condition
    if (cl + cr - g7*(ur - ul) <= MyFloat(0.0)) {
      pStar = (MyFloat)1.0e-30;  //0.0;
      uStar = MyFloat(0.0);
      return;
    }

    // Compute guess for pStar from Primitive Variable Riemann Solver
    ppv = MyFloat(0.5)*(pl + pr) + MyFloat(0.5)*(ul - ur)*cup;
    ppv = std::max(MyFloat(0.0), ppv);

    // Select PVRS Riemann Solver
    if (pmax/pmin <= quser && pmin <= ppv && ppv <= pmax) {
      pStar = ppv;
    }
    // Select two-rarefaction Riemann Solver
    else if (ppv < pmin) {
      const MyFloat pq  = pow(pl/pr, g1);
      const MyFloat um  = (pq*ul/cl + ur/cr + g4*(pq - MyFloat(1.0)))/(pq/cl + MyFloat(1.0)/cr);
      const MyFloat ptl = MyFloat(1.0) + g7*(ul - um)/cl;
      const MyFloat ptr = MyFloat(1.0) + g7*(um - ur)/cr;
      pStar = MyFloat(0.5)*(pl*pow(ptl, g3) + pr*pow(ptr, g3));
    }
    // Select two-shock Riemann Solver with PVRS as estimate
    else {
      const MyFloat gel = sqrt((g5/dl)/(g6*pl + ppv));
      const MyFloat ger = sqrt((g5/dr)/(g6*pr + ppv));
      pStar = (gel*pl + ger*pr - (ur - ul))/(gel + ger);
    }

    pStar = std::max(pStar, (MyFloat)1.0e-30);
    assert(std::isnormal(pStar));


    // Main iteration loop
    //---------------------------------------------------------------------------------------------
    do {
      iteration++;

      // Calculate contribution to f and fprime for LHS
      if (pStar > pl) {
        fl = (pStar - pl)*sqrt(Al/(pStar + Bl));
        flprime = sqrt(Al/(pStar + Bl))*(MyFloat(1.0) - MyFloat(0.5)*(pStar - pl)/(pStar + Bl));
      }
      else {
        fl = g4*cl*(pow(pStar/pl, g1) - MyFloat(1.0));
        flprime = pow(pStar/pl, -g2)/(dl*cl);
      }

      // Calculate contribution to f and fprime for RHS
      if (pStar > pr) {
        fr = (pStar - pr)*sqrt(Ar/(pStar + Br));
        frprime = sqrt(Ar/(pStar + Br))*(MyFloat(1.0) - MyFloat(0.5)*(pStar - pr)/(pStar + Br));
      }
      else {
        fr = g4*cr*(pow(pStar/pr, g1) - MyFloat(1.0));
        frprime = pow(pStar/pr, -g2)/(dr*cr);
      }

      // Perform Newton-Raphson iteration
      pold = pStar;
      pStar = pStar - (fl + fr + ur - ul)/(flprime + frprime);

      // Check if convergence has been achieved
      if (pStar < MyFloat(1.0e-30)) pStar = MyFloat(1.0e-30); //this could also be smaller than zero
      else if ((MyFloat) 2.0*fabs(pStar - pold)/(pStar + pold) < tolerance) break;

      assert(std::isnormal(pStar));

    } while (iteration < maxNumIterations);
    //---------------------------------------------------------------------------------------------

    // Compute velocity of star region
    uStar = MyFloat(0.5)*(ul + ur) + MyFloat(0.5)*(fr - fl);
    assert(std::isfinite(uStar));

    return;
  }


  //===============================================================================================
  /// Samples the solution of the Exact Riemann solver at the given
  //===============================================================================================
  template <int ndim>
  void ExactRiemannSolver<ndim>::SampleExactSolution
   (const MyFloat s,                               ///< [in] (Dimensionless) position in shock
    const MyFloat pStar,                           ///< [in] Pressure of star region
    const MyFloat uStar,                           ///< [in] Velocity of star region
    const MyFloat pl,                              ///< [in] LHS pressure
    const MyFloat pr,                              ///< [in] RHS pressure
    const MyFloat dl,                              ///< [in] LHS density
    const MyFloat dr,                              ///< [in] RHS density
    const MyFloat cl,                              ///< [in] LHS sound speed
    const MyFloat cr,                              ///< [in] RHS sound speed
    const MyFloat ul,                              ///< [in] LHS velocity
    const MyFloat ur,                              ///< [in] RHS velocity
    MyFloat &p,                                    ///< [out] Pressure at s
    MyFloat &d,                                    ///< [out] Density at s
    MyFloat &u) const                              ///< [out] Velocity at s
  {

    // If sampling point lies to the left of the contact discontinuity
    //---------------------------------------------------------------------------------------------
    if (s <= uStar) {

      // Left rarefaction
      //-------------------------------------------------------------------------------------------
      if (pStar <= pl) {
        const MyFloat sh = ul - cl;

        // Sampled point is left data state
        if (s <= sh) {
          d = dl;
          u = ul;
          p = pl;
        }
        else {
          const MyFloat cm = cl*pow(pStar/pl, g1);
          const MyFloat st = uStar - cm;

          // Sample point is star left state
          if (s > st) {
            d = dl*pow(pStar/pl, invgamma);
            u = uStar;
            p = pStar;
          }
          // Sampled point is inside left fan
          else {
            const MyFloat c = g5*(cl + g7*(ul - s));
            u = g5*(cl + g7*ul + s);
            d = dl*pow(c/cl, g4);
            p = pl*pow(c/cl, g3);
          }
        }
      }

      // Left shock
      //-------------------------------------------------------------------------------------------
      else {
        const MyFloat pml = pStar/pl;
        const MyFloat sl = ul - cl*sqrt(g2*pml + g1);

        // Sampled point is left data state
        if (s <= sl) {
          d = dl;
          u = ul;
          p = pl;
        }
        // Sampled point is star left state
        else {
          d = dl*(pml + g6)/(pml*g6 + MyFloat(1.0));
          u = uStar;
          p = pStar;
        }
      }
      //-------------------------------------------------------------------------------------------

    }
    // Sampling point lies to the right of the contact discontinuity
    //---------------------------------------------------------------------------------------------
    else {

      // Right shock
      //-------------------------------------------------------------------------------------------
      if (pStar >= pr) {
        const MyFloat pmr = pStar/pr;
        const MyFloat sr = ur + cr*sqrt(g2*pmr + g1);

        // Sampled point is right data state
        if (s >= sr) {
          d = dr;
          u = ur;
          p = pr;
        }
        // Sampled point is star right state
        else {
          d = dr*(pmr + g6)/(pmr*g6 + MyFloat(1.0));
          u = uStar;
          p = pStar;
        }

      }
      // Right rarefaction
      //-------------------------------------------------------------------------------------------
      else {
        const MyFloat sh = ur + cr;

        // Sampled point is left data state
        if (s >= sh) {
          d = dr;
          u = ur;
          p = pr;
        }
        else {
          const MyFloat cm = cr*pow(pStar/pr,g1);
          const MyFloat st = uStar + cm;

          // Sample point is star right state
          if (s <= st) {
            d = dr*pow(pStar/pr, invgamma);
            u = uStar;
            p = pStar;
          }
          // Sampled point is inside left fan
          else {
            const MyFloat c = g5*(cr - g7*(ur - s));
            u = g5*(-cr + g7*ur + s);
            d = dr*pow(c/cr,g4);
            p = pr*pow(c/cr,g3);
          }
        }
      }

    }
    //---------------------------------------------------------------------------------------------

    assert(std::isnormal(p));
    assert(std::isfinite(d));
    assert(std::isfinite(u));

    return;
  }


#if defined(ONEDIM)
  template class ExactRiemannSolver<1>;
#elif defined(TWODIMS)
  template class ExactRiemannSolver<2>;
#else
  template class ExactRiemannSolver<3>;
#endif

}
//-------------------------------------------------------------------------------------------------
