#include "FluidVector.hpp"


//-------------------------------------------------------------------------------------------------
namespace OpenGadget3
{

  //===============================================================================================
  /// Constructor to convert primitive fluid vector (W) to conserved fluid vector (Q).
  //===============================================================================================
  template <int ndim>
  QFluidVector<ndim>::QFluidVector
   (const WFluidVector<ndim> &W,                 ///< [in] Primtive fluid vector to convert to Q
    const MyFloat vol,                           ///< [in] Volume of particle
    const MyFloat gamma,                         ///< [in] Ratio of specific heats
    const MyFloat a)                             ///< [in] Scale factor
  {
    assert(std::isnormal(vol));
    assert(std::isnormal(gamma));
    assert(std::isnormal(W[irho]));
    assert(std::isnormal(W[ipress]));
    const MyFloat mass = W[irho]*vol;
    const MyFloat gammam1 = gamma - (MyFloat)1.0;
    MyFloat vSqd = (MyFloat)0.0;
    for (int k=0; k<ndim; k++) vSqd += (W[k]*W[k]);
    for (int k=0; k<ndim; k++) data[k] = W[k]*mass;
    data[irho]   = mass;
    data[ipress] = W[ipress]*vol/gammam1/pow(a,3*gammam1) + (MyFloat)0.5*mass*vSqd;
  }


  //===============================================================================================
  /// Constructor to convert primitive fluid vector (W) to conserved fluid vector (Q).
  //===============================================================================================
#if NUMDIMS == 1
  template <>
  MyFloat QFluidVector<1>::GetSpecificInternalEnergy() const
  {
    assert(std::isnormal(data[irho]));
    const MyFloat invm = (MyFloat) 1.0 / data[irho];
    return data[ipress]*invm - (MyFloat) 0.5*invm*invm*data[ivx]*data[ivx];
  }
#endif
#if NUMDIMS == 2
  template <>
  MyFloat QFluidVector<2>::GetSpecificInternalEnergy() const
  {
    assert(std::isnormal(data[irho]));
    const MyFloat invm = (MyFloat) 1.0 / data[irho];
    return data[ipress]*invm - (MyFloat) 0.5*invm*invm*(data[ivx]*data[ivx] + data[ivy]*data[ivy]);
  }
#endif
#if NUMDIMS == 3
  template <>
  MyFloat QFluidVector<3>::GetSpecificInternalEnergy() const
  {
    assert(std::isnormal(data[irho]));
    const MyFloat invm = (MyFloat) 1.0 / data[irho];
    return data[ipress]*invm - (MyFloat) 0.5*invm*invm*(data[ivx]*data[ivx] + data[ivy]*data[ivy] + data[ivz]*data[ivz]);
  }
#endif


  //===============================================================================================
  /// Constructor to convert primitive fluid vector (W) to conserved fluid vector (U).
  //===============================================================================================
  template <int ndim>
  UFluidVector<ndim>::UFluidVector
   (const WFluidVector<ndim> &W,                 ///< [in] Primtive fluid vector to convert to U
    const MyFloat gamma)                         ///< [in] Ratio of specific heats
  {
    assert(std::isnormal(gamma));
    const MyFloat gammam1 = gamma - (MyFloat)1.0;
    MyFloat vSqd = (MyFloat)0.0;
    for (int k=0; k<ndim; k++) vSqd += (W[k]*W[k]);
    for (int k=0; k<ndim; k++) data[k] = W[irho]*W[k];
    data[irho]   = W[irho];
    data[ipress] = W[ipress]/gammam1 + (MyFloat)0.5*W[irho]*vSqd;
  }


  //===============================================================================================
  /// Constructor to convert conserved fluid vector (Q) to primitive fluid vector (W).
  //===============================================================================================
  template <int ndim>
  WFluidVector<ndim>::WFluidVector
   (const QFluidVector<ndim> &Q,                 ///< [in] Conserved fluid vector to be converted
    const MyFloat vol,                           ///< [in] Volume of particle
    const MyFloat gamma,                         ///< [in] Ratio of specific heats
    const MyFloat a)
  {
    assert(std::isnormal(gamma));
    assert(std::isnormal(vol));
    assert(std::isnormal(Q[irho]));
    const MyFloat gammam1 = gamma - (MyFloat) 1.0;
    const MyFloat invMass = (MyFloat) 1.0 / Q[irho];
    const MyFloat invVol = (MyFloat) 1.0 / vol;
    MyFloat pSqd = (MyFloat) 0.0;
    for (int k=0; k<ndim; k++) pSqd += (Q[k]*Q[k]);
    for (int k=0; k<ndim; k++) data[k] = Q[k]*invMass;
    data[irho] = Q[irho]*invVol;
    data[ipress] = (Q[ipress] - (MyFloat) 0.5*pSqd*invMass)*gammam1*invVol*pow(a,3*gammam1);
    //data[ipress] = (Q[ipress] - (MyFloat) 0.5*pSqd*invMass)*invMass*pow(data[irho],gamma);
  }


  //===============================================================================================
  /// Constructor to convert conserved fluid vector (U) to primitive fluid vector (W).
  //===============================================================================================
  template <int ndim>
  WFluidVector<ndim>::WFluidVector
   (const UFluidVector<ndim> &U,                 ///< [in] Conserved fluid vector to be converted
    const MyFloat gamma)                         ///< [in] Ratio of specific heats
  {
    assert(std::isnormal(U[irho]));
    const MyFloat gammam1 = gamma - (MyFloat) 1.0;
    const MyFloat invRho = (MyFloat) 1.0 / U[irho];
    MyFloat pSqd = (MyFloat) 0.0;
    for (int k=0; k<ndim; k++) pSqd += (U[k]*U[k]);
    for (int k=0; k<ndim; k++) data[k] = U[k]*invRho;
    data[irho] = U[irho];
    data[ipress] = (U[ipress] - (MyFloat) 0.5*pSqd*invRho)*gammam1;
  }

  //===============================================================================================
  /// Constructor to convert (Sph)P variables to primitive fluid vector (W).
  //===============================================================================================
  template <int ndim>
  WFluidVector<ndim>::WFluidVector
   (const MyFloat *Vel,                          ///< [in] velocity
    const MyLongDouble Density,                  ///< [in] density
    const MyFloat Pressure)                      ///< [in] pressure
  {
    assert(std::isnormal(Density));
    for (int k=0; k<ndim; k++) data[k] = Vel[k];
    data[irho] = Density;
    data[ipress] = Pressure;
  }


  //===============================================================================================
  /// Definitions of functions to compute the time derivative of the primitive fluid vector.
  /// Used in the MUSCL scheme to predict primitive properties at the half timestep.
  //===============================================================================================
#if defined(ONEDIM)
  template <>
  WFluidVector<1> WFluidVector<1>::ComputePrimitiveTimeDerivative
   (const MyFloat cSoundSqd,                     ///< [in] Sound speed squared
    const WFluidVector<1> &W,                    ///< [in] Primitive fluid vector
    const Vector<1> *gradW,                      ///< [in] (Slope-limited) gradients
    const Vector<1> Vel_smooth)                  ///< [in] smmothed velocity over neighbours
  {
    WFluidVector<1> Wdot;
    const MyFloat divV = gradW[ivx][ivx];
    const MyFloat invrho = (MyFloat) 1.0 / W[irho];
    Wdot[ivx]    = -(W[ivx]-Vel_smooth[0])*gradW[ivx][0] - gradW[ipress][0]*invrho;
    Wdot[irho]   = -(W[ivx]-Vel_smooth[0])*gradW[irho][0] - W[irho]*divV;
    Wdot[ipress] = -(W[ivx]-Vel_smooth[0])*gradW[ipress][0] - W[irho]*cSoundSqd*divV;
    return Wdot;
  }
#elif defined(TWODIMS)
  template <>
  WFluidVector<2> WFluidVector<2>::ComputePrimitiveTimeDerivative
   (const MyFloat cSoundSqd,                     ///< [in] Sound speed squared
    const WFluidVector<2> &W,                    ///< [in] Primitive fluid vector
    const Vector<2> *gradW,                      ///< [in] (Slope-limited) gradients
    const Vector<2> Vel_smooth)                  ///< [in] smmothed velocity over neighbours
  {
    WFluidVector<2> Wdot;
    const MyFloat divV = gradW[ivx][ivx] + gradW[ivy][ivy];
    const MyFloat invrho = (MyFloat) 1.0 / W[irho];
    Wdot[ivx]    = -(W[ivx]-Vel_smooth[0])*gradW[ivx][0] - (W[ivy]-Vel_smooth[1])*gradW[ivx][1] - gradW[ipress][0]*invrho;
    Wdot[ivy]    = -(W[ivx]-Vel_smooth[0])*gradW[ivy][0] - (W[ivy]-Vel_smooth[1])*gradW[ivy][1] - gradW[ipress][1]*invrho;
    Wdot[irho]   = -(W[ivx]-Vel_smooth[0])*gradW[irho][0] - (W[ivy]-Vel_smooth[1])*gradW[irho][1] - W[irho]*divV;
    Wdot[ipress] = -(W[ivx]-Vel_smooth[0])*gradW[ipress][0] - (W[ivy]-Vel_smooth[1])*gradW[ipress][1] - W[irho]*cSoundSqd*divV;
    return Wdot;
  }
#else
  template <>
  WFluidVector<3> WFluidVector<3>::ComputePrimitiveTimeDerivative
   (const MyFloat cSoundSqd,                     ///< [in] Sound speed squared
    const WFluidVector<3> &W,                    ///< [in] Primitive fluid vector
    const Vector<3> *gradW,                      ///< [in] (Slope-limited) gradients
    const Vector<3> Vel_smooth)                  ///< [in] smmothed velocity over neighbours
  {
    WFluidVector<3> Wdot;
    const MyFloat divV = gradW[ivx][ivx] + gradW[ivy][ivy] + gradW[ivz][ivz];
    const MyFloat invrho = (MyFloat) 1.0 / W[irho];
    Wdot[ivx]    = -(W[ivx]-Vel_smooth[0])*gradW[ivx][0] - (W[ivy]-Vel_smooth[1])*gradW[ivx][1] - (W[ivz]-Vel_smooth[2])*gradW[ivx][2] - gradW[ipress][0]*invrho;
    Wdot[ivy]    = -(W[ivx]-Vel_smooth[0])*gradW[ivy][0] - (W[ivy]-Vel_smooth[1])*gradW[ivx][1] - (W[ivz]-Vel_smooth[2])*gradW[ivy][2] - gradW[ipress][1]*invrho;
    Wdot[ivz]    = -(W[ivx]-Vel_smooth[0])*gradW[ivz][0] - (W[ivy]-Vel_smooth[1])*gradW[ivz][1] - (W[ivz]-Vel_smooth[2])*gradW[ivz][2] - gradW[ipress][2]*invrho;
    Wdot[irho]   = -(W[ivx]-Vel_smooth[0])*gradW[irho][0] - (W[ivy]-Vel_smooth[1])*gradW[irho][1] - (W[ivz]-Vel_smooth[2])*gradW[irho][2] - W[irho]*divV;
    Wdot[ipress] = -(W[ivx]-Vel_smooth[0])*gradW[ipress][0] - (W[ivy]-Vel_smooth[1])*gradW[ipress][1] - (W[ivz]-Vel_smooth[2])*gradW[ipress][2] - W[irho]*cSoundSqd*divV;
    return Wdot;
  }
#endif
  

#if defined(ONEDIM)
  template class QFluidVector<1>;
  template class UFluidVector<1>;
  template class WFluidVector<1>;
#elif defined(TWODIMS)
  template class QFluidVector<2>;
  template class UFluidVector<2>;
  template class WFluidVector<2>;
#else
  template class QFluidVector<3>;
  template class UFluidVector<3>;
  template class WFluidVector<3>;
#endif

}
//-------------------------------------------------------------------------------------------------
