#ifndef _OPENGADGET3_EOS_HPP_
#define _OPENGADGET3_EOS_HPP_


//#include "../CodeBase/allvars.h"
#include "../System/assert.h"
#include <cmath>
#include "../CodeBase/precision.h"


//-------------------------------------------------------------------------------------------------
namespace OpenGadget3
{

  //===============================================================================================
  /// \brief   Base class for defining all equations-of-state classes in hydro computations.
  /// \author  D. A. Hubber
  /// \date    12/03/2018
  //===============================================================================================
  template <int ndim>
  class Eos
  {
  protected:

    const MyFloat gamma;                         ///< Ratio of specific heats
    const MyFloat gammam1;                       ///< gamma - 1
    const MyFloat invgammam1;                    ///< 1 / (gamma - 1)


  public:

    Eos(const MyFloat _gamma) :
     gamma(_gamma),
     gammam1(_gamma - MyFloat(1.0)),
     invgammam1(MyFloat(1.0) / (_gamma - MyFloat(1.0)))
    {
      assert(_gamma > MyFloat(1.0));
    };
    virtual ~Eos() {};

    virtual MyFloat Pressure(const MyFloat rho, const MyFloat u, const MyFloat a3inv) const = 0;
    virtual MyFloat SoundSpeed(const MyFloat rho, const MyFloat u, const MyFloat a3inv) const = 0;
    virtual MyFloat SpecificInternalEnergy(const MyFloat rho, const MyFloat u) const = 0;
    virtual MyFloat Temperature(const MyFloat rho, const MyFloat u) const = 0;

    inline MyFloat GetGamma() const {return gamma;}
    inline MyFloat GetGammaMinusOne() const {return gammam1;}
    inline MyFloat GetInvGammaMinusOne() const {return invgammam1;}

    inline MyFloat GetUFromEntropy(const MyFloat rho, const MyFloat entropy, const MyFloat a3inv) const
    {
      return invgammam1*entropy*pow(rho*a3inv, gammam1);
    }
    inline MyFloat GetEntropyFromU(const MyFloat rho, const MyFloat u, const MyFloat a3inv) const
    {
      return gammam1*u/pow(rho*a3inv, gammam1);
    }
    
    inline MyFloat GetDtQFromDtEntropy(const MyFloat dtentropy,
				       const MyFloat entropy, const MyFloat u,
				       const MyFloat mass,
				       const MyFloat divvel, const MyFloat hubble_a) const
    {
      return u * (dtentropy / entropy - gammam1*divvel) * hubble_a * mass;
    }
    inline MyFloat GetDtUFromDtEntropy(const MyFloat dtentropy,
				       const MyFloat entropy, const MyFloat u,
				       const MyFloat divvel, const MyFloat hubble_a,
				       const int comoving_integration_on) const
    {
      MyFloat DtU = u * (dtentropy / entropy - gammam1*divvel);
      if (comoving_integration_on)
	DtU -= u * 3*gammam1// *hubble_a
	  ;
      return DtU;
    }
    inline MyFloat GetDtEntropyFromDtU(const MyFloat dtu,
				       const MyFloat entropy, const MyFloat u,
				       const MyFloat divvel, const MyFloat hubble_a,
				       const int comoving_integration_on) const
    {
      MyFloat DtEntropy = entropy * dtu / u + gammam1*divvel;
      if (comoving_integration_on)
	DtEntropy += entropy * 3*gammam1// *hubble_a
	  ;
      return DtEntropy;
    }

  };



  //===============================================================================================
  /// \brief   Adiabatic equation-of-state class.
  /// \author  D. A. Hubber
  /// \date    12/03/2018
  //===============================================================================================
  template <int ndim>
  class AdiabaticEos final : public Eos<ndim>
  {
  protected:

    using Eos<ndim>::gamma;
    using Eos<ndim>::gammam1;


  public:

    // Constructor for setting all adiabatic parameters
    AdiabaticEos(const MyFloat _gamma) : Eos<ndim>(_gamma) {};
    ~AdiabaticEos() {};

    virtual MyFloat Pressure(const MyFloat rho, const MyFloat u, const MyFloat a3inv) const override
    {
      return gammam1*rho*u/pow(a3inv,gammam1);
    }
    virtual MyFloat SoundSpeed(const MyFloat rho, const MyFloat u, const MyFloat a3inv) const override
    {
      return sqrt(gamma*gammam1*u/pow(a3inv,gammam1));
    }
    virtual MyFloat SpecificInternalEnergy(const MyFloat rho, const MyFloat u) const override
    {
      return u;
    }
    virtual MyFloat Temperature(const MyFloat rho, const MyFloat u) const override
    {
      return gammam1*u;
    }

  };



  //===============================================================================================
  /// \brief   Isothermal equation-of-state class definition.
  /// \author  D. A. Hubber
  /// \date    06/06/2018
  //===============================================================================================
  template <int ndim>
  class IsothermalEos final : public Eos<ndim>
  {
  protected:

    using Eos<ndim>::gamma;
    using Eos<ndim>::gammam1;
    using Eos<ndim>::invgammam1;

    const MyFloat temp0;                         ///< Isothermal temperature
    const MyFloat muBar;                         ///< Mean gas particle mass
    const MyFloat invMuBar;                      ///< 1 / muBar
    const MyFloat cSound;                        ///< Isothermal sound speed


  public:

    // Constructor for setting all isothermal parameters
    IsothermalEos(const MyFloat _gamma, const MyFloat _temp0, const MyFloat _muBar) :
     Eos<ndim>(_gamma),
     temp0(_temp0),
     muBar(_muBar),
     invMuBar(MyFloat(1.0) / _muBar),
     cSound(sqrt(_temp0/_muBar))
    {
      assert(std::isnormal(_temp0));
      assert(std::isnormal(_muBar));
    }

    ~IsothermalEos() {};

    virtual MyFloat Pressure(const MyFloat rho, const MyFloat u, const MyFloat a3inv) const override
    {
      return rho*temp0*invMuBar;
    }
    virtual MyFloat SoundSpeed(const MyFloat rho, const MyFloat u, const MyFloat a3inv) const override
    {
      return cSound;
    }
    virtual MyFloat SpecificInternalEnergy(const MyFloat rho, const MyFloat u) const override
    {
      return temp0*invgammam1*invMuBar;
    }
    virtual MyFloat Temperature(const MyFloat rho, const MyFloat u) const override
    {
      return temp0;
    }

  };

}
//-------------------------------------------------------------------------------------------------
#endif
