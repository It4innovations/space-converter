//#include "Mfm/RiemannSolver.hpp"
#include <math.h>
#include "RiemannSolver.hpp"


//-------------------------------------------------------------------------------------------------
namespace OpenGadget3
{

  //===============================================================================================
  /// Base constructor for Riemann solver class.  Sets all gamma const variables.
  //===============================================================================================
  template <int ndim>
  RiemannSolver<ndim>::RiemannSolver(Eos<ndim>* const _eos) :
    gamma(_eos->GetGamma()),
    invgamma(MyFloat(1.0)/_eos->GetGamma()),
    g1(MyFloat(0.5)*(_eos->GetGamma() - MyFloat(1.0))/_eos->GetGamma()),
    g2(MyFloat(0.5)*(_eos->GetGamma() + MyFloat(1.0))/_eos->GetGamma()),
    g3(MyFloat(2.0)*_eos->GetGamma()/(_eos->GetGamma() - MyFloat(1.0))),
    g4(MyFloat(2.0)/(_eos->GetGamma() - MyFloat(1.0))),
    g5(MyFloat(2.0)/(_eos->GetGamma() + MyFloat(1.0))),
    g6((_eos->GetGamma() - 1.0)/(_eos->GetGamma() + MyFloat(1.0))),
    g7(MyFloat(0.5)*(_eos->GetGamma() - MyFloat(1.0))),
    g8(_eos->GetGamma() - MyFloat(1.0)),
    g9(MyFloat(1.0)/(_eos->GetGamma() - MyFloat(1.0))),
    eos(_eos)
  {
    assert(std::isnormal(gamma));
  }

#if defined(ONEDIM)
  template class RiemannSolver<1>;
#elif defined(TWODIMS)
  template class RiemannSolver<2>;
#else
  template class RiemannSolver<3>;
#endif

}
//-------------------------------------------------------------------------------------------------
