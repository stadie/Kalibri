//
// $Id: CalibData.cc,v 1.33 2010/05/19 13:34:48 stadie Exp $
//
#include "CalibData.h"


#include <map>
#include <cmath>

#include "Parameters.h"

double (*Event::scaleResidual)(double z2) = &Event::scaleNone;

//!  Scale the normalized, squared residual
//!  \f$ z^{2} = \chi^{2}/\textrm{weight} \f$
//!  using the Cauchy-Function
//!  \f$ z^{2} \rightarrow c^{2}\ln( 1 + (z/c)^{2} ) \f$
//!
//!  \param z2 Normalized and squared residual
//!  \return Scaled residual
double Event::scaleCauchy(double const z2)
{
  double const c = 2.3849;
  return (c*c) * log( 1 + z2*(1.0/(c*c)) );
}

//!  Scale the normalized, squared residual
//! \f$ z^{2} = \chi^{2}/\textrm{weight} \f$
//!  using the Huber-Function
//!  \f[
//!  z^{2} \rightarrow 
//!  \left\{
//!     \begin{array}{ll}
//!        z & \textrm{for } |z| <= c \\ c ( 2|z| - c ) & \textrm{for } |z| > c
//!      \end{array}
//!    \right.
//!  \f]
//!
//!  \param z2 Normalized and squared residual
//!  \return Scaled residual
double Event::scaleHuber(double const z2)
{
  static double const c = 1.345;
  double const z = sqrt(z2);
  return (  std::abs(z) <= c  ?  z2  :  c*(2.*std::abs(z) - c)  );
}


//!  \brief Cut on residuals
//!
//!  discards events with $|residual| > sqrt(c2) \sigma$
//!
//!  \param z2 Normalized and squared residual
//!  \return Scaled residual
double Event::scaleTukey(const double z2)
{
  const double c2 = 16;
  if(z2 > c2) return 0;
  double w = 1-z2/c2;
  return w*w * z2;
}

