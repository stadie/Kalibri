// $Id: SmearPhotonJet.cc,v 1.1 2009/06/11 17:29:25 mschrode Exp $

#include "SmearPhotonJet.h"



//!  \brief Get the negative log-likelihood of this event
//!  \return The negative log-likelihood of this event
// --------------------------------------------------
double SmearPhotonJet::chi2() const {
  TJet mess(GetMess());
  mess.pt = GetMess()->pt / GetTruth();
  
  return -1.*log( respPDF_(&mess) / GetTruth() ); // Need to divide by _truth to have probability (!= density)
}



//!  \brief Get the negative log-likelihood and the
//!         contributions to the 1. and 2. derivatives
//!
//!  Calculates the negative log-likelihood \f$ -\ln L \f$
//!  of this event. Moreover the contribution to the 
//!  first and second derivative ('temp_derivative1',
//!  'temp_derivative2') of the global \f$ -\ln L \f$
//!  function is calculated numerically and returned
//!  by reference, where
//!  \f[
//!   \frac{\partial (-\ln L) }{\partial p}
//!   = \sum \frac{\textrm{temp\_derivative1}}{2\epsilon}
//!  \f]
//!  and
//!  \f[
//!   \frac{\partial^{2} (-\ln L) }{\partial p^{2}}
//!   = \sum \frac{\textrm{temp\_derivative2}}{\epsilon^{2}}
//!  \f]
//!
//!  \param temp_derivative1 Pointer to first derivative contribution
//!  \param temp_derivative2 Pointer to second derivative contribution
//!  \param epsilon Step size \f$ \epsilon \f$ for derivative calculation
//!  \return The negative log-likelihood of this event
// --------------------------------------------------
double SmearPhotonJet::chi2_fast(double * temp_derivative1,
			     double * temp_derivative2,
			     double const epsilon) const {
  double f = chi2();

  int      idx;
  double * par;
  double   oldpar;
  double   temp1;
  double   temp2;

  // Vary parameters of response pdf
  idx = respPDF_.parIndex();
  par = respPDF_.firstPar();
  for(int i = 0; i < respPDF_.nPars(); i++) {
    oldpar = par[i];

    par[i] += epsilon;
    temp1   = chi2();

    par[i] -= 2.*epsilon;
    temp2   = chi2();

    par[i]  = oldpar;

    temp_derivative1[idx+i] += temp1 - temp2;
    temp_derivative2[idx+i] += temp1 + temp2 - 2*f;
  }

  return f;
}
