//
//    Class for all events with one jet and truth informatio
//
//    first version: Hartmut Stadie 2008/12/14
//    $Id: JetWidthEvent.cc,v 1.3 2010/12/13 10:38:28 stadie Exp $
//   

#include "JetWidthEvent.h"
#include "Parameters.h"
#include "Parametrization.h"

JetWidthEvent::~JetWidthEvent() 
{ 
  delete jet_;
}

double JetWidthEvent::chi2() const
{
  const MeanWidthParametrization* mwp = dynamic_cast<const MeanWidthParametrization*>(jet_->f()->parametrization());
  
  double mu = mwp->widthMean(jet_,jet_->f()->firstPar());
  double sigma = mwp->widthSigma(jet_,jet_->f()->firstPar());
  double jw = 0.5 * (jet_->momentPhiPhi()+jet_->momentEtaEta());
  //std::cout << "mu:" << mu <<" sigma:" << sigma << '\n';
  double diff = (jw-mu)/sigma;
  return weight_ * Event::scaleResidual(diff*diff);
}

double JetWidthEvent::chi2_fast(double * temp_derivative1, 
				double * temp_derivative2,  
				double * temp_derivative3, 
				double * temp_derivative4,
				const double* epsilon) const
{
  static const double log2pi = log(2*M_PI)+20;
  const Function* f = jet_->f();
  const MeanWidthParametrization* mwp = dynamic_cast<const MeanWidthParametrization*>(f->parametrization());
  
  double mu = mwp->widthMean(jet_,f->firstPar());
  double sigma = mwp->widthSigma(jet_,f->firstPar());
  double jw = 0.5 * (jet_->momentPhiPhi()+jet_->momentEtaEta());
  //std::cout << "mu:" << mu <<" sigma:" << sigma << '\n';
  double diff = (jw-mu)/sigma;
  double chi2 = weight_ * Event::scaleResidual(log2pi + 2*log(sigma)+diff*diff);
  //std::cout << "diff:" << diff << "  " << chi2 << '\n';
  assert(chi2 == chi2);
  double temp1,temp2;
  for(int i = 0 ; i < 10 ; ++i) {
    double orig = f->firstPar()[i];
    int parid = f->parIndex() + i;
    f->firstPar()[i] += epsilon[parid];
    mu = mwp->widthMean(jet_,f->firstPar());
    sigma = mwp->widthSigma(jet_,f->firstPar());
    temp1 = (jw-mu)/sigma;
    temp1*= temp1;   
    temp1 = weight_ * Event::scaleResidual(log2pi+2*log(sigma)+temp1);
    if(std::isinf(temp1)) {
      std::cout << "Grrh:" << mu << ", " << sigma << ", " << jw << '\n';
    }
    assert(temp1 == temp1);
    f->firstPar()[i] = orig - epsilon[parid];
    mu = mwp->widthMean(jet_,f->firstPar());
    sigma = mwp->widthSigma(jet_,f->firstPar());
    temp2 = (jw-mu)/sigma; 
    temp2*= temp2; 
    temp2 = weight_ * Event::scaleResidual(log2pi+2*log(sigma)+temp2);
    if(std::isinf(temp2)) {
      std::cout << "Grrh:" << mu << ", " << sigma << ", " << jw << '\n';
      std::cout << "i = " << i << " orig = " << orig << "  eps = " << epsilon[parid] << '\n';
    }
    assert(temp2 == temp2);
    f->firstPar()[i] = orig;
    double td1 = temp1 - temp2;
    if(td1 != td1) {
      std:: cout << "td = " << td1 << " temp1 = " << temp1 << " temp2 = " 
		 << temp2 << '\n';
    }
    assert(td1 == td1);
    temp_derivative1[parid] += td1; // for 1st derivative
    temp_derivative2[parid] += (temp1 - chi2 + temp2 - chi2); // for 2nd derivative
  }
  return chi2plots_ = chi2;
}
