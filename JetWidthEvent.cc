//
//    Class for all events with one jet and truth informatio
//
//    first version: Hartmut Stadie 2008/12/14
//    $Id: JetWidthEvent.cc,v 1.24 2010/11/01 15:47:41 stadie Exp $
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

double JetWidthEvent::chi2_fast(long double * temp_derivative1, 
				long double * temp_derivative2, 
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
    assert(temp1 == temp1);
    f->firstPar()[i] = orig - epsilon[parid];
    mu = mwp->widthMean(jet_,f->firstPar());
    sigma = mwp->widthSigma(jet_,f->firstPar());
    temp2 = (jw-mu)/sigma;
    temp2*= temp2;
    temp2 = weight_ * Event::scaleResidual(log2pi+2*log(sigma)+temp2);
    assert(temp2 == temp2);
    f->firstPar()[i] = orig;
    temp_derivative1[parid] += (temp1 - temp2); // for 1st derivative
    temp_derivative2[parid] += (temp1 - chi2 + temp2 - chi2); // for 2nd derivative
  }
  return chi2plots_ = chi2;
}
