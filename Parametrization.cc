//
//  $Id: Parametrization.h,v 1.50 2009/11/06 11:59:51 mschrode Exp $
//
#include "Parametrization.h"


#include "TH1D.h"
#include "TMath.h"


SmearStepGaussInter::SmearStepGaussInter(double tMin, double tMax, double rMin, double rMax, int rNBins, double ptDijetMin, double ptDijetMax, const std::vector<double>& rParScales, const std::vector<double>& gaussPar)
  : Parametrization(0,rNBins+4,0,1),
    tMin_(tMin),
    tMax_(tMax),
    rMin_(rMin),
    rMax_(rMax),
    ptDijetMin_(ptDijetMin),
    ptDijetMax_(ptDijetMax),
    nStepPar_(rNBins),
    binWidth_((rMax_ - rMin_)/nStepPar_),
    respParScales_(rParScales),
    gaussPar_(gaussPar)
{
  for(int i = 0; i < nStepPar_+1; i++) {
    binCenters_.push_back( rMin_ + (0.5+i)*binWidth_ );
  }
  assert( 0.0 <= tMin_ && tMin_ < tMax_ );
  assert( 0.0 <= rMin_ && rMin_ < rMax_ );
  assert( 0.0 <= ptDijetMin_ && ptDijetMin_ < ptDijetMax_ );
  assert( respParScales_.size() >= nJetPars() );
  assert( gaussPar_.size() >= 1 );
  
  print();
  
  // Integral over dijet resolution for truth pdf
  ptDijetCutInt_ = new TH1D("norm","",400,tMin_,tMax_);
}

SmearStepGaussInter::~SmearStepGaussInter() 
{ 
  binCenters_.clear(); 
  delete ptDijetCutInt_; 
}
  
void SmearStepGaussInter::update(const double * par)
{
  std::cout << "'" << name() << "': Updating ptDijet cut integral... ";
  ptDijetCutInt_->Reset();
  
  for(int bin = 1; bin < ptDijetCutInt_->GetNbinsX(); bin++) {
    Measurement x;
    x.pt = ptDijetCutInt_->GetBinCenter(bin);
    double integral = 0.;
    int nSteps = 400;
    double dPtDijet = (ptDijetMax_ - ptDijetMin_) / nSteps;
    for(int i = 0; i < nSteps; i++) {
      double ptDijet = ptDijetMin_ + i*dPtDijet;
      x.E = ptDijet / x.pt;
      double prob = correctedJetEt(&x,par);
      prob *= x.pt;
	integral += prob;
    }
    integral /= sqrt(2.);
    integral *= dPtDijet;
    ptDijetCutInt_->SetBinContent(bin,integral);
  }
  std::cout << "ok\n";
}


double SmearStepGaussInter::correctedGlobalJetEt(const Measurement *x,const double *par) const {
  // Norm of probability of dijet event configuration
  double norm = 0.;
  for(int bin = 1; bin < ptDijetCutInt_->GetNbinsX(); bin++) {
      double pt = ptDijetCutInt_->GetBinCenter(bin);
      norm += pow(pt,2-par[0]) * ptDijetCutInt_->GetBinContent(bin) * ptDijetCutInt_->GetXaxis()->GetBinWidth(1);
  }
  
  double p = 0.;
  if( norm ) {
    p = truthPDF(x->pt,par[0]);
    p /= norm;
  }
  
  return p;
}
 
double  SmearStepGaussInter::truthPDF(double pt, double n) const
{
  double prob = 0.;
  if( tMin_ < pt && pt < tMax_ ) {
    prob = 1. / pow( pt, n );
    int bin = ptDijetCutInt_->FindBin(pt);
    prob *= ptDijetCutInt_->GetBinContent(bin);
  }
  
  return prob;
}
