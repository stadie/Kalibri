//
// Original Author:  Hartmut Stadie
//         Created:  Mon Jun 30 11:00:00 CEST 2008
// $Id: toy.cc,v 1.4 2008/07/16 12:20:14 auterman Exp $
//
#include "ToyMC.h"


int main() {
  ToyMC* mc = new ToyMC();
  mc->mMinEta         = -1.5;
  mc->mMaxEta         =  1.5;
  mc->mMinPt          = 20;
  mc->mMaxPt          = 200;
  //PtSpectrum:
  //  - Uniform: flat in pt between minPt and maxPt
  //  - PowerLaw: falling pt spectrum from minPt with p_T^{-2.5}
  mc->mPtSpectrum     = ToyMC::PowerLaw;
  mc->mMaxPi0Frac     =  0.5;
  mc->mMaxEmf         =  0.0;
  mc->mTowConst       =  1.3;
  mc->mResoStochastic =  1.3;
  mc->mResoNoise      =  0.056;
  //mc->mTowConst       =  1.00;
  //mc->mResoStochastic =  0.0;
  //mc->mResoNoise      =  0.0;
  //simulated out-of-cone correction factor 
  //R = 1/(1-exp(-0.5(A+bE))) = 1 + exp(-0.5(A+BE)) + exp(-(A+BE)) +...
  mc->mJetSpreadA     =  -2 * log(1-0.995);//99.5% in 0.5 cone
  mc->mJetSpreadB     =  0.0;
  mc->mNoOutOfCone    =  true;
  //settings for symmetric distributions 
  mc->mModel          = ToyMC::Gauss;
  //setting for flat distribution (noise)
  //mc->mModel          = ToyMC::Flat;
  //settings for asymmetric distributions (noise)
  //mc->mModel          = ToyMC::Exp;
  //mc->mModel          = ToyMC::Slope;
  //mc->makeTrackCluster("trackcluster.root", 50000);
  mc->makePhotonJet("input/toy_photonjet.root",20000);
  
  return 0;
}

