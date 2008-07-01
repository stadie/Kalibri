//
// Original Author:  Hartmut Stadie
//         Created:  Mon Jun 30 11:00:00 CEST 2008
// $Id: toy.cc,v 1.1 2008/06/30 13:15:27 stadie Exp $
//
#include "ToyMC.h"


int main() {
  ToyMC* mc = new ToyMC();
  mc->mMinEta         = -2.5;
  mc->mMaxEta         =  2.5;
  mc->mMinPt          = 30;
  mc->mMaxPt          = 400;
  mc->mTowConst       =  1.00;
  mc->mResoStochastic =  1.25;
  mc->mResoNoise      =  0.056;
  mc->mJetSpread      =  0.10;
  mc->mNoOutOfCone    = true;
  //settings for symmetric distributions 
  //mc->mModel          = ToyMC::Gauss;
  //setting for flat distribution (noise)
  mc->mModel          = ToyMC::Flat;
  /*
    //settings for asymmetric distributions
    mc->mModel          = ToyMC::Landau;
    mc->mChunks         = 1;
  */
  //mc->makeTrackCluster("trackcluster.root", 50000);
  mc->makePhotonJet("input/toy_photonjet.root",5000);
  
  return 0;
}

