//
// Original Author:  Hartmut Stadie
//         Created:  Mon Jun 30 11:00:00 CEST 2008
// $Id: ToyMC.h,v 1.1 2008/06/30 13:15:27 stadie Exp $
//
#ifndef TOYMC_H
#define TOYMC_H

#include "TLorentzVector.h"

class TRandom;
class TTree;

class ToyMC {
public:
  enum Model { Gauss, Landau, Flat };
  double mMinEta, mMaxEta;
  double mMinPt, mMaxPt;
  double mTowConst;
  double mResoStochastic,mResoNoise;
  double mJetSpread;
  bool   mNoOutOfCone;
  Model  mModel;
  int    mChunks;
private: 
  TRandom* mRandom;
  TLorentzVector mPinput;

  void genInput();
  void calIds(float& eta, float &phi, int& ieta, int& iphi);
  void smearTower(double e, float& te, float& tem, float& thad, float& tout);  
  int splitJet(const TLorentzVector& jet ,float* et,float* eta,float * phi, int* ieta,int* iphi);
 public:
  ToyMC();
  ~ToyMC() {}
  int generatePhotonJetTree(TTree *tree, int nevents);
  int generateTrackClusterTree(TTree *tree, int nevents);
  int makeTrackCluster(const char* filename, int nevents);
  int makePhotonJet(const char* filename, int nevents);
};

#endif
