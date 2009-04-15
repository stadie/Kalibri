#ifndef TOYMC_H
#define TOYMC_H

#include "TLorentzVector.h"

class TRandom;
class TTree;


//!  \brief Generate toy MC data
//!
//!  \author Hartmut Stadie
//!  \date   Mon Jun 30 11:00:00 CEST 2008
//!  $Id: ToyMC.h,v 1.10 2009/04/06 14:51:05 mschrode Exp $
// ----------------------------------------------------------------  
class ToyMC {

 private:
  //!  \brief Resolution model
  // ----------------------------------------------------------------  
  enum Model { Gauss, Landau, Flat, Exp, Slope };


  //!  \brief Response for generation
  //!
  //!  The generated response
  //!  \f[
  //!   R(E^{\textrm{true}}_{T}) = \frac{E^{\textrm{meas}}_{T}}{E^{\textrm{true}}_{T}}
  //!  \f]
  //!  depends on the parameters \f$ A_{i} \f$ as defined in the config file
  //!  via the field 'ToyMC tower const'. The parameter values are
  //!  internally stored in the array mTowConst.
  //!
  //!  Possible response models are:
  //! 
  //!  - 'Constant': Apply a constant response factor to the
  //!    hadronic part of the jet Et:
  //!    \f[ R = \frac{1}{A_{0}} \f]
  //!    This is the default.
  //!    The appropriate correction functions are
  //!    - ToyParametrization
  //!    - ToyJetParametrization
  //!    - ToyStepParametrization
  //!    - ToyStepJetParametrization
  //!
  //!  - 'L3': Level 3 response from JetMET group. In each event,
  //!    the response factor is calculated from the total photon pt
  //!    (i.e. \f$ E^{\textrm{true}}_{T} = P^{\gamma}_{T} \f$)
  //!    but applied only to the hadronic part of the jet Et:
  //!    \f[ R = A_{0}
  //!          - \frac{A_{1}}{\log^{A_{2}}_{10}(E^{\textrm{true}}_{T}) + A_{3}}
  //!          +\frac{A_{4}}{E^{\textrm{true}}_{T}}  \f]
  //!    The appropriate correction functions are
  //!    - L2L3JetParametrization
  //!    - L2L3JetParametrization2
  //! 
  //!  - 'SimpleInverse': In each event,
  //!    the response factor is calculated from the total photon pt
  //!    (i.e. \f$ E^{\textrm{true}}_{T} = P^{\gamma}_{T} \f$)
  //!    but applied only to the hadronic part of the jet Et:
  //!    \f[ R = 1 - \frac{A_{0}}{E^{\textrm{true}}_{T} + A_{1}} \f]
  //!    The appropriate correction function is
  //!    - ToySimpleInverseParametrization
  // ----------------------------------------------------------------  
  enum Response { Constant, L3, SimpleInverse };


  //!  \brief Photon (truth) spectrum
  // ----------------------------------------------------------------  
  enum Spectrum { Uniform, PowerLaw };


  double          mMinEta;             //!< Minimum photon eta
  double          mMaxEta;             //!< Maximum photon eta
  double          mMinPt;              //!< Minimum photon pt
  double          mMaxPt;              //!< Maximum photon pt
  Spectrum        mPtSpectrum;         //!< Photon pt spectrum
  double          mTowConst[5];        //!< Parameters for Response
  double          mResoStochastic;     //!< Constant of stochastic term in resolution
  double          mResoNoise;          //!< Constant of noise term in resolution
  double          mJetSpreadA;
  double          mJetSpreadB;
  bool            mNoOutOfCone;
  Model           mModel;
  int             mChunks;
  double          mMaxPi0Frac;
  double          mMaxEmf;
  TLorentzVector  mPinput;             //!< Stores the truth lorentz vector of the current event
  TRandom*        mRandom;             //!< Random generator
  Response        mResponse;           //!< Response model


  void genInput();
  void calIds(float& eta, float &phi, int& ieta, int& iphi);
  void smearTower(double e, float& te, float& tem, float& thad, float& tout, float& temtrue, 
		  float& thadtrue, float& touttrue);  
  int  splitJet(const TLorentzVector& jet ,float* et,float* eta,float * phi, int* ieta,int* iphi);


 public:
  ToyMC();
  ~ToyMC() {}
  int generatePhotonJetTree(TTree *tree, int nevents);
  int generateTrackClusterTree(TTree *tree, int nevents);
  int generateDiJetTree(TTree* CalibTree, int nevents);
  int makeTrackCluster(const char* filename, int nevents);
  int makePhotonJet(const char* filename, int nevents);
  int makeDiJet(const char* filename, int nevents);
  void init(const std::string& configfile);
  void print() const;
};

#endif
