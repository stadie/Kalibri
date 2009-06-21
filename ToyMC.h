#ifndef TOYMC_H
#define TOYMC_H

#include <vector>

#include "TH1F.h"
#include "TLorentzVector.h"
#include "TRandom.h"

class TTree;


//!  \brief Generate toy MC data
//!
//!  \author Hartmut Stadie
//!  \date   Mon Jun 30 11:00:00 CEST 2008
//!  $Id: ToyMC.h,v 1.12 2009/06/11 17:32:15 mschrode Exp $
// ----------------------------------------------------------------  
class ToyMC {

 private:
  //!  \brief Resolution model
  //! 
  //!  - For calibration
  //!    - 'Gauss': An energy dependent Gaussian resolution
  //!      \f[
  //!       \frac{\sigma}{p_{T}} = \frac{b_{0}}{\sqrt(p_{T})} \oplus b_{1}
  //!      \f]
  //!      
  //!    - 'Landau': An energy dependent Landau resolution
  //!
  //!  - For jet smearing:
  //!    - 'GaussUniform': A Gaussian resolution with a uniform
  //!      low energy tail.
  //!      The appropriate parametrizations are
  //!      - SmearParametrizationFermiTail
  //!  
  //!    - 'TwoGauss': Resolution given by a central Gaussian around
  //!      1 and a second Gaussian to model tails.
  //!      The appropriate parametrizations are
  //!      - SmearParametrizationStepGauss
  // ----------------------------------------------------------------  
  enum ResolutionModel { Gauss, Landau, GaussUniform, TwoGauss };


  //!  \brief Response for generation
  //!
  //!  The generated response
  //!  \f[
  //!   R(E^{\textrm{true}}_{T}) = \frac{E^{\textrm{meas}}_{T}}{E^{\textrm{true}}_{T}}
  //!  \f]
  //!  depends on the parameters \f$ A_{i} \f$ as defined in the config file
  //!  via the field 'ToyMC response parameters'. The response model is definded
  //!  via the field 'ToyMC response model'.
  //!
  //!  Possible response models are:
  //!    - 'Constant': Apply a constant response factor to the
  //!      hadronic part of the jet Et:
  //!      \f[ R = \frac{1}{A_{0}} \f]
  //!      This is the default.
  //!      The appropriate correction functions are
  //!      - ToyParametrization
  //!      - ToyJetParametrization
  //!      - ToyStepParametrization
  //!      - ToyStepJetParametrization
  //!
  //!    - 'L3': Level 3 response from JetMET group. In each event,
  //!      the response factor is calculated from the total photon pt
  //!      (i.e. \f$ E^{\textrm{true}}_{T} = P^{\gamma}_{T} \f$)
  //!      but applied only to the hadronic part of the jet Et:
  //!      \f[ R = A_{0}
  //!            - \frac{A_{1}}{\log^{A_{2}}_{10}(E^{\textrm{true}}_{T}) + A_{3}}
  //!            +\frac{A_{4}}{E^{\textrm{true}}_{T}}  \f]
  //!      The appropriate correction functions are
  //!      - L2L3JetParametrization
  //!      - L2L3JetParametrization2
  //!   
  //!    - 'SimpleInverse': In each event,
  //!      the response factor is calculated from the total photon pt
  //!      (i.e. \f$ E^{\textrm{true}}_{T} = P^{\gamma}_{T} \f$)
  //!      but applied only to the hadronic part of the jet Et:
  //!      \f[ R = 1 - \frac{A_{0}}{E^{\textrm{true}}_{T} + A_{1}} \f]
  //!      The appropriate correction function is
  //!      - ToySimpleInverseParametrization
  //!
  //!    - 'Flat'
  //!    - 'Exp'
  //!    - 'Slope'
  // ----------------------------------------------------------------  
  enum ResponseModel { Constant, L3, SimpleInverse, Flat, Slope, Exp };


  //!  \brief Truth pt spectrum
  // ----------------------------------------------------------------  
  enum TruthSpectrum { Uniform, PowerLaw };


  // Global variables
  TRandom*        mRandom;             //!< Random generator
  int             mType;               //!< Event type: Photonjet (1), Dijet (2)

  // Parameters for truth
  double          mMinEta;             //!< Minimum truth eta
  double          mMaxEta;             //!< Maximum truth eta
  double          mMinPt;              //!< Minimum truth pt
  double          mMaxPt;              //!< Maximum truth pt
  TruthSpectrum   mPtSpectrum;         //!< Truth pt spectrum
  TLorentzVector  mPinput;             //!< Stores the truth lorentz vector of the current event

  // Parameters for measurement 
  int             mChunks;
  double          mJetSpreadA;
  double          mJetSpreadB;
  bool            mNoOutOfCone;
  double          mMaxPi0Frac;
  double          mMaxEmf;

  ResponseModel   mResponseModel;      //!< Response model
  std::vector<double> mParResp;        //!< Parameters for Response
  TH1F          * mHistResp;           //!< For histogramed response

  ResolutionModel mResolutionModel;    //!< Resolution model
  std::vector<double> mParReso;        //!< Parameters for Respolution

  double          mSmearFactor;        //!< Combined smear factor from response and resolution
  bool            mSmearTowersIndividually;  //!< If true, mSmearTowersIndividually is determined individually for each tower, else for each jet



  void genInput();
  void calIds(float& eta, float &phi, int& ieta, int& iphi);
  void smearTower(double e, bool calcSmearFactor, float& te, float& tem, float& thad, float& tout, float& temtrue, 
		  float& thadtrue, float& touttrue);  
  int  splitJet(const TLorentzVector& jet ,float* et,float* eta,float * phi, int* ieta,int* iphi);
  void CalculateSmearFactor(double pt);


 public:
  ToyMC();
  ~ToyMC() {
    delete mRandom;
    if( mHistResp ) delete mHistResp;
  }
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
