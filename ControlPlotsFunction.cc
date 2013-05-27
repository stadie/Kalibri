 // $Id: ControlPlotsFunction.cc,v 1.47 2013/05/21 13:49:37 kirschen Exp $

#include "ControlPlotsFunction.h"

#include "CalibMath.h"
#include "CorFactors.h"
#include "Jet.h"
#include "JetTruthEvent.h"
#include "TwoJetsPtBalanceEvent.h"
#include "TVector2.h"

/*
//!  Will return 0 if \p type does not exist in \p types_
// ----------------------------------------------------------------   
double ControlPlotsFunction::yValue(const Event * evt, ControlPlotsConfig::CorrectionType type) const {
  double yValue = 0.;
  std::map<ControlPlotsConfig::CorrectionType,Function>::const_iterator it = yFuncs_.find(type);
  if( it != yFuncs_.end() ) yValue = (this->*(it->second))(evt); 

  return yValue;
}
*/


// ----------------------------------------------------------------   
void ControlPlotsFunction::addYFunction(ControlPlotsConfig::CorrectionType type, Function func) {
  // Store correction types in a vector for easy access
  types_.push_back(type);
  // Store correction function
  std::map<ControlPlotsConfig::CorrectionType,Function>::const_iterator it = yFuncs_.find(type);
  if( it == yFuncs_.end() ) {
    yFuncs_[type] = func;
  } else {
    //std::cerr << "WARNING: Adding function for already existing CorrectionType '" << type << "'\n";
  }
}

//!  \brief Returns the number of PU events 
//!
//!  The \p Event \p evt has to be of type \p JetTruthEvent.
//!  Implements \p Function.
// ----------------------------------------------------------------   
double ControlPlotsFunction::jetTruthEventNPU(const Event *evt) const {
  const JetTruthEvent * jte = static_cast<const JetTruthEvent*>(evt);
  return jte->nPU();
}


//!  \brief Returns the number of truth PU events 
//!
//!  The \p Event \p evt has to be of type \p JetTruthEvent.
//!  Implements \p Function.
// ----------------------------------------------------------------   
double ControlPlotsFunction::jetTruthEventNPUTruth(const Event *evt) const {
  const JetTruthEvent * jte = static_cast<const JetTruthEvent*>(evt);
  return jte->nPUTruth();
}

//!  \brief Returns the number of PU events 
//!
//!  The \p Event \p evt has to be of type \p JetTruthEvent.
//!  Implements \p Function.
// ----------------------------------------------------------------   
double ControlPlotsFunction::jetTruthEventRho(const Event *evt) const {
  const JetTruthEvent * jte = static_cast<const JetTruthEvent*>(evt);
  return jte->rho();
}


//!  \brief Returns #eta of the jet
//!
//!  The \p Event \p evt has to be of type \p JetTruthEvent.
//!  Implements \p Function.
// ----------------------------------------------------------------   
double ControlPlotsFunction::jetTruthEventJetEta(const Event *evt) const {
  const JetTruthEvent * jte = static_cast<const JetTruthEvent*>(evt);
  return jte->jet()->eta();
}

//!  \brief Returns |#eta| of the jet
//!
//!  The \p Event \p evt has to be of type \p JetTruthEvent.
//!  Implements \p Function.
// ----------------------------------------------------------------   
double ControlPlotsFunction::jetTruthEventJetAbsEta(const Event *evt) const {
  const JetTruthEvent * jte = static_cast<const JetTruthEvent*>(evt);
  return std::abs(jte->jet()->eta());
}

//!  \brief Returns p_{T} of the jet
//!
//!  The \p Event \p evt has to be of type \p JetTruthEvent.
//!  Implements \p Function.
// ----------------------------------------------------------------   
double ControlPlotsFunction::jetTruthEventJetPt(const Event *evt) const {
  const JetTruthEvent * jte = static_cast<const JetTruthEvent*>(evt);
  return jte->jet()->pt();
}

//!  \brief Returns ECal fraction of the jet
//!
//!  The \p Event \p evt has to be of type \p JetTruthEvent.
//!  Implements \p Function.
// ----------------------------------------------------------------   
double ControlPlotsFunction::jetTruthEventJetEMF(const Event *evt) const {
  const JetTruthEvent * jte = static_cast<const JetTruthEvent*>(evt);
  return jte->jet()->emf();
}


//!  \brief Returns the #phi #phi moment of the jet
//!
//!  The \p Event \p evt has to be of type \p JetTruthEvent.
//!  Implements \p Function.
// ----------------------------------------------------------------   
double ControlPlotsFunction::jetTruthEventJetMomentPhiPhi(const Event *evt) const {
  const JetTruthEvent * jte = static_cast<const JetTruthEvent*>(evt);
  return jte->jet()->momentPhiPhi();
}


//!  \brief Returns the #eta #eta moment of the jet
//!
//!  The \p Event \p evt has to be of type \p JetTruthEvent.
//!  Implements \p Function.
// ----------------------------------------------------------------   
double ControlPlotsFunction::jetTruthEventJetMomentEtaEta(const Event *evt) const {
  const JetTruthEvent * jte = static_cast<const JetTruthEvent*>(evt);
  return jte->jet()->momentEtaEta();
}

//!  \brief Returns the #eta #eta moment of the jet
//!
//!  The \p Event \p evt has to be of type \p JetTruthEvent.
//!  Implements \p Function.
// ----------------------------------------------------------------   
double ControlPlotsFunction::jetTruthEventJetMeanMoment(const Event *evt) const {
  const JetTruthEvent * jte = static_cast<const JetTruthEvent*>(evt);
  return 0.5 * (jte->jet()->momentEtaEta()+jte->jet()->momentPhiPhi());
}

//!  \brief Returns truth pt
//!
//!  The \p Event \p evt has to be of type \p JetTruthEvent.
//!  Implements \p Function.
// ----------------------------------------------------------------   
double ControlPlotsFunction::jetTruthEventTruthPt(const Event *evt) const {
  const JetTruthEvent * jte = static_cast<const JetTruthEvent*>(evt);
  return jte->truth();
}

//!  \brief Returns flavor of jet
//!
//!  The \p Event \p evt has to be of type \p JetTruthEvent.
//!  Implements \p Function.
// ----------------------------------------------------------------   
double ControlPlotsFunction::jetTruthEventJetFlavor(const Event *evt) const {
  const JetTruthEvent * jte = static_cast<const JetTruthEvent*>(evt);
  return jte->jet()->flavor();
}


//!  \brief Returns the DeltaR to the closest reco jet
//!
//!  The \p Event \p evt has to be of type  \p JetTruthEvent
//!  DeltaR to closest other reco jet is returned
//!  
//!  
// ----------------------------------------------------------------   
double ControlPlotsFunction::jetTruthEventJetClosestJetdR(const Event * evt) const {
  const JetTruthEvent * jte = static_cast<const JetTruthEvent*>(evt);
  assert(jte->jet()->closestJetdR()>0);
  return jte->jet()->closestJetdR();
}




//!  \brief Returns the jet response
//!
//!  The \p Event \p evt has to be of type \p JetTruthEvent.
//!  The response is defined as
//!  \f[  p^{jet}_{T} / p^{true}_{T}\f].
//!  Implements \p Function.
// ----------------------------------------------------------------   
double ControlPlotsFunction::jetTruthEventResponse(const Event * evt) const {
  const JetTruthEvent * jte = static_cast<const JetTruthEvent*>(evt);
  Jet * jet = static_cast<Jet*>(jte->mess());

  return jet->pt() / jte->truth();
}



//!  \brief Returns the corrected jet response
//!
//!  The \p Event \p evt has to be of type \p JetTruthEvent.
//!  The response is defined as
//!  \f[  p^{jet'}_{T} / p^{true}_{T}\f],
//!  where \f$ p^{jet'}_{T} \f$ is the jet's pt corrected by
//!  the Kalibri JEC.
//!  Implements \p Function.
// ----------------------------------------------------------------   
double ControlPlotsFunction::jetTruthEventResponseKalibriCorrected(const Event * evt) const {
  const JetTruthEvent * jte = static_cast<const JetTruthEvent*>(evt);
  Jet * jet = static_cast<Jet*>(jte->mess());

  return jet->correctedEt() / jte->truth();
}



//!  \brief Returns the corrected jet response
//!
//!  The \p Event \p evt has to be of type \p JetTruthEvent.
//!  The response is defined as
//!  \f[  p^{jet'}_{T} / p^{true}_{T}\f],
//!  where \f$ p^{jet'}_{T} \f$ is the jet's pt corrected by
//!  the JetMET L2L3 JEC.
//!  Implements \p Function.
// ----------------------------------------------------------------   
double ControlPlotsFunction::jetTruthEventResponseL2L3Corrected(const Event * evt) const {
  const JetTruthEvent * jte = static_cast<const JetTruthEvent*>(evt);
  Jet * jet = static_cast<Jet*>(jte->mess());

  return jet->corFactors().getL2L3() * jet->pt() / jte->truth();
}

//!  \brief Returns the corrected jet response
//!
//!  The \p Event \p evt has to be of type \p JetTruthEvent.
//!  The response is defined as
//!  \f[  p^{jet'}_{T} / p^{true}_{T}\f],
//!  where \f$ p^{jet'}_{T} \f$ is the jet's pt corrected by
//!  the JetMET L2L3Res JEC.
//!  Implements \p Function.
// ----------------------------------------------------------------   
double ControlPlotsFunction::jetTruthEventResponseL2L3ResCorrected(const Event * evt) const {
  const JetTruthEvent * jte = static_cast<const JetTruthEvent*>(evt);
  Jet * jet = static_cast<Jet*>(jte->mess());

  return jet->corFactors().getL2L3Res() * jet->pt() / jte->truth();
}

//!  \brief Returns the corrected jet response
//!
//!  The \p Event \p evt has to be of type \p JetTruthEvent.
//!  The response is defined as
//!  \f[  p^{jet'}_{T} / p^{true}_{T}\f],
//!  where \f$ p^{jet'}_{T} \f$ is the jet's pt corrected by
//!  the JetMET L2L3L4 JEC.
//!  Implements \p Function.
// ----------------------------------------------------------------   
double ControlPlotsFunction::jetTruthEventResponseL2L3L4Corrected(const Event * evt) const {
  const JetTruthEvent * jte = static_cast<const JetTruthEvent*>(evt);
  Jet * jet = static_cast<Jet*>(jte->mess());

  return jet->corFactors().getL2L3() * jet->corFactors().getL4() * jet->pt() / jte->truth();
}


//!  \brief Returns the corrected jet response
//!
//!  The \p Event \p evt has to be of type \p JetTruthEvent.
//!  The response is defined as
//!  \f[  p^{jet'}_{T} / p^{true}_{T}\f],
//!  where \f$ p^{jet'}_{T} \f$ is the jet's pt corrected by
//!  the JetMET L2L3ResL4 JEC.
//!  Implements \p Function.
// ----------------------------------------------------------------   
double ControlPlotsFunction::jetTruthEventResponseL2L3ResL4Corrected(const Event * evt) const {
  const JetTruthEvent * jte = static_cast<const JetTruthEvent*>(evt);
  Jet * jet = static_cast<Jet*>(jte->mess());

  return jet->corFactors().getL2L3Res() * jet->corFactors().getL4() * jet->pt() / jte->truth();
}


//!  \brief Returns the corrected jet response
//!
//!  The \p Event \p evt has to be of type \p JetTruthEvent.
//!  The response is defined as
//!  \f[  p^{jet'}_{T} / p^{true}_{T}\f],
//!  where \f$ p^{jet'}_{T} \f$ is the jet's pt corrected by
//!  the JetMET L1L2L3 JEC.
//!  Implements \p Function.
// ----------------------------------------------------------------   
double ControlPlotsFunction::jetTruthEventResponseL1L2L3Corrected(const Event * evt) const {
  const JetTruthEvent * jte = static_cast<const JetTruthEvent*>(evt);
  Jet * jet = static_cast<Jet*>(jte->mess());

  return jet->corFactors().getL1()*jet->corFactors().getL2L3() * jet->pt() / jte->truth();
}


//!  \brief Returns the corrected jet response
//!
//!  The \p Event \p evt has to be of type \p JetTruthEvent.
//!  The response is defined as
//!  \f[  p^{jet'}_{T} / p^{true}_{T}\f],
//!  where \f$ p^{jet'}_{T} \f$ is the jet's pt corrected by
//!  the JetMET L5 JEC.
//!  Implements \p Function.
// ----------------------------------------------------------------   
double ControlPlotsFunction::jetTruthEventResponseL5Corrected(const Event * evt) const {
  const JetTruthEvent * jte = static_cast<const JetTruthEvent*>(evt);
  Jet * jet = static_cast<Jet*>(jte->mess());
  //  std::cout << "L1:" << jet->corFactors().getL1() << " l2l3:" << jet->corFactors().getL2L3() << " l2:" << jet->corFactors().getL2() << " l3:" << jet->corFactors().getL3() << " L5:" << jet->corFactors().getL5() << std::endl;
  return jet->corFactors().getL5() * jet->pt() / jte->truth();
}


//!  \brief Returns run number uf the Event
//!
//!  The \p Event \p evt has to be of type \p TwoJetsPtBalanceEvent.
//!  Implements \p Function.
// ----------------------------------------------------------------   
double ControlPlotsFunction::twoJetsPtBalanceEventRunNumber(const Event *evt) const {
  const TwoJetsPtBalanceEvent * jte = static_cast<const TwoJetsPtBalanceEvent*>(evt);
  return jte->runNumber();
}


//!  \brief Returns #phi of the jet
//!
//!  The \p Event \p evt has to be of type \p TwoJetsPtBalanceEvent.
//!  Implements \p Function.
// ----------------------------------------------------------------   
double ControlPlotsFunction::twoJetsPtBalanceEventJetPhi(const Event *evt) const {
  const TwoJetsPtBalanceEvent * jte = static_cast<const TwoJetsPtBalanceEvent*>(evt);
  return jte->getJet1()->phi();
}

//!  \brief Returns #eta of the jet
//!
//!  The \p Event \p evt has to be of type \p TwoJetsPtBalanceEvent.
//!  Implements \p Function.
// ----------------------------------------------------------------   
double ControlPlotsFunction::twoJetsPtBalanceEventJetEta(const Event *evt) const {
  const TwoJetsPtBalanceEvent * jte = static_cast<const TwoJetsPtBalanceEvent*>(evt);
  return jte->getJet1()->eta();
}

//!  \brief Returns |#eta| of the jet
//!
//!  The \p Event \p evt has to be of type \p TwoJetsPtBalanceEvent.
//!  Implements \p Function.
// ----------------------------------------------------------------   
double ControlPlotsFunction::twoJetsPtBalanceEventJetAbsEta(const Event *evt) const {
  const TwoJetsPtBalanceEvent * jte = static_cast<const TwoJetsPtBalanceEvent*>(evt);
  return std::abs(jte->getJet1()->eta());
}

//!  \brief Returns #eta of the jet2
//!
//!  The \p Event \p evt has to be of type \p TwoJetsPtBalanceEvent.
//!  Implements \p Function.
// ----------------------------------------------------------------   
double ControlPlotsFunction::twoJetsPtBalanceEventJet2Eta(const Event *evt) const {
  const TwoJetsPtBalanceEvent * jte = static_cast<const TwoJetsPtBalanceEvent*>(evt);
  return jte->getJet2()->eta();
}

//!  \brief Returns |#eta| of the jet2
//!
//!  The \p Event \p evt has to be of type \p TwoJetsPtBalanceEvent.
//!  Implements \p Function.
// ----------------------------------------------------------------   
double ControlPlotsFunction::twoJetsPtBalanceEventJet2AbsEta(const Event *evt) const {
  const TwoJetsPtBalanceEvent * jte = static_cast<const TwoJetsPtBalanceEvent*>(evt);
  return std::abs(jte->getJet2()->eta());
}


//!  \brief Returns absolute value of MET 
//!
//!  The \p Event \p evt has to be of type \p TwoJetsPtBalanceEvent.
//!  Implements \p Function.
// ----------------------------------------------------------------   
double ControlPlotsFunction::twoJetsPtBalanceEventMET(const Event *evt) const {
  const TwoJetsPtBalanceEvent * tjpbe = static_cast<const TwoJetsPtBalanceEvent*>(evt);
  return tjpbe->MET();
}


//!  \brief Returns absolute value of MET (with residual corrections)
//!
//!  The \p Event \p evt has to be of type \p TwoJetsPtBalanceEvent.
//!  Implements \p Function.
// ----------------------------------------------------------------   
double ControlPlotsFunction::twoJetsPtBalanceEventMETT1(const Event *evt) const {
  const TwoJetsPtBalanceEvent * tjpbe = static_cast<const TwoJetsPtBalanceEvent*>(evt);
  return tjpbe->METT1();
}
//!  \brief Returns absolute value of MET (with residual corrections)
//!
//!  The \p Event \p evt has to be of type \p TwoJetsPtBalanceEvent.
//!  Implements \p Function.
// ----------------------------------------------------------------   
double ControlPlotsFunction::twoJetsPtBalanceEventMETT1Res(const Event *evt) const {
  const TwoJetsPtBalanceEvent * tjpbe = static_cast<const TwoJetsPtBalanceEvent*>(evt);
  return tjpbe->METT1Res();
}

//!  \brief Returns absolute value of MET (with residual corrections)
//!
//!  The \p Event \p evt has to be of type \p TwoJetsPtBalanceEvent.
//!  Implements \p Function.
// ----------------------------------------------------------------   
double ControlPlotsFunction::twoJetsPtBalanceEventMETT2(const Event *evt) const {
  const TwoJetsPtBalanceEvent * tjpbe = static_cast<const TwoJetsPtBalanceEvent*>(evt);
  return tjpbe->METT2();
}
//!  \brief Returns absolute value of MET (with residual corrections)
//!
//!  The \p Event \p evt has to be of type \p TwoJetsPtBalanceEvent.
//!  Implements \p Function.
// ----------------------------------------------------------------   
double ControlPlotsFunction::twoJetsPtBalanceEventMETT2Res(const Event *evt) const {
  const TwoJetsPtBalanceEvent * tjpbe = static_cast<const TwoJetsPtBalanceEvent*>(evt);
  return tjpbe->METT2Res();
}

//!  \brief Returns value of MET projected on second jet axis 
//!
//!  The \p Event \p evt has to be of type \p TwoJetsPtBalanceEvent.
//!  Implements \p Function.
// ----------------------------------------------------------------   
double ControlPlotsFunction::twoJetsPtBalanceEventMETProj(const Event *evt) const {
  const TwoJetsPtBalanceEvent * tjpbe = static_cast<const TwoJetsPtBalanceEvent*>(evt);
  Jet * jet2 = tjpbe->getJet2();
  return  TMath::Abs(tjpbe->MET()*cos(deltaPhi(tjpbe->METphi(), jet2->phi())));
}

//!  \brief Returns value of MET (with residual corrections) projected on second jet axis 
//!
//!  The \p Event \p evt has to be of type \p TwoJetsPtBalanceEvent.
//!  Implements \p Function.
// ----------------------------------------------------------------   
double ControlPlotsFunction::twoJetsPtBalanceEventMETT1Proj(const Event *evt) const {
  const TwoJetsPtBalanceEvent * tjpbe = static_cast<const TwoJetsPtBalanceEvent*>(evt);
  Jet * jet2 = tjpbe->getJet2();
  return  TMath::Abs(tjpbe->METT1()*cos(deltaPhi(tjpbe->METT1phi(), jet2->phi())));
}

//!  \brief Returns value of MET (with residual corrections) projected on second jet axis 
//!
//!  The \p Event \p evt has to be of type \p TwoJetsPtBalanceEvent.
//!  Implements \p Function.
// ----------------------------------------------------------------   
double ControlPlotsFunction::twoJetsPtBalanceEventMETT1ResProj(const Event *evt) const {
  const TwoJetsPtBalanceEvent * tjpbe = static_cast<const TwoJetsPtBalanceEvent*>(evt);
  Jet * jet2 = tjpbe->getJet2();
  return  TMath::Abs(tjpbe->METT1Res()*cos(deltaPhi(tjpbe->METT1Resphi(), jet2->phi())));
}

//!  \brief Returns absolute value of METPhi 
//!
//!  The \p Event \p evt has to be of type \p TwoJetsPtBalanceEvent.
//!  Implements \p Function.
// ----------------------------------------------------------------   
double ControlPlotsFunction::twoJetsPtBalanceEventMETPhi(const Event *evt) const {
  const TwoJetsPtBalanceEvent * tjpbe = static_cast<const TwoJetsPtBalanceEvent*>(evt);
  return tjpbe->METphi();
}


//!  \brief Returns number of reconstructed vertices 
//!
//!  The \p Event \p evt has to be of type \p TwoJetsPtBalanceEvent.
//!  Implements \p Function.
// ----------------------------------------------------------------   
double ControlPlotsFunction::twoJetsPtBalanceEventVtxN(const Event *evt) const {
  const TwoJetsPtBalanceEvent * tjpbe = static_cast<const TwoJetsPtBalanceEvent*>(evt);
  return tjpbe->nVtx();
}


//!  \brief Returns number of MC PU vertices 
//!
//!  The \p Event \p evt has to be of type \p TwoJetsPtBalanceEvent.
//!  Implements \p Function.
// ----------------------------------------------------------------   
double ControlPlotsFunction::twoJetsPtBalanceEventMCNPUVtx(const Event *evt) const {
  const TwoJetsPtBalanceEvent * tjpbe = static_cast<const TwoJetsPtBalanceEvent*>(evt);
  return tjpbe->nPU();
}

//!  \brief Returns number of true MC PU vertices 
//!
//!  The \p Event \p evt has to be of type \p TwoJetsPtBalanceEvent.
//!  Implements \p Function.
// ----------------------------------------------------------------   
double ControlPlotsFunction::twoJetsPtBalanceEventMCNPUTruth(const Event *evt) const {
  const TwoJetsPtBalanceEvent * tjpbe = static_cast<const TwoJetsPtBalanceEvent*>(evt);
  return tjpbe->nPUTruth();
}


//!  \brief Returns the PF Charged hadron fraction (of jet1)
//!
//!  The \p Event \p evt has to be of type \p TwoJetsPtBalanceEvent.
//!  PF Charged hadron fraction (of jet1) is returned
//!  PF_CH_Fraction
//!  
// ----------------------------------------------------------------   
double ControlPlotsFunction::twoJetsPtBalanceEventPF_CH_Fraction(const Event * evt) const {
  const TwoJetsPtBalanceEvent * jte = static_cast<const TwoJetsPtBalanceEvent*>(evt);

  return jte->getJet1()->fCH();
}


//!  \brief Returns the PF neutral hadron fraction (of jet1)
//!
//!  The \p Event \p evt has to be of type \p TwoJetsPtBalanceEvent.
//!  PF neutral hadron fraction (of jet1) is returned
//!  PF_NH_Fraction
//!  
// ----------------------------------------------------------------   
double ControlPlotsFunction::twoJetsPtBalanceEventPF_NH_Fraction(const Event * evt) const {
  const TwoJetsPtBalanceEvent * jte = static_cast<const TwoJetsPtBalanceEvent*>(evt);

  return jte->getJet1()->fNH();
}


//!  \brief Returns the PF photon fraction (of jet1)
//!
//!  The \p Event \p evt has to be of type \p TwoJetsPtBalanceEvent.
//!  PF photon fraction (of jet1) is returned
//!  PF_PH_Fraction
//!  
// ----------------------------------------------------------------   
double ControlPlotsFunction::twoJetsPtBalanceEventPF_PH_Fraction(const Event * evt) const {
  const TwoJetsPtBalanceEvent * jte = static_cast<const TwoJetsPtBalanceEvent*>(evt);

  return jte->getJet1()->fPH();
}


//!  \brief Returns the PF electron fraction (of jet1)
//!
//!  The \p Event \p evt has to be of type \p TwoJetsPtBalanceEvent.
//!  PF electron fraction (of jet1) is returned
//!  PF_EL_Fraction
//!  
// ----------------------------------------------------------------   
double ControlPlotsFunction::twoJetsPtBalanceEventPF_EL_Fraction(const Event * evt) const {
  const TwoJetsPtBalanceEvent * jte = static_cast<const TwoJetsPtBalanceEvent*>(evt);

  return jte->getJet1()->fEL();
}

//!  \brief Returns the PF HFHad fraction (of jet1)
//!
//!  The \p Event \p evt has to be of type \p TwoJetsPtBalanceEvent.
//!  PF Charged hadron fraction (of jet1) is returned
//!  PF_CH_Fraction
//!  
// ----------------------------------------------------------------   
double ControlPlotsFunction::twoJetsPtBalanceEventPF_HFHad_Fraction(const Event * evt) const {
  const TwoJetsPtBalanceEvent * jte = static_cast<const TwoJetsPtBalanceEvent*>(evt);

  return jte->getJet1()->fHFHad();
}

//!  \brief Returns the PF HFEm fraction (of jet1)
//!
//!  The \p Event \p evt has to be of type \p TwoJetsPtBalanceEvent.
//!  PF Charged hadron fraction (of jet1) is returned
//!  PF_CH_Fraction
//!  
// ----------------------------------------------------------------   
double ControlPlotsFunction::twoJetsPtBalanceEventPF_HFEm_Fraction(const Event * evt) const {
  const TwoJetsPtBalanceEvent * jte = static_cast<const TwoJetsPtBalanceEvent*>(evt);

  return jte->getJet1()->fHFEm();
}



//!  \brief Returns the PF Charged hadron fraction (of jet1) multiplied 
//!  with the per-event L2L3-corrected response pt(jet1)/pt(jet2)
//!
//!  The \p Event \p evt has to be of type \p TwoJetsPtBalanceEvent.
//!  PF Charged hadron fraction (of jet1) multiplied 
//!  with the per-event L2L3-corrected response pt(jet1)/pt(jet2) is returned
//!  PF_CH_RespCorrFraction
//!  
// ----------------------------------------------------------------   
double ControlPlotsFunction::twoJetsPtBalanceEventPF_CH_RespCorrFraction(const Event * evt) const {
  const TwoJetsPtBalanceEvent * jte = static_cast<const TwoJetsPtBalanceEvent*>(evt);

  return jte->getJet1()->fCH()*(jte->getJet1()->pt() * jte->getJet1()->corFactors().getL2L3())/(jte->getJet2()->pt() * jte->getJet2()->corFactors().getL2L3());
}


//!  \brief Returns the PF neutral hadron fraction (of jet1) multiplied 
//!  with the per-event L2L3-corrected response pt(jet1)/pt(jet2)
//!
//!  The \p Event \p evt has to be of type \p TwoJetsPtBalanceEvent.
//!  PF neutral hadron fraction (of jet1) multiplied 
//!  with the per-event L2L3-corrected response pt(jet1)/pt(jet2) is returned
//!  PF_NH_RespCorrFraction
//!  
// ----------------------------------------------------------------   
double ControlPlotsFunction::twoJetsPtBalanceEventPF_NH_RespCorrFraction(const Event * evt) const {
  const TwoJetsPtBalanceEvent * jte = static_cast<const TwoJetsPtBalanceEvent*>(evt);

  return jte->getJet1()->fNH()*(jte->getJet1()->pt() * jte->getJet1()->corFactors().getL2L3())/(jte->getJet2()->pt() * jte->getJet2()->corFactors().getL2L3());
}


//!  \brief Returns the PF photon fraction (of jet1) multiplied 
//!  with the per-event L2L3-corrected response pt(jet1)/pt(jet2)
//!
//!  The \p Event \p evt has to be of type \p TwoJetsPtBalanceEvent.
//!  PF photon fraction (of jet1) multiplied 
//!  with the per-event L2L3-corrected response pt(jet1)/pt(jet2) is returned
//!  PF_PH_RespCorrFraction
//!  
// ----------------------------------------------------------------   
double ControlPlotsFunction::twoJetsPtBalanceEventPF_PH_RespCorrFraction(const Event * evt) const {
  const TwoJetsPtBalanceEvent * jte = static_cast<const TwoJetsPtBalanceEvent*>(evt);

  return jte->getJet1()->fPH()*(jte->getJet1()->pt() * jte->getJet1()->corFactors().getL2L3())/(jte->getJet2()->pt() * jte->getJet2()->corFactors().getL2L3());
}


//!  \brief Returns the PF electron fraction (of jet1) multiplied 
//!  with the per-event L2L3-corrected response pt(jet1)/pt(jet2)
//!
//!  The \p Event \p evt has to be of type \p TwoJetsPtBalanceEvent.
//!  PF electron fraction (of jet1) multiplied 
//!  with the per-event L2L3-corrected response pt(jet1)/pt(jet2) is returned
//!  PF_EL_RespCorrFraction
//!  
// ----------------------------------------------------------------   
double ControlPlotsFunction::twoJetsPtBalanceEventPF_EL_RespCorrFraction(const Event * evt) const {
  const TwoJetsPtBalanceEvent * jte = static_cast<const TwoJetsPtBalanceEvent*>(evt);

  return jte->getJet1()->fEL()*(jte->getJet1()->pt() * jte->getJet1()->corFactors().getL2L3())/(jte->getJet2()->pt() * jte->getJet2()->corFactors().getL2L3());
}

//!  \brief Returns the PF HFHad fraction (of jet1) multiplied 
//!  with the per-event L2L3-corrected response pt(jet1)/pt(jet2)
//!
//!  The \p Event \p evt has to be of type \p TwoJetsPtBalanceEvent.
//!  PF Charged hadron fraction (of jet1) multiplied 
//!  with the per-event L2L3-corrected response pt(jet1)/pt(jet2) is returned
//!  PF_CH_RespCorrFraction
//!  
// ----------------------------------------------------------------   
double ControlPlotsFunction::twoJetsPtBalanceEventPF_HFHad_RespCorrFraction(const Event * evt) const {
  const TwoJetsPtBalanceEvent * jte = static_cast<const TwoJetsPtBalanceEvent*>(evt);

  return jte->getJet1()->fHFHad()*(jte->getJet1()->pt() * jte->getJet1()->corFactors().getL2L3())/(jte->getJet2()->pt() * jte->getJet2()->corFactors().getL2L3());
}

//!  \brief Returns the PF HFEm fraction (of jet1) multiplied 
//!  with the per-event L2L3-corrected response pt(jet1)/pt(jet2)
//!
//!  The \p Event \p evt has to be of type \p TwoJetsPtBalanceEvent.
//!  PF Charged hadron fraction (of jet1) multiplied 
//!  with the per-event L2L3-corrected response pt(jet1)/pt(jet2) is returned
//!  PF_CH_RespCorrFraction
//!  
// ----------------------------------------------------------------   
double ControlPlotsFunction::twoJetsPtBalanceEventPF_HFEm_RespCorrFraction(const Event * evt) const {
  const TwoJetsPtBalanceEvent * jte = static_cast<const TwoJetsPtBalanceEvent*>(evt);

  return jte->getJet1()->fHFEm()*(jte->getJet1()->pt() * jte->getJet1()->corFactors().getL2L3())/(jte->getJet2()->pt() * jte->getJet2()->corFactors().getL2L3());
}




//!  \brief Returns flavor of jet
//!
//!  The \p Event \p evt has to be of type \p TwoJetsPtBalanceEvent.
//!  Implements \p Function.
// ----------------------------------------------------------------   
double ControlPlotsFunction::twoJetsPtBalanceEventJetFlavor(const Event *evt) const {
  const TwoJetsPtBalanceEvent * jte = static_cast<const TwoJetsPtBalanceEvent*>(evt);
  return jte->getJet1()->flavor();
}





//!  \brief Returns fraction of momentum of third jet projected to dijet axis 
//!
//!  The \p Event \p evt has to be of type \p TwoJetsPtBalanceEvent.
//!  Implements \p Function.
// ----------------------------------------------------------------   
double ControlPlotsFunction::twoJetsPtBalanceEventThirdJetFraction(const Event *evt) const {
  const TwoJetsPtBalanceEvent * tjpbe = static_cast<const TwoJetsPtBalanceEvent*>(evt);
  return tjpbe->relPtJet3ProjectionCorrL2L3();
}

//!  \brief Returns fraction of momentum of third jet and average dijet momentum  
//!
//!  The \p Event \p evt has to be of type \p TwoJetsPtBalanceEvent.
//!  Implements \p Function.
// ----------------------------------------------------------------   
double ControlPlotsFunction::twoJetsPtBalanceEventThirdJetFractionPlain(const Event *evt) const {
  const TwoJetsPtBalanceEvent * tjpbe = static_cast<const TwoJetsPtBalanceEvent*>(evt);
  return tjpbe->relPtJet3CorrL2L3();
}

//!  \brief Returns momentum of third jet  
//!
//!  The \p Event \p evt has to be of type \p TwoJetsPtBalanceEvent.
//!  Implements \p Function.
// ----------------------------------------------------------------   
double ControlPlotsFunction::twoJetsPtBalanceEventThirdJetPt(const Event *evt) const {
  const TwoJetsPtBalanceEvent * tjpbe = static_cast<const TwoJetsPtBalanceEvent*>(evt);
  if (! tjpbe->hasJet3()) return 0;
  return tjpbe->getJet3()->corFactors().getL2L3() * tjpbe->getJet3()->pt();
}

//!  \brief Returns mean p_{T} of the first two jets
//!
//!  The \p Event \p evt has to be of type \p TwoJetsPtBalanceEvent.
//!  Implements \p Function.
// ----------------------------------------------------------------   
double ControlPlotsFunction::twoJetsPtBalanceEventMeanPt(const Event *evt) const {
  const TwoJetsPtBalanceEvent * jte = static_cast<const TwoJetsPtBalanceEvent*>(evt);
  return jte->ptDijetCorrL2L3();
}


//!  \brief Returns p_{T} of the jet
//!
//!  The \p Event \p evt has to be of type \p TwoJetsPtBalanceEvent.
//!  Implements \p Function.
// ----------------------------------------------------------------   
double ControlPlotsFunction::twoJetsPtBalanceEventJetPt(const Event *evt) const {
  const TwoJetsPtBalanceEvent* jte = static_cast<const TwoJetsPtBalanceEvent*>(evt);
  return jte->getJet1()->pt();
}


//!  \brief Returns p_{T} of the jet (L2L3-corrected)
//!
//!  The \p Event \p evt has to be of type \p TwoJetsPtBalanceEvent.
//!  Implements \p Function.
// ----------------------------------------------------------------   
double ControlPlotsFunction::twoJetsPtBalanceEventJetPtL2L3Corrected(const Event *evt) const {
  const TwoJetsPtBalanceEvent* jte = static_cast<const TwoJetsPtBalanceEvent*>(evt);
  return jte->getJet1()->pt() * jte->getJet1()->corFactors().getL2L3() ;
}

//!  \brief Returns p_{T} of the jet (L2L3Res-corrected)
//!
//!  The \p Event \p evt has to be of type \p TwoJetsPtBalanceEvent.
//!  Implements \p Function.
// ----------------------------------------------------------------   
double ControlPlotsFunction::twoJetsPtBalanceEventJetPtL2L3ResCorrected(const Event *evt) const {
  const TwoJetsPtBalanceEvent* jte = static_cast<const TwoJetsPtBalanceEvent*>(evt);
  return jte->getJet1()->pt() * jte->getJet1()->corFactors().getL2L3Res() ;
}

//!  \brief Returns p_{T} of the second jet
//!
//!  The \p Event \p evt has to be of type \p TwoJetsPtBalanceEvent.
//!  Implements \p Function.
// ----------------------------------------------------------------   
double ControlPlotsFunction::twoJetsPtBalanceEventJet2Pt(const Event *evt) const {
  const TwoJetsPtBalanceEvent* jte = static_cast<const TwoJetsPtBalanceEvent*>(evt);
  return jte->getJet2()->pt();
}

//!  \brief Returns p_{T} of the second jet (L2L3-corrected)
//!
//!  The \p Event \p evt has to be of type \p TwoJetsPtBalanceEvent.
//!  Implements \p Function.
// ----------------------------------------------------------------   
double ControlPlotsFunction::twoJetsPtBalanceEventJet2PtL2L3Corrected(const Event *evt) const {
  const TwoJetsPtBalanceEvent* jte = static_cast<const TwoJetsPtBalanceEvent*>(evt);
  return jte->getJet2()->pt() * jte->getJet2()->corFactors().getL2L3() ;
}

//!  \brief Returns p_{T} of the second jet (L2L3Res-corrected)
//!
//!  The \p Event \p evt has to be of type \p TwoJetsPtBalanceEvent.
//!  Implements \p Function.
// ----------------------------------------------------------------   
double ControlPlotsFunction::twoJetsPtBalanceEventJet2PtL2L3ResCorrected(const Event *evt) const {
  const TwoJetsPtBalanceEvent* jte = static_cast<const TwoJetsPtBalanceEvent*>(evt);
  return jte->getJet2()->pt() * jte->getJet2()->corFactors().getL2L3Res() ;
}

//!  \brief Returns p_{T} of the leading of the two leading jets
//!
//!  The \p Event \p evt has to be of type \p TwoJetsPtBalanceEvent.
//!  Implements \p Function.
// ----------------------------------------------------------------   
double ControlPlotsFunction::twoJetsPtBalanceEventJetLeadPt(const Event *evt) const {
  const TwoJetsPtBalanceEvent* jte = static_cast<const TwoJetsPtBalanceEvent*>(evt);
  return twoJetsPtBalanceEventJet2Pt(jte) >= twoJetsPtBalanceEventJetPt(jte)? twoJetsPtBalanceEventJet2Pt(jte) : twoJetsPtBalanceEventJetPt(jte);
}

//!  \brief Returns p_{T} of the leading of the two leading  jet (L2L3-corrected)
//!
//!  The \p Event \p evt has to be of type \p TwoJetsPtBalanceEvent.
//!  Implements \p Function.
// ----------------------------------------------------------------   
double ControlPlotsFunction::twoJetsPtBalanceEventJetLeadPtL2L3Corrected(const Event *evt) const {
  const TwoJetsPtBalanceEvent* jte = static_cast<const TwoJetsPtBalanceEvent*>(evt);
  return twoJetsPtBalanceEventJet2PtL2L3Corrected(jte) >= twoJetsPtBalanceEventJetPtL2L3Corrected(jte)? twoJetsPtBalanceEventJet2PtL2L3Corrected(jte) : twoJetsPtBalanceEventJetPtL2L3Corrected(jte);
}

//!  \brief Returns p_{T} of the leading of the two leading  jet (L2L3Res-corrected)
//!
//!  The \p Event \p evt has to be of type \p TwoJetsPtBalanceEvent.
//!  Implements \p Function.
// ----------------------------------------------------------------   
double ControlPlotsFunction::twoJetsPtBalanceEventJetLeadPtL2L3ResCorrected(const Event *evt) const {
  const TwoJetsPtBalanceEvent* jte = static_cast<const TwoJetsPtBalanceEvent*>(evt);
  return twoJetsPtBalanceEventJet2PtL2L3ResCorrected(jte) >= twoJetsPtBalanceEventJetPtL2L3ResCorrected(jte)? twoJetsPtBalanceEventJet2PtL2L3ResCorrected(jte) : twoJetsPtBalanceEventJetPtL2L3ResCorrected(jte);
}
//////////////////////////////////////////////////////////////////////////////////
double ControlPlotsFunction::twoJetsPtBalanceEventJetLead2Pt(const Event *evt) const {
  const TwoJetsPtBalanceEvent* jte = static_cast<const TwoJetsPtBalanceEvent*>(evt);
  if(std::abs(jte->getJet1()->eta())<1.3)
   {
     if(jte->getJet1()->pt() > jte->getJet2()->pt())return jte->getJet1()->pt(); 	 
     else if(std::abs(jte->getJet2()->eta())<1.3)return jte->getJet2()->pt();
     else return jte->getJet1()->pt();
   }
  else return jte->getJet2()->pt();
}

//!  \brief Returns p_{T} of the leading of the two leading  barrel jets  (L2L3-corrected)
//!
//!  The \p Event \p evt has to be of type \p TwoJetsPtBalanceEvent.
//!  Implements \p Function.
// ----------------------------------------------------------------   
double ControlPlotsFunction::twoJetsPtBalanceEventJetLead2PtL2L3Corrected(const Event *evt) const {
  const TwoJetsPtBalanceEvent* jte = static_cast<const TwoJetsPtBalanceEvent*>(evt);
  if(std::abs(jte->getJet1()->eta())<1.3)
   {
     if(jte->getJet1()->pt() * jte->getJet1()->corFactors().getL2L3() > jte->getJet2()->pt() * jte->getJet2()->corFactors().getL2L3())return jte->getJet1()->pt() * jte->getJet1()->corFactors().getL2L3(); 	 
     else if(std::abs(jte->getJet2()->eta())<1.3)return jte->getJet2()->pt() * jte->getJet2()->corFactors().getL2L3();
     else return jte->getJet1()->pt() * jte->getJet1()->corFactors().getL2L3();
   }
  else return jte->getJet2()->pt() * jte->getJet2()->corFactors().getL2L3();
}

//!  \brief Returns p_{T} of the leading of the two leading barrel jets (L2L3Res-corrected)
//!
//!  The \p Event \p evt has to be of type \p TwoJetsPtBalanceEvent.
//!  Implements \p Function.
// ----------------------------------------------------------------   
double ControlPlotsFunction::twoJetsPtBalanceEventJetLead2PtL2L3ResCorrected(const Event *evt) const {
  const TwoJetsPtBalanceEvent* jte = static_cast<const TwoJetsPtBalanceEvent*>(evt);
  if(std::abs(jte->getJet1()->eta())<1.3)
   {
     if(jte->getJet1()->pt() * jte->getJet1()->corFactors().getL2L3Res() > jte->getJet2()->pt() *
     jte->getJet2()->corFactors().getL2L3Res())return jte->getJet1()->pt() * jte->getJet1()->corFactors().getL2L3Res(); 	 
     else if(std::abs(jte->getJet2()->eta())<1.3)return jte->getJet2()->pt() * jte->getJet2()->corFactors().getL2L3Res();
     else return jte->getJet1()->pt() * jte->getJet1()->corFactors().getL2L3Res();
   }
  else return jte->getJet2()->pt() * jte->getJet2()->corFactors().getL2L3Res();
}
////////////////////////////////////////////////////////////////////////////////
//!  \brief Returns ECal fraction of the jet
//!
//!  The \p Event \p evt has to be of type \p TwoJetsPtBalanceEvent.
//!  Implements \p Function.
// ----------------------------------------------------------------   
double ControlPlotsFunction::twoJetsPtBalanceEventJetEMF(const Event *evt) const {
  const TwoJetsPtBalanceEvent * jte = static_cast<const TwoJetsPtBalanceEvent*>(evt);
  return jte->getJet1()->emf();
}


//!  \brief Returns the #phi #phi moment of the jet
//!
//!  The \p Event \p evt has to be of type \p TwoJetsPtBalanceEvent.
//!  Implements \p Function.
// ----------------------------------------------------------------   
double ControlPlotsFunction::twoJetsPtBalanceEventJetMomentPhiPhi(const Event *evt) const {
  const TwoJetsPtBalanceEvent * jte = static_cast<const TwoJetsPtBalanceEvent*>(evt);
  return jte->getJet1()->momentPhiPhi();
}


//!  \brief Returns the #eta #eta moment of the jet
//!
//!  The \p Event \p evt has to be of type \p TwoJetsPtBalanceEvent.
//!  Implements \p Function.
// ----------------------------------------------------------------   
double ControlPlotsFunction::twoJetsPtBalanceEventJetMomentEtaEta(const Event *evt) const {
  const TwoJetsPtBalanceEvent * jte = static_cast<const TwoJetsPtBalanceEvent*>(evt);
  return jte->getJet1()->momentEtaEta();
}

//!  \brief Returns the mean of the #phi #phi and #eta #eta moments of the jet
//!
//!  The \p Event \p evt has to be of type \p TwoJetsPtBalanceEvent.
//!  Implements \p Function.
// ----------------------------------------------------------------   
double ControlPlotsFunction::twoJetsPtBalanceEventJetMeanMoment(const Event *evt) const {
  const TwoJetsPtBalanceEvent * jte = static_cast<const TwoJetsPtBalanceEvent*>(evt);
  return 0.5 * (jte->getJet1()->momentEtaEta() + jte->getJet1()->momentPhiPhi());
}

//!  \brief Returns the jet deltaphi
//!
//!  The \p Event \p evt has to be of type \p TwoJetsPtBalanceEvent.
//!  The deltaphi is defined as
//!  \f[  p^{jet}_{T} / p^{true}_{T}\f].
//!  Implements \p Function.
// ----------------------------------------------------------------   
double ControlPlotsFunction::twoJetsPtBalanceEventDeltaPhi(const Event * evt) const {
  const TwoJetsPtBalanceEvent * jte = static_cast<const TwoJetsPtBalanceEvent*>(evt);
  Jet * jet1 = jte->getJet1();
  Jet * jet2 = jte->getJet2();
  return (std::abs(TVector2::Phi_mpi_pi(jet1->phi() - jet2->phi())) );
}

//!  \brief Returns the DeltaR to the closest reco jet (of jet1)
//!
//!  The \p Event \p evt has to be of type \p TwoJetsPtBalanceEvent.
//!  DeltaR (of jet1) to closest other reco jet is returned
//!  
//!  
// ----------------------------------------------------------------   
double ControlPlotsFunction::twoJetsPtBalanceEventClosestJetdR(const Event * evt) const {
  const TwoJetsPtBalanceEvent * jte = static_cast<const TwoJetsPtBalanceEvent*>(evt);
  Jet * jet1 = jte->getJet1();
  assert(jet1->closestJetdR()>0);
  //DEBUG:  if(jet1->closestJetdR()<0.1)std::cout << jet1->closestJetdR() << " jet1->pt() " << jet1->pt() << std::endl;
  return jet1->closestJetdR();
}





//!  \brief Returns the jet asymmetry
//!
//!  The \p Event \p evt has to be of type \p TwoJetsPtBalanceEvent.
//!  The asymmetry is defined as
//!  \f[  p^{jet}_{T} / p^{true}_{T}\f].
//!  Implements \p Function.
// ----------------------------------------------------------------   
double ControlPlotsFunction::twoJetsPtBalanceEventAsymmetry(const Event * evt) const {
  const TwoJetsPtBalanceEvent * jte = static_cast<const TwoJetsPtBalanceEvent*>(evt);
  Jet * jet1 = jte->getJet1();
  Jet * jet2 = jte->getJet2();
  return (jet1->pt()-jet2->pt())/(jet1->pt()+jet2->pt());
}



//!  \brief Returns the corrected jet asymmetry
//!
//!  The \p Event \p evt has to be of type \p TwoJetsPtBalanceEvent.
//!  The asymmetry is defined as
//!  \f[  p^{jet'}_{T} / p^{true}_{T}\f],
//!  where \f$ p^{jet'}_{T} \f$ is the jet's pt corrected by
//!  the Kalibri JEC.
//!  Implements \p Function.
// ----------------------------------------------------------------   
double ControlPlotsFunction::twoJetsPtBalanceEventAsymmetryKalibriCorrected(const Event * evt) const {
  const TwoJetsPtBalanceEvent * jte = static_cast<const TwoJetsPtBalanceEvent*>(evt);
  Jet * jet1 = jte->getJet1();
  Jet * jet2 = jte->getJet2();
  return (jet1->correctedEt()-jet2->correctedEt())/(jet1->correctedEt()+jet2->correctedEt());
}



//!  \brief Returns the corrected jet asymmetry
//!
//!  The \p Event \p evt has to be of type \p TwoJetsPtBalanceEvent.
//!  The asymmetry is defined as
//!  \f[  p^{jet'}_{T} / p^{true}_{T}\f],
//!  where \f$ p^{jet'}_{T} \f$ is the jet's pt corrected by
//!  the JetMET L2L3 JEC.
//!  Implements \p Function.
// ----------------------------------------------------------------   
double ControlPlotsFunction::twoJetsPtBalanceEventAsymmetryL2L3Corrected(const Event * evt) const {
  const TwoJetsPtBalanceEvent * jte = static_cast<const TwoJetsPtBalanceEvent*>(evt);

  Jet * jet1 = jte->getJet1();
  Jet * jet2 = jte->getJet2();
  return (jet1->corFactors().getL2L3() * jet1->pt() - jet2->corFactors().getL2L3() * jet2->pt())/(jet1->corFactors().getL2L3() * jet1->pt()+jet2->corFactors().getL2L3() * jet2->pt());
}

//!  \brief Returns the corrected jet asymmetry
//!
//!  The \p Event \p evt has to be of type \p TwoJetsPtBalanceEvent.
//!  The asymmetry is defined as
//!  \f[  p^{jet'}_{T} / p^{true}_{T}\f],
//!  where \f$ p^{jet'}_{T} \f$ is the jet's pt corrected by
//!  the JetMET L2L3Res JEC.
//!  Implements \p Function.
// ----------------------------------------------------------------   
double ControlPlotsFunction::twoJetsPtBalanceEventAsymmetryL2L3ResCorrected(const Event * evt) const {
  const TwoJetsPtBalanceEvent * jte = static_cast<const TwoJetsPtBalanceEvent*>(evt);

  Jet * jet1 = jte->getJet1();
  Jet * jet2 = jte->getJet2();
  return (jet1->corFactors().getL2L3Res() * jet1->pt() - jet2->corFactors().getL2L3Res() * jet2->pt())/(jet1->corFactors().getL2L3Res() * jet1->pt()+jet2->corFactors().getL2L3Res() * jet2->pt());
}


//!  \brief Returns the corrected jet asymmetry
//!
//!  The \p Event \p evt has to be of type \p TwoJetsPtBalanceEvent.
//!  The asymmetry is defined as
//!  \f[  p^{jet'}_{T} / p^{true}_{T}\f],
//!  where \f$ p^{jet'}_{T} \f$ is the jet's pt corrected by
//!  the JetMET L2L3L4 JEC.
//!  Implements \p Function.
// ----------------------------------------------------------------   
double ControlPlotsFunction::twoJetsPtBalanceEventAsymmetryL2L3L4Corrected(const Event * evt) const {
  const TwoJetsPtBalanceEvent * jte = static_cast<const TwoJetsPtBalanceEvent*>(evt);

  Jet * jet1 = jte->getJet1();
  Jet * jet2 = jte->getJet2();
  return (jet1->corFactors().getL2L3L4() * jet1->pt() - jet2->corFactors().getL2L3L4() * jet2->pt())/(jet1->corFactors().getL2L3L4() * jet1->pt()+jet2->corFactors().getL2L3L4() * jet2->pt());
}


//!  \brief Returns the corrected jet asymmetry
//!
//!  The \p Event \p evt has to be of type \p TwoJetsPtBalanceEvent.
//!  The asymmetry is defined as
//!  \f[  p^{jet'}_{T} / p^{true}_{T}\f],
//!  where \f$ p^{jet'}_{T} \f$ is the jet's pt corrected by
//!  the JetMET L2L3ResL4 JEC.
//!  Implements \p Function.
// ----------------------------------------------------------------   
double ControlPlotsFunction::twoJetsPtBalanceEventAsymmetryL2L3ResL4Corrected(const Event * evt) const {
  const TwoJetsPtBalanceEvent * jte = static_cast<const TwoJetsPtBalanceEvent*>(evt);

  Jet * jet1 = jte->getJet1();
  Jet * jet2 = jte->getJet2();
  return ( jet1->corFactors().getL2L3Res() * jet1->corFactors().getL4() * jet1->pt() - jet2->corFactors().getL2L3Res() * jet2->corFactors().getL4() * jet2->pt())/( jet1->corFactors().getL2L3Res() * jet1->corFactors().getL4() * jet1->pt()+ jet2->corFactors().getL2L3Res() * jet2->corFactors().getL4() * jet2->pt());
}


//!  \brief Returns the jet genasymmetry
//!
//!  The \p Event \p evt has to be of type \p TwoJetsPtBalanceEvent.
//!  The genasymmetry is defined as
//!  \f[  (p^{true}_{T,1}-p^{true}_{T,2}) / (p^{true}_{T,1}+p^{true}_{T,2})\f].
//!  Implements \p Function.
// ----------------------------------------------------------------   
double ControlPlotsFunction::twoJetsPtBalanceEventGenAsymmetry(const Event * evt) const {
  const TwoJetsPtBalanceEvent * jte = static_cast<const TwoJetsPtBalanceEvent*>(evt);
  Jet * jet1 = jte->getJet1();
  Jet * jet2 = jte->getJet2();
  return 2*(jet1->genPt()-jet2->genPt())/(jet1->genPt()+jet2->genPt());
}

double ControlPlotsFunction::twoJetsPtBalanceEventGenLeadAsymmetry(const Event * evt) const {
  const TwoJetsPtBalanceEvent * jte = static_cast<const TwoJetsPtBalanceEvent*>(evt);
  Jet * jet1 = jte->getJet1();
  Jet * jet2 = jte->getJet2();
  Jet * jetL;
  if(jet1->genPt()>jet2->genPt()) jetL=jet1;
  else jetL=jet2;
  return (jet1->genPt()-jet2->genPt())/jetL->genPt();
}

double ControlPlotsFunction::twoJetsPtBalanceEventGenBarrAsymmetry(const Event * evt) const {
  const TwoJetsPtBalanceEvent * jte = static_cast<const TwoJetsPtBalanceEvent*>(evt);
  Jet * jet1 = jte->getJet1();
  Jet * jet2 = jte->getJet2();
  return (jet1->genPt()-jet2->genPt())/jet2->genPt();
}

//!  \brief Returns the jet B
//!
//!  The \p Event \p evt has to be of type \p TwoJetsPtBalanceEvent.
//!  The B is defined as
//!  \f[  p^{jet}_{T} / p^{true}_{T}\f].
//!  Implements \p Function.
// ----------------------------------------------------------------   
double ControlPlotsFunction::twoJetsPtBalanceEventB(const Event * evt) const {
  const TwoJetsPtBalanceEvent * jte = static_cast<const TwoJetsPtBalanceEvent*>(evt);
  Jet * jet1 = jte->getJet1();
  Jet * jet2 = jte->getJet2();
  return (jet1->pt()-jet2->pt())/(jet1->pt()+jet2->pt())*2;
}



//!  \brief Returns the corrected jet B
//!
//!  The \p Event \p evt has to be of type \p TwoJetsPtBalanceEvent.
//!  The B is defined as
//!  \f[  p^{jet'}_{T} / p^{true}_{T}\f],
//!  where \f$ p^{jet'}_{T} \f$ is the jet's pt corrected by
//!  the Kalibri JEC.
//!  Implements \p Function.
// ----------------------------------------------------------------   
double ControlPlotsFunction::twoJetsPtBalanceEventBKalibriCorrected(const Event * evt) const {
  const TwoJetsPtBalanceEvent * jte = static_cast<const TwoJetsPtBalanceEvent*>(evt);
  Jet * jet1 = jte->getJet1();
  Jet * jet2 = jte->getJet2();
  return (jet1->correctedEt()-jet2->correctedEt())/(jet1->correctedEt()+jet2->correctedEt())*2;
}



//!  \brief Returns the corrected jet B
//!
//!  The \p Event \p evt has to be of type \p TwoJetsPtBalanceEvent.
//!  The B is defined as
//!  \f[  p^{jet'}_{T} / p^{true}_{T}\f],
//!  where \f$ p^{jet'}_{T} \f$ is the jet's pt corrected by
//!  the JetMET L2L3L4 JEC.
//!  Implements \p Function.
// ----------------------------------------------------------------   
double ControlPlotsFunction::twoJetsPtBalanceEventBL2L3L4Corrected(const Event * evt) const {
  const TwoJetsPtBalanceEvent * jte = static_cast<const TwoJetsPtBalanceEvent*>(evt);

  Jet * jet1 = jte->getJet1();
  Jet * jet2 = jte->getJet2();
  return (jet1->corFactors().getL2L3L4() * jet1->pt() - jet2->corFactors().getL2L3L4() * jet2->pt())/(jet1->corFactors().getL2L3L4() * jet1->pt()+jet2->corFactors().getL2L3L4() * jet2->pt())*2;
}


//!  \brief Returns the corrected jet B
//!
//!  The \p Event \p evt has to be of type \p TwoJetsPtBalanceEvent.
//!  The B is defined as
//!  \f[  p^{jet'}_{T} / p^{true}_{T}\f],
//!  where \f$ p^{jet'}_{T} \f$ is the jet's pt corrected by
//!  the JetMET L2L3ResL4 JEC.
//!  Implements \p Function.
// ----------------------------------------------------------------   
double ControlPlotsFunction::twoJetsPtBalanceEventBL2L3ResL4Corrected(const Event * evt) const {
  const TwoJetsPtBalanceEvent * jte = static_cast<const TwoJetsPtBalanceEvent*>(evt);

  Jet * jet1 = jte->getJet1();
  Jet * jet2 = jte->getJet2();
  return (jet1->corFactors().getL2L3Res() * jet1->corFactors().getL4() * jet1->pt() - jet2->corFactors().getL2L3Res() * jet2->corFactors().getL4() * jet2->pt())/(jet1->corFactors().getL2L3Res() * jet1->corFactors().getL4() * jet1->pt()+jet2->corFactors().getL2L3Res() * jet2->corFactors().getL4() * jet2->pt())*2;
}

//!  \brief Returns the corrected jet B
//!
//!  The \p Event \p evt has to be of type \p TwoJetsPtBalanceEvent.
//!  The B is defined as
//!  \f[  p^{jet'}_{T} / p^{true}_{T}\f],
//!  where \f$ p^{jet'}_{T} \f$ is the jet's pt corrected by
//!  the JetMET L2L3 JEC.
//!  Implements \p Function.
// ----------------------------------------------------------------   
double ControlPlotsFunction::twoJetsPtBalanceEventBL2L3Corrected(const Event * evt) const {
  const TwoJetsPtBalanceEvent * jte = static_cast<const TwoJetsPtBalanceEvent*>(evt);

  Jet * jet1 = jte->getJet1();
  Jet * jet2 = jte->getJet2();
  return (jet1->corFactors().getL2L3() * jet1->pt() - jet2->corFactors().getL2L3() * jet2->pt())/(jet1->corFactors().getL2L3() * jet1->pt()+jet2->corFactors().getL2L3() * jet2->pt())*2;
}

//!  \brief Returns the corrected jet B
//!
//!  The \p Event \p evt has to be of type \p TwoJetsPtBalanceEvent.
//!  The B is defined as
//!  \f[  p^{jet'}_{T} / p^{true}_{T}\f],
//!  where \f$ p^{jet'}_{T} \f$ is the jet's pt corrected by
//!  the JetMET L2L3Res JEC.
//!  Implements \p Function.
// ----------------------------------------------------------------   
double ControlPlotsFunction::twoJetsPtBalanceEventBL2L3ResCorrected(const Event * evt) const {
  const TwoJetsPtBalanceEvent * jte = static_cast<const TwoJetsPtBalanceEvent*>(evt);

  Jet * jet1 = jte->getJet1();
  Jet * jet2 = jte->getJet2();
  return (jet1->corFactors().getL2L3Res() * jet1->pt() - jet2->corFactors().getL2L3Res() * jet2->pt())/(jet1->corFactors().getL2L3Res() * jet1->pt()+jet2->corFactors().getL2L3Res() * jet2->pt())*2;
}





//!  \brief Returns the jet MPF response
//!
//!  The \p Event \p evt has to be of type \p TwoJetsPtBalanceEvent.
//!  The MPF response is defined as 
//!  \f[ mpf = 1 + met*cos(delta_phi(metphi, jtphi[iref])) / p^{jet'}_{T} \f].
//!  Implements \p Function.
// ----------------------------------------------------------------   
double ControlPlotsFunction::twoJetsPtBalanceEventMPFResponse(const Event * evt) const {
  const TwoJetsPtBalanceEvent * jte = static_cast<const TwoJetsPtBalanceEvent*>(evt);
  Jet * jet2 = jte->getJet2();
  return (1 + jte->MET()*cos(deltaPhi<double>(jte->METphi(), jet2->phi())) / jet2->pt());
}

//!  \brief Returns the corrected jet MPF response
//!
//!  The \p Event \p evt has to be of type \p TwoJetsPtBalanceEvent.
//!  The MPF response is defined as
//!  \f[ mpf = 1 + met*cos(delta_phi(metphi, jtphi[iref])) / p^{jet'}_{T}\f],
//!  where \f$ p^{jet'}_{T} \f$ is the jet's pt corrected by
//!  the JetMET L2L3 JEC.
//!  Implements \p Function.
// ----------------------------------------------------------------   
double ControlPlotsFunction::twoJetsPtBalanceEventMPFResponseL2L3Corrected(const Event * evt) const {
  const TwoJetsPtBalanceEvent * jte = static_cast<const TwoJetsPtBalanceEvent*>(evt);
  Jet * jet2 = jte->getJet2();
  return (1 + jte->MET()*cos(deltaPhi<double>(jte->METphi(), jet2->phi())) / (jet2->corFactors().getL2L3() *jet2->pt()));
}

//!  \brief Returns the corrected jet MPF response
//!
//!  The \p Event \p evt has to be of type \p TwoJetsPtBalanceEvent.
//!  The MPF response is defined as
//!  \f[ mpf = 1 + met*cos(delta_phi(metphi, jtphi[iref])) / p^{jet'}_{T}\f],
//!  where \f$ p^{jet'}_{T} \f$ is the jet's pt corrected by
//!  the JetMET L2L3Res JEC.
//!  Implements \p Function.
// ----------------------------------------------------------------   
double ControlPlotsFunction::twoJetsPtBalanceEventMPFResponseL2L3ResCorrected(const Event * evt) const {
  const TwoJetsPtBalanceEvent * jte = static_cast<const TwoJetsPtBalanceEvent*>(evt);
  Jet * jet2 = jte->getJet2();
  return (1 + jte->MET()*cos(deltaPhi<double>(jte->METphi(), jet2->phi())) / (jet2->corFactors().getL2L3Res() *jet2->pt()));
}


//!  \brief Returns the jet MPFMETT1 response (with type-1 corrected MET)
//!
//!  The \p Event \p evt has to be of type \p TwoJetsPtBalanceEvent.
//!  The MPFMETT1 response is defined as 
//!  \f[ mpf = 1 + met*cos(delta_phi(metphi, jtphi[iref])) / p^{jet'}_{T} \f].
//!  Implements \p Function.
// ----------------------------------------------------------------   
double ControlPlotsFunction::twoJetsPtBalanceEventMPFMETT1Response(const Event * evt) const {
  const TwoJetsPtBalanceEvent * jte = static_cast<const TwoJetsPtBalanceEvent*>(evt);
  Jet * jet2 = jte->getJet2();
  return (1 + jte->METT1()*cos(deltaPhi<double>(jte->METT1phi(), jet2->phi())) / jet2->pt());
}



//!  \brief Returns the corrected jet MPFMETT1 response (with type-1 corrected MET)
//!
//!  The \p Event \p evt has to be of type \p TwoJetsPtBalanceEvent.
//!  The MPFMETT1 response is defined as
//!  \f[ mpf = 1 + met*cos(delta_phi(metphi, jtphi[iref])) / p^{jet'}_{T}\f],
//!  where \f$ p^{jet'}_{T} \f$ is the jet's pt corrected by
//!  the JetMET L2L3 JEC.
//!  Implements \p Function.
// ----------------------------------------------------------------   
double ControlPlotsFunction::twoJetsPtBalanceEventMPFMETT1ResponseL2L3Corrected(const Event * evt) const {
  const TwoJetsPtBalanceEvent * jte = static_cast<const TwoJetsPtBalanceEvent*>(evt);
  Jet * jet2 = jte->getJet2();
  return (1 + jte->METT1()*cos(deltaPhi<double>(jte->METT1phi(), jet2->phi())) / (jet2->corFactors().getL2L3() *jet2->pt()));
}

//!  \brief Returns the corrected jet MPFMETT1 response (with type-1 corrected MET)
//!
//!  The \p Event \p evt has to be of type \p TwoJetsPtBalanceEvent.
//!  The MPFMETT1 response is defined as
//!  \f[ mpf = 1 + met*cos(delta_phi(metphi, jtphi[iref])) / p^{jet'}_{T}\f],
//!  where \f$ p^{jet'}_{T} \f$ is the jet's pt corrected by
//!  the JetMET L2L3Res JEC.
//!  Implements \p Function.
// ----------------------------------------------------------------   
double ControlPlotsFunction::twoJetsPtBalanceEventMPFMETT1ResponseL2L3ResCorrected(const Event * evt) const {
  const TwoJetsPtBalanceEvent * jte = static_cast<const TwoJetsPtBalanceEvent*>(evt);
  Jet * jet2 = jte->getJet2();
  return (1 + jte->METT1Res()*cos(deltaPhi<double>(jte->METT1Resphi(), jet2->phi())) / (jet2->corFactors().getL2L3Res() *jet2->pt()));
}

//!  \brief Returns the corrected jet MPFMETT1 response (with type-1 corrected MET)
//!
//!  The \p Event \p evt has to be of type \p TwoJetsPtBalanceEvent.
//!  The MPFMETT1 response is defined as
//!  \f[ mpf = 1 + met*cos(delta_phi(metphi, jtphi[iref])) / p^{jet'}_{T}\f],
//!  where \f$ p^{jet'}_{T} \f$ is the jet's pt corrected by
//!  the JetMET L2L3 JEC.
//!  Implements \p Function.
// ----------------------------------------------------------------   

double ControlPlotsFunction::twoJetsPtBalanceEventMPFMETT2Response(const Event * evt) const {
  const TwoJetsPtBalanceEvent * jte = static_cast<const TwoJetsPtBalanceEvent*>(evt);
  Jet * jet2 = jte->getJet2();
  return (1 + jte->METT2()*cos(deltaPhi<double>(jte->METT2phi(), jet2->phi())) / jet2->pt());
}

double ControlPlotsFunction::twoJetsPtBalanceEventMPFMETT2ResponseL2L3Corrected(const Event * evt) const {
  const TwoJetsPtBalanceEvent * jte = static_cast<const TwoJetsPtBalanceEvent*>(evt);
  Jet * jet2 = jte->getJet2();
  return (1 + jte->METT2()*cos(deltaPhi<double>(jte->METT2phi(), jet2->phi())) / (jet2->corFactors().getL2L3() *jet2->pt()));
}

//!  \brief Returns the corrected jet MPFMETT1 response (with type-1 corrected MET)
//!
//!  The \p Event \p evt has to be of type \p TwoJetsPtBalanceEvent.
//!  The MPFMETT1 response is defined as
//!  \f[ mpf = 1 + met*cos(delta_phi(metphi, jtphi[iref])) / p^{jet'}_{T}\f],
//!  where \f$ p^{jet'}_{T} \f$ is the jet's pt corrected by
//!  the JetMET L2L3Res JEC.
//!  Implements \p Function.
// ----------------------------------------------------------------   
double ControlPlotsFunction::twoJetsPtBalanceEventMPFMETT2ResponseL2L3ResCorrected(const Event * evt) const {
  const TwoJetsPtBalanceEvent * jte = static_cast<const TwoJetsPtBalanceEvent*>(evt);
  Jet * jet2 = jte->getJet2();
  return (1 + jte->METT2Res()*cos(deltaPhi<double>(jte->METT2Resphi(), jet2->phi())) / (jet2->corFactors().getL2L3Res() *jet2->pt()));
}

//!  \brief Returns the jet1 response
//!
//!  The \p Event \p evt has to be of type \p TwoJetsPtBalanceEvent and "data" and MC must indeed be MC.
//!  The response is defined as
//!  \f[  p^{jet}_{T} / p^{true}_{T}\f].
//!  Implements \p Function.
// ----------------------------------------------------------------   
double ControlPlotsFunction::twoJetsPtBalanceEventMCTruthJet1Response(const Event * evt) const {
  const TwoJetsPtBalanceEvent * jte = static_cast<const TwoJetsPtBalanceEvent*>(evt);
  Jet * jet1 = jte->getJet1();

  return jet1->pt() / jet1->genPt();
}


//!  \brief Returns the L2L3-corrected jet1 response
//!
//!  The \p Event \p evt has to be of type \p TwoJetsPtBalanceEvent and "data" and MC must indeed be MC.
//!  The response is defined as
//!  \f[  p^{jet}_{T} / p^{true}_{T}\f].
//!  Implements \p Function.
// ----------------------------------------------------------------   
double ControlPlotsFunction::twoJetsPtBalanceEventMCTruthJet1L2L3Response(const Event * evt) const {
  const TwoJetsPtBalanceEvent * jte = static_cast<const TwoJetsPtBalanceEvent*>(evt);
  Jet * jet1 = jte->getJet1();

  return jet1->corFactors().getL2L3()*jet1->pt() / jet1->genPt();
}


