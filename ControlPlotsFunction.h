// $Id: ControlPlotsFunction.h,v 1.34 2013/03/14 13:16:03 kirschen Exp $

#ifndef CONTROLPLOTS_FUNCTION_H
#define CONTROLPLOTS_FUNCTION_H

#include "ControlPlotsConfig.h"

class Event;


//!  \brief Functions for the profile control plots
//!
//!  Stores the functions for the profile control plots. This class
//!  is used by \p ControlPlotsProfile and returns for a given \p Event
//!  the x value, the y value and the binning value for the profile
//!  plots. There are different y values depending on the applied
//!  correction.
//!
//!  For example, for a "Response vs GenJetPt" plot in bins of "Eta",
//!  this class stores the functions to get Response, GenJetPt and
//!  Eta from an \p Event.
//!
//!  \sa \p ControlPlotsConfig, \p ControlPlotsProfile  
//!
//!  \author Matthias Schroeder
//!  \date 2009/12/18
//!  $Id: ControlPlotsFunction.h,v 1.34 2013/03/14 13:16:03 kirschen Exp $
// ----------------------------------------------------------------   
class ControlPlotsFunction {
 public:
  //! The function's signature
  typedef double (ControlPlotsFunction::*Function)(const Event *evt) const;

  //! Constructor
  ControlPlotsFunction()
    : binFunc_(0), xFunc_(0), cutFunc_(0) {};

  //! Check if all functions are initialised
  bool isInit() const { return binFunc_ && xFunc_ && yFuncs_.size(); }

  //! Return iterator to first \p CorrectionType
  ControlPlotsConfig::CorrectionTypeIt correctionTypeBegin() const { return types_.begin(); }
  //! Return iterator to last \p CorrectionType
  ControlPlotsConfig::CorrectionTypeIt correctionTypeEnd() const { return types_.end(); }

  //! Interface to the profile: return the value of the binning quantity from \p evt
  double binValue(const Event * evt) const { return (this->*binFunc_)(evt); }
  //! Interface to the profile: return the value of the x quantity from \p evt
  double xValue(const Event * evt) const { return (this->*xFunc_)(evt); }  
  //! Interface to the profile: return the value of the cut quantity from \p evt
  double cutValue(const Event * evt) const { return cutFunc_ ? (this->*cutFunc_)(evt) : 0; }
  //! Interface to the profile: return the value of the y quantity from \p evt for the correction type \p type
  double yValue(const Event * evt, ControlPlotsConfig::CorrectionType type) const {
    return (this->*(yFuncs_.find(type)->second))(evt);
  }
  //! Set the function returning the binning value from an event
  void setBinFunction(Function func) { binFunc_ = func; }
  //! Set the function returning the x value from an event
  void setXFunction(Function func) { xFunc_ = func; } 
  //! Set the function returning the cut value from an event
  void setCutFunction(Function func) { cutFunc_ = func; }
  //! Set the functions returning the y value for different corrections from an event
  void addYFunction(ControlPlotsConfig::CorrectionType type, Function func);

  double jetTruthEventNPU(const Event *evt) const;
  double jetTruthEventRho(const Event *evt) const;
  double jetTruthEventJetEta(const Event *evt) const;
  double jetTruthEventJetAbsEta(const Event *evt) const;
  double jetTruthEventJetPt(const Event *evt) const;
  double jetTruthEventJetEMF(const Event *evt) const;
  double jetTruthEventJetMomentPhiPhi(const Event *evt) const;
  double jetTruthEventJetMomentEtaEta(const Event *evt) const;
  double jetTruthEventJetMeanMoment(const Event *evt) const;
  double jetTruthEventTruthPt(const Event *evt) const;
  double jetTruthEventJetFlavor(const Event *evt) const;
  double jetTruthEventResponse(const Event * evt) const;
  double jetTruthEventResponseKalibriCorrected(const Event * evt) const;
  double jetTruthEventResponseL2L3Corrected(const Event * evt) const;
  double jetTruthEventResponseL2L3ResCorrected(const Event * evt) const;
  double jetTruthEventResponseL2L3L4Corrected(const Event * evt) const;
  double jetTruthEventResponseL2L3ResL4Corrected(const Event * evt) const;
  double jetTruthEventResponseL1L2L3Corrected(const Event * evt) const;
  double jetTruthEventResponseL5Corrected(const Event * evt) const;

  double twoJetsPtBalanceEventRunNumber(const Event *evt) const; 
  double twoJetsPtBalanceEventJetPhi(const Event *evt) const; 
  double twoJetsPtBalanceEventJetEta(const Event *evt) const; 
  double twoJetsPtBalanceEventJetAbsEta(const Event *evt) const;
  double twoJetsPtBalanceEventJet2Eta(const Event *evt) const; 
  double twoJetsPtBalanceEventJet2AbsEta(const Event *evt) const;
  double twoJetsPtBalanceEventJetPt(const Event *evt) const;
  double twoJetsPtBalanceEventJetPtL2L3Corrected(const Event *evt) const;
  double twoJetsPtBalanceEventJetPtL2L3ResCorrected(const Event *evt) const;
  double twoJetsPtBalanceEventJet2Pt(const Event *evt) const;
  double twoJetsPtBalanceEventJet2PtL2L3Corrected(const Event *evt) const;
  double twoJetsPtBalanceEventJet2PtL2L3ResCorrected(const Event *evt) const;
  double twoJetsPtBalanceEventJetLeadPt(const Event *evt) const;
  double twoJetsPtBalanceEventJetLeadPtL2L3Corrected(const Event *evt) const;
  double twoJetsPtBalanceEventJetLeadPtL2L3ResCorrected(const Event *evt) const;
  double twoJetsPtBalanceEventJetLead2Pt(const Event *evt) const;
  double twoJetsPtBalanceEventJetLead2PtL2L3Corrected(const Event *evt) const;
  double twoJetsPtBalanceEventJetLead2PtL2L3ResCorrected(const Event *evt) const;  
  double twoJetsPtBalanceEventMeanPt(const Event *evt) const;
  double twoJetsPtBalanceEventJetEMF(const Event *evt) const; 
  double twoJetsPtBalanceEventMET(const Event *evt) const;
  double twoJetsPtBalanceEventVtxN(const Event *evt) const;
  double twoJetsPtBalanceEventMCNPUVtx(const Event *evt) const;
  double twoJetsPtBalanceEventMCNPUTruth(const Event *evt) const;
  double twoJetsPtBalanceEventPF_CH_Fraction(const Event * evt) const;
  double twoJetsPtBalanceEventPF_NH_Fraction(const Event * evt) const;
  double twoJetsPtBalanceEventPF_PH_Fraction(const Event * evt) const;
  double twoJetsPtBalanceEventPF_EL_Fraction(const Event * evt) const;
  double twoJetsPtBalanceEventPF_HFHad_Fraction(const Event * evt) const;
  double twoJetsPtBalanceEventPF_HFEm_Fraction(const Event * evt) const;
  double twoJetsPtBalanceEventJetFlavor(const Event *evt) const;
  double twoJetsPtBalanceEventThirdJetFraction(const Event *evt) const;
  double twoJetsPtBalanceEventThirdJetFractionPlain(const Event *evt) const;
  double twoJetsPtBalanceEventThirdJetPt(const Event *evt) const;
  double twoJetsPtBalanceEventJetMomentPhiPhi(const Event *evt) const;
  double twoJetsPtBalanceEventJetMomentEtaEta(const Event *evt) const; 
  double twoJetsPtBalanceEventJetMeanMoment(const Event *evt) const;
  double twoJetsPtBalanceEventDeltaPhi(const Event * evt) const;
  double twoJetsPtBalanceEventClosestJetdR(const Event * evt) const;
  double twoJetsPtBalanceEventAsymmetry(const Event * evt) const;
  double twoJetsPtBalanceEventAsymmetryKalibriCorrected(const Event * evt) const;
  double twoJetsPtBalanceEventAsymmetryL2L3Corrected(const Event * evt) const;
  double twoJetsPtBalanceEventAsymmetryL2L3ResCorrected(const Event * evt) const;
  double twoJetsPtBalanceEventAsymmetryL2L3L4Corrected(const Event * evt) const;
  double twoJetsPtBalanceEventAsymmetryL2L3ResL4Corrected(const Event * evt) const;
  double twoJetsPtBalanceEventGenAsymmetry(const Event * evt) const ;
  double twoJetsPtBalanceEventB(const Event * evt) const;
  double twoJetsPtBalanceEventBKalibriCorrected(const Event * evt) const;
  double twoJetsPtBalanceEventBL2L3Corrected(const Event * evt) const;
  double twoJetsPtBalanceEventBL2L3ResCorrected(const Event * evt) const;
  double twoJetsPtBalanceEventBL2L3L4Corrected(const Event * evt) const;
  double twoJetsPtBalanceEventBL2L3ResL4Corrected(const Event * evt) const;
  double twoJetsPtBalanceEventMPFResponse(const Event * evt) const;
  double twoJetsPtBalanceEventMPFResponseL2L3Corrected(const Event * evt) const;
  double twoJetsPtBalanceEventMPFResponseL2L3ResCorrected(const Event * evt) const;
  double twoJetsPtBalanceEventMCTruthJet1Response(const Event * evt) const;
  double twoJetsPtBalanceEventMCTruthJet1L2L3Response(const Event * evt) const;

 private:
  //! The different correction types of the y quantity
  std::vector<ControlPlotsConfig::CorrectionType> types_;
  //! The binning function
  Function binFunc_;
  //! The x value function
  Function xFunc_; 
  //! The cut value function
  Function cutFunc_;
  //! The y value functions for the different correction types
  std::map<ControlPlotsConfig::CorrectionType,Function> yFuncs_;    
};
#endif
