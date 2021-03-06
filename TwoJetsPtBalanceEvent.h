#ifndef TWO_JETS_PT_BALANCE_EVENT_H
#define TWO_JETS_PT_BALANCE_EVENT_H

#include <cmath>

#include "CalibData.h"
#include "Jet.h"
#include "CorFactors.h"
#include "TMath.h"

//!
//!  \brief Class for relative calibration in pseudorapidity
//!         using dijet events
//!
//!  \author Matthias Schroeder
//!  \date Mon Oct 26 21:03:43 CET 2009 
//!  $Id: TwoJetsPtBalanceEvent.h,v 1.19 2012/09/10 15:44:05 kirschen Exp $
// --------------------------------------------------
class TwoJetsPtBalanceEvent : public Event {
 public:
 TwoJetsPtBalanceEvent(Jet *j1, Jet *j2, Jet *j3, double ptHat, double w, short npu, float nputruth, short nvtx, float metraw, float metrawphi, float metT1, float metT1phi, float metT1res, float metT1resphi,int runNumber, float PUMCHighestSumPt)
   : Event(w,ptHat,npu,nputruth,nvtx,metraw,metrawphi,metT1,metT1phi,metT1res,metT1resphi,runNumber,PUMCHighestSumPt),
    jet1_(j1),
    jet2_(j2),
    jet3_(j3),
    residual_(0),
    varresidual_(0),
    flaggedBad_(false),
    chi2Plots_(1000.) {
    error1_ = jet1_ ? jet1_->error() : 0;
    error2_ = jet2_ ? jet2_->error() : 0;
  }

  virtual ~TwoJetsPtBalanceEvent() { delete jet1_; delete jet2_; if( hasJet3() ) delete jet3_; }

  virtual Measurement *mess() const { return jet1_; }
  virtual double parametrizedMess() const { return jet1_->correctedEt();}
  
  Measurement *mess2() const { return jet2_; }
  double parametrizedMess2() const { return jet2_->correctedEt();}

  Measurement *mess3() const { return jet3_; }
  double parametrizedMess3() const { return hasJet3() ? jet3_->correctedEt() : 0.;}

  Jet *getJet1() const { return jet1_; }
  Jet *getJet2() const { return jet2_; }
  Jet *getJet3() const { return jet3_; }
  Jet *jet1() const { return jet1_; }
  Jet *jet2() const { return jet2_; }
  Jet *jet3() const { return jet3_; }
  
  bool hasJet3() const { return jet3_ != 0 ? true : false; }

  virtual void setParameters(Parameters* param) {
    jet1_->setParameters(param);
    jet2_->setParameters(param);
    if( hasJet3() ) jet3_->setParameters(param);
  }

  virtual double truth() const { return 0.; }
  virtual DataType type() const { return PtBalance; } 

  double correctedMass() const;  
  virtual double chi2() const { return chi2_fast(0, 0, 0, 0, 0); }
  virtual double chi2_plots() const { return chi2Plots_; }
  virtual double chi2_fast(double * temp_derivative1, double * temp_derivative2, double * temp_derivative3, double * temp_derivative4, const double *epsilon) const { 
    //chi2Plots_ = chi2_fast_balance(temp_derivative1,temp_derivative2,epsilon);
    chi2Plots_ = chi2_relative(temp_derivative1,temp_derivative2,epsilon);
    return chi2Plots_;
  }
  virtual void updateError() {
    error1_ = jet1_->expectedError(ptDijetCorr());
    error2_ = jet2_->expectedError(ptDijetCorr());
  }

  double ptDijet() const { return 0.5*(jet1_->pt()+jet2_->pt()); }
  double ptDijetGen() const { return 0.5*(jet1_->genPt()+jet2_->genPt()); }
  double ptDijetCorr() const { return 0.5*(parametrizedMess()+parametrizedMess2()); }
  double ptDijetCorrL2L3() const { 
   if(TMath::IsNaN(0.5*( jet1_->corFactors().getL2L3() * jet1_->pt() + jet2_->corFactors().getL2L3() * jet2_->pt()))){
      std::cout << "L2L3corr: " <<jet1_->corFactors().getL2L3() << " L1-corrected jetpt: " << jet1_->pt() << " jeteta: " << jet1_->eta() << " Jet2...: "  << jet2_->corFactors().getL2L3() << " " << jet2_->pt()  << " " << jet2_->eta() << std::endl;
       }
       return 0.5*( jet1_->corFactors().getL2L3() * jet1_->pt() + jet2_->corFactors().getL2L3() * jet2_->pt() ); }

  double ptBalance() const { return (jet1_->pt() - jet2_->pt()) / ptDijet(); }
  double ptBalanceGen() const { return (jet1_->genPt()-jet2_->genPt()) / ptDijetGen(); }
  double ptBalanceCorr() const { return (parametrizedMess()-parametrizedMess2()) / ptDijetCorr(); }
  double ptBalanceCorrL2L3() const { return ( jet1_->corFactors().getL2L3() * jet1_->pt() - jet2_->corFactors().getL2L3() * jet2_->pt() ) / ptDijetCorrL2L3(); }

  double triggerPtVariableL2L3(bool useSingleJetTriggers) const;

  double ptSumAbs(double pt1, double pt2) const;
  double ptSumAbs() const;
  double ptSumAbsCorr() const;
  double ptSumAbsGen() const;
  double ptSumAbsCorrL2L3() const;

  virtual double relPtJet3() const ;
  virtual double relPtJet3CorrL2L3() const ;
  virtual double relPtJet3Projection() const ;
  virtual double relPtJet3ProjectionCorrL2L3() const ;

  bool flaggedBad() const { return flaggedBad_; }  //!< Status

protected:
  Jet *jet1_;
  Jet *jet2_;
  Jet *jet3_;

 private:
  double error1_;		//!< Store jet1_ error during major iteration
  double error2_;		//!< Store jet2_ error during major iteration
  double residual_;
  double varresidual_;

  mutable bool flaggedBad_;
  mutable double chi2Plots_;   //!< Store chi2 value from last iteration for plots

  double chi2_fast_simple(double * temp_derivative1, double * temp_derivative2, const double* epsilon) const;
  double chi2_fast_simple_res(double pt1, double pt2) const;
  double chi2_fast_simple_dRes2(double pt1, double pt2) const;

  double chi2_fast_balance(double * temp_derivative1, double * temp_derivative2, const double* epsilon) const;
  double chi2_fast_balance_res(double pt1, double pt2) const;
  double chi2_fast_balance_dRes2(double pt1, double pt2) const;  
  double chi2_relative(double * temp_derivative1, double * temp_derivative2, const double* epsilon) const;
  static const float lowThresholdJet1_;
  static const float lowThresholdJet3_;
};

#endif
