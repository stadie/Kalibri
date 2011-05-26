// $Id:  $

#ifndef DIJET_RESOLUTION_EVENT_H
#define DIJET_RESOLUTION_EVENT_H

#include "CalibData.h"
#include "Jet.h"
#include "ResolutionFunction.h"

//!  \brief A dijet event for resolution measurement
//!  \author Matthias Schroeder
//!  \date Tue Jun  9 15:24:49 CEST 2009
//!  $Id: $
// --------------------------------------------------
class DiJetResolutionEvent : public Event {
public:
  DiJetResolutionEvent(Jet* jet1, Jet* jet2, double deltaPhi12, double pPhi,
		       double ptJet3, double ptJet4, double pJ3, double pSJ, double ptRef,
		       double ptHat, double weight, const ResolutionFunction& pdf,
		       double min, double max, double eps, int niter);
  ~DiJetResolutionEvent();

  // Dummy methods to implement event interface
  void setParameters(Parameters* param);
  Measurement * mess() const { return 0; }
  double truth() const { return 0.; }
  double chi2_plots() const { return 0.; }
  double parametrizedMess() const { return 0.; }
  void updateError() { };

  // Implements event interface
  DataType type() const { return DiJetResolution; }
  double chi2() const { return chi2Spectrum(); }
  double chi2_fast(double * temp_derivative1, double * temp_derivative2, const double* epsilon) const;
  void printInitStats() const;

  // DiJetResolutionEvent specific methods
  unsigned int ptBin() const { return pdf_->ptBin(); }
  const Jet* jet1() const { return jet1_; }
  const Jet* jet2() const { return jet2_; }
  const Jet* jet(int i) const { 
    const Jet *jet = 0;
    if( i == 0 ) jet = jet1();
    else if( i == 1 ) jet = jet2();
    return jet;
  }
  double deltaPhi12() const { return deltaPhi12_; }
  double pPhi() const { return pPhi_; }
  double ptJet3() const { return ptJet3_; }
  double ptJet4() const { return ptJet4_; }
  double pJ3() const { return pJ3_; }
  double pSJ() const { return pSJ_; }
  double ptRef() const { return ptRef_; }
  double avePt() const { return 0.5 * (jet1()->pt() + jet2()->pt()); } 
  double avePtGen() const { return 0.5 * (jet1()->genPt() + jet2()->genPt()); } 
  double relJet3Pt() const { return ptJet3() / avePt(); }

  double pdfPtMeas(double ptMeas1, double ptMeas2, double ptTrue) const { return pdf_->pdfPtMeas(ptMeas1,ptMeas2,ptTrue); }
  double pdfPtTrue(double ptTrue) const { return pdf_->pdfPtTrue(ptTrue); }
  double pdfResp(double r, double ptTrue) const { return pdf_->pdfResp(r,ptTrue); }
  double pdfDijetAsym(double a, double ptTrue) const { return pdf_->pdfDijetAsym(a,ptTrue); }



 private:
  const ResolutionFunction* pdf_; 
  Jet* jet1_;
  Jet* jet2_;
  const double deltaPhi12_;
  const double pPhi_;
  const double ptRef_;
  const double ptJet3_;
  const double ptJet4_;
  const double pJ3_;
  const double pSJ_;

  const int    kMaxNIter_;   //!< Max number of iterations in integration
  const double kEps_;        //!< Integration precision for convergence
  const double kMin_;        //!< Minimum of truth pdf
  const double kMax_;        //!< Maximum of truth pdf
  
  double chi2Simple() const;
  double chi2Spectrum() const;
};

#endif
