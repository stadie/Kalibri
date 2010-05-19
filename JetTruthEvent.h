//
//    Class for all events with one jet and truth informatio
//
//    first version: Hartmut Stadie 2008/12/14
//    $Id: JetTruthEvent.h,v 1.15 2010/05/19 13:34:49 stadie Exp $
//   
#ifndef JETTRUTHEVENT_H
#define JETTRUTHEVENT_H

#include"CalibData.h"

#include "Jet.h"

//interface to Data
class JetTruthEvent : public Event
{
public:
 JetTruthEvent(Jet *j, double t, double w, bool binned = false) : Event(w),jet_(j),truth_(t),chi2plots_(1000.),flagged_bad_(false), binned_(binned) {}
  ~JetTruthEvent();

  //interface from Event
  Measurement *mess() const {return jet_;}
  double truth() const { return truth_;}
  double parametrizedMess() const { return jet_->correctedEt(jet_->Et());}

  void changeParAddress(double* oldpar, double* newpar) { jet_->changeParAddress(oldpar,newpar);}
  DataType type() const { return GammaJet;} 
  Jet* jet() const {return jet_;}
  
  double chi2() const;
  double chi2_plots() const { return chi2plots_; }
  double chi2_fast(double * temp_derivative1, double * temp_derivative2, const double* epsilon) const { 
    if(! binned_) {
      chi2plots_ = chi2_log_fast_invert(temp_derivative1,temp_derivative2,epsilon);
    } else {
      chi2plots_ = chi2_fast_scaled(temp_derivative1,temp_derivative2,epsilon);
    }
    return chi2plots_;
  }
  double chi2_fast_blobel(double * temp_derivative1, double * temp_derivative2, const double* epsilon) const;
  double chi2_fast_scaled(double * temp_derivative1, double * temp_derivative2, const double* epsilon) const;
  double chi2_fast_simple_scaled(double * temp_derivative1, double * temp_derivative2, const double* epsilon) const;
  double chi2_fast_simple(double * temp_derivative1, double * temp_derivative2, const double* epsilon) const;
  double chi2_fast_invert(double * temp_derivative1, double * temp_derivative2, const double* epsilon) const;
  double chi2_log_fast_invert(double * temp_derivative1, double * temp_derivative2, const double* epsilon) const;
  void updateError() { }
  bool flaggedBad() const { return flagged_bad_; }  //!< Status from inversion procedure

  static int nFlagged() { return nflagged_;}
  static void printStats();
 private:
  Jet* jet_;
  double truth_;
  mutable double chi2plots_;   //!< Store chi2 value from last iteration for plots
  mutable bool flagged_bad_;
  bool binned_;
  static int nflagged_;
};

#endif
