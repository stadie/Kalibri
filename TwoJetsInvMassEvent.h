//
//    Class for all events with one jet and truth informatio
//
//    first version: Hartmut Stadie 2008/12/14
//
//    $Id: TwoJetsInvMassEvent.h,v 1.12 2010/11/01 15:47:40 stadie Exp $
//   

#ifndef TWOJETSINVMASSEVENT_H
#define TWOJETSINVMASSEVENT_H

#include"CalibData.h"

#include "Jet.h"

//interface to Data
class TwoJetsInvMassEvent : public Event
{
public:
  TwoJetsInvMassEvent(Jet *j1, Jet *j2, double t, double w, double* p) 
    : jet1_(j1), jet2_(j2),truth_(t),flagged_bad_(false),chi2plots_(1000.),par_(p) {}
  ~TwoJetsInvMassEvent() { delete jet1_; delete jet2_;}

  //interface from Event
  Measurement *mess() const {return jet1_;}
  double truth() const { return truth_;}
  double parametrizedMess() const { return jet1_->correctedEt(jet1_->Et());}
  
  Measurement *mess2() const {return jet2_;}
  double parametrizedMess2() const { return jet2_->correctedEt(jet2_->Et());}

  Jet *jet1() const {return jet1_;}
  Jet *jet2() const {return jet2_;}

  void setParameters(Parameters* param);
  DataType type() const { return InvMass;} 
  double correctedMass() const;
  
  double chi2() const;
  double chi2_plots() const { return chi2plots_; }
  double chi2_fast(double * temp_derivative1, double * temp_derivative2,  double * temp_derivative3, double * temp_derivative4, const double* epsilon) const { 
    //chi2plots = chi2_fast_simple(temp_derivative1,temp_derivative2,epsilon);
    //chi2plots = chi2_fast_scaled(temp_derivative1,temp_derivative2,epsilon);
    //chi2plots = chi2_fast_const_error(temp_derivative1,temp_derivative2,epsilon);
    chi2plots_ = chi2_fast_inv(temp_derivative1,temp_derivative2,epsilon);
    return chi2plots_;
  }
  double chi2_fast_simple(double * temp_derivative1, double * temp_derivative2, const double* epsilon) const;
  double chi2_fast_const_error(double * temp_derivative1, double * temp_derivative2, const double* epsilon) const;
  double chi2_fast_scaled(double * temp_derivative1, double * temp_derivative2, const double* epsilon) const;
  double chi2_fast_inv(double * temp_derivative1, double * temp_derivative2, const double* epsilon) const;
  void updateError() {}

 private:
  Jet *jet1_,*jet2_;
  double truth_;
  mutable bool flagged_bad_;
  mutable double chi2plots_;   //!< Store chi2 value from last iteration for plots
  double *par_;
};

#endif
