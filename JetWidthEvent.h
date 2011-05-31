//
//    Class for all events with one jet and truth informatio
//
//    first version: Hartmut Stadie 2008/12/14
//    $Id: JetWidthEvent.h,v 1.3 2011/05/18 15:58:35 stadie Exp $
//   
#ifndef JETWIDTHEVENT_H
#define JETWIDTHEVENT_H

#include "CalibData.h"

#include "Jet.h"

class Parameters;

//interface to Data
class JetWidthEvent : public Event
{
public:
 JetWidthEvent(Jet *j, double w, double pthat = 0, short npu = 0) : Event(w,pthat,npu),jet_(j),chi2plots_(1000.) {}
  ~JetWidthEvent();

  //interface from Event
  Measurement *mess() const {return jet_;}
  double truth() const { return 0.0;}
  double parametrizedMess() const { return jet_->correctedEt(jet_->Et());}

  void setParameters(Parameters* param) { jet_->setParameters(param);}
  DataType type() const { return JWFit;} 
  Jet* jet() const {return jet_;}
  
  double chi2() const;
  double chi2_plots() const { return chi2plots_; }
  double chi2_fast(double * temp_derivative1, double * temp_derivative2,  double * temp_derivative3, double * temp_derivative4, const double* epsilon) const; 
  void updateError() {}
 private:
  Jet* jet_;
  mutable double chi2plots_;   //!< Store chi2 value from last iteration for plots
};

#endif
