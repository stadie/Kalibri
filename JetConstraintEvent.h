//
//    Class for constraints on the jet correction
//
//    first version: Hartmut Stadie 2009/07/23
//    $Id: JetConstraintEvent.h,v 1.8 2010/05/19 13:34:48 stadie Exp $
//   
#ifndef JETCONSTRAINTEVENT_H
#define JETCONSTRAINTEVENT_H

#include <map>
#include <vector>

#include"CalibData.h"
#include "Jet.h"

//interface to Data
class JetConstraintEvent : public Event
{
 public:
  JetConstraintEvent(double minpt, double maxpt, double mineta, double maxeta, 
		     double w) 
    : Event(w), minpt_(minpt),maxpt_(maxpt),mineta_(mineta),maxeta_(maxeta),trusum_(0),
    trusum2_(0),error_(0) {}
  ~JetConstraintEvent();
    
  void addJet(double truePt, const Jet* j, const Function* globalFunc = 0);
  

  //interface from TData
  Measurement *mess() const {return jets_[0];}
  double truth() const { return jets_.size() ? trusum_ / jets_.size() : 0;}
  double parametrizedMess() const { return jets_[0]->correctedEt(jets_[0]->Et());}

  void changeParAddress(double* oldpar, double* newpar);
  DataType type() const { return JetConstraint;} 
  
  double chi2() const;
  double chi2_plots() const { return chi2plots_; }
  double chi2_fast(double * temp_derivative1, double * temp_derivative2, const double* epsilon) const;
  void updateError() { } 
  double minEta() const { return mineta_;}
  double maxEta() const { return maxeta_;}
  double minPt() const { return minpt_;}
  double maxPt() const { return maxpt_;}
  int nJets() const { return jets_.size();}
 private:
  struct Variation{
    double uppersum_;
    double lowersum_;
    Variation() : uppersum_(0), lowersum_(0) {}
  };
  double minpt_,maxpt_;
  double mineta_,maxeta_;
  typedef std::map<int, Variation> VarMap;
  mutable VarMap varmap_;
  std::vector<Jet*> jets_;
  double trusum_;
  double trusum2_;
  double error_;
  mutable double chi2plots_;   //!< Store chi2 value from last iteration for plots
};

#endif
