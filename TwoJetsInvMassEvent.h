//
//    Class for all events with one jet and truth informatio
//
//    first version: Hartmut Stadie 2008/12/14
//    $Id: TwoJetsInvMassEvent.h,v 1.3 2009/01/22 17:48:10 stadie Exp $
//   
#ifndef TWOJETSINVMASSEVENT_H
#define TWOJETSINVMASSEVENT_H

#include"CalibData.h"

#include "Jet.h"

//interface to Data
class TwoJetsInvMassEvent : public TData
{
public:
  TwoJetsInvMassEvent(Jet *j1, Jet *j2, double t, double w) 
    : jet1(j1), jet2(j2),truth(t),weight(w),flagged_bad(false) {}
  ~TwoJetsInvMassEvent() { delete jet1; delete jet2;}

  //interface from TData
  TMeasurement *GetMess() const {return jet1;}
  double GetTruth() const { return truth;}
  double GetParametrizedMess() const { return jet1->correctedEt(jet1->Et());}
  
  TMeasurement *GetMess2() const {return jet2;}
  double GetParametrizedMess2() const { return jet2->correctedEt(jet2->Et());}

  void ChangeParAddress(double* oldpar, double* newpar) { 
    jet1->ChangeParAddress(oldpar,newpar);
    jet2->ChangeParAddress(oldpar,newpar);
  }
  DataType GetType() const { return InvMass;} 
  double GetWeight() const { return weight;}

  double correctedMass() const;
  
  double chi2() const;
  double chi2_fast(double * temp_derivative1, double * temp_derivative2, double const epsilon) const { 
    return chi2_fast_simple(temp_derivative1,temp_derivative2,epsilon);
  }
  double chi2_fast_simple(double * temp_derivative1, double * temp_derivative2, double const epsilon) const;
  void UpdateError() {}

 private:
  Jet *jet1,*jet2;
  double truth;
  double weight;
  bool flagged_bad;
};

#endif
