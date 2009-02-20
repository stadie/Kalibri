//
//    Class for all events with one jet and truth informatio
//
//    first version: Hartmut Stadie 2008/12/14
//    $Id: JetTruthEvent.h,v 1.3 2009/01/22 17:48:10 stadie Exp $
//   
#ifndef JETTRUTHEVENT_H
#define JETTRUTHEVENT_H

#include"CalibData.h"

#include "Jet.h"

//interface to Data
class JetTruthEvent : public TData
{
public:
  JetTruthEvent(Jet *j, double t, double w) : jet(j),truth(t),weight(w),flagged_bad(false) {}
  ~JetTruthEvent() { delete jet;}

  //interface from TData
  TMeasurement *GetMess() const {return jet;}
  double GetTruth() const { return truth;}
  double GetParametrizedMess() const { return jet->correctedEt(jet->Et());}

  void ChangeParAddress(double* oldpar, double* newpar) { jet->ChangeParAddress(oldpar,newpar);}
  DataType GetType() const { return GammaJet;} 
  double GetWeight() const { return weight;}
  
  double chi2() const;
  double chi2_fast(double * temp_derivative1, double * temp_derivative2, double const epsilon) const { 
    return chi2_fast_invert(temp_derivative1,temp_derivative2,epsilon);
  }
  double chi2_fast_blobel(double * temp_derivative1, double * temp_derivative2, double const epsilon) const;
  double chi2_fast_simple_scaled(double * temp_derivative1, double * temp_derivative2, double const epsilon) const;
  double chi2_fast_simple(double * temp_derivative1, double * temp_derivative2, double const epsilon) const;
  double chi2_fast_invert(double * temp_derivative1, double * temp_derivative2, double const epsilon) const;
  void UpdateError() {}

  static void printStats();
 private:
  Jet* jet;
  double truth;
  double weight;
  mutable bool flagged_bad;
  static int nflagged;
};

#endif
