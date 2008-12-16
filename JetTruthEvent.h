//
//    Class for all events with one jet and truth informatio
//
//    first version: Hartmut Stadie 2008/12/14
//    $Id: EventReader.h,v 1.1 2008/12/12 13:43:15 stadie Exp $
//   
#ifndef JETTRUTHEVENT_H
#define JETTRUTHEVENT_H

#include"CalibData.h"

#include "Jet.h"

//interface to Data
class JetTruthEvent : public TData
{
public:
  JetTruthEvent(Jet *j, double t, double w) : jet(j),truth(t),weight(w) {}
  ~JetTruthEvent() { delete jet;}

  //interface from TData
  TMeasurement *GetMess() const {return jet;}
  double GetTruth() const { return truth;}
  double GetParametrizedMess() const { return jet->correctedEt(jet->Et());}

  void ChangeParAddress(double* oldpar, double* newpar) { jet->ChangeParAddress(oldpar,newpar);}
  DataType GetType() const { return GammaJet;} 
  double GetWeight() const { return weight;}
  
  double chi2() const;
  double chi2_fast(double * temp_derivative1, double * temp_derivative2, double const epsilon) const;
  void UpdateError() {}

 private:
  Jet* jet;
  double truth;
  double weight;
};

#endif
