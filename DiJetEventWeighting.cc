//
//    Class for existing event processors
//    
//    Right now we only have the spectrum correctors.
//    Thus they are implemented directly in this class
//
//    first version: Hartmut Stadie 2008/12/14
//    $Id: DiJetEventWeighting.cc,v 1.5 2012/01/24 16:26:47 kirschen Exp $
//   
#include "DiJetEventWeighting.h"

#include "CalibData.h"
#include "ConfigFile.h"
#include "Parameters.h"
#include "TwoJetsPtBalanceEvent.h"


#include <iostream>

DiJetEventWeighting::DiJetEventWeighting(const std::string& configfile, Parameters* param)
  : EventProcessor("DiJetEventWeighting",configfile,param)
{  
  ConfigFile config(configfile.c_str()); 
    
  std::vector<double> trigthresholds = bag_of<double>(config.read<std::string>("Di-Jet trigger thresholds",""));

 
  for(int i = 0, l = trigthresholds.size() ; i < l ; ++i) {
    ndata_[trigthresholds[i]] = 0;
    ncontrol_[trigthresholds[i]] = 0;
  }
}
 
DiJetEventWeighting::~DiJetEventWeighting()
{
}
  

int DiJetEventWeighting::preprocess(std::vector<Event*>& data,
			     std::vector<Event*>& control1,
			     std::vector<Event*>& control2)
{

  std::cout << "start DiJetEventWeighting with:" << std::endl;
  std::cout << "  " << data.size() << " events in data" << std::endl;
  std::cout << "  " << control1.size() << " events in control1" << std::endl;
  //  std::cout << data.size() << " events in data" << std::endl;

  //count events in data and control sample

  //delete events that do not pass trigger thresholds (after JEC)
  CheckBadEventDiJetEventWeighting passCheckBadEventDiJetEventWeighting(this);
  std::vector<Event*>::iterator bound;
  bound= partition(data.begin(),data.end(),passCheckBadEventDiJetEventWeighting);
  data.erase(bound,data.end());

  std::cout << "after checking for bad events:" << std::endl;
  std::cout << "  " << data.size() << " events in data" << std::endl;


  //delete events that have a strange response (response <0.2 or >2.0), large deviation from pthat or a bad dR reco-gen-match.
  CheckStrangeEventDiJetEventWeighting passCheckStrangeEventDiJetEventWeighting(this);
  bound= partition(control1.begin(),control1.end(),passCheckStrangeEventDiJetEventWeighting);
  control1.erase(bound,control1.end());

  std::cout << "after checking for strange response:" << std::endl;
  std::cout << "  " << control1.size() << " events in control1" << std::endl;

  //compute new weight factors
  for(std::map<double,double>::const_iterator i = ndata_.begin() ;
      i != ndata_.end() ; ++i) {
    weights_[i->first] = ncontrol_[i->first] ? i->second/ncontrol_[i->first]: 0.0;
  }
  
  //set weights and delete events that would be weigthed to zero (due to trigger thresholds)
  SetWeightDiJetEventWeighting passSetWeightDiJetEventWeighting(this);
  bound= partition(control1.begin(),control1.end(),passSetWeightDiJetEventWeighting);
  control1.erase(bound,control1.end());

  std::cout << "end DiJetEventWeighting with:" << std::endl;
  std::cout << "  " << data.size() << " events in data" << std::endl;
  std::cout << "  " << control1.size() << " events in control1" << std::endl;


  return ndata_.size();
}
 

bool DiJetEventWeighting::passCheckBadEventDiJetEventWeighting(Event* event)
{
  if(event->type() != PtBalance) std::cout << "Warning: No TwoJetsPtBalanceEvent! ";
  TwoJetsPtBalanceEvent* tje = dynamic_cast<TwoJetsPtBalanceEvent*>(event);
  std::map<double,double>::iterator it = ndata_.lower_bound(tje->ptDijetCorrL2L3());
  if(!(it == ndata_.begin())){
    assert(it != ndata_.begin());
    --it;
    it->second += tje->weight();
    return true;
  }
  else {
    //	std::cout << "bad event ptDijetCorrL2L3: " << tje->ptDijetCorrL2L3() << " ptDijet(): " << tje->ptDijet() << " ptDijetGen: " << tje->ptDijetGen() <<  " lowest threshold: " <<it->first << '\n';
    return false;
  }
}

bool DiJetEventWeighting::passCheckStrangeEventDiJetEventWeighting(Event* event)
{
  if(event->type() != PtBalance) std::cout << "Warning: No TwoJetsPtBalanceEvent! ";
  TwoJetsPtBalanceEvent* tje = dynamic_cast<TwoJetsPtBalanceEvent*>(event);
  Jet * j1 = tje->getJet1();
  double res1 = j1->pt()/j1->genPt();
  double dr1 = j1->dR();
  Jet * j2 = tje->getJet2();
  double res2 = j2->pt()/j2->genPt();
  double dr2 = j2->dR();
  //      //cut on pthat introduced according to JES-mail by Mikko 11 Aug 2011
  //      //NOW: Introduce additional deltaR-matching cut...
  if(!( (res1 > 0.2) && (res1 < 2.0) && (res2 > 0.2) && (res2 < 2.0) && (j1->pt()< 2.0 * tje->ptHat())  && (j2->pt()< 2.0 * tje->ptHat())&& (dr1<0.25) && (dr2 <0.25) )){
    //      std::cout << "strange response: " << res1 << ", " <<  res2 << " pt:" <<  j1->pt() << ", " << j2->pt() << ", pthat: "  << tje->ptHat() << '\n';
    return false;
  }
  std::map<double,double>::iterator it = ncontrol_.lower_bound(tje->ptDijetCorrL2L3());
  if(!(it == ncontrol_.begin())){
    assert(it != ncontrol_.begin());
    --it;
    it->second += tje->weight();
  }
  return true;
}
 

bool DiJetEventWeighting::passSetWeightDiJetEventWeighting(Event* event)
{
  if(event->type() != PtBalance) std::cout << "Warning: No TwoJetsPtBalanceEvent! ";
  TwoJetsPtBalanceEvent* tje = dynamic_cast<TwoJetsPtBalanceEvent*>(event);
  std::map<double,double>::iterator it = weights_.lower_bound(tje->ptDijetCorrL2L3());
  if(it == weights_.begin()) {
    return false;
  }
  else {
    --it;
    tje->setWeight(tje->weight() * it->second);
    return true;
  }
  
}
