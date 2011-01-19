//
//    Class for existing event processors
//    
//    Right now we only have the spectrum correctors.
//    Thus they are implemented directly in this class
//
//    first version: Hartmut Stadie 2008/12/14
//    $Id: DiJetEventWeighting.cc,v 1.1 2010/12/20 11:08:13 stadie Exp $
//   
#include "DiJetEventWeighting.h"

#include "CalibData.h"
#include "ConfigFile.h"
#include "Parameters.h"
#include "TwoJetsPtBalanceEvent.h"


#include <iostream>

DiJetEventWeighting::DiJetEventWeighting(const std::string& configfile, Parameters* param)
  : EventProcessor(configfile,param)
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
  //count events in data and control sample
  for(std::vector<Event*>::iterator i = data.begin() ; i != data.end() ; ++i) {
    if((*i)->type() != PtBalance) continue;
    TwoJetsPtBalanceEvent* tje = dynamic_cast<TwoJetsPtBalanceEvent*>(*i);
    std::map<double,double>::iterator it = ndata_.lower_bound(tje->ptDijetCorrL2L3());
    if(it == ndata_.begin()) {
      std::cout << "bad event " << tje->ptDijetCorrL2L3() << " " << it->first << '\n';
      --i;
      delete tje;
      data.erase(i+1);
      continue;
    }
    assert(it != ndata_.begin());
    --it;
    it->second += tje->weight();
  }
  for(std::vector<Event*>::const_iterator i = control1.begin() ; i != control1.end() ; ++i) {
    if((*i)->type() != PtBalance) continue;
    TwoJetsPtBalanceEvent* tje = dynamic_cast<TwoJetsPtBalanceEvent*>(*i);
    std::map<double,double>::iterator it = ncontrol_.lower_bound(tje->ptDijetCorrL2L3());
    if(it == ncontrol_.begin()) continue;
    --it;
    it->second += tje->weight();
  }
  
  //compute new weight factors
  std::map<double,double> weights;
  for(std::map<double,double>::const_iterator i = ndata_.begin() ;
      i != ndata_.end() ; ++i) {
    weights[i->first] = ncontrol_[i->first] ? i->second/ncontrol_[i->first]: 0.0;
  }
  
  //set weights
  for(std::vector<Event*>::iterator i = control1.begin() ; i != control1.end() ; ++i) {
    if((*i)->type() != PtBalance) continue;
    TwoJetsPtBalanceEvent* tje = dynamic_cast<TwoJetsPtBalanceEvent*>(*i);
    std::map<double,double>::iterator it = weights.lower_bound(tje->ptDijetCorrL2L3());
    if(it == weights.begin()) {
      --i;
      delete tje;
      control1.erase(i+1);
      continue;
    };
    --it;
    //std::cout << tje->ptDijetCorrL2L3() << " " << it->first << "  w=" << it->second << '\n';
    tje->setWeight(tje->weight() * it->second);
  }
  
  return ndata_.size();
}
 

