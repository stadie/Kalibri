//
//    Class for existing event processors
//    
//    Right now we only have the spectrum correctors.
//    Thus they are implemented directly in this class
//
//    first version: Hartmut Stadie 2011/05/30
//    $Id: DiJetEventWeighting.cc,v 1.2 2011/01/19 15:05:21 stadie Exp $
//   
#include "DiJetEventBinning.h"

#include "Parameters.h"
#include "TwoJetsPtBalanceEvent.h"
#include "TwoJetsPtBalanceBinnedEvent.h"

#include <iostream>

DiJetEventBinning::DiJetEventBinning(const std::string& configfile, Parameters* param) 
  : EventProcessor("DiJetEventBinning",configfile,param), minEvents_(100) 
{
  ConfigFile config(configfile.c_str());
  etavalues_ = bag_of<double>(config.read<std::string>(name() +" eta bins","0"));
  for(std::vector<double>::const_iterator i = etavalues_.begin() ; 
      i !=  etavalues_.end() ; ++i) {
    for(int j = 0 ; j < Nalpha_ ; ++j) {
      databins_[j].insert(std::make_pair(EtaPtKey(*i,0),new TwoJetsPtBalanceBinnedEvent((j+1) *0.05)));  
      controlbins_[j].insert(std::make_pair(EtaPtKey(*i,0),new TwoJetsPtBalanceBinnedEvent((j+1) *0.05)));    
    }
  } 
}

DiJetEventBinning::~DiJetEventBinning() {}
  

DiJetEventBinning::EventMap::const_iterator DiJetEventBinning::getIterator(double eta, double pt, 
					       const EventMap& map) const {
  std::vector<double>::const_iterator i = std::lower_bound(etavalues_.begin(), 
							   etavalues_.end(), eta);
  if(i == etavalues_.end()) {
    std::cout << "cannot find lower bound for " << eta
	      << '\n';
    return map.end();
  }
  --i;
  EventMap::const_iterator j = map.lower_bound(EtaPtKey(*i,pt));
  --j;
  //std::cout << "returning key for " << eta << ", " << pt << " key:" << j->first.eta() 
  //	    << ", " << j->first.pt() << " low eta:" << *i << '\n'; 
  return j;
}


int DiJetEventBinning::preprocess(std::vector<Event*>& data,
			 std::vector<Event*>& control1,
			 std::vector<Event*>& control2) 
{ 
  std::vector<Event*> kept;
  for(std::vector<Event*>::iterator i = data.begin() ; i != data.end() ; ++i) {
    if((*i)->type() != PtBalance) {
      kept.push_back(*i);
      continue;
    }
    TwoJetsPtBalanceEvent* tje = dynamic_cast<TwoJetsPtBalanceEvent*>(*i);
    const double alpha = tje->relPtJet3();
    for(int j = 0 ; j < Nalpha_ ; ++j) {
      if(alpha < (j+1) *0.05) {
	TwoJetsPtBalanceBinnedEvent *bin = getEvent(tje->jet1()->eta(),tje->ptDijet(),databins_[j]);
	bin->addEvent(tje);
	if((j == 0) && (bin->size() == minEvents_*2)) {
	  split(tje->jet1()->eta(),tje->ptDijet());
	}
      }
    }
  }
  data = kept;
  kept.clear();
  for(std::vector<Event*>::iterator i = control1.begin() ; i != control1.end() ; ++i) {
    if((*i)->type() != PtBalance) {
      kept.push_back(*i);
      continue;
    }
    TwoJetsPtBalanceEvent* tje = dynamic_cast<TwoJetsPtBalanceEvent*>(*i);
    const double alpha = tje->relPtJet3();
    for(int j = 0 ; j < Nalpha_ ; ++j) {
      if(alpha < (j+1) *0.05) {
	TwoJetsPtBalanceBinnedEvent *bin = getEvent(tje->jet1()->eta(),tje->ptDijet(),controlbins_[j]);
	bin->addEvent(tje);
      }
    }
  }
  control1 =  kept;
  //add new events
  for(int i = 0 ; i < Nalpha_ ; ++i) {
    for(EventMap::iterator j = databins_[i].begin() ; 
	j != databins_[i].end() ; ++j) {
      TwoJetsPtBalanceBinnedEvent *bin = j->second;
      //std::cout << j->first.eta() << ", " << j->first.pt() << ":" << bin->size() << '\n';
      if(bin->size() > 0) {
	bin->computeJets();
	data.push_back(bin);
      }
    }
    for(EventMap::iterator j = controlbins_[i].begin() ; 
	j != controlbins_[i].end() ; ++j) {
      TwoJetsPtBalanceBinnedEvent *bin = j->second;
      //std::cout << j->first.eta() << ", " << j->first.pt() << ":" << bin->size() << '\n';;
      if(bin->size() > 0) {
	bin->computeJets();
	control1.push_back(bin);
      }
    }
  }
  return data.size();
}


void DiJetEventBinning::split(double eta, double pt) {
  double minpt = 0;
  for(int i = 0 ; i < Nalpha_ ; ++i) {
    EventMap::const_iterator it = getIterator(eta,pt,databins_[i]);
    TwoJetsPtBalanceBinnedEvent *bin = it->second;
    if(i == 0) minpt = bin->ptMedian();
    TwoJetsPtBalanceBinnedEvent *newbin = bin->split(minpt);
    //if(i == 0) {
    //  std::cout << "splitted event at " << it->first.eta() << ", " << minpt << " " 
    //		<< bin->size() << ", " << newbin->size() << '\n';
    //}
    databins_[i].insert(std::make_pair(EtaPtKey(it->first.eta(),minpt),newbin));
    it = getIterator(eta,pt,controlbins_[i]);
    bin = it->second;
    newbin = bin->split(minpt);
    controlbins_[i].insert(std::make_pair(EtaPtKey(it->first.eta(),minpt),newbin));
  }
}

