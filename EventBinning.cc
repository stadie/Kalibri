//
//    Class for existing event processors
//    
//    Right now we only have the spectrum correctors.
//    Thus they are implemented directly in this class
//
//    first version: Hartmut Stadie 2008/12/14
//    $Id: EventBinning.cc,v 1.1 2010/12/13 10:55:09 stadie Exp $
//   
#include "EventBinning.h"

#include "CalibData.h"
#include "ConfigFile.h"
#include "Parameters.h"
#include "Binning.h"
#include "Jet.h"
#include "JetTruthEvent.h"
#include "JetWidthEvent.h"

#include <iostream>
#include <algorithm>

EventBinning::EventBinning(const std::string& configfile, Parameters* param)
  : EventProcessor(configfile,param)
{  
  ConfigFile config(configfile.c_str()); 
  binning_ = new Binning(&config);
 
  binEvents_ = config.read<bool>("bin events",false);
  eventsperbin_ = config.read<int>("events per bin",1000);
  for(int i = 0 ; i < 4 ; ++i) jetfunc_[i] = 0;
  std::vector<std::string> vars = bag_of_string(config.read<std::string>("jet binning variables","")); 
  int j = 0;
  for(std::vector<std::string>::const_iterator i = vars.begin() ; 
      i != vars.end() ; ++i) {
    if(*i == "pt") {
      jetfunc_[j] = &Jet::genPt; 
    } else if(*i == "eta") {
      jetfunc_[j] = &Jet::eta;
    } else if(*i == "sigmaphi") {
      jetfunc_[j] = &Jet::momentPhiPhi;    
    } else if(*i == "sigmaeta") {
      jetfunc_[j] = &Jet::momentEtaEta;   
    } else if(*i == "emf") {
      jetfunc_[j] = &Jet::emf;
    } else if(*i == "meanMoment") {
      jetfunc_[j] = &Jet::meanMoment;   
    } else {
      std::cerr << "unknown binning variable: " << *i << '\n';
      exit(3);
    }
    ++j;
  }
}
 
EventBinning::~EventBinning()
{
  delete binning_;
}
  

int EventBinning::preprocess(std::vector<Event*>& data,
			     std::vector<Event*>& control1,
			     std::vector<Event*>& control2)
{
  if(! binEvents_) return data.size();
  float x[4];
  for(DataIter i = data.begin() ; i != data.end() ; ++i) {
    if(((*i)->type() != GammaJet) && ((*i)->type() != JWFit))continue;
    const Jet *jet = 0;
    if((*i)->type() == GammaJet) {
      JetTruthEvent* jte = dynamic_cast<JetTruthEvent*>(*i);
      jet = jte->jet();
    } 
    if((*i)->type() == JWFit) {
      JetWidthEvent* jwe = dynamic_cast<JetWidthEvent*>(*i);
      jet = jwe->jet();
    }
    
    for(int j  = 0 ; j < 4 ; ++j) {
      x[j] = (jetfunc_[j]) ? (jet->*(jetfunc_[j]))() : 0;
    }
    int id = binning_->findBin(x[0],x[1],x[2],x[3]);
    //std::cout << "add event to bin " <<  id << " for " << x[0] << ", " << x[1] << ", " << x[2] << '\n';
    bins_[id].push_back(*i);
    --i;
    data.erase(i+1);
  }
  for(BinMapIter i = bins_.begin() ; i != bins_.end() ; ++i) {
    std::sort(i->second.begin(),i->second.end(),ltgenjetpt());
    for(unsigned int j = 0, l = i->second.size() ; j < l ; j+= eventsperbin_) {
      unsigned int k = j + eventsperbin_;
      if(k > l) k = l;
      std::cout << j << " to " << k << " in " << l << '\n';
      if(i->second[j]->type() == GammaJet) {
	data.push_back(createBinnedJetTruthEvent(i->second.begin()+j,i->second.begin()+k));
      }
      if(i->second[j]->type() == JWFit) {
	data.push_back(createBinnedJetWidthEvent(i->second.begin()+j,i->second.begin()+k));
      }
    }
  }

  return data.size();
}
 
int EventBinning::postprocess(std::vector<Event*>& data,
			      std::vector<Event*>& control1,
			      std::vector<Event*>& control2)
{
  if(! binEvents_) return data.size();
  for(DataIter i = data.begin() ; i != data.end() ; ++i) {
    if(((*i)->type() != GammaJet) && ((*i)->type() != JWFit))continue;
    --i;
    data.erase(i+1);
  }
  for(BinMapIter i = bins_.begin() ; i != bins_.end() ; ++i) {
    data.insert(data.end(),i->second.begin(),i->second.end());
  }
  return data.size();
}
 

Event* EventBinning::createBinnedJetTruthEvent(DataIter begin, DataIter end)
{
  std::cout << "new JTE for " << end-begin << " entries.\n";
  return 0;
}

Event* EventBinning::createBinnedJetWidthEvent(DataIter begin,DataIter end)
{
  std::cout << "new JWE for " << end-begin << " entries.\n";
  return 0; 
}

bool EventBinning::ltgenjetpt::operator()(const Event* e1, const Event* e2) const
{
  const Jet *jet1 = 0, *jet2 = 0;
  if(e1->type() == GammaJet) {
    const JetTruthEvent* jte = dynamic_cast<const JetTruthEvent*>(e1);
    jet1 = jte->jet();
    jte = dynamic_cast<const JetTruthEvent*>(e2);
    jet2 = jte->jet();
    
  } 
  if(e1->type() == JWFit) {
    const JetWidthEvent* jwe = dynamic_cast<const JetWidthEvent*>(e1);
    jet1 = jwe->jet();
    jwe = dynamic_cast<const JetWidthEvent*>(e2);
    jet2 = jwe->jet();
  }
  return jet1->genPt() < jet2->genPt();
}
