//
//    Class for existing event processors
//    
//    Right now we only have the spectrum correctors.
//    Thus they are implemented directly in this class
//
//    first version: Hartmut Stadie 2008/12/14
//    $Id: DiJetEventWeighting.cc,v 1.10 2012/11/20 16:33:46 kirschen Exp $
//   
#include "DiJetEventWeighting.h"

#include "CalibData.h"
#include "ConfigFile.h"
#include "Parameters.h"
#include "TwoJetsPtBalanceEvent.h"


#include <iostream>
#include <iomanip>

DiJetEventWeighting::DiJetEventWeighting(const std::string& configfile, Parameters* param)
  : EventProcessor("DiJetEventWeighting",configfile,param)
{  
  ConfigFile config(configfile.c_str()); 
    
  std::vector<double> trigthresholds = bag_of<double>(config.read<std::string>("Di-Jet trigger thresholds",""));
  useSingleJetTriggers_ = config.read<bool>("Use single jet triggers",false);

 
  for(int i = 0, l = trigthresholds.size() ; i < l ; ++i) {
    ndata_[trigthresholds[i]] = 0;
    ncontrol_[trigthresholds[i]] = 0;
    nCountsData_[trigthresholds[i]] = 0;
    nCountsControl_[trigthresholds[i]] = 0;
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



  // do MC-cleaning on "data" as well, in case the data event vector actually contains MC events, test done on first event in data vector
  
  if(data.front()->type() != PtBalance) std::cout << "Warning: No TwoJetsPtBalanceEvent! ";
  TwoJetsPtBalanceEvent* tje = dynamic_cast<TwoJetsPtBalanceEvent*>(data.front());
  bool dataIsActuallyMC = tje->ptHat() >0 ? true : false;

  std::vector<Event*>::iterator bound;
  if(dataIsActuallyMC){
  CheckStrangeEventDiJetEventWeighting MCDATApassCheckStrangeEventDiJetEventWeighting(this);
  bound= partition(data.begin(),data.end(),MCDATApassCheckStrangeEventDiJetEventWeighting);
  data.erase(bound,data.end());
  std::cout << "after checking for strange response in data (which actually is MC):" << std::endl;
  std::cout << "  " << data.size() << " events in MC labelled as data" << std::endl;
  }



  //delete events that do not pass trigger thresholds (after JEC)
  CheckBadEventDiJetEventWeighting passCheckBadEventDiJetEventWeighting(this);
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



  ///////////////////////////////////////////////////////
  ////count events in control sample for reweighting/////
  ///////////////////////////////////////////////////////
  for(std::vector<Event*>::iterator i = control1.begin() ; i != control1.end() ; ++i) {
    if((*i)->type() != PtBalance) continue;
    TwoJetsPtBalanceEvent* tje = dynamic_cast<TwoJetsPtBalanceEvent*>(*i);
    std::map<double,double>::iterator it = ncontrol_.lower_bound(tje->triggerPtVariableL2L3(useSingleJetTriggers_));
    if(!(it == ncontrol_.begin())){
      assert(it != ncontrol_.begin());
      --it;
      if(tje->relPtJet3CorrL2L3()<0.2){
	it->second += tje->weight();
	(--nCountsControl_.lower_bound(tje->triggerPtVariableL2L3(useSingleJetTriggers_)))->second+=1;
      }
    }
  }
  ///////////////////////////////////////////////////////
  ////count events in data sample for reweighting   /////
  ///////////////////////////////////////////////////////
  for(std::vector<Event*>::iterator i = data.begin() ; i != data.end() ; ++i) {
    if((*i)->type() != PtBalance) continue;
    TwoJetsPtBalanceEvent* tje = dynamic_cast<TwoJetsPtBalanceEvent*>(*i);
    std::map<double,double>::iterator it = ndata_.lower_bound(tje->ptDijetCorrL2L3());
    assert(it != ndata_.begin());
    --it;
    if(tje->relPtJet3CorrL2L3()<0.2){
      it->second += tje->weight();
      (--nCountsData_.lower_bound(tje->triggerPtVariableL2L3(useSingleJetTriggers_)))->second+=1;
    }
  }


  std::cout << "------------------------------Summary table of pt-dependent reweighting counts------------------------------" << std::endl;
  std::cout <<std::setw(5) <<"thres" << " => " <<std::setw(12)<< "Weight sum D" << " => " <<std::setw(8)<< "Counts D" << " => " <<std::setw(12)<<  "Avg. w. data"<< " => " <<std::setw(5) <<"thres" << " => " <<std::setw(12)<< "Weight sum C" << " => " <<std::setw(8)<< "Counts C" << " => " <<std::setw(12)<<  "Avg. w. control"<< std::endl;
  for(std::map<double,double>::const_iterator it=ndata_.begin(); it != ndata_.end(); ++it ){
    std::map<double,double>::const_iterator control_it = ncontrol_.find((*it).first);
    std::map<double,int>::const_iterator CountsD_it = nCountsData_.find((*it).first);
    std::map<double,int>::const_iterator CountsC_it = nCountsControl_.find((*it).first);
    std::cout <<std::setw(5) <<(*it).first << " => " <<std::setw(12)<< (*it).second << " => " <<std::setw(8)<< (*CountsD_it).second << " => " <<std::setw(12)<<  (*CountsD_it).second>0 ? (*it).second/(*CountsD_it).second : 0<< " => " 
	      <<std::setw(5)<< (*control_it).first <<" => " <<std::setw(12)<< (*control_it).second  << " => " <<std::setw(8)<< (*CountsC_it).second << " => " <<std::setw(12)<<  (*CountsD_it).second > 0 ?(*control_it).second/(*CountsD_it).second : 0<<std::endl;

  }
  std::cout << "------------------------------------------------------------------------------------------------------------" << std::endl;


  //compute new weight factors
  for(std::map<double,double>::const_iterator i = ndata_.begin() ;
      i != ndata_.end() ; ++i) {
    weights_[i->first] = ncontrol_[i->first] ? i->second/ncontrol_[i->first]: 0.0;
  }
  
  //set weights and delete MC events that would be weigthed to zero (due to trigger thresholds)
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
  std::map<double,double>::iterator it = ndata_.lower_bound(tje->triggerPtVariableL2L3(useSingleJetTriggers_));
  if(!(it == ndata_.begin())){
    assert(it != ndata_.begin());
    --it;
    return true;
  }
  else {
    std::cout << "bad event ptDijetCorrL2L3: " << tje->triggerPtVariableL2L3(useSingleJetTriggers_) << " ptDijet(): " << tje->ptDijet() << " ptDijetGen: " << tje->ptDijetGen() <<  " lowest threshold: " <<it->first << '\n';
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
  std::map<double,double>::iterator it = ncontrol_.lower_bound(tje->triggerPtVariableL2L3(useSingleJetTriggers_));
  if(!(it == ncontrol_.begin())){
    assert(it != ncontrol_.begin());
    --it;
  }
  return true;
}
 

bool DiJetEventWeighting::passSetWeightDiJetEventWeighting(Event* event)
{
  if(event->type() != PtBalance) std::cout << "Warning: No TwoJetsPtBalanceEvent! ";
  TwoJetsPtBalanceEvent* tje = dynamic_cast<TwoJetsPtBalanceEvent*>(event);
  std::map<double,double>::iterator it = weights_.lower_bound(tje->triggerPtVariableL2L3(useSingleJetTriggers_));
  if(it == weights_.begin()) {
    return false;
  }
  else {
    --it;
    tje->setWeight(tje->weight() * it->second);
    return true;
  }
  
}

