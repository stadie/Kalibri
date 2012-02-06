//
//    Event processor to perform additional (harder) cuts
//    on the preselected dijet events
//
//    first version: Hartmut Stadie 2008/12/14
//   

#include "DiJetEventCuts.h"

#include "CalibData.h"
#include "ConfigFile.h"
#include "Parameters.h"
#include "TwoJetsPtBalanceEvent.h"
#include "CutFlow.h"

#include <iostream>

DiJetEventCuts::DiJetEventCuts(const std::string& configfile, Parameters* param)
  : EventProcessor("DiJetEventCuts",configfile,param)
{  

  //  configfile_=configfile;

}
 
DiJetEventCuts::~DiJetEventCuts()
{
}
  

int DiJetEventCuts::preprocess(std::vector<Event*>& data,
			     std::vector<Event*>& control1,
			     std::vector<Event*>& control2)
{
  std::cout << "processing data: " <<std::endl; 
  CutFlow CutFlow_DATA(this->configName());
  CutFlow_DATA.setAllSuppDiJetCuts();

  //  CutFlow CutFlow_DATA;
  //count events in data and control sample
  for(std::vector<Event*>::iterator i = data.begin() ; i != data.end() ; ++i) {
    if((*i)->type() != PtBalance) continue;
    TwoJetsPtBalanceEvent* tje = dynamic_cast<TwoJetsPtBalanceEvent*>(*i);
    CutFlow_DATA.setDiJetEvent(tje);
    if(CutFlow_DATA.doAllSuppDiJetCuts()==false){
      --i;
      delete tje;
      data.erase(i+1);
      continue;
    }
  }
  CutFlow_DATA.printCutFlow();
  
  std::cout << "processing control1: " <<std::endl; 
  CutFlow CutFlow_CONTROL1(this->configName());
  CutFlow_CONTROL1.setAllSuppDiJetCuts();

  //  CutFlow CutFlow_CONTROL1;
  //count events in data and control sample
  for(std::vector<Event*>::iterator i = control1.begin() ; i != control1.end() ; ++i) {
    if((*i)->type() != PtBalance) continue;
    TwoJetsPtBalanceEvent* tje = dynamic_cast<TwoJetsPtBalanceEvent*>(*i);
    CutFlow_CONTROL1.setDiJetEvent(tje);
    if(CutFlow_CONTROL1.doAllSuppDiJetCuts()==false){
      --i;
      delete tje;
      control1.erase(i+1);
      continue;
    }
  }
  CutFlow_CONTROL1.printCutFlow();
  
  
  return (data.size()+control1.size());
}
 

