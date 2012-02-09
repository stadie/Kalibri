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
#include "progressbar.h"

#include <algorithm>
#include <iostream>

DiJetEventCuts::DiJetEventCuts(const std::string& configfile, Parameters* param)
  : EventProcessor("DiJetEventCuts",configfile,param)
{  

  //  configfile_=configfile;

}
 
DiJetEventCuts::~DiJetEventCuts()
{
}
  
typedef  bool (CutFlow::*CutFlowMemFn)(Event* event);

int DiJetEventCuts::preprocess(std::vector<Event*>& data,
			     std::vector<Event*>& control1,
			     std::vector<Event*>& control2)
{
  std::cout << "processing data: " <<std::endl; 
  CutFlow CutFlow_DATA(this->configName());
  CutFlow_DATA.setAllSuppDiJetCuts();
  CutFlow_DATA.setNExpectedEvents(data.size());
  int l = data.size()/100;

  std::cout << "start with "  << data.size() << " events... check with cutflow:"<< std::endl;
 
  CheckDiJetCuts passCutsPredicate_DATA(&CutFlow_DATA);
  std::vector<Event*>::iterator bound;
  bound= partition(data.begin(),data.end(),passCutsPredicate_DATA);
  data.erase(bound,data.end());
  
  std::cout << "kept "  << data.size() << " events... check with cutflow:"<< std::endl;
  CutFlow_DATA.printCutFlow();
  
  std::cout << "processing control1: " <<std::endl; 
  CutFlow CutFlow_CONTROL1(this->configName());
  CutFlow_CONTROL1.setAllSuppDiJetCuts();
  CutFlow_CONTROL1.setNExpectedEvents(control1.size());
  CheckDiJetCuts passCutsPredicate_CONTROL1(&CutFlow_CONTROL1);

  bound= partition(control1.begin(),control1.end(),passCutsPredicate_CONTROL1);
  control1.erase(bound,control1.end());
  CutFlow_CONTROL1.printCutFlow();
  std::cout << "kept "  << control1.size() << " events... check with cutflow:"<< std::endl;
  
  return (data.size()+control1.size());
}
 

