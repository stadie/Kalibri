//
//    Base class for event processors
//    
//    first version: Hartmut Stadie 2008/12/14
//    $Id: EventProcessor.cc,v 1.11 2012/01/24 16:25:32 kirschen Exp $
//   
#include "EventProcessor.h"

#include "Parameters.h"
#include "ControlPlots.h"

#include <iostream>

EventProcessor::EventProcessor(const std::string& name, const std::string& configfile, Parameters* param)
  : par_(param), name_(name), configName_(configfile), active_(false)
{  
  ConfigFile config(configfile.c_str());
  config_= config;
  active_ = config.read<bool>(name_,false);
  std::cout << name_.c_str() << " is turned " << (active_ ? "on.\n" : "off.\n"); 
}

EventProcessor::~EventProcessor()
{
}
  

void EventProcessor::produceControlPlots(const std::vector<std::vector<Event*>* >& samples){
  std::cout << name_.c_str() << std::endl;
  // Make control plots
  std::cout << "****Plotting:****\n";
  ControlPlots * plots = new ControlPlots(&config_,samples,this);
  plots->makePlots();
  delete plots;
      
}

