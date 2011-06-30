//
//    Base class for event processors
//    
//    first version: Hartmut Stadie 2008/12/14
//    $Id: EventProcessor.cc,v 1.9 2011/06/06 14:58:58 mschrode Exp $
//   
#include "EventProcessor.h"

#include "Parameters.h"
#include "ConfigFile.h"

#include <iostream>

EventProcessor::EventProcessor(const std::string& name, const std::string& configfile, Parameters* param)
  : par_(param), name_(name), active_(false)
{  
  ConfigFile config(configfile.c_str());
  active_ = config.read<bool>(name_,false);
  std::cout << name_.c_str() << " is turned " << (active_ ? "on.\n" : "off.\n"); 
}

EventProcessor::~EventProcessor()
{
}
  



