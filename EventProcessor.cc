//
//    Class for existing event processors
//    
//    Right now we only have the spectrum correctors.
//    Thus they are implemented directly in this class
//
//    first version: Hartmut Stadie 2008/12/14
//    $Id: EventProcessor.cc,v 1.8 2010/12/13 10:55:09 stadie Exp $
//   
#include "EventProcessor.h"

#include "Parameters.h"



EventProcessor::EventProcessor(const std::string& configfile, Parameters* param)
  : par_(param)
{  
}
 
EventProcessor::~EventProcessor()
{
}
  



