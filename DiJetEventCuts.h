#ifndef DIJETEVENTCUTS_H
#define DIJETEVENTCUTS_H

#include "EventProcessor.h"

#include <map>


class Parameters;

//! \brief Event processor to perform additional (harder) cuts on the preselected dijet events
//! Prints separate cut-flow and reads in harder cuts from config file
//! 
//! 
// -----------------------------------------------------------------
class DiJetEventCuts : public EventProcessor
{
 public:
  DiJetEventCuts(const std::string& configfile, Parameters* param);
  virtual ~DiJetEventCuts();
protected:
  virtual int preprocess(std::vector<Event*>& data,
			 std::vector<Event*>& control1,
			 std::vector<Event*>& control2);
  virtual int postprocess(std::vector<Event*>& data,
			  std::vector<Event*>& control1,
			  std::vector<Event*>& control2) { return data.size();}
  
 private:
  //  std::string& configfile_;
};


#endif
