//
// $Id: EventProcessor.h,v 1.10 2012/01/24 16:25:33 kirschen Exp $
//
#ifndef EVENTPROCESSOR_H
#define EVENTPROCESSOR_H

#include <vector>
#include <string>
#include "ConfigFile.h"

class Event;
class Parameters;

// -----------------------------------------------------------------
class EventProcessor
{
 public:
  EventProcessor(const std::string& name, const std::string& configfile, 
		 Parameters* param);
  virtual ~EventProcessor();
  int process(std::vector<Event*>& data,
	      std::vector<Event*>& control1,
	      std::vector<Event*>& control2) {
    if(active_) return preprocess(data,control1,control2);
    return data.size();
  }
  int revert(std::vector<Event*>& data,
	     std::vector<Event*>& control1,
	     std::vector<Event*>& control2) {
    if(active_) return postprocess(data,control1,control2);
    return data.size();
  }
  const std::string& name() const { return name_;}
  const std::string& configName() const { return configName_;}
  void produceControlPlots(const std::vector<std::vector<Event*>* >& samples);
  
protected:
  virtual int preprocess(std::vector<Event*>& data,
			 std::vector<Event*>& control1,
			 std::vector<Event*>& control2) = 0;
  virtual int postprocess(std::vector<Event*>& data,
			  std::vector<Event*>& control1,
			  std::vector<Event*>& control2) = 0;

 private:
  Parameters* par_;
  const std::string name_;
  const std::string configName_;
  bool active_;
  ConfigFile config_;
};


#endif
