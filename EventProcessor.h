//
// $Id: EventProcessor.h,v 1.8 2011/04/28 16:05:34 stadie Exp $
//
#ifndef EVENTPROCESSOR_H
#define EVENTPROCESSOR_H

#include <vector>
#include <string>

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
  bool active_;
};


#endif
