//
// $Id: EventProcessor.h,v 1.7 2010/12/20 11:08:13 stadie Exp $
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
  EventProcessor(const std::string& configfile, Parameters* param);
  virtual ~EventProcessor();
  virtual int preprocess(std::vector<Event*>& data,
			 std::vector<Event*>& control1,
			 std::vector<Event*>& control2) = 0;
  virtual int postprocess(std::vector<Event*>& data,
			  std::vector<Event*>& control1,
			  std::vector<Event*>& control2) = 0;

 private:
  Parameters* par_;
};


#endif
