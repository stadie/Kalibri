//
// $Id: EventProcessor.h,v 1.5 2010/10/20 11:28:17 stadie Exp $
//
#ifndef EVENTPROCESSOR_H
#define EVENTPROCESSOR_H

#include <vector>

class Event;
class Parameters;

// -----------------------------------------------------------------
class EventProcessor
{
 public:
  EventProcessor(const std::string& configfile, Parameters* param);
  virtual ~EventProcessor();
  virtual int preprocess(std::vector<Event*>& data) = 0;
  virtual int postprocess(std::vector<Event*>& data) = 0;

 private:
  Parameters* par_;
};


#endif
