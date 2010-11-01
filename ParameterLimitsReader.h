//
//    Reader for Parameter Limits
//
//    This class add user defined parameter limits
//
//    first version: Hartmut Stadie 2008/12/12
//    $Id: ParameterLimitsReader.h,v 1.3 2010/10/20 11:28:19 stadie Exp $
//   
#ifndef PARAMETERLIMITSREADER_H
#define PARAMETERLIMITSREADER_H

#include "EventReader.h"

#include <string>

class ParameterLimitsReader : public EventReader{
 public:
  ParameterLimitsReader(const std::string& configfile, Parameters *p);
  virtual ~ParameterLimitsReader();
  int readEvents(std::vector<Event*>& data);
 private:
  class ParLimit {
  public:
    int index;
    double min;
    double max;
    double k;
    ParLimit(int index, double min, double max, double k) 
      : index(index), min(min), max(max), k(k) {}
  };
  std::vector<ParLimit> par_limits;
};


#endif
