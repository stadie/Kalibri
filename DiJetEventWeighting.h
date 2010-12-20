//
// $Id: DiJetEventWeighting.h,v 1.1 2010/12/13 10:55:09 stadie Exp $
//
#ifndef DIJETEVENTWEIGHTING_H
#define DIJETEVENTWEIGHTING_H

#include "EventProcessor.h"

#include <map>


class Parameters;

// -----------------------------------------------------------------
class DiJetEventWeighting : public EventProcessor
{
 public:
  DiJetEventWeighting(const std::string& configfile, Parameters* param);
  virtual ~DiJetEventWeighting();
  virtual int preprocess(std::vector<Event*>& data,
			 std::vector<Event*>& control1,
			 std::vector<Event*>& control2);
  virtual int postprocess(std::vector<Event*>& data,
			  std::vector<Event*>& control1,
			  std::vector<Event*>& control2) { return data.size();}
  
 private:
  std::map<double,double> ndata_,ncontrol_;
};


#endif
