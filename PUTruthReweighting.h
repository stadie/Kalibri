//
#ifndef PUTRUTHREWEIGHTING_H
#define PUTRUTHREWEIGHTING_H

#include "EventProcessor.h"
#include "TwoJetsPtBalanceEvent.h"
#include "CalibData.h"

#include <map>


class Parameters;

// -----------------------------------------------------------------
class PUTruthReweighting : public EventProcessor
{
 public:
  PUTruthReweighting(const std::string& configfile, Parameters* param);
  virtual ~PUTruthReweighting();
  double TriggerPtVariable(Event* event);

protected:
  virtual int preprocess(std::vector<Event*>& data,
			 std::vector<Event*>& control1,
			 std::vector<Event*>& control2);
  virtual int postprocess(std::vector<Event*>& data,
			  std::vector<Event*>& control1,
			  std::vector<Event*>& control2) { return data.size();}
  
 private:
  std::map<double,double> ndata_,ncontrol_;
  std::map<double,double> weights_;
  bool useSingleJetTriggers_;
  std::vector<std::string> trignames_;
  std::vector<double> trigthresholds_;

};



#endif
