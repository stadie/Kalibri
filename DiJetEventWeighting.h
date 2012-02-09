//
// $Id: DiJetEventWeighting.h,v 1.2 2011/06/30 14:27:14 stadie Exp $
//
#ifndef DIJETEVENTWEIGHTING_H
#define DIJETEVENTWEIGHTING_H

#include "EventProcessor.h"
#include "TwoJetsPtBalanceEvent.h"
#include "CalibData.h"

#include <map>


class Parameters;

// -----------------------------------------------------------------
class DiJetEventWeighting : public EventProcessor
{
 public:
  DiJetEventWeighting(const std::string& configfile, Parameters* param);
  virtual ~DiJetEventWeighting();
  bool passCheckBadEventDiJetEventWeighting(Event* event);
  bool passCheckStrangeEventDiJetEventWeighting(Event* event);
  bool passSetWeightDiJetEventWeighting(Event* event);

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
};

class CheckBadEventDiJetEventWeighting{
public:
  CheckBadEventDiJetEventWeighting(DiJetEventWeighting* DJEWptr) : DJEWptr_(DJEWptr){};
  bool operator() (Event* event)
  {
    return DJEWptr_->passCheckBadEventDiJetEventWeighting(event);
  }
private:
  DiJetEventWeighting* DJEWptr_;
  
};

class CheckStrangeEventDiJetEventWeighting{
public:
  CheckStrangeEventDiJetEventWeighting(DiJetEventWeighting* DJEWptr) : DJEWptr_(DJEWptr){};
  bool operator() (Event* event)
  {
    return DJEWptr_->passCheckStrangeEventDiJetEventWeighting(event);
  }
private:
  DiJetEventWeighting* DJEWptr_;
  
};

class SetWeightDiJetEventWeighting{
public:
  SetWeightDiJetEventWeighting(DiJetEventWeighting* DJEWptr) : DJEWptr_(DJEWptr){};
  bool operator() (Event* event)
  {
    return DJEWptr_->passSetWeightDiJetEventWeighting(event);
  }
private:
  DiJetEventWeighting* DJEWptr_;
  
};


#endif
