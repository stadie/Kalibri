//
// $Id: EventBinning.h,v 1.1 2010/12/13 10:55:09 stadie Exp $
//
#ifndef EVENTBINNING_H
#define EVENTBINNING_H

#include "EventProcessor.h"

#include <vector>
#include <map>
#include <set>

class Parameters;
class Binning;
class Event;
class Jet;
// -----------------------------------------------------------------
class EventBinning : public EventProcessor
{
 public:
  EventBinning(const std::string& configfile, Parameters* param);
  virtual ~EventBinning();
  virtual int preprocess(std::vector<Event*>& data,
			 std::vector<Event*>& control1,
			 std::vector<Event*>& control2);
  virtual int postprocess(std::vector<Event*>& data,
			  std::vector<Event*>& control1,
			  std::vector<Event*>& control2);
  
 private:
  typedef std::vector<Event*>::iterator DataIter;
  typedef std::vector<Event*>::const_iterator DataConstIter;
  Event* createBinnedJetTruthEvent(DataIter begin, DataIter end);
  Event* createBinnedJetWidthEvent(DataIter begin, DataIter end);
  bool binEvents_;
  int eventsperbin_;
  Binning* binning_;
  typedef float (Jet::*JetFunction)() const; 
  JetFunction  jetfunc_[4];
  struct ltgenjetpt
  {
    bool operator()(const Event* e1, const Event* e2) const;
  };

  typedef std::map<int,std::vector<Event*> > BinMap;
  typedef std::map<int,std::vector<Event*> >::iterator  BinMapIter;
  BinMap bins_;
};


#endif
