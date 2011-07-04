//
// $Id: DiJetEventWeighting.h,v 1.1 2010/12/20 11:08:13 stadie Exp $
//
#ifndef DIJETEVENTBINNING_H
#define DIJETEVENTBINNING_H

#include "EventProcessor.h"

#include <map>


class Parameters;
class TwoJetsPtBalanceEvent;
class TwoJetsPtBalanceBinnedEvent;

// -----------------------------------------------------------------
class DiJetEventBinning : public EventProcessor
{
public:
  DiJetEventBinning(const std::string& configfile, Parameters* param);
  virtual ~DiJetEventBinning();
protected:
  virtual int preprocess(std::vector<Event*>& data,
			 std::vector<Event*>& control1,
			 std::vector<Event*>& control2);
  virtual int postprocess(std::vector<Event*>& data,
			  std::vector<Event*>& control1,
			  std::vector<Event*>& control2) { return data.size();}
  
private: 
  class EtaPtKey {
  private:
    const double mineta_,minpt_;
  public:
    EtaPtKey(double mineta = 0, double minpt = 0) : mineta_(mineta),minpt_(minpt) {} 
    double eta() const {return mineta_;}
    double pt() const {return minpt_;}    
    bool operator<(const EtaPtKey& k) const {
      if(eta() < k.eta()) return true;
      if(eta() > k.eta()) return false;
      if(pt() < k.pt()) return true;
      return false;
    }
  };
  typedef std::map<EtaPtKey,TwoJetsPtBalanceBinnedEvent*> EventMap;
  EventMap::const_iterator getIterator(double eta, double pt, 
				       const EventMap& map) const;
  TwoJetsPtBalanceBinnedEvent* getEvent(double eta, double pt, 
					const EventMap& map) const {
    EventMap::const_iterator i = getIterator(eta,pt,map);
    if(i != map.end()) return i->second;
    return 0;
  }
  
  void split(double eta, double pt);

  static const int Nalpha_ = 4;
  unsigned int minEvents_;
  std::vector<double> etavalues_;
  EventMap  databins_[Nalpha_],controlbins_[Nalpha_];
};


#endif
