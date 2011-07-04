#ifndef TWO_JETS_PT_BALANCE_BINNED_EVENT_H
#define TWO_JETS_PT_BALANCE_BINNED_EVENT_H
//!
//!  \brief Class for relative calibration in pseudorapidity
//!         using binned dijet events
//!
//!  \author Hartmut Stadie
//!  \date Mon Oct 26 21:03:43 CET 2009 
//!  $Id: TwoJetsPtBalanceEvent.h,v 1.12 2010/11/01 15:47:40 stadie Exp $
// --------------------------------------------------


#include <cmath>

#include "TwoJetsPtBalanceEvent.h"


class TwoJetsPtBalanceBinnedEvent : public TwoJetsPtBalanceEvent {
 public:
  TwoJetsPtBalanceBinnedEvent(double maxAlpha); 
  ~TwoJetsPtBalanceBinnedEvent();

  int addEvent(TwoJetsPtBalanceEvent* tjpbe);
  unsigned int size() const { return events_.size();}
  TwoJetsPtBalanceBinnedEvent* split(double pt);
  void computeJets();
  // from TwoJetsPtBalanceEvent
  double ptMedian();
  double relPtJet3() const { return maxAlpha_ - 0.000001;}
 private:
  //should be a list sorted in pt....for ptMedian
  std::vector<TwoJetsPtBalanceEvent*> events_;
  double maxAlpha_;
  struct lessPtAve
  {
    bool operator()(const TwoJetsPtBalanceEvent* e1, const TwoJetsPtBalanceEvent* e2) const {
      return e1->ptDijet() < e2->ptDijet();
    }
  };
  struct lessThanPtAve
  {
  private:
    double pt_;
  public:
    lessThanPtAve(double pt) : pt_(pt) {}
    bool operator()(const TwoJetsPtBalanceEvent* e1) const {
      return e1->ptDijet() < pt_;
    }
  };
};

#endif
