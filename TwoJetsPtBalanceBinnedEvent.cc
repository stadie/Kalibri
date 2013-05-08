// $Id: TwoJetsPtBalanceBinnedEvent.cc,v 1.4 2012/04/27 12:50:13 kirschen Exp $

#include "TwoJetsPtBalanceBinnedEvent.h"
#include "JetBin.h"

#include <algorithm>
#include <functional>

TwoJetsPtBalanceBinnedEvent::TwoJetsPtBalanceBinnedEvent(double maxAlpha) :  
  TwoJetsPtBalanceEvent(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0), maxAlpha_(maxAlpha)
{
}

TwoJetsPtBalanceBinnedEvent::~TwoJetsPtBalanceBinnedEvent()
{
}

void TwoJetsPtBalanceBinnedEvent::computeJets()
{
  JetBin* jetbins[3];
  TwoJetsPtBalanceEvent* tjpbe = events_.front();
  Jet* jet = tjpbe->jet1();
  jetbins[0] = new JetBin(*jet->f(), jet->errFunc(), *jet->gf());
  jet = tjpbe->jet2();
  jetbins[1] = new JetBin(*jet->f(), jet->errFunc(), *jet->gf());
  jet = tjpbe->jet1(); 
  jetbins[2] = new JetBin(*jet->f(), jet->errFunc(), *jet->gf());
  for(std::vector<TwoJetsPtBalanceEvent*>::const_iterator i = events_.begin();
      i != events_.end() ; ++i) {
    TwoJetsPtBalanceEvent* tjpbe = *i;
    jetbins[0]->addJet(tjpbe->jet1());
    jetbins[1]->addJet(tjpbe->jet2());
    if(tjpbe->hasJet3()) {
      jetbins[2]->addJet(tjpbe->jet3());
    }
  }
  delete jet1_;
  jet1_ = jetbins[0]->jet();
  delete jet2_;
  jet2_ = jetbins[1]->jet();
  delete jet3_;
  jet3_ = jetbins[2]->jet();

  delete jetbins[0];
  delete jetbins[1];
  delete jetbins[2];
}




int TwoJetsPtBalanceBinnedEvent::addEvent(TwoJetsPtBalanceEvent* tjpbe)
{
  events_.push_back(tjpbe);
  return events_.size();
}
double TwoJetsPtBalanceBinnedEvent::ptMedian()
{
  assert(! (events_.size() % 2));
  std::vector<TwoJetsPtBalanceEvent*>::iterator middle = events_.begin() + events_.size()/2;
  std::nth_element(events_.begin(),middle,events_.end(),lessPtAve());
  return (*middle)->ptDijet();
}

// returns the newly created Event
TwoJetsPtBalanceBinnedEvent* TwoJetsPtBalanceBinnedEvent::split(double pt)
{
  TwoJetsPtBalanceBinnedEvent* e2 = new TwoJetsPtBalanceBinnedEvent(maxAlpha_);
  //rearrange events so that all events before middle have less pt
  std::vector<TwoJetsPtBalanceEvent*>::iterator middle = 
    std::partition(events_.begin(),events_.end(),lessThanPtAve(pt));
  for(std::vector<TwoJetsPtBalanceEvent*>::const_iterator i = middle;
      i != events_.end() ; ++i) {
    e2->addEvent(*i);
  }
  events_.erase(middle,events_.end());
  return e2;
}
