#ifndef CUTFLOW_H
#define CUTFLOW_H

#include "ConfigFile.h"
#include "TwoJetsPtBalanceEvent.h"
#include "progressbar.h"

class Event;

//!  \brief Class for performing cuts on DiJetEvents
//!         and keeping track of them
//!
//!  \author Henning Kirschenmann
//!  \date 2012/01/27
// ----------------------------------------------------------------   
class CutFlow{
 public:
  //  std::map<>
  //CutFlow();
  //CutFlow(const ConfigFile *configFile);
  CutFlow(const std::string& configfile);
  virtual ~CutFlow();

  void setEvent(Event* event);
  void setDiJetEvent(Event* event);
  bool doMaxEtaCut();
  bool doAllSuppDiJetCuts();
  bool survivesAllSuppDiJetCuts(Event* event);
  void setNExpectedEvents(int nExpectedEvents) {nExpectedEvents_=nExpectedEvents;}
  void setAllSuppDiJetCuts();
  void printCutFlow();
  bool operator() (Event* event)
    {
      //      if(nExpectedEvents_-nEvents_ % 10 == 0)progressbar(nEvents_*100/nExpectedEvents_);
      return survivesAllSuppDiJetCuts(event);
    }

 private:
  const ConfigFile *config_;
  int nExpectedEvents_;
  int nEvents_;
  int nMaxEtaCut_;
  Event* event_;
  TwoJetsPtBalanceEvent* diJetEvent_;
  int    nMinDeltaPhi_;         //!< Number of events rejected by maxDeltaPhi_ cut
  int    nMinCutOn3rdJet_;         //!< Number of events rejected by minRel3rdJetEt_ cut
  int    nMaxCutOn3rdJet_;         //!< Number of events rejected by maxRel3rdJetEt_ cut
  int    nMaxCutOnAsymmetry_;      //!< Number of events rejected by cut on asymmetry
  double minDeltaPhi_;          //!< Minimum DeltaPhi for 0 < DeltaPhi < Pi
  double maxRel3rdJetEt_;       //!< Maximum relative pt of 3rd jet in dijet event
  double minRel3rdJetEt_;       //!< Minimum relative pt of 3rd jet in dijet event
  double maxCutOnAsymmetry_;       //!< Maximum asymmetry to consider event
  bool useMinDeltaPhi_;         //!< Apply cut on Minimum DeltaPhi for 0 < DeltaPhi < Pi
  bool useMinRel3rdJetEt_;      //!< Apply cut on Minimum relative pt of 3rd jet in dijet event
  bool useMaxRel3rdJetEt_;      //!< Apply cut on Maximum relative pt of 3rd jet in dijet event
  bool useMaxCutOnAsymmetry_;      //!< Apply cut on Maximum asymmetry in dijet event

};


class CheckDiJetCuts{
public:
  CheckDiJetCuts(CutFlow* cutflowptr) : cutflowptr_(cutflowptr){};
  bool operator() (Event* event)
  {
    return (*cutflowptr_)(event);
  }
private:
  CutFlow* cutflowptr_;
  
};








#endif
