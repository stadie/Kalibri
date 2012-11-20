#ifndef DIJETREADER_H
#define DIJETREADER_H
//!
//!  \brief Reader for dijet events
//!
//!  This class reads dijet events from a ROOT-tree. They are read
//!  by calling the method readEvents(std::vector<Event*>& data).
//!  The data is stored in a format derived from Event depending on the
//!  'Di-Jet data class' field in the config file:
//!    - class 0: Event_PtBalance
//!
//!  It is also possible to store the first two jets of the dijet
//!  event as a JetTruthEvent, where the \f$ p^{\textrm{gen}}_{T} \f$
//!  is used as truth. This is useful for Monte Carlo based calibration:
//!    - class 11: JetTruthEvent of Jet
//!    - class 12: JetTruthEvent of JetWithTowers
//!  In this case, the following cuts as defined in the config file
//!  are applied on each jet in this order:
//!    -# Number of jets greater than 2
//!    -# \f$ p^{\textrm{gen}}_{T} > \texttt{Et genJet min both Jets} \f$
//!    -# \f$ p^{\textrm{gen}}_{T} < \texttt{Et genJet max both Jets} \f$
//!    -# \f$ \Delta R < \texttt{DeltaR cut on jet matching} \f$
//!    -# \f$ p^{\textrm{jet}}_{T} < \texttt{Et cut on jet} \f$
//!    -# \f$ |\eta| > \texttt{Eta cut on jet} \f$
//!    -# \f$ \textrm{Hadronic fraction} < \texttt{Max had fraction} \f$
//!    -# \f$ \textrm{Hadronic fraction} > \texttt{Min had fraction} \f$
//!  
//!
//!  \author Hartmut Stadie
//!  \date 2008/12/12
//!  $Id: DiJetReader.h,v 1.40 2012/11/19 14:27:09 kirschen Exp $
// ----------------------------------------------------------------   



#include "EventReader.h"

#include <string>
#include <memory>
#include <iterator>
#include <vector>

#include "Jet.h"
#include "NJetSel.h"
#include <boost/thread/mutex.hpp>



class NJetSel;
class TRandom;
class JetBin;
class TTree;
class TwoJetsPtBalanceEvent;

class DiJetReader : public EventReader{
 public:
  DiJetReader(const std::string& configfile, Parameters *p);
  virtual ~DiJetReader();
  virtual int readEvents(std::vector<Event*>& data);
  virtual int readControlEvents(std::vector<Event*>& control, int id);
 protected:
  TwoJetsPtBalanceEvent* createTwoJetsPtBalanceEvent();
  void twoJetsPtBalanceSmearJetsJER();
  Event* createDiJetResolutionEvent();
  Event* createDiJetResolutionEventRecoOrdered();
  Event* createDiJetResolutionEventGenOrdered();
  Event* createDiJetResolutionEventFromSkim();
  int createJetTruthEvents(std::vector<Event*>& data);
  CorFactors* createCorFactors(int jetid) const;
  bool passesJetId(int idx) const;
  int readEventsFromTree(std::vector<Event*>& data);
  void printCutFlow();
  double deltaR(double eta1, double eta2, double phi1, double phi2) const;
  double deltaRGenJetCol(int jet1Idx, int jet2Idx) const {
    return deltaR(nJet_->GenJetColEta[jet1Idx],nJet_->GenJetColEta[jet2Idx],nJet_->GenJetColPhi[jet1Idx],nJet_->GenJetColPhi[jet2Idx]);
  }
  bool passesHltChainANDOR() const;
  bool passesHlt(const std::string &trigger) const;
  double eventWeight() const;


  std::auto_ptr<NJetSel> nJet_;                //!< Njet Selector
  TRandom* rand_;             //!< Random number generator

  int    dataClass_;            //!< Data class, see also Event
  int    nDijetEvents_;         //!< Maximum number of read dijet events
  int    prescale_;             //!< only read every nth event
  bool   weights_eq_one_;       //!< force weight for each event to one
  bool   fire_all_dijet_triggers_;       //!< Set all trigger btis to one (branches need to be revised...)
  bool   JERReadInJ1J2SameEtaBin_;       //!< Read in events in same eta bin (not central + any eta) for JER studies; thresholds defined in JEREtaMap_
  std::map<float,int> JEREtaMap_;
  bool   HFAsReferenceRegion_;       //!< Use HF as reference region for TwoJetsPtBalanceEvent  (probe jet in eta>3 )
  bool   correctJECandScaleJER_;     //!< update JEC and scale JER in control events
  bool   readingInControlEvents_;    //!< bool to determine whether this reader is used for reading in cotrol events (default false, altered by readControlEvents(std::vector<Event*>& control, int id))

  double minJetEt_;             //!< Minimum pt of jet
  double maxJetEt_;             //!< Maximum pt of jet
  double minDijetEt_;           //!< Minimum dijet pt
  double maxDijetEt_;           //!< Maximum dijet pt
  double max3rdJetEt_;          //!< Maximum pt of 3rd jet in dijet event
  double minRel3rdJetEt_;       //!< Minimum relative pt of 3rd jet in dijet event
  double maxRel3rdJetEt_;       //!< Maximum relative pt of 3rd jet in dijet event
  double minRelSoftJetEt_;
  double maxRelSoftJetEt_;
  double minDeltaPhi_;          //!< Minimum DeltaPhi for 0 < DeltaPhi < Pi
  double minJetEta_;            //!< Minimum absolute jet eta
  double maxJetEta_;            //!< Maximum absolute jet eta
  double minJetHadFraction_;    //!< Minimum jet Had/(Had+EMF)
  double maxJetHadFraction_;    //!< Maximum jet Had/(Had+EMF)
  double minGenJetEt_;          //!< Minimum pt of genJets of dijets
  double maxGenJetEt_;          //!< Maximum pt of genJets of dijets
  double maxDeltaR_;            //!< Maximum DeltaR
  double minRunNumber_;         //!< Minimum run number allowed for selection
  double maxRunNumber_;         //!< Maximum run number allowed for selection
  std::vector<std::string> hltChainANDOR_;

  int    nReadEvts_;            //!< Number of read events
  int    nGoodEvts_;            //!< Number of events passing all cuts
  int    nDiJetCut_;            //!< Number of events with less than 2 jets
  int    nHlt_;
  int    nVtx_;
  int    nMinDijetEt_;          //!< Number of events rejected by minDijetEt_ cut
  int    nMaxDijetEt_;          //!< Number of events rejected by maxDijetEt_ cut
  int    nMinJetEt_;
  int    nMaxJetEt_;
  int    nCutOn3rdJet_;         //!< Number of events rejected by max3rdJetEt_ or maxRelJetEt_ cut
  int    nCutOnSoftJets_;
  int    nMinDeltaPhi_;         //!< Number of events rejected by maxDeltaPhi_ cut
  int    nMinJetEta_;           //!< Number of events rejected by minJetEta_ cut
  int    nMaxJetEta_;           //!< Number of events rejected by maxJetEta_ cut
  int    nMaxGenJetEt_;         //!< Number of events rejected by maxGenJetEt_ cut
  int    nMinGenJetEt_;         //!< Number of events rejected by minGenJetEt_ cut
  int    nJetIDCut_;
  int    nMinJetHadFraction_;   //!< Number of events rejected by minJetHadFraction_ cut
  int    nMaxJetHadFraction_;   //!< Number of events rejected by maxJetHadFraction_ cut
  int    nMaxDeltaR_;           //!< Number of events rejected by maxDeltaR_ cut
  int    nCutOnMinRunNumber_;   //!< Number of events rejected by minRunNumber_ cut
  int    nCutOnMaxRunNumber_;   //!< Number of events rejected by maxRunNumber_ cut
  int    nTriggerSel_;          //!< Number of events not passing the trigger selection
  int    nMaxMCPU_;             //!< Max number of mixed-in PU events
  int    nMaxVtxN_;             //!< Max number of reconstructed vertices (easy PU-cut)
  int    nMinVtxN_;             //!< Min number of reconstructed vertices (easy PU-cut)
  int    maxNIter_;             //!< Max number of iterations in integration
  double eps_;                  //!< Integration precision for convergence
  double truthSpecExp_;         //!< Exponent of truth spectrum
  double genjetpt_,genjete_,jeteta_,sigmaphi_,sigmaeta_,sumsigmaetaphi_,emf_,meanMoment_; //!< possible binning variables

  int minJetN90Hits_;
  double maxJetFHPD_;
  double maxJetFRBX_;

  const double* vars_[4];             //!< Jet binning variables
  const double zero_;           //!< just null

  std::vector<Jet::JetIndex*> jetIndices_;

  friend class ThreadedDiJetReader; 
  static boost::mutex dijetmutex_;

  //for DiJetAve trigger
  bool hltdijetave15incl_,hltdijetave30incl_,hltdijetave50incl_,hltdijetave70incl_,hltdijetave100incl_,hltdijetave140incl_,hltdijetave180incl_,hltdijetave300incl_;
  bool hltdijetavec30incl_,hltdijetavec60incl_,hltdijetavec80incl_,hltdijetavec110incl_,hltdijetavec150incl_,hltdijetavec190incl_,hltdijetavec240incl_,hltdijetavec300incl_,hltdijetavec370incl_;

  //for SingleJet trigger
  bool useSingleJetTriggers_;
  bool hltjetc30incl_,hltjetc60incl_,hltjetc80incl_,hltjetc110incl_,hltjetc150incl_,hltjetc190incl_,hltjetc240incl_,hltjetc300incl_,hltjetc370incl_;

  //for 2012 PF trigger 
  bool hltdiPFjetc40incl_,hltdiPFjetc80incl_,hltdiPFjetc140incl_,hltdiPFjetc200incl_,hltdiPFjetc260incl_,hltdiPFjetc320incl_,hltdiPFjetc400incl_,hltPFjetc40incl_,hltPFjetc80incl_,hltPFjetc140incl_,hltPFjetc200incl_,hltPFjetc260incl_,hltPFjetc320incl_,hltPFjetc400incl_;


  std::map<double,bool*> trigmap_;
  std::map<double,double> mcweightmap_[2];
  bool requireTrigger_, mcweightmapid_;
};


#endif
