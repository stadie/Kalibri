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
//!  $Id: DiJetReader.h,v 1.17 2010/04/13 13:44:09 mschrode Exp $
// ----------------------------------------------------------------   



#include "EventReader.h"

#include <string>
#include <memory>
#include <map>
#include <set>
#include <iterator>

class Jet;
class NJetSel;
class TRandom;
class JetBin;

class DiJetReader : public EventReader{
 public:
  DiJetReader(const std::string& configfile, TParameters *p);
  ~DiJetReader();
  int readEvents(std::vector<Event*>& data);


 private:
  Event* createTwoJetsPtBalanceEvent();
  Event* createSmearEvent(int callIdx = 0);
  int createJetTruthEvents(std::vector<Event*>& data);
  CorFactors* createCorFactors(int jetid) const;
  std::vector<Jet*> readCaloJets(int nJets) const;
  std::vector<Jet*> readGenJetSortedJets(int nJets) const;

  std::auto_ptr<NJetSel> nJet_;                //!< Njet Selector
  TRandom* rand_;             //!< Random number generator

  int    dataClass_;            //!< Data class, see also Event
  int    nDijetEvents_;         //!< Maximum number of read dijet events
  int    prescale_;             //!< only read every nth event

  double minJetEt_;             //!< Minimum pt of jet
  double maxJetEt_;             //!< Maximum pt of jet
  double minDijetEt_;           //!< Minimum dijet pt
  double maxDijetEt_;           //!< Maximum dijet pt
  double max3rdJetEt_;          //!< Maximum pt of 3rd jet in dijet event
  double maxRel3rdJetEt_;       //!< Maximum relative pt of 3rd jet in dijet event
  double minDeltaPhi_;          //!< Minimum DeltaPhi for 0 < DeltaPhi < Pi
  double maxJetEta_;            //!< Maximum absolute jet eta
  double minJetHadFraction_;    //!< Minimum jet Had/(Had+EMF)
  double maxJetHadFraction_;    //!< Maximum jet Had/(Had+EMF)
  double minGenJetEt_;          //!< Minimum pt of genJets of dijets
  double maxGenJetEt_;          //!< Maximum pt of genJets of dijets
  double maxDeltaR_;            //!< Maximum DeltaR

  int    nDiJetCut_;            //!< Number of events with less than 2 jets
  int    nMinJetEt_;            //!< Number of events rejected by minJetEt_ cut
  int    nMaxJetEt_;            //!< Number of events rejected by maxJetEt_ cut
  int    nMinDijetEt_;          //!< Number of events rejected by minDijetEt_ cut
  int    nMaxDijetEt_;          //!< Number of events rejected by maxDijetEt_ cut
  int    nCutOn3rdJet_;         //!< Number of events rejected by max3rdJetEt_ or maxRelJetEt_ cut
  int    nMinDeltaPhi_;         //!< Number of events rejected by maxDeltaPhi_ cut
  int    nMaxJetEta_;           //!< Number of events rejected by maxJetEta_ cut
  int    nMinJetHadFraction_;   //!< Number of events rejected by minJetHadFraction_ cut
  int    nMaxJetHadFraction_;   //!< Number of events rejected by maxJetHadFraction_ cut
  int    nMaxGenJetEt_;         //!< Number of events rejected by maxGenJetEt_ cut
  int    nMinGenJetEt_;         //!< Number of events rejected by minGenJetEt_ cut
  int    nMaxDeltaR_;           //!< Number of events rejected by maxDeltaR_ cut

  int    maxNIter_;             //!< Max number of iterations in integration
  double eps_;                  //!< Integration precision for convergence
  double min_;                  //!< Minimum of truth spectrum in integration
  double max_;                  //!< Maximum of truth spectrum in integration
  double truthSpecExp_;         //!< Exponent of truth spectrum

  //handle binned events
  int findBin(double eta, double pt) {
    //from CMSSW/JetMETCorrections/MCJet/test/Settings.h
    static const double Pt[33] = {5,10,12,15,18,22,26,30,35,40,45,51,57,64,72,80,90,105,120,135,150,175,200,250,300,350,400,500,650,800,1000,1500,5000};
    //static const double Pt[1] = {7000};
    static const double eta_boundaries[83] = {-5.191,-4.889,-4.716,-4.538,-4.363,-4.191,-4.013,-3.839,-3.664,-3.489,-3.314,-3.139,-2.964,-2.853,-2.650,-2.500,-2.322,-2.172,-2.043,-1.930,-1.830,-1.740,-1.653,-1.566,-1.479,-1.392,-1.305,-1.218,-1.131,-1.044,-0.957,-0.879,-0.783,-0.696,-0.609,-0.522,-0.435,-0.348,-0.261,-0.174,-0.087,0.000,0.087,0.174,0.261,0.348,0.435,0.522,0.609,0.696,0.783,0.879,0.957,1.044,1.131,1.218,1.305,1.392,1.479,1.566,1.653,1.740,1.830,1.930,2.043,2.172,2.322,2.500,2.650,2.853,2.964,3.139,3.314,3.489,3.664,3.839,4.013,4.191,4.363,4.538,4.716,4.889,5.191};
    static const std::set<double> ptbins(Pt,Pt+33);
    static const std::set<double> etabins(eta_boundaries,eta_boundaries+83);
    int ipt = std::distance(ptbins.begin(),ptbins.lower_bound(pt));
    int ieta = std::distance(etabins.begin(),etabins.lower_bound(eta));
    return ipt * 100 + ieta;
  }
  std::map<int,JetBin*> jetbins_;

};


#endif
