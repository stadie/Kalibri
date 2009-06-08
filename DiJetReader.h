#ifndef DIJETREADER_H
#define DIJETREADER_H

#include "EventReader.h"

#include <string>

#include "NJetSel.h"


//!
//!  \brief Reader for dijet events
//!
//!  This class reads dijet events from a ROOT-tree. They are read
//!  by calling the method readEvents(std::vector<TData*>& data).
//!  The data is stored in a format derived from TData depending on the
//!  'Di-Jet data class' field in the config file:
//!    - class 0: TData_PtBalance
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
//!    -# \f$ p^{\textrm{jet}}_{T} < \texttt{Et cut on jet}
//!    -# \f$ |\eta| > \texttt{Eta cut on jet} \f$
//!    -# \f$ \textrm{Hadronic fraction} < \texttt{Max had fraction} \f$
//!    -# \f$ \textrm{Hadronic fraction} > \texttt{Min had fraction} \f$
//!  
//!
//!  \author Hartmut Stadie
//!  \date 2008/12/12
//!  $Id: DiJetReader.h,v 1.6 2009/06/05 15:44:20 mschrode Exp $
// ----------------------------------------------------------------   
class DiJetReader : public EventReader{
 public:
  DiJetReader(const std::string& configfile, TParameters *p);
  virtual ~DiJetReader();
  int readEvents(std::vector<TData*>& data);

 private:
  TData* createPtBalanceEvent();
  int createJetTruthEvents(std::vector<TData*>& data);

  NJetSel njet;                //!< Njet Selector

  int    dataClass;            //!< Data class, see also TData
  int    n_dijet_events;       //!< Maximum number of read dijet events

  double Et_cut_on_jet;        //!< Minimum pt of jet
  double Et_cut_nplus1Jet;     //!< Maximum pt of other than leading 2 jets in dijet event
  double Rel_cut_on_nJet;      //!< Maximum relative pt  of other than leading 2 jets in dijet event
  double GenJetCutLow;         //!< Minimum pt of genJets of dijets
  double GenJetCutUp;          //!< Maximum pt of genJets of dijets
  double DeltaRMatchingCut;    //!< Maximum DeltaR
  double Eta_cut_on_jet;       //!< Maximum absolute jet eta
  double Had_cut_min;          //!< Minimum jet Had/(Had+EMF)
  double Had_cut_max;          //!< Maximum jet Had/(Had+EMF)

  int    nNjet_cut;            //!< Number of events with less than 2 jets
  int    nEt_cut_on_jet;       //!< Number of events rejected by Et_cut_on_jet
  int    nEt_cut_nplus1Jet;    //!< Number of events rejected by Et_cut_nplus1Jet cut
  int    nRel_cut_on_nJet;     //!< Number of events rejected by Rel_cut_on_nJet cut
  int    nGenJetCutLow;        //!< Number of events rejected by GenJetCutLow cut
  int    nGenJetCutUp;         //!< Number of events rejected by GenJetCutUp
  int    nDeltaRMatchingCut;   //!< Number of events rejected by DeltaRMatchingCut
  int    nEta_cut_on_jet;      //!< Number of events rejected by Eta_cut_on_jet
  int    nHad_cut_min;         //!< Number of events rejected by Had_cut_min
  int    nHad_cut_max;         //!< Number of events rejected by Had_cut_max
};


#endif
