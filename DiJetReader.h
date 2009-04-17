#ifndef DIJETREADER_H
#define DIJETREADER_H

#include "EventReader.h"

#include <string>

#include "NJetSel.h"


//!
//!    Reader for Di-Jet Events
//!
//!    This class reads events according to the NJetSel
//!
//!    \author Hartmut Stadie
//!    \date 2008/12/12
//!    $Id: DiJetReader.h,v 1.3 2009/03/04 17:26:51 thomsen Exp $
// ----------------------------------------------------------------   
class DiJetReader : public EventReader{
 public:
  DiJetReader(const std::string& configfile, TParameters *p);
  virtual ~DiJetReader();
  int readEvents(std::vector<TData*>& data);

 private:
  TData* createPtBalanceEvent();
  int createJetTruthEvents(std::vector<TData*>& data);

  NJetSel njet;
  double Et_cut_nplus1Jet;     //!< Maximum pt of other than leading 2 jets in dijet event
  double Rel_cut_on_nJet;      //!< Maximum relative pt  of other than leading 2 jets in dijet event
  double GenJetCut;            //!< Minimum pt of genJets of dijets
  double Eta_cut_on_jet;       //!< Maximum absolute jet eta
  double Had_cut_min;          //!< Minimum jet Had/(Had+EMF)
  double Had_cut_max;          //!< Maximum jet Had/(Had+EMF)
  int    n_dijet_events;       //!< Maximum number of read dijet events
  int    dataClass;            //!< Data class, see also TData
};


#endif
