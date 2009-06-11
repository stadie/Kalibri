#ifndef PHOTONJETREADER_H
#define PHOTONJETREADER_H

#include "EventReader.h"

#include <string>

#include "GammaJetSel.h"

//!
//!  \brief Reader for Photon Jet Events
//!
//!  This class reads events according to the GammaJetSel.
//!
//!  \author Hartmut Stadie
//!  \date 2008/12/12
//!  $Id: PhotonJetReader.h,v 1.4 2009/04/17 14:28:08 mschrode Exp $
// ----------------------------------------------------------------   
class PhotonJetReader : public EventReader{
 public:
  PhotonJetReader(const std::string& configfile, TParameters *p);
  virtual ~PhotonJetReader() {};
  int readEvents(std::vector<TData*>& data);

 private:
  TData* createTruthMultMessEvent();
  TData* createJetTruthEvent();
  TData* createSmearEvent();


  GammaJetSel gammajet;        //!< Gamma-jet Selector

  int    dataClass;            //!< Data class, see also TData
  int    n_gammajet_events;    //!< Maximum number of read photon jet events

  double Et_cut_on_jet;        //!< Minimum pt of jet
  double Et_cut_on_gamma;      //!< Maximum pt of photon
  double Rel_cut_on_gamma;     //!< Maximum relative pt of non-leading jet
  double GenJetCutLow;         //!< Minimum pt of genJet
  double GenJetCutUp;          //!< Maximum pt of genJet
  double DeltaRMatchingCut;    //!< Maximum DeltaR
  double Eta_cut_on_jet;       //!< Maximum absolute jet eta
  double Had_cut_min;          //!< Minimum jet Had/(Had+EMF)
  double Had_cut_max;          //!< Maximum jet Had/(Had+EMF)

  int    nEt_cut_on_jet;       //!< Number of events rejected by Et_cut_on_jet
  int    nEt_cut_on_gamma;     //!< Number of events rejected by Et_cut_on_gamma cut
  int    nRel_cut_on_gamma;    //!< Number of events rejected by nRel_cut_on_gamma cut
  int    nGenJetCutLow;        //!< Number of events rejected by GenJetCutLow cut
  int    nGenJetCutUp;         //!< Number of events rejected by GenJetCutUp
  int    nDeltaRMatchingCut;   //!< Number of events rejected by DeltaRMatchingCut
  int    nEta_cut_on_jet;      //!< Number of events rejected by Eta_cut_on_jet
  int    nHad_cut_min;         //!< Number of events rejected by Had_cut_min
  int    nHad_cut_max;         //!< Number of events rejected by Had_cut_max
};


#endif
