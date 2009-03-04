//
//    Reader for Photon Jet Events
//
//    This class reads events according to the GammaJetSel
//
//    first version: Hartmut Stadie 2008/12/12
//    $Id: PhotonJetReader.h,v 1.2 2008/12/16 15:21:26 stadie Exp $
//   
#ifndef PHOTONJETREADER_H
#define PHOTONJETREADER_H

#include "EventReader.h"

#include <string>

#include "GammaJetSel.h"

class PhotonJetReader : public EventReader{
 public:
  PhotonJetReader(const std::string& configfile, TParameters *p);
  virtual ~PhotonJetReader();
  int readEvents(std::vector<TData*>& data);
 private:
  TData* createTruthMultMessEvent();
  TData* createJetTruthEvent();

  GammaJetSel gammajet;
  double Et_cut_on_gamma,Et_cut_on_jet,Eta_cut_on_jet,Et_cut_on_genJet,Rel_cut_on_gamma; 
  int n_gammajet_events;
  int dataClass;
};


#endif
