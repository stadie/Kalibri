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
//!  $Id: PhotonJetReader.h,v 1.3 2009/03/04 17:26:51 thomsen Exp $
// ----------------------------------------------------------------   
class PhotonJetReader : public EventReader{
 public:
  PhotonJetReader(const std::string& configfile, TParameters *p);
  virtual ~PhotonJetReader();
  int readEvents(std::vector<TData*>& data);

 private:
  TData* createTruthMultMessEvent();
  TData* createJetTruthEvent();

  GammaJetSel gammajet;
  double Et_cut_on_gamma;      //!< Minimum photon Et
  double Et_cut_on_jet;        //!< Minimum jet Et
  double Eta_cut_on_jet;       //!< Maximum absolute jet eta
  double Et_cut_on_genJet;     //!< Minimum genJet Et
  double Rel_cut_on_gamma;     //!< Minimum fraction of Et of non-leading jet to photon Et
  double Had_cut_min;          //!< Minimum jet Had/(Had+EMF)
  double Had_cut_max;          //!< Maximum jet Had/(Had+EMF)
  int    n_gammajet_events;    //!< Maximum number of read photon jet events
  int    dataClass;            //!< Data class, see also TData
};


#endif
