#ifndef ZJETREADER_H
#define ZJETREADER_H
//!
//!    \brief Reader for Z Jet Events
//!
//!    This class reads events according to the ZJetSel
//!
//!    \author Hartmut Stadie
//!    \date 2008/12/12
//!    $Id: ZJetReader.h,v 1.4 2009/04/17 14:28:08 mschrode Exp $
// ----------------------------------------------------------------   

#include "EventReader.h"

#include <string>
#include <memory>

class ZJetSel;

class ZJetReader : public EventReader{
 public:
  ZJetReader(const std::string& configfile, TParameters *p);
  virtual ~ZJetReader();
  int readEvents(std::vector<Event*>& data);
 private:
  Event* createJetTruthEvent();

  std::auto_ptr<ZJetSel> zjet;
  double Et_cut_on_Z;          //!< Minimum Z Et
  double Et_cut_on_jet;        //!< Minimum jet Et
  double Eta_cut_on_jet;       //!< Maximum absolute jet eta
  double Et_cut_on_genJet;     //!< Minimum genJet Et
  double Had_cut_min;          //!< Minimum jet Had/(Had+EMF)
  double Had_cut_max;          //!< Maximum jet Had/(Had+EMF)
  int    n_zjet_events;        //!< Maximum number of read photon jet events
  int    dataClass;            //!< Data class, see also Event
};


#endif
