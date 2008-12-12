//
//    Reader for Di-Jet Events
//
//    This class reads events according to the NJetSel
//
//    first version: Hartmut Stadie 2008/12/12
//    $Id: caliber.h,v 1.33 2008/11/20 16:38:03 stadie Exp $
//   
#ifndef DIJETREADER_H
#define DIJETREADER_H

#include "EventReader.h"

#include <string>

#include "NJetSel.h"

class DiJetReader : public EventReader{
 public:
  DiJetReader(const std::string& configfile, TParameters *p);
  virtual ~DiJetReader();
  int readEvents(std::vector<TData*>& data);
 private:
  NJetSel njet;
  double Et_cut_nplus1Jet,Rel_cut_on_nJet;
  int n_dijet_events;
};


#endif
