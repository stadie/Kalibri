//
//    Reader for Z Jet Events
//
//    This class reads events according to the ZJetSel
//
//    first version: Hartmut Stadie 2008/12/12
//    $Id: ZJetReader.h,v 1.1 2008/12/12 17:06:00 stadie Exp $
//   
#ifndef ZJETREADER_H
#define ZJETREADER_H

#include "EventReader.h"

#include <string>

#include "ZJetSel.h"

class ZJetReader : public EventReader{
 public:
  ZJetReader(const std::string& configfile, TParameters *p);
  virtual ~ZJetReader();
  int readEvents(std::vector<TData*>& data);
 private:
  TData* createTruthMultMessEvent();
  TData* createJetTruthEvent();

  ZJetSel zjet;
  double Et_cut_on_Z,Et_cut_on_jet;
  int n_zjet_events;
  int dataClass;
};


#endif
