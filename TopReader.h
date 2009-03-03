//
//    Reader for ttbar Events
//
//    This class reads events according to the TopSel
//
//    first version: Hartmut Stadie 2008/12/12
//    $Id: TopReader.h,v 1.1 2008/12/12 17:06:00 stadie Exp $
//   
#ifndef TOPREADER_H
#define TOPREADER_H

#include "EventReader.h"

#include <string>

#include "TopSel.h"


class TData;

class TopReader : public EventReader{
 public:
  TopReader(const std::string& configfile, TParameters *p);
  virtual ~TopReader();
  int readEvents(std::vector<TData*>& data);
 private:
  TData* createTwoJetsInvMassEvents();

  TopSel top;
  double Et_cut_on_jet;
  bool useMassConstraintW, useMassConstraintTop;
  double massConstraint_W, massConstraint_Top; 
  int n_top_events;
  int dataClass;
};


#endif
