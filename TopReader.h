//
//    Reader for ttbar Events
//
//    This class reads events according to the TopSel
//
//    first version: Hartmut Stadie 2008/12/12
//    $Id: TopReader.h,v 1.5 2009/08/04 15:37:41 snaumann Exp $
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

  TopSel top_;
  double minJetEt_;
  double maxJetEta_;
  double minJetHadFrac_;
  double maxJetHadFrac_;
  bool useToL3CorrectedJets_;
  bool useMassConstraintW_;
  bool useMassConstraintTop_;
  bool useGenJetInvMass_;
  double massConstraintW_;
  double massConstraintTop_; 
  int nTopEvents_;
  int dataClass_;
};


#endif
