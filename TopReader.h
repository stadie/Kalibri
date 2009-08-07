//
//    Reader for ttbar Events
//
//    This class reads events according to the TopSel
//
//    first version: Hartmut Stadie 2008/12/12
//    $Id: TopReader.h,v 1.6 2009/08/05 12:16:30 stadie Exp $
//   
#ifndef TOPREADER_H
#define TOPREADER_H

#include <string>

#include "TH2F.h"

#include "EventReader.h"
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

  bool createGenWHist_;
  TH2F* genWPtEta_;
};


#endif
