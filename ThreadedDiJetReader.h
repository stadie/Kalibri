#ifndef THREADEDDIJETREADER_H
#define THREADEDDIJETREADER_H
//!
//!  \brief Reader for dijet events
//!
//!  This class reads dijet events from a ROOT-tree. They are read
//!  by calling the method readEvents(std::vector<Event*>& data).
//!  The data is stored in a format derived from Event depending on the
//!  'Di-Jet data class' field in the config file:
//!    - class 0: Event_PtBalance
//!
//!  It is also possible to store the first two jets of the dijet
//!  event as a JetTruthEvent, where the \f$ p^{\textrm{gen}}_{T} \f$
//!  is used as truth. This is useful for Monte Carlo based calibration:
//!    - class 11: JetTruthEvent of Jet
//!    - class 12: JetTruthEvent of JetWithTowers
//!  In this case, the following cuts as defined in the config file
//!  are applied on each jet in this order:
//!    -# Number of jets greater than 2
//!    -# \f$ p^{\textrm{gen}}_{T} > \texttt{Et genJet min both Jets} \f$
//!    -# \f$ p^{\textrm{gen}}_{T} < \texttt{Et genJet max both Jets} \f$
//!    -# \f$ \Delta R < \texttt{DeltaR cut on jet matching} \f$
//!    -# \f$ p^{\textrm{jet}}_{T} < \texttt{Et cut on jet} \f$
//!    -# \f$ |\eta| > \texttt{Eta cut on jet} \f$
//!    -# \f$ \textrm{Hadronic fraction} < \texttt{Max had fraction} \f$
//!    -# \f$ \textrm{Hadronic fraction} > \texttt{Min had fraction} \f$
//!  
//!
//!  \author Hartmut Stadie
//!  \date 2008/12/12
//!  $Id: ThreadedDiJetReader.h,v 1.3 2010/11/01 15:47:40 stadie Exp $
// ----------------------------------------------------------------   
#include "DiJetReader.h"

#include <vector>

#include <boost/shared_ptr.hpp>
#include <boost/thread/thread.hpp>
class TChain;

class ThreadedDiJetReader : public DiJetReader{
 public:
  ThreadedDiJetReader(const std::string& configfile, Parameters *p, int niot);
  ~ThreadedDiJetReader();
  int readEvents(std::vector<Event*>& data);
  int readControlEvents(std::vector<Event*>& control, int id);

 private:
  class ReadThread {
  private:
    DiJetReader* reader_;
    std::vector<Event*> data_;
    struct read_events
    {
    private:
      ReadThread *parent_;
    public:
      read_events(ReadThread *parent) : parent_(parent) {}
      void operator()();
    };
    boost::shared_ptr<boost::thread> thread_;
    friend class read_events;  
  public:
    ReadThread(const std::string& configfile, Parameters* p);
    ~ReadThread();
    void start();
    bool isDone();
    int addEvents(std::vector<Event*>& data) const;
    void reset();
    int nEvents() const { return data_.size();}
    DiJetReader* reader() { return reader_;}
  };
  std::vector<ReadThread*> readers_;
  TTree* tree_;
};


#endif
