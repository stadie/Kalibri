//
//    first version: Hartmut Stadie 2008/12/12
//    $Id: ThreadedDiJetReader.cc,v 1.6 2011/01/19 15:07:14 stadie Exp $
//   
#include "ThreadedDiJetReader.h"

#include "ConfigFile.h"
#include "Parameters.h"
#include "NJetSel.h"
#include "TChain.h"
#include "TFile.h"
#include "TChainElement.h"
#include "JetBin.h"
#include "Binning.h"
#include "JetTruthEvent.h"

#include "CorFactorsFactory.h"
#include <boost/thread/thread.hpp>
#include <boost/thread/mutex.hpp>

//!  \brief Constructor
//!
//!  Reads data from ROOT trees and stores them in an \p NJetSel selector.
//!  The data can be stored in a format derived from \p Event (as specified
//!  in the 'Di-Jet data class' field in the config file) by calling the
//!  method readEvents(<tt>std::vector<Event*>& data</tt>). Additionally,
//!  the cut thresholds are read from the configfile.
//!
//!  \param configfile Name of configfile
//!  \param p Pointer to \p Parameters object
// ----------------------------------------------------------------   
ThreadedDiJetReader::ThreadedDiJetReader(const std::string& configfile, 
					 Parameters* p, int niot)
  : DiJetReader(configfile,p)
{  
  for(int i = 0 ; i < niot ; ++i) {
    readers_.push_back(new ReadThread(configfile, p));
  }
  tree_ = createTree("Di-Jet");
};

ThreadedDiJetReader::~ThreadedDiJetReader()
{
  for(std::vector<ReadThread*>::iterator i = readers_.begin() ;
      i != readers_.end() ; ++i) {
    delete *i;
  }
  readers_.clear();
};

//!  \brief Read dijet data and store in format as specified
//!         in the config file
//!  \param data Read data objects are appended to data
//!  \return Number of appended objects
// ----------------------------------------------------------------   
int ThreadedDiJetReader::readEvents(std::vector<Event*>& data)
{
  if(nDijetEvents_ == 0) return 0; 
  std::cout << "Reading events of type ";
  if(dataClass_ == 1) {
    std::cout << "'TwoJetsPtBalanceEvent'";
  } else if((dataClass_ == 11)  || (dataClass_ == 12) || (dataClass_ == 21)) {
    std::cout << "'JetTruthEvent'";
  } else if(dataClass_ == 5) {
    std::cout << "'SmearData'";
    if( !correctToL3_ && !correctL2L3_ ) {
      std::cerr << "WARNING: Jets are not corrected!\n";
      exit(9);
    }
  } else if (dataClass_ == 31) {
    std::cout << "'JetWidthEvent'";
  } else {
    std::cerr << "Unknown data class " << dataClass_ << '\n';
    exit(9);
  }
  std::cout << " (data class " << dataClass_ << "):\n";
  
  nDiJetCut_          = 0;
  nMinJetEt_          = 0;
  nMaxJetEt_          = 0;
  nTriggerSel_        = 0;
  nMinDijetEt_        = 0;
  nMaxDijetEt_        = 0;
  nTriggerSel_        = 0;
  nCutOn3rdJet_       = 0;
  nCutOnSoftJets_     = 0;
  nMinGenJetEt_       = 0;    
  nMaxGenJetEt_       = 0;     
  nMaxDeltaR_         = 0;
  nMinJetEta_         = 0;
  nMaxJetEta_         = 0;  
  nMinJetHadFraction_ = 0;     
  nMaxJetHadFraction_ = 0;     
  nMinDeltaPhi_       = 0;
  nReadEvts_          = 0;
  nGoodEvts_          = 0;

  TChain* chain = dynamic_cast<TChain*>(tree_);
  unsigned int id = 0;
  std::vector<TFile*> files;
  if(! chain) { 
    id = 1;
    readers_[0]->reader()->nJet_->Init(tree_);
    readers_[0]->start();
  } else {
    TObjArray *fileElements=chain->GetListOfFiles();
    TIter next(fileElements);
    TChainElement *chEl=0;
    while (( chEl=(TChainElement*)next() )) {
      TFile* f = TFile::Open(chEl->GetTitle());
      //std::cout << "opened file " << f->GetName() << '\n';
      files.push_back(f);
      TTree* tree = (TTree*)f->Get(chain->GetName()); 
      //tree->GetEntriesFast();
      //std::cout << "adding " << f->GetName() << " with " <<  tree->GetEntriesFast() << " entries\n";
      readers_[id]->reader()->nJet_->Init(tree);
      //readers_[id]->reader()->nJet_->GetEntry(0);
      ++id;
      if(id == readers_.size()) {
	for(std::vector<ReadThread*>::iterator i = readers_.begin() ;
	    i != readers_.end() ; ++i) {
	  (*i)->start();
	}
	for(std::vector<ReadThread*>::iterator i = readers_.begin() ;
	    i != readers_.end() ; ++i) {
	  //std::cout << "waiting for reader " << i - readers_.begin() << '\n';
	  if((*i)->isDone()) {
	    DiJetReader* djr = (*i)->reader();
	    //std::cout << "adding for reader " << i - readers_.begin() 
	    //	    << " with " <<  djr->nReadEvts_ << " events.\n";
	    nReadEvts_          += djr->nReadEvts_;
	    nDiJetCut_          += djr->nDiJetCut_;
	    nMinJetEt_          += djr->nMinJetEt_;
	    nMaxJetEt_          += djr->nMaxJetEt_;
	    nTriggerSel_        += djr->nTriggerSel_;
	    nMinDijetEt_        += djr->nMinDijetEt_;
	    nMaxDijetEt_        += djr->nMaxDijetEt_;
	    nCutOn3rdJet_       += djr->nCutOn3rdJet_;
	    nCutOnSoftJets_     += djr->nCutOnSoftJets_;
	    nMinGenJetEt_       += djr->nMinGenJetEt_;
	    nMaxGenJetEt_       += djr->nMaxGenJetEt_;  
	    nMaxDeltaR_         += djr->nMaxDeltaR_;
	    nMinJetEta_         += djr->nMinJetEta_;
	    nMaxJetEta_         += djr->nMaxJetEta_;
	    nMinJetHadFraction_ += djr->nMinJetHadFraction_;    
	    nMaxJetHadFraction_ += djr->nMaxJetHadFraction_;    
	    nMinDeltaPhi_       += djr->nMinDeltaPhi_; 
	    //std::cout << *i << ":events " << (*i)->nEvents() << ";" << djr->nReadEvts_<< '\n';
	  }
	}
	std::cout << nReadEvts_ << " events read\n";
	id = 0;
	if((nDijetEvents_ > 0) && (nReadEvts_ >= nDijetEvents_)) break; 
	for(std::vector<TFile*>::const_iterator f = files.begin() ;
	    f!= files.end() ; ++f){
	  (*f)->Close();
	  delete *f;
	}
	files.clear();
      }
    } 
    //read last batch
    for(unsigned int  i = 0 ; i < id ; ++i) {
      readers_[i]->start();
    }
  }
  for(unsigned int  i = 0 ; i < id ; ++i) {
    if(readers_[i]->isDone()) {
      DiJetReader* djr = readers_[i]->reader(); 
      nReadEvts_          += djr->nReadEvts_;
      nDiJetCut_          += djr->nDiJetCut_;
      nMinJetEt_          += djr->nMinJetEt_;
      nMaxJetEt_          += djr->nMaxJetEt_;
      nMinDijetEt_        += djr->nMinDijetEt_;
      nMaxDijetEt_        += djr->nMaxDijetEt_;
      nCutOn3rdJet_       += djr->nCutOn3rdJet_;
      nCutOnSoftJets_     += djr->nCutOnSoftJets_;
      nMinGenJetEt_       += djr->nMinGenJetEt_;
      nMaxGenJetEt_       += djr->nMaxGenJetEt_;  
      nMaxDeltaR_         += djr->nMaxDeltaR_;
      nMinJetEta_         += djr->nMinJetEta_;
      nMaxJetEta_         += djr->nMaxJetEta_;
      nMinJetHadFraction_ += djr->nMinJetHadFraction_;    
      nMaxJetHadFraction_ += djr->nMaxJetHadFraction_;    
      nMinDeltaPhi_       += djr->nMinDeltaPhi_;  
      //std::cout << readers_[i] << ":events " << readers_[i]->nEvents() << ";" << djr->nReadEvts_<< '\n';
    }
  }
  std::cout << nReadEvts_ << " events read\n";
    
  for(std::vector<TFile*>::const_iterator f = files.begin() ;
      f!= files.end() ; ++f){
    (*f)->Close();
    delete *f;
  }
  files.clear();
  for(std::vector<ReadThread*>::iterator i = readers_.begin() ;
      i != readers_.end() ; ++i) {
    nGoodEvts_ += (*i)->addEvents(data);
    (*i)->reset();
  }  
  if(dataClass_ == 21) {
    for(Binning::BinMap::const_iterator ijb = binning_->bins().begin() ; 
	ijb != binning_->bins().end() ; ++ijb) {
      if(ijb->second->nJets() > 10) {
	Jet *jet = ijb->second->jet();
	
	if(corFactorsFactory_) {
	  jet->updateCorFactors(corFactorsFactory_->create(jet));
	}
	if(correctToL3_) {
	  jet->correctToL3();
	} else if(correctL2L3_) {
	  jet->correctL2L3();
	} 
	data.push_back(new JetTruthEvent(jet,ijb->second->genPt(),ijb->second->nJets(),true));
	++nGoodEvts_; 
      }      
      delete ijb->second;
    }
    binning_->clear();
  }
  printCutFlow();
  std::cout << "Stored " << nGoodEvts_ << " dijet events.\n";
  delete tree_;
  return nGoodEvts_;
};


//!  \brief Read dijet data and store in format as specified
//!         in the config file
//!  \param data Read data objects are appended to data
//!  \return Number of appended objects
// ----------------------------------------------------------------   
int ThreadedDiJetReader::readControlEvents(std::vector<Event*>& control, int id)
{
  std::ostringstream name;
  name << "Di-Jet Control" << id;
  nDijetEvents_ = config_->read<int>("use "+name.str()+" events",0);
  //hack disable trigger
  requireTrigger_ = false;
  for(std::vector<ReadThread*>::iterator i = readers_.begin() ;
	    i != readers_.end() ; ++i) {
    (*i)->reader()->requireTrigger_ = false;
  }
  if(nDijetEvents_ == 0) return 0;
  tree_ = createTree(name.str());
  if(tree_->GetEntries() == 0) return 0;
  delete corFactorsFactory_;
  corFactorsFactory_ = 0;
  int nev = readEvents(control);
  std::cout << "Will use events for control plots.\n";
  return nev;
}

ThreadedDiJetReader::ReadThread::ReadThread(const std::string& configfile, 
					    Parameters* p) 
{
  reader_ = new DiJetReader(configfile,p);
  delete reader_->nJet_->fChain;
}
  
ThreadedDiJetReader::ReadThread::~ReadThread() 
{ 
  delete reader_;
}
    
void  ThreadedDiJetReader::ReadThread::start() {
  thread_ = new boost::thread(read_events(this)); 
}


bool ThreadedDiJetReader::ReadThread::isDone() { 
  thread_->join(); 
  delete thread_; 
  return true;
}

int ThreadedDiJetReader::ReadThread::addEvents(std::vector<Event*>& data) const
{
  data.insert(data.end(),data_.begin(),data_.end());
  return data_.size();
}

void ThreadedDiJetReader::ReadThread::reset() {
  data_.clear();
}

void ThreadedDiJetReader::ReadThread::read_events::operator()() {
  parent_->reader_->readEventsFromTree(parent_->data_);
}

