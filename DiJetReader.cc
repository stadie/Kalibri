//
//    first version: Hartmut Stadie 2008/12/12
//    $Id: DiJetReader.cc,v 1.62 2010/10/20 11:28:07 stadie Exp $
//   
#include "DiJetReader.h"

#include "CalibData.h"
#include "SmearDiJet.h"
#include "ConfigFile.h"
#include "Parameters.h"
#include "JetTruthEvent.h"
#include "JetWithTowers.h"
#include "TwoJetsPtBalanceEvent.h"
#include "JetConstraintEvent.h"
#include "NJetSel.h"
#include "CorFactors.h"
#include "CorFactorsFactory.h"
#include "Function.h"
#include "SmearFunction.h"
#include "JetBin.h"
#include "Binning.h"

#include "TVector2.h"
#include "TRandom3.h"

#include <algorithm>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <boost/thread/thread.hpp>
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
DiJetReader::DiJetReader(const std::string& configfile, Parameters* p)
  : EventReader(configfile,p), nJet_(new NJetSel()), 
    rand_(new TRandom3(0)), zero_(0)
{
  // Maximum number of read events
  nDijetEvents_ = config_->read<int>("use Di-Jet events",-1);
  prescale_ = config_->read<int>("Di-Jet prescale",1);
  weightRelToNtuple_ = config_->read<double>("Di-Jet weight relative to ntuple weight",1.);

  // Cuts
  ptRef_             = config_->read<double>("Reference pt",1.);
  minJetEt_          = config_->read<double>("Et min cut on jet",0.0);
  maxJetEt_          = config_->read<double>("Et max cut on jet",10000.0);
  minDijetEt_        = config_->read<double>("Et min cut on dijet",0.0); 
  maxDijetEt_        = config_->read<double>("Et max cut on dijet",10000.0); 
  max3rdJetEt_       = config_->read<double>("Et cut on n+1 Jet",10000.0);
  minRel3rdJetEt_    = config_->read<double>("Min cut on relative n+1 Jet Et",0.);
  maxRel3rdJetEt_    = config_->read<double>("Max cut on relative n+1 Jet Et",1.);
  maxRelSoftJetEt_   = config_->read<double>("Max cut on relative Soft Jet Et",1.);
  minJetEta_         = config_->read<double>("Eta min cut on jet",0.0);
  maxJetEta_         = config_->read<double>("Eta max cut on jet",5.0);
  minJetHadFraction_ = config_->read<double>("Min had fraction",0.07);
  maxJetHadFraction_ = config_->read<double>("Max had fraction",0.95);
  minDeltaPhi_       = config_->read<double>("Min Delta Phi",2.5);
  minGenJetEt_       = config_->read<double>("Et genJet min",0.0);
  maxGenJetEt_       = config_->read<double>("Et genJet max",10000.0);
  maxDeltaR_         = config_->read<double>("DeltaR cut on jet matching",0.25);

  // Loose JetID for barrel
  // Replace by ntuple variabel in newer Kalibri-trees
  minJetN90Hits_ = 2;
  maxJetFHPD_ = 0.98;
  maxJetFRBX_ = 0.98;

  // Counter for cutflow
  nDiJetCut_          = 0;
  nMinJetEt_          = 0;
  nMaxJetEt_          = 0;
  nMinDijetEt_        = 0;
  nMaxDijetEt_        = 0;
  nCutOn3rdJet_       = 0;
  nCutOnSoftJets_     = 0;
  nMinJetEta_         = 0;  
  nMaxJetEta_         = 0;  
  nMinDeltaPhi_       = 0;
  nMinJetHadFraction_ = 0;     
  nMaxJetHadFraction_ = 0;     
  nMinGenJetEt_       = 0;    
  nMaxGenJetEt_       = 0;     
  nMaxDeltaR_         = 0;
  // Integration parameter for SmearData
  maxNIter_  = config_->read<int>("DiJet integration number of iterations",5);
  eps_       = config_->read<double>("DiJet integration epsilon",1.E-5);
  min_       = config_->read<double>("DiJet integration min",0.);
  max_       = config_->read<double>("DiJet integration max",1.);
  // Data class
  dataClass_ = config_->read<int>("Di-Jet data class", 0);
  if((dataClass_ != 1)&&(dataClass_ != 11)&&(dataClass_ != 12)&&(dataClass_ != 21)&&(dataClass_ != 5)) {
    std::cerr << "DiJetReader: Unknown data class " << dataClass_ << ". Aborting!" << std::endl;
    exit(9);
  }
  const std::string name = "jet binning";
  std::vector<std::string> vars = bag_of_string(config_->read<std::string>(name+" variables","")); 
  int j = 0;
  for(std::vector<std::string>::const_iterator i = vars.begin() ; 
      i != vars.end() ; ++i) {
    if(*i == "pt") {
      vars_[j] = &genjetpt_;
    } else if(*i == "eta") {
      vars_[j] = &jeteta_;
    } else if(*i == "sigmaphi") {
      vars_[j] = &sigmaphi_;    
    } else if(*i == "sigmaeta") {
      vars_[j] = &sigmaeta_;   
    } else if(*i == "sumsigmaetaphi") {
      vars_[j] = &sumsigmaetaphi_;
    } else if(*i == "emf") {
      vars_[j] = &emf_; } 
    else if(*i == "meanMoment") {
      vars_[j] = &meanMoment_;   
    } else {
      std::cerr << "unknown binning varible: " << *i << '\n';
      exit(3);
    }
    ++j;
  }
  for(; j < 4 ; ++j) vars_[j] = &zero_;
  
  jetIndices_.resize(8,0);

  //std::cout << "size:" << sizeof(JetTruthEvent) << ", " << sizeof(Jet) << ", " << sizeof(Measurement) << '\n';
}

DiJetReader::~DiJetReader() {
  delete rand_;
}



//!  \brief Read dijet data and store in format as specified
//!         in the config file
//!  \param data Read data objects are appended to data
//!  \return Number of appended objects
// ----------------------------------------------------------------   
int DiJetReader::readEvents(std::vector<Event*>& data)
{
  // Input files
  nJet_->Init(createTree("Di-Jet")); 
  // Some informative output for the interested calibrator
  // Check of correct data class
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
  } else {
    std::cerr << "Unknown data class " << dataClass_ << '\n';
    exit(9);
  }
  std::cout << " (data class " << dataClass_ << "):\n";

  int nev = readEventsFromTree(data);  
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
      }      
      delete ijb->second;
    }
    binning_->clear();
  }
  printCutFlow();
  std::cout << "Stored " << nGoodEvts_ << " dijet events for analysis.\n";
  delete nJet_->fChain;
  return nev;
}

//!  \brief Read dijet data and store in format as specified
//!         in the config file
//!  \param data Read data objects are appended to data
//!  \return Number of appended objects
// ----------------------------------------------------------------   
int DiJetReader::readEventsFromTree(std::vector<Event*>& data)
{
  if(nDijetEvents_ == 0) return 0;
  
  // Reset counters of rejected events
  nReadEvts_          = 0;
  nGoodEvts_          = 0;
  nDiJetCut_          = 0;
  nMinJetEt_          = 0;
  nMaxJetEt_          = 0;
  nMinDijetEt_        = 0;
  nMaxDijetEt_        = 0;
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

  //Run jet-Jet stuff  
  int nevent    = nJet_->fChain->GetEntries();  // Number of events in chain


  Int_t cachesize = 10000000; //10 MBytes
  nJet_->fChain->SetCacheSize(cachesize); 
  if((dataClass_ == 11)||(dataClass_ == 21)||(dataClass_ == 1)) { 
    nJet_->fChain->SetBranchStatus("Track*",0);
    nJet_->fChain->SetBranchStatus("Tow*",0);
    nJet_->fChain->SetBranchStatus("Vtx*",0);
    nJet_->fChain->SetBranchStatus("GenPart*",0);
    nJet_->fChain->SetBranchStatus("GenPartId*",1);
  } else if(dataClass_ == 12) {
    nJet_->fChain->SetBranchStatus("Track*",0);  
    nJet_->fChain->SetBranchStatus("Vtx*",0);
    nJet_->fChain->SetBranchStatus("GenPart*",0);
    nJet_->fChain->SetBranchStatus("GenPartId*",1);
  } else if(dataClass_ == 5) {
    nJet_->fChain->SetBranchStatus("Track*",0); 
    nJet_->fChain->SetBranchStatus("Vtx*",0);
    nJet_->fChain->SetBranchStatus("GenPart*",0);
  }
  // Read the events
  for (int i=0 ; i < nevent ; i+= prescale_) {
    if((i+1)%10000==0) std::cout << "  " << i+1 << std::endl;
    nJet_->fChain->GetEvent(i); 
    if (nJet_->NobjTow>10000 || nJet_->NobjJet>100) {
      std::cerr << "ERROR: Increase array sizes in NJet_Selector; NobjTow="
		<< nJet_->NobjTow<<", NobjJet="<<nJet_->NobjJet<<"!"<<std::endl;
      exit(9);
    }

//     if( nJet_->GenJetColPt[1] > nJet_->GenJetColPt[0] 
// 	|| nJet_->GenJetColPt[2] > nJet_->GenJetColPt[1] ) {
//     std::cout << "\nEvent " << i << std::endl;
//     std::cout << nJet_->JetPt[0] << std::endl;
//     std::cout << nJet_->JetPt[1] << std::endl;
//     std::cout << nJet_->JetPt[2] << std::endl;
//     std::cout << nJet_->GenJetColPt[0] << std::endl;
//     std::cout << nJet_->GenJetColPt[1] << std::endl;
//     std::cout << nJet_->GenJetColPt[2] << std::endl;
//     }



    if(dataClass_ == 1) {
      nReadEvts_++;
      TwoJetsPtBalanceEvent* td = createTwoJetsPtBalanceEvent(); 
      if(td) {
	//second jet should be central...
	if(std::abs(td->getJet1()->eta()) < 1.3) {
	  data.push_back(new TwoJetsPtBalanceEvent(td->getJet2()->clone(),td->getJet1()->clone(),
						   td->getJet3() ? td->getJet3()->clone():0,
						   td->ptHat(),td->weight()));
	  ++nGoodEvts_;
	}
	if(std::abs(td->getJet2()->eta()) < 1.3) {
	  data.push_back(td ); 
	  ++nGoodEvts_;
	} else {
	  delete td;
	}
      } 
    } else if((dataClass_ == 11)  || (dataClass_ == 12) || (dataClass_ == 21)) {
      nReadEvts_++;
      int nAddedJets = createJetTruthEvents(data);
      if( nAddedJets ) nGoodEvts_ += nAddedJets;    
    } else if(dataClass_ == 5) {
      for(int calls = 0; calls < 1; calls++) {
	nReadEvts_++;
	Event* td = createSmearEvent(calls); 
	if(td) {
	  nGoodEvts_++;
	  data.push_back(td ); 
	} 
      }
    }
    if(nReadEvts_>=nDijetEvents_ && nDijetEvents_>=0 ) break;
  }
  //printCutFlow(); 
  return nGoodEvts_;
}

void DiJetReader::printCutFlow()
{
  // Print cut flow
  std::cout << "Read " << nReadEvts_ << " dijet events:\n";
  if( dataClass_ == 1 ) {
    std::cout << "  " << (nReadEvts_-=nDiJetCut_) << std::flush;
    std::cout << " events with 2 or more jets\n"; 
    std::cout << "  " << (nReadEvts_-=nMinDeltaPhi_) << std::flush;
    std::cout << " dijet events with DeltaPhi > " << minDeltaPhi_ << '\n';;
    std::cout << "  " << (nReadEvts_-=nMinGenJetEt_) << std::flush;
    std::cout << " dijet events with ptgen > " << minGenJetEt_ << " GeV\n";
    std::cout << "  " << (nReadEvts_-=nMaxGenJetEt_) << std::flush;
    std::cout << " dijet events with ptgen < " << maxGenJetEt_ << " GeV\n";
    std::cout << "  " << (nReadEvts_-=nMaxDeltaR_) << std::flush;
    std::cout << " dijet events with DeltaR < " << maxDeltaR_ << "\n";
    std::cout << "  " << (nReadEvts_-=nMinJetEt_) << std::flush;
    std::cout << " dijet events Et > " << minJetEt_ << " GeV\n";
    std::cout << "  " << (nReadEvts_-=nMaxJetEt_) << std::flush;
    std::cout << " dijet events Et < " << maxJetEt_ << " GeV\n";
    std::cout << "  " << (nReadEvts_-=(nMinJetEta_+nMaxJetEta_)) << std::flush;
    std::cout << " dijet events with " << minJetEta_ << " < |eta| < " << maxJetEta_ << "\n";
    std::cout << "  " << (nReadEvts_-=nMinJetHadFraction_) << std::flush;
    std::cout << " dijet events with hadronic fraction > " << minJetHadFraction_ << "\n";
    std::cout << "  " << (nReadEvts_-=nMaxJetHadFraction_) << std::flush;
    std::cout << " dijet events with jet id and hadronic fraction < " << maxJetHadFraction_ << "\n";   
    std::cout << "  " << (nReadEvts_-=nCutOn3rdJet_) << std::flush;
    std::cout << " dijet events with " << minRel3rdJetEt_ << " < pt(jet3) / ptAve < " << maxRel3rdJetEt_ << "\n";
    std:: cout << "  " << (nReadEvts_-=nCutOnSoftJets_) << std::flush;
    std::cout << " dijet events with pt(jet>3) / ptAve < " << maxRelSoftJetEt_ << "\n";
  } else if( dataClass_ == 11 || dataClass_ == 12 || dataClass_ == 21 ) {
    std::cout << "  " << (nReadEvts_-=nDiJetCut_) << std::flush;
    std::cout << " events with 2 or more jets\n";
    std::cout << "  That are " << (nReadEvts_*=2) << " jet-truth events:\n";
    std::cout << "    " << (nReadEvts_-=nMinGenJetEt_) << std::flush;
    std::cout << " jet-truth events with ptgen > " << minGenJetEt_ << "\n";
    std::cout << "    " << (nReadEvts_-=nMaxGenJetEt_) << std::flush;
    std::cout << " jet-truth events with ptgen < " << maxGenJetEt_ << "\n";
    std::cout << "    " << (nReadEvts_-=nMaxDeltaR_) << std::flush;
    std::cout << " jet-truth events with DeltaR < " << maxDeltaR_ << "\n";
    std::cout << "    " << (nReadEvts_-=nMinJetEt_) << std::flush;
    std::cout << " jet-truth events Et > " << minJetEt_ << "\n";
    std::cout << "    " << (nReadEvts_-=nMaxJetEt_) << std::flush;
    std::cout << " jet-truth events Et < " << maxJetEt_ << " GeV\n";
    std::cout << "  " << (nReadEvts_-=(nMinJetEta_+nMaxJetEta_)) << std::flush;
    std::cout << " dijet events with " << minJetEta_ << " < |eta| < " << maxJetEta_ << "\n";
    std::cout << "    " << (nReadEvts_-=nMinJetHadFraction_) << std::flush;
    std::cout << " jet-truth events with hadronic fraction > " << minJetHadFraction_ << "\n";
    std::cout << "    " << (nReadEvts_-=nMaxJetHadFraction_) << std::flush;
    std::cout << " jet-truth events with hadronic fraction < " << maxJetHadFraction_ << "\n";
  } else if( dataClass_ == 5 ) {
    std::cout << "  " << (nReadEvts_-=nDiJetCut_) << std::flush;
    std::cout << " events with 2 or more jets\n";
    std::cout << "  " << (nReadEvts_-=nMinGenJetEt_) << std::flush;
    std::cout << " dijet events with ptgen > " << minGenJetEt_ << " GeV\n";
    std::cout << "  " << (nReadEvts_-=nMaxGenJetEt_) << std::flush;
    std::cout << " dijet events with ptgen < " << maxGenJetEt_ << " GeV\n";
    std::cout << "  " << (nReadEvts_-=nMinJetEt_) << std::flush;
    std::cout << " dijet events pt(1) > " << minJetEt_ << " GeV\n";
    std::cout << "  " << (nReadEvts_-=nMaxJetEt_) << std::flush;
    std::cout << " dijet events pt(1) < " << maxJetEt_ << " GeV\n";
    std::cout << "  " << (nReadEvts_-=nMinDijetEt_) << std::flush;
    std::cout << " dijet events pt(ave) > " << minDijetEt_ << " GeV\n";
    std::cout << "  " << (nReadEvts_-=nMaxDijetEt_) << std::flush;
    std::cout << " dijet events pt(ave) < " << maxDijetEt_ << " GeV\n";
    std::cout << "  " << (nReadEvts_-=(nMinJetEta_+nMaxJetEta_)) << std::flush;
    std::cout << " dijet events with " << minJetEta_ << " < |eta| < " << maxJetEta_ << "\n";
    std::cout << "  " << (nReadEvts_-=nMinDeltaPhi_) << std::flush;
    std::cout << " dijet events with DeltaPhi > " << minDeltaPhi_ << '\n';;
    std::cout << "  " << (nReadEvts_-=nCutOn3rdJet_) << std::flush;
    std::cout << " dijet events with " << minRel3rdJetEt_ << " < pt(jet3) / (" << ptRef_ << "GeV) < " << maxRel3rdJetEt_ << "\n";
    std::cout << "  " << (nReadEvts_-=nCutOnSoftJets_) << std::flush;
    std::cout << " dijet events with pt(jet>3) / (" << ptRef_ << " GeV) < " << maxRelSoftJetEt_ << "\n";
    std::cout << "  " << (nReadEvts_-=nMinJetHadFraction_) << std::flush;
    std::cout << " dijet events with hadronic fraction > " << minJetHadFraction_ << "\n";
    std::cout << "  " << (nReadEvts_-=nMaxJetHadFraction_) << std::flush;
    std::cout << " dijet events with hadronic fraction < " << maxJetHadFraction_ << "\n";
    std::cout << "  " << (nReadEvts_-=nMaxDeltaR_) << std::flush;
    std::cout << " dijet events with DeltaR < " << maxDeltaR_ << "\n";
  }
}
  
//!  \brief Read dijet data and store in format as specified
//!         in the config file
//!  \param data Read data objects are appended to data
//!  \return Number of appended objects
// ----------------------------------------------------------------   
int DiJetReader::readControlEvents(std::vector<Event*>& control, int id)
{ 
  std::ostringstream name;
  name << "Di-Jet Control" << id;
  nDijetEvents_ = config_->read<int>("use "+name.str()+" events",0);
  if(nDijetEvents_ == 0) return 0;
  TTree* tree = createTree(name.str());
  if(tree->GetEntries() == 0) return 0;  
  delete corFactorsFactory_;
  corFactorsFactory_ = 0;
  nJet_->Init(tree);
  int nev = readEventsFromTree(control);
  printCutFlow();
  std::cout << "Stored " << nev << " dijet events for control plots.\n";
  delete nJet_->fChain;
  return nev;
}


//!  \brief Use first two jets of dijet event as two JetTruthEvent
//!         objects where genjet Et is truth
//!  \note The dijets are ordered in genjet Et
//!  \param data Read JetTruthEvent objects are appended to data
//!  \return Number of appended JetTruthEvent objects (0 - 2)
// ----------------------------------------------------------------   
int DiJetReader::createJetTruthEvents(std::vector<Event*>& data)
{
  int injet = 2;
  if( nJet_->NobjGenJet < injet ) {
    nDiJetCut_++;
    return 0;
  }

  int     njets = 0;   // Number of stored JetTruthEvents; 0 or 2
  double * terr = new double[nJet_->NobjTow];


  // Loop over two jets with highest genjet Et
  for(int genJetIdx = 0; genJetIdx < 2; genJetIdx++) {
    int calJetIdx = nJet_->GenJetColJetIdx[genJetIdx]; // Closest (DeltaR) calo jet
    if( nJet_->NobjJet <= calJetIdx ) {
      nDiJetCut_++;
      return 0;
    }

    // Cuts
    if( nJet_->GenJetColPt[genJetIdx] < minGenJetEt_ ) {
      nMinGenJetEt_++;
      continue;
    } else if( nJet_->GenJetColPt[genJetIdx] > maxGenJetEt_ ) {
      nMaxGenJetEt_++;
      continue;
    }
    double dphi        = TVector2::Phi_mpi_pi( nJet_->JetPhi[calJetIdx] - nJet_->GenJetColPhi[genJetIdx] );
    double deta        = nJet_->JetEta[calJetIdx] - nJet_->GenJetColEta[genJetIdx];
    double drJetGenjet = sqrt( deta*deta + dphi*dphi );
    if( drJetGenjet > maxDeltaR_ ) {
      nMaxDeltaR_++;
      continue;
    } else if( nJet_->JetPt[calJetIdx] < minJetEt_ ) {
      nMinJetEt_++;
      continue;
    } else if( std::abs(nJet_->JetEta[calJetIdx]) > maxJetEta_ ) {
      nMaxJetEta_++;
      continue;
    }
//     if(nJet_->JetEtWeightedSigmaPhi[calJetIdx] < 0) {
//       std::cout << "warning: weighted sigma_phi < 0 : " << nJet_->JetEtWeightedSigmaPhi[calJetIdx] << std::endl;
//       continue;
//     }   
//     if(nJet_->JetEtWeightedSigmaPhi[calJetIdx] > 1.0) {
//       std::cout << "warning: weighted sigma_phi < 1.0 : " << nJet_->JetEtWeightedSigmaPhi[calJetIdx] << " jet Et: " << nJet_->JetEt[calJetIdx] << std::endl;
//       continue;
//     }
//     if(nJet_->JetEtWeightedSigmaPhi[calJetIdx] != nJet_->JetEtWeightedSigmaPhi[calJetIdx]) {
//       std::cout << "warning: weighted sigma_phi is nan: " << nJet_->JetEtWeightedSigmaPhi[calJetIdx] << " jet Et: " << nJet_->JetEt[calJetIdx] << std::endl;
//       continue;
//     }
    // Construct event   
    Measurement tower;
    double err2 = 0;
    if(dataClass_ == 12) {
      for(int n=0; n<nJet_->NobjTow; ++n){
	if(nJet_->Tow_jetidx[n] != calJetIdx) continue;//look for ij-jet's towers
	tower.pt     = nJet_->TowEt[n];
	double scale = nJet_->TowEt[n]/nJet_->TowE[n];
	tower.EMF    = nJet_->TowEm[n]*scale;
	tower.HadF   = nJet_->TowHad[n]*scale;
	tower.OutF   = nJet_->TowOE[n]*scale;
	tower.eta    = nJet_->TowEta[n];
	tower.phi    = nJet_->TowPhi[n];
	tower.E      = nJet_->TowE[n];
	terr[n]      = tower_error_param(&tower.pt,&tower,0); 
	if(terr[n] == 0) {
	  //assume toy MC???
	  terr[n] = Parameters::toy_tower_error_parametrization(&tower.pt,&tower);
	}
	terr[n]  *= terr[n];
	err2     += terr[n];
      }
    }
    // Cuts on hadronic fraction 
    if(1 - nJet_->JetEMF[calJetIdx] < minJetHadFraction_) {
      nMinJetHadFraction_++;
      continue;
    } else if(1 - nJet_->JetEMF[calJetIdx] > maxJetHadFraction_) { 
      nMaxJetHadFraction_++;
      continue;
    }
    tower.pt      = nJet_->JetPt[calJetIdx];
    tower.EMF     = tower.pt * nJet_->JetEMF[calJetIdx];
    tower.HadF    = tower.pt * (1.0 - nJet_->JetEMF[calJetIdx]);
    tower.OutF    = 0;
    tower.eta     = nJet_->JetEta[calJetIdx];
    tower.phi     = nJet_->JetPhi[calJetIdx];
    tower.E       = nJet_->JetE[calJetIdx];
    double err    = jet_error_param(&tower.pt,&tower,0);
    err2         += err * err;

    if(dataClass_ == 21) {
      
      jeteta_ = nJet_->JetEta[calJetIdx];
      genjetpt_ = nJet_->GenJetColPt[genJetIdx];
      sigmaeta_ = nJet_->JetEtWeightedSigmaEta[calJetIdx];
      sigmaphi_ = nJet_->JetEtWeightedSigmaPhi[calJetIdx];
      sumsigmaetaphi_ = sigmaeta_ + sigmaphi_;
      meanMoment_ = 0.5 * sumsigmaetaphi_;
      emf_ = nJet_->JetEMF[calJetIdx];
      boost::mutex::scoped_lock lock(dijetmutex);
      JetBin* bin = (*binning_)(*vars_[0],*vars_[1],*vars_[2],*vars_[3]);
      // std::cout << "adding jet to bin:" << ibin << " with pt=" << nJet_->GenJetColPt[genJetIdx] << ", eta = " << nJet_->JetEta[calJetIdx] << std::endl;
      if(! bin) { 
	bin = binning_->setBin(new JetBin(par_->jet_function(nJet_->JetIEta[calJetIdx],nJet_->JetIPhi[calJetIdx]),jet_error_param,par_->global_jet_function()),*vars_[0],*vars_[1],*vars_[2],*vars_[3]);
      }
      CorFactors* cf = createCorFactors(calJetIdx); 
      bin->addJet(nJet_->JetPt[calJetIdx],tower.EMF,tower.HadF,
		  tower.OutF,nJet_->JetE[calJetIdx],
		  nJet_->JetEta[calJetIdx],nJet_->JetPhi[calJetIdx],
		  nJet_->JetEtWeightedSigmaPhi[calJetIdx],
		  nJet_->JetEtWeightedSigmaEta[calJetIdx],
		  nJet_->GenJetColPt[genJetIdx],drJetGenjet,
		  *cf);
      delete cf;
      ++njets;
      continue;
    }

    Jet *jet;
    if(dataClass_ == 12) {
      JetWithTowers *jt = 
	new JetWithTowers(nJet_->JetPt[calJetIdx],tower.EMF,tower.HadF,
			  tower.OutF,nJet_->JetE[calJetIdx],
			  nJet_->JetEta[calJetIdx],nJet_->JetPhi[calJetIdx],
			  nJet_->JetEtWeightedSigmaPhi[calJetIdx],
			  nJet_->JetEtWeightedSigmaEta[calJetIdx], Jet::uds,
			  nJet_->GenJetColPt[genJetIdx],drJetGenjet,
			  createCorFactors(calJetIdx),
			  par_->jet_function(nJet_->JetIEta[calJetIdx],
					     nJet_->JetIPhi[calJetIdx]),
			  jet_error_param,par_->global_jet_function());
      for(int j = 0 ; j < nJet_->NobjTow ; ++j) {
	if (nJet_->Tow_jetidx[j]!= calJetIdx) continue;//look for ij-jet's towers
	double scale = nJet_->TowEt[j]/nJet_->TowE[j];
	jt->addTower(nJet_->TowEt[j],nJet_->TowEm[j]*scale,
		     nJet_->TowHad[j]*scale,nJet_->TowOE[j]*scale,
		     nJet_->TowE[j],nJet_->TowEta[j],nJet_->TowPhi[j],
		     par_->tower_function(nJet_->TowId_eta[j],nJet_->TowId_phi[j]),
		     tower_error_param);
      }
      jet = jt;
    }
    else { 
      jet = new Jet(nJet_->JetPt[calJetIdx],tower.EMF,tower.HadF,
		    tower.OutF,nJet_->JetE[calJetIdx],
		    nJet_->JetEta[calJetIdx],nJet_->JetPhi[calJetIdx],
		    nJet_->JetEtWeightedSigmaPhi[calJetIdx],
		    nJet_->JetEtWeightedSigmaEta[calJetIdx],
		    Jet::flavorFromPDG(nJet_->GenPartId_algo[calJetIdx]),
		    nJet_->GenJetColPt[genJetIdx],drJetGenjet,
		    createCorFactors(calJetIdx),
		    par_->jet_function(nJet_->JetIEta[calJetIdx],nJet_->JetIPhi[calJetIdx]),
		    jet_error_param,par_->global_jet_function());    
    }
    if(corFactorsFactory_) {
      jet->updateCorFactors(corFactorsFactory_->create(jet));
    }
    if(correctToL3_) {
      jet->correctToL3();
    } else if(correctL2L3_) {
      jet->correctL2L3();
    }
    JetTruthEvent* jte = new JetTruthEvent(jet,nJet_->GenJetColPt[genJetIdx],nJet_->Weight);
    data.push_back(jte);
    ++njets;
    //add jet to constraints
    for(unsigned int i = 0 ; i < constraints_.size() ; ++i) {
      constraints_[i]->addJet(nJet_->GenJetColPt[genJetIdx],jet,new Function(&Parametrization::correctedJetEt,0,0,0,0,cp_));
    }
  }
  delete [] terr;
  return njets;
}



//!  \brief Create \p SmearDiJet event for jet smearing
//!
//!  Uses the three jets with the highest corrected
//!  calo pt. For creation of a \p SmearDiJet event,
//!  the JetMET L2L3 correction is applied.
//!
//!  \return A \p SmearDiJet if all cuts are passed,
//!          else 0
// ----------------------------------------------------------------   
Event* DiJetReader::createSmearEvent(int callIdx) {
  //return createSmearEventGenOrdered();
  return createSmearEventCaloOrdered();
}


// ----------------------------------------------------------------   
Event* DiJetReader::createSmearEventCaloOrdered()
{
  SmearDiJet *dijet = 0;
  if( jetIndices_.size() != 8 ) std::cerr << "%%%%%%%%%%%% WARNING %%%%%%%%%%%%%%%%%\n";

  // There should be at least two calo jets in the event
  if( nJet_->NobjJet < 2 ) {
    nDiJetCut_++;
  } else {
    // Find indices of 2(3) leading L2L3 corrected calo jets
    size_t nJets = nJet_->NobjJet;
    if( jetIndices_.size() < nJets ) nJets = jetIndices_.size();
    for(size_t i = 0; i < nJets; ++i) {
      // jetIndices_[i] = new Jet::JetIndex(i,nJet_->JetPt[i]*nJet_->JetCorrL2L3[i]);
      if( i < 2 ) jetIndices_[i] = new Jet::JetIndex(i,nJet_->JetPt[i]*nJet_->JetCorrL2L3[i]);
      else jetIndices_[i] = new Jet::JetIndex(i,nJet_->JetPt[i]);
    }
    //std::sort(jetIndices_.begin(),jetIndices_.begin()+nJets,Jet::JetIndex::ptGreaterThan);

    double deltaPhi12 = TVector2::Phi_mpi_pi(nJet_->JetPhi[jetIndices_[0]->idx_]-nJet_->JetPhi[jetIndices_[1]->idx_]);

    // Phi of dijet axis
    double pPhi = TVector2::Phi_mpi_pi(nJet_->JetPhi[jetIndices_[0]->idx_]-0.5*deltaPhi12+M_PI/2.);
    double pJ3 = 0.;
    double ptJet3 = 0.;
    if( nJet_->NobjJet > 2 ) {
      pJ3 = jetIndices_[2]->pt_*cos(TVector2::Phi_mpi_pi(pPhi-nJet_->JetPhi[jetIndices_[2]->idx_]));
      ptJet3 = jetIndices_[2]->pt_;
    }
    double pSJ = 0.;
    for(int i = 3; i < nJet_->NobjJet; ++i) {
      int idx = i;
      if( i < static_cast<int>(nJets) ) idx = jetIndices_[i]->idx_;
      pSJ += nJet_->JetPt[idx]*cos(TVector2::Phi_mpi_pi(pPhi-nJet_->JetPhi[idx]));
    }
    double ptJet4 = nJets < 4 ? 0. : jetIndices_[3]->pt_;
    double pUCE = 0.;

    // Check if event is selected
    bool isGoodEvt = true;
    if( nJet_->GenJetPt[jetIndices_[0]->idx_] < minGenJetEt_ || 
	nJet_->GenJetPt[jetIndices_[1]->idx_] < minGenJetEt_ ) {
      nMinGenJetEt_++;
      isGoodEvt = false;
    }
    else if( nJet_->GenJetPt[jetIndices_[0]->idx_] > maxGenJetEt_ || 
	     nJet_->GenJetPt[jetIndices_[1]->idx_] > maxGenJetEt_ ) {
      nMaxGenJetEt_++;
      isGoodEvt = false;
    }
    else if( jetIndices_[0]->pt_+jetIndices_[1]->pt_ < 2.*minDijetEt_ ) {
      nMinDijetEt_++;
      isGoodEvt = false;
    }
    else if( jetIndices_[0]->pt_+jetIndices_[1]->pt_ > 2.*maxDijetEt_ ) {
      nMaxDijetEt_++;
      isGoodEvt = false;
    }
    else if( std::abs(nJet_->JetEta[jetIndices_[0]->idx_]) < minJetEta_  ||
	     std::abs(nJet_->JetEta[jetIndices_[1]->idx_]) < minJetEta_ ) {
      nMinJetEta_++;
      isGoodEvt = false;
    }
    else if( std::abs(nJet_->JetEta[jetIndices_[0]->idx_]) > maxJetEta_  ||
	     std::abs(nJet_->JetEta[jetIndices_[1]->idx_]) > maxJetEta_ ) {
      nMaxJetEta_++;
      isGoodEvt = false;
    }
    else if( std::abs(deltaPhi12) < minDeltaPhi_ ) {
      nMinDeltaPhi_++;
      isGoodEvt = false;
    }
    else if( std::abs(pJ3) > maxRel3rdJetEt_*ptRef_ ) {
      nCutOn3rdJet_++;
      isGoodEvt = false;
    }
    else if( std::abs(pSJ) > 0.75*maxRelSoftJetEt_*ptRef_ ) {
      nCutOnSoftJets_++;
      isGoodEvt = false;
    }
    else if( !passesJetId(jetIndices_[0]->idx_) || !passesJetId(jetIndices_[1]->idx_) ) {
      nMaxJetHadFraction_++;
      isGoodEvt = false;
    }
    
    if( isGoodEvt ) {
      // Create Jets
      std::vector<Jet*> jets(3);
      for(size_t i = 0; i < jets.size(); ++i) {
	if( i < nJets ) {
	  double dphi = TVector2::Phi_mpi_pi(nJet_->JetPhi[jetIndices_[i]->idx_] - nJet_->GenJetPhi[jetIndices_[i]->idx_]);
	  double deta         = nJet_->JetEta[jetIndices_[i]->idx_] - nJet_->GenJetEta[jetIndices_[i]->idx_];
	  double drJetGenjet  = sqrt( deta*deta + dphi*dphi );
	  
	  jets[i] = new Jet(nJet_->JetPt[jetIndices_[i]->idx_],
			    nJet_->JetEMF[jetIndices_[i]->idx_]*nJet_->JetEt[jetIndices_[i]->idx_],
			    (1.-nJet_->JetEMF[jetIndices_[i]->idx_])*nJet_->JetEt[jetIndices_[i]->idx_],
			    0,
			    nJet_->JetE[jetIndices_[i]->idx_],
			    nJet_->JetEta[jetIndices_[i]->idx_],
			    nJet_->JetPhi[jetIndices_[i]->idx_],
			    nJet_->JetEtWeightedSigmaPhi[jetIndices_[i]->idx_],
			    nJet_->JetEtWeightedSigmaEta[jetIndices_[i]->idx_],
			    Jet::uds,
			    nJet_->GenJetPt[jetIndices_[i]->idx_],
			    drJetGenjet,
			    createCorFactors(jetIndices_[i]->idx_),
			    par_->jet_function(nJet_->JetIEta[jetIndices_[i]->idx_],
					       nJet_->JetIPhi[jetIndices_[i]->idx_]),
			    jet_error_param,
			    par_->global_jet_function());
	  // Read external correction factors
	  if(corFactorsFactory_) {
	    jets[i]->updateCorFactors(corFactorsFactory_->create(jets[i]));
	  }
	  // Correct measurement to L3 (L1*L2*L3)
	  if(correctToL3_) {
	    jets[i]->correctToL3();
	  }
	  // Correct measurement with L2*L3
	  else if(correctL2L3_) {
	    jets[i]->correctL2L3();
	  }
	} else {
	  jets[i] = 0;
	}
      }
      // Create SmearDiJet event from three jets
      dijet = new SmearDiJet(jets[0],                    // First jet
			     jets[1],                    // Second jet
			     jets[2],                    // Third jet
			     std::abs(deltaPhi12),pPhi,
			     ptJet3,ptJet4,
			     pJ3,pSJ,pUCE,ptRef_,
			     nJet_->GenEvtScale,
			     weightRelToNtuple_*nJet_->Weight,
			     par_->resolutionFitPDF(1,1),
			     min_,                       // Integration minimum
			     max_,                       // Integration maximum
			     eps_,                       // Integration step length
			     maxNIter_);                 // Integration n iterations
    }
    for(size_t i = 0; i < nJets; ++i) {
      delete jetIndices_[i];
    }
  }
  
  return dijet;
}


// ----------------------------------------------------------------   
Event* DiJetReader::createSmearEventGenOrdered() {
  SmearDiJet *dijet = 0;

  // There should be at least two gen jets in the event
  if( nJet_->NobjGenJet < 2 ) {
    nDiJetCut_++;
  } else {
    // Check if event is selected
    bool isGoodEvt = true;

    double deltaPhi12 = TVector2::Phi_mpi_pi(nJet_->JetPhi[nJet_->GenJetColJetIdx[0]]-nJet_->JetPhi[nJet_->GenJetColJetIdx[1]]);
    double pt0 = nJet_->JetCorrL2L3[nJet_->GenJetColJetIdx[0]]*nJet_->JetPt[nJet_->GenJetColJetIdx[0]];
    double pt1 = nJet_->JetCorrL2L3[nJet_->GenJetColJetIdx[1]]*nJet_->JetPt[nJet_->GenJetColJetIdx[1]];
    double ptAve = 0.5*(pt0+pt1);

    // Phi of dijet axis
    double pPhi = TVector2::Phi_mpi_pi(nJet_->GenJetColPhi[0]-0.5*deltaPhi12+M_PI/2.);
    double pJ3 = 0.;
    if( nJet_->NobjGenJet > 2 ) {
      pJ3 = nJet_->GenJetColPt[2]*cos(TVector2::Phi_mpi_pi(pPhi-nJet_->GenJetColPhi[2]));
    }
    double pSJ = 0.;
    for(int i = 3; i < nJet_->NobjGenJet; ++i) {
      pSJ += nJet_->GenJetColPt[i]*cos(TVector2::Phi_mpi_pi(pPhi-nJet_->GenJetColPhi[i]));
    }
    double pUCE = 0.;


    if( nJet_->GenJetColPt[0] < minGenJetEt_ || 
  	nJet_->GenJetColPt[1] < minGenJetEt_ ) {
      nMinGenJetEt_++;
      isGoodEvt = false;
    }
    else if( nJet_->GenJetColPt[0] > maxGenJetEt_ || 
	     nJet_->GenJetColPt[1] > maxGenJetEt_ ) {
      nMaxGenJetEt_++;
      isGoodEvt = false;
    }
    else if( ptAve < minDijetEt_ ) {
      nMinDijetEt_++;
      isGoodEvt = false;
    }
    else if( ptAve > maxDijetEt_ ) {
       nMaxDijetEt_++;
       isGoodEvt = false;
    }
    else if( std::abs(nJet_->JetEta[nJet_->GenJetColJetIdx[0]]) < minJetEta_  ||
	     std::abs(nJet_->JetEta[nJet_->GenJetColJetIdx[1]]) < minJetEta_ ) {
      nMinJetEta_++;
      isGoodEvt = false;
    }
    else if( std::abs(nJet_->JetEta[nJet_->GenJetColJetIdx[0]]) > maxJetEta_  ||
	     std::abs(nJet_->JetEta[nJet_->GenJetColJetIdx[1]]) > maxJetEta_ ) {
      nMaxJetEta_++;
      isGoodEvt = false;
    }
     else if( std::abs(deltaPhi12) < minDeltaPhi_ ) {
       nMinDeltaPhi_++;
       isGoodEvt = false;
     }
     else if( std::abs(pJ3) > maxRel3rdJetEt_*ptRef_ ) {
       nCutOn3rdJet_++;
       isGoodEvt = false;
     }
     else if( std::abs(pSJ) > 0.75*maxRelSoftJetEt_*ptRef_ ) {
       nCutOnSoftJets_++;
       isGoodEvt = false;
     }
     else if( 1. - nJet_->JetEMF[nJet_->GenJetColJetIdx[0]] < minJetHadFraction_ ||
	      1. - nJet_->JetEMF[nJet_->GenJetColJetIdx[1]] < minJetHadFraction_ ) {
       nMinJetHadFraction_++;
       isGoodEvt = false;
     }
     else if( 1. - nJet_->JetEMF[nJet_->GenJetColJetIdx[0]] > maxJetHadFraction_ ||
	      1. - nJet_->JetEMF[nJet_->GenJetColJetIdx[1]] > maxJetHadFraction_ ) {
       nMaxJetHadFraction_++;
       isGoodEvt = false;
     }
    
    if( isGoodEvt ) {
      // Create Jets
      std::vector<Jet*> jets(3);
      
      for(size_t i = 0; i < jets.size(); ++i) {
	if( i < static_cast<size_t>(nJet_->NobjGenJet) ) {
	  double dphi = TVector2::Phi_mpi_pi(nJet_->JetPhi[nJet_->GenJetColJetIdx[i]]-nJet_->GenJetColPhi[i]);
	  double deta         = nJet_->JetEta[nJet_->GenJetColJetIdx[i]] - nJet_->GenJetColEta[i];
	  double drJetGenjet  = sqrt( deta*deta + dphi*dphi );
	  
	  jets[i] = new Jet(nJet_->JetPt[nJet_->GenJetColJetIdx[i]],
			    nJet_->JetEMF[nJet_->GenJetColJetIdx[i]]*nJet_->JetEt[nJet_->GenJetColJetIdx[i]],
			    (1.-nJet_->JetEMF[nJet_->GenJetColJetIdx[i]])*nJet_->JetEt[nJet_->GenJetColJetIdx[i]],
			    0,
			    nJet_->JetE[nJet_->GenJetColJetIdx[i]],
			    nJet_->JetEta[nJet_->GenJetColJetIdx[i]],
			    nJet_->JetPhi[nJet_->GenJetColJetIdx[i]],
			    nJet_->JetEtWeightedSigmaPhi[nJet_->GenJetColJetIdx[i]],
			    nJet_->JetEtWeightedSigmaEta[nJet_->GenJetColJetIdx[i]],
			    Jet::uds,
			    nJet_->GenJetColPt[i],
			    drJetGenjet,
			    createCorFactors(nJet_->GenJetColJetIdx[i]),
			    par_->jet_function(nJet_->JetIEta[nJet_->GenJetColJetIdx[i]],
					       nJet_->JetIPhi[nJet_->GenJetColJetIdx[i]]),
			    jet_error_param,
			    par_->global_jet_function());
	  // Read external correction factors
	  if(corFactorsFactory_) {
	    jets[i]->updateCorFactors(corFactorsFactory_->create(jets[i]));
	  }
	  // Correct measurement to L3 (L1*L2*L3)
	  if(correctToL3_) {
	    jets[i]->correctToL3();
	  }
	  // Correct measurement with L2*L3
	  else if(correctL2L3_) {
	    jets[i]->correctL2L3();
	  }
	} else {
	  jets[i] = 0;
	}
      }
      if( jets[0]->dR() > maxDeltaR_ || jets[1]->dR() > maxDeltaR_ ) {
	nMaxDeltaR_++;
	isGoodEvt = false;
      } else {
 	// Create SmearDiJet event
	dijet = new SmearDiJet(jets[0],                    // First jet
			       jets[1],                    // Second jet
			       jets[2],                    // Third jet
			       std::abs(deltaPhi12),pPhi,
			       (nJet_->NobjGenJet<2 ? 0. : nJet_->GenJetColPt[2]),
			       (nJet_->NobjGenJet<3 ? 0. : nJet_->GenJetColPt[3]),
			       pJ3,pSJ,pUCE,ptRef_,
			       nJet_->GenEvtScale,
			       weightRelToNtuple_*nJet_->Weight,
			       par_->resolutionFitPDF(1,1),
			       min_,                       // Integration minimum
			       max_,                       // Integration maximum
			       eps_,                       // Integration step length
			       maxNIter_);                 // Integration n iterations
      }
    }
  }
  
  return dijet;
}



//!  \brief Create an event of type \p TwoJetsPtBalanceEvent
//!
//!  The two jets leading in uncorrected calo pt are randomly
//!  assigned jet 1 and 2.
// ----------------------------------------------------------------   
TwoJetsPtBalanceEvent* DiJetReader::createTwoJetsPtBalanceEvent()
{
  // Number of jets in the event
  int nJets = nJet_->NobjJet;

  // There should be at least two jets
  // and at the most three jets
  if( nJets < 2 ) {
  nDiJetCut_++;
  return 0;
  } else if(nJets > 3) {
    nJets = 3;
  }
  if( std::abs(TVector2::Phi_mpi_pi(nJet_->JetPhi[0] - nJet_->JetPhi[1])) < minDeltaPhi_ ) {
    nMinDeltaPhi_++;
    return 0;
  }
  if ( nJet_->GenJetPt[0] < minGenJetEt_ || nJet_->GenJetPt[1] < minGenJetEt_ ) {
    nMinGenJetEt_++;
    return 0;
  }
  if( nJet_->GenJetPt[0] > maxGenJetEt_ || nJet_->GenJetPt[1] > maxGenJetEt_ ) {
    nMaxGenJetEt_++;
    return 0;
  }  
  double dphi        = TVector2::Phi_mpi_pi( nJet_->JetPhi[0] - nJet_->GenJetPhi[0] );
  double deta        = nJet_->JetEta[0] - nJet_->GenJetColEta[0];
  double dR1 = sqrt( deta*deta + dphi*dphi );
  dphi        = TVector2::Phi_mpi_pi( nJet_->JetPhi[1] - nJet_->GenJetPhi[1] );
  deta        = nJet_->JetEta[1] - nJet_->GenJetColEta[1];
  double dR2 = sqrt( deta*deta + dphi*dphi );
  if( dR1 > maxDeltaR_ || dR2 > maxDeltaR_ ) {
    nMaxDeltaR_++;
    return 0;
  }
  if( (nJet_->JetPt[0] < minJetEt_) || (nJet_->JetPt[1] < minJetEt_) ) {
    nMinJetEt_++;
    return 0;
  }
  if( nJet_->JetPt[0] > maxJetEt_ || nJet_->JetPt[1] > maxJetEt_ ) {
    nMaxJetEt_++;
    return 0;
  }
  if( std::abs(nJet_->JetEta[0]) > maxJetEta_ || std::abs(nJet_->JetEta[1]) > maxJetEta_ ) {
    nMaxJetEta_++;
    return 0;
  }
  if(1 - nJet_->JetEMF[0] < minJetHadFraction_ || 
	  1 - nJet_->JetEMF[1] < minJetHadFraction_ ) {
    nMinJetHadFraction_++;
    return 0;
  }
  if( 1 - nJet_->JetEMF[0] > maxJetHadFraction_ || 
	   1 - nJet_->JetEMF[1] > maxJetHadFraction_ ) {
    nMaxJetHadFraction_++;
    return 0;
  }
  //loose jet id 
  /* 
     if(! (nJet_->JetIDLoose[0] && nJet_->JetIDLoose[1])) {
     nMaxJetHadFraction_++;
    return 0;
  }
  */
  if(((( nJet_->JetEMF[0] <= 0.01) && (std::abs(nJet_->JetEta[0]) < 2.6) )) ||
     (( nJet_->JetEMF[1] <= 0.01) && (std::abs(nJet_->JetEta[1]) < 2.6) )) {
    nMaxJetHadFraction_++;
    return 0;
  }

  
  if( nJet_->JetN90Hits[0] <= 1 || nJet_->JetN90Hits[1] <= 1) {
    nMaxJetHadFraction_++;
    return 0;
  }
  if( nJet_->JetFHPD[0] >= 0.98 || nJet_->JetFHPD[1] >= 0.98) {
    nMaxJetHadFraction_++;
    return 0;
  }
  
  double ptAve = 0.5 * (nJet_->JetPt[0]+ nJet_->JetPt[1]);
  if( nJets > 2) {
    //compute dijet kin
    double deltaPhi12 = TVector2::Phi_mpi_pi(nJet_->JetPhi[0]-nJet_->JetPhi[1]);
    
    // Phi of dijet axis
    double pPhi = TVector2::Phi_mpi_pi(nJet_->JetPhi[0]-0.5*deltaPhi12+M_PI/2.);
    double pJ3 = 0.;
    double ptJet3 = 0.;
    if( nJet_->NobjJet > 2 ) {
      pJ3 = nJet_->JetPt[2]*cos(TVector2::Phi_mpi_pi(pPhi-nJet_->JetPhi[2]));
      ptJet3 = nJet_->JetPt[2];
    }
    double pSJ = 0.;
    for(int i = 3; i < nJet_->NobjJet; ++i) {
      pSJ += nJet_->JetPt[i]*cos(TVector2::Phi_mpi_pi(pPhi-nJet_->JetPhi[i]));
    }
    //double pUCE = 0.;
    
    if( std::abs(pJ3) > maxRel3rdJetEt_*ptAve ) {
      nCutOn3rdJet_++;
      return 0;
    }
    if( std::abs(pSJ) > 0.75*maxRelSoftJetEt_*ptAve ) {
      nCutOnSoftJets_++;
      return 0;
    }
  }
  
  // Pointer to the three jets leading in pt calo
  Jet * jet1 = 0;
  Jet * jet2 = 0;
  Jet * jet3 = 0;

  // Loop over the two or, if existing, three jets
  // leading in calo pt; the first two jets are
  // assigned randomly to have an unbiased balance
  // measure
  int calJetIdx[3];
  calJetIdx[0] = 0;
  calJetIdx[1] = 1;
  calJetIdx[2] = 2;
  for(int i = 0; i < nJets; i++) {
    // Jet - GenJet matching
    double dphi         = TVector2::Phi_mpi_pi( nJet_->JetPhi[calJetIdx[i]] - nJet_->GenJetPhi[calJetIdx[i]] );
    double deta         = nJet_->JetEta[calJetIdx[i]] - nJet_->GenJetEta[calJetIdx[i]];
    double drJetGenjet  = sqrt( deta*deta + dphi*dphi );
    
    // Create jet
    Jet * jet = new Jet(nJet_->JetPt[calJetIdx[i]],
			nJet_->JetPt[calJetIdx[i]] * nJet_->JetEMF[calJetIdx[i]],
			nJet_->JetPt[calJetIdx[i]] * (1.0 - nJet_->JetEMF[calJetIdx[i]]),
			0.0,
			nJet_->JetE[calJetIdx[i]],
			nJet_->JetEta[calJetIdx[i]],
			nJet_->JetPhi[calJetIdx[i]], 
			nJet_->JetEtWeightedSigmaPhi[calJetIdx[i]],
			nJet_->JetEtWeightedSigmaEta[calJetIdx[i]],
			Jet::uds,
			nJet_->GenJetPt[calJetIdx[i]],
			drJetGenjet,createCorFactors(calJetIdx[i]),
			par_->jet_function(nJet_->JetIEta[calJetIdx[i]],
					   nJet_->JetIPhi[calJetIdx[i]]),
			jet_error_param,
			par_->global_jet_function());    
    
    // Store jet
    if( i == 0 ) {
      jet1 = jet;
    } else if( i == 1 ) {
      jet2 = jet;
    } else if( i == 2 ) {
      jet3 = jet;
    }
  }  // End of loop over jets 
  if(corFactorsFactory_) {
    jet1->updateCorFactors(corFactorsFactory_->create(jet1));
    jet2->updateCorFactors(corFactorsFactory_->create(jet2));
    if(jet3) jet3->updateCorFactors(corFactorsFactory_->create(jet3));    
  }
  // Correct measurement to L3 (L1*L2*L3)
  if(correctToL3_) {
    jet1->correctToL3();
    jet2->correctToL3();
    if(jet3) jet3->correctToL3();
  }
  // Correct measurement with L2*L3
  if(correctL2L3_) {
    jet1->correctL2L3();
    jet2->correctL2L3();
    if(jet3) jet3->correctL2L3();
  }
  // Create TwoJetsInvMassEvent
  TwoJetsPtBalanceEvent * evt
    = new TwoJetsPtBalanceEvent(jet1,jet2,jet3,nJet_->GenEvtScale,nJet_->Weight) ;
  
  return evt;
}


// ----------------------------------------------------------------   
CorFactors* DiJetReader::createCorFactors(int jetid) const 
{
  return new CorFactors(nJet_->JetCorrZSP[jetid], // L1
			nJet_->JetCorrL2[jetid],  // L2
			nJet_->JetCorrL3[jetid],  // L3
			nJet_->JetCorrL4JW[jetid],  // L4
			1.,                         // L5
			nJet_->JetCorrJPT[jetid],
			nJet_->JetCorrL2L3JPT[jetid]); //JPTL2L3
}



// ----------------------------------------------------------------   
bool DiJetReader::passesJetId(int idx) const {
  bool ok = true;
  if( nJet_->JetN90Hits[idx] < minJetN90Hits_ )
    ok = false;
  else if( nJet_->JetFHPD[idx] > maxJetFHPD_ )
    ok = false;
  else if( (std::abs(nJet_->JetEta[idx]) < 2.55 && nJet_->JetEMF[idx] < 0.01) )
    ok = false;
  else if( (std::abs(nJet_->JetEta[idx]) > 2.55 && nJet_->JetEMF[idx] < -0.9) ||
	   (std::abs(nJet_->JetEta[idx]) > 2.55 && nJet_->JetPt[idx] > 80. && nJet_->JetEMF[idx] < -0.9) )
    ok = false;

  return ok;
}

boost::mutex DiJetReader::dijetmutex;
