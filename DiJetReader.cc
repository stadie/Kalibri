//
//    first version: Hartmut Stadie 2008/12/12
//    $Id: DiJetReader.cc,v 1.54 2010/07/22 17:38:32 mschrode Exp $
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



//!  \brief Constructor
//!
//!  Reads data from ROOT trees and stores them in an \p NJetSel selector.
//!  The data can be stored in a format derived from \p Event (as specified
//!  in the 'Di-Jet data class' field in the config file) by calling the
//!  method readEvents(<tt>std::vector<Event*>& data</tt>). Additionally,
//!  the cut thresholds are read from the configfile.
//!
//!  \param configfile Name of configfile
//!  \param p Pointer to \p TParameters object
// ----------------------------------------------------------------   
DiJetReader::DiJetReader(const std::string& configfile, TParameters* p)
  : EventReader(configfile,p), nJet_(new NJetSel()), rand_(new TRandom3(0)),
    zero_(0)
{
  // Maximum number of read events
  nDijetEvents_ = config_->read<int>("use Di-Jet events",-1);
  if(nDijetEvents_ == 0) return;
  prescale_ = config_->read<int>("Di-Jet prescale",1);

  // Cuts
  minJetEt_          = config_->read<double>("Et min cut on jet",0.0);
  maxJetEt_          = config_->read<double>("Et max cut on jet",10000.0);
  minDijetEt_        = config_->read<double>("Et min cut on dijet",0.0); 
  maxDijetEt_        = config_->read<double>("Et max cut on dijet",10000.0); 
  max3rdJetEt_       = config_->read<double>("Et cut on n+1 Jet",10.0);
  maxRel3rdJetEt_    = config_->read<double>("Relative n+1 Jet Et Cut",0.2);
  minJetEta_         = config_->read<double>("Eta min cut on jet",0.0);
  maxJetEta_         = config_->read<double>("Eta max cut on jet",5.0);
  minJetHadFraction_ = config_->read<double>("Min had fraction",0.07);
  maxJetHadFraction_ = config_->read<double>("Max had fraction",0.95);
  minDeltaPhi_       = config_->read<double>("Min Delta Phi",2.5);
  minGenJetEt_       = config_->read<double>("Et genJet min",0.0);
  maxGenJetEt_       = config_->read<double>("Et genJet max",10000.0);
  maxDeltaR_         = config_->read<double>("DeltaR cut on jet matching",0.25);
  // Counter for cutflow
  nDiJetCut_          = 0;
  nMinJetEt_          = 0;
  nMaxJetEt_          = 0;
  nMinDijetEt_        = 0;
  nMaxDijetEt_        = 0;
  nCutOn3rdJet_       = 0;
  nMinJetEta_         = 0;  
  nMaxJetEta_         = 0;  
  nMinDeltaPhi_       = 0;
  nMinJetHadFraction_ = 0;     
  nMaxJetHadFraction_ = 0;     
  nMinGenJetEt_       = 0;    
  nMaxGenJetEt_       = 0;     
  nMaxDeltaR_         = 0;
  // Integration parameter for SmearData
  maxNIter_          = config_->read<int>("DiJet integration number of iterations",5);
  eps_               = config_->read<double>("DiJet integration epsilon",1.E-5);
  min_               = config_->read<double>("DiJet integration min",0.);
  max_               = config_->read<double>("DiJet integration max",1.);
  // Data class
  dataClass_ = config_->read<int>("Di-Jet data class", 0);
  if((dataClass_ != 1)&&(dataClass_ != 11)&&(dataClass_ != 12)&&(dataClass_ != 21)&&(dataClass_ != 5)) {
    std::cerr << "DiJetReader: Unknown data class " << dataClass_ << ". Aborting!" << std::endl;
    exit(9);
  }
  // Input files
  nJet_->Init(createTree("dijet")); 

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
  if(nDijetEvents_ == 0) return 0;

  // Reset counters of rejected events
  nDiJetCut_          = 0;
  nMinJetEt_          = 0;
  nMaxJetEt_          = 0;
  nMinDijetEt_        = 0;
  nMaxDijetEt_        = 0;
  nCutOn3rdJet_       = 0;
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
  int nReadEvts = 0;                          // Number of read events
  int nGoodEvts = 0;                          // Number of events passing all cuts

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

  if((dataClass_ == 11)||(dataClass_ == 21)) { 
    nJet_->fChain->SetBranchStatus("Track*",0);
    //nJet_->fChain->SetBranchStatus("Tow*",0);
  } else if(dataClass_ == 12) {
    nJet_->fChain->SetBranchStatus("Track*",0);
  } else if(dataClass_ == 5) {
    nJet_->fChain->SetBranchStatus("Track*",0);
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
      nReadEvts++;
      TwoJetsPtBalanceEvent* td = createTwoJetsPtBalanceEvent(); 
      if(td) {
	//seconed jet should be central...
	if(std::abs(td->getJet1()->eta()) < 1.3) {
	  data.push_back(new TwoJetsPtBalanceEvent(td->getJet2()->clone(),td->getJet1()->clone(),
						   td->getJet3() ? td->getJet3()->clone():0,
						   td->ptHat(),td->weight()));
	  ++nGoodEvts;
	}
	if(std::abs(td->getJet2()->eta()) < 1.3) {
	  data.push_back(td ); 
	  ++nGoodEvts;
	} else {
	  delete td;
	}
      } 
    } else if((dataClass_ == 11)  || (dataClass_ == 12) || (dataClass_ == 21)) {
      nReadEvts++;
      int nAddedJets = createJetTruthEvents(data);
      if( nAddedJets ) nGoodEvts += nAddedJets;    
    } else if(dataClass_ == 5) {
      for(int calls = 0; calls < 2; calls++) {
	nReadEvts++;
	Event* td = createSmearEvent(calls); 
	if(td) {
	  nGoodEvts++;
	  data.push_back(td ); 
	} 
      }
    }
    if(nReadEvts>=nDijetEvents_ && nDijetEvents_>=0 ) break;
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
      }      
      delete ijb->second;
    }
    binning_->clear();
  }
  // Print cut flow
  std::cout << "Read " << nReadEvts << " dijet events:\n";
  if( dataClass_ == 1 ) {
    std::cout << "  " << (nReadEvts-=nDiJetCut_) << std::flush;
    std::cout << " events with 2 or more jets\n"; 
    std::cout << "  " << (nReadEvts-=nMinDeltaPhi_) << std::flush;
    std::cout << " dijet events with DeltaPhi > " << minDeltaPhi_ << '\n';;
    std::cout << "  " << (nReadEvts-=nMinGenJetEt_) << std::flush;
    std::cout << " dijet events with ptgen > " << minGenJetEt_ << " GeV\n";
    std::cout << "  " << (nReadEvts-=nMaxGenJetEt_) << std::flush;
    std::cout << " dijet events with ptgen < " << maxGenJetEt_ << " GeV\n";
    std::cout << "  " << (nReadEvts-=nMaxDeltaR_) << std::flush;
    std::cout << " dijet events with DeltaR < " << maxDeltaR_ << "\n";
    std::cout << "  " << (nReadEvts-=nMinJetEt_) << std::flush;
    std::cout << " dijet events Et > " << minJetEt_ << " GeV\n";
    std::cout << "  " << (nReadEvts-=nMaxJetEt_) << std::flush;
    std::cout << " dijet events Et < " << maxJetEt_ << " GeV\n";
    std::cout << "  " << (nReadEvts-=(nMinJetEta_+nMaxJetEta_)) << std::flush;
    std::cout << " dijet events with " << minJetEta_ << " < |eta| < " << maxJetEta_ << "\n";
    std::cout << "  " << (nReadEvts-=nMinJetHadFraction_) << std::flush;
    std::cout << " dijet events with hadronic fraction > " << minJetHadFraction_ << "\n";
    std::cout << "  " << (nReadEvts-=nMaxJetHadFraction_) << std::flush;
    std::cout << " dijet events with jet id and hadronic fraction < " << maxJetHadFraction_ << "\n";
    std::cout << "  " << (nReadEvts-=nCutOn3rdJet_) << std::flush;
    std::cout << " dijet events with pt(jet3) / pt(dijet) < " << maxRel3rdJetEt_ << " or ";
    std::cout << "pt(jet3) < " << max3rdJetEt_ << " GeV\n";
  } else if( dataClass_ == 11 || dataClass_ == 12 || dataClass_ == 21 ) {
    std::cout << "  " << (nReadEvts-=nDiJetCut_) << std::flush;
    std::cout << " events with 2 or more jets\n";
    std::cout << "  That are " << (nReadEvts*=2) << " jet-truth events:\n";
    std::cout << "    " << (nReadEvts-=nMinGenJetEt_) << std::flush;
    std::cout << " jet-truth events with ptgen > " << minGenJetEt_ << "\n";
    std::cout << "    " << (nReadEvts-=nMaxGenJetEt_) << std::flush;
    std::cout << " jet-truth events with ptgen < " << maxGenJetEt_ << "\n";
    std::cout << "    " << (nReadEvts-=nMaxDeltaR_) << std::flush;
    std::cout << " jet-truth events with DeltaR < " << maxDeltaR_ << "\n";
    std::cout << "    " << (nReadEvts-=nMinJetEt_) << std::flush;
    std::cout << " jet-truth events Et > " << minJetEt_ << "\n";
    std::cout << "    " << (nReadEvts-=nMaxJetEt_) << std::flush;
    std::cout << " jet-truth events Et < " << maxJetEt_ << " GeV\n";
    std::cout << "  " << (nReadEvts-=(nMinJetEta_+nMaxJetEta_)) << std::flush;
    std::cout << " dijet events with " << minJetEta_ << " < |eta| < " << maxJetEta_ << "\n";
    std::cout << "    " << (nReadEvts-=nMinJetHadFraction_) << std::flush;
    std::cout << " jet-truth events with hadronic fraction > " << minJetHadFraction_ << "\n";
    std::cout << "    " << (nReadEvts-=nMaxJetHadFraction_) << std::flush;
    std::cout << " jet-truth events with hadronic fraction < " << maxJetHadFraction_ << "\n";
  } else if( dataClass_ == 5 ) {
    std::cout << "  " << (nReadEvts-=nDiJetCut_) << std::flush;
    std::cout << " events with more than 3 jets\n";
    std::cout << "  " << (nReadEvts-=nMinGenJetEt_) << std::flush;
    std::cout << " dijet events with ptgen > " << minGenJetEt_ << " GeV\n";
    std::cout << "  " << (nReadEvts-=nMaxGenJetEt_) << std::flush;
    std::cout << " dijet events with ptgen < " << maxGenJetEt_ << " GeV\n";
    std::cout << "  " << (nReadEvts-=nMinJetEt_) << std::flush;
    std::cout << " dijet events Et > " << minJetEt_ << " GeV\n";
    std::cout << "  " << (nReadEvts-=nMaxJetEt_) << std::flush;
    std::cout << " dijet events Et < " << maxJetEt_ << " GeV\n";
    std::cout << "  " << (nReadEvts-=(nMinJetEta_+nMaxJetEta_)) << std::flush;
    std::cout << " dijet events with " << minJetEta_ << " < |eta| < " << maxJetEta_ << "\n";
    std::cout << "  " << (nReadEvts-=nCutOn3rdJet_) << std::flush;
    std::cout << " dijet events with pt(jet3) / pt(dijet) < " << maxRel3rdJetEt_ << "\n";
    std::cout << "  " << (nReadEvts-=nMinJetHadFraction_) << std::flush;
    std::cout << " dijet events with hadronic fraction > " << minJetHadFraction_ << "\n";
    std::cout << "  " << (nReadEvts-=nMaxJetHadFraction_) << std::flush;
    std::cout << " dijet events with hadronic fraction < " << maxJetHadFraction_ << "\n";
    std::cout << "  " << (nReadEvts-=nMaxDeltaR_) << std::flush;
    std::cout << " dijet events with DeltaR < " << maxDeltaR_ << "\n";
  }
  std::cout << "Stored " << nGoodEvts << " dijet events for analysis.\n";
  return nGoodEvts;
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
    double err2 = 0;
    Measurement tower;
    double dR        = 10;
    int closestTower = 0;
    for(int n=0; n<nJet_->NobjTow; ++n){
      if(nJet_->Tow_jetidx[n] != calJetIdx) continue;//look for ij-jet's towers
      if( dataClass_ == 12) {
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
	  terr[n] = TParameters::toy_tower_error_parametrization(&tower.pt,&tower);
	}
	terr[n]  *= terr[n];
	err2     += terr[n];
      }
      dphi      = TVector2::Phi_mpi_pi(nJet_->JetPhi[calJetIdx]-tower.phi);
      deta      = nJet_->JetEta[calJetIdx]-tower.eta;
      double dr = sqrt( deta*deta + dphi*dphi );     
      if(dr < dR) {
	dR = dr;
	closestTower = n;
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
			  par_->jet_function(nJet_->TowId_eta[closestTower],
					     nJet_->TowId_phi[closestTower]),
			  jet_error_param,par_->global_jet_function(),minJetEt_);
      for(int j = 0 ; j < nJet_->NobjTow ; ++j) {
	if (nJet_->Tow_jetidx[j]!= calJetIdx) continue;//look for ij-jet's towers
	double scale = nJet_->TowEt[j]/nJet_->TowE[j];
	jt->addTower(nJet_->TowEt[j],nJet_->TowEm[j]*scale,
		     nJet_->TowHad[j]*scale,nJet_->TowOE[j]*scale,
		     nJet_->TowE[j],nJet_->TowEta[j],nJet_->TowPhi[j],
		     par_->tower_function(nJet_->TowId_eta[calJetIdx],nJet_->TowId_phi[calJetIdx]),
		     tower_error_param);
      }
      jet = jt;
    }
    else { 
      jet = new Jet(nJet_->JetPt[calJetIdx],tower.EMF,tower.HadF,
		    tower.OutF,nJet_->JetE[calJetIdx],
		    nJet_->JetEta[calJetIdx],nJet_->JetPhi[calJetIdx],
		    0.,//nJet_->JetEtWeightedSigmaPhi[calJetIdx],
		    0.,//nJet_->JetEtWeightedSigmaEta[calJetIdx],
		    Jet::flavorFromPDG(nJet_->GenPartId_algo[calJetIdx]),
		    nJet_->GenJetColPt[genJetIdx],drJetGenjet,
		    createCorFactors(calJetIdx),
		    par_->jet_function(1,//nJet_->TowId_eta[closestTower],
				       1),//nJet_->TowId_phi[closestTower]),
		    jet_error_param,par_->global_jet_function(),minJetEt_);    
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
  //return createSmearEventGenOrdered(callIdx);
  return createSmearEventCaloOrdered(callIdx);
}


// ----------------------------------------------------------------   
Event* DiJetReader::createSmearEventCaloOrdered(int callIdx)
{
  SmearDiJet *dijet = 0;

  // There should be at least three calo jets in the event
  if( nJet_->NobjJet < 3 ) {
    nDiJetCut_++;
  } else {
    // Find indices of 3 leading L2L3 corrected calo jets
    for(size_t i = 0; i < jetIndices_.size(); ++i) {
      jetIndices_[i] = new Jet::JetIndex(i,nJet_->JetPt[i]*nJet_->JetCorrL2[i]*nJet_->JetCorrL3[i]);
    }
    std::sort(jetIndices_.begin(),jetIndices_.end(),Jet::JetIndex::ptGreaterThan);
    if( callIdx == 1 ) std::swap(jetIndices_[0],jetIndices_[1]);
    jetIndices_.erase(jetIndices_.begin()+3,jetIndices_.end());
    
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
    else if( jetIndices_[0]->pt_ < minJetEt_ ) {
      nMinJetEt_++;
      isGoodEvt = false;
    }
    else if( jetIndices_[0]->pt_ > maxJetEt_ ) {
      nMaxJetEt_++;
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
    else if( 2.*jetIndices_[2]->pt_ / (jetIndices_[0]->pt_+jetIndices_[1]->pt_ ) > maxRel3rdJetEt_ ) {
      nCutOn3rdJet_++;
      isGoodEvt = false;
    }
    else if( 1. - nJet_->JetEMF[jetIndices_[0]->idx_] < minJetHadFraction_ ||
	     1. - nJet_->JetEMF[jetIndices_[1]->idx_] < minJetHadFraction_ ) {
      nMinJetHadFraction_++;
      isGoodEvt = false;
    }
    else if( 1. - nJet_->JetEMF[jetIndices_[0]->idx_] > maxJetHadFraction_ ||
	     1. - nJet_->JetEMF[jetIndices_[1]->idx_] > maxJetHadFraction_ ) {
      nMaxJetHadFraction_++;
      isGoodEvt = false;
    }
  
    if( isGoodEvt ) {
      // Create Jets
      std::vector<Jet*> jets(3);
      for(size_t i = 0; i < jets.size(); ++i) {
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
      }
    
      if( jets[0]->dR() > maxDeltaR_ || jets[1]->dR() > maxDeltaR_ ) {
	nMaxDeltaR_++;
	isGoodEvt = false;
      } else {
	// Create SmearDiJet event from three jets
	dijet = new SmearDiJet(jets[0],                    // First jet
			       jets[1],                    // Second jet
			       jets[2],                    // Third jet
			       nJet_->GenEvtScale,
			       nJet_->Weight,
			       par_->resolutionFitPDF(1,1),
			       min_,                       // Integration minimum
			       max_,                       // Integration maximum
			       eps_,                       // Integration step length
			       maxNIter_);                 // Integration n iterations
      }
    }
    for(size_t i = 0; i < jetIndices_.size(); ++i) {
      delete jetIndices_[i];
    }
  }
  
  return dijet;
}


// ----------------------------------------------------------------   
Event* DiJetReader::createSmearEventGenOrdered(int callIdx) {
  SmearDiJet *dijet = 0;

  // There should be at least three gen jets in the event
  if( nJet_->NobjGenJet < 3 ) {
    nDiJetCut_++;
  } else {
    std::vector<int> calJetIdx(3);
    calJetIdx[0] = nJet_->GenJetColJetIdx[callIdx==0 ? 0 : 1];
    calJetIdx[1] = nJet_->GenJetColJetIdx[callIdx==0 ? 0 : 1];
    calJetIdx[2] = nJet_->GenJetColJetIdx[2];
    // Check if event is selected
    bool isGoodEvt = true;
    if( nJet_->GenJetColPt[0] < minGenJetEt_ || 
	nJet_->GenJetColPt[1] < minGenJetEt_ ) {
      nMinGenJetEt_++;
      isGoodEvt = false;
    }
    else if( nJet_->GenJetColPt[1] > maxGenJetEt_ || 
	     nJet_->GenJetColPt[1] > maxGenJetEt_ ) {
      nMaxGenJetEt_++;
      isGoodEvt = false;
    }
    else if( nJet_->JetCorrL2[calJetIdx[0]]*nJet_->JetCorrL3[calJetIdx[0]]*nJet_->JetPt[calJetIdx[0]] < minJetEt_ ) {
      nMinJetEt_++;
      isGoodEvt = false;
    }
    else if( nJet_->JetCorrL2[calJetIdx[0]]*nJet_->JetCorrL3[calJetIdx[0]]*nJet_->JetPt[calJetIdx[0]] > maxJetEt_ ) {
      nMaxJetEt_++;
      isGoodEvt = false;
    }
    else if( std::abs(nJet_->JetEta[calJetIdx[0]]) < minJetEta_  ||
	     std::abs(nJet_->JetEta[calJetIdx[1]]) < minJetEta_ ) {
      nMinJetEta_++;
      isGoodEvt = false;
    }
    else if( std::abs(nJet_->JetEta[calJetIdx[1]]) > maxJetEta_  ||
	     std::abs(nJet_->JetEta[calJetIdx[1]]) > maxJetEta_ ) {
      nMaxJetEta_++;
      isGoodEvt = false;
    }
    else if( 2.*nJet_->JetCorrL2[calJetIdx[2]]*nJet_->JetCorrL3[calJetIdx[2]]*nJet_->JetPt[calJetIdx[2]] /
	     (nJet_->JetCorrL2[calJetIdx[0]]*nJet_->JetCorrL3[calJetIdx[0]]*nJet_->JetPt[calJetIdx[0]]+
	      nJet_->JetCorrL2[calJetIdx[1]]*nJet_->JetCorrL3[calJetIdx[1]]*nJet_->JetPt[calJetIdx[1]])
	     > maxRel3rdJetEt_ ) {
      nCutOn3rdJet_++;
      isGoodEvt = false;
    }
    else if( 1. - nJet_->JetEMF[calJetIdx[0]] < minJetHadFraction_ ||
	     1. - nJet_->JetEMF[calJetIdx[1]] < minJetHadFraction_ ) {
      nMinJetHadFraction_++;
      isGoodEvt = false;
    }
    else if( 1. - nJet_->JetEMF[calJetIdx[0]] > maxJetHadFraction_ ||
	     1. - nJet_->JetEMF[calJetIdx[1]] > maxJetHadFraction_ ) {
      nMaxJetHadFraction_++;
      isGoodEvt = false;
    }
  
    if( isGoodEvt ) {
      // Create Jets
      std::vector<Jet*> jets(3);
      std::vector<int> genJetIdx(3);
      genJetIdx[0] = callIdx==0 ? 0 : 1;
      genJetIdx[1] = callIdx==0 ? 0 : 1;
      genJetIdx[2] = 2;

      for(size_t i = 0; i < jets.size(); ++i) {
	double dphi = TVector2::Phi_mpi_pi(nJet_->JetPhi[calJetIdx[i]] - nJet_->GenJetColPhi[genJetIdx[i]]);
	double deta         = nJet_->JetEta[calJetIdx[i]] - nJet_->GenJetColEta[genJetIdx[i]];
	double drJetGenjet  = sqrt( deta*deta + dphi*dphi );
    
	jets[i] = new Jet(nJet_->JetPt[calJetIdx[i]],
			  nJet_->JetEMF[calJetIdx[i]]*nJet_->JetEt[calJetIdx[i]],
			  (1.-nJet_->JetEMF[calJetIdx[i]])*nJet_->JetEt[calJetIdx[i]],
			  0,
			  nJet_->JetE[calJetIdx[i]],
			  nJet_->JetEta[calJetIdx[i]],
			  nJet_->JetPhi[calJetIdx[i]],
			  nJet_->JetEtWeightedSigmaPhi[calJetIdx[i]],
			  nJet_->JetEtWeightedSigmaEta[calJetIdx[i]],
			  Jet::uds,
			  nJet_->GenJetColPt[genJetIdx[i]],
			  drJetGenjet,
			  createCorFactors(calJetIdx[i]),
 			  par_->jet_function(nJet_->JetIEta[calJetIdx[i]],
 					     nJet_->JetIPhi[calJetIdx[i]]),
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
      }
    
      if( jets[0]->dR() > maxDeltaR_ || jets[1]->dR() > maxDeltaR_ ) {
	nMaxDeltaR_++;
	isGoodEvt = false;
      } else {
	// Create SmearDiJet event from three jets
	dijet = new SmearDiJet(jets[0],                    // First jet
			       jets[1],                    // Second jet
			       jets[2],                    // Third jet
			       nJet_->GenEvtScale,
			       nJet_->Weight,
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
  } else if( nJets >= 3 ) {
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
  if( nJets > 2 &&(  2 * nJet_->JetPt[2]/(nJet_->JetPt[0] + nJet_->JetPt[1]) > maxRel3rdJetEt_ ||  nJet_->JetPt[2] > max3rdJetEt_) ) {
    nCutOn3rdJet_++;
    return 0;
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
			par_->global_jet_function(),
			minJetEt_ );    
    
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
			1.,                         // L4
			1.,                         // L5
			nJet_->JetCorrJPT[jetid],
			nJet_->JetCorrL2L3JPT[jetid]); //JPTL2L3
}
