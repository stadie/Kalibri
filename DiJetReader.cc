//
//    first version: Hartmut Stadie 2008/12/12
//    $Id: DiJetReader.cc,v 1.21 2009/07/23 13:49:55 mschrode Exp $
//   
#include "DiJetReader.h"

#include "CalibData.h"
#include "SmearDiJet.h"
#include "ConfigFile.h"
#include "ToyMC.h"
#include "Parameters.h"
#include "JetTruthEvent.h"
#include "Jet.h"
#include "JetWithTowers.h"
#include "TVector2.h"

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <fstream>



//!  \brief Constructor
//!
//!  Reads data from ROOT trees and stores them in an NJetSel selector.
//!  The data can be stored in a format derived from TData (as specified
//!  in the 'Di-Jet data class' field in the config file) by calling the
//!  method readEvents(std::vector<TData*>& data). Additionally, the cut
//!  thresholds are read from the configfile.
//!
//!  \param configfile Name of configfile
//!  \param p Pointer to TParameters object
// ----------------------------------------------------------------   
DiJetReader::DiJetReader(const std::string& configfile, TParameters* p)
  : EventReader(configfile,p)
{
  // Maximum number of read events
  nDijetEvents_ = config->read<int>("use Di-Jet events",-1);
  if(nDijetEvents_ == 0) {
    delete config;
    config = 0;
    return;
  }

  // Cuts
  minJetEt_          = config->read<double>("Et cut on jet",0.0);
  minDijetEt_        = config->read<double>("Et min cut on dijet",0.0); 
  maxDijetEt_        = config->read<double>("Et max cut on dijet",100.0); 
  max3rdJetEt_       = config->read<double>("Et cut on n+1 Jet",10.0);
  maxRel3rdJetEt_    = config->read<double>("Relative n+1 Jet Et Cut",0.2);
  maxJetEta_         = config->read<double>("Eta cut on jet",5.0);
  minJetHadFraction_ = config->read<double>("Min had fraction",0.07);
  maxJetHadFraction_ = config->read<double>("Max had fraction",0.95);
  minDeltaPhi_       = config->read<double>("Min Delta Phi",2.5);
  minGenJetEt_       = config->read<double>("Et genJet min",0.0);
  maxGenJetEt_       = config->read<double>("Et genJet max",10000.0);
  maxDeltaR_         = config->read<double>("DeltaR cut on jet matching",0.25);
  // Counter for cutflow
  nDiJetCut_          = 0;
  nMinJetEt_          = 0;
  nMinDijetEt_        = 0;
  nMaxDijetEt_        = 0;
  nCutOn3rdJet_       = 0;
  nMaxJetEta_         = 0;  
  nMinDeltaPhi_       = 0;
  nMinJetHadFraction_ = 0;     
  nMaxJetHadFraction_ = 0;     
  nMinGenJetEt_       = 0;    
  nMaxGenJetEt_       = 0;     
  nMaxDeltaR_         = 0;
  // Integration parameter for SmearData
  maxNIter_          = config->read<int>("DiJet integration number of iterations",5);
  eps_               = config->read<double>("DiJet integration epsilon",1.E-5);
  std::vector<double> tBinEdges = bag_of<double>(config->read<string>("DiJet integration pt bin edges","0. 1."));
  min_               = tBinEdges.front();
  max_               = tBinEdges.back();
  // Data class
  dataClass_ = config->read<int>("Di-Jet data class", 0);
  if((dataClass_ != 0)&&(dataClass_ != 11)&&(dataClass_ != 12)&&(dataClass_ != 5)) {
    std::cout << "DiJetReader: Unknown data class " << dataClass_ << ". Using data class 0." << std::endl;
    dataClass_ = 0;
  }

  // Input files
  string default_tree_name   = config->read<string>("Default Tree Name","CalibTree");
  string treename_dijet      = config->read<string>("Di-Jet tree",default_tree_name);
  TTree * tchain_dijet;
  vector<string> input_dijet = bag_of_string(config->read<string>("Di-Jet input file","input/dijet.root"));  
  if(input_dijet[0] == "toy") { // Generate Toy MC sample
    std::cout << "generating " << nDijetEvents_ << " Di-Jet events\n";
    ToyMC* mc = new ToyMC();
    mc->init(configfile);
    mc->print();
    tchain_dijet = new TTree(treename_dijet.c_str(),"Di-Jet events");
    mc->generateDiJetTree(tchain_dijet,nDijetEvents_);
    delete mc;
  } else if(input_dijet[0] == "input/dijetlist") { // Open all files listed in "input/dijetlist"
    TChain* chain = new TChain(treename_dijet.c_str()); 
    std::ifstream filelist;
    filelist.open("input/dijetlist");
    std::string name = "";
    while( !filelist.eof() ) {
      filelist >> name;
      if( filelist.eof() ) break;
      cout << "...opening root-file " << name << " for Di-Jet analysis." << endl;
      chain->Add( name.c_str() );
    }
    filelist.close();
    tchain_dijet = chain;
  }
  else { // Open all files listed in configfile
    TChain* chain = new TChain(treename_dijet.c_str()); 
    for (bag_of_string::const_iterator it = input_dijet.begin(); it!=input_dijet.end(); ++it){
      cout << "...opening root-file " << (*it) << " for Di-Jet analysis." << endl;
      chain->Add( it->c_str() );
    }  
    tchain_dijet = chain;
  }
  nJet_.Init( tchain_dijet );
  
  delete config;
  config = 0;
}

DiJetReader::~DiJetReader()
{
  
}



//!  \brief Read dijet data and store in format as specified
//!         in the config file
//!  \param data Read data objects are appended to data
//!  \return Number of appended objects
// ----------------------------------------------------------------   
int DiJetReader::readEvents(std::vector<TData*>& data)
{
  if(nDijetEvents_ == 0) return 0;

  // Reset counters of rejected events
  nDiJetCut_          = 0;
  nMinJetEt_          = 0;
  nMinDijetEt_        = 0;
  nMaxDijetEt_        = 0;
  nCutOn3rdJet_       = 0;
  nMinGenJetEt_       = 0;    
  nMaxGenJetEt_       = 0;     
  nMaxDeltaR_         = 0;
  nMaxJetEta_         = 0;  
  nMinJetHadFraction_ = 0;     
  nMaxJetHadFraction_ = 0;     
  nMinDeltaPhi_       = 0;

  //Run jet-Jet stuff  
  int injet     = 2;
  int nevent    = nJet_.fChain->GetEntries();  // Number of events in chain
  int nReadEvts = 0;                          // Number of read events
  int nGoodEvts = 0;                          // Number of events passing all cuts

  cout << "\nReading " << injet << "-jet events...\n";
  for (int i=0;i<nevent;i++) {
    if((i+1)%10000==0) cout << i+1 << endl;
    nJet_.fChain->GetEvent(i); 
    if (nJet_.NobjTow>10000 || nJet_.NobjJet>100) {
      cerr << "ERROR: Increase array sizes in NJet_Selector; NobjTow="
	   << nJet_.NobjTow<<", NobjJet="<<nJet_.NobjJet<<"!"<<endl;
      exit(9);
    }
    if(dataClass_ == 0) {
      nReadEvts++;
      TData* td = createPtBalanceEvent(); 
      if(td) {
	nGoodEvts++;    
	data.push_back(td ); 
      } 
    } else if((dataClass_ == 11)  || (dataClass_ == 12)) {
      nReadEvts++;
      int nAddedJets = createJetTruthEvents(data);
      if( nAddedJets ) nGoodEvts += nAddedJets;    
    } else if(dataClass_ == 5) {
      nReadEvts++;
      TData* td = createSmearEvent(); 
      if(td) {
	nGoodEvts++;
	data.push_back(td ); 
      } 
    } else {
      std::cerr << "unknown data class:" << dataClass_ << '\n';
      exit(9);
    }
    if(nReadEvts>=nDijetEvents_ && nDijetEvents_>=0 ) break;
  }

  // Print cut flow
  std::cout << "Read " << nReadEvts << " " << injet << "-jet events:\n";
  if( dataClass_ == 11 || dataClass_ == 12 ) {
    std::cout << "  " << (nReadEvts-=nDiJetCut_) << std::flush;
    std::cout << " events with " << injet << " or more jets\n";
    std::cout << "  That are " << (nReadEvts*=2) << " jet-truth events:\n";
    std::cout << "    " << (nReadEvts-=nMinGenJetEt_) << std::flush;
    std::cout << " jet-truth events with ptgen > " << minGenJetEt_ << "\n";
    std::cout << "    " << (nReadEvts-=nMaxGenJetEt_) << std::flush;
    std::cout << " jet-truth events with ptgen < " << maxGenJetEt_ << "\n";
    std::cout << "    " << (nReadEvts-=nMaxDeltaR_) << std::flush;
    std::cout << " jet-truth events with DeltaR < " << maxDeltaR_ << "\n";
    std::cout << "    " << (nReadEvts-=nMinJetEt_) << std::flush;
    std::cout << " jet-truth events Et > " << minJetEt_ << "\n";
    std::cout << "    " << (nReadEvts-=nMaxJetEta_) << std::flush;
    std::cout << " jet-truth events with |eta| < " << maxJetEta_ << "\n";
    std::cout << "    " << (nReadEvts-=nMinJetHadFraction_) << std::flush;
    std::cout << " jet-truth events with hadronic fraction > " << minJetHadFraction_ << "\n";
    std::cout << "    " << (nReadEvts-=nMaxJetHadFraction_) << std::flush;
    std::cout << " jet-truth events with hadronic fraction < " << maxJetHadFraction_ << "\n";
  } else if( dataClass_ == 5 ) {
    std::cout << "  " << (nReadEvts-=nDiJetCut_) << std::flush;
    std::cout << " events with more than 3 jets\n";
    std::cout << "  " << (nReadEvts-=nMinGenJetEt_) << std::flush;
    std::cout << " dijet-jet events with ptgen > " << minGenJetEt_ << " GeV\n";
    std::cout << "  " << (nReadEvts-=nMaxGenJetEt_) << std::flush;
    std::cout << " dijet-jet events with ptgen < " << maxGenJetEt_ << " GeV\n";
    std::cout << "  " << (nReadEvts-=nMaxDeltaR_) << std::flush;
    std::cout << " dijet-jet events with DeltaR < " << maxDeltaR_ << "\n";
    std::cout << "  " << (nReadEvts-=nMinJetEt_) << std::flush;
    std::cout << " dijet-jet events Et > " << minJetEt_ << " GeV\n";
    std::cout << "  " << (nReadEvts-=nMaxJetEta_) << std::flush;
    std::cout << " dijet-jet events with |eta| < " << maxJetEta_ << "\n";
    std::cout << "  " << (nReadEvts-=nMinJetHadFraction_) << std::flush;
    std::cout << " dijet-jet events dijet pt > " << minDijetEt_ << " GeV\n";
    std::cout << "  " << (nReadEvts-=nMinDijetEt_) << std::flush;
    std::cout << " dijet-jet events dijet pt < " << maxDijetEt_ << " GeV\n";
    std::cout << "  " << (nReadEvts-=nMaxDijetEt_) << std::flush;
    std::cout << " dijet-jet events with hadronic fraction > " << minJetHadFraction_ << "\n";
    std::cout << "  " << (nReadEvts-=nMaxJetHadFraction_) << std::flush;
    std::cout << " dijet-jet events with hadronic fraction < " << maxJetHadFraction_ << "\n";
    std::cout << "  " << (nReadEvts-=nCutOn3rdJet_) << std::flush;
    std::cout << " dijet-jet events with pt(jet3) / pt(dijet) < " << maxRel3rdJetEt_ << " or ";
    std::cout << "pt(jet3) < " << max3rdJetEt_ << " GeV\n";
    std::cout << "  " << (nReadEvts-=nMinDeltaPhi_) << std::flush;
    std::cout << " dijet-jet events with DeltaPhi > " << minDeltaPhi_ << "\n";
  }
  std::cout << "Stored " << nGoodEvts << " events for analysis.\n";
  return nGoodEvts;
}
  



//!  \brief Create TData_PtBalance event from dijet data
//!  \return Pointer to TData_PtBalance event (0 if cuts are not passed)
// ----------------------------------------------------------------   
TData* DiJetReader::createPtBalanceEvent()
{
  //--------------
  //  n - Jet
  //-------------- 
  int injet = 2;
  TData_PtBalance * jj_data[nJet_.NobjJet];
  jj_data[0] = 0;
  //std::cout << "reading " << nJet_.NobjJet << " jets\n";
  
  int nstoredjets = 0;
  for (unsigned int ij = 0; (int)ij<nJet_.NobjJet; ++ij){
    if(nJet_.JetPt[ij] < max3rdJetEt_) continue;
    //Find the jets eta & phi index using the nearest tower to jet axis:
    int jet_index=-1;
    double min_tower_dr = 10.0;
    double em = 0;
    double had = 0;
    double out = 0;
    TLorentzVector LJet(0,0,0,0);
    LJet.SetPtEtaPhiE(nJet_.JetPt[ij],nJet_.JetEta[ij],nJet_.JetPhi[ij],nJet_.JetE[ij]);
    for (int n=0; n<nJet_.NobjTow; ++n){
      if (nJet_.Tow_jetidx[n]!=(int)ij) continue;//look for ij-jet's towers
      em += nJet_.TowEm[n];
      had += nJet_.TowHad[n];
      out += nJet_.TowOE[n];
      TLorentzVector Ltower(0,0,0,0);
      Ltower.SetPtEtaPhiE(nJet_.TowEt[n],nJet_.TowEta[n],nJet_.TowPhi[n],nJet_.TowE[n]);
      double dr = Ltower.DeltaR(LJet);
      if (dr<min_tower_dr) {
	jet_index = p->GetJetBin(p->GetJetEtaBin(nJet_.TowId_eta[n]),
				 p->GetJetPhiBin(nJet_.TowId_phi[n]));
	min_tower_dr = dr;
      }
    }
    if (jet_index<0){ 
      cerr<<"WARNING: JJ jet_index = " << jet_index << endl; 
      continue; 
    }
    double * direction = new double[2];
    direction[0] = sin(nJet_.JetPhi[ij]);
    direction[1] = cos(nJet_.JetPhi[ij]);
    TJet* jetp  = new TJet;
    jetp->pt  = nJet_.JetEt[ij];
    jetp->eta = nJet_.JetEta[ij];
    jetp->phi = nJet_.JetPhi[ij];
    jetp->E   = nJet_.JetE[ij];
    jetp->genPt =nJet_.GenJetPt[ij];
    jetp->corFactors = TJet::CorFactors(nJet_.JetCorrZSP[ij], // L1
					nJet_.JetCorrL2[ij],  // L2
					nJet_.JetCorrL3[ij],  // L3
					1.,              // L4
					1.,              // L5
					nJet_.JetCorrJPT[ij],
					nJet_.JetCorrL2L3JPT[ij]);
    //the following is not quite correct, as this factor is different for all towers. These values should be in the n-tupel as well
      double factor =  nJet_.JetEt[ij] /  nJet_.JetE[ij];
      jetp->HadF = had * factor;
      jetp->EMF = em * factor;
      jetp->OutF = out * factor;
      //Create an jet/Jet TData event
      jj_data[nstoredjets] = new TData_PtBalance( 
          jet_index * p->GetNumberOfJetParametersPerBin() + p->GetNumberOfTowerParameters(),
	  direction,                                     //p_T direction of this jet
	  0.0,                                           //truth//
	  sqrt(pow(0.5,2)+pow(0.10*nJet_.JetPt[ij],2)),   //error//
	  //nJet_.Weight,                                   //weight//
	  1.,                                          //weight//
	  p->GetJetParRef( jet_index ),                  //params
	  p->GetNumberOfJetParametersPerBin(),           //number of free jet param. p. bin
	  p->jet_parametrization,                        //function
	  //p->dummy_parametrization,
          jet_error_param,                               //error param. function
	  jetp                                           //jet momentum for plotting and scale
        );
//cout << "jet "<<nstoredjets<<"'s E="<<nJet_.JetE[ij]
//     << ", ntower:"<<endl;
      //Add the jet's towers to "jj_data":
      for (int n=0; n<nJet_.NobjTow; ++n){
        if (nJet_.Tow_jetidx[n]!=(int)ij) continue;//look for ij-jet's towers
	//if (nJet_.TowEt[n]<0.01) continue;

	int index = p->GetBin(p->GetEtaBin(nJet_.TowId_eta[n]),
			      p->GetPhiBin(nJet_.TowId_phi[n]));
//std::cout << "jet:" << ij << ", towid=" << n << ", bin index:" << index << "\n";
	if (index<0){ cerr<<"WARNING: JJ tower_index = " << index << endl; continue; }

	double relativEt = nJet_.TowEt[n]/nJet_.JetEt[ij];  
	//if (relativEt<=0) cerr << "relEt = " <<relativEt << endl; //continue;
	//This relativeE is used *only* for plotting! Therefore no cuts on this var!
	//create array with multidimensional measurement
	TMeasurement * mess = new TTower;
	mess->pt = double(nJet_.TowEt[n]);
	double scale = nJet_.TowEt[n]/nJet_.TowE[n];
	mess->EMF = double(nJet_.TowEm[n]*scale);
	mess->HadF = double(nJet_.TowHad[n]*scale);
	mess->OutF = double(nJet_.TowOE[n]*scale);
	mess->eta = double(nJet_.TowEta[n]);
	mess->phi = double(nJet_.TowPhi[n]);
	mess->E = double(nJet_.TowE[n]);
	//mess[7] = double( cos( nJet_.JetCalPhi-nJet_.TowPhi[n] ) ); // Projection factor for summing tower Pt

	jj_data[nstoredjets]->AddMess(new TData_TruthMess(
	    index,
	    mess,                                                   //mess//
	    nJet_.JetPt[ij] * relativEt,                             //truth//
	    sqrt(pow(0.5,2)+pow(0.1*nJet_.JetPt[ij]*relativEt,2)),   //error//
            //1.,                                                   //weight//
	    nJet_.Weight,                                            //weight//
	    p->GetTowerParRef( index ),                             //parameter//
	    p->GetNumberOfTowerParametersPerBin(),                  //number of free tower param. p. bin//
	    p->tower_parametrization,                               //function//
	    tower_error_param                                      //error param. function//
	  ));
      }
      //Add the jet's tracks to "gj_data":
      double* EfficiencyMap = p->GetEffMap();
      int track_index;
      for (int n=0; n<nJet_.NobjTrack; ++n){
        if (nJet_.Track_jetidx[n]!=(int)ij) continue;//look for ij-jet's tracks

	if((nJet_.TrackTowIdEta[n] == 0) || (nJet_.TrackTowIdPhi[n] == 0)) {
	  if(nJet_.TrackPt[n] > 2){
	    std::cerr << "WARNING: eta or phi id of track is zero!\n";
	    continue;
	  }
	  else track_index = 0; //bent low momentum tracks with no HCAL hit
	}
	else
	  track_index = p->GetTrackBin(p->GetTrackEtaBin(nJet_.TrackTowIdEta[n]),
					   p->GetTrackPhiBin(nJet_.TrackTowIdPhi[n]));
	if (track_index<0){ cerr<<"WARNING: JJ track_index = " << track_index << endl; continue; }
	//create array with multidimensional measurement
	//TMeasurement * Tmess = new TTrack;
	TTrack * Tmess = new TTrack;
	Tmess->TrackId = int(nJet_.TrackId[n]);
	Tmess->TowerId = int(nJet_.TrackTowId[n]);
	Tmess->pt = double(nJet_.TrackPt[n]);
	double scale = nJet_.TrackP[n]/nJet_.TrackPt[n];
	Tmess->EM1 = double(nJet_.TrackEMC1[n]*scale);
	Tmess->EMF = double(nJet_.TrackEMC3[n]*scale);
	Tmess->EM5 = double(nJet_.TrackEMC5[n]*scale);
	Tmess->Had1 = double(nJet_.TrackHAC1[n]*scale);
	Tmess->HadF = double(nJet_.TrackHAC3[n]*scale);
	Tmess->Had5 = double(nJet_.TrackHAC5[n]*scale);
	Tmess->OutF = 0;
	Tmess->DR = double(nJet_.TrackDR[n]);
	Tmess->DRout = double(nJet_.TrackDROut[n]);
	Tmess->eta = double(nJet_.TrackEta[n]);
	Tmess->etaOut = double(nJet_.TrackEtaOut[n]);
	Tmess->phi = double(nJet_.TrackPhi[n]);
	Tmess->phiOut = double(nJet_.TrackPhiOut[n]);
	Tmess->E = double(nJet_.TrackP[n]);
	Tmess->TrackChi2 = double(nJet_.TrackChi2[n]);
	Tmess->NValidHits = int(nJet_.TrackNHits[n]);
	Tmess->TrackQualityT = bool(nJet_.TrackQualityT[n]);
	Tmess->MuDR = double(nJet_.MuDR[n]);
	Tmess->MuDE = double(nJet_.MuDE[n]);
	int TrackEffBin = p->GetTrackEffBin(nJet_.TrackPt[n],nJet_.TrackEta[n]);
	Tmess->Efficiency = EfficiencyMap[TrackEffBin];
	//mess[7] = double( cos( nJet_.JetCalPhi-nJet_.TowPhi[n] ) ); // Projection factor for summing tower Pt
	//EM+=mess->EMF;
	//F+=mess->pt;
	jj_data[nstoredjets]->AddTrack(new TData_TruthMess(
					      track_index  * p->GetNumberOfTrackParametersPerBin() + p->GetNumberOfTowerParameters() + p->GetNumberOfJetParameters() ,
					      Tmess,                                                    //mess//
					      0,                           //truth//
					      0.05 + 0.00015 * nJet_.TrackPt[n], //error//
					      1.,                                                      //weight//
					      p->GetTrackParRef( track_index ),                              //parameter//
					      p->GetNumberOfTrackParametersPerBin(),                   //number of free tower param. p. bin//
					      p->track_parametrization,                                //function//
					      track_error_param                                        //error param.func.//
					      ));
      }
      jj_data[nstoredjets]->UseTracks(useTracks);   //check if track information is sufficient to use Track Parametrization  

    if(nstoredjets> 0)  
      jj_data[0]->AddNewMultMess( jj_data[nstoredjets] );
    ++nstoredjets;
  }//loop over all n-jets
  bool goodevent=true;
  if (nstoredjets < injet) goodevent = false;
  if (nstoredjets > injet){
    /*
      for (int i=0; i < nstoredjets; ++i){
      cout<<i<<"-ter Jet Pt: "<<jj_data[i]->GetMess()->pt<<endl;
      }
    */
    //relative Pt cut only works if jets are Pt sorted
    double scale=0;
    for (int i=0; i < injet; ++i){
      scale += jj_data[i]->GetMess()->pt;
    }
    scale /= injet;
    //cout<<"scale: "<<scale<<endl;
    if ( jj_data[injet]->GetMess()->pt > scale*maxRel3rdJetEt_ ) goodevent = false;
  }
  if ( (((TJet*)(jj_data[0]->GetMess()))->genPt < minGenJetEt_) || (((TJet*)(jj_data[1]->GetMess()))->genPt  < minGenJetEt_ ))  goodevent = false;
  if ( (std::abs(jj_data[0]->GetMess()->eta) > maxJetEta_) ||  (std::abs(jj_data[1]->GetMess()->eta) > maxJetEta_) )  goodevent = false;
  
  /*
  //sort jets. 1st is barrel, 2nd is probe
  if( nstoredjets ==  2) {
  if(std::abs(jj_data[0]->GetMess()[1]) > 1.2) {
  if(std::abs(jj_data[1]->GetMess()[1]) > 1.2) {
  delete jj_data[0];
  continue;
  } else {
  jj_data[0]->ClearMultMess();
  jj_data[1]->AddNewMultMess(jj_data[0]);
  TData_PtBalance* tmp = jj_data[1];
  jj_data[1] = jj_data[0];
  jj_data[0] = tmp; 
  }
  } else if(std::abs(jj_data[1]->GetMess()[1]) < 1.2) {
  //both jets central, roll the dice and swap
  if(rand()/(RAND_MAX+1.0) > 0.5) {
  jj_data[0]->ClearMultMess();
  jj_data[1]->AddNewMultMess(jj_data[0]);
  TData_PtBalance* tmp = jj_data[1];
  jj_data[1] = jj_data[0];
  jj_data[0] = tmp; 
  }
  }
  }    
  */
  if(! goodevent) {
    delete jj_data[0];
    return 0;
  }
  return jj_data[0];
}



//!  \brief Use first two jets of dijet event as two JetTruthEvent
//!         objects where genjet Et is truth
//!  \note The dijets are ordered in genjet Et
//!  \param data Read JetTruthEvent objects are appended to data
//!  \return Number of appended JetTruthEvent objects (0 - 2)
// ----------------------------------------------------------------   
int DiJetReader::createJetTruthEvents(std::vector<TData*>& data)
{
  int injet = 2;
  if( nJet_.NobjGenJet < injet ) {
    nDiJetCut_++;
    return 0;
  }

  int     njets = 0;   // Number of stored JetTruthEvents; 0 or 2
  double * terr = new double[nJet_.NobjTow];


  // Loop over two jets with highest genjet Et
  for(int genJetIdx = 0; genJetIdx < 2; genJetIdx++) {
    int calJetIdx = nJet_.GenJetColJetIdx[genJetIdx]; // Closest (DeltaR) calo jet
    if( nJet_.NobjJet <= calJetIdx ) {
      nDiJetCut_++;
      return 0;
    }

    // Cuts
    if( nJet_.GenJetColEt[genJetIdx] < minGenJetEt_ ) {
      nMinGenJetEt_++;
      continue;
    } else if( nJet_.GenJetColEt[genJetIdx] > maxGenJetEt_ ) {
      nMaxGenJetEt_++;
      continue;
    }
    double dphi        = TVector2::Phi_mpi_pi( nJet_.JetPhi[calJetIdx] - nJet_.GenJetColPhi[genJetIdx] );
    double deta        = nJet_.JetEta[calJetIdx] - nJet_.GenJetColEta[genJetIdx];
    double drJetGenjet = sqrt( deta*deta + dphi*dphi );
    if( drJetGenjet > maxDeltaR_ ) {
      nMaxDeltaR_++;
      continue;
    } else if( nJet_.JetEt[calJetIdx] < minJetEt_ ) {
      nMinJetEt_++;
      continue;
    } else if( std::abs(nJet_.JetEta[calJetIdx]) > maxJetEta_ ) {
      nMaxJetEta_++;
      continue;
    }

    // Construct event
    double em   = 0;
    double had  = 0;
    double out  = 0;
    double err2 = 0;
    TMeasurement tower;
    double dR        = 10;
    int closestTower = 0; 
    for(int n=0; n<nJet_.NobjTow; ++n){
      if(nJet_.Tow_jetidx[n] != calJetIdx) continue;//look for ij-jet's towers

      em          += nJet_.TowEm[n];
      had         += nJet_.TowHad[n];
      out         += nJet_.TowOE[n];  
      tower.pt     = nJet_.TowEt[n];
      double scale = nJet_.TowEt[n]/nJet_.TowE[n];
      tower.EMF    = nJet_.TowEm[n]*scale;
      tower.HadF   = nJet_.TowHad[n]*scale;
      tower.OutF   = nJet_.TowOE[n]*scale;
      tower.eta    = nJet_.TowEta[n];
      tower.phi    = nJet_.TowPhi[n];
      tower.E      = nJet_.TowE[n];
      terr[n]      = tower_error_param(&tower.pt,&tower,0); 
      if(terr[n] == 0) {
	//assume toy MC???
	terr[n] = TParameters::toy_tower_error_parametrization(&tower.pt,&tower);
      }
      terr[n]  *= terr[n];
      err2     += terr[n];
      dphi      = TVector2::Phi_mpi_pi(nJet_.JetPhi[calJetIdx]-tower.phi);
      deta      = nJet_.JetEta[calJetIdx]-tower.eta;
      double dr = sqrt( deta*deta + dphi*dphi );     
      if(dr < dR) {
	dR = dr;
	closestTower = n;
      }
    }
    // Cuts on hadronic fraction 
    if(had/(had + em) < minJetHadFraction_) {
      nMinJetHadFraction_++;
      continue;
    } else if(had/(had + em) > maxJetHadFraction_) { 
      nMaxJetHadFraction_++;
      continue;
    }
    double factor = nJet_.JetEt[calJetIdx] /  nJet_.JetE[calJetIdx];
    tower.pt      = nJet_.JetEt[calJetIdx];
    tower.EMF     = em * factor;
    tower.HadF    = had * factor;
    tower.OutF    = out * factor;
    tower.eta     = nJet_.JetEta[calJetIdx];
    tower.phi     = nJet_.JetPhi[calJetIdx];
    tower.E       = nJet_.JetE[calJetIdx];
    double err    = jet_error_param(&tower.pt,&tower,0);
    err2         += err * err;

    Jet *jet;
    if(dataClass_ == 12) {
      JetWithTowers *jt = 
	new JetWithTowers(nJet_.JetEt[calJetIdx],em * factor,had * factor,
			  out * factor,nJet_.JetE[calJetIdx],nJet_.JetEta[calJetIdx],
			  nJet_.JetPhi[calJetIdx],TJet::uds,nJet_.GenJetColEt[genJetIdx],drJetGenjet,
			  TJet::CorFactors(nJet_.JetCorrZSP[calJetIdx], // L1
					   nJet_.JetCorrL2[calJetIdx],  // L2
					   nJet_.JetCorrL3[calJetIdx],  // L3
					   1.,                         // L4
					   1.,                         // L5
					   nJet_.JetCorrJPT[calJetIdx],
					   nJet_.JetCorrL2[calJetIdx]*nJet_.JetCorrL3[calJetIdx]), //not the JPT specific L2L3 factors?
			  p->jet_function(nJet_.TowId_eta[closestTower],
					  nJet_.TowId_phi[closestTower]),
			  jet_error_param,p->global_jet_function(),minJetEt_);
      for(int j = 0 ; j < nJet_.NobjTow ; ++j) {
	if (nJet_.Tow_jetidx[j]!= calJetIdx) continue;//look for ij-jet's towers
	double scale = nJet_.TowEt[j]/nJet_.TowE[j];
	jt->addTower(nJet_.TowEt[j],nJet_.TowEm[j]*scale,
		     nJet_.TowHad[j]*scale,nJet_.TowOE[j]*scale,
		     nJet_.TowE[j],nJet_.TowEta[j],nJet_.TowPhi[j],
		     p->tower_function(nJet_.TowId_eta[calJetIdx],nJet_.TowId_phi[calJetIdx]),
		     tower_error_param);
      }
      jet = jt;
    }
    else { 
      jet = new Jet(nJet_.JetEt[calJetIdx],em * factor,had * factor,out * factor,
		    nJet_.JetE[calJetIdx],nJet_.JetEta[calJetIdx],nJet_.JetPhi[calJetIdx],
		    TJet::uds,nJet_.GenJetColEt[genJetIdx],drJetGenjet,
		    TJet::CorFactors(nJet_.JetCorrZSP[calJetIdx], // L1
				     nJet_.JetCorrL2[calJetIdx],  // L2
				     nJet_.JetCorrL3[calJetIdx],  // L3
				     1.,                         // L4
				     1.,                         // L5
				     nJet_.JetCorrJPT[calJetIdx],
				     nJet_.JetCorrL2[calJetIdx]*nJet_.JetCorrL3[calJetIdx]), //not the JPT specific L2L3 factors?
		    p->jet_function(nJet_.TowId_eta[closestTower],
				    nJet_.TowId_phi[closestTower]),
		    jet_error_param,p->global_jet_function(),minJetEt_);    
    }

    JetTruthEvent* jte = new JetTruthEvent(jet,nJet_.GenJetColEt[genJetIdx],1.);//nJet_.Weight);
    data.push_back(jte);
    ++njets;
  }     
  delete [] terr;
  return njets;
}



//!  \brief Create \p SmearDiJet \p event for jet smearing
//!
//!  Uses the three jets with the highest uncorrected
//!  calo pt. For creation of a \p SmearDiJet \p event,
//!  the JetMET L2L3 correction is applied.
//!
//!  \return A \p SmearDiJet \p if all cuts are passed,
//!          else 0
// ----------------------------------------------------------------   
TData* DiJetReader::createSmearEvent()
{
  if( nJet_.NobjJet < 3 ) {
    nDiJetCut_++;
    return 0;
  }

  SmearDiJet * jj_data = 0;
  
  TJet   * jet1        = 0;
  TJet   * jet2        = 0;

  // Loop over three jets with highest calojet Et
  for(int jetIdx = 0; jetIdx < 3; jetIdx++) {

    double dphi         = TVector2::Phi_mpi_pi( nJet_.JetPhi[jetIdx] - nJet_.GenJetPhi[jetIdx] );
    double deta         = nJet_.JetEta[jetIdx] - nJet_.GenJetEta[jetIdx];
    double drJetGenjet  = sqrt( deta*deta + dphi*dphi );
    double min_tower_dr = 10.;
    double emf          = 0;
    double had          = 0;
    double out          = 0;
    int    closestTower = 0; 

    // Loop over towers
    for (int n=0; n<nJet_.NobjTow; ++n) {
      if (nJet_.Tow_jetidx[n]!=(int)jetIdx) continue;//look for jetIdx-jet's towers
      emf += nJet_.TowEm[n];
      had += nJet_.TowHad[n];
      out += nJet_.TowOE[n];
      dphi = TVector2::Phi_mpi_pi( nJet_.JetPhi[jetIdx] - nJet_.TowPhi[n] );
      deta = nJet_.JetEta[jetIdx] - nJet_.TowEta[n];
      double dr = sqrt( deta*deta + dphi*dphi );     
      if (dr < min_tower_dr) {
	min_tower_dr = dr;
	closestTower = n;
      }
    } // End of loop over towers


    // Projection factor E --> Et
    // The following is not quite correct, as this factor is different for all towers
    // These values should be in the n-tupel as well
    double projFac   = nJet_.JetEt[jetIdx] /  nJet_.JetE[jetIdx];

    // Want L2L3 corrected jets
    double corrFac   = nJet_.JetCorrL2[jetIdx] * nJet_.JetCorrL3[jetIdx];

    // Set up measurement
    TJet * jetp      = new TJet;
    jetp->pt         = corrFac * nJet_.JetPt[jetIdx];
    jetp->E          = corrFac * nJet_.JetE[jetIdx];
    jetp->HadF       = corrFac * had * projFac;
    jetp->EMF        = corrFac * emf * projFac;
    jetp->OutF       = corrFac * out * projFac;
    jetp->eta        = nJet_.JetEta[jetIdx];
    jetp->phi        = nJet_.JetPhi[jetIdx];
    jetp->genPt      = nJet_.GenJetPt[jetIdx];
    jetp->dR         = drJetGenjet;
    jetp->ptHat      = nJet_.GenEvtScale;

    // All corrections initialised with 1. as jets are already L2L3 corrected
    jetp->corFactors = TJet::CorFactors();


    if     ( jetIdx == 0 ) { // Store first jet
      jet1 = jetp;
    }
    else if( jetIdx == 1 ) { // Store second jet
      jet2 = jetp;
    }
    else if( jetIdx == 2) { // Create a SmearDiJet event
      int etaBin = nJet_.TowId_eta[closestTower]; // This is the 3rd jet, needs to be adjusted!
      int phiBin = nJet_.TowId_phi[closestTower]; // This is the 3rd jet, needs to be adjusted!
      jj_data = new SmearDiJet(jet1,                       // First jet
			       jet2,                       // Second jet
			       jetp,                       // Third jet
			       1.,                         // Weights from EventProcessor
			       p->jet_function(etaBin,
					       phiBin),     // Response pdf
			       p->global_jet_function(),   // Truth pdf
			       min_,                       // Integration minimum
			       max_,                       // Integration maximum
			       eps_,                       // Integration step length
			       maxNIter_);                 // Integration n iterations
    }
  }  // End of loop over jets


  // Check if event is ok and return
  bool isGoodEvt = true;
  TJet * j1 = static_cast<TJet*>(jj_data->GetMess());
  TJet * j2 = static_cast<TJet*>(jj_data->GetSecondMess());
  TJet * j3 = static_cast<TJet*>(jj_data->GetThirdMess());
  if     ( j1->genPt < minGenJetEt_ || j2->genPt < minGenJetEt_ ) {
    nMinGenJetEt_++;
    isGoodEvt = false;
  }
  else if( j1->genPt > maxGenJetEt_ || j2->genPt > maxGenJetEt_ ) {
    nMaxGenJetEt_++;
    isGoodEvt = false;
  }
  else if( j1->dR > maxDeltaR_ || j2->dR > maxDeltaR_ ) {
    nMaxDeltaR_++;
    isGoodEvt = false;
  }
  else if( j1->pt < minJetEt_ || j2->pt < minJetEt_ ) {
    nMinJetEt_++;
    isGoodEvt = false;
  }
  else if( jj_data->dijetPt() < minDijetEt_ ) {
    nMinDijetEt_++;
    isGoodEvt = false;
  }
  else if( jj_data->dijetPt() > maxDijetEt_ ) {
    nMaxDijetEt_++;
    isGoodEvt = false;
  }
  else if( std::abs(j1->eta) > maxJetEta_ || std::abs(j2->eta) > maxJetEta_ ) {
    nMaxJetEta_++;
    isGoodEvt = false;
  }
  else if( j1->HadF/(j1->HadF + j1->EMF) < minJetHadFraction_ ||
	   j2->HadF/(j2->HadF + j2->EMF) < minJetHadFraction_ ) {
    nMinJetHadFraction_++;
    isGoodEvt = false;
  }
  else if( j1->HadF/(j1->HadF + j1->EMF) > maxJetHadFraction_ ||
	   j2->HadF/(j2->HadF + j2->EMF) > maxJetHadFraction_ ) {
    nMaxJetHadFraction_++;
    isGoodEvt = false;
  }
  else if( j3->pt / jj_data->dijetPt() > maxRel3rdJetEt_ && j3->pt > max3rdJetEt_ ) {
    nCutOn3rdJet_++;
    isGoodEvt = false;
  }
  else if( std::abs(TVector2::Phi_mpi_pi(j1->phi - j2->phi)) < minDeltaPhi_ ) {
    nMinDeltaPhi_++;
    isGoodEvt = false;
  }

  if(! isGoodEvt) {
    if( jj_data ) delete jj_data;
    jj_data = 0;
  }
  return jj_data;
}
