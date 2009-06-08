//
//    first version: Hartmut Stadie 2008/12/12
//    $Id: DiJetReader.cc,v 1.12 2009/06/05 15:44:20 mschrode Exp $
//   
#include "DiJetReader.h"

#include "CalibData.h"
#include "ConfigFile.h"
#include "ToyMC.h"
#include "Parameters.h"
#include "JetTruthEvent.h"
#include "Jet.h"
#include "JetWithTowers.h"

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
  n_dijet_events = config->read<int>("use Di-Jet events",-1);
  if(n_dijet_events == 0) {
    delete config;
    config = 0;
    return;
  }

  // Cuts
  Et_cut_on_jet     = config->read<double>("Et cut on jet",0.0); 
  Et_cut_nplus1Jet  = config->read<double>("Et cut on n+1 Jet",10.0);
  Rel_cut_on_nJet   = config->read<double>("Relative n+1 Jet Et Cut",0.2);
  GenJetCutLow      = config->read<double>("Et genJet min",0.0);
  GenJetCutUp       = config->read<double>("Et genJet max",10000.0);
  DeltaRMatchingCut = config->read<double>("DeltaR cut on jet matching",0.25);
  Eta_cut_on_jet    = config->read<double>("Eta cut on jet",5.0);
  Had_cut_min       = config->read<double>("Min had fraction",0.07);
  Had_cut_max       = config->read<double>("Max had fraction",0.95);
  // Counter for cutflow
  nEt_cut_on_jet     = 0;
  nNjet_cut          = 0;
  nEt_cut_nplus1Jet  = 0;
  nRel_cut_on_nJet   = 0; 
  nGenJetCutLow      = 0;    
  nGenJetCutUp       = 0;     
  nDeltaRMatchingCut = 0;
  nEta_cut_on_jet    = 0;  
  nHad_cut_min       = 0;     
  nHad_cut_max       = 0;     

  // Data class
  dataClass = config->read<int>("Di-Jet data class", 0);
  if((dataClass != 0)&&(dataClass != 11)&&(dataClass != 12)) {
    std::cout << "DiJetReader: Unknown data class " << dataClass << ". Using data class 0." << std::endl;
    dataClass = 0;
  }

  // Input files
  string default_tree_name   = config->read<string>("Default Tree Name","CalibTree");
  string treename_dijet      = config->read<string>("Di-Jet tree",default_tree_name);
  TTree * tchain_dijet;
  vector<string> input_dijet = bag_of_string(config->read<string>("Di-Jet input file","input/dijet.root"));  
  if(input_dijet[0] == "toy") { // Generate Toy MC sample
    std::cout << "generating " << n_dijet_events << " Di-Jet events\n";
    ToyMC* mc = new ToyMC();
    mc->init(configfile);
    mc->print();
    tchain_dijet = new TTree(treename_dijet.c_str(),"Di-Jet events");
    mc->generateDiJetTree(tchain_dijet,n_dijet_events);
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
  njet.Init( tchain_dijet );
  
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
  if(n_dijet_events == 0) return 0;

  // Reset counters of rejected events
  nNjet_cut          = 0;
  nEt_cut_on_jet     = 0;
  nEt_cut_nplus1Jet  = 0;
  nRel_cut_on_nJet   = 0; 
  nGenJetCutLow      = 0;    
  nGenJetCutUp       = 0;     
  nDeltaRMatchingCut = 0;
  nEta_cut_on_jet    = 0;  
  nHad_cut_min       = 0;     
  nHad_cut_max       = 0;     

  //Run jet-Jet stuff  
  int injet     = 2;
  int nevent    = njet.fChain->GetEntries();  // Number of events in chain
  int nReadEvts = 0;                          // Number of read events
  int nGoodEvts = 0;                          // Number of events passing all cuts

  cout << "\nReading " << injet << "-jet events...\n";
  for (int i=0;i<nevent;i++) {
    if((i+1)%10000==0) cout << i+1 << endl;
    njet.fChain->GetEvent(i); 
    if (njet.NobjTow>10000 || njet.NobjJet>100) {
      cerr << "ERROR: Increase array sizes in NJetSelector; NobjTow="
	   << njet.NobjTow<<", NobjJet="<<njet.NobjJet<<"!"<<endl;
      exit(9);
    }
    if(dataClass == 0) {
      nReadEvts++;
      TData* td = createPtBalanceEvent(); 
      if(td) {
	nGoodEvts++;    
	data.push_back(td ); 
      } 
    } else if((dataClass == 11)  || (dataClass == 12)) {
      nReadEvts++;
      int nAddedJets = createJetTruthEvents(data);
      if( nAddedJets ) nGoodEvts += nAddedJets;    
    } else {
      std::cerr << "unknown data class:" << dataClass << '\n';
      exit(9);
    }
    if(nReadEvts>=n_dijet_events && n_dijet_events>=0 ) break;
  }

  // Print cut flow
  std::cout << "Read " << nReadEvts << " " << injet << "-jet events:\n";
  if( dataClass == 11 || dataClass == 12 ) {
    std::cout << "  " << (nReadEvts-=nNjet_cut) << std::flush;
    std::cout << " events with more than " << injet << " jets\n";
    std::cout << "  That are " << (nReadEvts*=2) << " jet-truth events:\n";
    std::cout << "    " << (nReadEvts-=nGenJetCutLow) << std::flush;
    std::cout << " jet-truth events with ptgen > " << GenJetCutLow << "\n";
    std::cout << "    " << (nReadEvts-=nGenJetCutUp) << std::flush;
    std::cout << " jet-truth events with ptgen < " << GenJetCutUp << "\n";
    std::cout << "    " << (nReadEvts-=nDeltaRMatchingCut) << std::flush;
    std::cout << " jet-truth events with DeltaR < " << DeltaRMatchingCut << "\n";
    std::cout << "    " << (nReadEvts-=nEt_cut_on_jet) << std::flush;
    std::cout << " jet-truth events Et > " << Et_cut_on_jet << "\n";
    std::cout << "    " << (nReadEvts-=nEta_cut_on_jet) << std::flush;
    std::cout << " jet-truth events with |eta| < " << Eta_cut_on_jet << "\n";
    std::cout << "    " << (nReadEvts-=nHad_cut_min) << std::flush;
    std::cout << " jet-truth events with hadronic fraction > " << Had_cut_min << "\n";
    std::cout << "    " << (nReadEvts-=nHad_cut_max) << std::flush;
    std::cout << " jet-truth events with hadronic fraction < " << Had_cut_max << "\n";
  }
  std::cout << "Stored " << nGoodEvts << " jet-truth events for calibration.\n";
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
  TData_PtBalance * jj_data[njet.NobjJet];
  jj_data[0] = 0;
  //std::cout << "reading " << njet.NobjJet << " jets\n";
  
  int nstoredjets = 0;
  for (unsigned int ij = 0; (int)ij<njet.NobjJet; ++ij){
    if(njet.JetPt[ij] < Et_cut_nplus1Jet) continue;
    //Find the jets eta & phi index using the nearest tower to jet axis:
    int jet_index=-1;
    double min_tower_dr = 10.0;
    double em = 0;
    double had = 0;
    double out = 0;
    TLorentzVector Ljet(0,0,0,0);
    Ljet.SetPtEtaPhiE(njet.JetPt[ij],njet.JetEta[ij],njet.JetPhi[ij],njet.JetE[ij]);
    for (int n=0; n<njet.NobjTow; ++n){
      if (njet.Tow_jetidx[n]!=(int)ij) continue;//look for ij-jet's towers
      em += njet.TowEm[n];
      had += njet.TowHad[n];
      out += njet.TowOE[n];
      TLorentzVector Ltower(0,0,0,0);
      Ltower.SetPtEtaPhiE(njet.TowEt[n],njet.TowEta[n],njet.TowPhi[n],njet.TowE[n]);
      double dr = Ltower.DeltaR(Ljet);
      if (dr<min_tower_dr) {
	jet_index = p->GetJetBin(p->GetJetEtaBin(njet.TowId_eta[n]),
				 p->GetJetPhiBin(njet.TowId_phi[n]));
	min_tower_dr = dr;
      }
    }
    if (jet_index<0){ 
      cerr<<"WARNING: JJ jet_index = " << jet_index << endl; 
      continue; 
    }
    double * direction = new double[2];
    direction[0] = sin(njet.JetPhi[ij]);
    direction[1] = cos(njet.JetPhi[ij]);
    TJet* jetp  = new TJet;
    jetp->pt  = njet.JetEt[ij];
    jetp->eta = njet.JetEta[ij];
    jetp->phi = njet.JetPhi[ij];
    jetp->E   = njet.JetE[ij];
    jetp->genPt =njet.GenJetPt[ij];
    jetp->ZSPcor =njet.JetCorrZSP[ij]; 
    jetp->JPTcor =njet.JetCorrJPT[ij]; 
    jetp->L2cor =njet.JetCorrL2[ij]; 
    jetp->L3cor =njet.JetCorrL3[ij]; 
    jetp->L2L3cor =njet.JetCorrL2L3[ij]; 
    jetp->L2L3JPTcor =njet.JetCorrL2L3JPT[ij]; 
    //the following is not quite correct, as this factor is different for all towers. These values should be in the n-tupel as well
      double factor =  njet.JetEt[ij] /  njet.JetE[ij];
      jetp->HadF = had * factor;
      jetp->EMF = em * factor;
      jetp->OutF = out * factor;
      //Create an jet/Jet TData event
      jj_data[nstoredjets] = new TData_PtBalance( 
          jet_index * p->GetNumberOfJetParametersPerBin() + p->GetNumberOfTowerParameters(),
	  direction,                                     //p_T direction of this jet
	  0.0,                                           //truth//
	  sqrt(pow(0.5,2)+pow(0.10*njet.JetPt[ij],2)),   //error//
	  //njet.Weight,                                   //weight//
	  1.,                                          //weight//
	  p->GetJetParRef( jet_index ),                  //params
	  p->GetNumberOfJetParametersPerBin(),           //number of free jet param. p. bin
	  p->jet_parametrization,                        //function
	  //p->dummy_parametrization,
          jet_error_param,                               //error param. function
	  jetp                                           //jet momentum for plotting and scale
        );
//cout << "jet "<<nstoredjets<<"'s E="<<njet.JetE[ij]
//     << ", ntower:"<<endl;
      //Add the jet's towers to "jj_data":
      for (int n=0; n<njet.NobjTow; ++n){
        if (njet.Tow_jetidx[n]!=(int)ij) continue;//look for ij-jet's towers
	//if (njet.TowEt[n]<0.01) continue;

	int index = p->GetBin(p->GetEtaBin(njet.TowId_eta[n]),
			      p->GetPhiBin(njet.TowId_phi[n]));
//std::cout << "jet:" << ij << ", towid=" << n << ", bin index:" << index << "\n";
	if (index<0){ cerr<<"WARNING: JJ tower_index = " << index << endl; continue; }

	double relativEt = njet.TowEt[n]/njet.JetEt[ij];  
	//if (relativEt<=0) cerr << "relEt = " <<relativEt << endl; //continue;
	//This relativeE is used *only* for plotting! Therefore no cuts on this var!
	//create array with multidimensional measurement
	TMeasurement * mess = new TTower;
	mess->pt = double(njet.TowEt[n]);
	double scale = njet.TowEt[n]/njet.TowE[n];
	mess->EMF = double(njet.TowEm[n]*scale);
	mess->HadF = double(njet.TowHad[n]*scale);
	mess->OutF = double(njet.TowOE[n]*scale);
	mess->eta = double(njet.TowEta[n]);
	mess->phi = double(njet.TowPhi[n]);
	mess->E = double(njet.TowE[n]);
	//mess[7] = double( cos( njet.JetCalPhi-njet.TowPhi[n] ) ); // Projection factor for summing tower Pt

	jj_data[nstoredjets]->AddMess(new TData_TruthMess(
	    index,
	    mess,                                                   //mess//
	    njet.JetPt[ij] * relativEt,                             //truth//
	    sqrt(pow(0.5,2)+pow(0.1*njet.JetPt[ij]*relativEt,2)),   //error//
            //1.,                                                   //weight//
	    njet.Weight,                                            //weight//
	    p->GetTowerParRef( index ),                             //parameter//
	    p->GetNumberOfTowerParametersPerBin(),                  //number of free tower param. p. bin//
	    p->tower_parametrization,                               //function//
	    tower_error_param                                      //error param. function//
	  ));
      }
      //Add the jet's tracks to "gj_data":
      double* EfficiencyMap = p->GetEffMap();
      int track_index;
      for (int n=0; n<njet.NobjTrack; ++n){
        if (njet.Track_jetidx[n]!=(int)ij) continue;//look for ij-jet's tracks

	if((njet.TrackTowIdEta[n] == 0) || (njet.TrackTowIdPhi[n] == 0)) {
	  if(njet.TrackPt[n] > 2){
	    std::cerr << "WARNING: eta or phi id of track is zero!\n";
	    continue;
	  }
	  else track_index = 0; //bent low momentum tracks with no HCAL hit
	}
	else
	  track_index = p->GetTrackBin(p->GetTrackEtaBin(njet.TrackTowIdEta[n]),
					   p->GetTrackPhiBin(njet.TrackTowIdPhi[n]));
	if (track_index<0){ cerr<<"WARNING: JJ track_index = " << track_index << endl; continue; }
	//create array with multidimensional measurement
	//TMeasurement * Tmess = new TTrack;
	TTrack * Tmess = new TTrack;
	Tmess->TrackId = int(njet.TrackId[n]);
	Tmess->TowerId = int(njet.TrackTowId[n]);
	Tmess->pt = double(njet.TrackPt[n]);
	double scale = njet.TrackP[n]/njet.TrackPt[n];
	Tmess->EM1 = double(njet.TrackEMC1[n]*scale);
	Tmess->EMF = double(njet.TrackEMC3[n]*scale);
	Tmess->EM5 = double(njet.TrackEMC5[n]*scale);
	Tmess->Had1 = double(njet.TrackHAC1[n]*scale);
	Tmess->HadF = double(njet.TrackHAC3[n]*scale);
	Tmess->Had5 = double(njet.TrackHAC5[n]*scale);
	Tmess->OutF = 0;
	Tmess->DR = double(njet.TrackDR[n]);
	Tmess->DRout = double(njet.TrackDROut[n]);
	Tmess->eta = double(njet.TrackEta[n]);
	Tmess->etaOut = double(njet.TrackEtaOut[n]);
	Tmess->phi = double(njet.TrackPhi[n]);
	Tmess->phiOut = double(njet.TrackPhiOut[n]);
	Tmess->E = double(njet.TrackP[n]);
	Tmess->TrackChi2 = double(njet.TrackChi2[n]);
	Tmess->NValidHits = int(njet.TrackNHits[n]);
	Tmess->TrackQualityT = bool(njet.TrackQualityT[n]);
	Tmess->MuDR = double(njet.MuDR[n]);
	Tmess->MuDE = double(njet.MuDE[n]);
	int TrackEffBin = p->GetTrackEffBin(njet.TrackPt[n],njet.TrackEta[n]);
	Tmess->Efficiency = EfficiencyMap[TrackEffBin];
	//mess[7] = double( cos( njet.JetCalPhi-njet.TowPhi[n] ) ); // Projection factor for summing tower Pt
	//EM+=mess->EMF;
	//F+=mess->pt;
	jj_data[nstoredjets]->AddTrack(new TData_TruthMess(
					      track_index  * p->GetNumberOfTrackParametersPerBin() + p->GetNumberOfTowerParameters() + p->GetNumberOfJetParameters() ,
					      Tmess,                                                    //mess//
					      0,                           //truth//
					      0.05 + 0.00015 * njet.TrackPt[n], //error//
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
    if ( jj_data[injet]->GetMess()->pt > scale*Rel_cut_on_nJet ) goodevent = false;
  }
  if ( (((TJet*)(jj_data[0]->GetMess()))->genPt < GenJetCutLow) || (((TJet*)(jj_data[1]->GetMess()))->genPt  < GenJetCutLow ))  goodevent = false;
  if ( (fabs(jj_data[0]->GetMess()->eta) > Eta_cut_on_jet) ||  (fabs(jj_data[1]->GetMess()->eta) > Eta_cut_on_jet) )  goodevent = false;
  
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
  if( njet.NobjJet < injet ) {
    nNjet_cut++;
    return 0;
  }

  int     njets = 0;  
  double * terr = new double[njet.NobjTow];

  // Order dijets in genjet Et
  int genjetidx[2] = { 0, 1 }; // Store index of two genjets with highest Et
  for(int i = 0; i < njet.NobjJet; ++i) {
    if( njet.GenJetEt[i] > njet.GenJetEt[genjetidx[0]] ) genjetidx[0] = i;
  }
  if( genjetidx[0] == 1 ) genjetidx[1] = 0;
  for(int i = 0; i < njet.NobjJet; ++i) {
    if( i != genjetidx[0] ) {
      if( njet.GenJetEt[i] > njet.GenJetEt[genjetidx[1]] ) genjetidx[1] = i;
    }
  }
  assert( njet.GenJetEt[genjetidx[0]] >= njet.GenJetEt[genjetidx[1]] );
  
  // Loop over two jets with highest genjet Et
  for(int idx = 0; idx < 2; idx++) {
    int i = genjetidx[idx];

    // Cuts
    if( njet.GenJetEt[i] < GenJetCutLow ) {
      nGenJetCutLow++;
      continue;
    } else if( njet.GenJetEt[i] > GenJetCutUp ) {
      nGenJetCutUp++;
      continue;
    } else if( pow(njet.JetEta[i] - njet.GenJetEta[i],2)
	+ pow(njet.JetPhi[i] - njet.GenJetPhi[i],2) > pow(DeltaRMatchingCut,2) ) {
      nDeltaRMatchingCut++;
      continue;
    } else if( fabs(njet.JetEt[i]) < Et_cut_on_jet ) {
      nEt_cut_on_jet++;
      continue;
    } else if( fabs(njet.JetEta[i]) > Eta_cut_on_jet ) {
      nEta_cut_on_jet++;
      continue;
    }

    // Construct event
    double em = 0;
    double had = 0;
    double out = 0;
    double err2 = 0;
    TMeasurement tower;
    double dR = 10;
    int closestTower = 0; 
    for(int n=0; n<njet.NobjTow; ++n){
      if(njet.Tow_jetidx[n] != i) continue;//look for ij-jet's towers

      em += njet.TowEm[n];
      had +=  njet.TowHad[n];
      out +=  njet.TowOE[n];  
      tower.pt = njet.TowEt[n];
      double scale = njet.TowEt[n]/njet.TowE[n];
      tower.EMF = njet.TowEm[n]*scale;
      tower.HadF = njet.TowHad[n]*scale;
      tower.OutF = njet.TowOE[n]*scale;
      tower.eta = njet.TowEta[n];
      tower.phi = njet.TowPhi[n];
      tower.E = njet.TowE[n];
      terr[n] = tower_error_param(&tower.pt,&tower,0); 
      if(terr[n] == 0) {
	//assume toy MC???
	terr[n] = TParameters::toy_tower_error_parametrization(&tower.pt,&tower);
      }
      terr[n] *= terr[n];
      err2 += terr[n];
      double dphi = TVector2::Phi_mpi_pi(njet.JetPhi[i]-tower.phi);
      double dr = sqrt((njet.JetEta[i]-tower.eta)*(njet.JetEta[i]-tower.eta)+
		       dphi*dphi);     
      if(dr < dR) {
	dR = dr;
	closestTower = n;
      }
    }
    // Cuts on hadronic fraction 
    if(had/(had + em) < Had_cut_min) {
      nHad_cut_min++;
      continue;
    } else if(had/(had + em) > Had_cut_max) { 
      nHad_cut_max++;
      continue;
    }
    double factor =  njet.JetEt[i] /  njet.JetE[i];
    tower.pt = njet.JetEt[i];
    tower.EMF = em * factor;
    tower.HadF = had * factor;
    tower.OutF = out * factor;
    tower.eta = njet.JetEta[i];
    tower.phi = njet.JetPhi[i];
    tower.E   = njet.JetE[i];
    double err =  jet_error_param(&tower.pt,&tower,0);
    err2 += err * err;
    Jet *jet;
    if(dataClass == 12) {
      JetWithTowers *jt = 
	new JetWithTowers(njet.JetEt[i],em * factor,had * factor,
			  out * factor,njet.JetE[i],njet.JetEta[i],
			  njet.JetPhi[i],TJet::uds,njet.GenJetEt[i],
			  njet.JetCorrZSP[i],njet.JetCorrJPT[i],
			  njet.JetCorrL2[i],njet.JetCorrL3[i],
			  njet.JetCorrL2[i]*njet.JetCorrL3[i],
			  njet.JetCorrL2[i]*njet.JetCorrL3[i],
			  p->jet_function(njet.TowId_eta[closestTower],
					  njet.TowId_phi[closestTower]),
			  jet_error_param,p->global_jet_function(),Et_cut_on_jet);
      for(int j = 0 ; j < njet.NobjTow ; ++j) {
	if (njet.Tow_jetidx[j]!= i) continue;//look for ij-jet's towers
	double scale = njet.TowEt[j]/njet.TowE[j];
	jt->addTower(njet.TowEt[j],njet.TowEm[j]*scale,
		     njet.TowHad[j]*scale,njet.TowOE[j]*scale,
		     njet.TowE[j],njet.TowEta[j],njet.TowPhi[j],
		     p->tower_function(njet.TowId_eta[i],njet.TowId_phi[i]),
		     tower_error_param);
      }
      jet = jt;
    }
    else { 
      jet = new Jet(njet.JetEt[i],em * factor,had * factor,out * factor,
		    njet.JetE[i],njet.JetEta[i],njet.JetPhi[i],
		    TJet::uds,njet.GenJetEt[i],njet.JetCorrZSP[i],
		    njet.JetCorrJPT[i],njet.JetCorrL2[i],njet.JetCorrL3[i],
		    njet.JetCorrL2[i]*njet.JetCorrL3[i],
		    njet.JetCorrL2[i]*njet.JetCorrL3[i],
		    p->jet_function(njet.TowId_eta[closestTower],
				    njet.TowId_phi[closestTower]),
		    jet_error_param,p->global_jet_function(),Et_cut_on_jet);    
    }
    JetTruthEvent* jte = new JetTruthEvent(jet,njet.GenJetEt[i],1.0);//njet.Weight);
    data.push_back(jte);
    ++njets;
  }     
  delete [] terr;
  return njets;
}
