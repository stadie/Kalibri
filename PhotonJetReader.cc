//  $Id: PhotonJetReader.cc,v 1.18 2009/06/11 17:32:15 mschrode Exp $

#include "PhotonJetReader.h"

#include "CalibData.h"
#include "SmearPhotonJet.h"
#include "JetTruthEvent.h"
#include "Jet.h"
#include "JetWithTowers.h"
#include "ConfigFile.h"
#include "ToyMC.h"
#include "Parameters.h"

#include <iostream>
#include <cstdlib>

#include <TVector2.h>
#include <TLorentzVector.h>



// ----------------------------------------------------------------   
PhotonJetReader::PhotonJetReader(const std::string& configfile, TParameters* p) :
  EventReader(configfile,p)
{
  // Maximum number of read events
  n_gammajet_events     = config->read<int>("use Gamma-Jet events",-1); 
  if(n_gammajet_events == 0) {
    delete config;
    config = 0;
    return ;
  }
  // Cuts
  Et_cut_on_jet     = config->read<double>("Et cut on jet",0.0); 
  Et_cut_on_gamma   = config->read<double>("Et cut on gamma",0.0);
  Rel_cut_on_gamma  = config->read<double>("Relative Rest Jet Cut",0.2);
  GenJetCutLow      = config->read<double>("Et genJet min",0.0);
  GenJetCutUp       = config->read<double>("Et genJet max",10000.0);
  DeltaRMatchingCut = config->read<double>("DeltaR cut on jet matching",0.25);
  Eta_cut_on_jet    = config->read<double>("Eta cut on jet",5.0);
  Had_cut_min       = config->read<double>("Min had fraction",0.07);
  Had_cut_max       = config->read<double>("Max had fraction",0.95);
  // Counter for cutflow
  nEt_cut_on_jet     = 0;
  nEt_cut_on_gamma   = 0;
  nRel_cut_on_gamma  = 0;
  nGenJetCutLow      = 0;    
  nGenJetCutUp       = 0;     
  nDeltaRMatchingCut = 0;
  nEta_cut_on_jet    = 0;  
  nHad_cut_min       = 0;     
  nHad_cut_max       = 0;     

  // Data class
  dataClass = config->read<int>("Gamma-Jet data class", 0);
  if( !( dataClass == 0 || dataClass == 1 || dataClass == 2 || dataClass == 5) ) {
    std::cout << "PhotonJetReader: Unknown data class " << dataClass << ". Using data class 0." << std::endl;
    dataClass = 0;
  }

  // Input files
  string default_tree_name = config->read<string>("Default Tree Name","CalibTree");
  string treename_gammajet = config->read<string>("Gamma-Jet tree", default_tree_name);
  TTree* tchain_gammajet;
  vector<string> input_gammajet = 
    bag_of_string(config->read<string>( "Gamma-Jet input file", "input/gammajet.root" ));
  if(input_gammajet[0] == "toy") {
    std::cout << "generating " << n_gammajet_events << " Gamma-Jet events\n";
    ToyMC* mc = new ToyMC();
    mc->init(configfile);
    mc->print();
    tchain_gammajet = new TTree(treename_gammajet.c_str(),"Gamma Jet events");
    mc->generatePhotonJetTree(tchain_gammajet,n_gammajet_events);
    delete mc;
  } else {
    TChain* chain = new TChain(treename_gammajet.c_str());
    for (bag_of_string::const_iterator it = input_gammajet.begin(); it!=input_gammajet.end(); ++it){
      cout << "...opening root-file " << (*it) << " for Gamma-Jet analysis." << endl;
      chain->Add( it->c_str() );
    }
    tchain_gammajet = chain;
  }
  gammajet.Init( tchain_gammajet );
  
  delete config;
  config = 0;
}



// ----------------------------------------------------------------   
int PhotonJetReader::readEvents(std::vector<TData*>& data)
{
  if(n_gammajet_events == 0) return 0;

  // Reset counters of rejected events
  nEt_cut_on_jet     = 0;
  nEt_cut_on_gamma   = 0;
  nRel_cut_on_gamma  = 0;
  nGenJetCutLow      = 0;    
  nGenJetCutUp       = 0;     
  nDeltaRMatchingCut = 0;
  nEta_cut_on_jet    = 0;  
  nHad_cut_min       = 0;     
  nHad_cut_max       = 0;     

  int nevent    = gammajet.fChain->GetEntries();  // Number of events in chain
  int nReadEvts = 0;                              // Number of read events
  int nGoodEvts = 0;                              // Number of events passing all cuts

  cout << "\nReading Gamma-jet events for data class " << dataClass << "...\n";
  for (int i=0;i<nevent;i++) {
    nReadEvts++;

    if(i%10000==0) cout<<"Gamma-Jet Event: "<<i<<endl;
    gammajet.fChain->GetEvent(i); 
    if (gammajet.NobjTowCal>200) {
      cerr<<"ERROR: Increase array sizes in GammaJetSelector; NobjTowCal="
	  <<gammajet.NobjTowCal<<"!"<<endl;
      exit(8);
    }
 
    // Trivial cuts
    bool goodEvent = true;
    if( gammajet.JetGenEt < GenJetCutLow ) {
      nGenJetCutLow++;
      goodEvent = false;
    } else if( gammajet.JetGenEt > GenJetCutUp ) {
      nGenJetCutUp++;
      goodEvent = false;
    } else if( pow(gammajet.JetCalEta - gammajet.JetGenEta,2)
	       + pow(TVector2::Phi_mpi_pi(gammajet.JetCalPhi - gammajet.JetGenPhi),2)
	       > pow(DeltaRMatchingCut,2) ) {
      nDeltaRMatchingCut++;
      goodEvent = false;
    } else if( gammajet.PhotonEt < Et_cut_on_gamma ) {
      nEt_cut_on_gamma++;
      goodEvent = false;
    } else if( gammajet.JetCalEt < Et_cut_on_jet ) {
      nEt_cut_on_jet++;
      goodEvent = false;
    } else if( gammajet.NonLeadingJetPt / gammajet.PhotonPt > Rel_cut_on_gamma) {
      nRel_cut_on_gamma++;
      goodEvent = false;
    } else if( fabs(gammajet.JetCalEta) > Eta_cut_on_jet ) {
      nEta_cut_on_jet++;
      goodEvent = false;
    }

    if( goodEvent ) {
      
      TData*                                     ev = 0;
      if(dataClass == 0)                         ev = createTruthMultMessEvent();
      else if(dataClass == 1 || dataClass == 2)  ev = createJetTruthEvent();
      else if(dataClass == 5)                    ev = createSmearEvent();
      
      if(ev) {
	data.push_back(ev); 
	nGoodEvts++;
      }
    }

    if((n_gammajet_events >= 0) && (nReadEvts >= n_gammajet_events)) break;
  }

  // Print cut flow
  std::cout << "Read Gamma-jet events:\n";
  std::cout << "  " << (nReadEvts-=nGenJetCutLow) << std::flush;
  std::cout << " gamma-jet events with ptgen > " << GenJetCutLow << "\n";
  std::cout << "  " << (nReadEvts-=nGenJetCutUp) << std::flush;
  std::cout << " gamma-jet events with ptgen < " << GenJetCutUp << "\n";
  std::cout << "  " << (nReadEvts-=nDeltaRMatchingCut) << std::flush;
  std::cout << " gamma-jet events with DeltaR < " << DeltaRMatchingCut << "\n";
  std::cout << "  " << (nReadEvts-=nEt_cut_on_gamma) << std::flush;
  std::cout << " gamma-jet events photon Et > " << Et_cut_on_gamma << "\n";
  std::cout << "  " << (nReadEvts-=nEt_cut_on_jet) << std::flush;
  std::cout << " gamma-jet events jet Et > " << Et_cut_on_jet << "\n";
  std::cout << "  " << (nReadEvts-=nRel_cut_on_gamma) << std::flush;
  std::cout << " gamma-jet events (non-leading jet Et) / (photon Et) > " << Rel_cut_on_gamma << "\n";
  std::cout << "  " << (nReadEvts-=nEta_cut_on_jet) << std::flush;
  std::cout << " gamma-jet events with |eta| < " << Eta_cut_on_jet << "\n";
  std::cout << "  " << (nReadEvts-=nHad_cut_min) << std::flush;
  std::cout << " gamma-jet events with hadronic fraction > " << Had_cut_min << "\n";
  std::cout << "  " << (nReadEvts-=nHad_cut_max) << std::flush;
  std::cout << " gamma-jet events with hadronic fraction < " << Had_cut_max << "\n";
  std::cout << "Stored " << nGoodEvts << " events for analysis.\n";

  return nGoodEvts;
}



// ----------------------------------------------------------------   
TData* PhotonJetReader::createJetTruthEvent()
{
  double em        = 0;
  double had       = 0;
  double out       = 0;
  double dR        = 10;
  int closestTower = 0;
  TMeasurement tower;

  TLorentzVector LJet(0,0,0,0);
  LJet.SetPtEtaPhiE(gammajet.JetCalPt,gammajet.JetCalEta,gammajet.JetCalPhi,gammajet.JetCalE);
  TLorentzVector LGenJet(0,0,0,0);
  LGenJet.SetPtEtaPhiE(gammajet.JetGenPt,gammajet.JetGenEta,gammajet.JetGenPhi,gammajet.JetGenE);

  // Loop over towers, find closest tower to jet axis,
  // and sum up emf, hadf, outf
  for(int n = 0; n < gammajet.NobjTowCal; ++n) {
    em          += gammajet.TowEm[n];
    had         +=  gammajet.TowHad[n];
    out         +=  gammajet.TowOE[n];  
    
    double dphi  = TVector2::Phi_mpi_pi(gammajet.JetCalPhi-tower.phi);
    double dr    = sqrt((gammajet.JetCalEta-tower.eta)*(gammajet.JetCalEta-tower.eta)+
			dphi*dphi);     
    if(dr < dR) {
      dR = dr;
      closestTower = n;
    }
  } // End of loop over towers
  
  // Cuts on hadronic fraction 
  if(had/(had + em) < Had_cut_min) {
    nHad_cut_min++;
    return 0;
  } else if(had/(had + em) > Had_cut_max) { 
    nHad_cut_max++;
    return 0;
  }

  double factor = gammajet.JetCalEt /  gammajet.JetCalE;

  Jet *j;
  if(dataClass == 2) {
    JetWithTowers *jt = 
      new JetWithTowers(gammajet.JetCalEt,em * factor,had * factor,
			out * factor,gammajet.JetCalE,gammajet.JetCalEta,
			gammajet.JetCalPhi,TJet::uds,gammajet.JetGenEt,LJet.DeltaR(LGenJet),
			gammajet.JetCorrZSP,gammajet.JetCorrJPT,
			gammajet.JetCorrL2,gammajet.JetCorrL3,
			gammajet.JetCorrL2L3,gammajet.JetCorrL2L3JPT,
			p->jet_function(gammajet.TowId_eta[closestTower],
					gammajet.TowId_phi[closestTower]),
			jet_error_param,p->global_jet_function(),Et_cut_on_jet);
    for(int i = 0; i < gammajet.NobjTowCal; ++i) {
      double scale = gammajet.TowEt[i]/gammajet.TowE[i];
      jt->addTower(gammajet.TowEt[i],gammajet.TowEm[i]*scale,
		   gammajet.TowHad[i]*scale,gammajet.TowOE[i]*scale,
		   gammajet.TowE[i],gammajet.TowEta[i],gammajet.TowPhi[i],
		   p->tower_function(gammajet.TowId_eta[i],gammajet.TowId_phi[i]),
		   tower_error_param);
    }
    j = jt;
  }
  else { 
    j = new Jet(gammajet.JetCalEt,em * factor,had * factor,out * factor,
		gammajet.JetCalE,gammajet.JetCalEta,gammajet.JetCalPhi,
		TJet::uds,gammajet.JetGenEt,LJet.DeltaR(LGenJet),gammajet.JetCorrZSP,
		gammajet.JetCorrJPT,gammajet.JetCorrL2,gammajet.JetCorrL3,
		gammajet.JetCorrL2L3,gammajet.JetCorrL2L3JPT,
		p->jet_function(gammajet.TowId_eta[closestTower],
				gammajet.TowId_phi[closestTower]),
		jet_error_param,p->global_jet_function(),Et_cut_on_jet);
  }
  JetTruthEvent * jte = new JetTruthEvent(j,gammajet.PhotonEt,gammajet.EventWeight);
  
  return jte;
}



//!  \brief Create SmearPhotonJet event for jet smearing
//!  \note Measured pt is L2L3 corrected
// ----------------------------------------------------------------   
TData* PhotonJetReader::createSmearEvent()
{
  //Find the jets eta & phi index using the nearest tower to jet axis:
  int    jet_index    = -1;
  double min_tower_dr = 10.;
  double em           = 0;
  double had          = 0;
  double out          = 0;
  int    closestTower = 0;

  TLorentzVector LJet(0,0,0,0);
  LJet.SetPtEtaPhiE(gammajet.JetCalEt,gammajet.JetCalEta,gammajet.JetCalPhi,gammajet.JetCalE);
  TLorentzVector LGenJet(0,0,0,0);
  LGenJet.SetPtEtaPhiE(gammajet.JetGenPt,gammajet.JetGenEta,gammajet.JetGenPhi,gammajet.JetGenE);

  for (int n=0; n<gammajet.NobjTowCal; ++n) {
    em  += gammajet.TowEm[n];
    had += gammajet.TowHad[n];
    out += gammajet.TowOE[n];
    TLorentzVector LTower(0,0,0,0);
    LTower.SetPtEtaPhiE(gammajet.TowEt[n],gammajet.TowEta[n],gammajet.TowPhi[n],gammajet.TowE[n]);
    double dr = LTower.DeltaR(LJet);
    if (dr<min_tower_dr) {
      min_tower_dr = dr;
      closestTower = n;
    }
  }

  if (jet_index<0){ cerr<<"WARNING: jet_index = " << jet_index << endl; return 0; }
  if(had/(had + em) < Had_cut_min) { return 0;}
  if(had/(had + em) > Had_cut_max) { return 0;}

  // Set up measurement
  TJet * jet      = new TJet;
  jet->pt         = gammajet.JetCorrL2 * gammajet.JetCorrL3 * gammajet.JetCalEt;
  jet->eta        = gammajet.JetCalEta;
  jet->phi        = gammajet.JetCalPhi;
  jet->E          = gammajet.JetCalE;
  jet->genPt      = gammajet.JetGenPt;
  jet->dR         = LJet.DeltaR(LGenJet);
  jet->ZSPcor     = gammajet.JetCorrZSP; 
  jet->JPTcor     = gammajet.JetCorrJPT; 
  jet->L2cor      = gammajet.JetCorrL2; 
  jet->L3cor      = gammajet.JetCorrL3; 
  jet->L2L3cor    = gammajet.JetCorrL2 * gammajet.JetCorrL3; 
  jet->L2L3JPTcor = 1.;//gammajet.JetCorrL2L3JPT[ij]; 
  //the following is not quite correct, as this factor is different for all towers. These values should be in the n-tupel as well
  double factor    = gammajet.JetCalEt /  gammajet.JetCalE;
  jet->HadF       = had * factor;
  jet->EMF        = em * factor;
  jet->OutF       = out * factor;

  SmearPhotonJet * pje = new SmearPhotonJet(jet,gammajet.PhotonEt,1.,
					    p->jet_function(gammajet.TowId_eta[closestTower],
							    gammajet.TowId_phi[closestTower]));

  return pje;
}



TData* PhotonJetReader::createTruthMultMessEvent() 
{
    //Find the jets eta & phi index using the nearest tower to jet axis:
    int jet_index=-1;
    double min_tower_dr = 10.0;
    double em = 0;
    double had = 0;
    double out = 0;

    TLorentzVector LJet(0,0,0,0);
    LJet.SetPtEtaPhiE(gammajet.JetCalEt,gammajet.JetCalEta,gammajet.JetCalPhi,gammajet.JetCalE);
    TLorentzVector LGenJet(0,0,0,0);
    LGenJet.SetPtEtaPhiE(gammajet.JetGenPt,gammajet.JetGenEta,gammajet.JetGenPhi,gammajet.JetGenE);

    /*
    for (int n=0; n<gammajet.NobjTowCal; ++n){
      TLorentzVector Ltower(0,0,0,0);
      Ltower.SetPtEtaPhiE(gammajet.TowEt[n],gammajet.TowEta[n],gammajet.TowPhi[n],gammajet.TowE[n]);
      Ljet += Ltower;
    }
    */
    //Ljet.SetPtEtaPhiE(gammajet.JetCalEt,gammajet.JetCalEta,gammajet.JetCalPhi,gammajet.JetCalE);
    for (int n=0; n<gammajet.NobjTowCal; ++n){
      em += gammajet.TowEm[n];
      had +=  gammajet.TowHad[n];
      out +=  gammajet.TowOE[n];
      TLorentzVector Ltower(0,0,0,0);
      Ltower.SetPtEtaPhiE(gammajet.TowEt[n],gammajet.TowEta[n],gammajet.TowPhi[n],gammajet.TowE[n]);
      double dr = Ltower.DeltaR(LJet);
      if (dr<min_tower_dr) {
	jet_index = p->GetJetBin(p->GetJetEtaBin(gammajet.TowId_eta[n]),
				 p->GetJetPhiBin(gammajet.TowId_phi[n]));
	min_tower_dr = dr;
      }
    }
    if (jet_index<0){ cerr<<"WARNING: jet_index = " << jet_index << endl; return 0; }
    if(had/(had + em) < Had_cut_min) { return 0;}
    if(had/(had + em) > Had_cut_max) { return 0;}
    //jet_index: p->eta_granularity*p->phi_granularity*p->GetNumberOfTowerParametersPerBin()
    //           has to be added for a correct reference to k[...].

    TJet* jetp  = new TJet;
    jetp->pt  = gammajet.JetCalEt;
    jetp->eta = gammajet.JetCalEta;
    jetp->phi = gammajet.JetCalPhi;
    jetp->E   = gammajet.JetCalE;
    jetp->genPt =gammajet.JetGenPt;
    jetp->dR    = LJet.DeltaR(LGenJet);
    jetp->ZSPcor =gammajet.JetCorrZSP; 
    jetp->JPTcor =gammajet.JetCorrJPT; 
    jetp->L2cor =gammajet.JetCorrL2; 
    jetp->L3cor =gammajet.JetCorrL3; 
    jetp->L2L3cor =gammajet.JetCorrL2L3; 
    jetp->L2L3JPTcor =gammajet.JetCorrL2L3JPT; 
    //the following is not quite correct, as this factor is different for all towers. These values should be in the n-tupel as well
    double factor =  gammajet.JetCalEt /  gammajet.JetCalE;
    jetp->HadF = had * factor;
    jetp->EMF = em * factor;
    jetp->OutF = out * factor;

    //Create an Gamma/Jet TData event
    TData_TruthMultMess * gj_data = new TData_TruthMultMess
      (
       jet_index  * p->GetNumberOfJetParametersPerBin() + p->GetNumberOfTowerParameters(),
       //gammajet.PhotonEt,				    //truth//
       gammajet.JetGenPt,
       sqrt(pow(0.5,2)+pow(0.10*gammajet.PhotonEt,2)),   //error//
       //0.10*gammajet.PhotonEt.//error//				    
       gammajet.EventWeight,                             //weight//
       //1.0,                                            //weight//
       p->GetJetParRef( jet_index ),                     //params
       p->GetNumberOfJetParametersPerBin(),              //number of free jet param. p. bin
       p->jet_parametrization,                           //function
       jet_error_param,                                  //error param. function
       jetp                                              //measurement
       );

    //double EM=0.,F=0.;
    //Add the jet's towers to "gj_data":
    for (int n=0; n<gammajet.NobjTowCal; ++n){
      //if (gammajet.TowEt[n]<0.01) continue;
   
      int index = p->GetBin(p->GetEtaBin(gammajet.TowId_eta[n]),
			    p->GetPhiBin(gammajet.TowId_phi[n]));
      if (index<0){ cerr<<"WARNING: towewer_index = " << index << endl; continue; }

      //double dR = deltaR(gammajet.JetCalEta, gammajet.JetCalPhi, gammajet.TowEta[n], gammajet.TowPhi[n]);
	      
      double relativEt = gammajet.TowEt[n]/gammajet.JetCalEt;  
      //if (relativEt<=0) cerr << "relEt = " <<relativEt << endl; //continue;
      //This relativeE is used *only* for plotting! Therefore no cuts on this var!
      //create array with multidimensional measurement
      TMeasurement * mess = new TTower;
      mess->pt = double(gammajet.TowEt[n]);
      double scale = gammajet.TowEt[n]/gammajet.TowE[n];
      mess->EMF = double(gammajet.TowEm[n]*scale);
      mess->HadF = double(gammajet.TowHad[n]*scale);
      mess->OutF = double(gammajet.TowOE[n]*scale);
      mess->eta = double(gammajet.TowEta[n]);
      mess->phi = double(gammajet.TowPhi[n]);
      mess->E = double(gammajet.TowE[n]);
      //mess[7] = double( cos( gammajet.JetCalPhi-gammajet.TowPhi[n] ) ); // Projection factor for summing tower Pt
      //EM+=mess->EMF;
      //F+=mess->pt;
      gj_data->AddMess(new TData_TruthMess(index,
					   mess,                                                    //mess//
					   gammajet.PhotonEt * relativEt,                           //truth//
					   //sqrt(1.3 * 1.3/gammajet.TowHad[n] + 0.056 * 0.056) * mess->HadF,
					   sqrt(pow(0.5,2)+pow(0.1*gammajet.PhotonEt*relativEt,2)), //error//
					   1.,                                                      //weight//
					   p->GetTowerParRef( index ),                              //parameter//
					   p->GetNumberOfTowerParametersPerBin(),                   //number of free tower param. p. bin//
					   p->tower_parametrization,                                //function//
					   tower_error_param                                       //error param.func.//
					   ));
    } 

    //Add the jet's tracks to "gj_data":
    double* EfficiencyMap = p->GetEffMap();
    int track_index;
    for (int n=0; n<gammajet.NobjTrack; ++n){
      if((gammajet.TrackTowIdEta[n] == 0) || (gammajet.TrackTowIdPhi[n] == 0)) {
	if(gammajet.TrackPt[n] > 2){
	  std::cerr << "WARNING: eta or phi id of track is zero!\n";
	  continue;
	}
	else track_index = 0; //bent low momentum tracks with no HCAL hit
      }
      else
	//one trackindex for all tracks in jet = track_index
	track_index = p->GetTrackBin(p->GetTrackEtaBin(gammajet.TrackTowIdEta[n]),
					 p->GetTrackPhiBin(gammajet.TrackTowIdPhi[n]));
      if (track_index<0){ cerr<<"WARNING: track_index = " << track_index << endl; continue; }
      //create array with multidimensional measurement
      //TMeasurement * Tmess = new TTrack;
      TTrack * Tmess = new TTrack;
      Tmess->TrackId = int(gammajet.TrackId[n]);
      Tmess->TowerId = int(gammajet.TrackTowId[n]);
      Tmess->pt = double(gammajet.TrackPt[n]);
      double scale = gammajet.TrackP[n]/gammajet.TrackPt[n];
      Tmess->EM1 = double(gammajet.TrackEMC1[n]*scale);
      Tmess->EMF = double(gammajet.TrackEMC3[n]*scale);
      Tmess->EM5 = double(gammajet.TrackEMC5[n]*scale);
      Tmess->Had1 = double(gammajet.TrackHAC1[n]*scale);
      Tmess->HadF = double(gammajet.TrackHAC3[n]*scale);
      Tmess->Had5 = double(gammajet.TrackHAC5[n]*scale);
      Tmess->OutF = 0;
      Tmess->DR = double(gammajet.TrackDR[n]);
      Tmess->DRout = double(gammajet.TrackDROut[n]);
      Tmess->eta = double(gammajet.TrackEta[n]);
      Tmess->etaOut = double(gammajet.TrackEtaOut[n]);
      Tmess->phi = double(gammajet.TrackPhi[n]);
      Tmess->phiOut = double(gammajet.TrackPhiOut[n]);
      Tmess->E = double(gammajet.TrackP[n]);
      Tmess->TrackChi2 = double(gammajet.TrackChi2[n]);
      Tmess->NValidHits = int(gammajet.TrackNHits[n]);
      Tmess->TrackQualityT = bool(gammajet.TrackQualityT[n]);
      Tmess->MuDR = double(gammajet.MuDR[n]);
      Tmess->MuDE = double(gammajet.MuDE[n]);
      int TrackEffBin = p->GetTrackEffBin(gammajet.TrackPt[n],gammajet.TrackEta[n]);
      Tmess->Efficiency = EfficiencyMap[TrackEffBin];
      //mess[7] = double( cos( gammajet.JetCalPhi-gammajet.TowPhi[n] ) ); // Projection factor for summing tower Pt
      //EM+=mess->EMF;
      //F+=mess->pt;
      gj_data->AddTrack(new TData_TruthMess(
					    track_index  * p->GetNumberOfTrackParametersPerBin() + p->GetNumberOfTowerParameters() + p->GetNumberOfJetParameters() ,
					   Tmess,                                                    //mess//
					   0,                           //truth//
					   0.05 + 0.00015 * gammajet.TrackPt[n], //error//
					   1.,                                                      //weight//
					   p->GetTrackParRef( track_index ),                              //parameter//
					   p->GetNumberOfTrackParametersPerBin(),                   //number of free tower param. p. bin//
					   p->track_parametrization,                                //function//
					   track_error_param                                        //error param.func.//
					   ));
    }
    gj_data->UseTracks(useTracks);   //check if track information is sufficient to use Track Parametrization
    return  gj_data;
}
