//
//    Reader for Photon Jet Events
//
//    This class reads events according fo the GammaJetSel
//
//    first version: Hartmut Stadie 2008/12/12
//    $Id: PhotonJetReader.cc,v 1.14 2009/03/04 17:26:51 thomsen Exp $
//   
#include "PhotonJetReader.h"

#include "CalibData.h"
#include "JetTruthEvent.h"
#include "Jet.h"
#include "JetWithTowers.h"
#include "ConfigFile.h"
#include "ToyMC.h"
#include "Parameters.h"

#include <iostream>
#include <cstdlib>

PhotonJetReader::PhotonJetReader(const std::string& configfile, TParameters* p) :
  EventReader(configfile,p), Et_cut_on_gamma(0),Et_cut_on_jet(0),
  Rel_cut_on_gamma(10),n_gammajet_events(0),dataClass(0)
{
  n_gammajet_events     = config->read<int>("use Gamma-Jet events",-1); 
  if(n_gammajet_events == 0) {
    delete config;
    config = 0;
    return ;
  }
  string default_tree_name = config->read<string>("Default Tree Name","CalibTree");
  string treename_gammajet = config->read<string>("Gamma-Jet tree", default_tree_name);
  
  Et_cut_on_jet   = config->read<double>("Et cut on jet",0.0); 
  Eta_cut_on_jet  = config->read<double>("Eta cut on jet",5.0); 
  Et_cut_on_gamma = config->read<double>("Et cut on gamma",0.0); 
  Et_cut_on_genJet   = config->read<double>("Et cut on genJet",0.0);
  Rel_cut_on_gamma =  config->read<double>("Relative Rest Jet Cut",0.2);
  
  dataClass = config->read<int>("Gamma-Jet data class", 0);
  if((dataClass < 0) || (dataClass > 2)) dataClass = 0;

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

PhotonJetReader::~PhotonJetReader()
{
  
}

int PhotonJetReader::readEvents(std::vector<TData*>& data)
{
  if(n_gammajet_events == 0) return 0;
  int nevent = gammajet.fChain->GetEntries();
  int nevents_added = 0;
  for (int i=0;i<nevent;i++) {
    if(i%1000==0) cout<<"Gamma-Jet Event: "<<i<<endl;
    gammajet.fChain->GetEvent(i); 
    if (gammajet.NobjTowCal>200) {
      cerr<<"ERROR: Increase array sizes in GammaJetSelector; NobjTowCal="
	  <<gammajet.NobjTowCal<<"!"<<endl;
      exit(8);
    }
 

    //trivial cuts
    if (gammajet.PhotonEt<Et_cut_on_gamma || 
        //gammajet.JetGenPt<Et_cut_on_gamma ||
	gammajet.JetGenPt<Et_cut_on_genJet ||
        gammajet.JetCalPt<Et_cut_on_jet) continue;

    if(fabs(gammajet.JetCalEta) > Eta_cut_on_jet) continue;

    if (gammajet.NonLeadingJetPt   / gammajet.PhotonPt >  Rel_cut_on_gamma)  continue;    //fraction of unwanted stuff

    TData* ev = 0;
    if(dataClass == 0) ev = createTruthMultMessEvent();
    else if(dataClass > 0) ev = createJetTruthEvent();
    
    if(ev)  
      {
	data.push_back(ev); 
	++nevents_added;
	if((n_gammajet_events >= 0) && (nevents_added >= n_gammajet_events-1))
	  break;
      }
  }
  return nevents_added;
}

TData* PhotonJetReader::createJetTruthEvent()
{
  double em = 0;
  double had = 0;
  double out = 0;
  double err2 = 0;
  TMeasurement tower;
  double* terr = new double[gammajet.NobjTowCal];
  double dR = 10;
  int closestTower = 0;
  for(int n = 0; n < gammajet.NobjTowCal; ++n) {
    em += gammajet.TowEm[n];
    had +=  gammajet.TowHad[n];
    out +=  gammajet.TowOE[n];  
    tower.pt = gammajet.TowEt[n];
    double scale = gammajet.TowEt[n]/gammajet.TowE[n];
    tower.EMF = gammajet.TowEm[n]*scale;
    tower.HadF = gammajet.TowHad[n]*scale;
    tower.OutF = gammajet.TowOE[n]*scale;
    tower.eta = gammajet.TowEta[n];
    tower.phi = gammajet.TowPhi[n];
    tower.E = gammajet.TowE[n];
    terr[n] = tower_error_param(&tower.pt,&tower,0); 
    if(terr[n] == 0) {
      //assume toy MC???
      terr[n] = TParameters::toy_tower_error_parametrization(&tower.pt,&tower);
    }
    terr[n] *= terr[n];
    err2 += terr[n];
    double dphi = TVector2::Phi_mpi_pi(gammajet.JetCalPhi-tower.phi);
    double dr = sqrt((gammajet.JetCalEta-tower.eta)*(gammajet.JetCalEta-tower.eta)+
		     dphi*dphi);     
    if(dr < dR) {
      dR = dr;
      closestTower = n;
    }
  }  //calc jet error
  if(had/(had + em) < 0.07) { return 0;}
  if(had/(had + em) > 0.92) { return 0;}
  double factor =  gammajet.JetCalEt /  gammajet.JetCalE;
  tower.pt = gammajet.JetCalEt;
  tower.EMF = em * factor;
  tower.HadF = had * factor;
  tower.OutF = out * factor;
  tower.eta = gammajet.JetCalEta;
  tower.phi = gammajet.JetCalPhi;
  tower.E   = gammajet.JetCalE;
  double err =  jet_error_param(&tower.pt,&tower,0);
  err2 += err * err;
  //use first tower, as the towers should be sorted in E or Et 
  Jet *j;
  if(dataClass == 2) {
    JetWithTowers *jt = 
      new JetWithTowers(gammajet.JetCalEt,em * factor,had * factor,
			out * factor,gammajet.JetCalE,gammajet.JetCalEta,
			gammajet.JetCalPhi,TJet::uds,gammajet.JetGenEt,
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
		TJet::uds,gammajet.JetGenEt,gammajet.JetCorrZSP,
		gammajet.JetCorrJPT,gammajet.JetCorrL2,gammajet.JetCorrL3,
		gammajet.JetCorrL2L3,gammajet.JetCorrL2L3JPT,
		p->jet_function(gammajet.TowId_eta[closestTower],
				gammajet.TowId_phi[closestTower]),
		jet_error_param,p->global_jet_function(),Et_cut_on_jet);
  }
  JetTruthEvent* jte = new JetTruthEvent(j,gammajet.PhotonEt,gammajet.EventWeight);
  delete [] terr;
  return jte;
}

TData* PhotonJetReader::createTruthMultMessEvent() 
{
    //Find the jets eta & phi index using the nearest tower to jet axis:
    int jet_index=-1;
    double min_tower_dr = 10.0;
    double em = 0;
    double had = 0;
    double out = 0;
    TLorentzVector Ljet(0,0,0,0);
    Ljet.SetPtEtaPhiE(gammajet.JetCalEt,gammajet.JetCalEta,gammajet.JetCalPhi,gammajet.JetCalE);
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
      double dr = Ltower.DeltaR(Ljet);
      if (dr<min_tower_dr) {
	jet_index = p->GetJetBin(p->GetJetEtaBin(gammajet.TowId_eta[n]),
				 p->GetJetPhiBin(gammajet.TowId_phi[n]));
	min_tower_dr = dr;
      }
    }
    if (jet_index<0){ cerr<<"WARNING: jet_index = " << jet_index << endl; return 0; }
    if(had/(had + em) < 0.07) { return 0;}
    //if(had/(had + em) > 0.92) { continue;}
    //jet_index: p->eta_granularity*p->phi_granularity*p->GetNumberOfTowerParametersPerBin()
    //           has to be added for a correct reference to k[...].

    TJet* jetp  = new TJet;
    jetp->pt  = gammajet.JetCalEt;
    jetp->eta = gammajet.JetCalEta;
    jetp->phi = gammajet.JetCalPhi;
    jetp->E   = gammajet.JetCalE;
    jetp->genPt =gammajet.JetGenPt;
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
