//
//    Reader for Z Jet Events
//
//    This class reads events according fo the ZJetSel
//
//    first version: Hartmut Stadie 2008/12/12
//    $Id: ZJetReader.cc,v 1.3 2009/01/06 13:35:21 stadie Exp $
//   
#include "ZJetReader.h"

#include "CalibData.h"
#include "ConfigFile.h"
#include "Parameters.h"
#include "TLorentzVector.h"

#include <cstdlib>
#include <iostream>

ZJetReader::ZJetReader(const std::string& configfile, TParameters* p) :
  EventReader(configfile,p), Et_cut_on_Z(0),Et_cut_on_jet(0)
{
  n_zjet_events     = config->read<int>("use Z-Jet events",-1); 
  if(n_zjet_events == 0) {
    delete config;
    config = 0;
    return ;
  }

  Et_cut_on_Z     = config->read<double>("Et cut on Z",0.0); 
  Et_cut_on_jet   = config->read<double>("Et cut on jet",0.0);
  
  string default_tree_name = config->read<string>("Default Tree Name","CalibTree");
  string treename_zjet      = config->read<string>( "Z-Jet tree", default_tree_name );
  TChain * tchain_zjet      = new TChain( treename_zjet.c_str() );
  vector<string> input_zjet = bag_of_string( config->read<string>( "Z-Jet input file", "input/zjet.root" ) );
  for (bag_of_string::const_iterator it = input_zjet.begin(); it!=input_zjet.end(); ++it){
    cout << "...opening root-file " << (*it) << " for Z-Jet analysis." << endl;
    tchain_zjet->Add( it->c_str() );
  }  
  zjet.Init( tchain_zjet );
  
  delete config;
  config = 0;
}

ZJetReader::~ZJetReader()
{
}

//calculates from Z energy a truth value for one calo tower of the jet.
int ZJetReader::readEvents(std::vector<TData*>& data)
{
  if(n_zjet_events == 0) return 0;
  //Run Z-Jet stuff  
  int nevent = zjet.fChain->GetEntries();
  int nevents_added = 0;
  for (int i=0;i<nevent;i++) {
    if(i%1000==0) cout<<"Z-Jet Event: "<<i<<endl;
    zjet.fChain->GetEvent(i); 
    if (zjet.NobjTowCal>200) {
      cerr<<"ERROR: Increase array sizes in ZJetSelector; NobjTowCal="
	  <<zjet.NobjTowCal<<"!"<<endl;
      exit(8);
    }
    if (zjet.NobjETowCal>200) {
      cerr<<"ERROR: Increase array sizes in ZJetSelector; NobjETowCal="
	  <<zjet.NobjETowCal<<"!"<<endl;
      exit(8);
    }   
    if (zjet.NobjTrack>200) {
      cerr<<"ERROR: Increase array sizes in ZJetSelector; NobjTrackCal="
	  <<zjet.NobjTrack<<"!"<<endl;
      exit(8);
    }
    //trivial cuts
    if (zjet.ZPt<Et_cut_on_Z || zjet.JetCalPt<Et_cut_on_jet) continue;
     
    //Find the jets eta & phi index using the nearest tower to jet axis:
    int jet_index=-1;
    double min_tower_dr = 10.0;
    double em = 0;
    double had = 0;
    TLorentzVector Ljet(0,0,0,0);
    Ljet.SetPtEtaPhiE(zjet.JetCalEt,zjet.JetCalEta,zjet.JetCalPhi,zjet.JetCalE);
    for (int n=0; n<zjet.NobjTowCal; ++n){
      em += zjet.TowEm[n];
      had +=  zjet.TowHad[n];
      TLorentzVector Ltower(0,0,0,0);
      Ltower.SetPtEtaPhiE(zjet.TowEt[n],zjet.TowEta[n],zjet.TowPhi[n],zjet.TowE[n]);
      double dr = Ltower.DeltaR(Ljet);
      if (dr<min_tower_dr) {
	jet_index = p->GetJetBin(p->GetJetEtaBin(zjet.TowId_eta[n]),
				 p->GetJetPhiBin(zjet.TowId_phi[n]));
	min_tower_dr = dr;
      }
    }

    if (jet_index<0){ cerr<<"WARNING: jet_index = " << jet_index << endl; continue; }
    if(em == 0) { continue;}
    //jet_index: p->eta_granularity*p->phi_granularity*p->GetNumberOfTowerParametersPerBin()
    //           has to be added for a correct reference to k[...].

    TMeasurement* jetp  = new TJet;
    jetp->pt  = zjet.JetCalEt;
    jetp->eta = zjet.JetCalEta;
    jetp->phi = zjet.JetCalPhi;
    jetp->E   = zjet.JetCalE;
    //Create an Z/Jet TData event
    TData_TruthMultMess * gj_data = new 
      TData_TruthMultMess(jet_index  * p->GetNumberOfJetParametersPerBin() + p->GetNumberOfTowerParameters(),
			  zjet.ZPt,				    //truth//
			  //zjet.JetGenPt,
			  sqrt(pow(0.5,2)+pow(0.10*zjet.ZEt,2)),    //error//
			  //zjet.EventWeight,                       //weight//
			  1.0,                                      //weight//
			  p->GetJetParRef( jet_index ),             //params
			  p->GetNumberOfJetParametersPerBin(),      //number of free jet param. p. bin
			  p->jet_parametrization,                   //function
			  jet_error_param,                          //error param. function
			  jetp
			  );
    
    //Add the jet's towers to "gj_data":
    for (int n=0; n<zjet.NobjTowCal; ++n){
      //if (zjet.TowEt[n]<0.01) continue;
      
      int index = p->GetBin(p->GetEtaBin(zjet.TowId_eta[n]),
			    p->GetPhiBin(zjet.TowId_phi[n]));
      if (index<0){ cerr<<"WARNING: towewer_index = " << index << endl; continue; }
      
      //double dR = deltaR(zjet.JetCalEta, zjet.JetCalPhi, zjet.TowEta[n], zjet.TowPhi[n]);
      
      double relativEt = zjet.TowEt[n]/zjet.JetCalEt;  
      //if (relativEt<=0) cerr << "relEt = " <<relativEt << endl; //continue;
      //This relativeE is used *only* for plotting! Therefore no cuts on this var!
      //create array with multidimensional measurement
      TMeasurement * mess = new TTower;
      mess->pt = double(zjet.TowEt[n]);
      double scale = zjet.TowEt[n]/zjet.TowE[n];
      mess->EMF = double(zjet.TowEm[n]*scale);
      mess->HadF = double(zjet.TowHad[n]*scale);
      mess->OutF = double(zjet.TowOE[n]*scale);
      mess->eta = double(zjet.TowEta[n]);
      mess->phi = double(zjet.TowPhi[n]);
      mess->E = double(zjet.TowE[n]);
      //mess[7] = double( cos( zjet.JetCalPhi-zjet.TowPhi[n] ) ); // Projection factor for summing tower Pt
      gj_data->AddMess(new TData_TruthMess(index,
					   mess,                                           //mess//
					   zjet.ZEt * relativEt,                           //truth//
					   sqrt(pow(0.5,2)+pow(0.1*zjet.ZEt*relativEt,2)), //error//
					   1.,                                             //weight//
					   p->GetTowerParRef( index ),                     //parameter//
					   p->GetNumberOfTowerParametersPerBin(),          //number of free tower param. p. bin//
					   p->tower_parametrization,                       //function//
					   tower_error_param                                        //error param.func.//
					   ));
    } 

    
    //Add the jet's tracks to "gj_data":
    for (int n=0; n<zjet.NobjTrack; ++n){
      if((zjet.TrackTowIdEta[n] == 0) || (zjet.TrackTowIdPhi[n] == 0)) {
	std::cerr << "WARNING: eta or phi id of track is zero!\n";
	continue;
      }
      int index = p->GetTrackBin(p->GetTrackEtaBin(zjet.TrackTowIdEta[n]),
				 p->GetTrackPhiBin(zjet.TrackTowIdPhi[n]));
      if (index<0){ cerr<<"WARNING: track_index = " << index << endl; continue; }
      //create array with multidimensional measurement
      //TMeasurement * Tmess = new TTrack;
      TTrack * Tmess = new TTrack;
      Tmess->TrackId = int(zjet.TrackId[n]);
      Tmess->TowerId = int(zjet.TrackTowId[n]);
      Tmess->pt = double(zjet.TrackPt[n]);
      double scale = zjet.TrackP[n]/zjet.TrackPt[n];
      Tmess->EM1 = double(zjet.TrackEMC1[n]*scale);
      Tmess->EMF = double(zjet.TrackEMC3[n]*scale);
      Tmess->EM5 = double(zjet.TrackEMC5[n]*scale);
      Tmess->Had1 = double(zjet.TrackHAC1[n]*scale);
      Tmess->HadF = double(zjet.TrackHAC3[n]*scale);
      Tmess->Had5 = double(zjet.TrackHAC5[n]*scale);
      Tmess->OutF = 0;
      Tmess->DR = double(zjet.TrackDR[n]);
      Tmess->DRout = double(zjet.TrackDROut[n]);
      Tmess->eta = double(zjet.TrackEta[n]);
      Tmess->etaOut = double(zjet.TrackEtaOut[n]);
      Tmess->phi = double(zjet.TrackPhi[n]);
      Tmess->phiOut = double(zjet.TrackPhiOut[n]);
      Tmess->E = double(zjet.TrackP[n]);
      Tmess->TrackChi2 = double(zjet.TrackChi2[n]);
      Tmess->NValidHits = int(zjet.TrackNHits[n]);
      Tmess->MuDR = double(zjet.MuDR[n]);
      Tmess->MuDE = double(zjet.MuDE[n]);
      //mess[7] = double( cos( zjet.JetCalPhi-zjet.TowPhi[n] ) ); // Projection factor for summing tower Pt
      //EM+=mess->EMF;
      //F+=mess->pt;
       gj_data->AddTrack(new TData_TruthMess(index  * p->GetNumberOfTrackParametersPerBin() + p->GetNumberOfTowerParameters() + p->GetNumberOfJetParameters() ,
					   Tmess,                                                    //mess//
					   0,                           //truth//
					   0.015 * zjet.TrackPt[n], //error//
					   1.,                                                      //weight//
					   p->GetTrackParRef( index ),                              //parameter//
					   p->GetNumberOfTrackParametersPerBin(),                   //number of free tower param. p. bin//
					   p->track_parametrization,                                //function//
					   track_error_param                                        //error param.func.//
					   ));
    } 
    gj_data->UseTracks(useTracks);   //check if track information is sufficient to use Track Parametrization
    data.push_back( gj_data ); 
    ++nevents_added;
    if((n_zjet_events>=0) && (nevents_added >= n_zjet_events-1))
      break;
  }
  return nevents_added;
}

