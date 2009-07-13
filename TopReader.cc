//
//    Reader for ttbar events
//
//    This class reads events according fo the TopSel
//
//    first version: Hartmut Stadie 2008/12/12
//    $Id: TopReader.cc,v 1.6 2009/06/11 17:32:15 mschrode Exp $
//   
#include "TopReader.h"

#include "CalibData.h"
#include "ConfigFile.h"
#include "Parameters.h"
#include "TwoJetsInvMassEvent.h"
#include "Jet.h"
#include "JetWithTowers.h"

#include "TLorentzVector.h"

#include <iostream>
#include <cstdlib>

TopReader::TopReader(const std::string& configfile, TParameters* p) :
  EventReader(configfile,p), Et_cut_on_jet(0),useMassConstraintW(false),
  useMassConstraintTop(false),massConstraint_W(80.4),massConstraint_Top(172.4),
  dataClass(0)
{
  n_top_events     = config->read<int>("use Top events",-1); 
  if(n_top_events == 0) {
    delete config;
    config = 0;
    return ;
  }

  Et_cut_on_jet   = config->read<double>("Et cut on jet",0.0);
  useMassConstraintW   = config->read<bool>("use W mass constraint",true);
  useMassConstraintTop = config->read<bool>("use Top mass constraint",false);
  massConstraint_W   = config->read<double>("W mass",80.4);
  massConstraint_Top = config->read<double>("Top mass",172.4);
  if(!useMassConstraintW && !useMassConstraintTop) {
    std::cout << "W mass constraint and Top mass constraint were both turned off," << std::endl
	      << "you should enable at least one of these constraints if you want to run with Top...!" << std::endl;
  }
  dataClass = config->read<int>("Top data class", 0);
  if((dataClass != 0)&&(dataClass != 1)&&(dataClass != 2)) dataClass = 0;
  
  string default_tree_name = config->read<string>("Default Tree Name","CalibTree");
  string treename_top    = config->read<string>( "Top tree", default_tree_name );
  TChain * tchain_top = new TChain( treename_top.c_str() );
  vector<string> input_top = bag_of_string( config->read<string>( "Top input file", "input/top.root" ) );
  for (bag_of_string::const_iterator it = input_top.begin(); it!=input_top.end(); ++it){
    cout << "...opening root-file " << (*it) << " for Top analysis." << endl;
    tchain_top->Add( it->c_str() );
  }  
  top.Init( tchain_top );
  
  delete config;
  config = 0;
}

TopReader::~TopReader()
{
}

int TopReader::readEvents(std::vector<TData*>& data)
{
  if(n_top_events == 0) return 0;
  //Run Top stuff  
  int nevent = top.fChain->GetEntries();
  int evt=0;
  
  for (int i=0;i<nevent;i++) {
    if((i+1)%1000==0) 
      cout<<"Top Event: "<<i+1<<endl;
    top.fChain->GetEvent(i); 
    if (top.NobjTow>1000 || top.NobjJet>8) {
      cerr << "ERROR: Increase array sizes in topSelector; NobjTow="
	   << top.NobjTow<<", NobjBJet="<<top.NobjJet<<"!"<<endl;
      exit(10);
    }
    //--------------
    
    if(useMassConstraintW && (dataClass == 0)) {
      bool goodevent = false;
      TData_InvMass2 * top_data[2]; //two W-jets
      top_data[0] = 0;
      int nstoredjets = 0;
      for (unsigned int ij = 0; ij<3; ++ij){
	if(top.JetPt[ij] < Et_cut_on_jet) continue;
	if((TJet::Flavor)top.JetFlavor[ij] != TJet::uds) continue;
	
	//Find the jets eta & phi index using the nearest tower to jet axis:
	int jet_index=-1;
	double min_tower_dr = 10.0;
	double em = 0;
	double had = 0;
	double out = 0;
	TLorentzVector Ljet(0,0,0,0);
	
	Ljet.SetPtEtaPhiE(top.JetPt[ij],top.JetEta[ij],top.JetPhi[ij],top.JetE[ij]);
	for (int n=0; n<top.NobjTow; ++n){
	  if (top.Tow_jetidx[n]!=(int)ij) continue;//look for ij-jet's towers
	  em += top.TowEm[n];
	  had += top.TowHad[n];
	  out += top.TowOE[n];
	  TLorentzVector Ltower(0,0,0,0);
	  Ltower.SetPtEtaPhiE(top.TowEt[n],top.TowEta[n],top.TowPhi[n],top.TowE[n]);
	  double dr = Ltower.DeltaR(Ljet);
	  if (dr<min_tower_dr) {
	    jet_index = p->GetJetBin(p->GetJetEtaBin(top.TowId_eta[n]),
				     p->GetJetPhiBin(top.TowId_phi[n]));
	    min_tower_dr = dr;
	  }
	}
	if (jet_index<0){
	  cerr<<"WARNING: Top jet_index = " << jet_index << endl;
	  continue;
	}
	
	double * direction = new double[3];
	direction[0] = sin(top.JetPhi[ij]);
	direction[1] = cos(top.JetPhi[ij]);
	direction[2] = sqrt(top.JetE[ij]*top.JetE[ij] - top.JetEt[ij]*top.JetEt[ij]);
	if(top.JetEta[ij]<0) direction[2] *= -1.;
	TJet* jetp  = new TJet;
	jetp->pt  = top.JetEt[ij];
	jetp->eta = top.JetEta[ij];
	jetp->phi = top.JetPhi[ij];
	jetp->E   = top.JetE[ij];
	jetp->flavor = (TJet::Flavor)top.JetFlavor[ij];
	//the following is not quite correct, as this factor is different for all towers.
	//These values should be in the n-tupel as well
	double factor =  top.JetEt[ij] /  top.JetE[ij];
	jetp->HadF = had * factor;
	jetp->EMF = em * factor;
	jetp->OutF = out * factor;
	//Create an Top TData event
	top_data[nstoredjets] = 
	  new TData_InvMass2(jet_index * p->GetNumberOfJetParametersPerBin() + p->GetNumberOfTowerParameters(),
			     direction,                                     //p_T direction of this jet
			     massConstraint_W,                                         //truth//
			     sqrt(pow(0.5,2)+pow(0.10*top.JetPt[ij],2)),    //error//
			     top.Weight,                                    //weight//
			     //1.,                                          //weight//
			     p->GetJetParRef( jet_index ),                  //params
			     p->GetNumberOfJetParametersPerBin(),           //number of free jet param. p. bin
			     p->jet_parametrization,                        //function
			     //p->dummy_parametrization,
			     jet_error_param,                               //error param. function
			     jetp                                           //jet momentum for plotting and scale
			     );
	
	//Add the jet's towers to "top_data":
	for (int n=0; n<top.NobjTow; ++n){
	  if (top.Tow_jetidx[n]!=(int)ij) continue;//look for ij-jet's towers
	  //if (top.TowEt[n]<0.01) continue;
	  
	  int index = p->GetBin(p->GetEtaBin(top.TowId_eta[n]),
				p->GetPhiBin(top.TowId_phi[n]));
	  //std::cout << "jet:" << ij << "bin index:" << index << "\n";
	  if (index<0){ cerr<<"WARNING: Top tower_index = " << index << endl; continue; }
	  
	  double relativEt = top.TowEt[n]/top.JetEt[ij];
	  //if (relativEt<=0) cerr << "relEt = " <<relativEt << endl; //continue;
	  //This relativeE is used *only* for plotting! Therefore no cuts on this var!
	  //create array with multidimensional measurement
	  TMeasurement * mess = new TTower;
	  mess->pt = double(top.TowEt[n]);
	  double scale = 0.;
	  if (top.TowE[n]!=0.) scale = top.TowEt[n]/top.TowE[n];
	  mess->EMF = double(top.TowEm[n]*scale);
	  mess->HadF = double(top.TowHad[n]*scale);
	  mess->OutF = double(top.TowOE[n]*scale);
	  mess->eta = double(top.TowEta[n]);
	  mess->phi = double(top.TowPhi[n]);
	  mess->E = double(top.TowE[n]);
	  //mess[7] = double( cos( top.JetCalPhi-top.TowPhi[n] ) ); // Projection factor for summing tower Pt
	  
	  top_data[nstoredjets]->AddMess(new TData_TruthMess(index,
							     mess,                                                   //mess//
							     top.JetPt[ij] * relativEt,                             //truth//
							     sqrt(pow(0.5,2)+pow(0.1*top.JetPt[ij]*relativEt,2)),   //error//
							     //1., //weight//
							     top.Weight,                                            //weight//
							     p->GetTowerParRef( index ),                             //parameter//
							     p->GetNumberOfTowerParametersPerBin(),                  //number of free tower param. p. bin//
							     p->tower_parametrization,                               //function//
							     tower_error_param                                       //error param. function//
							     ));
	}
	if(nstoredjets> 0) {
          top_data[0]->AddNewMultMess( top_data[nstoredjets] );
	}
	++nstoredjets;
	goodevent = (nstoredjets==2 ? true : false);
      }//loop over all top-jets
      
      if (goodevent) {
	data.push_back( top_data[0] );
      }
    }

    if(useMassConstraintTop && (dataClass == 0)) {
      bool goodevent = false;
      TData_InvMass2 * top_data[3]; //two W-jets and one b-jet
      top_data[0] = 0;
      int nstoredjets = 0;
      for (unsigned int ij = 0; ij<3; ++ij){
	if(top.JetPt[ij] < Et_cut_on_jet) continue;
	
	//Find the jets eta & phi index using the nearest tower to jet axis:
	int jet_index=-1;
	double min_tower_dr = 10.0;
	double em = 0;
	double had = 0;
	double out = 0;
	TLorentzVector Ljet(0,0,0,0);
	
	Ljet.SetPtEtaPhiE(top.JetPt[ij],top.JetEta[ij],top.JetPhi[ij],top.JetE[ij]);
	for (int n=0; n<top.NobjTow; ++n){
	  if (top.Tow_jetidx[n]!=(int)ij) continue;//look for ij-jet's towers
	  em += top.TowEm[n];
	  had += top.TowHad[n];
	  out += top.TowOE[n];
	  TLorentzVector Ltower(0,0,0,0);
	  Ltower.SetPtEtaPhiE(top.TowEt[n],top.TowEta[n],top.TowPhi[n],top.TowE[n]);
	  double dr = Ltower.DeltaR(Ljet);
	  if (dr<min_tower_dr) {
	    jet_index = p->GetJetBin(p->GetJetEtaBin(top.TowId_eta[n]),
				     p->GetJetPhiBin(top.TowId_phi[n]));
	    min_tower_dr = dr;
	  }
	}
	if (jet_index<0){
	  cerr<<"WARNING: Top jet_index = " << jet_index << endl;
	  continue;
	}
	
	double * direction = new double[3];
	direction[0] = sin(top.JetPhi[ij]);
	direction[1] = cos(top.JetPhi[ij]);
	direction[2] = sqrt(top.JetE[ij]*top.JetE[ij] - top.JetEt[ij]*top.JetEt[ij]);
	TJet* jetp  = new TJet;
	jetp->pt  = top.JetEt[ij];
	jetp->eta = top.JetEta[ij];
	jetp->phi = top.JetPhi[ij];
	jetp->E   = top.JetE[ij];
	jetp->flavor = (TJet::Flavor)top.JetFlavor[ij];
	//the following is not quite correct, as this factor is different for all towers.
	//These values should be in the n-tupel as well
	double factor =  top.JetEt[ij] /  top.JetE[ij];
	jetp->HadF = had * factor;
	jetp->EMF = em * factor;
	jetp->OutF = out * factor;
	//Create an Top TData event
	top_data[nstoredjets] = 
	  new TData_InvMass2(jet_index * p->GetNumberOfJetParametersPerBin() + p->GetNumberOfTowerParameters(),
			     direction,                                     //p_T direction of this jet
			     massConstraint_Top,                                         //truth//
			     sqrt(pow(0.5,2)+pow(0.10*top.JetPt[ij],2)),    //error//
			     top.Weight,                                    //weight//
			     //1.,                                          //weight//
			     p->GetJetParRef( jet_index ),                  //params
			     p->GetNumberOfJetParametersPerBin(),           //number of free jet param. p. bin
			     p->jet_parametrization,                        //function
			     //p->dummy_parametrization,
			     jet_error_param,                               //error param. function
			     jetp                                           //jet momentum for plotting and scale
			     );
	
	//Add the jet's towers to "top_data":
	for (int n=0; n<top.NobjTow; ++n){
	  if (top.Tow_jetidx[n]!=(int)ij) continue;//look for ij-jet's towers
	  //if (top.TowEt[n]<0.01) continue;
	  
	  int index = p->GetBin(p->GetEtaBin(top.TowId_eta[n]),
				p->GetPhiBin(top.TowId_phi[n]));
	  //std::cout << "jet:" << ij << "bin index:" << index << "\n";
	  if (index<0){ cerr<<"WARNING: Top tower_index = " << index << endl; continue; }
	  
	  double relativEt = top.TowEt[n]/top.JetEt[ij];
	  //if (relativEt<=0) cerr << "relEt = " <<relativEt << endl; //continue;
	  //This relativeE is used *only* for plotting! Therefore no cuts on this var!
	  //create array with multidimensional measurement
	  TMeasurement * mess = new TTower;
	  mess->pt = double(top.TowEt[n]);
	  double scale = 0.;
	  if (top.TowE[n]!=0.) scale = top.TowEt[n]/top.TowE[n];
	  mess->EMF = double(top.TowEm[n]*scale);
	  mess->HadF = double(top.TowHad[n]*scale);
	  mess->OutF = double(top.TowOE[n]*scale);
	  mess->eta = double(top.TowEta[n]);
	  mess->phi = double(top.TowPhi[n]);
	  mess->E = double(top.TowE[n]);
	  //mess[7] = double( cos( top.JetCalPhi-top.TowPhi[n] ) ); // Projection factor for summing tower Pt
	  
	  top_data[nstoredjets]->AddMess(new TData_TruthMess(index,
							     mess,                                                   //mess//
							     top.JetPt[ij] * relativEt,                             //truth//
							     sqrt(pow(0.5,2)+pow(0.1*top.JetPt[ij]*relativEt,2)),   //error//
							     //1., //weight//
							     top.Weight,                                            //weight//
							     p->GetTowerParRef( index ),                             //parameter//
							     p->GetNumberOfTowerParametersPerBin(),                  //number of free tower param. p. bin//
							     p->tower_parametrization,                               //function//
							     tower_error_param                                       //error param. function//
							     ));
	}
	if(nstoredjets> 0) {
          top_data[0]->AddNewMultMess( top_data[nstoredjets] );
	}
	++nstoredjets;
	goodevent = (nstoredjets==3 ? true : false);
      }//loop over all top-jets
      
      if (goodevent) {
	data.push_back( top_data[0] );
      }
    }
    if(dataClass > 0) {
      TData* ev = createTwoJetsInvMassEvents();
      if(ev) data.push_back(ev);
      else continue;
    }
    ++evt;
    if (evt>=n_top_events && n_top_events>=0)
      break;
  }
  return evt;
}

TData* TopReader::createTwoJetsInvMassEvents()
{
  Jet *jets[2] = {0,0};
  double* terr = new double[top.NobjTow];
  for(int i = 0; i < 3; ++i) {
    if(top.JetPt[i] < Et_cut_on_jet) continue;
    if((TJet::Flavor)top.JetFlavor[i] != TJet::uds) continue;
    
    double em = 0;
    double had = 0;
    double out = 0;
    double err2 = 0;
    TMeasurement tower;
    double dR = 10;
    int closestTower = 0; 
    for(int n=0; n<top.NobjTow; ++n){
      if(top.Tow_jetidx[n] != i) continue;//look for ij-jet's towers
      
      em += top.TowEm[n];
      had +=  top.TowHad[n];
      out +=  top.TowOE[n];  
      tower.pt = top.TowEt[n];
      double scale = top.TowEt[n]/top.TowE[n];
      tower.EMF = top.TowEm[n]*scale;
      tower.HadF = top.TowHad[n]*scale;
      tower.OutF = top.TowOE[n]*scale;
      tower.eta = top.TowEta[n];
      tower.phi = top.TowPhi[n];
      tower.E = top.TowE[n];
      terr[n] = tower_error_param(&tower.pt,&tower,0); 
      if(terr[n] == 0) {
	//assume toy MC???
	terr[n] = TParameters::toy_tower_error_parametrization(&tower.pt,&tower);
      }
      terr[n] *= terr[n];
      err2 += terr[n];
      double dphi = TVector2::Phi_mpi_pi(top.JetPhi[i]-tower.phi);
      double dr = sqrt((top.JetEta[i]-tower.eta)*(top.JetEta[i]-tower.eta)+
		       dphi*dphi);     
      if(dr < dR) {
	dR = dr;
	closestTower = n;
      }
    } 
    if(had/(had + em) < 0.07) { return 0;}
    if(had/(had + em) > 0.92) { return 0;}
    double factor =  top.JetEt[i] /  top.JetE[i];
    tower.pt = top.JetEt[i];
    tower.EMF = em * factor;
    tower.HadF = had * factor;
    tower.OutF = out * factor;
    tower.eta = top.JetEta[i];
    tower.phi = top.JetPhi[i];
    tower.E   = top.JetE[i];
    double err =  jet_error_param(&tower.pt,&tower,0);
    err2 += err * err;
    Jet **jet = jets[0] ? &jets[1] : &jets[0];
    if(dataClass == 2) {
      JetWithTowers *jt = 
	new JetWithTowers(top.JetEt[i], em * factor, had * factor,
			  out * factor, top.JetE[i], top.JetEta[i],
			  top.JetPhi[i], TJet::uds, top.GenJetPt[i], 0.,
			  TJet::CorFactors(top.JetCorrL1[i],
					   top.JetCorrL2[i],
					   top.JetCorrL3[i],
					   top.JetCorrL4[i],            
					   top.JetCorrL5[i]),
			  p->jet_function(top.TowId_eta[closestTower],
					  top.TowId_phi[closestTower]),
			  jet_error_param, p->global_jet_function(), Et_cut_on_jet);
      for(int j = 0 ; j < top.NobjTow ; ++j) {
	if (top.Tow_jetidx[j]!= i) continue;//look for ij-jet's towers
	double scale = top.TowEt[j]/top.TowE[j];
	jt->addTower(top.TowEt[j],top.TowEm[j]*scale,
		     top.TowHad[j]*scale,top.TowOE[j]*scale,
		     top.TowE[j],top.TowEta[j],top.TowPhi[j],
		     p->tower_function(top.TowId_eta[i],top.TowId_phi[i]),
		     tower_error_param);
      }
      *jet = jt;
    }
    else { 
      *jet = new Jet(top.JetEt[i], em * factor, had * factor, out * factor,
		     top.JetE[i], top.JetEta[i], top.JetPhi[i],
		     TJet::uds, top.GenJetPt[i], 0.,
		     TJet::CorFactors(top.JetCorrL1[i],
				      top.JetCorrL2[i],
				      top.JetCorrL3[i],
				      top.JetCorrL4[i],            
				      top.JetCorrL5[i]),
		     p->jet_function(top.TowId_eta[closestTower],
				     top.TowId_phi[closestTower]),
		     jet_error_param, p->global_jet_function(), Et_cut_on_jet);    
    }
  }
  delete [] terr;
  if(jets[1]) {
    //    if(0.5*(jets[0]->pt+jets[1]->pt)>40.)
    return new TwoJetsInvMassEvent(jets[0],jets[1],massConstraint_W,1.0);
  } else {
    delete jets[0];
    delete jets[1];
  }
  return 0;
}
