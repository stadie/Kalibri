
//User libs
#include "ControlPlots.h"
#include "ConfigFile.h"
#include "CalibData.h"
//C++ libs
#include <iostream>
#include <cmath>
//Root libs
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TLatex.h"
#include "TCanvas.h"
#include "TPostScript.h"
#include "TStyle.h"
#include "TH1.h"
#include "TH2F.h"
#include "TF1.h"
#include "TLegend.h"
#include "TProfile.h"
#include "TROOT.h"

using namespace std;

void TControlPlots::ReadConfigFile(string file)
{
  ConfigFile config( file.c_str() );
  
  _doPlots = config.read<bool>("do plots",true); 

}

void TControlPlots::FitControlPlots()  // Fit Control Histograms w.r.t. towers
{
  TCanvas * c1 = new TCanvas("c1","",600,600);
  TPostScript ps("plots.ps",111);
  std::vector<TData*>::const_iterator data_it,it;
  TLatex latex;
  latex.SetTextSize(0.035);
  TH2F * constants = new TH2F("calib_constants","Calibration constants vs. Et and eta-bin with EMF=OUF=0",100,0.5,100.,p->GetEtaGranularity(),1,p->GetEtaGranularity());    
  double * testmess = new double[4];
  //int d = 0;
  for (int eta=0; eta<p->GetEtaGranularity();++eta){
    for (int phi=0; phi<p->GetPhiGranularity();++phi){
      int i = p->GetBin(eta,phi);
      double * val = p->GetTowerParRef(i);
      char * name = new char[100];
      TH1F * plot[4]; 
      sprintf(name, "h%d_eta%d_phi%d",i,eta+1,phi+1);
      plot[0] = new TH1F(name,";uncalibrated tower E_{T} [GeV];k-factor",100,0.5,100.);    
      sprintf(name, "htt%d_eta%d_phi%d",i,eta+1,phi+1);
      plot[1] = new TH1F(name,";uncalibrated tower E_{T} [GeV];k-factor",100,0.5,100.);    
      sprintf(name, "hgj%d_eta%d_phi%d",i,eta+1,phi+1);
      plot[2] = new TH1F(name,";uncalibrated tower E_{T} [GeV];k-factor",100,0.5,100.);    
      sprintf(name, "htc%d_eta%d_phi%d",i,eta+1,phi+1);
      plot[3] = new TH1F(name,";uncalibrated tower E_{T} [GeV];k-factor",100,0.5,100.);    
      TH1F * plot_had[4]; 
      sprintf(name, "hhad%d_eta%d_phi%d",i,eta+1,phi+1);
      plot_had[0] = new TH1F(name,";uncalibrated hadronic fraction of tower E_{T} [GeV];k-factor",100,0.5,100.);    
      sprintf(name, "hhadtt%d_eta%d_phi%d",i,eta+1,phi+1);
      plot_had[1] = new TH1F(name,";uncalibrated hadronic fraction of ttower E_{T} [GeV];k-factor",100,0.5,100.);    
      sprintf(name, "hhadgj%d_eta%d_phi%d",i,eta+1,phi+1);
      plot_had[2] = new TH1F(name,";uncalibrated hadronic fraction of ttower E_{T} [GeV];k-factor",100,0.5,100.);    
      sprintf(name, "hhadtc%d_eta%d_phi%d",i,eta+1,phi+1);
      plot_had[3] = new TH1F(name,";uncalibrated hadronic fraction of ttower E_{T} [GeV];k-factor",100,0.5,100.);    

      TH1F * k[4];
      sprintf(name, "k%d_eta%d_phi%d",i,eta+1,phi+1);
      k[0] = new TH1F(name,";uncalibrated tower E_{T} [GeV];k-factor",100,0.5,100.);    
      sprintf(name, "ktt%d_eta%d_phi%d",i,eta+1,phi+1);
      k[1] = new TH1F(name,";uncalibrated tower E_{T} [GeV];k-factor",100,0.5,100.);    
      sprintf(name, "kgj%d_eta%d_phi%d",i,eta+1,phi+1);
      k[2] = new TH1F(name,";uncalibrated tower E_{T} [GeV];k-factor",100,0.5,100.);    
      sprintf(name, "ktc%d_eta%d_phi%d",i,eta+1,phi+1);
      k[3] = new TH1F(name,";uncalibrated tower E_{T} [GeV];k-factor",100,0.5,100.);    
      TH1F * em[4];
      sprintf(name, "em%d_eta%d_phi%d",i,eta+1,phi+1);
      em[0] = new TH1F(name,"",100,0.5,100.);    
      sprintf(name, "emtt%d_eta%d_phi%d",i,eta+1,phi+1);
      em[1] = new TH1F(name,"",100,0.5,100.);    
      sprintf(name, "emgj%d_eta%d_phi%d",i,eta+1,phi+1);
      em[2] = new TH1F(name,"",100,0.5,100.);    
      sprintf(name, "emtc%d_eta%d_phi%d",i,eta+1,phi+1);
      em[3] = new TH1F(name,"",100,0.5,100.);    
      TH1F * had[4];
      sprintf(name, "had%d_eta%d_phi%d",i,eta+1,phi+1);
      had[0] = new TH1F(name,"",100,0.5,100.);    
      sprintf(name, "hadtt%d_eta%d_phi%d",i,eta+1,phi+1);
      had[1] = new TH1F(name,"",100,0.5,100.);    
      sprintf(name, "hadgj%d_eta%d_phi%d",i,eta+1,phi+1);
      had[2] = new TH1F(name,"",100,0.5,100.);    
      sprintf(name, "hadtc%d_eta%d_phi%d",i,eta+1,phi+1);
      had[3] = new TH1F(name,"",100,0.5,100.);    
      TH1F * au[4];
      sprintf(name, "au%d_eta%d_phi%d",i,eta+1,phi+1);
      au[0] = new TH1F(name,"",100,0.5,100.);    
      sprintf(name, "autt%d_eta%d_phi%d",i,eta+1,phi+1);
      au[1] = new TH1F(name,"",100,0.5,100.);    
      sprintf(name, "augj%d_eta%d_phi%d",i,eta+1,phi+1);
      au[2] = new TH1F(name,"",100,0.5,100.);    
      sprintf(name, "autc%d_eta%d_phi%d",i,eta+1,phi+1);
      au[3] = new TH1F(name,"",100,0.5,100.);    

      TH1F * khad[4];
      sprintf(name, "khad%d_eta%d_phi%d",i,eta+1,phi+1);
      khad[0] = new TH1F(name,";uncalibrated tower E_{T} [GeV];k-factor",100,0.5,100.);    
      sprintf(name, "khadtt%d_eta%d_phi%d",i,eta+1,phi+1);
      khad[1] = new TH1F(name,";uncalibrated tower E_{T} [GeV];k-factor",100,0.5,100.);    
      sprintf(name, "khadgj%d_eta%d_phi%d",i,eta+1,phi+1);
      khad[2] = new TH1F(name,";uncalibrated tower E_{T} [GeV];k-factor",100,0.5,100.);    
      sprintf(name, "khadtc%d_eta%d_phi%d",i,eta+1,phi+1);
      khad[3] = new TH1F(name,";uncalibrated tower E_{T} [GeV];k-factor",100,0.5,100.);    
      TH1F * et_vs_had[4];
      sprintf(name, "et_vs_had%d_eta%d_phi%d",i,eta+1,phi+1);
      et_vs_had[0] = new TH1F(name,"",100,0.5,100.);    
      sprintf(name, "et_vs_hadtt%d_eta%d_phi%d",i,eta+1,phi+1);
      et_vs_had[1] = new TH1F(name,"",100,0.5,100.);    
      sprintf(name, "et_vs_hadgj%d_eta%d_phi%d",i,eta+1,phi+1);
      et_vs_had[2] = new TH1F(name,"",100,0.5,100.);    
      sprintf(name, "et_vs_hadtc%d_eta%d_phi%d",i,eta+1,phi+1);
      et_vs_had[3] = new TH1F(name,"",100,0.5,100.);    
      TH1F * em_vs_had[4];
      sprintf(name, "em_vs_had%d_eta%d_phi%d",i,eta+1,phi+1);
      em_vs_had[0] = new TH1F(name,"",100,0.5,100.);    
      sprintf(name, "em_vs_hadtt%d_eta%d_phi%d",i,eta+1,phi+1);
      em_vs_had[1] = new TH1F(name,"",100,0.5,100.);    
      sprintf(name, "em_vs_hadgj%d_eta%d_phi%d",i,eta+1,phi+1);
      em_vs_had[2] = new TH1F(name,"",100,0.5,100.);    
      sprintf(name, "em_vs_hadtc%d_eta%d_phi%d",i,eta+1,phi+1);
      em_vs_had[3] = new TH1F(name,"",100,0.5,100.);    
      TH1F * au_vs_had[4];
      sprintf(name, "au_vs_had%d_eta%d_phi%d",i,eta+1,phi+1);
      au_vs_had[0] = new TH1F(name,"",100,0.5,100.);    
      sprintf(name, "au_vs_hadtt%d_eta%d_phi%d",i,eta+1,phi+1);
      au_vs_had[1] = new TH1F(name,"",100,0.5,100.);    
      sprintf(name, "au_vs_hadgj%d_eta%d_phi%d",i,eta+1,phi+1);
      au_vs_had[2] = new TH1F(name,"",100,0.5,100.);    
      sprintf(name, "au_vs_hadtc%d_eta%d_phi%d",i,eta+1,phi+1);
      au_vs_had[3] = new TH1F(name,"",100,0.5,100.);    
      TH1F * norm_vs_had[4];
      sprintf(name, "norm_vs_had%d_eta%d_phi%d",i,eta+1,phi+1);
      norm_vs_had[0] = new TH1F(name,"",100,0.5,100.);    
      sprintf(name, "norm_vs_hadtt%d_eta%d_phi%d",i,eta+1,phi+1);
      norm_vs_had[1] = new TH1F(name,"",100,0.5,100.);    
      sprintf(name, "norm_vs_hadgj%d_eta%d_phi%d",i,eta+1,phi+1);
      norm_vs_had[2] = new TH1F(name,"",100,0.5,100.);    
      sprintf(name, "norm_vs_hadtc%d_eta%d_phi%d",i,eta+1,phi+1);
      norm_vs_had[3] = new TH1F(name,"",100,0.5,100.);    

      TH1F * norm[4];
      sprintf(name, "hnorm%d",i);
      norm[0] = new TH1F(name,"",100,0.5,100.);    
      sprintf(name, "hnormtt%d",i);
      norm[1] = new TH1F(name,"",100,0.5,100.);    
      sprintf(name, "hnormgj%d",i);
      norm[2] = new TH1F(name,"",100,0.5,100.);    
      sprintf(name, "hnormtc%d",i);
      norm[3] = new TH1F(name,"",100,0.5,100.);    
      sprintf(name, "h2dgj%d_eta%d_phi%d",i,eta+1,phi+1);
      TH2F * plot2dgj = new TH2F(name,";uncalibrated tower E [GeV];k-factor",100,0.5,100.,100,0.0,5.0);    
      sprintf(name, "hdgj_weight%d_eta%d_phi%d",i,eta+1,phi+1);
      TH1F * plotgj_weight = new TH1F(name,";uncalibrated tower E [GeV];k-factor",100,0.5,100.);
      sprintf(name, "norm_weight%d_eta%d_phi%d",i,eta+1,phi+1);
      TH1F * norm_weight = new TH1F(name,";uncalibrated tower E [GeV];k-factor",100,0.5,100.);    
      sprintf(name, "h2dtt%d_eta%d_phi%d",i,eta+1,phi+1);
      TH2F * plot2dtt = new TH2F(name,";uncalibrated tower E [GeV];k-factor",100,0.5,100.,100,0.0,5.0);    
      sprintf(name, "h2dtc%d_eta%d_phi%d",i,eta+1,phi+1);
      TH2F * plot2dtc = new TH2F(name,";uncalibrated tower E [GeV];k-factor",100,0.5,100.,100,0.0,5.0);    
      TH1F * chi2[4];
      //TH1F * chi2red[4];
      sprintf(name, "hgj_chi2_%d_eta%d_phi%d",i,eta+1,phi+1);
      chi2[0] = new TH1F(name,";chi^{2};N",100,0.0,50.);    
      sprintf(name, "htt_chi2_%d_eta%d_phi%d",i,eta+1,phi+1);
      chi2[1] = new TH1F(name,";chi^{2};N",100,0.0,50.);    
      sprintf(name, "htc_chi2_%d_eta%d_phi%d",i,eta+1,phi+1);
      chi2[2] = new TH1F(name,";chi^{2};N",100,0.0,50.);    
      sprintf(name, "hjj_chi2_%d_eta%d_phi%d",i,eta+1,phi+1);
      chi2[3] = new TH1F(name,";chi^{2};N",100,0.0,50.);    
      
      data_it = data->begin();

      double mess, error;
      //double p[p->free_pars_per_bin];
      //int    index = -999, ndof=0;
      int thisIndexJet=0;
      //loop over all fit-events
      for (; data_it != data->end();++data_it){
        //if one fit event is composed of multiple towers, than loop over all
	mess=0.0; error=0.0;
	const std::vector<TData*>& data_ref = (*data_it)->GetRef();
	//double JetCorr = (*data_it)->GetParametrizedMess();
	double Jet=0.;
	for (it =data_ref.begin(); it!=data_ref.end(); ++it){
          Jet += (*it)->GetParametrizedMess();
	}
	double maxTowerET = 0.0;
	
	//first loop over towers
	for (it =data_ref.begin(); it!=data_ref.end(); ++it){
	  double * m = (*it)->GetMess();
	  int thisIndex = (*it)->GetIndex();
	  if (m[0]>maxTowerET) {
	    thisIndexJet = thisIndex;
	    maxTowerET = m[0];
	  }
	  if (thisIndex !=i)
	    continue; //tower belongs to a wrong bin

	  if (m[0]!=0.0){
	    norm[0]->Fill(m[0] );
	    em[ 0]->Fill( m[0], m[1] );
	    had[0]->Fill( m[0], m[2] );
	    au[ 0]->Fill( m[0], m[3] );
	    if ((*data_it)->GetType()==TypeGammaJet) {
  	      norm[2]->Fill(m[0] );
	      em[ 2]->Fill( m[0], m[1] );
	      had[2]->Fill( m[0], m[2] );
	      au[ 2]->Fill( m[0], m[3] );
	    }  
	    else if ((*data_it)->GetType()==TypeTrackTower) {
  	      norm[1]->Fill(m[0] );
	      em[ 1]->Fill( m[0], m[1] );
	      had[1]->Fill( m[0], m[2] );
	      au[ 1]->Fill( m[0], m[3] );
	    }  
	    else if ((*data_it)->GetType()==TypeTrackCluster) {
	      norm[3]->Fill(m[0] );
	      em[ 3]->Fill( m[0], m[1] );
	      had[3]->Fill( m[0], m[2] );
	      au[ 3]->Fill( m[0], m[3] );
	    }  
	  }
	  if (m[2]!=0.0){
	    norm_vs_had[0]->Fill( m[2] );
	    et_vs_had[  0]->Fill( m[2], m[0] );
	    em_vs_had[  0]->Fill( m[2], m[1] );
	    au_vs_had[  0]->Fill( m[2], m[3] );
	    if ((*data_it)->GetType()==TypeGammaJet) {
	      norm_vs_had[2]->Fill( m[2] );
	      et_vs_had[  2]->Fill( m[2], m[0] );
	      em_vs_had[  2]->Fill( m[2], m[1] );
	      au_vs_had[  2]->Fill( m[2], m[3] );
	    }  
	    else if ((*data_it)->GetType()==TypeTrackTower) {
	      norm_vs_had[1]->Fill( m[2] );
	      et_vs_had[  1]->Fill( m[2], m[0] );
	      em_vs_had[  1]->Fill( m[2], m[1] );
	      au_vs_had[  1]->Fill( m[2], m[3] );
	    }  
	    else if ((*data_it)->GetType()==TypeTrackCluster) {
	      norm_vs_had[3]->Fill( m[2] );
	      et_vs_had[  3]->Fill( m[2], m[0] );
	      em_vs_had[  3]->Fill( m[2], m[1] );
	      au_vs_had[  3]->Fill( m[2], m[3] );
	    }  
	  }
	}
      }
      em[0]->Divide(norm[0]);
      em[1]->Divide(norm[1]);
      em[2]->Divide(norm[2]);
      em[3]->Divide(norm[3]);
      had[0]->Divide(norm[0]);
      had[1]->Divide(norm[1]);
      had[2]->Divide(norm[2]);
      had[3]->Divide(norm[3]);
      au[0]->Divide(norm[0]);
      au[1]->Divide(norm[1]);
      au[2]->Divide(norm[2]);
      au[3]->Divide(norm[3]);

      em_vs_had[0]->Divide(norm_vs_had[0]);
      em_vs_had[1]->Divide(norm_vs_had[1]);
      em_vs_had[2]->Divide(norm_vs_had[2]);
      em_vs_had[3]->Divide(norm_vs_had[3]);
      et_vs_had[0]->Divide(norm_vs_had[0]);
      et_vs_had[1]->Divide(norm_vs_had[1]);
      et_vs_had[2]->Divide(norm_vs_had[2]);
      et_vs_had[3]->Divide(norm_vs_had[3]);
      au_vs_had[0]->Divide(norm_vs_had[0]);
      au_vs_had[1]->Divide(norm_vs_had[1]);
      au_vs_had[2]->Divide(norm_vs_had[2]);
      au_vs_had[3]->Divide(norm_vs_had[3]);

      //second loop over all fit-events
      for (data_it = data->begin(); data_it != data->end();++data_it){
        //if one fit event is composed of multiple towers, than loop over all
	mess=0.0; error=0.0;
	const std::vector<TData*>& data_ref = (*data_it)->GetRef();
	double JetCorr = (*data_it)->GetParametrizedMess();
	double Jet=0.;
	for (it =data_ref.begin(); it!=data_ref.end(); ++it){
          Jet += (*it)->GetParametrizedMess();
	}
	double maxTowerET = 0.0;
	
         //second loop over towers
	for (it =data_ref.begin(); it!=data_ref.end(); ++it){
	  double m = (*it)->GetMess()[0];
	  double mhad = (*it)->GetMess()[2];
	  double t   = (*it)->GetTruth()*(Jet/JetCorr);
          double tmp = (*it)->GetParametrizedMess();
          mess  += tmp;
	  error +=  (*it)->GetParametrizedErr(&tmp);
	  int thisIndex = (*it)->GetIndex();
	  if (m>maxTowerET) {
	    thisIndexJet = thisIndex;
	    maxTowerET = m;
	  }
	  if (thisIndex !=i)
	    continue; //tower belongs to a wrong bin

          if (mhad!=0.0){
	    testmess[0] = m;
	    testmess[1] = em[0]->GetBinContent(em[0]->GetXaxis()->FindBin(m));
	    testmess[2] = had[0]->GetBinContent(had[0]->GetXaxis()->FindBin(m));
	    testmess[3] = au[0]->GetBinContent(au[0]->GetXaxis()->FindBin(m));
	    plot_had[0]->Fill( mhad, t/m );
	    khad[0]->Fill(     mhad, p->plot_parametrization(testmess,val) );
	    if ((*data_it)->GetType()==TypeGammaJet) {
	      plot_had[2]->Fill( mhad, t/m );
	      khad[2]->Fill(     mhad, p->plot_parametrization(testmess,val) );
	    }  
	    else if ((*data_it)->GetType()==TypeTrackTower) {
	      khad[1]->Fill(     mhad, p->plot_parametrization(testmess,val) );
	      plot_had[1]->Fill( mhad, t/m );
	    }  
	    else if ((*data_it)->GetType()==TypeTrackCluster) {
  	      plot_had[3]->Fill( mhad, t/m );
	      khad[3]->Fill(     mhad, p->plot_parametrization(testmess,val) );
	    }  
	  }

	  if (m!=0.0){
	    plot[0]->Fill( m, t/m );
	    
	    testmess[0] = m;
	    testmess[1] = em[0]->GetBinContent(em[0]->GetXaxis()->FindBin(m));
	    testmess[2] = had[0]->GetBinContent(had[0]->GetXaxis()->FindBin(m));
	    testmess[3] = au[0]->GetBinContent(au[0]->GetXaxis()->FindBin(m));
	    k[0]->Fill( m, p->plot_parametrization(testmess,val) );

	    if ((*data_it)->GetType()==TypeGammaJet) {
	      plot[2]->Fill( m, t/m );
  	      plot2dgj->Fill( m, t/m );
	      double weight = t/Jet;
	      norm_weight->Fill( m, weight );
  	      plotgj_weight->Fill( m, t/m );
	      testmess[1] = em[2]->GetBinContent(em[2]->GetXaxis()->FindBin(m));
	      testmess[2] = had[2]->GetBinContent(had[2]->GetXaxis()->FindBin(m));
	      testmess[3] = au[2]->GetBinContent(au[2]->GetXaxis()->FindBin(m));
	      k[2]->Fill( m, p->plot_parametrization(testmess,val) );
	    }  
	    else if ((*data_it)->GetType()==TypeTrackTower) {
	      plot[1]->Fill( m, t/m );
  	      plot2dtt->Fill( m, t/m );
	      testmess[1] = em[1]->GetBinContent(em[1]->GetXaxis()->FindBin(m));
	      testmess[2] = had[1]->GetBinContent(had[1]->GetXaxis()->FindBin(m));
	      testmess[3] = au[1]->GetBinContent(au[1]->GetXaxis()->FindBin(m));
	      k[1]->Fill( m, p->plot_parametrization(testmess,val) );
	    }  
	    else if ((*data_it)->GetType()==TypeTrackCluster) {
	      plot[3]->Fill( m, t/m );
  	      plot2dtc->Fill( m, t/m );
	      testmess[1] = em[3]->GetBinContent(em[3]->GetXaxis()->FindBin(m));
	      testmess[2] = had[3]->GetBinContent(had[3]->GetXaxis()->FindBin(m));
	      testmess[3] = au[3]->GetBinContent(au[3]->GetXaxis()->FindBin(m));
	      k[3]->Fill( m, p->plot_parametrization(testmess,val) );
	    }  
	  }
	}


	if (thisIndexJet!=i)
	  continue; //event (jet or tower) belongs to a wrong bin

	switch ( (*data_it)->GetType()) {
	case TypeTrackTower://track-tower
	  if (mess != 0.0) chi2[0]->Fill( (*data_it)->chi2()/(*data_it)->GetWeight() ); break;
	case TypeGammaJet://gamma-jet
	  if (mess != 0.0) chi2[1]->Fill( (*data_it)->chi2()/(*data_it)->GetWeight() ); break;
	case TypeTrackCluster://track-cluster
	  if (mess != 0.0) chi2[2]->Fill( (*data_it)->chi2()/(*data_it)->GetWeight() ); 
	  break;
	case TypePtBalance://jet-jet
	  if (mess != 0.0) chi2[3]->Fill( (*data_it)->chi2()/(*data_it)->GetWeight() ); break;
	}
      }
      if (norm[0]->GetEntries()==0) continue;

      plot[0]->Divide(norm[0]);
      plot[1]->Divide(norm[1]);
      plot[2]->Divide(norm[2]);
      plot[3]->Divide(norm[3]);
      plot_had[0]->Divide(norm_vs_had[0]);
      plot_had[1]->Divide(norm_vs_had[1]);
      plot_had[2]->Divide(norm_vs_had[2]);
      plot_had[3]->Divide(norm_vs_had[3]);
      k[0]->Divide(norm[0]);
      k[1]->Divide(norm[1]);
      k[2]->Divide(norm[2]);
      k[3]->Divide(norm[3]);
      khad[0]->Divide(norm_vs_had[0]);
      khad[1]->Divide(norm_vs_had[1]);
      khad[2]->Divide(norm_vs_had[2]);
      khad[3]->Divide(norm_vs_had[3]);
      plotgj_weight->Divide(norm_weight);
      for (int j=0; j<4; ++j)//number of diff samples, i.e. gj, tt, tc, jj,...
	for (int b=0; b<100; ++b){//number of bins
	  if (norm[j]->GetBinContent(b)>0) {
	    //plot->SetBinContent(b, plot->GetBinContent(b) / norm->GetBinContent(b));
	    plot[j]->SetBinError(  b, 1./sqrt(
					      norm[j]->GetBinContent(b)     //stat
					      )*plot[j]->GetBinContent(b) );
	    plot_had[j]->SetBinError(  b, 1./sqrt(
					      norm_vs_had[j]->GetBinContent(b)     //stat
					      )*plot_had[j]->GetBinContent(b) );
	  }			      
	}

      plot[0]->SetMarkerStyle(8);
      plot[0]->SetMaximum(6.0);
      plot[0]->SetMinimum(0.0);

      plot[0]->Draw("pe");//average
      plot[3]->SetLineColor(3);
      plot[3]->SetMarkerColor(3);
      plot[3]->Draw("pe,same");//track-cluster
      plot[2]->SetLineColor(2);//gamma-jet
      plot[2]->SetMarkerColor(2);
      plot[2]->Draw("pe,same");
      plot[1]->SetLineColor(4);
      plot[1]->SetMarkerColor(4);
      plot[1]->Draw("pe,same");//track-tower

      k[0]->SetLineColor(1);
      k[0]->SetLineWidth(3);
      k[0]->Draw("l,same");//average
      k[3]->SetLineWidth(3);
      k[3]->SetLineColor(3);
      k[3]->Draw("l,same");//track-cluster
      k[2]->SetLineWidth(3);
      k[2]->SetLineColor(2);//gamma-jet
      k[2]->Draw("l,same");
      k[1]->SetLineWidth(3);
      k[1]->SetLineColor(4);
      k[1]->Draw("l,same");//track-tower

      latex.DrawLatex( 73,4.2,"Average");
      latex.DrawLatex( 73,3.8,"#color[2]{Gamma-Jet}");
      latex.DrawLatex( 73,3.4,"#color[3]{Track-Cluster}");
      latex.DrawLatex( 73,3.0,"#color[4]{Track-Tower}");
      
      //sprintf(name,"f(c_{i}) = %2.2f + %1.3f/#sqrt{E} + %1.4f/E",val[0],val[1],val[2]);
      //latex.DrawLatex( 0.3*(plot[0]->GetXaxis()->GetXmax()-plot[0]->GetXaxis()->GetXmin()),
      //		       0.4*(plot[0]->GetMaximum()-plot[0]->GetMinimum()),
      //		       name);
      c1->Draw(); 
      ps.NewPage();
    

      sprintf(name, "h_k_hadonly_chi2_%d_eta%d_phi%d",i,eta+1,phi+1);
      TH1F * khadonly = new TH1F(name,";k vs. had. E_{T} with EMF=OUF=0     [GeV];",100,0.5,100.);    
      sprintf(name, "h_k_Efrac02_chi2_%d_eta%d_phi%d",i,eta+1,phi+1);
      TH1F * kEfrac02 = new TH1F(name,";k vs. had. E_{T} with EMF=0.2 OUF=0 [GeV];",100,0.5,100.);    
      sprintf(name, "h_k_Efrac05_chi2_%d_eta%d_phi%d",i,eta+1,phi+1);
      TH1F * kEfrac05 = new TH1F(name,";k vs. had. E_{T} with EMF=0.5 OUF=0 [GeV];",100,0.5,100.);    
      testmess[3] = 0.0;
      for (int b=1; b<=100; ++b){
	testmess[0] = (double)b; // -> EMFrac = 0.0
	testmess[1] = 0.0;
	testmess[2] = (double)b;
	khadonly->SetBinContent(khadonly->GetXaxis()->FindBin(b), 
	                        p->plot_parametrization(testmess,val) );
	constants->SetBinContent(constants->GetXaxis()->FindBin(b),constants->GetYaxis()->FindBin(eta+1),
				 p->plot_parametrization(testmess,val));
	testmess[1] = (double)b*0.2; // -> EMFrac = 0.25
	testmess[2] = (double)b*0.8;
	kEfrac02->SetBinContent(kEfrac02->GetXaxis()->FindBin(b), 
	                        p->plot_parametrization(testmess,val) );       
	testmess[1] = (double)b*0.5; // -> EMFrac = 1.0
	testmess[2] = (double)b*0.5;
	kEfrac05->SetBinContent(kEfrac05->GetXaxis()->FindBin(b), 
	                        p->plot_parametrization(testmess,val) );       
      }
      khadonly->Draw("l"); 
      c1->Draw(); 
      ps.NewPage();

      kEfrac02->Draw("l"); 
      c1->Draw(); 
      ps.NewPage();

      kEfrac05->Draw("l"); 
      c1->Draw(); 
      ps.NewPage();

      plot_had[0]->SetMarkerStyle(8);
      plot_had[0]->SetMaximum(6.0);
      plot_had[0]->SetMinimum(0.0);
      plot_had[0]->Draw("pe");//average
      plot_had[3]->SetLineColor(3);
      plot_had[3]->SetMarkerColor(3);
      plot_had[3]->Draw("pe,same");//track-cluster
      plot_had[2]->SetLineColor(2);//gamma-jet
      plot_had[2]->SetMarkerColor(2);
      plot_had[2]->Draw("pe,same");
      plot_had[1]->SetLineColor(4);
      plot_had[1]->SetMarkerColor(4);
      plot_had[1]->Draw("pe,same");//track-tower

      khad[0]->SetLineColor(1);
      khad[0]->SetLineWidth(3);
      khad[0]->Draw("l,same");//average
      khad[3]->SetLineWidth(3);
      khad[3]->SetLineColor(3);
      khad[3]->Draw("l,same");//track-cluster
      khad[2]->SetLineWidth(3);
      khad[2]->SetLineColor(2);//gamma-jet
      khad[2]->Draw("l,same");
      khad[1]->SetLineWidth(3);
      khad[1]->SetLineColor(4);
      khad[1]->Draw("l,same");//track-tower

      latex.DrawLatex( 73,4.2,"Average");
      latex.DrawLatex( 73,3.8,"#color[2]{Gamma-Jet}");
      latex.DrawLatex( 73,3.4,"#color[3]{Track-Cluster}");
      latex.DrawLatex( 73,3.0,"#color[4]{Track-Tower}");
      
      //sprintf(name,"f(c_{i}) = %2.2f + %1.3f/#sqrt{E} + %1.4f/E",val[0],val[1],val[2]);
      //latex.DrawLatex( 0.3*(plot[0]->GetXaxis()->GetXmax()-plot[0]->GetXaxis()->GetXmin()),
      //		       0.4*(plot[0]->GetMaximum()-plot[0]->GetMinimum()),
      //		       name);
      c1->Draw(); 
      ps.NewPage();

      
      plot2dgj->SetLineColor(2);
      plot2dtc->SetLineColor(3);
      plot2dgj->Draw("BOX");	      
      plot2dtt->Draw("BOX,same");	      
      plot2dtc->Draw("BOX,same");	      
      latex.DrawLatex( 0.3*(plot2dgj->GetXaxis()->GetXmax()-plot2dgj->GetXaxis()->GetXmin()),
		       0.4*(plot2dgj->GetYaxis()->GetXmax()-plot2dgj->GetYaxis()->GetXmin()),
		       "#color[2]{Gamma-Jet}");
      latex.DrawLatex( 0.3*(plot2dgj->GetXaxis()->GetXmax()-plot2dgj->GetXaxis()->GetXmin()),
		       0.3*(plot2dgj->GetYaxis()->GetXmax()-plot2dgj->GetYaxis()->GetXmin()),
		       "#color[1]{TrackTower}");
      latex.DrawLatex( 0.3*(plot2dgj->GetXaxis()->GetXmax()-plot2dgj->GetXaxis()->GetXmin()),
		       0.2*(plot2dgj->GetYaxis()->GetXmax()-plot2dgj->GetYaxis()->GetXmin()),
		       "#color[3]{TrackCluster}");
      c1->Draw(); 
      ps.NewPage();
      
      //chi2[0]->SetMarkerStyle( 8 );
      //show overflow in the last bin
      double maximum = 0.0;
      for (int k=0; k<2; ++k){
        if (chi2[k]->GetMaximum()>maximum)
	  maximum = chi2[k]->GetMaximum(); 
	chi2[k]->SetBinContent(100,chi2[k]->GetBinContent(100)+chi2[k]->GetBinContent(101));
	chi2[k]->SetBinContent(1,chi2[k]->GetBinContent(0)+chi2[k]->GetBinContent(1));
      }
      chi2[0]->SetMaximum( maximum+sqrt(maximum) );
      chi2[0]->SetLineColor(1);
      chi2[1]->SetLineColor(2);
      chi2[2]->SetLineColor(3);
      chi2[0]->GetYaxis()->SetTitleOffset(1.3);
      chi2[0]->Draw("h");
      chi2[1]->Draw("h,same");
      chi2[2]->Draw("h,same");
      latex.DrawLatex( 0.3*(chi2[0]->GetXaxis()->GetXmax()-chi2[0]->GetXaxis()->GetXmin()),
		       0.4*(chi2[0]->GetMaximum()-chi2[0]->GetMinimum()),
		       "#color[2]{Gamma-Jet}");
      latex.DrawLatex( 0.3*(chi2[0]->GetXaxis()->GetXmax()-chi2[0]->GetXaxis()->GetXmin()),
		       0.3*(chi2[0]->GetMaximum()-chi2[0]->GetMinimum()),
		       "#color[1]{TrackTower}");
      latex.DrawLatex( 0.3*(chi2[0]->GetXaxis()->GetXmax()-chi2[0]->GetXaxis()->GetXmin()),
		       0.2*(chi2[0]->GetMaximum()-chi2[0]->GetMinimum()),
		       "#color[3]{TrackCluster}");
      c1->Draw(); 
      ps.NewPage();
    }
  }
  gStyle->SetPalette(1);
  constants->GetXaxis()->SetTitle("Et [GeV]");
  constants->GetYaxis()->SetTitle("Eta bin");
  constants->Draw("COLZ"); 
  c1->Draw(); 
  ps.NewPage();
  
  ps.Close();
}
  
void TControlPlots::GammaJetControlPlots()  // Gamma-Jet Control Histograms
{
  TCanvas * c1 = new TCanvas("gj1","",600,600);
  TPostScript ps("gammajet_plots_per_towerbin.ps",111);

  std::vector<TData*>::const_iterator data_it,it;
  TLatex latex;
  latex.SetTextSize(0.035);
  for (int eta=0; eta<p->GetEtaGranularity();++eta){
    for (int phi=0; phi<p->GetPhiGranularity();++phi){
      int i = p->GetBin(eta,phi);
      char * name = new char[100];
      sprintf(name, "hjes_gj%d_eta%d_phi%d",i,eta+1,phi+1);
      TH1F * plot_jes = new TH1F(name,";#sum calibrated tower E_{T} [GeV]; JES: ( E_{T}^{#gamma} / #sum E_{T}^{calib. tower})",100,0.0,400.);    
      sprintf(name, "h_gj%d_eta%d_phi%d",i,eta+1,phi+1);
      TH1F * plot = new TH1F(name,";uncalibrated jet E_{T} [GeV];average of ( E_{T}^{#gamma} / E_{T}^{uncalib. jet})",100,0.0,400.);    
      sprintf(name, "gj_fit%d",i);
      TH1F * fit  = new TH1F(name,"",100,0.0,400.);    
      sprintf(name, "gj_norm%d",i);
      TH1F * norm = new TH1F(name,"",100,0.0,400.);    
      sprintf(name, "gjjes_norm%d",i);
      TH1F * norm_jes = new TH1F(name,"",100,0.0,400.);    
      //sprintf(name, "gj_etaphi%d",i);
      //TH2F * etaphi = new TH2F(name,";eta;phi",eta_ntwr+1,-eta_ntwr/2-0.5,eta_ntwr/2+0.5, phi_ntwr, 0.5, phi_ntwr+0.5);
      int indexJet=0, ijets=0;      
      data_it = data->begin();
      //loop over all fit-events
      for (; data_it != data->end();++data_it){
        if ( (*data_it)->GetType()!=TypeGammaJet) continue;

	int indexTower=0;
	double Etmax=0, calib_tower_sum=0.0;
	double tower_sum = 0.0; //is equivalent to (*data_it)->GetMess(),
	                        //but since we need the index too, this is faster
	const std::vector<TData*>& data_ref = (*data_it)->GetRef();
	for (it =data_ref.begin(); it!=data_ref.end(); ++it){
	  double * tow_et = (*it)->GetMess();
	  if (tow_et[0]>Etmax){
	    indexTower=(*it)->GetIndex();
	    Etmax=tow_et[0];
	  }
	  tower_sum += tow_et[0];
	  calib_tower_sum += (*it)->GetParametrizedMess();
	}

	if (indexTower!=i)
	  continue; //event belongs to a wrong bin
        indexJet += (*data_it)->GetIndex();
	++ijets;

        double JetCorr = (*data_it)->GetParametrizedMess();

        fit->Fill(      tower_sum, JetCorr/tower_sum );
	plot->Fill(     tower_sum, (*data_it)->GetTruth()/tower_sum );
        norm->Fill(     tower_sum ); 
	plot_jes->Fill( calib_tower_sum, (*data_it)->GetTruth()/calib_tower_sum );
        norm_jes->Fill( calib_tower_sum ); 
	//etaphi->Fill( etaJet, phiJet );
      }
      if (norm->GetEntries()==0) continue;
      
      plot->Divide(norm);
      plot_jes->Divide(norm_jes);
      fit->Divide(norm);
      for (int b=0; b<100; ++b){
	if (norm->GetBinContent(b)>0) {
	  plot->SetBinError(  b, 1./sqrt(
					 norm->GetBinContent(b)     //stat
					 )*plot->GetBinContent(b) );
	}			      
	if (norm_jes->GetBinContent(b)>0) {
	  plot_jes->SetBinError(  b, 1./sqrt(
					     norm_jes->GetBinContent(b)     //stat
					     )*plot_jes->GetBinContent(b) );
	}			      
      }
	
      fit->SetLineColor( 2 );
      fit->SetLineWidth( 4 );
      plot->GetYaxis()->SetTitleOffset( 1.4 );
      plot->SetMarkerStyle( 8 );
      plot->SetMinimum(0.1);
      
      //c1->SetLogy(1);
      plot->Draw("pe");
      fit->Draw("h,same");
      TLatex latex;
      latex.SetTextSize(0.035);
      latex.DrawLatex( 0.3*(plot->GetXaxis()->GetXmax()-plot->GetXaxis()->GetXmin()),
		       0.6*(plot->GetMaximum()-plot->GetMinimum()),
		       "#color[2]{--  tower and jet corrections}");
      c1->Draw(); 
      ps.NewPage();

      plot_jes->Draw("pe");
      TF1 * res2 = new TF1("res2",p->jes_plot_parametrization, 0.5, 400., 3);
      i = indexJet/ijets - p->GetNumberOfTowerParameters();
      double * val = p->GetJetParRef(i);
      res2->SetParameters(val[0],val[1]);
      res2->SetLineWidth( 3 );
      res2->SetLineColor( 2 );
      res2->Draw("same");
      c1->Draw(); 
      ps.NewPage();
    }
  }
  ps.Close();
}

void TControlPlots::GammaJetControlPlotsJetBin()  // Gamma-Jet Control Histograms
{
  TCanvas * c1 = new TCanvas("gj2","",600,600);
  TPostScript ps("gammajet_plots_per_jetbin.ps",111);

  std::vector<TData*>::const_iterator data_it,it;
  TLatex latex;
  latex.SetTextSize(0.035);
  //one plot per *JET* bin!
  for (int eta=0; eta<p->GetEtaGranularityJet();++eta){
    for (int phi=0; phi<p->GetPhiGranularityJet();++phi){
      int i = p->GetJetBin(eta,phi) + p->GetNumberOfTowerParameters();
      char * name = new char[100];
      sprintf(name, "h2jes_gj%d_eta%d_phi%d",i,eta+1,phi+1);
      TH1F * plot_jes = new TH1F(name,";#sum calibrated tower E_{T} [GeV]; JES: ( E_{T}^{#gamma} / #sum E_{T}^{calib. tower})",100,0.0,400.);    
      sprintf(name, "h2_gj%d_eta%d_phi%d",i,eta+1,phi+1);
      TH1F * plot = new TH1F(name,";uncalibrated jet E_{T} [GeV];average of ( E_{T}^{#gamma} / E_{T}^{uncalib. jet})",100,0.0,400.);    
      sprintf(name, "gj2_fit%d",i);
      TH1F * fit  = new TH1F(name,"",100,0.0,400.);    
      sprintf(name, "gj2_norm%d",i);
      TH1F * norm = new TH1F(name,"",100,0.0,400.);    
      sprintf(name, "gj2jes_norm%d",i);
      TH1F * norm_jes = new TH1F(name,"",100,0.0,400.);    
      data_it = data->begin();
      //loop over all fit-events
      for (; data_it != data->end();++data_it){
        if ( (*data_it)->GetType()!=TypeGammaJet) continue;
	if ( (*data_it)->GetIndex()!=i)	  continue; //event belongs to a wrong bin
	
	double calib_tower_sum=0.0;
        double JetCorr = (*data_it)->GetParametrizedMess();

	double tower_sum = 0.0; //is equivalent to (*data_it)->GetMess(),
	const std::vector<TData*>& data_ref = (*data_it)->GetRef();
	for (it =data_ref.begin(); it!=data_ref.end(); ++it){
	  tower_sum += (*it)->GetMess()[0];
	  calib_tower_sum += (*it)->GetParametrizedMess();
	}


        fit->Fill(      tower_sum, JetCorr/tower_sum );
	plot->Fill(     tower_sum, (*data_it)->GetTruth()/tower_sum );
        norm->Fill(     tower_sum ); 
	plot_jes->Fill( calib_tower_sum, (*data_it)->GetTruth()/calib_tower_sum );
        norm_jes->Fill( calib_tower_sum ); 
      }
      if (norm->GetEntries()==0) continue;
      
      plot->Divide(norm);
      plot_jes->Divide(norm_jes);
      fit->Divide(norm);
      for (int b=0; b<100; ++b){
	if (norm->GetBinContent(b)>0) {
	  plot->SetBinError(  b, 1./sqrt(
					 norm->GetBinContent(b)     //stat
					 )*plot->GetBinContent(b) );
	}			      
	if (norm_jes->GetBinContent(b)>0) {
	  plot_jes->SetBinError(  b, 1./sqrt(
					     norm_jes->GetBinContent(b)     //stat
					     )*plot_jes->GetBinContent(b) );
	}			      
      }
	
      fit->SetLineColor( 2 );
      fit->SetLineWidth( 4 );
      plot->GetYaxis()->SetTitleOffset( 1.4 );
      plot->SetMarkerStyle( 8 );
      plot->SetMinimum(0.1);
      
      //c1->SetLogy(1);
      plot->Draw("pe");
      fit->Draw("h,same");
      TLatex latex;
      latex.SetTextSize(0.035);
      latex.DrawLatex( 0.3*(plot->GetXaxis()->GetXmax()-plot->GetXaxis()->GetXmin()),
		       0.6*(plot->GetMaximum()-plot->GetMinimum()),
		       "#color[2]{--  tower and jet corrections}");
      c1->Draw(); 
      ps.NewPage();

      plot_jes->Draw("pe");
      TF1 * res2 = new TF1("res2",p->jes_plot_parametrization, 0.5, 400., 3);
      i = p->GetJetBin(eta, phi);
      double * val = p->GetJetParRef(i);
      res2->SetParameters(val[0],val[1]);
      res2->SetLineWidth( 3 );
      res2->SetLineColor( 2 );
      res2->Draw("same");
      c1->Draw(); 
      ps.NewPage();
    }
  }
  ps.Close();
}

  
void TControlPlots::TrackTowerControlPlots()  // Track-Tower Control Histograms
{
  TCanvas * c1 = new TCanvas("c1","",600,600);
  TPostScript ps("tracktower_plots.ps",111);
  std::vector<TData*>::const_iterator data_it,it;
  TLatex latex;
  latex.SetTextSize(0.035);
//int d = 0;
  for (int eta=0; eta<p->GetEtaGranularity();++eta){
    for (int phi=0; phi<p->GetPhiGranularity();++phi){
      int i = p->GetBin(eta,phi);
      char * name = new char[100];
      TH1F * plot_tt, * norm_tt;
      sprintf(name, "h_tt%d_eta%d_phi%d",i,eta+1,phi+1);
      plot_tt = new TH1F(name,";uncalibrated tower E [GeV];k-factor",100,0.5,100.);    
      sprintf(name, "norm_tt%d",i);
      norm_tt = new TH1F(name,"",100,0.5,100.);
      sprintf(name, "h2d_tt%d_eta%d_phi%d",i,eta+1,phi+1);
      TH2F * plot2d_tt = new TH2F(name,";uncalibrated tower E [GeV];k-factor",100,0.5,100.,100,0.0,5.0);    
      data_it = data->begin();

      double mess, error;
      //double p[p->free_pars_per_bin];
      //int    index = -999, ndof=0;
      int thisIndexJet;
      //loop over all fit-events
      for (; data_it != data->end();++data_it){
        //if one fit event is composed of multiple towers, than loop over all
	mess=0.0; error=0.0;
	const std::vector<TData*>& data_ref = (*data_it)->GetRef();
	double JetCorr = (*data_it)->GetParametrizedMess();
	double Jet=0.;
	for (it =data_ref.begin(); it!=data_ref.end(); ++it){
          Jet += (*it)->GetParametrizedMess();
	}
	double maxTowerET = 0.0;
	for (it =data_ref.begin(); it!=data_ref.end(); ++it){
	  double m   = (*it)->GetMess()[0];
	  double t   = (*it)->GetTruth()*(Jet/JetCorr);
          double tmp = (*it)->GetParametrizedMess();
          mess  += tmp;
	  error +=  (*it)->GetParametrizedErr(&tmp);
	  int thisIndex = (*it)->GetIndex();
	  if (m>maxTowerET) {
	    thisIndexJet = thisIndex;
	    maxTowerET = m;
	  }
	  
	  if (thisIndex !=i)
	    continue; //tower belongs to a wrong bin
	  if (m!=0.0){
	    if ((*data_it)->GetType()==TypeTrackTower) {
	      plot_tt->Fill( m, t/m );
	      norm_tt->Fill( m );
  	      plot2d_tt->Fill( m, t/m );
	    } 
	  }
	}
	if (thisIndexJet!=i)
	  continue; //event (jet or tower) belongs to a wrong bin

      }
	  
      if (norm_tt->GetEntries()==0) continue;

      plot_tt->Divide(norm_tt);
      for (int b=0; b<100; ++b){
	if (norm_tt->GetBinContent(b)>0) {
	  //plot->SetBinContent(b, plot->GetBinContent(b) / norm->GetBinContent(b));
	  plot_tt->SetBinError(  b, 1./sqrt(
				      //pow(plot2c->GetBinContent(b)/plot2norm->GetBinContent(b),2)  //syst
				      //+  
				      norm_tt->GetBinContent(b)     //stat
				      )*plot_tt->GetBinContent(b) );
	}
      }

      TF1 * res1_tt = new TF1("res1",p->plot_parametrization, 0.5, 100., 3);
      double * val = p->GetTowerParRef(i);
      res1_tt->SetParameters(val[0],val[1]);
      res1_tt->SetLineWidth( 3 );
      res1_tt->SetLineColor( 2 );
      plot_tt->SetMarkerStyle(8);
      plot_tt->SetMarkerColor(4);
      plot_tt->SetLineColor(4);
      plot_tt->Draw("l");//track-tower
      res1_tt->Draw("same"); 
      
      //Check if fit values are equal to the start values, in that case use red text color
//      int r=1;
//      vector<double>::const_iterator st = start_values.begin();
//      for (CalibVal::const_iterator it = val.begin(); it!=val.end(), 
//           st!=start_values.end(); ++it, ++st)
//        r *= ((double)(*it)==(*st));
//      if (r==0)
	sprintf(name,"f(c_{i}) = %2.2f + %1.3f/#sqrt{E} + %1.4f/E",val[0],val[1],val[2]);
//      else
//	sprintf(name,"#color[2]{f(c_{i}) = %2.2f + %1.3f/#sqrt{E} + %1.4f/E}",val[0],val[1],val[2]);
      latex.DrawLatex( 0.3*(plot_tt->GetXaxis()->GetXmax()-plot_tt->GetXaxis()->GetXmin()),
		       0.4*(plot_tt->GetMaximum()-plot_tt->GetMinimum()),
		       name);
      c1->Draw(); 
      ps.NewPage();
      
      plot2d_tt->Draw("BOX");	      
      latex.DrawLatex( 0.3*(plot2d_tt->GetXaxis()->GetXmax()-plot2d_tt->GetXaxis()->GetXmin()),
		       0.3*(plot2d_tt->GetYaxis()->GetXmax()-plot2d_tt->GetYaxis()->GetXmin()),
		       "#color[1]{TrackTower}");
      c1->Draw(); 
      ps.NewPage();
      
    }
  }
  
  ps.Close();
  /*
  TCanvas * c1 = new TCanvas("trkttow1","",600,600);
  TPostScript ps("tracktower_plots.ps",111);

  //ps.Range(xsize,ysize);
  //c1->Clear();
  std::vector<TData*>::const_iterator it;
  gStyle->SetOptStat(0);
  const static unsigned ptbins = 20;
  const static double   maxpt  = 400;
  int events = 0;
  for (int eta=0; eta<p->GetEtaGranularity();++eta){
    for (int phi=0; phi<p->GetPhiGranularity();++phi){
      int i = eta*p->GetPhiGranularity()*p->free_pars_per_bin+phi*p->free_pars_per_bin;
      char * name = new char[100];
      TH1F * plota[ ptbins ];
      TH1F * plotb[ ptbins ];
      //sprintf(name, "gj_etaphi%d eta:%d phi:%d;eta;phi",i,eta,phi);
      //TH2F * etaphi = new TH2F(name,name,eta_ntwr,0,eta_ntwr, phi_ntwr1, 0, phi_ntwr1);
      for (unsigned p=0; p<ptbins; ++p){
	sprintf(name, "h%da_trkt%d_eta%d_phi%d",p,i,eta+1,phi+1);
	plota[p] = new TH1F(name,";(E^{tower}-E^{track})/E^{track};events",100,-20.0,20.);    
	sprintf(name, "h%db_trkt%d_eta%d_phi%d",p,i,eta+1,phi+1);
	plotb[p] = new TH1F(name,"",100,-20.0,20.);    
      }
      it = data->begin();
      //loop over all fit-events
      for (; it != data->end();++it){
        if ( (*it)->GetType()!=1) continue;
        //if one fit event is composed of multiple towers, than loop ovver all
	int thisIndex=(*it)->GetID();
	if (thisIndex!=i || (*it)->GetMess()==0.0)
	  continue; //event belongs to a wrong bin
	++events;

        CalibVal val = k.find((*it)->GetID())->second;
	double p[] = {val[0],val[1],val[2]};
	double tow_et = (*it)->GetMess();
	parametrization( &tow_et, p);
        for (unsigned k=0; k<ptbins; ++k){
	  if ( k*maxpt/ptbins < (*it)->GetTruth() && (*it)->GetTruth() <= (k+1)*maxpt/ptbins) {
	    plota[k]->Fill( (parametrization(&tow_et,p) - (*it)->GetTruth())/(*it)->GetTruth() );
	    plotb[k]->Fill( (tow_et - (*it)->GetTruth())/(*it)->GetTruth() );
	  }
	}
	//etaphi->Fill((*it)->GetID().first, (*it)->GetID().second);
      }
      //if (norm->GetEntries()==0) continue;
      double maximum = 0.0;
      TLegend leg(0.6,0.7,0.89,0.89);
      leg.SetFillColor(0);
      for (unsigned k=0; k<ptbins; ++k){
        if (plota[k]->GetEntries()<10) continue;
        //show overflow in the last bin
	plota[k]->SetBinContent(100,plota[k]->GetBinContent(100)+plota[k]->GetBinContent(101));
	plotb[k]->SetBinContent(100,plotb[k]->GetBinContent(100)+plotb[k]->GetBinContent(101));
        //show underflow in the first bin
	plota[k]->SetBinContent(1,plota[k]->GetBinContent(0)+plota[k]->GetBinContent(1));
	plotb[k]->SetBinContent(1,plotb[k]->GetBinContent(0)+plotb[k]->GetBinContent(1));
	
        if (plota[k]->GetMaximum()>maximum) maximum=plota[k]->GetMaximum();
        if (plotb[k]->GetMaximum()>maximum) maximum=plotb[k]->GetMaximum();
	plota[k]->GetYaxis()->SetTitleOffset( 1.4 );
	plota[k]->SetMinimum(0.1);
	plota[k]->SetLineColor( 2+k );
	plota[k]->SetLineWidth( 4 );
	plota[k]->SetLineStyle( 1 );
	plotb[k]->SetLineColor( 2+k );
	plotb[k]->SetLineWidth( 4 );
	plotb[k]->SetLineStyle( 2 );
	sprintf(name, "%d < E < %d GeV; calibrated", (int)(k*maxpt/ptbins),(int)((k+1)*maxpt/ptbins) );
        leg.AddEntry(plota[k],name, "l");
	sprintf(name, "%d < E < %d GeV; un-calib.", (int)(k*maxpt/ptbins),(int)((k+1)*maxpt/ptbins) );
        leg.AddEntry(plotb[k],name, "l");
      }
      plota[0]->SetMaximum(maximum+sqrt(maximum));
      plota[0]->Draw("h");
      plotb[0]->Draw("ha,same");
      for (unsigned k=1; k<ptbins; ++k){
        if (plota[k]->GetEntries()<10) continue;
	plota[k]->Draw("ha,same");
	plotb[k]->Draw("ha,same");
      }
      leg.Draw("same");
      //c1->SetLogy(1);
      c1->Draw(); 
      ps.NewPage();

      //etaphi->Draw("BOX");
      //c1->Draw(); 
      //ps.NewPage();
    }
  }
  cout << "...used " << events << " track-tower events." << endl;
  gStyle->SetOptStat(1);  
  ps.Close();
  */
}

void TControlPlots::TrackClusterControlPlots()  // Track-Cluster Control Histograms
{
  TCanvas * c1 = new TCanvas("tc1","",600,600);
  TPostScript ps("trackcluster_plots.ps",111);

  std::vector<TData*>::const_iterator data_it,it;
  TLatex latex;
  latex.SetTextSize(0.035);
  for (int eta=0; eta<p->GetEtaGranularity();++eta){
    for (int phi=0; phi<p->GetPhiGranularity();++phi){
      int i = p->GetBin(eta,phi);
      char * name = new char[100];
      sprintf(name, "hjes_tc%d_eta%d_phi%d",i,eta+1,phi+1);
//      TH1F * plot_jes = new TH1F(name,";#sum calibrated tower E_{T} [GeV]; JES: ( E_{T}^{#gamma} / #sum E_{T}^{calib. tower})",100,0.0,400.);    
      sprintf(name, "h_tc%d_eta%d_phi%d",i,eta+1,phi+1);
      TH1F * plot = new TH1F(name,";uncalibrated jet E_{T} [GeV];average of ( E_{T}^{#gamma} / E_{T}^{uncalib. jet})",100,0.0,400.);    
      sprintf(name, "tc_fit%d",i);
      TH1F * fit  = new TH1F(name,"",100,0.0,400.);    
      sprintf(name, "tc_norm%d",i);
      TH1F * norm = new TH1F(name,"",100,0.0,400.);    
      data_it = data->begin();
      //loop over all fit-events
      for (; data_it != data->end();++data_it){
        if ( (*data_it)->GetType()!=TypeTrackCluster) continue;

	int indexJet=0;
	double Etmax=0, calib_tower_sum=0.0;
	double tower_sum = 0.0; //is equivalent to (*data_it)->GetMess(),
	                        //but since we need the index too, this is faster
	const  std::vector<TData*>& data_ref = (*data_it)->GetRef();
	for (it =data_ref.begin(); it!=data_ref.end(); ++it){
	  double tow_et = (*it)->GetMess()[0];
	  if (tow_et>Etmax){
	    indexJet=(*it)->GetIndex();
	    Etmax=tow_et;
	  }
	  tower_sum += tow_et;
	  calib_tower_sum += (*it)->GetParametrizedMess();
	}

	if (indexJet!=i)
	  continue; //event belongs to a wrong bin

        double JetCorr = (*data_it)->GetParametrizedMess();

        fit->Fill(      tower_sum, JetCorr/tower_sum );
	plot->Fill(     tower_sum, (*data_it)->GetTruth()/tower_sum );
        norm->Fill(     tower_sum ); 
	//etaphi->Fill( etaJet, phiJet );
      }
      if (norm->GetEntries()==0) continue;
      
      plot->Divide(norm);
      fit->Divide(norm);
      for (int b=0; b<100; ++b){
	if (norm->GetBinContent(b)>0) {
	  plot->SetBinError(  b, 1./sqrt(
				      norm->GetBinContent(b)     //stat
				      )*plot->GetBinContent(b) );
	}			      
      }
	
      fit->SetLineColor( 2 );
      fit->SetLineWidth( 4 );
      plot->GetYaxis()->SetTitleOffset( 1.4 );
      plot->SetMarkerStyle( 8 );
      plot->SetMinimum(0.1);
      
      //c1->SetLogy(1);
      plot->Draw("pe");
      fit->Draw("h,same");
      TLatex latex;
      latex.SetTextSize(0.035);
      latex.DrawLatex( 0.5*(plot->GetXaxis()->GetXmax()-plot->GetXaxis()->GetXmin()),
		       0.6*(plot->GetMaximum()-plot->GetMinimum()),
		       "#color[2]{--  tower corrections}");
      c1->Draw(); 
      ps.NewPage();

    }
  }
  ps.Close();
}


void TControlPlots::GammaJetControlPlotsJetJEC()
{
  TCanvas * c1 = new TCanvas("controlplots","",600,600);
  TPostScript ps("controlplots.ps",111);

  //book hists
  TH2F* heta[12];
  heta[0] = new TH2F("heta","#gamma-jet;#eta",100,-5,5,100,0,4);
  for(int i = 1 ; i < 12 ; ++i) heta[i] = (TH2F*)heta[0]->Clone();
  heta[3]->SetTitle("#gamma-jet 10 < E_{T}^{#gamma} < 35 GeV;#eta");
  heta[6]->SetTitle("#gamma-jet 35 < E_{T}^{#gamma} < 90 GeV;#eta");
  heta[9]->SetTitle("#gamma-jet 90 < E_{T}^{#gamma} < 300 GeV;#eta");

  TH2F* hpt[3];
  hpt[0] = new TH2F("hpt","#gamma-jet;p_{T} [GeV]",100,20,220,100,0,4);
  hpt[1] = (TH2F*)hpt[0]->Clone();
  hpt[2] = (TH2F*)hpt[0]->Clone();
  
  TH2F* hemf[3];
  hemf[0] = new TH2F("hemf","#gamma-jet;EMF",100,0,1,100,0,4);
  hemf[1] = (TH2F*)hemf[0]->Clone();
  hemf[2] = (TH2F*)hemf[0]->Clone();
 
  TH1F* hptGamma = new TH1F("hptGamma","#gamma-jet",100,20,220);
  hptGamma->SetXTitle("p_{T} [GeV]");
  TH1F* hptGammaW = new TH1F("hptGammaW","#gamma-jet",100,20,220);
  hptGammaW->SetXTitle("p_{T} [GeV]");
  

  double bins[101];
  for(int i = 0; i < 101 ; ++i) {
    bins[i] = pow(10,(i+32)/40.0);
  }
  TH2F* hptlog[3];
  hptlog[0] = new TH2F("hptlog","#gamma-jet;p_{T} [GeV]",100,bins,100,0,4);
  hptlog[1] = (TH2F*)hptlog[0]->Clone();
  hptlog[2] = (TH2F*)hptlog[0]->Clone();
  
  //loop over all fit-events
  for ( std::vector<TData*>::iterator i = data->begin() ; i != data->end() ; ++i )  {
    TData* jg = *i;
    if( jg->GetType() != TypeGammaJet ) continue;
    double etjet = jg->GetMess()[0];
    double etjetcor = jg->GetParametrizedMess();
    double etajet = jg->GetMess()[1];
    //double phijet = jg->GetMess()[2];
    heta[0]->Fill(etajet,etjet/ jg->GetTruth(),jg->GetWeight());
    heta[1]->Fill(etajet,etjetcor/ jg->GetTruth(),jg->GetWeight());
    heta[2]->Fill(etajet,etjet/etjetcor,jg->GetWeight());
    if (jg->GetTruth() > 10 && jg->GetTruth() < 35) {
      heta[3]->Fill(etajet,etjet/ jg->GetTruth(),jg->GetWeight());
      heta[4]->Fill(etajet,etjetcor/ jg->GetTruth(),jg->GetWeight());
      heta[5]->Fill(etajet,etjet/etjetcor,jg->GetWeight());
    } else if (jg->GetTruth() > 35 && jg->GetTruth() < 90) {
      heta[6]->Fill(etajet,etjet/ jg->GetTruth(),jg->GetWeight());
      heta[7]->Fill(etajet,etjetcor/ jg->GetTruth(),jg->GetWeight());
      heta[8]->Fill(etajet,etjet/etjetcor,jg->GetWeight());
    } else if (jg->GetTruth() > 90 && jg->GetTruth() < 300) {
      heta[9]->Fill(etajet,etjet/ jg->GetTruth(),jg->GetWeight());
      heta[10]->Fill(etajet,etjetcor/ jg->GetTruth(),jg->GetWeight());
      heta[11]->Fill(etajet,etjet/etjetcor,jg->GetWeight());
    } 
    hpt[0]->Fill(jg->GetTruth(),etjet/ jg->GetTruth(),jg->GetWeight());
    hpt[1]->Fill(jg->GetTruth(),etjetcor/jg->GetTruth(),jg->GetWeight());
    hpt[2]->Fill(etjetcor,etjet/etjetcor,jg->GetWeight());    
    hptlog[0]->Fill(jg->GetTruth(),etjet/ jg->GetTruth(),jg->GetWeight());
    hptlog[1]->Fill(jg->GetTruth(),etjetcor/jg->GetTruth(),jg->GetWeight());
    hptlog[2]->Fill(etjetcor,etjet/etjetcor,jg->GetWeight());
    hptGamma->Fill(jg->GetTruth());
    hptGammaW->Fill(jg->GetTruth(),jg->GetWeight());
    //em fraction plots     
    double em = 0;
    double had = 0;
    for(std::vector<TData*>::const_iterator t = jg->GetRef().begin(); t != jg->GetRef().end(); ++t) {
      TData* tt = *t;
      em  += tt->GetMess()[1];
      had += tt->GetMess()[2];
      had += tt->GetMess()[3];
    }
    hemf[0]->Fill(em/(em+had),etjet/jg->GetTruth(),jg->GetWeight());
    hemf[1]->Fill(em/(em+had),etjetcor/jg->GetTruth(),jg->GetWeight());
    hemf[2]->Fill(em/(em+had),etjet/etjetcor,jg->GetWeight());
  } 
  TH1F* hists[12][6];
  TLegend* leg = new TLegend(0.7,0.96,0.96,0.72);
  leg->AddEntry(heta[0],"p^{jet}_{T}/ E_{T}^{#gamma}","p");
  leg->AddEntry(heta[2],"p_{T}^{jet}/p_{T}^{cor. jet}","p");
  leg->AddEntry(heta[1],"p_{T}^{cor. jet}/E_{T}^{#gamma}","p");
  for(int i = 0 ; i < 12 ; i+= 3) {
    heta[i]->SetMarkerStyle(20);
    heta[i]->SetMarkerColor(1);
    heta[i]->SetMinimum(0.5);
    heta[i]->SetMaximum(1.2);
    heta[i+2]->SetMarkerStyle(4);
    heta[i+2]->SetMarkerColor(4);  
    heta[i+1]->SetMarkerStyle(22);
    heta[i+1]->SetMarkerColor(2);
    heta[i]->SetMinimum(0.5);
    heta[i]->SetMaximum(1.2);
    leg->Draw();
    Fit2D(heta[i],hists[i]);
    Fit2D(heta[i+1],hists[i+1]);
    Fit2D(heta[i+2],hists[i+2]);
    for(int j = 0 ; j < 6 ; ++j) {
      hists[i][j]->Draw();
      hists[i][j]->SetStats(0);
      hists[i+1][j]->Draw("SAME");
      hists[i+2][j]->Draw("SAME");
      leg->Draw();
      c1->SetGrid();
      c1->Draw();   
      ps.NewPage(); 
    }
  } 
  TLegend* leg2 = new TLegend(0.7,0.96,0.96,0.72);
  leg2->AddEntry(hists[3][2],"10 < E_{T}^{#gamma} < 35 GeV","p");
  leg2->AddEntry(hists[6][2],"35 < E_{T}^{#gamma} < 90 GeV","p");
  leg2->AddEntry(hists[9][2],"90 < E_{T}^{#gamma} < 300 GeV","p");
  hists[3][2]->SetTitle("p^{jet}_{T}/ E_{T}^{#gamma}   mean;#eta");
  hists[3][3]->SetTitle("p^{jet}_{T}/ E_{T}^{#gamma}   rel. width;#eta");
  hists[4][2]->SetTitle("p^{cor. jet}_{T}/ E_{T}^{#gamma}   mean;#eta");
  hists[4][3]->SetTitle("p^{cor. jet}_{T}/ E_{T}^{#gamma}    rel. width;#eta");
  for(int j = 3 ; j < 5 ; ++j) {
    for(int i = 2 ; i < 4 ; ++i) {
      hists[j][i]->SetMarkerStyle(20);
      hists[j][i]->SetMarkerColor(1);
      hists[j+3][i]->SetMarkerStyle(4);
      hists[j+3][i]->SetMarkerColor(4);  
      hists[j+6][i]->SetMarkerStyle(22);
      hists[j+6][i]->SetMarkerColor(2);
      hists[j][i]->Draw();
      hists[j][i]->SetStats(0);
      hists[j+3][i]->Draw("SAME");
      hists[j+6][i]->Draw("SAME");
      leg2->Draw();
      c1->SetGrid();
      c1->Draw();   
      ps.NewPage();
    }
  }
  delete leg2;
  for(int i = 0 ; i < 12 ; ++i) {
    for(int j = 0 ; j < 6 ; ++j) {
      delete hists[i][j];
    }	
  }
  hpt[0]->SetMarkerStyle(20);
  hpt[0]->SetMarkerColor(1);
  hpt[2]->SetMarkerStyle(4);
  hpt[2]->SetMarkerColor(4);  
  hpt[1]->SetMarkerStyle(22);
  hpt[1]->SetMarkerColor(2);
  Fit2D(hpt[0],hists[0]);
  Fit2D(hpt[1],hists[1]);
  Fit2D(hpt[2],hists[2]);
  for(int i = 0 ; i < 6 ; ++i) {
    hists[0][i]->Draw();
    hists[0][i]->SetStats(0);
    hists[1][i]->Draw("SAME");
    hists[2][i]->Draw("SAME");
    leg->Draw();
    c1->SetGrid();
    c1->Draw();   
    ps.NewPage(); 
  }
  for(int i = 0 ; i < 6 ; ++i) {
    delete hists[0][i];
    delete hists[1][i];
    delete hists[2][i];
  }
  hemf[0]->SetMarkerStyle(20);
  hemf[0]->SetMarkerColor(1);
  hemf[2]->SetMarkerStyle(4);
  hemf[2]->SetMarkerColor(4);  
  hemf[1]->SetMarkerStyle(22);
  hemf[1]->SetMarkerColor(2);
  Fit2D(hemf[0],hists[0]);
  Fit2D(hemf[1],hists[1]);
  Fit2D(hemf[2],hists[2]);
  for(int i = 0 ; i < 6 ; ++i) {
    hists[0][i]->Draw();
    hists[0][i]->SetStats(0);
    hists[1][i]->Draw("SAME");
    hists[2][i]->Draw("SAME");
    leg->Draw();
    c1->SetGrid();
    c1->Draw();   
    ps.NewPage(); 
  }
  for(int i = 0 ; i < 6 ; ++i) {
    delete hists[0][i];
    delete hists[1][i];
    delete hists[2][i];
  }
  
  hptlog[0]->SetMarkerStyle(20);
  hptlog[0]->SetMarkerColor(1);
  hptlog[0]->SetMinimum(0.2);
  hptlog[0]->SetMaximum(1.8);
  hptlog[2]->SetMarkerStyle(4);
  hptlog[2]->SetMarkerColor(4);  
  hptlog[1]->SetMarkerStyle(22);
  hptlog[1]->SetMarkerColor(2);  
  Fit2D(hptlog[0],hists[0]);
  Fit2D(hptlog[1],hists[1]);
  Fit2D(hptlog[2],hists[2]);
  for(int i = 0 ; i < 6 ; ++i) {
    hists[0][i]->Draw();
    hists[0][i]->SetStats(0);
    hists[1][i]->Draw("SAME");
    hists[2][i]->Draw("SAME");
    c1->SetLogx(1);
    c1->SetGrid();
    leg->Draw();   
    c1->Draw();   
    ps.NewPage(); 
  }
  for(int i = 0 ; i < 6 ; ++i) {
    delete hists[0][i];
    delete hists[1][i];
    delete hists[2][i];
  }
  ps.NewPage(); 
  delete leg;
  hptGamma->SetMarkerStyle(20);
  hptGamma->SetMarkerColor(1);
  hptGammaW->SetMarkerStyle(22);
  hptGammaW->SetMarkerColor(2);
  hptGammaW->Draw("p");
  hptGammaW->SetStats(0);
  hptGamma->Draw("pSAME");
  c1->SetLogx(0);  
  c1->SetLogy(1);  
  c1->SetGrid();
  leg = new TLegend(0.7,0.96,0.96,0.72);
  leg->AddEntry(hptGamma,"no weights","p");
  leg->AddEntry(hptGammaW,"with weights","p");
  leg->Draw();
  c1->Draw(); 
  ps.NewPage(); 
  delete leg;
  for(int i = 0 ; i < 12 ; ++i) delete heta[i];
  delete hpt[0];
  delete hpt[1];
  delete hpt[2];
  delete hptlog[0];
  delete hptlog[1];
  delete hptlog[2];
  delete hptGamma;
  delete hptGammaW;  
  c1->SetLogx(0);  
  c1->SetLogy(0);   
  c1->SetGrid(0);
  //tower plots
  TH1F* htow[3];
  htow[0] = new TH1F("htow0","tower response;id_{#eta}",82,0,82);
  htow[1] = (TH1F*)htow[0]->Clone();
  htow[2] = (TH1F*)htow[0]->Clone();
  for (int eta=-41; eta < 41;++eta){
    int i = p->GetEtaBin(eta >= 0 ? eta +1 : eta);
    double x[4] = {0,0,0,0};
    double* par = p->GetTowerParRef(p->GetBin(i,0));
    htow[0]->Fill(eta + 41,TParameters::tower_parametrization(x,par));
    x[2] = 5;
    htow[1]->Fill(eta + 41,TParameters::tower_parametrization(x,par));
    x[2] = 50;
    htow[2]->Fill(eta + 41,TParameters::tower_parametrization(x,par));
  }
  htow[0]->SetMinimum(-5);
  htow[0]->SetMaximum(80);
  htow[0]->SetStats(0);  
  htow[0]->SetMarkerStyle(20);
  htow[0]->SetMarkerColor(1);
  htow[1]->SetMarkerStyle(22);
  htow[1]->SetMarkerColor(2);
  htow[2]->SetMarkerStyle(24);
  htow[2]->SetMarkerColor(4);
  htow[0]->Draw("p");
  htow[1]->Draw("p SAME");
  htow[2]->Draw("p SAME"); 
  leg = new TLegend(0.4,0.96,0.66,0.72);
  leg->AddEntry(htow[0]," 0 GeV","p");
  leg->AddEntry(htow[1]," 5 GeV","p");
  leg->AddEntry(htow[2],"50 GeV","p");
  leg->Draw();
  c1->Draw(); 
  ps.Close();
  delete htow[0];
  delete htow[1];
  delete htow[2];
  delete leg;
}

void TControlPlots::Fit2D(TH2F* hist, TH1F* hresults[6]) 
{
  //book hists
  TString s(hist->GetName());
  s.Append("_res");
  if( hist->GetXaxis()->GetXbins()->GetSize() == hist->GetNbinsX() +1) {
    hresults[0] = new TH1F(s,hist->GetTitle(),hist->GetNbinsX(),hist->GetXaxis()->GetXbins()->GetArray());
  } else {
    hresults[0] = new TH1F(s,hist->GetTitle(),hist->GetNbinsX(),hist->GetXaxis()->GetXmin(),
			   hist->GetXaxis()->GetXmax());
  }
  hresults[0]->SetXTitle(hist->GetXaxis()->GetTitle());
  hresults[0]->SetMarkerStyle(hist->GetMarkerStyle());
  hresults[0]->SetMarkerColor(hist->GetMarkerColor());
  hresults[0]->SetMarkerSize(hist->GetMarkerSize());
  for(int i = 1; i < 6 ; ++i) {
    hresults[i] = (TH1F*)hresults[0]->Clone();
  }
  s = hist->GetTitle();
  hresults[0]->SetTitle(s.Append(" mean")); 
  s = hist->GetTitle();
  hresults[1]->SetTitle(s.Append(" standard deviation")); 
  s = hist->GetTitle();
  hresults[2]->SetTitle(s.Append(" mean of Gauss fit")); 
  s = hist->GetTitle();
  hresults[3]->SetTitle(s.Append(" width of Gauss fit"));
  s = hist->GetTitle();
  hresults[4]->SetTitle(s.Append(" median"));    \
  for(int i = 0 ; i < 6 ; ++i) {
    hresults[i]->SetMinimum(0.2);
    hresults[i]->SetMaximum(1.8);
    ++i;
    hresults[i]->SetMinimum(0.0);
    hresults[i]->SetMaximum(0.3);
  }
  TH1F* htemp = new TH1F("htemp","",hist->GetNbinsY(),hist->GetYaxis()->GetXmin(),
		          hist->GetYaxis()->GetXmax());
  htemp->Sumw2();
  const int nq = 2;
  double yq[2],xq[2];
  xq[0] = 0.5;
  xq[1] = 0.90;
  for(int i = 1 ; i <= hist->GetNbinsX() ; ++i) {
    htemp->Reset();
    for(int j = 1 ; j <= hist->GetNbinsY() ; ++j) {
      htemp->Fill(htemp->GetBinCenter(j),hist->GetBinContent(hist->GetBin(i,j)));
    }  
    if(htemp->GetSumOfWeights() <= 0) continue;
    htemp->Fit("gaus","LLQNO","");
    TF1 *f = (TF1*)gROOT->GetFunction("gaus")->Clone();
    double mean = f->GetParameter(1);
    double meanerror = f->GetParError(1);
    double width = f->GetParameter(2);
    if(width < 0.2) width = 0.2;
    if( (htemp->Fit(f,"LLQNO","goff",mean - 2 * width, mean + 2 * width) == 0) && (f->GetProb() > 0.01)) {
      mean = f->GetParameter(1);
      meanerror = f->GetParError(1);
      width = f->GetParameter(2);
      hresults[2]->SetBinContent(i,mean);
      hresults[2]->SetBinError(i,meanerror);
      hresults[3]->SetBinContent(i,width/mean);
      hresults[3]->SetBinError(i, f->GetParError(2));
    }
    mean = htemp->GetMean();
    meanerror = htemp->GetMeanError();
    width = htemp->GetRMS();
    hresults[0]->SetBinContent(i,mean);
    hresults[0]->SetBinError(i,meanerror);
    hresults[1]->SetBinContent(i,width/mean); 
    hresults[1]->SetBinError(i,htemp->GetRMSError());
    htemp->GetQuantiles(nq,yq,xq);
    hresults[4]->SetBinContent(i,yq[0]);
    hresults[4]->SetBinError(i,0.0001);
    hresults[5]->SetBinContent(i,yq[1]/yq[0]-1);
    hresults[5]->SetBinError(i,0.0001);
    delete f;
  }
  delete htemp;
}
