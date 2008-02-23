
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
      TH1F * chi2red[4];
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
      int thisIndexJet;
      //loop over all fit-events
      for (; data_it != data->end();++data_it){
        //if one fit event is composed of multiple towers, than loop over all
	mess=0.0; error=0.0;
	std::vector<TData*> data_ref = (*data_it)->GetRef();
	double JetCorr = (*data_it)->GetParametrizedMess();
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
	std::vector<TData*> data_ref = (*data_it)->GetRef();
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
	  if (mess != 0.0) chi2[0]->Fill( (*data_it)->chi2() ); break;
	case TypeGammaJet://gamma-jet
	  if (mess != 0.0) chi2[1]->Fill( (*data_it)->chi2() ); break;
	case TypeTrackCluster://track-cluster
	  if (mess != 0.0) chi2[2]->Fill( (*data_it)->chi2() ); 
	  break;
	case TypePtBalance://jet-jet
	  if (mess != 0.0) chi2[3]->Fill( (*data_it)->chi2() ); break;
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

	int indexTower;
	double Etmax=0, calib_tower_sum=0.0;
	double tower_sum = 0.0; //is equivalent to (*data_it)->GetMess(),
	                        //but since we need the index too, this is faster
	std::vector<TData*> data_ref = (*data_it)->GetRef();
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
	std::vector<TData*> data_ref = (*data_it)->GetRef();
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
	std::vector<TData*> data_ref = (*data_it)->GetRef();
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
      int r=1;
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
      TH1F * plot_jes = new TH1F(name,";#sum calibrated tower E_{T} [GeV]; JES: ( E_{T}^{#gamma} / #sum E_{T}^{calib. tower})",100,0.0,400.);    
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

	int indexJet;
	double Etmax=0, calib_tower_sum=0.0;
	double tower_sum = 0.0; //is equivalent to (*data_it)->GetMess(),
	                        //but since we need the index too, this is faster
	std::vector<TData*> data_ref = (*data_it)->GetRef();
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
  TCanvas * c1 = new TCanvas("gj3","",600,600);
  TPostScript ps("gammajet_plots_ala_JEC.ps",111);

  //book hists
  TProfile* heta = new TProfile("heta","#gamma-jet",100,-5,5,0.2,1.8);
  heta->SetXTitle("#eta");
  //heta->SetYTitle("< #frac{p_{T,jet}}{E_{T,#gamma}}>");
  TProfile* hetacor = new TProfile("hetacor","#eta dependence",100,-5,5,0.2,1.8);
  hetacor->SetXTitle("#eta_{jet}^{cal}");
  //hetacor->SetYTitle("< #frac{p_{T,jet}}{E_{T,#gamma}}>");
  TProfile* hetapar = new TProfile("hetapar","#eta dependence",100,-5,5,0.2,1.8);
  hetapar->SetXTitle("#eta_{jet}^{cal}");
  //hetapar->SetYTitle("< #frac{p_{T,jet}}{E_{T,cor. jet}}>");

  TProfile* hpt = new TProfile("hpt","#gamma-jet",100,20,220,0.2,1.8);
  hpt->SetXTitle("p_{T} [GeV]");
  //hpt->SetYTitle("< #frac{p_{T,jet}}{E_{T,#gamma}}>");
  TProfile* hptcor = new TProfile("hptcor","#p_{T} dependence",100,20,220,0.2,1.8);
  hptcor->SetXTitle("p_{T}");
  //hptcor->SetYTitle("< #frac{p_{T,jet}}{E_{T,#gamma}}>");
  TProfile* hptpar = new TProfile("hptpar","p_{T} dependence",100,20,220,0.2,1.8);
  hptpar->SetXTitle("p_{T}");
  //hptpar->SetYTitle("< #frac{p_{T,jet}}{E_{T,cor. jet}}>");
 
  double bins[101];
  for(int i = 0; i < 101 ; ++i) {
    bins[i] = pow(10,(i+32)/40.0);
  }
  TProfile* hptlog = new TProfile("hptlog","#gamma-jet",100,bins,0,4);
  hptlog->SetXTitle("p_{T} [GeV]");
  //hptlog->SetYTitle("< #frac{p_{T,jet}}{E_{T,#gamma}}>");
  TProfile* hptlogcor = new TProfile("hptlogcor","#p_{T} dependence",100,bins,0,4);
  hptlogcor->SetXTitle("p_{T}");
  //hptlogcor->SetYTitle("< #frac{p_{T,jet}}{E_{T,#gamma}}>");
  TProfile* hptlogpar = new TProfile("hptlogpar","p_{T} dependence",100,bins,0,4);
  hptlogpar->SetXTitle("p_{T}");
  //hptlogpar->SetYTitle("< #frac{p_{T,jet}}{E_{T,cor. jet}}>");
  
  //loop over all fit-events
  for ( std::vector<TData*>::iterator i = data->begin() ; i != data->end() ; ++i )  {
    TData* jg = *i;
    if( jg->GetType() != TypeGammaJet ) continue;
    double etjet = jg->GetMess()[0];
    double etjetcor = jg->GetParametrizedMess();
    double etajet = jg->GetMess()[1];
    //double phijet = jg->GetMess()[2];
    heta->Fill(etajet,etjet/ jg->GetTruth());
    hetacor->Fill(etajet,etjetcor/ jg->GetTruth());
    hetapar->Fill(etajet,etjet/etjetcor);
    hpt->Fill(jg->GetTruth(),etjet/ jg->GetTruth());
    hptcor->Fill(jg->GetTruth(),etjetcor/jg->GetTruth());
    hptpar->Fill(etjetcor,etjet/etjetcor);    
    hptlog->Fill(jg->GetTruth(),etjet/ jg->GetTruth());
    hptlogcor->Fill(jg->GetTruth(),etjetcor/jg->GetTruth());
    hptlogpar->Fill(etjetcor,etjet/etjetcor);
  }
  heta->SetMarkerStyle(20);
  heta->SetMarkerColor(1);
  heta->SetMinimum(0.5);
  heta->SetMaximum(1.2);
  hetapar->SetMarkerStyle(4);
  hetapar->SetMarkerColor(4);  
  hetacor->SetMarkerStyle(22);
  hetacor->SetMarkerColor(2);
  heta->Draw();
  heta->SetStats(0);
  hetacor->Draw("SAME");
  hetapar->Draw("SAME");
  TLegend* leg = new TLegend(0.7,0.96,0.96,0.72);
  leg->AddEntry(heta,"< p^{jet}_{T}/ E_{T}^{#gamma}>","p");
  leg->AddEntry(hetapar,"<p_{T}^{jet}/p_{T}^{cor. jet}>","p");
  leg->AddEntry(hetacor,"<p_{T}^{cor. jet}/E_{T}^{#gamma}>","p");
  leg->Draw();
  c1->SetGrid();
  c1->Draw(); 
  ps.NewPage();
  delete leg;
  hpt->SetMarkerStyle(20);
  hpt->SetMarkerColor(1);
  hpt->SetMinimum(0.2);
  hpt->SetMaximum(1.8);
  hptpar->SetMarkerStyle(4);
  hptpar->SetMarkerColor(4);  
  hptcor->SetMarkerStyle(22);
  hptcor->SetMarkerColor(2);
  hpt->Draw();
  hpt->SetStats(0);
  hptcor->Draw("SAME");
  hptpar->Draw("SAME"); 
  leg = new TLegend(0.7,0.96,0.96,0.72);
  leg->AddEntry(hpt,"< p^{jet}_{T}/ E_{T}^{#gamma}>","p");
  leg->AddEntry(hptpar,"<p_{T}^{jet}/p_{T}^{cor. jet}>","p");
  leg->AddEntry(hptcor,"<p_{T}^{cor. jet}/E_{T}^{#gamma}>","p");
  leg->Draw();
  c1->SetGrid();
  c1->Draw();   
  ps.NewPage(); 
  delete leg;
  hptlog->SetMarkerStyle(20);
  hptlog->SetMarkerColor(1);
  hptlog->SetMinimum(0.2);
  hptlog->SetMaximum(1.8);
  hptlogpar->SetMarkerStyle(4);
  hptlogpar->SetMarkerColor(4);  
  hptlogcor->SetMarkerStyle(22);
  hptlogcor->SetMarkerColor(2);
  hptlog->Draw();
  hptlog->SetStats(0);
  hptlogcor->Draw("SAME");
  hptlogpar->Draw("SAME");
  c1->SetLogx(1);  
  c1->SetGrid();
  leg = new TLegend(0.7,0.96,0.96,0.72);
  leg->AddEntry(hptlog,"< p^{jet}_{T}/ E_{T}^{#gamma}>","p");
  leg->AddEntry(hptlogpar,"<p_{T}^{jet}/p_{T}^{cor. jet}>","p");
  leg->AddEntry(hptlogcor,"<p_{T}^{cor. jet}/E_{T}^{#gamma}>","p");
  leg->Draw();
  c1->Draw(); 
  ps.Close();
  delete leg;
  delete heta;
  delete hetacor;
  delete hetapar;
  delete hpt;
  delete hptcor;
  delete hptpar;
  delete hptlog;
  delete hptlogcor;
  delete hptlogpar;
}
