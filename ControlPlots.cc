
//User libs
#include "ControlPlots.h"
#include "ConfigFile.h"
#include "CalibData.h"
#include "CalibMath.h"
//C++ libs
#include <iostream>
#include <cmath>
//Root libs
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TLatex.h"
#include "TLine.h"
#include "TCanvas.h"
#include "TPostScript.h"
#include "TStyle.h"
#include "TH1.h"
#include "TH2F.h"
//#include "TF1.h"
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



  /////////////// OUTLIER CONTROL PLOTS //////////////////

  // Distribution of chi2 summands (normalized residuals)
  TH1F *h_chi2 = new TH1F("h_chi2","Scaled residuals z^{2} = 1/w #chi^{2};f(z^{2});dN / df(z^{2})",50,0,200);
  h_chi2->SetLineColor(1);
  h_chi2->SetLineStyle(3);

  // Distribution of Cauchy scaled normalized residuals
  TH1F *h_cauchy = static_cast<TH1F*>(h_chi2->Clone("h_cauchy"));
  h_cauchy->SetLineColor(2);
  h_cauchy->SetLineStyle(2);

  // Distribution of Huber scaled normalized residuals
  TH1F *h_huber = static_cast<TH1F*>(h_chi2->Clone("h_huber"));
  h_huber->SetLineColor(4);
  h_huber->SetLineStyle(1);

  // Cauchy-scaled versus no scaling
  TH2F *h_none_cauchy = new TH2F("h_none_cauchy","Residuals z^{2};z^{2};f(z^{2})",50,0,10,50,0,10);
  h_none_cauchy->SetLineColor(2);

  // Huber-scaled versus no scaling
  TH2F *h_none_huber = static_cast<TH2F*>(h_none_cauchy->Clone("h_none_huber"));
  h_none_huber->SetLineColor(4);

  for(  std::vector<TData*>::const_iterator it = data->begin();  it < data->end();  ++it )
    {
      double weight = (*it)->GetWeight();

      TData::ScaleResidual = &TData::ScaleNone;
      double res = ( (*it)->chi2() ) / weight;
      TData::ScaleResidual = &TData::ScaleCauchy;
      double res_cauchy = ( (*it)->chi2() ) / weight;
      TData::ScaleResidual = &TData::ScaleHuber;
      double res_huber = ( (*it)->chi2() ) / weight;

      h_chi2->Fill(res);
      h_cauchy->Fill(res_cauchy);
      h_huber->Fill(res_huber);
      h_none_cauchy->Fill(res,res_cauchy);
      h_none_huber->Fill(res,res_huber);
    }

  h_cauchy->Draw();
  h_huber->Draw("same");
  h_chi2->Draw("same");
  c1->SetLogy(1);

  TLegend *l_res = new TLegend(0.35,0.68,0.7,0.88);
  l_res->SetFillColor(0);
  l_res->SetBorderSize(0);
  l_res->SetHeader("Scaling function f");
  l_res->AddEntry(h_chi2,"None","L");
  l_res->AddEntry(h_huber,"Huber","L");
  l_res->AddEntry(h_cauchy,"Cauchy","L");
  l_res->Draw("same");

  c1->Draw();
  ps.NewPage();


  h_none_cauchy->Draw("box");
  h_none_huber->Draw("boxsame");
  c1->SetLogy(0);

  TLine *line1 = new TLine(0,0,7,7);
  line1->SetLineStyle(2);
  line1->SetLineColor(1);
  line1->SetLineWidth(1);
  line1->Draw("same");

  TLegend *l_res2 = new TLegend(0.35,0.68,0.7,0.88);
  l_res2->SetFillColor(0);
  l_res2->SetBorderSize(0);
  l_res2->SetHeader("Scaling function f");
  l_res2->AddEntry(line1,"None","L");
  l_res2->AddEntry(h_none_huber,"Huber","L");
  l_res2->AddEntry(h_none_cauchy,"Cauchy","L");
  l_res2->Draw("same");

  c1->Draw();
  ps.NewPage();

  delete h_chi2;
  delete h_cauchy;
  delete h_huber;
  delete l_res;
  delete h_none_cauchy;
  delete h_none_huber;
  delete line1;
  delete l_res2;

  /////////////// END: OUTLIER CONTROL PLOTS //////////////////

  
  ps.Close();
}
  
void TControlPlots::GammaJetControlPlots() // Gamma-Jet Control Histograms
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
      int i = p->GetJetBin(eta,phi)*p->GetNumberOfJetParametersPerBin() + p->GetNumberOfTowerParameters();
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
	if ( (*data_it)->GetIndex()!= i )	  continue; //event belongs to a wrong bin
	
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
      res2->SetParameters(val[0],val[1],val[2]);
      res2->SetLineWidth( 3 );
      res2->SetLineColor( 2 );
      res2->Draw("same");
      c1->Draw(); 
      ps.NewPage();
    }
  }
  ps.Close();
}

  
void TControlPlots::TrackTowerControlPlots() // Track-Tower Control Histograms
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

void TControlPlots::TrackClusterControlPlots() // Track-Cluster Control Histograms
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
  TCanvas * c2 = new TCanvas("controlplotsGauss","",600,600);
  c2->Divide(1,3);
  TPostScript ps("controlplotsGammaJet.ps",111);
  //int phi_granularity = p->GetPhiGranularity();

  //book hists
  TH1I* towerinjet[4];
  towerinjet[0]= new TH1I("No. Tower in Jet","No. Tower in Jet",40,0,40);
  for(int i =1;i<4;++i)   towerinjet[i] = (TH1I*) towerinjet[0]->Clone();
  towerinjet[1]->SetTitle("#gamma-jet 10 < E_{T}^{#gamma} < 35 GeV;No. Tower in Jet");
  towerinjet[2]->SetTitle("#gamma-jet 35 < E_{T}^{#gamma} < 90 GeV;No. Tower in Jet");
  towerinjet[3]->SetTitle("#gamma-jet 90 < E_{T}^{#gamma} < 300 GeV;No. Tower in Jet");
  TH2F* respvstet[4];
  respvstet[0] = new TH2F("Response vs. Tower Et","Response vs. Tower Et;raw tower Et;Correction factor",100,0,100,100,-4,10);
  for(int i = 1 ; i < 4 ; ++i) respvstet[i] = (TH2F*)respvstet[0]->Clone();
  respvstet[1]->SetTitle("Response vs. Tower Et (leading Tower);raw tower Et;Correction factor");
  respvstet[2]->SetTitle("Response vs. Tower Et (neighbouring leading Tower);raw tower Et;Correction factor");
  respvstet[3]->SetTitle("Response vs. Tower Et (2 towers next to leading Tower);raw tower Et;Correction factor");


  TH2F* testNegResponse[4];
  testNegResponse[0] = new TH2F("Response vs. Hadronic Fraction","Response vs. Hadronic Fraction; Had/EM+Had; Response",120,-0.1,1.1,100,-4,12);
  testNegResponse[1] = new TH2F("Response vs. Hadronic Fraction eta>3","Response vs. Hadronic Fraction #eta > 3; Had/EM+Had; Response",120,-0.1,1.1,100,-4,12);
  testNegResponse[2] = new TH2F("Response vs. Eta(jet)","Response vs. #eta(Jet);#eta(jet); Respone",100,-5,5,100,-4,12);
  testNegResponse[3] = new TH2F("Response vs. Eta_index","Response vs. #eta indx;#eta index; Respone",45,0,45,100,-4,12);
 
  int bining_plot = 9;
  TH1F* leadToNext[21];  //0+1 rings,2-12 eta dependence ,13+14 Relrings, 15-20 Pt depend (not rel)
  leadToNext[0] = new  TH1F("Et in ring / raw jetEt","Et in ring (raw = red) / raw jet Et;#Delta R;Ring Et / Jet Et",bining_plot,0,0.6);
  leadToNext[1] = (TH1F*) leadToNext[0]->Clone();   //raw tower
  leadToNext[2] = new  TH1F("Sum*(ET in ring)/Sum(ET in ring(raw))","Sum(ET in ring) / Sum(ET in ring(raw));#Delta R;Ring Et / Ring Et(raw)",bining_plot,0,0.6);
  for(int i=3;i<21;++i)
    {
      leadToNext[i] = (TH1F*) leadToNext[2]->Clone();
    }
  leadToNext[3]->SetTitle("Sum*(ET in ring)/Sum(ET in ring(raw)) abs#eta < 0.5");
  leadToNext[4]->SetTitle("Sum*(ET in ring)/Sum(ET in ring(raw)) 0.5 < abs#eta < 1.0");
  leadToNext[5]->SetTitle("Sum*(ET in ring)/Sum(ET in ring(raw)) 1.0 < abs#eta < 1.5");
  leadToNext[6]->SetTitle("Sum*(ET in ring)/Sum(ET in ring(raw)) 1.5 < abs#eta < 2.0");
  leadToNext[7]->SetTitle("Sum*(ET in ring)/Sum(ET in ring(raw)) 2.0 < abs#eta < 2.5");
  leadToNext[8]->SetTitle("Sum*(ET in ring)/Sum(ET in ring(raw)) 2.5 < abs#eta < 3.0");
  leadToNext[9]->SetTitle("Sum*(ET in ring)/Sum(ET in ring(raw)) 3.0 < abs#eta < 3.5");
  leadToNext[10]->SetTitle("Sum*(ET in ring)/Sum(ET in ring(raw)) 3.5 < abs#eta < 4.0");
  leadToNext[11]->SetTitle("Sum*(ET in ring)/Sum(ET in ring(raw)) 4.0 < abs#eta < 4.5");
  leadToNext[12]->SetTitle("Sum*(ET in ring)/Sum(ET in ring(raw)) 4.5 < abs#eta < 5.0");

  leadToNext[13]->SetTitle("Rel Et (divided by leading tower Et) in ring (raw = red);#Delta R;Rel Et in Ring");
  leadToNext[15]->SetTitle("ET in ring / raw Jet Et   10 < P_{T}^{Jet} < 35 GeV (raw = red);#Delta R;Ring Et / raw jet Et");
  leadToNext[17]->SetTitle("ET in ring / raw Jet Et   35 < P_{T}^{Jet} < 90 GeV (raw = red);#Delta R;Ring Et / raw jet Et");
  leadToNext[19]->SetTitle("ET in ring / raw Jet Et   90 < P_{T}^{Jet} < 300 GeV (raw = red);#Delta R;Ring Et / raw jet Et");

  TH2F* EtaPhiMap = new TH2F("Eta-Phi hit map","#eta-#Phi hit map;#eta;#Phi",200,-5,5,128,-3.2,3.2);


  TH2F* heta[12];
  heta[0] = new TH2F("heta","#gamma-jet;#eta",100,-5,5,100,0,4);
  for(int i = 1 ; i < 12 ; ++i) heta[i] = (TH2F*)heta[0]->Clone();
  heta[3]->SetTitle("#gamma-jet 10 < E_{T}^{#gamma} < 35 GeV;#eta");
  heta[6]->SetTitle("#gamma-jet 35 < E_{T}^{#gamma} < 90 GeV;#eta");
  heta[9]->SetTitle("#gamma-jet 90 < E_{T}^{#gamma} < 300 GeV;#eta");

  TH2F* hpt[3];
  hpt[0] = new TH2F("hpt","#gamma-jet;truth p_{T} [GeV]",400,0,400,100,0,4);
  hpt[1] = (TH2F*)hpt[0]->Clone();
  hpt[2] = (TH2F*)hpt[0]->Clone();

  TH2F* henergy[3];
  henergy[0] = new TH2F("energy","#gamma-jet;Energy [GeV]",600,0,600,100,0,4);
  henergy[1] = (TH2F*)henergy[0]->Clone();
  henergy[2] = (TH2F*)henergy[0]->Clone();
  
  TH2F* hpt_uncorr[12];
  hpt_uncorr[0] = new TH2F("hpt_uncorr","#gamma-jet;uncorr. jet p_{T} [GeV]",400,0,400,100,0,4);
  for(int i = 1 ; i < 12 ; ++i)  hpt_uncorr[i] = (TH2F*)hpt_uncorr[0]->Clone();
  hpt_uncorr[3]->SetTitle("#gamma-jet |#eta| < 1.4;uncorr. jet p_{T} [GeV]");
  hpt_uncorr[6]->SetTitle("#gamma-jet 1.4 < |#eta| < 3.0;uncorr. jet p_{T} [GeV]");
  hpt_uncorr[9]->SetTitle("#gamma-jet 3.0 < |#eta|;uncorr. jet p_{T} [GeV]");
  
  TH2F* hemf[3];
  hemf[0] = new TH2F("hemf","#gamma-jet;EMF",100,0,1,100,0,4);
  hemf[1] = (TH2F*)hemf[0]->Clone();
  hemf[2] = (TH2F*)hemf[0]->Clone();





 
  TH1F* hptGamma[3];
  TH1F* hptGammaW[3];
  hptGamma[0]  = new TH1F("hptGamma","#gamma-jet;p_{T} [GeV]",100,0,400);
  hptGammaW[0] = new TH1F("hptGammaW","#gamma-jet;p_{T} [GeV]",100,0,400);
  hptGamma[1]  = new TH1F("hetaGamma","#gamma-jet;#eta",100,-5,5);
  hptGammaW[1] = new TH1F("hetaGammaW","#gamma-jet;#eta",100,-5,5);
  hptGamma[2]  = new TH1F("hemfGamma","#gamma-jet;EMF",100,0,1);
  hptGammaW[2] = new TH1F("hemfGammaW","#gamma-jet;EMF",100,0,1);
  TH2F* hptGamma2D = new TH2F("hGamma2D","#gamma-jet;p_{T};EMF",500,0,500,100,0,1);
  TH2F* hptGamma2DW = new TH2F("hGamma2DW","#gamma-jet weights;p_{T};EMF",500,0,500,100,0,1);

  double bins[101];
  for(int i = 0; i < 101 ; ++i) {
    bins[i] = pow(10,(i+32)/40.0);
  }
  TH2F* hptlog[3];
  hptlog[0] = new TH2F("hptlog","#gamma-jet;p_{T} [GeV]",100,bins,100,0,4);
  hptlog[1] = (TH2F*)hptlog[0]->Clone();
  hptlog[2] = (TH2F*)hptlog[0]->Clone();


  for(int i=0; i<3;++i)
    {
	  heta[4*i]->Sumw2();
	  heta[4*i+1]->Sumw2();
	  heta[4*i+2]->Sumw2();
	  heta[4*i+3]->Sumw2();
	  hpt_uncorr[4*i]->Sumw2();
	  hpt_uncorr[4*i+1]->Sumw2();
	  hpt_uncorr[4*i+2]->Sumw2();
	  hpt_uncorr[4*i+3]->Sumw2();
	  henergy[i]->Sumw2();
	  hpt[i]->Sumw2();
	  hemf[i]->Sumw2();
	  hptlog[i]->Sumw2();
    }

  double ringsSum[11][bining_plot];
  double ringsRawSum[11][bining_plot];
    for(int b=0; b < 11; ++b)
      {
	for(int a=0; a < bining_plot; ++a)
	  {
	    ringsSum[b][a]=0;
	    ringsRawSum[b][a]=0;
	  }
      }
  
  //loop over all fit-events
  for ( std::vector<TData*>::iterator i = data->begin() ; i != data->end() ; ++i )  {
    TData* jg = *i;
    
    if( jg->GetType() != TypeGammaJet ) continue;
    double etjet = jg->GetMess()[0];
    double energyjet = jg->GetMess()[0];  //Et -> E see below
    double etjetcor = jg->GetParametrizedMess();
    double etajet = jg->GetMess()[1];
    double phijet = jg->GetMess()[2];
    int noTower=0;
    double maxTowerET=0;
    double maxTowerEnergy=0;
    double maxTowerETraw=0;
    int thisIndexJet=0;
    double maxTowerEta=0;
    double maxTowerPhi=0;
    const std::vector<TData*>& data_ref =(*i)->GetRef();
    for(std::vector<TData*>::const_iterator it = data_ref.begin();it != data_ref.end(); ++it)
      {
	double * m = (*it)->GetMess();
	double  pm = (*it)->GetParametrizedMess();
	respvstet[0]->Fill(m[0],pm/m[0]);
	int thisIndex = (*it)->GetIndex();
	if (pm>maxTowerET) {
	  thisIndexJet = thisIndex;
	  maxTowerET = pm;
	  maxTowerETraw = m[0];
	  maxTowerEta = m[4];
	  maxTowerPhi = m[5];
	  maxTowerEnergy = m[6];
	  }
	testNegResponse[0]->Fill(m[2]/(m[1]+m[2]),pm/m[0]);
	if(fabs(etajet) > 3)	testNegResponse[1]->Fill(m[2]/(m[1]+m[2]),pm/m[0]);
	testNegResponse[2]->Fill(etajet,pm/m[0]);
	testNegResponse[3]->Fill(thisIndex,pm/m[0]);

	EtaPhiMap->Fill(m[4],m[5]);

	++noTower;
      }
    energyjet *= maxTowerEnergy / maxTowerETraw;
    respvstet[1]->Fill(maxTowerETraw,maxTowerET/maxTowerETraw);
    double rings[bining_plot];
    double ringsRaw[bining_plot];
    for(int a=0; a < bining_plot; ++a)
      {
	rings[a]=0;
	ringsRaw[a]=0;
      }
    //second loop
    for(std::vector<TData*>::const_iterator it = data_ref.begin();it != data_ref.end(); ++it)
      {
	//next to leading towers:
	double * m  = (*it)->GetMess();
	double  pm = (*it)->GetParametrizedMess();
	if(m[0] != 0)
	  {
	    for(int a=0; a < bining_plot; ++a)
	      {
		double DeltaR = sqrt((deltaPhi(phijet,m[5]) * deltaPhi(phijet,m[5])) + ((etajet - m[4]) * (etajet - m[4])));
		if(DeltaR <= ((0.6 / bining_plot) * a))
		  {
		    rings[a] += pm / etjet;
		    ringsRaw[a] += m[0] / etjet;
		  }
	      }	    
	    
	    int index = (*it)->GetIndex();
	    if((abs(index - thisIndexJet) == 1)      )// || (abs(index - thisIndexJet + phi_granularity) < 2)  || (abs(index - thisIndexJet - phi_granularity) < 2))
	      respvstet[2]->Fill(m[0], (*it)->GetParametrizedMess()/m[0]);
	    
	    if(abs(index - thisIndexJet) == 2)  // different with proper Phi granularity
	      respvstet[3]->Fill(m[0], (*it)->GetParametrizedMess()/m[0]);
	  }
      }

    for(int a=0; a < bining_plot; ++a)
      {
	if(a == 0)
	  {
	    leadToNext[0]->Fill(a * (0.6 / bining_plot),rings[0]);
	    leadToNext[1]->Fill(a * (0.6 / bining_plot),ringsRaw[0]);
	    leadToNext[13]->Fill(a * (0.6 / bining_plot),1);
	    leadToNext[14]->Fill(a * (0.6 / bining_plot),1);
	    if((etjet > 10) && (etjet < 35))
	      {
		leadToNext[15]->Fill(a * (0.6 / bining_plot),rings[0]);
		leadToNext[16]->Fill(a * (0.6 / bining_plot),ringsRaw[0]);
	      }
	    if((etjet > 35) && (etjet < 90))
	      {
		leadToNext[17]->Fill(a * (0.6 / bining_plot),rings[0]);
		leadToNext[18]->Fill(a * (0.6 / bining_plot),ringsRaw[0]);
	      }
	    if((etjet > 90) && (etjet < 300))
	      {
		leadToNext[19]->Fill(a * (0.6 / bining_plot),rings[0]);
		leadToNext[20]->Fill(a * (0.6 / bining_plot),ringsRaw[0]);
	      }
	  }
	else
	  {
	    leadToNext[0]->Fill(a * (0.6 / bining_plot),rings[a]-rings[a-1]);
	    leadToNext[1]->Fill(a * (0.6 / bining_plot),ringsRaw[a]-ringsRaw[a-1]);
	    if((etjet > 10) && (etjet < 35))
	      {
		leadToNext[15]->Fill(a * (0.6 / bining_plot),rings[a]-rings[a-1]);
		leadToNext[16]->Fill(a * (0.6 / bining_plot),ringsRaw[a]-ringsRaw[a-1]);
	      }
	    if((etjet > 35) && (etjet < 90))
	      {
		leadToNext[17]->Fill(a * (0.6 / bining_plot),rings[a]-rings[a-1]);
		leadToNext[18]->Fill(a * (0.6 / bining_plot),ringsRaw[a]-ringsRaw[a-1]);
	      }
	    if((etjet > 90) && (etjet < 300))
	      {
		leadToNext[19]->Fill(a * (0.6 / bining_plot),rings[a]-rings[a-1]);
		leadToNext[20]->Fill(a * (0.6 / bining_plot),ringsRaw[a]-ringsRaw[a-1]);
	      }
	    if(rings[0]!=0)
	      {
		leadToNext[13]->Fill(a * (0.6 / bining_plot),(rings[a]-rings[a-1])/rings[0]);
		leadToNext[14]->Fill(a * (0.6 / bining_plot),(ringsRaw[a]-ringsRaw[a-1])/ringsRaw[0]);
	      }
	  }
	ringsSum[0][0] +=rings[0];
	ringsRawSum[0][0] +=ringsRaw[0];
	for(int b=1;b<11;b++)
	  {
	    if((fabs(etajet) > (b-1)*0.5) && (fabs(etajet) < b*0.5))
	      {
		ringsSum[b][0] +=rings[0];
		ringsRawSum[b][0] +=ringsRaw[0];
	      }
	  }
	if(a>0)
	  {
	    ringsSum[0][a] +=rings[a] - rings[a-1];
	    ringsRawSum[0][a] +=ringsRaw[a] - ringsRaw[a-1];
	    for(int b=1;b<11;b++)
	      {
		if((fabs(etajet) > (b-1)*0.5) && (fabs(etajet) < b*0.5))
		  {
		    ringsSum[b][a] +=rings[a] - rings[a-1];
		    ringsRawSum[b][a] +=ringsRaw[a] - ringsRaw[a-1];
		  }
	      }
	  }
      }
    towerinjet[0]->Fill(noTower);
    
    heta[0]->Fill(etajet,etjet/ jg->GetTruth(),jg->GetWeight());
    heta[1]->Fill(etajet,etjetcor/ jg->GetTruth(),jg->GetWeight());
    heta[2]->Fill(etajet,etjet/etjetcor,jg->GetWeight());
    if (jg->GetTruth() > 10 && jg->GetTruth() < 35) {
      heta[3]->Fill(etajet,etjet/ jg->GetTruth(),jg->GetWeight());
      heta[4]->Fill(etajet,etjetcor/ jg->GetTruth(),jg->GetWeight());
      heta[5]->Fill(etajet,etjet/etjetcor,jg->GetWeight());
      towerinjet[1]->Fill(noTower);
    } else if (jg->GetTruth() > 35 && jg->GetTruth() < 90) {
      heta[6]->Fill(etajet,etjet/ jg->GetTruth(),jg->GetWeight());
      heta[7]->Fill(etajet,etjetcor/ jg->GetTruth(),jg->GetWeight());
      heta[8]->Fill(etajet,etjet/etjetcor,jg->GetWeight());
      towerinjet[2]->Fill(noTower);
    } else if (jg->GetTruth() > 90 && jg->GetTruth() < 300) {
      heta[9]->Fill(etajet,etjet/ jg->GetTruth(),jg->GetWeight());
      heta[10]->Fill(etajet,etjetcor/ jg->GetTruth(),jg->GetWeight());
      heta[11]->Fill(etajet,etjet/etjetcor,jg->GetWeight());
      towerinjet[3]->Fill(noTower);
    } 





    hpt[0]->Fill(jg->GetTruth(),etjet/ jg->GetTruth(),jg->GetWeight());
    hpt[1]->Fill(jg->GetTruth(),etjetcor/jg->GetTruth(),jg->GetWeight());
    hpt[2]->Fill(jg->GetTruth(),etjet/etjetcor,jg->GetWeight());   




    henergy[0]->Fill(energyjet,etjet/ jg->GetTruth(),jg->GetWeight());
    henergy[1]->Fill(energyjet,etjetcor/jg->GetTruth(),jg->GetWeight());
    henergy[2]->Fill(energyjet,etjet/etjetcor,jg->GetWeight());    
    hptlog[0]->Fill(jg->GetTruth(),etjet/ jg->GetTruth(),jg->GetWeight());
    hptlog[1]->Fill(jg->GetTruth(),etjetcor/jg->GetTruth(),jg->GetWeight());
    hptlog[2]->Fill(etjetcor,etjet/etjetcor,jg->GetWeight());

    hpt_uncorr[0]->Fill(etjet,etjet/ jg->GetTruth(),jg->GetWeight());
    hpt_uncorr[1]->Fill(etjet,etjetcor/jg->GetTruth(),jg->GetWeight());
    hpt_uncorr[2]->Fill(etjet,etjet/etjetcor,jg->GetWeight());    
    if(fabs(etajet) < 1.4)      {
      hpt_uncorr[3]->Fill(etjet,etjet/ jg->GetTruth(),jg->GetWeight());
      hpt_uncorr[4]->Fill(etjet,etjetcor/jg->GetTruth(),jg->GetWeight());
      hpt_uncorr[5]->Fill(etjet,etjet/etjetcor,jg->GetWeight());    
    } else if(fabs(etajet) > 1.4  && fabs(etajet) < 3.0) {
      hpt_uncorr[6]->Fill(etjet,etjet/ jg->GetTruth(),jg->GetWeight());
      hpt_uncorr[7]->Fill(etjet,etjetcor/jg->GetTruth(),jg->GetWeight());
      hpt_uncorr[8]->Fill(etjet,etjet/etjetcor,jg->GetWeight());    
    } else if(fabs(etajet) > 3.0) {
      hpt_uncorr[9]->Fill(etjet,etjet/ jg->GetTruth(),jg->GetWeight());
      hpt_uncorr[10]->Fill(etjet,etjetcor/jg->GetTruth(),jg->GetWeight());
      hpt_uncorr[11]->Fill(etjet,etjet/etjetcor,jg->GetWeight());    
    }

/*
    hptlog[0]->Fill(jg->GetMess()[0],etjet/ jg->GetTruth(),jg->GetWeight());
    hptlog[1]->Fill(jg->GetMess()[0],etjetcor/jg->GetTruth(),jg->GetWeight());
    hptlog[2]->Fill(etjetcor,etjet/etjetcor,jg->GetWeight());
*/
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
    hptGamma[0]->Fill(jg->GetTruth());
    hptGammaW[0]->Fill(jg->GetTruth(),jg->GetWeight());
    hptGamma[1]->Fill(etajet);
    hptGammaW[1]->Fill(etajet,jg->GetWeight());
    hptGamma[2]->Fill(em/(em+had));
    hptGammaW[2]->Fill(em/(em+had),jg->GetWeight());
    hptGamma2D->Fill(jg->GetTruth(),em/(em+had));
    hptGamma2DW->Fill(jg->GetTruth(),em/(em+had),jg->GetWeight());
  } 

  TH1F* hists[12][8];
  c1->cd();
  for(int i=0;i<4;++i)
    {
      towerinjet[i]->Draw();
      c1->Draw();
      ps.NewPage();
    }
  /*
  for(int i=0;i<4;++i)
    {
      respvstet[i]->Draw("Box");
      c1->Draw();
      ps.NewPage();
    }
  */
  /*
  for(int i=0;i<3;++i)
    {
      testNegResponse[i]->Draw("Box");
      c1->Draw();
      ps.NewPage();
    }
  */
  for(int b=0; b< 11; ++b)
    {
      for(int a=0; a < bining_plot; ++a)
	{
	  if(ringsRawSum[b][a] != 0)
	    leadToNext[b+2]->Fill(a * (0.6 / bining_plot),ringsSum[b][a]/ringsRawSum[b][a]);
	}
    }
  leadToNext[0]->Draw();
  leadToNext[1]->SetLineColor(2);
  leadToNext[1]->Draw("same");
  c1->Draw();
  ps.NewPage();
  /*
  for(int b=2; b< 15; ++b)
    {
      if(b==13)
	{
	  leadToNext[13]->Draw();
	  leadToNext[14]->SetLineColor(2);
	  leadToNext[14]->Draw("same");
	  c1->Draw();
	  ps.NewPage();
	}
      //if((b!=13) && (b!=14))
	//{
	  //leadToNext[b]->Draw();
	  //c1->Draw();
	  //ps.NewPage();
	//}
    }
*/
  for(int i=0; i<3;++i)
    {
      leadToNext[15 + (2*i)]->Draw();
      leadToNext[16 + (2*i)]->SetLineColor(2);
      leadToNext[16 + (2*i)]->Draw("same");
      c1->Draw();
      ps.NewPage();
    }

  EtaPhiMap->Draw("Col");//,Palette");
  c1->Draw();
  ps.NewPage();
  EtaPhiMap->Draw("Palette");
  c1->Draw();
  ps.NewPage();

  TH1F* gaussplots[3][4];  
  TF1* gf[3][4];          

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

    Fit2D(heta[i],hists[i],gaussplots[0], gf[0]);
    Fit2D(heta[i+1],hists[i+1],gaussplots[1], gf[1]);
    Fit2D(heta[i+2],hists[i+2],gaussplots[2], gf[2]); 
    for(int a = 0; a<3;++a)
      {
	for(int b = 0 ; b < 4 ; ++b) {
	  hists[a+i][b]->SetMinimum(0.2);
	  hists[a+i][b]->SetMaximum(1.8);
	  ++b;
	  hists[a+i][b]->SetMinimum(0.0);
	  hists[a+i][b]->SetMaximum(0.5);
	}
      }


  for(int k=0;k<3;++k)
    {
    gaussplots[0][k]->SetMarkerStyle(20);
    gaussplots[0][k]->SetMarkerColor(1);
    gaussplots[1][k]->SetMarkerStyle(22);
    gaussplots[1][k]->SetMarkerColor(2);
    gaussplots[2][k]->SetMarkerStyle(4);
    gaussplots[2][k]->SetMarkerColor(4);
    }
  /*
  for(int a=0; a<1;++a) // 1: pt(jet)/Et(gamma), 2: pt(jet)/pt(corJet), 3: pt(corJet)/Et(gamma)
      {
	for(int b=0;b<3;++b) 
	  {
	    if(i==0)    
	      gaussplots[a][b]->SetTitle("#gamma-jet full energy range;#eta");
	    if(i==3)    
	      gaussplots[a][b]->SetTitle("#gamma-jet 10 < E_{T}^{#gamma} < 35 GeV;#eta");
	    if(i==6)    
	      gaussplots[a][b]->SetTitle("#gamma-jet 35 < E_{T}^{#gamma} < 90 GeV;#eta");
	    if(i==9)    
	      gaussplots[a][b]->SetTitle("#gamma-jet 90 < E_{T}^{#gamma} < 300 GeV;#eta");
	    c2->cd(b);
	    gaussplots[a][b]->Draw();
	    gf[a][b]->Draw("same");
	    c2->Update();
	  }
	c2->Draw();
	ps.NewPage();
      }
    c1->cd();
  */
    for(int j = 0 ; j < 8 ; ++j) {
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
    for(int i = 2 ; i < 8 ; ++i) {
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
    for(int j = 0 ; j < 8 ; ++j) {
      delete hists[i][j];
    }	
  }


  for(int i = 0 ; i < 12 ; i+= 3) {
    hpt_uncorr[i]->SetMarkerStyle(20);
    hpt_uncorr[i]->SetMarkerColor(1);
    hpt_uncorr[i]->SetMinimum(0.5);
    hpt_uncorr[i]->SetMaximum(1.2);
    hpt_uncorr[i+2]->SetMarkerStyle(4);
    hpt_uncorr[i+2]->SetMarkerColor(4);  
    hpt_uncorr[i+1]->SetMarkerStyle(22);
    hpt_uncorr[i+1]->SetMarkerColor(2);
    hpt_uncorr[i]->SetMinimum(0.5);
    hpt_uncorr[i]->SetMaximum(1.2);
    leg->Draw();

    Fit2D(hpt_uncorr[i],hists[i],gaussplots[0], gf[0]);
    Fit2D(hpt_uncorr[i+1],hists[i+1],gaussplots[1], gf[1]);
    Fit2D(hpt_uncorr[i+2],hists[i+2],gaussplots[2], gf[2]); 
    for(int a = 0; a<3;++a)
      {
	for(int b = 0 ; b < 4 ; ++b) {
	  hists[a+i][b]->SetMinimum(0.2);
	  hists[a+i][b]->SetMaximum(1.8);
	  ++b;
	  hists[a+i][b]->SetMinimum(0.0);
	  hists[a+i][b]->SetMaximum(0.5);
	}
      }


  for(int k=0;k<3;++k)
    {
    gaussplots[0][k]->SetMarkerStyle(20);
    gaussplots[0][k]->SetMarkerColor(1);
    gaussplots[1][k]->SetMarkerStyle(22);
    gaussplots[1][k]->SetMarkerColor(2);
    gaussplots[2][k]->SetMarkerStyle(4);
    gaussplots[2][k]->SetMarkerColor(4);
    }
  /*
  for(int a=0; a<1;++a) // 1: pt(jet)/Et(gamma), 2: pt(jet)/pt(corJet), 3: pt(corJet)/Et(gamma)
      {
	for(int b=0;b<3;++b) 
	  {
	    if(i==0)    
	      gaussplots[a][b]->SetTitle("#gamma-jet full energy range;#eta");
	    if(i==3)    
	      gaussplots[a][b]->SetTitle("#gamma-jet 10 < E_{T}^{#gamma} < 35 GeV;#eta");
	    if(i==6)    
	      gaussplots[a][b]->SetTitle("#gamma-jet 35 < E_{T}^{#gamma} < 90 GeV;#eta");
	    if(i==9)    
	      gaussplots[a][b]->SetTitle("#gamma-jet 90 < E_{T}^{#gamma} < 300 GeV;#eta");
	    c2->cd(b);
	    gaussplots[a][b]->Draw();
	    gf[a][b]->Draw("same");
	    c2->Update();
	  }
	c2->Draw();
	ps.NewPage();
      }
    c1->cd();
  */
    for(int j = 0 ; j < 8 ; ++j) {
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
  for(int i = 0 ; i < 12 ; ++i) {
    for(int j = 0 ; j < 8 ; ++j) {
      delete hists[i][j];
    }	
  }




  hpt[0]->SetMarkerStyle(20);
  hpt[0]->SetMarkerColor(1);
  hpt[2]->SetMarkerStyle(4);
  hpt[2]->SetMarkerColor(4);  
  hpt[1]->SetMarkerStyle(22);
  hpt[1]->SetMarkerColor(2);
  Fit2D(hpt[0],hists[0],gaussplots[0], gf[0]);
  Fit2D(hpt[1],hists[1],gaussplots[1], gf[1]);
  Fit2D(hpt[2],hists[2],gaussplots[2], gf[2]);
    for(int a = 0; a<3;++a)
      {
	for(int b = 0 ; b < 4 ; ++b) {
	  hists[a][b]->SetMinimum(0.2);
	  hists[a][b]->SetMaximum(1.8);
	  ++b;
	  hists[a][b]->SetMinimum(0.0);
	  hists[a][b]->SetMaximum(0.5);
	}
      }
  for(int k=0;k<3;++k)
    {
    gaussplots[0][k]->SetMarkerStyle(20);
    gaussplots[0][k]->SetMarkerColor(1);
    gaussplots[1][k]->SetMarkerStyle(22);
    gaussplots[1][k]->SetMarkerColor(2);
    gaussplots[2][k]->SetMarkerStyle(4);
    gaussplots[2][k]->SetMarkerColor(4);
    }
  /*
    for(int a=0; a<3;++a)
      {
	for(int b=0;b<3;++b)
	  { 
	    gaussplots[a][b]->SetTitle("#gamma-jet;p_{T} [GeV]");
	    c2->cd(b);
	    gaussplots[a][b]->Draw();
	    gf[a][b]->Draw("same");
	    c2->Update();
	  }
      }
    c2->Draw();
    ps.NewPage();
    c1->cd();
  */


  for(int i = 0 ; i < 8 ; ++i) {
    hists[0][i]->Draw();
    hists[0][i]->SetStats(0);
    hists[1][i]->Draw("SAME");
    hists[2][i]->Draw("SAME");
    leg->Draw();
    c1->SetGrid();
    c1->Draw();   
    ps.NewPage(); 
  }
  for(int i = 0 ; i < 8 ; ++i) {
    delete hists[0][i];
    delete hists[1][i];
    delete hists[2][i];
  }



  henergy[0]->SetMarkerStyle(20);
  henergy[0]->SetMarkerColor(1);
  henergy[2]->SetMarkerStyle(4);
  henergy[2]->SetMarkerColor(4);  
  henergy[1]->SetMarkerStyle(22);
  henergy[1]->SetMarkerColor(2);
  Fit2D(henergy[0],hists[0],gaussplots[0], gf[0]);
  Fit2D(henergy[1],hists[1],gaussplots[1], gf[1]);
  Fit2D(henergy[2],hists[2],gaussplots[2], gf[2]);
    for(int a = 0; a<3;++a)
      {
	for(int b = 0 ; b < 4 ; ++b) {
	  hists[a][b]->SetMinimum(0.2);
	  hists[a][b]->SetMaximum(1.8);
	  ++b;
	  hists[a][b]->SetMinimum(0.0);
	  hists[a][b]->SetMaximum(0.5);
	}
      }
  for(int k=0;k<3;++k)
    {
    gaussplots[0][k]->SetMarkerStyle(20);
    gaussplots[0][k]->SetMarkerColor(1);
    gaussplots[1][k]->SetMarkerStyle(22);
    gaussplots[1][k]->SetMarkerColor(2);
    gaussplots[2][k]->SetMarkerStyle(4);
    gaussplots[2][k]->SetMarkerColor(4);
    }
  /*
    for(int a=0; a<3;++a)
      {
	for(int b=0;b<3;++b)
	  { 
	    gaussplots[a][b]->SetTitle("#gamma-jet;Energy [GeV]");
	    c2->cd(b);
	    gaussplots[a][b]->Draw();
	    gf[a][b]->Draw("same");
	    c2->Update();
	  }
      }
    c2->Draw();
    ps.NewPage();
    c1->cd();
  */

  for(int i = 0 ; i < 8 ; ++i) {
    hists[0][i]->Draw();
    hists[0][i]->SetStats(0);
    hists[1][i]->Draw("SAME");
    hists[2][i]->Draw("SAME");
    leg->Draw();
    c1->SetGrid();
    c1->Draw();   
    ps.NewPage(); 
  }
  for(int i = 0 ; i < 8 ; ++i) {
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
  Fit2D(hemf[0],hists[0],gaussplots[0], gf[0]);
  Fit2D(hemf[1],hists[1],gaussplots[1], gf[1]);
  Fit2D(hemf[2],hists[2],gaussplots[2], gf[2]);
    for(int a = 0; a<3;++a)
      {
	for(int b = 0 ; b < 4 ; ++b) {
	  hists[a][b]->SetMinimum(0.2);
	  hists[a][b]->SetMaximum(1.8);
	  ++b;
	  hists[a][b]->SetMinimum(0.0);
	  hists[a][b]->SetMaximum(0.5);
	}
      }
  for(int k=0;k<3;++k)
    {
    gaussplots[0][k]->SetMarkerStyle(20);
    gaussplots[0][k]->SetMarkerColor(1);
    gaussplots[1][k]->SetMarkerStyle(22);
    gaussplots[1][k]->SetMarkerColor(2);
    gaussplots[2][k]->SetMarkerStyle(4);
    gaussplots[2][k]->SetMarkerColor(4);
    }
  /*
    for(int a=0; a<1;++a) 
      {
	for(int b=0;b<3;++b)
	  {
	    gaussplots[a][b]->SetTitle("#gamma-jet;EMF");
	    c2->cd(b);
	    gaussplots[a][b]->Draw();
	    gf[a][b]->Draw("same");
	    c2->Update();
	  }
      }
    c2->Draw();
    ps.NewPage();
    c1->cd();
  */

  for(int i = 0 ; i < 8 ; ++i) {
    hists[0][i]->Draw();
    hists[0][i]->SetStats(0);
    hists[1][i]->Draw("SAME");
    hists[2][i]->Draw("SAME");
    leg->Draw();
    c1->SetGrid();
    c1->Draw();   
    ps.NewPage(); 
  }
  for(int i = 0 ; i < 8 ; ++i) {
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
  Fit2D(hptlog[0],hists[0],gaussplots[0], gf[0]);
  Fit2D(hptlog[1],hists[1],gaussplots[1], gf[1]);
  Fit2D(hptlog[2],hists[2],gaussplots[2], gf[2]);
    for(int a = 0; a<3;++a)
      {
	for(int b = 0 ; b < 4 ; ++b) {
	  hists[a][b]->SetMinimum(0.2);
	  hists[a][b]->SetMaximum(1.8);
	  ++b;
	  hists[a][b]->SetMinimum(0.0);
	  hists[a][b]->SetMaximum(0.5);
	}
      }
  for(int k=0;k<3;++k)
    {
    gaussplots[0][k]->SetMarkerStyle(20);
    gaussplots[0][k]->SetMarkerColor(1);
    gaussplots[1][k]->SetMarkerStyle(22);
    gaussplots[1][k]->SetMarkerColor(2);
    gaussplots[2][k]->SetMarkerStyle(4);
    gaussplots[2][k]->SetMarkerColor(4);
    }
  /*
    for(int a=0; a<1;++a)
      {
	for(int b=0;b<3;++b) 
	  {
	    gaussplots[a][b]->SetTitle("#gamma-jet;p_{T} (log)[GeV]");
	    c2->cd(b);
	    gaussplots[a][b]->Draw();
	    gf[a][b]->Draw("same");
	    c2->Update();
	  }
      }
    c2->Draw();
    ps.NewPage();
    c1->cd();
  */
    for(int i = 0 ; i < 3 ; ++i) 
      {
	for(int j=0; j<3;++j)
	  {
	    delete gaussplots[i][j];
	    delete gf[i][j];
	  }
      }

  for(int i = 0 ; i < 8 ; ++i) {
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
  for(int i = 0 ; i < 8 ; ++i) {
    delete hists[0][i];
    delete hists[1][i];
    delete hists[2][i];
  }
  ps.NewPage(); 
  delete leg;
  for(int i = 0 ; i < 3 ; ++i) {
    hptGamma[i]->SetMarkerStyle(20);
    hptGamma[i]->SetMarkerColor(1);
    hptGammaW[i]->SetMarkerStyle(22);
    hptGammaW[i]->SetMarkerColor(2);
    if (hptGamma[i]->GetMaximum()>hptGammaW[i]->GetMaximum())
      hptGammaW[i]->SetMaximum( hptGamma[i]->GetMaximum() );

    hptGammaW[i]->Draw("p");
    hptGammaW[i]->SetStats(0);
    hptGamma[i]->Draw("pSAME");
    c1->SetLogx(0);  
    c1->SetLogy(1);  
    c1->SetGrid();
    leg = new TLegend(0.7,0.96,0.96,0.72);
    leg->AddEntry(hptGamma[i],"no weights","p");
    leg->AddEntry(hptGammaW[i],"with weights","p");
    leg->Draw();
    c1->Draw(); 
  }
  c1->SetLogx(0);  
  c1->SetLogy(0);   
  c1->SetGrid(0);
  hptGamma2D->Draw("hist");
  c1->Draw(); 
   ps.NewPage(); 
   
  hptGamma2DW->Draw("hist");
  c1->Draw(); 
  
  ps.NewPage(); 
  delete leg;
  for(int i = 0 ; i < 12 ; ++i) delete heta[i];
  for(int i = 0 ; i < 12 ; ++i) delete hpt_uncorr[i];
  delete hpt[0];
  delete hpt[1];
  delete hpt[2];
  delete henergy[0];
  delete henergy[1];
  delete henergy[2];
  delete hptlog[0];
  delete hptlog[1];
  delete hptlog[2];
  delete hptGamma[0];
  delete hptGammaW[0];  
  delete hptGamma[1];
  delete hptGammaW[1];  
  delete hptGamma[2];
  delete hptGammaW[2];  
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



void TControlPlots::DiJetControlPlots()
{
  TCanvas * c1 = new TCanvas("controlplots","",600,600);
  TCanvas * c2 = new TCanvas("controlplotsGauss","",600,600);
  c2->Divide(1,3);
  TPostScript ps("controlplotsDiJet.ps",111);

  //book hists
  TH2F* combmean[5];
  TH2F* difmean[5];
  TH2F* Bvsdphi;
  TH2F* Difvscomb;
  TH1F* eta[2];
  TH1F* ptspec[2];
  TH1F* dphi[11];
  
  ptspec[0] = new TH1F("Pt spectrum probe jet","Pt spectrum probe jet;Pt",300,0,2000);
  ptspec[1] = new TH1F("Pt spectrum barrel jet","Pt spectrum barrel jet;Pt",300,0,2000);
  eta[0] = new TH1F("eta probe jet","#eta probe jet;#eta",100,-5,5);
  eta[1] = new TH1F("eta barrel jet","#eta barrel jet;#eta",100,-5,5);
  dphi[0] = new TH1F("Delta Phi","#Delta Phi;#Delta #Phi",70,2.8,3.5);
  for(int i=1;i<5;++i)
    dphi[i] = (TH1F*)dphi[0]->Clone();
  dphi[1]->SetTitle("#Delta Phi (10-35 GeV);#Delta #Phi");
  dphi[2]->SetTitle("#Delta Phi (35-90 GeV);#Delta #Phi");
  dphi[3]->SetTitle("#Delta Phi (90-300 GeV);#Delta #Phi");
  dphi[4]->SetTitle("#Delta Phi (300+ GeV);#Delta #Phi");
  dphi[5] = new TH1F("Delta Phi off","#Delta Phi;#Delta #Phi",120,-3.4,3.4);
  dphi[6] = new TH1F("Delta Phi off+","#Delta Phi;#Delta #Phi",120,2.8,3.4);
  dphi[7] = new TH1F("abs Delta Phi off","#Delta Phi;#Delta #Phi",120,2.8,3.4);
  dphi[8] = new TH1F("Delta Phi wo abs","#Delta Phi;#Delta #Phi",140,-3.5,3.5);
  dphi[9] = new TH1F("Delta Phi wo abs +","#Delta Phi;#Delta #Phi",120,2.8,3.5);
  dphi[10] = new TH1F("Delta Phi wo abs -","#Delta Phi;#Delta #Phi",120,-3.5,-2.8);
  Bvsdphi = new TH2F("B vs DeltaPhi","B vs #Delta #Phi;#Delta#Phi;B",70,2.8,3.5,100,-1,1);
  Difvscomb = new TH2F("Delta Et vs Combined Jet","#Delta Et vs. Combined Jet;Combined Jet;#Delta Et",100,0,100,100,0,100);
  combmean[0] = new TH2F("DiJet","di-jet controlplot;scale (Pt);Pt of combined jet",100,0,2000,100,0,200);
  difmean[0] = new TH2F("DiJet difference","di-jet controlplot;scale (Pt);abs. difference in Pt",100,0,2000,100,0,200);
  for(int i=1;i<5;++i)
    {
      combmean[i] =  (TH2F*)combmean[0]->Clone();
      difmean[i] =  (TH2F*)difmean[0]->Clone();
    }
  combmean[1]->SetTitle("abs(#eta) < 1");
  combmean[2]->SetTitle("1 < abs(#eta) < 2");
  combmean[3]->SetTitle("2 < abs(#eta) < 3");
  combmean[4]->SetTitle("3 < abs(#eta) < 4");
  difmean[1]->SetTitle("abs(#eta) < 1");
  difmean[2]->SetTitle("1 < abs(#eta) < 2");
  difmean[3]->SetTitle("2 < abs(#eta) < 3");
  difmean[4]->SetTitle("3 < abs(#eta) < 4");
  //combmean->SetMarkerStyle(20);
  //combmean->SetMarkerColor(1);

  //TLegend* leg = new TLegend(0.7,0.96,0.96,0.72);
  //leg->AddEntry(combmean,"Difference Vs. mean Pt","p");


  //book hists
  TH2F* Beta[8];
  Beta[0] = new TH2F("Beta","di-jet;#eta",100,-5,5,100,-1,1);
  for(int i = 1 ; i < 8 ; ++i) Beta[i] = (TH2F*)Beta[0]->Clone();
  Beta[2]->SetTitle("di-jet 10 < E_{T}^{barrel jet} < 35 GeV;#eta");
  Beta[4]->SetTitle("di-jet 35 < E_{T}^{barrel jet} < 90 GeV;#eta");
  Beta[6]->SetTitle("di-jet 90 < E_{T}^{barrel jet} < 300 GeV;#eta");

  TH2F* Bpt[2];
  Bpt[0] = new TH2F("Bpt","di-jet;p_{T} [GeV]",400,0,400,100,-1,1);
  Bpt[1] = (TH2F*)Bpt[0]->Clone();

  TH2F* Benergy[2];
  Benergy[0] = new TH2F("Benergy","di-jet;Energy [GeV]",400,0,400,100,-1,1);
  Benergy[1] = (TH2F*)Benergy[0]->Clone();
  
  TH2F* Bemf[2];
  Bemf[0] = new TH2F("Bemf","di-jet;EMF (probe jet)",100,0,1,100,-1,1);
  Bemf[1] = (TH2F*)Bemf[0]->Clone();

  double bins[101];
  for(int i = 0; i < 101 ; ++i) {
    bins[i] = pow(10,(i+32)/40.0);
  }

  TH2F* Bptlog[2];
  Bptlog[0] = new TH2F("Bptlog","di-jet;p_{T} [GeV]",100,bins,100,-1,1); 
  Bptlog[1] = (TH2F*)Bptlog[0]->Clone();


  for(int i=0; i<2;++i)
    {
	  Beta[4*i]->Sumw2();
	  Beta[4*i+1]->Sumw2();
	  Beta[4*i+2]->Sumw2();
	  Beta[4*i+3]->Sumw2();
	  Bpt[i]->Sumw2();
	  Benergy[i]->Sumw2();
	  Bemf[i]->Sumw2();
	  Bptlog[i]->Sumw2();
    }

  //loop over all fit-events
  for ( std::vector<TData*>::iterator i = data->begin() ; i != data->end() ; ++i )  
    {
      TData* jj = *i;
      if(jj->GetType() != TypePtBalance) continue;
      TData_MessMess* jm = (TData_MessMess*) jj;
      double etscale = jm->GetScale();
      double etajet1 = jm->GetMultMess(0)[1];
      double etajet2 = jm->GetMultMess(1)[1];
      double etjetcomb = jm->GetMessCombination();
      double etjet1 = jm->GetMultParametrizedMess(0);      //Probe
      double etjet2 = jm->GetMultParametrizedMess(1);      //Barrel
      double etjet1uncor = jm->GetMultMess(0)[0];      //Probe
      double etjet2uncor = jm->GetMultMess(1)[0];      //Barrel
      double phijet1 = jm->GetMultMess(0)[2];      //Probe
      double phijet2 = jm->GetMultMess(1)[2];      //Barrel
      double B = (etjet1 - etjet2) / etscale;
      double Buncor = (etjet1uncor - etjet2uncor) * 2 / (etjet1uncor + etjet2uncor);
      double etaprobe = etajet1;
      double phiprobe = phijet1;
      double etprobe = etjet1;
      if(fabs(etajet1) < fabs(etajet2))  //unbias if both jets in barrel
	{
	  B *= -1;
	  Buncor *= -1;
	  etprobe = etjet2; 
	  etjet2 = etjet1;  
	  etaprobe = etajet2; 
	  etajet2 = etajet1;  
	  phiprobe = phijet2; 
	  phijet2 = phijet1;  
	}
      double deltaphi = fabs(phiprobe - phijet2);
      //double etjet3uncor = jm->GetMultMess(2)[0];
      double deltaphioff = deltaPhi(phiprobe,phijet2);
      //Test->Fill(etjet3uncor,deltaphi);
      //combmean->Fill(etmean, etjetdif, jj->GetWeight());
      //ptspec[0]

      ptspec[0]->Fill(etprobe);
      ptspec[1]->Fill(etjet2);
      eta[0]->Fill(etaprobe);
      eta[1]->Fill(etajet2);
      dphi[0]->Fill(deltaphi);
      dphi[5]->Fill(deltaphioff);
      dphi[6]->Fill(deltaphioff);
      dphi[7]->Fill(fabs(deltaphioff));
      dphi[8]->Fill(phiprobe - phijet2);             //
      if((phiprobe - phijet2) > 0)
	dphi[9]->Fill(phiprobe - phijet2);             //
      else
	dphi[10]->Fill(phiprobe - phijet2);             //
      Bvsdphi->Fill(deltaphi,B); 
      Difvscomb->Fill(fabs(etprobe - etjet2),etjetcomb);
      combmean[0]->Fill(etscale, etjetcomb);
      difmean[0]->Fill(etscale, fabs(etprobe - etjet2));

      Beta[0]->Fill(etaprobe, B,jj->GetWeight());
      Beta[1]->Fill(etaprobe, Buncor,jj->GetWeight());
      if (etscale > 10 && etscale < 35) {
      Beta[2]->Fill(etaprobe, B,jj->GetWeight());
      Beta[3]->Fill(etaprobe, Buncor,jj->GetWeight());
      dphi[1]->Fill(deltaphi);
      } else if (etscale > 35 && etscale < 90) {
      Beta[4]->Fill(etaprobe, B,jj->GetWeight());
      Beta[5]->Fill(etaprobe, Buncor,jj->GetWeight());
      dphi[2]->Fill(deltaphi);
      } else if (etscale > 90 && etscale < 300) {
      Beta[6]->Fill(etaprobe, B,jj->GetWeight());
      Beta[7]->Fill(etaprobe, Buncor,jj->GetWeight());
      dphi[3]->Fill(deltaphi);
      }  else if (etscale > 300)
	dphi[4]->Fill(deltaphi);
      Bpt[0]->Fill(etscale,B,jj->GetWeight());
      Bpt[1]->Fill(etscale,Buncor,jj->GetWeight());

      double theta1 = 2 * atan(exp(-etaprobe));
      double theta2 = 2 * atan(exp(-etajet2));
      double energy1 = etprobe * sin(theta1);
      double energy2 = etjet2 * sin(theta2);
      Benergy[0]->Fill((energy1 + energy2) /2,B,jj->GetWeight());
      Benergy[1]->Fill((energy1 + energy2) /2,Buncor,jj->GetWeight());
      Bptlog[0]->Fill(etscale,B,jj->GetWeight());
      Bptlog[1]->Fill(etscale,Buncor,jj->GetWeight());


      //em fraction plots     
      double em = 0;
      double had = 0;
      for(std::vector<TData*>::const_iterator t = jj->GetRef().begin(); t != jj->GetRef().end(); ++t) {
	TData* tt = *t;
	em  += tt->GetMess()[1];
	had += tt->GetMess()[2];
	had += tt->GetMess()[3];
      }
      Bemf[0]->Fill(em/(em+had),B,jj->GetWeight());
      Bemf[1]->Fill(em/(em+had),Buncor,jj->GetWeight());
            
      
      for(int i=0;i<4;++i)
	{
	  if((fabs(etaprobe) > i) && (fabs(etaprobe) < i+1)) 
	    {
	      combmean[i+1]->Fill(etscale, etjetcomb);
	      difmean[i+1]->Fill(etscale, fabs(etprobe - etjet2));
	    }
	}
    }


  c1->cd();
  for(int i=0;i<2;++i)
    {
      ptspec[i]->Draw();
      c1->Draw();
      ps.NewPage(); 
    }
  for(int i=0;i<2;++i)
    {
      eta[i]->Draw();
      c1->Draw();
      ps.NewPage(); 
    }
  for(int i=0;i<5;++i)
    {
      dphi[i]->Draw();
      c1->Draw();
      ps.NewPage(); 
    }
  Bvsdphi->Draw();
  c1->Draw();
  ps.NewPage();  
  Difvscomb->Draw();
  c1->Draw();
  ps.NewPage();  

  TH1F* hists[8][8];

  TH1F* gaussplots[2][4];  
  TF1* gf[2][4];     

  TLegend* leg = new TLegend(0.7,0.96,0.96,0.72);
  leg->AddEntry(Beta[0],"B = Pt^{probe} - Pt^{barrel} / scale");
  leg->AddEntry(Beta[1],"B before fit");
  for(int i = 0 ; i < 8 ; i+=2) {
    Beta[i]->SetMarkerStyle(20);
    Beta[i]->SetMarkerColor(1);
    Beta[i+1]->SetMarkerStyle(22);
    Beta[i+1]->SetMarkerColor(2);
    leg->Draw();

    Fit2D(Beta[i],hists[i],gaussplots[0], gf[0]);
    Fit2D(Beta[i+1],hists[i+1],gaussplots[1], gf[1]);
    for(int a = 0; a<2;++a)
      {
	for(int b = 0 ; b < 4 ; ++b) {
	  hists[a+i][b]->SetMinimum(-1.);
	  hists[a+i][b]->SetMaximum(1.);
	  ++b;
	  hists[a+i][b]->SetMinimum(0.0);
	  hists[a+i][b]->SetMaximum(1.);
	}
      }
    
    for(int k=0;k<3;++k)
      {
	gf[1][k]->SetLineColor(2);
	gaussplots[0][k]->SetMarkerStyle(20);
	gaussplots[0][k]->SetMarkerColor(1);
	gaussplots[1][k]->SetMarkerStyle(22);
	gaussplots[1][k]->SetMarkerColor(2);
      }
    /*
    for(int b=0;b<3;++b) 
      {
	for(int a=0; a<2;++a) // 1: B after correction, 2: B before correction
	  {
	    if(i==0)    
	      gaussplots[a][b]->SetTitle("di-jet full energy range;#eta");
	    if(i==2)    
	      gaussplots[a][b]->SetTitle("di-jet 10 < P_{T}^{scale} < 35 GeV;#eta");
	    if(i==4)    
	      gaussplots[a][b]->SetTitle("di-jet 35 < P_{T}^{scale} < 90 GeV;#eta");
	    if(i==6)    
	      gaussplots[a][b]->SetTitle("di-jet 90 < P_{T}^{scale} < 300 GeV;#eta");
	    c2->cd(b);
	    if(a==0)	    gaussplots[a][b]->Draw();
	    else            gaussplots[a][b]->Draw("same");
	    gf[a][b]->Draw("same");
	  }
	c2->Update();
	c2->Draw();
	ps.NewPage();
      }
    c1->cd();
    */
    for(int j = 0 ; j < 8 ; ++j) 
      {
	hists[i][j]->Draw();
	hists[i][j]->SetStats(0);
	hists[i+1][j]->Draw("SAME");
	leg->Draw();
	c1->SetGrid();
	c1->Draw();   
	ps.NewPage(); 
      }
  } 

  for(int i = 0 ; i < 8 ; ++i) {
    for(int j = 0 ; j < 8 ; ++j) {
      delete hists[i][j];
    }	
  }
  
  Bpt[0]->SetMarkerStyle(20);
  Bpt[0]->SetMarkerColor(1);
  Bpt[1]->SetMarkerStyle(22);
  Bpt[1]->SetMarkerColor(2);
  Fit2D(Bpt[0],hists[0],gaussplots[0], gf[0]);
  Fit2D(Bpt[1],hists[1],gaussplots[1], gf[1]);
    for(int a = 0; a<2;++a)
      {
	for(int b = 0 ; b < 4 ; ++b) {
	  hists[a][b]->SetMinimum(-1.);
	  hists[a][b]->SetMaximum(1.);
	  ++b;
	  hists[a][b]->SetMinimum(0.0);
	  hists[a][b]->SetMaximum(1.);
	}
      }

  for(int k=0;k<3;++k)
    {
      gf[1][k]->SetLineColor(2);
      gaussplots[0][k]->SetMarkerStyle(20);
      gaussplots[0][k]->SetMarkerColor(1);
      gaussplots[1][k]->SetMarkerStyle(22);
      gaussplots[1][k]->SetMarkerColor(2);
    }
  /*
  for(int b=0;b<3;++b)
    { 
      for(int a=0; a<2;++a)
	{
	  gaussplots[a][b]->SetTitle("di-jet;p_{T} [GeV]");
	  c2->cd(b);
	  if(a==0)	    gaussplots[a][b]->Draw();
	  else            gaussplots[a][b]->Draw("same");
	  gf[a][b]->Draw("same");
	}
      c2->Update();
      c2->Draw();
      ps.NewPage();
    }
  c1->cd();
  */
  for(int i = 0 ; i < 8 ; ++i) {
    hists[0][i]->Draw();
    hists[0][i]->SetStats(0);
    hists[1][i]->Draw("same");
    leg->Draw();
    c1->SetGrid();
    c1->Draw();   
    ps.NewPage(); 
  }
  for(int i = 0 ; i < 8 ; ++i) {
    delete hists[0][i];
    delete hists[1][i];
  }
  

  
  Benergy[0]->SetMarkerStyle(20);
  Benergy[0]->SetMarkerColor(1);
  Benergy[1]->SetMarkerStyle(22);
  Benergy[1]->SetMarkerColor(2);
  Fit2D(Benergy[0],hists[0],gaussplots[0], gf[0]);
  Fit2D(Benergy[1],hists[1],gaussplots[1], gf[1]);
    for(int a = 0; a<2;++a)
      {
	for(int b = 0 ; b < 4 ; ++b) {
	  hists[a][b]->SetMinimum(-1.);
	  hists[a][b]->SetMaximum(1.);
	  ++b;
	  hists[a][b]->SetMinimum(0.0);
	  hists[a][b]->SetMaximum(1.);
	}
      }

  for(int k=0;k<3;++k)
    {
      gf[1][k]->SetLineColor(2);
      gaussplots[0][k]->SetMarkerStyle(20);
      gaussplots[0][k]->SetMarkerColor(1);
      gaussplots[1][k]->SetMarkerStyle(22);
      gaussplots[1][k]->SetMarkerColor(2);
    }
  /*
  for(int b=0;b<3;++b)
    { 
      for(int a=0; a<2;++a)
	{
	  gaussplots[a][b]->SetTitle("di-jet;p_{T} [GeV]");
	  c2->cd(b);
	  if(a==0)	    gaussplots[a][b]->Draw();
	  else            gaussplots[a][b]->Draw("same");
	  gf[a][b]->Draw("same");
	}
      c2->Update();
      c2->Draw();
      ps.NewPage();
    }
  c1->cd();
  */
  for(int i = 0 ; i < 8 ; ++i) {
    hists[0][i]->Draw();
    hists[0][i]->SetStats(0);
    hists[1][i]->Draw("same");
    leg->Draw();
    c1->SetGrid();
    c1->Draw();   
    ps.NewPage(); 
  }
  for(int i = 0 ; i < 8 ; ++i) {
    delete hists[0][i];
    delete hists[1][i];
  }

  Bemf[0]->SetMarkerStyle(20);
  Bemf[0]->SetMarkerColor(1);
  Bemf[1]->SetMarkerStyle(22);
  Bemf[1]->SetMarkerColor(2);
  Fit2D(Bemf[0],hists[0],gaussplots[0], gf[0]);
  Fit2D(Bemf[1],hists[1],gaussplots[1], gf[1]);

    for(int a = 0; a<2;++a)
      {
	for(int b = 0 ; b < 4 ; ++b) {
	  hists[a][b]->SetMinimum(-1.);
	  hists[a][b]->SetMaximum(1.);
	  ++b;
	  hists[a][b]->SetMinimum(0.0);
	  hists[a][b]->SetMaximum(1.);
	}
      }
  for(int k=0;k<3;++k)
    {
      gf[1][k]->SetLineColor(2);
      gaussplots[0][k]->SetMarkerStyle(20);
      gaussplots[0][k]->SetMarkerColor(1);
      gaussplots[1][k]->SetMarkerStyle(22);
      gaussplots[1][k]->SetMarkerColor(2);
    }
  /*
  for(int b=0;b<3;++b)
    { 
      for(int a=0; a<2;++a) 
	{
	  gaussplots[a][b]->SetTitle("di-jet;p_{T} [GeV]");
	  c2->cd(b);
	  if(a==0)	    gaussplots[a][b]->Draw();
	  else            gaussplots[a][b]->Draw("same");
	  gf[a][b]->Draw("same");
	}
      c2->Update();
      c2->Draw();
      ps.NewPage();
    }
  c1->cd();
  */

  for(int i = 0 ; i < 8 ; ++i) {
    hists[0][i]->Draw();
    hists[0][i]->SetStats(0);
    hists[1][i]->Draw("SAME");
    leg->Draw();
    c1->SetGrid();
    c1->Draw();   
    ps.NewPage(); 
  }
  for(int i = 0 ; i < 8 ; ++i) {
    delete hists[0][i];
    delete hists[1][i];
  }
  Bptlog[0]->SetMarkerStyle(20);
  Bptlog[0]->SetMarkerColor(1);
  Bptlog[0]->SetMinimum(0.2);
  Bptlog[0]->SetMaximum(1.8);
  Bptlog[1]->SetMarkerStyle(22);
  Bptlog[1]->SetMarkerColor(2);  
  Fit2D(Bptlog[0],hists[0],gaussplots[0], gf[0]);
  Fit2D(Bptlog[1],hists[1],gaussplots[1], gf[1]);

    for(int a = 0; a<2;++a)
      {
	for(int b = 0 ; b < 4 ; ++b) {
	  hists[a][b]->SetMinimum(-1.);
	  hists[a][b]->SetMaximum(1.);
	  ++b;
	  hists[a][b]->SetMinimum(0.0);
	  hists[a][b]->SetMaximum(1.);
	}
      }
  for(int k=0;k<3;++k)
    {
      gf[1][k]->SetLineColor(2);
      gaussplots[0][k]->SetMarkerStyle(20);
      gaussplots[0][k]->SetMarkerColor(1);
      gaussplots[1][k]->SetMarkerStyle(22);
      gaussplots[1][k]->SetMarkerColor(2);
    }
  /*
  for(int b=0;b<3;++b)
    { 
      for(int a=0; a<2;++a)
	{
	  gaussplots[a][b]->SetTitle("di-jet;p_{T} [GeV]");
	  c2->cd(b);
	  if(a==0)	    gaussplots[a][b]->Draw();
	  else            gaussplots[a][b]->Draw("same");
	  gf[a][b]->Draw("same");
	}
      c2->Update();
      c2->Draw();
      ps.NewPage();
    }
  c1->cd();
  */
  for(int i = 0 ; i < 2 ; ++i) 
    {
      for(int j=0; j<3;++j)
	{
	  delete gaussplots[i][j];
	  delete gf[i][j];
	}
    }

  for(int i = 0 ; i < 8 ; ++i) {
    hists[0][i]->Draw();
    hists[0][i]->SetStats(0);
    hists[1][i]->Draw("SAME");
    c1->SetLogx(1);
    c1->SetGrid(); 
    leg->Draw();
    c1->SetGrid();
    c1->Draw();   
    ps.NewPage(); 
  }
  for(int i = 0 ; i < 8 ; ++i) {
    delete hists[0][i];
    delete hists[1][i];
  }
  ps.NewPage();
  delete leg;
  for(int i = 0 ; i < 8 ; ++i)  delete Beta[i];
  for(int i = 0 ; i < 2 ; ++i){
    delete Bpt[i];
    delete Benergy[i];
    delete Bptlog[i];
    delete Bemf[i];
  }
  
  /*
  c1->SetLogx(0);
  TF1* line = new TF1("line","x",0,200);
  for(int i=0;i<5;++i)
    {
      combmean[i]->Draw("Box");
      line->Draw("same");
      c1->Draw();   
      ps.NewPage();
    }
  for(int i=0;i<5;++i)
    {
      difmean[i]->Draw("Box");
      line->Draw("same");
      c1->Draw();   
      ps.NewPage();
    }
  */
  ps.Close();

  for(int i=0;i<5;++i)
    {
      delete combmean[i];
      delete difmean[i];
    }
  delete eta[0];
  delete eta[1];
  delete ptspec[0];
  delete ptspec[1];
  //delete line;
}


void TControlPlots::Fit2D(TH2F* hist, TH1F* hresults[8], TH1F* gaussplots[4], TF1* gf[4] )
{
  //book hists
  TString s(hist->GetName());
  s.Append("_res");
  s.Append(iplot++);
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
  for(int i = 1; i < 8 ; ++i) {
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
  hresults[4]->SetTitle(s.Append(" median"));  
  s = hist->GetTitle();
  hresults[5]->SetTitle(s.Append(" chisquared / no. free parameters")); 
  s = hist->GetTitle();
  hresults[6]->SetTitle(s.Append(" probability"));    
  /* for(int i = 0 ; i < 6 ; ++i) {
    hresults[i]->SetMinimum(0.2);
    hresults[i]->SetMaximum(1.8);
    ++i;
    hresults[i]->SetMinimum(0.0);
    hresults[i]->SetMaximum(0.3);
    }
*/
    hresults[5]->SetMinimum(0.0);
    hresults[5]->SetMaximum(100);
    hresults[6]->SetMinimum(0.0);
    hresults[6]->SetMaximum(1.05);
  TH1F* htemp = new TH1F("htemp","",hist->GetNbinsY(),hist->GetYaxis()->GetXmin(),
		          hist->GetYaxis()->GetXmax());
  htemp->Sumw2();
  for(int i=0;i<3;++i) 
    {
      gaussplots[i] = (TH1F*)htemp->Clone("hnew");
      gf[i] = new TF1("dummy","0",0,1);
    }

  const int nq = 2;
  double yq[2],xq[2];
  xq[0] = 0.5;
  xq[1] = 0.90;
  for(int i = 1 ; i <= hist->GetNbinsX() ; ++i) {
    htemp->Reset();
    for(int j = 1 ; j <= hist->GetNbinsY() ; ++j) {
      htemp->SetBinContent(j,hist->GetBinContent(hist->GetBin(i,j)));
      //htemp->Fill(htemp->GetBinCenter(j),hist->GetBinContent(hist->GetBin(i,j)));
      htemp->SetBinError(j,hist->GetBinError(i,j));
    }  

    if(htemp->GetSumOfWeights() <= 0) continue;
    htemp->Fit("gaus","LLQNO","");
    TF1 *f = (TF1*)gROOT->GetFunction("gaus")->Clone();
    double mean = f->GetParameter(1);
    double meanerror = f->GetParError(1);
    double width = f->GetParameter(2);
    if(width < 0.2) width = 0.2;
    if( (htemp->Fit(f,"LLQNO","goff",mean - 2 * width, mean + 2 * width) == 0) && (f->GetProb() > 0.01)) 
      {
	mean = f->GetParameter(1);
	meanerror = f->GetParError(1);
	width = f->GetParameter(2);
	hresults[2]->SetBinContent(i,mean);
	hresults[2]->SetBinError(i,meanerror);
	//hresults[3]->SetBinContent(i,width/mean);
	hresults[3]->SetBinContent(i,width);
	hresults[3]->SetBinError(i, f->GetParError(2));
      }
    hresults[5]->SetBinContent(i, f->GetChisquare() / f->GetNumberFreeParameters());
    hresults[5]->SetBinError(i, 0.01);
    hresults[6]->SetBinContent(i, f->GetProb());
    hresults[6]->SetBinError(i, 0.01);
    if(i ==  int(hist->GetNbinsX()/6))       
      {
	gaussplots[0] = (TH1F*)htemp->Clone("hnew"); 
	gf[0] = (TF1*)f->Clone();
      }

    if(i == int(hist->GetNbinsX()/3))  
      {
	gaussplots[1] = (TH1F*)htemp->Clone("hnew"); 
	gf[1] = (TF1*)f->Clone();
      }

    if(i == int(hist->GetNbinsX()/2))  
      {
	gaussplots[2] = (TH1F*)htemp->Clone("hnew"); 
	gf[2] = (TF1*)f->Clone();
      }
    mean = htemp->GetMean();
    meanerror = htemp->GetMeanError();
    width = htemp->GetRMS();
    hresults[0]->SetBinContent(i,mean);
    hresults[0]->SetBinError(i,meanerror);
    //hresults[1]->SetBinContent(i,width/mean); 
    hresults[1]->SetBinContent(i,width); 
    hresults[1]->SetBinError(i,htemp->GetRMSError());
    htemp->GetQuantiles(nq,yq,xq);
    hresults[4]->SetBinContent(i,yq[0]);
    hresults[4]->SetBinError(i,0.0001);
    hresults[7]->SetBinContent(i,yq[1]/yq[0]-1);
    hresults[7]->SetBinError(i,0.0001);
    delete f;
  }
  delete htemp;
}


double gauss_step(double *x, double *par);
//defined in caliber.cc

void TControlPlots::GammaJetSigmas()
{
  TCanvas * c1 = new TCanvas("controlplots","",600,600);
  TPostScript ps("sigmas_gammajet.ps",111);
  char * name = new char[100];

  TH1F * gauss_forpt[200];
  TH1F * gauss_forptcorr[200];
  gauss_forpt[0] = new TH1F("hgauss","pT bin[0..1GeV];#frac{pT jet - pT truth}{pT jet}",600,-3,3);
  gauss_forptcorr[0] = new TH1F("hgausscorr","corrected jet pT bin[0..1GeV];#frac{pT jet - pT truth}{pT jet}",600,-3,3);
  for(int i = 1 ; i < 200 ; ++i) {
    gauss_forpt[i] = (TH1F*)gauss_forpt[0]->Clone();
    sprintf(name,"pT bin[%d..%dGeV]",i,i+1);
    gauss_forpt[i]->SetTitle(name);
    gauss_forptcorr[i] = (TH1F*)gauss_forptcorr[0]->Clone();
    sprintf(name,"corrected jet pT bin[%d..%dGeV]",i,i+1);
    gauss_forptcorr[i]->SetTitle(name);
  }
  
  //loop over all fit-events
  for ( std::vector<TData*>::iterator i = data->begin(); 
        i != data->end() ; ++i )  {
    TData* jg = *i;
    double etjetcor = jg->GetParametrizedMess();
    if( jg->GetType() != TypeGammaJet ) continue;
    if(jg->GetMess()[0]>0 && jg->GetMess()[0]<200)
      gauss_forpt[(int)jg->GetMess()[0]]->Fill( (jg->GetMess()[0]-jg->GetTruth())/jg->GetMess()[0], jg->GetWeight() );
    if(etjetcor>0 && etjetcor<200)
      gauss_forptcorr[(int)etjetcor]->Fill( (etjetcor-jg->GetTruth())/etjetcor, jg->GetWeight() );
  }

  double edge;
  TText * text = new TText();
  text->SetTextSize(0.03);
  text->SetTextColor(2);
  //TF1 * f = new TF1("gauss_step",gauss_step,-10,10,5);
  for(int i = 0 ; i < 200 ; ++i) {
    //edge = 1.0-20./(((double)i)+0.5);
    //f->SetParameters(-1.,2.0,3.0, edge, 0.01);
    //f->FixParameter(3, edge);
    //f->FixParameter(4, 0.01);
    
    gauss_forpt[i]->Fit("gaus","LLQNO","");
    TF1 * f = (TF1*)gROOT->GetFunction("gaus")->Clone();
    //cout<<"bin "<<i
    //    <<": mean="<<f->GetParameter(0)
    //	<<", sigma="<<f->GetParameter(1)
    //	<<endl;
    
    //double mean = f->GetParameter(1);
    //double meanerror = f->GetParError(1);
    //double width = f->GetParameter(2);
    
    //gauss_forpt[i]->Fit("gauss_step","LLQNO","");
    //cout<<"bin "<<i<<": mean="<<f->GetParameter(0)
    //	<<", sigma="<<f->GetParameter(1)
    //	<<", height="<<f->GetParameter(2)
    //	<<", edge("<<edge<<")="<<f->GetParameter(3)
    //	<<", width-edge="<<f->GetParameter(4)
    //	<<endl;
    
    gauss_forpt[i]->Draw("h");
    f->SetLineColor(2);
    f->Draw("same");
    sprintf(name,"mean %f",f->GetParameter(1));
    text->DrawText(1.4,0.7*gauss_forpt[i]->GetMaximum(),name);
    //func->Draw("same");
    c1->Draw();
    delete f;
/*
    TF1 *g=0;
    //Fit1D(gauss_forptcorr[i],g);
    gauss_forptcorr[i]->Fit("gaus","LLQNO","");
    g = (TF1*)gROOT->GetFunction("gaus")->Clone();
    gauss_forptcorr[i]->Draw("h");
    g->SetLineColor(4);
    g->Draw("same");
    c1->Draw();
    delete g;
*/
  }
  //delete f;
  ps.Close();
  
  for(int i = 0 ; i < 200 ; ++i){
    delete gauss_forpt[i];  
    delete gauss_forptcorr[i];  
  }  
  delete name;
}

void TControlPlots::Fit1D(TH1F* hist, TF1* result)
{
   hist->Fit("gaus","LLQNO","");
   result = (TF1*)gROOT->GetFunction("gaus")->Clone();
}


