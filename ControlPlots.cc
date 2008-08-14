#include "ControlPlots.h"

#include <iostream>
using namespace std;

#include "TCanvas.h"
#include "TDirectory.h"
#include "TH1I.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TPostScript.h"
#include "TString.h"
#include "TROOT.h"

#include "CalibData.h"
#include "CalibMath.h"
#include "Parameters.h"


//---------------------------------------------------------------
//  Constructor
// 
//  Objects of 'TControlPlots' created by this constructor
//  can create control plots via the 'MakeControlPlots...()'
//  methods from several TData objects. The output is in
//  .ps or both ps. and .root format.
//  
//  Parameters:
//  data           TData objects to create controlplots from
//  par            Parameters
//  outputFormat   Specify what output format to use:
//                 0  .ps and .root (default)
//                 1  .ps
//---------------------------------------------------------------
TControlPlots::TControlPlots(const std::vector<TData*> *data, TParameters *par, int outputFormat)
  : _data(data), _par(par), _outFile( outputFormat==0 ? new TFile("controlplots.root","RECREATE","Cal calib control plots") : 0 )
{ 
  ptRatioName[0] = "p^{jet}_{T}/ E_{T}^{#gamma}";
  ptRatioName[1] = "p_{T}^{cor. jet}/E_{T}^{#gamma}";
  ptRatioName[2] = "p_{T}^{jet}/p_{T}^{cor. jet}";

  controlQuantityName[0] = "mean"; 
  controlQuantityName[1] = "standard deviation"; 
  controlQuantityName[2] = "mean of Gauss fit"; 
  controlQuantityName[3] = "width of Gauss fit"; 
  controlQuantityName[4] = "median"; 
  controlQuantityName[5] = "#chi^{2} / n.d.f."; 
  controlQuantityName[6] = "probability"; 
  controlQuantityName[7] = "quantiles"; 

  SetGStyle();

  if( outputFormat == 1 ) _outputROOT = false;
  else _outputROOT = true;
}


TControlPlots::~TControlPlots()
{
  if( _outFile !=0 )
    {
      if( _outFile->IsOpen() ) _outFile->Close();
      delete _outFile;
    }
}




//---------------------------------------------------------------
//   Tower control plots
//---------------------------------------------------------------
void TControlPlots::MakeControlPlotsTowers()
{
  std::vector<TObject*> objToBeWritten;
  std::vector<TData*>::const_iterator data_it,it;

  TCanvas * const c1 = new TCanvas("c1","",600,600);
  TPostScript * const ps = new TPostScript("controlplotsTowers.ps",111);


  //book objects
  TLatex *latex = new TLatex();
  latex->SetTextSize(0.03);
  latex->SetTextFont(42);

  TMeasurement * testmess = new TMeasurement();

  char name[100];
  int markerColor[5] = { 1,4,2,3,7 };
  int markerStyle[5] = { 8,1,1,1,1 };

  TH2F * constants = new TH2F("hCalibConstants","Calibration constants vs. E_{T} and #eta-bin with f_{em} = OUF = 0;E_{T} [GeV];#eta-bin",100,0.5,100.,_par->GetEtaGranularity(),1,_par->GetEtaGranularity()); 
   
  for (int eta=0; eta<_par->GetEtaGranularity();++eta) // Loop over eta bins
    {
      for (int phi=0; phi<_par->GetPhiGranularity();++phi) // Loop over phi bins
	{
	  int i = _par->GetBin(eta,phi);
	  double * val = _par->GetTowerParRef(i);

	  TH1F * plot[5]; 
	  sprintf(name, "h%d_eta%d_phi%d",i,eta+1,phi+1);
	  plot[0] = new TH1F(name,";uncorrected tower E_{T} [GeV];tower k-factor",100,0.5,100.);    
	  sprintf(name, "htt%d_eta%d_phi%d",i,eta+1,phi+1);
	  plot[1] = new TH1F(name,";uncorrected tower E_{T} [GeV];tower k-factor",100,0.5,100.);    
	  sprintf(name, "hgj%d_eta%d_phi%d",i,eta+1,phi+1);
	  plot[2] = new TH1F(name,";uncorrected tower E_{T} [GeV];tower k-factor",100,0.5,100.);    
	  sprintf(name, "htc%d_eta%d_phi%d",i,eta+1,phi+1);
	  plot[3] = new TH1F(name,";uncorrected tower E_{T} [GeV];tower k-factor",100,0.5,100.); 
	  sprintf(name, "hjj%d_eta%d_phi%d",i,eta+1,phi+1);
	  plot[4] = new TH1F(name,";uncorrected tower E_{T} [GeV];tower k-factor",100,0.5,100.); 
   
	  TH1F * plot_had[5]; 
	  sprintf(name, "hhad%d_eta%d_phi%d",i,eta+1,phi+1);
	  plot_had[0] = new TH1F(name,";uncorrected hadronic fraction of tower E_{T} [GeV];tower k-factor",100,0.5,100.);    
	  sprintf(name, "hhadtt%d_eta%d_phi%d",i,eta+1,phi+1);
	  plot_had[1] = new TH1F(name,";uncorrected hadronic fraction of ttower E_{T} [GeV];tower k-factor",100,0.5,100.);    
	  sprintf(name, "hhadgj%d_eta%d_phi%d",i,eta+1,phi+1);
	  plot_had[2] = new TH1F(name,";uncorrected hadronic fraction of ttower E_{T} [GeV];tower k-factor",100,0.5,100.);    
	  sprintf(name, "hhadtc%d_eta%d_phi%d",i,eta+1,phi+1);
	  plot_had[3] = new TH1F(name,";uncorrected hadronic fraction of ttower E_{T} [GeV];tower k-factor",100,0.5,100.);    
	  sprintf(name, "hhadjj%d_eta%d_phi%d",i,eta+1,phi+1);
	  plot_had[4] = new TH1F(name,";uncorrected hadronic fraction of ttower E_{T} [GeV];tower k-factor",100,0.5,100.);    

	  TH1F * k[5];
	  sprintf(name, "k%d_eta%d_phi%d",i,eta+1,phi+1);
	  k[0] = new TH1F(name,";uncorrected tower E_{T} [GeV];tower k-factor",100,0.5,100.);    
	  sprintf(name, "ktt%d_eta%d_phi%d",i,eta+1,phi+1);
	  k[1] = new TH1F(name,";uncorrected tower E_{T} [GeV];tower k-factor",100,0.5,100.);    
	  sprintf(name, "kgj%d_eta%d_phi%d",i,eta+1,phi+1);
	  k[2] = new TH1F(name,";uncorrected tower E_{T} [GeV];tower k-factor",100,0.5,100.);    
	  sprintf(name, "ktc%d_eta%d_phi%d",i,eta+1,phi+1);
	  k[3] = new TH1F(name,";uncorrected tower E_{T} [GeV];tower k-factor",100,0.5,100.); 
	  sprintf(name, "kjj%d_eta%d_phi%d",i,eta+1,phi+1);
	  k[4] = new TH1F(name,";uncorrected tower E_{T} [GeV];tower k-factor",100,0.5,100.); 
   
	  TH1F * em[5];
	  sprintf(name, "em%d_eta%d_phi%d",i,eta+1,phi+1);
	  em[0] = new TH1F(name,"",100,0.5,100.);    
	  sprintf(name, "emtt%d_eta%d_phi%d",i,eta+1,phi+1);
	  em[1] = new TH1F(name,"",100,0.5,100.);    
	  sprintf(name, "emgj%d_eta%d_phi%d",i,eta+1,phi+1);
	  em[2] = new TH1F(name,"",100,0.5,100.);    
	  sprintf(name, "emtc%d_eta%d_phi%d",i,eta+1,phi+1);
	  em[3] = new TH1F(name,"",100,0.5,100.);    
	  sprintf(name, "emjj%d_eta%d_phi%d",i,eta+1,phi+1);
	  em[4] = new TH1F(name,"",100,0.5,100.);    

	  TH1F * had[5];
	  sprintf(name, "had%d_eta%d_phi%d",i,eta+1,phi+1);
	  had[0] = new TH1F(name,"",100,0.5,100.);    
	  sprintf(name, "hadtt%d_eta%d_phi%d",i,eta+1,phi+1);
	  had[1] = new TH1F(name,"",100,0.5,100.);    
	  sprintf(name, "hadgj%d_eta%d_phi%d",i,eta+1,phi+1);
	  had[2] = new TH1F(name,"",100,0.5,100.);    
	  sprintf(name, "hadtc%d_eta%d_phi%d",i,eta+1,phi+1);
	  had[3] = new TH1F(name,"",100,0.5,100.);    
	  sprintf(name, "hadjj%d_eta%d_phi%d",i,eta+1,phi+1);
	  had[4] = new TH1F(name,"",100,0.5,100.);    

	  TH1F * au[5];
	  sprintf(name, "au%d_eta%d_phi%d",i,eta+1,phi+1);
	  au[0] = new TH1F(name,"",100,0.5,100.);    
	  sprintf(name, "autt%d_eta%d_phi%d",i,eta+1,phi+1);
	  au[1] = new TH1F(name,"",100,0.5,100.);    
	  sprintf(name, "augj%d_eta%d_phi%d",i,eta+1,phi+1);
	  au[2] = new TH1F(name,"",100,0.5,100.);    
	  sprintf(name, "autc%d_eta%d_phi%d",i,eta+1,phi+1);
	  au[3] = new TH1F(name,"",100,0.5,100.);    
	  sprintf(name, "aujj%d_eta%d_phi%d",i,eta+1,phi+1);
	  au[4] = new TH1F(name,"",100,0.5,100.);    

	  TH1F * khad[5];
	  sprintf(name, "khad%d_eta%d_phi%d",i,eta+1,phi+1);
	  khad[0] = new TH1F(name,";uncorrected tower E_{T} [GeV];tower k-factor",100,0.5,100.);    
	  sprintf(name, "khadtt%d_eta%d_phi%d",i,eta+1,phi+1);
	  khad[1] = new TH1F(name,";uncorrected tower E_{T} [GeV];tower k-factor",100,0.5,100.);    
	  sprintf(name, "khadgj%d_eta%d_phi%d",i,eta+1,phi+1);
	  khad[2] = new TH1F(name,";uncorrected tower E_{T} [GeV];tower k-factor",100,0.5,100.);    
	  sprintf(name, "khadtc%d_eta%d_phi%d",i,eta+1,phi+1);
	  khad[3] = new TH1F(name,";uncorrected tower E_{T} [GeV];tower k-factor",100,0.5,100.);    
	  sprintf(name, "khadjj%d_eta%d_phi%d",i,eta+1,phi+1);
	  khad[4] = new TH1F(name,";uncorrected tower E_{T} [GeV];tower k-factor",100,0.5,100.);    

	  TH1F * et_vs_had[5];
	  sprintf(name, "et_vs_had%d_eta%d_phi%d",i,eta+1,phi+1);
	  et_vs_had[0] = new TH1F(name,"",100,0.5,100.);    
	  sprintf(name, "et_vs_hadtt%d_eta%d_phi%d",i,eta+1,phi+1);
	  et_vs_had[1] = new TH1F(name,"",100,0.5,100.);    
	  sprintf(name, "et_vs_hadgj%d_eta%d_phi%d",i,eta+1,phi+1);
	  et_vs_had[2] = new TH1F(name,"",100,0.5,100.);    
	  sprintf(name, "et_vs_hadtc%d_eta%d_phi%d",i,eta+1,phi+1);
	  et_vs_had[3] = new TH1F(name,"",100,0.5,100.);    
	  sprintf(name, "et_vs_hadjj%d_eta%d_phi%d",i,eta+1,phi+1);
	  et_vs_had[4] = new TH1F(name,"",100,0.5,100.);    

	  TH1F * em_vs_had[5];
	  sprintf(name, "em_vs_had%d_eta%d_phi%d",i,eta+1,phi+1);
	  em_vs_had[0] = new TH1F(name,"",100,0.5,100.);    
	  sprintf(name, "em_vs_hadtt%d_eta%d_phi%d",i,eta+1,phi+1);
	  em_vs_had[1] = new TH1F(name,"",100,0.5,100.);    
	  sprintf(name, "em_vs_hadgj%d_eta%d_phi%d",i,eta+1,phi+1);
	  em_vs_had[2] = new TH1F(name,"",100,0.5,100.);    
	  sprintf(name, "em_vs_hadtc%d_eta%d_phi%d",i,eta+1,phi+1);
	  em_vs_had[3] = new TH1F(name,"",100,0.5,100.);    
	  sprintf(name, "em_vs_hadjj%d_eta%d_phi%d",i,eta+1,phi+1);
	  em_vs_had[4] = new TH1F(name,"",100,0.5,100.);    

	  TH1F * au_vs_had[5];
	  sprintf(name, "au_vs_had%d_eta%d_phi%d",i,eta+1,phi+1);
	  au_vs_had[0] = new TH1F(name,"",100,0.5,100.);    
	  sprintf(name, "au_vs_hadtt%d_eta%d_phi%d",i,eta+1,phi+1);
	  au_vs_had[1] = new TH1F(name,"",100,0.5,100.);    
	  sprintf(name, "au_vs_hadgj%d_eta%d_phi%d",i,eta+1,phi+1);
	  au_vs_had[2] = new TH1F(name,"",100,0.5,100.);    
	  sprintf(name, "au_vs_hadtc%d_eta%d_phi%d",i,eta+1,phi+1);
	  au_vs_had[3] = new TH1F(name,"",100,0.5,100.);    
	  sprintf(name, "au_vs_hadjj%d_eta%d_phi%d",i,eta+1,phi+1);
	  au_vs_had[4] = new TH1F(name,"",100,0.5,100.);    

	  TH1F * norm_vs_had[5];
	  sprintf(name, "norm_vs_had%d_eta%d_phi%d",i,eta+1,phi+1);
	  norm_vs_had[0] = new TH1F(name,"",100,0.5,100.);    
	  sprintf(name, "norm_vs_hadtt%d_eta%d_phi%d",i,eta+1,phi+1);
	  norm_vs_had[1] = new TH1F(name,"",100,0.5,100.);    
	  sprintf(name, "norm_vs_hadgj%d_eta%d_phi%d",i,eta+1,phi+1);
	  norm_vs_had[2] = new TH1F(name,"",100,0.5,100.);    
	  sprintf(name, "norm_vs_hadtc%d_eta%d_phi%d",i,eta+1,phi+1);
	  norm_vs_had[3] = new TH1F(name,"",100,0.5,100.);    
	  sprintf(name, "norm_vs_hadjj%d_eta%d_phi%d",i,eta+1,phi+1);
	  norm_vs_had[4] = new TH1F(name,"",100,0.5,100.);    

	  TH1F * norm[5];
	  sprintf(name, "hnorm%d",i);
	  norm[0] = new TH1F(name,"",100,0.5,100.);    
	  sprintf(name, "hnormtt%d",i);
	  norm[1] = new TH1F(name,"",100,0.5,100.);    
	  sprintf(name, "hnormgj%d",i);
	  norm[2] = new TH1F(name,"",100,0.5,100.);    
	  sprintf(name, "hnormtc%d",i);
	  norm[3] = new TH1F(name,"",100,0.5,100.);    
	  sprintf(name, "hnormjj%d",i);
	  norm[4] = new TH1F(name,"",100,0.5,100.);    

	  TH1F * chi2[4];
	  sprintf(name, "hgj_chi2_%d_eta%d_phi%d",i,eta+1,phi+1);
	  chi2[0] = new TH1F(name,";#chi^{2};N",100,0.0,50.);    
	  sprintf(name, "htt_chi2_%d_eta%d_phi%d",i,eta+1,phi+1);
	  chi2[1] = new TH1F(name,";#chi^{2};N",100,0.0,50.);    
	  sprintf(name, "htc_chi2_%d_eta%d_phi%d",i,eta+1,phi+1);
	  chi2[2] = new TH1F(name,";#chi^{2};N",100,0.0,50.);    
	  sprintf(name, "hjj_chi2_%d_eta%d_phi%d",i,eta+1,phi+1);
	  chi2[3] = new TH1F(name,";#chi^{2};N",100,0.0,50.);    

	  sprintf(name, "h2dgj%d_eta%d_phi%d",i,eta+1,phi+1);
	  TH2F * plot2dgj = new TH2F(name,";uncorrected tower E_{T} [GeV];tower k-factor",100,0.5,100.,100,0.0,5.0);    
	  sprintf(name, "hdgj_weight%d_eta%d_phi%d",i,eta+1,phi+1);
	  TH1F * plotgj_weight = new TH1F(name,";uncorrected tower E_{T} [GeV];tower k-factor",100,0.5,100.);

	  sprintf(name, "norm_weight%d_eta%d_phi%d",i,eta+1,phi+1);
	  TH1F * norm_weight = new TH1F(name,";uncorrected tower E_{T} [GeV];tower k-factor",100,0.5,100.);    

	  sprintf(name, "h2dtt%d_eta%d_phi%d",i,eta+1,phi+1);
	  TH2F * plot2dtt = new TH2F(name,";uncorrected tower E_{T} [GeV];tower k-factor",100,0.5,100.,100,0.0,5.0);    

	  sprintf(name, "h2dtc%d_eta%d_phi%d",i,eta+1,phi+1);
	  TH2F * plot2dtc = new TH2F(name,";uncorrected tower E_{T} [GeV];tower k-factor",100,0.5,100.,100,0.0,5.0);    
	  sprintf(name, "h2djj%d_eta%d_phi%d",i,eta+1,phi+1);
	  TH2F * plot2djj = new TH2F(name,";uncorrected tower E_{T} [GeV];tower k-factor",100,0.5,100.,100,0.0,5.0);    


	  // Fill histos
	  double Jet, JetCorr;

	  double maxTowerET;	      
	  int thisIndexJet;

	  double mess, error;

	  data_it = _data->begin();
	  for (; data_it != _data->end();++data_it) // 1. loop over all fit events
	    {
	      //if one fit event is composed of multiple towers, than loop over all
	      const std::vector<TData*>& data_ref = (*data_it)->GetRef();

	      // Get sum of corrected tower Et
	      Jet = 0.;
	      for(it = data_ref.begin(); it != data_ref.end(); ++it) //1. loop over towers
		{
		  Jet += (*it)->GetParametrizedMess();
		}

	      // Find tower with max Et
	      maxTowerET = 0.;	      
	      thisIndexJet = 0;
	      for (it = data_ref.begin(); it != data_ref.end(); ++it) //2. loop over towers
		{
		  TMeasurement * m = (*it)->GetMess();
		  int thisIndex = (*it)->GetIndex();
		  if ( m->pt > maxTowerET ) 
		    {
		      thisIndexJet = thisIndex;
		      maxTowerET = m->pt;
		    }
		  if ( thisIndex != i ) continue; //tower belongs to a wrong bin

		  if (m->pt != 0.0)
		    {
		      norm[0]->Fill(m->pt );
		      em[ 0]->Fill( m->pt, m->EMF );
		      had[0]->Fill( m->pt, m->HadF );
		      au[ 0]->Fill( m->pt, m->OutF );
		      if ( (*data_it)->GetType()==GammaJet ) 
			{
			  norm[2]->Fill(m->pt );
			  em[ 2]->Fill( m->pt, m->EMF );
			  had[2]->Fill( m->pt, m->HadF );
			  au[ 2]->Fill( m->pt, m->OutF );
			}  
		      else if ( (*data_it)->GetType()==TrackTower ) 
			{
			  norm[1]->Fill(m->pt );
			  em[ 1]->Fill( m->pt, m->EMF );
			  had[1]->Fill( m->pt, m->HadF );
			  au[ 1]->Fill( m->pt, m->OutF );
			}  
		      else if ( (*data_it)->GetType()==TrackCluster)
			{
			  norm[3]->Fill(m->pt );
			  em[ 3]->Fill( m->pt, m->EMF );
			  had[3]->Fill( m->pt, m->HadF );
			  au[ 3]->Fill( m->pt, m->OutF );
			}  
		      else if ( (*data_it)->GetType()==PtBalance)
			{
			  norm[4]->Fill(m->pt );
			  em[ 4]->Fill( m->pt, m->EMF );
			  had[4]->Fill( m->pt, m->HadF );
			  au[ 4]->Fill( m->pt, m->OutF );
			}  
		    }

		  if (m->HadF!=0.0)
		    {
		      norm_vs_had[0]->Fill( m->HadF );
		      et_vs_had[  0]->Fill( m->HadF, m->pt );
		      em_vs_had[  0]->Fill( m->HadF, m->EMF );
		      au_vs_had[  0]->Fill( m->HadF, m->OutF );
		      if ((*data_it)->GetType()==GammaJet) 
			{
			  norm_vs_had[2]->Fill( m->HadF );
			  et_vs_had[  2]->Fill( m->HadF, m->pt );
			  em_vs_had[  2]->Fill( m->HadF, m->EMF );
			  au_vs_had[  2]->Fill( m->HadF, m->OutF );
			}  
		      else if ((*data_it)->GetType()==TrackTower) 
			{
			  norm_vs_had[1]->Fill( m->HadF );
			  et_vs_had[  1]->Fill( m->HadF, m->pt );
			  em_vs_had[  1]->Fill( m->HadF, m->EMF );
			  au_vs_had[  1]->Fill( m->HadF, m->OutF );
			}  
		      else if ((*data_it)->GetType()==TrackCluster) 
			{
			  norm_vs_had[3]->Fill( m->HadF );
			  et_vs_had[  3]->Fill( m->HadF, m->pt );
			  em_vs_had[  3]->Fill( m->HadF, m->EMF );
			  au_vs_had[  3]->Fill( m->HadF, m->OutF );
			}  
		      else if ((*data_it)->GetType()==PtBalance) 
			{
			  norm_vs_had[4]->Fill( m->HadF );
			  et_vs_had[  4]->Fill( m->HadF, m->pt );
			  em_vs_had[  4]->Fill( m->HadF, m->EMF );
			  au_vs_had[  4]->Fill( m->HadF, m->OutF );
			}  
		    }
		} // End of 2. loop over towers
	    } // End of 1. loop over all fit events

	  for(int a = 0; a < 5; a++)
	    {
	      em[a]->Divide(norm[a]);
	      had[a]->Divide(norm[a]);
	      au[a]->Divide(norm[a]);
	      em_vs_had[a]->Divide(norm_vs_had[a]);
	      et_vs_had[a]->Divide(norm_vs_had[a]);
	      au_vs_had[a]->Divide(norm_vs_had[a]);
	    }
	  
	  for (data_it = _data->begin(); data_it != _data->end();++data_it) // 2. loop over all fit-events
	    {
	      //if one fit event is composed of multiple towers, than loop over all
	      mess = 0.;
	      error=0.;

	      const std::vector<TData*>& data_ref = (*data_it)->GetRef();
	      JetCorr = (*data_it)->GetParametrizedMess(); // Corrected jet Pt
	      Jet = 0.;		// Sum of corrected tower Et
	      for (it =data_ref.begin(); it!=data_ref.end(); ++it) // 3. loop over towers
		{
		  Jet += (*it)->GetParametrizedMess();
		}

	      maxTowerET = 0.;
	      thisIndexJet = 0;
	      for (it =data_ref.begin(); it!=data_ref.end(); ++it) // 4. loop over towers
		{
		  double m = (*it)->GetMess()->pt;
		  double mhad = (*it)->GetMess()->HadF;
		  double t   = (*it)->GetTruth()*(Jet/JetCorr);
		  double tmp = (*it)->GetParametrizedMess();
		  mess  += tmp;
		  error +=  (*it)->GetParametrizedErr(&tmp);
		  int thisIndex = (*it)->GetIndex();
		  if (m>maxTowerET)
		    {
		      thisIndexJet = thisIndex;
		      maxTowerET = m;
		    }
		  if (thisIndex !=i) continue; //tower belongs to a wrong bin

		  if (mhad!=0.0)
		    {
		      testmess->pt = m;
		      testmess->EMF = em[0]->GetBinContent(em[0]->GetXaxis()->FindBin(m));
		      testmess->HadF = had[0]->GetBinContent(had[0]->GetXaxis()->FindBin(m));
		      testmess->OutF = au[0]->GetBinContent(au[0]->GetXaxis()->FindBin(m));
		      plot_had[0]->Fill( mhad, t/m );
		      khad[0]->Fill(     mhad, _par->plot_parametrization(testmess,val) );
		      if ((*data_it)->GetType()==GammaJet)
			{
			  plot_had[2]->Fill( mhad, t/m );
			  khad[2]->Fill(     mhad, _par->plot_parametrization(testmess,val) );
			}  
		      else if ((*data_it)->GetType()==TrackTower)
			{
			  khad[1]->Fill(     mhad, _par->plot_parametrization(testmess,val) );
			  plot_had[1]->Fill( mhad, t/m );
			}  
		      else if ((*data_it)->GetType()==TrackCluster)
			{
			  plot_had[3]->Fill( mhad, t/m );
			  khad[3]->Fill(     mhad, _par->plot_parametrization(testmess,val) );
			}  
		      else if ((*data_it)->GetType()==PtBalance)
			{
			  plot_had[4]->Fill( mhad, t/m );
			  khad[4]->Fill(     mhad, _par->plot_parametrization(testmess,val) );
			}  
		    }

		  if (m!=0.0)
		    {
		      plot[0]->Fill( m, t/m );
	    
		      testmess->pt = m;
		      testmess->EMF = em[0]->GetBinContent(em[0]->GetXaxis()->FindBin(m));
		      testmess->HadF = had[0]->GetBinContent(had[0]->GetXaxis()->FindBin(m));
		      testmess->OutF = au[0]->GetBinContent(au[0]->GetXaxis()->FindBin(m));
		      k[0]->Fill( m, _par->plot_parametrization(testmess,val) );

		      if ((*data_it)->GetType()==GammaJet)
			{
			  plot[2]->Fill( m, t/m );
			  plot2dgj->Fill( m, t/m );
			  double weight = t/Jet;
			  norm_weight->Fill( m, weight );
			  plotgj_weight->Fill( m, t/m );
			  testmess->EMF = em[2]->GetBinContent(em[2]->GetXaxis()->FindBin(m));
			  testmess->HadF = had[2]->GetBinContent(had[2]->GetXaxis()->FindBin(m));
			  testmess->OutF = au[2]->GetBinContent(au[2]->GetXaxis()->FindBin(m));
			  k[2]->Fill( m, _par->plot_parametrization(testmess,val) );
			}  
		      else if ((*data_it)->GetType()==TrackTower)
			{
			  plot[1]->Fill( m, t/m );
			  plot2dtt->Fill( m, t/m );
			  testmess->EMF = em[1]->GetBinContent(em[1]->GetXaxis()->FindBin(m));
			  testmess->HadF = had[1]->GetBinContent(had[1]->GetXaxis()->FindBin(m));
			  testmess->OutF = au[1]->GetBinContent(au[1]->GetXaxis()->FindBin(m));
			  k[1]->Fill( m, _par->plot_parametrization(testmess,val) );
			}  
		      else if ((*data_it)->GetType()==TrackCluster) 
			{
			  plot[3]->Fill( m, t/m );
			  plot2dtc->Fill( m, t/m );
			  testmess->EMF = em[3]->GetBinContent(em[3]->GetXaxis()->FindBin(m));
			  testmess->HadF = had[3]->GetBinContent(had[3]->GetXaxis()->FindBin(m));
			  testmess->OutF = au[3]->GetBinContent(au[3]->GetXaxis()->FindBin(m));
			  k[3]->Fill( m, _par->plot_parametrization(testmess,val) );
			}  
		      else if ((*data_it)->GetType()==PtBalance) 
			{
			  plot[4]->Fill( m, t/m );
			  plot2djj->Fill( m, t/m );
			  testmess->EMF = em[3]->GetBinContent(em[3]->GetXaxis()->FindBin(m));
			  testmess->HadF = had[3]->GetBinContent(had[3]->GetXaxis()->FindBin(m));
			  testmess->OutF = au[3]->GetBinContent(au[3]->GetXaxis()->FindBin(m));
			  k[4]->Fill( m, _par->plot_parametrization(testmess,val) );
			}  
		    }
		} // End of 4. loop over towers


	      if (thisIndexJet!=i) continue; //event (jet or tower) belongs to a wrong bin

	      switch ( (*data_it)->GetType())
		{
		case TrackTower://track-tower
		  if (mess != 0.0) chi2[0]->Fill( (*data_it)->chi2()/(*data_it)->GetWeight() ); break;
		case GammaJet://gamma-jet
		  if (mess != 0.0) chi2[1]->Fill( (*data_it)->chi2()/(*data_it)->GetWeight() ); break;
		case TrackCluster://track-cluster
		  if (mess != 0.0) chi2[2]->Fill( (*data_it)->chi2()/(*data_it)->GetWeight() ); 
		  break;
		case PtBalance://jet-jet
		  if (mess != 0.0) chi2[3]->Fill( (*data_it)->chi2()/(*data_it)->GetWeight() ); break;
		default:// do nothing
		  break;
		}
	    }// End of 2. loop over all fit-events

	  if (norm[0]->GetEntries()==0) continue;

	  for(int a = 0; a < 5; a++)
	    {
	      plot[a]->Divide(norm[a]);
	      plot_had[a]->Divide(norm_vs_had[a]);
	      k[a]->Divide(norm[a]);
	      khad[a]->Divide(norm_vs_had[a]);
	    }
	  plotgj_weight->Divide(norm_weight);

	  for (int j=0; j<5; ++j)//number of diff samples, i.e. gj, tt, tc, jj,...
	    {
	      for (int b=1; b<=100; ++b) //number of bins
		{
		  if (norm[j]->GetBinContent(b)>0)
		    {
		      plot[j]->SetBinError(  b, 1./sqrt(
							norm[j]->GetBinContent(b)     //stat
							)*plot[j]->GetBinContent(b) );
		      plot_had[j]->SetBinError(  b, 1./sqrt(
							    norm_vs_had[j]->GetBinContent(b)     //stat
							    )*plot_had[j]->GetBinContent(b) );
		    }			      
		}
	    }


	  for(int a = 0; a < 5; a++)
	    {
	      plot[a]->SetMarkerStyle(markerStyle[a]);
	      plot[a]->SetMarkerColor(markerColor[a]);
	      plot[a]->SetLineColor(markerColor[a]);
	      plot[a]->GetYaxis()->SetRangeUser(0.,6.);
	      if( a == 0 ) plot[a]->Draw("PE");
	      else plot[a]->Draw("PE SAME");

	      k[a]->SetLineColor(markerColor[a]);
	      k[a]->SetLineWidth(3);
	      k[a]->GetYaxis()->SetRangeUser(0.,6.);
	      k[a]->Draw("L SAME");

	      objToBeWritten.push_back(plot[a]->Clone());
	      objToBeWritten.push_back(k[a]->Clone());
 	    }
	  latex->DrawLatex( 66,3.8,"Average");
	  latex->DrawLatex( 66,3.5,"#color[2]{Gamma-Jet}");
	  latex->DrawLatex( 66,3.2,"#color[3]{Track-Cluster}");
	  latex->DrawLatex( 66,2.9,"#color[4]{Track-Tower}");
	  latex->DrawLatex( 66,2.6,"#color[7]{p_{T}-Balance}");

	  c1->Draw(); 
	  ps->NewPage();

    
	  sprintf(name, "h_k_hadonly_chi2_%d_eta%d_phi%d",i,eta+1,phi+1);
	  TH1F * khadonly = new TH1F(name,"f_{em} = OUF = 0 [GeV];uncorrected hadronic tower E_{T} [GeV];tower k-factor",100,0.5,100.);    

	  sprintf(name, "h_k_Efrac02_chi2_%d_eta%d_phi%d",i,eta+1,phi+1);
	  TH1F * kEfrac02 = new TH1F(name,"f_{em} = 0.2, OUF = 0 [GeV];uncorrected hadronic tower E_{T} [GeV];tower k-factork",100,0.5,100.);    

	  sprintf(name, "h_k_Efrac05_chi2_%d_eta%d_phi%d",i,eta+1,phi+1);
	  TH1F * kEfrac05 = new TH1F(name,"f_{em} = 0.5, OUF = 0 [GeV];uncorrected hadronic tower E_{T} [GeV];tower k-factork",100,0.5,100.);    

	  testmess->OutF = 0.0;
	  for (int b=1; b<=100; ++b)
	    {
	      testmess->pt = (double)b; // -> EMFrac = 0.0
	      testmess->EMF = 0.0;
	      testmess->HadF = (double)b;
	      khadonly->SetBinContent(khadonly->GetXaxis()->FindBin(b), 
				      _par->plot_parametrization(testmess,val) );
	      constants->SetBinContent(constants->GetXaxis()->FindBin(b),constants->GetYaxis()->FindBin(eta+1),
				       _par->plot_parametrization(testmess,val));
	      testmess->EMF = (double)b*0.2; // -> EMFrac = 0.25
	      testmess->HadF = (double)b*0.8;
	      kEfrac02->SetBinContent(kEfrac02->GetXaxis()->FindBin(b), 
				      _par->plot_parametrization(testmess,val) );       
	      testmess->EMF = (double)b*0.5; // -> EMFrac = 1.0
	      testmess->HadF = (double)b*0.5;
	      kEfrac05->SetBinContent(kEfrac05->GetXaxis()->FindBin(b), 
				      _par->plot_parametrization(testmess,val) );       
	    }

	  khadonly->GetYaxis()->SetRangeUser(0.9,1.1);
	  khadonly->Draw("l"); 
	  c1->Draw(); 
	  ps->NewPage();
	  objToBeWritten.push_back(khadonly->Clone());

	  kEfrac02->GetYaxis()->SetRangeUser(0.9,1.1);
	  kEfrac02->Draw("l"); 
	  c1->Draw(); 
	  ps->NewPage();
	  objToBeWritten.push_back(kEfrac02->Clone());

	  kEfrac05->GetYaxis()->SetRangeUser(0.9,1.1);
	  kEfrac05->Draw("l"); 
	  c1->Draw(); 
	  ps->NewPage();
	  objToBeWritten.push_back(kEfrac05->Clone());

	  for(int a = 0; a < 5; a++)
	    {
	      plot_had[a]->SetMarkerStyle(markerStyle[a]);
	      plot_had[a]->SetMarkerColor(markerColor[a]);
	      plot_had[a]->SetLineColor(markerColor[a]);
	      plot_had[a]->GetYaxis()->SetRangeUser(0.,6.);
	      if( a == 0 ) plot_had[a]->Draw("PE");
	      else plot_had[a]->Draw("PE SAME");

	      khad[a]->SetLineColor(markerColor[a]);
	      khad[a]->SetLineWidth(3);
	      khad[a]->GetYaxis()->SetRangeUser(0.,6.);
	      khad[a]->Draw("L SAME");

	      objToBeWritten.push_back(plot_had[a]->Clone());
	      objToBeWritten.push_back(khad[a]->Clone());
	    }
	  latex->DrawLatex( 66,3.8,"Average");
	  latex->DrawLatex( 66,3.5,"#color[2]{Gamma-Jet}");
	  latex->DrawLatex( 66,3.2,"#color[3]{Track-Cluster}");
	  latex->DrawLatex( 66,2.9,"#color[4]{Track-Tower}");
	  latex->DrawLatex( 66,2.6,"#color[7]{p_{T}-Balance}");

	  c1->Draw(); 
	  ps->NewPage();


	  plot2dgj->SetLineColor(2);
	  plot2dtc->SetLineColor(3);
	  plot2dtt->SetLineColor(4);
	  plot2djj->SetLineColor(7);	      
	  plot2dgj->Draw("BOX");	      
	  plot2dtt->Draw("BOX,same");	      
	  plot2dtc->Draw("BOX,same");	      
	  plot2djj->Draw("BOX,same");	      

	  objToBeWritten.push_back(plot2dgj->Clone());
	  objToBeWritten.push_back(plot2dtt->Clone());
	  objToBeWritten.push_back(plot2dtc->Clone());
	  objToBeWritten.push_back(plot2djj->Clone());
			      
	  latex->DrawLatex( 0.66*(plot2dgj->GetXaxis()->GetXmax()-plot2dgj->GetXaxis()->GetXmin()),
			    0.5*(plot2dgj->GetYaxis()->GetXmax()-plot2dgj->GetYaxis()->GetXmin()),
			    "#color[2]{Gamma-Jet}");
	  latex->DrawLatex( 0.66*(plot2dgj->GetXaxis()->GetXmax()-plot2dgj->GetXaxis()->GetXmin()),
			    0.45*(plot2dgj->GetYaxis()->GetXmax()-plot2dgj->GetYaxis()->GetXmin()),
			    "#color[3]{TrackCluster}");
	  latex->DrawLatex( 0.66*(plot2dgj->GetXaxis()->GetXmax()-plot2dgj->GetXaxis()->GetXmin()),
			    0.4*(plot2dgj->GetYaxis()->GetXmax()-plot2dgj->GetYaxis()->GetXmin()),
			    "#color[4]{TrackTower}");
	  latex->DrawLatex( 0.66*(plot2dgj->GetXaxis()->GetXmax()-plot2dgj->GetXaxis()->GetXmin()),
			    0.35*(plot2dgj->GetYaxis()->GetXmax()-plot2dgj->GetYaxis()->GetXmin()),
			    "#color[7]{p_{T}-Balance}");
	  c1->Draw(); 
	  ps->NewPage();
      

	  //show overflow in the last bin
	  double maximum = 0.0;
	  for (int k=0; k<4; ++k)
	    {
	      if (chi2[k]->GetMaximum()>maximum)
		maximum = chi2[k]->GetMaximum(); 
	      chi2[k]->SetBinContent(100,chi2[k]->GetBinContent(100)+chi2[k]->GetBinContent(101));
	      chi2[k]->SetBinContent(1,chi2[k]->GetBinContent(0)+chi2[k]->GetBinContent(1));
	    }
	  for(int k = 0; k < 4; k++)
	    {
	      chi2[k]->SetMaximum( maximum+sqrt(maximum) );
	      chi2[k]->SetLineColor(markerColor[k+1]);
	      if( k == 0 ) chi2[k]->Draw("h");
	      else chi2[k]->Draw("h,same");
	      
	      objToBeWritten.push_back(chi2[k]->Clone());
	    }

	  latex->DrawLatex( 0.66*(chi2[0]->GetXaxis()->GetXmax()-chi2[0]->GetXaxis()->GetXmin()),
			    0.65*(chi2[0]->GetMaximum()-chi2[0]->GetMinimum()),
			    "#color[2]{Gamma-Jet}");
	  latex->DrawLatex( 0.66*(chi2[0]->GetXaxis()->GetXmax()-chi2[0]->GetXaxis()->GetXmin()),
			    0.6*(chi2[0]->GetMaximum()-chi2[0]->GetMinimum()),
			    "#color[1]{TrackTower}");
	  latex->DrawLatex( 0.66*(chi2[0]->GetXaxis()->GetXmax()-chi2[0]->GetXaxis()->GetXmin()),
			    0.55*(chi2[0]->GetMaximum()-chi2[0]->GetMinimum()),
			    "#color[3]{TrackCluster}");
	  latex->DrawLatex( 0.66*(chi2[0]->GetXaxis()->GetXmax()-chi2[0]->GetXaxis()->GetXmin()),
			    0.5*(chi2[0]->GetMaximum()-chi2[0]->GetMinimum()),
			    "#color[7]{p_{T}-Balance}");
	  c1->Draw(); 
	  ps->NewPage();

	  for(int a = 0; a < 5; a++)
	    {
	      delete plot[a];
	      delete plot_had[a];
	      delete k[a];
	      delete em[a];
	      delete had[a];
	      delete au[a];
	      delete khad[a];
	      delete et_vs_had[a];
	      delete em_vs_had[a];
	      delete au_vs_had[a];
	      delete norm_vs_had[a];
	      delete norm[a];
	      if(a < 4) delete chi2[a];
	    }
	  delete plot2dgj;
	  delete plotgj_weight;
	  delete norm_weight;
	  delete plot2dtt;
	  delete plot2dtc;
	  delete plot2djj;
	  delete khadonly;
	  delete kEfrac02; 
	  delete kEfrac05; 

	} // End of loop over phi bins
    } // End of loop over eta bins

  gStyle->SetPalette(1);
  constants->Draw("COLZ"); 
  c1->Draw(); 
  ps->NewPage();
  objToBeWritten.push_back(constants);


  // Distribution of chi2 summands (normalized residuals)
  TH1F *h_chi2 = new TH1F("h_chi2","Scaled residuals z^{2} = 1/w #chi^{2};f(z^{2});dN / df(z^{2})",50,0,200);
  h_chi2->SetLineColor(1);
  h_chi2->SetLineStyle(3);
  objToBeWritten.push_back(h_chi2);

  // Distribution of Cauchy scaled normalized residuals
  TH1F *h_cauchy = static_cast<TH1F*>(h_chi2->Clone("h_cauchy"));
  h_cauchy->SetLineColor(2);
  h_cauchy->SetLineStyle(2);
  objToBeWritten.push_back(h_cauchy);

  // Distribution of Huber scaled normalized residuals
  TH1F *h_huber = static_cast<TH1F*>(h_chi2->Clone("h_huber"));
  h_huber->SetLineColor(4);
  h_huber->SetLineStyle(1);
  objToBeWritten.push_back(h_huber);

  // Cauchy-scaled versus no scaling
  TH2F *h_none_cauchy = new TH2F("h_none_cauchy","Residuals z^{2};z^{2};f(z^{2})",50,0,10,50,0,10);
  h_none_cauchy->SetLineColor(2);
  objToBeWritten.push_back(h_none_cauchy);

  // Huber-scaled versus no scaling
  TH2F *h_none_huber = static_cast<TH2F*>(h_none_cauchy->Clone("h_none_huber"));
  h_none_huber->SetLineColor(4);
  objToBeWritten.push_back(h_none_huber);

  for(  std::vector<TData*>::const_iterator it = _data->begin();  it < _data->end();  ++it )
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
  ps->NewPage();


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


  if( _outputROOT ) WriteToRootFile(objToBeWritten, "Towers");
    
  ps->Close();

  delete h_chi2;
  delete h_cauchy;
  delete h_huber;
  delete l_res;
  delete h_none_cauchy;
  delete h_none_huber;
  delete line1;
  delete l_res2;
  delete constants;
  delete latex;
  delete testmess;
  delete c1;
  delete ps;
  objToBeWritten.clear();
}



//---------------------------------------------------------------
//   Gamma-Jet Control Histograms
//---------------------------------------------------------------
void TControlPlots::MakeControlPlotsGammaJet()
{
  std::vector<TObject*> objToBeWritten;


  TCanvas * const c1 = new TCanvas("c1","",600,600);
  TCanvas * const c2 = new TCanvas("c2","",600,600);
  c2->Divide(2,2);

  TPostScript * const ps = new TPostScript("controlplotsGammaJet.ps",111);


  //book hists
  char name[100];
  char title[200];
  int etLimit[4] = { 10, 35, 90, 300 };

  // Number towers in jet
  TH1I* towerinjet[4];
  towerinjet[0]= new TH1I("hNTowInJet0","#gamma-jet;# tower in jet",40,0,40);
  objToBeWritten.push_back(towerinjet[0]);
  for(int i =1;i<4;++i)
    {
      sprintf(name,"hNTowInJet%i",i);
      sprintf(title,"#gamma-jet,  %i < E_{T}^{#gamma} < %i GeV;# tower in jet",etLimit[i-1],etLimit[i]);
      towerinjet[i] = (TH1I*) towerinjet[0]->Clone(name);
      towerinjet[i]->SetTitle(title);
      objToBeWritten.push_back(towerinjet[i]);
    }


  // Energy in rings around jet axis
  int nRings = 9;
  TH1F* leadToNext[21];  //0+1 rings,2-12 eta dependence ,13+14 Relrings, 15-20 Pt depend (not rel)
  leadToNext[0] = new  TH1F("hLeadToNext0","#gamma-jet,  Et in ring (raw = red) / raw jet Et;#Delta R;Ring Et / Jet Et",nRings,0,0.6);
  objToBeWritten.push_back(leadToNext[0]);
  leadToNext[1] = (TH1F*) leadToNext[0]->Clone("hLeadToNext1");   //raw tower
  objToBeWritten.push_back(leadToNext[1]);
  leadToNext[2] = new  TH1F("hLeadToNext2","#gamma-jet,  Sum(ET in ring) / Sum(ET in ring(raw));#Delta R;Ring Et / Ring Et(raw)",nRings,0,0.6);
  objToBeWritten.push_back(leadToNext[2]);
  for(int i=3;i<=14;++i)
    {
      sprintf(name,"hLeadToNext%i",i);
      leadToNext[i] = (TH1F*) leadToNext[2]->Clone(name);
      float etaMin = 0.5*(i-3);
      float etaMax = etaMin + 0.5;
      sprintf(title,"Sum*(ET in ring)/Sum(ET in ring(raw)) %.1f < abs#eta < %.1f",etaMin,etaMax);
      leadToNext[i]->SetTitle(title);
      objToBeWritten.push_back(leadToNext[i]);
    }
  leadToNext[13]->SetTitle("Rel Et (divided by leading tower Et) in ring (raw = red);#Delta R;Rel Et in Ring");
  leadToNext[14]->SetTitle("Rel Et (divided by leading tower Et) in ring (raw = red);#Delta R;Rel Et in Ring");
  for(int i=14;i<21;++i)
    {
      sprintf(name,"hLeadToNext%i",i);
      leadToNext[i] = (TH1F*) leadToNext[2]->Clone(name);
      int etReg = (i-14)/2;
      sprintf(title,"#gamma-jet, E_{T} in ring / raw Jet E{T}   %i < P_{T}^{Jet} < %i GeV (raw = red);#Delta R;Ring E_{T} / raw jet E_{T}",etLimit[etReg],etLimit[etReg+1]);
      leadToNext[i]->SetTitle(title);
      objToBeWritten.push_back(leadToNext[i]);
    }

  double ringsSum[11][nRings];
  double ringsRawSum[11][nRings];
    for(int b=0; b < 11; ++b)
      {
	for(int a=0; a < nRings; ++a)
	  {
	    ringsSum[b][a]=0;
	    ringsRawSum[b][a]=0;
	  }
      }


  // Hits in (eta,phi)
  TH2F* EtaPhiMap = new TH2F("hEtaPhiHitMap","#gamma-jet,  #eta-#Phi hit map;#eta;#Phi",200,-5,5,128,-3.2,3.2);
  objToBeWritten.push_back(EtaPhiMap);


  // Response vs tower Et
  TH2F* respvstet[4];
  respvstet[0] = new TH2F("hResVsTowerEt0","#gamma-jet,  Response vs. Tower Et;raw tower Et;Correction factor",100,0,100,100,-4,10);
  objToBeWritten.push_back(respvstet[0]);
  for(int i = 1 ; i < 4 ; ++i)
    {
      sprintf(name,"hResVsTowerEt%i",i);
      respvstet[i] = (TH2F*)respvstet[0]->Clone(name);
      objToBeWritten.push_back(respvstet[i]);
    }
  respvstet[1]->SetTitle("#gamma-jet,  Response vs. Tower Et (leading Tower);raw tower Et;Correction factor");
  respvstet[2]->SetTitle("#gamma-jet,  Response vs. Tower Et (neighbouring leading Tower);raw tower Et;Correction factor");
  respvstet[3]->SetTitle("#gamma-jet,  Response vs. Tower Et (2 towers next to leading Tower);raw tower Et;Correction factor");



  // The following 2D histograms contain different Pt ratios:
  //   ptjet/etgamma
  //   ptjetcorr/etgamma 
  //   ptjet/ptjetcorr
  // versus different quantities:
  //   eta (for all Etgamma and binned in different Etgamma)
  //   uncorrected jet pt (for all eta and binned in different eta)
  //   Etgamma = true jet pt
  //   log(Etgamma)
  //   jet energy
  //   electromagnetic fraction emf
  //
  // For arrays of size 12:
  //   base index i = 0,3,6,9 indicates Etgamma or eta bin:
  //     0: All Et / all eta
  //     3: 10 < Et < 35   /  |eta| < 1.4
  //     6: 35 < Et < 90   /  1.4 < |eta| < 3.0
  //     9: 90 < Et < 300  /  3.0 < |eta|
  //   sub-index i+0, i+1, i+2 indicates plotted Pt ratio:
  //     base+0: ptjet/etgamma
  //     base+1: ptjetcorr/etgamma
  //     base+2: ptjet/ptjetcorr
  //
  // For arrays of size 3:
  //   index indicates plotted Pt ratio:
  //     0: ptjet/etgamma
  //     1: ptjetcorr/etgamma
  //     2: ptjet/ptjetcorr
  TH2F* heta[12];
  heta[0] = new TH2F("heta0","#gamma-jet;#eta",400,-5,5,100,0,4);
  objToBeWritten.push_back(heta[0]);

  TH2F* hpt_uncorr[12];
  hpt_uncorr[0] = new TH2F("hpt_uncorr0","#gamma-jet;p^{jet}_{T} [GeV]",400,0,400,100,0,4);
  objToBeWritten.push_back(hpt_uncorr[0]);

  for(int i = 1 ; i < 12 ; ++i)
    {
      sprintf(name,"heta%i",i);
      heta[i] = (TH2F*)heta[0]->Clone(name);
      heta[i]->Sumw2();

      sprintf(name,"hpt_uncorr%i",i);
      hpt_uncorr[i] = (TH2F*) hpt_uncorr[0]->Clone(name);
      hpt_uncorr[i]->Sumw2();

      if( i > 2 )
	{
	  int etReg = (i-3)/3;
	  sprintf(title,"#gamma-jet,  %i < E_{T}^{#gamma} < %i GeV",etLimit[etReg],etLimit[etReg+1]);
	  heta[i]->SetTitle(title);

	  if( etReg == 0 ) sprintf(title,"#gamma-jet,  |#eta| < 1.4");
	  if( etReg == 1 ) sprintf(title,"#gamma-jet,  1.4 < |#eta| < 3.0");
	  if( etReg == 2 ) sprintf(title,"#gamma-jet,  3.0 < |#eta|");

	  hpt_uncorr[i]->SetTitle(title);
	}
      objToBeWritten.push_back(heta[i]);
    }

  TH2F* hpt[3];
  hpt[0] = new TH2F("hpt0","#gamma-jet;E^{#gamma}_{T} [GeV]",400,0,400,100,0,4);
  hpt[1] = (TH2F*)hpt[0]->Clone("hpt1");
  hpt[2] = (TH2F*)hpt[0]->Clone("hpt2");

  TH2F* henergy[3];
  henergy[0] = new TH2F("henergy0","#gamma-jet;E^{jet} [GeV]",600,0,600,100,0,4);
  henergy[1] = (TH2F*)henergy[0]->Clone("henergy1");
  henergy[2] = (TH2F*)henergy[0]->Clone("henergy2");
  
  TH2F* hemf[3];
  hemf[0] = new TH2F("hemf0","#gamma-jet;electromagnetic fraction f_{em}",100,0,1,100,0,4);
  hemf[1] = (TH2F*)hemf[0]->Clone("hemf1");
  hemf[2] = (TH2F*)hemf[0]->Clone("hemf2");

  double bins[101];
  for(int i = 0; i < 101 ; ++i)
    {
      bins[i] = pow(10,(i+32)/40.0);
    }
  TH2F* hptlog[3];
  hptlog[0] = new TH2F("hptlog0","#gamma-jet;E^{#gamma}_{T} [GeV]",100,bins,100,0,4);
  hptlog[1] = (TH2F*) hptlog[0]->Clone("hptlog1");
  hptlog[2] = (TH2F*) hptlog[0]->Clone("hptlog2");

  for(int i = 0; i < 3; i++)
    {
      hpt[i]->Sumw2();
      henergy[i]->Sumw2();
      hemf[i]->Sumw2();
      hptlog[i]->Sumw2();

      objToBeWritten.push_back(hpt[i]);
      objToBeWritten.push_back(henergy[i]);
      objToBeWritten.push_back(hemf[i]);
      objToBeWritten.push_back(hptlog[i]);
    }


  // Number of events per pt, eta and emf bin
  // to test TCaliber::FlattenSpectra() 
  TH1F* hptGamma[3];
  TH1F* hptGammaW[3];
  hptGamma[0]  = new TH1F("hptGamma","#gamma-jet;p_{T} [GeV]",100,0,400);
  hptGammaW[0] = new TH1F("hptGammaW","#gamma-jet;p_{T} [GeV]",100,0,400);
  hptGamma[1]  = new TH1F("hetaGamma","#gamma-jet;#eta",100,-5,5);
  hptGammaW[1] = new TH1F("hetaGammaW","#gamma-jet;#eta",100,-5,5);
  hptGamma[2]  = new TH1F("hemfGamma","#gamma-jet;f_{em}",100,0,1);
  hptGammaW[2] = new TH1F("hemfGammaW","#gamma-jet;f_{em}",100,0,1);
  for(int i = 0; i < 3; i++)
    {
      objToBeWritten.push_back(hptGamma[i]);
      objToBeWritten.push_back(hptGammaW[i]);
    }

  // emf vs pt before and after weighting
  TH2F* hptGamma2D = new TH2F("hGamma2D","#gamma-jet,  w/o weights;p_{T} [GeV];f_{em}",500,0,500,100,0,1);
  objToBeWritten.push_back(hptGamma2D);
  TH2F* hptGamma2DW = new TH2F("hGamma2DW","#gamma-jet,  with weights;p_{T} [GeV];f_{em}",500,0,500,100,0,1);
  objToBeWritten.push_back(hptGamma2DW);



  // Fill histos

  //loop over all fit-events
  for ( std::vector<TData*>::const_iterator i = _data->begin() ; i != _data->end() ; ++i )
    {
      TData* jg = *i;
      if( jg->GetType() != GammaJet ) continue;

      double etjet = jg->GetMess()->pt;
      double energyjet = jg->GetMess()->pt;  //Et -> E see below
      double etjetcor = jg->GetParametrizedMess();
      double etajet = jg->GetMess()->eta;
      double phijet = jg->GetMess()->phi;
      int noTower=0;
      double maxTowerET=0;
      double maxTowerEnergy=0;
      double maxTowerETraw=0;
      int thisIndexJet=0;
      double maxTowerEta=0;
      double maxTowerPhi=0;

      // first loop over tower
      const std::vector<TData*>& data_ref = jg->GetRef();
      for(std::vector<TData*>::const_iterator it = data_ref.begin();it != data_ref.end(); ++it)
	{
	  TMeasurement * m = (*it)->GetMess();
	  double  pm = (*it)->GetParametrizedMess();
	  respvstet[0]->Fill(m->pt,pm/m->pt);
	  int thisIndex = (*it)->GetIndex();
	  if (pm>maxTowerET) {
	    thisIndexJet = thisIndex;
	    maxTowerET = pm;
	    maxTowerETraw = m->pt;
	    maxTowerEta = m->eta;
	    maxTowerPhi = m->phi;
	    maxTowerEnergy = m->E;
	  }

	  EtaPhiMap->Fill(m->eta,m->phi);

	  ++noTower;
	} // end of first loop over tower
      energyjet *= maxTowerEnergy / maxTowerETraw;
      respvstet[1]->Fill(maxTowerETraw,maxTowerET/maxTowerETraw);
      double rings[nRings];
      double ringsRaw[nRings];
      for(int a=0; a < nRings; ++a)
	{
	  rings[a]=0;
	  ringsRaw[a]=0;
	}

      // second loop over tower
      for(std::vector<TData*>::const_iterator it = data_ref.begin();it != data_ref.end(); ++it)
	{
	  //next to leading towers:
	  TMeasurement * m  = (*it)->GetMess();
	  double  pm = (*it)->GetParametrizedMess();
	  if(m->pt != 0)
	    {
	      for(int a=0; a < nRings; ++a)
		{
		  double DeltaR = sqrt((deltaPhi(phijet,m->phi) * deltaPhi(phijet,m->phi)) + ((etajet - m->eta) * (etajet - m->eta)));
		  if(DeltaR <= ((0.6 / nRings) * a))
		    {
		      rings[a] += pm / etjet;
		      ringsRaw[a] += m->pt / etjet;
		    }
		}	    
	    
	      int index = (*it)->GetIndex();
	      if((abs(index - thisIndexJet) == 1)      )// || (abs(index - thisIndexJet + phi_granularity) < 2)  || (abs(index - thisIndexJet - phi_granularity) < 2))
		respvstet[2]->Fill(m->pt, (*it)->GetParametrizedMess()/m->pt);
	    
	      if(abs(index - thisIndexJet) == 2)  // different with proper Phi granularity
		respvstet[3]->Fill(m->pt, (*it)->GetParametrizedMess()/m->pt);
	    }
	}  // end of second loop over tower

      // loop over rings
      for(int a=0; a < nRings; ++a)
	{
	  if(a == 0)
	    {
	      leadToNext[0]->Fill(a * (0.6 / nRings),rings[0]);
	      leadToNext[1]->Fill(a * (0.6 / nRings),ringsRaw[0]);
	      leadToNext[13]->Fill(a * (0.6 / nRings),1);
	      leadToNext[14]->Fill(a * (0.6 / nRings),1);
	      if((etjet > 10) && (etjet < 35))
		{
		  leadToNext[15]->Fill(a * (0.6 / nRings),rings[0]);
		  leadToNext[16]->Fill(a * (0.6 / nRings),ringsRaw[0]);
		}
	      if((etjet > 35) && (etjet < 90))
		{
		  leadToNext[17]->Fill(a * (0.6 / nRings),rings[0]);
		  leadToNext[18]->Fill(a * (0.6 / nRings),ringsRaw[0]);
		}
	      if((etjet > 90) && (etjet < 300))
		{
		  leadToNext[19]->Fill(a * (0.6 / nRings),rings[0]);
		  leadToNext[20]->Fill(a * (0.6 / nRings),ringsRaw[0]);
		}
	    }
	  else
	    {
	      leadToNext[0]->Fill(a * (0.6 / nRings),rings[a]-rings[a-1]);
	      leadToNext[1]->Fill(a * (0.6 / nRings),ringsRaw[a]-ringsRaw[a-1]);
	      if((etjet > 10) && (etjet < 35))
		{
		  leadToNext[15]->Fill(a * (0.6 / nRings),rings[a]-rings[a-1]);
		  leadToNext[16]->Fill(a * (0.6 / nRings),ringsRaw[a]-ringsRaw[a-1]);
		}
	      if((etjet > 35) && (etjet < 90))
		{
		  leadToNext[17]->Fill(a * (0.6 / nRings),rings[a]-rings[a-1]);
		  leadToNext[18]->Fill(a * (0.6 / nRings),ringsRaw[a]-ringsRaw[a-1]);
		}
	      if((etjet > 90) && (etjet < 300))
		{
		  leadToNext[19]->Fill(a * (0.6 / nRings),rings[a]-rings[a-1]);
		  leadToNext[20]->Fill(a * (0.6 / nRings),ringsRaw[a]-ringsRaw[a-1]);
		}
	      if(rings[0]!=0)
		{
		  leadToNext[13]->Fill(a * (0.6 / nRings),(rings[a]-rings[a-1])/rings[0]);
		  leadToNext[14]->Fill(a * (0.6 / nRings),(ringsRaw[a]-ringsRaw[a-1])/ringsRaw[0]);
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
	} // end of loop over rings
      towerinjet[0]->Fill(noTower);


      heta[0]->Fill(etajet,etjet/ jg->GetTruth(),jg->GetWeight());
      heta[1]->Fill(etajet,etjetcor/ jg->GetTruth(),jg->GetWeight());
      heta[2]->Fill(etajet,etjet/etjetcor,jg->GetWeight());
      if (jg->GetTruth() > 10 && jg->GetTruth() < 35)
	{
	  heta[3]->Fill(etajet,etjet/ jg->GetTruth(),jg->GetWeight());
	  heta[4]->Fill(etajet,etjetcor/ jg->GetTruth(),jg->GetWeight());
	  heta[5]->Fill(etajet,etjet/etjetcor,jg->GetWeight());
	  towerinjet[1]->Fill(noTower);
	}
      else if (jg->GetTruth() > 35 && jg->GetTruth() < 90)
	{
	  heta[6]->Fill(etajet,etjet/ jg->GetTruth(),jg->GetWeight());
	  heta[7]->Fill(etajet,etjetcor/ jg->GetTruth(),jg->GetWeight());
	  heta[8]->Fill(etajet,etjet/etjetcor,jg->GetWeight());
	  towerinjet[2]->Fill(noTower);
	}
      else if (jg->GetTruth() > 90 && jg->GetTruth() < 300)
	{
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
      hptlog[2]->Fill(jg->GetTruth(),etjet/etjetcor,jg->GetWeight());

      hpt_uncorr[0]->Fill(etjet,etjet/ jg->GetTruth(),jg->GetWeight());
      hpt_uncorr[1]->Fill(etjet,etjetcor/jg->GetTruth(),jg->GetWeight());
      hpt_uncorr[2]->Fill(etjet,etjet/etjetcor,jg->GetWeight());    
      if(fabs(etajet) < 1.4)
	{
	  hpt_uncorr[3]->Fill(etjet,etjet/ jg->GetTruth(),jg->GetWeight());
	  hpt_uncorr[4]->Fill(etjet,etjetcor/jg->GetTruth(),jg->GetWeight());
	  hpt_uncorr[5]->Fill(etjet,etjet/etjetcor,jg->GetWeight());    
	}
      else if(fabs(etajet) > 1.4  && fabs(etajet) < 3.0)
	{
	  hpt_uncorr[6]->Fill(etjet,etjet/ jg->GetTruth(),jg->GetWeight());
	  hpt_uncorr[7]->Fill(etjet,etjetcor/jg->GetTruth(),jg->GetWeight());
	  hpt_uncorr[8]->Fill(etjet,etjet/etjetcor,jg->GetWeight());    
	}
      else if(fabs(etajet) > 3.0)
	{
	  hpt_uncorr[9]->Fill(etjet,etjet/ jg->GetTruth(),jg->GetWeight());
	  hpt_uncorr[10]->Fill(etjet,etjetcor/jg->GetTruth(),jg->GetWeight());
	  hpt_uncorr[11]->Fill(etjet,etjet/etjetcor,jg->GetWeight());    
	}

      /*
	hptlog[0]->Fill(jg->GetMess()->pt,etjet/ jg->GetTruth(),jg->GetWeight());
	hptlog[1]->Fill(jg->GetMess()->pt,etjetcor/jg->GetTruth(),jg->GetWeight());
	hptlog[2]->Fill(etjetcor,etjet/etjetcor,jg->GetWeight());
      */
      //em fraction plots     
      double em = 0;
      double had = 0;
      for(std::vector<TData*>::const_iterator t = jg->GetRef().begin(); t != jg->GetRef().end(); ++t)
	{
	  TData* tt = *t;
	  em  += tt->GetMess()->EMF;
	  had += tt->GetMess()->HadF;
	  had += tt->GetMess()->OutF;
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
    } // end of loop over all fit-events

  // Fill energy in rings around jet axis
  for(int b=0; b< 11; ++b)
    {
      for(int a=0; a < nRings; ++a)
	{
	  if(ringsRawSum[b][a] != 0)
	    leadToNext[b+2]->Fill(a * (0.6 / nRings),ringsSum[b][a]/ringsRawSum[b][a]);
	}
    }

  // Draw histos
  c1->cd();

  // Number of towers in jets
  for(int i=0;i<4;++i)
    {
      towerinjet[i]->Draw();
      c1->Draw();
      ps->NewPage();
    }
  /*
    for(int i=0;i<4;++i)
    {
    respvstet[i]->Draw("Box");
    c1->Draw();
    ps->NewPage();
    }
  */

  // Energy in rings around jet axis
  leadToNext[0]->Draw();
  leadToNext[1]->SetLineColor(2);
  leadToNext[1]->Draw("same");
  c1->Draw();
  ps->NewPage();

  for(int i=0; i<3;++i)
    {
      leadToNext[15 + (2*i)]->Draw();
      leadToNext[16 + (2*i)]->SetLineColor(2);
      leadToNext[16 + (2*i)]->Draw("same");
      c1->Draw();
      ps->NewPage();
    }

  EtaPhiMap->Draw("Col");//,Palette");
  c1->Draw();
  ps->NewPage();
  EtaPhiMap->Draw("Palette");
  c1->Draw();
  ps->NewPage();




  // From the above specified 2D histograms containing different
  // Pt ratios versus different quantities, projections along
  // the x-axis are made per x bin using TControlPlots::Fit2D(...).
  // Some properties of these projected distributions, in the
  // following called 'control quantities', are plotted:
  //   0: Mean value
  //   1: Standard deviation
  //   2: Mean of Gauss fit
  //   3: Width of Gauss fit
  //   4: Median 
  //   5: chi2 / n.d.f.
  //   6: Probability of Gauss fit
  //   7: Quantiles Q0.9 / (Q0.9 - 1)
  //
  // Also, the projected distributions are plotted for some
  // example x-bins together with the Gauss fit


  // Control quantities from the Pt ratio vs eta plots.
  // First dimension 12:
  //   base index i = 0,3,6,9 indicates Etgamma bin:
  //     0: All Et / all eta
  //     3: 10 < Et < 35
  //     6: 35 < Et < 90
  //     9: 90 < Et < 300
  //   sub-index i+0, i+1, i+2 indicates plotted Pt ratio:
  //     base+0: ptjet/etgamma
  //     base+1: ptjetcorr/etgamma
  //     base+2: ptjet/ptjetcorr
  // Second dimension 8 is the above specified control
  // quantity
  TH1F* hists_eta[12][8];

  // Projected distributions for 4 example eta-bins
  // and the Gauss fits
  TH1F* gp_eta[12][4];	
  TF1* gf_eta[12][4];
  
  int markerStyle[3] = { 20, 22, 21 };
  int markerColor[3] = { 1, 2, 4 };

  TLegend* leg = new TLegend(0.7,0.68,0.96,0.9);
  leg->SetFillColor(0);
  leg->AddEntry(heta[0],ptRatioName[0],"p");
  leg->AddEntry(heta[2],ptRatioName[2],"p");
  leg->AddEntry(heta[1],ptRatioName[1],"p");

  for(int i = 0 ; i < 12 ; i += 3) // Loop over Etgamma bins
    {
      for(int a = 0; a < 3; a++) // Loop over Pt ratios
	{
	  heta[i+a]->SetMinimum(0.5);
	  heta[i+a]->SetMaximum(1.2);
	  heta[i+a]->SetMarkerStyle(markerStyle[a]);
	  heta[i+a]->SetMarkerColor(markerColor[a]);
	  heta[i+a]->SetLineColor(markerColor[a]);
	}

      // Do projections and determine control quantities
      Fit2D(heta[i],hists_eta[i],gp_eta[i], gf_eta[i]);
      Fit2D(heta[i+1],hists_eta[i+1],gp_eta[i+1], gf_eta[i+1]);
      Fit2D(heta[i+2],hists_eta[i+2],gp_eta[i+2], gf_eta[i+2]); 

      // Set axis ranges
      for(int a = 0; a<3;++a)
	{
	  for(int b = 0 ; b < 4 ; ++b)
	    {
	      hists_eta[a+i][b]->SetMinimum(0.2);
	      hists_eta[a+i][b]->SetMaximum(1.8);
	      ++b;
	      hists_eta[a+i][b]->SetMinimum(0.0);
	      hists_eta[a+i][b]->SetMaximum(0.5);
	    }
	}

      // Draw gaussplots for example eta bins
      // on multi-canvas
      for(int a = 0; a < 3; a++) // Loop over ptratios
	{
	  for(int b = 0; b < 3; b++) // Loop over example eta bins
	    {
	      // Find eta bin of gaussplot
	      int bin = 0;
	      if(  b == 0  )  bin = int(heta[i+a]->GetNbinsX()/6);
	      else if(  b == 1  )  bin = int(heta[i+a]->GetNbinsX()/3);
	      else if(  b == 2  )  bin = int(heta[i+a]->GetNbinsX()/2);
	      float etaMin = heta[i+a]->GetXaxis()->GetBinLowEdge(bin);
	      float etaMax = etaMin + heta[i+a]->GetXaxis()->GetBinWidth(bin);

	      // Set title according to energy and eta bin
	      if( i == 0 )  sprintf(title,"#gamma-jet, %.2f < #eta < %.2f",etaMin,etaMax);
	      else sprintf(title,"#gamma-jet, %i < E_{T}^{#gamma} < %i GeV, %.2f < #eta < %.2f",
			   etLimit[int(i/3)-1],etLimit[int(i/3)],etaMin,etaMax);
	      gp_eta[i+a][b]->SetTitle(title);
	      gp_eta[i+a][b]->SetXTitle(ptRatioName[a]);

	      // Set style and line color according to ptRatioName
	      gp_eta[i+a][b]->SetMarkerStyle(markerStyle[a]);
	      gp_eta[i+a][b]->SetMarkerColor(markerColor[a]);
	      gp_eta[i+a][b]->SetLineColor(markerColor[a]);
	      gf_eta[i+a][b]->SetLineColor(markerColor[a]);

	      // Plot gaussplots_eta
	      c2->cd(1+b);
	      gp_eta[i+a][b]->Draw();
	      gf_eta[i+a][b]->Draw("same");

	      objToBeWritten.push_back(gp_eta[i+a][b]);
	      objToBeWritten.push_back(gf_eta[i+a][b]);
	    } // End of loop over example eta bins
	  c2->Draw();
	  ps->NewPage();
	} // End of loop over ptratios

      // Draw control quantities
      c1->cd();
      for(int j = 0 ; j < 8 ; ++j) // Loop over control quantities
	{
	  hists_eta[i][j]->Draw("P");
	  hists_eta[i][j]->SetStats(0);
	  hists_eta[i+1][j]->Draw("P SAME");
	  hists_eta[i+2][j]->Draw("P SAME");
	  leg->Draw("SAME");
	  c1->SetGrid();
	  c1->Draw();   
	  ps->NewPage(); 

	  objToBeWritten.push_back(hists_eta[i][j]);
	  objToBeWritten.push_back(hists_eta[i+1][j]);
	  objToBeWritten.push_back(hists_eta[i+2][j]);
	}
    } // End loop over Etgamma bins


  //Comparison of control quantities in Etgamma bins
  //and Pt ratio bins
  TLegend* leg2 = new TLegend(0.7,0.68,0.96,0.9);
  leg2->SetFillColor(0);
  for(int i = 0; i < 3; i++)
    {
      sprintf(title,"%i < E_{T}^{#gamma} < %i GeV",etLimit[i],etLimit[i+1]);
      leg2->AddEntry(hists_eta[i][0],title,"P");
    }

  for(int a = 3 ; a < 5 ; ++a) // Loop over ptratio bins
    {
      for(int b = 0 ; b < 8 ; ++b) // Loop over control quantites
	{
	  for(int c = 0; c < 7; c += 3) // Loop over etgamma bins
	    {
	      TH1F *h = static_cast<TH1F*>(hists_eta[a+c][b]->Clone("h"));
	      h->SetTitle("#gamma-jet,  " + ptRatioName[a-3] + " " + controlQuantityName[b]);
	      h->SetMarkerStyle(markerStyle[c/3]);
	      h->SetMarkerColor(markerColor[c/3]);
	      h->SetLineColor(markerColor[c/3]);
	      if( c == 0 ) h->Draw("P");
	      else h->Draw("P SAME");
	      h->SetStats(0);
	    }
	  leg2->Draw("SAME");
	  c1->SetGrid();
	  c1->Draw();   
	  ps->NewPage();
	}
    }
  delete leg2;



  // Control quantities from the Pt ratio vs uncorrected
  // jet Pt plots.
  // First dimension 12:
  //   base index i = 0,3,6,9 indicates eta bin:
  //     0: All eta
  //     3: |eta| < 1.4
  //     6: 1.4 < |eta| < 3.0
  //     9: 3.0 < |eta|
  //   sub-index i+0, i+1, i+2 indicates plotted Pt ratio:
  //     base+0: ptjet/etgamma
  //     base+1: ptjetcorr/etgamma
  //     base+2: ptjet/ptjetcorr
  // Second dimension 8 is the above specified control
  // quantity
  TH1F* hists_ptuncorr[12][8];

  // Projected distributions for 4 example eta-bins
  // and the Gauss fits
  TH1F* gp_ptuncorr[12][4];  
  TF1* gf_ptuncorr[12][4];

  for(int i = 0 ; i < 12 ; i+= 3)  // Loop over eta-bins
    {
      for(int a = 0; a < 3; a++) // Loop over Pt ratios
	{
	  hpt_uncorr[i+a]->SetMinimum(0.5);
	  hpt_uncorr[i+a]->SetMaximum(1.2);
	  hpt_uncorr[i+a]->SetMarkerStyle(markerStyle[a]);
	  hpt_uncorr[i+a]->SetMarkerColor(markerColor[a]);
	  hpt_uncorr[i+a]->SetLineColor(markerColor[a]);
	}

      // Do projections and determine control quantities
      Fit2D(hpt_uncorr[i],hists_ptuncorr[i],gp_ptuncorr[i], gf_ptuncorr[i]);
      Fit2D(hpt_uncorr[i+1],hists_ptuncorr[i+1],gp_ptuncorr[i+1], gf_ptuncorr[i+1]);
      Fit2D(hpt_uncorr[i+2],hists_ptuncorr[i+2],gp_ptuncorr[i+2], gf_ptuncorr[i+2]); 
      for(int a = 0; a<3;++a)
	{
	  for(int b = 0 ; b < 4 ; ++b)
	    {
	      hists_ptuncorr[a+i][b]->SetMinimum(0.2);
	      hists_ptuncorr[a+i][b]->SetMaximum(1.8);
	      ++b;
	      hists_ptuncorr[a+i][b]->SetMinimum(0.0);
	      hists_ptuncorr[a+i][b]->SetMaximum(0.5);
	    }
	}

      // Draw gaussplots for example pt bins
      // on multi-canvas
      for(int a = 0; a < 3; a++) // Loop over ptratios
	{
	  for(int b = 0; b < 3; b++) // Loop over example pt bins
	    {
	      // Find pt bin of gaussplot
	      int bin = 0;
	      if(  b == 0  )  bin = int(hpt_uncorr[i+a]->GetNbinsX()/6);
	      else if(  b == 1  )  bin = int(hpt_uncorr[i+a]->GetNbinsX()/3);
	      else if(  b == 2  )  bin = int(hpt_uncorr[i+a]->GetNbinsX()/2);
	      float min = hpt_uncorr[i+a]->GetXaxis()->GetBinLowEdge(bin);
	      float max = min + hpt_uncorr[i+a]->GetXaxis()->GetBinWidth(bin);

	      // Set title according to energy and pt bin
	      if( i == 0 )      sprintf(title,"#gamma-jet, %.1f < p^{jet}_{T} < %.1f GeV",min,max);
	      else if( i == 3 ) sprintf(title,
					"#gamma-jet, |#eta| < 1.4, %.1f < p^{jet}_{T} < %.1f GeV",
					min,max);
	      else if( i == 6 ) sprintf(title,
					"#gamma-jet, 1.4 < |#eta| < 3.0, %.1f < p^{jet}_{T} < %.1f GeV",
					min,max);
	      else if( i == 9 ) sprintf(title,
					"#gamma-jet, 3.0 < |#eta|, %.1f < p^{jet}_{T} < %.1f GeV",
					min,max);
	      gp_ptuncorr[i+a][b]->SetTitle(title);
	      gp_ptuncorr[i+a][b]->SetXTitle(ptRatioName[a]);

	      // Set style and line color according to ptRatioName
	      gp_ptuncorr[i+a][b]->SetMarkerStyle(markerStyle[a]);
	      gp_ptuncorr[i+a][b]->SetMarkerColor(markerColor[a]);
	      gp_ptuncorr[i+a][b]->SetLineColor(markerColor[a]);
	      gf_ptuncorr[i+a][b]->SetLineColor(markerColor[a]);

	      // Plot gaussplots
	      c2->cd(1+b);
	      gp_ptuncorr[i+a][b]->Draw();
	      gf_ptuncorr[i+a][b]->Draw("same");

	      objToBeWritten.push_back(gp_ptuncorr[i+a][b]);
	      objToBeWritten.push_back(gf_ptuncorr[i+a][b]);
	    } // End of loop over example ptbins
	  c2->Draw();
	  ps->NewPage();
	} // End of loop over ptratios

      c1->cd();
      for(int j = 0 ; j < 8 ; ++j)
	{
	  hists_ptuncorr[i][j]->Draw("P");
	  hists_ptuncorr[i][j]->SetStats(0);
	  hists_ptuncorr[i+1][j]->Draw("P SAME");
	  hists_ptuncorr[i+2][j]->Draw("P SAME");
	  leg->Draw("SAME");
	  c1->SetGrid();
	  c1->Draw();   
	  ps->NewPage(); 

	  objToBeWritten.push_back(hists_ptuncorr[i][j]);
	  objToBeWritten.push_back(hists_ptuncorr[i+1][j]);
	  objToBeWritten.push_back(hists_ptuncorr[i+2][j]);
	}
    }   // End of loop over eta-bins




  // Control quantities from the Pt ratio vs Et gamma.
  // First dimension 12:
  //   index i = 0,1,2 indicates plotted pt-ratio
  //     0: ptjet/etgamma
  //     1: ptjetcorr/etgamma
  //     2: ptjet/ptjetcorr
  // Second dimension 8 is the above specified control
  // quantity
  TH1F* hists_pttrue[3][8];

  // Projected distributions for 4 example eta-bins
  // and the Gauss fits
  TH1F* gp_pttrue[3][4];  
  TF1* gf_pttrue[3][4];

  for(int a = 0; a < 3; a++) // Loop over Pt ratios
    {
      hpt[a]->SetMinimum(0.5);
      hpt[a]->SetMaximum(1.2);
      hpt[a]->SetMarkerStyle(markerStyle[a]);
      hpt[a]->SetMarkerColor(markerColor[a]);
      hpt[a]->SetLineColor(markerColor[a]);
    }

  // Do projections and determine control quantities
  Fit2D(hpt[0],hists_pttrue[0],gp_pttrue[0], gf_pttrue[0]);
  Fit2D(hpt[1],hists_pttrue[1],gp_pttrue[1], gf_pttrue[1]);
  Fit2D(hpt[2],hists_pttrue[2],gp_pttrue[2], gf_pttrue[2]); 
  for(int a = 0; a<3;++a)
    {
      for(int b = 0 ; b < 4 ; ++b)
	{
	  hists_pttrue[a][b]->SetMinimum(0.2);
	  hists_pttrue[a][b]->SetMaximum(1.8);
	  ++b;
	  hists_pttrue[a][b]->SetMinimum(0.0);
	  hists_pttrue[a][b]->SetMaximum(0.5);
	}
    }

  // Draw gaussplots for example pt bins
  // on multi-canvas
  for(int a = 0; a < 3; a++) // Loop over ptratios
    {
      for(int b = 0; b < 3; b++) // Loop over example pt bins
	{
	  // Find pt bin of gaussplot
	  int bin = 0;
	  if(  b == 0  )  bin = int(hpt[a]->GetNbinsX()/6);
	  else if(  b == 1  )  bin = int(hpt[a]->GetNbinsX()/3);
	  else if(  b == 2  )  bin = int(hpt[a]->GetNbinsX()/2);
	  float min = hpt[a]->GetXaxis()->GetBinLowEdge(bin);
	  float max = min + hpt[a]->GetXaxis()->GetBinWidth(bin);

	  // Set title according to pt bin
	  sprintf(title,"#gamma-jet, %.1f < E^{#gamma}_{T} < %.1f GeV",min,max);
	  gp_pttrue[a][b]->SetTitle(title);
	  gp_pttrue[a][b]->SetXTitle(ptRatioName[a]);

	  // Set style and line color according to ptRatioName
	  gp_pttrue[a][b]->SetMarkerStyle(markerStyle[a]);
	  gp_pttrue[a][b]->SetMarkerColor(markerColor[a]);
	  gp_pttrue[a][b]->SetLineColor(markerColor[a]);
	  gf_pttrue[a][b]->SetLineColor(markerColor[a]);

	  // Plot gaussplots
	  c2->cd(1+b);
	  gp_pttrue[a][b]->Draw();
	  gf_pttrue[a][b]->Draw("same");

	  objToBeWritten.push_back(gp_pttrue[a][b]);
	  objToBeWritten.push_back(gf_pttrue[a][b]);
	} // End of loop over example ptbins
      c2->Draw();
      ps->NewPage();
    } // End of loop over ptratios

  c1->cd();
  for(int j = 0 ; j < 8 ; ++j) // Loop over control quantities
    {
      hists_pttrue[0][j]->Draw("P");
      hists_pttrue[0][j]->SetStats(0);
      hists_pttrue[1][j]->Draw("P SAME");
      hists_pttrue[2][j]->Draw("P SAME");
      leg->Draw("SAME");
      c1->SetGrid();
      c1->Draw();   
      ps->NewPage(); 

      objToBeWritten.push_back(hists_pttrue[0][j]);
      objToBeWritten.push_back(hists_pttrue[1][j]);
      objToBeWritten.push_back(hists_pttrue[2][j]);
    }



  // Control quantities from the Pt ratio vs Et gamma
  // in log scale.
  // First dimension 12:
  //   index i = 0,1,2 indicates plotted pt-ratio
  //     0: ptjet/etgamma
  //     1: ptjetcorr/etgamma
  //     2: ptjet/ptjetcorr
  // Second dimension 8 is the above specified control
  // quantity
  TH1F* hists_ptlog[3][8];

  // Projected distributions for 4 example eta-bins
  // and the Gauss fits
  TH1F* gp_ptlog[3][4];  
  TF1* gf_ptlog[3][4];

  for(int a = 0; a < 3; a++) // Loop over Pt ratios
    {
      hptlog[a]->SetMinimum(0.5);
      hptlog[a]->SetMaximum(1.2);
      hptlog[a]->SetMarkerStyle(markerStyle[a]);
      hptlog[a]->SetMarkerColor(markerColor[a]);
      hptlog[a]->SetLineColor(markerColor[a]);
    }

  // Do projections and determine control quantities
  Fit2D(hptlog[0],hists_ptlog[0],gp_ptlog[0], gf_ptlog[0]);
  Fit2D(hptlog[1],hists_ptlog[1],gp_ptlog[1], gf_ptlog[1]);
  Fit2D(hptlog[2],hists_ptlog[2],gp_ptlog[2], gf_ptlog[2]); 
  for(int a = 0; a<3;++a)
    {
      for(int b = 0 ; b < 4 ; ++b)
	{
	  hists_ptlog[a][b]->SetMinimum(0.2);
	  hists_ptlog[a][b]->SetMaximum(1.8);
	  ++b;
	  hists_ptlog[a][b]->SetMinimum(0.0);
	  hists_ptlog[a][b]->SetMaximum(0.5);
	}
    }

  c1->cd();
  for(int j = 0 ; j < 8 ; ++j) // Loop over control quantities
    {
      hists_ptlog[0][j]->Draw("P");
      hists_ptlog[0][j]->SetStats(0);
      hists_ptlog[1][j]->Draw("P SAME");
      hists_ptlog[2][j]->Draw("P SAME");
      leg->Draw("SAME");
      c1->SetGrid();
      c1->SetLogx(1);
      c1->Draw();   
      ps->NewPage(); 

      objToBeWritten.push_back(hists_ptlog[0][j]);
      objToBeWritten.push_back(hists_ptlog[1][j]);
      objToBeWritten.push_back(hists_ptlog[2][j]);
    }
  c1->SetLogx(0);



  // Control quantities from the Pt ratio vs uncorrected
  // jet energy.
  // First dimension 12:
  //   index i = 0,1,2 indicates plotted pt-ratio
  //     0: ptjet/etgamma
  //     1: ptjetcorr/etgamma
  //     2: ptjet/ptjetcorr
  // Second dimension 8 is the above specified control
  // quantity
  TH1F* hists_energy[3][8];

  // Projected distributions for 4 example eta-bins
  // and the Gauss fits
  TH1F* gp_energy[3][4];  
  TF1* gf_energy[3][4];

  for(int a = 0; a < 3; a++) // Loop over Pt ratios
    {
      henergy[a]->SetMinimum(0.5);
      henergy[a]->SetMaximum(1.2);
      henergy[a]->SetMarkerStyle(markerStyle[a]);
      henergy[a]->SetMarkerColor(markerColor[a]);
      henergy[a]->SetLineColor(markerColor[a]);
    }

  // Do projections and determine control quantities
  Fit2D(henergy[0],hists_energy[0],gp_energy[0], gf_energy[0]);
  Fit2D(henergy[1],hists_energy[1],gp_energy[1], gf_energy[1]);
  Fit2D(henergy[2],hists_energy[2],gp_energy[2], gf_energy[2]); 
  for(int a = 0; a<3;++a)
    {
      for(int b = 0 ; b < 4 ; ++b)
	{
	  hists_energy[a][b]->SetMinimum(0.2);
	  hists_energy[a][b]->SetMaximum(1.8);
	  ++b;
	  hists_energy[a][b]->SetMinimum(0.0);
	  hists_energy[a][b]->SetMaximum(0.5);
	}
    }

  // Draw gaussplots for example energy bins
  // on multi-canvas
  for(int a = 0; a < 3; a++) // Loop over ptratios
    {
      for(int b = 0; b < 3; b++) // Loop over example energy bins
	{
	  // Find pt bin of gaussplot
	  int bin = 0;
	  if(  b == 0  )  bin = int(henergy[a]->GetNbinsX()/6);
	  else if(  b == 1  )  bin = int(henergy[a]->GetNbinsX()/3);
	  else if(  b == 2  )  bin = int(henergy[a]->GetNbinsX()/2);
	  float min = henergy[a]->GetXaxis()->GetBinLowEdge(bin);
	  float max = min + henergy[a]->GetXaxis()->GetBinWidth(bin);

	  // Set title according to energy bin
	  sprintf(title,"#gamma-jet, %.1f < E^{jet} < %.1f GeV",min,max);
	  gp_energy[a][b]->SetTitle(title);
	  gp_energy[a][b]->SetXTitle(ptRatioName[a]);

	  // Set style and line color according to ptRatioName
	  gp_energy[a][b]->SetMarkerStyle(markerStyle[a]);
	  gp_energy[a][b]->SetMarkerColor(markerColor[a]);
	  gp_energy[a][b]->SetLineColor(markerColor[a]);
	  gf_energy[a][b]->SetLineColor(markerColor[a]);

	  // Plot gaussplots
	  c2->cd(1+b);
	  gp_energy[a][b]->Draw();
	  gf_energy[a][b]->Draw("same");

	  objToBeWritten.push_back(gp_energy[a][b]);
	  objToBeWritten.push_back(gf_energy[a][b]);
	} // End of loop over example ptbins
      c2->Draw();
      ps->NewPage();
    } // End of loop over ptratios

  c1->cd();
  for(int j = 0 ; j < 8 ; ++j) // Loop over control quantities
    {
      hists_energy[0][j]->Draw("P");
      hists_energy[0][j]->SetStats(0);
      hists_energy[1][j]->Draw("P SAME");
      hists_energy[2][j]->Draw("P SAME");
      leg->Draw("SAME");
      c1->SetGrid();
      c1->Draw();   
      ps->NewPage(); 

      objToBeWritten.push_back(hists_energy[0][j]);
      objToBeWritten.push_back(hists_energy[1][j]);
      objToBeWritten.push_back(hists_energy[2][j]);
    }




  // Control quantities from the Pt ratio vs emf.
  // First dimension 12:
  //   index i = 0,1,2 indicates plotted pt-ratio
  //     0: ptjet/etgamma
  //     1: ptjetcorr/etgamma
  //     2: ptjet/ptjetcorr
  // Second dimension 8 is the above specified control
  // quantity
  TH1F* hists_emf[3][8];

  // Projected distributions for 4 example eta-bins
  // and the Gauss fits
  TH1F* gp_emf[3][4];  
  TF1* gf_emf[3][4];

  for(int a = 0; a < 3; a++) // Loop over Pt ratios
    {
      hemf[a]->SetMinimum(0.5);
      hemf[a]->SetMaximum(1.2);
      hemf[a]->SetMarkerStyle(markerStyle[a]);
      hemf[a]->SetMarkerColor(markerColor[a]);
      hemf[a]->SetLineColor(markerColor[a]);
    }

  // Do projections and determine control quantities
  Fit2D(hemf[0],hists_emf[0],gp_emf[0], gf_emf[0]);
  Fit2D(hemf[1],hists_emf[1],gp_emf[1], gf_emf[1]);
  Fit2D(hemf[2],hists_emf[2],gp_emf[2], gf_emf[2]); 
  for(int a = 0; a<3;++a)
    {
      for(int b = 0 ; b < 4 ; ++b)
	{
	  hists_emf[a][b]->SetMinimum(0.2);
	  hists_emf[a][b]->SetMaximum(1.8);
	  ++b;
	  hists_emf[a][b]->SetMinimum(0.0);
	  hists_emf[a][b]->SetMaximum(0.5);
	}
    }

  // Draw gaussplots for example emf bins
  // on multi-canvas
  for(int a = 0; a < 3; a++) // Loop over ptratios
    {
      for(int b = 0; b < 3; b++) // Loop over example emf bins
	{
	  // Find emf bin of gaussplot
	  int bin = 0;
	  if(  b == 0  )  bin = int(hemf[a]->GetNbinsX()/6);
	  else if(  b == 1  )  bin = int(hemf[a]->GetNbinsX()/3);
	  else if(  b == 2  )  bin = int(hemf[a]->GetNbinsX()/2);
	  float min = hemf[a]->GetXaxis()->GetBinLowEdge(bin);
	  float max = min + hemf[a]->GetXaxis()->GetBinWidth(bin);

	  // Set title according to emf bin
	  sprintf(title,"#gamma-jet, %.2f < f_{em} < %.2f",min,max);
	  gp_emf[a][b]->SetTitle(title);
	  gp_emf[a][b]->SetXTitle(ptRatioName[a]);

	  // Set style and line color according to ptRatioName
	  gp_emf[a][b]->SetMarkerStyle(markerStyle[a]);
	  gp_emf[a][b]->SetMarkerColor(markerColor[a]);
	  gp_emf[a][b]->SetLineColor(markerColor[a]);
	  gf_emf[a][b]->SetLineColor(markerColor[a]);

	  // Plot gaussplots
	  c2->cd(1+b);
	  gp_emf[a][b]->Draw();
	  gf_emf[a][b]->Draw("same");

	  objToBeWritten.push_back(gp_emf[a][b]);
	  objToBeWritten.push_back(gf_emf[a][b]);
	} // End of loop over example emf bins
      c2->Draw();
      ps->NewPage();
    } // End of loop over ptratios

  c1->cd();
  for(int j = 0 ; j < 8 ; ++j) // Loop over control quantities
    {
      hists_emf[0][j]->Draw("P");
      hists_emf[0][j]->SetStats(0);
      hists_emf[1][j]->Draw("P SAME");
      hists_emf[2][j]->Draw("P SAME");
      leg->Draw("SAME");
      c1->SetGrid();
      c1->Draw();   
      ps->NewPage(); 

      objToBeWritten.push_back(hists_emf[0][j]);
      objToBeWritten.push_back(hists_emf[1][j]);
      objToBeWritten.push_back(hists_emf[2][j]);
    }
  ps->NewPage(); 





  // Test TCaliber::FlattenSpectra()
  // Plot number of events per pt, eta, and emf bin
  for(int i = 0 ; i < 3 ; ++i)
    {
      hptGamma[i]->SetMarkerStyle(20);
      hptGamma[i]->SetMarkerColor(1);
      hptGammaW[i]->SetMarkerStyle(22);
      hptGammaW[i]->SetMarkerColor(2);
      if( hptGamma[i]->GetMaximum()>hptGammaW[i]->GetMaximum() )
	hptGammaW[i]->SetMaximum( hptGamma[i]->GetMaximum() );

      hptGammaW[i]->Draw("p");
      hptGammaW[i]->SetStats(0);
      hptGamma[i]->Draw("pSAME");
      c1->SetLogx(0);  
      c1->SetLogy(1);  
      c1->SetGrid();
      leg->Clear();
      leg->AddEntry(hptGamma[i],"no weights","p");
      leg->AddEntry(hptGammaW[i],"with weights","p");
      leg->Draw("SAME");
      c1->Draw(); 
    }


  // emf vs pt before and after TCaliber::FlattenSpectra()
  c1->SetLogx(0);  
  c1->SetLogy(0);   
  c1->SetGrid(0);
  hptGamma2D->SetMarkerStyle(7);
  hptGamma2D->Draw("hist");
  c1->Draw(); 
  ps->NewPage(); 

  hptGamma2DW->SetMarkerStyle(7);
  hptGamma2DW->Draw("hist");
  c1->Draw(); 
  
  c1->SetLogx(0);  
  c1->SetLogy(0);   
  c1->SetGrid(0);


  // tower response
  // 0: measured HCAL E = 0 GeV
  // 1: measured HCAL E = 5 GeV
  // 2: measured HCAL E = 50 GeV
  TH1F* htow[3];
  htow[0] = new TH1F("htow0","#gamma-jet,  tower response;id_{#eta}",82,0,82);
  htow[1] = (TH1F*)htow[0]->Clone("htow1");
  htow[2] = (TH1F*)htow[0]->Clone("htow2");
  
  for (int eta=-41; eta < 41;++eta)
    {
      int i = _par->GetEtaBin(eta >= 0 ? eta +1 : eta);
      TMeasurement x;
      double* par = _par->GetTowerParRef(_par->GetBin(i,0));
      htow[0]->Fill(eta + 41,TParameters::tower_parametrization(&x,par));
      x.HadF = 5;
      htow[1]->Fill(eta + 41,TParameters::tower_parametrization(&x,par));
      x.HadF = 50;
      htow[2]->Fill(eta + 41,TParameters::tower_parametrization(&x,par));
    }

  for(int i = 0; i < 3; i++)
    {
      htow[i]->SetMinimum(-5);
      htow[i]->SetMaximum(80);
      htow[i]->SetStats(0);  
      htow[i]->SetMarkerStyle(markerStyle[i]);
      htow[i]->SetMarkerColor(markerColor[i]);
      htow[i]->SetLineColor(markerColor[i]);
      if( i == 0 ) htow[i]->Draw("p");
      else htow[i]->Draw("p SAME");
      objToBeWritten.push_back(htow[i]);
    }
  leg->Clear();
  leg->SetHeader("E_{had} [GeV]");
  leg->AddEntry(htow[0]," 0","p");
  leg->AddEntry(htow[1]," 5","p");
  leg->AddEntry(htow[2],"50","p");
  leg->Draw("SAME");
  c1->Draw(); 
  ps->Close();


  // Closing ps file and writing objects to .root file
  ps->Close();
  if( _outputROOT ) WriteToRootFile(objToBeWritten, "GammaJet");


  // free memory
  for(int i = 0; i < 4; i++)
    {
      delete towerinjet[i];
      delete respvstet[i];
    }
  for(int i = 0; i < 21; i++) delete leadToNext[i];
  for(int i = 0; i < 12; i++)
    {
      delete heta[i];
      delete hpt_uncorr[i];
    }
  delete EtaPhiMap;
  for(int i = 0; i < 3; i++)
    {
      delete hpt[i];
      delete henergy[i];
      delete hemf[i];
      delete hptlog[i];
      delete hptGamma[i];
      delete hptGammaW[i];
    }
  delete hptGamma2D;
  delete hptGamma2DW;

  for(int i = 0; i < 12; i++)
    {
      for(int j = 0; j < 4; j++)
	{
	  delete hists_eta[i][j];
	  delete hists_eta[i][j+4];
	  delete gp_eta[i][j];
	  delete gf_eta[i][j];

	  delete hists_ptuncorr[i][j];
	  delete hists_ptuncorr[i][j+4];
	  delete gp_ptuncorr[i][j];
	  delete gf_ptuncorr[i][j];
	}
    }          
  for(int i = 0; i < 3; i++)
    {
      for(int j = 0; j < 4; j++)
	{
	  delete hists_pttrue[i][j];
	  delete hists_pttrue[i][j+4];
	  delete gp_pttrue[i][j];
	  delete gf_pttrue[i][j];

	  delete hists_ptlog[i][j];
	  delete hists_ptlog[i][j+4];
	  delete gp_ptlog[i][j];
	  delete gf_ptlog[i][j];

	  delete hists_energy[i][j];
	  delete hists_energy[i][j+4];
	  delete gp_energy[i][j];
	  delete gf_energy[i][j];

	  delete hists_emf[i][j];
	  delete hists_emf[i][j+4];
	  delete gp_emf[i][j];
	  delete gf_emf[i][j];
	}
    }
  delete htow[0];
  delete htow[1];
  delete htow[2];
          
  delete leg;

  delete ps;
  delete c1;
  delete c2;
}



//---------------------------------------------------------------
//   Gamma-Jet Control Histograms per tower bin
//   orig name: gammajet_plots_per_towerbin.ps
//---------------------------------------------------------------
void TControlPlots::MakeControlPlotsGammaJetPerTowerBin()
{
  std::vector<TObject*> objToBeWritten;
  std::vector<TObject*> objToBeDeleted;
  std::vector<TData*>::const_iterator data_it, it;

  TCanvas * const c1 = new TCanvas("c1","",600,600);
  objToBeDeleted.push_back(c1);
  TPostScript * const ps = new TPostScript("controlplotsGammaJetPerTowerBin.ps",111);
  objToBeDeleted.push_back(ps);

  //one plot per *TOWER* bin!
  for (int eta=0; eta<_par->GetEtaGranularity();++eta) // Loop over eta bins
    {
      for (int phi=0; phi<_par->GetPhiGranularity();++phi) // Loop over phi bins
	{
	  int i = _par->GetBin(eta,phi);

	  // Initialize histos
	  char * name = new char[100];
	  sprintf(name, "hjes_gj%d_eta%d_phi%d",i,eta+1,phi+1);
	  TH1F * plot_jes = new TH1F(name,";#sum calibrated tower P_{T} [GeV]; JES: ( P_{T}^{#gamma} / #sum P_{T}^{calib. tower})",100,0.0,400.);    
	  objToBeDeleted.push_back(plot_jes);

	  sprintf(name, "h_gj%d_eta%d_phi%d",i,eta+1,phi+1);
	  TH1F * plot = new TH1F(name,";uncorrected jet P_{T} [GeV];average of ( P_{T}^{#gamma} / P_{T}^{uncalib. jet})",100,0.0,400.);    
	  objToBeDeleted.push_back(plot);

	  sprintf(name, "gj_fit%d",i);
	  TH1F * fit  = new TH1F(name,"",100,0.0,400.);    
	  objToBeDeleted.push_back(fit);

	  sprintf(name, "gj_norm%d",i);
	  TH1F * norm = new TH1F(name,"",100,0.0,400.);    
	  objToBeDeleted.push_back(norm);

	  sprintf(name, "gjjes_norm%d",i);
	  TH1F * norm_jes = new TH1F(name,"",100,0.0,400.);    
	  objToBeDeleted.push_back(norm_jes);


	  int indexJet=0, ijets=0;      
	  data_it = _data->begin();
	  for (; data_it != _data->end();++data_it) // loop over all fit-events
	    {
	      if ( (*data_it)->GetType()!=GammaJet ) continue;

	      int indexTower = 0; // Index of max tower
	      double Etmax = 0.;
	      double calib_tower_sum = 0.;
	      double tower_sum = 0.; //is equivalent to (*data_it)->GetMess(),
	                             //but since we need the index too, this is faster
	      const std::vector<TData*>& data_ref = (*data_it)->GetRef();
	      for (it = data_ref.begin(); it != data_ref.end(); ++it) // Loop over towers
		{
		  double tow_et = (*it)->GetMess()->pt;
		  if (tow_et>Etmax)
		    {
		      indexTower = (*it)->GetIndex();
		      Etmax = tow_et;
		    }
		  tower_sum += tow_et; //*tow_et[7];
		  calib_tower_sum += (*it)->GetParametrizedMess(); //*tow_et[7];
		} // End of loop over towers

	      if ( indexTower!=i ) continue; //event belongs to a wrong bin

	      indexJet += (*data_it)->GetIndex();
	      ++ijets;

	      double JetCorr = (*data_it)->GetParametrizedMess();

	      fit->Fill(      tower_sum, JetCorr/tower_sum );
	      plot->Fill(     tower_sum, (*data_it)->GetTruth()/tower_sum );
	      norm->Fill(     tower_sum ); 
	      plot_jes->Fill( calib_tower_sum, (*data_it)->GetTruth()/calib_tower_sum );
	      norm_jes->Fill( calib_tower_sum ); 
	    } // End of loop over all fit-events

	  if ( norm->GetEntries()==0 ) continue;
      
	  // Normalize histos
	  plot->Divide(norm);
	  plot_jes->Divide(norm_jes);
	  fit->Divide(norm);
	  for (int b = 1; b <= norm->GetNbinsX(); ++b)
	    {
	      if (norm->GetBinContent(b)>0)
		{
		  plot->SetBinError(  b, 1./sqrt(
						 norm->GetBinContent(b)     //stat
						 )*plot->GetBinContent(b) );
		}			      
	      if (norm_jes->GetBinContent(b)>0)
		{
		  plot_jes->SetBinError(  b, 1./sqrt(
						     norm_jes->GetBinContent(b)     //stat
						     )*plot_jes->GetBinContent(b) );
		}			      
	    }
	
	  fit->SetLineColor( 2 );
	  fit->SetLineWidth( 4 );
	  plot->SetMarkerStyle( 8 );
	  plot->SetMinimum(0.1);
      
	  //c1->SetLogy(1);
	  plot->Draw("pe");
	  objToBeWritten.push_back(plot);
	  fit->Draw("h,same");
	  objToBeWritten.push_back(fit);

	  TLatex latex;
	  latex.SetTextSize(0.035);
	  latex.DrawLatex( 0.3*(plot->GetXaxis()->GetXmax()-plot->GetXaxis()->GetXmin()),
		       0.6*(plot->GetMaximum()-plot->GetMinimum()),
		       "#color[2]{--  tower and jet corrections}");

	  c1->Draw(); 
	  ps->NewPage();

	  plot_jes->Draw("pe");
	  objToBeWritten.push_back(plot_jes);

	  sprintf(name,"res2_%i",i);
	  //TF1 * res2 = new TF1(name,_par->jes_plot_parametrization, 0.5, 400., 3);
	  //objToBeDeleted.push_back(res2);
	  //i = indexJet/ijets - _par->GetNumberOfTowerParameters();
	  //double * val = _par->GetJetParRef(i);
	  //res2->SetParameters(val[0],val[1]);
	  //res2->SetLineWidth( 3 );
	  //res2->SetLineColor( 2 );
	  //res2->Draw("same");
	  //objToBeWritten.push_back(res2);
	  c1->Draw(); 
	  ps->NewPage();
    }
  }

  ps->Close();
  if( _outputROOT ) WriteToRootFile( objToBeWritten, "GammaJetPerTowerBin" );
  objToBeDeleted.clear();
}




//---------------------------------------------------------------
//   Gamma-Jet Control Histograms per jet bin
//   orig name: gammajet_plots_per_jetbin.ps
//---------------------------------------------------------------
void TControlPlots::MakeControlPlotsGammaJetPerJetBin()
{
  std::vector<TObject*> objToBeDeleted;
  std::vector<TObject*> objToBeWritten;
  std::vector<TData*>::const_iterator data_it, it;

  TCanvas * const c1 = new TCanvas("c1","",600,600);
  objToBeDeleted.push_back(c1);
  TPostScript * const ps = new TPostScript("controlplotsGammaJetPerJetBin.ps",111);
  objToBeDeleted.push_back(ps);

  //one plot per *JET* bin!
  for (int eta=0; eta<_par->GetEtaGranularityJet();++eta)
    {
      for (int phi=0; phi<_par->GetPhiGranularityJet();++phi)
	{
	  int i = _par->GetJetBin(eta,phi)*_par->GetNumberOfJetParametersPerBin() + _par->GetNumberOfTowerParameters();
	  char * name = new char[100];
	  sprintf(name, "h2jes_gj%d_eta%d_phi%d",i,eta+1,phi+1);
	  TH1F * plot_jes = new TH1F(name,";#sum calibrated tower E_{T} [GeV]; JES: ( E_{T}^{#gamma} / #sum E_{T}^{calib. tower})",100,0.0,400.);
	  objToBeDeleted.push_back(plot_jes);
	  sprintf(name, "h2_gj%d_eta%d_phi%d",i,eta+1,phi+1);
	  TH1F * plot = new TH1F(name,";uncorrected jet E_{T} [GeV];average of ( E_{T}^{#gamma} / E_{T}^{uncalib. jet})",100,0.0,400.);    
	  objToBeDeleted.push_back(plot);
	  sprintf(name, "gj2_fit%d",i);
	  TH1F * fit  = new TH1F(name,"",100,0.0,400.);    
	  objToBeDeleted.push_back(fit);
	  sprintf(name, "gj2_norm%d",i);
	  TH1F * norm = new TH1F(name,"",100,0.0,400.);    
	  objToBeDeleted.push_back(norm);
	  sprintf(name, "gj2jes_norm%d",i);
	  TH1F * norm_jes = new TH1F(name,"",100,0.0,400.);    
	  objToBeDeleted.push_back(norm_jes);

	  //loop over all fit-events
	  data_it = _data->begin();
	  for (; data_it != _data->end();++data_it)
	    {
	      if ( (*data_it)->GetType()!=GammaJet ) continue;
	      if ( (*data_it)->GetIndex()!= i )	  continue; //event belongs to a wrong bin
	
	      double JetCorr = (*data_it)->GetParametrizedMess();

	      double tower_sum = 0.0;
	      double calib_tower_sum = 0.0;

	      const std::vector<TData*>& data_ref = (*data_it)->GetRef(); // Tower
	      for (it = data_ref.begin(); it!=data_ref.end(); ++it)
		{
		  tower_sum += (*it)->GetMess()->pt; // * (*it)->GetMess()[7];
		  calib_tower_sum += (*it)->GetParametrizedMess(); // * (*it)->GetMess()[7];
		}

	      fit->Fill( tower_sum, JetCorr/tower_sum );
	      plot->Fill( tower_sum, (*data_it)->GetTruth()/tower_sum );
	      norm->Fill( tower_sum ); 

	      plot_jes->Fill( calib_tower_sum, (*data_it)->GetTruth()/calib_tower_sum );
	      norm_jes->Fill( calib_tower_sum ); 
	    }

	  if (norm->GetEntries()==0) continue;
	  plot->Divide(norm);
	  fit->Divide(norm);
	  plot_jes->Divide(norm_jes);

	  for(int b = 1; b <= norm->GetNbinsX(); b++)
	    {
	      if (norm->GetBinContent(b)>0)
		{
		  plot->SetBinError(  b, 1./sqrt(
						 norm->GetBinContent(b)     //stat
						 )*plot->GetBinContent(b) );
		}			      
	      if (norm_jes->GetBinContent(b)>0)
		{
		  plot_jes->SetBinError(  b, 1./sqrt(
						     norm_jes->GetBinContent(b)     //stat
						     )*plot_jes->GetBinContent(b) );
		}			      
	    }
	
	  fit->SetLineColor( 2 );
	  fit->SetLineWidth( 4 );
	  plot->SetMarkerStyle( 8 );
	  plot->SetMinimum(0.1);
      
	  plot->Draw("pe");
	  objToBeWritten.push_back( plot );
	  fit->Draw("h,same");
	  objToBeWritten.push_back( fit );

	  TLatex latex;
	  latex.SetTextSize(0.035);
	  latex.DrawLatex( 0.3*(plot->GetXaxis()->GetXmax()-plot->GetXaxis()->GetXmin()),
			   0.6*(plot->GetMaximum()-plot->GetMinimum()),
			   "#color[2]{---  tower and jet corrections}");

	  c1->Draw(); 
	  ps->NewPage();


	  plot_jes->Draw("pe");
	  objToBeWritten.push_back( plot_jes );

	  sprintf(name,"res2_%i",i);
	  //TF1 * res2 = new TF1(name,_par->jes_plot_parametrization, 0.5, 400., 3);
	  //objToBeDeleted.push_back(res2);
	  //i = _par->GetJetBin(eta, phi);
	  //double * val = _par->GetJetParRef(i);
	  //res2->SetParameters(val[0],val[1],val[2]);
	  //res2->SetLineWidth( 3 );
	  //res2->SetLineColor( 2 );
	  //res2->Draw("same");
	  //objToBeWritten.push_back( res2 );

	  c1->Draw(); 
	  ps->NewPage();
	}
    }
  ps->Close();
  if( _outputROOT ) WriteToRootFile( objToBeWritten, "GammaJetPerJetBin" );
  objToBeDeleted.clear();
}



//---------------------------------------------------------------
//   Gamma-Jet Control Histograms
//   orig name: sigmas_gammajet.ps
//---------------------------------------------------------------
void TControlPlots::MakeControlPlotsGammaJetSigmas()
{
  std::vector<TObject*> objToBeWritten;

  TCanvas * const c1 = new TCanvas("c1","",600,600);
  TPostScript * const ps = new TPostScript("controlplotsGammaJetSigmas.ps",111);

  int nPtBins = 200;

  TH1F * gauss_forpt[nPtBins];
  TH1F * gauss_forptcorr[nPtBins];
  gauss_forpt[0] = new TH1F("hgauss0","p^{jet}_{T} bin [0..1GeV];(p^{jet}_{T} - p^{#gamma}_{T}) / p^{jet}_{T}",600,-3,3);
  gauss_forptcorr[0] = new TH1F("hgausscorr0","p^{corr. jet}_{T} bin [0..1GeV];(p^{corr. jet}_{T} - p^{#gamma}_{T}) / p^{corr. jet}_{T}",600,-3,3);

  char name[100];
  for(int i = 1 ; i < nPtBins ; ++i)
    {
      sprintf(name,"hgauss%i",i);
      gauss_forpt[i] = (TH1F*)gauss_forpt[0]->Clone(name);
      sprintf(name,"p^{jet}_{T} bin [%d..%dGeV]",i,i+1);
      gauss_forpt[i]->SetTitle(name);

      sprintf(name,"hgausscorr%i",i);
      gauss_forptcorr[i] = (TH1F*)gauss_forptcorr[0]->Clone(name);
      sprintf(name,"p^{corr. jet}_{T} bin [%d..%dGeV]",i,i+1);
      gauss_forptcorr[i]->SetTitle(name);
  }
  
  //loop over all fit-events
  for ( std::vector<TData*>::const_iterator i = _data->begin(); i != _data->end() ; ++i )
    {
      if( (*i)->GetType() != GammaJet ) continue;

      double etjetcor = (*i)->GetParametrizedMess();
      if((*i)->GetMess()->pt>0 && (*i)->GetMess()->pt<200)
	gauss_forpt[(int)(*i)->GetMess()->pt]->Fill( ((*i)->GetMess()->pt - (*i)->GetTruth())/(*i)->GetMess()->pt, (*i)->GetWeight() );
      if(etjetcor>0 && etjetcor<200)
	gauss_forptcorr[(int)etjetcor]->Fill( (etjetcor - (*i)->GetTruth())/etjetcor, (*i)->GetWeight() );
    }


  TText * text = new TText();
  text->SetTextFont(42);
  text->SetTextSize(0.03);
  text->SetTextColor(2);
  
  TF1 * f[nPtBins];

  //TF1 * f = new TF1("gauss_step",gauss_step,-10,10,5);
  //double edge;
  for(int i = 0 ; i < nPtBins ; ++i) // Loop over pt bins
    {
      //edge = 1.0-20./(((double)i)+0.5);
      //f->SetParameters(-1.,2.0,3.0, edge, 0.01);
      //f->FixParameter(3, edge);
      //f->FixParameter(4, 0.01);
    
      gauss_forpt[i]->Fit("gaus","LLQNO","");
      sprintf(name,"fit_gausscorr%i",i);
      f[i] = (TF1*)gROOT->GetFunction("gaus")->Clone(name);
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
      objToBeWritten.push_back(gauss_forpt[i]);
      f[i]->SetLineColor(2);
      f[i]->Draw("same");
      objToBeWritten.push_back(f[i]);
      sprintf(name,"mean %f",f[i]->GetParameter(1));
      text->DrawText(1,0.7*gauss_forpt[i]->GetMaximum(),name);
      //func->Draw("same");
      c1->Draw();
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
    } // End of loop over pt bins
  
  ps->Close();

  if( _outputROOT ) WriteToRootFile( objToBeWritten, "GammaJetSigmas" );

  for(int i = 0 ; i < nPtBins ; ++i)
    {
      delete gauss_forpt[i];  
      delete gauss_forptcorr[i];  
      //delete f[i];
    }  
  delete text;
}





void TControlPlots::MakeControlPlotsDiJet()
{
  std::vector<TObject*> objToBeWritten;

  TCanvas * const c1 = new TCanvas("1","",600,600);
  TCanvas * const c2 = new TCanvas("c2","",600,600);
  c2->Divide(2,2);
  TPostScript * const ps = new TPostScript("controlplotsDiJet.ps",111);


  //book hists
  char name[100];

  TH2F* Scale[2];
  Scale[0] = new TH2F("hScaleDiff","Scale;Scale;Scale After Fit - Scale",100,0,1000,100,-100,100);
  Scale[1] = new TH2F("hScaleAfter","Scale;Scale; Scale After Fit",100,0,1000,100,0,1000);

  TH1F* ptspec[4];
  ptspec[0] = new TH1F("hPtSpecProbeJet","p_{T} spectrum probe jet;p_{T} [GeV]",300,0,2000);
  ptspec[1] = new TH1F("hPtSpecBarrelJet","p_{T} spectrum barrel jet;p_{T} [GeV]",300,0,2000);
  ptspec[2] = new TH1F("PtSpecProbeJetBeforeFit","p_{T} spectrum probe jet before fit;p_{T} [GeV]",300,0,2000);
  ptspec[3] = new TH1F("PtSpecBarrelJetBeforeFit","p_{T} spectrum barrel jet before fit;p_{T} [GeV]",300,0,2000);

  TH1F* eta[2];
  eta[0] = new TH1F("hEtaProbeJet","#eta probe jet;#eta",100,-5,5);
  eta[1] = new TH1F("hEtaBarrelJet","#eta barrel jet;#eta",100,-5,5);

  TH1F* dphi[11];
  dphi[0] = new TH1F("hDeltaPhi0","#Delta Phi;#Delta #Phi",70,2.8,3.5);
  for(int i = 1; i < 5; ++i)
    {
      sprintf(name,"hDeltaPhi%i",i);
      dphi[i] = (TH1F*)dphi[0]->Clone(name);
    }
  dphi[1]->SetTitle("#Delta #Phi (P^{scale}_{T} 10-35 GeV);#Delta #Phi");
  dphi[2]->SetTitle("#Delta #Phi (P^{scale}_{T} 35-90 GeV);#Delta #Phi");
  dphi[3]->SetTitle("#Delta #Phi (P^{scale}_{T} 90-300 GeV);#Delta #Phi");
  dphi[4]->SetTitle("#Delta #Phi (P^{scale}_{T} 300+ GeV);#Delta #Phi");
  dphi[5] = new TH1F("hDeltaPhiOff","#Delta Phi;#Delta #Phi",120,-3.4,3.4);
  dphi[6] = new TH1F("hDeltaPhiOff+","#Delta Phi;#Delta #Phi",120,2.8,3.4);
  dphi[7] = new TH1F("hAbsDeltaPhiOff","#Delta Phi;#Delta #Phi",120,2.8,3.4);
  dphi[8] = new TH1F("hDeltaPhiwoAbs","#Delta Phi;#Delta #Phi",140,-3.5,3.5);
  dphi[9] = new TH1F("hDeltaPhiwoAbs+","#Delta Phi;#Delta #Phi",120,2.8,3.5);
  dphi[10] = new TH1F("hDeltaPhiwoAbs-","#Delta Phi;#Delta #Phi",120,-3.5,-2.8);

  TH2F* Bvsdphi = new TH2F("hBvsDeltaPhi","B vs #Delta #Phi;#Delta#Phi;B",70,2.8,3.5,100,-1,1);

  TH2F* Difvscomb = new TH2F("hDeltaEtvsCombinedJet","#Delta Et vs. Combined Jet;Combined Jet;#Delta E_{T}",100,0,100,100,0,100);

  TH2F* combmean[5];
  combmean[0] = new TH2F("hDiJet0","di-jet controlplot;scale (Pt);Pt of combined jet",100,0,2000,100,0,200);

  TH2F* difmean[5];
  difmean[0] = new TH2F("hDiJetDif0","di-jet controlplot;scale (Pt);abs. difference in Pt",100,0,2000,100,0,200);
  for(int i=1;i<5;++i)
    {
      sprintf(name,"hDiJet%i",i);
      combmean[i] =  (TH2F*)combmean[0]->Clone(name);

      sprintf(name,"hDiJetDif%i",i);
      difmean[i] =  (TH2F*)difmean[0]->Clone(name);
    }
  combmean[1]->SetTitle("abs(#eta) < 1");
  combmean[2]->SetTitle("1 < abs(#eta) < 2");
  combmean[3]->SetTitle("2 < abs(#eta) < 3");
  combmean[4]->SetTitle("3 < abs(#eta) < 4");

  difmean[1]->SetTitle("abs(#eta) < 1");
  difmean[2]->SetTitle("1 < abs(#eta) < 2");
  difmean[3]->SetTitle("2 < abs(#eta) < 3");
  difmean[4]->SetTitle("3 < abs(#eta) < 4");


  TH2F* Beta[8];
  Beta[0] = new TH2F("hBeta0","di-jet;#eta",100,-5,5,100,-0.7,0.7);
  for(int i = 1 ; i < 8 ; ++i)
    {
      sprintf(name,"hBeta%i",i);
      Beta[i] = (TH2F*)Beta[0]->Clone(name);
    }
  Beta[2]->SetTitle("di-jet 10 < E_{T}^{scale} < 35 GeV;#eta");
  Beta[3]->SetTitle("di-jet 10 < E_{T}^{scale} < 35 GeV;#eta");
  Beta[4]->SetTitle("di-jet 35 < E_{T}^{scale} < 90 GeV;#eta");
  Beta[5]->SetTitle("di-jet 35 < E_{T}^{scale} < 90 GeV;#eta");
  Beta[6]->SetTitle("di-jet 90 < E_{T}^{scale} < 300 GeV;#eta");
  Beta[7]->SetTitle("di-jet 90 < E_{T}^{scale} < 300 GeV;#eta");

  TH2F* Bpt[2];
  Bpt[0] = new TH2F("hBpt0","di-jet;p_{T} [GeV]",400,0,400,100,-0.7,0.7);
  Bpt[1] = (TH2F*)Bpt[0]->Clone("hBpt1");

  TH2F* Benergy[2];
  Benergy[0] = new TH2F("hBenergy0","di-jet;E [GeV]",400,0,400,100,-0.7,0.7);
  Benergy[1] = (TH2F*)Benergy[0]->Clone("hBenergy1");
  
  TH2F* Bemf[2];
  Bemf[0] = new TH2F("hBemf0","di-jet;f_{em} (probe jet)",100,0,1,100,-0.7,0.7);
  Bemf[1] = (TH2F*)Bemf[0]->Clone("hBemf1");

  double bins[101];
  for(int i = 0; i < 101 ; ++i)
    {
      bins[i] = pow(10,(i+32)/40.0);
    }
  TH2F* Bptlog[2];
  Bptlog[0] = new TH2F("hBptlog0","di-jet;p_{T} [GeV]",100,bins,100,-0.7,0.7); 
  Bptlog[1] = (TH2F*)Bptlog[0]->Clone("hBptlog1");

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
  for( std::vector<TData*>::const_iterator i = _data->begin() ; i != _data->end() ; ++i )  
    {
      TData* jj = *i;
      if(jj->GetType() != PtBalance) continue;

      TData_MessMess* jm = (TData_MessMess*) jj;
      double etscale = jm->GetScale();


      double etparascale = 0.;
      for(std::vector<TData*>::const_iterator t = jm->GetRef().begin(); t != jm->GetRef().end(); ++t)
	{
	  etparascale += (*t)->GetParametrizedMess();
	}
      etparascale = ( etparascale + jm->GetParametrizedMess() )/2.;
      double etajet1 = jm->GetMess()->eta;
      double etajet2 = (*jm->GetSecondaryJets())[0]->GetMess()->eta;
      double etjetcomb = jm->GetMessCombination();
      double etjet1 = jm->GetParametrizedMess();      //Probe
      double etjet2 = (*jm->GetSecondaryJets())[0]->GetParametrizedMess();      //Barrel
      double etjet1uncor = jm->GetMess()->pt;      //Probe
      double etjet2uncor = (*jm->GetSecondaryJets())[0]->GetMess()->pt;      //Barrel
      double phijet1 = jm->GetMess()->phi;      //Probe
      double phijet2 = (*jm->GetSecondaryJets())[0]->GetMess()->phi;      //Barrel
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
	  double temp = etjet2uncor;
	  etjet2uncor = etjet1uncor;
	  etjet1uncor = temp;
	}
      double deltaphi = fabs(phiprobe - phijet2);
      double deltaphioff = deltaPhi(phiprobe,phijet2);

      Scale[0]->Fill(etscale,etparascale - etscale);
      Scale[1]->Fill(etscale,etparascale);
      ptspec[0]->Fill(etprobe);
      ptspec[1]->Fill(etjet2);
      ptspec[2]->Fill(etjet1uncor);
      ptspec[3]->Fill(etjet2uncor);
      eta[0]->Fill(etaprobe);
      eta[1]->Fill(etajet2);
      dphi[0]->Fill(deltaphi);
      dphi[5]->Fill(deltaphioff);
      dphi[6]->Fill(deltaphioff);
      dphi[7]->Fill(fabs(deltaphioff));
      dphi[8]->Fill(phiprobe - phijet2);             //
      if((phiprobe - phijet2) > 0) dphi[9]->Fill(phiprobe - phijet2);             //
      else dphi[10]->Fill(phiprobe - phijet2);             //
      Bvsdphi->Fill(deltaphi,B); 
      Difvscomb->Fill(fabs(etprobe - etjet2),etjetcomb);
      combmean[0]->Fill(etscale, etjetcomb);
      difmean[0]->Fill(etscale, fabs(etprobe - etjet2));

      Beta[0]->Fill(etaprobe, B,jj->GetWeight());
      Beta[1]->Fill(etaprobe, Buncor,jj->GetWeight());
      if (etscale > 10 && etscale < 35)
	{
	  Beta[2]->Fill(etaprobe, B,jj->GetWeight());
	  Beta[3]->Fill(etaprobe, Buncor,jj->GetWeight());
	  dphi[1]->Fill(deltaphi);
	}
      else if (etscale > 35 && etscale < 90)
	{
	  Beta[4]->Fill(etaprobe, B,jj->GetWeight());
	  Beta[5]->Fill(etaprobe, Buncor,jj->GetWeight());
	  dphi[2]->Fill(deltaphi);
	}
      else if (etscale > 90 && etscale < 300)
	{
	  Beta[6]->Fill(etaprobe, B,jj->GetWeight());
	  Beta[7]->Fill(etaprobe, Buncor,jj->GetWeight());
	  dphi[3]->Fill(deltaphi);
	}
      else if (etscale > 300) dphi[4]->Fill(deltaphi);

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
      for(std::vector<TData*>::const_iterator t = jj->GetRef().begin(); t != jj->GetRef().end(); ++t)
	{
	  TData* tt = *t;
	  em  += tt->GetMess()->EMF;
	  had += tt->GetMess()->HadF;
	  had += tt->GetMess()->OutF;
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
    }  //End of loop over all fit-events



  c1->cd();
  Scale[0]->Draw("box");
  objToBeWritten.push_back(Scale[0]);
  c1->Draw();
  ps->NewPage(); 

  Scale[1]->Draw("box");
  objToBeWritten.push_back(Scale[1]);
  c1->Draw();
  ps->NewPage(); 

  for(int i=0;i<4;++i)
    {
      ptspec[i]->Draw();
      objToBeWritten.push_back(ptspec[i]);
      c1->Draw();
      ps->NewPage(); 
    }
  for(int i=0;i<2;++i)
    {
      eta[i]->Draw();
      objToBeWritten.push_back(eta[i]);
      c1->Draw();
      ps->NewPage(); 
    }
  for(int i=0;i<5;++i)
    {
      dphi[i]->Draw();
      objToBeWritten.push_back(dphi[i]);
      c1->Draw();
      ps->NewPage(); 
    }

  Bvsdphi->Draw();
  objToBeWritten.push_back(Bvsdphi);
  c1->Draw();
  ps->NewPage();  

  Difvscomb->Draw();
  objToBeWritten.push_back(Difvscomb);
  c1->Draw();
  ps->NewPage();  


  // Control quantities vs eta
  TH1F* hists_beta[8][8];
  TH1F* gp_beta[8][4];  
  TF1* gf_beta[8][4];     

  TLegend* leg = new TLegend(0.7,0.7,0.96,0.9);
  leg->SetFillColor(0);
  leg->AddEntry(Beta[1],"B before fit");
  leg->AddEntry(Beta[0],"B after fit");

  int markerColor[2] = { 2,1 };
  int markerStyle[2] = { 22,20 };
  int etLimit[4] = { 10,35,90,300 };

  for(int i = 0 ; i < 8 ; i+=2) // Loop over balance plots
    {
      for(int a = 0; a < 2; a++)
	{
	  Beta[i+a]->SetMarkerStyle(markerStyle[a]);
	  Beta[i+a]->SetMarkerColor(markerColor[a]);
	  Beta[i+a]->SetLineColor(markerColor[a]);
	  objToBeWritten.push_back(Beta[i+a]);
	  
	  Fit2D(Beta[i+a],hists_beta[i+a],gp_beta[i+a], gf_beta[i+a]);

	  for(int b = 0 ; b < 4 ; ++b)
	    {
	      hists_beta[i+a][b]->SetMinimum(-0.5);
	      hists_beta[i+a][b]->SetMaximum(0.5);
	      ++b;
	      hists_beta[i+a][b]->SetMinimum(0.0);
	      hists_beta[i+a][b]->SetMaximum(1.);
	    }
	}

      for(int a = 1; a >=0; a--) // Loop over correction
	{
	  for(int b = 0; b < 3; b++) // Loop over example eta bins
	    {
	      // Find eta bin of gaussplot
	      int bin = 0;
	      if(  b == 0  )  bin = int(Beta[i+a]->GetNbinsX()/6);
	      else if(  b == 1  )  bin = int(Beta[i+a]->GetNbinsX()/3);
	      else if(  b == 2  )  bin = int(Beta[i+a]->GetNbinsX()/2);
	      float min = Beta[i+a]->GetXaxis()->GetBinLowEdge(bin);
	      float max = min + Beta[i+a]->GetXaxis()->GetBinWidth(bin);

	      // Set title according to eta bin
	      if( i == 0 )  sprintf(name,"di-jet, %.2f < #eta < %.2f",min,max);
	      else sprintf(name,"di-jet, %i < E_{T}^{scale} < %i GeV, %.2f < #eta < %.2f",
			   etLimit[int(i/2)-1],etLimit[int(i/2)],min,max);
	      gp_beta[a][b]->SetTitle(name);

	      if( a == 0 ) gp_beta[a][b]->SetXTitle("B after fit");
	      else gp_beta[a][b]->SetXTitle("B before fit");

	      // Set style and line color according to ptRatioName
	      gp_beta[a][b]->SetMarkerStyle(markerStyle[a]);
	      gp_beta[a][b]->SetMarkerColor(markerColor[a]);
	      gp_beta[a][b]->SetLineColor(markerColor[a]);
	      gf_beta[a][b]->SetLineColor(markerColor[a]);

	      // Plot gaussplots
	      c2->cd(1+b);
	      gp_beta[a][b]->Draw();
	      gf_beta[a][b]->Draw("same");

	      objToBeWritten.push_back(gp_beta[a][b]);
	      objToBeWritten.push_back(gf_beta[a][b]);
	    } // End of loop over example ptbins
	  c2->Draw();
	  ps->NewPage();
	} // End of loop over correction
      c1->cd();

      for(int j = 0 ; j < 8 ; ++j) 
	{
	  hists_beta[i][j]->Draw("P");
	  objToBeWritten.push_back(hists_beta[i][j]);
	  hists_beta[i][j]->SetStats(0);
	  hists_beta[i+1][j]->Draw("P SAME");
	  objToBeWritten.push_back(hists_beta[i+1][j]);
	  leg->Draw("SAME");
	  c1->SetGrid();
	  c1->Draw();   
	  ps->NewPage(); 
	}
    } // End of loop over balance plots 



  // Control quantities vs pt
  TH1F* hists_pt[2][8];
  TH1F* gp_pt[2][4];  
  TF1* gf_pt[2][4];     
  for(int a = 0; a < 2; a++)
    {
      Bpt[a]->SetMarkerStyle(markerStyle[a]);
      Bpt[a]->SetMarkerColor(markerColor[a]);
      Bpt[a]->SetLineColor(markerColor[a]);

      Fit2D(Bpt[a],hists_pt[a],gp_pt[a], gf_pt[a]);

      for(int b = 0 ; b < 4 ; ++b)
	{
	  hists_pt[a][b]->SetMinimum(-0.5);
	  hists_pt[a][b]->SetMaximum(0.5);
	  ++b;
	  hists_pt[a][b]->SetMinimum(0.0);
	  hists_pt[a][b]->SetMaximum(1.);
	}
    }

  for(int a = 1; a >=0; a--) // Loop over correction
    {
      for(int b = 0; b < 3; b++) // Loop over example pt bins
	{
	  // Find pt bin of gaussplot
	  int bin = 0;
	  if(  b == 0  )  bin = int(Bpt[a]->GetNbinsX()/6);
	  else if(  b == 1  )  bin = int(Bpt[a]->GetNbinsX()/3);
	  else if(  b == 2  )  bin = int(Bpt[a]->GetNbinsX()/2);
	  float min = Bpt[a]->GetXaxis()->GetBinLowEdge(bin);
	  float max = min + Bpt[a]->GetXaxis()->GetBinWidth(bin);

	  // Set title according to pt bin
	  sprintf(name,"di-jet, %.1f < p_{T} < %.1f [GeV]",min,max);
	  gp_pt[a][b]->SetTitle(name);

	  if( a == 0 ) gp_pt[a][b]->SetXTitle("B after fit");
	  else gp_pt[a][b]->SetXTitle("B before fit");

	  // Set style and line color according to ptRatioName
	  gp_pt[a][b]->SetMarkerStyle(markerStyle[a]);
	  gp_pt[a][b]->SetMarkerColor(markerColor[a]);
	  gp_pt[a][b]->SetLineColor(markerColor[a]);
	  gf_pt[a][b]->SetLineColor(markerColor[a]);

	  // Plot gaussplots
	  c2->cd(1+b);
	  gp_pt[a][b]->Draw();
	  gf_pt[a][b]->Draw("same");

	  objToBeWritten.push_back(gp_pt[a][b]);
	  objToBeWritten.push_back(gf_pt[a][b]);
	} // End of loop over example ptbins
      c2->Draw();
      ps->NewPage();
    } // End of loop over correction

  c1->cd();
  for(int j = 0 ; j < 8 ; ++j) 
    {
      hists_pt[0][j]->Draw("P");
      objToBeWritten.push_back(hists_pt[0][j]);
      hists_pt[0][j]->SetStats(0);
      hists_pt[1][j]->Draw("P SAME");
      objToBeWritten.push_back(hists_pt[1][j]);
      leg->Draw("SAME");
      c1->SetGrid();
      c1->Draw();   
      ps->NewPage(); 
    }



  // Control quantities vs log pt
  TH1F* hists_ptlog[2][8];
  TH1F* gp_ptlog[2][4];  
  TF1* gf_ptlog[2][4];     
  for(int a = 0; a < 2; a++)
    {
      Bptlog[a]->SetMarkerStyle(markerStyle[a]);
      Bptlog[a]->SetMarkerColor(markerColor[a]);
      Bptlog[a]->SetLineColor(markerColor[a]);

      Fit2D(Bptlog[a],hists_ptlog[a],gp_ptlog[a], gf_ptlog[a]);

      for(int b = 0 ; b < 4 ; ++b)
	{
	  hists_ptlog[a][b]->SetMinimum(-0.5);
	  hists_ptlog[a][b]->SetMaximum(0.5);
	  ++b;
	  hists_ptlog[a][b]->SetMinimum(0.0);
	  hists_ptlog[a][b]->SetMaximum(1.);
	}
    }

  c1->cd();
  for(int j = 0 ; j < 8 ; ++j) 
    {
      hists_ptlog[0][j]->Draw("P");
      objToBeWritten.push_back(hists_ptlog[0][j]);
      hists_ptlog[0][j]->SetStats(0);
      hists_ptlog[1][j]->Draw("P SAME");
      objToBeWritten.push_back(hists_ptlog[1][j]);
      leg->Draw("SAME");
      c1->SetGrid();
      c1->SetLogx(1);
      c1->Draw();   
      ps->NewPage(); 
    }

  c1->SetLogx(0);

  // Control quantities vs energy
  TH1F* hists_energy[2][8];
  TH1F* gp_energy[2][4];  
  TF1* gf_energy[2][4];     
  for(int a = 0; a < 2; a++)
    {
      Benergy[a]->SetMarkerStyle(markerStyle[a]);
      Benergy[a]->SetMarkerColor(markerColor[a]);
      Benergy[a]->SetLineColor(markerColor[a]);

      Fit2D(Benergy[a],hists_energy[a],gp_energy[a], gf_energy[a]);

      for(int b = 0 ; b < 4 ; ++b)
	{
	  hists_energy[a][b]->SetMinimum(-0.5);
	  hists_energy[a][b]->SetMaximum(0.5);
	  ++b;
	  hists_energy[a][b]->SetMinimum(0.0);
	  hists_energy[a][b]->SetMaximum(1.);
	}
    }

  for(int a = 1; a >=0; a--) // Loop over correction
    {
      for(int b = 0; b < 3; b++) // Loop over example energy bins
	{
	  // Find energy bin of gaussplot
	  int bin = 0;
	  if(  b == 0  )  bin = int(Benergy[a]->GetNbinsX()/6);
	  else if(  b == 1  )  bin = int(Benergy[a]->GetNbinsX()/3);
	  else if(  b == 2  )  bin = int(Benergy[a]->GetNbinsX()/2);
	  float min = Benergy[a]->GetXaxis()->GetBinLowEdge(bin);
	  float max = min + Benergy[a]->GetXaxis()->GetBinWidth(bin);

	  // Set title according to energy bin
	  sprintf(name,"di-jet, %.1f < E < %.1f [GeV]",min,max);
	  gp_energy[a][b]->SetTitle(name);

	  if( a == 0 ) gp_energy[a][b]->SetXTitle("B after fit");
	  else gp_energy[a][b]->SetXTitle("B before fit");

	  // Set style and line color according to ptRatioName
	  gp_energy[a][b]->SetMarkerStyle(markerStyle[a]);
	  gp_energy[a][b]->SetMarkerColor(markerColor[a]);
	  gp_energy[a][b]->SetLineColor(markerColor[a]);
	  gf_energy[a][b]->SetLineColor(markerColor[a]);

	  // Plot gaussplots
	  c2->cd(1+b);
	  gp_energy[a][b]->Draw();
	  gf_energy[a][b]->Draw("same");

	  objToBeWritten.push_back(gp_energy[a][b]);
	  objToBeWritten.push_back(gf_energy[a][b]);
	} // End of loop over example energy bins
      c2->Draw();
      ps->NewPage();
    } // End of loop over correction

  c1->cd();
  for(int j = 0 ; j < 8 ; ++j) 
    {
      hists_energy[0][j]->Draw("P");
      objToBeWritten.push_back(hists_energy[0][j]);
      hists_energy[0][j]->SetStats(0);
      hists_energy[1][j]->Draw("P SAME");
      objToBeWritten.push_back(hists_energy[1][j]);
      leg->Draw("SAME");
      c1->SetGrid();
      c1->Draw();   
      ps->NewPage(); 
    }


  // Control quantities vs emf
  TH1F* hists_emf[2][8];
  TH1F* gp_emf[2][4];  
  TF1* gf_emf[2][4];     
  for(int a = 0; a < 2; a++)
    {
      Bemf[a]->SetMarkerStyle(markerStyle[a]);
      Bemf[a]->SetMarkerColor(markerColor[a]);
      Bemf[a]->SetLineColor(markerColor[a]);

      Fit2D(Bemf[a],hists_emf[a],gp_emf[a], gf_emf[a]);

      for(int b = 0 ; b < 4 ; ++b)
	{
	  hists_emf[a][b]->SetMinimum(-0.5);
	  hists_emf[a][b]->SetMaximum(0.5);
	  ++b;
	  hists_emf[a][b]->SetMinimum(0.0);
	  hists_emf[a][b]->SetMaximum(1.);
	}
    }

  for(int a = 1; a >=0; a--) // Loop over correction
    {
      for(int b = 0; b < 3; b++) // Loop over example emf bins
	{
	  // Find emf bin of gaussplot
	  int bin = 0;
	  if(  b == 0  )  bin = int(Bemf[a]->GetNbinsX()/6);
	  else if(  b == 1  )  bin = int(Bemf[a]->GetNbinsX()/3);
	  else if(  b == 2  )  bin = int(Bemf[a]->GetNbinsX()/2);
	  float min = Bemf[a]->GetXaxis()->GetBinLowEdge(bin);
	  float max = min + Bemf[a]->GetXaxis()->GetBinWidth(bin);

	  // Set title according to energy bin
	  sprintf(name,"di-jet, %.2f < f_{em} < %.2f",min,max);
	  gp_emf[a][b]->SetTitle(name);

	  if( a == 0 ) gp_emf[a][b]->SetXTitle("B after fit");
	  else gp_emf[a][b]->SetXTitle("B before fit");

	  // Set style and line color according to ptRatioName
	  gp_emf[a][b]->SetMarkerStyle(markerStyle[a]);
	  gp_emf[a][b]->SetMarkerColor(markerColor[a]);
	  gp_emf[a][b]->SetLineColor(markerColor[a]);
	  gf_emf[a][b]->SetLineColor(markerColor[a]);

	  // Plot gaussplots
	  c2->cd(1+b);
	  gp_emf[a][b]->Draw();
	  gf_emf[a][b]->Draw("same");

	  objToBeWritten.push_back(gp_emf[a][b]);
	  objToBeWritten.push_back(gf_emf[a][b]);
	} // End of loop over example energy bins
      c2->Draw();
      ps->NewPage();
    } // End of loop over correction

  c1->cd();
  for(int j = 0 ; j < 8 ; ++j) 
    {
      hists_emf[0][j]->Draw("P");
      objToBeWritten.push_back(hists_emf[0][j]);
      hists_emf[0][j]->SetStats(0);
      hists_emf[1][j]->Draw("P SAME");
      objToBeWritten.push_back(hists_emf[1][j]);
      leg->Draw("SAME");
      c1->SetGrid();
      c1->Draw();   
      ps->NewPage(); 
    }
  

//   TF1* line = new TF1("line","x",0,200);
//   for(int i=0;i<5;++i)
//     {
//       combmean[i]->Draw("Box");
//       line->Draw("same");
//       c1->Draw();   
//       ps->NewPage();
//     }
//   for(int i=0;i<5;++i)
//     {
//       difmean[i]->Draw("Box");
//       line->Draw("same");
//       c1->Draw();   
//       ps->NewPage();
//     }
//   delete line;



  // Clean up
  ps->Close();

  if( _outputROOT ) WriteToRootFile( objToBeWritten, "DiJet" );

  delete c1;
  delete c2;
  delete ps;
  delete Scale[0];
  delete Scale[1];
  for(int i = 0; i < 5; i++)
    {
      if( i < 4 ) delete ptspec[i];
      delete combmean[i];
      delete difmean[i];
    }
  for(int i = 0; i <11; i++)
    {
      delete dphi[i];
    }
  delete Bvsdphi;
  delete Difvscomb;
  for(int i = 0; i < 8; i++)
    {
      delete Beta[i];
    }
  for(int i = 0; i < 2; i++)
    {
      delete Bpt[i];
      delete Benergy[i];
      delete Bemf[i];
      delete Bptlog[i];
    }
  for(int i = 0 ; i < 8 ; ++i)
    {
      for(int j = 0 ; j < 4 ; ++j)
	{
	  delete hists_beta[i][j];
	  delete hists_beta[i][j+4];
	  delete gp_beta[i][j];
	  delete gf_beta[i][j];
	}	
    }
  for(int i = 0 ; i < 2 ; ++i)
    {
      for(int j = 0 ; j < 4 ; ++j)
	{
	  delete hists_pt[i][j];
	  delete hists_pt[i][j+4];
	  delete gp_pt[i][j];
	  delete gf_pt[i][j];

	  delete hists_ptlog[i][j];
	  delete hists_ptlog[i][j+4];
	  delete gp_ptlog[i][j];
	  delete gf_ptlog[i][j];

	  delete hists_energy[i][j];
	  delete hists_energy[i][j+4];
	  delete gp_energy[i][j];
	  delete gf_energy[i][j];

	  delete hists_emf[i][j];
	  delete hists_emf[i][j+4];
	  delete gp_emf[i][j];
	  delete gf_emf[i][j];
	}	
    }
  delete leg;
}





//---------------------------------------------------------------
//   Takes a 2D histogram 'hist' and creates projections along
//   the x-axis per x bin. Some properties of these projected
//   distributions are filled, per x-bin, into 8 1D histograms
//   'hresuslts':
//      0: Mean value
//      1: Standard deviation
//      2: Mean of Gauss fit
//      3: Width of Gauss fit
//      4: Median 
//      5: chi2 / n.d.f.
//      6: Probability of Gauss fit
//      7: Quantiles Q0.9 / (Q0.9 - 1)
//   'hresuslts' are newly created (take care of deleting them!);
//   their object names are set to:
//      "<hist-name>_result<X>",
//   where <hist-name> = hist->GetName() and <X> is the index of
//   the above specified property (i.e. 0 for "mean").
//
//   Also, the projected distributions of 3 example x-bins:
//      0: hist->GetNbinsX() / 6
//      1: hist->GetNbinsX() / 3
//      2: hist->GetNbinsX() / 2
//   are filled into 'gaussplots' and the corresponding
//   Gauss fits are filled into 'gf'. Both 'gaussplots' and 'gf'
//   are newly created (take care of deleting them!) and their
//   names are set to:
//      "<hist-name>_gaussplot<X>",
//      "<hist-name>_gaussfit<X>"
//   respectively.
//---------------------------------------------------------------
void TControlPlots::Fit2D(const TH2F* hist, TH1F* hresults[8], TH1F* gaussplots[4], TF1* gf[4] ) const
{
  //book hists
  TString s = hist->GetName();
  s += "_result0";
  if( hist->GetXaxis()->GetXbins()->GetSize() == hist->GetNbinsX() +1)
    {
      hresults[0] = new TH1F(s,"",hist->GetNbinsX(),hist->GetXaxis()->GetXbins()->GetArray());
    }
  else
    {
      hresults[0] = new TH1F(s,"",hist->GetNbinsX(),hist->GetXaxis()->GetXmin(),
			   hist->GetXaxis()->GetXmax());
    }
  hresults[0]->SetXTitle(hist->GetXaxis()->GetTitle());
  hresults[0]->SetMarkerStyle(hist->GetMarkerStyle());
  hresults[0]->SetMarkerColor(hist->GetMarkerColor());
  hresults[0]->SetLineColor(hist->GetLineColor());
  hresults[0]->SetMarkerSize(hist->GetMarkerSize());
  for(int i = 1; i < 8 ; ++i)
    {
      s = hist->GetName();
      s += "_result";
      s += i;
      hresults[i] = (TH1F*) hresults[0]->Clone(s);
      s = hist->GetTitle();
      hresults[i]->SetTitle(s + ",  " + controlQuantityName[i]); 
    }
  s = hist->GetTitle();
  hresults[0]->SetTitle(s + ",  " + controlQuantityName[0]); 

  hresults[5]->SetMinimum(0.0);
  hresults[5]->SetMaximum(100);
  hresults[6]->SetMinimum(0.0);
  hresults[6]->SetMaximum(1.05);


  TH1F* htemp = new TH1F("htemp","",hist->GetNbinsY(),hist->GetYaxis()->GetXmin(),
		          hist->GetYaxis()->GetXmax());
  htemp->Sumw2();

  for(int i=0;i<4;++i) 
    {
      s = hist->GetName();
      s += "_gaussplot";
      s += i;
      gaussplots[i] = (TH1F*)htemp->Clone(s);
      
      s = hist->GetName();
      s += "_gaussfit";
      s += i;
      gf[i] = new TF1(s,"0");
    }


  const int nq = 2;
  double yq[2],xq[2];
  xq[0] = 0.5;
  xq[1] = 0.90;
  int index = 0;		// Counting index used for gaussplots
  for(int i = 1 ; i <= hist->GetNbinsX() ; ++i)
    {
      htemp->Reset();
      for(int j = 1 ; j <= hist->GetNbinsY() ; ++j)
	{
	  htemp->SetBinContent(j,hist->GetBinContent(hist->GetBin(i,j)));
	  htemp->SetBinError(j,hist->GetBinError(i,j));
	}  
      if(htemp->GetSumOfWeights() <= 0) continue;
      htemp->Fit("gaus","LLQNO","");
      TF1 *f = (TF1*)gROOT->GetFunction("gaus")->Clone();
      double mean = f->GetParameter(1);
      double meanerror = f->GetParError(1);
      double width = f->GetParameter(2);
      if(width < 0.2) width = 0.2;
      if( (htemp->Fit(f,"LLQNO","goff",mean - 2 * width, mean + 2 * width) == 0) && (f->GetProb() > 0.01) ) 
	{
	  mean = f->GetParameter(1);
	  meanerror = f->GetParError(1);
	  width = f->GetParameter(2);

	  hresults[2]->SetBinContent(i,mean);
	  hresults[2]->SetBinError(i,meanerror);
	  hresults[3]->SetBinContent(i,width);
	  hresults[3]->SetBinError(i, f->GetParError(2));
	}
      hresults[5]->SetBinContent(i, f->GetChisquare() / f->GetNumberFreeParameters());
      hresults[5]->SetBinError(i, 0.01);
      hresults[6]->SetBinContent(i, f->GetProb());
      hresults[6]->SetBinError(i, 0.01);

      if(  i == int(hist->GetNbinsX()/6)
	   || i == int(hist->GetNbinsX()/3)
	   ||  i == int(hist->GetNbinsX()/2)  )       
	{
	  gaussplots[index] = (TH1F*)htemp->Clone(gaussplots[index]->GetName());
	  gf[index] = (TF1*)f->Clone(gf[index]->GetName());
	  index++;
	}

      mean = htemp->GetMean();
      meanerror = htemp->GetMeanError();
      width = htemp->GetRMS();
      hresults[0]->SetBinContent(i,mean);
      hresults[0]->SetBinError(i,meanerror);
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




//---------------------------------------------------------------
// Write all TObjects in 'obj' to the file '_outFile' into the
// directory '_outFile:/dir'. If '_outFile:/dir' does not exist,
// it is created first.
//---------------------------------------------------------------
void TControlPlots::WriteToRootFile(std::vector<TObject*> obj, std::string dir)
{
  std::string directory = _outFile->GetName();
  directory += ":";
  gDirectory->cd(directory.c_str());
  directory += "/";
  directory += dir;
  bool dirExists = gDirectory->GetDirectory(directory.c_str());
  if( !dirExists )
    {
      gDirectory->mkdir(dir.c_str());
    }
  gDirectory->cd(directory.c_str());
  for(std::vector<TObject*>::const_iterator it = obj.begin(); it < obj.end(); it++)
    {
      int ok = gDirectory->WriteTObject( *it );
      if( !ok ) std::cerr << "Error writing object '" << (*it)->GetName() << "' to file." << std::endl;
    }
}



//---------------------------------------------------------------
// Set style option for ps output.
//---------------------------------------------------------------
void TControlPlots::SetGStyle() const
{
  // For the canvas:
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetCanvasColor(kWhite);
  gStyle->SetCanvasDefH(800); //Height of canvas
  gStyle->SetCanvasDefW(800); //Width of canvas
  gStyle->SetCanvasDefX(0);   //Position on screen
  gStyle->SetCanvasDefY(0);

  // For the frame
  gStyle->SetFrameBorderMode(0);
  gStyle->SetFrameBorderSize(1);
  gStyle->SetFrameFillColor(kBlack);
  gStyle->SetFrameFillStyle(0);
  gStyle->SetFrameLineColor(kBlack);
  gStyle->SetFrameLineStyle(0);
  gStyle->SetFrameLineWidth(1);

  // For the Pad:
  gStyle->SetPadBorderMode(0);
  gStyle->SetPadColor(kWhite);
  gStyle->SetPadGridX(false);
  gStyle->SetPadGridY(false);
  gStyle->SetGridColor(0);
  gStyle->SetGridStyle(3);
  gStyle->SetGridWidth(1);

  // For the histo:
  gStyle->SetHistLineColor(kBlack);
  gStyle->SetHistLineStyle(0);
  gStyle->SetHistLineWidth(1);

  // For the statistics box:
  gStyle->SetOptStat("neMR");
  gStyle->SetStatColor(kWhite);
  gStyle->SetStatFont(42);
  gStyle->SetStatFontSize(0.03);
  gStyle->SetStatTextColor(1);
  gStyle->SetStatFormat("6.4g");
  gStyle->SetStatBorderSize(1);
  gStyle->SetStatX(0.92);              
  gStyle->SetStatY(0.86);              
  gStyle->SetStatH(0.16);
  gStyle->SetStatW(0.22);

  // For the legend
  gStyle->SetLegendBorderSize(1);

  // Margins:
  gStyle->SetPadTopMargin(0.10);
  gStyle->SetPadBottomMargin(0.14);
  gStyle->SetPadLeftMargin(0.2);
  gStyle->SetPadRightMargin(0.04);

  // For the Global title:
  gStyle->SetOptTitle(1);
  gStyle->SetTitleFont(42,"");
  gStyle->SetTitleColor(1);
  gStyle->SetTitleTextColor(1);
  gStyle->SetTitleFillColor(0);
  gStyle->SetTitleFontSize(0.1);
  gStyle->SetTitleAlign(23);
  gStyle->SetTitleX(0.58);
  gStyle->SetTitleH(0.05);
  gStyle->SetTitleXOffset(0);
  gStyle->SetTitleYOffset(0);
  gStyle->SetTitleBorderSize(0);

  // For the axis titles:
  gStyle->SetTitleColor(1,"XYZ");
  gStyle->SetTitleFont(42,"XYZ");
  gStyle->SetTitleSize(0.04,"XYZ");
  gStyle->SetTitleXOffset(1.5);
  gStyle->SetTitleYOffset(2.0);

  // For the axis labels:
  gStyle->SetLabelColor(1,"XYZ");
  gStyle->SetLabelFont(42,"XYZ");
  gStyle->SetLabelOffset(0.007,"XYZ");
  gStyle->SetLabelSize(0.04,"XYZ");

  // For the axis:
  gStyle->SetAxisColor(1,"XYZ");
  gStyle->SetStripDecimals(kTRUE);
  gStyle->SetTickLength(0.03,"XYZ");
  gStyle->SetNdivisions(510,"XYZ");
  gStyle->SetPadTickX(1);  // To get tick marks on the opposite side of the frame
  gStyle->SetPadTickY(1);
}


