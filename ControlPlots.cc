#include "ControlPlots.h"

#include <algorithm>
#include <iostream>
using namespace std;

#include "TCanvas.h"
#include "TDirectory.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TH1I.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TLorentzVector.h"
#include "TPaveText.h"
#include "TPostScript.h"
#include "TProfile.h"
#include "TString.h"
#include "TROOT.h"

#include "CalibMath.h"

#include "JetTruthEvent.h"
#include "TwoJetsInvMassEvent.h"


//!  \brief Constructor
//! 
//!  \param configfile  Path to config file
//!  \param data        TData objects to create controlplots from
//!  \param par         Parameters
// -------------------------------------------------------------
TControlPlots::TControlPlots(const std::string& configfile, const std::vector<TData*> *data, TParameters *par)
  : mData(data), mPar(par), mConfig(new ConfigFile(configfile.c_str())), mOutFile(0)
{ 
  mPtRatioName[0] = "p^{jet}_{T}/ E_{T}^{#gamma}";
  mPtRatioName[1] = "p_{T}^{cor. jet}/E_{T}^{#gamma}";
  mPtRatioName[2] = "p_{T}^{jet}/p_{T}^{cor. jet}";

  mControlQuantityName[0] = "mean"; 
  mControlQuantityName[1] = "standard deviation"; 
  mControlQuantityName[2] = "mean of Gauss fit"; 
  mControlQuantityName[3] = "width of Gauss fit"; 
  mControlQuantityName[4] = "median"; 
  mControlQuantityName[5] = "#chi^{2} / n.d.f."; 
  mControlQuantityName[6] = "probability"; 
  mControlQuantityName[7] = "quantiles"; 

  SetGStyle();

  if(mConfig->read<bool>("plot output format",0)) {
    mOutputROOT = false;
  }
  else {
    mOutputROOT = true;
    mOutFile = new TFile("controlplots.root","RECREATE","Cal calib control plots");
  }    
}



// -------------------------------------------------------------
TControlPlots::~TControlPlots()
{
  if( mOutFile !=0 )
    {
      if( mOutFile->IsOpen() ) mOutFile->Close();
      delete mOutFile;
    }
}



//!  \brief Create all plots as specified in the config file
// -------------------------------------------------------------
void TControlPlots::MakePlots()
{
  cout << endl << "Writing control plots in .ps " << flush;
  if( mOutputROOT ) cout << "and .root " << flush;
  cout << "format:" << endl;
  if( mConfig->read<bool>("create jet-truth event plots",false) ) {
    cout << "Creating jet-truth event control plots... " << flush;
    MakeControlPlotsJetTruthEventResponse();
    cout << "ok" << endl;
  }
  if( mConfig->read<bool>("create binned response plots",false) ) {
    cout << "Creating binned control plots... " << flush;
    MakeControlPlotsBinnedResponse();
    cout << "ok" << endl;
  }
  if( mConfig->read<bool>("create L2L3 MC truth plots",false) ) {
    cout << "Creating L2L3 MC truth control plots... " << flush;
    MakeControlPlotsL2L3MCTruth();
    cout << "ok" << endl;
  }
  if( mConfig->read<bool>("create chi2 plots",false) ) {
    cout << "Creating chi2 control plots... " << flush;
    MakeControlPlotsChi2();
    cout << "ok" << endl;
  }
  if( mConfig->read<bool>("create correction function plots",false) ) {
    cout << "Creating correction function control plots... " << flush;
    MakeControlPlotsCorrectionFunction();
    cout << "ok" << endl;
  }
  if( mConfig->read<bool>("create tower plots" ,false) ) {
    cout << "Creating tower control plots... " << flush;
    MakeControlPlotsTowers();
    cout << "ok" << endl;
  }
  if( mConfig->read<bool>("create gamma jet plots" ,false) ) {
    cout << "Creating more gamma jet control plots... " << flush;
    MakeControlPlotsGammaJet();
    cout << "ok" << endl;
  }
  if( mConfig->read<bool>("create more gamma jet plots" ,false) ) {
    cout << "Creating gamma jet (tower bin) control plots... " << flush;
    MakeControlPlotsGammaJetPerTowerBin();
    cout << "ok" << endl;
      
    cout << "Creating gamma jet (jet bin) control plots... " << flush;
    MakeControlPlotsGammaJetPerJetBin();
    cout << "ok" << endl;
    
    cout << "Creating even more gamma jet control plots (sigmas)... " << flush;
    MakeControlPlotsGammaJetSigmas();
    cout << "ok" << endl;
  }
  if( mConfig->read<bool>("create dijet plots",false) ) {
    cout << "Creating di-jet control plots... " << flush;
    MakeControlPlotsDiJet();
    cout << "ok" << endl;
  }
  if( mConfig->read<bool>("create top plots",false) ) {
    cout << "Creating top control plots... " << flush;
    MakeControlPlotsTop();
    cout << "ok" << endl;
  }
  if( mConfig->read<bool>("create parameter scan plots" ,false) ) {
    cout << "Creating parameter scan control plots... " << flush;
    MakeControlPlotsParameterScan();
    cout << "ok" << endl;
  }
}



//!  \brief Response distribution in bins of truth Et and eta
//!
//!  Creates distributions of the response
//!  \f$ E^{jet}_{T} / E^{true}_{T} \f$
//!  in bins of true Et \f$ E^{true}_{T} \f$ and measured
//!  pseudorapidity \f$ \eta \f$, where \f$ E^{jet}_{T} \f$
//!  is the corrected jet Et.
//!
//!  There are two sets of distributions:
//!   - Correction from global fit
//!   - Comparison of correction from global fit
//!     and JetMET L2L3 correction
//!
//!  The plots are written to the file
//!  "controlplotsBinnedResponse.ps" and (if enabled)
//!  to the directory "BinnedResponse" in the file
//!  "controlplots.root".
// -------------------------------------------------------------
void TControlPlots::MakeControlPlotsBinnedResponse()
{
  std::vector<TObject*>                objToBeWritten;
  std::vector<TData*>::const_iterator  data_it;


  // Create a Et and eta binning
  // Et = x, eta = y
  std::vector<double> binEdgesEt
    = bag_of<double>(mConfig->read<std::string>("Control plots pt bin edges","0 100"));
  std::vector<double> binEdgesEta
    = bag_of<double>(mConfig->read<std::string>("Control plots eta bin edges","-5 5"));
  Binning bins(binEdgesEt,binEdgesEta);


  // Create histograms of response in
  // Et and eta bins
  std::vector<TH1F*> hResp;
  std::vector<TH1F*> hRespJetMet;
  for(int i = 0; i < bins.NBins(); i++)
    {
      char name[50];
      sprintf(name,"hResponse_pt%i_eta%i",bins.IX(i),bins.IY(i));
      TH1F * h = new TH1F(name,"",51,0,2);
      h->SetLineColor(2);
      h->GetXaxis()->SetTitleSize(0.05);
      h->GetXaxis()->SetTitleOffset(1.2);
      h->GetYaxis()->SetTitleSize(0.05);
      h->GetYaxis()->SetTitleOffset(1.6);
      hResp.push_back(h);
      objToBeWritten.push_back(h);

      sprintf(name,"hResponseJetMET_pt%i_eta%i",bins.IX(i),bins.IY(i));
      h = static_cast<TH1F*>(h->Clone(name));
      h->SetLineColor(1);
      hRespJetMet.push_back(h);
      objToBeWritten.push_back(h);
    }


  // Loop over data and fill response into
  // the corresponding histogram
  for(  data_it = mData->begin();  data_it != mData->end();  data_it++ )
    {
      double weight = (*data_it)->GetWeight();
      double etmeas = (*data_it)->GetMess()->pt;
      double etcorr = (*data_it)->GetParametrizedMess();
      double eta    = (*data_it)->GetMess()->eta;
      double ettrue = (*data_it)->GetTruth();

      // Discard events flagged bad
      JetTruthEvent *jte = dynamic_cast<JetTruthEvent*>(*data_it);
      if( jte )
	{
	  if( jte->FlaggedBad() ) continue;
	}

      int bin = bins.Bin(ettrue,eta);
      if( 0 <= bin && bin < bins.NBins() )
	{
	  hResp.at(bin)->Fill(etcorr/ettrue,weight);

	  // For comparison reasons, correct with JetMET
	  // correction; this works only for data class > 0
	  TJet *jet = dynamic_cast<TJet*>((*data_it)->GetMess());
	  if( jet )
	    {
	      double cjetmet = jet->L2L3cor;
	      hRespJetMet.at(bin)->Fill(cjetmet*etmeas/ettrue,weight);
	    }
	}
    }


  // Draw histograms into ps file, put 
  // at the most 6 histograms per page
  TPostScript * const ps = new TPostScript("controlplotsBinnedResponse.ps",112);
  ps->Range(25,1);

  double cw  = 300;
  double mw  = 70;

  double w   = 3*cw + mw;
  double h   = 2*cw + mw;
  
  double crw = cw/w;
  double mrw = mw/w;
  double crh = cw/h;
  double mrh = mw/h;

  TCanvas     * const c1 = new TCanvas("c1","",(int)w,(int)h);

  // Pads for the histograms
  std::vector<TPad*> cPads;
  cPads.push_back(new TPad("cPad0","",mrw,mrh+crh,mrw+crw,1.));
  cPads.push_back(new TPad("cPad1","",mrw+crw,mrh+crh,mrw+2*crw,1.));
  cPads.push_back(new TPad("cPad2","",mrw+2*crw,mrh+crh,1.,1.));
  cPads.push_back(new TPad("cPad3","",mrw,mrh,mrw+crw,mrh+crh));
  cPads.push_back(new TPad("cPad4","",mrw+crw,mrh,mrw+2*crw,mrh+crh));
  cPads.push_back(new TPad("cPad5","",mrw+2*crw,mrh,1.,mrh+crh));
  for(unsigned int i = 0; i < cPads.size(); i++)
    {
      cPads.at(i)->SetFillStyle(1001);
      cPads.at(i)->SetFrameFillColor(10);
      cPads.at(i)->SetFrameBorderMode(0);
      cPads.at(i)->SetTopMargin(0.04);
      cPads.at(i)->SetBottomMargin(0.1);
      cPads.at(i)->SetLeftMargin(0.1);
      cPads.at(i)->SetRightMargin(0.04);

      c1->cd();
      cPads.at(i)->Draw();
    }

  // Pads for margins holding titles
  TPad * mbPad = new TPad("mbPad","",0.,0.,1.,mrh);
  mbPad->SetFillStyle(1001);
  mbPad->SetFrameFillColor(10);
  mbPad->SetFrameBorderMode(0);
  TPaveText * bLabel = new TPaveText(0.8,0.5,0.95,1.,"NDC");
  bLabel->SetFillColor(0);
  bLabel->SetTextFont(42);
  bLabel->SetTextSize(0.7);
  bLabel->SetBorderSize(0);
  bLabel->AddText("p^{jet}_{T} / p^{true}_{T}");
  c1->cd();
  mbPad->Draw();
  mbPad->cd();
  bLabel->Draw();

  //   TPad * mlPad = new TPad("mlPad","",0.,mrh,mrw,1.);
  //   mlPad->SetFillStyle(1001);
  //   mlPad->SetFrameFillColor(10);
  //   mlPad->SetFrameBorderMode(0);
  //   TPaveText * lLabel = new TPaveText(0.1,0.5,0.95,1.,"NDC");
  //   lLabel->SetFillColor(0);
  //   lLabel->SetTextFont(42);
  //   lLabel->SetTextSize(0.7);
  //   lLabel->SetBorderSize(0);
  //   lLabel->SetTextAngle(90);
  //   lLabel->AddText("dN / d( E^{jet}_{T} / E^{true}_{T} )");
  //   c1->cd();
  //   mlPad->Draw();
  //   mlPad->cd();
  //   lLabel->Draw();

  // For some reason, this prevents the first ps-page
  // from looking weird...
  c1->Draw();
  ps->NewPage();
  for(unsigned int i = 0; i < hResp.size(); i++)
    {
      // Plot histogram
      int c = ( i % 6 );
      cPads.at(c)->cd();
      TH1F *h = hResp.at(i);
      h->GetYaxis()->SetRangeUser(0,1.6*(h->GetMaximum()));
      h->Draw();

      // Fit a gaussian in central part
      h->Fit("gaus","0QIL","",h->GetMean() - 1.5*(h->GetRMS()),h->GetMean() + 1.5*(h->GetRMS()));
      TF1 *fit = h->GetFunction("gaus");
      fit->Draw("same");

      // Draw fit parameters
      char label[100];
      TPaveText * fitstat = new TPaveText(0.13,0.66,0.93,0.93,"NDC");
      fitstat->SetFillColor(0);
      fitstat->SetTextFont(42);
      fitstat->SetTextAlign(12);
      sprintf(label,"%.0f < p^{true}_{T} < %.0f GeV, %.1f < #eta < %.1f",
	      bins.XLow(i),bins.XUp(i),bins.YLow(i),bins.YUp(i));
      fitstat->AddText(label);
      sprintf(label,"#mu = %.3f #pm %.3f",fit->GetParameter(1),fit->GetParError(1));
      fitstat->AddText(label);
      sprintf(label,"#sigma = %.3f #pm %.3f",fit->GetParameter(2),fit->GetParError(2));
      fitstat->AddText(label);
      fitstat->Draw("same");

      if( c == 5 )
	{
	  c1->Draw();
	  ps->NewPage();
	}
    }

  // Now the comparison plots to JetMET 
  TLegend *legComp = new TLegend(0.13,0.66,0.6,0.84);
  legComp->SetBorderSize(0);
  legComp->SetFillColor(0);
  legComp->SetTextFont(42);
  legComp->AddEntry(hResp.at(0),"Global fit","L");
  legComp->AddEntry(hRespJetMet.at(0),"JetMET L2L3","L");
  for(unsigned int i = 0; i < hResp.size(); i++)
    {
      int c = ( i % 6 );
      c1->cd();
      cPads.at(c)->cd();

      // Set maximum from histo with larger bin content
      double max = (hResp.at(i)->GetMaximum())/1.6; // Compensate for already setting this maximum above!
      if( hRespJetMet.at(i)->GetMaximum() > max ) max = hRespJetMet.at(i)->GetMaximum();
      hResp.at(i)->GetYaxis()->SetRangeUser(0,1.6*max);
      hResp.at(i)->Draw();
      hRespJetMet.at(i)->Draw("same");

      // Draw a vertical line at 1
      TLine *line = new TLine(1.,0.,1.,max);
      line->SetLineWidth(1);
      line->SetLineColor(4);
      line->SetLineStyle(2);
      line->Draw("same");

      // Label bin
      char label[100];
      TPaveText * fitstat = new TPaveText(0.13,0.84,0.93,0.93,"NDC");
      fitstat->SetFillColor(0);
      fitstat->SetTextFont(42);
      fitstat->SetTextAlign(12);
      sprintf(label,"%.0f < p^{true}_{T} < %.0f GeV, %.1f < #eta < %.1f",
	      bins.XLow(i),bins.XUp(i),bins.YLow(i),bins.YUp(i));
      fitstat->AddText(label);
      fitstat->Draw("same");

      if( c == 5 )
	{
	  legComp->Draw("same");
	  c1->Draw();
	  ps->NewPage();
	}
    }

  if( mOutputROOT ) WriteToRootFile(objToBeWritten,"BinnedResponse");


  // Clean up
  ps->Close();
  for(std::vector<TH1F*>::iterator it = hResp.begin();
      it != hResp.end(); it++)
    {
      delete *it;
    }
  hResp.clear();
  for(std::vector<TH1F*>::iterator it = hRespJetMet.begin();
      it != hRespJetMet.end(); it++)
    {
      delete *it;
    }
  hRespJetMet.clear();
  objToBeWritten.clear();
  delete c1;
  delete ps;
}



//!  \brief Chi2 distribution
//!
//!  Creates the following distributions of
//!  \f$ \chi^{2}_{i} \f$ from the sum
//!  \f[
//!   \chi^{2} = \sum_{i} w_{i}\chi^{2}_{i}
//!  \f]
//!  where \f$ w_{i} \f$ is the event weight
//!  factor:
//!   - The residual scaling scheme from the last 
//!     iteration is used.
//!   - All three residual scaling schemes are
//!     compared:
//!      - Distributions of scaled residuals
//!      - Scaled residuals vs residuals
//!   - The \f$ \chi^{2}_{i} \f$ probability
//!     prob(chi2,1).
//!
//!  The plots are written to the file
//!  "controlplotsChi2.ps" and (if enabled)
//!  to the directory "Chi2" in the file
//!  "controlplots.root".
// -------------------------------------------------------------
void TControlPlots::MakeControlPlotsChi2()
{
  std::vector<TObject*>                objToBeWritten;
  std::vector<TData*>::const_iterator  data_it;

  // Distribution of chi2 summands with last scaling
  // configuration scheme
  TH1F *h_chi2_last = new TH1F("h_chi2last",";#chi^{2}_{i};dN / d#chi^{2}_{i}",50,0,20);
  h_chi2_last->SetLineColor(1);
  h_chi2_last->SetLineStyle(1);
  objToBeWritten.push_back(h_chi2_last);

  // Distribution of chi2 probability with last scaling
  // configuration scheme
  TH1F *h_prob = new TH1F("h_prob",";prob(#chi^{2}_{i},1)",50,0,1);
  h_prob->SetMarkerStyle(7);
  objToBeWritten.push_back(h_prob);

  // Distribution of chi2 summands (normalized residuals)
  TH1F *h_chi2 = new TH1F("h_chi2","Scaled residuals z^{2} = 1/w #chi^{2};f(z^{2});dN / df(z^{2})",50,0,20);
  h_chi2->SetLineColor(1);
  h_chi2->SetLineStyle(1);
  objToBeWritten.push_back(h_chi2);

  // Distribution of Cauchy scaled normalized residuals
  TH1F *h_cauchy = static_cast<TH1F*>(h_chi2->Clone("h_cauchy"));
  h_cauchy->SetLineColor(2);
  h_cauchy->SetLineStyle(2);
  objToBeWritten.push_back(h_cauchy);

  // Distribution of Huber scaled normalized residuals
  TH1F *h_huber = static_cast<TH1F*>(h_cauchy->Clone("h_huber"));
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


  // Loop over data and fill histograms
  for(  data_it = mData->begin();  data_it != mData->end();  data_it++ )
    {
      double weight   = (*data_it)->GetWeight();
      double lastchi2 = ((*data_it)->chi2_plots())/weight;

      h_chi2_last->Fill(lastchi2);
      h_prob->Fill(TMath::Prob(lastchi2,1));

      TData::ScaleResidual = &TData::ScaleNone;
      double res = ( (*data_it)->chi2() / weight );
      TData::ScaleResidual = &TData::ScaleCauchy;
      double res_cauchy = ( (*data_it)->chi2() / weight );
      TData::ScaleResidual = &TData::ScaleHuber;
      double res_huber = ( (*data_it)->chi2() / weight );

      h_chi2->Fill(res);
      h_cauchy->Fill(res_cauchy);
      h_huber->Fill(res_huber);
      h_none_cauchy->Fill(res,res_cauchy);
      h_none_huber->Fill(res,res_huber);
    } // End loop over data


  // Draw histograms
  TPostScript * const ps = new TPostScript("controlplotsChi2.ps",111);
  TCanvas     * const c1 = new TCanvas("c1","",500,500);
  c1->cd();

  h_chi2_last->Draw();
  c1->SetLogy(1);
  c1->Draw();
  ps->NewPage();

  h_prob->GetYaxis()->SetRangeUser(0.,1.2*(h_prob->GetMaximum()));
  h_prob->Draw("PE1");
  c1->SetLogy(0);
  c1->Draw();
  ps->NewPage();

  h_cauchy->Draw();
  h_huber->Draw("same");
  h_chi2->Draw("same");
  c1->SetLogy(1);

  TLegend *l_res = new TLegend(0.35,0.68,0.7,0.88);
  l_res->SetFillColor(0);
  l_res->SetBorderSize(0);
  l_res->SetTextFont(42);
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
  l_res2->SetTextFont(42);
  l_res2->SetHeader("Scaling function f");
  l_res2->AddEntry(line1,"None","L");
  l_res2->AddEntry(h_none_huber,"Huber","L");
  l_res2->AddEntry(h_none_cauchy,"Cauchy","L");
  l_res2->Draw("same");

  c1->Draw();

  if( mOutputROOT ) WriteToRootFile(objToBeWritten, "Chi2");
  
  ps->Close();


  // Clean up
  delete l_res;
  delete line1;
  delete l_res2;
  delete c1;
  delete ps;
  for(std::vector<TObject*>::iterator it = objToBeWritten.begin();
      it != objToBeWritten.end(); it++)
    {
      delete *it;
    }
  objToBeWritten.clear();
}



//!  \brief Control plots for comparison with L2L3 MC truth correction
//!
//!  The plots are written to the file
//!  "controlplotsL2L3MCTruth.ps" and (if enabled)
//!  to the directory "L2L3MCTruth" in the file
//!  "controlplots.root".
// -------------------------------------------------------------
void TControlPlots::MakeControlPlotsL2L3MCTruth() {
  std::vector<TObject*>                objToBeWritten;
  std::vector<TData*>::const_iterator  data_it;


  // Absolute eta bins
  std::vector<double> etaBinEdge = bag_of<double>(mConfig->read<std::string>("Control plots eta bin edges",""));

  std::vector<double>::iterator itNeg = etaBinEdge.end();
  std::vector<double>::iterator it = etaBinEdge.begin();
  for(; it != etaBinEdge.end(); it++) {
    if( *it < 0 ) itNeg = it;
  }
  if( itNeg != etaBinEdge.end() ) etaBinEdge.erase(etaBinEdge.begin(),itNeg+1); // Want symmetry in eta

  // Find histogram ranges
  double ptGenMin = 10000.;
  double ptGenMax = 0.;
  for( data_it = mData->begin(); data_it != mData->end(); data_it++ ) {
    if( (*data_it)->GetTruth() < ptGenMin ) ptGenMin = (*data_it)->GetTruth();
    if( (*data_it)->GetTruth() > ptGenMax ) ptGenMax = (*data_it)->GetTruth();
  }
  ptGenMax *= 1.1;
  ptGenMin *= 0.9;


  // Create histograms:
  // Gen Pt in different eta bins
  std::vector<TH1F*> hGenPtSpec;
  for(unsigned int i = 1; i < etaBinEdge.size(); i++) {
    char name[50];
    sprintf(name,"hGenPtSpec_eta%i",i-1);
    TH1F * h = new TH1F(name,";p^{gen}_{T} (GeV);dN_{jet} / dp^{gen}_{T}  1 / (GeV)",
			50,ptGenMin,ptGenMax);
    h->GetXaxis()->SetNdivisions(505);
    h->SetMarkerStyle(19+i);
    int color = i;
    if( i >= 3 ) color++; // Avoid yellow
    h->SetMarkerColor(color);
    h->SetLineColor(color);
    hGenPtSpec.push_back(h);
    objToBeWritten.push_back(h);
  }


  // Loop over data and fill histograms
  for( data_it = mData->begin(); data_it != mData->end(); data_it++ ) {
    JetTruthEvent *jte = dynamic_cast<JetTruthEvent*>(*data_it);
    if( jte ) {
      if( jte->FlaggedBad() ) continue;    // Discard events flagged bad

      double weight = jte->GetWeight();
      double ptgen  = jte->GetTruth();
      double eta    = jte->GetMess()->eta;

      unsigned int etaBin = 0;
      if( fabs(eta) > etaBinEdge.back() ) continue;
      while( fabs(eta) > etaBinEdge.at(etaBin+1) ) etaBin++;
      hGenPtSpec.at(etaBin)->Fill( weight * ptgen );
    }
  }


  // Draw histograms
  TPostScript * const ps = new TPostScript("controlplotsL2L3MCTruth.ps",111);
  TCanvas     * const c1 = new TCanvas("c1","",500,500);
  c1->cd();

  double legminx = 0.6;
  double legminy = 0.65;
  TLegend * leg = new TLegend(legminx,legminy,0.9,0.85);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextFont(42);
  for(unsigned int i = 1; i < etaBinEdge.size(); i++) {
    char entry[50];
    sprintf(entry,"%.1f < |#eta| < %.1f",etaBinEdge.at(i-1),etaBinEdge.at(i));
    leg->AddEntry(hGenPtSpec.at(i-1),entry,"P");
  }

  hGenPtSpec.at(0)->GetYaxis()->SetRangeUser(0.1,15*(hGenPtSpec.at(0)->GetMaximum()));

  for(unsigned int i = 0; i < hGenPtSpec.size(); i++) {
    if( i == 0 ) hGenPtSpec.at(i)->Draw("PE1");
    else         hGenPtSpec.at(i)->Draw("PE1same");
  }
  c1->SetLogy(1);
  leg->Draw("same");
  c1->Draw();
  ps->NewPage();


  if( mOutputROOT ) WriteToRootFile(objToBeWritten, "L2L3MCTruth");
  
  ps->Close();


  // Clean up
  delete c1;
  delete ps;
  delete leg;
  std::vector<TObject*>::iterator objit = objToBeWritten.begin();
  for(; objit != objToBeWritten.end(); objit++) {
    delete *objit;
  }
  objToBeWritten.clear();
  hGenPtSpec.clear();
}



//!  \brief Response and resolution distribution in bins of
//!         truth Pt and eta
//!
//!  The plots are written to the file
//!  "controlplotsJetTruthEventResponse.ps",
//!  "controlplotsJetTruthEventResolution.ps", and (if enabled)
//!  to the directory "JetTruthEventResponse" in the file
//!  "controlplots.root".
// -------------------------------------------------------------
void TControlPlots::MakeControlPlotsJetTruthEventResponse() {
  std::vector<TObject*>                objToBeWritten;
  std::vector<TData*>::const_iterator  data_it;

  // -- Create 2D histograms of response vs eta, pt --------------

  // Create a pttrue and eta binning
  // pttrue = x, eta = y
  std::vector<double> binEdgesPt = bag_of<double>(mConfig->read<std::string>("Control plots pt bin edges",""));

  // Absolute eta bins
  std::vector<double> binEdgesEta = bag_of<double>(mConfig->read<std::string>("Control plots eta bin edges",""));
  Binning bins(binEdgesPt,binEdgesEta);
  //  bins.Print();


  // Create 2D histograms of response vs eta in
  // Pt and eta bins
  std::vector<TH2F*> h2EtaUncorr;    // Uncorrected response vs eta
  std::vector<TH2F*> h2EtaCorr;      // Response corrected by kalibri fit vs eta
  std::vector<TH2F*> h2EtaCorrL2;    // Response corrected by JetMET L2 correction vs eta
  std::vector<TH2F*> h2EtaCorrL2L3;  // Response corrected by JetMET L2L3 correction vs eta
  for(int ptbin = 0; ptbin < bins.NBinsX(); ptbin++) { // Loop over pttrue bins
    char name[50];
    sprintf(name,"h2EtaUncorr_pttrue%i",ptbin);
    TH2F * h2 = new TH2F(name,";#eta;< p^{jet}_{T} / p^{true}_{T} >",20,-5,5,51,0,2);
    h2EtaUncorr.push_back(h2);
    
    sprintf(name,"h2EtaCorr_pttrue%i",ptbin);
    h2EtaCorr.push_back(static_cast<TH2F*>(h2->Clone(name)));
    
    sprintf(name,"h2EtaCorrL2_pttrue%i",ptbin);
    h2EtaCorrL2.push_back(static_cast<TH2F*>(h2->Clone(name)));
    
    sprintf(name,"h2EtaCorrL2L3_pttrue%i",ptbin);
    h2EtaCorrL2L3.push_back(static_cast<TH2F*>(h2->Clone(name)));
  } // End of loop over pttrue bins

  // Create 2D histograms of response vs pttrue in
  // eta bins
  std::vector<TH2F*> h2PttrueUncorr;    // Uncorrected response vs pttrue
  std::vector<TH2F*> h2PttrueCorr;      // Response corrected by kalibri fit vs pttrue
  std::vector<TH2F*> h2PttrueCorrL2;    // Response corrected by JetMET L2 correction vs pttrue
  std::vector<TH2F*> h2PttrueCorrL2L3;  // Response corrected by JetMET L2L3 correction vs pttru
  // Logarithmic binning
  const int nLogBins = 15;
  double logBins[nLogBins+1];
  EquidistLogBins(logBins,nLogBins,bins.XLow(0),bins.XUp(bins.NBinsX()-1));
  for(int etabin = 0; etabin < bins.NBinsY(); etabin++) { // Loop over eta bins
    char name[50];
    sprintf(name,"h2PttrueUncorr_pttrue%i",etabin);
    TH2F * h2 = new TH2F(name,";p^{true}_{T} (GeV);< p^{jet}_{T} / p^{true}_{T} >",nLogBins,logBins,51,0,2);
    h2PttrueUncorr.push_back(h2);
    
    sprintf(name,"h2PttrueCorr_eta%i",etabin);
    h2PttrueCorr.push_back(static_cast<TH2F*>(h2->Clone(name)));
    
    sprintf(name,"h2PttrueCorrL2_eta%i",etabin);
    h2PttrueCorrL2.push_back(static_cast<TH2F*>(h2->Clone(name)));
    
    sprintf(name,"h2PttrueCorrL2L3_eta%i",etabin);
    h2PttrueCorrL2L3.push_back(static_cast<TH2F*>(h2->Clone(name)));
  } // End of loop over eta bins



  // Loop over data and fill response into
  // the corresponding 2D histogram
  for( data_it = mData->begin(); data_it != mData->end(); data_it++ ) {
    JetTruthEvent *jte = dynamic_cast<JetTruthEvent*>(*data_it);
    if( jte ) {
      if( jte->FlaggedBad() ) continue;   // Discard events flagged bad

      TJet * jet        = static_cast<TJet*>(jte->GetMess());
      double weight     = jte->GetWeight();
      double ptmeas     = jet->pt;
      double ptcorr     = jte->GetParametrizedMess();
      double ptcorrL2   = jet->L2cor * jet->pt;
      double ptcorrL2L3 = jet->L2L3cor * jet->pt;
      double eta        = jte->GetMess()->eta;
      double pttrue     = jte->GetTruth();
    
      // Find pttrue bin for response vs eta plot
      int ptbin = bins.IX(pttrue);
      if( 0 <= ptbin && ptbin < bins.NBinsX() ) {
	h2EtaUncorr.at(ptbin)->Fill(eta,ptmeas/pttrue,weight);
	h2EtaCorr.at(ptbin)->Fill(eta,ptcorr/pttrue,weight);
	h2EtaCorrL2.at(ptbin)->Fill(eta,ptcorrL2/pttrue,weight);
	h2EtaCorrL2L3.at(ptbin)->Fill(eta,ptcorrL2L3/pttrue,weight);
      }

      // Find pttrue bin for response vs eta plot
      int etabin = bins.IY(eta);
      if( 0 <= etabin && etabin < bins.NBinsY() ) {
	h2PttrueUncorr.at(etabin)->Fill(pttrue,ptmeas/pttrue,weight);
	h2PttrueCorr.at(etabin)->Fill(pttrue,ptcorr/pttrue,weight);
	h2PttrueCorrL2.at(etabin)->Fill(pttrue,ptcorrL2/pttrue,weight);
	h2PttrueCorrL2L3.at(etabin)->Fill(pttrue,ptcorrL2L3/pttrue,weight);
      }
    }
  } // End of loop over data


  // -- Mean response vs eta, pt ---------------------------------
  // Indices:
  //  0 - Mean
  //  1 - Sigma
  //  2 - Mean of Gauss
  //  3 - Width of Gauss
  // Suffix:
  //  Uncorr   - Uncorrected mean response
  //  Corr     - Mean response corrected by kalibri fit
  //  CorrL2   - Mean response corrected by JetMET L2 correction
  //  CorrL2L3 - Mean response corrected by JetMET L2L3 correction
  std::vector< std::vector<TH1F*> > hRespEtaUncorr(4,std::vector<TH1F*>(bins.NBinsX()));
  std::vector< std::vector<TH1F*> > hRespEtaCorr(4,std::vector<TH1F*>(bins.NBinsX()));
  std::vector< std::vector<TH1F*> > hRespEtaCorrL2(4,std::vector<TH1F*>(bins.NBinsX()));
  std::vector< std::vector<TH1F*> > hRespEtaCorrL2L3(4,std::vector<TH1F*>(bins.NBinsX()));

  std::vector< std::vector<TH1F*> > hRespPttrueUncorr(4,std::vector<TH1F*>(bins.NBinsY()));
  std::vector< std::vector<TH1F*> > hRespPttrueCorr(4,std::vector<TH1F*>(bins.NBinsY()));
  std::vector< std::vector<TH1F*> > hRespPttrueCorrL2(4,std::vector<TH1F*>(bins.NBinsY()));
  std::vector< std::vector<TH1F*> > hRespPttrueCorrL2L3(4,std::vector<TH1F*>(bins.NBinsY()));

  for(int ptbin = 0; ptbin < bins.NBinsX(); ptbin++) { // Loop over pt bins
    char title[50];
    sprintf(title,"%.1f < p^{true}_{T} < %.1f GeV",bins.XLow(bins.Bin(ptbin,0)),bins.XUp(bins.Bin(ptbin,0)));
    TH1F * hProjection[8];
    TH1F * gaussplots[4];
    TF1 * gaussfits[4];
    // Uncorrected response
    Fit2D(h2EtaUncorr.at(ptbin),hProjection,gaussplots,gaussfits);
    for(int i = 0; i < 8; i++) {
      if( i < 4 ) {
	hProjection[i]->SetTitle(title);
	hProjection[i]->SetMarkerStyle(20);
	hProjection[i]->SetMarkerColor(1);
	hProjection[i]->SetLineColor(1);
	hRespEtaUncorr.at(i).at(ptbin) = hProjection[i];
	objToBeWritten.push_back(hProjection[i]);
	delete gaussplots[i];
	delete gaussfits[i];
      } else {
	delete hProjection[i];
      }
    }
    
    // Kalibri corrected response
    Fit2D(h2EtaCorr.at(ptbin),hProjection,gaussplots,gaussfits);
    for(int i = 0; i < 8; i++) {
      if( i < 4 ) {
	hProjection[i]->SetTitle(title);
	hProjection[i]->SetMarkerStyle(21);
	hProjection[i]->SetMarkerColor(2);
	hProjection[i]->SetLineColor(2);
	hRespEtaCorr.at(i).at(ptbin) = hProjection[i];
	objToBeWritten.push_back(hProjection[i]);
	delete gaussplots[i];
	delete gaussfits[i];
      } else {
	delete hProjection[i];
      }
    }

    // L2 corrected response
    Fit2D(h2EtaCorrL2.at(ptbin),hProjection,gaussplots,gaussfits);
    for(int i = 0; i < 8; i++) {
      if( i < 4 ) {
	hProjection[i]->SetTitle(title);
	hProjection[i]->SetMarkerStyle(22);
	hProjection[i]->SetMarkerColor(3);
	hProjection[i]->SetLineColor(3);
	hRespEtaCorrL2.at(i).at(ptbin) = hProjection[i];
	objToBeWritten.push_back(hProjection[i]);
	delete gaussplots[i];
	delete gaussfits[i];
      } else {
	delete hProjection[i];
      }
    }

    // L2L3 corrected response
    Fit2D(h2EtaCorrL2L3.at(ptbin),hProjection,gaussplots,gaussfits);
    for(int i = 0; i < 8; i++) {
      if( i < 4 ) {
	hProjection[i]->SetTitle(title);
	hProjection[i]->SetMarkerStyle(23);
	hProjection[i]->SetMarkerColor(4);
	hProjection[i]->SetLineColor(4);
	hRespEtaCorrL2L3.at(i).at(ptbin) = hProjection[i];
	objToBeWritten.push_back(hProjection[i]);
	delete gaussplots[i];
	delete gaussfits[i];
      } else {
	delete hProjection[i];
      }
    }
  } // End of loop over pt bins

  //Delete 2D histograms response vs eta
  for(int i = 0; i < bins.NBinsX(); i++) {
    delete h2EtaUncorr.at(i);
    delete h2EtaCorr.at(i);
    delete h2EtaCorrL2.at(i);
    delete h2EtaCorrL2L3.at(i);
  }

  for(int etabin = 0; etabin < bins.NBinsY(); etabin++) { // Loop over eta bins
    char title[50];
    sprintf(title,"%.1f <  #eta < %.1f",bins.YLow(bins.Bin(0,etabin)),bins.YUp(bins.Bin(0,etabin)));
    TH1F * hProjection[8];
    TH1F * gaussplots[4];
    TF1 * gaussfits[4];
    // Uncorrected response
    Fit2D(h2PttrueUncorr.at(etabin),hProjection,gaussplots,gaussfits);
    for(int i = 0; i < 8; i++) {
      if( i < 4 ) {
	hProjection[i]->SetTitle(title);
	hProjection[i]->SetMarkerStyle(20);
	hProjection[i]->SetMarkerColor(1);
	hProjection[i]->SetLineColor(1);
	hRespPttrueUncorr.at(i).at(etabin) = hProjection[i];
	objToBeWritten.push_back(hProjection[i]);
	delete gaussplots[i];
	delete gaussfits[i];
      } else {
	delete hProjection[i];
      }
    }
    
    // Kalibri corrected response
    Fit2D(h2PttrueCorr.at(etabin),hProjection,gaussplots,gaussfits);
    for(int i = 0; i < 8; i++) {
      if( i < 4 ) {
	hProjection[i]->SetTitle(title);
	hProjection[i]->SetMarkerStyle(21);
	hProjection[i]->SetMarkerColor(2);
	hProjection[i]->SetLineColor(2);
	hRespPttrueCorr.at(i).at(etabin) = hProjection[i];
	objToBeWritten.push_back(hProjection[i]);
	delete gaussplots[i];
	delete gaussfits[i];
      } else {
	delete hProjection[i];
      }
    }

    // L2 corrected response
    Fit2D(h2PttrueCorrL2.at(etabin),hProjection,gaussplots,gaussfits);
    for(int i = 0; i < 8; i++) {
      if( i < 4 ) {
	hProjection[i]->SetTitle(title);
	hProjection[i]->SetMarkerStyle(22);
	hProjection[i]->SetMarkerColor(3);
	hProjection[i]->SetLineColor(3);
	hRespPttrueCorrL2.at(i).at(etabin) = hProjection[i];
	objToBeWritten.push_back(hProjection[i]);
	delete gaussplots[i];
	delete gaussfits[i];
      } else {
	delete hProjection[i];
      }
    }

    // L2L3 corrected response
    Fit2D(h2PttrueCorrL2L3.at(etabin),hProjection,gaussplots,gaussfits);
    for(int i = 0; i < 8; i++) {
      if( i < 4 ) {
	hProjection[i]->SetTitle(title);
	hProjection[i]->SetMarkerStyle(23);
	hProjection[i]->SetMarkerColor(4);
	hProjection[i]->SetLineColor(4);
	hRespPttrueCorrL2L3.at(i).at(etabin) = hProjection[i];
	objToBeWritten.push_back(hProjection[i]);
	delete gaussplots[i];
	delete gaussfits[i];
      } else {
	delete hProjection[i];
      }
    }
  } // End of loop over pt bins


  //Delete 2D histograms response vs pttrue
  for(int i = 0; i < bins.NBinsY(); i++) {
    delete h2PttrueUncorr.at(i);
    delete h2PttrueCorr.at(i);
    delete h2PttrueCorrL2.at(i);
    delete h2PttrueCorrL2L3.at(i);
  }


  // -- Draw histograms into ps file ----------------------------
  TLegend * leg = new TLegend(0.4,0.65,0.93,0.85);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextFont(42);
  leg->AddEntry(hRespEtaUncorr.at(0).at(0),"Uncorrected","P");
  leg->AddEntry(hRespEtaCorr.at(0).at(0),"Corrected Kalibri","P");
  //leg->AddEntry(hRespEtaCorrL2.at(0).at(0),"Corrected L2","P");
  leg->AddEntry(hRespEtaCorrL2L3.at(0).at(0),"Corrected L2L3","P");

  // Mean response vs eta per pttrue bin
  TPostScript       * ps = new TPostScript("controlplotsJetTruthEventResponse.ps",111);
  TCanvas     * const c1 = new TCanvas("c1","",500,500);
  c1->cd();

  TH1F * h = hRespEtaUncorr.at(0).at(0);
  TLine *leta = new TLine(h->GetXaxis()->GetBinLowEdge(1),1,
			  h->GetXaxis()->GetBinUpEdge(h->GetNbinsX()),1);
  leta->SetLineColor(1);
  leta->SetLineWidth(2);
  leta->SetLineStyle(2);

  for(int ptbin = 0; ptbin < bins.NBinsX(); ptbin++) {
    // Mean response vs eta
    hRespEtaUncorr.at(0).at(ptbin)->GetYaxis()->SetRangeUser(0,2);
    hRespEtaUncorr.at(0).at(ptbin)->GetYaxis()->SetTitle("< p^{jet}_{T} / p^{true}_{T} >");
    hRespEtaUncorr.at(0).at(ptbin)->Draw("PE1");
    leta->Draw("same");
    hRespEtaCorr.at(0).at(ptbin)->Draw("PE1same");
    //hRespEtaCorrL2.at(0).at(ptbin)->Draw("PE1same");
    hRespEtaCorrL2L3.at(0).at(ptbin)->Draw("PE1same");
    leg->Draw("same");
    c1->Draw();
    ps->NewPage();

    // Mean of Gauss fit on response vs eta
    hRespEtaUncorr.at(2).at(ptbin)->GetYaxis()->SetRangeUser(0,2);
    hRespEtaUncorr.at(2).at(ptbin)->GetYaxis()->SetTitle("Gauss fit < p^{jet}_{T} / p^{true}_{T} >");
    hRespEtaUncorr.at(2).at(ptbin)->Draw("PE1");
    leta->Draw("same");
    hRespEtaCorr.at(2).at(ptbin)->Draw("PE1same");
    //hRespEtaCorrL2.at(2).at(ptbin)->Draw("PE1same");
    hRespEtaCorrL2L3.at(2).at(ptbin)->Draw("PE1same");
    leg->Draw("same");
    c1->Draw();
    ps->NewPage();
  }
  for(int ptbin = 0; ptbin < bins.NBinsX(); ptbin++) {
    // Zoom: Mean response vs eta
    hRespEtaUncorr.at(0).at(ptbin)->GetYaxis()->SetRangeUser(0.8,1.2);
    hRespEtaUncorr.at(0).at(ptbin)->Draw("PE1");
    leta->Draw("same");
    hRespEtaCorr.at(0).at(ptbin)->Draw("PE1same");
    //hRespEtaCorrL2.at(0).at(ptbin)->Draw("PE1same");
    hRespEtaCorrL2L3.at(0).at(ptbin)->Draw("PE1same");
    leg->Draw("same");
    c1->Draw();
    ps->NewPage();
    // Zoom: Mean of Gauss fit on response vs eta
    hRespEtaUncorr.at(2).at(ptbin)->GetYaxis()->SetRangeUser(0.8,1.2);
    hRespEtaUncorr.at(2).at(ptbin)->GetYaxis()->SetTitle("Gauss fit < p^{jet}_{T} / p^{true}_{T} >");
    hRespEtaUncorr.at(2).at(ptbin)->Draw("PE1");
    leta->Draw("same");
    hRespEtaCorr.at(2).at(ptbin)->Draw("PE1same");
    //hRespEtaCorrL2.at(2).at(ptbin)->Draw("PE1same");
    hRespEtaCorrL2L3.at(2).at(ptbin)->Draw("PE1same");
    leg->Draw("same");
    c1->Draw();
    ps->NewPage();
  }

  // Mean response vs pttrue per eta bin
  h = hRespPttrueUncorr.at(0).at(0);
  TLine *lpttrue = new TLine(h->GetXaxis()->GetBinLowEdge(1),1,
			  h->GetXaxis()->GetBinUpEdge(h->GetNbinsX()),1);
  lpttrue->SetLineColor(1);
  lpttrue->SetLineWidth(2);
  lpttrue->SetLineStyle(2);
  for(int etabin = 0; etabin < bins.NBinsY(); etabin++) {
    // Mean response vs pttrue
    hRespPttrueUncorr.at(0).at(etabin)->GetYaxis()->SetRangeUser(0,2);
    hRespPttrueUncorr.at(0).at(etabin)->GetYaxis()->SetTitle("< p^{jet}_{T} / p^{true}_{T} >");
    hRespPttrueUncorr.at(0).at(etabin)->Draw("PE1");
    lpttrue->Draw("same");
    hRespPttrueCorr.at(0).at(etabin)->Draw("PE1same");
    //hRespPttrueCorrL2.at(0).at(etabin)->Draw("PE1same");
    hRespPttrueCorrL2L3.at(0).at(etabin)->Draw("PE1same");
    leg->Draw("same");
    c1->SetLogx(1);
    c1->Draw();
    ps->NewPage();
    // Mean response of Gauss fit vs pttrue
    hRespPttrueUncorr.at(2).at(etabin)->GetYaxis()->SetRangeUser(0,2);
    hRespPttrueUncorr.at(2).at(etabin)->GetYaxis()->SetTitle("Gauss fit < p^{jet}_{T} / p^{true}_{T} >");
    hRespPttrueUncorr.at(2).at(etabin)->Draw("PE1");
    lpttrue->Draw("same");
    hRespPttrueCorr.at(2).at(etabin)->Draw("PE1same");
    //hRespPttrueCorrL2.at(2).at(etabin)->Draw("PE1same");
    hRespPttrueCorrL2L3.at(2).at(etabin)->Draw("PE1same");
    leg->Draw("same");
    c1->Draw();
    ps->NewPage();
  }
  for(int etabin = 0; etabin < bins.NBinsY(); etabin++) {
    // Zoom: Mean response vs pttrue
    hRespPttrueUncorr.at(0).at(etabin)->GetYaxis()->SetRangeUser(0.8,1.2);
    hRespPttrueUncorr.at(0).at(etabin)->Draw("PE1");
    lpttrue->Draw("same");
    hRespPttrueCorr.at(0).at(etabin)->Draw("PE1same");
    //hRespPttrueCorrL2.at(0).at(etabin)->Draw("PE1same");
    hRespPttrueCorrL2L3.at(0).at(etabin)->Draw("PE1same");
    leg->Draw("same");
    c1->Draw();
    ps->NewPage();
    // Zoom: Mean response of Gauss fit vs pttrue
    hRespPttrueUncorr.at(2).at(etabin)->GetYaxis()->SetRangeUser(0.8,1.2);
    hRespPttrueUncorr.at(2).at(etabin)->Draw("PE1");
    lpttrue->Draw("same");
    hRespPttrueCorr.at(2).at(etabin)->Draw("PE1same");
    //hRespPttrueCorrL2.at(2).at(etabin)->Draw("PE1same");
    hRespPttrueCorrL2L3.at(2).at(etabin)->Draw("PE1same");
    leg->Draw("same");
    c1->Draw();
    ps->NewPage();
  }
  delete ps;
  delete leta;
  delete lpttrue;

  // Resolution vs eta, pttrue
  ps = new TPostScript("controlplotsJetTruthEventResolution.ps",111);
  c1->cd();
  c1->SetLogx(0);

  // Resolution vs eta per pttrue bin
  for(int ptbin = 0; ptbin < bins.NBinsX(); ptbin++) {
    // Sigma
    hRespEtaUncorr.at(1).at(ptbin)->GetYaxis()->SetRangeUser(0,0.4);
    hRespEtaUncorr.at(1).at(ptbin)->GetYaxis()->SetTitle("#sigma( p^{jet}_{T} / p^{true}_{T} )  /  < p^{jet}_{T} / p^{true}_{T} >");
    hRespEtaUncorr.at(1).at(ptbin)->Draw("PE1");
    hRespEtaCorr.at(1).at(ptbin)->Draw("PE1same");
    //hRespEtaCorrL2.at(1).at(ptbin)->Draw("PE1same");
    hRespEtaCorrL2L3.at(1).at(ptbin)->Draw("PE1same");
    leg->Draw("same");
    c1->Draw();
    ps->NewPage();
    // Width of Gauss fit
    hRespEtaUncorr.at(3).at(ptbin)->GetYaxis()->SetRangeUser(0,0.4);
    hRespEtaUncorr.at(3).at(ptbin)->GetYaxis()->SetTitle("Gauss fit #sigma( p^{jet}_{T} / p^{true}_{T} )  /  < p^{jet}_{T} / p^{true}_{T} >");
    hRespEtaUncorr.at(3).at(ptbin)->Draw("PE1");
    hRespEtaCorr.at(3).at(ptbin)->Draw("PE1same");
    //hRespEtaCorrL2.at(3).at(ptbin)->Draw("PE1same");
    hRespEtaCorrL2L3.at(3).at(ptbin)->Draw("PE1same");
    leg->Draw("same");
    c1->Draw();
    ps->NewPage();
  }
  // Resolution vs pttrue per eta bin
  for(int etabin = 0; etabin < bins.NBinsY(); etabin++) {
    // Sigma
    hRespPttrueUncorr.at(1).at(etabin)->GetYaxis()->SetRangeUser(0,0.4);
    hRespPttrueUncorr.at(1).at(etabin)->GetYaxis()->SetTitle("#sigma( p^{jet}_{T} / p^{true}_{T} )  /  < p^{jet}_{T} / p^{true}_{T} >");
    hRespPttrueUncorr.at(1).at(etabin)->Draw("PE1");
    hRespPttrueCorr.at(1).at(etabin)->Draw("PE1same");
    //hRespPttrueCorrL2.at(1).at(etabin)->Draw("PE1same");
    hRespPttrueCorrL2L3.at(1).at(etabin)->Draw("PE1same");
    leg->Draw("same");
    c1->SetLogx(1);
    c1->Draw();
    ps->NewPage();
    // Width of Gauss fit
    hRespPttrueUncorr.at(3).at(etabin)->GetYaxis()->SetRangeUser(0,0.4);
    hRespPttrueUncorr.at(3).at(etabin)->GetYaxis()->SetTitle("Gauss fit #sigma( p^{jet}_{T} / p^{true}_{T} )  /  < p^{jet}_{T} / p^{true}_{T} >");
    hRespPttrueUncorr.at(3).at(etabin)->Draw("PE1");
    hRespPttrueCorr.at(3).at(etabin)->Draw("PE1same");
    //hRespPttrueCorrL2.at(3).at(etabin)->Draw("PE1same");
    hRespPttrueCorrL2L3.at(3).at(etabin)->Draw("PE1same");
    leg->Draw("same");
    c1->Draw();
    ps->NewPage();
  }

  if( mOutputROOT ) WriteToRootFile(objToBeWritten,"JetTruthEventResponse");

  // Clean up
  ps->Close();
  std::vector<TObject*>::iterator it = objToBeWritten.begin();
  for(; it != objToBeWritten.end(); it++) {
    delete *it;
  }
  delete c1;
  delete ps;
  delete leg;
}



//!  \brief Fitted correction functions
//!
//!  The plots are written to the file
//!  "controlplotsCorrectionFunction.ps" and (if enabled)
//!  to the directory "CorrectionFunction" in the file
//!  "controlplots.root".
// -------------------------------------------------------------
void TControlPlots::MakeControlPlotsCorrectionFunction() {
  std::vector<TObject*>                objToBeWritten;
  std::vector<TData*>::const_iterator  data_it;

  std::vector<double> binEdgesPt = bag_of<double>(mConfig->read<std::string>("Control plots pt bin edges",""));


  // -- Global jet correctin function ----------
  // Create histogram
  TH1F * hGlobalJetCorr = new TH1F("hGlobalJetCorr",
				   "Global jet correction;p^{jet}_{T} (GeV);Correction factor",
				   500,binEdgesPt.front(),binEdgesPt.back());
  hGlobalJetCorr->SetLineWidth(2);
  objToBeWritten.push_back(hGlobalJetCorr);

  // Loop over pt bins and fill histo with correction factor
//   std::cout << ">> Global jet parameters:\n";
//   for(int i = 0; i < 4; i++) {
//     std::cout << i << ": " << mPar->GetGlobalJetParRef()[i] << std::endl;
//   }
  TMeasurement meas;
  for(int ptbin = 1; ptbin <= hGlobalJetCorr->GetNbinsX(); ptbin++) {
    meas.pt = hGlobalJetCorr->GetBinCenter(ptbin);
    double ptcorr = TParameters::global_jet_parametrization(&meas,mPar->GetGlobalJetParRef());
    hGlobalJetCorr->SetBinContent(ptbin,ptcorr/meas.pt);
    hGlobalJetCorr->GetYaxis()->SetRangeUser(1.0,3.0);
  }


  // -- Local jet correctin function -----------
  // Create one histogram per jet-bin
  std::vector<TH1F*> hLocalJetCorr(mPar->GetEtaGranularityJet());
  for(int jetbin = 0; jetbin < mPar->GetEtaGranularityJet(); jetbin++) {
    char name[50];
    sprintf(name,"hLocalJetCorr_bin%i",jetbin);
    char title[100];
    sprintf(title,"Jet correction - Bin %i;Correction factor",jetbin);
    TH1F * h = new TH1F(name,title,500,binEdgesPt.front(),binEdgesPt.back());
    h->SetLineWidth(2);

    // Loop over ptbins and fill histo with correction factor
    for(int ptbin = 1; ptbin <= h->GetNbinsX(); ptbin++) {
    meas.pt = h->GetBinCenter(ptbin);
    double ptcorr = TParameters::jet_parametrization(&meas,mPar->GetJetParRef(jetbin));
    h->SetBinContent(ptbin,ptcorr/meas.pt);
    }
    h->GetYaxis()->SetRangeUser(0.5,2.0);
    hLocalJetCorr.at(jetbin) = h;
    objToBeWritten.push_back(h);
  }
  

  // -- Draw histograms into ps file -----------
  TPostScript * const ps = new TPostScript("controlplotsCorrectionFunction.ps",111);
  TCanvas     * const c1 = new TCanvas("c1","",500,500);
  c1->cd();
  c1->SetLogx(1);

  // Global jet correction
  hGlobalJetCorr->Draw("L");
  c1->Draw();
  ps->NewPage();

  // Local jet correction
  for(unsigned int i = 0; i < hLocalJetCorr.size(); i++) {
    hLocalJetCorr.at(i)->Draw("L");
    c1->Draw();
    ps->NewPage();
  }


  if( mOutputROOT ) WriteToRootFile(objToBeWritten,"CorrectionFunction");


  // Clean up
  ps->Close();
  std::vector<TObject*>::iterator it = objToBeWritten.begin();
  for(; it != objToBeWritten.end(); it++) {
    delete *it;
  }
  delete c1;
  delete ps;
}





//---------------------------------------------------------------
//   Tower control plots
//---------------------------------------------------------------
void TControlPlots::MakeControlPlotsTowers()
{
  std::vector<TObject*> objToBeWritten;
  std::vector<TData*>::const_iterator data_it;

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

  TH2F * constants = new TH2F("hCalibConstants","Calibration constants vs. E_{T} and #eta-bin with f_{em} = OUF = 0;E_{T} [GeV];#eta-bin",100,0.5,100.,mPar->GetEtaGranularity(),1,mPar->GetEtaGranularity()); 
   
  for (int eta=0; eta<mPar->GetEtaGranularity();++eta) // Loop over eta bins
    {
      for (int phi=0; phi<mPar->GetPhiGranularity();++phi) // Loop over phi bins
	{
	  int i = mPar->GetBin(eta,phi);
	  double * val = mPar->GetTowerParRef(i);

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

	  data_it = mData->begin();
	  for (; data_it != mData->end();++data_it) // 1. loop over all fit events
	    {
	      //if one fit event is composed of multiple towers, than loop over all
	      TAbstractData* ad = dynamic_cast<TAbstractData*>(*data_it);
	      if(! ad) continue;
	      const std::vector<TAbstractData*>& data_ref = ad->GetRef();

	      // Get sum of corrected tower Et
	      Jet = 0.;
	      for(std::vector<TAbstractData*>::const_iterator it = data_ref.begin(); 
		  it != data_ref.end(); ++it) //1. loop over towers
		{
		  Jet += (*it)->GetParametrizedMess();
		}

	      // Find tower with max Et
	      maxTowerET = 0.;	      
	      thisIndexJet = 0;
	      for(std::vector<TAbstractData*>::const_iterator it = data_ref.begin(); 
		  it != data_ref.end(); ++it) //2. loop over towers
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
	  
	  for (data_it = mData->begin(); data_it != mData->end();++data_it) // 2. loop over all fit-events
	    {
	      //if one fit event is composed of multiple towers, than loop over all
	      mess = 0.;
	      error=0.;


	      JetCorr = (*data_it)->GetParametrizedMess(); // Corrected jet Pt
	      Jet = 0.;		// Sum of corrected tower Et
	      
	      TAbstractData* ad = dynamic_cast<TAbstractData*>(*data_it);
	      if(! ad) continue;
	      const std::vector<TAbstractData*>& data_ref = ad->GetRef();
	      for(std::vector<TAbstractData*>::const_iterator it =data_ref.begin(); 
		  it!=data_ref.end(); ++it) // 3. loop over towers
		{
		  Jet += (*it)->GetParametrizedMess();
		}
	      maxTowerET = 0.;
	      thisIndexJet = 0;
	      for(std::vector<TAbstractData*>::const_iterator it =data_ref.begin(); 
		  it!=data_ref.end(); ++it) // 4. loop over towers
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
		      khad[0]->Fill(     mhad, mPar->plot_parametrization(testmess,val) );
		      if ((*data_it)->GetType()==GammaJet)
			{
			  plot_had[2]->Fill( mhad, t/m );
			  khad[2]->Fill(     mhad, mPar->plot_parametrization(testmess,val) );
			}  
		      else if ((*data_it)->GetType()==TrackTower)
			{
			  khad[1]->Fill(     mhad, mPar->plot_parametrization(testmess,val) );
			  plot_had[1]->Fill( mhad, t/m );
			}  
		      else if ((*data_it)->GetType()==TrackCluster)
			{
			  plot_had[3]->Fill( mhad, t/m );
			  khad[3]->Fill(     mhad, mPar->plot_parametrization(testmess,val) );
			}  
		      else if ((*data_it)->GetType()==PtBalance)
			{
			  plot_had[4]->Fill( mhad, t/m );
			  khad[4]->Fill(     mhad, mPar->plot_parametrization(testmess,val) );
			}  
		    }

		  if (m!=0.0)
		    {
		      plot[0]->Fill( m, t/m );
	    
		      testmess->pt = m;
		      testmess->EMF = em[0]->GetBinContent(em[0]->GetXaxis()->FindBin(m));
		      testmess->HadF = had[0]->GetBinContent(had[0]->GetXaxis()->FindBin(m));
		      testmess->OutF = au[0]->GetBinContent(au[0]->GetXaxis()->FindBin(m));
		      k[0]->Fill( m, mPar->plot_parametrization(testmess,val) );

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
			  k[2]->Fill( m, mPar->plot_parametrization(testmess,val) );
			}  
		      else if ((*data_it)->GetType()==TrackTower)
			{
			  plot[1]->Fill( m, t/m );
			  plot2dtt->Fill( m, t/m );
			  testmess->EMF = em[1]->GetBinContent(em[1]->GetXaxis()->FindBin(m));
			  testmess->HadF = had[1]->GetBinContent(had[1]->GetXaxis()->FindBin(m));
			  testmess->OutF = au[1]->GetBinContent(au[1]->GetXaxis()->FindBin(m));
			  k[1]->Fill( m, mPar->plot_parametrization(testmess,val) );
			}  
		      else if ((*data_it)->GetType()==TrackCluster) 
			{
			  plot[3]->Fill( m, t/m );
			  plot2dtc->Fill( m, t/m );
			  testmess->EMF = em[3]->GetBinContent(em[3]->GetXaxis()->FindBin(m));
			  testmess->HadF = had[3]->GetBinContent(had[3]->GetXaxis()->FindBin(m));
			  testmess->OutF = au[3]->GetBinContent(au[3]->GetXaxis()->FindBin(m));
			  k[3]->Fill( m, mPar->plot_parametrization(testmess,val) );
			}  
		      else if ((*data_it)->GetType()==PtBalance) 
			{
			  plot[4]->Fill( m, t/m );
			  plot2djj->Fill( m, t/m );
			  testmess->EMF = em[3]->GetBinContent(em[3]->GetXaxis()->FindBin(m));
			  testmess->HadF = had[3]->GetBinContent(had[3]->GetXaxis()->FindBin(m));
			  testmess->OutF = au[3]->GetBinContent(au[3]->GetXaxis()->FindBin(m));
			  k[4]->Fill( m, mPar->plot_parametrization(testmess,val) );
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
				      mPar->plot_parametrization(testmess,val) );
	      constants->SetBinContent(constants->GetXaxis()->FindBin(b),constants->GetYaxis()->FindBin(eta+1),
				       mPar->plot_parametrization(testmess,val));
	      testmess->EMF = (double)b*0.2; // -> EMFrac = 0.25
	      testmess->HadF = (double)b*0.8;
	      kEfrac02->SetBinContent(kEfrac02->GetXaxis()->FindBin(b), 
				      mPar->plot_parametrization(testmess,val) );       
	      testmess->EMF = (double)b*0.5; // -> EMFrac = 1.0
	      testmess->HadF = (double)b*0.5;
	      kEfrac05->SetBinContent(kEfrac05->GetXaxis()->FindBin(b), 
				      mPar->plot_parametrization(testmess,val) );       
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

  if( mOutputROOT ) WriteToRootFile(objToBeWritten, "Towers");
    
  ps->Close();

  delete latex;
  delete constants;
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
  std::vector<std::string> tmpPlottedQuant = bag_of_string(mConfig->read<std::string>("gamma jet plotted quantities",""));
  std::set<std::string> plottedQuant;
  //do we really have to copy the vector and not just use mPlottedQuant in the first place???
  std::vector<std::string>::const_iterator it = tmpPlottedQuant.begin();
  for(; it < tmpPlottedQuant.end(); it++) {
    plottedQuant.insert(*it);
  }

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
  TH1F* Chi2Plot[4];
  Chi2Plot[0] = new TH1F("chi2","#chi2;#chi2",200,0,100);
  Chi2Plot[1] = (TH1F*)Chi2Plot[0]->Clone("Chi2TrackOnly");
  Chi2Plot[2] = (TH1F*)Chi2Plot[0]->Clone("Chi2CaloOnly");
  Chi2Plot[3] = (TH1F*)Chi2Plot[0]->Clone("Chi2TrackOnlyNoWeight");
  TH2F* Chi2Pt[2];
  Chi2Pt[0] = new TH2F("chi2Pt","#chi2Pt;Pt;#chi2",100,0,200,100,0,30);
  Chi2Pt[0]->SetMarkerStyle(22);
  Chi2Pt[1] = (TH2F*)Chi2Pt[0]->Clone();
  TH2F* Chi2Eta =  new TH2F("chi2Eta","#chi2Eta;Eta;#chi2",500,-5,5,100,0,30);
  TH2F* Chi2NoTrack =  new TH2F("chi2NoTrack","#chi2NoTrack;NoTrack;#chi2",25,0,25,100,0,30);
  TH2F* Chi2Error =  new TH2F("chi2Error","#chi2Error;Error;#chi2",50,0,50,100,0,30);
  TH2F* Diff2Pt[2];
  Diff2Pt[0] = new TH2F("diff2PtTrackOnly","diff2Pt;Pt;diff2",100,0,1000,100,0,100);
  Diff2Pt[0]->SetMarkerStyle(22);
  Diff2Pt[1] = (TH2F*)Diff2Pt[0]->Clone("diff2PtCaloOnly");
  TH2F* TrackCalo[3];
  TrackCalo[0] = new TH2F("TrackPtOverCaloPt",";CaloPt; ",100,0,1000,75,0,1.5);
  TrackCalo[1] = (TH2F*)TrackCalo[0]->Clone("TrackResponstOverCaloPt");
  TrackCalo[2] = (TH2F*)TrackCalo[0]->Clone("Part of Energy from TrackPt > 50Gev");

  TH2F *RelResPt[8];
  RelResPt[0] = new TH2F("hRelResPt0","Z-Jet;p_{T}(gen) [GeV];Jet/genJet",20,0,150,50,0.5,1.5);
  RelResPt[0]->Sumw2();  
  RelResPt[1] = (TH2F*)RelResPt[0]->Clone("hRelResPt1"); 
  RelResPt[2] = (TH2F*)RelResPt[0]->Clone("hRelResPt2");
  RelResPt[3] = (TH2F*)RelResPt[0]->Clone("hRelResPt3");
  RelResPt[4] = new TH2F("hRelResPt4","Z-Jet;p_{T}(gen) [GeV];Jet/genJet",20,0,500,50,0.5,1.5);  
  RelResPt[5] = (TH2F*)RelResPt[4]->Clone("hRelResPt5"); 
  RelResPt[6] = (TH2F*)RelResPt[4]->Clone("hRelResPt6");
  RelResPt[7] = (TH2F*)RelResPt[4]->Clone("hRelResPt7");

  TH2F *RelResEMF[4];
  RelResEMF[0] = new TH2F("hRelResEMF0","Z-Jet;EMF;Jet/genJet",20,0,1,50,0.5,1.5);  
  RelResEMF[1] = (TH2F*)RelResEMF[0]->Clone("hRelResEMF1"); 
  RelResEMF[2] = (TH2F*)RelResEMF[0]->Clone("hRelResEMF2");
  RelResEMF[3] = (TH2F*)RelResEMF[0]->Clone("hRelResEMF3");

  TH2F *RelResEta[4];
  RelResEta[0] = new TH2F("hRelResEta0","Z-Jet;Eta;Jet/genJet",20,-2.5,2.5,50,0.5,1.5);  
  RelResEta[1] = (TH2F*)RelResEta[0]->Clone("hRelResEta1"); 
  RelResEta[2] = (TH2F*)RelResEta[0]->Clone("hRelResEta2"); 
  RelResEta[3] = (TH2F*)RelResEta[0]->Clone("hRelResEta3");

  TH2F *RelResTrackMult[4];
  RelResTrackMult[0] = new TH2F("hRelResTrackMult0","Z-Jet;TrackMult;Jet/genJet",20,0,50,50,0.5,1.5);  
  RelResTrackMult[1] = (TH2F*)RelResTrackMult[0]->Clone("hRelResTrackMult1"); 
  RelResTrackMult[2] = (TH2F*)RelResTrackMult[0]->Clone("hRelResTrackMult2");
  RelResTrackMult[3] = (TH2F*)RelResTrackMult[0]->Clone("hRelResTrackMult3");


  // The following 2D histograms contain different Pt ratios:
  //   ptjet/etgamma
  //   ptjetcorr/etgamma 
  //   ptjet/ptjetcorr
  // versus different quantities (controlled via argument 'plottedQuant'):
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
  if( plottedQuant.count("eta") > 0 )
    {
      /*
      double etaBinEdge[83];
      for(int e = 0; e < 82; e++)
	{
	  int etaBin = 0;
	  if( e < 41 ) etaBin = -41 + e;
	  else etaBin = -40 + e;
	  etaBinEdge[e] = mPar->EtaLowerEdge(etaBin);
	}
      etaBinEdge[82] = mPar->EtaUpperEdge(41);
      heta[0] = new TH2F("heta0","#gamma-jet;#eta",82,etaBinEdge,100,0,4);
      */
      heta[0] = new TH2F("heta0","#gamma-jet;#eta",20,-5,5,100,0,4);
      objToBeWritten.push_back(heta[0]);
      for(int i = 1 ; i < 12 ; ++i)
	{
	  sprintf(name,"heta%i",i);
	  heta[i] = (TH2F*)heta[0]->Clone(name);
	  heta[i]->Sumw2();
	  if( i > 2 )
	    {
	      int etReg = (i-3)/3;
	      sprintf(title,"#gamma-jet,  %i < E_{T}^{#gamma} < %i GeV",etLimit[etReg],etLimit[etReg+1]);
	      heta[i]->SetTitle(title);
	    }
	  objToBeWritten.push_back(heta[i]);
	}
    }

  TH2F* hpt_uncorr[12];
  if( plottedQuant.count("uncorrected jet pt") > 0 )
    {
      hpt_uncorr[0] = new TH2F("hpt_uncorr0","#gamma-jet;p^{jet}_{T} [GeV]",80,0,400,100,0,4);
      objToBeWritten.push_back(hpt_uncorr[0]);
      for(int i = 1 ; i < 12 ; ++i)
	{
	  sprintf(name,"hpt_uncorr%i",i);
	  hpt_uncorr[i] = (TH2F*) hpt_uncorr[0]->Clone(name);
	  hpt_uncorr[i]->Sumw2();

	  if( i > 2 )
	    {
	      int etReg = (i-3)/3;
	      if( etReg == 0 ) sprintf(title,"#gamma-jet,  |#eta| < 1.4");
	      if( etReg == 1 ) sprintf(title,"#gamma-jet,  1.4 < |#eta| < 3.0");
	      if( etReg == 2 ) sprintf(title,"#gamma-jet,  3.0 < |#eta|");
	      hpt_uncorr[i]->SetTitle(title);
	    }
	  objToBeWritten.push_back(heta[i]);
	}
    }

  TH2F* hpt[3];
  if( plottedQuant.count("true jet pt") > 0 )
    {
      hpt[0] = new TH2F("hpt0","#gamma-jet;E^{#gamma}_{T} [GeV]",80,0,400,100,0,4);
      hpt[1] = (TH2F*)hpt[0]->Clone("hpt1");
      hpt[2] = (TH2F*)hpt[0]->Clone("hpt2");
      for(int i = 0; i < 3; i++)
	{
	  hpt[i]->Sumw2();
	  objToBeWritten.push_back(hpt[i]);
	}
    }

  TH2F* hetaTrack[3];
  if( plottedQuant.count("eta") > 0 )
    {
      hetaTrack[0] = new TH2F("hetaTrack0","#gamma-jet Track only;#eta",20,-5,5,100,0,4);
      hetaTrack[1] = (TH2F*)hetaTrack[0]->Clone("hetaTrack1");
      hetaTrack[2] = (TH2F*)hetaTrack[0]->Clone("hetaTrack2");
      for(int i = 0; i < 3; i++)
	{
	  hetaTrack[i]->Sumw2();
	  objToBeWritten.push_back(hetaTrack[i]);
	}
    }

  TH2F* hptTrack[3];
  if( plottedQuant.count("true jet pt") > 0 )
    {
      hptTrack[0] = new TH2F("hptTrack0","#gamma-jet Track only;E^{#gamma}_{T} [GeV]",80,0,400,100,0,4);
      hptTrack[1] = (TH2F*)hptTrack[0]->Clone("hptTrack1");
      hptTrack[2] = (TH2F*)hptTrack[0]->Clone("hptTrack2");
      for(int i = 0; i < 3; i++)
	{
	  hptTrack[i]->Sumw2();
	  objToBeWritten.push_back(hptTrack[i]);
	}
    }

  TH2F* henergy[3];
  if( plottedQuant.count("uncorrected jet energy") > 0 )
    {
      henergy[0] = new TH2F("henergy0","#gamma-jet;E^{jet} [GeV]",120,0,600,100,0,4);
      henergy[1] = (TH2F*)henergy[0]->Clone("henergy1");
      henergy[2] = (TH2F*)henergy[0]->Clone("henergy2");
      for(int i = 0; i < 3; i++)
	{
	  henergy[i]->Sumw2();
	  objToBeWritten.push_back(henergy[i]);
	}
    }

  TH2F* hemf[3];
  if( plottedQuant.count("emf") > 0 )
    {
      hemf[0] = new TH2F("hemf0","#gamma-jet;electromagnetic fraction f_{em}",50,0,1,100,0,4);
      hemf[1] = (TH2F*)hemf[0]->Clone("hemf1");
      hemf[2] = (TH2F*)hemf[0]->Clone("hemf2");
      for(int i = 0; i < 3; i++)
	{
	  hemf[i]->Sumw2();
	  objToBeWritten.push_back(hemf[i]);
	}
    }

  double bins[101];
  for(int i = 0; i < 101 ; ++i)
    {
      bins[i] = pow(10,(i+32)/40.0);
    }
  TH2F* hptlog[3];
  if( plottedQuant.count("log true jet pt") > 0 )
    {
      hptlog[0] = new TH2F("hptlog0","#gamma-jet;E^{#gamma}_{T} [GeV]",100,bins,100,0,4);
      hptlog[1] = (TH2F*) hptlog[0]->Clone("hptlog1");
      hptlog[2] = (TH2F*) hptlog[0]->Clone("hptlog2");
      for(int i = 0; i < 3; i++)
	{
	  hptlog[i]->Sumw2();
	  objToBeWritten.push_back(hptlog[i]);
	}
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

  int track=0;
  int noTrack=0;

  // Fill histos

  //loop over all fit-events
  for ( std::vector<TData*>::const_iterator i = mData->begin() ; i != mData->end() ; ++i )
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
      double SumTrack=0;
      double SumTrackOverFifty=0;
      double SumTrackResp=0;
      int NoTracks = 0; //number of Tracks in Jet
      // first loop over tower
      TAbstractData* ad = dynamic_cast<TAbstractData*>(jg);

      if(ad) {
	if(ad->GetTrackuse()) track++;
	else noTrack++;

	const std::vector<TAbstractData*>& data_refT = ad->GetRefTrack();
	for(std::vector<TAbstractData*>::const_iterator it = data_refT.begin();it != data_refT.end(); ++it)
	  {
	    TTrack * m = (TTrack*)(*it)->GetMess();
	    double  pm = (*it)->GetParametrizedMess();
	    if(m->DRout<0.5){
	      SumTrack += m->pt;
	      if( m->pt > 50)
		SumTrackOverFifty += m->pt;
	      SumTrackResp += pm;
	      NoTracks++;
	    }
	  }
	TrackCalo[0]->Fill(etjet,SumTrack/etjet);
	TrackCalo[1]->Fill(etjet,SumTrackResp/etjet);
	TrackCalo[2]->Fill(etjet,SumTrackOverFifty/SumTrack);

	const std::vector<TAbstractData*>& data_ref = ad->GetRef();
	for(std::vector<TAbstractData*>::const_iterator it = data_ref.begin();it != data_ref.end(); ++it)
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
	for(std::vector<TAbstractData*>::const_iterator it = data_ref.begin();it != data_ref.end(); ++it)
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
      }
      Chi2Plot[0]->Fill(jg->chi2());
      //if(ad->GetTrackuse()) {
      if(true) {
	Chi2Plot[1]->Fill(jg->chi2());
	Chi2Plot[3]->Fill(jg->chi2()/jg->GetWeight());
	Chi2Pt[0]->Fill(etjet,jg->chi2());
	Chi2Eta->Fill(etajet,jg->chi2());
	Chi2NoTrack->Fill(NoTracks,jg->chi2());
	if(ad) Chi2Error->Fill(ad->GetParametrizedErr(),jg->chi2());
	Diff2Pt[0]->Fill(etjet,std::abs(etjetcor - jg->GetTruth()));

	double em1 = jg->GetMess()->EMF;
	double had1 = jg->GetMess()->HadF+jg->GetMess()->OutF;
	TJet* jet = (TJet*)(jg->GetMess());
	double genJet = jet->genPt;
	RelResPt[0]->Fill(genJet,(jg->GetParametrizedMess() )/genJet,jg->GetWeight());
	RelResPt[1]->Fill(genJet,(jet->L2L3cor * jg->GetMess()->pt )/genJet,jg->GetWeight());
	RelResPt[2]->Fill(genJet,( jet->L2L3JPTcor * jet->ZSPcor * jet->JPTcor * jg->GetMess()->pt )/genJet,jg->GetWeight());
	RelResPt[3]->Fill(genJet,(jg->GetMess()->pt )/genJet,jg->GetWeight());
	RelResPt[4]->Fill(genJet,(jg->GetParametrizedMess() )/genJet,jg->GetWeight());
	RelResPt[5]->Fill(genJet,(jet->L2L3cor * jg->GetMess()->pt )/genJet,jg->GetWeight());
	RelResPt[6]->Fill(genJet,(jet->L2L3JPTcor * jet->ZSPcor * jet->JPTcor * jg->GetMess()->pt )/genJet,jg->GetWeight());
	RelResPt[7]->Fill(genJet,(jg->GetMess()->pt )/genJet,jg->GetWeight());
	RelResEMF[0]->Fill(em1/(em1+had1),(jg->GetParametrizedMess() )/genJet,jg->GetWeight());
	RelResEMF[1]->Fill(em1/(em1+had1),(jet->L2L3cor * jg->GetMess()->pt )/genJet,jg->GetWeight());
	RelResEMF[2]->Fill(em1/(em1+had1),(jet->L2L3JPTcor * jet->ZSPcor * jet->JPTcor * jg->GetMess()->pt )/genJet,jg->GetWeight());
	RelResEMF[3]->Fill(em1/(em1+had1),(jg->GetMess()->pt )/genJet,jg->GetWeight());
	RelResEta[0]->Fill(jg->GetMess()->eta,(jg->GetParametrizedMess() )/genJet,jg->GetWeight());
	RelResEta[1]->Fill(jg->GetMess()->eta,(jet->L2L3cor * jg->GetMess()->pt )/genJet,jg->GetWeight());
	RelResEta[2]->Fill(jg->GetMess()->eta,(jet->L2L3JPTcor * jet->ZSPcor * jet->JPTcor * jg->GetMess()->pt )/genJet,jg->GetWeight());
	RelResEta[3]->Fill(jg->GetMess()->eta,(jg->GetMess()->pt )/genJet,jg->GetWeight());
	
	if(ad) {
	  const std::vector<TAbstractData*>& data_refT = ad->GetRefTrack();
	  int TrackMult =  data_refT.size();
	  RelResTrackMult[0]->Fill(TrackMult,(jg->GetParametrizedMess() )/genJet,jg->GetWeight());
	  RelResTrackMult[1]->Fill(TrackMult,(jet->L2L3cor * jg->GetMess()->pt )/genJet,jg->GetWeight());
	  RelResTrackMult[2]->Fill(TrackMult,(jet->L2L3JPTcor * jet->ZSPcor * jet->JPTcor * jg->GetMess()->pt )/genJet,jg->GetWeight());
	  RelResTrackMult[3]->Fill(TrackMult,(jg->GetMess()->pt )/genJet,jg->GetWeight());
	}
      }
      else {
	Chi2Plot[2]->Fill(jg->chi2());
	Chi2Pt[1]->Fill(etjet,jg->chi2());
	Diff2Pt[1]->Fill(etjet,std::abs(etjetcor - jg->GetTruth()));
      }
      towerinjet[0]->Fill(noTower);
      if (jg->GetTruth() > 10 && jg->GetTruth() < 35)
	{
	  towerinjet[1]->Fill(noTower);
	}
      else if (jg->GetTruth() > 35 && jg->GetTruth() < 90)
	{
	  towerinjet[2]->Fill(noTower);
	}
      else if (jg->GetTruth() > 90 && jg->GetTruth() < 300)
	{
	  towerinjet[3]->Fill(noTower);
	} 

      if( plottedQuant.count("eta") > 0 )
	{
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
	  if(ad && ad->GetTrackuse()) 
	    {
	      hetaTrack[0]->Fill(etajet,etjet/ jg->GetTruth(),jg->GetWeight());
	      hetaTrack[1]->Fill(etajet,etjetcor/ jg->GetTruth(),jg->GetWeight());
	      hetaTrack[2]->Fill(etajet,etjet/etjetcor,jg->GetWeight());
	    }
	}

      if( plottedQuant.count("true jet pt") > 0 )
	{
	  hpt[0]->Fill(jg->GetTruth(),etjet/jg->GetTruth(),jg->GetWeight());
	  hpt[1]->Fill(jg->GetTruth(),etjetcor/jg->GetTruth(),jg->GetWeight());
	  hpt[2]->Fill(jg->GetTruth(),etjet/etjetcor,jg->GetWeight());
	  if(ad && ad->GetTrackuse()) 
	    {   
	      hptTrack[0]->Fill(jg->GetTruth(),etjet/ jg->GetTruth(),jg->GetWeight());
	      hptTrack[1]->Fill(jg->GetTruth(),etjetcor/jg->GetTruth(),jg->GetWeight());
	      hptTrack[2]->Fill(jg->GetTruth(),etjet/etjetcor,jg->GetWeight());
	    }
	}

      if( plottedQuant.count("uncorrected jet energy") > 0 )
	{
	  henergy[0]->Fill(energyjet,etjet/ jg->GetTruth(),jg->GetWeight());
	  henergy[1]->Fill(energyjet,etjetcor/jg->GetTruth(),jg->GetWeight());
	  henergy[2]->Fill(energyjet,etjet/etjetcor,jg->GetWeight());    
	}

      if( plottedQuant.count("log true jet pt") > 0 )
      {
	hptlog[0]->Fill(jg->GetTruth(),etjet/ jg->GetTruth(),jg->GetWeight());
	hptlog[1]->Fill(jg->GetTruth(),etjetcor/jg->GetTruth(),jg->GetWeight());
	hptlog[2]->Fill(jg->GetTruth(),etjet/etjetcor,jg->GetWeight());
      }

      if( plottedQuant.count("uncorrected jet pt") > 0 )
	{
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
	}

      //em fraction plots     
      if( plottedQuant.count("emf") > 0 )
	{
	  double em = jg->GetMess()->EMF;
	  double had = jg->GetMess()->HadF + jg->GetMess()->OutF;
	  /*
	    for(std::vector<TData*>::const_iterator t = jg->GetRef().begin(); t != jg->GetRef().end(); ++t)
	    {
	    TData* tt = *t;
	    em  += tt->GetMess()->EMF;
	    had += tt->GetMess()->HadF;
	    had += tt->GetMess()->OutF;
	    }
	  */
	  hemf[0]->Fill(em/(em+had),etjet/jg->GetTruth(),jg->GetWeight());
	  hemf[1]->Fill(em/(em+had),etjetcor/jg->GetTruth(),jg->GetWeight());
	  hemf[2]->Fill(em/(em+had),etjet/etjetcor,jg->GetWeight());
	  hptGamma[2]->Fill(em/(em+had));
	  hptGammaW[2]->Fill(em/(em+had),jg->GetWeight());
	  hptGamma2D->Fill(jg->GetTruth(),em/(em+had));
	  hptGamma2DW->Fill(jg->GetTruth(),em/(em+had),jg->GetWeight());
	}

      hptGamma[0]->Fill(jg->GetTruth());
      hptGammaW[0]->Fill(jg->GetTruth(),jg->GetWeight());
      hptGamma[1]->Fill(etajet);
      hptGammaW[1]->Fill(etajet,jg->GetWeight());
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
  c1->SetLogy(1);
  for(int i=3;i<4;++i)
    {
      Chi2Plot[i]->Draw();
      c1->Draw();
      ps->NewPage();
    }
  /*
  Chi2Pt[0]->SetMarkerColor(2);
  Chi2Pt[0]->Draw();
  Chi2Pt[1]->Draw("Same");
  c1->Draw();
  ps->NewPage();
  Chi2Eta->Draw();
  c1->Draw();
  ps->NewPage();
  Chi2NoTrack->Draw();
  c1->Draw();
  ps->NewPage();
  Chi2Error->Draw();
  c1->Draw();
  ps->NewPage();
  Diff2Pt[0]->SetMarkerColor(2);
  Diff2Pt[0]->Draw();
  c1->Draw();
  ps->NewPage();
  Diff2Pt[1]->Draw();
  c1->Draw();
  ps->NewPage();
  */
  c1->SetLogy(0);

  //Basic Track Plots
  for(int i=0;i<3;++i)
    {
      TrackCalo[i]->Draw("box");
      c1->Draw();
      ps->NewPage();
    }

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
  leg->AddEntry(heta[0],mPtRatioName[0],"p");
  leg->AddEntry(heta[2],mPtRatioName[2],"p");
  leg->AddEntry(heta[1],mPtRatioName[1],"p");

  if( plottedQuant.count("eta") > 0 )
    {
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
	      hists_eta[a+i][4]->SetMinimum(0.4);
	      hists_eta[a+i][4]->SetMaximum(1.6);
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
		  gp_eta[i+a][b]->SetXTitle(mPtRatioName[a]);

		  // Set style and line color according to mPtRatioName
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
		  h->SetTitle("#gamma-jet,  " + mPtRatioName[a-3] + " " + mControlQuantityName[b]);
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
    }


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

  if( plottedQuant.count("uncorrected jet pt") > 0 )
    {
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
	      hists_ptuncorr[a+i][4]->SetMinimum(0.4);
	      hists_ptuncorr[a+i][4]->SetMaximum(1.6);
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
		  gp_ptuncorr[i+a][b]->SetXTitle(mPtRatioName[a]);

		  // Set style and line color according to mPtRatioName
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
    }



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

  if( plottedQuant.count("true jet pt") > 0 )
    {
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
	  hists_pttrue[a][4]->SetMinimum(0.4);
	  hists_pttrue[a][4]->SetMaximum(1.6);
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
	      gp_pttrue[a][b]->SetXTitle(mPtRatioName[a]);

	      // Set style and line color according to mPtRatioName
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
    }

  // Control quantities from the Pt ratio vs Et gamma.
  // TracksOnly!!!!!!!!!
  // First dimension 12:
  //   index i = 0,1,2 indicates plotted pt-ratio
  //     0: ptjet/etgamma
  //     1: ptjetcorr/etgamma
  //     2: ptjet/ptjetcorr
  // Second dimension 8 is the above specified control
  // quantity
  TH1F* hists_ptTrack[3][8];

  // Projected distributions for 4 example eta-bins
  // and the Gauss fits
  TH1F* gp_ptTrack[3][4];  
  TF1* gf_ptTrack[3][4];

  if( plottedQuant.count("true jet pt") > 0 )
    {
      for(int a = 0; a < 3; a++) // Loop over Pt ratios
	{
	  hptTrack[a]->SetMinimum(0.5);
	  hptTrack[a]->SetMaximum(1.2);
	  hptTrack[a]->SetMarkerStyle(markerStyle[a]);
	  hptTrack[a]->SetMarkerColor(markerColor[a]);
	  hptTrack[a]->SetLineColor(markerColor[a]);
	}

      // Do projections and determine control quantities
      Fit2D(hptTrack[0],hists_ptTrack[0],gp_ptTrack[0], gf_ptTrack[0]);
      Fit2D(hptTrack[1],hists_ptTrack[1],gp_ptTrack[1], gf_ptTrack[1]);
      Fit2D(hptTrack[2],hists_ptTrack[2],gp_ptTrack[2], gf_ptTrack[2]); 
      for(int a = 0; a<3;++a)
	{
	  for(int b = 0 ; b < 4 ; ++b)
	    {
	      hists_ptTrack[a][b]->SetMinimum(0.2);
	      hists_ptTrack[a][b]->SetMaximum(1.8);
	      ++b;
	      hists_ptTrack[a][b]->SetMinimum(0.0);
	      hists_ptTrack[a][b]->SetMaximum(0.5);
	    }
	  hists_ptTrack[a][4]->SetMinimum(0.4);
	  hists_ptTrack[a][4]->SetMaximum(1.6);
	}

      // Draw gaussplots for example pt bins
      // on multi-canvas
      for(int a = 0; a < 3; a++) // Loop over ptratios
	{
	  for(int b = 0; b < 3; b++) // Loop over example pt bins
	    {
	      // Find pt bin of gaussplot
	      int bin = 0;
	      if(  b == 0  )  bin = int(hptTrack[a]->GetNbinsX()/6);
	      else if(  b == 1  )  bin = int(hptTrack[a]->GetNbinsX()/3);
	      else if(  b == 2  )  bin = int(hptTrack[a]->GetNbinsX()/2);
	      float min = hptTrack[a]->GetXaxis()->GetBinLowEdge(bin);
	      float max = min + hptTrack[a]->GetXaxis()->GetBinWidth(bin);

	      // Set title according to pt bin
	      sprintf(title,"#gamma-jet Track only, %.1f < E^{#gamma}_{T} < %.1f GeV",min,max);
	      gp_ptTrack[a][b]->SetTitle(title);
	      gp_ptTrack[a][b]->SetXTitle(mPtRatioName[a]);

	      // Set style and line color according to mPtRatioName
	      gp_ptTrack[a][b]->SetMarkerStyle(markerStyle[a]);
	      gp_ptTrack[a][b]->SetMarkerColor(markerColor[a]);
	      gp_ptTrack[a][b]->SetLineColor(markerColor[a]);
	      gf_ptTrack[a][b]->SetLineColor(markerColor[a]);

	      // Plot gaussplots
	      c2->cd(1+b);
	      gp_ptTrack[a][b]->Draw();
	      gf_ptTrack[a][b]->Draw("same");

	      objToBeWritten.push_back(gp_ptTrack[a][b]);
	      objToBeWritten.push_back(gf_ptTrack[a][b]);
	    } // End of loop over example ptbins
	  c2->Draw();
	  ps->NewPage();
	} // End of loop over ptratios

      c1->cd();
      for(int j = 0 ; j < 8 ; ++j) // Loop over control quantities
	{
	  hists_ptTrack[0][j]->Draw("P");
	  hists_ptTrack[0][j]->SetStats(0);
	  hists_ptTrack[1][j]->Draw("P SAME");
	  hists_ptTrack[2][j]->Draw("P SAME");
	  leg->Draw("SAME");
	  c1->SetGrid();
	  c1->Draw();   
	  ps->NewPage(); 

	  objToBeWritten.push_back(hists_ptTrack[0][j]);
	  objToBeWritten.push_back(hists_ptTrack[1][j]);
	  objToBeWritten.push_back(hists_ptTrack[2][j]);
	}
    }


  // Control quantities from the Pt ratio vs Eta.
  // TracksOnly!!!!!!!!!
  // First dimension 12:
  //   index i = 0,1,2 indicates plotted pt-ratio
  //     0: ptjet/etgamma
  //     1: ptjetcorr/etgamma
  //     2: ptjet/ptjetcorr
  // Second dimension 8 is the above specified control
  // quantity
  TH1F* hists_etaTrack[3][8];

  // Projected distributions for 4 example eta-bins
  // and the Gauss fits
  TH1F* gp_etaTrack[3][4];  
  TF1* gf_etaTrack[3][4];

  if( plottedQuant.count("eta") > 0 )
    {
      for(int a = 0; a < 3; a++) // Loop over Pt ratios
	{
	  hetaTrack[a]->SetMinimum(0.5);
	  hetaTrack[a]->SetMaximum(1.2);
	  hetaTrack[a]->SetMarkerStyle(markerStyle[a]);
	  hetaTrack[a]->SetMarkerColor(markerColor[a]);
	  hetaTrack[a]->SetLineColor(markerColor[a]);
	}

      // Do projections and determine control quantities
      Fit2D(hetaTrack[0],hists_etaTrack[0],gp_etaTrack[0], gf_etaTrack[0]);
      Fit2D(hetaTrack[1],hists_etaTrack[1],gp_etaTrack[1], gf_etaTrack[1]);
      Fit2D(hetaTrack[2],hists_etaTrack[2],gp_etaTrack[2], gf_etaTrack[2]); 
      for(int a = 0; a<3;++a)
	{
	  for(int b = 0 ; b < 4 ; ++b)
	    {
	      hists_etaTrack[a][b]->SetMinimum(0.2);
	      hists_etaTrack[a][b]->SetMaximum(1.8);
	      ++b;
	      hists_etaTrack[a][b]->SetMinimum(0.0);
	      hists_etaTrack[a][b]->SetMaximum(0.5);
	    }
	  hists_etaTrack[a][4]->SetMinimum(0.4);
	  hists_etaTrack[a][4]->SetMaximum(1.6);
	}

      // Draw gaussplots for example eta bins
      // on multi-canvas
      for(int a = 0; a < 3; a++) // Loop over ptratios
	{
	  for(int b = 0; b < 3; b++) // Loop over example pt bins
	    {
	      // Find eta bin of gaussplot
	      int bin = 0;
	      if(  b == 0  )  bin = int(hetaTrack[a]->GetNbinsX()/6);
	      else if(  b == 1  )  bin = int(hetaTrack[a]->GetNbinsX()/3);
	      else if(  b == 2  )  bin = int(hetaTrack[a]->GetNbinsX()/2);
	      float min = hetaTrack[a]->GetXaxis()->GetBinLowEdge(bin);
	      float max = min + hetaTrack[a]->GetXaxis()->GetBinWidth(bin);

	      // Set title according to eta bin
	      sprintf(title,"#gamma-jet Track only, %.1f < E^{#gamma}_{T} < %.1f GeV",min,max);
	      gp_etaTrack[a][b]->SetTitle(title);
	      gp_etaTrack[a][b]->SetXTitle(mPtRatioName[a]);

	      // Set style and line color according to mEtaRatioName
	      gp_etaTrack[a][b]->SetMarkerStyle(markerStyle[a]);
	      gp_etaTrack[a][b]->SetMarkerColor(markerColor[a]);
	      gp_etaTrack[a][b]->SetLineColor(markerColor[a]);
	      gf_etaTrack[a][b]->SetLineColor(markerColor[a]);

	      // Plot gaussplots
	      c2->cd(1+b);
	      gp_etaTrack[a][b]->Draw();
	      gf_etaTrack[a][b]->Draw("same");

	      objToBeWritten.push_back(gp_etaTrack[a][b]);
	      objToBeWritten.push_back(gf_etaTrack[a][b]);
	    } // End of loop over example etabins
	  c2->Draw();
	  ps->NewPage();
	} // End of loop over ptratios

      c1->cd();
      for(int j = 0 ; j < 8 ; ++j) // Loop over control quantities
	{
	  hists_etaTrack[0][j]->Draw("P");
	  hists_etaTrack[0][j]->SetStats(0);
	  hists_etaTrack[1][j]->Draw("P SAME");
	  hists_etaTrack[2][j]->Draw("P SAME");
	  leg->Draw("SAME");
	  c1->SetGrid();
	  c1->Draw();   
	  ps->NewPage(); 

	  objToBeWritten.push_back(hists_etaTrack[0][j]);
	  objToBeWritten.push_back(hists_etaTrack[1][j]);
	  objToBeWritten.push_back(hists_etaTrack[2][j]);
	}
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

  if( plottedQuant.count("log true jet pt") > 0 )
    {
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
	  hists_ptlog[a][4]->SetMinimum(0.4);
	  hists_ptlog[a][4]->SetMaximum(1.6);
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
    }


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

  if( plottedQuant.count("uncorrected jet energy") > 0 )
    {
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
	  hists_energy[a][4]->SetMinimum(0.4);
	  hists_energy[a][4]->SetMaximum(1.6);
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
	      gp_energy[a][b]->SetXTitle(mPtRatioName[a]);

	      // Set style and line color according to mPtRatioName
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

  if( plottedQuant.count("emf") > 0 )
    {
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
	  hists_emf[a][4]->SetMinimum(0.4);
	  hists_emf[a][4]->SetMaximum(1.6);
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
	      gp_emf[a][b]->SetXTitle(mPtRatioName[a]);

	      // Set style and line color according to mPtRatioName
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
    }



  int markerColorRes[4] = { 2,1,4,9 };
  int markerStyleRes[4] = { 22,20,23,21 };
  //JPT(3), JetMET(2) & MyCalibration(1) & Calo(4)
  TLegend* legRes = new TLegend(0.7,0.7,0.96,0.9);
  legRes->SetFillColor(0);
  legRes->AddEntry(RelResPt[0],"My Calibration");
  legRes->AddEntry(RelResPt[1],"JetMET Calibration");
  legRes->AddEntry(RelResPt[2],"JPT Calibration");
  legRes->AddEntry(RelResPt[3],"Calo Jet");

  TH1F* hists_RelResPt[8][8];
  TH1F* gp_RelResPt[8][4];  
  TF1* gf_RelResPt[8][4];     
  for(int a = 0; a < 4; a++)
    {
      RelResPt[a]->SetMarkerStyle(markerStyleRes[a]);
      RelResPt[a]->SetMarkerColor(markerColorRes[a]);
      RelResPt[a]->SetLineColor(markerColorRes[a]);
      RelResPt[a+4]->SetMarkerStyle(markerStyleRes[a]);
      RelResPt[a+4]->SetMarkerColor(markerColorRes[a]);
      RelResPt[a+4]->SetLineColor(markerColorRes[a]);

      Fit2DRes(RelResPt[a],hists_RelResPt[a],gp_RelResPt[a], gf_RelResPt[a]);
      Fit2DRes(RelResPt[a+4],hists_RelResPt[a+4],gp_RelResPt[a+4], gf_RelResPt[a+4]);

      for(int b = 0 ; b < 4 ; ++b)
	{
	  hists_RelResPt[a][b]->SetMinimum(0.5);
	  hists_RelResPt[a][b]->SetMaximum(1.3);
	  hists_RelResPt[a+4][b]->SetMinimum(0.5);
	  hists_RelResPt[a+4][b]->SetMaximum(1.3);
	  ++b;
	  hists_RelResPt[a][b]->SetMinimum(0.0);
	  hists_RelResPt[a][b]->SetMaximum(0.25);
	  hists_RelResPt[a+4][b]->SetMinimum(0.0);
	  hists_RelResPt[a+4][b]->SetMaximum(0.25);
	}
      hists_RelResPt[a][4]->SetMinimum(0.4);
      hists_RelResPt[a][4]->SetMaximum(1.6);
      hists_RelResPt[a+4][4]->SetMinimum(0.4);
      hists_RelResPt[a+4][4]->SetMaximum(1.6);
    }

  c1->cd();
  for(int j = 0 ; j < 5 ; ++j) 
    {
      hists_RelResPt[0][j]->Draw("P");
      objToBeWritten.push_back(hists_RelResPt[0][j]);
      hists_RelResPt[0][j]->SetStats(0);
      hists_RelResPt[1][j]->Draw("P SAME");
      hists_RelResPt[2][j]->Draw("P SAME");
      hists_RelResPt[3][j]->Draw("P SAME");
      objToBeWritten.push_back(hists_RelResPt[1][j]);
      objToBeWritten.push_back(hists_RelResPt[2][j]);
      objToBeWritten.push_back(hists_RelResPt[3][j]);
      legRes->Draw("SAME");
      c1->SetGrid();
      c1->Draw();   
      ps->NewPage(); 
    }
  for(int j = 0 ; j < 5 ; ++j) 
    {
      hists_RelResPt[4][j]->Draw("P");
      objToBeWritten.push_back(hists_RelResPt[4][j]);
      hists_RelResPt[4][j]->SetStats(0);
      hists_RelResPt[5][j]->Draw("P SAME");
      hists_RelResPt[6][j]->Draw("P SAME");
      hists_RelResPt[7][j]->Draw("P SAME");
      objToBeWritten.push_back(hists_RelResPt[5][j]);
      objToBeWritten.push_back(hists_RelResPt[6][j]);
      objToBeWritten.push_back(hists_RelResPt[7][j]);
      legRes->Draw("SAME");
      c1->SetGrid();
      c1->Draw();   
      ps->NewPage(); 
    }


  TH1F* hists_RelResEMF[4][8];
  TH1F* gp_RelResEMF[4][4];  
  TF1* gf_RelResEMF[4][4];     
  for(int a = 0; a < 4; a++)
    {
      RelResEMF[a]->SetMarkerStyle(markerStyleRes[a]);
      RelResEMF[a]->SetMarkerColor(markerColorRes[a]);
      RelResEMF[a]->SetLineColor(markerColorRes[a]);

      Fit2DRes(RelResEMF[a],hists_RelResEMF[a],gp_RelResEMF[a], gf_RelResEMF[a]);

      for(int b = 0 ; b < 4 ; ++b)
	{
	  hists_RelResEMF[a][b]->SetMinimum(0.5);
	  hists_RelResEMF[a][b]->SetMaximum(1.3);
	  ++b;
	  hists_RelResEMF[a][b]->SetMinimum(0.0);
	  hists_RelResEMF[a][b]->SetMaximum(0.3);
	}
      hists_RelResEMF[a][4]->SetMinimum(0.4);
      hists_RelResEMF[a][4]->SetMaximum(1.6);
    }

  c1->cd();
  for(int j = 0 ; j < 5 ; ++j) 
    {
      hists_RelResEMF[0][j]->Draw("P");
      objToBeWritten.push_back(hists_RelResEMF[0][j]);
      hists_RelResEMF[0][j]->SetStats(0);
      hists_RelResEMF[1][j]->Draw("P SAME");
      hists_RelResEMF[2][j]->Draw("P SAME");
      hists_RelResEMF[3][j]->Draw("P SAME");
      objToBeWritten.push_back(hists_RelResEMF[1][j]);
      objToBeWritten.push_back(hists_RelResEMF[2][j]);
      objToBeWritten.push_back(hists_RelResEMF[3][j]);
      legRes->Draw("SAME");
      c1->SetGrid();
      c1->Draw();   
      ps->NewPage(); 
    }

  TH1F* hists_RelResEta[4][8];
  TH1F* gp_RelResEta[4][4];  
  TF1* gf_RelResEta[4][4];     
  for(int a = 0; a < 4; a++)
    {
      RelResEta[a]->SetMarkerStyle(markerStyleRes[a]);
      RelResEta[a]->SetMarkerColor(markerColorRes[a]);
      RelResEta[a]->SetLineColor(markerColorRes[a]);

      Fit2DRes(RelResEta[a],hists_RelResEta[a],gp_RelResEta[a], gf_RelResEta[a]);

      for(int b = 0 ; b < 4 ; ++b)
	{
	  hists_RelResEta[a][b]->SetMinimum(0.5);
	  hists_RelResEta[a][b]->SetMaximum(1.3);
	  ++b;
	  hists_RelResEta[a][b]->SetMinimum(0.0);
	  hists_RelResEta[a][b]->SetMaximum(0.3);
	}
      hists_RelResEta[a][4]->SetMinimum(0.4);
      hists_RelResEta[a][4]->SetMaximum(1.6);
    }

  c1->cd();
  for(int j = 0 ; j < 5 ; ++j) 
    {
      hists_RelResEta[0][j]->Draw("P");
      objToBeWritten.push_back(hists_RelResEta[0][j]);
      hists_RelResEta[0][j]->SetStats(0);
      hists_RelResEta[1][j]->Draw("P SAME");
      hists_RelResEta[2][j]->Draw("P SAME");
      hists_RelResEta[3][j]->Draw("P SAME");
      objToBeWritten.push_back(hists_RelResEta[1][j]);
      objToBeWritten.push_back(hists_RelResEta[2][j]);
      objToBeWritten.push_back(hists_RelResEta[3][j]);
      legRes->Draw("SAME");
      c1->SetGrid();
      c1->Draw();   
      ps->NewPage(); 
    }

  TH1F* hists_RelResTrackMult[4][8];
  TH1F* gp_RelResTrackMult[4][4];  
  TF1* gf_RelResTrackMult[4][4];     
  for(int a = 0; a < 4; a++)
    {
      RelResTrackMult[a]->SetMarkerStyle(markerStyleRes[a]);
      RelResTrackMult[a]->SetMarkerColor(markerColorRes[a]);
      RelResTrackMult[a]->SetLineColor(markerColorRes[a]);

      Fit2DRes(RelResTrackMult[a],hists_RelResTrackMult[a],gp_RelResTrackMult[a], gf_RelResTrackMult[a]);

      for(int b = 0 ; b < 4 ; ++b)
	{
	  hists_RelResTrackMult[a][b]->SetMinimum(0.5);
	  hists_RelResTrackMult[a][b]->SetMaximum(1.3);
	  ++b;
	  hists_RelResTrackMult[a][b]->SetMinimum(0.0);
	  hists_RelResTrackMult[a][b]->SetMaximum(0.3);
	}
      hists_RelResTrackMult[a][4]->SetMinimum(0.4);
      hists_RelResTrackMult[a][4]->SetMaximum(1.6);
    }

  c1->cd();
  for(int j = 0 ; j < 5 ; ++j) 
    {
      hists_RelResTrackMult[0][j]->Draw("P");
      objToBeWritten.push_back(hists_RelResTrackMult[0][j]);
      hists_RelResTrackMult[0][j]->SetStats(0);
      hists_RelResTrackMult[1][j]->Draw("P SAME");
      hists_RelResTrackMult[2][j]->Draw("P SAME");
      hists_RelResTrackMult[3][j]->Draw("P SAME");
      objToBeWritten.push_back(hists_RelResTrackMult[1][j]);
      objToBeWritten.push_back(hists_RelResTrackMult[2][j]);
      objToBeWritten.push_back(hists_RelResTrackMult[3][j]);
      legRes->Draw("SAME");
      c1->SetGrid();
      c1->Draw();   
      ps->NewPage(); 
    }




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
      int i = mPar->GetEtaBin(eta >= 0 ? eta +1 : eta);
      TMeasurement x;
      double* par = mPar->GetTowerParRef(mPar->GetBin(i,0));
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


  // Closing ps file and writing objects to .root file
  ps->Close();
  if( mOutputROOT ) WriteToRootFile(objToBeWritten, "GammaJet");


  // free memory
  for(int i = 0; i < 4; i++)
    {
      delete towerinjet[i];
      delete respvstet[i];
    }
  for(int i = 0; i < 21; i++) delete leadToNext[i];
  for(int i = 0; i < 3; i++)  delete Chi2Plot[i];
  for(int i = 0; i < 2; i++)  delete Chi2Pt[i];
  delete Chi2Eta;
  delete Chi2NoTrack;
  delete Chi2Error;
  for(int i = 0; i < 2; i++)  delete Diff2Pt[i];
  if( plottedQuant.count("eta") > 0 )
    {
      delete hetaTrack[0];
      delete hetaTrack[1];
      delete hetaTrack[2];
      for(int i = 0; i < 3; i++)
	{
	  for(int j = 0; j < 4; j++)
	    {
	      delete hists_etaTrack[i][j];
	      delete hists_etaTrack[i][j+4];
	      delete gp_etaTrack[i][j];
	      delete gf_etaTrack[i][j];
	    }
	}
      for(int i = 0; i < 12; i++)
	{
	  delete heta[i];
	  for(int j = 0; j < 4; j++)
	    {
	      delete hists_eta[i][j];
	      delete hists_eta[i][j+4];
	      delete gp_eta[i][j];
	      delete gf_eta[i][j];
	    }
	}
    } 

  if( plottedQuant.count("uncorrected jet pt") > 0 )
    {
      for(int i = 0; i < 12; i++)
	{
	  delete hpt_uncorr[i];
	  for(int j = 0; j < 4; j++)
	    {
	      delete hists_ptuncorr[i][j];
	      delete hists_ptuncorr[i][j+4];
	      delete gp_ptuncorr[i][j];
	      delete gf_ptuncorr[i][j];
	    }
	}
    }
  delete EtaPhiMap;
  if( plottedQuant.count("true jet pt") > 0 )
    {
      delete hpt[0];
      delete hpt[1];
      delete hpt[2];
      delete hptTrack[0];
      delete hptTrack[1];
      delete hptTrack[2];
      for(int i = 0; i < 3; i++)
	{
	  for(int j = 0; j < 4; j++)
	    {
	      delete hists_pttrue[i][j];
	      delete hists_pttrue[i][j+4];
	      delete gp_pttrue[i][j];
	      delete gf_pttrue[i][j];
	      delete hists_ptTrack[i][j];
	      delete hists_ptTrack[i][j+4];
	      delete gp_ptTrack[i][j];
	      delete gf_ptTrack[i][j];
	    }
	}
    }

  if( plottedQuant.count("uncorrected jet energy") > 0 )
    {
      delete henergy[0];
      delete henergy[1];
      delete henergy[2];
      for(int i = 0; i < 3; i++)
	{
	  for(int j = 0; j < 4; j++)
	    {
	      delete hists_energy[i][j];
	      delete hists_energy[i][j+4];
	      delete gp_energy[i][j];
	      delete gf_energy[i][j];
	    }
	}
    }

  if( plottedQuant.count("log true jet pt") > 0 )
    {
      delete hptlog[0];
      delete hptlog[1];
      delete hptlog[2];
      for(int i = 0; i < 3; i++)
	{
	  for(int j = 0; j < 4; j++)
	    {
	      delete hists_ptlog[i][j];
	      delete hists_ptlog[i][j+4];
	      delete gp_ptlog[i][j];
	      delete gf_ptlog[i][j];
	    }
	}
    }

  if( plottedQuant.count("emf") > 0 )
    {
      delete hemf[0];
      delete hemf[1];
      delete hemf[2];
      for(int i = 0; i < 3; i++)
	{
	  for(int j = 0; j < 4; j++)
	    {
	      delete hists_emf[i][j];
	      delete hists_emf[i][j+4];
	      delete gp_emf[i][j];
	      delete gf_emf[i][j];
	    }
	}
    }

  cout<<"Percentage of Jets with Track Calibration: "<<(double)(track*100)/(track+noTrack)<<"%"<<endl;	  
	  
  delete hptGamma[0];
  delete hptGamma[1];
  delete hptGamma[2];
  delete hptGammaW[0];
  delete hptGammaW[1];
  delete hptGammaW[2];
    
  delete hptGamma2D;
  delete hptGamma2DW;
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
  std::vector<TData*>::const_iterator data_it;
  std::vector<TAbstractData*>::const_iterator it;

  TCanvas * const c1 = new TCanvas("c1","",600,600);
  objToBeDeleted.push_back(c1);
  TPostScript * const ps = new TPostScript("controlplotsGammaJetPerTowerBin.ps",111);
  objToBeDeleted.push_back(ps);

  //one plot per *TOWER* bin!
  for (int eta=0; eta<mPar->GetEtaGranularity();++eta) // Loop over eta bins
    {
      for (int phi=0; phi<mPar->GetPhiGranularity();++phi) // Loop over phi bins
	{
	  int i = mPar->GetBin(eta,phi);

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
	  data_it = mData->begin();
	  for (; data_it != mData->end();++data_it) // loop over all fit-events
	    {
	      if ( (*data_it)->GetType()!=GammaJet ) continue;
	      TAbstractData* ad = dynamic_cast<TAbstractData*>(*data_it);
	      if(! ad) continue;

	      int indexTower = 0; // Index of max tower
	      double Etmax = 0.;
	      double calib_tower_sum = 0.;
	      double tower_sum = 0.; //is equivalent to (*data_it)->GetMess(),
	                             //but since we need the index too, this is faster
	      const std::vector<TAbstractData*>& data_ref = ad->GetRef();
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

	      indexJet += ad->GetIndex();
	      ++ijets;

	      double JetCorr = ad->GetParametrizedMess();

	      fit->Fill(      tower_sum, JetCorr/tower_sum );
	      plot->Fill(     tower_sum, ad->GetTruth()/tower_sum );
	      norm->Fill(     tower_sum ); 
	      plot_jes->Fill( calib_tower_sum, ad->GetTruth()/calib_tower_sum );
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
	  //TF1 * res2 = new TF1(name,mPar->jes_plot_parametrization, 0.5, 400., 3);
	  //objToBeDeleted.push_back(res2);
	  //i = indexJet/ijets - mPar->GetNumberOfTowerParameters();
	  //double * val = mPar->GetJetParRef(i);
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
  if( mOutputROOT ) WriteToRootFile( objToBeWritten, "GammaJetPerTowerBin" );
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
  std::vector<TData*>::const_iterator data_it;
  std::vector<TAbstractData*>::const_iterator it;

  TCanvas * const c1 = new TCanvas("c1","",600,600);
  objToBeDeleted.push_back(c1);
  TPostScript * const ps = new TPostScript("controlplotsGammaJetPerJetBin.ps",111);
  objToBeDeleted.push_back(ps);

  //one plot per *JET* bin!
  for (int eta=0; eta<mPar->GetEtaGranularityJet();++eta)
    {
      for (int phi=0; phi<mPar->GetPhiGranularityJet();++phi)
	{
	  int i = mPar->GetJetBin(eta,phi)*mPar->GetNumberOfJetParametersPerBin() + mPar->GetNumberOfTowerParameters();
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
	  data_it = mData->begin();
	  for (; data_it != mData->end();++data_it)
	    {
	      TAbstractData* ad = dynamic_cast<TAbstractData*>(*data_it);
	      if(! ad) continue;
	      if ( ad->GetType()!=GammaJet ) continue;
	      if ( ad->GetIndex()!= i )	  continue; //event belongs to a wrong bin
	
	      double JetCorr = ad->GetParametrizedMess();

	      double tower_sum = 0.0;
	      double calib_tower_sum = 0.0;

	      const std::vector<TAbstractData*>& data_ref = ad->GetRef(); // Tower
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
	  //TF1 * res2 = new TF1(name,mPar->jes_plot_parametrization, 0.5, 400., 3);
	  //objToBeDeleted.push_back(res2);
	  //i = mPar->GetJetBin(eta, phi);
	  //double * val = mPar->GetJetParRef(i);
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
  if( mOutputROOT ) WriteToRootFile( objToBeWritten, "GammaJetPerJetBin" );
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
  for ( std::vector<TData*>::const_iterator i = mData->begin(); i != mData->end() ; ++i )
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

  if( mOutputROOT ) WriteToRootFile( objToBeWritten, "GammaJetSigmas" );

  for(int i = 0 ; i < nPtBins ; ++i)
    {
      delete gauss_forpt[i];  
      delete gauss_forptcorr[i];  
      //delete f[i];
    }  
  delete text;
}



//---------------------------------------------------------------
//   Dijet-Jet Control Histograms
//---------------------------------------------------------------
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
  TH1F* JetResolutionPlot[12];
  JetResolutionPlot[0] = new TH1F("CaloJetResolution","CaloJetResolution;genJet-CaloJet",200,-300,500);
  JetResolutionPlot[1] = new TH1F("JetMetJetResolution","JetMetJetResolution;genJet-JetMetJet",200,-300,500);
  JetResolutionPlot[2] = new TH1F("JPTJetResolution","JPTJetResolution;genJet-JPTJet",200,-300,500);
  JetResolutionPlot[3] = new TH1F("MyJetResolution","MyJetResolution;genJet-MyJet",200,-300,500);
  JetResolutionPlot[4] = new TH1F("CaloJetRelativeResolution","CaloJetRelativeResolution;(genJet-CaloJet)/genJet",100,-1.5,1.5);
  JetResolutionPlot[5] = new TH1F("JetMetJetRelativeResolution","JetMetJetRelativeResolution;(genJet-JetMetJet)/genJet",100,-1.5,1.5);
  JetResolutionPlot[6] = new TH1F("JPTJetRelativeResolution","JPTJetRelativeResolution;(genJet-JPTJet)/genJet",100,-1.5,1.5);
  JetResolutionPlot[7] = new TH1F("MyJetRelativeResolution","MyJetRelativeResolution;(genJet-MyJet)/genJet",100,-1.5,1.5);
  JetResolutionPlot[8] = new TH1F("CaloJetResponse","CaloJetResponse;CaloJet/genJet",100,0.5,1.5);
  JetResolutionPlot[9] = new TH1F("JetMetJetResponse","JetMetJetResponse;JetMetJet/genJet",100,0.5,1.5);
  JetResolutionPlot[10] = new TH1F("JPTJetResponsen","JPTJetResponse;JPTJet/genJet",100,0.5,1.5);
  JetResolutionPlot[11] = new TH1F("MyJetResponse","MyJetResponse;MyJet/genJet",100,0.5,1.5);
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
  Beta[0] = new TH2F("hBeta0","di-jet;#eta",50,-5,5,100,-0.7,0.7);
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
  Bpt[0] = new TH2F("hBpt0","di-jet;p_{T} [GeV]",80,0,400,100,-0.7,0.7);
  Bpt[1] = (TH2F*)Bpt[0]->Clone("hBpt1");

  TH2F* Benergy[2];
  Benergy[0] = new TH2F("hBenergy0","di-jet;E [GeV]",100,0,400,100,-0.7,0.7);
  Benergy[1] = (TH2F*)Benergy[0]->Clone("hBenergy1");
  
  TH2F* Bemf[2];
  Bemf[0] = new TH2F("hBemf0","di-jet;f_{em} (probe jet)",50,0,1,100,-0.7,0.7);
  Bemf[1] = (TH2F*)Bemf[0]->Clone("hBemf1");

  double bins[101];
  for(int i = 0; i < 101 ; ++i)
    {
      bins[i] = pow(10,(i+32)/40.0);
    }
  TH2F* Bptlog[2];
  Bptlog[0] = new TH2F("hBptlog0","di-jet;p_{T} [GeV]",100,bins,100,-0.7,0.7); 
  Bptlog[1] = (TH2F*)Bptlog[0]->Clone("hBptlog1");

  TH2F *RelResPt[8];
  RelResPt[0] = new TH2F("hRelResPt0","di-jet;p_{T}(gen) [GeV];Jet/genJet",20,0,150,50,0.5,1.5);  
  RelResPt[1] = (TH2F*)RelResPt[0]->Clone("hRelResPt1"); 
  RelResPt[2] = (TH2F*)RelResPt[0]->Clone("hRelResPt2");
  RelResPt[3] = (TH2F*)RelResPt[0]->Clone("hRelResPt3");
  RelResPt[4] = new TH2F("hRelResPt4","di-jet;p_{T}(gen) [GeV];Jet/genJet",20,0,800,50,0.5,1.5);  
  RelResPt[5] = (TH2F*)RelResPt[4]->Clone("hRelResPt5"); 
  RelResPt[6] = (TH2F*)RelResPt[4]->Clone("hRelResPt6");
  RelResPt[7] = (TH2F*)RelResPt[4]->Clone("hRelResPt7");

  TH2F *RelResEMF[4];
  RelResEMF[0] = new TH2F("hRelResEMF0","di-jet;EMF;Jet/genJet",20,0,1,50,0.5,1.5);  
  RelResEMF[1] = (TH2F*)RelResEMF[0]->Clone("hRelResEMF1"); 
  RelResEMF[2] = (TH2F*)RelResEMF[0]->Clone("hRelResEMF2");
  RelResEMF[3] = (TH2F*)RelResEMF[0]->Clone("hRelResEMF3");

  TH2F *RelResEta[4];
  RelResEta[0] = new TH2F("hRelResEta0","di-jet;Eta;Jet/genJet",20,-2.5,2.5,50,0.5,1.5);  
  RelResEta[1] = (TH2F*)RelResEta[0]->Clone("hRelResEta1"); 
  RelResEta[2] = (TH2F*)RelResEta[0]->Clone("hRelResEta2"); 
  RelResEta[3] = (TH2F*)RelResEta[0]->Clone("hRelResEta3");

  TH2F *RelResTrackMult[4];
  RelResTrackMult[0] = new TH2F("hRelResTrackMult0","di-jet;TrackMult;Jet/genJet",20,0,50,50,0.5,1.5);  
  RelResTrackMult[1] = (TH2F*)RelResTrackMult[0]->Clone("hRelResTrackMult1"); 
  RelResTrackMult[2] = (TH2F*)RelResTrackMult[0]->Clone("hRelResTrackMult2");
  RelResTrackMult[3] = (TH2F*)RelResTrackMult[0]->Clone("hRelResTrackMult3");

  //TH1F *TrackOOC[2];
  //TrackOOC[0] = new TH1F("TrackOOC","Out Of Cone Part;JetGenPt; mean OOC part",50,0,300);
  //TrackOOC[1] = (TH1F*)TrackOOC[0]->Clone();
  //TH2F *TrackOOC2D = new TH2F("TrackOOC2D","Out Of Cone Part;JetGenPt; OOC part of genJet",50,0,300,60,0,1.2);

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

  for(int i=0; i<4;++i)
    {
      RelResPt[i]->Sumw2();
      RelResPt[i+4]->Sumw2();
      RelResEMF[i]->Sumw2();
      RelResEta[i]->Sumw2();
      RelResTrackMult[i]->Sumw2();
    }
  int track=0;
  int noTrack=0;

  //loop over all fit-events
  for( std::vector<TData*>::const_iterator i = mData->begin() ; i != mData->end() ; ++i )  
    {
      TData* jj = *i;
      if(jj->GetType() != PtBalance) continue;

      TData_MessMess* jm = (TData_MessMess*) jj;
      double etscale = jm->GetScale();

      //em fraction plots     
      double em = jj->GetMess()->EMF;
      double had = jj->GetMess()->HadF+jj->GetMess()->OutF;

      double etparascale = 0.;
      for(std::vector<TAbstractData*>::const_iterator t = jm->GetRef().begin(); t != jm->GetRef().end(); ++t)
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



      //if(jm->GetTrackuse()) {
	if(true) { 
	TAbstractData* ad = dynamic_cast<TAbstractData*>(jm);
	const std::vector<TAbstractData*>& data_refT = ad->GetRefTrack();
	int TrackMult =  data_refT.size();
	double em1 = jj->GetMess()->EMF;
	double had1 = jj->GetMess()->HadF+jj->GetMess()->OutF;

	track++;
	TJet* jet = (TJet*)(jm->GetMess());
	double genJet = jet->genPt;
	//TrackOOC[0]->Fill(genJet,jm->GetOOC()/genJet);
	//TrackOOC[1]->Fill(genJet);
	//TrackOOC2D->Fill(genJet,jm->GetOOC()/genJet);
	JetResolutionPlot[0]->Fill(genJet- jm->GetMess()->pt);
	JetResolutionPlot[1]->Fill(genJet- jet->L2L3cor * jm->GetMess()->pt);
	JetResolutionPlot[2]->Fill(genJet- jet->L2L3JPTcor * jet->ZSPcor * jet->JPTcor * jm->GetMess()->pt);
	JetResolutionPlot[3]->Fill(genJet- jm->GetParametrizedMess());
	JetResolutionPlot[4]->Fill((genJet- jm->GetMess()->pt)/genJet);
	JetResolutionPlot[5]->Fill((genJet- jet->L2L3cor * jm->GetMess()->pt)/genJet);
	JetResolutionPlot[6]->Fill((genJet- jet->L2L3JPTcor * jet->ZSPcor * jet->JPTcor * jm->GetMess()->pt)/genJet);
	JetResolutionPlot[7]->Fill((genJet- jm->GetParametrizedMess())/genJet);
	JetResolutionPlot[8]->Fill(jm->GetMess()->pt/genJet);
	JetResolutionPlot[9]->Fill( jet->L2L3cor * jm->GetMess()->pt/genJet);
	JetResolutionPlot[10]->Fill( jet->L2L3JPTcor * jet->ZSPcor * jet->JPTcor * jm->GetMess()->pt/genJet);
	JetResolutionPlot[11]->Fill(jm->GetParametrizedMess()/genJet);

	RelResPt[0]->Fill(genJet,(jm->GetParametrizedMess() )/genJet,jj->GetWeight());
	RelResPt[1]->Fill(genJet,(jet->L2L3cor * jm->GetMess()->pt )/genJet,jj->GetWeight());
	RelResPt[2]->Fill(genJet,( jet->L2L3JPTcor * jet->ZSPcor * jet->JPTcor * jm->GetMess()->pt )/genJet,jj->GetWeight());
	RelResPt[3]->Fill(genJet,(jm->GetMess()->pt )/genJet,jj->GetWeight());
	RelResPt[4]->Fill(genJet,(jm->GetParametrizedMess() )/genJet,jj->GetWeight());
	RelResPt[5]->Fill(genJet,(jet->L2L3cor * jm->GetMess()->pt )/genJet,jj->GetWeight());
	RelResPt[6]->Fill(genJet,(jet->L2L3JPTcor * jet->ZSPcor * jet->JPTcor * jm->GetMess()->pt )/genJet,jj->GetWeight());
	RelResPt[7]->Fill(genJet,(jm->GetMess()->pt )/genJet,jj->GetWeight());
	RelResEMF[0]->Fill(em1/(em1+had1),(jm->GetParametrizedMess() )/genJet,jj->GetWeight());
	RelResEMF[1]->Fill(em1/(em1+had1),(jet->L2L3cor * jm->GetMess()->pt )/genJet,jj->GetWeight());
	RelResEMF[2]->Fill(em1/(em1+had1),(jet->L2L3JPTcor * jet->ZSPcor * jet->JPTcor * jm->GetMess()->pt )/genJet,jj->GetWeight());
	RelResEMF[3]->Fill(em1/(em1+had1),(jm->GetMess()->pt )/genJet,jj->GetWeight());
	RelResEta[0]->Fill(jm->GetMess()->eta,(jm->GetParametrizedMess() )/genJet,jj->GetWeight());
	RelResEta[1]->Fill(jm->GetMess()->eta,(jet->L2L3cor * jm->GetMess()->pt )/genJet,jj->GetWeight());
	RelResEta[2]->Fill(jm->GetMess()->eta,(jet->L2L3JPTcor * jet->ZSPcor * jet->JPTcor * jm->GetMess()->pt )/genJet,jj->GetWeight());
	RelResEta[3]->Fill(jm->GetMess()->eta,(jm->GetMess()->pt )/genJet,jj->GetWeight());
	RelResTrackMult[0]->Fill(TrackMult,(jm->GetParametrizedMess() )/genJet,jj->GetWeight());
	RelResTrackMult[1]->Fill(TrackMult,(jet->L2L3cor * jm->GetMess()->pt )/genJet,jj->GetWeight());
	RelResTrackMult[2]->Fill(TrackMult,(jet->L2L3JPTcor * jet->ZSPcor * jet->JPTcor * jm->GetMess()->pt )/genJet,jj->GetWeight());
	RelResTrackMult[3]->Fill(TrackMult,(jm->GetMess()->pt )/genJet,jj->GetWeight());
      }
      else noTrack++;
	//if((*jm->GetSecondaryJets())[0]->GetTrackuse()) {
      if(true) {
	track++;
	TJet* jet = (TJet*)((*jm->GetSecondaryJets())[0]->GetMess());
	double genJet = jet->genPt;
	double parametrizedMess = (*jm->GetSecondaryJets())[0]->GetParametrizedMess();
	double messPt =  (*jm->GetSecondaryJets())[0]->GetMess()->pt;
	TAbstractData* ad = dynamic_cast<TAbstractData*>((*jm->GetSecondaryJets())[0]);
	const std::vector<TAbstractData*>& data_refT = ad->GetRefTrack();
	int TrackMult =  data_refT.size();
	double em2 = (*jm->GetSecondaryJets())[0]->GetMess()->EMF;
	double had2 =(*jm->GetSecondaryJets())[0]->GetMess()->HadF+(*jm->GetSecondaryJets())[0]->GetMess()->OutF;
	double etaJet2 =  (*jm->GetSecondaryJets())[0]->GetMess()->eta;

	//TrackOOC[0]->Fill(genJet,(*jm->GetSecondaryJets())[0]->GetOOC()/genJet);
	//TrackOOC[1]->Fill(genJet);
	//TrackOOC2D->Fill(genJet,(*jm->GetSecondaryJets())[0]->GetOOC()/genJet);
	JetResolutionPlot[0]->Fill(genJet- messPt);
	JetResolutionPlot[1]->Fill(genJet- jet->L2L3cor * messPt);
	JetResolutionPlot[2]->Fill(genJet- jet->L2L3JPTcor * jet->ZSPcor * jet->JPTcor * messPt);
	JetResolutionPlot[3]->Fill(genJet- parametrizedMess);
	JetResolutionPlot[4]->Fill((genJet- messPt)/genJet);
	JetResolutionPlot[5]->Fill((genJet- jet->L2L3cor * messPt)/genJet);
	JetResolutionPlot[6]->Fill((genJet- jet->L2L3JPTcor * jet->ZSPcor * jet->JPTcor * messPt)/genJet);
	JetResolutionPlot[7]->Fill((genJet- parametrizedMess)/genJet);
	JetResolutionPlot[8]->Fill(messPt/genJet);
	JetResolutionPlot[9]->Fill( jet->L2L3cor * messPt/genJet);
	JetResolutionPlot[10]->Fill(jet->L2L3JPTcor * jet->ZSPcor * jet->JPTcor * messPt/genJet);
	JetResolutionPlot[11]->Fill(parametrizedMess/genJet);

	RelResPt[0]->Fill(genJet,(parametrizedMess )/genJet,jj->GetWeight());
	RelResPt[1]->Fill(genJet,(jet->L2L3cor * messPt )/genJet,jj->GetWeight());
	RelResPt[2]->Fill(genJet,(jet->L2L3JPTcor * jet->ZSPcor * jet->JPTcor * messPt )/genJet,jj->GetWeight());
	RelResPt[3]->Fill(genJet,(messPt )/genJet,jj->GetWeight());
	RelResPt[4]->Fill(genJet,(parametrizedMess )/genJet,jj->GetWeight());
	RelResPt[5]->Fill(genJet,(jet->L2L3cor * messPt )/genJet,jj->GetWeight());
	RelResPt[6]->Fill(genJet,(jet->L2L3JPTcor * jet->ZSPcor * jet->JPTcor * messPt )/genJet,jj->GetWeight());
	RelResPt[7]->Fill(genJet,(messPt )/genJet,jj->GetWeight());
	RelResEMF[0]->Fill(em2/(em2+had2),(parametrizedMess )/genJet,jj->GetWeight());
	RelResEMF[1]->Fill(em2/(em2+had2),(jet->L2L3cor * messPt )/genJet,jj->GetWeight());
	RelResEMF[2]->Fill(em2/(em2+had2),(jet->L2L3JPTcor * jet->ZSPcor * jet->JPTcor * messPt )/genJet,jj->GetWeight());
	RelResEMF[3]->Fill(em2/(em2+had2),(messPt )/genJet,jj->GetWeight());
	RelResEta[0]->Fill(etaJet2,(parametrizedMess )/genJet,jj->GetWeight());
	RelResEta[1]->Fill(etaJet2,(jet->L2L3cor * messPt )/genJet,jj->GetWeight());
	RelResEta[2]->Fill(etaJet2,(jet->L2L3JPTcor * jet->ZSPcor * jet->JPTcor * messPt )/genJet,jj->GetWeight());
	RelResEta[3]->Fill(etaJet2,(messPt )/genJet,jj->GetWeight());
	RelResTrackMult[0]->Fill(TrackMult,(parametrizedMess )/genJet,jj->GetWeight());
	RelResTrackMult[1]->Fill(TrackMult,(jet->L2L3cor * messPt )/genJet,jj->GetWeight());
	RelResTrackMult[2]->Fill(TrackMult,(jet->L2L3JPTcor * jet->ZSPcor * jet->JPTcor * messPt )/genJet,jj->GetWeight());
	RelResTrackMult[3]->Fill(TrackMult,(messPt )/genJet,jj->GetWeight());
      }
      else noTrack++;

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


      /*
	for(std::vector<TAbstractData*>::const_iterator t = jj->GetRef().begin(); t != jj->GetRef().end(); ++t)
	{
	TData* tt = *t;
	em  += tt->GetMess()->EMF;
	had += tt->GetMess()->HadF;
	had += tt->GetMess()->OutF;
	}
      */
      
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
  //TrackOOC[0]->Divide(TrackOOC[1]);
  /*
  c1->cd();
  TrackOOC[0]->Draw();
  objToBeWritten.push_back(TrackOOC[0]);
  c1->Draw();
  ps->NewPage(); 

  c1->cd();
  TrackOOC2D->Draw("box");
  objToBeWritten.push_back(TrackOOC2D);
  c1->Draw();
  ps->NewPage(); */

  c1->cd();
  Scale[0]->Draw("box");
  objToBeWritten.push_back(Scale[0]);
  c1->Draw();
  ps->NewPage(); 

  Scale[1]->Draw("box");
  objToBeWritten.push_back(Scale[1]);
  c1->Draw();
  ps->NewPage(); 
  /*
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
  */
  /*
  Bvsdphi->Draw();
  objToBeWritten.push_back(Bvsdphi);
  c1->Draw();
  ps->NewPage();  

  Difvscomb->Draw();
  objToBeWritten.push_back(Difvscomb);
  c1->Draw();
  ps->NewPage();  
  */
  for(int i=0;i<12;++i)
    {
      JetResolutionPlot[i]->Draw();
      objToBeWritten.push_back(JetResolutionPlot[i]);
      c1->Draw();
      ps->NewPage(); 
    }

  // Control quantities vs eta
  TH1F* hists_beta[8][8];
  TH1F* gp_beta[8][4];  
  TF1* gf_beta[8][4];     

  TLegend* leg = new TLegend(0.7,0.7,0.96,0.9);
  leg->SetFillColor(0);
  leg->AddEntry(Beta[1],"B before fit");
  leg->AddEntry(Beta[0],"B after fit");

  int markerColor[4] = { 2,1,4,9 };
  int markerStyle[4] = { 22,20,23,21 };
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
	  hists_beta[i+a][4]->SetMinimum(-0.6);
	  hists_beta[i+a][4]->SetMaximum(0.6);
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

	      // Set style and line color according to mPtRatioName
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
      hists_pt[a][4]->SetMinimum(-0.6);
      hists_pt[a][4]->SetMaximum(0.6);
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

	  // Set style and line color according to mPtRatioName
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
      hists_ptlog[a][4]->SetMinimum(-0.6);
      hists_ptlog[a][4]->SetMaximum(0.6);
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
      hists_energy[a][4]->SetMinimum(-0.6);
      hists_energy[a][4]->SetMaximum(0.6);
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

	  // Set style and line color according to mPtRatioName
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
      hists_emf[a][4]->SetMinimum(-0.6);
      hists_emf[a][4]->SetMaximum(0.6);
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

	  // Set style and line color according to mPtRatioName
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
  


  //JPT(3), JetMET(2) & MyCalibration(1) & Calo(4)
  TLegend* legRes = new TLegend(0.7,0.7,0.96,0.9);
  legRes->SetFillColor(0);
  legRes->AddEntry(RelResPt[0],"My Calibration");
  legRes->AddEntry(RelResPt[1],"JetMET Calibration");
  legRes->AddEntry(RelResPt[2],"JPT Calibration");
  legRes->AddEntry(RelResPt[3],"Calo Jet");

  TH1F* hists_RelResPt[8][8];
  TH1F* gp_RelResPt[8][4];  
  TF1* gf_RelResPt[8][4];     
  for(int a = 0; a < 4; a++)
    {
      RelResPt[a]->SetMarkerStyle(markerStyle[a]);
      RelResPt[a]->SetMarkerColor(markerColor[a]);
      RelResPt[a]->SetLineColor(markerColor[a]);
      RelResPt[a+4]->SetMarkerStyle(markerStyle[a]);
      RelResPt[a+4]->SetMarkerColor(markerColor[a]);
      RelResPt[a+4]->SetLineColor(markerColor[a]);

      Fit2DRes(RelResPt[a],hists_RelResPt[a],gp_RelResPt[a], gf_RelResPt[a]);
      Fit2DRes(RelResPt[a+4],hists_RelResPt[a+4],gp_RelResPt[a+4], gf_RelResPt[a+4]);

      for(int b = 0 ; b < 4 ; ++b)
	{
	  hists_RelResPt[a][b]->SetMinimum(0.5);
	  hists_RelResPt[a][b]->SetMaximum(1.3);
	  hists_RelResPt[a+4][b]->SetMinimum(0.5);
	  hists_RelResPt[a+4][b]->SetMaximum(1.3);
	  ++b;
	  hists_RelResPt[a][b]->SetMinimum(0.0);
	  hists_RelResPt[a][b]->SetMaximum(0.25);
	  hists_RelResPt[a+4][b]->SetMinimum(0.0);
	  hists_RelResPt[a+4][b]->SetMaximum(0.25);
	}
      hists_RelResPt[a][4]->SetMinimum(0.4);
      hists_RelResPt[a][4]->SetMaximum(1.6);
      hists_RelResPt[a+4][4]->SetMinimum(0.4);
      hists_RelResPt[a+4][4]->SetMaximum(1.6);
    }

  c1->cd();
  for(int j = 0 ; j < 5 ; ++j) 
    {
      hists_RelResPt[0][j]->Draw("P");
      objToBeWritten.push_back(hists_RelResPt[0][j]);
      hists_RelResPt[0][j]->SetStats(0);
      hists_RelResPt[1][j]->Draw("P SAME");
      hists_RelResPt[2][j]->Draw("P SAME");
      hists_RelResPt[3][j]->Draw("P SAME");
      objToBeWritten.push_back(hists_RelResPt[1][j]);
      objToBeWritten.push_back(hists_RelResPt[2][j]);
      objToBeWritten.push_back(hists_RelResPt[3][j]);
      legRes->Draw("SAME");
      c1->SetGrid();
      c1->Draw();   
      ps->NewPage(); 
    }
  for(int j = 0 ; j < 5 ; ++j) 
    {
      hists_RelResPt[4][j]->Draw("P");
      objToBeWritten.push_back(hists_RelResPt[4][j]);
      hists_RelResPt[4][j]->SetStats(0);
      hists_RelResPt[5][j]->Draw("P SAME");
      hists_RelResPt[6][j]->Draw("P SAME");
      hists_RelResPt[7][j]->Draw("P SAME");
      objToBeWritten.push_back(hists_RelResPt[5][j]);
      objToBeWritten.push_back(hists_RelResPt[6][j]);
      objToBeWritten.push_back(hists_RelResPt[7][j]);
      legRes->Draw("SAME");
      c1->SetGrid();
      c1->Draw();   
      ps->NewPage(); 
    }

  TH1F* hists_RelResEMF[4][8];
  TH1F* gp_RelResEMF[4][4];  
  TF1* gf_RelResEMF[4][4];     
  for(int a = 0; a < 4; a++)
    {
      RelResEMF[a]->SetMarkerStyle(markerStyle[a]);
      RelResEMF[a]->SetMarkerColor(markerColor[a]);
      RelResEMF[a]->SetLineColor(markerColor[a]);

      Fit2DRes(RelResEMF[a],hists_RelResEMF[a],gp_RelResEMF[a], gf_RelResEMF[a]);

      for(int b = 0 ; b < 4 ; ++b)
	{
	  hists_RelResEMF[a][b]->SetMinimum(0.5);
	  hists_RelResEMF[a][b]->SetMaximum(1.3);
	  ++b;
	  hists_RelResEMF[a][b]->SetMinimum(0.0);
	  hists_RelResEMF[a][b]->SetMaximum(0.3);
	}
      hists_RelResEMF[a][4]->SetMinimum(0.4);
      hists_RelResEMF[a][4]->SetMaximum(1.6);
    }

  c1->cd();
  for(int j = 0 ; j < 5 ; ++j) 
    {
      hists_RelResEMF[0][j]->Draw("P");
      objToBeWritten.push_back(hists_RelResEMF[0][j]);
      hists_RelResEMF[0][j]->SetStats(0);
      hists_RelResEMF[1][j]->Draw("P SAME");
      hists_RelResEMF[2][j]->Draw("P SAME");
      hists_RelResEMF[3][j]->Draw("P SAME");
      objToBeWritten.push_back(hists_RelResEMF[1][j]);
      objToBeWritten.push_back(hists_RelResEMF[2][j]);
      objToBeWritten.push_back(hists_RelResEMF[3][j]);
      legRes->Draw("SAME");
      c1->SetGrid();
      c1->Draw();   
      ps->NewPage(); 
    }

  TH1F* hists_RelResEta[4][8];
  TH1F* gp_RelResEta[4][4];  
  TF1* gf_RelResEta[4][4];     
  for(int a = 0; a < 4; a++)
    {
      RelResEta[a]->SetMarkerStyle(markerStyle[a]);
      RelResEta[a]->SetMarkerColor(markerColor[a]);
      RelResEta[a]->SetLineColor(markerColor[a]);

      Fit2DRes(RelResEta[a],hists_RelResEta[a],gp_RelResEta[a], gf_RelResEta[a]);

      for(int b = 0 ; b < 4 ; ++b)
	{
	  hists_RelResEta[a][b]->SetMinimum(0.5);
	  hists_RelResEta[a][b]->SetMaximum(1.3);
	  ++b;
	  hists_RelResEta[a][b]->SetMinimum(0.0);
	  hists_RelResEta[a][b]->SetMaximum(0.3);
	}
      hists_RelResEta[a][4]->SetMinimum(0.4);
      hists_RelResEta[a][4]->SetMaximum(1.6);
    }

  c1->cd();
  for(int j = 0 ; j < 5 ; ++j) 
    {
      hists_RelResEta[0][j]->Draw("P");
      objToBeWritten.push_back(hists_RelResEta[0][j]);
      hists_RelResEta[0][j]->SetStats(0);
      hists_RelResEta[1][j]->Draw("P SAME");
      hists_RelResEta[2][j]->Draw("P SAME");
      hists_RelResEta[3][j]->Draw("P SAME");
      objToBeWritten.push_back(hists_RelResEta[1][j]);
      objToBeWritten.push_back(hists_RelResEta[2][j]);
      objToBeWritten.push_back(hists_RelResEta[3][j]);
      legRes->Draw("SAME");
      c1->SetGrid();
      c1->Draw();   
      ps->NewPage(); 
    }

  TH1F* hists_RelResTrackMult[4][8];
  TH1F* gp_RelResTrackMult[4][4];  
  TF1* gf_RelResTrackMult[4][4];     
  for(int a = 0; a < 4; a++)
    {
      RelResTrackMult[a]->SetMarkerStyle(markerStyle[a]);
      RelResTrackMult[a]->SetMarkerColor(markerColor[a]);
      RelResTrackMult[a]->SetLineColor(markerColor[a]);

      Fit2DRes(RelResTrackMult[a],hists_RelResTrackMult[a],gp_RelResTrackMult[a], gf_RelResTrackMult[a]);

      for(int b = 0 ; b < 4 ; ++b)
	{
	  hists_RelResTrackMult[a][b]->SetMinimum(0.5);
	  hists_RelResTrackMult[a][b]->SetMaximum(1.3);
	  ++b;
	  hists_RelResTrackMult[a][b]->SetMinimum(0.0);
	  hists_RelResTrackMult[a][b]->SetMaximum(0.3);
	}
      hists_RelResTrackMult[a][4]->SetMinimum(0.4);
      hists_RelResTrackMult[a][4]->SetMaximum(1.6);
    }

  c1->cd();
  for(int j = 0 ; j < 5 ; ++j) 
    {
      hists_RelResTrackMult[0][j]->Draw("P");
      objToBeWritten.push_back(hists_RelResTrackMult[0][j]);
      hists_RelResTrackMult[0][j]->SetStats(0);
      hists_RelResTrackMult[1][j]->Draw("P SAME");
      hists_RelResTrackMult[2][j]->Draw("P SAME");
      hists_RelResTrackMult[3][j]->Draw("P SAME");
      objToBeWritten.push_back(hists_RelResTrackMult[1][j]);
      objToBeWritten.push_back(hists_RelResTrackMult[2][j]);
      objToBeWritten.push_back(hists_RelResTrackMult[3][j]);
      legRes->Draw("SAME");
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

  if( mOutputROOT ) WriteToRootFile( objToBeWritten, "DiJet" );

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
  for(int i = 0; i < 3; i++)
    {
      delete RelResPt[i];
      delete RelResPt[i+3];
      delete RelResEMF[i];
      delete RelResEta[i];
      delete RelResTrackMult[i];
    }
  for(int i = 0 ; i < 3 ; ++i)
    {
      for(int j = 0 ; j < 4 ; ++j)
	{
	  delete hists_RelResPt[i][j];
	  delete hists_RelResPt[i][j+4];
	  delete gp_RelResPt[i][j];
	  delete gf_RelResPt[i][j];
	  delete hists_RelResPt[i+3][j];
	  delete hists_RelResPt[i+3][j+4];
	  delete gp_RelResPt[i+3][j];
	  delete gf_RelResPt[i+3][j];

	  delete hists_RelResEMF[i][j];
	  delete hists_RelResEMF[i][j+4];
	  delete gp_RelResEMF[i][j];
	  delete gf_RelResEMF[i][j];

	  delete hists_RelResEta[i][j];
	  delete hists_RelResEta[i][j+4];
	  delete gp_RelResEta[i][j];
	  delete gf_RelResEta[i][j];

	  delete hists_RelResTrackMult[i][j];
	  delete hists_RelResTrackMult[i][j+4];
	  delete gp_RelResTrackMult[i][j];
	  delete gf_RelResTrackMult[i][j];
	}	
    }


  cout<<"Percentage of Jets with Track Calibration: "<<(double)(track*100)/(track+noTrack)<<"%"<<endl;	

  delete leg;
  delete legRes;
}



//---------------------------------------------------------------
//   Top Control Histograms
//---------------------------------------------------------------
void TControlPlots::MakeControlPlotsTop()
{
  std::vector<TObject*> objToBeWritten;

  TCanvas * const c = new TCanvas("c","",600,600);

  TPostScript * const ps = new TPostScript("controlplotsTop.ps",111);

  bool printEps = false;  // just as a temporary solution. don't we want to support this in general?

  // book hists

  TH1F* scale  = new TH1F("scale" , "Scale" , 200,    0, 200);
  TH1F* weight = new TH1F("weight", "Weight", 200,    0, 200);
  TH1F* truth  = new TH1F("truth" , "Truth" , 200,    0, 200);
  TH1F* pt     = new TH1F("pt"    , "Pt"    , 200,    0, 200);
  TH1F* eta    = new TH1F("eta"   , "Eta"   ,  80,  -4.,  4.);
  TH1F* phi    = new TH1F("phi"   , "Phi"   ,  68, -3.4, 3.4);

  TH1F* meanPt  = new TH1F("meanPt" , "MeanPt" , 200,    0, 200);
  TH1F* meanEta = new TH1F("meanEta", "MeanEta",  80,  -4.,  4.);

  TH1F* invMass         [2];
  TH1F* messTruth       [2];
  TProfile* messTruthPt [2];
  TProfile* messTruthEta[2];
  TProfile* responsePt  [2];
  TProfile* responseEta [2];
  TString suffix[2] = { "Before", "After" };
  for(unsigned a=0; a<2; a++){
    invMass     [a] = new TH1F("invMass"  +suffix[a], "",  40, 0., 200.);
    messTruth   [a] = new TH1F("messTruth"+suffix[a], "",  40, 0.,   2.);
    messTruthPt [a] = new TProfile("messTruthPt" +suffix[a], "", 35,    0, 140);
    messTruthEta[a] = new TProfile("messTruthEta"+suffix[a], "", 40,  -4.,  4.);
    responsePt  [a] = new TProfile("responsePt"  +suffix[a], "", 35,    0, 140);
    responseEta [a] = new TProfile("responseEta" +suffix[a], "", 40,  -4.,  4.);
  }

  //loop over all fit-events and fill hists

  for( std::vector<TData*>::const_iterator i = mData->begin() ; i != mData->end() ; ++i )  
    {

      TData* data = *i;
      if(data->GetType() != InvMass) continue;

      std::vector<TJet*> jets;
      
      TData_InvMass2 *invM2 = dynamic_cast<TData_InvMass2*>(data);	
      TwoJetsInvMassEvent* ev = (TwoJetsInvMassEvent*)data;
      if(invM2) {	
	jets.push_back( (TJet*) invM2->GetMess() );
	for(unsigned j=0; j<invM2->MultMessSize(); j++)
	  jets.push_back( (TJet*) (*invM2->GetSecondaryJets())[j]->GetMess() );
      } 
      else if(ev) {
	jets.push_back( ev->GetJet1() );
	jets.push_back( ev->GetJet2() );
      } else { 
      continue;
      }

      TLorentzVector jet4Vec, combined4Vec;
      for(unsigned j=0; j<jets.size(); j++) {
	jet4Vec.SetPtEtaPhiE( jets[j]->pt, jets[j]->eta, jets[j]->phi, jets[j]->E );
	if(j==0) combined4Vec = jet4Vec;
	else combined4Vec += jet4Vec;
      }
      double s,t,w;
      
      if(invM2) {
	s = invM2->GetScale();
	w = invM2->GetWeight();
	t = invM2->GetTruth();
      }
      if(ev) {
	s = ev->GetTruth();
	w = ev->GetWeight();
	t = ev->GetTruth();
      }
      scale ->Fill(s);
      weight->Fill(w);
      truth ->Fill(t);
      
      invMass  [0]->Fill( combined4Vec.M()   );
      messTruth[0]->Fill( combined4Vec.M()/t );

      invMass  [1]->Fill( invM2 ? invM2->GetMessCombination()   : ev->correctedMass()   );
      messTruth[1]->Fill( invM2 ? invM2->GetMessCombination()/t : ev->correctedMass()/t );

      double mPt  = 0.;
      double mEta = 0.;

      for(unsigned j=0; j<jets.size(); j++) {
	pt ->Fill( jets[j]->pt  );
	eta->Fill( jets[j]->eta );
	phi->Fill( jets[j]->phi );

	mPt  += jets[j]->pt;
	mEta += jets[j]->eta;

	messTruthPt [0]->Fill( jets[j]->pt  , combined4Vec.M()/t );
	messTruthEta[0]->Fill( jets[j]->eta , combined4Vec.M()/t );
	
	messTruthPt [1]->Fill( jets[j]->pt  , invM2 ? invM2->GetMessCombination()/t : ev->correctedMass()/t );
	messTruthEta[1]->Fill( jets[j]->eta , invM2 ? invM2->GetMessCombination()/t : ev->correctedMass()/t );

	double response[2];
	if(j==0) {
	  response[0] = invM2 ? invM2->GetMess()->pt / jets[j]->genPt : 
	                           ev->GetMess()->pt / jets[j]->genPt;
	  response[1] = invM2 ? invM2->GetParametrizedMess() / jets[j]->genPt : 
	                           ev->GetParametrizedMess() / jets[j]->genPt;
	}
	else {
	  if(invM2) {
	    const TData_MessMess* mm = (*invM2->GetSecondaryJets())[j-1];
	    response[0] = mm->GetMess()->pt         / jets[j]->genPt;
	    response[1] = mm->GetParametrizedMess() / jets[j]->genPt;
	  } else {
	    response[0] = ev->GetMess2()->pt         / jets[j]->genPt;
	    response[1] = ev->GetParametrizedMess2() / jets[j]->genPt;
	  }
	}
	for(unsigned a=0; a<2; a++){
	  responsePt [a]->Fill( jets[j]->pt  , response[a] );
	  responseEta[a]->Fill( jets[j]->eta , response[a] );
	}
      }

      mPt  /= jets.size();
      mEta /= jets.size();

      meanPt ->Fill( mPt  );
      meanEta->Fill( mEta );

    }  //End of loop over all fit-events

  // configure hists

  int markerColor[2] = { 2, 1 };
  int markerStyle[2] = { 22, 20 };

  TPaveText *paveText[2];
  paveText[0] = new TPaveText(0.7, 0.6 , 0.91, 0.75, "NDC");
  paveText[1] = new TPaveText(0.7, 0.45, 0.91, 0.6 , "NDC");

  for(unsigned a=0; a<2; a++){

    invMass     [a]->Fit( "gaus", "Q0" );
    messTruth   [a]->Fit( "gaus", "Q0" );

    double invMassMu      = invMass[a]->GetFunction("gaus")->GetParameter(1);
    double invMassSigma   = invMass[a]->GetFunction("gaus")->GetParameter(2);

    double messTruthMu    = messTruth[a]->GetFunction("gaus")->GetParameter(1);
    double messTruthSigma = messTruth[a]->GetFunction("gaus")->GetParameter(2);

    invMass     [a]->Fit( "gaus", "Q", "", invMassMu  -2*invMassSigma  , invMassMu  +2*invMassSigma   );
    messTruth   [a]->Fit( "gaus", "Q", "", messTruthMu-2*messTruthSigma, messTruthMu+2*messTruthSigma );

    invMass     [a]->SetLineColor( markerColor[a] );
    messTruth   [a]->SetLineColor( markerColor[a] );

    invMass     [a]->GetFunction("gaus")->SetLineColor( markerColor[a] );
    messTruth   [a]->GetFunction("gaus")->SetLineColor( markerColor[a] );

    invMass     [a]->GetFunction("gaus")->SetLineWidth( 1 );
    messTruth   [a]->GetFunction("gaus")->SetLineWidth( 1 );

    invMass     [a]->SetMarkerColor( markerColor[a] );
    messTruth   [a]->SetMarkerColor( markerColor[a] );
    messTruthPt [a]->SetMarkerColor( markerColor[a] );
    messTruthEta[a]->SetMarkerColor( markerColor[a] );

    invMass     [a]->SetMarkerStyle( markerStyle[a] );
    messTruth   [a]->SetMarkerStyle( markerStyle[a] );
    messTruthPt [a]->SetMarkerStyle( markerStyle[a] );
    messTruthEta[a]->SetMarkerStyle( markerStyle[a] );

    invMass     [a]->SetMarkerSize( 1.5 );
    messTruth   [a]->SetMarkerSize( 1.5 );
    messTruthPt [a]->SetMarkerSize( 1.5 );
    messTruthEta[a]->SetMarkerSize( 1.5 );

    paveText[a]->SetFillColor( 0 );
    paveText[a]->SetBorderSize( 1 );
    paveText[a]->SetTextColor( markerColor[a] );
    paveText[a]->SetTextAlign( 12 );

    double mu    = invMass[a]->GetFunction("gaus")->GetParameter(1);
    double sigma = invMass[a]->GetFunction("gaus")->GetParameter(2);
    double relSigma = sigma / mu;

    char *tmpTxt = new char[100];

    sprintf(tmpTxt, "#mu = %4.1f GeV", mu);
    paveText[a]->AddText(tmpTxt);
    sprintf(tmpTxt, "#sigma = %4.1f GeV", sigma);
    paveText[a]->AddText(tmpTxt);
    sprintf(tmpTxt, "#sigma/#mu = %4.2f", relSigma);
    paveText[a]->AddText(tmpTxt);

    invMass     [a]->SetStats( 0 );
    messTruth   [a]->SetStats( 0 );
    messTruthPt [a]->SetStats( 0 );
    messTruthEta[a]->SetStats( 0 );
    
    invMass     [a]->SetXTitle( "invariant mass [GeV]" );
    messTruth   [a]->SetXTitle( "measurement/truth" );
    messTruthPt [a]->SetXTitle( "p_{T} [GeV]" );
    messTruthEta[a]->SetXTitle( "#eta" );

    invMass     [a]->SetYTitle( "events" );
    messTruth   [a]->SetYTitle( "events" );
    messTruthPt [a]->SetYTitle( "measurement/truth" );
    messTruthEta[a]->SetYTitle( "measurement/truth" );

    messTruthPt [a]->SetMinimum( 0.4 );
    messTruthEta[a]->SetMinimum( 0.4 );

    messTruthPt [a]->SetMaximum( 1.6 );
    messTruthEta[a]->SetMaximum( 1.6 );

    responsePt [a]->SetMarkerColor( markerColor[a] );
    responseEta[a]->SetMarkerColor( markerColor[a] );

    responsePt [a]->SetMarkerStyle( markerStyle[a] );
    responseEta[a]->SetMarkerStyle( markerStyle[a] );

    responsePt [a]->SetMarkerSize( 1.5 );
    responseEta[a]->SetMarkerSize( 1.5 );
    
    responsePt [a]->SetStats( 0 );
    responseEta[a]->SetStats( 0 );
    
    responsePt [a]->SetXTitle( "p_{T} [GeV]" );
    responseEta[a]->SetXTitle( "#eta" );

    responsePt [a]->SetYTitle( "p_{T} (rec) / p_{T} (gen)" );
    responseEta[a]->SetYTitle( "p_{T} (rec) / p_{T} (gen)" );

    responsePt [a]->SetMinimum( 0. );
    responseEta[a]->SetMinimum( 0. );

    responsePt [a]->SetMaximum( 1.2 );
    responseEta[a]->SetMaximum( 1.2 );

  }

  // create a legend

  TLegend* legend = new TLegend(0.7, 0.75, 0.91, 0.85);
  legend->SetFillColor(0);
  legend->AddEntry(invMass[0],"before fit");
  legend->AddEntry(invMass[1],"after fit");

  // create a line

  TLine* line = new TLine();
  line->SetLineStyle(2);

  // draw hists

  c->cd(1);

  scale->Draw();
  c->Draw();
  if(printEps) c->Print("top_scale.eps");
  ps->NewPage();

  weight->Draw();
  c->Draw();
  if(printEps) c->Print("top_weight.eps");
  ps->NewPage();

  truth->Draw();
  c->Draw();
  if(printEps) c->Print("top_truth.eps");
  ps->NewPage();

  pt->Draw();
  c->Draw();
  if(printEps) c->Print("top_pt.eps");
  ps->NewPage();

  eta->Draw();
  c->Draw();
  if(printEps) c->Print("top_eta.eps");
  ps->NewPage();

  phi->Draw();
  c->Draw();
  if(printEps) c->Print("top_phi.eps");
  ps->NewPage();

  meanPt->Draw();
  c->Draw();
  if(printEps) c->Print("top_meanPt.eps");
  ps->NewPage();

  meanEta->Draw();
  c->Draw();
  if(printEps) c->Print("top_meanEta.eps");
  ps->NewPage();

  invMass[0]->Draw("p");
  invMass[1]->Draw("p same");
  legend->Draw("same");
  paveText[0]->Draw();
  paveText[1]->Draw();
  c->Draw();
  if(printEps) c->Print("top_invMass.eps");
  ps->NewPage();

  messTruth[0]->Draw("p");
  messTruth[1]->Draw("p same");
  legend->Draw("same");
  c->Draw();
  if(printEps) c->Print("top_messTruth.eps");
  ps->NewPage();

  messTruthPt[0]->Draw("p");
  messTruthPt[1]->Draw("p same");
  legend->Draw("same");
  line->DrawLine(messTruthPt[0]->GetXaxis()->GetXmin(), 1., messTruthPt[0]->GetXaxis()->GetXmax(), 1.);
  c->Draw();
  if(printEps) c->Print("top_messTruthPt.eps");
  ps->NewPage();

  messTruthEta[0]->Draw("p");
  messTruthEta[1]->Draw("p same");
  legend->Draw("same");
  line->DrawLine(messTruthEta[0]->GetXaxis()->GetXmin(), 1., messTruthEta[0]->GetXaxis()->GetXmax(), 1.);
  c->Draw();
  if(printEps) c->Print("top_messTruthEta.eps");
  ps->NewPage();

  responsePt[0]->Draw("p");
  responsePt[1]->Draw("p same");
  legend->Draw("same");
  line->DrawLine(responsePt[0]->GetXaxis()->GetXmin(), 1., responsePt[0]->GetXaxis()->GetXmax(), 1.);
  c->Draw();
  if(printEps) c->Print("top_responsePt.eps");
  ps->NewPage();

  responseEta[0]->Draw("p");
  responseEta[1]->Draw("p same");
  legend->Draw("same");
  line->DrawLine(responseEta[0]->GetXaxis()->GetXmin(), 1., responseEta[0]->GetXaxis()->GetXmax(), 1.);
  c->Draw();
  if(printEps) c->Print("top_responseEta.eps");

  ps->Close();

  // write to root-file

  objToBeWritten.push_back( scale  );
  objToBeWritten.push_back( weight );
  objToBeWritten.push_back( truth  );
  objToBeWritten.push_back( pt  );
  objToBeWritten.push_back( eta );
  objToBeWritten.push_back( phi );
  objToBeWritten.push_back( meanPt  );
  objToBeWritten.push_back( meanEta );
  for(unsigned a=0; a<2; a++){
    objToBeWritten.push_back( invMass     [a] );
    objToBeWritten.push_back( messTruth   [a] );
    objToBeWritten.push_back( messTruthPt [a] );
    objToBeWritten.push_back( messTruthEta[a] );
    objToBeWritten.push_back( responsePt  [a] );
    objToBeWritten.push_back( responseEta [a] );
  }

  if( mOutputROOT ) WriteToRootFile( objToBeWritten, "Top" );

  // clean up

  delete scale;
  delete weight;
  delete truth;
  delete pt;
  delete eta;
  delete phi;
  delete meanPt;
  delete meanEta;
  for(unsigned a=0; a<2; a++){
    delete paveText    [a];
    delete invMass     [a];
    delete messTruth   [a];
    delete messTruthPt [a];
    delete messTruthEta[a];
    delete responsePt  [a];
    delete responseEta [a];
  }

  delete legend;
  delete line;

}



//---------------------------------------------------------------
//   Vary parameters around fitted value and plot
//   chi2 profile. The method generates one 1-dim plot
//   for each parameter in which that parameter is varied
//   and all other parameters are left as in TParameters::GetPars().
//   Additionally, correlations between parameters are
//   shown in 2-dim plots, where two parameters are varied.
//   So far, only correlations between the 1st and 2nd, the
//   3rd and 4th and so on are plotted. This has to be extended
//   in a clever way...
//
//   NOTE: Chi2 is calculated w/o scaling of residuals
//   ( = scaling scheme '0' )
//---------------------------------------------------------------
void TControlPlots::MakeControlPlotsParameterScan()
{
  // Put scaling scheme to no scaling for the time being
  TData::ScaleResidual = &TData::ScaleNone;

  std::vector<TObject*> objToBeWritten;

  TCanvas * const c1 = new TCanvas("c1","",600,600);
  TPostScript * const ps = new TPostScript("controlplotsParameterScan.ps",111);


  // Loop over parameters and vary one parameter
  // while the others are left as in TParameter::GetPars()
  TGraph *gParScan[mPar->GetNumberOfParameters()];
  TH1F *hFrame[mPar->GetNumberOfParameters()];
  TLine *line[mPar->GetNumberOfParameters()];
  for(int i = 0; i < mPar->GetNumberOfParameters(); i++)
    {
      // Store original value of parameter i
      double origPar = mPar->GetPars()[i];

      // Vary parameter and get chi2
      double x_param[21];
      double y_chi2[21];
      for(int a = 0; a < 21; a++)
	{
	  x_param[a] = 0.;
	  y_chi2[a] = 0.;
	}
      for(int a = 0; a < 21; a++)
	{
	  double variedPar = origPar - 0.1 + 0.01*a;
	  x_param[a] = variedPar;
	  mPar->GetPars()[i] = variedPar;
	  for( std::vector<TData*>::const_iterator it = mData->begin();  it < mData->end();  ++it )
	    {
	      y_chi2[a] += (*it)->chi2();
	    }
	}

      // Reset original parameter
      mPar->GetPars()[i] = origPar;


      // Draw graphs
      TString name = "gParScan";
      name += i;
      TString title = "Parameter ";
      title += i;
      if( i < mPar->GetNumberOfTowerParameters() )
	{
	  title += " (tower parameter ";
	  title += i;
	}
      else
	{
	  title += " (jet parameter ";
	  title += i - (mPar->GetNumberOfTowerParameters());
	}
      title += ");Parameter p_{";
      title += i;
      title += "};#chi^{2}";
      gParScan[i] = new TGraph(21,x_param,y_chi2);
      gParScan[i]->SetName(name);
      gParScan[i]->SetTitle(title);
      gParScan[i]->SetMarkerStyle(20);
      objToBeWritten.push_back(gParScan[i]);

      name = "hParScanFrame";
      name += i;
      hFrame[i] = new TH1F(name,title,1,x_param[0]-0.01,x_param[20]+0.01);
      double y_min = *min_element(y_chi2,y_chi2+21);
      double y_max = *max_element(y_chi2,y_chi2+21);
      double range = y_max - y_min;
      hFrame[i]->GetYaxis()->SetRangeUser(y_min - range/20, y_max + range/20);

      line[i] = new TLine(origPar,hFrame[i]->GetMinimum(),origPar,hFrame[i]->GetMaximum());
      line[i]->SetLineStyle(2);
      line[i]->SetLineWidth(1);

      c1->cd();
      hFrame[i]->Draw();
      gParScan[i]->Draw("Psame");
      line[i]->Draw("same");
      c1->Draw();   
    }



  // Correlations between 2 parameters
  int nPlots = mPar->GetNumberOfParameters()/2;
  TGraph2D *gParScan2D[nPlots];
  for(int i = 0; i < nPlots; i++)
    {
      int parIdx = 2*i;

      // Store original value of parameters i and i+1
      double origParX = mPar->GetPars()[parIdx];
      double origParY = mPar->GetPars()[parIdx+1];

      // Vary parameter and get chi2
      double x_param[441];
      double y_param[441];
      double z_chi2[441];
      for(int a = 0; a < 441; a++)
	{
	  x_param[a] = 0.;
	  y_param[a] = 0.;
	  z_chi2[a] = 0.;
	}
      for(int a = 0; a < 21; a++)
	{
	  double variedParX = origParX - 0.1 + 0.01*a;
	  mPar->GetPars()[parIdx] = variedParX;
	  for(int b = 0; b < 21; b++)
	    {
	      double variedParY = origParY - 0.1 + 0.01*b;
	      mPar->GetPars()[parIdx+1] = variedParY;

	      int pointIdx = b + 21*a;
	      x_param[pointIdx] = variedParX;
	      y_param[pointIdx] = variedParY;

	      for( std::vector<TData*>::const_iterator it = mData->begin();  it < mData->end();  ++it )
		{
		  z_chi2[pointIdx] += (*it)->chi2();
		}
	    }
	}

      // Reset original parameter
      mPar->GetPars()[parIdx] = origParX;
      mPar->GetPars()[parIdx+1] = origParY;


      // Draw graphs
      TString name = "gParScan2D";
      name += i;
      TString title = "Parameter ";
      title += parIdx;
      title += " and ";
      title += parIdx+1;
      title += ";Parameter p_{";
      title += parIdx;
      title += "};Parameter p_{";
      title += parIdx+1;
      title += "};#chi^{2}";
      gParScan2D[i] = new TGraph2D(441,x_param,y_param,z_chi2);
      gParScan2D[i]->SetName(name);
      gParScan2D[i]->SetTitle(title);
      objToBeWritten.push_back(gParScan2D[i]);

      c1->cd();
      gParScan2D[i]->Draw("cont3");
      c1->SetGrid();
      c1->Draw();   
      if( i < nPlots-1 ) ps->NewPage();
    }




  // Clean up
  ps->Close();

  if( mOutputROOT ) WriteToRootFile( objToBeWritten, "ParScan" );

  delete c1;
  delete ps;

  for(int i = 0; i < mPar->GetNumberOfParameters(); i++)
    {
      delete gParScan[i];
      delete hFrame[i];
      delete line[i];
    }
  for(int i = 0; i < nPlots; i++)
    {
      delete gParScan2D[i];
    }
 }




//!  \brief Takes a 2D histogram 'hist' and creates projections along
//!         the x-axis per x bin
//!
//!   Some properties of these projected
//!   distributions are filled, per x-bin, into 8 1D histograms
//!   'hresuslts':
//!    - 0: Mean value
//!    - 1: Standard deviation
//!    - 2: Mean of Gauss fit
//!    - 3: Width of Gauss fit
//!    - 4: Median 
//!    - 5: chi2 / n.d.f.
//!    - 6: Probability of Gauss fit
//!    - 7: Quantiles Q0.9 / (Q0.9 - 1)
//!   'hresuslts' are newly created (take care of deleting them!);
//!   their object names are set to:
//!      "(hist-name)_result(X)",
//!   where (hist-name) = hist->GetName() and (X) is the index of
//!   the above specified property (i.e. 0 for "mean").
//!
//!   Also, the projected distributions of 3 example x-bins:
//!    - 0: hist->GetNbinsX() / 6
//!    - 1: hist->GetNbinsX() / 3
//!    - 2: hist->GetNbinsX() / 2
//!   are filled into 'gaussplots' and the corresponding
//!   Gauss fits are filled into 'gf'. Both 'gaussplots' and 'gf'
//!   are newly created (take care of deleting them!) and their
//!   names are set to:
//!      "(hist-name)_gaussplot(X)",
//!      "(hist-name)_gaussfit(X)"
//!   respectively.
//! ---------------------------------------------------------------
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
      hresults[i]->SetTitle(s + ",  " + mControlQuantityName[i]); 
    }
  s = hist->GetTitle();
  hresults[0]->SetTitle(s + ",  " + mControlQuantityName[0]); 

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
      double mean = htemp->GetMean(); 
      double meanerror = htemp->GetMeanError();
      double width = htemp->GetRMS();
      if(width < 0.1) width = 0.1;
      if(htemp->GetSumOfWeights() <= 0) continue; 
      htemp->Fit("gaus","QNO","", mean - 3 * width,mean + 3 * width);
      TF1 *f = (TF1*)gROOT->GetFunction("gaus")->Clone();
      mean = f->GetParameter(1);
      meanerror = f->GetParError(1);
      width = f->GetParameter(2);
      if(width < 0.05) width = 0.05;
      if( (htemp->Fit(f,"LLQNO","goff",mean - 1.5 * width, mean + 1.5 * width) == 0) && (f->GetProb() > 0.01) ) {
	mean = f->GetParameter(1);
	meanerror = f->GetParError(1);
	width = f->GetParameter(2);
	
	hresults[2]->SetBinContent(i,mean);
	hresults[2]->SetBinError(i,meanerror);
	hresults[3]->SetBinContent(i,width/mean);
	hresults[3]->SetBinError(i, f->GetParError(2)/mean);
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
      hresults[1]->SetBinContent(i,width/mean); 
      hresults[1]->SetBinError(i,htemp->GetRMSError()/mean);
      htemp->GetQuantiles(nq,yq,xq);
      hresults[4]->SetBinContent(i,yq[0]);
      hresults[4]->SetBinError(i,0.0001);
      hresults[7]->SetBinContent(i,yq[1]/yq[0]-1);
      hresults[7]->SetBinError(i,0.0001);
      delete f;
    }
  delete htemp;
}


void TControlPlots::Fit2DRes(const TH2F* hist, TH1F* hresults[8], TH1F* gaussplots[4], TF1* gf[4] ) const
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
      hresults[i]->SetTitle(s + ",  " + mControlQuantityName[i]);
    }
  hresults[0]->SetYTitle("Mittelwert von Jet/genJet"); 
  hresults[2]->SetYTitle("Gauss Peak von Jet/genJet"); 
  hresults[1]->SetYTitle("(Breite von Jet/genJet) / Mittelwert"); 
  hresults[3]->SetYTitle("(Gauss Breite von Jet/genJet) / Peak Value"); 
  s = hist->GetTitle();
  hresults[0]->SetTitle(s + ",  " + mControlQuantityName[0]); 

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
      //if( (htemp->Fit(f,"LLQNO","goff",mean - 2 * width, mean + 2 * width) == 0) && (f->GetProb() > 0.01) ) 
      if( htemp->Fit(f,"LLQNO","goff",mean - 2 * width, mean + 2 * width) == 0) 
	{
	  mean = f->GetParameter(1);
	  meanerror = f->GetParError(1);
	  width = f->GetParameter(2);

	  hresults[2]->SetBinContent(i,mean);
	  hresults[2]->SetBinError(i,meanerror);
	  hresults[3]->SetBinContent(i,width/mean);
	  hresults[3]->SetBinError(i, f->GetParError(2)/mean);
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
      hresults[1]->SetBinContent(i,width/mean); 
      hresults[1]->SetBinError(i,htemp->GetRMSError()/mean);
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
// Write all TObjects in 'obj' to the file 'mOutFile' into the
// directory 'mOutFile:/dir'. If 'mOutFile:/dir' does not exist,
// it is created first.
//---------------------------------------------------------------
void TControlPlots::WriteToRootFile(std::vector<TObject*> obj, std::string dir)
{
  std::string directory = mOutFile->GetName();
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
  gStyle->SetOptStat(0);
  //  gStyle->SetOptStat("neMR");
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
  gStyle->SetPadTopMargin(0.11);
  gStyle->SetPadBottomMargin(0.18);
  gStyle->SetPadLeftMargin(0.25);
  gStyle->SetPadRightMargin(0.04);

  // For the Global title:
  gStyle->SetOptTitle(1);
  gStyle->SetTitleFont(42,"");
  gStyle->SetTitleColor(1);
  gStyle->SetTitleTextColor(1);
  gStyle->SetTitleFillColor(0);
  gStyle->SetTitleFontSize(0.1);
  gStyle->SetTitleAlign(23);
  gStyle->SetTitleX(0.6);
  gStyle->SetTitleH(0.05);
  gStyle->SetTitleBorderSize(0);

  // For the axis titles:
  gStyle->SetTitleColor(1,"XYZ");
  gStyle->SetLabelColor(1,"XYZ");
  // For the axis labels:
  gStyle->SetLabelFont(42,"XYZ");
  gStyle->SetLabelOffset(0.007,"XYZ");
  gStyle->SetLabelSize(0.045,"XYZ");
  // For the axis titles:
  gStyle->SetTitleFont(42,"XYZ");
  gStyle->SetTitleSize(0.06,"XYZ");
  gStyle->SetTitleXOffset(1.2);
  gStyle->SetTitleYOffset(2.0);

  // For the axis:
  gStyle->SetAxisColor(1,"XYZ");
  gStyle->SetStripDecimals(kTRUE);
  gStyle->SetTickLength(0.03,"XYZ");
  gStyle->SetNdivisions(510,"XYZ");
  gStyle->SetPadTickX(1);  // To get tick marks on the opposite side of the frame
  gStyle->SetPadTickY(1);
}



//!  Filling 'bins' with borders of 'nBins' bins between 'first' and 'last'
//!  that are equidistant when viewed in log scale,
//!  so 'bins' must have length nBins+1;
//!  If 'first', 'last' or 'nBins' are not positive, failure is reported.
// -------------------------------------------------------------
bool TControlPlots::EquidistLogBins(double * bins, int nBins, double first, double last) const {
  if( nBins < 1 || first <= 0. || last <= 0. ) return false;
  bins[0]     = first;
  bins[nBins] = last;
  const double firstLog = log10(bins[0]);
  const double lastLog  = log10(bins[nBins]);
  for (int i = 1; i < nBins; ++i) {
    bins[i] = pow(10., firstLog + i*(lastLog-firstLog)/(nBins));
  }

  return true;
}



// -------------------------------------------------------------
TControlPlots::Binning::Binning(const std::vector<double>& binEdgesX, const std::vector<double>& binEdgesY)
{
  assert( binEdgesX.size() > 1 );
  assert( binEdgesY.size() > 1 );

  for(unsigned int i = 0; i < binEdgesX.size(); i++)
    {
      if( i > 0 ) assert( binEdgesX.at(i) > mEdgesX.at(i-1) );
      mEdgesX.push_back(binEdgesX.at(i));
    }
  for(unsigned int i = 0; i < binEdgesY.size(); i++)
    {
      if( i > 0 ) assert( binEdgesY.at(i) > mEdgesY.at(i-1) );
      mEdgesY.push_back(binEdgesY.at(i));
    }
}



// -------------------------------------------------------------
int TControlPlots::Binning::Bin(int ix, int iy) const
{
  int bin = -1;
  if( ix >= 0 && ix < NBinsX() && iy >= 0 && iy < NBinsY() )
    {
      bin = ix + iy*NBinsX();
    }
  return bin;
}



// -------------------------------------------------------------
int TControlPlots::Binning::IX(double x) const
{
  int ix = -1;                     // Underflow
  if( x > mEdgesX.at(NBinsX()) )   // Overflow
    {
      ix = NBinsX();
    }
  else if( x > mEdgesX.at(0) )
    {
      ix = 0;
      while( x > mEdgesX.at(ix+1) ) ix++;
    }
  return ix;
}



// -------------------------------------------------------------
int TControlPlots::Binning::IY(double y) const
{
  int iy = -1;                     // Underflow
  if( y > mEdgesY.at(NBinsY()) )   // Overflow
    {
      iy = NBinsY();
    }
  else if( y > mEdgesY.at(0) )
    {
      iy = 0;
      while( y > mEdgesY.at(iy+1) ) iy++;
    }
  return iy;
}



// -------------------------------------------------------------
void TControlPlots::Binning::Print() const
{
  for(int i = 0; i < NBins(); i++)
    {
      std::cout << i << ":  " << IX(i) << " (" << XLow(i) << ", " << XUp(i) << "),  " << IY(i) << " (" << YLow(i) << ", " << YUp(i) << ")" << std::endl;
    }
}

