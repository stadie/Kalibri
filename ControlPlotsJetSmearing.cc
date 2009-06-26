// $Id: ControlPlotsJetSmearing.cc,v 1.1 2009/06/11 17:34:05 mschrode Exp $

#include "ControlPlotsJetSmearing.h"

#include <cmath>
#include <iostream>

#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TLine.h"
#include "TPostScript.h"
#include "TStyle.h"

#include "SmearData.h"
#include "SmearDiJet.h"
#include "SmearPhotonJet.h"


// --------------------------------------------------
ControlPlotsJetSmearing::ControlPlotsJetSmearing(const std::string& configfile, const std::vector<TData*> * data, TParameters * param)
  : mData(data), mConfig(new ConfigFile(configfile.c_str())), mParam(param),
    mDijetNBins(100), mDijetMin(50), mDijetMax(1000),
    mPhotonJetNBins(100), mPhotonJetMin(50), mPhotonJetMax(1000),
    mRespNBins(150), mRespMin(0.), mRespMax(6.),
    mDir("./controlplots")
{
  SetGStyle();
}



// --------------------------------------------------
void ControlPlotsJetSmearing::PlotDijets() const {
  std::cout << "Creating dijet control plots..." << std::endl;

  // Create histograms
  TH1F * hDijetPtTrue       = new TH1F("hDijetPtTrue",";p^{true}_{T} (GeV)",mDijetNBins,mDijetMin,mDijetMax);
  TH1F * hDijetPtMeas       = new TH1F("hDijetPtMeas",";p^{jet}_{T} (GeV)",mDijetNBins,mDijetMin,mDijetMax);
  TH1F * hDijetDeltaPtMeas  = new TH1F("hDijetDeltaPtMeas",
				       ";(p^{jet}_{T,0} - p^{jet}_{T,1}) / (p^{jet}_{T,0} + p^{jet}_{T,1})",
				       mDijetNBins,-0.5,0.5);
  TH1F * hDijetDeltaPtTrue  = new TH1F("hDijetDeltaPtTrue",
				       ";(p^{true}_{T,0} - p^{true}_{T,1}) / (p^{true}_{T,0} + p^{true}_{T,1})",
				       mDijetNBins,-0.5,0.5);
  TH1F * hDijetRelPt        = new TH1F("hDijetRelPt",";p^{jet}_{T} / p^{true}_{T}",mDijetNBins,0,2);


  // Fill histograms
  for(DataIt datait = mData->begin(); datait != mData->end(); datait++) {  // Loop over Data
    if( (*datait)->GetType() == TypeSmearDiJet )  {                        // Select DiJet events
      SmearDiJet * dijet = static_cast<SmearDiJet*>(*datait);  
      for(int i = 0; i < 2; i++) {        // Loop over both jets
	TJet * jet = static_cast<TJet*>(dijet->GetMess());
	if( i == 1 ) jet = static_cast<TJet*>(dijet->GetSecondMess());
	if( jet ) {
	  double t = jet->genPt;          // truth
	  double m = jet->pt;             // measurement
	    
	  hDijetPtTrue->Fill(t);
	  hDijetPtMeas->Fill(m);
	  hDijetRelPt->Fill(m/t);		    
	}
      }

      TJet * jet1 = static_cast<TJet*>(dijet->GetMess());
      TJet * jet2 = static_cast<TJet*>(dijet->GetSecondMess());

      if( jet1 && jet2 ) {
	double dif = jet1->genPt - jet2->genPt;
	double sum = jet1->genPt + jet2->genPt;
	hDijetDeltaPtTrue->Fill( dif/sum );

	dif = jet1->pt - jet2->pt;
	sum = jet1->pt + jet2->pt;
	hDijetDeltaPtMeas->Fill( dif/sum );
      }
    }
  } // End loop over Data


    // Write histos to ps file
  std::vector<TH1F*> hists;
  hists.push_back(hDijetPtTrue);
  hists.push_back(hDijetPtMeas);
  hists.push_back(hDijetDeltaPtMeas);
  hists.push_back(hDijetDeltaPtTrue);
  hists.push_back(hDijetRelPt);

  TPostScript * const ps = new TPostScript((mDir+"/Dijets.ps").c_str(),111);
  TCanvas *c1 = new TCanvas("c1","Dijets",0,0,600,600);
  c1->cd();

  std::vector<TH1F*>::iterator it = hists.begin();
  for( ; it != hists.end(); it++) {
    (*it)->Draw();
    c1->Draw();
    ps->NewPage();
  }
  ps->Close();


  // Clean up
  for( ; it != hists.end(); it++) {
    delete *it;
  }
  hists.clear();
  delete ps;
  delete c1;

  std::cout << "Done" << std::endl;
}



// --------------------------------------------------
void ControlPlotsJetSmearing::PlotResponse() const
{
  std::cout << "Creating response control plots..." << std::endl;

  // Create histograms of response
  TH1D* hRespMeas = new TH1D("hRespMeas",";R = p^{jet}_{T} / p^{true}_{T};1/(Nw)  dN / dR",mRespNBins,mRespMin,mRespMax);
  hRespMeas->Sumw2();

  TH1D* hRespFit = new TH1D("hRespFit",";R = p^{jet}_{T} / p^{true}_{T};1/(Nw)  dN / dR",5*mRespNBins,mRespMin,mRespMax);
  hRespFit->SetLineColor(2);
  hRespFit->SetLineWidth(2);
  hRespFit->Sumw2();

  TH1D* hRespFitStep = new TH1D("hRespFitStep",";R = p^{jet}_{T} / p^{true}_{T};1 / (Nw)  dN / dR",
				mConfig->read<int>("Response pdf nsteps",10),
				mConfig->read<double>("Response pdf min",0.),
				mConfig->read<double>("Response pdf max",1.8));
  hRespFitStep->Sumw2();
  hRespFitStep->SetLineColor(9);
  hRespFitStep->SetLineWidth(2);

  TH1D* hRespFitGaus = new TH1D("hRespFitGaus",";R = p^{jet}_{T} / p^{true}_{T};1 / (Nw)  dN / dR",5*mRespNBins,mRespMin,mRespMax);
  hRespFitGaus->Sumw2();
  hRespFitGaus->SetLineColor(8);
  hRespFitGaus->SetLineWidth(2);

  TH1D* hRespFitSum = new TH1D("hRespFitSum",";R = p^{jet}_{T} / p^{true}_{T};1 / (Nw)  dN / dR",5*mRespNBins,mRespMin,mRespMax);
  hRespFitSum->Sumw2();
  hRespFitSum->SetLineColor(1);
  hRespFitSum->SetLineWidth(2);


  // Create histograms of genjet pdf
  double mingpt = mConfig->read<double>("DiJet integration min truth",0);
  double maxgpt = mConfig->read<double>("DiJet integration max truth",1);

  TH1D * hDijetGenJetPt     = new TH1D("hDijetGenJetPt",
				       ";p^{true}_{T} (GeV);1 / (Nw)  dN / dp^{true}_{T}",
				       25,0.9*mingpt,1.1*maxgpt);
  hDijetGenJetPt->GetXaxis()->SetNdivisions(505);
  hDijetGenJetPt->Sumw2();

  TH1D * hDijetTruthPDF     = new TH1D("hDijetTruthPDF",
				       ";p^{true}_{T} (GeV);1 / (Nw)  dN / dp^{true}_{T}",
				       5*mRespNBins,0.9*mingpt,1.1*maxgpt);
  hDijetTruthPDF->SetLineColor(2);
  hDijetTruthPDF->SetLineWidth(2);

  // Fill histogram of measured response
  for(DataIt datait = mData->begin(); datait != mData->end(); datait++) {
    // Select DiJet events
    if( (*datait)->GetType() == TypeSmearDiJet )  {
      SmearDiJet * dijet = static_cast<SmearDiJet*>(*datait);  
      for(int i = 0; i < 2; i++) {        // Loop over both jets
	TJet * jet = static_cast<TJet*>(dijet->GetMess());
	if( i == 1 ) jet = static_cast<TJet*>(dijet->GetSecondMess());
	hDijetGenJetPt->Fill(jet->genPt);
	hRespMeas->Fill(jet->pt / jet->genPt);
      }
    }

    // Select PhotonJet events
    if( (*datait)->GetType() == TypeSmearPhotonJet ) {
      SmearPhotonJet * pjet = static_cast<SmearPhotonJet*>(*datait);  
      double photonpt = pjet->GetTruth();
      double jetpt    = pjet->GetMess()->pt;
      hRespMeas->Fill( jetpt / photonpt );		    
    }
  } // End of loop over data

  if( hRespMeas->Integral("width") ) hRespMeas->Scale(1./hRespMeas->Integral("width"));
  if( hDijetGenJetPt->Integral("width") ) hDijetGenJetPt->Scale(1./hDijetGenJetPt->Integral("width"));


  // Fill histogram of fitted response
  SmearData * smeardata = dynamic_cast<SmearData*>(mData->front());
  if( smeardata ) {
    // Interpolated fit function
    for(int bin = 1; bin <= hRespFit->GetNbinsX(); bin++) {
      double r = hRespFit->GetBinCenter(bin);
      hRespFit->SetBinContent(bin,smeardata->RespPDF(r));
    }

    // In case of step + gauss parametrizations
    std::string param = mConfig->read<std::string>("Parametrization Class","");
    if( param == "SmearParametrizationStepGauss" || param == "SmearParametrizationStepGaussInter" ) {
      std::vector<double> scale = bag_of<double>(mConfig->read<string>("Jet parameter scales",""));

      // Step part of fit function
      for(int bin = 1; bin <= hRespFitStep->GetNbinsX(); bin++) {
	double val  = scale.at(bin+2)*(smeardata->GetRespPar()[bin+2]);
	hRespFitStep->SetBinContent(bin,val);
      }
      hRespFitStep->Scale(1./hRespFitStep->Integral("width"));
      hRespFitStep->Scale(1. - scale.at(0)*(smeardata->GetRespPar()[0]));
      
      // Gauss part of fit function
      for(int bin = 1; bin <= hRespFitGaus->GetNbinsX(); bin++) {
	double c     = scale.at(0)*(smeardata->GetRespPar()[0]);
	double mu    = scale.at(1)*(smeardata->GetRespPar()[1]);
	double sigma = scale.at(2)*(smeardata->GetRespPar()[2]);
	double r     = hRespFitGaus->GetBinCenter(bin);
	double val   = c * exp( -pow((mu-r)/sigma,2) / 2. ) / sqrt(2.*M_PI) / sigma;
	hRespFitGaus->SetBinContent(bin,val);
      }
      
      // Sum
      for(int binGaus = 1; binGaus <= hRespFitGaus->GetNbinsX(); binGaus++) {
	int    binStep = hRespFitStep->FindBin(hRespFitGaus->GetBinCenter(binGaus));
	double val     = hRespFitStep->GetBinContent(binStep) + hRespFitGaus->GetBinContent(binGaus);
	hRespFitSum->SetBinContent(binGaus,val);
      }
    }

    // In case of two gauss parametrizations
    else if( param == "SmearParametrizationTwoGauss" ) {
    
      std::vector<double> scale = bag_of<double>(mConfig->read<string>("Jet parameter scales",""));
      
      // Central Gauss of fit function
      for(int bin = 1; bin <= hRespFitGaus->GetNbinsX(); bin++) {
	double c     = scale.at(0)*(smeardata->GetRespPar()[0]);
	double mu    = scale.at(1)*(smeardata->GetRespPar()[1]);
	double sigma = scale.at(2)*(smeardata->GetRespPar()[2]);
	double r     = hRespFitGaus->GetBinCenter(bin);
	double val   = c * exp( -pow((mu-r)/sigma,2) / 2. ) / sqrt(2.*M_PI) / sigma;
	hRespFitGaus->SetBinContent(bin,val);
      }
      
      // Tail Gauss of fit function
      for(int bin = 1; bin <= hRespFitSum->GetNbinsX(); bin++) {
	double c     = scale.at(0)*(smeardata->GetRespPar()[0]);
	double mu    = scale.at(3)*(smeardata->GetRespPar()[3]);
	double sigma = scale.at(4)*(smeardata->GetRespPar()[4]);
	double r     = hRespFitGaus->GetBinCenter(bin);
	double val   = (1.-c) * exp( -pow((mu-r)/sigma,2) / 2. ) / sqrt(2.*M_PI) / sigma;
	hRespFitSum->SetBinContent(bin,val);
      }
    }
  }

  // Fill histogram of assumed dijet truth pdf
  // and fit with 1/x^n function
  double n      = 6.5;
  DataIt datait = mData->begin();
  while( (*datait)->GetType() != TypeSmearDiJet  &&  datait != mData->end() ) datait++;
  if( datait != mData->end() ) {
    SmearDiJet * dijet = static_cast<SmearDiJet*>(*datait);  
    for(int bin = 1; bin <= hDijetTruthPDF->GetNbinsX(); bin++) {
      double t = hDijetTruthPDF->GetBinCenter(bin);
      hDijetTruthPDF->SetBinContent(bin,dijet->TruthPDF(t));
    }

    n            = (dijet->GetTruthPar()[0]);
    double k     = n - 3.;
    double m     = n - 1.;
    double norm1 = k / ( pow(mingpt,-k) - pow(maxgpt,-k) );
    double norm2 = m / ( pow(mingpt,-m) - pow(maxgpt,-m) );
    hDijetTruthPDF->Scale(norm2/norm1);

    //   double norm = ( pow(maxgpt,3) + pow(mingpt,3) ) / 3.;
    //   hDijetTruthPDF->Scale(norm/(maxgpt-mingpt));
  }

  // Find populated x-axis range
  int maxBin = 0;
  int minBin = 1;
  for(int bin = 1; bin <= hRespMeas->GetNbinsX(); bin++) {
    if( hRespMeas->GetBinContent(bin) > 0  ) maxBin = bin;
    if( minBin > maxBin ) minBin = bin;
  }
  if( maxBin < hRespMeas->GetNbinsX() ) maxBin++;
  if( minBin > 1 )                      minBin--;
  hRespMeas->GetXaxis()->SetRange(minBin,maxBin);


  // Write histos to eps file
  TCanvas *c1 = new TCanvas("c1","Jet Response",0,0,600,600);
  c1->cd();
  hRespMeas->Draw();
  hRespFit->Draw("Lsame");
  c1->SetLogy();
  c1->SaveAs((mDir+"/JetResponse.eps").c_str());

  c1->cd();
  hRespMeas->Draw();
  std::string param = mConfig->read<std::string>("Parametrization Class","");
  if( param == "SmearParametrizationStepGauss" || param == "SmearParametrizationStepGaussInter" ) {
    hRespFitStep->Draw("same");
    hRespFitGaus->Draw("same");
    hRespFitSum->Draw("same");
  } else if( param == "SmearParametrizationTwoGauss" ) {
    hRespFitGaus->Draw("same");
    hRespFitSum->Draw("same");
  }
  hRespFit->Draw("Lsame");
  c1->SetLogy();
  c1->SaveAs((mDir+"/JetResponseDetail.eps").c_str());

  c1->cd();
  hDijetGenJetPt->Draw();
  hDijetTruthPDF->Draw("SAME");

  TLegend *leg = new TLegend(0.4,0.75,0.89,0.89);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextFont(42);
  char entry[50];
  sprintf(entry,"Fit #propto 1 / (p^{true}_{T})^{%.1f}",n);
  leg->AddEntry(hDijetTruthPDF,entry,"L");
  leg->Draw("same");

  c1->SetLogy();
  c1->SaveAs((mDir+"/DiJetTruth.eps").c_str());

  // Write histos to root file
  TFile rootfile((mDir+"/JetResponse.root").c_str(),"RECREATE");
  rootfile.WriteTObject(hRespMeas);
  rootfile.WriteTObject(hRespFit);
  rootfile.WriteTObject(hRespFitStep);
  rootfile.WriteTObject(hRespFitGaus);
  rootfile.WriteTObject(hRespFitSum);
  rootfile.WriteTObject(hDijetGenJetPt);
  rootfile.WriteTObject(hDijetTruthPDF);
  rootfile.Close();

  // Clean up
  delete c1;
  delete hRespMeas;
  delete hRespFit;
  delete hRespFitStep;
  delete hRespFitGaus;
  delete hRespFitSum;
  delete hDijetGenJetPt;
  delete hDijetTruthPDF;
  delete leg;

  std::cout << "Done" << std::endl;
}



// --------------------------------------------------
void ControlPlotsJetSmearing::PlotParameterScan(const std::vector<unsigned int>& pars) const {
  std::cout << "Creating parameter scan control plots..." << std::endl;

  // Create one histogram of likelihood per parameter
  std::vector<TH1D*> hLikelihood;
  for(int i = 0; i < mParam->GetNumberOfParameters(); i++) {
    char name[50];
    sprintf(name,"hNLogL%i",i);
    char title[50];
    sprintf(title,";Parameter %i;-ln(L)",i);
    double min = mParam->GetPars()[i] - 0.5;
    if( min < 0. ) min = 0.;
    double max = mParam->GetPars()[i] + 0.5;
    TH1D * h = new TH1D(name,title,10,min,max);
    h->SetMarkerStyle(24);
    hLikelihood.push_back(h);
  }

  // Do parameter scan
  for(unsigned int i = 0; i < pars.size(); i++) {
    int idx = pars.at(i);
    if( idx < mParam->GetNumberOfParameters() ) {
      std::cout << "  Par " << idx << "... " << std::flush;
      double oldpar = mParam->GetPars()[idx];
      TH1D * h = hLikelihood.at(idx);
      for(int bin = 1; bin <= h->GetNbinsX(); bin++) {
	mParam->GetPars()[idx] = h->GetBinCenter(bin);
	double fsum = 0.;
	for(DataIt datait = mData->begin(); datait != mData->end(); datait++) {
	  fsum += (*datait)->chi2();
	}
	h->SetBinContent(bin,fsum);
      }
      mParam->GetPars()[idx] = oldpar;
      std::cout << "ok" << std::endl;
    }
  }

  // Plot histograms
  TPostScript * const ps = new TPostScript((mDir+"/ParameterScan.ps").c_str(),111);
  TCanvas *c1 = new TCanvas("c1","Parameter Scan",0,0,600,600);
  ps->NewPage();
  c1->cd();
  for(unsigned int i = 0; i < pars.size(); i++) {
    int idx = pars.at(i);
    if( idx < mParam->GetNumberOfParameters() ) {
      TH1D * h = hLikelihood.at(idx);
      h->Draw("P");
      TLine * line = new TLine(mParam->GetPars()[idx],h->GetMinimum(),
			       mParam->GetPars()[idx],h->GetMaximum());
      line->SetLineStyle(2);
      line->SetLineWidth(2);
      line->Draw("same");
      c1->Draw();
      ps->NewPage();
      delete line;
    }
  }
  ps->Close();

  // Clean up
  for(std::vector<TH1D*>::iterator it = hLikelihood.begin(); it != hLikelihood.end(); it++) {
    delete *it;
  }
  hLikelihood.clear();
  delete c1;
  delete ps;

  std::cout << "Done" << std::endl;
}



// --------------------------------------------------
void ControlPlotsJetSmearing::SetGStyle() const
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

  //  Margins
  // -------------------------------------------
  gStyle->SetPadTopMargin(0.05);
  gStyle->SetPadBottomMargin(0.18);
  gStyle->SetPadLeftMargin(0.19);
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

  // For the axis labels:
  //  For the axis labels and titles
  // -------------------------------------------
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
  gStyle->SetTitleYOffset(1.5);

  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
}

