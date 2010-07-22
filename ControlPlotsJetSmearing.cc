// $Id: ControlPlotsJetSmearing.cc,v 1.17 2010/07/22 13:58:30 mschrode Exp $

#include "ControlPlotsJetSmearing.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <sstream>

#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraphAsymmErrors.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TLine.h"
#include "TPad.h"
#include "TPaveText.h"
#include "TPostScript.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TVector2.h"
#include "TH1.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TH2.h"
#include "TError.h"
#include "TH2F.h"
#include "TMath.h"
#include "TPaveText.h"
#include "TRandom3.h"
#include "TROOT.h"


#include "SmearData.h"
#include "SmearDiJet.h"
#include "SmearPhotonJet.h"
#include "Jet.h"



//!  \brief Constructor
//!
//!  By default, all controlplots are written to the
//!  directory \p ./controlplots/
//!
//!  \param configfile Name of the configuration file
//!  \param data The data
//!  \param param The parametrization
// --------------------------------------------------
ControlPlotsJetSmearing::ControlPlotsJetSmearing(const std::string& configfile, const std::vector<Event*> * data, TParameters * param, const std::string &outDir)
  : data_(data),
    config_(new ConfigFile(configfile.c_str())),
    param_(param),
    respNBins_(60),
    respMin_(0.),
    respMax_(2.),
    dir_(outDir)
{
  // Do not print ROOT message if eps file has been created
  gErrorIgnoreLevel = 1001;

  // Override possible existing root file
  TFile rootfile((dir_+"/jsResponse.root").c_str(),"RECREATE");
  rootfile.Close();
  setGStyle();
  rand_ = new TRandom3(0);

  parClass_ = config_->read<std::string>("Parametrization Class","");
  startParJet_ = bag_of<double>(config_->read<string>("jet start values",""));
  scale_ = bag_of<double>(config_->read<string>("jet parameter scales",""));
  truthPar_ = bag_of<double>(config_->read<string>("plots true resolution parameters",""));

  // Create pt binning
  std::string binning = config_->read<std::string>("plots pt binning","");
  ptBinningVar_ = "ptGen";
  if( binning.find("ptDijet") != std::string::npos ) ptBinningVar_ = "ptDijet";

  ptBinEdges_ = std::vector<double>(2,0.);
  if( binning.find("binning") != std::string::npos ) {
    ptBinEdges_.clear();
    ptBinEdges_ = bag_of<double>(config_->read<std::string>("plots pt bin edges","100 500"));
  } else if( binning.find("cuts") != std::string::npos ) {
    if( ptBinningVar_ == "ptGen" ) {
      ptBinEdges_.clear();
      ptBinEdges_.push_back(config_->read<double>("Et genJet min",100));
      ptBinEdges_.push_back(config_->read<double>("Et genJet max",500));
    } else if( ptBinningVar_ == "pt" ) {
      ptBinEdges_.clear();
      ptBinEdges_.push_back(config_->read<double>("Et min cut on jet",100));
      ptBinEdges_.push_back(config_->read<double>("Et max cut on jet",500));
    } else if( ptBinningVar_ == "ptDijet" ) {
      ptBinEdges_.clear();
      ptBinEdges_.push_back(config_->read<double>("Et min cut on dijet",100));
      ptBinEdges_.push_back(config_->read<double>("Et max cut on dijet",500));
    }
  }
  ptBinCenters_ = std::vector<double>(nPtBins());
  for(int i = 0; i < nPtBins(); i++) {
    ptBinCenters_[i] = 0.5 * ( ptBinEdges_[i] + ptBinEdges_[i+1] );
  }

  minJetPt_ = config_->read<double>("Et min cut on jet",0.);
  maxJetPt_ = config_->read<double>("Et max cut on jet",1.);

  outNamePrefix_ = dir_+"/"+config_->read<std::string>("plots name prefix","JS");

  saveAsEps_ = config_->read<bool>("plots save as eps",false);
}

ControlPlotsJetSmearing::~ControlPlotsJetSmearing() {
  ptBinEdges_.clear();
  ptBinCenters_.clear();
  delete rand_;
}



// --------------------------------------------------
void ControlPlotsJetSmearing::makePlots() const {
  if( config_->read<bool>("create 3rd jet plots",false) )
    plot3rdJet();
  if( config_->read<bool>("create dijet plots",false) )
    plotDijets();
  if( config_->read<bool>("create logP plots",false) )
    plotLogP();
  if( config_->read<bool>("create mean response and resolution plots",false) )
    plotMeanResponseAndResolution();
  if( config_->read<bool>("create parameter error plots",false) )
    plotParameters();
  if( config_->read<bool>("create parameter scan plots",false) )
    plotParameterScan();
  if( config_->read<bool>("create response plots",false) )
    plotResponse();
  if( config_->read<bool>("create simulated asymmetry plots",false) )
    plotAsymmetrySimulation();
}



//!  \brief Draw response control plots for events
//!         of type \p SmearDiJet
// --------------------------------------------------
void ControlPlotsJetSmearing::plotResponse() const
{
  std::cout << "Creating response control plots\n";

  // --- Create histograms of response and spectrum ---------------------
  std::string param = config_->read<std::string>("Parametrization Class","");
  double tMin = config_->read<double>("DiJet integration min",0.);
  double tMax = config_->read<double>("DiJet integration max",1.);
  double rMin = config_->read<double>("Response pdf min",0.);
  double rMax = config_->read<double>("Response pdf max",2.);

  std::vector<TH1*> hRespMeasAbs(nPtBins());   // The response ptJet / ptGen absolute entries
  std::vector<TH1*> hRespMeas(nPtBins());      // The response ptJet / ptGen
  std::vector<TH1*> hRespMeasJet1(nPtBins());      // The response ptJet1 / ptGen
  std::vector<TH1*> hRespMeasJet2(nPtBins());      // The response ptJet2 / ptGen
  std::vector<TH1*> hRespMCPtHat(nPtBins());      // The response ptJet / ptGen
  std::vector<TH1*> hRespFitStart(nPtBins());  // The response pdf with start values
  std::vector<TH1*> hRespFit(nPtBins());       // The fitted response pdf
  std::vector<TH1*> hRespFitErrStat(nPtBins());// The fitted response pdf with fitted errors
  std::vector<TH1*> hRespFitStep(nPtBins());   // Step function part of the response pdf
  std::vector<TH1*> hRespFitGaus(nPtBins());   // Gauss part of the response pdf
  std::vector<TH1*> hRespFitSum(nPtBins());    // Sum of step and Gauss part
  std::vector<TH1*> hRespFitBins(nPtBins());       // The fitted response pdf in bins as MC truth
  std::vector<TH1*> hRespRatio(nPtBins());
  std::vector<TH1*> hRespRatioJet1(nPtBins());
  std::vector<TH1*> hRespRatioJet2(nPtBins());
  std::vector<TH1*> hPtGenAbsBins(nPtBins());
  std::vector<TH1*> hPtGenAsym(nPtBins());
  std::vector<TH1*> hPtAsym(nPtBins());
  std::vector<TH1*> hPtAsymBiased(nPtBins());
  std::vector<TH1*> hFitPtAsym(nPtBins());

  TH1 * hRespRatioFrame = 0;
  TH1 * hRespRatioFrameJet1 = 0;
  TH1 * hRespRatioFrameJet2 = 0;
  TH1 * hTruthPDF = 0;      // Truth pdf
  TH1 * hTruthPDFErrStat = 0;      // Truth pdf
  TH1 * hPtGenAbs = 0;         // PtGen spectrum
  TH1 * hPtGen = 0;         // PtGen spectrum
  TH1 * hPtGenJet1 = 0;         // PtGen spectrum
  TH1 * hPtGenJet2 = 0;         // PtGen spectrum
  TH1 * hPtHat = 0;         // PtHat spectrum
  TH1 * hPtDijet = 0;       // Dijet spectrum
  TH1 * hGaussWidth = 0;
  TH1 * hGaussWidthErr = 0;
  TH1 * hGaussWidthTruth = 0;
  TH1 * hGaussWidthMC = 0;
  TH1 * hGaussWidthPseudo = 0;
  TH1 * hGaussWidthRatio = 0;
  TH1 * hGaussWidthRatioErr = 0;
  TH1 * hGaussWidthRatioMC = 0;		  
  TH1 * hGaussWidthRatioPseudo = 0;
  TH1 * hJESFrame = 0;
  TH1 * hJES = 0;
  TH1 * hJESJet1 = 0;
  TH1 * hJESJet2 = 0;

  for(int ptBin = 0; ptBin < nPtBins(); ptBin++) {
    std::string name = "hRespMeasAbs_" + toString(ptBin);
    hRespMeasAbs[ptBin] = new TH1D(name.c_str(),";R = p^{jet}_{T} / p^{gen}_{T};dN / dR",
				   respNBins_,respMin_,respMax_);
    hRespMeasAbs[ptBin]->SetLineWidth(2);

    name = "hRespMeas_" + toString(ptBin);
    hRespMeas[ptBin] = new TH1D(name.c_str(),";R = p^{jet}_{T} / p^{gen}_{T};1 / N  dN / dR",
				respNBins_,respMin_,respMax_);
    hRespMeas[ptBin]->Sumw2();
    hRespMeas[ptBin]->SetMarkerStyle(20);
    hRespMeas[ptBin]->SetLineWidth(2);

    name = "hRespMeasJet1_" + toString(ptBin);
    hRespMeasJet1[ptBin] = static_cast<TH1D*>(hRespMeas[ptBin]->Clone(name.c_str()));
    hRespMeasJet1[ptBin]->GetXaxis()->SetTitle("R_{1} = p^{jet}_{T,1} / p^{gen}_{T,1}");

    name = "hRespMeasJet2_" + toString(ptBin);
    hRespMeasJet2[ptBin] = static_cast<TH1D*>(hRespMeas[ptBin]->Clone(name.c_str()));
    hRespMeasJet2[ptBin]->GetXaxis()->SetTitle("R_{2} = p^{jet}_{T,2} / p^{gen}_{T,2}");

    name = "hRespMCPtHat_" + toString(ptBin);
    hRespMCPtHat[ptBin] = static_cast<TH1D*>(hRespMeas[ptBin]->Clone(name.c_str()));
    hRespMCPtHat[ptBin]->SetTitle(";R = p^{jet}_{T} / #hat{p}_{T};1 / N  dN / dR");

    name = "hRespFit_" + toString(ptBin);
    hRespFit[ptBin] = new TH1D(name.c_str(),";R = p^{jet}_{T} / p^{true}_{T};1 / N  dN / dR",
			       5*respNBins_,respMin_,respMax_);
    hRespFit[ptBin]->SetLineColor(2);
    hRespFit[ptBin]->SetLineWidth(2);

    name = "hRespFitErrStat_" + toString(ptBin);
    hRespFitErrStat[ptBin] = static_cast<TH1D*>(hRespFit[ptBin]->Clone(name.c_str()));
    hRespFitErrStat[ptBin]->SetFillColor(45);

    name = "hRespFitStart_" + toString(ptBin);
    hRespFitStart[ptBin] = static_cast<TH1D*>(hRespFit[ptBin]->Clone(name.c_str()));
    hRespFitStart[ptBin]->SetLineStyle(2);

    name = "hRespFitStep_" + toString(ptBin);
    hRespFitStep[ptBin] = new TH1D(name.c_str(),";R = p^{jet}_{T} / p^{true}_{T};1 / N  dN / dR",
				   config_->read<int>("Response pdf nsteps",10),
				   config_->read<double>("Response pdf min",0.),
				   config_->read<double>("Response pdf max",1.8));
    hRespFitStep[ptBin]->Sumw2();
    hRespFitStep[ptBin]->SetLineColor(9);
    hRespFitStep[ptBin]->SetLineWidth(2);

    name = "hRespFitGaus_" + toString(ptBin);
    hRespFitGaus[ptBin] = static_cast<TH1D*>(hRespFit[ptBin]->Clone(name.c_str()));
    hRespFitGaus[ptBin]->SetLineColor(8);

    name = "hRespFitSum_" + toString(ptBin);
    hRespFitSum[ptBin] = static_cast<TH1D*>(hRespFit[ptBin]->Clone(name.c_str()));
    hRespFitSum[ptBin]->SetLineColor(1);

    name = "hRespFitBins_" + toString(ptBin);
    hRespFitBins[ptBin] = static_cast<TH1D*>(hRespMeas[ptBin]->Clone(name.c_str()));
    hRespFitBins[ptBin]->SetLineColor(2);
    hRespFitBins[ptBin]->SetLineWidth(2);

    name = "hRespRatio_" + toString(ptBin);
    hRespRatio[ptBin] = static_cast<TH1D*>(hRespMeas[ptBin]->Clone(name.c_str()));
    hRespRatio[ptBin]->SetYTitle("Fit / MC truth");

    name = "hRespRatioJet1_" + toString(ptBin);
    hRespRatioJet1[ptBin] = static_cast<TH1D*>(hRespRatio[ptBin]->Clone(name.c_str()));
    hRespRatioJet1[ptBin]->SetXTitle(";R_{1} = p^{jet}_{T,1} / p^{true}_{T,1}");

    name = "hRespRatioJet2_" + toString(ptBin);
    hRespRatioJet2[ptBin] = static_cast<TH1D*>(hRespRatio[ptBin]->Clone(name.c_str()));
    hRespRatioJet2[ptBin]->SetXTitle(";R_{2} = p^{jet}_{T,2} / p^{true}_{T,2}");

    name = "hPtGenAbs_" + toString(ptBin);
    hPtGenAbsBins[ptBin] = new TH1D(name.c_str(),";p^{gen}_{T} (GeV);dN / dp^{gen}_{T}  1 / (GeV)"
				    ,25,0.8*ptBinEdges_[ptBin],1.1*ptBinEdges_[ptBin+1]);
    hPtGenAbsBins[ptBin]->GetXaxis()->SetNdivisions(505);
    hPtGenAbsBins[ptBin]->SetLineWidth(2);

    name = "hPtGenAsym_" + toString(ptBin);
    hPtGenAsym[ptBin] = new TH1D(name.c_str(),
				 ";p^{gen}_{T} asymmetry;",
				 30,-0.3,0.3);
    hPtGenAsym[ptBin]->Sumw2();
    hPtGenAsym[ptBin]->SetLineWidth(2);

    name = "hPtAsym_" + toString(ptBin);
    hPtAsym[ptBin] = new TH1D(name.c_str(),";p_{T} asymmetry;",75,-1.,1.);
    hPtAsym[ptBin]->Sumw2();
    hPtAsym[ptBin]->SetMarkerStyle(20);
    hPtAsym[ptBin]->SetLineWidth(2);

    name = "hPtAsymBiased_" + toString(ptBin);
    hPtAsymBiased[ptBin] = new TH1D(name.c_str(),";Biased p_{T} asymmetry;",32,0.,1.);
    hPtAsymBiased[ptBin]->Sumw2();
    hPtAsymBiased[ptBin]->SetMarkerStyle(20);
    hPtAsymBiased[ptBin]->SetLineWidth(2);

    name = "hFitPtAsym_" + toString(ptBin);
    hFitPtAsym[ptBin] = new TH1D(name.c_str(),
				 ";p_{T} asymmetry;",
				 5*respNBins_,-0.4,0.4);
    hFitPtAsym[ptBin]->Sumw2();
    hFitPtAsym[ptBin]->SetLineWidth(2);
    hFitPtAsym[ptBin]->SetLineColor(2);
  }

  hRespRatioFrame = static_cast<TH1D*>(hRespRatio[0]->Clone("hRespRatioFrame"));
  for(int bin = 1; bin <= hRespRatioFrame->GetNbinsX(); ++bin) {
    hRespRatioFrame->SetBinContent(bin,1.);
  }
  hRespRatioFrame->SetLineStyle(2);
  hRespRatioFrame->GetYaxis()->SetRangeUser(0,2.8);
  hRespRatioFrame->GetXaxis()->SetRangeUser(rMin,rMax);
  hRespRatioFrameJet1 = static_cast<TH1D*>(hRespRatioFrame->Clone("hRespRatioFrameJet1"));
  hRespRatioFrameJet1->SetXTitle("R_{1} = p^{jet}_{T,1} / p^{true}_{T,1}");
  hRespRatioFrameJet2 = static_cast<TH1D*>(hRespRatioFrame->Clone("hRespRatioFrameJet2"));
  hRespRatioFrameJet2->SetXTitle("R_{2} = p^{jet}_{T,2} / p^{true}_{T,2}");

  hPtGenAbs = new TH1D("hPtGenAbs",";p^{gen}_{T} (GeV);dN / dp^{gen}_{T}  1 / (GeV)",
		       60,0.8*ptBinEdges_.front(),1.1*ptBinEdges_.back());
  hPtGenAbs->SetMarkerStyle(20);
  hPtGenAbs->GetXaxis()->SetNdivisions(505);
  hPtGenAbs->SetLineWidth(2);

  hPtGen = static_cast<TH1D*>(hPtGenAbs->Clone("hPtGen"));
  hPtGen->SetTitle(";p^{gen}_{T} (GeV);1 / N  dN / dp^{gen}_{T}  1 / (GeV)");
  hPtGen->Sumw2();

  hPtGenJet1 = static_cast<TH1D*>(hPtGen->Clone("hPtGenJet1"));
  hPtGenJet1->SetTitle(";p^{gen}_{T,1} (GeV);1 / N  dN / dp^{gen}_{T,1}  1 / (GeV)");

  hPtGenJet2 = static_cast<TH1D*>(hPtGen->Clone("hPtGenJet2"));
  hPtGenJet2->SetTitle(";p^{gen}_{T,2} (GeV);1 / N  dN / dp^{gen}_{T,2}  1 / (GeV)");
  
  hPtHat = static_cast<TH1D*>(hPtGen->Clone("hPtHat"));
  hPtHat->SetTitle(";#hat{p}_{T} (GeV);1 / N  dN / d#hat{p}_{T}  1 / (GeV)");

  hPtDijet = static_cast<TH1D*>(hPtGen->Clone("hPtDijet"));
  hPtDijet->SetTitle(";p^{dijet}_{T} (GeV);1 / N  dN / dp^{dijet}_{T}  1 / (GeV)");

  hTruthPDF = new TH1D("hTruthPDF",";p^{true}_{T} (GeV);1 / N  dN / dp^{true}_{T}  1 /  (GeV)",
		       5*hPtGen->GetNbinsX(),tMin,tMax);
  hTruthPDF->SetLineColor(2);
  hTruthPDF->SetLineWidth(2);
  hTruthPDF->Sumw2();

  hTruthPDFErrStat = static_cast<TH1D*>(hTruthPDF->Clone("hTruthPDFErrStat"));
  hTruthPDFErrStat->SetFillColor(45);

  hGaussWidth = new TH1D("hGaussWidth",";p_{T} (GeV);#sigma(p_{T}) / p_{T}",
			 500,ptBinEdges_.front(),ptBinEdges_.back());
  hGaussWidth->SetNdivisions(505);
  hGaussWidth->SetMarkerStyle(1);
  hGaussWidth->SetLineWidth(2);
  hGaussWidth->SetLineColor(2);
  hGaussWidthErr = static_cast<TH1D*>(hGaussWidth->Clone("hGaussWidthErr"));
  hGaussWidthErr->SetMarkerStyle(1);
  hGaussWidthErr->SetLineWidth(0);
  hGaussWidthErr->SetLineColor(43);
  hGaussWidthErr->SetFillColor(43);
  hGaussWidthTruth = static_cast<TH1D*>(hGaussWidth->Clone("hGaussWidthTruth"));
  hGaussWidthTruth->SetMarkerStyle(1);
  hGaussWidthTruth->SetLineColor(4);
  hGaussWidthTruth->SetLineStyle(2);
  hGaussWidthRatio = static_cast<TH1D*>(hGaussWidth->Clone("hGaussWidthRatio"));
  hGaussWidthRatio->SetTitle(";p_{T} (GeV);#sigma_{fit} / #sigma_{true}");
  hGaussWidthRatioErr = static_cast<TH1D*>(hGaussWidthRatio->Clone("hGaussWidthRatioErr"));
  hGaussWidthRatioErr->SetMarkerStyle(1);
  hGaussWidthRatioErr->SetLineWidth(0);
  hGaussWidthRatioErr->SetLineColor(43);
  hGaussWidthRatioErr->SetFillColor(43);

  hGaussWidthMC = new TH1D("hGaussWidthMC",";p_{T} (GeV);#sigma(p_{T}) / p_{T}",
			   nPtBins(),&(ptBinEdges_.front()));
  hGaussWidthMC->SetMarkerStyle(20);
  hGaussWidthMC->SetMarkerColor(1);
  hGaussWidthMC->SetLineColor(1);
  hGaussWidthMC->SetLineStyle(1);
  hGaussWidthRatioMC = static_cast<TH1D*>(hGaussWidthMC->Clone("hGaussWidthRatioMC"));
  hGaussWidthRatioMC->SetMarkerStyle(20);
  hGaussWidthRatioMC->SetMarkerColor(1);
  hGaussWidthRatioMC->SetLineColor(1);
  hGaussWidthRatioMC->SetLineStyle(1);
  hGaussWidthPseudo = static_cast<TH1D*>(hGaussWidthMC->Clone("hGaussWidthPseudo"));
  hGaussWidthPseudo->SetMarkerStyle(21);
  hGaussWidthPseudo->SetMarkerColor(2);
  hGaussWidthPseudo->SetLineColor(2);
  hGaussWidthRatioPseudo = static_cast<TH1D*>(hGaussWidthPseudo->Clone("hGaussWidthRatioPseudo"));

  hJESFrame = new TH1D("hJESFrame",";p_{T} (GeV);JES",
		       nPtBins(),&(ptBinEdges_.front()));
  for(int bin = 1; bin <= hJESFrame->GetNbinsX(); ++bin) {
    hJESFrame->SetBinContent(bin,1.);
  }
  hJESFrame->GetYaxis()->SetRangeUser(0.95,1.1);
  hJESFrame->SetLineStyle(2);
  hJES = static_cast<TH1D*>(hJESFrame->Clone("hJES"));
  hJES->SetMarkerStyle(20);
  hJES->SetLineStyle(1);
  hJESJet1 = static_cast<TH1D*>(hJES->Clone("hJESJet1"));
  hJESJet1->SetMarkerStyle(21);
  hJESJet1->SetMarkerColor(2);
  hJESJet1->SetLineColor(2);
  hJESJet2 = static_cast<TH1D*>(hJES->Clone("hJESJet2"));
  hJESJet2->SetMarkerStyle(22);
  hJESJet2->SetMarkerColor(4);
  hJESJet2->SetLineColor(4);


  // --- Fill histograms of measured response --------------
  for(DataIt datait = data_->begin(); datait != data_->end(); datait++) {
    // Select DiJet events
    if( (*datait)->type() == TypeSmearDiJet )  {
      SmearDiJet * dijet = static_cast<SmearDiJet*>(*datait);  

      hPtHat->Fill( dijet->ptHat(), dijet->weight() );

      const Jet *j1 = dijet->jet1();
      const Jet *j2 = dijet->jet2();
      if( rand_->Uniform() > 0.5 ) {
	j1 = dijet->jet2();
	j2 = dijet->jet1();
      }
      double ptGenAsym = j1->genPt() - j2->genPt();
      ptGenAsym /= j1->genPt() + j2->genPt();
      double ptAsym = j1->pt() - j2->pt();
      ptAsym /= j1->pt() + j2->pt();
      double ptAsymBiased = j1->pt() - j2->pt();
      if( j2->pt() > j1->pt() ) ptAsymBiased = j2->pt() - j1->pt();
      ptAsymBiased /= j1->pt() + j2->pt();

      for(int i = 0; i < 2; i++) {        // Loop over both jets
	const Jet * jet = dijet->jet1();
	if( i == 1 ) jet = dijet->jet2();

	hPtGenAbs->Fill( jet->genPt(), dijet->weight() );
	hPtGen->Fill( jet->genPt(), dijet->weight() );
	if( i == 0 ) hPtGenJet1->Fill( jet->genPt(), dijet->weight() );
	if( i == 1 ) hPtGenJet2->Fill( jet->genPt(), dijet->weight() );

	for(int bin = 0; bin < nPtBins(); bin++) {
	  double var = 0.;
	  if( ptBinningVar_ == "ptGen" ) var = jet->genPt();
	  else if( ptBinningVar_ == "pt" ) var = jet->pt();
	  else if( ptBinningVar_ == "ptDijet" ) var = dijet->dijetPt();
	  if( ptBinEdges_[bin] <= var && var < ptBinEdges_[bin+1] ) {
	    hRespMeasAbs[bin]->Fill( jet->pt() / jet->genPt(), dijet->weight() );
	    hRespMeas[bin]->Fill( jet->pt() / jet->genPt(), dijet->weight() );
	    if( i == 0 ) {
	      hRespMeasJet1[bin]->Fill( jet->pt() / jet->genPt(), dijet->weight() );
	      hRespRatioJet1[bin]->Fill( jet->pt() / jet->genPt(), dijet->weight() );
	    } else if ( i == 1 ) {
	      hRespMeasJet2[bin]->Fill( jet->pt() / jet->genPt(), dijet->weight() );
	      hRespRatioJet2[bin]->Fill( jet->pt() / jet->genPt(), dijet->weight() );
	    }
	    hPtGenAbsBins[bin]->Fill( jet->genPt(), dijet->weight() );
	    hPtGenAsym[bin]->Fill( ptGenAsym, dijet->weight() );
	    hPtAsym[bin]->Fill( ptAsym, dijet->weight() );
	    hPtAsymBiased[bin]->Fill( ptAsymBiased, dijet->weight() );
	    hRespMCPtHat[bin]->Fill( jet->pt()/dijet->ptHat(), dijet->weight() );
	    continue;
	  }
	}
      }
      hPtDijet->Fill( dijet->dijetPt(), dijet->weight() );
    }
  } // End of loop over data
  for(int ptBin = 0; ptBin < nPtBins(); ptBin++) {
    normHist(hRespMeas[ptBin],"width");
    normHist(hRespMeasJet1[ptBin],"width");
    normHist(hRespMeasJet2[ptBin],"width");
    normHist(hRespRatioJet1[ptBin],"width");
    normHist(hRespRatioJet2[ptBin],"width");
    normHist(hRespMCPtHat[ptBin],"width");
    normHist(hPtGenAsym[ptBin],"width");
    normHist(hPtAsym[ptBin],"width");
    normHist(hPtAsymBiased[ptBin],"width");
  }
  normHist(hPtGen,tMin,tMax,"width");
  normHist(hPtGenJet1,tMin,tMax,"width");
  normHist(hPtGenJet2,tMin,tMax,"width");
  normHist(hPtHat,tMin,tMax,"width");
  normHist(hPtDijet,tMin,tMax,"width");


  // --- Fill histograms of fitted response ----------------
  // Get parameters
  std::vector<double> fittedPar(param_->GetNumberOfParameters());
  for(int i = 0; i < param_->GetNumberOfParameters(); i++) {
    fittedPar.at(i) = param_->GetPars()[i];
  }

  std::vector<double> auxPar = bag_of<double>(config_->read<string>("mean response parameters","1 0"));
  SmearData * smearData = dynamic_cast<SmearData*>(data_->front());
  if( smearData ) {
    // Loop over ptBins
    for(int ptBin = 0; ptBin < nPtBins(); ptBin++) {
      // Interpolated response function
      for(int bin = 1; bin <= hRespFit[ptBin]->GetNbinsX(); bin++) {
	double r = hRespFit[ptBin]->GetBinCenter(bin);
	double pt = hPtGenAbsBins[ptBin]->GetMean();
	double val = smearData->pdfResp(r,pt);
	hRespFit[ptBin]->SetBinContent(bin,val);
	hRespFitErrStat[ptBin]->SetBinContent(bin,val);
	hRespFitErrStat[ptBin]->SetBinError(bin,smearData->pdfRespError(r,pt));

	hRespFitBins[ptBin]->Fill(r,val);
	hRespRatio[ptBin]->Fill(r,val);
	hRespRatioJet1[ptBin]->Fill(r,val);
	hRespRatioJet2[ptBin]->Fill(r,val);
      }
      hRespFitBins[ptBin]->Scale(1.*hRespFitBins[ptBin]->GetNbinsX()/hRespFit[ptBin]->GetNbinsX());
      hRespRatio[ptBin]->Scale(1.*hRespRatio[ptBin]->GetNbinsX()/hRespFit[ptBin]->GetNbinsX());
      hRespRatioJet1[ptBin]->Scale(1.*hRespRatioJet1[ptBin]->GetNbinsX()/hRespFit[ptBin]->GetNbinsX());
      hRespRatioJet2[ptBin]->Scale(1.*hRespRatioJet2[ptBin]->GetNbinsX()/hRespFit[ptBin]->GetNbinsX());

      // Ratio plots
      hRespRatio[ptBin]->Divide(hRespMeas[ptBin]);
      hRespRatioJet1[ptBin]->Divide(hRespMeasJet1[ptBin]);
      hRespRatioJet2[ptBin]->Divide(hRespMeasJet2[ptBin]);

      // Pt asymmetry
      for(int bin = 1; bin <= hFitPtAsym[ptBin]->GetNbinsX(); bin++) {
	double a = hFitPtAsym[ptBin]->GetBinCenter(bin);
	double pt = hPtGenAbsBins[ptBin]->GetMean();
	hFitPtAsym[ptBin]->SetBinContent(bin,smearData->pdfDijetAsym(a,pt));
      }

      // Interpolated fit function with start values
      // Copy start values into parameter array
      for(int i = 0; i < param_->GetNumberOfParameters(); i++) {
	param_->GetPars()[i] = startParJet_.at(i);
      }
      // Plot response function
      for(int bin = 1; bin <= hRespFitStart[ptBin]->GetNbinsX(); bin++) {
	double pt = hPtGenAbsBins[ptBin]->GetMean();
	double r = hRespFitStart[ptBin]->GetBinCenter(bin);
	hRespFitStart[ptBin]->SetBinContent(bin,smearData->pdfResp(r,pt));
      }
      // Copy back fitted values into parameter array
      for(int i = 0; i < param_->GetNumberOfParameters(); i++) {
	param_->GetPars()[i] = fittedPar.at(i);
      }
    } // End of loop over ptBins
  } // End if( smearData )


  // --- Fill histograms of fitted Gaussian width -----------
  if( param == "SmearParametrizationGauss" 
      || param == "SmearParametrizationGaussExtrapolation"
      || param == "SmearParametrizationGaussImbalance" ) {
    std::vector<double> truthPar = bag_of<double>(config_->read<string>("plots true resolution parameters",""));
    for(int bin = 1; bin <= hGaussWidth->GetNbinsX(); ++bin) {
      double pt = hGaussWidth->GetBinCenter(bin);
      double sigma = gaussianWidth(pt) / pt;
      double err = gaussianWidthError(pt) / pt;
      hGaussWidth->SetBinContent(bin,sigma);
      hGaussWidthErr->SetBinContent(bin,sigma);
      hGaussWidthErr->SetBinError(bin,err);
      double sigmaTrue = gaussianWidthTruth(pt) / pt;
      hGaussWidthTruth->SetBinContent(bin,sigmaTrue);
      hGaussWidthRatio->SetBinContent(bin,sigma/sigmaTrue);
      hGaussWidthRatioErr->SetBinContent(bin,hGaussWidthRatio->GetBinContent(bin));
      hGaussWidthRatioErr->SetBinError(bin,err/sigmaTrue);
    }

    for(int bin = 1; bin <= hJESFrame->GetNbinsX(); ++bin) {
      for(int i = 0; i < 3; ++i) {
	TH1 *hScale = 0;
	TH1 *hReso = 0;
	TH1 *hMC = 0;
	if( i == 0 ) {
	  hScale = hJES;
	  hMC = hRespMeas[bin-1];
	  hReso = hGaussWidthMC;
	} else if( i == 1 ) {
	  hScale = hJESJet1;
	  hMC = hRespMeasJet1[bin-1];
	} else if( i == 2 ) {
	  hScale = hJESJet2;
	  hMC = hRespMeasJet2[bin-1];
	}

	// Fit central Gaussian
	double scale = hMC->GetMean();
	double scaleErr = hMC->GetMeanError();
	double sigma = hMC->GetRMS();
	double sigmaErr = hMC->GetRMSError();
	if( scale ) {
	  hMC->Fit("gaus","0QIR","",scale-1.5*sigma,scale+1.5*sigma);
	  TF1 *fit = hMC->GetFunction("gaus");
	  scale = fit->GetParameter(1);
	  scaleErr = fit->GetParError(1);
	  sigma = fit->GetParameter(2);
	  sigmaErr = fit->GetParError(2);
	}
	hScale->SetBinContent(bin,scale);
	hScale->SetBinError(bin,scaleErr);
	if( hReso ) {
	  hReso->SetBinContent(bin,sigma);
	  hReso->SetBinError(bin,sigmaErr);
	}
	// Perform pseudo fit (scale influence)
	if( i == 0 && scale ) {						    
	  TF1 *fit = new TF1("PseudoFit","gaus",0.,2.);
	  fit->SetParameter(0,1./sqrt(2.*M_PI)/sigma);
	  fit->FixParameter(1,1.);
	  fit->SetParameter(2,sigma);
	  hMC->Fit(fit,"0QIRB");

	  hGaussWidthPseudo->SetBinContent(bin,fit->GetParameter(2));
	  hGaussWidthPseudo->SetBinError(bin,fit->GetParError(2));

	  delete fit;
	}
      }

      double pt = hPtGenAbsBins[bin-1]->GetMean();
      double sigmaTrue = gaussianWidthTruth(pt) / pt;
      hGaussWidthRatioMC->SetBinContent(bin,hGaussWidthMC->GetBinContent(bin)/sigmaTrue);
      hGaussWidthRatioMC->SetBinError(bin,hGaussWidthMC->GetBinError(bin)/sigmaTrue);
      hGaussWidthRatioPseudo->SetBinContent(bin,hGaussWidthPseudo->GetBinContent(bin)/sigmaTrue);
      hGaussWidthRatioPseudo->SetBinError(bin,hGaussWidthPseudo->GetBinError(bin)/sigmaTrue);
    }
  }


  // --- Fill histograms of fitted truth spectrum -----------

  // Fill histogram of assumed dijet truth pdf
  for(DataIt datait = data_->begin(); datait != data_->end(); datait++) {
    if( (*datait)->type() == TypeSmearDiJet ) {
      SmearDiJet * dijet = static_cast<SmearDiJet*>(*datait);  
      
      for(int bin = 1; bin <= hTruthPDF->GetNbinsX(); bin++) {
	double t = hTruthPDF->GetBinCenter(bin);
	hTruthPDF->SetBinContent(bin,dijet->pdfPtTrue(t));
      }
      for(int bin = 1; bin <= hTruthPDF->GetNbinsX(); bin++) {
	hTruthPDFErrStat->SetBinContent(bin,hTruthPDF->GetBinContent(bin));
	double t = hTruthPDFErrStat->GetBinCenter(bin);
	hTruthPDFErrStat->SetBinError(bin,dijet->pdfPtTrueError(t));
      }
      break;
    }
  }


  // --- Set x-axis ranges ----------------------------------
  for(int ptBin = 0; ptBin < nPtBins(); ptBin++) {
    hRespMeasAbs[ptBin]->GetXaxis()->SetRangeUser(rMin,rMax);   
    hRespMeas[ptBin]->GetXaxis()->SetRangeUser(rMin,rMax);      
    hRespMeasJet1[ptBin]->GetXaxis()->SetRangeUser(rMin,rMax);      
    hRespMeasJet2[ptBin]->GetXaxis()->SetRangeUser(rMin,rMax);      
    hRespMCPtHat[ptBin]->GetXaxis()->SetRangeUser(rMin,rMax);  
    hRespFitStart[ptBin]->GetXaxis()->SetRangeUser(rMin,rMax);  
    hRespFit[ptBin]->GetXaxis()->SetRangeUser(rMin,rMax);       
    hRespFitErrStat[ptBin]->GetXaxis()->SetRangeUser(rMin,rMax);
    hRespFitStep[ptBin]->GetXaxis()->SetRangeUser(rMin,rMax);   
    hRespFitGaus[ptBin]->GetXaxis()->SetRangeUser(rMin,rMax);   
    hRespFitSum[ptBin]->GetXaxis()->SetRangeUser(rMin,rMax);    
  }

  // --- Set y-axis ranges ----------------------------------
  for(int ptBin = 0; ptBin < nPtBins(); ptBin++) {
    double minFac = 0.;
    double maxFac = 1.6;
    
    double min = 0.;
    double max = 0.;
    findYRange(hRespMeas[ptBin],min,max);
    min *= minFac;
    max *= maxFac;
    hRespMeas[ptBin]->GetYaxis()->SetRangeUser(min,max);
    hRespMeasJet2[ptBin]->GetYaxis()->SetRangeUser(min,max);
    hRespMeasJet1[ptBin]->GetYaxis()->SetRangeUser(min,max);

    findYRange(hRespMCPtHat[ptBin],min,max);
    min *= minFac;
    max *= maxFac;
    hRespMCPtHat[ptBin]->GetYaxis()->SetRangeUser(min,max);

    findYRange(hPtAsym[ptBin],min,max);
    min *= minFac;
    max *= maxFac;
    hPtAsym[ptBin]->GetYaxis()->SetRangeUser(min,max);

    findYRange(hPtAsymBiased[ptBin],min,max);
    min *= minFac;
    max *= maxFac;
    hPtAsymBiased[ptBin]->GetYaxis()->SetRangeUser(min,max);

    findYRange(hPtGenAsym[ptBin],min,max);
    min *= minFac;
    max *= maxFac;
    hPtGenAsym[ptBin]->GetYaxis()->SetRangeUser(min,max);

    setYRange(hRespMeasAbs[ptBin],0.5,50.);
  }
  setYRange(hPtDijet, 0.5, 100.);
  setYRange(hPtGen, 0.5, 100.);
  setYRange(hPtGenJet1, 0.5, 100.);
  setYRange(hPtGenJet2, 0.5, 100.);
  setYRange(hPtHat, 0.5, 100.);


  // --- Plot histograms -----------------------------------
  // Label bins
  std::vector<TLegend*> legPtRange(nPtBins());
  std::vector<TLegend*> legPtRangeAndCenters(nPtBins());
  for(int ptBin = 0; ptBin < nPtBins(); ptBin++) {
    legPtRange[ptBin] = new TLegend(0.23,0.72,0.78,0.8);
    legPtRange[ptBin]->SetBorderSize(0);
    legPtRange[ptBin]->SetFillColor(0);
    legPtRange[ptBin]->SetTextFont(42);

    legPtRangeAndCenters[ptBin] = new TLegend(0.23,0.65,0.8,0.8);
    legPtRangeAndCenters[ptBin]->SetBorderSize(0);
    legPtRangeAndCenters[ptBin]->SetFillColor(0);
    legPtRangeAndCenters[ptBin]->SetTextFont(42);

    std::string binVar;
    if( ptBinningVar_ == "ptDijet" ) binVar = "p^{dijet}_{T}";
    else if( ptBinningVar_ == "pt" ) binVar = "p_{T}";
    else if( ptBinningVar_ == "ptGen" ) binVar = "p^{gen}_{T}";

    std::string label = toString(ptBinEdges_[ptBin])
      + " < " + binVar + " < "
      + toString(ptBinEdges_[ptBin+1])
      + " GeV";
    legPtRange[ptBin]->AddEntry(hRespMeas[ptBin],label.c_str(),"P");
    legPtRangeAndCenters[ptBin]->AddEntry(hRespMeas[ptBin],label.c_str(),"P");
    label = "p_{T} = " + toString(hPtGenAbsBins[ptBin]->GetMean()) + " GeV";
    legPtRangeAndCenters[ptBin]->AddEntry(hRespFit[ptBin],"Fit","L");
  }

  // Write histos to ps file
  TPostScript * const ps = new TPostScript((dir_+"/jsResponse.ps").c_str(),111);
  TCanvas *c1 = new TCanvas("c1","Jet Response",0,0,600,600);

  TLegend *legFitStart = new TLegend(0.23,0.5,0.5,0.65);
  legFitStart->SetBorderSize(0);
  legFitStart->SetFillColor(0);
  legFitStart->SetTextFont(42);
  legFitStart->AddEntry(hRespFitStart[0],"At start","L");
  legFitStart->AddEntry(hRespFit[0],"After fit","L");
  
  int logy = 0;
  for(int ptBin = 0; ptBin < nPtBins(); ptBin++) {
    // MC truth and fitted response
    //     ps->NewPage();
    //     c1->cd();
    //     hRespMeasAbs[ptBin]->Draw("PE1");
    //     legPtRange[ptBin]->Draw("same");
    //     c1->SetLogy();
    //     c1->Draw();

    //     ps->NewPage();
    //     c1->cd();
    //     hRespMeas[ptBin]->Draw("PE1");
    //     legPtRange[ptBin]->Draw("same");
    //     c1->SetLogy();
    //     c1->Draw();

    //     ps->NewPage();
    //     c1->cd();
    //     hRespMeas[ptBin]->Draw("PE1");
    //     hRespFitErrStat[ptBin]->Draw("E3 same");
    //     hRespMeas[ptBin]->Draw("same");
    //     hRespFit[ptBin]->Draw("Lsame");
    //     gPad->RedrawAxis();
    //     legPtRangeAndCenters[ptBin]->Draw("same");
    //     c1->SetLogy(logy);
    //     c1->Draw();

    ps->NewPage();
    c1->cd();
    hRespMeas[ptBin]->Draw("PE1");
    hRespFit[ptBin]->Draw("Lsame");
    hRespFitStart[ptBin]->Draw("Lsame");
    gPad->RedrawAxis();
    c1->SetLogy(logy);
    legPtRangeAndCenters[ptBin]->Draw("same");
    legFitStart->Draw("same");
    c1->Draw();

    ps->NewPage();
    c1->cd();
    hRespMeas[ptBin]->Draw("PE1");
    if( param == "SmearParametrizationStepGaussInter" ) {
      hRespFitStep[ptBin]->Draw("same");
      hRespFitGaus[ptBin]->Draw("same");
      //hRespFitSum[ptBin]->Draw("same");
    }
    hRespFit[ptBin]->Draw("Lsame");
    gPad->RedrawAxis();
    legPtRangeAndCenters[ptBin]->Draw("same");
    c1->SetLogy(logy);
    c1->Draw();

    ps->NewPage();
    c1->cd();
    hRespMeas[ptBin]->Draw("PE1");
    hRespFitBins[ptBin]->Draw("Hsame");
    gPad->RedrawAxis();
    legPtRangeAndCenters[ptBin]->Draw("same");
    c1->SetLogy(0);
    c1->Draw();
    c1->SetLogy(0);

    ps->NewPage();
    c1->cd();
    hRespRatioFrame->Draw("L");
    hRespRatio[ptBin]->Draw("PE1same");
    c1->SetLogy(0);
    c1->Draw();

    ps->NewPage();
    c1->cd();
    hRespMeas[ptBin]->GetYaxis()->SetRangeUser(5E-5,5E3);
    hRespMeas[ptBin]->Draw("PE1");
    hRespFit[ptBin]->Draw("Lsame");
    gPad->RedrawAxis();
    legPtRangeAndCenters[ptBin]->Draw("same");
    c1->SetLogy(1);
    c1->Draw();
    c1->SetLogy(0);

    ps->NewPage();
    c1->cd();
    hRespMeasJet1[ptBin]->Draw("PE1");
    hRespFit[ptBin]->Draw("Lsame");
    gPad->RedrawAxis();
    legPtRangeAndCenters[ptBin]->Draw("same");
    c1->Draw();

    ps->NewPage();
    c1->cd();
    hRespMeasJet1[ptBin]->Draw("PE1");
    hRespFitBins[ptBin]->Draw("Hsame");
    gPad->RedrawAxis();
    legPtRangeAndCenters[ptBin]->Draw("same");
    c1->SetLogy(0);
    c1->Draw();

    ps->NewPage();
    c1->cd();
    hRespRatioFrameJet1->Draw("L");
    hRespRatioJet1[ptBin]->Draw("PE1same");
    c1->SetLogy(0);
    c1->Draw();

    ps->NewPage();
    c1->cd();
    hRespMeasJet2[ptBin]->Draw("PE1");
    hRespFit[ptBin]->Draw("Lsame");
    gPad->RedrawAxis();
    legPtRangeAndCenters[ptBin]->Draw("same");
    c1->Draw();

    ps->NewPage();
    c1->cd();
    hRespMeasJet2[ptBin]->Draw("PE1");
    hRespFitBins[ptBin]->Draw("Hsame");
    gPad->RedrawAxis();
    legPtRangeAndCenters[ptBin]->Draw("same");
    c1->SetLogy(0);
    c1->Draw();

    ps->NewPage();
    c1->cd();
    hRespRatioFrameJet2->Draw("L");
    hRespRatioJet2[ptBin]->Draw("PE1same");
    c1->SetLogy(0);
    c1->Draw();

    ps->NewPage();
    c1->cd();
    std::vector<TF1*> fitResp(2);
    std::vector<TLine*> fitRespLine(2);
    for(int jetIdx = 0; jetIdx < 2; ++jetIdx) {
      TH1 *hMC = 0;
      if( jetIdx == 0 ) {
	hMC = hRespMeasJet1[ptBin];
	hMC->SetMarkerColor(2*(jetIdx+1));
	hMC->SetLineColor(2*(jetIdx+1));
	hMC->Draw("PE1");
      } else if( jetIdx == 1 ) {
	hMC = hRespMeasJet2[ptBin];
	hMC->SetMarkerColor(2*(jetIdx+1));
	hMC->SetLineColor(2*(jetIdx+1));
	hMC->Draw("PE1same");
      }
      TString fitName = "fitResp";
      fitName += jetIdx;
      fitResp[jetIdx] = new TF1(fitName,"gaus",0.3,1.7);
      hMC->Fit(fitResp[jetIdx],"0QIR");
      fitResp[jetIdx]->SetLineStyle(2);
      fitResp[jetIdx]->SetLineWidth(2);
      fitResp[jetIdx]->SetLineColor(2*(jetIdx+1));
      fitResp[jetIdx]->Draw("same");

      fitRespLine[jetIdx] = new TLine(fitResp[jetIdx]->GetParameter(1),0.,
				      fitResp[jetIdx]->GetParameter(1),
				      fitResp[jetIdx]->Eval(fitResp[jetIdx]->GetParameter(1)));
      fitRespLine[jetIdx]->SetLineColor(2*(jetIdx+1));
      fitRespLine[jetIdx]->SetLineStyle(2);
      fitRespLine[jetIdx]->Draw("same");
    }
    hRespFit[ptBin]->Draw("Lsame");
    gPad->RedrawAxis();
    legPtRangeAndCenters[ptBin]->Draw("same");
    c1->Draw();
    for(int jetIdx = 0; jetIdx < 2; ++jetIdx) {
      delete fitResp[jetIdx];
      delete fitRespLine[jetIdx];
    }

    //     ps->NewPage();
    //     c1->cd();
    //     hRespMCPtHat[ptBin]->Draw();
    //     hRespFit[ptBin]->Draw("Lsame");
    //     gPad->RedrawAxis();
    //     legPtRangeAndCenters[ptBin]->Draw("same");
    //     c1->SetLogy(logy);
    //     c1->Draw();
//     ps->NewPage();
//     c1->cd();
//     hPtGenAsym[ptBin]->Draw("PE1");
//     legPtRange[ptBin]->Draw("same");
//     c1->SetLogy(logy);
//     c1->Draw();


    ps->NewPage();
    c1->cd();
    hPtAsym[ptBin]->Draw("PE1");
    TF1 *fitAsym = new TF1("fitAsym","gaus",-0.5,0.5);
    fitAsym->SetParameter(0,1./sqrt(2.*M_PI)/0.1);
    fitAsym->FixParameter(1,0.);
    fitAsym->SetParameter(2,0.1);
    hPtAsym[ptBin]->Fit(fitAsym,"Q0IRB");
    fitAsym->SetLineStyle(2);
    fitAsym->SetLineWidth(2);
    fitAsym->Draw("same");
    hFitPtAsym[ptBin]->Draw("Lsame");
    legPtRangeAndCenters[ptBin]->Draw("same");
    gPad->RedrawAxis();
    c1->SetLogy(logy);
    c1->Draw();

    ps->NewPage();
    c1->cd();
    hPtAsymBiased[ptBin]->Draw("PE1");
    legPtRange[ptBin]->Draw("same");
    gPad->RedrawAxis();
    c1->SetLogy(logy);
    c1->Draw();

    delete fitAsym;
  }

  // Truth spectrum
  delete legPtRange[0];
  legPtRange[0] = new TLegend(0.23,0.72,0.78,0.8);
  legPtRange[0]->SetBorderSize(0);
  legPtRange[0]->SetFillColor(0);
  legPtRange[0]->SetTextFont(42);
  std::string binVar = config_->read<std::string>("plots pt binning","");
  if( binVar.find("ptDijet") != std::string::npos ) binVar = "p^{dijet}_{T}";
  else if( binVar.find("ptGen") != std::string::npos ) binVar = "p^{gen}_{T}";
  else binVar = "p^{gen}_{T}";
  std::string label = toString(ptBinEdges_.front())
    + " < " + binVar + " < "
    + toString(ptBinEdges_.back())
    + " GeV";
  legPtRange[0]->AddEntry(hPtGen,label.c_str(),"L");

  ps->NewPage();
  c1->cd();
  hPtGenAbs->Draw("PE1");
  legPtRange[0]->Draw("same");
  c1->SetLogy();
  c1->Draw();

  ps->NewPage();
  c1->cd();
  hPtGen->Draw("PE1");
  //hTruthPDFErrStat->Draw("E3same");
  //hPtGen->Draw("same");
  hTruthPDF->Draw("Lsame");
  gPad->RedrawAxis();
  legPtRange[0]->Draw("same");
  c1->SetLogy();
  c1->Draw();

  ps->NewPage();
  c1->cd();
  double tmpMin = 0.;
  double tmpMax = 0.;
  findYRange(hPtGen,tmpMin,tmpMax);
  tmpMax *= 1.4;
  hPtGen->GetYaxis()->SetRangeUser(0.,tmpMax);
  hPtGen->Draw("PE1");
  hTruthPDF->Draw("Lsame");
  gPad->RedrawAxis();
  legPtRange[0]->Draw("same");
  c1->SetLogy(0);
  c1->SetLogx(0);
  c1->Draw();
  c1->SetLogx(0);

  ps->NewPage();
  c1->cd();
  hPtGenJet1->Draw("PE1");
  hTruthPDF->Draw("Lsame");
  gPad->RedrawAxis();
  legPtRange[0]->Draw("same");
  c1->SetLogy();
  c1->Draw();

  ps->NewPage();
  c1->cd();
  tmpMin = 0.;
  tmpMax = 0.;
  findYRange(hPtGenJet1,tmpMin,tmpMax);
  tmpMax *= 1.4;
  hPtGenJet1->GetYaxis()->SetRangeUser(0.,tmpMax);
  hPtGenJet1->Draw("PE1");
  hTruthPDF->Draw("Lsame");
  gPad->RedrawAxis();
  legPtRange[0]->Draw("same");
  c1->SetLogy(0);
  c1->SetLogx(0);
  c1->Draw();
  c1->SetLogx(0);

  std::vector<TObject*> objs;
  objs.clear();
  objs.push_back(hPtDijet);
  objs.push_back(hTruthPDF);
  drawPSPage(ps,c1,objs,"",true);

  // Gaussian width
  if( param == "SmearParametrizationGauss"
      || param == "SmearParametrizationGaussExtrapolation"
      || param == "SmearParametrizationGaussImbalance" ) {
    TLegend *legSig = createLegend(4,0.8);
    legSig->AddEntry(hGaussWidthTruth,"True width","L");
    legSig->AddEntry(hGaussWidthMC,"MC truth distributions","P");
    legSig->AddEntry(hGaussWidth,"Fitted width","L");
    legSig->AddEntry(hGaussWidthErr,"Statistical uncertainty","F");
    
    ps->NewPage();
    c1->cd();
    hGaussWidth->Draw("L");
    //hGaussWidthErr->Draw("LE3same");
    hGaussWidthTruth->Draw("Lsame");
    hGaussWidth->Draw("Lsame");
    hGaussWidthMC->Draw("PE1same");
    legSig->Draw("same");
    c1->SetLogx(1);
    c1->SetLogy(0);
    c1->Draw();

    ps->NewPage();
    c1->cd();
    hGaussWidthRatio->GetYaxis()->SetRangeUser(0.8,1.4);
    hGaussWidthRatio->Draw("L");
    hGaussWidthRatioErr->Draw("LE3same");
    hGaussWidthRatio->Draw("Lsame");
    hGaussWidthRatioMC->Draw("PE1same");
    TLine line(hGaussWidthRatio->GetXaxis()->GetBinLowEdge(1),1.,
	       hGaussWidthRatio->GetXaxis()->GetBinLowEdge(hGaussWidthRatio->GetNbinsX()),1.);
    line.SetLineStyle(2);
    line.Draw("same");
    legSig->Draw("same");
    c1->SetLogx(1);
    c1->SetLogy(0);
    c1->Draw();

    delete legSig;
    legSig = createLegend(3);
    legSig->AddEntry(hGaussWidthTruth,"MC truth","L");
    legSig->AddEntry(hGaussWidthMC,"MC truth distributions","P");
    legSig->AddEntry(hGaussWidth,"Pseudo fit","L");

    ps->NewPage();
    c1->cd();
    hGaussWidthPseudo->Draw("PE1");
    hGaussWidthMC->Draw("PE1same");
    hGaussWidthTruth->Draw("Lsame");
    legSig->Draw("same");
    c1->SetLogx(1);
    c1->SetLogy(0);
    c1->Draw();

    ps->NewPage();
    c1->cd();
    hGaussWidthRatioPseudo->GetYaxis()->SetRangeUser(0.8,1.4);
    hGaussWidthRatioPseudo->Draw("PE1");
    hGaussWidthRatioMC->Draw("PE1same");
    line.Draw("same");
    legSig->Draw("same");
    c1->SetLogx(1);
    c1->SetLogy(0);
    c1->Draw();
  }

  TLegend *legJES = createLegend(3,0.5);
  legJES->AddEntry(hJESJet1,"Jet 1","P");
  legJES->AddEntry(hJESJet2,"Jet 2","P");
  legJES->AddEntry(hJES,"Both","P");

  ps->NewPage();
  c1->cd();
  hJESFrame->Draw("L");
  hJES->Draw("PE1same");
  hJESJet1->Draw("PE1same");
  hJESJet2->Draw("PE1same");
  legJES->Draw("same");
  c1->SetLogx(1);
  c1->SetLogy(0);
  c1->Draw();



  // Write histos to root file
  TFile rootfile((dir_+"/jsResponse.root").c_str(),"UPDATE");
  for(int ptBin = 0; ptBin < nPtBins(); ptBin++) {
    rootfile.WriteTObject(hRespMeasAbs[ptBin]);
    rootfile.WriteTObject(hRespMeas[ptBin]);
    rootfile.WriteTObject(hRespMeasJet1[ptBin]);
    rootfile.WriteTObject(hRespMeasJet2[ptBin]);
    rootfile.WriteTObject(hRespMCPtHat[ptBin]);
    rootfile.WriteTObject(hRespFit[ptBin]);
    rootfile.WriteTObject(hRespFitErrStat[ptBin]);
    rootfile.WriteTObject(hRespFitStart[ptBin]);
    rootfile.WriteTObject(hRespFitStep[ptBin]);
    rootfile.WriteTObject(hRespFitGaus[ptBin]);
    rootfile.WriteTObject(hRespFitSum[ptBin]);
    rootfile.WriteTObject(hRespFitBins[ptBin]);
    rootfile.WriteTObject(hRespRatio[ptBin]);
    rootfile.WriteTObject(hRespRatioJet1[ptBin]);
    rootfile.WriteTObject(hRespRatioJet2[ptBin]);
    rootfile.WriteTObject(hPtGenAsym[ptBin]);
    rootfile.WriteTObject(hPtAsym[ptBin]);
    rootfile.WriteTObject(hPtAsymBiased[ptBin]);
    rootfile.WriteTObject(hFitPtAsym[ptBin]);
  }
  rootfile.WriteTObject(hPtGenAbs);
  rootfile.WriteTObject(hPtGen);
  rootfile.WriteTObject(hPtGenJet1);
  rootfile.WriteTObject(hPtGenJet2);
  rootfile.WriteTObject(hPtHat);
  rootfile.WriteTObject(hPtDijet);
  rootfile.WriteTObject(hTruthPDF);
  rootfile.WriteTObject(hTruthPDFErrStat);
  rootfile.WriteTObject(hGaussWidth);
  rootfile.WriteTObject(hGaussWidthMC);
  rootfile.WriteTObject(hGaussWidthErr);
  rootfile.WriteTObject(hGaussWidthTruth);
  rootfile.WriteTObject(hGaussWidthPseudo);
  rootfile.WriteTObject(hGaussWidthRatio);
  rootfile.WriteTObject(hGaussWidthRatioErr);
  rootfile.WriteTObject(hGaussWidthRatioMC);
  rootfile.WriteTObject(hGaussWidthRatioPseudo);

  rootfile.Close();


  // --- Clean up ------------------------------------------
  for(int ptBin = 0; ptBin < nPtBins(); ptBin++) {
    delete hRespMeasAbs[ptBin];
    delete hRespMeas[ptBin];
    delete hRespMeasJet1[ptBin];
    delete hRespMeasJet2[ptBin];
    delete hRespMCPtHat[ptBin];
    delete hRespFit[ptBin];
    delete hRespFitErrStat[ptBin];
    delete hRespFitStart[ptBin];
    delete hRespFitStep[ptBin];
    delete hRespFitGaus[ptBin];
    delete hRespFitSum[ptBin];
    delete hRespFitBins[ptBin];
    delete legPtRangeAndCenters[ptBin];
    delete hRespRatio[ptBin];
    delete hRespRatioJet1[ptBin];
    delete hRespRatioJet2[ptBin];
    delete hPtGenAbsBins[ptBin];
    delete hPtGenAsym[ptBin];
    delete hPtAsym[ptBin];
    delete hPtAsymBiased[ptBin];
    delete hFitPtAsym[ptBin];
  }
  delete legFitStart;
  delete hPtGenAbs;
  delete hPtGen;
  delete hPtGenJet1;
  delete hPtGenJet2;
  delete hPtHat;
  delete hPtDijet;
  delete hTruthPDF;
  delete hTruthPDFErrStat;
  delete hGaussWidth;
  delete hGaussWidthErr;
  delete hGaussWidthMC;
  delete hGaussWidthTruth;
  delete hGaussWidthPseudo;
  delete hGaussWidthRatio;
  delete hGaussWidthRatioErr;
  delete hGaussWidthRatioMC;
  delete hGaussWidthRatioPseudo;
  delete legJES;
  delete c1;
  delete ps;
}



// --------------------------------------------------
void ControlPlotsJetSmearing::plotAsymmetrySimulation() const {
  std::cout << "Plotting asymmetry simulations" << std::endl;


  // ----- Create histograms ----------------------------------------
  std::vector< std::vector<TH1*> > hResp(nPtBins());
  std::vector< std::vector<TH1*> > hAsymGenBins(nPtBins());  // [ptBin][response]
  std::vector< std::vector<TH1*> > hAsymPtBins(nPtBins());  // [ptBin][response]
  std::vector< std::vector<TH1*> > hAsymPtAveBins(nPtBins());  // [ptBin][response]
  std::vector< std::vector<TH1*> > hAsymBiasGenBins(nPtBins());  // [ptBin][response]
  std::vector< std::vector<TH1*> > hAsymBiasPtBins(nPtBins());  // [ptBin][response]
  std::vector< std::vector<TH1*> > hAsymBiasPtAveBins(nPtBins());  // [ptBin][response]
  for(int ptBin = 0; ptBin < nPtBins(); ++ptBin) {

    for(int r = 0; r < 5; ++r) {
      TString name = "plotAsymmetrySimulation:hResp_pt";
      name += ptBin;
      name += "_r";
      name += r;
      TH1 *h = new TH1D(name,";Response;Events / #DeltaR",100,0.,2.);
      h->SetLineColor(color(r));
      hResp[ptBin].push_back(h);


      name = "plotAsymmetrySimulation:hAsymGenBins_pt";
      name += ptBin;
      name += "_r";
      name += r;
      h = new TH1D(name,";A;Events / #DeltaA",21,-1.,1);
      h->Sumw2();
      h->SetMarkerStyle(20+r);
      h->SetMarkerColor(color(r));
      h->SetLineColor(color(r));
      hAsymGenBins[ptBin].push_back(h);

      name = "plotAsymmetrySimulation:hAsymPtBins_pt";
      name += ptBin;
      name += "_r";
      name += r;
      h = static_cast<TH1D*>(h->Clone(name));
      hAsymPtBins[ptBin].push_back(h);

      name = "plotAsymmetrySimulation:hAsymPtAveBins_pt";
      name += ptBin;
      name += "_r";
      name += r;
      h = static_cast<TH1D*>(h->Clone(name));
      hAsymPtAveBins[ptBin].push_back(h);


      name = "plotAsymmetrySimulation:hAsymBiasGenBins_pt";
      name += ptBin;
      name += "_r";
      name += r;
      h = new TH1D(name,";A_{bias};Events / #DeltaA_{bias}",21,0.,1);
      h->Sumw2();
      h->SetMarkerStyle(20+r);
      h->SetMarkerColor(color(r));
      h->SetLineColor(color(r));
      hAsymBiasGenBins[ptBin].push_back(h);

      name = "plotAsymmetrySimulation:hAsymBiasPtBins_pt";
      name += ptBin;
      name += "_r";
      name += r;
      h = static_cast<TH1D*>(h->Clone(name));
      hAsymBiasPtBins[ptBin].push_back(h);

      name = "plotAsymmetrySimulation:hAsymBiasPtAveBins_pt";
      name += ptBin;
      name += "_r";
      name += r;
      h = static_cast<TH1D*>(h->Clone(name));
      hAsymBiasPtAveBins[ptBin].push_back(h);
    }
  }
  
  for(size_t r = 0; r < hAsymGenBins[0].size(); ++r) {
  }



  // ----- Simulating asymmetries -----------------------------------
  std::cout << "  Simulating asymmetries\n";

  std::vector<double> m1(hAsymGenBins.at(0).size());
  std::vector<double> m2(hAsymGenBins.at(0).size());
  Parametrization::CrystalBallFunction cbFunc;
  // Loop over data and fill histograms of asymmetry
  for(DataIt datait = data_->begin(); datait != data_->end(); datait++) {
    // Select DiJet events
    if( (*datait)->type() == TypeSmearDiJet )  {
      SmearDiJet * dijet = static_cast<SmearDiJet*>(*datait);  
      double t1 = dijet->jet1()->genPt();
      double t2 = dijet->jet2()->genPt();

      int nSimEvts = 10;
      for(int n = 0; n < nSimEvts; ++n) {
	// Gaussian response
	m1[0] = rand_->Gaus(t1,gaussianWidthTruth(t1));
	m2[0] = rand_->Gaus(t2,gaussianWidthTruth(t2));

	// Crystal Ball responses
	m1[1] = t1*cbFunc.random(1.,gaussianWidthTruth(t1)/t1,2.,2.);
	m2[1] = t2*cbFunc.random(1.,gaussianWidthTruth(t2)/t2,2.,2.);

	m1[2] = t1*cbFunc.random(1.,gaussianWidthTruth(t1)/t1,2.,4.);
	m2[2] = t2*cbFunc.random(1.,gaussianWidthTruth(t2)/t2,2.,4.);

	m1[3] = t1*cbFunc.random(1.,gaussianWidthTruth(t1)/t1,3.,2.);
	m2[3] = t2*cbFunc.random(1.,gaussianWidthTruth(t2)/t2,3.,2.);

	m1[4] = t1*cbFunc.truncRandom(1.,gaussianWidthTruth(t1)/t1,2.,2.,0.4);
	m2[4] = t2*cbFunc.truncRandom(1.,gaussianWidthTruth(t2)/t2,2.,2.,0.4);

	for(size_t r = 0; r < m1.size(); ++r) {
	  // Binning by mean ptGen
	  int bin = findPtBin(0.5*(t1+t2));
	  if( bin > -1 ) {
	    if( m1[r] > 0 && m2[r] > 0 ) {
	      // Fill asymmetry
	      if( rand_->Uniform() > 0.5 ) hAsymGenBins[bin][r]->Fill( (m1[r]-m2[r])/(m1[r]+m2[r]), dijet->weight() );
	      else                         hAsymGenBins[bin][r]->Fill( (m2[r]-m1[r])/(m1[r]+m2[r]), dijet->weight() );
	      // Fill pt ordered asymmetry
	      if( m1[r] > m2[r] ) {
		hAsymBiasGenBins[bin][r]->Fill( (m1[r]-m2[r])/(m1[r]+m2[r]), dijet->weight() );
	      } else {
		hAsymBiasGenBins[bin][r]->Fill( (m2[r]-m1[r])/(m1[r]+m2[r]), dijet->weight() );
	      }
	    }
	  }

	  // Binning by pt
	  if( rand_->Uniform() > 0.5 ) {
	    bin = findPtBin(m1[r]);
	  } else {
	    bin = findPtBin(m1[r]);
	  }
	  if( bin > -1 ) {
	    if( m1[r] > 0 && m2[r] > 0 ) {
	      // Fill asymmetry
	      if( rand_->Uniform() > 0.5 ) hAsymPtBins[bin][r]->Fill( (m1[r]-m2[r])/(m1[r]+m2[r]), dijet->weight() );
	      else                         hAsymPtBins[bin][r]->Fill( (m2[r]-m1[r])/(m1[r]+m2[r]), dijet->weight() );
	      // Fill pt ordered asymmetry
	      if( m1[r] > m2[r] ) {
		hAsymBiasPtBins[bin][r]->Fill( (m1[r]-m2[r])/(m1[r]+m2[r]), dijet->weight() );
	      } else {
		hAsymBiasPtBins[bin][r]->Fill( (m2[r]-m1[r])/(m1[r]+m2[r]), dijet->weight() );
	      }
	    }
	  }

	  // Binning by pt average
	  bin = findPtBin(0.5*(m1[r]+m2[r]));
	  if( bin > -1 ) {
	    if( m1[r] > 0 && m2[r] > 0 ) {
	      // Fill asymmetry
	      if( rand_->Uniform() > 0.5 ) hAsymPtAveBins[bin][r]->Fill( (m1[r]-m2[r])/(m1[r]+m2[r]), dijet->weight() );
	      else                         hAsymPtAveBins[bin][r]->Fill( (m2[r]-m1[r])/(m1[r]+m2[r]), dijet->weight() );
	      // Fill pt ordered asymmetry
	      if( m1[r] > m2[r] ) {
		hAsymBiasPtAveBins[bin][r]->Fill( (m1[r]-m2[r])/(m1[r]+m2[r]), dijet->weight() );
	      } else {
		hAsymBiasPtAveBins[bin][r]->Fill( (m2[r]-m1[r])/(m1[r]+m2[r]), dijet->weight() );
	      }
	    }
	  }
	}
      }
    }
  } // End of loop over data

  // Fill response
  for(int ptBin = 0; ptBin < nPtBins(); ptBin++) {
    double t = 0.5*(ptBinEdges_[ptBin]+ptBinEdges_[ptBin+1]);
    for(int rBin = 1; rBin <= hResp[ptBin][0]->GetNbinsX(); ++rBin) {
      double r = hResp[ptBin][0]->GetBinCenter(rBin);
      // Gaussian response
      hResp[ptBin][0]->SetBinContent(rBin,TMath::Gaus(r,1.,gaussianWidthTruth(t)/t,kTRUE));
      // Crystal Ball responses
      hResp[ptBin][1]->SetBinContent(rBin,cbFunc.pdf(r,1.,gaussianWidthTruth(t)/t,2.,2.));
      hResp[ptBin][2]->SetBinContent(rBin,cbFunc.pdf(r,1.,gaussianWidthTruth(t)/t,2.,4.));
      hResp[ptBin][3]->SetBinContent(rBin,cbFunc.pdf(r,1.,gaussianWidthTruth(t)/t,3.,2.));
      hResp[ptBin][4]->SetBinContent(rBin,cbFunc.truncPdf(r,1.,gaussianWidthTruth(t)/t,2.,2.,0.4));
    }
  }


  // ----- Write histos to file -------------------------------------
  std::cout << "  Writing histograms to file" << std::endl;
  std::vector<TPaveText*> labelGenBins(nPtBins());
  std::vector<TPaveText*> labelPtAveBins(nPtBins());
  std::vector<TPaveText*> labelPtBins(nPtBins());
  for(int ptBin = 0; ptBin < nPtBins(); ptBin++) {
    char label[50];

    labelGenBins[ptBin] = new TPaveText(0.2,0.7,0.5,0.8,"NDC");
    labelGenBins[ptBin]->SetBorderSize(0);
    labelGenBins[ptBin]->SetFillColor(0);
    labelGenBins[ptBin]->SetTextFont(42);
    labelGenBins[ptBin]->SetTextAlign(12);
    sprintf(label,"%.0f < p^{gen}_{T} < %.0f GeV",ptBinEdges_[ptBin],ptBinEdges_[ptBin+1]);
    labelGenBins[ptBin]->AddText(label);
    labelGenBins[ptBin]->AddText("L = 500 pb^{-1}");

    labelPtAveBins[ptBin] = new TPaveText(0.2,0.7,0.5,0.8,"NDC");
    labelPtAveBins[ptBin]->SetBorderSize(0);
    labelPtAveBins[ptBin]->SetFillColor(0);
    labelPtAveBins[ptBin]->SetTextFont(42);
    labelPtAveBins[ptBin]->SetTextAlign(12);
    sprintf(label,"%.0f < p^{ave}_{T} < %.0f GeV",ptBinEdges_[ptBin],ptBinEdges_[ptBin+1]);
    labelPtAveBins[ptBin]->AddText(label);
    labelPtAveBins[ptBin]->AddText("L = 500 pb^{-1}");

    labelPtBins[ptBin] = new TPaveText(0.2,0.7,0.5,0.8,"NDC");
    labelPtBins[ptBin]->SetBorderSize(0);
    labelPtBins[ptBin]->SetFillColor(0);
    labelPtBins[ptBin]->SetTextFont(42);
    labelPtBins[ptBin]->SetTextAlign(12);
    sprintf(label,"%.0f < p_{T} < %.0f GeV",ptBinEdges_[ptBin],ptBinEdges_[ptBin+1]);
    labelPtBins[ptBin]->AddText(label);
    labelPtBins[ptBin]->AddText("L = 500 pb^{-1}");
  }
  TLegend *leg = new TLegend(0.52,0.8-0.04*hAsymGenBins[0].size(),0.82,0.8);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextFont(42);
  leg->AddEntry(hAsymGenBins[0][0],"Gaussian","P");
  leg->AddEntry(hAsymGenBins[0][1],"CB #alpha = 2, n = 2","P");
  leg->AddEntry(hAsymGenBins[0][2],"CB #alpha = 2, n = 4","P");
  leg->AddEntry(hAsymGenBins[0][3],"CB #alpha = 3, n = 2","P");
  leg->AddEntry(hAsymGenBins[0][4],"CB #alpha = 3, n = 2, min = 0.4","P");

  // Write histos to ps file
  TString outName = outNamePrefix_+"SimAsym_";
  TPostScript * const ps = 0;
  if( !saveAsEps_ ) new TPostScript((dir_+"/jsSimulatedAsymmetry.ps").c_str(),111);
  TCanvas *c1 = new TCanvas("c1","Simulated Asymmetry",0,0,600,600);
  c1->SetLogy();

  // Responses
  for(int ptBin = 0; ptBin < nPtBins(); ptBin++) {
    if( !saveAsEps_ ) ps->NewPage();
    c1->cd();
    hResp[ptBin][0]->GetYaxis()->SetRangeUser(3E-5,4200.);
    hResp[ptBin][0]->Draw("L");
    for(size_t r = 1; r < hResp[ptBin].size(); ++r) {
      hResp[ptBin][r]->Draw("Lsame");
    }
    labelGenBins[ptBin]->Draw("same");
    leg->Draw("same");
    if( !saveAsEps_ ) {
      c1->Draw();
    } else {
      TString name = outName;
      name += "hResp_PtBin";
      name += ptBin;
      name += ".eps";
      c1->SaveAs(name,"eps");
    }
  }

  // PtGen binning
  for(int ptBin = 0; ptBin < nPtBins(); ptBin++) {
    if( !saveAsEps_ ) ps->NewPage();
    c1->cd();
    hAsymGenBins[ptBin][0]->SetMaximum(250.*hAsymGenBins[ptBin][0]->GetMaximum());
    hAsymGenBins[ptBin][0]->Draw("PE1");
    for(size_t r = 1; r < hAsymGenBins[ptBin].size(); ++r) {
      hAsymGenBins[ptBin][r]->Draw("PE1same");
    }
    labelGenBins[ptBin]->Draw("same");
    leg->Draw("same");
    if( !saveAsEps_ ) {
      c1->Draw();
    } else {
      TString name = outName;
      name += "hAsymGenBins_PtBin";
      name += ptBin;
      name += ".eps";
      c1->SaveAs(name,"eps");
    }
  }
  for(int ptBin = 0; ptBin < nPtBins(); ptBin++) {
    if( !saveAsEps_ ) ps->NewPage();
    c1->cd();
    hAsymBiasGenBins[ptBin][0]->SetMaximum(40.*hAsymBiasGenBins[ptBin][0]->GetMaximum());
    hAsymBiasGenBins[ptBin][0]->Draw("PE1");
    for(size_t r = 1; r < hAsymBiasGenBins[ptBin].size(); ++r) {
      hAsymBiasGenBins[ptBin][r]->Draw("PE1same");
    }
    labelGenBins[ptBin]->Draw("same");
    leg->Draw("same");
    if( !saveAsEps_ ) {
      c1->Draw();
    } else {
      TString name = outName;
      name += "hAsymBiasedGenBins_PtBin";
      name += ptBin;
      name += ".eps";
      c1->SaveAs(name,"eps");
    }
  }

  // Pt binning
  for(int ptBin = 0; ptBin < nPtBins(); ptBin++) {
    if( !saveAsEps_ ) ps->NewPage();
    c1->cd();
    hAsymPtBins[ptBin][0]->SetMaximum(250.*hAsymPtBins[ptBin][0]->GetMaximum());
    hAsymPtBins[ptBin][0]->Draw("PE1");
    for(size_t r = 1; r < hAsymPtBins[ptBin].size(); ++r) {
      hAsymPtBins[ptBin][r]->Draw("PE1same");
    }
    labelPtBins[ptBin]->Draw("same");
    leg->Draw("same");
    if( !saveAsEps_ ) {
      c1->Draw();
    } else {
      TString name = outName;
      name += "hAsymPtBins_PtBin";
      name += ptBin;
      name += ".eps";
      c1->SaveAs(name,"eps");
    }
  }
  for(int ptBin = 0; ptBin < nPtBins(); ptBin++) {
    if( !saveAsEps_ ) ps->NewPage();
    c1->cd();
    hAsymBiasPtBins[ptBin][0]->SetMaximum(40.*hAsymBiasPtBins[ptBin][0]->GetMaximum());
    hAsymBiasPtBins[ptBin][0]->Draw("PE1");
    for(size_t r = 1; r < hAsymBiasPtBins[ptBin].size(); ++r) {
      hAsymBiasPtBins[ptBin][r]->Draw("PE1same");
    }
    labelPtBins[ptBin]->Draw("same");
    leg->Draw("same");
    if( !saveAsEps_ ) {
      c1->Draw();
    } else {
      TString name = outName;
      name += "hAsymBiasedPtBins_PtBin";
      name += ptBin;
      name += ".eps";
      c1->SaveAs(name,"eps");
    }
  }

  // PtAve binning
  for(int ptBin = 0; ptBin < nPtBins(); ptBin++) {
    if( !saveAsEps_ ) ps->NewPage();
    c1->cd();
    hAsymPtAveBins[ptBin][0]->SetMaximum(250.*hAsymPtAveBins[ptBin][0]->GetMaximum());
    hAsymPtAveBins[ptBin][0]->Draw("PE1");
    for(size_t r = 1; r < hAsymPtAveBins[ptBin].size(); ++r) {
      hAsymPtAveBins[ptBin][r]->Draw("PE1same");
    }
    labelPtAveBins[ptBin]->Draw("same");
    leg->Draw("same");
    if( !saveAsEps_ ) {
      c1->Draw();
    } else {
      TString name = outName;
      name += "hAsymPtAveBins_PtBin";
      name += ptBin;
      name += ".eps";
      c1->SaveAs(name,"eps");
    }
  }
  for(int ptBin = 0; ptBin < nPtBins(); ptBin++) {
    if( !saveAsEps_ ) ps->NewPage();
    c1->cd();
    hAsymBiasPtAveBins[ptBin][0]->SetMaximum(40.*hAsymBiasPtAveBins[ptBin][0]->GetMaximum());
    hAsymBiasPtAveBins[ptBin][0]->Draw("PE1");
    for(size_t r = 1; r < hAsymBiasPtAveBins[ptBin].size(); ++r) {
      hAsymBiasPtAveBins[ptBin][r]->Draw("PE1same");
    }
    labelPtAveBins[ptBin]->Draw("same");
    leg->Draw("same");
    if( !saveAsEps_ ) {
      c1->Draw();
    } else {
      TString name = outName;
      name += "hAsymBiasedPtAveBins_PtBin";
      name += ptBin;
      name += ".eps";
      c1->SaveAs(name,"eps");
    }
  }
  if( !saveAsEps_ ) ps->Close();
    
  // ----- Clean up -------------------------------------------------
  for(size_t pt = 0; pt < hAsymGenBins.size(); ++pt) {
    delete labelGenBins[pt];
    delete labelPtBins[pt];
    delete labelPtAveBins[pt];
    for(size_t r = 0; r < hAsymGenBins[pt].size(); ++r) {
      delete hResp[pt][r];
      delete hAsymGenBins[pt][r];  
      delete hAsymPtBins[pt][r];  
      delete hAsymPtAveBins[pt][r];
      delete hAsymBiasGenBins[pt][r];
      delete hAsymBiasPtBins[pt][r]; 
      delete hAsymBiasPtAveBins[pt][r];
    }
  }
  delete leg;
  delete c1;
  if( !saveAsEps_ ) delete ps;
}



//!  \brief Draw control plots for events
//!         of type \p SmearDiJet
// --------------------------------------------------
void ControlPlotsJetSmearing::plotDijets() const
{
  std::cout << "Creating dijet control plots\n";

  // --- Create histograms --------------------------

  // Find pt ranges
  double minPtHat    = 10000.;
  double maxPtHat    = 0.;
  double minGenJetPt = 10000.;
  double maxGenJetPt = 0.;
  double minCalJetPt = 10000.;
  double maxCalJetPt = 0.;
  double minDijetPt  = 10000.;
  double maxDijetPt  = 0.;
  double min3rdJetPt = 10000.;
  double max3rdJetPt = 0.;
  for(DataIt datait = data_->begin(); datait != data_->end(); datait++) {
    // Select DiJet events
    if( (*datait)->type() == TypeSmearDiJet )  {
      SmearDiJet * dijet = static_cast<SmearDiJet*>(*datait);  

      // Loop over both jets
      for(int i = 0; i < 2; i++) {        
	const Jet * jet = dijet->jet1();
	if( i == 1 ) jet = dijet->jet2();

	if( jet->genPt() < minGenJetPt ) minGenJetPt = jet->genPt();
	if( jet->genPt() > maxGenJetPt ) maxGenJetPt = jet->genPt();

	if( jet->pt() < minCalJetPt ) minCalJetPt = jet->pt();
	if( jet->pt() > maxCalJetPt ) maxCalJetPt = jet->pt();
      }

      if( dijet->ptHat() < minPtHat ) minPtHat = dijet->ptHat();
      if( dijet->ptHat() > maxPtHat ) maxPtHat = dijet->ptHat();

      if( dijet->dijetPt() < minDijetPt ) minDijetPt = dijet->dijetPt();
      if( dijet->dijetPt() > maxDijetPt ) maxDijetPt = dijet->dijetPt();

      const Jet * jet3 = dijet->jet3();
      if( jet3->pt() < min3rdJetPt ) min3rdJetPt = jet3->pt();
      if( jet3->pt() > max3rdJetPt ) max3rdJetPt = jet3->pt();
    }
  }

  // Pt distributions
  TH1D * hPtHat = new TH1D("hPtHat",";#hat{p}_{T} (GeV);dN / d#hat{p}_{T}  1 / (GeV)",
			   50,0.9*minPtHat, 1.1*maxPtHat);
  hPtHat->SetLineWidth(2);
  hPtHat->SetNdivisions(505);
  hPtHat->Sumw2();

  std::vector<TH1*> hGenJetPt;
  std::vector<TH1*> hCalJetPt;
  std::string genJetNames[3] = { "hGenJetPtBothJets",
				 "hGenJetPtJet1", 
				 "hGenJetPtJet2" };
  std::string calJetNames[3] = { "hCalJetPtBothJets",
				 "hCalJetPtJet1", 
				 "hCalJetPtJet2" };
  int color[3] = { 1, 2, 4 };
  for(int i = 0; i < 3; i++) {
    TH1 * h = 0;

    h = new TH1D(genJetNames[i].c_str(),";p^{gen}_{T} (GeV);dN / dp^{gen}_{T}  1 / (GeV)",
		 50,0.9*minGenJetPt,1.1*maxGenJetPt);
    h->Sumw2();
    h->SetLineWidth(2);
    h->SetLineColor(color[i]);
    h->SetNdivisions(505);
    hGenJetPt.push_back(h);

    h = new TH1D(calJetNames[i].c_str(),";p^{jet}_{T} (GeV);dN / dp^{jet}_{T}  1 / (GeV)",
		 50,0.9*minCalJetPt,1.1*maxCalJetPt);
    h->Sumw2();
    h->SetLineWidth(2);
    h->SetLineColor(color[i]);
    h->SetNdivisions(505);
    hCalJetPt.push_back(h);
  }

  TH2 * hCalJet2vsCalJet1Pt = new TH2D("hCalJet2vsCalJet1Pt",
					";p^{jet1}_{T} (GeV);p^{jet2}_{T} (GeV)",
					50,0.9*minCalJetPt,1.1*maxCalJetPt,
					50,0.9*minCalJetPt,1.1*maxCalJetPt);

  TH1 * hDijetPt = new TH1D("hDijetPt",";p^{dijet}_{T} (GeV);1 / N  dN / dp^{dijet}_{T}  1 / (GeV)",
			     50,0.9*minDijetPt, 1.1*maxDijetPt);
  hDijetPt->SetLineWidth(2);
  hDijetPt->Sumw2();

  TH1 * h3rdJetPt = new TH1D("h3rdJetPt",";p^{jet3}_{T} (GeV);1 / N  dN / dp^{jet3}_{T}  1 / (GeV)",
			      50,0.9*min3rdJetPt, 1.1*max3rdJetPt);
  h3rdJetPt->SetLineWidth(2);
  h3rdJetPt->Sumw2();

  TH2 * h3rdJetvsDijetPt = new TH2D("h3rdJetvsDijetPt",
				     ";p^{dijet}_{T} (GeV);p^{jet3}_{T} (GeV)",
				     50,0.9*minDijetPt, 1.1*maxDijetPt,
				     50,0.9*min3rdJetPt, 1.1*max3rdJetPt);

  TH1 * hRel3rdJetPt = new TH1D("hRel3rdJetPt",
				 ";p^{jet3}_{T,rel} = p^{jet3}_{T} / p^{dijet}_{T};1 / N  dN / dp^{jet3}_{T,rel}",
				 50,0,1.4);
  hRel3rdJetPt->SetLineWidth(2);
  hRel3rdJetPt->Sumw2();

  TH1 * hDeltaPhi = new TH1D("hDeltaPhi",";#Delta#phi;1 / N  dN / d#Delta#phi",
			      25,1.5,M_PI);
  hDeltaPhi->SetLineWidth(2);
  hDeltaPhi->Sumw2();

  TH2 * hDeltaPhivsRel3rdJetPt = new TH2D("hDeltaPhivsRel3rdJetPt",
					   ";p^{jet3}_{T,rel} = p^{jet3}_{T} / p^{dijet}_{T};#Delta#phi",
					   25,0,1,25,1.5,M_PI);

  // Response correlations
  TH2 * hRvsDeltaPhi = new TH2D("hRvsDeltaPhi",";#Delta#phi;p^{jet}_{T} / p^{gen}_{T}",
				 50,1.5,M_PI,50,0,2);
  hRvsDeltaPhi->SetMarkerStyle(7);
  TH2 * hRvsRel3rdJetPt = new TH2D("hRvsRel3rdJetPt",
				    ";p^{jet3}_{T,rel} = p^{jet3}_{T} / p^{dijet}_{T};p^{jet}_{T} / p^{gen}_{T}",
				    50,0,1.4,50,0,2);
  hRvsRel3rdJetPt->SetMarkerStyle(7);
  TH2 * hRvs3rdJetPt = new TH2D("hRvs3rdJetPt",
				 ";p^{jet3} (GeV);p^{jet}_{T} / p^{gen}_{T}",
				 50,0.9*min3rdJetPt, 1.1*max3rdJetPt,50,0,2);
  hRvs3rdJetPt->SetMarkerStyle(7);
  TH2 * hRvsEMF = new TH2D("hRvsEMF",";EMF;p^{jet}_{T} / p^{gen}_{T}",
			    50,0,1,50,0,2);
  hRvsEMF->SetMarkerStyle(7);
  TH2 * hRvsDeltaR = new TH2D("hRvsDeltaR",";#Delta R(jet,genJet);p^{jet}_{T} / p^{gen}_{T}",
			       25,0,0.4,50,0,2);
  hRvsDeltaR->SetMarkerStyle(7);
  hRvsDeltaR->GetXaxis()->SetNdivisions(505);
    



  // --- Fill histograms ----------------------------
  for(DataIt datait = data_->begin(); datait != data_->end(); datait++) {
    if( (*datait)->type() == TypeSmearDiJet )  { // Select DiJet events
      SmearDiJet * dijet = static_cast<SmearDiJet*>(*datait);  

      const Jet * jet1 = dijet->jet1();
      const Jet * jet2 = dijet->jet2();
      const Jet * jet3 = dijet->jet3();

      double weight = dijet->weight();

      double dPhi = std::abs(TVector2::Phi_mpi_pi( jet1->phi() - jet2->phi() ));
      hDeltaPhi->Fill( dPhi, weight );

      hRvsDeltaPhi->Fill( dPhi, jet1->pt() / jet1->genPt(), weight );
      hRvsDeltaPhi->Fill( dPhi, jet2->pt() / jet2->genPt(), weight );

      hRvsEMF->Fill( jet1->EmEt() / jet1->pt(), jet1->pt() / jet1->genPt(), weight );
      hRvsEMF->Fill( jet2->EmEt() / jet2->pt(), jet2->pt() / jet2->genPt(), weight );

      hRvsRel3rdJetPt->Fill( jet3->pt() / dijet->dijetPt(), jet1->pt() / jet1->genPt(), weight );
      hRvsRel3rdJetPt->Fill( jet3->pt() / dijet->dijetPt(), jet2->pt() / jet2->genPt(), weight );

      hRvs3rdJetPt->Fill( jet3->pt(), jet1->pt() / jet1->genPt(), weight );
      hRvs3rdJetPt->Fill( jet3->pt(), jet2->pt() / jet2->genPt(), weight );

      hRvsDeltaR->Fill( jet1->dR(), jet1->pt() / jet1->genPt(), weight );
      hRvsDeltaR->Fill( jet2->dR(), jet2->pt() / jet2->genPt(), weight );

      hPtHat->Fill( dijet->ptHat(), weight );

      hGenJetPt.at(0)->Fill( jet1->genPt(), weight );
      hGenJetPt.at(0)->Fill( jet2->genPt(), weight );
      hGenJetPt.at(1)->Fill( jet1->genPt(), weight );
      hGenJetPt.at(2)->Fill( jet2->genPt(), weight );

      hCalJetPt.at(0)->Fill( jet1->pt(), weight );
      hCalJetPt.at(0)->Fill( jet2->pt(), weight );
      hCalJetPt.at(1)->Fill( jet1->pt(), weight );
      hCalJetPt.at(2)->Fill( jet2->pt(), weight );

      hCalJet2vsCalJet1Pt->Fill( jet1->pt(), jet2->pt(), weight );
      hDijetPt->Fill( dijet->dijetPt(), weight );

      h3rdJetPt->Fill( jet3->pt(), weight );
      h3rdJetvsDijetPt->Fill( dijet->dijetPt(), jet3->pt(), weight );
      hRel3rdJetPt->Fill( jet3->pt() / dijet->dijetPt(), weight );
      hDeltaPhivsRel3rdJetPt->Fill( jet3->pt() / dijet->dijetPt(), dPhi, weight );
    }
  }


  // Normalizing histograms
  //   for(size_t i = 0; i < hGenJetPt.size(); i++) {
  //     normHist( hGenJetPt.at(i) );
  //     normHist( hCalJetPt.at(i) );
  //   }
  //  normHist( hPtHat );
  //  normHist( hDijetPt );
  //  normHist( h3rdJetPt );
  //  normHist( hRel3rdJetPt );
  //  normHist( hDeltaPhi );


  // --- Plot histograms ----------------------------
  TPostScript * const ps = new TPostScript((dir_+"/jsDijets.ps").c_str(),111);
  TCanvas           * c1 = new TCanvas("c1","Dijets",0,0,600,600);

  // Applied cuts
  TPaveText * appCuts = new TPaveText(0.1,0.1,0.9,0.9,"NDC");
  appCuts->SetFillColor(0);
  appCuts->SetTextFont(42);
  //  appCuts->SetTextSize(0.1);
  appCuts->SetTextAlign(12);
  appCuts->SetBorderSize(0);
  appCuts->AddText(("E^{jet}_{T} > "
		    + toString(config_->read<double>("Et cut on jet",0.)) ).c_str());
  appCuts->AddText(("|#eta| < "
		    + toString(config_->read<double>("Eta cut on jet",0.)) ).c_str());
  appCuts->AddText((toString(config_->read<double>("Et min cut on dijet",0.))
		    + " < E^{dijet}_{T} < "
		    + toString(config_->read<double>("Et max cut on dijet",0.))).c_str());
  appCuts->AddText(("E^{jet3}_{T} < "
		    + toString(config_->read<double>("Et cut on n+1 Jet",0.))
		    + " or E^{jet3}_{T} / E^{dijet}_{T} < "
		    + toString(config_->read<double>("Relative n+1 Jet Et Cut",0.)) ).c_str());
  appCuts->AddText(("#Delta#phi > "
		    + toString(config_->read<double>("Min Delta Phi",0.)) ).c_str());
  appCuts->AddText((toString(config_->read<double>("Min had fraction",0.))
		    + " < E^{had}_{T} / E^{jet}_{T} < "
		    + toString(config_->read<double>("Max had fraction")) ).c_str());
  appCuts->AddText((toString(config_->read<double>("Et genJet min",0.))
		    + " < E^{gen}_{T} < "
		    + toString(config_->read<double>("Et genJet max",0.)) ).c_str());
  appCuts->AddText(("#Delta R(jet,genJet) < "
		    + toString(config_->read<double>("DeltaR cut on jet matching")) ).c_str());
  drawPSPage(ps,c1,appCuts,"",false);

  drawPSPage(ps,c1,hPtHat,"",true);
  
  TLegend * leg = new TLegend(0.7,0.75,0.9,0.9);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextFont(42);
  leg->AddEntry(hGenJetPt.at(1),"Jet 1","L");
  leg->AddEntry(hGenJetPt.at(2),"Jet 2","L");

  std::vector<TObject*> objs;
  for(size_t i = 1; i < hGenJetPt.size(); i++) {
    objs.push_back(hGenJetPt.at(i));
  }
  objs.push_back(leg);
  drawPSPage(ps,c1,hGenJetPt.at(0),"",true);
  drawPSPage(ps,c1,objs,"",true);

  objs.clear();
  for(size_t i = 1; i < hCalJetPt.size(); i++) {
    objs.push_back(hCalJetPt.at(i));
  }
  objs.push_back(leg);
  drawPSPage(ps,c1,hCalJetPt.at(0),"",true);
  drawPSPage(ps,c1,objs,"",true);

  drawPSPage(ps,c1,hCalJet2vsCalJet1Pt);
  drawPSPage(ps,c1,hDijetPt,"",true);
  drawPSPage(ps,c1,h3rdJetPt,"",true);
  drawPSPage(ps,c1,h3rdJetvsDijetPt);
  drawPSPage(ps,c1,hRel3rdJetPt,"",true);
  drawPSPage(ps,c1,hDeltaPhi,"",true);
  drawPSPage(ps,c1,hDeltaPhivsRel3rdJetPt);
  drawPSPage(ps,c1,hRvsRel3rdJetPt,"COLZ",true);
  drawPSPage(ps,c1,hRvs3rdJetPt,"COLZ",true);
  drawPSPage(ps,c1,hRvsDeltaPhi,"COLZ",true);
  drawPSPage(ps,c1,hRvsEMF,"COLZ",true);
  drawPSPage(ps,c1,hRvsDeltaR,"COLZ",true);

  ps->Close();

  // Clean up
  for(size_t i = 0; i < hGenJetPt.size(); i++) {
    delete hGenJetPt.at(i);
    delete hCalJetPt.at(i);
  }
  delete appCuts;
  delete leg;
  delete hCalJet2vsCalJet1Pt;
  delete hPtHat;
  delete hDijetPt;
  delete h3rdJetPt;
  delete h3rdJetvsDijetPt;
  delete hRel3rdJetPt;
  delete hDeltaPhi;
  delete hDeltaPhivsRel3rdJetPt;
  delete hRvsRel3rdJetPt;
  delete hRvs3rdJetPt;
  delete hRvsDeltaPhi;
  delete hRvsEMF;
  delete hRvsDeltaR;
  delete c1;
  delete ps;
}



// --------------------------------------------------
void ControlPlotsJetSmearing::plotParameters() const {
  std::cout << "Creating parameter control plots\n";

  // ----- Quantities defining the parametrization -----
  std::vector<double> scale = bag_of<double>(config_->read<string>("jet parameter scales",""));
  int nPar = param_->GetNumberOfParameters();


  // ----- Create histograms -----
  // Fitted parameter values with errors
  TH1D *hPars = new TH1D("hParameters",";Parameter index;Parameter value",nPar,-0.5,nPar-0.5);
  hPars->SetMarkerStyle(20);
  hPars->SetNdivisions(nPar);
  // Fitted absolute parameter values (par*scale) with errors
  TH1D *hAbsPars = static_cast<TH1D*>(hPars->Clone("hAbsoluteParameters"));
  hAbsPars->SetYTitle("Absolute parameter value");
  // Relative parameter errors
  TH1D *hRelParErrors = static_cast<TH1D*>(hPars->Clone("hRelativeParameterErrors"));
  hRelParErrors->SetYTitle("Relative parameter error");
  // Global parameter correlations
  TH1D *hGlobalCorr = static_cast<TH1D*>(hPars->Clone("hGlobalParameterCorrelations"));
  hGlobalCorr->SetYTitle("Global parameter correlation");
  // Parameter correlations
  TH2D *hParCorr = new TH2D("hParameterCorrelations",";Parameter index;Parameter index",
			    nPar,-0.5,nPar-0.5,nPar,-0.5,nPar-0.5);
  hParCorr->SetNdivisions(nPar,"XY");


  // ----- Fill histograms -----
  for(int i = 0; i < nPar; i++) {
    int bin = 1+i;

    hPars->SetBinContent(bin,param_->GetPars()[i]);
    hPars->SetBinError(bin,param_->GetErrors()[i]);

    hAbsPars->SetBinContent(bin,scale_[i]*param_->GetPars()[i]);
    hAbsPars->SetBinError(bin,scale_[i]*param_->GetErrors()[i]);
  
    hRelParErrors->SetBinContent(bin,param_->GetErrors()[i]/param_->GetPars()[i]);

    hGlobalCorr->SetBinContent(bin,param_->GetGlobalCorrCoeff()[i]);
  }
  for(int i = 0; i < nPar; i++) {
    for(int j = 0; j < i+1; j++) {
      int idx = (i*i + i)/2 + j;
      double corr = 0.;
      if( param_->GetErrors()[i] && param_->GetErrors()[j] ) {
	corr = param_->GetCovCoeff()[idx] / param_->GetErrors()[i] / param_->GetErrors()[j];
      }
      hParCorr->Fill(i,j,corr);
      if( i != j ) {
	hParCorr->Fill(j,i,corr);
      }
    }
  }


  // ----- Plot histograms -----
  TPostScript * const ps = new TPostScript((dir_+"/jsParameters.ps").c_str(),111);
  TCanvas           * c1 = new TCanvas("c1","Parameters",0,0,600,600);

  drawPSPage(ps,c1,hPars,"PE1");
  drawPSPage(ps,c1,hAbsPars,"PE1");
  drawPSPage(ps,c1,hRelParErrors,"P");
  drawPSPage(ps,c1,hGlobalCorr,"P");
  drawPSPage(ps,c1,hParCorr,"COLZ");

  ps->Close();

  TFile rootfile((dir_+"/jsResponse.root").c_str(),"UPDATE");
  rootfile.WriteTObject(hPars);
  rootfile.WriteTObject(hAbsPars);
  rootfile.WriteTObject(hRelParErrors);
  rootfile.WriteTObject(hGlobalCorr);
  rootfile.WriteTObject(hParCorr);


  // ----- Clean up -----
  delete hPars;
  delete hAbsPars;
  delete hRelParErrors;
  delete hGlobalCorr;
  delete hParCorr;
  delete c1;
  delete ps;
}



//! For each pair (i,j) of free parameters in the fit the
//! likelihood is plotted in the (i,j) plane. The parameters
//! are varied in 3 steps of 2*sigma (error from the fit)
//! below and above the fitted parameter value.
//!
//! The plots are written to the files "jsResponse.root"
//! and "jsParameterScans.ps".
// --------------------------------------------------
void ControlPlotsJetSmearing::plotParameterScan() const {
  std::cout << "Creating parameter scan control plots\n";

  // ----- Set up quantities -----
  int n = 4;
  int nSteps = 2*n+1;
  double nSigma = 2;
  int nPar = param_->GetNumberOfParameters();
  std::vector<TH2D*> hParScans2D;
  std::vector<TLine*> lines;

  // Store likelihood for original parameter values
  double offset = 0.;
  for(DataIt dataIt = data_->begin(); dataIt != data_->end(); dataIt++) {
    if( (*dataIt)->type() == TypeSmearDiJet )  { // Select DiJet events
      offset += (*dataIt)->chi2();
    }
  }     

  // ----- Vary parameters and calculate likelihood -----
  // Store original parameter values
  std::vector<double> origPars(nPar);
  for(int i = 0; i < nPar; i++) {
    origPars[i] = param_->GetPars()[i];
  }
  // Outer loop over parameters
  for(int i = 0; i < nPar; i++) {
    if( param_->isFixedPar(i) ) continue;
    double idVal = nSigma*param_->GetErrors()[i];
    if( idVal == 0 ) idVal = 0.1;
			
    // Inner loop over parameters
    for(int j = 0; j < i; j++) {
      if( param_->isFixedPar(j) ) continue;
      double jdVal = nSigma*param_->GetErrors()[j];
      if( jdVal == 0 ) jdVal = 0.046;
      // Create histogram of likelihood from i vs j
      TString name = "hParScan2D";
      name += hParScans2D.size();
      TString title = "- #Deltaln(L);Parameter ";
      title += i;
      if( param_->parName(i) != "" ) title += " (" + param_->parName(i) + ")";
      title += ";Parameter ";
      title += j;
      if( param_->parName(j) != "" ) title += " (" + param_->parName(j) + ")";
      hParScans2D.push_back(new TH2D(name,title,
				     nSteps,origPars[i]-(n+0.5)*idVal,origPars[i]+(n+0.5)*idVal,
				     nSteps,origPars[j]-(n+0.5)*jdVal,origPars[j]+(n+0.5)*jdVal));
			    
      // Vary parameters i and j
      for(int is = 0; is < nSteps; is++) {
	double iPar = origPars[i] + (is-n)*idVal; 
	param_->GetPars()[i] = iPar;
	for(int js = 0; js < nSteps; js++) {
	  double jPar = origPars[j] + (js-n)*jdVal; 
	  param_->GetPars()[j] = jPar;
	  double deltaLkh = 0.;
	  // Calculate likelihood for varied parameters
	  if( is != n || js != n ) {
	    for(DataIt dataIt = data_->begin(); dataIt != data_->end(); dataIt++) {
	      if( (*dataIt)->type() == TypeSmearDiJet )  { // Select DiJet events
		deltaLkh += (*dataIt)->chi2();
	      }
	    }    
	    deltaLkh -= offset;
	  }
	  hParScans2D.back()->Fill(iPar,jPar,deltaLkh);
	}
      }
      // Reset parameters to original values
      param_->GetPars()[i] = origPars[i];
      param_->GetPars()[j] = origPars[j];
    } // End of inner loop over parameters
  } // End of outer loop over parameters

  // Project out 1D parameter scans
  std::vector<TH1D*> hParScansTmp(nPar);
  int parScan2Didx = -1;
  int parScan1Didx = -1;

  // Create histograms of likelihood from i
  for(int i = 0; i < nPar; i++) {
    hParScansTmp[i] = 0;
    if( param_->isFixedPar(i) ) continue;

    // Inner loop over parameters
    for(int j = 0; j < i; j++) {
      if( param_->isFixedPar(j) ) continue;
      parScan2Didx++;
      const TH2D *h2 = hParScans2D[parScan2Didx];

      // Does 1D hist for parameter j exist?
      if( hParScansTmp[j] == 0 ) {
	TString name = "hParScan";
	name += ++parScan1Didx;
	TString title = ";";
	title += h2->GetYaxis()->GetTitle();
	title += ";- #Deltaln(L)";
	TH1D *h = new TH1D(name,title,
			   h2->GetNbinsY(),
			   h2->GetYaxis()->GetXmin(),
			   h2->GetYaxis()->GetXmax());
	// Copy y bin content in the central x bin
	for(int yBin = 1; yBin <= h2->GetNbinsY(); yBin++) {
	  h->SetBinContent(yBin,h2->GetBinContent(h2->GetBin(n+1,yBin)));
	}
	hParScansTmp[j] = h;
      }
      // Does 1D hist for parameter i exist?
      if( hParScansTmp[i] == 0 ) {
	TString name = "hParScan";
	name += ++parScan1Didx;
	TString title = ";";
	title += h2->GetXaxis()->GetTitle();
	title += ";- #Deltaln(L)";
	TH1D *h = new TH1D(name,title,
			   h2->GetNbinsX(),
			   h2->GetXaxis()->GetXmin(),
			   h2->GetXaxis()->GetXmax());
	// Copy x bin content in the central y bin
	for(int xBin = 1; xBin <= h2->GetNbinsX(); xBin++) {
	  h->SetBinContent(xBin,h2->GetBinContent(h2->GetBin(xBin,n+1)));
	}
	hParScansTmp[i] = h;
      }
    }
  }
  std::vector<TH1D*> hParScans;
  for(int i = 0; i < nPar; i++) {
    if( param_->isFixedPar(i) ) continue;
    hParScans.push_back(hParScansTmp[i]);
  }
  hParScansTmp.clear();


  // ----- Plot histograms -----
  TPostScript * const ps = new TPostScript((dir_+"/jsParameterScans.ps").c_str(),111);
  TCanvas           * c1 = new TCanvas("c1","ParameterScans",0,0,600,600);

  for(size_t i = 0; i < hParScans.size(); i++) {
    hParScans[i]->SetMarkerStyle(20);
    drawPSPage(ps,c1,hParScans[i],"P",true);
  }
  for(size_t i = 0; i < hParScans2D.size(); i++) {
    drawPSPage(ps,c1,hParScans2D[i],"COLZ",true);
  }

  ps->Close();

  TFile rootfile((dir_+"/jsResponse.root").c_str(),"UPDATE");
  for(size_t i = 0; i < hParScans.size(); i++) {
    rootfile.WriteTObject(hParScans[i]);
  }
  for(size_t i = 0; i < hParScans2D.size(); i++) {
    rootfile.WriteTObject(hParScans2D[i]);
  }


  // ----- Clean up -----
  for(size_t i = 0; i < hParScans2D.size(); i++) {
    delete hParScans2D[i];
  }
  for(size_t i = 0; i < hParScans.size(); i++) {
    delete hParScans[i];
  }
  for(size_t i = 0; i < lines.size(); i++) {
    delete lines[i];
  }
}  



//! These are the distributions of the negative logarithm
//! of the probability density of each event multiplied
//! by the event weight. That are the summands each
//! event adds to the negative log-likelihood of the fit.
//! (For a \f$\chi^{2}\f$-fit these would be the pull
//! distributions.)
//!
//! The plots are written to the files "jsResponse.root"
//! and "jsLogP.ps".
// --------------------------------------------------
void ControlPlotsJetSmearing::plotLogP() const {
  std::cout << "Creating -log(P) control plots\n";

  // ----- Fill vectors of -log(P) -----
  std::vector<double> logPstart;
  std::vector<double> logPend;
  std::vector<double> logPWstart;
  std::vector<double> logPWend;

  // Loop over data and fill -log(P) with fitted
  // parameter values
  for(DataIt dataIt = data_->begin(); dataIt != data_->end(); dataIt++) {
    // Select DiJet events
    if( (*dataIt)->type() == TypeSmearDiJet )  {
      logPWend.push_back((*dataIt)->chi2());
      logPend.push_back((*dataIt)->chi2()/(*dataIt)->weight());
    }
  }

  // Store fitted parameters
  std::vector<double> fittedPar(param_->GetNumberOfParameters());
  for(int i = 0; i < param_->GetNumberOfParameters(); i++) {
    fittedPar[i] = param_->GetPars()[i];
  }
  // Copy start values into parameter array
  std::vector<double> startParGlobal = bag_of<double>(config_->read<string>("global jet start values",""));
  for(int i = 0; i < param_->GetNumberOfParameters(); i++) {
    if( i < param_->GetNumberOfJetParameters() )
      param_->GetPars()[i] = startParJet_.at(i);
    else
      param_->GetPars()[i] = startParGlobal[i-param_->GetNumberOfJetParameters()];
  }
  // Loop over data and fill -log(P) with start
  // parameter values
  for(DataIt dataIt = data_->begin(); dataIt != data_->end(); dataIt++) {
    // Select DiJet events
    if( (*dataIt)->type() == TypeSmearDiJet )  {
      logPWstart.push_back((*dataIt)->chi2());
      logPstart.push_back((*dataIt)->chi2()/(*dataIt)->weight());
    }
  }
  // Copy back fitted values into parameter array
  for(int i = 0; i < param_->GetNumberOfParameters(); i++) {
    param_->GetPars()[i] = fittedPar.at(i);
  }


  // ----- Fill histograms of -log(P) -----
  std::sort(logPstart.begin(),logPstart.end());
  std::sort(logPend.begin(),logPend.end());
  double max = logPstart.back() > logPend.back() ? logPstart.back() : logPend.back();
  TH1F *hLogPstart = new TH1F("hLogPstart",";- ln(P)",100,0.,1.1*max);
  hLogPstart->SetLineWidth(2);
  hLogPstart->SetLineColor(4);
  TH1F *hLogPend = static_cast<TH1F*>(hLogPstart->Clone("hLogPend"));
  hLogPend->SetLineColor(2);

  std::sort(logPWstart.begin(),logPWstart.end());
  std::sort(logPWend.begin(),logPWend.end());
  max = logPWstart.back() > logPWend.back() ? logPWstart.back() : logPWend.back();
  TH1F *hLogPWstart = new TH1F("hLogPWstart",";- w #upoint ln(P)",100,0.,1.1*max);
  hLogPWstart->SetLineWidth(2);
  hLogPWstart->SetLineColor(4);
  TH1F *hLogPWend = static_cast<TH1F*>(hLogPWstart->Clone("hLogPWend"));
  hLogPWend->SetLineColor(2);

  for(size_t i = 0; i < logPend.size(); i++) {
    hLogPstart->Fill(logPstart[i]);
    hLogPend->Fill(logPend[i]);
    hLogPWstart->Fill(logPWstart[i]);
    hLogPWend->Fill(logPWend[i]);
  }
  logPstart.clear();
  logPend.clear();
  logPWstart.clear();
  logPWend.clear();

  max = hLogPend->GetMaximum() > hLogPstart->GetMaximum() ?
    hLogPend->GetMaximum() : hLogPstart->GetMaximum();
  hLogPstart->GetYaxis()->SetRangeUser(5E-2,5.*max);
  hLogPend->GetYaxis()->SetRangeUser(5E-2,5.*max);

  max = hLogPWend->GetMaximum() > hLogPWstart->GetMaximum() ?
    hLogPWend->GetMaximum() : hLogPWstart->GetMaximum();
  hLogPWstart->SetMaximum(5.*max);
  hLogPWend->SetMaximum(5.*max);


  // ----- Plot histograms -----
  TPostScript * const ps = new TPostScript((dir_+"/jsLogP.ps").c_str(),111);
  TCanvas * c1 = new TCanvas("c1","LogP",0,0,600,600);

  TLegend * leg = new TLegend(0.23,0.64,0.78,0.8);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextFont(42);
  leg->AddEntry(hLogPstart,"Before fit","L");
  leg->AddEntry(hLogPend,"After fit","L");
  
  ps->NewPage();
  c1->cd()->SetLogy();
  hLogPstart->Draw();
  hLogPend->Draw("same");
  gPad->RedrawAxis();
  leg->Draw("same");
  c1->Draw();

  ps->NewPage();
  c1->cd()->SetLogy();
  hLogPWstart->Draw();
  hLogPWend->Draw("same");
  gPad->RedrawAxis();
  leg->Draw("same");
  c1->Draw();
  ps->Close();

  TFile rootfile((dir_+"/jsResponse.root").c_str(),"UPDATE");
  rootfile.WriteTObject(hLogPstart);
  rootfile.WriteTObject(hLogPend);
  rootfile.WriteTObject(hLogPWstart);
  rootfile.WriteTObject(hLogPWend);

  // ----- Clean up -----
  delete hLogPstart;
  delete hLogPend;
  delete hLogPWstart;
  delete hLogPWend;
  delete leg;
  delete c1;
  delete ps;
}



// --------------------------------------------------
void ControlPlotsJetSmearing::plotMeanResponseAndResolution() const {
  std::cout << "Creating response and resolution fits\n";

  // ----- Create histograms -----

  // Response vs ptgen
  TH2* hRespVsPtGen = new TH2D("MeanResp_hRespVsPtGen",
			       ";p^{gen}_{T} (GeV);p^{jet}_{T} / p^{gen}_{T}",
			       nPtBins(),&(ptBinEdges_.front()),51,0.,2.);
  hRespVsPtGen->SetNdivisions(505);
  hRespVsPtGen->Sumw2();

  std::vector<double> bins(nPtBins()+1);
  equidistLogBins(bins,bins.size()-1,ptBinsMin(),ptBinsMax());
  TH2* hRespVsPtGenLog = new TH2D("MeanResp_hRespVsPtGenLog",
				  ";p^{gen}_{T} (GeV);p^{jet}_{T} / p^{gen}_{T}",
				  bins.size()-1,&(bins.front()),51,0.,2.);
  hRespVsPtGenLog->SetNdivisions(505);
  hRespVsPtGenLog->Sumw2();

  TH1* hPtGen = new TH1D("MeanResp_hPtGen",";p^{gen}_{T} (GeV);dN / dp^{gen}_{T}  1 / (GeV)",
			 60,ptBinsMin(),ptBinsMax());
  hPtGen->SetNdivisions(505);
  hPtGen->Sumw2();

  TH1* hPtGenNorm = 0;
  

  // ----- Fill histograms -----
  for(DataIt datait = data_->begin(); datait != data_->end(); datait++) {
    // Select DiJet events
    if( (*datait)->type() == TypeSmearDiJet )  {
      SmearDiJet * dijet = static_cast<SmearDiJet*>(*datait); 

      hPtGen->Fill(dijet->jet1()->genPt(),dijet->weight());
      hPtGen->Fill(dijet->jet2()->genPt(),dijet->weight());

      hRespVsPtGen->Fill(dijet->jet1()->genPt(),dijet->jet1()->pt()/dijet->jet1()->genPt(),dijet->weight());
      hRespVsPtGen->Fill(dijet->jet2()->genPt(),dijet->jet2()->pt()/dijet->jet2()->genPt(),dijet->weight());

      hRespVsPtGenLog->Fill(dijet->jet1()->genPt(),dijet->jet1()->pt()/dijet->jet1()->genPt(),dijet->weight());
      hRespVsPtGenLog->Fill(dijet->jet2()->genPt(),dijet->jet2()->pt()/dijet->jet2()->genPt(),dijet->weight());
    }
  }
  hPtGenNorm = static_cast<TH1D*>(hPtGen->Clone("MeanResp_hPtGenNorm"));
  hPtGenNorm->SetYTitle("1 / N  dN / dp^{gen}_{T}  1 / (GeV)");
  normHist(hPtGenNorm,"width");

  // ----- Fit profiles -----
  std::vector<TH1*> hResp(2);
  hResp[0] = new TH1D("MeanResp_hResp",";p^{gen}_{T} (GeV);< p^{jet}_{T} / p^{gen}_{T} >",
		      hRespVsPtGen->GetNbinsX(),
		      hRespVsPtGen->GetXaxis()->GetXbins()->GetArray());
  hResp[1] = new TH1D("MeanResp_hRespLog",";p^{gen}_{T} (GeV);< p^{jet}_{T} / p^{gen}_{T} >",
		      hRespVsPtGenLog->GetNbinsX(),
		      hRespVsPtGenLog->GetXaxis()->GetXbins()->GetArray());

  std::vector<TH1*> hReso(2);
  hReso[0] = new TH1D("MeanResp_hReso",
		      ";p^{gen}_{T} (GeV);#sigma(p^{jet}_{T} / p^{gen}_{T}) / < p^{jet}_{T} / p^{gen}_{T} >",
		      hRespVsPtGen->GetNbinsX(),
		      hRespVsPtGen->GetXaxis()->GetXbins()->GetArray());
  hReso[1] = new TH1D("MeanResp_hResoLog",
		      ";p^{gen}_{T} (GeV);#sigma(p^{jet}_{T} / p^{gen}_{T}) / < p^{jet}_{T} / p^{gen}_{T} >",
		      hRespVsPtGenLog->GetNbinsX(),
		      hRespVsPtGenLog->GetXaxis()->GetXbins()->GetArray());

  // Get 1D slices and get mean, sigma
  for(size_t i = 0; i < 2; ++i) {
    hResp[i]->SetMarkerStyle(20);
    hResp[i]->SetLineWidth(2);
    hResp[i]->SetNdivisions(505);
    hResp[i]->Sumw2();
    
    hReso[i]->SetMarkerStyle(20);
    hReso[i]->SetLineWidth(2);
    hReso[i]->SetNdivisions(505);
    hReso[i]->Sumw2();

    TH2* h = hRespVsPtGen;
    if( i == 1 ) h = hRespVsPtGenLog;
    TH1* hSlice = new TH1D("hSlice","",h->GetNbinsY(),h->GetYaxis()->GetXmin(),
			   h->GetYaxis()->GetXmax());
    hSlice->Sumw2();
    for(int xBin = 1; xBin <= h->GetNbinsX(); xBin++) {
      hSlice->Reset();
      for(int yBin = 1; yBin <= h->GetNbinsY(); yBin++) {
	hSlice->SetBinContent(yBin,h->GetBinContent(h->GetBin(xBin,yBin)));
	hSlice->SetBinError(yBin,h->GetBinError(h->GetBin(xBin,yBin)));
      }  

      double mean       = hSlice->GetMean();
      double width      = hSlice->GetRMS();
      //      if( width < 0.1 ) width = 0.1;
      hSlice->Fit("gaus","QNLI0","",mean-1.5*width,mean+1.5*width);
      TF1 *f = static_cast<TF1*>(gROOT->GetFunction("gaus"));
      mean = f->GetParameter(1);
      double meanError = f->GetParError(1);
      width = f->GetParameter(2);
      double widthError = f->GetParError(2);
      
      hResp[i]->SetBinContent(xBin,mean);
      hResp[i]->SetBinError(xBin,meanError);
      if( mean ) {
	hReso[i]->SetBinContent(xBin,width/mean);
	hReso[i]->SetBinError(xBin,widthError/mean);
      }
    }
    delete hSlice;
    hResp[i]->GetYaxis()->SetRangeUser(0.9,1.2);
    hReso[i]->GetYaxis()->SetRangeUser(0.,0.4);
  }


  // ----- Fit resolution -----
  std::vector<TF1*> fReso(2);
  fReso[0] = new TF1("MeanResp_fitReso","sqrt([0]*[0]/x/x + [1]*[1]/x + [2]*[2])",
		     ptBinEdges_.front(),ptBinEdges_.back());
  fReso[1] = static_cast<TF1*>(fReso[0]->Clone("MeanResp_fitResoLog"));
  std::vector<TPaveText*> fitLabel(2);

  for(size_t i = 0; i < 2; ++i) {
    fReso[i]->SetLineWidth(2);
    fReso[i]->SetLineColor(2);
    fReso[i]->FixParameter(0,0.);
    fReso[i]->SetParameter(1,1.3);
    fReso[i]->SetParameter(2,0.04);
    hReso[i]->Fit(fReso[i],"ILQ0R");

    fitLabel[i] = createPaveText(2);
    char txt[50];
    sprintf(txt,"a_{0} = %.2f #pm %.2f #sqrt{GeV}",fReso[i]->GetParameter(0),fReso[i]->GetParError(0));
    fitLabel[i]->AddText(txt);
    sprintf(txt,"a_{1} = %.2f #pm %.2f #sqrt{GeV}",fReso[i]->GetParameter(1),fReso[i]->GetParError(1));
    fitLabel[i]->AddText(txt);
    sprintf(txt,"a_{2} = %.3f #pm %.3f GeV",fReso[i]->GetParameter(2),fReso[i]->GetParError(2));
    fitLabel[i]->AddText(txt);    
  }


  // ----- Fitting spectrum -----
  TF1 *fitSpectrum = new TF1("fitSpectrum",spectrum,ptBinsMin(),ptBinsMax(),8);
  fitSpectrum->FixParameter(0,ptBinsMin());
  fitSpectrum->FixParameter(1,ptBinsMax());
  fitSpectrum->SetParameter(2,3.85);
  fitSpectrum->SetParameter(3,0.03);
  fitSpectrum->SetParameter(4,6.73);
  fitSpectrum->SetParameter(5,0.02);
  fitSpectrum->SetParameter(6,9.31);
  fitSpectrum->SetParameter(7,0.009);
  fitSpectrum->SetLineColor(2);
  hPtGenNorm->Fit(fitSpectrum,"0ILR");

  std::cout << "Int(" << ptBinsMin() << "," << ptBinsMax() << ") = ";
  std::cout << fitSpectrum->Integral(ptBinsMin(),ptBinsMax()) << std::endl;
  std::cout << std::endl;
  for(int i = 0; i < fitSpectrum->GetNpar(); i++) {
    std::cout << "$" << i-2 << "$ & $" << std::flush;
    std::cout << fitSpectrum->GetParameter(i) << " \\pm " << std::flush;
    std::cout << fitSpectrum->GetParError(i) << "$ \\\\" << std::endl;
  }

  TPaveText *sampleLabel = createPaveText(3,0.7);
  sampleLabel->AddText("CMS QCD simulation");
  sampleLabel->AddText("#sqrt{s} = 7 TeV,  L = 10 pb^{-1}");
  sampleLabel->AddText("Anti-k_{T} R = 0.5 dijets");

  TLegend *fitLeg = createLegend(1,0.7,lineHeight(),3.5*lineHeight());
  fitLeg->AddEntry(fitSpectrum,"Fit: N #sum^{2}_{i=0} e^{(-a_{2i}-a_{2i+1}#upoint p_{T})}","L");
  

  // ----- Plot histograms -----
  TPostScript * const ps = new TPostScript((dir_+"/jsMeanRespAndReso.ps").c_str(),111);
  TCanvas * c1 = new TCanvas("c1","MeanRespAndReso",0,0,600,600);

  drawPSPage(ps,c1,hPtGen,"PE1");

  std::vector<TObject*> objs;
  objs.push_back(hPtGenNorm);
  objs.push_back(fitSpectrum);
  objs.push_back(sampleLabel);
  objs.push_back(fitLeg);
  drawPSPage(ps,c1,objs,"PE1",true);

  for(size_t i = 0; i < 2; ++i) {
    drawPSPage(ps,c1,hResp[i],"PE1");
    
    ps->NewPage();
    c1->cd();
    if( i == 1 ) gPad->SetLogx();
    hReso[i]->Draw("PE1");
    fReso[i]->Draw("same");
    fitLabel[i]->Draw("same");
    gPad->RedrawAxis();
    c1->Draw();
  }
  
  TFile rootfile((dir_+"/jsResponse.root").c_str(),"UPDATE");
  rootfile.WriteTObject(hResp[0]);
  rootfile.WriteTObject(hResp[1]);
  rootfile.WriteTObject(hReso[0]);
  rootfile.WriteTObject(hReso[1]);
  rootfile.WriteTObject(hPtGen);


  // ----- Clean up -----
  for(size_t i = 0; i < 2; ++i) {
    delete hResp[i];
    delete hReso[i];
    delete fReso[i];
    delete fitLabel[i];
  }
  delete hPtGen;
  delete hPtGenNorm;
  delete fitSpectrum;
  delete sampleLabel;
  delete fitLeg;
  delete c1;
  delete ps;
}


// --------------------------------------------------
double ControlPlotsJetSmearing::spectrum(double *x, double *par) {
  double min = par[0];
  double max = par[1];
  double f = 0.;
  if( min < x[0] && x[0] < max ) {
    double norm = exp(-par[2])*(exp(-par[3]*min)-exp(-par[3]*max))/par[3];
    norm += exp(-par[4])*(exp(-par[5]*min)-exp(-par[5]*max))/par[5];
    norm += exp(-par[6])*(exp(-par[7]*min)-exp(-par[7]*max))/par[7];
    f = (exp(-par[2] - par[3]*x[0]) + exp(-par[4] - par[5]*x[0]) + exp(-par[6] - par[7]*x[0]))/norm;
  }
  return f;
}



// --------------------------------------------------
void ControlPlotsJetSmearing::plot3rdJet() const {
  std::cout << "Creating 3rd jet control plots\n";

  // Binning in pt3rel
  std::vector<double> pt3BinEdges;
  pt3BinEdges.push_back(0.);
  pt3BinEdges.push_back(0.06);
  pt3BinEdges.push_back(0.08);
  pt3BinEdges.push_back(0.10);
  pt3BinEdges.push_back(0.12);
  pt3BinEdges.push_back(0.15);
  pt3BinEdges.push_back(0.20);
  int nPt3Bins = static_cast<int>(pt3BinEdges.size()-1);

  // ----- Create histograms -----

  // -- Relative pt of 2nd and 3rd jet --
  // Pt of (i+1)th jet relative to ith jet
  TH1 *hPtGen21 = new TH1D("hPtGen21",
			   ";p^{gen}_{T,2} / p^{gen}_{T,1};Number of events",
			   50,0.,1.5);
  TH1 *hPtGen31 = new TH1D("hPtGen31",
			   ";p^{gen}_{T,3} / p^{gen}_{T,1};Number of events",
			   25,0.,0.3);
  TH1 *hRelJet3Pt = new TH1D("hRelJet3Pt",
			     ";p^{rel}_{T,3};Number of events",
			     25,0.,0.3);
  TH2 *hRelJet3PtVsPtGen31 = new TH2D("hRelJet3PtVsPtGen31",
				      ";p^{rel}_{T,3};p^{gen}_{T,3} / p^{gen}_{T,1}",
				      50,0.,0.2,50,0,0.2);
  hRelJet3PtVsPtGen31->SetNdivisions(505,"XY");
  TGraphErrors *gDeltaPtGen21Width = 0;
  TH1 *hDeltaPtGen21WidthFrame = new TH1D("hDeltaPtGen21WidthFrame",
					  ";p^{gen}_{T,1} (GeV);#sigma / p^{gen}_{T,1}",
					  100,ptBinsMin(),ptBinsMax());
  TF1 *fDeltaPtGen21Width = new TF1("fDeltaPtGen21Width","[0]/sqrt(x)",ptBinsMin(),ptBinsMax());
  fDeltaPtGen21Width->SetParameter(0,1.);
  fDeltaPtGen21Width->SetLineWidth(2);
  fDeltaPtGen21Width->SetLineColor(2);

  TH1 *hAngularScale2 = new TH1D("hAngularScale2",";(Scale2)  /  (p^{gen}_{T,2} / p^{gen}_{T,1})",50,0,2);
  TH1 *hAngularScale3 = static_cast<TH1D*>(hAngularScale2->Clone("hAngularScale3"));
  hAngularScale3->SetXTitle("(Scale3)  /  (p^{gen}_{T,3} / p^{gen}_{T,1})");

  TH2 *hAngularScale2vsGenMHT = new TH2D("hAngularScale2vsGenMHT",
					 ";GenJet MHT / HT;(Scale2)  /  (p^{gen}_{T,2} / p^{gen}_{T,1})",
					 500,0,0.25,500,0,2);
  hAngularScale2vsGenMHT->SetNdivisions(505,"X");
  TH2 *hAngularScale3vsGenMHT = static_cast<TH2D*>(hAngularScale2vsGenMHT->Clone("hAngularScale3vsGenMHT"));
  hAngularScale3vsGenMHT->SetYTitle("(Scale3)  /  (p^{gen}_{T,3} / p^{gen}_{T,1})");

  TH1 *hNAllVsPtRec = new TH1D("hNAllVsPtRec",";p^{reco,ave}_{T} (GeV);N(all)",
			       nPtBins(),&(ptBinEdges_.front()));
  TH1 *hNMisorderedVsPtRec = static_cast<TH1D*>(hNAllVsPtRec->Clone("hNMisorderedVsPtRec"));
  hNMisorderedVsPtRec->SetYTitle("N(incorrect)");
  TH1 *hOrderingVsPtRec = static_cast<TH1D*>(hNAllVsPtRec->Clone("hOrderingVsPtRec"));
  hOrderingVsPtRec->SetYTitle("N(incorrect) / N(all)");
  hOrderingVsPtRec->GetYaxis()->SetRangeUser(0.,1.);
  TGraphAsymmErrors *gOrderingVsPtRec = new TGraphAsymmErrors();
  gOrderingVsPtRec->SetMarkerStyle(20);
  TH1 *hSameGenJetVsPtRec = static_cast<TH1D*>(hOrderingVsPtRec->Clone("hSameGenJetVsPtRec"));
  hSameGenJetVsPtRec->SetYTitle("N(same GenJet) / N(all)");
  std::vector<TH1*> hNSameGenJetVsPtRec;
  std::vector<TGraphAsymmErrors*> gSameGenJetVsPtRec;
  for(int i = 0; i < 3; ++i) {
    int color = i>0 ? 2*i : 1;
    TString name = "hNSameGenJetVsPtRec";
    name += i;
    hNSameGenJetVsPtRec.push_back(static_cast<TH1D*>(hNMisorderedVsPtRec->Clone(name)));
    hNSameGenJetVsPtRec.back()->SetYTitle("N(same GenJet)");
    hNSameGenJetVsPtRec.back()->SetMarkerStyle(20+i);
    hNSameGenJetVsPtRec.back()->SetMarkerColor(color);
    hNSameGenJetVsPtRec.back()->SetLineColor(color);
    gSameGenJetVsPtRec.push_back(new TGraphAsymmErrors());
    gSameGenJetVsPtRec.back()->SetMarkerStyle(20+i);
    gSameGenJetVsPtRec.back()->SetMarkerColor(color);
    gSameGenJetVsPtRec.back()->SetLineColor(color);
  }

  TH1 *hNAllVsPtGen = static_cast<TH1D*>(hNAllVsPtRec->Clone("hNAllVsPtGen"));
  hNAllVsPtGen->SetXTitle("p^{gen,ave}_{T} (GeV)");
  TH1 *hNMisorderedVsPtGen = static_cast<TH1D*>(hNAllVsPtGen->Clone("hNMisorderedVsPtGen"));
  hNMisorderedVsPtGen->SetYTitle("N(incorrect)");
  TH1 *hOrderingVsPtGen = static_cast<TH1D*>(hOrderingVsPtRec->Clone("hOrderingVsPtGen"));
  hOrderingVsPtGen->SetXTitle("p^{gen,ave}_{T} (GeV)");
  TGraphAsymmErrors *gOrderingVsPtGen = static_cast<TGraphAsymmErrors*>(gOrderingVsPtRec->Clone());
  TH1 *hSameGenJetVsPtGen = static_cast<TH1D*>(hOrderingVsPtGen->Clone("hSameGenJetVsPtGen"));
  hSameGenJetVsPtGen->SetYTitle("N(same GenJet) / N(all)");
  std::vector<TH1*> hNSameGenJetVsPtGen;
  std::vector<TGraphAsymmErrors*> gSameGenJetVsPtGen;
  for(int i = 0; i < 3; ++i) {
    int color = i>0 ? 2*i : 1;
    TString name = "hNSameGenJetVsPtGen";
    name += i;
    hNSameGenJetVsPtGen.push_back(static_cast<TH1D*>(hNMisorderedVsPtGen->Clone(name)));
    hNSameGenJetVsPtGen.back()->SetYTitle("N(same GenJet) / N(all)");
    hNSameGenJetVsPtGen.back()->SetMarkerStyle(20+i);
    hNSameGenJetVsPtGen.back()->SetMarkerColor(color);
    hNSameGenJetVsPtGen.back()->SetLineColor(color);
    gSameGenJetVsPtGen.push_back(new TGraphAsymmErrors());
    gSameGenJetVsPtGen.back()->SetMarkerStyle(20);
    gSameGenJetVsPtGen.back()->SetMarkerColor(color);
    gSameGenJetVsPtGen.back()->SetLineColor(color);
  }

  // Pt of (i+1)th jet relative to ith jet
  // in different pt bins
  std::vector<TH1*> hPtGen21Bins(nPtBins());
  std::vector<TH1*> hDeltaPtGen21(nPtBins());
  std::vector<TF1*> fDeltaPtGen21(nPtBins());
  std::vector<TH1*> hPtGen31Bins(nPtBins());
  std::vector<TH1*> hRelJet3PtBins(nPtBins());
  std::vector<TH2*> hRelJet3PtVsPtGen31Bins(nPtBins());
  std::vector<TH2*> hDeltaPtGenVsRelJet3Pt(nPtBins());
  std::vector< std::vector<TH1*> > hRelJet3Bins(nPtBins());
  std::vector<TH1*> hMeanRelJet3Pt(nPt3Bins);
  for(int ptBin = 0; ptBin < nPtBins(); ++ptBin) {
    TString name = "hPtGen21Bin";
    name += ptBin;
    hPtGen21Bins[ptBin] = static_cast<TH1D*>(hPtGen21->Clone(name));
    hPtGen21Bins[ptBin]->SetYTitle("Relative number of events");
    hPtGen21Bins[ptBin]->Sumw2();

    name = "hDeltaPtGen21";
    name += ptBin;
    hDeltaPtGen21[ptBin] = new TH1D(name,";p^{gen}_{T,2} - p^{gen}_{T,1} (GeV);Number of events;",
				    25,-0.4*ptBinEdges_[ptBin],0.4*ptBinEdges_[ptBin+1]);

    name = "fDeltaPtGen21";
    name += ptBin;
    fDeltaPtGen21[ptBin] = new TF1(name,"gaus",-0.3*ptBinEdges_[ptBin],0.3*ptBinEdges_[ptBin+1]);
    fDeltaPtGen21[ptBin]->SetLineWidth(2);
    fDeltaPtGen21[ptBin]->SetLineColor(2);

    name = "hPtGen31Bin";
    name += ptBin;
    hPtGen31Bins[ptBin] = static_cast<TH1D*>(hPtGen31->Clone(name));
    hPtGen31Bins[ptBin]->SetYTitle("Relative number of events");
    hPtGen31Bins[ptBin]->Sumw2();
    
    name = "hRelJet3PtBin";
    name += ptBin;
    hRelJet3PtBins[ptBin] = static_cast<TH1D*>(hRelJet3Pt->Clone(name));
    hRelJet3PtBins[ptBin]->SetYTitle("Relative number of events");
    hRelJet3PtBins[ptBin]->Sumw2();

    name = "hRelJet3PtVsPtGen31Bins";
    name += ptBin;
    hRelJet3PtVsPtGen31Bins[ptBin] = static_cast<TH2D*>(hRelJet3PtVsPtGen31->Clone(name));

    name = "hDeltaPtGenVsRelJet3Pt";
    name += ptBin;
    hDeltaPtGenVsRelJet3Pt[ptBin] = new TH2D(name,";p^{rel}_{T};p^{gen}_{T,1} - p^{gen}_{T,2} (GeV)",
					     60,0.,0.4,50,-0.4*ptBinEdges_.at(ptBin),
					     0.4*ptBinEdges_.at(ptBin+1));

    hRelJet3Bins[ptBin] = std::vector<TH1*>(nPt3Bins);
    for(int pt3Bin = 0; pt3Bin < nPt3Bins; ++pt3Bin) {
      name = "hRelJet3Bins_";
      name += "Pt";
      name += ptBin;
      name += "_Pt3";
      name += pt3Bin;
      hRelJet3Bins[ptBin][pt3Bin] = new TH1D(name,";p^{rel}_{T,3};Number of events",
					     50,pt3BinEdges.front(),pt3BinEdges.back());
      hRelJet3Bins[ptBin][pt3Bin]->SetMarkerStyle(20);
      hRelJet3Bins[ptBin][pt3Bin]->SetMarkerColor(color(pt3Bin));
      hRelJet3Bins[ptBin][pt3Bin]->SetLineColor(color(pt3Bin));
      hRelJet3Bins[ptBin][pt3Bin]->SetNdivisions(505);
    }
  }
  for(int pt3Bin = 0; pt3Bin < nPt3Bins; ++pt3Bin) {
    TString name = "hMeanRelJet3PtBin";
    name += pt3Bin;
    hMeanRelJet3Pt[pt3Bin] = new TH1D(name,";p_{T} (GeV);< p^{rel}_{T,3} >",
				      nPtBins(),&(ptBinEdges_.front()));
    hMeanRelJet3Pt[pt3Bin]->SetMarkerStyle(20);
    double deltaPt3 = pt3BinEdges[pt3Bin+1] - pt3BinEdges[pt3Bin];
    hMeanRelJet3Pt[pt3Bin]->GetYaxis()->SetRangeUser(pt3BinEdges[pt3Bin]-0.2*deltaPt3,
						     pt3BinEdges[pt3Bin+1]+0.5*deltaPt3);
  }


  // -- Sigma depending on 3rd jet pt --

  // Response and mean sigma histograms in bins of pt
  std::vector<TH1*> hResMC(nPtBins());		// MC truth response from ptgen
  std::vector< std::vector<TH1*> > hResMean(nPt3Bins); // Response from mean pttrue

  TH1 *hSigmaMC = new TH1D("hSigmaMC",
			   ";p_{T} (GeV);#sigma(p_{T}) / p_{T}",
			   nPtBins(),&(ptBinEdges_.front()));
  hSigmaMC->SetNdivisions(505);
  hSigmaMC->SetMarkerStyle(20);
  std::vector<TH1*> hSigmaMean(nPt3Bins);

  TH1 *hSigmaVsPt3Frame	= new TH1D("hSigmaVsPt3Frame",
				   ";p^{rel}_{T,3};#sigma(p_{T}) / p_{T}",
				   500,0.,1.5*pt3BinEdges.back());
  hSigmaVsPt3Frame->SetNdivisions(505);
  hSigmaVsPt3Frame->GetYaxis()->SetRangeUser(0.,0.3);
  std::vector<TH1*> hSigmaMCVsPt3(nPtBins());
  std::vector<TGraphErrors*> gSigmaVsPt3(nPtBins());
  std::vector<TF1*> fSigmaVsPt3(nPtBins());

  for(int pt3Bin = 0; pt3Bin < nPt3Bins; ++pt3Bin) { // Loop over rel pt3 bins
    hResMean[pt3Bin] = std::vector<TH1*>(nPtBins());

    TString name = "hSigmaMean";
    name += pt3Bin;
    hSigmaMean[pt3Bin] = static_cast<TH1D*>(hSigmaMC->Clone(name));
    hSigmaMean[pt3Bin]->SetMarkerColor(color(pt3Bin));
    hSigmaMean[pt3Bin]->SetLineColor(color(pt3Bin));

    for(int ptBin = 0; ptBin < nPtBins(); ++ptBin) { // Loop over ptBins
      if( pt3Bin == 0 ) {
	name = "hResMC";
	name += ptBin;
	hResMC[ptBin] = new TH1D(name,"",51,0.,2.);

	name = "hSigmaMCVsPt3";
	name += ptBin;
	hSigmaMCVsPt3[ptBin] = static_cast<TH1D*>(hSigmaVsPt3Frame->Clone(name));
	hSigmaMCVsPt3[ptBin]->SetLineColor(4);
	hSigmaMCVsPt3[ptBin]->SetLineStyle(2);

	name = "fSigmaVsPt3";
	name += ptBin;
	fSigmaVsPt3[ptBin] = new TF1(name,"pol1",0.,1.5*pt3BinEdges.back());
	fSigmaVsPt3[ptBin]->SetLineColor(2);
	fSigmaVsPt3[ptBin]->SetLineWidth(2);
      }
      name = "hResMean";
      name += pt3Bin;
      name += "_";
      name += ptBin;
      hResMean[pt3Bin][ptBin] = static_cast<TH1D*>(hResMC[ptBin]->Clone(name));
    } //End of loop over ptBins
  } // End of loop over pt3Bins

  TH1 *hGaussWidth = new TH1D("hGaussWidth3rdJet",";p_{T} (GeV);#sigma(p_{T}) / p_{T}",
			      500,ptBinsMin(),ptBinsMax());
  hGaussWidth->SetMarkerStyle(1);
  hGaussWidth->SetLineWidth(2);
  hGaussWidth->SetLineStyle(2);
  TH1 *hGaussWidthTruth = static_cast<TH1D*>(hGaussWidth->Clone("hGaussWidthTruth3rdJet"));
  hGaussWidthTruth->SetLineStyle(1);


  // ----- Fill histograms -----
  for(DataIt datait = data_->begin(); datait != data_->end(); datait++) {
    // Select DiJet events
    if( (*datait)->type() == TypeSmearDiJet )  {
      SmearDiJet * dijet = static_cast<SmearDiJet*>(*datait);  

      const Jet *j1 = dijet->jet1();
      const Jet *j2 = dijet->jet2();
      const Jet *j3 = dijet->jet3();

      hAngularScale2->Fill(j1->genPt()*dijet->scalePt2()/j2->genPt());
      hAngularScale3->Fill(j1->genPt()*dijet->scalePt3()/j3->genPt());
      hAngularScale2vsGenMHT->Fill(dijet->relGenMet(),j1->genPt()*dijet->scalePt2()/j2->genPt());
      hAngularScale3vsGenMHT->Fill(dijet->relGenMet(),j1->genPt()*dijet->scalePt3()/j3->genPt());

      hNAllVsPtRec->Fill(dijet->avePt());
      hNAllVsPtGen->Fill(dijet->avePtGen());
      double eps = 1E-5;
      if( j3->genPt() > j1->genPt()+eps || j3->genPt() > j2->genPt()+eps ) {
	hNMisorderedVsPtRec->Fill(dijet->avePt());
	hNMisorderedVsPtGen->Fill(dijet->avePtGen());
      }
      if( j2->genPt() > j1->genPt()-eps && j2->genPt() < j1->genPt()+eps ) {
	hNSameGenJetVsPtRec[0]->Fill(dijet->avePt());
	hNSameGenJetVsPtGen[0]->Fill(dijet->avePtGen());
      }
      if( j3->genPt() > j1->genPt()-eps && j3->genPt() < j1->genPt()+eps ) {
	hNSameGenJetVsPtRec[1]->Fill(dijet->avePt());
	hNSameGenJetVsPtGen[1]->Fill(dijet->avePtGen());
      }
      if( j3->genPt() > j2->genPt()-eps && j3->genPt() < j2->genPt()+eps ) {
	hNSameGenJetVsPtRec[2]->Fill(dijet->avePt());
	hNSameGenJetVsPtGen[2]->Fill(dijet->avePtGen());
      }


//       if( rand_->Uniform() > 0.5 ) {
// 	j1 = dijet->jet2();
// 	j2 = dijet->jet1();
//       }

      hPtGen21->Fill(j2->genPt()/j1->genPt(),dijet->weight());
      hPtGen31->Fill(j3->genPt()/j1->genPt(),dijet->weight());
      hRelJet3Pt->Fill(dijet->relJet3Pt(),dijet->weight());
      hRelJet3PtVsPtGen31->Fill(j3->genPt()/j1->genPt(),dijet->relJet3Pt(),dijet->weight());

      int ptBin = findPtBin(j1->genPt());
      if( ptBin >= 0 ) {
	hPtGen21Bins[ptBin]->Fill(j2->genPt()/j1->genPt(),dijet->weight());
	hDeltaPtGen21[ptBin]->Fill(j2->genPt() - j1->genPt(),dijet->weight());
	hPtGen31Bins[ptBin]->Fill(j3->genPt()/j1->genPt(),dijet->weight());
	hResMC[ptBin]->Fill(j1->pt()/j1->genPt(),dijet->weight());
	hResMC[ptBin]->Fill(j2->pt()/j2->genPt(),dijet->weight());
	hRelJet3PtBins[ptBin]->Fill(dijet->relJet3Pt(),dijet->weight());
	hRelJet3PtVsPtGen31Bins[ptBin]->Fill(j3->genPt()/j1->genPt(),dijet->relJet3Pt(),dijet->weight());
	hDeltaPtGenVsRelJet3Pt[ptBin]->Fill(dijet->relJet3Pt(),
					    dijet->jet1()->genPt() - dijet->jet2()->genPt(),
					    dijet->weight());

	int pt3Bin = findBin(dijet->relJet3Pt(),pt3BinEdges);
	if( pt3Bin >= 0 && pt3Bin < nPt3Bins ) {
	  hRelJet3Bins[ptBin][pt3Bin]->Fill(dijet->relJet3Pt(),dijet->weight());
	}
      }
      double meanPtGen = 0.5*(j1->genPt()+j2->genPt());
      ptBin = findPtBin(meanPtGen);
      if( ptBin >=0 ) {
	int pt3bin = findBin(dijet->relJet3Pt(),pt3BinEdges);
	if( pt3bin >= 0 ) {
	  hResMean[pt3bin][ptBin]->Fill(j1->pt()/meanPtGen,dijet->weight());
	  hResMean[pt3bin][ptBin]->Fill(j2->pt()/meanPtGen,dijet->weight());
	}
      }
    }
  } // End of loop over data
  gOrderingVsPtRec->BayesDivide(hNMisorderedVsPtRec,hNAllVsPtRec);
  for(int i = 0; i < gOrderingVsPtRec->GetN(); ++i) {
    gOrderingVsPtRec->SetPointEXlow(i,0.);
    gOrderingVsPtRec->SetPointEXhigh(i,0.);
  }
  gOrderingVsPtGen->BayesDivide(hNMisorderedVsPtGen,hNAllVsPtGen);
  for(int i = 0; i < gOrderingVsPtGen->GetN(); ++i) {
    gOrderingVsPtGen->SetPointEXlow(i,0.);
    gOrderingVsPtGen->SetPointEXhigh(i,0.);
  }
  for(size_t k = 0; k < gSameGenJetVsPtRec.size(); ++k) {
    gSameGenJetVsPtRec[k]->BayesDivide(hNSameGenJetVsPtRec[k],hNAllVsPtRec);
    gSameGenJetVsPtGen[k]->BayesDivide(hNSameGenJetVsPtGen[k],hNAllVsPtGen);
    for(int i = 0; i < gSameGenJetVsPtGen[k]->GetN(); ++i) {
      gSameGenJetVsPtRec[k]->SetPointEXlow(i,0.);
      gSameGenJetVsPtRec[k]->SetPointEXhigh(i,0.);
      gSameGenJetVsPtGen[k]->SetPointEXlow(i,0.);
      gSameGenJetVsPtGen[k]->SetPointEXhigh(i,0.);
    }
  }

  // Fit Gaussian width per bin
  for(int pt3bin = 0; pt3bin < nPt3Bins; ++pt3bin) {
    for(int ptBin = 0; ptBin < nPtBins(); ++ptBin) {
      if( pt3bin == 0 ) {
	hResMC[ptBin]->Fit("gaus","IL0Q");
	TF1 *fit = hResMC[ptBin]->GetFunction("gaus");
	hSigmaMC->SetBinContent(1+ptBin,fit->GetParameter(2));
	hSigmaMC->SetBinError(1+ptBin,fit->GetParError(2));
      }
      hResMean[pt3bin][ptBin]->Fit("gaus","IL0Q");
      TF1 *fit = hResMean[pt3bin][ptBin]->GetFunction("gaus");
      hSigmaMean[pt3bin]->SetBinContent(1+ptBin,fit->GetParameter(2));
      hSigmaMean[pt3bin]->SetBinError(1+ptBin,fit->GetParError(2));
    }
  }

  // Fill histogram of true and fitted Gaussian width
  if( parClass_ == "SmearParametrizationGauss" ) {
    for(int bin = 1; bin <= hGaussWidth->GetNbinsX(); ++bin) {
      double pt = hGaussWidth->GetBinCenter(bin);

      double sigma = gaussianWidth(pt) / pt;
      hGaussWidth->SetBinContent(bin,sigma);
      hGaussWidth->SetBinError(bin,0.);
      
      sigma = gaussianWidthTruth(pt) / pt;
      hGaussWidthTruth->SetBinContent(bin,sigma);
      hGaussWidthTruth->SetBinError(bin,0.);
    }
  }

  // Create graph of sigma vs pt3rel and fit with linear function
  std::vector<double> xVal(nPt3Bins);
  std::vector<double> yVal(nPt3Bins);
  std::vector<double> xErr(nPt3Bins,0.);
  std::vector<double> yErr(nPt3Bins);
  for(int ptBin = 0; ptBin < nPtBins(); ++ptBin) {
    for(int pt3bin = 0; pt3bin < nPt3Bins; ++pt3bin) {
      xVal[pt3bin] = 0.5*(pt3BinEdges[pt3bin]+pt3BinEdges[pt3bin+1]);
      yVal[pt3bin] = hSigmaMean[pt3bin]->GetBinContent(ptBin+1);
      yErr[pt3bin] = hSigmaMean[pt3bin]->GetBinError(ptBin+1);
    }
    gSigmaVsPt3[ptBin] = new TGraphErrors(xVal.size(),&(xVal.front()),&(yVal.front()),
					  &(xErr.front()),&(yErr.front()));
    gSigmaVsPt3[ptBin]->Fit(fSigmaVsPt3[ptBin],"0QR");

    for(int bin = 1; bin <= hSigmaMCVsPt3[ptBin]->GetNbinsX(); ++bin) {
      hSigmaMCVsPt3[ptBin]->SetBinContent(bin,hSigmaMC->GetBinContent(ptBin+1));
    }
  }

  // Fit ptGen2-ptGen1 distributions
  xVal.clear();
  yVal.clear();
  xErr.clear();
  yErr.clear();
  for(int ptBin = 0; ptBin < nPtBins(); ++ptBin) {
    hDeltaPtGen21[ptBin]->Fit(fDeltaPtGen21[ptBin],"IR0Q");

    xVal.push_back(0.5*(ptBinEdges_[ptBin]+ptBinEdges_[ptBin+1]));
    yVal.push_back(fDeltaPtGen21[ptBin]->GetParameter(2)/xVal[ptBin]);
    xErr.push_back(0.);
    yErr.push_back(fDeltaPtGen21[ptBin]->GetParError(2)/xVal[ptBin]);
  }
  gDeltaPtGen21Width = new TGraphErrors(xVal.size(),&(xVal.front()),&(yVal.front()),
					&(xErr.front()),&(yErr.front()));
  gDeltaPtGen21Width->Fit(fDeltaPtGen21Width,"R0Q");
  double maxWidthDeltaPtGen21 = *(std::max_element(yVal.begin(),yVal.end()));

  // Mean rel 3rd jet pt
  for(int pt3Bin = 0; pt3Bin < nPt3Bins; ++pt3Bin) {
    for(int ptBin = 0; ptBin < nPtBins(); ++ptBin) {
      hMeanRelJet3Pt[pt3Bin]->SetBinContent(ptBin+1,hRelJet3Bins[ptBin][pt3Bin]->GetMean());
      hMeanRelJet3Pt[pt3Bin]->SetBinError(ptBin+1,hRelJet3Bins[ptBin][pt3Bin]->GetMeanError());
    }
  }


  // ----- Edit histograms -----
  for(int ptBin = 0; ptBin < nPtBins(); ++ptBin) {
    normHist(hPtGen21Bins[ptBin]);
    normHist(hPtGen31Bins[ptBin]);
    normHist(hRelJet3PtBins[ptBin]);

    setYRange(hPtGen21Bins[ptBin],0.,1.4);
    setYRange(hDeltaPtGen21[ptBin],0.,1.4);
    setYRange(hPtGen31Bins[ptBin],0.,1.4);
    setYRange(hRelJet3PtBins[ptBin],0.,1.4);
    setYRange(hSigmaMC,0.5,2.);
    double max = hRelJet3Bins[ptBin][0]->GetMaximum();
    for(int pt3Bin = 1; pt3Bin < nPt3Bins; ++pt3Bin) {
      if( hRelJet3Bins[ptBin][pt3Bin]->GetMaximum() > max ) {
	max = hRelJet3Bins[ptBin][pt3Bin]->GetMaximum();
      }      
    }
    for(int pt3Bin = 0; pt3Bin < nPt3Bins; ++pt3Bin) {
      hRelJet3Bins[ptBin][pt3Bin]->GetYaxis()->SetRangeUser(0.,1.5*max);
    }
    
  }
  hDeltaPtGen21WidthFrame->GetYaxis()->SetRangeUser(0.,1.5*maxWidthDeltaPtGen21);
  


  // ----- Plot histograms -----
  TPostScript * const ps = new TPostScript((dir_+"/js3rdJet.ps").c_str(),111);
  TCanvas * c1 = new TCanvas("c1","3rd Jet",0,0,600,600);
  std::vector<TObject*> objs;

  ps->NewPage();
  c1->cd();
  hNMisorderedVsPtGen->Draw("PE1");
  gPad->SetLogx(1);
  gPad->SetLogy(1);
  c1->Draw();
  if( saveAsEps_ ) c1->SaveAs((outNamePrefix_+"_NMisorderedVsPtGen.eps").c_str(),"eps");

  ps->NewPage();
  c1->cd();
  hNMisorderedVsPtRec->Draw("PE1");
  gPad->SetLogx(1);
  gPad->SetLogy(1);
  if( saveAsEps_ ) c1->SaveAs((outNamePrefix_+"_NMisorderedVsPtRec.eps").c_str(),"eps");
  c1->Draw();

  ps->NewPage();
  c1->cd();
  gPad->SetLogx(1);
  gPad->SetLogy(0);
  hOrderingVsPtGen->Draw("H");
  gOrderingVsPtGen->Draw("PE1same");
  if( saveAsEps_ ) c1->SaveAs((outNamePrefix_+"_OrderingVsPtGen.eps").c_str(),"eps");
  c1->Draw();

  ps->NewPage();
  c1->cd();
  gPad->SetLogx(1);
  gPad->SetLogy(0);
  hOrderingVsPtRec->Draw("H");
  gOrderingVsPtRec->Draw("PE1same");
  if( saveAsEps_ ) c1->SaveAs((outNamePrefix_+"_OrderingVsPtRec.eps").c_str(),"eps");
  c1->Draw();

  int maxIdx = 0;
  if( hNSameGenJetVsPtGen[1]->GetMaximum() > hNSameGenJetVsPtGen[maxIdx]->GetMaximum() ) maxIdx = 1;
  if( hNSameGenJetVsPtGen[2]->GetMaximum() > hNSameGenJetVsPtGen[maxIdx]->GetMaximum() ) maxIdx = 2;
  TLegend *leg = createLegend(3);
  leg->AddEntry(hNSameGenJetVsPtGen[0],"GenJet2 = GenJet1","P");
  leg->AddEntry(hNSameGenJetVsPtGen[1],"GenJet3 = GenJet1","P");
  leg->AddEntry(hNSameGenJetVsPtGen[2],"GenJet3 = GenJet2","P");
  ps->NewPage();
  c1->cd();
  gPad->SetLogx(1);
  gPad->SetLogy(1);
  hNSameGenJetVsPtGen[maxIdx]->GetYaxis()->SetRangeUser(0.5,100.*hNSameGenJetVsPtGen[maxIdx]->GetMaximum());
  hNSameGenJetVsPtGen[maxIdx]->Draw("PE1");
  for(size_t k = 0; k < hNSameGenJetVsPtGen.size(); ++k) {
    hNSameGenJetVsPtGen[k]->Draw("PE1same");
  }
  leg->Draw("same");
  if( saveAsEps_ ) c1->SaveAs((outNamePrefix_+"_NSameGenJetVsPtGen.eps").c_str(),"eps");
  c1->Draw();

  maxIdx = 0;
  if( hNSameGenJetVsPtRec[1]->GetMaximum() > hNSameGenJetVsPtRec[maxIdx]->GetMaximum() ) maxIdx = 1;
  if( hNSameGenJetVsPtRec[2]->GetMaximum() > hNSameGenJetVsPtRec[maxIdx]->GetMaximum() ) maxIdx = 2;
  ps->NewPage();
  c1->cd();
  gPad->SetLogx(1);
  gPad->SetLogy(1);
  hNSameGenJetVsPtRec[maxIdx]->GetYaxis()->SetRangeUser(0.5,100.*hNSameGenJetVsPtRec[maxIdx]->GetMaximum());
  hNSameGenJetVsPtRec[maxIdx]->Draw("PE1");
  for(size_t k = 0; k < hNSameGenJetVsPtRec.size(); ++k) {
    hNSameGenJetVsPtRec[k]->Draw("PE1same");
  }
  leg->Draw("same");
  if( saveAsEps_ ) c1->SaveAs((outNamePrefix_+"_NSameGenJetVsPtRec.eps").c_str(),"eps");
  c1->Draw();

  ps->NewPage();
  c1->cd();
  gPad->SetLogx(1);
  gPad->SetLogy(0);
  hSameGenJetVsPtGen->Draw();
  for(size_t k = 0; k < gSameGenJetVsPtGen.size(); ++k) {
    gSameGenJetVsPtGen[k]->Draw("PE1same");
  }
  leg->Draw("same");
  if( saveAsEps_ ) c1->SaveAs((outNamePrefix_+"_SameGenJetVsPtGen.eps").c_str(),"eps");
  c1->Draw();

  ps->NewPage();
  c1->cd();
  gPad->SetLogx(1);
  hSameGenJetVsPtRec->Draw();
  for(size_t k = 0; k < gSameGenJetVsPtRec.size(); ++k) {
    gSameGenJetVsPtRec[k]->Draw("PE1same");
  }
  leg->Draw("same");
  if( saveAsEps_ ) c1->SaveAs((outNamePrefix_+"_SameGenJetVsPtRec.eps").c_str(),"eps");
  c1->Draw();

  gPad->SetLogx(0);
  gPad->SetLogy(0);

  drawPSPage(ps,c1,hAngularScale2);
  drawPSPage(ps,c1,hAngularScale3);
  drawPSPage(ps,c1,hAngularScale2vsGenMHT,"COLZ");
  drawPSPage(ps,c1,hAngularScale3vsGenMHT,"COLZ");

  std::vector<TPaveText*> ptBinLabels(nPtBins());
  for(int ptBin = 0; ptBin < nPtBins(); ++ptBin) {
    ptBinLabels[ptBin] = createPaveText(1,0.8);
    char txt[50];
    sprintf(txt,"%.0f < p^{gen}_{T,1} < %.0f GeV",ptBinEdges_[ptBin],ptBinEdges_[ptBin+1]);
    ptBinLabels[ptBin]->AddText(txt);
  }

  // Mean 3rd jet pt per pt3rel bin
  for(int ptBin = 0; ptBin < nPtBins(); ++ptBin) {
    objs.clear();
    for(int pt3Bin = 0; pt3Bin < nPt3Bins; ++pt3Bin) {
      objs.push_back(hRelJet3Bins[ptBin][pt3Bin]);
    }
    objs.push_back(ptBinLabels[ptBin]);
    drawPSPage(ps,c1,objs,"PE1");
    TString name = outNamePrefix_;
    name += "_Rel3rdJetPt_PtBin";
    name += ptBin;
    name += ".eps";
    if( saveAsEps_ ) c1->SaveAs(name,"eps");
  }
  for(int pt3Bin = 0; pt3Bin < nPt3Bins; ++pt3Bin) {
    hMeanRelJet3Pt[pt3Bin]->Fit("pol0","0Q");
    TF1 *fitMean = hMeanRelJet3Pt[pt3Bin]->GetFunction("pol0");
    fitMean->SetLineColor(2);
    fitMean->SetLineWidth(2);
    
    std::cout << "  " << pt3BinEdges[pt3Bin] << " < pt3rel < " << pt3BinEdges[pt3Bin+1] << ":  " << std::flush;
    std::cout << fitMean->GetParameter(0) << " +/- " << fitMean->GetParError(0) << std::endl;

    TPaveText *txt = createPaveText(2);
    char entry[50];
    sprintf(entry,"%.2f < p^{rel}_{T,3} < %.2f GeV",pt3BinEdges[pt3Bin],pt3BinEdges[pt3Bin+1]);
    txt->AddText(entry);
    sprintf(entry,"Fit: < p^{rel}_{T,3} > = %.4f #pm %.4f",fitMean->GetParameter(0),fitMean->GetParError(0));
    txt->AddText(entry);

    objs.clear();
    objs.push_back(hMeanRelJet3Pt[pt3Bin]);
    objs.push_back(fitMean);
    objs.push_back(txt);
    drawPSPage(ps,c1,objs,"PE1");
    TString name = outNamePrefix_;
    name += "_MeanRel3rdJetPt_Pt3Bin";
    name += pt3Bin;
    name += ".eps";
    if( saveAsEps_ ) c1->SaveAs(name,"eps");
    delete fitMean;
  }

  drawPSPage(ps,c1,hPtGen21,"PE1");
  drawPSPage(ps,c1,hPtGen31,"PE1");
  //  drawPSPage(ps,c1,hPtGen21,"PE1",true);
  //  drawPSPage(ps,c1,hPtGen31,"PE1",true);
  drawPSPage(ps,c1,hRelJet3Pt,"PE1");
  drawPSPage(ps,c1,hRelJet3PtVsPtGen31,"COLZ",true);

  for(int ptBin = 0; ptBin < nPtBins(); ++ptBin) {
    objs.clear();
    objs.push_back(hDeltaPtGenVsRelJet3Pt[ptBin]);
    objs.push_back(ptBinLabels[ptBin]);
    drawPSPage(ps,c1,objs,"COLZ",true);
    TString name = outNamePrefix_;
    name += "_DeltaPtGenVsJet3Pt_PtBin";
    name += ptBin;
    name += ".eps";
    //if( saveAsEps_ ) c1->SaveAs(name,"eps");
  }
  for(int ptBin = 0; ptBin < nPtBins(); ++ptBin) {
    objs.clear();
    objs.push_back(hPtGen21Bins[ptBin]);
    objs.push_back(ptBinLabels[ptBin]);
    drawPSPage(ps,c1,objs,"PE1");
  }
  for(int ptBin = 0; ptBin < nPtBins(); ++ptBin) {
    objs.clear();
    objs.push_back(hPtGen31Bins[ptBin]);
    objs.push_back(ptBinLabels[ptBin]);
    drawPSPage(ps,c1,objs,"PE1");
  }
  for(int ptBin = 0; ptBin < nPtBins(); ++ptBin) {
    objs.clear();
    objs.push_back(hDeltaPtGen21[ptBin]);
    objs.push_back(fDeltaPtGen21[ptBin]);
    objs.push_back(ptBinLabels[ptBin]);
    drawPSPage(ps,c1,objs,"PE1");
  }

  TLegend *legDeltaPtGen21Width = createLegend(1);
  char txt[50];
  sprintf(txt,"Fit: (%.3f #pm %0.3f) #sqrt{GeV} / #sqrt{p^{gen}_{T,1}}",
	  fDeltaPtGen21Width->GetParameter(0),fDeltaPtGen21Width->GetParError(0));
  legDeltaPtGen21Width->AddEntry(fDeltaPtGen21Width,txt,"L");
  objs.clear();
  objs.push_back(hDeltaPtGen21WidthFrame);
  objs.push_back(gDeltaPtGen21Width);
  objs.push_back(fDeltaPtGen21Width);
  objs.push_back(legDeltaPtGen21Width);
  drawPSPage(ps,c1,objs,"PE1");

  for(int ptBin = 0; ptBin < nPtBins(); ++ptBin) {
    objs.clear();
    objs.push_back(hRelJet3PtBins[ptBin]);
    objs.push_back(ptBinLabels[ptBin]);
    drawPSPage(ps,c1,objs,"PE1");
  }
  for(int ptBin = 0; ptBin < nPtBins(); ++ptBin) {
    objs.clear();
    objs.push_back(hRelJet3PtVsPtGen31Bins[ptBin]);
    objs.push_back(ptBinLabels[ptBin]);
    drawPSPage(ps,c1,objs,"COLZ",true);
  }

  for(int ptBin = 0; ptBin < nPtBins(); ++ptBin) {
    ps->NewPage();
    c1->cd();
    hSigmaVsPt3Frame->Draw();
    hSigmaMCVsPt3[ptBin]->Draw("Hsame");
    fSigmaVsPt3[ptBin]->Draw("same");
    gSigmaVsPt3[ptBin]->Draw("PE1same");
    ptBinLabels[ptBin]->Draw("same");
    gPad->RedrawAxis();
    c1->Draw();
  }

  ps->NewPage();
  c1->cd();
  delete leg;
  leg = createLegend(nPt3Bins+1,0.6);
  leg->AddEntry(hSigmaMC,"MC truth","P");
  hSigmaMC->Draw("PE1");
  hGaussWidthTruth->Draw("Hsame");
  hGaussWidth->Draw("Hsame");
  for(int pt3bin = 0; pt3bin < nPt3Bins; ++pt3bin) {
    char txt[50];
    sprintf(txt,"%.2f < p^{rel}_{T,3} < %.2f",pt3BinEdges[pt3bin],pt3BinEdges[pt3bin+1]);
    leg->AddEntry(hSigmaMean[pt3bin],txt,"P");
    hSigmaMean[pt3bin]->Draw("PE1same");
  }
  leg->Draw("same");
  gPad->RedrawAxis();
  c1->Draw();

  // ----- Write histos to file -----
  TFile rootfile((dir_+"/jsResponse.root").c_str(),"UPDATE");
  for(int ptBin = 0; ptBin < nPtBins(); ++ptBin) {
    rootfile.WriteTObject(hPtGen21Bins[ptBin]);
    rootfile.WriteTObject(hDeltaPtGen21[ptBin]);
    rootfile.WriteTObject(fDeltaPtGen21[ptBin]);
    rootfile.WriteTObject(hPtGen31Bins[ptBin]);
    rootfile.WriteTObject(hRelJet3PtBins[ptBin]);
    rootfile.WriteTObject(hRelJet3PtVsPtGen31Bins[ptBin]);
    rootfile.WriteTObject(gSigmaVsPt3[ptBin]);
    rootfile.WriteTObject(fSigmaVsPt3[ptBin]);
    rootfile.WriteTObject(hDeltaPtGenVsRelJet3Pt[ptBin]);
  }
  rootfile.WriteTObject(hPtGen21);
  rootfile.WriteTObject(hPtGen31);
  rootfile.WriteTObject(hRelJet3PtVsPtGen31);
  rootfile.WriteTObject(gDeltaPtGen21Width);
  rootfile.WriteTObject(fDeltaPtGen21Width);


  // ----- Clean up -----
  delete hPtGen21;
  delete hPtGen31;
  delete hRelJet3PtVsPtGen31;
  delete hAngularScale2;
  delete hAngularScale3;
  delete hAngularScale2vsGenMHT;
  delete hAngularScale3vsGenMHT;
  for(int ptBin = 0; ptBin < nPtBins(); ++ptBin) {
    delete hPtGen21Bins[ptBin];
    delete hDeltaPtGen21[ptBin];
    delete fDeltaPtGen21[ptBin];
    delete hPtGen31Bins[ptBin];
    delete hRelJet3PtVsPtGen31Bins[ptBin];
    delete hResMC[ptBin];
    delete ptBinLabels[ptBin];
    delete hSigmaMCVsPt3[ptBin];
    delete gSigmaVsPt3[ptBin];
    delete fSigmaVsPt3[ptBin];
    delete hDeltaPtGenVsRelJet3Pt[ptBin];
    for(int pt3Bin = 0; pt3Bin < nPt3Bins; ++pt3Bin) {
      delete hRelJet3Bins[ptBin][pt3Bin];
    }
  }
  delete hSigmaVsPt3Frame;
  delete hSigmaMC;
  for(int pt3Bin = 0; pt3Bin < nPt3Bins; ++pt3Bin) {
    delete hSigmaMean[pt3Bin];
    delete hMeanRelJet3Pt[pt3Bin];
    for(int ptBin = 0; ptBin < nPtBins(); ++ptBin) {
      delete hResMean[pt3Bin][ptBin];
    }
  }
  delete hGaussWidth;
  delete hGaussWidthTruth;
  delete hDeltaPtGen21WidthFrame;
  delete gDeltaPtGen21Width;
  delete fDeltaPtGen21Width;
  delete legDeltaPtGen21Width;

  delete hNAllVsPtRec;
  delete hNAllVsPtGen;
  delete hNMisorderedVsPtRec;
  delete hNMisorderedVsPtGen;
  delete hOrderingVsPtRec;
  delete hOrderingVsPtGen;
  delete gOrderingVsPtRec;
  delete gOrderingVsPtGen;
  delete hSameGenJetVsPtRec;
  delete hSameGenJetVsPtGen;
  for(int k = 0; k < 3; ++k) {
    delete hNSameGenJetVsPtRec[k];
    delete hNSameGenJetVsPtGen[k];
    delete gSameGenJetVsPtRec[k];
    delete gSameGenJetVsPtGen[k];
  }

  delete c1;
  delete ps;
}


// --------------------------------------------------
double ControlPlotsJetSmearing::gaussianWidth(double pt) const {
  double a[3];
  for(int i = 0; i < 3; i++) {
    a[i] = scale_[i]*param_->GetPars()[i];
  }
  return sqrt( a[0]*a[0] + a[1]*a[1]*pt + a[2]*a[2]*pt*pt );
}


// ------------------------------------------------------------------------
double ControlPlotsJetSmearing::gaussianWidthError(double pt) const {
  // Calculate derivatives
  std::vector<double> ds(3);
  double s = gaussianWidth(pt);
  for(int i = 0; i < 3; i++) {
    ds[i] = scale_[i]*param_->GetPars()[i]/s;
    if( i == 1 ) ds[i] *= pt;
    if( i == 2 ) ds[i] *= pt*pt;
  }

  // Calculate variance
  double var = 0.;
  for(int i = 0; i < 3; i++) { // Outer loop over parameters
    for(int j = 0; j < i+1; j++) { // Inner loop over parameters
      int idx = (i*i + i)/2 + j; // Index of (i,j) in covariance vector
      if( i == j ) { // Diagonal terms
	var += ds[i]*ds[i]*scale_[i]*scale_[i]*param_->GetCovCoeff()[idx];
      } else { // Off-diagonal terms
	var += 2*ds[i]*ds[j]*scale_[i]*scale_[j]*param_->GetCovCoeff()[idx];
      }
    } // End of inner loop over parameters
  } // End of outer loop over parameters

  return sqrt(var);
}


// --------------------------------------------------
double ControlPlotsJetSmearing::gaussianWidthTruth(double pt) const {
  double s = 0.;
  if( truthPar_.size() >= 3 ) {
    s = sqrt( truthPar_[0]*truthPar_[0]
	      +truthPar_[1]*truthPar_[1]*pt
	      +truthPar_[2]*truthPar_[2]*pt*pt );
  }

  return s;
}


//!  \brief Draw TObject on one page of a ps file
//!
//!  Opens a new page in the PostScript file \p ps
//!  and draws the TObject \p obj on it.
//!
//!  \param ps     A new page is added to this file
//!                containing the TObject \p obj
//!  \param can    The \p obj is drawn on this canvas
//!  \param obj    The TObject to be drawn
//!  \param option Draw options
//!  \param log    Sets log-scale on last axis if true
// --------------------------------------------------
void ControlPlotsJetSmearing::drawPSPage(TPostScript * ps, TCanvas * can, TObject * obj, const std::string &option, bool log) const {
  std::vector<TObject*> objs;
  objs.push_back(obj);
  drawPSPage(ps,can,objs,option,log);
}



//!  \brief Draw TObjects on one page of a ps file
//!
//!  Opens a new page in the PostScript file \p ps
//!  and draws all TObjects in \p obs on it. If \p objs
//!  contains more than one TObject, the draw option
//!  "same" is automatically used.
//!
//!  \param ps     A new page is added to this file
//!                containing all TObjects in \p objs
//!  \param can    The \p objs are drawn on this canvas
//!  \param objs   Contains the objects that are drawn
//!  \param option Draw options (except for "same",
//!                see above)
//!  \param log    Sets log-scale on last axis if true
// --------------------------------------------------
void ControlPlotsJetSmearing::drawPSPage(TPostScript * ps, TCanvas * can, std::vector<TObject*> objs, const std::string &option, bool log) const {
  ps->NewPage();
  can->cd();
  for( size_t i = 0; i < objs.size(); i++ ) {
    std::string opt = option;
    if( i ) opt.append("same");
    objs.at(i)->Draw(opt.c_str());
  }
  //gPad->RedrawAxis();
  std::string hist = objs.at(0)->ClassName();
  if( hist.find("TH1") != std::string::npos ) {
    if( log ) can->SetLogy(1);
    else      can->SetLogy(0);
  } else if( hist.find("TH2") != std::string::npos ) {
    can->SetLogy(0);
    if( log ) can->SetLogz(1);
    else      can->SetLogz(0);
  }

  can->Draw();
}



//!  \brief Find y-axis range
//!
//!  Sets \p min and \p max to the minimum (non-zero) and
//!  maximum bin content of \p h, respectively.
// --------------------------------------------------
void ControlPlotsJetSmearing::findYRange(const TH1 * h, double& min, double& max) const {
  min = 10000.;
  max = 0.;
  for(int bin = 1; bin <= h->GetNbinsX(); bin++) {
    double val = h->GetBinContent(bin);
    if( val > 0. && val < min ) min = val;
    if( val > 0. && val > max ) max = val;
  }
  if( min > max ) {
    min = 1E-3;
    max = 1;
  }
}



//!  \brief Set default \p gStyle options
// --------------------------------------------------
void ControlPlotsJetSmearing::setGStyle() const
{
  gStyle->SetErrorX(0);
  gStyle->SetPalette(1);

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
  //  gStyle->SetMarkerStyle(20);

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
  gStyle->SetPadTopMargin(0.16);
  gStyle->SetPadBottomMargin(0.18);
  gStyle->SetPadLeftMargin(0.18);
  gStyle->SetPadRightMargin(0.16);

  // For the Global title:
  gStyle->SetOptTitle(1);
  gStyle->SetTitleFont(42,"");
  gStyle->SetTitleColor(1);
  gStyle->SetTitleTextColor(1);
  gStyle->SetTitleFillColor(0);
  gStyle->SetTitleFontSize(0.12);
  gStyle->SetTitleAlign(23);
  gStyle->SetTitleX(0.515);
  gStyle->SetTitleH(0.06);
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



// --------------------------------------------------
TLegend *ControlPlotsJetSmearing::createLegend(int nEntries, double width, double lineHgt, double yOffset) const {
  double margin = 0.04;
  double x0 = gStyle->GetPadLeftMargin()+margin;
  double x1 = 1.-(gStyle->GetPadRightMargin()+margin);
  x0 = x0 + (1.-width)*(x1-x0);
  double y1 = 1.-(gStyle->GetPadTopMargin()+margin+0.02+yOffset);
  double height = lineHeight();
  if( lineHgt > 0 ) height = lineHgt;
  double y0 = y1 - nEntries*height;
  TLegend *leg = new TLegend(x0,y0,x1,y1);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextFont(42);
  return leg;
}


// --------------------------------------------------
TPaveText *ControlPlotsJetSmearing::createPaveText(int nEntries, double width, double lineHgt) const {
  double margin = 0.04;
  double x0 = gStyle->GetPadLeftMargin()+margin;
  double x1 = 1.-(gStyle->GetPadRightMargin()+margin);
  x0 = x0 + (1.-width)*(x1-x0);
  double y1 = 1.-(gStyle->GetPadTopMargin()+margin+0.02);
  double height = lineHeight();
  if( lineHgt > 0 ) height = lineHgt;
  double y0 = y1 - nEntries*height;
  TPaveText *txt = new TPaveText(x0,y0,x1,y1,"NDC");
  txt->SetBorderSize(0);
  txt->SetFillColor(0);
  txt->SetTextFont(42);
  txt->SetTextAlign(12);
  return txt;
}


//!  \brief Adjust y-axis range
//!
//!  Sets the y-axis range of \p h from
//!  <tt> c1 * min</tt> to <tt> c2 * max</tt>,
//!  where \p min and \p max are the minimal and
//!  the maximal bin non-zero content, respectively.
//!  If <tt>min < minLimit</tt>, \p minLimit is used
//!  instead as minimum.
// --------------------------------------------------
void ControlPlotsJetSmearing::setYRange(TH1 * h, double c1, double c2, double minLimit) const {
  double min = 0.;
  double max = 0.;
  findYRange(h,min,max);
  min *= c1;
  max *= c2;
  if( min < minLimit ) min = minLimit;
  h->GetYaxis()->SetRangeUser( min, max );
}

void ControlPlotsJetSmearing::normHist(TH1 *h, std::string option) const { 
  if( h->Integral(option.c_str()) ) h->Scale(1./h->Integral(option.c_str())); 
}

void ControlPlotsJetSmearing::normHist(TH1 *h, double min, double max, std::string option) const { 
  double norm = h->Integral(h->FindBin(min),h->FindBin(max),option.c_str());
  if( norm ) h->Scale(1./norm);
}

//! \brief Convert to std::string
//!
//! Converts a given object to an std::string
//! using std::stringstream.
//!
//! \param t Object to be converted to a string
//! \return String representation of t
// --------------------------------------------------
template <class T> std::string ControlPlotsJetSmearing::toString(const T& t) const {
  std::stringstream ss;
  ss << t;
  return ss.str();
}
 
 
// --------------------------------------------------
int ControlPlotsJetSmearing::findPtBin(const Jet *jet) const {
  double var = 0.;
  if( ptBinningVar_ == "ptGen" ) var = jet->genPt();
  else if( ptBinningVar_ == "pt" ) var = jet->pt();

  return findBin(var,ptBinEdges_);
}


// --------------------------------------------------
int ControlPlotsJetSmearing::findPtBin(double pt) const {
  return findBin(pt,ptBinEdges_);
}
 

// --------------------------------------------------
int ControlPlotsJetSmearing::findBin(double x, const std::vector<double> &binEdges) const {
  int bin = -1; // Under- / overflow bin
  if( x >= binEdges.front() && x <= binEdges.back() ) {
    for(int i = 0; i < static_cast<int>(binEdges.size()-1); ++i) {
      if( x > binEdges[i] ) bin = i;
      else break;
    }
  }

  return bin;
}
  

//!  Filling \p bins with borders of \p nBins bins between \p first
//!  and \p last that are equidistant when viewed in log scale,
//!  so \p bins must have length \p nBins+1. If \p first, \p last
//!  or \p nBins are not positive, failure is reported.
// -------------------------------------------------------------
bool ControlPlotsJetSmearing::equidistLogBins(std::vector<double>& bins, int 
					      nBins, double first, double last) const {
  if( static_cast<int>(bins.size()) != nBins+1 ) return false;
  if( nBins < 1 || first <= 0. || last <= 0. || first >= last ) return false;

  bins[0]     = first;
  bins[nBins] = last;
  const double firstLog = log10(bins[0]);
  const double lastLog  = log10(bins[nBins]);
  for (int i = 1; i < nBins; ++i) {
    bins[i] = pow(10., firstLog + i*(lastLog-firstLog)/(nBins));
  }

  return true;
}


int ControlPlotsJetSmearing::color(int i) const {
  int               col = 1;
  if( i == 1 )      col = 2;
  else if( i == 2 ) col = 4;
  else if( i == 3 ) col = 6;
  else if( i == 4 ) col = 8;
  else if( i == 5 ) col = 9;

  return col;
}
