// $Id: ControlPlotsJetSmearing.cc,v 1.19 2010/09/22 13:29:44 mschrode Exp $

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
  std::cout << "PtBinning variable: " << ptBinningVar_ << std::endl;

  minJetPt_ = config_->read<double>("Et min cut on jet",0.);
  maxJetPt_ = config_->read<double>("Et max cut on jet",1.);

  // Override possible existing root file
  TFile rootfile((dir_+"/jsResponse.root").c_str(),"RECREATE");
  rootfile.Close();

  outNamePrefix_ = dir_+"/"+config_->read<std::string>("plots name prefix","JS")+"_";
  saveAsEps_ = config_->read<bool>("plots save as eps",false);
}

ControlPlotsJetSmearing::~ControlPlotsJetSmearing() {
  ptBinEdges_.clear();
  ptBinCenters_.clear();
  delete rand_;
}



// --------------------------------------------------
void ControlPlotsJetSmearing::makePlots() const {
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
  if( config_->read<bool>("create parallel components plots",false) )
    plotParallelComponents();
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

  TH1 * hDeltaPtJet12;
  TH1 * hDeltaPtJet12Fit;
  TH1 * hPJet3;
  TH1 * hPJet3Rel;
  TH1 * hPJet3GenRel;
  TH1 * hPSJ;
  TH1 * hPSJRel;
  TH1 * hPSJGenRel;
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
  TH1 * hPtJet1 = 0;
  TH1 * hPtJet2 = 0;
  TH1 * hPtJet3 = 0;
  TH1 * hPtJet4 = 0;
  TH2 * hPSJvsPtJet4 = 0;
  TH2 * hPtJet1vs2 = 0;
  TH1 * hEta = 0;
  TH1 * hDeltaPhi12 = 0;
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
    hPtGenAsym[ptBin] = new TH1D(name.c_str(),";p^{gen}_{T} asymmetry;",45,-0.2,0.2);
    hPtGenAsym[ptBin]->Sumw2();
    hPtGenAsym[ptBin]->SetLineWidth(2);

    name = "hPtAsym_" + toString(ptBin);
    hPtAsym[ptBin] = new TH1D(name.c_str(),";p_{T} asymmetry;",175,-1.,1.);
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

  hDeltaPtJet12 = new TH1D("hDeltaPtJet12",";0.5#upoint|p_{T,1} - p_{T,2}| (GeV);",
			   50,0.,0.3*ptBinEdges_.back());
  hDeltaPtJet12->SetMarkerStyle(20);
  hDeltaPtJet12->SetNdivisions(505);
  hDeltaPtJet12->Sumw2();

  hDeltaPtJet12Fit = new TH1D("hDeltaPtJet12Fit",";0.5#upoint|p_{T,1} - p_{T,2}| (GeV);",
			      500,0.,0.5*ptBinEdges_.back());

  hPJet3 = new TH1D("hPJet3",";p_{||,3} (GeV);Events",
		    50,-0.3*ptBinEdges_.front(),0.3*ptBinEdges_.back());
  hPJet3->SetMarkerStyle(20);
  
  hPJet3Rel = new TH1D("hPJet3Rel",";p_{||,3} / p^{ref}_{T};Events",50,-0.25,0.25);
  hPJet3Rel->SetMarkerStyle(20);

  hPJet3GenRel = static_cast<TH1D*>(hPJet3Rel->Clone("hPJet3GenRel"));
  hPJet3GenRel->SetXTitle(";p_{||,3} / p^{gen,ave}_{T}");

  
  hPSJ = new TH1D("hPSJ",";p_{||,>3} (GeV);Events",
		  50,-0.3*ptBinEdges_.front(),0.3*ptBinEdges_.back());
  hPSJ->SetMarkerStyle(20);
  
  hPSJRel = new TH1D("hPSJRel",";p_{||,>3} / p^{ref}_{T};Events",50,-0.1,0.1);
  hPSJRel->SetMarkerStyle(20);

  hPSJGenRel = static_cast<TH1D*>(hPJet3Rel->Clone("hPSJGenRel"));
  hPSJGenRel->SetXTitle(";p_{||,>3} / p^{gen,ave}_{T}");

  hPSJvsPtJet4 = new TH2D("hPSJvsPtJet4",";p_{t,4};p_{||,>3} (GeV) (GeV)",
			  50,0.,0.3*ptBinEdges_.back(),
			  50,-0.1*ptBinEdges_.front(),0.1*ptBinEdges_.back());


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
		       120,0.8*ptBinEdges_.front(),1.1*ptBinEdges_.back());
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
  hPtDijet->SetTitle(";p^{dijet}_{T} (GeV);Events");

  hPtJet1 = new TH1D("hPtJet1",";p_{T,1} (GeV);Jets",70,0.5*ptBinEdges_.front(),1.3*ptBinEdges_.back());
  hPtJet1->SetMarkerStyle(20);
  hPtJet1->GetXaxis()->SetNdivisions(505);
  hPtJet1->SetLineWidth(2);
  setColor(hPtJet1,2);

  hPtJet2 = static_cast<TH1D*>(hPtJet1->Clone("hPtJet2"));
  hPtJet2->SetXTitle("p_{T,2} (GeV)");
  setColor(hPtJet2,4);

  hPtJet3 = new TH1D("hPtJet3",";p_{T,3} (GeV);Events",50,0.,0.2*ptBinEdges_.back());

  hPtJet4 = static_cast<TH1D*>(hPtJet3->Clone("hPtJet4"));
  hPtJet4->SetXTitle("p_{T,4} (GeV)");

  hEta = new TH1D("hEta",";#eta;Jets",101,-5.1,5.1);
  hDeltaPhi12 = new TH1D("hDeltaPhi12",";|#Delta#phi_{12}|;Events",100,0.,3.2);

  hPtJet1vs2 = new TH2D("hPtJet1vs2",";p_{T,1} (GeV);p_{T,2} (GeV);",70,0.5*ptBinEdges_.front(),1.3*ptBinEdges_.back(),70,0.5*ptBinEdges_.front(),1.3*ptBinEdges_.back());
  hPtJet1vs2->GetXaxis()->SetNdivisions(505);
  hPtJet1vs2->GetYaxis()->SetNdivisions(505);

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
    if( (*datait)->type() == TypeSmearDiJet ) {
      SmearDiJet * dijet = static_cast<SmearDiJet*>(*datait);  
      double weight = dijet->weight();
      double ptAve = dijet->avePt();
      double ptHat = dijet->ptHat();
            
      // Event wise spectra
      hPtHat->Fill(ptHat,weight);
      hPtDijet->Fill(ptAve,weight);
      hPJet3->Fill(dijet->pJ3(),weight);
      hPJet3Rel->Fill(dijet->pJ3()/dijet->ptRef(),weight);
      hPSJ->Fill(dijet->pSJ()/weight);
      hPSJRel->Fill(dijet->pSJ()/dijet->ptRef(),weight);
      hPSJvsPtJet4->Fill(dijet->ptJet4(),dijet->pSJ(),weight);
      double ptGenAve = 0.5*(dijet->jet1()->genPt()+dijet->jet2()->genPt());
      if( ptGenAve > 0 ) {
	hPJet3GenRel->Fill(dijet->pJ3()/ptGenAve,weight);
	hPSJGenRel->Fill(dijet->pSJ()/ptGenAve,weight);
      }
      hPtGenJet1->Fill(dijet->jet1()->genPt(),weight);
      hPtJet1->Fill(dijet->jet1()->pt(),weight);
      hPtGenJet2->Fill(dijet->jet2()->genPt(),weight);
      hPtJet2->Fill(dijet->jet2()->pt(),weight);
      hPtJet3->Fill(dijet->ptJet3(),weight);
      hPtJet4->Fill(dijet->ptJet4(),weight);
      hDeltaPhi12->Fill(dijet->deltaPhi12(),weight);
      hDeltaPtJet12->Fill(0.5*std::abs(dijet->jet1()->pt()-dijet->jet2()->pt()),weight);
      
      // Pt asymmetry variables
      const Jet *j1 = dijet->jet1();
      const Jet *j2 = dijet->jet2();
      hPtJet1vs2->Fill(j1->pt(),j2->pt(),weight);
      double ht = j1->pt() + j2->pt();
      if( ht > 0 ) {
	double ptAsym = (j1->pt() - j2->pt())/ht;
	hPtAsymBiased[0]->Fill(ptAsym,weight);
	hPtAsym[0]->Fill(ptAsym,weight);
	hPtAsym[0]->Fill(-1.*ptAsym,weight);
      }
      double ptGenAsym = j1->genPt() - j2->genPt();
      ht = j1->genPt() + j2->genPt();
      if( ht > 0 ) {
	ptGenAsym /= ht;
      } else {
	ptGenAsym = -2.;
      }
      
      for(int i = 0; i < 2; i++) {        // Loop over both jets
	const Jet * jet = dijet->jet(i);
	
	// Jet wise spectra
	hPtGenAbs->Fill(jet->genPt(),weight);
	hPtGen->Fill(jet->genPt(),weight);
	hEta->Fill(jet->eta(),weight);
	
	// Response distributions
	double var = 0.;
	if( ptBinningVar_ == "ptGen" )        var = jet->genPt();
	else if( ptBinningVar_ == "pt" )      var = jet->pt();
	else if( ptBinningVar_ == "ptDijet" ) var = ptAve;
	int bin = findPtBin(var);
	if( bin >= 0 && bin < nPtBins() ) {
	  if( jet->genPt() > 0 ) {
	    hRespMeasAbs[bin]->Fill( jet->pt() / jet->genPt(),weight);
	    hRespMeas[bin]->Fill( jet->pt() / jet->genPt(),weight);
	    if( i == 0 ) {
	      hRespMeasJet1[bin]->Fill( jet->pt() / jet->genPt(),weight);
	    } else if ( i == 1 ) {
	      hRespMeasJet2[bin]->Fill( jet->pt() / jet->genPt(),weight);
	    }
	  }
	  hPtGenAbsBins[bin]->Fill(jet->genPt(),weight);
	  hRespMCPtHat[bin]->Fill(jet->pt()/ptHat,weight);
	  hPtGenAsym[bin]->Fill(ptGenAsym,weight);
	  hPtGenAsym[bin]->Fill(-1.*ptGenAsym,weight);
	}
      }
    }
  } // End of loop over data
  for(int ptBin = 0; ptBin < nPtBins(); ptBin++) {
    normHist(hRespMeas[ptBin],"width");
    normHist(hRespMeasJet1[ptBin],"width");
    normHist(hRespMeasJet2[ptBin],"width");
    normHist(hRespMCPtHat[ptBin],"width");
    normHist(hPtAsym[ptBin],"width");
    normHist(hPtAsymBiased[ptBin],"width");
  }
  normHist(hDeltaPtJet12,"width");
  hDeltaPtJet12->Scale(0.5);
  normHist(hPtGen,tMin,tMax,"width");
  normHist(hPtGenJet1,tMin,tMax,"width");
  normHist(hPtGenJet2,tMin,tMax,"width");
  normHist(hPtHat,tMin,tMax,"width");


  // --- Fill histograms of fitted response ----------------
  // Get parameters
  std::vector<double> fittedPar(param_->GetNumberOfParameters());
  for(int i = 0; i < param_->GetNumberOfParameters(); i++) {
    fittedPar.at(i) = param_->GetPars()[i];
  }

  std::vector<double> auxPar = bag_of<double>(config_->read<string>("mean response parameters","1 0"));
  SmearData * smearData = dynamic_cast<SmearData*>(data_->front());
  if( smearData ) {
    for(int bin = 1; bin <= hDeltaPtJet12Fit->GetNbinsX(); bin++) {
        double s = param_->GetPars()[0]/sqrt(2.);
        double d = hDeltaPtJet12Fit->GetBinCenter(bin);
        if( d < 2.*s ) {
	  double u = d/s;
	  double norm = s*sqrt(M_PI*2.);//*erf(1./sqrt(2.));
	  hDeltaPtJet12Fit->SetBinContent(bin,exp(-0.5*u*u)/norm);
 	}
//       SmearDiJet *dijet = dynamic_cast<SmearDiJet*>(smearData);
//       if( dijet ) {
//  	double pdf = dijet->pdfPtMeas(hDeltaPtJet12Fit->GetBinCenter(bin),0.,dijet->ptRef());
//    	hDeltaPtJet12Fit->SetBinContent(bin,pdf);
//       }
    }

    // Loop over ptBins
    for(int ptBin = 0; ptBin < nPtBins(); ptBin++) {
      // Interpolated response function
      double ptMean = 1.;
      if( ptBinningVar_ == "ptGen" ) ptMean = hPtGenAbsBins[ptBin]->GetMean();
      else if( ptBinningVar_ == "ptDijet" ) ptMean = hPtDijet->GetMean();
      for(int bin = 1; bin <= hRespFit[ptBin]->GetNbinsX(); bin++) {
	double r = hRespFit[ptBin]->GetBinCenter(bin);
	double val = smearData->pdfResp(r,ptMean);
	hRespFit[ptBin]->SetBinContent(bin,val);
	hRespFitErrStat[ptBin]->SetBinContent(bin,val);
	hRespFitErrStat[ptBin]->SetBinError(bin,smearData->pdfRespError(r,ptMean));

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
	hFitPtAsym[ptBin]->SetBinContent(bin,smearData->pdfDijetAsym(a,ptMean));
      }

      // Interpolated fit function with start values
      // Copy start values into parameter array
      for(int i = 0; i < param_->GetNumberOfParameters(); i++) {
	param_->GetPars()[i] = startParJet_.at(i);
      }
      // Plot response function
      for(int bin = 1; bin <= hRespFitStart[ptBin]->GetNbinsX(); bin++) {
	double r = hRespFitStart[ptBin]->GetBinCenter(bin);
	hRespFitStart[ptBin]->SetBinContent(bin,smearData->pdfResp(r,ptMean));
      }
      // Copy back fitted values into parameter array
      for(int i = 0; i < param_->GetNumberOfParameters(); i++) {
	param_->GetPars()[i] = fittedPar.at(i);
      }
    } // End of loop over ptBins
  } // End if( smearData )


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

    setYRange(hRespMeasAbs[ptBin],0.5,50.);
  }
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
    if( ptBinningVar_ == "ptDijet" ) binVar = "p^{ave}_{T}";
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
    hDeltaPtJet12->Draw("PE1");
    hDeltaPtJet12Fit->Draw("Lsame");
    TF1 *func = new TF1("func","gaus",0.,2.*hDeltaPtJet12->GetRMS());
    func->SetParameter(0,1./sqrt(2.*M_PI)/hDeltaPtJet12->GetRMS());
    func->FixParameter(1,0.);
    func->SetParameter(2,hDeltaPtJet12->GetRMS());
    hDeltaPtJet12->Fit(func,"I0BR");
    std::cout << "\nHistogram Fit 0.5|x1 - x2|: " << std::abs(func->GetParameter(2))*sqrt(2.) << " +/- " << func->GetParError(2)*sqrt(2.) << std::endl;
    func->SetLineColor(2);
    func->SetLineWidth(1);
    func->SetLineStyle(2);
    func->Draw("same");
    gPad->RedrawAxis();
    c1->Draw();

    ps->NewPage();
    c1->cd();
    hDeltaPtJet12->Draw("PE1");
    hDeltaPtJet12Fit->Draw("Lsame");
    func->SetRange(0,hDeltaPtJet12->GetXaxis()->GetBinUpEdge(hDeltaPtJet12->GetNbinsX()));
    func->Draw("same");
    std::vector<TLine*> lines;
    for(int k = 0; k < 6; ++k) {
      double x = (2.+k)*std::abs(func->GetParameter(2));//*sqrt(2.);
      if( x > hDeltaPtJet12->GetXaxis()->GetBinUpEdge(hDeltaPtJet12->GetNbinsX()) ) break;
      TLine *line = new TLine(x,3E-5,x,hDeltaPtJet12->GetMaximum());
      line->SetLineStyle(2);
      line->SetLineWidth(1);
      line->SetLineColor(4);
      line->Draw("same");
      lines.push_back(line);
    }
    gPad->RedrawAxis();
    c1->SetLogy();
    c1->Draw();

    for(size_t k = 0; k < lines.size(); ++k) {
      delete lines[k];
    }      

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
    hRespFitBins[ptBin]->Draw("HISTsame");
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
    hRespFitBins[ptBin]->Draw("HISTsame");
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
    hRespFitBins[ptBin]->Draw("HISTsame");
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
      fitResp[jetIdx] = new TF1(fitName,"gaus",0.,2.);
      hMC->Fit(fitResp[jetIdx],"Q0IR");
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

    //     ps->NewPage();
    //     c1->cd();
    //     hRespMCPtHat[ptBin]->Draw();
    //     hRespFit[ptBin]->Draw("Lsame");
    //     gPad->RedrawAxis();
    //     legPtRangeAndCenters[ptBin]->Draw("same");
    //     c1->SetLogy(logy);
    //     c1->Draw();
    ps->NewPage();
    c1->cd();
    hPtGenAsym[ptBin]->Draw("PE1");
    legPtRange[ptBin]->Draw("same");
    c1->SetLogy(logy);
    c1->Draw();


    ps->NewPage();
    c1->cd();
    //TF1 *fitAsym = new TF1("fitAsym","gaus",-0.5,0.5);
    TF1 *fitAsym = new TF1("fitAsym","gaus",hPtAsym[ptBin]->GetMean()-1.8*hPtAsym[ptBin]->GetRMS(),hPtAsym[ptBin]->GetMean()+1.8*hPtAsym[ptBin]->GetRMS());
    fitAsym->SetParameter(0,1./sqrt(2.*M_PI)/0.1);
    fitAsym->SetParameter(1,0.);
    fitAsym->SetParameter(2,0.1);
    hPtAsym[ptBin]->Fit(fitAsym,"0QIRB");
    fitAsym->SetLineStyle(2);
    fitAsym->SetLineWidth(2);
    hPtAsym[ptBin]->GetXaxis()->SetRangeUser(-0.2,0.2);
    hPtAsym[ptBin]->Draw("PE1");
    fitAsym->Draw("same");
    hFitPtAsym[ptBin]->Draw("Lsame");
    legPtRangeAndCenters[ptBin]->Draw("same");
    gPad->RedrawAxis();
    c1->SetLogy(logy);
    c1->Draw();

    double ptMean = 1.;
    if( ptBinningVar_ == "ptGen" ) ptMean = hPtGenAbsBins[ptBin]->GetMean();
    else if( ptBinningVar_ == "ptDijet" ) ptMean = hPtDijet->GetMean();
    double par = scale(0)*param_->GetPars()[0]/ptMean;
    double parErr = scale(0)*param_->GetErrors()[0]/ptMean;
    double parMCTruth = 0.5*(fitResp[0]->GetParameter(2)+fitResp[1]->GetParameter(2));
    double parMCTruthErr = 0.5*(fitResp[0]->GetParError(2)+fitResp[1]->GetParError(2));

    std::cout << std::endl;
    std::cout << "Pt          = " << ptMean << std::endl;
    std::cout << "Resp        = " << parMCTruth << " +/- " << parMCTruthErr << std::endl;
    std::cout << "Fit (Resp)  = " << par << " +/- " << parErr << std::endl;
    std::cout << "Dev         = " << par / parMCTruth - 1. << std::endl;
    std::cout << "Asym        = " << std::abs(fitAsym->GetParameter(2)) << " +/- " << fitAsym->GetParError(2) << std::endl;
    std::cout << "Fit (Asym)  = " << par/sqrt(2.) << " +/- " << parErr/sqrt(2.) << std::endl;
    std::cout << "Dev         = " << par/sqrt(2.)/std::abs(fitAsym->GetParameter(2)) - 1. << std::endl;

    for(int jetIdx = 0; jetIdx < 2; ++jetIdx) {
      delete fitResp[jetIdx];
      delete fitRespLine[jetIdx];
    }
    

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

  ps->NewPage();
  c1->cd();
  hPtJet1->Draw("PE1");
  hPtJet2->Draw("PE1same");
  c1->Draw();
 
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
  rootfile.WriteTObject(hDeltaPtJet12Fit);
  rootfile.WriteTObject(hDeltaPtJet12);
  rootfile.WriteTObject(hPJet3);
  rootfile.WriteTObject(hPJet3Rel);
  rootfile.WriteTObject(hPJet3GenRel);
  rootfile.WriteTObject(hPSJ);
  rootfile.WriteTObject(hPSJRel);
  rootfile.WriteTObject(hPSJGenRel);
  rootfile.WriteTObject(hPtGenAbs);
  rootfile.WriteTObject(hPtGen);
  rootfile.WriteTObject(hPtGenJet1);
  rootfile.WriteTObject(hPtGenJet2);
  rootfile.WriteTObject(hPtHat);
  rootfile.WriteTObject(hPtDijet);
  rootfile.WriteTObject(hPtJet1);
  rootfile.WriteTObject(hPtJet2);
  rootfile.WriteTObject(hPtJet3);
  rootfile.WriteTObject(hPtJet4);
  rootfile.WriteTObject(hPSJvsPtJet4);
  rootfile.WriteTObject(hPtJet1vs2);
  rootfile.WriteTObject(hEta);
  rootfile.WriteTObject(hDeltaPhi12);
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
  delete hDeltaPtJet12;
  delete hDeltaPtJet12Fit;
  delete hPJet3;
  delete hPJet3Rel;
  delete hPJet3GenRel;
  delete hPSJ;
  delete hPSJRel;
  delete hPSJGenRel;
  delete legFitStart;
  delete hPtGenAbs;
  delete hPtGen;
  delete hPtGenJet1;
  delete hPtGenJet2;
  delete hPtHat;
  delete hPtDijet;
  delete hPtJet1;
  delete hPtJet2;
  delete hPtJet3;
  delete hPtJet4;
  delete hPSJvsPtJet4;
  delete hPtJet1vs2;
  delete hEta;
  delete hDeltaPhi12;
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
void ControlPlotsJetSmearing::plotParallelComponents() const {
  std::cout << "Plotting parallel components" << std::endl;

  std::cout << "  Creating histograms" << std::endl;
  std::vector<TH1*> hPtGenAve(nPtBins());
  std::vector< std::vector<TH1*> > hPtGen(nPtBins());
  std::vector< std::vector<TH2*> > hPpOverPtGenVsDeltaPhi12(nPtBins());
  std::vector<TH1*> hDeltaPhi(nPtBins());
  std::vector<TH1*> hPPhi(nPtBins());

  std::vector<TH1*> hResp(nPtBins());
  std::vector<TH1*> hPtAsym(nPtBins());
  std::vector<TH1*> hPtGenAsym(nPtBins());
  std::vector<TH1*> hPJet3(nPtBins());
  std::vector<TH1*> hPSJ(nPtBins());
  std::vector<TH1*> hPUCE(nPtBins());
  std::vector<TH1*> hPJet3Rel(nPtBins());
  std::vector<TH1*> hPSJRel(nPtBins());
  std::vector<TH1*> hPUCERel(nPtBins());

  std::vector<double> pt3Limits;
  pt3Limits.push_back(0.06);
  pt3Limits.push_back(0.1);
  pt3Limits.push_back(0.15);
  std::vector< std::vector<TH1*> > hPtAsymPt3Cuts(nPtBins());
  std::vector< std::vector<TH1*> > hPtAsymPp3Cuts(nPtBins());

  for(int ptBin = 0; ptBin < nPtBins(); ++ptBin) {
    TString binLabelGen = toString(ptBinEdges_[ptBin])+" < p^{gen,ave}_{T} < "+toString(ptBinEdges_[ptBin+1])+" GeV";

    TString name = "ParallelComponents:hPtAveGen"+toTString(ptBin);
    TH1 *h1 = new TH1D(name,binLabelGen+";p^{gen,ave}_{T} (GeV);Events",50,0.9*ptBinEdges_[ptBin],1.1*ptBinEdges_[ptBin+1]);
    h1->SetMarkerStyle(20);
    hPtGenAve[ptBin] = h1;

    for(int j = 0; j < 3; ++j) {
      name = "ParallelComponents:hPtGenJet"+toTString(j)+"_"+toTString(ptBin);
      h1 = new TH1D(name,binLabelGen+";p^{gen}_{T,"+toString(1+j)+"} (GeV);Events",
		    50,j<2 ? 0.8*ptBinEdges_[ptBin] : 0.,1.2*ptBinEdges_[ptBin+1]);
      h1->SetMarkerStyle(20);
      hPtGen[ptBin].push_back(h1);

      name = "ParallelComponents:hPpOverPtGenVsDeltaPhi12";
      name += ptBin;
      name += "_Jet";
      name += j;
      TH2 *h2 = new TH2D(name,binLabelGen+";#Delta#Phi_{12};p^{gen}_{||,"+toString(1+j)+"} / p^{gen}_{T,"+toString(1+j)+"};Events",
			100,0.,3.5,100,0.6,1.05);
      hPpOverPtGenVsDeltaPhi12[ptBin].push_back(h2);
    }

    name = "ParallelComponents:hDeltaPhi"+toTString(ptBin);
    h1 = new TH1D(name,binLabelGen+";#Delta#Phi_{12};Events",50,0.,M_PI);
    h1->SetMarkerStyle(20);
    hDeltaPhi[ptBin] = h1;

    name = "ParallelComponents:hPPhi"+toTString(ptBin);
    h1 = new TH1D(name,binLabelGen+";#Phi_{||};Events",50,-6.,6.);
    h1->SetMarkerStyle(20);
    hPPhi[ptBin] = h1;

    name = "ParallelComponents:hPJet3"+toTString(ptBin);
    h1 = new TH1D(name,binLabelGen+";p^{gen}_{||,3} (GeV);Events",50,0.,0.8*ptBinEdges_[ptBin]);
    h1->SetMarkerStyle(20);
    hPJet3[ptBin] = h1;

    name = "ParallelComponents:hPSJ"+toTString(ptBin);
    h1 = new TH1D(name,binLabelGen+";p^{gen}_{||,SJ} (GeV);Events",50,0.,0.8*ptBinEdges_[ptBin]);
    h1->SetMarkerStyle(20);
    hPSJ[ptBin] = h1;

    name = "ParallelComponents:hPUCE"+toTString(ptBin);
    h1 = new TH1D(name,binLabelGen+";p^{gen}_{||,UCE} (GeV);Events",50,0.,0.8*ptBinEdges_[ptBin]);
    h1->SetMarkerStyle(20);
    hPUCE[ptBin] = h1;

    name = "ParallelComponents:hPJet3Rel"+toTString(ptBin);
    h1 = new TH1D(name,binLabelGen+";p^{gen}_{||,3} / <p^{gen,ave}_{T}>;Events",50,0.,0.8);
    h1->SetMarkerStyle(20);
    hPJet3Rel[ptBin] = h1;

    name = "ParallelComponents:hPSJRel"+toTString(ptBin);
    h1 = new TH1D(name,binLabelGen+";p^{gen}_{||,SJ} / <p^{gen,ave}_{T}>;Events",50,0.,0.8);
    h1->SetMarkerStyle(20);
    hPSJRel[ptBin] = h1;

    name = "ParallelComponents:hPUCERel"+toTString(ptBin);
    h1 = new TH1D(name,binLabelGen+";p^{gen}_{||,UCE} / <p^{gen,ave}_{T}>;Events",50,0.,0.8);
    h1->SetMarkerStyle(20);
    hPUCERel[ptBin] = h1;

    name = "ParallelComponents:hResp"+toTString(ptBin);
    h1 = new TH1D(name,binLabelGen+";Response;Events",51,0.,2.);
    h1->SetMarkerStyle(20);
    hResp[ptBin] = h1;

    name = "ParallelComponents:hPtAsym"+toTString(ptBin);
    h1 = new TH1D(name,binLabelGen+";Asymmetry;Events",51,-1.,1.);
    h1->SetMarkerStyle(20);
    hPtAsym[ptBin] = h1;

    name = "ParallelComponents:hPtGenAsym"+toTString(ptBin);
    h1 = new TH1D(name,binLabelGen+";p^{gen}_{T} Asymmetry;Events",51,-1.,1.);
    h1->SetMarkerStyle(20);
    hPtGenAsym[ptBin] = h1;

    for(size_t c = 0; c < pt3Limits.size(); ++c) {
      name = "ParallelComponents:hPtAsymPt3Cut"+toTString(pt3Limits[c])+"_"+toTString(ptBin);
      h1 = new TH1D(name,binLabelGen+";Asymmetry",51,-1.,1.);
      h1->SetYTitle("1 / N  Events / "+toTString(h1->GetBinWidth(1)));
      h1->Sumw2();
      h1->SetMarkerStyle(20);
      setColor(h1,color(c));
      hPtAsymPt3Cuts[ptBin].push_back(h1);

      name = "ParallelComponents:hPtAsymPp3Cut"+toTString(pt3Limits[c])+"_"+toTString(ptBin);
      h1 = new TH1D(name,binLabelGen+";Asymmetry",51,-1.,1.);
      h1->SetYTitle("1 / N  Events / "+toTString(h1->GetBinWidth(1)));
      h1->Sumw2();
      h1->SetMarkerStyle(20);
      setColor(h1,color(c));
      hPtAsymPp3Cuts[ptBin].push_back(h1);
    }
  }



  // Fill histograms
  std::cout << "  Filling histograms" << std::endl;

  // 1. loop over data
  for(DataIt datait = data_->begin(); datait != data_->end(); datait++) {
    // Select DiJet events
    if( (*datait)->type() == TypeSmearDiJet )  {
      SmearDiJet * dijet = static_cast<SmearDiJet*>(*datait);  
      double ptGenAve = 0.5*(dijet->jet1()->genPt()+dijet->jet2()->genPt());
      int ptGenAveBin = findPtBin(ptGenAve);
      if( ptGenAveBin >= 0 ) {
	hPtGenAve[ptGenAveBin]->Fill(ptGenAve,dijet->weight());
	for(int i = 0; i < (dijet->jet3() ? 3 : 2); ++i) {
	  hPtGen[ptGenAveBin][i]->Fill(dijet->jet(i)->genPt(),dijet->weight());
	  hPpOverPtGenVsDeltaPhi12[ptGenAveBin][i]->Fill(dijet->deltaPhi12(),cos(TVector2::Phi_mpi_pi(dijet->jet(i)->phi()-dijet->pPhi())),dijet->weight());
	}
	hDeltaPhi[ptGenAveBin]->Fill(dijet->deltaPhi12(),dijet->weight());
	hPPhi[ptGenAveBin]->Fill(dijet->pPhi(),dijet->weight());
	
	hResp[ptGenAveBin]->Fill(dijet->jet1()->pt()/dijet->jet1()->genPt(),dijet->weight());
	hResp[ptGenAveBin]->Fill(dijet->jet2()->pt()/dijet->jet2()->genPt(),dijet->weight());

	// Pt asymmetry
	double ptAsym = dijet->jet1()->pt() + dijet->jet2()->pt();
	if( ptAsym > 0. ) {
	  ptAsym = (dijet->jet1()->pt() - dijet->jet2()->pt()) / ptAsym;
	  if( rand_->Uniform(0.,1.) > 0.5 ) ptAsym *= -1.;
	} else {
	  std::cerr << "WARNING: pt1 + pt2 <= 0" << std::endl;
	  ptAsym = -2.;
	}
	hPtAsym[ptGenAveBin]->Fill(ptAsym,dijet->weight());
	double pt3 = dijet->jet3() ? dijet->jet3()->genPt() : 0.;
	for(size_t c = 0; c < pt3Limits.size(); ++c) {
	  double maxPt3 = pt3Limits[c]*ptGenAve;
	  if( pt3 < maxPt3 ) hPtAsymPt3Cuts[ptGenAveBin][c]->Fill(ptAsym,dijet->weight());
	  if( dijet->pJ3() < maxPt3 ) hPtAsymPp3Cuts[ptGenAveBin][c]->Fill(ptAsym,dijet->weight());
	}

	// PtGen asymmetry
	double ptGenAsym = dijet->jet1()->genPt() + dijet->jet2()->genPt();
	if( ptGenAsym > 0. ) {
	  ptGenAsym = (dijet->jet1()->genPt() - dijet->jet2()->genPt()) / ptGenAsym;
	  if( rand_->Uniform(0.,1.) > 0.5 ) ptGenAsym *= -1.;
	} else {
	  std::cerr << "WARNING: pt1Gen + pt2Gen <= 0" << std::endl;
	  ptGenAsym = -2.;
	}
	hPtGenAsym[ptGenAveBin]->Fill(ptGenAsym,dijet->weight());

      }
    }
  }
  for(int ptBin = 0; ptBin < nPtBins(); ++ptBin) {
    for(size_t c = 0; c < pt3Limits.size(); ++c) {
      normHist(hPtAsymPt3Cuts[ptBin][c],"width");
      normHist(hPtAsymPp3Cuts[ptBin][c],"width");
    }
  }  

  std::vector<double> meanPtGenAve(nPtBins(),1.);
  for(size_t ptBin = 0; ptBin < meanPtGenAve.size(); ++ptBin) {
    meanPtGenAve[ptBin] = hPtGenAve[ptBin]->GetMean();
  }
  // 2. loop over data
  for(DataIt datait = data_->begin(); datait != data_->end(); datait++) {
    // Select DiJet events
    if( (*datait)->type() == TypeSmearDiJet )  {
      SmearDiJet * dijet = static_cast<SmearDiJet*>(*datait);  
      double ptGenAve = 0.5*(dijet->jet1()->genPt()+dijet->jet2()->genPt());
      int ptGenAveBin = findPtBin(ptGenAve);
      if( ptGenAveBin >= 0 ) {
	hPJet3[ptGenAveBin]->Fill(dijet->pJ3(),dijet->weight());
	hPSJ[ptGenAveBin]->Fill(dijet->pSJ(),dijet->weight());
	hPUCE[ptGenAveBin]->Fill(dijet->pUCE(),dijet->weight());

	hPJet3Rel[ptGenAveBin]->Fill(std::abs(dijet->pJ3())/meanPtGenAve[ptGenAveBin],dijet->weight());
	hPSJRel[ptGenAveBin]->Fill(std::abs(dijet->pSJ())/meanPtGenAve[ptGenAveBin],dijet->weight());
	hPUCERel[ptGenAveBin]->Fill(std::abs(dijet->pUCE())/meanPtGenAve[ptGenAveBin],dijet->weight());
      }
    }
  }



  // Contributions to asymmetry
  std::vector<double> meanPtGenAveErr(meanPtGenAve.size());
  std::vector<double> stdDevPtAsym(meanPtGenAve.size());
  std::vector<double> stdDevPtAsymErr(meanPtGenAve.size());
  std::vector<double> stdDevPtGenAsym(meanPtGenAve.size());
  std::vector<double> stdDevPtGenAsymErr(meanPtGenAve.size());
  std::vector<double> stdDevResp(meanPtGenAve.size());
  std::vector<double> stdDevRespErr(meanPtGenAve.size());
  std::vector<double> stdDevPJet3(meanPtGenAve.size());
  std::vector<double> stdDevPJet3Err(meanPtGenAve.size());
  std::vector<double> stdDevPSJ(meanPtGenAve.size());
  std::vector<double> stdDevPSJErr(meanPtGenAve.size());
  std::vector<double> stdDevPUCE(meanPtGenAve.size());
  std::vector<double> stdDevPUCEErr(meanPtGenAve.size());
  for(size_t i = 0; i < meanPtGenAve.size(); ++i) {
    meanPtGenAveErr[i] = hPtGenAve[i]->GetMeanError();

    // Standard deviation
    stdDevPtAsym[i] = sqrt(2.)*hPtAsym[i]->GetRMS();
    stdDevPtAsymErr[i] = sqrt(2.)*hPtAsym[i]->GetRMSError();

    stdDevPtGenAsym[i] = sqrt(2.)*hPtGenAsym[i]->GetRMS();
    stdDevPtGenAsymErr[i] = sqrt(2.)*hPtGenAsym[i]->GetRMSError();
    
    stdDevResp[i] = hResp[i]->GetRMS();
    stdDevRespErr[i] = hResp[i]->GetRMSError();
    
    stdDevPJet3[i] = hPJet3Rel[i]->GetRMS();
    stdDevPJet3Err[i] = hPJet3Rel[i]->GetRMSError();
    stdDevPSJ[i] = hPSJRel[i]->GetRMS();
    stdDevPSJErr[i] = hPSJRel[i]->GetRMSError();
    stdDevPUCE[i] = hPUCERel[i]->GetRMS();
    stdDevPUCEErr[i] = hPUCERel[i]->GetRMSError();
  }

  TGraphErrors *gStdDevPtAsym
    = new TGraphErrors(meanPtGenAve.size(),&(meanPtGenAve.front()),&(stdDevPtAsym.front()),
		       &(meanPtGenAveErr.front()),&(stdDevPtAsymErr.front()));
  gStdDevPtAsym->SetMarkerStyle(20);
  gStdDevPtAsym->SetMarkerColor(color(0));
  gStdDevPtAsym->SetLineColor(color(0));

  TGraphErrors *gStdDevPtGenAsym
    = new TGraphErrors(meanPtGenAve.size(),&(meanPtGenAve.front()),&(stdDevPtGenAsym.front()),
		       &(meanPtGenAveErr.front()),&(stdDevPtGenAsymErr.front()));
  gStdDevPtGenAsym->SetMarkerStyle(24);
  gStdDevPtGenAsym->SetMarkerColor(color(4));
  gStdDevPtGenAsym->SetLineColor(color(4));

  TGraphErrors *gStdDevResp
    = new TGraphErrors(meanPtGenAve.size(),&(meanPtGenAve.front()),&(stdDevResp.front()),
		       &(meanPtGenAveErr.front()),&(stdDevRespErr.front()));
  gStdDevResp->SetMarkerStyle(21);
  gStdDevResp->SetMarkerColor(color(1));
  gStdDevResp->SetLineColor(color(1));

  TGraphErrors *gStdDevPJet3
    = new TGraphErrors(meanPtGenAve.size(),&(meanPtGenAve.front()),&(stdDevPJet3.front()),
		       &(meanPtGenAveErr.front()),&(stdDevPJet3Err.front()));
  gStdDevPJet3->SetMarkerStyle(22);
  gStdDevPJet3->SetMarkerColor(color(2));
  gStdDevPJet3->SetLineColor(color(2));

  TGraphErrors *gStdDevPSJ
    = new TGraphErrors(meanPtGenAve.size(),&(meanPtGenAve.front()),&(stdDevPSJ.front()),
		       &(meanPtGenAveErr.front()),&(stdDevPSJErr.front()));
  gStdDevPSJ->SetMarkerStyle(23);
  gStdDevPSJ->SetMarkerColor(color(3));
  gStdDevPSJ->SetLineColor(color(3));

  TGraphErrors *gStdDevPUCE
    = new TGraphErrors(meanPtGenAve.size(),&(meanPtGenAve.front()),&(stdDevPUCE.front()),
		       &(meanPtGenAveErr.front()),&(stdDevPUCEErr.front()));
  gStdDevPUCE->SetMarkerStyle(24);
  gStdDevPUCE->SetMarkerColor(color(4));
  gStdDevPUCE->SetLineColor(color(4));



  // Dependence of asymmetry on pt3
  std::vector<TGraphErrors*> gResPtAsymPt3Cuts(pt3Limits.size());
  std::vector<TGraphErrors*> gResPtAsymPp3Cuts(pt3Limits.size());
  for(size_t c = 0; c < pt3Limits.size(); ++c) {
    std::vector<double> resPtAsymPt3Cuts(meanPtGenAve.size());
    std::vector<double> resPtAsymPt3CutsErr(meanPtGenAve.size());
    std::vector<double> resPtAsymPp3Cuts(meanPtGenAve.size());
    std::vector<double> resPtAsymPp3CutsErr(meanPtGenAve.size());
    for(size_t i = 0; i < meanPtGenAve.size(); ++i) {
      resPtAsymPt3Cuts[i] = sqrt(2.)*hPtAsymPt3Cuts[i][c]->GetRMS();
      resPtAsymPt3CutsErr[i] = sqrt(2.)*hPtAsymPt3Cuts[i][c]->GetRMSError();
      resPtAsymPp3Cuts[i] = sqrt(2.)*hPtAsymPp3Cuts[i][c]->GetRMS();
      resPtAsymPp3CutsErr[i] = sqrt(2.)*hPtAsymPp3Cuts[i][c]->GetRMSError();

    }
    gResPtAsymPt3Cuts[c] = new TGraphErrors(meanPtGenAve.size(),&(meanPtGenAve.front()),&(resPtAsymPt3Cuts.front()),&(meanPtGenAveErr.front()),&(resPtAsymPt3CutsErr.front()));
    gResPtAsymPt3Cuts[c]->SetMarkerStyle(24+c);
    gResPtAsymPt3Cuts[c]->SetLineColor(color(c));
    gResPtAsymPt3Cuts[c]->SetMarkerColor(color(c));

    gResPtAsymPp3Cuts[c] = new TGraphErrors(meanPtGenAve.size(),&(meanPtGenAve.front()),&(resPtAsymPp3Cuts.front()),&(meanPtGenAveErr.front()),&(resPtAsymPp3CutsErr.front()));
    gResPtAsymPp3Cuts[c]->SetMarkerStyle(20+c);
    gResPtAsymPp3Cuts[c]->SetLineColor(color(c));
    gResPtAsymPp3Cuts[c]->SetMarkerColor(color(c));
  }



  // Plot
  std::cout << "  Plotting histograms" << std::endl;
  TString outName = outNamePrefix_+"ParallelComponent_";
  TPostScript *ps = 0;
  if( !saveAsEps_ ) ps = new TPostScript((dir_+"/jsParallelComponent.ps").c_str(),111);
  TCanvas *c1 = new TCanvas("c1","Parallel Component",0,0,600,600);

  for(int ptBin = 0; ptBin < nPtBins(); ptBin++) {
    if( !saveAsEps_ ) ps->NewPage();
    c1->cd();
    hPtGenAve[ptBin]->Draw("PE1");
    if( !saveAsEps_ ) c1->Draw();
    else c1->SaveAs(outName+"hPtGenAve_PtBin"+toString(ptBin)+".eps","eps");

    if( !saveAsEps_ ) ps->NewPage();
    c1->cd();
    hDeltaPhi[ptBin]->Draw("PE1");
    if( !saveAsEps_ ) c1->Draw();
    else c1->SaveAs(outName+"hDeltaPhi_PtBin"+toString(ptBin)+".eps","eps");

    if( !saveAsEps_ ) ps->NewPage();
    c1->cd();
    hPPhi[ptBin]->Draw("PE1");
    if( !saveAsEps_ ) c1->Draw();
    else c1->SaveAs(outName+"hPPhi_PtBin"+toString(ptBin)+".eps","eps");

    for(size_t j = 0; j < hPpOverPtGenVsDeltaPhi12[ptBin].size(); ++j) {
      if( !saveAsEps_ ) ps->NewPage();
      c1->cd();
      hPtGen[ptBin][j]->Draw("PE1");
      if( !saveAsEps_ ) c1->Draw();
      else c1->SaveAs(outName+"hPtGenJet"+toTString(1+j)+"_PtBin"+toString(ptBin)+".eps","eps");

      if( !saveAsEps_ ) ps->NewPage();
      c1->cd();
      hPpOverPtGenVsDeltaPhi12[ptBin][j]->Draw("COLZ");
      if( !saveAsEps_ ) c1->Draw();
      else c1->SaveAs(outName+"hPpOverPtGenJet"+toTString(1+j)+"VsDeltaPhi12_PtBin"+toTString(ptBin)+".eps","eps");
    }

    if( !saveAsEps_ ) ps->NewPage();
    c1->cd();
    hResp[ptBin]->GetXaxis()->SetRangeUser(0.3,1.7);
    hResp[ptBin]->Draw("PE1");
    if( !saveAsEps_ ) c1->Draw();
    else c1->SaveAs(outName+"hResp_PtBin"+toString(ptBin)+".eps","eps");

    if( !saveAsEps_ ) ps->NewPage();
    c1->cd();
    hPtAsym[ptBin]->GetXaxis()->SetRangeUser(-0.7,0.7);
    hPtAsym[ptBin]->Draw("PE1");
    if( !saveAsEps_ ) c1->Draw();
    else c1->SaveAs(outName+"hPtAsym_PtBin"+toString(ptBin)+".eps","eps");

    if( !saveAsEps_ ) ps->NewPage();
    c1->cd();
    hPtGenAsym[ptBin]->GetXaxis()->SetRangeUser(-0.7,0.7);
    hPtGenAsym[ptBin]->Draw("PE1");
    if( !saveAsEps_ ) c1->Draw();
    else c1->SaveAs(outName+"hPtGenAsym_PtBin"+toString(ptBin)+".eps","eps");

    if( !saveAsEps_ ) ps->NewPage();
    c1->cd();
    hPtAsymPt3Cuts[ptBin][0]->GetXaxis()->SetRangeUser(-0.7,0.7);
    hPtAsymPt3Cuts[ptBin][0]->Draw("PE1");
    for(size_t c = 1; c < pt3Limits.size(); ++c) {
      hPtAsymPt3Cuts[ptBin][c]->Draw("PE1same");
    }
    if( !saveAsEps_ ) c1->Draw();
    else c1->SaveAs(outName+"hPtAsymPt3Cuts_PtBin"+toString(ptBin)+".eps","eps");

    if( !saveAsEps_ ) ps->NewPage();
    c1->cd();
    hPtAsymPp3Cuts[ptBin][0]->GetXaxis()->SetRangeUser(-0.7,0.7);
    hPtAsymPp3Cuts[ptBin][0]->Draw("PE1");
    for(size_t c = 1; c < pt3Limits.size(); ++c) {
      hPtAsymPp3Cuts[ptBin][c]->Draw("PE1same");
    }
    if( !saveAsEps_ ) c1->Draw();
    else c1->SaveAs(outName+"hPpAsymPt3Cuts_PtBin"+toString(ptBin)+".eps","eps");

    if( !saveAsEps_ ) ps->NewPage();
    c1->cd();
    hPJet3[ptBin]->Draw("PE1");
    if( !saveAsEps_ ) c1->Draw();
    else c1->SaveAs(outName+"hPJet3_PtBin"+toString(ptBin)+".eps","eps");

    if( !saveAsEps_ ) ps->NewPage();
    c1->cd();
    hPSJ[ptBin]->Draw("PE1");
    if( !saveAsEps_ ) c1->Draw();
    else c1->SaveAs(outName+"hPSJ_PtBin"+toString(ptBin)+".eps","eps");

    if( !saveAsEps_ ) ps->NewPage();
    c1->cd();
    hPUCE[ptBin]->Draw("PE1");
    if( !saveAsEps_ ) c1->Draw();
    else c1->SaveAs(outName+"hPUCE_PtBin"+toString(ptBin)+".eps","eps");

    if( !saveAsEps_ ) ps->NewPage();
    c1->cd();
    hPJet3Rel[ptBin]->Draw("PE1");
    if( !saveAsEps_ ) c1->Draw();
    else c1->SaveAs(outName+"hPJet3Rel_PtBin"+toString(ptBin)+".eps","eps");

    if( !saveAsEps_ ) ps->NewPage();
    c1->cd();
    hPSJRel[ptBin]->Draw("PE1");
    if( !saveAsEps_ ) c1->Draw();
    else c1->SaveAs(outName+"hPSJRel_PtBin"+toString(ptBin)+".eps","eps");

    if( !saveAsEps_ ) ps->NewPage();
    c1->cd();
    hPUCERel[ptBin]->Draw("PE1");
    if( !saveAsEps_ ) c1->Draw();
    else c1->SaveAs(outName+"hPUCERel_PtBin"+toString(ptBin)+".eps","eps");
  }

  if( !saveAsEps_ ) ps->NewPage();
  c1->cd();

  TF1 *fPtGenAsym = new TF1("fPtGenAsym","sqrt([0]*[0]/x[0]/x[0]+[1]*[1]/x[0]+[2]*[2])",
			      ptBinEdges_[0],ptBinEdges_[nPtBins()]);
  fPtGenAsym->SetParameter(0,2.54877);
  fPtGenAsym->SetParameter(1,0.149045);
  fPtGenAsym->SetParameter(0,0.0109168);
  fPtGenAsym->SetLineWidth(1);
  fPtGenAsym->SetLineStyle(2);
  fPtGenAsym->SetLineColor(color(4));


  TH1 *hSumContributions = new TH1D("hSumContributions",";p^{gen}_{T} (GeV);Standard Deviation",
				    nPtBins(),&(ptBinEdges_.front()));
  for(int i = 0; i < nPtBins(); ++i) {
    double resp = gStdDevResp->GetY()[i];
    double pJ3 = gStdDevPJet3->GetY()[i];
    double pSJ = gStdDevPSJ->GetY()[i];
    double ptGenAsym = fPtGenAsym->Eval(meanPtGenAve[i]);
    //double pUCE = gStdDevPUCE->GetY()[i];
    hSumContributions->SetBinContent(1+i,sqrt( resp*resp + pJ3*pJ3 + pSJ*pSJ + ptGenAsym*ptGenAsym ));
  }
  hSumContributions->GetYaxis()->SetRangeUser(0,0.8);
  hSumContributions->GetXaxis()->SetMoreLogLabels();
  hSumContributions->Draw();

  TLegend *leg = createLegend(5,0.78);
  TF1 *fMCTruthResp = new TF1("fMCTruthResp","sqrt([0]*[0]/x[0]/x[0]+[1]*[1]/x[0]+[2]*[2])",
			      ptBinEdges_[0],ptBinEdges_[nPtBins()]);
  for(int i = 0; i < 3; ++i) {
    fMCTruthResp->SetParameter(i,truthPar_[i]);
  }
  fMCTruthResp->SetLineWidth(1);
  fMCTruthResp->SetLineStyle(2);
  fMCTruthResp->SetLineColor(color(4));
  //  fMCTruthResp->Draw("same");

  fPtGenAsym->Draw("same");
  gStdDevResp->Draw("PE1same");

  gStdDevPtAsym->Draw("PE1same");
  gStdDevPJet3->Draw("PE1same");
  gStdDevPSJ->Draw("PE1same");
  //  gStdDevPtGenAsym->Draw("PE1same");

  leg->AddEntry(gStdDevPtAsym,"#sqrt{2} #sigma(Asymmetry)","P");
  leg->AddEntry(hSumContributions,"Total","L");
  leg->AddEntry(gStdDevResp,"MC Truth #sigma(R_{MC})","P");
  leg->AddEntry(gStdDevPJet3,"#sigma(p_{||,3})","P");
  leg->AddEntry(gStdDevPSJ,"#sigma(p_{||,Soft})","P");
  leg->AddEntry(fPtGenAsym,"#sqrt{2} #sigma(gen-Asymmetry)","L");

  leg->Draw("same");
  c1->SetLogx();
  if( !saveAsEps_ ) c1->Draw();
  else c1->SaveAs(outName+"hParallelContributions.eps","eps");


  if( !saveAsEps_ ) ps->NewPage();
  c1->cd();
  TH1 *hFramePtAsym = new TH1D("hFramePtAsym",";<p^{gen,ave}> (GeV);Resolution",
			       1000,ptBinEdges_[0],ptBinEdges_[nPtBins()]);
  hFramePtAsym->GetYaxis()->SetRangeUser(0,0.6);
  hFramePtAsym->Draw();
  fMCTruthResp->Draw("same");
  TLegend *legPtAsym = createLegend(5,0.5);
  for(size_t c = 0; c < pt3Limits.size(); ++c) {
    gResPtAsymPt3Cuts[c]->Draw("PE1same");
    legPtAsym->AddEntry(gResPtAsymPt3Cuts[c],"p_{T,3} < "+toTString(pt3Limits[c])+" p^{gen,ave}","P");
    gResPtAsymPp3Cuts[c]->Draw("PE1same");
    legPtAsym->AddEntry(gResPtAsymPp3Cuts[c],"p_{||,3} < "+toTString(pt3Limits[c])+" p^{gen,ave}","P");    
  }
  legPtAsym->Draw("same");
  c1->SetLogx();
  if( !saveAsEps_ ) c1->Draw();
  else c1->SaveAs(outName+"hPtAsymmetryResolutions.eps","eps");



  if( !saveAsEps_ ) ps->Close();


  // Delete
  std::cout << "  Cleaning up" << std::endl;
  for(size_t i = 0; i < hPPhi.size(); ++i) {
    delete hPtGenAve[i];
    delete hPPhi[i];
    delete hPtAsym[i];
    delete hPJet3[i];
    delete hPSJ[i];
    delete hPUCE[i];
    delete hPJet3Rel[i];
    delete hPSJRel[i];
    delete hPUCERel[i];
    delete hResp[i];
    for(size_t j = 0; j < hPpOverPtGenVsDeltaPhi12[i].size(); ++j) {
      delete hPtGen[i][j];
      delete hPpOverPtGenVsDeltaPhi12[i][j];
    }
    for(size_t c = 0; c < pt3Limits.size(); ++c) {
      delete hPtAsymPt3Cuts[i][c];
      delete hPtAsymPp3Cuts[i][c];
    }
  }
  for(size_t c = 0; c < pt3Limits.size(); ++c) {
    delete gResPtAsymPt3Cuts[c];
    delete gResPtAsymPp3Cuts[c];
  }
  delete hSumContributions;
  delete fMCTruthResp;
  delete legPtAsym;
  delete leg;
  delete gStdDevPtAsym;
  delete gStdDevResp;
  delete gStdDevPJet3;
  delete gStdDevPSJ;
  delete gStdDevPUCE;
  delete c1;
  if( ps ) delete ps;
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
  TPostScript *ps = 0;
  if( !saveAsEps_ ) ps = new TPostScript((dir_+"/jsSimulatedAsymmetry.ps").c_str(),111);
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

    hAbsPars->SetBinContent(bin,scale(i)*param_->GetPars()[i]);
    hAbsPars->SetBinError(bin,scale(i)*param_->GetErrors()[i]);
  
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
  std::vector<double> bins(35+1);
  equidistLogBins(bins,35,ptBinsMin(),ptBinsMax());
  TH2* hRespVsPtGen = new TH2D("MeanResp_hRespVsPtGen",
			       ";p^{gen}_{T} (GeV);p^{jet}_{T} / p^{gen}_{T}",
			       bins.size()-1,&(bins.front()),51,0.,2.);
  hRespVsPtGen->SetNdivisions(505);
  hRespVsPtGen->Sumw2();
  

  // ----- Fill histograms -----
  for(DataIt datait = data_->begin(); datait != data_->end(); datait++) {
    // Select DiJet events
    if( (*datait)->type() == TypeSmearDiJet )  {
      SmearDiJet * dijet = static_cast<SmearDiJet*>(*datait); 
      hRespVsPtGen->Fill(dijet->jet1()->genPt(),dijet->jet1()->pt()/dijet->jet1()->genPt(),dijet->weight());
      hRespVsPtGen->Fill(dijet->jet2()->genPt(),dijet->jet2()->pt()/dijet->jet2()->genPt(),dijet->weight());
    }
  }


  // ----- Fit profiles -----
  std::vector<TH1*> hResp(2);
  hResp[0] = new TH1D("MeanResp_hRespMean",";p^{gen}_{T} (GeV);< p^{jet}_{T} / p^{gen}_{T} >",
		      hRespVsPtGen->GetNbinsX(),hRespVsPtGen->GetXaxis()->GetXbins()->GetArray());
  hResp[0]->SetNdivisions(505);
  hResp[0]->SetMarkerStyle(20);
  hResp[0]->SetMarkerColor(4);
  hResp[0]->SetLineColor(hResp[0]->GetMarkerColor());
  hResp[1] = static_cast<TH1D*>(hResp[0]->Clone("MeanResp_hRespGauss"));
  hResp[1]->SetMarkerColor(2);
  hResp[1]->SetLineColor(hResp[1]->GetMarkerColor());


  std::vector<TH1*> hReso(2);
  hReso[0] = static_cast<TH1D*>(hResp[0]->Clone("MeanResp_hResoStd"));
  hReso[0]->SetMarkerColor(1);
  hReso[0]->SetYTitle("Standard Deviation");
  hReso[0]->SetLineColor(hReso[0]->GetMarkerColor());
  hReso[1] = static_cast<TH1D*>(hReso[0]->Clone("MeanResp_hResoGauss"));
  hReso[1]->SetYTitle("Gauss Fit #sigma / #mu");

  // Get 1D slices and get mean, sigma
  TH1* hSlice = new TH1D("hSlice","",hRespVsPtGen->GetNbinsY(),hRespVsPtGen->GetYaxis()->GetXmin(),
			 hRespVsPtGen->GetYaxis()->GetXmax());
  hSlice->Sumw2();
  for(int xBin = 1; xBin <= hRespVsPtGen->GetNbinsX(); xBin++) {
    hSlice->Reset();
    for(int yBin = 1; yBin <= hRespVsPtGen->GetNbinsY(); yBin++) {
      hSlice->SetBinContent(yBin,hRespVsPtGen->GetBinContent(hRespVsPtGen->GetBin(xBin,yBin)));
      hSlice->SetBinError(yBin,hRespVsPtGen->GetBinError(hRespVsPtGen->GetBin(xBin,yBin)));
    }  

    double mean = hSlice->GetMean();
    double width = hSlice->GetRMS();

    if( width < 0.1 ) width = 0.1;
    if( hSlice->GetSumOfWeights() <= 0 ) continue;
    hSlice->Fit("gaus","QNO","",mean-3.*width,mean+3.*width);
    TF1 *f = static_cast<TF1*>(gROOT->GetFunction("gaus")->Clone());
    mean = f->GetParameter(1);
    double meanErr = f->GetParError(1);
    width = f->GetParameter(2);
    if( width < 0.05 ) width = 0.05;
    if( (hSlice->Fit(f,"LLQNO","goff",mean-1.5*width,mean+1.5*width) == 0) ) {
      mean = f->GetParameter(1);
      meanErr = f->GetParError(1);
      width = f->GetParameter(2);
      hResp[1]->SetBinContent(xBin,mean);
      hResp[1]->SetBinError(xBin,meanErr);
      hReso[1]->SetBinContent(xBin,width);///mean);
      hReso[1]->SetBinError(xBin,f->GetParError(2));///mean);
    }
    delete f;

//     if( width < 0.1 ) width = 0.1;
//     if( hSlice->GetSumOfWeights() <= 0 ) continue;
//     TF1 *f = new TF1("fit",gaussian,mean-3.*width,mean+3.*width,2);
//     f->SetParameter(0,hSlice->Integral());
//     f->SetParameter(1,0.1);
//     for(int bin = 1; bin <= hSlice->GetNbinsX(); ++bin) {
//       std::cout << f->Eval(hSlice->GetBinCenter(bin)) << std::endl;
//     }
//     hSlice->Fit(f,"0R");
//     width = f->GetParameter(1);
//     if( width < 0.05 ) width = 0.05;
//     if( (hSlice->Fit(f,"LLQNO","goff",mean-2.5*width,mean+2.5*width) == 0) ) {
//       width = f->GetParameter(1);
//       hResp[1]->SetBinContent(xBin,mean);
//       hResp[1]->SetBinError(xBin,0.);
//       hReso[1]->SetBinContent(xBin,width);
//       hReso[1]->SetBinError(xBin,f->GetParError(1));
//     }
//     delete f;

    mean = hSlice->GetMean();
    meanErr = hSlice->GetMeanError();
    width = hSlice->GetRMS();
    hResp[0]->SetBinContent(xBin,mean);
    hResp[0]->SetBinError(xBin,meanErr);
    hReso[0]->SetBinContent(xBin,width);///mean);
    hReso[0]->SetBinError(xBin,hSlice->GetRMSError());///mean);
  }
  delete hSlice;


  // ----- Fit resolution -----
  std::vector<TF1*> fReso(2);
  fReso[0] = new TF1("MeanResp_fitResoStd","sqrt([0]*[0]/x/x + [1]*[1]/x + [2]*[2])",
		     ptBinEdges_.front(),ptBinEdges_.back());
  fReso[1] = static_cast<TF1*>(fReso[0]->Clone("MeanResp_fitResoGauss"));
  std::vector<TPaveText*> fitLabel(2);

  for(size_t i = 0; i < 2; ++i) {
    fReso[i]->SetLineWidth(1);
    fReso[i]->SetLineColor(2);
    fReso[i]->SetParameter(0,2.);
    fReso[i]->SetParameter(1,1.3);
    fReso[i]->SetParameter(2,0.04);
    hReso[i]->Fit(fReso[i],"INQR");

    fitLabel[i] = createPaveText(2);
    char txt[50];
    sprintf(txt,"a_{0} = %.2f #pm %.2f #sqrt{GeV}",fReso[i]->GetParameter(0),fReso[i]->GetParError(0));
    fitLabel[i]->AddText(txt);
    sprintf(txt,"a_{1} = %.2f #pm %.2f #sqrt{GeV}",fReso[i]->GetParameter(1),fReso[i]->GetParError(1));
    fitLabel[i]->AddText(txt);
    sprintf(txt,"a_{2} = %.3f #pm %.3f GeV",fReso[i]->GetParameter(2),fReso[i]->GetParError(2));
    fitLabel[i]->AddText(txt);    
  }


  // ----- Plot histograms -----
  TPostScript * const ps = new TPostScript((dir_+"/jsMeanRespAndReso.ps").c_str(),111);
  TCanvas * c1 = new TCanvas("c1","MeanRespAndReso",0,0,600,600);

  ps->NewPage();
  c1->cd();
  hResp[0]->GetYaxis()->SetRangeUser(0.9,1.2);
  hResp[0]->Draw("PE1");
  hResp[1]->Draw("PE1same");
  TLegend *leg  = createLegend(2);
  leg->AddEntry(hResp[0],"Mean","P");
  leg->AddEntry(hResp[1],"Gauss Fit Mean","P");
  leg->Draw("same");
  c1->SetLogx();
  c1->Draw();

  for(size_t i = 0; i < 2; ++i) {
    ps->NewPage();
    c1->cd();
    hReso[i]->GetYaxis()->SetRangeUser(0.,0.2);
    hReso[i]->Draw("PE1");
    fReso[i]->Draw("same");
    fitLabel[i]->Draw("same");
    c1->SetLogx();
    c1->Draw();
  }    
  
  TFile rootfile((dir_+"/jsResponse.root").c_str(),"UPDATE");
  rootfile.WriteTObject(hResp[0]);
  rootfile.WriteTObject(hResp[1]);
  rootfile.WriteTObject(hReso[0]);
  rootfile.WriteTObject(hReso[1]);


  // ----- Clean up -----
  for(size_t i = 0; i < 2; ++i) {
    delete hResp[i];
    delete hReso[i];
    delete fReso[i];
    delete fitLabel[i];
  }
  delete leg;
  delete c1;
  delete ps;
}



// // --------------------------------------------------
// void ControlPlotsJetSmearing::plotFlavors() const {
//   std::cout << "Creating flavor control plots\n";


//   // ----- Create histograms -----

//   std::vector<double> bins(35);
//   equidistLogBins(bins,bins.size()-1,ptBinsMin(),ptBinsMax());

//   // 0: Unknown
//   // 1: g
//   // 2: uds
//   // 3: c
//   // 4: b
//   std::vector<TH2*> hRespVsPtGen(5);
//   for(size_t f = 0; f < hRespVsPtGen.size(); ++f) {
//     hRespVsPtGen[f] = new TH2D("Flavor"+toTString(f)+"_hRespVsPtGen",
// 			       ";p^{gen}_{T} (GeV);p^{calo}_{T} / p^{gen}_{T}",
// 			       bins.size()-1,&(bins.front()),51,0.,2.);
//     hRespVsPtGen[f]->SetNdivisions(505);
//     hRespVsPtGen[f]->Sumw2();
//   }
//   // 0: Everything else
//   // 1: One b + anything else
//   // 2: Two b
//   for(size_t f = 0; f < hAsymVsPtGen.size(); ++f) {
//     std::vector<TH2*> hAsymVsPtGen(3);
//     hAsymVsPtGen[f] = new TH2D("Flavor"+toTString(f)+"_hAsymVsPtGen",
// 			       ";p^{gen,ave}_{T} (GeV);Asymmetry",
// 			       bins.size()-1,&(bins.front()),51,-1.,1.);
//     hAsymVsPtGen[f]->SetNdivisions(505);
//     hAsymVsPtGen[f]->Sumw2();
//   }

  

//   // ----- Fill histograms -----
//   for(DataIt datait = data_->begin(); datait != data_->end(); datait++) {
//     // Select DiJet events
//     if( (*datait)->type() == TypeSmearDiJet )  {
//       SmearDiJet * dijet = static_cast<SmearDiJet*>(*datait);

//       int bIdx = 0;
//       for(int i = 0; i < 2; ++i) {
// 	int fIdx = 0;
// 	if( dijet->jet(i)->flavor() == Jet::gluon ) fIdx = 1;
// 	else if( dijet->jet(i)->flavor() == Jet::uds ) fIdx = 2;
// 	else if( dijet->jet(i)->flavor() == Jet::c ) fIdx = 3;
// 	else if( dijet->jet(i)->flavor() == Jet::b ) {
// 	  fIdx = 4;
// 	  ++bIdx;
// 	}
// 	hRespVsPtGen[fIdx]->Fill(dijet->jet(i)->genPt(),dijet->jet(i)->pt()/dijet->jet(i)->genPt(),dijet->weight());
//       }
//       double ptAsym = dijet->jet1()->pt() + dijet->jet2()->pt();
//       if( ptAsym > 0 ) {
// 	ptAsym = (dijet->jet1()->pt() - dijet->jet2()->pt())/ptAsym;
// 	if( rand_->Uniform(0,1) > 0.5 ) ptAsym *= -1.;
// 	hAsymVsPtGen[bIdx]->Fill(0.5*(dijet->jet1()->genPt()+dijet->jet2()->genPt()),ptAsym,dijet->weight());
//       }
//     }
//   }



//   // ----- Store and fit profiles -----

//   std::vector<TH1*> hScale(hRespVsPtGen.size());
//   std::vector<TH1*> hReso(hRespVsPtGen.size());
//   std::vector< std::vector<TH1*> > hRespBins(hRespVsPtGen[0]->GetNbinsX());
//   for(size_t ptBin = 0; ptBin < hRespBins.size(); ++ptBin) {
//     hRespBins[ptBin] = std::vector<TH1*>(hRespVsPtGen.size());
//     for(size_t f = 0; f < hRespVsPtGen.size(); ++f) {
//       TH1 *h = new TH1D("hResp_PtBin"+toTString(ptBin)+"_Flavor"+toTString(f),"",
// 			hRespVsPtGen->GetNbinsY(),hRespVsPtGen->GetYaxis()->GetXmin(),
// 			hRespVsPtGen->GetYaxis()->GetXmax());
//       h->Sumw2();
//       h->SetLineColor(color(f));
//       for(int yBin = 1; yBin <= hRespVsPtGen->GetNbinsY(); yBin++) {
// 	h->SetBinContent(yBin,hRespVsPtGen->GetBinContent(hRespVsPtGen->GetBin(ptBin,yBin)));
// 	h->SetBinError(yBin,hRespVsPtGen->GetBinError(hRespVsPtGen->GetBin(ptBin,yBin)));
//       }  
//       hRespBins[ptBin][f] = h;

//       if( ptBin == 0 ) {
// 	hScale[f] = new TH1D("hScale_Flavor"+toTString(f),";p^{gen}_{T} (GeV);Scale",
// 			     hRespVsPtGen->GetNbinsX(),hRespVsPtGen->GetXaxis()->GetXbins()->GetArray());
// 	hScale[f]->SetNdivisions(505);
// 	hScale[f]->SetMarkerStyle(20+f);
// 	hScale[f]->SetMarkerColor(color(f));
// 	hScale[f]->SetLineColor(hScale[f]->GetMarkerColor());

// 	hReso[f] = new TH1D("hReso_Flavor"+toTString(f),";p^{gen}_{T} (GeV);Resolution",
// 			     hRespVsPtGen->GetNbinsX(),hRespVsPtGen->GetXaxis()->GetXbins()->GetArray());
// 	hReso[f]->SetNdivisions(505);
// 	hReso[f]->SetMarkerStyle(20+f);
// 	hReso[f]->SetMarkerColor(color(f));
// 	hReso[f]->SetLineColor(hReso[f]->GetMarkerColor());
//       }
//       double mean = h->GetMean();
//       double width = h->GetRMS();
//       if( width < 0.1 ) width = 0.1;
//       if( h->GetSumOfWeights() <= 0 ) continue;
//       h->Fit("gaus","QNO","",mean-3.*width,mean+3.*width);
//       TF1 *fit = static_cast<TF1*>(gROOT->GetFunction("gaus")->Clone());
//       mean = fit->GetParameter(1);
//       double meanErr = fit->GetParError(1);
//       width = fit->GetParameter(2);
//       if( width < 0.05 ) width = 0.05;
//       if( (h->Fit(fit,"LLQNO","goff",mean-1.5*width,mean+1.5*width) == 0) ) {
// 	mean = fit->GetParameter(1);
// 	meanErr = fit->GetParError(1);
// 	width = fit->GetParameter(2);
// 	hScale[f]->SetBinContent(ptBin,mean);
// 	hScale[f]->SetBinError(ptBin,meanErr);
// 	hReso[f]->SetBinContent(ptBin,width);
// 	hReso[f]->SetBinError(ptBin,fit->GetParError(2));
//       }
//     } // End of loop over flavors
//   } // End of loop over pt bins


//   std::vector< std::vector<TH1*> > hAsymBins(hAsymVsPtGen[0]->GetNbinsX());
//   for(size_t ptBin = 0; ptBin < hAsymBins.size(); ++ptBin) {
//     hAsymBins[ptBin] = std::vector<TH1*>(hAsymVsPtGen.size());
//     for(size_t f = 0; f < hAsymVsPtGen.size(); ++f) {
//       TH1 *h = new TH1D("hAsym_PtBin"+toTString(ptBin)+"_Flavor"+toTString(f),"",
// 			hAsymVsPtGen->GetNbinsY(),hAsymVsPtGen->GetYaxis()->GetXmin(),
// 			hAsymVsPtGen->GetYaxis()->GetXmax());
//       h->Sumw2();
//       h->SetLineColor(color(f));
//       for(int yBin = 1; yBin <= hAsymVsPtGen->GetNbinsY(); yBin++) {
// 	h->SetBinContent(yBin,hAsymVsPtGen->GetBinContent(hAsymVsPtGen->GetBin(ptBin,yBin)));
// 	h->SetBinError(yBin,hAsymVsPtGen->GetBinError(hAsymVsPtGen->GetBin(ptBin,yBin)));
//       }  
//       hAsymBins[ptBin][f] = h;
//     }
//   }



//   // ----- Plot histograms -----
//   TPostScript * const ps = new TPostScript((dir_+"/jsFlavors.ps").c_str(),111);
//   TCanvas * c1 = new TCanvas("c1","Flavors",0,0,600,600);

//   ps->NewPage();
//   c1->cd();
//   hResp[0]->GetYaxis()->SetRangeUser(0.9,1.2);
//   hResp[0]->Draw("PE1");
//   hResp[1]->Draw("PE1same");
//   TLegend *leg  = createLegend(2);
//   leg->AddEntry(hResp[0],"Mean","P");
//   leg->AddEntry(hResp[1],"Gauss Fit Mean","P");
//   leg->Draw("same");
//   c1->SetLogx();
//   c1->Draw();

//   for(size_t i = 0; i < 2; ++i) {
//     ps->NewPage();
//     c1->cd();
//     hReso[i]->GetYaxis()->SetRangeUser(0.,0.2);
//     hReso[i]->Draw("PE1");
//     fReso[i]->Draw("same");
//     fitLabel[i]->Draw("same");
//     c1->SetLogx();
//     c1->Draw();
//   }    
  


//   // ----- Clean up -----
//   TFile rootfile((dir_+"/jsResponse.root").c_str(),"UPDATE");

//   for(size_t i = 0; i < hRespVsPtGen.size(); ++i) {
//     rootfile.WriteTObject(hScale[i]);
//     rootfile.WriteTObject(hReso[i]);

//     delete hScale[i];
//     delete hReso[i];
//     delete hRespVsPtGen[i];
//   }
//   for(size_t i = 0; i < hRespBins.size(); ++i) {
//     for(size_t j = 0; j < hRespBins[i].size(); ++j) {
//       rootfile.WriteTObject(hRespBins[i][j]);

//       delete hRespBins[i][j];
//     }
//   }
//   for(size_t i = 0; i < hAsymVsPtGen.size(); ++i) {
//     delete hAsymVsPtGen[i];
//   }
//   for(size_t i = 0; i < hAsymBins.size(); ++i) {
//     for(size_t j = 0; j < hAsymBins[i].size(); ++j) {
//       rootfile.WriteTObject(hAsymBins[i][j]);

//       delete hAsymBins[i][j];
//     }
//   }
//   delete leg;
//   delete c1;
//   delete ps;
// }



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


// -------------------------------------------------------------------------------------
double ControlPlotsJetSmearing::gaussian(double *x, double *par) {
  double u = (x[0]-1.)/par[1];
  return par[0]*exp(-0.5*u*u);
}



// --------------------------------------------------
double ControlPlotsJetSmearing::gaussianWidth(double pt) const {
  double a[3];
  for(int i = 0; i < 3; i++) {
    a[i] = scale(i)*param_->GetPars()[i];
  }
  return sqrt( a[0]*a[0] + a[1]*a[1]*pt + a[2]*a[2]*pt*pt );
}


// ------------------------------------------------------------------------
double ControlPlotsJetSmearing::gaussianWidthError(double pt) const {
  // Calculate derivatives
  std::vector<double> ds(3);
  double s = gaussianWidth(pt);
  for(int i = 0; i < 3; i++) {
    ds[i] = scale(i)*param_->GetPars()[i]/s;
    if( i == 1 ) ds[i] *= pt;
    if( i == 2 ) ds[i] *= pt*pt;
  }

  // Calculate variance
  double var = 0.;
  for(int i = 0; i < 3; i++) { // Outer loop over parameters
    for(int j = 0; j < i+1; j++) { // Inner loop over parameters
      int idx = (i*i + i)/2 + j; // Index of (i,j) in covariance vector
      if( i == j ) { // Diagonal terms
	var += ds[i]*ds[i]*scale(i)*scale(i)*param_->GetCovCoeff()[idx];
      } else { // Off-diagonal terms
	var += 2*ds[i]*ds[j]*scale(i)*scale(j)*param_->GetCovCoeff()[idx];
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
  // Zero horizontal error bars
  gStyle->SetErrorX(0);
  
  //  For 'colz' TH2
  gStyle->SetPalette(1);
  
  //  For the canvas
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetCanvasColor(kWhite);
  gStyle->SetCanvasDefH(800); //Height of canvas
  gStyle->SetCanvasDefW(800); //Width of canvas
  gStyle->SetCanvasDefX(0);   //Position on screen
  gStyle->SetCanvasDefY(0);
  
  //  For the frame
  gStyle->SetFrameBorderMode(0);
  gStyle->SetFrameBorderSize(1);
  gStyle->SetFrameFillColor(kBlack);
  gStyle->SetFrameFillStyle(0);
  gStyle->SetFrameLineColor(kBlack);
  gStyle->SetFrameLineStyle(0);
  gStyle->SetFrameLineWidth(1);
  
  //  For the Pad
  gStyle->SetPadBorderMode(0);
  gStyle->SetPadColor(kWhite);
  gStyle->SetPadGridX(false);
  gStyle->SetPadGridY(false);
  gStyle->SetGridColor(0);
  gStyle->SetGridStyle(3);
  gStyle->SetGridWidth(1);
      
  //  Margins
  gStyle->SetPadTopMargin(0.05);
  gStyle->SetPadBottomMargin(0.18);
  gStyle->SetPadLeftMargin(0.19);
  gStyle->SetPadRightMargin(0.04);

  //  For the histo:
  gStyle->SetHistLineColor(kBlack);
  gStyle->SetHistLineStyle(0);
  gStyle->SetHistLineWidth(1);
  
  //  For the statistics box:
  gStyle->SetOptStat(0);
    
  //  For the Global title:
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
  
  //  For the axis
  gStyle->SetAxisColor(1,"XYZ");
  gStyle->SetTickLength(0.03,"XYZ");
  gStyle->SetNdivisions(510,"XYZ");
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetStripDecimals(kFALSE);
    
  //  For the axis labels and titles
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




//   gStyle->SetErrorX(0);
//   gStyle->SetPalette(1);

//   // For the canvas:
//   gStyle->SetCanvasBorderMode(0);
//   gStyle->SetCanvasColor(kWhite);
//   gStyle->SetCanvasDefH(800); //Height of canvas
//   gStyle->SetCanvasDefW(800); //Width of canvas
//   gStyle->SetCanvasDefX(0);   //Position on screen
//   gStyle->SetCanvasDefY(0);

//   // For the frame
//   gStyle->SetFrameBorderMode(0);
//   gStyle->SetFrameBorderSize(1);
//   gStyle->SetFrameFillColor(kBlack);
//   gStyle->SetFrameFillStyle(0);
//   gStyle->SetFrameLineColor(kBlack);
//   gStyle->SetFrameLineStyle(0);
//   gStyle->SetFrameLineWidth(1);

//   // For the Pad:
//   gStyle->SetPadBorderMode(0);
//   gStyle->SetPadColor(kWhite);
//   gStyle->SetPadGridX(false);
//   gStyle->SetPadGridY(false);
//   gStyle->SetGridColor(0);
//   gStyle->SetGridStyle(3);
//   gStyle->SetGridWidth(1);

//   // For the histo:
//   gStyle->SetHistLineColor(kBlack);
//   gStyle->SetHistLineStyle(0);
//   gStyle->SetHistLineWidth(1);
//   //  gStyle->SetMarkerStyle(20);

//   // For the statistics box:
//   gStyle->SetOptStat(0);
//   //  gStyle->SetOptStat("neMR");
//   gStyle->SetStatColor(kWhite);
//   gStyle->SetStatFont(42);
//   gStyle->SetStatFontSize(0.03);
//   gStyle->SetStatTextColor(1);
//   gStyle->SetStatFormat("6.4g");
//   gStyle->SetStatBorderSize(1);
//   gStyle->SetStatX(0.92);              
//   gStyle->SetStatY(0.86);              
//   gStyle->SetStatH(0.16);
//   gStyle->SetStatW(0.22);

//   // For the legend
//   gStyle->SetLegendBorderSize(1);

//   //  Margins
//   // -------------------------------------------
//   gStyle->SetPadTopMargin(0.16);
//   gStyle->SetPadBottomMargin(0.18);
//   gStyle->SetPadLeftMargin(0.18);
//   gStyle->SetPadRightMargin(0.16);

//   // For the Global title:
//   gStyle->SetOptTitle(1);
//   gStyle->SetTitleFont(42,"");
//   gStyle->SetTitleColor(1);
//   gStyle->SetTitleTextColor(1);
//   gStyle->SetTitleFillColor(0);
//   gStyle->SetTitleFontSize(0.12);
//   gStyle->SetTitleAlign(23);
//   gStyle->SetTitleX(0.515);
//   gStyle->SetTitleH(0.06);
//   gStyle->SetTitleXOffset(0);
//   gStyle->SetTitleYOffset(0);
//   gStyle->SetTitleBorderSize(0);

//   // For the axis labels:
//   //  For the axis labels and titles
//   // -------------------------------------------
//   gStyle->SetTitleColor(1,"XYZ");
//   gStyle->SetLabelColor(1,"XYZ");
//   // For the axis labels:
//   gStyle->SetLabelFont(42,"XYZ");
//   gStyle->SetLabelOffset(0.007,"XYZ");
//   gStyle->SetLabelSize(0.045,"XYZ");
  
//   // For the axis titles:
//   gStyle->SetTitleFont(42,"XYZ");
//   gStyle->SetTitleSize(0.06,"XYZ");
//   gStyle->SetTitleXOffset(1.2);
//   gStyle->SetTitleYOffset(1.5);

//   gStyle->SetPadTickX(1);
//   gStyle->SetPadTickY(1);
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


// -------------------------------------------------------------------------------------
void ControlPlotsJetSmearing::setAxisTitles(TH1 *h, const std::string &xTitle, const std::string &xUnit, const std::string &yTitle, bool norm) const {
  // x axis label
  h->SetXTitle(xUnit.length()==0 ? xTitle.c_str() : (xTitle+" ("+xUnit+")").c_str());

  // y axis label
  std::string yAxisTitle;
  if( norm ) yAxisTitle += "1 / N  ";
  if( yTitle == "jets" || yTitle == "events" ) {
    if( yTitle == "jets" )        yAxisTitle += "Jets / ";
    else if( yTitle == "events" ) yAxisTitle += "Events / ";
    yAxisTitle += toString(h->GetBinWidth(1));
    if( xUnit.length() ) yAxisTitle += " "+xUnit;
  } else {
    yAxisTitle = yTitle;
  }
  h->SetYTitle(yAxisTitle.c_str());
}


// -------------------------------------------------------------------------------------
void ControlPlotsJetSmearing::setColor(TH1 *h, int color) const {
  h->SetMarkerColor(color);
  h->SetLineColor(color);
}
