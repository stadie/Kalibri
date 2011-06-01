// $Id: ControlPlotsResolution.cc,v 1.1 2011/05/26 07:42:52 mschrode Exp $

#include "ControlPlotsResolution.h"

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
#include "TH1D.h"
#include "TH2.h"
#include "TH2D.h"
#include "TError.h"
#include "TMath.h"
#include "TPaveText.h"
#include "TRandom3.h"
#include "TROOT.h"

#include "ResolutionFunction.h"
#include "DiJetResolutionEvent.h"
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
ControlPlotsResolution::ControlPlotsResolution(const std::string& configfile, const std::vector<Event*> * data, Parameters * param, const std::string &outDir)
  : data_(data),
    config_(new ConfigFile(configfile.c_str())),
    param_(param),
    respNBins_(120),
    respMin_(0.),
    respMax_(2.),
    dir_(outDir)
{
  // Do not print ROOT message if eps file has been created
  gErrorIgnoreLevel = 1001;

  setGStyle();
  rand_ = new TRandom3(0);

  // Parametrization
  parClass_ = config_->read<std::string>("Parametrization Class","");
  scale_ = bag_of<double>(config_->read<string>("jet parameter scales",""));
  while( static_cast<int>(scale_.size()) < param_->numberOfParameters() ) scale_.push_back(1.);    

  // MC Truth resolution
  truthPar_ = bag_of<double>(config_->read<string>("plots true resolution parameters","3.8663  0.728714  0.  0.224013"));
  std::string formula = config_->read<std::string>("plots true resolution formula","sqrt(((TMath::Sign(1,[0])*sq([0]/x))+(sq([1])*(x^([3]-1))))+sq([2]))");
  truthRes_ = new TF1("ControlPlotsResolution::truthRes",formula.c_str(),0.1,1500.);
  if( truthRes_->GetNpar() != static_cast<int>(truthPar_.size()) ) {
    std::cerr << "ERROR: Number of truth parameters does not match parametrization" << std::endl;
    exit(-1);
  }
  for(size_t i = 0; i < truthPar_.size(); ++i) {
    truthRes_->SetParameter(i,truthPar_.at(i));
  }
  truthRes_->SetLineWidth(1);
  truthRes_->SetLineColor(kRed);

  // Create pt binning
  std::string binning = config_->read<std::string>("plots pt binning","");
  ptBinningVar_ = "ptGen";
  if( binning.find("ptAve") != std::string::npos ) ptBinningVar_ = "ptAve";

  std::cout << "PtBinning variable: " << ptBinningVar_ << std::endl;

  etaMin_ = config_->read<double>("Eta min cut on jet",0.);
  etaMax_  = config_->read<double>("Eta max cut on jet",0.);
  std::string binVar;
  if( ptBinningVar_ == "ptAve" ) binVar = "p^{ave}_{T}";
  else if( ptBinningVar_ == "ptGen" ) binVar = "p^{gen}_{T}";
  titleBins_ = std::vector<TString>(param_->nPtBins());
  for(unsigned int ptBin = 0; ptBin < param_->nPtBins(); ptBin++) {
    titleBins_[ptBin] = toString(param_->ptMin(ptBin))+" < "+binVar+" < "+toString(param_->ptMax(ptBin))+" GeV,  "+toString(etaMin_)+" < |#eta| < "+toString(etaMax_);
  }
  title_ = toString(param_->ptMin())+" < "+binVar+" < "+toString(param_->ptMax())+" GeV,  "+toString(etaMin_)+" < |#eta| < "+toString(etaMax_);

  // Override possible existing root file
  TFile rootfile((dir_+"/jsResponse.root").c_str(),"RECREATE");
  rootfile.Close();

  outNamePrefix_ = dir_+"/"+config_->read<std::string>("plots name prefix","JS")+"_";
  saveAsEps_ = config_->read<bool>("plots save as eps",false);
}

ControlPlotsResolution::~ControlPlotsResolution() {
  delete rand_;
  delete truthRes_;
}



// --------------------------------------------------
void ControlPlotsResolution::makePlots() const {
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
  if( config_->read<bool>("create parallel components plots",false) )
    plotParallelComponents();
  if( config_->read<bool>("create tail plots",false) )
    plotTails();
}



//!  \brief Draw response control plots for events
//!         of type \p DiJetResolutionEvent
// --------------------------------------------------
void ControlPlotsResolution::plotResponse() const
{
  std::cout << "Creating response control plots\n";

  // --- Create histograms of response and spectrum ---------------------
  std::cout << "  Reading parameters" << std::endl;

  std::string param = config_->read<std::string>("Parametrization Class","");
  double rMin = 0.6;
  double rMax = 1.4;

  std::cout << "  Creating histograms" << std::endl;
  std::vector<TH1*> hRespMeasAbs(param_->nPtBins());   // The response ptJet / ptGen absolute entries
  std::vector<TH1*> hRespMeas(param_->nPtBins());      // The response ptJet / ptGen
  std::vector<TH1*> hRespMeasJet1(param_->nPtBins());      // The response ptJet1 / ptGen
  std::vector<TH1*> hRespMeasJet2(param_->nPtBins());      // The response ptJet2 / ptGen
  std::vector<TH1*> hRespSymAbs(param_->nPtBins());   // Symmetrized response
  std::vector<TH1*> hRespMCPtHat(param_->nPtBins());      // The response ptJet / ptGen
  std::vector<TH1*> hRespFitStart(param_->nPtBins());  // The response pdf with start values
  std::vector<TH1*> hRespFit(param_->nPtBins());       // The fitted response pdf
  std::vector<TH1*> hRespFitErrStat(param_->nPtBins());// The fitted response pdf with fitted errors
  std::vector<TH1*> hRespFitStep(param_->nPtBins());   // Step function part of the response pdf
  std::vector<TH1*> hRespFitGaus(param_->nPtBins());   // Gauss part of the response pdf
  std::vector<TH1*> hRespFitSum(param_->nPtBins());    // Sum of step and Gauss part
  std::vector<TH1*> hRespFitBins(param_->nPtBins());       // The fitted response pdf in bins as MC truth
  std::vector<TH1*> hRespRatio(param_->nPtBins());
  std::vector<TH1*> hRespRatioJet1(param_->nPtBins());
  std::vector<TH1*> hRespRatioJet2(param_->nPtBins());
  std::vector<TH1*> hRespRatioFrame(param_->nPtBins());
  std::vector<TH1*> hRespRatioFrameJet1(param_->nPtBins());
  std::vector<TH1*> hRespRatioFrameJet2(param_->nPtBins());
  std::vector<TH1*> hPtGenBins(param_->nPtBins());
  std::vector<TH1*> hPtGenJet1(param_->nPtBins());         // PtGen spectrum
  std::vector<TH1*> hPtGenJet2(param_->nPtBins());         // PtGen spectrum
  std::vector<TH1*> hPtGenAbsBins(param_->nPtBins());
  std::vector<TH1*> hTruthPDF(param_->nPtBins());      // Truth pdf
  std::vector<TH1*> hPtAveAbsBins(param_->nPtBins());
  std::vector<TH1*> hPtGenAsym(param_->nPtBins());
  std::vector<TH1*> hPtAsym(param_->nPtBins());
  std::vector<TH1*> hPtAbsAsym(param_->nPtBins());
  std::vector<TH1*> hPtAsymBiased(param_->nPtBins());
  std::vector<TH1*> hFitPtAsym(param_->nPtBins());
  std::vector<TH1*> hDeltaPtJet12(param_->nPtBins());
  std::vector<TH1*> hDeltaPtJet12Fit(param_->nPtBins());
  std::vector<TH1*> hPJet3(param_->nPtBins());
  std::vector<TH1*> hPJet3Rel(param_->nPtBins());
  std::vector<TH1*> hPJet3GenRel(param_->nPtBins());
  std::vector<TH1*> hPSJ(param_->nPtBins());
  std::vector<TH1*> hPSJRel(param_->nPtBins());
  std::vector<TH1*> hPSJGenRel(param_->nPtBins());
  std::vector<TH1*> hPtJet1(param_->nPtBins());
  std::vector<TH1*> hPtJet2(param_->nPtBins());
  std::vector<TH1*> hPtJet3(param_->nPtBins());
  std::vector<TH1*> hPtJet4(param_->nPtBins());
  std::vector<TH1*> hPtJet3Rel(param_->nPtBins());
  std::vector<TH2*> hPSJvsPtJet4(param_->nPtBins());
  std::vector<TH2*> hPtJet1vs2(param_->nPtBins());
  std::vector<TH1*> hEta(param_->nPtBins());
  std::vector<TH1*> hDeltaPhi12(param_->nPtBins());
  std::vector<TH1*> hNumPU(param_->nPtBins());
  std::vector<TH1*> hWeights(param_->nPtBins());

  TH1* hTruthPDFErrStat = 0;      // Truth pdf
  TH1* hPtGenAbs = 0;         // PtGen spectrum
  TH1* hPtGen = 0;         // PtGen spectrum
  TH1* hPtHat = 0;         // PtHat spectrum
  TH1* hPtAveAbs = 0;       // Dijet spectrum
  TH1* hJESFrame = 0;
  TH1* hJES = 0;
  TH1* hJESJet1 = 0;
  TH1* hJESJet2 = 0;
  TH1* hPtAveCombined = 0;


  for(unsigned int ptBin = 0; ptBin < param_->nPtBins(); ptBin++) {
    std::string name = "hRespMeasAbs_" + toString(ptBin);
    hRespMeasAbs[ptBin] = new TH1D(name.c_str(),titleBins_[ptBin]+";R = p^{jet}_{T} / p^{gen}_{T};dN / dR",
				   respNBins_,respMin_,respMax_);
    hRespMeasAbs[ptBin]->SetLineWidth(2);

    name = "hRespSymAbs_" + toString(ptBin);
    hRespSymAbs[ptBin] = new TH1D(name.c_str(),titleBins_[ptBin]+";R = p^{jet}_{T} / p^{gen}_{T};dN / dR",
				  81,0.,2.);
    

    name = "hRespMeas_" + toString(ptBin);
    hRespMeas[ptBin] = new TH1D(name.c_str(),titleBins_[ptBin]+";R = p^{jet}_{T} / p^{gen}_{T};1 / N  dN / dR",
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
    hRespMCPtHat[ptBin]->SetTitle(titleBins_[ptBin]+";R = p^{jet}_{T} / #hat{p}_{T};1 / N  dN / dR");

    name = "hRespFit_" + toString(ptBin);
    hRespFit[ptBin] = new TH1D(name.c_str(),titleBins_[ptBin]+";R = p^{jet}_{T} / p^{true}_{T};1 / N  dN / dR",
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
    hRespFitStep[ptBin] = new TH1D(name.c_str(),titleBins_[ptBin]+";R = p^{jet}_{T} / p^{true}_{T};1 / N  dN / dR",
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
    hRespRatioJet1[ptBin]->SetXTitle(titleBins_[ptBin]+";R_{1} = p^{jet}_{T,1} / p^{true}_{T,1}");

    name = "hRespRatioJet2_" + toString(ptBin);
    hRespRatioJet2[ptBin] = static_cast<TH1D*>(hRespRatio[ptBin]->Clone(name.c_str()));
    hRespRatioJet2[ptBin]->SetXTitle(titleBins_[ptBin]+";R_{2} = p^{jet}_{T,2} / p^{true}_{T,2}");

    name = "hRespRatioFrame_" + toString(ptBin);
    hRespRatioFrame[ptBin] = static_cast<TH1D*>(hRespRatio[0]->Clone(name.c_str()));
    for(int bin = 1; bin <= hRespRatioFrame[ptBin]->GetNbinsX(); ++bin) {
      hRespRatioFrame[ptBin]->SetBinContent(bin,1.);
    }
    hRespRatioFrame[ptBin]->SetLineStyle(2);
    hRespRatioFrame[ptBin]->GetYaxis()->SetRangeUser(0,2.8);
    hRespRatioFrame[ptBin]->GetXaxis()->SetRangeUser(rMin,rMax);
    name = "hRespRatioFrameJet1_" + toString(ptBin);
    hRespRatioFrameJet1[ptBin] = static_cast<TH1D*>(hRespRatioFrame[ptBin]->Clone(name.c_str()));
    hRespRatioFrameJet1[ptBin]->SetXTitle("R_{1} = p^{jet}_{T,1} / p^{true}_{T,1}");
    name = "hRespRatioFrameJet2_" + toString(ptBin);
    hRespRatioFrameJet2[ptBin] = static_cast<TH1D*>(hRespRatioFrame[ptBin]->Clone(name.c_str()));
    hRespRatioFrameJet2[ptBin]->SetXTitle("R_{2} = p^{jet}_{T,2} / p^{true}_{T,2}");

    name = "hPtGenAbs_" + toString(ptBin);
    hPtGenAbsBins[ptBin] = new TH1D(name.c_str(),titleBins_[ptBin]+";p^{gen}_{T} (GeV);dN / dp^{gen}_{T}  1 / (GeV)",75,param_->ptTrueMin(ptBin),param_->ptTrueMax(ptBin));
    hPtGenAbsBins[ptBin]->SetMarkerStyle(20);
    hPtGenAbsBins[ptBin]->GetXaxis()->SetNdivisions(505);
    hPtGenAbsBins[ptBin]->SetLineWidth(2);

    name = "hPtAveAbs_" + toString(ptBin);
    hPtAveAbsBins[ptBin] = new TH1D(name.c_str(),titleBins_[ptBin]+";p^{ave}_{T} (GeV);dN / dp^{ave}_{T}  1 / (GeV)" ,75,0.8*param_->ptMin(ptBin),1.1*param_->ptMax(ptBin));
    hPtAveAbsBins[ptBin]->SetMarkerStyle(20);
    hPtAveAbsBins[ptBin]->GetXaxis()->SetNdivisions(505);
    hPtAveAbsBins[ptBin]->SetLineWidth(2);
    hPtAveAbsBins[ptBin]->Sumw2();

    name = "hPtGen_" + toString(ptBin);
    hPtGenBins[ptBin] = static_cast<TH1D*>(hPtGenAbsBins[ptBin]->Clone(name.c_str()));
    hPtGenBins[ptBin]->Sumw2();
    hPtGenBins[ptBin]->SetYTitle("1 / N  dN / dp^{gen}_{T}  1 / (GeV)");

    name = "hPtGenJet1_" + toString(ptBin);
    hPtGenJet1[ptBin] = static_cast<TH1D*>(hPtGenBins[ptBin]->Clone(name.c_str()));
    hPtGenJet1[ptBin]->SetTitle(titleBins_[ptBin]+";p^{gen}_{T,1} (GeV);1 / N  dN / dp^{gen}_{T,1}  1 / (GeV)");
    
    name = "hPtGenJet2_" + toString(ptBin);
    hPtGenJet2[ptBin] = static_cast<TH1D*>(hPtGenBins[ptBin]->Clone(name.c_str()));
    hPtGenJet2[ptBin]->SetTitle(titleBins_[ptBin]+";p^{gen}_{T,2} (GeV);1 / N  dN / dp^{gen}_{T,2}  1 / (GeV)");

    name = "hTruthPDF_" + toString(ptBin);
    hTruthPDF[ptBin] = new TH1D(name.c_str(),title_+";p^{true}_{T} (GeV);1 / N  dN / dp^{true}_{T}  1 /  (GeV)",
				5*hPtGenBins[ptBin]->GetNbinsX(),param_->ptTrueMin(ptBin),param_->ptTrueMax(ptBin));
    hTruthPDF[ptBin]->SetLineColor(2);
    hTruthPDF[ptBin]->SetLineWidth(2);
    hTruthPDF[ptBin]->Sumw2();

    name = "hPtGenAsym_" + toString(ptBin);
    hPtGenAsym[ptBin] = new TH1D(name.c_str(),titleBins_[ptBin]+";p^{gen}_{T} asymmetry;",45,-0.2,0.2);
    hPtGenAsym[ptBin]->Sumw2();
    hPtGenAsym[ptBin]->SetLineWidth(2);

    name = "hPtAsym_" + toString(ptBin);
    hPtAsym[ptBin] = new TH1D(name.c_str(),titleBins_[ptBin]+";p_{T} asymmetry;",101,-1.,1.);
    hPtAsym[ptBin]->Sumw2();
    hPtAsym[ptBin]->SetMarkerStyle(20);
    hPtAsym[ptBin]->SetLineWidth(2);

    name = "hPtAbsAsym_" + toString(ptBin);
    hPtAbsAsym[ptBin] = new TH1D(name.c_str(),titleBins_[ptBin]+";p_{T} asymmetry;events",101,-1.,1.);
    hPtAbsAsym[ptBin]->SetMarkerStyle(20);
    hPtAbsAsym[ptBin]->SetLineWidth(2);
    hPtAbsAsym[ptBin]->Sumw2();

    name = "hPtAsymBiased_" + toString(ptBin);
    hPtAsymBiased[ptBin] = new TH1D(name.c_str(),titleBins_[ptBin]+";Biased p_{T} asymmetry;",32,0.,1.);
    hPtAsymBiased[ptBin]->Sumw2();
    hPtAsymBiased[ptBin]->SetMarkerStyle(20);
    hPtAsymBiased[ptBin]->SetLineWidth(2);

    name = "hFitPtAsym_" + toString(ptBin);
    hFitPtAsym[ptBin] = new TH1D(name.c_str(),
				 titleBins_[ptBin]+";p_{T} asymmetry;",
				 5*respNBins_,-0.4,0.4);
    hFitPtAsym[ptBin]->Sumw2();
    hFitPtAsym[ptBin]->SetLineWidth(2);
    hFitPtAsym[ptBin]->SetLineColor(2);

    name = "hDeltaPtJet12_" + toString(ptBin);
    hDeltaPtJet12[ptBin] = new TH1D(name.c_str(),titleBins_[ptBin]+";0.5#upoint|p_{T,1} - p_{T,2}| (GeV);",
				    50,0.,0.5*param_->ptMax(ptBin));
    hDeltaPtJet12[ptBin]->SetMarkerStyle(20);
    hDeltaPtJet12[ptBin]->SetNdivisions(505);
    hDeltaPtJet12[ptBin]->Sumw2();

    name = "hDeltaPtJet12Fit_" + toString(ptBin);
    hDeltaPtJet12Fit[ptBin] = new TH1D(name.c_str(),titleBins_[ptBin]+";0.5#upoint|p_{T,1} - p_{T,2}| (GeV);",
				       500,0.,0.1*param_->ptMax(ptBin));
    hDeltaPtJet12Fit[ptBin]->SetLineColor(2);

    name = "hPJet3_" + toString(ptBin);
    hPJet3[ptBin] = new TH1D(name.c_str(),titleBins_[ptBin]+";p_{||,3} (GeV);Events",
		      50,0.,0.3*param_->ptMax(ptBin));
    hPJet3[ptBin]->SetMarkerStyle(20);

    name = "hPJet3Rel_" + toString(ptBin);
    hPJet3Rel[ptBin] = new TH1D(name.c_str(),titleBins_[ptBin]+";p_{||,3} / p^{ref}_{T};Events",200,0.,0.3);
    hPJet3Rel[ptBin]->SetMarkerStyle(20);

    name = "hPJet3GenRel_" + toString(ptBin);
    hPJet3GenRel[ptBin] = static_cast<TH1D*>(hPJet3Rel[ptBin]->Clone(name.c_str()));
    hPJet3GenRel[ptBin]->SetXTitle(titleBins_[ptBin]+";p_{||,3} / p^{gen,ave}_{T}");

    name = "hPSJ_" + toString(ptBin);
    hPSJ[ptBin] = new TH1D(name.c_str(),titleBins_[ptBin]+";p_{||,>3} (GeV);Events",
		    50,-0.3*param_->ptMin(),0.3*param_->ptMax(ptBin));
    hPSJ[ptBin]->SetMarkerStyle(20);

    name = "hPSJRel_" + toString(ptBin);
    hPSJRel[ptBin] = new TH1D(name.c_str(),titleBins_[ptBin]+";p_{||,>3} / p^{ref}_{T};Events",200,0.,0.1);
    hPSJRel[ptBin]->SetMarkerStyle(20);

    name = "hPSJGenRel_" + toString(ptBin);
    hPSJGenRel[ptBin] = static_cast<TH1D*>(hPSJRel[ptBin]->Clone(name.c_str()));
    hPSJGenRel[ptBin]->SetXTitle(titleBins_[ptBin]+";p_{||,>3} / p^{gen,ave}_{T}");

    name = "hPSJvsPtJet4_" + toString(ptBin);
    hPSJvsPtJet4[ptBin] = new TH2D(name.c_str(),titleBins_[ptBin]+";p_{t,4};p_{||,>3} (GeV) (GeV)",
				   50,0.,0.3*param_->ptMax(ptBin),
				   50,-0.1*param_->ptMin(),0.1*param_->ptMax(ptBin));

    name = "hPtJet1_" + toString(ptBin);
    hPtJet1[ptBin] = new TH1D(name.c_str(),titleBins_[ptBin]+";p_{T,1} (GeV);Jets",70,0.5*param_->ptMin(ptBin),1.3*param_->ptMax(ptBin));
    hPtJet1[ptBin]->SetMarkerStyle(20);
    hPtJet1[ptBin]->GetXaxis()->SetNdivisions(505);
    hPtJet1[ptBin]->SetLineWidth(2);
    hPtJet1[ptBin]->Sumw2();
    setColor(hPtJet1[ptBin],2);

    name = "hPtJet2_" + toString(ptBin);
    hPtJet2[ptBin] = static_cast<TH1D*>(hPtJet1[ptBin]->Clone(name.c_str()));
    hPtJet2[ptBin]->SetXTitle("p_{T,2} (GeV)");
    setColor(hPtJet2[ptBin],4);

    name = "hPtJet3_" + toString(ptBin);
    hPtJet3[ptBin] = new TH1D(name.c_str(),titleBins_[ptBin]+";p_{T,3} (GeV);Events",50,0.,0.2*param_->ptMax(ptBin));

    name = "hPtJet4_" + toString(ptBin);
    hPtJet4[ptBin] = static_cast<TH1D*>(hPtJet3[ptBin]->Clone(name.c_str()));
    hPtJet4[ptBin]->SetXTitle("p_{T,4} (GeV)");

    name = "hPtJet3Rel_" + toString(ptBin);
    hPtJet3Rel[ptBin] = new TH1D(name.c_str(),titleBins_[ptBin]+";p^{rel}_{T,3};Events",50,0.,0.4);

    name = "hPtJet1vs2_" + toString(ptBin);
    hPtJet1vs2[ptBin] = new TH2D(name.c_str(),titleBins_[ptBin]+";p_{T,1} (GeV);p_{T,2} (GeV);",70,0.5*param_->ptMin(ptBin),1.3*param_->ptMax(ptBin),70,0.5*param_->ptMin(ptBin),1.3*param_->ptMax(ptBin));
    hPtJet1vs2[ptBin]->GetXaxis()->SetNdivisions(505);
    hPtJet1vs2[ptBin]->GetYaxis()->SetNdivisions(505);

    name = "hEta_" + toString(ptBin);
    hEta[ptBin] = new TH1D(name.c_str(),titleBins_[ptBin]+";#eta;Jets",101,-5.1,5.1);

    name = "hDeltaPhi12_" + toString(ptBin);
    hDeltaPhi12[ptBin] = new TH1D(name.c_str(),titleBins_[ptBin]+";|#Delta#phi_{12}|;Events",100,0.,3.2);

    name = "hNumPU_" + toString(ptBin);
    hNumPU[ptBin] = new TH1D(name.c_str(),titleBins_[ptBin]+";N(PU Vertices);Events",40,0,40);

    name = "hWeights_" + toString(ptBin);
    hWeights[ptBin] = new TH1D(name.c_str(),titleBins_[ptBin]+";Weight;Events",1000,0.,50.);
  }


  hPtGenAbs = new TH1D("hPtGenAbs",title_+";p^{gen}_{T} (GeV);dN / dp^{gen}_{T}  1 / (GeV)",
		       120,0.8*param_->ptTrueMin(),1.1*param_->ptTrueMax());
  //500,0.8*param_->ptTrueMin(),1.1*param_->ptTrueMax());
  hPtGenAbs->SetMarkerStyle(20);
  hPtGenAbs->GetXaxis()->SetNdivisions(505);
  hPtGenAbs->SetLineWidth(2);

  hPtGen = static_cast<TH1D*>(hPtGenAbs->Clone("hPtGen"));
  hPtGen->SetTitle(title_+";p^{gen}_{T} (GeV);1 / N  dN / dp^{gen}_{T}  1 / (GeV)");
  hPtGen->Sumw2();

  hPtHat = static_cast<TH1D*>(hPtGen->Clone("hPtHat"));
  hPtHat->SetTitle(title_+";#hat{p}_{T} (GeV);1 / N  dN / d#hat{p}_{T}  1 / (GeV)");

  int nbins = static_cast<int>((param_->ptMax()-param_->ptMin())/10.);
  hPtAveAbs = new TH1D("hPtAveAbs",title_+";p^{ave}_{T} (GeV);Events",
		       nbins,param_->ptMin(),param_->ptMax());
  hPtAveAbs->Sumw2();
  hPtAveAbs->SetMarkerStyle(20);

  std::vector<double> edges = findBinEdgesForCombinedSpectrum();
  hPtAveCombined = new TH1D("hPtAveCombined",title_+";p^{ave}_{T} (GeV);Events",
			    edges.size()-1,&(edges.front()));
  hPtAveCombined->Sumw2();
  hPtAveCombined->SetMarkerStyle(20);

  hTruthPDFErrStat = static_cast<TH1D*>(hTruthPDF[0]->Clone("hTruthPDFErrStat"));
  hTruthPDFErrStat->SetFillColor(45);

  hJESFrame = new TH1D("hJESFrame",title_+";p_{T} (GeV);JES",
		       param_->nPtBins(),&(param_->ptBinEdges().front()));
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


  // --- Fill histograms with data quantities --------------
  std::cout << "  Filling histograms with Event data" << std::endl;

  for(DataIt datait = data_->begin(); datait != data_->end(); datait++) {
    if( (*datait)->type() == DiJetResolution ) {
      DiJetResolutionEvent* dijet = static_cast<DiJetResolutionEvent*>(*datait);  
      double weight = dijet->weight();
      double ptAve = dijet->avePt();
      double ptHat = dijet->ptHat();

      // Event wise spectra, all pt
      hPtHat->Fill(ptHat,weight);
      hPtAveAbs->Fill(ptAve,weight);
      hPtAveCombined->Fill(ptAve,weight);

      // Event wise spectra, pt bins
      unsigned int bin = dijet->ptBin();
      hPtAveAbsBins[bin]->Fill(ptAve,weight);
      hPJet3[bin]->Fill(std::abs(dijet->pJ3()),weight);
      hPJet3Rel[bin]->Fill(std::abs(dijet->pJ3()/dijet->ptRef()),weight);
      hPSJ[bin]->Fill(dijet->pSJ()/weight);
      hPSJRel[bin]->Fill(dijet->pSJ()/dijet->ptRef(),weight);
      hPSJvsPtJet4[bin]->Fill(dijet->ptJet4(),dijet->pSJ(),weight);
      double ptGenAve = 0.5*(dijet->jet1()->genPt()+dijet->jet2()->genPt());
      if( ptGenAve > 0 ) {
	hPJet3GenRel[bin]->Fill(dijet->pJ3()/ptGenAve,weight);
	hPSJGenRel[bin]->Fill(dijet->pSJ()/ptGenAve,weight);
      }
      hPtGenJet1[bin]->Fill(dijet->jet1()->genPt(),weight);
      hPtJet1[bin]->Fill(dijet->jet1()->pt(),weight);
      hPtGenJet2[bin]->Fill(dijet->jet2()->genPt(),weight);
      hPtJet2[bin]->Fill(dijet->jet2()->pt(),weight);
      hPtJet3[bin]->Fill(dijet->ptJet3(),weight);
      hPtJet4[bin]->Fill(dijet->ptJet4(),weight);
      hPtJet3Rel[bin]->Fill(dijet->ptJet3()/dijet->ptRef(),weight);
      hDeltaPtJet12[bin]->Fill(0.5*std::abs(dijet->jet1()->pt()-dijet->jet2()->pt()),weight);
      hDeltaPhi12[bin]->Fill(dijet->deltaPhi12(),weight);
      hNumPU[bin]->Fill(dijet->nPU());
      hWeights[bin]->Fill(weight);
      
      // Pt asymmetry variables
      const Jet *j1 = dijet->jet1();
      const Jet *j2 = dijet->jet2();
      hPtJet1vs2[bin]->Fill(j1->pt(),j2->pt(),weight);
      double ht = j1->pt() + j2->pt();
      if( ht > 0 ) {
	double ptAsym = (j1->pt() - j2->pt())/ht;
	hPtAsymBiased[bin]->Fill(ptAsym,weight);
	hPtAsym[bin]->Fill(ptAsym,weight);
	hPtAsym[bin]->Fill(-1.*ptAsym,weight);
	hPtAbsAsym[bin]->Fill(ptAsym,weight);
	hPtAbsAsym[bin]->Fill(-1.*ptAsym,weight);
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
	
	// Response distributions
	hPtGenBins[bin]->Fill(jet->genPt(),weight);
	if( jet->genPt() > 0 ) {
	  hRespMeasAbs[bin]->Fill( jet->pt() / jet->genPt(),weight);
	  hRespMeas[bin]->Fill( jet->pt() / jet->genPt(),weight);
	  if( i == 0 ) {
	    hRespMeasJet1[bin]->Fill( jet->pt() / jet->genPt(),weight);
	  } else if ( i == 1 ) {
	    hRespMeasJet2[bin]->Fill( jet->pt() / jet->genPt(),weight);
	  }

	  // Symmetrized repsonse
	  double sr = (jet->pt()-jet->genPt())/jet->genPt();
	  hRespSymAbs[bin]->Fill(1.+sr,weight);
	  hRespSymAbs[bin]->Fill(1.-sr,weight);
	}
	hPtGenAbsBins[bin]->Fill(jet->genPt(),weight);
	hRespMCPtHat[bin]->Fill(jet->pt()/ptHat,weight);
	hPtGenAsym[bin]->Fill(ptGenAsym,weight);
	hPtGenAsym[bin]->Fill(-1.*ptGenAsym,weight);
	hEta[bin]->Fill(jet->eta(),weight);
      }
    }
  } // End of loop over data

  std::cout << "  Normalising histograms" << std::endl;
  for(unsigned int ptBin = 0; ptBin < param_->nPtBins(); ptBin++) {
    normHist(hRespMeas[ptBin],"width");
    normHist(hRespMeasJet1[ptBin],"width");
    normHist(hRespMeasJet2[ptBin],"width");
    normHist(hRespMCPtHat[ptBin],"width");
    normHist(hPtAsym[ptBin],"width");
    normHist(hPtAsymBiased[ptBin],"width");
    normHist(hDeltaPtJet12[ptBin],"width");
    hDeltaPtJet12[ptBin]->Scale(0.5);
    normHist(hPtGenBins[ptBin],param_->ptTrueMin(),param_->ptTrueMax(),"width");
    normHist(hPtGenJet1[ptBin],param_->ptTrueMin(),param_->ptTrueMax(),"width");
    normHist(hPtGenJet2[ptBin],param_->ptTrueMin(),param_->ptTrueMax(),"width");
  }
  normHist(hPtGen,param_->ptTrueMin(),param_->ptTrueMax(),"width");
  normHist(hPtHat,param_->ptTrueMin(),param_->ptTrueMax(),"width");


  // --- Fill histograms of fitted response ----------------
  std::cout << "  Filling fitted response and truth spectrum" << std::endl;

  // Get parameters
  std::vector<double> fittedPar(param_->numberOfParameters());
  for(int i = 0; i < param_->numberOfParameters(); i++) {
    fittedPar.at(i) = param_->parameters()[i];
  }
  // Mean pt in bin
  std::vector<double> ptMean(param_->nPtBins());
  //  std::vector<double> auxPar = bag_of<double>(config_->read<string>("mean response parameters","1 0"));
  // Loop over ptBins
  for(unsigned int ptBin = 0; ptBin < param_->nPtBins(); ptBin++) {
    const ResolutionFunction& resFunc = param_->resolutionFitPDF(ptBin,1,1);
    ptMean[ptBin] = hPtAveAbsBins[ptBin]->GetMean();

    for(int bin = 1; bin <= hDeltaPtJet12Fit[ptBin]->GetNbinsX(); bin++) {
      double diff = hDeltaPtJet12Fit[ptBin]->GetBinCenter(bin);
//       if( d < dMeas[ptBin] ) {
// 	double u = d/s;
// 	double norm = s*sqrt(M_PI*2.);//*erf(1./sqrt(2.));
// 	hDeltaPtJet12Fit[ptBin]->SetBinContent(bin,exp(-0.5*u*u)/norm);
//       }
      //std::cout << "    Bin " << bin << ": 6.1" << std::endl;
      double pdf = resFunc.pdfPtMeas(2.*diff,0.,ptMean[ptBin]);
      hDeltaPtJet12Fit[ptBin]->SetBinContent(bin,pdf);
    }

    // Interpolated response function
    for(int bin = 1; bin <= hRespFit[ptBin]->GetNbinsX(); bin++) {
      double r = hRespFit[ptBin]->GetBinCenter(bin);
      double val = resFunc.pdfResp(r,ptMean[ptBin]);
      hRespFit[ptBin]->SetBinContent(bin,val);
      hRespFitErrStat[ptBin]->SetBinContent(bin,val);
      //	hRespFitErrStat[ptBin]->SetBinError(bin,smearData->pdfRespError(r,ptMean[ptBin]));

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
      hFitPtAsym[ptBin]->SetBinContent(bin,resFunc.pdfDijetAsym(a,hPtGenBins[ptBin]->GetMean()));
    }

//     for(int ptAveBin = 1; ptAveBin <= hPtAveAbsBins[ptBin]->GetNbinsX(); ++ptAveBin) {
//       for(int bin = 1; bin <= hFitPtAsym[ptBin]->GetNbinsX(); bin++) {
// 	double a = hFitPtAsym[ptBin]->GetBinCenter(bin);
// 	hFitPtAsym[ptBin]->Fill(a,hPtAveAbsBins[ptBin]->GetBinContent(ptAveBin)*resFunc.pdfDijetAsym(a,hPtAveAbsBins[ptBin]->GetBinCenter(ptAveBin)));
//       }
//     }
//     normHist(hFitPtAsym[ptBin],"width");

    // Interpolated fit function with start values
    // Copy start values into parameter array
    for(int i = 0; i < param_->numberOfParameters(); i++) {
      param_->parameters()[i] = param_->jetStartPar(i);
    }
    // Plot response function
    for(int bin = 1; bin <= hRespFitStart[ptBin]->GetNbinsX(); bin++) {
      double r = hRespFitStart[ptBin]->GetBinCenter(bin);
      hRespFitStart[ptBin]->SetBinContent(bin,resFunc.pdfResp(r,ptMean[ptBin]));
    }

    // Copy back fitted values into parameter array
    for(int i = 0; i < param_->numberOfParameters(); i++) {
      param_->parameters()[i] = fittedPar.at(i);
    }

    
    // --- Fill histograms of fitted truth spectrum -----------
    
    // Fill histogram of assumed dijet truth pdf
    for(int bin = 1; bin <= hTruthPDF[ptBin]->GetNbinsX(); bin++) {
      double t = hTruthPDF[ptBin]->GetBinCenter(bin);
      hTruthPDF[ptBin]->SetBinContent(bin,resFunc.pdfPtTrue(t));
    }
  } // End of loop over ptBins


  std::cout << "  Setting axis ranges" << std::endl;

  // --- Set x-axis ranges ----------------------------------
  for(unsigned int ptBin = 0; ptBin < param_->nPtBins(); ptBin++) {
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
  for(unsigned int ptBin = 0; ptBin < param_->nPtBins(); ptBin++) {
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

    setYRange(hPtGenBins[ptBin], 0.5, 100.);
    setYRange(hPtGenJet1[ptBin], 0.5, 100.);
    setYRange(hPtGenJet2[ptBin], 0.5, 100.);
  }
  setYRange(hPtGen, 0.5, 100.);
  setYRange(hPtHat, 0.5, 100.);



  // --- Plot histograms -----------------------------------
  std::cout << "  Plotting histograms" << std::endl;

  // Label bins
  std::vector<TLegend*> legPtRangeAndCenters(param_->nPtBins());
  for(unsigned int ptBin = 0; ptBin < param_->nPtBins(); ptBin++) {
    legPtRangeAndCenters[ptBin] = new TLegend(0.23,0.75,0.8,0.9);
    legPtRangeAndCenters[ptBin]->SetBorderSize(0);
    legPtRangeAndCenters[ptBin]->SetFillColor(0);
    legPtRangeAndCenters[ptBin]->SetTextFont(42);

    std::string binVar;
    if( ptBinningVar_ == "ptAve" ) binVar = "p^{ave}_{T}";
    else if( ptBinningVar_ == "pt" ) binVar = "p_{T}";
    else if( ptBinningVar_ == "ptGen" ) binVar = "p^{gen}_{T}";

    std::string label = toString(param_->ptMin(ptBin))
      + " < " + binVar + " < "
      + toString(param_->ptMax(ptBin))
      + " GeV";
    legPtRangeAndCenters[ptBin]->AddEntry(hRespMeas[ptBin],label.c_str(),"P");
    if( ptBinningVar_ == "ptAve" ) label = "p_{T} = " + toString(hPtAveAbsBins[ptBin]->GetMean()) + " GeV";
    else if( ptBinningVar_ == "ptGen" ) label = "p_{T} = " + toString(hPtGenAbsBins[ptBin]->GetMean()) + " GeV";    legPtRangeAndCenters[ptBin]->AddEntry(hRespFit[ptBin],"Fit","L");
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
  for(unsigned int ptBin = 0; ptBin < param_->nPtBins(); ptBin++) {
    ps->NewPage();
    c1->cd();
    hPtAveAbsBins[ptBin]->Draw("PE1");
    c1->SetLogy(1);
    c1->Draw();

    ps->NewPage();
    c1->cd();
    hPtGenBins[ptBin]->Draw("PE1");
    c1->SetLogy(1);
    c1->Draw();

    ps->NewPage();
    c1->cd();
    double tmpMin = 0.;
    double tmpMax = 0.;
    findYRange(hPtGenBins[ptBin],tmpMin,tmpMax);
    tmpMax *= 1.4;
    hPtGenBins[ptBin]->GetYaxis()->SetRangeUser(0.,tmpMax);
    hPtGenBins[ptBin]->Draw("PE1");
    hTruthPDF[ptBin]->Draw("Lsame");
    gPad->RedrawAxis();
    c1->SetLogy(0);
    c1->Draw();

    ps->NewPage();
    c1->cd();
    hPtGenJet1[ptBin]->Draw("PE1");
    hTruthPDF[ptBin]->Draw("Lsame");
    gPad->RedrawAxis();
    c1->SetLogy(1);
    c1->Draw();

    ps->NewPage();
    c1->cd();
    tmpMin = 0.;
    tmpMax = 0.;
    findYRange(hPtGenJet1[ptBin],tmpMin,tmpMax);
    tmpMax *= 1.4;
    hPtGenJet1[ptBin]->GetYaxis()->SetRangeUser(0.,tmpMax);
    hPtGenJet1[ptBin]->Draw("PE1");
    hTruthPDF[ptBin]->Draw("Lsame");
    gPad->RedrawAxis();
    c1->SetLogy(0);
    c1->SetLogx(0);
    c1->Draw();
    c1->SetLogx(0);

    ps->NewPage();
    c1->cd();
    hPtGenJet2[ptBin]->Draw("PE1");
    hTruthPDF[ptBin]->Draw("Lsame");
    gPad->RedrawAxis();
    c1->SetLogy();
    c1->Draw();

    ps->NewPage();
    c1->cd();
    tmpMin = 0.;
    tmpMax = 0.;
    findYRange(hPtGenJet2[ptBin],tmpMin,tmpMax);
    tmpMax *= 1.4;
    hPtGenJet2[ptBin]->GetYaxis()->SetRangeUser(0.,tmpMax);
    hPtGenJet2[ptBin]->Draw("PE1");
    hTruthPDF[ptBin]->Draw("Lsame");
    gPad->RedrawAxis();
    c1->SetLogy(0);
    c1->SetLogx(0);
    c1->Draw();
    c1->SetLogx(0);
    
    ps->NewPage();
    c1->cd();
    hPtJet1[ptBin]->Draw("PE1");
    hPtJet2[ptBin]->Draw("PE1same");
    c1->Draw();

    ResolutionFunction resFunc = param_->resolutionFitPDF(ptBin,1,1);
    if( hDeltaPtJet12[ptBin]->Integral() ) {
      int maxBin = hDeltaPtJet12[ptBin]->FindBin(resFunc.dMeasMax());
      double norm = hDeltaPtJet12[ptBin]->Integral(1,maxBin,"width");
      hDeltaPtJet12[ptBin]->Scale(1./norm);
    }
    TF1 *func = new TF1("func","gaus",0.,resFunc.dMeasMax());
    func->FixParameter(1,0.);
    func->SetParameter(2,hDeltaPtJet12[ptBin]->GetRMS());
    func->SetParameter(0,1./sqrt(2.*M_PI)/func->GetParameter(2));
    hDeltaPtJet12[ptBin]->Fit(func,"I0QBR");
    func->SetLineColor(kBlue);
    func->SetLineWidth(1);
    func->SetLineStyle(2);

    hDeltaPtJet12Fit[ptBin]->Scale(func->Integral(0.,resFunc.dMeasMax()));

//     std::cout << "Int(pdf)   " << hDeltaPtJet12Fit[ptBin]->Integral("width") << std::endl;
//     std::cout << "Int(fit)   " << func->Integral(0.,resFunc.dMeasMax()) << std::endl;
//     std::cout << "Int(hist)  " << hDeltaPtJet12[ptBin]->Integral("width") << std::endl;
//     std::cout << std::endl;

    std::cout << "    DeltaMeas(pdf) \t" << (param_->parameters()[ptBin])/sqrt(2.) << " +/- " << (param_->errors()[ptBin])/sqrt(2.) << std::endl;
    std::cout << "    DeltaMeas(fit) \t" << func->GetParameter(2) << " +/- " << func->GetParError(2) << std::endl;
    

    ps->NewPage();
    c1->cd();
    hDeltaPtJet12[ptBin]->Draw("PE1");
    hDeltaPtJet12Fit[ptBin]->Draw("Lsame");
    func->Draw("same");
    gPad->RedrawAxis();
    c1->SetLogy(0);
    c1->Draw();

    ps->NewPage();
    c1->cd();
    hDeltaPtJet12[ptBin]->Draw("PE1");
    hDeltaPtJet12Fit[ptBin]->Draw("Lsame");
    func->SetRange(0,hDeltaPtJet12[ptBin]->GetXaxis()->GetBinUpEdge(hDeltaPtJet12[ptBin]->GetNbinsX()));
    func->Draw("same");
    std::vector<TLine*> lines;
    for(int k = 0; k < 6; ++k) {
      double x = (2.+k)*std::abs(func->GetParameter(2));//*sqrt(2.);
      if( x > hDeltaPtJet12[ptBin]->GetXaxis()->GetBinUpEdge(hDeltaPtJet12[ptBin]->GetNbinsX()) ) break;
      TLine *line = new TLine(x,3E-5,x,hDeltaPtJet12[ptBin]->GetMaximum());
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
    TF1 *fitAsym = new TF1("fitAsym","gaus",hPtAsym[ptBin]->GetMean()-2.*resFunc.dMeasMax()/ptMean[ptBin],hPtAsym[ptBin]->GetMean()+2.*resFunc.dMeasMax()/ptMean[ptBin]);
    fitAsym->SetParameter(0,1./sqrt(2.*M_PI)/0.1);
    fitAsym->SetParameter(1,0.);
    fitAsym->SetParameter(2,0.1);
    hPtAsym[ptBin]->Fit(fitAsym,"0QIRB");
    fitAsym->SetLineStyle(2);
    fitAsym->SetLineWidth(2);
    hPtAsym[ptBin]->GetXaxis()->SetRangeUser(-0.2,0.2);
    hPtAsym[ptBin]->Draw("PE1");
    fitAsym->Draw("same");
    hFitPtAsym[ptBin]->Draw("HISTsame");
    legPtRangeAndCenters[ptBin]->Draw("same");
    gPad->RedrawAxis();
    c1->SetLogy(logy);
    c1->Draw();

    std::cout << "    Pdf  \t" << (param_->parameters()[ptBin])/ptMean[ptBin] << " +/- " << (param_->errors()[ptBin])/ptMean[ptBin] << std::endl;
    std::cout << "    Asym \t" << sqrt(2.)*fitAsym->GetParameter(2) << " +/- " << sqrt(2.)*fitAsym->GetParError(2) << std::endl;
    delete fitAsym;

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
    if( hRespMeas[ptBin]->Fit("gaus","I0Q") == 0 ) {
      hRespMeas[ptBin]->GetFunction("gaus")->Draw("same");
      //std::cout << "Res  " << hRespMeas[ptBin]->GetFunction("gaus")->GetParameter(2) << " +/- " << hRespMeas[ptBin]->GetFunction("gaus")->GetParError(2) << std::endl;
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
    hRespRatioFrame[ptBin]->Draw("L");
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

//     ps->NewPage();
//     c1->cd();
//     hRespMeasJet1[ptBin]->Draw("PE1");
//     hRespFitBins[ptBin]->Draw("HISTsame");
//     gPad->RedrawAxis();
//     legPtRangeAndCenters[ptBin]->Draw("same");
//     c1->SetLogy(0);
//     c1->Draw();

//     ps->NewPage();
//     c1->cd();
//     hRespRatioFrameJet1[ptBin]->Draw("L");
//     hRespRatioJet1[ptBin]->Draw("PE1same");
//     c1->SetLogy(0);
//     c1->Draw();

//     ps->NewPage();
//     c1->cd();
//     hRespMeasJet2[ptBin]->Draw("PE1");
//     hRespFit[ptBin]->Draw("Lsame");
//     gPad->RedrawAxis();
//     legPtRangeAndCenters[ptBin]->Draw("same");
//     c1->Draw();

//     ps->NewPage();
//     c1->cd();
//     hRespMeasJet2[ptBin]->Draw("PE1");
//     hRespFitBins[ptBin]->Draw("HISTsame");
//     gPad->RedrawAxis();
//     legPtRangeAndCenters[ptBin]->Draw("same");
//     c1->SetLogy(0);
//     c1->Draw();

//     ps->NewPage();
//     c1->cd();
//     hRespRatioFrameJet2[ptBin]->Draw("L");
//     hRespRatioJet2[ptBin]->Draw("PE1same");
//     c1->SetLogy(0);
//     c1->Draw();


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
    c1->SetLogy(logy);
    c1->Draw();

//     double ptMean = 1.;
//     if( ptBinningVar_ == "ptGen" ) ptMean = hPtGenAbsBins[ptBin]->GetMean();
//     else if( ptBinningVar_ == "ptDijet" ) ptMean = hPtAveAbs->GetMean();
//     double par = scale(0)*param_->parameters()[0]/ptMean;
//     double parErr = scale(0)*param_->errors()[0]/ptMean;
//     double parMCTruth = 0.5*(fitResp[0]->GetParameter(2)+fitResp[1]->GetParameter(2));
//     double parMCTruthErr = 0.5*(fitResp[0]->GetParError(2)+fitResp[1]->GetParError(2));

//     std::cout << std::endl;
//     std::cout << "Pt          = " << ptMean << std::endl;
//     std::cout << "Resp        = " << parMCTruth << " +/- " << parMCTruthErr << std::endl;
//     std::cout << "Fit (Resp)  = " << par << " +/- " << parErr << std::endl;
//     std::cout << "Dev         = " << par / parMCTruth - 1. << std::endl;
//     std::cout << "Asym        = " << std::abs(fitAsym->GetParameter(2)) << " +/- " << fitAsym->GetParError(2) << std::endl;
//     std::cout << "Fit (Asym)  = " << par/sqrt(2.) << " +/- " << parErr/sqrt(2.) << std::endl;
//     std::cout << "Dev         = " << par/sqrt(2.)/std::abs(fitAsym->GetParameter(2)) - 1. << std::endl;

    for(int jetIdx = 0; jetIdx < 2; ++jetIdx) {
      delete fitResp[jetIdx];
      delete fitRespLine[jetIdx];
    }
    
    ps->NewPage();
    c1->cd();
    hPtAsymBiased[ptBin]->Draw("PE1");
    gPad->RedrawAxis();
    c1->SetLogy(logy);
    c1->Draw();
  }

  ps->NewPage();
  c1->cd();
  hPtGenAbs->Draw("PE1");
  c1->SetLogy();
  c1->Draw();

  ps->NewPage();
  c1->cd();
  hPtGen->Draw("PE1");
  gPad->RedrawAxis();
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
  gPad->RedrawAxis();
  c1->SetLogy(0);
  c1->SetLogx(0);
  c1->Draw();
  c1->SetLogx(0);

  ps->NewPage();
  c1->cd();
  hPtAveAbs->Draw("PE1");
  c1->SetLogy();
  c1->Draw();


  std::cout << "  Writing histos to ROOT file" << std::endl;

  // Write histos to root file
  TFile rootfile((dir_+"/jsResponse.root").c_str(),"UPDATE");
  for(unsigned int ptBin = 0; ptBin < param_->nPtBins(); ptBin++) {
    rootfile.WriteTObject(hRespMeasAbs[ptBin]);
    rootfile.WriteTObject(hRespMeas[ptBin]);
    rootfile.WriteTObject(hRespMeasJet1[ptBin]);
    rootfile.WriteTObject(hRespMeasJet2[ptBin]);
    rootfile.WriteTObject(hRespSymAbs[ptBin]);
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
    rootfile.WriteTObject(hPtGenBins[ptBin]);
    rootfile.WriteTObject(hPtAsym[ptBin]);
    rootfile.WriteTObject(hPtAbsAsym[ptBin]);
    rootfile.WriteTObject(hPtAsymBiased[ptBin]);
    rootfile.WriteTObject(hFitPtAsym[ptBin]);

    rootfile.WriteTObject(hDeltaPtJet12Fit[ptBin]);
    rootfile.WriteTObject(hDeltaPtJet12[ptBin]);
    rootfile.WriteTObject(hPJet3[ptBin]);
    rootfile.WriteTObject(hPJet3Rel[ptBin]);
    rootfile.WriteTObject(hPJet3GenRel[ptBin]);
    rootfile.WriteTObject(hPSJ[ptBin]);
    rootfile.WriteTObject(hPSJRel[ptBin]);
    rootfile.WriteTObject(hPSJGenRel[ptBin]);
    rootfile.WriteTObject(hPtJet1[ptBin]);
    rootfile.WriteTObject(hPtJet2[ptBin]);
    rootfile.WriteTObject(hPtJet3[ptBin]);
    rootfile.WriteTObject(hPtJet4[ptBin]);
    rootfile.WriteTObject(hPtJet3Rel[ptBin]);
    rootfile.WriteTObject(hPSJvsPtJet4[ptBin]);
    rootfile.WriteTObject(hPtJet1vs2[ptBin]);
    rootfile.WriteTObject(hPtGenAbsBins[ptBin]);
    rootfile.WriteTObject(hPtAveAbsBins[ptBin]);
    rootfile.WriteTObject(hPtGenBins[ptBin]);
    rootfile.WriteTObject(hPtGenJet1[ptBin]);
    rootfile.WriteTObject(hPtGenJet2[ptBin]);
    rootfile.WriteTObject(hTruthPDF[ptBin]);
    rootfile.WriteTObject(hEta[ptBin]);
    rootfile.WriteTObject(hDeltaPhi12[ptBin]);
    rootfile.WriteTObject(hWeights[ptBin]);
    rootfile.WriteTObject(hNumPU[ptBin]);
  }
  rootfile.WriteTObject(hPtGenAbs);
  rootfile.WriteTObject(hPtGen);
  rootfile.WriteTObject(hPtHat);
  rootfile.WriteTObject(hPtAveAbs);
  rootfile.WriteTObject(hPtAveCombined);
  rootfile.WriteTObject(hTruthPDFErrStat);

  rootfile.Close();


  std::cout << "  Cleaning up" << std::endl;

  // --- Clean up ------------------------------------------
  for(unsigned int ptBin = 0; ptBin < param_->nPtBins(); ptBin++) {
    delete hRespMeasAbs[ptBin];
    delete hRespMeas[ptBin];
    delete hRespMeasJet1[ptBin];
    delete hRespMeasJet2[ptBin];
    delete hRespSymAbs[ptBin];
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
    delete hRespRatioFrame[ptBin];
    delete hRespRatioFrameJet1[ptBin];
    delete hRespRatioFrameJet2[ptBin];
    delete hPtGenAbsBins[ptBin];
    delete hPtAveAbsBins[ptBin];
    delete hPtGenBins[ptBin];
    delete hPtGenJet1[ptBin];
    delete hPtGenJet2[ptBin];
    delete hPtGenAsym[ptBin];
    delete hPtAsym[ptBin];
    delete hPtAbsAsym[ptBin];
    delete hPtAsymBiased[ptBin];
    delete hFitPtAsym[ptBin];
    delete hDeltaPtJet12[ptBin];
    delete hDeltaPtJet12Fit[ptBin];
    delete hPJet3[ptBin];
    delete hPJet3Rel[ptBin];
    delete hPJet3GenRel[ptBin];
    delete hPSJ[ptBin];
    delete hPSJRel[ptBin];
    delete hPSJGenRel[ptBin];
    delete hPtJet1[ptBin];
    delete hPtJet2[ptBin];
    delete hPtJet3[ptBin];
    delete hPtJet4[ptBin];
    delete hPtJet3Rel[ptBin];
    delete hPSJvsPtJet4[ptBin];
    delete hPtJet1vs2[ptBin];
    delete hTruthPDF[ptBin];
    delete hEta[ptBin];
    delete hDeltaPhi12[ptBin];
    delete hWeights[ptBin];
    delete hNumPU[ptBin];
  }
  delete hPtGenAbs;
  delete hPtGen;
  delete hPtHat;
  delete hPtAveAbs;
  delete hPtAveCombined;
  delete hTruthPDFErrStat;
  delete legFitStart;
  delete c1;
  delete ps;
}


// --------------------------------------------------
std::vector<double> ControlPlotsResolution::findBinEdgesForCombinedSpectrum() const {
  // Divide bins
  std::vector<double> binEdges;
  double edge = 45.;
  while( edge < 1000. ) {
    binEdges.push_back(edge);
    edge += 2.5;
  }

  return binEdges;
}



// --------------------------------------------------
void ControlPlotsResolution::plotTails() const {
  std::cout << "Plotting tails" << std::endl;

  std::cout << "  Creating response histograms" << std::endl;


  // --- Creating histograms --------------
  std::cout << "  Creating histograms" << std::endl;

  // Response histograms
  std::vector<TH1*> hRespPtBins(param_->nPtBins());

  for(unsigned int ptBin = 0; ptBin < param_->nPtBins(); ++ptBin) {
    TString name = "hResp_Pt"+toString(ptBin);
    hRespPtBins[ptBin] = new TH1D(name,titleBins_[ptBin]+";R = p_{T} / p^{gen}_{T};Probability",51,0.,2.);
    hRespPtBins[ptBin]->Sumw2();
    hRespPtBins[ptBin]->SetMarkerStyle(20);
  }


  std::cout << "  Filling histograms" << std::endl;

  for(DataIt datait = data_->begin(); datait != data_->end(); datait++) {
    // Select DiJet events
    if( (*datait)->type() == DiJetResolution ) {
      DiJetResolutionEvent* dijet = static_cast<DiJetResolutionEvent*>(*datait);  
      double weight = dijet->weight();
      unsigned int ptBin = dijet->ptBin();
      
      for(int i = 0; i < 2; i++) {        // Loop over both jets
	const Jet * jet = dijet->jet(i);
	if( jet->genPt() > 0 ) {
	  hRespPtBins[ptBin]->Fill( jet->pt() / jet->genPt(),weight);
	}
      }
    }
  } // End of loop over data

  for(unsigned int ptBin = 0; ptBin < param_->nPtBins(); ptBin++) {
    normHist(hRespPtBins[ptBin],"width");
  }
  


//   std::cout << "  Subtracting Gaussian" << std::endl;
//   std::vector<TH1*> hTails(param_->nPtBins());
//   std::vector<TH1*> hTailsClean(param_->nPtBins());
//   std::vector<TF1*> fGauss(param_->nPtBins());
//   std::vector<TH1*> hTailsInter(param_->nPtBins());
//   for(unsigned int ptBin = 0; ptBin < param_->nPtBins(); ++ptBin) {
//     subtractGaussian(hResp[ptBin],hTails[ptBin],hTailsClean[ptBin],fGauss[ptBin]);
//     interpolateTails(hTailsClean[ptBin],hTailsInter[ptBin]);
//   }


  std::cout << "  Writing histos to file" << std::endl;
  // Write histos to file
  //  TString outName = outNamePrefix_+"Tails_";
  TPostScript * const ps = new TPostScript((dir_+"/jsTails.ps").c_str(),111);
  TCanvas *c1 = new TCanvas("c1","Tails",0,0,600,600);

  for(unsigned int ptBin = 0; ptBin < param_->nPtBins(); ptBin++) {
    ps->NewPage();
    c1->cd();
    hRespPtBins[ptBin]->Draw("PE1");
    c1->SetLogy();
    c1->Draw();
  }
  ps->Close();



//   // Full response distributions (all jets)
//   for(unsigned int ptBin = 0; ptBin < param_->nPtBins(); ++ptBin) {
//     ps->NewPage();
//     c1->cd();
//     hResp[ptBin]->Draw("PE1");
//     c1->SetLogy(1);
//     if( saveAsEps_ ) c1->SaveAs(outName+"Response_"+toString(ptBin)+".eps","eps");
//     else c1->Draw();
//   }
//   // Response distributions, Gaussian fit, tails (all jets)
//   for(unsigned int ptBin = 0; ptBin < param_->nPtBins(); ++ptBin) {
//     ps->NewPage();
//     c1->cd();
//     hResp[ptBin]->Draw("PE1");
//     fGauss[ptBin]->Draw("same");
//     hTails[ptBin]->Draw("same");
//     hTailsClean[ptBin]->Draw("same");
//     c1->SetLogy(1);
//     if( saveAsEps_ ) c1->SaveAs(outName+"GaussAndTails_"+toString(ptBin)+".eps","eps");
//     else c1->Draw();
//   }
//   // Tails (linear) (all jets)
//   for(unsigned int ptBin = 0; ptBin < param_->nPtBins(); ++ptBin) {
//     ps->NewPage();
//     c1->cd();
//     hTails[ptBin]->GetYaxis()->SetRangeUser(0,1.1*hTails[ptBin]->GetMaximum());
//     hTails[ptBin]->Draw("HIST");
//     hTailsClean[ptBin]->Draw("HISTsame");
//     c1->SetLogy(0);
//     if( saveAsEps_ ) c1->SaveAs(outName+"Tails_"+toString(ptBin)+".eps","eps");
//     else c1->Draw();
//   }
//   // Cleaned tails (linear) (all jets)
//   for(unsigned int ptBin = 0; ptBin < param_->nPtBins(); ++ptBin) {
//     ps->NewPage();
//     c1->cd();
//     hTailsClean[ptBin]->GetYaxis()->SetRangeUser(0.,1.1*hTailsClean[ptBin]->GetMaximum());
//     hTailsClean[ptBin]->Draw("HIST");
//     c1->SetLogy(0);
//     if( saveAsEps_ ) c1->SaveAs(outName+"TailsClean_"+toString(ptBin)+".eps","eps");
//     else c1->Draw();
//   }
//   // Cleaned tails (linear, interpolated) (all jets)
//   for(unsigned int ptBin = 0; ptBin < param_->nPtBins(); ++ptBin) {
//     ps->NewPage();
//     c1->cd();
//     hTailsInter[ptBin]->GetYaxis()->SetRangeUser(0.,1.1*hTailsInter[ptBin]->GetMaximum());
//     hTailsInter[ptBin]->Draw("HIST");
//     c1->SetLogy(0);
//     if( saveAsEps_ ) c1->SaveAs(outName+"TailsInter_"+toString(ptBin)+".eps","eps");
//     else c1->Draw();
//   }

//   ps->Close();



  std::cout << "  Writing histos to ROOT file" << std::endl;

  // Write histos to root file
  TFile rootfile((dir_+"/jsResponse.root").c_str(),"UPDATE");
  for(unsigned int ptBin = 0; ptBin < param_->nPtBins(); ptBin++) {
    rootfile.WriteTObject(hRespPtBins[ptBin]);
  }
  rootfile.Close();
  

  std::cout << "  Cleaning up" << std::endl;

  // --- Clean up ------------------------------------------
  for(unsigned int ptBin = 0; ptBin < param_->nPtBins(); ptBin++) {
    delete hRespPtBins[ptBin];
  }
  delete c1;
  delete ps;


//   // Write histos to root file
//   TFile rootfile((dir_+"/jsResponse.root").c_str(),"UPDATE");
//   for(unsigned int ptBin = 0; ptBin < param_->nPtBins(); ptBin++) {
//     rootfile.WriteTObject(hResp[ptBin]);
//     rootfile.WriteTObject(hTails[ptBin]);
//     rootfile.WriteTObject(hTailsClean[ptBin]);
//     rootfile.WriteTObject(hTailsInter[ptBin]);
//     rootfile.WriteTObject(fGauss[ptBin]);
//   }
//   rootfile.Close();

  
//   for(unsigned int ptBin = 0; ptBin < param_->nPtBins(); ++ptBin) {
//     delete hResp[ptBin];
//     delete hTails[ptBin];
//     delete hTailsClean[ptBin];
//     delete hTailsInter[ptBin];
//     delete fGauss[ptBin];
//   }
//   delete c1;
//   delete ps;
}



// --------------------------------------------------
void ControlPlotsResolution::subtractGaussian(const TH1* hResp, TH1* &hTail, TH1* &hTailClean, TF1* &gauss) const {

  TString name = hResp->GetName();
  name += "_Tails";
  hTail = static_cast<TH1D*>(hResp->Clone(name));
  hTail->SetLineColor(kBlue);
  hTail->SetMarkerColor(kBlue);
  hTail->SetMarkerStyle(hResp->GetMarkerStyle()+1);

  name = hResp->GetName();
  name += "_TailsClean";
  hTailClean = static_cast<TH1D*>(hResp->Clone(name));
  hTailClean->Reset();
  hTailClean->SetLineColor(kBlue+3);
  hTailClean->SetMarkerColor(kBlue+3);
  hTailClean->SetMarkerStyle(hResp->GetMarkerStyle()+2);

  name = hResp->GetName();
  name += "_GaussFit";
  gauss = new TF1(name,"gaus",0.,2.);
  gauss->SetLineWidth(1);
  gauss->SetLineColor(kRed);

   double mean = hTail->GetMean();
   double sigma = hTail->GetRMS();
   if( hTail->Fit(gauss,"I0Q","",mean-2.*sigma,mean+2.*sigma) == 0 ) {
     mean = gauss->GetParameter(1);
     sigma = gauss->GetParameter(2);
     int binMin = hTail->FindBin(mean-5.*sigma);
     int binMax = hTail->FindBin(mean+5.*sigma);
     for(int bin = binMin; bin <= binMax; ++bin) {
       double min = hTail->GetXaxis()->GetBinLowEdge(bin);
       double max = hTail->GetXaxis()->GetBinUpEdge(bin);
       double gaussPdf = gauss->Integral(min,max)/hTail->GetBinWidth(1);
       double tailPdf = hTail->GetBinContent(bin) - gaussPdf;
       hTail->SetBinContent(bin,tailPdf);
     }
     for(int bin = 1; bin <= hTailClean->GetNbinsX(); ++bin) {
       if( hTail->GetBinContent(bin) > 0.4*gauss->Eval(hTail->GetBinCenter(bin)) ) {
	 hTailClean->SetBinContent(bin,hTail->GetBinContent(bin));
	 hTailClean->SetBinError(bin,hTail->GetBinError(bin));
       }
     }
   }
}



// --------------------------------------------------
void ControlPlotsResolution::interpolateTails(TH1* hTails, TH1* &hInter) const {
  TString name = hTails->GetName();
  name += "Interpolated";
  hInter = new TH1D(name,"",10*hTails->GetNbinsX(),hTails->GetXaxis()->GetXmin(),hTails->GetXaxis()->GetXmax());
  hInter->SetTitle(hTails->GetTitle());
  hInter->SetXTitle(hTails->GetXaxis()->GetTitle());
  hInter->SetYTitle(hTails->GetYaxis()->GetTitle());
  hInter->SetLineColor(hTails->GetLineColor());
  for(int bin = 1; bin <= hInter->GetNbinsX(); ++bin) {
    hInter->SetBinContent(bin,hTails->Interpolate(hInter->GetBinCenter(bin)));
  }
}



// --------------------------------------------------
void ControlPlotsResolution::plotParallelComponents() const {
  std::cout << "Plotting parallel components" << std::endl;

  std::cout << "  Creating histograms" << std::endl;
  std::vector<TH1*> hPtGenAve(param_->nPtBins());
  std::vector< std::vector<TH1*> > hPtGen(param_->nPtBins());
  std::vector< std::vector<TH2*> > hPpOverPtGenVsDeltaPhi12(param_->nPtBins());
  std::vector<TH1*> hDeltaPhi(param_->nPtBins());
  std::vector<TH1*> hPPhi(param_->nPtBins());

  std::vector<TH1*> hResp(param_->nPtBins());
  std::vector<TH1*> hPtAsym(param_->nPtBins());
  std::vector<TH1*> hPtGenAsym(param_->nPtBins());
  std::vector<TH1*> hPJet3(param_->nPtBins());
  std::vector<TH1*> hPSJ(param_->nPtBins());
  std::vector<TH1*> hPJet3Rel(param_->nPtBins());
  std::vector<TH1*> hPSJRel(param_->nPtBins());

  std::vector<double> pt3Limits;
  pt3Limits.push_back(0.06);
  pt3Limits.push_back(0.1);
  pt3Limits.push_back(0.15);
  std::vector< std::vector<TH1*> > hPtAsymPt3Cuts(param_->nPtBins());
  std::vector< std::vector<TH1*> > hPtAsymPp3Cuts(param_->nPtBins());

  for(unsigned int ptBin = 0; ptBin < param_->nPtBins(); ++ptBin) {
    TString binLabelGen = toString(param_->ptMin(ptBin))+" < p^{gen,ave}_{T} < "+toString(param_->ptMax(ptBin))+" GeV";

    TString name = "ParallelComponents:hPtAveGen"+toTString(ptBin);
    TH1 *h1 = new TH1D(name,binLabelGen+";p^{gen,ave}_{T} (GeV);Events",50,0.9*param_->ptMin(ptBin),1.1*param_->ptMax(ptBin));
    h1->SetMarkerStyle(20);
    hPtGenAve[ptBin] = h1;

    for(int j = 0; j < 3; ++j) {
      name = "ParallelComponents:hPtGenJet"+toTString(j)+"_"+toTString(ptBin);
      h1 = new TH1D(name,binLabelGen+";p^{gen}_{T,"+toString(1+j)+"} (GeV);Events",
		    50,j<2 ? 0.8*param_->ptMin(ptBin) : 0.,1.2*param_->ptMax(ptBin));
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
    h1 = new TH1D(name,binLabelGen+";p^{gen}_{||,3} (GeV);Events",50,0.,0.8*param_->ptMin(ptBin));
    h1->SetMarkerStyle(20);
    hPJet3[ptBin] = h1;

    name = "ParallelComponents:hPSJ"+toTString(ptBin);
    h1 = new TH1D(name,binLabelGen+";p^{gen}_{||,SJ} (GeV);Events",50,0.,0.8*param_->ptMin(ptBin));
    h1->SetMarkerStyle(20);
    hPSJ[ptBin] = h1;

    name = "ParallelComponents:hPJet3Rel"+toTString(ptBin);
    h1 = new TH1D(name,binLabelGen+";p^{gen}_{||,3} / <p^{gen,ave}_{T}>;Events",50,0.,0.8);
    h1->SetMarkerStyle(20);
    hPJet3Rel[ptBin] = h1;

    name = "ParallelComponents:hPSJRel"+toTString(ptBin);
    h1 = new TH1D(name,binLabelGen+";p^{gen}_{||,SJ} / <p^{gen,ave}_{T}>;Events",50,0.,0.8);
    h1->SetMarkerStyle(20);
    hPSJRel[ptBin] = h1;

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
    if( (*datait)->type() == DiJetResolution )  {
      DiJetResolutionEvent* dijet = static_cast<DiJetResolutionEvent*>(*datait);  
      double ptGenAve = dijet->avePtGen();
      unsigned int ptGenAveBin = 0;
      if( param_->findPtBin(ptGenAve,ptGenAveBin) ) {
 	hPtGenAve[ptGenAveBin]->Fill(ptGenAve,dijet->weight());
 	for(int i = 0; i < 2; ++i) {
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
 	double pt3 = 0.;////// Should be pt3Gen: dijet->jet3() ? dijet->jet3()->genPt() : 0.;
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
  for(unsigned int ptBin = 0; ptBin < param_->nPtBins(); ++ptBin) {
    for(size_t c = 0; c < pt3Limits.size(); ++c) {
      normHist(hPtAsymPt3Cuts[ptBin][c],"width");
      normHist(hPtAsymPp3Cuts[ptBin][c],"width");
    }
  }  

  std::vector<double> meanPtGenAve(param_->nPtBins(),1.);
  for(size_t ptBin = 0; ptBin < meanPtGenAve.size(); ++ptBin) {
    meanPtGenAve[ptBin] = hPtGenAve[ptBin]->GetMean();
  }
  // 2. loop over data
  for(DataIt datait = data_->begin(); datait != data_->end(); datait++) {
    // Select DiJet events
    if( (*datait)->type() == DiJetResolution )  {
      DiJetResolutionEvent* dijet = static_cast<DiJetResolutionEvent*>(*datait);  
      double ptGenAve = dijet->avePtGen();
      unsigned int ptGenAveBin = 0;
      if( param_->findPtBin(ptGenAve,ptGenAveBin) ) {
	hPJet3[ptGenAveBin]->Fill(dijet->pJ3(),dijet->weight());
 	hPSJ[ptGenAveBin]->Fill(dijet->pSJ(),dijet->weight());
	
  	hPJet3Rel[ptGenAveBin]->Fill(std::abs(dijet->pJ3())/meanPtGenAve[ptGenAveBin],dijet->weight());
  	hPSJRel[ptGenAveBin]->Fill(std::abs(dijet->pSJ())/meanPtGenAve[ptGenAveBin],dijet->weight());
//  	hPJet3Rel[ptGenAveBin]->Fill(std::abs(dijet->pJ3())/dijet->ptRef(),dijet->weight());
//  	hPSJRel[ptGenAveBin]->Fill(std::abs(dijet->pSJ())/dijet->ptRef(),dijet->weight());
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
  }

  TGraphErrors *gStdDevPtAsym
    = new TGraphErrors(meanPtGenAve.size(),&(meanPtGenAve.front()),&(stdDevPtAsym.front()),
		       &(meanPtGenAveErr.front()),&(stdDevPtAsymErr.front()));
  gStdDevPtAsym->SetMarkerStyle(20);
  gStdDevPtAsym->SetMarkerColor(color(0));
  gStdDevPtAsym->SetLineColor(color(0));
  gStdDevPtAsym->SetName("gStdDevPtAsym");

  TGraphErrors *gStdDevPtGenAsym
    = new TGraphErrors(meanPtGenAve.size(),&(meanPtGenAve.front()),&(stdDevPtGenAsym.front()),
		       &(meanPtGenAveErr.front()),&(stdDevPtGenAsymErr.front()));
  gStdDevPtGenAsym->SetMarkerStyle(24);
  gStdDevPtGenAsym->SetMarkerColor(color(4));
  gStdDevPtGenAsym->SetLineColor(color(4));
  gStdDevPtGenAsym->SetName("gStdDevPtGenAsym");

  TGraphErrors *gStdDevResp
    = new TGraphErrors(meanPtGenAve.size(),&(meanPtGenAve.front()),&(stdDevResp.front()),
		       &(meanPtGenAveErr.front()),&(stdDevRespErr.front()));
  gStdDevResp->SetMarkerStyle(21);
  gStdDevResp->SetMarkerColor(color(1));
  gStdDevResp->SetLineColor(color(1));
  gStdDevResp->SetName("gStdDevResp");

  TGraphErrors *gStdDevPJet3
    = new TGraphErrors(meanPtGenAve.size(),&(meanPtGenAve.front()),&(stdDevPJet3.front()),
		       &(meanPtGenAveErr.front()),&(stdDevPJet3Err.front()));
  gStdDevPJet3->SetMarkerStyle(22);
  gStdDevPJet3->SetMarkerColor(color(2));
  gStdDevPJet3->SetLineColor(color(2));
  gStdDevPJet3->SetName("gStdDevPJet3");

  TGraphErrors *gStdDevPSJ
    = new TGraphErrors(meanPtGenAve.size(),&(meanPtGenAve.front()),&(stdDevPSJ.front()),
		       &(meanPtGenAveErr.front()),&(stdDevPSJErr.front()));
  gStdDevPSJ->SetMarkerStyle(23);
  gStdDevPSJ->SetMarkerColor(color(3));
  gStdDevPSJ->SetLineColor(color(3));
  gStdDevPSJ->SetName("gStdDevPSJ");


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

  for(unsigned int ptBin = 0; ptBin < param_->nPtBins(); ptBin++) {
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
    hPJet3Rel[ptBin]->Draw("PE1");
    if( !saveAsEps_ ) c1->Draw();
    else c1->SaveAs(outName+"hPJet3Rel_PtBin"+toString(ptBin)+".eps","eps");

    if( !saveAsEps_ ) ps->NewPage();
    c1->cd();
    hPSJRel[ptBin]->Draw("PE1");
    if( !saveAsEps_ ) c1->Draw();
    else c1->SaveAs(outName+"hPSJRel_PtBin"+toString(ptBin)+".eps","eps");
  }

  if( !saveAsEps_ ) ps->NewPage();
  c1->cd();

  TF1 *fPtGenAsym = new TF1("fPtGenAsym","sqrt([0]*[0]/x[0]/x[0]+[1]*[1]/x[0]+[2]*[2])",
			    param_->ptMin(),param_->ptMax());
  fPtGenAsym->SetParameter(0,2.54877);
  fPtGenAsym->SetParameter(1,0.149045);
  fPtGenAsym->SetParameter(0,0.0109168);
  fPtGenAsym->SetLineWidth(1);
  fPtGenAsym->SetLineStyle(2);
  fPtGenAsym->SetLineColor(color(4));


  TH1 *hSumContributions = new TH1D("hSumContributions",";p^{gen}_{T} (GeV);Standard Deviation",
				    param_->nPtBins(),&(param_->ptBinEdges().front()));
  for(unsigned int i = 0; i < param_->nPtBins(); ++i) {
    double resp = gStdDevResp->GetY()[i];
    double pJ3 = gStdDevPJet3->GetY()[i];
    double pSJ = gStdDevPSJ->GetY()[i];
    double ptGenAsym = fPtGenAsym->Eval(meanPtGenAve[i]);
    hSumContributions->SetBinContent(1+i,sqrt( resp*resp + pJ3*pJ3 + pSJ*pSJ + ptGenAsym*ptGenAsym ));
  }
  hSumContributions->GetYaxis()->SetRangeUser(0,0.8);
  hSumContributions->GetXaxis()->SetMoreLogLabels();
  hSumContributions->Draw();

  TLegend *leg = createLegend(5,0.78);

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
  c1->SaveAs(outName+"hParallelContributions.C");


  if( !saveAsEps_ ) ps->NewPage();
  c1->cd();
  TH1 *hFramePtAsym = new TH1D("hFramePtAsym",";<p^{gen,ave}> (GeV);Resolution",
			       1000,param_->ptMin(),param_->ptMax());
  hFramePtAsym->GetYaxis()->SetRangeUser(0,0.6);
  hFramePtAsym->Draw();
  truthRes_->Draw("same");
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

  TFile rootfile((dir_+"/jsResponse.root").c_str(),"UPDATE");
  rootfile.WriteTObject(hSumContributions);
  rootfile.WriteTObject(fPtGenAsym);
  rootfile.WriteTObject(gStdDevResp);
  rootfile.WriteTObject(gStdDevPtAsym);
  rootfile.WriteTObject(gStdDevPJet3);
  rootfile.WriteTObject(gStdDevPSJ);
  rootfile.Close();


  // Delete
  std::cout << "  Cleaning up" << std::endl;
  for(size_t i = 0; i < hPPhi.size(); ++i) {
    delete hPtGenAve[i];
    delete hPPhi[i];
    delete hPtAsym[i];
    delete hPJet3[i];
    delete hPSJ[i];
    delete hPJet3Rel[i];
    delete hPSJRel[i];
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
  delete legPtAsym;
  delete leg;
  delete gStdDevPtAsym;
  delete gStdDevResp;
  delete gStdDevPJet3;
  delete gStdDevPSJ;
  delete c1;
  if( ps ) delete ps;
}



// --------------------------------------------------
void ControlPlotsResolution::plotParameters() const {
  std::cout << "Creating parameter control plots\n";

  // ----- Quantities defining the parametrization -----
  int nPar = param_->numberOfParameters();


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

    hPars->SetBinContent(bin,param_->parameters()[i]);
    hPars->SetBinError(bin,param_->errors()[i]);

    hAbsPars->SetBinContent(bin,scale(i)*param_->parameters()[i]);
    hAbsPars->SetBinError(bin,scale(i)*param_->errors()[i]);
  
    hRelParErrors->SetBinContent(bin,param_->errors()[i]/param_->parameters()[i]);

    hGlobalCorr->SetBinContent(bin,param_->globalCorrCoeff()[i]);
  }
  for(int i = 0; i < nPar; i++) {
    for(int j = 0; j < i+1; j++) {
      int idx = (i*i + i)/2 + j;
      double corr = 0.;
      if( param_->errors()[i] && param_->errors()[j] ) {
	corr = param_->covCoeff()[idx] / param_->errors()[i] / param_->errors()[j];
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
void ControlPlotsResolution::plotParameterScan() const {
  std::cout << "Creating parameter scan control plots\n";

  // ----- Set up quantities -----
  int n = 4;
  int nSteps = 2*n+1;
  double nSigma = 2;
  int nPar = param_->numberOfParameters();
  std::vector<TH2D*> hParScans2D;
  std::vector<TLine*> lines;

  // Store likelihood for original parameter values
  double offset = 0.;
  for(DataIt dataIt = data_->begin(); dataIt != data_->end(); dataIt++) {
    if( (*dataIt)->type() == DiJetResolution )  { // Select DiJet events
      offset += (*dataIt)->chi2();
    }
  }     

  // ----- Vary parameters and calculate likelihood -----
  // Store original parameter values
  std::vector<double> origPars(nPar);
  for(int i = 0; i < nPar; i++) {
    origPars[i] = param_->parameters()[i];
  }
  // Outer loop over parameters
  for(int i = 0; i < nPar; i++) {
    if( param_->isFixedPar(i) ) continue;
    double idVal = nSigma*param_->errors()[i];
    if( idVal == 0 ) idVal = 0.1;
			
    // Inner loop over parameters
    for(int j = 0; j < i; j++) {
      if( param_->isFixedPar(j) ) continue;
      double jdVal = nSigma*param_->errors()[j];
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
	param_->parameters()[i] = iPar;
	for(int js = 0; js < nSteps; js++) {
	  double jPar = origPars[j] + (js-n)*jdVal; 
	  param_->parameters()[j] = jPar;
	  double deltaLkh = 0.;
	  // Calculate likelihood for varied parameters
	  if( is != n || js != n ) {
	    for(DataIt dataIt = data_->begin(); dataIt != data_->end(); dataIt++) {
	      if( (*dataIt)->type() == DiJetResolution )  { // Select DiJet events
		deltaLkh += (*dataIt)->chi2();
	      }
	    }    
	    deltaLkh -= offset;
	  }
	  hParScans2D.back()->Fill(iPar,jPar,deltaLkh);
	}
      }
      // Reset parameters to original values
      param_->parameters()[i] = origPars[i];
      param_->parameters()[j] = origPars[j];
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
void ControlPlotsResolution::plotLogP() const {
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
    if( (*dataIt)->type() == DiJetResolution )  {
      logPWend.push_back((*dataIt)->chi2());
      logPend.push_back((*dataIt)->chi2()/(*dataIt)->weight());
    }
  }

  // Store fitted parameters
  std::vector<double> fittedPar(param_->numberOfParameters());
  for(int i = 0; i < param_->numberOfParameters(); i++) {
    fittedPar[i] = param_->parameters()[i];
  }
  // Copy start values into parameter array
  std::vector<double> startParGlobal = bag_of<double>(config_->read<string>("global jet start values",""));
  for(int i = 0; i < param_->numberOfParameters(); i++) {
    if( i < param_->numberOfJetParameters() )
      param_->parameters()[i] = param_->jetStartPar(i);
    else
      param_->parameters()[i] = startParGlobal[i-param_->numberOfJetParameters()];
  }
  // Loop over data and fill -log(P) with start
  // parameter values
  for(DataIt dataIt = data_->begin(); dataIt != data_->end(); dataIt++) {
    // Select DiJet events
    if( (*dataIt)->type() == DiJetResolution )  {
      logPWstart.push_back((*dataIt)->chi2());
      logPstart.push_back((*dataIt)->chi2()/(*dataIt)->weight());
    }
  }
  // Copy back fitted values into parameter array
  for(int i = 0; i < param_->numberOfParameters(); i++) {
    param_->parameters()[i] = fittedPar.at(i);
  }


  // ----- Fill histograms of -log(P) -----
  std::sort(logPstart.begin(),logPstart.end());
  std::sort(logPend.begin(),logPend.end());
  double max = logPstart.back() > logPend.back() ? logPstart.back() : logPend.back();
  TH1* hLogPstart = new TH1D("hLogPstart",";- ln(P)",100,0.,1.1*max);
  hLogPstart->SetLineWidth(2);
  hLogPstart->SetLineColor(4);
  TH1* hLogPend = static_cast<TH1*>(hLogPstart->Clone("hLogPend"));
  hLogPend->SetLineColor(2);

  std::sort(logPWstart.begin(),logPWstart.end());
  std::sort(logPWend.begin(),logPWend.end());
  max = logPWstart.back() > logPWend.back() ? logPWstart.back() : logPWend.back();
  TH1* hLogPWstart = new TH1D("hLogPWstart",";- w #upoint ln(P)",100,0.,1.1*max);
  hLogPWstart->SetLineWidth(2);
  hLogPWstart->SetLineColor(4);
  TH1* hLogPWend = static_cast<TH1*>(hLogPWstart->Clone("hLogPWend"));
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
void ControlPlotsResolution::plotMeanResponseAndResolution() const {
  std::cout << "Creating response and resolution fits\n";

  // ----- Create histograms -----

  // Response vs ptgen
  // Default binning for MC truth plots
  std::vector<double> bins(35+1);
  equidistLogBins(bins,35,4.,2700.);

  TH2* hRespVsPtGen = new TH2D("MeanResp_hRespVsPtGen",
			       ";p^{gen}_{T} (GeV);p^{rec}_{T} / p^{gen}_{T}",
			       bins.size()-1,&(bins.front()),201,0.,2.);
  hRespVsPtGen->SetNdivisions(505);
  hRespVsPtGen->Sumw2();

  TH2* hRespVsPtGen_FineBins = new TH2D("MeanResp_hRespVsPtGen_FineBins",
					";p^{gen}_{T} (GeV);p^{rec}_{T} / p^{gen}_{T}",
					500,0.,1000.,501,0.,50.);
  hRespVsPtGen_FineBins->SetNdivisions(505);
  hRespVsPtGen_FineBins->Sumw2();

  TH2* hRespVsPtRec_FineBins = static_cast<TH2*>(hRespVsPtGen_FineBins->Clone("hRespVsPtRec_FineBins"));
  hRespVsPtRec_FineBins->GetXaxis()->SetTitle("p^{rec}_{T} (GeV)");

  TH2* hWeightVsResp = new TH2D("MeanResp_hWeightVsResp",
				";p^{rec}_{T} / p^{gen}_{T};Event Weight",
				1000,0.,100.,1000,0.,10000.);
  hWeightVsResp->SetNdivisions(505);
  hWeightVsResp->Sumw2();


  // ----- Fill histograms -----
  for(DataIt datait = data_->begin(); datait != data_->end(); datait++) {
    // Select DiJet events
    if( (*datait)->type() == DiJetResolution )  {
      DiJetResolutionEvent* dijet = static_cast<DiJetResolutionEvent*>(*datait);

      for(int i = 0; i < 2; ++i) {
	double ptrec = dijet->jet(i)->pt();
	double ptgen = dijet->jet(i)->genPt();
	if( ptgen > 0. ) {
	  hRespVsPtGen->Fill(ptgen,ptrec/ptgen);//,dijet->weight());
	  hRespVsPtGen_FineBins->Fill(ptgen,ptrec/ptgen);//,dijet->weight());
	  hRespVsPtRec_FineBins->Fill(ptrec,ptrec/ptgen);//,dijet->weight());
	  hWeightVsResp->Fill(ptrec/ptgen,dijet->weight());
	}
      } 
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
		     param_->ptMin(),param_->ptMax());
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
  hRespVsPtGen_FineBins->Draw("COLZ");
  c1->SetLogz();
  c1->Draw();

  ps->NewPage();
  c1->cd();
  hRespVsPtRec_FineBins->Draw("COLZ");
  c1->SetLogz();
  c1->Draw();

  ps->NewPage();
  c1->cd();
  hWeightVsResp->Draw("COLZ");
  c1->SetLogz();
  c1->Draw();

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
  rootfile.WriteTObject(hRespVsPtGen);
  rootfile.WriteTObject(hRespVsPtGen_FineBins);
  rootfile.WriteTObject(hRespVsPtRec_FineBins);
  rootfile.WriteTObject(hWeightVsResp);


  // ----- Clean up -----
  for(size_t i = 0; i < 2; ++i) {
    delete hResp[i];
    delete hReso[i];
    delete fReso[i];
    delete fitLabel[i];
  }
  delete hRespVsPtGen;
  delete hRespVsPtGen_FineBins;
  delete hRespVsPtRec_FineBins;
  delete hWeightVsResp;
  delete leg;
  delete c1;
  delete ps;
}



// --------------------------------------------------
double ControlPlotsResolution::spectrum(double *x, double *par) {
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
double ControlPlotsResolution::gaussian(double *x, double *par) {
  double u = (x[0]-1.)/par[1];
  return par[0]*exp(-0.5*u*u);
}



// --------------------------------------------------
double ControlPlotsResolution::gaussianWidth(double pt) const {
  double a[3];
  for(int i = 0; i < 3; i++) {
    a[i] = scale(i)*param_->parameters()[i];
  }
  return sqrt( a[0]*a[0] + a[1]*a[1]*pt + a[2]*a[2]*pt*pt );
}


// ------------------------------------------------------------------------
double ControlPlotsResolution::gaussianWidthError(double pt) const {
  // Calculate derivatives
  std::vector<double> ds(3);
  double s = gaussianWidth(pt);
  for(int i = 0; i < 3; i++) {
    ds[i] = scale(i)*param_->parameters()[i]/s;
    if( i == 1 ) ds[i] *= pt;
    if( i == 2 ) ds[i] *= pt*pt;
  }

  // Calculate variance
  double var = 0.;
  for(int i = 0; i < 3; i++) { // Outer loop over parameters
    for(int j = 0; j < i+1; j++) { // Inner loop over parameters
      int idx = (i*i + i)/2 + j; // Index of (i,j) in covariance vector
      if( i == j ) { // Diagonal terms
	var += ds[i]*ds[i]*scale(i)*scale(i)*param_->covCoeff()[idx];
      } else { // Off-diagonal terms
	var += 2*ds[i]*ds[j]*scale(i)*scale(j)*param_->covCoeff()[idx];
      }
    } // End of inner loop over parameters
  } // End of outer loop over parameters

  return sqrt(var);
}


// --------------------------------------------------
double ControlPlotsResolution::gaussianWidthTruth(double pt) const {
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
void ControlPlotsResolution::drawPSPage(TPostScript * ps, TCanvas * can, TObject * obj, const std::string &option, bool log) const {
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
void ControlPlotsResolution::drawPSPage(TPostScript * ps, TCanvas * can, std::vector<TObject*> objs, const std::string &option, bool log) const {
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
void ControlPlotsResolution::findYRange(const TH1 * h, double& min, double& max) const {
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
void ControlPlotsResolution::setGStyle() const
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
  gStyle->SetPadTopMargin(0.11);
  gStyle->SetPadBottomMargin(0.16);
  gStyle->SetPadLeftMargin(0.21);
  gStyle->SetPadRightMargin(0.06);

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
  //  gStyle->SetTitleFontSize(0.1);
  gStyle->SetTitleAlign(13);
  gStyle->SetTitleX(0.18);
//   gStyle->SetTitleY(1.);
//   gStyle->SetTitleH(0.15);
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
}



// --------------------------------------------------
TLegend *ControlPlotsResolution::createLegend(int nEntries, double width, double lineHgt, double yOffset) const {
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
TPaveText *ControlPlotsResolution::createPaveText(int nEntries, double width, double lineHgt) const {
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
void ControlPlotsResolution::setYRange(TH1 * h, double c1, double c2, double minLimit) const {
  double min = 0.;
  double max = 0.;
  findYRange(h,min,max);
  min *= c1;
  max *= c2;
  if( min < minLimit ) min = minLimit;
  h->GetYaxis()->SetRangeUser( min, max );
}

void ControlPlotsResolution::normHist(TH1 *h, std::string option) const { 
  if( h->Integral(option.c_str()) ) h->Scale(1./h->Integral(option.c_str())); 
}

void ControlPlotsResolution::normHist(TH1 *h, double min, double max, std::string option) const { 
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
template <class T> std::string ControlPlotsResolution::toString(const T& t) const {
  std::stringstream ss;
  ss << t;
  return ss.str();
}
   

//!  Filling \p bins with borders of \p nBins bins between \p first
//!  and \p last that are equidistant when viewed in log scale,
//!  so \p bins must have length \p nBins+1. If \p first, \p last
//!  or \p nBins are not positive, failure is reported.
// -------------------------------------------------------------
bool ControlPlotsResolution::equidistLogBins(std::vector<double>& bins, int 
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


int ControlPlotsResolution::color(int i) const {
  int               col = 1;
  if( i == 1 )      col = 2;
  else if( i == 2 ) col = 4;
  else if( i == 3 ) col = 6;
  else if( i == 4 ) col = 8;
  else if( i == 5 ) col = 9;

  return col;
}


// -------------------------------------------------------------------------------------
void ControlPlotsResolution::setAxisTitles(TH1 *h, const std::string &xTitle, const std::string &xUnit, const std::string &yTitle, bool norm) const {
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
void ControlPlotsResolution::setColor(TH1 *h, int color) const {
  h->SetMarkerColor(color);
  h->SetLineColor(color);
}
