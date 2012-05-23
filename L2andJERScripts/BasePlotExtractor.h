#ifndef BasePlotExtractor_h
#define BasePlotExtractor_h

#include <vector>
#include <iostream>
#include <algorithm>
//#include "../progressbar.h"
#include "../ConfigFile.cc"
#include "../ControlPlotsConfig.cc"
#include "../ControlPlotsProfile.cc"
#include "../ControlPlotsFunction.cc"
#include "DefaultStyles.cc"
//#include "../ControlPlots.cc"
//#include "../Bin.cc"
#include "TString.h"
#include "TH1.h"
#include "TF1.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TList.h"

#include "TFitResult.h"
#include "THelpers_2.h"
#include "MakeDateDir.h"
//#include "TString.h"
#define UTILS_AS_HEADER_FILE
#include "util/utils.h"
#include "util/HistOps.h"
#include "util/LabelFactory.h"
#include "util/StyleSettings.h"


typedef std::vector<TH1D*> TH1vec_t;
typedef std::vector<TH1vec_t > VecOfTH1vec_t;
typedef std::vector<VecOfTH1vec_t > VecOfVecOfTH1vec_t;

//!  \brief Reads (PF-)fraction plots and produces stacked
//!         histograms of the input
//!
//!  Makes use of Kalibri classes used for ControlPlots to get binning
//!  names and other stuff right.
//!  Default is to 
//!
//!
//!  
//!
//!  \author Henning Kirschenmann
//!  \date 2012/03/07
// ----------------------------------------------------------------   
class BasePlotExtractor {
 public :
  BasePlotExtractor(TString plotsnames="AbsPFFractionVsPt",TString kalibriPlotsPath="KalibriPlots.root");
  void init(TString profileType="Mean");
  void makeRatioVsBinVarHistos();
  void refreshRatiosDataMC();
  void drawConfidenceIntervals(TH1D* histo);
  void addFunctionLabelsToLegend(TH1D* histo, TLegend* leg);
  void outputTable(TString label, TH1D* histo);
  void readInExtraInfo();
  TString kalibriPlotsPath() {return kalibriPlotsPath_;};
  TString profileType() {return profileType_;};
  std::vector<double> yProfileMinMax() {return yProfileMinMax_;};
  TString yProfileTitle() {return yProfileTitle_;};
  std::vector<double> yRatioMinMax() {return yRatioMinMax_;};
  TString yRatioTitle() {return yRatioTitle_;};
  std::vector<double> yDifferenceMinMax() {return yDifferenceMinMax_;};
  TString yDifferenceTitle() {return yDifferenceTitle_;};

  //  void Plot();
  // private:
  ConfigFile* getConfig() {return &config_;}
  ConfigFile config_;
  ConfigFile ExternalConfig_;
  std::vector <TH1vec_t> FractionPlots_; 
  std::vector<std::vector<Event*>* > dummy_;       //!< The data
  ControlPlotsFunction::Function findTwoJetsPtBalanceEventFunction(const std::string& varName, ControlPlotsConfig::CorrectionType type = ControlPlotsConfig::Uncorrected) const;
  std::vector<std::string> names_;
  std::vector<ControlPlotsConfig*> configs_;
  std::vector<ControlPlotsFunction*> functions_;
  std::vector<ControlPlotsProfile*> profiles_;
  VecOfVecOfTH1vec_t AllPlots_;
  VecOfTH1vec_t AllRatiosDataMC_;
  VecOfTH1vec_t AllDifferencesDataMC_;
  TH1vec_t RatioVsBinVarHistos_;
  TString kalibriPlotsPath_;
  TString plotsnames_;
  TString jetType_;
  TString jetLabel_;
  TString profileType_;
  TString binningSelection_;
  std::vector<double> yProfileMinMax_;
  TString yProfileTitle_;
  std::vector<double> yRatioMinMax_;
  TString yRatioTitle_;
  std::vector<double> yDifferenceMinMax_;
  TString yDifferenceTitle_;
};

#endif 


