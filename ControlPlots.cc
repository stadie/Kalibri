#include "ControlPlots.h"

#include <iostream>
#include <string>
#include <vector>

#include "TError.h"
#include "TStyle.h"

#include "CalibData.h"
#include "ConfigFile.h"
#include "ControlPlotsConfig.h"
#include "ControlPlotsFunction.h"
#include "ControlPlotsProfile.h"
#include "Function.h"
#include "Jet.h"
#include "JetTruthEvent.h"
#include "TwoJetsPtBalanceEvent.h"

//!  \brief Constructor
//! 
//!  \param configFile  Configuration file
//!  \param data        Plotted data
// -------------------------------------------------------------
ControlPlots::ControlPlots(const ConfigFile *configFile, const std::vector<Event*> *data)
  : config_(configFile), data_(data) {
  gErrorIgnoreLevel = 1001;
  setGStyle();
}



//!  \brief Create all plots as specified in the config file
// -------------------------------------------------------------
void ControlPlots::makePlots() const {
  if( config_->read<bool>("create JetTruthEvent plots",false) )
    createJetTruthEventPlots();
  if( config_->read<bool>("create TwoJetsPtBalanceEvent plots",false) )
    createTwoJetsPtBalanceEventPlots();
  
}



//!  \brief Control plots for \p JetTruthEvent data
// -------------------------------------------------------------
void ControlPlots::createJetTruthEventPlots() const {
  std::cout << "Creating JetTruthEvent plots\n";
    
  // Read different control plot names
  std::vector<std::string> names = bag_of_string(config_->read<std::string>("JetTruthEvent plots names",""));
  // Loop over names
  std::vector<ControlPlotsConfig*> configs(names.size());
  std::vector<ControlPlotsFunction*> functions(names.size());
  std::vector<ControlPlotsProfile*> profiles(names.size());
  for(size_t i = 0; i < names.size(); i++) {
    std::cout << " Creating plots '" << names.at(i) << "'\n";
  
    // Create ControlPlotsConfig    
    ControlPlotsConfig *pConfig = new ControlPlotsConfig(config_,names.at(i));
    configs.at(i) = pConfig;

    // Create functions
    ControlPlotsFunction *func = new ControlPlotsFunction();
    func->setBinFunction(findJetTruthEventFunction(pConfig->binVariable()));
    func->setXFunction(findJetTruthEventFunction(pConfig->xVariable()));
    func->setCutFunction(findJetTruthEventFunction(pConfig->cutVariable()));
    ControlPlotsConfig::CorrectionTypeIt corrTypeIt = pConfig->correctionTypesBegin();
    for(; corrTypeIt != pConfig->correctionTypesEnd(); corrTypeIt++) {
      func->addYFunction(*corrTypeIt,findJetTruthEventFunction(pConfig->yVariable(),*corrTypeIt));
    }
    functions.at(i) = func;

    // Create profile
    profiles.at(i) = new ControlPlotsProfile(pConfig,func);
  } // End of loop over names

  
  // Fill histograms
  std::cout << "  Filling plots\n";	  
  for( DataIt evt = data_->begin(); evt != data_->end(); evt++ ) {
    JetTruthEvent *jte = dynamic_cast<JetTruthEvent*>(*evt);
    if( jte ) {
      for(size_t i = 0; i < configs.size(); i++) {
	profiles.at(i)->fill(jte);
      }
    }
  }

  // Fitting profiles and writing plots to file
  std::cout << "  Fitting profiles and writing plots to file" << std::endl;
	  
  for(size_t i = 0; i < configs.size(); i++) {
    profiles.at(i)->fitProfiles();
    profiles.at(i)->draw();
  }

  // Cleaning up
  for(size_t i = 0; i < configs.size(); i++) {
    delete configs.at(i);
    delete functions.at(i);
    delete profiles.at(i);
  }
}

//!  \brief Control plots for \p TwoJetPtBalanceEvent data
// -------------------------------------------------------------
void ControlPlots::createTwoJetsPtBalanceEventPlots() const {
  std::cout << "Creating TwoJetsPtBalanceEvent plots\n";
    
  // Read different control plot names
  std::vector<std::string> names = bag_of_string(config_->read<std::string>("TwoJetsPtBalanceEvent plots names",""));
  // Loop over names
  std::vector<ControlPlotsConfig*> configs(names.size());
  std::vector<ControlPlotsFunction*> functions(names.size());
  std::vector<ControlPlotsProfile*> profiles(names.size());
  for(size_t i = 0; i < names.size(); i++) {
    std::cout << " Creating plots '" << names.at(i) << "'\n";
  
    // Create ControlPlotsConfig    
    ControlPlotsConfig *pConfig = new ControlPlotsConfig(config_,names.at(i));
    configs.at(i) = pConfig;

    // Create functions
    ControlPlotsFunction *func = new ControlPlotsFunction();
    func->setBinFunction(findTwoJetsPtBalanceEventFunction(pConfig->binVariable()));
    func->setXFunction(findTwoJetsPtBalanceEventFunction(pConfig->xVariable()));
    func->setCutFunction(findTwoJetsPtBalanceEventFunction(pConfig->cutVariable()));
    ControlPlotsConfig::CorrectionTypeIt corrTypeIt = pConfig->correctionTypesBegin();
    for(; corrTypeIt != pConfig->correctionTypesEnd(); corrTypeIt++) {
      func->addYFunction(*corrTypeIt,findTwoJetsPtBalanceEventFunction(pConfig->yVariable(),*corrTypeIt));
    }
    functions.at(i) = func;

    // Create profile
    profiles.at(i) = new ControlPlotsProfile(pConfig,func);
  } // End of loop over names

  
  // Fill histograms
  std::cout << "  Filling plots\n";	  
  for( DataIt evt = data_->begin(); evt != data_->end(); evt++ ) {
    TwoJetsPtBalanceEvent *jte = dynamic_cast<TwoJetsPtBalanceEvent*>(*evt);
    if( jte ) {
      for(size_t i = 0; i < configs.size(); i++) {
	profiles.at(i)->fill(jte);
      }
    }
  }

  // Fitting profiles and writing plots to file
  std::cout << "  Fitting profiles and writing plots to file" << std::endl;
	  
  for(size_t i = 0; i < configs.size(); i++) {
    profiles.at(i)->fitProfiles();
    profiles.at(i)->draw();
  }

  // Cleaning up
  for(size_t i = 0; i < configs.size(); i++) {
    delete configs.at(i);
    delete functions.at(i);
    delete profiles.at(i);
  }
}



//!  \brief Helper method for \p createJetTruthEventPlots()
//!
//!  Returns the \p ControlPlotsFunction::Function for a variable
//!  with name \p varName and a \p ControlPlotsConfig::CorrectionType
//!  \p type.
// -------------------------------------------------------------
ControlPlotsFunction::Function ControlPlots::findJetTruthEventFunction(const std::string& varName, ControlPlotsConfig::CorrectionType type) const {
  if( varName == "Eta" )
    return  &ControlPlotsFunction::jetTruthEventJetEta;
  if( varName == "AbsEta" )
    return  &ControlPlotsFunction::jetTruthEventJetAbsEta;
  if( varName == "Pt" )
   return  &ControlPlotsFunction::jetTruthEventJetPt;
  if( varName == "EMF" )
    return  &ControlPlotsFunction::jetTruthEventJetEMF;
  if( varName == "momentEtaEta" )
    return  &ControlPlotsFunction::jetTruthEventJetMomentEtaEta;
  if( varName == "scaledEtaEta" )
    return  &ControlPlotsFunction::jetTruthEventJetScaledEtaEta;
  if( varName == "scaledPhiPhi" )
    return  &ControlPlotsFunction::jetTruthEventJetScaledPhiPhi;
  if( varName == "GenJetPt" )
    return  &ControlPlotsFunction::jetTruthEventTruthPt;
  if( varName == "momentPhiPhi" )
    return &ControlPlotsFunction::jetTruthEventJetMomentPhiPhi; 
  if( varName == "meanMoment" )
    return &ControlPlotsFunction::jetTruthEventJetMeanMoment;
  if( varName == "Flavor" )
    return &ControlPlotsFunction::jetTruthEventJetFlavor;
  if( varName == "GenJetResponse" && type == ControlPlotsConfig::Uncorrected )
   return &ControlPlotsFunction::jetTruthEventResponse;
  if( varName == "GenJetResponse" && type == ControlPlotsConfig::Kalibri )
    return &ControlPlotsFunction::jetTruthEventResponseKalibriCorrected;
  if( varName == "GenJetResponse" && type == ControlPlotsConfig::L2L3 )
    return  &ControlPlotsFunction::jetTruthEventResponseL2L3Corrected;
  if( varName == "") {
    return 0;
  }
  std::cerr << "ControlPlots: unknown variable " << varName << std::endl;
  return 0;
}


//!  \brief Helper method for \p createTwoJetsPtBalanceEventPlots()
//!
//!  Returns the \p ControlPlotsFunction::Function for a variable
//!  with name \p varName and a \p ControlPlotsConfig::CorrectionType
//!  \p type.
// -------------------------------------------------------------
ControlPlotsFunction::Function ControlPlots::findTwoJetsPtBalanceEventFunction(const std::string& varName, ControlPlotsConfig::CorrectionType type) const {
  if( varName == "Eta" )
    return  &ControlPlotsFunction::twoJetsPtBalanceEventJetEta;
  if( varName == "AbsEta" )
    return  &ControlPlotsFunction::twoJetsPtBalanceEventJetAbsEta;
  if( varName == "Pt" )
   return  &ControlPlotsFunction::twoJetsPtBalanceEventJetPt; 
  if( varName == "MeanPt" )
    return  &ControlPlotsFunction::twoJetsPtBalanceEventMeanPt;
  if( varName == "EMF" )
    return  &ControlPlotsFunction::twoJetsPtBalanceEventJetEMF;
  if( varName == "momentEtaEta" )
    return  &ControlPlotsFunction::twoJetsPtBalanceEventJetMomentEtaEta;
  if( varName == "momentPhiPhi" )
    return &ControlPlotsFunction::twoJetsPtBalanceEventJetMomentPhiPhi; 
  if( varName == "meanMoment" )
    return &ControlPlotsFunction::twoJetsPtBalanceEventJetMeanMoment;
  if( varName == "Asymmetry" && type == ControlPlotsConfig::Uncorrected )
    return &ControlPlotsFunction::twoJetsPtBalanceEventAsymmetry;
  if( varName == "Asymmetry" && type == ControlPlotsConfig::Kalibri )
    return &ControlPlotsFunction::twoJetsPtBalanceEventAsymmetryKalibriCorrected;
  if( varName == "Asymmetry" && type == ControlPlotsConfig::L2L3 )
    return  &ControlPlotsFunction::twoJetsPtBalanceEventAsymmetryL2L3Corrected;
  if( varName == "B" && type == ControlPlotsConfig::Uncorrected )
    return &ControlPlotsFunction::twoJetsPtBalanceEventB;
  if( varName == "B" && type == ControlPlotsConfig::Kalibri )
    return &ControlPlotsFunction::twoJetsPtBalanceEventBKalibriCorrected;
  if( varName == "B" && type == ControlPlotsConfig::L2L3 )
    return  &ControlPlotsFunction::twoJetsPtBalanceEventBL2L3Corrected;
  if( varName == "") {
    return 0;
  }
  
  std::cerr << "ControlPlots: unknown variable " << varName << std::endl;
  return 0;
}



//!  Set style option for the output.
//---------------------------------------------------------------
void ControlPlots::setGStyle() const {
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
  gStyle->SetPadLeftMargin(0.19);
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

