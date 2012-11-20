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
#include "progressbar.h"


//!  \brief Constructor
//! 
//!  \param configFile  Configuration file
//!  \param data        Plotted data
// -------------------------------------------------------------
ControlPlots::ControlPlots(const ConfigFile *configFile, const std::vector<std::vector<Event*>* >& samples, const EventProcessor* eventProcessor)
  : config_(configFile), samples_(samples), eventProcessor_(eventProcessor){
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
    
    if(eventProcessor_!=0) {
      pConfig->determineOutPlotSuffix(eventProcessor_->name());
    }
    else pConfig->setOutPlotSuffix("");
    //    pConfig->setOutPlotSuffix(root_filename_suffix);
    pConfig->setOutRootFileName("KalibriPlots"+pConfig->outPlotSuffix()+".root");
    configs.at(i) = pConfig;

    // Create functions
    ControlPlotsFunction *func = new ControlPlotsFunction();
    func->setBinFunction(findJetTruthEventFunction(pConfig->binVariable()));
    func->setXFunction(findJetTruthEventFunction(pConfig->xVariable()));
    func->setCutFunction(findJetTruthEventFunction(pConfig->cutVariable()));
    for(ControlPlotsConfig::InputTagsIterator it = pConfig->inputTagsBegin() ; 
	it != pConfig->inputTagsEnd(); ++it) {
      func->addYFunction(it->second,findJetTruthEventFunction(pConfig->yVariable(),it->second));
    }
    functions.at(i) = func;

    // Create profile
    profiles.at(i) = new ControlPlotsProfile(pConfig,func);
  } // End of loop over names

  
  // Fill histograms
  std::cout << "  Filling plots\n";
  for(unsigned int id = 0 ; id < samples_.size() ; ++id) {
    int l = samples_[id]->size()/100;
    for(DataIt evt = samples_[id]->begin() ; evt != samples_[id]->end() ; ++evt)  {
      if((evt-samples_[id]->begin() +1) % 10000 == 0)
	progressbar((evt-samples_[id]->begin()+1)/l);
      JetTruthEvent *jte = dynamic_cast<JetTruthEvent*>(*evt); 
      if( jte ) {
        for(size_t i = 0; i < configs.size(); i++) {
	  profiles.at(i)->fill(jte,id);
	}
      }
    }
  std::cout << '\n';
  }  

  // Fitting profiles and writing plots to file
  std::cout << "  Fitting profiles and writing plots to file" << std::endl;
	  
  for(size_t i = 0, l = configs.size(); i < l; i++) {
    // Create / open ROOT file for output
    configs.at(i)->openRootFile();
    profiles.at(i)->fitProfiles();
    progressbar((i+1)*100/l);
    profiles.at(i)->draw();
    configs.at(i)->closeRootFile();

  }
  std::cout << '\n';
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
  std::vector<std::string> names;
  std::string root_filename_suffix="";  

  //  
  if(eventProcessor_!=0) {

    std::cout << "plots done before event processor execution... "<< eventProcessor_->name()<< std::endl;
    //    std::cout << "JFGIJHFGIJHFGI COTEIOJTJIRG CONFIG " << eventProcessor_->name()<< std::endl;
    names = bag_of_string(config_->read<std::string>("TwoJetsPtBalanceEvent " +eventProcessor_->name()+" plots names",""));

  }
  else
    {
  // Read different control plot names
  names = bag_of_string(config_->read<std::string>("TwoJetsPtBalanceEvent plots names",""));
    }
  // Loop over names
  std::vector<ControlPlotsConfig*> configs(names.size());
  std::vector<ControlPlotsFunction*> functions(names.size());
  std::vector<ControlPlotsProfile*> profiles(names.size());
  for(size_t i = 0; i < names.size(); i++) {
    std::cout << " Creating plots '" << names.at(i) << "'\n";
  
    // Create ControlPlotsConfig    
    ControlPlotsConfig *pConfig = new ControlPlotsConfig(config_,names.at(i));
    configs.at(i) = pConfig;
    if(eventProcessor_!=0) {
      pConfig->determineOutPlotSuffix(eventProcessor_->name());
    }
    else pConfig->setOutPlotSuffix("");
    //    pConfig->setOutPlotSuffix(root_filename_suffix);
    pConfig->setOutRootFileName("KalibriPlots"+pConfig->outPlotSuffix()+".root");
    //    std::cout << pConfig->outRootFileName()  << std::endl;
    // Create functions
    ControlPlotsFunction *func = new ControlPlotsFunction();
    func->setBinFunction(findTwoJetsPtBalanceEventFunction(pConfig->binVariable()));
    func->setXFunction(findTwoJetsPtBalanceEventFunction(pConfig->xVariable()));
    func->setCutFunction(findTwoJetsPtBalanceEventFunction(pConfig->cutVariable()));
    for(ControlPlotsConfig::InputTagsIterator it = pConfig->inputTagsBegin() ; 
	it != pConfig->inputTagsEnd(); ++it) {
      func->addYFunction(it->second,findTwoJetsPtBalanceEventFunction(pConfig->yVariable(),it->second));
    }
    functions.at(i) = func;

    // Create profile
    profiles.at(i) = new ControlPlotsProfile(pConfig,func);
  } // End of loop over names

  
  // Fill histograms
  std::cout << "  Filling plots\n"; 
  for(unsigned int id = 0 ; id < samples_.size() ; ++id) {
    int l = samples_[id]->size()/100;
    for( DataIt evt = samples_[id]->begin(); evt != samples_[id]->end(); ++evt ) { 
      if((evt-samples_[id]->begin() + 1) % 10000 == 0) 
	progressbar((evt-samples_[id]->begin())/l);
      TwoJetsPtBalanceEvent *jte = dynamic_cast<TwoJetsPtBalanceEvent*>(*evt);
      if( jte ) {
	for(size_t i = 0; i < configs.size(); i++) {
	  profiles.at(i)->fill(jte,id);
	}
      }
    }
    std::cout << '\n';
  }

  // Fitting profiles and writing plots to file
  std::cout << "  Fitting profiles and writing plots to file" << std::endl;
	  
  for(size_t i = 0, l = configs.size(); i < l; i++) {
    // Create / open ROOT file for output
    configs.at(i)->openRootFile();
    profiles.at(i)->fitProfiles();
    progressbar((i+1)*100/l);
    profiles.at(i)->draw();
    configs.at(i)->closeRootFile();

  }
  std::cout << '\n';
  // Cleaning up
  for(size_t i = 0; i < configs.size(); i++) {
    //std::cout << 'works\n';
    delete configs.at(i);
    //std::cout << 'worksfunctions\n';
    delete functions.at(i);
    //std::cout << 'worksprofiles\n';
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
  if( varName == "NPU" )
    return  &ControlPlotsFunction::jetTruthEventNPU;
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
  if( varName == "GenJetResponse" && type == ControlPlotsConfig::L2L3Res )
    return  &ControlPlotsFunction::jetTruthEventResponseL2L3ResCorrected;
  if( varName == "GenJetResponse" && type == ControlPlotsConfig::L2L3L4 )
    return  &ControlPlotsFunction::jetTruthEventResponseL2L3L4Corrected;
  if( varName == "GenJetResponse" && type == ControlPlotsConfig::L2L3ResL4 )
    return  &ControlPlotsFunction::jetTruthEventResponseL2L3ResL4Corrected;
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
  if( varName == "RunNumber" )
    return  &ControlPlotsFunction::twoJetsPtBalanceEventRunNumber;
  if( varName == "Phi" )
    return  &ControlPlotsFunction::twoJetsPtBalanceEventJetPhi;
  if( varName == "Eta" )
    return  &ControlPlotsFunction::twoJetsPtBalanceEventJetEta;
  if( varName == "AbsEta" )
    return  &ControlPlotsFunction::twoJetsPtBalanceEventJetAbsEta;
  if( varName == "Jet2Eta")
   return  &ControlPlotsFunction::twoJetsPtBalanceEventJet2Eta; 
  if( varName == "Jet2AbsEta")
   return  &ControlPlotsFunction::twoJetsPtBalanceEventJet2AbsEta; 
  if( varName == "Pt" && type == ControlPlotsConfig::Uncorrected )
   return  &ControlPlotsFunction::twoJetsPtBalanceEventJetPt; 
  if( varName == "Pt" && type == ControlPlotsConfig::L2L3  )
   return  &ControlPlotsFunction::twoJetsPtBalanceEventJetPtL2L3Corrected; 
  if( varName == "Pt" && type == ControlPlotsConfig::L2L3Res  )
   return  &ControlPlotsFunction::twoJetsPtBalanceEventJetPtL2L3ResCorrected; 
  if( varName == "Jet2Pt" && type == ControlPlotsConfig::Uncorrected  )
   return  &ControlPlotsFunction::twoJetsPtBalanceEventJet2Pt; 
  if( varName == "Jet2Pt" && type == ControlPlotsConfig::L2L3  )
   return  &ControlPlotsFunction::twoJetsPtBalanceEventJet2PtL2L3Corrected; 
  if( varName == "Jet2Pt" && type == ControlPlotsConfig::L2L3Res  )
   return  &ControlPlotsFunction::twoJetsPtBalanceEventJet2PtL2L3ResCorrected; 
  if( varName == "JetLeadPt" && type == ControlPlotsConfig::Uncorrected  )
   return  &ControlPlotsFunction::twoJetsPtBalanceEventJetLeadPt; 
  if( varName == "JetLeadPt" && type == ControlPlotsConfig::L2L3  )
   return  &ControlPlotsFunction::twoJetsPtBalanceEventJetLeadPtL2L3Corrected; 
  if( varName == "JetLeadPt" && type == ControlPlotsConfig::L2L3Res  )
   return  &ControlPlotsFunction::twoJetsPtBalanceEventJetLeadPtL2L3ResCorrected; 
  if( varName == "JetLead2Pt" && type == ControlPlotsConfig::Uncorrected  )
   return  &ControlPlotsFunction::twoJetsPtBalanceEventJetLead2Pt; 
  if( varName == "JetLead2Pt" && type == ControlPlotsConfig::L2L3  )
   return  &ControlPlotsFunction::twoJetsPtBalanceEventJetLead2PtL2L3Corrected; 
  if( varName == "JetLead2Pt" && type == ControlPlotsConfig::L2L3Res  )
   return  &ControlPlotsFunction::twoJetsPtBalanceEventJetLead2PtL2L3ResCorrected;   
//  if( varName == "MeanPt")
//   return  &ControlPlotsFunction::twoJetsPtBalanceEventJet2PtL2L3Corrected; 
//  //    return  &ControlPlotsFunction::twoJetsPtBalanceEventMeanPt;
  if( varName == "MeanPt")
    return  &ControlPlotsFunction::twoJetsPtBalanceEventMeanPt;
  if( varName == "EMF" )
    return  &ControlPlotsFunction::twoJetsPtBalanceEventJetEMF;
  if( varName == "momentEtaEta" )
    return  &ControlPlotsFunction::twoJetsPtBalanceEventJetMomentEtaEta;
  if( varName == "momentPhiPhi" )
    return &ControlPlotsFunction::twoJetsPtBalanceEventJetMomentPhiPhi; 
  if( varName == "meanMoment" )
    return &ControlPlotsFunction::twoJetsPtBalanceEventJetMeanMoment;
  if( varName == "VtxN" )
    return &ControlPlotsFunction::twoJetsPtBalanceEventVtxN;
  if( varName == "MCNPUVtx" )
    return &ControlPlotsFunction::twoJetsPtBalanceEventMCNPUVtx;
  if( varName == "MCNPUTruth" )
    return &ControlPlotsFunction::twoJetsPtBalanceEventMCNPUTruth;
  if( varName == "PF_CH_Fraction" )
    return &ControlPlotsFunction::twoJetsPtBalanceEventPF_CH_Fraction;
  if( varName == "PF_NH_Fraction" )
    return &ControlPlotsFunction::twoJetsPtBalanceEventPF_NH_Fraction;
  if( varName == "PF_PH_Fraction" )
    return &ControlPlotsFunction::twoJetsPtBalanceEventPF_PH_Fraction;
  if( varName == "PF_EL_Fraction" )
    return &ControlPlotsFunction::twoJetsPtBalanceEventPF_EL_Fraction;
  if( varName == "PF_HFHad_Fraction" )
    return &ControlPlotsFunction::twoJetsPtBalanceEventPF_HFHad_Fraction;
  if( varName == "PF_HFEm_Fraction" )
    return &ControlPlotsFunction::twoJetsPtBalanceEventPF_HFEm_Fraction;
  if( varName == "Flavor" )
    return &ControlPlotsFunction::twoJetsPtBalanceEventJetFlavor;
  if( varName == "DeltaPhi" )
    return &ControlPlotsFunction::twoJetsPtBalanceEventDeltaPhi;
  if( varName == "ThirdJetFraction" )
    return &ControlPlotsFunction::twoJetsPtBalanceEventThirdJetFraction;
  if( varName == "ThirdJetFractionPlain" )
    return &ControlPlotsFunction::twoJetsPtBalanceEventThirdJetFractionPlain;
  if( varName == "ThirdJetPt" )
    return &ControlPlotsFunction::twoJetsPtBalanceEventThirdJetPt;
  if( varName == "Asymmetry" && type == ControlPlotsConfig::Uncorrected )
    return &ControlPlotsFunction::twoJetsPtBalanceEventAsymmetry;
  if( varName == "Asymmetry" && type == ControlPlotsConfig::Kalibri )
    return &ControlPlotsFunction::twoJetsPtBalanceEventAsymmetryKalibriCorrected;
  if( varName == "Asymmetry" && type == ControlPlotsConfig::L2L3 )
    return  &ControlPlotsFunction::twoJetsPtBalanceEventAsymmetryL2L3Corrected;
  if( varName == "Asymmetry" && type == ControlPlotsConfig::L2L3Res )
    return  &ControlPlotsFunction::twoJetsPtBalanceEventAsymmetryL2L3ResCorrected;
  if( varName == "Asymmetry" && type == ControlPlotsConfig::L2L3L4 )
    return  &ControlPlotsFunction::twoJetsPtBalanceEventAsymmetryL2L3L4Corrected;
  if( varName == "Asymmetry" && type == ControlPlotsConfig::L2L3ResL4 )
    return  &ControlPlotsFunction::twoJetsPtBalanceEventAsymmetryL2L3ResL4Corrected;
  if( varName == "B" && type == ControlPlotsConfig::Uncorrected )
    return &ControlPlotsFunction::twoJetsPtBalanceEventB;
  if( varName == "B" && type == ControlPlotsConfig::Kalibri )
    return &ControlPlotsFunction::twoJetsPtBalanceEventBKalibriCorrected;
  if( varName == "B" && type == ControlPlotsConfig::L2L3 )
    return  &ControlPlotsFunction::twoJetsPtBalanceEventBL2L3Corrected;
  if( varName == "B" && type == ControlPlotsConfig::L2L3Res )
    return  &ControlPlotsFunction::twoJetsPtBalanceEventBL2L3ResCorrected;
  if( varName == "B" && type == ControlPlotsConfig::L2L3L4 )
    return  &ControlPlotsFunction::twoJetsPtBalanceEventBL2L3L4Corrected;
  if( varName == "B" && type == ControlPlotsConfig::L2L3ResL4 )
    return  &ControlPlotsFunction::twoJetsPtBalanceEventBL2L3ResL4Corrected;
  if( varName == "MPFResponse" && type == ControlPlotsConfig::Uncorrected )
    return &ControlPlotsFunction::twoJetsPtBalanceEventMPFResponse;
  if( varName == "MPFResponse" && type == ControlPlotsConfig::L2L3 )
    return  &ControlPlotsFunction::twoJetsPtBalanceEventMPFResponseL2L3Corrected;
  if( varName == "MPFResponse" && type == ControlPlotsConfig::L2L3Res )
    return  &ControlPlotsFunction::twoJetsPtBalanceEventMPFResponseL2L3ResCorrected;
  if( varName == "") {
    return 0;
  }
  
  std::cerr << "ControlPlots: unknown variable " << varName << std::endl;
  return 0;
}



//!  Set style option for the output.
//---------------------------------------------------------------
void ControlPlots::setGStyle() const {
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

