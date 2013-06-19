#ifndef BasePlotExtractor_cc
#define BasePlotExtractor_cc

#include "BasePlotExtractor.h"
#include "../scripts/tdrstyle_mod.C"
#include <string>
#include <fstream>
#include <iostream>
#include <sstream>
#include <iomanip>

//! Default constructor, reads in information about the 
//! plots. These plots have to be defined in config_
//! and ExternalConfig_ read in here
//!
//!
//!  \author Henning Kirschenmann
//!  \date 2012/03/07
// ----------------------------------------------------------------   
BasePlotExtractor::BasePlotExtractor(TString plotsnames,TString kalibriPlotsShortName, TString ExternalConfigPath) {
  kalibriPlotsShortName_=kalibriPlotsShortName;
  outputPathROOT_=gSystem->pwd()+(TString)"/"+MakeDateDir();//+"/Output_"+plotsnames+".root";
  std::cout << outputPathROOT_ << std::endl;
  if(chdir(kalibriPlotsShortName) != 0){ 
    mkdir(kalibriPlotsShortName, S_IRWXU|S_IRWXG|S_IRWXO); 
    chdir(kalibriPlotsShortName); 
  } 
  if(chdir(plotsnames) != 0){ 
    mkdir(plotsnames, S_IRWXU|S_IRWXG|S_IRWXO); 
    chdir(plotsnames); 
  } 
  //  chdir("..");

  plotsnames_=plotsnames;
  //read in config files in kalibri-style
  ExternalConfig_=ConfigFile(ExternalConfigPath.Data());
  //  ExternalConfig_=ConfigFile("/afs/naf.desy.de/user/k/kirschen/scratch/2012_05_L2L3ResidualsFinal/L2andJERScripts/ExternalConfigs.cfg");
  //  ExternalConfig_=ConfigFile("/afs/naf.desy.de/user/k/kirschen/scratch/2012_05_L2L3ResidualsFinal/L2andJERScripts/ExternalConfigsFine.cfg");
  //  ExternalConfig_=ConfigFile("/afs/naf.desy.de/user/k/kirschen/public/ExternalConfigs.cfg");
  readInExtraInfo();
  std::cout << "readInExtraInfo(); executed " << yProfileTitle()<< std::endl;
  std::cout << "readInExtraInfo(); executed " << binningSelection_<< std::endl;
  std::cout << "opening L2L3.cfg: " << pathToConfig_ << std::endl;
  std::cout << "opening KalibriPlots.root from: " <<   kalibriPlotsPath_ << std::endl;
  config_=ConfigFile((std::string)pathToConfig_);//"/afs/naf.desy.de/user/k/kirschen/scratch/2012_05_L2L3ResidualsFinal/L2andJERScripts/AllPlots_L2L3.cfg");


  DeviationTypes_.push_back(std::make_pair("RatioDev"," "+mcLabel_+"/"+dataLabel_+" ratio"));
  DeviationTypes_.push_back(std::make_pair("MCDev"," in "+mcLabel_+""));
  DeviationTypes_.push_back(std::make_pair("DataDev"," in "+dataLabel_+""));

}

//! Manually read in extra information from
//! ExternalConfig_ to have nice x,y axis titles
//! and consistent labels (store in class members)
//!
//!  \author Henning Kirschenmann
//!  \date 2012/03/07
// ----------------------------------------------------------------   
void BasePlotExtractor::readInExtraInfo() {
  TString fileNames ="AK5PF";
  jetType_ = util::LabelFactory::jetAlgo(fileNames)+util::LabelFactory::jetType(fileNames);
  jetLabel_ = util::LabelFactory::labelJet(fileNames);
  std::cout << "import " <<(std::string)plotsnames_+" plots profile yMinMax" << std::endl;
  dataLabel_ = ExternalConfig_.read<std::string>((std::string)kalibriPlotsShortName_+" dataLabel",ExternalConfig_.read<std::string>("Default dataLabel","Data"));
  mcLabel_   = ExternalConfig_.read<std::string>((std::string)kalibriPlotsShortName_+" mcLabel",ExternalConfig_.read<std::string>("Default mcLabel","MC"));

  yProfileMinMax_ = bag_of<double>(ExternalConfig_.read<std::string>((std::string)plotsnames_+" plots profile yMinMax","0 0.5"));
  assert(yProfileMinMax_.size()==2);
  yProfileTitle_ = ExternalConfig_.read<std::string>((std::string)plotsnames_+" plots profile yTitle","DUMMYResolution");
  yProfileTitle_.ReplaceAll("Data",dataLabel_);
  yProfileTitle_.ReplaceAll("MC",mcLabel_);
  yRatioMinMax_ = bag_of<double>(ExternalConfig_.read<std::string>((std::string)plotsnames_+" plots ratio yMinMax","0.5 1.5"));
  assert(yRatioMinMax_.size()==2);
  yMCDataRatioMinMax_ = bag_of<double>(ExternalConfig_.read<std::string>((std::string)plotsnames_+" plots MCData ratio yMinMax","0.5 1.5"));
  assert(yMCDataRatioMinMax_.size()==2);
  yRatioTitle_ = ExternalConfig_.read<std::string>((std::string)plotsnames_+" plots ratio yTitle","DUMMY Ratio");
  yRatioTitle_.ReplaceAll("Data",dataLabel_);
  yRatioTitle_.ReplaceAll("MC",mcLabel_);
  yDifferenceMinMax_ = bag_of<double>(ExternalConfig_.read<std::string>((std::string)plotsnames_+" plots difference yMinMax","-15 10"));
  assert(yDifferenceMinMax_.size()==2);
  yDifferenceTitle_ = ExternalConfig_.read<std::string>((std::string)plotsnames_+" plots difference yTitle","DUMMY Difference");
  yDifferenceTitle_.ReplaceAll("Data",dataLabel_);
  yDifferenceTitle_.ReplaceAll("MC",mcLabel_);
  profileType_ = ExternalConfig_.read<std::string>((std::string)plotsnames_+" plots profile type","Mean");

  kalibriPlotsPath_ = ExternalConfig_.read<std::string>((std::string)kalibriPlotsShortName_+" kalibriPlotsPath","KalibriPlots.root");
  pathToConfig_ = ExternalConfig_.read<std::string>((std::string)kalibriPlotsShortName_+" pathToConfig","/afs/naf.desy.de/user/k/kirschen/scratch/2012_05_L2L3ResidualsFinal/L2andJERScripts/AllPlots_L2L3.cfg");
  makeEquiDistHistos_ = ExternalConfig_.read<bool>((std::string)kalibriPlotsShortName_+" makeEquiDistHistos",ExternalConfig_.read<bool>("Default makeEquiDistHistos",1));
  intLumi_ = ExternalConfig_.read<double>((std::string)kalibriPlotsShortName_+" intLumi",ExternalConfig_.read<double>("Default intLumi",0));
  sqrtS_ = ExternalConfig_.read<int>((std::string)kalibriPlotsShortName_+" SqrtS",ExternalConfig_.read<int>("Default SqrtS",7));

  runNumbers_ = bag_of<double>(ExternalConfig_.read<std::string>((std::string)kalibriPlotsShortName_+" RunNumbers",ExternalConfig_.read<std::string>("TimeDependence RunNumbers","")));
  runNumbersLabels_ = bag_of_string(ExternalConfig_.read<std::string>((std::string)kalibriPlotsShortName_+" RunNumbersLabels",ExternalConfig_.read<std::string>("TimeDependence RunNumbersLabels","")));

  runNumbersMap_.clear();

  for(unsigned int i=0;i<runNumbers_.size();i++){
    if(DEBUG)    std::cout << runNumbers_.size() << " and " << runNumbersLabels_.size()<<std::endl;
    if(DEBUG)    std::cout << runNumbers_.at(i) << " and " << runNumbersLabels_.at(i) <<std::endl;
  }
  assert(runNumbers_.size()==runNumbersLabels_.size());
  for(unsigned int i=0;i<runNumbers_.size();i++){
    runNumbersMap_[runNumbers_.at(i)]=runNumbersLabels_.at(i);
  }

  //  binningSelection_ = ExternalConfig_.read<std::string>("binning selection","kostas");
  if(kalibriPlotsPath_.Contains("kostas"))binningSelection_ = "kostas";
  else if(kalibriPlotsPath_.Contains("k_HFfix"))binningSelection_ = "k_HFfix";

}

BasePlotExtractor::~BasePlotExtractor() {
  std::cout << "calling the destructor of BasePlotExtractor " << kalibriPlotsShortName() << " and " << plotsnames_ << std::endl;
  //  assert(AllPlots_.size()==AllTwoDPlots_.size());
  for(unsigned int i=0;i<AllPlots_.size();i++){
    //    assert(AllPlots_.at(i).size()==AllTwoDPlots_.at(i).size());
    for(unsigned int j=0;j<AllPlots_.at(i).size();j++){
      for(unsigned int k=0;k<AllPlots_.at(i).at(j).size();k++){
	if(DEBUG)	std::cout << "trying to delete i "  << i << " j " << j << " k " << k << std::endl;
	delete AllPlots_.at(i).at(j).at(k);
      }
    }
  }
  std::cout << "deleted allplots vector" << std::endl;

  for(unsigned int i=0;i<AllTwoDPlots_.size();i++){
    for(unsigned int j=0;j<AllTwoDPlots_.at(i).size();j++){
      for(unsigned int k=0;k<AllTwoDPlots_.at(i).at(j).size();k++){
	delete AllTwoDPlots_.at(i).at(j).at(k);
      }
    }
  }
  std::cout << "deleted alltwodplots vector" << std::endl;

  std::cout << "AllRatiosDataMC_.size() " << AllRatiosDataMC_.size() << " AllDifferencesDataMC_.size() " << AllDifferencesDataMC_.size() << std::endl;
  assert(AllRatiosDataMC_.size()==AllDifferencesDataMC_.size());
  std::cout << "deleting AllRatiosDataMC_ and AllDifferencesDataMC_ now" << std::endl;
  for(unsigned int i=0;i<AllRatiosDataMC_.size();i++){
    for(unsigned int j=0;j<AllRatiosDataMC_.at(i).size();j++){
      //      std::cout << "deleting AllRatiosDataMC_ and AllDifferencesDataMC_ now" << std::endl;
      delete AllRatiosDataMC_.at(i).at(j);
      delete AllDifferencesDataMC_.at(i).at(j);
    }
  }
  std::cout << "configs_.size() " << configs_.size() << " functions_.size() " << functions_.size() << " profiles_.size() " << profiles_.size() << std::endl;
  std::cout << "deleted AllRatiosDataMC and AllDifferencesDataMC vector" << std::endl;

  assert (configs_.size()==functions_.size());
  assert (configs_.size()==profiles_.size());
  for(unsigned int i=0;i<configs_.size();i++){
    delete configs_.at(i);
    delete functions_.at(i);
    delete profiles_.at(i);
  }
  std::cout << "deleted configs functions profiles vector" << std::endl;
  for(unsigned int i=0;i<RatioVsBinVarHistos_.size();i++){
    delete RatioVsBinVarHistos_.at(i);
    delete MCVsBinVarHistos_.at(i);
    delete DataVsBinVarHistos_.at(i);
  }
  std::cout << "deleted RatioVsBinVarHistos vector" << std::endl;

  for(unsigned int i=0;i<AllDeviationsVsBinVarHistos_.size();i++){
    for(unsigned int j=0;j<AllDeviationsVsBinVarHistos_.at(i).size();j++){
      delete AllDeviationsVsBinVarHistos_.at(i).at(j);
    }
  }
  std::cout << "deleted AllDeviationsVsBinVarHistos vector" << std::endl;

  delete inf_;  //  inf_->Close();

}

//! Read in the basic plots for the profile types
//! and store them in class members
//!   AllRatiosDataMC_ contains all ratio plots (ratio->Divide(DataMCProfiles.at(0),DataMCProfiles.at(1));)
//!   AllDifferencesDataMC_ contains differences between Data and MC as histo, scaled to some %-scale ( 100 * DataMCProfiles.at(0) - 100 * DataMCProfiles.at(1));)
//!   AllPlots_        contains all profile plots (best: 2, one Data and one MC, but depends on definition of correction types); 
//!    - OneBinAsymmetryVsPt30 correction types  =  L2L3###
//!    - OneBinAsymmetryVsPt30 1 correction types  =  L2L3
//!    - OneBinAsymmetryVsPt30 input samples     =  0:data;1:MC
//!    - would be a good definition resulting in DataMCProfiles to have two entries
//!    - OneBinAsymmetryVsPt30 correction types  =  L2L3 L2L3res
//!    - would result in three entries
//!
//!   Following for Kalibri bookkeeping of plots (can be used later to extract e.g. labels)
//!   names_	   = names;	
//!   configs_	   = configs;	
//!   functions_	   = functions;	
//!   profiles_        = profiles;    
//! 
//!  \author Henning Kirschenmann
//!  \date 2012/03/07
// ----------------------------------------------------------------   
void BasePlotExtractor::init(TString profileType) {
  //  TFile *inf_;
  //      inf_ = new TFile("KalibriPlots_withOldSmear.root");
  //       inf_ = new TFile("KalibriPlots_noSmear.root");
  //    inf_ = new TFile("KalibriPlots.root");
  inf_ = new TFile(kalibriPlotsPath());
  std::vector<std::vector<Event*>* > samples;
  setTDRStyle();


  std::vector<std::string> names;
  names = bag_of_string(ExternalConfig_.read<std::string>((std::string)plotsnames_+" plots names",""));
  names_	          = names;	
  std::vector<ControlPlotsConfig*> configs(names.size());
  std::vector<ControlPlotsFunction*> functions(names.size());
  std::vector<ControlPlotsProfile*> profiles(names.size());
  VecOfVecOfTH2vec_t AllTwoDPlots;
  VecOfVecOfTH1vec_t AllPlots;
  VecOfTH1vec_t AllRatiosDataMC;
  VecOfTH1vec_t AllDifferencesDataMC;

  // Read different control plot names
  for(size_t i = 0; i < names.size(); i++) {
    std::cout << " Creating plots '" << names.at(i) << "'\n";
    
    // Create ControlPlotsConfig    
    ControlPlotsConfig *pConfig = new ControlPlotsConfig(getConfig(),names.at(i));
    configs.at(i) = pConfig;
    configs.at(i)->setOutDirName((std::string)outputPathROOT_);
    configs.at(i)->setOutRootFileName((std::string)("Output"+kalibriPlotsShortName()+".root"));
    //    if(i==0)configs.at(i)->openRootFile() ;
    // Create functions
    ControlPlotsFunction *func = new ControlPlotsFunction();
    std::cout << "test" << std::endl;
    func->setBinFunction(findTwoJetsPtBalanceEventFunction(pConfig->binVariable()));
    func->setXFunction(findTwoJetsPtBalanceEventFunction(pConfig->xVariable()));
    func->setCutFunction(findTwoJetsPtBalanceEventFunction(pConfig->cutVariable()));
    func->setCut2Function(findTwoJetsPtBalanceEventFunction(pConfig->cut2Variable()));
    for(ControlPlotsConfig::InputTagsIterator it = pConfig->inputTagsBegin() ; 
	it != pConfig->inputTagsEnd(); ++it) {
      func->addYFunction(it->second,findTwoJetsPtBalanceEventFunction(pConfig->yVariable(),it->second));
    }
    functions.at(i) = func;
    std::cout << "test2" << std::endl;
    
    // Create profile
    profiles.at(i) = new ControlPlotsProfile(pConfig,func);
    
    //    std::cout << pConfig->nBins()<< std::endl;
    VecOfTH1vec_t BinsDMCF;
    VecOfTH2vec_t BinsTwoDDMCF;
    TH1vec_t DataMCRatios; //Vector to contain Data/MC-ratios
    TH1vec_t DataMCDifferences; //Vector to contain Data-MC differences
    for(int bin_i=0;bin_i<pConfig->nBins();bin_i++){
      //      std::cout << bin_i << " of " << pConfig->nBins() << std::endl;
      TH1vec_t DataMCProfiles;//Vector to contain Data and MC TH1
      TH2vec_t DataMCTwoDHistos;//Vector to contain Data and MC TH2
      for(ControlPlotsConfig::InputTagsIterator i = pConfig->inputTagsBegin();
	  i != pConfig->inputTagsEnd(); ++i) {
	TString name = pConfig->name();
	name += "/" +pConfig->name();
	name += "_"+pConfig->yVariable()+"Vs"+pConfig->xVariable();
	name += "_"+pConfig->sampleName(i->first);
	name += "_"+pConfig->correctionTypeName(i->second);
	name += "_"+pConfig->binName(bin_i);
	TString TwoDname = name;
	name += "_"+profileType;
	//Import the profile plots (TH1D) as defined in the config (usually two, depending on number of correction types)
	//Replace histos with simple equidistant binning for runnumber histograms (depending on default in ExternalConfigs.cfg)
	if(makeEquiDistHistos_&&pConfig->xVariable()=="RunNumber")DataMCProfiles.push_back(replaceHistosWithEquiDistBins((TH1D*) inf_->Get(name)));
	else{
	  DataMCProfiles.push_back((TH1D*) inf_->Get(name));
       	}
	DataMCTwoDHistos.push_back((TH2D*) inf_->Get(TwoDname));
	if(DEBUG)	std::cout << DataMCProfiles.back()->GetName() << std::endl;
      }
      DataMCProfiles.at(0)->Sumw2();
      DataMCProfiles.at(1)->Sumw2();
      // Do ratios of first and second (0 should be Data, 1 should be MC, can be checked in output above)
      TH1D* ratio = (TH1D*) DataMCProfiles.at(0)->Clone();
      ratio->Sumw2();
      ratio->Divide(DataMCProfiles.at(0),DataMCProfiles.at(1));
      DataMCRatios.push_back(ratio);
      TH1D* difference = (TH1D*) DataMCProfiles.at(0)->Clone();
      difference->Sumw2();
      difference->Add(difference, DataMCProfiles.at(1), 100, -100);
      DataMCDifferences.push_back(difference);
      //      DataMCProfiles.at(1)->Draw("hist");
      //      DataMCProfiles.at(0)->Draw("same");
      //      c->SaveAs(DataMCProfiles.back()->GetName()+(TString)".eps");
      BinsDMCF.push_back(DataMCProfiles);
      BinsTwoDDMCF.push_back(DataMCTwoDHistos);
    }// (TH1D*)inf_->Get("kFSR_vs_Abseta_histo_res1");
    AllPlots.push_back(BinsDMCF);
    AllTwoDPlots.push_back(BinsTwoDDMCF);
    AllRatiosDataMC.push_back(DataMCRatios);
    AllDifferencesDataMC.push_back(DataMCDifferences);
  }


  AllRatiosDataMC_        = AllRatiosDataMC;
  AllDifferencesDataMC_   = AllDifferencesDataMC;
  AllPlots_               = AllPlots;
  AllTwoDPlots_           = AllTwoDPlots;
  //  names_	          = names;	
  configs_	          = configs;	
  functions_	          = functions;	
  profiles_               = profiles;    
  
  makeRatioVsBinVarHistos();
}

//! Refresh DataMCRatios vector after performing e.g. an  
//! extrapolation (see Extrapolation.cc /.h) 
//! This is needed to get the ratios and any following plots with the 
//! (e.g. radiation) extrapolation appplied
//! 
//!  \author Henning Kirschenmann
//!  \date 2012/03/07
// ----------------------------------------------------------------   
void BasePlotExtractor::refreshRatiosDataMC(){
  //delete old stuff
  std::cout << "AllRatiosDataMC_.size() " << AllRatiosDataMC_.size() << " AllDifferencesDataMC_.size() " << AllDifferencesDataMC_.size() << std::endl;
  assert(AllRatiosDataMC_.size()==AllDifferencesDataMC_.size());
  std::cout << "deleting AllRatiosDataMC_ and AllDifferencesDataMC_ now" << std::endl;
  for(unsigned int i=0;i<AllRatiosDataMC_.size();i++){
    for(unsigned int j=0;j<AllRatiosDataMC_.at(i).size();j++){
      //      std::cout << "deleting AllRatiosDataMC_ and AllDifferencesDataMC_ now" << std::endl;
      delete AllRatiosDataMC_.at(i).at(j);
      delete AllDifferencesDataMC_.at(i).at(j);
    }
  }
  std::cout << "configs_.size() " << configs_.size() << " functions_.size() " << functions_.size() << " profiles_.size() " << profiles_.size() << std::endl;
  std::cout << "deleted AllRatiosDataMC and AllDifferencesDataMC vector" << std::endl;




  VecOfTH1vec_t AllRatiosDataMC;
  VecOfTH1vec_t AllDifferencesDataMC;

  // Read different control plot names
  for(size_t i = 0; i < names_.size(); i++) {
    TH1vec_t DataMCRatios; //Vector to contain Data/MC-ratios
    TH1vec_t DataMCDifferences; //Vector to contain Data-MC differences
    for(int bin_i=0;bin_i<configs_.at(i)->nBins();bin_i++){
      // Do ratios of first and second (0 should be Data, 1 should be MC, can be checked in output above)
      TH1D* ratio = (TH1D*) AllPlots_.at(i).at(bin_i).at(0)->Clone();
      ratio->Sumw2();
      ratio->Divide(AllPlots_.at(i).at(bin_i).at(0),AllPlots_.at(i).at(bin_i).at(1));
      DataMCRatios.push_back(ratio);
      TH1D* difference = (TH1D*) AllPlots_.at(i).at(bin_i).at(0)->Clone();
      difference->Sumw2();
      difference->Add(difference, AllPlots_.at(i).at(bin_i).at(1), 100, -100);
      DataMCDifferences.push_back(difference);
    }
    AllRatiosDataMC.push_back(DataMCRatios);
    AllDifferencesDataMC.push_back(DataMCDifferences);
  }
  //save to class instance
  AllRatiosDataMC_ = AllRatiosDataMC;
  AllDifferencesDataMC_   = AllDifferencesDataMC;
  std::cout << "done with refreshing ratio-plots" <<std::endl;
}

//! Fits const and loglin function to DataMC-ratios (added to list of functions of ratio histogram)
//! then creates a RatioVsBinVarHisto using the const fit
//! All these histos are saved to RatioVsBinVarHistos_
//! 
//!  \author Henning Kirschenmann
//!  \date 2012/03/07
// ----------------------------------------------------------------   
void BasePlotExtractor::makeRatioVsBinVarHistos(){
  TH1vec_t RatioVsBinVarHistos;
  TH1vec_t MCVsBinVarHistos;
  TH1vec_t DataVsBinVarHistos;

  VecOfTH1vec_t AllDeviationsVsBinVarHistos;
  //  TH1vec_t DeviationsOfRatioVsBinVarHistos;
  //insert deviation plots here
  for(int conf_i=0;conf_i<configs_.size();conf_i++){
    //    std::cout << names_.size() << std::endl;
    //    std::cout << names_.at(conf_i) << std::endl;
    TH1D* RatioVsBinVar = new TH1D(("RatioVsBinVar"+names_.at(conf_i)).c_str(),"",configs_.at(conf_i)->nBins(),&(configs_.at(conf_i)->binEdges()->front()));
    RatioVsBinVar->Sumw2();
    RatioVsBinVar->GetXaxis()->SetTitle(configs_.at(conf_i)->binAxisTitle().c_str());
    RatioVsBinVar->GetYaxis()->SetTitle(dataLabel_+"/"+mcLabel_+" ratio (const fit)");

    TH1D* MCVsBinVar = new TH1D(("MCVsBinVar"+names_.at(conf_i)).c_str(),"",configs_.at(conf_i)->nBins(),&(configs_.at(conf_i)->binEdges()->front()));
    MCVsBinVar->Sumw2();
    MCVsBinVar->GetXaxis()->SetTitle(configs_.at(conf_i)->binAxisTitle().c_str());
    MCVsBinVar->GetYaxis()->SetTitle(mcLabel_+" resp.estimator (const fit)");

    TH1D* DataVsBinVar = new TH1D(("DataVsBinVar"+names_.at(conf_i)).c_str(),"",configs_.at(conf_i)->nBins(),&(configs_.at(conf_i)->binEdges()->front()));
    DataVsBinVar->Sumw2();
    DataVsBinVar->GetXaxis()->SetTitle(configs_.at(conf_i)->binAxisTitle().c_str());
    DataVsBinVar->GetYaxis()->SetTitle(dataLabel_ + " resp. estimator (const fit)");

    TH1vec_t DeviationsVsBinVarHistos;
    for(size_t dev_i =0; dev_i<DeviationTypes_.size();dev_i++){
      DeviationsVsBinVarHistos.push_back(new TH1D(("DeviationsOfRatioVsBinVar"+names_.at(conf_i)+"_"+DeviationTypes_.at(dev_i).first)/*.c_str()*/,(";"+configs_.at(conf_i)->binAxisTitle()+";Weighted standard deviation [%]").c_str(),configs_.at(conf_i)->nBins(),&(configs_.at(conf_i)->binEdges()->front())));
      //    TH1D* DeviationsOfRatioVsBinVar = new TH1D(("DeviationsOfRatioVsBinVar"+names_.at(conf_i)).c_str(),"",configs_.at(conf_i)->nBins(),&(configs_.at(conf_i)->binEdges()->front()));
//      DeviationsVsBinVarHistos.back()->GetXaxis()->SetTitle(configs_.at(conf_i)->binAxisTitle().c_str());
//      DeviationsVsBinVarHistos.back()->GetYaxis()->SetTitle("Weighted standard deviation [%]");
    }


    //  for(int bin_i=0;bin_i<configs_.at(0)->nBins()-1;bin_i++){
  for(int bin_i=0;bin_i<configs_.at(conf_i)->nBins();bin_i++){
    fitFunctionsToPlot(AllRatiosDataMC_.at(conf_i).at(bin_i));
    fillRatioVsBinVarPlot(RatioVsBinVar,AllRatiosDataMC_.at(conf_i).at(bin_i),bin_i,"fit_const");
    fitFunctionsToPlot(AllPlots_.at(conf_i).at(bin_i).at(1)); //mc
    fitFunctionsToPlot(AllPlots_.at(conf_i).at(bin_i).at(0)); //data

    fillRatioVsBinVarPlot(MCVsBinVar,  AllPlots_.at(conf_i).at(bin_i).at(1),bin_i,"fit_const");
    fillRatioVsBinVarPlot(DataVsBinVar,AllPlots_.at(conf_i).at(bin_i).at(0),bin_i,"fit_const");

    for(size_t dev_i =0; dev_i<DeviationTypes_.size();dev_i++){
      if(DeviationTypes_.at(dev_i).first=="RatioDev")fillDeviationsOfRatioVsBinVarPlot(DeviationsVsBinVarHistos.at(dev_i),AllRatiosDataMC_.at(conf_i).at(bin_i),bin_i,"fit_const");
      else if(DeviationTypes_.at(dev_i).first=="MCDev")fillDeviationsOfRatioVsBinVarPlot(DeviationsVsBinVarHistos.at(dev_i),AllPlots_.at(conf_i).at(bin_i).at(1),bin_i,"fit_const");
      else if(DeviationTypes_.at(dev_i).first=="DataDev")fillDeviationsOfRatioVsBinVarPlot(DeviationsVsBinVarHistos.at(dev_i),AllPlots_.at(conf_i).at(bin_i).at(0),bin_i,"fit_const");
      
    }
  }
//  for(int bin_i=0;bin_i<configs_.at(0)->nBins()-1;bin_i++){
//    std::cout << RatioVsBinVar->GetBinContent(bin_i+1) << ", " << std::endl;
//  }
  RatioVsBinVarHistos.push_back(RatioVsBinVar);
  MCVsBinVarHistos  .push_back(MCVsBinVar  );
  DataVsBinVarHistos.push_back(DataVsBinVar);


  AllDeviationsVsBinVarHistos.push_back(DeviationsVsBinVarHistos);
  //  DeviationsOfRatioVsBinVarHistos.push_back(DeviationsOfRatioVsBinVar);  
  }
  RatioVsBinVarHistos_=RatioVsBinVarHistos;
  MCVsBinVarHistos_=MCVsBinVarHistos;
  DataVsBinVarHistos_=DataVsBinVarHistos;
  AllDeviationsVsBinVarHistos_=AllDeviationsVsBinVarHistos;
  //DeviationsOfRatioVsBinVarHistos_=DeviationsOfRatioVsBinVarHistos;
}


//! Fit functions to ratio plot
//! 
//! 
//!  \author Henning Kirschenmann
//!  \date 2012/07/11
// ----------------------------------------------------------------   
void BasePlotExtractor::fitFunctionsToPlot(TH1D* histo){
  TF1 *fit_const = new TF1("fit_const","[0]",histo->GetXaxis()->GetXmin(),histo->GetXaxis()->GetXmax()); //was used before...
  fit_const->SetParameters(1,1);
  fit_const->SetParName(0,"const");
  histo->Fit("fit_const","EMQ","same");

  if(!plotsnames_.Contains("Resolution")){//ignore pt-dependence of
    //    resolution ratios for my purposes... might change this to a
    //  config-file option in the future
    if(plotsnames_.Contains("VsPt")||plotsnames_.Contains("VsMeanPt")||plotsnames_.Contains("VsJetLeadPt")||plotsnames_.Contains("VsJetLead2Pt")||plotsnames_.Contains("VsJet2Pt")){
      TF1 *fit_loglin = new TF1("fit_loglin","[0]+[1]*TMath::Log(x)",histo->GetXaxis()->GetXmin(),histo->GetXaxis()->GetXmax()); //was used before...
      fit_loglin->SetParameters(1,1);
      fit_loglin->SetParName(0,"const");
      fit_loglin->SetParName(1,"slope");
      histo->Fit("fit_loglin","EMQ+","same");
    }
    else{
      TF1 *fit_lin = new TF1("fit_lin","[0]+[1]*x",histo->GetXaxis()->GetXmin(),histo->GetXaxis()->GetXmax()); //was used before...
      fit_lin->SetParameters(1,1);
      fit_lin->SetParName(0,"const");
      fit_lin->SetParName(1,"slope");
      histo->Fit("fit_lin","EMQ+","same");
    }
  }
}

//! Fill RatioVsBinVarPlotWith fit results
//! 
//! 
//!  \author Henning Kirschenmann
//!  \date 2012/07/11
// ----------------------------------------------------------------   
void BasePlotExtractor::fillRatioVsBinVarPlot(TH1D* RatioVsBinVarHisto, TH1D* HistoOfBin_i, Int_t bin_i, TString func_name){
  if(HistoOfBin_i->GetFunction(func_name)){
    RatioVsBinVarHisto->SetBinContent(bin_i+1,HistoOfBin_i->GetFunction(func_name)->GetParameter(0));
    RatioVsBinVarHisto->SetBinError(bin_i+1,HistoOfBin_i->GetFunction(func_name)->GetParError(0));
    //    std::cout <<RatioVsBinVarHisto->GetBinContent(bin_i+1)<<std::endl;
  }
  else {
    std::cout << RatioVsBinVarHisto->GetName() << "WARNING: NO FIT FUNCTIONS ADDED... SETTING BIN TO ZERO: " << bin_i<< std::endl;
    RatioVsBinVarHisto->SetBinContent(bin_i+1,0.);
    RatioVsBinVarHisto->SetBinError(bin_i+1,0.);
  }


}

//! Fill DeviationsOfRatioVsBinVarPlotWith fit results
//! Current implemenation: determine weighted standard deviation from (constant) fit, set bin error to zero (was: paerror(0) - error of constant fit parameter)
//! 
//!  \author Henning Kirschenmann
//!  \date 2012/07/11
// ----------------------------------------------------------------   
void BasePlotExtractor::fillDeviationsOfRatioVsBinVarPlot(TH1D* DeviationsOfRatioVsBinVarHisto, TH1D* HistoOfBin_i, Int_t bin_i, TString func_name){
  if(HistoOfBin_i->GetFunction(func_name)){

    std::vector<Jet::JetIndex*> DeviationIndices;


    std::vector <Double_t> RatioValues;
    std::vector <Double_t> RatioErrorValues;
    std::vector <Double_t> DeviationsFromFit;
    for(int i=1;i<HistoOfBin_i->GetNbinsX()+1;i++){
      if(HistoOfBin_i->GetBinContent(i)!=0){
	RatioValues.push_back(HistoOfBin_i->GetBinContent(i));
	RatioErrorValues.push_back(HistoOfBin_i->GetBinError(i));
	DeviationsFromFit.push_back(TMath::Abs(HistoOfBin_i->GetBinContent(i)-HistoOfBin_i->GetFunction(func_name)->GetParameter(0)));
	if(DEBUG)std::cout << HistoOfBin_i->GetName() << " " << HistoOfBin_i->GetBinContent(i) << "\t " <<HistoOfBin_i->GetFunction(func_name)->GetParameter(0) << "\t " <<  DeviationsFromFit.back() << std::endl;
      }
    }
    Double_t SumDeviations=0;
    for(int i=1;i<DeviationsFromFit.size();i++){
      if(DEBUG)std::cout << "DeviationsFromFit.at(i): " << DeviationsFromFit.at(i) << std::endl;
      SumDeviations+=DeviationsFromFit.at(i);
    }
    if(DeviationsFromFit.size()>0){
    assert(DeviationsFromFit.size()>0);
    SumDeviations = SumDeviations/DeviationsFromFit.size();
    SumDeviations = TMath::Sqrt(SumDeviations);
    }
    else std::cout << "WARNING: apparently no entries... something went wrong?" << std::endl;

    if(DEBUG)std::cout << "about to set bincontent for bin_i " << bin_i << std::endl;
    DeviationsOfRatioVsBinVarHisto->SetBinError(bin_i+1,0);
    Double_t fitResultConst = HistoOfBin_i->GetFunction(func_name)->GetParameter(0);
    //weighted standard deviation
    DeviationsOfRatioVsBinVarHisto->SetBinContent(bin_i+1,getWeightedStandardDeviation(RatioValues,RatioErrorValues)/fitResultConst*100);//DeviationsFromFit.size()>5 ? SumDeviations : 0);
    if(DEBUG)std::cout << "set bincontent for bin_i " << bin_i << std::endl;
	
	for(size_t i = 0; i < DeviationIndices.size(); ++i) {
	  delete DeviationIndices[i];
	}
    if(DEBUG)std::cout << "deleted DeviationIndices" << std::endl;
    
  }
  else {
    std::cout << DeviationsOfRatioVsBinVarHisto->GetName() << "WARNING: NO FIT FUNCTIONS ADDED... SETTING BIN TO ZERO: " << bin_i<< std::endl;
    DeviationsOfRatioVsBinVarHisto->SetBinContent(bin_i+1,0.);
    DeviationsOfRatioVsBinVarHisto->SetBinError(bin_i+1,0.);
  }

  
}

Double_t BasePlotExtractor::getSampleMean(std::vector <Double_t> sampleVector){
  Int_t size = sampleVector.size();
  if(size>1){
    Double_t Sum=0;
    for(int i=0;i<size;i++){
      Sum+=sampleVector.at(i);
    }
    if(DEBUG)    std::cout << "sampleMean: " << Sum/size << std::endl;
    return Sum/size;
  }
  else{
    std::cout << "single or no elements for mean calculation. Quitting" <<std::endl;
    //    assert(size>1);
    return -1;
  }
  
  
}

Double_t BasePlotExtractor::getSampleStandardDeviation(std::vector <Double_t> sampleVector){
  Int_t size = sampleVector.size();
  if(size>1){
    Double_t mean = getSampleMean(sampleVector);
    Double_t SumDiffSquared=0;
    for(int i=0;i<size;i++){
      SumDiffSquared+= TMath::Power(sampleVector.at(i)-mean,2);
    }
    if(DEBUG)    std::cout << "SumDiffSquared: " << SumDiffSquared << std::endl;
    if(DEBUG)    std::cout << "size: " << size << std::endl;
    Double_t sampleStandardDeviation=TMath::Sqrt(SumDiffSquared/(size-1));
    if(DEBUG)    std::cout << "sampleStandardDeviation: " << sampleStandardDeviation << std::endl;
    return sampleStandardDeviation;
    
  }
  else{
    return -1;
  }
  
  
}

Double_t BasePlotExtractor::getWeightedStandardDeviation(std::vector <Double_t> sampleVector, std::vector <Double_t> sampleVectorErrors){
  Int_t size = sampleVector.size();
  assert(sampleVector.size()==sampleVectorErrors.size());
  if(size>1){
    Double_t s0=0;
    Double_t s1=0;
    Double_t s2=0;
    for(int i=0;i<size;i++){
      Double_t weight = 1/sampleVectorErrors.at(i);
      s0 += weight;
      s1 += weight * TMath::Power(sampleVector.at(i),1);
      s2 += weight * TMath::Power(sampleVector.at(i),2);
    }

    Double_t WeightedStandardDeviation = TMath::Sqrt(s0*s2-TMath::Power(s1,2))/s0;

    return WeightedStandardDeviation;
    
  }
  else{
    return -1;
  }
  
  
}



//!  Get JER values from histos
//! 
//!  \author Kristin Heine/Henning Kirschenmann
//!  \date 2012/05/23
// ----------------------------------------------------------------
void BasePlotExtractor::outputTable(TString label, TH1D* histo){

   ofstream myfile;
   myfile.open("Residual_results_"+label+".txt");

   for(int i = 1; i < histo->GetNbinsX()+1; i++) {
     myfile << "Eta Bin: " <<std::setw(8)<< histo->GetXaxis()->GetBinLowEdge(i) <<std::setw(8) <<  histo->GetXaxis()->GetBinUpEdge(i) << " : " << "Residual: " << std::setw(12)<<histo->GetBinContent(i) << " +- " << std::setw(12)<<histo->GetBinError(i) << "\n"; 
   }
   myfile.close();

   myfile.open("ResidualCorrectionFormat_"+label+".txt");
   myfile << "{ 1 JetEta 1 JetPt [0] Correction L2Relative} \n";
   for(int i = 1; i < histo->GetNbinsX()+1; i++) {
     myfile << histo->GetXaxis()->GetBinLowEdge(i) <<"\t" <<  histo->GetXaxis()->GetBinUpEdge(i) << "\t        3              3           3500      "<< histo->GetBinContent(i) <<"\n";
   }
   myfile.close();

   myfile.open("TimeDepTableFormat_"+label+".txt");
   for(int i = 1; i < histo->GetNbinsX()+1; i++) {
     myfile  << histo->GetXaxis()->GetBinLowEdge(i) <<"\t" <<  histo->GetXaxis()->GetBinUpEdge(i) <<":\t" << std::setw(12) << histo->GetBinContent(i) << " \t+- " << std::setw(12)<<histo->GetBinError(i)<<"\n";
   }
   myfile.close();


}


//! Fills legend with information from ExternalConfig_, 
//! also adds chi^2/ndf to labels
//!
//!  \author Henning Kirschenmann
//!  \date 2012/03/07
// ----------------------------------------------------------------   
void BasePlotExtractor::addFunctionLabelsToLegend(TH1D* histo, TLegend* leg){
  DefaultStyles style;
  style.setStyle("Confidence");
  TIter next(histo->GetListOfFunctions());
  TObject *obj;
  Int_t i=0;
  while ((obj = next())){
    //    obj->Print();//Draw(next.GetOption());
    TF1* verytemp_tf1 = (TF1*) obj;
    verytemp_tf1->SetLineColor(style.getColor(i));
    verytemp_tf1->SetFillColor(style.getColor(i));

    //introduced memory leak here for plotting?!
    TF1* temp_tf1 = (TF1*) verytemp_tf1->Clone();
    temp_tf1->SetFillStyle(1001);

    TFitResultPtr r = histo->Fit(temp_tf1->GetName(),"NSEMQ");
    Double_t chisquared1;
    unsigned int ndf1;
    char buffer [50];
    chisquared1= r->Chi2();
    ndf1= r->Ndf();

    sprintf (buffer, " (#chi^{2}/ndf = %.1f / %u )", chisquared1, ndf1);
    std::cout << buffer << std::endl;

    std::string legendName = ExternalConfig_.read<std::string>(temp_tf1->GetName()+(std::string) " legend label","DUMMYlabel");
    //      yProfileTitle_ = ExternalConfig_.read<std::string>((std::string)plotsnames_+" plots profile yTitle","DUMMYResolution");
    leg->AddEntry(temp_tf1,(legendName+buffer).c_str(),"FL");
    i++;
  }

}


//! Draws confidence intervals for all functions in
//! the listoffunctions of the histogram
//! For this, it repeats the fit and then extracts the band
//!
//!  \author Henning Kirschenmann
//!  \date 2012/03/07
// ----------------------------------------------------------------   
void BasePlotExtractor::drawConfidenceIntervals(TH1D* histo){
  DefaultStyles style;
  style.setStyle("Confidence");
  
  //  std::cout << "List of functions: " << std::endl;
  //	   histo->GetListOfFunctions()->Print();
  TIter next(histo->GetListOfFunctions());
  TObject *obj;
  Int_t i=0;
  while ((obj = next())){
    //    obj->Print();//Draw(next.GetOption());
    TF1* temp_tf1 = (TF1*) obj;

    histo->Fit(temp_tf1->GetName(),"NEMQ");
    //Create a histogram to hold the confidence intervals
    TH1D *hint = new TH1D("hint","Fitted function with .95 conf.band", 100, histo->GetXaxis()->GetXmin(), histo->GetXaxis()->GetXmax());
    hint->Sumw2();
    (TVirtualFitter::GetFitter())->GetConfidenceIntervals(hint);
    //Now the "hint" histogram has the fitted function values as the 
    //bin contents and the confidence intervals as bin errors
    hint->SetStats(kFALSE);
    hint->SetMarkerStyle(1);
    //    if(i<fill_colors_.size()){
      hint->SetMarkerColor(style.getColor(i));
      hint->SetFillColor(style.getColor(i));
      temp_tf1->SetLineColor(style.getColor(i));
      temp_tf1->SetFillColor(style.getColor(i));
      //    }
      //    else{
      //      hint->SetMarkerColor(38);
      //      hint->SetFillColor(38);
      //    }
    i++;
    hint->DrawClone("e4 same");
    delete hint;
  }
  
}


//! Draws confidence intervals for all functions in
//! the listoffunctions of the histogram
//! For this, it repeats the fit and then extracts the band
//!
//!  \author Henning Kirschenmann
//!  \date 2012/03/07
// ----------------------------------------------------------------   
void BasePlotExtractor::drawConfidenceIntervals(TGraphErrors* histo){
  DefaultStyles style;
  style.setStyle("Confidence");
  
  //  std::cout << "List of functions: " << std::endl;
  //	   histo->GetListOfFunctions()->Print();
  TIter next(histo->GetListOfFunctions());
  TObject *obj;
  Int_t i=0;
  while ((obj = next())){
    //    obj->Print();//Draw(next.GetOption());
    TF1* temp_tf1 = (TF1*) obj;
    //    std::cout << "drawing confidence intervals" << std::endl;
    histo->Fit(temp_tf1,"NEMQ");
    //Create a histogram to hold the confidence intervals
    TH1D *hint = new TH1D("hint","Fitted function with .95 conf.band", 100, histo->GetXaxis()->GetXmin(), histo->GetXaxis()->GetXmax());
    hint->Sumw2();
    (TVirtualFitter::GetFitter())->GetConfidenceIntervals(hint);
    //Now the "hint" histogram has the fitted function values as the 
    //bin contents and the confidence intervals as bin errors
    hint->SetStats(kFALSE);
    hint->SetMarkerStyle(1);
    //    if(i<fill_colors_.size()){
      hint->SetMarkerColor(style.getColor(i));
      hint->SetFillColor(style.getColor(i));
      temp_tf1->SetLineColor(style.getColor(i));
      temp_tf1->SetFillColor(style.getColor(i));
      //    }
      //    else{
      //      hint->SetMarkerColor(38);
      //      hint->SetFillColor(38);
      //    }
    i++;
    hint->DrawClone("e4 same");
    delete hint;
  }
  
}

//!  \brief modify plots to have simply numbered equidistant binning
//!
//!  Useful e.g. for time dependence studies where the runnumbers 
//!  
//!  
//!  
// -------------------------------------------------------------

void BasePlotExtractor::drawRunNumberLabels(TH1D* histo, ControlPlotsConfig* config){
  
  DefaultStyles style;
  std::cout << "temp" << std::endl;
  std::vector<double> xBinEdges = *config->xBinEdges();

  TStyle *tdrStyle = (TStyle*)gROOT->FindObject("tdrStyle");

  for(unsigned int i=0;i<runNumbers_.size();i++){
    //    std::cout << runNumbers_.size() << " and " << runNumbersLabels_.size()<<std::endl;
    if(DEBUG)    if(xBinEdges.size()>i)std::cout << runNumbers_.at(i) << " and " << runNumbersLabels_.at(i) << " and " << xBinEdges.at(i)<< std::endl;
  }

  std::vector <TLine*> lines;
  std::vector <TText*> labels;
  std::vector <int> binNumbers;
  binNumbers.push_back(1);
  //  std::cout << histo->GetYaxis()->GetFirst() << " and " << histo->GetYaxis()->GetLast() << " " << histo->GetYaxis()->GetBinLowEdge(1) << " " << histo->GetYaxis()->GetBinUpEdge(1) << " " << histo->GetMinimum() << " " << histo->GetMaximum() <<std::endl;
  for(unsigned int xbin_i=0;xbin_i<xBinEdges.size();xbin_i++){
    //    std::cout << xBinEdges.at(xbin_i) <<std::endl;
    if(runNumbersMap_[xBinEdges.at(xbin_i)]!=""){
      //      std::cout << xBinEdges.at(xbin_i) << " also bin nummer " << xbin_i<< std::endl;
      binNumbers.push_back(xbin_i);
      lines.push_back(new TLine(histo->GetXaxis()->GetBinUpEdge(xbin_i),histo->GetMinimum(),histo->GetXaxis()->GetBinUpEdge(xbin_i),histo->GetMaximum()));
      labels.push_back(new TText(histo->GetXaxis()->GetBinUpEdge(xbin_i),histo->GetMinimum()+0.05*(histo->GetMaximum()-histo->GetMinimum()),(runNumbersMap_[xBinEdges.at(xbin_i)]).c_str()));

    }
  }
  for(unsigned int line_i=0;line_i<lines.size();line_i++){
    lines.at(line_i)->SetLineColor(style.getColor(line_i));
    lines.at(line_i)->SetLineStyle(2);
    lines.at(line_i)->DrawClone();
    labels.at(line_i)->SetX(histo->GetXaxis()->GetBinUpEdge(binNumbers.at(line_i)));
    labels.at(line_i)->SetTextColor(style.getColor(line_i));
    labels.at(line_i)->SetTextSize(  tdrStyle->GetTitleFontSize()/2);
    labels.at(line_i)->DrawClone();

  }

  for(unsigned int line_i=0;line_i<lines.size();line_i++){
    delete lines.at(line_i);
    delete labels.at(line_i);
  }
}


//!  \brief Encapsulate cmsPrelDraw from tdrstyle_mod.C
//!
//!  Chooses correct centre of mass energy and integrated luminosity
//!  
// -------------------------------------------------------------

void BasePlotExtractor::drawCMSPrel(bool wide){

  cmsPrel(intLumi_,wide, sqrtS_);

}


//!  \brief modify plots to have simply numbered equidistant binning
//!
//!  Useful e.g. for time dependence studies where the runnumbers 
//!  are widely spread and sparse
//!  
//!  
// -------------------------------------------------------------
TH1D* BasePlotExtractor::replaceHistosWithEquiDistBins(TH1D* histo){


  //continue here...

  Int_t nbins = histo->GetNbinsX();
  std::vector <Double_t> binEdges;
  for(int i=0;i<nbins+1;i++){
    binEdges.push_back(i);
    if(i<nbins+1&&i>1&&TString(histo->GetName()).Contains("MC")){
      histo->SetBinContent(i,histo->GetBinContent(i-1));
      histo->SetBinError(i,histo->GetBinError(i-1));
    }
  }
  histo->GetXaxis()->Set(nbins,&binEdges[0]);

  return histo;
}


//!  \brief Helper method for \p createTwoJetsPtBalanceEventPlots()
//!
//!  A copy from CalibCore, needs to be updated manually
//!  Returns the \p ControlPlotsFunction::Function for a variable
//!  with name \p varName and a \p ControlPlotsConfig::CorrectionType
//!  \p type.
// -------------------------------------------------------------
ControlPlotsFunction::Function BasePlotExtractor::findTwoJetsPtBalanceEventFunction(const std::string& varName, ControlPlotsConfig::CorrectionType type) const {
  if(! TString(names_.at(0)).Contains("MCTruth")){
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
    if( varName == "MET" )
      return &ControlPlotsFunction::twoJetsPtBalanceEventMET;
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
    if( varName == "PF_CH_RespCorrFraction" )
      return &ControlPlotsFunction::twoJetsPtBalanceEventPF_CH_RespCorrFraction;
    if( varName == "PF_NH_RespCorrFraction" )
      return &ControlPlotsFunction::twoJetsPtBalanceEventPF_NH_RespCorrFraction;
    if( varName == "PF_PH_RespCorrFraction" )
      return &ControlPlotsFunction::twoJetsPtBalanceEventPF_PH_RespCorrFraction;
    if( varName == "PF_EL_RespCorrFraction" )
      return &ControlPlotsFunction::twoJetsPtBalanceEventPF_EL_RespCorrFraction;
    if( varName == "PF_HFHad_RespCorrFraction" )
      return &ControlPlotsFunction::twoJetsPtBalanceEventPF_HFHad_RespCorrFraction;
    if( varName == "PF_HFEm_RespCorrFraction" )
      return &ControlPlotsFunction::twoJetsPtBalanceEventPF_HFEm_RespCorrFraction;
    if( varName == "Flavor" )
      return &ControlPlotsFunction::twoJetsPtBalanceEventJetFlavor;
    if( varName == "DeltaPhi" )
      return &ControlPlotsFunction::twoJetsPtBalanceEventDeltaPhi;
    if( varName == "ClosestJetdR" )
      return &ControlPlotsFunction::twoJetsPtBalanceEventClosestJetdR;
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
    if( varName == "GenAsymmetry" )
      return &ControlPlotsFunction::twoJetsPtBalanceEventGenAsymmetry;
    if( varName == "GenLeadAsymmetry" )
      return &ControlPlotsFunction::twoJetsPtBalanceEventGenLeadAsymmetry;  
    if( varName == "GenBarrAsymmetry" )
      return &ControlPlotsFunction::twoJetsPtBalanceEventGenBarrAsymmetry;       
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
    if( varName == "MCTruthResponse" && type == ControlPlotsConfig::Uncorrected )
      return &ControlPlotsFunction::twoJetsPtBalanceEventMCTruthJet1Response;
    if( varName == "MCTruthResponse" && type == ControlPlotsConfig::L2L3 )
      return &ControlPlotsFunction::twoJetsPtBalanceEventMCTruthJet1L2L3Response;
    if( varName == "MCTruthResponse" && type == ControlPlotsConfig::L2L3Res )
      return &ControlPlotsFunction::twoJetsPtBalanceEventMCTruthJet1L2L3Response;
  }

  //add truth functions in addition for compatibility with JEC validation plots
  if( varName == "NPU" )
    return  &ControlPlotsFunction::jetTruthEventNPU;
  if( varName == "Rho" )
    return  &ControlPlotsFunction::jetTruthEventRho;
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
  if( varName == "ClosestJetdR" )
    return &ControlPlotsFunction::jetTruthEventJetClosestJetdR;
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
  if( varName == "GenJetResponse" && type == ControlPlotsConfig::L1L2L3 )
    return  &ControlPlotsFunction::jetTruthEventResponseL1L2L3Corrected;
  if( varName == "GenJetResponse" && type == ControlPlotsConfig::L5 )
    return  &ControlPlotsFunction::jetTruthEventResponseL5Corrected;
  if( varName == "") {
    return 0;
  }
  
  std::cerr << "ControlPlots: unknown variable " << varName << std::endl;
  return 0;
}

#endif 

