#include "BasePlotExtractor.h"
#include "../scripts/tdrstyle_mod.C"
#include <string>

//Default constructor
BasePlotExtractor::BasePlotExtractor(TString plotsnames,TString kalibriPlotsPath) {
  //init config_file;
  //  config_=ConfigFile("L2L3.cfg");
  //  plotsnames_=plotsnames;
  //  init();
  kalibriPlotsPath_=kalibriPlotsPath;
   config_=ConfigFile("/afs/naf.desy.de/user/k/kirschen/scratch/2012_04_L2L3ResidualsSummary/scripts/L2L3AllPlots.cfg");
  //  config_=ConfigFile("L2L3.cfg");
  plotsnames_=plotsnames;
  ExternalConfig_=ConfigFile("/afs/naf.desy.de/user/k/kirschen/scratch/2012_04_L2L3ResidualsSummary/scripts/ExternalConfigs.cfg");
  readInExtraInfo();
  std::cout << "readInExtraInfo(); executed " << yProfileTitle()<< std::endl;

}

void BasePlotExtractor::readInExtraInfo() {
  TString fileNames ="AK5PF";
  jetType_ = util::LabelFactory::jetAlgo(fileNames)+util::LabelFactory::jetType(fileNames);
  jetLabel_ = util::LabelFactory::labelJet(fileNames);
  std::cout << "import " <<(std::string)plotsnames_+" plots profile yMinMax" << std::endl;
  yProfileMinMax_ = bag_of<double>(ExternalConfig_.read<std::string>((std::string)plotsnames_+" plots profile yMinMax","0 0.5"));
  assert(yProfileMinMax_.size()==2);
  yProfileTitle_ = ExternalConfig_.read<std::string>((std::string)plotsnames_+" plots profile yTitle","DUMMYResolution");
  yRatioMinMax_ = bag_of<double>(ExternalConfig_.read<std::string>((std::string)plotsnames_+" plots ratio yMinMax","0.5 1.5"));
  assert(yRatioMinMax_.size()==2);
  yRatioTitle_ = ExternalConfig_.read<std::string>((std::string)plotsnames_+" plots ratio yTitle","DUMMY Ratio");
  profileType_ = ExternalConfig_.read<std::string>((std::string)plotsnames_+" plots profile type","Mean");


}


void BasePlotExtractor::init(TString profileType) {
  TFile *inf;
  //      inf = new TFile("KalibriPlots_withOldSmear.root");
  //       inf = new TFile("KalibriPlots_noSmear.root");
  //    inf = new TFile("KalibriPlots.root");
  inf = new TFile(kalibriPlotsPath());
  std::vector<std::vector<Event*>* > samples;
  setTDRStyle();

  std::vector<std::string> names;
  names = bag_of_string(ExternalConfig_.read<std::string>((std::string)plotsnames_+" plots names",""));
  std::vector<ControlPlotsConfig*> configs(names.size());
  std::vector<ControlPlotsFunction*> functions(names.size());
  std::vector<ControlPlotsProfile*> profiles(names.size());
  VecOfVecOfTH1vec_t AllPlots;
  VecOfTH1vec_t AllRatiosDataMC;

  // Read different control plot names
  for(size_t i = 0; i < names.size(); i++) {
    std::cout << " Creating plots '" << names.at(i) << "'\n";
    
    // Create ControlPlotsConfig    
    ControlPlotsConfig *pConfig = new ControlPlotsConfig(getConfig(),names.at(i));
    configs.at(i) = pConfig;
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
    
    //    std::cout << pConfig->nBins()<< std::endl;
    VecOfTH1vec_t BinsDMCF;
    TH1vec_t DataMCRatios; //Vector to contain Data/MC-ratios
    for(int bin_i=0;bin_i<pConfig->nBins();bin_i++){
      //      std::cout << bin_i << " of " << pConfig->nBins() << std::endl;
      TH1vec_t DataMCFractions;//Vector to contain Data and MC TH1
      for(ControlPlotsConfig::InputTagsIterator i = pConfig->inputTagsBegin();
	  i != pConfig->inputTagsEnd(); ++i) {
	TString name = pConfig->name();
	name += "/" +pConfig->name();
	name += "_"+pConfig->yVariable()+"Vs"+pConfig->xVariable();
	name += "_"+pConfig->sampleName(i->first);
	name += "_"+pConfig->correctionTypeName(i->second);
	name += "_"+pConfig->binName(bin_i);
	name += "_"+profileType;
	//	std::cout << name << std::endl;
	DataMCFractions.push_back((TH1D*) inf->Get(name));
	std::cout << DataMCFractions.back()->GetName() << std::endl;
      }
      TH1D* ratio = (TH1D*) DataMCFractions.at(0)->Clone();
      ratio->Divide(DataMCFractions.at(0),DataMCFractions.at(1));
      DataMCRatios.push_back(ratio);
      //      DataMCFractions.at(1)->Draw("hist");
      //      DataMCFractions.at(0)->Draw("same");
      //      c->SaveAs(DataMCFractions.back()->GetName()+(TString)".eps");
      BinsDMCF.push_back(DataMCFractions);
    }// (TH1D*)inf->Get("kFSR_vs_Abseta_histo_res1");
    AllPlots.push_back(BinsDMCF);
    AllRatiosDataMC.push_back(DataMCRatios);
  }


  AllRatiosDataMC_ = AllRatiosDataMC;
  AllPlots_        = AllPlots;
  names_	   = names;	
  configs_	   = configs;	
  functions_	   = functions;	
  profiles_        = profiles;    
  
  makeRatioVsBinVarHistos();
}

void BasePlotExtractor::refreshRatiosDataMC(){
  VecOfTH1vec_t AllRatiosDataMC;

  // Read different control plot names
  for(size_t i = 0; i < names_.size(); i++) {
    TH1vec_t DataMCRatios; //Vector to contain Data/MC-ratios
    for(int bin_i=0;bin_i<configs_.at(i)->nBins();bin_i++){
      TH1D* ratio = (TH1D*) AllPlots_.at(i).at(bin_i).at(0)->Clone();
      //      TH1D* ratio = (TH1D*) DataMCFractions.at(0)->Clone();
      ratio->Divide(AllPlots_.at(i).at(bin_i).at(0),AllPlots_.at(i).at(bin_i).at(1));
      //      ratio->Divide(DataMCFractions.at(0),DataMCFractions.at(1));
      DataMCRatios.push_back(ratio);
      //      DataMCFractions.at(1)->Draw("hist");
      //      DataMCFractions.at(0)->Draw("same");
      //      c->SaveAs(DataMCFractions.back()->GetName()+(TString)".eps");
      //      BinsDMCF.push_back(DataMCFractions);
    }// (TH1D*)inf->Get("kFSR_vs_Abseta_histo_res1");
    //    AllPlots.push_back(BinsDMCF);
    AllRatiosDataMC.push_back(DataMCRatios);
  }

  AllRatiosDataMC_ = AllRatiosDataMC;
}

void BasePlotExtractor::makeRatioVsBinVarHistos(){
  TH1vec_t RatioVsBinVarHistos;
  for(int conf_i=0;conf_i<configs_.size();conf_i++){
    std::cout << "done with refreshing ratio-plots" << conf_i<<std::endl;
    std::cout << names_.size() << std::endl;
    std::cout << names_.at(conf_i) << std::endl;
    TH1D* RatioVsBinVar = new TH1D(("RatioVsBinVar"+names_.at(conf_i)).c_str(),"",configs_.at(conf_i)->nBins(),&(configs_.at(conf_i)->binEdges()->front()));
    std::cout << names_.at(conf_i) << std::endl;
    RatioVsBinVar->GetXaxis()->SetTitle(configs_.at(conf_i)->binAxisTitle().c_str());
    RatioVsBinVar->GetYaxis()->SetTitle("Data/MC ratio (const fit)");
  TF1 *fit_const = new TF1("fit_const","[0]",RatioVsBinVar->GetXaxis()->GetXmin(),RatioVsBinVar->GetXaxis()->GetXmax()); //was used before...
  fit_const->SetParameters(1,1);
  fit_const->SetParName(0,"const");
  
  TF1 *fit_loglin = new TF1("fit_loglin","[0]+[1]*TMath::Log(x)",RatioVsBinVar->GetXaxis()->GetXmin(),RatioVsBinVar->GetXaxis()->GetXmax()); //was used before...
  fit_loglin->SetParameters(1,1);
  fit_loglin->SetParName(0,"const");
  fit_loglin->SetParName(1,"slope");

     std::cout << names_.at(conf_i) << "test4" << std::endl;
      //    hXSpectrum_[id] = new TH1D(name.c_str(),"",
      //			       config_->nXBins(),&(config_->xBinEdges()->front()));

    for(int bin_i=0;bin_i<configs_.at(0)->nBins()-1;bin_i++){
      //      std::cout << names_.at(conf_i) << "test5"<< std::endl;
      AllRatiosDataMC_.at(conf_i).at(bin_i)->Fit("fit_const","EMQ","same");
      AllRatiosDataMC_.at(conf_i).at(bin_i)->Fit("fit_loglin","EMQ+","same");
      //      std::cout << names_.at(conf_i) << "test5"<< std::endl;
      if(AllRatiosDataMC_.at(conf_i).at(bin_i)->GetFunction("fit_const")){
	//      TList *lfits = AllRatiosDataMC_.at(conf_i).at(bin_i)->GetListOfFunctions();
	//      lfits->ls();
	//      std::cout << AllRatiosDataMC_.at(conf_i).at(bin_i)->GetFunction("fit_const")->GetParameter(0) << std::endl;
	RatioVsBinVar->SetBinContent(bin_i+1,AllRatiosDataMC_.at(conf_i).at(bin_i)->GetFunction("fit_const")->GetParameter(0));
	//      std::cout << names_.at(conf_i) << "test5"<< std::endl;
	RatioVsBinVar->SetBinError(bin_i+1,AllRatiosDataMC_.at(conf_i).at(bin_i)->GetFunction("fit_const")->GetParError(0));

      }
      else {
	std::cout << names_.at(conf_i) << "WARNING: NO FIT FUNCTIONS ADDED... SETTING TO ZERO BIN: " << bin_i<< std::endl;
	RatioVsBinVar->SetBinContent(bin_i+1,0.);
	RatioVsBinVar->SetBinError(bin_i+1,0.);
      }
      std::cout << names_.at(conf_i) << "test6"<< std::endl;
      
    }
        for(int bin_i=0;bin_i<configs_.at(0)->nBins()-1;bin_i++){
	  std::cout << RatioVsBinVar->GetBinContent(bin_i+1) << ", " << std::endl;
    }

    RatioVsBinVarHistos.push_back(RatioVsBinVar);

    std::cout << "done with refreshing ratio-plots" << conf_i<<std::endl;
  }

  
  RatioVsBinVarHistos_=RatioVsBinVarHistos;
}

void BasePlotExtractor::addFunctionLabelsToLegend(TH1D* histo, TLegend* leg){
  DefaultStyles style;
  style.setStyle("Confidence");
  TIter next(histo->GetListOfFunctions());
  TObject *obj;
  Int_t i=0;
  while ((obj = next())){
    //    obj->Print();//Draw(next.GetOption());
    TF1* temp_tf1 = (TF1*) obj;
    temp_tf1->SetLineColor(style.getColor(i));
    temp_tf1->SetFillColor(style.getColor(i));

    TFitResultPtr r = histo->Fit(temp_tf1->GetName(),"NSEMQ");
    Double_t chisquared1;
    unsigned int ndf1;
    char buffer [50];

//    import_kFSR_vs_Abseta_histo_res1_Dir1->Fit(kFSR_fit,"N 0  S E");
//    TMatrixDSym cov = r->GetCovarianceMatrix();  //  to access the covariance matrix
    chisquared1= r->Chi2();
    ndf1= r->Ndf();


    //    histo->Fit(temp_tf1->GetName(),"NS");
    sprintf (buffer, " (#chi^{2}/ndf = %.1f / %u )", chisquared1, ndf1);
    std::cout << buffer << std::endl;

    std::string legendName = ExternalConfig_.read<std::string>(temp_tf1->GetName()+(std::string) " legend label","DUMMYlabel");
    //      yProfileTitle_ = ExternalConfig_.read<std::string>((std::string)plotsnames_+" plots profile yTitle","DUMMYResolution");

    leg->AddEntry(temp_tf1,(legendName+buffer).c_str(),"F");
    i++;
  }

}


void BasePlotExtractor::drawConfidenceIntervals(TH1D* histo){
  DefaultStyles style;
  style.setStyle("Confidence");
  
  std::cout << "List of functions: " << std::endl;
  //	   histo->GetListOfFunctions()->Print();
  TIter next(histo->GetListOfFunctions());
  TObject *obj;
  Int_t i=0;
  while ((obj = next())){
    //    obj->Print();//Draw(next.GetOption());
    TF1* temp_tf1 = (TF1*) obj;

    histo->Fit(temp_tf1->GetName(),"NEMQ");
    //Create a histogram to hold the confidence intervals
    TH1D *hint = new TH1D("hint","Fitted gaussian with .95 conf.band", 100, histo->GetXaxis()->GetXmin(), histo->GetXaxis()->GetXmax());
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
    hint->Draw("e4 same");
    //	     hint->Print("all"); 
  }
  
}

//!  \brief Helper method for \p createTwoJetsPtBalanceEventPlots()
//!
//!  Returns the \p ControlPlotsFunction::Function for a variable
//!  with name \p varName and a \p ControlPlotsConfig::CorrectionType
//!  \p type.
// -------------------------------------------------------------
ControlPlotsFunction::Function BasePlotExtractor::findTwoJetsPtBalanceEventFunction(const std::string& varName, ControlPlotsConfig::CorrectionType type) const {
  if( varName == "Eta" )
    return  &ControlPlotsFunction::twoJetsPtBalanceEventJetEta;
  if( varName == "AbsEta" )
    return  &ControlPlotsFunction::twoJetsPtBalanceEventJetAbsEta;
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
  if( varName == "PF_CH_Fraction" )
    return &ControlPlotsFunction::twoJetsPtBalanceEventPF_CH_Fraction;
  if( varName == "PF_NH_Fraction" )
    return &ControlPlotsFunction::twoJetsPtBalanceEventPF_NH_Fraction;
  if( varName == "PF_PH_Fraction" )
    return &ControlPlotsFunction::twoJetsPtBalanceEventPF_PH_Fraction;
  if( varName == "PF_EL_Fraction" )
    return &ControlPlotsFunction::twoJetsPtBalanceEventPF_EL_Fraction;
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
