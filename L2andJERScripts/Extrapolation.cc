#include "Extrapolation.h"


//! Default constructor, see BasePlotExtractor constructor
//! for more information and initialization.
//!
//!  \author Henning Kirschenmann
//!  \date 2012/03/07
// ----------------------------------------------------------------   
Extrapolation::Extrapolation(TString plotsnames,TString kalibriPlotsShortName) : BasePlotExtractor(plotsnames,kalibriPlotsShortName){
  init(profileType_); //baseplotextracotr method to read in plots
  std::cout << "still initializing..." <<std::endl;
  extrapolInit();
  std::cout << "creating extrapolation plots" <<std::endl;
  createPtRelExtrapol();

}

//! Reads in extra information for extrapolation
//! As of now: read in the cuts on pt3rel to read in the correct
//! histograms
//!
//!  \author Henning Kirschenmann
//!  \date 2012/03/07
// ----------------------------------------------------------------   
void Extrapolation::extrapolInit() {
  cutNames_ = bag_of_string(ExternalConfig_.read<std::string>("TwoJetsPtBalanceEvent plots cut_list",""));
  cutNumbers_ = bag_of<double>(ExternalConfig_.read<std::string>("TwoJetsPtBalanceEvent plots cut_no_list",""));
  if(kalibriPlotsShortName_.Contains("TimePtDependence")){
  cutNames_ = bag_of_string(ExternalConfig_.read<std::string>("TimePtDependence plots cut_list",""));
  cutNumbers_ = bag_of<double>(ExternalConfig_.read<std::string>("TimePtDependence plots cut_no_list",""));
  }
  if(plotsnames_.Contains("VsClosestJetdRPtCut")){//2012PFCHSPtDependence")){
  cutNames_ = bag_of_string(ExternalConfig_.read<std::string>("2012PFCHSPtDependence plots cut_list",""));
  cutNumbers_ = bag_of<double>(ExternalConfig_.read<std::string>("2012PFCHSPtDependence plots cut_no_list",""));
  }
  cutNamesValueToNormalize_ = ExternalConfig_.read<string>((std::string)plotsnames_+" cut_listValueToNormalize",ExternalConfig_.read<string>("Default cut_listValueToNormalize","20"));
  indexToNormalizeTo_ = -1;
  for(unsigned int i=0;i<cutNames_.size();i++){
    std::cout << "cutNamesValue: " << cutNames_.at(i) << " -  " << cutNamesValueToNormalize_ << " - " << indexToNormalizeTo_ << std::endl;
    if(cutNames_.at(i)==cutNamesValueToNormalize_)indexToNormalizeTo_=i;
  }
  assert(indexToNormalizeTo_>-1);
  doPlotExtrapol_ = ExternalConfig_.read<bool>((std::string)plotsnames_+" doPlotExtrapol",ExternalConfig_.read<bool>("Default doPlotExtrapol",1));
  exportOnlyLinearExtrapolation_ = ExternalConfig_.read<bool>((std::string)kalibriPlotsShortName_+" ExportOnlyLinearExtrapolation",ExternalConfig_.read<bool>("Default ExportOnlyLinearExtrapolation",true));

  //  yProfileTitle_ = ExternalConfig_.read<std::string>((std::string)plotsnames_+" plots profile yTitle","DUMMYResolution");
  //  sqrtS_ = ExternalConfig_.read<int>((std::string)kalibriPlotsShortName_+" SqrtS",ExternalConfig_.read<int>("Default SqrtS",7));
  std::cout << "doPlotExtrapol: " << doPlotExtrapol_ <<std::endl;
  if(kalibriPlotsShortName_.Contains("PhiDependence"))doPlotExtrapol_=false;


  plotToRootFileSetup_ = ExternalConfig_.read<string>((std::string)plotsnames_+" plotToRootFileSetup",ExternalConfig_.read<string>("Default plotToRootFileSetup","default"));

  if(plotToRootFileSetup_=="default"){
    saveDeviationPlots_=true;
    saveVsXVariablePlotsAndRatios_=true;
    saveVsBinVarPlots_=true;
    saveExtrapolPlots_=true;

  }
  else if(plotToRootFileSetup_=="MinimalFlavor"){
    saveDeviationPlots_=false;
    saveVsXVariablePlotsAndRatios_=true;
    saveVsBinVarPlots_=false;
    saveExtrapolPlots_=false;
  }
  else std::cerr << "No valid plot to root setup chosen" << std::endl;
     



}



//! Plot and save all the histograms
//!
//!  \author Henning Kirschenmann
//!  \date 2012/03/07
// ----------------------------------------------------------------   
void Extrapolation::Plot() {
  //  MakeDateDir();
  if(chdir("Extrapol") != 0){ 
    mkdir("Extrapol", S_IRWXU|S_IRWXG|S_IRWXO); 
    chdir("Extrapol"); 
  } 



  setTDRStyle();
  TCanvas* c = new TCanvas("c","",600,600);
  c->SetLogx(configs_.at(0)->logX());
  DefaultStyles style;
  //  style.setStyle();

  //  createRatios();
  //  for(int bin_i=0;bin_i<pConfig->nBins();bin_i++){
  for(int bin_i=0;bin_i<configs_.at(0)->nBins();bin_i++){
    int nEntries =2;//names_.size();
    TLegend *leg = new TLegend(0.2,0.9-nEntries*0.07,0.4,0.9);
    leg->SetBorderSize(0);
    leg->SetFillColor(0);
    leg->SetFillStyle(0);
    leg->SetTextFont(42);
    leg->SetTextSize(0.04);
    
    for(int conf_i=0;conf_i<configs_.size();conf_i++){
      for(int j=0;j<2;j++){
	AllPlots_.at(conf_i).at(bin_i).at(j)->SetLineColor(style.getColor(0));
	AllPlots_.at(conf_i).at(bin_i).at(j)->SetMarkerColor(style.getColor(0));
	AllPlots_.at(conf_i).at(bin_i).at(j)->SetMarkerStyle(style.getMarker(0));
      }
    }
//    leg->AddEntry(AllPlots_.at(0).at(bin_i).at(0),(/*configs_.at(0)->yTitle()+*/" Data")/*.c_str()*/,"P");
//    leg->AddEntry(AllPlots_.at(0).at(bin_i).at(1),(/*configs_.at(0)->yTitle()+*/" MC")/*.c_str()*/,"L");
    leg->AddEntry(AllPlots_.at(0).at(bin_i).at(0),dataLabel()/*+" Data"*/,"P");
    leg->AddEntry(AllPlots_.at(0).at(bin_i).at(1),mcLabel()/*+" MC"*/,"L");
    
    

    for(int i=0;i<configs_.size();i++){
      //Draw res in data and MC
      
      AllPlots_.at(i).at(bin_i).at(1)->Draw("hist");
      AllPlots_.at(i).at(bin_i).at(1)->GetYaxis()->SetRangeUser(yProfileMinMax().at(0),yProfileMinMax().at(1));
      AllPlots_.at(i).at(bin_i).at(1)->GetYaxis()->SetTitle(yProfileTitle());
      TStyle *tdrStyle = (TStyle*)gROOT->FindObject("tdrStyle"); 
      AllPlots_.at(i).at(bin_i).at(1)->GetYaxis()->SetTitleOffset(tdrStyle->GetTitleYOffset());
      AllPlots_.at(i).at(bin_i).at(1)->GetXaxis()->SetRangeUser(configs_.at(0)->xMin(),configs_.at(0)->xMax());
      
      AllPlots_.at(i).at(bin_i).at(0)->Draw("p e0 hist same");
      
      AllPlots_.at(i).at(bin_i).at(1)->Draw("hist same");
      drawRunNumberLabels(AllPlots_.at(i).at(bin_i).at(1),configs_.at(i));
      leg->SetHeader((configs_.at(i)->binTitle(configs_.at(i)->binEdges()->at(bin_i),
					       configs_.at(i)->binEdges()->at(bin_i+1))).c_str());
      leg->Draw();
      drawCMSPrel();
      TString outname = "ResolutionPlots_"+plotsnames_+"_"+cutNames_.at(i)+"_"+configs_.at(0)->binName(bin_i);
      c->RedrawAxis();
      c->SaveAs(outname+".pdf");
      if(saveVsXVariablePlotsAndRatios_){
	configs_.at(0)->safelyToRootFile(AllPlots_.at(i).at(bin_i).at(0),outname+"_"+dataLabel()+"_hist");
	configs_.at(0)->safelyToRootFile(AllPlots_.at(i).at(bin_i).at(1),outname+"_"+mcLabel()+"_hist");
      }
     
      //Draw DataMC-ratios
      std::cout << AllRatiosDataMC_.size() << " and " << std::endl;
      AllRatiosDataMC_.at(i).at(bin_i)->Draw();
      drawConfidenceIntervals(AllRatiosDataMC_.at(i).at(bin_i));
      AllRatiosDataMC_.at(i).at(bin_i)->Draw("same");
      AllRatiosDataMC_.at(i).at(bin_i)->GetYaxis()->SetRangeUser(yRatioMinMax().at(0),yRatioMinMax().at(1));
      AllRatiosDataMC_.at(i).at(bin_i)->GetYaxis()->SetTitle(yRatioTitle());
      //    AllRatiosDataMC_.at(i).at(bin_i)->GetYaxis()->SetRangeUser(0.5,1.5);
      //    AllRatiosDataMC_.at(i).at(bin_i)->GetYaxis()->SetTitle("Resolution ratio");
      AllRatiosDataMC_.at(i).at(bin_i)->GetYaxis()->SetTitleOffset(tdrStyle->GetTitleYOffset());
      AllRatiosDataMC_.at(i).at(bin_i)->GetXaxis()->SetRangeUser(configs_.at(0)->xMin(),configs_.at(0)->xMax());
      TPaveText *label = util::LabelFactory::createPaveTextWithOffset(2,1.0,0.5);
      label->AddText(jetLabel_);//+  |#eta_{1,2}| > "+util::toTString(1.4)+",  L = "+util::StyleSettings::luminosity(4.6));
      label->AddText((configs_.at(i)->binTitle(configs_.at(i)->binEdges()->at(bin_i),
					       configs_.at(i)->binEdges()->at(bin_i+1))).c_str());//(configs_.at(0)->BinTitle(xBin_i,configs_.at(0)->binEdges()->at(bin_i),configs_.at(0)->binEdges()->at(bin_i+1),0)).c_str()/*+(TString)" and "*/);
      drawRunNumberLabels( AllRatiosDataMC_.at(i).at(bin_i),configs_.at(i));
      TLegend* leg1 = util::LabelFactory::createLegendWithOffset(2,0.6);
      addFunctionLabelsToLegend( AllRatiosDataMC_.at(i).at(bin_i),leg1);
      //  leg1->AddEntry(MCExtrapols_.at(xBin_i),"Extrapolation (MC)","LP");
      //  leg1->AddEntry(DataExtrapols_.at(xBin_i),"Extrapolation (data)","LP");
      label->Draw("same");
      leg1->Draw();
      drawCMSPrel();
      outname = "ResolutionPlots_"+plotsnames_+"_ratio"+cutNames_.at(i)+"_"+configs_.at(0)->binName(bin_i);
      //  outname+=bin_i;
      //      outname+=".pdf";
      c->RedrawAxis();
      c->SaveAs(outname+".pdf");
      //      if(exportOnlyLinearExtrapolation_){
      if(saveVsXVariablePlotsAndRatios_)configs_.at(0)->safelyToRootFile(AllRatiosDataMC_.at(i).at(bin_i),outname+"_hist");

    }
  }
  
//  //Save RatioVsBinVar plots
//  for(int conf_i=0;conf_i<configs_.size();conf_i++){
//  c->SetLogx(0);
//  //  RatioVsBinVarHistos_.at(conf_i)
//  //    RatioVsBinVarHistos_.at(conf_i)->GetYaxis()->SetRangeUser(0.7,1.3);
//    RatioVsBinVarHistos_.at(conf_i)->GetYaxis()->SetRangeUser(yRatioMinMax().at(0),yRatioMinMax().at(1));
//
//    //    RatioVsBinVarHistos_.at(conf_i)->GetYaxis()->SetTitle("Resolution ratio");
//    TStyle *tdrStyle = (TStyle*)gROOT->FindObject("tdrStyle"); 
//    RatioVsBinVarHistos_.at(conf_i)->GetYaxis()->SetTitleOffset(tdrStyle->GetTitleYOffset());
//    RatioVsBinVarHistos_.at(conf_i)->Draw();
//    RatioVsBinVarHistos_.at(conf_i)->Draw("histsame");
//    drawCMSPrel();
//    TPaveText *label = util::LabelFactory::createPaveTextWithOffset(2,1.0,0.5);
//    label->AddText(jetLabel_);//+  |#eta_{1,2}| > "+util::toTString(1.4)+",  L = "+util::StyleSettings::luminosity(4.6));
//    label->SetFillStyle(0);
//    label->AddText(yProfileTitle()/*plotsnames_*/);
//    label->Draw("same");
//    TString outname = "ResolutionPlots_"+plotsnames_+"_RatioVsBinVar"+cutNames_.at(conf_i)+"_"+names_.at(conf_i);
//    //  outname+=bin_i;
//    RatioVsBinVarHistos_.at(conf_i)->SetName(outname);
//    outname+=".pdf";
//    c->RedrawAxis();
//    c->SaveAs(outname);
//    configs_.at(0)->safelyToRootFile(RatioVsBinVarHistos_.at(conf_i));
//
//
//  }
//


  //Save Ratio/Data/MCVsBinVar plots
  for(int conf_i=0;conf_i<configs_.size();conf_i++){
    c->SetLogx(0);
    TH1vec_t Ratio_Data_MC_histos;
    std::vector <TString> Ratio_Data_MC_labels;
    
    Ratio_Data_MC_histos.push_back(RatioVsBinVarHistos_.at(conf_i));
    Ratio_Data_MC_labels.push_back("Ratio");
    Ratio_Data_MC_histos.push_back(MCVsBinVarHistos_.at(conf_i));
    Ratio_Data_MC_labels.push_back("MC");
    Ratio_Data_MC_histos.push_back(DataVsBinVarHistos_.at(conf_i));
    Ratio_Data_MC_labels.push_back("Data");
    assert(Ratio_Data_MC_histos.size()==Ratio_Data_MC_labels.size());
    
    for(int histo_i=0;histo_i<Ratio_Data_MC_histos.size();histo_i++){
      if(Ratio_Data_MC_labels.at(histo_i).Contains("Ratio"))Ratio_Data_MC_histos.at(histo_i)->GetYaxis()->SetRangeUser(yRatioMinMax().at(0),yRatioMinMax().at(1));
      else Ratio_Data_MC_histos.at(histo_i)->GetYaxis()->SetRangeUser(yProfileMinMax().at(0),yProfileMinMax().at(1));
      TStyle *tdrStyle = (TStyle*)gROOT->FindObject("tdrStyle"); 
      Ratio_Data_MC_histos.at(histo_i)->GetYaxis()->SetTitleOffset(tdrStyle->GetTitleYOffset());
      Ratio_Data_MC_histos.at(histo_i)->Draw();
      Ratio_Data_MC_histos.at(histo_i)->Draw("histsame");
      drawCMSPrel();
      TPaveText *label = util::LabelFactory::createPaveTextWithOffset(2,1.0,0.5);
      label->AddText(jetLabel_);//+  |#eta_{1,2}| > "+util::toTString(1.4)+",  L = "+util::StyleSettings::luminosity(4.6));
      label->SetFillStyle(0);
      label->AddText(yProfileTitle()/*plotsnames_*/);
      label->Draw("same");
      TString outname = "ResolutionPlots_"+plotsnames_+"_"+Ratio_Data_MC_labels.at(histo_i)+"VsBinVar"+cutNames_.at(conf_i)+"_"+names_.at(conf_i);
      //  outname+=bin_i;
      Ratio_Data_MC_histos.at(histo_i)->SetName(outname);
      //      outname+=".pdf";
      c->RedrawAxis();
      c->SaveAs(outname+".pdf");
      if(saveVsBinVarPlots_)configs_.at(0)->safelyToRootFile(Ratio_Data_MC_histos.at(histo_i));
      
    }
  }




  //Save DeviationsOfRatioVsBinVar plot
  for(int conf_i=0;conf_i<configs_.size();conf_i++){
  c->SetLogx(0);
  //  AllDeviationsVsBinVarHistos_.at(conf_i).at(0)->GetYaxis()->SetRangeUser(-0.05,1.2);

    TStyle *tdrStyle = (TStyle*)gROOT->FindObject("tdrStyle"); 
    for(size_t dev_i =0; dev_i<DeviationTypes_.size();dev_i++){
      
      AllDeviationsVsBinVarHistos_.at(conf_i).at(dev_i)->GetYaxis()->SetRangeUser(0,10);
      AllDeviationsVsBinVarHistos_.at(conf_i).at(dev_i)->GetYaxis()->SetTitleOffset(tdrStyle->GetTitleYOffset());
      AllDeviationsVsBinVarHistos_.at(conf_i).at(dev_i)->SetFillColor(kYellow+1);
      AllDeviationsVsBinVarHistos_.at(conf_i).at(dev_i)->SetLineColor(kYellow+1);
      AllDeviationsVsBinVarHistos_.at(conf_i).at(dev_i)->SetMarkerColor(kYellow+1);
      //    AllDeviationsVsBinVarHistos_.at(conf_i).at(dev_i)->Draw();
      AllDeviationsVsBinVarHistos_.at(conf_i).at(dev_i)->Draw("hist");
      drawCMSPrel();
      TPaveText *label = util::LabelFactory::createPaveTextWithOffset(2,1.0,0.05);
      label->AddText(jetLabel_);//+  |#eta_{1,2}| > "+util::toTString(1.4)+",  L = "+util::StyleSettings::luminosity(4.6));
      label->SetFillStyle(0);
      label->AddText("Weighted standard deviation of " + yProfileTitle() + DeviationTypes_.at(dev_i).second);
      label->AddText((TString)"evaluated as function of " + configs_.at(dev_i)->xTitle());
      label->Draw("same");
      TString outname = "ResolutionPlots_"+plotsnames_+"_DeviationsOfRatioVsBinVar"+cutNames_.at(conf_i)+"_"+names_.at(conf_i)+"_"+DeviationTypes_.at(dev_i).first;
      //  outname+=bin_i;
      AllDeviationsVsBinVarHistos_.at(conf_i).at(dev_i)->SetName(outname);
      outname+=".pdf";
      c->RedrawAxis();
      c->SaveAs(outname);
      if(saveDeviationPlots_)configs_.at(0)->safelyToRootFile(AllDeviationsVsBinVarHistos_.at(conf_i).at(dev_i));
    }
  }

  if(DEBUG)std::cout << "configs_.at(0)->nBins()" <<std::endl;
  if(DEBUG)std::cout << "configs_.at(0)->nBins()" << configs_.at(0)->nBins() << std::endl;
  //Save RatioVsBinVar plots for MCDataRatioVsBinVarHistos_
  for(size_t ratio_i=0;ratio_i<MCDataRatioVsBinVarHistos_.size();ratio_i++){
    c->SetLogx(0);
    MCDataRatioVsBinVarHistos_.at(ratio_i)->GetYaxis()->SetRangeUser(yMCDataRatioMinMax().at(0),yMCDataRatioMinMax().at(1));
    if(TString(MCDataRatioVsBinVarHistos_.at(ratio_i)->GetName()).Contains("Normalized"))MCDataRatioVsBinVarHistos_.at(ratio_i)->GetYaxis()->SetRangeUser(0.9, 1.08);
    //    MCDataRatioVsBinVarHistos_.at(ratio_i)->GetYaxis()->SetRangeUser(0.7,1.3);
    TStyle *tdrStyle = (TStyle*)gROOT->FindObject("tdrStyle"); 
    MCDataRatioVsBinVarHistos_.at(ratio_i)->GetYaxis()->SetTitleOffset(tdrStyle->GetTitleYOffset());
    //    MCDataRatioVsBinVarHistos_.at(ratio_i)->Dump();
    MCDataRatioVsBinVarHistos_.at(ratio_i)->Draw();
    MCDataRatioVsBinVarHistos_.at(ratio_i)->Draw("histsame");
    drawCMSPrel();
    TPaveText *label = util::LabelFactory::createPaveTextWithOffset(2,1.0,0.5);
    label->AddText(jetLabel_);//+  |#eta_{1,2}| > "+util::toTString(1.4)+",  L = "+util::StyleSettings::luminosity(4.6));
    label->SetFillStyle(0);
    label->AddText(yProfileTitle()/*plotsnames_*/);
    label->Draw("same");
    TString tempSaveHistoName = MCDataRatioVsBinVarHistos_.at(ratio_i)->GetName();
    TString outname = (TString)"ResolutionPlots_"+plotsnames_+"_RatioVsBinVar_"+MCDataRatioVsBinVarHistos_.at(ratio_i)->GetName();
    MCDataRatioVsBinVarHistos_.at(ratio_i)->SetName(outname);
    outname+=".pdf";
    c->RedrawAxis();
    if(!(TString(MCDataRatioVsBinVarHistos_.at(ratio_i)->GetName()).Contains("Quad")&&exportOnlyLinearExtrapolation_)){
      c->SaveAs(outname);
      if(saveVsBinVarPlots_)configs_.at(0)->safelyToRootFile(MCDataRatioVsBinVarHistos_.at(ratio_i));
    }
    MCDataRatioVsBinVarHistos_.at(ratio_i)->SetName(tempSaveHistoName);

  }

  c->SetLogx(configs_.at(0)->logX());
  for(size_t ratio_i=0;ratio_i<MCDataRatioVsBinVarHistos_.size();ratio_i++){
    for(int bin_i=0;bin_i<configs_.at(0)->nBins();bin_i++){
      TH1D* DRAWHISTO = All_CollectExtrapolatedAllMCDataRatios_.at(ratio_i).at(bin_i);
      DRAWHISTO->GetYaxis()->SetRangeUser(yMCDataRatioMinMax().at(0),yMCDataRatioMinMax().at(1));
      if(TString(MCDataRatioVsBinVarHistos_.at(ratio_i)->GetName()).Contains("Normalized"))DRAWHISTO->GetYaxis()->SetRangeUser(0.9, 1.08);
      //      DRAWHISTO->GetYaxis()->SetRangeUser(0.5,1.1);
      DRAWHISTO->GetYaxis()->SetTitle(MCDataRatioVsBinVarHistos_.at(ratio_i)->GetYaxis()->GetTitle());
      TStyle *tdrStyle = (TStyle*)gROOT->FindObject("tdrStyle"); 
      DRAWHISTO->GetYaxis()->SetTitleOffset(tdrStyle->GetTitleYOffset());
      DRAWHISTO->Draw();
      drawConfidenceIntervals(DRAWHISTO);
      DRAWHISTO->Draw("same");
      drawRunNumberLabels(DRAWHISTO,configs_.at(0));
      drawCMSPrel();
      TPaveText *label = util::LabelFactory::createPaveTextWithOffset(2,1.0,0.5);
      label->AddText(jetLabel_);//+  |#eta_{1,2}| > "+util::toTString(1.4)+",  L = "+util::StyleSettings::luminosity(4.6));
      label->AddText((configs_.at(0)->binTitle(configs_.at(0)->binEdges()->at(bin_i),
					       configs_.at(0)->binEdges()->at(bin_i+1),false)).c_str());
      TLegend* leg1 = util::LabelFactory::createLegendWithOffset(2,0.6);
      addFunctionLabelsToLegend( DRAWHISTO,leg1);
      label->Draw("same");
      leg1->Draw();

      TString outname = (TString)"ResolutionPlots_"+plotsnames_+"_ExtrapolatedMCDataRatios_"+MCDataRatioVsBinVarHistos_.at(ratio_i)->GetName()+"_"+configs_.at(0)->binName(bin_i);
      c->SetName(outname);
      //      outname+=".pdf";
      c->RedrawAxis();
      if(!(TString(MCDataRatioVsBinVarHistos_.at(ratio_i)->GetName()).Contains("Quad")&&exportOnlyLinearExtrapolation_)){
	c->SaveAs(outname+".pdf");
	if(saveVsBinVarPlots_){
	  configs_.at(0)->safelyToRootFile(c);
	  configs_.at(0)->safelyToRootFile(All_CollectExtrapolatedAllMCDataRatios_.at(ratio_i).at(bin_i),outname+"_hist");
	}
      }
    }
  }      

  
  chdir(outputPathROOT()+"/..");
  // chdir("../../."); 
  // configs_.at(0)->closeRootFile();
}

//! Do the radiation extrapolation in each bin (for resolutions mainly bins in eta (and then xbins in pt)
//! Then save the result of the extrapolation and refresh all previous plots (e.g. Data/MC-ratios)
//! 
//!  \author Henning Kirschenmann
//!  \date 2012/03/07
// ----------------------------------------------------------------   
void Extrapolation::createPtRelExtrapol() {
  if(DEBUG)std::cout << "configs_.at(0)->nBins()" <<std::endl;
  if(DEBUG)std::cout << "configs_.at(0)->nBins()" << configs_.at(0)->nBins() << std::endl;
  VecOfTH1vec_t AllExtrapolatedRes;
  //loop over all bins (e.g. eta)
  for(int bin_i=0;bin_i<configs_.at(0)->nBins();bin_i++){
    //create instance of helper class ExtrapolateBin that actually does the extrapolation
    ExtrapolateBin ExtrapolationBin(this);
    if(DEBUG)std::cout << "DEBUG: created ExtrapolateBin " << bin_i << std::endl;
    
    //add all MC and Data histos, i.e. loop over the cuts and then add a MC and data histo for each cut
    for(int conf_i=0;conf_i<configs_.size();conf_i++){
      if(DEBUG)std::cout << "DEBUG: " << conf_i << " of " << configs_.size()<< ". Adding MC and data histo wit X entries MC: " << AllPlots_.at(conf_i).at(bin_i).at(1)->GetEntries()<< " data: " <<AllPlots_.at(conf_i).at(bin_i).at(0)->GetEntries() << std::endl;
      ExtrapolationBin.addMCHisto(AllPlots_.at(conf_i).at(bin_i).at(1));
      ExtrapolationBin.addDataHisto(AllPlots_.at(conf_i).at(bin_i).at(0));

    }
    ExtrapolationBin.calculateAndAddMCDataRatio();
    ExtrapolationBin.calculateAndAddMCDataRatiosNormalizeToSecondCut();
    std::cout << "added all MC and Data histos for current bin" << std::endl;
    // now create extrapolation tgrapherrors for each xbin (e.g. pt)
    // and plot this extrapolation plot
    std::cout << "So many bins... " <<  configs_.at(0)->nXBins() << std::endl;
    for(Int_t xbin_i=0;xbin_i<configs_.at(0)->nXBins();xbin_i++){
      if(DEBUG)std::cout << "TEST "<< xbin_i << " of " << configs_.at(0)->nXBins() <<std::endl;
      ExtrapolationBin.createExtrapolationTGraphErrors(xbin_i+1);
      if(DEBUG)std::cout << "created extrapolation tgrapherrors for bin " << xbin_i<< std::endl;
      if(doPlotExtrapol_)ExtrapolationBin.plotExtrapol(xbin_i,bin_i);
    }
    std::cout << "created extrapolation tgrapherrors for each bin" << std::endl;

    //create extrapolated histograms (e.g. resolution vs pt) for that bin
    ExtrapolationBin.produceExtrapolatedRes();
    //collect those extrapolated histos
    TH1vec_t CollectExtrapolatedRes;
    CollectExtrapolatedRes.push_back(ExtrapolationBin.ExtrapolatedResData());
    CollectExtrapolatedRes.push_back(ExtrapolationBin.ExtrapolatedResMC());
    AllExtrapolatedRes.push_back(CollectExtrapolatedRes);
    CollectExtrapolatedMCDataRatios_.push_back(ExtrapolationBin.ExtrapolatedMCDataRatio());
    CollectExtrapolatedNormalizedMCDataRatios_.push_back(ExtrapolationBin.ExtrapolatedNormalizedMCDataRatio());
    CollectExtrapolatedQuadMCDataRatios_.push_back(ExtrapolationBin.ExtrapolatedQuadMCDataRatio());
    CollectExtrapolatedQuadNormalizedMCDataRatios_.push_back(ExtrapolationBin.ExtrapolatedQuadNormalizedMCDataRatio());
    CollectExtrapolatedLinQuadMCDataRatios_.push_back(ExtrapolationBin.ExtrapolatedLinQuadMCDataRatio());
    CollectExtrapolatedLinQuadNormalizedMCDataRatios_.push_back(ExtrapolationBin.ExtrapolatedLinQuadNormalizedMCDataRatio());
  }
  All_CollectExtrapolatedAllMCDataRatios_.push_back(CollectExtrapolatedMCDataRatios_);
  All_CollectExtrapolatedAllMCDataRatios_.push_back(CollectExtrapolatedNormalizedMCDataRatios_);
  All_CollectExtrapolatedAllMCDataRatios_.push_back(CollectExtrapolatedQuadMCDataRatios_);
  All_CollectExtrapolatedAllMCDataRatios_.push_back(CollectExtrapolatedQuadNormalizedMCDataRatios_);
  All_CollectExtrapolatedAllMCDataRatios_.push_back(CollectExtrapolatedLinQuadMCDataRatios_);
  All_CollectExtrapolatedAllMCDataRatios_.push_back(CollectExtrapolatedLinQuadNormalizedMCDataRatios_);


  // Append these extrapolated histograms as if they were read in right from the start,
  // add the corresponding "virtual" cut of 0.0 to the cutlist and repeat the same "appending"
  // for all other vectors in question
  AllPlots_.push_back(AllExtrapolatedRes);
  cutNumbers_.push_back(0.0);
  cutNames_.push_back("00");
  configs_.push_back(new ControlPlotsConfig(*configs_.at(0)));
  configs_.back()->setCutMax(cutNumbers_.back());
  functions_.push_back(new ControlPlotsFunction(*functions_.at(0)));
  //  profiles_.push_back(new ControlPlotsProfile(*profiles_.at(0)));
  profiles_.push_back(new ControlPlotsProfile(configs_.back(),functions_.back()));
  names_.push_back((((TString)names_.at(0)).ReplaceAll((TString)cutNames_.at(0),(TString)cutNames_.back())/*+"_Extrapol"*/).Data());
  //redo the calculation of data/MC-ratios to propagate the additional extrapolated histograms (defined in BasePlotExtractor)
  refreshRatiosDataMC();
  // fits const and loglin function to dataMC-ratio and produces RatioVsBinVar-plots (e.g. Data/MC-ratio of resolution vs. eta) (defined in BasePlotExtractor)
  makeRatioVsBinVarHistos();



  //MCDataRatios_
  makeMCDataRatioAndNormalizedMCDataRatioVsBinVarHistos();


}


//! Fits const and loglin function to MCData-ratios (added to list of functions of ratio histogram)
//! then creates a RatioVsBinVarHisto using the const fit
//! All these histos are saved to RatioVsBinVarHistos_
//! 
//!  \author Henning Kirschenmann
//!  \date 2012/03/07
// ----------------------------------------------------------------   
void Extrapolation::makeMCDataRatioAndNormalizedMCDataRatioVsBinVarHistos(){
  MCDataRatioVsBinVarHistos_.clear();
  MCDataRatioVsBinVarHistos_.push_back(new TH1D("ExtrapolatedMCDataRatioVsBinVar","",configs_.at(0)->nBins(),&(configs_.at(0)->binEdges()->front())));
  MCDataRatioVsBinVarHistos_.back()->GetYaxis()->SetTitle(""+mcLabel_+"/"+dataLabel_+" ratio (const fit)");
  MCDataRatioVsBinVarHistos_.back()->GetYaxis()->SetRangeUser(yMCDataRatioMinMax().at(0),yMCDataRatioMinMax().at(1));
      //  MCDataRatioVsBinVarHistos_.back()->GetYaxis()->SetRangeUser(0.5,1.1);
  MCDataRatioVsBinVarHistos_.push_back(new TH1D("ExtrapolatedNormalizedMCDataRatioVsBinVar","",configs_.at(0)->nBins(),&(configs_.at(0)->binEdges()->front())));
  MCDataRatioVsBinVarHistos_.back()->GetYaxis()->SetTitle("Radiation correction");
  MCDataRatioVsBinVarHistos_.back()->GetYaxis()->SetRangeUser(0.9,1.08);

  MCDataRatioVsBinVarHistos_.push_back(new TH1D("ExtrapolatedQuadMCDataRatioVsBinVar","",configs_.at(0)->nBins(),&(configs_.at(0)->binEdges()->front())));
  MCDataRatioVsBinVarHistos_.back()->GetYaxis()->SetTitle(""+mcLabel_+"/"+dataLabel_+" ratio (const fit, quadr.)");
  MCDataRatioVsBinVarHistos_.back()->GetYaxis()->SetRangeUser(yMCDataRatioMinMax().at(0),yMCDataRatioMinMax().at(1));
  MCDataRatioVsBinVarHistos_.push_back(new TH1D("ExtrapolatedQuadNormalizedMCDataRatioVsBinVar","",configs_.at(0)->nBins(),&(configs_.at(0)->binEdges()->front())));
  MCDataRatioVsBinVarHistos_.back()->GetYaxis()->SetTitle("Radiation correction (quadr.)");
  MCDataRatioVsBinVarHistos_.back()->GetYaxis()->SetRangeUser(0.9,1.08);


  MCDataRatioVsBinVarHistos_.push_back(new TH1D("ExtrapolatedLinQuadMCDataRatioVsBinVar","",configs_.at(0)->nBins(),&(configs_.at(0)->binEdges()->front())));
  MCDataRatioVsBinVarHistos_.back()->GetYaxis()->SetTitle(""+mcLabel_+"/"+dataLabel_+" ratio (const fit, lin./quadr.)");
  MCDataRatioVsBinVarHistos_.back()->GetYaxis()->SetRangeUser(yMCDataRatioMinMax().at(0),yMCDataRatioMinMax().at(1));
  MCDataRatioVsBinVarHistos_.push_back(new TH1D("ExtrapolatedLinQuadNormalizedMCDataRatioVsBinVar","",configs_.at(0)->nBins(),&(configs_.at(0)->binEdges()->front())));
  MCDataRatioVsBinVarHistos_.back()->GetYaxis()->SetTitle("Radiation correction (lin./quadr.)");
  MCDataRatioVsBinVarHistos_.back()->GetYaxis()->SetRangeUser(0.9,1.08);

  for(unsigned int i=0;i<MCDataRatioVsBinVarHistos_.size();i++){
    MCDataRatioVsBinVarHistos_.at(i)->Sumw2();
    MCDataRatioVsBinVarHistos_.at(i)->GetXaxis()->SetTitle(configs_.at(0)->binAxisTitle().c_str());
    //    for(int bin_i=0;bin_i<configs_.at(0)->nBins()-1;bin_i++){
    for(int bin_i=0;bin_i<configs_.at(0)->nBins();bin_i++){
      //      std::cout <<"TEST"<<std::endl;
      if((TString)MCDataRatioVsBinVarHistos_.at(i)->GetName()=="ExtrapolatedMCDataRatioVsBinVar"){
	fitFunctionsToPlot(CollectExtrapolatedMCDataRatios_.at(bin_i));
	fillRatioVsBinVarPlot(MCDataRatioVsBinVarHistos_.at(i),CollectExtrapolatedMCDataRatios_.at(bin_i),bin_i,"fit_const");
      }
      else if((TString)MCDataRatioVsBinVarHistos_.at(i)->GetName()=="ExtrapolatedNormalizedMCDataRatioVsBinVar"){
	fitFunctionsToPlot(CollectExtrapolatedNormalizedMCDataRatios_.at(bin_i));
	fillRatioVsBinVarPlot(MCDataRatioVsBinVarHistos_.at(i),CollectExtrapolatedNormalizedMCDataRatios_.at(bin_i),bin_i,"fit_const");
      }
      else if((TString)MCDataRatioVsBinVarHistos_.at(i)->GetName()=="ExtrapolatedQuadMCDataRatioVsBinVar"){
	fitFunctionsToPlot(CollectExtrapolatedQuadMCDataRatios_.at(bin_i));
	fillRatioVsBinVarPlot(MCDataRatioVsBinVarHistos_.at(i),CollectExtrapolatedQuadMCDataRatios_.at(bin_i),bin_i,"fit_const");
      }
      else if((TString)MCDataRatioVsBinVarHistos_.at(i)->GetName()=="ExtrapolatedQuadNormalizedMCDataRatioVsBinVar"){
	fitFunctionsToPlot(CollectExtrapolatedQuadNormalizedMCDataRatios_.at(bin_i));
	fillRatioVsBinVarPlot(MCDataRatioVsBinVarHistos_.at(i),CollectExtrapolatedQuadNormalizedMCDataRatios_.at(bin_i),bin_i,"fit_const");
      }
      else if((TString)MCDataRatioVsBinVarHistos_.at(i)->GetName()=="ExtrapolatedLinQuadMCDataRatioVsBinVar"){
	fitFunctionsToPlot(CollectExtrapolatedLinQuadMCDataRatios_.at(bin_i));
	fillRatioVsBinVarPlot(MCDataRatioVsBinVarHistos_.at(i),CollectExtrapolatedLinQuadMCDataRatios_.at(bin_i),bin_i,"fit_const");
      }
      else if((TString)MCDataRatioVsBinVarHistos_.at(i)->GetName()=="ExtrapolatedLinQuadNormalizedMCDataRatioVsBinVar"){
	fitFunctionsToPlot(CollectExtrapolatedLinQuadNormalizedMCDataRatios_.at(bin_i));
	fillRatioVsBinVarPlot(MCDataRatioVsBinVarHistos_.at(i),CollectExtrapolatedLinQuadNormalizedMCDataRatios_.at(bin_i),bin_i,"fit_const");
      }
      //      else 	     std::cout <<"TEST4"<<std::endl;

    }
//    for(int bin_i=0;bin_i<configs_.at(0)->nBins()-1;bin_i++){
//      std::cout << MCDataRatioVsBinVarHistos_.at(i)->GetBinContent(bin_i+1) << ", " << std::endl;
//    }

  }
}




//!  Put JER values into table
//! 
//!  \author Kristin Heine/Henning Kirschenmann
//!  \date 2012/05/23
// ----------------------------------------------------------------   
void Extrapolation::ExportTables() {
  MakeDateDir();
   if(chdir("TXT_tables") != 0){ 
      mkdir("TXT_tables", S_IRWXU|S_IRWXG|S_IRWXO); 
      chdir("TXT_tables"); 
  } 
  if(chdir(kalibriPlotsShortName_) != 0){ 
    mkdir(kalibriPlotsShortName_, S_IRWXU|S_IRWXG|S_IRWXO); 
    chdir(kalibriPlotsShortName_); 
  } 

   for(int conf_i=0;conf_i<configs_.size();conf_i++){
      TString outname = "ResolutionPlots_"+plotsnames_+"_RatioVsBinVar"+cutNames_.at(conf_i)+"_"+names_.at(conf_i);
      outputTable(outname,RatioVsBinVarHistos_.at(conf_i));
   }
   for(size_t ratio_i=0;ratio_i<MCDataRatioVsBinVarHistos_.size();ratio_i++){
     TString outname = (TString)"ResolutionPlots_"+plotsnames_+"_RatioVsBinVar_"+MCDataRatioVsBinVarHistos_.at(ratio_i)->GetName();
     outputTable(outname,MCDataRatioVsBinVarHistos_.at(ratio_i));
   }


  for(int bin_i=0;bin_i<configs_.at(0)->nBins();bin_i++){
    for(int i=0;i<configs_.size();i++){
//      std::cout << "test to setup outname" << std::endl;
//      std::cout <<plotsnames_ << std::endl;
//      std::cout << cutNames_.at(i)<< std::endl;
//      std::cout <<configs_.at(0)->binName(bin_i) << std::endl;
//      std::cout << << std::endl;
//      std::cout << << std::endl;
      TString outname = "ResolutionPlots_"+plotsnames_+"_ratio"+cutNames_.at(i)+"_"+configs_.at(0)->binName(bin_i);
      //      std::cout << "trying to export " << outname << std::endl;
      outputTable(outname,AllRatiosDataMC_.at(i).at(bin_i));
    }
  }
   chdir("../../../."); 
} 


//! Default constructor taking a pointer to the encapsulating object 
//! (in order to access all information saved there)
//! 
//!  \author Henning Kirschenmann
//!  \date 2012/03/07
// ----------------------------------------------------------------   
Extrapolation::ExtrapolateBin::ExtrapolateBin(Extrapolation* Outer) {
  Outer_=Outer;

}

 
 
Extrapolation::ExtrapolateBin::~ExtrapolateBin() {
  for(unsigned int vec_i; vec_i< MCHistos_.size(); vec_i++)delete MCHistos_.at(vec_i);
  for(unsigned int vec_i; vec_i< DataHistos_.size(); vec_i++)delete DataHistos_.at(vec_i);
  for(unsigned int vec_i; vec_i< MCDataRatiosHistos_.size(); vec_i++)delete MCDataRatiosHistos_.at(vec_i);
  for(unsigned int vec_i; vec_i< NormalizedMCDataRatiosHistos_.size(); vec_i++)delete NormalizedMCDataRatiosHistos_.at(vec_i);

//  for(unsigned int vec_i; vec_i< MCExtrapols_.size(); vec_i++)delete MCExtrapols_.at(vec_i);
//  for(unsigned int vec_i; vec_i< DataExtrapols_.size(); vec_i++)delete DataExtrapols_.at(vec_i);
//  for(unsigned int vec_i; vec_i< MCDataRatioExtrapols_.size(); vec_i++)delete MCDataRatioExtrapols_.at(vec_i);
//  for(unsigned int vec_i; vec_i< NormalizedMCDataRatioExtrapols_.size(); vec_i++)delete NormalizedMCDataRatioExtrapols_.at(vec_i);

}


void Extrapolation::ExtrapolateBin::addMCHisto(TH1D* MCHisto){
  MCHistos_.push_back((TH1D*)MCHisto->Clone());
}


void Extrapolation::ExtrapolateBin::addDataHisto(TH1D* DataHisto){
  DataHistos_.push_back((TH1D*)DataHisto->Clone());
}

void Extrapolation::ExtrapolateBin::calculateAndAddMCDataRatio(){

  for(int conf_i=0;conf_i<Outer_->configs_.size();conf_i++){
    TH1D* temp = (TH1D*) (DataHistos_.at(conf_i)->Clone());
    temp->Divide(MCHistos_.at(conf_i),DataHistos_.at(conf_i));
    MCDataRatiosHistos_.push_back(temp);
  }
}

// CONTINUE WITH THIS STUFF...
void Extrapolation::ExtrapolateBin::calculateAndAddMCDataRatiosNormalizeToSecondCut(){
  for(int conf_i=0;conf_i<Outer_->configs_.size();conf_i++){
    TH1D* temp = (TH1D*) (MCDataRatiosHistos_.at(conf_i)->Clone());
//    //hard coded second cut
//    //    jdghjsdhghsjg change here to adapt to different third jet cut
//      //for 0.1,0.2,0.3,0.4
//    //    temp->Divide(MCDataRatiosHistos_.at(conf_i),MCDataRatiosHistos_.at(1));
//   //for fine
    //determine indextonormalize to from config file (using '(std::string)plotsnames_+" cut_listValueToNormalize"')
    temp->Divide(MCDataRatiosHistos_.at(conf_i),MCDataRatiosHistos_.at(Outer_->indexToNormalizeTo_));
    NormalizedMCDataRatiosHistos_.push_back(temp);
  }

}
  //  TH1D* temp = (TH1D*) (DataHistos_.back()->Clone());
//  temp->Divide(MCHistos_.back(),DataHistos_.back());
//
//  NormalizedToCutNo2MCDataRatiosHistos_.push_back
//  TH1D* temp = (TH1D*) (DataHistos_.back()->Clone());
//(DataHistos_.back()->Clone()->Divide(MCHistos_.back(),DataHistos_.back()));
//}

//! create extrapolation tgrapherrors for a xbin
//! 
//! 
//!  \author Henning Kirschenmann
//!  \date 2012/03/07
// ----------------------------------------------------------------   
void Extrapolation::ExtrapolateBin::createExtrapolationTGraphErrors(Int_t xBin_i){

  std::vector<double> x,x_e,MCy,MCy_e,Datay,Datay_e,MCDataRatioy,MCDataRatioy_e,NormalizedMCDataRatioy,NormalizedMCDataRatioy_e;


  //  std::cout << "started to createExtrapolationTGraphError for bin "<<xBin_i <<std::endl;
  //  loop over cuts and fill vectors
  if(DEBUG)std::cout << "DEBUG: " << Outer_->cutNumbers_.size() <<std::endl;

  for (Int_t cut_i =0;cut_i<Outer_->cutNumbers_.size();cut_i++){
    x.push_back(Outer_->cutNumbers_.at(cut_i));
    x_e.push_back(0.);
    MCy.push_back(MCHistos_.at(cut_i)->GetBinContent(xBin_i));
    //    std::cout << "TESTINGMCHistos " << MCHistos_.at(cut_i)->GetBinContent(xBin_i) << std::endl;
    //    std::cout << "TESTINGDataHistos " << DataHistos_.at(cut_i)->GetBinContent(xBin_i) << std::endl;
    MCy_e.push_back(MCHistos_.at(cut_i)->GetBinError(xBin_i));
    Datay.push_back(DataHistos_.at(cut_i)->GetBinContent(xBin_i));
    Datay_e.push_back(DataHistos_.at(cut_i)->GetBinError(xBin_i));
    //  std::cout << "started to createExtrapolationTGraphError for bin "<<xBin_i <<std::endl;

    MCDataRatioy.push_back(MCDataRatiosHistos_.at(cut_i)->GetBinContent(xBin_i));
    MCDataRatioy_e.push_back(MCDataRatiosHistos_.at(cut_i)->GetBinError(xBin_i));
    NormalizedMCDataRatioy.push_back(NormalizedMCDataRatiosHistos_.at(cut_i)->GetBinContent(xBin_i));
    NormalizedMCDataRatioy_e.push_back(NormalizedMCDataRatiosHistos_.at(cut_i)->GetBinError(xBin_i));

    //       std::cout << Outer_->cutNumbers_.at(cut_i) << std::endl;
       //       std::cout<< MCy.back() << " and error: " << MCy_e.back() << std::endl;
  }

  
  //create TGraphErrors from previously defined vectors
  TGraphErrors *extrapol_MC = new TGraphErrors(Outer_->cutNumbers_.size(),&x[0],&MCy[0],&x_e[0],&MCy_e[0]);
  TGraphErrors *extrapol_Data = new TGraphErrors(Outer_->cutNumbers_.size(),&x[0],&Datay[0],&x_e[0],&Datay_e[0]);
  TGraphErrors *extrapol_MCDataRatio = new TGraphErrors(Outer_->cutNumbers_.size(),&x[0],&MCDataRatioy[0],&x_e[0],&MCDataRatioy_e[0]);
  TGraphErrors *extrapol_NormalizedMCDataRatio = new TGraphErrors(Outer_->cutNumbers_.size(),&x[0],&NormalizedMCDataRatioy[0],&x_e[0],&NormalizedMCDataRatioy_e[0]);
  //  TGraphErrors *gr_res1 = new TGraphErrors(n,&x_ptthree_[0],&y_mean_ratio_res1_[0],&ex_ptthree_[0],&ey_mean_ratio_res1_[0]);
  TF1 *lin_extrapol = new TF1("lin_extrapol","[0]+[1]*x",0,Outer_->cutNumbers_.back()+0.05); //was used before...
  TF1 *quad_extrapol = new TF1("quad_extrapol","[0]+[1]*x+[2]*x*x",0,Outer_->cutNumbers_.back()+0.05); //was used before...
  lin_extrapol->SetParameters(1,-0.1);
  lin_extrapol->SetParName(0,"ResZero");
  lin_extrapol->SetParName(1,"slope");
  
  quad_extrapol->SetParameters(100,-0.1,0.1);
  quad_extrapol->SetParName(0,"ResZero");
  quad_extrapol->SetParName(1,"slope");
  quad_extrapol->SetParName(1,"quadraticComponent");
  
  //fit a linear extrapolation function to the TGraphErrors for data and MC
  //and fit a quadratic extrapolation function to the TGraphErrors
  extrapol_MC->Fit("lin_extrapol","Q","same",0,Outer_->cutNumbers_.back()+0.05);
  extrapol_Data->Fit("lin_extrapol","Q","same",0,Outer_->cutNumbers_.back()+0.05);
  extrapol_MCDataRatio->Fit("lin_extrapol","Q","same",0,Outer_->cutNumbers_.back()+0.05);
  extrapol_NormalizedMCDataRatio->Fit("lin_extrapol","Q","same",0,Outer_->cutNumbers_.back()+0.05);

  if(Outer_->exportOnlyLinearExtrapolation_){
    extrapol_MC->GetListOfFunctions()->Add(quad_extrapol);
    extrapol_Data->GetListOfFunctions()->Add(quad_extrapol);
    extrapol_MCDataRatio->GetListOfFunctions()->Add(quad_extrapol);
    extrapol_NormalizedMCDataRatio->GetListOfFunctions()->Add(quad_extrapol);
    //    h->GetListOfFunctions()->Add(func);
  }
  else{
    extrapol_MC->Fit("quad_extrapol","Q+","same",0,Outer_->cutNumbers_.back()+0.05);
    extrapol_Data->Fit("quad_extrapol","Q+","same",0,Outer_->cutNumbers_.back()+0.05);
    extrapol_MCDataRatio->Fit("quad_extrapol","Q+","same",0,Outer_->cutNumbers_.back()+0.05);
    extrapol_NormalizedMCDataRatio->Fit("quad_extrapol","Q+","same",0,Outer_->cutNumbers_.back()+0.05);
  }


  //collect the extrapolation TGraphErrors
  MCExtrapols_.push_back(extrapol_MC);
  DataExtrapols_.push_back(extrapol_Data);
  MCDataRatioExtrapols_.push_back(extrapol_MCDataRatio);
  NormalizedMCDataRatioExtrapols_.push_back(extrapol_NormalizedMCDataRatio);


}

//! Helper function to determine min/max values of the TGraphErrors for extrapolation
//! 
//! 
//!  \author Henning Kirschenmann
//!  \date 2012/03/07
// ----------------------------------------------------------------   
std::pair <float,float> Extrapolation::determineMinMax(TGraphErrors* graph){
  std::pair <float,float> minMaxPair(graph->GetY()[0],graph->GetY()[0]);
  for(Int_t i=0;i<graph->GetN();i++){
    if(graph->GetY()[i]<minMaxPair.first)minMaxPair.first=graph->GetY()[i];
    if(graph->GetY()[i]>minMaxPair.second)minMaxPair.second=graph->GetY()[i];
  }
  return minMaxPair;
}


//! Plot and save the Extrapolation TGraphErrors
//! 
//! 
//!  \author Henning Kirschenmann
//!  \date 2012/03/07
// ----------------------------------------------------------------   
void Extrapolation::ExtrapolateBin::plotExtrapol(Int_t xBin_i, Int_t bin_i){
  //  std::cout << Outer_->cutNumbers_.at(1) << std::endl;
  if(chdir("Resolution") != 0){ 
    mkdir("Resolution", S_IRWXU|S_IRWXG|S_IRWXO); 
    chdir("Resolution"); 
  } 
  
  setTDRStyle();
  DefaultStyles style;
  //  style.setStyle();
  TCanvas* c = new TCanvas("c","",600,600);
    std::pair <float,float> minMaxPair = Outer_->determineMinMax(DataExtrapols_.at(xBin_i));
  //  std::pair <float,float> minMaxPair = std::make_pair(-0.1,0.2);
  c->DrawFrame(0,minMaxPair.first*0.5-0.05,Outer_->cutNumbers_.back()+0.05,minMaxPair.second*1.2,(";cut on "+Outer_->configs_.at(0)->cutAxisTitle()+";"+Outer_->yProfileTitle()/*"#sqrt{2} #sigma"*/)/*.c_str()*/);
  MCExtrapols_.at(xBin_i)->Draw("P");
  Outer_->drawConfidenceIntervals(MCExtrapols_.at(xBin_i));
  Outer_->drawConfidenceIntervals(DataExtrapols_.at(xBin_i));
  MCExtrapols_.at(xBin_i)->Draw("Psame");
  MCExtrapols_.at(xBin_i)->SetLineColor(style.getColor(0));
  MCExtrapols_.at(xBin_i)->SetMarkerColor(style.getColor(0));
  MCExtrapols_.at(xBin_i)->SetMarkerStyle(style.getMarker(0));
  MCExtrapols_.at(xBin_i)->GetFunction("lin_extrapol")->SetLineColor(style.getColor(0));
  MCExtrapols_.at(xBin_i)->GetFunction("lin_extrapol")->SetLineStyle(2);
  TF1* MCTemp=(TF1*) MCExtrapols_.at(xBin_i)->GetFunction("lin_extrapol")->Clone();
  MCTemp->SetRange(0.1,1);
  MCTemp->SetLineStyle(1);
  MCTemp->Draw("same");
  MCExtrapols_.at(xBin_i)->GetFunction("quad_extrapol")->SetLineColor(style.getColor(0));
  MCExtrapols_.at(xBin_i)->GetFunction("quad_extrapol")->SetLineStyle(2);
  TF1* MCQuadTemp=(TF1*) MCExtrapols_.at(xBin_i)->GetFunction("quad_extrapol")->Clone();
  MCQuadTemp->SetRange(0.1,1);
  MCQuadTemp->SetLineStyle(1);
  MCQuadTemp->Draw("same");
  //  TH1F* DrawHist = (TH1F*) MCExtrapols_.at(xBin_i)->GetHistogram();
//  DrawHist->GetXaxis()->SetRangeUser(0.0,Outer_->cutNumbers_.back()+0.05);
//  DrawHist->GetYaxis()->SetRangeUser(0.0,0.15);
  DataExtrapols_.at(xBin_i)->Draw("Psame");
  DataExtrapols_.at(xBin_i)->SetLineColor(style.getColor(1));
  DataExtrapols_.at(xBin_i)->SetMarkerColor(style.getColor(1));
  DataExtrapols_.at(xBin_i)->SetMarkerStyle(style.getMarker(1));
  DataExtrapols_.at(xBin_i)->GetFunction("lin_extrapol")->SetLineColor(style.getColor(1));
  DataExtrapols_.at(xBin_i)->GetFunction("lin_extrapol")->SetLineStyle(2);
  TF1* DataTemp=(TF1*) DataExtrapols_.at(xBin_i)->GetFunction("lin_extrapol")->Clone();
  DataTemp->SetRange(0.1,1);
  DataTemp->SetLineStyle(1);
  DataTemp->Draw("same");
  DataExtrapols_.at(xBin_i)->GetFunction("quad_extrapol")->SetLineColor(style.getColor(0));
  DataExtrapols_.at(xBin_i)->GetFunction("quad_extrapol")->SetLineStyle(2);
  TF1* DataQuadTemp=(TF1*) DataExtrapols_.at(xBin_i)->GetFunction("quad_extrapol")->Clone();
  DataQuadTemp->SetRange(0.1,1);
  DataQuadTemp->SetLineStyle(1);
  DataQuadTemp->Draw("same");
  int nEntries =2;


  TPaveText *label = util::LabelFactory::createPaveTextWithOffset(2,1.0,0.5);
  label->AddText(Outer_->jetLabel_);//+  |#eta_{1,2}| > "+util::toTString(1.4)+",  L = "+util::StyleSettings::luminosity(4.6));
  label->AddText((Outer_->configs_.at(0)->xBinTitle(xBin_i,Outer_->configs_.at(0)->binEdges()->at(bin_i),Outer_->configs_.at(0)->binEdges()->at(bin_i+1),0)).c_str()/*+(TString)" and "*/);
  label->Draw("same");
  TLegend* leg1 = util::LabelFactory::createLegendWithOffset(2,0.6);
 
  leg1->AddEntry(MCExtrapols_.at(xBin_i),"Extrapolation ("+Outer_->mcLabel_+")","LP");
  leg1->AddEntry(DataExtrapols_.at(xBin_i),"Extrapolation ("+Outer_->dataLabel_+")","LP");

  leg1->Draw();
  Outer_->drawCMSPrel();


  TString outname = "ResolutionPlots_ExtrapolPtThree_"+Outer_->plotsnames_+"_"+Outer_->configs_.at(0)->binName(bin_i)+"_"+Outer_->configs_.at(0)->xBinName(xBin_i);
  //  outname+=bin_i;
  outname+=".pdf";
  c->RedrawAxis();
  c->SaveAs(outname);



  //Plot Extrapol of MCDataRatio
  TCanvas* c2 = new TCanvas("c2","",600,600);
  minMaxPair = Outer_->determineMinMax(MCDataRatioExtrapols_.at(xBin_i));
  //  std::pair <float,float> minMaxPair = std::make_pair(-0.1,0.2);
  c2->DrawFrame(0,minMaxPair.first*0.85-0.05,Outer_->cutNumbers_.back()+0.05,minMaxPair.second*1.1,(";cut on "+Outer_->configs_.at(0)->cutAxisTitle()+";"+Outer_->mcLabel_+"/"+Outer_->dataLabel_+" Ratio"/*"#sqrt{2} #sigma"*/));
  MCDataRatioExtrapols_.at(xBin_i)->Draw("P");
  Outer_->drawConfidenceIntervals(MCDataRatioExtrapols_.at(xBin_i));
  MCDataRatioExtrapols_.at(xBin_i)->Draw("Psame");
  MCDataRatioExtrapols_.at(xBin_i)->SetLineColor(style.getColor(0));
  MCDataRatioExtrapols_.at(xBin_i)->SetMarkerColor(style.getColor(0));
  MCDataRatioExtrapols_.at(xBin_i)->SetMarkerStyle(style.getMarker(0));
  MCDataRatioExtrapols_.at(xBin_i)->GetFunction("lin_extrapol")->SetLineColor(style.getColor(0));
  MCDataRatioExtrapols_.at(xBin_i)->GetFunction("lin_extrapol")->SetLineStyle(2);
  TF1* MCDataRatioTemp=(TF1*) MCDataRatioExtrapols_.at(xBin_i)->GetFunction("lin_extrapol")->Clone();
  MCDataRatioTemp->SetRange(0.1,1);
  MCDataRatioTemp->SetLineStyle(1);
  MCDataRatioTemp->Draw("same");
  MCDataRatioExtrapols_.at(xBin_i)->GetFunction("quad_extrapol")->SetLineColor(style.getColor(0));
  MCDataRatioExtrapols_.at(xBin_i)->GetFunction("quad_extrapol")->SetLineStyle(2);
  TF1* MCDataRatioQuadTemp=(TF1*) MCDataRatioExtrapols_.at(xBin_i)->GetFunction("quad_extrapol")->Clone();
  MCDataRatioQuadTemp->SetRange(0.1,1);
  MCDataRatioQuadTemp->SetLineStyle(1);
  MCDataRatioQuadTemp->Draw("same");
  //  TGraphErrors *GraphTemp = (TGraphErrors*)  MCDataRatioExtrapols_.at(xBin_i)->Clone();
  //  Outer_->drawConfidenceIntervals(GraphTemp);
  label = util::LabelFactory::createPaveTextWithOffset(2,1.0,0.5);
  label->AddText(Outer_->jetLabel_);//+  |#eta_{1,2}| > "+util::toTString(1.4)+",  L = "+util::StyleSettings::luminosity(4.6));
  label->AddText((Outer_->configs_.at(0)->xBinTitle(xBin_i,Outer_->configs_.at(0)->binEdges()->at(bin_i),Outer_->configs_.at(0)->binEdges()->at(bin_i+1),0)).c_str()/*+(TString)" and "*/);
  label->Draw("same");
  leg1 = util::LabelFactory::createLegendWithOffset(1,0.6);
  leg1->AddEntry(MCDataRatioExtrapols_.at(xBin_i),"Extrapolation " +Outer_->yProfileTitle(),"LP");
  leg1->Draw();
  Outer_->drawCMSPrel();


  outname = "ResolutionPlots_ExtrapolPtThree_MCDataRatio_"+Outer_->plotsnames_+"_"+Outer_->configs_.at(0)->binName(bin_i)+"_"+Outer_->configs_.at(0)->xBinName(xBin_i);
  //  outname+=bin_i;
  MCDataRatioExtrapols_.at(xBin_i)->SetName(outname);
  outname+=".pdf";
  c2->SaveAs(outname);
  if(Outer_->saveExtrapolPlots_)Outer_->configs_.at(0)->safelyToRootFile(MCDataRatioExtrapols_.at(xBin_i));




  //Plot Extrapol of Normalized MCDataRatio
  TCanvas* c3 = new TCanvas("c3","",600,600);
  minMaxPair = Outer_->determineMinMax(NormalizedMCDataRatioExtrapols_.at(xBin_i));
  //  std::pair <float,float> minMaxPair = std::make_pair(-0.1,0.2);
  c3->DrawFrame(0,minMaxPair.first*0.85-0.05,Outer_->cutNumbers_.back()+0.05,minMaxPair.second*1.1,(";cut on "+Outer_->configs_.at(0)->cutAxisTitle()+";Normalized "+Outer_->mcLabel_+"/"+Outer_->dataLabel_+" Ratio"/*"#sqrt{2} #sigma"*/));
  NormalizedMCDataRatioExtrapols_.at(xBin_i)->Draw("P");
  Outer_->drawConfidenceIntervals(NormalizedMCDataRatioExtrapols_.at(xBin_i));
  NormalizedMCDataRatioExtrapols_.at(xBin_i)->Draw("Psame");
  NormalizedMCDataRatioExtrapols_.at(xBin_i)->SetLineColor(style.getColor(0));
  NormalizedMCDataRatioExtrapols_.at(xBin_i)->SetMarkerColor(style.getColor(0));
  NormalizedMCDataRatioExtrapols_.at(xBin_i)->SetMarkerStyle(style.getMarker(0));
  NormalizedMCDataRatioExtrapols_.at(xBin_i)->GetFunction("lin_extrapol")->SetLineColor(style.getColor(0));
  NormalizedMCDataRatioExtrapols_.at(xBin_i)->GetFunction("lin_extrapol")->SetLineStyle(2);
  TF1* NormalizedMCDataRatioTemp=(TF1*) NormalizedMCDataRatioExtrapols_.at(xBin_i)->GetFunction("lin_extrapol")->Clone();
  NormalizedMCDataRatioTemp->SetRange(0.1,1);
  NormalizedMCDataRatioTemp->SetLineStyle(1);
  NormalizedMCDataRatioTemp->Draw("same");
  NormalizedMCDataRatioExtrapols_.at(xBin_i)->GetFunction("quad_extrapol")->SetLineColor(style.getColor(0));
  NormalizedMCDataRatioExtrapols_.at(xBin_i)->GetFunction("quad_extrapol")->SetLineStyle(2);
  TF1* NormalizedMCDataRatioQuadTemp=(TF1*) NormalizedMCDataRatioExtrapols_.at(xBin_i)->GetFunction("quad_extrapol")->Clone();
  NormalizedMCDataRatioQuadTemp->SetRange(0.1,1);
  NormalizedMCDataRatioQuadTemp->SetLineStyle(1);
  NormalizedMCDataRatioQuadTemp->Draw("same");

  label = util::LabelFactory::createPaveTextWithOffset(2,1.0,0.5);
  label->AddText(Outer_->jetLabel_);//+  |#eta_{1,2}| > "+util::toTString(1.4)+",  L = "+util::StyleSettings::luminosity(4.6));
  label->AddText((Outer_->configs_.at(0)->xBinTitle(xBin_i,Outer_->configs_.at(0)->binEdges()->at(bin_i),Outer_->configs_.at(0)->binEdges()->at(bin_i+1),0)).c_str()/*+(TString)" and "*/);
  label->Draw("same");
  leg1 = util::LabelFactory::createLegendWithOffset(1,0.6);
  leg1->AddEntry(NormalizedMCDataRatioExtrapols_.at(xBin_i),"Extrapolation " +Outer_->yProfileTitle(),"LP");
  leg1->Draw();
  Outer_->drawCMSPrel();


  outname = "ResolutionPlots_ExtrapolPtThree_NormalizedMCDataRatio_"+Outer_->plotsnames_+"_"+Outer_->configs_.at(0)->binName(bin_i)+"_"+Outer_->configs_.at(0)->xBinName(xBin_i);
  //  outname+=bin_i;
  outname+=".pdf";
  c3->SaveAs(outname);




  chdir(".."); 


}

//! Create and save the extrapolated histograms for each bin
//! (e.g. resolutions)
//! 
//!  \author Henning Kirschenmann
//!  \date 2012/03/07
// ----------------------------------------------------------------   
void Extrapolation::ExtrapolateBin::produceExtrapolatedRes(){

  TH1D* tempExtrapolatedResMC = (TH1D*) MCHistos_.at(0)->Clone();
  TH1D* tempExtrapolatedResData = (TH1D*) DataHistos_.at(0)->Clone();
  TH1D* tempExtrapolatedMCDataRatio = (TH1D*) MCDataRatiosHistos_.at(0)->Clone();
  TH1D* tempExtrapolatedNormalizedMCDataRatio = (TH1D*) NormalizedMCDataRatiosHistos_.at(0)->Clone();
  TH1D* tempExtrapolatedQuadMCDataRatio = (TH1D*) MCDataRatiosHistos_.at(0)->Clone();
  TH1D* tempExtrapolatedQuadNormalizedMCDataRatio = (TH1D*) NormalizedMCDataRatiosHistos_.at(0)->Clone();
  TH1D* tempExtrapolatedLinQuadMCDataRatio = (TH1D*) MCDataRatiosHistos_.at(0)->Clone();
  TH1D* tempExtrapolatedLinQuadNormalizedMCDataRatio = (TH1D*) NormalizedMCDataRatiosHistos_.at(0)->Clone();


  for(Int_t xbin_i=0;xbin_i<Outer_->configs_.at(0)->nXBins();xbin_i++){
    tempExtrapolatedResMC->SetBinContent(xbin_i+1,MCExtrapols_.at(xbin_i)->GetFunction("lin_extrapol")->GetParameter(0)>0.01 ? MCExtrapols_.at(xbin_i)->GetFunction("lin_extrapol")->GetParameter(0) : 0.0);
    tempExtrapolatedResData->SetBinContent(xbin_i+1,DataExtrapols_.at(xbin_i)->GetFunction("lin_extrapol")->GetParameter(0)>0.01 ? DataExtrapols_.at(xbin_i)->GetFunction("lin_extrapol")->GetParameter(0) : 0.0);
    tempExtrapolatedMCDataRatio->SetBinContent(xbin_i+1,MCDataRatioExtrapols_.at(xbin_i)->GetFunction("lin_extrapol")->GetParameter(0)>0.01 ? MCDataRatioExtrapols_.at(xbin_i)->GetFunction("lin_extrapol")->GetParameter(0) : 0.0);
    tempExtrapolatedNormalizedMCDataRatio->SetBinContent(xbin_i+1,NormalizedMCDataRatioExtrapols_.at(xbin_i)->GetFunction("lin_extrapol")->GetParameter(0)>0.01 ? NormalizedMCDataRatioExtrapols_.at(xbin_i)->GetFunction("lin_extrapol")->GetParameter(0) : 0.0);
    tempExtrapolatedQuadMCDataRatio->SetBinContent(xbin_i+1,MCDataRatioExtrapols_.at(xbin_i)->GetFunction("quad_extrapol")->GetParameter(0)>0.01 ? MCDataRatioExtrapols_.at(xbin_i)->GetFunction("quad_extrapol")->GetParameter(0) : 0.0);
    tempExtrapolatedQuadNormalizedMCDataRatio->SetBinContent(xbin_i+1,NormalizedMCDataRatioExtrapols_.at(xbin_i)->GetFunction("quad_extrapol")->GetParameter(0)>0.01 ? NormalizedMCDataRatioExtrapols_.at(xbin_i)->GetFunction("quad_extrapol")->GetParameter(0) : 0.0);
    tempExtrapolatedLinQuadMCDataRatio->SetBinContent(xbin_i+1,(tempExtrapolatedMCDataRatio->GetBinContent(xbin_i+1)+tempExtrapolatedQuadMCDataRatio->GetBinContent(xbin_i+1))/2 );
    tempExtrapolatedLinQuadNormalizedMCDataRatio->SetBinContent(xbin_i+1,(tempExtrapolatedNormalizedMCDataRatio->GetBinContent(xbin_i+1)+tempExtrapolatedQuadNormalizedMCDataRatio->GetBinContent(xbin_i+1))/2 );
//    tempExtrapolatedLinQuadMCDataRatio->SetBinContent(xbin_i+1,1.02 );
//    tempExtrapolatedLinQuadNormalizedMCDataRatio->SetBinContent(xbin_i+1,1.01 );

    tempExtrapolatedResMC->SetBinError(xbin_i+1,MCExtrapols_.at(xbin_i)->GetFunction("lin_extrapol")->GetParameter(0)>0.01 ? MCExtrapols_.at(xbin_i)->GetFunction("lin_extrapol")->GetParError(0) : 0.0);
    tempExtrapolatedResData->SetBinError(xbin_i+1,DataExtrapols_.at(xbin_i)->GetFunction("lin_extrapol")->GetParameter(0)>0.01 ? DataExtrapols_.at(xbin_i)->GetFunction("lin_extrapol")->GetParError(0) : 0.0);
    tempExtrapolatedMCDataRatio->SetBinError(xbin_i+1,MCDataRatioExtrapols_.at(xbin_i)->GetFunction("lin_extrapol")->GetParameter(0)>0.01 ? MCDataRatioExtrapols_.at(xbin_i)->GetFunction("lin_extrapol")->GetParError(0) : 0.0);
    tempExtrapolatedNormalizedMCDataRatio->SetBinError(xbin_i+1,NormalizedMCDataRatioExtrapols_.at(xbin_i)->GetFunction("lin_extrapol")->GetParameter(0)>0.01 ? NormalizedMCDataRatioExtrapols_.at(xbin_i)->GetFunction("lin_extrapol")->GetParError(0) : 0.0);
    tempExtrapolatedQuadMCDataRatio->SetBinError(xbin_i+1,MCDataRatioExtrapols_.at(xbin_i)->GetFunction("quad_extrapol")->GetParameter(0)>0.01 ? MCDataRatioExtrapols_.at(xbin_i)->GetFunction("quad_extrapol")->GetParError(0) : 0.0);
    tempExtrapolatedQuadNormalizedMCDataRatio->SetBinError(xbin_i+1,NormalizedMCDataRatioExtrapols_.at(xbin_i)->GetFunction("quad_extrapol")->GetParameter(0)>0.01 ? NormalizedMCDataRatioExtrapols_.at(xbin_i)->GetFunction("quad_extrapol")->GetParError(0) : 0.0);

    double LinQuadDiff = TMath::Abs(tempExtrapolatedMCDataRatio->GetBinContent(xbin_i+1)-tempExtrapolatedQuadMCDataRatio->GetBinContent(xbin_i+1));
    double NormalizedLinQuadDiff = TMath::Abs(tempExtrapolatedNormalizedMCDataRatio->GetBinContent(xbin_i+1)-tempExtrapolatedQuadNormalizedMCDataRatio->GetBinContent(xbin_i+1));
    double LinQuadSYS=LinQuadDiff/2;
    double NormalizedLinQuadSYS= NormalizedLinQuadDiff/2;
    if(tempExtrapolatedMCDataRatio->GetBinError(xbin_i+1)>LinQuadSYS)LinQuadSYS=tempExtrapolatedMCDataRatio->GetBinError(xbin_i+1);
    if(tempExtrapolatedNormalizedMCDataRatio->GetBinError(xbin_i+1)>LinQuadSYS)LinQuadSYS=tempExtrapolatedNormalizedMCDataRatio->GetBinError(xbin_i+1);
    if(tempExtrapolatedQuadMCDataRatio->GetBinError(xbin_i+1)>LinQuadSYS)LinQuadSYS=tempExtrapolatedQuadMCDataRatio->GetBinError(xbin_i+1);
    if(tempExtrapolatedQuadNormalizedMCDataRatio->GetBinError(xbin_i+1)>LinQuadSYS)LinQuadSYS=tempExtrapolatedQuadNormalizedMCDataRatio->GetBinError(xbin_i+1);

    tempExtrapolatedLinQuadMCDataRatio->SetBinError(xbin_i+1,LinQuadDiff);
    tempExtrapolatedLinQuadNormalizedMCDataRatio->SetBinError(xbin_i+1,NormalizedLinQuadDiff);

  }
  ExtrapolatedResMC_=(TH1D*) tempExtrapolatedResMC;
  ExtrapolatedResData_=(TH1D*) tempExtrapolatedResData;
  ExtrapolatedMCDataRatio_=(TH1D*) tempExtrapolatedMCDataRatio;
  ExtrapolatedNormalizedMCDataRatio_=(TH1D*) tempExtrapolatedNormalizedMCDataRatio;
  ExtrapolatedQuadMCDataRatio_              =(TH1D*) tempExtrapolatedQuadMCDataRatio;	       
  ExtrapolatedQuadNormalizedMCDataRatio_    = (TH1D*) tempExtrapolatedQuadNormalizedMCDataRatio;   
  ExtrapolatedLinQuadMCDataRatio_           = (TH1D*) tempExtrapolatedLinQuadMCDataRatio;	       
  ExtrapolatedLinQuadNormalizedMCDataRatio_ =  (TH1D*) tempExtrapolatedLinQuadNormalizedMCDataRatio;

}
