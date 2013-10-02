#include "Extrapolation.h"
#include "TInterpreter.h"
#include "Rtypes.h"
#include "TMatrixD.h"

#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <sstream>
#include <iomanip>

// any data to be used by chi2_linear for evaluationg the function:
struct fit_data{
    // x, y values:
   std::vector<double> x_val, y_val;
    // variance-covariance matrix for y values:
    TMatrixD y_cov;
    // inverted cov matrix; calculated by chi2_linear "on demand".
    TMatrixD y_cov_inv;
    
    void reset(){
        x_val.clear();
        y_val.clear();
        y_cov.ResizeTo(0,0);
        y_cov_inv.ResizeTo(0,0);
    }
};

fit_data data;

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
  doPlotExtrapol_ = ExternalConfig_.read<bool>((std::string)plotsnames_+" doPlotExtrapol",ExternalConfig_.read<bool>("Default doPlotExtrapol",1));
  //  yProfileTitle_ = ExternalConfig_.read<std::string>((std::string)plotsnames_+" plots profile yTitle","DUMMYResolution");
  //  sqrtS_ = ExternalConfig_.read<int>((std::string)kalibriPlotsShortName_+" SqrtS",ExternalConfig_.read<int>("Default SqrtS",7));
  std::cout << "doPlotExtrapol: " << doPlotExtrapol_ <<std::endl;
  if(kalibriPlotsShortName_.Contains("PhiDependence"))doPlotExtrapol_=false;
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

  std::cout << "creating extrapolation plots" <<std::endl;
  createPtRelExtrapol();


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
    leg->AddEntry(AllPlots_.at(0).at(bin_i).at(0),(/*configs_.at(0)->yTitle()+*/" Data")/*.c_str()*/,"P");
    leg->AddEntry(AllPlots_.at(0).at(bin_i).at(1),(/*configs_.at(0)->yTitle()+*/" MC")/*.c_str()*/,"L");
    
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
      AllPlots_.at(i).at(bin_i).at(0)->SetName(outname + "_Data");
      AllPlots_.at(i).at(bin_i).at(1)->SetName(outname + "_MC");
      //  outname+=bin_i;
      outname+=".eps";
      c->RedrawAxis();
      c->SaveAs(outname);
      configs_.at(0)->safelyToRootFile(AllPlots_.at(i).at(bin_i).at(0));
      configs_.at(0)->safelyToRootFile(AllPlots_.at(i).at(bin_i).at(1));

      
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
      outname+=".eps";
      c->RedrawAxis();
      c->SaveAs(outname);
    }
  }
  
  //Save RatioVsBinVar plots
  for(int conf_i=0;conf_i<configs_.size();conf_i++){
  c->SetLogx(0);
  //  RatioVsBinVarHistos_.at(conf_i)

  // RatioVsBinVarHistos_.at(conf_i)->GetYaxis()->SetRangeUser(0.4,1.6); //0.7 1.3
  //    RatioVsBinVarHistos_.at(conf_i)->GetYaxis()->SetRangeUser(0.7,1.3);
    RatioVsBinVarHistos_.at(conf_i)->GetYaxis()->SetRangeUser(yRatioMinMax().at(0),yRatioMinMax().at(1));

    //    RatioVsBinVarHistos_.at(conf_i)->GetYaxis()->SetTitle("Resolution ratio");
    TStyle *tdrStyle = (TStyle*)gROOT->FindObject("tdrStyle"); 
    RatioVsBinVarHistos_.at(conf_i)->GetYaxis()->SetTitleOffset(tdrStyle->GetTitleYOffset());
    RatioVsBinVarHistos_.at(conf_i)->Draw();
    RatioVsBinVarHistos_.at(conf_i)->Draw("histsame");
    drawCMSPrel();
    TPaveText *label = util::LabelFactory::createPaveTextWithOffset(2,1.0,0.5);
    label->AddText(jetLabel_);//+  |#eta_{1,2}| > "+util::toTString(1.4)+",  L = "+util::StyleSettings::luminosity(4.6));
    label->SetFillStyle(0);
    label->AddText(yProfileTitle()/*plotsnames_*/);
    label->Draw("same");
    TString outname = "ResolutionPlots_"+plotsnames_+"_RatioVsBinVar"+cutNames_.at(conf_i)+"_"+names_.at(conf_i);
    //  outname+=bin_i;
    RatioVsBinVarHistos_.at(conf_i)->SetName(outname);
    outname+=".eps";
    c->RedrawAxis();
    c->SaveAs(outname);
    configs_.at(0)->safelyToRootFile(RatioVsBinVarHistos_.at(conf_i));


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
      outname+=".eps";
      c->RedrawAxis();
      c->SaveAs(outname);
      configs_.at(0)->safelyToRootFile(AllDeviationsVsBinVarHistos_.at(conf_i).at(dev_i));
    }
  }



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
    outname+=".eps";
    c->RedrawAxis();
    c->SaveAs(outname);
    configs_.at(0)->safelyToRootFile(MCDataRatioVsBinVarHistos_.at(ratio_i));
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
      outname+=".eps";
      c->RedrawAxis();
      c->SaveAs(outname);
      configs_.at(0)->safelyToRootFile(c);
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
      ExtrapolationBin.createExtrapolationTGraphErrors(xbin_i+1, bin_i);
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
  }
  All_CollectExtrapolatedAllMCDataRatios_.push_back(CollectExtrapolatedMCDataRatios_);
  All_CollectExtrapolatedAllMCDataRatios_.push_back(CollectExtrapolatedNormalizedMCDataRatios_);


  // Append these extrapolated histograms as if they were read in right from the start,
  // add the corresponding "virtual" cut of 0.0 to the cutlist and repeat the same "appending"
  // for all other vectors in question
  AllPlots_.push_back(AllExtrapolatedRes);
  cutNumbers_.push_back(0.0);
  cutNames_.push_back("00");
  configs_.push_back(new ControlPlotsConfig(*configs_.at(0)));
  configs_.back()->setCutMax(cutNumbers_.back());
  functions_.push_back(new ControlPlotsFunction(*functions_.at(0)));
  profiles_.push_back(new ControlPlotsProfile(*profiles_.at(0)));
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
  MCDataRatioVsBinVarHistos_.back()->GetYaxis()->SetTitle("MC/Data ratio (const fit)");
  MCDataRatioVsBinVarHistos_.back()->GetYaxis()->SetRangeUser(yMCDataRatioMinMax().at(0),yMCDataRatioMinMax().at(1));
      //  MCDataRatioVsBinVarHistos_.back()->GetYaxis()->SetRangeUser(0.5,1.1);
  MCDataRatioVsBinVarHistos_.push_back(new TH1D("ExtrapolatedNormalizedMCDataRatioVsBinVar","",configs_.at(0)->nBins(),&(configs_.at(0)->binEdges()->front())));
  MCDataRatioVsBinVarHistos_.back()->GetYaxis()->SetTitle("Radiation correction");
  MCDataRatioVsBinVarHistos_.back()->GetYaxis()->SetRangeUser(0.9,1.08);

  for(unsigned int i=0;i<MCDataRatioVsBinVarHistos_.size();i++){
    MCDataRatioVsBinVarHistos_.at(i)->Sumw2();
    MCDataRatioVsBinVarHistos_.at(i)->GetXaxis()->SetTitle(configs_.at(0)->binAxisTitle().c_str());
    //    for(int bin_i=0;bin_i<configs_.at(0)->nBins()-1;bin_i++){
    for(int bin_i=0;bin_i<configs_.at(0)->nBins();bin_i++){
      //      std::cout <<"TEST"<<std::endl;
      if((TString)MCDataRatioVsBinVarHistos_.at(i)->GetName()=="ExtrapolatedMCDataRatioVsBinVar"){
	//	     std::cout <<"TEST2"<<std::endl;
	fitFunctionsToPlot(CollectExtrapolatedMCDataRatios_.at(bin_i));
	fillRatioVsBinVarPlot(MCDataRatioVsBinVarHistos_.at(i),CollectExtrapolatedMCDataRatios_.at(bin_i),bin_i,"fit_const");
      }
      else if((TString)MCDataRatioVsBinVarHistos_.at(i)->GetName()=="ExtrapolatedNormalizedMCDataRatioVsBinVar"){
	//	     std::cout <<"TEST3"<<std::endl;
	fitFunctionsToPlot(CollectExtrapolatedNormalizedMCDataRatios_.at(bin_i));
	fillRatioVsBinVarPlot(MCDataRatioVsBinVarHistos_.at(i),CollectExtrapolatedNormalizedMCDataRatios_.at(bin_i),bin_i,"fit_const");
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


void Extrapolation::ExtrapolateBin::addMCHisto(TH1D* MCHisto){
  MCHistos_.push_back(MCHisto);
}


void Extrapolation::ExtrapolateBin::addDataHisto(TH1D* DataHisto){
  DataHistos_.push_back(DataHisto);
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
    //hard coded second cut
    temp->Divide(MCDataRatiosHistos_.at(conf_i),MCDataRatiosHistos_.at(1));
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

// the chi^2 to minimize for fitting a linear function
//   y = p[0]*x + p[1]
// with fit parameters p[0], p[1] to data with known x and y and covariance
// matrix for y.
void chi2_linear(Int_t& npar, Double_t* grad, Double_t& fval, Double_t* p, Int_t status){
    if(data.y_cov_inv.GetNcols()==0){
        double dummy;
        int ncols = data.y_cov.GetNcols();
        data.y_cov_inv.ResizeTo(ncols, ncols);
        data.y_cov_inv = data.y_cov.Invert(&dummy);
    }
    const size_t ndata = data.x_val.size(); // number of data points in x,y graph to fit to
    std::vector<double> delta_y(ndata);
    for(size_t i=0; i<ndata; ++i){
        delta_y[i] = data.x_val[i]*p[0] + p[1] - data.y_val[i];
    }
    // now calculate the chi2, i.e.
    //  dy^T * C^{-1} * dy
    // where C is the variance--covariance matrix and dy = (y_data - y_pred)
    // This could probably be implemented in ROOT, but it's so simple, we just do it here:
    fval = 0.0;
    for(size_t i=0; i<ndata; ++i){
        for(size_t j=0; j<ndata; ++j){
            fval += delta_y[i] * delta_y[j] * data.y_cov_inv(i,j);
        }
    }
}

void make_lin_fit(double & slope, double & d_slope, double & offset, double & d_offset){
    TMinuit min;
    min.SetPrintLevel(-1);
    int err = min.DefineParameter(0, "slope", slope, d_slope, -10.0, 10.0);
    assert(err==0);
    err = min.DefineParameter(1, "offset", offset, d_offset, 0.0, 1e3);
    assert(err==0);
    min.SetFCN(chi2_linear);
    min.mnmigr();
    min.GetParameter(0, slope, d_slope);
    min.GetParameter(1, offset, d_offset);
}

//! create extrapolation tgrapherrors for a xbin
//! 
//! 
//!  \author Henning Kirschenmann
//!  \date 2012/03/07
// ----------------------------------------------------------------   
void Extrapolation::ExtrapolateBin::createExtrapolationTGraphErrors(Int_t xBin_i, Int_t bin_i){

   doLinExtrapol_ = true;

   // std::vector<double> x,x_e,MCy,MCy_e,Datay,Datay_e,MCDataRatioy,MCDataRatioy_e,NormalizedMCDataRatioy,NormalizedMCDataRatioy_e;
   std::vector<double> x,x_e,x_mc,x_e_mc,x_data,x_e_data,MCy,MCy_e,Datay,Datay_e,MCDataRatioy,MCDataRatioy_e,NormalizedMCDataRatioy,NormalizedMCDataRatioy_e;

  //  std::cout << "started to createExtrapolationTGraphError for bin "<<xBin_i <<std::endl;
  //  loop over cuts and fill vectors
  if(DEBUG)std::cout << "DEBUG: " << Outer_->cutNumbers_.size() <<std::endl;

  // Covariance matrices needed for fitting 
  TMatrixD y_cov_mc;
  y_cov_mc.ResizeTo(Outer_->cutNumbers_.size(), Outer_->cutNumbers_.size());
  TMatrixD y_cov_data;
  y_cov_data.ResizeTo(Outer_->cutNumbers_.size(), Outer_->cutNumbers_.size());

  std::cout << "Bin: " << Outer_->configs_.at(0)->xBinTitle(xBin_i-1,Outer_->configs_.at(0)->binEdges()->at(bin_i),Outer_->configs_.at(0)->binEdges()->at(bin_i+1),0) << std::endl;

  for (Int_t cut_i =0;cut_i<Outer_->cutNumbers_.size();cut_i++){
     if( doLinExtrapol_ ) {
        x.push_back(Outer_->cutNumbers_.at(cut_i));
        x_e.push_back(0.);
        x_mc.push_back(Outer_->cutNumbers_.at(cut_i));
        x_e_mc.push_back(0.);
        x_data.push_back(Outer_->cutNumbers_.at(cut_i));
        x_e_data.push_back(0.);
     }
     else {
        x.push_back(Outer_->cutNumbers_.at(cut_i));
        x_e.push_back(0.);
        x_mc.push_back(getMeanAlphaMC(xBin_i, bin_i).at(cut_i));
        x_e_mc.push_back(getMeanAlphaMCErr(xBin_i, bin_i).at(cut_i));
        x_data.push_back(getMeanAlphaData(xBin_i, bin_i).at(cut_i));
        x_e_data.push_back(getMeanAlphaDataErr(xBin_i, bin_i).at(cut_i));
     }
    MCy.push_back(MCHistos_.at(cut_i)->GetBinContent(xBin_i));
    //    std::cout << "TESTINGMCHistos " << MCHistos_.at(cut_i)->GetBinContent(xBin_i) << std::endl;
    //    std::cout << "TESTINGDataHistos " << DataHistos_.at(cut_i)->GetBinContent(xBin_i) << std::endl;
    MCy_e.push_back(MCHistos_.at(cut_i)->GetBinError(xBin_i));
    Datay.push_back(DataHistos_.at(cut_i)->GetBinContent(xBin_i));
    Datay_e.push_back(DataHistos_.at(cut_i)->GetBinError(xBin_i));

    // fill covariance matrix for data and mc
    for (Int_t cut_j =0;cut_j<Outer_->cutNumbers_.size();cut_j++){
       if( cut_i <= cut_j ) {
          double n1_mc = pow(MCHistos_.at(cut_i)->GetBinContent(xBin_i),2)/(2*pow(MCHistos_.at(cut_i)->GetBinError(xBin_i),2));
          double n2_mc = pow(MCHistos_.at(cut_j)->GetBinContent(xBin_i),2)/(2*pow(MCHistos_.at(cut_j)->GetBinError(xBin_i),2));

          double n1_data = pow(DataHistos_.at(cut_i)->GetBinContent(xBin_i),2)/(2*pow(DataHistos_.at(cut_i)->GetBinError(xBin_i),2));
          double n2_data = pow(DataHistos_.at(cut_j)->GetBinContent(xBin_i),2)/(2*pow(DataHistos_.at(cut_j)->GetBinError(xBin_i),2));

          y_cov_mc(cut_i, cut_j) = pow(MCHistos_.at(cut_i)->GetBinError(xBin_i),2)*
             pow((n1_mc/n2_mc),2)*
             (MCHistos_.at(cut_i)->GetBinContent(xBin_i)/MCHistos_.at(cut_j)->GetBinContent(xBin_i));
          y_cov_data(cut_i, cut_j) = pow(DataHistos_.at(cut_i)->GetBinError(xBin_i),2)*
             pow((n1_data/n2_data),2)*
             (DataHistos_.at(cut_i)->GetBinContent(xBin_i)/DataHistos_.at(cut_j)->GetBinContent(xBin_i));
       }
       else {
          double n1_mc = pow(MCHistos_.at(cut_j)->GetBinContent(xBin_i),2)/(2*pow(MCHistos_.at(cut_j)->GetBinError(xBin_i),2));
          double n2_mc = pow(MCHistos_.at(cut_i)->GetBinContent(xBin_i),2)/(2*pow(MCHistos_.at(cut_i)->GetBinError(xBin_i),2));

          double n1_data = pow(DataHistos_.at(cut_j)->GetBinContent(xBin_i),2)/(2*pow(DataHistos_.at(cut_j)->GetBinError(xBin_i),2));
          double n2_data = pow(DataHistos_.at(cut_i)->GetBinContent(xBin_i),2)/(2*pow(DataHistos_.at(cut_i)->GetBinError(xBin_i),2));

          y_cov_mc(cut_i, cut_j) = pow(MCHistos_.at(cut_j)->GetBinError(xBin_i),2)*
             pow((n1_mc/n2_mc),2)*
             (MCHistos_.at(cut_j)->GetBinContent(xBin_i)/MCHistos_.at(cut_i)->GetBinContent(xBin_i));
          y_cov_data(cut_i, cut_j) = pow(DataHistos_.at(cut_j)->GetBinError(xBin_i),2)*
             pow((n1_data/n2_data),2)*
             (DataHistos_.at(cut_j)->GetBinContent(xBin_i)/DataHistos_.at(cut_i)->GetBinContent(xBin_i));
       }
    }
  
       //  std::cout << "started to createExtrapolationTGraphError for bin "<<xBin_i <<std::endl;

    MCDataRatioy.push_back(MCDataRatiosHistos_.at(cut_i)->GetBinContent(xBin_i));
    MCDataRatioy_e.push_back(MCDataRatiosHistos_.at(cut_i)->GetBinError(xBin_i));
    NormalizedMCDataRatioy.push_back(NormalizedMCDataRatiosHistos_.at(cut_i)->GetBinContent(xBin_i));
    NormalizedMCDataRatioy_e.push_back(NormalizedMCDataRatiosHistos_.at(cut_i)->GetBinError(xBin_i));

    std::cout << Outer_->cutNumbers_.at(cut_i) << std::endl;
    std::cout<< "MC y values: " << MCy.back() << " and error: " << MCy_e.back() << std::endl;
    std::cout<< "Data y values: " << Datay.back() << " and error: " << Datay_e.back() << std::endl;
  }
  
  //create TGraphErrors from previously defined vectors
  TGraphErrors *extrapol_MC = new TGraphErrors(Outer_->cutNumbers_.size(),&x_mc[0],&MCy[0],&x_e_mc[0],&MCy_e[0]);
  TGraphErrors *extrapol_Data = new TGraphErrors(Outer_->cutNumbers_.size(),&x_data[0],&Datay[0],&x_e_data[0],&Datay_e[0]);
  TGraphErrors *extrapol_MCDataRatio = new TGraphErrors(Outer_->cutNumbers_.size(),&x[0],&MCDataRatioy[0],&x_e[0],&MCDataRatioy_e[0]);
  TGraphErrors *extrapol_NormalizedMCDataRatio = new TGraphErrors(Outer_->cutNumbers_.size(),&x[0],&NormalizedMCDataRatioy[0],&x_e[0],&NormalizedMCDataRatioy_e[0]);
  //  TGraphErrors *gr_res1 = new TGraphErrors(n,&x_ptthree_[0],&y_mean_ratio_res1_[0],&ex_ptthree_[0],&ey_mean_ratio_res1_[0]);

  // fit linear extrapolation function
  TF1 *lin_extrapol = new TF1("lin_extrapol","[0]+[1]*x",0,Outer_->cutNumbers_.back()+0.05); //was used before...
  TF1 *lin_extrapol_mc = new TF1("lin_extrapol_mc","[0]+[1]*x",0,Outer_->cutNumbers_.back()+0.05); //was used before...
  TF1 *lin_extrapol_data = new TF1("lin_extrapol_data","[0]+[1]*x",0,Outer_->cutNumbers_.back()+0.05);
  // lin_extrapol->SetParameters(1,-0.1);
  //lin_extrapol->SetParName(0,"ResZero");
  //lin_extrapol->SetParName(1,"slope");
  
  // fit sqrt of quadratic function 
  TF1 *quad_extrapol = new TF1("quad_extrapol","TMath::Sqrt(pow([0],2)+pow([1],2)*pow(x,2))",0,Outer_->cutNumbers_.back()+0.05); //was used before...
  //quad_extrapol->SetParLimits(0, 0.01, 0.3);
 //  quad_extrapol->FixParameter(0, 0.005);
//   quad_extrapol->FixParameter(1, 0.4);
  quad_extrapol->SetParameters(0, 0.01);
  quad_extrapol->SetParameters(1, 0.4);
  
  //fit extrapolation function to the TGraphErrors for data and MC  
  if(doLinExtrapol_) {
     // extrapol_MC->Fit("lin_extrapol","Q","same",0,Outer_->cutNumbers_.back()+0.05);
     // extrapol_Data->Fit("lin_extrapol","Q","same",0,Outer_->cutNumbers_.back()+0.05);
    
     // fit mc
     data.reset();
     data.x_val = x;
     data.y_val = MCy;
     data.y_cov.ResizeTo(Outer_->cutNumbers_.size(),Outer_->cutNumbers_.size());
     data.y_cov = y_cov_mc;

     double slope = 0.04;
     double d_slope = 0.04;
     double offset = MCHistos_.at(0)->GetBinContent(xBin_i);
     double d_offset = MCHistos_.at(0)->GetBinError(xBin_i);
     make_lin_fit(slope, d_slope, offset, d_offset);

     std::cout << "MC fit values: " << "slope: " << slope << " offset: " << offset << std::endl; 
     lin_extrapol_mc->SetParameter(0, offset);
     lin_extrapol_mc->SetParError(0, d_offset);
     lin_extrapol_mc->SetParameter(1, slope);
     lin_extrapol_mc->SetParError(1, d_slope);
     extrapol_MC->GetListOfFunctions()->Add(lin_extrapol_mc);

     data.reset();

     // fit data
     data.x_val = x;
     data.y_val = Datay;
     data.y_cov.ResizeTo(Outer_->cutNumbers_.size(),Outer_->cutNumbers_.size());
     data.y_cov = y_cov_data;

     slope = 0.04;
     d_slope = 0.04;
     offset = DataHistos_.at(0)->GetBinContent(xBin_i);
     d_offset = DataHistos_.at(0)->GetBinError(xBin_i);
     make_lin_fit(slope, d_slope, offset, d_offset);
     std::cout << "Data fit values: " << "slope: " << slope << " offset: " << offset << std::endl; 

     lin_extrapol_data->SetParameter(0, offset);
     lin_extrapol_data->SetParError(0, d_offset);
     lin_extrapol_data->SetParameter(1, slope);
     lin_extrapol_data->SetParError(1, d_slope);
     extrapol_Data->GetListOfFunctions()->Add(lin_extrapol_data);

     data.reset();
  }
  else {
     extrapol_MC->Fit("quad_extrapol","Q","same",0,Outer_->cutNumbers_.back()+0.05);
     // std::cout << "MC Par 1: " << quad_extrapol->GetParameter(1) << std::endl;
     extrapol_Data->Fit("quad_extrapol","Q","same",0,Outer_->cutNumbers_.back()+0.05);
     // std::cout << "Data Par 1: " << quad_extrapol->GetParameter(1) << std::endl;
     //extrapol_MC->GetListOfFunctions()->Add(quad_extrapol);
     //extrapol_Data->GetListOfFunctions()->Add(quad_extrapol);
  }
  
  extrapol_MCDataRatio->Fit("lin_extrapol","Q","same",0,Outer_->cutNumbers_.back()+0.05);
  extrapol_NormalizedMCDataRatio->Fit("lin_extrapol","Q","same",0,Outer_->cutNumbers_.back()+0.05);

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
  //  std::pair <float,float> minMaxPair = Outer_->determineMinMax(DataExtrapols_.at(xBin_i));
  std::pair <float,float> minMaxPair = Outer_->determineMinMax(MCExtrapols_.at(xBin_i));
  //  std::pair <float,float> minMaxPair = std::make_pair(-0.1,0.2);
  //c->DrawFrame(0,minMaxPair.first*0.5-0.05,Outer_->cutNumbers_.back()+0.05,minMaxPair.second*1.2,(";cut on "+Outer_->configs_.at(0)->cutAxisTitle()+";"+Outer_->yProfileTitle()/*"#sqrt{2} #sigma"*/)/*.c_str()*/);
  c->DrawFrame(0,minMaxPair.first*0.5-0.05,Outer_->cutNumbers_.back()+0.05,minMaxPair.second*1.47,(";cut on "+Outer_->configs_.at(0)->cutAxisTitle()+";"+Outer_->yProfileTitle()/*"#sqrt{2} #sigma"*/)/*.c_str()*/);
  MCExtrapols_.at(xBin_i)->Draw("P");
  MCExtrapols_.at(xBin_i)->SetLineColor(style.getColor(0));
  MCExtrapols_.at(xBin_i)->SetMarkerColor(style.getColor(0));
  MCExtrapols_.at(xBin_i)->SetMarkerStyle(style.getMarker(0));
  DataExtrapols_.at(xBin_i)->Draw("Psame");
  DataExtrapols_.at(xBin_i)->SetLineColor(style.getColor(1));
  DataExtrapols_.at(xBin_i)->SetMarkerColor(style.getColor(1));
  DataExtrapols_.at(xBin_i)->SetMarkerStyle(style.getMarker(1));
  TF1* MCTemp = new TF1();
  TF1* DataTemp = new TF1();
  if(doLinExtrapol_) {
     MCExtrapols_.at(xBin_i)->GetFunction("lin_extrapol_mc")->SetLineColor(style.getColor(0));
     MCExtrapols_.at(xBin_i)->GetFunction("lin_extrapol_mc")->SetLineStyle(2);
     DataExtrapols_.at(xBin_i)->GetFunction("lin_extrapol_data")->SetLineColor(style.getColor(1));
     DataExtrapols_.at(xBin_i)->GetFunction("lin_extrapol_data")->SetLineStyle(2);
     MCTemp=(TF1*) MCExtrapols_.at(xBin_i)->GetFunction("lin_extrapol_mc")->Clone();
     DataTemp=(TF1*) DataExtrapols_.at(xBin_i)->GetFunction("lin_extrapol_data")->Clone();
  }
  else {
     MCExtrapols_.at(xBin_i)->GetFunction("quad_extrapol")->SetLineColor(style.getColor(0));
     MCExtrapols_.at(xBin_i)->GetFunction("quad_extrapol")->SetLineStyle(2);
     DataExtrapols_.at(xBin_i)->GetFunction("quad_extrapol")->SetLineColor(style.getColor(1));
     DataExtrapols_.at(xBin_i)->GetFunction("quad_extrapol")->SetLineStyle(2);
     MCTemp=(TF1*) MCExtrapols_.at(xBin_i)->GetFunction("quad_extrapol")->Clone();
     DataTemp=(TF1*) DataExtrapols_.at(xBin_i)->GetFunction("quad_extrapol")->Clone();
  }

  TPaveText *label_chi2 = util::LabelFactory::createPaveTextWithOffset(2,1.0,0.001);
  //TPaveText *label_chi2 = util::LabelFactory::createPaveTextWithOffset(1,1.0,0.001);
  label_chi2->AddText(Form("#chi2 / ndf Data = %4.2f/%i", DataTemp->GetChisquare(), DataTemp->GetNDF()));
  label_chi2->AddText(Form("#chi2 / ndf MC = %4.2f/%i", MCTemp->GetChisquare(), MCTemp->GetNDF()));
  label_chi2->Draw("same");

  MCTemp->SetRange(0.1,1);
  MCTemp->SetLineStyle(1);
  MCTemp->Draw("same");
  DataTemp->SetRange(0.1,1);
  DataTemp->SetLineStyle(1);
  DataTemp->Draw("same");
 
  //  TH1F* DrawHist = (TH1F*) MCExtrapols_.at(xBin_i)->GetHistogram();
//  DrawHist->GetXaxis()->SetRangeUser(0.0,Outer_->cutNumbers_.back()+0.05);
//  DrawHist->GetYaxis()->SetRangeUser(0.0,0.15);

  int nEntries =2;

  TPaveText *label = util::LabelFactory::createPaveTextWithOffset(2,1.0,0.5);
  label->AddText(Outer_->jetLabel_);//+  |#eta_{1,2}| > "+util::toTString(1.4)+",  L = "+util::StyleSettings::luminosity(4.6));
  label->AddText((Outer_->configs_.at(0)->xBinTitle(xBin_i,Outer_->configs_.at(0)->binEdges()->at(bin_i),Outer_->configs_.at(0)->binEdges()->at(bin_i+1),0)).c_str()/*+(TString)" and "*/);
  label->Draw("same");

  TLegend* leg1 = util::LabelFactory::createLegendWithOffset(2,0.6);
  leg1->AddEntry(MCExtrapols_.at(xBin_i),"Extrapolation (MC)","LP");
  leg1->AddEntry(DataExtrapols_.at(xBin_i),"Extrapolation (data)","LP");

  leg1->Draw();
  Outer_->drawCMSPrel();


  TString outname = "ResolutionPlots_ExtrapolPtThree_"+Outer_->plotsnames_+"_"+Outer_->configs_.at(0)->binName(bin_i)+"_"+Outer_->configs_.at(0)->xBinName(xBin_i);
  //  outname+=bin_i;
  MCExtrapols_.at(xBin_i)->SetName(outname);
  outname+=".eps";
  c->RedrawAxis();
  c->SaveAs(outname);
  Outer_->configs_.at(0)->safelyToRootFile(MCExtrapols_.at(xBin_i));


  //Plot Extrapol of MCDataRatio
  TCanvas* c2 = new TCanvas("c2","",600,600);
  minMaxPair = Outer_->determineMinMax(MCDataRatioExtrapols_.at(xBin_i));
  //  std::pair <float,float> minMaxPair = std::make_pair(-0.1,0.2);
  c2->DrawFrame(0,minMaxPair.first*0.85-0.05,Outer_->cutNumbers_.back()+0.05,minMaxPair.second*1.1,(";cut on "+Outer_->configs_.at(0)->cutAxisTitle()+";MC/Data Ratio"/*"#sqrt{2} #sigma"*/).c_str());
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
  outname+=".eps";
  c2->SaveAs(outname);
  Outer_->configs_.at(0)->safelyToRootFile(MCDataRatioExtrapols_.at(xBin_i));




  //Plot Extrapol of Normalized MCDataRatio
  TCanvas* c3 = new TCanvas("c3","",600,600);
  minMaxPair = Outer_->determineMinMax(NormalizedMCDataRatioExtrapols_.at(xBin_i));
  //  std::pair <float,float> minMaxPair = std::make_pair(-0.1,0.2);
  c3->DrawFrame(0,minMaxPair.first*0.85-0.05,Outer_->cutNumbers_.back()+0.05,minMaxPair.second*1.1,(";cut on "+Outer_->configs_.at(0)->cutAxisTitle()+";Normalized MC/Data Ratio"/*"#sqrt{2} #sigma"*/).c_str());
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
  outname+=".eps";
  c3->SaveAs(outname);




  chdir(".."); 


}

//! Get mean of alpha MC to use for extrapolation graphs with exclusive alpha bins
//! 
//! 
//!  \author Kristin Heine
//!  \date 2012/07/17
// ----------------------------------------------------------------   
std::vector<Double_t> Extrapolation::ExtrapolateBin::getMeanAlphaMC(Int_t xBin_i, Int_t bin_i){

   //std::cout << "bin_i: " << bin_i << " xBin_i: " << xBin_i << std::endl;

   // TFile* file = new TFile("/afs/naf.desy.de/user/k/kriheine/scratch/Kalibri/20132013ABCD_ReReco_CORR2013SummerV1_AK5_MC_Su12Z253_kostas_ChangeAlphaRangeMoreFitPointsAlphaExclusive_v2/dijetsFall10_TuneZ2_AK5PF_weighted_residuals_JERMatt/plots/KalibriPlots.root", "READ");

   //TFile* file = new TFile("/afs/naf.desy.de/user/k/kriheine/scratch/Kalibri/20132013ABCD_ReReco_CORR2013SummerV1_AK5_MC_Su12Z253_kostas_ChangeAlphaRangeMoreFitPointsAlphaExclusive_v5_fineAsymm/dijetsFall10_TuneZ2_AK5PF_weighted_residuals_JERMatt/plots/KalibriPlots.root", "READ");

   TFile* file = new TFile("/afs/naf.desy.de/user/k/kriheine/scratch/Kalibri/20132013ABCD_ReReco_CORR2013SummerV1_AK5_MC_Su12Z253_kostas_ChangeAlphaRangeAlphaExclusive_ThirdJetFraction_v1/dijetsFall10_TuneZ2_AK5PF_weighted_residuals_JERMatt/plots/KalibriPlots.root", "READ");

   TH1D* temp_0_MC;
   TH1D* temp_1_MC;
   TH1D* temp_2_MC;
   TH1D* temp_3_MC;
   TH1D* temp_4_MC;
   TH1D* temp_5_MC;
   TH1D* temp_6_MC;
   TH1D* temp_7_MC;
    
   TString name0_MC = Form("ThirdJetFractionPlainVsPt5_ThirdJetFractionVsMeanPt_MC_L2L3_AbsEta%i_Mean", bin_i);
   TString name1_MC = Form("ThirdJetFractionPlainVsPt10_ThirdJetFractionVsMeanPt_MC_L2L3_AbsEta%i_Mean", bin_i);
   TString name2_MC = Form("ThirdJetFractionPlainVsPt12.5_ThirdJetFractionVsMeanPt_MC_L2L3_AbsEta%i_Mean", bin_i);
   TString name3_MC = Form("ThirdJetFractionPlainVsPt15_ThirdJetFractionVsMeanPt_MC_L2L3_AbsEta%i_Mean", bin_i);
   TString name4_MC = Form("ThirdJetFractionPlainVsPt17.5_ThirdJetFractionVsMeanPt_MC_L2L3_AbsEta%i_Mean", bin_i);
   TString name5_MC = Form("ThirdJetFractionPlainVsPt20_ThirdJetFractionVsMeanPt_MC_L2L3_AbsEta%i_Mean", bin_i);
   TString name6_MC = Form("ThirdJetFractionPlainVsPt22.5_ThirdJetFractionVsMeanPt_MC_L2L3_AbsEta%i_Mean", bin_i);
   TString name7_MC = Form("ThirdJetFractionPlainVsPt25_ThirdJetFractionVsMeanPt_MC_L2L3_AbsEta%i_Mean", bin_i);

   //  file->cd("ThirdJetFractionPlainVsPt5");
   // temp_0_MC = (TH1D*) gDirectory->FindObjectAny(name0_MC);
   file->cd("ThirdJetFractionPlainVsPt10");
   temp_1_MC = (TH1D*) gDirectory->FindObjectAny(name1_MC);
   //file->cd("ThirdJetFractionPlainVsPt12.5");
   //temp_2_MC = (TH1D*) gDirectory->FindObjectAny(name2_MC);
   file->cd("ThirdJetFractionPlainVsPt15");
   temp_3_MC = (TH1D*) gDirectory->FindObjectAny(name3_MC);
   //file->cd("ThirdJetFractionPlainVsPt17.5");
   //temp_4_MC = (TH1D*) gDirectory->FindObjectAny(name4_MC);
   file->cd("ThirdJetFractionPlainVsPt20");
   temp_5_MC = (TH1D*) gDirectory->FindObjectAny(name5_MC);  
   //file->cd("ThirdJetFractionPlainVsPt22.5");
   //temp_6_MC = (TH1D*) gDirectory->FindObjectAny(name6_MC);
   file->cd("ThirdJetFractionPlainVsPt25");
   temp_7_MC = (TH1D*) gDirectory->FindObjectAny(name7_MC);
    
   std::vector<Double_t> temp_vec_MC;   
   //temp_vec_MC.push_back(temp_0_MC->GetBinContent(xBin_i));
   temp_vec_MC.push_back(temp_1_MC->GetBinContent(xBin_i));
   //temp_vec_MC.push_back(temp_2_MC->GetBinContent(xBin_i));
   temp_vec_MC.push_back(temp_3_MC->GetBinContent(xBin_i));
   //temp_vec_MC.push_back(temp_4_MC->GetBinContent(xBin_i));
   temp_vec_MC.push_back(temp_5_MC->GetBinContent(xBin_i));
   //temp_vec_MC.push_back(temp_6_MC->GetBinContent(xBin_i));
   temp_vec_MC.push_back(temp_7_MC->GetBinContent(xBin_i)); 

   file->Close();

   return  temp_vec_MC;
}

//! Get error on mean of alpha MC to use for extrapolation graphs with exclusive alpha bins
//! 
//! 
//!  \author Kristin Heine
//!  \date 2012/07/17
// ----------------------------------------------------------------   
std::vector<Double_t> Extrapolation::ExtrapolateBin::getMeanAlphaMCErr(Int_t xBin_i, Int_t bin_i){

   // TFile* file = new TFile("/afs/naf.desy.de/user/k/kriheine/scratch/Kalibri/20132013ABCD_ReReco_CORR2013SummerV1_AK5_MC_Su12Z253_kostas_ChangeAlphaRangeMoreFitPointsAlphaExclusive_v2/dijetsFall10_TuneZ2_AK5PF_weighted_residuals_JERMatt/plots/KalibriPlots.root", "READ");

   //TFile* file = new TFile("/afs/naf.desy.de/user/k/kriheine/scratch/Kalibri/20132013ABCD_ReReco_CORR2013SummerV1_AK5_MC_Su12Z253_kostas_ChangeAlphaRangeMoreFitPointsAlphaExclusive_v5_fineAsymm/dijetsFall10_TuneZ2_AK5PF_weighted_residuals_JERMatt/plots/KalibriPlots.root", "READ");

   TFile* file = new TFile("/afs/naf.desy.de/user/k/kriheine/scratch/Kalibri/20132013ABCD_ReReco_CORR2013SummerV1_AK5_MC_Su12Z253_kostas_ChangeAlphaRangeAlphaExclusive_ThirdJetFraction_v1/dijetsFall10_TuneZ2_AK5PF_weighted_residuals_JERMatt/plots/KalibriPlots.root", "READ");
  
   TH1D* temp_0_MC;
   TH1D* temp_1_MC;
   TH1D* temp_2_MC;
   TH1D* temp_3_MC;
   TH1D* temp_4_MC;
   TH1D* temp_5_MC;
   TH1D* temp_6_MC;
   TH1D* temp_7_MC;
    
   TString name0_MC = Form("ThirdJetFractionPlainVsPt5_ThirdJetFractionVsMeanPt_MC_L2L3_AbsEta%i_Mean", bin_i);
   TString name1_MC = Form("ThirdJetFractionPlainVsPt10_ThirdJetFractionVsMeanPt_MC_L2L3_AbsEta%i_Mean", bin_i);
   TString name2_MC = Form("ThirdJetFractionPlainVsPt12.5_ThirdJetFractionVsMeanPt_MC_L2L3_AbsEta%i_Mean", bin_i);
   TString name3_MC = Form("ThirdJetFractionPlainVsPt15_ThirdJetFractionVsMeanPt_MC_L2L3_AbsEta%i_Mean", bin_i);
   TString name4_MC = Form("ThirdJetFractionPlainVsPt17.5_ThirdJetFractionVsMeanPt_MC_L2L3_AbsEta%i_Mean", bin_i);
   TString name5_MC = Form("ThirdJetFractionPlainVsPt20_ThirdJetFractionVsMeanPt_MC_L2L3_AbsEta%i_Mean", bin_i);
   TString name6_MC = Form("ThirdJetFractionPlainVsPt22.5_ThirdJetFractionVsMeanPt_MC_L2L3_AbsEta%i_Mean", bin_i);
   TString name7_MC = Form("ThirdJetFractionPlainVsPt25_ThirdJetFractionVsMeanPt_MC_L2L3_AbsEta%i_Mean", bin_i);

   //  file->cd("ThirdJetFractionPlainVsPt5");
   //  temp_0_MC = (TH1D*) gDirectory->FindObjectAny(name0_MC);
   file->cd("ThirdJetFractionPlainVsPt10");
   temp_1_MC = (TH1D*) gDirectory->FindObjectAny(name1_MC);
   // file->cd("ThirdJetFractionPlainVsPt12.5");
   //temp_2_MC = (TH1D*) gDirectory->FindObjectAny(name2_MC);
   file->cd("ThirdJetFractionPlainVsPt15");
   temp_3_MC = (TH1D*) gDirectory->FindObjectAny(name3_MC);
   //file->cd("ThirdJetFractionPlainVsPt17.5");
   //temp_4_MC = (TH1D*) gDirectory->FindObjectAny(name4_MC);
   file->cd("ThirdJetFractionPlainVsPt20");
   temp_5_MC = (TH1D*) gDirectory->FindObjectAny(name5_MC);
   //file->cd("ThirdJetFractionPlainVsPt22.5");
   //temp_6_MC = (TH1D*) gDirectory->FindObjectAny(name6_MC);
   file->cd("ThirdJetFractionPlainVsPt25");
   temp_7_MC = (TH1D*) gDirectory->FindObjectAny(name7_MC);
     
   std::vector<Double_t> temp_vec_MC_err;
   //temp_vec_MC_err.push_back(temp_0_MC->GetBinError(xBin_i));
   temp_vec_MC_err.push_back(temp_1_MC->GetBinError(xBin_i));
   //temp_vec_MC_err.push_back(temp_2_MC->GetBinError(xBin_i));
   temp_vec_MC_err.push_back(temp_3_MC->GetBinError(xBin_i));
   //temp_vec_MC_err.push_back(temp_4_MC->GetBinError(xBin_i));
   temp_vec_MC_err.push_back(temp_5_MC->GetBinError(xBin_i));
   //temp_vec_MC_err.push_back(temp_6_MC->GetBinError(xBin_i));
   temp_vec_MC_err.push_back(temp_7_MC->GetBinError(xBin_i)); 

   file->Close();

   return temp_vec_MC_err;
}

//! Get mean of alpha data to use for extrapolation graphs with exclusive alpha bins
//! 
//! 
//!  \author Kristin Heine
//!  \date 2012/07/17
// ----------------------------------------------------------------   
std::vector<Double_t> Extrapolation::ExtrapolateBin::getMeanAlphaData(Int_t xBin_i, Int_t bin_i){

   //TFile* file = new TFile("/afs/naf.desy.de/user/k/kriheine/scratch/Kalibri/20132013ABCD_ReReco_CORR2013SummerV1_AK5_MC_Su12Z253_kostas_ChangeAlphaRangeMoreFitPointsAlphaExclusive_v2/dijetsFall10_TuneZ2_AK5PF_weighted_residuals_JERMatt/plots/KalibriPlots.root", "READ");

   // TFile* file = new TFile("/afs/naf.desy.de/user/k/kriheine/scratch/Kalibri/20132013ABCD_ReReco_CORR2013SummerV1_AK5_MC_Su12Z253_kostas_ChangeAlphaRangeMoreFitPointsAlphaExclusive_v5_fineAsymm/dijetsFall10_TuneZ2_AK5PF_weighted_residuals_JERMatt/plots/KalibriPlots.root", "READ");

   TFile* file = new TFile("/afs/naf.desy.de/user/k/kriheine/scratch/Kalibri/20132013ABCD_ReReco_CORR2013SummerV1_AK5_MC_Su12Z253_kostas_ChangeAlphaRangeAlphaExclusive_ThirdJetFraction_v1/dijetsFall10_TuneZ2_AK5PF_weighted_residuals_JERMatt/plots/KalibriPlots.root", "READ");
       
   TH1D* temp_0_Data;
   TH1D* temp_1_Data;
   TH1D* temp_2_Data;
   TH1D* temp_3_Data;
   TH1D* temp_4_Data;
   TH1D* temp_5_Data;
   TH1D* temp_6_Data;
   TH1D* temp_7_Data;

   TString name0_Data = Form("ThirdJetFractionPlainVsPt5_ThirdJetFractionVsMeanPt_data_L2L3res_AbsEta%i_Mean", bin_i);
   TString name1_Data = Form("ThirdJetFractionPlainVsPt10_ThirdJetFractionVsMeanPt_data_L2L3res_AbsEta%i_Mean", bin_i);
   TString name2_Data = Form("ThirdJetFractionPlainVsPt12.5_ThirdJetFractionVsMeanPt_data_L2L3res_AbsEta%i_Mean", bin_i);
   TString name3_Data = Form("ThirdJetFractionPlainVsPt15_ThirdJetFractionVsMeanPt_data_L2L3res_AbsEta%i_Mean", bin_i);
   TString name4_Data = Form("ThirdJetFractionPlainVsPt17.5_ThirdJetFractionVsMeanPt_data_L2L3res_AbsEta%i_Mean", bin_i);
   TString name5_Data = Form("ThirdJetFractionPlainVsPt20_ThirdJetFractionVsMeanPt_data_L2L3res_AbsEta%i_Mean", bin_i);
   TString name6_Data = Form("ThirdJetFractionPlainVsPt22.5_ThirdJetFractionVsMeanPt_data_L2L3res_AbsEta%i_Mean", bin_i);
   TString name7_Data = Form("ThirdJetFractionPlainVsPt25_ThirdJetFractionVsMeanPt_data_L2L3res_AbsEta%i_Mean", bin_i);
      
   // file->cd("ThirdJetFractionPlainVsPt5");
   // temp_0_Data = (TH1D*) gDirectory->FindObjectAny(name0_Data);
   file->cd("ThirdJetFractionPlainVsPt10");
   temp_1_Data = (TH1D*) gDirectory->FindObjectAny(name1_Data);
   //file->cd("ThirdJetFractionPlainVsPt12.5");
   //temp_2_Data = (TH1D*) gDirectory->FindObjectAny(name2_Data);
   file->cd("ThirdJetFractionPlainVsPt15");
   temp_3_Data = (TH1D*) gDirectory->FindObjectAny(name3_Data);
   //file->cd("ThirdJetFractionPlainVsPt17.5");
   //temp_4_Data = (TH1D*) gDirectory->FindObjectAny(name4_Data);
   file->cd("ThirdJetFractionPlainVsPt20");
   temp_5_Data = (TH1D*) gDirectory->FindObjectAny(name5_Data);
   //file->cd("ThirdJetFractionPlainVsPt22.5");
   //temp_6_Data = (TH1D*) gDirectory->FindObjectAny(name6_Data);
   file->cd("ThirdJetFractionPlainVsPt25");
   temp_7_Data = (TH1D*) gDirectory->FindObjectAny(name7_Data);
    
   std::vector<Double_t> temp_vec_Data;
   //temp_vec_Data.push_back(temp_0_Data->GetBinContent(xBin_i));
   temp_vec_Data.push_back(temp_1_Data->GetBinContent(xBin_i));
   // temp_vec_Data.push_back(temp_2_Data->GetBinContent(xBin_i));
   temp_vec_Data.push_back(temp_3_Data->GetBinContent(xBin_i));
   // temp_vec_Data.push_back(temp_4_Data->GetBinContent(xBin_i));
   temp_vec_Data.push_back(temp_5_Data->GetBinContent(xBin_i));
   // temp_vec_Data.push_back(temp_6_Data->GetBinContent(xBin_i));
   temp_vec_Data.push_back(temp_7_Data->GetBinContent(xBin_i));   

   file->Close();

   return temp_vec_Data;
}

//! Get error on mean of alpha data to use for extrapolation graphs with exclusive alpha bins
//! 
//! 
//!  \author Kristin Heine
//!  \date 2012/07/17
// ----------------------------------------------------------------   
std::vector<Double_t> Extrapolation::ExtrapolateBin::getMeanAlphaDataErr(Int_t xBin_i, Int_t bin_i){

   //TFile* file = new TFile("/afs/naf.desy.de/user/k/kriheine/scratch/Kalibri/20132013ABCD_ReReco_CORR2013SummerV1_AK5_MC_Su12Z253_kostas_ChangeAlphaRangeMoreFitPointsAlphaExclusive_v2/dijetsFall10_TuneZ2_AK5PF_weighted_residuals_JERMatt/plots/KalibriPlots.root", "READ");

   //TFile* file = new TFile("/afs/naf.desy.de/user/k/kriheine/scratch/Kalibri/20132013ABCD_ReReco_CORR2013SummerV1_AK5_MC_Su12Z253_kostas_ChangeAlphaRangeMoreFitPointsAlphaExclusive_v5_fineAsymm/dijetsFall10_TuneZ2_AK5PF_weighted_residuals_JERMatt/plots/KalibriPlots.root", "READ");

   TFile* file = new TFile("/afs/naf.desy.de/user/k/kriheine/scratch/Kalibri/20132013ABCD_ReReco_CORR2013SummerV1_AK5_MC_Su12Z253_kostas_ChangeAlphaRangeAlphaExclusive_ThirdJetFraction_v1/dijetsFall10_TuneZ2_AK5PF_weighted_residuals_JERMatt/plots/KalibriPlots.root", "READ");
    
   TH1D* temp_0_Data;
   TH1D* temp_1_Data;
   TH1D* temp_2_Data;
   TH1D* temp_3_Data;
   TH1D* temp_4_Data;
   TH1D* temp_5_Data;
   TH1D* temp_6_Data;
   TH1D* temp_7_Data;

   TString name0_Data = Form("ThirdJetFractionPlainVsPt5_ThirdJetFractionVsMeanPt_data_L2L3res_AbsEta%i_Mean", bin_i);
   TString name1_Data = Form("ThirdJetFractionPlainVsPt10_ThirdJetFractionVsMeanPt_data_L2L3res_AbsEta%i_Mean", bin_i);
   TString name2_Data = Form("ThirdJetFractionPlainVsPt12.5_ThirdJetFractionVsMeanPt_data_L2L3res_AbsEta%i_Mean", bin_i);
   TString name3_Data = Form("ThirdJetFractionPlainVsPt15_ThirdJetFractionVsMeanPt_data_L2L3res_AbsEta%i_Mean", bin_i);
   TString name4_Data = Form("ThirdJetFractionPlainVsPt17.5_ThirdJetFractionVsMeanPt_data_L2L3res_AbsEta%i_Mean", bin_i);
   TString name5_Data = Form("ThirdJetFractionPlainVsPt20_ThirdJetFractionVsMeanPt_data_L2L3res_AbsEta%i_Mean", bin_i);
   TString name6_Data = Form("ThirdJetFractionPlainVsPt22.5_ThirdJetFractionVsMeanPt_data_L2L3res_AbsEta%i_Mean", bin_i);
   TString name7_Data = Form("ThirdJetFractionPlainVsPt25_ThirdJetFractionVsMeanPt_data_L2L3res_AbsEta%i_Mean", bin_i);
      
   // file->cd("ThirdJetFractionPlainVsPt5");
   // temp_0_Data = (TH1D*) gDirectory->FindObjectAny(name0_Data);
   file->cd("ThirdJetFractionPlainVsPt10");
   temp_1_Data = (TH1D*) gDirectory->FindObjectAny(name1_Data);
   //file->cd("ThirdJetFractionPlainVsPt12.5");
   //temp_2_Data = (TH1D*) gDirectory->FindObjectAny(name2_Data);
   file->cd("ThirdJetFractionPlainVsPt15");
   temp_3_Data = (TH1D*) gDirectory->FindObjectAny(name3_Data);
   //file->cd("ThirdJetFractionPlainVsPt17.5");
   //temp_4_Data = (TH1D*) gDirectory->FindObjectAny(name4_Data);
   file->cd("ThirdJetFractionPlainVsPt20");
   temp_5_Data = (TH1D*) gDirectory->FindObjectAny(name5_Data);
   //file->cd("ThirdJetFractionPlainVsPt22.5");
   //temp_6_Data = (TH1D*) gDirectory->FindObjectAny(name6_Data);
   file->cd("ThirdJetFractionPlainVsPt25");
   temp_7_Data = (TH1D*) gDirectory->FindObjectAny(name7_Data);
       
   std::vector<double> temp_vec_Data_err;
   //temp_vec_Data_err.push_back(temp_0_Data->GetBinError(xBin_i));
   temp_vec_Data_err.push_back(temp_1_Data->GetBinError(xBin_i));
   //temp_vec_Data_err.push_back(temp_2_Data->GetBinError(xBin_i));
   temp_vec_Data_err.push_back(temp_3_Data->GetBinError(xBin_i));
   //temp_vec_Data_err.push_back(temp_4_Data->GetBinError(xBin_i));
   temp_vec_Data_err.push_back(temp_5_Data->GetBinError(xBin_i));
   //temp_vec_Data_err.push_back(temp_6_Data->GetBinError(xBin_i));
   temp_vec_Data_err.push_back(temp_7_Data->GetBinError(xBin_i)); 

   file->Close();

   return temp_vec_Data_err;
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

  for(Int_t xbin_i=0;xbin_i<Outer_->configs_.at(0)->nXBins();xbin_i++){
     if(doLinExtrapol_) {
        tempExtrapolatedResMC->SetBinContent(xbin_i+1,MCExtrapols_.at(xbin_i)->GetFunction("lin_extrapol_mc")->GetParameter(0)>0.005 ? MCExtrapols_.at(xbin_i)->GetFunction("lin_extrapol_mc")->GetParameter(0) : 0.0);
        tempExtrapolatedResData->SetBinContent(xbin_i+1,DataExtrapols_.at(xbin_i)->GetFunction("lin_extrapol_data")->GetParameter(0)>0.005 ? DataExtrapols_.at(xbin_i)->GetFunction("lin_extrapol_data")->GetParameter(0) : 0.0);
     }
     else {
        tempExtrapolatedResMC->SetBinContent(xbin_i+1,MCExtrapols_.at(xbin_i)->GetFunction("quad_extrapol")->GetParameter(0)>0.005 ? MCExtrapols_.at(xbin_i)->GetFunction("quad_extrapol")->GetParameter(0) : 0.0);
        tempExtrapolatedResData->SetBinContent(xbin_i+1,DataExtrapols_.at(xbin_i)->GetFunction("quad_extrapol")->GetParameter(0)>0.005 ? DataExtrapols_.at(xbin_i)->GetFunction("quad_extrapol")->GetParameter(0) : 0.0);
     }
    tempExtrapolatedMCDataRatio->SetBinContent(xbin_i+1,MCDataRatioExtrapols_.at(xbin_i)->GetFunction("lin_extrapol")->GetParameter(0)>0.005 ? MCDataRatioExtrapols_.at(xbin_i)->GetFunction("lin_extrapol")->GetParameter(0) : 0.0);
    tempExtrapolatedNormalizedMCDataRatio->SetBinContent(xbin_i+1,NormalizedMCDataRatioExtrapols_.at(xbin_i)->GetFunction("lin_extrapol")->GetParameter(0)>0.005 ? NormalizedMCDataRatioExtrapols_.at(xbin_i)->GetFunction("lin_extrapol")->GetParameter(0) : 0.0);
    if(doLinExtrapol_) {
       tempExtrapolatedResMC->SetBinError(xbin_i+1,MCExtrapols_.at(xbin_i)->GetFunction("lin_extrapol_mc")->GetParameter(0)>0.005 ? MCExtrapols_.at(xbin_i)->GetFunction("lin_extrapol_mc")->GetParError(0) : 0.0);
       tempExtrapolatedResData->SetBinError(xbin_i+1,DataExtrapols_.at(xbin_i)->GetFunction("lin_extrapol_data")->GetParameter(0)>0.005 ? DataExtrapols_.at(xbin_i)->GetFunction("lin_extrapol_data")->GetParError(0) : 0.0);
    }
    else {
       tempExtrapolatedResMC->SetBinError(xbin_i+1,MCExtrapols_.at(xbin_i)->GetFunction("quad_extrapol")->GetParameter(0)>0.005 ? MCExtrapols_.at(xbin_i)->GetFunction("quad_extrapol")->GetParError(0) : 0.0);
       tempExtrapolatedResData->SetBinError(xbin_i+1,DataExtrapols_.at(xbin_i)->GetFunction("quad_extrapol")->GetParameter(0)>0.005 ? DataExtrapols_.at(xbin_i)->GetFunction("quad_extrapol")->GetParError(0) : 0.0);
    }
    tempExtrapolatedMCDataRatio->SetBinError(xbin_i+1,MCDataRatioExtrapols_.at(xbin_i)->GetFunction("lin_extrapol")->GetParameter(0)>0.005 ? MCDataRatioExtrapols_.at(xbin_i)->GetFunction("lin_extrapol")->GetParError(0) : 0.0);
    tempExtrapolatedNormalizedMCDataRatio->SetBinError(xbin_i+1,NormalizedMCDataRatioExtrapols_.at(xbin_i)->GetFunction("lin_extrapol")->GetParameter(0)>0.005 ? NormalizedMCDataRatioExtrapols_.at(xbin_i)->GetFunction("lin_extrapol")->GetParError(0) : 0.0);
  }
  ExtrapolatedResMC_=(TH1D*) tempExtrapolatedResMC;
  ExtrapolatedResData_=(TH1D*) tempExtrapolatedResData;
  ExtrapolatedMCDataRatio_=(TH1D*) tempExtrapolatedMCDataRatio;
  ExtrapolatedNormalizedMCDataRatio_=(TH1D*) tempExtrapolatedNormalizedMCDataRatio;
    
}
