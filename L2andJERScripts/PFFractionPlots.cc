#include "PFFractionPlots.h"
#include "Extrapolation.h"

//Default constructor, see BasePlotExtractor constructor
//! for more information and initialization.
//!
//!  \author Henning Kirschenmann
//!  \date 2012/03/07
// ----------------------------------------------------------------   
PFFractionPlots::PFFractionPlots(TString plotsnames,TString kalibriPlotsShortName) : BasePlotExtractor(plotsnames,kalibriPlotsShortName), drawResponseCorrected_("0"){
  //init config_file;
  //  config_=ConfigFile("L2L3.cfg");
  //  plotsnames_=plotsnames;
  init("Mean");
}


//! Plot and save all the histograms
//!
//!  \author Henning Kirschenmann
//!  \date 2012/03/07
// ----------------------------------------------------------------   
void PFFractionPlots::Plot() {
  setTDRStyle();
  TCanvas* c = new TCanvas("c","",600,600);
  c->SetLogx(configs_.at(0)->logX());
  DefaultStyles style;
  style.setStyle("PFComp");
  //  MakeDateDir();
  if(chdir("PFComp") != 0){ 
    mkdir("PFComp", S_IRWXU|S_IRWXG|S_IRWXO); 
    chdir("PFComp"); 
  } 

  TString PFFractionTitle= "PF energy fraction";
  if(drawResponseCorrected_.Contains("MPF")){
    PFFractionTitle = "PF energy fraction (MPF corr.)";
    yRatioTitle_.ReplaceAll("fractions","fractions (MPF corr.)");
    yDifferenceTitle_.ReplaceAll("fractions","fractions (MPF corr.)");
  }
  if(drawResponseCorrected_.Contains("RR")){
    PFFractionTitle = "PF energy fraction (RR corr.)";
    yRatioTitle_.ReplaceAll("fractions","fractions (RR corr.)");
    yDifferenceTitle_.ReplaceAll("fractions","fractions (RR corr.)");
    yRatioTitle_.ReplaceAll("(MPF corr.)","");
    yDifferenceTitle_.ReplaceAll("(MPF corr.)","");
  }

//
//  for(int i=0;i<names_.size();i++){
//    std::cout << names_.at(i) << std::endl;
//  }
//
  //  for(int bin_i=0;bin_i<pConfig->nBins();bin_i++){
  for(int bin_i=0;bin_i<configs_.at(0)->nBins();bin_i++){
  int nEntries =names_.size();
  TLegend *leg = new TLegend(0.2,0.45-(nEntries/2)*0.07,0.8,0.45);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetTextFont(42);
  leg->SetTextSize(0.04);
  leg->SetNColumns(2);

    for(int conf_i=0;conf_i<configs_.size();conf_i++){
      for(int j=0;j<2;j++){
	if(j==1){
	  AllPlots_.at(conf_i).at(bin_i).at(j)->SetLineColor(style.getColor(conf_i));
	  AllPlots_.at(conf_i).at(bin_i).at(j)->SetMarkerColor(1);
	}
	else {
	  AllPlots_.at(conf_i).at(bin_i).at(j)->SetLineColor(1);
	  AllPlots_.at(conf_i).at(bin_i).at(j)->SetMarkerColor(1);
	  
	}
	  AllPlots_.at(conf_i).at(bin_i).at(j)->SetFillColor(style.getColor(conf_i));
	  AllPlots_.at(conf_i).at(bin_i).at(j)->SetFillStyle(1001);
	  AllPlots_.at(conf_i).at(bin_i).at(j)->SetMarkerStyle(style.getMarker(conf_i));
	  AllPlots_.at(conf_i).at(bin_i).at(j)->GetListOfFunctions()->Delete();

      }
      leg->AddEntry(AllPlots_.at(conf_i).at(bin_i).at(1),(configs_.at(conf_i)->yTitle()).c_str(),"FP");
    }


  

    //    std::cout << (configs_.at(0)->binTitle(configs_.at(0)->binEdges()->at(bin_i),
    //					   configs_.at(0)->binEdges()->at(bin_i+1))).c_str() << std::endl;

    AllPlots_.at(0).at(bin_i).at(1)->Draw();
    //    AllPlots_.at(0).at(bin_i).at(1)->SetTitle((configs_.at(0)->xBinTitle(bin_i,20,2000)).c_str());
    AllPlots_.at(0).at(bin_i).at(1)->GetYaxis()->SetRangeUser(0,1.1);
    AllPlots_.at(0).at(bin_i).at(1)->GetYaxis()->SetTitle(PFFractionTitle);
    //    AllPlots_.at(0).at(bin_i).at(1)->GetYaxis()->SetTitleOffset(1.1);
    TStyle *tdrStyle = (TStyle*)gROOT->FindObject("tdrStyle"); 
    AllPlots_.at(0).at(bin_i).at(1)->GetYaxis()->SetTitleOffset(tdrStyle->GetTitleYOffset());
    AllPlots_.at(0).at(bin_i).at(1)->GetXaxis()->SetRangeUser(configs_.at(0)->xMin(),configs_.at(0)->xMax());
  THStack *hs = new THStack("hs","Stacked 1D histograms MC ");
  THStack *hs2 = new THStack("hs2","Stacked 1D histograms DATA");
  for(int i=0;i<configs_.size();i++){
  hs->Add(AllPlots_.at(i).at(bin_i).at(1));
  hs2->Add(AllPlots_.at(i).at(bin_i).at(0));
  }
  hs->Draw("histsame");
  hs2->Draw("same");
  drawRunNumberLabels(AllPlots_.at(0).at(bin_i).at(1),configs_.at(0));
  leg->SetHeader((configs_.at(0)->binTitle(configs_.at(0)->binEdges()->at(bin_i),
					   configs_.at(0)->binEdges()->at(bin_i+1))).c_str());
  leg->Draw();
  drawCMSPrel();
  TString outname = "FractionPlots_"+plotsnames_+"_"+configs_.at(0)->binName(bin_i);
  outname+="RespCorr";
  outname+=drawResponseCorrected_;
  outname+=".pdf";
  c->RedrawAxis();
  c->SaveAs(outname);

  //Ratio plots 
  //Ratio plots 
  //Ratio plots 

  AllRatiosDataMC_.at(0).at(bin_i)->GetListOfFunctions()->Delete();
  AllRatiosDataMC_.at(0).at(bin_i)->Draw();
  for(int conf_i=0;conf_i<configs_.size();conf_i++){
    TStyle *tdrStyle = (TStyle*)gROOT->FindObject("tdrStyle"); 
    TList *lfits = AllRatiosDataMC_.at(conf_i).at(bin_i)->GetListOfFunctions();
    lfits->Delete();
    //Draw DataMC-ratios
    if(DEBUG)std::cout << AllRatiosDataMC_.size() << " and " << std::endl;
    AllRatiosDataMC_.at(conf_i).at(bin_i)->Draw("same");
    drawConfidenceIntervals(AllRatiosDataMC_.at(conf_i).at(bin_i));
    AllRatiosDataMC_.at(conf_i).at(bin_i)->Draw("same");
    AllRatiosDataMC_.at(conf_i).at(bin_i)->GetYaxis()->SetRangeUser(yRatioMinMax().at(0),yRatioMinMax().at(1));
    AllRatiosDataMC_.at(conf_i).at(bin_i)->GetYaxis()->SetTitle(yRatioTitle());
    //    AllRatiosDataMC_.at(i).at(bin_i)->GetYaxis()->SetRangeUser(0.5,1.5);
    //    AllRatiosDataMC_.at(i).at(bin_i)->GetYaxis()->SetTitle("Resolution ratio");
    AllRatiosDataMC_.at(conf_i).at(bin_i)->GetYaxis()->SetTitleOffset(tdrStyle->GetTitleYOffset());
    AllRatiosDataMC_.at(conf_i).at(bin_i)->GetXaxis()->SetRangeUser(configs_.at(0)->xMin(),configs_.at(0)->xMax());


    AllRatiosDataMC_.at(conf_i).at(bin_i)->SetLineColor(style.getColor(conf_i));
    AllRatiosDataMC_.at(conf_i).at(bin_i)->SetMarkerColor(style.getColor(conf_i));
    AllRatiosDataMC_.at(conf_i).at(bin_i)->SetMarkerStyle(style.getMarker(conf_i));
    AllRatiosDataMC_.at(conf_i).at(bin_i)->SetFillColor(style.getColor(conf_i));
    AllRatiosDataMC_.at(conf_i).at(bin_i)->SetFillStyle(0);
  }

  //  label->Draw("same");
  drawRunNumberLabels(AllRatiosDataMC_.at(0).at(bin_i),configs_.at(0));
  leg->Draw();
     drawCMSPrel();
    outname = "FractionPlots_"+plotsnames_+"_ratio"+"_"+configs_.at(0)->binName(bin_i);
    outname+="RespCorr";
  outname+=drawResponseCorrected_;
    outname+=".pdf";
    c->RedrawAxis();
    c->SaveAs(outname);


  //Difference plots 
  //Difference plots 
  //Difference plots 

  AllDifferencesDataMC_.at(0).at(bin_i)->GetListOfFunctions()->Delete();
  AllDifferencesDataMC_.at(0).at(bin_i)->Draw();
  for(int conf_i=0;conf_i<configs_.size();conf_i++){
    TStyle *tdrStyle = (TStyle*)gROOT->FindObject("tdrStyle"); 
    TList *lfits = AllDifferencesDataMC_.at(conf_i).at(bin_i)->GetListOfFunctions();
    lfits->Delete();
    //Draw DataMC-differences
    if(DEBUG)std::cout << AllDifferencesDataMC_.size() << " and " << std::endl;
    AllDifferencesDataMC_.at(conf_i).at(bin_i)->Draw("same");
    drawConfidenceIntervals(AllDifferencesDataMC_.at(conf_i).at(bin_i));
    AllDifferencesDataMC_.at(conf_i).at(bin_i)->Draw("same");
    AllDifferencesDataMC_.at(conf_i).at(bin_i)->GetYaxis()->SetRangeUser(yDifferenceMinMax().at(0),yDifferenceMinMax().at(1));
    AllDifferencesDataMC_.at(conf_i).at(bin_i)->GetYaxis()->SetTitle(yDifferenceTitle());
    AllDifferencesDataMC_.at(conf_i).at(bin_i)->GetYaxis()->SetTitleOffset(tdrStyle->GetTitleYOffset());
    AllDifferencesDataMC_.at(conf_i).at(bin_i)->GetXaxis()->SetRangeUser(configs_.at(0)->xMin(),configs_.at(0)->xMax());


    AllDifferencesDataMC_.at(conf_i).at(bin_i)->SetLineColor(style.getColor(conf_i));
    AllDifferencesDataMC_.at(conf_i).at(bin_i)->SetMarkerColor(style.getColor(conf_i));
    AllDifferencesDataMC_.at(conf_i).at(bin_i)->SetMarkerStyle(style.getMarker(conf_i));
    AllDifferencesDataMC_.at(conf_i).at(bin_i)->SetFillColor(style.getColor(conf_i));
    AllDifferencesDataMC_.at(conf_i).at(bin_i)->SetFillStyle(0);
  }

  //  label->Draw("same");
  drawRunNumberLabels(AllDifferencesDataMC_.at(0).at(bin_i),configs_.at(0));
  leg->Draw();
     drawCMSPrel();
    outname = "FractionPlots_"+plotsnames_+"_difference"+"_"+configs_.at(0)->binName(bin_i);
    outname+="RespCorr";
    outname+=drawResponseCorrected_;
    outname+=".pdf";
    c->RedrawAxis();
    c->SaveAs(outname);









  }

  chdir(outputPathROOT()+"/..");

}


void PFFractionPlots::PlotDifferencesWithResponseCorrection() {



  Extrapolation* MPFResponseResults = new Extrapolation("MPFResponseVsRunNumber",kalibriPlotsShortName());
  chdir(outputPathROOT()+"/"+kalibriPlotsShortName_+"/"+plotsnames_);
  std::cout << "indexToNormalizeTo_ " << MPFResponseResults->indexToNormalizeTo_  << std::endl;
  for(int bin_i=0;bin_i<configs_.at(0)->nBins();bin_i++){
    for(int conf_i=0;conf_i<configs_.size();conf_i++){
      for(int j=0;j<2;j++){
	for(int xbin_i=0;xbin_i<AllPlots_.at(conf_i).at(bin_i).at(j)->GetNbinsX();xbin_i++){
	  Double_t FractionValue = AllPlots_.at(conf_i).at(bin_i).at(j)->GetBinContent(xbin_i);
	  Double_t FractionError = AllPlots_.at(conf_i).at(bin_i).at(j)->GetBinError(xbin_i);
	  Double_t ResponseValue = MPFResponseResults->AllPlots_.at(MPFResponseResults->indexToNormalizeTo_).at(bin_i).at(j)->GetBinContent(xbin_i);
	  Double_t ResponseError = MPFResponseResults->AllPlots_.at(MPFResponseResults->indexToNormalizeTo_).at(bin_i).at(j)->GetBinError(xbin_i);
	  if(DEBUGALL)std::cout << "conf_i " << conf_i << " bin_i " << bin_i << " j " << j << " xbin_i " << xbin_i << " AllPlots_.at(conf_i).at(bin_i).at(j)->GetBinContent(xbin_i) " << FractionValue << " MPFResponseResults->AllPlots_.at(MPFResponseResults->indexToNormalizeTo_).at(bin_i).at(j)->GetBinContent(xbin_i) " << ResponseValue << std::endl;
	  AllPlots_.at(conf_i).at(bin_i).at(j)->SetBinContent(xbin_i,FractionValue*ResponseValue);
	  // plain f(x) = u * v gaussian error propagation:
	  AllPlots_.at(conf_i).at(bin_i).at(j)->SetBinError(xbin_i,TMath::Sqrt(TMath::Power(FractionValue*ResponseError,2)+TMath::Power(ResponseValue*FractionError,2)));
	}
      }
    }
  }
  refreshRatiosDataMC();
  drawResponseCorrected_="MPF";
  Plot();



  Extrapolation* RRResponseResults = new Extrapolation("RelResponseVsRunNumber",kalibriPlotsShortName());
  chdir(outputPathROOT()+"/"+kalibriPlotsShortName_+"/"+plotsnames_);
  std::cout << "indexToNormalizeTo_ " << RRResponseResults->indexToNormalizeTo_  << std::endl;
  for(int bin_i=0;bin_i<configs_.at(0)->nBins();bin_i++){
    for(int conf_i=0;conf_i<configs_.size();conf_i++){
      for(int j=0;j<2;j++){
	for(int xbin_i=0;xbin_i<AllPlots_.at(conf_i).at(bin_i).at(j)->GetNbinsX();xbin_i++){
	  Double_t FractionValue = AllPlots_.at(conf_i).at(bin_i).at(j)->GetBinContent(xbin_i);
	  Double_t FractionError = AllPlots_.at(conf_i).at(bin_i).at(j)->GetBinError(xbin_i);
	  Double_t ResponseValue = RRResponseResults->AllPlots_.at(RRResponseResults->indexToNormalizeTo_).at(bin_i).at(j)->GetBinContent(xbin_i);
	  Double_t ResponseError = RRResponseResults->AllPlots_.at(RRResponseResults->indexToNormalizeTo_).at(bin_i).at(j)->GetBinError(xbin_i);
	  if(DEBUGALL)std::cout << "conf_i " << conf_i << " bin_i " << bin_i << " j " << j << " xbin_i " << xbin_i << " AllPlots_.at(conf_i).at(bin_i).at(j)->GetBinContent(xbin_i) " << FractionValue << " RRResponseResults->AllPlots_.at(RRResponseResults->indexToNormalizeTo_).at(bin_i).at(j)->GetBinContent(xbin_i) " << ResponseValue << std::endl;
	  AllPlots_.at(conf_i).at(bin_i).at(j)->SetBinContent(xbin_i,FractionValue*ResponseValue);
	  // plain f(x) = u * v gaussian error propagation:
	  AllPlots_.at(conf_i).at(bin_i).at(j)->SetBinError(xbin_i,TMath::Sqrt(TMath::Power(FractionValue*ResponseError,2)+TMath::Power(ResponseValue*FractionError,2)));
	}
      }
    }
  }
  refreshRatiosDataMC();
  drawResponseCorrected_="RR";
  Plot();

//  // somehow RR/MPFResponseResults seems to interfer with the PFFRactionPlots instance, deleting ResponseResults crashes destructor of PFFractionPlots instance. Any reason for this?
//  // keep memory leak as a workaround (as this is just a plotting macro)
//  // might be linked to BasePlotExtractor::refreshRatiosDataMC() ?
//  delete MPFResponseResults;
//  delete RRResponseResults;

}
