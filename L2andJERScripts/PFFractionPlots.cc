#include "PFFractionPlots.h"

//Default constructor, see BasePlotExtractor constructor
//! for more information and initialization.
//!
//!  \author Henning Kirschenmann
//!  \date 2012/03/07
// ----------------------------------------------------------------   
PFFractionPlots::PFFractionPlots(TString plotsnames,TString kalibriPlotsShortName) : BasePlotExtractor(plotsnames,kalibriPlotsShortName){
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
	if(j==0){
	  AllPlots_.at(conf_i).at(bin_i).at(j)->SetLineColor(style.getColor(conf_i));
	}
	else {
	  AllPlots_.at(conf_i).at(bin_i).at(j)->SetLineColor(1);
	  AllPlots_.at(conf_i).at(bin_i).at(j)->SetMarkerColor(1);
	  
	}
	  AllPlots_.at(conf_i).at(bin_i).at(j)->SetFillColor(style.getColor(conf_i));
	  AllPlots_.at(conf_i).at(bin_i).at(j)->SetFillStyle(3001);
	  AllPlots_.at(conf_i).at(bin_i).at(j)->SetMarkerStyle(style.getMarker(conf_i));
      }
      leg->AddEntry(AllPlots_.at(conf_i).at(bin_i).at(1),(configs_.at(conf_i)->yTitle()).c_str(),"FP");
    }


  

    //    std::cout << (configs_.at(0)->binTitle(configs_.at(0)->binEdges()->at(bin_i),
    //					   configs_.at(0)->binEdges()->at(bin_i+1))).c_str() << std::endl;

    AllPlots_.at(0).at(bin_i).at(0)->Draw();
    //    AllPlots_.at(0).at(bin_i).at(0)->SetTitle((configs_.at(0)->xBinTitle(bin_i,20,2000)).c_str());
    AllPlots_.at(0).at(bin_i).at(0)->GetYaxis()->SetRangeUser(0,1.1);
    AllPlots_.at(0).at(bin_i).at(0)->GetYaxis()->SetTitle("PF energy fraction");
    //    AllPlots_.at(0).at(bin_i).at(0)->GetYaxis()->SetTitleOffset(1.1);
    TStyle *tdrStyle = (TStyle*)gROOT->FindObject("tdrStyle"); 
    AllPlots_.at(0).at(bin_i).at(0)->GetYaxis()->SetTitleOffset(tdrStyle->GetTitleYOffset());
    AllPlots_.at(0).at(bin_i).at(0)->GetXaxis()->SetRangeUser(configs_.at(0)->xMin(),configs_.at(0)->xMax());
  THStack *hs = new THStack("hs","Stacked 1D histograms MC ");
  THStack *hs2 = new THStack("hs2","Stacked 1D histograms DATA");
  for(int i=0;i<configs_.size();i++){
  hs->Add(AllPlots_.at(i).at(bin_i).at(0));
  hs2->Add(AllPlots_.at(i).at(bin_i).at(1));
  }
  hs->Draw("histsame");
  hs2->Draw("same");
  leg->SetHeader((configs_.at(0)->binTitle(configs_.at(0)->binEdges()->at(bin_i),
					   configs_.at(0)->binEdges()->at(bin_i+1))).c_str());
  leg->Draw();
  cmsPrel();
  TString outname = "FractionPlots_"+plotsnames_+"_"+configs_.at(0)->binName(bin_i);
  //  outname+=bin_i;
  outname+=".eps";
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
    std::cout << AllRatiosDataMC_.size() << " and " << std::endl;
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
  leg->Draw();
     cmsPrel();
    outname = "FractionPlots_"+plotsnames_+"_ratio"+"_"+configs_.at(0)->binName(bin_i);
    //  outname+=bin_i;
    outname+=".eps";
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
    std::cout << AllDifferencesDataMC_.size() << " and " << std::endl;
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
  leg->Draw();
     cmsPrel();
    outname = "FractionPlots_"+plotsnames_+"_difference"+"_"+configs_.at(0)->binName(bin_i);
    //  outname+=bin_i;
    outname+=".eps";
    c->SaveAs(outname);









  }

  chdir(outputPathROOT()+"/..");

}

