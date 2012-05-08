/*************************************************
  Plot histograms together, then plot
  their ratio on a histogram below.

  Must provied at least 1 histogram 
  to be the "denomintator", and
  up to 3 histograms to be "numerators".

  All the plots will be shown together in a
  big plot, and below that
  will be a plot showing the ratio of all
  the "numerators" to the "denominator".

  
  Example use (first step is to prevent a bug in ROOT crashing):
   root [1] .O0
   root [2] .L plotHistswSomeUncAndAbsUncertainty.C
   root [3] plotHistswSomeUncAndAbsUncertainty(photonEt_BarrelPIDcut,photonEt_Barrel,"Photon E_{T}","E_{T} (GeV)","entries/20 GeV bin")  

  Michael Anderson
  March 18, 2009
*************************************************/

#include "TH1D.h"
#include "TLine.h"
#include <vector>
#include <cmath>
#include "THelpers_2.h"

TCanvas* plotHistswSomeUncAndAbsUncertainty(std::vector<TH1D*> fit_func_list, TH1D* baseHist, TString title="", TString xTitle="", TString yTitle="", TH1D* pt_uncertainty_band=0) {

  int numberOfFitFuncs = fit_func_list.size();
  if (numberOfFitFuncs>5) {
    cout << "Too many histograms for numerator (currently only supports up to 5)" << endl;
    exit;
  }
  if (!baseHist) {
    cout << "baseHist provided does not exist" << endl;
    exit;
  }

  //*************************************************
  // Variables
  bool topPlotLogY = 0;      // 0 = no log; 1= log
  TString yTitle2 = "ratio"; // bottom plot y axis title

  std::vector<int> histColors; 
  histColors.push_back(kBlue);  // change colors as you like
  histColors.push_back(kRed);
  histColors.push_back(kGreen-1);
  histColors.push_back(8);
  histColors.push_back(33);
  //  histColors.push_back(kGreen-1);

  std::vector<int> fillStyles; 
  fillStyles.push_back(3002);
  fillStyles.push_back(3004);
  //  fillStyles.push_back(3305);
  //  fillStyles.push_back(3490);
  fillStyles.push_back(3004);
  fillStyles.push_back(3005);
  fillStyles.push_back(3006);

  std::vector<int> good_markers; 
   good_markers.push_back(20);
   good_markers.push_back(25);
   good_markers.push_back(21);
   good_markers.push_back(24);
   good_markers.push_back(22);
   good_markers.push_back(26);
   good_markers.push_back(29);
   good_markers.push_back(30);
   good_markers.push_back(23);
   good_markers.push_back(28);
   good_markers.push_back(34);


  int histDenominatorColor = kBlack;

  float defaultRatioYmin = 1.02;
  float defaultRatioYmax = 0.60;
  // END of Variables
  //*************************************************

  setTDRStyle();
  TCanvas* c = new TCanvas("c","",600,600);

  TH1D* baseHistogram = (TH1D*)baseHist->Clone();


  // Create ratio histograms
  std::vector<TH1D*> hists_with_fitted_uncertainties;
  for (int i=0; i<numberOfFitFuncs; i++) {
    ////    hists_with_fitted_uncertainties.push_back( (TH1D*)fit_func_list.at(i)->Clone() );
//    hists_with_fitted_uncertainties.push_back( (TH1D*)baseHist->Clone() );
//    hists_with_fitted_uncertainties[i]->GetTitle();   
//    for(unsigned int bin_i=1;bin_i<= hists_with_fitted_uncertainties[i]->GetNbinsX();bin_i++){
//      //Get uncertainty from fitted function
//      Double_t x_bc=baseHist->GetBinCenter(bin_i);
//      Double_t uncertainty=fit_func_list.at(i)->Eval(x_bc);
//      hists_with_fitted_uncertainties[i]->SetBinError(bin_i,uncertainty);
//    }
    hists_with_fitted_uncertainties.push_back( (TH1D*)baseHist->Clone() );
    hists_with_fitted_uncertainties[i]->GetTitle();   
    for(unsigned int bin_i=1;bin_i<= hists_with_fitted_uncertainties[i]->GetNbinsX();bin_i++){
      Double_t uncertainty=fit_func_list.at(i)->GetBinContent(bin_i);
      hists_with_fitted_uncertainties[i]->SetBinError(bin_i,uncertainty);
    }
  }

  for(unsigned int bin_i=1;bin_i<= baseHist->GetNbinsX();bin_i++){
    Double_t uncertainty= pt_uncertainty_band->GetBinError(bin_i)*baseHist->GetBinContent(bin_i);
    pt_uncertainty_band->SetBinError(bin_i,uncertainty);
    pt_uncertainty_band->SetBinContent(bin_i,baseHist->GetBinContent(bin_i));
  }
 
  TH1D* abs_uncertainty_otf = (TH1D*)baseHist->Clone();
  for(unsigned int bin_i=1;bin_i<= abs_uncertainty_otf->GetNbinsX();bin_i++){
    Double_t sum_x2=0;
    for (int i=0; i<numberOfFitFuncs; i++) {
      sum_x2+=TMath::Power(hists_with_fitted_uncertainties[i]->GetBinError(bin_i),2);
    }
    abs_uncertainty_otf->SetBinError(bin_i,TMath::Sqrt(sum_x2));
  }
  abs_uncertainty_otf->SetFillStyle(3001);
  abs_uncertainty_otf->SetFillColor(kGray);
  abs_uncertainty_otf->SetMarkerStyle(1);
 

  baseHistogram->SetLineWidth(2);
  baseHistogram->SetLineColor(histDenominatorColor);
  baseHistogram->SetMarkerColor(histDenominatorColor);
  baseHistogram->Draw();
  abs_uncertainty_otf->Draw("same e3");
  pt_uncertainty_band->Draw("same e3");
  baseHistogram->Draw("same");
  //      baseHistogram->SetLabelSize(0.0);
  //      //  baseHistogram->GetYaxis()->SetNdivisions(10);
//      baseHistogram->GetYaxis()->SetRangeUser(0.95,1.25);
//      baseHistogram->GetXaxis()->SetTitleSize(0.00);
//      baseHistogram->GetYaxis()->SetLabelSize(0.07);
//      baseHistogram->GetYaxis()->SetTitleSize(0.08);
//      baseHistogram->GetYaxis()->SetTitleOffset(0.76);
  baseHistogram->SetTitle(title+";"+xTitle+";"+yTitle);
  //      
  for (int i=0; i<numberOfFitFuncs; i++) {
    hists_with_fitted_uncertainties[i]->SetLineWidth(2);
    hists_with_fitted_uncertainties[i]->SetLineColor(histColors[i]);
    hists_with_fitted_uncertainties[i]->SetFillColor(histColors[i]);
    hists_with_fitted_uncertainties[i]->SetFillStyle(fillStyles[i]);
    hists_with_fitted_uncertainties[i]->SetMarkerColor(histColors[i]);
    //    hists_with_fitted_uncertainties[i]->SetMarkerStyle(good_markers[i+1]);
    hists_with_fitted_uncertainties[i]->SetMarkerStyle(1);
    hists_with_fitted_uncertainties[i]->Draw("same e3");
  }
  c->SetLogy(topPlotLogY);
  return c;


}


