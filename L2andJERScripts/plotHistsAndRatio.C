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
   root [2] .L plotHistsAndRatio.C
   root [3] plotHistsAndRatio(photonEt_BarrelPIDcut,photonEt_Barrel,"Photon E_{T}","E_{T} (GeV)","entries/20 GeV bin")  

  Michael Anderson
  March 18, 2009
*************************************************/

#include "TH1D.h"
#include "TLine.h"
#include <vector>
#include <cmath>
#include "THelpers_2.h"

TCanvas* plotHistsAndRatio(std::vector<TH1D*> numeratorHistograms, TH1D* denominatorHist, TString title="", TString xTitle="", TString yTitle="", Bool_t plot_ratio=true) {

  int numberOfNumeratorHists = numeratorHistograms.size();
  if (numberOfNumeratorHists>5) {
    cout << "Too many histograms for numerator (currently only supports up to 5)" << endl;
    exit;
  }
  if (!denominatorHist) {
    cout << "denominatorHist provided does not exist" << endl;
    exit;
  }

  //*************************************************
  // Variables
  bool topPlotLogY = 0;      // 0 = no log; 1= log
  TString yTitle2 = "ratio"; // bottom plot y axis title
  //2010
  //  Double_t y_min=0.85;
  //  Double_t y_max=1.15;
  //2011
  Double_t y_min=0.95;
  Double_t y_max=1.25;


  std::vector<int> histColors; 
  histColors.push_back(kBlue);  // change colors as you like
  histColors.push_back(kRed);
  histColors.push_back(kGreen-1);
  histColors.push_back(8);
  histColors.push_back(33);
  //  histColors.push_back(kGreen-1);

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

  TCanvas *c1 = new TCanvas("c1", "c1",0,0,600,500);
  c1->Range(0,0,1,1);

  std::vector<TH1D*> hists;
  for (int i=0; i<numberOfNumeratorHists; i++) {
    hists.push_back( (TH1D*)numeratorHistograms[i] );
  }
  TH1D* denominatorHistogram = (TH1D*)denominatorHist->Clone();


  // Create ratio histograms
  std::vector<TH1D*> hist_over_denomHist;
  for (int i=0; i<numberOfNumeratorHists; i++) {
    hist_over_denomHist.push_back( (TH1D*)numeratorHistograms[i]->Clone() );
    hist_over_denomHist[i]->GetTitle();   
    hist_over_denomHist[i]->Divide(denominatorHistogram);
  }

  if(!plot_ratio)
    {

      setTDRStyle();
      TCanvas* c = new TCanvas("c","",600,600);


      denominatorHistogram->SetLineWidth(2);
      denominatorHistogram->SetLineColor(histDenominatorColor);
      denominatorHistogram->SetMarkerColor(histDenominatorColor);
      denominatorHistogram->Draw();
//      denominatorHistogram->SetLabelSize(0.0);
//      //  denominatorHistogram->GetYaxis()->SetNdivisions(10);
      denominatorHistogram->GetYaxis()->SetRangeUser(y_min,y_max);
//      denominatorHistogram->GetXaxis()->SetTitleSize(0.00);
//      denominatorHistogram->GetYaxis()->SetLabelSize(0.07);
//      denominatorHistogram->GetYaxis()->SetTitleSize(0.08);
//      denominatorHistogram->GetYaxis()->SetTitleOffset(0.76);
      denominatorHistogram->SetTitle(title+";"+xTitle+";"+yTitle);
//      
      for (int i=0; i<numberOfNumeratorHists; i++) {
	hists[i]->SetLineWidth(2);
	hists[i]->SetLineColor(histColors[i]);
	hists[i]->SetMarkerColor(histColors[i]);
	hists[i]->SetMarkerStyle(good_markers[i+1]);
	hists[i]->Draw("same");
      }

      c->SetLogy(topPlotLogY);
      return c;
    }

  c1->cd();

  //*************************************************
  // Bottom plot
  TPad *c1_1 = new TPad("c1_1", "newpad",0.01,0.01,0.99,0.32);
  c1_1->Draw();
  c1_1->cd();
  c1_1->SetTopMargin(0.01);
  c1_1->SetBottomMargin(0.3);
  c1_1->SetRightMargin(0.1);
  c1_1->SetFillStyle(0);


  hist_over_denomHist[0]->Draw();
  drawUncertainty(1.0, 0.02, hist_over_denomHist[0]->GetXaxis()->GetXmin() , hist_over_denomHist[0]->GetXaxis()->GetXmax());
  drawUncertainty(1.0, 0.005, hist_over_denomHist[0]->GetXaxis()->GetXmin() , hist_over_denomHist[0]->GetXaxis()->GetXmax(),41);
  hist_over_denomHist[0]->Draw("same");
  hist_over_denomHist[0]->SetLineWidth(1);
  hist_over_denomHist[0]->SetLineColor(histColors[0]);
  hist_over_denomHist[0]->SetMarkerColor(histColors[0]);
  hist_over_denomHist[0]->SetMarkerStyle(good_markers[0]);
  hist_over_denomHist[0]->SetMinimum(defaultRatioYmin);
  hist_over_denomHist[0]->SetMaximum(defaultRatioYmax);
  hist_over_denomHist[0]->GetYaxis()->SetNdivisions(5);
  hist_over_denomHist[0]->GetYaxis()->SetRangeUser(0.95,1.05);
  hist_over_denomHist[0]->SetTitle(";"+xTitle+";"+yTitle2);
  hist_over_denomHist[0]->GetXaxis()->SetTitleSize(0.14);
  hist_over_denomHist[0]->GetXaxis()->SetLabelSize(0.14);
  hist_over_denomHist[0]->GetYaxis()->SetLabelSize(0.11);
  hist_over_denomHist[0]->GetYaxis()->SetTitleSize(0.14);
  hist_over_denomHist[0]->GetYaxis()->SetTitleOffset(0.48);
  for (int i=1; i<numberOfNumeratorHists; i++) {
    hist_over_denomHist[i]->SetLineWidth(1);
    hist_over_denomHist[i]->SetLineColor(histColors[i]);
    hist_over_denomHist[i]->SetMarkerColor(histColors[i]);
    hist_over_denomHist[i]->SetMarkerStyle(good_markers[i]);
    hist_over_denomHist[i]->Draw("same");
  }
  TLine *line_eta = new TLine(hist_over_denomHist[0]->GetXaxis()->GetXmin(),1.,hist_over_denomHist[0]->GetXaxis()->GetXmax(),1.);
  line_eta->SetLineStyle(2);
  line_eta->SetLineColor(1);
  line_eta->Draw();

  // End bottom plot
  //*************************************************
 

  //*************************************************
  // Top Plot
  c1->cd();
  TPad *c1_2 = new TPad("c1_2", "newpad",0.01,0.33,0.99,0.99);
  c1_2->Draw(); 
  c1_2->cd();
  c1_2->SetTopMargin(0.1);
  c1_2->SetBottomMargin(0.01);
  c1_2->SetRightMargin(0.1);
  c1_1->SetFillStyle(0);


  denominatorHistogram->SetLineWidth(2);
  denominatorHistogram->SetLineColor(histDenominatorColor);
  denominatorHistogram->SetMarkerColor(histDenominatorColor);
  denominatorHistogram->Draw();
  denominatorHistogram->SetLabelSize(0.0);
  //  denominatorHistogram->GetYaxis()->SetNdivisions(10);
  denominatorHistogram->GetYaxis()->SetRangeUser(y_min,y_max);
  denominatorHistogram->GetXaxis()->SetTitleSize(0.00);
  denominatorHistogram->GetYaxis()->SetLabelSize(0.07);
  denominatorHistogram->GetYaxis()->SetTitleSize(0.08);
  denominatorHistogram->GetYaxis()->SetTitleOffset(0.76);
  denominatorHistogram->SetTitle(title+";;"+yTitle);

  for (int i=0; i<numberOfNumeratorHists; i++) {
    hists[i]->SetLineWidth(2);
    hists[i]->SetLineColor(histColors[i]);
    hists[i]->SetMarkerColor(histColors[i]);
    hists[i]->SetMarkerStyle(good_markers[i]);
    hists[i]->Draw("same");
  }

  c1_2->SetLogy(topPlotLogY);
  // End bottom plot
  //*************************************************

  return c1;
}



TCanvas* plotHistsAndRatio(TH1D* numeratorHist, TH1D* denominatorHist, TString title="", TString xTitle="", TString yTitle="", Bool_t plot_ratio=true) {  
  std::vector<TH1D*> numeratorHistograms;
  numeratorHistograms.push_back( numeratorHist );
  TCanvas* c1 = plotHistsAndRatio(numeratorHistograms, denominatorHist, title, xTitle, yTitle);
  return c1;
}

TCanvas* plotHistsAndRatio(TH1D* numeratorHist1, TH1D* numeratorHist2, TH1D* denominatorHist, TString title="", TString xTitle="", TString yTitle="", Bool_t plot_ratio=true) {
  std::vector<TH1D*> numeratorHistograms;
  numeratorHistograms.push_back( numeratorHist1 );
  numeratorHistograms.push_back( numeratorHist2 );
  TCanvas* c1 = plotHistsAndRatio(numeratorHistograms, denominatorHist, title, xTitle, yTitle);
  return c1;
}

TCanvas* plotHistsAndRatio(TH1D* numeratorHist1, TH1D* numeratorHist2, TH1D* numeratorHist3, TH1D* denominatorHist, TString title="", TString xTitle="", TString yTitle="", Bool_t plot_ratio=true) {
  std::vector<TH1D*> numeratorHistograms;
  numeratorHistograms.push_back( numeratorHist1 );
  numeratorHistograms.push_back( numeratorHist2 );
  numeratorHistograms.push_back( numeratorHist3 );
  TCanvas *c1 = plotHistsAndRatio(numeratorHistograms, denominatorHist, title, xTitle, yTitle);
  return c1;
}
