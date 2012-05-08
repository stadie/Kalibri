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
   root [2] .L plotRatios.C
   root [3] plotRatios(photonEt_BarrelPIDcut,photonEt_Barrel,"Photon E_{T}","E_{T} (GeV)","entries/20 GeV bin")  

  Michael Anderson
  March 18, 2009
*************************************************/

#include "TH1D.h"
#include "TF1.h"
#include "TLine.h"
#include "MakeDateDir.h"
#include <vector>
#include <cmath>
#include "THelpers_2.h"

TH1D* get_absolute_uncertainty_band_from_histos(TH1D* draw_hist, std::vector<TF1*> functions_, std::vector<TH1D*> hist_over_denomHist){

  //  Int_t no_x_bins = 100;
  Int_t no_x_bins = draw_hist->GetNbinsX();

  TH1D* abs_uncertainty_band = (TH1D*) draw_hist->Clone();
  //new TH1D("abs_uncertainty_band", "", no_x_bins, draw_hist->GetXaxis()->GetXmin(), draw_hist->GetXaxis()->GetXmax());

    abs_uncertainty_band->SetName("abs_uncertainty_band");
    abs_uncertainty_band->SetTitle("abs_uncertainty_band");
    //new TH1D("res_hist_pt_"+(Long_t)pt, "res_hist_pt_"+(Long_t)pt,no_x_bins, trad_res_hist->GetXaxis()->GetXmin() , trad_res_hist->GetXaxis()->GetXmax());

    for(Int_t i=1;i< no_x_bins+1;i++){
      Double_t BinCenter = abs_uncertainty_band->GetXaxis()->GetBinCenter(i);
      Double_t sum_x2 = 0;
      for(unsigned int f_i=0;f_i<functions_.size();f_i++){
	//	if(functions_.at(f_i)->GetName()!="NOFIT")sum_x2+= TMath::Abs(TMath::Power(functions_.at(f_i)->Eval(BinCenter),2));
	if(functions_.at(f_i)->GetName()!="NOFIT")sum_x2+=TMath::Abs(TMath::Power(hist_over_denomHist.at(f_i)->GetBinContent(i),2));
      }
      abs_uncertainty_band->SetBinContent(i,0);
      abs_uncertainty_band->SetBinError(i,TMath::Sqrt(sum_x2));
    }
    return abs_uncertainty_band;
}



TCanvas* plotRatios_wAbsUncFit(std::vector<TH1D*> numeratorHistograms, TH1D* denominatorHist, TString title="", TString xTitle="", TString yTitle="", TH1D* pt_uncertainty_band=0) {

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
  bool showHists = 1;     
  TString yTitle2 = "|C_{X}/C_{base}-1|"; // bottom plot y axis title

  std::vector<int> histColors; 
  histColors.push_back(kBlue);  // change colors as you like
  histColors.push_back(kRed);
  histColors.push_back(kGreen-1);
  histColors.push_back(8);
  histColors.push_back(33);
  //  histColors.push_back(kGreen-1);

  std::vector<int> fillStyles; 
  fillStyles.push_back(3305);
  fillStyles.push_back(3350);
  fillStyles.push_back(3004);
  fillStyles.push_back(3005);
  fillStyles.push_back(3006);


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


  TF1 *systematic_fit = new TF1("systematic_fit","[0]+[1]*TMath::CosH(x)+[2]*TMath::CosH(x)",0,10);
  TF1 *systematic_fit_2 = new TF1("systematic_fit_2","[0]+[1]*TMath::CosH(x)+[2]*TMath::CosH(x)",0,10);
  TF1 *kFSR_fit = new TF1("kFSR_fit","(x<2.9)*([0]+[1]*cosh(x)/(1+cosh(x)*[2]))+(x>2.9)*[3]",0,10); //was used before...
  kFSR_fit->SetParameters(0.9,0.014,0.26,0.01);
  kFSR_fit->SetParName(0,"const");
  kFSR_fit->SetParName(1,"par1");
  kFSR_fit->SetParName(2,"par2");
  kFSR_fit->SetParName(3,"const2");
  TF1 *kFSR_fit_2 = new TF1("kFSR_fit_2","(x<2.9)*([0]+[1]*cosh(x)/(1+cosh(x)*[2]))+(x>=2.9)*[3]",0,10); //was used before...
  kFSR_fit_2->SetParameters(0.9,0.014,0.26,0.01);
  TF1 *errorfct_fit = new TF1("errorfct_fit","0.5*[2]*(TMath::Erf([0]*(x-[1]))+1)",0,10);
  TF1 *errorfct_fit_2 = new TF1("errorfct_fit_2","0.5*[2]*(TMath::Erf([0]*(x-[1]))+1)",0,10);
  TF1 *pol5_fit = new TF1("pol5_fit","pol5",0,10);
  TF1 *pol5_fit_2 = new TF1("pol5_fit_2","pol5",0,10);
  TF1 *NOFIT = new TF1("NOFIT","[0]",0,10);
  NOFIT->SetParameters(0.1,0.1);
  std::vector<TString> fit_functions_names_;
  TString fit_name;
  // Create ratio histograms
  std::vector<TH1D*> hist_over_denomHist;
  for (int i=0; i<numberOfNumeratorHists; i++) {
    hist_over_denomHist.push_back( (TH1D*)numeratorHistograms[i]->Clone() );
    hist_over_denomHist[i]->GetTitle();   
    hist_over_denomHist[i]->Divide(denominatorHistogram);
    for(int bin_i = 1; bin_i <= hist_over_denomHist[i]->GetNbinsX(); bin_i++) {
      //HACK!!! IN ORDER TO GET REASONABLE (AND CONSERVATIVE) FITS
      hist_over_denomHist[i]->SetBinContent(bin_i,TMath::Abs(hist_over_denomHist[i]->GetBinContent(bin_i)-1));
    }
    TString plotname(hist_over_denomHist[i]->GetName());
    if(plotname.Contains("JER-syst")){
      fit_functions_names_.push_back("kFSR_fit");
      fit_name=fit_functions_names_.back();
    }
    else if(plotname.Contains("cons-triggers")){
      fit_functions_names_.push_back("systematic_fit");
      fit_name=fit_functions_names_.back();
    }
    else if(plotname.Contains("Radiation")){
      //      fit_functions_names_.push_back("pol5_fit");
            fit_functions_names_.push_back("kFSR_fit");
      //      fit_functions_names_.push_back("errorfct_fit");
      fit_name=fit_functions_names_.back();
    }
    else if(plotname.Contains("Time")){
      fit_functions_names_.push_back("NOFIT");
      fit_name=fit_functions_names_.back();
    }
    else if(plotname.Contains("MPF")){
      fit_functions_names_.push_back("pol5_fit");
      //      fit_functions_names_.push_back("kFSR_fit");
      //      fit_functions_names_.push_back("errorfct_fit");
      fit_name=fit_functions_names_.back();
    }
    else if(plotname.Contains("RR")){
      //fit_functions_names_.push_back("NOFIT");
      //      fit_functions_names_.push_back("pol5_fit");
            fit_functions_names_.push_back("kFSR_fit");
      //      fit_functions_names_.push_back("errorfct_fit");
      fit_name=fit_functions_names_.back();
    }
    else{
      fit_functions_names_.push_back("kFSR_fit");
      //    fit_functions_names_.push_back("systematic_fit");
    fit_name=fit_functions_names_.back();
    }

    if(fit_name!="NOFIT"){
    hist_over_denomHist[i]->Fit(fit_name,"","same");
    hist_over_denomHist[i]->Fit(fit_name+"_2","+","same");
    //    hist_over_denomHist[i]->GetFunction(fit_name)->SetFillStyle(fillStyles.at(i));
    //    hist_over_denomHist[i]->GetFunction(fit_name)->SetFillColor(histColors.at(i));
    hist_over_denomHist[i]->GetFunction(fit_name)->SetLineColor(histColors.at(i));
    hist_over_denomHist[i]->GetFunction(fit_name)->SetLineWidth(2);
    //    hist_over_denomHist[i]->GetFunction(fit_name+"_2")->SetFillStyle(fillStyles.at(i));
    //    hist_over_denomHist[i]->GetFunction(fit_name+"_2")->SetFillColor(histColors.at(i));
    hist_over_denomHist[i]->GetFunction(fit_name+"_2")->SetLineColor(histColors.at(i));
    hist_over_denomHist[i]->GetFunction(fit_name+"_2")->SetLineWidth(2);
       if(fit_name=="systematic_fit"){
    hist_over_denomHist[i]->GetFunction(fit_name+"_2")->SetParameter(0,-hist_over_denomHist[i]->GetFunction(fit_name+"_2")->GetParameter(0));
    hist_over_denomHist[i]->GetFunction(fit_name+"_2")->SetParameter(1,-hist_over_denomHist[i]->GetFunction(fit_name+"_2")->GetParameter(1));
    hist_over_denomHist[i]->GetFunction(fit_name+"_2")->SetParameter(2,-hist_over_denomHist[i]->GetFunction(fit_name+"_2")->GetParameter(2));
       }
    }
    else{
      hist_over_denomHist[i]->SetFillStyle(fillStyles.at(i));
      hist_over_denomHist[i]->SetFillColor(histColors.at(i));
    }
  }


  std::vector<TF1*> fit_functions_;
  cout << "works here.." << endl;
  for (int i=0; i<numberOfNumeratorHists; i++) {
    if(hist_over_denomHist[i]->GetFunction(fit_functions_names_.at(i))){
	fit_functions_.push_back(hist_over_denomHist[i]->GetFunction(fit_functions_names_.at(i)));
      }
      else{
	fit_functions_.push_back(NOFIT);
      }
  }
  cout << "works here.." << fit_functions_.size() <<endl;
  TH1D* abs_uncertainty_band = get_absolute_uncertainty_band_from_histos(hist_over_denomHist[0],fit_functions_,hist_over_denomHist);
  cout << "works here.." << endl;
  abs_uncertainty_band->SetFillStyle(3001);
  abs_uncertainty_band->SetFillColor(kGray);
  abs_uncertainty_band->SetMarkerStyle(1);

//  std::vector<TF1*> fits_hist_over_denomHist_;
//  for (int i=0; i<numberOfNumeratorHists; i++) {
//      fits_hist_over_denomHist_.push_back(new TF1("fa1","sin(x)/x",0,10));
//      //   
//
//    
//  }

  //*************************************************
  // Bottom plot
//  TPad *c1_1 = new TPad("c1_1", "newpad",0.01,0.01,0.99,0.32);
//  c1_1->Draw();
//  c1_1->cd();
//  c1_1->SetTopMargin(0.01);
//  c1_1->SetBottomMargin(0.3);
//  c1_1->SetRightMargin(0.1);
//  c1_1->SetFillStyle(0);
//

  /*if(showHists)*/hist_over_denomHist[0]->Draw("hist");
  //  else hist_over_denomHist[0]->GetFunction(fit_functions_names_.at(0))->Draw();
  drawUncertainty(0.0, 0.02, hist_over_denomHist[0]->GetXaxis()->GetXmin() , hist_over_denomHist[0]->GetXaxis()->GetXmax());
  drawUncertainty(0.0, 0.005, hist_over_denomHist[0]->GetXaxis()->GetXmin() , hist_over_denomHist[0]->GetXaxis()->GetXmax(),41);
  abs_uncertainty_band->Draw("same E3");
  if(pt_uncertainty_band){
    cout << "found uncertainty band..." << endl;
    pt_uncertainty_band->Draw("same E3");
  }
  if(fit_functions_names_.at(0)=="NOFIT"){hist_over_denomHist[0]->Draw("same hist");
  abs_uncertainty_band->Draw("same E3");
  }
  else if(showHists)hist_over_denomHist[0]->Draw("same");
  else hist_over_denomHist[0]->GetFunction(fit_functions_names_.at(0))->Draw("same");
  hist_over_denomHist[0]->SetLineWidth(1);
  hist_over_denomHist[0]->SetLineColor(histColors[0]);
  hist_over_denomHist[0]->SetMarkerColor(histColors[0]);
  hist_over_denomHist[0]->SetMinimum(defaultRatioYmin);
  hist_over_denomHist[0]->SetMaximum(defaultRatioYmax);
  //  hist_over_denomHist[0]->GetYaxis()->SetNdivisions(5);
  hist_over_denomHist[0]->GetYaxis()->SetRangeUser(0.0,0.15);
  hist_over_denomHist[0]->SetTitle(";"+xTitle+";"+yTitle2);
//  hist_over_denomHist[0]->GetXaxis()->SetTitleSize(0.14);
//  hist_over_denomHist[0]->GetXaxis()->SetLabelSize(0.14);
//  hist_over_denomHist[0]->GetYaxis()->SetLabelSize(0.11);
//  hist_over_denomHist[0]->GetYaxis()->SetTitleSize(0.14);
//  hist_over_denomHist[0]->GetYaxis()->SetTitleOffset(0.48);
  for (int i=1; i<numberOfNumeratorHists; i++) {
    hist_over_denomHist[i]->SetLineWidth(1);
    hist_over_denomHist[i]->SetLineColor(histColors[i]);
    hist_over_denomHist[i]->SetMarkerColor(histColors[i]);
    if(showHists)hist_over_denomHist[i]->Draw("same");
    else if(fit_functions_names_.at(i)=="NOFIT")hist_over_denomHist[i]->Draw("same hist");
    else hist_over_denomHist[i]->GetFunction(fit_functions_names_.at(i))->Draw("same");
  }
  TLine *line_eta = new TLine(hist_over_denomHist[0]->GetXaxis()->GetXmin(),0.,hist_over_denomHist[0]->GetXaxis()->GetXmax(),0.);
  line_eta->SetLineStyle(2);
  line_eta->SetLineColor(1);
  line_eta->Draw();


  TFile *outf = new TFile(GetDateDir()+"/"+((TString)"Systematic_fit_funcs_HIST_ABSUNC_rel_to_"+denominatorHist->GetName())+".root","RECREATE");

  for(unsigned int i=0; i< fit_functions_names_.size();i++){
    TF1* tempptr=  hist_over_denomHist[i]->GetFunction(fit_functions_names_.at(i));
    tempptr->SetName((TString)"Fit_func_"+hist_over_denomHist[i]->GetName());
    tempptr->Write();
    hist_over_denomHist[i]->GetYaxis()->SetRangeUser(0.0,0.05);
    hist_over_denomHist[i]->Write();
    //    abs_uncertainty_band->SetName(abs_uncertainty_band->GetName()+(TString)"_"+hist_over_denomHist[i]->GetName());
  }

  abs_uncertainty_band->Write();


  outf->Close();

  return c1;
}



TCanvas* plotRatios_wAbsUncFit(TH1D* numeratorHist, TH1D* denominatorHist, TString title="", TString xTitle="", TString yTitle="") {  
  std::vector<TH1D*> numeratorHistograms;
  numeratorHistograms.push_back( numeratorHist );
  TCanvas* c1 = plotRatios(numeratorHistograms, denominatorHist, title, xTitle, yTitle);
  return c1;
}

TCanvas* plotRatios_wAbsUncFit(TH1D* numeratorHist1, TH1D* numeratorHist2, TH1D* denominatorHist, TString title="", TString xTitle="", TString yTitle="") {
  std::vector<TH1D*> numeratorHistograms;
  numeratorHistograms.push_back( numeratorHist1 );
  numeratorHistograms.push_back( numeratorHist2 );
  TCanvas* c1 = plotRatios(numeratorHistograms, denominatorHist, title, xTitle, yTitle);
  return c1;
}

TCanvas* plotRatios_wAbsUncFit(TH1D* numeratorHist1, TH1D* numeratorHist2, TH1D* numeratorHist3, TH1D* denominatorHist, TString title="", TString xTitle="", TString yTitle="") {
  std::vector<TH1D*> numeratorHistograms;
  numeratorHistograms.push_back( numeratorHist1 );
  numeratorHistograms.push_back( numeratorHist2 );
  numeratorHistograms.push_back( numeratorHist3 );
  TCanvas *c1 = plotRatios(numeratorHistograms, denominatorHist, title, xTitle, yTitle);
  return c1;
}
