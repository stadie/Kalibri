#include "TStyle.h"
#include "TFile.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TImage.h"
#include "TPad.h"
#include "TCanvas.h"
#include "TLine.h"

#include <iostream>

void setStyle()
{
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
  gStyle->SetPadRightMargin(0.18);

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

void OverlayPlots(const char* mc, const char* data, const char* hnbeg, const char* hnend, const char* cor1, const char* cor2, double min, double max,const char* lab1,const char* lab2) 
{

  TString cname("c");
  static int nc = 0;
  cname += nc;
  cname +=".";
  TCanvas* c = new TCanvas(cname,"",500,500);
  ++nc;
  c->SetRightMargin(0.04);
  c->SetTopMargin(0.13);
  
  TPad *p1 = new TPad("i1", "i1",0.79, 0.71,0.94,0.84);
  p1->SetFillStyle(4000);  
  TImage* img = TImage::Open("kalibriLogoSmall.gif");
  p1->cd();
  img->Draw("XZ");
 
  TFile* fmc = TFile::Open(mc);
  TFile* fdata = TFile::Open(data);
  
  if(! fmc) {
    std::cout << "file " << mc << " not found!\n";
    return;
  }
  if(! fdata) {
    std::cout << "file " << data << " not found!\n";
    fmc->Close();
    return;
  }
  TString hist1(hnbeg);
  hist1 += "_";
  TString hist2(hist1);
  hist1 += cor1;
  hist2 += cor2;
  hist1 += "_";
  hist2 += "_";
  hist1 += hnend;
  hist2 += hnend;
  TH1D* hunmc  = (TH1D*)fmc->Get(hist1);
  TH1D* hcormc = (TH1D*)fmc->Get(hist2);
  TH1D* hundata  = (TH1D*)fdata->Get(hist1);
  TH1D* hcordata = (TH1D*)fdata->Get(hist2);

    
  TLine *line = new TLine(hunmc->GetXaxis()->GetXmin(),1.0,hunmc->GetXaxis()->GetXmax(),1.0);
  line->SetLineStyle(2);
  line->SetLineColor(1);

  c->cd();
  hunmc->Draw("HIST");
  hunmc->SetMinimum(min);
  hunmc->SetMaximum(max);
  hcormc->Draw("HISTSAME");  
  hcormc->SetMinimum(min);
  hcormc->SetMaximum(max);
  hundata->Draw("SAME");
  hcordata->Draw("SAME");
  
  c->cd();
  //p1->DrawClone();  
  int nEntries = 4;
  TLegend * leg = new TLegend(0.5,0.85-nEntries*0.07,0.8,0.85);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextFont(42);
  leg->SetTextSize(0.04);

  TString label1(lab1);
  TString label2(lab2);
  leg->AddEntry(hundata,label1,"P");
  leg->AddEntry(hcordata,label2,"P");
  leg->AddEntry(hundata,label1+" (MC)","L");
  leg->AddEntry(hcordata,label2+" (MC)","L");

  leg->Draw();
  if(hist1.Contains("VsPt")) c->SetLogx();

  line->Draw();
  c->Print(0,"pdf");
  delete c;
  fmc->Close();
  fdata->Close();
}



void Overlay() 
{
  setStyle();
  /*
  const char* f1 = "/afs/naf.desy.de/user/s/stadie/scratch/dijetsMCSpring10/plots/KalibriPlots.root";
  const char* f2 = "/afs/naf.desy.de/user/s/stadie/scratch/dijetsRun2010A-DCS-RES/plots/KalibriPlots.root";

  const char* cor1 = "Uncorrected";
  const char* cor2 = "L2L3";
  const char* lab1 = "raw";
  const char* lab2 = "L2L3";
  */
  const char* f1 = "/afs/naf.desy.de/user/s/stadie/scratch/dijetsMCSpring10-CorTest/plots/KalibriPlots.root";
  const char* f2 = "/afs/naf.desy.de/user/s/stadie/scratch/dijetsRun2010A-DCS-CorRes/plots/KalibriPlots.root";

  const char* cor1 = "Uncorrected";
  const char* cor2 = "Kalibri";
  const char* lab1 = "L2L3";
  const char* lab2 = "L2L3L4JW";
  
  OverlayPlots(f1,f2,"AsymmetryVsPt/AsymmetryVsPt_AsymmetryVsMeanPt","AbsEta0_RatioOfMeans",cor1,cor2,0.7,1.7,lab1,lab2);
  OverlayPlots(f1,f2,"AsymmetryVsPt/AsymmetryVsPt_AsymmetryVsMeanPt","AbsEta1_RatioOfMeans",cor1,cor2,0.7,1.7,lab1,lab2);
  OverlayPlots(f1,f2,"AsymmetryVsPt/AsymmetryVsPt_AsymmetryVsMeanPt","AbsEta2_RatioOfMeans",cor1,cor2,0.5,2.5,lab1,lab2);
  OverlayPlots(f1,f2,"AsymmetryVsPt/AsymmetryVsPt_AsymmetryVsMeanPt","AbsEta3_RatioOfMeans",cor1,cor2,0.5,2.5,lab1,lab2);
 
  OverlayPlots(f1,f2,"AsymmetryVsEta/AsymmetryVsEta_AsymmetryVsEta","MeanPt2_RatioOfMeans",cor1,cor2,0.5,2.2,lab1,lab2);  
  OverlayPlots(f1,f2,"AsymmetryVsEta/AsymmetryVsEta_AsymmetryVsEta","MeanPt3_RatioOfMeans",cor1,cor2,0.5,2.2,lab1,lab2);  
  OverlayPlots(f1,f2,"AsymmetryVsEta/AsymmetryVsEta_AsymmetryVsEta","MeanPt4_RatioOfMeans",cor1,cor2,0.5,2.2,lab1,lab2);  
  OverlayPlots(f1,f2,"AsymmetryVsEta/AsymmetryVsEta_AsymmetryVsEta","MeanPt5_RatioOfMeans",cor1,cor2,0.5,2.2,lab1,lab2);  
  OverlayPlots(f1,f2,"AsymmetryVsEta/AsymmetryVsEta_AsymmetryVsEta","MeanPt6_RatioOfMeans",cor1,cor2,0.5,2.2,lab1,lab2);   
  OverlayPlots(f1,f2,"AsymmetryVsEta/AsymmetryVsEta_AsymmetryVsEta","MeanPt7_RatioOfMeans",cor1,cor2,0.5,2.2,lab1,lab2);   
  OverlayPlots(f1,f2,"AsymmetryVsEta/AsymmetryVsEta_AsymmetryVsEta","MeanPt8_RatioOfMeans",cor1,cor2,0.5,2.2,lab1,lab2);   

  OverlayPlots(f1,f2,"AsymmetryVsPt/AsymmetryVsPt_AsymmetryVsMeanPt","AbsEta0_StandardDeviation",cor1,cor2,0,0.3,lab1,lab2);
  OverlayPlots(f1,f2,"AsymmetryVsPt/AsymmetryVsPt_AsymmetryVsMeanPt","AbsEta1_StandardDeviation",cor1,cor2,0,0.3,lab1,lab2);
  OverlayPlots(f1,f2,"AsymmetryVsPt/AsymmetryVsPt_AsymmetryVsMeanPt","AbsEta2_StandardDeviation",cor1,cor2,0,0.3,lab1,lab2);
  OverlayPlots(f1,f2,"AsymmetryVsPt/AsymmetryVsPt_AsymmetryVsMeanPt","AbsEta3_StandardDeviation",cor1,cor2,0,0.3,lab1,lab2); 

  OverlayPlots(f1,f2,"AsymmetryVsMeanMoment/AsymmetryVsMeanMoment_AsymmetryVsmeanMoment","AbsEta0_RatioOfMeans",cor1,cor2,0.5,2.0,lab1,lab2);
  OverlayPlots(f1,f2,"AsymmetryVsMeanMoment/AsymmetryVsMeanMoment_AsymmetryVsmeanMoment","AbsEta1_RatioOfMeans",cor1,cor2,0.5,2.0,lab1,lab2);
  OverlayPlots(f1,f2,"AsymmetryVsMeanMoment/AsymmetryVsMeanMoment_AsymmetryVsmeanMoment","AbsEta2_RatioOfMeans",cor1,cor2,0.5,2.5,lab1,lab2);
  OverlayPlots(f1,f2,"AsymmetryVsMeanMoment/AsymmetryVsMeanMoment_AsymmetryVsmeanMoment","AbsEta3_RatioOfMeans",cor1,cor2,0.5,2.5,lab1,lab2);
  
  OverlayPlots(f1,f2,"AsymmetryVsMeanMomentMeanPt/AsymmetryVsMeanMomentMeanPt_AsymmetryVsmeanMoment","MeanPt2_RatioOfMeans",cor1,cor2,0.5,2.5,lab1,lab2);
  OverlayPlots(f1,f2,"AsymmetryVsMeanMomentMeanPt/AsymmetryVsMeanMomentMeanPt_AsymmetryVsmeanMoment","MeanPt3_RatioOfMeans",cor1,cor2,0.5,2.5,lab1,lab2);
  OverlayPlots(f1,f2,"AsymmetryVsMeanMomentMeanPt/AsymmetryVsMeanMomentMeanPt_AsymmetryVsmeanMoment","MeanPt4_RatioOfMeans",cor1,cor2,0.5,2.5,lab1,lab2);
  OverlayPlots(f1,f2,"AsymmetryVsMeanMomentMeanPt/AsymmetryVsMeanMomentMeanPt_AsymmetryVsmeanMoment","MeanPt5_RatioOfMeans",cor1,cor2,0.5,2.5,lab1,lab2);
  OverlayPlots(f1,f2,"AsymmetryVsMeanMomentMeanPt/AsymmetryVsMeanMomentMeanPt_AsymmetryVsmeanMoment","MeanPt6_RatioOfMeans",cor1,cor2,0.5,2.5,lab1,lab2); 
  OverlayPlots(f1,f2,"AsymmetryVsMeanMomentMeanPt/AsymmetryVsMeanMomentMeanPt_AsymmetryVsmeanMoment","MeanPt7_RatioOfMeans",cor1,cor2,0.5,2.5,lab1,lab2); 
  OverlayPlots(f1,f2,"AsymmetryVsMeanMomentMeanPt/AsymmetryVsMeanMomentMeanPt_AsymmetryVsmeanMoment","MeanPt8_RatioOfMeans",cor1,cor2,0.5,2.5,lab1,lab2); 

  OverlayPlots(f1,f2,"AsymmetryVsMeanMomentMeanPt2/AsymmetryVsMeanMomentMeanPt2_AsymmetryVsmeanMoment","MeanPt2_RatioOfMeans",cor1,cor2,0.5,2.5,lab1,lab2);
  OverlayPlots(f1,f2,"AsymmetryVsMeanMomentMeanPt2/AsymmetryVsMeanMomentMeanPt2_AsymmetryVsmeanMoment","MeanPt3_RatioOfMeans",cor1,cor2,0.5,2.5,lab1,lab2);
  OverlayPlots(f1,f2,"AsymmetryVsMeanMomentMeanPt2/AsymmetryVsMeanMomentMeanPt2_AsymmetryVsmeanMoment","MeanPt4_RatioOfMeans",cor1,cor2,0.5,2.5,lab1,lab2);
  OverlayPlots(f1,f2,"AsymmetryVsMeanMomentMeanPt2/AsymmetryVsMeanMomentMeanPt2_AsymmetryVsmeanMoment","MeanPt5_RatioOfMeans",cor1,cor2,0.5,2.5,lab1,lab2);
  OverlayPlots(f1,f2,"AsymmetryVsMeanMomentMeanPt2/AsymmetryVsMeanMomentMeanPt2_AsymmetryVsmeanMoment","MeanPt6_RatioOfMeans",cor1,cor2,0.5,2.5,lab1,lab2);
  OverlayPlots(f1,f2,"AsymmetryVsMeanMomentMeanPt2/AsymmetryVsMeanMomentMeanPt2_AsymmetryVsmeanMoment","MeanPt7_RatioOfMeans",cor1,cor2,0.5,2.5,lab1,lab2);
  OverlayPlots(f1,f2,"AsymmetryVsMeanMomentMeanPt2/AsymmetryVsMeanMomentMeanPt2_AsymmetryVsmeanMoment","MeanPt8_RatioOfMeans",cor1,cor2,0.5,2.5,lab1,lab2);
}
