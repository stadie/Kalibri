#include "TStyle.h"
#include "TFile.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TImage.h"
#include "TPad.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TF1.h"

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

void OverlayPlots(const char* mc, const char* data, const char* hn, double min, double max,const char* lab1,const char* lab2) 
{
  TF1 resolmc("resolmc","sqrt(sign([0])*[0]*[0]/x/x+[1]*[1]*pow(x,[3]-1)+[2]*[2])",3); 
  resolmc.SetParameters(1.0,0.7,0,0.05);
  
  TF1 resoldata("resoldata","sqrt(sign([0])*[0]*[0]/x/x+[1]*[1]*pow(x,[3]-1)+[2]*[2])",3); 
  resoldata.SetParameters(1.0,0.7,0,0.05);
 

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
  TString hist1(hn);
  TH1D* hunmc  =   (TH1D*)fmc->Get(hist1);
  TH1D* hundata  = (TH1D*)fdata->Get(hist1);

  hundata->SetMarkerStyle(3);
  hundata->SetMarkerColor(2);
  hundata->SetLineColor(2);

  TLine *line = new TLine(hunmc->GetXaxis()->GetXmin(),1.0,hunmc->GetXaxis()->GetXmax(),1.0);
  line->SetLineStyle(2);
  line->SetLineColor(1);

  resolmc.SetLineColor(4);
  resoldata.SetLineColor(2);
  

  c->cd();
  hunmc->Fit("resolmc","","PE");  hunmc->Fit("resolmc","M","PE");
  //hunmc->Draw("PE");
  hunmc->SetMinimum(min);
  hunmc->SetMaximum(max);
  //hundata->Draw("PE SAME");
  hundata->Fit("resoldata","","PE SAME");  hundata->Fit("resoldata","M","PE SAME");
  
  c->cd();
  //p1->DrawClone();  
  int nEntries = 2;
  TLegend * leg = new TLegend(0.3,0.85-nEntries*0.07,0.7,0.85);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextFont(42);
  leg->SetTextSize(0.04);

  TString label1(lab1);
  TString label2(lab2);
  leg->AddEntry(hunmc,label1,"PL");
  leg->AddEntry(hundata,label2,"PL");

  leg->Draw();
  if(hist1.Contains("VsPt")) c->SetLogx();
  if(hist1.Contains("VsGenJetPt")) c->SetLogx();

  //line->Draw();
  p1->Draw("XZ");
  c->Print(0,"eps");
  delete c;
  fmc->Close();
  fdata->Close();
}



void OverlayResol() 
{
  setStyle();

  const char* f1 = "~/scratch/plots/Fall11plotsak5FastPF/KalibriPlots.root";
  const char* f2 = "~/scratch/plots/2012V7plotsak5FastPFS6/KalibriPlots.root";

  const char* lab1 = "2011 Fall11   42x";
  const char* lab2 = "2012 Summer12 52x";
  
  OverlayPlots(f1,f2,"MCTruthResolVsGenJetPt3/MCTruthResolVsGenJetPt3_GenJetResponseVsGenJetPt__L2L3_AbsEta0_GaussFitWidth",0,0.4,lab1,lab2);
  OverlayPlots(f1,f2,"MCTruthResolVsGenJetPt3/MCTruthResolVsGenJetPt3_GenJetResponseVsGenJetPt__L2L3_AbsEta1_GaussFitWidth",0,0.4,lab1,lab2);
  OverlayPlots(f1,f2,"MCTruthResolVsGenJetPt3/MCTruthResolVsGenJetPt3_GenJetResponseVsGenJetPt__L2L3_AbsEta2_GaussFitWidth",0,0.4,lab1,lab2);
  OverlayPlots(f1,f2,"MCTruthResolVsGenJetPt3/MCTruthResolVsGenJetPt3_GenJetResponseVsGenJetPt__L2L3_AbsEta3_GaussFitWidth",0,0.4,lab1,lab2);
  OverlayPlots(f1,f2,"MCTruthResolVsGenJetPt3/MCTruthResolVsGenJetPt3_GenJetResponseVsGenJetPt__L2L3_AbsEta4_GaussFitWidth",0,0.4,lab1,lab2);
  OverlayPlots(f1,f2,"MCTruthResolVsGenJetPt3/MCTruthResolVsGenJetPt3_GenJetResponseVsGenJetPt__L2L3_AbsEta5_GaussFitWidth",0,0.4,lab1,lab2);
  OverlayPlots(f1,f2,"MCTruthResolVsGenJetPt3/MCTruthResolVsGenJetPt3_GenJetResponseVsGenJetPt__L2L3_AbsEta6_GaussFitWidth",0,0.4,lab1,lab2);

  OverlayPlots(f1,f2,"MCTruthResolVsGenJetPt3/MCTruthResolVsGenJetPt3_GenJetResponseVsGenJetPt__L2L3_AbsEta0_StandardDeviation",0,0.4,lab1,lab2);
  OverlayPlots(f1,f2,"MCTruthResolVsGenJetPt3/MCTruthResolVsGenJetPt3_GenJetResponseVsGenJetPt__L2L3_AbsEta1_StandardDeviation",0,0.4,lab1,lab2);
  OverlayPlots(f1,f2,"MCTruthResolVsGenJetPt3/MCTruthResolVsGenJetPt3_GenJetResponseVsGenJetPt__L2L3_AbsEta2_StandardDeviation",0,0.4,lab1,lab2);
  OverlayPlots(f1,f2,"MCTruthResolVsGenJetPt3/MCTruthResolVsGenJetPt3_GenJetResponseVsGenJetPt__L2L3_AbsEta3_StandardDeviation",0,0.4,lab1,lab2);
  OverlayPlots(f1,f2,"MCTruthResolVsGenJetPt3/MCTruthResolVsGenJetPt3_GenJetResponseVsGenJetPt__L2L3_AbsEta4_StandardDeviation",0,0.4,lab1,lab2);
  OverlayPlots(f1,f2,"MCTruthResolVsGenJetPt3/MCTruthResolVsGenJetPt3_GenJetResponseVsGenJetPt__L2L3_AbsEta5_StandardDeviation",0,0.4,lab1,lab2);
  OverlayPlots(f1,f2,"MCTruthResolVsGenJetPt3/MCTruthResolVsGenJetPt3_GenJetResponseVsGenJetPt__L2L3_AbsEta6_StandardDeviation",0,0.4,lab1,lab2);

}
