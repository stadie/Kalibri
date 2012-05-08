#include "MakeDateDir.h"
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
  gStyle->SetPadRightMargin(0.16);

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

void plotMPF_vs_relresponse(TString root_file, TString label_for_filename, TString algo) 
{
  setStyle();

  TFile* fone = TFile::Open(root_file);
    if (fone->IsZombie()) {
       cout << "Error opening file" << endl;
       exit(-1);
    }
    else cout << "ROOT-Datei erfolgreich geladen. " << endl;


  TH1D* MPF_data= (TH1D*)fone->Get("MPFVsEta_all_pt/MPFVsEta_all_pt_MPFResponseVsEta_data_L2L3_Jet2Pt0_Mean");
  TH1D* MPF_mc= (TH1D*)fone->Get("MPFVsEta_all_pt/MPFVsEta_all_pt_MPFResponseVsEta_MC_L2L3_Jet2Pt0_Mean");
  TH1D* Asymmetry_data= (TH1D*)fone->Get("AsymmetryVsEta_all_pt/AsymmetryVsEta_all_pt_AsymmetryVsEta_data_L2L3_MeanPt0_RatioOfMeans");
  TH1D* Asymmetry_mc= (TH1D*)fone->Get("AsymmetryVsEta_all_pt/AsymmetryVsEta_all_pt_AsymmetryVsEta_MC_L2L3_MeanPt0_RatioOfMeans");

  TH1D* MPF_spectrum_data= (TH1D*)fone->Get("MPFVsEta_all_pt/MPFVsEta_all_pt_data_EtaSpectrum");
  TH1D* MPF_spectrum_mc= (TH1D*)fone->Get("MPFVsEta_all_pt/MPFVsEta_all_pt_MC_EtaSpectrum");


  TLatex latex;
  latex.SetNDC();
  latex.SetTextSize(0.05);
  latex.SetTextAlign(13);  //align at top


  TCanvas* c= new TCanvas("c","",500,500); 
  c->SetRightMargin(0.02);
  c->SetTopMargin(0.13);
  //c->SetLogx();
  //  c->SetGridy(); 
  //  c->SetGridx();
  TPad *c1_1 = new TPad("c1_1", "newpad",0.01,0.0.01,0.99,0.32);
  c1_1->Draw();
  c1_1->cd();
  c1_1->SetTopMargin(0.01);
  c1_1->SetBottomMargin(0.3);
  c1_1->SetRightMargin(0.1);
  c1_1->SetFillStyle(0);



  TH1D* hresidual_MPF = (TH1D*) MPF_data->Clone();
  hresidual_MPF->SetTitle("");
  hresidual_MPF->SetYTitle("X_{data}/X_{MC}");
  hresidual_MPF->SetXTitle("#eta");
  hresidual_MPF->Divide(MPF_mc,MPF_data);
  //void Divide(const TH1* h1, const TH1* h2, Double_t c1 = 1, Double_t c2 = 1, Option_t* option = "")

  TH1D* hresidual_Asymmetry = (TH1D*) Asymmetry_data->Clone();
  hresidual_Asymmetry->SetTitle("");
  hresidual_Asymmetry->SetYTitle("Residual correction");
  hresidual_Asymmetry->SetXTitle("#eta");
  hresidual_Asymmetry->Divide(Asymmetry_mc,Asymmetry_data);
  //void Divide(const TH1* h1, const TH1* h2, Double_t c1 = 1, Double_t c2 = 1, Option_t* option = "")


  TH1D* hresidual_MPF_Asymmetry_ratio = (TH1D*) MPF_data->Clone();
  hresidual_MPF_Asymmetry_ratio->Divide(hresidual_MPF,hresidual_Asymmetry);

  hresidual_MPF_Asymmetry_ratio->Draw();
  hresidual_MPF_Asymmetry_ratio->SetLineWidth(1);
  hresidual_MPF_Asymmetry_ratio->SetTitle(";#eta;ratio");
  hresidual_MPF_Asymmetry_ratio->GetXaxis()->SetTitleSize(0.14);
  hresidual_MPF_Asymmetry_ratio->GetXaxis()->SetLabelSize(0.14);
  hresidual_MPF_Asymmetry_ratio->GetYaxis()->SetLabelSize(0.11);
  hresidual_MPF_Asymmetry_ratio->GetYaxis()->SetTitleSize(0.14);
  hresidual_MPF_Asymmetry_ratio->GetYaxis()->SetTitleOffset(0.6);
  hresidual_MPF_Asymmetry_ratio->GetYaxis()->SetNdivisions(5);

  hresidual_MPF_Asymmetry_ratio->GetYaxis()->SetRangeUser(0.95,1.05);
  TLine *line_eta = new TLine(hresidual_MPF->GetXaxis()->GetXmin(),1.,hresidual_MPF->GetXaxis()->GetXmax(),1.);
  line_eta->SetLineStyle(2);
  line_eta->SetLineColor(1);
  line_eta->Draw();

  c->cd();
  c1_2 = new TPad("c1_2", "newpad",0.01,0.33,0.99,0.99);
  c1_2->Draw(); 
  c1_2->cd();
  c1_2->SetTopMargin(0.1);
  c1_2->SetBottomMargin(0.01);
  c1_2->SetRightMargin(0.1);
  c1_1->SetFillStyle(0);

  hresidual_MPF->SetLineColor(1);
  hresidual_MPF->SetMarkerColor(1);
  hresidual_MPF->SetMarkerStyle(20);
  hresidual_Asymmetry->SetLineColor(2);
  hresidual_Asymmetry->SetMarkerColor(2);
  hresidual_Asymmetry->SetMarkerStyle(21);
  hresidual_MPF->GetYaxis()->SetRangeUser(0.9,1.2);
  hresidual_MPF->Draw("C E0 X0 L");
  //  hresidual_->GetXaxis()->SetRange(-1,25);
  hresidual_Asymmetry->Draw("C E0 X0 L SAME"); 
  TPad *p1 = new TPad("i1", "i1",0.8, 0.71,0.9,0.84);
  p1->SetFillStyle(4000);  
  TImage* img = TImage::Open("../kalibriLogoSmall.gif");
  //img->Scale(img->GetWidth(),img->GetHeight()); 
  p1->cd();
  img->Draw("XZ");
  c1_2->cd();
  p1->DrawClone();  
  int nEntries =2;
  TLegend * leg = leg = new TLegend(0.4,0.85-nEntries*0.07,0.8,0.85);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextFont(42);
  leg->SetTextSize(0.04);

  leg->AddEntry(hresidual_MPF,"MPF response","PL");
  leg->AddEntry(hresidual_Asymmetry,"relative response","PL");
  leg->Draw();

  latex.DrawLatex(.25,.8,algo+" Jets");
  
  line_eta->Draw();
  
  c->SaveAs(GetDateDir()+"/MPF_asymm_comp_"+label_for_filename+"_"+ algo + ".eps");

  TCanvas* c2= new TCanvas("c2","",500,500); 
  c2->SetRightMargin(0.02);
  c2->SetTopMargin(0.13);
 c2->SetLogy();

  TH1D* MPF_spectrum_data= (TH1D*)fone->Get("MPFVsEta_all_pt/MPFVsEta_all_pt_data_EtaSpectrum");
  TH1D* MPF_spectrum_mc= (TH1D*)fone->Get("MPFVsEta_all_pt/MPFVsEta_all_pt_MC_EtaSpectrum");


  MPF_spectrum_data->Draw("C E0 X0 L");
  MPF_spectrum_mc->Draw("same hist");

  //  p1->DrawClone();  

  TLegend * leg2 = leg = new TLegend(0.55,0.45-nEntries*0.07,0.8,0.45);
  leg2->SetBorderSize(0);
  leg2->SetFillColor(0);
  leg2->SetTextFont(42);
  leg2->SetTextSize(0.04);
     
  leg2->AddEntry(MPF_spectrum_data,"data","P");
  leg2->AddEntry(MPF_spectrum_mc,"MC","L");
  leg2->Draw();

  latex.DrawLatex(.35,.4,algo+" Jets");
  c2->SaveAs(GetDateDir()+"/MPF_asymm_comp_"+label_for_filename+"_Eta_Spectrum_"+ algo + ".eps");

}
