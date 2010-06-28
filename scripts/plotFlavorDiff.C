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

void plotFlavorDiff() 
{
  setStyle();
  TH1D* hglu1= (TH1D*)gDirectory->Get("MCTruthRespFlavorVsGenJetPt/MCTruthRespFlavorVsGenJetPt_GenJetResponseVsGenJetPt_Uncorrected_Flavor0_GaussFitMean");
  TH1D* huds1= (TH1D*)gDirectory->Get("MCTruthRespFlavorVsGenJetPt/MCTruthRespFlavorVsGenJetPt_GenJetResponseVsGenJetPt_Uncorrected_Flavor1_GaussFitMean");
  TH1D* hglu2= (TH1D*)gDirectory->Get("MCTruthRespFlavorVsGenJetPt/MCTruthRespFlavorVsGenJetPt_GenJetResponseVsGenJetPt_Kalibri_Flavor0_GaussFitMean");
  TH1D* huds2= (TH1D*)gDirectory->Get("MCTruthRespFlavorVsGenJetPt/MCTruthRespFlavorVsGenJetPt_GenJetResponseVsGenJetPt_Kalibri_Flavor1_GaussFitMean");

  TCanvas* c= new TCanvas("c","",500,500); 
  c->SetRightMargin(0.02);
  c->SetTopMargin(0.13);
  //c->SetLogx();
  c->SetGridy(); 
  c->SetGridx();
  TH1D* hdiff1 = (TH1D*) huds1->Clone();
  hdiff1->SetTitle("");
  hdiff1->SetYTitle("R_{uds}-R_{gluon}");
  hdiff1->SetXTitle("p_{T}^{gen} (GeV)");
  hdiff1->Add(hglu1,-1);
  TH1D* hdiff2 = (TH1D*) huds2->Clone();
  hdiff2->Add(hglu2,-1);
  hdiff1->SetMaximum(0.15);
  hdiff1->SetMinimum(0);
  hdiff1->SetBinContent(1,10);
  hdiff1->SetBinContent(2,10);
  hdiff1->SetBinContent(3,10);
  hdiff1->SetBinContent(4,10);  
  hdiff1->SetBinContent(5,10);  
  hdiff2->SetBinContent(1,10);
  hdiff2->SetBinContent(2,10);
  hdiff2->SetBinContent(3,10);
  hdiff2->SetBinContent(4,10);
  hdiff2->SetBinContent(5,10);

  hdiff1->Draw("C E0 X0 L");
  hdiff1->GetXaxis()->SetRange(-1,25);
  hdiff2->Draw("C E0 X0 L SAME"); 
  TPad *p1 = new TPad("i1", "i1",0.79, 0.71,0.94,0.84);
  p1->SetFillStyle(4000);  
  TImage* img = TImage::Open("kalibriLogoSmall.gif");
  //img->Scale(img->GetWidth(),img->GetHeight()); 
  p1->cd();
  img->Draw("XZ");
  c->cd();
  p1->DrawClone();  
  int nEntries =2;
  TLegend * leg = leg = new TLegend(0.5,0.85-nEntries*0.07,0.8,0.85);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextFont(42);
  leg->SetTextSize(0.04);

  leg->AddEntry(hdiff1,"L2L3","PL");
  leg->AddEntry(hdiff2,"Kalibri","PL");
  leg->Draw();
  TCanvas* c2= new TCanvas("c2","",500,500);
  TH1D* hres = (TH1D*)gDirectory->Get("MCTruthResponseVsGenJetPt/MCTruthResponseVsGenJetPt_GenJetResponseVsGenJetPt_Uncorrected_Eta2_GaussFitWidth"); 
  TH1D* hres2 = (TH1D*)gDirectory->Get("MCTruthResponseVsGenJetPt/MCTruthResponseVsGenJetPt_GenJetResponseVsGenJetPt_Kalibri_Eta2_GaussFitWidth");
  hres->SetYTitle("#sigma/p_{T}");
  hres->SetXTitle("p_{T}^{gen} (GeV)");
  TF1 resol("resol","sqrt([0]*[0]/x + [1] * [1])");
  resol.SetParameters(1.2,0.05);
  //hres->Fit("resol");
  //hres2->Fit("resol");  
  c2->SetGridy(); 
  c2->SetGridx();  
  c2->SetRightMargin(0.02);
  c2->SetTopMargin(0.13);
  hres->SetFillColor(hres->GetMarkerColor());
  hres->Draw("PE0 X0");
  hres->SetMaximum(0.3);
  hres2->SetMaximum(0.3);
  hres->GetXaxis()->SetRange(-1,25);
  hres2->SetFillColor(hres2->GetMarkerColor());
  hres2->Draw("PE0 X0 SAME");
  p1->DrawClone();
  leg->DrawClone();


  TCanvas* c3= new TCanvas("c3","",500,500);

  TH1D* rel_hres = hres->Clone();
  TH1D* rel_hres2 = hres2->Clone();
  rel_hres->Divide(hres,hres);
  rel_hres2->Divide(hres2,hres);
  //void Divide(const TH1* h1, const TH1* h2, Double_t  c1 = 1, Double_t  c2 = 1, Option_t* option = "")

  rel_hres->SetYTitle("(#sigma/p_{T})/(#sigma_{Uncorrected}/p_{T})");
  rel_hres->SetXTitle("p_{T}^{gen} (GeV)");

  c3->SetGridy(); 
  c3->SetGridx();  
  c3->SetRightMargin(0.02);
  c3->SetTopMargin(0.13);
  rel_hres->SetFillColor(rel_hres->GetMarkerColor());
  rel_hres->Draw("PE0 X0");
  rel_hres->GetYaxis()->SetRangeUser(0.6,1.2);
  rel_hres->SetMaximum(1.2);
  rel_hres2->SetMaximum(1.1);
  rel_hres->GetXaxis()->SetRange(-1,25);
  rel_hres2->SetFillColor(rel_hres2->GetMarkerColor());
  rel_hres2->Draw("PE0 X0 SAME");
  p1->DrawClone();
  leg->DrawClone();

}
