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
  gStyle->SetTitleY(0.999);
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


void TH1style(TH1D* histo, std::vector<Int_t> line_styles_, std::vector<Int_t> colours_, std::vector<Int_t> markers_, Int_t line_style, Int_t colour, Int_t marker, TString ytitle)
{
  histo->SetMarkerStyle(markers_[marker]);
  histo->SetMarkerColor(colours_[colour]);
  histo->SetLineColor(colours_[colour]);
  histo->SetLineStyle(line_styles_[line_style]);
  histo->GetYaxis()->SetTitle(ytitle);

}


void TH1style_plus_fit(TH1D* histo, std::vector<Int_t> line_styles_, std::vector<Int_t> colours_, std::vector<Int_t> markers_, Int_t line_style, Int_t colour, Int_t marker, TString ytitle, TString fit_func_name)
{
  histo->SetMarkerStyle(markers_[marker]);
  histo->SetMarkerColor(colours_[colour]);
  histo->SetLineColor(colours_[colour]);
  histo->SetLineStyle(line_styles_[line_style]);
  histo->GetYaxis()->SetTitle(ytitle);
  if( histo->GetFunction(fit_func_name)){
    histo->GetFunction(fit_func_name)->SetLineColor(colours_[colour]);
    histo->GetFunction(fit_func_name)->SetLineStyle(line_styles_[line_style]);
  }
}

void TGraphErrorsstyle(TGraphErrors* histo, std::vector<Int_t> line_styles_, std::vector<Int_t> colours_, std::vector<Int_t> markers_, Int_t line_style, Int_t colour, Int_t marker, TString ytitle)
{
  histo->SetMarkerStyle(markers_[marker]);
  histo->SetMarkerColor(colours_[colour]);
  histo->SetLineColor(colours_[colour]);
  histo->SetLineStyle(line_styles_[line_style]);
  histo->GetYaxis()->SetTitle(ytitle);

}


void pm_eta_TGraph_draw(TCanvas* draw_canvas, TGraphErrors* all_eta, TGraphErrors* abs_eta,std::vector<Int_t> line_styles_, std::vector<Int_t> colours_, std::vector<Int_t> markers_)
{

  Int_t no_abseta_points=abs_eta->GetN();
  if(all_eta->GetN()!=2*no_abseta_points)cout << "da laeuft was schief alleta N: "<< all_eta->GetN() << " and " <<no_abseta_points  << endl;

  Double_t* abs_eta_x=abs_eta->GetX();
  Double_t* abs_eta_ex=abs_eta->GetEX();
  Double_t* all_eta_y=all_eta->GetY();
  Double_t* all_eta_ey=all_eta->GetEY();

  //  Double_t* all_eta_x=all_eta->GetX();
  //  Double_t* all_eta_ex=all_eta->GetEX();


  std::vector < Double_t > plus_eta_y;
  std::vector < Double_t > plus_eta_ey;

  std::vector < Double_t > minus_eta_y;
  std::vector < Double_t > minus_eta_ey;
  for(Int_t array_i=0;array_i<no_abseta_points;array_i++){
    plus_eta_y.push_back(all_eta_y[no_abseta_points+array_i]);
    plus_eta_ey.push_back(all_eta_ey[no_abseta_points+array_i]);
    minus_eta_y.push_back(all_eta_y[no_abseta_points-(array_i+1)]);
    minus_eta_ey.push_back(all_eta_ey[no_abseta_points-(array_i+1)]);
    //    cout << "Abseta eta: " << abs_eta_x[array_i] << " alleta plus_eta: " << all_eta_x[no_abseta_points+array_i] << " alleta minus_eta: " << all_eta_x[no_abseta_points-(array_i+1)] << endl;
  }

TGraphErrors *plus_eta = new TGraphErrors(no_abseta_points,abs_eta_x,&plus_eta_y[0],abs_eta_ex,&plus_eta_ey[0]);
TGraphErrors *minus_eta = new TGraphErrors(no_abseta_points,abs_eta_x,&minus_eta_y[0],abs_eta_ex,&minus_eta_ey[0]);
 TGraphErrorsstyle(abs_eta,line_styles_, colours_, markers_, 0,0,0,"bla");
 TGraphErrorsstyle(plus_eta,line_styles_, colours_, markers_, 0,1,1,"bla");
 TGraphErrorsstyle(minus_eta,line_styles_, colours_, markers_, 0,2,2,"bla");

  abs_eta->Draw("ALP");
  abs_eta->GetXaxis()->SetTitle("|#eta|");
  abs_eta->GetYaxis()->SetTitle("Correction factor");
  abs_eta->GetYaxis()->SetRangeUser(0.8,1.2);
  //p1->DrawClone();  
   cout << "was ist hier...?" << endl;
   plus_eta->Draw("P same");
   minus_eta->Draw("P same");

   TLegend *leg_selected;
   leg_selected = new TLegend(0.2,0.70,0.55,0.85);
   leg_selected->SetFillColor(kWhite);
   //   leg->SetHeader("Legende");
   leg_selected->AddEntry(abs_eta,"(|#eta|)","lep");
   leg_selected->AddEntry(plus_eta,"(+#eta)","lep");
   leg_selected->AddEntry(minus_eta,"(-#eta)","lep");

   leg_selected->Draw();
   draw_canvas->Print((TString)"PM_eta_draw"+abs_eta->GetName()+(TString) ".eps");


}


void pm_eta_TGraph_draw(TCanvas* draw_canvas, TGraphErrors* all_eta, TGraphErrors* abs_eta, TH1D* abs_eta_histo, std::vector<Int_t> line_styles_, std::vector<Int_t> colours_, std::vector<Int_t> markers_, TString ytitle, Double_t y_low=0.8, Double_t y_hig=1.2)
{

  Int_t no_abseta_points=abs_eta->GetN();
  if(all_eta->GetN()!=2*no_abseta_points)cout << "da laeuft was schief alleta N: "<< all_eta->GetN() << " and " <<no_abseta_points  << endl;

  Double_t* abs_eta_x=abs_eta->GetX();
  Double_t* abs_eta_ex=abs_eta->GetEX();
  Double_t* all_eta_y=all_eta->GetY();
  Double_t* all_eta_ey=all_eta->GetEY();

  //  Double_t* all_eta_x=all_eta->GetX();
  //  Double_t* all_eta_ex=all_eta->GetEX();


  std::vector < Double_t > plus_eta_y;
  std::vector < Double_t > plus_eta_ey;

  std::vector < Double_t > minus_eta_y;
  std::vector < Double_t > minus_eta_ey;

  std::vector < Double_t > diff_eta_y;
  std::vector < Double_t > diff_eta_ey;
  std::vector < Double_t > diff_eta_offset_y;
  std::vector < Double_t > diff_eta_offset_ey;
  for(Int_t array_i=0;array_i<no_abseta_points;array_i++){
    plus_eta_y.push_back(all_eta_y[no_abseta_points+array_i]);
    plus_eta_ey.push_back(all_eta_ey[no_abseta_points+array_i]);
    minus_eta_y.push_back(all_eta_y[no_abseta_points-(array_i+1)]);
    minus_eta_ey.push_back(all_eta_ey[no_abseta_points-(array_i+1)]);
    diff_eta_y.push_back(plus_eta_y.back()-minus_eta_y.back());
    diff_eta_ey.push_back(TMath::Sqrt(TMath::Power(plus_eta_ey.back(),2)+TMath::Power(minus_eta_ey.back(),2)));
    diff_eta_offset_y.push_back(plus_eta_y.back()-minus_eta_y.back()+y_low);
    diff_eta_offset_ey.push_back(TMath::Sqrt(TMath::Power(plus_eta_ey.back(),2)+TMath::Power(minus_eta_ey.back(),2)));
    //    cout << "Abseta eta: " << abs_eta_x[array_i] << " alleta plus_eta: " << all_eta_x[no_abseta_points+array_i] << " alleta minus_eta: " << all_eta_x[no_abseta_points-(array_i+1)] << endl;
  }

TGraphErrors *plus_eta = new TGraphErrors(no_abseta_points,abs_eta_x,&plus_eta_y[0],abs_eta_ex,&plus_eta_ey[0]);
TGraphErrors *minus_eta = new TGraphErrors(no_abseta_points,abs_eta_x,&minus_eta_y[0],abs_eta_ex,&minus_eta_ey[0]);
TGraphErrors *diff_eta = new TGraphErrors(no_abseta_points,abs_eta_x,&diff_eta_y[0],abs_eta_ex,&diff_eta_ey[0]);
TGraphErrors *diff_offset_eta = new TGraphErrors(no_abseta_points,abs_eta_x,&diff_eta_offset_y[0],abs_eta_ex,&diff_eta_offset_ey[0]);
 TGraphErrorsstyle(abs_eta,line_styles_, colours_, markers_, 0,0,0,"bla");
 TGraphErrorsstyle(plus_eta,line_styles_, colours_, markers_, 0,1,1,"bla");
 TGraphErrorsstyle(minus_eta,line_styles_, colours_, markers_, 0,2,2,"bla");

 TGraphErrorsstyle(diff_eta,line_styles_, colours_, markers_, 4,4,4,"bla");
 TGraphErrorsstyle(diff_offset_eta,line_styles_, colours_, markers_, 4,4,4,"bla");

 TH1style(abs_eta_histo,line_styles_, colours_, markers_, 0,0,0,"bla");
  abs_eta_histo->Draw("hist");
  abs_eta_histo->GetXaxis()->SetTitle("|#eta|");
  abs_eta_histo->GetYaxis()->SetTitle(ytitle);
  abs_eta_histo->GetYaxis()->SetRangeUser(y_low,y_hig);
  //p1->DrawClone();  
   cout << "was ist hier...?" << endl;
   plus_eta->Draw("P same");
   minus_eta->Draw("P same");

   TLegend *leg_selected;
   leg_selected = new TLegend(0.2,0.70,0.55,0.85);
   leg_selected->SetFillColor(kWhite);
   //   leg->SetHeader("Legende");
   leg_selected->AddEntry(abs_eta,"(|#eta|)","le");
   leg_selected->AddEntry(plus_eta,"(+#eta)","lep");
   leg_selected->AddEntry(minus_eta,"(-#eta)","lep");

   leg_selected->Draw();
   draw_canvas->Print((TString)"PM_eta_draw_hist"+abs_eta->GetName()+(TString) ".eps");

   TLine *line_offset = new TLine(abs_eta_x[0],y_low,abs_eta_x[no_abseta_points-1]+2,y_low);
   line_offset->SetLineStyle(2);
   abs_eta_histo->GetYaxis()->SetRangeUser(y_low-0.1,y_hig);
   diff_offset_eta->Draw("P same");
   line_offset->Draw();
   TPaveText pt(.01,.25,.1,.35,"NDC ARC br");
   pt.SetFillColor(kWhite);
   pt.SetBorderSize(5);
   pt.SetLineColor(kBlack);
   //   pt.SetLabel("Diff:");
   pt.AddText("Diff:");
   pt.Draw();
   draw_canvas->Print((TString)"PM_eta_draw_hist_overlay_diff_"+abs_eta->GetName()+(TString) ".eps");

   TLine *line = new TLine(abs_eta_x[0],0.0,abs_eta_x[no_abseta_points-1]+2,0.0);
   line->SetLineStyle(2);

   diff_eta->Draw("ALP");
   diff_eta->GetYaxis()->SetRangeUser(-0.025,+0.025);
   line->Draw();
   draw_canvas->Print((TString)"PM_eta_draw_hist_difference"+abs_eta->GetName()+(TString) ".eps");

}

