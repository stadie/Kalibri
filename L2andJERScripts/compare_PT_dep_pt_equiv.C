#include "tdrstyle_mod.C"
#include "MakeDateDir.h"

TH1D* get_res_hist_at_pt(TH1D* trad_res_hist, TH1D* res_CONST_hist, TH1D* res_SLOPE_hist, Double_t rel_pt_dev){

    TF1 *fit_loglin = new TF1("fit_loglin","[0]+[1]*TMath::Log(x)",10,2000); //was used before...
    fit_loglin->SetParameters(1,1);
    fit_loglin->SetParName(0,"const");
    fit_loglin->SetParName(1,"slope");
    Int_t no_x_bins = trad_res_hist->GetNbinsX();

    TH1D* res_hist_at_pt = trad_res_hist->Clone();
    res_hist_at_pt->SetName("res_hist_pt_"+(Long_t)rel_pt_dev);
    res_hist_at_pt->SetTitle("res_hist_pt_"+(Long_t)rel_pt_dev);
    //new TH1D("res_hist_pt_"+(Long_t)pt, "res_hist_pt_"+(Long_t)pt,no_x_bins, trad_res_hist->GetXaxis()->GetXmin() , trad_res_hist->GetXaxis()->GetXmax());

    for(Int_t i=1;i< no_x_bins+1;i++){
      //      cout << "ACHTIUNG" << res_CONST_hist->GetBinContent(i) << endl;
      fit_loglin->SetParameter(0,res_CONST_hist->GetBinContent(i));
      fit_loglin->SetParError(0,res_CONST_hist->GetBinError(i));
      fit_loglin->SetParameter(1,res_SLOPE_hist->GetBinContent(i));
      fit_loglin->SetParError(1,res_SLOPE_hist->GetBinError(i));
      //      fit_loglin->Print();
      Double_t pt = fit_loglin->GetX(trad_res_hist->GetBinContent(i));//,50, 500, 1.E-10, 10000, true);
      //      cout << "x_bin: " << i<< " pt: " <<pt << " y value " << rel_pt_dev/100*trad_res_hist->GetBinContent(i) << " test: " << fit_loglin->Eval(300)<<endl;
      res_hist_at_pt->SetBinContent(i,fit_loglin->Eval(rel_pt_dev/100*pt));
      //      cout << fit_loglin->Eval(pt) << endl;
    }
    return res_hist_at_pt;
}

TH1D* get_central_pt(TH1D* trad_res_hist, TH1D* res_CONST_hist, TH1D* res_SLOPE_hist){

    TF1 *fit_loglin = new TF1("fit_loglin","[0]+[1]*TMath::Log(x)",10,2000); //was used before...
    fit_loglin->SetParameters(1,1);
    fit_loglin->SetParName(0,"const");
    fit_loglin->SetParName(1,"slope");
    Int_t no_x_bins = trad_res_hist->GetNbinsX();

    TH1D* central_pt = trad_res_hist->Clone();
    central_pt->GetYaxis()->SetRangeUser(100,1000);
    central_pt->SetName("central_pt_");
    central_pt->SetTitle("central_pt_");
    //new TH1D("res_hist_pt_"+(Long_t)pt, "res_hist_pt_"+(Long_t)pt,no_x_bins, trad_res_hist->GetXaxis()->GetXmin() , trad_res_hist->GetXaxis()->GetXmax());

    for(Int_t i=1;i< no_x_bins+1;i++){
      //      cout << "ACHTIUNG" << res_CONST_hist->GetBinContent(i) << endl;
      fit_loglin->SetParameter(0,res_CONST_hist->GetBinContent(i));
      fit_loglin->SetParError(0,res_CONST_hist->GetBinError(i));
      fit_loglin->SetParameter(1,res_SLOPE_hist->GetBinContent(i));
      fit_loglin->SetParError(1,res_SLOPE_hist->GetBinError(i));
      //      fit_loglin->Print();
      Double_t pt = fit_loglin->GetX(trad_res_hist->GetBinContent(i));//,50, 500, 1.E-10, 10000, true);
      central_pt->SetBinContent(i,pt);
      //      cout << fit_loglin->Eval(pt) << endl;
    }
    return central_pt;
}

TH1D* get_pseudo_uncertainty_band(TH1D* high_res_hist, TH1D* low_res_hist){

    Int_t no_x_bins = high_res_hist->GetNbinsX();

    TH1D* pt_uncertainty_band = high_res_hist->Clone();
    pt_uncertainty_band->SetName("pt_uncertainty_band");
    pt_uncertainty_band->SetTitle("pt_uncertainty_band");
    //new TH1D("res_hist_pt_"+(Long_t)pt, "res_hist_pt_"+(Long_t)pt,no_x_bins, trad_res_hist->GetXaxis()->GetXmin() , trad_res_hist->GetXaxis()->GetXmax());

    for(Int_t i=0;i< no_x_bins+1;i++){
      //      cout << "ACHTIUNG" << res_CONST_hist->GetBinContent(i) << endl;
      pt_uncertainty_band->SetBinContent(i,(high_res_hist->GetBinContent(i)+low_res_hist->GetBinContent(i))/2);
      pt_uncertainty_band->SetBinError(i,(high_res_hist->GetBinContent(i)-low_res_hist->GetBinContent(i))/2);
    }
    return pt_uncertainty_band;
}


void compare_PT_dep_pt_equiv(TString dir1, TString algo="PF", TString dir_prefix="KOSTAS_L1_on_pt_plain_", TString label_dir1="", TString binning_select ="kostas", Bool_t pub_style=true){
  setTDRStyle();
  TCanvas* c = new TCanvas("c","",600,600);
  //  c->SetRightMargin(0.12);
   TPad *pad1 = new TPad("pad1","",0,0,1,1);
   TPad *pad2 = new TPad("pad2","",0,0,1,1);
   pad1->SetRightMargin(0.15);
   pad2->SetFillStyle(4000); //will be transparent
   pad1->SetTicky(0);
   pad1->Draw();
   pad1->cd();
  TLatex latex;
  latex.SetNDC();
  latex.SetTextSize(0.05);
  latex.SetTextAlign(13);  //align at top

  //  Double_t rel_lower_pt=20;
  //  Double_t rel_upper_pt=500;
  Double_t rel_lower_pt=50;
  Double_t rel_upper_pt=200;

    TFile *inf_Dir1;
    inf_Dir1 = new TFile(dir_prefix+dir1+"/res_corrections_histos.root","OPEN");
    if (inf_Dir1->IsZombie()) {
       cout << "Error opening file" << endl;
       exit(-1);
    }
    else cout << "ROOT-Datei erfolgreich geladen. " << endl;

        if(label_dir1!="")dir1=label_dir1;

    TH1D* import_FSRcorr_residuals_eta_res1_Dir1;
    TH1D* import_FSRcorr_residuals_Abseta_res1_Dir1;
    //    cout << "kFSR_eq_one_"+binning_select+"_"+algo+"_eta_res_hist" << endl;
    import_FSRcorr_residuals_eta_res1_Dir1  = (TH1D*)inf_Dir1->Get("kFSR_eq_one_"+binning_select+"_TuneZ2_"+algo+"_eta_res_hist");
    import_FSRcorr_residuals_eta_CONST_FIT_Dir1  = (TH1D*)inf_Dir1->Get("kFSR_eq_one_"+binning_select+"_TuneZ2_"+algo+"_eta_res_const_slope_hist");
    import_FSRcorr_residuals_eta_SLOPE_FIT_Dir1  = (TH1D*)inf_Dir1->Get("kFSR_eq_one_"+binning_select+"_TuneZ2_"+algo+"_eta_res_slope_hist");

    import_FSRcorr_residuals_Abseta_res1_Dir1  = (TH1D*)inf_Dir1->Get("kFSR_eq_one_"+binning_select+"_TuneZ2_"+algo+"_Abseta_res_hist");
    import_FSRcorr_residuals_Abseta_CONST_FIT_Dir1  = (TH1D*)inf_Dir1->Get("kFSR_eq_one_"+binning_select+"_TuneZ2_"+algo+"_Abseta_res_const_slope_hist");
    import_FSRcorr_residuals_Abseta_SLOPE_FIT_Dir1  = (TH1D*)inf_Dir1->Get("kFSR_eq_one_"+binning_select+"_TuneZ2_"+algo+"_Abseta_res_slope_hist");

    TLine *line_eta = new TLine(import_FSRcorr_residuals_eta_res1_Dir1->GetXaxis()->GetXmin(),1.,import_FSRcorr_residuals_eta_res1_Dir1->GetXaxis()->GetXmax(),1.);
    line_eta->Draw();
    line_eta->SetLineStyle(2);
    line_eta->SetLineColor(1);
    TLine *line_Abseta = new TLine(import_FSRcorr_residuals_Abseta_res1_Dir1->GetXaxis()->GetXmin(),1.,import_FSRcorr_residuals_Abseta_res1_Dir1->GetXaxis()->GetXmax(),1.);
    line_Abseta->Draw();
    line_Abseta->SetLineStyle(2);
    line_Abseta->SetLineColor(1);


    //kFSR_eq_one_k_HFfix_TuneZ2_PF_eta_res_slope_hist

    import_FSRcorr_residuals_eta_res1_Dir1->SetStats(0);
    import_FSRcorr_residuals_eta_res1_Dir1->GetYaxis()->SetTitle("Residual Correction");
    import_FSRcorr_residuals_eta_res1_Dir1->GetXaxis()->SetTitle("#eta");

    import_FSRcorr_residuals_Abseta_res1_Dir1->SetStats(0);
    import_FSRcorr_residuals_Abseta_res1_Dir1->GetYaxis()->SetTitle("Residual Correction");
    import_FSRcorr_residuals_Abseta_res1_Dir1->GetXaxis()->SetTitle("|#eta|");



    //ETA

    import_FSRcorr_residuals_eta_res1_Dir1->SetLineColor(1);
    import_FSRcorr_residuals_eta_res1_Dir1->SetMarkerColor(1);
    import_FSRcorr_residuals_eta_res1_Dir1->SetMarkerStyle(24);
    import_FSRcorr_residuals_eta_res1_Dir1->GetYaxis()->SetRangeUser(0.85,1.2);
    import_FSRcorr_residuals_eta_res1_Dir1->Draw("");

    TH1D* low_pt_eta_hist = get_res_hist_at_pt(import_FSRcorr_residuals_eta_res1_Dir1, import_FSRcorr_residuals_eta_CONST_FIT_Dir1, import_FSRcorr_residuals_eta_SLOPE_FIT_Dir1, rel_lower_pt);
    TH1D* high_pt_eta_hist = get_res_hist_at_pt(import_FSRcorr_residuals_eta_res1_Dir1, import_FSRcorr_residuals_eta_CONST_FIT_Dir1, import_FSRcorr_residuals_eta_SLOPE_FIT_Dir1, rel_upper_pt);
    TH1D* central_pt_eta = get_central_pt(import_FSRcorr_residuals_eta_res1_Dir1, import_FSRcorr_residuals_eta_CONST_FIT_Dir1, import_FSRcorr_residuals_eta_SLOPE_FIT_Dir1);

    low_pt_eta_hist->SetLineColor(2);	 
    low_pt_eta_hist->SetMarkerColor(2); 
    low_pt_eta_hist->SetMarkerStyle(21);
    low_pt_eta_hist->Draw("same");      
    high_pt_eta_hist->SetLineColor(4);	 
    high_pt_eta_hist->SetMarkerColor(4); 
    high_pt_eta_hist->SetMarkerStyle(20);
    high_pt_eta_hist->Draw("same");      

    TH1D* pseudo_uncertainty_eta= get_pseudo_uncertainty_band(high_pt_eta_hist,low_pt_eta_hist);
    pseudo_uncertainty_eta->SetFillColor(5);
    pseudo_uncertainty_eta->SetMarkerColor(5);
    pseudo_uncertainty_eta->SetMarkerStyle(1);
    pseudo_uncertainty_eta->Draw("same e5");
    import_FSRcorr_residuals_eta_res1_Dir1->Draw("same");
    low_pt_eta_hist->Draw("same");      
    high_pt_eta_hist->Draw("same");      



   Double_t ymin = 100;
   Double_t ymax = 500;
   Double_t dy = (ymax-ymin)/0.8; //5+15 per cent margins top and bottom
   Double_t xmin = pseudo_uncertainty_eta->GetXaxis()->GetXmin();
   Double_t xmax = pseudo_uncertainty_eta->GetXaxis()->GetXmax();
   Double_t dx = (xmax-xmin)/0.7; //2* 15 per cent margins left and right
   //   pad2->Range(xmin-0.1*dx,ymin-0.1*dy,xmax+0.1*dx,ymax+0.1*dy);
      pad2->Range(xmin-0.15*dx,ymin-0.15*dy,xmax+0.15*dx,ymax+0.05*dy);
   //   pad2->Range(xmin,ymin,xmax,ymax);
   pad2->SetRightMargin(0.15);
   pad2->SetLeftMargin(0.15);
   pad2->SetTopMargin(0.05);
   pad2->SetBottomMargin(0.15);
   pad2->Draw();
   pad2->cd();
    //    central_pt_Abseta->GetYaxis()->SetRangeUser(10,2000);

   central_pt_eta->SetLineColor(kGreen+3);
   central_pt_eta->Draw("][sames hist");
   //   central_pt_Abseta->Draw("][sames");
   pad2->Update();
   // draw axis on the right side of the pad
   TGaxis *axis = new TGaxis(xmax,ymin,xmax,ymax,ymin,ymax,50510,"+L");
   axis->SetLabelColor(kGreen+3);
   axis->Draw();
   axis->SetTitleColor(kGreen+3);
   axis->SetTitleFont(42);
   axis->SetTitleSize(0.06);
   axis->SetTitleOffset(1.25);
   //   axis->SetTickSize(0.03);
   axis->SetNdivisions(510);
   axis->SetTitle("p_{T} of intersection [GeV]");

   pad1->cd();
 

    TLegend *leg_eta_comb;
    leg_eta_comb = new TLegend(0.25,0.70,0.7,0.85);
    //    leg_eta_comb->SetFillColor(kWhite);
    leg_eta_comb->SetFillStyle(kNone);
    leg_eta_comb->SetTextFont(42);
      //   leg->SetHeader("Legende");
    if(pub_style){
    leg_eta_comb->AddEntry(import_FSRcorr_residuals_eta_res1_Dir1,"2010","p");
    }
    else{
    leg_eta_comb->AddEntry(import_FSRcorr_residuals_eta_res1_Dir1,dir1,"p");
    leg_eta_comb->AddEntry(low_pt_eta_hist,TString("at p_{T}^{is} #upoint")+(Long_t)rel_lower_pt+" %","p");
    leg_eta_comb->AddEntry(high_pt_eta_hist,TString("at p_{T}^{is} #upoint")+(Long_t)rel_upper_pt+" %","p");
    //    leg_eta_comb->AddEntry(high_pt_eta_hist,TString("+ ")+(Long_t)rel_upper_pt+" %","p");
    }
    leg_eta_comb->Draw();
    //"ResComp_"+
    cmsPrel(intLumi=36, false);
    latex.DrawLatex(.25,.9,algo+" Jets");
    line_eta->Draw();

    c->SaveAs(GetDateDir()+"/ResPtDepComp_impr_"+dir1+"_FSRcorr_residuals_eta_"+ algo +".eps");
    central_pt_eta->Draw("][hist");
    central_pt_eta->GetYaxis()->SetTitle("p_{T} of intersection [GeV]");
    central_pt_eta->GetXaxis()->SetTitle("#eta");
    central_pt_eta->GetYaxis()->SetRangeUser(100,400);
    c->SaveAs(GetDateDir()+"/ResPtDepComp_impr_"+dir1+"_FSRcorr_residuals_eta_pt_is_"+ algo +".eps");

  //ABSETA
    //    pad2->Delete();
   c->cd();
   pad1->Draw();
   pad1->cd();
   TPad *pad3 = new TPad("pad3","",0,0,1,1);
   pad3->SetFillStyle(4000); //will be transparent


    import_FSRcorr_residuals_Abseta_res1_Dir1->SetLineColor(1);
    import_FSRcorr_residuals_Abseta_res1_Dir1->SetMarkerColor(1);
    import_FSRcorr_residuals_Abseta_res1_Dir1->SetMarkerStyle(24);
    import_FSRcorr_residuals_Abseta_res1_Dir1->GetYaxis()->SetRangeUser(0.85,1.2);
    import_FSRcorr_residuals_Abseta_res1_Dir1->Draw("");

    TH1D* low_pt_Abseta_hist = get_res_hist_at_pt(import_FSRcorr_residuals_Abseta_res1_Dir1, import_FSRcorr_residuals_Abseta_CONST_FIT_Dir1, import_FSRcorr_residuals_Abseta_SLOPE_FIT_Dir1, rel_lower_pt);
    TH1D* high_pt_Abseta_hist = get_res_hist_at_pt(import_FSRcorr_residuals_Abseta_res1_Dir1, import_FSRcorr_residuals_Abseta_CONST_FIT_Dir1, import_FSRcorr_residuals_Abseta_SLOPE_FIT_Dir1, rel_upper_pt);
    TH1D* central_pt_Abseta = get_central_pt(import_FSRcorr_residuals_Abseta_res1_Dir1, import_FSRcorr_residuals_Abseta_CONST_FIT_Dir1, import_FSRcorr_residuals_Abseta_SLOPE_FIT_Dir1);

    low_pt_Abseta_hist->SetLineColor(2);	 
    low_pt_Abseta_hist->SetMarkerColor(2); 
    low_pt_Abseta_hist->SetMarkerStyle(21);
    low_pt_Abseta_hist->Draw("same");      
    high_pt_Abseta_hist->SetLineColor(4);	 
    high_pt_Abseta_hist->SetMarkerColor(4); 
    high_pt_Abseta_hist->SetMarkerStyle(20);
    high_pt_Abseta_hist->Draw("same");      

    TH1D* pseudo_uncertainty_Abseta= get_pseudo_uncertainty_band(high_pt_Abseta_hist,low_pt_Abseta_hist);
    pseudo_uncertainty_Abseta->SetFillColor(5);
    pseudo_uncertainty_Abseta->SetMarkerColor(5);
    pseudo_uncertainty_Abseta->SetMarkerStyle(1);
    pseudo_uncertainty_Abseta->Draw("same e5");
    import_FSRcorr_residuals_Abseta_res1_Dir1->Draw("same");
    low_pt_Abseta_hist->Draw("same");      
    high_pt_Abseta_hist->Draw("same");      
 
    cout << "hier funzt es noch.. " << endl;

 //scale hint1 to the pad coordinates
//    Float_t rightmax = 500;//1.1*central_pt_Abseta->GetMaximum();
//   Float_t scale = gPad->GetUymax()/rightmax;
//   central_pt_Abseta->SetLineColor(kRed);
//   central_pt_Abseta->Scale(scale);
//   central_pt_Abseta->Draw("same");
//   
//   //draw an axis on the right side
//   TGaxis *axis = new TGaxis(gPad->GetUxmax(),gPad->GetUymin(),
//         gPad->GetUxmax(), gPad->GetUymax(),0,rightmax,510,"+L");
//   axis->SetLineColor(kRed);
//   axis->SetLabelColor(kRed);
//   axis->Draw();
//

//    pad2->cd();
   //Double_t ymin = 100;
   //Double_t ymax = 500;
   //Double_t dy = (ymax-ymin)/0.8; //5+15 per cent margins top and bottom
   //Double_t xmin = pseudo_uncertainty_Abseta->GetXaxis()->GetXmin();
   //Double_t xmax = pseudo_uncertainty_Abseta->GetXaxis()->GetXmax();
   //Double_t dx = (xmax-xmin)/0.7; //2* 15 per cent margins left and right
    xmin = pseudo_uncertainty_Abseta->GetXaxis()->GetXmin();
    xmax = pseudo_uncertainty_Abseta->GetXaxis()->GetXmax();
    dx = (xmax-xmin)/0.7; //2* 15 per cent margins left and right
   //   pad3->Range(xmin-0.1*dx,ymin-0.1*dy,xmax+0.1*dx,ymax+0.1*dy);
   cout << "hier funzt es noch.. " << endl;
     pad3->Range(xmin-0.15*dx,ymin-0.15*dy,xmax+0.15*dx,ymax+0.05*dy);
    cout << "hier funzt es noch.. " << endl;
    //            pad3->Range(xmin,ymin,xmax,ymax);
   pad3->SetRightMargin(0.15);
   pad3->SetLeftMargin(0.15);
   pad3->SetTopMargin(0.05);
   pad3->SetBottomMargin(0.15);
   pad3->Draw();
   pad3->cd();
    //    central_pt_Abseta->GetYaxis()->SetRangeUser(10,2000);

   central_pt_Abseta->SetLineColor(kGreen+3);
   central_pt_Abseta->Draw("][sames hist");
   pad3->Update();
   // draw axis on the right side of the pad
   //   TGaxis *axis = new TGaxis(xmax,ymin,xmax,ymax,ymin,ymax,50510,"+L");
   axis->SetLabelColor(kGreen+3);
   axis->Draw();
   axis->SetTitleColor(kGreen+3);
   axis->SetTitleFont(42);
   axis->SetTitleSize(0.06);
   axis->SetTitleOffset(1.25);
   //   axis->SetTickSize(0.03);
   axis->SetNdivisions(510);
   axis->SetTitle("p_{T} of intersection [GeV]");

   pad1->cd();
 


    cmsPrel(intLumi=36, false);
    latex.DrawLatex(.25,.9,algo+" Jets");
    line_Abseta->Draw();
     leg_eta_comb->Draw();
    c->SaveAs(GetDateDir()+"/ResPtDepComp_impr_"+dir1+"_FSRcorr_residuals_Abseta_"+ algo +".eps");

    central_pt_Abseta->Draw("][hist");
    central_pt_Abseta->GetYaxis()->SetTitle("p_{T} of intersection [GeV]");
    central_pt_Abseta->GetXaxis()->SetTitle("|#eta|");
    central_pt_Abseta->GetYaxis()->SetRangeUser(100,400);
    c->SaveAs(GetDateDir()+"/ResPtDepComp_impr_"+dir1+"_FSRcorr_residuals_Abseta_pt_is_"+ algo +".eps");


     //RATIO-plots eta
     low_pt_eta_hist->Divide(low_pt_eta_hist,import_FSRcorr_residuals_eta_res1_Dir1);
     high_pt_eta_hist->Divide(high_pt_eta_hist,import_FSRcorr_residuals_eta_res1_Dir1);

    TH1D*  pseudo_uncertainty_ratio_eta= get_pseudo_uncertainty_band(high_pt_eta_hist,low_pt_eta_hist);
    pseudo_uncertainty_ratio_eta->SetFillColor(5);
    pseudo_uncertainty_ratio_eta->SetMarkerColor(5);
    pseudo_uncertainty_ratio_eta->SetMarkerStyle(1);

     pseudo_uncertainty_ratio_eta->Draw("e5");
     low_pt_eta_hist->Draw("same");      
     high_pt_eta_hist->Draw("same");      

    leg_eta_comb->DrawClone();

     cmsPrel(intLumi=36, false);
    latex.DrawLatex(.25,.9,algo+" Jets");
    line_eta->Draw();

    pseudo_uncertainty_ratio_eta->GetYaxis()->SetRangeUser(0.9,1.2);
    pseudo_uncertainty_ratio_eta->GetYaxis()->SetTitle("Deviation from nominal residuals");
 
    line_eta->Draw();
  
    cmsPrel(intLumi=36, false);
     c->SaveAs(GetDateDir()+"/ResPtDepComp_impr_"+dir1+"_FSRcorr_residuals_ratio_eta_"+ algo +".eps");

     //RATIO-plots Abseta
     low_pt_Abseta_hist->Divide(low_pt_Abseta_hist,import_FSRcorr_residuals_Abseta_res1_Dir1);
     high_pt_Abseta_hist->Divide(high_pt_Abseta_hist,import_FSRcorr_residuals_Abseta_res1_Dir1);

     TH1D* pseudo_uncertainty_ratio_Abseta= get_pseudo_uncertainty_band(high_pt_Abseta_hist,low_pt_Abseta_hist);
    pseudo_uncertainty_ratio_Abseta->SetFillColor(5);
    pseudo_uncertainty_ratio_Abseta->SetMarkerColor(5);
    pseudo_uncertainty_ratio_Abseta->SetMarkerStyle(1);

     pseudo_uncertainty_ratio_Abseta->Draw("e5");
     low_pt_Abseta_hist->Draw("same");      
     high_pt_Abseta_hist->Draw("same");      

    leg_eta_comb->DrawClone();

     cmsPrel(intLumi=36, false);
    latex.DrawLatex(.25,.9,algo+" Jets");
    line_Abseta->Draw();

    pseudo_uncertainty_ratio_Abseta->GetYaxis()->SetRangeUser(0.9,1.2);
    pseudo_uncertainty_ratio_Abseta->GetYaxis()->SetTitle("Deviation from nominal residuals");
 
    line_Abseta->Draw();
  
    cmsPrel(intLumi=36, false);
     c->SaveAs(GetDateDir()+"/ResPtDepComp_impr_"+dir1+"_FSRcorr_residuals_ratio_Abseta_"+ algo +".eps");

     //     TString dir1, TString algo="PF", TString dir_prefix="KOSTAS_L1_on_pt_plain_", TString label_dir1="", TString binning_select ="kostas", Bool_t pub_style=true){

    TFile *outf = new TFile(GetDateDir()+"/PTDEPENDENCE_"+binning_select+"_"+algo+".root","UPDATE");
    pseudo_uncertainty_eta->SetName(label_dir1+"_pseudo_uncertainty_eta");
    pseudo_uncertainty_ratio_eta->SetName(label_dir1+"_pseudo_uncertainty_ratio_eta");
    pseudo_uncertainty_Abseta->SetName(label_dir1+"_pseudo_uncertainty_Abseta");
    pseudo_uncertainty_ratio_Abseta->SetName(label_dir1+"_pseudo_uncertainty_ratio_Abseta");
    pseudo_uncertainty_eta->Write();
    pseudo_uncertainty_ratio_eta->Write();
    pseudo_uncertainty_Abseta->Write();
    pseudo_uncertainty_ratio_Abseta->Write();
    outf->Close();

}

