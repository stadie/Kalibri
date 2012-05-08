#include "tdrstyle_mod.C"
#include "MakeDateDir.h"

TH1D* get_res_hist_at_pt(TH1D* trad_res_hist, TH1D* res_CONST_hist, TH1D* res_SLOPE_hist, Double_t pt){

    TF1 *fit_loglin = new TF1("fit_loglin","[0]+[1]*TMath::Log(x)",trad_res_hist->GetXaxis()->GetXmin(),trad_res_hist->GetXaxis()->GetXmax()); //was used before...
    fit_loglin->SetParameters(1,1);
    fit_loglin->SetParName(0,"const");
    fit_loglin->SetParName(1,"slope");
    Int_t no_x_bins = trad_res_hist->GetNbinsX();

    TH1D* res_hist_at_pt = trad_res_hist->Clone();
    res_hist_at_pt->SetName("res_hist_pt_"+(Long_t)pt);
    res_hist_at_pt->SetTitle("res_hist_pt_"+(Long_t)pt);
    //new TH1D("res_hist_pt_"+(Long_t)pt, "res_hist_pt_"+(Long_t)pt,no_x_bins, trad_res_hist->GetXaxis()->GetXmin() , trad_res_hist->GetXaxis()->GetXmax());

    for(Int_t i=1;i< no_x_bins+1;i++){
      //      cout << "ACHTIUNG" << res_CONST_hist->GetBinContent(i) << endl;
      fit_loglin->SetParameter(0,res_CONST_hist->GetBinContent(i));
      fit_loglin->SetParError(0,res_CONST_hist->GetBinError(i));
      fit_loglin->SetParameter(1,res_SLOPE_hist->GetBinContent(i));
      fit_loglin->SetParError(1,res_SLOPE_hist->GetBinError(i));
      //      fit_loglin->Print();
      res_hist_at_pt->SetBinContent(i,fit_loglin->Eval(pt));
      //      cout << fit_loglin->Eval(pt) << endl;
    }
    return res_hist_at_pt;
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


void compare_PT_dep(TString dir1, TString algo="PF", TString dir_prefix="KOSTAS_L1_on_pt_plain_", TString label_dir1="", TString binning_select ="kostas", Bool_t pub_style=true){
  setTDRStyle();
  TCanvas* c = new TCanvas("c","",600,600);
  TLatex latex;
  latex.SetNDC();
  latex.SetTextSize(0.05);
  latex.SetTextAlign(13);  //align at top

  Double_t abs_lower_pt=50;
  Double_t abs_upper_pt=500;

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
    cout << "kFSR_eq_one_"+binning_select+"_"+algo+"_eta_res_hist" << endl;
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

    TH1D* low_pt_eta_hist = get_res_hist_at_pt(import_FSRcorr_residuals_eta_res1_Dir1, import_FSRcorr_residuals_eta_CONST_FIT_Dir1, import_FSRcorr_residuals_eta_SLOPE_FIT_Dir1, abs_lower_pt);
    TH1D* high_pt_eta_hist = get_res_hist_at_pt(import_FSRcorr_residuals_eta_res1_Dir1, import_FSRcorr_residuals_eta_CONST_FIT_Dir1, import_FSRcorr_residuals_eta_SLOPE_FIT_Dir1, abs_upper_pt);

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
    leg_eta_comb->AddEntry(low_pt_eta_hist,TString("at ")+(Long_t)abs_lower_pt+" GeV","p");
    leg_eta_comb->AddEntry(high_pt_eta_hist,TString("at ")+(Long_t)abs_upper_pt+" GeV","p");
    }
    leg_eta_comb->Draw();
    //"ResComp_"+
    cmsPrel(intLumi=36, false);
    latex.DrawLatex(.25,.9,algo+" Jets");
    line_eta->Draw();

    c->SaveAs(GetDateDir()+"/ResPtDepComp_"+dir1+"_FSRcorr_residuals_eta_"+ algo +".eps");

  //ABSETA



    import_FSRcorr_residuals_Abseta_res1_Dir1->SetLineColor(1);
    import_FSRcorr_residuals_Abseta_res1_Dir1->SetMarkerColor(1);
    import_FSRcorr_residuals_Abseta_res1_Dir1->SetMarkerStyle(24);
    import_FSRcorr_residuals_Abseta_res1_Dir1->GetYaxis()->SetRangeUser(0.85,1.2);
    import_FSRcorr_residuals_Abseta_res1_Dir1->Draw("");

    TH1D* low_pt_Abseta_hist = get_res_hist_at_pt(import_FSRcorr_residuals_Abseta_res1_Dir1, import_FSRcorr_residuals_Abseta_CONST_FIT_Dir1, import_FSRcorr_residuals_Abseta_SLOPE_FIT_Dir1, abs_lower_pt);
    TH1D* high_pt_Abseta_hist = get_res_hist_at_pt(import_FSRcorr_residuals_Abseta_res1_Dir1, import_FSRcorr_residuals_Abseta_CONST_FIT_Dir1, import_FSRcorr_residuals_Abseta_SLOPE_FIT_Dir1, abs_upper_pt);

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




    leg_eta_comb->DrawClone();

     cmsPrel(intLumi=36, false);
    latex.DrawLatex(.25,.9,algo+" Jets");
    line_Abseta->Draw();
     c->SaveAs(GetDateDir()+"/ResPtDepComp_"+dir1+"_FSRcorr_residuals_Abseta_"+ algo +".eps");



     //RATIO-plots eta
     low_pt_eta_hist->Divide(low_pt_eta_hist,import_FSRcorr_residuals_eta_res1_Dir1);
     high_pt_eta_hist->Divide(high_pt_eta_hist,import_FSRcorr_residuals_eta_res1_Dir1);

     pseudo_uncertainty_eta= get_pseudo_uncertainty_band(high_pt_eta_hist,low_pt_eta_hist);
    pseudo_uncertainty_eta->SetFillColor(5);
    pseudo_uncertainty_eta->SetMarkerColor(5);
    pseudo_uncertainty_eta->SetMarkerStyle(1);

     pseudo_uncertainty_eta->Draw("e5");
     low_pt_eta_hist->Draw("same");      
     high_pt_eta_hist->Draw("same");      

    leg_eta_comb->DrawClone();

     cmsPrel(intLumi=36, false);
    latex.DrawLatex(.25,.9,algo+" Jets");
    line_eta->Draw();

    pseudo_uncertainty_eta->GetYaxis()->SetRangeUser(0.9,1.2);
    pseudo_uncertainty_eta->GetYaxis()->SetTitle("Deviation from nominal residuals");
 
    line_eta->Draw();
  
    cmsPrel(intLumi=36, false);
     c->SaveAs(GetDateDir()+"/ResPtDepComp_"+dir1+"_FSRcorr_residuals_ratio_eta_"+ algo +".eps");

     //RATIO-plots Abseta
     low_pt_Abseta_hist->Divide(low_pt_Abseta_hist,import_FSRcorr_residuals_Abseta_res1_Dir1);
     high_pt_Abseta_hist->Divide(high_pt_Abseta_hist,import_FSRcorr_residuals_Abseta_res1_Dir1);

     pseudo_uncertainty_Abseta= get_pseudo_uncertainty_band(high_pt_Abseta_hist,low_pt_Abseta_hist);
    pseudo_uncertainty_Abseta->SetFillColor(5);
    pseudo_uncertainty_Abseta->SetMarkerColor(5);
    pseudo_uncertainty_Abseta->SetMarkerStyle(1);

     pseudo_uncertainty_Abseta->Draw("e5");
     low_pt_Abseta_hist->Draw("same");      
     high_pt_Abseta_hist->Draw("same");      

    leg_eta_comb->DrawClone();

     cmsPrel(intLumi=36, false);
    latex.DrawLatex(.25,.9,algo+" Jets");
    line_Abseta->Draw();

    pseudo_uncertainty_Abseta->GetYaxis()->SetRangeUser(0.9,1.2);
    pseudo_uncertainty_Abseta->GetYaxis()->SetTitle("Deviation from nominal residuals");
 
    line_Abseta->Draw();
  
    cmsPrel(intLumi=36, false);
     c->SaveAs(GetDateDir()+"/ResPtDepComp_"+dir1+"_FSRcorr_residuals_ratio_Abseta_"+ algo +".eps");


}

