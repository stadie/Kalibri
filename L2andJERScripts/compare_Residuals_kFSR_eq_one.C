#include "tdrstyle_mod.C"
#include "MakeDateDir.h"

void compare_Residuals_kFSR_eq_one(TString dir1, TString dir2, TString
  algo="PF", TString dir_prefix="KOSTAS_L1_on_pt_plain_", TString
  label_dir1="", TString label_dir2="" ,TString binning_select ="kostas"){
  setTDRStyle();
  TCanvas* c = new TCanvas("c","",600,600);
  TLatex latex;
  latex.SetNDC();
  latex.SetTextSize(0.05);
  latex.SetTextAlign(13);  //align at top

    TFile *inf_Dir1;
    inf_Dir1 = new TFile(dir_prefix+dir1+"/res_corrections_histos.root","OPEN");

    TFile *inf_Dir2;
    inf_Dir2 = new TFile(dir_prefix+dir2+"/res_corrections_histos.root","OPEN");


    TH1D* import_FSRcorr_residuals_eta_res1_Dir1;
    TH1D* import_FSRcorr_residuals_eta_res1_Dir2;
    TH1D* import_FSRcorr_residuals_Abseta_res1_Dir1;
    TH1D* import_FSRcorr_residuals_Abseta_res1_Dir2;

    import_FSRcorr_residuals_eta_res1_Dir1  = (TH1D*)inf_Dir1->Get("kFSR_eq_one_"+binning_select+"_use_easy_mean_TuneZ2_"+algo+"_eta_res_hist");
    import_FSRcorr_residuals_eta_res1_Dir2  = (TH1D*)inf_Dir2->Get("kFSR_eq_one_"+binning_select+"_use_easy_mean_TuneZ2_"+algo+"_eta_res_hist");

    import_FSRcorr_residuals_Abseta_res1_Dir1  = (TH1D*)inf_Dir1->Get("kFSR_eq_one_"+binning_select+"_use_easy_mean_TuneZ2_"+algo+"_Abseta_res_hist");
    import_FSRcorr_residuals_Abseta_res1_Dir2  = (TH1D*)inf_Dir2->Get("kFSR_eq_one_"+binning_select+"_use_easy_mean_TuneZ2_"+algo+"_Abseta_res_hist");

    //    #LT r_{Sample1}/r_{Sample2}  #GT_{#alpha = 0.2}

        if(label_dir1!="")dir1=label_dir1;
        if(label_dir2!="")dir2=label_dir2;


    import_FSRcorr_residuals_eta_res1_Dir1->SetStats(0);
    import_FSRcorr_residuals_eta_res1_Dir2->SetStats(0);
    import_FSRcorr_residuals_eta_res1_Dir1->GetYaxis()->SetTitle(" #LT r_{MC}/r_{data}  #GT_{#alpha #leq 0.2}");
    import_FSRcorr_residuals_eta_res1_Dir2->GetYaxis()->SetTitle(" #LT r_{MC}/r_{data}  #GT_{#alpha #leq 0.2}");
    import_FSRcorr_residuals_eta_res1_Dir1->GetXaxis()->SetTitle("#eta");
    import_FSRcorr_residuals_eta_res1_Dir2->GetXaxis()->SetTitle("#eta");

    import_FSRcorr_residuals_Abseta_res1_Dir1->SetStats(0);
    import_FSRcorr_residuals_Abseta_res1_Dir2->SetStats(0);
    import_FSRcorr_residuals_Abseta_res1_Dir1->GetYaxis()->SetTitle(" #LT r_{MC}/r_{data}  #GT_{#alpha #leq 0.2}");
    import_FSRcorr_residuals_Abseta_res1_Dir2->GetYaxis()->SetTitle(" #LT r_{MC}/r_{data}  #GT_{#alpha #leq 0.2}");
    import_FSRcorr_residuals_Abseta_res1_Dir1->GetXaxis()->SetTitle("|#eta|");
    import_FSRcorr_residuals_Abseta_res1_Dir2->GetXaxis()->SetTitle("|#eta|");



    //ETA

    import_FSRcorr_residuals_eta_res1_Dir1->SetLineColor(1);
    import_FSRcorr_residuals_eta_res1_Dir1->SetMarkerColor(1);
    import_FSRcorr_residuals_eta_res1_Dir1->SetMarkerStyle(24);
    import_FSRcorr_residuals_eta_res1_Dir1->GetYaxis()->SetRangeUser(0.95,1.25);
    import_FSRcorr_residuals_eta_res1_Dir1->Draw("");

    import_FSRcorr_residuals_eta_res1_Dir2->SetLineColor(2);
    import_FSRcorr_residuals_eta_res1_Dir2->SetMarkerColor(2);
    import_FSRcorr_residuals_eta_res1_Dir2->SetMarkerStyle(21);
    import_FSRcorr_residuals_eta_res1_Dir2->Draw("same");

    

    TLegend *leg_eta_comb;
    leg_eta_comb = new TLegend(0.25,0.75,0.7,0.85);
    leg_eta_comb->SetFillColor(kWhite);
    leg_eta_comb->SetTextFont(42);
      //   leg->SetHeader("Legende");
    leg_eta_comb->AddEntry(import_FSRcorr_residuals_eta_res1_Dir1," "+dir1+"","p");
    leg_eta_comb->AddEntry(import_FSRcorr_residuals_eta_res1_Dir2," "+dir2+"","p");
    leg_eta_comb->Draw();
    //"ResComp_kFSR_eq_one_"+
    cmsPrel(intLumi=36, false);
     latex.DrawLatex(.25,.9,algo+" Jets");

    c->SaveAs(GetDateDir()+"/ResComp_kFSR_eq_one_"+dir1+"_"+dir2+"_FSRcorr_residuals_eta_"+ algo +".eps");


    import_FSRcorr_residuals_eta_res1_Dir1->Divide(import_FSRcorr_residuals_eta_res1_Dir1,import_FSRcorr_residuals_eta_res1_Dir2);

    import_FSRcorr_residuals_eta_res1_Dir1->Draw();
    import_FSRcorr_residuals_eta_res1_Dir1->GetYaxis()->SetRangeUser(0.95,1.05);
    import_FSRcorr_residuals_eta_res1_Dir1->GetYaxis()->SetTitle("R("+dir1+")/R("+dir2+")");
 
   TLine *line_eta = new TLine(import_FSRcorr_residuals_eta_res1_Dir1->GetXaxis()->GetXmin(),1.,import_FSRcorr_residuals_eta_res1_Dir1->GetXaxis()->GetXmax(),1.);
    line_eta->Draw();
    line_eta->SetLineStyle(2);
    line_eta->SetLineColor(1);
 
    cmsPrel(intLumi=36, false);
     latex.DrawLatex(.25,.9,algo+" Jets");
    c->SaveAs(GetDateDir()+"/ResComp_kFSR_eq_one_"+dir1+"_"+dir2+"_FSRcorr_residuals_ratio_eta_"+ algo +".eps");



    //ABSETA


    import_FSRcorr_residuals_Abseta_res1_Dir1->SetLineColor(1);
    import_FSRcorr_residuals_Abseta_res1_Dir1->SetMarkerColor(1);
    import_FSRcorr_residuals_Abseta_res1_Dir1->SetMarkerStyle(24);
    import_FSRcorr_residuals_Abseta_res1_Dir1->GetYaxis()->SetRangeUser(0.95,1.25);
    import_FSRcorr_residuals_Abseta_res1_Dir1->Draw("");

    import_FSRcorr_residuals_Abseta_res1_Dir2->SetLineColor(2);
    import_FSRcorr_residuals_Abseta_res1_Dir2->SetMarkerColor(2);
    import_FSRcorr_residuals_Abseta_res1_Dir2->SetMarkerStyle(21);
    import_FSRcorr_residuals_Abseta_res1_Dir2->Draw("same");

    

    TLegend *leg_Abseta_comb;
    leg_Abseta_comb = new TLegend(0.25,0.75,0.7,0.85);
    leg_Abseta_comb->SetFillColor(kWhite);
    leg_Abseta_comb->SetTextFont(42);
      //   leg->SetHeader("Legende");
    leg_Abseta_comb->AddEntry(import_FSRcorr_residuals_Abseta_res1_Dir1," "+dir1+"","p");
    leg_Abseta_comb->AddEntry(import_FSRcorr_residuals_Abseta_res1_Dir2," "+dir2+"","p");
    leg_Abseta_comb->Draw();

    cmsPrel(intLumi=36, false);
     latex.DrawLatex(.25,.9,algo+" Jets");
    c->SaveAs(GetDateDir()+"/ResComp_kFSR_eq_one_"+dir1+"_"+dir2+"_FSRcorr_residuals_Abseta_"+ algo +".eps");


    import_FSRcorr_residuals_Abseta_res1_Dir1->Divide(import_FSRcorr_residuals_Abseta_res1_Dir1,import_FSRcorr_residuals_Abseta_res1_Dir2);

    import_FSRcorr_residuals_Abseta_res1_Dir1->Draw();
    import_FSRcorr_residuals_Abseta_res1_Dir1->GetYaxis()->SetRangeUser(0.95,1.05);
    import_FSRcorr_residuals_Abseta_res1_Dir1->GetYaxis()->SetTitle("R("+dir1+")/R("+dir2+")");
 
   TLine *line_Abseta = new TLine(import_FSRcorr_residuals_Abseta_res1_Dir1->GetXaxis()->GetXmin(),1.,import_FSRcorr_residuals_Abseta_res1_Dir1->GetXaxis()->GetXmax(),1.);
    line_Abseta->Draw();
    line_Abseta->SetLineStyle(2);
    line_Abseta->SetLineColor(1);
 
    cmsPrel(intLumi=36, false);
     latex.DrawLatex(.25,.9,algo+" Jets");
    c->SaveAs(GetDateDir()+"/ResComp_kFSR_eq_one_"+dir1+"_"+dir2+"_FSRcorr_residuals_ratio_Abseta_"+ algo +".eps");




}

