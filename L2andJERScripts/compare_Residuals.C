#include "MakeDateDir.h"
#include "tdrstyle_mod.C"
#include "compare_PT_dep_pt_equiv.C"

void compare_Residuals(TString dir1, TString dir2, TString algo="PF", TString dir_prefix="KOSTAS_L1_on_pt_plain_", TString label_dir1="", TString label_dir2="",TString binning_select ="kostas", Bool_t pub_style=true){
  setTDRStyle();
  TCanvas* c = new TCanvas("c","",600,600);
  TLatex latex;
  latex.SetNDC();
  latex.SetTextSize(0.05);
  latex.SetTextAlign(13);  //align at top
  Bool_t up_down_exist=false;

    TFile *inf_Dir1;
    inf_Dir1 = new TFile(dir_prefix+dir1+"/res_corrections_histos.root","OPEN");

    TFile *inf_Dir2;
    inf_Dir2 = new TFile(dir_prefix+dir2+"/res_corrections_histos.root","OPEN");

    TFile *inf_Dir2_u;
    inf_Dir2_u = new TFile(dir_prefix+dir2+"_u/res_corrections_histos.root","OPEN");
    TFile *inf_Dir2_d;
    inf_Dir2_d = new TFile(dir_prefix+dir2+"_d/res_corrections_histos.root","OPEN");

    if (inf_Dir2_u->IsZombie() || inf_Dir2_d->IsZombie()) {
       cout << "Error opening file" << endl;
       //       exit(-1);
    }
    else{
      cout << "ROOT-Datei erfolgreich geladen. " << endl;
      up_down_exist=true;
    }


        if(label_dir1!="")dir1=label_dir1;
        if(label_dir2!="")dir2=label_dir2;

    TH1D* import_FSRcorr_residuals_eta_res1_Dir1;
    TH1D* import_FSRcorr_residuals_eta_res1_Dir2;
    TH1D* import_FSRcorr_residuals_Abseta_res1_Dir1;
    TH1D* import_FSRcorr_residuals_Abseta_res1_Dir2;


    TH1D* import_FSRcorr_residuals_eta_res1_Dir2_u;
    TH1D* import_FSRcorr_residuals_Abseta_res1_Dir2_u;
    TH1D* import_FSRcorr_residuals_eta_res1_Dir2_d;
    TH1D* import_FSRcorr_residuals_Abseta_res1_Dir2_d;

    import_FSRcorr_residuals_eta_res1_Dir1  = (TH1D*)inf_Dir1->Get(""+binning_select+"_use_coarse_kFSRAbs_TuneZ2_"+algo+"_eta_res_hist");
    import_FSRcorr_residuals_eta_res1_Dir2  = (TH1D*)inf_Dir2->Get(""+binning_select+"_use_coarse_kFSRAbs_TuneZ2_"+algo+"_eta_res_hist");

    import_FSRcorr_residuals_Abseta_res1_Dir1  = (TH1D*)inf_Dir1->Get(""+binning_select+"_use_coarse_kFSRAbs_TuneZ2_"+algo+"_Abseta_res_hist");
    import_FSRcorr_residuals_Abseta_res1_Dir2  = (TH1D*)inf_Dir2->Get(""+binning_select+"_use_coarse_kFSRAbs_TuneZ2_"+algo+"_Abseta_res_hist");

    import_FSRcorr_residuals_eta_res1_Dir2_u  = (TH1D*)inf_Dir2_u->Get(""+binning_select+"_use_coarse_kFSRAbs_TuneZ2_"+algo+"_eta_res_hist");
    import_FSRcorr_residuals_Abseta_res1_Dir2_u  = (TH1D*)inf_Dir2_u->Get(""+binning_select+"_use_coarse_kFSRAbs_TuneZ2_"+algo+"_Abseta_res_hist");
    import_FSRcorr_residuals_eta_res1_Dir2_d  = (TH1D*)inf_Dir2_d->Get(""+binning_select+"_use_coarse_kFSRAbs_TuneZ2_"+algo+"_eta_res_hist");
    import_FSRcorr_residuals_Abseta_res1_Dir2_d  = (TH1D*)inf_Dir2_d->Get(""+binning_select+"_use_coarse_kFSRAbs_TuneZ2_"+algo+"_Abseta_res_hist");

//
// standard until 10 Nov 2011
//    import_FSRcorr_residuals_eta_res1_Dir1  = (TH1D*)inf_Dir1->Get(""+binning_select+"_TuneZ2_"+algo+"_eta_res_hist");
//    import_FSRcorr_residuals_eta_res1_Dir2  = (TH1D*)inf_Dir2->Get(""+binning_select+"_TuneZ2_"+algo+"_eta_res_hist");
//
//    import_FSRcorr_residuals_Abseta_res1_Dir1  = (TH1D*)inf_Dir1->Get(""+binning_select+"_TuneZ2_"+algo+"_Abseta_res_hist");
//    import_FSRcorr_residuals_Abseta_res1_Dir2  = (TH1D*)inf_Dir2->Get(""+binning_select+"_TuneZ2_"+algo+"_Abseta_res_hist");
//
//    import_FSRcorr_residuals_eta_res1_Dir2_u  = (TH1D*)inf_Dir2_u->Get(""+binning_select+"_TuneZ2_"+algo+"_eta_res_hist");
//    import_FSRcorr_residuals_Abseta_res1_Dir2_u  = (TH1D*)inf_Dir2_u->Get(""+binning_select+"_TuneZ2_"+algo+"_Abseta_res_hist");
//    import_FSRcorr_residuals_eta_res1_Dir2_d  = (TH1D*)inf_Dir2_d->Get(""+binning_select+"_TuneZ2_"+algo+"_eta_res_hist");
//    import_FSRcorr_residuals_Abseta_res1_Dir2_d  = (TH1D*)inf_Dir2_d->Get(""+binning_select+"_TuneZ2_"+algo+"_Abseta_res_hist");


//
//    Int_t nbinsx_eta = import_FSRcorr_residuals_eta_res1_Dir2->GetNbinsX();
//    for(Int_t i =0;i< nbinsx_eta;i++){
//      import_FSRcorr_residuals_eta_res1_Dir2->SetBinContent(i,import_FSRcorr_residuals_eta_res1_Dir2->GetBinContent(i)*import_FSRcorr_residuals_eta_res1_Dir1->GetBinContent(i));
//    }
//
//    Int_t nbinsx_Abseta = import_FSRcorr_residuals_Abseta_res1_Dir2->GetNbinsX();
//    for(Int_t i =0;i< nbinsx_Abseta;i++){
//      import_FSRcorr_residuals_Abseta_res1_Dir2->SetBinContent(i,import_FSRcorr_residuals_Abseta_res1_Dir2->GetBinContent(i)*import_FSRcorr_residuals_Abseta_res1_Dir1->GetBinContent(i));
//    }

//
//    TF1 *add1 = new TF1("add1","1",-10,10);
//    import_FSRcorr_residuals_eta_res1_Dir2->Add(add1,1);
//    import_FSRcorr_residuals_Abseta_res1_Dir2->Add(add1,1);
//
//    //void Add(TF1* h1, Double_t c1 = 1, Option_t* option = "")
//    import_FSRcorr_residuals_eta_res1_Dir2->Add( import_FSRcorr_residuals_eta_res1_Dir2,  import_FSRcorr_residuals_eta_res1_Dir1, 1, -1);
//    import_FSRcorr_residuals_Abseta_res1_Dir2->Add( import_FSRcorr_residuals_Abseta_res1_Dir2,  import_FSRcorr_residuals_Abseta_res1_Dir1, 1, -1);
//    //void Add(const TH1* h, const TH1* h2, Double_t c1 = 1, Double_t c2 = 1)
//
//

    import_FSRcorr_residuals_eta_res1_Dir1->SetStats(0);
    import_FSRcorr_residuals_eta_res1_Dir2->SetStats(0);
    import_FSRcorr_residuals_eta_res1_Dir1->GetYaxis()->SetTitle("Residual Correction");
    import_FSRcorr_residuals_eta_res1_Dir2->GetYaxis()->SetTitle("Residual Correction");
    import_FSRcorr_residuals_eta_res1_Dir1->GetXaxis()->SetTitle("#eta");
    import_FSRcorr_residuals_eta_res1_Dir2->GetXaxis()->SetTitle("#eta");

    import_FSRcorr_residuals_Abseta_res1_Dir1->SetStats(0);
    import_FSRcorr_residuals_Abseta_res1_Dir2->SetStats(0);
    import_FSRcorr_residuals_Abseta_res1_Dir1->GetYaxis()->SetTitle("Residual Correction");
    import_FSRcorr_residuals_Abseta_res1_Dir2->GetYaxis()->SetTitle("Residual Correction");
    import_FSRcorr_residuals_Abseta_res1_Dir1->GetXaxis()->SetTitle("|#eta|");
    import_FSRcorr_residuals_Abseta_res1_Dir2->GetXaxis()->SetTitle("|#eta|");



    //ETA

    import_FSRcorr_residuals_eta_res1_Dir1->SetLineColor(1);
    import_FSRcorr_residuals_eta_res1_Dir1->SetMarkerColor(1);
    import_FSRcorr_residuals_eta_res1_Dir1->SetMarkerStyle(24);
    import_FSRcorr_residuals_eta_res1_Dir1->GetYaxis()->SetRangeUser(0.85,1.2);
    import_FSRcorr_residuals_eta_res1_Dir1->Draw("");

    import_FSRcorr_residuals_eta_res1_Dir2->SetLineColor(2);
    import_FSRcorr_residuals_eta_res1_Dir2->SetMarkerColor(2);
    import_FSRcorr_residuals_eta_res1_Dir2->SetMarkerStyle(21);
    import_FSRcorr_residuals_eta_res1_Dir2->Draw("same");

    TH1D* pseudo_uncertainty_eta;
    if(up_down_exist){
      pseudo_uncertainty_eta= get_pseudo_uncertainty_band(import_FSRcorr_residuals_eta_res1_Dir2_u,import_FSRcorr_residuals_eta_res1_Dir2_d);
      pseudo_uncertainty_eta->SetFillColor(5);
      pseudo_uncertainty_eta->SetMarkerColor(5);
      pseudo_uncertainty_eta->SetMarkerStyle(1);
      pseudo_uncertainty_eta->Draw("same e5");
    }
    import_FSRcorr_residuals_eta_res1_Dir1->Draw("same");
    import_FSRcorr_residuals_eta_res1_Dir2->Draw("same");

    TLegend *leg_eta_comb;
    leg_eta_comb = new TLegend(0.25,0.70,0.7,0.85);
    //    leg_eta_comb->SetFillColor(kWhite);
    leg_eta_comb->SetFillStyle(kNone);
    leg_eta_comb->SetTextFont(42);
      //   leg->SetHeader("Legende");
    if(pub_style){
    leg_eta_comb->AddEntry(import_FSRcorr_residuals_eta_res1_Dir1,"2010","p");
    leg_eta_comb->AddEntry(import_FSRcorr_residuals_eta_res1_Dir2,"2011","p");
    }
    else{
      //    leg_eta_comb->AddEntry(import_FSRcorr_residuals_eta_res1_Dir1,"FSR-corrected residuals ("+dir1+")","p");
      //    leg_eta_comb->AddEntry(import_FSRcorr_residuals_eta_res1_Dir2,"FSR-corrected residuals ("+dir2+")","p");
    leg_eta_comb->AddEntry(import_FSRcorr_residuals_eta_res1_Dir1,dir1,"p");
    leg_eta_comb->AddEntry(import_FSRcorr_residuals_eta_res1_Dir2,dir2,"p");
    }
    leg_eta_comb->Draw();
    //"ResComp_"+
    cmsPrel(intLumi=36, false);
    latex.DrawLatex(.25,.9,algo+" Jets");

    c->SaveAs(GetDateDir()+"/ResComp_"+dir1+"_"+dir2+"_FSRcorr_residuals_eta_"+ algo +".eps");

    TH1D* ratio_eta_Dir2_Dir1 = import_FSRcorr_residuals_eta_res1_Dir1->Clone();
    TH1D* ratio_eta_Dir2_Dir1_u = import_FSRcorr_residuals_eta_res1_Dir1->Clone();
    TH1D* ratio_eta_Dir2_Dir1_d = import_FSRcorr_residuals_eta_res1_Dir1->Clone();
    ratio_eta_Dir2_Dir1->Divide(import_FSRcorr_residuals_eta_res1_Dir2,import_FSRcorr_residuals_eta_res1_Dir1);

    if(up_down_exist){
    ratio_eta_Dir2_Dir1_u->Divide(import_FSRcorr_residuals_eta_res1_Dir2_u,import_FSRcorr_residuals_eta_res1_Dir1);
    ratio_eta_Dir2_Dir1_d->Divide(import_FSRcorr_residuals_eta_res1_Dir2_d,import_FSRcorr_residuals_eta_res1_Dir1);
    }
    ratio_eta_Dir2_Dir1->Draw();
    ratio_eta_Dir2_Dir1->GetYaxis()->SetRangeUser(0.95,1.05);
    ratio_eta_Dir2_Dir1->GetYaxis()->SetTitle("R("+dir2+")/R("+dir1+")");


    TH1D* pseudo_uncertainty_ratio_eta;
    if(up_down_exist){
      pseudo_uncertainty_ratio_eta= get_pseudo_uncertainty_band(ratio_eta_Dir2_Dir1_u,ratio_eta_Dir2_Dir1_d);
      pseudo_uncertainty_ratio_eta->SetFillColor(5);
      pseudo_uncertainty_ratio_eta->SetMarkerColor(5);
      pseudo_uncertainty_ratio_eta->SetMarkerStyle(1);
      pseudo_uncertainty_ratio_eta->Draw("same e5");
    }
    ratio_eta_Dir2_Dir1->Draw("same");
  
   TLine *line_eta = new TLine(import_FSRcorr_residuals_eta_res1_Dir1->GetXaxis()->GetXmin(),1.,import_FSRcorr_residuals_eta_res1_Dir1->GetXaxis()->GetXmax(),1.);
    line_eta->Draw();
    line_eta->SetLineStyle(2);
    line_eta->SetLineColor(1);
 
    cmsPrel(intLumi=36, false);
     c->SaveAs(GetDateDir()+"/ResComp_"+dir1+"_"+dir2+"_FSRcorr_residuals_ratio_eta_"+ algo +".eps");

    import_FSRcorr_residuals_eta_res1_Dir2->SetLineColor(1);
    import_FSRcorr_residuals_eta_res1_Dir2->SetMarkerColor(1);
    import_FSRcorr_residuals_eta_res1_Dir2->SetMarkerStyle(24);
    import_FSRcorr_residuals_eta_res1_Dir2->GetYaxis()->SetRangeUser(0.95,1.35);
    import_FSRcorr_residuals_eta_res1_Dir2->Draw("");
    import_FSRcorr_residuals_eta_res1_Dir2->Draw("");
    TLegend *leg_eta;
    leg_eta = new TLegend(0.25,0.75,0.7,0.85);
    //    leg_eta->SetFillColor(kWhite);
    leg_eta->SetFillStyle(kNone);
    leg_eta->SetTextFont(42);
      //   leg->SetHeader("Legende");
    if(pub_style){
    leg_eta->AddEntry(import_FSRcorr_residuals_eta_res1_Dir2,"2010","p");
    }
    else{
      //    leg_eta->AddEntry(import_FSRcorr_residuals_eta_res1_Dir1,"FSR-corrected residuals ("+dir1+")","p");
      //    leg_eta->AddEntry(import_FSRcorr_residuals_eta_res1_Dir2,"FSR-corrected residuals ("+dir2+")","p");
    leg_eta->AddEntry(import_FSRcorr_residuals_eta_res1_Dir2,dir2,"p");
    }
    leg_eta->Draw();
    //"ResComp_"+
    cmsPrel(intLumi=36, false);
    latex.DrawLatex(.25,.9,algo+" Jets");

    c->SaveAs(GetDateDir()+"/ResComp_"+dir2+"_FSRcorr_residuals_eta_"+ algo +".eps");

    //ABSETA


    import_FSRcorr_residuals_Abseta_res1_Dir1->SetLineColor(1);
    import_FSRcorr_residuals_Abseta_res1_Dir1->SetMarkerColor(1);
    import_FSRcorr_residuals_Abseta_res1_Dir1->SetMarkerStyle(24);
    import_FSRcorr_residuals_Abseta_res1_Dir1->GetYaxis()->SetRangeUser(0.85,1.2);
    import_FSRcorr_residuals_Abseta_res1_Dir1->Draw("");

    import_FSRcorr_residuals_Abseta_res1_Dir2->SetLineColor(2);
    import_FSRcorr_residuals_Abseta_res1_Dir2->SetMarkerColor(2);
    import_FSRcorr_residuals_Abseta_res1_Dir2->SetMarkerStyle(21);
    import_FSRcorr_residuals_Abseta_res1_Dir2->Draw("same");


    TH1D* pseudo_uncertainty_Abseta;
    if(up_down_exist){
    pseudo_uncertainty_Abseta= get_pseudo_uncertainty_band(import_FSRcorr_residuals_Abseta_res1_Dir2_u,import_FSRcorr_residuals_Abseta_res1_Dir2_d);
    pseudo_uncertainty_Abseta->SetFillColor(5);
    pseudo_uncertainty_Abseta->SetMarkerColor(5);
    pseudo_uncertainty_Abseta->SetMarkerStyle(1);
    pseudo_uncertainty_Abseta->Draw("same e5");
    }
        
    import_FSRcorr_residuals_Abseta_res1_Dir1->Draw("same");
    import_FSRcorr_residuals_Abseta_res1_Dir2->Draw("same");
   

    TLegend *leg_Abseta_comb;
    leg_Abseta_comb = new TLegend(0.2,0.80,0.7,0.9);
    //    leg_Abseta_comb->SetFillColor(kWhite);
    leg_Abseta_comb->SetFillStyle(kNone);
    leg_Abseta_comb->SetTextFont(42);
      //   leg->SetHeader("Legende");
    if(pub_style){
    leg_Abseta_comb->AddEntry(import_FSRcorr_residuals_Abseta_res1_Dir1,"2010","l");
    leg_Abseta_comb->AddEntry(import_FSRcorr_residuals_Abseta_res1_Dir2,"2011","l");
    }
    else{
      leg_Abseta_comb->AddEntry(import_FSRcorr_residuals_Abseta_res1_Dir1,/*"FSR-corrected residuals ("+*/dir1/*+")"*/,"p");
    leg_Abseta_comb->AddEntry(import_FSRcorr_residuals_Abseta_res1_Dir2,/*"FSR-corrected residuals ("+*/dir2/*+")"*/,"p");
    }
    leg_Abseta_comb->Draw();

    cmsPrel(intLumi=36, false);
     c->SaveAs(GetDateDir()+"/ResComp_"+dir1+"_"+dir2+"_FSRcorr_residuals_Abseta_"+ algo +".eps");



    TH1D* ratio_Abseta_Dir2_Dir1 = import_FSRcorr_residuals_Abseta_res1_Dir1->Clone();
    TH1D* ratio_Abseta_Dir2_Dir1_u = import_FSRcorr_residuals_Abseta_res1_Dir1->Clone();
    TH1D* ratio_Abseta_Dir2_Dir1_d = import_FSRcorr_residuals_Abseta_res1_Dir1->Clone();
    ratio_Abseta_Dir2_Dir1->Divide(import_FSRcorr_residuals_Abseta_res1_Dir2,import_FSRcorr_residuals_Abseta_res1_Dir1);

    if(up_down_exist){
    ratio_Abseta_Dir2_Dir1_u->Divide(import_FSRcorr_residuals_Abseta_res1_Dir2_u,import_FSRcorr_residuals_Abseta_res1_Dir1);
    ratio_Abseta_Dir2_Dir1_d->Divide(import_FSRcorr_residuals_Abseta_res1_Dir2_d,import_FSRcorr_residuals_Abseta_res1_Dir1);
    }
    ratio_Abseta_Dir2_Dir1->Draw();
    ratio_Abseta_Dir2_Dir1->GetYaxis()->SetRangeUser(0.95,1.05);
    ratio_Abseta_Dir2_Dir1->GetYaxis()->SetTitle("R("+dir2+")/R("+dir1+")");


    TH1D* pseudo_uncertainty_ratio_Abseta;
    if(up_down_exist){
    pseudo_uncertainty_ratio_Abseta= get_pseudo_uncertainty_band(ratio_Abseta_Dir2_Dir1_u,ratio_Abseta_Dir2_Dir1_d);
    pseudo_uncertainty_ratio_Abseta->SetFillColor(5);
    pseudo_uncertainty_ratio_Abseta->SetMarkerColor(5);
    pseudo_uncertainty_ratio_Abseta->SetMarkerStyle(1);
    pseudo_uncertainty_ratio_Abseta->Draw("same e5");
    }
   ratio_Abseta_Dir2_Dir1->Draw("same");
 

   TLine *line_Abseta = new TLine(import_FSRcorr_residuals_Abseta_res1_Dir1->GetXaxis()->GetXmin(),1.,import_FSRcorr_residuals_Abseta_res1_Dir1->GetXaxis()->GetXmax(),1.);
    line_Abseta->Draw();
    line_Abseta->SetLineStyle(2);
    line_Abseta->SetLineColor(1);
 
    cmsPrel(intLumi=36, false);
     c->SaveAs(GetDateDir()+"/ResComp_"+dir1+"_"+dir2+"_FSRcorr_residuals_ratio_Abseta_"+ algo +".eps");



    if(up_down_exist){
//      pseudo_uncertainty_eta            ->SetName("pseudo_uncertainty_eta"); 
//      pseudo_uncertainty_ratio_eta      ->SetName("pseudo_uncertainty_ratio_eta"); 
//      pseudo_uncertainty_Abseta         ->SetName("pseudo_uncertainty_Abseta");
//      pseudo_uncertainty_ratio_Abseta   ->SetName("pseudo_uncertainty_ratio_Abseta");

    TFile *outf = new TFile("Additional_up_down_uncertaintybands_"+binning_select+"_"+algo+".root","UPDATE");
    pseudo_uncertainty_eta->SetName(label_dir2+"_pseudo_uncertainty_eta");
    pseudo_uncertainty_ratio_eta->SetName(label_dir2+"_"+label_dir1+"_pseudo_uncertainty_ratio_eta");
    pseudo_uncertainty_Abseta->SetName(label_dir2+"_pseudo_uncertainty_Abseta");
    pseudo_uncertainty_ratio_Abseta->SetName(label_dir2+"_"+label_dir1+"_pseudo_uncertainty_ratio_Abseta");
    pseudo_uncertainty_eta->Write();
    pseudo_uncertainty_ratio_eta->Write();
    pseudo_uncertainty_Abseta->Write();
    pseudo_uncertainty_ratio_Abseta->Write();



    import_FSRcorr_residuals_eta_res1_Dir2_u    ->SetName(label_dir2+"_import_FSRcorr_residuals_eta_res1_Dir2_u");        
    import_FSRcorr_residuals_Abseta_res1_Dir2_u ->SetName(label_dir2+"_import_FSRcorr_residuals_Abseta_res1_Dir2_u");           
    import_FSRcorr_residuals_eta_res1_Dir2_d    ->SetName(label_dir2+"_import_FSRcorr_residuals_eta_res1_Dir2_d");        
    import_FSRcorr_residuals_Abseta_res1_Dir2_d ->SetName(label_dir2+"_import_FSRcorr_residuals_Abseta_res1_Dir2_d");           

    import_FSRcorr_residuals_eta_res1_Dir2_u    ->Write();
    import_FSRcorr_residuals_Abseta_res1_Dir2_u ->Write();
    import_FSRcorr_residuals_eta_res1_Dir2_d    ->Write();
    import_FSRcorr_residuals_Abseta_res1_Dir2_d ->Write();





    outf->Close();

    }


}

