#include "tdrstyle_mod.C"
#include "MakeDateDir.h"

void compare_RES_VAL_Data_MC(TString dir1, TString algo="PF", TString dir_prefix="KOSTAS_L1_on_pt_plain_", TString label_dir1="",TString binning_select ="kostas"){
  setTDRStyle();
  TCanvas* c = new TCanvas("c","",600,600);
  TLatex latex;
  latex.SetNDC();
  latex.SetTextSize(0.05);
  latex.SetTextAlign(13);  //align at top
  
    TFile *inf_Dir1;
    inf_Dir1 = new TFile(dir_prefix+dir1+"/res_corrections_histos.root","OPEN");
        if(label_dir1!="")dir1=label_dir1;


    TH1D* import_FSRcorr_residuals_eta_res1_Dir1;
    TH1D* import_FSRcorr_residuals_eta_val1_Dir1;
    TH1D* import_FSRcorr_residuals_Abseta_res1_Dir1;
    TH1D* import_FSRcorr_residuals_Abseta_val1_Dir1;
    import_FSRcorr_residuals_eta_res1_Dir1  = (TH1D*)inf_Dir1->Get(""+binning_select+"_use_easy_mean_TuneZ2_"+algo+"_eta_res_hist");
    import_FSRcorr_residuals_eta_val1_Dir1  = (TH1D*)inf_Dir1->Get(""+binning_select+"_use_easy_mean_TuneZ2_"+algo+"_eta_val_hist");

    import_FSRcorr_residuals_Abseta_res1_Dir1  = (TH1D*)inf_Dir1->Get(""+binning_select+"_use_easy_mean_TuneZ2_"+algo+"_Abseta_res_hist");
    import_FSRcorr_residuals_Abseta_val1_Dir1  = (TH1D*)inf_Dir1->Get(""+binning_select+"_use_easy_mean_TuneZ2_"+algo+"_Abseta_val_hist");


//    import_FSRcorr_residuals_eta_res1_Dir1  = (TH1D*)inf_Dir1->Get(""+binning_select+"_use_coarse_kFSRAbs_use_easy_mean_TuneZ2_"+algo+"_eta_res_hist");
//    import_FSRcorr_residuals_eta_val1_Dir1  = (TH1D*)inf_Dir1->Get(""+binning_select+"_use_coarse_kFSRAbs_use_easy_mean_TuneZ2_"+algo+"_eta_val_hist");
//
//    import_FSRcorr_residuals_Abseta_res1_Dir1  = (TH1D*)inf_Dir1->Get(""+binning_select+"_use_coarse_kFSRAbs_use_easy_mean_TuneZ2_"+algo+"_Abseta_res_hist");
//    import_FSRcorr_residuals_Abseta_val1_Dir1  = (TH1D*)inf_Dir1->Get(""+binning_select+"_use_coarse_kFSRAbs_use_easy_mean_TuneZ2_"+algo+"_Abseta_val_hist");
//   
    Int_t no_binsx = import_FSRcorr_residuals_eta_res1_Dir1->GetNbinsX();
    for (Int_t binx_i=0;binx_i<=no_binsx+1;binx_i++){
      //      cout << "content:" <<import_FSRcorr_residuals_eta_res1_Dir1->GetBinContent(binx_i) << endl;
      if(import_FSRcorr_residuals_eta_res1_Dir1->GetBinContent(binx_i)!=0){import_FSRcorr_residuals_eta_res1_Dir1->SetBinContent(binx_i,1/import_FSRcorr_residuals_eta_res1_Dir1->GetBinContent(binx_i));
	//	cout << "error" << import_FSRcorr_residuals_eta_res1_Dir1->GetBinError(binx_i) << endl;
	import_FSRcorr_residuals_eta_res1_Dir1->SetBinError(binx_i,TMath::Abs(import_FSRcorr_residuals_eta_res1_Dir1->GetBinError(binx_i)/import_FSRcorr_residuals_eta_res1_Dir1->GetBinContent(binx_i)));
	//	cout << "new error" <<import_FSRcorr_residuals_eta_res1_Dir1->GetBinError(binx_i) << endl;
      }
      //      cout << "new bincontent" << import_FSRcorr_residuals_eta_res1_Dir1->GetBinContent(binx_i) << endl;
      if(import_FSRcorr_residuals_eta_val1_Dir1->GetBinContent(binx_i)!=0){import_FSRcorr_residuals_eta_val1_Dir1->SetBinContent(binx_i,1/import_FSRcorr_residuals_eta_val1_Dir1->GetBinContent(binx_i));
	import_FSRcorr_residuals_eta_val1_Dir1->SetBinError(binx_i,TMath::Abs(import_FSRcorr_residuals_eta_val1_Dir1->GetBinError(binx_i)/import_FSRcorr_residuals_eta_val1_Dir1->GetBinContent(binx_i)));
      }
      if(import_FSRcorr_residuals_Abseta_res1_Dir1->GetBinContent(binx_i)!=0){import_FSRcorr_residuals_Abseta_res1_Dir1->SetBinContent(binx_i,1/import_FSRcorr_residuals_Abseta_res1_Dir1->GetBinContent(binx_i));
	import_FSRcorr_residuals_Abseta_res1_Dir1->SetBinError(binx_i,TMath::Abs(import_FSRcorr_residuals_Abseta_res1_Dir1->GetBinError(binx_i)/import_FSRcorr_residuals_Abseta_res1_Dir1->GetBinContent(binx_i)));
      }
      if(import_FSRcorr_residuals_Abseta_val1_Dir1->GetBinContent(binx_i)!=0){import_FSRcorr_residuals_Abseta_val1_Dir1->SetBinContent(binx_i,1/import_FSRcorr_residuals_Abseta_val1_Dir1->GetBinContent(binx_i));
	import_FSRcorr_residuals_Abseta_val1_Dir1->SetBinError(binx_i,TMath::Abs(import_FSRcorr_residuals_Abseta_val1_Dir1->GetBinError(binx_i)/import_FSRcorr_residuals_Abseta_val1_Dir1->GetBinContent(binx_i)));
      }
    }

    import_FSRcorr_residuals_eta_res1_Dir1->GetYaxis()->SetTitle("Data / Monte Carlo");
    import_FSRcorr_residuals_eta_val1_Dir1->GetYaxis()->SetTitle("Data / Monte Carlo");
    import_FSRcorr_residuals_Abseta_res1_Dir1->GetYaxis()->SetTitle("Data / Monte Carlo");
    import_FSRcorr_residuals_Abseta_val1_Dir1->GetYaxis()->SetTitle("Data / Monte Carlo");



    import_FSRcorr_residuals_eta_res1_Dir1->SetStats(0);
    import_FSRcorr_residuals_eta_val1_Dir1->SetStats(0);
    import_FSRcorr_residuals_eta_res1_Dir1->GetXaxis()->SetTitle("#eta");
    import_FSRcorr_residuals_eta_val1_Dir1->GetXaxis()->SetTitle("#eta");

    import_FSRcorr_residuals_Abseta_res1_Dir1->SetStats(0);
    import_FSRcorr_residuals_Abseta_val1_Dir1->SetStats(0);
    import_FSRcorr_residuals_Abseta_res1_Dir1->GetXaxis()->SetTitle("|#eta|");
    import_FSRcorr_residuals_Abseta_val1_Dir1->GetXaxis()->SetTitle("|#eta|");



    //ETA
   TLine *line_eta = new TLine(import_FSRcorr_residuals_eta_res1_Dir1->GetXaxis()->GetXmin(),1.,import_FSRcorr_residuals_eta_res1_Dir1->GetXaxis()->GetXmax(),1.);
    line_eta->SetLineStyle(2);
    line_eta->SetLineColor(1);
 

    import_FSRcorr_residuals_eta_res1_Dir1->SetLineColor(1);
    import_FSRcorr_residuals_eta_res1_Dir1->SetMarkerColor(1);
    import_FSRcorr_residuals_eta_res1_Dir1->SetMarkerStyle(24);
    import_FSRcorr_residuals_eta_res1_Dir1->GetYaxis()->SetRangeUser(0.8,1.2);
    import_FSRcorr_residuals_eta_res1_Dir1->Draw("");

    import_FSRcorr_residuals_eta_val1_Dir1->SetLineColor(2);
    import_FSRcorr_residuals_eta_val1_Dir1->SetMarkerColor(2);
    import_FSRcorr_residuals_eta_val1_Dir1->SetMarkerStyle(21);
    import_FSRcorr_residuals_eta_val1_Dir1->Draw("same");

    

    TLegend *leg_eta_comb;
    leg_eta_comb = new TLegend(0.2,0.75,0.7,0.85);
    leg_eta_comb->SetFillColor(kWhite);
    leg_eta_comb->SetTextFont(42);
    //for Paper
       leg_eta_comb->SetTextSize(0.035);
      //   leg->SetHeader("Legende");
    leg_eta_comb->AddEntry(import_FSRcorr_residuals_eta_res1_Dir1,"Before residual correction","p");
    leg_eta_comb->AddEntry(import_FSRcorr_residuals_eta_val1_Dir1,"After residual correction","p");
    //    leg_eta_comb->AddEntry(import_FSRcorr_residuals_eta_res1_Dir1,"Before residual correction ("+dir1+")","p");
    //    leg_eta_comb->AddEntry(import_FSRcorr_residuals_eta_val1_Dir1,"After residual correction ("+dir1+")","p");
    leg_eta_comb->Draw();
    //"ResVal_Comp_Data_MC_"+
    line_eta->Draw();
    cmsPrel(intLumi=36, false);
    latex.DrawLatex(.25,.9,algo+" Jets");
   
    c->SaveAs(GetDateDir()+"/ResVal_Comp_Data_MC_"+dir1+"_FSRcorr_residuals_eta_"+ algo +".eps");


    import_FSRcorr_residuals_eta_res1_Dir1->Divide(import_FSRcorr_residuals_eta_res1_Dir1,import_FSRcorr_residuals_eta_val1_Dir1);

    import_FSRcorr_residuals_eta_res1_Dir1->Draw();
    import_FSRcorr_residuals_eta_res1_Dir1->GetYaxis()->SetRangeUser(0.95,1.05);
    import_FSRcorr_residuals_eta_res1_Dir1->GetYaxis()->SetTitle("R("+dir1+",L2L3)/R("+dir1+",L2L3res)");
 
    line_eta->Draw();
    cmsPrel(intLumi=36, false);
    latex.DrawLatex(.25,.9,algo+" Jets");
    c->SaveAs(GetDateDir()+"/ResVal_Comp_Data_MC_"+dir1+"_FSRcorr_residuals_ratio_eta_"+ algo +".eps");



    //ABSETA

   TLine *line_Abseta = new TLine(import_FSRcorr_residuals_Abseta_res1_Dir1->GetXaxis()->GetXmin(),1.,import_FSRcorr_residuals_Abseta_res1_Dir1->GetXaxis()->GetXmax(),1.);
    line_Abseta->SetLineStyle(2);
    line_Abseta->SetLineColor(1);

    import_FSRcorr_residuals_Abseta_res1_Dir1->SetLineColor(1);
    import_FSRcorr_residuals_Abseta_res1_Dir1->SetMarkerColor(1);
    import_FSRcorr_residuals_Abseta_res1_Dir1->SetMarkerStyle(24);
    import_FSRcorr_residuals_Abseta_res1_Dir1->GetYaxis()->SetRangeUser(0.8,1.2);
    import_FSRcorr_residuals_Abseta_res1_Dir1->Draw("");

    import_FSRcorr_residuals_Abseta_val1_Dir1->SetLineColor(2);
    import_FSRcorr_residuals_Abseta_val1_Dir1->SetMarkerColor(2);
    import_FSRcorr_residuals_Abseta_val1_Dir1->SetMarkerStyle(21);
    import_FSRcorr_residuals_Abseta_val1_Dir1->Draw("same");

    

    TLegend *leg_Abseta_comb;
    leg_Abseta_comb = new TLegend(0.2,0.75,0.7,0.85);
    leg_Abseta_comb->SetFillColor(kWhite);
    leg_Abseta_comb->SetTextFont(42);
      //   leg->SetHeader("Legende");
    leg_Abseta_comb->AddEntry(import_FSRcorr_residuals_Abseta_res1_Dir1,"Before residual correction","p");
    leg_Abseta_comb->AddEntry(import_FSRcorr_residuals_Abseta_val1_Dir1,"After residual correction","p");
    //    leg_Abseta_comb->AddEntry(import_FSRcorr_residuals_Abseta_res1_Dir1,"Before residual correction ("+dir1+")","p");
    //    leg_Abseta_comb->AddEntry(import_FSRcorr_residuals_Abseta_val1_Dir1,"After residual correction ("+dir1+")","p");
    leg_Abseta_comb->Draw();

    line_Abseta->Draw();
    cmsPrel(intLumi=36, false);
    latex.DrawLatex(.25,.9,algo+" Jets");
    c->SaveAs(GetDateDir()+"/ResVal_Comp_Data_MC_"+dir1+"_FSRcorr_residuals_Abseta_"+ algo +".eps");


    import_FSRcorr_residuals_Abseta_res1_Dir1->Divide(import_FSRcorr_residuals_Abseta_res1_Dir1,import_FSRcorr_residuals_Abseta_val1_Dir1);

    import_FSRcorr_residuals_Abseta_res1_Dir1->Draw();
    import_FSRcorr_residuals_Abseta_res1_Dir1->GetYaxis()->SetRangeUser(0.95,1.05);
    import_FSRcorr_residuals_Abseta_res1_Dir1->GetYaxis()->SetTitle("R("+dir1+",L2L3)/R("+dir1+",L2L3res)");
 
 
    line_Abseta->Draw();
    cmsPrel(intLumi=36, false);
    latex.DrawLatex(.25,.9,algo + " Jets");
    c->SaveAs(GetDateDir()+"/ResVal_Comp_Data_MC_"+dir1+"_FSRcorr_residuals_ratio_Abseta_"+ algo +".eps");





    import_FSRcorr_residuals_eta_res1_Dir1->GetYaxis()->SetTitle("Data / Monte Carlo");
    import_FSRcorr_residuals_eta_val1_Dir1->GetYaxis()->SetTitle("Data / Monte Carlo");
    import_FSRcorr_residuals_Abseta_res1_Dir1->GetYaxis()->SetTitle("Data / Monte Carlo");
    import_FSRcorr_residuals_Abseta_val1_Dir1->GetYaxis()->SetTitle("Data / Monte Carlo");


}

