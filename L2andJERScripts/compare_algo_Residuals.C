#include <iostream>
#include "TFile.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TLatex.h"
#include "TLegend.h"
#include "MakeDateDir.h"

#include "tdrstyle_mod.C"
#include "plotHistsAndRatio.C"
#include "plotRatios.C"

void do_the_comparison(std::vector <TString> algo_list, TString dir, TString dirlabel, TString dir_prefix="KOSTAS_L1_on_pt_plain_",TString binning_select ="kostas", Bool_t pub_style=true){
  //  do_the_comparison(ALGO_list,dir1,label_dir1,dir_prefix,binning_select,pub_style);

  setTDRStyle();
  Double_t intLumi=36;
  //  TCanvas* c = new TCanvas("c","",600,600);
  TCanvas* c;
  TLatex latex;
  latex.SetNDC();
  latex.SetTextSize(0.08);
  latex.SetTextAlign(13);  //align at top
  //  TList *FileList = new TList();

  TFile* file = new TFile(dir_prefix+dir+"/res_corrections_histos.root","OPEN");

  std::vector < TH1D * > plot_list_;
  std::vector < TH1D * > plot_list_abs_;


  if (file->IsZombie()){
       cout << "Error opening file" << endl;
       //       exit(-1);
    }
    else{
      cout << "ROOT-Datei erfolgreich geladen. " << endl;
      std::cout << dir_prefix+dir+"/res_corrections_histos.root" << std::endl;
    }


  for(unsigned int i=0;i<algo_list.size();i++){
//    plot_list_.push_back((TH1D*)file->Get(""+binning_select+"_TuneZ2_"+algo_list.at(i)+"_eta_res_hist"));
//    plot_list_abs_.push_back((TH1D*)file->Get(""+binning_select+"_TuneZ2_"+algo_list.at(i)+"_Abseta_res_hist"));
//    plot_list_.push_back((TH1D*)file->Get(""+binning_select+"_use_coarse_kFSRAbs_TuneZ2_"+algo_list.at(i)+"_eta_res_hist"));
//    plot_list_abs_.push_back((TH1D*)file->Get(""+binning_select+"_use_coarse_kFSRAbs_TuneZ2_"+algo_list.at(i)+"_Abseta_res_hist"));
    //no fit
    //    plot_list_.push_back((TH1D*)file->Get(""+binning_select+"_use_easy_mean_TuneZ2_"+algo_list.at(i)+"_eta_res_hist"));
    //    plot_list_abs_.push_back((TH1D*)file->Get(""+binning_select+"_use_easy_mean_TuneZ2_"+algo_list.at(i)+"_Abseta_res_hist"));
    //no fit use easy mean
    plot_list_.push_back((TH1D*)file->Get(""+binning_select+"_use_coarse_kFSRAbs_use_easy_mean_TuneZ2_"+algo_list.at(i)+"_eta_res_hist"));
    plot_list_abs_.push_back((TH1D*)file->Get(""+binning_select+"_use_coarse_kFSRAbs_use_easy_mean_TuneZ2_"+algo_list.at(i)+"_Abseta_res_hist"));

    std::cout << "works2" << std::endl;
  }
  

  TH1D* denom_plot = plot_list_.at(0);
  plot_list_.erase (plot_list_.begin());

  c = plotHistsAndRatio(plot_list_,denom_plot,"blab","#eta","Residual Correction",true);

  plot_list_.insert(plot_list_.begin(),denom_plot);
  TString algo_names="";
    TLegend *leg_eta_comb;
    leg_eta_comb = new TLegend(0.2,0.55,0.6,0.80);
    leg_eta_comb->SetFillColor(kWhite);
    leg_eta_comb->SetFillStyle(kNone);
    leg_eta_comb->SetTextFont(42);
    //for Paper
    //    leg_eta_comb->SetTextSize(0.035);
    //         leg_eta_comb->SetHeader(dirlabel);
  for(unsigned int i=0;i<algo_list.size();i++){
    //    leg_eta_comb->AddEntry(denom_plot,"Before residual correction","p");
    leg_eta_comb->AddEntry(plot_list_.at(i),algo_list.at(i)+" Jets","p");
    algo_names=algo_names+algo_list.at(i)+"_";
  }
    leg_eta_comb->Draw();
    latex.DrawLatex(.25,.85,dirlabel);
  c->SaveAs(GetDateDir()+"/ALGO_comp_residuals_w_ratio_"+dirlabel+"_"+ algo_names + ".eps");

  plot_list_.erase (plot_list_.begin());

  c = plotHistsAndRatio(plot_list_,denom_plot,"blab","#eta","Residual Correction",false);

  plot_list_.insert(plot_list_.begin(),denom_plot);
    leg_eta_comb->Draw();
    latex.DrawLatex(.25,.85,dirlabel);
  c->SaveAs(GetDateDir()+"/ALGO_comp_residuals_"+dirlabel+"_"+ algo_names + ".eps");

//  plot_list_.erase (plot_list_.begin());
//
//  c = plotRatios(plot_list_,denom_plot,"blab","#eta","Residual Correction");
//
//  plot_list_.insert(plot_list_.begin(),denom_plot);
//  leg_eta_comb->Draw();
//  latex.DrawLatex(.35,.85,algo+" Jets");
//  c->SaveAs(GetDateDir()+"/SYSTEMATICS_comp_RATIO_"+sample_names+"_"+ algo + ".eps");
//
//
//
//
 TH1D* denom_plot_abs = plot_list_abs_.at(0);
  plot_list_abs_.erase (plot_list_abs_.begin());

  c = plotHistsAndRatio(plot_list_abs_,denom_plot_abs,"blab","|#eta|","Residual Correction",true);
  plot_list_abs_.insert(plot_list_abs_.begin(),denom_plot_abs);

    leg_eta_comb->Draw();
    //    latex.DrawLatex(.35,.85,algo+" Jets");



    latex.DrawLatex(.25,.85,dirlabel);
  c->SaveAs(GetDateDir()+"/ALGO_comp_residuals_w_ratio_Abseta_"+dirlabel+"_"+ algo_names + ".eps");

  plot_list_abs_.erase (plot_list_abs_.begin());

  c = plotHistsAndRatio(plot_list_abs_,denom_plot_abs,"blab","|#eta|","Residual Correction",false);
  plot_list_abs_.insert(plot_list_abs_.begin(),denom_plot_abs);

    leg_eta_comb->Draw();
    //    latex.DrawLatex(.35,.85,algo+" Jets");



    latex.DrawLatex(.25,.85,dirlabel);
  c->SaveAs(GetDateDir()+"/ALGO_comp_residuals_Abseta_"+dirlabel+"_"+ algo_names + ".eps");


//  //  plot_list_.insert(plot_list_.begin(),denom_plot);
//
//  plot_list_abs_.erase (plot_list_abs_.begin());
//
//  c = plotRatios(plot_list_abs_,denom_plot_abs,"blab","|#eta|","Ratio");
//  plot_list_abs_.insert(plot_list_abs_.begin(),denom_plot_abs);
//
//    leg_eta_comb->Draw();
//    latex.DrawLatex(.35,.85,algo+" Jets");
//
//
//
//  c->SaveAs(GetDateDir()+"/SYSTEMATICS_comp_Abseta_RATIO_"+sample_names+"_"+ algo + ".eps");
//  //  plot_list_.insert(plot_list_.begin(),denom_plot);
//
//  for(unsigned int i=0;i<dirlabel_list.size();i++){
//    //    cout << dirlabel_list.at(i) << endl;
//    if(dirlabel_list.at(i)!="")dir_list.at(i)=dirlabel_list.at(i);
//
//  }
}




void compare_algo_Residuals(TString ALGO1, TString ALGO2, TString dir1, TString label_dir1, TString dir_prefix="KOSTAS_L1_on_pt_plain_",TString binning_select ="kostas", Bool_t pub_style=true){
  std::vector <TString> ALGO_list;
  ALGO_list.push_back(ALGO1);
  ALGO_list.push_back(ALGO2);
  do_the_comparison(ALGO_list,dir1,label_dir1,dir_prefix,binning_select,pub_style);


}

void compare_algo_Residuals(TString ALGO1, TString ALGO2, TString ALGO3, TString dir1, TString label_dir1, TString dir_prefix,TString binning_select , Bool_t pub_style){

  std::vector <TString> ALGO_list;
  ALGO_list.push_back(ALGO1);
  ALGO_list.push_back(ALGO2);
  ALGO_list.push_back(ALGO3);
  do_the_comparison(ALGO_list,dir1,label_dir1,dir_prefix,binning_select,pub_style);

}

void compare_algo_Residuals(TString ALGO1, TString ALGO2, TString ALGO3, TString ALGO4, TString dir1, TString label_dir1, TString dir_prefix,TString binning_select, Bool_t pub_style){
  std::vector <TString> ALGO_list;
  ALGO_list.push_back(ALGO1);
  ALGO_list.push_back(ALGO2);
  ALGO_list.push_back(ALGO3);
  ALGO_list.push_back(ALGO4);
  do_the_comparison(ALGO_list,dir1,label_dir1,dir_prefix,binning_select,pub_style);


}

void compare_algo_Residuals(TString ALGO1, TString ALGO2, TString ALGO3, TString ALGO4, TString ALGO5, TString dir1, TString label_dir1, TString dir_prefix,TString binning_select, Bool_t pub_style){
  std::vector <TString> ALGO_list;
  ALGO_list.push_back(ALGO1);
  ALGO_list.push_back(ALGO2);
  ALGO_list.push_back(ALGO3);
  ALGO_list.push_back(ALGO4);
  ALGO_list.push_back(ALGO5);
  do_the_comparison(ALGO_list,dir1,label_dir1,dir_prefix,binning_select,pub_style);
}
