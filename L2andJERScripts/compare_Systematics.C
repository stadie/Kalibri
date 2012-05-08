#include <iostream>
#include "MakeDateDir.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TF1.h"
#include "TLatex.h"
#include "TLegend.h"
#include "THelpers_2.h"

#include "tdrstyle_mod.C"
#include "plotHistsAndRatio.C"
#include "plotRatios.C"
#include "plotRatios_wFit.C"
#include "plotRatios_wAbsUncFit.C"
#include "plotHistswithFittedSystematics.C"
#include "plotHistswSomeUncAndAbsUncertainty.C"
void do_the_comparison(std::vector <TString> dir_list, std::vector <TString> dirlabel_list, TString algo="PF", TString dir_prefix="KOSTAS_L1_on_pt_plain_",TString binning_select ="kostas", Bool_t pub_style=true){
  setTDRStyle();
  Double_t intLumi=36;
  //  TCanvas* c = new TCanvas("c","",600,600);
  TCanvas* c;
  TLatex latex;
  latex.SetNDC();
  latex.SetTextSize(0.05);
  latex.SetTextAlign(13);  //align at top
  //  TList *FileList = new TList();
  std::vector <TFile*> file_list;
  for(unsigned int i=0;i<dir_list.size();i++){
    std::cout << dir_prefix <<dir_list.at(i)<<"/res_corrections_histos.root" << std::endl;
    //    FileList->Add( TFile::Open(dir_prefix+dir_list.at(i)+"/res_corrections_histos.root") );

        file_list.push_back(new TFile(dir_prefix+dir_list.at(i)+"/res_corrections_histos.root","OPEN"));
	std::cout << file_list.at(i) << std::endl;
	//        file_list.at(i)->ls();
        if (file_list.at(i)->IsZombie()) {
	  std::cout << "Error opening file" << std::endl;
           exit(-1);
        }
  }
//
//  TFile *first_source = (TFile*)FileList->First();
//

  std::vector < TH1D * > plot_list_;
  std::vector < TH1D * > plot_list_abs_;


  for(unsigned int i=0;i<dir_list.size();i++){
    //    cout << "works1 " << file_list.size() << " " <<i << endl;
    //    std::cout << ""+binning_select+"_use_coarse_kFSRAbs_TuneZ2_"+algo+"_eta_res_hist" << std::endl;
    //first_source->ls();
    //    file_list.at(i)->ls();
    //use fitted    
    //plot_list_.push_back((TH1D*)file_list.at(i)->Get(""+binning_select+"_use_coarse_kFSRAbs_TuneZ2_"+algo+"_eta_res_hist"));
    //plot_list_abs_.push_back((TH1D*)file_list.at(i)->Get(""+binning_select+"_use_coarse_kFSRAbs_TuneZ2_"+algo+"_Abseta_res_hist"));
    //no fit
        plot_list_.push_back((TH1D*)file_list.at(i)->Get(""+binning_select+"_use_easy_mean_TuneZ2_"+algo+"_eta_res_hist"));
        plot_list_abs_.push_back((TH1D*)file_list.at(i)->Get(""+binning_select+"_use_easy_mean_TuneZ2_"+algo+"_Abseta_res_hist"));
    plot_list_.back()->SetName(dirlabel_list.at(i)+"_"+algo+"_eta");
    plot_list_abs_.back()->SetName(dirlabel_list.at(i)+"_"+algo+"_Abseta");
    //    cout << dir_prefix <<dir_list.at(i)<<"/res_corrections_histos.root" << endl;
    //    file_list.push_back(new TFile(dir_prefix+dir_list.at(i)+"/res_corrections_histos.root","OPEN"));
    //    plot_list_.at(i)->Draw();
    //    c2->cd();
    std::cout << "works2" << std::endl;
  }
  
  TFile *inf = new TFile(GetDateDir()+"/PTDEPENDENCE_"+binning_select+"_"+algo+".root","OPEN");
  TH1D* ratio_abseta_pt_dependence;
  inf->GetObject(dirlabel_list.at(0)+"_pseudo_uncertainty_ratio_Abseta",ratio_abseta_pt_dependence);
  //      fone->GetObject((easy_mean_prefix+XVsPt+ptthreecuts[cut_i])+("/"+easy_mean_prefix+XVsPt+ptthreecuts[cut_i])+"_"+XVsPtType+"_MC_L2L3_Eta"+(Long_t)eta_i+ratio_of_mean_or_GM,hMC1);
  if(!ratio_abseta_pt_dependence)cout << "FEHLER - NICHT IMPORTIERT!" << endl;
  for(int bin_i = 1; bin_i <= ratio_abseta_pt_dependence->GetNbinsX(); bin_i++)ratio_abseta_pt_dependence->SetBinContent(bin_i,ratio_abseta_pt_dependence->GetBinContent(bin_i)-1);
  TH1D* abseta_pt_dependence;
  inf->GetObject(dirlabel_list.at(0)+"_pseudo_uncertainty_Abseta",abseta_pt_dependence);
  if(!abseta_pt_dependence)cout << "FEHLER - NICHT IMPORTIERT!" << endl;


  TH1D* denom_plot = plot_list_.at(0);
  plot_list_.erase (plot_list_.begin());

  c = plotHistsAndRatio(plot_list_,denom_plot,"blab","#eta","Residual Correction");

  plot_list_.insert(plot_list_.begin(),denom_plot);
  TString sample_names="";

    TLegend *leg_eta_comb;
    leg_eta_comb = new TLegend(0.15,0.55,0.75,0.80);
    leg_eta_comb->SetFillColor(kWhite);
    leg_eta_comb->SetFillStyle(kNone);
    leg_eta_comb->SetTextFont(42);
    //for Paper
    //       leg_eta_comb->SetTextSize(0.035);
      //   leg->SetHeader("Legende");
  for(unsigned int i=0;i<dirlabel_list.size();i++){
    //    leg_eta_comb->AddEntry(denom_plot,"Before residual correction","p");
    leg_eta_comb->AddEntry(plot_list_.at(i),dirlabel_list.at(i),"p");
    sample_names=sample_names+dirlabel_list.at(i)+"_";
  }
      //    leg_eta_comb->AddEntry(import_FSRcorr_residuals_eta_res1_Dir1,"Before residual correction ("+dir1+")","p");
    //    leg_eta_comb->AddEntry(import_FSRcorr_residuals_eta_val1_Dir1,"After residual correction ("+dir1+")","p");
    leg_eta_comb->Draw();
    //"ResVal_Comp_Data_MC_"+
    //    line_eta->Draw();

    //    cmsPrel(intLumi=36, false);
    latex.DrawLatex(.20,.85,algo+" Jets");


    TLegend *leg_eta;
    leg_eta = new TLegend(0.15,0.55,0.75,0.80);
    leg_eta->SetFillColor(kWhite);
    leg_eta->SetFillStyle(kNone);
    leg_eta->SetTextFont(42);
    //for Paper
    //       leg_eta->SetTextSize(0.035);
    //   leg->SetHeader("Legende");
    for(unsigned int i=1;i<dirlabel_list.size();i++){
      leg_eta->AddEntry(plot_list_.at(i),dirlabel_list.at(i)+"/"+dirlabel_list.at(0),"p");
    }





    Double_t xmin = denom_plot->GetXaxis()->GetXmin(); 
    Double_t xmax = denom_plot->GetXaxis()->GetXmax(); 
    Double_t y=1.0;
    Double_t ye=0.02;


  c->SaveAs(GetDateDir()+"/SYSTEMATICS_comp_"+sample_names+"_"+ algo + ".eps");
  plot_list_.erase (plot_list_.begin());

  c = plotRatios(plot_list_,denom_plot,"blab","#eta","Residual Correction");

  plot_list_.insert(plot_list_.begin(),denom_plot);
  latex.DrawLatex(.20,.85,algo+" Jets");
  leg_eta->Draw();
  c->SaveAs(GetDateDir()+"/SYSTEMATICS_comp_RATIO_"+sample_names+"_"+ algo + ".eps");




  TH1D* denom_plot_abs = plot_list_abs_.at(0);
  plot_list_abs_.erase (plot_list_abs_.begin());

  c = plotHistsAndRatio(plot_list_abs_,denom_plot_abs,"blab","|#eta|","Residual Correction");
  plot_list_abs_.insert(plot_list_abs_.begin(),denom_plot_abs);

    leg_eta_comb->Draw();
    latex.DrawLatex(.20,.85,algo+" Jets");

    Double_t xmin_abs = denom_plot_abs->GetXaxis()->GetXmin(); 
    Double_t xmax_abs = denom_plot_abs->GetXaxis()->GetXmax(); 

  c->SaveAs(GetDateDir()+"/SYSTEMATICS_comp_Abseta_"+sample_names+"_"+ algo + ".eps");
  //  plot_list_.insert(plot_list_.begin(),denom_plot);

  plot_list_abs_.erase (plot_list_abs_.begin());

  c = plotRatios(plot_list_abs_,denom_plot_abs,"blab","|#eta|","Ratio");
  plot_list_abs_.insert(plot_list_abs_.begin(),denom_plot_abs);

    leg_eta->Draw();
    latex.DrawLatex(.20,.85,algo+" Jets");



  c->SaveAs(GetDateDir()+"/SYSTEMATICS_comp_Abseta_RATIO_"+sample_names+"_"+ algo + ".eps");
  //  plot_list_.insert(plot_list_.begin(),denom_plot);


  std::vector<int> fillStyles; 
  //  fillStyles.push_back(3002);
  //  fillStyles.push_back(3004);
  //  fillStyles.push_back(3305);
  //  fillStyles.push_back(3490);
  fillStyles.push_back(3305);
  fillStyles.push_back(3350);
  fillStyles.push_back(3004);
  fillStyles.push_back(3005);
  fillStyles.push_back(3006);

    leg_eta = new TLegend(0.15,0.55,0.75,0.80);
    leg_eta->SetFillColor(kWhite);
    leg_eta->SetFillStyle(kNone);
    leg_eta->SetTextFont(42);
    //for Paper
    //       leg_eta->SetTextSize(0.035);
    //   leg->SetHeader("Legende");
    for(unsigned int i=1;i<dirlabel_list.size();i++){
      if(dirlabel_list.at(i).Contains("Time")){
	cout << i << " i index" << endl;
	plot_list_.at(i)->SetFillStyle(fillStyles.at(i-1));
	plot_list_.at(i)->SetFillColor(plot_list_.at(i)->GetLineColor());
	leg_eta->AddEntry(plot_list_.at(i),dirlabel_list.at(i)+"/"+dirlabel_list.at(0),"f");
      }
      else leg_eta->AddEntry(plot_list_.at(i),dirlabel_list.at(i)+"/"+dirlabel_list.at(0),"l");
    }

  plot_list_abs_.erase (plot_list_abs_.begin());
  ratio_abseta_pt_dependence->SetLineColor(41);
  ratio_abseta_pt_dependence->SetFillColor(41);
  ratio_abseta_pt_dependence->SetFillStyle(3359);
  ratio_abseta_pt_dependence->SetMarkerStyle(1);
  c = plotRatios_wFit(plot_list_abs_,denom_plot_abs,"blab","|#eta|","Ratio",ratio_abseta_pt_dependence);
  plot_list_abs_.insert(plot_list_abs_.begin(),denom_plot_abs);
  //  ratio_abseta_pt_dependence->Draw("same e3");
  denom_plot_abs->SetFillColor(kGray);
  denom_plot_abs->SetFillStyle(3001);
  leg_eta->AddEntry(denom_plot_abs,"absolute uncertainty (w.o. pt)","f");
  
  TLegend* leg_eta_wo_pt = (TLegend*) leg_eta->Clone();
  leg_eta->AddEntry(ratio_abseta_pt_dependence,"pt-dependence","f");
  
    leg_eta->Draw();
    latex.DrawLatex(.20,.85,algo+" Jets");



  c->SaveAs(GetDateDir()+"/SYSTEMATICS_comp_Abseta_RATIO_wFit_"+sample_names+"_"+ algo + ".eps");
  //  plot_list_.insert(plot_list_.begin(),denom_plot);





  plot_list_abs_.erase (plot_list_abs_.begin());
  c = plotRatios_wAbsUncFit(plot_list_abs_,denom_plot_abs,"blab","|#eta|","Ratio");
  plot_list_abs_.insert(plot_list_abs_.begin(),denom_plot_abs);
  //  ratio_abseta_pt_dependence->Draw("same e3");
  denom_plot_abs->SetFillColor(kGray);
  denom_plot_abs->SetFillStyle(3001);
  
    leg_eta_wo_pt->Draw();
    latex.DrawLatex(.20,.85,algo+" Jets");



  c->SaveAs(GetDateDir()+"/SYSTEMATICS_comp_Abseta_RATIO_wABSUNCFit_"+sample_names+"_"+ algo + ".eps");
  //  plot_list_.insert(plot_list_.begin(),denom_plot);









  //    plot_list_abs_.back()->SetName(dirlabel_list.at(i)+"_"+algo+"_Abseta");

  TFile *read_systematics = new TFile(GetDateDir()+"/"+((TString)"Systematic_fit_funcs_rel_to_"+denom_plot_abs->GetName())+".root","OPEN");

  std::vector <TF1*> fit_func_list_abs_;
  cout << dir_list.size() << endl;
  for(unsigned int i=1;i<dir_list.size();i++){
    cout << (TString)"Fit_func_"+plot_list_abs_.at(i)->GetName() << endl;
    fit_func_list_abs_.push_back((TF1*)read_systematics->Get((TString)"Fit_func_"+plot_list_abs_.at(i)->GetName()));
    if(!fit_func_list_abs_.back())cout << "FEHLER - NICHT IMPORTIERT!" << endl;
    else cout << "it worked" << endl;
  }

  c = plotHistswithFittedSystematics(fit_func_list_abs_,denom_plot_abs,"blab","|#eta|","Residual correction",ratio_abseta_pt_dependence);
    leg_eta->Draw();
    latex.DrawLatex(.20,.85,algo+" Jets");

  c->SaveAs(GetDateDir()+"/SYSTEMATICS_comp_Abseta_fitted_uncertainties_"+sample_names+"_"+ algo + ".eps");
  //  plot_list_.insert(plot_list_.begin(),denom_plot);






  TFile *read_HIST_systematics = new TFile(GetDateDir()+"/"+((TString)"Systematic_fit_funcs_HIST_ABSUNC_rel_to_"+denom_plot_abs->GetName())+".root","OPEN");

  TH1D* abs_uncertainty_not_fitted;
  read_HIST_systematics->GetObject("abs_uncertainty_band",abs_uncertainty_not_fitted);

  std::vector <TF1*> some_fit_func_list_abs_;
  std::vector <TH1D*> some_ratio_unc_list_abs_;
  std::vector <TH1D*> dummy_hists_for_legend;
  cout << dir_list.size() << endl;
  for(unsigned int i=1;i<dir_list.size();i++){
    Bool_t import = false;
    if(dirlabel_list.at(i).Contains("JER-syst"))import=true;
    if(dirlabel_list.at(i).Contains("Radiation"))import=true;
    if(import){
      cout << (TString)"Fit_func_"+plot_list_abs_.at(i)->GetName() << endl;
      dummy_hists_for_legend.push_back((TH1D*)plot_list_abs_.at(i)->Clone());
      dummy_hists_for_legend.back()->SetName(dirlabel_list.at(i));
      some_fit_func_list_abs_.push_back((TF1*)read_systematics->Get((TString)"Fit_func_"+plot_list_abs_.at(i)->GetName()));
      some_ratio_unc_list_abs_.push_back((TH1D*)read_systematics->Get((TString)plot_list_abs_.at(i)->GetName()));
      if(!some_fit_func_list_abs_.back())cout << "FEHLER - NICHT IMPORTIERT!" << endl;
      else cout << "it worked" << endl;
    }
  }

  //  c = plotHistswSomeUncAndAbsUncertainty(some_fit_func_list_abs_,denom_plot_abs,"blab","|#eta|","Residual correction",abs_uncertainty_not_fitted);
  c = plotHistswSomeUncAndAbsUncertainty(some_ratio_unc_list_abs_,denom_plot_abs,"blab","|#eta|","Residual correction",abs_uncertainty_not_fitted);

    leg_eta = new TLegend(0.15,0.55,0.75,0.80);
    leg_eta->SetFillColor(kWhite);
    leg_eta->SetFillStyle(kNone);
    leg_eta->SetTextFont(42);
    //for Paper
    //       leg_eta->SetTextSize(0.035);
    //   leg->SetHeader("Legende");
    for(unsigned int i=0;i<dummy_hists_for_legend.size();i++){
      //      if(dirlabel_list.at(match_for_legend.at(i)).Contains("Radiation")){
	dummy_hists_for_legend.at(i)->SetFillStyle(fillStyles.at(i));
	dummy_hists_for_legend.at(i)->SetLineColor(plot_list_abs_.at(i+1)->GetLineColor());
	dummy_hists_for_legend.at(i)->SetFillColor(plot_list_abs_.at(i+1)->GetLineColor());
	leg_eta->AddEntry(dummy_hists_for_legend.at(i),dummy_hists_for_legend.at(i)->GetName(),"f");
	//      }
      //      else leg_eta->AddEntry(dummy_hists_for_legend.at(i),dirlabel_list.at(i)+"/"+dirlabel_list.at(0),"l");
    }
    leg_eta->AddEntry(denom_plot_abs,"absolute uncertainty","f");
    leg_eta->Draw();
    latex.DrawLatex(.20,.85,algo+" Jets");



  c->SaveAs(GetDateDir()+"/SYSTEMATICS_comp_Abseta_SomeUncAndAbsUncertainty_"+sample_names+"_"+ algo + ".eps");
  //  plot_list_.insert(plot_list_.begin(),denom_plot);


  for(unsigned int i=0;i<dirlabel_list.size();i++){
    //    cout << dirlabel_list.at(i) << endl;
    if(dirlabel_list.at(i)!="")dir_list.at(i)=dirlabel_list.at(i);

  }
}




void compare_Systematics(TString dir1, TString dir2, TString label_dir1="", TString label_dir2="", TString algo="PF", TString dir_prefix="KOSTAS_L1_on_pt_plain_",TString binning_select ="kostas", Bool_t pub_style=true){
  std::vector <TString> dir_list;
  dir_list.push_back(dir1);
  dir_list.push_back(dir2);
  std::vector <TString> dirlabel_list;
  dirlabel_list.push_back(label_dir1);
  dirlabel_list.push_back(label_dir2);
  do_the_comparison(dir_list,dirlabel_list,algo,dir_prefix,binning_select,pub_style);


}

void compare_Systematics(TString dir1, TString dir2, TString dir3, TString label_dir1, TString label_dir2, TString label_dir3, TString algo, TString dir_prefix,TString binning_select , Bool_t pub_style){

  std::vector <TString> dir_list;
  dir_list.push_back(dir1);
  dir_list.push_back(dir2);
  dir_list.push_back(dir3);
  std::vector <TString> dirlabel_list;
  dirlabel_list.push_back(label_dir1);
  dirlabel_list.push_back(label_dir2);
  dirlabel_list.push_back(label_dir3);
  do_the_comparison(dir_list,dirlabel_list,algo,dir_prefix,binning_select,pub_style);

}

void compare_Systematics(TString dir1, TString dir2, TString dir3, TString dir4, TString label_dir1, TString label_dir2, TString label_dir3, TString label_dir4, TString algo, TString dir_prefix,TString binning_select, Bool_t pub_style){
  std::vector <TString> dir_list;
  dir_list.push_back(dir1);
  dir_list.push_back(dir2);
  dir_list.push_back(dir3);
  dir_list.push_back(dir4);
  std::vector <TString> dirlabel_list;
  dirlabel_list.push_back(label_dir1);
  dirlabel_list.push_back(label_dir2);
  dirlabel_list.push_back(label_dir3);
  dirlabel_list.push_back(label_dir4);
  do_the_comparison(dir_list,dirlabel_list,algo,dir_prefix,binning_select,pub_style);


}

void compare_Systematics(TString dir1, TString dir2, TString dir3, TString dir4, TString dir5, TString label_dir1, TString label_dir2, TString label_dir3, TString label_dir4, TString label_dir5, TString algo, TString dir_prefix,TString binning_select, Bool_t pub_style){
  std::vector <TString> dir_list;
  dir_list.push_back(dir1);
  dir_list.push_back(dir2);
  dir_list.push_back(dir3);
  dir_list.push_back(dir4);
  dir_list.push_back(dir5);
  std::vector <TString> dirlabel_list;
  dirlabel_list.push_back(label_dir1);
  dirlabel_list.push_back(label_dir2);
  dirlabel_list.push_back(label_dir3);
  dirlabel_list.push_back(label_dir4);
  dirlabel_list.push_back(label_dir5);
  do_the_comparison(dir_list,dirlabel_list,algo,dir_prefix,binning_select,pub_style);
}

void compare_Systematics(TString dir1, TString dir2, TString dir3, TString dir4, TString dir5, TString dir6, TString label_dir1, TString label_dir2, TString label_dir3, TString label_dir4, TString label_dir5, TString label_dir6, TString algo, TString dir_prefix,TString binning_select, Bool_t pub_style){
  std::vector <TString> dir_list;
  dir_list.push_back(dir1);
  dir_list.push_back(dir2);
  dir_list.push_back(dir3);
  dir_list.push_back(dir4);
  dir_list.push_back(dir5);
  dir_list.push_back(dir6);
  std::vector <TString> dirlabel_list;
  dirlabel_list.push_back(label_dir1);
  dirlabel_list.push_back(label_dir2);
  dirlabel_list.push_back(label_dir3);
  dirlabel_list.push_back(label_dir4);
  dirlabel_list.push_back(label_dir5);
  dirlabel_list.push_back(label_dir6);
  do_the_comparison(dir_list,dirlabel_list,algo,dir_prefix,binning_select,pub_style);
}
