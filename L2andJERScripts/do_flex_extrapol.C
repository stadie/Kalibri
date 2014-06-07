#include <iostream>
#include "TMath.h"
#include "TString.h"
#include "TFile.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TImage.h"
#include "TGraphErrors.h"
#include "TH1.h"
#include "THelpers_2.h"
#include <fstream>
#include "extrapol_helpers.h"
#include <cmath>
#include "tdrstyle_mod.C"
#include "TFitResult.h"
#include "TFitResultPtr.h"
#include "TMatrixD.h"
#include "do_flex_extrapol.h"



Double_t get_error(Double_t x, TMatrixDSym cov, TF1 *kFSR_fit){

  //    Double_t par_a=kFSR_fit->GetParameter(0);
    Double_t par_b=kFSR_fit->GetParameter(1);
    Double_t par_c=kFSR_fit->GetParameter(2);

    Double_t dqdz [3];
    dqdz[0]=1;
    dqdz[1]=cosh(x)/(1+par_c*cosh(x));
    dqdz[2]=(-par_b*cosh(x)*cosh(x))/TMath::Power((1+par_c*cosh(x)),2);

    TMatrixD result(1,1);
    TMatrixD part_left(1,3);
    TMatrixD part_right(3,1);
    for(Int_t j=0;j<3;j++){
      part_left[0][j]=dqdz[j];
      part_right[j][0]=dqdz[j];
    }
    cov.Print();
    part_left.Print();
    part_right.Print();

    part_right=cov*part_right;
    result= part_left*part_right;
    Double_t ey = TMath::Sqrt(result[0][0]);

    return ey;

}


//! Central method to call for determination of residuals
//!
//! Includes all neccessary steps to determine the residuals such as extrapolations for the radiation correction
//! and can be adapted by a large number of parameters. Plots and correction .txt files are saved in sub-folders 
//! and root-files in accordance to the different parameters. 
//!
//!
//!
//! Parameters:
//!      jet_type: PF, Calo, JPT, PFchs, ...
//!      generatorone_/generatortwo_:
//!        - Allows to compare different generators (not used by batch steering script) 
//!          directly within the script (you will find e.g. res1 and res2-plots)
//!          Sample use:
//!           You have two subfolders in your CalobCore/SAMPLENAME directory 
//!           dijetsFall10_TunePYTHIA_AK5PF_weighted_residuals_k_HFfix --> generatorone = "TunePYTHIA"
//!           dijetsFall10_TuneHerwig_AK5PF_weighted_residuals_k_HFfix --> generatortwo = "TuneHerwig"
//!          The script will then create almost all plots in parallel and creat comparison plots
//!      image_ext: can e.g. be .eps or .pdf (eps recommmended by CMS), if image_ext="" --> no export to individual files
//!      root_export: decide whether plots and (important when scripting) extrapolations/krads are saved to root-files
//!      use_imported_kFSRAbs_: switch whether or not to use an imported kRad
//!      fine_coarse: as of now mainly a switch between "kostas" and "k_HFfix"-binning, but other binnings are still available
//!      use_easy_mean: switch that decides whether "OneBin"/pt-independent histograms are read in or the pt-dependent variants are chosen
//!      use_fitted_kFSR: switch whether or not to use the fit to the used krad (either imported or from the same run of the script, depending on use_imported_kFSRAbs_)
//!      corr_generation_: SAMPLENAME of the path to the CalobCore/SAMPLENAME directory, e.g. 2011Full2011_CORRF11DB_He_AK5_MC_F11Z2wPUsm_Y_f_kostas_MPF_AK5
//!      ratio_of_mean_or_GM_: switch whether to use the mean, the Gaussian mean or the interquartile mean (where the IQmean needs a very fine binning of the response to work)
//!      export_all_plots_: switch that saves some space if set to false
//!      kFSR_eq_one_: switch to force kRAD-correction to be =1 for all eta
//!      MPF_or_rel_response_: switch for choosing between relative response and MPF-response histograms when reading in.
//!
//!
//! Interesting plots that are good to know that they are there:
//!      - (+eta)/(-eta)/|eta| comparison plots
//!      - plots with response vs. no. of reconstructed vertices (needs coarsely binned histogram to give useful information)
//!
//!
//!
//!
//!
// ----------------------------------------------------------------   
void do_flex_extrapol::Loop(TString jet_type_, TString generatorone_, TString generatortwo_, TString image_ext_, TString root_export_, TString use_imported_kFSRAbs_, TString fine_coarse_, TString use_easy_mean_, TString use_fitted_kFSR_, TString corr_generation_, TString ratio_of_mean_or_GM_, Bool_t export_all_plots_, TString kFSR_eq_one_, TString MPF_or_rel_response_){


jet_type 	      =  jet_type_		;              
generatorone 	      =  generatorone_	        ;              
generatortwo 	      =  generatortwo_	        ;              
image_ext 	      =  image_ext_		;              
root_export 	      =  root_export_		;              
use_imported_kFSRAbs  =  use_imported_kFSRAbs_  ;              
fine_coarse 	      =  fine_coarse_		;              
use_easy_mean 	      =  use_easy_mean_	        ;              
use_fitted_kFSR       =	 use_fitted_kFSR_	;              
corr_generation       =	 corr_generation_	;              
ratio_of_mean_or_GM   =  ratio_of_mean_or_GM_	;              
export_all_plots      =  export_all_plots_	;              
kFSR_eq_one           =  kFSR_eq_one_           ;              
MPF_or_rel_response   =  MPF_or_rel_response_   ;

  Bool_t no_entries_went_wrong=false;

  setTDRStyle();
     if(chdir(corr_generation) != 0){ 
       mkdir(corr_generation, S_IRWXU|S_IRWXG|S_IRWXO); 
       chdir(corr_generation); 
     } 



  TString easy_mean_prefix="";
  if(use_easy_mean.Contains("use_easy_mean"))
    easy_mean_prefix="OneBin";

  if(use_imported_kFSRAbs.Contains("true"))do_flex_extrapol::import_plots();
  cov.Print();


  TCanvas* c = new TCanvas("c","",600,600);
//  c->SetRightMargin(0.04);
//  c->SetTopMargin(0.099);
//  c->SetBottomMargin(0.16);

  do_flex_extrapol::define_cosmetics_and_cuts();

  TH1D* dummy_histo = new TH1D("dummy_histo", "dummy_histo", 10, 20, 2000);
  for(Int_t i = 0; i<10;i++){
    dummy_histo->SetBinContent(i,1.0);
    dummy_histo->SetBinError(i,1.0);
  }

  std::vector <TString> faulty_values_dummy_histos;


  TString dir_prefix("../../"+corr_generation+"/");//../../2011_01_new_Kalibri_L2L3res_and_JWPF/scripts/");

  TString Residual_plots_one(dir_prefix + "dijetsFall10_"+generatorone+"_AK5"+jet_type+"_weighted_residuals_"+ fine_coarse+"/plots/KalibriPlots.root");
  TString Residual_plots_two(dir_prefix + "dijetsFall10_"+generatortwo+"_AK5"+jet_type+"_weighted_residuals_"+ fine_coarse+"/plots/KalibriPlots.root");


  do_flex_extrapol::define_eta_bins_and_labels();
cout << "labels defined" << endl;
  TString XVsPt, XVsPtType;
  if(MPF_or_rel_response.Contains("rel_response")){
    ratio_of_mean_or_GM="_RatioOf"+ ratio_of_mean_or_GM+"s";
    XVsPt = "AsymmetryVsPt";
    XVsPtType = "AsymmetryVsMeanPt";
  }
  else if(MPF_or_rel_response.Contains("MPF")){
    ratio_of_mean_or_GM="_"+ ratio_of_mean_or_GM;
    XVsPt = "MPFT1VsPtAve";
    XVsPtType = "MPFMETT1ResponseVsMeanPt";
  }
cout << "open one" << endl;
  TFile* fone = TFile::Open(Residual_plots_one);
cout << "open two" << endl;
  TFile* ftwo = TFile::Open(Residual_plots_two);
cout << "Zombie 1" << endl;
    if (fone->IsZombie()) {
       cout << "Error opening file" << endl;
       exit(-1);
    }
    else cout << "ROOT-Datei erfolgreich geladen. " << endl;

    if (ftwo->IsZombie()) {
       cout << "Error opening file" << endl;
       exit(-1);
    }
    else cout << "ROOT-Datei erfolgreich geladen. " << endl;


  if(use_imported_kFSRAbs.Contains("true"))fine_coarse=fine_coarse+"_use_coarse_kFSRAbs";
  if(use_fitted_kFSR.Contains("use_fitted_kFSR"))fine_coarse=fine_coarse+"_fit";
  if(use_easy_mean.Contains("use_easy_mean"))fine_coarse=fine_coarse+"_"+use_easy_mean;
  if(kFSR_eq_one.Contains("kFSR_eq_one"))fine_coarse="kFSR_eq_one_"+fine_coarse;


TH1D* MCone_eta_spectrum_MC =  (TH1D*)fone->Get("AsymmetryVsEta/AsymmetryVsEta_MC_EtaSpectrum");
TH1D* MCone_eta_spectrum_data =  (TH1D*)fone->Get("AsymmetryVsEta/AsymmetryVsEta_data_EtaSpectrum");
TH1D* MCtwo_eta_spectrum_MC =  (TH1D*)ftwo->Get("AsymmetryVsEta/AsymmetryVsEta_MC_EtaSpectrum");
TH1D* MCtwo_eta_spectrum_data =  (TH1D*)ftwo->Get("AsymmetryVsEta/AsymmetryVsEta_data_EtaSpectrum");


TH1D* MCone_pt_spectrum_MC =  (TH1D*)fone->Get("AsymmetryVsPt20/AsymmetryVsPt20_MC_MeanPtSpectrum");
TH1D* MCone_pt_spectrum_data =  (TH1D*)fone->Get("AsymmetryVsPt20/AsymmetryVsPt20_data_MeanPtSpectrum");
TH1D* MCtwo_pt_spectrum_MC =  (TH1D*)ftwo->Get("AsymmetryVsPt20/AsymmetryVsPt20_MC_MeanPtSpectrum");
TH1D* MCtwo_pt_spectrum_data =  (TH1D*)ftwo->Get("AsymmetryVsPt20/AsymmetryVsPt20_data_MeanPtSpectrum");

  cout << Residual_plots_one << endl;

  std::vector < std::vector <TH1D*> > all_ptthree_all_eta_MC1_L2L3_;	 
  std::vector < std::vector <TH1D*> > all_ptthree_all_eta_MC2_L2L3_;	 
  std::vector < std::vector <TH1D*> > all_ptthree_all_eta_D_L2L3_;	 
  std::vector < std::vector <TH1D*> > all_ptthree_all_eta_D_L2L3res_;	 
  std::vector < std::vector <TH1D*> > all_ptthree_all_eta_ratio_val1_;	 
  std::vector < std::vector <TH1D*> > all_ptthree_all_eta_ratio_val2_;	 
  std::vector < std::vector <TH1D*> > all_ptthree_all_eta_ratio_res1_;	 
  std::vector < std::vector <TH1D*> > all_ptthree_all_eta_ratio_res2_;	 

  std::vector < std::vector <TH1D*> > all_ptthree_all_Abseta_MC1_L2L3_;	 
  std::vector < std::vector <TH1D*> > all_ptthree_all_Abseta_MC2_L2L3_;	 
  std::vector < std::vector <TH1D*> > all_ptthree_all_Abseta_D_L2L3_;	 
  std::vector < std::vector <TH1D*> > all_ptthree_all_Abseta_D_L2L3res_;	 
  std::vector < std::vector <TH1D*> > all_ptthree_all_Abseta_ratio_val1_;	 
  std::vector < std::vector <TH1D*> > all_ptthree_all_Abseta_ratio_val2_;	 
  std::vector < std::vector <TH1D*> > all_ptthree_all_Abseta_ratio_res1_;	 
  std::vector < std::vector <TH1D*> > all_ptthree_all_Abseta_ratio_res2_;	 

  std::vector < std::vector <std::pair < Double_t, Double_t> > > all_ptthree_all_Abseta_res1_ptreach_;	 

  for(unsigned int cut_i=0;cut_i<ptthreecuts.size();cut_i++){
    //cout << ptthreecuts[cut_i] << endl;

    std::vector <TH1D*> all_eta_MC1_L2L3_;	 
    std::vector <TH1D*> all_eta_MC2_L2L3_;	 
    std::vector <TH1D*> all_eta_D_L2L3_;	 
    std::vector <TH1D*> all_eta_D_L2L3res_;	 
    std::vector <TH1D*> all_eta_ratio_val1_;	 
    std::vector <TH1D*> all_eta_ratio_val2_;	 
    std::vector <TH1D*> all_eta_ratio_res1_;	 
    std::vector <TH1D*> all_eta_ratio_res2_;	 


    std::vector <TH1D*> all_Abseta_MC1_L2L3_;	 
    std::vector <TH1D*> all_Abseta_MC2_L2L3_;	 
    std::vector <TH1D*> all_Abseta_D_L2L3_;	 
    std::vector <TH1D*> all_Abseta_D_L2L3res_;	 
    std::vector <TH1D*> all_Abseta_ratio_val1_;	 
    std::vector <TH1D*> all_Abseta_ratio_val2_;	 
    std::vector <TH1D*> all_Abseta_ratio_res1_;	 
    std::vector <TH1D*> all_Abseta_ratio_res2_;	 

    std::vector <std::pair < Double_t, Double_t> > all_Abseta_res1_ptreach_;
    for(Int_t eta_i=0;eta_i<no_eta_bins;eta_i++){


      //      MyClass *obj;
      //      directory->GetObject("some object",obj);
      //      if (obj) { ... the object exist and inherits from MyClass ... }

      TH1D* hMC1;
      fone->GetObject((easy_mean_prefix+XVsPt+ptthreecuts[cut_i])+("/"+easy_mean_prefix+XVsPt+ptthreecuts[cut_i])+"_"+XVsPtType+"_MC_L2L3_Eta"+(Long_t)eta_i+ratio_of_mean_or_GM,hMC1);
      if(!hMC1)cout << "FEHLER - NICHT IMPORTIERT!" << endl;
      TH1D* hMC2;
      ftwo->GetObject((easy_mean_prefix+XVsPt+ptthreecuts[cut_i])+("/"+easy_mean_prefix+XVsPt+ptthreecuts[cut_i])+"_"+XVsPtType+"_MC_L2L3_Eta"+(Long_t)eta_i+ratio_of_mean_or_GM,hMC2);
      if(!hMC2)cout << "FEHLER - NICHT IMPORTIERT!" << endl;
      TH1D* hL2L3;
      fone->GetObject((easy_mean_prefix+XVsPt+ptthreecuts[cut_i])+("/"+easy_mean_prefix+XVsPt+ptthreecuts[cut_i])+"_"+XVsPtType+"_data_L2L3_Eta"+(Long_t)eta_i+ratio_of_mean_or_GM,hL2L3);
      if(!hL2L3)cout << "FEHLER - NICHT IMPORTIERT!" << endl;
      TH1D* hL2L3res;
      fone->GetObject((easy_mean_prefix+XVsPt+ptthreecuts[cut_i])+("/"+easy_mean_prefix+XVsPt+ptthreecuts[cut_i])+"_"+XVsPtType+"_data_L2L3res_Eta"+(Long_t)eta_i+ratio_of_mean_or_GM,hL2L3res);
      if(!hL2L3res)cout << "FEHLER - NICHT IMPORTIERT!" << endl;


      
      if ( hMC1->GetEntries()==0||hMC2->GetEntries()==0||hL2L3->GetEntries()==0||hL2L3res->GetEntries()==0){
	std::cout << "DAS IST NICHT GUT!!! FEHLER!!!" << std::endl;
	//	hMC1->Dump();
      	std::cout << (easy_mean_prefix+XVsPt+ptthreecuts[cut_i])+("/"+easy_mean_prefix+XVsPt+ptthreecuts[cut_i])+"_"+XVsPtType+"_MC_L2L3_Eta"+(Long_t)eta_i+ratio_of_mean_or_GM << std::endl;


	faulty_values_dummy_histos.push_back("hMC-hL2L3: " +(easy_mean_prefix+XVsPt+ptthreecuts[cut_i])+("/"+easy_mean_prefix+XVsPt+ptthreecuts[cut_i])+"_"+XVsPtType+"_MC_L2L3_Eta"+(Long_t)eta_i+ratio_of_mean_or_GM+  " hMC1: " +  (Long_t) hMC1->GetEntries() + " hMC2: " + (Long_t) hMC2->GetEntries()+ " hL2L3: " + (Long_t) hL2L3->GetEntries()+ " hL2L3res: " + (Long_t) hL2L3res->GetEntries());	
	hMC1=dummy_histo;
	hMC2=dummy_histo;
	hL2L3=dummy_histo;
	hL2L3res=dummy_histo;

	no_entries_went_wrong=true;
      }

//      if(all_eta_MC1_L2L3_.back()->GetEntries()==0)std::cout << "this cannot work..." << std::endl;
//      if(all_eta_MC2_L2L3_.back()->GetEntries()==0)std::cout << "this cannot work..." << std::endl;
//      if(all_eta_D_L2L3_.back()->GetEntries()==0)std::cout << "this cannot work..." << std::endl;
//      if(all_eta_D_L2L3res_.back()->GetEntries()==0)std::cout << "this cannot work..." << std::endl;

      all_eta_MC1_L2L3_.push_back(hMC1);
      all_eta_MC2_L2L3_.push_back(hMC2);
      all_eta_D_L2L3_.push_back(hL2L3);
      all_eta_D_L2L3res_.push_back(hL2L3res);

      //      cout << "test" << endl;
      all_eta_MC1_L2L3_ .back()->GetYaxis()->SetRangeUser(0.5,1.5);
      all_eta_MC2_L2L3_ .back()->GetYaxis()->SetRangeUser(0.5,1.5);
      all_eta_D_L2L3_   .back()->GetYaxis()->SetRangeUser(0.5,1.5);
      all_eta_D_L2L3res_.back()->GetYaxis()->SetRangeUser(0.5,1.5);
      //      cout << "test2" << endl;


      TH1D* temp_res1 = (TH1D*) all_eta_D_L2L3res_.back()->Clone();
      TH1D* temp_val1 = (TH1D*) all_eta_D_L2L3_.back()->Clone();
      TH1D* temp_res2 = (TH1D*) all_eta_D_L2L3res_.back()->Clone();
      TH1D* temp_val2 = (TH1D*) all_eta_D_L2L3_.back()->Clone();

      temp_res1->Sumw2();
      temp_val1->Sumw2();
      temp_res2->Sumw2();
      temp_val2->Sumw2();
      temp_val1->Divide(all_eta_MC1_L2L3_.back(),all_eta_D_L2L3res_.back());
      temp_res1->Divide(all_eta_MC1_L2L3_.back(),all_eta_D_L2L3_.back());
      temp_val2->Divide(all_eta_MC2_L2L3_.back(),all_eta_D_L2L3res_.back());
      temp_res2->Divide(all_eta_MC2_L2L3_.back(),all_eta_D_L2L3_.back());

      TH1style(temp_val1,line_styles_, colours_, markers_, 0,0,0,"R(MC)/R(data)");
      TH1style(temp_val2,line_styles_, colours_, markers_, 1,1,1,"R(MC)/R(data)");
      TH1style(temp_res1,line_styles_, colours_, markers_, 2,2,2,"R(MC)/R(data)");
      TH1style(temp_res2,line_styles_, colours_, markers_, 3,3,3,"R(MC)/R(data)");

      TF1 *fit_const = new TF1("fit_const","[0]",temp_val1->GetXaxis()->GetXmin(),temp_val1->GetXaxis()->GetXmax()); //was used before...
      fit_const->SetParameters(1,1);
      fit_const->SetParName(0,"const");

      TF1 *fit_loglin = new TF1("fit_loglin","[0]+[1]*TMath::Log(x)",temp_val1->GetXaxis()->GetXmin(),temp_val1->GetXaxis()->GetXmax()); //was used before...
      fit_loglin->SetParameters(1,1);
      fit_loglin->SetParName(0,"const");
      fit_loglin->SetParName(1,"slope");
      //      fit_loglin->SetParName(2,"x0");


      if (  temp_res1->GetEntries()==0|| temp_res2->GetEntries()==0|| temp_val1->GetEntries()==0||temp_val2->GetEntries()==0){
	std::cout << "DAS IST NICHT GUT!!! FEHLER!!!" << std::endl;
	//	hMC1->Dump();
	faulty_values_dummy_histos.push_back("divided (ratios): " +(easy_mean_prefix+XVsPt+ptthreecuts[cut_i])+("/"+easy_mean_prefix+XVsPt+ptthreecuts[cut_i])+"_"+XVsPtType+"_MC_L2L3_Eta"+(Long_t)eta_i+ratio_of_mean_or_GM+ " temp_res1: " +  (Long_t)  temp_res1->GetEntries()+ " temp_val1: " +  (Long_t)  temp_val1->GetEntries()+ " temp_res2: " +  (Long_t)  temp_res2->GetEntries()+ " temp_val2: " +  (Long_t)  temp_val2->GetEntries());
	temp_res1 = dummy_histo;
	temp_val1 = dummy_histo;
	temp_res2 = dummy_histo;
	temp_val2 = dummy_histo;
      }

      cout << hMC1->GetName() << endl;
      cout << "eta: " << eta_i << " and test " <<  eta_bins_labels_[eta_i].first << " " << eta_bins_labels_[eta_i].second << endl;
      fit_const->SetLineColor(temp_res1->GetLineColor());
      fit_const->SetLineStyle(temp_res1->GetLineStyle());
      temp_res1->Fit("fit_const","","same");
      temp_res1->Fit("fit_loglin","+","same");

      cout << "eta: " << eta_i << " and test " <<  eta_bins_labels_[eta_i].first << " " << eta_bins_labels_[eta_i].second << endl;
      fit_const->SetLineColor(temp_res2->GetLineColor());
      fit_const->SetLineStyle(temp_res2->GetLineStyle());
      temp_res2->Fit("fit_const","","same");

      cout << "eta: " << eta_i << " and test " <<  eta_bins_labels_[eta_i].first << " " << eta_bins_labels_[eta_i].second << endl;
      fit_const->SetLineColor(temp_val1->GetLineColor());
      fit_const->SetLineStyle(temp_val1->GetLineStyle());
      temp_val1->Fit("fit_const","","same");

      cout << "eta: " << eta_i << " and test " <<  eta_bins_labels_[eta_i].first << " " << eta_bins_labels_[eta_i].second << endl;
      fit_const->SetLineColor(temp_val2->GetLineColor());
      fit_const->SetLineStyle(temp_val2->GetLineStyle());
      temp_val2->Fit("fit_const","","same");



      all_eta_ratio_val1_.push_back(temp_val1);
      all_eta_ratio_val2_.push_back(temp_val2);
      all_eta_ratio_res1_.push_back(temp_res1);
      all_eta_ratio_res2_.push_back(temp_res2);

    }

    std::vector <std::vector<TH1D* > > overlay_eta_valres1_;
    overlay_eta_valres1_.push_back(all_eta_ratio_val1_);
    for(unsigned int style_i=0;style_i<all_eta_ratio_val1_.size();style_i++)TH1style_plus_fit(all_eta_ratio_val1_[style_i],line_styles_, colours_, markers_, 0,0,0,"R(MC)/R(data)","fit_const");
    overlay_eta_valres1_.push_back(all_eta_ratio_res1_);
    for(unsigned int style_i=0;style_i<all_eta_ratio_res1_.size();style_i++)TH1style_plus_fit(all_eta_ratio_res1_[style_i],line_styles_, colours_, markers_, 1,1,1,"R(MC)/R(data)","fit_const");
    
    TLegend *leg_valres1;
    leg_valres1 = new TLegend(0.25,0.70,0.55,0.85);
    leg_valres1->SetFillColor(kWhite);
    leg_valres1->SetFillStyle(0);
      //   leg->SetHeader("Legende");
    leg_valres1->AddEntry(all_eta_ratio_val1_.back(),"r("+generatorone+")/r(L2L3res)","lep");
    leg_valres1->AddEntry(all_eta_ratio_res1_.back(),"r("+generatorone+")/r(L2L3)","lep");

    if(export_all_plots)draw_Overlay_TH1D_save_PS(leg_valres1,overlay_eta_valres1_, "Overlay_Eta_"+generatorone+"_"+generatortwo+"_res1_and_val1_vs_pt_with_TJF_"+ptthreecuts[cut_i], "nice", "legend", "x1_y0_z0", 0,-1, 0.85, 1.35, 0,-1, image_ext, fine_coarse+"_"+ generatorone+"_"+generatortwo+"_"+jet_type+"_PS",1.);

    std::vector <std::vector<TH1D* > > overlay_eta_valres2_;
    overlay_eta_valres2_.push_back(all_eta_ratio_val2_);
    for(unsigned int style_i=0;style_i<all_eta_ratio_val2_.size();style_i++)TH1style_plus_fit(all_eta_ratio_val2_[style_i],line_styles_, colours_, markers_, 0,0,0,"R(MC)/R(data)","fit_const");
    overlay_eta_valres2_.push_back(all_eta_ratio_res2_);
    for(unsigned int style_i=0;style_i<all_eta_ratio_res2_.size();style_i++)TH1style_plus_fit(all_eta_ratio_res2_[style_i],line_styles_, colours_, markers_, 1,1,1,"R(MC)/R(data)","fit_const");
    
    TLegend *leg_valres2;
    leg_valres2 = new TLegend(0.25,0.70,0.55,0.85);
    leg_valres2->SetFillColor(kWhite);
    leg_valres2->SetFillStyle(0);
      //   leg->SetHeader("Legende");
    leg_valres2->AddEntry(all_eta_ratio_val2_.back(),"r("+generatortwo+")/r(L2L3res)","lep");
    leg_valres2->AddEntry(all_eta_ratio_res2_.back(),"r("+generatortwo+")/r(L2L3)","lep");
    
    // save memory
    //    if(export_all_plots)draw_Overlay_TH1D_save_PS(leg_valres2,overlay_eta_valres2_, "Overlay_Eta_"+generatorone+"_"+generatortwo+"_res2_and_val2_vs_pt_with_TJF_"+ptthreecuts[cut_i], "nice", "legend", "x1_y0_z0", 0,-1, 0.85, 1.35, 0,-1, image_ext, fine_coarse+"_"+ generatorone+"_"+generatortwo+"_"+jet_type+"_PS",1.);


    std::vector <std::vector<TH1D* > > overlay_eta_res1res2_;
    overlay_eta_res1res2_.push_back(all_eta_ratio_res1_);
    for(unsigned int style_i=0;style_i<all_eta_ratio_res1_.size();style_i++)TH1style_plus_fit(all_eta_ratio_val2_[style_i],line_styles_, colours_, markers_, 0,0,0,"R(MC)/R(data)","fit_const");
    overlay_eta_res1res2_.push_back(all_eta_ratio_res2_);
    for(unsigned int style_i=0;style_i<all_eta_ratio_res2_.size();style_i++)TH1style_plus_fit(all_eta_ratio_res2_[style_i],line_styles_, colours_, markers_, 1,1,1,"R(MC)/R(data)","fit_const");
    
    TLegend *leg_res1res2;
    leg_res1res2 = new TLegend(0.25,0.70,0.55,0.85);
    leg_res1res2->SetFillColor(kWhite);
    leg_res1res2->SetFillStyle(0);
      //   leg->SetHeader("Legende");
    leg_res1res2->AddEntry(all_eta_ratio_val2_.back(),"r("+generatortwo+")/r(L2L3res)","lep");
    leg_res1res2->AddEntry(all_eta_ratio_res2_.back(),"r("+generatortwo+")/r(L2L3)","lep");
    cout << "hier noch o.k." << endl;
    // save memory
    //if(export_all_plots)draw_Overlay_TH1D_save_PS(leg_res1res2,overlay_eta_res1res2_, "Overlay_Eta_"+generatorone+"_"+generatortwo+"_res1_and_res2_vs_pt_with_TJF_"+ptthreecuts[cut_i], "nice", "legend", "x1_y0_z0", 0,-1, 0.85, 1.35, 0,-1, image_ext, fine_coarse+"_"+ generatorone+"_"+generatortwo+"_"+jet_type+"_PS",1.);
    cout << "hier noch o.k." << endl;

 all_ptthree_all_eta_MC1_L2L3_      .push_back(all_eta_MC1_L2L3_   );	        	
 all_ptthree_all_eta_MC2_L2L3_      .push_back(all_eta_MC2_L2L3_   );	        	
 all_ptthree_all_eta_D_L2L3_        .push_back(all_eta_D_L2L3_ 	   );   	
 all_ptthree_all_eta_D_L2L3res_     .push_back(all_eta_D_L2L3res_  );	         	
 all_ptthree_all_eta_ratio_val1_    .push_back(all_eta_ratio_val1_ );	          	
 all_ptthree_all_eta_ratio_val2_    .push_back(all_eta_ratio_val2_ );	          	
 all_ptthree_all_eta_ratio_res1_    .push_back(all_eta_ratio_res1_ );	          	
 all_ptthree_all_eta_ratio_res2_    .push_back(all_eta_ratio_res2_ );

    cout << "hier noch o.k." << endl;



    for(Int_t Abseta_i=0;Abseta_i<no_Abseta_bins;Abseta_i++){

     //      directory->GetObject("some object",obj);
      //      if (obj) { ... the object exist and inherits from MyClass ... }

      cout <<
	(easy_mean_prefix+"Abs"+XVsPt+ptthreecuts[cut_i])+("/"+easy_mean_prefix+"Abs"+XVsPt+ptthreecuts[cut_i])+"_"+XVsPtType+"_MC_L2L3_AbsEta"+(Long_t)Abseta_i+ratio_of_mean_or_GM
	 << endl;

      TH1D* hMC1;
      fone->GetObject((easy_mean_prefix+"Abs"+XVsPt+ptthreecuts[cut_i])+("/"+easy_mean_prefix+"Abs"+XVsPt+ptthreecuts[cut_i])+"_"+XVsPtType+"_MC_L2L3_AbsEta"+(Long_t)Abseta_i+ratio_of_mean_or_GM,hMC1);
      if(!hMC1)cout << "FEHLER - NICHT IMPORTIERT!" << endl;
      TH1D* hMC2;
      ftwo->GetObject((easy_mean_prefix+"Abs"+XVsPt+ptthreecuts[cut_i])+("/"+easy_mean_prefix+"Abs"+XVsPt+ptthreecuts[cut_i])+"_"+XVsPtType+"_MC_L2L3_AbsEta"+(Long_t)Abseta_i+ratio_of_mean_or_GM,hMC2);
      if(!hMC2)cout << "FEHLER - NICHT IMPORTIERT!" << endl;
      TH1D* hL2L3;
      fone->GetObject((easy_mean_prefix+"Abs"+XVsPt+ptthreecuts[cut_i])+("/"+easy_mean_prefix+"Abs"+XVsPt+ptthreecuts[cut_i])+"_"+XVsPtType+"_data_L2L3_AbsEta"+(Long_t)Abseta_i+ratio_of_mean_or_GM,hL2L3);
      if(!hL2L3)cout << "FEHLER - NICHT IMPORTIERT!" << endl;
      TH1D* hL2L3res;
      fone->GetObject((easy_mean_prefix+"Abs"+XVsPt+ptthreecuts[cut_i])+("/"+easy_mean_prefix+"Abs"+XVsPt+ptthreecuts[cut_i])+"_"+XVsPtType+"_data_L2L3res_AbsEta"+(Long_t)Abseta_i+ratio_of_mean_or_GM,hL2L3res);
      if(!hL2L3res)cout << "FEHLER - NICHT IMPORTIERT!" << endl;


      
      if ( hMC1->GetEntries()==0|| hMC2->GetEntries()==0||hL2L3->GetEntries()==0||hL2L3res->GetEntries()==0){
	std::cout << "DAS IST NICHT GUT!!! FEHLER!!!" << std::endl;
	//	hMC1->Dump();

	faulty_values_dummy_histos.push_back("hMC-hL2L3: " +(easy_mean_prefix+"Abs"+XVsPt+ptthreecuts[cut_i])+("/"+easy_mean_prefix+"Abs"+XVsPt+ptthreecuts[cut_i])+"_"+XVsPtType+"_MC_L2L3_AbsEta"+(Long_t)Abseta_i+ratio_of_mean_or_GM +  " hMC1: " +  (Long_t) hMC1->GetEntries() + " hMC2: " + (Long_t) hMC2->GetEntries()+ " hL2L3: " + (Long_t) hL2L3->GetEntries()+ " hL2L3res: " + (Long_t) hL2L3res->GetEntries());
	hMC1=dummy_histo;
	hMC2=dummy_histo;
	hL2L3=dummy_histo;
	hL2L3res=dummy_histo;
	no_entries_went_wrong=true;
      }
    cout << "hier noch o.k. Abs" << endl;

      cout << ((easy_mean_prefix+"Abs"+XVsPt+ptthreecuts[cut_i])+("/AbsAsymmetryVsPt"+ptthreecuts[cut_i])+"_"+XVsPtType+"_MC_L2L3_AbsEta"+(Long_t)Abseta_i+ratio_of_mean_or_GM) << endl;
      all_Abseta_MC1_L2L3_.push_back(hMC1);
      all_Abseta_MC2_L2L3_.push_back(hMC2);
      all_Abseta_D_L2L3_.push_back(hL2L3);
      all_Abseta_D_L2L3res_.push_back(hL2L3res);

      all_Abseta_MC1_L2L3_ .back()->GetYaxis()->SetRangeUser(0.5,1.5);
      all_Abseta_MC2_L2L3_ .back()->GetYaxis()->SetRangeUser(0.5,1.5);
      all_Abseta_D_L2L3_   .back()->GetYaxis()->SetRangeUser(0.5,1.5);
      all_Abseta_D_L2L3res_.back()->GetYaxis()->SetRangeUser(0.5,1.5);


      TH1D* temp_res1 = (TH1D*) all_Abseta_D_L2L3res_.back()->Clone();
      TH1D* temp_val1 = (TH1D*) all_Abseta_D_L2L3_.back()->Clone();
      TH1D* temp_res2 = (TH1D*) all_Abseta_D_L2L3res_.back()->Clone();
      TH1D* temp_val2 = (TH1D*) all_Abseta_D_L2L3_.back()->Clone();

      temp_res1->Sumw2();
      temp_val1->Sumw2();
      temp_res2->Sumw2();
      temp_val2->Sumw2();

      temp_val1->Divide(all_Abseta_MC1_L2L3_.back(),all_Abseta_D_L2L3res_.back());
      temp_res1->Divide(all_Abseta_MC1_L2L3_.back(),all_Abseta_D_L2L3_.back());
      temp_val2->Divide(all_Abseta_MC2_L2L3_.back(),all_Abseta_D_L2L3res_.back());
      temp_res2->Divide(all_Abseta_MC2_L2L3_.back(),all_Abseta_D_L2L3_.back());

      TH1style(temp_val1,line_styles_, colours_, markers_, 0,0,0,"R(MC)/R(data)");
      TH1style(temp_val2,line_styles_, colours_, markers_, 1,1,1,"R(MC)/R(data)");
      TH1style(temp_res1,line_styles_, colours_, markers_, 2,2,2,"R(MC)/R(data)");
      TH1style(temp_res2,line_styles_, colours_, markers_, 3,3,3,"R(MC)/R(data)");

      if (  temp_res1->GetEntries()==0|| temp_res2->GetEntries()==0|| temp_val1->GetEntries()==0||temp_val2->GetEntries()==0){
	std::cout << "DAS IST NICHT GUT!!! FEHLER!!!" << std::endl;
	//	hMC1->Dump();
	faulty_values_dummy_histos.push_back("divided (ratios): " +(easy_mean_prefix+"Abs"+XVsPt+ptthreecuts[cut_i])+("/"+easy_mean_prefix+"Abs"+XVsPt+ptthreecuts[cut_i])+"_"+XVsPtType+"_MC_L2L3_AbsEta"+(Long_t)Abseta_i+ratio_of_mean_or_GM+  " temp_res1: " +  (Long_t)  temp_res1->GetEntries()+ " temp_val1: " +  (Long_t)  temp_val1->GetEntries()+ " temp_res2: " +  (Long_t)  temp_res2->GetEntries()+ " temp_val2: " +  (Long_t)  temp_val2->GetEntries());
	temp_res1 = dummy_histo;
	temp_val1 = dummy_histo;
	temp_res2 = dummy_histo;
	temp_val2 = dummy_histo;
      }

      TF1 *fit_const = new TF1("fit_const","[0]",temp_val1->GetXaxis()->GetXmin(),temp_val1->GetXaxis()->GetXmax()); //was used before...
      fit_const->SetParameters(1,1);
      fit_const->SetParName(0,"const");

      TF1 *fit_loglin = new TF1("fit_loglin","[0]+[1]*TMath::Log(x)",temp_val1->GetXaxis()->GetXmin(),temp_val1->GetXaxis()->GetXmax()); //was used before...
      fit_loglin->SetParameters(1,1);
      fit_loglin->SetParName(0,"const");
      fit_loglin->SetParName(1,"slope");

      fit_const->SetLineColor(temp_res1->GetLineColor());
      fit_const->SetLineStyle(temp_res1->GetLineStyle());
      temp_res1->Fit("fit_const","","same");
      temp_res1->Fit("fit_loglin","+","same");

      fit_const->SetLineColor(temp_res2->GetLineColor());
      fit_const->SetLineStyle(temp_res2->GetLineStyle());
      temp_res2->Fit("fit_const","","same");

      fit_const->SetLineColor(temp_val1->GetLineColor());
      fit_const->SetLineStyle(temp_val1->GetLineStyle());
      temp_val1->Fit("fit_const","","same");

      fit_const->SetLineColor(temp_val2->GetLineColor());
      fit_const->SetLineStyle(temp_val2->GetLineStyle());
      temp_val2->Fit("fit_const","","same");


      temp_val1->SetTitle(Abseta_bins_labels_[Abseta_i].first +" < |#eta| < " + Abseta_bins_labels_[Abseta_i].second);
      temp_val2->SetTitle(Abseta_bins_labels_[Abseta_i].first +" < |#eta| < " + Abseta_bins_labels_[Abseta_i].second);
      temp_res1->SetTitle(Abseta_bins_labels_[Abseta_i].first +" < |#eta| < " + Abseta_bins_labels_[Abseta_i].second);
      temp_res2->SetTitle(Abseta_bins_labels_[Abseta_i].first +" < |#eta| < " + Abseta_bins_labels_[Abseta_i].second);

      all_Abseta_ratio_val1_.push_back(temp_val1);
      all_Abseta_ratio_val2_.push_back(temp_val2);
      all_Abseta_ratio_res1_.push_back(temp_res1);
      all_Abseta_ratio_res2_.push_back(temp_res2);

      Int_t bin_low = temp_res1->FindFirstBinAbove(0.1) ;
      Int_t bin_hig = temp_res1->FindLastBinAbove(0.1) ;
      cout << "lowest bin " << bin_low  << " highest bin " << bin_hig << endl; 
      cout << "lowest bin xlow " << temp_res1->GetXaxis()->GetBinLowEdge(bin_low)  << " highest bin xhig" << temp_res1->GetXaxis()->GetBinUpEdge(bin_hig) << endl; 
      all_Abseta_res1_ptreach_.push_back(make_pair( temp_res1->GetXaxis()->GetBinLowEdge(bin_low),temp_res1->GetXaxis()->GetBinUpEdge(bin_hig)));
    }

      cout << "test4" << endl;
      cout << "test4" << endl;
      cout << "test4" << endl;
    std::vector <std::vector<TH1D* > > overlay_Abseta_valres1_;
    overlay_Abseta_valres1_.push_back(all_Abseta_ratio_val1_);
    for(unsigned int style_i=0;style_i<all_Abseta_ratio_val1_.size();style_i++)TH1style_plus_fit(all_Abseta_ratio_val1_[style_i],line_styles_, colours_, markers_, 0,0,0,"R(MC)/R(data)","fit_const");
    overlay_Abseta_valres1_.push_back(all_Abseta_ratio_res1_);
    for(unsigned int style_i=0;style_i<all_Abseta_ratio_res1_.size();style_i++)TH1style_plus_fit(all_Abseta_ratio_res1_[style_i],line_styles_, colours_, markers_, 1,1,1,"R(MC)/R(data)","fit_const");
    
    TLegend *leg_Abseta_valres1;
    leg_Abseta_valres1 = new TLegend(0.25,0.70,0.55,0.85);
    leg_Abseta_valres1->SetFillColor(kWhite);
      //   leg->SetHeader("Legende");
    leg_Abseta_valres1->AddEntry(all_Abseta_ratio_val1_.back(),"r("+generatorone+")/r(L2L3res)","lep");
    leg_Abseta_valres1->AddEntry(all_Abseta_ratio_res1_.back(),"r("+generatorone+")/r(L2L3)","lep");
    
    if(export_all_plots)draw_Overlay_TH1D_save_PS(leg_Abseta_valres1,overlay_Abseta_valres1_, "Overlay_Abseta_"+generatorone+"_"+generatortwo+"_res1_and_val1_vs_pt_with_TJF_"+ptthreecuts[cut_i], "nice", "legend", "x1_y0_z0", 0,-1, 0.85, 1.35, 0,-1, image_ext, fine_coarse+"_"+ generatorone+"_"+generatortwo+"_"+jet_type+"_PS",1.);

    std::vector <std::vector<TH1D* > > overlay_Abseta_valres2_;
    overlay_Abseta_valres2_.push_back(all_Abseta_ratio_val2_);
    for(unsigned int style_i=0;style_i<all_Abseta_ratio_val2_.size();style_i++)TH1style_plus_fit(all_Abseta_ratio_val2_[style_i],line_styles_, colours_, markers_, 0,0,0,"R(MC)/R(data)","fit_const");
    overlay_Abseta_valres2_.push_back(all_Abseta_ratio_res2_);
    for(unsigned int style_i=0;style_i<all_Abseta_ratio_res2_.size();style_i++)TH1style_plus_fit(all_Abseta_ratio_res2_[style_i],line_styles_, colours_, markers_, 1,1,1,"R(MC)/R(data)","fit_const");
    
    TLegend *leg_Abseta_valres2;
    leg_Abseta_valres2 = new TLegend(0.25,0.70,0.55,0.85);
    leg_Abseta_valres2->SetFillColor(kWhite);
      //   leg->SetHeader("Legende");
    leg_Abseta_valres2->AddEntry(all_Abseta_ratio_val2_.back(),"r("+generatortwo+")/r(L2L3res)","lep");
    leg_Abseta_valres2->AddEntry(all_Abseta_ratio_res2_.back(),"r("+generatortwo+")/r(L2L3)","lep");
    
    // save memory
    //    if(export_all_plots)draw_Overlay_TH1D_save_PS(leg_Abseta_valres2,overlay_Abseta_valres2_, "Overlay_Abseta_"+generatorone+"_"+generatortwo+"_res2_and_val2_vs_pt_with_TJF_"+ptthreecuts[cut_i], "nice", "legend", "x1_y0_z0", 0,-1, 0.85, 1.35, 0,-1, image_ext, fine_coarse+"_"+ generatorone+"_"+generatortwo+"_"+jet_type+"_PS",1.);


    std::vector <std::vector<TH1D* > > overlay_Abseta_res1res2_;
    overlay_Abseta_res1res2_.push_back(all_Abseta_ratio_res1_);
    for(unsigned int style_i=0;style_i<all_Abseta_ratio_res1_.size();style_i++)TH1style_plus_fit(all_Abseta_ratio_res1_[style_i],line_styles_, colours_, markers_, 0,0,0,"R(MC)/R(data)","fit_const");
    overlay_Abseta_res1res2_.push_back(all_Abseta_ratio_res2_);
    for(unsigned int style_i=0;style_i<all_Abseta_ratio_res2_.size();style_i++)TH1style_plus_fit(all_Abseta_ratio_res2_[style_i],line_styles_, colours_, markers_, 1,1,1,"R(MC)/R(data)","fit_const");
    
    TLegend *leg_Abseta_res1res2;
    leg_Abseta_res1res2 = new TLegend(0.25,0.70,0.55,0.85);
    leg_Abseta_res1res2->SetFillColor(kWhite);
      //   leg->SetHeader("Legende");
    leg_Abseta_res1res2->AddEntry(all_Abseta_ratio_res1_.back(),"r("+generatorone+")/r(L2L3)","lep");
    leg_Abseta_res1res2->AddEntry(all_Abseta_ratio_res2_.back(),"r("+generatortwo+")/r(L2L3)","lep");
    
    // save memory
    //    if(export_all_plots)draw_Overlay_TH1D_save_PS(leg_Abseta_res1res2,overlay_Abseta_res1res2_, "Overlay_Abseta_"+generatorone+"_"+generatortwo+"_res1_and_res2_vs_pt_with_TJF_"+ptthreecuts[cut_i], "nice", "legend", "x1_y0_z0", 0,-1, 0.85, 1.35, 0,-1, image_ext, fine_coarse+"_"+ generatorone+"_"+generatortwo+"_"+jet_type+"_PS",1.);


    
    TLegend *leg_Abseta_res1;
    leg_Abseta_res1 = new TLegend(0.25,0.70,0.55,0.85);
    leg_Abseta_res1->SetFillColor(kWhite);
      //   leg->SetHeader("Legende");
    leg_Abseta_res1->AddEntry(all_Abseta_ratio_res1_.back(),"r("+generatorone+")/r(L2L3)","lep");
    std::vector <std::vector<TH1D* > > overlay_Abseta_res1_;
    overlay_Abseta_res1_.push_back(all_Abseta_ratio_res1_);
    if(export_all_plots)draw_Overlay_TH1D_save_PS(leg_Abseta_res1,overlay_Abseta_res1_, "Overlay_Abseta_"+generatorone+"_"+generatortwo+"_res1_vs_pt_with_TJF_"+ptthreecuts[cut_i], "tdr", "legend", "x1_y0_z0", 0,-1, 0.85, 1.35, 0,-1, image_ext, fine_coarse+"_"+ generatorone+"_"+generatortwo+"_"+jet_type+"_PS",1.,true);

    //     if(export_all_plots)draw_TGraphErrors_save_PS(fine_coarse+"_"+ generatorone+"_"+generatortwo+"_"+jet_type+"_PS_kFSR_extrapols",image_ext,extrapol_inclusive_Abseta_res1_,"testen_Abseta_res1","tdr","P","no","x0",0.,0.45,0.95,1.05,"#alpha","Normalized residual correction");//"#alpha","#frac{<R(MC)/R(data)>_{#alpha}}{<R(MC)/R(data)>_{#alpha=0.2}}");


 all_ptthree_all_Abseta_MC1_L2L3_      .push_back(all_Abseta_MC1_L2L3_   );	        	
 all_ptthree_all_Abseta_MC2_L2L3_      .push_back(all_Abseta_MC2_L2L3_   );	        	
 all_ptthree_all_Abseta_D_L2L3_        .push_back(all_Abseta_D_L2L3_ 	   );   	
 all_ptthree_all_Abseta_D_L2L3res_     .push_back(all_Abseta_D_L2L3res_  );	         	
 all_ptthree_all_Abseta_ratio_val1_    .push_back(all_Abseta_ratio_val1_ );	          	
 all_ptthree_all_Abseta_ratio_val2_    .push_back(all_Abseta_ratio_val2_ );	          	
 all_ptthree_all_Abseta_ratio_res1_    .push_back(all_Abseta_ratio_res1_ );	          	
 all_ptthree_all_Abseta_ratio_res2_    .push_back(all_Abseta_ratio_res2_ );

 all_ptthree_all_Abseta_res1_ptreach_  .push_back(all_Abseta_res1_ptreach_);


	          	

  }



    std::vector <TH1D*> NPV_all_eta_MC1_L2L3_;	 
    std::vector <TH1D*> NPV_all_eta_MC2_L2L3_;	 
    std::vector <TH1D*> NPV_all_eta_D_L2L3_;	 
    std::vector <TH1D*> NPV_all_eta_D_L2L3res_;	 
    std::vector <TH1D*> NPV_all_eta_ratio_val1_;	 
    std::vector <TH1D*> NPV_all_eta_ratio_val2_;	 
    std::vector <TH1D*> NPV_all_eta_ratio_res1_;	 
    std::vector <TH1D*> NPV_all_eta_ratio_res2_;	 

    std::vector <TH1D*> NPV_all_Abseta_MC1_L2L3_;	 
    std::vector <TH1D*> NPV_all_Abseta_MC2_L2L3_;	 
    std::vector <TH1D*> NPV_all_Abseta_D_L2L3_;	 
    std::vector <TH1D*> NPV_all_Abseta_D_L2L3res_;	 
    std::vector <TH1D*> NPV_all_Abseta_ratio_val1_;	 
    std::vector <TH1D*> NPV_all_Abseta_ratio_val2_;	 
    std::vector <TH1D*> NPV_all_Abseta_ratio_res1_;	 
    std::vector <TH1D*> NPV_all_Abseta_ratio_res2_;	 

    for(Int_t eta_i=0;eta_i<no_eta_bins;eta_i++){


      //      MyClass *obj;
      //      directory->GetObject("some object",obj);
      //      if (obj) { ... the object exist and inherits from MyClass ... }
      cout << ((TString) "AsymmetryVsNPV")+("/AsymmetryVsNPV")+"_AsymmetryVsVtxN_MC_L2L3_Eta"+(Long_t)eta_i+ratio_of_mean_or_GM << endl;
      TH1D* hMC1;
      fone->GetObject(((TString) "AsymmetryVsNPV")+("/AsymmetryVsNPV")+"_AsymmetryVsVtxN_MC_L2L3_Eta"+(Long_t)eta_i+ratio_of_mean_or_GM,hMC1);
      if(!hMC1)cout << "FEHLER - NICHT IMPORTIERT!" << endl;
      TH1D* hMC2;
      ftwo->GetObject(((TString) "AsymmetryVsNPV")+("/AsymmetryVsNPV")+"_AsymmetryVsVtxN_MC_L2L3_Eta"+(Long_t)eta_i+ratio_of_mean_or_GM,hMC2);
      if(!hMC2)cout << "FEHLER - NICHT IMPORTIERT!" << endl;
      TH1D* hL2L3;
      fone->GetObject(((TString) "AsymmetryVsNPV")+("/AsymmetryVsNPV")+"_AsymmetryVsVtxN_data_L2L3_Eta"+(Long_t)eta_i+ratio_of_mean_or_GM,hL2L3);
      if(!hL2L3)cout << "FEHLER - NICHT IMPORTIERT!" << endl;
      TH1D* hL2L3res;
      fone->GetObject(((TString) "AsymmetryVsNPV")+("/AsymmetryVsNPV")+"_AsymmetryVsVtxN_data_L2L3res_Eta"+(Long_t)eta_i+ratio_of_mean_or_GM,hL2L3res);
      if(!hL2L3res)cout << "FEHLER - NICHT IMPORTIERT!" << endl;


      
      if ( hMC1->GetEntries()==0||hMC2->GetEntries()==0||hL2L3->GetEntries()==0||hL2L3res->GetEntries()==0){
	std::cout << "DAS IST NICHT GUT!!! FEHLER!!!" << std::endl;
	//	hMC1->Dump();
      	std::cout <<  "AsymmetryVsNPV/AsymmetryVsNPV_AsymmetryVsVtxN_MC_L2L3_Eta"+(Long_t)eta_i+(ratio_of_mean_or_GM) << std::endl;


	faulty_values_dummy_histos.push_back("hMC-hL2L3: " +( (TString)"AsymmetryVsNPV")+("/AsymmetryVsNPV")+"_AsymmetryVsVtxN_MC_L2L3_Eta"+(Long_t)eta_i+ratio_of_mean_or_GM+  " hMC1: " +  (Long_t) hMC1->GetEntries() + " hMC2: " + (Long_t) hMC2->GetEntries()+ " hL2L3: " + (Long_t) hL2L3->GetEntries()+ " hL2L3res: " + (Long_t) hL2L3res->GetEntries());	
	hMC1=dummy_histo;
	hMC2=dummy_histo;
	hL2L3=dummy_histo;
	hL2L3res=dummy_histo;

	no_entries_went_wrong=true;
      }

//      if(NPV_all_eta_MC1_L2L3_.back()->GetEntries()==0)std::cout << "this cannot work..." << std::endl;
//      if(NPV_all_eta_MC2_L2L3_.back()->GetEntries()==0)std::cout << "this cannot work..." << std::endl;
//      if(NPV_all_eta_D_L2L3_.back()->GetEntries()==0)std::cout << "this cannot work..." << std::endl;
//      if(NPV_all_eta_D_L2L3res_.back()->GetEntries()==0)std::cout << "this cannot work..." << std::endl;

      NPV_all_eta_MC1_L2L3_.push_back(hMC1);
      NPV_all_eta_MC2_L2L3_.push_back(hMC2);
      NPV_all_eta_D_L2L3_.push_back(hL2L3);
      NPV_all_eta_D_L2L3res_.push_back(hL2L3res);

      //      cout << "test" << endl;
      NPV_all_eta_MC1_L2L3_ .back()->GetYaxis()->SetRangeUser(0.5,1.5);
      NPV_all_eta_MC2_L2L3_ .back()->GetYaxis()->SetRangeUser(0.5,1.5);
      NPV_all_eta_D_L2L3_   .back()->GetYaxis()->SetRangeUser(0.5,1.5);
      NPV_all_eta_D_L2L3res_.back()->GetYaxis()->SetRangeUser(0.5,1.5);
      //      cout << "test2" << endl;


      TH1D* temp_res1 = (TH1D*) NPV_all_eta_D_L2L3res_.back()->Clone();
      TH1D* temp_val1 = (TH1D*) NPV_all_eta_D_L2L3_.back()->Clone();
      TH1D* temp_res2 = (TH1D*) NPV_all_eta_D_L2L3res_.back()->Clone();
      TH1D* temp_val2 = (TH1D*) NPV_all_eta_D_L2L3_.back()->Clone();

      temp_res1->Sumw2();
      temp_val1->Sumw2();
      temp_res2->Sumw2();
      temp_val2->Sumw2();
      temp_val1->Divide(NPV_all_eta_MC1_L2L3_.back(),NPV_all_eta_D_L2L3res_.back());
      temp_res1->Divide(NPV_all_eta_MC1_L2L3_.back(),NPV_all_eta_D_L2L3_.back());
      temp_val2->Divide(NPV_all_eta_MC2_L2L3_.back(),NPV_all_eta_D_L2L3res_.back());
      temp_res2->Divide(NPV_all_eta_MC2_L2L3_.back(),NPV_all_eta_D_L2L3_.back());

      TH1style(temp_val1,line_styles_, colours_, markers_, 0,0,0,"R(MC)/R(data)");
      TH1style(temp_val2,line_styles_, colours_, markers_, 1,1,1,"R(MC)/R(data)");
      TH1style(temp_res1,line_styles_, colours_, markers_, 2,2,2,"R(MC)/R(data)");
      TH1style(temp_res2,line_styles_, colours_, markers_, 3,3,3,"R(MC)/R(data)");

      TF1 *fit_const = new TF1("fit_const","[0]",temp_val1->GetXaxis()->GetXmin(),temp_val1->GetXaxis()->GetXmax()); //was used before...
      fit_const->SetParameters(1,1);
      fit_const->SetParName(0,"const");


      if (  temp_res1->GetEntries()==0|| temp_res2->GetEntries()==0|| temp_val1->GetEntries()==0||temp_val2->GetEntries()==0){
	std::cout << "DAS IST NICHT GUT!!! FEHLER!!!" << std::endl;
	//	hMC1->Dump();
	faulty_values_dummy_histos.push_back("divided (ratios): " +( (TString)"AsymmetryVsNPV")+("/AsymmetryVsNPV")+"_AsymmetryVsVtxN_MC_L2L3_Eta"+(Long_t)eta_i+ratio_of_mean_or_GM+ " temp_res1: " +  (Long_t)  temp_res1->GetEntries()+ " temp_val1: " +  (Long_t)  temp_val1->GetEntries()+ " temp_res2: " +  (Long_t)  temp_res2->GetEntries()+ " temp_val2: " +  (Long_t)  temp_val2->GetEntries());
	temp_res1 = dummy_histo;
	temp_val1 = dummy_histo;
	temp_res2 = dummy_histo;
	temp_val2 = dummy_histo;
      }

      cout << hMC1->GetName() << endl;
      cout << "eta: " << eta_i << " and test " <<  eta_bins_labels_[eta_i].first << " " << eta_bins_labels_[eta_i].second << endl;
      fit_const->SetLineColor(temp_res1->GetLineColor());
      fit_const->SetLineStyle(temp_res1->GetLineStyle());
      temp_res1->Fit("fit_const","","same");

      cout << "eta: " << eta_i << " and test " <<  eta_bins_labels_[eta_i].first << " " << eta_bins_labels_[eta_i].second << endl;
      fit_const->SetLineColor(temp_res2->GetLineColor());
      fit_const->SetLineStyle(temp_res2->GetLineStyle());
      temp_res2->Fit("fit_const","","same");

      cout << "eta: " << eta_i << " and test " <<  eta_bins_labels_[eta_i].first << " " << eta_bins_labels_[eta_i].second << endl;
      fit_const->SetLineColor(temp_val1->GetLineColor());
      fit_const->SetLineStyle(temp_val1->GetLineStyle());
      temp_val1->Fit("fit_const","","same");

      cout << "eta: " << eta_i << " and test " <<  eta_bins_labels_[eta_i].first << " " << eta_bins_labels_[eta_i].second << endl;
      fit_const->SetLineColor(temp_val2->GetLineColor());
      fit_const->SetLineStyle(temp_val2->GetLineStyle());
      temp_val2->Fit("fit_const","","same");



      NPV_all_eta_ratio_val1_.push_back(temp_val1);
      NPV_all_eta_ratio_val2_.push_back(temp_val2);
      NPV_all_eta_ratio_res1_.push_back(temp_res1);
      NPV_all_eta_ratio_res2_.push_back(temp_res2);

    }


    for(Int_t Abseta_i=0;Abseta_i<no_Abseta_bins;Abseta_i++){


      //      MyClass *obj;
      //      directory->GetObject("some object",obj);
      //      if (obj) { ... the object exist and inherits from MyClass ... }

      TH1D* hMC1;
      cout <<((TString) "AbsAsymmetryVsNPV")+("/AbsAsymmetryVsNPV")+"_AsymmetryVsVtxN_MC_L2L3_AbsEta"+(Long_t)Abseta_i+ratio_of_mean_or_GM << endl;
      fone->GetObject(((TString) "AbsAsymmetryVsNPV")+("/AbsAsymmetryVsNPV")+"_AsymmetryVsVtxN_MC_L2L3_AbsEta"+(Long_t)Abseta_i+ratio_of_mean_or_GM,hMC1);
      if(!hMC1)cout << "FEHLER - NICHT IMPORTIERT!" << endl;
      TH1D* hMC2;
      ftwo->GetObject(((TString) "AbsAsymmetryVsNPV")+("/AbsAsymmetryVsNPV")+"_AsymmetryVsVtxN_MC_L2L3_AbsEta"+(Long_t)Abseta_i+ratio_of_mean_or_GM,hMC2);
      if(!hMC2)cout << "FEHLER - NICHT IMPORTIERT!" << endl;
      TH1D* hL2L3;
      fone->GetObject(((TString) "AbsAsymmetryVsNPV")+("/AbsAsymmetryVsNPV")+"_AsymmetryVsVtxN_data_L2L3_AbsEta"+(Long_t)Abseta_i+ratio_of_mean_or_GM,hL2L3);
      if(!hL2L3)cout << "FEHLER - NICHT IMPORTIERT!" << endl;
      TH1D* hL2L3res;
      fone->GetObject(((TString) "AbsAsymmetryVsNPV")+("/AbsAsymmetryVsNPV")+"_AsymmetryVsVtxN_data_L2L3res_AbsEta"+(Long_t)Abseta_i+ratio_of_mean_or_GM,hL2L3res);
      if(!hL2L3res)cout << "FEHLER - NICHT IMPORTIERT!" << endl;


      
      if ( hMC1->GetEntries()==0||hMC2->GetEntries()==0||hL2L3->GetEntries()==0||hL2L3res->GetEntries()==0){
	std::cout << "DAS IST NICHT GUT!!! FEHLER!!!" << std::endl;
	//	hMC1->Dump();
      	std::cout << ( "AbsAsymmetryVsNPV")+((TString)"/AbsAsymmetryVsNPV")+"_AsymmetryVsVtxN_MC_L2L3_AbsEta"+(Long_t)Abseta_i+ratio_of_mean_or_GM << std::endl;


	faulty_values_dummy_histos.push_back("hMC-hL2L3: " +((TString) "AbsAsymmetryVsNPV")+("/AbsAsymmetryVsNPV")+"_AsymmetryVsVtxN_MC_L2L3_AbsEta"+(Long_t)Abseta_i+ratio_of_mean_or_GM+  " hMC1: " +  (Long_t) hMC1->GetEntries() + " hMC2: " + (Long_t) hMC2->GetEntries()+ " hL2L3: " + (Long_t) hL2L3->GetEntries()+ " hL2L3res: " + (Long_t) hL2L3res->GetEntries());	
	hMC1=dummy_histo;
	hMC2=dummy_histo;
	hL2L3=dummy_histo;
	hL2L3res=dummy_histo;

	no_entries_went_wrong=true;
      }

//      if(NPV_all_Abseta_MC1_L2L3_.back()->GetEntries()==0)std::cout << "this cannot work..." << std::endl;
//      if(NPV_all_Abseta_MC2_L2L3_.back()->GetEntries()==0)std::cout << "this cannot work..." << std::endl;
//      if(NPV_all_Abseta_D_L2L3_.back()->GetEntries()==0)std::cout << "this cannot work..." << std::endl;
//      if(NPV_all_Abseta_D_L2L3res_.back()->GetEntries()==0)std::cout << "this cannot work..." << std::endl;

      NPV_all_Abseta_MC1_L2L3_.push_back(hMC1);
      NPV_all_Abseta_MC2_L2L3_.push_back(hMC2);
      NPV_all_Abseta_D_L2L3_.push_back(hL2L3);
      NPV_all_Abseta_D_L2L3res_.push_back(hL2L3res);

      //      cout << "test" << endl;
      NPV_all_Abseta_MC1_L2L3_ .back()->GetYaxis()->SetRangeUser(0.5,1.5);
      NPV_all_Abseta_MC2_L2L3_ .back()->GetYaxis()->SetRangeUser(0.5,1.5);
      NPV_all_Abseta_D_L2L3_   .back()->GetYaxis()->SetRangeUser(0.5,1.5);
      NPV_all_Abseta_D_L2L3res_.back()->GetYaxis()->SetRangeUser(0.5,1.5);
      //      cout << "test2" << endl;


      TH1D* temp_res1 = (TH1D*) NPV_all_Abseta_D_L2L3res_.back()->Clone();
      TH1D* temp_val1 = (TH1D*) NPV_all_Abseta_D_L2L3_.back()->Clone();
      TH1D* temp_res2 = (TH1D*) NPV_all_Abseta_D_L2L3res_.back()->Clone();
      TH1D* temp_val2 = (TH1D*) NPV_all_Abseta_D_L2L3_.back()->Clone();

      temp_res1->Sumw2();
      temp_val1->Sumw2();
      temp_res2->Sumw2();
      temp_val2->Sumw2();
      temp_val1->Divide(NPV_all_Abseta_MC1_L2L3_.back(),NPV_all_Abseta_D_L2L3res_.back());
      temp_res1->Divide(NPV_all_Abseta_MC1_L2L3_.back(),NPV_all_Abseta_D_L2L3_.back());
      temp_val2->Divide(NPV_all_Abseta_MC2_L2L3_.back(),NPV_all_Abseta_D_L2L3res_.back());
      temp_res2->Divide(NPV_all_Abseta_MC2_L2L3_.back(),NPV_all_Abseta_D_L2L3_.back());

      TH1style(temp_val1,line_styles_, colours_, markers_, 0,0,0,"R(MC)/R(data)");
      TH1style(temp_val2,line_styles_, colours_, markers_, 1,1,1,"R(MC)/R(data)");
      TH1style(temp_res1,line_styles_, colours_, markers_, 2,2,2,"R(MC)/R(data)");
      TH1style(temp_res2,line_styles_, colours_, markers_, 3,3,3,"R(MC)/R(data)");

      TF1 *fit_const = new TF1("fit_const","[0]",temp_val1->GetXaxis()->GetXmin(),temp_val1->GetXaxis()->GetXmax()); //was used before...
      fit_const->SetParameters(1,1);
      fit_const->SetParName(0,"const");


      if (  temp_res1->GetEntries()==0|| temp_res2->GetEntries()==0|| temp_val1->GetEntries()==0||temp_val2->GetEntries()==0){
	std::cout << "DAS IST NICHT GUT!!! FEHLER!!!" << std::endl;
	//	hMC1->Dump();
	faulty_values_dummy_histos.push_back("divided (ratios): " +( (TString)"AbsAsymmetryVsNPV")+("/AbsAsymmetryVsNPV")+"_AsymmetryVsVtxN_MC_L2L3_Abseta"+(Long_t)Abseta_i+ratio_of_mean_or_GM+ " temp_res1: " +  (Long_t)  temp_res1->GetEntries()+ " temp_val1: " +  (Long_t)  temp_val1->GetEntries()+ " temp_res2: " +  (Long_t)  temp_res2->GetEntries()+ " temp_val2: " +  (Long_t)  temp_val2->GetEntries());
	temp_res1 = dummy_histo;
	temp_val1 = dummy_histo;
	temp_res2 = dummy_histo;
	temp_val2 = dummy_histo;
      }

      cout << hMC1->GetName() << endl;
      cout << "Abseta: " << Abseta_i << " and test " <<  Abseta_bins_labels_[Abseta_i].first << " " << Abseta_bins_labels_[Abseta_i].second << endl;
      fit_const->SetLineColor(temp_res1->GetLineColor());
      fit_const->SetLineStyle(temp_res1->GetLineStyle());
      temp_res1->Fit("fit_const","","same");

      cout << "Abseta: " << Abseta_i << " and test " <<  Abseta_bins_labels_[Abseta_i].first << " " << Abseta_bins_labels_[Abseta_i].second << endl;
      fit_const->SetLineColor(temp_res2->GetLineColor());
      fit_const->SetLineStyle(temp_res2->GetLineStyle());
      temp_res2->Fit("fit_const","","same");

      cout << "Abseta: " << Abseta_i << " and test " <<  Abseta_bins_labels_[Abseta_i].first << " " << Abseta_bins_labels_[Abseta_i].second << endl;
      fit_const->SetLineColor(temp_val1->GetLineColor());
      fit_const->SetLineStyle(temp_val1->GetLineStyle());
      temp_val1->Fit("fit_const","","same");

      cout << "Abseta: " << Abseta_i << " and test " <<  Abseta_bins_labels_[Abseta_i].first << " " << Abseta_bins_labels_[Abseta_i].second << endl;
      fit_const->SetLineColor(temp_val2->GetLineColor());
      fit_const->SetLineStyle(temp_val2->GetLineStyle());
      temp_val2->Fit("fit_const","","same");



      NPV_all_Abseta_ratio_val1_.push_back(temp_val1);
      NPV_all_Abseta_ratio_val2_.push_back(temp_val2);
      NPV_all_Abseta_ratio_res1_.push_back(temp_res1);
      NPV_all_Abseta_ratio_res2_.push_back(temp_res2);

    }



    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////Extrapolation section/////////////////////////////////////////////////////////////////////
    //////////////////////////////Extrapolation section/////////////////////////////////////////////////////////////////////
    //////////////////////////////Extrapolation section/////////////////////////////////////////////////////////////////////
    //////////////////////////////Extrapolation section/////////////////////////////////////////////////////////////////////
    //////////////////////////////Extrapolation section/////////////////////////////////////////////////////////////////////
    //////////////////////////////Extrapolation section/////////////////////////////////////////////////////////////////////
    //////////////////////////////Extrapolation section/////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  std::vector <TGraphErrors*> extrapol_inclusive_Abseta_val1_;
  std::vector <TGraphErrors*> extrapol_inclusive_Abseta_res1_;
  std::vector <TGraphErrors*> extrapol_inclusive_Abseta_val2_;
  std::vector <TGraphErrors*> extrapol_inclusive_Abseta_res2_;
    for(Int_t Abseta_i=0;Abseta_i<no_Abseta_bins;Abseta_i++){

      //      TH1D *hextrapol_val1 = new TH1D("hextrapol_val1","",ptthreecuts_Double_.size(),&ptthreecuts_Double_[0]);


      std::vector <Double_t> x_ptthree_;
      std::vector <Double_t> ex_ptthree_;
      std::vector <Double_t> y_mean_ratio_res1_;
      std::vector <Double_t> ey_mean_ratio_res1_;
      std::vector <Double_t> y_mean_ratio_val1_;
      std::vector <Double_t> ey_mean_ratio_val1_;

      std::vector <Double_t> y_mean_ratio_res2_;
      std::vector <Double_t> ey_mean_ratio_res2_;
      std::vector <Double_t> y_mean_ratio_val2_;
      std::vector <Double_t> ey_mean_ratio_val2_;

      Int_t cut_20=1;
      
      for(unsigned int cut_j=0;cut_j<ptthreecuts.size();cut_j++){
	if(ptthreecuts[cut_j]=="20"){
	  cut_20=cut_j;}
      }

      for(unsigned int cut_i=0;cut_i<ptthreecuts.size();cut_i++){
	x_ptthree_.push_back(ptthreecuts_Double_.at(cut_i));
	ex_ptthree_.push_back(0.00);


	Double_t y_val1=all_ptthree_all_Abseta_ratio_val1_  .at(cut_i).at(Abseta_i)->GetFunction("fit_const")->GetParameter(0);
	Double_t ey_val1=all_ptthree_all_Abseta_ratio_val1_  .at(cut_i).at(Abseta_i)->GetFunction("fit_const")->GetParError(0);
	Double_t y_val1_20=all_ptthree_all_Abseta_ratio_val1_  .at(cut_20).at(Abseta_i)->GetFunction("fit_const")->GetParameter(0);
	Double_t ey_val1_20=all_ptthree_all_Abseta_ratio_val1_  .at(cut_20).at(Abseta_i)->GetFunction("fit_const")->GetParError(0); 
	Double_t y_val2=all_ptthree_all_Abseta_ratio_val2_  .at(cut_i).at(Abseta_i)->GetFunction("fit_const")->GetParameter(0);
	Double_t ey_val2=all_ptthree_all_Abseta_ratio_val2_  .at(cut_i).at(Abseta_i)->GetFunction("fit_const")->GetParError(0);
	Double_t y_val2_20=all_ptthree_all_Abseta_ratio_val2_  .at(cut_20).at(Abseta_i)->GetFunction("fit_const")->GetParameter(0);
	Double_t ey_val2_20=all_ptthree_all_Abseta_ratio_val2_  .at(cut_20).at(Abseta_i)->GetFunction("fit_const")->GetParError(0); 
	Double_t y_res1=all_ptthree_all_Abseta_ratio_res1_  .at(cut_i).at(Abseta_i)->GetFunction("fit_const")->GetParameter(0);
	Double_t ey_res1=all_ptthree_all_Abseta_ratio_res1_  .at(cut_i).at(Abseta_i)->GetFunction("fit_const")->GetParError(0);
	Double_t y_res1_20=all_ptthree_all_Abseta_ratio_res1_  .at(cut_20).at(Abseta_i)->GetFunction("fit_const")->GetParameter(0);
	Double_t ey_res1_20=all_ptthree_all_Abseta_ratio_res1_  .at(cut_20).at(Abseta_i)->GetFunction("fit_const")->GetParError(0); 
	Double_t y_res2=all_ptthree_all_Abseta_ratio_res2_  .at(cut_i).at(Abseta_i)->GetFunction("fit_const")->GetParameter(0);
	Double_t ey_res2=all_ptthree_all_Abseta_ratio_res2_  .at(cut_i).at(Abseta_i)->GetFunction("fit_const")->GetParError(0);
	Double_t y_res2_20=all_ptthree_all_Abseta_ratio_res2_  .at(cut_20).at(Abseta_i)->GetFunction("fit_const")->GetParameter(0);
	Double_t ey_res2_20=all_ptthree_all_Abseta_ratio_res2_  .at(cut_20).at(Abseta_i)->GetFunction("fit_const")->GetParError(0); 

	y_mean_ratio_val1_.push_back(y_val1/y_val1_20);
	y_mean_ratio_val2_.push_back(y_val2/y_val2_20);
	y_mean_ratio_res1_.push_back(y_res1/y_res1_20);
	y_mean_ratio_res2_.push_back(y_res2/y_res2_20);


	//classic method: assume uncorrelated errors:
	ey_mean_ratio_val1_.push_back(TMath::Sqrt((ey_val1*ey_val1/(y_val1_20*y_val1_20))+((y_val1*y_val1*ey_val1_20*ey_val1_20)/TMath::Power(y_val1_20,4))));
	ey_mean_ratio_val2_.push_back(TMath::Sqrt((ey_val2*ey_val2/(y_val2_20*y_val2_20))+((y_val2*y_val2*ey_val2_20*ey_val2_20)/TMath::Power(y_val2_20,4))));

	ey_mean_ratio_res1_.push_back(TMath::Sqrt((ey_res1*ey_res1/(y_res1_20*y_res1_20))+((y_res1*y_res1*ey_res1_20*ey_res1_20)/TMath::Power(y_res1_20,4))));
	ey_mean_ratio_res2_.push_back(TMath::Sqrt((ey_res2*ey_res2/(y_res2_20*y_res2_20))+((y_res2*y_res2*ey_res2_20*ey_res2_20)/TMath::Power(y_res2_20,4))));



	cout << "DAS KOMMT RAUS: " << y_mean_ratio_res1_.back() << " und Fehler: " << ey_mean_ratio_res1_.back() << endl;
	cout << "DETAILS.. y_res1: " << y_res1 << " ey_res1: " << ey_res1 << " y_res1_20: " << y_res1_20 << " ey_res1_20: " << y_res1_20 << endl;
	cout << "In der Wurzel:  " << 1/(y_res1_20*y_res1_20)*ey_res1*ey_res1+ey_res1_20*ey_res1_20*(TMath::Power(y_res1-y_res1_20,2)/TMath::Power(y_res1_20,4)-1/TMath::Power(y_res1_20,2)) << " Teil 1: " << TMath::Power(y_res1-y_res1_20,2)/TMath::Power(y_res1_20,4) << " Teil2: "<< 1/TMath::Power(y_res1_20,2) << endl;

	
      }


      Int_t n = y_mean_ratio_val1_.size();



      TGraphErrors *gr_val1 = new TGraphErrors(n,&x_ptthree_[0],&y_mean_ratio_val1_[0],&ex_ptthree_[0],&ey_mean_ratio_val1_[0]);
      gr_val1->SetTitle(Abseta_bins_labels_[Abseta_i].first +" < |#eta| < " + Abseta_bins_labels_[Abseta_i].second);
      gr_val1->SetName((TString)"kFSR_extrapol_val1_Abseta_"+(Long_t)Abseta_i);
      //      TGraphErrorsstyle(gr_val1,line_styles_, colours_, markers_, 0,0,0,"bla");


      TGraphErrors *gr_res1 = new TGraphErrors(n,&x_ptthree_[0],&y_mean_ratio_res1_[0],&ex_ptthree_[0],&ey_mean_ratio_res1_[0]);
      gr_res1->Draw();
      gr_res1->SetTitle(Abseta_bins_labels_[Abseta_i].first +" < |#eta| < " + Abseta_bins_labels_[Abseta_i].second);
      gr_res1->SetName((TString)"kFSR_extrapol_res1_Abseta_"+(Long_t)Abseta_i);
      //           TGraphErrorsstyle(gr_res1,line_styles_, colours_, markers_, 0,0,0,"bla");

      TGraphErrors *gr_val2 = new TGraphErrors(n,&x_ptthree_[0],&y_mean_ratio_val2_[0],&ex_ptthree_[0],&ey_mean_ratio_val2_[0]);
      gr_val2->SetTitle(Abseta_bins_labels_[Abseta_i].first +" < |#eta| < " + Abseta_bins_labels_[Abseta_i].second);
      gr_val2->SetName((TString)"kFSR_extrapol_val2_Abseta_"+(Long_t)Abseta_i);
      //      TGraphErrorsstyle(gr_val2,line_styles_, colours_, markers_, 0,0,0,"bla");


      TGraphErrors *gr_res2 = new TGraphErrors(n,&x_ptthree_[0],&y_mean_ratio_res2_[0],&ex_ptthree_[0],&ey_mean_ratio_res2_[0]);
      gr_res2->SetTitle(Abseta_bins_labels_[Abseta_i].first +" < |#eta| < " + Abseta_bins_labels_[Abseta_i].second);
      gr_res2->SetName((TString)"kFSR_extrapol_res2_Abseta_"+(Long_t)Abseta_i);
      //      TGraphErrorsstyle(gr_res2,line_styles_, colours_, markers_, 0,0,0,"bla");


     TF1 *lin_extrapol = new TF1("lin_extrapol","[0]+[1]*x",0,0.5); //was used before...
      lin_extrapol->SetParameters(1,-0.1);
      lin_extrapol->SetParName(0,"k_{FSR}");
      lin_extrapol->SetParName(1,"slope");

      //deactivate because JINST
      //      lin_extrapol->SetLineColor(gr_val1->GetLineColor());
      lin_extrapol->SetLineColor(2);
      gr_val1->Fit("lin_extrapol","","same",0,0.5);

      //      lin_extrapol->SetLineColor(gr_res1->GetLineColor());
      gr_res1->Fit("lin_extrapol","","same",0,0.5);

      extrapol_inclusive_Abseta_val1_.push_back(gr_val1);
      extrapol_inclusive_Abseta_res1_.push_back(gr_res1);

      //      lin_extrapol->SetLineColor(gr_val2->GetLineColor());
      gr_val2->Fit("lin_extrapol","","same",0,0.5);

      //      lin_extrapol->SetLineColor(gr_res2->GetLineColor());
      gr_res2->Fit("lin_extrapol","","same",0,0.5);

      extrapol_inclusive_Abseta_val2_.push_back(gr_val2);
      extrapol_inclusive_Abseta_res2_.push_back(gr_res2);

    }


    if(export_all_plots)draw_TGraphErrors_save_PS(fine_coarse+"_"+ generatorone+"_"+generatortwo+"_"+jet_type+"_PS_kFSR_extrapols",image_ext,extrapol_inclusive_Abseta_res1_,"testen_Abseta_res1","tdr","P","no","x0",0.,0.45,0.95,1.05,"cut on #alpha","Normalized residual correction",true);//"#alpha","#frac{<R(MC)/R(data)>_{#alpha}}{<R(MC)/R(data)>_{#alpha=0.2}}");
    // save memory
    //     if(export_all_plots)draw_TGraphErrors_save_PS(fine_coarse+"_"+ generatorone+"_"+generatortwo+"_"+jet_type+"_PS_kFSR_extrapols",image_ext,extrapol_inclusive_Abseta_res2_,"testen_Abseta_res2","tdr","P","no","x0",0.,0.45,0.95,1.05,"#alpha","Normalized residual correction");//"#alpha","#frac{<R(MC)/R(data)>_{#alpha}}{<R(MC)/R(data)>_{#alpha=0.2}}");
    if(export_all_plots)draw_TGraphErrors_save_PS(fine_coarse+"_"+ generatorone+"_"+generatortwo+"_"+jet_type+"_PS_kFSR_extrapols",image_ext,extrapol_inclusive_Abseta_val1_,"testen_Abseta_val1","tdr","P","no","x0",0.,0.45,0.95,1.05,"cut on #alpha","Normalized residual correction",true);//"#alpha","#frac{<R(MC)/R(data)>_{#alpha}}{<R(MC)/R(data)>_{#alpha=0.2}}");
    // save memory
    //     if(export_all_plots)draw_TGraphErrors_save_PS(fine_coarse+"_"+ generatorone+"_"+generatortwo+"_"+jet_type+"_PS_kFSR_extrapols",image_ext,extrapol_inclusive_Abseta_val2_,"testen_Abseta_val2","tdr","P","no","x0",0.,0.45,0.95,1.05,"#alpha","Normalized residual correction");//"#alpha","#frac{<R(MC)/R(data)>_{#alpha}}{<R(MC)/R(data)>_{#alpha=0.2}}");






    /////////////////////////
    ////////////////////////
    ////////////same with eta
    /////////////////////////
    /////////////////////////



  std::vector <TGraphErrors*> extrapol_inclusive_eta_val1_;
  std::vector <TGraphErrors*> extrapol_inclusive_eta_res1_;
  std::vector <TGraphErrors*> extrapol_inclusive_eta_val2_;
  std::vector <TGraphErrors*> extrapol_inclusive_eta_res2_;
    for(Int_t eta_i=0;eta_i<no_eta_bins;eta_i++){

      //      TH1D *hextrapol_val1 = new TH1D("hextrapol_val1","",ptthreecuts_Double_.size(),&ptthreecuts_Double_[0]);


      std::vector <Double_t> x_ptthree_;
      std::vector <Double_t> ex_ptthree_;
      std::vector <Double_t> y_mean_ratio_res1_;
      std::vector <Double_t> ey_mean_ratio_res1_;
      std::vector <Double_t> y_mean_ratio_val1_;
      std::vector <Double_t> ey_mean_ratio_val1_;

      std::vector <Double_t> y_mean_ratio_res2_;
      std::vector <Double_t> ey_mean_ratio_res2_;
      std::vector <Double_t> y_mean_ratio_val2_;
      std::vector <Double_t> ey_mean_ratio_val2_;

      Int_t cut_20=1;
      
      for(unsigned int cut_j=0;cut_j<ptthreecuts.size();cut_j++){
	if(ptthreecuts[cut_j]=="20"){
	  cut_20=cut_j;}
      }


      for(unsigned int cut_i=0;cut_i<ptthreecuts.size();cut_i++){
	x_ptthree_.push_back(ptthreecuts_Double_.at(cut_i));
	ex_ptthree_.push_back(0.00);



	Double_t y_val1=all_ptthree_all_eta_ratio_val1_  .at(cut_i).at(eta_i)->GetFunction("fit_const")->GetParameter(0);
	Double_t ey_val1=all_ptthree_all_eta_ratio_val1_  .at(cut_i).at(eta_i)->GetFunction("fit_const")->GetParError(0);
	Double_t y_val1_20=all_ptthree_all_eta_ratio_val1_  .at(cut_20).at(eta_i)->GetFunction("fit_const")->GetParameter(0);
	Double_t ey_val1_20=all_ptthree_all_eta_ratio_val1_  .at(cut_20).at(eta_i)->GetFunction("fit_const")->GetParError(0); 
	Double_t y_val2=all_ptthree_all_eta_ratio_val2_  .at(cut_i).at(eta_i)->GetFunction("fit_const")->GetParameter(0);
	Double_t ey_val2=all_ptthree_all_eta_ratio_val2_  .at(cut_i).at(eta_i)->GetFunction("fit_const")->GetParError(0);
	Double_t y_val2_20=all_ptthree_all_eta_ratio_val2_  .at(cut_20).at(eta_i)->GetFunction("fit_const")->GetParameter(0);
	Double_t ey_val2_20=all_ptthree_all_eta_ratio_val2_  .at(cut_20).at(eta_i)->GetFunction("fit_const")->GetParError(0); 
	Double_t y_res1=all_ptthree_all_eta_ratio_res1_  .at(cut_i).at(eta_i)->GetFunction("fit_const")->GetParameter(0);
	Double_t ey_res1=all_ptthree_all_eta_ratio_res1_  .at(cut_i).at(eta_i)->GetFunction("fit_const")->GetParError(0);
	Double_t y_res1_20=all_ptthree_all_eta_ratio_res1_  .at(cut_20).at(eta_i)->GetFunction("fit_const")->GetParameter(0);
	Double_t ey_res1_20=all_ptthree_all_eta_ratio_res1_  .at(cut_20).at(eta_i)->GetFunction("fit_const")->GetParError(0); 
	Double_t y_res2=all_ptthree_all_eta_ratio_res2_  .at(cut_i).at(eta_i)->GetFunction("fit_const")->GetParameter(0);
	Double_t ey_res2=all_ptthree_all_eta_ratio_res2_  .at(cut_i).at(eta_i)->GetFunction("fit_const")->GetParError(0);
	Double_t y_res2_20=all_ptthree_all_eta_ratio_res2_  .at(cut_20).at(eta_i)->GetFunction("fit_const")->GetParameter(0);
	Double_t ey_res2_20=all_ptthree_all_eta_ratio_res2_  .at(cut_20).at(eta_i)->GetFunction("fit_const")->GetParError(0); 

	y_mean_ratio_val1_.push_back(y_val1/y_val1_20);
	y_mean_ratio_val2_.push_back(y_val2/y_val2_20);
	y_mean_ratio_res1_.push_back(y_res1/y_res1_20);
	y_mean_ratio_res2_.push_back(y_res2/y_res2_20);


	//classic method: assume uncorrelated errors:
	ey_mean_ratio_val1_.push_back(TMath::Sqrt((ey_val1*ey_val1/(y_val1_20*y_val1_20))+((y_val1*y_val1*ey_val1_20*ey_val1_20)/TMath::Power(y_val1_20,4))));
	ey_mean_ratio_val2_.push_back(TMath::Sqrt((ey_val2*ey_val2/(y_val2_20*y_val2_20))+((y_val2*y_val2*ey_val2_20*ey_val2_20)/TMath::Power(y_val2_20,4))));

	ey_mean_ratio_res1_.push_back(TMath::Sqrt((ey_res1*ey_res1/(y_res1_20*y_res1_20))+((y_res1*y_res1*ey_res1_20*ey_res1_20)/TMath::Power(y_res1_20,4))));
	ey_mean_ratio_res2_.push_back(TMath::Sqrt((ey_res2*ey_res2/(y_res2_20*y_res2_20))+((y_res2*y_res2*ey_res2_20*ey_res2_20)/TMath::Power(y_res2_20,4))));


	
      }


      Int_t n = y_mean_ratio_val1_.size();
      //      cout << "wie oft...: " << n << endl;


      TGraphErrors *gr_val1 = new TGraphErrors(n,&x_ptthree_[0],&y_mean_ratio_val1_[0],&ex_ptthree_[0],&ey_mean_ratio_val1_[0]);
      gr_val1->SetTitle(eta_bins_labels_[eta_i].first +" < #eta < " + eta_bins_labels_[eta_i].second);
      gr_val1->SetName((TString)"kFSR_extrapol_val1_eta_"+(Long_t)eta_i);
      //      TGraphErrorsstyle(gr_val1,line_styles_, colours_, markers_, 0,0,0,"bla");


      TGraphErrors *gr_res1 = new TGraphErrors(n,&x_ptthree_[0],&y_mean_ratio_res1_[0],&ex_ptthree_[0],&ey_mean_ratio_res1_[0]);
      gr_res1->SetTitle(eta_bins_labels_[eta_i].first +" < #eta < " + eta_bins_labels_[eta_i].second);
      gr_res1->SetName((TString)"kFSR_extrapol_res1_eta_"+(Long_t)eta_i);
      //      TGraphErrorsstyle(gr_res1,line_styles_, colours_, markers_, 1,1,1,"bla");

      TGraphErrors *gr_val2 = new TGraphErrors(n,&x_ptthree_[0],&y_mean_ratio_val2_[0],&ex_ptthree_[0],&ey_mean_ratio_val2_[0]);
      gr_val2->SetTitle(eta_bins_labels_[eta_i].first +" < #eta < " + eta_bins_labels_[eta_i].second);
      gr_val2->SetName((TString)"kFSR_extrapol_val2_eta_"+(Long_t)eta_i);
      //      TGraphErrorsstyle(gr_val2,line_styles_, colours_, markers_, 0,0,0,"bla");


      TGraphErrors *gr_res2 = new TGraphErrors(n,&x_ptthree_[0],&y_mean_ratio_res2_[0],&ex_ptthree_[0],&ey_mean_ratio_res2_[0]);
      gr_res2->SetTitle(eta_bins_labels_[eta_i].first +" < #eta < " + eta_bins_labels_[eta_i].second);
      gr_res2->SetName((TString)"kFSR_extrapol_res2_eta_"+(Long_t)eta_i);
      //      TGraphErrorsstyle(gr_res2,line_styles_, colours_, markers_, 1,1,1,"bla");


     TF1 *lin_extrapol = new TF1("lin_extrapol","[0]+[1]*x",0,0.5); //was used before...
      lin_extrapol->SetParameters(1,-0.1);
      lin_extrapol->SetParName(0,"k_{FSR}");
      lin_extrapol->SetParName(1,"slope");

      //deactivate because JINST
      //      lin_extrapol->SetLineColor(gr_val1->GetLineColor());
      lin_extrapol->SetLineColor(2);
      gr_val1->Fit("lin_extrapol","","same",0,0.5);

      //      lin_extrapol->SetLineColor(gr_res1->GetLineColor());
      gr_res1->Fit("lin_extrapol","","same",0,0.5);

      extrapol_inclusive_eta_val1_.push_back(gr_val1);
      extrapol_inclusive_eta_res1_.push_back(gr_res1);

      //      lin_extrapol->SetLineColor(gr_val2->GetLineColor());
      gr_val2->Fit("lin_extrapol","","same",0,0.5);

      //      lin_extrapol->SetLineColor(gr_res2->GetLineColor());
      gr_res2->Fit("lin_extrapol","","same",0,0.5);

      extrapol_inclusive_eta_val2_.push_back(gr_val2);
      extrapol_inclusive_eta_res2_.push_back(gr_res2);


    }


     if(export_all_plots)draw_TGraphErrors_save_PS(fine_coarse+"_"+ generatorone+"_"+generatortwo+"_"+jet_type+"_PS_kFSR_extrapols",image_ext,extrapol_inclusive_eta_res1_,"testen_eta_res1","tdr","P","no","x0",0.,0.45,0.95,1.05,"cut on #alpha","Normalized residual correction",true);//"#alpha","#frac{<R(MC)/R(data)>_{#alpha}}{<R(MC)/R(data)>_{#alpha=0.2}}");
    // save memory
    //     if(export_all_plots)draw_TGraphErrors_save_PS(fine_coarse+"_"+ generatorone+"_"+generatortwo+"_"+jet_type+"_PS_kFSR_extrapols",image_ext,extrapol_inclusive_eta_res2_,"testen_eta_res2","tdr","P","no","x0",0.,0.45,0.95,1.05,"#alpha","Normalized residual correction");//"#alpha","#frac{<R(MC)/R(data)>_{#alpha}}{<R(MC)/R(data)>_{#alpha=0.2}}");
     if(export_all_plots)draw_TGraphErrors_save_PS(fine_coarse+"_"+ generatorone+"_"+generatortwo+"_"+jet_type+"_PS_kFSR_extrapols",image_ext,extrapol_inclusive_eta_val1_,"testen_eta_val1","tdr","P","no","x0",0.,0.45,0.95,1.05,"cut on #alpha","Normalized residual correction",true);//"#alpha","#frac{<R(MC)/R(data)>_{#alpha}}{<R(MC)/R(data)>_{#alpha=0.2}}");
    // save memory
    //     if(export_all_plots)draw_TGraphErrors_save_PS(fine_coarse+"_"+ generatorone+"_"+generatortwo+"_"+jet_type+"_PS_kFSR_extrapols",image_ext,extrapol_inclusive_eta_val2_,"testen_eta_val2","tdr","P","no","x0",0.,0.45,0.95,1.05,"#alpha","Normalized residual correction");//"#alpha","#frac{<R(MC)/R(data)>_{#alpha}}{<R(MC)/R(data)>_{#alpha=0.2}}");





    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////Extract kRAD*MC/DATA results//////////////////////////////////////////////////////////////
    //////////////////////////////Extract kRAD*MC/DATA results//////////////////////////////////////////////////////////////
    //////////////////////////////Extract kRAD*MC/DATA results//////////////////////////////////////////////////////////////
    //////////////////////////////Extract kRAD*MC/DATA results//////////////////////////////////////////////////////////////
    //////////////////////////////Extract kRAD*MC/DATA results//////////////////////////////////////////////////////////////
    //////////////////////////////Extract kRAD*MC/DATA results//////////////////////////////////////////////////////////////
    //////////////////////////////Extract kRAD*MC/DATA results//////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


     std::cout << "still works" << std::endl;

  std::vector <Double_t> SLOPE_Abseta_;
  std::vector <Double_t> eSLOPE_Abseta_;
  std::vector <Double_t> CONST_SLOPE_Abseta_;
  std::vector <Double_t> eCONST_SLOPE_Abseta_;

  std::vector <Double_t> x_Abseta_;
  std::vector <Double_t> ex_Abseta_;
  std::vector <Double_t> trad_x_Abseta_;
  std::vector <Double_t> trad_ex_Abseta_;
  std::vector <Double_t> y_kFSR_res1_;
  std::vector <Double_t> ey_kFSR_res1_;
  std::vector <Double_t> y_kFSR_val1_;
  std::vector <Double_t> ey_kFSR_val1_;
  std::vector <Double_t> y_residual_res1_;
  std::vector <Double_t> ey_residual_res1_;
  std::vector <Double_t> y_validation_val1_;
  std::vector <Double_t> ey_validation_val1_;
  std::vector <Double_t> y_kFSR_res2_;
  std::vector <Double_t> ey_kFSR_res2_;
  std::vector <Double_t> y_kFSR_val2_;
  std::vector <Double_t> ey_kFSR_val2_;
  std::vector <Double_t> y_residual_res2_;
  std::vector <Double_t> ey_residual_res2_;
  std::vector <Double_t> y_validation_val2_;
  std::vector <Double_t> ey_validation_val2_;

  for(Int_t Abseta_i=0;Abseta_i<trad_no_Abseta_bins;Abseta_i++){
    trad_x_Abseta_.push_back((trad_Abseta_bins_[Abseta_i].second+trad_Abseta_bins_[Abseta_i].first)/2);
    trad_ex_Abseta_.push_back((trad_Abseta_bins_[Abseta_i].second-trad_Abseta_bins_[Abseta_i].first)/2);
  }
  for(Int_t Abseta_i=0;Abseta_i<no_Abseta_bins;Abseta_i++){

    x_Abseta_.push_back((Abseta_bins_[Abseta_i].second+Abseta_bins_[Abseta_i].first)/2);
    ex_Abseta_.push_back((Abseta_bins_[Abseta_i].second-Abseta_bins_[Abseta_i].first)/2);

    Int_t cut_i=1;

      for(unsigned int cut_j=0;cut_j<ptthreecuts.size();cut_j++){
	if(ptthreecuts[cut_j]=="20"){
	  cut_i=cut_j;}
	}
std::cout << "cut_i: " << cut_i << "abseta_i" << Abseta_i << std::endl;
      SLOPE_Abseta_.push_back(all_ptthree_all_Abseta_ratio_res1_  .at(cut_i).at(Abseta_i)->GetFunction("fit_loglin")->GetParameter(1));
      eSLOPE_Abseta_.push_back(all_ptthree_all_Abseta_ratio_res1_  .at(cut_i).at(Abseta_i)->GetFunction("fit_loglin")->GetParError(1));
      CONST_SLOPE_Abseta_.push_back(all_ptthree_all_Abseta_ratio_res1_  .at(cut_i).at(Abseta_i)->GetFunction("fit_loglin")->GetParameter(0));
      eCONST_SLOPE_Abseta_.push_back(all_ptthree_all_Abseta_ratio_res1_  .at(cut_i).at(Abseta_i)->GetFunction("fit_loglin")->GetParError(0));

//     std::cout << "still works" << std::endl;
//     std::cout << Abseta_i <<std::endl;
    //all_ptthree_all_Abseta_ratio_val1_  all_ptthree_all_Abseta_MC1_L2L3_
    Double_t corr_residual_res1 = all_ptthree_all_Abseta_ratio_res1_  .at(cut_i).at(Abseta_i)->GetFunction("fit_const")->GetParameter(0);
    Double_t corr_validation_val1 = all_ptthree_all_Abseta_ratio_val1_  .at(cut_i).at(Abseta_i)->GetFunction("fit_const")->GetParameter(0);
    Double_t kFSR_res1 = extrapol_inclusive_Abseta_res1_.at(Abseta_i)->GetFunction("lin_extrapol")->GetParameter(0);
    Double_t kFSR_val1 = extrapol_inclusive_Abseta_val1_.at(Abseta_i)->GetFunction("lin_extrapol")->GetParameter(0);

    Double_t e_corr_residual_res1 = all_ptthree_all_Abseta_ratio_res1_  .at(cut_i).at(Abseta_i)->GetFunction("fit_const")->GetParError(0);
    Double_t e_corr_validation_val1 = all_ptthree_all_Abseta_ratio_val1_  .at(cut_i).at(Abseta_i)->GetFunction("fit_const")->GetParError(0);
    Double_t e_kFSR_res1 = extrapol_inclusive_Abseta_res1_.at(Abseta_i)->GetFunction("lin_extrapol")->GetParError(0);
    //    cout << "which is bigger... direkt kFSR error: " << extrapol_inclusive_Abseta_res1_.at(Abseta_i)->GetFunction("lin_extrapol")->GetParError(0) << " or 0.2 times slopeerror: " << 0.2*extrapol_inclusive_Abseta_res1_.at(Abseta_i)->GetFunction("lin_extrapol")->GetParError(1) << endl;;
    Double_t e_kFSR_val1 = extrapol_inclusive_Abseta_val1_.at(Abseta_i)->GetFunction("lin_extrapol")->GetParError(0);

    Double_t corr_residual_res2 = all_ptthree_all_Abseta_ratio_res2_  .at(cut_i).at(Abseta_i)->GetFunction("fit_const")->GetParameter(0);
    Double_t corr_validation_val2 = all_ptthree_all_Abseta_ratio_val2_  .at(cut_i).at(Abseta_i)->GetFunction("fit_const")->GetParameter(0);
    Double_t kFSR_res2 = extrapol_inclusive_Abseta_res2_.at(Abseta_i)->GetFunction("lin_extrapol")->GetParameter(0);
    Double_t kFSR_val2 = extrapol_inclusive_Abseta_val2_.at(Abseta_i)->GetFunction("lin_extrapol")->GetParameter(0);

    Double_t e_corr_residual_res2 = all_ptthree_all_Abseta_ratio_res2_  .at(cut_i).at(Abseta_i)->GetFunction("fit_const")->GetParError(0);
    Double_t e_corr_validation_val2 = all_ptthree_all_Abseta_ratio_val2_  .at(cut_i).at(Abseta_i)->GetFunction("fit_const")->GetParError(0);
    Double_t e_kFSR_res2 = extrapol_inclusive_Abseta_res2_.at(Abseta_i)->GetFunction("lin_extrapol")->GetParError(0);
    Double_t e_kFSR_val2 = extrapol_inclusive_Abseta_val2_.at(Abseta_i)->GetFunction("lin_extrapol")->GetParError(0);


    //    if(use_imported_kFSRAbs.Contains("true")){

    if(use_imported_kFSRAbs.Contains("true")){
      Double_t eta_search= std::abs(x_Abseta_.back());
      Double_t import_kFSR= import_kFSR_vs_Abseta_histo_res1->GetBinContent(import_kFSR_vs_Abseta_histo_res1->FindBin(eta_search));
      Double_t import_ekFSR= import_kFSR_vs_Abseta_histo_res1->GetBinError(import_kFSR_vs_Abseta_histo_res1->FindBin(eta_search));
      kFSR_res1=import_kFSR;
      kFSR_res2=import_kFSR;
      kFSR_val1=import_kFSR;
      kFSR_val2=import_kFSR;

      e_kFSR_res1=import_ekFSR;
      e_kFSR_res2=import_ekFSR;
      e_kFSR_val1=import_ekFSR;
      e_kFSR_val2=import_ekFSR;
    }

    //
    if(use_fitted_kFSR.Contains("use_fitted_kFSR")){
     
      cout << "imported_kFSR_func test..."<< x_Abseta_.back() <<endl;
      cout << "bla value: " << import_kFSR_fit->Eval(x_Abseta_.back()) << endl;
      Double_t import_kFSR = import_kFSR_fit->Eval(x_Abseta_.back());
      cout << "imported_kFSR_func: " << import_kFSR <<endl;

      Double_t import_ekFSR= get_error(x_Abseta_.back(), cov, import_kFSR_fit);
      cout << "import_kFSR: " << import_kFSR << " import_ekFSR: " << import_ekFSR << endl;
      //           Double_t import_ekFSR=0.001;
      kFSR_res1=import_kFSR;
      kFSR_res2=import_kFSR;
      kFSR_val1=import_kFSR;
      kFSR_val2=import_kFSR;

      e_kFSR_res1=import_ekFSR;
      e_kFSR_res2=import_ekFSR;
      e_kFSR_val1=import_ekFSR;
      e_kFSR_val2=import_ekFSR;
    }
    if(kFSR_eq_one.Contains("kFSR_eq_one")){
      kFSR_res1=1.0;
      kFSR_res2=1.0;
      kFSR_val1=1.0;
      kFSR_val2=1.0;
    }

      y_residual_res1_.push_back(corr_residual_res1*kFSR_res1);
      ey_residual_res1_.push_back(TMath::Sqrt(TMath::Power(e_corr_residual_res1,2)+TMath::Power(e_kFSR_res1,2)));

      y_validation_val1_.push_back(corr_validation_val1*kFSR_val1);
      ey_validation_val1_.push_back(TMath::Sqrt(TMath::Power(e_corr_validation_val1,2)+TMath::Power(e_kFSR_val1,2)));

      y_kFSR_res1_.push_back(kFSR_res1);
      y_kFSR_val1_.push_back(kFSR_val1);

      ey_kFSR_res1_.push_back(e_kFSR_res1);
      ey_kFSR_val1_.push_back(e_kFSR_val1);


      y_residual_res2_.push_back(corr_residual_res2*kFSR_res2);
      ey_residual_res2_.push_back(TMath::Sqrt(TMath::Power(e_corr_residual_res2,2)+TMath::Power(e_kFSR_res2,2)));

      y_validation_val2_.push_back(corr_validation_val2*kFSR_val2);
      ey_validation_val2_.push_back(TMath::Sqrt(TMath::Power(e_corr_validation_val2,2)+TMath::Power(e_kFSR_val2,2)));

      y_kFSR_res2_.push_back(kFSR_res2);
      y_kFSR_val2_.push_back(kFSR_val2);

      ey_kFSR_res2_.push_back(e_kFSR_res2);
      ey_kFSR_val2_.push_back(e_kFSR_val2);


      }


      Int_t n = y_residual_res1_.size();
      //      Int_t abs_n = y_residual_res1_.size();
      //      cout << "wie oft...: " << n << endl;

      TH1D* Residual_slope_Abseta_correction_histo_res1=new TH1D("Residual_slope_Abseta_correction_histo_res1", "Residual_slope_Abseta_correction_histo_res1", zero_eta , &eta_binning[zero_eta]);
      TH1D* Residual_const_slope_Abseta_correction_histo_res1=new TH1D("Residual_const_slope_Abseta_correction_histo_res1", "Residual_const_slope_Abseta_correction_histo_res1", zero_eta , &eta_binning[zero_eta]);

      TH1D* Residual_Abseta_correction_histo_res1=new TH1D("Residual_Abseta_correction_histo_res1", "Residual_Abseta_correction_histo_res1", zero_eta , &eta_binning[zero_eta]);
      TH1D* Residual_Abseta_correction_histo_val1=new TH1D("Residual_Abseta_correction_histo_val1", "Residual_Abseta_correction_histo_val1", zero_eta , &eta_binning[zero_eta]);
      TH1D* Residual_Abseta_correction_histo_res2=new TH1D("Residual_Abseta_correction_histo_res2", "Residual_Abseta_correction_histo_res2", zero_eta , &eta_binning[zero_eta]);
      TH1D* Residual_Abseta_correction_histo_val2=new TH1D("Residual_Abseta_correction_histo_val2", "Residual_Abseta_correction_histo_val2", zero_eta , &eta_binning[zero_eta]);

      //eta_binning[]={-6.0,-4.0,-3.5,-3.2,-3.0,-2.8,-2.6,-2.4,-2.2,-2.0,-1.8,-1.5,-1.4,-1.3,-1.2,-1.1,-0.9,-0.6,-0.3,0.0,0.3,0.6,0.9,1.1,1.2,1.3,1.4,1.5,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.2,3.5,4.0,6.0};

      //  Int_t  no_eta =38;
      //  Int_t  zero_eta =19;


      //      Residual_Abseta_correction_histo_res1->Print("all");
  for(Int_t Abseta_i=0;Abseta_i<zero_eta;Abseta_i++){
//     std::cout << "still works" << std::endl;
//     std::cout << Abseta_i <<std::endl;
    Residual_slope_Abseta_correction_histo_res1->SetBinContent(Abseta_i+1,SLOPE_Abseta_[Abseta_i]);
    //    cout << "x: " << x_Abseta_[Abseta_i] << " and y: "<< SLOPE_Abseta_[Abseta_i] << " and error: " << eSLOPE_Abseta_[Abseta_i] << endl;
    Residual_slope_Abseta_correction_histo_res1->SetBinError(Abseta_i+1,eSLOPE_Abseta_[Abseta_i]);
    Residual_slope_Abseta_correction_histo_res1->SetTitle(Abseta_bins_labels_[Abseta_i].first +" < |#eta| < " + Abseta_bins_labels_[Abseta_i].second);
    Residual_const_slope_Abseta_correction_histo_res1->SetBinContent(Abseta_i+1,CONST_SLOPE_Abseta_[Abseta_i]);
    //    cout << "x: " << x_Abseta_[Abseta_i] << " and y: "<< CONST_SLOPE_Abseta_[Abseta_i] << " and error: " << eCONST_SLOPE_Abseta_[Abseta_i] << endl;
    Residual_const_slope_Abseta_correction_histo_res1->SetBinError(Abseta_i+1,eCONST_SLOPE_Abseta_[Abseta_i]);
    Residual_const_slope_Abseta_correction_histo_res1->SetTitle(Abseta_bins_labels_[Abseta_i].first +" < |#eta| < " + Abseta_bins_labels_[Abseta_i].second);

    Residual_Abseta_correction_histo_res1->SetBinContent(Abseta_i+1,y_residual_res1_[Abseta_i]);
    //    cout << "x: " << x_Abseta_[Abseta_i] << " and y: "<< y_residual_res1_[Abseta_i] << " and error: " << ey_residual_res1_[Abseta_i] << endl;
    Residual_Abseta_correction_histo_res1->SetBinError(Abseta_i+1,ey_residual_res1_[Abseta_i]);
    Residual_Abseta_correction_histo_val1->SetBinContent(Abseta_i+1,y_validation_val1_[Abseta_i]);
    //    cout << "x: " << x_Abseta_[Abseta_i] << " and y: "<< y_validation_val1_[Abseta_i] << " and error: " << ey_validation_val1_[Abseta_i] << endl;
    Residual_Abseta_correction_histo_val1->SetBinError(Abseta_i+1,ey_validation_val1_[Abseta_i]);

    Residual_Abseta_correction_histo_res2->SetBinContent(Abseta_i+1,y_residual_res2_[Abseta_i]);
    //    cout << "x: " << x_Abseta_[Abseta_i] << " and y: "<< y_residual_res2_[Abseta_i] << " and error: " << ey_residual_res2_[Abseta_i] << endl;
    Residual_Abseta_correction_histo_res2->SetBinError(Abseta_i+1,ey_residual_res2_[Abseta_i]);
    Residual_Abseta_correction_histo_val2->SetBinContent(Abseta_i+1,y_validation_val2_[Abseta_i]);
    //    cout << "x: " << x_Abseta_[Abseta_i] << " and y: "<< y_validation_val2_[Abseta_i] << " and error: " << ey_validation_val2_[Abseta_i] << endl;
    Residual_Abseta_correction_histo_val2->SetBinError(Abseta_i+1,ey_validation_val2_[Abseta_i]);
  }

  TGraphErrors* Residual_Abseta_correction_res1 = new TGraphErrors(n,&x_Abseta_[0],&y_residual_res1_[0],&ex_Abseta_[0],&ey_residual_res1_[0]);
  Residual_Abseta_correction_res1->SetTitle("Residual correction (res1 constants)");
  Residual_Abseta_correction_res1->SetName("Residual_Abseta_correction_res1");
  TGraphErrorsstyle(Residual_Abseta_correction_res1,line_styles_, colours_, markers_, 1,1,1,"#eta(+) - #eta(-)");
  TGraphErrors* Residual_Abseta_correction_val1 = new TGraphErrors(n,&x_Abseta_[0],&y_validation_val1_[0],&ex_Abseta_[0],&ey_validation_val1_[0]);
  Residual_Abseta_correction_val1->SetTitle("Residual correction (on top of L2L3res)");
  Residual_Abseta_correction_val1->SetName("Residual_Abseta_correction_val1");
  TGraphErrorsstyle(Residual_Abseta_correction_val1,line_styles_, colours_, markers_, 0,0,0,"#eta(+) - #eta(-)");

  TH1D* kFSR_vs_Abseta_histo_res1=new TH1D("kFSR_vs_Abseta_histo_res1", "kFSR_vs_Abseta_histo_res1", zero_eta , &eta_binning[zero_eta]);
  TH1D* kFSR_vs_Abseta_histo_res2=new TH1D("kFSR_vs_Abseta_histo_res2", "kFSR_vs_Abseta_histo_res2", zero_eta , &eta_binning[zero_eta]);
  //      Residual_Abseta_correction_histo_res1->Print("all");
  for(Int_t Abseta_i=0;Abseta_i<zero_eta;Abseta_i++){
     std::cout << "still works" << std::endl;
     std::cout << Abseta_i <<std::endl;

    kFSR_vs_Abseta_histo_res1->SetBinContent(Abseta_i+1,y_kFSR_res1_[Abseta_i]);
    kFSR_vs_Abseta_histo_res1->SetBinError(Abseta_i+1,ey_kFSR_res1_[Abseta_i]);
    kFSR_vs_Abseta_histo_res2->SetBinContent(Abseta_i+1,y_kFSR_res2_[Abseta_i]);
    kFSR_vs_Abseta_histo_res2->SetBinError(Abseta_i+1,ey_kFSR_res2_[Abseta_i]);
  }


  TGraphErrors* kFSR_vs_Abseta_res1 = new TGraphErrors(n,&x_Abseta_[0],&y_kFSR_res1_[0],&ex_Abseta_[0],&ey_kFSR_res1_[0]);
  kFSR_vs_Abseta_res1->SetTitle("kFSR (res1)");
  kFSR_vs_Abseta_res1->SetName("kFSR_vs_Abseta_res1");
  TGraphErrorsstyle(kFSR_vs_Abseta_res1,line_styles_, colours_, markers_, 1,1,1,"#eta(+) - #eta(-)");
  TGraphErrors* kFSR_vs_Abseta_val1 = new TGraphErrors(n,&x_Abseta_[0],&y_kFSR_val1_[0],&ex_Abseta_[0],&ey_kFSR_val1_[0]);
  kFSR_vs_Abseta_val1->SetTitle("kFSR (val1idation)");
  kFSR_vs_Abseta_val1->SetName("kFSR_vs_Abseta_val1");
  TGraphErrorsstyle(kFSR_vs_Abseta_val1,line_styles_, colours_, markers_, 0,0,0,"bla");
  TGraphErrors* Residual_Abseta_correction_res2 = new TGraphErrors(n,&x_Abseta_[0],&y_residual_res2_[0],&ex_Abseta_[0],&ey_residual_res2_[0]);
  Residual_Abseta_correction_res2->SetTitle("Residual correction (res2 constants)");
  Residual_Abseta_correction_res2->SetName("Residual_Abseta_correction_res2");
  TGraphErrorsstyle(Residual_Abseta_correction_res2,line_styles_, colours_, markers_, 1,1,1,"bla");
  TGraphErrors* Residual_Abseta_correction_val2 = new TGraphErrors(n,&x_Abseta_[0],&y_validation_val2_[0],&ex_Abseta_[0],&ey_validation_val2_[0]);
  Residual_Abseta_correction_val2->SetTitle("Residual correction (on top of L2L3res)");
  Residual_Abseta_correction_val2->SetName("Residual_Abseta_correction_val2");
  TGraphErrorsstyle(Residual_Abseta_correction_val2,line_styles_, colours_, markers_, 0,0,0,"bla");
  TGraphErrors* kFSR_vs_Abseta_res2 = new TGraphErrors(n,&x_Abseta_[0],&y_kFSR_res2_[0],&ex_Abseta_[0],&ey_kFSR_res2_[0]);
  kFSR_vs_Abseta_res2->SetTitle("kFSR (res2)");
  kFSR_vs_Abseta_res2->SetName("kFSR_vs_Abseta_res2");
  TGraphErrorsstyle(kFSR_vs_Abseta_res2,line_styles_, colours_, markers_, 1,1,1,"#eta(+) - #eta(-)");
  TGraphErrors* kFSR_vs_Abseta_val2 = new TGraphErrors(n,&x_Abseta_[0],&y_kFSR_val2_[0],&ex_Abseta_[0],&ey_kFSR_val2_[0]);
  kFSR_vs_Abseta_val2->SetTitle("kFSR (val1idation)");
  kFSR_vs_Abseta_val2->SetName("kFSR_vs_Abseta_val2");
  TGraphErrorsstyle(kFSR_vs_Abseta_val2,line_styles_, colours_, markers_, 0,0,0,"#eta(+) - #eta(-)");


  std::vector <Double_t> SLOPE_eta_;
  std::vector <Double_t> eSLOPE_eta_;
  std::vector <Double_t> CONST_SLOPE_eta_;
  std::vector <Double_t> eCONST_SLOPE_eta_;

  std::vector <Double_t> x_eta_;
  std::vector <Double_t> ex_eta_;
  std::vector <Double_t> trad_x_eta_;
  std::vector <Double_t> trad_ex_eta_;
  y_kFSR_res1_.clear();
  ey_kFSR_res1_.clear();
  y_kFSR_val1_.clear();
  ey_kFSR_val1_.clear();
  y_residual_res1_.clear();
  ey_residual_res1_.clear();
  y_validation_val1_.clear();
  ey_validation_val1_.clear();
  y_kFSR_res2_.clear();
  ey_kFSR_res2_.clear();
  y_kFSR_val2_.clear();
  ey_kFSR_val2_.clear();
  y_residual_res2_.clear();
  ey_residual_res2_.clear();
  y_validation_val2_.clear();
  ey_validation_val2_.clear();





  for(Int_t eta_i=0;eta_i<trad_eta_bins;eta_i++){
    trad_x_eta_.push_back((trad_eta_bins_[eta_i].second+trad_eta_bins_[eta_i].first)/2);
    trad_ex_eta_.push_back((trad_eta_bins_[eta_i].second-trad_eta_bins_[eta_i].first)/2);
  }
  for(Int_t eta_i=0;eta_i<eta_bins;eta_i++){

    x_eta_.push_back((eta_bins_[eta_i].second+eta_bins_[eta_i].first)/2);
    ex_eta_.push_back((eta_bins_[eta_i].second-eta_bins_[eta_i].first)/2);

    Int_t cut_i=1;

      for(unsigned int cut_j=0;cut_j<ptthreecuts.size();cut_j++){
	if(ptthreecuts[cut_j]=="20"){
	  cut_i=cut_j;}
	}

      std::cout << "still works" << std::endl;
     std::cout << eta_i <<std::endl;
     if(all_ptthree_all_eta_ratio_res1_  .at(cut_i).at(eta_i)->GetFunction("fit_loglin")){
     SLOPE_eta_.push_back(all_ptthree_all_eta_ratio_res1_  .at(cut_i).at(eta_i)->GetFunction("fit_loglin")->GetParameter(1));
      eSLOPE_eta_.push_back(all_ptthree_all_eta_ratio_res1_  .at(cut_i).at(eta_i)->GetFunction("fit_loglin")->GetParError(1));
      CONST_SLOPE_eta_.push_back(all_ptthree_all_eta_ratio_res1_  .at(cut_i).at(eta_i)->GetFunction("fit_loglin")->GetParameter(0));
      eCONST_SLOPE_eta_.push_back(all_ptthree_all_eta_ratio_res1_  .at(cut_i).at(eta_i)->GetFunction("fit_loglin")->GetParError(0));
     }
     else {
       SLOPE_eta_.push_back(0);
       eSLOPE_eta_.push_back(0);
       CONST_SLOPE_eta_.push_back(0);
       eCONST_SLOPE_eta_.push_back(0);
       
     }
      //    Double_t corr_residual_const_slope_res1 = all_ptthree_all_eta_ratio_res1_  .at(cut_i).at(eta_i)->GetFunction("fit_loglin")->GetParameter(1);
      //    Double_t e_corr_residual_const_slope_res1 = all_ptthree_all_eta_ratio_res1_  .at(cut_i).at(eta_i)->GetFunction("fit_loglin")->GetParError(1);

    //all_ptthree_all_eta_ratio_val1_  all_ptthree_all_eta_MC1_L2L3_
      std::cout << "still works" << std::endl;
     std::cout << eta_i <<std::endl;
    Double_t corr_residual_res1 = all_ptthree_all_eta_ratio_res1_  .at(cut_i).at(eta_i)->GetFunction("fit_const")->GetParameter(0);
    Double_t corr_validation_val1 = all_ptthree_all_eta_ratio_val1_  .at(cut_i).at(eta_i)->GetFunction("fit_const")->GetParameter(0);
    Double_t kFSR_res1 = extrapol_inclusive_eta_res1_.at(eta_i)->GetFunction("lin_extrapol")->GetParameter(0);
    Double_t kFSR_val1 = extrapol_inclusive_eta_val1_.at(eta_i)->GetFunction("lin_extrapol")->GetParameter(0);

    Double_t e_corr_residual_res1 = all_ptthree_all_eta_ratio_res1_  .at(cut_i).at(eta_i)->GetFunction("fit_const")->GetParError(0);
    Double_t e_corr_validation_val1 = all_ptthree_all_eta_ratio_val1_  .at(cut_i).at(eta_i)->GetFunction("fit_const")->GetParError(0);
    Double_t e_kFSR_res1 = extrapol_inclusive_eta_res1_.at(eta_i)->GetFunction("lin_extrapol")->GetParError(0);
    Double_t e_kFSR_val1 = extrapol_inclusive_eta_val1_.at(eta_i)->GetFunction("lin_extrapol")->GetParError(0);

    Double_t corr_residual_res2 = all_ptthree_all_eta_ratio_res2_  .at(cut_i).at(eta_i)->GetFunction("fit_const")->GetParameter(0);
    Double_t corr_validation_val2 = all_ptthree_all_eta_ratio_val2_  .at(cut_i).at(eta_i)->GetFunction("fit_const")->GetParameter(0);
    Double_t kFSR_res2 = extrapol_inclusive_eta_res2_.at(eta_i)->GetFunction("lin_extrapol")->GetParameter(0);
    Double_t kFSR_val2 = extrapol_inclusive_eta_val2_.at(eta_i)->GetFunction("lin_extrapol")->GetParameter(0);

    Double_t e_corr_residual_res2 = all_ptthree_all_eta_ratio_res2_  .at(cut_i).at(eta_i)->GetFunction("fit_const")->GetParError(0);
    Double_t e_corr_validation_val2 = all_ptthree_all_eta_ratio_val2_  .at(cut_i).at(eta_i)->GetFunction("fit_const")->GetParError(0);
    Double_t e_kFSR_res2 = extrapol_inclusive_eta_res2_.at(eta_i)->GetFunction("lin_extrapol")->GetParError(0);
    Double_t e_kFSR_val2 = extrapol_inclusive_eta_val2_.at(eta_i)->GetFunction("lin_extrapol")->GetParError(0);


    if(use_imported_kFSRAbs.Contains("true")){
      Double_t eta_search= std::abs(x_eta_.back());
      Double_t import_kFSR= import_kFSR_vs_Abseta_histo_res1->GetBinContent(import_kFSR_vs_Abseta_histo_res1->FindBin(eta_search));
      Double_t import_ekFSR= import_kFSR_vs_Abseta_histo_res1->GetBinError(import_kFSR_vs_Abseta_histo_res1->FindBin(eta_search));
      kFSR_res1=import_kFSR;
      kFSR_res2=import_kFSR;
      kFSR_val1=import_kFSR;
      kFSR_val2=import_kFSR;

      e_kFSR_res1=import_ekFSR;
      e_kFSR_res2=import_ekFSR;
      e_kFSR_val1=import_ekFSR;
      e_kFSR_val2=import_ekFSR;
    }

    if(use_fitted_kFSR.Contains("use_fitted_kFSR")){
     
      Double_t import_kFSR = import_kFSR_fit->Eval(std::abs(x_eta_.back()));
      cout << "imported_kFSR_func: " << import_kFSR <<endl;
      Double_t import_ekFSR= get_error(std::abs(x_eta_.back()), cov, import_kFSR_fit);
      //            Double_t import_ekFSR=0.001;
      kFSR_res1=import_kFSR;
      kFSR_res2=import_kFSR;
      kFSR_val1=import_kFSR;
      kFSR_val2=import_kFSR;

      e_kFSR_res1=import_ekFSR;
      e_kFSR_res2=import_ekFSR;
      e_kFSR_val1=import_ekFSR;
      e_kFSR_val2=import_ekFSR;
    }

    if(kFSR_eq_one.Contains("kFSR_eq_one")){
      kFSR_res1=1.0;
      kFSR_res2=1.0;
      kFSR_val1=1.0;
      kFSR_val2=1.0;
    }

      y_residual_res1_.push_back(corr_residual_res1*kFSR_res1);
      ey_residual_res1_.push_back(TMath::Sqrt(TMath::Power(e_corr_residual_res1,2)+TMath::Power(e_kFSR_res1,2)));

      y_validation_val1_.push_back(corr_validation_val1*kFSR_val1);
      ey_validation_val1_.push_back(TMath::Sqrt(TMath::Power(e_corr_validation_val1,2)+TMath::Power(e_kFSR_val1,2)));

      y_kFSR_res1_.push_back(kFSR_res1);
      y_kFSR_val1_.push_back(kFSR_val1);

      ey_kFSR_res1_.push_back(e_kFSR_res1);
      ey_kFSR_val1_.push_back(e_kFSR_val1);


      y_residual_res2_.push_back(corr_residual_res2*kFSR_res2);
      ey_residual_res2_.push_back(TMath::Sqrt(TMath::Power(e_corr_residual_res2,2)+TMath::Power(e_kFSR_res2,2)));

      y_validation_val2_.push_back(corr_validation_val2*kFSR_val2);
      ey_validation_val2_.push_back(TMath::Sqrt(TMath::Power(e_corr_validation_val2,2)+TMath::Power(e_kFSR_val2,2)));

      y_kFSR_res2_.push_back(kFSR_res2);
      y_kFSR_val2_.push_back(kFSR_val2);

      ey_kFSR_res2_.push_back(e_kFSR_res2);
      ey_kFSR_val2_.push_back(e_kFSR_val2);


      
}

      n = y_residual_res1_.size();
      //      cout << "wie oft...: " << n << endl;


      TH1D* Residual_slope_eta_correction_histo_res1=new TH1D("Residual_slope_eta_correction_histo_res1", "Residual_slope_eta_correction_histo_res1", no_eta , &eta_binning[0]);
      TH1D* Residual_const_slope_eta_correction_histo_res1=new TH1D("Residual_const_slope_eta_correction_histo_res1", "Residual_const_slope_eta_correction_histo_res1", no_eta , &eta_binning[0]);


      TH1D* Residual_eta_correction_histo_res1=new TH1D("Residual_eta_correction_histo_res1", "Residual_eta_correction_histo_res1", no_eta , &eta_binning[0]);
      TH1D* Residual_eta_correction_histo_val1=new TH1D("Residual_eta_correction_histo_val1", "Residual_eta_correction_histo_val1", no_eta , &eta_binning[0]);
      TH1D* Residual_eta_correction_histo_res2=new TH1D("Residual_eta_correction_histo_res2", "Residual_eta_correction_histo_res2", no_eta , &eta_binning[0]);
      TH1D* Residual_eta_correction_histo_val2=new TH1D("Residual_eta_correction_histo_val2", "Residual_eta_correction_histo_val2", no_eta , &eta_binning[0]);

      //eta_binning[]={-6.0,-4.0,-3.5,-3.2,-3.0,-2.8,-2.6,-2.4,-2.2,-2.0,-1.8,-1.5,-1.4,-1.3,-1.2,-1.1,-0.9,-0.6,-0.3,0.0,0.3,0.6,0.9,1.1,1.2,1.3,1.4,1.5,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.2,3.5,4.0,6.0};

      //  Int_t  no_eta =38;
      //  Int_t  zero_eta =19;


      //      Residual_eta_correction_histo_res1->Print("all");
  for(Int_t eta_i=0;eta_i<no_eta;eta_i++){
    Residual_slope_eta_correction_histo_res1->SetBinContent(eta_i+1,SLOPE_eta_[eta_i]);
    cout << "x: " << x_eta_[eta_i] << " and y: "<< SLOPE_eta_[eta_i] << " and error: " << eSLOPE_eta_[eta_i] << endl;
    Residual_slope_eta_correction_histo_res1->SetBinError(eta_i+1,eSLOPE_eta_[eta_i]);
    Residual_slope_eta_correction_histo_res1->SetTitle(eta_bins_labels_[eta_i].first +" < #eta < " + eta_bins_labels_[eta_i].second);
    Residual_const_slope_eta_correction_histo_res1->SetBinContent(eta_i+1,CONST_SLOPE_eta_[eta_i]);
    cout << "x: " << x_eta_[eta_i] << " and y: "<< CONST_SLOPE_eta_[eta_i] << " and error: " << eCONST_SLOPE_eta_[eta_i] << endl;
    Residual_const_slope_eta_correction_histo_res1->SetBinError(eta_i+1,eCONST_SLOPE_eta_[eta_i]);
    Residual_const_slope_eta_correction_histo_res1->SetTitle(eta_bins_labels_[eta_i].first +" < #eta < " + eta_bins_labels_[eta_i].second);


    Residual_eta_correction_histo_res1->SetBinContent(eta_i+1,y_residual_res1_[eta_i]);
    cout << "x: " << x_eta_[eta_i] << " and y: "<< y_residual_res1_[eta_i] << " and error: " << ey_residual_res1_[eta_i] << endl;
    Residual_eta_correction_histo_res1->SetBinError(eta_i+1,ey_residual_res1_[eta_i]);
    Residual_eta_correction_histo_val1->SetBinContent(eta_i+1,y_validation_val1_[eta_i]);
    cout << "x: " << x_eta_[eta_i] << " and y: "<< y_validation_val1_[eta_i] << " and error: " << ey_validation_val1_[eta_i] << endl;
    Residual_eta_correction_histo_val1->SetBinError(eta_i+1,ey_validation_val1_[eta_i]);

    Residual_eta_correction_histo_res2->SetBinContent(eta_i+1,y_residual_res2_[eta_i]);
    cout << "x: " << x_eta_[eta_i] << " and y: "<< y_residual_res2_[eta_i] << " and error: " << ey_residual_res2_[eta_i] << endl;
    Residual_eta_correction_histo_res2->SetBinError(eta_i+1,ey_residual_res2_[eta_i]);
    Residual_eta_correction_histo_val2->SetBinContent(eta_i+1,y_validation_val2_[eta_i]);
    cout << "x: " << x_eta_[eta_i] << " and y: "<< y_validation_val2_[eta_i] << " and error: " << ey_validation_val2_[eta_i] << endl;
    Residual_eta_correction_histo_val2->SetBinError(eta_i+1,ey_validation_val2_[eta_i]);
  }




  //  std::vector <Double_t> SLOPE_eta_;
  //  std::vector <Double_t> eSLOPE_eta_;


  TGraphErrors* Residual_eta_correction_res1 = new TGraphErrors(n,&x_eta_[0],&y_residual_res1_[0],&ex_eta_[0],&ey_residual_res1_[0]);
  Residual_eta_correction_res1->SetTitle("Residual correction (res1 constants)");
  Residual_eta_correction_res1->SetName("Residual_eta_correction_res1");
  TGraphErrorsstyle(Residual_eta_correction_res1,line_styles_, colours_, markers_, 1,1,1,"bla");
  TGraphErrors* Residual_eta_correction_val1 = new TGraphErrors(n,&x_eta_[0],&y_validation_val1_[0],&ex_eta_[0],&ey_validation_val1_[0]);
  Residual_eta_correction_val1->SetTitle("Residual correction (on top of L2L3res)");
  Residual_eta_correction_val1->SetName("Residual_eta_correction_val1");
  TGraphErrorsstyle(Residual_eta_correction_val1,line_styles_, colours_, markers_, 0,0,0,"bla");
  TGraphErrors* kFSR_vs_eta_res1 = new TGraphErrors(n,&x_eta_[0],&y_kFSR_res1_[0],&ex_eta_[0],&ey_kFSR_res1_[0]);
  kFSR_vs_eta_res1->SetTitle("kFSR (res1)");
  kFSR_vs_eta_res1->SetName("kFSR_vs_eta_res1");
  TGraphErrorsstyle(kFSR_vs_eta_res1,line_styles_, colours_, markers_, 1,1,1,"bla");
  TGraphErrors* kFSR_vs_eta_val1 = new TGraphErrors(n,&x_eta_[0],&y_kFSR_val1_[0],&ex_eta_[0],&ey_kFSR_val1_[0]);
  kFSR_vs_eta_val1->SetTitle("kFSR (val1idation)");
  kFSR_vs_eta_val1->SetName("kFSR_vs_eta_val1");
  TGraphErrorsstyle(kFSR_vs_eta_val1,line_styles_, colours_, markers_, 0,0,0,"bla");
  TGraphErrors* Residual_eta_correction_res2 = new TGraphErrors(n,&x_eta_[0],&y_residual_res2_[0],&ex_eta_[0],&ey_residual_res2_[0]);
  Residual_eta_correction_res2->SetTitle("Residual correction (res2 constants)");
  Residual_eta_correction_res2->SetName("Residual_eta_correction_res2");
  TGraphErrorsstyle(Residual_eta_correction_res2,line_styles_, colours_, markers_, 1,1,1,"bla");
  TGraphErrors* Residual_eta_correction_val2 = new TGraphErrors(n,&x_eta_[0],&y_validation_val2_[0],&ex_eta_[0],&ey_validation_val2_[0]);
  Residual_eta_correction_val2->SetTitle("Residual correction (on top of L2L3res)");
  Residual_eta_correction_val2->SetName("Residual_eta_correction_val2");
  TGraphErrorsstyle(Residual_eta_correction_val2,line_styles_, colours_, markers_, 0,0,0,"bla");
  TGraphErrors* kFSR_vs_eta_res2 = new TGraphErrors(n,&x_eta_[0],&y_kFSR_res2_[0],&ex_eta_[0],&ey_kFSR_res2_[0]);
  kFSR_vs_eta_res2->SetTitle("kFSR (res2)");
  kFSR_vs_eta_res2->SetName("kFSR_vs_eta_res2");
  TGraphErrorsstyle(kFSR_vs_eta_res2,line_styles_, colours_, markers_, 1,1,1,"bla");
  TGraphErrors* kFSR_vs_eta_val2 = new TGraphErrors(n,&x_eta_[0],&y_kFSR_val2_[0],&ex_eta_[0],&ey_kFSR_val2_[0]);
  kFSR_vs_eta_val2->SetTitle("kFSR (val1idation)");
  kFSR_vs_eta_val2->SetName("kFSR_vs_eta_val2");
  TGraphErrorsstyle(kFSR_vs_eta_val2,line_styles_, colours_, markers_, 0,0,0,"bla");


  // TF1 *kFSR_fit = new TF1("kFSR_fit","[0]+[1]*x+[2]*(2*x*x-1)",kFSR_vs_Abseta_res1->GetXaxis()->GetXmin(),kFSR_vs_Abseta_res1->GetXaxis()->GetXmax()); //was used before...
   TF1 *kFSR_fit = new TF1("kFSR_fit","[0]+[1]*cosh(x)/(1+cosh(x)*[2])",kFSR_vs_Abseta_res1->GetXaxis()->GetXmin(),kFSR_vs_Abseta_res1->GetXaxis()->GetXmax()); //was used before...
  kFSR_fit->SetParameters(0.9,0.014,0.26);
  kFSR_fit->SetParName(0,"const");
  kFSR_fit->SetParName(1,"par1");
  kFSR_fit->SetParName(2,"par2");
  //[0] + cosh(eta)*([1] + cosh(eta)*[2])


  TFitResultPtr kFSR_fit_result = kFSR_vs_Abseta_histo_res1->Fit("kFSR_fit","S");
  kFSR_vs_Abseta_histo_res2->Fit("kFSR_fit");
  kFSR_vs_Abseta_res1->Fit("kFSR_fit");
  kFSR_vs_Abseta_res2->Fit("kFSR_fit");
  kFSR_vs_Abseta_val1->Fit("kFSR_fit");
  kFSR_vs_Abseta_val2->Fit("kFSR_fit");


//  kFSR_vs_Abseta_histo_res1->Fit("kFSR_fit","","",0.01,3.5);
//  kFSR_vs_Abseta_histo_res2->Fit("kFSR_fit","","",0.01,3.5);
//  kFSR_vs_Abseta_res1->Fit("kFSR_fit","","",0.01,3.5);
//  kFSR_vs_Abseta_res2->Fit("kFSR_fit","","",0.01,3.5);
//  kFSR_vs_Abseta_val1->Fit("kFSR_fit","","",0.01,3.5);
//  kFSR_vs_Abseta_val2->Fit("kFSR_fit","","",0.01,3.5);
//
//  kFSR_vs_Abseta_histo_res1->Fit("kFSR_fit","","",0.01,5.1);
//  kFSR_vs_Abseta_histo_res2->Fit("kFSR_fit","","",0.01,5.1);
//  kFSR_vs_Abseta_res1->Fit("kFSR_fit","","",0.01,5.1);
//  kFSR_vs_Abseta_res2->Fit("kFSR_fit","","",0.01,5.1);
//  kFSR_vs_Abseta_val1->Fit("kFSR_fit","","",0.01,5.1);
//  kFSR_vs_Abseta_val2->Fit("kFSR_fit","","",0.01,5.1);


  Double_t residual_calo[]={0.993637,1.0024,1.01285,1.02215,1.02115,1.02437,0.99169,0.995611,0.986447,0.941971,0.927278,0.971096,0.977086,0.929923};
  Double_t residual_pf[]={1.00014,1.00137,1.00926,1.0181,1.01773,1.03351,1.01248,1.00628,0.990574,0.936734,0.945384,1.02474,1.04382,1.02495};
  Double_t residual_jpt[]={1.1,1.1,1.1,1.1,1.1,1.1,1.1,1.1,1.1,1.1,1.1,1.1,1.1,1.1};
  Double_t e_residuals_existing[]={0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  
  for(unsigned int L3_i=0; L3_i<trad_x_Abseta_.size();L3_i++){
    residual_calo[L3_i]=residual_calo[L3_i]/1.00725;
    residual_pf[L3_i]=residual_pf[L3_i]/1.00741;
  }

  TGraphErrors* Residual_Abseta_correction_existing;
  if(jet_type=="PF")Residual_Abseta_correction_existing = new
    TGraphErrors(trad_x_Abseta_.size(),&trad_x_Abseta_[0],&residual_pf[0],&trad_ex_Abseta_[0],&e_residuals_existing[0]);
  if(jet_type=="PFCHS")Residual_Abseta_correction_existing = new
    TGraphErrors(trad_x_Abseta_.size(),&trad_x_Abseta_[0],&residual_pf[0],&trad_ex_Abseta_[0],&e_residuals_existing[0]);
  if(jet_type=="Calo")Residual_Abseta_correction_existing = new
    TGraphErrors(trad_x_Abseta_.size(),&trad_x_Abseta_[0],&residual_calo[0],&trad_ex_Abseta_[0],&e_residuals_existing[0]);
  if(jet_type=="JPT")Residual_Abseta_correction_existing = new
    TGraphErrors(trad_x_Abseta_.size(),&trad_x_Abseta_[0],&residual_jpt[0],&trad_ex_Abseta_[0],&e_residuals_existing[0]);



//vector<Double_t>::iterator p;
//
//for (p = trad_x_Abseta_.begin(); p != trad_x_Abseta_.end(); p++)
//   cout << *p << " ";
//  trad_x_Abseta_.Print();
  Residual_Abseta_correction_existing->SetTitle("Residual correction (existing)");
  Residual_Abseta_correction_existing->SetName("Residual_Abseta_correction_existing");
  TGraphErrorsstyle(Residual_Abseta_correction_existing,line_styles_, colours_, markers_, 2,2,2,"bla");

  std::cout<<"it should die right here"<<std::endl;
  TH1D* Residual_Abseta_correction_histo_existing=new TH1D("Residual_Abseta_correction_histo_existing", "Residual_Abseta_correction_histo_existing", trad_zero_eta , &trad_eta_binning[trad_zero_eta]);

      //      Residual_Abseta_correction_histo_existing->Print("all");
  for(Int_t Abseta_i=0;Abseta_i<trad_zero_eta;Abseta_i++){
    if(jet_type=="PF"){
    Residual_Abseta_correction_histo_existing->SetBinContent(Abseta_i+1,residual_pf[Abseta_i]);
    Residual_Abseta_correction_histo_existing->SetBinError(Abseta_i+1,e_residuals_existing[Abseta_i]);
    }
    else if(jet_type=="Calo"){
    Residual_Abseta_correction_histo_existing->SetBinContent(Abseta_i+1,residual_calo[Abseta_i]);
    Residual_Abseta_correction_histo_existing->SetBinError(Abseta_i+1,e_residuals_existing[Abseta_i]);
    }
  }
  TH1style(Residual_Abseta_correction_histo_existing,line_styles_, colours_, markers_, 2,2,2,"bla");

     if(chdir(fine_coarse+"_"+ generatorone+"_"+generatortwo+"_"+jet_type) != 0){ 
       mkdir(fine_coarse+"_"+ generatorone+"_"+generatortwo+"_"+jet_type, S_IRWXU|S_IRWXG|S_IRWXO); 
       chdir(fine_coarse+"_"+ generatorone+"_"+generatortwo+"_"+jet_type); 
     } 
  std::cout<<"sanity"<<std::endl;


  pm_eta_TGraph_draw(c, Residual_eta_correction_val1, Residual_Abseta_correction_val1, line_styles_, colours_, markers_);
  Residual_Abseta_correction_histo_val1->SetTitle("");
  pm_eta_TGraph_draw(c, Residual_eta_correction_val1, Residual_Abseta_correction_val1, Residual_Abseta_correction_histo_val1, line_styles_, colours_, markers_, "Correction Factor");

  Residual_Abseta_correction_res1->SetTitle("");
  pm_eta_TGraph_draw(c, Residual_eta_correction_res1, Residual_Abseta_correction_res1, line_styles_, colours_, markers_);
  Residual_Abseta_correction_histo_res1->SetTitle("");
  pm_eta_TGraph_draw(c, Residual_eta_correction_res1, Residual_Abseta_correction_res1, Residual_Abseta_correction_histo_res1, line_styles_, colours_, markers_, "Correction Factor");
  pm_eta_TGraph_draw(c, kFSR_vs_eta_res1, kFSR_vs_Abseta_res1, line_styles_, colours_, markers_);
  kFSR_vs_Abseta_histo_res1->SetTitle("");
  pm_eta_TGraph_draw(c, kFSR_vs_eta_res1, kFSR_vs_Abseta_res1, kFSR_vs_Abseta_histo_res1, line_styles_, colours_, markers_, "k_{FSR}",0.9,1.1);




  kFSR_vs_Abseta_histo_res1->Draw("P");
  kFSR_vs_Abseta_histo_res1->SetTitle("");
  kFSR_vs_Abseta_histo_res1->GetYaxis()->SetTitle("k_{FSR}");
  kFSR_vs_Abseta_histo_res1->GetYaxis()->SetRangeUser(0.9,1.1);
  c->Print("kFSR_res1_Abseta"+image_ext);
  //  TH1style_plus_fit(kFSR_vs_Abseta_histo_res2,line_styles_, colours_, markers_, 1,1,1,"bla","kFSR_fit");
  kFSR_vs_Abseta_histo_res2->Draw("P same");
  TLegend *leg_kFSR;
  leg_kFSR = new TLegend(0.25,0.20,0.55,0.40);
  leg_kFSR->SetFillColor(kWhite);
  //   leg->SetHeader("Legende");
  leg_kFSR->AddEntry(kFSR_vs_Abseta_histo_res1,generatorone,"lep");
  leg_kFSR->AddEntry(kFSR_vs_Abseta_histo_res2,generatortwo,"lep");
  leg_kFSR->Draw();
  c->Print("kFSR_res1_res2_Abseta"+image_ext);


  Residual_Abseta_correction_res1->Draw("AP");
  Residual_Abseta_correction_res1->GetXaxis()->SetTitle("|#eta|");
  Residual_Abseta_correction_res1->GetYaxis()->SetRangeUser(0.8,1.1);
  Residual_Abseta_correction_res1->GetYaxis()->SetTitle("Correction Factor");

  Residual_Abseta_correction_res2->Draw("P same");
  Residual_Abseta_correction_histo_existing->Draw("hist same");
  //  Residual_Abseta_correction_existing->Draw("LP same");

  TLegend *leg_selected;
  leg_selected = new TLegend(0.25,0.20,0.55,0.40);
  leg_selected->SetFillColor(kWhite);
  //   leg->SetHeader("Legende");
  leg_selected->AddEntry(Residual_Abseta_correction_res1,"("+ generatorone+ ")","lep");
  leg_selected->AddEntry(Residual_Abseta_correction_res2,"("+ generatortwo+ ")","lep");
  leg_selected->AddEntry(Residual_Abseta_correction_existing,"Residual correction (existing)","lep");
  leg_selected->Draw();
  c->Print(fine_coarse+"_"+ generatorone+"_"+generatortwo+"_"+jet_type+"Corr_factor_res1_res2"+image_ext);

  Residual_Abseta_correction_existing->Draw("ALP");
  c->Print(fine_coarse+"_"+ generatorone+"_"+generatortwo+"_"+jet_type+"Corr_factor_existing"+image_ext);


  Residual_Abseta_correction_val1->Draw("AP");
  Residual_Abseta_correction_val1->GetXaxis()->SetTitle("|#eta|");
  Residual_Abseta_correction_val1->GetYaxis()->SetRangeUser(0.8,1.1);
  Residual_Abseta_correction_val1->GetYaxis()->SetTitle("Correction Factor");

  TGraphErrorsstyle(Residual_Abseta_correction_val2,line_styles_, colours_, markers_, 1,1,1,"bla");
  Residual_Abseta_correction_val2->Draw("P same");
  Residual_Abseta_correction_histo_existing->Draw("hist same");
  //  Residual_Abseta_correction_existing->Draw("LP same");

  TLegend *leg_val;
  leg_val = new TLegend(0.25,0.20,0.55,0.40);
  leg_val->SetFillColor(kWhite);
  //   leg->SetHeader("Legende");
  leg_val->AddEntry(Residual_Abseta_correction_val1,"L2L3Res on top ("+ generatorone+ ")","lep");
  leg_val->AddEntry(Residual_Abseta_correction_val2,"L2L3Res on top ("+ generatortwo+ ")","lep");
  leg_val->AddEntry(Residual_Abseta_correction_existing,"Residual correction (existing)","lep");
  leg_val->Draw();
  c->Print(fine_coarse+"_"+ generatorone+"_"+generatortwo+"_"+jet_type+"Corr_factor_val1_val2"+image_ext);




  c->SetLogy();
MCone_eta_spectrum_MC->Draw("hist");
MCone_eta_spectrum_data->Draw("pe same");
c->Print("Eta_Spectrum_"+ generatorone+"_overlay"+image_ext);


MCtwo_eta_spectrum_MC->Draw("hist");
MCtwo_eta_spectrum_data->Draw("pe same");
c->Print("Eta_Spectrum_"+ generatortwo+"_overlay"+image_ext);


 c->SetLogx();
MCone_pt_spectrum_MC->Draw("hist");
MCone_pt_spectrum_data->Draw("pe same");
c->Print("Pt_Spectrum_"+ generatorone+"_overlay"+image_ext);


MCtwo_pt_spectrum_MC->Draw("hist");
MCtwo_pt_spectrum_data->Draw("pe same");
c->Print("Pt_Spectrum_"+ generatortwo+"_overlay"+image_ext);


  chdir("..");




  //////////////NPV
 Int_t Nbinsx = NPV_all_eta_MC1_L2L3_.at(0)->GetNbinsX();
 
 std::vector <TGraphErrors*> Residuals_Abseta_NPV_;
 std::vector <TGraphErrors*> Residuals_eta_NPV_;
 std::vector <TGraphErrors*> Residuals_norm_Abseta_NPV_;
 std::vector <TGraphErrors*> Residuals_norm_eta_NPV_;
 for(Int_t bins_i=1;bins_i<=Nbinsx;bins_i++){

   std::vector <Double_t> Abseta_res_y_;
   std::vector <Double_t> Abseta_norm_res_y_;
   std::vector <Double_t> Abseta_res_ey_;
   for(Int_t Abseta_i=0;Abseta_i<no_Abseta_bins;Abseta_i++){
     cout << "Abseta_i: " << Abseta_i << " bins_i: " << bins_i << " " <<NPV_all_Abseta_ratio_res1_.at(Abseta_i)->GetBinContent(bins_i) << endl;
     Abseta_res_y_.push_back(NPV_all_Abseta_ratio_res1_.at(Abseta_i)->GetBinContent(bins_i));
     Abseta_norm_res_y_.push_back(NPV_all_Abseta_ratio_res1_.at(Abseta_i)->GetBinContent(bins_i)/NPV_all_Abseta_ratio_res1_.at(Abseta_i)->GetBinContent(3));
     Abseta_res_ey_.push_back(NPV_all_Abseta_ratio_res1_.at(Abseta_i)->GetBinError(bins_i));
   }
   TGraphErrors* collect_Abseta_NPV = new TGraphErrors(no_Abseta_bins,&x_Abseta_[0],&Abseta_res_y_[0],&ex_Abseta_[0],&Abseta_res_ey_[0]);
   collect_Abseta_NPV->SetName((TString)"NPVAbseta_res_"+(Long_t)bins_i);
   collect_Abseta_NPV->SetTitle((TString)"NPVAbseta_res_"+(Long_t)bins_i+";|#eta|;MC/Data");
   Residuals_Abseta_NPV_.push_back(collect_Abseta_NPV);
   TGraphErrors* collect_norm_Abseta_NPV = new TGraphErrors(no_Abseta_bins,&x_Abseta_[0],&Abseta_norm_res_y_[0],&ex_Abseta_[0],&Abseta_res_ey_[0]);
   collect_norm_Abseta_NPV->SetName((TString)"NPVAbseta_norm_res_"+(Long_t)bins_i);
   collect_norm_Abseta_NPV->SetTitle((TString)"NPVAbseta_norm_res_"+(Long_t)bins_i+";|#eta|;deviation from 4-5");
   Residuals_norm_Abseta_NPV_.push_back(collect_norm_Abseta_NPV);

   std::vector <Double_t> eta_res_y_;
   std::vector <Double_t> eta_norm_res_y_;
   std::vector <Double_t> eta_res_ey_;
   for(Int_t eta_i=0;eta_i<no_eta_bins;eta_i++){
     cout << "eta_i: " << eta_i << " bins_i: " << bins_i << " " <<NPV_all_eta_ratio_res1_.at(eta_i)->GetBinContent(bins_i) << endl;
     eta_res_y_.push_back(NPV_all_eta_ratio_res1_.at(eta_i)->GetBinContent(bins_i));
     eta_norm_res_y_.push_back(NPV_all_eta_ratio_res1_.at(eta_i)->GetBinContent(bins_i)/NPV_all_eta_ratio_res1_.at(eta_i)->GetBinContent(3));
     eta_res_ey_.push_back(NPV_all_eta_ratio_res1_.at(eta_i)->GetBinError(bins_i));
   }
   TGraphErrors* collect_eta_NPV = new TGraphErrors(no_eta_bins,&x_eta_[0],&eta_res_y_[0],&ex_eta_[0],&eta_res_ey_[0]);
   collect_eta_NPV->SetName((TString)"NPVeta_res_"+(Long_t)bins_i);
   collect_eta_NPV->SetTitle((TString)"NPVeta_res_"+(Long_t)bins_i+";#eta;MC/Data");
   Residuals_eta_NPV_.push_back(collect_eta_NPV);
   TGraphErrors* collect_norm_eta_NPV = new TGraphErrors(no_eta_bins,&x_eta_[0],&eta_norm_res_y_[0],&ex_eta_[0],&eta_res_ey_[0]);
   collect_norm_eta_NPV->SetName((TString)"NPVeta_norm_res_"+(Long_t)bins_i);
   collect_norm_eta_NPV->SetTitle((TString)"NPVeta_norm_res_"+(Long_t)bins_i+";#eta;dev. from 4-5");
   Residuals_norm_eta_NPV_.push_back(collect_norm_eta_NPV);
 }

    TLegend *NPV_leg;
    NPV_leg = new TLegend(0.25,0.70,0.55,0.85);
    NPV_leg->SetFillColor(kWhite);
      //   leg->SetHeader("Legende");
    for(Int_t bins_i=0;bins_i<Nbinsx;bins_i++){
      NPV_leg->AddEntry(Residuals_eta_NPV_.at(bins_i),(TString)"NVtx-reg "+(Long_t)bins_i,"lep");
    }

    draw_Overlay_TGraphErrors(NPV_leg,Residuals_eta_NPV_, "Overlay_Eta_"+generatorone+"_Residuals_eta_NPV","nice", "ALPlegend", "x0_y0_z0", 0,-1, 0.97, 1.17, 0,-1, image_ext, ".");//,1.);
    draw_Overlay_TGraphErrors(NPV_leg,Residuals_Abseta_NPV_, "Overlay_Abseta_"+generatorone+"_Residuals_Abseta_NPV","nice", "ALPlegend", "x0_y0_z0", 0,-1, 0.97, 1.17, 0,-1, image_ext, ".");//,1.);
    draw_Overlay_TGraphErrors(NPV_leg,Residuals_norm_eta_NPV_, "Overlay_Eta_"+generatorone+"_Residuals_norm_eta_NPV","nice", "ALPlegend", "x0_y0_z0", 0,-1, 0.9, 1.1, 0,-1, image_ext, ".");//,1.);
    draw_Overlay_TGraphErrors(NPV_leg,Residuals_norm_Abseta_NPV_, "Overlay_Abseta_"+generatorone+"_Residuals_norm_Abseta_NPV","nice", "ALPlegend", "x0_y0_z0", 0,-1, 0.9, 1.1, 0,-1, image_ext, ".");//,1.);





  //draw_Overlay_TGraphErrors(TLegend *leg,std::vector < TGraphErrors* > histos_, TString PS_name ="DEFAULT_PS", TString format = "six", TString draw_options ="", TString Logxyz="x0_y0_z0", Double_t ru_xlow=0, Double_t ru_xhig=-1, Double_t ru_ylow=0, Double_t ru_yhig=-1, Double_t ru_zlow=0, Double_t ru_zhig=-1, TString img_exp="ENTER_WITH_POINT", TString SaveTo="PS", Double_t h_line=-55.)
 ////////////NPV


  if(root_export.Contains("true")){


    TFile *outf = new TFile(fine_coarse+"_"+ generatorone+"_"+generatortwo+"_"+jet_type+"_kFSR_histos.root","RECREATE");
    
    for(Int_t bins_i=0;bins_i<Nbinsx;bins_i++){
      Residuals_Abseta_NPV_.at(bins_i)->Write();
      Residuals_eta_NPV_.at(bins_i)->Write();
    }
    kFSR_vs_eta_val1   ->Write();
    kFSR_vs_eta_val2   ->Write();
    kFSR_vs_eta_res1   ->Write();
    kFSR_vs_eta_res2   ->Write();
    kFSR_vs_Abseta_val1->Write();
    kFSR_vs_Abseta_val2->Write();
    kFSR_vs_Abseta_res1->Write();
    kFSR_vs_Abseta_res2->Write();
    kFSR_vs_Abseta_histo_res1->Write();
    kFSR_vs_Abseta_histo_res2->Write();
    kFSR_fit_result->Write();   
  for(unsigned int cut_i=0;cut_i<ptthreecuts.size();cut_i++){
    if(ptthreecuts[cut_i]=="20"){
    for(Int_t Abseta_i=0;Abseta_i<no_Abseta_bins;Abseta_i++){

      all_ptthree_all_Abseta_ratio_res1_.at(cut_i).at(Abseta_i)->Write();
      all_ptthree_all_Abseta_ratio_val1_.at(cut_i).at(Abseta_i)->Write();
    }}
  }




    outf->Close();
  }
  
  TFile *collectf = new TFile("res_corrections_histos.root","UPDATE");

  Residual_slope_eta_correction_histo_res1->SetName(fine_coarse+"_"+ generatorone+"_"+jet_type+"_eta_res_slope_hist");
  Residual_slope_eta_correction_histo_res1->Write();
  Residual_slope_Abseta_correction_histo_res1->SetName(fine_coarse+"_"+ generatorone+"_"+jet_type+"_Abseta_res_slope_hist");
  Residual_slope_Abseta_correction_histo_res1->Write();

  Residual_const_slope_eta_correction_histo_res1->SetName(fine_coarse+"_"+ generatorone+"_"+jet_type+"_eta_res_const_slope_hist");
  Residual_const_slope_eta_correction_histo_res1->Write();
  Residual_const_slope_Abseta_correction_histo_res1->SetName(fine_coarse+"_"+ generatorone+"_"+jet_type+"_Abseta_res_const_slope_hist");
  Residual_const_slope_Abseta_correction_histo_res1->Write();

  Residual_eta_correction_histo_res1->SetName(fine_coarse+"_"+ generatorone+"_"+jet_type+"_eta_res_hist");
  Residual_eta_correction_histo_res1->Write();
  Residual_eta_correction_histo_res2->SetName(fine_coarse+"_"+ generatortwo+"_"+jet_type+"_eta_res_hist");
  Residual_eta_correction_histo_res2->Write();

  Residual_eta_correction_histo_val1->SetName(fine_coarse+"_"+ generatorone+"_"+jet_type+"_eta_val_hist");
  Residual_eta_correction_histo_val1->Write();
  Residual_eta_correction_histo_val2->SetName(fine_coarse+"_"+ generatortwo+"_"+jet_type+"_eta_val_hist");
  Residual_eta_correction_histo_val2->Write();

  Residual_Abseta_correction_histo_res1->SetName(fine_coarse+"_"+ generatorone+"_"+jet_type+"_Abseta_res_hist");
  Residual_Abseta_correction_histo_res1->Write();
  Residual_Abseta_correction_histo_res2->SetName(fine_coarse+"_"+ generatortwo+"_"+jet_type+"_Abseta_res_hist");
  Residual_Abseta_correction_histo_res2->Write();

  Residual_Abseta_correction_histo_val1->SetName(fine_coarse+"_"+ generatorone+"_"+jet_type+"_Abseta_val_hist");
  Residual_Abseta_correction_histo_val1->Write();
  Residual_Abseta_correction_histo_val2->SetName(fine_coarse+"_"+ generatortwo+"_"+jet_type+"_Abseta_val_hist");
  Residual_Abseta_correction_histo_val2->Write();

  collectf->Close();

   ofstream myfile;
  myfile.open (fine_coarse+"_"+ generatorone+"_L2L3Residual_AK5"+jet_type+".txt");
  myfile << "{ 1 JetEta 1 JetPt [0] Correction L2Relative}\n";
  //  Double_t L3_abs_offset=1.010;
  //Double_t L3_abs_offset=1.025;
  Double_t L3_abs_offset=1.014;
  //  myfile << L3_abs_offset  << "\n";

  for(Int_t eta_i=0;eta_i<no_eta;eta_i++){

    myfile << setw(15) << eta_binning[eta_i] << "       " << setw(10)  << eta_binning[eta_i+1] << "    3     3     3500        "  << setw(10) << Residual_eta_correction_histo_res1->GetBinContent(eta_i+1)*L3_abs_offset << "\n";


    //      TH1D* Residual_eta_correction_histo_res1=new TH1D("Residual_eta_correction_histo_res1", "Residual_eta_correction_histo_res1", no_eta , &eta_binning[0]);

    //eta_binning[]={-6.0,-4.0,-3.5,-3.2,-3.0,-2.8,-2.6,-2.4,-2.2,-2.0,-1.8,-1.5,-1.4,-1.3,-1.2,-1.1,-0.9,-0.6,-0.3,0.0,0.

    //  -5.191         -3.489              3              3           3500       0.997653


  }
  myfile.close();


  myfile.open (fine_coarse+"_"+ generatorone+"_L2L3Residual_AK5"+jet_type+"_PTDEP.txt");
  myfile << "{ 1 JetEta 1 JetPt [0]*([1]+[2]*TMath::Log(x)) Correction L2Relative}\n";
  //  myfile << L3_abs_offset  << "\n";

  Int_t cut_20=20;
    for(unsigned int cut_j=0;cut_j<ptthreecuts.size();cut_j++){
      if(ptthreecuts[cut_j]=="20"){
	cut_20=cut_j;}
    }
    cout << cut_20 << endl;

  for(Int_t eta_i=0;eta_i<no_eta;eta_i++){

    Int_t Abseta_i=TMath::Abs(TMath::Abs(eta_i-zero_eta+0.75));
    cout << Abseta_i << endl;
    Double_t kFSR_value = kFSR_vs_Abseta_histo_res1->GetBinContent(kFSR_vs_Abseta_histo_res1->FindBin( std::abs( (eta_binning[eta_i]+ eta_binning[eta_i+1])/2 )));
    
    myfile << setw(15) << eta_binning[eta_i] << "       " << setw(10)  << eta_binning[eta_i+1] << "    5         "  << setw(10)  <<  all_ptthree_all_Abseta_res1_ptreach_.at(cut_20).at(Abseta_i).first << setw(15) <<all_ptthree_all_Abseta_res1_ptreach_.at(cut_20).at(Abseta_i).second << setw(15) << kFSR_value*L3_abs_offset << setw(15) << Residual_const_slope_eta_correction_histo_res1->GetBinContent(eta_i+1) << setw(15) << Residual_slope_eta_correction_histo_res1->GetBinContent(eta_i+1) << "\n";

    //      TH1D* Residual_eta_correction_histo_res1=new TH1D("Residual_eta_correction_histo_res1", "Residual_eta_correction_histo_res1", no_eta , &eta_binning[0]);

    //eta_binning[]={-6.0,-4.0,-3.5,-3.2,-3.0,-2.8,-2.6,-2.4,-2.2,-2.0,-1.8,-1.5,-1.4,-1.3,-1.2,-1.1,-0.9,-0.6,-0.3,0.0,0.

    //  -5.191         -3.489              3              3           3500       0.997653


  }
  myfile.close();
 
  myfile.open (fine_coarse+"_"+ generatorone+"_Abseta_L2L3Residual_AK5"+jet_type+".txt");
  myfile << "{ 1 JetEta 1 JetPt [0] Correction L2Relative}\n";

  for(Int_t eta_i=0;eta_i<no_eta;eta_i++){
    Int_t Abseta_i=TMath::Abs(TMath::Abs(eta_i-zero_eta+0.75));
    myfile << setw(15) << eta_binning[eta_i] << "       " << setw(10)  << eta_binning[eta_i+1] << "    3     3     3500        "  << setw(10) << Residual_Abseta_correction_histo_res1->GetBinContent(kFSR_vs_Abseta_histo_res1->FindBin( std::abs( (eta_binning[eta_i]+ eta_binning[eta_i+1])/2 )))*L3_abs_offset << "\n";


  }
  myfile.close();




 chdir(".."); 

 if(no_entries_went_wrong)cout << "DA GING ETWAS SCHIEF (Z.B. PF jets bei -3 in eta) und zwar bei:" << endl;

 for(unsigned int i=0;i<faulty_values_dummy_histos.size();i++)std::cout << faulty_values_dummy_histos.at(i) << endl;








}


