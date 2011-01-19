#ifndef flex_res_h
#define flex_res_h

#include "THelpers.h"
#include <TMatrixDSym.h>
#include <TDecompSVD.h>

class flex_res : public base_corr{

public :
  //  std::vector < std::vector < TH1D* > > tlj_X_counts_all_;

  Double_t  PDG_id;

  
  TF1* Fit_B_all;
  TF1* Fit_C_all;     
  TF1* Fit_X0_all;    
  TF1* Fit_Sigma_all; 
  TF1* get_correction;


  std::vector < std::vector < TH2D* > > tlj_X_response_2D_all_;
  std::vector < std::vector < TProfile* > > tlj_X_response_prof_all_;
  std::vector < std::vector < TH1D* > > tlj_X_response_GMP_mean_all_;
  std::vector < std::vector < TH1D* > > tlj_X_response_all_corrections_all_;


  ////RESOLUTION PART
  //  std::vector < TMatrixDSym > covariance_helper_matrices;
	const static unsigned int n = 2;	// number of variables
  std::vector< std::vector < std::vector < TH2D* > > >tlj_Sel_Correlations_2D_counts_all_;
  std::vector< std::vector < std::vector < TH2D* > > >tlj_Sel_Correlations_2D_counts_helper_all_;
  std::vector< std::vector < std::vector < TH2D* > > >tlj_Sel_Correlations_2D_response_all_;
  std::vector< std::vector < std::vector < TH2D* > > >tlj_Sel_Correlations_2D_response2_all_;
  std::vector< std::vector < std::vector < TH2D* > > >tlj_Sel_Correlations_2D_mean_response_all_;

  std::vector < std::vector < TH1D* > > tlj_corrected_response_barrel_all_;
  std::vector < TGraphErrors* > tlj_Response_Graphs_mean_all_;
  std::vector < TGraphErrors* > tlj_Response_Graphs_rel_sigma_all_;
  std::vector < TGraphErrors* > tlj_Response_Graphs_rel_resol_all_;

  std::vector < TGraphErrors* > tlj_Response_Graphs_mean_selec_;
  std::vector < TGraphErrors* > tlj_Response_Graphs_rel_sigma_selec_;
  std::vector < TGraphErrors* > tlj_Response_Graphs_rel_resol_selec_;

  TString root_resol_name_binning;


  virtual void Loop(Bool_t setvalues=0, Int_t par_bin_choice=0, Int_t par_X_choice=0, Int_t par_eta_choice=0, TString par_Corr_choice="12", TString par_img_choice=".pdf", Int_t par_PDG_choice=5);
  virtual void Import_Histos();
  virtual void Book_Histos();
  virtual void Write_Histos();
  virtual void Write_TGraphErrors();
  Double_t parametrized_correction(Double_t pt, Double_t sigma_phi);



};

#endif
