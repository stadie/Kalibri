#ifndef flex_corr_h
#define flex_corr_h

#include "THelpers.h"

class flex_corr : public base_corr{

public :

std::vector < std::vector < TH1D* > > tlj_Corr_Vars_counts_all_;

std::vector < std::vector < TH2D* > > tlj_X_response_2D_all_;
std::vector < std::vector < TProfile* > > tlj_X_response_prof_all_;
std::vector < std::vector < TH1D* > > tlj_X_response_GMP_mean_all_;
std::vector < std::vector < TH1D* > > tlj_X_response_GMP_sigma_all_;
std::vector < std::vector < TH1D* > > tlj_X_response_GMP_rel_sigma_all_;
std::vector < std::vector < TH2D* > > tlj_X_raw_response_2D_all_;
std::vector < std::vector < TProfile* > > tlj_X_raw_response_prof_all_;
std::vector < std::vector < TH1D* > > tlj_X_raw_response_GMP_mean_all_;
std::vector < std::vector < TH2D* > > tlj_X_L2_response_2D_all_;
std::vector < std::vector < TProfile* > > tlj_X_L2_response_prof_all_;
std::vector < std::vector < TH1D* > > tlj_X_L2_response_GMP_mean_all_;



  virtual void     Loop(Bool_t setvalues=0, Int_t par_bin_choice=0, Int_t par_eta_choice=0, TString par_img_choice=".pdf");
 virtual void Book_Histos();
 virtual void Write_Histos();

};

#endif
