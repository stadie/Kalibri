#ifndef flex_param_h
#define flex_param_h

#include "THelpers.h"
#include <TMatrixDSym.h>
#include <TDecompSVD.h>
//#include "Fit/FitResult.h"
//#include "Fit/FitResultPtr.h"
#include "TFitResult.h"

//#include <TFitResultPtr.h>


class flex_param : public base_corr{

public :

  const static Int_t dof_in_fit_ = 5; //no. of fit parameters

  std::vector < std::vector < TH1D* > > tlj_X_response_all_corrections_all_;
  std::vector < std::vector < TH1D* > > tlj_Corr_Vars_counts_all_;
  std::vector < std::vector < TH2D* > > tlj_X_response_2D_all_;
  std::vector < std::vector < TProfile* > > tlj_X_response_prof_all_;
  std::vector < std::vector < TH1D* > > tlj_X_response_GMP_mean_all_;


  std::vector < std::vector < TGraphErrors* > > fit_params_in_X_all_;

  std::vector < TGraphErrors* > Correlation_B_C_;
  std::vector < std::vector < std::vector < Double_t > > > Corr_coeffs_all_;
  std::vector <  TString > correlation_labels_;

  Int_t offset_for_X;

  virtual void     Loop(Bool_t setvalues=0, Int_t par_bin_choice=0, Int_t par_X_choice=0, Int_t par_eta_choice=0, TString par_Corr_choice="12", TString par_img_choice=".pdf");
  virtual void Import_Histos();
  virtual void Book_Histos();
  virtual void Fit_Histos();
  virtual void Create_TGraphErrors();
  virtual void Write_Histos();
  virtual void Write_TGraphErrors();


};

#endif
