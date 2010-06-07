#ifndef flex_param_h
#define flex_param_h

#include "THelpers.h"
#include <TMatrixDSym.h>
#include <TDecompSVD.h>

class flex_param : public base_corr{

public :

  const static Int_t dof_in_fit_ = 4; //no. of fit parameters

  std::vector < std::vector < TH1D* > > tlj_X_response_all_corrections_all_;
  std::vector < std::vector < TH1D* > > tlj_Corr_Vars_counts_all_;
  std::vector < std::vector < TH2D* > > tlj_X_response_2D_all_;
  std::vector < std::vector < TProfile* > > tlj_X_response_prof_all_;
  std::vector < std::vector < TH1D* > > tlj_X_response_GMP_mean_all_;

  std::vector < std::vector < TF1* > > fit_functions_all_;

  std::vector < std::vector < TGraphErrors* > > fit_params_in_X_all_;



  virtual void     Loop();
  virtual void Import_Histos();
  virtual void Book_Histos();
  virtual void Fit_Histos();
  virtual void Create_TGraphErrors();
  virtual void Write_Histos();
  virtual void Write_TGraphErrors();


};

#endif
