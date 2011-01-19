#include "flex_param.C"

void run_param(Bool_t setvalues=0, Int_t par_bin_choice=0, Int_t par_X_choice=0, Int_t par_eta_choice=0, TString par_Corr_choice="12", TString par_img_choice=".pdf"){


   flex_param t;
t.Loop(setvalues, par_bin_choice, par_X_choice, par_eta_choice, par_Corr_choice, par_img_choice);
}
