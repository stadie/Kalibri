#include "flex_corr.C"

void run_corr(Bool_t setvalues=0, Int_t par_bin_choice=0, Int_t par_eta_choice=0, TString par_img_choice=".pdf"){


   flex_corr t;
t.Loop(setvalues, par_bin_choice, par_eta_choice, par_img_choice);
}
