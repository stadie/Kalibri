#include "flex_flav_diff.C"

void run_flav_diff(Bool_t setvalues=0, Int_t par_bin_choice=0, Int_t par_X_choice=0, TString par_eta_choice="0", TString par_Corr_choice="12", TString par_img_choice=".pdf", TString par_PDG_choice="01"){


   flex_flav_diff t;
t.Loop(setvalues, par_bin_choice, par_X_choice, par_eta_choice, par_Corr_choice, par_img_choice, par_PDG_choice);
}
