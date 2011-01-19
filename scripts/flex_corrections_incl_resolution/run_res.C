#include "flex_res.C"

void run_res(Bool_t setvalues=0, Int_t par_bin_choice=0, Int_t par_X_choice=0, Int_t par_eta_choice=0, TString par_Corr_choice="12", TString par_img_choice=".pdf", Int_t par_PDG_choice=4){


   flex_res t;
t.Loop(setvalues, par_bin_choice, par_X_choice, par_eta_choice, par_Corr_choice, par_img_choice, par_PDG_choice);
}
