#include "do_flex_extrapol.C"

void run_do_flex_extrapol(TString jet_type, TString generatorone, TString
generatortwo, TString image_ext, TString root_export, TString
use_imported_kFSRAbs, TString fine_coarse, TString use_easy_mean,
TString use_fitted_kFSR, TString corr_generation="Spring10", TString
		      ratio_of_mean_or_GM="Means", Bool_t
			  export_all_plots=true, TString kFSR_eq_one ="kFSR_no_eq_one", TString MPF_or_rel_response ="rel_response"){


   do_flex_extrapol t;
t.Loop(jet_type, generatorone,generatortwo,image_ext, root_export,use_imported_kFSRAbs, fine_coarse, use_easy_mean,use_fitted_kFSR,corr_generation, ratio_of_mean_or_GM, export_all_plots, kFSR_eq_one, MPF_or_rel_response);
}
