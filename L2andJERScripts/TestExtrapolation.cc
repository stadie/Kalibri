#include "Extrapolation.cc"

void TestExtrapolation(){
  //2011 resolutions
  //  Extrapolation test("ResolutionVsPt","/afs/naf.desy.de/user/k/kirschen/scratch/2012_04_L2L3ResidualsSummary/2011Full2011_CORRF11DB_He_AK5_MC_F11Z2wPUsm_Y_f_kostas_AK5/dijetsFall10_TuneZ2_AK5PF_weighted_residuals_kostas/plots/KalibriPlots.root");
  //  test.Plot();

  //2012 resolutions
  Extrapolation test("ResolutionVsPt","/afs/naf.desy.de/user/k/kirschen/scratch/2012_04_L2L3ResidualsSummary/2012TEST_CORRFinal2011_AK5_MC_Su12Z2Star_kostas_TrueReweighting_AK5/dijetsFall10_TuneZ2_AK5PF_weighted_residuals_kostas/plots/KalibriPlots.root");
  test.Plot();



  //plots for pt-dependent krad-correction and resulting residuals
//  Extrapolation test2("RelResponseVsPt");
//  test2.Plot();
//  
//  Extrapolation test3("MPFResponseVsPt");
//  test3.Plot();
  


};


# ifndef __CINT__
int main()
{
  TestExtrapolation();

  return 0;
}
# endif
