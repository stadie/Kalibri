#include "Extrapolation.cc"

void TestExtrapolation(){
  //2011 resolutions with smearing
  //    Extrapolation test("ResolutionVsPt","/afs/naf.desy.de/user/k/kirschen/public/2011withJERcorr_kostas_KalibriPlots.root");
  //  test.Plot();


  //2011 resolutions without JER smearing
  //    Extrapolation test("ResolutionVsPt","/afs/naf.desy.de/user/k/kirschen/public/2011_kostas_KalibriPlots.root");
  //  test.Plot();

  //2012 resolutions
  Extrapolation test("ResolutionVsPt","/afs/naf.desy.de/user/k/kriheine/scratch/Kalibri/2011Full2011_CORRFinal2011_AK5_MC_F1144Z2wPU_kostas_DefaultTestTEst/dijetsFall10_TuneZ2_AK5PF_weighted_residuals_JER_common/plots/KalibriPlots.root");
  test.Plot();
  test.ExportTables();


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
