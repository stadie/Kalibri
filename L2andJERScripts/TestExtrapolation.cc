#include "Extrapolation.cc"

void TestExtrapolation(){
  //2011 resolutions with smearing
  //    Extrapolation test("ResolutionVsPt","/afs/naf.desy.de/user/k/kirschen/public/2011withJERcorr_KalibriPlots.root");
  //  test.Plot();


  //2011 resolutions without JER smearing
  //    Extrapolation test("ResolutionVsPt","/afs/naf.desy.de/user/k/kirschen/public/2011_KalibriPlots.root");
  //  test.Plot();

  //2012 resolutions
  Extrapolation test("ResolutionVsPt","/afs/naf.desy.de/user/k/kirschen/public/2012_KalibriPlots.root");
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
