#include "PFFractionPlots.cc"

void TestPFPlots(){
  //  TString pathToRootFile ="../2012TEST_CORRFinal2011_AK5_MC_Su12Z2Star_kostas_TEMP_MOREPLOTS_TEMP_TEMP_AK5/dijetsFall10_TuneZ2_AK5PF_weighted_residuals_kostas/plots/KalibriPlots.root";
  //hffix - needs to be set in externalconfigs as well
  TString pathToRootFile ="../2012TEST_CORRFinal2011_AK5_MC_Su12Z2Star_kostas_TEMP_MOREPLOTS_TEMP_TEMP_AK5/dijetsFall10_TuneZ2_AK5PF_weighted_residuals_k_HFfix/plots/KalibriPlots.root";

  PFFractionPlots test("AbsPFFractionVsPt",pathToRootFile);
  test.Plot();

  //  PFFractionPlots test2("AbsPFFractionVsPt");
  //  test2.Plot();

//  PFFractionPlots test2("PFFractionVsEta");
//  test2.Plot();
//

  PFFractionPlots test3("AbsPFFractionVsEta",pathToRootFile);
  test3.Plot();

  PFFractionPlots test4("AbsPFFractionVsPt_alleta",pathToRootFile);
  test4.Plot();

  PFFractionPlots test5("AbsPFFractionVsEta_allpt",pathToRootFile);
  test5.Plot();

  //  test.setBinning("kostas",true);
  //  test.setBinning("kostas",true);
  //  test.print();
};


# ifndef __CINT__
int main()
{
  TestPFPlots();

  return 0;
}
# endif
