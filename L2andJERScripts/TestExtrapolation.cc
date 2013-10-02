#include "Extrapolation.cc"
#include "TMinuit.h"

void TestExtrapolation(){

  Extrapolation test("ResolutionVsPt","Resolutions2012");
  test.Plot();
  test.ExportTables();

};


# ifndef __CINT__
int main()
{
  TestExtrapolation();

  return 0;
}
# endif
