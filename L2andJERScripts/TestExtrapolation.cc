#include "Extrapolation.cc"

void TestExtrapolation(){

  Extrapolation test("ResolutionVsPt","Resolutions44Z2");
  test.Plot();


  std::vector<TString> samplelist;
  samplelist.push_back("42Z2");
  samplelist.push_back("44Z2");
  samplelist.push_back("44Hpp");


  for(size_t i=0;i<samplelist.size();i++){
    Extrapolation test1("MPFResponseVsRunNumber",samplelist.at(i));
    test1.Plot();

    Extrapolation test2("MPFResponseVsJetWidth",samplelist.at(i));
    test2.Plot();
    
    Extrapolation test3("MPFResponseVsPF_CH_Fraction",samplelist.at(i));
    test3.Plot();

    Extrapolation test4("MPFResponseVsJetWidth",samplelist.at(i));
    test4.Plot();
    
    Extrapolation test5("MPFResponseVsNPV",samplelist.at(i));
    test5.Plot();
    
    Extrapolation test6("RelResponseVsNPV",samplelist.at(i));
    test6.Plot();
    
    Extrapolation test7("MPFResponseVsPt",samplelist.at(i));
    test7.Plot();
    
    Extrapolation test8("RelResponseVsPt",samplelist.at(i));
    test8.Plot();
    
    Extrapolation test9("OneBinRelResponseVsPt",samplelist.at(i));
    test9.Plot();
    
//    Extrapolation test10("OneBinMPFVsPt",samplelist.at(i));
//    test10.Plot();
  }



};


# ifndef __CINT__
int main()
{
  TestExtrapolation();

  return 0;
}
# endif
