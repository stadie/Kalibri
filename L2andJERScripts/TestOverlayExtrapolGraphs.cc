//#include "TimeDependence.cc"
#include "OverlayExtrapolGraphs.cc"

void TestOverlayExtrapolGraphs(){

  std::vector<TString> samplelist;

  samplelist.push_back("44Z2");
  samplelist.push_back("44Hpp");
  samplelist.push_back("42Z2");

  for(size_t i=0;i<samplelist.size();i++){
//    Extrapolation test1b("RelResponseVsRunNumber",samplelist.at(i));
//    test1b.Plot();

    OverlayExtrapolGraphs test1a;
    test1a.addPlots("OneBinRelResponseVsPt",samplelist.at(i));
    //    test1a.addPlots("OneBinMPFResponseVsPt",samplelist.at(i));
    test1a.plotAll();

    OverlayExtrapolGraphs test2;
    test2.addPlots("MPFResponseVsMeanPt",samplelist.at(i));
    test2.addPlots("RelResponseVsPt",samplelist.at(i));
    test2.plotAll();
  }

  OverlayExtrapolGraphs test1b;
  test1b.addPlots("OneBinRelResponseVsPt","44Z2");
  test1b.addPlots("OneBinMPFResponseVsPt","44Z2");
  test1b.addPlots("OneBinRelResponseVsPt","44Hpp");
  test1b.addPlots("OneBinMPFResponseVsPt","44Hpp");
  test1b.plotAll();


};


# ifndef __CINT__
int main()
{
  TestOverlayExtrapolGraphs();

  return 0;
}
# endif
