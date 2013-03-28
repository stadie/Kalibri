//#include "TimeDependence.cc"
#include "OverlayPlainProfileAnd2DPlots.cc"

void TestOverlayPlainProfileAnd2DPlots(){


  std::vector<TString> flavor_list;
  flavor_list.push_back("Unmatched");
  flavor_list.push_back("Gluon");
  flavor_list.push_back("uds");
  flavor_list.push_back("Charm");
  flavor_list.push_back("Bottom");

  assert(flavor_list.size()==5);//needs to be consistent with binning in ControlPlots

  std::vector<TString> InputTag_list;
  InputTag_list.push_back("Z2star");
  InputTag_list.push_back("Z2");
  InputTag_list.push_back("Herwig");

  std::vector<TString> samplelist;
  samplelist.push_back("PFCHSwithNu(phys)");
  samplelist.push_back("PFCHS(algo)");
  samplelist.push_back("PFCHS(phys)");
  samplelist.push_back("Calo(phys)");

  for(size_t i=0;i<samplelist.size();i++){

	//comparing all flavor responses of Herwig to Z2star 
    OverlayPlainProfileAnd2DPlots test10a("|#eta|<1.3; " + samplelist.at(i)+"; (Herwig/Z2star)","Flavor Response (Herwig/Z2star)");
    test10a.ConfigureSetRangeUser(0.95,1.15);
    for(size_t f_i=0;f_i<flavor_list.size();f_i++){
      test10a.addPlots("FlavorResponse",samplelist.at(i),"MCTruthRespFlavorVsGenJetPt","Z2star",f_i,flavor_list.at(f_i));
      test10a.addPlots("FlavorResponse",samplelist.at(i),"MCTruthRespFlavorVsGenJetPt","Herwig",f_i,flavor_list.at(f_i));
    }
    test10a.plotNormalizedComboRatiosAll();


//
//    //comparing the normalized flavor response for each MC generator (e.g. b/all for Herwig, Z2, Z2star)
//    OverlayPlainProfileAnd2DPlots test9a("|#eta|<1.3; " + samplelist.at(i)+"","Herwig++/PYTHIA Z2star (Response)");
//    for(size_t f_i=0;f_i<flavor_list.size();f_i++){
//      //	for(size_t i_i=0;i_i<InputTag_list.size();i_i++){
//      test9a.addPlots("DoubleRatioFlavorResponse",samplelist.at(i),"MCTruthResponseVsGenJetPt","Z2star",2,flavor_list.at(f_i));
//      test9a.addPlots("DoubleRatioFlavorResponse",samplelist.at(i),"MCTruthRespFlavorVsGenJetPt","Z2star",f_i,flavor_list.at(f_i));
//      test9a.addPlots("DoubleRatioFlavorResponse",samplelist.at(i),"MCTruthResponseVsGenJetPt","Herwig",2,flavor_list.at(f_i));
//      test9a.addPlots("DoubleRatioFlavorResponse",samplelist.at(i),"MCTruthRespFlavorVsGenJetPt","Herwig",f_i,flavor_list.at(f_i));
//      //	}
//    }
//    test9a.plotNormalizedDoubleComboRatiosAll();
//
//
//
//
//
//      for(size_t f_i=0;f_i<flavor_list.size();f_i++){
//	//comparing the normalized flavor response for each MC generator (e.g. b/all for Herwig, Z2, Z2star)
//	OverlayPlainProfileAnd2DPlots test8a("|#eta|<1.3; " + samplelist.at(i)+"; "+flavor_list.at(f_i)+"(norm.)","Flavor Response/All Response");
//	for(size_t i_i=0;i_i<InputTag_list.size();i_i++){
//	  test8a.addPlots("FlavorResponse",samplelist.at(i),"MCTruthResponseVsGenJetPt",InputTag_list.at(i_i),2,InputTag_list.at(i_i));
//	  test8a.addPlots("FlavorResponse",samplelist.at(i),"MCTruthRespFlavorVsGenJetPt",InputTag_list.at(i_i),f_i,flavor_list.at(f_i));
//	}
//	test8a.plotNormalizedComboRatiosAll();
//      }
//
//
//
//      //plot all flavor response for different generators
//      OverlayPlainProfileAnd2DPlots test7a("|#eta|<1.3; " + samplelist.at(i)+"; All flavor response","Mean Response");
//      for(size_t i_i=0;i_i<InputTag_list.size();i_i++)      test7a.addPlots("FlavorResponse",samplelist.at(i),"MCTruthResponseVsGenJetPt",InputTag_list.at(i_i),2,InputTag_list.at(i_i));
//      test7a.plotAll();
//
//
//      //plot all flavor response and response of each flavor for different generators
//      for(size_t i_i=0;i_i<InputTag_list.size();i_i++){
//	OverlayPlainProfileAnd2DPlots test5a("|#eta|<1.3; " + samplelist.at(i)+"; "+InputTag_list.at(i_i)+"(All)","Mean Response");
//	test5a.addPlots("FlavorResponse",samplelist.at(i),"MCTruthResponseVsGenJetPt",InputTag_list.at(i_i),2,"All Flavor");
//	for(size_t f_i=0;f_i<flavor_list.size();f_i++){
//	  test5a.addPlots("FlavorResponse",samplelist.at(i),"MCTruthRespFlavorVsGenJetPt",InputTag_list.at(i_i),f_i,flavor_list.at(f_i));
//	}
//	test5a.plotAll();
//      }
//
//      //plot response of each generator for different flavors
//      for(size_t f_i=0;f_i<flavor_list.size();f_i++){
//	OverlayPlainProfileAnd2DPlots test6a("|#eta|<1.3; " + samplelist.at(i)+"; "+flavor_list.at(f_i)+"(All)","Mean Response");
//	test6a.addPlots("FlavorResponse",samplelist.at(i),"MCTruthResponseVsGenJetPt",InputTag_list.at(0),2,"All Flavor("+InputTag_list.at(0) +")");
//	for(size_t i_i=0;i_i<InputTag_list.size();i_i++)test6a.addPlots("FlavorResponse",samplelist.at(i),"MCTruthRespFlavorVsGenJetPt",InputTag_list.at(i_i),f_i,InputTag_list.at(i_i));
//	test6a.plotAll();
//      }
//
//////    Extrapolation test1b("RelResponseVsRunNumber",samplelist.at(i));
//////    test1b.Plot();
//
//      //plot response of each flavor for different generators
//      for(size_t i_i=0;i_i<InputTag_list.size();i_i++){
//	OverlayPlainProfileAnd2DPlots test4a("|#eta|<1.3; " + samplelist.at(i)+"; "+InputTag_list.at(i_i),"Mean Response");
//	for(size_t f_i=0;f_i<flavor_list.size();f_i++){
//	  test4a.addPlots("FlavorResponse",samplelist.at(i),"MCTruthRespFlavorVsGenJetPt",InputTag_list.at(i_i),f_i,flavor_list.at(f_i));
//	}
//	test4a.plotAll();
//      }
//    
//      //plot response of each generator for different flavors (in barrel, endcap, HF) 
//    for(size_t f_i=0;f_i<flavor_list.size();f_i++){
//      OverlayPlainProfileAnd2DPlots test1a("|#eta|<1.3; " + samplelist.at(i)+"; "+flavor_list.at(f_i),"Mean Response");
//      for(size_t i_i=0;i_i<InputTag_list.size();i_i++)test1a.addPlots("FlavorResponse",samplelist.at(i),"MCTruthRespFlavorVsGenJetPt",InputTag_list.at(i_i),f_i,InputTag_list.at(i_i));
//      test1a.plotAll();
//
//      OverlayPlainProfileAnd2DPlots test1b("|#eta|<1.3 (G); " + samplelist.at(i)+"; "+flavor_list.at(f_i),"Gaussian Mean Response");
//      for(size_t i_i=0;i_i<InputTag_list.size();i_i++)test1b.addPlots("GaussFlavorResponse",samplelist.at(i),"MCTruthRespFlavorVsGenJetPt",InputTag_list.at(i_i),f_i,InputTag_list.at(i_i));
//      test1b.plotAll();
//
//      OverlayPlainProfileAnd2DPlots test2a("1.3<|#eta|<3; " + samplelist.at(i)+"; "+flavor_list.at(f_i),"Mean Response");
//      for(size_t i_i=0;i_i<InputTag_list.size();i_i++)test2a.addPlots("FlavorResponse",samplelist.at(i),"MCTruthEndCapRespFlavorVsGenJetPt",InputTag_list.at(i_i),f_i,InputTag_list.at(i_i));
//      test2a.plotAll();
//
//      OverlayPlainProfileAnd2DPlots test2b("3<|#eta|<5; " + samplelist.at(i)+"; "+flavor_list.at(f_i),"Mean Response");
//      for(size_t i_i=0;i_i<InputTag_list.size();i_i++)test2b.addPlots("FlavorResponse",samplelist.at(i),"MCTruthHFRespFlavorVsGenJetPt",InputTag_list.at(i_i),f_i,InputTag_list.at(i_i));
//      test2b.plotAll();
//
//
//    }


 
  }





//  OverlayPlainProfileAnd2DPlots test1b;
//  test1b.addPlots("OneBinRelResponseVsPt","44Z2");
//  test1b.addPlots("OneBinMPFResponseVsPt","44Z2");
//  test1b.addPlots("OneBinRelResponseVsPt","44Hpp");
//  test1b.addPlots("OneBinMPFResponseVsPt","44Hpp");
//  test1b.plotAll();


};


# ifndef __CINT__
int main()
{
  TestOverlayPlainProfileAnd2DPlots();

  return 0;
}
# endif
