#include "TString.h"
#include "TFile.h"
#include "TChain.h"
#include "TH1D.h"
#include <iostream>

void getMCPUDists(){

  std::vector <TString> paths, shortNames;
  paths.push_back("/scratch/hh/current/cms/user/stadie/2011/v8/QCD_Pt-15to3000_TuneZ2_Flat_7TeV_pythia6_Fall11-PU_S6_START42_V14B-v1/merged/ak5FastPF_*.root");
  shortNames.push_back("42XFall11");

  paths.push_back("/scratch/hh/current/cms/user/kirschen/v6_QCD_Pt-15to3000_TuneZ2_Flat_7TeV_pythia6_Summer11-PU_S3_START42_V11-v2_smeared_Matthias_with_METcorr_nominal/ak5FastPF_*.root");
  shortNames.push_back("Summer11");

  paths.push_back("/scratch/hh/current/cms/user/kirschen/2011_Jets_v10/44X_Final/QCD_Pt-15to3000_TuneZ2_Flat_7TeV_pythia6_Fall11-PU_S6_START44_V9B-v1/merged/ak5FastPF_*.root");
  shortNames.push_back("44XFall11");

  paths.push_back("/scratch/hh/current/cms/user/kirschen/2011_Jets_v10/44X_Final/QCD_Pt-15to3000_Tune23_Flat_7TeV_herwigpp_Fall11-PU_S6_START44_V9B-v1/merged/ak5FastPF_*.root");
  shortNames.push_back("44XHerwigFall11");

  paths.push_back("/scratch/hh/current/cms/user/kirschen/2012_Jets_v2/QCD_Pt-15to3000_TuneZ2star_Flat_8TeV_pythia6_Summer12-PU_S7_START50_V15-v1/merged/ak5FastPF_*.root");
  shortNames.push_back("Summer12");

  paths.push_back("/scratch/hh/current/cms/user/kirschen/2012_Jets_v3/Merge_PUS6_PUS7_QCD_Pt-15to3000_TuneZ2_Flat_8TeV_pythia6_Summer12-_START52_V9-v1/ak5FastPF_*.root");
  shortNames.push_back("Summer12S6PlusS7");

  paths.push_back("/scratch/hh/current/cms/user/kirschen/2012_Jets_v3/Merge_PUS6Z2star1mio_PUS7Z210mio_QCD_Pt-15to3000_TuneZ2_Flat_8TeV_pythia6_Summer12-_START52_V9-v1/ak5FastPF_*.root");
  shortNames.push_back("Summer12S6Plus10MioS7");


  for(unsigned int i=0; i<paths.size();i++){
      TChain * chain = new TChain("DiJetTree","");
      chain->Add(paths.at(i));   
      //Observed
      TH1D* pileup = new TH1D("pileup","pileup",60,0,60);
      chain->Draw("PUMCNumVtx>>pileup");
      pileup->DrawNormalized();
      pileup->Scale(1/pileup->Integral());
      pileup->GetYaxis()->SetRangeUser(0,0.08);
      pileup->Draw();
      std::cout << "Double_t Observed"<<shortNames.at(i)<< "[60] = {\n";
      for(unsigned int j=1; j<60; j++){
	std::cout << pileup->GetBinContent(j) << ",";// <<std::endl;
      }
      std::cout << pileup->GetBinContent(60) << "}\n" << std::endl;
      TFile *outf = new TFile(shortNames.at(i)+"_ObservedMCPUDistributions.root","RECREATE");
      pileup->Write();
      outf->Close();

      //Truth
      chain->Draw("PUMCNumTruth>>pileup");
      pileup->DrawNormalized();
      pileup->Scale(1/pileup->Integral());
      pileup->GetYaxis()->SetRangeUser(0,0.08);
      pileup->Draw();
      std::cout << "Double_t True"<<shortNames.at(i)<< "[60] = {\n";
      for(unsigned int j=1; j<60; j++){
	std::cout << pileup->GetBinContent(j) << ",";// <<std::endl;
      }
      std::cout << pileup->GetBinContent(60) << "}\n" << std::endl;
      outf = new TFile(shortNames.at(i)+"_TrueMCPUDistributions.root","RECREATE");
      pileup->Write();
      outf->Close();

  }

}
