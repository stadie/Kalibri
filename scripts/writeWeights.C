// $Id: writeWeights.C,v 1.1 2010/07/20 17:13:38 stadie Exp $

#include "TChain.h"
#include "TBranch.h"
#include "TChainElement.h"
#include "TFile.h"

#include <iostream>

void writeWeightsToTree(const char* path, float lumi, int nev = -1) {
  TChain* c = new TChain("DiJetTree","DiJetTree");
  c->Add(path);
  
  Float_t weight;
  Float_t xsec;
  TBranch* bCrossSection;
  c->SetBranchAddress("CrossSection", &xsec,&bCrossSection);
  
  Long64_t nentries = c->GetEntries();
  if(nev < 0) nev = nentries;
  bCrossSection->GetEntry(0);
  Float_t newweight = (xsec * lumi) / nev;  
  std::cout << path << " and n events:" << nentries <<'\n';
  std::cout << "gen events:" << nev << " x-sec:" << xsec << " pb\n";
  std::cout << "new weight:" << newweight << '\n';
  TObjArray *fileElements=c->GetListOfFiles();
  TIter next(fileElements);
  TChainElement *chEl=0;
  while (( chEl=(TChainElement*)next() )) {
    //std::cout << "opening file " << chEl->GetTitle() << '\n';
    TFile f(chEl->GetTitle(),"update");
    TTree* djt = (TTree*)f.Get("DiJetTree");
    djt->SetBranchStatus("*",1);
    djt->SetBranchStatus("Weight",0);
    TTree *ndjt = djt->CloneTree(-1,"fast");
    delete djt;
    TBranch* bWeight = ndjt->Branch("Weight", &weight, "Weight/F");
    weight = newweight;
    //read the number of entries in the t3
    nentries = ndjt->GetEntries();
    for (Long64_t i = 0; i < nentries; i++){
      bWeight->Fill();
      //std::cout << bWeight->Fill() << '\n';
    }
    // save only the new version of the tree
    //f.Delete("DiJetTree;1");
    ndjt->Write("", TObject::kOverwrite);
    f.Close();
  }
  delete c;
}

void writeWeights() {
  writeWeightsToTree("/scratch/hh/current/cms/user/stadie/QCDDiJetSummer10-START36_V9_S09-v1A/Pt0to15*ak5Calo.root",0.1,2197029);
  writeWeightsToTree("/scratch/hh/current/cms/user/stadie/QCDDiJetSummer10-START36_V9_S09-v1A/Pt1000to1400*ak5Calo.root",0.1);
  writeWeightsToTree("/scratch/hh/current/cms/user/stadie/QCDDiJetSummer10-START36_V9_S09-v1A/Pt120to170*ak5Calo.root",0.1,58888);
  writeWeightsToTree("/scratch/hh/current/cms/user/stadie/QCDDiJetSummer10-START36_V9_S09-v1A/Pt1400to1800*ak5Calo.root",0.1);
  writeWeightsToTree("/scratch/hh/current/cms/user/stadie/QCDDiJetSummer10-START36_V9_S09-v1A/Pt15to20*ak5Calo.root",0.1,2256430);
  writeWeightsToTree("/scratch/hh/current/cms/user/stadie/QCDDiJetSummer10-START36_V9_S09-v1A/Pt170to230*ak5Calo.root",0.1,51680);
  writeWeightsToTree("/scratch/hh/current/cms/user/stadie/QCDDiJetSummer10-START36_V9_S09-v1A/Pt1800to2200*ak5Calo.root",0.1);
  writeWeightsToTree("/scratch/hh/current/cms/user/stadie/QCDDiJetSummer10-START36_V9_S09-v1A/Pt20to30*ak5Calo.root",0.1,1037110);
  writeWeightsToTree("/scratch/hh/current/cms/user/stadie/QCDDiJetSummer10-START36_V9_S09-v1A/Pt2200to2600*ak5Calo.root",0.1);
  writeWeightsToTree("/scratch/hh/current/cms/user/stadie/QCDDiJetSummer10-START36_V9_S09-v1A/Pt230to300*ak5Calo.root",0.1,52894);
  writeWeightsToTree("/scratch/hh/current/cms/user/stadie/QCDDiJetSummer10-START36_V9_S09-v1A/Pt2600to3000*ak5Calo.root",0.1);
  writeWeightsToTree("/scratch/hh/current/cms/user/stadie/QCDDiJetSummer10-START36_V9_S09-v1A/Pt3000to3500*ak5Calo.root",0.1);
  writeWeightsToTree("/scratch/hh/current/cms/user/stadie/QCDDiJetSummer10-START36_V9_S09-v1A/Pt300to380*ak5Calo.root",0.1,64265);
  writeWeightsToTree("/scratch/hh/current/cms/user/stadie/QCDDiJetSummer10-START36_V9_S09-v1A/Pt30to50*ak5Calo*.root",0.1,1161768);
  writeWeightsToTree("/scratch/hh/current/cms/user/stadie/QCDDiJetSummer10-START36_V9_S09-v1A/Pt380to470*ak5Calo.root",0.1);
  writeWeightsToTree("/scratch/hh/current/cms/user/stadie/QCDDiJetSummer10-START36_V9_S09-v1A/Pt470to600*ak5Calo.root",0.1);
  writeWeightsToTree("/scratch/hh/current/cms/user/stadie/QCDDiJetSummer10-START36_V9_S09-v1A/Pt50to80*ak5Calo.root",0.1,111289);
  writeWeightsToTree("/scratch/hh/current/cms/user/stadie/QCDDiJetSummer10-START36_V9_S09-v1A/Pt600to800*ak5Calo.root",0.1);
  writeWeightsToTree("/scratch/hh/current/cms/user/stadie/QCDDiJetSummer10-START36_V9_S09-v1A/Pt800to1000*ak5Calo.root",0.1);
  writeWeightsToTree("/scratch/hh/current/cms/user/stadie/QCDDiJetSummer10-START36_V9_S09-v1A/Pt80to120*ak5Calo.root",0.1,606771);
}

void testWeights() 
{
  TChain* c = new TChain("DiJetTree","DiJetTree");
  c->Add("/scratch/hh/current/cms/user/stadie/QCDDiJetSummer10-START36_V9_S09-v1A/Pt*ak5Calo.root");
  c->Draw("JetPt[0] >>hpt(700,0,3500)","Weight");
}
