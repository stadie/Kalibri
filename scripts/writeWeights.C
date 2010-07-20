// $Id: firstDataAnalysis.C,v 1.3 2010/01/18 15:55:04 mschrode Exp $

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
  writeWeightsToTree("/scratch/hh/current/cms/user/stadie/Spring10QCDDiJetPt0to15/ak5Calo*.root",0.1,2197029);
  writeWeightsToTree("/scratch/hh/current/cms/user/stadie/Spring10QCDDiJetPt1000to1400/ak5Calo*.root",0.1);
  writeWeightsToTree("/scratch/hh/current/cms/user/stadie/Spring10QCDDiJetPt120to170/ak5Calo*.root",0.1);
  writeWeightsToTree("/scratch/hh/current/cms/user/stadie/Spring10QCDDiJetPt1400to1800/ak5Calo*.root",0.1);
  writeWeightsToTree("/scratch/hh/current/cms/user/stadie/Spring10QCDDiJetPt15to20/ak5Calo*.root",0.1,2256430);
  writeWeightsToTree("/scratch/hh/current/cms/user/stadie/Spring10QCDDiJetPt170to230/ak5Calo*.root",0.1);
  writeWeightsToTree("/scratch/hh/current/cms/user/stadie/Spring10QCDDiJetPt1800to2200/ak5Calo*.root",0.1);
  writeWeightsToTree("/scratch/hh/current/cms/user/stadie/Spring10QCDDiJetPt20to30/ak5Calo*.root",0.1,1034680);
  writeWeightsToTree("/scratch/hh/current/cms/user/stadie/Spring10QCDDiJetPt2200to2600/ak5Calo*.root",0.1);
  writeWeightsToTree("/scratch/hh/current/cms/user/stadie/Spring10QCDDiJetPt230to300/ak5Calo*.root",0.1);
  writeWeightsToTree("/scratch/hh/current/cms/user/stadie/Spring10QCDDiJetPt2600to3000/ak5Calo*.root",0.1);
  writeWeightsToTree("/scratch/hh/current/cms/user/stadie/Spring10QCDDiJetPt3000to3500/ak5Calo*.root",0.1);
  writeWeightsToTree("/scratch/hh/current/cms/user/stadie/Spring10QCDDiJetPt300to380/ak5Calo*.root",0.1);
  writeWeightsToTree("/scratch/hh/current/cms/user/stadie/Spring10QCDDiJetPt30to50/ak5Calo*.root",0.1,1161768);
  writeWeightsToTree("/scratch/hh/current/cms/user/stadie/Spring10QCDDiJetPt380to470/ak5Calo*.root",0.1);
  writeWeightsToTree("/scratch/hh/current/cms/user/stadie/Spring10QCDDiJetPt470to600/ak5Calo*.root",0.1);
  writeWeightsToTree("/scratch/hh/current/cms/user/stadie/Spring10QCDDiJetPt50to80/ak5Calo*.root",0.1,111289);
  writeWeightsToTree("/scratch/hh/current/cms/user/stadie/Spring10QCDDiJetPt600to800/ak5Calo*.root",0.1);
  writeWeightsToTree("/scratch/hh/current/cms/user/stadie/Spring10QCDDiJetPt800to1000/ak5Calo*.root",0.1);
  writeWeightsToTree("/scratch/hh/current/cms/user/stadie/Spring10QCDDiJetPt80to120/ak5Calo*.root",0.1,606771);
}
