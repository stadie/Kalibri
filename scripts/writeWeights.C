// $Id: writeWeights.C,v 1.2 2010/10/14 17:28:27 stadie Exp $

#include "TChain.h"
#include "TBranch.h"
#include "TChainElement.h"
#include "TFile.h"

#include <iostream>

void writeWeightsToTree(const char* path, float lumi, int nev = -1, double ptHatExpo = 0., double ptHatNorm = 1.) {
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
  std::cout << "new weight:" << newweight;
  if( ptHatExpo < 0. ) std::cout << " * ptHat^{" << ptHatExpo << "}";
  std::cout << '\n';
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
    Float_t ptHat = 0.;
    ndjt->SetBranchAddress("GenEvtScale",&ptHat);
    //read the number of entries in the t3
    nentries = ndjt->GetEntries();
    for (Long64_t i = 0; i < nentries; i++){
      if( ptHatExpo < 0. ) {
	ndjt->GetEntry(i);
	weight = newweight*pow(static_cast<double>(ptHat/ptHatNorm),ptHatExpo);
	//	if( i < 10 ) std::cout << "PtHat " << ptHat << " -- > " << weight << std::endl;
      }
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
  writeWeightsToTree("/scratch/hh/current/cms/user/mschrode/mc/tmp/QCD_Pt_15to3000_TuneZ2_Flat_7TeV_pythia6job_*_ak5JPT.root",1.,-1,-4.5,15.);
}

void testWeights() 
{
  TChain* c = new TChain("DiJetTree","DiJetTree");
  c->Add("/scratch/hh/current/cms/user/stadie/QCDDiJetSummer10-START36_V9_S09-v1A/Pt*ak5Calo.root");
  c->Draw("JetPt[0] >>hpt(700,0,3500)","Weight");
}
