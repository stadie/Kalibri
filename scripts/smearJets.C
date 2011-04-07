// $Id: smearJets.C,v 1.3 2011/02/04 10:35:44 mschrode Exp $

#include "TChain.h"
#include "TBranch.h"
#include "TFile.h"

#include <cmath>
#include <iostream>

#include "scripts/ptresolution.h"

void smearJetsInFile(const char* path,const char* outpath, int smeartype) {
  Int_t           NobjGenJet;
  Float_t         GenJetColPt[100];   //[NobjGenJet]
  Int_t           GenJetColJetIdx[100];   //[NobjGenJet]
  Int_t           NobjJet;
  Float_t         JetPtold[100];   //[NobjJet]
  Float_t         JetEtold[100];   //[NobjJet]
  Float_t         JetEold[100];   //[NobjJet]
  Float_t         JetPt[100];   //[NobjJet]
  Float_t         JetEt[100];   //[NobjJet]
  Float_t         JetEta[100];   //[NobjJet]
  Float_t         JetE[100];   //[NobjJet]
  Float_t         JetCorrL2L3[100];   //[NobjJet]
  //according to JME-10-14 table 1
  float smearFactors[4][4] = { 
    { 1.079, 1.054, 1.061, 1.174},//Calo
    { 1.059, 1.115, 0.993, 1.062},//JPT
    { 1.04,  1.02,  0.91,  1.04 },//PF
    { 1.07,  1.10,  1.07,  1.18 } //PF from Gamma
  };

  TFile *fin = TFile::Open(path);
  TTree *tin = (TTree*)fin->Get("DiJetTree");
  
  tin->SetBranchStatus("*",1);
  tin->SetBranchStatus("JetPt",0);
  tin->SetBranchStatus("JetEt",0);
  tin->SetBranchStatus("JetE",0);

  Long64_t nentries = tin->GetEntries();

  TFile *fout = TFile::Open(outpath,"UPDATE");
  TTree *tout =  tin->CloneTree(nentries,"fast");


  tin->SetBranchStatus("*",0); 
  tin->SetBranchStatus("NobjJet",1);
  tin->SetBranchStatus("JetPt",1);
  tin->SetBranchStatus("JetEt",1);
  tin->SetBranchStatus("JetEta",1);
  tin->SetBranchStatus("JetE",1);
  tin->SetBranchStatus("NobjGenJet",1);
  tin->SetBranchStatus("GenJetColPt",1);
  tin->SetBranchStatus("GenJetColJetIdx",1); 
  tin->SetBranchStatus("JetCorrL2L3",1); 

  tin->SetBranchAddress("NobjJet", &NobjJet);
  tin->SetBranchAddress("JetPt", JetPtold);
  tin->SetBranchAddress("JetEt", JetEtold);
  tin->SetBranchAddress("JetEta", JetEta);
  tin->SetBranchAddress("JetE", JetEold);
  tin->SetBranchAddress("NobjGenJet", &NobjGenJet);
  tin->SetBranchAddress("GenJetColPt", GenJetColPt);
  tin->SetBranchAddress("GenJetColJetIdx", GenJetColJetIdx);
  tin->SetBranchAddress("JetCorrL2L3", JetCorrL2L3);

  
  TBranch* b_JetPt = tout->Branch("JetPt",JetPt, "JetPt[NobjJet]/F");
  TBranch* b_JetEt = tout->Branch("JetEt",JetEt, "JetEt[NobjJet]/F");
  TBranch* b_JetE  = tout->Branch("JetE",JetE, "JetE[NobjJet]/F");
  //tout->SetBranchStatus("*",1);

  for (Long64_t i = 0; i < nentries; ++i) {
    tin->GetEntry(i);
    for(int j = 0 ; j <  NobjJet ; ++j) {
      JetPt[j] = JetPtold[j];
      JetEt[j] = JetPtold[j];
      JetE[j] = JetPtold[j];
    }
    for(int j = NobjGenJet -1  ; j >= 0 ; --j) {
      int id = GenJetColJetIdx[j];
      if((id >=  NobjJet)|| (id < 0)) continue;
      float diff = JetCorrL2L3[id] * JetPtold[id] - GenJetColPt[j];
      float aeta = std::abs(JetEta[id]);
      int ieta = 3;
      if (aeta < 1.1) ieta = 0;
      else if (aeta < 1.7) ieta = 1;
      else if (aeta < 2.3) ieta = 2;
      float smearfactor = 1.0;
      if(smeartype < 4) smearfactor = smearFactors[smeartype][ieta];
      else {
	if(smeartype == 4) {
	  _ismcjer = true;
	  float resmc = ptresolution_calo(GenJetColPt[j],JetEta[id]);
	  _ismcjer = false;
	  float resdata = ptresolution_calo(GenJetColPt[j],JetEta[id]);
	  smearfactor = resdata/resmc;	
	}
	if(smeartype == 5) {
	  _ismcjer = true;
	  float resmc = ptresolution_jpt(GenJetColPt[j],JetEta[id]);
	  _ismcjer = false;
	  float resdata = ptresolution_jpt(GenJetColPt[j],JetEta[id]);
	  smearfactor = resdata/resmc;
	}
	if(smeartype == 6) {
	  _ismcjer = true;
	  float resmc = ptresolution(GenJetColPt[j],JetEta[id]);
	  _ismcjer = false;
	  float resdata = ptresolution(GenJetColPt[j],JetEta[id]);
	  smearfactor = resdata/resmc;
	}
      }
      //std::cout << GenJetColPt[j] << ", " << JetEta[id] << ": " <<  smearfactor
      //		<< '\n';
      JetPt[id] =  JetPtold[id] + (smearfactor - 1.0) * diff /JetCorrL2L3[id];
      float k =  JetPt[id] /  JetPtold[id];
      //float k =  (1-smearfactor) * GenJetColPt[j]/JetPtold[id] + smearfactor;
      //JetPt[id] = k * JetPtold[id];
      JetEt[id] = k * JetEtold[id];
      JetE[id] = k * JetEold[id];
    }
    b_JetPt->Fill();
    b_JetEt->Fill();
    b_JetE->Fill();
  }
  fin->Close();
  tout->Write("", TObject::kOverwrite);
  fout->Close(); 
}

void smearJets() {
  TString sin("/scratch/hh/current/cms/user/stadie/QCD_Pt_15to3000_TuneZ2_Flat_7TeV_pythia6_Fall10-START38_V12-v1Dmerged/ak5Calo");
  TString son("/scratch/hh/current/cms/user/stadie/QCD_Pt_15to3000_TuneZ2_Flat_7TeV_pythia6_Fall10-START38_V12-v1Dsmeared3/ak5Calo");
  for(int i = 0 ; i < 10 ; ++i) {
    TString sif = sin + "_";
    sif += i;
    sif +=".root";
    std::cout << "smearing " << sif << '\n';
    TString sof = son + "_";
    sof += i;
    sof += ".root";
    smearJetsInFile(sif,sof,4);
  }
}
