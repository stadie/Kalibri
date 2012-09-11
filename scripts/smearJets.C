// $Id: smearJets.C,v 1.3 2012/09/11 12:57:39 kirschen Exp $

#include "TChain.h"
#include "TString.h"
#include "TBranch.h"
#include "TFile.h"
#include <TVector2.h>

#include <cmath>
#include <iostream>

#include "scripts/ptresolution.h"

void smearJetsInFile(const char* path,const char* outpath) {
  Int_t           NobjGenJet;
  Float_t         GenJetColPt[100];   //[NobjGenJet]
  Int_t           GenJetColJetIdx[100];   //[NobjGenJet]
  Int_t           NobjJet;
  Float_t         JetPtold[100];   //[NobjJet]
  Float_t         JetEtold[100];   //[NobjJet]
  Float_t         JetEold[100];   //[NobjJet]
  Float_t         JetPhiold[100];   //[NobjJet]
  Float_t         JetPt[100];   //[NobjJet]
  Float_t         JetEt[100];   //[NobjJet]
  Float_t         JetEta[100];   //[NobjJet]
  Float_t         JetE[100];   //[NobjJet] 
  Float_t         JetCorrL2L3[100];   //[NobjJet]
  Float_t         JetCorrL1[100];   //[NobjJet]
  Float_t         MetOld;
  Float_t         MetPhiOld;
  Float_t         Met;
  Float_t         MetPhi;

  Bool_t do_MET_smearing=true;
  Int_t count_nans_for_MET_smearing=1;

  TFile *fin = TFile::Open(path);
  TTree *tin = (TTree*)fin->Get("DiJetTree");
  
  tin->SetBranchStatus("*",1);
  tin->SetBranchStatus("JetPt",0);
  tin->SetBranchStatus("JetEt",0);
  tin->SetBranchStatus("JetE",0);
  tin->SetBranchStatus("Met",0); 
  tin->SetBranchStatus("MetPhi",0); 

  Long64_t nentries = tin->GetEntries();

  //  nentries=1000;

  TFile *fout = TFile::Open(outpath,"UPDATE");
  TTree *tout =  tin->CloneTree(nentries,"fast");


  tin->SetBranchStatus("*",0); 
  tin->SetBranchStatus("NobjJet",1);
  tin->SetBranchStatus("JetPt",1);
  tin->SetBranchStatus("JetEt",1);
  tin->SetBranchStatus("JetEta",1);
  tin->SetBranchStatus("JetE",1);
  tin->SetBranchStatus("JetPhi",1);
  tin->SetBranchStatus("NobjGenJet",1);
  tin->SetBranchStatus("GenJetColPt",1);
  tin->SetBranchStatus("GenJetColJetIdx",1); 
  tin->SetBranchStatus("JetCorrL2L3",1); 
  tin->SetBranchStatus("JetCorrL1",1); 
  tin->SetBranchStatus("Met",1); 
  tin->SetBranchStatus("MetPhi",1); 

  tin->SetBranchAddress("NobjJet", &NobjJet);
  tin->SetBranchAddress("JetPt", JetPtold);
  tin->SetBranchAddress("JetEt", JetEtold);
  tin->SetBranchAddress("JetEta", JetEta);
  tin->SetBranchAddress("JetE", JetEold);
  tin->SetBranchAddress("JetPhi", JetPhiold);
  tin->SetBranchAddress("NobjGenJet", &NobjGenJet);
  tin->SetBranchAddress("GenJetColPt", GenJetColPt);
  tin->SetBranchAddress("GenJetColJetIdx", GenJetColJetIdx);
  tin->SetBranchAddress("JetCorrL2L3", JetCorrL2L3);
  tin->SetBranchAddress("JetCorrL1", JetCorrL1);
  tin->SetBranchAddress("Met", &MetOld);
  tin->SetBranchAddress("MetPhi", &MetPhiOld);

  
  TBranch* b_JetPt  = tout->Branch("JetPt",JetPt, "JetPt[NobjJet]/F");
  TBranch* b_JetEt  = tout->Branch("JetEt",JetEt, "JetEt[NobjJet]/F");
  TBranch* b_JetE   = tout->Branch("JetE",JetE, "JetE[NobjJet]/F");
  TBranch* b_Met    = tout->Branch("Met",&Met, "Met/F");
  TBranch* b_MetPhi = tout->Branch("MetPhi",&MetPhi, "MetPhi/F");
  //tout->SetBranchStatus("*",1);

  for (Long64_t i = 0; i < nentries; ++i) {
    if(i % 100000 == 0)printf("event %d\n", i); //Alle tausend Events was sagen
    tin->GetEntry(i);
    for(int j = 0 ; j <  NobjJet ; ++j) {
      JetPt[j] = JetPtold[j];
      JetEt[j] = JetPtold[j];
      JetE[j] = JetPtold[j];
    }
    Met=MetOld;
    MetPhi=MetPhiOld;
    TVector2* MET = new TVector2(1,1);
    MET->SetMagPhi(Met,MetPhi);
    //    cout << "MAG: " << Met << " phi: " << MetPhi << endl;
    for(int j = NobjGenJet -1  ; j >= 0 ; --j) {
      int id = GenJetColJetIdx[j];
      if((id >=  NobjJet)|| (id < 0)) continue;
      //do not let JetCorrL1 be negative...
      float diff = max(JetCorrL1[id],(Float_t)0.0001) * JetCorrL2L3[id] * JetPtold[id] - GenJetColPt[j];
      float smearfactor = 1.0;
      smearfactor = getScaleFactor(JetEta[id]);
      if(TMath::IsNaN(smearfactor)){
	cout << "DEBUG... GenJetColPt[j]: " << GenJetColPt[j] << " JetEta[id]: " << JetEta[id] <<  endl;
      }
//      if(i % 100000 == 0){
//	printf("event %d\n", i); //Alle tausend Events was sagen
//	cout << smearfactor << endl;
//	cout << "DEBUG... GenJetColPt[j]: " << GenJetColPt[j] << " JetEta[id]: " << JetEta[id] <<  endl;
//      }

      // let L1JetCorr not be negative and make sure JetPt is above
      // some threshold
      if(GenJetColPt[j]>10.){
	JetPt[id] =  JetPtold[id] + (smearfactor - 1.0) * diff / (max(JetCorrL1[id],(Float_t)0.0001) * JetCorrL2L3[id]);
      }
      float k =  JetPt[id] /  JetPtold[id];
      //float k =  (1-smearfactor) * GenJetColPt[j]/JetPtold[id] + smearfactor;
      //JetPt[id] = k * JetPtold[id];
      JetEt[id] = k * JetEtold[id];
      JetE[id] = k * JetEold[id];
      if(do_MET_smearing){
	TVector2* MET_alteration = new TVector2(1,1);
	TVector2* pt_old = new TVector2(1,1);
	TVector2* pt_smeared = new TVector2(1,1);
	
	pt_old->SetMagPhi(JetPtold[id],JetPhiold[id]);
	pt_smeared->SetMagPhi(JetPt[id],JetPhiold[id]);
	*MET_alteration= *pt_old-*pt_smeared;
	//    MET_alteration->Print();
	//    MET->Print();
	if(TMath::IsNaN(MET_alteration->Mod())){
	  count_nans_for_MET_smearing++;
	  cout << "DEBUG... Ptold: " << JetPtold[id] << " JetPt: " << JetPt[id]<< " JetPhi: " << JetPhiold[id]  << " smearfactor: " << smearfactor << " JetCorrL2L3[id]: " << JetCorrL2L3[id] << " JetCorrL1[id]: " << JetCorrL1[id] << " MET-alteration failed: " << count_nans_for_MET_smearing << endl;
	}
	else{
	  *MET+=*MET_alteration;
	}
	//    MET->Print();
	MET_alteration->Delete();
	pt_old->Delete();
	pt_smeared->Delete();
      }
    }
    //        cout << "BEFORE... MAG: " << Met << " phi: " << MetPhi << endl;
    Met=MET->Mod();
    MetPhi=MET->Phi_mpi_pi(MET->Phi());
    //    cout << "AFTER.... MAG: " << Met << " phi: " << MetPhi << endl;

    b_JetPt->Fill();
    b_JetEt->Fill();
    b_JetE->Fill();
    b_Met->Fill();
    b_MetPhi->Fill();
  }
  fin->Close();
  tout->Write("", TObject::kOverwrite);
  fout->Close(); 
  cout << "MET-alterations failing due to NANs: " << count_nans_for_MET_smearing << endl;

}


void smearJets() {
  std::vector<TString> sInPath;
  std::vector<TString> sOutPath;
  std::vector<TString> scaleFactors;

  sInPath.push_back("/scratch/hh/current/cms/user/rathjd/Calibration/MCSummer12S10DX53/ak5FastPF");
  sOutPath.push_back("/scratch/hh/current/cms/user/kirschen/2012_Jets_v6/MCSummer12S10DX53_Smeared/ak5FastPF");
  scaleFactors.push_back("PF_Matthias");
  sInPath.push_back("/scratch/hh/current/cms/user/rathjd/Calibration/MCSummer12S10DX53/ak5PFCHS");
  sOutPath.push_back("/scratch/hh/current/cms/user/kirschen/2012_Jets_v6/MCSummer12S10DX53_Smeared/ak5PFCHS");
  scaleFactors.push_back("PF_Matthias");
  sInPath.push_back("/scratch/hh/current/cms/user/rathjd/Calibration/MCSummer12S10DX53/ak5FastPF");
  sOutPath.push_back("/scratch/hh/current/cms/user/kirschen/2012_Jets_v6/MCSummer12S10DX53_Smeared_u/ak5FastPF");
  scaleFactors.push_back("PF_Matthias_u");
  sInPath.push_back("/scratch/hh/current/cms/user/rathjd/Calibration/MCSummer12S10DX53/ak5PFCHS");
  sOutPath.push_back("/scratch/hh/current/cms/user/kirschen/2012_Jets_v6/MCSummer12S10DX53_Smeared_u/ak5PFCHS");
  scaleFactors.push_back("PF_Matthias_u");
  sInPath.push_back("/scratch/hh/current/cms/user/rathjd/Calibration/MCSummer12S10DX53/ak5FastPF");
  sOutPath.push_back("/scratch/hh/current/cms/user/kirschen/2012_Jets_v6/MCSummer12S10DX53_Smeared_d/ak5FastPF");
  scaleFactors.push_back("PF_Matthias_d");
  sInPath.push_back("/scratch/hh/current/cms/user/rathjd/Calibration/MCSummer12S10DX53/ak5PFCHS");
  sOutPath.push_back("/scratch/hh/current/cms/user/kirschen/2012_Jets_v6/MCSummer12S10DX53_Smeared_d/ak5PFCHS");
  scaleFactors.push_back("PF_Matthias_d");

  sInPath.push_back("/scratch/hh/current/cms/user/rathjd/Calibration/MCSummer12S10DX53/ak5Calo");
  sOutPath.push_back("/scratch/hh/current/cms/user/kirschen/2012_Jets_v6/MCSummer12S10DX53_Smeared/ak5Calo");
  scaleFactors.push_back("Calo_JME");




  assert(sInPath.size()==sOutPath.size());
  assert(sInPath.size()==scaleFactors.size());


  for(unsigned int sample_i = 0 ; sample_i < sOutPath.size() ; ++sample_i) {
    configureSmearfactor(scaleFactors.at(sample_i));
    for(int i = 0 ; i < 10 ; ++i) {
      TString sif = sInPath.at(sample_i)  + "_";
      sif += i;
      sif +="_sam0.root";
      std::cout << "smearing " << sif << '\n';
      TString sof = sOutPath.at(sample_i) + "_";
      sof += i;
      sof += "_sam0.root";
      smearJetsInFile(sif,sof);
    }
  }
}
