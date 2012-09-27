#include "TH1.h"
#include "TH2.h"
#include "TMath.h"
#include "TF1.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TChain.h"
#include "cmath" 
#include "TPad.h"
#include "TStyle.h"
#include "TSpectrum.h"
#include "TString.h"
#include "TVector.h"
#include "TVector2.h"
#include "TLorentzVector.h"
#include "TFile.h"
#include <iostream>
#include <vector>
#include <list>

static inline double round(double d, int decPlaces) {
    d *= pow(10.,1.*decPlaces);
    d = std::floor(d+0.5);
    d /= pow(10.,1.*decPlaces);

    return d;
  }

void CaloPFMatching(int low, int high)
{
  //////////////////////////////////////////////////////////
//   This file has been automatically generated 
//     (Mon Sep 13 16:00:13 2010 by ROOT version5.26/00b)
//   from TTree DiJetTree/
//   found on file: /rdata/cms-data01/proof/calibration/JetMET_Run2010A-PromptReco-v4_DCSONLY_132440-143328_DiJetAve50U/Dijet-ak5Calo_101_1.root
//////////////////////////////////////////////////////////

//Reset ROOT and connect tree file
  gROOT->Reset();
  gROOT->SetStyle("Plain");
  bool verbose=false; //activation of debug mode
  if(verbose)std::cout<<"define eta binning"<<endl; 
  double etabins[8]={0.,0.783,1.305,1.93,2.5,2.964,3.2,5.191};
  if(verbose)std::cout<<"define histograms"<<endl;
  std::vector<TH2F*> Histograms;
  std::vector<TH2F*> HistsBalance;
  std::vector<TH2F*> HistsBalanceSame;
  std::vector<TH2F*> HistsBalanceRR;
  for(int number=0; number<7; number++)
    {
      Histograms.push_back(new TH2F(TString::Format("Eta%d",number),TString::Format("Matching %.3f<|#eta_{PF}|<%.3f",etabins[number],etabins[number+1]), 100, -100, 100, 100, 0, 1));
      Histograms[number]->Sumw2();
      Histograms[number]->GetXaxis()->SetTitle("p_{T}^{PF, jet1}-p_{T}^{Calo, jet1}");
      Histograms[number]->GetYaxis()->SetTitle("dR_{PF, Calo}^{jet1}");    
      Histograms[number]->SetDrawOption("colz");
      if(verbose)std::cout<<"hist "<<number<<" done"<<endl; 
      HistsBalance.push_back(new TH2F(TString::Format("BEta%d",number),TString::Format("Balance %.3f<|#eta_{PF}|<%.3f",etabins[number],etabins[number+1]), 300, 0, 3, 100, 0, 1));
      HistsBalance[number]->Sumw2();
      HistsBalance[number]->GetXaxis()->SetTitle("p_{T}^{Calo, jet1}/p_{T}^{PF, jet2}");
      HistsBalance[number]->GetYaxis()->SetTitle("dR_{PF, Calo}^{jet1}");    
      HistsBalance[number]->SetDrawOption("colz");
      HistsBalanceSame.push_back(new TH2F(TString::Format("BEta%dS",number),TString::Format("Balance %.3f<|#eta_{PF}|<%.3f",etabins[number],etabins[number+1]), 300, 0, 3, 100, 0, 1));
      HistsBalanceSame[number]->Sumw2();
      HistsBalanceSame[number]->GetXaxis()->SetTitle("p_{T}^{Calo, jet1}/p_{T}^{PF, jet1}");
      HistsBalanceSame[number]->GetYaxis()->SetTitle("dR_{PF, Calo}^{jet1}");    
      HistsBalanceSame[number]->SetDrawOption("colz");
      HistsBalanceRR.push_back(new TH2F(TString::Format("BEta%dRR",number),TString::Format("Balance %.3f<|#eta_{PF}|<%.3f",etabins[number],etabins[number+1]), 300, 0, 3, 100, 0, 1));
      HistsBalanceRR[number]->Sumw2();
      HistsBalanceRR[number]->GetXaxis()->SetTitle("2#frac{p_{T}^{Calo, jet1}-p_{T}^{PF, jet2}}{p_{T}^{Calo, jet1}+p_{T}^{PF, jet2}}");
      HistsBalanceRR[number]->GetYaxis()->SetTitle("dR_{PF, Calo}^{jet1}");    
      HistsBalanceRR[number]->SetDrawOption("colz");     
    }
  for(Int_t RunRange=low; RunRange<high; RunRange++) //min 0, max 13
    {
      //Vector-Definitionen
      if(verbose)std::cout<<"define and clean vectors"<<endl;
      std::vector<TLorentzVector*> PFJetVec1;
      int size=PFJetVec1.size();
      if(verbose)std::cout<<"size PFJV1: "<<size<<endl;
      if(size>0) for(int i=0; i<size;i++)
		   {
		     PFJetVec1[i]->Delete();
		   }		
      PFJetVec1.erase(PFJetVec1.begin(), PFJetVec1.end());
      std::vector<TLorentzVector*> PFJetVec2;
      size=PFJetVec2.size();
      if(verbose)std::cout<<"size PFJV2: "<<size<<endl;
      if(size>0) for(int i=0; i<size;i++)
		   {
		     PFJetVec2[i]->Delete();
		   }		
      PFJetVec2.erase(PFJetVec2.begin(), PFJetVec2.end());
      std::vector<UInt_t> RunNumberVec;
      size=RunNumberVec.size();
      if(verbose)std::cout<<"size RNV: "<<size<<endl;	
      RunNumberVec.erase(RunNumberVec.begin(),RunNumberVec.end());
      std::vector<UInt_t> EventNumberVec;
      size=EventNumberVec.size();
      if(verbose)std::cout<<"size ENV: "<<size<<endl;	
      EventNumberVec.erase(EventNumberVec.begin(), EventNumberVec.end());
      if(verbose)std::cout<<"start data readout"<<endl; 
      for(Int_t method=1; method<3;method++) //1 for PF, 2 for Calo
	{
	  if(verbose && method==1)std::cout<<"start PF readout"<<endl; 
	  if(verbose && method==2)std::cout<<"start Calo readout"<<endl;       
	  TChain* DiJetTree = new TChain("DiJetTree");
	  if(method==1) DiJetTree->Add("/scratch/hh/current/cms/user/rathjd/Calibration/2012_Jets_v10_noLowPu/combine_2012AJet_2012BJetMon-JetHT13JulyReReco_2012CJetMonJetHTv2_190456-202305/ak5FastPF_*_sam*.root");
	  if(method==2) DiJetTree->Add("/scratch/hh/current/cms/user/rathjd/Calibration/2012_Jets_v10_noLowPu/combine_2012AJet_2012BJetMon-JetHT13JulyReReco_2012CJetMonJetHTv2_190456-202305/ak5Calo_*_sam*.root");
	  if(verbose)std::cout<<"TChain configured"<<endl;
	  //Declaration of leaves types
	  UInt_t          RunNumber;
	  UInt_t          LumiBlockNumber;
	  UInt_t          EventNumber;
	  Bool_t          HltPhysicsDeclared;
	  Bool_t          HltDiPFJetAve40;
	  Bool_t          HltDiPFJetAve80;	  
	  Bool_t          HltDiPFJetAve140;
	  Bool_t          HltDiPFJetAve200;
	  Bool_t          HltDiPFJetAve260;
	  Bool_t          HltDiPFJetAve320;
	  Bool_t          HltDiPFJetAve400;
	  Int_t           VtxNTracks;
	  Float_t         VtxPosX;
	  Float_t         VtxPosY;
	  Float_t         VtxPosZ;
	  Float_t         VtxNormalizedChi2;
	  Float_t         VtxNDof;
	  Bool_t          VtxIsFake;
	  Int_t           NobjTow;
	  Int_t           TowId[1000];
	  Int_t           TowId_phi[1000];
	  Int_t           TowId_eta[1000];
	  Float_t         TowEt[1000];
	  Float_t         TowEta[1000];
	  Float_t         TowPhi[1000];
	  Float_t         TowE[1000];
	  Float_t         TowEm[1000];
	  Float_t         TowHad[1000];
	  Float_t         TowOE[1000];
	  Int_t           Tow_jetidx[1000];
	  UInt_t          TowNumBadEcalCells[1000];
	  UInt_t          TowNumBadHcalCells[1000];
	  UInt_t          TowNumProblematicEcalCells[1000];
	  UInt_t          TowNumProblematicHcalCells[1000];
	  UInt_t          TowNumRecoveredEcalCells[1000];
	  UInt_t          TowNumRecoveredHcalCells[1000];
	  Int_t           NobjTrack;
	  Int_t           TrackTowId[1000];
	  Int_t           TrackTowIdPhi[1000];
	  Int_t           TrackTowIdEta[1000];
	  Int_t           TrackId[1000];
	  Int_t           TrackNHits[1000];
	  Bool_t          TrackQualityL[1000];
	  Bool_t          TrackQualityT[1000];
	  Bool_t          TrackQualityHP[1000];
	  Float_t         TrackChi2[1000];
	  Float_t         TrackPt[1000];
	  Float_t         TrackEta[1000];
	  Float_t         TrackPhi[1000];
	  Float_t         TrackP[1000];
	  Float_t         TrackDR[1000];
	  Float_t         TrackPhiOut[1000];
	  Float_t         TrackEtaOut[1000];
	  Float_t         TrackDROut[1000];
	  Float_t         TrackEMC1[1000];
	  Float_t         TrackEMC3[1000];
	  Float_t         TrackEMC5[1000];
	  Float_t         TrackHAC1[1000];
	  Float_t         TrackHAC3[1000];
	  Float_t         TrackHAC5[1000];
	  Int_t           Track_jetidx[1000];
	  Float_t         MuDR[1000];
	  Float_t         MuDE[1000];
	  Float_t         TrackD0[1000];
	  Float_t         TrackZ0[1000];
	  Int_t           NobjJet;
	  Float_t         JetPt[200];
	  Float_t         JetPhi[200];
	  Float_t         JetEta[200];
	  Float_t         JetEt[200];
	  Float_t         JetE[200];
	  Int_t           JetN90Hits[200];
	  Float_t         JetHad[200];
	  Float_t         JetEMF[200];
	  Float_t         JetFHPD[200];
	  Float_t         JetFRBX[200];
	  Bool_t          JetIDLoose[200];
	  Bool_t          JetIDTight[200];
	  Float_t         JetEtWeightedSigmaPhi[200];
	  Float_t         JetEtWeightedSigmaEta[200];
	  Float_t         JetCorrZSP[200];
	  Float_t         JetCorrL1[200];
	  Float_t         JetCorrL2[200];
	  Float_t         JetCorrL3[200];
	  Float_t         JetCorrJPT[200];
	  Float_t         JetCorrL2L3[200];
	  Float_t         JetCorrL2L3JPT[200];
	  Float_t         JetCorrL4JW[200];
	  Int_t           JetIEta[200];
	  Int_t           JetIPhi[200];
	  Float_t         JetGenJetDeltaR[200];
	  Float_t         GenJetPt[200];
	  Float_t         GenJetPhi[200];
	  Float_t         GenJetEta[200];
	  Float_t         GenJetEt[200];
	  Float_t         GenJetE[200];
	  Int_t           NobjGenJet;
	  Float_t         GenJetColPt[0];
	  Float_t         GenJetColPhi[0];
	  Float_t         GenJetColEta[0];
	  Float_t         GenJetColEt[0];
	  Float_t         GenJetColE[0];
	  Float_t         GenJetColEmE[0];
	  Float_t         GenJetColHadE[0];
	  Float_t         GenJetColInvE[0];
	  Float_t         GenJetColAuxE[0];
	  Int_t           GenJetColJetIdx[0];
	  Float_t         GenPartPt_algo[200];
	  Float_t         GenPartPhi_algo[200];
	  Float_t         GenPartEta_algo[200];
	  Float_t         GenPartEt_algo[200];
	  Float_t         GenPartE_algo[200];
	  Float_t         GenPartM_algo[200];
	  Int_t           GenPartId_algo[200];
	  Float_t         GenPartPt_phys[200];
	  Float_t         GenPartPhi_phys[200];
	  Float_t         GenPartEta_phys[200];
	  Float_t         GenPartEt_phys[200];
	  Float_t         GenPartE_phys[200];
	  Float_t         GenPartM_phys[200];
	  Int_t           GenPartId_phys[200];
	  Float_t         GenEvtScale;
	  Float_t         Met;
	  Float_t         MetPhi;
	  Float_t         MetSum;
	  Float_t         Weight;
	  Float_t         CrossSection;
	  
	  // Set branch addresses.
	  DiJetTree->SetBranchAddress("RunNumber",&RunNumber);
	  DiJetTree->SetBranchAddress("LumiBlockNumber",&LumiBlockNumber);
	  DiJetTree->SetBranchAddress("EventNumber",&EventNumber);
	  DiJetTree->SetBranchAddress("HltPhysicsDelcared",&HltPhysicsDeclared);
	  DiJetTree->SetBranchAddress("HltDiPFJetAve40",&HltDiPFJetAve40);
	  DiJetTree->SetBranchAddress("HltDiPFJetAve80",&HltDiPFJetAve80);
	  DiJetTree->SetBranchAddress("HltDiPFJetAve140",&HltDiPFJetAve140);
	  DiJetTree->SetBranchAddress("HltDiPFJetAve200",&HltDiPFJetAve200);
	  DiJetTree->SetBranchAddress("HltDiPFJetAve260",&HltDiPFJetAve260);
	  DiJetTree->SetBranchAddress("HltDiPFJetAve320",&HltDiPFJetAve320);
	  DiJetTree->SetBranchAddress("HltDiPFJetAve400",&HltDiPFJetAve400);
	  DiJetTree->SetBranchAddress("VtxNTracks",&VtxNTracks);
	  DiJetTree->SetBranchAddress("VtxPosX",&VtxPosX);
	  DiJetTree->SetBranchAddress("VtxPosY",&VtxPosY);
	  DiJetTree->SetBranchAddress("VtxPosZ",&VtxPosZ);
	  DiJetTree->SetBranchAddress("VtxNormalizedChi2",&VtxNormalizedChi2);
	  DiJetTree->SetBranchAddress("VtxNDof",&VtxNDof);
	  DiJetTree->SetBranchAddress("VtxIsFake",&VtxIsFake);
	  DiJetTree->SetBranchAddress("NobjTow",&NobjTow);
	  DiJetTree->SetBranchAddress("TowId",TowId);
	  DiJetTree->SetBranchAddress("TowId_phi",TowId_phi);
	  DiJetTree->SetBranchAddress("TowId_eta",TowId_eta);
	  DiJetTree->SetBranchAddress("TowEt",TowEt);
	  DiJetTree->SetBranchAddress("TowEta",TowEta);
	  DiJetTree->SetBranchAddress("TowPhi",TowPhi);
	  DiJetTree->SetBranchAddress("TowE",TowE);
	  DiJetTree->SetBranchAddress("TowEm",TowEm);
	  DiJetTree->SetBranchAddress("TowHad",TowHad);
	  DiJetTree->SetBranchAddress("TowOE",TowOE);
	  DiJetTree->SetBranchAddress("Tow_jetidx",Tow_jetidx);
	  DiJetTree->SetBranchAddress("TowNumBadEcalCells",TowNumBadEcalCells);
	  DiJetTree->SetBranchAddress("TowNumBadHcalCells",TowNumBadHcalCells);
	  DiJetTree->SetBranchAddress("TowNumProblematicEcalCells",TowNumProblematicEcalCells);
	  DiJetTree->SetBranchAddress("TowNumProblematicHcalCells",TowNumProblematicHcalCells);
	  DiJetTree->SetBranchAddress("TowNumRecoveredEcalCells",TowNumRecoveredEcalCells);
	  DiJetTree->SetBranchAddress("TowNumRecoveredHcalCells",TowNumRecoveredHcalCells);
	  DiJetTree->SetBranchAddress("NobjTrack",&NobjTrack);
	  DiJetTree->SetBranchAddress("TrackTowId",TrackTowId);
	  DiJetTree->SetBranchAddress("TrackTowIdPhi",TrackTowIdPhi);
	  DiJetTree->SetBranchAddress("TrackTowIdEta",TrackTowIdEta);
	  DiJetTree->SetBranchAddress("TrackId",TrackId);
	  DiJetTree->SetBranchAddress("TrackNHits",TrackNHits);
	  DiJetTree->SetBranchAddress("TrackQualityL",TrackQualityL);
	  DiJetTree->SetBranchAddress("TrackQualityT",TrackQualityT);
	  DiJetTree->SetBranchAddress("TrackQualityHP",TrackQualityHP);
	  DiJetTree->SetBranchAddress("TrackChi2",TrackChi2);
	  DiJetTree->SetBranchAddress("TrackPt",TrackPt);
	  DiJetTree->SetBranchAddress("TrackEta",TrackEta);
	  DiJetTree->SetBranchAddress("TrackPhi",TrackPhi);
	  DiJetTree->SetBranchAddress("TrackP",TrackP);
	  DiJetTree->SetBranchAddress("TrackDR",TrackDR);
	  DiJetTree->SetBranchAddress("TrackPhiOut",TrackPhiOut);
	  DiJetTree->SetBranchAddress("TrackEtaOut",TrackEtaOut);
	  DiJetTree->SetBranchAddress("TrackDROut",TrackDROut);
	  DiJetTree->SetBranchAddress("TrackEMC1",TrackEMC1);
	  DiJetTree->SetBranchAddress("TrackEMC3",TrackEMC3);
	  DiJetTree->SetBranchAddress("TrackEMC5",TrackEMC5);
	  DiJetTree->SetBranchAddress("TrackHAC1",TrackHAC1);
	  DiJetTree->SetBranchAddress("TrackHAC3",TrackHAC3);
	  DiJetTree->SetBranchAddress("TrackHAC5",TrackHAC5);
	  DiJetTree->SetBranchAddress("Track_jetidx",Track_jetidx);
	  DiJetTree->SetBranchAddress("MuDR",MuDR);
	  DiJetTree->SetBranchAddress("MuDE",MuDE);
	  DiJetTree->SetBranchAddress("TrackD0",TrackD0);
	  DiJetTree->SetBranchAddress("TrackZ0",TrackZ0);
	  DiJetTree->SetBranchAddress("NobjJet",&NobjJet);
	  DiJetTree->SetBranchAddress("JetPt",JetPt);
	  DiJetTree->SetBranchAddress("JetPhi",JetPhi);
	  DiJetTree->SetBranchAddress("JetEta",JetEta);
	  DiJetTree->SetBranchAddress("JetEt",JetEt);
	  DiJetTree->SetBranchAddress("JetE",JetE);
	  DiJetTree->SetBranchAddress("JetN90Hits",JetN90Hits);
	  DiJetTree->SetBranchAddress("JetHad",JetHad);
	  DiJetTree->SetBranchAddress("JetEMF",JetEMF);
	  DiJetTree->SetBranchAddress("JetFHPD",JetFHPD);
	  DiJetTree->SetBranchAddress("JetFRBX",JetFRBX);
	  DiJetTree->SetBranchAddress("JetIDLoose",JetIDLoose);
	  DiJetTree->SetBranchAddress("JetIDTight",JetIDTight);
	  DiJetTree->SetBranchAddress("JetEtWeightedSigmaPhi",JetEtWeightedSigmaPhi);
	  DiJetTree->SetBranchAddress("JetEtWeightedSigmaEta",JetEtWeightedSigmaEta);
	  DiJetTree->SetBranchAddress("JetCorrZSP",JetCorrZSP);
	  DiJetTree->SetBranchAddress("JetCorrL1",JetCorrL1);
	  DiJetTree->SetBranchAddress("JetCorrL2",JetCorrL2);
	  DiJetTree->SetBranchAddress("JetCorrL3",JetCorrL3);
	  DiJetTree->SetBranchAddress("JetCorrJPT",JetCorrJPT);
	  DiJetTree->SetBranchAddress("JetCorrL2L3",JetCorrL2L3);
	  DiJetTree->SetBranchAddress("JetCorrL2L3JPT",JetCorrL2L3JPT);
	  DiJetTree->SetBranchAddress("JetCorrL4JW",JetCorrL4JW);
	  DiJetTree->SetBranchAddress("JetIEta",JetIEta);
	  DiJetTree->SetBranchAddress("JetIPhi",JetIPhi);
	  DiJetTree->SetBranchAddress("JetGenJetDeltaR",JetGenJetDeltaR);
	  DiJetTree->SetBranchAddress("GenJetPt",GenJetPt);
	  DiJetTree->SetBranchAddress("GenJetPhi",GenJetPhi);
	  DiJetTree->SetBranchAddress("GenJetEta",GenJetEta);
	  DiJetTree->SetBranchAddress("GenJetEt",GenJetEt);
	  DiJetTree->SetBranchAddress("GenJetE",GenJetE);
	  DiJetTree->SetBranchAddress("NobjGenJet",&NobjGenJet);
	  DiJetTree->SetBranchAddress("GenJetColPt",&GenJetColPt);
	  DiJetTree->SetBranchAddress("GenJetColPhi",&GenJetColPhi);
	  DiJetTree->SetBranchAddress("GenJetColEta",&GenJetColEta);
	  DiJetTree->SetBranchAddress("GenJetColEt",&GenJetColEt);
	  DiJetTree->SetBranchAddress("GenJetColE",&GenJetColE);
	  DiJetTree->SetBranchAddress("GenJetColEmE",&GenJetColEmE);
	  DiJetTree->SetBranchAddress("GenJetColHadE",&GenJetColHadE);
	  DiJetTree->SetBranchAddress("GenJetColInvE",&GenJetColInvE);
	  DiJetTree->SetBranchAddress("GenJetColAuxE",&GenJetColAuxE);
	  DiJetTree->SetBranchAddress("GenJetColJetIdx",&GenJetColJetIdx);
	  DiJetTree->SetBranchAddress("GenPartPt_algo",GenPartPt_algo);
	  DiJetTree->SetBranchAddress("GenPartPhi_algo",GenPartPhi_algo);
	  DiJetTree->SetBranchAddress("GenPartEta_algo",GenPartEta_algo);
	  DiJetTree->SetBranchAddress("GenPartEt_algo",GenPartEt_algo);
	  DiJetTree->SetBranchAddress("GenPartE_algo",GenPartE_algo);
	  DiJetTree->SetBranchAddress("GenPartM_algo",GenPartM_algo);
	  DiJetTree->SetBranchAddress("GenPartId_algo",GenPartId_algo);
	  DiJetTree->SetBranchAddress("GenPartPt_phys",GenPartPt_phys);
	  DiJetTree->SetBranchAddress("GenPartPhi_phys",GenPartPhi_phys);
	  DiJetTree->SetBranchAddress("GenPartEta_phys",GenPartEta_phys);
	  DiJetTree->SetBranchAddress("GenPartEt_phys",GenPartEt_phys);
	  DiJetTree->SetBranchAddress("GenPartE_phys",GenPartE_phys);
	  DiJetTree->SetBranchAddress("GenPartM_phys",GenPartM_phys);
	  DiJetTree->SetBranchAddress("GenPartId_phys",GenPartId_phys);
	  DiJetTree->SetBranchAddress("GenEvtScale",&GenEvtScale);
	  DiJetTree->SetBranchAddress("Met",&Met);
	  DiJetTree->SetBranchAddress("MetPhi",&MetPhi);
	  DiJetTree->SetBranchAddress("MetSum",&MetSum);
	  DiJetTree->SetBranchAddress("Weight",&Weight);
	  DiJetTree->SetBranchAddress("CrossSection",&CrossSection);
	  
	  //     This is the loop skeleton
	  //       To read only selected branches, Insert statements like:
	  DiJetTree->SetBranchStatus("*",0);  // disable all branches
	  DiJetTree->SetBranchStatus("Jet*",1);  // enable Jet branc
	  DiJetTree->SetBranchStatus("Event*",1);  // enable Event branc
	  DiJetTree->SetBranchStatus("Run*",1);  // enable Run branc
	  DiJetTree->SetBranchStatus("Hlt*",1);
	  if(verbose)std::cout<<"start event readout"<<endl;  
	  Long64_t nentries=DiJetTree->GetEntries();
	  //nentries = 1000; //switch on for testruns
	  Long64_t nbytes=0;  
	  std::cout<<"RunRange: "<<190000+RunRange*1000<<"-"<<(190000+(RunRange+1)*1000)<<endl;
	  for(Long64_t j=0; j<nentries; j++)
	    {
	      if(j %1000000==0)std::cout<<"event: "<<j<<endl;
	      nbytes += DiJetTree->GetEntry(j);
	      if(RunNumber<(190000+RunRange*1000) || RunNumber>=(190000+(RunRange+1)*1000)) continue;
	      int firstjet=-1;
	      int secondjet=-1;
	      bool dijets=false;
	      double ptave=0;
	      for(int i=0; i<NobjJet-1; i++) //find dijet events and save jet index
		{
		  if(!JetIDTight) continue;
		  if(JetPt[i]*JetCorrL2L3[i]*JetCorrL1[i]<20) continue;
		  if(firstjet < 0 && std::abs(JetEta[i]) > 2.5) continue;
		  else if(firstjet < 0) firstjet=i;
		  //if(verbose)std::cout<<"firstjet: "<<firstjet<<endl;
		  Double_t dphi = TVector2::Phi_mpi_pi(JetPhi[i]-JetPhi[firstjet]);
		  if(firstjet>=0 && std::abs(dphi)>2.7 && JetPt[i+1]*JetCorrL2L3[i+1]*JetCorrL1[i+1]/(JetPt[firstjet]*JetCorrL2L3[firstjet]*JetCorrL1[firstjet])<=0.2) secondjet=i;
		  //if(verbose)std::cout<<"secondjet: "<<secondjet<<endl;
		  if(firstjet>=0 && secondjet>=0)
		    {
		      dijets=true;
		      ptave=(JetPt[firstjet]*JetCorrL2L3[firstjet]*JetCorrL1[firstjet]+JetPt[secondjet]*JetCorrL2L3[secondjet]*JetCorrL1[secondjet])/2;
		      if(verbose)std::cout<<"found dijet event with ptave="<<ptave<<endl;
		      break;
		    }
		}
	      if(!dijets) continue;
	      //checking for trigger tresholds
	      bool triggered=false;
	      if(method==1)
		{
		  if(HltDiPFJetAve40 && ptave>=60) triggered=true;
		  if(HltDiPFJetAve80 && ptave>=105) triggered=true;
		  if(HltDiPFJetAve140 && ptave>=174) triggered=true;
		  if(HltDiPFJetAve200 && ptave>=242) triggered=true;
		  if(HltDiPFJetAve260 && ptave>=311) triggered=true;
		  if(HltDiPFJetAve320 && ptave>=380) triggered=true;
		  if(HltDiPFJetAve400 && ptave>=468) triggered=true;
		}
	      else if(method==2)
		{
		  if(HltDiPFJetAve40 && ptave>=65) triggered=true;
		  if(HltDiPFJetAve80 && ptave>=108) triggered=true;
		  if(HltDiPFJetAve140 && ptave>=183) triggered=true;
		  if(HltDiPFJetAve200 && ptave>=253) triggered=true;
		  if(HltDiPFJetAve260 && ptave>=324) triggered=true;
		  if(HltDiPFJetAve320 && ptave>=395) triggered=true;
		  if(HltDiPFJetAve400 && ptave>=482) triggered=true;
		}
	      if(verbose) std::cout<<HltDiPFJetAve40<<HltDiPFJetAve80<<HltDiPFJetAve140<<HltDiPFJetAve200<<HltDiPFJetAve200<<HltDiPFJetAve260<<HltDiPFJetAve320<<HltDiPFJetAve400<<endl;
	      if(verbose) std::cout<<"triggered: "<<triggered<<endl;
	      if(!triggered) continue;
	      if(method==1) //save dijet-events for further use
		{		  
		  if(verbose)std::cout<<"saving fourvectors for PF"<<endl;
		  TLorentzVector *One=new TLorentzVector(0.,0.,0.,0.);
		  One->SetPtEtaPhiE(JetPt[firstjet]*JetCorrL2L3[firstjet]*JetCorrL1[firstjet], JetEta[firstjet], JetPhi[firstjet], JetE[firstjet]);
		  TLorentzVector *Two=new TLorentzVector(0.,0.,0.,0.);
		  Two->SetPtEtaPhiE(JetPt[secondjet]*JetCorrL2L3[secondjet]*JetCorrL1[secondjet], JetEta[secondjet], JetPhi[secondjet], JetE[secondjet]);
		  bool found=false;
		  UInt_t redge=0;
		  bool first=false;
		  if(EventNumberVec.size()>0)redge=EventNumberVec.size()-1;
		  else first=true;
		  UInt_t ledge=0.;
		  UInt_t mid=0;
		  bool linear=false;
		  if(first)
		    {
		      if(verbose)std::cout<<"enter first vector elements"<<endl;
		      found=true;
		      PFJetVec1.push_back(One);
		      PFJetVec2.push_back(Two);
		      RunNumberVec.push_back(RunNumber);
		      EventNumberVec.push_back(EventNumber);
		    }
		  //binary ordering into the right positions by event and then runnumber
		  if(verbose)std::cout<<"enter sorted vector insertion"<<endl;
		  while(!found)
		    {
		      if(redge<ledge || ledge>EventNumberVec.size()-1 || redge>EventNumberVec.size()-1)
			{
			  found=true;
			  std::cout<<"no save position found"<<endl;
			  break;
			}
		      if(verbose)std::cout<<"l: "<<ledge<<", r: "<<redge<<", m: "<<mid<<", N-1:"<<EventNumberVec.size()-1<<endl;
		      if(verbose)std::cout<<round((redge+ledge)*0.5,0)<<endl;
		      if(!linear)mid=round((redge+ledge)*0.5,0);
		      if(mid>0 && mid<EventNumberVec.size()-1) //check for last bins
			{
			  if(EventNumberVec[mid]<EventNumber)
			    {
			      ledge=mid+1;//search right
			    }
			  else if(EventNumberVec[mid]==EventNumber)
			    {
			      if(RunNumberVec[mid]<RunNumber)
				{
				  if(EventNumberVec[mid+1]==EventNumber)//check event to the right
				    {
				      if(RunNumberVec[mid+1]<RunNumber)
					{
					  mid=mid+1;//search linear right
					  linear=true;
					}
				      else if(RunNumberVec[mid+1]==RunNumber)
					{
					  std::cout<<"doublecount error, skip event!"<<endl;
					  found=true;
					}
				      else
					{
					  PFJetVec1.insert(PFJetVec1.begin()+mid+1, One);
					  PFJetVec2.insert(PFJetVec2.begin()+mid+1, Two);
					  RunNumberVec.insert(RunNumberVec.begin()+mid+1, RunNumber);
					  EventNumberVec.insert(EventNumberVec.begin()+mid+1, EventNumber);
					  found=true;
					}
				    }
				  else
				    {
				      PFJetVec1.insert(PFJetVec1.begin()+mid+1, One);
				      PFJetVec2.insert(PFJetVec2.begin()+mid+1, Two);
				      RunNumberVec.insert(RunNumberVec.begin()+mid+1, RunNumber);
				      EventNumberVec.insert(EventNumberVec.begin()+mid+1, EventNumber);
				      found=true; 
				    }//end right event check
				}
			      else if(RunNumberVec[mid]==RunNumber)
				{
				  std::cout<<"doublecount error, skip event!"<<endl;
				  found=true;
				}
			      else
				{
				  if(EventNumberVec[mid-1]==EventNumber)//check event to the left
				    {
				      if(RunNumberVec[mid-1]<RunNumber)
					{
					  PFJetVec1.insert(PFJetVec1.begin()+mid, One);
					  PFJetVec2.insert(PFJetVec2.begin()+mid, Two);
					  RunNumberVec.insert(RunNumberVec.begin()+mid, RunNumber);
					  EventNumberVec.insert(EventNumberVec.begin()+mid, EventNumber);
					  found=true; 
					}
				      else if(RunNumberVec[mid-1]==RunNumber)
					{
					  std::cout<<"doublecount error, skip event!"<<endl;
					  found=true;
					}
				      else
					{
					  mid=mid-1;//search linear left
					  linear=true; 
					}
				    }
				  else
				    {
				      PFJetVec1.insert(PFJetVec1.begin()+mid, One);
				      PFJetVec2.insert(PFJetVec2.begin()+mid, Two);
				      RunNumberVec.insert(RunNumberVec.begin()+mid, RunNumber);
				      EventNumberVec.insert(EventNumberVec.begin()+mid, EventNumber);
				      found=true; 
				    }//end left event check
				}
			    }
			  else
			    {
			      if(EventNumberVec[mid-1]<EventNumber)
				{
				  PFJetVec1.insert(PFJetVec1.begin()+mid, One);
				  PFJetVec2.insert(PFJetVec2.begin()+mid, Two);
				  RunNumberVec.insert(RunNumberVec.begin()+mid, RunNumber);
				  EventNumberVec.insert(EventNumberVec.begin()+mid, EventNumber);
				  found=true;   
				}
			      else if(EventNumberVec[mid-1]==EventNumber)
				{
				  if(RunNumberVec[mid-1]<RunNumber)
				    {
				      PFJetVec1.insert(PFJetVec1.begin()+mid, One);
				      PFJetVec2.insert(PFJetVec2.begin()+mid, Two);
				      RunNumberVec.insert(RunNumberVec.begin()+mid, RunNumber);
				      EventNumberVec.insert(EventNumberVec.begin()+mid, EventNumber);
				      found=true;    
				    }
				  else if(RunNumberVec[mid-1]==RunNumber)
				    {
				      std::cout<<"doublecount error, skip event!"<<endl;
				      found=true; 
				    }
				  else
				    {
				      mid=mid-1;//search linear left
				      linear=true; 
				    }
				}
			      else
				{
				  redge=mid-1;//search right
				}
			    }
			}
		      else if(mid==0)
			{
			  if(EventNumberVec[mid]>EventNumber)
			    {
			      PFJetVec1.insert(PFJetVec1.begin()+mid, One);
			      PFJetVec2.insert(PFJetVec2.begin()+mid, Two);
			      RunNumberVec.insert(RunNumberVec.begin()+mid, RunNumber);
			      EventNumberVec.insert(EventNumberVec.begin()+mid, EventNumber);
			      found=true;   
			    }
			  else if(EventNumberVec[mid]==EventNumber)
			    {
			      if(RunNumberVec[mid]<RunNumber)
				{
				  PFJetVec1.insert(PFJetVec1.begin()+mid+1, One);
				  PFJetVec2.insert(PFJetVec2.begin()+mid+1, Two);
				  RunNumberVec.insert(RunNumberVec.begin()+mid+1, RunNumber);
				  EventNumberVec.insert(EventNumberVec.begin()+mid+1, EventNumber);
				  found=true;   
				}
			      else if(RunNumberVec[mid]==RunNumber)
				{
				  std::cout<<"doublecount error, skip event!"<<endl;
				  found=true; 
				}
			      else
				{
				  PFJetVec1.insert(PFJetVec1.begin()+mid, One);
				  PFJetVec2.insert(PFJetVec2.begin()+mid, Two);
				  RunNumberVec.insert(RunNumberVec.begin()+mid, RunNumber);
				  EventNumberVec.insert(EventNumberVec.begin()+mid, EventNumber);
				  found=true;   
				}
			    }
			  else
			    {
			      if(EventNumberVec[mid+1]<EventNumber)
				{
				  std::cout<<"missort error, skip event!"<<endl;
				  found=true; 
				}
			      if(EventNumberVec[mid+1]==EventNumber)
				{
				  if(RunNumberVec[mid+1]<RunNumber)
				    {
				      std::cout<<"missort error, skip event!"<<endl;
				      found=true; 
				    }
				  else if(RunNumberVec[mid+1]==RunNumber)
				    {
				      std::cout<<"doublecount error, skip event!"<<endl;
				      found=true; 
				    }
				  else
				    {
				      PFJetVec1.insert(PFJetVec1.begin()+mid+1, One);
				      PFJetVec2.insert(PFJetVec2.begin()+mid+1, Two);
				      RunNumberVec.insert(RunNumberVec.begin()+mid+1, RunNumber);
				      EventNumberVec.insert(EventNumberVec.begin()+mid+1, EventNumber);
				      found=true;   
				    }
				}
			      else
				{
				  PFJetVec1.insert(PFJetVec1.begin()+mid+1, One);
				  PFJetVec2.insert(PFJetVec2.begin()+mid+1, Two);
				  RunNumberVec.insert(RunNumberVec.begin()+mid+1, RunNumber);
				  EventNumberVec.insert(EventNumberVec.begin()+mid+1, EventNumber);
				  found=true;   
				}
			    }
			}
		      else if(mid==EventNumberVec.size()-1)
			{
			  if(EventNumberVec[mid]<EventNumber)
			    {
			      PFJetVec1.push_back(One);
			      PFJetVec2.push_back(Two);
			      RunNumberVec.push_back(RunNumber);
			      EventNumberVec.push_back(EventNumber);
			      found=true;   
			    }
			  else if(EventNumberVec[mid]==EventNumber)
			    {
			      if(RunNumberVec[mid]<RunNumber)
				{
				  PFJetVec1.push_back(One);
				  PFJetVec2.push_back(Two);
				  RunNumberVec.push_back(RunNumber);
				  EventNumberVec.push_back(EventNumber);
				  found=true;  
				}
			      else if(RunNumberVec[mid]==RunNumber)
				{
				  std::cout<<"doublecount error, skip event!"<<endl;
				  found=true; 
				}
			      else
				{
				  PFJetVec1.insert(PFJetVec1.begin()+mid, One);
				  PFJetVec2.insert(PFJetVec2.begin()+mid, Two);
				  RunNumberVec.insert(RunNumberVec.begin()+mid, RunNumber);
				  EventNumberVec.insert(EventNumberVec.begin()+mid, EventNumber);
				  found=true;    
				}
			    }
			  else
			    {
			      PFJetVec1.insert(PFJetVec1.begin()+mid, One);
			      PFJetVec2.insert(PFJetVec2.begin()+mid, Two);
			      RunNumberVec.insert(RunNumberVec.begin()+mid, RunNumber);
			      EventNumberVec.insert(EventNumberVec.begin()+mid, EventNumber);
			      found=true;   
			    }
			}
		    }//end while
		  if(verbose)std::cout<<"exit sorted vector filling"<<endl;
		}//end method 1
	      if(method==2 && EventNumberVec.size()>0) //check for matching events and match jets
		{
		  if(verbose)std::cout<<"matching PF to Calo"<<endl; 
		  bool found=false;
		  bool linear=false;
		  UInt_t redge=EventNumberVec.size()-1.;
		  UInt_t ledge=0.;
		  size_t browse;
		  UInt_t mid;
		  bool abort=false;
		  //binary search for the right position by event and then runnumber
		  if(verbose)std::cout<<"enter binary search for run and event match"<<endl;
		  while(!found)
		    {
		      if(verbose)std::cout<<"l: "<<ledge<<", r: "<<redge<<", m: "<<mid<<", N-1:"<<EventNumberVec.size()-1<<endl;
		      if(verbose)std::cout<<round((redge+ledge)*0.5,0)<<endl;
		      if(!linear)mid=round((redge+ledge)*0.5,0);
		      if(linear && EventNumberVec[mid]!=EventNumber)
			{
			  found=true;
			  abort=true;
			  std::cout<<"no matching element found in linear search<<"<<endl;
			}
		      if(redge<ledge || ledge>EventNumberVec.size()-1 || redge>EventNumberVec.size()-1)
			{
			  found=true;
			  abort=true;
			  std::cout<<"no matching element found"<<endl;
			  break;
			}
		      if(EventNumberVec[mid]<EventNumber)
			{
			  ledge=mid+1;//search right
			}
		      else if(EventNumberVec[mid]==EventNumber)
			{
			  if(RunNumberVec[mid]<RunNumber)
			    {
			      linear=true;
			      mid=mid+1;//linear search right
			    }
			  else if(RunNumberVec[mid]==RunNumber)
			    {
			      browse=mid;
			      found=true;
			    }
			  else
			    {
			      linear=true;
			      mid=mid-1;
			    }
			}
		      else
			{
			  redge=mid-1;//search left
			}
		    }
		  if(verbose)std::cout<<"exit binary search for run and event match"<<endl;
		  if(verbose && abort)std::cout<<"aborted due to no found element"<<endl;
		  if(abort) continue;
		  Double_t dphi11= TVector2::Phi_mpi_pi(PFJetVec1[browse]->Phi()-JetPhi[firstjet]); 
		  Double_t dR11=TMath::Sqrt(dphi11*dphi11+TMath::Power(PFJetVec1[browse]->Eta()-JetEta[firstjet],2));
		  Double_t dphi12= TVector2::Phi_mpi_pi(PFJetVec1[browse]->Phi()-JetPhi[secondjet]); 
		  Double_t dR12=TMath::Sqrt(dphi12*dphi12+TMath::Power(PFJetVec1[browse]->Eta()-JetEta[secondjet],2));
		  Double_t dphi21= TVector2::Phi_mpi_pi(PFJetVec2[browse]->Phi()-JetPhi[firstjet]); 
		  Double_t dR21=TMath::Sqrt(dphi21*dphi21+TMath::Power(PFJetVec2[browse]->Eta()-JetEta[firstjet],2));
		  Double_t dphi22= TVector2::Phi_mpi_pi(PFJetVec2[browse]->Phi()-JetPhi[secondjet]); 
		  Double_t dR22=TMath::Sqrt(dphi22*dphi22+TMath::Power(PFJetVec2[browse]->Eta()-JetEta[secondjet],2));
		  if(dR11<dR12 && dR11<dR21 && dR11<dR22)
		    {
		      if(verbose)std::cout<<"dR11: "<<dR11<<endl;
		      for(size_t number=0; number<7; number++)
			{
			  if(PFJetVec1[browse]->Eta()>=etabins[number] && PFJetVec1[browse]->Eta()<etabins[number+1])
			    {
			      Histograms[number]->Fill(PFJetVec1[browse]->Pt()-JetPt[firstjet]*JetCorrL2L3[firstjet]*JetCorrL1[firstjet], dR11);
			      HistsBalance[number]->Fill(JetPt[firstjet]*JetCorrL2L3[firstjet]*JetCorrL1[firstjet]/PFJetVec2[browse]->Pt(), dR11);
			      HistsBalanceSame[number]->Fill(JetPt[firstjet]*JetCorrL2L3[firstjet]*JetCorrL1[firstjet]/PFJetVec1[browse]->Pt(), dR11);
			      HistsBalanceRR[number]->Fill(2*(JetPt[firstjet]*JetCorrL2L3[firstjet]*JetCorrL1[firstjet]-PFJetVec2[browse]->Pt())/(PFJetVec1[browse]->Pt()+PFJetVec2[browse]->Pt()), dR11);
			    }
			  if(PFJetVec2[browse]->Eta()>=etabins[number] && PFJetVec2[browse]->Eta()<etabins[number+1])
			    {
			      Histograms[number]->Fill(PFJetVec2[browse]->Pt()-JetPt[secondjet]*JetCorrL2L3[secondjet]*JetCorrL1[secondjet], dR22);
			      HistsBalance[number]->Fill(JetPt[secondjet]*JetCorrL2L3[secondjet]*JetCorrL1[secondjet]/PFJetVec1[browse]->Pt(), dR22);
			      HistsBalanceSame[number]->Fill(JetPt[secondjet]*JetCorrL2L3[secondjet]*JetCorrL1[secondjet]/PFJetVec2[browse]->Pt(), dR22);
			      HistsBalanceRR[number]->Fill(2*(JetPt[secondjet]*JetCorrL2L3[secondjet]*JetCorrL1[secondjet]-PFJetVec1[browse]->Pt())/(PFJetVec1[browse]->Pt()+PFJetVec2[browse]->Pt()), dR22);			      
			    }
			}
		    }
		  else if(dR12<dR21 && dR12<dR22)
		    {
		      if(verbose)std::cout<<"dR12: "<<dR12<<endl;
		      for(size_t number=0; number<7; number++)
			{
			  if(PFJetVec1[browse]->Eta()>=etabins[number] && PFJetVec1[browse]->Eta()<etabins[number+1])
			    {
			      Histograms[number]->Fill(PFJetVec1[browse]->Pt()-JetPt[secondjet]*JetCorrL2L3[secondjet]*JetCorrL1[secondjet], dR12);
			      HistsBalance[number]->Fill(JetPt[secondjet]*JetCorrL2L3[secondjet]*JetCorrL1[secondjet]/PFJetVec2[browse]->Pt(), dR12);
			      HistsBalanceSame[number]->Fill(JetPt[secondjet]*JetCorrL2L3[secondjet]*JetCorrL1[secondjet]/PFJetVec1[browse]->Pt(), dR12);
			      HistsBalanceRR[number]->Fill(2*(JetPt[secondjet]*JetCorrL2L3[secondjet]*JetCorrL1[secondjet]-PFJetVec2[browse]->Pt())/(PFJetVec1[browse]->Pt()+PFJetVec2[browse]->Pt()), dR12);			      
			    }
			  if(PFJetVec2[browse]->Eta()>=etabins[number] && PFJetVec2[browse]->Eta()<etabins[number+1])
			    {
			      Histograms[number]->Fill(PFJetVec2[browse]->Pt()-JetPt[firstjet]*JetCorrL2L3[firstjet]*JetCorrL1[firstjet], dR21);
			      HistsBalance[number]->Fill(JetPt[firstjet]*JetCorrL2L3[firstjet]*JetCorrL1[firstjet]/PFJetVec1[browse]->Pt(), dR21);
			      HistsBalanceSame[number]->Fill(JetPt[firstjet]*JetCorrL2L3[firstjet]*JetCorrL1[firstjet]/PFJetVec2[browse]->Pt(), dR21);
			      HistsBalanceRR[number]->Fill(2*(JetPt[firstjet]*JetCorrL2L3[firstjet]*JetCorrL1[firstjet]-PFJetVec1[browse]->Pt())/(PFJetVec1[browse]->Pt()+PFJetVec2[browse]->Pt()), dR21);				      
			    }
			}
		    }
		  else if(dR21<dR22)
		    {
		      if(verbose)std::cout<<"dR21: "<<dR21<<endl;
		      for(size_t number=0; number<7; number++)
			{
			  if(PFJetVec1[browse]->Eta()>=etabins[number] && PFJetVec1[browse]->Eta()<etabins[number+1])
			    {
			      Histograms[number]->Fill(PFJetVec1[browse]->Pt()-JetPt[secondjet]*JetCorrL2L3[secondjet]*JetCorrL1[secondjet], dR12);
			      HistsBalance[number]->Fill(JetPt[secondjet]*JetCorrL2L3[secondjet]*JetCorrL1[secondjet]/PFJetVec2[browse]->Pt(), dR12);
			      HistsBalanceSame[number]->Fill(JetPt[secondjet]*JetCorrL2L3[secondjet]*JetCorrL1[secondjet]/PFJetVec1[browse]->Pt(), dR12);
			      HistsBalanceRR[number]->Fill(2*(JetPt[secondjet]*JetCorrL2L3[secondjet]*JetCorrL1[secondjet]-PFJetVec2[browse]->Pt())/(PFJetVec1[browse]->Pt()+PFJetVec2[browse]->Pt()), dR12);			      
			    }
			  if(PFJetVec2[browse]->Eta()>=etabins[number] && PFJetVec2[browse]->Eta()<etabins[number+1])
			    {
			      Histograms[number]->Fill(PFJetVec2[browse]->Pt()-JetPt[firstjet]*JetCorrL2L3[firstjet]*JetCorrL1[firstjet], dR21);
			      HistsBalance[number]->Fill(JetPt[firstjet]*JetCorrL2L3[firstjet]*JetCorrL1[firstjet]/PFJetVec1[browse]->Pt(), dR21);
			      HistsBalanceSame[number]->Fill(JetPt[firstjet]*JetCorrL2L3[firstjet]*JetCorrL1[firstjet]/PFJetVec2[browse]->Pt(), dR21);
			      HistsBalanceRR[number]->Fill(2*(JetPt[firstjet]*JetCorrL2L3[firstjet]*JetCorrL1[firstjet]-PFJetVec1[browse]->Pt())/(PFJetVec1[browse]->Pt()+PFJetVec2[browse]->Pt()), dR21);			      
			    }
			}
		    }
		  else
		    {  
		      if(verbose)std::cout<<"dR22: "<<dR22<<endl;
		      for(size_t number=0; number<7; number++)
			{
			  if(PFJetVec1[browse]->Eta()>=etabins[number] && PFJetVec1[browse]->Eta()<etabins[number+1])
			    {
			      Histograms[number]->Fill(PFJetVec1[browse]->Pt()-JetPt[firstjet]*JetCorrL2L3[firstjet]*JetCorrL1[firstjet], dR11);
			      HistsBalance[number]->Fill(JetPt[firstjet]*JetCorrL2L3[firstjet]*JetCorrL1[firstjet]/PFJetVec2[browse]->Pt(), dR11);
			      HistsBalanceSame[number]->Fill(JetPt[firstjet]*JetCorrL2L3[firstjet]*JetCorrL1[firstjet]/PFJetVec1[browse]->Pt(), dR11);
			      HistsBalanceRR[number]->Fill(2*(JetPt[firstjet]*JetCorrL2L3[firstjet]*JetCorrL1[firstjet]-PFJetVec2[browse]->Pt())/(PFJetVec1[browse]->Pt()+PFJetVec2[browse]->Pt()), dR11);			      
			    }
			  if(PFJetVec2[browse]->Eta()>=etabins[number] && PFJetVec2[browse]->Eta()<etabins[number+1])
			    {
			      Histograms[number]->Fill(PFJetVec2[browse]->Pt()-JetPt[secondjet]*JetCorrL2L3[secondjet]*JetCorrL1[secondjet], dR22);
			      HistsBalance[number]->Fill(JetPt[secondjet]*JetCorrL2L3[secondjet]*JetCorrL1[secondjet]/PFJetVec1[browse]->Pt(), dR22);
			      HistsBalanceSame[number]->Fill(JetPt[secondjet]*JetCorrL2L3[secondjet]*JetCorrL1[secondjet]/PFJetVec2[browse]->Pt(), dR22);
			      HistsBalanceRR[number]->Fill(2*(JetPt[secondjet]*JetCorrL2L3[secondjet]*JetCorrL1[secondjet]-PFJetVec1[browse]->Pt())/(PFJetVec1[browse]->Pt()+PFJetVec2[browse]->Pt()), dR22);				      
			    }
			}
		    }
		  /*if(verbose)std::cout<<"clean vectors"<<endl; 
		  //erase obsolete entries
		  PFJetVec1[browse]->Delete();
		  PFJetVec1.erase(PFJetVec1.begin()+browse);
		  PFJetVec2[browse]->Delete();		  
		  PFJetVec2.erase(PFJetVec2.begin()+browse);
		  RunNumberVec.erase(RunNumberVec.begin()+browse);
		  EventNumberVec.erase(EventNumberVec.begin()+browse);
		  if(verbose)std::cout<<"cleaned vectors"<<endl; */
		}//end method 2
	    }//end loop over tree
	}//end loop over samples
      //save data
      std::cout<<"PF1: "<<PFJetVec1.size()<<", PF2: "<<PFJetVec2.size()<<", Run: "<<RunNumberVec.size()<<", Evt: "<<EventNumberVec.size()<<endl;
    //}//end loop over runranges
  if(verbose)std::cout<<"save histograms"<<endl;
  TFile *f=new TFile(TString::Format("RunRange_%d-%d.root",190000+RunRange*1000,190000+(RunRange+1)*1000), "recreate", "MatchPFCalo", 1);
  for(size_t number=0; number<7; number++)
    {
      Histograms[number]->Write();
      HistsBalance[number]->Write();
      HistsBalanceSame[number]->Write();
      HistsBalanceRR[number]->Write();
    }
  f->Close(); 
  }//end loop over runranges
    for(size_t number=0; number<7; number++)
    {
      Histograms[number]->Delete();
      HistsBalance[number]->Delete();
      HistsBalanceSame[number]->Delete();
      HistsBalanceRR[number]->Delete();
    }
}

