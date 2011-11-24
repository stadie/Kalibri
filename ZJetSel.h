//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Oct 20 13:34:03 2011 by ROOT version 5.30/01
// from TTree ZJetTree/
// found on file: /scratch/hh/current/cms/user/stadie/2011/v8/DoubleMu2011BPromptRecoV1_Cert_160404-178078/merged/ak5FastPF_1.root
//////////////////////////////////////////////////////////

#ifndef ZJetSel_h
#define ZJetSel_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>

class ZJetSel : public TSelector {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain

   // Declaration of leaf types
   UInt_t          RunNumber;
   UInt_t          LumiBlockNumber;
   UInt_t          EventNumber;
   Bool_t          HltPhysicsDelcared;
   Bool_t          HltL1Jet6U;
   Bool_t          HltDiJetAve15U;
   Bool_t          HltDiJetAve30U;
   Bool_t          HltDiJetAve50U;
   Bool_t          HltDiJetAve70U;
   Bool_t          HltDiJetAve100U;
   Bool_t          HltDiJetAve140U;
   Bool_t          HltDiJetAve180U;
   Bool_t          HltDiJetAve300U;
   Bool_t          HltDiJetAve30;
   Bool_t          HltDiJetAve60;
   Bool_t          HltDiJetAve80;
   Bool_t          HltDiJetAve110;
   Bool_t          HltDiJetAve150;
   Bool_t          HltDiJetAve190;
   Bool_t          HltDiJetAve240;
   Bool_t          HltDiJetAve300;
   Bool_t          HltDiJetAve370;
   Bool_t          HltJet30;
   Bool_t          HltJet60;
   Bool_t          HltJet80;
   Bool_t          HltJet110;
   Bool_t          HltJet150;
   Bool_t          HltJet190;
   Bool_t          HltJet240;
   Bool_t          HltJet300;
   Bool_t          HltJet370;
   Int_t           VtxN;
   Int_t           VtxNTracks;
   Float_t         VtxPosX;
   Float_t         VtxPosY;
   Float_t         VtxPosZ;
   Float_t         VtxNormalizedChi2;
   Float_t         VtxNDof;
   Bool_t          VtxIsFake;
   Int_t           PUMCNumVtx;
   Int_t           PUMCNumVtxOOT;
   Float_t         Rho;
   Int_t           NobjTow;
   Int_t           TowId[1000];   //[NobjTow]
   Int_t           TowId_phi[1000];   //[NobjTow]
   Int_t           TowId_eta[1000];   //[NobjTow]
   Float_t         TowEt[1000];   //[NobjTow]
   Float_t         TowEta[1000];   //[NobjTow]
   Float_t         TowPhi[1000];   //[NobjTow]
   Float_t         TowE[1000];   //[NobjTow]
   Float_t         TowEm[1000];   //[NobjTow]
   Float_t         TowHad[1000];   //[NobjTow]
   Float_t         TowOE[1000];   //[NobjTow]
   Int_t           Tow_jetidx[1000];   //[NobjTow]
   UInt_t          TowNumBadEcalCells[1000];   //[NobjTow]
   UInt_t          TowNumBadHcalCells[1000];   //[NobjTow]
   UInt_t          TowNumProblematicEcalCells[1000];   //[NobjTow]
   UInt_t          TowNumProblematicHcalCells[1000];   //[NobjTow]
   UInt_t          TowNumRecoveredEcalCells[1000];   //[NobjTow]
   UInt_t          TowNumRecoveredHcalCells[1000];   //[NobjTow]
   Int_t           NobjTrack;
   Int_t           TrackTowId[1000];   //[NobjTrack]
   Int_t           TrackTowIdPhi[1000];   //[NobjTrack]
   Int_t           TrackTowIdEta[1000];   //[NobjTrack]
   Int_t           TrackId[1000];   //[NobjTrack]
   Int_t           TrackNHits[1000];   //[NobjTrack]
   Bool_t          TrackQualityL[1000];   //[NobjTrack]
   Bool_t          TrackQualityT[1000];   //[NobjTrack]
   Bool_t          TrackQualityHP[1000];   //[NobjTrack]
   Float_t         TrackChi2[1000];   //[NobjTrack]
   Float_t         TrackPt[1000];   //[NobjTrack]
   Float_t         TrackEta[1000];   //[NobjTrack]
   Float_t         TrackPhi[1000];   //[NobjTrack]
   Float_t         TrackP[1000];   //[NobjTrack]
   Float_t         TrackDR[1000];   //[NobjTrack]
   Float_t         TrackPhiOut[1000];   //[NobjTrack]
   Float_t         TrackEtaOut[1000];   //[NobjTrack]
   Float_t         TrackDROut[1000];   //[NobjTrack]
   Float_t         TrackEMC1[1000];   //[NobjTrack]
   Float_t         TrackEMC3[1000];   //[NobjTrack]
   Float_t         TrackEMC5[1000];   //[NobjTrack]
   Float_t         TrackHAC1[1000];   //[NobjTrack]
   Float_t         TrackHAC3[1000];   //[NobjTrack]
   Float_t         TrackHAC5[1000];   //[NobjTrack]
   Int_t           Track_jetidx[1000];   //[NobjTrack]
   Float_t         MuDR[1000];   //[NobjTrack]
   Float_t         MuDE[1000];   //[NobjTrack]
   Float_t         TrackD0[1000];   //[NobjTrack]
   Float_t         TrackZ0[1000];   //[NobjTrack]
   Int_t           NobjJet;
   Float_t         JetPt[100];   //[NobjJet]
   Float_t         JetPhi[100];   //[NobjJet]
   Float_t         JetEta[100];   //[NobjJet]
   Float_t         JetEt[100];   //[NobjJet]
   Float_t         JetE[100];   //[NobjJet]
   Int_t           JetN90Hits[100];   //[NobjJet]
   Float_t         JetHad[100];   //[NobjJet]
   Float_t         JetEMF[100];   //[NobjJet]
   Float_t         JetFHPD[100];   //[NobjJet]
   Float_t         JetFRBX[100];   //[NobjJet]
   Float_t         JetFChargedHadron[100];   //[NobjJet]
   Float_t         JetNeutralHadrons[100];   //[NobjJet]
   Float_t         JetFPhotons[100];   //[NobjJet]
   Float_t         JetFElectrons[100];   //[NobjJet]
   Bool_t          JetIDLoose[100];   //[NobjJet]
   Bool_t          JetIDTight[100];   //[NobjJet]
   Float_t         JetEtWeightedSigmaPhi[100];   //[NobjJet]
   Float_t         JetEtWeightedSigmaEta[100];   //[NobjJet]
   Float_t         JetArea[100];   //[NobjJet]
   Float_t         JetCorrZSP[100];   //[NobjJet]
   Float_t         JetCorrL1[100];   //[NobjJet]
   Float_t         JetCorrL2[100];   //[NobjJet]
   Float_t         JetCorrL3[100];   //[NobjJet]
   Float_t         JetCorrJPT[100];   //[NobjJet]
   Float_t         JetCorrL2L3[100];   //[NobjJet]
   Float_t         JetCorrL2L3JPT[100];   //[NobjJet]
   Float_t         JetCorrL4JW[100];   //[NobjJet]
   Int_t           JetIEta[100];   //[NobjJet]
   Int_t           JetIPhi[100];   //[NobjJet]
   Int_t           JetNChargedHadrons[100];   //[NobjJet]
   Float_t         JetGenJetDeltaR[100];   //[NobjJet]
   Float_t         GenJetPt[100];   //[NobjJet]
   Float_t         GenJetPhi[100];   //[NobjJet]
   Float_t         GenJetEta[100];   //[NobjJet]
   Float_t         GenJetEt[100];   //[NobjJet]
   Float_t         GenJetE[100];   //[NobjJet]
   Int_t           NobjGenJet;
   Float_t         GenJetColPt[100];   //[NobjGenJet]
   Float_t         GenJetColPhi[100];   //[NobjGenJet]
   Float_t         GenJetColEta[100];   //[NobjGenJet]
   Float_t         GenJetColEt[100];   //[NobjGenJet]
   Float_t         GenJetColE[100];   //[NobjGenJet]
   Float_t         GenJetColEmE[100];   //[NobjGenJet]
   Float_t         GenJetColHadE[100];   //[NobjGenJet]
   Float_t         GenJetColInvE[100];   //[NobjGenJet]
   Float_t         GenJetColAuxE[100];   //[NobjGenJet]
   Int_t           GenJetColJetIdx[100];   //[NobjGenJet]
   Float_t         GenPartPt_algo[100];   //[NobjJet]
   Float_t         GenPartPhi_algo[100];   //[NobjJet]
   Float_t         GenPartEta_algo[100];   //[NobjJet]
   Float_t         GenPartEt_algo[100];   //[NobjJet]
   Float_t         GenPartE_algo[100];   //[NobjJet]
   Float_t         GenPartM_algo[100];   //[NobjJet]
   Int_t           GenPartId_algo[100];   //[NobjJet]
   Float_t         GenPartPt_phys[100];   //[NobjJet]
   Float_t         GenPartPhi_phys[100];   //[NobjJet]
   Float_t         GenPartEta_phys[100];   //[NobjJet]
   Float_t         GenPartEt_phys[100];   //[NobjJet]
   Float_t         GenPartE_phys[100];   //[NobjJet]
   Float_t         GenPartM_phys[100];   //[NobjJet]
   Int_t           GenPartId_phys[100];   //[NobjJet]
   Float_t         GenEvtScale;
   Float_t         Met;
   Float_t         MetPhi;
   Float_t         MetSum;
   Float_t         Weight;
   Float_t         CrossSection;
   Bool_t          HltDoubleMu3;
   Bool_t          HltDoubleMu5;
   Bool_t          HltDoubleMu6;
   Bool_t          HltDoubleMu7;
   Bool_t          HltDoubleMu45;
   Bool_t          HltMu5IsoMu5;
   Bool_t          HltMu13Mu8;
   Bool_t          HltMu17Mu8;
   Int_t           NobjMuon;
   Float_t         MuonPt[10];   //[NobjMuon]
   Float_t         MuonPhi[10];   //[NobjMuon]
   Float_t         MuonEta[10];   //[NobjMuon]
   Float_t         MuonE[10];   //[NobjMuon]
   Float_t         MuonTrackChi2[10];   //[NobjMuon]
   Float_t         MuonD0[10];   //[NobjMuon]
   Float_t         MuonIsoR03Trk[10];   //[NobjMuon]
   Float_t         MuonIsoR03Em[10];   //[NobjMuon]
   Float_t         MuonIsoR03Had[10];   //[NobjMuon]
   Float_t         MuonIsoR05Trk[10];   //[NobjMuon]
   Float_t         MuonIsoR05Em[10];   //[NobjMuon]
   Float_t         MuonIsoR05Had[10];   //[NobjMuon]
   Bool_t          MuonGlobal[10];   //[NobjMuon]
   Bool_t          MuonTracker[10];   //[NobjMuon]
   Bool_t          MuonStandalone[10];   //[NobjMuon]
   Bool_t          MuonCalo[10];   //[NobjMuon]
   Int_t           MuonNChambers[10];   //[NobjMuon]
   Int_t           MuonNTrackerHits[10];   //[NobjMuon]
   Int_t           MuonNPixelHits[10];   //[NobjMuon]
   Int_t           MuonNStations[10];   //[NobjMuon]
   Int_t           MuonTrackNDoF[10];   //[NobjMuon]
   Int_t           MuonCharge[10];   //[NobjMuon]
   Int_t           NobjGenMuon;
   Float_t         GenMuonPt[10];   //[NobjGenMuon]
   Float_t         GenMuonPhi[10];   //[NobjGenMuon]
   Float_t         GenMuonEta[10];   //[NobjGenMuon]
   Float_t         GenMuonE[10];   //[NobjGenMuon]

   // List of branches
   TBranch        *b_RunNumber;   //!
   TBranch        *b_LumiBlockNumber;   //!
   TBranch        *b_EventNumber;   //!
   TBranch        *b_HltPhysicsDelcared;   //!
   TBranch        *b_HltL1Jet6U;   //!
   TBranch        *b_HltDiJetAve15U;   //!
   TBranch        *b_HltDiJetAve30U;   //!
   TBranch        *b_HltDiJetAve50U;   //!
   TBranch        *b_HltDiJetAve70U;   //!
   TBranch        *b_HltDiJetAve100U;   //!
   TBranch        *b_HltDiJetAve140U;   //!
   TBranch        *b_HltDiJetAve180U;   //!
   TBranch        *b_HltDiJetAve300U;   //!
   TBranch        *b_HltDiJetAve30;   //!
   TBranch        *b_HltDiJetAve60;   //!
   TBranch        *b_HltDiJetAve80;   //!
   TBranch        *b_HltDiJetAve110;   //!
   TBranch        *b_HltDiJetAve150;   //!
   TBranch        *b_HltDiJetAve190;   //!
   TBranch        *b_HltDiJetAve240;   //!
   TBranch        *b_HltDiJetAve300;   //!
   TBranch        *b_HltDiJetAve370;   //!
   TBranch        *b_HltJet30;   //!
   TBranch        *b_HltJet60;   //!
   TBranch        *b_HltJet80;   //!
   TBranch        *b_HltJet110;   //!
   TBranch        *b_HltJet150;   //!
   TBranch        *b_HltJet190;   //!
   TBranch        *b_HltJet240;   //!
   TBranch        *b_HltJet300;   //!
   TBranch        *b_HltJet370;   //!
   TBranch        *b_VtxN;   //!
   TBranch        *b_VtxNTracks;   //!
   TBranch        *b_VtxPosX;   //!
   TBranch        *b_VtxPosY;   //!
   TBranch        *b_VtxPosZ;   //!
   TBranch        *b_VtxNormalizedChi2;   //!
   TBranch        *b_VtxNDof;   //!
   TBranch        *b_VtxIsFake;   //!
   TBranch        *b_PUMCNumVtx;   //!
   TBranch        *b_PUMCNumVtxOOT;   //!
   TBranch        *b_Rho;   //!
   TBranch        *b_NobjTow;   //!
   TBranch        *b_TowId;   //!
   TBranch        *b_TowId_phi;   //!
   TBranch        *b_TowId_eta;   //!
   TBranch        *b_TowEt;   //!
   TBranch        *b_TowEta;   //!
   TBranch        *b_TowPhi;   //!
   TBranch        *b_TowE;   //!
   TBranch        *b_TowEm;   //!
   TBranch        *b_TowHad;   //!
   TBranch        *b_TowOE;   //!
   TBranch        *b_Tow_jetidx;   //!
   TBranch        *b_TowNumBadEcalCells;   //!
   TBranch        *b_TowNumBadHcalCells;   //!
   TBranch        *b_TowNumProblematicEcalCells;   //!
   TBranch        *b_TowNumProblematicHcalCells;   //!
   TBranch        *b_TowNumRecoveredEcalCells;   //!
   TBranch        *b_TowNumRecoveredHcalCells;   //!
   TBranch        *b_NobjTrack;   //!
   TBranch        *b_TrackTowId;   //!
   TBranch        *b_TrackTowIdPhi;   //!
   TBranch        *b_TrackTowIdEta;   //!
   TBranch        *b_TrackId;   //!
   TBranch        *b_TrackNHits;   //!
   TBranch        *b_TrackQualityL;   //!
   TBranch        *b_TrackQualityT;   //!
   TBranch        *b_TrackQualityHP;   //!
   TBranch        *b_TrackChi2;   //!
   TBranch        *b_TrackPt;   //!
   TBranch        *b_TrackEta;   //!
   TBranch        *b_TrackPhi;   //!
   TBranch        *b_TrackP;   //!
   TBranch        *b_TrackDR;   //!
   TBranch        *b_TrackPhiOut;   //!
   TBranch        *b_TrackEtaOut;   //!
   TBranch        *b_TrackDROut;   //!
   TBranch        *b_TrackEMC1;   //!
   TBranch        *b_TrackEMC3;   //!
   TBranch        *b_TrackEMC5;   //!
   TBranch        *b_TrackHAC1;   //!
   TBranch        *b_TrackHAC3;   //!
   TBranch        *b_TrackHAC5;   //!
   TBranch        *b_Track_jetidx;   //!
   TBranch        *b_MuDR;   //!
   TBranch        *b_MuDE;   //!
   TBranch        *b_TrackD0;   //!
   TBranch        *b_TrackZ0;   //!
   TBranch        *b_NobjJet;   //!
   TBranch        *b_JetPt;   //!
   TBranch        *b_JetPhi;   //!
   TBranch        *b_JetEta;   //!
   TBranch        *b_JetEt;   //!
   TBranch        *b_JetE;   //!
   TBranch        *b_JetN90Hits;   //!
   TBranch        *b_JetHad;   //!
   TBranch        *b_JetEMF;   //!
   TBranch        *b_JetFHPD;   //!
   TBranch        *b_JetFRBX;   //!
   TBranch        *b_JetFChargedHadron;   //!
   TBranch        *b_JetNeutralHadrons;   //!
   TBranch        *b_JetFPhotons;   //!
   TBranch        *b_JetFElectrons;   //!
   TBranch        *b_JetIDLoose;   //!
   TBranch        *b_JetIDTight;   //!
   TBranch        *b_JetEtWeightedSigmaPhi;   //!
   TBranch        *b_JetEtWeightedSigmaEta;   //!
   TBranch        *b_JetArea;   //!
   TBranch        *b_JetCorrZSP;   //!
   TBranch        *b_JetCorrL1;   //!
   TBranch        *b_JetCorrL2;   //!
   TBranch        *b_JetCorrL3;   //!
   TBranch        *b_JetCorrJPT;   //!
   TBranch        *b_JetCorrL2L3;   //!
   TBranch        *b_JetCorrL2L3JPT;   //!
   TBranch        *b_JetCorrL4JW;   //!
   TBranch        *b_JetIEta;   //!
   TBranch        *b_JetIPhi;   //!
   TBranch        *b_JetNChargedHadrons;   //!
   TBranch        *b_JetGenJetDeltaR;   //!
   TBranch        *b_GenJetPt;   //!
   TBranch        *b_GenJetPhi;   //!
   TBranch        *b_GenJetEta;   //!
   TBranch        *b_GenJetEt;   //!
   TBranch        *b_GenJetE;   //!
   TBranch        *b_NobjGenJet;   //!
   TBranch        *b_GenJetColPt;   //!
   TBranch        *b_GenJetColPhi;   //!
   TBranch        *b_GenJetColEta;   //!
   TBranch        *b_GenJetColEt;   //!
   TBranch        *b_GenJetColE;   //!
   TBranch        *b_GenJetColEmE;   //!
   TBranch        *b_GenJetColHadE;   //!
   TBranch        *b_GenJetColInvE;   //!
   TBranch        *b_GenJetColAuxE;   //!
   TBranch        *b_GenJetColJetIdx;   //!
   TBranch        *b_GenPartPt_algo;   //!
   TBranch        *b_GenPartPhi_algo;   //!
   TBranch        *b_GenPartEta_algo;   //!
   TBranch        *b_GenPartEt_algo;   //!
   TBranch        *b_GenPartE_algo;   //!
   TBranch        *b_GenPartM_algo;   //!
   TBranch        *b_GenPartId_algo;   //!
   TBranch        *b_GenPartPt_phys;   //!
   TBranch        *b_GenPartPhi_phys;   //!
   TBranch        *b_GenPartEta_phys;   //!
   TBranch        *b_GenPartEt_phys;   //!
   TBranch        *b_GenPartE_phys;   //!
   TBranch        *b_GenPartM_phys;   //!
   TBranch        *b_GenPartId_phys;   //!
   TBranch        *b_GenEvtScale;   //!
   TBranch        *b_Met;   //!
   TBranch        *b_MetPhi;   //!
   TBranch        *b_MetSum;   //!
   TBranch        *b_Weight;   //!
   TBranch        *b_CrossSection;   //!
   TBranch        *b_HltDoubleMu3;   //!
   TBranch        *b_HltDoubleMu5;   //!
   TBranch        *b_HltDoubleMu6;   //!
   TBranch        *b_HltDoubleMu7;   //!
   TBranch        *b_HltDoubleMu45;   //!
   TBranch        *b_HltMu5IsoMu5;   //!
   TBranch        *b_HltMu13Mu8;   //!
   TBranch        *b_HltMu17Mu8;   //!
   TBranch        *b_NobjMuon;   //!
   TBranch        *b_MuonPt;   //!
   TBranch        *b_MuonPhi;   //!
   TBranch        *b_MuonEta;   //!
   TBranch        *b_MuonE;   //!
   TBranch        *b_MuonTrackChi2;   //!
   TBranch        *b_MuonD0;   //!
   TBranch        *b_MuonIsoR03Trk;   //!
   TBranch        *b_MuonIsoR03Em;   //!
   TBranch        *b_MuonIsoR03Had;   //!
   TBranch        *b_MuonIsoR05Trk;   //!
   TBranch        *b_MuonIsoR05Em;   //!
   TBranch        *b_MuonIsoR05Had;   //!
   TBranch        *b_MuonGlobal;   //!
   TBranch        *b_MuonTracker;   //!
   TBranch        *b_MuonStandalone;   //!
   TBranch        *b_MuonCalo;   //!
   TBranch        *b_MuonNChambers;   //!
   TBranch        *b_MuonNTrackerHits;   //!
   TBranch        *b_MuonNPixelHits;   //!
   TBranch        *b_MuonNStations;   //!
   TBranch        *b_MuonTrackNDoF;   //!
   TBranch        *b_MuonCharge;   //!
   TBranch        *b_NobjGenMuon;   //!
   TBranch        *b_GenMuonPt;   //!
   TBranch        *b_GenMuonPhi;   //!
   TBranch        *b_GenMuonEta;   //!
   TBranch        *b_GenMuonE;   //!

   ZJetSel(TTree * /*tree*/ =0) { }
   virtual ~ZJetSel() { }
   virtual Int_t   Version() const { return 2; }
   virtual void    Begin(TTree *tree);
   virtual void    SlaveBegin(TTree *tree);
   virtual void    Init(TTree *tree);
   virtual Bool_t  Notify();
   virtual Bool_t  Process(Long64_t entry);
   virtual Int_t   GetEntry(Long64_t entry, Int_t getall = 0) { return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0; }
   virtual void    SetOption(const char *option) { fOption = option; }
   virtual void    SetObject(TObject *obj) { fObject = obj; }
   virtual void    SetInputList(TList *input) { fInput = input; }
   virtual TList  *GetOutputList() const { return fOutput; }
   virtual void    SlaveTerminate();
   virtual void    Terminate();

  //ClassDef(ZJetSel,0);
};

#endif

#ifdef ZJetSel_cxx
void ZJetSel::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("RunNumber", &RunNumber, &b_RunNumber);
   fChain->SetBranchAddress("LumiBlockNumber", &LumiBlockNumber, &b_LumiBlockNumber);
   fChain->SetBranchAddress("EventNumber", &EventNumber, &b_EventNumber);
   fChain->SetBranchAddress("HltPhysicsDelcared", &HltPhysicsDelcared, &b_HltPhysicsDelcared);
   fChain->SetBranchAddress("HltL1Jet6U", &HltL1Jet6U, &b_HltL1Jet6U);
   fChain->SetBranchAddress("HltDiJetAve15U", &HltDiJetAve15U, &b_HltDiJetAve15U);
   fChain->SetBranchAddress("HltDiJetAve30U", &HltDiJetAve30U, &b_HltDiJetAve30U);
   fChain->SetBranchAddress("HltDiJetAve50U", &HltDiJetAve50U, &b_HltDiJetAve50U);
   fChain->SetBranchAddress("HltDiJetAve70U", &HltDiJetAve70U, &b_HltDiJetAve70U);
   fChain->SetBranchAddress("HltDiJetAve100U", &HltDiJetAve100U, &b_HltDiJetAve100U);
   fChain->SetBranchAddress("HltDiJetAve140U", &HltDiJetAve140U, &b_HltDiJetAve140U);
   fChain->SetBranchAddress("HltDiJetAve180U", &HltDiJetAve180U, &b_HltDiJetAve180U);
   fChain->SetBranchAddress("HltDiJetAve300U", &HltDiJetAve300U, &b_HltDiJetAve300U);
   fChain->SetBranchAddress("HltDiJetAve30", &HltDiJetAve30, &b_HltDiJetAve30);
   fChain->SetBranchAddress("HltDiJetAve60", &HltDiJetAve60, &b_HltDiJetAve60);
   fChain->SetBranchAddress("HltDiJetAve80", &HltDiJetAve80, &b_HltDiJetAve80);
   fChain->SetBranchAddress("HltDiJetAve110", &HltDiJetAve110, &b_HltDiJetAve110);
   fChain->SetBranchAddress("HltDiJetAve150", &HltDiJetAve150, &b_HltDiJetAve150);
   fChain->SetBranchAddress("HltDiJetAve190", &HltDiJetAve190, &b_HltDiJetAve190);
   fChain->SetBranchAddress("HltDiJetAve240", &HltDiJetAve240, &b_HltDiJetAve240);
   fChain->SetBranchAddress("HltDiJetAve300", &HltDiJetAve300, &b_HltDiJetAve300);
   fChain->SetBranchAddress("HltDiJetAve370", &HltDiJetAve370, &b_HltDiJetAve370);
   fChain->SetBranchAddress("HltJet30", &HltJet30, &b_HltJet30);
   fChain->SetBranchAddress("HltJet60", &HltJet60, &b_HltJet60);
   fChain->SetBranchAddress("HltJet80", &HltJet80, &b_HltJet80);
   fChain->SetBranchAddress("HltJet110", &HltJet110, &b_HltJet110);
   fChain->SetBranchAddress("HltJet150", &HltJet150, &b_HltJet150);
   fChain->SetBranchAddress("HltJet190", &HltJet190, &b_HltJet190);
   fChain->SetBranchAddress("HltJet240", &HltJet240, &b_HltJet240);
   fChain->SetBranchAddress("HltJet300", &HltJet300, &b_HltJet300);
   fChain->SetBranchAddress("HltJet370", &HltJet370, &b_HltJet370);
   fChain->SetBranchAddress("VtxN", &VtxN, &b_VtxN);
   fChain->SetBranchAddress("VtxNTracks", &VtxNTracks, &b_VtxNTracks);
   fChain->SetBranchAddress("VtxPosX", &VtxPosX, &b_VtxPosX);
   fChain->SetBranchAddress("VtxPosY", &VtxPosY, &b_VtxPosY);
   fChain->SetBranchAddress("VtxPosZ", &VtxPosZ, &b_VtxPosZ);
   fChain->SetBranchAddress("VtxNormalizedChi2", &VtxNormalizedChi2, &b_VtxNormalizedChi2);
   fChain->SetBranchAddress("VtxNDof", &VtxNDof, &b_VtxNDof);
   fChain->SetBranchAddress("VtxIsFake", &VtxIsFake, &b_VtxIsFake);
   fChain->SetBranchAddress("PUMCNumVtx", &PUMCNumVtx, &b_PUMCNumVtx);
   fChain->SetBranchAddress("PUMCNumVtxOOT", &PUMCNumVtxOOT, &b_PUMCNumVtxOOT);
   fChain->SetBranchAddress("Rho", &Rho, &b_Rho);
   fChain->SetBranchAddress("NobjTow", &NobjTow, &b_NobjTow);
   fChain->SetBranchAddress("TowId", &TowId, &b_TowId);
   fChain->SetBranchAddress("TowId_phi", &TowId_phi, &b_TowId_phi);
   fChain->SetBranchAddress("TowId_eta", &TowId_eta, &b_TowId_eta);
   fChain->SetBranchAddress("TowEt", &TowEt, &b_TowEt);
   fChain->SetBranchAddress("TowEta", &TowEta, &b_TowEta);
   fChain->SetBranchAddress("TowPhi", &TowPhi, &b_TowPhi);
   fChain->SetBranchAddress("TowE", &TowE, &b_TowE);
   fChain->SetBranchAddress("TowEm", &TowEm, &b_TowEm);
   fChain->SetBranchAddress("TowHad", &TowHad, &b_TowHad);
   fChain->SetBranchAddress("TowOE", &TowOE, &b_TowOE);
   fChain->SetBranchAddress("Tow_jetidx", &Tow_jetidx, &b_Tow_jetidx);
   fChain->SetBranchAddress("TowNumBadEcalCells", &TowNumBadEcalCells, &b_TowNumBadEcalCells);
   fChain->SetBranchAddress("TowNumBadHcalCells", &TowNumBadHcalCells, &b_TowNumBadHcalCells);
   fChain->SetBranchAddress("TowNumProblematicEcalCells", &TowNumProblematicEcalCells, &b_TowNumProblematicEcalCells);
   fChain->SetBranchAddress("TowNumProblematicHcalCells", &TowNumProblematicHcalCells, &b_TowNumProblematicHcalCells);
   fChain->SetBranchAddress("TowNumRecoveredEcalCells", &TowNumRecoveredEcalCells, &b_TowNumRecoveredEcalCells);
   fChain->SetBranchAddress("TowNumRecoveredHcalCells", &TowNumRecoveredHcalCells, &b_TowNumRecoveredHcalCells);
   fChain->SetBranchAddress("NobjTrack", &NobjTrack, &b_NobjTrack);
   fChain->SetBranchAddress("TrackTowId", &TrackTowId, &b_TrackTowId);
   fChain->SetBranchAddress("TrackTowIdPhi", &TrackTowIdPhi, &b_TrackTowIdPhi);
   fChain->SetBranchAddress("TrackTowIdEta", &TrackTowIdEta, &b_TrackTowIdEta);
   fChain->SetBranchAddress("TrackId", &TrackId, &b_TrackId);
   fChain->SetBranchAddress("TrackNHits", &TrackNHits, &b_TrackNHits);
   fChain->SetBranchAddress("TrackQualityL", &TrackQualityL, &b_TrackQualityL);
   fChain->SetBranchAddress("TrackQualityT", &TrackQualityT, &b_TrackQualityT);
   fChain->SetBranchAddress("TrackQualityHP", &TrackQualityHP, &b_TrackQualityHP);
   fChain->SetBranchAddress("TrackChi2", &TrackChi2, &b_TrackChi2);
   fChain->SetBranchAddress("TrackPt", &TrackPt, &b_TrackPt);
   fChain->SetBranchAddress("TrackEta", &TrackEta, &b_TrackEta);
   fChain->SetBranchAddress("TrackPhi", &TrackPhi, &b_TrackPhi);
   fChain->SetBranchAddress("TrackP", &TrackP, &b_TrackP);
   fChain->SetBranchAddress("TrackDR", &TrackDR, &b_TrackDR);
   fChain->SetBranchAddress("TrackPhiOut", &TrackPhiOut, &b_TrackPhiOut);
   fChain->SetBranchAddress("TrackEtaOut", &TrackEtaOut, &b_TrackEtaOut);
   fChain->SetBranchAddress("TrackDROut", &TrackDROut, &b_TrackDROut);
   fChain->SetBranchAddress("TrackEMC1", &TrackEMC1, &b_TrackEMC1);
   fChain->SetBranchAddress("TrackEMC3", &TrackEMC3, &b_TrackEMC3);
   fChain->SetBranchAddress("TrackEMC5", &TrackEMC5, &b_TrackEMC5);
   fChain->SetBranchAddress("TrackHAC1", &TrackHAC1, &b_TrackHAC1);
   fChain->SetBranchAddress("TrackHAC3", &TrackHAC3, &b_TrackHAC3);
   fChain->SetBranchAddress("TrackHAC5", &TrackHAC5, &b_TrackHAC5);
   fChain->SetBranchAddress("Track_jetidx", &Track_jetidx, &b_Track_jetidx);
   fChain->SetBranchAddress("MuDR", &MuDR, &b_MuDR);
   fChain->SetBranchAddress("MuDE", &MuDE, &b_MuDE);
   fChain->SetBranchAddress("TrackD0", &TrackD0, &b_TrackD0);
   fChain->SetBranchAddress("TrackZ0", &TrackZ0, &b_TrackZ0);
   fChain->SetBranchAddress("NobjJet", &NobjJet, &b_NobjJet);
   fChain->SetBranchAddress("JetPt", JetPt, &b_JetPt);
   fChain->SetBranchAddress("JetPhi", JetPhi, &b_JetPhi);
   fChain->SetBranchAddress("JetEta", JetEta, &b_JetEta);
   fChain->SetBranchAddress("JetEt", JetEt, &b_JetEt);
   fChain->SetBranchAddress("JetE", JetE, &b_JetE);
   fChain->SetBranchAddress("JetN90Hits", JetN90Hits, &b_JetN90Hits);
   fChain->SetBranchAddress("JetHad", JetHad, &b_JetHad);
   fChain->SetBranchAddress("JetEMF", JetEMF, &b_JetEMF);
   fChain->SetBranchAddress("JetFHPD", JetFHPD, &b_JetFHPD);
   fChain->SetBranchAddress("JetFRBX", JetFRBX, &b_JetFRBX);
   fChain->SetBranchAddress("JetFChargedHadron", JetFChargedHadron, &b_JetFChargedHadron);
   fChain->SetBranchAddress("JetNeutralHadrons", JetNeutralHadrons, &b_JetNeutralHadrons);
   fChain->SetBranchAddress("JetFPhotons", JetFPhotons, &b_JetFPhotons);
   fChain->SetBranchAddress("JetFElectrons", JetFElectrons, &b_JetFElectrons);
   fChain->SetBranchAddress("JetIDLoose", JetIDLoose, &b_JetIDLoose);
   fChain->SetBranchAddress("JetIDTight", JetIDTight, &b_JetIDTight);
   fChain->SetBranchAddress("JetEtWeightedSigmaPhi", JetEtWeightedSigmaPhi, &b_JetEtWeightedSigmaPhi);
   fChain->SetBranchAddress("JetEtWeightedSigmaEta", JetEtWeightedSigmaEta, &b_JetEtWeightedSigmaEta);
   fChain->SetBranchAddress("JetArea", JetArea, &b_JetArea);
   fChain->SetBranchAddress("JetCorrZSP", JetCorrZSP, &b_JetCorrZSP);
   fChain->SetBranchAddress("JetCorrL1", JetCorrL1, &b_JetCorrL1);
   fChain->SetBranchAddress("JetCorrL2", JetCorrL2, &b_JetCorrL2);
   fChain->SetBranchAddress("JetCorrL3", JetCorrL3, &b_JetCorrL3);
   fChain->SetBranchAddress("JetCorrJPT", JetCorrJPT, &b_JetCorrJPT);
   fChain->SetBranchAddress("JetCorrL2L3", JetCorrL2L3, &b_JetCorrL2L3);
   fChain->SetBranchAddress("JetCorrL2L3JPT", JetCorrL2L3JPT, &b_JetCorrL2L3JPT);
   fChain->SetBranchAddress("JetCorrL4JW", JetCorrL4JW, &b_JetCorrL4JW);
   fChain->SetBranchAddress("JetIEta", JetIEta, &b_JetIEta);
   fChain->SetBranchAddress("JetIPhi", JetIPhi, &b_JetIPhi);
   fChain->SetBranchAddress("JetNChargedHadrons", JetNChargedHadrons, &b_JetNChargedHadrons);
   fChain->SetBranchAddress("JetGenJetDeltaR", JetGenJetDeltaR, &b_JetGenJetDeltaR);
   fChain->SetBranchAddress("GenJetPt", GenJetPt, &b_GenJetPt);
   fChain->SetBranchAddress("GenJetPhi", GenJetPhi, &b_GenJetPhi);
   fChain->SetBranchAddress("GenJetEta", GenJetEta, &b_GenJetEta);
   fChain->SetBranchAddress("GenJetEt", GenJetEt, &b_GenJetEt);
   fChain->SetBranchAddress("GenJetE", GenJetE, &b_GenJetE);
   fChain->SetBranchAddress("NobjGenJet", &NobjGenJet, &b_NobjGenJet);
   fChain->SetBranchAddress("GenJetColPt", &GenJetColPt, &b_GenJetColPt);
   fChain->SetBranchAddress("GenJetColPhi", &GenJetColPhi, &b_GenJetColPhi);
   fChain->SetBranchAddress("GenJetColEta", &GenJetColEta, &b_GenJetColEta);
   fChain->SetBranchAddress("GenJetColEt", &GenJetColEt, &b_GenJetColEt);
   fChain->SetBranchAddress("GenJetColE", &GenJetColE, &b_GenJetColE);
   fChain->SetBranchAddress("GenJetColEmE", &GenJetColEmE, &b_GenJetColEmE);
   fChain->SetBranchAddress("GenJetColHadE", &GenJetColHadE, &b_GenJetColHadE);
   fChain->SetBranchAddress("GenJetColInvE", &GenJetColInvE, &b_GenJetColInvE);
   fChain->SetBranchAddress("GenJetColAuxE", &GenJetColAuxE, &b_GenJetColAuxE);
   fChain->SetBranchAddress("GenJetColJetIdx", &GenJetColJetIdx, &b_GenJetColJetIdx);
   fChain->SetBranchAddress("GenPartPt_algo", GenPartPt_algo, &b_GenPartPt_algo);
   fChain->SetBranchAddress("GenPartPhi_algo", GenPartPhi_algo, &b_GenPartPhi_algo);
   fChain->SetBranchAddress("GenPartEta_algo", GenPartEta_algo, &b_GenPartEta_algo);
   fChain->SetBranchAddress("GenPartEt_algo", GenPartEt_algo, &b_GenPartEt_algo);
   fChain->SetBranchAddress("GenPartE_algo", GenPartE_algo, &b_GenPartE_algo);
   fChain->SetBranchAddress("GenPartM_algo", GenPartM_algo, &b_GenPartM_algo);
   fChain->SetBranchAddress("GenPartId_algo", GenPartId_algo, &b_GenPartId_algo);
   fChain->SetBranchAddress("GenPartPt_phys", GenPartPt_phys, &b_GenPartPt_phys);
   fChain->SetBranchAddress("GenPartPhi_phys", GenPartPhi_phys, &b_GenPartPhi_phys);
   fChain->SetBranchAddress("GenPartEta_phys", GenPartEta_phys, &b_GenPartEta_phys);
   fChain->SetBranchAddress("GenPartEt_phys", GenPartEt_phys, &b_GenPartEt_phys);
   fChain->SetBranchAddress("GenPartE_phys", GenPartE_phys, &b_GenPartE_phys);
   fChain->SetBranchAddress("GenPartM_phys", GenPartM_phys, &b_GenPartM_phys);
   fChain->SetBranchAddress("GenPartId_phys", GenPartId_phys, &b_GenPartId_phys);
   fChain->SetBranchAddress("GenEvtScale", &GenEvtScale, &b_GenEvtScale);
   fChain->SetBranchAddress("Met", &Met, &b_Met);
   fChain->SetBranchAddress("MetPhi", &MetPhi, &b_MetPhi);
   fChain->SetBranchAddress("MetSum", &MetSum, &b_MetSum);
   fChain->SetBranchAddress("Weight", &Weight, &b_Weight);
   fChain->SetBranchAddress("CrossSection", &CrossSection, &b_CrossSection);
   fChain->SetBranchAddress("HltDoubleMu3", &HltDoubleMu3, &b_HltDoubleMu3);
   fChain->SetBranchAddress("HltDoubleMu5", &HltDoubleMu5, &b_HltDoubleMu5);
   fChain->SetBranchAddress("HltDoubleMu6", &HltDoubleMu6, &b_HltDoubleMu6);
   fChain->SetBranchAddress("HltDoubleMu7", &HltDoubleMu7, &b_HltDoubleMu7);
   fChain->SetBranchAddress("HltDoubleMu45", &HltDoubleMu45, &b_HltDoubleMu45);
   fChain->SetBranchAddress("HltMu5IsoMu5", &HltMu5IsoMu5, &b_HltMu5IsoMu5);
   fChain->SetBranchAddress("HltMu13Mu8", &HltMu13Mu8, &b_HltMu13Mu8);
   fChain->SetBranchAddress("HltMu17Mu8", &HltMu17Mu8, &b_HltMu17Mu8);
   fChain->SetBranchAddress("NobjMuon", &NobjMuon, &b_NobjMuon);
   fChain->SetBranchAddress("MuonPt", MuonPt, &b_MuonPt);
   fChain->SetBranchAddress("MuonPhi", MuonPhi, &b_MuonPhi);
   fChain->SetBranchAddress("MuonEta", MuonEta, &b_MuonEta);
   fChain->SetBranchAddress("MuonE", MuonE, &b_MuonE);
   fChain->SetBranchAddress("MuonTrackChi2", MuonTrackChi2, &b_MuonTrackChi2);
   fChain->SetBranchAddress("MuonD0", MuonD0, &b_MuonD0);
   fChain->SetBranchAddress("MuonIsoR03Trk", MuonIsoR03Trk, &b_MuonIsoR03Trk);
   fChain->SetBranchAddress("MuonIsoR03Em", MuonIsoR03Em, &b_MuonIsoR03Em);
   fChain->SetBranchAddress("MuonIsoR03Had", MuonIsoR03Had, &b_MuonIsoR03Had);
   fChain->SetBranchAddress("MuonIsoR05Trk", MuonIsoR05Trk, &b_MuonIsoR05Trk);
   fChain->SetBranchAddress("MuonIsoR05Em", MuonIsoR05Em, &b_MuonIsoR05Em);
   fChain->SetBranchAddress("MuonIsoR05Had", MuonIsoR05Had, &b_MuonIsoR05Had);
   fChain->SetBranchAddress("MuonGlobal", MuonGlobal, &b_MuonGlobal);
   fChain->SetBranchAddress("MuonTracker", MuonTracker, &b_MuonTracker);
   fChain->SetBranchAddress("MuonStandalone", MuonStandalone, &b_MuonStandalone);
   fChain->SetBranchAddress("MuonCalo", MuonCalo, &b_MuonCalo);
   fChain->SetBranchAddress("MuonNChambers", MuonNChambers, &b_MuonNChambers);
   fChain->SetBranchAddress("MuonNTrackerHits", MuonNTrackerHits, &b_MuonNTrackerHits);
   fChain->SetBranchAddress("MuonNPixelHits", MuonNPixelHits, &b_MuonNPixelHits);
   fChain->SetBranchAddress("MuonNStations", MuonNStations, &b_MuonNStations);
   fChain->SetBranchAddress("MuonTrackNDoF", MuonTrackNDoF, &b_MuonTrackNDoF);
   fChain->SetBranchAddress("MuonCharge", MuonCharge, &b_MuonCharge);
   fChain->SetBranchAddress("NobjGenMuon", &NobjGenMuon, &b_NobjGenMuon);
   fChain->SetBranchAddress("GenMuonPt", &GenMuonPt, &b_GenMuonPt);
   fChain->SetBranchAddress("GenMuonPhi", &GenMuonPhi, &b_GenMuonPhi);
   fChain->SetBranchAddress("GenMuonEta", &GenMuonEta, &b_GenMuonEta);
   fChain->SetBranchAddress("GenMuonE", &GenMuonE, &b_GenMuonE);
}

Bool_t ZJetSel::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

#endif // #ifdef ZJetSel_cxx
