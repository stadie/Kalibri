//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Apr 15 16:21:23 2010 by ROOT version 5.24/00
// from TTree DiJetTree/
// found on file: /rdata/cms-data01/proof/calibration/Summer09QCDFlat-MC_31X_V9_7TeV-v1_E/Summer09QCDFlat_Pt15to3000MC_31X_V9_7TeV-v1_92.root
//////////////////////////////////////////////////////////

#ifndef base_corr_h
#define base_corr_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH2.h>
#include <TH1.h>
#include <TProfile.h>
#include <TProfile.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TMath.h>
#include <iostream>
#include <iomanip>
#include <cmath>

class base_corr {

public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain


   //////////EIGENER TEIL

  std::stringstream testout;
  std::string s;

  
  const static Bool_t debug=0;
  const static Bool_t img_exp=1;
  const static Bool_t plots=1;
  std::stringstream img_extension;

  const static Int_t pt_bins_x = 100;
  const static Int_t s_phi_bins_x = 200;
  const static Int_t response_bins = 100;

  const static Double_t s_phi_xlow = -.5;
  const static Double_t s_phi_xhigh = 1.05;
  const static Double_t response_low = -0.5;
  const static Double_t response_high = 2.0;

  const static Double_t iso_max       = 0.5;

  Int_t no_pt_bins_;
  std::vector< std::pair <Double_t,Double_t> > pt_bins_;

  Int_t no_X_labels_;
  Int_t no_Corr_labels_;
  Int_t no_all_Corr_labels_;
  Int_t bin_choice;
  TString img_choice;
  Int_t X_choice;
  TString Corr_choice;

  std::vector <  TString > X_labels_;
  std::vector < std::pair < TString,TString > > Corr_labels_;
  std::vector <  TString > double_gauss_labels_;
  std::vector < Int_t > Corr_selected_labels_;


  TCanvas *test;

      Double_t _GenJetColE;	    
      Double_t _GenJetColEt;	    
      Double_t _GenJetColPt;	    
      Double_t _L2L3JetResponse;	    
      Double_t _L2L3JetE;	    
      Double_t _L2L3JetPt;	    
      Double_t _JetResponse;	    
      Double_t _JetE;	    
      Double_t _JetPt;	    
      Double_t _JetEt;	    
      Double_t _JetEtWeightedSigmaPhi;	    
      Double_t _JetEtWeightedSigmaEta;	    
      Double_t _JetEMF;         
      Double_t _JetEmE;
      Double_t _JetCorrEmE;
      Double_t _JetEMFCorr;




    /////base functions for variables, labels and histo-definitions 
      virtual void base_corr::Fill_obvious_vars(Int_t genjet_i, Int_t match);
      virtual void base_corr::Declare_Labels();
   TH1D* base_corr::define_pt_histo(TString Name, TString Title, std::vector< std::pair <Double_t,Double_t> > pt_bins_, Int_t pt_i,Int_t nbinsx, Double_t xlow, Double_t xup);
   TH2D* base_corr::define_pt_histo(TString Name, TString Title, std::vector< std::pair <Double_t,Double_t> > pt_bins_, Int_t pt_i,Int_t nbinsx, Double_t xlow, Double_t xup, Int_t nbinsy, Double_t ylow, Double_t yup);
   TProfile* base_corr::define_pt_histo(TString Name, TString Title, std::vector< std::pair <Double_t,Double_t> > pt_bins_, Int_t pt_i,Int_t nbinsx, Double_t xlow, Double_t xup, Double_t ylow, Double_t yup);
   TString base_corr::define_pt_histo_name(TString Name, std::vector< std::pair <Double_t,Double_t> > pt_bins_, Int_t pt_i);


   std::vector < Double_t > base_corr::get_EMF_vars_(Int_t genjet_i, Int_t match);
   std::vector < Double_t > base_corr::get_new_sigma_phi_vars_(Int_t genjet_i, Int_t match);
   std::vector < Double_t > base_corr::get_X_Vars_(Int_t genjet_i, Int_t match);
   std::vector < Double_t > base_corr::get_Correction_Vars_(Int_t genjet_i, Int_t match);
   std::vector<Double_t> base_corr::doublefit_gaus(TH1D *histo);

 

   ////////ENDE EIGENER TEIL





   // Declaration of leaf types
   UInt_t          RunNumber;
   UInt_t          LumiBlockNumber;
   UInt_t          EventNumber;
   Bool_t          HltPhysicsDelcared;
   Int_t           VtxNTracks;
   Float_t         VtxPosX;
   Float_t         VtxPosY;
   Float_t         VtxPosZ;
   Float_t         VtxNormalizedChi2;
   Float_t         VtxNDof;
   Bool_t          VtxIsFake;
   Int_t           NobjTow;
   Int_t           TowId[1500];   //[NobjTow]
   Int_t           TowId_phi[1500];   //[NobjTow]
   Int_t           TowId_eta[1500];   //[NobjTow]
   Float_t         TowEt[1500];   //[NobjTow]
   Float_t         TowEta[1500];   //[NobjTow]
   Float_t         TowPhi[1500];   //[NobjTow]
   Float_t         TowE[1500];   //[NobjTow]
   Float_t         TowEm[1500];   //[NobjTow]
   Float_t         TowHad[1500];   //[NobjTow]
   Float_t         TowOE[1500];   //[NobjTow]
   Int_t           Tow_jetidx[1500];   //[NobjTow]
   UInt_t          TowNumBadEcalCells[1500];   //[NobjTow]
   UInt_t          TowNumBadHcalCells[1500];   //[NobjTow]
   UInt_t          TowNumProblematicEcalCells[1500];   //[NobjTow]
   UInt_t          TowNumProblematicHcalCells[1500];   //[NobjTow]
   UInt_t          TowNumRecoveredEcalCells[1500];   //[NobjTow]
   UInt_t          TowNumRecoveredHcalCells[1500];   //[NobjTow]
   Int_t           NobjTrack;
   Int_t           TrackTowId[1500];   //[NobjTrack]
   Int_t           TrackTowIdPhi[1500];   //[NobjTrack]
   Int_t           TrackTowIdEta[1500];   //[NobjTrack]
   Int_t           TrackId[1500];   //[NobjTrack]
   Int_t           TrackNHits[1500];   //[NobjTrack]
   Bool_t          TrackQualityL[1500];   //[NobjTrack]
   Bool_t          TrackQualityT[1500];   //[NobjTrack]
   Bool_t          TrackQualityHP[1500];   //[NobjTrack]
   Float_t         TrackChi2[1500];   //[NobjTrack]
   Float_t         TrackPt[1500];   //[NobjTrack]
   Float_t         TrackEta[1500];   //[NobjTrack]
   Float_t         TrackPhi[1500];   //[NobjTrack]
   Float_t         TrackP[1500];   //[NobjTrack]
   Float_t         TrackDR[1500];   //[NobjTrack]
   Float_t         TrackPhiOut[1500];   //[NobjTrack]
   Float_t         TrackEtaOut[1500];   //[NobjTrack]
   Float_t         TrackDROut[1500];   //[NobjTrack]
   Float_t         TrackEMC1[1500];   //[NobjTrack]
   Float_t         TrackEMC3[1500];   //[NobjTrack]
   Float_t         TrackEMC5[1500];   //[NobjTrack]
   Float_t         TrackHAC1[1500];   //[NobjTrack]
   Float_t         TrackHAC3[1500];   //[NobjTrack]
   Float_t         TrackHAC5[1500];   //[NobjTrack]
   Int_t           Track_jetidx[1500];   //[NobjTrack]
   Float_t         MuDR[1500];   //[NobjTrack]
   Float_t         MuDE[1500];   //[NobjTrack]
   Float_t         TrackD0[1500];   //[NobjTrack]
   Float_t         TrackZ0[1500];   //[NobjTrack]
   Int_t           NobjJet;
   Float_t         JetPt[150];   //[NobjJet]
   Float_t         JetPhi[150];   //[NobjJet]
   Float_t         JetEta[150];   //[NobjJet]
   Float_t         JetEt[150];   //[NobjJet]
   Float_t         JetE[150];   //[NobjJet]
   Int_t           JetN90Hits[150];   //[NobjJet]
   Float_t         JetEMF[150];   //[NobjJet]
   Float_t         JetFHPD[150];   //[NobjJet]
   Float_t         JetFRBX[150];   //[NobjJet]
   Float_t         JetEtWeightedSigmaPhi[150];   //[NobjJet]
   Float_t         JetEtWeightedSigmaEta[150];   //[NobjJet]
   Float_t         JetCorrZSP[150];   //[NobjJet]
   Float_t         JetCorrL2[150];   //[NobjJet]
   Float_t         JetCorrL3[150];   //[NobjJet]
   Float_t         JetCorrJPT[150];   //[NobjJet]
   Float_t         JetCorrL2L3[150];   //[NobjJet]
   Float_t         JetCorrL2L3JPT[150];   //[NobjJet]
   Float_t         JetGenJetDeltaR[150];   //[NobjJet]
   Float_t         GenJetPt[150];   //[NobjJet]
   Float_t         GenJetPhi[150];   //[NobjJet]
   Float_t         GenJetEta[150];   //[NobjJet]
   Float_t         GenJetEt[150];   //[NobjJet]
   Float_t         GenJetE[150];   //[NobjJet]
   Int_t           NobjGenJet;
   Float_t         GenJetColPt[150];   //[NobjGenJet]
   Float_t         GenJetColPhi[150];   //[NobjGenJet]
   Float_t         GenJetColEta[150];   //[NobjGenJet]
   Float_t         GenJetColEt[150];   //[NobjGenJet]
   Float_t         GenJetColE[150];   //[NobjGenJet]
   Float_t         GenJetColEmE[150];   //[NobjGenJet]
   Float_t         GenJetColHadE[150];   //[NobjGenJet]
   Float_t         GenJetColInvE[150];   //[NobjGenJet]
   Float_t         GenJetColAuxE[150];   //[NobjGenJet]
   Int_t           GenJetColJetIdx[150];   //[NobjGenJet]
   Int_t           NobjGenJetPart;
   Float_t         GenJetPartE[1500];   //[NobjGenJetPart]
   Float_t         GenJetPartPt[1500];   //[NobjGenJetPart]
   Float_t         GenJetPartEta[1500];   //[NobjGenJetPart]
   Float_t         GenJetPartPhi[1500];   //[NobjGenJetPart]
   Int_t           GenJetPartPDG[1500];   //[NobjGenJetPart]
   Int_t           GenJetPartGenJetColIdx[1500];   //[NobjGenJetPart]
   Float_t         GenPartPt_algo[150];   //[NobjJet]
   Float_t         GenPartPhi_algo[150];   //[NobjJet]
   Float_t         GenPartEta_algo[150];   //[NobjJet]
   Float_t         GenPartEt_algo[150];   //[NobjJet]
   Float_t         GenPartE_algo[150];   //[NobjJet]
   Float_t         GenPartM_algo[150];   //[NobjJet]
   Int_t           GenPartId_algo[150];   //[NobjJet]
   Float_t         GenPartPt_phys[150];   //[NobjJet]
   Float_t         GenPartPhi_phys[150];   //[NobjJet]
   Float_t         GenPartEta_phys[150];   //[NobjJet]
   Float_t         GenPartEt_phys[150];   //[NobjJet]
   Float_t         GenPartE_phys[150];   //[NobjJet]
   Float_t         GenPartM_phys[150];   //[NobjJet]
   Int_t           GenPartId_phys[150];   //[NobjJet]
   Float_t         GenEvtScale;
   Float_t         Met;
   Float_t         MetPhi;
   Float_t         MetSum;
   Float_t         Weight;

   // List of branches
   TBranch        *b_RunNumber;   //!
   TBranch        *b_LumiBlockNumber;   //!
   TBranch        *b_EventNumber;   //!
   TBranch        *b_HltPhysicsDelcared;   //!
   TBranch        *b_VtxNTracks;   //!
   TBranch        *b_VtxPosX;   //!
   TBranch        *b_VtxPosY;   //!
   TBranch        *b_VtxPosZ;   //!
   TBranch        *b_VtxNormalizedChi2;   //!
   TBranch        *b_VtxNDof;   //!
   TBranch        *b_VtxIsFake;   //!
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
   TBranch        *b_JetEMF;   //!
   TBranch        *b_JetFHPD;   //!
   TBranch        *b_JetFRBX;   //!
   TBranch        *b_JetEtWeightedSigmaPhi;   //!
   TBranch        *b_JetEtWeightedSigmaEta;   //!
   TBranch        *b_JetCorrZSP;   //!
   TBranch        *b_JetCorrL2;   //!
   TBranch        *b_JetCorrL3;   //!
   TBranch        *b_JetCorrJPT;   //!
   TBranch        *b_JetCorrL2L3;   //!
   TBranch        *b_JetCorrL2L3JPT;   //!
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
   TBranch        *b_NobjGenJetPart;   //!
   TBranch        *b_GenJetPartE;   //!
   TBranch        *b_GenJetPartPt;   //!
   TBranch        *b_GenJetPartEta;   //!
   TBranch        *b_GenJetPartPhi;   //!
   TBranch        *b_GenJetPartPDG;   //!
   TBranch        *b_GenJetPartGenJetColIdx;   //!
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

   base_corr(TTree *tree=0);
   virtual ~base_corr();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);


};

#endif

#ifdef base_corr_cxx
base_corr::base_corr(TTree *tree)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.

   if (tree == 0) {
     /*
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/rdata/cms-data01/proof/calibration/MinBiasSpring10-START3X_V26A_356ReReco-v1/MinBiasSpring10-START3X_V26A_356ReReco-v1_99.root");
      if (!f) {
         f = new TFile("/rdata/cms-data01/proof/calibration/MinBiasSpring10-START3X_V26A_356ReReco-v1/MinBiasSpring10-START3X_V26A_356ReReco-v1_99.root");
      }
      tree = (TTree*)gDirectory->Get("DiJetTree");
*/

      TChain * chain = new TChain("DiJetTree","");
 

        chain->Add("/rdata/cms-data01/proof/calibration/Summer09QCDFlat-MC_31X_V9_7TeV-v1_E/Summer09QCDFlat_Pt15to3000MC_31X_V9_7TeV-v1_*.root");   
	//      chain->Add("/rdata/cms-data01/proof/calibration/MinBiasSpring10-START3X_V26A_356ReReco-v1/MinBiasSpring10-START3X_V26A_356ReReco-v1_*.root");

      tree = chain;



   }
   Init(tree);
test = new TCanvas("TEST", "TEST",341,360,700,530);
}

base_corr::~base_corr()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t base_corr::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t base_corr::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (!fChain->InheritsFrom(TChain::Class()))  return centry;
   TChain *chain = (TChain*)fChain;
   if (chain->GetTreeNumber() != fCurrent) {
      fCurrent = chain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void base_corr::Init(TTree *tree)
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
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("RunNumber", &RunNumber, &b_RunNumber);
   fChain->SetBranchAddress("LumiBlockNumber", &LumiBlockNumber, &b_LumiBlockNumber);
   fChain->SetBranchAddress("EventNumber", &EventNumber, &b_EventNumber);
   fChain->SetBranchAddress("HltPhysicsDelcared", &HltPhysicsDelcared, &b_HltPhysicsDelcared);
   fChain->SetBranchAddress("VtxNTracks", &VtxNTracks, &b_VtxNTracks);
   fChain->SetBranchAddress("VtxPosX", &VtxPosX, &b_VtxPosX);
   fChain->SetBranchAddress("VtxPosY", &VtxPosY, &b_VtxPosY);
   fChain->SetBranchAddress("VtxPosZ", &VtxPosZ, &b_VtxPosZ);
   fChain->SetBranchAddress("VtxNormalizedChi2", &VtxNormalizedChi2, &b_VtxNormalizedChi2);
   fChain->SetBranchAddress("VtxNDof", &VtxNDof, &b_VtxNDof);
   fChain->SetBranchAddress("VtxIsFake", &VtxIsFake, &b_VtxIsFake);
   fChain->SetBranchAddress("NobjTow", &NobjTow, &b_NobjTow);
   fChain->SetBranchAddress("TowId", TowId, &b_TowId);
   fChain->SetBranchAddress("TowId_phi", TowId_phi, &b_TowId_phi);
   fChain->SetBranchAddress("TowId_eta", TowId_eta, &b_TowId_eta);
   fChain->SetBranchAddress("TowEt", TowEt, &b_TowEt);
   fChain->SetBranchAddress("TowEta", TowEta, &b_TowEta);
   fChain->SetBranchAddress("TowPhi", TowPhi, &b_TowPhi);
   fChain->SetBranchAddress("TowE", TowE, &b_TowE);
   fChain->SetBranchAddress("TowEm", TowEm, &b_TowEm);
   fChain->SetBranchAddress("TowHad", TowHad, &b_TowHad);
   fChain->SetBranchAddress("TowOE", TowOE, &b_TowOE);
   fChain->SetBranchAddress("Tow_jetidx", Tow_jetidx, &b_Tow_jetidx);
   fChain->SetBranchAddress("TowNumBadEcalCells", TowNumBadEcalCells, &b_TowNumBadEcalCells);
   fChain->SetBranchAddress("TowNumBadHcalCells", TowNumBadHcalCells, &b_TowNumBadHcalCells);
   fChain->SetBranchAddress("TowNumProblematicEcalCells", TowNumProblematicEcalCells, &b_TowNumProblematicEcalCells);
   fChain->SetBranchAddress("TowNumProblematicHcalCells", TowNumProblematicHcalCells, &b_TowNumProblematicHcalCells);
   fChain->SetBranchAddress("TowNumRecoveredEcalCells", TowNumRecoveredEcalCells, &b_TowNumRecoveredEcalCells);
   fChain->SetBranchAddress("TowNumRecoveredHcalCells", TowNumRecoveredHcalCells, &b_TowNumRecoveredHcalCells);
   fChain->SetBranchAddress("NobjTrack", &NobjTrack, &b_NobjTrack);
   fChain->SetBranchAddress("TrackTowId", TrackTowId, &b_TrackTowId);
   fChain->SetBranchAddress("TrackTowIdPhi", TrackTowIdPhi, &b_TrackTowIdPhi);
   fChain->SetBranchAddress("TrackTowIdEta", TrackTowIdEta, &b_TrackTowIdEta);
   fChain->SetBranchAddress("TrackId", TrackId, &b_TrackId);
   fChain->SetBranchAddress("TrackNHits", TrackNHits, &b_TrackNHits);
   fChain->SetBranchAddress("TrackQualityL", TrackQualityL, &b_TrackQualityL);
   fChain->SetBranchAddress("TrackQualityT", TrackQualityT, &b_TrackQualityT);
   fChain->SetBranchAddress("TrackQualityHP", TrackQualityHP, &b_TrackQualityHP);
   fChain->SetBranchAddress("TrackChi2", TrackChi2, &b_TrackChi2);
   fChain->SetBranchAddress("TrackPt", TrackPt, &b_TrackPt);
   fChain->SetBranchAddress("TrackEta", TrackEta, &b_TrackEta);
   fChain->SetBranchAddress("TrackPhi", TrackPhi, &b_TrackPhi);
   fChain->SetBranchAddress("TrackP", TrackP, &b_TrackP);
   fChain->SetBranchAddress("TrackDR", TrackDR, &b_TrackDR);
   fChain->SetBranchAddress("TrackPhiOut", TrackPhiOut, &b_TrackPhiOut);
   fChain->SetBranchAddress("TrackEtaOut", TrackEtaOut, &b_TrackEtaOut);
   fChain->SetBranchAddress("TrackDROut", TrackDROut, &b_TrackDROut);
   fChain->SetBranchAddress("TrackEMC1", TrackEMC1, &b_TrackEMC1);
   fChain->SetBranchAddress("TrackEMC3", TrackEMC3, &b_TrackEMC3);
   fChain->SetBranchAddress("TrackEMC5", TrackEMC5, &b_TrackEMC5);
   fChain->SetBranchAddress("TrackHAC1", TrackHAC1, &b_TrackHAC1);
   fChain->SetBranchAddress("TrackHAC3", TrackHAC3, &b_TrackHAC3);
   fChain->SetBranchAddress("TrackHAC5", TrackHAC5, &b_TrackHAC5);
   fChain->SetBranchAddress("Track_jetidx", Track_jetidx, &b_Track_jetidx);
   fChain->SetBranchAddress("MuDR", MuDR, &b_MuDR);
   fChain->SetBranchAddress("MuDE", MuDE, &b_MuDE);
   fChain->SetBranchAddress("TrackD0", TrackD0, &b_TrackD0);
   fChain->SetBranchAddress("TrackZ0", TrackZ0, &b_TrackZ0);
   fChain->SetBranchAddress("NobjJet", &NobjJet, &b_NobjJet);
   fChain->SetBranchAddress("JetPt", JetPt, &b_JetPt);
   fChain->SetBranchAddress("JetPhi", JetPhi, &b_JetPhi);
   fChain->SetBranchAddress("JetEta", JetEta, &b_JetEta);
   fChain->SetBranchAddress("JetEt", JetEt, &b_JetEt);
   fChain->SetBranchAddress("JetE", JetE, &b_JetE);
   fChain->SetBranchAddress("JetN90Hits", JetN90Hits, &b_JetN90Hits);
   fChain->SetBranchAddress("JetEMF", JetEMF, &b_JetEMF);
   fChain->SetBranchAddress("JetFHPD", JetFHPD, &b_JetFHPD);
   fChain->SetBranchAddress("JetFRBX", JetFRBX, &b_JetFRBX);
   fChain->SetBranchAddress("JetEtWeightedSigmaPhi", JetEtWeightedSigmaPhi, &b_JetEtWeightedSigmaPhi);
   fChain->SetBranchAddress("JetEtWeightedSigmaEta", JetEtWeightedSigmaEta, &b_JetEtWeightedSigmaEta);
   fChain->SetBranchAddress("JetCorrZSP", JetCorrZSP, &b_JetCorrZSP);
   fChain->SetBranchAddress("JetCorrL2", JetCorrL2, &b_JetCorrL2);
   fChain->SetBranchAddress("JetCorrL3", JetCorrL3, &b_JetCorrL3);
   fChain->SetBranchAddress("JetCorrJPT", JetCorrJPT, &b_JetCorrJPT);
   fChain->SetBranchAddress("JetCorrL2L3", JetCorrL2L3, &b_JetCorrL2L3);
   fChain->SetBranchAddress("JetCorrL2L3JPT", JetCorrL2L3JPT, &b_JetCorrL2L3JPT);
   fChain->SetBranchAddress("JetGenJetDeltaR", JetGenJetDeltaR, &b_JetGenJetDeltaR);
   fChain->SetBranchAddress("GenJetPt", GenJetPt, &b_GenJetPt);
   fChain->SetBranchAddress("GenJetPhi", GenJetPhi, &b_GenJetPhi);
   fChain->SetBranchAddress("GenJetEta", GenJetEta, &b_GenJetEta);
   fChain->SetBranchAddress("GenJetEt", GenJetEt, &b_GenJetEt);
   fChain->SetBranchAddress("GenJetE", GenJetE, &b_GenJetE);
   fChain->SetBranchAddress("NobjGenJet", &NobjGenJet, &b_NobjGenJet);
   fChain->SetBranchAddress("GenJetColPt", GenJetColPt, &b_GenJetColPt);
   fChain->SetBranchAddress("GenJetColPhi", GenJetColPhi, &b_GenJetColPhi);
   fChain->SetBranchAddress("GenJetColEta", GenJetColEta, &b_GenJetColEta);
   fChain->SetBranchAddress("GenJetColEt", GenJetColEt, &b_GenJetColEt);
   fChain->SetBranchAddress("GenJetColE", GenJetColE, &b_GenJetColE);
   fChain->SetBranchAddress("GenJetColEmE", GenJetColEmE, &b_GenJetColEmE);
   fChain->SetBranchAddress("GenJetColHadE", GenJetColHadE, &b_GenJetColHadE);
   fChain->SetBranchAddress("GenJetColInvE", GenJetColInvE, &b_GenJetColInvE);
   fChain->SetBranchAddress("GenJetColAuxE", GenJetColAuxE, &b_GenJetColAuxE);
   fChain->SetBranchAddress("GenJetColJetIdx", GenJetColJetIdx, &b_GenJetColJetIdx);
   fChain->SetBranchAddress("NobjGenJetPart", &NobjGenJetPart, &b_NobjGenJetPart);
   fChain->SetBranchAddress("GenJetPartE", GenJetPartE, &b_GenJetPartE);
   fChain->SetBranchAddress("GenJetPartPt", GenJetPartPt, &b_GenJetPartPt);
   fChain->SetBranchAddress("GenJetPartEta", GenJetPartEta, &b_GenJetPartEta);
   fChain->SetBranchAddress("GenJetPartPhi", GenJetPartPhi, &b_GenJetPartPhi);
   fChain->SetBranchAddress("GenJetPartPDG", GenJetPartPDG, &b_GenJetPartPDG);
   fChain->SetBranchAddress("GenJetPartGenJetColIdx", GenJetPartGenJetColIdx, &b_GenJetPartGenJetColIdx);
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
   Notify();
}

Bool_t base_corr::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void base_corr::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t base_corr::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef base_corr_cxx
