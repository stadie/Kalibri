//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu May  8 15:44:54 2008 by ROOT version 5.18/00
// from TTree DiJetTree/
// found on file: NJet_test.root
//////////////////////////////////////////////////////////

#ifndef NJetSel_h
#define NJetSel_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>

class NJetSel : public TSelector {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain

   // Declaration of leaf types
   Int_t           NobjTow;
   Int_t           NobjTrack;
   Int_t           TowId[10000];   //[NobjTow]
   Int_t           TowId_phi[10000];   //[NobjTow]
   Int_t           TowId_eta[10000];   //[NobjTow]
   Float_t         TowEt[10000];   //[NobjTow]
   Float_t         TowEta[10000];   //[NobjTow]
   Float_t         TowPhi[10000];   //[NobjTow]
   Float_t         TowE[10000];   //[NobjTow]
   Float_t         TowEm[10000];   //[NobjTow]
   Float_t         TowHad[10000];   //[NobjTow]
   Float_t         TowOE[10000];   //[NobjTow]
   Int_t           Tow_jetidx[10000];   //[NobjTow]
   Int_t           TrackId[200];   //[NobjTrack]
   Int_t           TrackTowId[200];   //[NobjTrack]
   Int_t           TrackTowIdPhi[200];   //[NobjTrack]
   Int_t           TrackTowIdEta[200];   //[NobjTrack]
   Float_t         TrackPt[200];   //[NobjTrack]
   Float_t         TrackEta[200];   //[NobjTrack]
   Float_t         TrackEtaOut[200];   //[NobjTrack]
   Float_t         TrackPhi[200];   //[NobjTrack]
   Float_t         TrackPhiOut[200];   //[NobjTrack]
   Float_t         TrackDR[200];   //[NobjTrack]
   Float_t         TrackDROut[200];   //[NobjTrack]
   Float_t         TrackP[200];   //[NobjTrack]
   Float_t         TrackEMC1[200];   //[NobjTrack]
   Float_t         TrackEMC3[200];   //[NobjTrack]
   Float_t         TrackEMC5[200];   //[NobjTrack]
   Float_t         TrackHAC1[200];   //[NobjTrack]
   Float_t         TrackHAC3[200];   //[NobjTrack]
   Float_t         TrackHAC5[200];   //[NobjTrack]
   Float_t         TrackChi2[200];   //[NobjTrack]
   Int_t           TrackNHits[200];   //[NobjTrack]
   Float_t         MuDE[200];   //[NobjTrack]
   Float_t         MuDR[200];   //[NobjTrack]
   Int_t           Track_jetidx[10000];   //[NobjTow]
   Int_t           NobjJet;
   Float_t         JetPt[100];   //[NobjJet]
   Float_t         JetPhi[100];   //[NobjJet]
   Float_t         JetEta[100];   //[NobjJet]
   Float_t         JetEt[100];   //[NobjJet]
   Float_t         JetE[100];   //[NobjJet]
   Float_t         Met;
   Float_t         MetPhi;
   Float_t         MetSum;
   Float_t         Weight;

   // List of branches
   TBranch        *b_NobjTow;   //!
   TBranch        *b_NobjTrack;   //!
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
   TBranch        *b_TrackId;   //
   TBranch        *b_TrackTowId;   //
   TBranch        *b_TrackTowIdPhi;   //
   TBranch        *b_TrackTowIdEta;   //
   TBranch        *b_TrackPt;   //
   TBranch        *b_TrackEta;   //
   TBranch        *b_TrackEtaOut;   //
   TBranch        *b_TrackPhi;   //
   TBranch        *b_TrackPhiOut;   //
   TBranch        *b_TrackDR;   //
   TBranch        *b_TrackDROut;   //
   TBranch        *b_TrackP;   //
   TBranch        *b_TrackEMC1;   //
   TBranch        *b_TrackEMC3;   //
   TBranch        *b_TrackEMC5;   //
   TBranch        *b_TrackHAC1;   //
   TBranch        *b_TrackHAC3;   //
   TBranch        *b_TrackHAC5;   //
   TBranch        *b_TrackChi2;   //
   TBranch        *b_TrackNHits;   //
   TBranch        *b_MuDE;   //
   TBranch        *b_MuDR;   //
   TBranch        *b_Track_jetidx;   //!
   TBranch        *b_NobjJet;   //!
   TBranch        *b_JetPt;   //!
   TBranch        *b_JetPhi;   //!
   TBranch        *b_JetEta;   //!
   TBranch        *b_JetEt;   //!
   TBranch        *b_JetE;   //!
   TBranch        *b_Met;   //!
   TBranch        *b_MetPhi;   //!
   TBranch        *b_MetSum;   //!
   TBranch        *b_Weight;   //!

   NJetSel(TTree * /*tree*/ =0) { }
   virtual ~NJetSel() { }
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

   //ClassDef(NJetSel,0);
};

#endif

#ifdef NJetSel_cxx
void NJetSel::Init(TTree *tree)
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

   fChain->SetBranchAddress("NobjTow", &NobjTow, &b_NobjTow);
   fChain->SetBranchAddress("NobjTrack", &NobjTrack, &b_NobjTrack);
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
   fChain->SetBranchAddress("TrackId", TrackId, &b_TrackId);
   fChain->SetBranchAddress("TrackTowId", TrackTowId, &b_TrackTowId);
   fChain->SetBranchAddress("TrackTowIdPhi", TrackTowIdPhi, &b_TrackTowIdPhi);
   fChain->SetBranchAddress("TrackTowIdEta", TrackTowIdEta, &b_TrackTowIdEta);
   fChain->SetBranchAddress("TrackPt", TrackPt, &b_TrackPt);
   fChain->SetBranchAddress("TrackEta", TrackEta, &b_TrackEta);
   fChain->SetBranchAddress("TrackEtaOut", TrackEtaOut, &b_TrackEtaOut);
   fChain->SetBranchAddress("TrackPhi", TrackPhi, &b_TrackPhi);
   fChain->SetBranchAddress("TrackPhiOut", TrackPhiOut, &b_TrackPhiOut);
   fChain->SetBranchAddress("TrackDR", TrackDR, &b_TrackDR);
   fChain->SetBranchAddress("TrackDROut", TrackDROut, &b_TrackDROut);
   fChain->SetBranchAddress("TrackP", TrackP, &b_TrackP);
   fChain->SetBranchAddress("TrackEMC1", TrackEMC1, &b_TrackEMC1);
   fChain->SetBranchAddress("TrackEMC3", TrackEMC3, &b_TrackEMC3);
   fChain->SetBranchAddress("TrackEMC5", TrackEMC5, &b_TrackEMC5);
   fChain->SetBranchAddress("TrackHAC1", TrackHAC1, &b_TrackHAC1);
   fChain->SetBranchAddress("TrackHAC3", TrackHAC3, &b_TrackHAC3);
   fChain->SetBranchAddress("TrackHAC5", TrackHAC5, &b_TrackHAC5);
   fChain->SetBranchAddress("TrackChi2", TrackChi2, &b_TrackChi2);
   fChain->SetBranchAddress("TrackNHits", TrackNHits, &b_TrackNHits);
   fChain->SetBranchAddress("MuDR", MuDR, &b_MuDR);
   fChain->SetBranchAddress("MuDE", MuDE, &b_MuDE);
   fChain->SetBranchAddress("Track_jetidx", Track_jetidx, &b_Track_jetidx);
   fChain->SetBranchAddress("NobjJet", &NobjJet, &b_NobjJet);
   fChain->SetBranchAddress("JetPt", JetPt, &b_JetPt);
   fChain->SetBranchAddress("JetPhi", JetPhi, &b_JetPhi);
   fChain->SetBranchAddress("JetEta", JetEta, &b_JetEta);
   fChain->SetBranchAddress("JetEt", JetEt, &b_JetEt);
   fChain->SetBranchAddress("JetE", JetE, &b_JetE);
   fChain->SetBranchAddress("Met", &Met, &b_Met);
   fChain->SetBranchAddress("MetPhi", &MetPhi, &b_MetPhi);
   fChain->SetBranchAddress("MetSum", &MetSum, &b_MetSum);
   fChain->SetBranchAddress("Weight", &Weight, &b_Weight);
}

Bool_t NJetSel::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

#endif // #ifdef NJetSel_cxx
