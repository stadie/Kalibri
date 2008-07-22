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
