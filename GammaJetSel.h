//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Aug 14 17:47:13 2007 by ROOT version 5.14/00b
// from TTree CalibTree/
// found on file: calib.root
//////////////////////////////////////////////////////////

#ifndef GammaJetSel_h
#define GammaJetSel_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>

#include <iostream>
class GammaJetSel : public TSelector {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain

   // Declaration of leave types
   Int_t           NobjTowCal;
   Int_t           TowId[200];   //[NobjTowCal]
   Int_t           TowId_phi[200];   //[NobjTowCal]
   Int_t           TowId_eta[200];   //[NobjTowCal]
   Float_t         TowEt[200];   //[NobjTowCal]
   Float_t         TowEta[200];   //[NobjTowCal]
   Float_t         TowPhi[200];   //[NobjTowCal]
   Float_t         TowE[200];   //[NobjTowCal]
   Float_t         TowEm[200];   //[NobjTowCal]
   Float_t         TowHad[200];   //[NobjTowCal]
   Float_t         TowOE[200];   //[NobjTowCal]
   Float_t         JetCalPt;
   Float_t         JetCalPhi;
   Float_t         JetCalEta;
   Float_t         JetCalEt;
   Float_t         JetCalE; 
   Float_t         JetGenPt;
   Float_t         JetGenPhi;
   Float_t         JetGenEta;
   Float_t         JetGenEt;
   Float_t         JetGenE;
   Float_t         MetCal;
   Float_t         MetCalPhi;
   Float_t         MetCalSum;
   Float_t         PhotonPt;
   Float_t         PhotonPhi;
   Float_t         PhotonEta;
   Float_t         PhotonEt;
   Float_t         PhotonE;
   Float_t         GenPhotonPt;
   Float_t         GenPhotonPhi;
   Float_t         GenPhotonEta;
   Float_t         GenPhotonEt;
   Float_t         GenPhotonE;
   //Int_t           ProcessID;
   Float_t         EventWeight;   

   // List of branches
   TBranch        *b_NobjTowCal;   //!
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
   TBranch        *b_JetCalPt;   //!
   TBranch        *b_JetCalPhi;   //!
   TBranch        *b_JetCalEta;   //!
   TBranch        *b_JetCalEt;   //!
   TBranch        *b_JetCalE;   //!  
   TBranch        *b_JetGenPt;   //!
   TBranch        *b_JetGenPhi;   //!
   TBranch        *b_JetGenEta;   //!
   TBranch        *b_JetGenEt;   //!
   TBranch        *b_JetGenE;   //!
   TBranch        *b_MetCal;   //!
   TBranch        *b_MetCalPhi;   //!
   TBranch        *b_MetCalSum;   //!
   TBranch        *b_PhotonPt;   //!
   TBranch        *b_PhotonPhi;   //!
   TBranch        *b_PhotonEta;   //!
   TBranch        *b_PhtonEt;   //!
   TBranch        *b_PhotonE;   //!
   TBranch        *b_GenPhotonPt;   //!
   TBranch        *b_GenPhotonPhi;   //!
   TBranch        *b_GenPhotonEta;   //!
   TBranch        *b_GenPhtonEt;   //!
   TBranch        *b_GenPhotonE;   //!
   //TBranch        *b_ProcessID;   //!
   TBranch        *b_EventWeight;   //!

   //GammaJetSel(TTree * /*tree*/ =0) : ProcessID(-1), EventWeight(1) { }
   GammaJetSel(TTree * /*tree*/ =0) : EventWeight(1) { }
   virtual ~GammaJetSel() { }
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

   //   ClassDef(GammaJetSel,0);
};

#endif

#ifdef GammaJetSel_cxx
void GammaJetSel::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normaly not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).
   // Set branch addresses and branch pointers
   if (!tree) { 
     std::cerr<<"tree is bad!"<<std::endl;
     return;
   }
   fChain = tree;
   fChain->SetMakeClass(1);
   fChain->SetBranchAddress("NobjTowCal", &NobjTowCal, &b_NobjTowCal);
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
   fChain->SetBranchAddress("JetCalPt", &JetCalPt, &b_JetCalPt);
   fChain->SetBranchAddress("JetCalPhi", &JetCalPhi, &b_JetCalPhi);
   fChain->SetBranchAddress("JetCalEta", &JetCalEta, &b_JetCalEta);
   fChain->SetBranchAddress("JetCalEt", &JetCalEt, &b_JetCalEt);
   fChain->SetBranchAddress("JetCalE", &JetCalE, &b_JetCalE); 
   fChain->SetBranchAddress("JetGenPt", &JetGenPt, &b_JetGenPt);
   fChain->SetBranchAddress("JetGenPhi", &JetGenPhi, &b_JetGenPhi);
   fChain->SetBranchAddress("JetGenEta", &JetGenEta, &b_JetGenEta);
   fChain->SetBranchAddress("JetGenEt", &JetGenEt, &b_JetGenEt);
   fChain->SetBranchAddress("JetGenE", &JetGenE, &b_JetGenE);
   fChain->SetBranchAddress("MetCal", &MetCal, &b_MetCal);
   fChain->SetBranchAddress("MetCalPhi", &MetCalPhi, &b_MetCalPhi);
   fChain->SetBranchAddress("MetCalSum", &MetCalSum, &b_MetCalSum);
   fChain->SetBranchAddress("PhotonPt", &PhotonPt, &b_PhotonPt);
   fChain->SetBranchAddress("PhotonPhi", &PhotonPhi, &b_PhotonPhi);
   fChain->SetBranchAddress("PhotonEta", &PhotonEta, &b_PhotonEta);
   fChain->SetBranchAddress("PhotonEt", &PhotonEt, &b_PhtonEt);
   fChain->SetBranchAddress("PhotonE", &PhotonE, &b_PhotonE);
   fChain->SetBranchAddress("GenPhotonPt", &GenPhotonPt, &b_GenPhotonPt);
   fChain->SetBranchAddress("GenPhotonPhi", &GenPhotonPhi, &b_GenPhotonPhi);
   fChain->SetBranchAddress("GenPhotonEta", &GenPhotonEta, &b_GenPhotonEta);
   fChain->SetBranchAddress("GenPhotonEt", &GenPhotonEt, &b_GenPhtonEt);
   fChain->SetBranchAddress("GenPhotonE", &GenPhotonE, &b_GenPhotonE);
   //fChain->SetBranchAddress("ProcessID", &ProcessID, &b_ProcessID);
   fChain->SetBranchAddress("EventWeight", &EventWeight, &b_EventWeight);
}

Bool_t GammaJetSel::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normaly not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

#endif // #ifdef GammaJetSel_cxx
