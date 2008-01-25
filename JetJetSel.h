//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Jan 25 11:05:05 2008 by ROOT version 5.14/00
// from TTree JetJetTree/
// found on file: susyincl400400test.jetjet.root
//////////////////////////////////////////////////////////

#ifndef JetJetSel_h
#define JetJetSel_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>

class JetJetSel : public TSelector {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain

   // Declaration of leave types
   Int_t           NobjTowJ1Cal;
   Int_t           TowJ1Id[200];   //[NobjTowJ1Cal]
   Int_t           TowJ1Id_phi[200];   //[NobjTowJ1Cal]
   Int_t           TowJ1Id_eta[200];   //[NobjTowJ1Cal]
   Float_t         TowJ1Et[200];   //[NobjTowJ1Cal]
   Float_t         TowJ1Eta[200];   //[NobjTowJ1Cal]
   Float_t         TowJ1Phi[200];   //[NobjTowJ1Cal]
   Float_t         TowJ1E[200];   //[NobjTowJ1Cal]
   Float_t         TowJ1Em[200];   //[NobjTowJ1Cal]
   Float_t         TowJ1Had[200];   //[NobjTowJ1Cal]
   Float_t         TowJ1OE[200];   //[NobjTowJ1Cal]
   Int_t           NobjTowJ2Cal;
   Int_t           TowJ2Id[200];   //[NobjTowJ2Cal]
   Int_t           TowJ2Id_phi[200];   //[NobjTowJ2Cal]
   Int_t           TowJ2Id_eta[200];   //[NobjTowJ2Cal]
   Float_t         TowJ2Et[200];   //[NobjTowJ2Cal]
   Float_t         TowJ2Eta[200];   //[NobjTowJ2Cal]
   Float_t         TowJ2Phi[200];   //[NobjTowJ2Cal]
   Float_t         TowJ2E[200];   //[NobjTowJ2Cal]
   Float_t         TowJ2Em[200];   //[NobjTowJ2Cal]
   Float_t         TowJ2Had[200];   //[NobjTowJ2Cal]
   Float_t         TowJ2OE[200];   //[NobjTowJ2Cal]
   Float_t         FirstJetPt;
   Float_t         FirstJetPhi;
   Float_t         FirstJetEta;
   Float_t         FirstJetEt;
   Float_t         FirstJetE;
   Float_t         ScndJetPt;
   Float_t         ScndJetPhi;
   Float_t         ScndJetEta;
   Float_t         ScndJetEt;
   Float_t         ScndJetE;
   Float_t         FirstJetGenPt;
   Float_t         FirstJetGenPhi;
   Float_t         FirstJetGenEta;
   Float_t         FirstJetGenEt;
   Float_t         FirstJetGenE;
   Float_t         ScndJetGenPt;
   Float_t         ScndJetGenPhi;
   Float_t         ScndJetGenEta;
   Float_t         ScndJetGenEt;
   Float_t         ScndJetGenE;
   Float_t         MetCal;
   Float_t         MetCalPhi;
   Float_t         MetCalSum;

   // List of branches
   TBranch        *b_NobjTowJ1Cal;   //!
   TBranch        *b_TowJ1Id;   //!
   TBranch        *b_TowJ1Id_phi;   //!
   TBranch        *b_TowJ1Id_eta;   //!
   TBranch        *b_TowJ1Et;   //!
   TBranch        *b_TowJ1Eta;   //!
   TBranch        *b_TowJ1Phi;   //!
   TBranch        *b_TowJ1E;   //!
   TBranch        *b_TowJ1Em;   //!
   TBranch        *b_TowJ1Had;   //!
   TBranch        *b_TowJ1OE;   //!
   TBranch        *b_NobjTowJ2Cal;   //!
   TBranch        *b_TowJ2Id;   //!
   TBranch        *b_TowJ2Id_phi;   //!
   TBranch        *b_TowJ2Id_eta;   //!
   TBranch        *b_TowJ2Et;   //!
   TBranch        *b_TowJ2Eta;   //!
   TBranch        *b_TowJ2Phi;   //!
   TBranch        *b_TowJ2E;   //!
   TBranch        *b_TowJ2Em;   //!
   TBranch        *b_TowJ2Had;   //!
   TBranch        *b_TowJ2OE;   //!
   TBranch        *b_FirstJetPt;   //!
   TBranch        *b_FirstJetPhi;   //!
   TBranch        *b_FirstJetEta;   //!
   TBranch        *b_FirstJetEt;   //!
   TBranch        *b_FirstJetE;   //!
   TBranch        *b_ScndJetPt;   //!
   TBranch        *b_ScndJetPhi;   //!
   TBranch        *b_ScndJetEta;   //!
   TBranch        *b_ScndJetEt;   //!
   TBranch        *b_ScndJetE;   //!
   TBranch        *b_FirstJetGenPt;   //!
   TBranch        *b_FirstJetGenPhi;   //!
   TBranch        *b_FirstJetGenEta;   //!
   TBranch        *b_FirstJetGenEt;   //!
   TBranch        *b_FirstJetGenE;   //!
   TBranch        *b_ScndJetGenPt;   //!
   TBranch        *b_ScndJetGenPhi;   //!
   TBranch        *b_ScndJetGenEta;   //!
   TBranch        *b_ScndJetGenEt;   //!
   TBranch        *b_ScndJetGenE;   //!
   TBranch        *b_MetCal;   //!
   TBranch        *b_MetCalPhi;   //!
   TBranch        *b_MetCalSum;   //!

   JetJetSel(TTree * /*tree*/ =0) { }
   virtual ~JetJetSel() { }
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

   //ClassDef(JetJetSel,0);
};

#endif

#ifdef JetJetSel_cxx
void JetJetSel::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normaly not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("NobjTowJ1Cal", &NobjTowJ1Cal, &b_NobjTowJ1Cal);
   fChain->SetBranchAddress("TowJ1Id", TowJ1Id, &b_TowJ1Id);
   fChain->SetBranchAddress("TowJ1Id_phi", TowJ1Id_phi, &b_TowJ1Id_phi);
   fChain->SetBranchAddress("TowJ1Id_eta", TowJ1Id_eta, &b_TowJ1Id_eta);
   fChain->SetBranchAddress("TowJ1Et", TowJ1Et, &b_TowJ1Et);
   fChain->SetBranchAddress("TowJ1Eta", TowJ1Eta, &b_TowJ1Eta);
   fChain->SetBranchAddress("TowJ1Phi", TowJ1Phi, &b_TowJ1Phi);
   fChain->SetBranchAddress("TowJ1E", TowJ1E, &b_TowJ1E);
   fChain->SetBranchAddress("TowJ1Em", TowJ1Em, &b_TowJ1Em);
   fChain->SetBranchAddress("TowJ1Had", TowJ1Had, &b_TowJ1Had);
   fChain->SetBranchAddress("TowJ1OE", TowJ1OE, &b_TowJ1OE);
   fChain->SetBranchAddress("NobjTowJ2Cal", &NobjTowJ2Cal, &b_NobjTowJ2Cal);
   fChain->SetBranchAddress("TowJ2Id", TowJ2Id, &b_TowJ2Id);
   fChain->SetBranchAddress("TowJ2Id_phi", TowJ2Id_phi, &b_TowJ2Id_phi);
   fChain->SetBranchAddress("TowJ2Id_eta", TowJ2Id_eta, &b_TowJ2Id_eta);
   fChain->SetBranchAddress("TowJ2Et", TowJ2Et, &b_TowJ2Et);
   fChain->SetBranchAddress("TowJ2Eta", TowJ2Eta, &b_TowJ2Eta);
   fChain->SetBranchAddress("TowJ2Phi", TowJ2Phi, &b_TowJ2Phi);
   fChain->SetBranchAddress("TowJ2E", TowJ2E, &b_TowJ2E);
   fChain->SetBranchAddress("TowJ2Em", TowJ2Em, &b_TowJ2Em);
   fChain->SetBranchAddress("TowJ2Had", TowJ2Had, &b_TowJ2Had);
   fChain->SetBranchAddress("TowJ2OE", TowJ2OE, &b_TowJ2OE);
   fChain->SetBranchAddress("FirstJetPt", &FirstJetPt, &b_FirstJetPt);
   fChain->SetBranchAddress("FirstJetPhi", &FirstJetPhi, &b_FirstJetPhi);
   fChain->SetBranchAddress("FirstJetEta", &FirstJetEta, &b_FirstJetEta);
   fChain->SetBranchAddress("FirstJetEt", &FirstJetEt, &b_FirstJetEt);
   fChain->SetBranchAddress("FirstJetE", &FirstJetE, &b_FirstJetE);
   fChain->SetBranchAddress("ScndJetPt", &ScndJetPt, &b_ScndJetPt);
   fChain->SetBranchAddress("ScndJetPhi", &ScndJetPhi, &b_ScndJetPhi);
   fChain->SetBranchAddress("ScndJetEta", &ScndJetEta, &b_ScndJetEta);
   fChain->SetBranchAddress("ScndJetEt", &ScndJetEt, &b_ScndJetEt);
   fChain->SetBranchAddress("ScndJetE", &ScndJetE, &b_ScndJetE);
   fChain->SetBranchAddress("FirstJetGenPt", &FirstJetGenPt, &b_FirstJetGenPt);
   fChain->SetBranchAddress("FirstJetGenPhi", &FirstJetGenPhi, &b_FirstJetGenPhi);
   fChain->SetBranchAddress("FirstJetGenEta", &FirstJetGenEta, &b_FirstJetGenEta);
   fChain->SetBranchAddress("FirstJetGenEt", &FirstJetGenEt, &b_FirstJetGenEt);
   fChain->SetBranchAddress("FirstJetGenE", &FirstJetGenE, &b_FirstJetGenE);
   fChain->SetBranchAddress("ScndJetGenPt", &ScndJetGenPt, &b_ScndJetGenPt);
   fChain->SetBranchAddress("ScndJetGenPhi", &ScndJetGenPhi, &b_ScndJetGenPhi);
   fChain->SetBranchAddress("ScndJetGenEta", &ScndJetGenEta, &b_ScndJetGenEta);
   fChain->SetBranchAddress("ScndJetGenEt", &ScndJetGenEt, &b_ScndJetGenEt);
   fChain->SetBranchAddress("ScndJetGenE", &ScndJetGenE, &b_ScndJetGenE);
   fChain->SetBranchAddress("MetCal", &MetCal, &b_MetCal);
   fChain->SetBranchAddress("MetCalPhi", &MetCalPhi, &b_MetCalPhi);
   fChain->SetBranchAddress("MetCalSum", &MetCalSum, &b_MetCalSum);
}

Bool_t JetJetSel::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normaly not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

#endif // #ifdef JetJetSel_cxx
