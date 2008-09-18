#include <vector>
#include <iostream>

#include "TLorentzVector.h"
#include "TTree.h"
#include "TFile.h"
#include "TROOT.h"
#include "TVirtualFitter.h"
#include "TMinuit.h"
#include "Minuit2/FCNBase.h"
#include "Minuit2/FunctionMinimum.h"
#include "Minuit2/MnMigrad.h"
#include "Minuit2/MnHesse.h"
#include "Minuit2/MnUserParameters.h"
#include "Minuit2/MnPrint.h"
#include "Minuit2/MnContours.h"
#include "Minuit2/MnPlot.h"
#include "Minuit2/MnParameterScan.h"
#include "TH2.h"
#include "TStyle.h"

class Event;
std::vector<Event*> gData;
double gPar[] = { 1.2 , 1.0 , 1.0 , 1.0 , 1.0 };


TLorentzVector correctedTower(const TLorentzVector& v) {
  if(v.E() < 5) return v * gPar[0];
  if(v.E() < 10) return v * gPar[1];
  if(v.E() < 25) return v * gPar[2];
  if(v.E() < 50) return v * gPar[3];
  return v * gPar[4];
}


class Event {
public:
  double truept;
  std::vector<TLorentzVector> towers; 
  
  Event(double tp, int ntower) : truept(tp),towers(ntower) {}
  double chi2() {
    TLorentzVector jet(0,0,0,0); 
    double sigma2 = 0;
    for(unsigned int i = 0 ; i < towers.size() ; ++i) {
      TLorentzVector ct(correctedTower(towers[i]));
      if(ct.E() > 0) {
	jet += ct;
	sigma2+= (1.3 * 1.3 / ct.E() + 0.056 * 0.056) * ct.Pt() * ct.Pt();
      }
    }
    return (truept - jet.Pt()) * (truept - jet.Pt()) / sigma2;
  }
};

int fillData() {
   TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/scratch/current/cms/user/stadie/toy_photonjet2.root");
   if (!f) {
      f = new TFile("/scratch/current/cms/user/stadie/toy_photonjet2.root");
   }
   TTree *GammaJetTree = (TTree*)gDirectory->Get("GammaJetTree");

//Declaration of leaves types
   Int_t           NobjTowCal;
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
   Float_t         TowEmTrue[1000];
   Float_t         TowHadTrue[1000];
   Float_t         TowOETrue[1000];
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
   Float_t         PhtonEt;
   Float_t         PhotonE;
   Float_t         Weight;

   // Set branch addresses.
   GammaJetTree->SetBranchAddress("NobjTowCal",&NobjTowCal);
   GammaJetTree->SetBranchAddress("TowId",TowId);
   GammaJetTree->SetBranchAddress("TowId_phi",TowId_phi);
   GammaJetTree->SetBranchAddress("TowId_eta",TowId_eta);
   GammaJetTree->SetBranchAddress("TowEt",TowEt);
   GammaJetTree->SetBranchAddress("TowEta",TowEta);
   GammaJetTree->SetBranchAddress("TowPhi",TowPhi);
   GammaJetTree->SetBranchAddress("TowE",TowE);
   GammaJetTree->SetBranchAddress("TowEm",TowEm);
   GammaJetTree->SetBranchAddress("TowHad",TowHad);
   GammaJetTree->SetBranchAddress("TowOE",TowOE);
   GammaJetTree->SetBranchAddress("TowEmTrue",TowEmTrue);
   GammaJetTree->SetBranchAddress("TowHadTrue",TowHadTrue);
   GammaJetTree->SetBranchAddress("TowOETrue",TowOETrue);
   GammaJetTree->SetBranchAddress("JetCalPt",&JetCalPt);
   GammaJetTree->SetBranchAddress("JetCalPhi",&JetCalPhi);
   GammaJetTree->SetBranchAddress("JetCalEta",&JetCalEta);
   GammaJetTree->SetBranchAddress("JetCalEt",&JetCalEt);
   GammaJetTree->SetBranchAddress("JetCalE",&JetCalE);
   GammaJetTree->SetBranchAddress("JetGenPt",&JetGenPt);
   GammaJetTree->SetBranchAddress("JetGenPhi",&JetGenPhi);
   GammaJetTree->SetBranchAddress("JetGenEta",&JetGenEta);
   GammaJetTree->SetBranchAddress("JetGenEt",&JetGenEt);
   GammaJetTree->SetBranchAddress("JetGenE",&JetGenE);
   GammaJetTree->SetBranchAddress("MetCal",&MetCal);
   GammaJetTree->SetBranchAddress("MetCalPhi",&MetCalPhi);
   GammaJetTree->SetBranchAddress("MetCalSum",&MetCalSum);
   GammaJetTree->SetBranchAddress("PhotonPt",&PhotonPt);
   GammaJetTree->SetBranchAddress("PhotonPhi",&PhotonPhi);
   GammaJetTree->SetBranchAddress("PhotonEta",&PhotonEta);
   GammaJetTree->SetBranchAddress("PhotonEt",&PhtonEt);
   GammaJetTree->SetBranchAddress("PhotonE",&PhotonE);
   GammaJetTree->SetBranchAddress("Weight",&Weight);
   // GammaJetTree->SetBranchStatus("*",0);  // disable all branches
   // TTreePlayer->SetBranchStatus("branchname",1);  // activate branchname

   Long64_t nentries = GammaJetTree->GetEntries();
   //nentries = 1000;
   gData.reserve(nentries);
   for(Long64_t i = 0 ; i < nentries; ++i) {
     GammaJetTree->GetEntry(i);
     Event* ev = new Event(PhotonPt,NobjTowCal);
     TLorentzVector emjet(0,0,0,0);
     for(int j = 0 ; j < NobjTowCal ; ++j) {
       ev->towers[j].SetPtEtaPhiM(TowEt[j],TowEta[j],TowPhi[j],0);
       //ev->towers[j] *= ((TowEmTrue[j] + TowHadTrue[j])/ TowE[j]);
       double emf = TowEm[j]/TowE[j];
       emjet += emf * ev->towers[j];
       ev->towers[j] *= (1 - emf);
     }
     ev->truept -= emjet.Pt();
     //if(ev->chi2() > 1) GammaJetTree->Show(i);
     gData.push_back(ev);

   }
   return gData.size();
}

void clearData() {
  for(unsigned int i = 0 ; i < gData.size() ; ++i) {
    delete gData[i];
  }
  gData.clear();
}


class ParFCN : public ROOT::Minuit2::FCNBase {
public:
  
  ParFCN() {}
  ~ParFCN() {}

  double operator()(const std::vector<double>& par) const {
    for(unsigned int i = 0 ; i < par.size() ; ++i) {
      gPar[i] = par[i];
    }
    double chi2sum = 0;
    for(unsigned int i = 0 ; i < gData.size() ; ++i) {
      double c = gData[i]->chi2();
      chi2sum += c;
    }
    return chi2sum;
  }
  double Up() const {return 1.;}
  
private:
  
};


void fitData2() {
  using namespace ROOT::Minuit2;
  
  ParFCN fcn;
  MnUserParameters upar;
  upar.Add("par0", 0.8, 0.1);
  upar.Add("par1", 1.2, 0.1);
  upar.Add("par2", 1.3, 0.1);
  upar.Add("par3", 1.1, 0.1);
  upar.Add("par4", 0.9, 0.1);

  MnMigrad migrad(fcn, upar);
  FunctionMinimum min = migrad();
  std::cout<<"minimum: "<<min<<std::endl;  
 
  return;
  MnParameterScan parscan(fcn,min.UserParameters());  
  MnPlot plot;
  for(unsigned int i = 0 ; i < upar.VariableParameters() ; ++i) {
    std::cout << "scan for par" << i << "\n";
    std::vector<std::pair<double,double> > scan = parscan(i);
    plot(scan);
  }
  MnContours contours(fcn,min);
  std::cout << "contour: par0:par1\n";
  fcn.SetErrorDef(2.42);
  std::vector<std::pair<double,double> > cont = contours(0,1,40);
  plot(cont);

  // try to run hesse 
  //MnHesse hesse; 
  //hesse(fcn, min); 
  //std::cout<<"minimum after hesse: "<<min<<std::endl;
}

void make2DHist() {
  TH2F *h = new TH2F("h","",20,0.9,1.1,20,0.9,1.1);
  
  ParFCN fcn;
  std::vector<double> params;
  params.push_back(1.0);
  params.push_back(1.0);
  double step = 0.2/20;
  for(params[0] = 0.8 + step/2 ; params[0] < 1.2 ; params[0] += step) {
    for(params[1] = 0.8 + step/2 ; params[1] < 1.2 ; params[1] += step) {
      h->Fill(params[0],params[1],fcn(params));
    }
  }
  gStyle->SetPalette(1);
  h->DrawCopy("CONT1");
}



void simpleFit() {
  std::cout << "read " << fillData() << " events.\n";
  //for(unsigned int i = 0 ; i < gData.size() ; ++i) {
  //  std::cout << "chi2:" << gData[i]->chi2() << "\n";
  //}
  //fitData();
  //make2DHist();
  fitData2();
  clearData();
  
}
