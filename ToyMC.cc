//
// Original Author:  Hartmut Stadie
//         Created:  Mon Jun 30 11:00:00 CEST 2008
// $Id: ToyMC.cc,v 1.3 2008/06/30 13:15:27 stadie Exp $
//
#include "ToyMC.h"

#include <cmath>
#include <iostream>
#include <map>
#include <cassert> 

#include "TRandom.h"
#include "TRandom3.h"
#include "TTree.h"
#include "TFile.h"



ToyMC::ToyMC() : mMinEta(-2.5),mMaxEta(2.5),mMinPt(30), mMaxPt(400), mTowConst(1.25),mResoStochastic(1.20),mResoNoise(0.05),mJetSpread(0.07),mNoOutOfCone(true),mModel(Gauss),mChunks(200)
{
  mRandom = new TRandom3();
  mRandom->SetSeed(0);
}

void ToyMC::genInput() { 
  static double rand[3];

  mRandom->RndmArray(3,rand);
  mPinput.SetPtEtaPhiM(rand[0]*(mMaxPt - mMinPt)+mMinPt,
		       rand[1]*(mMaxEta - mMinEta)+mMinEta,
		       rand[2]*2 * M_PI - M_PI, 0);
}

void ToyMC::calIds(float& eta, float &phi, int& ieta, int& iphi) 
{
  const float dEta = 2.5/30;
  const float dPhi =  2 * M_PI / 72;
  
  if(eta > 0) {
    ieta = (int)(eta / dEta) + 1;
  } else {
    ieta = -((int)(std::abs(eta) / dEta) + 1);
  }
  if(ieta > 0) {
    eta = (ieta - 0.5) * dEta;
  } else {
    eta = (ieta + 0.5) * dEta;
  }
  
  if(phi < 0) phi += 2 * M_PI;
  iphi = (int)(phi / dPhi);
  iphi++;
  if(iphi > 72) iphi = 72; 
  phi = (iphi-0.5) * dPhi;
  assert(phi < 2 * M_PI);
  assert(ieta != 0);
  assert(iphi > 0);
  assert(iphi <= 72);
  assert(ieta <= 40);
  assert(ieta >= -40);
}

void ToyMC::smearTower(double e, float& te, float& tem, float& thad, float& tout) {
  float emf = mRandom->Uniform(0.5);
  tem = emf * e;
  thad = (1-emf) * e / mTowConst;
  tout = 0;
  if(mModel == Gauss) {
    thad = mRandom->Gaus(1.0,sqrt(mResoStochastic * mResoStochastic/ thad + 
				  mResoNoise * mResoNoise)) * thad;
  }  else if(mModel == Landau) {
    double smear;
    do {
      smear =  mRandom->Landau(1,sqrt(mResoStochastic * mResoStochastic/ thad + mResoNoise * mResoNoise));
    } while((smear < 0) || (smear > 2));
    smear = 2 - smear;
    thad *= smear;
  }
  if(thad < 0) thad = 0;
  te = tem + thad + tout;
}

int ToyMC::splitJet(const TLorentzVector& jet ,float* et,float* eta,float * phi, int* ieta,int* iphi) {
  typedef std::map<int,int> TowerMap;

  TowerMap towers;
  double jphi = jet.Phi();
  if(jphi < 0) jphi += 2 * M_PI;
  //std::cout << "jet: Pt:" << jet.Pt() << " Phi:" << jet.Phi() << " Eta:" << jet.Eta() << '\n';
  double de = jet.E() / mChunks;
  int ntowers = 0;
  TLorentzVector rec(0,0,0,0);
  TLorentzVector tow;
  for(int i = 0 ; i < mChunks ; ++i) {
    float teta = mRandom->Gaus(jet.Eta(), mJetSpread);
    float tphi = mRandom->Gaus(jet.Phi(), mJetSpread);
    if( tphi < 0) tphi += 2 * M_PI;
    int ie, ip;
    calIds(teta, tphi, ie, ip); 
    //std::cout << "vorher:" << teta << ", " << tphi << ", " << ie << ", " 
    //	      << ip << "\n";
    float dphi = TVector2::Phi_mpi_pi(tphi-jphi);
    float deta = teta-jet.Eta();
    if(sqrt(deta*deta + dphi*dphi) > 0.5) {
      //std::cout << "Out of cone:" << teta << ":" << jet.Eta() << " , " << dphi << '\n';
      if(mNoOutOfCone) --i;
      continue;
    }
    int id = towers[ie*1000 + ip];
    if(! id) {
      ++ntowers;
      id = ntowers;
      towers[ie*1000 + ip] = id;
      et[id-1] = 0;
    }
    --id;
    tow.SetPtEtaPhiM(1,teta,tphi,0);
    tow *= de/tow.E();
    et[id] += tow.Pt();
    eta[id] = teta;
    phi[id] = tphi;
    ieta[id] = ie;
    iphi[id] = ip;
    rec += tow;
  }
  //std::cout  << "Eta:" << jet.Eta() <<  "       : " << rec.Pt() << "," << rec.E() << "  == " << jet.Pt() << "," << jet.E() << '\n';
  //std::cout << "lost energy:" << lostE/jet.E() << '\n';
  //assert(lostE/jet.E() < 0.25);
  return ntowers;
}

int ToyMC::generateTrackClusterTree(TTree* CalibTree, int nevents) 
{
  //make tree
  const int kMaxTower = 1;
  int NobjTowCal;
  float towet[kMaxTower];
  float toweta[kMaxTower];
  float towphi[kMaxTower];
  float towen[kMaxTower];
  float towem[kMaxTower];
  float towhd[kMaxTower];
  float towoe[kMaxTower];
  int towid_phi[kMaxTower];
  int towid_eta[kMaxTower];
  int towid [kMaxTower];
  float tracket;
  float tracketerr;
  float tracketa;
  float trackphi;
  float tracken;
  CalibTree->Branch("NobjTowCal",&NobjTowCal,"NobjTowCal/I");
  CalibTree->Branch("TowId",towid,"TowId[NobjTowCal]/I");
  CalibTree->Branch("TowId_phi",towid_phi,"TowId_phi[NobjTowCal]/I");
  CalibTree->Branch("TowId_eta",towid_eta,"TowId_eta[NobjTowCal]/I");
  CalibTree->Branch("TowEt",towet,"TowEt[NobjTowCal]/F");
  CalibTree->Branch("TowEta",toweta,"TowEta[NobjTowCal]/F");
  CalibTree->Branch("TowPhi",towphi,"TowPhi[NobjTowCal]/F");
  CalibTree->Branch("TowE",towen,"TowE[NobjTowCal]/F");
  CalibTree->Branch("TowEm",towem,"TowEm[NobjTowCal]/F");
  CalibTree->Branch("TowHad",towhd,"TowHad[NobjTowCal]/F");
  CalibTree->Branch("TowOE",towoe,"TowOE[NobjTowCal]/F");
  CalibTree->Branch("TrackEt",&tracket,"TrackEt/F");
  CalibTree->Branch("TrackEterr",&tracketerr,"TrackEterr/F");
  CalibTree->Branch("TrackEta",&tracketa,"TrackEta/F");
  CalibTree->Branch("TrackPhi",&trackphi,"TrackPhi/F");
  CalibTree->Branch("TrackE",&tracken,"TrackE/F");

  for(int i = 0; i < nevents ; ++i) {
    genInput();
    tracket = mPinput.Pt();
    tracketerr = 0;
    tracketa = mPinput.Eta();
    trackphi = mPinput.Phi();
    tracken = mPinput.E();
    NobjTowCal = 1;
    towphi[0] = mPinput.Phi(); 
    toweta[0] = mPinput.Eta();
    towid[0] = 0;
    //std::cout << "vorher:" << toweta[0] << ", " << towphi[0] << ", " << towid_eta[0] << ", " 
    //	      << towid_phi[0] << "\n";
    calIds(toweta[0],towphi[0],towid_eta[0],towid_phi[0]);
    //std::cout << "nachher:" << toweta[0] << ", " << towphi[0] << ", " << towid_eta[0] << ", " 
    //	      << towid_phi[0] << "\n";
    smearTower(mPinput.E(),towen[0],towem[0],towhd[0],towoe[0]);
    towet[0] = towen[0]/mPinput.E() * mPinput.Pt();
    CalibTree->Fill();
    if(i % 1000 == 0) std::cout << "generated event " << i << '\n';
  }
  return CalibTree->GetEntriesFast(); 
}

int ToyMC::makeTrackCluster(const char* filename, int nevents) {
  TFile* file = new TFile(filename,"recreate");
  TTree* CalibTree = new TTree("TrackClusterTree","TrackClusterTre");
 
  nevents = generateTrackClusterTree(CalibTree, nevents);
  file->Write();
  file->Close();
  return  nevents;
}

int ToyMC::generatePhotonJetTree(TTree* CalibTree, int nevents)
{
  //make tree 
  const int kMaxTower = 1000;
  int NobjTowCal;
  float towet[kMaxTower];
  float toweta[kMaxTower];
  float towphi[kMaxTower];
  float towen[kMaxTower];
  float towem[kMaxTower];
  float towhd[kMaxTower];
  float towoe[kMaxTower];
  int towid_phi[kMaxTower];
  int towid_eta[kMaxTower];
  int towid [kMaxTower];
  float jcalpt,jcalphi,jcaleta,jcalet,jcale;
  float jgenpt,jgenphi,jgeneta,jgenet,jgene;
  float mcalmet,mcalphi,mcalsum;
  float photonpt,photonphi,photoneta,photonet,photone;
  CalibTree->Branch("NobjTowCal",&NobjTowCal,"NobjTowCal/I");
  CalibTree->Branch("TowId",towid,"TowId[NobjTowCal]/I");
  CalibTree->Branch("TowId_phi",towid_phi,"TowId_phi[NobjTowCal]/I");
  CalibTree->Branch("TowId_eta",towid_eta,"TowId_eta[NobjTowCal]/I");
  CalibTree->Branch("TowEt",towet,"TowEt[NobjTowCal]/F");
  CalibTree->Branch("TowEta",toweta,"TowEta[NobjTowCal]/F");
  CalibTree->Branch("TowPhi",towphi,"TowPhi[NobjTowCal]/F");
  CalibTree->Branch("TowE",towen,"TowE[NobjTowCal]/F");
  CalibTree->Branch("TowEm",towem,"TowEm[NobjTowCal]/F");
  CalibTree->Branch("TowHad",towhd,"TowHad[NobjTowCal]/F");
  CalibTree->Branch("TowOE",towoe,"TowOE[NobjTowCal]/F");
  // Jet- MEt-specific branches of the tree
  CalibTree->Branch("JetCalPt",&jcalpt,"JetCalPt/F");
  CalibTree->Branch("JetCalPhi",&jcalphi,"JetCalPhi/F");
  CalibTree->Branch("JetCalEta",&jcaleta,"JetCalEta/F");
  CalibTree->Branch("JetCalEt",&jcalet,"JetCalEt/F");
  CalibTree->Branch("JetCalE",&jcale,"JetCalE/F");
  // Gen- Jet- branches of the tree
  CalibTree->Branch("JetGenPt",&jgenpt,"JetGenPt/F");
  CalibTree->Branch("JetGenPhi",&jgenphi,"JetGenPhi/F");
  CalibTree->Branch("JetGenEta",&jgeneta,"JetGenEta/F");
  CalibTree->Branch("JetGenEt",&jgenet,"JetGenEt/F");
  CalibTree->Branch("JetGenE",&jgene,"JetGenE/F");
  //met
  CalibTree->Branch("MetCal",&mcalmet,"MetCal/F");
  CalibTree->Branch("MetCalPhi",&mcalphi,"MetCalPhi/F");
  CalibTree->Branch("MetCalSum",&mcalsum,"MetCalSum/F");
  //photons
  CalibTree->Branch("PhotonPt",&photonpt,"PhotonPt/F");
  CalibTree->Branch("PhotonPhi",&photonphi,"PhotonPhi/F");
  CalibTree->Branch("PhotonEta",&photoneta,"PhotonEta/F");
  CalibTree->Branch("PhotonEt",&photonet,"PhtonEt/F");
  CalibTree->Branch("PhotonE",&photone,"PhotonE/F");
 
  TLorentzVector jet, tower;
  for(int i = 0; i < nevents ; ++i) {
    genInput();
    photonpt = mPinput.Pt();
    photoneta = mPinput.Eta();
    photonphi = mPinput.Phi();
    photonet = mPinput.Pt();
    photone = mPinput.E();
    jgenpt = mRandom->Gaus(photonpt,0.04 * photonpt);
    jgeneta = mRandom->Gaus(photoneta,1.0);
    if((jgeneta > 2.5) || (jgeneta < -2.5)) {
      --i;
      continue;
    }
    jgenphi = photonphi +M_PI;
    if(jgenphi > M_PI) jgenphi -= 2 * M_PI;
    jgenet = jgenpt;
    jet.SetPtEtaPhiM(jgenpt,jgeneta,jgenphi,0);
    jgene = jet.E();
    NobjTowCal = splitJet(jet,towet,toweta,towphi,towid_eta,towid_phi);
    jet.SetPtEtaPhiM(0,0,0,0);
    double towmean = mTowConst;
    if(mModel == Flat) mTowConst = 1/mRandom->Uniform(1.5);
    for(int j = 0; j < NobjTowCal ; ++j) {
      tower.SetPtEtaPhiM(towet[j],toweta[j],towphi[j],0);
      towen[j] = tower.E();
      smearTower(tower.E(),towen[j],towem[j],towhd[j],towoe[j]); 
      tower *= towen[j]/tower.E();
      towet[j] = tower.Pt();
      jet += tower;
    }
    mTowConst = towmean;
    jcalpt = jet.Pt();
    jcaleta = jet.Eta();
    jcalphi = jet.Phi();
    jcalet = jet.Pt();
    jcale = jet.E();
    mcalmet = 0;
    mcalphi = 0;
    mcalsum = 0;
    CalibTree->Fill();
    if(i % 1000 == 0) std::cout << "writing event " << i << '\n';
  } 
  return CalibTree->GetEntriesFast();
}

int ToyMC::makePhotonJet(const char* filename, int nevents) {
  TFile* file = new TFile(filename,"recreate");
  TTree* CalibTree = new TTree("GammaJetTree","GammaJetTree");

  nevents = generatePhotonJetTree(CalibTree,nevents);
  //file->ls();
  file->Write();
  file->Close();
  return nevents;
}
