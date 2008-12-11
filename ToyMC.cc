//
// Original Author:  Hartmut Stadie
//         Created:  Mon Jun 30 11:00:00 CEST 2008
// $Id: ToyMC.cc,v 1.15 2008/11/20 16:38:03 stadie Exp $
//
#include "ToyMC.h"

#include <cmath>
#include <iostream>
#include <map>
#include <cassert> 
#include <ext/hash_map>

#include "TRandom.h"
#include "TRandom3.h"
#include "TTree.h"
#include "TFile.h"

#include "ConfigFile.h"

ToyMC::ToyMC() : mMinEta(-2.5),mMaxEta(2.5),mMinPt(30), mMaxPt(400),mPtSpectrum(Uniform),
		 mResoStochastic(1.20),mResoNoise(0.05),mJetSpreadA(0.5),mJetSpreadB(0),
		 mNoOutOfCone(true),mModel(Gauss),mChunks(200),mMaxPi0Frac(0.5),mMaxEmf(0.5)
{
  mTowConst[0] = 1.24;
  mTowConst[1] = 0;
  mTowConst[2] = 1;
  mTowConst[3] = 1;
  mTowConst[4] = 0;
  mRandom = new TRandom3();
  mRandom->SetSeed(0);
}

void ToyMC::genInput() { 
  static double rand[3];
  mRandom->RndmArray(3,rand);
  double pt = 0;  
  if(mPtSpectrum == Uniform) {
    mRandom->RndmArray(3,rand);
    pt = rand[0]*(mMaxPt - mMinPt)+mMinPt;
  } else if(mPtSpectrum == PowerLaw) {
    pt = mMinPt * pow(rand[0],-1.0/3.5);
  }
  mPinput.SetPtEtaPhiM(pt,
		       rand[1]*(mMaxEta - mMinEta)+mMinEta,
		       rand[2]*2 * M_PI - M_PI, 0);
}

void ToyMC::calIds(float& eta, float &phi, int& ieta, int& iphi) 
{
  const static float dEta = 2.5/30;
  const static float dPhi =  2 * M_PI / 72;
  
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

void ToyMC::smearTower(double e, float& te, float& tem, float& thad, float& tout, 
		       float& temtrue, float& thadtrue, float& touttrue) 
{
  float emf = mRandom->Uniform(mMaxEmf);
  temtrue = emf * e;
  tem = temtrue;
  thadtrue = (1-emf) * e;
  if(mTowConst[1] == 0) {
    thadtrue /= mTowConst[0];
  } else { 
    double c =  (thadtrue < 1.0) ? mTowConst[0] - mTowConst[1]/mTowConst[3] +  mTowConst[4] :  mTowConst[0] - mTowConst[1]/(pow(log(thadtrue),mTowConst[2]) + mTowConst[3]) + mTowConst[4]/thadtrue;
    //std::cout << thadtrue << ",  " << c << '\n';
    thadtrue /= c;
  }
  thad = thadtrue;
  touttrue = 0;
  tout = touttrue;
  if(mModel == Landau) {
    double smear;
    do {
      smear =  mRandom->Landau(1,sqrt(mResoStochastic * mResoStochastic/ thad + mResoNoise * mResoNoise));
    } while((smear < 0) || (smear > 2));
    smear = 2 - smear;
    thad *= smear;
  } else {
    double smear;
    do {
      smear = mRandom->Gaus(1.0,sqrt(mResoStochastic * mResoStochastic/ thad + 
				     mResoNoise * mResoNoise));
    } while((smear < 0) || (smear > 2));
    thad *= smear;
  }
  if(thad < 0) thad = 0;
  te = tem + thad + tout;
}

int ToyMC::splitJet(const TLorentzVector& jet ,float* et,float* eta,float * phi, int* ieta,int* iphi) {
  typedef __gnu_cxx::hash_map<int,int> TowerMap;

  TowerMap towers;
  double jphi = jet.Phi();
  if(jphi < 0) jphi += 2 * M_PI;
  //std::cout << "jet: Pt:" << jet.Pt() << " Phi:" << jet.Phi() << " Eta:" << jet.Eta() << '\n';
  //double de = jet.E() / mChunks;
  int ntowers = 0;
  TLorentzVector rec(0,0,0,0);
  TLorentzVector tow;
  double lostPt = 0;
  double dpt = jet.Pt() / mChunks;
  if(mChunks < 0) {
    dpt = 0.3;
    mChunks = (int)std::ceil(jet.Pt() / dpt);
  }
  for(int i = 0 ; i < mChunks ; ++i) {
    //float teta = mRandom->Gaus(jet.Eta(), jetspread);
    //float tphi = mRandom->Gaus(jet.Phi(), jetspread);
    float R = mRandom->Exp(1/(mJetSpreadA +mJetSpreadB * jet.E()));
    float PHI = mRandom->Uniform(2 * M_PI);
    //std::cout << "E:" << jet.E() << "  R:" << R << '\n';
    float teta = jet.Eta() + R * cos(PHI);
    float tphi = jet.Phi() + R * sin(PHI);
    
    tphi = TVector2::Phi_0_2pi(tphi);
    if(std::abs(teta) > 3.33333) {
      //std::cout << "chunk outside simulated calo\n";
      if(mNoOutOfCone) --i;
      else lostPt += dpt;
      continue;
    }
    double de = dpt/cos(tphi-jphi);
    int ie, ip;
    calIds(teta, tphi, ie, ip); 
    //std::cout << "vorher:" << teta << ", " << tphi << ", " << ie << ", " 
    //	      << ip << "\n";
    float dphi = TVector2::Phi_mpi_pi(tphi-jphi);
    float deta = teta-jet.Eta();
    if(sqrt(deta*deta + dphi*dphi) > 0.5) {
      //std::cout << "Out of cone:" << teta << ":" << jet.Eta() << " , " << dphi << '\n';
      if(mNoOutOfCone) --i;
      else lostPt += dpt;
      continue;
    }
    int id = 0;
    assert(ie*1000 + ip != 0);
    TowerMap::const_iterator towit = towers.find(ie*1000 + ip);
    if(towit != towers.end()) {
      id = towit->second;
    } else {
      ++ntowers;
      id = ntowers;
      towers[ie*1000 + ip] = id;
      et[id-1] = 0;
    }
    --id;
    tow.SetPtEtaPhiM(de,teta,tphi,0);
    //tow *= de/tow.E();
    et[id] += tow.Pt();
    eta[id] = teta;
    phi[id] = TVector2::Phi_mpi_pi(tphi);
    ieta[id] = ie;
    iphi[id] = ip;
    rec += tow;
  }
  //std::cout  << "Eta:" << jet.Eta() <<  "       : " << rec.Pt() << "," << rec.E() << "  == " << jet.Pt() << "," << jet.E() << '\n';
  //std::cout << "lost energy:" << lostE/jet.E() << '\n';
  if(mNoOutOfCone) {
    double scale = jet.Pt()/rec.Pt();
    assert(scale < 1.1); 
    TLorentzVector rec2(0,0,0,0);
    for(int i = 0 ; i < ntowers ; ++i) {
      et[i] *= scale; 
      tow.SetPtEtaPhiM(et[i],eta[i],phi[i],0);
      rec2 += tow;
    }
    assert(std::abs((rec2.Pt()-jet.Pt())/jet.Pt()) < 0.001); 
    //std::cout << " vorher:" << scale << "  nachher:" << rec2.E()/jet.E() << "\n";
  }
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
  float towemtrue[kMaxTower];
  float towhdtrue[kMaxTower];
  float towoetrue[kMaxTower];
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
  CalibTree->Branch("TowEmTrue",towemtrue,"TowEmTrue[NobjTowCal]/F");
  CalibTree->Branch("TowHadTrue",towhdtrue,"TowHadTrue[NobjTowCal]/F");
  CalibTree->Branch("TowOETrue",towoetrue,"TowOETrue[NobjTowCal]/F");
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
    smearTower(mPinput.E(),towen[0],towem[0],towhd[0],towoe[0],towemtrue[0],towhdtrue[0],towoetrue[0]);
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
  float towemtrue[kMaxTower];
  float towhdtrue[kMaxTower];
  float towoetrue[kMaxTower];
  int towid_phi[kMaxTower];
  int towid_eta[kMaxTower];
  int towid [kMaxTower];
  float jcalpt,jcalphi,jcaleta,jcalet,jcale;
  float jgenpt,jgenphi,jgeneta,jgenet,jgene;
  float mcalmet,mcalphi,mcalsum;
  float photonpt,photonphi,photoneta,photonet,photone;
  float weight = 1.0;
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
  CalibTree->Branch("TowEmTrue",towemtrue,"TowEmTrue[NobjTowCal]/F");
  CalibTree->Branch("TowHadTrue",towhdtrue,"TowHadTrue[NobjTowCal]/F");
  CalibTree->Branch("TowOETrue",towoetrue,"TowOETrue[NobjTowCal]/F");
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
  CalibTree->Branch("Weight",&weight,"Weight/F");
  

  TLorentzVector jet,genjet, tower;
  for(int i = 0; i < nevents ; ++i) {
    genInput();
    photonpt = mPinput.Pt();
    photoneta = mPinput.Eta();
    photonphi = mPinput.Phi();
    photonet = mPinput.Pt();
    photone = mPinput.E();
    //jgenpt = mRandom->Gaus(photonpt,0.04 * photonpt);
    jgenpt = photonet;
    jgeneta = mRandom->Gaus(photoneta,1.0);
    if((jgeneta > 3.0) || (jgeneta < -3.0)) {
      --i;
      continue;
    }
    jgenphi = photonphi +M_PI;
    genjet.SetPtEtaPhiM(jgenpt,jgeneta,jgenphi,0);
    NobjTowCal = splitJet(genjet,towet,toweta,towphi,towid_eta,towid_phi);
    jet.SetPtEtaPhiM(0,0,0,0);
    genjet.SetPtEtaPhiM(0,0,0,0);
    double towmean = mTowConst[0];
    if(mModel == Flat) mTowConst[0] = 1/mRandom->Uniform(1.5);
    else if(mModel == Exp) mTowConst[0] = 1/mRandom->Exp(0.5);
    else if(mModel == Slope) {
      double u1 = mRandom->Uniform(2);
      double u2 = mRandom->Uniform(2);
      mTowConst[0] = 1/(2 - std::max(u1,u2));
    }
    double p0frac = mRandom->Uniform(mMaxPi0Frac);
    for(int j = 0; j < NobjTowCal ; ++j) {
      tower.SetPtEtaPhiM(towet[j],toweta[j],towphi[j],0);
      towen[j] =  tower.E();
      genjet += tower;
      smearTower((1 - p0frac) * tower.E(),towen[j],towem[j],towhd[j],towoe[j],
		 towemtrue[j],towhdtrue[j],towoetrue[j]); 
      towen[j] += p0frac * tower.E();
      towem[j] += p0frac * tower.E();
      towemtrue[j] += p0frac * tower.E();
      tower *= towen[j]/tower.E();
      towet[j] = tower.Pt();
      jet += tower;
    }
    mTowConst[0] = towmean;
    jcalpt = jet.Pt();
    jcaleta = jet.Eta();
    jcalphi = jet.Phi();
    jcalet = jet.Pt();
    jcale = jet.E(); 
    jgenphi = genjet.Phi();
    jgenet = genjet.Pt();
    jgenpt = genjet.Pt();
    jgene = genjet.E();
    jgeneta = genjet.Eta();
    mcalmet = 0;
    mcalphi = 0;
    mcalsum = 0;
    CalibTree->Fill();
    if(i % 1000 == 0) std::cout << "generating event " << i << '\n';
  } 
  return CalibTree->GetEntriesFast();
}

int ToyMC::generateDiJetTree(TTree* CalibTree, int nevents)
{
  //make tree 
  const int kMAX = 1000;
  int NobjTow;
  float towet[kMAX];
  float toweta[kMAX];
  float towphi[kMAX];
  float towen[kMAX];
  float towem[kMAX];
  float towhd[kMAX];
  float towoe[kMAX];  
  float towemtrue[kMAX];
  float towhdtrue[kMAX];
  float towoetrue[kMAX];
  int towid_phi[kMAX];
  int towid_eta[kMAX];
  int towid[kMAX];
  int tow_jetidx[kMAX];
  
  float ttowet[kMAX];
  float ttoweta[kMAX];
  float ttowphi[kMAX];
  int ttowid_phi[kMAX];
  int ttowid_eta[kMAX]; 

  const int kjMAX = 2;
  int NobjJet = 2;
  float jetpt[kjMAX];
  float jetphi[kjMAX];
  float jeteta[kjMAX];
  float jetet[kjMAX];
  float jete[kjMAX]; 
  float jetgenpt[kjMAX];
  float jetgenphi[kjMAX];
  float jetgeneta[kjMAX];
  float jetgenet[kjMAX];
  float jetgene[kjMAX];
  float weight = 1;
  // CaloTower branches
  CalibTree->Branch("NobjTow",&NobjTow,"NobjTow/I");
  CalibTree->Branch("TowId",towid,"TowId[NobjTow]/I");
  CalibTree->Branch("TowId_phi",towid_phi,"TowId_phi[NobjTow]/I");
  CalibTree->Branch("TowId_eta",towid_eta,"TowId_eta[NobjTow]/I");
  CalibTree->Branch("TowEt",towet,"TowEt[NobjTow]/F");
  CalibTree->Branch("TowEta",toweta,"TowEta[NobjTow]/F");
  CalibTree->Branch("TowPhi",towphi,"TowPhi[NobjTow]/F");
  CalibTree->Branch("TowE",towen,"TowE[NobjTow]/F");
  CalibTree->Branch("TowEm",towem,"TowEm[NobjTow]/F");
  CalibTree->Branch("TowHad",towhd,"TowHad[NobjTow]/F");
  CalibTree->Branch("TowOE",towoe,"TowOE[NobjTow]/F"); 
  CalibTree->Branch("TowEmTrue",towemtrue,"TowEmTrue[NobjTowCal]/F");
  CalibTree->Branch("TowHadTrue",towhdtrue,"TowHadTrue[NobjTowCal]/F");
  CalibTree->Branch("TowOETrue",towoetrue,"TowOETrue[NobjTowCal]/F");
  CalibTree->Branch("Tow_jetidx",tow_jetidx,"Tow_jetidx[NobjTow]/I");
  // Jet-specific branches of the tree
  CalibTree->Branch("NobjJet",&NobjJet,"NobjJet/I"             );
  CalibTree->Branch("JetPt",jetpt,"JetPt[NobjJet]/F" );
  CalibTree->Branch("JetPhi",jetphi,"JetPhi[NobjJet]/F");
  CalibTree->Branch("JetEta",jeteta,"JetEta[NobjJet]/F");
  CalibTree->Branch("JetEt",jetet,"JetEt[NobjJet]/F" );
  CalibTree->Branch("JetE",jete,"JetE[NobjJet]/F"  );
  CalibTree->Branch("Weight",&weight,"Weight/F"  );


  TLorentzVector jet[2],genjet[2], tower;
  for(int n = 0; n < nevents ; ++n) {
    genInput();
    jetgenpt[0] = mPinput.Pt();
    jetgeneta[0] = mPinput.Eta();
    jetgenphi[0] = mPinput.Phi();
    jetgenet[0] = mPinput.Pt();
    jetgene[0] = mPinput.E();
    //jgenpt = mRandom->Gaus(photonpt,0.04 * photonpt);
    jetgenpt[1] = jetgenpt[0];
    jetgeneta[1] = mRandom->Gaus(jetgeneta[0],1.0);
    if((jetgeneta[1] > 3.0) || (jetgeneta[1] < -3.0)) {
      --n;
      continue;
    }
    jetgenphi[1] = jetgenphi[0] +M_PI;
    if(mModel == Flat) mTowConst[0] = 1/mRandom->Uniform(1.5);
    else if(mModel == Exp) mTowConst[0] = 1/mRandom->Exp(0.5);
    else if(mModel == Slope) {
      double u1 = mRandom->Uniform(2);
      double u2 = mRandom->Uniform(2);
      mTowConst[0] = 1/(2 - std::max(u1,u2));
    }
    NobjTow = 0;
    for(int i = 0 ; i < NobjJet ; ++i) {
      genjet[i].SetPtEtaPhiM(jetgenpt[i],jetgeneta[i],jetgenphi[i],0);
      int ntow = splitJet(genjet[i],ttowet,ttoweta,ttowphi,ttowid_eta,
			  ttowid_phi);  
      jet[i].SetPtEtaPhiM(0,0,0,0);
      genjet[i].SetPtEtaPhiM(0,0,0,0);
      double p0frac = mRandom->Uniform(mMaxPi0Frac);
      for(int j = 0 ; j < ntow ; ++j) {
	int k = NobjTow + j;
	towet[k] = ttowet[j];
	toweta[k] = ttoweta[j];
	towphi[k] = ttowphi[j];	
	towid[k] = 0;
	towid_eta[k] = ttowid_eta[j];
	towid_phi[k] = ttowid_phi[j];
	tow_jetidx[k] = i;
	tower.SetPtEtaPhiM(towet[k],toweta[k],towphi[k],0);
	towen[k] =  tower.E();
	genjet[i] += tower;
	smearTower((1 - p0frac) * tower.E(),towen[k],towem[k],towhd[k],
		   towoe[k],towemtrue[k],towhdtrue[k],towoetrue[k]); 
	towen[k] += p0frac * tower.E();
	towem[k] += p0frac * tower.E();
	towemtrue[k] += p0frac * tower.E();
	tower *= towen[k]/tower.E();
	towet[k] = tower.Pt();
	jet[i] += tower;
      }
      NobjTow += ntow; 
      jetpt[i] = jet[i].Pt();
      jetphi[i] = jet[i].Phi();
      jeteta[i] = jet[i].Eta();
      jetet[i] = jet[i].Pt();
      jete[i] = jet[i].E(); 
      jetgenpt[i] = genjet[i].Pt();
      jetgenphi[i] = genjet[i].Phi();
      jetgeneta[i] = genjet[i].Eta();
      jetgenet[i] = genjet[i].Pt();
      jetgene[i] = genjet[i].E();
    }
    if((jeteta[0] < mMinEta) || (jeteta[0] > mMaxEta)) {
      --n;
      continue;
    }
    if((jeteta[1] > mMinEta) && (jeteta[1] < mMaxEta)
       && (jetpt[1] > jetpt[0])) {
      //swap jets
      jetpt[0] = jet[1].Pt();
      jetphi[0] = jet[1].Phi();
      jeteta[0] = jet[1].Eta();
      jetet[0] = jet[1].Pt();
      jete[0] = jet[1].E(); 
      jetgenpt[0] = genjet[1].Pt();
      jetgenphi[0] = genjet[1].Phi();
      jetgeneta[0] = genjet[1].Eta();
      jetgenet[0] = genjet[1].Pt();
      jetgene[0] = genjet[1].E();
      jetpt[1] = jet[0].Pt();
      jetphi[1] = jet[0].Phi();
      jeteta[1] = jet[0].Eta();
      jetet[1] = jet[0].Pt();
      jete[1] = jet[0].E(); 
      jetgenpt[1] = genjet[0].Pt();
      jetgenphi[1] = genjet[0].Phi();
      jetgeneta[1] = genjet[0].Eta();
      jetgenet[1] = genjet[0].Pt();
      jetgene[1] = genjet[0].E();
      for(int j = 0 ; j < NobjTow ; ++j) {
	tow_jetidx[j] = (tow_jetidx[j] == 0) ? 1 : 0;
      }
    }
    CalibTree->Fill();
    if(n % 1000 == 0) std::cout << "generating event " << n << '\n';
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

int ToyMC::makeDiJet(const char* filename, int nevents) {
  TFile* file = new TFile(filename,"recreate");
  TTree* CalibTree = new TTree("DiJetTree","DiJetTree");

  nevents = generateDiJetTree(CalibTree,nevents);
  //file->ls();
  file->Write();
  file->Close();
  return nevents;
}

void ToyMC::init(const std::string& configfile) {
  ConfigFile config(configfile.c_str());

   mMinEta = config.read<double>("ToyMC min eta",-2.5);
   mMaxEta = config.read<double>("ToyMC max eta",2.5);
   mMinPt  = config.read<double>("ToyMC min pt",30);
   mMaxPt  = config.read<double>("ToyMC max pt",400);
   std::string spectrum = config.read<std::string>("ToyMC pt spectrum","uniform");
   if(spectrum == "powerlaw") {
     mPtSpectrum = PowerLaw; 
   } else if(spectrum == "uniform") {
     mPtSpectrum = Uniform;
   } else {
     std::cerr << "unknown ToyMC pt spectrum:" << spectrum << '\n';
     exit(1);
   }
   std::vector<double> auter = bag_of<double>(config.read<string>("ToyMC tower const","1.25 0.0 1.0 1.0 0.0"));
   while(auter.size() < 5) auter.push_back(0);
   assert(auter.size() == 5);
   for(unsigned int i = 0; i < auter.size() ; ++i) mTowConst[i] = auter[i];
   mResoStochastic = config.read<double>("ToyMC tower resolution stochastic",1.20);
   mResoNoise = config.read<double>("ToyMC tower resolution noise",0.05);
   mJetSpreadA = config.read<double>("ToyMC jet spread A",0.5);
   mJetSpreadB = config.read<double>("ToyMC jet spread B",0);
   mNoOutOfCone = config.read<bool>("ToyMC avoid out-of-cone",true);
   std::string model = config.read<std::string>("ToyMC model","gauss");
   if(model == "gauss") {
     mModel = Gauss;
   } else if(model  == "landau") {
     mModel = Landau;
   } else if(model == "flat") {
     mModel = Flat;
     mTowConst[1] = 0;
   } else if(model == "exp") {
     mModel = Exp;
     mTowConst[1] = 0;
   } else if(model == "slope") {
     mModel = Slope;
     mTowConst[1] = 0;
   } else {
     std::cerr << "unknown ToyMC model:" << model << '\n';
     exit(1);
   }
   mChunks = config.read<int>("ToyMC chunks",200);
   mMaxPi0Frac = config.read<double>("ToyMC max pi0 fraction",0.5);
   mMaxEmf = config.read<double>("ToyMC tower max EMF",0.5);
   int seed = config.read<int>("ToyMC seed",0); 
   mRandom->SetSeed(seed);
}

void ToyMC::print() const {
  std::cout << "ToyMC configuration:\n";
  std::cout << "  primary: " << mMinEta << " < eta < " << mMaxEta << '\n';
  std::cout << "           " << mMinPt << " < pt < " << mMaxPt 
	    << " spectrum = " << mPtSpectrum << '\n';
  std::cout << "    tower: c = (" << mTowConst[0] <<", " << mTowConst[1]
	    << ", " << mTowConst[2] << ", " << mTowConst[3] << ", "
	    << mTowConst[4] << ") stoch = " << mResoStochastic << " noise = "
	    << mResoNoise << " max EMF = " << mMaxEmf << " max piO = "
	    << mMaxPi0Frac << '\n';
  std::cout << "     jets: spread A = " << mJetSpreadA << " B = " 
	    << mJetSpreadB << " avoid out-of-cone = " << mNoOutOfCone
	    << '\n';
  std::cout << "    model: n chunks = " <<  mChunks << " model = "
	    << mModel << '\n';
}
