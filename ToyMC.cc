// $Id: ToyMC.cc,v 1.25 2009/06/21 18:16:00 mschrode Exp $

#include "ToyMC.h"

#include <cmath>
#include <iostream>
#include <map>
#include <cassert> 
#include <ext/hash_map>

#include "TF1.h"
#include "TFile.h"
#include "TRandom3.h"
#include "TTree.h"

#include "ConfigFile.h"


//!  \brief Default constructor
// -----------------------------------------------------------------
ToyMC::ToyMC() : mType(1), mMinEta(-2.5),mMaxEta(2.5),mMinPt(30), mMaxPt(400),mPtSpectrum(Uniform),
		 mChunks(200),mJetSpreadA(0.5),mJetSpreadB(0),mNoOutOfCone(true),mMaxPi0Frac(0.5),
		 mMaxEmf(0.5),mResponseModel(Constant),mHistResp(0),
		 mResolutionModel(Gauss), mSmearFactor(1.)
{
  mRandom      = new TRandom3();
  mRandom->SetSeed(0);
  mHistResp = 0;
}


//!  \brief Generates the input (truth) of an event
//!
//!  Generates a random input lorentz vector of an event 
//!  (e.g. the lorentz vector of the photon for photon-jet
//!  events) which is then smeared due to detector response
//!  and  resolution effects etc.
//!  The randomly generated components are
//!  - pt between mMinPt and mMaxPt according to mPtSpectrum
//!  - eta uniformly between mMinEta and mMaxEta
//!  - phi uniformly between 0 and 2Pi
//!  - a zero mass
// -----------------------------------------------------------------
void ToyMC::genInput() { 
  static double rand[3];
  mRandom->RndmArray(3,rand);
  double pt = 0;  
  if(mPtSpectrum == Uniform) {
    mRandom->RndmArray(3,rand);
    pt = rand[0]*(mMaxPt - mMinPt)+mMinPt;
  } else if(mPtSpectrum == PowerLaw) {
    pt = mMinPt * pow(rand[0],-1.0/5.5);
  }
  mPinput.SetPtEtaPhiM(pt,
		       rand[1]*(mMaxEta - mMinEta)+mMinEta,
		       rand[2]*2 * M_PI - M_PI, 0);
}


//!  \brief Transform eta and phi values of the center of
//!         the corresponding tower and find the tower indices
//!
//!  \param eta Eta, is transformed to eta of corresponding tower center
//!  \param phi Phi, is transformed to phi of corresponding tower center
//!  \param ieta Index of tower covering eta
//!  \param iphi Index of tower covering phi
// -----------------------------------------------------------------
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



//!  \brief Calculate pt smear factor due to
//!         response and resolution
//!  \param pt Pt for calculation of some smear factors
//----------------------------------------------------------
void ToyMC::CalculateSmearFactor(double pt) {
  // Reset smear factor
  mSmearFactor = 1.;

  // Apply resolution
  if( mResponseModel == Constant
      || (mResponseModel == Flat)
      || (mResponseModel == Exp)
      || (mResponseModel == Slope)) {
    mSmearFactor *= mParResp.at(0);
  }
  else if( mResponseModel == L3 ) { 
    if( mPinput.Pt() < 0.1 ) {
      mSmearFactor *= mParResp.at(0) - mParResp.at(1)/mParResp.at(3) +  mParResp.at(4);
    } else {
      mSmearFactor *= mParResp.at(0) - mParResp.at(1)/(pow(log10(pt),mParResp.at(2)) + mParResp.at(3)) + mParResp.at(4)/pt;
    }
  }
  else if( mResponseModel == SimpleInverse ) { 
    mSmearFactor *= 1. - mParResp.at(0)/(pt + mParResp.at(1));
  }

  // Apply resolution
  double smear = 1.;
  if( mResolutionModel == Landau) {
    do {
      smear =  mRandom->Landau(1,sqrt(mParReso.at(0)*mParReso.at(0)/pt + mParReso.at(1)*mParReso.at(1)));
    } while((smear < 0) || (smear > 2));
    smear = 2 - smear;
  }
  else if ( (mResolutionModel == Gauss) ) {
    do {
      smear = mRandom->Gaus(1.0,sqrt(mParReso.at(0)*mParReso.at(0)/pt + mParReso.at(1)*mParReso.at(1)));
    } while((smear < 0) || (smear > 2));
  }
  else if( mResolutionModel == GaussUniform ) {
    do{
      smear = mRandom->Gaus( 1.0, mParReso.at(0) );
    } while( smear < 0 || smear > 2 );
    if( mRandom->Uniform() < mParReso.at(1) )
      smear = mRandom->Uniform();
  }
  else if( mResolutionModel == TwoGauss ) {
    do {
      smear = mHistResp->GetRandom();
    } while ( smear < 0.3 );
  }

  mSmearFactor *= smear;
}




//!  \brief Calculate emf, scale hadronic tower energy with response factor,
//!         and smear with resolution
//!
//!  \param e True tower energy without pi0 part
//!  \param calcSmearFactor If true, the smear factor is calculated newly during this function call
//!  \param te Tower energy after scaling
//!  \param tem Measured em part of tower energy after scaling and smearing
//!  \param thad Measured had part of tower energy after scaling and smearing
//!  \param tout Measured HO part of tower energy after scaling and smearing
//!  \param temtrue True em part of tower energy after scaling
//!  \param thadtrue True had part of tower energy after scaling
//!  \param touttrue True HO part of tower energy after scaling
// -----------------------------------------------------------------
void ToyMC::smearTower(double e, bool calcSmearFactor, float& te, float& tem, float& thad, float& tout, 
		       float& temtrue, float& thadtrue, float& touttrue) 
{
  // Generate emf and set electromagnetic
  // and hadronic fraction; smear hadronic
  // fraction with response
  touttrue   = 0.;
  tout       = touttrue;
  float emf  = mRandom->Uniform(mMaxEmf);
  temtrue    = emf * e;
  tem        = temtrue;
  thadtrue   = (1-emf) * e;
  thad       = thadtrue;

  // Apply response and resolution to hadronic fraction
  if( calcSmearFactor ) {
    if( mSmearTowersIndividually ) CalculateSmearFactor(thad);
    else                           CalculateSmearFactor(mPinput.Pt());
  }
  thad      *= mSmearFactor;

  // Add up tower parts to total tower energy
  te = tem + thad + tout;
}



//!  \brief Split generated jet into towers
//!
//!  Splits generated true jet pt into 'mChunks' portions
//!  (particles), spreads them in eta and phi and sums up
//!  tower pt.
// -----------------------------------------------------------------
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



// -----------------------------------------------------------------
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

    // Calculate smear factor
    bool calcSmearFactor = false;
    if( i == 0 || mSmearTowersIndividually ) calcSmearFactor = true;

    smearTower(mPinput.E(),calcSmearFactor,towen[0],towem[0],towhd[0],towoe[0],towemtrue[0],towhdtrue[0],towoetrue[0]);
    towet[0] = towen[0]/mPinput.E() * mPinput.Pt();
    CalibTree->Fill();
    if(i % 1000 == 0) std::cout << "generated event " << i << '\n';
  }
  return CalibTree->GetEntriesFast(); 
}



// -----------------------------------------------------------------
int ToyMC::makeTrackCluster(const char* filename, int nevents) {
  TFile* file = new TFile(filename,"recreate");
  TTree* CalibTree = new TTree("TrackClusterTree","TrackClusterTre");
 
  nevents = generateTrackClusterTree(CalibTree, nevents);
  file->Write();
  file->Close();
  return  nevents;
}



//!  \brief Generate photon-jet events and write into tree
//!  \param CalibTree ROOT tree
//!  \param nevents Number of photon-jet events
//!  \return Number of generated events
// -----------------------------------------------------------------
int ToyMC::generatePhotonJetTree(TTree* CalibTree, int nevents)
{
  //make tree 
  const int kMaxTower = 1000;
  int NobjTowCal,NobjTrack = 0;
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
  const int kMaxTracks = 1;
  float trackpt[kMaxTracks];
  float tracketa[kMaxTracks];
  float trackphi[kMaxTracks];
  float trackp[kMaxTracks];
  float trackdr[kMaxTracks];
  float tracketaout[kMaxTracks];
  float trackphiout[kMaxTracks];
  float trackdrout[kMaxTracks];
  float trackemc1[kMaxTracks];
  float trackemc3[kMaxTracks];
  float trackemc5[kMaxTracks];
  float trackhac1[kMaxTracks];
  float trackhac3[kMaxTracks];
  float trackhac5[kMaxTracks];
  int tracktowid[kMaxTracks];
  int tracktowidphi[kMaxTracks];
  int tracktowideta[kMaxTracks];
  int trackid[kMaxTracks];
  int tracknhits[kMaxTracks];
  int trackQualityL[kMaxTracks];
  int trackQualityT[kMaxTracks];
  int trackQualityHP[kMaxTracks];
  float trackchi2[kMaxTracks];
  float muDR[kMaxTracks];
  float muDE[kMaxTracks];

  float jcalpt,jcalphi,jcaleta,jcalet,jcale;

  // All correction factors are 1 in ToyMC
  float jscaleZSP    = 1.;
  float jscalel2     = 1.;
  float jscalel3     = 1.;
  float jscalel23    = 1.;
  float jscaleJPT    = 1.;
  float jscalel23JPT = 1.;

  float jgenpt,jgenphi,jgeneta,jgenet,jgene;
  float mcalmet,mcalphi,mcalsum;
  float photonpt,photonphi,photoneta,photonet,photone; 
  float gphotonpt,gphotonphi,gphotoneta,gphotonet,gphotone;
  float nonleadingjetspt = 0.0;
  float weight = 1.0;
  CalibTree->Branch("NobjTowCal",&NobjTowCal,"NobjTowCal/I");   
  CalibTree->Branch("NobjTrack", &NobjTrack, "NobjTowCal/I");
  //tower branches
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
  //track branches
  CalibTree->Branch( "NobjTrack",  &NobjTrack, "NobjTrack/I"             );
  CalibTree->Branch( "TrackTowId", tracktowid, "TrackTowId[NobjTrack]/I" );
  CalibTree->Branch( "TrackTowIdPhi", tracktowidphi, "TrackTowIdPhi[NobjTrack]/I" );
  CalibTree->Branch( "TrackTowIdEta", tracktowideta, "TrackTowIdEta[NobjTrack]/I" );
  CalibTree->Branch( "TrackId",    trackid,    "TrackId[NobjTrack]/I"    );
  CalibTree->Branch( "TrackNHits", tracknhits, "TrackNHits[NobjTrack]/I" );
  CalibTree->Branch( "TrackQualityL",trackQualityL,"TrackQualityL[NobjTrack]/O");
  CalibTree->Branch( "TrackQualityT",trackQualityT,"TrackQualityT[NobjTrack]/O");
  CalibTree->Branch( "TrackQualityHP",trackQualityHP,"TrackQualityHP[NobjTrack]/O");
  CalibTree->Branch( "TrackChi2",  trackchi2,  "TrackChi2[NobjTrack]/F"  );
  CalibTree->Branch( "TrackPt",    trackpt,    "TrackPt[NobjTrack]/F"    );
  CalibTree->Branch( "TrackEta",   tracketa,   "TrackEta[NobjTrack]/F"   );
  CalibTree->Branch( "TrackPhi",   trackphi,   "TrackPhi[NobjTrack]/F"   );
  CalibTree->Branch( "TrackP" ,    trackp,     "TrackP[NobjTrack]/F"     );
  CalibTree->Branch( "TrackDR" ,   trackdr,    "TrackDR[NobjTrack]/F"    );
  CalibTree->Branch( "TrackPhiOut",trackphiout,"TrackPhiout[NobjTrack]/F");
  CalibTree->Branch( "TrackEtaOut",tracketaout,"TrackEtaout[NobjTrack]/F");
  CalibTree->Branch( "TrackDROut", trackdrout, "TrackDRout[NobjTrack]/F" );
  CalibTree->Branch( "TrackEMC1",  trackemc1,  "TrackEMC1[NobjTrack]/F"  );
  CalibTree->Branch( "TrackEMC3",  trackemc3,  "TrackEMC3[NobjTrack]/F"  );
  CalibTree->Branch( "TrackEMC5",  trackemc5,  "TrackEMC5[NobjTrack]/F"  );
  CalibTree->Branch( "TrackHAC1",  trackhac1,  "TrackHAC1[NobjTrack]/F"  );
  CalibTree->Branch( "TrackHAC3",  trackhac3,  "TrackHAC3[NobjTrack]/F"  );
  CalibTree->Branch( "TrackHAC5",  trackhac5,  "TrackHAC5[NobjTrack]/F"  );
  CalibTree->Branch( "MuDR", muDR,  "MuDR[NobjTrack]/F"  );
  CalibTree->Branch( "MuDE", muDE,  "MuDE[NobjTrack]/F"  );

  // Jet- MEt-specific branches of the tree
  CalibTree->Branch("JetCalPt",&jcalpt,"JetCalPt/F");
  CalibTree->Branch("JetCalPhi",&jcalphi,"JetCalPhi/F");
  CalibTree->Branch("JetCalEta",&jcaleta,"JetCalEta/F");
  CalibTree->Branch("JetCalEt",&jcalet,"JetCalEt/F");
  CalibTree->Branch("JetCalE",&jcale,"JetCalE/F");  
  CalibTree->Branch( "JetCorrZSP",&jscaleZSP, "JetCorrZSP/F" );
  CalibTree->Branch( "JetCorrL2", &jscalel2,  "JetCorrL2/F" );
  CalibTree->Branch( "JetCorrL3", &jscalel3,  "JetCorrL3/F" );
  CalibTree->Branch( "JetCorrL2L3", &jscalel23,  "JetCorrL2L3/F" );
  CalibTree->Branch( "JetCorrJPT",&jscaleJPT, "JetCorrJPT/F" );
  CalibTree->Branch( "JetCorrL2L3JPT", &jscalel23JPT,  "JetCorrL2L3JPT/F" );

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
  // GenPhotons branches
  CalibTree->Branch( "GenPhotonPt",  &gphotonpt,  "GenPhotonPt/F"  );
  CalibTree->Branch( "GenPhotonPhi", &gphotonphi, "GenPhotonPhi/F" );
  CalibTree->Branch( "GenPhotonEta", &gphotoneta, "GenPhotonEta/F" );
  CalibTree->Branch( "GenPhotonEt",  &gphotonet,  "GenPhotonEt/F"   );
  CalibTree->Branch( "GenPhotonE",   &gphotone,   "GenPhotonE/F"   );
  // NonLeadingJetPt branch
  CalibTree->Branch( "NonLeadingJetPt", &nonleadingjetspt,   "NonLeadingJetPt/F"   );


  CalibTree->Branch("EventWeight",&weight,"EventWeight/F");
  

  // Generate events
  TLorentzVector jet;
  TLorentzVector genjet;
  TLorentzVector tower;

  // Loop over events
  for(int i = 0; i < nevents ; ++i) {
    // Generate truth 4-momentum
    genInput();

    // Assign it's variables to genphoton
    // and genphoton (perfectly measured)
    photonpt   = mPinput.Pt();
    photoneta  = mPinput.Eta();
    photonphi  = mPinput.Phi();
    photonet   = mPinput.Pt();
    photone    = mPinput.E();
    gphotonpt  = photonpt;
    gphotoneta = photoneta;
    gphotonphi = photonphi;
    gphotonet  = photonet;
    gphotone   = photone;

    // Gen jet
    jgenpt     = gphotonpt;
    jgeneta    = mRandom->Gaus(photoneta,1.0);
    if((jgeneta > 3.0) || (jgeneta < -3.0)) {
      --i;
      continue;
    }
    jgenphi    = photonphi + M_PI;

    // Create 4-momentum of genjet (massless)
    genjet.SetPtEtaPhiM(jgenpt,jgeneta,jgenphi,0);

    // Split it into towers and set tower truth
    NobjTowCal = splitJet(genjet,towet,toweta,towphi,towid_eta,towid_phi);

    // Reset jet and genjet 4-momenta. They will
    // be repopulated by sum of tower 4-momenta
    jet.SetPtEtaPhiM(0,0,0,0);
    genjet.SetPtEtaPhiM(0,0,0,0);

    // Set random response parameters for this event
    // if required by response model
    if      (mResponseModel == Flat) mParResp.at(0) = mRandom->Uniform(1.5);
    else if (mResponseModel == Exp)  mParResp.at(0) = mRandom->Exp(0.5);
    else if (mResponseModel == Slope) {
      double u1 = mRandom->Uniform(2);
      double u2 = mRandom->Uniform(2);
      mParResp.at(0) = 2. - std::max(u1,u2);
    }

    // Generate pi0 fraction i.e. non-hadronic fraction
    double p0frac = mRandom->Uniform(mMaxPi0Frac);

    // Loop over towers and smear truth with response factor
    for(int j = 0; j < NobjTowCal ; ++j) {
      // Set tower truth 4-momentum
      tower.SetPtEtaPhiM(towet[j],toweta[j],towphi[j],0);
      towen[j]  =  tower.E();

      // Add tower truth to genjet
      genjet   += tower;

      // Calculate smear factor
      bool calcSmearFactor = false;
      if( j == 0 || mSmearTowersIndividually ) calcSmearFactor = true;

      // Smear hadronic par of tower energy
      smearTower((1 - p0frac) * tower.E(),calcSmearFactor,towen[j],towem[j],towhd[j],towoe[j],
		 towemtrue[j],towhdtrue[j],towoetrue[j]); 

      // Add remaining em part
      towen[j]     += p0frac * tower.E();
      towem[j]     += p0frac * tower.E();
      towemtrue[j] += p0frac * tower.E();

      // Multiply tower 4-momentum with response
      tower        *= towen[j]/tower.E();
      towet[j]      = tower.Pt();

      // Add tower to jet
      jet += tower;
    } // End loop over towers

    // Set jet and genjet variables
    jcalpt   = jet.Pt();
    jcaleta  = jet.Eta();
    jcalphi  = jet.Phi();
    jcalet   = jet.Pt();
    jcale    = jet.E(); 
    jgenphi  = genjet.Phi();
    jgenet   = genjet.Pt();
    jgenpt   = genjet.Pt();
    jgene    = genjet.E();
    jgeneta  = genjet.Eta();
    mcalmet  = 0;
    mcalphi  = 0;
    mcalsum  = 0;

    // Fill tree
    CalibTree->Fill();
    if(i % 1000 == 0) std::cout << "generating event " << i << '\n';
  } 
  return CalibTree->GetEntriesFast();
}



//!  \brief Generate dijet events and write into tree
//!  \param CalibTree ROOT tree
//!  \param nevents Number of dijet events
//!  \return Number of generated events
//----------------------------------------------------------
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

  // All correction factors are 1 in ToyMC
  float jscaleZSP[2] = { 1., 1. };
  float jscalel2[2] = { 1., 1. };
  float jscalel3[2] = { 1., 1. };
  float jscalel23[2] = { 1., 1. };
  float jscaleJPT[2] = { 1., 1. };
  float jscalel23JPT[2] = { 1., 1. };

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
  CalibTree->Branch("GenJetEt",jetgenet,"GenJetEt[NobjJet]/F" );
  CalibTree->Branch("GenJetPt",jetgenpt,"GenJetPt[NobjJet]/F" );

  CalibTree->Branch( "JetCorrZSP",     jscaleZSP, "JetCorrZSP[NobjJet]/F" );
  CalibTree->Branch( "JetCorrL2",      jscalel2,  "JetCorrL2[NobjJet]/F" );
  CalibTree->Branch( "JetCorrL3",      jscalel3,  "JetCorrL3[NobjJet]/F" );
  CalibTree->Branch( "JetCorrL2L3",    jscalel23,  "JetCorrL2L3[NobjJet]/F" );
  CalibTree->Branch( "JetCorrJPT",     jscaleJPT, "JetCorrJPT[NobjJet]/F" );
  CalibTree->Branch( "JetCorrL2L3JPT", jscalel23JPT,  "JetCorrL2L3JPT[NobjJet]/F" );


  // Generate events
  TLorentzVector jet[2];
  TLorentzVector genjet[2];
  TLorentzVector tower;

  // Loop over events
  for(int n = 0; n < nevents ; ++n) {
    // Generate truth 4-momentum
    genInput();

    // Assign it's variables to first genjet
    jetgenpt[0]  = mPinput.Pt();
    jetgeneta[0] = mPinput.Eta();
    jetgenphi[0] = mPinput.Phi();
    jetgenet[0]  = mPinput.Pt();
    jetgene[0]   = mPinput.E();

    // Second genjet gets random eta
    // between -3 and 3, and phi+PI
    jetgenpt[1]  = jetgenpt[0];
    jetgeneta[1] = mRandom->Gaus(jetgeneta[0],1.0);
    if((jetgeneta[1] > 3.0) || (jetgeneta[1] < -3.0)) {
      --n;
      continue;
    }
    jetgenphi[1] = jetgenphi[0] + M_PI;

    // Set random response paramters for this event
    // if required by response model
    if     (mResponseModel == Flat) mParResp.at(0) = mRandom->Uniform(1.5);
    else if(mResponseModel == Exp)  mParResp.at(0) = mRandom->Exp(0.5);
    else if(mResponseModel == Slope) {
      double u1 = mRandom->Uniform(2);
      double u2 = mRandom->Uniform(2);
      mParResp.at(0) = 2. - std::max(u1,u2);
    }

    // Loop over jets and
    // split genjets into towers
    NobjTow = 0;
    for(int i = 0 ; i < NobjJet ; ++i) {
      // Create 4-momentum of genjet (massless)
      genjet[i].SetPtEtaPhiM(jetgenpt[i],jetgeneta[i],jetgenphi[i],0);

      // Split it into towers and set truth of towers
      int ntow = splitJet(genjet[i],ttowet,ttoweta,ttowphi,ttowid_eta,ttowid_phi);  

      // Reset jet and genjet 4-momenta. They will
      // be repopulated by sum of tower 4-momenta
      jet[i].SetPtEtaPhiM(0,0,0,0);
      genjet[i].SetPtEtaPhiM(0,0,0,0);

      // Generate pi0 fraction i.e. non-hadronic fraction
      double p0frac = mRandom->Uniform(mMaxPi0Frac);

      // Loop over towers and smear truth with response
      // and resolution
      for(int j = 0 ; j < ntow ; ++j) {
	int k         = NobjTow + j;

	// Copy tower truth to tower
	towet[k]      = ttowet[j];
	toweta[k]     = ttoweta[j];
	towphi[k]     = ttowphi[j];	
	towid[k]      = 0;
	towid_eta[k]  = ttowid_eta[j];
	towid_phi[k]  = ttowid_phi[j];
	tow_jetidx[k] = i;
	tower.SetPtEtaPhiM(towet[k],toweta[k],towphi[k],0);
	towen[k]      = tower.E();

	// Add unsmeared tower to genjet
	genjet[i]    += tower;

	// Calculate smear factor
	bool calcSmearFactor = false;
	if( j == 0 || mSmearTowersIndividually ) calcSmearFactor = true;

	// Smear hadronic par of tower energy
	smearTower((1 - p0frac) * tower.E(),calcSmearFactor,towen[k],towem[k],towhd[k],
		   towoe[k],towemtrue[k],towhdtrue[k],towoetrue[k]); 

	// Add remaining em part
	towen[k]     += p0frac * tower.E();
	towem[k]     += p0frac * tower.E();
	towemtrue[k] += p0frac * tower.E();

	// Multiply tower 4-momentum with response
	tower        *= towen[k]/tower.E();
	towet[k]      = tower.Pt();

	// Add tower to jet
	jet[i]       += tower;
      } // End loop over towers

      // Set jet and genjet variables
      NobjTow     += ntow; 
      jetpt[i]     = jet[i].Pt();
      jetphi[i]    = jet[i].Phi();
      jeteta[i]    = jet[i].Eta();
      jetet[i]     = jet[i].Pt();
      jete[i]      = jet[i].E(); 
      jetgenpt[i]  = genjet[i].Pt();
      jetgenphi[i] = genjet[i].Phi();
      jetgeneta[i] = genjet[i].Eta();
      jetgenet[i]  = genjet[i].Pt();
      jetgene[i]   = genjet[i].E();
    } // End of loop over jets

    // Check generated eta measurement
    // of first jet
    if((jeteta[0] < mMinEta) || (jeteta[0] > mMaxEta)) {
      --n;
      continue;
    }

    // Order jets by pt
    if((jeteta[1] > mMinEta) && (jeteta[1] < mMaxEta)
       && (jetpt[1] > jetpt[0])) {
      //swap jets
      jetpt[0]     = jet[1].Pt();
      jetphi[0]    = jet[1].Phi();
      jeteta[0]    = jet[1].Eta();
      jetet[0]     = jet[1].Pt();
      jete[0]      = jet[1].E(); 
      jetgenpt[0]  = genjet[1].Pt();
      jetgenphi[0] = genjet[1].Phi();
      jetgeneta[0] = genjet[1].Eta();
      jetgenet[0]  = genjet[1].Pt();
      jetgene[0]   = genjet[1].E();
      jetpt[1]     = jet[0].Pt();
      jetphi[1]    = jet[0].Phi();
      jeteta[1]    = jet[0].Eta();
      jetet[1]     = jet[0].Pt();
      jete[1]      = jet[0].E(); 
      jetgenpt[1]  = genjet[0].Pt();
      jetgenphi[1] = genjet[0].Phi();
      jetgeneta[1] = genjet[0].Eta();
      jetgenet[1]  = genjet[0].Pt();
      jetgene[1]   = genjet[0].E();

      // Swap towers
      for(int j = 0 ; j < NobjTow ; ++j) {
	tow_jetidx[j] = (tow_jetidx[j] == 0) ? 1 : 0;
      }
    }


    // Fill tree
    CalibTree->Fill();

    if(n % 1000 == 0) std::cout << "generating event " << n << '\n';
  }  // End of loop over events

  return CalibTree->GetEntriesFast();
}



//!  \brief Generate photon-jet events and write into tree
//!  \param filename Name of output file
//!  \param nevents Number of photon-jet events
//!  \return Number of generated events
//----------------------------------------------------------
int ToyMC::makePhotonJet(const char* filename, int nevents) {
  TFile* file = new TFile(filename,"recreate");
  TTree* CalibTree = new TTree("GammaJetTree","GammaJetTree");

  nevents = generatePhotonJetTree(CalibTree,nevents);
  //file->ls();
  file->Write();
  file->Close();
  return nevents;
}



//!  \brief Generate dijet events
//!  \param filename Name of output file
//!  \param nevents Number of dijet events
//!  \return Number of generated events
//----------------------------------------------------------
int ToyMC::makeDiJet(const char* filename, int nevents) {
  TFile* file = new TFile(filename,"recreate");
  TTree* CalibTree = new TTree("DiJetTree","DiJetTree");

  nevents = generateDiJetTree(CalibTree,nevents);
  //file->ls();
  file->Write();
  file->Close();
  return nevents;
}



//----------------------------------------------------------
void ToyMC::init(const std::string& configfile) {
  ConfigFile config(configfile.c_str());

  // Ranges
  mMinEta = config.read<double>("ToyMC min eta",-2.5);
  mMaxEta = config.read<double>("ToyMC max eta",2.5);
  mMinPt  = config.read<double>("ToyMC min pt",30);
  mMaxPt  = config.read<double>("ToyMC max pt",400);

  // Truth spectrum
  std::string spectrum = config.read<std::string>("ToyMC pt spectrum","uniform");
  if(spectrum == "powerlaw") {
    mPtSpectrum = PowerLaw; 
  } else if(spectrum == "uniform") {
    mPtSpectrum = Uniform;
  } else {
    std::cerr << "unknown ToyMC pt spectrum:" << spectrum << '\n';
    exit(1);
  }

  // Response model
  mParResp             = bag_of<double>(config.read<string>("ToyMC response parameters","1"));
  std::string response = config.read<std::string>("ToyMC response model","Constant");
  if        ( response == "Constant" ) {
    mResponseModel = Constant;
    assert( mParResp.size() >= 1 );
  } else if ( response == "Flat" ) {
    mResponseModel = Flat;
    assert( mParResp.size() >= 1 );
  } else if ( response == "Slope" ) {
    mResponseModel = Slope;
    assert( mParResp.size() >= 1 );
  } else if ( response == "Exp" ) {
    mResponseModel = Exp;
    assert( mParResp.size() >= 1 );
  } else if ( response == "L3" ) {
    mResponseModel = L3;
    assert( mParResp.size() >= 4 );
  } else if( response == "SimpleInverse" ) {
    mResponseModel = SimpleInverse;
    assert( mParResp.size() >= 2 );
  } else {
    std::cerr << "unknown ToyMC response model: " << response << '\n';
    exit(1);
  }

  // Resolution model
  mParReso               = bag_of<double>(config.read<string>("ToyMC resolution parameters","1.2 0.05"));
  std::string resolution = config.read<std::string>("ToyMC resolution model","Gauss");
  if(resolution == "Gauss") {
    mResolutionModel = Gauss;
    assert( mParReso.size() >= 2 );
  } else if(resolution  == "Landau") {
    mResolutionModel = Landau;
    assert( mParReso.size() >= 2 );
  } else if( resolution == "GaussUniform" ) {
    mResolutionModel = GaussUniform; 
    assert( mParReso.size() >= 3 );
  } else if( resolution == "TwoGauss" ) {
    mResolutionModel = TwoGauss;
    assert( mParReso.size() >= 4 );
    
    // Set up sum of two Gaussians as pdf
    double c  = mParReso.at(0);  // Normalization
    double u0 = 1.;              // Mean of central Gaussian (scale)
    double s0 = mParReso.at(1);  // Width of central Gaussian
    double u1 = mParReso.at(2);  // Mean of second Gaussian
    double s1 = mParReso.at(3);  // Width of central Gaussian

    double minResp = 0.;
    double maxResp = 2.;

    TF1 * f = new TF1("f","gaus(0)+gaus(3)",minResp,maxResp);
    f->SetParameter(0,c/sqrt(2*M_PI)/s0);
    f->SetParameter(1,u0);
    f->SetParameter(2,s0);
    f->SetParameter(3,(1.-c)/sqrt(2*M_PI)/s1);
    f->SetParameter(4,u1);
    f->SetParameter(5,s1);

    // Fill response histogram according to f
    mHistResp = new TH1F("hHistResp",";p^{jet}_{T} / p^{true}_{T};1/(Nw) dN / d(p^{jet}_{T} / p^{true}_{T})",			     1000,minResp,maxResp);
    for(int bin = 1; bin <= mHistResp->GetNbinsX(); bin++) {
      double r = f->Eval(mHistResp->GetBinCenter(bin));
      mHistResp->SetBinContent(bin,r);
    }

    double norm = mHistResp->Integral("width");
    mHistResp->Scale(1./norm);
    delete f;
   } else {
     std::cerr << "unknown ToyMC resolution model: " << resolution << '\n';
     exit(1);
   }

  // Calculate smear factor for each tower or each jet
  mSmearTowersIndividually = config.read<bool>("ToyMC smear towers individually",false);

  // Jets
  mJetSpreadA  = config.read<double>("ToyMC jet spread A",0.5);
  mJetSpreadB  = config.read<double>("ToyMC jet spread B",0);
   mNoOutOfCone = config.read<bool>("ToyMC avoid out-of-cone",true);
   mChunks      = config.read<int>("ToyMC chunks",200);
   mMaxPi0Frac  = config.read<double>("ToyMC max pi0 fraction",0.5);
   mMaxEmf      = config.read<double>("ToyMC tower max EMF",0.5);

  // General
  int seed = config.read<int>("ToyMC seed",0); 
  mRandom->SetSeed(seed);
  mType = config.read<int>("ToyMC type",1);
  if( !( mType == 1 || mType == 2 ) ) {
    std::cout << "unknown ToyMC event type " << mType << std::endl;
    exit(1);
  }
}



//----------------------------------------------------------
void ToyMC::print() const {
  std::cout << "\n  ToyMC configuration:\n";
  std::cout << " -----------------------------------------\n";
  std::cout << "  primary:      " << mMinEta << " < eta < " << mMaxEta << '\n';
  std::cout << "                " << mMinPt << " < pt < " << mMaxPt << "\n";

  std::cout << "  spectrum:     ";
  if( mPtSpectrum == Uniform )
    std::cout << "Uniform\n";
  else if( mPtSpectrum == PowerLaw )
    std::cout << "PowerLaw\n";

  std::cout << "  response:     ";
  if( mResponseModel == Constant )
    std::cout << "Constant\n";
  else if( mResponseModel == Flat )
    std::cout << "Flat\n";
  else if( mResponseModel == Exp )
    std::cout << "Exp\n";
  else if( mResponseModel == Slope )
    std::cout << "Slope\n";
  else if( mResponseModel == L3 )
    std::cout << "L3\n";
  else if( mResponseModel == SimpleInverse )
    std::cout << "SimpleInverse\n";

  std::cout << "  parameters:   ";
  for(unsigned int i = 0; i < mParResp.size(); i++)
    std::cout << mParResp.at(i) << ",  ";
  std::cout << "\n";

  std::cout << "  resolution:   ";
  if( mResolutionModel == Gauss )
    std::cout << "Gauss\n";
  else if( mResolutionModel == Landau )
    std::cout << "Landau\n";
  else if( mResolutionModel == GaussUniform )
    std::cout << "GaussUniform\n";
  else if( mResolutionModel == TwoGauss )
    std::cout << "TwoGauss\n";

  std::cout << "  parameters:   ";
  for(unsigned int i = 0; i < mParReso.size(); i++)
    std::cout << mParReso.at(i) << ",  ";
  std::cout << "\n";

  std::cout << "  Smear factor is calculated for each ";
  if( mSmearTowersIndividually ) std::cout << "tower\n";
  else                           std::cout << "jet\n";

  std::cout << "  max EMF:      " << mMaxEmf << "\n";
  std::cout << "  max pi0:      " << mMaxPi0Frac << "\n";

  std::cout << "  jet spread:   ";
  std::cout << "A = " << mJetSpreadA << ",  B = " << mJetSpreadB << "\n";
  std::cout << "  n chunks:     " << mChunks << "\n";
  std::cout << "  out-of-cone:  ";
  if( mNoOutOfCone )
    std::cout << "no\n";
  else 
    std::cout << "yes\n";

  std::cout << "  type:         ";
  if( mType == 1 )
    std::cout << "Photon-jet events\n";
  else if( mType == 2 )
    std::cout << "Dijet events\n";

  std::cout << "\n";
}
