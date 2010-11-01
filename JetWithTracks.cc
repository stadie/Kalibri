//
//    Class for jets with tracks 
//
//    first version: Hartmut Stadie 2009/04/08
//    $Id: JetWithTracks.cc,v 1.14 2010/10/20 11:28:12 stadie Exp $
//   
#include"JetWithTracks.h"

#include "Parameters.h"
#include "TLorentzVector.h"

JetWithTracks::JetWithTracks(float Et, float EmEt, float HadEt ,float OutEt, float E,
			     float eta,float phi, float phiphi, float etaeta, Flavor flavor, float genPt, float dR,
			     CorFactors* corFactors, const Function& func, 
			     float (*errfunc)(const float *x, const Measurement *xorig, float err), 
			     const Function& gfunc) 
  :  Jet(Et,EmEt,HadEt,OutEt,E,eta,phi,phiphi,etaeta,flavor,genPt,dR,corFactors,
	 func,errfunc,gfunc),
     ntrackpars_(0), expectedCaloEt_(0), trackPt_(0)
{
}

JetWithTracks::JetWithTracks(const JetWithTracks& j) 
  :  Jet(j),ntrackpars_(0), expectedCaloEt_(0), trackPt_(0)
{
  for(TrackCollConstIter i = j.tracks_.begin() ; i != j.tracks_.end() ; ++i) {
    const Track *t = *i;
    addTrack(t->Et(),t->EmEt(),t->HadEt(),t->OutEt(),t->E(),t->eta(),t->phi(),
	     t->TrackId,t->TowerId,t->DR,t->DRout,t->etaOut,t->phiOut,
	     t->EM1,t->EM5,t->Had1,t->Had5,t->TrackChi2,t->NValidHits,
	     t->TrackQualityT,t->MuDR,t->MuDE,t->Efficiency,*(t->f_),t->errf_);
  }
}

JetWithTracks::~JetWithTracks() 
{
  for(TrackCollIter i = tracks_.begin() ; i != tracks_.end() ; ++i) {
    delete *i;
  }
}  

void JetWithTracks::setParameters(Parameters* param) 
{
  Jet::setParameters(param);
  for(TrackCollIter i = tracks_.begin() ; i != tracks_.end() ; ++i) {
    const Function& f = (*i)->setParameters(param);
    trackpars_[f.parIndex()] = f.firstPar();
  }
}

float JetWithTracks::correctedEt(float Et,bool fast) const
{
  const float ConeRadius = 0.5; 
  const float MIPsignal = 4;   
  int NoUsedTracks = 0;
  bool isMuon;
  trackPt_ = 0;
  expectedCaloEt_ = 0;
  for(TrackCollConstIter i = tracks_.begin() ; i != tracks_.end() ; ++i) {
    Track *t = (*i);
    if(! t->goodTrack()) continue; 
    if( t->pt > 100) continue;

    ++NoUsedTracks;
    isMuon = ( t->trackId() == 13);
  
    if(t->dR() < ConeRadius) {
      if(t->dRout() < ConeRadius) {
	if(!isMuon) {
	  trackPt_ += t->pt;
	  expectedCaloEt_ += t->expectedEt();
	}
      } else { // out-of-cone
	trackPt_ += t->pt;
      }
    } else { // subtract curl-ins
      if(isMuon) expectedCaloEt_ -= MIPsignal;
      else expectedCaloEt_ -= t->expectedEt();
    }
  }
  //correct neutral rest
  float res = (Et > expectedCaloEt_) ? Jet::correctedEt(Et -  expectedCaloEt_) : 0;
  //std::cout << "Et:" << Et << " track corrected:" << ccet << " jet cor:" << res << '\n';
  return res + trackPt_;
}

// varies all parameters for this jet by eps and returns a vector of the
// parameter id and the Et for the par + eps and par - eps variation
const Parameters::VariationColl& JetWithTracks::varyPars(const double* eps, float Et, float start)
{
  Jet::varyPars(eps,Et,start);
  int i = Jet::nPar();
  Parameters::VariationColl& varcoll = parameters_->cachedVariationColl();
  varcoll.resize(nPar());
  for(std::map<int,double*>::const_iterator iter = trackpars_.begin() ;
      iter != trackpars_.end() ; ++iter) {
    double *p = iter->second;
    int id = iter->first;
    //std::cout << p << "," << id << " track[0]:" << tracks[0]->Par() << '\n';
    for(int trkpar = 0 ; trkpar < ntrackpars_ ; ++trkpar) {
      //std::cout <<  "truth: " << Et << " alternating par:" << id + trkpar << "  = " << p[trkpar] << std::endl;
      double orig = p[trkpar]; 
      p[trkpar] += eps[id + trkpar];
      varcoll[i].upperEt = Jet::expectedEt(Et,start,varcoll[i].upperError);
      p[trkpar] = orig - eps[id + trkpar];
      varcoll[i].lowerEt = Jet::expectedEt(Et,start,varcoll[i].lowerError); 
      p[trkpar] = orig;
      varcoll[i].parid = id + trkpar;
      ++i;
    }
  }
  //std::cout << i << " parameters modified.\n";
  return varcoll;
}
// varies all parameters for this jet by eps and returns a vector of the
// parameter id and the Et for the par + eps and par - eps variation
const Parameters::VariationColl& JetWithTracks::varyParsDirectly(const double* eps, bool computeDeriv)
{
  Jet::varyParsDirectly(eps,computeDeriv);
  int i = Jet::nPar();

  Parameters::VariationColl& varcoll = parameters_->cachedVariationColl();
  varcoll.resize(nPar());
  const float deltaE = 0.1;
  for(std::map<int,double*>::const_iterator iter = trackpars_.begin() ;
      iter != trackpars_.end() ; ++iter) {
    double *p = iter->second;
    int id = iter->first;
    //std::cout << p << "," << id << " track[0]:" << tracks[0]->Par() << '\n';
    for(int trkpar = 0 ; trkpar < ntrackpars_ ; ++trkpar) {
      //std::cout <<  " alternating par:" << id + trkpar << "  = " << p[trkpar] 
      //		<< ", " << p << std::endl;
      double orig = p[trkpar]; 
      p[trkpar] += eps[id + trkpar];
      varcoll[i].upperEt = correctedEt(Measurement::pt);
      varcoll[i].upperError = expectedError(varcoll[i].upperEt);
      if(computeDeriv) {
	varcoll[i].upperEtDeriv =  (correctedEt(Measurement::pt+deltaE) -  correctedEt(Measurement::pt-deltaE))/2/deltaE;
      }
      p[trkpar] = orig - eps[id + trkpar];
      varcoll[i].lowerEt = correctedEt(Measurement::pt); 
      varcoll[i].lowerError = expectedError(varcoll[i].lowerEt);
      if(computeDeriv) {
	varcoll[i].lowerEtDeriv =  (correctedEt(Measurement::pt+deltaE) -  correctedEt(Measurement::pt-deltaE))/2/deltaE;
      }      
      p[trkpar] = orig;
      varcoll[i].parid = id + trkpar;
      //std::cout << "up:" << varcoll[i].upperEt << " low:" << varcoll[i].lowerEt << '\n'; 
      ++i;
    }
  }
  //std::cout << i << " parameters modified.\n";
  return varcoll;
}

float JetWithTracks::error() const {
  float var = 0, err;
  for(TrackCollConstIter i = tracks_.begin() ; i != tracks_.end() ; ++i) {
    err = (*i)->error();
    var += err * err;
  }   
  float jeterr = Jet::error();
  return sqrt(var + jeterr * jeterr);
}

float JetWithTracks::expectedError(float et) const
{
  float var = 0, err;
  for(TrackCollConstIter i = tracks_.begin() ; i != tracks_.end() ; ++i) {
    err = (*i)->error();
    var += err * err;
  }
  float jeterr = Jet::expectedError(et);
  return sqrt(var + jeterr * jeterr);
}
  
void JetWithTracks::addTrack(float Et, float EmEt, float HadEt ,
			     float OutEt, float E,float eta,float phi,
			     int TrackId, int TowerId, float DR, float DRout,
			     float etaOut, float phiOut, float EM1, float EM5, 
			     float Had1, float Had5, float TrackChi2, 
			     int NValidHits, bool TrackQualityT, float MuDR, 
			     float MuDE, float Efficiency, const Function& func,
			     float (*errfunc)(const float *x, const Measurement *xorig, float err))
{
  tracks_.push_back(new Track(Et,EmEt,HadEt,OutEt,E,eta,phi,TrackId,TowerId,DR,DRout,etaOut,phiOut,
			     EM1,EM5,Had1,Had5,TrackChi2,NValidHits,TrackQualityT,MuDR,MuDE,Efficiency,
			     func,errfunc)); 
  ntrackpars_ = func.nPars();
  trackpars_[func.parIndex()] = func.firstPar();
}



JetWithTracks::Track::Track(float Et, float EmEt, float HadEt ,
			    float OutEt, float E,float eta,float phi,
			    int TrackId, int TowerId, float DR, float DRout,
			    float etaOut, float phiOut, float EM1, float EM5, 
			    float Had1, float Had5, float TrackChi2, 
			    int NValidHits, bool TrackQualityT, float MuDR, 
			    float MuDE, float Efficiency, const Function& func,
			    float (*errfunc)(const float *x, const Measurement *xorig, float err))
  :  TTrack(Et,EmEt,HadEt,OutEt,E,eta,phi,TrackId,TowerId,DR,DRout,etaOut,phiOut,EM1,EM5,Had1,Had5,
	    TrackChi2,NValidHits,TrackQualityT,MuDR,MuDE,Efficiency), f_(&func), errf_(errfunc)
{ 
}

const Function& JetWithTracks::Track::setParameters(Parameters* param) {
  f_ = &(param->function(*f_));
  return *f_;
}


float JetWithTracks::Track::expectedEt() const
{
  float et = (*f_)(this);
  assert(et == et);
  assert(et >= 0);
  return et;
}

  
float JetWithTracks::expectedEt(float truth, float start, bool fast)
{
  correctedEt(start);
  //std::cout << truth << ", " << pt << ", " << cet << ", " 
  //	    <<  expectedCaloEt << ", " << trackPt << "\n";
  if(truth < trackPt_) return (expectedCaloEt_ > 0) ? expectedCaloEt_ : 1.0;
  if(expectedCaloEt_ >=  Measurement::pt) return expectedCaloEt_;
  float m = Jet::expectedEt(truth, start, fast);
  //std::cout << "expected ET:" << m << '\n';
  //assert(m > 0);
  return m;
}
