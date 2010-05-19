//
//    Class for jets with tracks 
//
//    first version: Hartmut Stadie 2009/04/08
//    $Id: JetWithTracks.cc,v 1.12 2010/05/19 13:34:48 stadie Exp $
//   
#include"JetWithTracks.h"

#include "TLorentzVector.h"

JetWithTracks::JetWithTracks(double Et, double EmEt, double HadEt ,double OutEt, double E,
			     double eta,double phi, double phiphi, double etaeta, Flavor flavor, double genPt, double dR,
			     CorFactors* corFactors, const Function& func, 
			     double (*errfunc)(const double *x, const Measurement *xorig, double err), 
			     const Function& gfunc, double Etmin) 
  :  Jet(Et,EmEt,HadEt,OutEt,E,eta,phi,phiphi,etaeta,flavor,genPt,dR,corFactors,
	 func,errfunc,gfunc,Etmin),
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
	     t->TrackQualityT,t->MuDR,t->MuDE,t->Efficiency,t->f_,t->errf_);
  }
}

JetWithTracks::~JetWithTracks() 
{
  for(TrackCollIter i = tracks_.begin() ; i != tracks_.end() ; ++i) {
    delete *i;
  }
}  

void JetWithTracks::changeParAddress(double* oldpar, double* newpar) 
{
  Jet::changeParAddress(oldpar,newpar);
  for(TrackCollIter i = tracks_.begin() ; i != tracks_.end() ; ++i) {
    (*i)->changeParAddress(oldpar,newpar);
  }
  for(std::map<int,double*>::iterator iter = trackpars_.begin() ;
      iter != trackpars_.end() ; ++iter) {
    iter->second += newpar - oldpar;
  }
}

double JetWithTracks::correctedEt(double Et,bool fast) const
{
  const double ConeRadius = 0.5; 
  const double MIPsignal = 4;   
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
  double res = (Et > expectedCaloEt_) ? Jet::correctedEt(Et -  expectedCaloEt_) : 0;
  //std::cout << "Et:" << Et << " track corrected:" << ccet << " jet cor:" << res << '\n';
  return res + trackPt_;
}

// varies all parameters for this jet by eps and returns a vector of the
// parameter id and the Et for the par + eps and par - eps variation
const Jet::VariationColl& JetWithTracks::varyPars(const double* eps, double Et, double start)
{
  Jet::varyPars(eps,Et,start);
  int i = Jet::nPar();

  for(std::map<int,double*>::const_iterator iter = trackpars_.begin() ;
      iter != trackpars_.end() ; ++iter) {
    double *p = iter->second;
    int id = iter->first;
    //std::cout << p << "," << id << " track[0]:" << tracks[0]->Par() << '\n';
    for(int trkpar = 0 ; trkpar < ntrackpars_ ; ++trkpar) {
      //std::cout <<  "truth: " << Et << " alternating par:" << id + trkpar << "  = " << p[trkpar] << std::endl;
      double orig = p[trkpar]; 
      p[trkpar] += eps[id + trkpar];
      varcoll_[i].upperEt = Jet::expectedEt(Et,start,varcoll_[i].upperError);
      p[trkpar] = orig - eps[id + trkpar];
      varcoll_[i].lowerEt = Jet::expectedEt(Et,start,varcoll_[i].lowerError); 
      p[trkpar] = orig;
      varcoll_[i].parid = id + trkpar;
      ++i;
    }
  }
  //std::cout << i << " parameters modified.\n";
  return varcoll_;
}
// varies all parameters for this jet by eps and returns a vector of the
// parameter id and the Et for the par + eps and par - eps variation
const Jet::VariationColl& JetWithTracks::varyParsDirectly(const double* eps, bool computeDeriv)
{
  Jet::varyParsDirectly(eps,computeDeriv);
  int i = Jet::nPar();

  const double deltaE = 0.1;
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
      varcoll_[i].upperEt = correctedEt(Measurement::pt);
      varcoll_[i].upperError = expectedError(varcoll_[i].upperEt);
      if(computeDeriv) {
	varcoll_[i].upperEtDeriv =  (correctedEt(Measurement::pt+deltaE) -  correctedEt(Measurement::pt-deltaE))/2/deltaE;
      }
      p[trkpar] = orig - eps[id + trkpar];
      varcoll_[i].lowerEt = correctedEt(Measurement::pt); 
      varcoll_[i].lowerError = expectedError(varcoll_[i].lowerEt);
      if(computeDeriv) {
	varcoll_[i].lowerEtDeriv =  (correctedEt(Measurement::pt+deltaE) -  correctedEt(Measurement::pt-deltaE))/2/deltaE;
      }      
      p[trkpar] = orig;
      varcoll_[i].parid = id + trkpar;
      //std::cout << "up:" << varcoll_[i].upperEt << " low:" << varcoll_[i].lowerEt << '\n'; 
      ++i;
    }
  }
  //std::cout << i << " parameters modified.\n";
  return varcoll_;
}

double JetWithTracks::error() const {
  double var = 0, err;
  for(TrackCollConstIter i = tracks_.begin() ; i != tracks_.end() ; ++i) {
    err = (*i)->error();
    var += err * err;
  }   
  double jeterr = Jet::error();
  return sqrt(var + jeterr * jeterr);
}

double JetWithTracks::expectedError(double et) const
{
  double var = 0, err;
  for(TrackCollConstIter i = tracks_.begin() ; i != tracks_.end() ; ++i) {
    err = (*i)->error();
    var += err * err;
  }
  double jeterr = Jet::expectedError(et);
  return sqrt(var + jeterr * jeterr);
}
  
void JetWithTracks::addTrack(double Et, double EmEt, double HadEt ,
			     double OutEt, double E,double eta,double phi,
			     int TrackId, int TowerId, double DR, double DRout,
			     double etaOut, double phiOut, double EM1, double EM5, 
			     double Had1, double Had5, double TrackChi2, 
			     int NValidHits, bool TrackQualityT, double MuDR, 
			     double MuDE, double Efficiency, const Function& func,
			     double (*errfunc)(const double *x, const Measurement *xorig, double err))
{
  tracks_.push_back(new Track(Et,EmEt,HadEt,OutEt,E,eta,phi,TrackId,TowerId,DR,DRout,etaOut,phiOut,
			     EM1,EM5,Had1,Had5,TrackChi2,NValidHits,TrackQualityT,MuDR,MuDE,Efficiency,
			     func,errfunc)); 
  ntrackpars_ = func.nPars();
  trackpars_[func.parIndex()] = func.firstPar();
  varcoll_.resize(Jet::nPar() + trackpars_.size() * ntrackpars_);
}



JetWithTracks::Track::Track(double Et, double EmEt, double HadEt ,
			    double OutEt, double E,double eta,double phi,
			    int TrackId, int TowerId, double DR, double DRout,
			    double etaOut, double phiOut, double EM1, double EM5, 
			    double Had1, double Had5, double TrackChi2, 
			    int NValidHits, bool TrackQualityT, double MuDR, 
			    double MuDE, double Efficiency, const Function& func,
			    double (*errfunc)(const double *x, const Measurement *xorig, double err))
  :  TTrack(Et,EmEt,HadEt,OutEt,E,eta,phi,TrackId,TowerId,DR,DRout,etaOut,phiOut,EM1,EM5,Had1,Had5,
	    TrackChi2,NValidHits,TrackQualityT,MuDR,MuDE,Efficiency), f_(func), errf_(errfunc)
{ 
}

double JetWithTracks::Track::expectedEt() const
{
  double et = f_(this);
  assert(et == et);
  assert(et >= 0);
  return et;
}

  
double JetWithTracks::expectedEt(double truth, double start, bool fast)
{
  correctedEt(start);
  //std::cout << truth << ", " << pt << ", " << cet << ", " 
  //	    <<  expectedCaloEt << ", " << trackPt << "\n";
  if(truth < trackPt_) return (expectedCaloEt_ > 0) ? expectedCaloEt_ : 1.0;
  if(expectedCaloEt_ >=  Measurement::pt) return expectedCaloEt_;
  double m = Jet::expectedEt(truth, start, fast);
  //std::cout << "expected ET:" << m << '\n';
  //assert(m > 0);
  return m;
}
