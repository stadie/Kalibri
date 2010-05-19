//
//    Class for jets with towers 
//
//    first version: Hartmut Stadie 2008/12/25
//    $Id: JetWithTowers.cc,v 1.25 2010/05/19 13:34:48 stadie Exp $
//   
#include"JetWithTowers.h"

#include "TLorentzVector.h"

JetWithTowers::JetWithTowers(double Et, double EmEt, double HadEt ,double OutEt, double E,
			     double eta,double phi, double phiphi, double etaeta, 
			     Flavor flavor, double genPt, double dR,
			     CorFactors* corFactors, const Function& func, 
			     double (*errfunc)(const double *x, const Measurement *xorig, double err), 
			     const Function& gfunc, double Etmin) 
  :  Jet(Et,EmEt,HadEt,OutEt,E,eta,phi,phiphi,etaeta,flavor,genPt,dR,corFactors,
	 func,errfunc,gfunc,Etmin),
     ntowerpars_(0)
{
}

JetWithTowers::JetWithTowers(const JetWithTowers& j) 
  : Jet(j),ntowerpars_(0)
{ 
  for(TowerCollConstIter i = j.towers_.begin() ; i != j.towers_.end() ; ++i) {
    const Tower *t = *i;
    addTower(t->Et(),t->EmEt(),t->HadEt(),t->OutEt(),t->E(),t->eta(),t->phi(),
	     t->f_,t->errf_);
  }
}

JetWithTowers::~JetWithTowers() 
{
  for(TowerCollIter i = towers_.begin() ; i != towers_.end() ; ++i) {
    delete *i;
  }
  towers_.clear();
}  

void JetWithTowers::changeParAddress(double* oldpar, double* newpar) 
{
  Jet::changeParAddress(oldpar,newpar);
  for(TowerCollIter i = towers_.begin() ; i != towers_.end() ; ++i) {
    (*i)->changeParAddress(oldpar,newpar);
  }
  for(std::map<int,double*>::iterator iter = towerpars_.begin() ;
      iter != towerpars_.end() ; ++iter) {
    iter->second += newpar - oldpar;
  }
}

double JetWithTowers::correctedEt(double Et,bool fast) const
{
  double ccet = EmEt() + OutEt();
  double HadEt = Et - ccet;
  if(! fast) {
    double chad = 0; 
    for(TowerCollConstIter i = towers_.begin() ; i != towers_.end() ; ++i) {
      chad += (*i)->projectionToJetAxis() * (*i)->correctedHadEt((*i)->HadEt());
    }  
    //std::cout << "jet ET:" << pt << " sum of tower:" << chad + EMF + OutF << "\n";
    for(TowerCollConstIter i = towers_.begin() ; i != towers_.end() ; ++i) {
      (*i)->setFractionOfJetHadEt((*i)->lastCorrectedHadEt()/chad);
      ccet += (*i)->projectionToJetAxis() * (*i)->correctedHadEt((*i)->fractionOfJetHadEt() * HadEt);
    }
  } else {
    for(TowerCollConstIter i = towers_.begin() ; i != towers_.end() ; ++i) {
      ccet += (*i)->projectionToJetAxis() * (*i)->correctedHadEt((*i)->fractionOfJetHadEt() * HadEt);
    }
  }
  //std::cout << "scale ET:" << Et << " cor. sum of tower:" << ccet << "\n";
  return Jet::correctedEt(ccet);
}

// varies all parameters for this jet by eps and returns a vector of the
// parameter id and the Et for the par + eps and par - eps variation
const Jet::VariationColl& JetWithTowers::varyPars(const double* eps, double Et, double start)
{
  Jet::varyPars(eps,Et,start);
  int i = Jet::nPar();

  for(std::map<int,double*>::const_iterator iter = towerpars_.begin() ;
      iter != towerpars_.end() ; ++iter) {
    double *p = iter->second;
    int id = iter->first;
    //std::cout << p << "," << id << " tower[0]:" << towers[0]->Par() << '\n';
    for(int towpar = 0 ; towpar < ntowerpars_ ; ++towpar) {
      //std::cout <<  "truth: " << Et << " alternating par:" << id + towpar << "  = " << p[towpar] 
      //		<< std::endl;
      double orig = p[towpar]; 
      p[towpar] += eps[id + towpar];
      varcoll_[i].upperEt = Jet::expectedEt(Et,start,varcoll_[i].upperError);
      if( varcoll_[i].upperEt < 0) varcoll_[i].upperEt = Measurement::pt;
      p[towpar] = orig - eps[id + towpar];
      varcoll_[i].lowerEt = Jet::expectedEt(Et,start,varcoll_[i].lowerError); 
      if( varcoll_[i].lowerEt < 0) varcoll_[i].lowerEt = Measurement::pt;
      p[towpar] = orig;
      varcoll_[i].parid = id + towpar;
      ++i;
    }
  }
  //std::cout << i << " parameters modified.\n";
  return varcoll_;
}
// varies all parameters for this jet by eps and returns a vector of the
// parameter id and the Et for the par + eps and par - eps variation
const Jet::VariationColl& JetWithTowers::varyParsDirectly(const double* eps, bool computeDeriv)
{
  Jet::varyParsDirectly(eps,computeDeriv);
  int i = Jet::nPar();
  
  const double deltaE = 0.1;
  for(std::map<int,double*>::const_iterator iter = towerpars_.begin() ;
      iter != towerpars_.end() ; ++iter) {
    double *p = iter->second;
    int id = iter->first;
    //std::cout << p << "," << id << " tower[0]:" << towers[0]->Par() << '\n';
    for(int towpar = 0 ; towpar < ntowerpars_ ; ++towpar) {
      //std::cout <<  " alternating par:" << id + towpar << "  = " << p[towpar] 
      //		<< ", " << p << std::endl;
      double orig = p[towpar]; 
      p[towpar] += eps[id + towpar];
      varcoll_[i].upperEt = correctedEt(Measurement::pt);
      varcoll_[i].upperError = expectedError(varcoll_[i].upperEt);
      if(computeDeriv) {
	varcoll_[i].upperEtDeriv =  (correctedEt(Measurement::pt+deltaE) -  correctedEt(Measurement::pt-deltaE))/2/deltaE;
      }
      p[towpar] = orig - eps[id + towpar];
      varcoll_[i].lowerEt = correctedEt(Measurement::pt); 
      varcoll_[i].lowerError = expectedError(varcoll_[i].lowerEt);
      if(computeDeriv) {
	varcoll_[i].lowerEtDeriv =  (correctedEt(Measurement::pt+deltaE) -  correctedEt(Measurement::pt-deltaE))/2/deltaE;
      }
      p[towpar] = orig;
      varcoll_[i].parid = id + towpar;
      //std::cout << "up:" << varcoll_[i].upperEt << " low:" << varcoll_[i].lowerEt << '\n'; 
      ++i;
    }
  }
  //std::cout << i << " parameters modified.\n";
  return varcoll_;
}

double JetWithTowers::error() const {
  double var = 0, err;
  for(TowerCollConstIter i = towers_.begin() ; i != towers_.end() ; ++i) {
    err = (*i)->projectionToJetAxis() * (*i)->error();
    var += err * err;
  }   
  double jeterr = Jet::error();
  return sqrt(var + jeterr * jeterr);
}

double JetWithTowers::expectedError(double et) const
{
  double var = 0, err;
  double HadEt = et - EmEt() - OutEt();
  //std::cout << "hadET:" << HadEt << '\n';
  for(TowerCollConstIter i = towers_.begin() ; i != towers_.end() ; ++i) {
    if((*i)->fractionOfJetHadEt() == 0) continue;
    //std::cout << "tower fraction:" << (*i)->fractionOfJetHadEt() << "   tower error:" << (*i)->expectedError((*i)->fractionOfJetHadEt() * HadEt + (*i)->EmEt() + (*i)->OutEt()) << '\n';
    err = (*i)->projectionToJetAxis() * (*i)->expectedError((*i)->fractionOfJetHadEt() * HadEt + (*i)->EmEt() + (*i)->OutEt());
    //assert(err == err);
    var += err * err;
  }
  double jeterr = Jet::expectedError(et);
  return sqrt(var + jeterr * jeterr);
}
  
void JetWithTowers::addTower(double Et, double EmEt, double HadEt ,
			     double OutEt, double E,double eta,double phi,
			     const Function& func,
			     double (*errfunc)(const double *x, const Measurement *xorig, double err))
{
  TLorentzVector jet, towp;
  double en = Measurement::HadF + Measurement::EMF + Measurement::OutF;
  en *= Measurement::E/Measurement::pt;
  double m = en > Measurement::E ? sqrt(en * en - Measurement::E * Measurement::E) : 0;
  //std::cout << "mass:" << m << '\n';
  jet.SetPtEtaPhiM(Measurement::pt,Measurement::eta,Measurement::phi,m);
  towp.SetPtEtaPhiM(Et,eta,phi,0);
  jet -= towp;
  double projection = (Measurement::pt - jet.Pt())/Et;
  //projection = 1.0;
  //std::cout << "Jet:" << jet.Eta() << ", " << jet.Phi()
  //	    << " tower:" << towp.Eta() << ", " << towp.Phi() << " :" 
  //	    << projection << std::endl;
  towers_.push_back(new Tower(Et,EmEt,HadEt,OutEt,E,eta,phi,projection,func,errfunc)); 
  ntowerpars_ = func.nPars();
  //std::cout << func.parIndex() << ", " << func.firstPar() << ", " <<  ntowerpars_ << '\n';
  towerpars_[func.parIndex()] = func.firstPar();
  varcoll_.resize(Jet::nPar() + towerpars_.size() * ntowerpars_);
  //std::cout << "tower:" << Et << " projected:" << projection * Et << std::endl;
}



JetWithTowers::Tower::Tower(double Et, double EmEt, double HadEt ,
			    double OutEt, double E,double eta,double phi, 
			    double alpha,const Function& func,
			    double (*errfunc)(const double *x, const Measurement *xorig, double err))
  :  Measurement(Et,EmEt,HadEt,OutEt,E,eta,phi), alpha_(alpha), lastCorHadEt_(0),
     fraction_(0), f_(func), errf_(errfunc)
{ 
  temp_ = *this;
}

JetWithTowers::Tower::Tower(double Et, double EmEt, double HadEt ,double OutEt, double E,
			    double EmEttrue, double HadEttrue, double OutEttrue,
			    double eta,double phi, double alpha, const Function& func,
			    double (*errfunc)(const double *x, const Measurement *xorig, double err))
  :  Measurement(Et,EmEt,HadEt,OutEt,E,eta,phi), alpha_(alpha), error_(0),
     mEmEttrue_(EmEttrue), mHadEttrue_(HadEttrue), mOutEttrue_(OutEttrue),
     lastCorHadEt_(0), fraction_(0), f_(func), errf_(errfunc)
{ 
  mEttrue_ = mEmEttrue_ + mHadEttrue_ + mOutEttrue_;
  temp_ = *this;
}


double JetWithTowers::Tower::correctedHadEt(double HadEt) const
{
  //assume that only the hadronic energy gets modified!
  temp_.pt   = HadEt + OutF + EMF;  
  temp_.HadF = HadEt;
  temp_.E    = Measurement::E * temp_.pt/pt;
  lastCorHadEt_ = f_(&temp_) - OutF - EMF;
  if(lastCorHadEt_ < 0) lastCorHadEt_ = 0;
  //std::cout << temp.HadF << ", " << HadEt << ":"  << lastCorHadEt << " par:" << f.firstPar() << std::endl;
  return lastCorHadEt_;
}
  
