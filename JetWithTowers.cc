//
//    Class for jets with towers 
//
//    first version: Hartmut Stadie 2008/12/25
//    $Id: JetWithTowers.cc,v 1.8 2009/02/12 19:38:51 stadie Exp $
//   
#include"JetWithTowers.h"

#include "TLorentzVector.h"

JetWithTowers::JetWithTowers(double Et, double EmEt, double HadEt,
			     double OutEt, double E,double eta,double phi, 
			     Flavor flavor, const Function& func,
			     double (*errfunc)(const double *x, const TMeasurement *xorig,double err), 
			     const Function& gfunc)
  : Jet(Et,EmEt,HadEt,OutEt,E,eta,phi,flavor,func,errfunc,gfunc),
    ntowerpars(0)
{
}

JetWithTowers::~JetWithTowers() 
{
  for(TowerCollIter i = towers.begin() ; i != towers.end() ; ++i) {
    delete *i;
  }
}  

void JetWithTowers::ChangeParAddress(double* oldpar, double* newpar) 
{
  Jet::ChangeParAddress(oldpar,newpar);
  for(TowerCollIter i = towers.begin() ; i != towers.end() ; ++i) {
    (*i)->ChangeParAddress(oldpar,newpar);
  }
  for(std::map<int,double*>::iterator iter = towerpars.begin() ;
      iter != towerpars.end() ; ++iter) {
    iter->second += newpar - oldpar;
  }
}

double JetWithTowers::correctedEt(double Et,bool fast) const
{
  double ccet = EmEt() + OutEt();
  double HadEt = Et - ccet;
  if(! fast) {
    double chad = 0; 
    for(TowerCollConstIter i = towers.begin() ; i != towers.end() ; ++i) {
      chad += (*i)->projectionToJetAxis() * (*i)->correctedHadEt((*i)->HadEt());
    }  
    //std::cout << "jet ET:" << pt << " sum of tower:" << chad + EMF + OutF << "\n";
    for(TowerCollConstIter i = towers.begin() ; i != towers.end() ; ++i) {
      (*i)->setFractionOfJetHadEt((*i)->lastCorrectedHadEt()/chad);
      ccet += (*i)->projectionToJetAxis() * (*i)->correctedHadEt((*i)->fractionOfJetHadEt() * HadEt);
    }
  } else {
    for(TowerCollConstIter i = towers.begin() ; i != towers.end() ; ++i) {
      ccet += (*i)->projectionToJetAxis() * (*i)->correctedHadEt((*i)->fractionOfJetHadEt() * HadEt);
    }
  }
  //std::cout << "scale ET:" << Et << " cor. sum of tower:" << ccet << "\n";
  return Jet::correctedEt(ccet);
}

// varies all parameters for this jet by eps and returns a vector of the
// parameter id and the Et for the par + eps and par - eps variation
const Jet::VariationColl& JetWithTowers::varyPars(double eps, double Et, double start)
{
  Jet::varyPars(eps,Et,start);
  int i = Jet::nPar();

  for(std::map<int,double*>::const_iterator iter = towerpars.begin() ;
      iter != towerpars.end() ; ++iter) {
    double *p = iter->second;
    int id = iter->first;
    //std::cout << p << "," << id << " tower[0]:" << towers[0]->Par() << '\n';
    for(int towpar = 0 ; towpar < ntowerpars ; ++towpar) {
      //std::cout <<  "truth: " << Et << " alternating par:" << id + towpar << "  = " << p[towpar] 
      //		<< std::endl;
      double orig = p[towpar]; 
      p[towpar] += eps;
      varcoll[i].upperEt = expectedEt(Et,start);
      if( varcoll[i].upperEt < 0) varcoll[i].upperEt = pt;
      varcoll[i].upperError = expectedError(varcoll[i].upperEt);
      p[towpar] = orig - eps;
      varcoll[i].lowerEt = expectedEt(Et,start); 
      if( varcoll[i].lowerEt < 0) varcoll[i].lowerEt = pt;
      varcoll[i].lowerError = expectedError(varcoll[i].lowerEt);
      p[towpar] = orig;
      varcoll[i].parid = id + towpar;
      ++i;
    }
  }
  //std::cout << i << " parameters modified.\n";
  return varcoll;
}
// varies all parameters for this jet by eps and returns a vector of the
// parameter id and the Et for the par + eps and par - eps variation
const Jet::VariationColl& JetWithTowers::varyParsDirectly(double eps)
{
  Jet::varyParsDirectly(eps);
  int i = Jet::nPar();

  for(std::map<int,double*>::const_iterator iter = towerpars.begin() ;
      iter != towerpars.end() ; ++iter) {
    double *p = iter->second;
    int id = iter->first;
    //std::cout << p << "," << id << " tower[0]:" << towers[0]->Par() << '\n';
    for(int towpar = 0 ; towpar < ntowerpars ; ++towpar) {
      //std::cout <<  " alternating par:" << id + towpar << "  = " << p[towpar] 
      //		<< ", " << p << std::endl;
      double orig = p[towpar]; 
      p[towpar] += eps;
      varcoll[i].upperEt = correctedEt(pt);
      varcoll[i].upperError = expectedError(varcoll[i].upperEt);
      p[towpar] = orig - eps;
      varcoll[i].lowerEt = correctedEt(pt); 
      varcoll[i].lowerError = expectedError(varcoll[i].lowerEt);
      p[towpar] = orig;
      varcoll[i].parid = id + towpar;
      //std::cout << "up:" << varcoll[i].upperEt << " low:" << varcoll[i].lowerEt << '\n'; 
      ++i;
    }
  }
  //std::cout << i << " parameters modified.\n";
  return varcoll;
}

double JetWithTowers::Error() const {
  double var = 0, err;
  for(TowerCollConstIter i = towers.begin() ; i != towers.end() ; ++i) {
    err = (*i)->projectionToJetAxis() * (*i)->Error();
    var += err * err;
  }   
  double jeterr = Jet::Error();
  return sqrt(var + jeterr * jeterr);
}

double JetWithTowers::expectedError(double truth) const
{
  double var = 0, err;
  double HadEt = truth - EmEt() - OutEt();
  //std::cout << "hadET:" << HadEt << '\n';
  for(TowerCollConstIter i = towers.begin() ; i != towers.end() ; ++i) {
    if((*i)->fractionOfJetHadEt() == 0) continue;
    //std::cout << "tower fraction:" << (*i)->fractionOfJetHadEt() << "   tower error:" << (*i)->expectedError((*i)->fractionOfJetHadEt() * HadEt + (*i)->EmEt() + (*i)->OutEt()) << '\n';
    err = (*i)->projectionToJetAxis() * (*i)->expectedError((*i)->fractionOfJetHadEt() * HadEt + (*i)->EmEt() + (*i)->OutEt());
    //assert(err == err);
    var += err * err;
  }
  double jeterr = Jet::expectedError(truth);
  return sqrt(var + jeterr * jeterr);
}
  
void JetWithTowers::addTower(double Et, double EmEt, double HadEt ,
			     double OutEt, double E,double eta,double phi,
			     const Function& func,
			     double (*errfunc)(const double *x, const TMeasurement *xorig, double err))
{
  TLorentzVector jet, towp;
  jet.SetPtEtaPhiM(TMeasurement::pt,TMeasurement::eta,TMeasurement::phi,0);
  towp.SetPtEtaPhiM(Et,eta,phi,0);
  jet -= towp;
  double projection = (TMeasurement::pt - jet.Pt())/Et;
  //projection = 1.0;
  //std::cout << "Jet:" << jet.Eta() << ", " << jet.Phi()
  //	    << " tower:" << towp.Eta() << ", " << towp.Phi() << " :" 
  //	    << projection << std::endl;
  towers.push_back(new Tower(Et,EmEt,HadEt,OutEt,E,eta,phi,projection,func,errfunc)); 
  ntowerpars = func.nPars();
  //std::cout << func.parIndex() << ", " << func.firstPar() << ", " <<  ntowerpars << '\n';
  towerpars[func.parIndex()] = func.firstPar();
  varcoll.resize(Jet::nPar() + towerpars.size() * ntowerpars);
  //std::cout << "tower:" << Et << " projected:" << projection * Et << std::endl;
}



JetWithTowers::Tower::Tower(double Et, double EmEt, double HadEt ,
			    double OutEt, double E,double eta,double phi, 
			    double alpha,const Function& func,
			    double (*errfunc)(const double *x, const TMeasurement *xorig, double err))
  :  TMeasurement(Et,EmEt,HadEt,OutEt,E,eta,phi), alpha(alpha), lastCorHadEt(0),
     fraction(0), f(func), errf(errfunc)
{ 
  temp = *this;
}

double JetWithTowers::Tower::correctedHadEt(double HadEt) const
{
  //assume that only the hadronic energy gets modified!
  temp.pt   = HadEt + OutF + EMF;  
  temp.HadF = HadEt;
  temp.E    = TMeasurement::E * temp.pt/pt;
  lastCorHadEt = f(&temp) - OutF - EMF;
  if(lastCorHadEt < 0) lastCorHadEt = 0;
  //std::cout << temp.HadF << ", " << HadEt << ":"  << lastCorHadEt << " par:" << f.firstPar() << std::endl;
  return lastCorHadEt;
}
  
