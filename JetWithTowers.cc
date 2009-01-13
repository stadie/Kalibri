//
//    Class for jets with towers 
//
//    first version: Hartmut Stadie 2008/12/25
//    $Id: JetWithTowers.cc,v 1.3 2009/01/09 18:09:58 stadie Exp $
//   
#include"JetWithTowers.h"

#include "TLorentzVector.h"

JetWithTowers::JetWithTowers(double Et, double EmEt, double HadEt,
			     double OutEt, double E,double eta,double phi, 
			     Flavor flavor,
			     double const(*func)(TMeasurement *const x, double *const par),
			     double err, double* firstpar, int id, int npars)
  : Jet(Et,EmEt,HadEt,OutEt,E,eta,phi,flavor,func,err,firstpar,id,npars),
    njetpars(npars),ntowerpars(0)
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


//varies the i'th parameter for this jet by eps and returns its overall 
// parameter id and sets the Et for the par + eps and par - eps result
int JetWithTowers::varyPar(int i, double eps, double Et, double scale, 
			   double& upperEt, double& lowerEt)
{
  if(i < njetpars) return Jet::varyPar(i,eps,Et,scale,upperEt,lowerEt);
  int towid = (i - njetpars)/ntowerpars;
  int towpar = (i - njetpars) % ntowerpars;
  //std::cout << "I: " << i << " tow:" << towid << "; " << towpar  << std::endl;
  std::map<int,double*>::const_iterator iter = towerpars.begin();
  for(int i  = 0 ; i < towid ; ++i) ++iter;
  double *p = iter->second;
  int id = iter->first;
  //std::cout << "truth: " << Et << "alternating par:" << id + towpar << "  = " << p[towpar] << std::endl;
  double orig = p[towpar]; 
  p[towpar] += eps;
  double s = scale;
  upperEt = expectedEt(Et,s,true);
  p[towpar] = orig - eps;
  s = scale;
  lowerEt = expectedEt(Et,s,true);
  p[towpar] = orig;
  return id + towpar; 
}

// varies all parameters for this jet by eps and returns a vector of the
// parameter id and the Et for the par + eps and par - eps variation
const Jet::VariationColl& JetWithTowers::varyPars(double eps, double Et, double scale)
{
  Jet::varyPars(eps,Et,scale);
  int i = njetpars;

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
      double s = scale;
      varcoll[i].upperEt = expectedEt(Et,s,true);
      p[towpar] = orig - eps;
      s = scale;
      varcoll[i].lowerEt = expectedEt(Et,s,true);
      p[towpar] = orig;
      varcoll[i].parid = id + towpar;
      ++i;
    }
  }
  //std::cout << i << " parameters modified.\n";
  return varcoll;
}

  
void JetWithTowers::addTower(double Et, double EmEt, double HadEt ,
			     double OutEt, double E,double eta,double phi,
			     double const(*func)(TMeasurement *const x, double *const par),
			     double err, double* firstpar, int id, int npars)
{
  TLorentzVector jet, towp;
  jet.SetPtEtaPhiM(TMeasurement::pt,TMeasurement::eta,TMeasurement::phi,0);
  towp.SetPtEtaPhiM(Et,eta,phi,0);
  jet -= towp;
  double projection = (TMeasurement::pt-jet.Et())/Et;
  projection = 1.0;
  //std::cout << "Jet:" << jet.Eta() << ", " << jet.Phi()
  //	    << " tower:" << towp.Eta() << ", " << towp.Phi() << " :" 
  //	    << projection << std::endl;
  towers.push_back(new Tower(Et,EmEt,HadEt,OutEt,E,eta,phi,projection,func,
			     err,firstpar,id,npars)); 
  ntowerpars = npars;
  towerpars[id] = firstpar;
  varcoll.resize(njetpars + towerpars.size() * ntowerpars);
  //std::cout << "tower:" << Et << " projected:" << projection * Et << std::endl;
}



JetWithTowers::Tower::Tower(double Et, double EmEt, double HadEt ,
			    double OutEt, double E,double eta,double phi, 
			    double alpha,double const(*func)(TMeasurement *const x, double *const par),
			    double err, double* firstpar, int id, int npars)
  :  TMeasurement(Et,EmEt,HadEt,OutEt,E,eta,phi), alpha(alpha), par(firstpar), 
     npar(npars), parid(id), error(err), lastCorHadEt(0), fraction(0), f(func)
{ 
  temp = *this;
}

double JetWithTowers::Tower::correctedHadEt(double HadEt) const
{
  //assume that only the hadronic energy gets modified!
  temp.pt   = HadEt + OutF + EMF;  
  temp.HadF = HadEt;
  temp.E    = TMeasurement::E * temp.pt/pt;
  lastCorHadEt = f(&temp,par) - OutF - EMF;
  //std::cout << pt << ", " << Et << ":"  << lastCorEt << " par:" << par[0] << std::endl;
  return lastCorHadEt;
}
  
