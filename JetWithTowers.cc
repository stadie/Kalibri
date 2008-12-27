//
//    Class for jets with towers 
//
//    first version: Hartmut Stadie 2008/12/25
//    $Id: Jet.h,v 1.1 2008/12/16 15:21:27 stadie Exp $
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
    towerpars[(*i)->FirstPar()] = (*i)->Par();
  }
}

double JetWithTowers::correctedEt(double Et) const
{
  double cet = 0; 
  for(TowerCollConstIter i = towers.begin() ; i != towers.end() ; ++i) {
    cet += (*i)->projectionToJetAxis() * (*i)->correctedEt((*i)->Et());
  }
  //std::cout << "jet ET:" << pt << " sum of tower:" << cet << "\n";
  double ccet = 0;
  for(TowerCollConstIter i = towers.begin() ; i != towers.end() ; ++i) {
    ccet += (*i)->projectionToJetAxis() * 
      (*i)->correctedEt((*i)->lastCorrectedEt()/cet * Et);
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
  //std::cout << "alternating par:" << id + towpar << "  = " << p[towpar] << std::endl;
  double orig = p[towpar];
  p[towpar] += eps;
  upperEt = expectedEt(Et,scale,true);
  p[towpar] = orig - eps;
  lowerEt = expectedEt(Et,scale,true);
  p[towpar] = orig;
  return id + towpar; 
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
  Jet::npar = njetpars + towerpars.size() * ntowerpars;
  //std::cout << "tower:" << Et << " projected:" << projection * Et << std::endl;
}



JetWithTowers::Tower::Tower(double Et, double EmEt, double HadEt ,
			    double OutEt, double E,double eta,double phi, 
			    double alpha,double const(*func)(TMeasurement *const x, double *const par),
			    double err, double* firstpar, int id, int npars)
  :  TMeasurement(Et,EmEt,HadEt,OutEt,E,eta,phi), alpha(alpha), par(firstpar), 
     npar(npars), parid(id), error(err),f(func)
{ 
  temp = *this;
}

double JetWithTowers::Tower::correctedEt(double Et) const
{
  //assume that only the hadronic energy gets modified!
  temp.pt   = Et;  
  temp.HadF = Et - OutF - EMF;
  temp.E    = TMeasurement::E * Et/pt;
  lastCorEt = f(&temp,par);
  //std::cout << pt << ", " << Et << ":"  << lastCorEt << " par:" << par[0] << std::endl;
  return lastCorEt;
}
  
