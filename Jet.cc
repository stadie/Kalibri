//
//    Class for basic jets 
//
//    first version: Hartmut Stadie 2008/12/14
//    $Id: EventReader.h,v 1.1 2008/12/12 13:43:15 stadie Exp $
//   
#include "Jet.h"  


Jet::Jet(double Et, double EmEt, double HadEt ,double OutEt, double E,
	 double eta,double phi, Flavor flavor,
	 double const(*func)(TMeasurement *const x, double *const par),
	 double err,double* firstpar, int id, int npars)
  : TJet(Et,EmEt,HadEt,OutEt,E,eta,phi,flavor), par(firstpar), npar(npars), parid(id),
    error(err),f(func)
{
  temp = *this;
}

//varies the i'th parameter for this jet by eps and returns its overall 
// parameter id and sets the Et for the par + eps and par - eps result
int Jet::varyPar(int i, double eps, double Et, double &upperEt, double& lowerEt) 
{
  double orig = par[i];
  par[i] += eps;
  upperEt = correctedEt(Et);
  //std::cout << "id: " << parid+i << "par:" << orig << "," << par[i] << ", " << upperEt << '\n';
  par[i] = orig - eps;
  lowerEt = correctedEt(Et);
  par[i] = orig;
  return  parid + i;
}

double Jet::correctedEt(double Et) {
  double s = Et/pt;
  //scale jet properties accordingly
  temp.pt   = Et;
  temp.EMF  = EMF * s;
  temp.HadF = HadF * s;
  temp.OutF = OutF * s;
  temp.E    = TJet::E * s;
  return f(&temp,par);
}
