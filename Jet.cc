//
//    Class for basic jets 
//
//    first version: Hartmut Stadie 2008/12/14
//    $Id: Jet.cc,v 1.2 2008/12/17 09:37:18 stadie Exp $
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
int Jet::varyPar(int i, double eps, double Et, double scale, double &upperEt, double& lowerEt) 
{
  double orig = par[i];
  par[i] += eps;
  upperEt = expectedEt(Et,scale,true);
  par[i] = orig - eps;
  lowerEt = expectedEt(Et,scale,true);
  par[i] = orig;
  return  parid + i;
}

double Jet::correctedEt(double Et) const {
  //assume that only the hadronic energy gets modified!
  temp.pt   = Et;  
  temp.HadF = Et - OutF - EMF;
  temp.E    = TJet::E * Et/pt;
  return f(&temp,par);
}

double Jet::expectedEt(double truth, double& scale, bool extrapolate)
{
  static const double eps = 1.0e-4;
  const double up = 4 * truth;
  const double low = 0.2 * truth;
  //find root of truth - jet->correctedEt(expectedEt)
  // x: expectedEt
  // y: truth -  jet->correctedEt(expectedEt)
  // f: jet->correctedEt(expectedEt)
  double x1 = scale;
  double f1 = correctedEt(x1);
  //get second point assuming a constant correction factor
  double x2 = (truth - EMF) * (x1 - EMF)/(f1 - EMF) + EMF;
  double f2 = correctedEt(x2);
  double y2 = truth - f2;
  if(extrapolate || (std::abs(y2) < eps)) {
    //std::cout << "extrapolated:" << x2 << ", " << y2 << " at scale " << scale << std::endl;
    return ((x2 < up )&&(x2 > low)) ? x2 : scale;
  }
  //use bisection method
  ++ncalls;
  double y1 = truth - f1;
  int i = 0;
  double dx = std::abs(x1-x2), dx1 = scale, dx2 = scale;
  while(dx > eps) {
    //std::cout << i << ":" << x1 << ", " << x2 << " : " << y1 << ", " << y2 << std::endl;
    double x3 = x1 + y1 * (x2-x1)/(f2 - f1);
    if((x3 > up )||(x3 < low)) {
      x3 = 0.5 *(x1 + x2);
    }
    dx2 = dx1;
    dx1 = dx;
    dx = std::abs(x2 - x3);
    if((i > 0) && (dx2 < dx)) {
      //std::cout << "Warning: fit alternating!\n";
      //std::cout << i << ": last three intervall sizes " << dx << ", " 
      //		<< dx1 << ", " << dx2 << std::endl;
      ++nwarns;
      x3 = 0.5 * (x2+x1);
      x1 = x2;
      x2 = x3;
      dx = scale;
      dx1 = scale;
      dx2 = scale;
      ++i;
      continue;
    }
    x1 = x2;
    f1 = f2;
    x2 = x3;
    f2 = correctedEt(x2);
    y1 = truth - f1;
    y2 = truth - f2;
    ++i;
    //std::cout << i << ":" << x1 << ", " << x2 << ":" << truth - f1 << "; " << truth - f2 << "\n";
    if(i > 10000) {
      //std::cout << "failed to find good root\n";
      //std::cout << i << ":" << x1 << ", " << x2 << ":" << truth - f2 << "\n";
      x2 = 0.5 * (x2+x1);
      ++nfails;
      break;
    }
  }
  ntries += i;
  scale = x2;
  f1 = correctedEt(scale);
  x2 = (truth - EMF) * (scale - EMF)/(f1 - EMF) + EMF;
  //std::cout << i << ": scale:" << scale << ", expected:" << (truth - EMF) * (scale - EMF)/(f1 - EMF) + EMF << "  dist for scale:" << truth - f1 << "\n";
  return ((x2 < up )&&(x2 > low)) ? x2 : scale;
}

int Jet::ncalls = 0;
int Jet::ntries = 0;
int Jet::nfails = 0;
int Jet::nwarns = 0;

void Jet::printInversionStats()
{
  std::cout << "Inversion statistics for expectedEt:\n";
  std::cout << "calls: " << ncalls << " average number of iterations:"
	    << (double)ntries/ncalls << " failures:" << (double)nfails/ntries
	    << "%    warnings:" << (double)nwarns/ntries << "%" <<std::endl;

}
