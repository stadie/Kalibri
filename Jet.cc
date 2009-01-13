//
//    Class for basic jets 
//
//    first version: Hartmut Stadie 2008/12/14
//    $Id: Jet.cc,v 1.6 2009/01/09 18:09:58 stadie Exp $
//   
#include "Jet.h"  


Jet::Jet(double Et, double EmEt, double HadEt ,double OutEt, double E,
	 double eta,double phi, Flavor flavor,
	 double const(*func)(TMeasurement *const x, double *const par),
	 double err,double* firstpar, int id, int npars)
  : TJet(Et,EmEt,HadEt,OutEt,E,eta,phi,flavor), par(firstpar), npar(npars), parid(id),
    error(err),f(func),varcoll(npars)
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
 
// varies all parameters for this jet by eps and returns a vector of the
// parameter id and the Et for the par + eps and par - eps variation
const Jet::VariationColl& Jet::varyPars(double eps, double Et, double scale)
{
  //scale = Et;
  double s = scale;
  for(int i = 0 ; i < npar ; ++i) {
    double orig = par[i];
    par[i] += eps;
    varcoll[i].upperEt = expectedEt(Et,s,true);
    //varcoll[i].upperEt = expectedEt(Et,s,false);
    s = scale;
    par[i] = orig - eps;;
    varcoll[i].lowerEt = expectedEt(Et,s,true);
    //varcoll[i].lowerEt = expectedEt(Et,s,false);
    s = scale;
    par[i] = orig;
    varcoll[i].parid = parid + i;
  }
  return varcoll;
}


double Jet::correctedEt(double Et, bool fast) const {
  //assume that only the hadronic energy gets modified!
  temp.pt   = Et;  
  temp.HadF = Et - OutF - EMF;
  temp.E    = TJet::E * Et/pt;
  return f(&temp,par);
}

double Jet::expectedEt(double truth, double& scale, bool extrapolate)
{
  static const double eps = 1.0e-5;
  const double up = 4 * truth;
  const double low = 0.2 * truth;
  //find root of truth - jet->correctedEt(expectedEt)
  // x: expectedEt
  // y: truth -  jet->correctedEt(expectedEt)
  // f: jet->correctedEt(expectedEt)
  double x1 = scale;
  double f1 = correctedEt(x1,extrapolate);
  //get second point assuming a constant correction factor
  double x2 = (truth - EMF - OutF) * (x1 - EMF - OutF)/(f1 - EMF - OutF) + EMF;
  if((x2 > up )||(x2 < low)) x2 = scale;
  double f2 = correctedEt(x2,true);
  double y2 = truth - f2;
  //std::cout << "truth:" << truth << " scale:" << scale << "  f1:" << f1 << " x2:" << x2 
  //	    << " f2:" << f2 << '\n';
  if(extrapolate || (std::abs(y2) < eps)) {
    //std::cout << "extrapolated:" << x2 << ", " << y2 << " at scale " << scale << std::endl;  
    return x2;
  }
  if(! secant(truth,x1,x2,eps)) return -1;
  scale = x2;
  f1 = correctedEt(scale,true);
  x2 = (truth - EMF - OutF) * (scale - EMF - OutF)/(f1 - EMF - OutF) + EMF + OutF;
  //std::cout << i << ": scale:" << scale << ", expected:" << (truth - EMF) * (scale - EMF)/(f1 - EMF) + EMF << "  dist for scale:" << truth - f1 << "\n";
  //assert(std::abs(correctedEt(x2)-truth)/truth < eps); 
  return ((x2 < up )&&(x2 > low)) ? x2 : scale;
}

bool Jet::falseposition(double truth, double& x1, double& x2,double eps)
{
  //x2 is the best estimate!
  double f1 = correctedEt(x1,true);
  double f2 = correctedEt(x2,true);
  double temp;
  double step = 0.1 * truth;
  ++ncalls;
  if(x1 > x2) {
    temp = x1;
    x1 = x2;
    x2 = temp;
    temp = f1;
    f1 = f2;
    f2 = temp;
  }  
  double y2 = truth - f2;
  double y1 = truth - f1;
  //std::cout << "x1,2:" << x1 << "," << x2 << " f1,2:" << f1 << ", " << f2 
  //	    << " y1,2:" << y1 << ", " << y2 << std::endl;
  int i = 0;
  while(y1 * y2 > 0) {
    //std::cout << "x1,2:" << x1 << "," << x2 << " f1,2:" << f1 << ", " << f2 
    //      << " y1,2:" << y1 << ", " << y2 << std::endl;
    if(f1 > truth) {
      x1 -= step;
      f1 = correctedEt(x1,true);
      y1 = truth - f1;
      ++ntries;
    }
    if(f2 < truth) {
      x2 += step;
      f2 = correctedEt(x2,true);
      y2 = truth - f2;
      ++ntries;
    }
    ++i;
    if(i > 5) break;
  } 
  i = 0;
  while(std::abs((x2-x1)/x1) > eps) {
    //std::cout << i << ":" << x1 << ", " << x2 << " : " << y1 << ", " << y2 << std::endl;
    double x3 = x1 + y1 * (x2-x1)/(f2 - f1);
    double f3 = correctedEt(x3,true);
    double y3 = truth - f3;
    ++i;
    if(y1 * y3 < 0) {
      x2 = x3;
      f2 = f3;
      y2 = y3;
    } else {
      x1 = x3;
      f1 = f3;
      y1 = y3;
    }
    if(i > 100) {
      ++nfails;
      ntries += i;
      return false;
    }
  } 
  ntries += i;
  x2 = 0.5*(x1 + x2);
  return true;
}


bool Jet::secant(double truth, double& x1, double& x2,double eps)
{
  //x2 is the best estimate!
  const double up = 4 * truth;
  const double low = 0.2 * truth;
  double f1 = correctedEt(x1,true);
  double f2 = correctedEt(x2,true);
  double y2 = truth - f2;
  double y1 = truth - f1;
  ++ncalls;
  int i = 0;
  double dx = std::abs(x1-x2), dx1 = truth, dx2 = truth;
  while((dx/x1 > eps)&&(i < 100)) {
    //std::cout << i << ":" << x1 << ", " << x2 << " : " << y1 << ", " << y2 << std::endl;
    double x3 = x1 + y1 * (x2-x1)/(f2 - f1);
    if((x3 > up )||(x3 < low)) {
      x3 = 0.5 *(x1 + x2);
    }
    dx2 = dx1;
    dx1 = dx;
    dx = std::abs(x2 - x3);
    if(dx2 < dx) {
      //std::cout << "Warning: fit alternating!\n";
      //std::cout << i << ": last three intervall sizes " << dx << ", " 
      //		<< dx1 << ", " << dx2 << std::endl;
      ++nwarns;
      x3 = 0.5 * (x2+x1);
      x1 = x2;
      x2 = x3;
      dx = truth;
      dx1 = truth;
      dx2 = truth;
      ++i;
      continue;
    }
    double f3 = correctedEt(x3,true);
    double y3 = truth - f3;
    //use false position if root is bracketed
    if(y1 * y3 < 0) {
      x2 = x3;
      f2 = f3;
      y2 = y3;
    } else {
      x1 = x2;
      f1 = f2;
      y1 = y2;
      x2 = x3;
      f2 = f3;
      y2 = y3;
    }
    ++i;
    //std::cout << i << ":" << x1 << ", " << x2 << ":" << truth - f1 << "; " << truth - f2 << "\n";
    //     if(i > 100) {
    //       //std::cout << "failed to find good root\n";
    //       //std::cout << i << ":" << x1 << ", " << x2 << ":" << truth - f2 << "\n";
    //       x2 = 0.5 * (x2+x1);
    //       ++nfails;
    //       ntries += i;
    //       return false;
    //     } 
  } 
  ntries += i;
  if(std::abs(y2) > 0.001 * truth) {
    //std::cout << "failed to find good root\n";
    //std::cout << i << ":" << x1 << ", " << x2 << ":" << truth - f2 << "\n";
    ++nfails;
    return false;
  }
  return true;
}

int Jet::ncalls = 0;
int Jet::ntries = 0;
int Jet::nfails = 0;
int Jet::nwarns = 0;

void Jet::printInversionStats()
{
  if(ncalls) {
    std::cout << "Inversion statistics for expectedEt:\n";
    std::cout << "calls: " << ncalls << " average number of iterations:"
	      << (double)ntries/ncalls << " failures:" << (double)nfails/ntries
	      << "%    warnings:" << (double)nwarns/ntries << "%" <<std::endl;
  }
}
