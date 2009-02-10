//
//    Class for basic jets 
//
//    first version: Hartmut Stadie 2008/12/14
//    $Id: Jet.cc,v 1.10 2009/01/22 15:30:30 stadie Exp $
//   
#include "Jet.h"  


Jet::Jet(double Et, double EmEt, double HadEt ,double OutEt, double E,
	 double eta,double phi, Flavor flavor,
	 double (*func)(const TMeasurement *x, const double *par),
	 double (*errfunc)(const double *x, const TMeasurement *xorig, double err), 
	 double* firstpar, int id, int npars)
  : TJet(Et,EmEt,HadEt,OutEt,E,eta,phi,flavor,0.0,1.0,1.0,1.0,1.0), par(firstpar), 
    npar(npars), parid(id), f(func),errf(errfunc),varcoll(npars)
{
  temp = *this;
}

Jet::Jet(double Et, double EmEt, double HadEt ,double OutEt, double E,
         double eta,double phi, Flavor flavor, double genPt, double ZSPcor, 
	 double JPTcor, double L2cor, double L3cor,
         double (*func)(const TMeasurement *x, const double *par),
         double (*errfunc)(const double *x, const TMeasurement *xorig, double err),
         double* firstpar, int id, int npars)
  : TJet(Et,EmEt,HadEt,OutEt,E,eta,phi,flavor,genPt,ZSPcor,JPTcor,L2cor,L3cor), 
    par(firstpar), npar(npars), parid(id), f(func),errf(errfunc),varcoll(npars)
{
  temp = *this;
}
 
// varies all parameters for this jet by eps and returns a vector of the
// parameter id and the Et for the par + eps and par - eps variation
const Jet::VariationColl& Jet::varyPars(double eps, double Et, double start)
{
  //start = Et;
  for(int i = 0 ; i < npar ; ++i) {
    double orig = par[i];
    par[i] += eps;
    varcoll[i].upperEt = expectedEt(Et,start);
    if( varcoll[i].upperEt < 0) varcoll[i].upperEt = 0.999 * pt;
    varcoll[i].upperError = expectedError(varcoll[i].upperEt);
    //varcoll[i].upperEt = expectedEt(Et,s,false);
    par[i] = orig - eps;;
    varcoll[i].lowerEt = expectedEt(Et,start); 
    if( varcoll[i].lowerEt < 0) varcoll[i].lowerEt = 1.001 * pt;
    varcoll[i].lowerError = expectedError(varcoll[i].lowerEt);
    //varcoll[i].lowerEt = expectedEt(Et,s,false);
    par[i] = orig;
    varcoll[i].parid = parid + i;
  }
  return varcoll;
}

// varies all parameters for this jet by eps and returns a vector of the
// parameter id and the Et for the par + eps and par - eps variation
const Jet::VariationColl& Jet::varyParsDirectly(double eps)
{
  for(int i = 0 ; i < npar ; ++i) {
    double orig = par[i];
    par[i] += eps;
    varcoll[i].upperEt = correctedEt(pt);
    varcoll[i].upperError = expectedError(varcoll[i].upperEt);
    par[i] = orig - eps;;
    varcoll[i].lowerEt = correctedEt(pt); 
    varcoll[i].lowerError = expectedError(varcoll[i].lowerEt);
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
  double corEt = f(&temp,par);
  if(corEt <  OutF + EMF) corEt = OutF + EMF;
  return corEt;
}

double Jet::expectedEt(double truth, double start, bool fast)
{
  static const double eps = 1.0e-8;
  //const double up = 4 * truth;
  //const double low = 0.2 * truth;
  double x1 = start,x2;
  //find root of truth - jet->correctedEt(expectedEt)
  // x: expectedEt
  // y: truth -  jet->correctedEt(expectedEt)
  // f: jet->correctedEt(expectedEt)
  double f1 = correctedEt(x1,fast);
  //get second point assuming a constant correction factor
  x2 = (truth - EMF - OutF) * (x1 - EMF - OutF)/(f1 - EMF - OutF) + EMF;
  //if((x2 > up )||(x2 < low)) x2 = x1;
  ///double f2 = correctedEt(x2,true);
  //double y2 = truth - f2;
  //std::cout << "truth:" << truth << " scale:" << scale << "  f1:" << f1 << " x2:" << x2 
  //	    << " f2:" << f2 << '\n';
  /*
  if(extrapolate || (std::abs(y2) < eps)) {
    //std::cout << "extrapolated:" << x2 << ", " << y2 << " at scale " << scale << std::endl;  
    return x2;
  }
  */
  if(! secant(truth,x2,x1,eps)) return -1;
  //f1 = correctedEt(scale,true);
  //x2 = (truth - EMF - OutF) * (scale - EMF - OutF)/(f1 - EMF - OutF) + EMF + OutF;
  //std::cout << i << ": scale:" << scale << ", expected:" << (truth - EMF) * (scale - EMF)/(f1 - EMF) + EMF << "  dist for scale:" << truth - f1 << "\n";
  //std::cout << "x1=" << x1 << "  x2=" << x2 << '\n';
  assert(std::abs(correctedEt(x2)-truth)/truth < eps); 
  return x2;
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


bool Jet::secant(double truth, double& x2, double& x1,double eps)
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
  if(dx < 1e-12) {
    x2 = 1.0001 * x1;
    dx = std::abs(x1-x2);
  }
  //std::cout << "first intervall size:" << dx/x1 << '\n';
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
  if(std::abs(y2) > eps * truth) {
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
    std::cout << "Inversion statistics for Jet::expectedEt:\n";
    std::cout << "calls: " << ncalls << " average number of iterations:"
	      << (double)ntries/ncalls << " failures:" << (double)nfails/ncalls*100
	      << "%    warnings:" << (double)nwarns/ntries*100 << "%" <<std::endl;
  }
}
