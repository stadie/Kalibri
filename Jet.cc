//
//    Class for basic jets 
//
//    first version: Hartmut Stadie 2008/12/14
//    $Id: Jet.cc,v 1.20 2009/04/08 14:46:18 stadie Exp $
//   
#include "Jet.h"  
#include "TMath.h"


Jet::Jet(double Et, double EmEt, double HadEt ,double OutEt, double E,
	 double eta,double phi, Flavor flavor,   
	 const Function& f, 
	 double (*errfunc)(const double *x, const TMeasurement *xorig, double err), 
	 const Function& gf, double Etmin) 
  : TJet(Et,EmEt,HadEt,OutEt,E,eta,phi,flavor,0.0,1.0,1.0,1.0,1.0,1.0,1.0), 
    f(f),gf(gf),errf(errfunc),etmin(Etmin)
{
  temp = *this;
  varcoll.resize(f.nPars() + gf.nPars());
}

Jet::Jet(double Et, double EmEt, double HadEt ,double OutEt, double E,
         double eta,double phi, Flavor flavor, double genPt, double ZSPcor, 
	 double JPTcor, double L2cor, double L3cor,double L2L3cor, 
	 double L2L3JPTcor, const Function& f, 
	 double (*errfunc)(const double *x, const TMeasurement *xorig, double err), 
	 const Function& gf, double Etmin) 
  : TJet(Et,EmEt,HadEt,OutEt,E,eta,phi,flavor,genPt,ZSPcor,JPTcor,L2cor,L3cor,L2L3cor,L2L3JPTcor), 
    f(f),gf(gf),errf(errfunc),etmin(Etmin)
{
  temp = *this;
  varcoll.resize(f.nPars() + gf.nPars());
}
 
// varies all parameters for this jet by eps and returns a vector of the
// parameter id and the Et for the par + eps and par - eps variation
const Jet::VariationColl& Jet::varyPars(double eps, double Et, double start)
{
  //start = Et;
  for(int i = 0 ; i < f.nPars() ; ++i) {
    double orig = f.firstPar()[i];
    f.firstPar()[i] += eps;
    varcoll[i].upperEt = expectedEt(Et,start,varcoll[i].upperError);
    //if( varcoll[i].upperEt < 0) return varcoll;
    //varcoll[i].upperEt = expectedEt(Et,s,false);
    f.firstPar()[i] = orig - eps;;
    varcoll[i].lowerEt = expectedEt(Et,start,varcoll[i].lowerError); 
    //if( varcoll[i].lowerEt < 0) return varcoll;
    //varcoll[i].lowerEt = expectedEt(Et,s,false);
    f.firstPar()[i] = orig;
    varcoll[i].parid = f.parIndex() + i;
  }
  for(int i = 0,j =  f.nPars(); i < gf.nPars() ; ++i,++j) {
    double orig = gf.firstPar()[i];
    gf.firstPar()[i] += eps;
    varcoll[j].upperEt = expectedEt(Et,start,varcoll[j].upperError);
    //if( varcoll[j].upperEt < 0) return varcoll;
    //varcoll[j].upperEt = expectedEt(Et,s,false);
    gf.firstPar()[i] = orig - eps;;
    varcoll[j].lowerEt = expectedEt(Et,start,varcoll[j].lowerError);
    //if( varcoll[j].lowerEt < 0) return varcoll;
    //varcoll[j].lowerEt = expectedEt(Et,s,false);
    gf.firstPar()[i] = orig;
    varcoll[j].parid = gf.parIndex() + i;
  }
  return varcoll;
}

// varies all parameters for this jet by eps and returns a vector of the
// parameter id and the Et for the par + eps and par - eps variation
const Jet::VariationColl& Jet::varyParsDirectly(double eps)
{
  for(int i = 0 ; i < f.nPars() ; ++i) {
    double orig = f.firstPar()[i];
    f.firstPar()[i] += eps;
    varcoll[i].upperEt = correctedEt(pt);
    varcoll[i].upperError = expectedError(varcoll[i].upperEt);
    f.firstPar()[i] = orig - eps;;
    varcoll[i].lowerEt = correctedEt(pt); 
    varcoll[i].lowerError = expectedError(varcoll[i].lowerEt);
    f.firstPar()[i] = orig;
    varcoll[i].parid = f.parIndex() + i;
  }  
  for(int i = 0, j =  f.nPars(); i < gf.nPars() ; ++i,++j) {
    double orig = gf.firstPar()[i];
    gf.firstPar()[i] += eps;
    varcoll[j].upperEt = correctedEt(pt);
    varcoll[j].upperError = expectedError(varcoll[j].upperEt);
    gf.firstPar()[i] = orig - eps;;
    varcoll[j].lowerEt = correctedEt(pt); 
    varcoll[j].lowerError = expectedError(varcoll[j].lowerEt);
    gf.firstPar()[i] = orig;
    varcoll[j].parid = gf.parIndex() + i;
  }
  return varcoll;
}


//!  \note Modifies only hadronic part of tower energy
double Jet::correctedEt(double Et, bool fast) const {
  
  //std::cout << "Pars:" << f.firstPar()[0] << ", " << f.firstPar()[1] << ", " << f.firstPar()[2]
  //	    << ", " << gf.firstPar()[0] << ", " << gf.firstPar()[1] << ", " << gf.firstPar()[2]
  //	    << ", " <<  gf.firstPar()[3] << '\n';
  // 
  //assume that only the hadronic energy gets modified!
  temp.pt   = Et;  
  temp.HadF = Et - OutF - EMF;
  if(temp.HadF < 0) temp.HadF = 0;
  temp.E    = TJet::E * Et/TJet::pt;
  double corEt = f(&temp);
  if(corEt != corEt) 
    std::cout << "Et:" << Et << "  orig Et:" << pt << " cor Et:" << corEt << "\n";
  assert(corEt == corEt);
  //if(corEt <  OutF + EMF) corEt = OutF + EMF;
  if(corEt < 0) {
    //    std::cout << "WARNING: jet cor. Et < 0.0 GeV:" << corEt << '\n';
    corEt = 1.0;
  }
  temp.pt   = corEt;  
  temp.HadF = corEt - OutF - EMF;
  if(temp.HadF < 0) temp.HadF = 0;
  temp.E    = TJet::E * corEt/TJet::pt;
  corEt = gf(&temp);
  if(corEt != corEt) 
    std::cout << "Et:" << Et << "  orig Et:" << pt << " cor Et:" << corEt << "\n";
  assert(corEt == corEt);
  //if(corEt <  OutF + EMF) corEt = OutF + EMF;
  if(corEt <= 1.0) {
    //    std::cout << "WARNING: global jet cor. Et < 1.0 GeV:" << corEt << '\n';
    corEt = 1.0;
  }
  return corEt;
}


//!  \brief Find mean measured Et from correction function for a given truth
//!
//!  
//!  \param truth The true Et
//!  \param start Start value for inversion procedure
//!  \param fast No functionality yet
//!  \return Mean measured Et for given truth
double Jet::expectedEt(double truth, double start, bool fast)
{
  static const double eps = 1.0e-12;
  //const double up = 4 * truth;
  //const double low = 0.2 * truth;
  double x1 = start,x2;
  //find root of truth - jet->correctedEt(expectedEt)
  // x: expectedEt
  // y: truth -  jet->correctedEt(expectedEt)
  // f: jet->correctedEt(expectedEt)
  double f1 = correctedEt(x1,false);
  if(std::abs(f1 - truth) < eps) {
    return x1;
  }
  //get second point assuming a constant correction factor
  //x2 = (truth - EMF - OutF) * (x1 - EMF - OutF)/(f1 - EMF - OutF) + EMF;
  x2 = truth * x1 / f1;
  //if((x2 > up )||(x2 < low)) x2 = x1;
  ///double f2 = correctedEt(x2,true);
  //double y2 = truth - f2;
  //std::cout << "truth:" << truth << " start:" << start << "  f1:" << f1 << " x2:" << x2	<< '\n';
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

double Jet::expectedEt(double truth, double start, double& error,bool fast)
{
  //truncate mean for jet min Et-cut
  double m = expectedEt(truth,start,fast);
  if(m < 0) return m;
  double s = expectedError(m);
  double x = (etmin - m)/s;
  if(x < -10) {
    error = s;
    return m;
  }
  //truncated mean:
  //m + (E^(-((a - m)^2/(2 s^2))) Sqrt[2/\[Pi]] s)/Erfc[(a - m)/(Sq[2] s)]
  //truncated RMS
  //m^2 + s^2 + (E^(-((a - m)^2/(2 s^2))) (a + m) Sqrt[2/\[Pi]] s)/Erfc[(a - m)/(Sqrt[2] s)]
  double l = exp(-x*x/2) * sqrt(2/M_PI) * s/TMath::Erfc(x/sqrt(2));
  m += l;
  s =  expectedError(m);
  error = sqrt(l*(etmin - m) + s * s);
  return m;
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
  const double low = 1; 
  if(x2 > up ) x2 = up;
  if(x2 < low) x2 = low;
 
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
    if(f1 == f2) {
      //std::cout << "Warning: no difference in corrected Et:" << f1 << "," << f2 << '\n';
      //print();
      x2 = 0.5 * (x1 + x2);
      ++nfails;
      return false;
    }
    double x3 = x1 + y1 * (x2-x1)/(f2 - f1);
    //std::cout << i << ":" << x1 << ", " << x2 << " : " << y1 << ", " << y2 << ", " << x3 << std::endl;
    if(x3 < low) x3 = low;
    if(x3 > up) x3 = up;

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
    //std::cout << x1 << ":" << f1 << " " << x2 << ":" << f2 << " " << x3 << ":" << f3 << '\n';
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
    //std::cout << i << ":" << x1 << ", " << x2 << ":" << truth - f2 << " truth:" << truth << "fs:" << f1 << "," << f2 << "\n";
    ++nfails;
    return false;
  }
  
  if(x2 != x2) {
//     std::cout << "failed to find good root\n";
//     std::cout << i << ":" << x1 << ", " << x2 << ":" << truth - f2 << "\n";
    ++nfails;
    return false;
  }
  
  return true;
}



// ------------------------------------------------------------
void Jet::print()
{
  std::cout << "Jet  Et: " << Et() << " GeV, eta: " << eta() << std::endl;
}




// ------------------------------------------------------------
long long Jet::ncalls = 0;
long long Jet::ntries = 0;
long long Jet::nfails = 0;
long long Jet::nwarns = 0;
void Jet::printInversionStats()
{
  if(ncalls) {
    std::cout << "Inversion statistics for Jet::expectedEt:\n";
    std::cout << "calls: " << ncalls << " average number of iterations:"
	      << (double)ntries/ncalls << " failures:" << (double)nfails/ncalls*100
	      << "%    warnings:" << (double)nwarns/ntries*100 << "%" <<std::endl;
  }
}
