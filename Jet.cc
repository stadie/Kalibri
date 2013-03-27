//
//    Class for basic jets 
//
//    first version: Hartmut Stadie 2008/12/14
//    $Id: Jet.cc,v 1.56 2013/03/27 12:14:13 kirschen Exp $
//   
#include "Jet.h"  

#include "CorFactors.h"
#include "Parameters.h"

#include <iostream>
#include <iomanip>

Jet::Jet(float Et, float EmEt, float HadEt ,float OutEt, float E,
         float eta,float phi, float phiphi, float etaeta, Flavor flavor, 
	 float fCH, float fNH, float fPH, float fEL, float fHFEm, float fHFHad, 
	 float genPt, float dR, CorFactors* corFactors, const Function& f, 
	 float (*errfunc)(const float *x, const Measurement *xorig, float err), 
	 const Function& gf, float closestJetdR) 
  : Measurement(Et,EmEt,HadEt,OutEt,E,eta,phi,phiphi,etaeta),flavor_(flavor), 
    fCH_(fCH), fNH_(fNH), fPH_(fPH), fEL_(fEL), fHFEm_(fHFEm), fHFHad_(fHFHad),
    genPt_(genPt),dR_(dR),closestJetdR_(closestJetdR),corFactors_(corFactors),f_(&f),gf_(&gf),
    errf_(errfunc),parameters_(0)
{
  //std::cout << "size:" << sizeof(Jet::GslImplementation) << ", " << sizeof(Jet) << ", " << sizeof(Function) << ", " << sizeof(Measurement) << '\n';
}

Jet::Jet(const Jet& j) 
  : Measurement(j), flavor_(j.flavor_), 
    fCH_(j.fCH_), fNH_(j.fNH_), fPH_(j.fPH_), fEL_(j.fEL_), fHFEm_(j.fHFEm_), fHFHad_(j.fHFHad_),
    genPt_(j.genPt_),dR_(j.dR_),closestJetdR_(j.closestJetdR_), 
    corFactors_(new CorFactors(*(j.corFactors_))),f_(j.f_),
    gf_(j.gf_),errf_(j.errf_),parameters_(0)
{
}

Jet::~Jet()
{
  delete corFactors_;
}
 
void Jet::setParameters(Parameters* param) {
  parameters_ = param;
  f_ = &(param->function(*f_));
  gf_ = &(param->function(*gf_));
}

void Jet::updateCorFactors(CorFactors *cor)
{
  delete corFactors_;
  corFactors_ = cor;
}

void Jet::correctL1()
{
  Measurement::pt   *= corFactors_->getL1();
  Measurement::E    *= corFactors_->getL1();
  Measurement::EMF  *= corFactors_->getL1();
  Measurement::HadF *= corFactors_->getL1();
  Measurement::OutF *= corFactors_->getL1();
}

void Jet::correctToL3()
{
  Measurement::pt   *= corFactors_->getToL3();
  Measurement::E    *= corFactors_->getToL3();
  Measurement::EMF  *= corFactors_->getToL3();
  Measurement::HadF *= corFactors_->getToL3();
  Measurement::OutF *= corFactors_->getToL3();
}

void Jet::correctL2L3()
{
  Measurement::pt   *= corFactors_->getL2L3();
  Measurement::E    *= corFactors_->getL2L3();
  Measurement::EMF  *= corFactors_->getL2L3();
  Measurement::HadF *= corFactors_->getL2L3();
  Measurement::OutF *= corFactors_->getL2L3();
}

void Jet::correctToLRes() {
  Measurement::pt   *= corFactors_->getToLRes();
  Measurement::E    *= corFactors_->getToLRes();
  Measurement::EMF  *= corFactors_->getToLRes();
  Measurement::HadF *= corFactors_->getToLRes();
  Measurement::OutF *= corFactors_->getToLRes();
}


//!  \brief Varies all parameters for this jet by +/-eps
//!
//!  The corrected Et and errors obtained by applying the
//!  correction function with the varied parameters are
//!  stored in a VariationColl.
//!
//!  \param eps Amount by which the parameters are varied
//!  \param Et Et which is to be corrected after parameter variation
//!  \param start Start value for expectedEt(float truth, float start, bool fast)
//!  \return Vector of ParameterVariation
//!  \sa varyParsDirectly
// -------------------------------------------------------
const Parameters::VariationColl& Jet::varyPars(const double* eps, float Et, float start)
{
  //start = Et; 
  Parameters::VariationColl& varcoll = parameters_->cachedVariationColl();
  varcoll.resize(nPar());
  for(int i = 0 ; i < f_->nPars() ; ++i) {
    double orig = f_->firstPar()[i];
    f_->firstPar()[i] += eps[f_->parIndex() + i];
    varcoll[i].upperEt = expectedEt(Et,start,varcoll[i].upperError);
    //if( varcoll[i].upperEt < 0) return varcoll;
    //varcoll[i].upperEt = expectedEt(Et,s,false);
    f_->firstPar()[i] = orig - eps[f_->parIndex() + i];
    varcoll[i].lowerEt = expectedEt(Et,start,varcoll[i].lowerError); 
    //if( varcoll[i].lowerEt < 0) return varcoll;
    //varcoll[i].lowerEt = expectedEt(Et,s,false); 
    f_->firstPar()[i] =  orig + 2 * eps[f_->parIndex() + i];
    varcoll[i].upperEt2 = expectedEt(Et,start,varcoll[i].upperError2);
    //if( varcoll[i].upperEt < 0) return varcoll;
    //varcoll[i].upperEt = expectedEt(Et,s,false);
    f_->firstPar()[i] = orig - 2 * eps[f_->parIndex() + i];
    varcoll[i].lowerEt2 = expectedEt(Et,start,varcoll[i].lowerError2); 
    //if( varcoll[i].lowerEt < 0) return varcoll;
    //varcoll[i].lowerEt = expectedEt(Et,s,false);
    f_->firstPar()[i] = orig;
    varcoll[i].parid = f_->parIndex() + i;
  }
  for(int i = 0,j =  f_->nPars(); i < gf_->nPars() ; ++i,++j) {
    double orig = gf_->firstPar()[i];
    gf_->firstPar()[i] += eps[gf_->parIndex() + i];
    varcoll[j].upperEt = expectedEt(Et,start,varcoll[j].upperError);
    //if( varcoll[j].upperEt < 0) return varcoll;
    //varcoll[j].upperEt = expectedEt(Et,s,false);
    gf_->firstPar()[i] = orig - eps[gf_->parIndex() + i];
    varcoll[j].lowerEt = expectedEt(Et,start,varcoll[j].lowerError);
    //if( varcoll[j].lowerEt < 0) return varcoll;
    //varcoll[j].lowerEt = expectedEt(Et,s,false);   
    gf_->firstPar()[i] = orig + 2*eps[gf_->parIndex() + i];
    varcoll[j].upperEt2 = expectedEt(Et,start,varcoll[j].upperError2);
    //if( varcoll[j].upperEt < 0) return varcoll;
    //varcoll[j].upperEt = expectedEt(Et,s,false);
    gf_->firstPar()[i] = orig - 2*eps[gf_->parIndex() + i];
    varcoll[j].lowerEt2 = expectedEt(Et,start,varcoll[j].lowerError2);
    //if( varcoll[j].lowerEt < 0) return varcoll;
    //varcoll[j].lowerEt = expectedEt(Et,s,false);
    gf_->firstPar()[i] = orig;
    varcoll[j].parid = gf_->parIndex() + i;
  }
  return varcoll;
}

//!  \brief Varies all parameters for this jet by +/-eps
//!
//!  The corrected original Et and errors obtained by applying the
//!  correction function with the varied parameters are
//!  stored in a VariationColl.
//!
//!  \note In contrast to varyPars, the originally measured
//!        Et is corrected.
//!  \param eps Amount by which the parameters are varied
//!  \return Vector of ParameterVariation
//!  \sa varyPars
// -------------------------------------------------------
const Parameters::VariationColl& Jet::varyParsDirectly(const double* eps, bool computeDeriv, float Et)
{ 
  Parameters::VariationColl& varcoll = parameters_->cachedVariationColl();
  varcoll.resize(nPar());
  if(Et == 0) Et = Measurement::pt;
  const float deltaE = 1e-05 * Et;
  for(int i = 0 ; i < f_->nPars() ; ++i) {
    double orig = f_->firstPar()[i];
    f_->firstPar()[i] += eps[f_->parIndex() + i];
    varcoll[i].upperEt = correctedEt(Et);
    varcoll[i].upperError = expectedError(varcoll[i].upperEt);
    if(computeDeriv) {
      varcoll[i].upperEtDeriv =  (correctedEt(Et+deltaE) -  correctedEt(Et-deltaE))/2/deltaE;
    }
    f_->firstPar()[i] = orig - eps[f_->parIndex() + i];
    varcoll[i].lowerEt = correctedEt(Et); 
    varcoll[i].lowerError = expectedError(varcoll[i].lowerEt);
    if(computeDeriv) {
      varcoll[i].lowerEtDeriv =  (correctedEt(Et+deltaE) -  correctedEt(Et-deltaE))/2/deltaE;
    }
    f_->firstPar()[i] = orig + 2*eps[f_->parIndex() + i];
    varcoll[i].upperEt2 = correctedEt(Et);
    varcoll[i].upperError2 = expectedError(varcoll[i].upperEt2);
    if(computeDeriv) {
      varcoll[i].upperEtDeriv2 =  (correctedEt(Et+deltaE) -  correctedEt(Et-deltaE))/2/deltaE;
    }
    f_->firstPar()[i] = orig - 2*eps[f_->parIndex() + i];
    varcoll[i].lowerEt2 = correctedEt(Et); 
    varcoll[i].lowerError2 = expectedError(varcoll[i].lowerEt2);
    if(computeDeriv) {
      varcoll[i].lowerEtDeriv2 =  (correctedEt(Et+deltaE) -  correctedEt(Et-deltaE))/2/deltaE;
    }
    f_->firstPar()[i] = orig;
    varcoll[i].parid = f_->parIndex() + i;
  }  
  for(int i = 0, j =  f_->nPars(); i < gf_->nPars() ; ++i,++j) {
    double orig = gf_->firstPar()[i];
    gf_->firstPar()[i] += eps[gf_->parIndex() + i];
    varcoll[j].upperEt = correctedEt(Et);
    varcoll[j].upperError = expectedError(varcoll[j].upperEt);
    if(computeDeriv) {
      varcoll[j].upperEtDeriv =  (correctedEt(Et+deltaE) -  correctedEt(Et-deltaE))/2/deltaE;
    }
    gf_->firstPar()[i] = orig - eps[gf_->parIndex() + i];
    varcoll[j].lowerEt = correctedEt(Et); 
    varcoll[j].lowerError = expectedError(varcoll[j].lowerEt);
    if(computeDeriv) {
      varcoll[j].lowerEtDeriv =  (correctedEt(Et+deltaE) -  correctedEt(Et-deltaE))/2/deltaE;
    }   
    gf_->firstPar()[i] = orig + 2*eps[gf_->parIndex() + i];
    varcoll[j].upperEt2 = correctedEt(Et);
    varcoll[j].upperError2 = expectedError(varcoll[j].upperEt2);
    if(computeDeriv) {
      varcoll[j].upperEtDeriv2 =  (correctedEt(Et+deltaE) -  correctedEt(Et-deltaE))/2/deltaE;
    }
    gf_->firstPar()[i] = orig - 2*eps[gf_->parIndex() + i];
    varcoll[j].lowerEt2 = correctedEt(Et); 
    varcoll[j].lowerError2 = expectedError(varcoll[j].lowerEt2);
    if(computeDeriv) {
      varcoll[j].lowerEtDeriv2 =  (correctedEt(Et+deltaE) -  correctedEt(Et-deltaE))/2/deltaE;
    }
    gf_->firstPar()[i] = orig;
    varcoll[j].parid = gf_->parIndex() + i;
  }
  return varcoll;
}


//!  \brief Correct a given jet Et
//!
//!  The given Et is corrected applying successively
//!  the local and the global jet correction functions
//!  f and gf (see also Parametrization).
//!
//!  \note Modifies only hadronic part of tower Et
//!  \param Et Jet Et which gets corrected
//!  \param fast No functionality yet
//!  \return Corrected jet Et
// -------------------------------------------------------
float Jet::correctedEt(float Et, bool fast) const {
  
  //std::cout << "Pars:" << f_->firstPar()[0] << ", " << f_->firstPar()[1] << ", " << f_->firstPar()[2]
  //	    << ", " << gf_->firstPar()[0] << ", " << gf_->firstPar()[1] << ", " << gf_->firstPar()[2]
  //	    << ", " <<  gf_->firstPar()[3] << '\n';
  // 
  //assume that only the hadronic energy gets modified!
  if(Et < 1.0) Et = 1.0;
  float oldpt = Measurement::pt;
  float oldE  = Measurement::E;
  Jet* self = const_cast<Jet*>(this);
  self->Measurement::pt   = Et;  
  //temp_.HadF = Et - OutF - EMF;
  // if(temp_.HadF < 0) temp_.HadF = 0;
  self->Measurement::E    = Et/oldpt;
  self->Measurement::pt = (*f_)(this);
  /*
    if(corEt != corEt) 
    std::cout << "Et:" << Et << "  orig Et:" << pt << " cor Et:" << corEt << "\n";
    assert(corEt == corEt);
    //if(corEt <  OutF + EMF) corEt = OutF + EMF;
  */
  if(Measurement::pt <= 0.1) {
    //std::cout << "WARNING: jet cor. Et <= 0.1 GeV:" << temp_.pt << " at eta " << Measurement::eta << '\n';
    self->Measurement::pt = 0.1;
  }
  //temp_.HadF = corEt - OutF - EMF;
  //if(temp_.HadF < 0) temp_.HadF = 0;
  self->Measurement::E = oldE * Measurement::pt/oldpt;  
  float ptcor = (*gf_)(this);
  //if(corEt != corEt) std::cout << "Et:" << Et << "  orig Et:" << pt << " cor Et:" << corEt << "\n";
  //assert(corEt == corEt);
  //if(corEt <  OutF + EMF) corEt = OutF + EMF;
  if(ptcor <= 1.0) {
    //std::cout << "WARNING: global jet cor. Et <= 0.1 GeV:" << temp_.pt << " at eta " << Measurement::eta << '\n';
    ptcor = 1.0;
  }
  self->Measurement::pt  = oldpt;
  self->Measurement::E   = oldE;
  return ptcor;
}


//!  \brief Find mean measured Et from correction function for a given truth
//!
//!  Finds the mean measured Et \f$ \bar{E}_{T} \f$ corresponding to a
//!  given truth from the correction function \f$ g_{p} \f$ by solving
//!  \f[
//!    0 = E^{\textrm{true}}_{T} - g_{p}(\bar{E}_{T})
//!  \f]
//!  Here, \f$ p \f$ is the set of current parameters.
//!  The mean Et is defined by the (unknown) response \f$ R \f$ by
//!  \f[
//!   \bar{E}_{T} = R(E^{\textrm{true}}_{T}) \cdot E^{\textrm{true}}_{T}.
//!  \f]
//!  The solution is found numerically using 
//!  secant(double truth, double& x1, double& x2, double eps) or
//!  falseposition(double truth, double& x1, double& x2, double eps)
//!  
//!  \param truth The true Et
//!  \param start Start value for inversion procedure
//!  \param fast No functionality yet
//!  \return Mean measured Et for given truth
// -------------------------------------------------------
float Jet::expectedEt(float truth, float start, bool fast)
{
  if(f_->hasInverse()) {
    float oldpt = Measurement::pt;
    float oldE  = Measurement::E;
    Measurement::pt   = truth;  
    Measurement::E    *= truth/Measurement::pt;
    float pt = f_->inverse(this);
    Measurement::pt = oldpt;
    Measurement::E  = oldE;
    return pt;
  }
  static const double eps = 1.0e-10;

  /*
  double x1 = root_;
  //find root of truth - jet->correctedEt(expectedEt)
  // x: expectedEt
  // y: truth -  jet->correctedEt(expectedEt)
  // f: jet->correctedEt(expectedEt)
  double f1 = correctedEt(x1,false);

  
  if(std::abs(f1 - truth) < eps * truth) {
    //std::cout << "this really happens!?\n";
    return x1;
  }
  */
  //bracket the root
  /*
    double x2 = 1.01 * x1;
    
    double f2 = correctedEt(x2,false);
    for( int i = 0 ; ; ++i) {
    if(i > 100) {
    ++nwarns;
    //std::cout << "Warning failed to bag: " << x1 << ", " << x2 << ":" << f1 << " < " << truth << " < " << f2 << std::endl;
    //assert(i < 10);
    x1 = 0.1 * (truth - 20);
    if(x1 < 0.1) x1 = 0.1;
    x2 = 10 * truth;
    break;
    //++nfails_;
    //return -1;
    }
    if(f1 >= f2) {
    if((f1 > truth) && (f2 < truth)) break; 
    double step = 0.5 * i;
    x1 -= 2 * step;
    x2 += 2 * step;
    f1 = correctedEt(x1,false);
    f2 = correctedEt(x2,false);
    } else {
    if(f1 > truth) {
    double step = (f1 - truth) * (x2-x1)/(f2-f1);
    if(step < 0.1) step = 0.1;
    x1 -=  2 * step;
    f1 = correctedEt(x1,false);
    ++ntries_;
    } else if(f2 < truth) {
    double step = (truth - f2) * (x2-x1)/(f2-f1);
    if(step < 0.1) step = 0.1;
    x2 += 2 * step;
    f2 = correctedEt(x2,false);
    ++ntries_;
    } else break;
    }
    }
  */
  //x1 = truth * Measurement::pt / correctedEt(Measurement::pt,false);
  
  double x1 = Measurement::pt,f1;
  double f2,x2 = x1 + error();
  //bracket root from numerical recipes 3rd edition, p447)
  const int ntries = 50;
  const double factor = 1.6;
  f1 = correctedEt(x1,false);
  f2 = correctedEt(x2,false);
  for(int i = 0 ; i < ntries ; ++i) {
    if(((f1 -truth) * (f2 - truth)) <= 0) {
      rf_par par(truth,this);
      if(! parameters_->findRoot(rf,&par,x1,x2,eps)) return -1;
      return x2;
    }
    ++ntries_;
    if(std::abs(f1 - truth) < std::abs(f2 - truth)) {
      f1 = correctedEt(x1 += factor*(x1-x2));
    } else {
      f2 = correctedEt(x2 += factor*(x2-x1));
    }
  }
  ++nfails_;
  return -1;

  //get second point assuming a constant correction factor
  //x2 = (truth - EMF - OutF) * (x1 - EMF - OutF)/(f1 - EMF - OutF) + EMF;
  //x2 = truth * x1 / f1;
  //if(! secant(truth,x2,x1,eps)) return -1;
  //assert(std::abs(correctedEt(x2)-truth) < eps * truth); 
  //root = x2;
  //return x2;
}


//!  \brief Find mean measured Et and error from correction function for a given truth
//!
//!  Finds the mean measured Et and the corresponding error
//!  corresponding to a given truth from the correction function
//!  (see expectedEt(double truth, double start, bool fast)).
//!  The error is calculated from the obtained mean Et using
//!  expectedError(double et) and returned by reference.
//!
//!  Both Et and the error are corrected for effects due to
//!  a cut on the Et spectrum (etmin).
//!  
//!  \param truth The true Et
//!  \param start Start value for inversion procedure
//!  \param error Will be filled with the error
//!  \param fast No functionality yet
//!  \return Mean measured Et for given truth
// -------------------------------------------------------
float Jet::expectedEt(float truth, float start, float& error,bool fast)
{
  //truncate mean for jet min Et-cut
  float m = expectedEt(truth,start,fast);
  if(m < 0) return m;
  float s = expectedError(m);

  // hack
  error = s;
  return m;

//   float x = (etmin - m)/s;
//   if(x < -10) {
//     error = s;
//     return m;
//   }
//   //truncated mean:
//   //m + (E^(-((a - m)^2/(2 s^2))) Sqrt[2/\[Pi]] s)/Erfc[(a - m)/(Sq[2] s)]
//   //truncated RMS
//   //m^2 + s^2 + (E^(-((a - m)^2/(2 s^2))) (a + m) Sqrt[2/\[Pi]] s)/Erfc[(a - m)/(Sqrt[2] s)]
//   float l = exp(-x*x/2) * sqrt(2/M_PI) * s/TMath::Erfc(x/sqrt(2));
//   m += l;
//   s =  expectedError(m);
//   error = sqrt(l*(etmin - m) + s * s);
//   return m;
}


bool Jet::falseposition(double truth, double& x1, double& x2,double eps)
{
  //x2 is the best estimate!
  double f1 = correctedEt(x1,true);
  double f2 = correctedEt(x2,true);
  double temp;
  double step = 0.1 * truth;
  ++ncalls_;
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
  while( 1 ) {
    //std::cout << "x1,2:" << x1 << "," << x2 << " f1,2:" << f1 << ", " << f2 
    //      << " y1,2:" << y1 << ", " << y2 << std::endl;
    if(f1 > truth) {
      x1 -= step;
      f1 = correctedEt(x1,true);
      y1 = truth - f1;
      ++ntries_;
    } else if(f2 < truth) {
      x2 += step;
      f2 = correctedEt(x2,true);
      y2 = truth - f2;
      ++ntries_;
    }
    else break;
    ++i;
    //if(i > 5) break;
    assert(i < 100);
  } 
  i = 0;
  while(std::abs(y2)/truth > eps) {
    //std::cout << i << ":" << x1 << ", " << x2 << " : " << y1 << ", " << y2 << std::endl;
    double x3 = x1 + y1 * (x2-x1)/(f2 - f1);
    double f3 = correctedEt(x3,true);
    double y3 = truth - f3;
    //std::cout << i << ":" << x3 << ":" << y3 << std::endl;
    ++i;
    if(std::abs(y3) / truth < eps) {
      x2 = x3;
      x1 = x1;
      return true;
    }
    if(y3 > 0) {
      x1 = x3;
      f1 = f3;
      y1 = y3;
    } else {
      x2 = x3;
      f2 = f3;
      y2 = y3;
    }
    if(i > 100) {
      ++nfails_;
      ntries_ += i;
      return false;
    }
  } 
  ntries_ += i;
  return true;
}


bool Jet::secant(double truth, double& x2, double& x1,double eps)
{
  //x2 is the best estimate!
  const double up = 5 * x2;
  const double low = 1.0; 
 
  double f1 = correctedEt(x1,true);
  double f2 = correctedEt(x2,true);
  double y2 = truth - f2;
  double y1 = truth - f1;
  ++ncalls_;
  int i = 0;
  double dx = std::abs(x1-x2), dx1 = truth, dx2 = truth;
  if(dx < 1e-12) {
    x2 = 1.0001 * x1;
    dx = std::abs(x1-x2);
  }
  //std::cout << "first intervall size:" << dx/x1 << '\n';
  while((std::abs(y2) > eps * truth)&&(i < 1000)) {
    if(f1 == f2) {
      std::cout << "Warning: no difference in corrected Et:" << f1 << "," << f2 << '\n';
      //print();
      x2 = 0.5 * (x1 + x2);
      ++nfails_;
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
      ++nwarns_;
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
    if(((y3 < 0) && (y1 > 0)) || ((y3 > 0) && (y1 < 0))) {
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
    //       ++nfails_;
    //       ntries_ += i;
    //       return false;
    //     } 
  } 
  ntries_ += i;
  if( i >= 1000) {
    //std::cout << "failed to find good root\n";
    //std::cout << i << ":" << x1 << ", " << x2 << ":" << truth - f2 << " truth:" << truth << "fs:" << f1 << "," << f2 << "\n";
    ++nfails_;
    return false;
  }
  
  if(x2 != x2) {
    std::cout << "failed to find good root\n";
    std::cout << i << ":" << x1 << ", " << x2 << ":" << truth - f2 << "\n";
    ++nfails_;
    return false;
  }
  
  return true;
}



// ------------------------------------------------------------
long long Jet::ncalls_ = 0;
long long Jet::ntries_ = 0;
long long Jet::nfails_ = 0;
long long Jet::nwarns_ = 0;

void Jet::printInversionStats()
{
  if(ncalls_) {
    std::cout << "Inversion statistics for Jet::expectedEt:\n";
    std::cout << "calls: " << ncalls_ << " average number of iterations:"
	      << (double)ntries_/ncalls_ << " failures:" << (double)nfails_/ncalls_*100
	      << "% warnings:" << (double)nwarns_/ntries_*100 << "%" <<std::endl;
  }
}
