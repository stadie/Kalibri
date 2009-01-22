//
//    Class for basic jets 
//
//    first version: Hartmut Stadie 2008/12/14
//    $Id: Jet.h,v 1.6 2009/01/16 08:46:40 stadie Exp $
//   
#ifndef JET_H
#define JET_H

#include"CalibData.h"

class Jet : public TJet
{
 public:
  Jet(double Et, double EmEt, double HadEt ,double OutEt, double E,
      double eta,double phi, Flavor flavor,
      double (*func)(const TMeasurement *x, const double *par),
      double (*errfunc)(const double *x, const TMeasurement *xorig, double err), 
      double* firstpar, int id, int npars);
  virtual ~Jet() {};
  double Et()     const {return pt;}
  double EmEt()   const {return EMF;}
  double HadEt()  const {return HadF;}
  double OutEt()  const {return OutF;}
  double E()      const {return TMeasurement::E;}
  double eta()    const {return TMeasurement::eta;}
  double phi()    const {return TMeasurement::phi;}
  Flavor flavor() const {return TJet::flavor;}
  virtual void ChangeParAddress(double* oldpar, double* newpar) {par += newpar - oldpar;}
  virtual double correctedEt(double Et, bool fast = false) const;
  double expectedEt(double truth, double& scale, bool extrapolate = false);
  virtual double Error() const {return errf(&(TMeasurement::pt),this,0);}
  virtual double expectedError(double truth) const { return  errf(&truth,this,0);}
  virtual int nPar() const {return npar;}
  //varies the i'th parameter for this jet by eps and returns its overall 
  // parameter id and sets the Et for the par + eps and par - eps result
  virtual int varyPar(int i, double eps, double Et, double scale, double& upperEt, double& lowerEt);
  struct ParameterVariation {
    int parid;
    double upperEt;
    double lowerEt;
    double upperError;
    double lowerError;
  };
  typedef std::vector<ParameterVariation> VariationColl;
  typedef std::vector<ParameterVariation>::const_iterator VariationCollIter;
  // varies all parameters for this jet by eps and returns a vector of the
  // parameter id and the Et for the par + eps and par - eps variation
  virtual const VariationColl& varyPars(double eps, double Et, double scale);
  virtual const VariationColl& varyParsDirectly(double eps);

  static void printInversionStats();
 private:
  double* par;//address to first parameter for this jet 
  int npar,parid;
  double error; 
  double (*f)(const TMeasurement *x, const double *par);
  double (*errf)(const double *x, const TMeasurement *xorig, double err);
 protected:
  mutable VariationColl varcoll;
 private:
  mutable TMeasurement temp;
  bool secant(double truth, double& x1, double& x2, double eps);
  bool falseposition(double truth, double& x1, double& x2, double eps);
  static int ncalls, ntries, nfails, nwarns;
};

#endif
