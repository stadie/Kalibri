//
//    Class for basic jets 
//
//    first version: Hartmut Stadie 2008/12/14
//    $Id: Jet.h,v 1.14 2009/04/08 14:46:18 stadie Exp $
//   
#ifndef JET_H
#define JET_H

#include"CalibData.h"
#include "Function.h"


class Jet : public TJet
{
 public:
  Jet(double Et, double EmEt, double HadEt ,double OutEt, double E,
      double eta,double phi, Flavor flavor, const Function& f, 
      double (*errfunc)(const double *x, const TMeasurement *xorig, double err), 
      const Function& gf, double Etmin = 0);
  Jet(double Et, double EmEt, double HadEt ,double OutEt, double E,
      double eta,double phi, Flavor flavor, double genPt, double ZSPcor,
      double JPTcor, double L2cor, double L3cor,double L2L3cor, 
      double L2L3JPTcor, const Function& f,
      double (*errfunc)(const double *x, const TMeasurement *xorig, double err), 
      const Function& gf, double Etmin = 0); 
  virtual ~Jet() {};
  double Et()     const {return pt;}
  double EmEt()   const {return EMF;}
  double HadEt()  const {return HadF;}
  double OutEt()  const {return OutF;}
  double E()      const {return TMeasurement::E;}
  double eta()    const {return TMeasurement::eta;}
  double phi()    const {return TMeasurement::phi;}
  Flavor flavor() const {return TJet::flavor;}
  virtual void ChangeParAddress(double* oldpar, double* newpar) {
    f.changeParBase(oldpar,newpar);
    gf.changeParBase(oldpar,newpar);
  }
  virtual double correctedEt(double Et, bool fast = false) const;
  double expectedEt(double truth, double start, double& error,
		    bool fast = false);
  virtual double Error() const {return errf(&(TMeasurement::pt),this,0);}  //!< Return error of original measurement
  virtual double expectedError(double truth) const { return errf(&truth,this,0);}
  double errorFunction(double pt) const { return errf(&pt,this,0);}
  virtual int nPar() const {return f.nPars() + gf.nPars();}
  struct ParameterVariation {
    int    parid;
    double upperEt;
    double lowerEt;
    double upperError;
    double lowerError;
    bool operator==(int b) const { return parid == b;}
  };
  typedef std::vector<ParameterVariation> VariationColl;
  typedef std::vector<ParameterVariation>::const_iterator VariationCollIter;
  // varies all parameters for this jet by eps and returns a vector of the
  // parameter id and the Et for the par + eps and par - eps variation
  virtual const VariationColl& varyPars(double eps, double Et, double start);
  virtual const VariationColl& varyParsDirectly(double eps);

  void print();                       //!< Print some jet members
  static void printInversionStats();  //!< Print some info on inversion

 protected:
  mutable VariationColl varcoll;
  virtual double expectedEt(double truth, double start, bool fast = false);

 private:
  double    error;                //!< Stores error for constant error mode
  Function  f;                    //!< Jet correction function
  Function  gf;                   //!< Global jet correction function
  double    (*errf)(const double *x, const TMeasurement *xorig, double err);   //!< Error function
  double    etmin;                //!< Lower cut on measured Et

  bool      secant(double truth, double& x1, double& x2, double eps);
  bool      falseposition(double truth, double& x1, double& x2, double eps);

  static long long ncalls;        //!< Number of calls of inversion methods secant and falseposition
  static long long ntries;        //!< Number of tries in iteration during inversion
  static long long nfails;        //!< Number of failed tries during inversion
  static long long nwarns;        //!< Number of warnings during inversion

  mutable TMeasurement temp;
};

#endif
