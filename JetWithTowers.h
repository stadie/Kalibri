//
//    Class for jets with towers 
//
//    first version: Hartmut Stadie 2008/12/25
//    $Id: JetWithTowers.h,v 1.4 2009/01/13 13:39:24 stadie Exp $
//   
#ifndef JETWITHTOWERS_H
#define JETWITHTOWERS_H

#include"Jet.h"

#include <vector>
#include <map>

class JetWithTowers : public Jet
{
 public:
  JetWithTowers(double Et, double EmEt, double HadEt ,double OutEt, double E,
		double eta,double phi, Flavor flavor,
		double (*func)(const TMeasurement *x, const double *par),
		double (*errfunc)(const double *x, const TMeasurement *xorig, double err),
		double* firstpar, int id, int npars);
  virtual ~JetWithTowers(); 
  virtual int nPar() const {return njetpars + towerpars.size() * ntowerpars;}
  virtual void ChangeParAddress(double* oldpar, double* newpar);
  virtual double correctedEt(double Et,bool fast = false) const; 
  virtual double Error() const;
  virtual double expectedError(double truth) const;
  //varies the i'th parameter for this jet by eps and returns its overall 
  // parameter id and sets the Et for the par + eps and par - eps result
  virtual int varyPar(int i, double eps, double Et, double scale, double& upperEt, double& lowerEt);
  // varies all parameters for this jet by eps and returns a vector of the
  // parameter id and the Et for the par + eps and par - eps variation
  virtual const VariationColl& varyPars(double eps, double Et, double scale);
  void addTower(double Et, double EmEt, double HadEt ,double OutEt, double E,
		double eta,double phi,
		double (*func)(const TMeasurement *x, const double *par),
		double (*errfunc)(const double *x, const TMeasurement *xorig, double err),
		double* firstpar, int id, int npars);
 private:
  class Tower : public TMeasurement {
  public:
    Tower(double Et, double EmEt, double HadEt ,double OutEt, double E,
	  double eta,double phi, double alpha,
	  double (*func)(const TMeasurement *x, const double *par),
	  double (*errfunc)(const double *x, const TMeasurement *xorig, double err), 
	  double* firstpar, int id, int npars);
    double Et()     const {return pt;}
    double EmEt()   const {return EMF;}
    double HadEt()  const {return HadF;}
    double OutEt()  const {return OutF;}
    double E()      const {return TMeasurement::E;}
    double eta()    const {return TMeasurement::eta;}
    double phi()    const {return TMeasurement::phi;}
    double projectionToJetAxis() const {return alpha;}
    double fractionOfJetHadEt() const { return fraction;}
    void setFractionOfJetHadEt(double frac) const { fraction = frac;}
    void ChangeParAddress(double* oldpar, double* newpar) {par += newpar - oldpar;}
    double correctedHadEt(double HadEt) const;
    double lastCorrectedHadEt() const { return lastCorHadEt;}  
    double Error() const {return errf(&(TMeasurement::pt),this,0);}
    double expectedError(double truth) const { return  errf(&truth,this,0);}
    int nPar() const {return npar;}
    int FirstPar() const {return parid;}
    double *Par() const {return par;}
  private:
    double alpha;
    double* par;//address to first parameter for this jet
    int npar,parid;
    double error; 
    mutable TMeasurement temp;
    mutable double lastCorHadEt;
    mutable double fraction;
    double (*f)(const TMeasurement *x, const double *par);
    double (*errf)(const double *x, const TMeasurement *xorig, double err);
  };
  typedef std::vector<Tower*> TowerColl;
  typedef TowerColl::iterator TowerCollIter;
  typedef TowerColl::const_iterator TowerCollConstIter;
  TowerColl towers;
  int njetpars,ntowerpars;
  std::map<int,double*> towerpars;
};

#endif
