//
//    Class for jets with towers 
//
//    first version: Hartmut Stadie 2008/12/25
//    $Id: JetWithTowers.h,v 1.1 2008/12/27 16:39:02 stadie Exp $
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
		double const(*func)(TMeasurement *const x, double *const par),
		double err, double* firstpar, int id, int npars);
  virtual ~JetWithTowers();
  virtual void ChangeParAddress(double* oldpar, double* newpar);
  virtual double correctedEt(double Et) const;
  //varies the i'th parameter for this jet by eps and returns its overall 
  // parameter id and sets the Et for the par + eps and par - eps result
  virtual int varyPar(int i, double eps, double Et, double scale, double& upperEt, double& lowerEt);
  // varies all parameters for this jet by eps and returns a vector of the
  // parameter id and the Et for the par + eps and par - eps variation
  virtual const VariationColl& varyPars(double eps, double Et, double scale);
  void addTower(double Et, double EmEt, double HadEt ,double OutEt, double E,
		double eta,double phi,
		double const(*func)(TMeasurement *const x, double *const par),
		double err, double* firstpar, int id, int npars);
 protected:
  class Tower : public TMeasurement {
  public:
    Tower(double Et, double EmEt, double HadEt ,double OutEt, double E,
	  double eta,double phi, double alpha,
	  double const(*func)(TMeasurement *const x, double *const par),
	  double err, double* firstpar, int id, int npars);
    double Et()     const {return pt;}
    double EmEt()   const {return EMF;}
    double HadEt()  const {return HadF;}
    double OutEt()  const {return OutF;}
    double E()      const {return TMeasurement::E;}
    double eta()    const {return TMeasurement::eta;}
    double phi()    const {return TMeasurement::phi;}
    double projectionToJetAxis() const {return alpha;}
    void ChangeParAddress(double* oldpar, double* newpar) {par += newpar - oldpar;}
    double correctedEt(double Et) const;
    double lastCorrectedEt() const { return lastCorEt;}
    double Error() const {return error;}
    int nPar() const {return npar;}
    int FirstPar() const {return parid;}
    double *Par() const {return par;}
  private:
    double alpha;
    double* par;//address to first parameter for this jet
    int npar,parid;
    double error; 
    mutable TMeasurement temp;
    mutable double lastCorEt;
    double const(*f)(TMeasurement *const x, double *const par);
  };
  typedef std::vector<Tower*> TowerColl;
  typedef TowerColl::iterator TowerCollIter;
  typedef TowerColl::const_iterator TowerCollConstIter;
  TowerColl towers;
  int njetpars,ntowerpars;
  std::map<int,double*> towerpars;
};

#endif
