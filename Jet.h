//
//    Class for basic jets 
//
//    first version: Hartmut Stadie 2008/12/14
//    $Id: EventReader.h,v 1.1 2008/12/12 13:43:15 stadie Exp $
//   
#ifndef JET_H
#define JET_H

#include"CalibData.h"

class Jet : public TJet
{
 public:
  Jet(double Et, double EmEt, double HadEt ,double OutEt, double E,
      double eta,double phi, Flavor flavor,
      double const(*func)(TMeasurement *const x, double *const par),
      double err, double* firstpar, int id, int npars);
  double Et()     const {return pt;}
  double EmEt()   const {return EMF;}
  double HadEt()  const {return HadF;}
  double OutEt()  const {return OutF;}
  double E()      const {return TMeasurement::E;}
  double eta()    const {return TMeasurement::eta;}
  double phi()    const {return TMeasurement::phi;}
  Flavor flavor() const {return TJet::flavor;}
  void ChangeParAddress(double* oldpar, double* newpar) {par += newpar - oldpar;}
  double correctedEt(double Et);
  double Error() const {return error;}
  int nPar() const {return npar;}
  //varies the i'th parameter for this jet by eps and returns its overall 
  // parameter id and sets the Et for the par + eps and par - eps result
  int varyPar(int i, double eps, double Et, double& upperEt, double& lowerEt);
 protected:
  double* par;//address to first parameter for this jet
  int npar,parid;
  double error; 
  mutable TMeasurement temp;
  double const(*f)(TMeasurement *const x, double *const par);
  //double const(*err)(double *const x, TMeasurement *const x_original, double const error);
};

#endif
