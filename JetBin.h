//!    \brief Class for  jet bins 
//!
//!    \author Hartmut Stadie
//!
//!    \date 2010/05/10
//!
//!    $Id: JetBin.h,v 1.2 2010/05/19 16:01:42 stadie Exp $
#ifndef JETBIN_H
#define JETBIN_H

#include "Function.h"
#include "CalibData.h"

class CorFactors;
class Jet;

class JetBin 
{
 public:
  JetBin(const Function& f,
	 double (*errfunc)(const double *x, const Measurement *xorig, double err), 
	 const Function& gf) 
    : sumMess_(),sumPt2_(0), sumGenPt_(0), sumGenPt2_(0), sumdR_(0), sumL1_(0), sumL2_(0), sumL3_(0),
    sumL4_(0), sumL5_(0), sumJPT_(0), sumJPTL2L3_(0), njets_(0), f_(f), gf_(gf), errf_(errfunc) 
  {}
    
  virtual ~JetBin() {}
  
  void addJet(double Et, double EmEt, double HadEt ,double OutEt, double E,
	      double eta,double phi, double phiphi, double etaeta,  
	      double genPt, double dR, const CorFactors& corFactors);
  
  Jet* jet() const { return createJet(); }
  double genPt() const { return njets_ ? sumGenPt_/njets_ : 0;}
  int nJets() const { return njets_;}

 protected:
  Jet* createJet() const;

  Measurement  sumMess_;//!< sums for Measurement
  double sumPt2_;
  double sumGenPt_,sumGenPt2_, sumdR_;        //!< sums for Jet
  double sumL1_, sumL2_, sumL3_, sumL4_, sumL5_, sumJPT_, sumJPTL2L3_;  //!< sums for CorFactor
  int njets_;                      //!< number of jets in this bin
  Function  f_;                    //!< Jet correction function
  Function  gf_;                   //!< Global jet correction function
  double    (*errf_)(const double *x, const Measurement *xorig, double err);   //!< Error function
};

#endif
