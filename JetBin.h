
//!    \brief Class for  jet bins 
//!
//!    \author Hartmut Stadie
//!
//!    \date 2010/05/10
//!
//!    $Id: JetBin.h,v 1.5 2010/11/01 15:47:40 stadie Exp $
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
	 float (*errfunc)(const float *x, const Measurement *xorig, float err), 
	 const Function& gf) 
    : sumMess_(),sumPt2_(0), sumGenPt_(0), sumGenPt2_(0), sumdR_(0), sumL1_(0), sumL2_(0), sumL3_(0),
    sumLres_(0),sumL4_(0), sumL5_(0), sumJPT_(0), sumJPTL2L3_(0), njets_(0), f_(&f), gf_(&gf),
    errf_(errfunc) 
  {}
    
  virtual ~JetBin() {}
  
  void addJet(float Et, float EmEt, float HadEt ,float OutEt, float E,
	      float eta,float phi, float phiphi, float etaeta,  
	      float genPt, float dR, const CorFactors& corFactors);
  
  Jet* jet() const { return createJet(); }
  float genPt() const { return njets_ ? sumGenPt_/njets_ : 0;}
  int nJets() const { return njets_;}

 protected:
  Jet* createJet() const;
  Measurement  sumMess_;//!< sums for Measurement
  float sumPt2_;
  float sumGenPt_,sumGenPt2_, sumdR_;        //!< sums for Jet
  float sumL1_, sumL2_, sumL3_, sumLres_, sumL4_, sumL5_, sumJPT_, sumJPTL2L3_;  //!< sums for CorFactor
  int njets_;                      //!< number of jets in this bin
  const Function*  f_;                    //!< Jet correction function
  const Function*  gf_;                   //!< Global jet correction function
  float    (*errf_)(const float *x, const Measurement *xorig, float err);   //!< Error function
};

#endif
