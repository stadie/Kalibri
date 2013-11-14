//!  \brief   Container class for jet correction factors
//
//    $Id: JetMETCorFactorsFactory.h,v 1.8 2011/04/01 10:23:50 kirschen Exp $
//   
#ifndef JETMETCORFACTORSFACTORY_H
#define JETMETCORFACTORSFACTORY_H

#include "CorFactorsFactory.h"
#include <boost/thread/mutex.hpp>
#include <vector>

class FactorizedJetCorrector;
class JetCorrectionUncertainty;
class JetCorrectorParameters;

class JetMETCorFactorsFactory : public CorFactorsFactory
{  
 public:
  enum Levels {L1L2L3,L2L3,L2L3res,L2L3L4,L1L2L3resL4,L1L2L3res,L1L2L3UP,L1L2L3DOWN};
  JetMETCorFactorsFactory(const std::string& name, std::string files,Levels type);
  ~JetMETCorFactorsFactory();

  CorFactors* create(const Jet* j, int nPV = 1, double rho=-9999., double jetA=0.);
  CorFactorsFactory* clone() const {
    return new JetMETCorFactorsFactory(*this);
  }
 private:
  FactorizedJetCorrector* cor_; 
  JetCorrectionUncertainty* corunc_; 
  std::vector<JetCorrectorParameters> vParam_;
  Levels type_;

  class Register {
  public:
    Register();
  private:
    JetMETCorFactorsFactory* create(const std::string& name, const std::string& files,Levels type) const;
  };
  static Register register_;
  static boost::mutex mutex_;
};



#endif
