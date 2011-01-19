//!  \brief   Container class for jet correction factors
//
//    $Id: JetMETCorFactorsFactory.h,v 1.5 2010/10/14 17:26:55 stadie Exp $
//   
#ifndef JETMETCORFACTORSFACTORY_H
#define JETMETCORFACTORSFACTORY_H

#include "CorFactorsFactory.h"
#include <vector>

class FactorizedJetCorrector;
class JetCorrectorParameters;

class JetMETCorFactorsFactory : public CorFactorsFactory
{  
 public:
  enum Levels {L1L2L3,L2L3,L2L3res,L2L3L4};
  JetMETCorFactorsFactory(const std::string& name, const std::string& files,Levels type);
  ~JetMETCorFactorsFactory();

  CorFactors* create(const Jet* j);
  CorFactorsFactory* clone() const {
    return new JetMETCorFactorsFactory(*this);
  }
 private:
  FactorizedJetCorrector* cor_;
  std::vector<JetCorrectorParameters> vParam_;
  Levels type_;

  class Register {
  public:
    Register();
  private:
    JetMETCorFactorsFactory* create(const std::string& name, const std::string& files,Levels type) const;
  };
  static Register register_;
};



#endif
