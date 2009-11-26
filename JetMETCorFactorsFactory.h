//!  \brief   Container class for jet correction factors
//
//    $Id: CorFactors.h,v 1.1 2009/11/25 13:07:45 stadie Exp $
//   
#ifndef JETMETCORFACTORSFACTORY_H
#define JETMETCORFACTORSFACTORY_H

#include "CorFactorsFactory.h"

class FactorizedJetCorrector;

class JetMETCorFactorsFactory : public CorFactorsFactory
{
 public:
  JetMETCorFactorsFactory(const std::string& name, const std::string& files);
  ~JetMETCorFactorsFactory();

  CorFactors* create(const Jet* j);
 private:
  FactorizedJetCorrector* cor_;

  static JetMETCorFactorsFactory *Summer09_7TeV_AK5Calo,*Summer09_AK5Calo;
  
};



#endif
