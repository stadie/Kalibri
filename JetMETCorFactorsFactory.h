//!  \brief   Container class for jet correction factors
//
//    $Id: JetMETCorFactorsFactory.h,v 1.1 2009/11/26 10:28:48 stadie Exp $
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

  static JetMETCorFactorsFactory *Summer09_7TeV_AK5Calo,*Summer09_AK5Calo,
    *Summer09_7TeV_ReReco332_AK5Calo;
  
};



#endif
