//
//    $Id: JetMETCorFactorsFactory.cc,v 1.2 2009/11/26 12:41:02 stadie Exp $
//   
#include "JetMETCorFactorsFactory.h"
#include "CorFactors.h"
#include "Jet.h"

#include "JetMETObjects/interface/FactorizedJetCorrector.h"

#include <iostream>

JetMETCorFactorsFactory::JetMETCorFactorsFactory(const std::string& name,
						 const std::string& files)
  : CorFactorsFactory(name)
{
  cor_ = new FactorizedJetCorrector("L2:L3",files);
  std::cout << "created JetMETCorFactorsFactory: " << name << '\n';
}

JetMETCorFactorsFactory::~JetMETCorFactorsFactory()
{
  delete cor_;
}


CorFactors* JetMETCorFactorsFactory::create(const Jet* j)
{
  cor_->setJetEta(j->eta());
  cor_->setJetPt(j->pt()); 
  cor_->setJetE(j->E());
  cor_->setJetPhi(j->phi());
  cor_->setJetEMF(j->EmEt()/j->pt()); 
  
  std::vector<float> levels = cor_->getSubCorrections();

  return new CorFactors(1.0,
			levels[0],
			levels[1]/levels[0],
			1.0,1.0,0.0,0.0);			
}
JetMETCorFactorsFactory::Register JetMETCorFactorsFactory::register_;

JetMETCorFactorsFactory::Register::Register() 
{
  try {
    new JetMETCorFactorsFactory("Summer09_7TeV_AK5Calo","JetMETObjects/data/Summer09_7TeV_L2Relative_AK5Calo.txt:JetMETObjects/data/Summer09_7TeV_L3Absolute_AK5Calo.txt");
  } 
  catch(std::exception) {
    std::cout << "... failed to create!\n";
  } 
  try {
    new JetMETCorFactorsFactory("Summer09_AK5Calo","JetMETObjects/data/Summer09_L2Relative_AK5Calo.txt:JetMETObjects/data/Summer09_L3Absolute_AK5Calo.txt");
  } 
  catch(std::exception) {
    std::cout << "... failed to create!\n";
  } 
  try {
    new JetMETCorFactorsFactory("Summer09_7TeV_ReReco332_AK5Calo",
				"Summer09_7TeV_ReReco332_L2Relative_AK5Calo.txt:Summer09_7TeV_ReReco332_L3Absolute_AK5Calo.txt");
  } 
  catch(std::exception) {
    std::cout << "... failed to create!\n";
  } 
}
