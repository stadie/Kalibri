//
//    $Id: JetMETCorFactorsFactory.cc,v 1.8 2010/05/01 09:35:34 stadie Exp $
//   
#include "JetMETCorFactorsFactory.h"
#include "CorFactors.h"
#include "Jet.h"

#include "JetMETObjects/interface/FactorizedJetCorrector.h"
#include "JetMETObjects/interface/JetCorrectorParameters.h"

#include <iostream>
#include <vector>
#include <string>
#include <sstream>

JetMETCorFactorsFactory::JetMETCorFactorsFactory(const std::string& name,
						 const std::string& files)
  : CorFactorsFactory(name)
{
  //split files by ":"
  std::vector<JetCorrectorParameters> vParam;
  std::stringstream ss(files);
  std::string file;
  while(std::getline(ss,file,':')) {
    vParam.push_back(JetCorrectorParameters(file));
  }
  cor_ = new FactorizedJetCorrector(vParam);
  
  //std::cout << "created JetMETCorFactorsFactory: " << name << '\n';
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
  cor_->setJetEMF(j->EmEt()/(j->EmEt() + j->HadEt())); 
  
  std::vector<float> levels = cor_->getSubCorrections();

  //std::cout << "eta:" << j->eta() << "  cor levels:" << levels.size() << " :";
  //std::cout << levels[0] << ", " << levels[1];
  //if(levels.size() == 3) std::cout << ", "<< levels[2] << '\n';
  //else std::cout << levels[2] << '\n';
  return new CorFactors(1.0,
			levels[0],
			(levels.size() == 3) ? levels[2]/levels[0] : levels[1]/levels[0],
			1.0,1.0,0.0,0.0);			
}
JetMETCorFactorsFactory::Register JetMETCorFactorsFactory::register_;

JetMETCorFactorsFactory::Register::Register() 
{
  //create("Summer09_7TeV_AK5Calo","JetMETObjects/data/Summer09_7TeV_L2Relative_AK5Calo.txt:JetMETObjects/data/Summer09_7TeV_L3Absolute_AK5Calo.txt");
  //create("Summer09_AK5Calo","JetMETObjects/data/Summer09_L2Relative_AK5Calo.txt:JetMETObjects/data/Summer09_L3Absolute_AK5Calo.txt");
  create("Summer09_7TeV_ReReco332_AK5Calo","JetMETObjects/data/Summer09_7TeV_ReReco332_L2Relative_AK5Calo.txt:JetMETObjects/data/Summer09_7TeV_ReReco332_L3Absolute_AK5Calo.txt");  
  create("Spring10_AK5Calo","JetMETObjects/data/Spring10_L2Relative_AK5Calo.txt:JetMETObjects/data/Spring10_L3Absolute_AK5Calo.txt");
  create("Spring10_AK5PF","JetMETObjects/data/Spring10_L2Relative_AK5PF.txt:JetMETObjects/data/Spring10_L3Absolute_AK5PF.txt");
  create("Spring10_AK5TRK","JetMETObjects/data/Spring10_L2Relative_AK5TRK.txt:JetMETObjects/data/Spring10_L3Absolute_AK5TRK.txt");
  create("Spring10_AK5JPT","JetMETObjects/data/Spring10_L2Relative_AK5JPT.txt:JetMETObjects/data/Spring10_L3Absolute_AK5JPT.txt");
  create("Spring10_AK5CaloData","JetMETObjects/data/Spring10_L2Relative_AK5Calo.txt:JetMETObjects/data/Spring10_L3Absolute_AK5Calo.txt:JetMETObjects/data/Spring10DataV1_L2L3Residual_AK5Calo.txt"); 
  create("Spring10_AK5PFData","JetMETObjects/data/Spring10_L2Relative_AK5PF.txt:JetMETObjects/data/Spring10_L3Absolute_AK5PF.txt:JetMETObjects/data/Spring10DataV1_L2L3Residual_AK5PF.txt"); 
  create("Spring10_AK5JPTData","JetMETObjects/data/Spring10_L2Relative_AK5JPT.txt:JetMETObjects/data/Spring10_L3Absolute_AK5JPT.txt:JetMETObjects/data/Spring10DataV1_L2L3Residual_AK5JPT.txt");
}

JetMETCorFactorsFactory* JetMETCorFactorsFactory::Register::create(const std::string& name, const std::string& files) const
{
  JetMETCorFactorsFactory* jmcff = 0;
  try {
    jmcff = new JetMETCorFactorsFactory(name,files);
  } 
  catch(std::exception& e) {
    std::cout << "...failed to create " << name << ":\n";
    std::cout << "     " << e.what() << std::endl;
    jmcff = 0;
  } 
  return jmcff;
}
