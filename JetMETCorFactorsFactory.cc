//
//    $Id: JetMETCorFactorsFactory.cc,v 1.17 2011/02/18 15:33:54 stadie Exp $
//   
#include "JetMETCorFactorsFactory.h"
#include "CorFactors.h"
#include "Jet.h"

#include "JetMETObjects/interface/FactorizedJetCorrector.h"
#include "JetMETObjects/interface/JetCorrectorParameters.h"

#include <iostream>
#include <string>
#include <sstream>

boost::mutex JetMETCorFactorsFactory::mutex_;

JetMETCorFactorsFactory::JetMETCorFactorsFactory(const std::string& name,
						 const std::string& files,
						 Levels type)
  : CorFactorsFactory(name), type_(type)
{
  //split files by ":"
  std::stringstream ss(files);
  std::string file;
  while(std::getline(ss,file,':')) {
    vParam_.push_back(JetCorrectorParameters(file));
    //vParam_.back().printScreen();
  }
  cor_ = new FactorizedJetCorrector(vParam_);
  //std::cout << "created JetMETCorFactorsFactory: " << name << '\n';
}

JetMETCorFactorsFactory::~JetMETCorFactorsFactory()
{
  //delete cor_;
}


CorFactors* JetMETCorFactorsFactory::create(const Jet* j,int nPV)
{
  boost::mutex::scoped_lock l(mutex_);
  if(nPV < 1) nPV = 1;
  cor_->setJetEta(j->eta());
  cor_->setJetPt(j->pt()); 
  cor_->setJetE(j->E());
  cor_->setJetPhi(j->phi());
  cor_->setJetEMF(j->EmEt()/(j->EmEt() + j->HadEt())); 
  cor_->setJetEtaEtaMoment(j->momentEtaEta());
  cor_->setJetPhiPhiMoment(j->momentPhiPhi());
  cor_->setNPV(nPV);

  std::vector<float> levels = cor_->getSubCorrections();

  //std::cout << "eta:" << j->eta() << " pt:" << j->pt() << " etaeta:" << j->momentEtaEta() << "  cor levels:" << levels.size() << " :";
  //std::cout << levels[0] << ", " << levels[1];
  //if(levels.size() == 3) std::cout << ", "<< levels[2] << '\n';
  //else std::cout << '\n';
  switch(type_) {
  case L2L3:
    return new CorFactors(1.0,
			  levels[0],
			  levels[1]/levels[0],
			  1.0,1.0,1.0,0.0,0.0);			
  case L1L2L3:
    return new CorFactors(levels[0],
			  levels[1]/levels[0],
			  levels[2]/levels[1],
			  1.0,1.0,1.0,0.0,0.0);
  case L2L3res:
    return new CorFactors(1.0,
			  levels[0],
			  levels[1]/levels[0],
			  levels[2]/levels[1],
			  1.0,1.0,0.0,0.0); 
  case L1L2L3res:
    return new CorFactors(levels[0],
			  levels[1]/levels[0],
			  levels[2]/levels[1],
			  levels[3]/levels[2],
			  1.0,1.0,0.0,0.0);
  case L2L3L4:
    return new CorFactors(1.0,
			  levels[0],
			  levels[1]/levels[0],
			  1.0,
			  levels[2]/levels[1],
			  1.0,0.0,0.0);
  };
  return new CorFactors(1.0,1.0,1.0,1.0,
			1.0,1.0,0.0,0.0);			
}

JetMETCorFactorsFactory::Register JetMETCorFactorsFactory::register_;

JetMETCorFactorsFactory::Register::Register() 
{
  //create("Summer09_7TeV_AK5Calo","JetMETObjects/data/Summer09_7TeV_L2Relative_AK5Calo.txt:JetMETObjects/data/Summer09_7TeV_L3Absolute_AK5Calo.txt");
  //create("Summer09_AK5Calo","JetMETObjects/data/Summer09_L2Relative_AK5Calo.txt:JetMETObjects/data/Summer09_L3Absolute_AK5Calo.txt");
 
  //  create("Summer09_7TeV_ReReco332_AK5Calo","JetMETObjects/data/Summer09_7TeV_ReReco332_L2Relative_AK5Calo.txt:JetMETObjects/data/Summer09_7TeV_ReReco332_L3Absolute_AK5Calo.txt",L2L3);  
  //  create("Spring10_AK5Calo","JetMETObjects/data/Spring10_L2Relative_AK5Calo.txt:JetMETObjects/data/Spring10_L3Absolute_AK5Calo.txt",L2L3);
  //  create("Spring10_AK5PF","JetMETObjects/data/Spring10_L2Relative_AK5PF.txt:JetMETObjects/data/Spring10_L3Absolute_AK5PF.txt",L2L3);
  //  create("Spring10_AK5TRK","JetMETObjects/data/Spring10_L2Relative_AK5TRK.txt:JetMETObjects/data/Spring10_L3Absolute_AK5TRK.txt",L2L3);
  //  create("Spring10_AK5JPT","JetMETObjects/data/Spring10_L2Relative_AK5JPT.txt:JetMETObjects/data/Spring10_L3Absolute_AK5JPT.txt",L2L3);
  //  create("Spring10_AK5CaloData","JetMETObjects/data/Spring10_L2Relative_AK5Calo.txt:JetMETObjects/data/Spring10_L3Absolute_AK5Calo.txt:JetMETObjects/data/Spring10DataV2_L2L3Residual_AK5Calo.txt",L2L3res); 
  //  create("Spring10_AK5PFData","JetMETObjects/data/Spring10_L2Relative_AK5PF.txt:JetMETObjects/data/Spring10_L3Absolute_AK5PF.txt:JetMETObjects/data/Spring10DataV2_L2L3Residual_AK5PF.txt",L2L3res); 
  //  create("Spring10_AK5JPTData","JetMETObjects/data/Spring10_L2Relative_AK5JPT.txt:JetMETObjects/data/Spring10_L3Absolute_AK5JPT.txt:JetMETObjects/data/Spring10DataV2_L2L3Residual_AK5JPT.txt",L2L3res);
  //  create("Spring10_AK5CaloJW","JetMETObjects/data/Spring10_L2Relative_AK5Calo.txt:JetMETObjects/data/Spring10_L3Absolute_AK5Calo.txt:JetMETObjects/data/L4JW_AK5Calo.txt",L2L3L4);  
  create("Fall10_AK5Calo","JetMETObjects/data/Fall10_L1Offset_AK5Calo.txt:JetMETObjects/data/Fall10_L2Relative_AK5Calo.txt:JetMETObjects/data/Fall10_L3Absolute_AK5Calo.txt:JetMETObjects/data/Jec10V1_L2L3Residual_AK5Calo.txt",L1L2L3res);
  create("Fall10_AK7Calo","JetMETObjects/data/Fall10_L1Offset_AK7Calo.txt:JetMETObjects/data/Fall10_L2Relative_AK7Calo.txt:JetMETObjects/data/Fall10_L3Absolute_AK7Calo.txt:JetMETObjects/data/Jec10V1_L2L3Residual_AK7Calo.txt",L1L2L3res);
  create("Fall10_AK5PF","JetMETObjects/data/Fall10_L1Offset_AK5PF.txt:JetMETObjects/data/Fall10_L2Relative_AK5PF.txt:JetMETObjects/data/Fall10_L3Absolute_AK5PF.txt:JetMETObjects/data/Jec10V1_L2L3Residual_AK5PF.txt",L1L2L3res);
  create("Fall10_AK7PF","JetMETObjects/data/Fall10_L1Offset_AK7PF.txt:JetMETObjects/data/Fall10_L2Relative_AK7PF.txt:JetMETObjects/data/Fall10_L3Absolute_AK7PF.txt:JetMETObjects/data/Jec10V1_L2L3Residual_AK7PF.txt",L1L2L3res);
  create("Fall10_AK5CaloNoOffset","JetMETObjects/data/Fall10_L2Relative_AK5Calo.txt:JetMETObjects/data/Fall10_L3Absolute_AK5Calo.txt:JetMETObjects/data/Jec10V1_L2L3Residual_AK5Calo.txt",L2L3res);
  create("Fall10_AK5PFNoOffset","JetMETObjects/data/Fall10_L2Relative_AK5PF.txt:JetMETObjects/data/Fall10_L3Absolute_AK5PF.txt:JetMETObjects/data/Jec10V1_L2L3Residual_AK5PF.txt",L2L3res);  
  create("Fall10_AK7CaloNoOffset","JetMETObjects/data/Fall10_L2Relative_AK7Calo.txt:JetMETObjects/data/Fall10_L3Absolute_AK7Calo.txt:JetMETObjects/data/Jec10V1_L2L3Residual_AK7Calo.txt",L2L3res);
  create("Fall10_AK7PFNoOffset","JetMETObjects/data/Fall10_L2Relative_AK7PF.txt:JetMETObjects/data/Fall10_L3Absolute_AK7PF.txt:JetMETObjects/data/Jec10V1_L2L3Residual_AK7PF.txt",L2L3res);
  create("Fall10_AK5JPT","JetMETObjects/data/Fall10_L2Relative_AK5JPT.txt:JetMETObjects/data/Fall10_L3Absolute_AK5JPT.txt:JetMETObjects/data/Jec10V1_L2L3Residual_AK5JPT.txt",L2L3res);
}

JetMETCorFactorsFactory* JetMETCorFactorsFactory::Register::create(const std::string& name, const std::string& files, Levels type) const
{
  JetMETCorFactorsFactory* jmcff = 0;
  try {
    jmcff = new JetMETCorFactorsFactory(name,files,type);
  } 
  catch(std::exception& e) {
    std::cout << "...failed to create " << name << ":\n";
    std::cout << "     " << e.what() << std::endl;
    jmcff = 0;
  } 
  return jmcff;
}
