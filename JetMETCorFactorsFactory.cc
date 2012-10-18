//
//    $Id: JetMETCorFactorsFactory.cc,v 1.31 2012/10/18 14:12:27 rathjd Exp $
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


CorFactors* JetMETCorFactorsFactory::create(const Jet* j,int nPV, double rho, double jetA)
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
  cor_->setRho(rho);
  cor_->setJetA(jetA);

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
!			  levels[1]/levels[0],
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
  case L1L2L3resL4:
    return new CorFactors(levels[0],
			  levels[1]/levels[0],
			  levels[2]/levels[1],
			  levels[3]/levels[2],
			  levels[4]/levels[3],
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
  create("Spring11_AK5Calo","JetMETObjects/data/Fall10_L1Offset_AK5Calo.txt:JetMETObjects/data/Spring11_L2Relative_AK5Calo.txt:JetMETObjects/data/Spring11_L3Absolute_AK5Calo.txt:JetMETObjects/data/Fall10_L2L3Residual_AK5Calo.txt",L1L2L3res);
  create("Spring11_AK5PF","JetMETObjects/data/Fall10_L1Offset_AK5PF.txt:JetMETObjects/data/Spring11_L2Relative_AK5PF.txt:JetMETObjects/data/Spring11_L3Absolute_AK5PF.txt:JetMETObjects/data/Fall10_L2L3Residual_AK5PF.txt",L1L2L3res);

  create("Fall10_Henning_AK5Calo","JetMETObjects/data/Fall10_L1Offset_AK5Calo.txt:JetMETObjects/data/Fall10_L2Relative_AK5Calo.txt:JetMETObjects/data/Fall10_L3Absolute_AK5Calo.txt:JetMETObjects/data/Fall10_1015absTuneZ2_L2L3Residual_AK5Calo.txt",L1L2L3res);
  create("Fall10_Henning_AK5PF","JetMETObjects/data/Fall10_L1Offset_AK5PF.txt:JetMETObjects/data/Fall10_L2Relative_AK5PF.txt:JetMETObjects/data/Fall10_L3Absolute_AK5PF.txt:JetMETObjects/data/Fall10_1015absTuneZ2_L2L3Residual_AK5PF.txt",L1L2L3res);

  create("Fall10_HenningJER_AK5Calo","JetMETObjects/data/Fall10_L1Offset_AK5Calo.txt:JetMETObjects/data/Fall10_L2Relative_AK5Calo.txt:JetMETObjects/data/Fall10_L3Absolute_AK5Calo.txt:JetMETObjects/data/ConstSmear1015absTuneZ2_L2L3Residual_AK5Calo.txt",L1L2L3res);
  create("Fall10_HenningJER_AK5PF","JetMETObjects/data/Fall10_L1Offset_AK5PF.txt:JetMETObjects/data/Fall10_L2Relative_AK5PF.txt:JetMETObjects/data/Fall10_L3Absolute_AK5PF.txt:JetMETObjects/data/ConstSmear1015absTuneZ2_L2L3Residual_AK5PF.txt",L1L2L3res);
  create("Fall10_HenningJER_AK5JPT","JetMETObjects/data/Fall10_L2Relative_AK5JPT.txt:JetMETObjects/data/Fall10_L3Absolute_AK5JPT.txt:JetMETObjects/data/ConstSmear1015absTuneZ2_L2L3Residual_AK5JPT.txt",L2L3res);

  create("Fall10_L4JW_AK5Calo","JetMETObjects/data/Fall10_L1Offset_AK5Calo.txt:JetMETObjects/data/Fall10_L2Relative_AK5Calo.txt:JetMETObjects/data/Fall10_L3Absolute_AK5Calo.txt:JetMETObjects/data/Fall10_L2L3Residual_AK5Calo.txt:JetMETObjects/data/L4JW_AK5Calo.txt",L1L2L3resL4);
  create("Fall10_AK5Calo","JetMETObjects/data/Fall10_L1Offset_AK5Calo.txt:JetMETObjects/data/Fall10_L2Relative_AK5Calo.txt:JetMETObjects/data/Fall10_L3Absolute_AK5Calo.txt:JetMETObjects/data/Fall10_L2L3Residual_AK5Calo.txt",L1L2L3res);
  create("Fall10_AK7Calo","JetMETObjects/data/Fall10_L1Offset_AK7Calo.txt:JetMETObjects/data/Fall10_L2Relative_AK7Calo.txt:JetMETObjects/data/Fall10_L3Absolute_AK7Calo.txt:JetMETObjects/data/Fall10_L2L3Residual_AK5Calo.txt",L1L2L3res);
  create("Fall10_AK5PF","JetMETObjects/data/Fall10_L1Offset_AK5PF.txt:JetMETObjects/data/Fall10_L2Relative_AK5PF.txt:JetMETObjects/data/Fall10_L3Absolute_AK5PF.txt:JetMETObjects/data/Fall10_L2L3Residual_AK5PF.txt",L1L2L3res);
  create("Fall10_AK7PF","JetMETObjects/data/Fall10_L1Offset_AK7PF.txt:JetMETObjects/data/Fall10_L2Relative_AK7PF.txt:JetMETObjects/data/Fall10_L3Absolute_AK7PF.txt:JetMETObjects/data/Fall10_L2L3Residual_AK5PF.txt",L1L2L3res); 
  create("Fall10_AK5CaloNoOffset","JetMETObjects/data/Fall10_L2Relative_AK5Calo.txt:JetMETObjects/data/Fall10_L3Absolute_AK5Calo.txt:JetMETObjects/data/Fall10_L2L3Residual_AK5Calo.txt",L2L3res);
  create("Fall10_AK5PFNoOffset","JetMETObjects/data/Fall10_L2Relative_AK5PF.txt:JetMETObjects/data/Fall10_L3Absolute_AK5PF.txt:JetMETObjects/data/Fall10_L2L3Residual_AK5PF.txt",L2L3res);
  create("Fall10_AK7CaloNoOffset","JetMETObjects/data/Fall10_L2Relative_AK7Calo.txt:JetMETObjects/data/Fall10_L3Absolute_AK7Calo.txt:JetMETObjects/data/Fall10_L2L3Residual_AK7Calo.txt",L2L3res);
  create("Fall10_AK7PFNoOffset","JetMETObjects/data/Fall10_L2Relative_AK7PF.txt:JetMETObjects/data/Fall10_L3Absolute_AK7PF.txt:JetMETObjects/data/Fall10_L2L3Residual_AK7PF.txt",L2L3res);
  create("Fall10_AK5JPT","JetMETObjects/data/Fall10_L2Relative_AK5JPT.txt:JetMETObjects/data/Fall10_L3Absolute_AK5JPT.txt:JetMETObjects/data/Fall10_L2L3Residual_AK5JPT.txt",L2L3res);
  create("Su11_He_AK5Calo","JetMETObjects/data/GR_R_42_V19_AK5Calo_L1Offset.txt:JetMETObjects/data/GR_R_42_V19_AK5Calo_L2Relative.txt:JetMETObjects/data/GR_R_42_V19_AK5Calo_L3Absolute.txt:JetMETObjects/data/GR_R_42_V19_AK5Calo_L2L3Residual.txt",L1L2L3res);
  create("Su11_He_AK5JPT","JetMETObjects/data/GR_R_42_V19_AK5JPT_L1Offset.txt:JetMETObjects/data/GR_R_42_V19_AK5JPT_L2Relative.txt:JetMETObjects/data/GR_R_42_V19_AK5JPT_L3Absolute.txt:JetMETObjects/data/GR_R_42_V19_AK5JPT_L2L3Residual.txt",L1L2L3res);
  //updated 14.11.2011 (see below)  create("Su11_He_AK5PF","JetMETObjects/data/GR_R_42_V19_AK5PF_L1FastJet.txt:JetMETObjects/data/GR_R_42_V19_AK5PF_L2Relative.txt:JetMETObjects/data/GR_R_42_V19_AK5PF_L3Absolute.txt:JetMETObjects/data/GR_R_42_V19_AK5PF_L2L3Residual.txt",L1L2L3res);
  create("Su11_He_AK5PF","JetMETObjects/data/GR_R_42_V19_AK5PF_L1FastJet.txt:JetMETObjects/data/GR_R_42_V19_AK5PF_L2Relative.txt:JetMETObjects/data/GR_R_42_V19_AK5PF_L3Absolute.txt:JetMETObjects/data/base_RR_PTDEP_pt_L2L3Residual_AK5PF.txt",L1L2L3res);
  //  create("Su11_He_AK5JPT","JetMETObjects/data/GR_R_42_V12_AK5JPT_L1Offset.txt:JetMETObjects/data/GR_R_42_V12_L2Relative_AK5JPT.txt:JetMETObjects/data/GR_R_42_V12_L3Absolute_AK5JPT.txt:JetMETObjects/data/11dcomb_Su11MC_HF_L2L3Residual_AK5JPT.txt",L1L2L3res);
  //  create("Su11_He_AK5PF","JetMETObjects/data/GR_R_42_V12_AK5PF_L1FastJet.txt:JetMETObjects/data/GR_R_42_V12_L2Relative_AK5PF.txt:JetMETObjects/data/GR_R_42_V12_L3Absolute_AK5PF.txt:JetMETObjects/data/11dcomb_Su11MC_HF_L2L3Residual_AK5PF.txt",L1L2L3res);
  create("F11_He_AK5Calo","JetMETObjects/data/Summer11_L1Offset_AK5Calo_SCALED.txt:JetMETObjects/data/Summer11_Jec11_V7_AK5Calo_L2Relative.txt:JetMETObjects/data/Summer11_Jec11_V7_AK5Calo_L3Absolute.txt:JetMETObjects/data/11_11_24_no_JERcorr_RR_L2L3Residual_AK5Calo.txt",L1L2L3res);
  create("F11_He_AK5PF","JetMETObjects/data/Summer11_L1FastJet_AK5PF.txt:JetMETObjects/data/Summer11_Jec11_V7_AK5PF_L2Relative.txt:JetMETObjects/data/Summer11_Jec11_V7_AK5PF_L3Absolute.txt:JetMETObjects/data/11_11_24_base_L2L3Residual_AK5PF.txt",L1L2L3res);
  //  create("F11_He_AK5JPT","JetMETObjects/data/GR_R_42_V19_AK5Calo_L1Offset.txt:JetMETObjects/data/Summer11_Jec11_V7_AK5JPT_L2Relative.txt:JetMETObjects/data/Summer11_Jec11_V7_AK5JPT_L3Absolute.txt:JetMETObjects/data/GR_R_42_V19_AK5JPT_L2L3Residual.txt",L1L2L3res);
    create("F11_He_AK5JPT","JetMETObjects/data/GR_R_42_V19_AK5Calo_L1Offset.txt:JetMETObjects/data/GR_R_42_V19_AK5JPT_L2Relative.txt:JetMETObjects/data/GR_R_42_V19_AK5JPT_L3Absolute.txt:JetMETObjects/data/GR_R_42_V19_AK5JPT_L2L3Residual.txt",L1L2L3res);

//  create("F11DB_He_AK5Calo","JetMETObjects/data/DB_Prelim_Jec11_V9_AK5Calo_L1Offset.txt:JetMETObjects/data/DB_Prelim_Jec11_V9_AK5Calo_L2Relative.txt:JetMETObjects/data/DB_Prelim_Jec11_V9_AK5Calo_L3Absolute.txt:JetMETObjects/data/Y_sm_fine_RR_ALL_EASY_L2L3Residual_AK5Calo.txt",L1L2L3res);
//  create("F11DB_He_AK5PF","JetMETObjects/data/DB_Prelim_Jec11_V9_AK5PF_L1FastJet.txt:JetMETObjects/data/DB_Prelim_Jec11_V9_AK5PF_L2Relative.txt:JetMETObjects/data/DB_Prelim_Jec11_V9_AK5PF_L3Absolute.txt:JetMETObjects/data/Y_sm_fine_RR_ALL_EASY_L2L3Residual_AK5PF.txt",L1L2L3res);
//  create("F11DB_He_AK5JPT","JetMETObjects/data/DB_Prelim_Jec11_V9_AK5JPT_L1Offset.txt:JetMETObjects/data/DB_Prelim_Jec11_V9_AK5JPT_L2Relative.txt:JetMETObjects/data/DB_Prelim_Jec11_V9_AK5JPT_L3Absolute.txt:JetMETObjects/data/Y_sm_fine_RR_ALL_EASY_L2L3Residual_AK5JPT.txt",L1L2L3res);

//updated 11/12/10
//  create("F11DB_He_AK5Calo","JetMETObjects/data/DB_Prelim_Jec11_V9_AK5Calo_L1Offset.txt:JetMETObjects/data/DB_Prelim_Jec11_V9_AK5Calo_L2Relative.txt:JetMETObjects/data/DB_Prelim_Jec11_V9_AK5Calo_L3Absolute.txt:JetMETObjects/data/11DB_MPF_imp_kFSRAbs_AbsetaFIX_L2L3Residual_AK5Calo.txt",L1L2L3res);
//  create("F11DB_He_AK5PF","JetMETObjects/data/DB_Prelim_Jec11_V9_AK5PF_L1FastJet.txt:JetMETObjects/data/DB_Prelim_Jec11_V9_AK5PF_L2Relative.txt:JetMETObjects/data/DB_Prelim_Jec11_V9_AK5PF_L3Absolute.txt:JetMETObjects/data/11DB_MPF_imp_kFSRAbs_AbsetaFIX_L2L3Residual_AK5PF.txt",L1L2L3res);
//  create("F11DB_He_AK5JPT","JetMETObjects/data/DB_Prelim_Jec11_V9_AK5JPT_L1Offset.txt:JetMETObjects/data/DB_Prelim_Jec11_V9_AK5JPT_L2Relative.txt:JetMETObjects/data/DB_Prelim_Jec11_V9_AK5JPT_L3Absolute.txt:JetMETObjects/data/11DB_MPF_imp_kFSRAbs_AbsetaFIX_L2L3Residual_AK5JPT.txt",L1L2L3res);

//updated 11/12/13
  create("F11DB_He_AK5Calo","JetMETObjects/data/Jec11_V11_sqlfile_AK5Calo_L1Offset.txt:JetMETObjects/data/Jec11_V11_sqlfile_AK5Calo_L2Relative.txt:JetMETObjects/data/Jec11_V11_sqlfile_AK5Calo_L3Absolute.txt:JetMETObjects/data/Jec11_V11_sqlfile_AK5Calo_L2L3Residual.txt",L1L2L3res);
  create("F11DB_He_AK5PF","JetMETObjects/data/Jec11_V11_sqlfile_AK5PF_L1FastJet.txt:JetMETObjects/data/Jec11_V11_sqlfile_AK5PF_L2Relative.txt:JetMETObjects/data/Jec11_V11_sqlfile_AK5PF_L3Absolute.txt:JetMETObjects/data/Jec11_V11_sqlfile_AK5PF_L2L3Residual.txt",L1L2L3res);
  //  create("F11DB_He_AK5JPT","JetMETObjects/data/Jec11_V11_sqlfile_AK5JPT_L1Offset.txt:JetMETObjects/data/Jec11_V11_sqlfile_AK5JPT_L2Relative.txt:JetMETObjects/data/Jec11_V11_sqlfile_AK5JPT_L3Absolute.txt:JetMETObjects/data/Jec11_V11_sqlfile_AK5JPT_L2L3Residual.txt",L1L2L3res);
  //FIX 11/12/18 as there was no L1Offset file in the sqlite-file
  create("F11DB_He_AK5JPT","JetMETObjects/data/DB_Prelim_Jec11_V9_AK5JPT_L1Offset.txt:JetMETObjects/data/Jec11_V11_sqlfile_AK5JPT_L2Relative.txt:JetMETObjects/data/Jec11_V11_sqlfile_AK5JPT_L3Absolute.txt:JetMETObjects/data/Jec11_V11_sqlfile_AK5JPT_L2L3Residual.txt",L1L2L3res);
  create("F11DB_He_AK7PF","JetMETObjects/data/Jec11_V11_sqlfile_AK5PF_L1FastJet.txt:JetMETObjects/data/Jec11_V11_sqlfile_AK7PF_L2Relative.txt:JetMETObjects/data/Jec11_V11_sqlfile_AK7PF_L3Absolute.txt:JetMETObjects/data/Jec11_V11_sqlfile_AK5PF_L2L3Residual.txt",L1L2L3res);
  create("F11DB_He_AK7Calo","JetMETObjects/data/Jec11_V11_sqlfile_AK5Calo_L1FastJet.txt:JetMETObjects/data/Jec11_V11_sqlfile_AK7Calo_L2Relative.txt:JetMETObjects/data/Jec11_V11_sqlfile_AK7Calo_L3Absolute.txt:JetMETObjects/data/Jec11_V11_sqlfile_AK5Calo_L2L3Residual.txt",L1L2L3res);

  //updated 29/03/2012
  create("Final2011_AK5Calo","JetMETObjects/data/GR_R_44_V14_AK5Calo_L1Offset.txt:JetMETObjects/data/GR_R_44_V14_AK5Calo_L2Relative.txt:JetMETObjects/data/GR_R_44_V14_AK5Calo_L3Absolute.txt:JetMETObjects/data/GR_R_44_V14_AK5Calo_L2L3Residual.txt",L1L2L3res);
  create("Final2011_AK5JPT","JetMETObjects/data/GR_R_44_V14_AK5JPT_L1Offset.txt:JetMETObjects/data/GR_R_44_V14_AK5JPT_L2Relative.txt:JetMETObjects/data/GR_R_44_V14_AK5JPT_L3Absolute.txt:JetMETObjects/data/GR_R_44_V14_AK5JPT_L2L3Residual.txt",L1L2L3res);
  create("Final2011_AK5PF","JetMETObjects/data/GR_R_44_V14_AK5PF_L1Offset.txt:JetMETObjects/data/GR_R_44_V14_AK5PF_L2Relative.txt:JetMETObjects/data/GR_R_44_V14_AK5PF_L3Absolute.txt:JetMETObjects/data/GR_R_44_V14_AK5PF_L2L3Residual.txt",L1L2L3res);
  create("Final2011_AK5PFchs","JetMETObjects/data/GR_R_44_V14_AK5PFchs_L1Offset.txt:JetMETObjects/data/GR_R_44_V14_AK5PFchs_L2Relative.txt:JetMETObjects/data/GR_R_44_V14_AK5PFchs_L3Absolute.txt:JetMETObjects/data/GR_R_44_V14_AK5PFchs_L2L3Residual.txt",L1L2L3res);
  create("Final2011_AK7Calo","JetMETObjects/data/GR_R_44_V14_AK7Calo_L1Offset.txt:JetMETObjects/data/GR_R_44_V14_AK7Calo_L2Relative.txt:JetMETObjects/data/GR_R_44_V14_AK7Calo_L3Absolute.txt:JetMETObjects/data/GR_R_44_V14_AK7Calo_L2L3Residual.txt",L1L2L3res);
  create("Final2011_AK7JPT","JetMETObjects/data/GR_R_44_V14_AK7JPT_L1Offset.txt:JetMETObjects/data/GR_R_44_V14_AK7JPT_L2Relative.txt:JetMETObjects/data/GR_R_44_V14_AK7JPT_L3Absolute.txt:JetMETObjects/data/GR_R_44_V14_AK7JPT_L2L3Residual.txt",L1L2L3res);
  create("Final2011_AK7PF","JetMETObjects/data/GR_R_44_V14_AK7PF_L1Offset.txt:JetMETObjects/data/GR_R_44_V14_AK7PF_L2Relative.txt:JetMETObjects/data/GR_R_44_V14_AK7PF_L3Absolute.txt:JetMETObjects/data/GR_R_44_V14_AK7PF_L2L3Residual.txt",L1L2L3res);


  create("Final2011_AK5FastCalo","JetMETObjects/data/GR_R_44_V14_AK5Calo_L1FastJet.txt:JetMETObjects/data/GR_R_44_V14_AK5Calo_L2Relative.txt:JetMETObjects/data/GR_R_44_V14_AK5Calo_L3Absolute.txt:JetMETObjects/data/GR_R_44_V14_AK5Calo_L2L3Residual.txt",L1L2L3res);
  create("Final2011_AK5FastJPT","JetMETObjects/data/GR_R_44_V14_AK5JPT_L1FastJet.txt:JetMETObjects/data/GR_R_44_V14_AK5JPT_L2Relative.txt:JetMETObjects/data/GR_R_44_V14_AK5JPT_L3Absolute.txt:JetMETObjects/data/GR_R_44_V14_AK5JPT_L2L3Residual.txt",L1L2L3res);
  create("Final2011_AK5FastPF","JetMETObjects/data/GR_R_44_V14_AK5PF_L1FastJet.txt:JetMETObjects/data/GR_R_44_V14_AK5PF_L2Relative.txt:JetMETObjects/data/GR_R_44_V14_AK5PF_L3Absolute.txt:JetMETObjects/data/GR_R_44_V14_AK5PF_L2L3Residual.txt",L1L2L3res);
  create("Final2011_AK5FastPFchs","JetMETObjects/data/GR_R_44_V14_AK5PFchs_L1FastJet.txt:JetMETObjects/data/GR_R_44_V14_AK5PFchs_L2Relative.txt:JetMETObjects/data/GR_R_44_V14_AK5PFchs_L3Absolute.txt:JetMETObjects/data/GR_R_44_V14_AK5PFchs_L2L3Residual.txt",L1L2L3res);
  create("Final2011_AK7FastCalo","JetMETObjects/data/GR_R_44_V14_AK7Calo_L1FastJet.txt:JetMETObjects/data/GR_R_44_V14_AK7Calo_L2Relative.txt:JetMETObjects/data/GR_R_44_V14_AK7Calo_L3Absolute.txt:JetMETObjects/data/GR_R_44_V14_AK7Calo_L2L3Residual.txt",L1L2L3res);
  create("Final2011_AK7FastJPT","JetMETObjects/data/GR_R_44_V14_AK7JPT_L1FastJet.txt:JetMETObjects/data/GR_R_44_V14_AK7JPT_L2Relative.txt:JetMETObjects/data/GR_R_44_V14_AK7JPT_L3Absolute.txt:JetMETObjects/data/GR_R_44_V14_AK7JPT_L2L3Residual.txt",L1L2L3res);
  create("Final2011_AK7FastPF","JetMETObjects/data/GR_R_44_V14_AK7PF_L1FastJet.txt:JetMETObjects/data/GR_R_44_V14_AK7PF_L2Relative.txt:JetMETObjects/data/GR_R_44_V14_AK7PF_L3Absolute.txt:JetMETObjects/data/GR_R_44_V14_AK7PF_L2L3Residual.txt",L1L2L3res);



  create("2012V7_AK5FastPF","/scratch/hh/current/cms/user/kirschen/JECTextFiles/Jec12_V7_AK5PFLocal_L1FastJet.txt:/scratch/hh/current/cms/user/kirschen/JECTextFiles/Jec12_V7_AK5PFLocal_L2Relative.txt:/scratch/hh/current/cms/user/kirschen/JECTextFiles/Jec12_V7_AK5PFLocal_L3Absolute.txt:/scratch/hh/current/cms/user/kirschen/JECTextFiles/2012_Abseta_L2L3Residual_AK5PF.txt",L1L2L3res);

  create("2012V8_AK5FastPF","/scratch/hh/current/cms/user/kirschen/JECTextFiles/Summer12_V1/Summer12_V1_DATA_AK5PFLocal_L1FastJet.txt:/scratch/hh/current/cms/user/kirschen/JECTextFiles/Summer12_V1/Summer12_V1_DATA_AK5PFLocal_L2Relative.txt:/scratch/hh/current/cms/user/kirschen/JECTextFiles/Summer12_V1/Summer12_V1_DATA_AK5PFLocal_L3Absolute.txt:/scratch/hh/current/cms/user/kirschen/JECTextFiles/Summer12_V1/Summer12_V1_DATA_AK5PFLocal_L2L3Residual.txt",L1L2L3res);
  create("2012V8_AK5FastPF_MC","/scratch/hh/current/cms/user/kirschen/JECTextFiles/Summer12_V1/Summer12_V1_MC_AK5PFLocal_L1FastJet.txt:/scratch/hh/current/cms/user/kirschen/JECTextFiles/Summer12_V1/Summer12_V1_MC_AK5PFLocal_L2Relative.txt:/scratch/hh/current/cms/user/kirschen/JECTextFiles/Summer12_V1/Summer12_V1_MC_AK5PFLocal_L3Absolute.txt",L1L2L3);
 
create("2012V8_AK5PFCHS","/scratch/hh/current/cms/user/kirschen/JECTextFiles/Summer12_V1/Summer12_V1_DATA_AK5PFchsLocal_L1FastJet.txt:/scratch/hh/current/cms/user/kirschen/JECTextFiles/Summer12_V1/Summer12_V1_DATA_AK5PFchsLocal_L2Relative.txt:/scratch/hh/current/cms/user/kirschen/JECTextFiles/Summer12_V1/Summer12_V1_DATA_AK5PFchsLocal_L3Absolute.txt:/scratch/hh/current/cms/user/kirschen/JECTextFiles/Summer12_V1/Summer12_V1_DATA_AK5PFchsLocal_L2L3Residual.txt",L1L2L3res);
 
create("2012V8_AK5PFCHS_MC","/scratch/hh/current/cms/user/kirschen/JECTextFiles/Summer12_V1/Summer12_V1_MC_AK5PFchsLocal_L1FastJet.txt:/scratch/hh/current/cms/user/kirschen/JECTextFiles/Summer12_V1/Summer12_V1_MC_AK5PFchsLocal_L2Relative.txt:/scratch/hh/current/cms/user/kirschen/JECTextFiles/Summer12_V1/Summer12_V1_MC_AK5PFchsLocal_L3Absolute.txt",L1L2L3);


  create("2012V8_AK5FastCalo","/scratch/hh/current/cms/user/kirschen/JECTextFiles/Summer12_V1/Summer12_V1_DATA_AK5CaloLocal_L1FastJet.txt:/scratch/hh/current/cms/user/kirschen/JECTextFiles/Summer12_V1/Summer12_V1_DATA_AK5CaloLocal_L2Relative.txt:/scratch/hh/current/cms/user/kirschen/JECTextFiles/Summer12_V1/Summer12_V1_DATA_AK5CaloLocal_L3Absolute.txt:/scratch/hh/current/cms/user/kirschen/JECTextFiles/Summer12_V1/Summer12_V1_DATA_AK5CaloLocal_L2L3Residual.txt",L1L2L3res);
  create("2012V8_AK5FastCalo_MC","/scratch/hh/current/cms/user/kirschen/JECTextFiles/Summer12_V1/Summer12_V1_MC_AK5CaloLocal_L1FastJet.txt:/scratch/hh/current/cms/user/kirschen/JECTextFiles/Summer12_V1/Summer12_V1_MC_AK5CaloLocal_L2Relative.txt:/scratch/hh/current/cms/user/kirschen/JECTextFiles/Summer12_V1/Summer12_V1_MC_AK5CaloLocal_L3Absolute.txt",L1L2L3);
  create("2012V8_AK5Calo","/scratch/hh/current/cms/user/kirschen/JECTextFiles/Summer12_V1/Summer12_V1_DATA_AK5CaloLocal_L1Offset.txt:/scratch/hh/current/cms/user/kirschen/JECTextFiles/Summer12_V1/Summer12_V1_DATA_AK5CaloLocal_L2Relative.txt:/scratch/hh/current/cms/user/kirschen/JECTextFiles/Summer12_V1/Summer12_V1_DATA_AK5CaloLocal_L3Absolute.txt:/scratch/hh/current/cms/user/kirschen/JECTextFiles/Summer12_V1/Summer12_V1_DATA_AK5CaloLocal_L2L3Residual.txt",L1L2L3res);
  create("2012V8_AK5Calo_MC","/scratch/hh/current/cms/user/kirschen/JECTextFiles/Summer12_V1/Summer12_V1_MC_AK5CaloLocal_L1Offset.txt:/scratch/hh/current/cms/user/kirschen/JECTextFiles/Summer12_V1/Summer12_V1_MC_AK5CaloLocal_L2Relative.txt:/scratch/hh/current/cms/user/kirschen/JECTextFiles/Summer12_V1/Summer12_V1_MC_AK5CaloLocal_L3Absolute.txt",L1L2L3);


  create("2012V8_AK5FastJPT","/scratch/hh/current/cms/user/kirschen/JECTextFiles/Summer12_V1/Summer12_V1_DATA_AK5JPTLocal_L1FastJet.txt:/scratch/hh/current/cms/user/kirschen/JECTextFiles/Summer12_V1/Summer12_V1_DATA_AK5JPTLocal_L2Relative.txt:/scratch/hh/current/cms/user/kirschen/JECTextFiles/Summer12_V1/Summer12_V1_DATA_AK5JPTLocal_L3Absolute.txt:/scratch/hh/current/cms/user/kirschen/JECTextFiles/Summer12_V1/Summer12_V1_DATA_AK5JPTLocal_L2L3Residual.txt",L1L2L3res);
  create("2012V8_AK5FastJPT_MC","/scratch/hh/current/cms/user/kirschen/JECTextFiles/Summer12_V1/Summer12_V1_MC_AK5JPTLocal_L1FastJet.txt:/scratch/hh/current/cms/user/kirschen/JECTextFiles/Summer12_V1/Summer12_V1_MC_AK5JPTLocal_L2Relative.txt:/scratch/hh/current/cms/user/kirschen/JECTextFiles/Summer12_V1/Summer12_V1_MC_AK5JPTLocal_L3Absolute.txt",L1L2L3);
  create("2012V8_AK5JPT","/scratch/hh/current/cms/user/kirschen/JECTextFiles/Summer12_V1/Summer12_V1_DATA_AK5JPTLocal_L1Offset.txt:/scratch/hh/current/cms/user/kirschen/JECTextFiles/Summer12_V1/Summer12_V1_DATA_AK5JPTLocal_L2Relative.txt:/scratch/hh/current/cms/user/kirschen/JECTextFiles/Summer12_V1/Summer12_V1_DATA_AK5JPTLocal_L3Absolute.txt:/scratch/hh/current/cms/user/kirschen/JECTextFiles/Summer12_V1/Summer12_V1_DATA_AK5JPTLocal_L2L3Residual.txt",L1L2L3res);
  create("2012V8_AK5JPT_MC","/scratch/hh/current/cms/user/kirschen/JECTextFiles/Summer12_V1/Summer12_V1_MC_AK5JPTLocal_L1Offset.txt:/scratch/hh/current/cms/user/kirschen/JECTextFiles/Summer12_V1/Summer12_V1_MC_AK5JPTLocal_L2Relative.txt:/scratch/hh/current/cms/user/kirschen/JECTextFiles/Summer12_V1/Summer12_V1_MC_AK5JPTLocal_L3Absolute.txt",L1L2L3);



  create("2012SQLV3_AK5FastPF","/scratch/hh/current/cms/user/kirschen/JECTextFiles/Summer12_V3/Summer12_V3_DATA_AK5PFLocal_L1FastJet.txt:/scratch/hh/current/cms/user/kirschen/JECTextFiles/Summer12_V3/Summer12_V3_DATA_AK5PFLocal_L2Relative.txt:/scratch/hh/current/cms/user/kirschen/JECTextFiles/Summer12_V3/Summer12_V3_DATA_AK5PFLocal_L3Absolute.txt:/scratch/hh/current/cms/user/kirschen/JECTextFiles/Summer12_V3/Summer12_V3_DATA_AK5PFLocal_L2L3Residual.txt",L1L2L3res);
  create("2012SQLV3_AK5FastPF_MC","/scratch/hh/current/cms/user/kirschen/JECTextFiles/Summer12_V3/Summer12_V3_MC_AK5PFLocal_L1FastJet.txt:/scratch/hh/current/cms/user/kirschen/JECTextFiles/Summer12_V3/Summer12_V3_MC_AK5PFLocal_L2Relative.txt:/scratch/hh/current/cms/user/kirschen/JECTextFiles/Summer12_V3/Summer12_V3_MC_AK5PFLocal_L3Absolute.txt",L1L2L3);

  create("2012SQLV3_AK5FastPFCHS","/scratch/hh/current/cms/user/kirschen/JECTextFiles/Summer12_V3/Summer12_V3_DATA_AK5PFchsLocal_L1FastJet.txt:/scratch/hh/current/cms/user/kirschen/JECTextFiles/Summer12_V3/Summer12_V3_DATA_AK5PFchsLocal_L2Relative.txt:/scratch/hh/current/cms/user/kirschen/JECTextFiles/Summer12_V3/Summer12_V3_DATA_AK5PFchsLocal_L3Absolute.txt:/scratch/hh/current/cms/user/kirschen/JECTextFiles/Summer12_V3/Summer12_V3_DATA_AK5PFchsLocal_L2L3Residual.txt",L1L2L3res);
  create("2012SQLV3_AK5FastPFCHS_MC","/scratch/hh/current/cms/user/kirschen/JECTextFiles/Summer12_V3/Summer12_V3_MC_AK5PFchsLocal_L1FastJet.txt:/scratch/hh/current/cms/user/kirschen/JECTextFiles/Summer12_V3/Summer12_V3_MC_AK5PFchsLocal_L2Relative.txt:/scratch/hh/current/cms/user/kirschen/JECTextFiles/Summer12_V3/Summer12_V3_MC_AK5PFchsLocal_L3Absolute.txt",L1L2L3);


  create("2012SQLV7_AK5FastPF","/scratch/hh/current/cms/user/kirschen/JECTextFiles/Summer12_V7/Summer12_V7_DATA_L1FastJet_AK5PF.txt:/scratch/hh/current/cms/user/kirschen/JECTextFiles/Summer12_V7/Summer12_V7_DATA_L2Relative_AK5PF.txt:/scratch/hh/current/cms/user/kirschen/JECTextFiles/Summer12_V7/Summer12_V7_DATA_L3Absolute_AK5PF.txt:/scratch/hh/current/cms/user/kirschen/JECTextFiles/Summer12_V7/Summer12_V7_DATA_L2L3Residual_AK5PF.txt",L1L2L3res);
  create("2012SQLV7_AK5FastPF_MC","/scratch/hh/current/cms/user/kirschen/JECTextFiles/Summer12_V7/Summer12_V7_MC_L1FastJet_AK5PF.txt:/scratch/hh/current/cms/user/kirschen/JECTextFiles/Summer12_V7/Summer12_V7_MC_L2Relative_AK5PF.txt:/scratch/hh/current/cms/user/kirschen/JECTextFiles/Summer12_V7/Summer12_V7_MC_L3Absolute_AK5PF.txt",L1L2L3);

  create("2012SQLV7_AK5FastPFCHS","/scratch/hh/current/cms/user/kirschen/JECTextFiles/Summer12_V7/Summer12_V7_DATA_L1FastJet_AK5PFchs.txt:/scratch/hh/current/cms/user/kirschen/JECTextFiles/Summer12_V7/Summer12_V7_DATA_L2Relative_AK5PFchs.txt:/scratch/hh/current/cms/user/kirschen/JECTextFiles/Summer12_V7/Summer12_V7_DATA_L3Absolute_AK5PFchs.txt:/scratch/hh/current/cms/user/kirschen/JECTextFiles/Summer12_V7/Summer12_V7_DATA_L2L3Residual_AK5PFchs.txt",L1L2L3res);
  create("2012SQLV7_AK5FastPFCHS_MC","/scratch/hh/current/cms/user/kirschen/JECTextFiles/Summer12_V7/Summer12_V7_MC_L1FastJet_AK5PFchs.txt:/scratch/hh/current/cms/user/kirschen/JECTextFiles/Summer12_V7/Summer12_V7_MC_L2Relative_AK5PFchs.txt:/scratch/hh/current/cms/user/kirschen/JECTextFiles/Summer12_V7/Summer12_V7_MC_L3Absolute_AK5PFchs.txt",L1L2L3);

  create("2012SQLV7_AK5Calo","/scratch/hh/current/cms/user/kirschen/JECTextFiles/Summer12_V7/Summer12_V7_DATA_L1Offset_AK5Calo.txt:/scratch/hh/current/cms/user/kirschen/JECTextFiles/Summer12_V7/Summer12_V7_DATA_L2Relative_AK5Calo.txt:/scratch/hh/current/cms/user/kirschen/JECTextFiles/Summer12_V7/Summer12_V7_DATA_L3Absolute_AK5Calo.txt:/scratch/hh/current/cms/user/kirschen/JECTextFiles/Summer12_V7/Summer12_V7_DATA_L2L3Residual_AK5Calo.txt",L1L2L3res);
  create("2012SQLV7_AK5Calo_MC","/scratch/hh/current/cms/user/kirschen/JECTextFiles/Summer12_V7/Summer12_V7_MC_L1Offset_AK5Calo.txt:/scratch/hh/current/cms/user/kirschen/JECTextFiles/Summer12_V7/Summer12_V7_MC_L2Relative_AK5Calo.txt:/scratch/hh/current/cms/user/kirschen/JECTextFiles/Summer12_V7/Summer12_V7_MC_L3Absolute_AK5Calo.txt",L1L2L3);

  create("2012SQLV7_AK5JPT","/scratch/hh/current/cms/user/kirschen/JECTextFiles/Summer12_V7/Summer12_V7_DATA_L1FastJet_AK5JPT.txt:/scratch/hh/current/cms/user/kirschen/JECTextFiles/Summer12_V7/Summer12_V7_DATA_L2Relative_AK5JPT.txt:/scratch/hh/current/cms/user/kirschen/JECTextFiles/Summer12_V7/Summer12_V7_DATA_L3Absolute_AK5JPT.txt:/scratch/hh/current/cms/user/kirschen/JECTextFiles/Summer12_V7/Summer12_V7_DATA_L2L3Residual_AK5JPT.txt",L1L2L3res);
  create("2012SQLV7_AK5JPT_MC","/scratch/hh/current/cms/user/kirschen/JECTextFiles/Summer12_V7/Summer12_V7_MC_L1FastJet_AK5JPT.txt:/scratch/hh/current/cms/user/kirschen/JECTextFiles/Summer12_V7/Summer12_V7_MC_L2Relative_AK5JPT.txt:/scratch/hh/current/cms/user/kirschen/JECTextFiles/Summer12_V7/Summer12_V7_MC_L3Absolute_AK5JPT.txt",L1L2L3);

 
 
  create("2012FallV1_AK7PF","/scratch/hh/dust/naf/cms/user/rathjd/Calibration/CorrectionFactors/Fall12_V1/Fall12_V1_DATA_L1FastJet_AK7PF.txt:/scratch/hh/dust/naf/cms/user/rathjd/Calibration/CorrectionFactors/Fall12_V1/Fall12_V1_DATA_L2Relative_AK7PF.txt:/scratch/hh/dust/naf/cms/user/rathjd/Calibration/CorrectionFactors/Fall12_V1/Fall12_V1_DATA_L3Absolute_AK7PF.txt:/scratch/hh/current/cms/user/kirschen/JECTextFiles/Summer12_V7/Summer12_V7_DATA_L2L3Residual_AK7PF.txt",L1L2L3res);
  create("2012FallV1_AK7PF_MC","/scratch/hh/dust/naf/cms/user/rathjd/Calibration/CorrectionFactors/Fall12_V1/Fall12_V1_MC_L1FastJet_AK7PF.txt:/scratch/hh/dust/naf/cms/user/rathjd/Calibration/CorrectionFactors/Fall12_V1/Fall12_V1_MC_L2Relative_AK7PF.txt:/scratch/hh/dust/naf/cms/user/rathjd/Calibration/CorrectionFactors/Fall12_V1/Fall12_V1_MC_L3Absolute_AK7PF.txt",L1L2L3);

  create("2012FallV1_AK5FastPF","/scratch/hh/dust/naf/cms/user/rathjd/Calibration/CorrectionFactors/Fall12_V1/Fall12_V1_DATA_L1FastJet_AK5PF.txt:/scratch/hh/dust/naf/cms/user/rathjd/Calibration/CorrectionFactors/Fall12_V1/Fall12_V1_DATA_L2Relative_AK5PF.txt:/scratch/hh/dust/naf/cms/user/rathjd/Calibration/CorrectionFactors/Fall12_V1/Fall12_V1_DATA_L3Absolute_AK5PF.txt://scratch/hh/current/cms/user/kirschen/JECTextFiles/Summer12_V7/Summer12_V7_DATA_L2L3Residual_AK5PF.txt",L1L2L3res);
  create("2012FallV1_AK5FastPF_MC","/scratch/hh/dust/naf/cms/user/rathjd/Calibration/CorrectionFactors/Fall12_V1/Fall12_V1_MC_L1FastJet_AK5PF.txt:/scratch/hh/dust/naf/cms/user/rathjd/Calibration/CorrectionFactors/Fall12_V1/Fall12_V1_MC_L2Relative_AK5PF.txt:/scratch/hh/dust/naf/cms/user/rathjd/Calibration/CorrectionFactors/Fall12_V1/Fall12_V1_MC_L3Absolute_AK5PF.txt",L1L2L3);

  create("2012FallV1_AK5FastPFCHS","/scratch/hh/dust/naf/cms/user/rathjd/Calibration/CorrectionFactors/Fall12_V1/Fall12_V1_DATA_L1FastJet_AK5PFchs.txt:/scratch/hh/dust/naf/cms/user/rathjd/Calibration/CorrectionFactors/Fall12_V1/Fall12_V1_DATA_L2Relative_AK5PFchs.txt:/scratch/hh/dust/naf/cms/user/rathjd/Calibration/CorrectionFactors/Fall12_V1/Fall12_V1_DATA_L3Absolute_AK5PFchs.txt:/scratch/hh/current/cms/user/kirschen/JECTextFiles/Summer12_V7/Summer12_V7_DATA_L2L3Residual_AK5PFchs.txt",L1L2L3res);
  create("2012FallV1_AK5FastPFCHS_MC","/scratch/hh/dust/naf/cms/user/rathjd/Calibration/CorrectionFactors/Fall12_V1/Fall12_V1_MC_L1FastJet_AK5PFchs.txt:/scratch/hh/dust/naf/cms/user/rathjd/Calibration/CorrectionFactors/Fall12_V1/Fall12_V1_MC_L2Relative_AK5PFchs.txt:/scratch/hh/dust/naf/cms/user/rathjd/Calibration/CorrectionFactors/Fall12_V1/Fall12_V1_MC_L3Absolute_AK5PFchs.txt",L1L2L3);

  create("2012FallV1_AK5Calo","/scratch/hh/dust/naf/cms/user/rathjd/Calibration/CorrectionFactors/Fall12_V1/Fall12_V1_DATA_L1Offset_AK5Calo.txt:/scratch/hh/dust/naf/cms/user/rathjd/Calibration/CorrectionFactors/Fall12_V1/Fall12_V1_DATA_L2Relative_AK5Calo.txt:/scratch/hh/dust/naf/cms/user/rathjd/Calibration/CorrectionFactors/Fall12_V1/Fall12_V1_DATA_L3Absolute_AK5Calo.txt:/scratch/hh/current/cms/user/kirschen/JECTextFiles/Summer12_V7/Summer12_V7_DATA_L2L3Residual_AK5Calo.txt",L1L2L3res);
  create("2012FallV1_AK5Calo_MC","/scratch/hh/dust/naf/cms/user/rathjd/Calibration/CorrectionFactors/Fall12_V1/Fall12_V1_MC_L1Offset_AK5Calo.txt:/scratch/hh/dust/naf/cms/user/rathjd/Calibration/CorrectionFactors/Fall12_V1/Fall12_V1_MC_L2Relative_AK5Calo.txt:/scratch/hh/dust/naf/cms/user/rathjd/Calibration/CorrectionFactors/Fall12_V1/Fall12_V1_MC_L3Absolute_AK5Calo.txt",L1L2L3);

  create("2012FallV1_AK5JPT","/scratch/hh/dust/naf/cms/user/rathjd/Calibration/CorrectionFactors/Fall12_V1/Fall12_V1_DATA_L1FastJet_AK5JPT.txt:/scratch/hh/dust/naf/cms/user/rathjd/Calibration/CorrectionFactors/Fall12_V1/Fall12_V1_DATA_L2Relative_AK5JPT.txt:/scratch/hh/dust/naf/cms/user/rathjd/Calibration/CorrectionFactors/Fall12_V1/Fall12_V1_DATA_L3Absolute_AK5JPT.txt:/scratch/hh/current/cms/user/kirschen/JECTextFiles/Summer12_V7/Summer12_V7_DATA_L2L3Residual_AK5JPT.txt",L1L2L3res);
  create("2012FallV1_AK5JPT_MC","/scratch/hh/dust/naf/cms/user/rathjd/Calibration/CorrectionFactors/Fall12_V1/Fall12_V1_MC_L1FastJet_AK5JPT.txt:/scratch/hh/dust/naf/cms/user/rathjd/Calibration/CorrectionFactors/Fall12_V1/Fall12_V1_MC_L2Relative_AK5JPT.txt:/scratch/hh/dust/naf/cms/user/rathjd/Calibration/CorrectionFactors/Fall12_V1/Fall12_V1_MC_L3Absolute_AK5JPT.txt",L1L2L3);
 
  create("2012FallV1_AK7PF","/scratch/hh/dust/naf/cms/user/rathjd/Calibration/CorrectionFactors/Fall12_V1/Fall12_V1_DATA_L1FastJet_AK7PF.txt:/scratch/hh/dust/naf/cms/user/rathjd/Calibration/CorrectionFactors/Fall12_V1/Fall12_V1_DATA_L2Relative_AK7PF.txt:/scratch/hh/dust/naf/cms/user/rathjd/Calibration/CorrectionFactors/Fall12_V1/Fall12_V1_DATA_L3Absolute_AK7PF.txt:/scratch/hh/current/cms/user/kirschen/JECTextFiles/Summer12_V7/Summer12_V7_DATA_L2L3Residual_AK7PF.txt",L1L2L3res); 
  create("2012FallV1_AK7PF_MC","/scratch/hh/dust/naf/cms/user/rathjd/Calibration/CorrectionFactors/Fall12_V1/Fall12_V1_MC_L1FastJet_AK7PF.txt:/scratch/hh/dust/naf/cms/user/rathjd/Calibration/CorrectionFactors/Fall12_V1/Fall12_V1_MC_L2Relative_AK7PF.txt:/scratch/hh/dust/naf/cms/user/rathjd/Calibration/CorrectionFactors/Fall12_V1/Fall12_V1_MC_L3Absolute_AK7PF.txt",L1L2L3);


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
    jmcff = 0;			// 
  } 
  return jmcff;
}
