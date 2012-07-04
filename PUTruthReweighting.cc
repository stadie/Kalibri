//
//    Class for reweighting PU according to
//    official recipes (2011/2012)
//    Truth PU
//   
#include "PUTruthReweighting.h"

#include "CalibData.h"
#include "ConfigFile.h"
#include "Parameters.h"
#include "TwoJetsPtBalanceEvent.h"
#include "PUReweighting/LumiReweightingStandAlone.h"

#include <iostream>

PUTruthReweighting::PUTruthReweighting(const std::string& configfile, Parameters* param)
  : EventProcessor("PUTruthReweighting",configfile,param)
{  
  ConfigFile config(configfile.c_str()); 
    
  trignames_ = bag_of_string(config.read<std::string>("Di-Jet trigger names",""));
  trigthresholds_ = bag_of<double>(config.read<std::string>("Di-Jet trigger thresholds",""));
  TruthWeightingDir_ = config.read<std::string>("PU TruthWeighting","Cert_2012_190456-191859");
  TruthWeightingMCDistribution_ = config.read<std::string>("PU TruthWeighting MC distribution","TrueSummer12");
  if(trignames_.size() != trigthresholds_.size()) {
    std::cerr << "DiJetReader: number of triggers and thresholds mismatch." 
	      << trignames_.size() << " != " << trigthresholds_.size() << std::endl;
    exit(9);
  }
  useSingleJetTriggers_ = config.read<bool>("Use single jet triggers",false);

 
  for(int i = 0, l = trigthresholds_.size() ; i < l ; ++i) {
    //    ndata_[trigthresholds_[i]] = 1;
    std::cout << i << " "<< trignames_.at(i) << " " << trigthresholds_.at(i) << std::endl;

    controlTrigger_[trigthresholds_[i]] = i;
  }
}
 
PUTruthReweighting::~PUTruthReweighting()
{
}
  

int PUTruthReweighting::preprocess(std::vector<Event*>& data,
			     std::vector<Event*>& control1,
			     std::vector<Event*>& control2)
{

  std::cout << "start PUTruthReweighting with:" << std::endl;
  std::cout << "  " << data.size() << " events in data" << std::endl;
  std::cout << "  " << control1.size() << " events in control1" << std::endl;
  //  std::cout << data.size() << " events in data" << std::endl;


//  std::cout <<"initialize reweight::LumiReWeighting LumiWeights_" << std::endl;
//  //      reweight::LumiReWeighting LumiWeights_= reweight::LumiReWeighting("/afs/naf.desy.de/user/k/kirschen/PUDistributions/TruePU_Distributions2012_04_24.root", "/afs/naf.desy.de/user/k/kirschen/scratch/2012_03_PUperHLT/CMSSW_5_2_3_patch4/src/PUData/2012A/MyDataPileupHistogramTruth_AllHLT.root", "PU_profile_TrueSummer12", "pileup");
//  //true
//  reweight::LumiReWeighting LumiWeights_= reweight::LumiReWeighting("/afs/naf.desy.de/user/k/kirschen/PUDistributions/TruePU_Distributions2012_04_25.root", "/afs/naf.desy.de/user/k/kirschen/scratch/2012_03_PUperHLT/CMSSW_5_2_3_patch4/src/MyDataPileupHistogramTrueAllHLT.root", "PU_profile_TrueSummer12", "pileup");
//  //observed
//  //      reweight::LumiReWeighting LumiWeights_= reweight::LumiReWeighting("/afs/naf.desy.de/user/k/kirschen/PUDistributions/TruePU_Distributions2012_04_25.root", "/afs/naf.desy.de/user/k/kirschen/scratch/2012_03_PUperHLT/CMSSW_5_2_3_patch4/src/MyDataPileupHistogramObservedAllHLT.root", "PU_profile_ObservedSummer12", "pileup");
//  //count events in data and control sample
//  std::cout <<"initialized reweight::LumiReWeighting LumiWeights_ successfully?" << std::endl;

  //initializing
  std::cout <<"initialize reweight::LumiReWeighting LumiWeights_" << std::endl;
  std::vector <reweight::LumiReWeighting> LumiWeightsPerTrigger;
  for(int i = 0, l = trigthresholds_.size() ; i < l ; ++i) {
    //    LumiWeightsPerTrigger.push_back(reweight::LumiReWeighting("/afs/naf.desy.de/user/k/kirschen/PUDistributions/TruePU_Distributions2012_04_25.root",("/afs/naf.desy.de/user/k/kirschen/scratch/2012_03_PUperHLT/CMSSW_5_2_3_patch4/src/PUData/2012AMaxBin60_191859/MyDataPileupHistogram"+(TString)trignames_.at(i)+".root").Data(), "PU_profile_TrueSummer12", "pileup"));
    LumiWeightsPerTrigger.push_back(reweight::LumiReWeighting("/scratch/hh/current/cms/user/kirschen/PUDistributions/TruePU_Distributions.root",("/scratch/hh/current/cms/user/"+TruthWeightingDir_+"/MyDataPileupHistogram"+(TString)trignames_.at(i)+".root").Data(), ("PU_profile_"+TruthWeightingMCDistribution_).c_str(), "pileup"));
    std::cout <<"initialized reweight::LumiReWeighting LumiWeights_ successfully..." << i <<std::endl;

  }

 
  int nProcEvts = 0; // Number of processed events
  double WeightSumAfter = 0; // Sum of Weights
  double WeightSumBefore = 0; // Sum of Weights
  std::vector<Event*>::iterator evt1 = control1.begin();
  for(; evt1 != control1.end(); ++evt1, ++nProcEvts) {
    if( (*evt1)->type() == ParLimit ) continue;
    float nputruth = (*evt1)->nPUTruth();
    if( (*evt1)->type() != PtBalance) std::cout << "Warning: No TwoJetsPtBalanceEvent! ";
    TwoJetsPtBalanceEvent* tje = dynamic_cast<TwoJetsPtBalanceEvent*>((*evt1));
    std::map<double,int>::iterator it = controlTrigger_.lower_bound(TriggerPtVariable(tje));
   if(!(it == controlTrigger_.begin())){
     assert(it != controlTrigger_.begin());
     --it;
   }
   //    std::cout << (*it).second <<" pt: " <<TriggerPtVariable(tje)<< std::endl;
    //    std::cout << (*controlTrigger_.lower_bound(600.)).second <<" pt: " <<"600"<< std::endl;
   double MyWeight = LumiWeightsPerTrigger.at((*it).second).ITweight3BX( nputruth );
    //    double MyWeight = LumiWeights_.ITweight3BX( nputruth );
    //    int nputruth = (*evt1)->nPU();
    //    double MyWeight = LumiWeights_.ITweight( nputruth );
    WeightSumBefore=+(*evt1)->weight();
    WeightSumAfter=+(MyWeight * ((*evt1)->weight()));
    if( nputruth < 60 ) {//WARNING: hard-coded 60
      (*evt1)->setWeight( MyWeight * ((*evt1)->weight()) );
    } else {
      std::cerr << "WARNING in PUTruthReweighting::preprocess: Number of PU vertices = " << nputruth << " out of histogram binning." << std::endl;
    }
  }

  std::cout << "  Applied weights for " << nProcEvts << " events in control1 and control2. \n";
  std::cout << "  Average weight was " << WeightSumAfter/WeightSumBefore << "\n";
  
  std::cout << "end PUTruthReweighting with:" << std::endl;
  std::cout << "  " << data.size() << " events in data" << std::endl;
  std::cout << "  " << control1.size() << " events in control1" << std::endl;


  return nProcEvts;
}
 

double PUTruthReweighting::TriggerPtVariable(Event* event)
{
  if(event->type() != PtBalance) std::cout << "Warning: No TwoJetsPtBalanceEvent! ";
  TwoJetsPtBalanceEvent* tje = dynamic_cast<TwoJetsPtBalanceEvent*>(event);
  if(!useSingleJetTriggers_)return tje->ptDijetCorrL2L3();
  else{
  Jet * j1 = tje->getJet1();
  Jet * j2 = tje->getJet2();

  double ptcorj1,ptcorj2;
  ptcorj1 = j1->corFactors().getL2L3() * j1->pt();
  ptcorj2 = j2->corFactors().getL2L3() * j2->pt();

  if(ptcorj1>ptcorj2)return ptcorj1;
  else return ptcorj2;
  }

}
