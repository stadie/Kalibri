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
  if(trignames_.size() != trigthresholds_.size()) {
    std::cerr << "DiJetReader: number of triggers and thresholds mismatch." 
	      << trignames_.size() << " != " << trigthresholds_.size() << std::endl;
    exit(9);
  }
  useSingleJetTriggers_ = config.read<bool>("Use single jet triggers",false);

 
  for(int i = 0, l = trigthresholds_.size() ; i < l ; ++i) {
    ndata_[trigthresholds_[i]] = 0;
    ncontrol_[trigthresholds_[i]] = 0;
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


  std::cout <<"initialize reweight::LumiReWeighting LumiWeights_" << std::endl;
  //      reweight::LumiReWeighting LumiWeights_= reweight::LumiReWeighting("/afs/naf.desy.de/user/k/kirschen/PUDistributions/TruePU_Distributions2012_04_24.root", "/afs/naf.desy.de/user/k/kirschen/scratch/2012_03_PUperHLT/CMSSW_5_2_3_patch4/src/PUData/2012A/MyDataPileupHistogramTruth_AllHLT.root", "PU_profile_TrueSummer12", "pileup");
  //true
  reweight::LumiReWeighting LumiWeights_= reweight::LumiReWeighting("/afs/naf.desy.de/user/k/kirschen/PUDistributions/TruePU_Distributions2012_04_25.root", "/afs/naf.desy.de/user/k/kirschen/scratch/2012_03_PUperHLT/CMSSW_5_2_3_patch4/src/MyDataPileupHistogramTrueAllHLT.root", "PU_profile_TrueSummer12", "pileup");
  //observed
  //      reweight::LumiReWeighting LumiWeights_= reweight::LumiReWeighting("/afs/naf.desy.de/user/k/kirschen/PUDistributions/TruePU_Distributions2012_04_25.root", "/afs/naf.desy.de/user/k/kirschen/scratch/2012_03_PUperHLT/CMSSW_5_2_3_patch4/src/MyDataPileupHistogramObservedAllHLT.root", "PU_profile_ObservedSummer12", "pileup");
  //count events in data and control sample
  std::cout <<"initialized reweight::LumiReWeighting LumiWeights_ successfully?" << std::endl;

 
  int nProcEvts = 0; // Number of processed events
  double WeightSum = 0; // Sum of Weights
  std::vector<Event*>::iterator evt1 = control1.begin();
  for(; evt1 != control1.end(); ++evt1, ++nProcEvts) {
    if( (*evt1)->type() == ParLimit ) continue;
    float nputruth = (*evt1)->nPUTruth();
    double MyWeight = LumiWeights_.ITweight3BX( nputruth );
    //    int nputruth = (*evt1)->nPU();
    //    double MyWeight = LumiWeights_.ITweight( nputruth );
    WeightSum=+(MyWeight * ((*evt1)->weight()));
    if( nputruth < 60 ) {//WARNING: hard-coded 60
      (*evt1)->setWeight( MyWeight * ((*evt1)->weight()) );
    } else {
      std::cerr << "WARNING in PUTruthReweighting::preprocess: Number of PU vertices = " << nputruth << " out of histogram binning." << std::endl;
    }
  }

  std::cout << "  Applied weights for " << nProcEvts << " events in control1 and control2. \n";
  std::cout << "  Average weight was " << WeightSum/nProcEvts << "\n";
  
  std::cout << "end PUTruthReweighting with:" << std::endl;
  std::cout << "  " << data.size() << " events in data" << std::endl;
  std::cout << "  " << control1.size() << " events in control1" << std::endl;


  return ndata_.size();
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