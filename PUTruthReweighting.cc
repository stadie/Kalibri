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
#include "JetTruthEvent.h"
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

  //initializing
  std::cout <<"initialize reweight::LumiReWeighting LumiWeights_" << std::endl;
  std::vector <reweight::LumiReWeighting> LumiWeightsPerTrigger;
  for(int i = 0, l = trigthresholds_.size() ; i < l ; ++i) {
    LumiWeightsPerTrigger.push_back(reweight::LumiReWeighting(("/scratch/hh/current/cms/user/"+TruthWeightingMCDistribution_+"_TrueMCPUDistributions.root").c_str(),("/scratch/hh/current/cms/user/"+TruthWeightingDir_+"/MyDataPileupHistogram"+(TString)trignames_.at(i)+".root").Data(), "pileup", "pileup"));
    std::cout <<"initialized reweight::LumiReWeighting LumiWeights_ successfully..." << i <<std::endl;

  }

 
  int nProcEvts = 0; // Number of processed events
  double WeightSumAfter = 0; // Sum of Weights
  double WeightSumBefore = 0; // Sum of Weights
  std::vector<Event*>::iterator evt1 = control1.begin();
  if( (*evt1)->type() != PtBalance)std::cout << "Warning: No TwoJetsPtBalanceEvent! Take special care when e.g. analyzing JetTruthEvents"; 
  for(; evt1 != control1.end(); ++evt1, ++nProcEvts) {
    if( (*evt1)->type() == ParLimit ) continue;
    float nputruth = (*evt1)->nPUTruth();
    std::map<double,int>::iterator it;
    if( (*evt1)->type() == PtBalance){
    TwoJetsPtBalanceEvent* tje = dynamic_cast<TwoJetsPtBalanceEvent*>((*evt1));
    it = controlTrigger_.lower_bound(tje->triggerPtVariableL2L3(useSingleJetTriggers_));
    }
    else if((*evt1)->type() == GammaJet){ //GammaJet is the data type used for JetTruthEvents
      JetTruthEvent * jte = dynamic_cast<JetTruthEvent*>((*evt1));
      it = controlTrigger_.lower_bound(jte->jet()->corFactors().getL2L3() * jte->jet()->pt());//WARNING, hack: use L2L3-corrected reco pt as pt-variable for trigger thresholds for trutevents
    }
    else{
      std::cout << "Warning: No TwoJetsPtBalanceEvent and no GammaJet/JetTruthEvent!"; 
    }
   if(!(it == controlTrigger_.begin())){
     assert(it != controlTrigger_.begin());
     --it;
   }
   //    std::cout << (*it).second <<" pt: " <<tje->triggerPtVariableL2L3(useSingleJetTriggers_)<< std::endl;
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
 
