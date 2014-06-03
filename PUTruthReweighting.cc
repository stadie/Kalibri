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

#include <iostream>

PUTruthReweighting::PUTruthReweighting(const std::string& configfile, Parameters* param)
  : EventProcessor("PUTruthReweighting",configfile,param)
{  
  if(! isActive()) return;

  ConfigFile config(configfile.c_str()); 
    

  trignames_ = bag_of_string(config.read<std::string>("Di-Jet trigger names",""));
  trigthresholds_ = bag_of<double>(config.read<std::string>("Di-Jet trigger thresholds",""));
  useMCReweightAll_ = config.read<bool>("PU TruthWeighting Reweight all eventvectors (for MC validation)",false);
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

  //initializing
  std::cout <<"initialize reweight::LumiReWeighting LumiWeights_" << std::endl;
  for(int i = 0, l = trigthresholds_.size() ; i < l ; ++i) {
    LumiWeightsPerTrigger_.push_back(reweight::LumiReWeighting(("/nfs/dust/cms/user/"+TruthWeightingMCDistribution_+"_TrueMCPUDistributions.root").c_str(),("/nfs/dust/cms/user/"+TruthWeightingDir_+"/MyDataPileupHistogram"+(TString)trignames_.at(i)+".root").Data(), "pileup", "pileup"));
    std::cout <<"initialized reweight::LumiReWeighting LumiWeights_ successfully..." << i <<std::endl;
    
  }

}
 
PUTruthReweighting::~PUTruthReweighting()
{
}
  

int PUTruthReweighting::reweightEventVector(std::vector<Event*>& evtVector){
  
  if(evtVector.size()==0){
    std::cout << "WARNING: evtVector empty" <<std::endl;
    return 0;
  }

  int nProcEvts = 0; // Number of processed events
  double WeightSumAfter = 0; // Sum of Weights
  double WeightSumBefore = 0; // Sum of Weights

  std::vector<Event*>::iterator evt1 = evtVector.begin();
  if( (*evt1)->type() != PtBalance)std::cout << "Warning: No TwoJetsPtBalanceEvent! Take special care when e.g. analyzing JetTruthEvents"; 
  for(; evt1 != evtVector.end(); ++evt1, ++nProcEvts) {
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
   double MyWeight = LumiWeightsPerTrigger_.at((*it).second).ITweight3BX( nputruth );
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

    std::cout << "  Applied weights for " << nProcEvts << " events in evtVector. \n";
    std::cout << "  Average weight was " << WeightSumAfter/WeightSumBefore << "\n";

    return nProcEvts;
}



int PUTruthReweighting::preprocess(std::vector<Event*>& data,
			     std::vector<Event*>& control1,
			     std::vector<Event*>& control2)
{
  int nProcEvts=0;

  std::cout << "start PUTruthReweighting with:" << std::endl;
  std::cout << "  " << data.size() << " events in data" << std::endl;
  std::cout << "  " << control1.size() << " events in control1" << std::endl;
  std::cout << "  " << control2.size() << " events in control2" << std::endl;


  if(useMCReweightAll_){
    std::cout << "reweighting all event vectors (data,control1, control2) with the same PU-histogram definition:" << std::endl;
    std::cout << "    TruthWeightingDir (target):   " << TruthWeightingDir_ << std::endl;
    std::cout << "    TruthWeightingMCDistribution: " << TruthWeightingMCDistribution_ << std::endl;
    
    if((*data.begin())->type() == GammaJet)nProcEvts += reweightEventVector(data);
    else std::cout << "WARNING: data is not of type 'GammaJet' which corresponds to JetTruthEvents. No reweighting will be applied!" << std::endl;
    if((*control1.begin())->type() == GammaJet)nProcEvts += reweightEventVector(control1);
    else std::cout << "WARNING: control1 is not of type 'GammaJet' which corresponds to JetTruthEvents. No reweighting will be applied!" << std::endl;
    if((*control2.begin())->type() == GammaJet)nProcEvts += reweightEventVector(control2);
    else std::cout << "WARNING: control2 is not of type 'GammaJet' which corresponds to JetTruthEvents. No reweighting will be applied!" << std::endl;
  }


  else{//default: only reweigh control (MC) accordingly
    nProcEvts +=  reweightEventVector(control1);
    nProcEvts +=  reweightEventVector(control2);
  }
  



  std::cout << "end PUTruthReweighting with:" << std::endl;
  std::cout << "  " << data.size() << " events in data" << std::endl;
  std::cout << "  " << control1.size() << " events in control1" << std::endl;
  std::cout << "  " << control2.size() << " events in control2" << std::endl;
  


  return nProcEvts;
}
 
