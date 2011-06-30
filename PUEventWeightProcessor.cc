#include "PUEventWeightProcessor.h"

#include <iostream>

#include "TFile.h"
#include "TH1.h"

#include "ConfigFile.h"
#include "CalibData.h"


// -----------------------------------------------------------------
PUEventWeightProcessor::PUEventWeightProcessor(const std::string& configfile, Parameters* param)
  : EventProcessor("PU weighting",configfile,param), weightEvents_(false) {

  ConfigFile config(configfile.c_str());

  std::string histFileName  = config.read<string>(name()+" histogram","");

  if( histFileName != "" ) {
    weightEvents_ = true;
    std::cout << "Weighting events to specified PU scenario" << std::endl;
    std::cout << "  Reading PU scenario from '" << histFileName << "'"  << std::endl;
    
    TFile file(histFileName.c_str(),"READ");
    TH1 *h = 0;
    file.GetObject("pileup",h);
    if( h ) {
      h->SetDirectory(0);
    } else {
      std::cerr << "ERROR in PUEventWeightProcessor: Histogram 'pileup' does not exist in file '" << histFileName << "'\n.";
      std::cerr << "See https://twiki.cern.ch/twiki/bin/view/CMS/HamburgWikiAnalysisCalibration#Pile_Up_Reweighting for available input distributions." << std::endl;
      exit(1);
    }
    file.Close();

    std::cout << "  Computing weights" << std::endl;
    weights_ = generate_flat10_weights(h);

    delete h;
  }
}



// -----------------------------------------------------------------
int PUEventWeightProcessor::preprocess(std::vector<Event*>& data,
				       std::vector<Event*>& control1,
				       std::vector<Event*>& control2) {

  if( !weightEvents_ ) return data.size();
  
  int nProcEvts = 0; // Number of processed events
  std::vector<Event*>::iterator evt = data.begin();
  for(; evt != data.end(); ++evt, ++nProcEvts) {
    if( (*evt)->type() == ParLimit ) continue;
    short int npu = (*evt)->nPU();
    if( npu < static_cast<short int>(weights_.size()) ) {
      (*evt)->setWeight( weights_.at(npu) * ((*evt)->weight()) );
    } else {
      std::cerr << "WARNING in PUEventWeightProcessor::preprocess: Number of PU vertices = " << npu << " out of histogram binning." << std::endl;
    }
  }
  std::cout << "  Applied weights for " << nProcEvts << " events\n";
  
  return nProcEvts;
}



// Generate weights for Flat10 PU scenario for given
// data PU distribution
// Code from: https://twiki.cern.ch/twiki/bin/viewauth/CMS/PileupReweighting
// --------------------------------------------------
std::vector<double> PUEventWeightProcessor::generate_flat10_weights(const TH1* data_npu_estimated) const {
  // see SimGeneral/MixingModule/python/mix_E7TeV_FlatDist10_2011EarlyData_inTimeOnly_cfi.py; copy and paste from there:
  const double npu_probs[25] = {0.0698146584, 0.0698146584, 0.0698146584,0.0698146584,0.0698146584,0.0698146584,0.0698146584,0.0698146584,0.0698146584,0.0698146584,0.0698146584 /* <-- 10*/,0.0630151648,0.0526654164,0.0402754482,0.0292988928,0.0194384503,0.0122016783,0.007207042,0.004003637,0.0020278322,0.0010739954,0.0004595759,0.0002229748,0.0001028162,4.58337152809607E-05 /* <-- 24 */};
  std::vector<double> result(25);
  double s = 0.0;
  for(int npu=0; npu<25; ++npu) {
    double npu_estimated = data_npu_estimated->GetBinContent(data_npu_estimated->GetXaxis()->FindBin(npu));                              
    result[npu] = npu_estimated / npu_probs[npu];
    s += npu_estimated;
  }
  // normalize weights such that the total sum of weights over thw whole sample is 1.0, i.e., sum_i  result[i] * npu_probs[i] should be 1.0 (!)
  for(int npu=0; npu<25; ++npu) {
    result[npu] /= s;
  }
  return result;
}
