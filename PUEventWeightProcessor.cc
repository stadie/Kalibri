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
  std::string PU_mixing_era  = config.read<string>(name()+" era","Flat10");

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
    if(PU_mixing_era=="Fall11"){
      weights_ = generate_fall11_weights(h);
      std::cout << "using Fall11 distribution" << std::endl;
    }
    else {
      weights_ = generate_flat10_weights(h);
      std::cout << "using Flat10 distribution" << std::endl;
    }
    delete h;
  }
}



// -----------------------------------------------------------------
int PUEventWeightProcessor::preprocess(std::vector<Event*>& data,
				       std::vector<Event*>& control1,
				       std::vector<Event*>& control2) {

  if( !weightEvents_ ) return data.size();
  
  int nProcEvts = 0; // Number of processed events
  std::vector<Event*>::iterator evt1 = control1.begin();
  for(; evt1 != control1.end(); ++evt1, ++nProcEvts) {
    if( (*evt1)->type() == ParLimit ) continue;
    short int npu = (*evt1)->nPU();
    if( npu < static_cast<short int>(weights_.size()) ) {
      (*evt1)->setWeight( weights_.at(npu) * ((*evt1)->weight()) );
    } else {
      std::cerr << "WARNING in PUEventWeightProcessor::preprocess: Number of PU vertices = " << npu << " out of histogram binning." << std::endl;
    }
  }
  std::vector<Event*>::iterator evt2 = control2.begin();
  for(; evt2 != control2.end(); ++evt2, ++nProcEvts) {
    if( (*evt2)->type() == ParLimit ) continue;
    short int npu = (*evt2)->nPU();
    if( npu < static_cast<short int>(weights_.size()) ) {
      (*evt2)->setWeight( weights_.at(npu) * ((*evt2)->weight()) );
    } else {
      std::cerr << "WARNING in PUEventWeightProcessor::preprocess: Number of PU vertices = " << npu << " out of histogram binning." << std::endl;
    }
  }
  std::cout << "  Applied weights for " << nProcEvts << " events in control1 and control2. \n";
  
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
  // normalize weights such that the total sum of weights over the whole sample is 1.0, i.e., sum_i  result[i] * npu_probs[i] should be 1.0 (!)
  for(int npu=0; npu<25; ++npu) {
    result[npu] /= s;
  }
  return result;
}



// Generate weights for Fall11 Distribution PU scenario for given
// data PU distribution
// Code from: https://twiki.cern.ch/twiki/bin/viewauth/CMS/PileupReweighting
// --------------------------------------------------
std::vector<double> PUEventWeightProcessor::generate_fall11_weights(const TH1* data_npu_estimated) const {
  //Distribution extracted from Fall11_PYTHIA_S6 MC
Double_t npu_probs_Fall2011[50] = {
  0.00867359,0.0188749,0.0311076,0.0424069,0.0527288,0.0586238,0.0623814,0.0618389,0.0605955,0.0569491,0.0529175,0.0493315,0.0454852,0.0423341,0.0389899,0.0362847,0.0334895,0.0309533,0.0282723,0.0258533,0.0233813,0.0209824,0.0185732,0.0163106,0.0142085,0.0122677,0.0105115,0.00880698,0.00739917,0.00610915,0.00499312,0.0040509,0.00325658,0.00258698,0.00202705,0.00158217,0.00123929,0.000932184,0.000730729,0.000534201,0.000398986,0.000301543,0.000216509,0.000162587,0.000110855,8.18409e-05,5.9305e-05,4.21522e-05,3.03824e-05,2.02549e-05};
  // see https://twiki.cern.ch/twiki/bin/viewauth/CMS/PileupInformation and mix_E7TeV_Fall2011_Reprocess_50ns_PoissonOOTPU_cfi.py copy and paste from there:
//  Double_t npu_probs_Fall2011[50] = { 0.003388501, 0.010357558, 0.024724258, 0.042348605, 0.058279812, 0.068851751, 0.072914824, 0.071579609, 0.066811668, 0.060672356, 0.054528356, 0.04919354, 0.044886042, 0.041341896, 0.0384679, 0.035871463, 0.03341952, 0.030915649, 0.028395374, 0.025798107, 0.023237445, 0.020602754, 0.0180688, 0.015559693, 0.013211063, 0.010964293, 0.008920993, 0.007080504, 0.005499239, 0.004187022, 0.003096474, 0.002237361, 0.001566428, 0.001074149, 0.000721755, 0.000470838, 0.00030268, 0.000184665, 0.000112883, 6.74043E-05, 3.82178E-05, 2.22847E-05, 1.20933E-05, 6.96173E-06, 3.4689E-06, 1.96172E-06, 8.49283E-07, 5.02393E-07, 2.15311E-07, 9.56938E-08};

  std::vector<double> result(50);
  double s = 0.0;
  for(int npu=0; npu<50; ++npu) {
    double npu_estimated = data_npu_estimated->GetBinContent(data_npu_estimated->GetXaxis()->FindBin(npu));                              
    result[npu] = npu_estimated / npu_probs_Fall2011[npu];
    s += npu_estimated;
  }
  // normalize weights such that the total sum of weights over the whole sample is 1.0, i.e., sum_i  result[i] * npu_probs_Fall2011[i] should be 1.0 (!)
  for(int npu=0; npu<50; ++npu) {
    result[npu] /= s;
  }
  return result;
}
