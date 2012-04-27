//
// $Id: PUEventWeightProcessor.h,v 1.2 2011/11/25 07:14:32 kirschen Exp $
//
#ifndef PU_EVENT_WEIGHT_PROCESSOR_H
#define PU_EVENT_WEIGHT_PROCESSOR_H

#include <vector>
#include <string>

#include "EventProcessor.h"

class Event;
class Parameters;
class TH1;



//! \brief Weight MC events to a PU scenario by number of added PU interactions
//! 
//! The weighting function is specified as a histogram from
//! a ROOT file.
//!
//! \see https://twiki.cern.ch/twiki/bin/viewauth/CMS/PileupReweighting
// -----------------------------------------------------------------
class PUEventWeightProcessor : public EventProcessor {
public:
  PUEventWeightProcessor(const std::string& configfile, Parameters* param);

  int preprocess(std::vector<Event*>& data,
		 std::vector<Event*>& control1,
		 std::vector<Event*>& control2);
  int postprocess(std::vector<Event*>& data,
		  std::vector<Event*>& control1,
		  std::vector<Event*>& control2) { return 0; }


private:
  bool weightEvents_;
  std::vector<double> weights_;

  std::vector<double> generate_flat10_weights(const TH1* data_npu_estimated) const;
  std::vector<double> generate_fall11_weights(const TH1* data_npu_estimated) const;
  std::vector<double> generate_summer12_weights(const TH1* data_npu_estimated) const;
};

#endif

