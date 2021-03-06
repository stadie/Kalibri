//!
//!  \brief Apply event weights
//!    
//!  A weighting factor \f$ w \f$ is calculated per
//!  \f$ \hat{p}_{T} \f$ bin \f$ i \f$ as
//!  \f[
//!   w = N_{ref} \frac{\sigma_{i}}{N_{i}}
//!  \f]
//!  where \f$ \sigma_{i} \f$ is the cross section and
//!  \f$ N_{i} \f$ the number of events in that bin.
//!  \f$ N_{ref} \f$ is the reference to which the weight
//!  is normalized. This is either the luminosity
//!  (\f$ N_{ref} = \mathcal{L} \f$) or else the number of
//!  events in a reference bin \p k  \p
//!  (\f$ N_{ref} = N_{k} / \sigma_{k}  \f$).
//!
//!  The cross sections, number of events and reference have
//!  to be specified in the config file.
//!
//!  \author Matthias Schroeder
//!  \date 2009/07/22
//!  $Id: EventWeightProcessor.h,v 1.6 2010/12/20 11:08:13 stadie Exp $
// -----------------------------------------------------------------

#ifndef EVENT_WEIGHT_PROCESSOR_H
#define EVENT_WEIGHT_PROCESSOR_H

#include <string>
#include <vector>

#include "CalibData.h"
#include "ConfigFile.h"
#include "EventProcessor.h"
#include "Parameters.h"

class EventWeightProcessor : public EventProcessor
{
 public:
  EventWeightProcessor(const std::string& configfile, Parameters* param);
  ~EventWeightProcessor();

protected:
  virtual int preprocess(std::vector<Event*>& data,
			 std::vector<Event*>& control1,
			 std::vector<Event*>& control2);
  virtual int postprocess(std::vector<Event*>& data,
			  std::vector<Event*>& control1,
			  std::vector<Event*>& control2) { return data.size();}
  
  bool applyWeights() const { return weightEvents_; }

 private:
  void calculateWeightsForBins(const std::vector<double>& xSection,
			       const std::vector<int>& nEvents,
			       double lumi, int refPtHatBin);

  bool weightEvents_;              //!< Apply weighting if true
  int type_;                       //!< Type of weighting method
  std::vector<double> minPtHat_;   //!< Minima of \f$ \hat{p}_{T} \f$ bins (type 0)
  std::vector<double> weights_;    //!< Weights of \f$ \hat{p}_{T} \f$ bins (type 0)
  double globalWeight_;            //!< \f$ L\sigma/N_{MC} \f$ (type 1)
  double expo_;                    //!< Exponent for pthat weighting (type 1)
};
#endif
