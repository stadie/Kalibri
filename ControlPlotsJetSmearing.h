// $Id: $

#ifndef JS_CONTROLPLOTS_JETSMEARING_H
#define JS_CONTROLPLOTS_JETSMEARING_H

#include <string>

#include "CalibData.h"
#include "ConfigFile.h"
#include "Parameters.h"


//!  \brief Generates validation plots for jet-smearing method
//!  \author Matthias Schroeder
//!  \date Thu May  7 11:30:28 CEST 2009 
//!  $Id: $
// --------------------------------------------------
class ControlPlotsJetSmearing {
 public:
  ControlPlotsJetSmearing(const std::string& configfile,const std::vector<TData*> * data, TParameters * param);
  ~ControlPlotsJetSmearing() {};
  
  void PlotDijets() const;
  void PlotResponse() const;
  void PlotParameterScan(const std::vector<unsigned int>& pars) const;
  void SetBinningResp(int nbins, double min, double max) { mRespNBins = nbins; mRespMin = min; mRespMax = max;}
  void SetBinningDiJet(int nbins, double min, double max) { mDijetNBins = nbins; mDijetMin = min; mDijetMax = max;}
  void SetBinningPhotonJet(int nbins, double min, double max) { mPhotonJetNBins = nbins; mPhotonJetMin = min; mPhotonJetMax = max;}


 private:
  typedef std::vector<TData*>::const_iterator DataIt;

  const std::vector<TData*> * mData;
  const ConfigFile          * mConfig;
  mutable TParameters       * mParam;
  
  int          mDijetNBins;
  double       mDijetMin;
  double       mDijetMax;
  int          mPhotonJetNBins;
  double       mPhotonJetMin;
  double       mPhotonJetMax;
  int          mRespNBins;
  double       mRespMin;
  double       mRespMax;
  std::string  mDir;
  
  void SetGStyle() const;
};
#endif
