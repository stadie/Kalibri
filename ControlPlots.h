#ifndef PLOTS_H
#define PLOTS_H

#include <string>
#include <vector>

#include "ControlPlotsConfig.h"
#include "ControlPlotsFunction.h"

class ConfigFile;
class Event;


//!  \brief Create control plots
//!
//!  Creates control plots via the \p makePlots() method from
//!  different \p Event types. The output is in .eps and .root format.
//!  The attributes of the control plots are  specified via 
//!  the configuration file.
//!
//!  $Id: ControlPlots.h,v 1.35 2010/06/09 22:30:53 stadie Exp $
// -------------------------------------------------------------
class ControlPlots {
 public:
  typedef std::vector<Event*>::const_iterator DataIt;

  ControlPlots(const ConfigFile *configFile, const std::vector<std::vector<Event*>* >& samples);
  ~ControlPlots() {};

  void makePlots() const;

 private:
  const ConfigFile *config_; //!< The configuration file
  const std::vector<std::vector<Event*>* >& samples_; //!< The plotted data

  void createJetTruthEventPlots() const;
  void createTwoJetsPtBalanceEventPlots() const;
  ControlPlotsFunction::Function findJetTruthEventFunction(const std::string& varName, ControlPlotsConfig::CorrectionType type = ControlPlotsConfig::Uncorrected) const;
  ControlPlotsFunction::Function findTwoJetsPtBalanceEventFunction(const std::string& varName, ControlPlotsConfig::CorrectionType type = ControlPlotsConfig::Uncorrected) const;
  void setGStyle() const;
};

#endif
