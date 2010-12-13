// $Id: ControlPlotsProfile.h,v 1.3 2010/11/19 10:11:40 stadie Exp $

#ifndef CONTROL_PLOTS_PROFILE_H
#define CONTROL_PLOTS_PROFILE_H


#include <map>
#include <string>
#include <vector>

#include "TH1D.h"
#include "TH2D.h"
#include "TLegend.h"
#include "TLine.h"

#include "ControlPlotsConfig.h"


class ControlPlotsFunction;
class Event;
class TGraph;

//! \brief Creates different profile plots from 2D histograms
//!
//! This class creates 2D histograms y vs x in different bins
//! of a binning variable for different jet energy correction
//! types. The values of y, x, and the binning variable are
//! found for each event from a \p ControlPlotsFunction.
//!
//! Different profiles along x are created from the 2D histograms
//! such as the mean value or the standard deviation.
//!
//! \author Matthias Schroeder
//! \date 2009/12/18
//! $Id: ControlPlotsProfile.h,v 1.3 2010/11/19 10:11:40 stadie Exp $
// ----------------------------------------------------------------   
class ControlPlotsProfile {
 public:
  //! Constructor
  ControlPlotsProfile(const ControlPlotsConfig *config, const ControlPlotsFunction *function);
  //! Destructor
  ~ControlPlotsProfile();

  //! Draws the histograms and profiles
  void draw();
  //! Fills the 2D histograms
  void fill(const Event *evt, int id);
  //! Fits the profile histograms from the 2D histograms in all bins
  void fitProfiles();

 private:
  //! \brief Helper class for \p ControlPlotsProfile
  //!
  //! Contains and fills the 2D histograms, profiles, and y distributions
  //! in one bin for the different correction and profile types.
  //! 
  //! \author Matthias Schroeder
  //! \date 2009/12/18
  //! $Id: ControlPlotsProfile.h,v 1.3 2010/11/19 10:11:40 stadie Exp $
  // ----------------------------------------------------------------   
  class Bin {
  public:
    //! Constructor
    Bin(int binIdx, double min, double max, const ControlPlotsConfig *config);
    //! Destructor
    ~Bin();

    //! Returns the bin's minimum
    double min() const { return min_; }
    //! Returns the bin's maximum
    double max() const { return max_; }

    //! Returns the 2D histogram of the correction type \p corrType
    TH2D *hYvsX(const ControlPlotsConfig::InputTag& tag);
    //! Returns the profile of type \p profType for the correction type \p corrType
    TH1D *hXProfile(const ControlPlotsConfig::InputTag& tag, ControlPlotsConfig::ProfileType profType);
    //! Returns the y distribution of the nth x bin for the correction type \p corrType
    TH1D *hYDistribution(int n, const ControlPlotsConfig::InputTag& tag);
    //! Returns the number of y distributions i.e. the number of x bins
    int nDistributions() const { return (hYDistributions_.begin())->second.size(); }
    //! Fills the 2D histogram y vs x with weight w for the correction type \p corrType
    int fill(double x, double y, double w,const ControlPlotsConfig::InputTag& tag);
    //! Fits the profiles from the 2D histograms
    int fitProfiles();

    //! Returns a horizontal line to be drawn into the profile 
    TLine *createHorizontalLine() const;
    //! Returns a legend to be drawn into the profile
    TLegend *createLegend();
    //! Returns the file name of the 2D histograms
    std::string hist2DFileName(const ControlPlotsConfig::InputTag& tag) const;
    //! Returns the file name of the profile histograms
    std::string profileFileName(ControlPlotsConfig::ProfileType type) const;
    //! Returns the file name of the y distributions
    std::string distributionFileName(int xBin, const ControlPlotsConfig::InputTag& tag) const;

    int id() const { return idx_;}
  private:
    const int idx_;
    const double min_;
    const double max_;
    const ControlPlotsConfig *config_;

    std::map< ControlPlotsConfig::InputTag, TH2D* > hYvxX_;
    std::map< ControlPlotsConfig::InputTag, std::map< ControlPlotsConfig::ProfileType, TH1D* > > hXProfile_;
    std::map< ControlPlotsConfig::InputTag, std::vector< TH1D* > > hYDistributions_;
 
    //! Counter of the number of calls of \p fitProfiles()
    int nCallsFitProfiles_;
  };

  const ControlPlotsConfig *config_;
  const ControlPlotsFunction *function_;


  std::vector<Bin*> bins_;
  TH1D *hXSpectrum_[3];

  //! Find the index of the \p Bin in which the event \p evt falls
  int findBin(const Event *evt) const;

  static TGraph* makeGraph(TH1* h);
};
#endif
