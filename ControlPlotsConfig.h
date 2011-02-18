// $Id: ControlPlotsConfig.h,v 1.8 2010/12/13 10:38:28 stadie Exp $

#ifndef CONTROLPLOTS_CONFIG_H
#define CONTROLPLOTS_CONFIG_H

#include <cmath>
#include <map>
#include <string>
#include <vector>

#include "TObject.h"

class ConfigFile;


//!  \brief Reads parameters for profile plots from configuration
//!         file and defines centrally labels and axes titles etc.
//!
//!  This class stores what variables are plotted in the profile
//!  plots, what jet energy correction types are applied and
//!  what profile types are created.
//!
//!  Furthermore, it defines key words which are to be used in the
//!  configuration file to define e.g.
//!  - what variables are plotted
//!  - what jet energy corrections are applied
//!  - what profile types are created.
//!  The class also defines for each key word string representations
//!  which can be used in histogram titles, legends etc.
//!
//!  Possible key words of axis and binning variables are
//!  - Eta
//!  - GenJetPt
//!  - GenJetResponse
//!
//!  Possible key words for jet energy correction types are
//!  - Uncorrected
//!  - Kalibri
//!  - L2L3
//!  
//!  Possible key words for profile types are
//!  - Mean
//!  - StandardDeviation
//!  - GaussFitMean
//!  - GaussFitWidth
//!  - Median
//!  - Chi2
//!  - Probability
//!  - Quantiles
//!
//!  \sa \p ControlPlotsProfile
//!
//!  \author Matthias Schroeder
//!  \date 2009/12/18
//!  $Id: ControlPlotsConfig.h,v 1.8 2010/12/13 10:38:28 stadie Exp $
// ----------------------------------------------------------------   
class ControlPlotsConfig {
 public:
  //! Different jet energy correction types
  enum CorrectionType { Uncorrected, Kalibri, L2L3, L2L3Res, L2L3L4 };
  typedef std::vector<CorrectionType>::const_iterator CorrectionTypeIt;  
  typedef std::pair<int,CorrectionType> InputTag;
  typedef std::vector<InputTag>::const_iterator InputTagsIterator;  

  //! Number of defined profile types
  static const int nProfileTypes = 10;
  //! Different types of profile histograms
  enum ProfileType { Mean=0, StandardDeviation, GaussFitMean, GaussFitWidth, Median, Chi2, Probability, Quantiles, RatioOfMeans, RatioOfGaussFitMeans};
  typedef std::vector<ProfileType>::const_iterator ProfileTypeIt;

  ControlPlotsConfig(const ConfigFile *configFile, const std::string &name);

  //! Returns the profile's name
  std::string name() const { return name_; }

  //! Returns the bin edges of the binning variable
  const std::vector<double> *binEdges() const { return &binEdges_; }
  //! Returns the number of bins
  int nBins() const { return nBins_; }
  //! Returns the minimum of the binning range
  double min() const { return binEdges_.front(); }
  //! Returns the maximum of the binning range
  double max() const { return binEdges_.back(); }
  //! Returns the name of the binning variable
  std::string binVariable() const { return binVar_; }
  //! Returns the name of the bin
  std::string binName(int binIdx) const;
  //! Returns the title of the bin
  std::string binTitle(double min, double max) const;
  //! Returns the name of the cut variable
  std::string cutVariable() const { return cutVar_; }
  //! Returns the minimum of the cut range
  double cutMin() const { return cutEdges_.first; }
  //! Returns the maximum of the cut range
  double cutMax() const { return cutEdges_.second; }
  
  //! Returns the bin edges of the x variable
  const std::vector<double> *xBinEdges() const { return &xBinEdges_; }
  //! Returns the number of x bins
  int nXBins() const { return nXBins_; }
  //! Returns the minimum of the x range
  double xMin() const { return xBinEdges_.front(); }
  //! Returns the maximum of the x range
  double xMax() const { return xBinEdges_.back(); }
  //! Returns the name of the x variable
  std::string xVariable() const { return xVar_; }
  //! Returns the title of the x axis
  std::string xTitle() const { return axisTitle(xVariable()); }
  //! Returns the name of the x bin
  std::string xBinName(int xBinIdx) const;
  //! Returns the title of the x bin
  std::string xBinTitle(int xBinIdx, double binMin, double binMax) const;
  //! Specifies whether there is logarithmic binning of the x axis
  bool logX() const { return logX_; }

  //! Returns the bin edges of the y variable
  const std::vector<double> *yBinEdges() const { return &yBinEdges_; }
  //! Returns the number of y bins
  int nYBins() const { return nYBins_; }
  //! Returns the minimum of the y range
  double yMin() const { return yBinEdges_.front(); }
  //! Returns the maximum of the y range
  double yMax() const { return yBinEdges_.back(); }
  //! Returns the zoomed minimum of the y range
  double yMinZoom(ProfileType profType) const { 
    std::map<ProfileType,double>::const_iterator i = yMinZoom_.find(profType);
    return i->second; 
  }
  //! Returns the zoomed maximum of the y range
  double yMaxZoom(ProfileType profType) const {
    std::map<ProfileType,double>::const_iterator i = yMaxZoom_.find(profType);
    return i->second; 
  }  
  //! Returns the name of the y variable
  std::string yVariable() const { return yVar_; }
  //! Returns the title of the y axis
  std::string yTitle() const { return axisTitle(yVariable()); }
  //! Returns the title of the profile's y axis depending on the \p ProfileType \p type
  std::string yProfileTitle(ProfileType type) const;

  //! Returns the marker and line color for the \p InputTag \p tag 
  int color(const InputTag& tag) const;
  //! Returns the marker style for the \p InputTag \p tag
  int markerStyle(const InputTag& tag) const;
  //! Returns a legend label for the \p InputTag \p tag
  std::string legendLabel(const InputTag& tag) const;

  //! Returns the correction types for which profiles are to be plotted
  const std::vector<InputTag> *inputTags() const { return &inputTags_; }
  //! Returns an iterator to the first correction type to be plotted
  InputTagsIterator inputTagsBegin() const { return inputTags_.begin(); }
  //! Returns an iterator to the last correction type to be plotted
  InputTagsIterator inputTagsEnd() const { return inputTags_.end(); }

  //! Specifies whether the y distributions per x bin are to be drawn
  bool drawDistributions() const { return ! inputTagsDistributions_.empty(); }
  //! Returns the input tags for which y distributions are to be plotted
  const std::vector<InputTag> *distributionInputTags() const { return &inputTagsDistributions_; }
  //! Returns an iterator to the first correction type to be plotted
  InputTagsIterator distributionInputTagsBegin() const { return inputTagsDistributions_.begin(); }
  //! Returns an iterator to the last correction type to be plotted
  InputTagsIterator distributionInputTagsEnd() const { return inputTagsDistributions_.end(); }

  //! Returns the \p CorrectionType for a given correction type name
  CorrectionType correctionType(const std::string &typeName) const;
   //! Returns the name of a given \p CorrectionType
  std::string correctionTypeName(CorrectionType corrType) const; 
  //! Returns the name of a given \p Sample
  std::string sampleName(int sample) const;
  //! Returns the name of a given \p InputTag
  std::string inputTagName(const InputTag& tag) const;

  //! Returns the profile types which are to be drawn
  const std::vector<ProfileType> *profileTypes() const { return &profTypes_; }
  //! Returns an iterator to the first profile type to be drawn
  ProfileTypeIt profileTypesBegin() const { return profTypes_.begin(); }
  //! Returns an iterator to the last profile type to be drawn
  ProfileTypeIt profileTypesEnd() const { return profTypes_.end(); }
  //! Returns the \p ProfileType for a given profile type name
  ProfileType profileType(const std::string &typeName) const;
  //! Returns the profile type name for a given \p ProfileType
  std::string profileTypeName(ProfileType profType) const;

  //! Returns the name of the directory in which the control plots are stored
  std::string outDirName() const { return outDirName_; }
  //! Returns the file ending ("eps") of the control plots
  std::string outFileType() const { return outFileType_; }
  //! Specifies whether the plots are only saved in a root-file.
  bool outOnlyRoot() const { return outOnlyRoot_; }
  //! Writes a \p obj to ROOT file 
  void toRootFile(TObject *obj) const;


 private:
  const ConfigFile *config_;
  const std::string name_;

  std::string binVar_;
  std::vector<double> binEdges_;
  int nBins_;

  std::string xVar_;
  std::vector<double> xBinEdges_;
  int nXBins_;
  bool logX_;

  std::string cutVar_;
  std::pair<double,double> cutEdges_;

  std::string yVar_;
  std::vector<double> yBinEdges_;
  int nYBins_;
  std::map<ProfileType,double> yMinZoom_;
  std::map<ProfileType,double> yMaxZoom_;

  std::map<int,std::string> sampleNames_;
  std::vector<InputTag> inputTags_;
  std::vector<InputTag> inputTagsDistributions_;
  std::vector<ProfileType> profTypes_;

  std::string outDirName_;
  std::string outFileType_;
  bool outOnlyRoot_;

  std::map<InputTag,int> colors_;
  std::map<InputTag,int> markerStyles_;
  std::map<InputTag,std::string> legendLabels_;

  //! Read parameter values from configuration file
  void init();
  //! Create a binning with equal bin sizes in logarithmic scale
  bool equidistLogBins(std::vector<double>& bins, int nBins, double first, double last) const;
  //! Returns the histogram axis title for a variable \p varName
  std::string axisTitle(const std::string &varName) const;
  //! Returns the unit for a variable \p varName
  std::string unitTitle(const std::string &varName) const;
  //! Returns the variable name as it appears in axis titles etc for \p varName
  std::string varTitle(const std::string &varName) const;
  
  double round(double x) const { return (x > 0.5) ? ceil(x) : floor(x); }
  template <class T> std::string toString(const T& t) const;
};
#endif
