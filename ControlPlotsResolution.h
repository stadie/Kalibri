//
// $Id: ControlPlotsResolution.h,v 1.12 2010/09/22 13:29:44 mschrode Exp $
//
#ifndef CONTROLPLOTS_RESOLUTION_H
#define CONTROLPLOTS_RESOLUTION_H

#include <string>
#include <vector>

#include "CalibData.h"
#include "ConfigFile.h"
#include "Parameters.h"

#include "TString.h"

class Jet;
class TCanvas;
class TF1;
class TH1;
class TLegend;
class TObject;
class TPaveText;
class TPostScript;
class TRandom3;


//!  \brief Generates validation plots for resolution measurement
//!  \author Matthias Schroeder
//!  \date Thu May  7 11:30:28 CEST 2009 
//!  $Id: ControlPlotsResolution.h,v 1.12 2010/09/22 13:29:44 mschrode Exp $
// --------------------------------------------------
class ControlPlotsResolution {
 public:
  ControlPlotsResolution(const std::string& configfile,const std::vector<Event*> * data, Parameters * param, const std::string &outDir = "./controlPlots");
  ~ControlPlotsResolution();

  void makePlots() const;
  void setBinningResp(int nbins, double min, double max) { respNBins_ = nbins; respMin_ = min; respMax_ = max;}


 private:
  static double spectrum(double *x, double *par);
  static double gaussian(double *x, double *par);

  typedef std::vector<Event*>::const_iterator DataIt;

  const std::vector<Event*> * data_;   //!< The data which is plotted
  const ConfigFile          * config_; //!< The configuration file
  mutable Parameters        * param_;  //!< The parametrization
  
  bool saveAsEps_;
  std::string outNamePrefix_;
  int          respNBins_;             //!< Number of bins in response control plots \p plotResponse()
  double       respMin_;               //!< Minimum of response control plots \p plotResponse()
  double       respMax_;               //!< Maximum of response control plots \p plotResponse()
  std::string  dir_;                   //!< Directory in which the control plots are written
  TRandom3 *rand_;
  std::string parClass_;
  std::vector<double> scale_;
  std::vector<double> truthPar_;
  TF1* truthRes_;
  std::string ptBinningVar_;
  std::vector<double> ptBinEdges_;
  double etaMin_;
  double etaMax_;
  std::vector<TString> titleBins_;
  TString title_;

  void plotResponse() const;
  void plotParameters() const;
  //! Plots the negative log-likelihood for different parameter values
  void plotParameterScan() const;
  //! Plots the distributions of the probability density of
  //! each event before and after the fit
  void plotLogP() const;
  void plotMeanResponseAndResolution() const;
  void plotParallelComponents() const;
  void plotTails() const;

  std::vector<double> findBinEdgesForCombinedSpectrum() const;
  void interpolateTails(TH1* hTails, TH1* &hInter) const;
  void subtractGaussian(const TH1* hResp, TH1* &hTail, TH1* &hTailClean, TF1* &gauss) const;

  double gaussianWidth(double pt) const;
  double gaussianWidthError(double pt) const;
  double gaussianWidthTruth(double pt) const;

  double scale(unsigned int i) const {
    return i < scale_.size() ? scale_[i] : 1.;
  }
  bool equidistLogBins(std::vector<double>& bins, int nBins, double first, double last) const;

  TLegend *createLegend(int nEntries, double width = 1., double lineHgt = -1., double yOffset = 0.) const;
  TPaveText *createPaveText(int nEntries, double width = 1., double lineHgt = -1.) const;
  void drawPSPage(TPostScript * ps, TCanvas * can, TObject * obj, const std::string &option = "", bool log = false) const;
  void drawPSPage(TPostScript * ps, TCanvas * can, std::vector<TObject*> objs, const std::string &option = "", bool log = false) const;
  void findYRange(const TH1 * h, double& min, double& max) const;
  double lineHeight() const { return 0.06; }
  void normHist(TH1 * h, std::string option = "") const;
  void normHist(TH1 *h, double min, double max, std::string option = "") const;
  void setGStyle() const;
  void setYRange(TH1 * h, double c1 = 0.9, double c2 = 1.1, double minLimit = 0.) const;
  template <class T> std::string toString(const T& t) const;
  template <class T> TString toTString(const T& t) const { return toString(t).c_str(); }
  int color(int i) const;
  void setAxisTitles(TH1 *h, const std::string &xTitle, const std::string &xUnit, const std::string &yTitle, bool norm = false) const;
  void setColor(TH1 *h, int color) const;
};
#endif
