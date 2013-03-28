#ifndef OverlayPlainProfileAnd2DPlots_h
#define OverlayPlainProfileAnd2DPlots_h

#include "BasePlotExtractor.cc"


class OverlayPlainProfileAnd2DPlots {
 public :
  OverlayPlainProfileAnd2DPlots(TString titleLabel="",TString yAxisLabel="");
  ~OverlayPlainProfileAnd2DPlots();
  void addPlots(TString plotsnames="AbsPFFractionVsPt",TString kalibriPlotsShortName="DEFAULT", TString specificPlotName="bla", TString specificSampleName = "PYTHIA", Int_t BinNumber=-1, TString LegendEntry ="def. legentry ");
  void Init();
  void checkConsistency();
  void plotOverlays();
  void plotAll();
  void plotNormalizedComboOverlays();
  void plotNormalizedComboRatiosAll();
  void plotNormalizedDoubleComboOverlays();
  void plotNormalizedDoubleComboRatiosAll();
  void ConfigureSetRangeUser(double ymin, double ymax);

 private:
  int IndexOfSpecificPlot(BasePlotExtractor* BPExtractor, TString specificPlotName);
  int IndexOfSpecificSample(BasePlotExtractor* BPExtractor, TString specificSampleName);
  std::vector<BasePlotExtractor*> BasePlotExtractors_;
  std::vector<TString> SpecificPlotNames_;
  std::vector<TString> SpecificSampleNames_;
  std::vector<TString> LegendEntries_;
  std::vector<Int_t> BinNumbers_;
  std::vector<Int_t> SpecificPlotIndices_;
  std::vector<Int_t> SpecificSampleIndices_;

  std::vector <TH1D*> collectOneDPlots_;
  std::vector < TH1 * > collectRatios_;
  std::vector <TH2D*> collectTwoDPlots_;

  double yminrange_;
  double ymaxrange_;

  int nBins_;
  int nXBins_;
  int nExtrapols_;
  TString titleLabel_;
  TString yAxisLabel_;
};


#endif 


