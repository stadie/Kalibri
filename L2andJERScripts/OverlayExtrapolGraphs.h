#ifndef OverlayExtrapolGraphs_h
#define OverlayExtrapolGraphs_h

#include "Extrapolation.cc"

class OverlayExtrapolGraphs {
 public :
  void addPlots(TString plotsnames="AbsPFFractionVsPt",TString kalibriPlotsShortName="DEFAULT");
  void readInGraphs();
  void checkConsistency();
  void plotOverlays();
  void plotAll();

 private:
  std::vector<Extrapolation> Extrapolations_;

  std::vector<VecOfTGErrvec_t> allCollectedGraphsXBins_; //contains vectors for each xBin with vectors of tgraphs containing all extrapolations for one bin.
  //  VecOfTGErrvec_t allCollectedGraphs_; //contains vectors of tgraphs containing all extrapolations for one bin.
  int nBins_;
  int nXBins_;
  int nExtrapols_;

};


#endif 


