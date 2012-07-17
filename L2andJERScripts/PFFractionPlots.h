#ifndef PFFractionPlots_h
#define PFFractionPlots_h

#include "BasePlotExtractor.cc"

//!  \brief Reads (PF-)fraction plots and produces stacked
//!         histograms of the input
//!
//!  Makes use of Kalibri classes used for ControlPlots to get binning
//!  names and other stuff right.
//!  Default is to 
//!
//!
//!  
//!
//!  \author Henning Kirschenmann
//!  \date 2012/03/07
// ----------------------------------------------------------------   
class PFFractionPlots : public BasePlotExtractor{
 public :
  PFFractionPlots(TString plotsnames="AbsPFFractionVsPt",TString kalibriPlotsShortName="KalibriPlots.root");
  void Plot();

};

#endif 


