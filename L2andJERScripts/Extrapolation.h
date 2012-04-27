#ifndef Extrapolation_h
#define Extrapolation_h

#include "BasePlotExtractor.cc"
#include "TGraphErrors.h"

typedef std::vector<TGraphErrors*> TGErrvec_t;
typedef std::vector<TGraphErrors*> VecOfTGErrvec_t;


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
class Extrapolation :public BasePlotExtractor{
 public :
  Extrapolation(TString plotsnames="ResolutionVsPt",TString kalibriPlotsPath="KalibriPlots.root");
  void Plot();
 private:
  void extrapolInit();
  void createPtRelExtrapol();
  std::vector<std::string> cutNames_;
  std::vector<double> cutNumbers_;
  //  VecOfVecOfTH1vec_t AllPlots_;
  //  void doExtrapol();  
  class ExtrapolateBin {
  public:
    ExtrapolateBin(Extrapolation* Outer);
    void addMCHisto(TH1D* MCHisto);
    void addDataHisto(TH1D* DataHisto);
    void createTGraphErrors(Int_t xBin_i);
    void plotExtrapol(Int_t xBin_i, Int_t bin_i);
    void produceExtrapolatedRes();
    TH1D* ExtrapolatedResMC() const {return ExtrapolatedResMC_;}
    TH1D* ExtrapolatedResData() const {return ExtrapolatedResData_;}
    std::pair <float,float> determineMinMax(TGraphErrors* graph);

  private:
    std::vector<TH1D*> MCHistos_;
    std::vector<TH1D*> DataHistos_;
    TGErrvec_t MCExtrapols_;
    TGErrvec_t DataExtrapols_;
    //    std::vector<double> cutNumbers_;
    //    std::vector<double> cutNumbers_;
    Extrapolation* Outer_;
    TH1D* ExtrapolatedResMC_;
    TH1D* ExtrapolatedResData_;
  };

};

#endif 


