#include "TFile.h"
#include <iostream>
#include "TH1D.h"
#include "TString.h"
#include "TCanvas.h"
#include "TAxis.h"

#define UTILS_AS_HEADER_FILE
#include "/afs/naf.desy.de/user/k/kirschen/public/util/utils.h"
#include "/afs/naf.desy.de/user/k/kirschen/public/util/FileOps.h"
#include "/afs/naf.desy.de/user/k/kirschen/public/util/HistOps.h"
#include "/afs/naf.desy.de/user/k/kirschen/public/util/LabelFactory.h"
#include "/afs/naf.desy.de/user/k/kirschen/public/util/StyleSettings.h"


//Script to extract bins with even event numbers from spectrum plots
//ExampleUse: root -b 'extractBinning.C++("KalibriPlots.root","AbsMPFVsRunNumber20/AbsMPFVsRunNumber20_data_RunNumberSpectrum",10)'
void extractBinning(TString kalibriplotsPath, TString SpectrumPlotName, Int_t numberOfDesiredBins=10){

  util::StyleSettings::setStyleJMEPaper();

  TH1D* SpectrumHisto = (TH1D*) util::FileOps::readTH1(kalibriplotsPath,SpectrumPlotName);


  Double_t xq[numberOfDesiredBins+1];  // position where to compute the quantiles in [0,1]
  Double_t yq[numberOfDesiredBins+1];  // array to contain the quantiles
  for (Int_t i=0;i<=numberOfDesiredBins;i++) xq[i] = Float_t(i)/numberOfDesiredBins;
  
  SpectrumHisto->GetQuantiles(numberOfDesiredBins+1,yq,xq);



  std::cout << "Extracting bin edges with similar event numbers for: " << SpectrumHisto->GetXaxis()->GetTitle() << std::endl;
  for (Int_t i=0;i<=numberOfDesiredBins;i++) std::cout << xq[i] << "\t ";// << 
  std::cout <<"\n";
  for (Int_t i=0;i<=numberOfDesiredBins;i++) std::cout << yq[i] << "\t ";// << std::endl;
  std::cout <<"\n";
  SpectrumHisto->Draw();
}
