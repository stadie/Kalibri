#ifndef TControlPlots_h
#define TControlPlots_h

//C++ libs
#include <string>
#include <cmath>
//User
#include "Parameters.h"
#include "TF1.h"

#include <iostream>

class TData;
class TH2F;
class TH1F;
class TControlPlots {
public:
  //TControlPlots(){this->ReadConfigFile("config/calibration.cfg");};
  TControlPlots(std::string config, std::vector<TData*> * d, TParameters * pars)
                {this->ReadConfigFile(config); data=d; p=pars; iplot=0;};
  //TControlPlots(std::string config, std::vector<TData*> * d, TStepEfracParameters * pars)
  //              {this->ReadConfigFile(config); data=d; p=pars;};
  ~TControlPlots(){};
  
  void GammaJetControlPlots();
  void GammaJetControlPlotsJetBin();
  void GammaJetControlPlotsJetJEC();
  void DiJetControlPlots();
  void TrackTowerControlPlots();
  void TrackClusterControlPlots();
  void FitControlPlots();
  void OutlierControlPlots() const;

private:  
  void ReadConfigFile(std::string config);
  void Fit2D(TH2F* hist, TH1F* hresults[8], TH1F* gaussplots[4], TF1* gf[4]);

  std::vector<TData*> * data;
  TParameters * p;
  int iplot;

  bool _doPlots;
  //... bins, pt ranges, etc...
  
};

#endif
