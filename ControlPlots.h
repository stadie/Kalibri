#ifndef TControlPlots_h
#define TControlPlots_h

//C++ libs
#include <string>
#include <cmath>
//User
#include "Parameters.h"

#include <iostream>

class TData;
  
class TControlPlots {
public:
  //TControlPlots(){this->ReadConfigFile("config/calibration.cfg");};
  TControlPlots(std::string config, std::vector<TData*> * d, TStepParameters * pars)
                {this->ReadConfigFile(config); data=d; p=pars;};
  ~TControlPlots(){};
  
  void GammaJetControlPlots();
  void GammaJetControlPlotsJetBin();
  void TrackTowerControlPlots();
  void TrackClusterControlPlots();
  void FitControlPlots();

  TStepParameters * p;

private:  
  void ReadConfigFile(std::string config);

  std::vector<TData*> * data;
 
  bool _doPlots;
  //... bins, pt ranges, etc...
  
};

#endif
