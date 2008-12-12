//
// Original Author:  Christian Autermann
//         Created:  Wed Jul 18 13:54:50 CEST 2007
// $Id: caliber.h,v 1.35 2008/12/12 17:06:00 stadie Exp $
//
#ifndef caliber_h
#define caliber_h


//C++ libs
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <fstream>
#include <iostream>
#include <cstdlib>
#include <utility> //pair
#include <vector>
#include <map>
#include <sstream>
#include <cmath>
#include <set>


class TParameters;
class TControlPlots;
class TData;
class TMeasurement;

class TCaliber {
public :
  TCaliber(const std::string& f)
    : configfile(f),p(0),plots(0)
 {};
  ~TCaliber(){};

  void Init();
  void Run();
  void Done();
  const char * GetOutputFile(){ return output_file.c_str(); };

protected:  
  //internal functions
  void Run_Lvmini();

  void FlattenSpectra();
  void BalanceSpectra();

private:
  int GetSpectraBin(double m1, double m2, double m3);
  
  //internal variables
  int fit_method, n_gammajet_events, n_dijet_events;
  int n_trijet_events,n_trackcluster_events, n_zjet_events, n_top_events;
  std::string configfile, output_file;              //input/output
  int use_GammaJetTowerMethod,use_DisplayMethod;    //plots
  double Et_cut_on_gamma, Et_cut_on_jet;
  bool useMassConstraintW;
  bool useMassConstraintTop;
  double RelWeight[7];//@@ Replace 7 by something meaningful
 

  std::vector<int> _residualScalingScheme;          // Iteration scheme of scaling of residuals
  double OutlierChi2Cut;                            // Cut on outlier when no scaling is chosen
  int nthreads;
  bool flatten_spectra;
  std::vector<TData*> data;
  
  TParameters * p;    //fit parameters, depend on number of bins & geometry

  TControlPlots * plots;  //the control plots
};

#endif
