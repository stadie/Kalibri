//
// Original Author:  Christian Autermann
//         Created:  Wed Jul 18 13:54:50 CEST 2007
// $Id: caliber.h,v 1.8 2008/02/25 10:07:45 stadie Exp $
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
//Root
#include "TMinuit.h"


//User libs
#include "GammaJetSel.h"
#include "TrackTowerSel.h"
#include "TrackClusterSel.h"
#include "JetJetSel.h"

class TParameters;
class TControlPlots;
class TData;

class TCaliber {
public :
  TCaliber() : plots(0) {};
  ~TCaliber(){};

  void Init(std::string f);
  void Run();
  void Done();
  const char * GetOutputFile(){ return output_file.c_str(); };

protected:  
  //TSelectors:
  GammaJetSel     gammajet;
  TrackTowerSel   tracktower;
  TrackClusterSel trackcluster;
  JetJetSel       jetjet;

  //internal functions
  void Run_Lvmini();
  void Run_Minuit();

  void Run_GammaJet();
  void Run_TrackTower();
  void Run_TrackCluster();
  void Run_JetJet();

private:
  void global_fit(int &npar, double *gin, double &f, 
                         double *allpar, int iflag); 
  double numeric_derivate( void (*func)(int&,double*,double&,double*,int),
                           double * pars, int npar, int index);
  double analytic_derivate( double * pars, int npar, int index);

  //internal variables
  int fit_method, n_gammajet_events, n_tracktower_events, 
      n_trackcluster_events, n_jetjet_events;
  std::string configfile, output_file;              //input/output
  int use_GammaJetTowerMethod,use_DisplayMethod;    //plots
  double Et_cut_on_jet, Et_cut_on_gamma,            //kin. cuts
         Et_cut_on_track, Et_cut_on_tower, Et_cut_on_cluster;

  int    OutlierIterationSteps;                     //outlier rejection
  int nthreads;
  double OutlierChi2Cut, OutlierChi2CutPresel;
  std::vector<TData*> data;
  
  TParameters * p;    //fit parameters, depend on number of bins & geometry

  TControlPlots * plots;  //the control plots
};

#endif
