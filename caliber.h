//
// Original Author:  Christian Autermann
//         Created:  Wed Jul 18 13:54:50 CEST 2007
// $Id: caliber.h,v 1.29 2008/09/29 10:15:16 mschrode Exp $
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

//Root
#include "TMinuit.h"


//User libs
#include "CalibData.h"
#include "GammaJetSel.h"
#include "TrackTowerSel.h"
#include "TrackClusterSel.h"
#include "NJetSel.h"
#include "ZJetSel.h"
#include "TopSel.h"

class TParameters;
class TControlPlots;
class TData;

class TCaliber {
public :
  TCaliber()
    : makeControlPlotsTowers(0), makeControlPlotsGammaJet(0), makeControlPlotsGammaJet2(0),
      makeControlPlotsDiJet(0), makeControlPlotsParScan(0), plots(0)
 {};
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
  NJetSel         dijet;
  NJetSel         trijet;
  ZJetSel         zjet;
  TopSel          top;

  //internal functions
  void Run_Lvmini();
  void Run_Minuit();

  void Run_GammaJet();
  void Run_TrackTower();
  void Run_TrackCluster();
  void Run_NJet(NJetSel & njet, int njet);
  void Run_ZJet();
  void Run_Top();

  void AddTowerConstraint();
  void AddParameterLimits();
  void FlattenSpectra();
  void BalanceSpectra();

private:
  double analytic_derivate( double * pars, int npar, int index);
  int GetSpectraBin(double m1, double m2, double m3);
  
  //internal variables
  int fit_method, n_gammajet_events, n_tracktower_events, 
      n_trackcluster_events, n_dijet_events, n_trijet_events, n_zjet_events,
      n_top_events;
  std::string configfile, output_file;              //input/output
  int use_GammaJetTowerMethod,use_DisplayMethod;    //plots
  double Et_cut_on_jet, Et_cut_on_gamma, Et_cut_nplus1Jet,     //kin. cuts
         Et_cut_on_track, Et_cut_on_tower, Et_cut_on_cluster, Et_cut_on_Z,
         Rel_cut_on_gamma, Rel_cut_on_nJet;
  double RelWeight[7];//@@ Replace 7 by something meaningful
  bool makeControlPlotsTowers;
  bool makeControlPlotsGammaJet;
  bool makeControlPlotsGammaJet2;
  bool makeControlPlotsDiJet;
  bool makeControlPlotsParScan;
  std::set<std::string> mPlottedQuant;

  std::vector<int> _residualScalingScheme;          // Iteration scheme of scaling of residuals
  double OutlierChi2Cut;                            // Cut on outlier when no scaling is chosen
  int nthreads;
  class TowerConstraint {
  public:
    int mineta;
    int maxeta;
    double hadEt;
    double emEt;
    double weight;
    TowerConstraint(int mineta,int maxeta,double hadEt, double emEt, double weight) :
      mineta(mineta),maxeta(maxeta),hadEt(hadEt),emEt(emEt),weight(weight) {}
  };
  class ParameterLimit {
  public:
    int index;
    double min;
    double max;
    double k;
    ParameterLimit(int index, double min, double max, double k) 
      : index(index), min(min), max(max), k(k) {}
  };
  
  std::vector<TowerConstraint> tower_constraints;
  std::vector<ParameterLimit> par_limits;
  bool flatten_spectra;
  std::vector<TData*> data;
  
  TParameters * p;    //fit parameters, depend on number of bins & geometry
  double const (*tower_error_param)(double *const x, TMeasurement *const xorig, double const err);
  double const (*jet_error_param)  (double *const x, TMeasurement *const xorig, double const err);

  TControlPlots * plots;  //the control plots
};

#endif
