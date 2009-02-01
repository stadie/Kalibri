//
// Original Author:  Christian Autermann
//         Created:  Wed Jul 18 13:54:50 CEST 2007
// $Id: caliber.h,v 1.39 2009/01/18 13:13:39 stadie Exp $
//
#ifndef caliber_h
#define caliber_h

#include <vector>
#include <string>

class TParameters;
class TControlPlots;
class TData;
class TMeasurement;

class TCaliber {
public :
  TCaliber(const std::string& f)
  : configfile(f),p(0),plots(0),deriv_step(1e-03),eps(1e-02),
  wlf1(1e-04),wlf2(0.9),print_parnderiv(false)
 {};
  ~TCaliber(){};

  void Init();
  void Run();
  void Done();
  const char * GetOutputFile(){ return output_file.c_str(); };

protected:  
  //internal functions
  void Run_Lvmini();

private:
  //internal variables
  int fit_method, n_gammajet_events, n_dijet_events;
  int n_trijet_events,n_trackcluster_events, n_zjet_events, n_top_events;
  std::string configfile, output_file;              //input/output
  //int use_GammaJetTowerMethod,use_DisplayMethod;    //plots
  //bool useMassConstraintW;
  //bool useMassConstraintTop;
 

  std::vector<int> _residualScalingScheme;          // Iteration scheme of scaling of residuals
  double OutlierChi2Cut;                            // Cut on outlier when no scaling is chosen
  int nthreads;
  bool flatten_spectra;
  std::vector<TData*> data;
  
  TParameters * p;    //fit parameters, depend on number of bins & geometry

  TControlPlots * plots;  //the control plots
  // control parameters of fit
  double deriv_step;
  float eps,wlf1,wlf2;
  bool print_parnderiv;
  std::vector<int> globaljetpars;
  std::vector<int> fixedpars;
};

#endif
