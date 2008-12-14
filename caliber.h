//
// Original Author:  Christian Autermann
//         Created:  Wed Jul 18 13:54:50 CEST 2007
// $Id: caliber.h,v 1.36 2008/12/12 17:52:14 stadie Exp $
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
};

#endif
