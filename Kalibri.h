//  $Id: Kalibri.h,v 1.12 2011/06/03 15:53:55 stadie Exp $

//!  \mainpage
//!
//!  \image html kalibriLogoSmall.jpg
//!  Package for data driven calibration using an unbinned fit (see also the related
//!  <A HREF="https://twiki.cern.ch/twiki/bin/view/CMS/HamburgWikiAnalysisCalibration">
//!  Twiki Page</A>).
//!
//!  \section label_sec_src Source Code
//!  The source code can be found
//!  <A HREF="http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/UserCode/Bromo/Calibration/CalibCore/">here</A>.
//!
//!  \section label_sec_workflow Workflow
//!  \image html kalibri_workflow.png
//!  (Graphic in <A HREF="../graphic/kalibri.eps">eps</A> format.)
//!
//!  \section label_sec_authors Authors
//!  - Christian Autermann
//!  - Ulla Gebbert
//!  - Robert Klanner
//!  - Bj&ouml;rn Kolodzey
//!  - Sebastian Naumann-Emme
//!  - Christian Sander
//!  - Matthias Schr&ouml;der
//!  - Torben Schum
//!  - Hartmut Stadie
//!  - Jan Thomsen
//!  - Roger Wolf

 

#ifndef caliber_h
#define caliber_h

#include <vector>
#include <string>

#include "include/lbfgs.h"
#include "Math/IFunction.h"
#include <iostream>

class Parameters;
class Controlplots;
class Event;
class Measurement;
class ComputeThread;

//!  \brief Main program
//!  \note  For profiling:
//!         To prevent gprof from missing the threads: 
//!         wget http://sam.zoy.org/writings/programming/gprof-helper.c
//!         gcc -shared -fPIC gprof-helper.c -o gprof-helper.so -lpthread -ldl 
//!         LD_PRELOAD=./gprof-helper.so ./junk
//!  \authors Christian Autermann, Hartmut Stadie, Matthias Schroeder
//!  \date Wed Jul 18 13:54:50 CEST 2007
//!  $Id: Kalibri.h,v 1.12 2011/06/03 15:53:55 stadie Exp $
// -----------------------------------------------------------------
class Kalibri {
public :
  //!  \brief Constructor
  //!  \param f Name of the configuration file
  // -----------------------------------------------------------------
  Kalibri(const std::string& f)
    : configFile_(f), par_(0), fitMethod_(1), minName_("Minuit2"),
      algoName_(""), nThreads_(1), nGammajetEvents_(0),
      nDijetEvents_(0), nTrijetEvents_(0), nTrackClusterEvents_(0),
      nZjetEvents_(0), nTopEvents_(0), printParNDeriv_(false), 
      derivStep_(1e-03), mvec_(6), nIter_(100), eps_(1e-02), wlf1_(1e-04),
      wlf2_(0.9), calcCov_(false), epsilon_(0), temp_derivative1_(0),
      temp_derivative2_(0),threads_(0)
  {};

  ~Kalibri(){};

  void init();   //!< Read parameters from configfile, read data
  void run();    //!< Run the fit
  void done();   //!< Make control plots, clean up
  const char * getOutputFile() { return outputFile_.c_str(); }; //!< Get the ouputfile name

protected:  
  //internal functions
  void run_Lvmini();  //!< Run the fit
  void run_lbfgs();   //!< Run the fit
  void run_Minimizer();   //!< Run the fit
  void stressTest(); //!< Check the math

  double eval(const double *x, double *f1 = 0, double *f2=0);
  class Funct : public ROOT::Math::IGradientFunctionMultiDim {
  public:
    Funct(Kalibri *k) : ROOT::Math::IGradientFunctionMultiDim(), k_(k) {}
    
    void FdF(const double* x, double& f, double* df) const {
      f = k_->eval(x,df,0);
    }
    ROOT::Math::IBaseFunctionMultiDim* Clone() const {
      return new Funct(k_);
    }
    void Gradient(const double *x, double * grad) const { 
      //std::cerr << "Kalibri::Funct::Gradient called. though it is inefficient!\n";
      k_->eval(x,grad,0);
    }
    unsigned int NDim() const; 
    
  private:
    Kalibri* const k_; 
    double DoEval(const double*x) const {
      return k_->eval(x,0,0); 
    }
    double DoDerivative(const double* x, unsigned int i) const {
      std::cerr << "Kalibri::Funct::DoDerivative called, but very slow\n";
      double *df = new double[NDim()];
      k_->eval(x,df,0);
      double res = df[i];
      delete [] df;
      return res;
    }
  };
  


  static lbfgsfloatval_t lbfgs_evaluate(void *instance, 
					const lbfgsfloatval_t *x,
					lbfgsfloatval_t *g, const int npar,
					const lbfgsfloatval_t step);
  static int lbfgs_progress(void *instance, const lbfgsfloatval_t *x,
			    const lbfgsfloatval_t *g, const lbfgsfloatval_t fx,
			    const lbfgsfloatval_t xnorm, 
			    const lbfgsfloatval_t gnorm,
			    const lbfgsfloatval_t step,
			    int n, int k, int ls);
private:
  //internal variables
  std::string configFile_;         //!< The configuration file name
  std::string outputFile_;         //!< The output file name
  Parameters * par_;              //!< Fit parameters, depend on number of bins & geometry
  int fitMethod_;                  //!< Running mode
  std::string  minName_,algoName_; //!< ROOT::Minimizer name and algo
  int nThreads_;                   //!< Number of threads
  std::vector<Event*> data_;       //!< The data
  std::vector<Event*> control_[2]; //!< control samples
  int nGammajetEvents_;            //!< Number of gamma-jet events
  int nDijetEvents_;               //!< Number of dijet events
  int nTrijetEvents_;              //!< Number of trijet events
  int nTrackClusterEvents_;        //!< Number of track-cluster events
  int nZjetEvents_;                //!< Number of Zjet events
  int nTopEvents_;                 //!< Number of top events
  int mode_;

  // control parameters of fit
  bool printParNDeriv_;     //!< Control whether to print derivatives in each iteration
  std::vector<int> residualScalingScheme_;    //!< Iteration scheme of scaling of residuals
  double outlierChi2Cut_;                     //!< Cut on outlier when no scaling is chosen
  std::vector<int> fixedJetPars_;             //!< List of fixed jet parameters
  std::vector<int> fixedGlobalJetPars_;       //!< List of fixed global jet parameters

  // LVMINI parameters
  double derivStep_;        //!< Step width for derivative calculation
  int mvec_;                //!< Number of stored vector pairs in LVMINI
  int nIter_;               //!< Number of iterations in LVMINI
  float eps_;               //!< Convergence parameter in LVMINI
  float wlf1_;              //!< Parameter 1 of strong Wolfe condition in LVMINI
  float wlf2_;              //!< Parameter 2 of strong Wolfe condition in LVMINI
  bool calcCov_;            //!< If true, calculate covariance matrix of fitted parameters
  double *epsilon_;
  double *temp_derivative1_;
  double *temp_derivative2_;
  double *temp_derivative3_;
  double *temp_derivative4_;
  ComputeThread **threads_;
};

#endif
