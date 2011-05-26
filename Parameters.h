//
// Original Authors:  Christian Autermann, Hartmut Stadie
//         Created:  Wed Jul 18 13:54:50 CEST 2007
// $Id: Parameters.h,v 1.65 2010/11/01 15:47:40 stadie Exp $
//
#ifndef Parameters_h
#define Parameters_h

//C++ libs
#include <vector>
#include <map>
#include <string>
#include <utility> 

#include <iostream>
#include <cmath>
#include <cstring>

#include "ConfigFile.h"
#include "Parametrization.h"
#include "Function.h"
#include "ResolutionParametrization.h"
#include "ResolutionFunction.h"
#include "gsl/gsl_errno.h"
#include "gsl/gsl_math.h"
#include "gsl/gsl_roots.h"

class TH1;

//!  \brief Connection between detector geometry and fit parameters,
//!         interface to response and error parametrizations
//!  \author Christian Autermann, Hartmut Stadie
//!  \date   Wed Jul 18 13:54:50 CEST 2007
//!  $Id: Parameters.h,v 1.65 2010/11/01 15:47:40 stadie Exp $
// -----------------------------------------------------------------
class Parameters {  
 public:
  static Parameters* createParameters(const ConfigFile& config);

  int etaBin(int const eta_id) const { return etaBin(eta_id, eta_granularity_, phi_granularity_, eta_symmetry_);}
  int phiBin(int const phi_id) const { return phiBin(phi_id, phi_granularity_);}
  int jetEtaBin(int const eta_id) const { return etaBin(eta_id, eta_granularity_jet_, phi_granularity_jet_, eta_symmetry_);}
  int jetPhiBin(int const phi_id) const { return phiBin(phi_id, phi_granularity_jet_);}
  int trackEtaBin(int const eta_id) const { return etaBin(eta_id, eta_granularity_track_, phi_granularity_track_, eta_symmetry_);}
  int trackPhiBin(int const phi_id) const { return phiBin(phi_id, phi_granularity_jet_);}
  int bin(unsigned const etabin, unsigned const phibin) const {if (etabin<0) return etabin; else return etabin*phi_granularity_ + phibin;}
  int jetBin(unsigned const etabin, unsigned const phibin) const { if (etabin<0) return etabin; else return etabin*phi_granularity_jet_ + phibin;}
  int trackBin(unsigned const etabin, unsigned const phibin) const {if (etabin<0) return etabin; else return etabin*phi_granularity_track_ + phibin;}

  int numberOfTowerParameters() const{return p_->nTowerPars() *eta_granularity_*phi_granularity_;}
  int numberOfJetParameters() const{return p_->nJetPars()*eta_granularity_jet_*phi_granularity_jet_;}
  int numberOfTrackParameters() const{return p_->nTrackPars()*eta_granularity_track_*phi_granularity_track_;}
  int numberOfGlobalJetParameters() const{return p_->nGlobalJetPars();}
  int numberOfFixedParameters() const {
    int n = 0;
    for(std::vector<bool>::const_iterator it = isFixedPar_.begin();
	it != isFixedPar_.end(); it++) {
      if( *it ) n++;
    }
    return n;
  }

  int numberOfParameters() const{return numberOfTowerParameters()+numberOfJetParameters() + numberOfTrackParameters()+ numberOfGlobalJetParameters();}
  int numberOfTowerParametersPerBin() const {return p_->nTowerPars();}
  int numberOfJetParametersPerBin() const {return p_->nJetPars();}
  int numberOfTrackParametersPerBin() const {return p_->nTrackPars();}
  int numberOfCovCoeffs() const { 
    return (numberOfParameters()*numberOfParameters()+numberOfParameters())/2;
  }

  int etaGranularity() const { return eta_granularity_;}
  int phiGranularity() const { return phi_granularity_;}
  int etaGranularityJet() const { return eta_granularity_jet_;}
  int phiGranularityJet() const { return phi_granularity_jet_;}
  int etaGranularityTrack() const { return eta_granularity_track_;}
  int phiGranularityTrack() const { return phi_granularity_track_;}

  unsigned int nPtBins() const { return ptBinEdges_.size() - 1; }
  unsigned int nParPerPtBin() const { return resParam_->nParPerPtBin(); }
  bool findPtBin(double pt, unsigned int &bin) const { return findBin(pt,ptBinEdges_,bin); }
  double ptMin(unsigned int bin) const { return ptBinEdges_.at(bin); }
  double ptMax(unsigned int bin) const { return ptBinEdges_.at(bin+1); }
  double ptMin() const { return ptBinEdges_.front(); }
  double ptMax() const { return ptBinEdges_.back(); }
  double ptBinEdge(unsigned int bin) const { return ptBinEdges_.at(bin); }
  const std::vector<double>& ptBinEdges() const { return ptBinEdges_; }
  double ptTrueMin(unsigned int bin) const { return ptTrueMin_.at(bin); }
  double ptTrueMax(unsigned int bin) const { return ptTrueMax_.at(bin); }
  double ptTrueMin() const { return ptTrueMin_.front(); }
  double ptTrueMax() const { return ptTrueMax_.back(); }

  void writeCalibrationTxt(const char* name); //!< write calibration constants to txt file
  void writeCalibrationTex(const char* name, const ConfigFile& config); //!< write calibration constants and some paraemters of the fit to tex file

  double* towerParRef(int bin) { return k_ + bin*p_->nTowerPars(); }
  double* jetParRef(int jetbin)  { return k_ + numberOfTowerParameters()+jetbin*p_->nJetPars();}
  double* trackParRef(int trackbin)  { return k_ + numberOfTowerParameters() + numberOfJetParameters() +trackbin*p_->nTrackPars();}
  double* globalJetParRef()  { return k_ + numberOfTowerParameters() + numberOfJetParameters() + numberOfTrackParameters();}

  double* towerParErrorRef(int bin) { 
    return parErrors_ + bin*p_->nTowerPars();
  }
  double* jetParErrorRef(int jetbin)  { 
    return parErrors_ + numberOfTowerParameters()+jetbin*p_->nJetPars();
  }
  double* trackParErrorRef(int trackbin)  {
    return parErrors_ + numberOfTowerParameters() + numberOfJetParameters() +trackbin*p_->nTrackPars();
  }
  double* globalJetParErrorRef()  { 
    return parErrors_ + numberOfTowerParameters() + numberOfJetParameters() + numberOfTrackParameters();
  }

  double* towerParGlobalCorrCoeffRef(int bin) { 
    return parGCorr_ + bin*p_->nTowerPars();
  }
  double* jetParGlobalCorrCoeffRef(int jetbin)  { 
    return parGCorr_ + numberOfTowerParameters()+jetbin*p_->nJetPars();
  }
  double* trackParGlobalCorrCoeffRef(int trackbin)  {
    return parGCorr_ + numberOfTowerParameters() + numberOfJetParameters() +trackbin*p_->nTrackPars();
  }
  double* globalJetParGlobalCorrCoeffRef()  { 
    return parGCorr_ + numberOfTowerParameters() + numberOfJetParameters() + numberOfTrackParameters();
  }

  bool isFixedPar(int i) const { 
    assert( i >= 0 && i < numberOfParameters() );
    return isFixedPar_[i];
  }
  std::string parName(int i) const {
    assert( i >= 0 && i < numberOfParameters() );
    return parNames_[i];
  }

  void setParameters(double *np) {
    std::memcpy(k_,np,numberOfParameters()*sizeof(double));
  }
  void setErrors(double *ne) {
    std::memcpy(parErrors_,ne,numberOfParameters()*sizeof(double));
  }  
  void setGlobalCorrCoeff(double *gcc) {
    std::memcpy(parGCorr_,gcc,numberOfParameters()*sizeof(double));
  }  
  void setCovCoeff(double *cov) {
    std::memcpy(parCov_,cov,numberOfCovCoeffs()*sizeof(double));
  }
  void fixPar(int i) {
    assert( i >= 0 && i < numberOfParameters() );
    isFixedPar_[i] = true;
  }
  void fillErrors(double* copy) const {
    std::memcpy(copy,parErrors_,numberOfParameters()*sizeof(double));
  }
  double* parameters() { return k_; }
  double* errors() { return parErrors_; }
  double* globalCorrCoeff() { return parGCorr_; }
  double* covCoeff() { return parCov_; }
  double* effMap() {return trackEff_;}
  int trackEffBin(double pt, double eta);

  double jetStartPar(unsigned int i) const { return jet_start_values_.at(i); }

  void print() const;


  void printFuncs() {
    std::cout << "funcs: " << &funcmap_ << " \n";
    for(FunctionMap::const_iterator i = funcmap_.begin() ; 
	i != funcmap_.end() ; ++i) {
      std::cout << " f:" << i->first.intVal() << " " << &(i->second) << ", " << i->second << '\n';
    }
  }


  bool needsUpdate() const { return resParam_->needsUpdate(); }
  void update() { resParam_->update(parameters()); }
  
  //Error parametrization functions:
  template<int Et> static float const_error(const float *x, const Measurement *xorig=0, float errorig=0) {
    return Et;
  }
  static float tower_error_parametrization(const float *x, const Measurement *xorig=0, float errorig=0) { 
    return (x[0]>0 ?  1.25 * sqrt( x[0])   :   1.25 * sqrt(-x[0]) );
  }
  static float jet_error_parametrization(const float *x, const Measurement *xorig=0, float errorig=0) {
    return (x[0]>0. ? 0.033*x[0] + 5.6   :   0.033*(-x[0]) + 5.6 ); 
  }

  static float track_error_parametrization(const float *x, const Measurement *xorig=0, float errorig=0) { 
    //for full error also see Grooms paper 0605164v4, p.25
    float error=0,error2=0;
    error =  (x[0]>0 ? x[0] *( 0.05 + 0.00015 * x[0])   : (-x[0]) *(  0.05 + 0.00015 * (-x[0]) )); //trackerror  to be checken and dependent on pt, eta, chi2, nohits, ....
    error2 = error * error;

    //Pi0 Fehler s.Clemens
    if(x[0] > 3)      error = x[0] * 0.15 + 3;              //p. 70
    else              error = x[0] * 1.15;
    error2 += error * error;

    //error2 += (1-1/1.48)*(1-1/1.48)*0.125*0.125*x[0]*x[0];   //*(x[0]/100)^(-0.076)          //1/1.48 = h/e
    //following term has to be checked!!!!
    //float a = 1/(1.48 * 1.48) * 1.25 * 1.25 * pow((fabs(x[0])* (xorig->E / xorig->pt) / 0.96),(0.816 - 1));  // 1- Pi0 * error(h)^2 (h/e)^2
    //error2 += (x[0]>0 ?  a * x[0]  : a * (-x[0]));    //intrinsic term (HCAL)
    error = sqrt(error2);
    return error;
  }


  static float jet_only_tower_error_parametrization(const float *x, const Measurement *xorig=0, float errorig=0) { 
    return 0;
  }



  //!  \brief Parameters from V. Chetluru's fit to L2L3 corrected jets
  //!
  //!  Use results from V. Chetluru's talk:
  //!  <A HREF="http://indico.cern.ch/getFile.py/access?contribId=1&resId=1&materialId=slides&confId=52598">
  //!  Jet energy resolution studies
  //!  </A>.
  //!
  //!  The absolute resolution is given by
  //!  \f[
  //!   \sigma^{2} = a^{2} + b^{2}p_{T} + c^{2}p^{2}_{T}
  //!  \f]
  //!  with the \f$ \eta \f$ dependent parameters
  //!  <TABLE>
  //!   <TR>
  //!    <TD>  </TD>
  //!    <TD> a </TD>
  //!    <TD> b </TD>
  //!    <TD> c </TD>
  //!   </TR>
  //!   <TR>
  //!    <TD> \f$ 0 < \eta < 0.8 \f$ </TD>
  //!    <TD> 4.44 </TD>
  //!    <TD> 1.11 </TD>
  //!    <TD> 0.03 </TD>
  //!   </TR>
  //!   <TR>
  //!    <TD> \f$ 0.8 < \eta < 1.5 \f$ </TD>
  //!    <TD> 4.35 </TD>
  //!    <TD> 1.17 </TD>
  //!    <TD> 0.04 </TD>
  //!   </TR>
  //!   <TR>
  //!    <TD> \f$ 1.5 < \eta < 2.4 \f$ </TD>
  //!    <TD> 4.34 </TD>
  //!    <TD> 0.85 </TD>
  //!    <TD> 0.03 </TD>
  //!   </TR>
  //!   <TR>
  //!    <TD> \f$ 2.4 < \eta < 3.2 \f$ </TD>
  //!    <TD> 4.08 </TD>
  //!    <TD> 0.45 </TD>
  //!    <TD> 0.04 </TD>
  //!   </TR>
  //!   <TR>
  //!    <TD> \f$ 3.2 < \eta \f$ </TD>
  //!    <TD> 3.90 </TD>
  //!    <TD> 0.29 </TD>
  //!    <TD> 0.09 </TD>
  //!   </TR>
  //!  </TABLE>
  //!
  //!  \return The absolute resolution
  // -----------------------------------------------------
  static float jet_only_jet_error_parametrization_et(const float *x, const Measurement *xorig=0, float errorig=0) {
    const static float a[5] = { 4.44 * 4.44, 4.35 * 4.35, 4.34 * 4.34 , 4.08 * 4.08, 3.90 * 3.90 };
    const static float b[5] = { 1.11 * 1.11, 1.17 * 1.17, 0.85 * 0.85, 0.45 * 0.45, 0.29 * 0.29};
    const static float c[5] = { 0.03 * 0.03, 0.04 * 0.04, 0.03 * 0.03, 0.04 * 0.04, 0.09 * 0.09};

    float abseta = std::abs(xorig->eta);
    int i = (abseta < 0.8) ? 0 : ((abseta < 1.5) ? 1 : ((abseta < 2.4) ? 2 : (abseta < 3.2) ? 3 : 4));
    return sqrt(a[i] + (b[i] + c[i] *x[0]) * x[0]);
  }

  static float jet_only_jet_error_parametrization_energy(const float *x, const Measurement *xorig=0, float errorig=0) {
    /*
    float pmess;
    if(std::abs(xorig->eta) < 3.0)  
      pmess =  x[0] * xorig->E / (xorig->pt * xorig->pt) * (xorig->HadF + xorig->OutF); //Et->E hadronic
    else
      pmess =  x[0] * (xorig->E / xorig->pt);  //Et->E 
    //constant before stochastic term is not properly knowen
    return (x[0]>0. ? 0.033*x[0] + 5.6 + 1.0 * sqrt(pmess)  :   0.033*(-x[0]) + 5.6 + 1.0 * sqrt(-pmess) ); 
    */
    float E = x[0] * xorig->E/xorig->pt;
    //float sqE = sqrt(E);
    return sqrt(1.3*1.3/E + 0.056 * 0.056) * x[0];
  }


  static float dummy_error_parametrization(const float *x, const Measurement *xorig=0, float errorig=0) {        
    return x[0];  
  }
  static float fast_error_parametrization(const float *x, const Measurement *xorig, float errorig)  {
    return (xorig->pt==0. ? errorig : errorig*x[0]/xorig->pt );  
  }
  static float jans_E_tower_error_parametrization(const float *x, const Measurement *xorig=0, float errorig=0)  {
    
    // E = x[0]*xorig[7];  x[0]=param. mess;    xorig == _mess
    float pmess;
    if(std::abs(xorig->eta) < 3.0)  
      pmess =  x[0] * xorig->E / (xorig->pt * xorig->pt) * (xorig->HadF + xorig->OutF); //Et->E hadronic
    else
      pmess =  x[0] * (xorig->E / xorig->pt);  //Et->E 
    return (xorig->E!=0. ? tower_error_parametrization(&pmess,xorig,errorig) * xorig->pt / xorig->E : 0.0);
    
    //return 0;
  }
  
  static float toy_tower_error_parametrization(const float *x, const Measurement *xorig=0, float errorig=0);
  
  static float toy_jet_error_parametrization(const float *x, const Measurement *xorig=0, float errorig=0);
  
  static float const_error_parametrization(const float *x, const Measurement *xorig, float errorig)  {
    return errorig;  
  }
  
  float etaEdge(int const etaBin, bool lowerEdge);
  //! return upper edge of bin in eta
  float etaUpperEdge(int const etaBin) { return etaEdge(etaBin, false); };
  //! return lower edge of bin in eta
  float etaLowerEdge(int const etaBin) { return etaEdge(etaBin, true ); };

  // Return parametrization functions
  const Function& tower_function(int etaid, int phiid);
  const Function& jet_function(int etaid, int phiid);
  const Function& track_function(int etaid, int phiid);
  const Function& global_jet_function();
  const Function& function(const Function& f);
  const ResolutionFunction& function(const ResolutionFunction& f);
  const ResolutionFunction& resolutionFitPDF(unsigned int ptBin, int etaid, int phiid);

  void readCalibrationCfi(const std::string& file);
  void readCalibrationTxt(const std::string& file);
  void readCalibrationJetMET(const std::vector<std::string>& inputFileNames);
  void readCalibrationJetMETL2(const std::string& inputFileName);
  void readCalibrationJetMETL3(const std::string& inputFileName);
  
  //static const Parametrization* parametrization() { return instance->p;}
  
  Parameters* clone() const;
  static void removeClone(Parameters* p);
  
  bool findRoot(double (* f) (double x, void * params), void* params,
		double& x1, double& x2, double eps); 
  
  struct Variation {
    int    parid;        //!< Id of varied parameter
    float upperEt;      //!< Expected Et if parameter is varied by +eps
    float lowerEt;      //!< Expected Et if parameter is varied by -eps
    float upperError;   //!< Expected error if parameter is varied by +eps
    float lowerError;   //!< Expected error if parameter is varied by -eps
    float upperEtDeriv; //!< Derivative of Et if parameter is  varied by +eps
    float lowerEtDeriv; //!< Derivative of Et if parameter is  varied by +eps
    bool operator==(int b) const { return parid == b;} //!< Two ParameterVariation are the same if they have the same parid
  };
  typedef std::vector<Variation> VariationColl;
  typedef std::vector<Variation>::const_iterator VariationCollIter;


  VariationColl& cachedVariationColl() const { return cachedvariationcoll_;}


 protected:
  Parameters();
  Parameters(const Parametrization& p) {};
  Parameters(Parametrization* p);
  Parameters(ResolutionParametrization *resParam); 
  virtual ~Parameters();
  Parameters& operator=(const Parameters& p);


 private:
  int etaBin(int phi_id, int etagranu, int phigranu, bool etasym) const;
  int phiBin(int phi_id, int phigranu) const;
  //! Return one line of LaTeX tabular containing the name and value of a given parameter from config file
  template<class T> std::string texTabularLine(const ConfigFile& config, const std::string& fieldname) const;
  //! Return submatrix of covariance matrix for \p nPar parameters from \p firstPar
  std::vector<int> findCovIndices(int firstPar, int nPar) const;
  //! Return stati (is fixed?) for \p nPar parameters from \p firstPar
  std::vector<bool> findParStatus(int firstPar, int nPar) const;

  //Towers in Eta-, Phi- direction (according to PTDR Vol I, p.201)
  static  const unsigned int eta_ntwr_=82, phi_ntwr_=72;
  unsigned int eta_ntwr_used_;
  bool eta_symmetry_;
  unsigned int eta_granularity_, phi_granularity_,eta_granularity_jet_, phi_granularity_jet_, eta_granularity_track_, phi_granularity_track_;
  std::vector<double> start_values_, jet_start_values_, track_start_values_, global_jet_start_values_;
  std::vector<std::string> parNames_;

  //The parametrization functions:
  static Parametrization* p_;
  static ResolutionParametrization* resParam_;

  double * k_; //!< all fit-parameters
  std::vector<bool> isFixedPar_;
  double * parErrors_; //!< all fit-parameter errors
  double * parGCorr_; //!< Global correlation coefficients of parameters
  double * parCov_;
  double * trackEff_; //!< track Efficiency 13eta X 13 ptbins;

  std::vector<double> ptBinEdges_;
  std::vector<double> ptBinCenters_;
  std::vector<double> ptTrueMin_;
  std::vector<double> ptTrueMax_;

  /// ------------------------------------------------------
  /// private functions

  void init(const ConfigFile& config);
  void readTrackEffTxt(const std::string& file);
  std::string trim(std::string const& source, char const* delims = " {}\t\r\n");
  bool findBin(double x, const std::vector<double> &binEdges, unsigned int &bin) const;

  static Parameters *instance_;
  static std::vector<Parameters*> clones_;

  static Parametrization* createParametrization(const std::string& name, const ConfigFile& config);
  static ResolutionParametrization* createResolutionParametrization(const std::string& name, const ConfigFile& config);

  
  class Cleaner
  {
  public:
    Cleaner() {}
    ~Cleaner()
    {
      //delete clones
      for(std::vector<Parameters*>::const_iterator i = clones_.begin();
	  i != clones_.end() ; ++i) {
	delete *i;
      }
      if(Parameters::instance_) { 
	delete Parameters::instance_; 
	Parameters::instance_ = 0; 
      }  
      delete p_;
    }
  };
  friend class Cleaner;
  
  enum FunctionType { Tower, Jet, Global, Track, Resolution };
  
  class FunctionID {
    FunctionType t_;
    const unsigned short int i_;
  public:
    FunctionID(const Function::ParametrizationFunction& f, unsigned short int i)
      : i_(i) 
      {
	if(f == &Parametrization::expectedResponse) t_ = Track;
	else if(f == &Parametrization::correctedGlobalJetEt) t_ = Global;
	else if(f == &Parametrization::correctedJetEt) t_ = Jet;
	else if(f == &Parametrization::correctedTowerEt) t_ = Tower;
	else {
	  exit(12);
	}
      }
      FunctionID(FunctionType t, unsigned short int i)
	: t_(t), i_(i) 
	{}  
	
	bool operator<(const FunctionID& r) const {
	  if(i_ < r.i_) return true;
	  if(i_ > r.i_) return false;
	  if(t_ < r.t_) return true;
	  return false;
	}
	int intVal() const {
	  return t_*10000 + i_;
	}
  };
  typedef std::map<FunctionID,Function*> FunctionMap;
  FunctionMap funcmap_;
  mutable VariationColl cachedvariationcoll_;
  //root finding
  gsl_root_fsolver* s_;
  gsl_function F_; 
  static long long ncalls_;        //!< Number of calls of inversion methods 
  static long long ntries_;        //!< Number of tries in iteration during inversion
  static long long nfails_;        //!< Number of failed tries during inversion
  static long long nwarns_;        //!< Number of warnings during inversion
};




#endif
