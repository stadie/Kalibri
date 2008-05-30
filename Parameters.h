//
// Original Author:  Christian Autermann
//         Created:  Wed Jul 18 13:54:50 CEST 2007
// $Id: Parameters.h,v 1.20 2008/05/08 17:13:17 auterman Exp $
//
#ifndef TParameters_h
#define TParameters_h

//C++ libs
#include <vector>
#include <map>
#include <string>
#include <utility> 

#include <iostream>
#include <cmath>

#include "Parametrization.h"

class ConfigFile;
class TParameters {  
public :
  
  static TParameters* CreateParameters(const std::string& configfile);

  std::string GetName() const;
  int GetEtaBin(int const eta_id) const;
  int GetPhiBin(int const phi_id) const; 
  int GetJetEtaBin(int const eta_id) const {return (abs(eta_id)-1)*2*eta_granularity_jet/eta_ntwr;}
  int GetJetPhiBin(int const phi_id) const {return (phi_id-1)*phi_granularity_jet/phi_ntwr;}
  //int GetJetEtaBin(double const eta_id) const {return fabs(eta_id)*eta_granularity_jet/4.0}
  //int GetJetPhiBin(double const phi_id) const {return phi_id*phi_granularity_jet/6.3;} 
  int GetBin(unsigned const etabin, unsigned const phibin) const {if (etabin<0) return etabin; else return etabin*phi_granularity + phibin;}
  //int GetJetBin(unsigned const etabin, unsigned const phibin) const { if (etabin<0) return etabin; else return eta_granularity*phi_granularity*free_pars_per_bin+etabin*phi_granularity_jet + phibin;}
  int GetJetBin(unsigned const etabin, unsigned const phibin) const { if (etabin<0) return etabin; else return etabin*phi_granularity_jet + phibin;}
  int GetNumberOfTowerParameters() const{return p->nTowerPars() *eta_granularity*phi_granularity;}
  int GetNumberOfJetParameters() const{return p->nJetPars()*eta_granularity_jet*phi_granularity_jet;}
  int GetNumberOfParameters() const{return GetNumberOfTowerParameters()+GetNumberOfJetParameters();}
  int GetNumberOfTowerParametersPerBin() const {return p->nTowerPars();}
  int GetNumberOfJetParametersPerBin() const {return p->nJetPars();}
  int GetEtaGranularity() const { return eta_granularity;}
  int GetPhiGranularity() const { return phi_granularity;}
  int GetEtaGranularityJet() const { return eta_granularity_jet;}
  int GetPhiGranularityJet() const { return phi_granularity_jet;}

  double* GetTowerParRef(int bin) { return k + bin*p->nTowerPars(); }
  double* GetJetParRef(int jetbin)  { return k + GetNumberOfTowerParameters()+jetbin*p->nJetPars();}
  void SetErrors(double *ne) { std::memcpy(e,ne,GetNumberOfParameters()*sizeof(double));}  
  void SetParameters(double *np) { std::memcpy(k,np,GetNumberOfParameters()*sizeof(double));}
  void SetFitChi2(double chi2) { fitchi2 = chi2;}
  double GetFitChi2() const { return fitchi2;}
  void FillErrors(double* copy) const {
    std::memcpy(copy,e,GetNumberOfParameters()*sizeof(double));
  }
  double* GetPars() { return k; }

  void Print() const;
  friend std::ostream& operator<<( std::ostream& os, const TParameters& c );
  
  static double tower_parametrization(double *x,double *par) {
    return instance->p->correctedTowerEt(x,par);
  }
  static double jet_parametrization(double *x,double *par) {
    return instance->p->correctedJetEt(x,par);
  }
  static double dummy_parametrization(double *x,double *par) {
    return x[0];
  }

  template<int Et> static double const_error(double * x) {
    return Et;
  }
  
  static double tower_error_parametrization(double * x) {
    return (x[0]>0 ? 0.0*x[0] + 1.0*sqrt( x[0]) + 0. : 
	    0.0*x[0] - 1.0*sqrt(-x[0]) + 0.); 
  }
  static double jet_error_parametrization(double * x) {
    return (x[0]>0. ? 0.05*x[0] + 1.0*sqrt( x[0]) + 0.:
	    0.05*x[0] - 1.0*sqrt(-x[0]) + 0.); 
  }
  
  static double plot_parametrization(double * x,double *par) {
    return tower_parametrization(x,par)/x[0]; 
  }

  static double jes_plot_parametrization(double * x,double *par) {
    return jet_parametrization(x,par)/x[0];
  }


protected:
  TParameters(Parametrization* p) 
    : p(p),k(0),e(0),fitchi2(0) {
  };
  virtual ~TParameters() {
    delete p;
    delete [] k;
    delete [] e;
  };
private:
  TParameters();
  TParameters(const TParameters&) {}

  //Towers in Eta-, Phi- direction (according to PTDR Vol I, p.201)
  unsigned const static eta_ntwr=82, phi_ntwr=72;
  unsigned eta_ntwr_used;
  bool eta_symmetry;
  unsigned int eta_granularity, phi_granularity,eta_granularity_jet, phi_granularity_jet;
  std::vector<double> start_values, jet_start_values;
  //The parametrization functions:
  Parametrization* p;

  double * k; //all fit-parameters
  double * e; //all fit-parameter errors
  double fitchi2;
  //private functions
  void Init(const ConfigFile& config);
  void Read_Calibration(const std::string& file);
  std::string trim(std::string const& source, char const* delims = " {}\t\r\n");
  std::string input_calibration;
  

  static TParameters *instance; 

  static Parametrization* CreateParametrization(const std::string& name);
  
  class Cleaner
  {
  public:
    Cleaner() {}
    ~Cleaner()
    {
      if(TParameters::instance) { 
	delete TParameters::instance; 
	TParameters::instance = 0; 
      }
    }
  };
  friend class Cleaner;
};

#endif
