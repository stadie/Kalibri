//
// Original Author:  Christian Autermann
//         Created:  Wed Jul 18 13:54:50 CEST 2007
// $Id: Parameters.h,v 1.28 2008/07/31 12:53:35 auterman Exp $
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
  
  int GetEtaBin(int const eta_id) const { return GetEtaBin(eta_id, eta_granularity, phi_granularity, eta_symmetry);}
  int GetPhiBin(int const phi_id) const { return GetPhiBin(phi_id, phi_granularity);}
  int GetJetEtaBin(int const eta_id) const { return GetEtaBin(eta_id, eta_granularity_jet, phi_granularity_jet, eta_symmetry);}
  int GetJetPhiBin(int const phi_id) const { return GetPhiBin(phi_id, phi_granularity_jet);}
  int GetBin(unsigned const etabin, unsigned const phibin) const {if (etabin<0) return etabin; else return etabin*phi_granularity + phibin;}
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
  
  static const double tower_parametrization(TMeasurement *const x,double *const par) {
    return instance->p->correctedTowerEt(x,par);
  }
  static const double jet_parametrization(TMeasurement *const x,double *const par) {
    return instance->p->correctedJetEt(x,par);
  }
  static const double dummy_parametrization(TMeasurement *const x,double *const par) {
    return x->pt;
  }

  //Error parametrization functions:
  template<int Et> static const double const_error(double *const x, TMeasurement *const xorig=0, double const errorig=0) {
    return Et;
  }
  static const double tower_error_parametrization(double *const x, TMeasurement *const xorig=0, double const errorig=0) {        
    return (x[0]>0 ?  1.25 * sqrt( x[0])   :   1.25 * sqrt(-x[0]) );  
  }
  static const double jet_error_parametrization(double *const x, TMeasurement *const xorig=0, double const errorig=0) {
    return (x[0]>0. ? 0.033*x[0] + 5.6   :   0.03*(-x[0]) + 5.6); 
  }
  static const double dummy_error_parametrization(double *const x, TMeasurement *const xorig=0, double const errorig=0) {        
    return x[0];  
  }
  static const double fast_error_parametrization(double *const x, TMeasurement *const xorig, double const errorig)  {
    return (xorig->pt=0. ? errorig : errorig*x[0]/xorig->pt );  
  }
  static const double jans_E_tower_error_parametrization(double *const x, TMeasurement *const xorig=0, double errorig=0)  {
    // E = x[0]*xorig[7];  x[0]=param. mess;    xorig == _mess
    double pmess;
    if(std::abs(xorig->eta) < 3.0)  
      pmess =  x[0] * xorig->E / (xorig->pt * xorig->pt) * (xorig->HadF + xorig->OutF); //Et->E hadronic
    else
      pmess =  x[0] * (xorig->E / xorig->pt);  //Et->E 
    return (xorig->E!=0. ? tower_error_parametrization(&pmess,xorig,errorig) * xorig->pt / xorig->E : 0.0);
  }


  //Plot paramterization stuff
  static const double plot_parametrization(TMeasurement *const x,double *const par) {
    return tower_parametrization(x,par)/x->pt; 
  }

  static const double jes_plot_parametrization(TMeasurement *const x,double *const par)  {
    return jet_parametrization(x,par)/x->pt;
  }

  //Limiting parameters
  static const double parameter_limit(TMeasurement *const x, double *const par) {
    double min = x->pt;
    double max = x->EMF;  //@@ Are you sure this is correct????
    if(par[0] < min) return (min-par[0]);
    if(par[0] > max) return (par[0]-max);
    return 0;
    //return 1e4/(1+exp(k* (par[0] - min))) + 1e4/(1+exp(-k* (par[0] - max));
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
  int GetEtaBin(int phi_id, int etagranu, int phigranu, bool etasym) const;
  int GetPhiBin(int phi_id, int phigranu) const;
  

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
