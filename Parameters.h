//
// Original Author:  Christian Autermann
//         Created:  Wed Jul 18 13:54:50 CEST 2007
// $Id: Parameters.h,v 1.38 2008/11/14 12:41:37 thomsen Exp $
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
#include <cstring>

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
  int GetTrackEtaBin(int const eta_id) const { return GetEtaBin(eta_id, eta_granularity_track, phi_granularity_track, eta_symmetry);}
  int GetTrackPhiBin(int const phi_id) const { return GetPhiBin(phi_id, phi_granularity_jet);}
  int GetBin(unsigned const etabin, unsigned const phibin) const {if (etabin<0) return etabin; else return etabin*phi_granularity + phibin;}
  int GetJetBin(unsigned const etabin, unsigned const phibin) const { if (etabin<0) return etabin; else return etabin*phi_granularity_jet + phibin;}
  int GetTrackBin(unsigned const etabin, unsigned const phibin) const {if (etabin<0) return etabin; else return etabin*phi_granularity_track + phibin;}

  int GetNumberOfTowerParameters() const{return p->nTowerPars() *eta_granularity*phi_granularity;}
  int GetNumberOfJetParameters() const{return p->nJetPars()*eta_granularity_jet*phi_granularity_jet;}
  int GetNumberOfTrackParameters() const{return p->nTrackPars()*eta_granularity_track*phi_granularity_track;}
  int GetNumberOfParameters() const{return GetNumberOfTowerParameters()+GetNumberOfJetParameters() + GetNumberOfTrackParameters();}
  int GetNumberOfTowerParametersPerBin() const {return p->nTowerPars();}
  int GetNumberOfJetParametersPerBin() const {return p->nJetPars();}
  int GetNumberOfTrackParametersPerBin() const {return p->nTrackPars();}

  int GetEtaGranularity() const { return eta_granularity;}
  int GetPhiGranularity() const { return phi_granularity;}
  int GetEtaGranularityJet() const { return eta_granularity_jet;}
  int GetPhiGranularityJet() const { return phi_granularity_jet;}
  int GetEtaGranularityTrack() const { return eta_granularity_track;}
  int GetPhiGranularityTrack() const { return phi_granularity_track;}

  /// write calibration constants to cfi file
  void Write_CalibrationCfi(const char* name); 
  /// write calibration constants to txt file
  void Write_CalibrationTxt(const char* name); 

  double* GetTowerParRef(int bin) { return k + bin*p->nTowerPars(); }
  double* GetJetParRef(int jetbin)  { return k + GetNumberOfTowerParameters()+jetbin*p->nJetPars();}
  double* GetTrackParRef(int trackbin)  { return k + GetNumberOfTowerParameters() + GetNumberOfJetParameters() +trackbin*p->nTrackPars();}
  void SetErrors(double *ne) { std::memcpy(e,ne,GetNumberOfParameters()*sizeof(double));}  
  void SetParameters(double *np) { std::memcpy(k,np,GetNumberOfParameters()*sizeof(double));}
  void SetFitChi2(double chi2) { fitchi2 = chi2;}
  double GetFitChi2() const { return fitchi2;}
  void FillErrors(double* copy) const {
    std::memcpy(copy,e,GetNumberOfParameters()*sizeof(double));
  }
  double* GetPars() { return k; }
  double* GetErrors() { return e; }

  void Print() const;
  
  static const double tower_parametrization(TMeasurement *const x,double *const par) {
    return instance->p->correctedTowerEt(x,par);
  }
  static const double jet_parametrization(TMeasurement *const x,double *const par) {
    return instance->p->correctedJetEt(x,par);
  }  
  static const double track_parametrization(TMeasurement *const x,double *const par) {
    return instance->p->GetExpectedResponse(x,par);
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
    return (x[0]>0. ? 0.033*x[0] + 5.6   :   0.033*(-x[0]) + 5.6 ); 
  }

  static const double track_error_parametrization(double *const x, TMeasurement *const xorig=0, double const errorig=0) { 
    double error=0,error2=0;
    error =  (x[0]>0 ? x[0] *( 0.05 + 0.00015 * x[0])   : (-x[0]) *(  0.05 + 0.00015 * (-x[0]) )); //trackerror
    //for full error see Grooms paper 0605164v4, p.25
    error2 = error * error;
    //Pi0 Fehler
    error2 += (1-1/1.48)*(1-1/1.48)*0.125*0.125*x[0]*x[0];   //*(x[0]/100)^(-0.076)          //1/1.48 = h/e
    //folling term has to be checked!!!!
    double a = 1/(1.48 * 1.48) * 1.25 * 1.25 * pow((fabs(x[0])* (xorig->E / xorig->pt) / 0.96),(0.816 - 1));  // 1- Pi0 * error(h)^2 (h/e)^2
    error2 += (x[0]>0 ?  a * x[0]  : a * (-x[0]));    //intrinsic term (HCAL)
    error = sqrt(error2);
    return error;
  }


  static const double jet_only_tower_error_parametrization(double *const x, TMeasurement *const xorig=0, double const errorig=0) { 
    return 0;
  }
  static const double jet_only_jet_error_parametrization_et(double *const x, TMeasurement *const xorig=0, double const errorig=0) {
    return (x[0]>0. ? 0.033*x[0] + 5.6 + 1.25 * sqrt( x[0])   :   0.033*(-x[0]) + 5.6 + 1.25 * sqrt( -x[0]) ); 
  }
  static const double jet_only_jet_error_parametrization_energy(double *const x, TMeasurement *const xorig=0, double const errorig=0) {
    double pmess;
    if(std::abs(xorig->eta) < 3.0)  
      pmess =  x[0] * xorig->E / (xorig->pt * xorig->pt) * (xorig->HadF + xorig->OutF); //Et->E hadronic
    else
      pmess =  x[0] * (xorig->E / xorig->pt);  //Et->E 
    //constant before stochastic term is not properly knowen
    return (x[0]>0. ? 0.033*x[0] + 5.6 + 1.0 * sqrt(pmess)  :   0.033*(-x[0]) + 5.6 + 1.0 * sqrt(-pmess) ); 
  }



  static const double dummy_error_parametrization(double *const x, TMeasurement *const xorig=0, double const errorig=0) {        
    return x[0];  
  }
  static const double fast_error_parametrization(double *const x, TMeasurement *const xorig, double const errorig)  {
    return (xorig->pt==0. ? errorig : errorig*x[0]/xorig->pt );  
  }
  static const double jans_E_tower_error_parametrization(double *const x, TMeasurement *const xorig=0, double errorig=0)  {
    
    // E = x[0]*xorig[7];  x[0]=param. mess;    xorig == _mess
    double pmess;
    if(std::abs(xorig->eta) < 3.0)  
      pmess =  x[0] * xorig->E / (xorig->pt * xorig->pt) * (xorig->HadF + xorig->OutF); //Et->E hadronic
    else
      pmess =  x[0] * (xorig->E / xorig->pt);  //Et->E 
    return (xorig->E!=0. ? tower_error_parametrization(&pmess,xorig,errorig) * xorig->pt / xorig->E : 0.0);
    
    //return 0;
  }

  static const double toy_tower_error_parametrization(double *const x, TMeasurement *const xorig=0, double const errorig=0) {        
    double hade = x[0] / (xorig->pt *xorig->pt )* xorig->HadF * xorig->E;
    return sqrt( 1.3 * 1.3 / hade + 0.056 * 0.056) * hade * xorig->pt / xorig->E;
  }
  
  static const double toy_jet_error_parametrization(double *const x, TMeasurement *const xorig=0, double const errorig=0) {
    return 0;
  }

  static const double const_error_parametrization(double *const x, 
						  TMeasurement *const xorig, 
						  double const errorig)  {
    return errorig;  
  }

  //Plot paramterization stuff
  static const double plot_parametrization(TMeasurement *const x,double *const par) {
    return tower_parametrization(x,par)/x->pt; 
  }

  static const double jes_plot_parametrization(double * x,double * par)  {
    //return jet_parametrization(x,par)/x->pt;
    return ( 1. + 0.295 * par[0] * exp(- 0.02566 * par[1] * x[0]));   

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

  /// return upper or lower eta eta edge
  float EtaEdge(int const etaBin, bool lowerEdge);
  /// return upper edge of bin in eta
  float EtaUpperEdge(int const etaBin) { return EtaEdge(etaBin, false); };
  /// return lower edge of bin in eta
  float EtaLowerEdge(int const etaBin) { return EtaEdge(etaBin, true ); };


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
  unsigned int eta_granularity, phi_granularity,eta_granularity_jet, phi_granularity_jet, eta_granularity_track, phi_granularity_track;
  std::vector<double> start_values, jet_start_values, track_start_values;
  //The parametrization functions:
  Parametrization* p;

  double * k; //all fit-parameters
  double * e; //all fit-parameter errors
  double fitchi2;

  /// ------------------------------------------------------
  /// private functions

  void Init(const ConfigFile& config);

  /// read predefined calibration constants from cfi file 
  void Read_CalibrationCfi(const std::string& file);
  /// read predefined calibration constants from txt file
  void Read_CalibrationTxt(const std::string& file);

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
