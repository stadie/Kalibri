//
// Original Author:  Christian Autermann
//         Created:  Wed Jul 18 13:54:50 CEST 2007
// $Id: Parameters.h,v 1.44 2009/02/18 17:51:38 stadie Exp $
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
#include "Function.h"

#include "TMath.h"

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
  int GetNumberOfGlobalJetParameters() const{return p->nGlobalJetPars();}

  int GetNumberOfParameters() const{return GetNumberOfTowerParameters()+GetNumberOfJetParameters() + GetNumberOfTrackParameters()+GetNumberOfGlobalJetParameters();}
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
  double* GetGlobalJetParRef()  { return k + GetNumberOfTowerParameters() + GetNumberOfJetParameters() + GetNumberOfTrackParameters();}

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
  
  static double tower_parametrization(const TMeasurement* x, const double* par) {
    return instance->p->correctedTowerEt(x,par);
  }
  static double jet_parametrization(const TMeasurement* x, const double* par) {
    return instance->p->correctedJetEt(x,par);
  }  
  static double track_parametrization(const TMeasurement* x, const double* par) {
    return instance->p->GetExpectedResponse(x,par);
  }
  static double global_jet_parametrization(const TMeasurement* x, const double* par) {
    return instance->p->correctedGlobalJetEt(x,par);
  }

  static double dummy_parametrization(const TMeasurement* x, const double* par) {
    return x->pt;
  }

  //Error parametrization functions:
  template<int Et> static double const_error(const double *x, const TMeasurement *xorig=0, double errorig=0) {
    return Et;
  }
  static double tower_error_parametrization(const double *x, const TMeasurement *xorig=0, double errorig=0) { 
    return (x[0]>0 ?  1.25 * sqrt( x[0])   :   1.25 * sqrt(-x[0]) );
  }
  static double jet_error_parametrization(const double *x, const TMeasurement *xorig=0, double errorig=0) {
    return (x[0]>0. ? 0.033*x[0] + 5.6   :   0.033*(-x[0]) + 5.6 ); 
  }

  static double track_error_parametrization(const double *x, const TMeasurement *xorig=0, double errorig=0) { 
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


  static double jet_only_tower_error_parametrization(const double *x, const TMeasurement *xorig=0, double errorig=0) { 
    return 0;
  }

  static double jet_only_jet_error_parametrization_et(const double *x, const TMeasurement *xorig=0, double errorig=0) {
    //use results from V. Chetluru
    //http://indico.cern.ch/getFile.py/access?contribId=1&resId=1&materialId=slides&confId=52598
    // rel. sigma^2 = a^2/pt^2 + b^2/pt + c^2
    const static double a[5] = { 4.44 , 4.35 , 4.34 , 4.08 , 3.90 };
    const static double b[5] = { 1.11 , 1.17 , 0.85 , 0.45 , 0.29 };
    const static double c[5] = { 0.03 , 0.04 , 0.03 , 0.04 , 0.09 };

    double abseta = std::abs(xorig->eta);
    int i = (abseta < 0.8) ? 0 : ((abseta < 1.5) ? 1 : ((abseta < 2.4) ? 2 : (abseta < 3.2) ? 3 : 4));
    return sqrt(a[i]*a[i]/x[0]/x[0] + b[i]*b[i]/x[0] + c[i]*c[i]) * x[0];
  }

  static double jet_only_jet_error_parametrization_energy(const double *x, const TMeasurement *xorig=0, double errorig=0) {
    /*
    double pmess;
    if(std::abs(xorig->eta) < 3.0)  
      pmess =  x[0] * xorig->E / (xorig->pt * xorig->pt) * (xorig->HadF + xorig->OutF); //Et->E hadronic
    else
      pmess =  x[0] * (xorig->E / xorig->pt);  //Et->E 
    //constant before stochastic term is not properly knowen
    return (x[0]>0. ? 0.033*x[0] + 5.6 + 1.0 * sqrt(pmess)  :   0.033*(-x[0]) + 5.6 + 1.0 * sqrt(-pmess) ); 
    */
    double E = x[0] * xorig->E/xorig->pt;
    //double sqE = sqrt(E);
    return sqrt(1.3*1.3/E + 0.056 * 0.056) * x[0];
  }


  static double dummy_error_parametrization(const double *x, const TMeasurement *xorig=0, double errorig=0) {        
    return x[0];  
  }
  static double fast_error_parametrization(const double *x, const TMeasurement *xorig, double errorig)  {
    return (xorig->pt==0. ? errorig : errorig*x[0]/xorig->pt );  
  }
  static double jans_E_tower_error_parametrization(const double *x, const TMeasurement *xorig=0, double errorig=0)  {
    
    // E = x[0]*xorig[7];  x[0]=param. mess;    xorig == _mess
    double pmess;
    if(std::abs(xorig->eta) < 3.0)  
      pmess =  x[0] * xorig->E / (xorig->pt * xorig->pt) * (xorig->HadF + xorig->OutF); //Et->E hadronic
    else
      pmess =  x[0] * (xorig->E / xorig->pt);  //Et->E 
    return (xorig->E!=0. ? tower_error_parametrization(&pmess,xorig,errorig) * xorig->pt / xorig->E : 0.0);
    
    //return 0;
  }

  static double toy_tower_error_parametrization(const double *x, const TMeasurement *xorig=0, double errorig=0) {        
    double hadet = x[0] - xorig->EMF - xorig->OutF;
    if(hadet < 0.001) hadet = 0.001;
    double hade = hadet * xorig->E / xorig->pt; 
    //std::cout << "had Et:" << hadet << " , " << "had E:" << hade << '\n';
    double var = 1.3 * 1.3/hade + 0.056 * 0.056;  
    //truncate variance accordingly
    double truncvar = - sqrt(var) * exp(-0.5/var) * sqrt(2/M_PI) + var * TMath::Erf(1/(sqrt(2 * var)));
    return sqrt(truncvar) * hadet;
  }
  
  static double toy_jet_error_parametrization(const double *x, const TMeasurement *xorig=0, double errorig=0) {
    return 0;
  }

  static double const_error_parametrization(const double *x, const TMeasurement *xorig, double errorig)  {
    return errorig;  
  }

  //Plot paramterization stuff
  static double plot_parametrization(const TMeasurement* x, const double* par) {
    return tower_parametrization(x,par)/x->pt; 
  }

  static double jes_plot_parametrization(double * x,double * par)  {
    //return jet_parametrization(x,par)/x->pt;
    return ( 1. + 0.295 * par[0] * exp(- 0.02566 * par[1] * x[0]));   

  }

  //Limiting parameters
  static double parameter_limit(const TMeasurement* x, const double *par) {
    double min = x->pt;
    double max = x->EMF;
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

  Function tower_function(int etaid, int phiid);
  Function jet_function(int etaid, int phiid);
  Function track_function(int etaid, int phiid);
  Function global_jet_function();


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
  std::vector<double> start_values, jet_start_values, track_start_values, global_jet_start_values;
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
