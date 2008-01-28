#ifndef TParameters_h
#define TParameters_h

//C++ libs
#include <vector>
#include <map>
#include <string>
#include <utility> //pair

#include <iostream>
#include <cmath>

class TParameters {  

public :
  TParameters(){};
  TParameters(std::string config):free_pars_per_bin(3),free_pars_per_bin_jet(2){
  this->ReadConfigFile(config);};
  virtual ~TParameters(){};
  friend std::ostream& operator<<( std::ostream& os, const TParameters& c );

  int GetEtaBin(int const eta_id) const;
  int GetPhiBin(int const phi_id) const; 
  int GetJetEtaBin(int const eta_id) const {return (abs(eta_id)-1)*2*eta_granularity_jet/eta_ntwr;};
  int GetJetPhiBin(int const phi_id) const {return (phi_id-1)*phi_granularity_jet/phi_ntwr;}; 
  //int GetJetEtaBin(double const eta_id) const {return fabs(eta_id)*eta_granularity_jet/4.0};
  //int GetJetPhiBin(double const phi_id) const {return phi_id*phi_granularity_jet/6.3;}; 
  int GetBin(unsigned const etabin, unsigned const phibin) const {if (etabin<0) return etabin; else return etabin*phi_granularity + phibin;};
  //int GetJetBin(unsigned const etabin, unsigned const phibin) const { if (etabin<0) return etabin; else return eta_granularity*phi_granularity*free_pars_per_bin+etabin*phi_granularity_jet + phibin;};
  int GetJetBin(unsigned const etabin, unsigned const phibin) const { if (etabin<0) return etabin; else return etabin*phi_granularity_jet + phibin;};
  virtual int GetNumberOfTowerParameters()const{return free_pars_per_bin*eta_granularity*phi_granularity;};
  virtual int GetNumberOfJetParameters()const{return   free_pars_per_bin_jet*eta_granularity_jet*phi_granularity_jet;};
  virtual int GetNumberOfParameters()const{return GetNumberOfTowerParameters()+GetNumberOfJetParameters();};

  double * GetTowerParRef(int const bin){return &k[bin*free_pars_per_bin]; };
  double * GetJetParRef(  int const jetbin){return &k[GetNumberOfTowerParameters()+jetbin*free_pars_per_bin_jet];};
  double * k; //all fit-parameters
  double * e; //all fit-parameter errors
  
  unsigned eta_granularity, phi_granularity,
           eta_granularity_jet, phi_granularity_jet;
  std::vector<double> start_values, jet_start_values;

  //The parametrization functions:
  unsigned free_pars_per_bin, free_pars_per_bin_jet;

//  static double tower_parametrization(double *x,double *par) {
//    return x[1] + par[0]*x[2] + par[1]*log(x[0]) + par[2] +x[3];};
//  static double jet_parametrization(double *x,double *par) {
//    return par[0]*x[0] + par[1];};
  static double dummy_parametrization(double *x,double *par) {
    return x[0];};

  static double tower_error_parametrization(double * x) {
    return (x[0]>0 ? 0.0*x[0] + 1.0*sqrt( x[0]) + 0. : 
	    0.0*x[0] - 1.0*sqrt(-x[0]) + 0.); };
  static double jet_error_parametrization(double * x) {
    return (x[0]>0. ? 0.05*x[0] + 1.0*sqrt( x[0]) + 0.:
	    0.05*x[0] - 1.0*sqrt(-x[0]) + 0.); };

//  static double plot_parametrization(double * x,double *par) {
//    return tower_parametrization(x,par)/x[0]; }

//  static double jes_plot_parametrization(double * x,double *par) {
//    return jet_parametrization(x,par)/x[0]; }


protected:
  void ReadConfigFile(std::string file);

private :
  //Towers in Eta-, Phi- direction (according to PTDR Vol I, p.201)
  unsigned const static eta_ntwr=82, phi_ntwr=72;
  unsigned eta_ntwr_used;
  bool eta_symmetry;

  //private functions
  void Read_Calibration(const std::string& file);
  std::string trim(std::string const& source, char const* delims = " {}\t\r\n");
  std::string input_calibration;
};

// Parametrization of hadronic response by a step function
class TStepParameters: public TParameters {
  public:
  TStepParameters(){};
  TStepParameters(std::string config){
  free_pars_per_bin=12;free_pars_per_bin_jet=2;ReadConfigFile(config);};
  
  static double tower_parametrization(double *x,double *par) {
    double result = 0;

    if      (x[2]>=0.0  && x[2]<=1.0)  result = x[1]+x[3] + par[0]*x[2];
    else if (x[2]>1.0   && x[2]<=2.0)  result = x[1]+x[3] + par[1]*x[2];
    else if (x[2]>2.0   && x[2]<=5.0)  result = x[1]+x[3] + par[2]*x[2];
    else if (x[2]>5.0   && x[2]<=10.0)  result = x[1]+x[3] + par[3]*x[2];
    else if (x[2]>10.0  && x[2]<=20.0)  result = x[1]+x[3] + par[4]*x[2];
    else if (x[2]>20.0  && x[2]<=40.0)  result = x[1]+x[3] + par[5]*x[2];
    else if (x[2]>40.0  && x[2]<=80.0) result = x[1]+x[3] + par[6]*x[2];
    else if (x[2]>80.0  && x[2]<=160.0) result = x[1]+x[3] + par[7]*x[2];
    else if (x[2]>160.0 && x[2]<=300.0) result = x[1]+x[3] + par[8]*x[2];
    else if (x[2]>300.0 && x[2]<=600.0) result = x[1]+x[3] + par[9]*x[2];
    else if (x[2]>600.0 && x[2]<=1000.0) result = x[1]+x[3] + par[10]*x[2];
    else if (x[2]>1000.0 )              result = x[1]+x[3] + par[11]*x[2];
    return result;};

  static double jet_parametrization(double *x,double *par) {
    return  par[0]*x[0] + par[1];};

  static double plot_parametrization(double * x,double *par) {
    return tower_parametrization(x,par)/x[0]; }

  static double jes_plot_parametrization(double * x,double *par) {
    return jet_parametrization(x,par)/x[0]; }

};

// Parametrization of hadronic response by a step function
// 3 Sets of Parameters for different EM fraction
class TStepEfracParameters: public TParameters {
  public:
  TStepEfracParameters(std::string config){
    free_pars_per_bin=36;free_pars_per_bin_jet=2;ReadConfigFile(config);};
  
  static double tower_parametrization(double *x,double *par) {
    double result=0;
    
    double Efrac = x[1]/(x[2]+x[3]);
    if (Efrac<0.1) {
      if      (x[2]>0.0   && x[2]<=1.0)    result = x[1]+x[3] + par[0]*x[2];
      else if (x[2]>1.0   && x[2]<=2.0)    result = x[1]+x[3] + par[1]*x[2];
      else if (x[2]>2.0   && x[2]<=5.0)    result = x[1]+x[3] + par[2]*x[2];
      else if (x[2]>5.0   && x[2]<=10.0)   result = x[1]+x[3] + par[3]*x[2];
      else if (x[2]>10.0  && x[2]<=20.0)   result = x[1]+x[3] + par[4]*x[2];
      else if (x[2]>20.0  && x[2]<=40.0)   result = x[1]+x[3] + par[5]*x[2];
      else if (x[2]>40.0  && x[2]<=80.0)   result = x[1]+x[3] + par[6]*x[2];
      else if (x[2]>80.0  && x[2]<=160.0)  result = x[1]+x[3] + par[7]*x[2];
      else if (x[2]>160.0 && x[2]<=300.0)  result = x[1]+x[3] + par[8]*x[2];
      else if (x[2]>300.0 && x[2]<=600.0)  result = x[1]+x[3] + par[9]*x[2];
      else if (x[2]>600.0 && x[2]<=1000.0) result = x[1]+x[3] + par[10]*x[2];
      else if (x[2]>1000.0 )               result = x[1]+x[3] + par[11]*x[2];
    } else if (Efrac<0.3) {
      if      (x[2]>0.0   && x[2]<=1.0)    result = x[1]+x[3] + par[12]*x[2];
      else if (x[2]>1.0   && x[2]<=2.0)    result = x[1]+x[3] + par[13]*x[2];
      else if (x[2]>2.0   && x[2]<=5.0)    result = x[1]+x[3] + par[14]*x[2];
      else if (x[2]>5.0   && x[2]<=10.0)   result = x[1]+x[3] + par[15]*x[2];
      else if (x[2]>10.0  && x[2]<=20.0)   result = x[1]+x[3] + par[16]*x[2];
      else if (x[2]>20.0  && x[2]<=40.0)   result = x[1]+x[3] + par[17]*x[2];
      else if (x[2]>40.0  && x[2]<=80.0)   result = x[1]+x[3] + par[18]*x[2];
      else if (x[2]>80.0  && x[2]<=160.0)  result = x[1]+x[3] + par[19]*x[2];
      else if (x[2]>160.0 && x[2]<=300.0)  result = x[1]+x[3] + par[20]*x[2];
      else if (x[2]>300.0 && x[2]<=600.0)  result = x[1]+x[3] + par[21]*x[2];
      else if (x[2]>600.0 && x[2]<=1000.0) result = x[1]+x[3] + par[22]*x[2];
      else if (x[2]>1000.0 )               result = x[1]+x[3] + par[23]*x[2];
    } else {
      if      (x[2]>0.0   && x[2]<=1.0)    result = x[1]+x[3] + par[24]*x[2];
      else if (x[2]>1.0   && x[2]<=2.0)    result = x[1]+x[3] + par[25]*x[2];
      else if (x[2]>2.0   && x[2]<=5.0)    result = x[1]+x[3] + par[26]*x[2];
      else if (x[2]>5.0   && x[2]<=10.0)   result = x[1]+x[3] + par[27]*x[2];
      else if (x[2]>10.0  && x[2]<=20.0)   result = x[1]+x[3] + par[28]*x[2];
      else if (x[2]>20.0  && x[2]<=40.0)   result = x[1]+x[3] + par[29]*x[2];
      else if (x[2]>40.0  && x[2]<=80.0)   result = x[1]+x[3] + par[30]*x[2];
      else if (x[2]>80.0  && x[2]<=160.0)  result = x[1]+x[3] + par[31]*x[2];
      else if (x[2]>160.0 && x[2]<=300.0)  result = x[1]+x[3] + par[32]*x[2];
      else if (x[2]>300.0 && x[2]<=600.0)  result = x[1]+x[3] + par[33]*x[2];
      else if (x[2]>600.0 && x[2]<=1000.0) result = x[1]+x[3] + par[34]*x[2];
      else if (x[2]>1000.0 )               result = x[1]+x[3] + par[35]*x[2];
    }
    return result;};
  
  static double jet_parametrization(double *x,double *par) {
    return  par[0]*x[0] + par[1];};

  static double plot_parametrization(double * x,double *par) {
    return tower_parametrization(x,par)/x[0]; }

  static double jes_plot_parametrization(double * x,double *par) {
    return jet_parametrization(x,par)/x[0]; }

};

// Parametrization of response by some "clever" function
class TMyParameters: public TParameters {
  public:
  TMyParameters(std::string config){
    free_pars_per_bin=3;free_pars_per_bin_jet=2;ReadConfigFile(config);};

  static double tower_parametrization(double *x,double *par) {
    return x[1] + par[0]*x[2] + par[1]*log(x[0]) + par[2];};

  static double jet_parametrization(double *x,double *par) {
    return par[0]*x[0] + par[1];};
  static double dummy_parametrization(double *x,double *par) {
    return x[0];};

  static double tower_error_parametrization(double * x) {
    return (x[0]>0 ? 0.0*x[0] + 1.0*sqrt( x[0]) + 0. : 
	    0.0*x[0] - 1.0*sqrt(-x[0]) + 0.); };
  static double jet_error_parametrization(double * x) {
    return (x[0]>0. ? 0.05*x[0] + 1.0*sqrt( x[0]) + 0.:
	    0.05*x[0] - 1.0*sqrt(-x[0]) + 0.); };

  static double plot_parametrization(double * x,double *par) {
    return tower_parametrization(x,par)/x[0]; }

  static double jes_plot_parametrization(double * x,double *par) {
    return jet_parametrization(x,par)/x[0]; }

};


#endif
