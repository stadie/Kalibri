
//!  \author Hartmut Stadie
//!  \date Thu Apr 03 17:09:50 CEST 2008
//!  $Id: Parametrization.h,v 1.34 2009/04/15 17:26:59 mschrode Exp $

#ifndef CALIBCORE_PARAMETRIZATION_H
#define CALIBCORE_PARAMETRIZATION_H

#include <cmath>
#include "CalibData.h"

//!  \brief Abstract base class for parametrizations of
//!         correction functions
//!
//!  Interface to different parametrizations for the 
//!  underlying hypothesis of the GlobalFit. Allows
//!  to correct a tower or jet measurement.
// -----------------------------------------------------------------
class Parametrization 
{
public:
  //!  \brief Constructor
  //!  \param ntowerpars Number of parameters for tower parametrization
  //!  \param njetpars Number of parameters for (eta,phi)-dependent jet parametrization
  //!  \param ntrackpars Number of parameters for track parametrization
  //!  \param nglobaljetpars Number of parameters for global jet parametrization
  Parametrization(unsigned int ntowerpars, unsigned int njetpars, 
		  unsigned int ntrackpars, unsigned int nglobaljetpars) : 
    ntowerpars_(ntowerpars), njetpars_(njetpars), ntrackpars_(ntrackpars), 
    nglobaljetpars_(nglobaljetpars) {}

  virtual ~Parametrization() {}


  //!  \brief Corrects the measured tower Et
  //!
  //!  The parameters of the tower correction function
  //!  may depend on the tower Et, eta, and phi.
  //!
  //!  \param x Tower Et measurement that is to be corrected.
  //!           The measurement contains the following
  //!           members (see also TMeasurement):
  //!           - x->pt   : et  of whole tower
  //!           - x->EMF  : et  of ECAL  part
  //!           - x->HadF : et  of HCAL  part
  //!           - x->OutF : et  of Outer part
  //!           - x->E    : en  of Outer part
  //!           - x->eta  : eta of tower
  //!           - x->phi  : phi of tower
  //!  \param par Parameters of the correction function of this tower
  //!  \return The corrected Et of a tower 
  // -----------------------------------------------------------------
  virtual double correctedTowerEt(const TMeasurement *x,const double *par) const = 0;


  //!  \brief Corrects the measured jet Et
  //!
  //!  The parameters of the jet correction function
  //!  may depend on the jet Et, eta, and phi.
  //!
  //!  \param x Jet Et measurement that is to be corrected.
  //!           The measurement contains the following
  //!           members (see also TMeasurement):
  //!           - x->pt   : et  of whole jet
  //!           - x->EMF  : et  of ECAL  part
  //!           - x->HadF : et  of HCAL  part
  //!           - x->OutF : et  of Outer part
  //!           - x->E    : en  of Outer part
  //!           - x->eta  : eta of jet
  //!           - x->phi  : phi of jet
  //!  \param par Parameters of the correction function of this jet
  //!  \return The corrected Et of a jet
  // -----------------------------------------------------------------
  virtual double correctedJetEt(const TMeasurement *x,const double *par) const = 0;

 
  //!  \brief Returns the expected signal of a track in the Calorimeter
  //!
  //!  \param x Track Et measurement from which the expected response
  //!           of the calorimeter is calculated (see also TTrack).
  //!  \param par Parameters of the response function of this track
  //!  \return The expected calorimeter response
  // -----------------------------------------------------------------
  virtual double GetExpectedResponse(const TMeasurement *x,const double *par) const { return x->pt;}


  //!  \brief Corrects the measured jet Et with global
  //!         correction function
  //!
  //!  The parameters of the global jet correction function
  //!  are independent of the jet Et, eta, and phi.
  //!
  //!  \param x Jet Et measurement that is to be corrected.
  //!           The measurement contains the following
  //!           members (see also TMeasurement):
  //!           - x->pt   : et  of whole jet
  //!           - x->EMF  : et  of ECAL  part
  //!           - x->HadF : et  of HCAL  part
  //!           - x->OutF : et  of Outer part
  //!           - x->E    : en  of Outer part
  //!           - x->eta  : eta of jet
  //!           - x->phi  : phi of jet
  //!  \param par Parameters of the global correction function
  //!  \return The corrected Et of a jet
  // -----------------------------------------------------------------
  virtual double correctedGlobalJetEt(const TMeasurement *x,const double *par) const { return x->pt;}


  //!  \brief Get the name of the parametrization class
  //!  \return Name of the parametrization class
  // -----------------------------------------------------------------
  virtual const char * name() const = 0;


  //!  \brief Get the number of parameters of the tower parametrization
  //!  \return Number of parameters of the tower parametrization
  // -----------------------------------------------------------------
  unsigned int nTowerPars() const { return ntowerpars_;}


  //!  \brief Get the number of parameters of the jet parametrization
  //!  \return Number of parameters of the jet parametrization
  // -----------------------------------------------------------------
  unsigned int nJetPars() const { return njetpars_;}


  //!  \brief Get the number of parameters of the track parametrization
  //!  \return Number of parameters of the track parametrization
  // -----------------------------------------------------------------
  unsigned int nTrackPars() const { return ntrackpars_;}


  //!  \brief Get the number of parameters of the global jet parametrization
  //!  \return Number of parameters of the global jet parametrization
  // -----------------------------------------------------------------
  unsigned int nGlobalJetPars() const { return nglobaljetpars_;}


private: 
  Parametrization();
  unsigned int ntowerpars_;      //!< Number of parameters of the tower parametrization
  unsigned int njetpars_;        //!< Number of parameters of the jet parametrization
  unsigned int ntrackpars_;      //!< Number of parameters of the track parametrization
  unsigned int nglobaljetpars_;  //!< Number of parameters of the global jet parametrization
};



//!  \brief Parametrization of the hadronic tower response 
//!         by a step function in Et
//!
//!  The total jet measurement remains unchanged.
//!
//!  This parametrization has 12 tower parameters.
//!
//!  \sa Parametrization
// -----------------------------------------------------------------
class StepParametrization : public Parametrization { 
public:
  StepParametrization() : Parametrization(12,0,0,0) {}
  const char* name() const { return "StepParametrization";}
  
  double correctedTowerEt(const TMeasurement *x,const double *par) const {
    double result = 0;
    
    if(x->HadF>=0.0  && x->HadF<=1.0)  result = x->EMF+x->OutF + par[0]*x->HadF; 
    else if (x->HadF>   1.0 && x->HadF<=   2.0) result = x->EMF+x->OutF+par[ 1]*x->HadF;
    else if (x->HadF>   2.0 && x->HadF<=   5.0) result = x->EMF+x->OutF+par[ 2]*x->HadF;
    else if (x->HadF>   5.0 && x->HadF<=  10.0) result = x->EMF+x->OutF+par[ 3]*x->HadF;
    else if (x->HadF>  10.0 && x->HadF<=  20.0) result = x->EMF+x->OutF+par[ 4]*x->HadF;
    else if (x->HadF>  20.0 && x->HadF<=  40.0) result = x->EMF+x->OutF+par[ 5]*x->HadF;
    else if (x->HadF>  40.0 && x->HadF<=  80.0) result = x->EMF+x->OutF+par[ 6]*x->HadF;
    else if (x->HadF>  80.0 && x->HadF<= 160.0) result = x->EMF+x->OutF+par[ 7]*x->HadF;
    else if (x->HadF> 160.0 && x->HadF<= 300.0) result = x->EMF+x->OutF+par[ 8]*x->HadF;
    else if (x->HadF> 300.0 && x->HadF<= 600.0) result = x->EMF+x->OutF+par[ 9]*x->HadF;
    else if (x->HadF> 600.0 && x->HadF<=1000.0) result = x->EMF+x->OutF+par[10]*x->HadF;
    else if (x->HadF>1000.0 )                   result = x->EMF+x->OutF+par[11]*x->HadF;
    return result;
  }
    
  double correctedJetEt(const TMeasurement *x,const double *par) const {
    return  x->pt;   
    //OutOfCone, Dominant, parametrized in Et since cone R lorenz invariant
    //return x->pt * ( 1. + 0.295 * par[0] * exp(- 0.02566 * par[1] * x->pt)); 
    /*
      double result = 0;
      if(x->pt>=0.0  && x->pt<=5.0)          result =  par[0]*x->pt + par[1];
      else if (x->pt>5.0   && x->pt<=20.0)   result =  par[2]*x->pt + par[3];
      else if (x->pt>20.0  && x->pt<=80.0)   result =  par[4]*x->pt + par[5];
      else if (x->pt>80.0 )                 result =  par[6]*x->pt + par[7];
      return result;
    */
  }
};



//!  \brief Parametrization of the hadronic tower response 
//!         by a step function in E
//!
//!  Additionally, the total jet measurement is corrected
//!  by an Et-dependent function.
//!
//!  This parametrization has 12 tower parameters and 2
//!  jet parameters.
//!
//!  \sa Parametrization
// -----------------------------------------------------------------
class StepParametrizationEnergy : public Parametrization { 
public:
  StepParametrizationEnergy() : Parametrization(12,2,0,0) {}
  const char* name() const { return "StepParametrizationEnergy";}
  
  double correctedTowerEt(const TMeasurement *x,const double *par) const {
    double result = 0;
    // reweight from et to en
    double e =  x->HadF * x->E / x->pt;
    
    if(e>=0.0  && e<=1.0)  result = x->EMF+x->OutF + par[0]*x->HadF;
    else if (e>   1.0 && e<=   2.0) result = x->EMF+x->OutF+par[ 1]*x->HadF;
    else if (e>   2.0 && e<=   5.0) result = x->EMF+x->OutF+par[ 2]*x->HadF;
    else if (e>   5.0 && e<=  10.0) result = x->EMF+x->OutF+par[ 3]*x->HadF;
    else if (e>  10.0 && e<=  20.0) result = x->EMF+x->OutF+par[ 4]*x->HadF;
    else if (e>  20.0 && e<=  40.0) result = x->EMF+x->OutF+par[ 5]*x->HadF;
    else if (e>  40.0 && e<=  80.0) result = x->EMF+x->OutF+par[ 6]*x->HadF;
    else if (e>  80.0 && e<= 160.0) result = x->EMF+x->OutF+par[ 7]*x->HadF;
    else if (e> 160.0 && e<= 300.0) result = x->EMF+x->OutF+par[ 8]*x->HadF;
    else if (e> 300.0 && e<= 600.0) result = x->EMF+x->OutF+par[ 9]*x->HadF;
    else if (e> 600.0 && e<=1000.0) result = x->EMF+x->OutF+par[10]*x->HadF;
    else if (e>1000.0 )             result = x->EMF+x->OutF+par[11]*x->HadF;
    
    return result;
  }
    
  double correctedJetEt(const TMeasurement *x,const double *par) const {
    //OutOfCone, Dominant, parametrized in Et since cone R lorenz invariant
    return x->pt * ( 1. + 0.295 * par[0] * exp(- 0.02566 * par[1] * x->pt));   
  }
};



//!  \brief Parametrization of the hadronic tower
//!         response by a step function in EMF
//!
//!  Parametrization of hadronic response of a tower
//!  by a step function with 3 sets of parametrizations for 
//!  different em fractions.
//!
//!  This parametrization has 36 tower parameters.
//!
//!  The total jet measurement remains unchanged.
//!
//!  \sa Parametrization
// -----------------------------------------------------------------
class StepEfracParametrization : public Parametrization {
public:
  StepEfracParametrization() : Parametrization(36,0,0,0) {}  //(36,2) {}
  const char* name() const { return "StepEfracParametrization";}

  double correctedTowerEt(const TMeasurement *x,const double *par) const {
    double result=0;
    
    double Efrac = x->EMF/(x->HadF+x->OutF+x->EMF);
    if( Efrac < 0.2 ) {
      if(x->HadF>=0.0 && x->HadF<=1.0)            result = x->EMF+x->OutF+par[ 0]*x->HadF;
      else if (x->HadF>   1.0 && x->HadF<=   2.0) result = x->EMF+x->OutF+par[ 1]*x->HadF;
      else if (x->HadF>   2.0 && x->HadF<=   5.0) result = x->EMF+x->OutF+par[ 2]*x->HadF;
      else if (x->HadF>   5.0 && x->HadF<=  10.0) result = x->EMF+x->OutF+par[ 3]*x->HadF;
      else if (x->HadF>  10.0 && x->HadF<=  20.0) result = x->EMF+x->OutF+par[ 4]*x->HadF;
      else if (x->HadF>  20.0 && x->HadF<=  40.0) result = x->EMF+x->OutF+par[ 5]*x->HadF;
      else if (x->HadF>  40.0 && x->HadF<=  80.0) result = x->EMF+x->OutF+par[ 6]*x->HadF;
      else if (x->HadF>  80.0 && x->HadF<= 160.0) result = x->EMF+x->OutF+par[ 7]*x->HadF;
      else if (x->HadF> 160.0 && x->HadF<= 300.0) result = x->EMF+x->OutF+par[ 8]*x->HadF;
      else if (x->HadF> 300.0 && x->HadF<= 600.0) result = x->EMF+x->OutF+par[ 9]*x->HadF;
      else if (x->HadF> 600.0 && x->HadF<=1000.0) result = x->EMF+x->OutF+par[10]*x->HadF;
      else if (x->HadF>1000.0 )                   result = x->EMF+x->OutF+par[11]*x->HadF;
    } else if (Efrac < 0.5) {
      if(x->HadF>=0.0 && x->HadF<=1.0)            result = x->EMF+x->OutF+par[12]*x->HadF;
      else if (x->HadF>   1.0 && x->HadF<=   2.0) result = x->EMF+x->OutF+par[13]*x->HadF;
      else if (x->HadF>   2.0 && x->HadF<=   5.0) result = x->EMF+x->OutF+par[14]*x->HadF;
      else if (x->HadF>   5.0 && x->HadF<=  10.0) result = x->EMF+x->OutF+par[15]*x->HadF;
      else if (x->HadF>  10.0 && x->HadF<=  20.0) result = x->EMF+x->OutF+par[16]*x->HadF;
      else if (x->HadF>  20.0 && x->HadF<=  40.0) result = x->EMF+x->OutF+par[17]*x->HadF;
      else if (x->HadF>  40.0 && x->HadF<=  80.0) result = x->EMF+x->OutF+par[18]*x->HadF;
      else if (x->HadF>  80.0 && x->HadF<= 160.0) result = x->EMF+x->OutF+par[19]*x->HadF;
      else if (x->HadF> 160.0 && x->HadF<= 300.0) result = x->EMF+x->OutF+par[20]*x->HadF;
      else if (x->HadF> 300.0 && x->HadF<= 600.0) result = x->EMF+x->OutF+par[21]*x->HadF;
      else if (x->HadF> 600.0 && x->HadF<=1000.0) result = x->EMF+x->OutF+par[22]*x->HadF;
      else if (x->HadF>1000.0 )                   result = x->EMF+x->OutF+par[23]*x->HadF;
    } else {
      if(x->HadF>=0.0 && x->HadF<=1.0)            result = x->EMF+x->OutF+par[24]*x->HadF;
      else if (x->HadF>   1.0 && x->HadF<=   2.0) result = x->EMF+x->OutF+par[25]*x->HadF;
      else if (x->HadF>   2.0 && x->HadF<=   5.0) result = x->EMF+x->OutF+par[26]*x->HadF;
      else if (x->HadF>   5.0 && x->HadF<=  10.0) result = x->EMF+x->OutF+par[27]*x->HadF;
      else if (x->HadF>  10.0 && x->HadF<=  20.0) result = x->EMF+x->OutF+par[28]*x->HadF;
      else if (x->HadF>  20.0 && x->HadF<=  40.0) result = x->EMF+x->OutF+par[29]*x->HadF;
      else if (x->HadF>  40.0 && x->HadF<=  80.0) result = x->EMF+x->OutF+par[30]*x->HadF;
      else if (x->HadF>  80.0 && x->HadF<= 160.0) result = x->EMF+x->OutF+par[31]*x->HadF;
      else if (x->HadF> 160.0 && x->HadF<= 300.0) result = x->EMF+x->OutF+par[32]*x->HadF;
      else if (x->HadF> 300.0 && x->HadF<= 600.0) result = x->EMF+x->OutF+par[33]*x->HadF;
      else if (x->HadF> 600.0 && x->HadF<=1000.0) result = x->EMF+x->OutF+par[34]*x->HadF;
      else if (x->HadF>1000.0 )                   result = x->EMF+x->OutF+par[35]*x->HadF;
    }
    return result;
  }
  
  double correctedJetEt(const TMeasurement *x,const double *par) const {
    return  x->pt;   
    //OutOfCone, Dominant, parametrized in Et since cone R lorenz invariant
    //return x->pt * ( 1. + 0.295 * par[0] * exp(- 0.02566 * par[1] * x->pt));   
  }
};



//!  \brief Parametrization of the hadronic jet response 
//!         by a step function in EMF
//!
//!  Parametrization of hadronic response of a jet by a step 
//!  function with 3 sets of parametrizations for 
//!  different em fractions; the tower response remains
//!  unchanged.
//!
//!  This parametrization has 65 jet parameters.
//!
//!  \sa Parametrization
// -----------------------------------------------------------------
class StepJetParametrization : public Parametrization { 
public:
  StepJetParametrization() : Parametrization(0,65,0,0) {}
  const char* name() const { return "StepJetParametrization";}

  double correctedTowerEt(const TMeasurement *x,const double *par) const {
    return x->pt;
  }
    
  double correctedJetEt(const TMeasurement *x,const double *par) const {
    double pt     = x->pt;
    double Efrac  = x->EMF/( x->EMF+x->HadF+x->OutF);
    double result = 0.;
    if(Efrac < 0.2) {
      if(pt>=0.0 && pt<=10.0)         result = x->EMF+x->OutF+par[ 0]*x->HadF;
      else if (pt> 10.0 && pt<= 20.0) result = x->EMF+x->OutF+par[ 1]*x->HadF;
      else if (pt> 20.0 && pt<= 30.0) result = x->EMF+x->OutF+par[ 2]*x->HadF;
      else if (pt> 30.0 && pt<= 40.0) result = x->EMF+x->OutF+par[ 3]*x->HadF;
      else if (pt> 40.0 && pt<= 60.0) result = x->EMF+x->OutF+par[ 4]*x->HadF;
      else if (pt> 60.0 && pt<= 80.0) result = x->EMF+x->OutF+par[ 5]*x->HadF;
      else if (pt> 80.0 && pt<=100.0) result = x->EMF+x->OutF+par[ 6]*x->HadF;
      else if (pt>100.0 && pt<=120.0) result = x->EMF+x->OutF+par[ 7]*x->HadF;
      else if (pt>120.0 && pt<=140.0) result = x->EMF+x->OutF+par[ 8]*x->HadF;
      else if (pt>140.0 && pt<=160.0) result = x->EMF+x->OutF+par[ 9]*x->HadF;
      else if (pt>160.0 && pt<=180.0) result = x->EMF+x->OutF+par[10]*x->HadF;
      else if (pt>180.0 && pt<=200.0) result = x->EMF+x->OutF+par[11]*x->HadF;
      else if (pt>200.0 && pt<=225.0) result = x->EMF+x->OutF+par[12]*x->HadF;
      else if (pt>225.0 && pt<=250.0) result = x->EMF+x->OutF+par[13]*x->HadF;
      else if (pt>250.0 && pt<=275.0) result = x->EMF+x->OutF+par[14]*x->HadF;
      else if (pt>275.0 && pt<=300.0) result = x->EMF+x->OutF+par[15]*x->HadF;
      else if (pt>300.0 && pt<=350.0) result = x->EMF+x->OutF+par[16]*x->HadF;
      else if (pt>350.0 && pt<=400.0) result = x->EMF+x->OutF+par[17]*x->HadF;
      else if (pt>400.0 && pt<=500.0) result = x->EMF+x->OutF+par[18]*x->HadF;
      else if (pt>500.0 && pt<=700.0) result = x->EMF+x->OutF+par[19]*x->HadF;
      else if (pt>700.0 )             result = x->EMF+x->OutF+par[20]*x->HadF;
    } else if (Efrac < 0.5) {
      if(pt>=0.0 && pt<=10.0)         result = x->EMF+x->OutF+par[21]*x->HadF;
      else if (pt> 10.0 && pt<= 20.0) result = x->EMF+x->OutF+par[22]*x->HadF;
      else if (pt> 20.0 && pt<= 30.0) result = x->EMF+x->OutF+par[23]*x->HadF;
      else if (pt> 30.0 && pt<= 40.0) result = x->EMF+x->OutF+par[24]*x->HadF;
      else if (pt> 40.0 && pt<= 60.0) result = x->EMF+x->OutF+par[25]*x->HadF;
      else if (pt> 60.0 && pt<= 80.0) result = x->EMF+x->OutF+par[26]*x->HadF;
      else if (pt> 80.0 && pt<=100.0) result = x->EMF+x->OutF+par[27]*x->HadF;
      else if (pt>100.0 && pt<=120.0) result = x->EMF+x->OutF+par[28]*x->HadF;
      else if (pt>120.0 && pt<=140.0) result = x->EMF+x->OutF+par[29]*x->HadF;
      else if (pt>140.0 && pt<=160.0) result = x->EMF+x->OutF+par[30]*x->HadF;
      else if (pt>160.0 && pt<=180.0) result = x->EMF+x->OutF+par[31]*x->HadF;
      else if (pt>180.0 && pt<=200.0) result = x->EMF+x->OutF+par[32]*x->HadF;
      else if (pt>200.0 && pt<=225.0) result = x->EMF+x->OutF+par[33]*x->HadF;
      else if (pt>225.0 && pt<=250.0) result = x->EMF+x->OutF+par[34]*x->HadF;
      else if (pt>250.0 && pt<=275.0) result = x->EMF+x->OutF+par[35]*x->HadF;
      else if (pt>275.0 && pt<=300.0) result = x->EMF+x->OutF+par[36]*x->HadF;
      else if (pt>300.0 && pt<=350.0) result = x->EMF+x->OutF+par[37]*x->HadF;
      else if (pt>350.0 && pt<=400.0) result = x->EMF+x->OutF+par[38]*x->HadF;
      else if (pt>400.0 && pt<=500.0) result = x->EMF+x->OutF+par[39]*x->HadF;
      else if (pt>500.0 && pt<=700.0) result = x->EMF+x->OutF+par[40]*x->HadF;
      else if (pt>700.0 )             result = x->EMF+x->OutF+par[41]*x->HadF;
    } else {
      if(pt>=0.0 && pt<=10.0)   result = x->EMF+x->OutF+par[42]*x->HadF;
      else if (pt> 10.0 && pt<= 20.0) result = x->EMF+x->OutF+par[43]*x->HadF;
      else if (pt> 20.0 && pt<= 30.0) result = x->EMF+x->OutF+par[44]*x->HadF;
      else if (pt> 30.0 && pt<= 40.0) result = x->EMF+x->OutF+par[45]*x->HadF;
      else if (pt> 40.0 && pt<= 60.0) result = x->EMF+x->OutF+par[46]*x->HadF;
      else if (pt> 60.0 && pt<= 80.0) result = x->EMF+x->OutF+par[47]*x->HadF;
      else if (pt> 80.0 && pt<=100.0) result = x->EMF+x->OutF+par[48]*x->HadF;
      else if (pt>100.0 && pt<=120.0) result = x->EMF+x->OutF+par[49]*x->HadF;
      else if (pt>120.0 && pt<=140.0) result = x->EMF+x->OutF+par[50]*x->HadF;
      else if (pt>140.0 && pt<=160.0) result = x->EMF+x->OutF+par[51]*x->HadF;
      else if (pt>160.0 && pt<=180.0) result = x->EMF+x->OutF+par[52]*x->HadF;
      else if (pt>180.0 && pt<=200.0) result = x->EMF+x->OutF+par[53]*x->HadF;
      else if (pt>200.0 && pt<=225.0) result = x->EMF+x->OutF+par[54]*x->HadF;
      else if (pt>225.0 && pt<=250.0) result = x->EMF+x->OutF+par[55]*x->HadF;
      else if (pt>250.0 && pt<=275.0) result = x->EMF+x->OutF+par[56]*x->HadF;
      else if (pt>275.0 && pt<=300.0) result = x->EMF+x->OutF+par[57]*x->HadF;
      else if (pt>300.0 && pt<=350.0) result = x->EMF+x->OutF+par[58]*x->HadF;
      else if (pt>350.0 && pt<=400.0) result = x->EMF+x->OutF+par[59]*x->HadF;
      else if (pt>400.0 && pt<=500.0) result = x->EMF+x->OutF+par[60]*x->HadF;
      else if (pt>500.0 && pt<=700.0) result = x->EMF+x->OutF+par[61]*x->HadF;
      else if (pt>700.0 )             result = x->EMF+x->OutF+par[62]*x->HadF;
    }

    //OutOfCone, Dominant, parametrized in Et since cone R lorenz invariant
    result *= ( 1. + 0.295 * par[63] * exp(- 0.02566 * par[64] * result));   

    return  result;
  }
};



//!  \brief Parametrization of tower and jet response
//!         by some "clever" function
//!
//!
//!  This parametrization has 2 jet parameters.
//!
//!  \sa Parametrization
// -----------------------------------------------------------------
class MyParametrization: public Parametrization {
 public:
  MyParametrization() : Parametrization(0,2,0,0) {}
  const char* name() const { return "MyParametrization";}
  
  double correctedTowerEt(const TMeasurement *x,const double *par) const {
    return x->EMF + x->HadF + x->OutF;
    //return x->EMF + par[0]*x->HadF + par[1]*log(x->pt) + par[2];
  }
  
  double correctedJetEt(const TMeasurement *x,const double *par) const {
    return x->EMF + (par[0] + x->HadF * par[1]/1000) * x->HadF + x->OutF;;
  }
};



//!  \brief Parametrization of tower and jet response
//!         with some ideas from the JetMET group
//!
//!  This parametrization has 3 tower parameters and 5
//!  jet parameters.
//!
//!  \sa Parametrization
// -----------------------------------------------------------------
class JetMETParametrization: public Parametrization {
public:
  JetMETParametrization() : Parametrization(3,5,0,0) {}
  const char* name() const { return "JetMETParametrization";}
  
  double correctedTowerEt(const TMeasurement *x,const double *par) const {
    return par[1] * x->HadF + par[2] * x->EMF + x->OutF + par[0];
  }
  
  double correctedJetEt(const TMeasurement *x,const double *par) const {
    double logx = log(x->HadF);
    if(logx < 0) logx = 0;
    return (par[0] - par[1]/(pow(logx,par[2]) + par[3]) + par[4]/x->HadF) * x->HadF;  
  }
};



//!  \brief Simple tower and jet parametrization
//!
//!  This parametrization has 3 tower parameters and 3
//!  jet parameters.
//!
//!  \sa Parametrization
// -----------------------------------------------------------------
class SimpleParametrization: public Parametrization {
public:
  SimpleParametrization() : Parametrization(3,3,0,0) {}
  const char* name() const { return "SimpleParametrization";}

  double correctedTowerEt(const TMeasurement *x,const double *par) const {
    return par[1] * x->EMF + par[2] * x->HadF + x->OutF + par[0];
  }

  double correctedJetEt(const TMeasurement *x,const double *par) const {
    return x->pt * ( par[2] + par[0] * exp( -par[1] * x->pt ) );  
  }
};



//!  \brief Parametrization for toy MC with constant response
//!
//!  This is the parametrization of the correction for ToyMC
//!  events with a constant tower response i.e. when specifying
//!  only one tower constant. In this parametrization, the
//!  hadronic part of the tower Et is corrected.
//!
//!  This parametrization has 1 tower parameter.
//!
//!  \sa Parametrization, ToyMC
// -----------------------------------------------------------------
class ToyParametrization: public Parametrization {
public:
  ToyParametrization() : Parametrization(1,0,0,0) {}
  const char* name() const { return "ToyParametrization";}

  double correctedTowerEt(const TMeasurement *x,const double *par) const {
    return par[0] * x->HadF + x->EMF + x->OutF;
  }

  double correctedJetEt(const TMeasurement *x,const double *par) const {
    return x->pt;  
  }
};



//!  \brief Parametrization for toy MC with constant response
//!
//!  This is the parametrization of the correction for ToyMC
//!  events with a constant tower response i.e. when specifying
//!  only one tower constant.  In this parametrization, the
//!  hadronic part of the jet Et is corrected.
//!
//!  This parametrization has 1 jet parameters.
//!
//!  \sa Parametrization, ToyMC
// -----------------------------------------------------------------
class ToyJetParametrization: public Parametrization {
public:
  ToyJetParametrization() : Parametrization(0,1,0,0) {}
  const char* name() const { return "ToyJetParametrization";}

  double correctedTowerEt(const TMeasurement *x,const double *par) const {
    return x->pt;
  }

  double correctedJetEt(const TMeasurement *x,const double *par) const {
    return par[0] * x->HadF + x->EMF + x->OutF;
  }
};



//!  \brief Parametrization for toy MC with step function of
//!         hadronic tower Et
//!
//!  In this parametrization, the hadronic part of
//!  tower Et is corrected by a step function. It is
//!  intended for a pt spectrum from 0 - 300 GeV
//!
//!  This parametrization has 15 tower parameters.
//!
//!  \note This parametrization is intended for studying
//!        cutoff and resolution effects on the minimization
//!        procedure.
//!
//!  \sa Parametrization, ToyMC
// -----------------------------------------------------------------
class ToyStepParametrization : public Parametrization { 
public:
  ToyStepParametrization() : Parametrization(15,0,0,0) {}
  const char* name() const { return "ToyStepParametrization";}
  
  double correctedTowerEt(const TMeasurement *x,const double *par) const {
    double result = 0.;
    double pt = x->HadF;
    
    if(pt < 2.)        result = x->EMF+x->OutF+par[ 0]*x->HadF;
    else if(pt <   5.) result = x->EMF+x->OutF+par[ 1]*x->HadF;
    else if(pt <  10.) result = x->EMF+x->OutF+par[ 2]*x->HadF;
    else if(pt <  20.) result = x->EMF+x->OutF+par[ 3]*x->HadF;
    else if(pt <  30.) result = x->EMF+x->OutF+par[ 4]*x->HadF;
    else if(pt <  40.) result = x->EMF+x->OutF+par[ 5]*x->HadF;
    else if(pt <  50.) result = x->EMF+x->OutF+par[ 6]*x->HadF;
    else if(pt <  60.) result = x->EMF+x->OutF+par[ 7]*x->HadF;
    else if(pt <  70.) result = x->EMF+x->OutF+par[ 8]*x->HadF;
    else if(pt <  80.) result = x->EMF+x->OutF+par[ 9]*x->HadF;
    else if(pt <  90.) result = x->EMF+x->OutF+par[10]*x->HadF;
    else if(pt < 100.) result = x->EMF+x->OutF+par[11]*x->HadF;
    else if(pt < 110.) result = x->EMF+x->OutF+par[12]*x->HadF;
    else if(pt < 120.) result = x->EMF+x->OutF+par[13]*x->HadF;
    else               result = x->EMF+x->OutF+par[14]*x->HadF;
    return result;
  }
    
  double correctedJetEt(const TMeasurement *x,const double *par) const {
    return x->pt;
  }
};



//!  \brief Parametrization for toy MC with step function of
//!         hadronic jet Et
//!
//!  In this parametrization, the hadronic part of the
//!  jet Et is corrected by a step function. It is
//!  intended for a pt spectrum from 0 - 300 GeV
//!
//!  This parametrization has 15 jet parameters.
//!
//!  \note This parametrization is intended for studying
//!        cutoff and resolution effects on the minimization
//!        procedure.
//!
//!  \sa Parametrization, ToyMC
// -----------------------------------------------------------------
class ToyStepJetParametrization : public Parametrization { 
 public:
  ToyStepJetParametrization() : Parametrization(0,15,0,0) {}
  const char* name() const { return "ToyStepJetParametrization";}
  
  double correctedTowerEt(const TMeasurement *x,const double *par) const {
    return x->pt;
  }
  
  double correctedJetEt(const TMeasurement *x,const double *par) const {
    double pt = x->HadF;
    double result = 0.;
    
    if(pt < 25. )      result = x->EMF+x->OutF+par[ 0]*x->HadF;
    else if(pt <  50.) result = x->EMF+x->OutF+par[ 1]*x->HadF;
    else if(pt <  75.) result = x->EMF+x->OutF+par[ 2]*x->HadF;
    else if(pt < 100.) result = x->EMF+x->OutF+par[ 3]*x->HadF;
    else if(pt < 125.) result = x->EMF+x->OutF+par[ 4]*x->HadF;
    else if(pt < 150.) result = x->EMF+x->OutF+par[ 5]*x->HadF;
    else if(pt < 175.) result = x->EMF+x->OutF+par[ 6]*x->HadF;
    else if(pt < 200.) result = x->EMF+x->OutF+par[ 7]*x->HadF;
    else if(pt < 225.) result = x->EMF+x->OutF+par[ 8]*x->HadF;
    else if(pt < 250.) result = x->EMF+x->OutF+par[ 9]*x->HadF;
    else if(pt < 275.) result = x->EMF+x->OutF+par[10]*x->HadF;
    else if(pt < 300.) result = x->EMF+x->OutF+par[11]*x->HadF;
    else if(pt < 325.) result = x->EMF+x->OutF+par[12]*x->HadF;
    else if(pt < 350.) result = x->EMF+x->OutF+par[13]*x->HadF;
    else               result = x->EMF+x->OutF+par[14]*x->HadF;
    return  result;
  }
};



//!  \brief Complete Track Parametrization
//!
//!  Same parametrization as StepEfracParametrization,
//!  if track outside tracker or track errors too large
//!
//!  This parametrization has 12 tower parameters,
//!  3 jet parameters, and 6 track parameters.
//!
//!  \sa Parametrization
// -----------------------------------------------------------------
class TrackParametrization : public Parametrization {
 public:
  TrackParametrization() : Parametrization(12,3,6,0) {}  //(36,3,3,0) {}
  const char* name() const { return "TrackParametrization";}

  double correctedTowerEt(const TMeasurement *x,const double *par) const 
    {
      double result=0;
      
      //double Efrac = x->EMF/(x->HadF+x->OutF+x->EMF);
      //if( Efrac < 0.2 ) {
      if(x->HadF>=0.0 && x->HadF<=1.0)            result = x->EMF+x->OutF+par[ 0]*x->HadF;
      else if (x->HadF>   1.0 && x->HadF<=   2.0) result = x->EMF+x->OutF+par[ 1]*x->HadF;
      else if (x->HadF>   2.0 && x->HadF<=   5.0) result = x->EMF+x->OutF+par[ 2]*x->HadF;
      else if (x->HadF>   5.0 && x->HadF<=  10.0) result = x->EMF+x->OutF+par[ 3]*x->HadF;
      else if (x->HadF>  10.0 && x->HadF<=  20.0) result = x->EMF+x->OutF+par[ 4]*x->HadF;
      else if (x->HadF>  20.0 && x->HadF<=  40.0) result = x->EMF+x->OutF+par[ 5]*x->HadF;
      else if (x->HadF>  40.0 && x->HadF<=  80.0) result = x->EMF+x->OutF+par[ 6]*x->HadF;
      else if (x->HadF>  80.0 && x->HadF<= 160.0) result = x->EMF+x->OutF+par[ 7]*x->HadF;
      else if (x->HadF> 160.0 && x->HadF<= 300.0) result = x->EMF+x->OutF+par[ 8]*x->HadF;
      else if (x->HadF> 300.0 && x->HadF<= 600.0) result = x->EMF+x->OutF+par[ 9]*x->HadF;
      else if (x->HadF> 600.0 && x->HadF<=1000.0) result = x->EMF+x->OutF+par[10]*x->HadF;
      else if (x->HadF>1000.0 )                   result = x->EMF+x->OutF+par[11]*x->HadF;

      if(result < x->EMF+x->OutF) result = x->EMF+x->OutF;
      return result;
    }
  
  double correctedJetEt(const TMeasurement *x,const double *par) const 
    {
    double result;
    /*
    //result = x->pt * ( 1. + 0.295 * par[3] * exp(- 0.02566 * par[4] * x->pt));
    result = x->pt; 
    */ 
    if(x->E < -500) //set to -1000 or -800 for track jets! Positive for all others
      {
	if(x->E < 900)  //set to -1000 for Calo Rest
	  result = x->pt; //*par[0];
	else //set to -800 for complete Track Jet
	  result = x->pt * par[2];
      }
    else 
      result = par[1] * x->pt;
    if(result < 0)   result = 0; 
    return result;
  }
   
  double GetExpectedResponse(const TMeasurement *x,const double *par) const   
    {
    double result=0;
    double PiFrac;
    //Groom
    double eh = 1.48  * par[0]; 
    //Wigman
    //double eh = 1.39  * par[0]; 

    //double ehECAL = 1.6 ;//* par[1]; 
    double TrackEMF = 0;
    bool LS = false;
    
    TTrack* temp = (TTrack*)(x);
    
    if(temp->EM1 < 1.2) LS = true;   //late showering particle
    
    //this is Groom's Parametrization:
    TrackEMF = temp->EM1 / (temp->Had1 + temp->EM1);  //bei 1X1 EMF ueberschaetzt, da shower schmaler, bei 3X3 zu viele andere Tracks
    if(TrackEMF > 1) TrackEMF = 1; 
    if(TrackEMF < 0) TrackEMF = 0;
    //JPT alg (2004/15):
    //TrackEMF = 0.4;  //mean value from Test Beam

    //Groom
    PiFrac = 1 - pow(((x->E + par[4]) /(fabs(par[2]) * 0.96) ),(par[3] * 0.816) - 1 );  
    //Wigman
    //PiFrac = 0.11*par[2] * log(par[3] * x->E); 

    if(PiFrac > 1)  PiFrac = 1;   //do not allow unphysical results
    if(PiFrac < 0)  PiFrac = 0;  //happens e.g. for TrackPt<2GeV
    if(eh < 0.1) eh=0.1; 

    
    double responseH = (1 + (eh - 1) * PiFrac) / eh;
    
    /*
    //double responseE = (1 + (ehECAL - 1) * PiFrac) / ehECAL;
    double responseE = 1;
    
    double resultE = x->pt * TrackEMF * responseE;
    //if(LS) resultE = temp->EM1;
    if((temp->EM5 > 0) && (temp->EM5  < resultE) )  resultE = temp->EM5;
    
    double resultH = x->pt * (1 - TrackEMF) * responseH;
    if((temp->Had5 > 0) && (temp->Had5  < resultH) )  resultH = temp->Had5;
    result = resultE + resultH;
    */

    result = x->pt * responseH; 
 
    return result;
  }
};



//!  \brief L2L3 JetMET parametrization
//!
//!  This parametrization uses the correction functions
//!  from the JetMET group:
//!  - There is no tower correction function
//!  - The jet Et is corrected eta-dependent with the
//!    L2 correction function
//!  - The jet Et is corrected globally with the L3
//!    correction function
//!
//!  This parametrization has 3 jet parameters and
//!  4 global jet parameters.
//!
//!  \sa Parametrization
// -----------------------------------------------------------------
class L2L3JetParametrization : public Parametrization { 
public:
  L2L3JetParametrization() : Parametrization(0,3,0,4) {}
  const char* name() const { return "L2L3JetParametrization";}
  
  double correctedTowerEt(const TMeasurement *x,const double *par) const {
    return x->pt;
  }

  //!  \brief Code from L2RelativeCorrector
  //!  \code
  //!  double pt = (fPt < p[0]) ? p[0] : (fPt > p[1]) ? p[1] : fPt;
  //!  double logpt = log10(pt);
  //!  double result = p[2]+logpt*(p[3]+logpt*(p[4]+logpt*(p[5]+logpt*(p[6]+logpt*p[7]))));
  //!  \endcode   
  double correctedJetEt(const TMeasurement *x,const double *par) const {
    double pt = (x->pt < 4.0) ? 4.0 : (x->pt > 2000.0) ? 2000.0 : x->pt; 
    double logpt = log10(pt);
    //double result = par[0]+logpt*(par[1]+logpt*(par[2]+logpt*(par[3]+logpt*(par[4]+logpt*par[5]))));
    double c1 = par[0]+logpt*(0.1 * par[1]+logpt * 0.1* par[2]);
    return c1 * x->pt;
  }

  //!  \brief Code from SimpleL3AbsoluteCorrector
  //!  \code
  //!  double pt = (fPt < p[0]) ? p[0] : (fPt > p[1]) ? p[1] : fPt;
  //!  double log10pt = log10(pt);
  //!  double result = p[2]+p[3]/(pow(log10pt,p[4])+p[5]);
  //!  \endcode
  double correctedGlobalJetEt(const TMeasurement *x,const double *par) const {
    double pt = (x->pt < 4.0) ? 4.0 : (x->pt > 2000.0) ? 2000.0 : x->pt; 
    double logpt = log10(pt);
    double c2 = par[0] + par[1]/(pow(logpt,par[2]) + par[3]);
    return  c2 * x->pt; 
  }
};



//!  \brief L2L3 JetMET parametrization
//!
//!  This parametrization uses the correction functions
//!  from the JetMET group:
//!  - There is no tower correction function
//!  - The jet Et is corrected eta-dependent with the
//!    product of the L2 and L3 correction functions.
//!
//!  This parametrization has 7 jet parameters and
//!
//!  \sa Parametrization
// -----------------------------------------------------------------
class L2L3JetParametrization2 : public Parametrization { 
public:
  L2L3JetParametrization2() : Parametrization(0,7,0,0) {}
  const char* name() const { return "L2L3JetParametrization";}
  
  double correctedTowerEt(const TMeasurement *x,const double *par) const {
    return x->pt;
  }
    
  double correctedJetEt(const TMeasurement *x,const double *par) const {
    //code from L2RelativeCorrector
    //double pt = (fPt < p[0]) ? p[0] : (fPt > p[1]) ? p[1] : fPt;
    //double logpt = log10(pt);
    //double result = p[2]+logpt*(p[3]+logpt*(p[4]+logpt*(p[5]+logpt*(p[6]+logpt*p[7]))));
    double pt = (x->pt < 4.0) ? 4.0 : (x->pt > 2000.0) ? 2000.0 : x->pt; 
    double logpt = log10(pt);
    //double result = par[0]+logpt*(par[1]+logpt*(par[2]+logpt*(par[3]+logpt*(par[4]+logpt*par[5]))));
    double c1 = par[0]+logpt*(par[1]*0.01+logpt*0.001*par[2]);
    pt = (c1 * x->pt < 4.0) ? 4.0 : (c1 * x->pt > 2000.0) ? 2000.0 : c1 * x->pt; 
    logpt = log10(pt);
    //result = par[6] + par[7]/(pow(logpt,par[8]) + par[9]);
    double c2 = par[3] + par[4]/(pow(logpt,par[5]) + par[6]);
    //std::cout << par[0] << ", " << par[1] << ", " << par[2] << ", " << par[3] << ", " 
    //	      <<  par[4] << " c2:" << c2 << " Et:" << x->pt << '\n';
    return  c2 * c1 * x->pt; 
  }
};



//!  \brief Complete Track Parametrization
//!
//!  Same parametrization as L2L3JetParametrization,
//!  if track outside tracker or track errors too large
//!
//!  This parametrization has 3 jet parameters, 5 track
//!  parameters, and 4 global jet parameters.
//!
//!  \sa Parametrization
// -----------------------------------------------------------------

class L2L3JetTrackParametrization : public Parametrization { 
public:
  L2L3JetTrackParametrization() : Parametrization(0,3,5,4) {}
  const char* name() const { return "L2L3JetTrackParametrization";}
  
  double correctedTowerEt(const TMeasurement *x,const double *par) const {
    return x->pt;
  }
    
  double correctedJetEt(const TMeasurement *x,const double *par) const {
    //code from L2RelativeCorrector
    //double pt = (fPt < p[0]) ? p[0] : (fPt > p[1]) ? p[1] : fPt;
    //double logpt = log10(pt);
    //double result = p[2]+logpt*(p[3]+logpt*(p[4]+logpt*(p[5]+logpt*(p[6]+logpt*p[7]))));
    double pt = (x->pt < 4.0) ? 4.0 : (x->pt > 2000.0) ? 2000.0 : x->pt; 
    double logpt = log10(pt);
    //double result = par[0]+logpt*(par[1]+logpt*(par[2]+logpt*(par[3]+logpt*(par[4]+logpt*par[5]))));
    double c1 = par[0]+logpt*(0.1 * par[1]+logpt * 0.1* par[2]);
    return c1 * x->pt;
  }

  double correctedGlobalJetEt(const TMeasurement *x,const double *par) const {
    // code from SimpleL3AbsoluteCorrector
    //double pt = (fPt < p[0]) ? p[0] : (fPt > p[1]) ? p[1] : fPt;
    //double log10pt = log10(pt);
    //double result = p[2]+p[3]/(pow(log10pt,p[4])+p[5]);
    double pt = (x->pt < 4.0) ? 4.0 : (x->pt > 2000.0) ? 2000.0 : x->pt; 
    double logpt = log10(pt);
    //result = par[6] + par[7]/(pow(logpt,par[8]) + par[9]);
    double c2 = par[0] + par[1]/(pow(logpt,par[2]) + par[3]);
    if(c2 < par[0]) c2 = par[0];
    //std::cout << par[0] << ", " << par[1] << ", " << par[2] << ", " << par[3] << ", " 
    //	      << " c2:" << c2 << " Et:" << x->pt << '\n';
    return  c2 * x->pt; 
  }
  
  double GetExpectedResponse(const TMeasurement *x,const double *par) const {
    double result=0;
    //Groom
    double eh = 1.48 * par[0]; 
    //Wigman
    //double eh = 1.39  * par[0]; 

    //double ehECAL = 1.6 ;//* par[1]; 
    double TrackEMF = 0;
    bool LS = false;
    
    TTrack* temp = (TTrack*)(x);
    
    if(temp->EM1 < 1.2) LS = true;   //late showering particle
    
    //this is Groom's Parametrization:
    TrackEMF = temp->EM1 / (temp->Had1 + temp->EM1);  //bei 1X1 EMF ueberschaetzt, da shower schmaler, bei 3X3 zu viele andere Tracks
    if(TrackEMF > 1) TrackEMF = 1; 
    if(TrackEMF < 0) TrackEMF = 0;
    //JPT alg (2004/15):
    //TrackEMF = 0.4;  //mean value from Test Beam
    
    //Groom
    double Ecor = (x->E + par[4] * 0.01) /par[2] /0.96;
    double PiFrac = Ecor > 0 ? 1 - pow(Ecor,par[3] * 0.816 - 1 ) : 1.0;  
    //PiFrac = 1 - pow(Ecor,-0.16);
    //Wigman
    //PiFrac = 0.11*par[2] * log(par[3] * x->E); 
    assert(PiFrac == PiFrac);
    if(PiFrac > 1)  PiFrac = 1;   //do not allow unphysical results
    if(PiFrac < 0)  PiFrac = 0;  //happens e.g. for TrackPt<2GeV
    if(eh < 0.1) eh=0.1; 
    
    
    double responseH = (1 + (eh - 1) * PiFrac) / eh;
    
    /*
    //double responseE = (1 + (ehECAL - 1) * PiFrac) / ehECAL;
    double responseE = 1;
    
    double resultE = x->pt * TrackEMF * responseE;
    //if(LS) resultE = temp->EM1;
    if((temp->EM5 > 0) && (temp->EM5  < resultE) )  resultE = temp->EM5;
    
    double resultH = x->pt * (1 - TrackEMF) * responseH;
    if((temp->Had5 > 0) && (temp->Had5  < resultH) )  resultH = temp->Had5;
    result = resultE + resultH;
    */

    result = x->pt * responseH; 

    return result;
  }
};


//!  \brief Parametrization for Toy MC with ToyMC::Response SimpleInverse
//!
//!  This parametrization contains the correction function
//!  for the ToyMC::Response SimpleInverse
//!  \f[ C(E_{T}) = \sqrt{ A^{2} + A_{1}E_{T}} - A \f]
//!  where \f$ A = -0.5(A_{1} - A_{0} - E_{T}) \f$.
//!  This function was derived by analytical inversion.
//!  - There is no tower correction function
//!  - There is no jet Et correction function
//!  - The hadronic part of the jet Et is corrected
//!    globally with the correction function given above
//!
//!  This parametrization 2 global jet parameters.
//!
//!  \sa Parametrization
// -----------------------------------------------------------------
class ToySimpleInverseParametrization : public Parametrization { 
public:
  ToySimpleInverseParametrization() : Parametrization(0,0,0,2) {}
  const char* name() const { return "ToySimpleInverseParametrization";}
  
  double correctedTowerEt(const TMeasurement *x,const double *par) const {
    return x->pt;
  }

  double correctedJetEt(const TMeasurement *x,const double *par) const {
    return x->pt;
  }

  double correctedGlobalJetEt(const TMeasurement *x,const double *par) const {
    double a = 0.5*(par[1] - par[0] - x->HadF);
    return x->EMF + x->OutF + sqrt( a*a + par[1]*x->HadF ) - a;
  }
};


#endif
