//
// Original Author:  Hartmut Stadie
//         Created:  Thu Apr 03 17:09:50 CEST 2008
// $Id: Parametrization.h,v 1.8 2008/07/31 12:51:15 thomsen Exp $
//
#ifndef CALIBCORE_PARAMETRIZATION_H
#define CALIBCORE_PARAMETRIZATION_H

#include <cmath>
#include "CalibData.h"

class Parametrization 
{
public:
  Parametrization(unsigned int ntowerpars, unsigned int njetpars) 
    : ntowerpars_(ntowerpars), njetpars_(njetpars) {}
  virtual ~Parametrization() {}

  /** correctedTowerEt(double *x,double *par)
      returns the corrected Et of a tower 
           input:            x->pt : Et of whole tower
                             x->EMF : Et of ECAL part
                             x->HadF : Et of HCAL part
                             x->OutF : Et of Outer part
           par:  the correction parameters of this tower
  **/
  virtual double correctedTowerEt(TMeasurement *x,double *par) const = 0;

  /** correctedJetEt(double *x,double *par)
      returns the corrected Et of a jet 
      input: x->OutF: x->pt : Et of uncorrected jet
                        x->EMF : eta of uncorrected jet(not used)
                        x->HadF : phi of uncorrected jet(not used)
           par:  the correction parameters of this jet
  **/
  virtual double correctedJetEt(TMeasurement *x,double *par) const = 0;
  virtual const char * name() const = 0;

  unsigned int nTowerPars() const { return ntowerpars_;}
  unsigned int nJetPars() const { return njetpars_;}

private: 
  Parametrization();
  unsigned int  ntowerpars_, njetpars_;
};

// Parametrization of hadronic response by a step function
class StepParametrization : public Parametrization { 
public:
  StepParametrization() : Parametrization(12,8) {}
  //StepParametrization() : Parametrization(12,24) {}
 
  const char* name() const { return "StepParametrization";}

  double correctedTowerEt(TMeasurement *x,double *par) const {
    double result = 0;
    
    if(x->HadF>=0.0  && x->HadF<=1.0)  result = x->EMF+x->OutF + par[0]*x->HadF;
    else if (x->HadF>1.0   && x->HadF<=2.0)  result = x->EMF+x->OutF + par[1]*x->HadF;
    else if (x->HadF>2.0   && x->HadF<=5.0)  result = x->EMF+x->OutF + par[2]*x->HadF;
    else if (x->HadF>5.0   && x->HadF<=10.0)  result = x->EMF+x->OutF + par[3]*x->HadF;
    else if (x->HadF>10.0  && x->HadF<=20.0)  result = x->EMF+x->OutF + par[4]*x->HadF;
    else if (x->HadF>20.0  && x->HadF<=40.0)  result = x->EMF+x->OutF + par[5]*x->HadF;
    else if (x->HadF>40.0  && x->HadF<=80.0) result = x->EMF+x->OutF + par[6]*x->HadF;
    else if (x->HadF>80.0  && x->HadF<=160.0) result = x->EMF+x->OutF + par[7]*x->HadF;
    else if (x->HadF>160.0 && x->HadF<=300.0) result = x->EMF+x->OutF + par[8]*x->HadF;
    else if (x->HadF>300.0 && x->HadF<=600.0) result = x->EMF+x->OutF + par[9]*x->HadF;
    else if (x->HadF>600.0 && x->HadF<=1000.0) result = x->EMF+x->OutF + par[10]*x->HadF;
    else if (x->HadF>1000.0 )              result = x->EMF+x->OutF + par[11]*x->HadF;
    return result;
  }
    
  double correctedJetEt(TMeasurement *x,double *par) const {
    double result = 0;
    /*
    if(x->pt>=0.0  && x->pt<=1.0)          result =  par[0]*x->pt + par[1];
    else if (x->pt>1.0   && x->pt<=2.0)    result =  par[2]*x->pt + par[3];
    else if (x->pt>2.0   && x->pt<=5.0)    result =  par[4]*x->pt + par[5];
    else if (x->pt>5.0   && x->pt<=10.0)   result =  par[6]*x->pt + par[7];
    else if (x->pt>10.0  && x->pt<=20.0)   result =  par[8]*x->pt + par[9];
    else if (x->pt>20.0  && x->pt<=40.0)   result =  par[10]*x->pt + par[11];
    else if (x->pt>40.0  && x->pt<=80.0)   result =  par[12]*x->pt + par[13];
    else if (x->pt>80.0  && x->pt<=160.0)  result =  par[14]*x->pt + par[15];
    else if (x->pt>160.0 && x->pt<=300.0)  result =  par[16]*x->pt + par[17];
    else if (x->pt>300.0 && x->pt<=600.0)  result =  par[18]*x->pt + par[19];
    else if (x->pt>600.0 && x->pt<=1000.0) result =  par[20]*x->pt + par[21];
    else if (x->pt>1000.0 )               result =  par[22]*x->pt + par[23];
    return result;
    */
    
    if(x->pt>=0.0  && x->pt<=5.0)          result =  par[0]*x->pt + par[1];
    else if (x->pt>5.0   && x->pt<=20.0)   result =  par[2]*x->pt + par[3];
    else if (x->pt>20.0  && x->pt<=80.0)   result =  par[4]*x->pt + par[5];
    else if (x->pt>80.0 )                 result =  par[6]*x->pt + par[7];
    return result;
 

    //return  x->pt;
    //return  par[0]*x->pt + par[1];
  }
};

// Parametrization of hadronic response by a step function
class StepParametrizationEnergy : public Parametrization { 
public:
  StepParametrizationEnergy() : Parametrization(14,2) {} //14,5
 
  const char* name() const { return "StepParametrizationEnergy";}

  double correctedTowerEt(TMeasurement *x,double *par) const {
    double result = 0;
    double e =  x->HadF * x->E / x->pt;

    if(e>=0.0  && e<=1.0)  result = x->EMF+x->OutF + par[0]*x->HadF;
    else if (e>1.0   && e<=2.0)  result = x->EMF+x->OutF + par[1]*x->HadF;
    else if (e>2.0   && e<=5.0)  result = x->EMF+x->OutF + par[2]*x->HadF;
    else if (e>5.0   && e<=10.0)  result = x->EMF+x->OutF + par[3]*x->HadF;
    else if (e>10.0  && e<=20.0)  result = x->EMF+x->OutF + par[4]*x->HadF;
    else if (e>20.0  && e<=40.0)  result = x->EMF+x->OutF + par[5]*x->HadF;
    else if (e>40.0  && e<=80.0) result = x->EMF+x->OutF + par[6]*x->HadF;
    else if (e>80.0  && e<=160.0) result = x->EMF+x->OutF + par[7]*x->HadF;
    else if (e>160.0 && e<=300.0) result = x->EMF+x->OutF + par[8]*x->HadF;
    else if (e>300.0 && e<=600.0) result = x->EMF+x->OutF + par[9]*x->HadF;
    else if (e>600.0 && e<=1000.0) result = x->EMF+x->OutF + par[10]*x->HadF;
    else if (e>1000.0 )              result = x->EMF+x->OutF + par[11]*x->HadF;
    return result;
  }
    
  double correctedJetEt(TMeasurement *x,double *par) const {
    double result = x->pt * ( 1. + par[0] * exp(-x->pt));   //Out of Cone, Dominant, parametrized in Et since cone R lorenz invariant
    if(x->OutF > exp(par[1]))        //punch through?
      result = result * (1. + par[2] * (log(x->OutF) - par[1]));

    /*
    if(par[2]<0) par[2] = 0;
    if(par[3]<0) par[3] = 0;
    if(x[3] > exp(par[2]))        //punch through?
      result = result * (1 + par[3] * (log(x[3]) - par[2]));
    */
    //result -= par[4];
    //if(result<0)   std::cout<<"gfdgfdfghk"<<std::endl;  // 
    //result *= -1;
    return result;
  }
};

// Parametrization of hadronic response by a step function
// 3 Sets of Parametrization for different EM fraction

class StepEfracParametrization : public Parametrization {
public:
  StepEfracParametrization() : Parametrization(36,2) {}
  
  const char* name() const { return "StepEfracParametrization";}

  double correctedTowerEt(TMeasurement *x,double *par) const {
    double result=0;
    
    //double Efrac = x->EMF/(x->HadF+x->OutF);
    if( x->EMF < 0.1 * (x->HadF+x->OutF) ) {
      if      (x->HadF>=0.0   && x->HadF<=1.0)   result = x->EMF+x->OutF + par[0]*x->HadF;
      else if (x->HadF>1.0   && x->HadF<=2.0)    result = x->EMF+x->OutF + par[1]*x->HadF;
      else if (x->HadF>2.0   && x->HadF<=5.0)    result = x->EMF+x->OutF + par[2]*x->HadF;
      else if (x->HadF>5.0   && x->HadF<=10.0)   result = x->EMF+x->OutF + par[3]*x->HadF;
      else if (x->HadF>10.0  && x->HadF<=20.0)   result = x->EMF+x->OutF + par[4]*x->HadF;
      else if (x->HadF>20.0  && x->HadF<=40.0)   result = x->EMF+x->OutF + par[5]*x->HadF;
      else if (x->HadF>40.0  && x->HadF<=80.0)   result = x->EMF+x->OutF + par[6]*x->HadF;
      else if (x->HadF>80.0  && x->HadF<=160.0)  result = x->EMF+x->OutF + par[7]*x->HadF;
      else if (x->HadF>160.0 && x->HadF<=300.0)  result = x->EMF+x->OutF + par[8]*x->HadF;
      else if (x->HadF>300.0 && x->HadF<=600.0)  result = x->EMF+x->OutF + par[9]*x->HadF;
      else if (x->HadF>600.0 && x->HadF<=1000.0) result = x->EMF+x->OutF + par[10]*x->HadF;
      else if (x->HadF>1000.0 )               result = x->EMF+x->OutF + par[11]*x->HadF;
    } else if (x->EMF<0.3*(x->HadF+x->OutF)) {
      if      (x->HadF>=0.0   && x->HadF<=1.0)   result = x->EMF+x->OutF + par[12]*x->HadF;
      else if (x->HadF>1.0   && x->HadF<=2.0)    result = x->EMF+x->OutF + par[13]*x->HadF;
      else if (x->HadF>2.0   && x->HadF<=5.0)    result = x->EMF+x->OutF + par[14]*x->HadF;
      else if (x->HadF>5.0   && x->HadF<=10.0)   result = x->EMF+x->OutF + par[15]*x->HadF;
      else if (x->HadF>10.0  && x->HadF<=20.0)   result = x->EMF+x->OutF + par[16]*x->HadF;
      else if (x->HadF>20.0  && x->HadF<=40.0)   result = x->EMF+x->OutF + par[17]*x->HadF;
      else if (x->HadF>40.0  && x->HadF<=80.0)   result = x->EMF+x->OutF + par[18]*x->HadF;
      else if (x->HadF>80.0  && x->HadF<=160.0)  result = x->EMF+x->OutF + par[19]*x->HadF;
      else if (x->HadF>160.0 && x->HadF<=300.0)  result = x->EMF+x->OutF + par[20]*x->HadF;
      else if (x->HadF>300.0 && x->HadF<=600.0)  result = x->EMF+x->OutF + par[21]*x->HadF;
      else if (x->HadF>600.0 && x->HadF<=1000.0) result = x->EMF+x->OutF + par[22]*x->HadF;
      else if (x->HadF>1000.0 )               result = x->EMF+x->OutF + par[23]*x->HadF;
    } else {
      if      (x->HadF>=0.0   && x->HadF<=1.0)   result = x->EMF+x->OutF + par[24]*x->HadF;
      else if (x->HadF>1.0   && x->HadF<=2.0)    result = x->EMF+x->OutF + par[25]*x->HadF;
      else if (x->HadF>2.0   && x->HadF<=5.0)    result = x->EMF+x->OutF + par[26]*x->HadF;
      else if (x->HadF>5.0   && x->HadF<=10.0)   result = x->EMF+x->OutF + par[27]*x->HadF;
      else if (x->HadF>10.0  && x->HadF<=20.0)   result = x->EMF+x->OutF + par[28]*x->HadF;
      else if (x->HadF>20.0  && x->HadF<=40.0)   result = x->EMF+x->OutF + par[29]*x->HadF;
      else if (x->HadF>40.0  && x->HadF<=80.0)   result = x->EMF+x->OutF + par[30]*x->HadF;
      else if (x->HadF>80.0  && x->HadF<=160.0)  result = x->EMF+x->OutF + par[31]*x->HadF;
      else if (x->HadF>160.0 && x->HadF<=300.0)  result = x->EMF+x->OutF + par[32]*x->HadF;
      else if (x->HadF>300.0 && x->HadF<=600.0)  result = x->EMF+x->OutF + par[33]*x->HadF;
      else if (x->HadF>600.0 && x->HadF<=1000.0) result = x->EMF+x->OutF + par[34]*x->HadF;
      else if (x->HadF>1000.0 )               result = x->EMF+x->OutF + par[35]*x->HadF;
    }
    return result;
  }
  
  double correctedJetEt(TMeasurement *x,double *par) const {
    return  par[0]*x->pt + par[1];
  }
};

// Parametrization of response by some "clever" function
class MyParametrization: public Parametrization {
public:
  MyParametrization() : Parametrization(3,2) {}
  const char* name() const { return "MyParametrization";}
  double correctedTowerEt(TMeasurement *x,double *par) const {
    return x->EMF + par[0]*x->HadF + par[1]*log(x->pt) + par[2];
  }
  double correctedJetEt(TMeasurement *x,double *par) const {
    return par[0]*x->pt + par[1];
  }
};

// Parametrization of response with some ideas from the JetMET group
class JetMETParametrization: public Parametrization {
public:
  JetMETParametrization() : Parametrization(3,5) {}
  const char* name() const { return "JetMETParametrization";}
  double correctedTowerEt(TMeasurement *x,double *par) const {
    if(par[0] < -10) par[0] = -10;
    if(par[1] < 0) par[1] = -par[1];
    if(par[2] < 0) par[2] = -par[2];
    return par[1] * x->HadF + par[2] * x->EMF + x->OutF + par[0];
  }
  double correctedJetEt(TMeasurement *x,double *par) const {
    double logx = log(x->pt);
    if(logx < 0) logx = 0;
    if(par[1] < 0) par[1] *= -1;
    if(par[2] < 0) par[2] *= -1;
    if(par[3] < 0) par[3] *= -1;
    if(par[4] < 0) par[4] *= -1;
    return (par[0] - par[1]/(pow(logx,par[2]) + par[3]) + par[4]/x->pt) * x->pt;  
  }
};

// Simple parametrization
class SimpleParametrization: public Parametrization {
public:
  SimpleParametrization() : Parametrization(3,3) {}
  const char* name() const { return "SimpleParametrization";}
  double correctedTowerEt(TMeasurement *x,double *par) const {
    if(par[0] < -10) par[0] = -10;
    if(par[1] < 0) par[1] = -par[1];
    if(par[2] < 0) par[2] = -par[2];
    return par[1] * x->EMF + par[2] * x->HadF + x->OutF + par[0];
  }
  double correctedJetEt(TMeasurement *x,double *par) const {
    if(par[0] < 0) par[0] *= -1;
    if(par[1] < 0) par[1] *= -1;
    if(par[2] < 0) par[2] *= -1;
    return x->pt * ( par[2] + par[0] * exp( -par[1] * x->pt ) );  
  }
};

// Parametrization for toy MC
class ToyParametrization: public Parametrization {
public:
  ToyParametrization() : Parametrization(3,0) {}
  const char* name() const { return "ToyParametrization";}
  double correctedTowerEt(TMeasurement *x,double *par) const {
    return par[0] * x->HadF + par[1] * x->EMF + x->OutF + par[2];
  }
  double correctedJetEt(TMeasurement *x,double *par) const {
    return x->pt;  
  }
};
#endif
