//
// Original Author:  Hartmut Stadie
//         Created:  Thu Apr 03 17:09:50 CEST 2008
// $Id: Parametrization.h,v 1.13 2008/09/19 14:00:39 thomsen Exp $
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
  virtual double correctedTowerEt(TMeasurement *const x,double *const par) const = 0;

  /** correctedJetEt(double *x,double *par)
      returns the corrected Et of a jet 
      input: x->OutF: x->pt : Et of uncorrected jet
                        x->EMF : eta of uncorrected jet(not used)
                        x->HadF : phi of uncorrected jet(not used)
           par:  the correction parameters of this jet
  **/
  virtual double correctedJetEt(TMeasurement *const x,double *const par) const = 0;
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
  StepParametrization() : Parametrization(12,0) {}
  //StepParametrization() : Parametrization(12,2) {}
 
  const char* name() const { return "StepParametrization";}

  double correctedTowerEt(TMeasurement *const x,double *par) const {
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
    
  double correctedJetEt(TMeasurement *const x,double *const par) const {
    double result = 0;
    return  x->pt;   
    //return x->pt * ( 1. + 0.295 * par[0] * exp(- 0.02566 * par[1] * x->pt));   //Out of Cone, Dominant, parametrized in Et since cone R lorenz invariant


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
    
    
    if(x->pt>=0.0  && x->pt<=5.0)          result =  par[0]*x->pt + par[1];
    else if (x->pt>5.0   && x->pt<=20.0)   result =  par[2]*x->pt + par[3];
    else if (x->pt>20.0  && x->pt<=80.0)   result =  par[4]*x->pt + par[5];
    else if (x->pt>80.0 )                 result =  par[6]*x->pt + par[7];
    return result;
 
    */
  }
};

// Parametrization of hadronic response by a step function
class StepParametrizationEnergy : public Parametrization { 
public:
  StepParametrizationEnergy() : Parametrization(13,2) {} //14,5
 
  const char* name() const { return "StepParametrizationEnergy";}

  double correctedTowerEt(TMeasurement *const x,double *const par) const {
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
    //result += par[12];
    return result;
  }
    
  double correctedJetEt(TMeasurement *const x,double *const par) const {
    return x->pt * ( 1. + 0.295 * par[0] * exp(- 0.02566 * par[1] * x->pt));   //Out of Cone, Dominant, parametrized in Et since cone R lorenz invariant
  }
};

// Parametrization of hadronic response by a step function
// 3 Sets of Parametrization for different EM fraction

class StepEfracParametrization : public Parametrization {
public:
  StepEfracParametrization() : Parametrization(36,0) {}  //(36,2) {}
  
  const char* name() const { return "StepEfracParametrization";}

  double correctedTowerEt(TMeasurement *const x,double *const par) const {
    double result=0;
    
    double Efrac = x->EMF/(x->HadF+x->OutF+x->EMF);
    if( Efrac < 0.2 ) {
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
    } else if (Efrac < 0.5)) {
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
  
  double correctedJetEt(TMeasurement *const x,double *const par) const {
    return  x->pt;   
    //return x->pt * ( 1. + 0.295 * par[0] * exp(- 0.02566 * par[1] * x->pt));   //Out of Cone, Dominant, parametrized in Et since cone R lorenz invariant
  }
};

// Parametrization of response by some "clever" function
class MyParametrization: public Parametrization {
public:
  MyParametrization() : Parametrization(3,2) {}
  const char* name() const { return "MyParametrization";}
  double correctedTowerEt(TMeasurement *const x,double *const par) const {
    return x->EMF + par[0]*x->HadF + par[1]*log(x->pt) + par[2];
  }
  double correctedJetEt(TMeasurement *const x,double *const par) const {
    return par[0]*x->pt + par[1];
  }
};

// Parametrization of response with some ideas from the JetMET group
class JetMETParametrization: public Parametrization {
public:
  JetMETParametrization() : Parametrization(3,5) {}
  const char* name() const { return "JetMETParametrization";}
  double correctedTowerEt(TMeasurement *const x,double *const par) const {
    if(par[0] < -10) par[0] = -10;
    if(par[1] < 0) par[1] = -par[1];
    if(par[2] < 0) par[2] = -par[2];
    return par[1] * x->HadF + par[2] * x->EMF + x->OutF + par[0];
  }
  double correctedJetEt(TMeasurement *const x,double *const par) const {
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
  double correctedTowerEt(TMeasurement *const x,double *const par) const {
    if(par[0] < -10) par[0] = -10;
    if(par[1] < 0) par[1] = -par[1];
    if(par[2] < 0) par[2] = -par[2];
    return par[1] * x->EMF + par[2] * x->HadF + x->OutF + par[0];
  }
  double correctedJetEt(TMeasurement *const x,double *const par) const {
    if(par[0] < 0) par[0] *= -1;
    if(par[1] < 0) par[1] *= -1;
    if(par[2] < 0) par[2] *= -1;
    return x->pt * ( par[2] + par[0] * exp( -par[1] * x->pt ) );  
  }
};

// Parametrization for toy MC
class ToyParametrization: public Parametrization {
public:
  ToyParametrization() : Parametrization(1,0) {}
  const char* name() const { return "ToyParametrization";}
  double correctedTowerEt(TMeasurement *const x,double *const par) const {
    //  if(fabs( x->HadF + x->EMF + x->OutF -   x->pt)  >  0.1)   std::cout<< x->HadF + x->EMF + x->OutF -   x->pt<<"   hjsdgfuyksdgf   "<< x->HadF + x->EMF + x->OutF -   x->E<<std::endl;
    //std::cout<<par[0]<<"       "<<par<<std::endl;
    return par[0] * x->HadF + x->EMF + x->OutF;
  }
  double correctedJetEt(TMeasurement *const x,double *const par) const {
    return x->pt;  
  }
};

// Parametrization of hadronic response by a step function
class ToyStepParametrizationEnergy : public Parametrization { 
public:
  ToyStepParametrizationEnergy() : Parametrization(80,0) {} //5,0
 
  const char* name() const { return "ToyStepParametrizationEnergy";}
  
  double correctedTowerEt(TMeasurement *const x,double *const par) const {
    double result = 0;
    double e = x->HadF ;//* x->E / x->pt;

    if(e<=5.0)        result =  x->EMF+x->OutF + par[0]*x->HadF;
    else if (e<=10.0) result =  x->EMF+x->OutF + par[1]*x->HadF;
    else if (e<=15.0) result =  x->EMF+x->OutF + par[2]*x->HadF;
    else if (e<=20)   result =  x->EMF+x->OutF + par[3]*x->HadF;
    else if (e<=25.0) result =  x->EMF+x->OutF + par[4]*x->HadF;
    else if (e<=30.0) result =  x->EMF+x->OutF + par[5]*x->HadF;
    else if (e<=35)   result =  x->EMF+x->OutF + par[6]*x->HadF;
    else if (e<=40.0) result =  x->EMF+x->OutF + par[7]*x->HadF;
    else if (e<=45.0) result =  x->EMF+x->OutF + par[8]*x->HadF;
    else if (e<=50.0) result =  x->EMF+x->OutF + par[9]*x->HadF;
    else if (e<=55.0) result =  x->EMF+x->OutF + par[10]*x->HadF;
    else if (e<=60)   result =  x->EMF+x->OutF + par[11]*x->HadF;
    else if (e<=65.0) result =  x->EMF+x->OutF + par[12]*x->HadF;
    else if (e<=70.0) result =  x->EMF+x->OutF + par[13]*x->HadF;
    else if (e<=75)   result =  x->EMF+x->OutF + par[14]*x->HadF;
    else if (e<=80.0) result =  x->EMF+x->OutF + par[15]*x->HadF;
    else if (e<=85.0) result =  x->EMF+x->OutF + par[16]*x->HadF;
    else if (e<=90.0) result =  x->EMF+x->OutF + par[17]*x->HadF;
    else if (e<=95)   result =  x->EMF+x->OutF + par[18]*x->HadF;
    else if (e<=100.0) result =  x->EMF+x->OutF + par[19]*x->HadF;
    else if (e<=105)   result =  x->EMF+x->OutF + par[20]*x->HadF;
    else if (e<=110.0) result =  x->EMF+x->OutF + par[21]*x->HadF;
    else if (e<=115.0) result =  x->EMF+x->OutF + par[22]*x->HadF;
    else if (e<=120)   result =  x->EMF+x->OutF + par[23]*x->HadF;
    else if (e<=125.0) result =  x->EMF+x->OutF + par[24]*x->HadF;
    else if (e<=130.0) result =  x->EMF+x->OutF + par[25]*x->HadF;
    else if (e<=135)   result =  x->EMF+x->OutF + par[26]*x->HadF;
    else if (e<=140.0) result =  x->EMF+x->OutF + par[27]*x->HadF;
    else if (e<=145.0) result =  x->EMF+x->OutF + par[28]*x->HadF;
    else if (e<=150.0) result =  x->EMF+x->OutF + par[29]*x->HadF;
    else if (e<=155.0) result =  x->EMF+x->OutF + par[30]*x->HadF;
    else if (e<=160)   result =  x->EMF+x->OutF + par[31]*x->HadF;
    else if (e<=165.0) result =  x->EMF+x->OutF + par[32]*x->HadF;
    else if (e<=170.0) result =  x->EMF+x->OutF + par[33]*x->HadF;
    else if (e<=175)   result =  x->EMF+x->OutF + par[34]*x->HadF;
    else if (e<=180.0) result =  x->EMF+x->OutF + par[35]*x->HadF;
    else if (e<=185.0) result =  x->EMF+x->OutF + par[36]*x->HadF;
    else if (e<=190.0) result =  x->EMF+x->OutF + par[37]*x->HadF;
    else if (e<=195)   result =  x->EMF+x->OutF + par[38]*x->HadF;
    else if (e<=200)   result =  x->EMF+x->OutF + par[39]*x->HadF;
    else if (e<=205)   result =  x->EMF+x->OutF + par[40]*x->HadF;
    else if (e<=210.0) result =  x->EMF+x->OutF + par[41]*x->HadF;
    else if (e<=215.0) result =  x->EMF+x->OutF + par[42]*x->HadF;
    else if (e<=220)   result =  x->EMF+x->OutF + par[43]*x->HadF;
    else if (e<=225.0) result =  x->EMF+x->OutF + par[44]*x->HadF;
    else if (e<=230.0) result =  x->EMF+x->OutF + par[45]*x->HadF;
    else if (e<=235)   result =  x->EMF+x->OutF + par[46]*x->HadF;
    else if (e<=240.0) result =  x->EMF+x->OutF + par[47]*x->HadF;
    else if (e<=245.0) result =  x->EMF+x->OutF + par[48]*x->HadF;
    else if (e<=250.0) result =  x->EMF+x->OutF + par[49]*x->HadF;
    else if (e<=255.0) result =  x->EMF+x->OutF + par[50]*x->HadF;
    else if (e<=260)   result =  x->EMF+x->OutF + par[51]*x->HadF;
    else if (e<=265.0) result =  x->EMF+x->OutF + par[52]*x->HadF;
    else if (e<=270.0) result =  x->EMF+x->OutF + par[53]*x->HadF;
    else if (e<=275)   result =  x->EMF+x->OutF + par[54]*x->HadF;
    else if (e<=280.0) result =  x->EMF+x->OutF + par[55]*x->HadF;
    else if (e<=285.0) result =  x->EMF+x->OutF + par[56]*x->HadF;
    else if (e<=290.0) result =  x->EMF+x->OutF + par[57]*x->HadF;
    else if (e<=295)   result =  x->EMF+x->OutF + par[58]*x->HadF;
    else if (e<=300)   result =  x->EMF+x->OutF + par[59]*x->HadF;
    else if (e<=305)   result =  x->EMF+x->OutF + par[60]*x->HadF;
    else if (e<=310.0) result =  x->EMF+x->OutF + par[61]*x->HadF;
    else if (e<=315.0) result =  x->EMF+x->OutF + par[62]*x->HadF;
    else if (e<=320)   result =  x->EMF+x->OutF + par[63]*x->HadF;
    else if (e<=325.0) result =  x->EMF+x->OutF + par[64]*x->HadF;
    else if (e<=330.0) result =  x->EMF+x->OutF + par[65]*x->HadF;
    else if (e<=335)   result =  x->EMF+x->OutF + par[66]*x->HadF;
    else if (e<=340.0) result =  x->EMF+x->OutF + par[67]*x->HadF;
    else if (e<=345.0) result =  x->EMF+x->OutF + par[68]*x->HadF;
    else if (e<=350.0) result =  x->EMF+x->OutF + par[69]*x->HadF;
    else if (e<=355.0) result =  x->EMF+x->OutF + par[70]*x->HadF;
    else if (e<=360)   result =  x->EMF+x->OutF + par[71]*x->HadF;
    else if (e<=365.0) result =  x->EMF+x->OutF + par[72]*x->HadF;
    else if (e<=370.0) result =  x->EMF+x->OutF + par[73]*x->HadF;
    else if (e<=375)   result =  x->EMF+x->OutF + par[74]*x->HadF;
    else if (e<=380.0) result =  x->EMF+x->OutF + par[75]*x->HadF;
    else if (e<=385.0) result =  x->EMF+x->OutF + par[76]*x->HadF;
    else if (e<=390.0) result =  x->EMF+x->OutF + par[77]*x->HadF;
    else if (e<=395)   result =  x->EMF+x->OutF + par[78]*x->HadF;
    else               result =  x->EMF+x->OutF + par[79]*x->HadF;
    return result;
  }
    
  double correctedJetEt(TMeasurement *const x,double *const par) const {
    return x->pt;
  }
};

// Parametrization of Jet hadronic response by a step function
class StepJetParametrization : public Parametrization { 
public:
  StepJetParametrization() : Parametrization(0,63) {}
 
  const char* name() const { return "StepJetParametrization";}

  double correctedTowerEt(TMeasurement *const x,double *const par) const {
    return x->pt;
  }
    
  double correctedJetEt(TMeasurement *const x,double *const par) const {
    double pt = x->pt;
    double Efrac = x->EMF / ( x->EMF + x->HadF + x->OutF);
    double result = 0;
    if(Efrac < 0.2)
      {
      if      (pt>=0.0   && pt<=10.0)   result = x->EMF+x->OutF + par[0]*x->HadF;
      else if (pt>10.0   && pt<=20.0)    result = x->EMF+x->OutF + par[1]*x->HadF;
      else if (pt>20.0   && pt<=30.0)    result = x->EMF+x->OutF + par[2]*x->HadF;
      else if (pt>30.0   && pt<=40.0)   result = x->EMF+x->OutF + par[3]*x->HadF;
      else if (pt>40.0  && pt<=60.0)   result = x->EMF+x->OutF + par[4]*x->HadF;
      else if (pt>60.0  && pt<=80.0)   result = x->EMF+x->OutF + par[5]*x->HadF;
      else if (pt>80.0  && pt<=100.0)   result = x->EMF+x->OutF + par[6]*x->HadF;
      else if (pt>100.0  && pt<=120.0)  result = x->EMF+x->OutF + par[7]*x->HadF;
      else if (pt>120.0 && pt<=140.0)  result = x->EMF+x->OutF + par[8]*x->HadF;
      else if (pt>140.0 && pt<=160.0)  result = x->EMF+x->OutF + par[9]*x->HadF;
      else if (pt>160.0 && pt<=180.0) result = x->EMF+x->OutF + par[10]*x->HadF;
      else if (pt>180.0 && pt<=200.0) result = x->EMF+x->OutF + par[11]*x->HadF;
      else if (pt>200.0 && pt<=225.0)  result = x->EMF+x->OutF + par[12]*x->HadF;
      else if (pt>225.0 && pt<=250.0)  result = x->EMF+x->OutF + par[13]*x->HadF;
      else if (pt>250.0 && pt<=275.0)  result = x->EMF+x->OutF + par[14]*x->HadF;
      else if (pt>275.0 && pt<=300.0) result = x->EMF+x->OutF + par[15]*x->HadF;
      else if (pt>300.0 && pt<=350.0) result = x->EMF+x->OutF + par[16]*x->HadF;
      else if (pt>350.0 && pt<=400.0)  result = x->EMF+x->OutF + par[17]*x->HadF;
      else if (pt>400.0 && pt<=500.0) result = x->EMF+x->OutF + par[18]*x->HadF;
      else if (pt>500.0 && pt<=700.0) result = x->EMF+x->OutF + par[19]*x->HadF;
      else if (pt>700.0 )               result = x->EMF+x->OutF + par[20]*x->HadF;
    } else if (Efrac < 0.5) {
      if      (pt>=0.0   && pt<=10.0)   result = x->EMF+x->OutF + par[21]*x->HadF;
      else if (pt>10.0   && pt<=20.0)    result = x->EMF+x->OutF + par[22]*x->HadF;
      else if (pt>20.0   && pt<=30.0)    result = x->EMF+x->OutF + par[23]*x->HadF;
      else if (pt>30.0   && pt<=40.0)   result = x->EMF+x->OutF + par[24]*x->HadF;
      else if (pt>40.0  && pt<=60.0)   result = x->EMF+x->OutF + par[25]*x->HadF;
      else if (pt>60.0  && pt<=80.0)   result = x->EMF+x->OutF + par[26]*x->HadF;
      else if (pt>80.0  && pt<=100.0)   result = x->EMF+x->OutF + par[27]*x->HadF;
      else if (pt>100.0  && pt<=120.0)  result = x->EMF+x->OutF + par[28]*x->HadF;
      else if (pt>120.0 && pt<=140.0)  result = x->EMF+x->OutF + par[29]*x->HadF;
      else if (pt>140.0 && pt<=160.0)  result = x->EMF+x->OutF + par[30]*x->HadF;
      else if (pt>160.0 && pt<=180.0) result = x->EMF+x->OutF + par[31]*x->HadF;
      else if (pt>180.0 && pt<=200.0) result = x->EMF+x->OutF + par[32]*x->HadF;
      else if (pt>200.0 && pt<=225.0)  result = x->EMF+x->OutF + par[33]*x->HadF;
      else if (pt>225.0 && pt<=250.0)  result = x->EMF+x->OutF + par[34]*x->HadF;
      else if (pt>250.0 && pt<=275.0)  result = x->EMF+x->OutF + par[35]*x->HadF;
      else if (pt>275.0 && pt<=300.0) result = x->EMF+x->OutF + par[36]*x->HadF;
      else if (pt>300.0 && pt<=350.0) result = x->EMF+x->OutF + par[37]*x->HadF;
      else if (pt>350.0 && pt<=400.0)  result = x->EMF+x->OutF + par[38]*x->HadF;
      else if (pt>400.0 && pt<=500.0) result = x->EMF+x->OutF + par[39]*x->HadF;
      else if (pt>500.0 && pt<=700.0) result = x->EMF+x->OutF + par[40]*x->HadF;
      else if (pt>700.0 )               result = x->EMF+x->OutF + par[41]*x->HadF;
    } else {
      if      (pt>=0.0   && pt<=10.0)   result = x->EMF+x->OutF + par[42]*x->HadF;
      else if (pt>10.0   && pt<=20.0)    result = x->EMF+x->OutF + par[43]*x->HadF;
      else if (pt>20.0   && pt<=30.0)    result = x->EMF+x->OutF + par[44]*x->HadF;
      else if (pt>30.0   && pt<=40.0)   result = x->EMF+x->OutF + par[45]*x->HadF;
      else if (pt>40.0  && pt<=60.0)   result = x->EMF+x->OutF + par[46]*x->HadF;
      else if (pt>60.0  && pt<=80.0)   result = x->EMF+x->OutF + par[47]*x->HadF;
      else if (pt>80.0  && pt<=100.0)   result = x->EMF+x->OutF + par[48]*x->HadF;
      else if (pt>100.0  && pt<=120.0)  result = x->EMF+x->OutF + par[49]*x->HadF;
      else if (pt>120.0 && pt<=140.0)  result = x->EMF+x->OutF + par[50]*x->HadF;
      else if (pt>140.0 && pt<=160.0)  result = x->EMF+x->OutF + par[51]*x->HadF;
      else if (pt>160.0 && pt<=180.0) result = x->EMF+x->OutF + par[52]*x->HadF;
      else if (pt>180.0 && pt<=200.0) result = x->EMF+x->OutF + par[53]*x->HadF;
      else if (pt>200.0 && pt<=225.0)  result = x->EMF+x->OutF + par[54]*x->HadF;
      else if (pt>225.0 && pt<=250.0)  result = x->EMF+x->OutF + par[55]*x->HadF;
      else if (pt>250.0 && pt<=275.0)  result = x->EMF+x->OutF + par[56]*x->HadF;
      else if (pt>275.0 && pt<=300.0) result = x->EMF+x->OutF + par[57]*x->HadF;
      else if (pt>300.0 && pt<=350.0) result = x->EMF+x->OutF + par[58]*x->HadF;
      else if (pt>350.0 && pt<=400.0)  result = x->EMF+x->OutF + par[59]*x->HadF;
      else if (pt>400.0 && pt<=500.0) result = x->EMF+x->OutF + par[60]*x->HadF;
      else if (pt>500.0 && pt<=700.0) result = x->EMF+x->OutF + par[61]*x->HadF;
      else if (pt>700.0 )               result = x->EMF+x->OutF + par[62]*x->HadF;
    }

    //result *= ( 1. + 0.295 * par[12] * exp(- 0.02566 * par[13] * result));   //Out of Cone, Dominant, parametrized in Et since cone R lorenz invariant
    return  result;
  }
};
#endif
