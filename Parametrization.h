//
// Original Author:  Hartmut Stadie
//         Created:  Thu Apr 03 17:09:50 CEST 2008
// $Id: Parametrization.h,v 1.4 2008/06/30 13:24:09 stadie Exp $
//
#ifndef CALIBCORE_PARAMETRIZATION_H
#define CALIBCORE_PARAMETRIZATION_H

#include <cmath>

class Parametrization 
{
public:
  Parametrization(unsigned int ntowerpars, unsigned int njetpars) 
    : ntowerpars_(ntowerpars), njetpars_(njetpars) {}
  virtual ~Parametrization() {}

  /** correctedTowerEt(double *x,double *par)
      returns the corrected Et of a tower 
           input: x[4]: x[0] : Et of whole tower
                             x[1] : Et of ECAL part
                             x[2] : Et of HCAL part
                             x[3] : Et of Outer part
           par:  the correction parameters of this tower
  **/
  virtual double correctedTowerEt(double *x,double *par) const = 0;

  /** correctedJetEt(double *x,double *par)
      returns the corrected Et of a jet 
      input: x[3]: x[0] : Et of uncorrected jet
                        x[1] : eta of uncorrected jet(not used)
                        x[2] : phi of uncorrected jet(not used)
           par:  the correction parameters of this jet
  **/
  virtual double correctedJetEt(double *x,double *par) const = 0;
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
  StepParametrization() : Parametrization(12,10) {}
 
  const char* name() const { return "StepParametrization";}

  double correctedTowerEt(double *x,double *par) const {
    double result = 0;
    
    if(x[2]>=0.0  && x[2]<=1.0)  result = x[1]+x[3] + par[0]*x[2];
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
    return result;
  }
    
  double correctedJetEt(double *x,double *par) const {
    //return  par[0]*x[0] + par[1];
    double result = 0;
    
    if(x[0]>=0.0        && x[0]<=5.0)    result = par[0]*x[0];
    else if (x[0]>5.0   && x[0]<=10.0)   result = par[1]*x[0];
    else if (x[0]>10.0  && x[0]<=20.0)   result = par[2]*x[0];
    else if (x[0]>20.0  && x[0]<=40.0)   result = par[3]*x[0];
    else if (x[0]>40.0  && x[0]<=80.0)   result = par[4]*x[0];
    else if (x[0]>80.0  && x[0]<=160.0)  result = par[5]*x[0];
    else if (x[0]>160.0 && x[0]<=300.0)  result = par[6]*x[0];
    else if (x[0]>300.0 && x[0]<=600.0)  result = par[7]*x[0];
    else if (x[0]>600.0 && x[0]<=1000.0) result = par[8]*x[0];
    else if (x[0]>1000.0 )               result = par[9]*x[0];
    return result;
  }
};

// Parametrization of hadronic response by a step function
// 3 Sets of Parametrization for different EM fraction

class StepEfracParametrization : public Parametrization {
public:
  StepEfracParametrization() : Parametrization(36,2) {}
  
  const char* name() const { return "StepEfracParametrization";}

  double correctedTowerEt(double *x,double *par) const {
    double result=0;
    
    //double Efrac = x[1]/(x[2]+x[3]);
    if( x[1] < 0.1 * (x[2]+x[3]) ) {
      if      (x[2]>=0.0   && x[2]<=1.0)   result = x[1]+x[3] + par[0]*x[2];
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
    } else if (x[1]<0.3*(x[2]+x[3])) {
      if      (x[2]>=0.0   && x[2]<=1.0)   result = x[1]+x[3] + par[12]*x[2];
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
      if      (x[2]>=0.0   && x[2]<=1.0)   result = x[1]+x[3] + par[24]*x[2];
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
    return result;
  }
  
  double correctedJetEt(double *x,double *par) const {
    return  par[0]*x[0] + par[1];
  }
};

// Parametrization of response by some "clever" function
class MyParametrization: public Parametrization {
public:
  MyParametrization() : Parametrization(3,2) {}
  const char* name() const { return "MyParametrization";}
  double correctedTowerEt(double *x,double *par) const {
    return x[1] + par[0]*x[2] + par[1]*log(x[0]) + par[2];
  }
  double correctedJetEt(double *x,double *par) const {
    return par[0]*x[0] + par[1];
  }
};

// Parametrization of response with some ideas from the JetMET group
class JetMETParametrization: public Parametrization {
public:
  JetMETParametrization() : Parametrization(3,5) {}
  const char* name() const { return "JetMETParametrization";}
  double correctedTowerEt(double *x,double *par) const {
    if(par[0] < -10) par[0] = -10;
    if(par[1] < 0) par[1] = -par[1];
    if(par[2] < 0) par[2] = -par[2];
    return par[1] * x[2] + par[2] * x[1] + x[3] + par[0];
  }
  double correctedJetEt(double *x,double *par) const {
    double logx = log(x[0]);
    if(logx < 0) logx = 0;
    if(par[1] < 0) par[1] *= -1;
    if(par[2] < 0) par[2] *= -1;
    if(par[3] < 0) par[3] *= -1;
    if(par[4] < 0) par[4] *= -1;
    return (par[0] - par[1]/(pow(logx,par[2]) + par[3]) + par[4]/x[0]) * x[0];  
  }
};

// Parametrization for toy MC
class ToyParametrization: public Parametrization {
public:
  ToyParametrization() : Parametrization(3,0) {}
  const char* name() const { return "ToyParametrization";}
  double correctedTowerEt(double *x,double *par) const {
    return par[0] * x[2] + par[1] * x[1] + x[3] + par[2];
  }
  double correctedJetEt(double *x,double *par) const {
    return x[0];  
  }
};
#endif
