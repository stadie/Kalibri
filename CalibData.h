//s
// Original Author:  Christian Autermann
//         Created:  Wed Jul 18 13:54:50 CEST 2007
// $Id: CalibData.h,v 1.20 2008/06/27 14:24:32 mschrode Exp $
//
#ifndef CalibData_h
#define CalibData_h

#include <TMath.h>

#include <iostream>
#include <vector> 
#include <map> 
#include <utility> //pair
#include <cmath>

#define __DimensionMeasurement 7

#define TypeDefault      0
#define TypeTrackTower   1
#define TypeGammaJet     2
#define TypeTrackCluster 3
#define TypeMessMess     4
#define TypePtBalance    5
#define TypeInvMass      6
#define TypeTowerConstraint 7
//#define __FastErrorCalculation

//virtual data base class -> not directly used!
class TData
{
public:
  TData(){_mess=0; _par=0; };
  TData(unsigned short int index, double * mess, double truth, double error, double weight, double * par, unsigned short int n_par,
        double(*func)(double*,double*), double(*err)(double*)){
    _index=index;_mess=mess;_truth=truth;_error=error;_par=par;_n_par=n_par;_func=func;_err=err;_weight=weight; };
  virtual ~TData(){
    delete [] _mess;
  };
  double *GetMess(){ return _mess;};//used only for plotting
  virtual double GetParametrizedMess(){return _func(_mess,_par);}
  virtual double GetParametrizedErr(double *paramess){ return _err(paramess);};
  double GetTruth(){ return _truth;};
  virtual double GetScale(){return GetTruth();};//flatten spectrum w.r.t. this
  double GetError(){ return _error;};
  double GetWeight(){ return _weight;};
  void   SetWeight(double weight){ _weight=weight;};
  short unsigned int GetType() {return _type;};
  void SetType(short unsigned int type) {_type=type;};
  unsigned short int GetIndex(){return _index;};
  virtual const std::vector<TData*>& GetRef() = 0;
  virtual double chi2() = 0;
  virtual double chi2_fast(double * temp_derivative1, double*  temp_derivative2, double epsilon) = 0;
  double * GetPar(){return _par;};
  unsigned short int GetNumberOfPars() {return _n_par;};
  virtual void ChangeParAddress(double* oldpar, double* newpar) { _par += newpar - oldpar;}

  static unsigned int total_n_pars;

  static double (*ScaleResidual)(double z2);         // Set to one of the following functions to scale the normalized residual z2 = chi2/weight in chi2() or chi2_fast(): 
  static double ScaleNone(double z2) { return z2; }  // No scaling of residuals in chi2() or chi2_fast() (default)
  static double ScaleCauchy(double z2);	             // Scaling of residuals with Cauchy-Function in chi2() or chi2_fast()
  static double ScaleHuber(double z2);               // Scaling of residuals with Huber-Function in chi2() or chi2_fast()

protected:
  unsigned short int _index, _n_par, _type; //limited from 0 to 65535
  double _truth, _error, _paramess, _weight;
  double *_mess, *_par;
  double (*_func)(double * x, double * par);
  double (*_err)(double * x);
};

//data class for data providing one truth and one messurement, 
//e.g. track-tower
class TData_TruthMess : public TData
{
public:
  TData_TruthMess(unsigned short int index, double * mess, double truth, double error, double weight, double * par, unsigned short int n_par,
        double(*func)(double*,double*),double(*err)(double*)):
    TData(index, mess, truth, error, weight, par, n_par, func, err){_type=TypeTrackTower;};

  virtual const std::vector<TData*>& GetRef() { 
    resultcache.clear();	
    resultcache.push_back( this );
    return resultcache;
  };
  virtual double chi2(){ 
    double new_mess  = GetParametrizedMess();
    double *mess = GetMess();
    //double new_error = GetParametrizedErr(&new_mess);
    double new_messE = new_mess * mess[6] / mess[0]; //Et->E  
    double new_error = GetParametrizedErr(&new_messE) * mess[0] / mess[6]; 
    double weight = GetWeight();
    //-- return weight*(_truth-new_mess)*(_truth-new_mess)/(new_error*new_error);
    return weight * (*TData::ScaleResidual)( (_truth-new_mess)*(_truth-new_mess)/(new_error*new_error) );
  };
  virtual double chi2_fast(double * temp_derivative1, double*  temp_derivative2, double epsilon);
  
private:
  static std::vector<TData*> resultcache;
};


//data class for data providing one truth and multiple messurements, 
//e.g. gamma-jet or track-cluster
class TData_TruthMultMess : public TData_TruthMess
{
public:
  TData_TruthMultMess(unsigned short int index, double truth, double error, 
		      double weight, double * par, unsigned short int n_par,
		      double(*func)(double*,double*), double(*err)(double*), double *mess) :
    TData_TruthMess(index, mess, truth, error, weight, par, n_par, func, err){_type=TypeGammaJet;};
  virtual ~TData_TruthMultMess() {
    for (std::vector<TData*>::const_iterator it=_vecmess.begin();
	 it!=_vecmess.end(); ++it)
      delete *it;
    _vecmess.clear();	
  };
  void AddMess(TData_TruthMess * m){ _vecmess.push_back(m);};
  // use GetMess from base class
/*   virtual double * GetMess(){//used only for plotting */
/*     double *dummy, *result=new double[__DimensionMeasurement]; */
/*     for (unsigned i=0; i<__DimensionMeasurement; ++i) result[i]=0.0; */
/*     for (std::vector<TData*>::const_iterator it=_vecmess.begin(); */
/*    	 it!=_vecmess.end(); ++it){ */
/*       dummy = (*it)->GetMess(); */
/*       for (int i=0; i<__DimensionMeasurement; ++i)	  */
/* 	result[i]+=dummy[i]; */
/*     }	  */
/*     return result;	   */
/*   }; */
  virtual double GetParametrizedMess(){
    double result=0.0;
    for (std::vector<TData*>::const_iterator it=_vecmess.begin();
  	 it!=_vecmess.end(); ++it){
      result+=(*it)->GetParametrizedMess();
    }
    return _func(&result, _par);	  
  };

  virtual double chi2(){ 
    double sum_mess=0.0, sum_error2=0.0, new_error, new_mess, new_messE;
    double weight = GetWeight();
    double *mess;
    for (std::vector<TData*>::const_iterator it=_vecmess.begin();
         it!=_vecmess.end(); ++it) {
       new_mess    = (*it)->GetParametrizedMess();
       mess = (*it)->GetMess();	 
       sum_mess   += new_mess;
       //new_error   = (*it)->GetParametrizedErr(&new_mess);
       new_messE =  new_mess  * mess[6] / mess[0]; //Et->E  
       new_error = (*it)->GetParametrizedErr(&new_messE) * mess[0] / mess[6];
       sum_error2 += new_error * new_error;
    }
    new_mess  = _func( &sum_mess, _par);  
    new_error = _err(  &new_mess );
    //-- return weight*(_truth-new_mess)*(_truth-new_mess)/(sum_error2 + new_error*new_error);
    return weight*(*TData::ScaleResidual)( (_truth-new_mess)*(_truth-new_mess)/(sum_error2 + new_error*new_error) );
  };
  virtual double chi2_fast(double * temp_derivative1, double*  temp_derivative2, double epsilon);
  virtual const std::vector<TData*>& GetRef() {return _vecmess;};
  virtual void ChangeParAddress(double* oldpar, double* newpar) { 
    TData::ChangeParAddress(oldpar,newpar);
    for (std::vector<TData*>::iterator it=_vecmess.begin();
         it !=_vecmess.end(); ++it) { (*it)->ChangeParAddress(oldpar,newpar);}
  }
protected:  
  std::vector<TData*> _vecmess;
};

//virtual data class for data providing only multiple messurements
//NOT DIRECTLY USED, SINCE "combine" FUNCTION IS NOT YET DEFINED HERE!
class TData_MessMess : public TData_TruthMultMess
{
public:
    TData_MessMess(unsigned short int index, double * dir, double truth, double error, 
                   double weight, double * par, unsigned short int n_par,
		   double(*func)(double*,double*), double(*err)(double*), double *mess):
      TData_TruthMultMess(index, truth, error, weight, par, n_par, func, err, mess){
      _direction=dir; _type=TypeMessMess;};
    virtual ~TData_MessMess() {
      for (std::vector<TData_MessMess*>::const_iterator it=_m2.begin();
	   it!=_m2.end(); ++it)
	delete *it;
      _m2.clear();
      delete [] _direction;	
    };
    virtual void AddNewMultMess(TData_MessMess * m2 ){_m2.push_back(m2);};
    virtual double GetMultParametrizedMess(int i) {
      if(i == 0) return GetParametrizedMess();
      return _m2[i-1]->GetParametrizedMess();
    }
    virtual double* GetMultMess(int i) {
      if(i == 0) return GetMess();
      return _m2[i-1]->GetMess();
    }
    virtual double GetMessCombination(){ return combine(); };
    virtual double * GetDirection(){ return _direction; };
    virtual double chi2(){ 
      double sum_error2=0.0, new_error, new_mess, new_messE; 
      double weight = GetWeight();
      for (std::vector<TData_MessMess*>::const_iterator it=_m2.begin();
           it!=_m2.end(); ++it) {
	 new_mess    = (*it)->GetParametrizedMess();	 
  	 //new_error   = (*it)->GetParametrizedErr(&new_mess);
	 double *mess = (*it)->GetMess();
	 new_messE =  new_mess  * mess[6] / mess[0]; //Et->E   
	 new_error = (*it)->GetParametrizedErr(&new_messE) * mess[0] / mess[6];
	 sum_error2 += new_error * new_error;
      }
      new_mess  = GetParametrizedMess();
      new_error = _err( &new_mess );
      new_mess  = GetMessCombination();

      //Combination -> B -> to response
      /*
      new_mess /= GetScale();
      new_mess = (2+new_mess)/(2-new_mess);
      new_error *= 4 / ( 4 - (4*new_mess) + (new_mess*new_mess));
      sum_error2 += new_error * new_error;
      */


      //std::cout << "new_mess:" << new_mess << " error:" << new_error << " and " << sqrt(sum_error2) 
      //		<< '\n';
      //-- return weight*(_truth-new_mess)*(_truth-new_mess)/(sum_error2 + new_error*new_error);
      return weight*(*TData::ScaleResidual)( (_truth-new_mess)*(_truth-new_mess)/(sum_error2 + new_error*new_error) );
    };
    virtual double chi2_fast(double * temp_derivative1, double*  temp_derivative2, double epsilon);
    virtual void ChangeParAddress(double* oldpar, double* newpar) { 
      TData::ChangeParAddress(oldpar,newpar);
      for (std::vector<TData*>::iterator it=_vecmess.begin();  it !=_vecmess.end(); ++it) 
	(*it)->ChangeParAddress(oldpar,newpar);
      for (std::vector<TData_MessMess*>::const_iterator it=_m2.begin(); it!=_m2.end(); ++it)
	(*it)->ChangeParAddress(oldpar,newpar);
    }
    virtual double GetScale(){
      double sum=0.;
      for (std::vector<TData_MessMess*>::const_iterator it= _m2.begin(); it !=_m2.end(); ++it) { 
	sum+=(*it)->GetMess()[0];
      }
      sum = (sum + GetMess()[0])/2.;	 
      return sum;
    }
protected:
  virtual double combine() = 0;
  std::vector<TData_MessMess*> _m2;
  double * _direction;
};

//e.g. jet-jet, jet-jet-jet, etc..
class TData_PtBalance : public TData_MessMess
{
public:
  TData_PtBalance(unsigned short int index, double * dir, double truth, double error, 
		  double weight, double * par, unsigned short int n_par,
		  double(*func)(double*,double*),	double(*err)(double*), double *mess) 
    : TData_MessMess(index, dir, truth*truth, error, weight, par, n_par, func, err, mess) {
    _type=TypePtBalance;
  };
  virtual ~TData_PtBalance(){};
  virtual double GetTruth(){ return sqrt(_truth);};
protected:
  virtual double combine();  
};

//e.g. top-Mass, W-mass, etc..
class TData_InvMass : public TData_MessMess
{
public:
  TData_InvMass(unsigned short int index, double * dir, double truth, double error, double weight,
		double * par, unsigned short int n_par,
		double(*func)(double*,double*),	double(*err)(double*), double *mess) 
    : TData_MessMess(index, dir, truth*truth, error, weight, par, n_par, func, err, mess) {
    _type=TypeInvMass;
  };
  virtual ~TData_InvMass(){};
  virtual double GetTruth(){ return sqrt(_truth);};
 protected:
  virtual double combine();  
};

#endif
