#ifndef CalibData_h
#define CalibData_h

#include <iostream>
#include <vector> 
#include <map> 
#include <utility> //pair
//#include <cmath>

#define TypeDefault      0
#define TypeTrackTower   1
#define TypeGammaJet     2
#define TypeTrackCluster 3
#define TypeJetJet       4

//virtual data base class -> not directly used!
class TData
{
public:
  TData(){};
  TData(unsigned short int index, double * mess, double truth, double error, double * par, unsigned short int n_par,
        double(*func)(double*,double*), double(*err)(double*)){
    _index=index;_mess=mess;_truth=truth;_error=error;_par=par;_n_par=n_par;_func=func;_err=err;};
  virtual ~TData(){};
  virtual double * GetMess(){  return _mess;};//used only for plotting
  virtual double GetParametrizedMess(){ return _func(_mess,_par);};
  virtual double GetParametrizedErr(double *_paramess){ return _err(_paramess);};
  virtual double GetTruth(){ return _truth;};
  virtual double GetError(){ return _error;};
  virtual short unsigned int GetType(){return _type;};
  virtual void SetType(short unsigned int type){_type=type;};
  unsigned short int GetIndex(){return _index;};
  virtual std::vector<TData*> GetRef(){};
  virtual double chi2(){};
  virtual double chi2_fast(){};
  virtual double * GetPar(){return _par;};
  unsigned short int GetNumberOfPars() {return _n_par;};
  void  AddToPar(unsigned short int const i, double const e){_par[i]+=e;};

  static unsigned int total_n_pars;
  double static * temp_derivative1;
  double static * temp_derivative2;
  double const static epsilon;

protected:
  unsigned short int _index, _n_par, _type; //limited from 0 to 65535
  double _truth, _error, _paramess;
  double *_mess, *_par;
  double (*_func)(double * x, double * par);
  double (*_err)(double * x);
};

//data class for data providing one truth and one messurement, 
//e.g. track-tower
class TData_TruthMess : public TData
{
public:
  TData_TruthMess(unsigned short int index, double * mess, double truth, double error, double * par, unsigned short int n_par,
        double(*func)(double*,double*),double(*err)(double*)):
    TData(index, mess, truth, error, par, n_par, func, err){_type=TypeTrackTower;};

  virtual std::vector<TData*> GetRef(){ 
    std::vector<TData*> result;
    result.push_back( this );
    return result;
  };
  double chi2(){ 
    double new_mess  = GetParametrizedMess();
    double new_error = GetParametrizedErr(&new_mess);
    return (_truth-new_mess)*(_truth-new_mess)/(new_error*new_error);
  };
  double chi2_fast();
};

//data class for data providing one truth and multiple messurements, 
//e.g. gamma-jet
class TData_TruthMultMess : public TData_TruthMess
{
public:
  TData_TruthMultMess(unsigned short int index, double truth, double error, double * par, unsigned short int n_par,
        double(*func)(double*,double*),double(*err)(double*)):
    TData_TruthMess(index, 0, truth, error, par, n_par, func, err){_type=TypeGammaJet;};
    virtual ~TData_TruthMultMess() {
      for (std::vector<TData*>::const_iterator it=_vecmess.begin();
	   it!=_vecmess.end(); ++it)
	delete *it;
    }
  void   AddMess(TData_TruthMess * m){_vecmess.push_back(m);};
  double * GetMess(){//used only for plotting
    double *dummy, *result=new double[4];
    result[0]=0.0;result[1]=0.0;result[2]=0.0;result[3]=0.0;
    for (std::vector<TData*>::const_iterator it=_vecmess.begin();
  	 it!=_vecmess.end(); ++it){
	 dummy = (*it)->GetMess();
       for (int i=0; i<4; ++i)	 
         result[i]+=dummy[i];
    }	 
    return result;	  
  };
  double GetParametrizedMess(){
    double result=0.0;
    for (std::vector<TData*>::const_iterator it=_vecmess.begin();
  	 it!=_vecmess.end(); ++it){
      result+=(*it)->GetParametrizedMess();
    }
    return _func(&result, _par);	  
  };

  double chi2(){ 
    double sum_mess=0.0, sum_error2=0.0, new_error, new_mess;
    for (std::vector<TData*>::const_iterator it=_vecmess.begin();
         it!=_vecmess.end(); ++it) {
       new_mess    = (*it)->GetParametrizedMess();	 
       sum_mess   += new_mess;
       new_error   = (*it)->GetParametrizedErr(&new_mess);
       sum_error2 += new_error * new_error;
    }
    new_mess  = _func( &sum_mess, _par);  
    new_error = _err(  &new_mess );
    return (_truth-new_mess)*(_truth-new_mess)/(sum_error2 + new_error*new_error);
  };
  double chi2_fast();
  virtual std::vector<TData*> GetRef(){return _vecmess;};

protected:  
  std::vector<TData*> _vecmess;
};

//data class for data providing only multiple messurements, no truth 
//e.g. jet-jet
class TData_MessMess : public TData
{
public:
  TData_MessMess(unsigned short int index, double * mess, double error, double * par, unsigned short int n_par,
        double(*func)(double*,double*),double(*err)(double*)):
     TData(index, mess,0,error,par,n_par,func,err){_m2=0;_type=TypeJetJet;};
  TData_MessMess(unsigned short int index, double * mess, double error, double * par, unsigned short int n_par,
        double(*func)(double*,double*), double(*err)(double*), TData_MessMess * m2):
     TData(index, mess,0,error,par,n_par,func,err){_m2=m2;_type=TypeJetJet;};

  TData_MessMess * GetMess2Ref(){return _m2;};
  void   SetMess2Ref(TData_MessMess * m2){_m2=m2;};

  virtual std::vector<TData*> GetRef(){ 
    std::vector<TData*> result;
    result.push_back(this);
    return result;
  };

protected:
  TData_MessMess * _m2;
};

#endif
