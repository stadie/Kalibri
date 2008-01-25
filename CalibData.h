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
#define TypePtBalance    5

//virtual data base class -> not directly used!
class TData
{
public:
  TData(){_mess=0; _par=0;};
  TData(unsigned short int index, double * mess, double truth, double error, double * par, unsigned short int n_par,
        double(*func)(double*,double*), double(*err)(double*)){
    _index=index;_mess=mess;_truth=truth;_error=error;_par=par;_n_par=n_par;_func=func;_err=err;};
  virtual ~TData(){
    //if (_mess) delete [] _mess;
  };
  virtual double * GetMess(){  return _mess;};//used only for plotting
  virtual double GetParametrizedMess(){ 
	       return _func(_mess,_par);}
  virtual double GetParametrizedErr(double *paramess){ return _err(paramess);};
  virtual double GetTruth(){ return _truth;};
  virtual double GetError(){ return _error;};
  virtual short unsigned int GetType(){return _type;};
  virtual void SetType(short unsigned int type){_type=type;};
  unsigned short int GetIndex(){return _index;};
  virtual const std::vector<TData*>& GetRef() = 0;
  virtual double chi2() = 0;
  virtual double chi2_fast() = 0;
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

  virtual const std::vector<TData*>& GetRef() { 
    //std::vector<TData*> result;
    resultcache.clear();	
    resultcache.push_back( this );
    return resultcache;
  };
  virtual double chi2(){ 
    double new_mess  = GetParametrizedMess();
    double new_error = GetParametrizedErr(&new_mess);
    return (_truth-new_mess)*(_truth-new_mess)/(new_error*new_error);
  };
  virtual double chi2_fast();
  
  private:
  	static std::vector<TData*> resultcache;
};


//data class for data providing one truth and multiple messurements, 
//e.g. gamma-jet or track-cluster
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
      _vecmess.clear();	
  };
  void   AddMess(TData_TruthMess * m){_vecmess.push_back(m);};
  virtual double * GetMess(){//used only for plotting
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
  virtual double GetParametrizedMess(){
    double result=0.0;
    for (std::vector<TData*>::const_iterator it=_vecmess.begin();
  	 it!=_vecmess.end(); ++it){
      result+=(*it)->GetParametrizedMess();
    }
    return _func(&result, _par);	  
  };

  virtual double chi2(){ 
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
  virtual double chi2_fast();
  virtual const std::vector<TData*>& GetRef() {return _vecmess;};

protected:  
  std::vector<TData*> _vecmess;
};

//data class for data providing only multiple messurements, no truth 
//e.g. jet-jet, jet-jet-jet, toptop, etc..
class TData_MessMess : public TData_TruthMultMess
{
public:
    TData_MessMess(unsigned short int index, double * dir, double truth, double error, 
                   double * par, unsigned short int n_par,
        double(*func)(double*,double*),	double(*err)(double*)):
    TData_TruthMultMess(index, truth, error, par, n_par, func, err){
      _direction=dir; _type=TypeJetJet;};
    virtual ~TData_MessMess() {
      for (std::vector<TData*>::const_iterator it=_vecmess.begin();
	   it!=_vecmess.end(); ++it)
	delete *it;
      _vecmess.clear();	
      for (std::vector<TData_MessMess*>::const_iterator it=_m2.begin();
	   it!=_m2.end(); ++it)
	(*it)->~TData_MessMess();
      _m2.clear();
      delete [] _direction;	
    };
    virtual void AddNewMultMess(TData_MessMess * m2 ){_m2.push_back(m2);};
    virtual double GetMessCombination(){ return combine(); };
    virtual double * GetDirection(){ return _direction; };
    virtual double chi2(){ 
      double sum_error2=0.0, new_error, new_mess;
      for (std::vector<TData_MessMess*>::const_iterator it=_m2.begin();
           it!=_m2.end(); ++it) {
	 new_mess    = (*it)->GetParametrizedMess();	 
  	 new_error   = (*it)->GetParametrizedErr(&new_mess);
	 sum_error2 += new_error * new_error;
      }
      //new_mess  = _func( _mess, _par);  
      new_mess  = GetParametrizedMess();
      new_error = _err(  &new_mess );
      new_mess  = GetMessCombination();  
      return (_truth-new_mess)*(_truth-new_mess)/(sum_error2 + new_error*new_error);
    };

protected:
  virtual double combine() = 0;  
  std::vector<TData_MessMess*> _m2;
  double * _direction;
};

//e.g. jet-jet, jet-jet-jet, etc..
class TData_PtBalance : public TData_MessMess
{
public:
    TData_PtBalance(unsigned short int index, double * dir, double truth, double error, double * par, unsigned short int n_par,
        double(*func)(double*,double*),	double(*err)(double*)):
    TData_MessMess(index, dir, truth, error, par, n_par, func, err){_type=TypePtBalance;};
    virtual ~TData_PtBalance(){};
    virtual double chi2_fast();
protected:
    virtual double combine();  
};

#endif
