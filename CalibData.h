//
// Original Author:  Christian Autermann
//         Created:  Wed Jul 18 13:54:50 CEST 2007
// $Id: CalibData.h,v 1.34 2008/08/01 14:56:51 mschrode Exp $
//
#ifndef CalibData_h
#define CalibData_h

//#include <iostream>//cout
#include <vector> 
#include <utility> //pair
#include <cmath>

//#define __DimensionMeasurement 8 //size of _mess[] array
// _mess[0] = pT
// _mess[1] = EMF * pT
// _mess[2] = HAD * pT
// _mess[3] = OutF * pT
// _mess[4] = ?
// _mess[5] = ?
// _mess[6] = ?
// _mess[7] = ?

enum DataType {Default, TrackTower, GammaJet, TrackCluster, MessMess, PtBalance,
               InvMass, typeTowerConstraint, ParLimit};
//#define TypeDefault      0
//#define TypeTrackTower   1
//#define TypeGammaJet     2
//#define TypeTrackCluster 3
//#define TypeMessMess     4
//#define TypePtBalance    5
//#define TypeInvMass      6
//#define TypeTowerConstraint 7
//#define TypeParLimit        8


//Note: The "measurement" of a tower or a jet will be stored 
//      in a multi-dimensional class below. However, the para-
//      metrized measurement, wich will be compared to the
//      'truth' will remain a single 'double' value (i.e. the 
//       same type as 'truth')!
class TMeasurement
{
public:
  TMeasurement():pt(0.),EMF(0.),HadF(0.),OutF(0.),E(0.),eta(0.),phi(0.){};
  TMeasurement(TMeasurement* m):pt(m->pt),EMF(m->EMF),HadF(m->HadF),OutF(m->OutF),
                                E(m->E),eta(m->eta),phi(m->phi){};
  //TMeasurement(double *x){pt=x[0];EMF=x[1];HadF=x[2];OutF=x[3];E=x[6];eta=x[4];phi=x[5];};
  //all common variables
  double pt;
  double EMF;
  double HadF;
  double OutF;
  double E;
  double eta;//necessary???
  double phi;//necessary???
};

class TTower : public TMeasurement
{ 
public:
  TTower(){};
  TTower(TMeasurement* t):TMeasurement(t){};
  TTower(TTower* t):TMeasurement(t){/*further initialization*/};
//variables specific only to towers
};

class TJet : public TMeasurement
{
public:
  TJet(){};
  TJet(TMeasurement* j):TMeasurement(j){};
  TJet(TJet* j):TMeasurement(j){/*further initialization*/};
//variables specific only to jets (i.e. mass)
};

//virtual data base class -> not directly used!
class TData
{
public:
  TData(){_mess=0; _par=0; };
  TData(unsigned short int index, TMeasurement * mess, double truth, double error, double weight, double * par, unsigned short int n_par,
        double(*func)(TMeasurement*,double*),double(*err)(double*,TMeasurement*,double))
  : _index(index),_mess(mess),_truth(truth),_error(error),_weight(weight),_par(par),_n_par(n_par),_func(func),_err(err){};
  virtual ~TData(){
    delete _mess;
  };
  TMeasurement *GetMess(){ return _mess;};
  virtual double GetParametrizedMess(){return _func(_mess,_par);}
  virtual double GetParametrizedErr(double *paramess){ return _err(paramess,_mess,_error);};
  virtual double GetParametrizedErr2(double *paramess){ 
    double error = GetParametrizedErr(paramess);
    return error *error;
  };
  double GetTruth(){ return _truth;};
  virtual double GetScale(){return GetTruth();};//flatten spectrum w.r.t. this
  double GetError(){ return _error;};
  double GetWeight(){ return _weight;};
  void   SetWeight(double weight){ _weight=weight;};
  DataType GetType() {return _type;};
  void SetType(DataType type) {_type=type;};
  unsigned short int GetIndex(){return _index;};
  virtual const std::vector<TData*>& GetRef() = 0;
  virtual double chi2() = 0;
  virtual double chi2_fast(double * temp_derivative1, double * temp_derivative2, double epsilon) = 0;
  double * GetPar(){return _par;};
  unsigned short int GetNumberOfPars() {return _n_par;};
  virtual void ChangeParAddress(double* oldpar, double* newpar) { _par += newpar - oldpar;}

  static unsigned int total_n_pars;

  static double (*ScaleResidual)(double z2);         // Set to one of the following functions to scale the normalized residual z2 = chi2/weight in chi2() or chi2_fast(): 
  static double ScaleNone(double z2) { return z2; }  // No scaling of residuals in chi2() or chi2_fast() (default)
  static double ScaleCauchy(double z2);	             // Scaling of residuals with Cauchy-Function in chi2() or chi2_fast()
  static double ScaleHuber(double z2);               // Scaling of residuals with Huber-Function in chi2() or chi2_fast()

protected:
  unsigned short int _index; //limited from 0 to 65535
  TMeasurement *_mess;
  double _truth, _error, _weight;
  double *_par;
  unsigned short int _n_par; //limited from 0 to 65535
  double (*_func)(TMeasurement * x, double * par);
  double (*_err)(double * x, TMeasurement * x_original, double error);
  DataType _type;
};

//data class for data providing one truth and one messurement, 
//e.g. track-tower
class TData_TruthMess : public TData
{
public:
  TData_TruthMess(unsigned short int index, TMeasurement * mess, double truth, double error, double weight, double * par, unsigned short int n_par,
        double(*func)(TMeasurement*,double*),double(*err)(double*,TMeasurement*,double))
  : TData(index, mess, truth, error, weight, par, n_par, func, err){_type=TrackTower;};

  virtual const std::vector<TData*>& GetRef() { 
    resultcache.clear();	
    resultcache.push_back( this );
    return resultcache;
  };

  virtual double chi2(){ 
    double new_mess  = GetParametrizedMess();
    double new_error = GetParametrizedErr(&new_mess);
    return GetWeight() * (*TData::ScaleResidual)( (_truth-new_mess)*(_truth-new_mess)/(new_error*new_error) );
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
		      double(*func)(TMeasurement*,double*), double(*err)(double*,TMeasurement*,double), 
		      TMeasurement *mess)
  : TData_TruthMess(index, mess, truth, error, weight, par, n_par, func, err){_type=GammaJet;};
  virtual ~TData_TruthMultMess() {
    for (std::vector<TData*>::const_iterator it=_vecmess.begin();
	 it!=_vecmess.end(); ++it)
      delete *it;
    _vecmess.clear();	
  };
  void AddMess(TData_TruthMess * m){ _vecmess.push_back(m);};

  virtual double GetParametrizedMess(){
    double tower_pt_sum=0.0;
    for (std::vector<TData*>::const_iterator it=_vecmess.begin();
  	 it!=_vecmess.end(); ++it){
      tower_pt_sum += (*it)->GetParametrizedMess(); // Sum of tower Pt
    }
    TJet jet(_mess);
    jet.pt = tower_pt_sum;
    return _func(&jet, _par);
  };
  
  virtual double chi2(){ 
    double weight = GetWeight(); 
    double new_mess  = GetParametrizedMess();             
    double new_error =  GetParametrizedErr( &new_mess );
    return (new_error!=0. ? weight*(*TData::ScaleResidual)(_truth-new_mess)*(_truth-new_mess)/(new_error*new_error):0.0);
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
//Base class for PT-balance and Inv-Mass !!!
class TData_MessMess : public TData_TruthMultMess
{
public:
    TData_MessMess(unsigned short int index, double * dir, double truth, double error, 
                   double weight, double * par, unsigned short int n_par,
		   double(*func)(TMeasurement*,double*),double(*err)(double*,TMeasurement*,double), 
		   TMeasurement *mess)
    : TData_TruthMultMess(index, truth, error, weight, par, n_par, func, err,mess),
      _direction(dir){_type=MessMess;};
    virtual ~TData_MessMess() {
      for (std::vector<TData_MessMess*>::const_iterator it=_m2.begin();
	   it!=_m2.end(); ++it)
	delete *it;
      _m2.clear();
      delete [] _direction;	
    };
    virtual void AddNewMultMess(TData_MessMess * m2 ){assert(m2->_m2.empty());_m2.push_back(m2);};
    void ClearMultMess() { _m2.clear();}
    virtual double chi2(){ 
      double sum_error2=0.0, new_error, new_mess;
      double weight = GetWeight();
      for (std::vector<TData_MessMess*>::const_iterator it=_m2.begin();
           it!=_m2.end(); ++it) {
	 new_mess    = (*it)->GetParametrizedMess();	 
  	 new_error   = (*it)->GetParametrizedErr(&new_mess);
	 sum_error2 += new_error * new_error;
      }
      new_mess  = GetParametrizedMess();
      new_error = GetParametrizedErr( &new_mess );
      new_mess  = GetMessCombination();
      sum_error2 = new_error*new_error;
      return (sum_error2!=0 ? weight*(_truth-new_mess)*(_truth-new_mess)/sum_error2 : 0.0);
    };
    virtual double chi2_fast(double * temp_derivative1, double*  temp_derivative2, double epsilon);
    virtual void ChangeParAddress(double* oldpar, double* newpar) { 
      TData::ChangeParAddress(oldpar,newpar);
      for (std::vector<TData*>::iterator it=_vecmess.begin();  it !=_vecmess.end(); ++it) 
	(*it)->ChangeParAddress(oldpar,newpar);
      for (std::vector<TData_MessMess*>::const_iterator it=_m2.begin(); it!=_m2.end(); ++it)
	(*it)->ChangeParAddress(oldpar,newpar);
    }
    virtual double GetMessCombination(){ return combine(); }; // for plotting
    virtual double * GetDirection(){ return _direction; };
    virtual std::vector<TData_MessMess*>* GetSecondaryJets(){return &_m2;};
    //virtual TMeasurement* GetMultMess(int i) { 
    //  if(i == 0) return GetMess();
    //  return _m2[i-1]->GetMess();
    //}
    //virtual double GetMultParametrizedMess(int i) { // for plotting
    //  if(i == 0) return GetParametrizedMess();
    //  return _m2[i-1]->GetParametrizedMess();
    //}
   
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
		  double(*func)(TMeasurement*,double*),double(*err)(double*,TMeasurement*,double), 
		  TMeasurement *mess)
  : TData_MessMess(index, dir, truth*truth, error, weight, par, n_par, func, err, mess){_type=PtBalance;};
  virtual ~TData_PtBalance(){};
  virtual double GetTruth(){ return sqrt(_truth);};
  virtual double GetScale(){
    double scale = GetMess()->pt;
    for (std::vector<TData_MessMess*>::const_iterator it=_m2.begin();it!=_m2.end(); ++it)
      scale += (*it)->GetMess()->pt;
    //@@ return scale/double(1+_m2.size());
    return scale/2.;
  }
protected:
  virtual double combine();  
};

//e.g. top-Mass, W-mass, etc..
class TData_InvMass2 : public TData_MessMess
{
public:
  TData_InvMass2(unsigned short int index, double * dir, double truth, double error, double weight,
		double * par, unsigned short int n_par,
		double(*func)(TMeasurement*,double*),	double(*err)(double*,TMeasurement*,double), 
		TMeasurement *mess)
  : TData_MessMess(index, dir, truth*truth, error, weight, par, n_par, func, err, mess) { _type=InvMass; };
  virtual ~TData_InvMass2(){};
  virtual double GetTruth(){ return sqrt(_truth);};
 protected:
  virtual double combine();  
};

//data class to limit a parameter
class TData_ParLimit : public TData
{
 public:
  TData_ParLimit(unsigned short int index, TMeasurement *mess, double error, double *par,double(*func)(TMeasurement*,double*)) 
  : TData(index,mess,0,error,1.0,par,1,func,0){ _type=ParLimit;
    };
    
    virtual const std::vector<TData*>& GetRef() { 
      _cache.clear();
      _cache.push_back(this);
      return _cache;
    };
    
    virtual double GetParametrizedErr(double *paramess){ 
      return _error;
    };
        
    virtual double chi2(){  
      double new_mess  = GetParametrizedMess();
      double new_error = GetParametrizedErr(&new_mess);
      return new_mess * new_mess / (new_error * new_error);
    };
    virtual double chi2_fast(double* temp_derivative1, 
			     double* temp_derivative2, double epsilon);
    
 private:
    static std::vector<TData*> _cache;
};

#endif
