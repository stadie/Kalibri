//
// Original Author:  Christian Autermann
//         Created:  Wed Jul 18 13:54:50 CEST 2007
// $Id: CalibData.h,v 1.46 2008/10/13 15:15:05 rwolf Exp $
//
#ifndef CalibData_h
#define CalibData_h

#include <iostream>
using namespace std;
#include <vector> 
#include <cmath>

enum DataType {Default, TrackTower, GammaJet, TrackCluster, MessMess, PtBalance,
               InvMass, typeTowerConstraint, ParLimit};

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
  TTower():TMeasurement(){};
  TTower(TMeasurement* t):TMeasurement(t){};
  TTower(TTower* t):TMeasurement(t){/*further initialization*/};
//variables specific only to towers (i.e. # EM cells)
};

class TJet : public TMeasurement
{
public:
  TJet():TMeasurement(){};
  TJet(TMeasurement* j):TMeasurement(j){};
  TJet(TJet* j):TMeasurement(j){/*further initialization*/};
//variables specific only to jets (i.e. mass)
};

//virtual data base class -> not directly used!
class TData
{
public:
  TData(){_mess=0; _par=0;};
  TData(unsigned short int index, TMeasurement * mess, double truth, double error, double weight, double * par, unsigned short int n_par,
        double const(*func)(TMeasurement *const,double *const),
	double const(*err)(double *const,TMeasurement *const,double const))
  : _index(index),_mess(mess),_truth(truth),_error(error),_weight(weight),_par(par),_n_par(n_par),_func(func),_err(err){};
  virtual ~TData(){
    delete _mess;
  };
  TMeasurement *GetMess() const { return _mess;};
  virtual double GetParametrizedMess() const {return _func(_mess,_par);};
  virtual double GetParametrizedMess(double *const paramess) const { // For derivative calculation
    TMeasurement m(_mess);
    m.pt = paramess[0];
    return _func(&m,_par);
  };
  virtual double GetParametrizedErr(double *const paramess) const { return _err(paramess,_mess,_error);};
  virtual double GetParametrizedErr2(double *const paramess){ 
    double error = GetParametrizedErr(paramess);
    return error *error;
  };
  double GetTruth() const { return _truth;};
  virtual double GetScale() const {return GetTruth();};//flatten spectrum w.r.t. this
  virtual void UpdateError(){};
  double GetError() const { return _error;};
  double GetWeight() const { return _weight;};
  void   SetWeight(const double & weight) { _weight=weight;};
  DataType GetType() const {return _type;};
  void SetType(DataType type) {_type=type;};
  unsigned short int GetIndex(){return _index;};
  virtual const std::vector<TData*>& GetRef() = 0;
  virtual double chi2() const {return 0.;};
  //virtual double chi2() const = 0;
  virtual double chi2_fast(double * temp_derivative1, double * temp_derivative2, double const epsilon) const = 0;
  double * GetPar(){return _par;};
  unsigned short int GetNumberOfPars() const {return _n_par;};
  virtual void ChangeParAddress(double* oldpar, double* newpar) { _par += newpar - oldpar;}

  static unsigned int total_n_pars;

  static double (*ScaleResidual)(double z2);         // Set to one of the following functions to scale the normalized residual z2 = chi2/weight in chi2() or chi2_fast(): 
  static double ScaleNone(double z2){ return z2; }  // No scaling of residuals in chi2() or chi2_fast() (default)
  static double ScaleCauchy(double z2);	             // Scaling of residuals with Cauchy-Function in chi2() or chi2_fast()
  static double ScaleHuber(double z2);               // Scaling of residuals with Huber-Function in chi2() or chi2_fast()


protected:
  unsigned short int _index; //limited from 0 to 65535
  TMeasurement *_mess;
  double _truth, _error, _weight;
  double *_par;
  unsigned short int _n_par; //limited from 0 to 65535
  double const(*_func)(TMeasurement *const x, double *const par);
  double const(*_err)(double *const x, TMeasurement *const x_original, double const error);
  DataType _type;
};

//data class for data providing one truth and one messurement, 
//e.g. track-tower
class TData_TruthMess : public TData
{
public:
  TData_TruthMess(unsigned short int index, TMeasurement * mess, double truth, double error, double weight, double * par, unsigned short int n_par,
        double const(*func)(TMeasurement *const,double *const),
		  double const(*err)(double *const,TMeasurement *const,double const))
  : TData(index, mess, truth, error, weight, par, n_par, func, err ){_type=TrackTower;};

  virtual const std::vector<TData*>& GetRef() { 
    resultcache.clear();	
    resultcache.push_back( this );
    return resultcache;
  };


//  virtual double GetParametrizedErr(double *const paramess) const{ 	 
//    double pmess; 	 
//    if(std::abs(_mess->eta) < 3.0) 	 
//      pmess =  paramess[0] * _mess->E / (_mess->pt * _mess->pt) * (_mess->HadF + _mess->OutF); //Et->E hadronic 	 
//    else 	 
//      pmess =  paramess[0] * (_mess->E / _mess->pt);  //Et->E 	 
//    return _err(&pmess,0,0) * _mess->pt / _mess->E; 	 
//  };     //search 	 
	 

  virtual double chi2() const{ 
    double new_mess  = GetParametrizedMess();
    double new_error = GetParametrizedErr(&new_mess);
    return GetWeight() * (*TData::ScaleResidual)( (_truth-new_mess)*(_truth-new_mess)/(new_error*new_error) );
  };
  virtual double chi2_fast(double * temp_derivative1, double*  temp_derivative2, double const epsilon) const;
  
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
        	      double const(*func)(TMeasurement *const,double *const),
		      double const(*err)(double *const,TMeasurement *const,double const),
		      TMeasurement *mess)
  : TData_TruthMess(index, mess, truth, error, weight, par, n_par, func, err){_type=GammaJet; };
  virtual ~TData_TruthMultMess() {
    for (std::vector<TData*>::const_iterator it=_vecmess.begin();
	 it!=_vecmess.end(); ++it)
      delete *it;
    _vecmess.clear();	
  };
  void AddMess(TData_TruthMess * m){ _vecmess.push_back(m);};

  virtual double GetParametrizedMess() const{
    double tower_pt_sum=0.0;
    for (std::vector<TData*>::const_iterator it=_vecmess.begin();
  	 it!=_vecmess.end(); ++it){
      tower_pt_sum += (*it)->GetParametrizedMess(); // Sum of tower Pt
    }
    TJet jet(_mess);
    jet.pt = tower_pt_sum;

    return _func(&jet, _par);
  };
  virtual double GetParametrizedErr() const { //returns total jet error and mot only the non-tower part like _err
    double error = 0, new_mess=0, sum_mess=0, new_error=0,sum_error2=0;
    for (std::vector<TData*>::const_iterator it=_vecmess.begin();
	 it!=_vecmess.end(); ++it) {
      new_mess    = (*it)->GetParametrizedMess();
      sum_mess   += new_mess;
      new_error   = (*it)->GetParametrizedErr(&new_mess);
      sum_error2 += new_error * new_error;
    }
    TJet jet(_mess);
    jet.pt    = sum_mess;
    new_mess  = _func(&jet, _par);
    new_error =  _err( &new_mess, _mess, _error );
    sum_error2 += new_error * new_error;
    
    error = sqrt(sum_error2);
    return error;
  };

  virtual double GetParametrizedMess(double *const paramess) const { // For derivative calculation
    TJet jet(_mess);
    jet.pt = paramess[0];
    return _func(&jet,_par);
  };

  virtual double chi2() const{ 
    double weight = GetWeight(); 
    double new_mess, new_error;
    double sum_mess = 0.;
    double sum_error2 = 0.;
    for (std::vector<TData*>::const_iterator it=_vecmess.begin();
	 it!=_vecmess.end(); ++it) {
      new_mess    = (*it)->GetParametrizedMess();
      sum_mess   += new_mess;
      new_error   = (*it)->GetParametrizedErr(&new_mess);
      sum_error2 += new_error * new_error;
    }
    TJet jet(_mess);
    jet.pt    = sum_mess;
    new_mess  = _func(&jet, _par);
    new_error =   _err( &new_mess, _mess, _error );
    return (new_error!=0. ? weight*(*TData::ScaleResidual)( (_truth-new_mess)*(_truth-new_mess)/(sum_error2 + new_error*new_error) ):0.0);
  };
  virtual double chi2_fast(double * temp_derivative1, double*  temp_derivative2, double const epsilon) const;
  virtual const std::vector<TData*>& GetRef() {return _vecmess;};
  virtual void ChangeParAddress(double* oldpar, double* newpar) { 
    TData::ChangeParAddress(oldpar,newpar);
    for (std::vector<TData*>::iterator it=_vecmess.begin();
	 it !=_vecmess.end(); ++it) { (*it)->ChangeParAddress(oldpar,newpar);}
  }
protected:  
  std::vector<TData*> _vecmess; 
  double JetError2;
};

//virtual data class for data providing only multiple messurements
//NOT DIRECTLY USED, SINCE "combine" FUNCTION IS NOT YET DEFINED HERE!
//Base class for PT-balance and Inv-Mass !!!
class TData_MessMess : public TData_TruthMultMess
{
public:
  TData_MessMess(unsigned short int index, double * dir, double truth, double error, 
                 double weight, double * par, unsigned short int n_par,
       	         double const(*func)(TMeasurement *const,double *const),
		 double const(*err)(double *const,TMeasurement *const,double const),
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
  virtual double chi2() const{ 
    double sum_error2=0.0, new_error, new_mess;
    double weight = GetWeight();
    
    for (std::vector<TData_MessMess*>::const_iterator it=_m2.begin();
         it!=_m2.end(); ++it) {
      new_error   = (*it)->GetParametrizedErr();
      sum_error2 += new_error * new_error;
    }
    new_error = GetParametrizedErr();
    new_mess  = GetMessCombination();

    return (sum_error2!=0 ? weight*(*TData::ScaleResidual)( (_truth-new_mess)*(_truth-new_mess)/(sum_error2 + new_error*new_error) ) : 0.0) ;
  };
  virtual double chi2_fast(double * temp_derivative1, double*  temp_derivative2, double const epsilon) const;
  virtual void ChangeParAddress(double* oldpar, double* newpar) { 
    TData::ChangeParAddress(oldpar,newpar);
    for (std::vector<TData*>::iterator it=_vecmess.begin();  it !=_vecmess.end(); ++it) 
      (*it)->ChangeParAddress(oldpar,newpar);
    for (std::vector<TData_MessMess*>::const_iterator it=_m2.begin(); it!=_m2.end(); ++it)
      (*it)->ChangeParAddress(oldpar,newpar);
  }
  virtual double GetMessCombination()const{ return combine(); }; // for plotting
  virtual double * GetDirection() const{ return _direction; };
  virtual std::vector<TData_MessMess*>const * GetSecondaryJets() const{return &_m2;};
   
protected:
  virtual double combine() const {return 0.;};
  std::vector<TData_MessMess*> _m2;
  double * _direction;
};

//e.g. jet-jet, jet-jet-jet, etc..
class TData_PtBalance : public TData_MessMess
{
public:
  TData_PtBalance(unsigned short int index, double * dir, double truth, double error, 
		  double weight, double * par, unsigned short int n_par,
        	  double const(*func)(TMeasurement *const,double *const),
		  double const(*err)(double *const,TMeasurement *const,double const),
		  TMeasurement *mess)
  : TData_MessMess(index, dir, truth, error, weight, par, n_par, func, err, mess){_type=PtBalance;};
  virtual ~TData_PtBalance(){};
  virtual double GetScale() const{
    double scale = GetMess()->pt;
    for (std::vector<TData_MessMess*>::const_iterator it=_m2.begin();it!=_m2.end(); ++it)
      scale += (*it)->GetMess()->pt;
    return scale/2.;
  }
  virtual void UpdateError() {
    double totalsum = GetParametrizedMess();
    for (std::vector<TData_MessMess*>::const_iterator it=_m2.begin();
	 it!=_m2.end(); ++it) {
      totalsum += (*it)->GetParametrizedMess();
    }
    
    double scale = totalsum / 2; //should be changed to ptsum (scalar) of the part projected to leading jet axis devided by 2 for n>2 n-jets     (also chi2() and chi2fast())

    double new_mess = GetParametrizedMess();
    double sum = totalsum - new_mess; 
    double new_error = GetParametrizedErr();
    new_error *= sum / (scale * scale);
    double sum_error2 = new_error * new_error;
    for (std::vector<TData_MessMess*>::const_iterator it=_m2.begin();
	 it!=_m2.end(); ++it) {
      new_mess    = (*it)->GetParametrizedMess();
      new_error   = (*it)->GetParametrizedErr();
      sum = totalsum - new_mess;
      new_error *= sum / (scale * scale);
      sum_error2 += new_error * new_error;
    }
    if( (GetParametrizedMess()>0) && (totalsum>GetParametrizedMess()))   //update only if new value is reasonable (Pt > 0)
    _error = sqrt(sum_error2);  
};

  virtual double chi2() const{ 
    double new_mess;
    double weight = GetWeight();
      
    double parascale = GetParametrizedMess();
    for (std::vector<TData_MessMess*>::const_iterator it=_m2.begin();
	 it!=_m2.end(); ++it) {
      parascale += (*it)->GetParametrizedMess();
    }
    
    //should be changed to ptsum (scalar) of the part projected to leading jet axis devided by 2 for n>2 n-jets           same for parascale and in ChiFast...
    parascale /= 2;

    /*
    for (std::vector<TData_MessMess*>::const_iterator it=_m2.begin();
         it!=_m2.end(); ++it) {
      //new_mess    = (*it)->GetParametrizedMess();
      new_mess    = (*it)->GetMess()->pt;	 //const error
      new_error   = (*it)->GetParametrizedErr(&new_mess);
      sum = totalsum -  (*it)->GetMess()->pt;   //new_mess?
      new_error *= sum / (scale * scale);
      sum_error2 += new_error * new_error;
    }
    new_mess  = GetMess()->pt;          //const error
    sum = totalsum - GetMess()->pt;   //new_mess?
    new_error = GetParametrizedErr( &new_mess );
    new_error *= sum / (scale * scale);   //
    */
    new_mess  = combine() / parascale;

    return (_error!=0 ? weight*(*TData::ScaleResidual)( (_truth-new_mess)*(_truth-new_mess)/(_error * _error) ) : 0.0) ;
  };
  virtual double chi2_fast(double * temp_derivative1, double*  temp_derivative2, double const epsilon) const;
protected:
  virtual double combine() const;  
};

//e.g. top-Mass, W-mass, etc..
class TData_InvMass2 : public TData_MessMess
{
public:
  TData_InvMass2(unsigned short int index, double * dir, double truth, double error, double weight,
		double * par, unsigned short int n_par,
        	double const(*func)(TMeasurement *const,double *const),
		double const(*err)(double *const,TMeasurement *const,double const),
		TMeasurement *mess)
  : TData_MessMess(index, dir, truth, error, weight, par, n_par, func, err, mess) { _type=InvMass; };
  virtual ~TData_InvMass2(){};
 protected:
  virtual double combine() const;  
};

//data class to limit a parameter
class TData_ParLimit : public TData
{
 public:
  TData_ParLimit(unsigned short int index, TMeasurement *mess, double error,
                 double *par,double const(*func)(TMeasurement *const,double *const)) 
  : TData(index,mess,0,error,1.0,par,1,func,0){ _type=ParLimit; };
    
    virtual const std::vector<TData*>& GetRef() { 
      _cache.clear();
      _cache.push_back(this);
      return _cache;
    };
    
    virtual double GetParametrizedErr(double *const paramess) const{ 
      return _error;
    };
        
    virtual double chi2() const{  
      double new_mess  = GetParametrizedMess();
      double new_error = GetParametrizedErr(&new_mess);
      return new_mess * new_mess / (new_error * new_error);
    };
    virtual double chi2_fast(double* temp_derivative1, 
			     double* temp_derivative2, double const epsilon) const;
    
 private:
    static std::vector<TData*> _cache;
};

#endif
