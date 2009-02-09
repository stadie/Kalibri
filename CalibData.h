//
// Original Author:  Christian Autermann
//         Created:  Wed Jul 18 13:54:50 CEST 2007
// $Id: CalibData.h,v 1.57 2009/01/16 08:46:40 stadie Exp $
//
#ifndef CalibData_h
#define CalibData_h

#include <iostream>
using namespace std;
#include <vector> 
#include <cmath>
#include <cassert>

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
  TMeasurement(double Et,double EmEt,double HadEt,double OutEt,double E,
	       double eta,double phi)
    : pt(Et),EMF(EmEt),HadF(HadEt),OutF(OutEt),E(E),eta(eta),phi(phi) {}
  TMeasurement(TMeasurement* m):pt(m->pt),EMF(m->EMF),HadF(m->HadF),OutF(m->OutF),
                                E(m->E),eta(m->eta),phi(m->phi){};
  //all common variables
  double pt;
  double EMF;
  double HadF;
  double OutF;
  double E;
  double eta;
  double phi;
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
  enum Flavor{ gluon=0, uds=1, c=2, b=3 };
  TJet():TMeasurement(){}; 
  TJet(double Et,double EmEt,double HadEt,double OutEt,double E,double eta,
       double phi, Flavor flavor, double genPt, double ZSPcor, double JPTcor, double L2cor, double L3cor)
    : TMeasurement(Et,EmEt,HadEt,OutEt,E,eta,phi),flavor(flavor), genPt(genPt) {};
  TJet(TMeasurement* j):TMeasurement(j){};
  TJet(TJet* j):TMeasurement(j){/*further initialization*/};
//variables specific only to jets (i.e. mass)
  Flavor flavor;
  double genPt;
  double ZSPcor;
  double JPTcor;
  double L2cor;
  double L3cor;
};

class TTrack : public TMeasurement
{
public:
  TTrack():TMeasurement(){};
  TTrack(TMeasurement* tr):TMeasurement(tr){};
  TTrack(TTrack* tr):TMeasurement(tr){/*further initialization*/};
//variables specific only to Tracks
  int TrackId;
  int TowerId;
  double DR;
  double DRout;
  double etaOut;
  double phiOut;
  double EM1;
  double EM5;
  double Had1;
  double Had5;
  double TrackChi2;
  int NValidHits;
  bool TrackQualityT;
  double MuDR;
  double MuDE;
};

//interface to Data
class TData
{
public:
  virtual ~TData() {}
  virtual TMeasurement *GetMess() const = 0;
  virtual double GetTruth() const = 0;
  virtual double GetParametrizedMess() const = 0;

  virtual void ChangeParAddress(double* oldpar, double* newpar) = 0;
  virtual DataType GetType() const = 0;
  virtual double GetWeight() const = 0;
  
  virtual double chi2() const = 0;
  virtual double chi2_fast(double * temp_derivative1, double * temp_derivative2, double const epsilon) const = 0;
  virtual void UpdateError() = 0;

  static double (*ScaleResidual)(double z2);         // Set to one of the following functions to scale the normalized residual z2 = chi2/weight in chi2() or chi2_fast(): 
  static double ScaleNone(double z2){ return z2; }  // No scaling of residuals in chi2() or chi2_fast() (default)
  static double ScaleCauchy(double z2);	             // Scaling of residuals with Cauchy-Function in chi2() or chi2_fast()
  static double ScaleHuber(double z2);               // Scaling of residuals with Huber-Function in chi2() or chi2_fast()
};

//virtual data base class -> not directly used!
class TAbstractData : public TData
{
public:
  TAbstractData() : TData() {_par=0;};
  TAbstractData(unsigned short int index, TMeasurement * mess, double truth, double error, double weight, double * par, unsigned short int n_par,
        double (*func)(const TMeasurement*, const double*),
	double (*err)(const double *,const TMeasurement *,double))
  : _index(index), _mess(mess),_truth(truth),_error(error),_weight(weight),_par(par),_n_par(n_par),_func(func),_err(err){};
  virtual ~TAbstractData(){
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
  virtual double GetParametrizedErr() const { return _error;};
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
  virtual const std::vector<TAbstractData*>& GetRef() = 0;
  virtual const std::vector<TAbstractData*>& GetRefTrack() = 0;
  virtual double chi2() const {return 0.;};
  //virtual double chi2() const = 0;
  virtual double chi2_fast(double * temp_derivative1, double * temp_derivative2, double const epsilon) const = 0;
  double * GetPar(){return _par;};
  unsigned short int GetNumberOfPars() const {return _n_par;};
  virtual void ChangeParAddress(double* oldpar, double* newpar) { _par += newpar - oldpar;}
  virtual bool GetTrackuse() {return true;}

  static unsigned int total_n_pars;

protected:
  unsigned short int _index; //limited from 0 to 65535
  TMeasurement *_mess;
  double _truth, _error, _weight;
  double *_par;
  double _invisError;
  unsigned short int _n_par; //limited from 0 to 65535
  double (*_func)(const TMeasurement *x, const double *par);
  double (*_err)(const double *x, const TMeasurement *x_original, double error);
  DataType _type;
};

//data class for data providing one truth and one messurement, 
//e.g. track-tower
class TData_TruthMess : public TAbstractData
{
public:
  TData_TruthMess(unsigned short int index,  TMeasurement * mess, double truth, double error, double weight, double * par, unsigned short int n_par, double (*func)(const TMeasurement *,const double*),  double (*err)(const double*,const TMeasurement*,double))
  : TAbstractData(index, mess, truth, error, weight, par, n_par, func, err ){_type=TrackTower;_invisError=0;};

  virtual const std::vector<TAbstractData*>& GetRef() { 
    resultcache.clear();	
    resultcache.push_back( this );
    return resultcache;
  };
  
  //only a dummy
  virtual const std::vector<TAbstractData*>& GetRefTrack()  { 
    resultcache.clear();	
    resultcache.push_back( this );
    return resultcache;
  };


  virtual double chi2() const{ 
    double new_mess  = GetParametrizedMess();
    double new_error = GetParametrizedErr(&new_mess);
    return GetWeight() * (*TData::ScaleResidual)( (_truth-new_mess)*(_truth-new_mess)/(new_error*new_error) );
  };
  virtual double chi2_fast(double * temp_derivative1, double*  temp_derivative2, double const epsilon) const;
  
private:
  static std::vector<TAbstractData*> resultcache;
};


//data class for data providing one truth and multiple messurements, 
//e.g. gamma-jet or track-cluster
class TData_TruthMultMess : public TData_TruthMess
{
public:
  TData_TruthMultMess(unsigned short int index, double truth, double error, 
		      double weight, double * par, unsigned short int n_par,
        	      double (*func)(const TMeasurement *,const double*),
		      double (*err)(const double*,const TMeasurement*,double),
		      TMeasurement *mess)
  : TData_TruthMess(index,  mess, truth, error, weight, par, n_par, func, err){
    _type=GammaJet; trackuse = false;_invisError=0;};
  virtual ~TData_TruthMultMess() {
    for (std::vector<TAbstractData*>::const_iterator it=_vecmess.begin();
	 it!=_vecmess.end(); ++it)
      delete *it;
    _vecmess.clear();	
  };
  void AddMess(TData_TruthMess * m){ _vecmess.push_back(m);};
  void AddTrack(TData_TruthMess * m){ _vectrack.push_back(m);};

  void UseTracks(bool use){   //check if tracks shall be used in this jet:
    if((fabs(_mess->eta) < 2.1) &&  ( _vectrack.size() > 0) && use)      //@ eta > 2.1or eta < -2.1 parts of the cone are outside the tracker. 
      trackuse = true;
  };
  virtual bool GetTrackuse(){return trackuse;};

  virtual void UpdateError() {
    if(trackuse)
      {
	double CaloRest = _mess->pt;
	bool IsMuon;
	double ConeRadius = 0.5;
	for (std::vector<TAbstractData*>::const_iterator it=_vectrack.begin();
	     it!=_vectrack.end(); ++it) {
	  //if( ((TTrack*)(*it)->GetMess())->TrackChi2 > 6  || ((TTrack*)(*it)->GetMess())->NValidHits < 9) // || TrackQuality != 1)
	  if( !((TTrack*)(*it)->GetMess())->TrackQualityT)
	    continue;  // qualityTracks = false;

	  if( (*it)->GetMess()->pt > 100) continue;

	  TTrack* temp = (TTrack*)(*it)->GetMess();
	  if(temp->TrackId == 13)  IsMuon = true;
	  else                     IsMuon = false;
	  if(temp->DR < ConeRadius) 
	    {
	      if(temp->DRout < ConeRadius)
		{
		  if(!IsMuon) {
		    CaloRest -= (*it)->GetParametrizedMess();   //this depends on start values. Error will be updated, but look for a more stable way  
		  }
		}
	    }
	}
	std::vector<TAbstractData*>::const_iterator tower=_vecmess.begin();
	if(CaloRest > 0)
	  _invisError = (*tower)->GetParametrizedErr(&CaloRest);  //Rest Error same as tower error, save this invisible energy error
      }
  }

  virtual double GetParametrizedMess() const{
    double result = 0;
    if(trackuse) 
      { 
	result = GetParametrizedTrackMess();
      }
    else
      {    //no tracks are used
	double tower_pt_sum=0.0;
	for (std::vector<TAbstractData*>::const_iterator it=_vecmess.begin();
	     it!=_vecmess.end(); ++it){
	  tower_pt_sum += (*it)->GetParametrizedMess(); // Sum of tower Pt
	}
	TJet jet(_mess);
	jet.pt = tower_pt_sum;
	result = _func(&jet, _par);
      }
    return result;
  };

  virtual double GetParametrizedErr() const { //returns total jet error and not only the non-tower part like _err
    double error = 0, new_mess=0, sum_mess=0, new_error=0,sum_error2=0;
    if(trackuse)
      {
	for (std::vector<TAbstractData*>::const_iterator it=_vectrack.begin();
	     it!=_vectrack.end(); ++it) {

	  if( !((TTrack*)(*it)->GetMess())->TrackQualityT)  continue;
	  if( (*it)->GetMess()->pt > 100) continue;

	  new_mess =  (*it)->GetMess()->pt;
	  new_error   = (*it)->GetParametrizedErr(&new_mess);
	  sum_error2 += new_error * new_error;
	}
	sum_error2 += _invisError * _invisError;
      }
    else
      {
      for (std::vector<TAbstractData*>::const_iterator it=_vecmess.begin();
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
	}
    error = sqrt(sum_error2);
    return error;
  };

  virtual double GetParametrizedMess(double *const paramess) const { // For derivative calculation
    double result = 0;
    //When tracks are used, GetParametrizedTrackMess() should have been taken, i.e. no tower correction should have been used before 
    if(trackuse) 	result = GetParametrizedTrackMess();
    else {
      TJet jet(_mess);
      jet.pt = paramess[0];
      result =  _func(&jet,_par);
    }
    return result;
  };

  virtual double GetParametrizedTrackMess() const{
    double JetPt = 0, CaloRest, CaloTrackPt;  
    int NoUsedTracks = 0; 
    bool IsMuon;
    const double ConeRadius = 0.5;      //Jet Cone Radius should not be hard coded or must be changed if Radius != 0.5
    const double MIPsignal = 4;      
    CaloRest = _mess->pt;
    TJet* myTJet =  (TJet*)(GetMess());
    CaloRest *= myTJet->ZSPcor;   //ZSP correction
    for (std::vector<TAbstractData*>::const_iterator it=_vectrack.begin();
	 it!=_vectrack.end(); ++it) {
      //if( ((TTrack*)(*it)->GetMess())->TrackChi2 > 6  || ((TTrack*)(*it)->GetMess())->NValidHits < 9) // || TrackQuality != 1)
      if( !((TTrack*)(*it)->GetMess())->TrackQualityT)
	continue;  // qualityTracks = false;

      if( (*it)->GetMess()->pt > 100) continue;

      NoUsedTracks++;
      TTrack* temp = (TTrack*)(*it)->GetMess();

      CaloTrackPt = (*it)->GetParametrizedMess();  //Expected Signal of track in Calorimeter

      //-Track Reco Efficiency betrachten

      if(temp->TrackId == 13)  IsMuon = true;
      else                     IsMuon = false;
      if(temp->DR < ConeRadius) 
	{
	  if(temp->DRout < ConeRadius)
	    {
	      if(!IsMuon) {
		JetPt += temp->pt;
		CaloRest -= CaloTrackPt;
	      }
	    }
	  else //Out of Cone
	    {//commet out if correction should be done to genJet level
	      JetPt += temp->pt;
	    }
	}
      else
	{
	  if(IsMuon)     CaloRest -= MIPsignal;
	  else           CaloRest -= CaloTrackPt;
	}
    }
    //cout<<"TrackPt: "<<JetPt<<"   Truth: "<<_truth<<"   rel Err: "<<(JetPt - _truth)/_truth<<"    Rest: "<<CaloRest<<"  relErrorIncRest: "<<(JetPt + CaloRest - _truth)/_truth<<"    #Tracks: "<<_vectrack.size()<<"   #UsedTracks: "<<NoUsedTracks<<endl;

    if(CaloRest > 0) {
      TJet jet(_mess);
      jet.pt = CaloRest;
      jet.E = -1000;  //signal that this is only the calo rest of a track jet!
      JetPt += _func(&jet, _par);  //ein anderer Parameter als bei der normalen korrektur nehmen (da hier nur neutr. Had + other rest)
    }

    return JetPt;
  };

  virtual double chi2() const{ 
    double weight = GetWeight(); 
    double new_mess, new_error;
    double sum_mess = 0.;
    double sum_error2 = 0.;
    if(trackuse) {
      new_error = GetParametrizedErr();
      new_mess = GetParametrizedTrackMess();
    }
    else{
      for (std::vector<TAbstractData*>::const_iterator it=_vecmess.begin();
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
    }
    return (new_error!=0. ? weight*(*TData::ScaleResidual)( (_truth-new_mess)*(_truth-new_mess)/(sum_error2 + new_error*new_error) ):0.0);
  };
  virtual double chi2_fast(double * temp_derivative1, double*  temp_derivative2, double const epsilon) const;
  virtual const std::vector<TAbstractData*>& GetRef() {return _vecmess;};
  virtual const std::vector<TAbstractData*>& GetRefTrack() {return _vectrack;};
  virtual void ChangeParAddress(double* oldpar, double* newpar) { 
    TAbstractData::ChangeParAddress(oldpar,newpar);
    for (std::vector<TAbstractData*>::iterator it=_vecmess.begin();
	 it !=_vecmess.end(); ++it) { (*it)->ChangeParAddress(oldpar,newpar);} 
    if(trackuse) {
      for (std::vector<TAbstractData*>::iterator it = _vectrack.begin() ;
	   it !=_vectrack.end() ; ++it) { (*it)->ChangeParAddress(oldpar,newpar);} 
    }
  }
protected:  
  std::vector<TAbstractData*> _vecmess; 
  std::vector<TAbstractData*> _vectrack; 
  bool trackuse;
  //double JetError2;
};

//virtual data class for data providing only multiple messurements
//NOT DIRECTLY USED, SINCE "combine" FUNCTION IS NOT YET DEFINED HERE!
//Base class for PT-balance and Inv-Mass !!!
class TData_MessMess : public TData_TruthMultMess
{
public:
  TData_MessMess(unsigned short int index, double * dir, double truth, double error, 
                 double weight, double * par, unsigned short int n_par,
       	         double (*func)(const TMeasurement *,const double*),
		 double (*err)(const double*,const TMeasurement*,double),
		 TMeasurement *mess)
  : TData_TruthMultMess(index, truth, error, weight, par, n_par, func, err,mess),
    _direction(dir){_type=MessMess; trackuse = false;_invisError=0;};
  virtual ~TData_MessMess() {
    for (std::vector<TData_MessMess*>::const_iterator it=_m2.begin();
	 it!=_m2.end(); ++it)
      delete *it;
    _m2.clear();
    delete [] _direction;	
  };
  virtual void AddNewMultMess(TData_MessMess * m2 ){assert(m2->_m2.empty());_m2.push_back(m2);};
  void ClearMultMess() { _m2.clear();}
  unsigned MultMessSize(){return _m2.size();};
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
    TData_TruthMultMess::ChangeParAddress(oldpar,newpar);
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
        	  double (*func)(const TMeasurement *,const double*),
		  double (*err)(const double*,const TMeasurement*,double),
		  TMeasurement *mess)
  : TData_MessMess(index, dir, truth, error, weight, par, n_par, func, err, mess){_type=PtBalance; trackuse = false;_invisError=0;};
  virtual ~TData_PtBalance(){};
  virtual double GetScale() const{
    double scale = GetMess()->pt;
    for (std::vector<TData_MessMess*>::const_iterator it=_m2.begin();it!=_m2.end(); ++it)
      scale += (*it)->GetMess()->pt;
    return scale/2.;
  }
  virtual void UpdateError() {
    if(trackuse) //update error on invisible energy 1st jet
      {
	double CaloRest = _mess->pt;
	bool IsMuon;
	double ConeRadius = 0.5;
	for (std::vector<TAbstractData*>::const_iterator it=_vectrack.begin();
	     it!=_vectrack.end(); ++it) {
	  //if( ((TTrack*)(*it)->GetMess())->TrackChi2 > 6  || ((TTrack*)(*it)->GetMess())->NValidHits < 9) // || TrackQuality != 1)
	  if( !((TTrack*)(*it)->GetMess())->TrackQualityT)
	    continue;  // qualityTracks = false;

	  if( (*it)->GetMess()->pt > 100) continue;

	  TTrack* temp = (TTrack*)(*it)->GetMess();
	  if(temp->TrackId == 13)  IsMuon = true;
	  else                     IsMuon = false;
	  if(temp->DR < ConeRadius) 
	    {
	      if(temp->DRout < ConeRadius)
		{
		  if(!IsMuon) {
		    CaloRest -= (*it)->GetParametrizedMess();   //this depends on start values. Error will be updated, but look for a more stable way  
		  }
		}
	    }
	}
	std::vector<TAbstractData*>::const_iterator tower=_vecmess.begin();
	if(CaloRest > 0)
	  _invisError = (*tower)->GetParametrizedErr(&CaloRest);  //Rest Error same as tower error, save this invisible energy error
      }

    double totalsum = GetParametrizedMess();
    for (std::vector<TData_MessMess*>::const_iterator it=_m2.begin();
	 it!=_m2.end(); ++it) {
      (*it)->UpdateError();  //update error on invisible energy other jets
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
        	double (*func)(const TMeasurement *,const double*),
		double (*err)(const double*,const TMeasurement*,double),
		TMeasurement *mess)
  : TData_MessMess(index, dir, truth, error, weight, par, n_par, func, err, mess) { _type=InvMass; trackuse = false;_invisError=0; };
  virtual ~TData_InvMass2(){};
 protected:
  virtual double combine() const;  
};

//data class to limit a parameter
class TData_ParLimit : public TAbstractData
{
 public:
  TData_ParLimit(unsigned short int index, TMeasurement *mess, double error,
                 double *par,double (*func)(const TMeasurement *,const double*))
  : TAbstractData(index,mess,0,error,1.0,par,1,func,0){ _type=ParLimit;_invisError=0;};
    
    virtual const std::vector<TAbstractData*>& GetRef() { 
      _cache.clear();
      _cache.push_back(this);
      return _cache;
    };
    //only a dummy
    virtual const std::vector<TAbstractData*>& GetRefTrack() { 
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
    static std::vector<TAbstractData*> _cache;
};

#endif
