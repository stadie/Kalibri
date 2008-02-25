//
// Original Author:  Christian Autermann
//         Created:  Wed Jul 18 13:54:50 CEST 2007
// $Id: CalibData.C,v 1.10 2008/02/25 10:07:45 stadie Exp $
//
#include "CalibData.h"
#include "map"
//#include <iostream>//cout


unsigned int TData::total_n_pars = 0;
//unsigned int TData::total_tower_pars = 0;
//unsigned int TData::total_jet_pars = 0;

std::vector<TData*> TData_TruthMess::resultcache = std::vector<TData*>(1);

double TData_TruthMess::chi2_fast(double* temp_derivative1, double*  temp_derivative2, 
			        double epsilon)
{ 
  double new_mess  = GetParametrizedMess();
#ifndef __FastErrorCalculation
  double new_error = GetParametrizedErr(&new_mess);
#else
  double new_error = _error;
  if (GetMess()[0]!=0.) new_error = _error*new_mess/GetMess()[0];
#endif   
  double weight = GetWeight();
  double new_chi2  = weight*(_truth-new_mess)*(_truth-new_mess)/(new_error*new_error);
  
  double dmess_dp, derror_dp;
  unsigned idx = _index*_n_par; //_index==bin; idx==bin*Free_parameters_per_bin
  for (unsigned i=0; i<idx; ++i){
    temp_derivative2[i]+=new_chi2;
    temp_derivative1[i]+=new_chi2;
  }
  for (unsigned i=idx; i<idx+_n_par; ++i){
    double oldpar = _par[i-idx];
    _par[i-idx]  += epsilon;
    dmess_dp  = GetParametrizedMess();
#ifndef __FastErrorCalculation
    derror_dp = GetParametrizedErr(&dmess_dp);
#else
    if (GetMess()[0]!=0.) derror_dp = _error*dmess_dp/GetMess()[0];
    else derror_dp = _error;
#endif   
    temp_derivative2[i]+= weight*(_truth-dmess_dp)*(_truth-dmess_dp)/(derror_dp*derror_dp);

    _par[i-idx]  = oldpar - epsilon;
    dmess_dp  = GetParametrizedMess();
#ifndef __FastErrorCalculation
    derror_dp = GetParametrizedErr(&dmess_dp);
#else
    if (GetMess()[0]!=0.) derror_dp = _error*dmess_dp/GetMess()[0];
    else derror_dp = _error;
#endif   
    temp_derivative1[i]+= weight*(_truth-dmess_dp)*(_truth-dmess_dp)/(derror_dp*derror_dp);
    _par[i-idx]  = oldpar;
  }
  for (unsigned i=idx+_n_par; i< total_n_pars; ++i){
    temp_derivative2[i]+=new_chi2;
    temp_derivative1[i]+=new_chi2;
  }
  return new_chi2;
};

double TData_TruthMultMess::chi2_fast(double* temp_derivative1, double*  temp_derivative2,
				      double epsilon) {
  double sum_mess=0.0, sum_error2=0.0, new_error, new_error2, new_mess, new_chi2,
         dmess_dp, derror_dp;
  double sm1[total_n_pars],sm2[total_n_pars],se1[total_n_pars],se2[total_n_pars];//store sum_m & sum_e2 for df/dp_i
  double weight = GetWeight();
  for (unsigned i=0; i<total_n_pars; ++i){
    sm1[i] = 0.0;
    sm2[i] = 0.0;
    se1[i] = 0.0;
    se2[i] = 0.0;
  }
  unsigned idx;
  for (std::vector<TData*>::const_iterator it=_vecmess.begin();
       it!=_vecmess.end(); ++it) {
    new_mess    = (*it)->GetParametrizedMess();	 
    sum_mess   += new_mess;
//#ifndef __FastErrorCalculation
    new_error   = (*it)->GetParametrizedErr(&new_mess);
//#else
//    if ((*it)->GetMess()[0]!=0.) new_error = (*it)->GetError()*new_mess/(*it)->GetMess()[0];
//    else new_error = (*it)->GetError();
//#endif   
    new_error2  = new_error * new_error;
    sum_error2 += new_error2;
    for (unsigned i=0; i<total_n_pars; ++i){
      idx = (*it)->GetIndex()*(*it)->GetNumberOfPars();
      if (i>=idx && i<idx+(*it)->GetNumberOfPars()) {   
	double oldpar = (*it)->GetPar()[i-idx];
	(*it)->GetPar()[i-idx]  += epsilon;
	dmess_dp  =(*it)->GetParametrizedMess();
	sm2[i]   += dmess_dp;
//#ifndef __FastErrorCalculation
	new_error = (*it)->GetParametrizedErr(&dmess_dp);
//#else
//        if ((*it)->GetMess()[0]!=0.) new_error = (*it)->GetError()*dmess_dp/(*it)->GetMess()[0];
//	else new_error = (*it)->GetError();
//#endif   
	se2[i]   += new_error * new_error;
	
	(*it)->GetPar()[i-idx]  = oldpar - epsilon;
	dmess_dp  =(*it)->GetParametrizedMess();
	sm1[i]   += dmess_dp;
//#ifndef __FastErrorCalculation
	new_error = (*it)->GetParametrizedErr(&dmess_dp);
//#else
//        if ((*it)->GetMess()[0]!=0.) new_error = (*it)->GetError()*dmess_dp/(*it)->GetMess()[0];
//	else new_error = (*it)->GetError();
//#endif   
	se1[i]   += new_error * new_error;
	(*it)->GetPar()[i-idx] = oldpar;
      } else {
	sm1[i] += new_mess;
	sm2[i] += new_mess;
	se1[i] += new_error2;
	se2[i] += new_error2;
      }
    }
  } 
  new_mess  = _func( &sum_mess, _par);  
//#ifndef __FastErrorCalculation
  new_error = _err(&new_mess); 
//#else
//  if (_mess!=0) if (_mess[0]!=0.) new_error = _error*new_mess/_mess[0];
//  else new_error = _error;
//#endif   
  new_chi2  = weight*(_truth-new_mess)*(_truth-new_mess)/(sum_error2 + new_error*new_error);

  idx = _index; //@@to be fixed -> introduce a eta-phi binning for JES
  for (unsigned i=0; i<total_n_pars; ++i){
    if (i>=idx && i<idx+_n_par) continue;//considered below
    new_mess  = sm2[i];  //the measurement with modified parameter "i"
    dmess_dp  = _func(&new_mess,_par);//calc. the jet's energy
//#ifndef __FastErrorCalculation
    derror_dp = _err(&dmess_dp);
//#else
//    if (new_mess!=0.) derror_dp = _error*dmess_dp/new_mess;
//    else derror_dp = _error;
//#endif   
    temp_derivative2[i]+= weight*(_truth-dmess_dp)*(_truth-dmess_dp)/(se2[i] + derror_dp*derror_dp);
    // same for p_i-epsilon:
    new_mess  = sm1[i];  
    dmess_dp  = _func(&new_mess,_par);
//#ifndef __FastErrorCalculation
    derror_dp = _err(&dmess_dp);
//#else
//    if (new_mess!=0.) derror_dp = _error*dmess_dp/new_mess;
//    else derror_dp = _error;
//#endif   
    derror_dp = _err(&dmess_dp);
    temp_derivative1[i]+= weight*(_truth-dmess_dp)*(_truth-dmess_dp)/(se1[i] + derror_dp*derror_dp);
  }
  	
  for (unsigned i=idx; i<idx+_n_par; ++i){
    //ok, we have to change the jet's parametrization:
    double oldpar =  _par[i-idx];
    _par[i-idx]  += epsilon;
    dmess_dp  = _func(&sum_mess,_par);
//#ifndef __FastErrorCalculation
    derror_dp = _err(&dmess_dp);
//#else
//    if (sum_mess!=0.) derror_dp = _error*dmess_dp/sum_mess;
//    else derror_dp = _error;
//#endif
    temp_derivative2[i]+= weight*(_truth-dmess_dp)*(_truth-dmess_dp)/(se2[i] + derror_dp*derror_dp);

    _par[i-idx]  = oldpar - epsilon;
    dmess_dp  = _func(&sum_mess,_par);
//#ifndef __FastErrorCalculation
    derror_dp = _err(&dmess_dp);
//#else
//    if (sum_mess!=0.) derror_dp = _error*dmess_dp/sum_mess;
//    else derror_dp = _error;
//#endif
    temp_derivative1[i]+= weight*(_truth-dmess_dp)*(_truth-dmess_dp)/(se1[i] + derror_dp*derror_dp);
    _par[i-idx]  = oldpar;
  }

  return new_chi2;
};

//#include <iostream>

double TData_MessMess::chi2_fast(double * temp_derivative1, double*  temp_derivative2, 
			        double epsilon)
{
  double new_chi2, new_mess, new_error, sum_error2=0.0;

  double weight = GetWeight();

  //Get all tower parameter used in this event:  
  std::map<int,double*> tpars;
  for (std::vector<TData*>::const_iterator it=_vecmess.begin();
       it!=_vecmess.end(); ++it) 
     for (unsigned i=(*it)->GetIndex(); i<(*it)->GetIndex()+(*it)->GetNumberOfPars(); ++i)
	 tpars[ i ]= &((*it)->GetPar()[i]);
  for (std::vector<TData_MessMess*>::const_iterator mit=_m2.begin();
       mit!=_m2.end();++mit){
     std::vector<TData*>::const_iterator mitend=(*mit)->GetRef().end();  
     for (std::vector<TData*>::const_iterator it=(*mit)->GetRef().begin();
	it!=mitend; ++it) {
       for (unsigned i=(*it)->GetIndex(); i<(*it)->GetIndex()+(*it)->GetNumberOfPars(); ++i)
	   tpars[ i ]= &((*it)->GetPar()[i]);
       }
     }  
  //Add all jet parameters in this event:
  for (unsigned i=0; i<_n_par; ++i)
    tpars[ i+_index ]= &(_par[i]);
  for (std::vector<TData_MessMess*>::const_iterator mit=_m2.begin();
       mit!=_m2.end();++mit)
    for (unsigned i=0; i<(*mit)->GetNumberOfPars(); ++i)
	{tpars[ i+(*mit)->GetIndex() ]= &((*mit)->GetPar()[i]);
	}

  //This event's chi^2 for the current (unchanged) parameters:
  new_chi2 = chi2();

  //Calc. & Cache derivatives w.r.t. tower parameters:
  for (unsigned i=0; i<total_n_pars; ++i){
     // in case the ith parameter is used for this event, calculate 
     // derivative dchi/dpar[i] and cache it:
     if ( tpars.find(i)!=tpars.end() ) {
       double oldpar =  *tpars[i];
       *tpars[i] += epsilon;
       sum_error2 = 0.0;
       for (std::vector<TData_MessMess*>::const_iterator mit=_m2.begin();
            mit!=_m2.end(); ++mit) {
	  new_mess    = (*mit)->GetParametrizedMess();	 
//#ifndef __FastErrorCalculation
  	  new_error   = (*mit)->GetParametrizedErr(&new_mess);
//#else
//  	  if ((*mit)->GetMess()) if ((*mit)->GetMess()[0]!=0.) new_error = (*mit)->GetError()*new_mess/(*mit)->GetMess()[0];
//	  else new_error   = (*mit)->GetError();
//#endif
	  sum_error2 += new_error * new_error;
       }
       new_mess  = GetParametrizedMess();
//#ifndef __FastErrorCalculation
       new_error = _err(  &new_mess );
//#else
//       new_error = _error;
//#endif
       new_mess  = GetMessCombination();  
       temp_derivative2[ i ]+=weight*(_truth-new_mess)*(_truth-new_mess)/(sum_error2 + new_error*new_error);
       sum_error2 = 0.0;

       *tpars[i] = oldpar - epsilon;
       for (std::vector<TData_MessMess*>::const_iterator mit=_m2.begin();
            mit!=_m2.end(); ++mit) {
	  new_mess    = (*mit)->GetParametrizedMess();	 
//#ifndef __FastErrorCalculation
  	  new_error   = (*mit)->GetParametrizedErr(&new_mess);
//#else
//  	  if ((*mit)->GetMess()) if ((*mit)->GetMess()[0]!=0.) new_error = (*mit)->GetError()*new_mess/(*mit)->GetMess()[0];
//	  else new_error   = (*mit)->GetError();
//#endif
	  sum_error2 += new_error * new_error;
       }
       new_mess  = GetParametrizedMess();
//#ifndef __FastErrorCalculation
       new_error = _err(  &new_mess );
//#else
//       new_error = _error;
//#endif
       new_mess  = GetMessCombination();  
       temp_derivative1[ i ]+=weight*(_truth-new_mess)*(_truth-new_mess)/(sum_error2 + new_error*new_error);
       *tpars[i] = oldpar;
    }
    else 
    {
      temp_derivative1[ i ]+=new_chi2;
      temp_derivative2[ i ]+=new_chi2;
    }
  }
  return new_chi2;
}

double TData_PtBalance::combine(){
  double x, y, dummy = GetParametrizedMess();
  //here _direction is a normalized vector in eta-phi plane in pT direction.
  
  x = dummy * _direction[0];
  y = dummy * _direction[1];
  for (std::vector<TData_MessMess*>::const_iterator it=_m2.begin();
       it!=_m2.end();++it){
    dummy = (*it)->GetParametrizedMess();
    x += dummy * (*it)->GetDirection()[0];
    y += dummy * (*it)->GetDirection()[1];  
  }
  return x*x+y*y;     
};

double TData_InvMass::combine(){
  double x, y, z=0., e, tx, ty, tz, dummy = GetParametrizedMess();
  //here _direction[1,2] is a normalized vector in eta-phi plane in pT direction.
  //_direction[3] is the p_z component of the measurement.
  
  x = dummy * _direction[0];
  y = dummy * _direction[1];
  if (_mess[0]!=0.) z = dummy/_mess[0] * _direction[2];
  e = sqrt( x*x + y*y + z*z );
  for (std::vector<TData_MessMess*>::const_iterator it=_m2.begin();
       it!=_m2.end();++it){
    dummy = (*it)->GetParametrizedMess();
    tx = dummy * (*it)->GetDirection()[0];
    ty = dummy * (*it)->GetDirection()[1];  
    if ((*it)->GetMess()[0]!=0.) tz = dummy/(*it)->GetMess()[0] * (*it)->GetDirection()[2];  
    else tz = 0.;
    x += tx;
    y += ty;
    z += tz;
    e += sqrt( tx*tx + ty*ty + tz*tz );
  }
  return e*e - x*x - y*y - z*z;     
};
