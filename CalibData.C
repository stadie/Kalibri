#include "CalibData.h"
#include "map"
//#include <iostream>//cout

unsigned int TData::total_n_pars = 0;
//unsigned int TData::total_tower_pars = 0;
//unsigned int TData::total_jet_pars = 0;

double * TData::temp_derivative1 = 0;
double * TData::temp_derivative2 = 0;
double const TData::epsilon = 1.E-3;

std::vector<TData*> TData_TruthMess::resultcache = std::vector<TData*>(1);

double TData_TruthMess::chi2_fast(){ 
  double new_mess  = GetParametrizedMess();
  double new_error = GetParametrizedErr(&new_mess);
  double new_chi2  = (_truth-new_mess)*(_truth-new_mess)/(new_error*new_error);
  
  double dmess_dp, derror_dp;
  unsigned idx = _index*_n_par; //_index==bin; idx==bin*Free_parameters_per_bin
  for (unsigned i=0; i<idx; ++i){
    temp_derivative2[i]+=new_chi2;
    temp_derivative1[i]+=new_chi2;
  }
  for (unsigned i=idx; i<idx+_n_par; ++i){
    _par[i-idx]  += epsilon;
    dmess_dp  = GetParametrizedMess();
    derror_dp = GetParametrizedErr(&dmess_dp);
    temp_derivative2[i]+= (_truth-dmess_dp)*(_truth-dmess_dp)/(derror_dp*derror_dp);

    _par[i-idx]  -= 2.0*epsilon;
    dmess_dp  = GetParametrizedMess();
    derror_dp = GetParametrizedErr(&dmess_dp);
    temp_derivative1[i]+= (_truth-dmess_dp)*(_truth-dmess_dp)/(derror_dp*derror_dp);
    _par[i-idx]  += epsilon;
  }
  for (unsigned i=idx+_n_par; i< total_n_pars; ++i){
    temp_derivative2[i]+=new_chi2;
    temp_derivative1[i]+=new_chi2;
  }
  return new_chi2;
};

double TData_TruthMultMess::chi2_fast(){
  double sum_mess=0.0, sum_error2=0.0, new_error, new_error2, new_mess, new_chi2,
         dmess_dp, derror_dp;
  double sm1[total_n_pars],sm2[total_n_pars],se1[total_n_pars],se2[total_n_pars];//store sum_m & sum_e2 for df/dp_i
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
    new_error   = (*it)->GetParametrizedErr(&new_mess);
    new_error2  = new_error * new_error;
    sum_error2 += new_error2;
    for (unsigned i=0; i<total_n_pars; ++i){
      idx = (*it)->GetIndex()*(*it)->GetNumberOfPars();
      if (i>=idx && i<idx+(*it)->GetNumberOfPars()) {
	(*it)->AddToPar(i-idx,epsilon);
	dmess_dp  =(*it)->GetParametrizedMess();
	sm2[i]   += dmess_dp;
	new_error = (*it)->GetParametrizedErr(&dmess_dp);
	se2[i]   += new_error * new_error;
	
	(*it)->AddToPar(i-idx,-2.0*epsilon);
	dmess_dp  =(*it)->GetParametrizedMess();
	sm1[i]   += dmess_dp;
	new_error = (*it)->GetParametrizedErr(&dmess_dp);
	se1[i]   += new_error * new_error;
	(*it)->AddToPar(i-idx,epsilon);
      } else {
	sm1[i] += new_mess;
	sm2[i] += new_mess;
	se1[i] += new_error2;
	se2[i] += new_error2;
      }
    }
  } 
  new_mess  = _func( &sum_mess, _par);  
  new_error = _err(&new_mess); 
  new_chi2  = (_truth-new_mess)*(_truth-new_mess)/(sum_error2 + new_error*new_error);

  idx = _index; //@@to be fixed -> introduce a eta-phi binning for JES
  for (unsigned i=0; i<total_n_pars; ++i){
    if (i>=idx && i<idx+_n_par) continue;//considered below
    new_mess  = sm2[i];  //the measurement with modified parameter "i"
    dmess_dp  = _func(&new_mess,_par);//calc. the jet's energy
    //derror_dp = _error*dmess_dp/new_mess;
    derror_dp = _err(&dmess_dp);
    temp_derivative2[i]+= (_truth-dmess_dp)*(_truth-dmess_dp)/(se2[i] + derror_dp*derror_dp);
    // same for p_i-epsilon:
    new_mess  = sm1[i];  
    dmess_dp  = _func(&new_mess,_par);
    //derror_dp = _error*dmess_dp/new_mess;
    derror_dp = _err(&dmess_dp);
    temp_derivative1[i]+= (_truth-dmess_dp)*(_truth-dmess_dp)/(se1[i] + derror_dp*derror_dp);
  }
  	
  for (unsigned i=idx; i<idx+_n_par; ++i){
    //ok, we have to change the jet's parametrization:
    _par[i-idx]  += epsilon;
    dmess_dp  = _func(&sum_mess,_par);
    //derror_dp = _error*dmess_dp/sum_mess;
    derror_dp = _err(&dmess_dp);
    temp_derivative2[i]+= (_truth-dmess_dp)*(_truth-dmess_dp)/(se2[i] + derror_dp*derror_dp);

    _par[i-idx]  -= 2.0*epsilon;
    dmess_dp  = _func(&sum_mess,_par);
    //derror_dp = _error*dmess_dp/sum_mess;
    derror_dp = _err(&dmess_dp);
    temp_derivative1[i]+= (_truth-dmess_dp)*(_truth-dmess_dp)/(se1[i] + derror_dp*derror_dp);
    _par[i-idx]  += epsilon;
  }

  return new_chi2;
};

#include <iostream>

double TData_PtBalance::chi2_fast()
{
  double new_mess, new_error, sum_error2=0.0;

  //Get all tower parameter used in this event:  
  std::map<int,double*> tpars;
  for (std::vector<TData*>::const_iterator it=_vecmess.begin();
       it!=_vecmess.end(); ++it) 
    tpars[ (*it)->GetIndex() ] = (*it)->GetPar();   
  for (std::vector<TData_MessMess*>::const_iterator mit=_m2.begin();
       mit!=_m2.end();++mit){
     std::vector<TData*>::const_iterator mitend=(*mit)->GetRef().end();  
     for (std::vector<TData*>::const_iterator it=(*mit)->GetRef().begin();
	it!=mitend; ++it) {
       tpars[ (*it)->GetIndex() ] = (*it)->GetPar();   
       }
     }  

  //Calc. & Cache derivatives w.r.t. tower parameters:
  for (std::map<int,double*>::const_iterator it=tpars.begin(); 
       it!=tpars.end(); ++it ) {
     for (unsigned i=0; i< (*_vecmess.begin())->GetNumberOfPars(); ++i){
       (it->second)[i] += epsilon;

       for (std::vector<TData_MessMess*>::const_iterator mit=_m2.begin();
            mit!=_m2.end(); ++mit) {
	  new_mess    = (*mit)->GetParametrizedMess();	 
  	  new_error   = (*mit)->GetParametrizedErr(&new_mess);
	  sum_error2 += new_error * new_error;
       }
       new_mess  = GetParametrizedMess();
       new_error = _err(  &new_mess );
       new_mess  = GetMessCombination();  
       temp_derivative2[ it->first+i ]+=(_truth-new_mess)*(_truth-new_mess)/(sum_error2 + new_error*new_error);

       (it->second)[i] -= 2.*epsilon;
       for (std::vector<TData_MessMess*>::const_iterator mit=_m2.begin();
            mit!=_m2.end(); ++mit) {
	  new_mess    = (*mit)->GetParametrizedMess();	 
  	  new_error   = (*mit)->GetParametrizedErr(&new_mess);
	  sum_error2 += new_error * new_error;
       }
       new_mess  = GetParametrizedMess();
       new_error = _err(  &new_mess );
       new_mess  = GetMessCombination();  
       temp_derivative1[ it->first+i ]+=(_truth-new_mess)*(_truth-new_mess)/(sum_error2 + new_error*new_error);
       (it->second)[i] += epsilon;
//std::cout << "par["<<(it->first)<<"+"<<i<<"] = " << (it->second)[i]<< std::endl;
    }
  }

  //Get all jet parameter used in this event:  
  tpars.clear();
  tpars[_index]=_par;
  for (std::vector<TData_MessMess*>::const_iterator mit=_m2.begin();
       mit!=_m2.end();++mit)
    tpars[ (*mit)->GetIndex() ] = (*mit)->GetPar();
  //Calc. & Cache derivatives w.r.t. jet parameters:
  for (std::map<int,double*>::const_iterator it=tpars.begin(); 
       it!=tpars.end(); ++it ) {
     for (unsigned i=0; i<_n_par; ++i){
       (it->second)[i] += epsilon;

       for (std::vector<TData_MessMess*>::const_iterator mit=_m2.begin();
            mit!=_m2.end(); ++mit) {
	  new_mess    = (*mit)->GetParametrizedMess();	 
  	  new_error   = (*mit)->GetParametrizedErr(&new_mess);
	  sum_error2 += new_error * new_error;
       }
       new_mess  = GetParametrizedMess();
       new_error = _err(  &new_mess );
       new_mess  = GetMessCombination();  
       temp_derivative2[ it->first+i ]+=(_truth-new_mess)*(_truth-new_mess)/(sum_error2 + new_error*new_error);

       (it->second)[i] -= 2.*epsilon;
       for (std::vector<TData_MessMess*>::const_iterator mit=_m2.begin();
            mit!=_m2.end(); ++mit) {
	  new_mess    = (*mit)->GetParametrizedMess();	 
  	  new_error   = (*mit)->GetParametrizedErr(&new_mess);
	  sum_error2 += new_error * new_error;
       }
       new_mess  = GetParametrizedMess();
       new_error = _err(  &new_mess );
       new_mess  = GetMessCombination();  
       temp_derivative1[ it->first+i ]+=(_truth-new_mess)*(_truth-new_mess)/(sum_error2 + new_error*new_error);
       (it->second)[i] += epsilon;
//std::cout << "par["<<(it->first)<<"+"<<i<<"] = " << (it->second)[i]<< std::endl;
    }
  }

  //std::cout<< "TData_PtBalance::chi2_fast tpars->size= " << tpars.size() << std::endl;
  return chi2();
}

/*
double TData_PtBalance::chi2_fast()
{
  double sum_mess=0.0, sum_error2=0.0, new_error, new_error2, 
         new_mess, new_chi2, dmess_dp, derror_dp;
  double sm1[total_n_pars],sm2[total_n_pars],se1[total_n_pars],se2[total_n_pars];//store sum_m & sum_e2 for df/dp_i
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
    new_error   = (*it)->GetParametrizedErr(&new_mess);
    new_error2  = new_error * new_error;
    sum_error2 += new_error2;
    for (unsigned i=0; i<total_n_pars; ++i){
      idx = (*it)->GetIndex()*(*it)->GetNumberOfPars();
      if (i>=idx && i<idx+(*it)->GetNumberOfPars()) {
	(*it)->AddToPar(i-idx,epsilon);
	dmess_dp  =(*it)->GetParametrizedMess();
	sm2[i]   += dmess_dp;
	new_error = (*it)->GetParametrizedErr(&dmess_dp);
	se2[i]   += new_error * new_error;
	
	(*it)->AddToPar(i-idx,-2.0*epsilon);
	dmess_dp  =(*it)->GetParametrizedMess();
	sm1[i]   += dmess_dp;
	new_error = (*it)->GetParametrizedErr(&dmess_dp);
	se1[i]   += new_error * new_error;
	(*it)->AddToPar(i-idx,epsilon);
      } else {
	sm1[i] += new_mess;
	sm2[i] += new_mess;
	se1[i] += new_error2;
	se2[i] += new_error2;
      }
    }
  } 
  new_mess  = _func( &sum_mess, _par);  
  new_error = _err(&new_mess); 
  new_mess  = GetMessCombination();  
  new_chi2  = (_truth-new_mess)*(_truth-new_mess)/(sum_error2 + new_error*new_error);

  //derivative calc. w.r.t. tower-parameters using param.mess from above
  idx = _index; //@@to be fixed -> introduce a eta-phi binning for JES
  for (unsigned i=0; i<total_n_pars; ++i){
    if (i>=idx && i<idx+_n_par) continue;//considered below
    new_mess  = sm2[i];  //the measurement with modified parameter "i"
    dmess_dp  = _func(&new_mess,_par);//calc. the jet's energy
    //derror_dp = _error*dmess_dp/new_mess;
    derror_dp = _err(&dmess_dp);
    dmess_dp  = GetMessCombination();  
    temp_derivative2[i]+= (_truth-dmess_dp)*(_truth-dmess_dp)/(se2[i] + derror_dp*derror_dp);
    // same for p_i-epsilon:
    new_mess  = sm1[i];  
    dmess_dp  = _func(&new_mess,_par);
    //derror_dp = _error*dmess_dp/new_mess;
    derror_dp = _err(&dmess_dp);
    dmess_dp  = GetMessCombination();  
    temp_derivative1[i]+= (_truth-dmess_dp)*(_truth-dmess_dp)/(se1[i] + derror_dp*derror_dp);
  }

  //derivatives w.r.t. to the jet's parametrization:
  for (unsigned i=idx; i<idx+_n_par; ++i){
    _par[i-idx]  += epsilon;
    dmess_dp  = _func(&sum_mess,_par);
    //derror_dp = _error*dmess_dp/sum_mess;
    derror_dp = _err(&dmess_dp);
    dmess_dp  = GetMessCombination();  
    temp_derivative2[i]+= (_truth-dmess_dp)*(_truth-dmess_dp)/(se2[i] + derror_dp*derror_dp);

    _par[i-idx]  -= 2.0*epsilon;
    dmess_dp  = _func(&sum_mess,_par);
    //derror_dp = _error*dmess_dp/sum_mess;
    derror_dp = _err(&dmess_dp);
    dmess_dp  = GetMessCombination();  
    temp_derivative1[i]+= (_truth-dmess_dp)*(_truth-dmess_dp)/(se1[i] + derror_dp*derror_dp);
    _par[i-idx]  += epsilon;
  }

  return new_chi2;
};
*/

double TData_PtBalance::combine(){
  double x, y, dummy = GetParametrizedMess();
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
