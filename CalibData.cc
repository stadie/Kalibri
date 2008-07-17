//
// Original Author:  Christian Autermann
//         Created:  Wed Jul 18 13:54:50 CEST 2007
// $Id: CalibData.cc,v 1.9 2008/07/16 14:36:04 thomsen Exp $
//
#include "CalibData.h"
#include "map"
//#include <iostream>//cout

unsigned int TData::total_n_pars = 0;
double (*TData::ScaleResidual)(double z2) = &TData::ScaleNone;

//unsigned int TData::total_tower_pars = 0;
//unsigned int TData::total_jet_pars = 0;

std::vector<TData*> TData_TruthMess::resultcache = std::vector<TData*>(1);

double TData_TruthMess::chi2_fast(double *temp_derivative1, double *temp_derivative2, double epsilon)
{ 
  double new_mess  = GetParametrizedMess();
#ifndef __FastErrorCalculation
  double new_error = GetParametrizedErr(&new_mess);
  
#else
  double new_error = _error;
  if (GetMess()[0]!=0.) new_error = _error*new_mess/GetMess()[0];
#endif   
  double weight = GetWeight();
  double new_chi2 = weight*(*TData::ScaleResidual)( (_truth-new_mess)*(_truth-new_mess)/(new_error*new_error) );

  double temp1 = 0.;		// Value of chi2 at par+epsilon
  double temp2 = 0.;		// Value of chi2 at par-epsilon
  
  double dmess_dp, derror_dp;   
  unsigned idx = _index*_n_par; //_index==bin; idx==bin*Free_parameters_per_bin
  for (unsigned i=idx; i<idx+_n_par; ++i){
    temp1 = 0.;
    temp2 = 0.;
    double oldpar = _par[i-idx];
    _par[i-idx]  += epsilon;
    dmess_dp  = GetParametrizedMess();
#ifndef __FastErrorCalculation
    derror_dp = GetParametrizedErr(&dmess_dp);
#else
    if (GetMess()[0]!=0.) derror_dp = _error*dmess_dp/GetMess()[0];
    else derror_dp = _error;
#endif   
    temp2 = weight*(*TData::ScaleResidual)( (_truth-dmess_dp)*(_truth-dmess_dp)/(derror_dp*derror_dp) );


    _par[i-idx]  = oldpar - epsilon;
    dmess_dp  = GetParametrizedMess();
#ifndef __FastErrorCalculation
    derror_dp = GetParametrizedErr(&dmess_dp);
#else
    if (GetMess()[0]!=0.) derror_dp = _error*dmess_dp/GetMess()[0];
    else derror_dp = _error;
#endif   
    temp1 = weight*(*TData::ScaleResidual)( (_truth-dmess_dp)*(_truth-dmess_dp)/(derror_dp*derror_dp) );

    // Difference of chi2 at par+epsilon and par-epsilon
    temp_derivative1[i] += (temp2 - temp1); // for 1st derivative
    temp_derivative2[i] += (temp2 + temp1 - 2*new_chi2); // for 2nd derivative

    _par[i-idx]  = oldpar;
  }
  return new_chi2;
};



double TData_TruthMultMess::chi2_fast(double* temp_derivative1, double* temp_derivative2, double epsilon) {
  // store sum_m & sum_e2 
  TLorentzVector sum_vec(0,0,0,0);
  double sum_mess = 0.;
  double sum_error2 = 0.;
  TLorentzVector single_vec(0,0,0,0);
  double single_mess = 0.;
  double single_error = 0.;
  double single_error2 = 0.;
  TLorentzVector dvec_dp(0,0,0,0);
  double new_mess = 0.;
  double new_error = 0.;
  double dmess_dp = 0.;
  double derror_dp = 0.;   
  double new_chi2 = 0;

  TLorentzVector sm1[total_n_pars];
  TLorentzVector sm2[total_n_pars];
  double se1[total_n_pars];
  double se2[total_n_pars];

  double weight = GetWeight();

  for (unsigned i=0; i<total_n_pars; ++i){
    sm1[i] = TLorentzVector(0,0,0,0);
    sm2[i] = TLorentzVector(0,0,0,0);
    se1[i] = 0.0;
    se2[i] = 0.0;
  }
  unsigned idx;

  for (std::vector<TData*>::const_iterator it=_vecmess.begin();
       it!=_vecmess.end(); ++it) {
    single_vec.SetPtEtaPhiM((*it)->GetParametrizedMess(),(*it)->GetMess()[4],(*it)->GetMess()[5],0);
    sum_vec += single_vec;
//#ifndef __FastErrorCalculation
    single_mess = single_vec.Et();
    single_error = (*it)->GetParametrizedErr(&single_mess);

//#else
//    if ((*it)->GetMess()[0]!=0.) new_error = (*it)->GetError()*new_mess/(*it)->GetMess()[0];
//    else new_error = (*it)->GetError();
//#endif   
    single_error2 = single_error * single_error;
    sum_error2 += single_error2;
    for (unsigned i=0; i<total_n_pars; ++i) {
      idx = (*it)->GetIndex()*(*it)->GetNumberOfPars(); // Sollte aus for-Schleife raus?
      if (i>=idx && i<idx+(*it)->GetNumberOfPars()) {   
	double oldpar = (*it)->GetPar()[i-idx];
	(*it)->GetPar()[i-idx]  += epsilon;
	dvec_dp.SetPtEtaPhiM((*it)->GetParametrizedMess(),(*it)->GetMess()[4],(*it)->GetMess()[5],0);
	sm2[i] += dvec_dp;
//#ifndef __FastErrorCalculation
	single_mess = dvec_dp.Et();
	single_error = (*it)->GetParametrizedErr(&single_mess);

//#else
//        if ((*it)->GetMess()[0]!=0.) new_error = (*it)->GetError()*dmess_dp/(*it)->GetMess()[0];
//	else new_error = (*it)->GetError();
//#endif   
	se2[i] += single_error * single_error;
	
	(*it)->GetPar()[i-idx]  = oldpar - epsilon;
	dvec_dp.SetPtEtaPhiM((*it)->GetParametrizedMess(),(*it)->GetMess()[4],(*it)->GetMess()[5],0);
	sm1[i] += dvec_dp;
//#ifndef __FastErrorCalculation
	single_mess = dvec_dp.Et();
	single_error = (*it)->GetParametrizedErr(&single_mess);

//#else
//        if ((*it)->GetMess()[0]!=0.) new_error = (*it)->GetError()*dmess_dp/(*it)->GetMess()[0];
//	else new_error = (*it)->GetError();
//#endif   
	se1[i] += single_error * single_error;
	(*it)->GetPar()[i-idx] = oldpar;
      } else {
	sm1[i] += single_vec;
	sm2[i] += single_vec;
	se1[i] += single_error2;
	se2[i] += single_error2;
      }
    }
  } 
  sum_mess = sum_vec.Et();
  new_mess  = GetJetCor( &sum_mess); 
  //new_mess  = GetParametrizedMess();               //jan
//#ifndef __FastErrorCalculation
  new_error = _err(&new_mess); 
//#else
//  if (_mess!=0) if (_mess[0]!=0.) new_error = _error*new_mess/_mess[0];
//  else new_error = _error;
//#endif   
  new_chi2  = weight*(*TData::ScaleResidual)( (_truth-new_mess)*(_truth-new_mess)/(sum_error2 + new_error*new_error) );

  double temp1 = 0.;		// Value of chi2 at par+epsilon
  double temp2 = 0.;		// Value of chi2 at par-epsilon
  idx = _index; //@@to be fixed -> introduce a eta-phi binning for JES
  for (unsigned i=0; i<total_n_pars; ++i){
    if (i>=idx && i<idx+_n_par) continue;//considered below
    temp1 = 0.;
    temp2 = 0.;
    new_mess  = sm2[i].Et();  //the measurement with modified parameter "i"
    dmess_dp  = GetJetCor(&new_mess);//calc. the jet's energy 
    //dmess_dp  = GetParametrizedMess();                //jan
//#ifndef __FastErrorCalculation
    derror_dp = _err(&dmess_dp);
//#else
//    if (new_mess!=0.) derror_dp = _error*dmess_dp/new_mess;
//    else derror_dp = _error;
//#endif   
    temp2 = weight*(*TData::ScaleResidual)( (_truth-dmess_dp)*(_truth-dmess_dp)/(se2[i] + derror_dp*derror_dp) );

    // same for p_i-epsilon:
    new_mess  = sm1[i].Et();  
    dmess_dp  = GetJetCor(&new_mess); 
    //dmess_dp  = GetParametrizedMess();              //jan
//#ifndef __FastErrorCalculation
    derror_dp = _err(&dmess_dp);
//#else
//    if (new_mess!=0.) derror_dp = _error*dmess_dp/new_mess;
//    else derror_dp = _error;
//#endif   
    temp1 = weight*(*TData::ScaleResidual)( (_truth-dmess_dp)*(_truth-dmess_dp)/(se1[i] + derror_dp*derror_dp) );

    // Difference of chi2 at par+epsilon and par-epsilon
    temp_derivative1[i] += (temp2 - temp1); // for 1st derivative
    temp_derivative2[i] += (temp2 + temp1 - 2*new_chi2); // for 2nd derivative
  }
  	
  for (unsigned i=idx; i<idx+_n_par; ++i){
    temp1 = 0.;
    temp2 = 0.;
    //ok, we have to change the jet's parametrization:
    double oldpar =  _par[i-idx];
    _par[i-idx]  += epsilon;
    new_mess = sum_vec.Et();
    dmess_dp  = GetJetCor(&new_mess); 
    //dmess_dp  = GetParametrizedMess();               //jan
//#ifndef __FastErrorCalculation
    derror_dp = _err(&dmess_dp);
//#else
//    if (sum_mess!=0.) derror_dp = _error*dmess_dp/sum_mess;
//    else derror_dp = _error;
//#endif
    temp2 = weight*(*TData::ScaleResidual)( (_truth-dmess_dp)*(_truth-dmess_dp)/(se2[i] + derror_dp*derror_dp) );

    _par[i-idx]  = oldpar - epsilon;
    new_mess = sum_vec.Et();
    dmess_dp  = GetJetCor(&new_mess);              //jan
    //dmess_dp  = GetParametrizedMess();  
//#ifndef __FastErrorCalculation
    derror_dp = _err(&dmess_dp);
//#else
//    if (sum_mess!=0.) derror_dp = _error*dmess_dp/sum_mess;
//    else derror_dp = _error;
//#endif
    temp1 = weight*(*TData::ScaleResidual)( (_truth-dmess_dp)*(_truth-dmess_dp)/(se1[i] + derror_dp*derror_dp) );
    // Difference of chi2 at par+epsilon and par-epsilon
    temp_derivative1[i] += (temp2 - temp1); // for 1st derivative
    temp_derivative2[i] += (temp2 + temp1 - 2*new_chi2); // for 2nd derivative

    _par[i-idx]  = oldpar;
  }

  return new_chi2;
};


double TData_MessMess::chi2_fast(double * temp_derivative1, double*  temp_derivative2, 
			        double epsilon)
{
  double new_chi2, new_mess,  sum_error2=0.0;
  double weight = GetWeight();

  //Get all tower parameter used in this event:  
  std::map<int,double*> tpars;
  for (std::vector<TData*>::const_iterator it=_vecmess.begin();
       it!=_vecmess.end(); ++it) {
    for (unsigned i= 0 ; i < (*it)->GetNumberOfPars(); ++i) {
      tpars[ (*it)->GetIndex() * (*it)->GetNumberOfPars() + i ] = &((*it)->GetPar()[i]);
    }
  }
  for (std::vector<TData_MessMess*>::const_iterator mit=_m2.begin();
       mit!=_m2.end();++mit) {
    std::vector<TData*>::const_iterator mitend=(*mit)->GetRef().end();  
    for (std::vector<TData*>::const_iterator it=(*mit)->GetRef().begin();
	 it!=mitend; ++it) { 
      for (unsigned i= 0 ; i < (*it)->GetNumberOfPars(); ++i) {
	tpars[ (*it)->GetIndex() * (*it)->GetNumberOfPars() + i ] = &((*it)->GetPar()[i]);
      }
    }  
  }
  //Add all jet parameters in this event:
  for (unsigned i=0; i<_n_par; ++i)
    tpars[ i+_index ]= &(_par[i]);
  for (std::vector<TData_MessMess*>::const_iterator mit=_m2.begin();
       mit!=_m2.end();++mit)
    for (unsigned i=0; i<(*mit)->GetNumberOfPars(); ++i) {
      tpars[ i+(*mit)->GetIndex() ]= &((*mit)->GetPar()[i]);
    }
  
  //This event's chi^2 for the current (unchanged) parameters:
  new_chi2 = chi2();

  //Calc. & Cache derivatives w.r.t. tower parameters:
  double temp1 = 0.;		//Value of chi2 at par+epsilon
  double temp2 = 0.;		//Value of chi2 at par-epsilon
  for (unsigned i=0; i<total_n_pars; ++i){
     // in case the ith parameter is used for this event, calculate 
     // derivative dchi/dpar[i] and cache it:
     if ( tpars.find(i)!=tpars.end() ) {
       double oldpar =  *tpars[i];
       *tpars[i] += epsilon;
       sum_error2 = 0.0;
       
       for (std::vector<TData_MessMess*>::const_iterator mit=_m2.begin();
            mit!=_m2.end(); ++mit) {
	 //#ifndef __FastErrorCalculation
	 sum_error2 += (*mit)->GetParametrizedErr2(); //all tower & jet errors
	 //#else
	 //  	  if ((*mit)->GetMess()) if ((*mit)->GetMess()[0]!=0.) new_error = (*mit)->GetError()*new_mess/(*mit)->GetMess()[0];
//	  else new_error   = (*mit)->GetError();
//#endif
       }
       sum_error2 += GetParametrizedErr2(); //total error first jet
       new_mess  = GetMessCombination();  
       temp2 = weight*(*TData::ScaleResidual)( (_truth-new_mess)*(_truth-new_mess)/sum_error2 );
       sum_error2 = 0.0;
       

       *tpars[i] = oldpar - epsilon;

       for (std::vector<TData_MessMess*>::const_iterator mit=_m2.begin();
            mit!=_m2.end(); ++mit) {
	 //#ifndef __FastErrorCalculation
	 sum_error2 += (*mit)->GetParametrizedErr2(); //all tower & jet errors
	 //#else
	 //  	  if ((*mit)->GetMess()) if ((*mit)->GetMess()[0]!=0.) new_error = (*mit)->GetError()*new_mess/(*mit)->GetMess()[0];
//	  else new_error   = (*mit)->GetError();
//#endif
       }
       sum_error2 += GetParametrizedErr2(); //total error first jet
       new_mess  = GetMessCombination();  
       temp1 = weight*(*TData::ScaleResidual)( (_truth-new_mess)*(_truth-new_mess)/sum_error2 );
      // Difference of chi2 at par+epsilon and par-epsilon
      temp_derivative1[i] += (temp2 - temp1); // for 1st derivative
      temp_derivative2[i] += (temp2 + temp1 - 2*new_chi2); // for 2nd derivative
      
      *tpars[i] = oldpar;
     }
  }
  return new_chi2;
}


double TData_PtBalance::combine(){
  double x, y, dummy = GetParametrizedMess();
  //here _direction is a normalized vector in eta-phi plane in pT direction.
  
  x = dummy * _direction[0];
  y = dummy * _direction[1];
  //x = 0;
  //y = 0;
  //dummy = 0;
  for (std::vector<TData_MessMess*>::const_iterator it=_m2.begin();
       it!=_m2.end();++it){
    dummy = (*it)->GetParametrizedMess(); 
    x += dummy * (*it)->GetDirection()[0];
    y += dummy * (*it)->GetDirection()[1];  
  }
  return sqrt(x*x+y*y);
  //return GetParametrizedMess() - dummy;
  //return GetMess()[0] - dummy;
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
  return sqrt(e*e - x*x - y*y - z*z);     
};

std::vector<TData*> TData_ParLimit::_cache = std::vector<TData*>(1);

double TData_ParLimit::chi2_fast(double * temp_derivative1, 
				 double* temp_derivative2, double epsilon) {
  double new_chi2 = chi2();
  //if((_par[0] < _mess[0]) || (_par[0] > _mess[1])) {
  //  std::cout << "WARNING: paramter " << _index << " at limit:" 
  //	      << _par[0] << '\n'; 
  //}
  //std::cout << "ParLimit: par=" << _par[0] << " chi2:" <<   new_chi2 << '\n';
  double oldpar = _par[0];
  _par[0] += epsilon;
  double temp2 = chi2();
  _par[0]  = oldpar - epsilon;
  double temp1 = chi2();
  // Difference of chi2 at par+epsilon and par-epsilon
  temp_derivative1[_index] += (temp2 - temp1); // for 1st derivative
  temp_derivative2[_index] += (temp2 + temp1 - 2*new_chi2); // for 2nd derivative
  _par[0]  = oldpar;
  return new_chi2;
}
    

/*
  Scale the normalized residual 'z^2 = chi^2/weight' using 
  the Cauchy-Function:  z^2 --> c^2 * Ln( 1 + (z/c)^2 )
*/
double TData::ScaleCauchy(double z2)
{
  double c = 2.3849;
  return c*c * log( 1 + z2/(c*c) );
}

/*
  Scale the normalized residual 'z^2 = chi^2/weight' using 
  the Huber-Function   z^2 --> z                 for |z| <= c
	                       c * ( 2*|z| - c ) for |z| > c
*/
double TData::ScaleHuber(double z2)
{
  double c = 1.345;
  double z = sqrt(z2);
  return (  fabs(z) <= c  ?  z2  :  c*(2*fabs(z) - c)  );
}
