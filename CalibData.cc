//
// Original Author:  Christian Autermann
//         Created:  Wed Jul 18 13:54:50 CEST 2007
// $Id: CalibData.cc,v 1.13 2008/07/30 15:19:38 auterman Exp $
//
#include "CalibData.h"
#include "map"
//#include <iostream>//cout

unsigned int TData::total_n_pars = 0;
double (*TData::ScaleResidual)(double z2) = &TData::ScaleNone;
std::vector<TData*> TData_TruthMess::resultcache = std::vector<TData*>(1);

double TData_TruthMess::chi2_fast(double *temp_derivative1, double *temp_derivative2, double epsilon)
{ 
  double new_mess  = GetParametrizedMess();
  double new_error = GetParametrizedErr(&new_mess);
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
    derror_dp = GetParametrizedErr(&dmess_dp);
    temp2 = weight*(*TData::ScaleResidual)( (_truth-dmess_dp)*(_truth-dmess_dp)/(derror_dp*derror_dp) );


    _par[i-idx]  = oldpar - epsilon;
    dmess_dp  = GetParametrizedMess();
    derror_dp = GetParametrizedErr(&dmess_dp);
    temp1 = weight*(*TData::ScaleResidual)( (_truth-dmess_dp)*(_truth-dmess_dp)/(derror_dp*derror_dp) );

    // Difference of chi2 at par+epsilon and par-epsilon
    temp_derivative1[i] += (temp2 - temp1); // for 1st derivative
    temp_derivative2[i] += (temp2 + temp1 - 2*new_chi2); // for 2nd derivative

    _par[i-idx]  = oldpar;
  }
  return new_chi2;
};



double TData_TruthMultMess::chi2_fast(double* temp_derivative1, double* temp_derivative2, double epsilon) {
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
    //@@ Why not replacing "tower pT" by "tower pT * (*it)->GetMess()[7]"
    //   directly when TData vector is filled ???
    sum_mess   += new_mess; // * (*it)->GetMess()[7]; // Sum of tower Pt
    new_error   = (*it)->GetParametrizedErr(&new_mess);
    new_error2  = new_error * new_error;
    sum_error2 +=  new_error2;
    for (unsigned i=0; i<total_n_pars; ++i){
      idx = (*it)->GetIndex()*(*it)->GetNumberOfPars(); // Sollte aus for-Schleife raus?
      if (i>=idx && i<idx+(*it)->GetNumberOfPars()) {   
	double oldpar = (*it)->GetPar()[i-idx];
	(*it)->GetPar()[i-idx]  += epsilon;
	dmess_dp  =(*it)->GetParametrizedMess();
	sm2[i]   += dmess_dp;// * (*it)->GetMess()[7]; // Sum of tower Pt
	new_error = (*it)->GetParametrizedErr(&dmess_dp);
	se2[i]   += new_error * new_error;
	
	(*it)->GetPar()[i-idx]  = oldpar - epsilon;
	dmess_dp  =(*it)->GetParametrizedMess();
	sm1[i]   += dmess_dp;// * (*it)->GetMess()[7]; // Sum of tower Pt
	new_error = (*it)->GetParametrizedErr(&dmess_dp);
	se1[i]   += new_error * new_error;
	(*it)->GetPar()[i-idx] = oldpar;
      } else {
	sm1[i] += new_mess;// * (*it)->GetMess()[7]; // Sum of tower Pt
	sm2[i] += new_mess;// * (*it)->GetMess()[7]; // Sum of tower Pt
	se1[i] += new_error2;
	se2[i] += new_error2;
      }
    }
  } 
  new_mess  = GetParametrizedMess();
  new_error = GetParametrizedErr(&new_mess); 
  new_chi2  = weight*(*TData::ScaleResidual)( (_truth-new_mess)*(_truth-new_mess)/(sum_error2 + new_error*new_error) );

  double temp1 = 0.;		// Value of chi2 at par+epsilon
  double temp2 = 0.;		// Value of chi2 at par-epsilon
  idx = _index; //@@to be fixed -> introduce a eta-phi binning for JES
  for (unsigned i=0; i<total_n_pars; ++i){
    if (i>=idx && i<idx+_n_par) continue;//considered below
    temp1 = 0.;
    temp2 = 0.;
    new_mess  = sm2[i];  //the measurement with modified parameter "i"

    dmess_dp  = GetParametrizedMess();
    derror_dp = GetParametrizedErr(&dmess_dp);
    temp2 = weight*(*TData::ScaleResidual)( (_truth-dmess_dp)*(_truth-dmess_dp)/(se2[i] + derror_dp*derror_dp) );

    // same for p_i-epsilon:
    new_mess  = sm1[i];  
    dmess_dp  = GetParametrizedMess();
    derror_dp = GetParametrizedErr(&dmess_dp);
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
    dmess_dp  = GetParametrizedMess();
    derror_dp = GetParametrizedErr(&dmess_dp);
    temp2 = weight*(*TData::ScaleResidual)( (_truth-dmess_dp)*(_truth-dmess_dp)/(se2[i] + derror_dp*derror_dp) );

    _par[i-idx]  = oldpar - epsilon;
    dmess_dp  = GetParametrizedMess();  
    derror_dp = GetParametrizedErr(&dmess_dp);
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
  double new_chi2, new_mess, new_error, sum_error2=0.0;
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
	  new_mess    = (*mit)->GetParametrizedMess();	 
  	  new_error   = (*mit)->GetParametrizedErr(&new_mess);
	  sum_error2 += new_error * new_error;
       }
       new_mess  = GetParametrizedMess();
       new_error = GetParametrizedErr( &new_mess );
       new_mess  = GetMessCombination();  
       temp2 += weight*(_truth-new_mess)*(_truth-new_mess)/(sum_error2 + new_error*new_error);
       sum_error2 = 0.0;

       *tpars[i] = oldpar - epsilon;
       for (std::vector<TData_MessMess*>::const_iterator mit=_m2.begin();
            mit!=_m2.end(); ++mit) {
	  new_mess    = (*mit)->GetParametrizedMess();	 
  	  new_error   = (*mit)->GetParametrizedErr(&new_mess);
	  sum_error2 += new_error * new_error;
       }
       new_mess  = GetParametrizedMess();
       new_error = GetParametrizedErr(  &new_mess );
       new_mess  = GetMessCombination();  
       temp1 += weight*(_truth-new_mess)*(_truth-new_mess)/(sum_error2 + new_error*new_error);
   
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
  //here _direction[0,1] is a normalized vector in eta-phi plane in pT direction.
  
  x = dummy * _direction[0];
  y = dummy * _direction[1];
  for (std::vector<TData_MessMess*>::const_iterator it=_m2.begin();
       it!=_m2.end();++it){
    dummy = (*it)->GetParametrizedMess(); 
    x += dummy * (*it)->GetDirection()[0];
    y += dummy * (*it)->GetDirection()[1];  
  }
  return sqrt(x*x+y*y);
};

double TData_InvMass2::combine(){
  double x, y, z=0., e, tx, ty, tz, dummy = GetParametrizedMess();
  //here _direction is three dimensional. dimension[0,1] is a normalized vector
  //in eta-phi plane in pT direction. _direction[2] is the p_z component of 
  //the measurement.
  
  x = dummy * _direction[0];
  y = dummy * _direction[1];
  if (GetMess()->pt!=0.) z = dummy/GetMess()->pt * _direction[2];
  e = sqrt( x*x + y*y + z*z );
  for (std::vector<TData_MessMess*>::const_iterator it=_m2.begin();
       it!=_m2.end();++it){
    dummy = (*it)->GetParametrizedMess();
    tx = dummy * (*it)->GetDirection()[0];
    ty = dummy * (*it)->GetDirection()[1];  
    if ((*it)->GetMess()->pt!=0.) tz = dummy/(*it)->GetMess()->pt * (*it)->GetDirection()[2];  
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
    

//  Scale the normalized residual 'z^2 = chi^2/weight' using 
//  the Cauchy-Function:  z^2 --> c^2 * Ln( 1 + (z/c)^2 )
double TData::ScaleCauchy(double z2)
{
  //double c = 2.3849;
  double c2 = 5.687748;
  //return c*c * log( 1 + z2/(c*c) );
  return c2 * log( 1. + z2/c2 );
}

//  Scale the normalized residual 'z^2 = chi^2/weight' using 
//  the Huber-Function   z^2 --> z                 for |z| <= c
//	                       c * ( 2*|z| - c ) for |z| > c
double TData::ScaleHuber(double z2)
{
  double c = 1.345;
  double z = sqrt(z2);
  return (  fabs(z) <= c  ?  z2  :  c*(2.*fabs(z) - c)  );
}
