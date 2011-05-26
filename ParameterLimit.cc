#include "ParameterLimit.h"
#include "Parameters.h"

double ParameterLimit::chi2() const 
{  
    double diff = 0;
    if(par_->parameters()[index_] < min_) diff = min_-par_->parameters()[index_];
    if(par_->parameters()[index_] > max_) diff = par_->parameters()[index_]-max_;
    //return 1e4/(1+exp(k* (par[0] - min))) + 1e4/(1+exp(-k* (par[0] - max));
    return diff * diff / error2_;
};

void ParameterLimit::setParameters(Parameters* param) 
{
  par_ = param;
}

double ParameterLimit::chi2_fast(double* temp_derivative1, 
				 double* temp_derivative2, 
				 const double* epsilon) const
{
  
  // Penalty term with current parameter values
  double new_chi2  = chi2();
 
  // Variation of parameters
  double oldpar              = par_->parameters()[index_];
  par_->parameters()[index_] += epsilon[index_];
  double temp2               = chi2();
  par_->parameters()[index_]  = oldpar - epsilon[index_];
  double temp1     = chi2();
 
  // Difference of chi2 at par+epsilon and par-epsilon
  temp_derivative1[index_] += (temp2 - temp1);                // for 1st derivative
  temp_derivative2[index_] += (temp2 + temp1 - 2.*new_chi2);  // for 2nd derivative
 
  // Reset original parameter value
  par_->parameters()[index_]  = oldpar;
 
  return new_chi2;
}
