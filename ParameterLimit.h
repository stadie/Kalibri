#ifndef PARAMETERLIMIT_H
#define PARAMETERLIMIT_H


#include "CalibData.h"

class Parameters;

//!  \brief Data class to limit a parameter
class ParameterLimit : public Event
{
 public:
  ParameterLimit(unsigned short int index, Parameters* parameters, double min, double max, double error = 1.0)
   : Event(1.0,0.0), index_(index),  min_(min), max_(max), 
    error2_(error*error), par_(parameters) 
    { 
    }
  
  double chi2() const;
  
  double chi2_fast(double* temp_derivative1, double* temp_derivative2,
		   double * temp_derivative3, double * temp_derivative4,
		   const double* epsilon) const;
  Measurement *mess() const { return 0;}
  double truth() const { return 0;}
  double parametrizedMess() const { return 0;}
  void setParameters(Parameters* param); 
  DataType type() const { return ParLimit;}
  double chi2_plots() const { return chi2();}
  void updateError() {}
  
 private:
  unsigned short int index_;
  double min_, max_, error2_;
  Parameters* par_;
};

#endif
