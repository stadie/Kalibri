// $Id: SmearData.cc,v 1.10 2010/05/19 13:34:48 stadie Exp $

#include "SmearData.h"

#include "CalibData.h"

#include "Parameters.h"
//!  \brief Constructor
//!  \param type Data type
//!  \param mess The measurement
//!  \param truth The truth of the measurement
//!  \param weight Event weight
//!  \param respPDF Response probability density
// --------------------------------------------------
SmearData::SmearData(DataType type, Measurement * mess, double truth, double ptHat, double weight, const SmearFunction& pdf)
  : Event(weight,ptHat),pdf_(pdf),mess_(mess),kTruth_(truth),kType_(type)
{
};


void SmearData::setParameters(Parameters* param) { 
  double *oldpar = pdf_.par() - pdf_.parIdx();
   pdf_.changeParBase(oldpar,param->parameters()); 
 }
