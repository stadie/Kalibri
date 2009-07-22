// $Id: SmearData.cc,v 1.2 2009/07/16 14:45:53 mschrode Exp $

#include "SmearData.h"



//!  \brief Constructor
//!  \param type Data type
//!  \param mess The measurement
//!  \param truth The truth of the measurement
//!  \param weight Event weight
//!  \param respPDF Response probability density
// --------------------------------------------------
SmearData::SmearData(DataType type, TMeasurement * mess, double truth, double weight, const Function& respPDF)
    : TData(),
    respPDF_(respPDF),
    kTruth_(truth),
    kType_(type),
    kWeight_(weight),
    mess_(mess) {};



//!  \brief Response pdf
//!  \param r Response
//!  \return The probability density of the response r
// --------------------------------------------------
double SmearData::RespPDF(double r) const {
  TJet jet;
  jet.pt = r;
  return respPDF_(&jet);
}
