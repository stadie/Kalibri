// $Id: SmearData.cc,v 1.3 2009/07/22 11:46:18 mschrode Exp $

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
    weight_(weight),
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
