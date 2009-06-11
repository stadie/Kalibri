// $Id: SmearData.cc,v 1.11 2009/06/10 14:19:56 mschrode Exp $

#include "SmearData.h"



//!  \brief Response pdf
//!  \param r Response
//!  \return The probability density of the response r
// --------------------------------------------------
double SmearData::RespPDF(double r) const {
  TJet jet;
  jet.pt = r;
  return mRespPDF(&jet);
}
