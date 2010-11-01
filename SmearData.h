// $Id: SmearData.h,v 1.12 2010/09/22 13:29:44 mschrode Exp $

#ifndef SmearData_h
#define SmearData_h

#include "CalibData.h"
#include "SmearFunction.h"


//!  \brief Abstract base class for jetsmearing method
//!  \author Matthias Schroeder
//!  \date Tue Jun  9 15:24:49 CEST 2009
//!  $Id: SmearData.h,v 1.12 2010/09/22 13:29:44 mschrode Exp $
// --------------------------------------------------
class SmearData : public Event {
 public:
  SmearData(DataType type, Measurement * mess, double truth, double ptHat, double weight, const SmearFunction& pdf);
  virtual ~SmearData() { delete mess_; }

  //!  \brief Get the negative log-likelihood of this event
  //!  \return The negative log-likelihood of this event
  // --------------------------------------------------
  virtual double chi2() const = 0;
  virtual double chi2_fast(double * temp_derivative1, double * temp_derivative2, const double* epsilon) const = 0;
  virtual void printInitStats() const = 0;

  virtual void setParameters(Parameters* param);
  virtual Measurement * mess() const { return mess_; }
  virtual double truth() const { return kTruth_; }
  virtual DataType type() const { return kType_; }
  double par(int i) { return pdf_.par(i); }
  double pdfPtMeas(double ptMeas1, double ptMeas2, double ptTrue) const { return pdf_.pdfPtMeas(ptMeas1,ptMeas2,ptTrue); }
  double pdfPtTrue(double ptTrue) const { return pdf_.pdfPtTrue(ptTrue); }
  double pdfPtTrueError(double ptTrue) const { return pdf_.pdfPtTrueError(ptTrue); }
  double pdfResp(double r, double ptTrue) const { return pdf_.pdfResp(r,ptTrue); }
  double pdfRespError(double r, double ptTrue) const { return pdf_.pdfRespError(r,ptTrue); }
  double pdfDijetAsym(double a, double ptTrue) const { return pdf_.pdfDijetAsym(a,ptTrue); }

  virtual double chi2_plots() const { return 0.; }                 //!< Dummy, no functionality
  virtual double parametrizedMess() const { return 0.; }        //!< Dummy, no functionality
  virtual void updateError() { }                                   //!< Dummy, no functionality


 protected:
  SmearFunction pdf_; 
  Measurement * mess_;


 private:
  const double    kTruth_;                     //!< Truth
  const DataType  kType_;                      //!< Event type
};

#endif
