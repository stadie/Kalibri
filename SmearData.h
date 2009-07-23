// $Id: SmearData.h,v 1.3 2009/07/22 11:46:18 mschrode Exp $

#ifndef SmearData_h
#define SmearData_h

#include "CalibData.h"
#include "Function.h"


//!  \brief Abstract base class for jetsmearing method
//!  \author Matthias Schroeder
//!  \date Tue Jun  9 15:24:49 CEST 2009
//!  $Id: SmearData.h,v 1.3 2009/07/22 11:46:18 mschrode Exp $
// --------------------------------------------------
class SmearData : public TData {
 public:
  SmearData(DataType type, TMeasurement * mess, double truth, double weight, const Function& respPDF);
  virtual ~SmearData() { delete mess_; }

  //!  \brief Get the negative log-likelihood of this event
  //!  \return The negative log-likelihood of this event
  // --------------------------------------------------
  virtual double chi2() const = 0;
  virtual double chi2_fast(double * temp_derivative1, double * temp_derivative2, double const epsilon) const = 0;
  virtual void printInitStats() const = 0;
  virtual double ptHat() const = 0;


  virtual void ChangeParAddress(double* oldpar, double* newpar) { respPDF_.changeParBase(oldpar,newpar); }
  virtual TMeasurement * GetMess() const { return mess_; }
  virtual double GetTruth() const { return kTruth_; }
  virtual DataType GetType() const { return kType_; }
  virtual double GetWeight() const { return weight_; }
  double * GetRespPar() { return respPDF_.firstPar(); }
  double RespPDF(double r) const;
  virtual void setWeight(double w) { weight_ = w; } 

  virtual double chi2_plots() const { return 0.; }                 //!< Dummy, no functionality
  virtual double GetParametrizedMess() const { return 0.; }        //!< Dummy, no functionality
  virtual void UpdateError() { }                                   //!< Dummy, no functionality


 protected:
  Function       respPDF_;                    //!< Response pdf


 private:
  const double    kTruth_;                     //!< Truth
  const DataType  kType_;                      //!< Event type
  double          weight_;                     //!< Event weight
  TMeasurement  * mess_;                       //!< The measurement
};

#endif
