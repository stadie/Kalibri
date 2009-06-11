// $Id: SmearData.h,v 1.7 2009/06/10 14:19:56 mschrode Exp $

#ifndef SmearData_h
#define SmearData_h

#include "CalibData.h"
#include "Function.h"


//!  \brief Abstract base class for jetsmearing method
//!  \author Matthias Schroeder
//!  \date Tue Jun  9 15:24:49 CEST 2009
//!  $Id: SmearData.h,v 1.7 2009/06/10 14:19:56 mschrode Exp $
// --------------------------------------------------
class SmearData : public TData {
 public:
  SmearData(DataType type, TMeasurement * mess, double truth, double weight, const Function& respPDF)
    : TData(),
    mRespPDF(respPDF),
    kTruth(truth),
    kType(type),
    kWeight(weight),
    mMess(mess) {};
  virtual ~SmearData() { delete mMess; }

  //!  \brief Get the negative log-likelihood of this event
  //!  \return The negative log-likelihood of this event
  // --------------------------------------------------
  virtual double chi2() const = 0;
  virtual double chi2_fast(double * temp_derivative1, double * temp_derivative2, double const epsilon) const = 0;
  virtual void PrintInitStats() const = 0;


  virtual void ChangeParAddress(double* oldpar, double* newpar) { mRespPDF.changeParBase(oldpar,newpar); }
  virtual TMeasurement * GetMess() const { return mMess; }
  virtual double GetTruth() const { return kTruth; }
  virtual DataType GetType() const { return kType; }
  virtual double GetWeight() const { return kWeight; }
  double * GetRespPar() { return mRespPDF.firstPar(); }
  double RespPDF(double r) const;

  virtual double chi2_plots() const { return 0.; }                 //!< Dummy, no functionality
  virtual double GetParametrizedMess() const { return 0.; }        //!< Dummy, no functionality
  virtual void UpdateError() { }                                   //!< Dummy, no functionality


 protected:
  Function       mRespPDF;                    //!< Response pdf


 private:
  const double    kTruth;                     //!< Truth
  const DataType  kType;                      //!< Event type
  const double    kWeight;                    //!< Event weight
  TMeasurement * mMess;                       //!< The measurement
};

#endif
