// $Id: SmearDiJet.h,v 1.1 2009/06/10 14:20:49 mschrode Exp $

#ifndef SmearDiJet_h
#define SmearDiJet_h

#include "CalibData.h"
#include "SmearData.h"


//!  \brief Dijet data for jetsmearing method
//!  \author Matthias Schroeder
//!  \date Tue Jun  9 18:23:44 CEST 2009
//!  $Id: SmearDiJet.h,v 1.1 2009/06/10 14:20:49 mschrode Exp $
// --------------------------------------------------
class SmearDiJet : public SmearData {
 public:
  SmearDiJet(TMeasurement * mess, TMeasurement * scndmess, double weight,
	     const Function& respPDF, const Function& truthPDF,
	     double min, double max, double eps, int niter)
    : SmearData(TypeSmearDiJet,mess,0,weight,respPDF),
    kMaxNIter(niter),
    kEps(eps),
    kMin(min),
    kMax(max),
    mScndMess(scndmess),
    mTruthPDF(truthPDF) {};
  ~SmearDiJet() { delete mScndMess; }

  virtual void ChangeParAddress(double* oldpar, double* newpar);
  virtual double chi2() const;
  virtual double chi2_fast(double * temp_derivative1, double * temp_derivative2, double const epsilon) const;
  virtual void PrintInitStats() const;

  TMeasurement * GetSecondMess() const { return mScndMess; }  //!< Get second jet
  double * GetTruthPar() { return mTruthPDF.firstPar(); }
  double TruthPDF(double t) const;


 private:
  const int    kMaxNIter;   //!< Max number of iterations in integration
  const double kEps;        //!< Integration precision for convergence
  const double kMin;        //!< Minimum of truth pdf
  const double kMax;        //!< Maximum of truth pdf

  TMeasurement * mScndMess; //!< Second jet
  Function       mTruthPDF; //!< Truth pdf
};
#endif
