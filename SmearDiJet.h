// $Id: SmearDiJet.h,v 1.1 2009/06/11 17:29:25 mschrode Exp $

#ifndef SmearDiJet_h
#define SmearDiJet_h

#include "CalibData.h"
#include "SmearData.h"


//!  \brief Dijet data for jetsmearing method
//!  \author Matthias Schroeder
//!  \date Tue Jun  9 18:23:44 CEST 2009
//!  $Id: SmearDiJet.h,v 1.1 2009/06/11 17:29:25 mschrode Exp $
// --------------------------------------------------
class SmearDiJet : public SmearData {
 public:
  SmearDiJet(TMeasurement * mess,
	     TMeasurement * secndMess,
	     TMeasurement * thirdMess,
	     double weight,
	     const Function& respPDF,
	     const Function& truthPDF,
	     double min,
	     double max,
	     double eps,
	     int niter);
  ~SmearDiJet() { delete secndMess_; delete thirdMess_; }

  virtual void ChangeParAddress(double* oldpar, double* newpar);
  virtual double chi2() const;
  virtual double chi2_fast(double * temp_derivative1, double * temp_derivative2, double const epsilon) const;
  virtual void PrintInitStats() const;

  double dijetPt() const { return 0.5 * (GetMess()->pt + GetSecondMess()->pt); }
  TMeasurement * GetSecondMess() const { return secndMess_; }  //!< Get second jet
  TMeasurement * GetThirdMess() const { return thirdMess_; }   //!< Get third jet
  double * GetTruthPar() { return truthPDF_.firstPar(); }
  double TruthPDF(double t) const;


 private:
  const int    kMaxNIter_;   //!< Max number of iterations in integration
  const double kEps_;        //!< Integration precision for convergence
  const double kMin_;        //!< Minimum of truth pdf
  const double kMax_;        //!< Maximum of truth pdf

  TMeasurement * secndMess_; //!< Second jet
  TMeasurement * thirdMess_; //!< Third jet
  Function       truthPDF_;  //!< Truth pdf
};
#endif
