// $Id: SmearDiJet.h,v 1.9 2010/05/19 13:34:49 stadie Exp $

#ifndef SmearDiJet_h
#define SmearDiJet_h

#include <vector>

#include "SmearData.h"
#include "SmearFunction.h"
#include "Jet.h"


//!  \brief Dijet data for jetsmearing method
//!  \author Matthias Schroeder
//!  \date Tue Jun  9 18:23:44 CEST 2009
//!  $Id: SmearDiJet.h,v 1.9 2010/05/19 13:34:49 stadie Exp $
// --------------------------------------------------
class SmearDiJet : public SmearData {
 public:
  SmearDiJet(Jet * jet1,
	     Jet * jet2,
	     Jet * jet3,
	     double ptHat,
	     double weight,
	     const SmearFunction& pdf,
	     double min,
	     double max,
	     double eps,
	     int niter);
  ~SmearDiJet();

  //  virtual void ChangeParAddress(double* oldpar, double* newpar);
  virtual double chi2() const;
  virtual double chi2_fast(double * temp_derivative1, double * temp_derivative2, const double* epsilon) const;
  virtual void printInitStats() const;

  const Jet * jet1() const { return static_cast<Jet*>(mess_); }
  const Jet * jet2() const { return jet2_; }
  const Jet * jet3() const { return jet3_; }
  const Jet * jet(int i) const { 
    const Jet *jet = 0;
    if( i == 0 ) jet = jet1();
    else if( i == 1 ) jet = jet2();
    else if( i == 2 ) jet = jet3();
    return jet;
  }
  double pdfPtMeasJet1(double ptMeas, double ptTrue, double pt3Rel) const {
    return pdf_.pdfPtMeasJet1(ptMeas,ptTrue,pt3Rel);
  }
  double pdfPtMeasJet2(double ptMeas, double ptTrue, double pt3Rel) const {
    return pdf_.pdfPtMeasJet2(ptMeas,ptTrue,pt3Rel);
  }

  //! Get dijet pt \f$ \frac{1}{2} (p^{1}_{T} + p^{2}_{T}) \f$
  double dijetPt() const { return 0.5 * (jet1()->pt() + jet2()->pt()); } 
  double avePt() const { return 0.5 * (jet1()->pt() + jet2()->pt()); } 
  double avePtGen() const { return 0.5 * (jet1()->genPt() + jet2()->genPt()); } 
  double relJet3Pt() const { return jet3()->pt() / dijetPt(); }
  
  double scalePt2() const { return scalePt2_; }
  double scalePt3() const { return scalePt3_; }
  double relGenMet() const { return relGenMet_; }


 private:
  const int    kMaxNIter_;   //!< Max number of iterations in integration
  const double kEps_;        //!< Integration precision for convergence
  const double kMin_;        //!< Minimum of truth pdf
  const double kMax_;        //!< Maximum of truth pdf

  Jet * jet2_; //!< Second jet
  Jet * jet3_; //!< Third jet
  
  double scalePt2_;
  double scalePt3_;
  double relGenMet_;
};
#endif
