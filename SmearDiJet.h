// $Id: SmearDiJet.h,v 1.10 2010/07/22 13:58:30 mschrode Exp $

#ifndef SmearDiJet_h
#define SmearDiJet_h

#include <vector>

#include "SmearData.h"
#include "SmearFunction.h"
#include "Jet.h"


//!  \brief Dijet data for jetsmearing method
//!  \author Matthias Schroeder
//!  \date Tue Jun  9 18:23:44 CEST 2009
//!  $Id: SmearDiJet.h,v 1.10 2010/07/22 13:58:30 mschrode Exp $
// --------------------------------------------------
class SmearDiJet : public SmearData {
 public:
  SmearDiJet(Jet * jet1, const Jet * jet2, const Jet * jet3,
	     double deltaPhi12, double pPhi, double ptJet3, double ptJet4, 
	     double pJ3, double pSJ, double pUCE, double ptRef,
	     double ptHat, double weight,
	     const SmearFunction& pdf,
	     double min, double max, double eps, int niter);
  ~SmearDiJet();

  //  virtual void ChangeParAddress(double* oldpar, double* newpar);
  virtual double chi2_fast(double * temp_derivative1, double * temp_derivative2, const double* epsilon) const;
  virtual double chi2() const { return chi2Spectrum(); }
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
  double deltaPhi12() const { return deltaPhi12_; }
  double pPhi() const { return pPhi_; }
  double pJ3() const { return pJ3_; }
  double pSJ() const { return pSJ_; }
  double pUCE() const { return pUCE_; }
  double ptRef() const { return ptRef_; }

  double dijetPt() const { return 0.5 * (jet1()->pt() + jet2()->pt()); } 
  double avePt() const { return 0.5 * (jet1()->pt() + jet2()->pt()); } 
  double avePtGen() const { return 0.5 * (jet1()->genPt() + jet2()->genPt()); } 
  double relJet3Pt() const { return jet3()->pt() / dijetPt(); }
  double ptJet3() const { return ptJet3_; }
  double ptJet4() const { return ptJet4_; }


 private:
  const int    kMaxNIter_;   //!< Max number of iterations in integration
  const double kEps_;        //!< Integration precision for convergence
  const double kMin_;        //!< Minimum of truth pdf
  const double kMax_;        //!< Maximum of truth pdf

  const Jet * jet2_; //!< Second jet
  const Jet * jet3_; //!< Third jet
  const double deltaPhi12_;
  const double pPhi_;
  const double pJ3_;
  const double pSJ_;
  const double pUCE_;
  const double ptRef_;
  const double ptJet3_;
  const double ptJet4_;
  
  double chi2Simple() const;
  double chi2Spectrum() const;
};
#endif
