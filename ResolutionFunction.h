// $ Id: $

#ifndef RESOLUTION_FUNCTION_H
#define RESOLUTION_FUNCTION_H

#include <vector>

#include "Function.h"

class ResolutionParametrization;

class ResolutionFunction : public Function {
  typedef double (ResolutionParametrization::*PdfPtMeas)(double ptMeas1, double ptMeas2, double ptTrue, unsigned int ptBin, const double*) const;
  typedef double (ResolutionParametrization::*PdfPtTrue)(double ptTrue, unsigned int ptBin) const;
  typedef double (ResolutionParametrization::*PdfResp)(double r, double ptTrue, unsigned int ptBin, const double*) const;
  typedef double (ResolutionParametrization::*PdfDijetAsym)(double a, double ptTrue, unsigned int ptBin, const double*) const;
  typedef double (ResolutionParametrization::*DMeasMax)(unsigned int ptBin) const;

 public:
  ResolutionFunction(unsigned int nPars, unsigned int ptBin, unsigned int parIdx, double *firstPar,
		     const std::vector<bool>& isFixedPar,
		     PdfPtMeas pdfPtMeas, PdfPtTrue pdfPtTrue, PdfResp pdfResp, PdfDijetAsym pdfDijetAsym,
		     DMeasMax dMeasMax, const ResolutionParametrization *p)
    : Function(0,0,firstPar,parIdx,nPars,0),
    ptBin_(ptBin), isFixedPar_(isFixedPar),
    pdfPtMeas_(pdfPtMeas), pdfPtTrue_(pdfPtTrue), pdfResp_(pdfResp), pdfDijetAsym_(pdfDijetAsym),
    dMeasMax_(dMeasMax) {
    param_ = p;
  }
  
  Function* clone() {
    return new ResolutionFunction(*this);
  }

  unsigned int ptBin() const { return ptBin_; }
  bool isFixedPar(unsigned int i) const { return isFixedPar_[i]; }
  
  double pdfPtMeas(double ptMeas1, double ptMeas2, double ptTrue) const {
    return (param_->*pdfPtMeas_)(ptMeas1,ptMeas2,ptTrue,ptBin_,firstpar_);
  }
  double pdfPtTrue(double ptTrue) const {
    return (param_->*pdfPtTrue_)(ptTrue,ptBin_);
  }
  double pdfResp(double r, double ptTrue) const {
    return (param_->*pdfResp_)(r,ptTrue,ptBin_,firstpar_);
  }
  double pdfDijetAsym(double a, double ptTrue) const {
    return (param_->*pdfDijetAsym_)(a,ptTrue,ptBin_,firstpar_);
  }
  double dMeasMax() const { return (param_->*dMeasMax_)(ptBin_); }
  

 private:
    static const ResolutionParametrization* param_;

    const unsigned int ptBin_;
    const std::vector<bool> isFixedPar_;

    const PdfPtMeas pdfPtMeas_;
    const PdfPtTrue pdfPtTrue_;
    const PdfResp pdfResp_;
    const PdfDijetAsym pdfDijetAsym_;
    const DMeasMax dMeasMax_;
};
#endif
