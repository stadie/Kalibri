//!    \brief Different resolution parametrizations for unbinned fit
//!
//!    \author Matthias Schroeder
//!
//!    \date 2010
//!
//! $Id: $

#ifndef RESOLUTION_PARAMETRIZATION_H
#define RESOLUTION_PARAMETRIZATION_H

#include <cmath>
#include <iostream>
#include <string>
#include <vector>


class TH1;
class TRandom;


//! Abstract base of response function parametrizations
//! for resolution fit
// ------------------------------------------------------------------------
class ResolutionParametrization {
 public:
  class CrystalBallFunction {
  public:
    CrystalBallFunction();
    ~CrystalBallFunction();

    double integral(double mean, double sigma, double alpha, double n, double min, double max) const;	
    double norm(double mean, double sigma, double alpha, double n) const;
    double pdf(double x, double mean, double sigma, double alpha, double n) const {
      return norm(mean,sigma,alpha,n)*value(x,mean,sigma,alpha,n);
    }
    double random(double mean, double sigma, double alpha, double n) const;
    double value(double x, double mean, double sigma, double alpha, double n) const;

    double truncPdf(double x, double mean, double sigma, double alpha, double n, double min) const;
    double truncRandom(double mean, double sigma, double alpha, double n, double min) const;
    double truncValue(double x, double mean, double sigma, double alpha, double n, double min) const {
      return x > min ? value(x,mean,sigma,alpha,n) : 0.;
    }
    
	
  private:
    TRandom *rand_;
    std::vector<double> par_;
  };


  ResolutionParametrization(unsigned int nPtBins, unsigned int nParPerPtBin);
    virtual ~ResolutionParametrization() {};
    
  virtual std::string name() const = 0;
  virtual bool needsUpdate() const = 0;
  virtual void update(const double *par) = 0;

  //! Returns probability density of measured jet pts given a true jet pt
  virtual double pdfPtMeas(double ptMeas1, double ptMeas2, double ptTrue, unsigned int ptBin, const double *par) const = 0;
  //! Returns probability density of true jet pt
  virtual double pdfPtTrue(double ptTrue, unsigned int ptBin) const = 0;
  virtual double pdfResp(double r, double ptTrue, unsigned int ptBin, const double *par) const = 0;
  virtual double pdfDijetAsym(double a, double ptTrue, unsigned int ptBin, const double *par) const = 0;

  //! Print some initialization details
  virtual void print() const { std::cout << "Parametrization class '" << name() << "'\n"; }  
  unsigned int nPtBins() const { return nPtBins_; }
  unsigned int nParPerPtBin() const { return nParPerPtBin_; }
  double dMeasMax(unsigned int ptBin) const { return dMeasMax_.at(ptBin); }


 protected:
  std::vector<double> dMeasMax_;


 private:
  const unsigned int nPtBins_;
  const unsigned int nParPerPtBin_;
};



// ------------------------------------------------------------------------
class ResolutionEmpty : public ResolutionParametrization {
 public:
  ResolutionEmpty(unsigned int nPtBins) : ResolutionParametrization(nPtBins,0) {};

  std::string name() const { return "ResolutionEmpty"; }
  bool needsUpdate() const { return false; }
  void update(const double *par) { }

  double pdfPtMeas(double ptMeas1, double ptMeas2, double ptTrue, unsigned int ptBin, const double *par) const { return 0.; }
  double pdfPtTrue(double ptTrue, unsigned int ptBin) const { return 0.; }
  double pdfResp(double r, double ptTrue, unsigned int ptBin, const double *par) const { return 0.; }
  double pdfDijetAsym(double a, double ptTrue, unsigned int ptBin, const double *par) const { return 0.; }
};



//! \brief Gaussian function with one sigma (not pt-dependent)
//!        and the assumption pttrue = ptave
// ------------------------------------------------------------------------
class ResolutionGaussAvePt : public ResolutionParametrization
{ 
 public:
  //! Constructor
  ResolutionGaussAvePt(unsigned int nPtBins);

  std::string name() const { return "ResolutionGaussAvePt"; }
  bool needsUpdate() const { return true; }
  void update(const double * par);

  double pdfPtMeas(double ptMeas1, double ptMeas2, double ptTrue, unsigned int ptBin, const double *par) const;
  double pdfPtTrue(double ptTrue, unsigned int ptBin) const { return 0.; }
  double pdfResp(double r, double ptTrue, unsigned int ptBin, const double *par) const;
  double pdfDijetAsym(double a, double ptTrue, unsigned int ptBin, const double *par) const;


 private:
  double sigma(const double *par) const {
    //return (par[0] > 1./sqrt(M_PI) ? par[0] : 1./sqrt(M_PI));
    return (par[0] > 2./sqrt(M_PI) ? par[0] : 2./sqrt(M_PI));
  }

  void print() const;
};



//! \brief Gaussian function with one sigma (not pt-dependent)
// ------------------------------------------------------------------------
class ResolutionGauss : public ResolutionParametrization
{ 
 public:
  ResolutionGauss(unsigned int nPtBins, const std::vector<double> &ptAveBinEdges, const std::vector<TH1*> &spectra);
  ~ResolutionGauss();

  std::string name() const { return "ResolutionGauss";}
  bool needsUpdate() const { return true; }
  void update(const double * par);

  double pdfPtMeas(double ptMeas1, double ptMeas2, double ptTrue, unsigned int ptBin, const double *par) const;
  double pdfPtTrue(double ptTrue, unsigned int ptBin) const;
  double pdfResp(double r, double ptTrue, unsigned int ptBin, const double *par) const;
  double pdfDijetAsym(double a, double ptTrue, unsigned int ptBin, const double *par) const;


 private:
  const std::vector<double> ptAveBinEdges_;
  const std::vector<TH1*> spectra_;

  std::vector<TH1*> hashTablePdfPtTrue_;

  double ptAveMin(unsigned int ptBin) const { return ptAveBinEdges_.at(ptBin); }
  double ptAveMax(unsigned int ptBin) const { return ptAveBinEdges_.at(ptBin+1); }
  double sigma(const double *par) const {
    return ( par[0] > 2.*sqrt(2./M_PI) ? par[0] : 2.*sqrt(2./M_PI) );
  }
  double sigma(const double *par, unsigned int ptBin) const {
    return ( par[ptBin] > 2.*sqrt(2./M_PI) ? par[ptBin] : 2.*sqrt(2./M_PI) );
  }
  double underlyingPdfPtTrue(double ptTrue, unsigned int ptBin) const;

  void print() const;
};



#endif
