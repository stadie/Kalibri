#include "ResolutionParametrization.h"

#include <cassert>

#include "TH1.h"
#include "TH1D.h"
#include "TRandom.h"
#include "TRandom3.h"




// ------------------------------------------------------------------------
ResolutionParametrization::ResolutionParametrization(unsigned int nPtBins, unsigned int nParPerPtBin)
  : nPtBins_(nPtBins), nParPerPtBin_(nParPerPtBin) {
  dMeasMax_ = std::vector<double>(nPtBins,0.);
}



// ------------------------------------------------------------------------
ResolutionGaussAvePt::ResolutionGaussAvePt(unsigned int nPtBins)
  : ResolutionParametrization(nPtBins,1) {
  dMeasMax_ = std::vector<double>(nPtBins);
  print();
}


// ------------------------------------------------------------------------
void ResolutionGaussAvePt::update(const double * par) {
  std::cout << name() << ": Updating hashed parameters" << std::endl;
  std::cout << "  Updating maximum dMeas " << std::endl;
  for(unsigned int i = 0; i < nPtBins(); ++i) {
    dMeasMax_[i] = 2.*par[i]/sqrt(2.); // Two sigma in dMeas
    std::cout << "    Bin " << i+1 << ": sigma = " << par[i] << " GeV  -->  dMeasMax = " << dMeasMax_[i] << " GeV\n";
  }
}


// ------------------------------------------------------------------------
double ResolutionGaussAvePt::pdfPtMeas(double ptMeas1, double ptMeas2, double ptTrue, unsigned int ptBin, const double *par) const {
//   double pdf = 0.;
//   double dMeas = 0.5*(ptMeas1 - ptMeas2);
//   if( std::abs(dMeas) < dMeasMax(ptBin) ) {
//     double s = sigma(par)/sqrt(2.);
//     double u = dMeas/s;
//     double norm = s*sqrt(M_PI*2.)*erf(dMeasMax(ptBin)/sqrt(2.)/s);
//     pdf = exp(-0.5*u*u)/norm;
//   }
//   return pdf;

  double pdf = 0.;
  double dMeas = 0.5*std::abs(ptMeas1 - ptMeas2);
  if( dMeas < dMeasMax(ptBin) ) {
    double s = sigma(par)/sqrt(2.);
    double u = dMeas/s;
    double norm = s*sqrt(M_PI/2.)*erf(dMeasMax(ptBin)/sqrt(2.)/s);
    pdf = exp(-0.5*u*u)/norm;
  }
  return pdf;

}


// ------------------------------------------------------------------------
double ResolutionGaussAvePt::pdfResp(double r, double ptTrue, unsigned int ptBin, const double *par) const {
  double s = sigma(par)/ptTrue;
  double u = (r - 1.)/s;
  return exp(-0.5*u*u)/sqrt(2.*M_PI)/s;
}


// ------------------------------------------------------------------------
double ResolutionGaussAvePt::pdfDijetAsym(double a, double ptTrue, unsigned int ptBin, const double *par) const {
  double s = sigma(par)/ptTrue/sqrt(2.);
  double u = a/s;
  return exp(-0.5*u*u)/sqrt(2.*M_PI)/s;
}


// ------------------------------------------------------------------------
void ResolutionGaussAvePt::print() const {
  std::cout << "Parametrization class '" << name() << "'\n";
  std::cout << "  " << nPtBins() << " ptAve bins\n";
  std::cout << "  " << nParPerPtBin() << " parameters per ptAve bin\n";
  std::cout << std::endl;
}





// ------------------------------------------------------------------------
ResolutionGauss::ResolutionGauss(unsigned int nPtBins, const std::vector<double> &ptAveBinEdges, const std::vector<TH1*> &spectra)
  : ResolutionParametrization(nPtBins,1), ptAveBinEdges_(ptAveBinEdges), spectra_(spectra) {

  assert( ptAveBinEdges_.size() == nPtBins+1 );
  assert( spectra_.size() == nPtBins );

  hashTablePdfPtTrue_ = std::vector<TH1*>(nPtBins);
  for(unsigned int i = 0; i < nPtBins; ++i) {
    TString name = "hashTablePdfPtTrue_";
    name += i;
    hashTablePdfPtTrue_[i] = new TH1D(name,"",5000,spectra_[i]->GetXaxis()->GetBinLowEdge(1),spectra_[i]->GetXaxis()->GetBinUpEdge(spectra_[i]->GetNbinsX()));
  }

  print();
}

// ------------------------------------------------------------------------
ResolutionGauss::~ResolutionGauss() { 
  for(unsigned int i = 0; i < nPtBins(); ++i) {
    delete spectra_[i];
    delete hashTablePdfPtTrue_[i];
  }
}



// ------------------------------------------------------------------------
void ResolutionGauss::update(const double * par) {
   std::cout << name() << ": Updating hashed parameters" << std::endl;

   // Updating dMeas
   std::cout << "  Updating maximum dMeas" << std::endl;
   for(unsigned int i = 0; i < nPtBins(); ++i) {
     dMeasMax_[i] = 2.*par[i]/sqrt(2.); // Two sigma in dMeas
     std::cout << "    Bin " << i+1 << ": sigma = " << par[i] << " GeV  -->  dMeasMax = " << dMeasMax_[i] << " GeV\n";
   }

   // Updating truth spectrum
   std::cout << "  Updating truth pdf" << std::endl;
   // Loop over pt bins
   for(unsigned int ptBin = 0; ptBin < nPtBins(); ++ptBin) {
     // Loop over tMin < ptTrue < tMax in this pt bin
     for(int tBin = 1; tBin <= hashTablePdfPtTrue_[ptBin]->GetNbinsX(); ++tBin) {
       double ptTrue = hashTablePdfPtTrue_[ptBin]->GetBinCenter(tBin);
       // Convolution with cuts on ptAve
       double s = sigma(par,ptBin)/sqrt(2.);
       double c = sqrt(M_PI/2.)*s*( erf((ptAveMax(ptBin)-ptTrue)/s/sqrt(2.)) - erf((ptAveMin(ptBin)-ptTrue)/s/sqrt(2.)) );
       double pdf = c*underlyingPdfPtTrue(ptTrue,ptBin);

       // Store (un-normalized) truth pdf for
       // this value of ptTrue in hash table
       hashTablePdfPtTrue_[ptBin]->SetBinContent(tBin,pdf);
     } // End of loop over ptTrue

     // Normalise values of truth pdf
     hashTablePdfPtTrue_[ptBin]->Scale(1./hashTablePdfPtTrue_[ptBin]->Integral("width"));
   } // End of loop over pt bins
}



// ------------------------------------------------------------------------
double ResolutionGauss::pdfPtMeas(double ptMeas1, double ptMeas2, double ptTrue, unsigned int ptBin, const double *par) const {
//   double pdf = 0.;
//    double dMeas = 0.5*(ptMeas1 - ptMeas2);
//    if( std::abs(dMeas) < dMeasMax(ptBin) ) {
//      double s = sigma(par)/sqrt(2.);
//      double u1 = dMeas/s;
//      double u2 = (0.5*(ptMeas1 + ptMeas2)-ptTrue)/s;
//      double norm = M_PI*s*s*erf(dMeasMax(ptBin)/sqrt(2.)/s)*( erf((ptAveMax(ptBin)-ptTrue)/sqrt(2.)/s) - erf((ptAveMin(ptBin)-ptTrue)/sqrt(2.)/s) );
//      if( norm < 1E-3 ) {
//        pdf = 0.;
//      } else {
//        pdf = exp(-0.5*u1*u1-0.5*u2*u2)/norm;
//      }
//    }
//   return pdf;

  double pdf = 0.;
  double dMeas = 0.5*std::abs(ptMeas1 - ptMeas2);
  if( dMeas < dMeasMax(ptBin) ) {
    double s = sigma(par)/sqrt(2.);
    double u1 = dMeas/s;
    double u2 = (0.5*(ptMeas1 + ptMeas2)-ptTrue)/s;
    double norm = M_PI*s*s*erf(dMeasMax(ptBin)/sqrt(2.)/s)/2.*( erf((ptAveMax(ptBin)-ptTrue)/sqrt(2.)/s) - erf((ptAveMin(ptBin)-ptTrue)/sqrt(2.)/s) );
    if( norm < 1E-3 ) {
      pdf = 0.;
    } else {
      pdf = exp(-0.5*u1*u1-0.5*u2*u2)/norm;
    }
  }
  return pdf;
}


// ------------------------------------------------------------------------
double ResolutionGauss::pdfResp(double r, double ptTrue, unsigned int ptBin, const double *par) const {
  double s = sigma(par)/ptTrue;
  double u = (r - 1.)/s;
  double cut = (1.+erf(ptTrue/sqrt(2.)/s))/2.;
  return exp(-0.5*u*u)/sqrt(2.*M_PI)/s/cut;
}


// ------------------------------------------------------------------------
double ResolutionGauss::pdfDijetAsym(double a, double ptTrue, unsigned int ptBin, const double *par) const {
  double s = sigma(par)/ptTrue/sqrt(2.);
  double u = a/s;
  return exp(-0.5*u*u)/sqrt(2.*M_PI)/s;
}


// ------------------------------------------------------------------------
double ResolutionGauss::pdfPtTrue(double ptTrue, unsigned int ptBin) const {
  return hashTablePdfPtTrue_.at(ptBin)->GetBinContent(hashTablePdfPtTrue_.at(ptBin)->FindBin(ptTrue));
}


// ------------------------------------------------------------------------
double ResolutionGauss::underlyingPdfPtTrue(double ptTrue, unsigned int ptBin) const {
  return spectra_.at(ptBin)->Interpolate(ptTrue);

//   // For ToyMC
//   return exp(-ptTrue/80.);
}


// ------------------------------------------------------------------------
void ResolutionGauss::print() const {
  std::cout << "Parametrization class '" << name() << "'\n";
  std::cout << "  PtAve binning (GeV): " << std::flush;
  for(unsigned int ptBin = 0; ptBin < nPtBins(); ++ptBin) {
    std::cout << ptAveMin(ptBin) << ", " << std::flush;
  }
  std::cout << ptAveMax(nPtBins()-1) << std::endl;
  std::cout << std::endl;
}




//!  A Crystal Ball function representation
// ------------------------------------------------------------------------
ResolutionParametrization::CrystalBallFunction::CrystalBallFunction()
  : rand_(new TRandom3(0)) {}
	
ResolutionParametrization::CrystalBallFunction::~CrystalBallFunction() {
  delete rand_;
}


//! The Crystal Ball function value (not normalized)
// ------------------------------------------------------------------------
double ResolutionParametrization::CrystalBallFunction::value(double x, double mean, double sigma, double alpha, double n) const {
  double f = 0.;
  if( x > 0. ) {
    double u = (x - mean)/sigma;
    if( u > -alpha ) {             // Gaussian part
      f = exp(-0.5*u*u);
    } else {                       // Powerlaw part
      f = exp(-0.5*alpha*alpha)/pow(1.-alpha*alpha/n-alpha*u/n,n); //A*pow(B-u,-n);
    }
  }

  return f;
}


//! The inverse of the integral over Crystal Ball function 
//! from 0 to infinity
// ------------------------------------------------------------------------
double ResolutionParametrization::CrystalBallFunction::norm(double mean, double sigma, double alpha, double n) const {
  double m = n-1.;
  double k = n*sigma*exp(-0.5*alpha*alpha)/m/alpha;
  double norm = sigma*sqrt(M_PI/2.)*( 1. + erf(alpha/sqrt(2)) );
  if( n == 1. ) {
    double B = n/alpha - alpha;
    norm += k/pow(1-alpha*alpha/n,m)*log( (B+mean/sigma)/(B+alpha) );
  } else {
    norm += k*( 1. - pow( 1 + alpha/n*( mean/sigma - alpha ),-m ) );
  }
  norm = 1./norm;
  if( norm < 1E-10 ) norm = 1E-10;
  return norm;
}


// ------------------------------------------------------------------------
double ResolutionParametrization::CrystalBallFunction::random(double mean, double sigma, double alpha, double n) const {
  double max = pdf(mean,mean,sigma,alpha,n);
  double x = -1.;
  do {
    x = rand_->Uniform(0.,2.);
  } while( pdf(x,mean,sigma,alpha,n) < rand_->Uniform(0.,1.)*max );
  return x;
}


// ------------------------------------------------------------------------
double ResolutionParametrization::CrystalBallFunction::integral(double mean, double sigma, double alpha, double n, double min, double max) const {
  double m = n - 1.;
  double c = mean - alpha*sigma;

  double in = 0.;
  if( min > c ) {
    // Integral from Gaussian part
    in = sigma*sqrt(M_PI/2.)*( erf((mean-min)/sqrt(2)/sigma) - erf((mean-max)/sqrt(2)/sigma) );
  } else if( max < c ) {
    // Integral from powerlaw part
    double k = n*sigma*exp(-0.5*alpha*alpha)/m/alpha;
    if( n == 1. ) {
      double B = n/alpha - alpha;
      in = k/pow(1-alpha*alpha/n,m)*log( (B+(mean-min)/sigma)/(B+(mean-max)/sigma) );
    } else {
      in = k*( pow( 1 + alpha/n*( (mean-max)/sigma - alpha ),-m )
	       -pow( 1 + alpha/n*( (mean-min)/sigma - alpha ),-m ) );
    }
  } else {
    // Integral from both parts
    double k = n*sigma*exp(-0.5*alpha*alpha)/m/alpha;
    in = sigma*sqrt(M_PI/2.)*( erf(alpha/sqrt(2)) - erf((mean-max)/sqrt(2)/sigma) );
    if( n == 1. ) {
      double B = n/alpha - alpha;
      in += k/pow(1-alpha*alpha/n,m)*log( (B+(mean-min)/sigma)/(B+alpha) );
    } else {
      in += k*( 1. - pow( 1 + alpha/n*( (mean-min)/sigma - alpha ),-m ) );
    }
  }

  return in;
}


// ------------------------------------------------------------------------
double ResolutionParametrization::CrystalBallFunction::truncPdf(double x, double mean, double sigma, double alpha, double n, double min) const {
  double p = 0.;
  double norm = integral(mean,sigma,alpha,n,min,10.);
  if( norm ) p = truncValue(x,mean,sigma,alpha,n,min)/norm;
  return p;
}


// ------------------------------------------------------------------------
double ResolutionParametrization::CrystalBallFunction::truncRandom(double mean, double sigma, double alpha, double n, double min) const {
  double max = truncPdf(mean,mean,sigma,alpha,n,min);
  double x = -1.;
  do {
    x = rand_->Uniform(0.,2.);
  } while( truncPdf(x,mean,sigma,alpha,n,min) < rand_->Uniform(0.,1.)*max );
  return x;
}
