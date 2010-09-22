//
//  $Id: Parametrization.cc,v 1.6 2010/07/22 13:58:30 mschrode Exp $
//
#include "Parametrization.h"


#include "TF1.h"
#include "TFile.h"
#include "TH1.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TMath.h"
#include "TRandom.h"
#include "TRandom3.h"



// ------------------------------------------------------------------------
SmearGaussAvePt::SmearGaussAvePt(double ptAveMin, double ptAveMax)
  : Parametrization(0,1,0,0),
    ptAveMin_(ptAveMin),
    ptAveMax_(ptAveMax) {
  assert( 0.0 <= ptAveMin_ && ptAveMin_ < ptAveMax_ );
  dMeasMax_ = 100.;
  print();
}


// ------------------------------------------------------------------------
void SmearGaussAvePt::update(const double * par) {
  std::cout << "Updating maximum dMeas to " << std::flush;
  dMeasMax_ = 2.*par[0]/sqrt(2.); // Two sigma in dMeas
  std::cout << dMeasMax_ << std::endl;
}

// ------------------------------------------------------------------------
double SmearGaussAvePt::pdfPtMeas(double ptMeas1, double ptMeas2, double ptTrue, const double *par) const {
  double pdf = 0.;
  double dMeas = 0.5*(ptMeas1 - ptMeas2);
  if( std::abs(dMeas) < dMeasMax_ ) {
    double s = sigma(par)/sqrt(2.);
    double u = dMeas/s;
    double norm = s*sqrt(M_PI*2.)*erf(dMeasMax_/sqrt(2.)/s);
    pdf = exp(-0.5*u*u)/norm;
  }
  return pdf;
}


// ------------------------------------------------------------------------
double SmearGaussAvePt::pdfResponse(double r, double ptTrue, const double *par) const {
  double s = sigma(par)/ptTrue;
  double u = (r - 1.)/s;
  return exp(-0.5*u*u)/sqrt(2.*M_PI)/s;
}



// ------------------------------------------------------------------------
double SmearGaussAvePt::pdfDijetAsym(double a, double ptTrue, const double *par) const {
  double s = sigma(par)/ptTrue/sqrt(2.);
  double u = a/s;
  return exp(-0.5*u*u)/sqrt(2.*M_PI)/s;
}



// ------------------------------------------------------------------------
void SmearGaussAvePt::print() const {
  std::cout << "Parametrization class '" << name() << "'\n";
  std::cout << "  " << ptAveMin_ << " < ptAve < " << ptAveMax_ << " GeV\n";
  std::cout << std::endl;
}



// ------------------------------------------------------------------------
SmearGaussPtBin::SmearGaussPtBin(double tMin, double tMax, double xMin, double xMax, const std::vector<double> &parScales, const std::vector<double> &startPar, const std::string &spectrum)
  : Parametrization(0,1,0,0),
    tMin_(tMin),
    tMax_(tMax),
    xMin_(xMin),
    xMax_(xMax),
    scale_(parScales) {
  assert( 0.0 <= tMin_ && tMin_ < tMax_ );
  assert( 0.0 <= xMin_ && xMin_ < xMax_ );
  assert( scale_.size() >= nJetPars() );
  assert( startPar.size() >= nJetPars() );
  dMeasMax_ = 1000.;
  print();

  TFile file(spectrum.c_str(),"READ");
  file.GetObject("hPtGen",hPdfPtTrue_);
  if( !hPdfPtTrue_ ) {
    std::cerr << "ERROR: No histogram found in file '" << file.GetName() << "'\n";
    exit(1);
  } else {
    std::cout << "Getting truth pdf from file '" << spectrum << "'\n";
    hPdfPtTrue_->SetDirectory(0);
    hPdfPtTrue_->SetName("hPdfPtTrue");
    int binMin = hPdfPtTrue_->FindBin(tMin_);
    int binMax = hPdfPtTrue_->FindBin(tMax_);
    if( hPdfPtTrue_->Integral(binMin,binMax,"width") )
      hPdfPtTrue_->Scale(1./hPdfPtTrue_->Integral(binMin,binMax,"width"));
  }
  file.Close();
  
  hashTablePdfPtTrue_ = new TH1D("hashTablePdfPtTrue_","",5000,tMin_-0.01,tMax_+0.01);
}

SmearGaussPtBin::~SmearGaussPtBin() { 
  delete hashTablePdfPtTrue_;
  if( hPdfPtTrue_ ) delete hPdfPtTrue_;
}


// ------------------------------------------------------------------------
void SmearGaussPtBin::update(const double * par) {
   std::cout << name() << ": Updating hashed parameters" << std::endl;
   std::cout << "  Updating maximum dMeas to " << std::flush;
   dMeasMax_ = 2.*scale_[0]*par[0]/sqrt(2.); // Two sigma in dMeas
   std::cout << dMeasMax_ << std::endl;
   hashPdfPtTrue(par);
}


// ------------------------------------------------------------------------
double SmearGaussPtBin::pdfPtMeas(double ptMeas1, double ptMeas2, double ptTrue, const double *par) const {
  double pdf = 0.;
   double dMeas = 0.5*(ptMeas1 - ptMeas2);
   if( std::abs(dMeas) < dMeasMax_ ) {
     double s = sigma(par)/sqrt(2.);
     double u1 = dMeas/s;
     double u2 = (0.5*(ptMeas1 + ptMeas2)-ptTrue)/s;
     double norm = M_PI*s*s*erf(dMeasMax_/sqrt(2.)/s)*( erf((xMax_-ptTrue)/sqrt(2.)/s) - erf((xMin_-ptTrue)/sqrt(2.)/s) );
     if( norm < 1E-3 ) {
       pdf = 0.;
     } else {
       pdf = exp(-0.5*u1*u1-0.5*u2*u2)/norm;
     }
   }
  return pdf;
}


// ------------------------------------------------------------------------
double SmearGaussPtBin::pdfPtTrue(double ptTrue, const double *par) const {
  return hashTablePdfPtTrue_->GetBinContent(hashTablePdfPtTrue_->FindBin(ptTrue));
}


// ------------------------------------------------------------------------
double SmearGaussPtBin::pdfResponse(double r, double ptTrue, const double *par) const {
  double s = sigma(par)/ptTrue;
  double u = (r - 1.)/s;
  double cut = (1.+erf(ptTrue/sqrt(2.)/s))/2.;
  return exp(-0.5*u*u)/sqrt(2.*M_PI)/s/cut;
}


// ------------------------------------------------------------------------
double SmearGaussPtBin::pdfResponseError(double r, double ptTrue, const double *par, const double *cov, const std::vector<int> &covIdx) const {
  // Store derivatives
  double s = sigma(par)/ptTrue;
  double u = (r-1.)/s;
  double df = pdfResponse(r,ptTrue,par)*(u*u - 1.)/s;

  // Calculate variance
  double var = df*df*scale_[0]*scale_[0]*cov[0]/ptTrue/ptTrue;

  // Return standard deviation
  return sqrt(var);
}


// ------------------------------------------------------------------------
double SmearGaussPtBin::pdfDijetAsym(double a, double ptTrue, const double *par) const {
  double s = sigma(par)/ptTrue/sqrt(2.);
  double u = a/s;
  return exp(-0.5*u*u)/sqrt(2.*M_PI)/s;
}


// ------------------------------------------------------------------------
void SmearGaussPtBin::hashPdfPtTrue(const double *par) const {
  std::cout << "  Hashing truth pdf with sigma = " << par[0] << std::endl;

  // Loop over tMin_ < ptTrue < tMax_ values
  for(int bin = 1; bin <= hashTablePdfPtTrue_->GetNbinsX(); bin++) {
    double ptTrue = hashTablePdfPtTrue_->GetBinCenter(bin);
    // Store (un-normalized) truth pdf for
    // this value of ptTrue in hash table
    hashTablePdfPtTrue_->SetBinContent(bin,pdfPtTrueNotNorm(ptTrue,par));
  }

  // Normalise values of truth pdf
  hashTablePdfPtTrue_->Scale(1./hashTablePdfPtTrue_->Integral("width"));
}



// ------------------------------------------------------------------------
double SmearGaussPtBin::pdfPtTrueNotNorm(double ptTrue, const double *par) const {
  double pdf = underlyingPdfPtTrue(ptTrue,par);

  // Convolution with cuts on ptAve
  double s = scale_[0]*par[0];
  s /= sqrt(2.);
  double c = sqrt(M_PI/2.)*s*( erf((xMax_-ptTrue)/s/sqrt(2.)) - erf((xMin_-ptTrue)/s/sqrt(2.)) );

  return c*pdf;
}


// ------------------------------------------------------------------------
double SmearGaussPtBin::underlyingPdfPtTrue(double ptTrue, const double *par) const {
  return hPdfPtTrue_->Interpolate(ptTrue);
}


// ------------------------------------------------------------------------
void SmearGaussPtBin::print() const {
  std::cout << "Parametrization class '" << name() << "'\n";
  std::cout << "  " << tMin_ << " < ptTrue < " << tMax_ << " GeV\n";
  std::cout << "  " << xMin_ << " < ptMeas < " << xMax_ << " GeV\n";
  std::cout << std::endl;
}



// ------------------------------------------------------------------------
SmearCrystalBallPtBin::SmearCrystalBallPtBin(double tMin, double tMax, double xMin, double xMax, double rMin,double rMax, const std::vector<double> &parScales, const std::vector<double> &startPar, const std::string &spectrum)
  : Parametrization(0,6,0,0),
    tMin_(tMin),
    tMax_(tMax),
    xMin_(xMin),
    xMax_(xMax),
    rMin_(rMin),
    rMax_(rMax),
    scale_(parScales),
    cb_(new CrystalBallFunction()) {
  assert( 0.0 <= tMin_ && tMin_ < tMax_ );
  assert( 0.0 <= xMin_ && xMin_ < xMax_ );
  assert( scale_.size() >= nJetPars() );
  
  print();

  TFile file(spectrum.c_str(),"READ");
  file.GetObject("hPtGen",hPdfPtTrue_);
  if( !hPdfPtTrue_ ) {
    std::cerr << "ERROR: No histogram found in file '" << file.GetName() << "'\n";
    exit(1);
  } else {
    std::cout << "Getting truth pdf from file '" << spectrum << "'\n";
    hPdfPtTrue_->SetDirectory(0);
    hPdfPtTrue_->SetName("hPdfPtTrue");
    int binMin = hPdfPtTrue_->FindBin(tMin_);
    int binMax = hPdfPtTrue_->FindBin(tMax_);
    if( hPdfPtTrue_->Integral(binMin,binMax,"width") )
      hPdfPtTrue_->Scale(1./hPdfPtTrue_->Integral(binMin,binMax,"width"));
  }
  file.Close();
  
  hashTablePdfPtTrue_ = new TH1D("hashTablePdfPtTrue_","",5000,tMin_,tMax_);
  hashPdfPtTrue(&(startPar.front()));
}

SmearCrystalBallPtBin::~SmearCrystalBallPtBin() {
  delete hashTablePdfPtTrue_;
  delete hPdfPtTrue_;
  delete cb_;
}


// ------------------------------------------------------------------------
double SmearCrystalBallPtBin::pdfPtMeas(double ptMeas1, double ptMeas2, double ptTrue, const double *par) const {
  double pdf = 0.;
//   if( xMin_ < ptMeas && ptMeas < xMax_ ) {
//     double norm = cb_->integral(ptTrue,sigma(par),alpha(par),n(par),xMin_,xMax_);
//     if( norm < 1E-10 ) norm = 1E-10;
//     pdf = cb_->value(ptMeas,ptTrue,sigma(par),alpha(par),n(par)) / norm;
//   }
  return pdf;
}


// // ------------------------------------------------------------------------
// double SmearCrystalBallPtBin::pdfPtMeasJet2(double ptMeas, double ptTrue, double pt3Rel, const double *par) const {
//   return cb_->pdf(ptMeas,ptTrue,sigma(par),alpha(par),n(par));
// }


// ------------------------------------------------------------------------
double SmearCrystalBallPtBin::pdfPtTrue(double ptTrue, const double *par) const {
  return hashTablePdfPtTrue_->GetBinContent(hashTablePdfPtTrue_->FindBin(ptTrue));
}


// ------------------------------------------------------------------------
double SmearCrystalBallPtBin::pdfResponse(double r, double ptTrue, const double *par) const {
  double pdf = 0.;
  double norm = cb_->integral(1.,sigma(par)/ptTrue,alpha(par),n(par),rMin_,rMax_);
  if( norm > 0. ) pdf = cb_->value(r,1.,sigma(par)/ptTrue,alpha(par),n(par))/norm;
  return pdf;
}


//! \note Only Gaussian core so far
// ------------------------------------------------------------------------
double SmearCrystalBallPtBin::pdfDijetAsym(double a, double ptTrue, const double *par) const {
  double s = sigma(par)/ptTrue/sqrt(2.);
  double u = a/s;
  return exp(-0.5*u*u)/sqrt(2.*M_PI)/s;
}


// ------------------------------------------------------------------------
void SmearCrystalBallPtBin::hashPdfPtTrue(const double *par) const {
  std::cout << "  Hashing truth pdf... " << std::flush;

  // Loop over tMin_ < ptTrue < tMax_ values
  for(int bin = 1; bin <= hashTablePdfPtTrue_->GetNbinsX(); bin++) {
    double ptTrue = hashTablePdfPtTrue_->GetBinCenter(bin);
    // Store (un-normalized) truth pdf for
    // this value of ptTrue in hash table
    hashTablePdfPtTrue_->SetBinContent(bin,pdfPtTrueNotNorm(ptTrue,par));
  }
  // Normalise values of truth pdf
  hashTablePdfPtTrue_->Scale(1./hashTablePdfPtTrue_->Integral("width"));

  std::cout << "ok" << std::endl;
}


// ------------------------------------------------------------------------
double SmearCrystalBallPtBin::pdfPtTrueNotNorm(double ptTrue, const double *par) const {
  double pdf = underlyingPdfPtTrue(ptTrue,par);

  // Description of cuts on 1. jet pt
  double s = sqrt( specSigmaPar(par,0)*specSigmaPar(par,0) +
		   specSigmaPar(par,1)*specSigmaPar(par,1)*ptTrue +
		   specSigmaPar(par,2)*specSigmaPar(par,2)*ptTrue*ptTrue );

  // Assuming Crystal Ball resolution
  double c = cb_->norm(ptTrue,s,alpha(par),n(par))*cb_->integral(ptTrue,s,alpha(par),n(par),xMin_,xMax_);

  return c*pdf;
}


// ------------------------------------------------------------------------
double SmearCrystalBallPtBin::underlyingPdfPtTrue(double ptTrue, const double *par) const {
  return hPdfPtTrue_->Interpolate(ptTrue);
}


// ------------------------------------------------------------------------
void SmearCrystalBallPtBin::print() const {
  std::cout << "Parametrization class '" << name() << "'\n";
  std::cout << "  " << tMin_ << " < ptTrue < " << tMax_ << " GeV\n";
  std::cout << "  " << xMin_ << " < ptMeas < " << xMax_ << " GeV\n";
  std::cout << std::endl;
}



//!  A Crystal Ball function representation
// ------------------------------------------------------------------------
Parametrization::CrystalBallFunction::CrystalBallFunction()
  : rand_(new TRandom3(0)) {}
	
Parametrization::CrystalBallFunction::~CrystalBallFunction() {
  delete rand_;
}


//! The Crystal Ball function value (not normalized)
// ------------------------------------------------------------------------
double Parametrization::CrystalBallFunction::value(double x, double mean, double sigma, double alpha, double n) const {
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
double Parametrization::CrystalBallFunction::norm(double mean, double sigma, double alpha, double n) const {
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
double Parametrization::CrystalBallFunction::random(double mean, double sigma, double alpha, double n) const {
  double max = pdf(mean,mean,sigma,alpha,n);
  double x = -1.;
  do {
    x = rand_->Uniform(0.,2.);
  } while( pdf(x,mean,sigma,alpha,n) < rand_->Uniform(0.,1.)*max );
  return x;
}


// ------------------------------------------------------------------------
double Parametrization::CrystalBallFunction::integral(double mean, double sigma, double alpha, double n, double min, double max) const {
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
double Parametrization::CrystalBallFunction::truncPdf(double x, double mean, double sigma, double alpha, double n, double min) const {
  double p = 0.;
  double norm = integral(mean,sigma,alpha,n,min,10.);
  if( norm ) p = truncValue(x,mean,sigma,alpha,n,min)/norm;
  return p;
}


// ------------------------------------------------------------------------
double Parametrization::CrystalBallFunction::truncRandom(double mean, double sigma, double alpha, double n, double min) const {
  double max = truncPdf(mean,mean,sigma,alpha,n,min);
  double x = -1.;
  do {
    x = rand_->Uniform(0.,2.);
  } while( truncPdf(x,mean,sigma,alpha,n,min) < rand_->Uniform(0.,1.)*max );
  return x;
}
