//
//  $Id: Parametrization.cc,v 1.5 2010/04/13 13:38:24 mschrode Exp $
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
SmearGauss::SmearGauss(double tMin, double tMax, double xMin, double xMax, const std::vector<double>& parScales, const std::vector<double> &startPar, const std::string &spectrum)
  : Parametrization(0,3,0,0),
    tMin_(tMin),
    tMax_(tMax),
    xMin_(xMin),
    xMax_(xMax),
    scale_(parScales),
    hPdfPtTrue_(0) {
  assert( 0.0 <= tMin_ && tMin_ < tMax_ );
  assert( 0.0 <= xMin_ && xMin_ < xMax_ );
  assert( scale_.size() >= nJetPars() );
  assert( startPar.size() >= nJetPars() );

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
  
  hashTablePdfPtTrue_ = new TH1D("hashTablePdfPtTrue_","",10000,tMin_,tMax_);
  hashPdfPtTrue(&(startPar.front()));

}


SmearGauss::~SmearGauss() { 
  delete hPdfPtTrue_;
  delete hashTablePdfPtTrue_;
 }


// ------------------------------------------------------------------------
double SmearGauss::pdfPtMeasJet1(double ptMeas, double ptTrue, double pt3Rel, const double *par) const {
  double pdf = 0.;
  if( xMin_ < ptMeas && ptMeas < xMax_ ) {
    double s = sigma(ptTrue,par);
    double u = (ptMeas - ptTrue)/s;
    double norm = sqrt(M_PI/2.)*s*( erf((xMax_-ptTrue)/sqrt(2.)/s) - erf((xMin_-ptTrue)/sqrt(2.)/s) );
    // This should be caught more cleverly
    if( norm < 1E-10 ) norm = 1E-10;
    pdf = exp(-0.5*u*u)/norm; 
  }
  return pdf;
}

// ------------------------------------------------------------------------
double SmearGauss::pdfPtMeasJet2(double ptMeas, double ptTrue, double pt3Rel, const double *par) const {
  double s = sigma(ptTrue,par);
  double u = (ptMeas - ptTrue)/s;
  double norm = sqrt(M_PI/2.)*s*( 1. + erf(ptTrue/sqrt(2.)/s) );
  // This should be caught more cleverly
  if( norm < 1E-10 ) norm = 1E-10;

  return exp(-0.5*u*u)/norm;
}


// ------------------------------------------------------------------------
double SmearGauss::pdfPtTrue(double ptTrue, const double *par) const {
  return hashTablePdfPtTrue_->GetBinContent(hashTablePdfPtTrue_->FindBin(ptTrue));
}


// ------------------------------------------------------------------------
double SmearGauss::pdfResponse(double r, double ptTrue, const double *par) const {
  double s = sigma(ptTrue,par)/ptTrue;
  double u = (r - 1.)/s;
  return exp(-0.5*u*u)/sqrt(2.*M_PI)/s;
}


// ------------------------------------------------------------------------
double SmearGauss::pdfDijetAsym(double a, double ptTrue, const double *par) const {
  double s = sigma(ptTrue,par)/ptTrue/sqrt(2.);
  double u = a/s;
  return exp(-0.5*u*u)/sqrt(2.*M_PI)/s;
}


// ------------------------------------------------------------------------
double SmearGauss::pdfResponseError(double r, double ptTrue, const double *par, const double *cov, const std::vector<int> &covIdx) const {
  // Store derivatives
  std::vector<double> df(nJetPars());
  for(size_t i = 0; i < nJetPars(); i++) {
    df[i] = pdfResponseDeriv(r,ptTrue,par,i);
  }

  // Calculate variance
  double var = 0.;
  for(int i = 0; i < static_cast<int>(nJetPars()); i++) { // Outer loop over parameters
    for(int j = 0; j < i+1; j++) { // Inner loop over parameters
      int idx = (i*i + i)/2 + j; // Index of (i,j) in covariance vector
      if( cov[idx] ) {
	if( i == j ) { // Diagonal terms
	  var += df[i]*df[i]*scale_[i]*scale_[i]*cov[idx];
	} else { // Off-diagonal terms
	  var += 2*df[i]*df[j]*scale_[i]*scale_[j]*cov[idx];
	}
      }
    } // End of inner loop over parameters
  } // End of outer loop over parameters
  // Return standard deviation
  return sqrt(var);
}


// ------------------------------------------------------------------------
double SmearGauss::pdfResponseDeriv(double r, double ptTrue, const double *par, int i) const {
  double df = 0.;
  if( i < 3 ) {
    double s = sigma(ptTrue,par);
    double u = ptTrue*(r-1.)/s;
    df = pdfResponse(r,ptTrue,par) * scale_[i]*par[i]/s/s * (u*u - 1.);
    if( i == 1 ) df *= ptTrue;
    if( i == 2 ) df *= ptTrue*ptTrue;
  }

  return df;
}


// ------------------------------------------------------------------------
void SmearGauss::hashPdfPtTrue(const double *par) const {
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
double SmearGauss::pdfPtTrueNotNorm(double ptTrue, const double *par) const {
  double pdf = underlyingPdfPtTrue(ptTrue,par);
  // Add description of cuts here

  return pdf;
}


// ------------------------------------------------------------------------
double SmearGauss::underlyingPdfPtTrue(double ptTrue, const double *par) const {
  return hPdfPtTrue_->Interpolate(ptTrue);
}




// ------------------------------------------------------------------------
void SmearGauss::print() const {
  std::cout << "Parametrization class '" << name() << "'\n";
  std::cout << "  Probability density of ptTrue spectrum:\n";
  std::cout << "    " << tMin_ << " < ptTrue < " << tMax_ << " GeV\n";
  std::cout << "    " << xMin_ << " < pt < " << xMax_ << " GeV\n";
  std::cout << std::endl;
}



// ------------------------------------------------------------------------
SmearGaussPtBin::SmearGaussPtBin(double tMin, double tMax, double xMin, double xMax, const std::vector<double> &parScales, const std::vector<double> &startPar, const std::string &spectrum)
  : Parametrization(0,4,0,0),
    tMin_(tMin),
    tMax_(tMax),
    xMin_(xMin),
    xMax_(xMax),
    scale_(parScales) {
  assert( 0.0 <= tMin_ && tMin_ < tMax_ );
  assert( 0.0 <= xMin_ && xMin_ < xMax_ );
  assert( scale_.size() >= nJetPars() );
  assert( startPar.size() >= nJetPars() );

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

SmearGaussPtBin::~SmearGaussPtBin() { 
  delete hashTablePdfPtTrue_;
  if( hPdfPtTrue_ ) delete hPdfPtTrue_;
}



// ------------------------------------------------------------------------
double SmearGaussPtBin::pdfPtMeasJet1(double ptMeas, double ptTrue, double pt3Rel, const double *par) const {
  double pdf = 0.;
  if( xMin_ < ptMeas && ptMeas < xMax_ ) {
    double s = sigma(par);
    double u = (ptMeas - ptTrue)/s;
    double norm = sqrt(M_PI/2.)*s*( erf((xMax_-ptTrue)/sqrt(2.)/s) - erf((xMin_-ptTrue)/sqrt(2.)/s) );
    // This should be caught more cleverly
    if( norm < 1E-10 ) norm = 1E-10;
   
    pdf = exp(-0.5*u*u)/norm; 
  }
  return pdf;
}



// ------------------------------------------------------------------------
double SmearGaussPtBin::pdfPtMeasJet2(double ptMeas, double ptTrue, double pt3Rel, const double *par) const {
  double s = sigma(par);
  double u = (ptMeas - ptTrue)/s;
  double norm = sqrt(M_PI/2.)*s*( 1. + erf(ptTrue/sqrt(2.)/s) );
  // This should be caught more cleverly
  if( norm < 1E-10 ) norm = 1E-10;
  return exp(-0.5*u*u)/norm;
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
double SmearGaussPtBin::pdfPtTrueNotNorm(double ptTrue, const double *par) const {
  double pdf = underlyingPdfPtTrue(ptTrue,par);
  // Convolution with cuts on 1. jet
  double s = sqrt( specSigmaPar(par,0)*specSigmaPar(par,0) +
		   specSigmaPar(par,1)*specSigmaPar(par,1)*ptTrue +
		   specSigmaPar(par,2)*specSigmaPar(par,2)*ptTrue*ptTrue );
  double c = 0.5*( erf((xMax_-ptTrue)/s/sqrt(2.)) - erf((xMin_-ptTrue)/s/sqrt(2.)) );

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
double SmearCrystalBallPtBin::pdfPtMeasJet1(double ptMeas, double ptTrue, double pt3Rel, const double *par) const {
  double pdf = 0.;
  if( xMin_ < ptMeas && ptMeas < xMax_ ) {
    double norm = cb_->integral(ptTrue,sigma(par),alpha(par),n(par),xMin_,xMax_);
    if( norm < 1E-10 ) norm = 1E-10;
    pdf = cb_->value(ptMeas,ptTrue,sigma(par),alpha(par),n(par)) / norm;
  }
  return pdf;
}


// ------------------------------------------------------------------------
double SmearCrystalBallPtBin::pdfPtMeasJet2(double ptMeas, double ptTrue, double pt3Rel, const double *par) const {
  return cb_->pdf(ptMeas,ptTrue,sigma(par),alpha(par),n(par));
}


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
